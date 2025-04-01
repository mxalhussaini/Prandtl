#include "DGSEMIntegrator.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{

DGSEMIntegrator::DGSEMIntegrator(
      std::shared_ptr<ParMesh> pmesh_,
      std::shared_ptr<ParFiniteElementSpace> fes0_,
      std::shared_ptr<ParGridFunction> alpha_,
      std::unique_ptr<LiftingScheme> liftingScheme_,
      NumericalFlux &rsolver_, int Np)
    : NonlinearFormIntegrator(), pmesh(pmesh_), fes0(fes0_), alpha(alpha_), liftingScheme(std::move(liftingScheme_)),
      rsolver(rsolver_), fluxFunction(rsolver_.GetFluxFunction()),
      Np_x(Np), Np_y(fluxFunction.dim > 1 ? Np : 1), Np_z(fluxFunction.dim > 2 ? Np : 1),
      num_equations(fluxFunction.num_equations), dim(num_equations - 2), num_elements(pmesh->GetNE()),
      GLIntRules(0, Quadrature1D::GaussLobatto)
{
    ir = &GLIntRules.Get(Geometry::SEGMENT, 2 * Np_x - 3);
    if (dim == 1)
    {
        ir_face = &GLIntRules.Get(Geometry::POINT, 2 * Np_x - 3);
        ir_vol = &GLIntRules.Get(Geometry::SEGMENT, 2 * Np_x - 3);
    }
    else if (dim == 2)
    {
        ir_face = &GLIntRules.Get(Geometry::SEGMENT, 2 * Np_x - 3);
        ir_vol = &GLIntRules.Get(Geometry::SQUARE, 2 * Np_x - 3);
    }
    else
    {
        ir_face = &GLIntRules.Get(Geometry::SQUARE, 2 * Np_x - 3);
        ir_vol = &GLIntRules.Get(Geometry::CUBE, 2 * Np_x - 3);
    }

    max_char_speed = -1.0;
    
    D_T.SetSize(Np_x);
    Dhat_T.SetSize(Np_x);
    Dhat2_T.SetSize(Np_x);

    Vector wBary(Np_x);
    wBary = 1.0;

    for (int i = 1; i < Np_x; i++)
    {
        for (int j = 0; j < i; j++)
        {
            wBary(j) *= (ir->IntPoint(j).x - ir->IntPoint(i).x);
            wBary(i) *= (ir->IntPoint(i).x - ir->IntPoint(j).x);
        }
    }
    
    wBary.Reciprocal();
    D_T = 0.0;
    for (int iL = 0; iL < Np_x; iL++)
    {
        for (int i = 0; i < Np_x; i++)
        {
            if (iL != i)
            {
                D_T(i, iL) = wBary(iL) / wBary(i) / (ir->IntPoint(i).x - ir->IntPoint(iL).x);
                D_T(i, i) -= D_T(i, iL);
            }
        }
    }
    
    Dhat_T = D_T;
    Dhat_T(0, 0) += 1.0 / ir->IntPoint(0).weight;
    Dhat_T(Np - 1, Np - 1) -= 1.0 / ir->IntPoint(Np - 1).weight;
    Dhat_T.Transpose();

    Dhat2_T = D_T;
    Dhat2_T *= 2.0;
    Dhat2_T(0, 0) += 1.0 / ir->IntPoint(0).weight;
    Dhat2_T(Np - 1, Np - 1) -= 1.0 / ir->IntPoint(Np - 1).weight;
    Dhat2_T.Transpose();

    D_T.Transpose();

    dof = dof1 = dof2 = ir_vol->GetNPoints();

    shape1.SetSize(ir_vol->GetNPoints());
    shape2.SetSize(ir_vol->GetNPoints());

    state1.SetSize(num_equations);
    state2.SetSize(num_equations);

    f.SetSize(num_equations);
    g.SetSize(num_equations);
    h.SetSize(num_equations);

    flux_num.SetSize(num_equations);
    flux_mat1.SetSize(num_equations, dim);
    flux_mat2.SetSize(num_equations, dim);
    flux_mat.SetSize(dim, num_equations);

    flux_mat1 = flux_mat2 = 0.0;

    dqdx.SetSize(num_equations);
    dqdy.SetSize(num_equations);
    dqdz.SetSize(num_equations);

    el_dudxi.SetSize(num_equations);
    el_dudeta.SetSize(num_equations);
    el_dudzeta.SetSize(num_equations);

    adj1.SetSize(dim); adj2.SetSize(dim);
    metric1.SetSize(dim); metric2.SetSize(dim);
    nor.SetSize(dim);

    F_inviscid.SetSize(num_equations, Np_x, Np_x * Np_y * Np_z);
    G_inviscid.SetSize(num_equations, Np_y, Np_x * Np_y * Np_z);
    H_inviscid.SetSize(num_equations, Np_z, Np_x * Np_y * Np_z);
    
    F_viscous.SetSize(num_equations, Np_x * Np_y * Np_z);
    G_viscous.SetSize(num_equations, Np_x * Np_y * Np_z);
    H_viscous.SetSize(num_equations, Np_x * Np_y * Np_z);

    D_row.SetSize(Np_x);

    dU_subcell.SetSize(num_equations);
    dU_inviscid.SetSize(num_equations);
    dU_viscous.SetSize(num_equations);
    dU_volume.SetSize(num_equations);
    dU_face1.SetSize(num_equations);
    dU_face2.SetSize(num_equations);
    dU.SetSize(num_equations);

    liftingScheme->SetLiftingParameters(ir, ir_face, ir_vol, D_T, num_equations, Np, dim);

#ifdef SUBCELL_FV_BLENDING
    SubcellMetricXi.SetSize(dim, Np_z * Np_y * (Np_x + 1), pmesh->GetNE());
    SubcellMetricEta.SetSize(dim, Np_z * (Np_y + 1) * Np_x, pmesh->GetNE());
    SubcellMetricZeta.SetSize(dim, (Np_z + 1) * Np_y * Np_x, pmesh->GetNE());
    ComputeSubcellMetrics();
#endif
}

void DGSEMIntegrator::AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                         FaceElementTransformations &Tr, const Vector &el_u,
                                         Vector &el_dudt)
{
    el_dudt.SetSize((dof1 + dof2) * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat1(el_u.GetData(), dof1, num_equations);
    const DenseMatrix el_u_mat2(el_u.GetData() + dof1 * num_equations, dof2, num_equations);
    
    DenseMatrix el_dudt_mat1(el_dudt.GetData(), dof1, num_equations);
    DenseMatrix el_dudt_mat2(el_dudt.GetData() + dof1 * num_equations, dof2, num_equations);

    for (int i = 0; i < ir_face->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir_face->IntPoint(i);
        Tr.SetAllIntPoints(&ip);
        J1 = Tr.GetElement1Transformation().Weight();
        J2 = Tr.GetElement2Transformation().Weight();
        el1.CalcShape(Tr.GetElement1IntPoint(), shape1);
        el2.CalcShape(Tr.GetElement2IntPoint(), shape2);

        el_u_mat1.MultTranspose(shape1, state1);
        el_u_mat2.MultTranspose(shape2, state2);        

        if (dim == 1)
        {
            nor(0) = (Tr.GetElement1IntPoint().x - 0.5) * 2.0;
        }
        else
        {
            CalcOrtho(Tr.Jacobian(), nor);
        }

        max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
        dU_face1 = dU_face2 = flux_num;
        dU_face1.Neg();

        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J1), shape1, dU_face1, el_dudt_mat1);
        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J2), shape2, dU_face2, el_dudt_mat2);
    }
}

void DGSEMIntegrator::AssembleElementVector(const FiniteElement &el,
        ElementTransformation &Tr, const Vector &el_u, Vector &el_dudt)
{
    fes0->GetElementDofs(Tr.ElementNo, alpha_indx);
    alpha->GetSubVector(alpha_indx, el_alpha);

    el_dudt.SetSize(dof * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat(el_u.GetData(), dof, num_equations);
    DenseMatrix el_dudt_mat(el_dudt.GetData(), dof, num_equations);
#ifdef SUBCELL_FV_BLENDING    
    ComputeFVFluxes(el_u_mat, el_alpha(0), Tr, el_dudt_mat);
#endif
    
    for (int k = 0; k < Np_z; k++)
    {
        for (int j = 0; j < Np_y; j++)
        {
            for (int i = 0; i < Np_x; i++)
            {
                id1 = k * Np_y * Np_x + j * Np_x + i;
                const IntegrationPoint &ip1 = ir_vol->IntPoint(id1);
                el_u_mat.GetRow(id1, state1);
                Tr.SetIntPoint(&ip1);
                J = Tr.Weight();
                adj1 = Tr.AdjugateJacobian();
                adj1.GetRow(0, metric1);
                f = 0.0;
                F_inviscid(id1).SetCol(i, f);

                for (int m = i + 1; m < Np_x; m++)
                {
                    id2 = k * Np_y * Np_x + j * Np_x + m;
                    const IntegrationPoint &ip2 = ir_vol->IntPoint(id2);
                    el_u_mat.GetRow(id2, state2);
                    Tr.SetIntPoint(&ip2);
                    adj2 = Tr.AdjugateJacobian();
                    adj2.GetRow(0, metric2);
                    max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, f));
                    F_inviscid(id1).SetCol(m, f);
                    F_inviscid(id2).SetCol(i, f);
                }

                if (dim > 1)
                {
                    adj1.GetRow(1, metric1);
                    g = 0.0;
                    G_inviscid(id1).SetCol(j, g);
                    for (int m = j + 1; m < Np_y; m++)
                    {
                        id2 = k * Np_y * Np_x + m * Np_x + i;
                        const IntegrationPoint &ip3 = ir_vol->IntPoint(id2);
                        el_u_mat.GetRow(id2, state2);
                        Tr.SetIntPoint(&ip3);
                        adj2 = Tr.AdjugateJacobian();
                        adj2.GetRow(1, metric2);
                        max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, g));
                        G_inviscid(id1).SetCol(m, g);
                        G_inviscid(id2).SetCol(j, g);
                    }
                    if (dim > 2)
                    {
                        adj1.GetRow(2, metric1);
                        h = 0.0;
                        H_inviscid(id1).SetCol(k, h);
                        for (int m = k + 1; m < Np_z; m++)
                        {
                            id2 = m * Np_y * Np_x + j * Np_x + i;
                            const IntegrationPoint &ip4 = ir_vol->IntPoint(id2);
                            el_u_mat.GetRow(id2, state2);
                            Tr.SetIntPoint(&ip4);
                            adj2 = Tr.AdjugateJacobian();
                            adj2.GetRow(2, metric2);
                            max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, h));
                            H_inviscid(id1).SetCol(m, h);
                            H_inviscid(id2).SetCol(k, g);
                        }
                    }
                }

                Dhat2_T.GetColumn(i, D_row); 
                F_inviscid(id1).Mult(D_row, dU_inviscid);

                if (dim > 1)
                {
                    Dhat2_T.GetColumn(j, D_row);
                    G_inviscid(id1).AddMult(D_row, dU_inviscid);
                    
                    if (dim > 2)
                    {
                        Dhat2_T.GetColumn(k, D_row);
                        H_inviscid(id1).AddMult(D_row, dU_inviscid);
                    }
                }
                
                dU_inviscid.Neg();
#ifdef SUBCELL_FV_BLENDING
                dU_inviscid *= (1.0 - el_alpha(0));
#endif
                dU_inviscid /= J;
                AddRow(el_dudt_mat, dU_inviscid, id1);
            }
        }
    }
}

void DGSEMIntegrator::AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                         FaceElementTransformations &Tr, const Vector &el_u,
                                         const Vector &el_dudx, const Vector &el_dudy,
                                         const Vector &el_dudz, Vector &el_dudt)
{
    el_dudt.SetSize((dof1 + dof2) * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat1(el_u.GetData(), dof1, num_equations);
    const DenseMatrix el_u_mat2(el_u.GetData() + dof1 * num_equations, dof2, num_equations);

    const DenseMatrix el_dudx_mat1(el_dudx.GetData(), dof1, num_equations);
    const DenseMatrix el_dudx_mat2(el_dudx.GetData() + dof1 * num_equations, dof2, num_equations);

    const DenseMatrix el_dudy_mat1(el_dudy.GetData(), dof1, num_equations);
    const DenseMatrix el_dudy_mat2(el_dudy.GetData() + dof1 * num_equations, dof2, num_equations);

    const DenseMatrix el_dudz_mat1(el_dudz.GetData(), dof1, num_equations);
    const DenseMatrix el_dudz_mat2(el_dudz.GetData() + dof1 * num_equations, dof2, num_equations);
    
    DenseMatrix el_dudt_mat1(el_dudt.GetData(), dof1, num_equations);
    DenseMatrix el_dudt_mat2(el_dudt.GetData() + dof1 * num_equations, dof2, num_equations);

    for (int i = 0; i < ir_face->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir_face->IntPoint(i);
        Tr.SetAllIntPoints(&ip);
        J1 = Tr.GetElement1Transformation().Weight();
        J2 = Tr.GetElement2Transformation().Weight();
        el1.CalcShape(Tr.GetElement1IntPoint(), shape1);
        el2.CalcShape(Tr.GetElement2IntPoint(), shape2);

        el_u_mat1.MultTranspose(shape1, state1);
        el_dudx_mat1.MultTranspose(shape1, dqdx);
        el_dudy_mat1.MultTranspose(shape1, dqdy);
        el_dudz_mat1.MultTranspose(shape1, dqdz);
        fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, dqdz, flux_mat1);

        el_u_mat2.MultTranspose(shape2, state2);
        el_dudx_mat2.MultTranspose(shape2, dqdx);
        el_dudy_mat2.MultTranspose(shape2, dqdy);
        el_dudz_mat2.MultTranspose(shape2, dqdz);   
        fluxFunction.ComputeViscousFlux(state2, dqdx, dqdy, dqdz, flux_mat2);

        CalcOrtho(Tr.Jacobian(), nor);

        max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
        dU_face1 = dU_face2 = flux_num;
        dU_face1.Neg();

        flux_mat1 += flux_mat2;
        flux_mat1 *= 0.5;
        flux_mat1.Mult(nor, flux_num);
        dU_face1 += flux_num;
        dU_face2 -= flux_num;  

        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J1), shape1, dU_face1, el_dudt_mat1);
        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J2), shape2, dU_face2, el_dudt_mat2);
    }
}

void DGSEMIntegrator::AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                         FaceElementTransformations &Tr, const Vector &el_u,
                                         const Vector &el_dudx, const Vector &el_dudy, Vector &el_dudt)
{
    el_dudt.SetSize((dof1 + dof2) * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat1(el_u.GetData(), dof1, num_equations);
    const DenseMatrix el_u_mat2(el_u.GetData() + dof1 * num_equations, dof2, num_equations);

    const DenseMatrix el_dudx_mat1(el_dudx.GetData(), dof1, num_equations);
    const DenseMatrix el_dudx_mat2(el_dudx.GetData() + dof1 * num_equations, dof2, num_equations);

    const DenseMatrix el_dudy_mat1(el_dudy.GetData(), dof1, num_equations);
    const DenseMatrix el_dudy_mat2(el_dudy.GetData() + dof1 * num_equations, dof2, num_equations);
    
    DenseMatrix el_dudt_mat1(el_dudt.GetData(), dof1, num_equations);
    DenseMatrix el_dudt_mat2(el_dudt.GetData() + dof1 * num_equations, dof2, num_equations);

    for (int i = 0; i < ir_face->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir_face->IntPoint(i);
        Tr.SetAllIntPoints(&ip);
        J1 = Tr.GetElement1Transformation().Weight();
        J2 = Tr.GetElement2Transformation().Weight();
        el1.CalcShape(Tr.GetElement1IntPoint(), shape1);
        el2.CalcShape(Tr.GetElement2IntPoint(), shape2);

        el_u_mat1.MultTranspose(shape1, state1);
        el_dudx_mat1.MultTranspose(shape1, dqdx);
        el_dudy_mat1.MultTranspose(shape1, dqdy);
        fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, flux_mat1);

        el_u_mat2.MultTranspose(shape2, state2);
        el_dudx_mat2.MultTranspose(shape2, dqdx);
        el_dudy_mat2.MultTranspose(shape2, dqdy);
        fluxFunction.ComputeViscousFlux(state2, dqdx, dqdy, flux_mat2);

        CalcOrtho(Tr.Jacobian(), nor);

        max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
        dU_face1 = dU_face2 = flux_num;
        dU_face1.Neg();

        flux_mat1 += flux_mat2;
        flux_mat1 *= 0.5;
        flux_mat1.Mult(nor, flux_num);
        dU_face1 += flux_num;
        dU_face2 -= flux_num;  

        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J1), shape1, dU_face1, el_dudt_mat1);
        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J2), shape2, dU_face2, el_dudt_mat2);
    }
}

void DGSEMIntegrator::AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                         FaceElementTransformations &Tr, const Vector &el_u,
                                         const Vector &el_dudx, Vector &el_dudt)
{
    el_dudt.SetSize((dof1 + dof2) * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat1(el_u.GetData(), dof1, num_equations);
    const DenseMatrix el_u_mat2(el_u.GetData() + dof1 * num_equations, dof2, num_equations);

    const DenseMatrix el_dudx_mat1(el_dudx.GetData(), dof1, num_equations);
    const DenseMatrix el_dudx_mat2(el_dudx.GetData() + dof1 * num_equations, dof2, num_equations);
    
    DenseMatrix el_dudt_mat1(el_dudt.GetData(), dof1, num_equations);
    DenseMatrix el_dudt_mat2(el_dudt.GetData() + dof1 * num_equations, dof2, num_equations);

    for (int i = 0; i < ir_face->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir_face->IntPoint(i);
        Tr.SetAllIntPoints(&ip);
        J1 = Tr.GetElement1Transformation().Weight();
        J2 = Tr.GetElement2Transformation().Weight();
        el1.CalcShape(Tr.GetElement1IntPoint(), shape1);
        el2.CalcShape(Tr.GetElement2IntPoint(), shape2);

        el_u_mat1.MultTranspose(shape1, state1);
        el_dudx_mat1.MultTranspose(shape1, dqdx);
        fluxFunction.ComputeViscousFlux(state1, dqdx, flux_mat1);

        el_u_mat2.MultTranspose(shape2, state2);
        el_dudx_mat2.MultTranspose(shape2, dqdx);
        fluxFunction.ComputeViscousFlux(state2, dqdx, flux_mat2);  

        nor(0) = (Tr.GetElement1IntPoint().x - 0.5) * 2.0;

        max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
        dU_face1 = dU_face2 = flux_num;
        dU_face1.Neg();

        flux_mat1 += flux_mat2;
        flux_mat1 *= 0.5;
        flux_mat1.Mult(nor, flux_num);
        dU_face1 += flux_num;
        dU_face2 -= flux_num;  

        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J1), shape1, dU_face1, el_dudt_mat1);
        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J2), shape2, dU_face2, el_dudt_mat2);
    }
}

void DGSEMIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, const Vector &el_dudx, const Vector &el_dudy, const Vector &el_dudz, Vector &el_dudt)
{
    fes0->GetElementDofs(Tr.ElementNo, alpha_indx);
    alpha->GetSubVector(alpha_indx, el_alpha);

    el_dudt.SetSize(dof * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat(el_u.GetData(), dof, num_equations);
    const DenseMatrix el_dudx_mat(el_dudx.GetData(), dof, num_equations);
    const DenseMatrix el_dudy_mat(el_dudy.GetData(), dof, num_equations);
    const DenseMatrix el_dudz_mat(el_dudz.GetData(), dof, num_equations);
    DenseMatrix el_dudt_mat(el_dudt.GetData(), dof, num_equations);
#ifdef SUBCELL_FV_BLENDING
    ComputeFVFluxes(el_u_mat, el_alpha(0), Tr, el_dudt_mat);
#endif
    for (int i = 0; i < ir_vol->GetNPoints(); i++)
    {
        const IntegrationPoint &ip1 = ir_vol->IntPoint(i);
        el_u_mat.GetRow(id1, state1);
        Tr.SetIntPoint(&ip1);
        adj1 = Tr.AdjugateJacobian();

        el_dudx_mat.GetRow(i, dqdx);

        el_dudy_mat.GetRow(i, dqdy);

        el_dudz_mat.GetRow(i, dqdz);      

        fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, dqdz, flux_mat1);

        mfem::MultABt(adj1, flux_mat1, flux_mat);
        flux_mat.GetRow(0, f);
        F_viscous.SetCol(i, f);

        flux_mat.GetRow(1, g);
        G_viscous.SetCol(i, g);

        flux_mat.GetRow(2, h);
        H_viscous.SetCol(i, h);        
    }
    
    for (int k = 0; k < Np_z; k++)
    {
        for (int j = 0; j < Np_y; j++)
        {
            for (int i = 0; i < Np_x; i++)
            {
                id1 = k * Np_y * Np_x + j * Np_x + i;
                const IntegrationPoint &ip1 = ir_vol->IntPoint(id1);
                el_u_mat.GetRow(id1, state1);
                Tr.SetIntPoint(&ip1);
                J = Tr.Weight();
                adj1 = Tr.AdjugateJacobian();

                adj1.GetRow(0, metric1);
                f = 0.0;
                F_inviscid(id1).SetCol(i, f);

                for (int m = i + 1; m < Np_x; m++)
                {
                    id2 = k * Np_y * Np_x + j * Np_x + m;
                    const IntegrationPoint &ip2 = ir_vol->IntPoint(id2);
                    el_u_mat.GetRow(id2, state2);
                    Tr.SetIntPoint(&ip2);
                    adj2 = Tr.AdjugateJacobian();
                    adj2.GetRow(0, metric2);
                    max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, f));
                    F_inviscid(id1).SetCol(m, f);
                    F_inviscid(id2).SetCol(i, f);
                }

                adj1.GetRow(1, metric1);
                g = 0.0;
                G_inviscid(id1).SetCol(j, g);
                for (int m = j + 1; m < Np_y; m++)
                {
                    id2 = k * Np_y * Np_x + m * Np_x + i;
                    const IntegrationPoint &ip3 = ir_vol->IntPoint(id2);
                    el_u_mat.GetRow(id2, state2);
                    Tr.SetIntPoint(&ip3);
                    adj2 = Tr.AdjugateJacobian();
                    adj2.GetRow(1, metric2);
                    max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, g));
                    G_inviscid(id1).SetCol(m, g);
                    G_inviscid(id2).SetCol(j, g);
                }

                adj1.GetRow(2, metric1);
                h = 0.0;
                H_inviscid(id1).SetCol(k, h);
                for (int m = k + 1; m < Np_z; m++)
                {
                    id2 = m * Np_y * Np_x + j * Np_x + i;
                    const IntegrationPoint &ip4 = ir_vol->IntPoint(id2);
                    el_u_mat.GetRow(id2, state2);
                    Tr.SetIntPoint(&ip4);
                    adj2 = Tr.AdjugateJacobian();
                    adj2.GetRow(2, metric2);
                    max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, h));
                    H_inviscid(id1).SetCol(m, h);
                    H_inviscid(id2).SetCol(k, g);
                }

                dU_viscous = 0.0;

                Dhat2_T.GetColumn(i, D_row); 
                F_inviscid(id1).Mult(D_row, dU_inviscid);

                Dhat_T.GetColumn(i, D_row);
                for (int l = 0; l < Np_x; l++)
                {
                    F_viscous.GetColumn(k * Np_x * Np_y + j * Np_x + l, dU);
                    dU *= D_row(l);
                    dU_viscous += dU;
                }

                Dhat2_T.GetColumn(j, D_row);
                G_inviscid(id1).AddMult(D_row, dU_inviscid);
                
                Dhat_T.GetColumn(j, D_row);
                for (int l = 0; l < Np_y; l++)
                {
                    G_viscous.GetColumn(k * Np_x * Np_y + l * Np_x + i, dU);
                    dU *= D_row(l);
                    dU_viscous += dU;
                }

                Dhat2_T.GetColumn(k, D_row);
                H_inviscid(id1).AddMult(D_row, dU_inviscid);

                Dhat_T.GetColumn(k, D_row);
                for (int l = 0; l < Np_z; l++)
                {
                    H_viscous.GetColumn(l * Np_x * Np_y + j * Np_x + i, dU);
                    dU *= D_row(l);
                    dU_viscous += dU;
                }
         
                dU_inviscid.Neg();
#ifdef SUBCELL_FV_BLENDING
                dU_inviscid *= (1.0 - el_alpha(0));
#endif
                add(dU_inviscid, dU_viscous, dU_volume);
                dU_volume /= J;
                AddRow(el_dudt_mat, dU_volume, id1);
            }
        }
    }    
}

void DGSEMIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, const Vector &el_dudx, const Vector &el_dudy, Vector &el_dudt)
{
    fes0->GetElementDofs(Tr.ElementNo, alpha_indx);
    alpha->GetSubVector(alpha_indx, el_alpha);

    el_dudt.SetSize(dof * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat(el_u.GetData(), dof, num_equations);
    const DenseMatrix el_dudx_mat(el_dudx.GetData(), dof, num_equations);
    const DenseMatrix el_dudy_mat(el_dudy.GetData(), dof, num_equations);
    DenseMatrix el_dudt_mat(el_dudt.GetData(), dof, num_equations);
#ifdef SUBCELL_FV_BLENDING
    ComputeFVFluxes(el_u_mat, el_alpha(0), Tr, el_dudt_mat);
#endif
    for (int i = 0; i < ir_vol->GetNPoints(); i++)
    {
        const IntegrationPoint &ip1 = ir_vol->IntPoint(i);
        el_u_mat.GetRow(id1, state1);
        Tr.SetIntPoint(&ip1);
        adj1 = Tr.AdjugateJacobian();

        el_dudx_mat.GetRow(i, dqdx);

        el_dudy_mat.GetRow(i, dqdy);

        fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, flux_mat1);

        mfem::MultABt(adj1, flux_mat1, flux_mat);
        flux_mat.GetRow(0, f);
        F_viscous.SetCol(i, f);

        flux_mat.GetRow(1, g);
        G_viscous.SetCol(i, g);
    }
    
    for (int j = 0; j < Np_y; j++)
    {
        for (int i = 0; i < Np_x; i++)
        {
            id1 = j * Np_x + i;
            const IntegrationPoint &ip1 = ir_vol->IntPoint(id1);
            el_u_mat.GetRow(id1, state1);
            Tr.SetIntPoint(&ip1);
            J = Tr.Weight();
            adj1 = Tr.AdjugateJacobian();
            adj1.GetRow(0, metric1);
            f = 0.0;
            F_inviscid(id1).SetCol(i, f);

            for (int m = i + 1; m < Np_x; m++)
            {
                id2 = j * Np_x + m;
                const IntegrationPoint &ip2 = ir_vol->IntPoint(id2);
                el_u_mat.GetRow(id2, state2);
                Tr.SetIntPoint(&ip2);
                adj2 = Tr.AdjugateJacobian();
                adj2.GetRow(0, metric2);
                max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, f));
                F_inviscid(id1).SetCol(m, f);
                F_inviscid(id2).SetCol(i, f);
            }

            adj1.GetRow(1, metric1);
            g = 0.0;
            G_inviscid(id1).SetCol(j, g);
            for (int m = j + 1; m < Np_y; m++)
            {
                id2 = m * Np_x + i;
                const IntegrationPoint &ip3 = ir_vol->IntPoint(id2);
                el_u_mat.GetRow(id2, state2);
                Tr.SetIntPoint(&ip3);
                adj2 = Tr.AdjugateJacobian();
                adj2.GetRow(1, metric2);
                max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, g));
                G_inviscid(id1).SetCol(m, g);
                G_inviscid(id2).SetCol(j, g);
            }

            dU_viscous = 0.0;

            Dhat2_T.GetColumn(i, D_row); 
            F_inviscid(id1).Mult(D_row, dU_inviscid);

            Dhat_T.GetColumn(i, D_row);
            for (int l = 0; l < Np_x; l++)
            {
                F_viscous.GetColumn(j * Np_x + l, dU);
                dU *= D_row(l);
                dU_viscous += dU;
            }

            Dhat2_T.GetColumn(j, D_row);
            G_inviscid(id1).AddMult(D_row, dU_inviscid);
            
            Dhat_T.GetColumn(j, D_row);
            for (int l = 0; l < Np_y; l++)
            {
                G_viscous.GetColumn(l * Np_x + i, dU);
                dU *= D_row(l);
                dU_viscous += dU;
            }
        
            dU_inviscid.Neg();
#ifdef SUBCELL_FV_BLENDING
            dU_inviscid *= (1.0 - el_alpha(0));
#endif
            add(dU_inviscid, dU_viscous, dU_volume);
            dU_volume /= J;
            AddRow(el_dudt_mat, dU_volume, id1);
        }
    }   
}

void DGSEMIntegrator::AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, const Vector &el_dudx, Vector &el_dudt)
{
    fes0->GetElementDofs(Tr.ElementNo, alpha_indx);
    alpha->GetSubVector(alpha_indx, el_alpha);

    el_dudt.SetSize(dof * num_equations);
    el_dudt = 0.0;

    const DenseMatrix el_u_mat(el_u.GetData(), dof, num_equations);
    const DenseMatrix el_dudx_mat(el_dudx.GetData(), dof, num_equations);
    DenseMatrix el_dudt_mat(el_dudt.GetData(), dof, num_equations);
#ifdef SUBCELL_FV_BLENDING
    ComputeFVFluxes(el_u_mat, el_alpha(0), Tr, el_dudt_mat);
#endif
    for (int i = 0; i < ir_vol->GetNPoints(); i++)
    {
        const IntegrationPoint &ip1 = ir_vol->IntPoint(i);
        el_u_mat.GetRow(id1, state1);
        Tr.SetIntPoint(&ip1);
        adj1 = Tr.AdjugateJacobian();

        el_dudx_mat.GetRow(i, dqdx);

        fluxFunction.ComputeViscousFlux(state1, dqdx, flux_mat1);

        mfem::MultABt(adj1, flux_mat1, flux_mat);
        flux_mat.GetRow(0, f);
        F_viscous.SetCol(i, f);
    }
    
    for (int i = 0; i < Np_x; i++)
    {
        const IntegrationPoint &ip1 = ir_vol->IntPoint(id1);
        el_u_mat.GetRow(id1, state1);
        Tr.SetIntPoint(&ip1);
        J = Tr.Weight();
        adj1 = Tr.AdjugateJacobian();
        adj1.GetRow(0, metric1);
        f = 0.0;
        F_inviscid(id1).SetCol(i, f);

        for (int m = i + 1; m < Np_x; m++)
        {
            const IntegrationPoint &ip2 = ir_vol->IntPoint(id2);
            el_u_mat.GetRow(id2, state2);
            Tr.SetIntPoint(&ip2);
            adj2 = Tr.AdjugateJacobian();
            adj2.GetRow(0, metric2);
            max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, f));
            F_inviscid(i).SetCol(m, f);
            F_inviscid(m).SetCol(i, f);
        }

        dU_viscous = 0.0;

        Dhat2_T.GetColumn(i, D_row); 
        F_inviscid(id1).Mult(D_row, dU_inviscid);

        Dhat_T.GetColumn(i, D_row);
        for (int l = 0; l < Np_x; l++)
        {
            F_viscous.GetColumn(l, dU);
            dU *= D_row(l);
            dU_viscous += dU;
        }
    
        dU_inviscid.Neg();
#ifdef SUBCELL_FV_BLENDING
        dU_inviscid *= (1.0 - el_alpha(0));
#endif
        add(dU_inviscid, dU_viscous, dU_volume);
        dU_volume /= J;
        AddRow(el_dudt_mat, dU_volume, id1);
    }   
}

void DGSEMIntegrator::AssembleLiftingFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy, Vector &el_dudz)
{
    liftingScheme->AssembleLiftingFaceVector(el1, el2, Tr, el_u, el_dudx, el_dudy, el_dudz);
}

void DGSEMIntegrator::AssembleLiftingFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy)
{
    liftingScheme->AssembleLiftingFaceVector(el1, el2, Tr, el_u, el_dudx, el_dudy);
}

void DGSEMIntegrator::AssembleLiftingFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, Vector &el_dudx)
{
    liftingScheme->AssembleLiftingFaceVector(el1, el2, Tr, el_u, el_dudx);
}

void DGSEMIntegrator::AssembleLiftingElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy, Vector &el_dudz)
{
    liftingScheme->AssembleLiftingElementVector(el, Tr, el_u, el_dudx, el_dudy, el_dudz);
}

void DGSEMIntegrator::AssembleLiftingElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy)
{
    liftingScheme->AssembleLiftingElementVector(el, Tr, el_u, el_dudx, el_dudy);
}

void DGSEMIntegrator::AssembleLiftingElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, Vector &el_dudx)
{
    liftingScheme->AssembleLiftingElementVector(el, Tr, el_u, el_dudx);
}

void DGSEMIntegrator::ComputeFVFluxes(const DenseMatrix &el_u_mat, real_t alpha_value, ElementTransformation &Tr, DenseMatrix &el_dudt_mat)
{
    for (int k = 0; k < Np_z; k++)
    {
        for (int j = 0; j < Np_y; j++)
        {
            dU_subcell = 0.0;
            id1 = k * Np_y * Np_x + j * Np_x;
            el_u_mat.GetRow(id1, state1);
            for (int i = 0; i < Np_x - 1; i++)
            {
                id2 = id1 + 1;
                el_u_mat.GetRow(id2, state2);
                SubcellMetricXi(Tr.ElementNo).GetColumn(id2, nor);
                max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
                Tr.SetIntPoint(&ir_vol->IntPoint(id1));

                dU_subcell -= flux_num;
                dU_subcell /= Tr.Weight() * (ir->IntPoint(i).weight);
                el_dudt_mat.SetRow(id1, dU_subcell);

                dU_subcell = flux_num;
                state1 = state2;
                id1 = id2;
            }
            Tr.SetIntPoint(&ir_vol->IntPoint(id1));
            dU_subcell /= Tr.Weight() * (ir->IntPoint(Np_x - 1).weight);
            el_dudt_mat.SetRow(id1, dU_subcell);
        }
    }

    if (dim > 1)
    {
        for (int k = 0; k < Np_z; k++)
        {
            for (int i = 0; i < Np_x; i++)
            {
                dU_subcell = 0.0;
                id1 = k * Np_y * Np_x + i;
                el_u_mat.GetRow(id1, state1);
                for (int j = 0; j < Np_y - 1; j++)
                {
                    id2 = k * Np_y * Np_x + (j + 1) * Np_x + i;
                    el_u_mat.GetRow(id2, state2);
                    SubcellMetricEta(Tr.ElementNo).GetColumn(id2, nor);
                    max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
                    Tr.SetIntPoint(&ir_vol->IntPoint(id1));

                    dU_subcell -= flux_num;
                    dU_subcell /= Tr.Weight() * (ir->IntPoint(j).weight);
                    AddRow(el_dudt_mat, dU_subcell, id1);

                    dU_subcell = flux_num;
                    state1 = state2;
                    id1 = id2;         
                }
                Tr.SetIntPoint(&ir_vol->IntPoint(id1));
                dU_subcell /= Tr.Weight() * (ir->IntPoint(Np_y - 1).weight);
                AddRow(el_dudt_mat, dU_subcell, id1);
            }
        }

        if (dim > 2)
        {
            for (int j = 0; j < Np_y; j++)
            {
                for (int i = 0; i < Np_x; i++)
                {
                    dU_subcell = 0.0;
                    id1 = j * Np_x + i;
                    el_u_mat.GetRow(id1, state1);
                    for (int k = 0; k < Np_z - 1; k++)
                    {
                        id2 = (k + 1) * Np_y * Np_x + j * Np_x + i;
                        el_u_mat.GetRow(id2, state2);
                        SubcellMetricZeta(Tr.ElementNo).GetColumn(id2, nor);
                        max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
                        Tr.SetIntPoint(&ir_vol->IntPoint(id1));

                        dU_subcell -= flux_num;
                        dU_subcell /= Tr.Weight() * (ir->IntPoint(k).weight);
                        AddRow(el_dudt_mat, dU_subcell, id1);

                        dU_subcell = flux_num;
                        state1 = state2;
                        id1 = id2;            
                    }
                    Tr.SetIntPoint(&ir_vol->IntPoint(id1));
                    dU_subcell /= Tr.Weight() * (ir->IntPoint(Np_z - 1).weight);
                    AddRow(el_dudt_mat, dU_subcell, id1);
                }
            }
        }
    }
    el_dudt_mat *= alpha_value;
}


void DGSEMIntegrator::ComputeSubcellMetrics()
{
    for (int el = 0; el < pmesh->GetNE(); el++)
    {
        ElementTransformation *Tr = pmesh->GetElementTransformation(el);
        DenseMatrix &nor_mat_xi = SubcellMetricXi(el);
        DenseMatrix &nor_mat_eta = SubcellMetricEta(el);
        DenseMatrix &nor_mat_zeta = SubcellMetricZeta(el);
        Vector tmp(dim);
        for (int k = 0; k < Np_z; k++)
        {
            for (int j = 0; j < Np_y; j++)
            {
                const IntegrationPoint &ip = ir_vol->IntPoint(k * Np_y * Np_x + j * Np_x); // left xi-face
                Tr->SetIntPoint(&ip);
                Tr->AdjugateJacobian().GetRow(0, metric1);
                for (int i = 0; i < Np_x + 1; i++)
                {
                    nor = metric1;
                    for (int l = 0; l < i; l++)
                    {
                        D_T.GetColumn(l, D_row);
                        tmp = 0.0;
                        real_t weight = ir->IntPoint(l).weight;
                        for (int m = 0; m < Np_x; m++)
                        {
                            Tr->SetIntPoint(&ir_vol->IntPoint(k * Np_y * Np_x + j * Np_x + m));
                            Tr->AdjugateJacobian().GetRow(0, metric2);
                            metric2 *= D_row(m);
                            tmp += metric2;
                        }
                        tmp *= weight;
                        nor += tmp;
                    }
                    nor_mat_xi.SetCol(k * Np_y * Np_x + j * Np_x + i, nor);
                }
            }
        }

        if (dim > 1)
        {
            for (int k = 0; k < Np_z; k++)
            {
                for (int i = 0; i < Np_x; i++)
                {
                    const IntegrationPoint &ip = ir_vol->IntPoint(k * Np_y * Np_x + i); // bottom eta-face
                    Tr->SetIntPoint(&ip);
                    Tr->AdjugateJacobian().GetRow(1, metric1);
                    for (int j = 0; j < Np_y + 1; j++)
                    {
                        nor = metric1;
                        for (int l = 0; l < j; l++)
                        {
                            D_T.GetColumn(l, D_row);
                            tmp = 0.0;
                            real_t weight = ir->IntPoint(l).weight;
                            for (int m = 0; m < Np_y; m++)
                            {
                                Tr->SetIntPoint(&ir_vol->IntPoint(k * Np_y * Np_x + m * Np_x + i));
                                Tr->AdjugateJacobian().GetRow(1, metric2);
                                metric2 *= D_row(m);
                                tmp += metric2;
                            }
                            tmp *= weight;
                            nor += tmp;
                        }
                        nor_mat_eta.SetCol(k * Np_y * Np_x + j * Np_x + i, nor);
                    }
                }
            }

            if (dim > 2)
            {
                for (int j = 0; j < Np_y; j++)
                {
                    for (int i = 0; i < Np_x; i++)
                    {
                        const IntegrationPoint &ip = ir_vol->IntPoint(j * Np_x + i); // bottom zeta-face
                        Tr->SetIntPoint(&ip);
                        Tr->AdjugateJacobian().GetRow(2, metric1);
                        for (int k = 0; k < Np_z + 1; k++)
                        {
                            nor = metric1;
                            for (int l = 0; l < k; l++)
                            {
                                D_T.GetColumn(l, D_row);
                                tmp = 0.0;
                                real_t weight = ir->IntPoint(l).weight;
                                for (int m = 0; m < Np_z; m++)
                                {
                                    Tr->SetIntPoint(&ir_vol->IntPoint(m * Np_y * Np_x + j * Np_x + i));
                                    Tr->AdjugateJacobian().GetRow(2, metric2);
                                    metric2 *= D_row(m);
                                    tmp += metric2;
                                }
                                tmp *= weight;
                                nor += tmp;
                            }
                            nor_mat_zeta.SetCol(k * Np_y * Np_x + j * Np_x + i, nor);
                        }
                    }
                }
            }
        }
    }   
}


}