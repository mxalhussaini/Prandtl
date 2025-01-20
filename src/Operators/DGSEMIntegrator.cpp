#include "DGSEMIntegrator.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{

DGSEMIntegrator::DGSEMIntegrator(
      std::shared_ptr<ParMesh> pmesh_,
      std::shared_ptr<ParFiniteElementSpace> vfes_, std::shared_ptr<ParFiniteElementSpace> fes0_,
      std::shared_ptr<ParGridFunction> grad_x_, std::shared_ptr<ParGridFunction> grad_y_,
      std::shared_ptr<ParGridFunction> grad_z_, std::shared_ptr<ParGridFunction> alpha_,
      NumericalFlux &rsolver_, int Np)
    : BilinearFormIntegrator(), pmesh(pmesh_), vfes(vfes_), fes0(fes0_),
      grad_x(grad_x_), grad_y(grad_y_), grad_z(grad_z_), alpha(alpha_),
      rsolver(rsolver_), fluxFunction(rsolver_.GetFluxFunction()),
      Np_x(Np), Np_y(fluxFunction.dim > 1 ? Np : 1), Np_z(fluxFunction.dim > 2 ? Np : 1),
      num_equations(fluxFunction.num_equations), dim(num_equations - 2), num_elements(pmesh->GetNE()),
      fe(new L2_SegmentElement(Np_x - 1, BasisType::GaussLobatto)),
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
    
    D.SetSize(Np_x);
    D_mod.SetSize(Np_x);
    D2_mod.SetSize(Np_x);

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
    D = 0.0;
    for (int iL = 0; iL < Np_x; iL++)
    {
        for (int i = 0; i < Np_x; i++)
        {
            if (iL != i)
            {
                D(i, iL) = wBary(iL) / wBary(i) / (ir->IntPoint(i).x - ir->IntPoint(iL).x);
                D(i, i) -= D(i, iL);
            }
        }
    }
    
    D_mod = D;
    D_mod(0, 0) += 1.0 / ir->IntPoint(0).weight;
    D_mod(Np - 1, Np - 1) -= 1.0 / ir->IntPoint(Np - 1).weight;
    D_mod.Transpose();

    D2_mod = D;
    D2_mod *= 2.0;
    D2_mod(0, 0) += 1.0 / ir->IntPoint(0).weight;
    D2_mod(Np - 1, Np - 1) -= 1.0 / ir->IntPoint(Np - 1).weight;
    D2_mod.Transpose();

    state1.SetSize(num_equations);
    state2.SetSize(num_equations);

    shape1.SetSize(ir_vol->GetNPoints());
    shape2.SetSize(ir_vol->GetNPoints());

    f.SetSize(num_equations);
    g.SetSize(num_equations);
    h.SetSize(num_equations);

    flux_num.SetSize(num_equations);
    flux1.SetSize(num_equations);
    flux2.SetSize(num_equations);
    flux_mat1.SetSize(num_equations, dim);
    flux_mat2.SetSize(num_equations, dim);
    flux_mat_ref.SetSize(dim, num_equations);

    prim.SetSize(num_equations);
    grad_state.SetSize(num_equations);
    grad_mat1.SetSize(num_equations, dim);
    grad_mat2.SetSize(num_equations, dim);

    adj1.SetSize(dim); adj2.SetSize(dim);
    metric1.SetSize(dim); metric2.SetSize(dim);
    nor.SetSize(dim); nor1.SetSize(dim); nor2.SetSize(dim); unit_nor.SetSize(dim);

    F_inviscid.SetSize(num_equations, Np_x, Np_x * Np_y * Np_z);
    G_inviscid.SetSize(num_equations, Np_y, Np_x * Np_y * Np_z);
    H_inviscid.SetSize(num_equations, Np_z, Np_x * Np_y * Np_z);
    
    F_viscous.SetSize(num_equations, Np_x * Np_y * Np_z);
    G_viscous.SetSize(num_equations, Np_x * Np_y * Np_z);
    H_viscous.SetSize(num_equations, Np_x * Np_y * Np_z);

    Dx.SetSize(Np_x);

    dU_sub.SetSize(num_equations);
    dU_inviscid.SetSize(num_equations);
    dU_viscous.SetSize(num_equations);
    dU.SetSize(num_equations);


    SubcellMetricXi.SetSize(dim, Np_z * Np_y * (Np_x + 1), pmesh->GetNE());
    SubcellMetricEta.SetSize(dim, Np_z * (Np_y + 1) * Np_x, pmesh->GetNE());
    SubcellMetricZeta.SetSize(dim, (Np_z + 1) * Np_y * Np_x, pmesh->GetNE());

    ComputeSubcellMetrics();
}

void DGSEMIntegrator::AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2,
                                         FaceElementTransformations &Tr, const Vector &elfun,
                                         Vector &elvect)
{
    const int dof1 = el1.GetDof();
    const int dof2 = el2.GetDof();

    shape1.SetSize(dof1);
    shape2.SetSize(dof2);

    elvect.SetSize((dof1 + dof2) * num_equations);
    elvect = 0.0;

    const DenseMatrix elfun1_mat(elfun.GetData(), dof1, num_equations);
    const DenseMatrix elfun2_mat(elfun.GetData() + dof1 * num_equations, dof2, num_equations);
    
    DenseMatrix elvect1_mat(elvect.GetData(), dof1, num_equations);
    DenseMatrix elvect2_mat(elvect.GetData() + dof1 * num_equations, dof2, num_equations);

#ifdef PARABOLIC
    if (currentMode == DIVERGENCE)
    {
        bool is_face_nbr = Tr.Elem2No >= num_elements;
        int elem2_no = is_face_nbr ? (Tr.Elem2No - num_elements) : Tr.Elem2No;
        vfes->GetElementVDofs(Tr.Elem1No, grad_indices1);
        if (!is_face_nbr)
        {
            vfes->GetElementVDofs(elem2_no, grad_indices2);
        }
        else
        {
            vfes->GetFaceNbrElementVDofs(elem2_no, grad_indices2);
        }

        auto extractGrads = [&] (std::shared_ptr<ParGridFunction> grad, Vector &grad_vdofs, const Array<int> &grad_indices, bool is_face_nbr)
        {
            if (!is_face_nbr)
            {
                grad->GetSubVector(grad_indices, grad_vdofs);
            }
            else
            {
                grad->FaceNbrData().GetSubVector(grad_indices, grad_vdofs);
            }
        };


        extractGrads(grad_x, grad_x_vdofs1, grad_indices1, false);
        extractGrads(grad_x, grad_x_vdofs2, grad_indices2, is_face_nbr);

        grad_x_mat1.UseExternalData(grad_x_vdofs1.GetData(), dof1, num_equations);
        grad_x_mat2.UseExternalData(grad_x_vdofs2.GetData(), dof2, num_equations);

        if (dim > 1)
        {
            extractGrads(grad_y, grad_y_vdofs1, grad_indices1, false);
            extractGrads(grad_y, grad_y_vdofs2, grad_indices2, is_face_nbr);
        
            grad_y_mat1.UseExternalData(grad_y_vdofs1.GetData(), dof1, num_equations);
            grad_y_mat2.UseExternalData(grad_y_vdofs2.GetData(), dof2, num_equations);

            if (dim > 2)
            {
                extractGrads(grad_z, grad_z_vdofs1, grad_indices1, false);
                extractGrads(grad_z, grad_z_vdofs2, grad_indices2, is_face_nbr);
            
                grad_z_mat1.UseExternalData(grad_z_vdofs1.GetData(), dof1, num_equations);
                grad_z_mat2.UseExternalData(grad_z_vdofs2.GetData(), dof2, num_equations);
            }
        }
    }
#endif

    for (int i = 0; i < ir_face->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir_face->IntPoint(i);
        Tr.SetAllIntPoints(&ip);
        J1 = Tr.GetElement1Transformation().Weight();
        J2 = Tr.GetElement2Transformation().Weight();
        el1.CalcShape(Tr.GetElement1IntPoint(), shape1);
        el2.CalcShape(Tr.GetElement2IntPoint(), shape2);

        elfun1_mat.MultTranspose(shape1, state1);
        elfun2_mat.MultTranspose(shape2, state2);        

        if (dim == 1)
        {
            nor(0) = (Tr.GetElement1IntPoint().x - 0.5) * 2.0;
        }
        else
        {
            CalcOrtho(Tr.Jacobian(), nor);
        }

        if (currentMode == DIVERGENCE)
        {
            max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor, flux_num));
            flux1 = flux_num;
            flux2 = flux_num;
            flux1.Neg();
#ifdef PARABOLIC
            grad_x_mat1.MultTranspose(shape1, grad_state);
            grad_mat1.SetCol(0, grad_state);

            grad_x_mat2.MultTranspose(shape2, grad_state);
            grad_mat2.SetCol(0, grad_state);

            if (dim > 1)
            {
                grad_y_mat1.MultTranspose(shape1, grad_state);
                grad_mat1.SetCol(1, grad_state);

                grad_y_mat2.MultTranspose(shape2, grad_state);
                grad_mat2.SetCol(1, grad_state);

                if (dim > 2)
                {
                    grad_z_mat1.MultTranspose(shape1, grad_state);
                    grad_mat1.SetCol(2, grad_state);

                    grad_z_mat2.MultTranspose(shape2, grad_state);
                    grad_mat2.SetCol(2, grad_state);
                }
            }

            fluxFunction.ComputeViscousFlux(state1, grad_mat1, flux_mat1);
            fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat2);
            // add viscous eigenvalue estimates
            flux_mat1 += flux_mat2;
            flux_mat1 *= 0.5;
            flux_mat1.Mult(nor, flux_num);
            flux1 += flux_num;
            flux2 -= flux_num;
#endif   
        }
        else
        {
            flux1 = state2;
            flux1 -= state1;
            flux1 *= 0.5 * nor(dir);
            flux2 = flux1;
        }

        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J1), shape1, flux1, elvect1_mat);
        AddMult_a_VWt(+1.0 / (ir->IntPoint(0).weight * J2), shape2, flux2, elvect2_mat);
    }
}

void DGSEMIntegrator::AssembleElementMatrix2(const FiniteElement &trial_fe, const FiniteElement &test_fe,
   ElementTransformation &Trans, DenseMatrix &elmat)
{
    int dof = trial_fe.GetDof();
    Vector d_col;

    DenseMatrix dshape(dof, dim);
    DenseMatrix gshape(dof, dim);
    DenseMatrix Jinv(dim);
    Vector shape(dof);
    elmat.SetSize(dim * dof, dof);
    DenseMatrix elmat_comp(dof, dof);

    elmat = 0.0;

    for (int i = 0; i < ir_vol->GetNPoints(); i++)
    {
        const IntegrationPoint &ip = ir_vol->IntPoint(i);
        Trans.SetIntPoint(&ip);

        CalcInverse(Trans.Jacobian(), Jinv);

        test_fe.CalcPhysShape(Trans, shape);
        trial_fe.CalcDShape(ip, dshape);

        Mult(dshape, Jinv, gshape);

        for (int d = 0; d < dim; ++d)
        {
            gshape.GetColumnReference(d, d_col);
            MultVWt(shape, d_col, elmat_comp);
            for (int jj = 0; jj < dof; ++jj)
            {
                for (int ii = 0; ii < dof; ++ii)
                {
                    elmat(d * dof + ii, jj) += elmat_comp(ii, jj);
                }
            }
        }
    }
}

void DGSEMIntegrator::AssembleElementVector(const FiniteElement &el,
        ElementTransformation &Tr, const Vector &elfun, Vector &elvect)
{
    const int dof = el.GetDof();

    fes0->GetElementDofs(Tr.ElementNo, alpha_indx);
    alpha->GetSubVector(alpha_indx, alpha_dof);

    elvect.SetSize(dof * num_equations);
    elvect = 0.0;

    const DenseMatrix elfun_mat(elfun.GetData(), dof, num_equations);
    DenseMatrix elvect_mat(elvect.GetData(), dof, num_equations);

    ComputeFVFluxes(elfun_mat, alpha_dof(0), Tr, elvect_mat);

#ifdef PARABOLIC
    vfes->GetElementVDofs(Tr.ElementNo, grad_indices1);

    grad_x->GetSubVector(grad_indices1, grad_x_vdofs1);
    grad_x_mat1.UseExternalData(grad_x_vdofs1.GetData(), dof, num_equations);

    if (dim > 1)
    {
        grad_y->GetSubVector(grad_indices1, grad_y_vdofs1);
        grad_y_mat1.UseExternalData(grad_y_vdofs1.GetData(), dof, num_equations);
        if (dim > 2)
        {
            grad_z->GetSubVector(grad_indices1, grad_z_vdofs1);
            grad_z_mat1.UseExternalData(grad_z_vdofs1.GetData(), dof, num_equations);
        }
    }

    for (int i = 0; i < ir_vol->GetNPoints(); i++)
    {
        const IntegrationPoint &ip1 = ir_vol->IntPoint(i);
        elfun_mat.GetRow(id1, state1);
        Tr.SetIntPoint(&ip1);
        adj1 = Tr.AdjugateJacobian();

        grad_x_mat1.GetRow(i, grad_state);
        grad_mat1.SetCol(0, grad_state);

        if (dim > 1)
        {
            grad_y_mat1.GetRow(i, grad_state);
            grad_mat1.SetCol(1, grad_state);

            if (dim > 2)
            {
                grad_z_mat1.GetRow(i, grad_state);
                grad_mat1.SetCol(2, grad_state);
            }
        }

        fluxFunction.ComputeViscousFlux(state1, grad_mat1, flux_mat1);
        // add viscous eigenvalue estimates

        mfem::MultABt(adj1, flux_mat1, flux_mat_ref);
        flux_mat_ref.GetRow(0, f);
        F_viscous.SetCol(i, f);
        if (dim > 1)
        {
            flux_mat_ref.GetRow(1, g);
            G_viscous.SetCol(i, g);
            if (dim > 2)
            {
                flux_mat_ref.GetRow(2, h);
                H_viscous.SetCol(i, h);
            }
        }
    }
    // grad_x_mat1.ClearExternalData();
    // if (dim > 1)
    // {
    //     grad_y_mat1.ClearExternalData();
    //     if (dim > 2)
    //     {
    //         grad_z_mat1.ClearExternalData();
    //     }
    // }
#endif
    
    for (int k = 0; k < Np_z; k++)
    {
        for (int j = 0; j < Np_y; j++)
        {
            for (int i = 0; i < Np_x; i++)
            {
                id1 = k * Np_y * Np_x + j * Np_x + i;
                const IntegrationPoint &ip1 = ir_vol->IntPoint(id1);
                elfun_mat.GetRow(id1, state1);
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
                    elfun_mat.GetRow(id2, state2);
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
                        elfun_mat.GetRow(id2, state2);
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
                            elfun_mat.GetRow(id2, state2);
                            Tr.SetIntPoint(&ip4);
                            adj2 = Tr.AdjugateJacobian();
                            adj2.GetRow(2, metric2);
                            max_char_speed = std::max(max_char_speed, rsolver.ComputeVolumeFlux(state1, state2, metric1, metric2, h));
                            H_inviscid(id1).SetCol(m, h);
                            H_inviscid(id2).SetCol(k, g);
                        }
                    }
                }

                elvect_mat.GetRow(id1, dU_sub);

                D2_mod.GetColumn(i, Dx); 
                F_inviscid(id1).Mult(Dx, dU_inviscid);
#ifdef PARABOLIC
                dU_viscous = 0.0;
                D_mod.GetColumn(i, Dx); //maybe transpose for better memory efficiency?
                for (int l = 0; l < Np_x; l++)
                {
                    F_viscous.GetColumn(k * Np_x * Np_y + j * Np_x + l, dU);
                    dU *= Dx(l);
                    dU_viscous += dU;
                }
#endif

                if (dim > 1)
                {
                    D2_mod.GetColumn(j, Dx);
                    G_inviscid(id1).Mult(Dx, dU);
                    dU_inviscid += dU;
#ifdef PARABOLIC
                    D_mod.GetColumn(j, Dx);
                    for (int l = 0; l < Np_y; l++)
                    {
                        G_viscous.GetColumn(k * Np_x * Np_y + l * Np_x + i, dU);
                        dU *= Dx(l);
                        dU_viscous += dU;
                    }
#endif
                    
                    if (dim > 2)
                    {
                        D2_mod.GetColumn(k, Dx);
                        H_inviscid(id1).Mult(Dx, dU);
                        dU_inviscid += dU;
#ifdef PARABOLIC
                        D_mod.GetColumn(k, Dx);
                        for (int l = 0; l < Np_z; l++)
                        {
                            H_viscous.GetColumn(l * Np_x * Np_y + j * Np_x + i, h);
                            dU = h;
                            dU *= Dx(l);
                            dU_viscous += dU;
                        }
#endif
                    }
                }
                
                dU_inviscid.Neg();
                dU_inviscid *= (1.0 - alpha_dof(0));

                dU_inviscid += dU_sub;
#ifdef PARABOLIC    
                dU_inviscid += dU_viscous;
#endif
                dU_inviscid /=J;
                elvect_mat.SetRow(id1, dU_inviscid);
            }
        }
    }
}

void DGSEMIntegrator::ComputeFVFluxes(const DenseMatrix &vdof_mat, const real_t &alpha_val, ElementTransformation &Tr, DenseMatrix &elvect_mat)
{
    for (int k = 0; k < Np_z; k++)
    {
        for (int j = 0; j < Np_y; j++)
        {
            flux1 = 0.0;
            for (int i = 0; i < Np_x - 1; i++)
            {
                int id1 = k * Np_y * Np_x + j * Np_x + i;
                int id2 = id1 + 1;
                vdof_mat.GetRow(id1, state1); // these two lines can be improved?
                vdof_mat.GetRow(id2, state2); // maybe no need to extract both states - one can be saved from previous iteration?
                SubcellMetricXi(Tr.ElementNo).GetColumn(id2, nor1);
                max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor1, flux_num));

                flux1 -= flux_num;
                flux1 *= alpha_val / (ir->IntPoint(i).weight);
                elvect_mat.SetRow(id1, flux1);

                flux1 = flux_num;
                // state1 = state2; like so?
            }
            id1 = k * Np_y * Np_x + j * Np_x + Np_x - 1;
            flux1 *= alpha_val / (ir->IntPoint(Np_x - 1).weight);
            elvect_mat.SetRow(id1, flux1);
        }
    }

    if (dim > 1)
    {
        for (int k = 0; k < Np_z; k++)
        {
            for (int i = 0; i < Np_x; i++)
            {
                flux1 = 0.0;
                for (int j = 0; j < Np_y - 1; j++)
                {
                    int id1 = k * Np_y * Np_x + j * Np_x + i;
                    int id2 = k * Np_y * Np_x + (j + 1) * Np_x + i;
                    vdof_mat.GetRow(id1, state1);
                    vdof_mat.GetRow(id2, state2);
                    SubcellMetricEta(Tr.ElementNo).GetColumn(id2, nor1);
                    max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor1, flux_num));

                    flux1 -= flux_num;
                    flux1 *= alpha_val / (ir->IntPoint(j).weight);
                    elvect_mat.GetRow(id1, flux2);
                    flux2 += flux1;
                    elvect_mat.SetRow(id1, flux2);

                    flux1 = flux_num;             
                }
                id1 = k * Np_y * Np_x + (Np_y - 1) * Np_x + i;
                flux1 *= alpha_val / (ir->IntPoint(Np_y - 1).weight);
                elvect_mat.GetRow(id1, flux2); // could this be imroved so that there is
                flux2 += flux1; //  no need to extract previously computed contribution?
                elvect_mat.SetRow(id1, flux2); // just directly adding the new contribution to elvect_mat?
            }
        }

        if (dim > 2)
        {
            for (int j = 0; j < Np_y; j++)
            {
                for (int i = 0; i < Np_x; i++)
                {
                    flux1 = 0.0;
                    for (int k = 0; k < Np_z - 1; k++)
                    {
                        int id1 = k * Np_y * Np_x + j * Np_x + i;
                        int id2 = (k + 1) * Np_y * Np_x + j * Np_x + i;
                        vdof_mat.GetRow(id1, state1);
                        vdof_mat.GetRow(id2, state2);
                        SubcellMetricZeta(Tr.ElementNo).GetColumn(id2, nor1);
                        max_char_speed = std::max(max_char_speed, rsolver.ComputeFaceFlux(state1, state2, nor1, flux_num));

                        flux1 -= flux_num;
                        flux1 *= alpha_val / (ir->IntPoint(k).weight);
                        elvect_mat.GetRow(id1, flux2);
                        flux2 += flux1;
                        elvect_mat.SetRow(id1, flux2);

                        flux1 = flux_num;              
                    }
                    id1 = (Np_z - 1) * Np_y * Np_x + j * Np_x + i;
                    flux1 *= alpha_val / (ir->IntPoint(Np_z - 1).weight);
                    elvect_mat.GetRow(id1, flux2);
                    flux2 += flux1;
                    elvect_mat.SetRow(id1, flux2);
                }
            }
        }
    }
}


void DGSEMIntegrator::ComputeSubcellMetrics()
{
    for (int el = 0; el < pmesh->GetNE(); el++)
    {
        ElementTransformation *Tr = pmesh->GetElementTransformation(el);
        DenseMatrix &nor_mat_xi = SubcellMetricXi(el);
        DenseMatrix &nor_mat_eta = SubcellMetricEta(el);
        DenseMatrix &nor_mat_zeta = SubcellMetricZeta(el);
        Vector nor_xi(dim), nor_eta(dim), nor_zeta(dim), tmp(dim);
        for (int k = 0; k < Np_z; k++)
        {
            for (int j = 0; j < Np_y; j++)
            {
                const IntegrationPoint &ip = ir_vol->IntPoint(k * Np_y * Np_x + j * Np_x); // left xi-face
                Tr->SetIntPoint(&ip);
                Tr->AdjugateJacobian().GetRow(0, nor1);
                for (int i = 0; i < Np_x + 1; i++)
                {
                    // nor_mat_xi.GetColumnReference(k * Np_y * Np_x + j * Np_x + i, nor_xi);
                    nor2 = nor1;
                    for (int l = 0; l < i; l++)
                    {
                        D.GetRow(l, Dx);
                        tmp = 0.0;
                        real_t weight = ir->IntPoint(l).weight;
                        for (int m = 0; m < Np_x; m++)
                        {
                            Tr->SetIntPoint(&ir_vol->IntPoint(k * Np_y * Np_x + j * Np_x + m));
                            Tr->AdjugateJacobian().GetRow(0, metric1);
                            metric1 *= Dx(m);
                            tmp += metric1;
                        }
                        tmp *= weight;
                        nor2 += tmp;
                    }
                    // Normalize(nor2);
                    nor_mat_xi.SetCol(k * Np_y * Np_x + j * Np_x + i, nor2);
                    // nor_xi = nor2;
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
                    Tr->AdjugateJacobian().GetRow(1, nor1);
                    for (int j = 0; j < Np_y + 1; j++)
                    {
                        // nor_mat_eta.GetColumnReference(k * Np_y * Np_x + j * Np_y + i, nor_eta);
                        nor2 = nor1;
                        for (int l = 0; l < j; l++)
                        {
                            D.GetRow(l, Dx);
                            tmp = 0.0;
                            real_t weight = ir->IntPoint(l).weight;
                            for (int m = 0; m < Np_y; m++)
                            {
                                Tr->SetIntPoint(&ir_vol->IntPoint(k * Np_y * Np_x + m * Np_x + i));
                                Tr->AdjugateJacobian().GetRow(1, metric1);
                                metric1 *= Dx(m);
                                tmp += metric1;
                            }
                            tmp *= weight;
                            nor2 += tmp;
                        }
                        // Normalize(nor2);
                        nor_mat_eta.SetCol(k * Np_y * Np_x + j * Np_x + i, nor2);
                        // nor_eta = nor2;
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
                        Tr->AdjugateJacobian().GetRow(2, nor1);
                        for (int k = 0; k < Np_z + 1; k++)
                        {
                            // nor_mat_zeta.GetColumnReference(k * Np_y * Np_x + j * Np_x + i, nor_zeta);
                            nor2 = nor1;
                            for (int l = 0; l < k; l++)
                            {
                                D.GetRow(l, Dx);
                                tmp = 0.0;
                                real_t weight = ir->IntPoint(l).weight;
                                for (int m = 0; m < Np_z; m++)
                                {
                                    Tr->SetIntPoint(&ir_vol->IntPoint(m * Np_y * Np_x + j * Np_x + i));
                                    Tr->AdjugateJacobian().GetRow(2, metric1);
                                    metric1 *= Dx(m);
                                    tmp += metric1;
                                }
                                tmp *= weight;
                                nor2 += tmp;
                            }
                            // Normalize(nor2);
                            nor_mat_zeta.SetCol(k * Np_y * Np_x + j * Np_x + i, nor2);
                            // nor_zeta = nor2;
                        }
                    }
                }
            }
        }
    }   
}


}