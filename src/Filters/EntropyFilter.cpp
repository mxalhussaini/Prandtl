#include "EntropyFilter.hpp"
#include "BasicOperations.hpp"
#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{

void EntropyFilter::ComputeMinEntropy(const Vector &state, real_t el_min_e)
{
    rho = state(0);
    p = ComputePressure(state, gamma);
    e = ComputeEntropy(rho, p, gamma);
    el_min_e = std::min(el_min_e, e);
}

void EntropyFilter::GetEntropyConstraints(const Vector &x)
{
    sol = x;
    sol.ExchangeFaceNbrData();

    Array<int> el2el, el2bdrel, vdof_indices, ent_indices;
    int indx, shifted_indx;

    for (int el = 0; el < vfes.GetNE(); el++)
    {
        element2element.GetRow(el, el2el);
        element2bdrelement.GetRow(el, el2bdrel);
        el2el.Prepend(el);
        fes0.GetElementVDofs(el, ent_indices);

        el_min_e = __DBL_MAX__;
        for (int i = 0; i < el2el.Size(); i++)
        {
            indx = el2el[i];
            const FiniteElement* fe2 = vfes.GetFE(indx);

            if (indx >= num_elements)
            {
                shifted_indx = indx - num_elements;
                vfes.GetFaceNbrElementVDofs(shifted_indx, vdof_indices);
                sol.FaceNbrData().GetSubVector(vdof_indices, el_vdofs);
                // Tr = vfes.GetFaceNbrElementTransformation(shifted_indx);
            }
            else
            {
                vfes.GetElementVDofs(indx, vdof_indices);
                sol.GetSubVector(vdof_indices, el_vdofs);
                // Tr = vfes.GetElementTransformation(indx);
            }

            vdof_mat2.UseExternalData(el_vdofs.GetData(), ndofs, num_equations);

            if (i == 0)
            {
                fe1 = fe2;
                vdof_mat1 = vdof_mat2;
            }

            for (int f = 0; f < Geometry::NumBdrArray[egeom]; f++)
            {
                pmesh.GetLocalFaceTransformation(ftype, etype, Tr_ip->Transf, 64 * f);
                Tr_ip->Transform(*face_ir, *mapped_face_ir);
                for (int j = 0; j < mapped_face_ir->GetNPoints(); j++)
                {
                    const IntegrationPoint& eip = mapped_face_ir->IntPoint(j);
                    fe2->CalcShape(eip, shape);
                    vdof_mat2.MultTranspose(shape, state2);
                    ComputeMinEntropy(state2, el_min_e);
                }
            }

            for (int k = 0; k < vol_ir->GetNPoints(); k++)
            {
                const IntegrationPoint &ip = vol_ir->IntPoint(k);
                fe2->CalcShape(ip, shape);
                vdof_mat2.MultTranspose(shape, state2);
                ComputeMinEntropy(state2, el_min_e);
            }

            for (int l = 0; l < ndofs; l++)
            {
                vdof_mat2.GetRow(l, state2);
                ComputeMinEntropy(state2, el_min_e);
            }
        }

        if (nonlinearForm.GetBdrFaceIntegrators().Size())
        {
            if (el2bdrel)
            {
                for (int m = 0; m < el2bdrel.Size(); m++)
                {
                    const int bdr_attr = pmesh.GetBdrAttribute(el2bdrel[m]);
                    if (bdr_attr_marker[bdr_attr-1] == 0) { continue; }

                    Tr_face = vfes.GetParMesh()->GetBdrFaceTransformations(el2bdrel[m]);

                    if (Tr_face != NULL)
                    {
                        for (int n = 0; n < bfnfi.Size(); n++)
                        {
                            if (bfnfi_marker[n] && (*bfnfi_marker[n])[bdr_attr-1] == 0)
                            {
                                continue;
                            }
                            for (int i = 0; i < bdr_face_ir->GetNPoints(); i++)
                            {
                                const IntegrationPoint &ip = bdr_face_ir->IntPoint(i);
                                Tr_face->SetAllIntPoints(&ip);
                                fe1->CalcShape(Tr_face->GetElement1IntPoint(), shape);
                                vdof_mat1.MultTranspose(shape, state1);

                                if (dim == 1) 
                                {
                                    nor(0) = (Tr_face->GetElement1IntPoint().x - 0.5) * 2.0;
                                }
                                else
                                {
                                    CalcOrtho(Tr_face->Jacobian(), nor);
                                }
                                Normalize(nor);
                                Tr_face->Transform(ip, phys_ip);
                                state2 = 0.0;
                                bfnfi[n]->ComputeOuterState(state1, state2, *Tr_face, ip);
                                ComputeMinEntropy(state2, el_min_e);
                            }
                        }
                    }

                }
            }
        }
        entropy.SetSubVector(ent_indices, el_min_e);
    }
}

void EntropyFilter::FilterSolution(Vector &x)
{
    std::cout << "111" << "\n";
}

}