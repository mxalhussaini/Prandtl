#include "PerssonPeraireIndicator.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{

PerssonPeraireIndicator::PerssonPeraireIndicator(std::shared_ptr<ParFiniteElementSpace> vfes,
                    std::shared_ptr<ParFiniteElementSpace> fes0,
                    std::shared_ptr<ParGridFunction> eta,
                    std::shared_ptr<ModalBasis> modalBasis)
                    : Indicator(vfes, fes0, eta), modalBasis(modalBasis),
                    ubdegs(modalBasis->GetPolyDegs())
{
    rho_p.SetSize(ndofs);
    modes.SetSize(ndofs);
    modesM1.SetSize(ndofs);
    modesM2.SetSize(ndofs);
    // ubdegs.SetSize(ndofs, vfes->GetMesh()->SpaceDimension());
    ubdegs_row.SetSize(vfes->GetMesh()->SpaceDimension());
}

void PerssonPeraireIndicator::CheckSmoothness(const Vector &x)
{
    for (int el = 0; el < vfes->GetNE(); el++)
    {
        vfes->GetElementVDofs(el, vdof_indices);
        fes0->GetElementDofs(el, ind_indx);
        x.GetSubVector(vdof_indices, el_vdofs);
        eta->GetSubVector(ind_indx, ind_dof);
        DenseMatrix vdof_mat(el_vdofs.GetData(), ndofs, num_equations);

        for (int i = 0; i < ndofs; i++)
        {
            vdof_mat.GetRow(i, state);
            rho_p(i) = state(0) * ComputePressure(state);
        }
        modalBasis->ComputeModes(rho_p);
        modalBasis->GetModes(modes);
        modesM1 = modes;
        modesM2 = modes;
        for (int i = 0; i < ndofs; i++)
        {
            ubdegs.GetRow(i, ubdegs_row);
            for (int j = 0; j < ubdegs_row.Size(); j++)
            {
                if (ubdegs_row[j] > order - 2)
                {
                    modesM2[i] = 0.0;
                    if (ubdegs_row[j] > order - 1)
                    {
                        modesM1[i] = 0.0;
                    }
                }
            }
        }
        // modes.Print(std::cout);
        // modesM1.Print(std::cout);
        // modesM2.Print(std::cout);
        // ubdegs.Print(std::cout);
        // std::cout << "----" << std::endl;
        ind_dof(0) = 1.0 - (modesM1 * modesM1) / (modes * modes);
        ind_dof(0) = std::max(ind_dof(0), 1.0 - (modesM2 * modesM2) / (modesM1 * modesM1));
        eta->SetSubVector(ind_indx, ind_dof);
    }
}

}