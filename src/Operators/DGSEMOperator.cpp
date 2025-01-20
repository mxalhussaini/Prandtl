#include "DGSEMOperator.hpp"
#include "Physics.hpp"

namespace Prandtl
{

DGSEMOperator::DGSEMOperator(std::shared_ptr<ParFiniteElementSpace> vfes_,
                             std::shared_ptr<ParFiniteElementSpace> fes0_,
                             std::shared_ptr<ParMesh> pmesh_,
                             std::shared_ptr<ParGridFunction> eta_,
                             std::shared_ptr<ParGridFunction> alpha_,
                             std::shared_ptr<ParGridFunction> grad_x_,
                             std::shared_ptr<ParGridFunction> grad_y_,
                             std::shared_ptr<ParGridFunction> grad_z_,
                             std::unique_ptr<DGSEMIntegrator> integrator_,
                             std::unique_ptr<Indicator> indicator_,
                             std::vector<BdrFaceIntegrator*> bfnfi_,
                             std::vector<Array<int>> bdr_marker_,
                             const real_t alpha_max, const real_t alpha_min)
                             : TimeDependentOperator(vfes_->GetTrueVSize()),
                             vfes(vfes_), fes0(fes0_), pmesh(pmesh_),
                             eta(eta_), alpha(alpha_), grad_x(grad_x_), grad_y(grad_y_), grad_z(grad_z_),
                             integrator(std::move(integrator_)), indicator(std::move(indicator_)), bfnfi(bfnfi_),
                             bdr_marker(bdr_marker_), num_equations(vfes->GetVDim()), dim(pmesh->SpaceDimension()),
                             order(vfes->GetElementOrder(0)),
                             num_elements(pmesh->GetNE()), Ndofs(vfes->GetFE(0)->GetDof()),
                             modalThreshold(0.5 * std::pow(10.0, -1.8 * std::pow(order, 0.25))),
                             alpha_max(alpha_max), alpha_min(alpha_min)
{
    nonlinearForm.reset(new ParNonlinearForm(vfes.get()));
#ifdef PARABOLIC
    nonlinearForm_Lifting.reset(new ParNonlinearForm(vfes.get()));
#endif
    nonlinearForm->AddDomainIntegrator(integrator.get());
    nonlinearForm->AddInteriorFaceIntegrator(integrator.get());

    std::vector<BdrFaceIntegrator*>::iterator it1 = bfnfi.begin();
    std::vector<Array<int>>::iterator it2 = bdr_marker.begin();

    for (; it1 != bfnfi.end() && it2 != bdr_marker.end(); ++it1, ++it2)
    {
        nonlinearForm->AddBdrFaceIntegrator(*it1, *it2);
#ifdef PARABOLIC
        nonlinearForm_Lifting->AddBdrFaceIntegrator(*it1, *it2);
#endif
    }
    nonlinearForm->UseExternalIntegrators();

#ifdef PARABOLIC
    nonlinearForm_Lifting->AddInteriorFaceIntegrator(integrator.get());
    nonlinearForm_Lifting->UseExternalIntegrators();
    global_entropy.SetSize(vfes->GetVSize());

    grad_mats.resize(pmesh->GetNE());
    for (int el = 0; el < num_elements; el++)
    {
        integrator->AssembleElementMatrix2(*vfes->GetFE(el), *vfes->GetFE(el),
                                           *vfes->GetElementTransformation(el),
                                           grad_mats[el]);  
    }
    grad_vdof_mats.SetSize(Ndofs * dim, num_equations);
#endif
}

void DGSEMOperator::ComputeBlendingCoefficient(const Vector &x) const
{
    indicator->CheckSmoothness(x);
    for (int el = 0; el < num_elements; el++)
    {
        fes0->GetElementDofs(el, ind_indx);
        eta->GetSubVector(ind_indx, ind_dof);
        alpha_dof = 1.0 / (1.0 + std::exp(-sharpness_fac * (ind_dof(0) - modalThreshold) / modalThreshold));
        if (alpha_dof < alpha_min)
        {
            alpha_dof = 0.0;
        }
        else if (alpha_dof > (1.0 - alpha_min))
        {
            alpha_dof = 1.0;
        }
        alpha_dof = std::min(alpha_dof, alpha_max);
        alpha->SetSubVector(ind_indx, alpha_dof);
    }
    
}

void DGSEMOperator::ComputeGlobalEntropyVector(const Vector &x, Vector &y) const
{
    DenseMatrix ent_mat(Ndofs, num_equations);
    for (int el = 0; el < num_elements; el++)
    {
        ElementTransformation *Tr = vfes->GetElementTransformation(el);
        vfes->GetElementVDofs(el, vdof_indices);
        x.GetSubVector(vdof_indices, el_vdofs);
        DenseMatrix vdof_mat(el_vdofs.GetData(), Ndofs, num_equations);
        Conserv2Entropy(vdof_mat, ent_mat);
        y.SetSubVector(vdof_indices, ent_mat.GetData());
    }
}

void DGSEMOperator::Mult(const Vector &x, Vector &y) const
{
#ifdef PARABOLIC
    ComputeGlobalEntropyVector(x, global_entropy);

    integrator->ChooseDirection(0);
    integrator->SetMode(DGSEMIntegrator::GRADIENT);

    for (int b = 0; b < bfnfi.size(); b++)
    {
        bfnfi[b]->ChooseDirection(0);
        bfnfi[b]->SetMode(BdrFaceIntegrator::GRADIENT);
    }
    nonlinearForm_Lifting->Mult(global_entropy, *grad_x);
    grad_x->ExchangeFaceNbrData();

    if (dim > 1)
    {
        integrator->ChooseDirection(1);
        for (int b = 0; b < bfnfi.size(); b++)
        {
            bfnfi[b]->ChooseDirection(1);
        }
        nonlinearForm_Lifting->Mult(global_entropy, *grad_y);
        grad_y->ExchangeFaceNbrData();

        if (dim > 2)
        {
            integrator->ChooseDirection(2);
            for (int b = 0; b < bfnfi.size(); b++)
            {
                bfnfi[b]->ChooseDirection(2);
            }
            nonlinearForm_Lifting->Mult(global_entropy, *grad_z);
            grad_z->ExchangeFaceNbrData();
        }
    }
    
    for (int el = 0; el < num_elements; el++)
    {
        ElementTransformation *Tr = vfes->GetElementTransformation(el);
        vfes->GetElementVDofs(el, vdof_indices);
        global_entropy.GetSubVector(vdof_indices, ent_vdofs);
        DenseMatrix ent_mat(ent_vdofs.GetData(), Ndofs, num_equations);
        x.GetSubVector(vdof_indices, el_vdofs);
        DenseMatrix vdof_mat(el_vdofs.GetData(), Ndofs, num_equations);

        mfem::Mult(grad_mats[el], ent_mat, grad_vdof_mats);

        grad_x->GetSubVector(vdof_indices, grad_vdofs);
        grad_mat1.UseExternalData(grad_vdofs.GetData(), Ndofs, num_equations);
        grad_vdof_mats.GetSubMatrix(0, Ndofs, 0, num_equations, grad_mat2);
        // grad_x->AddElementVector(vdof_indices, grad_mat2.GetData());
        grad_mat1 += grad_mat2;
        EntropyGrad2PrimGrad(vdof_mat, grad_mat1);
        grad_x->SetSubVector(vdof_indices, grad_mat1.GetData());

        if (dim > 1)
        {
            grad_y->GetSubVector(vdof_indices, grad_vdofs);
            grad_mat1.UseExternalData(grad_vdofs.GetData(), Ndofs, num_equations);
            grad_vdof_mats.GetSubMatrix(Ndofs, 2 * Ndofs, 0, num_equations, grad_mat2);
            grad_mat1 += grad_mat2;
            EntropyGrad2PrimGrad(vdof_mat, grad_mat1);
            grad_y->SetSubVector(vdof_indices, grad_mat1.GetData());

            if (dim > 2)
            {
                grad_z->GetSubVector(vdof_indices, grad_vdofs);
                grad_mat1.UseExternalData(grad_vdofs.GetData(), Ndofs, num_equations);
                grad_vdof_mats.GetSubMatrix(2 * Ndofs, 3 * Ndofs, 0, num_equations, grad_mat2);
                grad_mat1 += grad_mat2;
                EntropyGrad2PrimGrad(vdof_mat, grad_mat1);
                grad_z->SetSubVector(vdof_indices, grad_mat1.GetData());
            }
        }
        // vdof_mat.ClearExternalData();
        grad_mat1.ClearExternalData();
        // ent_mat.ClearExternalData();
    }
#endif
    integrator->SetMode(DGSEMIntegrator::DIVERGENCE);
    for (int b = 0; b < bfnfi.size(); b++)
    {
        bfnfi[b]->SetMode(BdrFaceIntegrator::DIVERGENCE);
    }

    ComputeBlendingCoefficient(x);
    nonlinearForm->Mult(x, y);

    max_char_speed = integrator->GetMaxCharSpeed();
    for (int b = 0; b < bfnfi.size(); b++)
    {
        max_char_speed = std::max(bfnfi[b]->GetMaxCharSpeed(), max_char_speed);
    }
}

void DGSEMOperator::AddBdrFaceIntegrator(BdrFaceIntegrator *bfi, Array<int> &bdr_marker)
{
    nonlinearForm->AddBdrFaceIntegrator(bfi, bdr_marker);
#ifdef PARABOLIC
    nonlinearForm_Lifting->AddBdrFaceIntegrator(bfi, bdr_marker);
    bfnfi.push_back(bfi);
#endif
}

}