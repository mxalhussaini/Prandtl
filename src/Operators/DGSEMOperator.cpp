#include "DGSEMOperator.hpp"

namespace Prandtl
{

DGSEMOperator::DGSEMOperator(std::shared_ptr<ParFiniteElementSpace> vfes_,
                             std::shared_ptr<ParFiniteElementSpace> fes0_,
                             std::shared_ptr<ParMesh> pmesh_,
                             std::shared_ptr<ParGridFunction> eta_,
                             std::shared_ptr<ParGridFunction> alpha_,
                             std::shared_ptr<ParGridFunction> dudx_,
                             std::shared_ptr<ParGridFunction> dudy_,
                             std::shared_ptr<ParGridFunction> dudz_,
                             std::unique_ptr<DGSEMIntegrator> integrator_,
                             std::unique_ptr<Indicator> indicator_,
                             const real_t gamma_,
                             const real_t alpha_max, const real_t alpha_min)
                             : TimeDependentOperator(vfes_->GetTrueVSize()),
                             vfes(vfes_), fes0(fes0_), pmesh(pmesh_),
                             eta(eta_), alpha(alpha_), dudx(dudx_), dudy(dudy_), dudz(dudz_),
                             integrator(std::move(integrator_)), indicator(std::move(indicator_)),
                             num_equations(vfes->GetVDim()), dim(pmesh->SpaceDimension()),
                             order(vfes->GetElementOrder(0)), num_elements(pmesh->GetNE()),
                             Ndofs(vfes->GetFE(0)->GetDof()),
                             modalThreshold(0.5 * std::pow(10.0, -1.8 * std::pow(order, 0.25))),
                             gamma(gamma_), gammaM1(gamma - 1.0), gammaM1Inverse(1.0 / gammaM1),
                             alpha_max(alpha_max), alpha_min(alpha_min)
{
    nonlinearForm.reset(new DGSEMNonlinearForm(vfes.get()));

    nonlinearForm->AddDomainIntegrator(integrator.get());
    nonlinearForm->AddInteriorFaceIntegrator(integrator.get());

    std::vector<BdrFaceIntegrator*>::iterator it1 = bfnfi.begin();
    std::vector<Array<int>>::iterator it2 = bdr_marker.begin();

    for (; it1 != bfnfi.end() && it2 != bdr_marker.end(); ++it1, ++it2)
    {
        nonlinearForm->AddBdrFaceIntegrator(*it1, *it2);
    }
    nonlinearForm->UseExternalIntegrators();

#ifdef PARABOLIC
    global_entropy.SetSize(vfes->GetVSize());
#endif
}

DGSEMOperator::~DGSEMOperator()
{
    for (auto ptr : bfnfi)
    {
        delete ptr;
    }
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

void DGSEMOperator::ComputeGlobalEntropyVector(const Vector &u, Vector &global_entropy) const
{
    DenseMatrix ent_mat(Ndofs, num_equations);
    for (int el = 0; el < num_elements; el++)
    {
        ElementTransformation *Tr = vfes->GetElementTransformation(el);
        vfes->GetElementVDofs(el, vdof_indices);
        u.GetSubVector(vdof_indices, el_vdofs);
        DenseMatrix vdof_mat(el_vdofs.GetData(), Ndofs, num_equations);
        Conserv2Entropy(vdof_mat, ent_mat, gamma, gammaM1, gammaM1Inverse);
        global_entropy.SetSubVector(vdof_indices, ent_mat.GetData());
    }
}

void DGSEMOperator::ComputeGlobalPrimitiveGradVector(const Vector &u, Vector &dudx) const
{
    for (int el = 0; el < num_elements; el++)
    {
        ElementTransformation *Tr = vfes->GetElementTransformation(el);
        vfes->GetElementVDofs(el, vdof_indices);
        u.GetSubVector(vdof_indices, el_vdofs);
        DenseMatrix vdof_mat(el_vdofs.GetData(), Ndofs, num_equations);

        dudx.GetSubVector(vdof_indices, grad_vdofs);
        DenseMatrix grad_mat(grad_vdofs.GetData(), Ndofs, num_equations);
        EntropyGrad2PrimGrad(vdof_mat, grad_mat, gammaM1, gammaM1Inverse);
        dudx.SetSubVector(vdof_indices, grad_mat.GetData());
    }
}

void DGSEMOperator::ComputeGlobalPrimitiveGradVector(const Vector &u, Vector &dudx, Vector &dudy) const
{
    for (int el = 0; el < num_elements; el++)
    {
        ElementTransformation *Tr = vfes->GetElementTransformation(el);
        vfes->GetElementVDofs(el, vdof_indices);
        u.GetSubVector(vdof_indices, el_vdofs);
        DenseMatrix vdof_mat(el_vdofs.GetData(), Ndofs, num_equations);

        dudx.GetSubVector(vdof_indices, grad_vdofs);
        DenseMatrix grad_mat1(grad_vdofs.GetData(), Ndofs, num_equations);
        EntropyGrad2PrimGrad(vdof_mat, grad_mat1, gammaM1, gammaM1Inverse);
        dudx.SetSubVector(vdof_indices, grad_mat1.GetData());

        dudy.GetSubVector(vdof_indices, grad_vdofs);
        DenseMatrix grad_mat2(grad_vdofs.GetData(), Ndofs, num_equations);
        EntropyGrad2PrimGrad(vdof_mat, grad_mat2, gammaM1, gammaM1Inverse);
        dudy.SetSubVector(vdof_indices, grad_mat2.GetData());    
        
    }
}

void DGSEMOperator::ComputeGlobalPrimitiveGradVector(const Vector &u, Vector &dudx, Vector &dudy, Vector &dudz) const
{
    for (int el = 0; el < num_elements; el++)
    {
        ElementTransformation *Tr = vfes->GetElementTransformation(el);
        vfes->GetElementVDofs(el, vdof_indices);
        u.GetSubVector(vdof_indices, el_vdofs);
        DenseMatrix vdof_mat(el_vdofs.GetData(), Ndofs, num_equations);

        dudx.GetSubVector(vdof_indices, grad_vdofs);
        DenseMatrix grad_mat1(grad_vdofs.GetData(), Ndofs, num_equations);
        EntropyGrad2PrimGrad(vdof_mat, grad_mat1, gammaM1, gammaM1Inverse);
        dudx.SetSubVector(vdof_indices, grad_mat1.GetData());

        dudy.GetSubVector(vdof_indices, grad_vdofs);
        DenseMatrix grad_mat2(grad_vdofs.GetData(), Ndofs, num_equations);
        EntropyGrad2PrimGrad(vdof_mat, grad_mat2, gammaM1, gammaM1Inverse);
        dudy.SetSubVector(vdof_indices, grad_mat2.GetData());

        dudz.GetSubVector(vdof_indices, grad_vdofs);
        DenseMatrix grad_mat3(grad_vdofs.GetData(), Ndofs, num_equations);
        EntropyGrad2PrimGrad(vdof_mat, grad_mat3, gammaM1, gammaM1Inverse);
        dudz.SetSubVector(vdof_indices, grad_mat3.GetData());      
    
    }
}

void DGSEMOperator::Mult(const Vector &u, Vector &dudt) const
{
#ifdef SUBCELL_FV_BLENDING
    ComputeBlendingCoefficient(u);
#endif

#ifdef PARABOLIC
    ComputeGlobalEntropyVector(u, global_entropy);

    if (dim == 1)
    {
        nonlinearForm->MultLifting(global_entropy, *dudx);
        ComputeGlobalPrimitiveGradVector(u, *dudx);
        nonlinearForm->Mult(u, *dudx, dudt);
    }
    else if (dim == 2)
    {
        nonlinearForm->MultLifting(global_entropy, *dudx, *dudy);
        ComputeGlobalPrimitiveGradVector(u, *dudx, *dudy);
        nonlinearForm->Mult(u, *dudx, *dudy, dudt);
    }
    else
    {
        ComputeGlobalPrimitiveGradVector(u, *dudx, *dudy, *dudz);
        nonlinearForm->Mult(u, *dudx, *dudy, *dudz, dudt);
        nonlinearForm->MultLifting(global_entropy, *dudx, *dudy, *dudz);
    }
#else
    nonlinearForm->Mult(u, dudt);
    max_char_speed = integrator->GetMaxCharSpeed();
    for (int b = 0; b < bfnfi.size(); b++)
    {
        max_char_speed = std::max(bfnfi[b]->GetMaxCharSpeed(), max_char_speed);
    }
#endif
}

void DGSEMOperator::AddBdrFaceIntegrator(BdrFaceIntegrator *bfi, Array<int> &bdr_marker)
{
    nonlinearForm->AddBdrFaceIntegrator(bfi, bdr_marker);
}

}