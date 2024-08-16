#include "DGOperator.hpp"

namespace Prandtl
{

// Implementation of class DGOperator
DGOperator::DGOperator(
   std::shared_ptr<ParFiniteElementSpace> vfes_,
   std::shared_ptr<ParFiniteElementSpace> fes0_,
   std::shared_ptr<ParMesh> mesh_,
   std::shared_ptr<ParGridFunction> sol_,
   std::unique_ptr<DGFormIntegrator> formIntegrator_,
   std::unique_ptr<Filter> filter_,
   bool preassembleWeakDivergence)
   :  TimeDependentOperator(vfes_->GetTrueVSize()),
      vfes(vfes_), fes0(fes0_), mesh(mesh_), sol(sol_),
      formIntegrator(std::move(formIntegrator_)), filter(std::move(filter_)),
      num_equations(formIntegrator_->num_equations),
      dim(mesh_->SpaceDimension()),
      z(vfes_->GetTrueVSize())
{
   // Standard local assembly and inversion for energy mass matrices.
   ComputeInvMass();
#ifndef MFEM_USE_MPI
   nonlinearForm.reset(new NonlinearForm(vfes));
#else
   nonlinearForm.reset(new ParNonlinearForm(vfes.get()));
#endif
   if (preassembleWeakDivergence)
   {
      ComputeWeakDivergence();
   }
   else
   {
      nonlinearForm->AddDomainIntegrator(formIntegrator.get());
   }
   nonlinearForm->AddInteriorFaceIntegrator(formIntegrator.get());
   nonlinearForm->UseExternalIntegrators();
}

void DGOperator::ComputeWeakDivergence()
{
   TransposeIntegrator weak_div(new GradientIntegrator());
   DenseMatrix weakdiv_bynodes;

   weak_div.SetIntRule(formIntegrator->GetVolumeIr());

   weakdiv.resize(vfes->GetNE());
   for (int i=0; i<vfes->GetNE(); i++)
   {
      int dof = vfes->GetFE(i)->GetDof();
      weakdiv_bynodes.SetSize(dof, dof*dim);
      weak_div.AssembleElementMatrix2(*vfes->GetFE(i), *vfes->GetFE(i),
                                      *vfes->GetElementTransformation(i),
                                      weakdiv_bynodes);
      weakdiv[i].SetSize(dof, dof*dim);
      // Reorder so that trial space is ByDim.
      // This makes applying weak divergence to flux value simpler.
      for (int j=0; j<dof; j++)
      {
         for (int d=0; d<dim; d++)
         {
            weakdiv[i].SetCol(j*dim + d, weakdiv_bynodes.GetColumn(d*dof + j));
         }
      }

   }
}

void DGOperator::ComputeInvMass()
{
   InverseIntegrator inv_mass(new MassIntegrator());

   inv_mass.SetIntRule(formIntegrator->GetVolumeIr());

   invmass.resize(vfes->GetNE());
   for (int i=0; i<vfes->GetNE(); i++)
   {
      int dof = vfes->GetFE(i)->GetDof();
      invmass[i].SetSize(dof);
      inv_mass.AssembleElementMatrix(*vfes->GetFE(i),
                                     *vfes->GetElementTransformation(i),
                                     invmass[i]);
   }
}


void DGOperator::Mult(const Vector &x, Vector &y) const
{
   // 0. Reset wavespeed computation before operator application.
   formIntegrator->ResetMaxCharSpeed();
   // 1. Apply Nonlinear form to obtain an auxiliary result
   //         z = - <F̂(u_h,n), [[v]]>_e
   //    If weak-divergence is not preassembled, we also have weak-divergence
   //         z = - <F̂(u_h,n), [[v]]>_e + (F(u_h), ∇v)

   nonlinearForm->Mult(x, z);
   if (!weakdiv.empty()) // if weak divergence is pre-assembled
   {
      // Apply weak divergence to F(u_h), and inverse mass to z_loc + weakdiv_loc
      Vector current_state; // view of current state at a node
      DenseMatrix current_flux; // flux of current state
      DenseMatrix flux; // element flux value. Whose column is ordered by dim.
      DenseMatrix current_xmat; // view of current states in an element, dof x num_eq
      DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
      DenseMatrix current_ymat; // view of element result, dof x num_eq
      const FluxFunction &fluxFunction = formIntegrator->GetFluxFunction();
      Array<int> vdofs;
      Vector xval, zval;
      for (int i=0; i<vfes->GetNE(); i++)
      {
         ElementTransformation* Tr = vfes->GetElementTransformation(i);
         int dof = vfes->GetFE(i)->GetDof();
         vfes->GetElementVDofs(i, vdofs);
         x.GetSubVector(vdofs, xval);
         current_xmat.UseExternalData(xval.GetData(), dof, num_equations);
         flux.SetSize(num_equations, dim*dof);
         for (int j=0; j<dof; j++) // compute flux for all nodes in the element
         {
            current_xmat.GetRow(j, current_state);
            current_flux.UseExternalData(flux.GetData() + num_equations*dim*j,
                                         num_equations, dof);
            fluxFunction.ComputeFlux(current_state, *Tr, current_flux);
         }
         // Compute weak-divergence and add it to auxiliary result, z
         // Recalling that weakdiv is reordered by dim, we can apply
         // weak-divergence to the transpose of flux.
         z.GetSubVector(vdofs, zval);
         current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
         mfem::AddMult_a_ABt(1.0, weakdiv[i], flux, current_zmat);
         // Apply inverse mass to auxiliary result to obtain the final result
         current_ymat.SetSize(dof, num_equations);
         mfem::Mult(invmass[i], current_zmat, current_ymat);
         y.SetSubVector(vdofs, current_ymat.GetData());
      }
   }
   else
   {
      // Apply block inverse mass
      Vector zval; // z_loc, dof*num_eq

      DenseMatrix current_zmat; // view of element auxiliary result, dof x num_eq
      DenseMatrix current_ymat; // view of element result, dof x num_eq
      Array<int> vdofs;
      for (int i=0; i<vfes->GetNE(); i++)
      {
         int dof = vfes->GetFE(i)->GetDof();
         vfes->GetElementVDofs(i, vdofs);
         z.GetSubVector(vdofs, zval);
         current_zmat.UseExternalData(zval.GetData(), dof, num_equations);
         current_ymat.SetSize(dof, num_equations);
         mfem::Mult(invmass[i], current_zmat, current_ymat);
         y.SetSubVector(vdofs, current_ymat.GetData());
      }
   }
   max_char_speed = formIntegrator->GetMaxCharSpeed();
}

void DGOperator::Update()
{
   nonlinearForm->Update();
   height = nonlinearForm->Height();
   width = height;
   z.SetSize(height);

   ComputeInvMass();
   if (!weakdiv.empty()) {ComputeWeakDivergence();}
}

void DGOperator::AddBdrFaceIntegrator(BdrFaceIntegrator *bfi, Array<int> &bdr_marker)
{
   nonlinearForm->AddBdrFaceIntegrator(bfi, bdr_marker);
   filter->AddBdrFaceIntegrator(bfi, bdr_marker);
}

// void DGOperator::AddBdrFaceIntegrators(Array<BdrFaceIntegrator*> bfnfi, Array<Array<int>> bfnfi_marker)
// {
//    for (int i = 0; i < bfnfi.Size(); i++)
//    {
//       nonlinearForm->AddBdrFaceIntegrator(bfnfi[i], bfnfi_marker[i]);
//    }
//    filter
// }

}

