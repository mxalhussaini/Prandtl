#include "BdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

BdrFaceIntegrator::BdrFaceIntegrator(const NumericalFlux &rsolver, int Np,
      std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
      std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
      const real_t &time, bool constant, bool t_dependent)
   : NonlinearFormIntegrator(),
     grad_x(grad_x), grad_y(grad_y), grad_z(grad_z), vfes(vfes),
     rsolver(rsolver), fluxFunction(rsolver.GetFluxFunction()),
     Np_x(Np), GLIntRules(0, Quadrature1D::GaussLobatto), 
     num_equations(fluxFunction.num_equations), dim(fluxFunction.dim),
     time(time), constant(constant), t_dependent(t_dependent)
{
   ir = &GLIntRules.Get(Geometry::SEGMENT, 2 * Np_x - 3);
   if (dim == 1)
   {
      ir_face = &GLIntRules.Get(Geometry::POINT, 2 * Np_x - 3);
   }
   else if (dim == 2)
   {
      ir_face = &GLIntRules.Get(Geometry::SEGMENT, 2 * Np_x - 3);
   }
   else
   {
      ir_face = &GLIntRules.Get(Geometry::SQUARE, 2 * Np_x - 3);
   }

   state1.SetSize(num_equations);
   state2.SetSize(num_equations);
   flux_num.SetSize(num_equations);
   flux1.SetSize(num_equations);
   flux2.SetSize(num_equations);
   flux_mat.SetSize(num_equations, dim);
   nor.SetSize(dim);
   grad_mat1.SetSize(num_equations, dim);
   grad_mat2.SetSize(num_equations, dim);
   grad_state.SetSize(num_equations);
}



void BdrFaceIntegrator::AssembleFaceVector( const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &elfun, Vector &elvect)
{
   const int dof1 = el1.GetDof();

   shape1.SetSize(dof1);

   elvect.SetSize(dof1 * num_equations);
   elvect = 0.0;

   const DenseMatrix elfun1_mat(elfun.GetData(), dof1, num_equations);

   DenseMatrix elvect1_mat(elvect.GetData(), dof1, num_equations);

#ifdef PARABOLIC
   if (currentMode == DIVERGENCE)
   {
      vfes->GetElementVDofs(Tr.Elem1No, vdof_indices);
      grad_x->GetSubVector(vdof_indices, grad_x_vdofs);
      grad_x_mat.UseExternalData(grad_x_vdofs.GetData(), dof1, num_equations);
      if (dim > 1)
      {
         grad_y->GetSubVector(vdof_indices, grad_y_vdofs);
         grad_y_mat.UseExternalData(grad_y_vdofs.GetData(), dof1, num_equations);
         if (dim > 2)
         {
            grad_z->GetSubVector(vdof_indices, grad_z_vdofs);
            grad_z_mat.UseExternalData(grad_z_vdofs.GetData(), dof1, num_equations);
         }
      }
   }
#endif
   
   for (int i = 0; i < ir_face->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir_face->IntPoint(i);

      Tr.SetAllIntPoints(&ip); // set face and element int. points
      J1 = Tr.GetElement1Transformation().Weight();
      el1.CalcShape(Tr.GetElement1IntPoint(), shape1);

      elfun1_mat.MultTranspose(shape1, state1);

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
         speed = ComputeBdrFaceInviscidFlux(state1, state2, flux1, nor, Tr, ip);
         flux1.Neg();
#ifdef PARABOLIC
         grad_x_mat.MultTranspose(shape1, grad_state);
         grad_mat1.SetCol(0, grad_state);
         if (dim > 1)
         {
            grad_y_mat.MultTranspose(shape1, grad_state);
            grad_mat1.SetCol(1, grad_state);
            if (dim > 2)
            {
               grad_z_mat.MultTranspose(shape1, grad_state);
               grad_mat1.SetCol(2, grad_state);
            }
         }

         ComputeBdrFaceViscousFlux(state1, state2, grad_mat1, flux_num, nor, Tr, ip);
         // add viscous eigenvalue estimates
         flux1 += flux_num;
#endif
      }
      else
      {
         ComputeBdrFaceLiftingFlux(state1, state2, flux1, Tr, ip);
         flux1 *= nor(dir);
      }
      // Update the global max char speed
      max_char_speed = std::max(speed, max_char_speed);

      // pre-multiply integration weight to flux
      AddMult_a_VWt(1.0 / (ir->IntPoint(0).weight * J1), shape1, flux1, elvect1_mat);
   }
}

real_t BdrFaceIntegrator::ComputeBdrFaceInviscidFlux(const Vector &state1, Vector &state2, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
   ComputeOuterInviscidState(state1, state2, Tr, ip);

   // Compute F(u+, x) and F(u-, x) with maximum characteristic speed
   // Compute hat(F) using evaluated quantities
   const real_t speed = rsolver.ComputeFaceFlux(state1, state2, nor, fluxN);
   
   return speed;
}

void BdrFaceIntegrator::ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
   Entropy2Conserv(state1, state2);
   ComputeOuterInviscidState(state2, fluxN, Tr, ip);
   Conserv2Entropy(fluxN, state2);
   fluxN = state2;
   fluxN -= state1;
   fluxN *= 0.5; 
}

}
