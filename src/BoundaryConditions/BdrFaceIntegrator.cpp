#include "BdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

BdrFaceIntegrator::BdrFaceIntegrator(
   const RiemannSolver &rsolver,
   const int IntOrderOffset)
   : NonlinearFormIntegrator(),
     rsolver(rsolver),
     fluxFunction(rsolver.GetFluxFunction()),
     IntOrderOffset(IntOrderOffset),
     num_equations(fluxFunction.num_equations)
{
#ifndef MFEM_THREAD_SAFE
   state1.SetSize(num_equations);
   state2.SetSize(num_equations);
   fluxN.SetSize(num_equations);
   nor.SetSize(fluxFunction.dim);
   phys_ip.SetSize(fluxFunction.dim);
   // unit_nor.SetSize(fluxFunction.dim);
#endif
}


void BdrFaceIntegrator::AssembleFaceVector(
   const FiniteElement &el1, const FiniteElement &el2,
   FaceElementTransformations &Tr, const Vector &elfun, Vector &elvect)
{
   // current elements' the number of degrees of freedom
   // does not consider the number of equations
   const int dof1 = el1.GetDof();

#ifdef MFEM_THREAD_SAFE
   // Local storage for element integration

   // shape function value at an integration point - first elem
   Vector shape1(dof1);
   // normal vector (usually not a unit vector)
   Vector nor(el1.GetDim());
   // unit normal vector
   // Vector unit_nor(el1.GetDim());  
   // integration point in physical space
   Vector phys_ip(el1.GetDim());  
   // state value at an integration point - first elem
   Vector state1(num_equations);
   // state value at an integration point - boundary state
   Vector state2(num_equations);
   // hat(F)(u,x)
   Vector fluxN(num_equations);
#else
   shape1.SetSize(dof1);
#endif

   elvect.SetSize(dof1 * num_equations);
   elvect = 0.0;

   const DenseMatrix elfun1_mat(elfun.GetData(), dof1, num_equations);

   DenseMatrix elvect1_mat(elvect.GetData(), dof1, num_equations);

   // Obtain integration rule. If integration is rule is given, then use it.
   // Otherwise, get (2*p + IntOrderOffset) order integration rule
   const IntegrationRule *ir = IntRule;
   if (!ir)
   {
      const int order = 2*std::max(el1.GetOrder(), el2.GetOrder()) + IntOrderOffset;
      ir = &IntRules.Get(Tr.GetGeometryType(), order);
   }
   // loop over integration points
   for (int i = 0; i < ir->GetNPoints(); i++)
   {
      const IntegrationPoint &ip = ir->IntPoint(i);

      Tr.SetAllIntPoints(&ip); // set face and element int. points

      // Calculate basis functions on both elements at the face
      el1.CalcShape(Tr.GetElement1IntPoint(), shape1);

      // Interpolate elfun at the point
      elfun1_mat.MultTranspose(shape1, state1);

      // Get the normal vector and the flux on the face
      if (nor.Size() == 1)  // if 1D, use 1 or -1.
      {
         // This assume the 1D integration point is in (0,1). This may not work
         // if this changes.
         nor(0) = (Tr.GetElement1IntPoint().x - 0.5) * 2.0;
      }
      else
      {
         CalcOrtho(Tr.Jacobian(), nor);
      }
      
      // unit_nor = nor;
      // Normalize(unit_nor);
      // Tr.Transform(ip, phys_ip);
      ComputeOuterState(state1, state2, Tr, ip);

      // Compute F(u+, x) and F(u-, x) with maximum characteristic speed
      // Compute hat(F) using evaluated quantities
      const real_t speed = rsolver.Eval(state1, state2, nor, Tr, fluxN);

      // Update the global max char speed
      max_char_speed = std::max(speed, max_char_speed);

      // pre-multiply integration weight to flux
      AddMult_a_VWt(-ip.weight, shape1, fluxN, elvect1_mat);
   }
}

}