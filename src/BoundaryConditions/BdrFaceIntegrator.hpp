#pragma once

#include "mfem.hpp"

namespace Prandtl

{

using namespace mfem;

/**
 * @brief Abstract boundary face integrator class
 *
 */
class BdrFaceIntegrator : public NonlinearFormIntegrator
{
private:
   // The maximum characteristic speed, updated during element/face vector assembly
   real_t max_char_speed;
   const RiemannSolver &rsolver;   // Numerical flux that maps F(u±,x) to hat(F)
   const FluxFunction &fluxFunction;
   // const int IntOrderOffset; // integration order offset, 2*p + IntOrderOffset.
#ifndef MFEM_THREAD_SAFE
   Vector shape1;  // shape function value at an integration point - first elem
   Vector state1;  // state value at an integration point - first elem
   Vector state2;  // state value at an integration point - boundary state
   Vector phys_ip; // integration point in physical space
   Vector fluxN;   // hat(F)(u,x)
#endif

protected:
   Vector nor;     // normal vector, @see CalcOrtho

public:
   const int num_equations;  // the number of equations
   /**
    * @brief Construct a new Hyperbolic Form Integrator object
    *
    * @param[in] rsolver numerical flux
    * @param[in] IntOrderOffset integration order offset
    */
   BdrFaceIntegrator(const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir);

   /**
    * @brief Reset the Max Char Speed 0
    *
    */
   void ResetMaxCharSpeed()
   {
      max_char_speed = 0.0;
   }

   real_t GetMaxCharSpeed()
   {
      return max_char_speed;
   }

   const FluxFunction &GetFluxFunction() { return fluxFunction; }

   virtual void ComputeOuterState(const Vector& state1,
      Vector& state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) = 0;
   virtual real_t ComputeBdrFaceFlux(const Vector &state1,
      Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip);

   /**
    * @brief implement <-hat(F)(u,x) n, [[v]]> with abstract hat(F) computed by
    * ComputeFluxDotN and numerical flux object
    *
    * @param[in] el1 finite element of the first element
    * @param[in] el2 finite element of the second element
    * @param[in] Tr face element transformations
    * @param[in] elfun local coefficient of basis from both elements
    * @param[out] elvect evaluated dual vector <-hat(F)(u,x) n, [[v]]>
    */
   void AssembleFaceVector(const FiniteElement &el1,
                           const FiniteElement &el2,
                           FaceElementTransformations &Tr,
                           const Vector &elfun, Vector &elvect) override;

};


}