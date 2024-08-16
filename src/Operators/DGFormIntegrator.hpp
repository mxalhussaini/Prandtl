#pragma once

#include "mfem.hpp"

namespace Prandtl
{

using namespace mfem;

/**
 * @brief Abstract hyperbolic form integrator, (F(u, x), ∇v) and (F̂(u±, x, n))
 *
 */
class DGFormIntegrator : public NonlinearFormIntegrator
{
private:
   // The maximum characteristic speed, updated during element/face vector assembly
   real_t max_char_speed;
   const RiemannSolver &rsolver;   // Numerical flux that maps F(u±,x) to hat(F)
   const FluxFunction &fluxFunction;
   const IntegrationRule *vol_ir, *face_ir;
#ifndef MFEM_THREAD_SAFE
   // Local storage for element integration
   Vector shape;              // shape function value at an integration point
   Vector state;              // state value at an integration point
   DenseMatrix flux;          // flux value at an integration point
   DenseMatrix dshape;  // derivative of shape function at an integration point

   Vector shape1;  // shape function value at an integration point - first elem
   Vector shape2;  // shape function value at an integration point - second elem
   Vector state1;  // state value at an integration point - first elem
   Vector state2;  // state value at an integration point - second elem
   Vector nor;     // normal vector, @see CalcOrtho
   Vector fluxN;   // hat(F)(u,x)
#endif

public:
   const int num_equations;  // the number of equations
   /**
    * @brief Construct a new Hyperbolic Form Integrator object
    *
    * @param[in] rsolver numerical flux
    * @param[in] IntOrderOffset integration order offset
    */
   DGFormIntegrator(const RiemannSolver &rsolver, const IntegrationRule *vol_ir_,
    const IntegrationRule *face_ir_);

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

   inline const IntegrationRule* GetVolumeIr()
   {
      return vol_ir;
   }

   inline const IntegrationRule* GetFaceIr()
   {
      return face_ir;
   }

   const FluxFunction &GetFluxFunction() { return fluxFunction; }

   /**
    * @brief implement (F(u), grad v) with abstract F computed by ComputeFlux
    *
    * @param[in] el local finite element
    * @param[in] Tr element transformation
    * @param[in] elfun local coefficient of basis
    * @param[out] elvect evaluated dual vector
    */
   void AssembleElementVector(const FiniteElement &el,
                              ElementTransformation &Tr,
                              const Vector &elfun, Vector &elvect) override;

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