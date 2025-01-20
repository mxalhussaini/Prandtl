#pragma once

#include "mfem.hpp"
#include "NumericalFlux.hpp"
#include "NavierStokesFlux.hpp"
#include "Physics.hpp"

namespace Prandtl

{

using namespace mfem;

/**
 * @brief Abstract boundary face integrator class
 *
 */
class BdrFaceIntegrator : public NonlinearFormIntegrator
{
public:
   enum Mode
   {
      GRADIENT,
      DIVERGENCE
   };
private:
   // The maximum characteristic speed, updated during element/face vector assembly
   real_t max_char_speed;
   const NumericalFlux &rsolver;   // Numerical flux that maps F(u±,x) to hat(F)
   Array<int> vdof_indices;
   Vector grad_state, grad_x_vdofs, grad_y_vdofs, grad_z_vdofs;
   DenseMatrix grad_mat1, grad_x_mat, grad_y_mat, grad_z_mat;
   int dir;

   Vector shape1;  // shape function value at an integration point - first elem
   Vector state1;  // state value at an integration point - first elem
   Vector state2;  // state value at an integration point - boundary state
   Vector flux_num, flux1, flux2;   // hat(F)(u,x)

   std::shared_ptr<ParGridFunction> grad_x, grad_y, grad_z;
   std::shared_ptr<ParFiniteElementSpace> vfes;

   const int Np_x;
   IntegrationRules GLIntRules;
   const IntegrationRule *ir, *ir_face;

protected:
   Mode currentMode;
   const NavierStokesFlux fluxFunction;
   Vector nor;     // normal vector, @see CalcOrtho
   const int num_equations, dim;
   real_t J1, speed;
   DenseMatrix flux_mat, grad_mat2;
   const real_t &time;
   bool constant, t_dependent;

public:

   BdrFaceIntegrator(const NumericalFlux &rsolver, int Np,
      std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
      std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
      const real_t &time, bool constant, bool t_dependent);
   
   void SetMode(Mode op)
   {
      currentMode = op;
   }
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

   void ChooseDirection(int dir_)
   {
      dir = dir_;
   }

   const NavierStokesFlux &GetFluxFunction() { return fluxFunction; }

   virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) {};
   virtual real_t ComputeBdrFaceInviscidFlux(const Vector &state1, Vector &state2, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip);
   virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) = 0;
   virtual void ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip);

   void AssembleFaceVector(const FiniteElement &el1,const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &elfun, Vector &elvect) override;

};


}