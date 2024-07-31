#pragma once

#include "mfem.hpp"

namespace Prandtl
{

using namespace mfem;

class DGOperator : public TimeDependentOperator
{
private:
   const int num_equations; // the number of equations
   const int dim;
   ParFiniteElementSpace &vfes; // vector finite element space
   ParMesh &mesh; // mesh object
   // Element integration form. Should contain ComputeFlux
   std::unique_ptr<HyperbolicFormIntegrator> formIntegrator;
   // Base Nonlinear Form
   std::unique_ptr<NonlinearForm> nonlinearForm;
   // element-wise inverse mass matrix
   std::vector<DenseMatrix> invmass; // local scalar inverse mass
   std::vector<DenseMatrix> weakdiv; // local weak divergence (trial space ByDim)
   // global maximum characteristic speed. Updated by form integrators
   mutable real_t max_char_speed;
   // auxiliary variable used in Mult
   mutable Vector z;

   // Compute element-wise inverse mass matrix
   void ComputeInvMass();
   // Compute element-wise weak-divergence matrix
   void ComputeWeakDivergence();
public:
   /**
    * @brief Construct a new DGHyperbolicConservationLaws object
    *
    * @param vfes_ vector finite element space. Only tested for DG [Pₚ]ⁿ
    * @param formIntegrator_ integrator (F(u,x), grad v)
    * @param preassembleWeakDivergence preassemble weak divergence for faster
    *                                  assembly
    */
   DGOperator(
      ParFiniteElementSpace &vfes_, ParMesh &mesh_,
      std::unique_ptr<HyperbolicFormIntegrator> formIntegrator_,
      bool preassembleWeakDivergence=true);
   /**
    * @brief Apply nonlinear form to obtain M⁻¹(DIVF + JUMP HAT(F))
    *
    * @param x current solution vector
    * @param y resulting dual vector to be used in an EXPLICIT solver
    */
   void Mult(const Vector &x, Vector &y) const override;
   // get global maximum characteristic speed to be used in CFL condition
   // where max_char_speed is updated during Mult.
   inline real_t GetMaxCharSpeed() { return max_char_speed; }
   void Update();

};

}