#pragma once

#include "DGSEMIntegrator.hpp"
#include "DGSEMNonlinearForm.hpp"
#include "BdrFaceIntegrator.hpp"
#include "ModalBasis.hpp"
#include "Indicator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

using namespace mfem;

class DGSEMOperator : public TimeDependentOperator
{
private:
    std::shared_ptr<ParFiniteElementSpace> vfes;
    std::shared_ptr<ParFiniteElementSpace> fes0;
    std::shared_ptr<ParMesh> pmesh;
    std::shared_ptr<ParGridFunction> eta;
    std::shared_ptr<ParGridFunction> alpha;
    std::shared_ptr<ParGridFunction> dudx, dudy, dudz;
    std::unique_ptr<DGSEMIntegrator> integrator;
    std::unique_ptr<Indicator> indicator;
    std::unique_ptr<DGSEMNonlinearForm> nonlinearForm;

    mutable Array<int> vdof_indices;
    mutable Vector el_vdofs, grad_vdofs;

    const int Ndofs;

    mutable Vector global_entropy;
    
    const int num_equations, dim, order, num_elements;
    const real_t sharpness_fac = 9.21024;
    const real_t modalThreshold;
    const real_t alpha_min;
    const real_t alpha_max;

    mutable real_t max_char_speed;
    
    std::vector<BdrFaceIntegrator*> bfnfi;
    std::vector<Array<int>> bdr_marker;
    mutable Array<int> ind_indx;
    mutable Vector ind_dof;
    mutable real_t alpha_dof;

    const real_t gamma;
    const real_t gammaM1;
    const real_t gammaM1Inverse;
    
    void ComputeGlobalEntropyVector(const Vector &u, Vector &global_entropy) const;
    void ComputeGlobalPrimitiveGradVector(const Vector &u, Vector &dudx) const;
    void ComputeGlobalPrimitiveGradVector(const Vector &u, Vector &dudx, Vector &dudy) const;
    void ComputeGlobalPrimitiveGradVector(const Vector &u, Vector &dudx, Vector &dudy, Vector &dudz) const;
    void ComputeBlendingCoefficient(const Vector &u) const;

public:
    DGSEMOperator(std::shared_ptr<ParFiniteElementSpace> vfes,
                  std::shared_ptr<ParFiniteElementSpace> fes0,
                  std::shared_ptr<ParMesh> pmesh,
                  std::shared_ptr<ParGridFunction> eta,
                  std::shared_ptr<ParGridFunction> alpha,
                  std::shared_ptr<ParGridFunction> dudx,
                  std::shared_ptr<ParGridFunction> dudy,
                  std::shared_ptr<ParGridFunction> dudz,
                  std::unique_ptr<DGSEMIntegrator> integrator,
                  std::unique_ptr<Indicator> indicator,
                  real_t gamma,
                  const real_t alpha_max = 0.5, const real_t alpha_min = 0.001);
    
    ~DGSEMOperator();
    
    void AddBdrFaceIntegrator(BdrFaceIntegrator *bfi, Array<int> &bdr_marker);
    
    void Mult(const Vector &u, Vector &dudt) const override;
    inline real_t GetMaxCharSpeed()
    {
        return max_char_speed;
    }

    inline real_t& GetTimeRef()
    {
        return t;
    }

};

}