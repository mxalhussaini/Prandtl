#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{

using namespace mfem;

class RiemannInvariantBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    VectorFunctionCoefficient prim_state_fun;
    Vector prim_state;
    real_t rho_o, p_o, p_i, p_b, a_o, a_i, a_b;
    Vector V_o;
    real_t Vn_o, Vn_i, Vn_b, dVn_o, dVn_i;
    Vector unit_nor;
    real_t R_plus, R_minus; // Riemann Invariants

    const real_t gammaInverse;

public:
    RiemannInvariantBdrFaceIntegrator(std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np, const real_t &time, real_t gamma, VectorFunctionCoefficient &prim_state_fun, bool t_dependent = false);
    RiemannInvariantBdrFaceIntegrator(std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np, const real_t &time, real_t gamma, const Vector &prim_state);

    virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, const Vector &dqdy, const Vector &dqdz, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, const Vector &dqdy, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};

}