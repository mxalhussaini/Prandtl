#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{

using namespace mfem;

class RiemannInvariantBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    FunctionCoefficient rho_out, p_out;
    VectorFunctionCoefficient V_out;
    real_t rho_o, p_o, p_i, p_b, a_o, a_i, a_b;
    Vector V_o;
    real_t Vn_o, Vn_i, Vn_b, dVn_o, dVn_i;
    Vector unit_nor;
    real_t R_plus, R_minus; // Riemann Invariants

public:
    RiemannInvariantBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                      std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                      std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                      const real_t &time, FunctionCoefficient &rho_out, FunctionCoefficient &p_out,
                                      VectorFunctionCoefficient &V_out, bool constant = true, bool t_dependent = false);
    virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};

}