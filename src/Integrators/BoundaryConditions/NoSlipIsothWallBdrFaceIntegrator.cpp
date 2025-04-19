#include "NoSlipIsothWallBdrFaceIntegrator.hpp"

namespace Prandtl
{

NoSlipIsothWallBdrFaceIntegrator::NoSlipIsothWallBdrFaceIntegrator(
    std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np,
    const real_t &time, real_t gamma, real_t R_gas_, FunctionCoefficient &T_wall_, VectorFunctionCoefficient &V_wall_,
    bool t_dependent)
    : SlipWallBdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, false, t_dependent),
    T_wall(T_wall_), V_wall(V_wall_), R_gas(R_gas_) {}

NoSlipIsothWallBdrFaceIntegrator::NoSlipIsothWallBdrFaceIntegrator(
    std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np,
    const real_t &time, real_t gamma, real_t R_gas_, real_t &T, const Vector &V)
    : SlipWallBdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, true, false),
    T(T), V(V), T_wall(std::function<real_t(const Vector&)>()), V_wall(dim, std::function<void(const Vector&, Vector&)>()), R_gas(R_gas_)
{
    v = 1.0 / (R_gas * T);
}

void NoSlipIsothWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, const Vector &dqdy, const Vector &dqdz, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, dqdz, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

void NoSlipIsothWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, const Vector &dqdy, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

void NoSlipIsothWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    fluxFunction.ComputeViscousFlux(state1, dqdx, flux_mat);
    flux_mat.Mult(nor, fluxN);
}


void NoSlipIsothWallBdrFaceIntegrator::ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        T_wall.SetTime(time);
        V_wall.SetTime(time);

        T = T_wall.Eval(Tr, ip);
        V_wall.Eval(V, Tr, ip);

        v = 1.0 / (R_gas * T);
    }
    fluxN(0) = state1(0);
    fluxN(1) = V(0) * v;
    if (dim > 1)
    {
        fluxN(2) = V(1) * v;
        if (dim > 2)
        {
            fluxN(3) = V(2) * v;
        }
    }
    fluxN(num_equations - 1) = -v;
    fluxN -= state1;
}

}