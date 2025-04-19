#include "NoSlipAdiabWallBdrFaceIntegrator.hpp"

namespace Prandtl
{

NoSlipAdiabWallBdrFaceIntegrator::NoSlipAdiabWallBdrFaceIntegrator(
    std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np,
    const real_t &time, real_t gamma, FunctionCoefficient &qn_wall_, VectorFunctionCoefficient &V_wall_,
    bool t_dependent)
    : SlipWallBdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, false, t_dependent),
    qn_wall(qn_wall_), V_wall(V_wall_) {}

NoSlipAdiabWallBdrFaceIntegrator::NoSlipAdiabWallBdrFaceIntegrator(
    std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np,
    const real_t &time, real_t gamma, real_t qn, const Vector &V)
    : SlipWallBdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, true, false),
    qn(qn), V(V),
    qn_wall(std::function<real_t(const Vector&)>()), V_wall(dim, std::function<void(const Vector &, Vector&)>()) {}

void NoSlipAdiabWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, const Vector &dqdy, const Vector &dqdz, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        qn_wall.SetTime(time);
        V_wall.SetTime(time);

        qn = qn_wall.Eval(Tr, ip);
        V_wall.Eval(V, Tr, ip);
    }

    qn *= std::sqrt(nor * nor);

    fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, dqdz,  flux_mat);
    flux_mat.Mult(nor, fluxN);
    fluxN(num_equations - 1) = V(0) * fluxN(1) + V(1) * fluxN(2) + V(2) * fluxN(3) + qn;
}

void NoSlipAdiabWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, const Vector &dqdy, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        qn_wall.SetTime(time);
        V_wall.SetTime(time);

        qn = qn_wall.Eval(Tr, ip);
        V_wall.Eval(V, Tr, ip);
    }

    qn *= std::sqrt(nor * nor);

    fluxFunction.ComputeViscousFlux(state1, dqdx, dqdy, flux_mat);
    flux_mat.Mult(nor, fluxN);
    fluxN(num_equations - 1) = V(0) * fluxN(1) + V(1) * fluxN(2) + qn;
}

void NoSlipAdiabWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        qn_wall.SetTime(time);
        V_wall.SetTime(time);

        qn = qn_wall.Eval(Tr, ip);
        V_wall.Eval(V, Tr, ip);
    }

    qn *= std::sqrt(nor * nor);

    fluxFunction.ComputeViscousFlux(state1, dqdx, flux_mat);
    flux_mat.Mult(nor, fluxN);
    fluxN(num_equations - 1) = V(0) * fluxN(1) + qn;
}

void NoSlipAdiabWallBdrFaceIntegrator::ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        V_wall.SetTime(time);
        V_wall.Eval(V, Tr, ip);        
    }
    v = -state1(num_equations - 1);
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