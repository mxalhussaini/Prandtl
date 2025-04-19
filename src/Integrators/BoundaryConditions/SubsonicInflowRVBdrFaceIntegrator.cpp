#include "SubsonicInflowRVBdrFaceIntegrator.hpp"

namespace Prandtl
{

SubsonicInflowRVBdrFaceIntegrator::SubsonicInflowRVBdrFaceIntegrator(
    std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np,
    const real_t &time, real_t gamma, FunctionCoefficient &rho_, VectorFunctionCoefficient &V_, bool t_dependent)
    : BdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, false, t_dependent),
    rho(rho_), V(V_) {}

SubsonicInflowRVBdrFaceIntegrator::SubsonicInflowRVBdrFaceIntegrator(
    std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, const int Np,
    const real_t &time, real_t gamma, real_t rho, const Vector &V)
    : BdrFaceIntegrator(liftingScheme,rsolver, Np, time, gamma, true, false),
    r(rho), u(V), rho(std::function<real_t(const Vector&)>()), V(dim, std::function<void(const Vector&, Vector&)>()) {}

void SubsonicInflowRVBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        if (t_dependent)
        {
            rho.SetTime(time);
            V.SetTime(time);
        }
        r = rho.Eval(Tr, ip);
        V.Eval(u, Tr, ip);
    }

    state2(0) = r;
    state2(1) = r * u(0);
    dke = -state1(1) * state1(1);
    state2(num_equations - 1) = state1(num_equations - 1);
    if (dim > 1)
    {
        state2(2) = r * u(1);
        dke -= state1(2) * state1(2);
        if (dim > 2)
        {
            state2(3) = r * u(2);
            dke -= state1(3) * state1(3);
        }
    }
    dke *= 0.5 / state1(0);
    dke += 0.5 * (u * u) * r;
    state2(num_equations - 1) += dke;
}

void SubsonicInflowRVBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx_, const Vector &dqdy_, const Vector &dqdz_, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    dqdx = dqdy = dqdz = 0.0;
    fluxFunction.ComputeViscousFlux(state2, dqdx, dqdy, dqdz, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

void SubsonicInflowRVBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx_, const Vector &dqdy_, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    dqdx = dqdy = 0.0;
    fluxFunction.ComputeViscousFlux(state2, dqdx, dqdy, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

void SubsonicInflowRVBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const Vector &dqdx_, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    dqdx = 0.0;
    fluxFunction.ComputeViscousFlux(state2, dqdx, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}