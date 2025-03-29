#include "SubsonicInflowPtTtAngBdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{

SubsonicInflowPtTtAngBdrFaceIntegrator::SubsonicInflowPtTtAngBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &pt, FunctionCoefficient &Tt,
    real_t theta, real_t phi, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, false, t_dependent),
    pt(pt), Tt(Tt), theta(theta), phi(phi)
{
    V_comps.SetSize(dim);

    if (dim > 1)
    {
        V_comps(0) = std::cos(theta * M_PI / 180.0);
        V_comps(1) = std::sin(theta * M_PI / 180.0);
        if (dim > 2)
        {
            V_comps(0) *= std::sin(phi * M_PI / 180.0);
            V_comps(1) *= std::sin(phi * M_PI / 180.0);
            V_comps(3) = std::cos(phi * M_PI / 180.0);
        }
    }
}

SubsonicInflowPtTtAngBdrFaceIntegrator::SubsonicInflowPtTtAngBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, real_t pt, real_t Tt,
    real_t theta, real_t phi)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false),
    p0(pt), T0(Tt), theta(theta), phi(phi),
    pt(std::function<real_t(const Vector&)>()), Tt(std::function<real_t(const Vector&)>())
{
    V_comps.SetSize(dim);

    if (dim > 1)
    {
        V_comps(0) = std::cos(theta * M_PI / 180.0);
        V_comps(1) = std::sin(theta * M_PI / 180.0);
        if (dim > 2)
        {
            V_comps(0) *= std::sin(phi * M_PI / 180.0);
            V_comps(1) *= std::sin(phi * M_PI / 180.0);
            V_comps(3) = std::cos(phi * M_PI / 180.0);
        }
    }
}

void SubsonicInflowPtTtAngBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        if (t_dependent)
        {
            pt.SetTime(time);
            Tt.SetTime(time);
        }
        p0 = pt.Eval(Tr, ip);
        T0 = Tt.Eval(Tr, ip);
    }

    p = ComputePressure(state1);
    V_sq = 2.0 * cp * T0 * (1.0 - std::pow(p / p0, gammaM1_gammaInverse));
    V_sq = std::max(0.0, V_sq);

    state2(0) = gamma_gammaM1Inverse * p / (cp * T0 - 0.5 * V_sq);
    state2(1) = state2(0) * std::sqrt(V_sq) * V_comps(0);
    if (dim > 1)
    {
        state2(2) = state2(0) * std::sqrt(V_sq) * V_comps(1);
        if (dim > 2)
        {
            state2(3) = state2(0) * std::sqrt(V_sq) * V_comps(2);
        }
    }
    state2(num_equations - 1) = p * gammaM1Inverse + 0.5 * state2(0) * V_sq;
}

void SubsonicInflowPtTtAngBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1,
        Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    fluxFunction.ComputeViscousFlux(state2, grad_mat1, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}