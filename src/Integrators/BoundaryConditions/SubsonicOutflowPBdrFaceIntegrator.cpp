#include "SubsonicOutflowPBdrFaceIntegrator.hpp"
#include "Physics.hpp"

namespace Prandtl
{

SubsonicOutflowPBdrFaceIntegrator::SubsonicOutflowPBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &p_fun, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, false, t_dependent),
    p_fun(p_fun) {}

SubsonicOutflowPBdrFaceIntegrator::SubsonicOutflowPBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, real_t p)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false),
    p(p), p_fun(std::function<real_t(const Vector&)>()) {}


void SubsonicOutflowPBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        if (t_dependent)
        {
            p_fun.SetTime(time);
        }
        p = p_fun.Eval(Tr, ip);
    }
    Vector vel(state1.GetData() + 1, dim);
    vel /= state1(0);
    state2 = state1;
    state2(num_equations - 1) = p * gammaM1Inverse + 0.5 * state1(0) * (vel * vel);
}

void SubsonicOutflowPBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1,
        Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}