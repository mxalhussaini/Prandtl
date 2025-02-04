#include "SubsonicInflowRVBdrFaceIntegrator.hpp"

namespace Prandtl
{

SubsonicInflowRVBdrFaceIntegrator::SubsonicInflowRVBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &rho_, VectorFunctionCoefficient V_, bool constant, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, constant, t_dependent),
    rho(rho_), V(V_)
{
    u.SetSize(dim);
    if (constant)
    {
        rho.SetTime(time);
        V.SetTime(time);
        ElementTransformation *Tr = vfes->GetElementTransformation(0);
        IntegrationPoint ip;
        ip.x = 0.0;
        if (dim > 1)
        {
            ip.y = 0.0;
            if (dim > 2)
            {
                ip.z = 0.0;
            }
        }
        r = rho.Eval(*Tr, ip);
        V.Eval(u, *Tr, ip);
    }
}

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

void SubsonicInflowRVBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1,
        Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}