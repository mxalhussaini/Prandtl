#include "SupersonicInflowBdrFaceIntegrator.hpp"
#include "Physics.hpp"

namespace Prandtl
{

SupersonicInflowBdrFaceIntegrator::SupersonicInflowBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &rho_in_, FunctionCoefficient &p_in_,
    VectorFunctionCoefficient &V_in_, bool constant, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, constant, t_dependent),
    rho_in(rho_in_), p_in(p_in_), V_in(V_in_)
{
    V.SetSize(dim);
    if (constant)
    {
        rho_in.SetTime(time);
        p_in.SetTime(time);
        V_in.SetTime(time);
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
        rho = rho_in.Eval(*Tr, ip);
        p = p_in.Eval(*Tr, ip);
        V_in.Eval(V, *Tr, ip);
    }
}

void SupersonicInflowBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        if (t_dependent)
        {
            rho_in.SetTime(time);
            p_in.SetTime(time);
            V_in.SetTime(time);
        }
        rho = rho_in.Eval(Tr, ip);
        p = p_in.Eval(Tr, ip);
        V_in.Eval(V, Tr, ip);
    }

    state2(0) = rho;
    state2(1) = rho * V(0);
    if (dim > 1)
    {
        state2(2) = rho * V(2);
        if (dim > 2)
        {
            state2(3) = rho * V(3);
        }
    }
    state2(num_equations - 1) = p * gammaM1Inverse + 0.5 * rho * (V * V);
}

void SupersonicInflowBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}