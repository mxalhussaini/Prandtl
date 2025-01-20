#include "SubsonicOutFlowPBdrFaceIntegrator.hpp"
#include "Physics.hpp"

namespace Prandtl
{

SubsonicOutFlowPBdrFaceIntegrator::SubsonicOutFlowPBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &p_out_,
    bool constant, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, constant, t_dependent),
    p_out(p_out_)
{
    if (constant)
    {
        p_out.SetTime(time);
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
        p = p_out.Eval(*Tr, ip);
    }
}

void SubsonicOutFlowPBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        if (t_dependent)
        {
            p_out.SetTime(time);
        }
        p = p_out.Eval(Tr, ip);
    }
    Vector vel(state1.GetData(), dim);
    vel /= state1(0);
    state2 = state1;
    state2(num_equations - 1) = p * gammaM1Inverse + 0.5 * state1(0) * (vel * vel);
}

void SubsonicOutFlowPBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1,
        Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}