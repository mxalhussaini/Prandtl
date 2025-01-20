#include "SpecifiedStateBdrFaceIntegrator.hpp"

namespace Prandtl
{

SpecifiedStateBdrFaceIntegrator::SpecifiedStateBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, VectorFunctionCoefficient &state_fun_, bool constant, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, constant, t_dependent),
    state_fun(state_fun_)
{
    state.SetSize(num_equations);
    if (constant)
    {
        state_fun.SetTime(time);
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
        state_fun.Eval(state, *Tr, ip);        
    }
}

void SpecifiedStateBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2,
    FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        if (t_dependent)
        {
            state_fun.SetTime(time);
        }
    }
    state_fun.Eval(state2, Tr, ip);
}

void SpecifiedStateBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}