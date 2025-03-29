#include "SpecifiedStateBdrFaceIntegrator.hpp"

namespace Prandtl
{

// Constructor for SpecifiedStateBdrfaceIntegrator with a variable (space- and/or time-dependent) conservative state
SpecifiedStateBdrFaceIntegrator::SpecifiedStateBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, VectorFunctionCoefficient &conserv_state_fun_, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, false, t_dependent),
    conserv_state_fun(conserv_state_fun_) {}


// Constructor for SpecifiedStateBdrfaceIntegrator with a constant conservative state
SpecifiedStateBdrFaceIntegrator::SpecifiedStateBdrFaceIntegrator(
    const NumericalFlux &rsolver, int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, const Vector &conserv_state)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false), 
    conserv_state(conserv_state), conserv_state_fun(num_equations, std::function<void(const Vector&, Vector&)>()) {}

void SpecifiedStateBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2,
    FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (constant)
    {
        state2 = conserv_state;
    }
    else
    {
        if (t_dependent)
        {
            conserv_state_fun.SetTime(time);
        }
        conserv_state_fun.Eval(state2, Tr, ip);
    }
}

void SpecifiedStateBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}