#include "SupersonicInflowBdrFaceIntegrator.hpp"

namespace Prandtl
{

// Constructor for SupersonicInflowBdrFaceIntegrator with a variable (space- and/or time-dependent) conservative state
SupersonicInflowBdrFaceIntegrator::SupersonicInflowBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, VectorFunctionCoefficient &conserv_state_fun, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, false, t_dependent),
    conserv_state_fun(conserv_state_fun) {}

// Constructor for SupersonicInflowBdrFaceIntegrator with a constant (space- and/or time-dependent) conservative state
SupersonicInflowBdrFaceIntegrator::SupersonicInflowBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, const Vector &conserv_state)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false),
    conserv_state(conserv_state), conserv_state_fun(num_equations, std::function<void(const Vector&, Vector&)>())
{
}

void SupersonicInflowBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
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

void SupersonicInflowBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}