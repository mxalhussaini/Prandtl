#include "RiemannInvariantBdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{
// Constructor for RiemannInvariantBdrFaceIntegrator with a variable (space- and/or time-dependent) primitive state
RiemannInvariantBdrFaceIntegrator::RiemannInvariantBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, VectorFunctionCoefficient &prim_state_fun, bool t_dependent)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, false, t_dependent),
    prim_state_fun(prim_state_fun)
{
    unit_nor.SetSize(dim);
    V_o.SetSize(dim);
    prim_state.SetSize(num_equations);
}

// Constructor for RiemannInvariantBdrFaceIntegrator with a constant primitive state
RiemannInvariantBdrFaceIntegrator::RiemannInvariantBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, const Vector &prim_state)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false),
    prim_state_fun(num_equations, std::function<void(const Vector&, Vector&)>())
{
    unit_nor.SetSize(dim);
    V_o.SetSize(dim);

    rho_o = prim_state(0);
    p_o = prim_state(num_equations - 1);
    V_o(0) = prim_state(1);
    if (dim > 1)
    {
        V_o(1) = prim_state(2);
        if (dim > 2)
        {
            V_o(2) = prim_state(3);
        }
    }
}


void RiemannInvariantBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2,
    FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    unit_nor = nor;
    Normalize(unit_nor);

    if (!constant)
    {
        if (t_dependent)
        {
            prim_state_fun.SetTime(time);
        }
  
        prim_state_fun.Eval(prim_state, Tr, ip);
        rho_o = prim_state(0);
        p_o = prim_state(num_equations - 1);
        V_o(0) = prim_state(1);
        if (dim > 1)
        {
            V_o(1) = prim_state(2);
            if (dim > 2)
            {
                V_o(2) = prim_state(3);
            }
        }
    }

    Vector V_i(state1.GetData() + 1, dim);
    V_i /= state1(0);
    Vn_i = V_i * unit_nor;
    Vn_o = V_o * unit_nor;

    p_i = ComputePressure(state1);

    a_o = ComputeSoundSpeed(p_o, rho_o);
    a_i = ComputeSoundSpeed(p_i, state1(0));

    R_minus = (std::abs(Vn_o) >= a_o && Vn_i >= 0.0) ? 
               Vn_i - 2.0 * a_i * gammaM1Inverse : Vn_o - 2.0 * a_o * gammaM1Inverse;
    R_plus =  (std::abs(Vn_o) >= a_o && Vn_i < 0.0) ? 
               Vn_o + 2.0 * a_o * gammaM1Inverse : Vn_i + 2.0 * a_i * gammaM1Inverse;

    Vn_b = 0.5 * (R_plus + R_minus);
    a_b = 0.25 * gammaM1 * (R_plus - R_minus);
    
    dVn_i = Vn_b - Vn_i;
    dVn_o = Vn_b - Vn_o;
    
    state2(0) = (Vn_i < 0.0) ?
                 std::pow(a_b * a_b * gammaInverse * rho_o / p_o, gammaM1Inverse) * rho_o :
                 std::pow(a_b * a_b * gammaInverse * state1(0) / p_i, gammaM1Inverse) * state1(0);
    state2(1) = (Vn_i >= 0.0) ?
                 state2(0) * (V_i(0) + dVn_i * unit_nor(0)) :
                 state2(0) * (V_o(0) + dVn_o * unit_nor(0));
    state2(num_equations - 1) = state2(1) * state2(1);
    if (dim > 1)
    {
        state2(2) = (Vn_i >= 0.0) ?
                    state2(0) * (V_i(1) + dVn_i * unit_nor(1)) :
                    state2(0) * (V_o(1) + dVn_o * unit_nor(1));
        state2(num_equations - 1) += state2(2) * state2(2);
        if (dim > 2)
        {
            state2(3) = (Vn_i >= 0.0) ?
                         state2(0) * (V_i(2) + dVn_i * unit_nor(2)) :
                         state2(0) * (V_o(2) + dVn_o * unit_nor(2));
            state2(num_equations - 1) += state2(3) * state2(3);
        }
    }
    p_b = a_b * a_b * state2(0) * gammaInverse;
    state2(num_equations - 1) *= (0.5 / state2(0));
    state2(num_equations - 1) += p_b * gammaM1Inverse;
}

void RiemannInvariantBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    grad_mat2 = 0.0;
    fluxFunction.ComputeViscousFlux(state2, grad_mat2, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}