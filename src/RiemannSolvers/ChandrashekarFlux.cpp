#include "ChandrashekarFlux.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"


namespace Prandtl
{

ChandrashekarFlux::ChandrashekarFlux(const NavierStokesFlux &fluxFunction)
    : NumericalFlux(fluxFunction), metric(fluxFunction.dim),
      metric1(fluxFunction.dim),metric2(fluxFunction.dim), unit_nor(fluxFunction.dim),
      state1(fluxFunction.num_equations), state2(fluxFunction.num_equations) {}

real_t ChandrashekarFlux::ComputeVolumeFlux(const Vector &state1, const Vector &state2,
                                            const Vector &metric1, const Vector &metric2,
                                            Vector &F_tilde)
{
    Vector prim1, prim2;
    Conserv2Prim(state1, prim1);
    Conserv2Prim(state2, prim2);
    Vector vel1(prim1.GetData() + 1, metric1.Size());
    Vector vel2(prim2.GetData() + 1, metric2.Size());
    
    lambda_max = std::sqrt(vel1 * vel1) + ComputeSoundSpeed(prim1(prim1.Size() - 1), prim1(0));
    lambda_max = std::max(lambda_max, std::sqrt(vel2 * vel2) + ComputeSoundSpeed(prim2(prim2.Size() - 1), prim2(0)));

    beta1 = 0.5 * state1(0) / prim1(prim1.Size() - 1);
    beta2 = 0.5 * state2(0) / prim2(prim2.Size() - 1);  

    p_hat = 0.5 * ComputeMean(state1(0), state2(0)) / ComputeMean(beta1, beta2);
    beta_ln = ComputeLogMean(beta1, beta2);
    rho_ln = ComputeLogMean(state1(0), state2(0));
    ComputeMean(metric1, metric2, metric);
    u_mean = ComputeMean(vel1(0), vel2(0));

    drho = ComputeJump(state1(0), state2(0));
    du = ComputeJump(vel1(0), vel2(0));


    h_hat = 0.5 / beta_ln * gammaM1Inverse -
            0.5 * ComputeMean(std::pow(vel1(0), 2.0), std::pow(vel2(0), 2.0))
                + std::pow(u_mean, 2.0) + p_hat / rho_ln;
    q = u_mean * metric(0);

    if (fluxFunction.dim > 1)
    {
        v_mean = ComputeMean(vel1(1), vel2(1));
        dv = ComputeJump(vel1(1), vel2(1));
        q += v_mean * metric(1);
        h_hat -= 0.5 * ComputeMean(std::pow(vel1(1), 2.0), std::pow(vel2(1), 2.0));
        h_hat += std::pow(v_mean, 2.0);
        if (fluxFunction.dim > 2)
        {
            w_mean = ComputeMean(vel1(2), vel2(2));
            dw = ComputeJump(vel1(2), vel2(2));
            q += w_mean * metric(2);
            h_hat -= 0.5 * ComputeMean(std::pow(vel1(2), 2.0), std::pow(vel2(2), 2.0));
            h_hat += std::pow(w_mean, 2.0);
        }
    }

    Vector momentum(fluxFunction.dim);
    momentum(0) = rho_ln * u_mean;
    if(fluxFunction.dim > 1)
    {
        momentum(1) = rho_ln * v_mean;
        if (fluxFunction.dim > 2)
        {
            momentum(2) = rho_ln * w_mean;
        }
    }
    
    F_tilde(0) = momentum * metric;
    const real_t velN = F_tilde(0) / rho_ln;
    for (int d = 0; d < fluxFunction.dim; d++)
    {
        F_tilde(1 + d) = velN * momentum(d) + p_hat * metric(d);
    }
    F_tilde(1 + fluxFunction.dim) = rho_ln * velN * h_hat;


    return lambda_max;
}

real_t ChandrashekarFlux::ComputeFaceFlux(const Vector &state1_, const Vector &state2_,
                                          const Vector &nor, Vector &flux) const
{
    unit_nor = nor;
    Normalize(unit_nor);

    state1 = state1_;
    state2 = state2_;

    // if (nor.Size() > 1)
    // {
    //     RotateState(state1, unit_nor);
    //     RotateState(state2, unit_nor);
    // }

    flux = 0.0;
    Vector prim1, prim2;
    Conserv2Prim(state1, prim1);
    Conserv2Prim(state2, prim2);
    Vector vel1(prim1.GetData() + 1, nor.Size());
    Vector vel2(prim2.GetData() + 1, nor.Size());
    
    lambda_max = std::sqrt(vel1 * vel1) + ComputeSoundSpeed(prim1(prim1.Size() - 1), prim1(0));
    lambda_max = std::max(lambda_max, std::sqrt(vel2 * vel2) + ComputeSoundSpeed(prim2(prim2.Size() - 1), prim2(0)));
    
    beta1 = 0.5 * state1(0) / prim1(prim1.Size() - 1);
    beta2 = 0.5 * state2(0) / prim2(prim2.Size() - 1);

    p_hat = 0.5 * ComputeMean(state1(0), state2(0)) / ComputeMean(beta1, beta2);
    beta_ln = ComputeLogMean(beta1, beta2);

    rho_ln = ComputeLogMean(state1(0), state2(0));
    rho_mean = ComputeMean(state1(0), state2(0));
    drho = ComputeJump(state1(0), state2(0));
    u_mean = ComputeMean(vel1(0), vel2(0));
    du = ComputeJump(vel1(0), vel2(0));
    h_hat = 0.5 / beta_ln * gammaM1Inverse -
            0.5 * ComputeMean(std::pow(vel1(0), 2.0), std::pow(vel2(0), 2.0))
                + p_hat / rho_ln + std::pow(u_mean, 2.0);


    // flux(0) = rho_ln * u_mean;
    // flux(1) = flux(0) * u_mean + p_hat;
    if (fluxFunction.dim > 1)
    {
        v_mean = ComputeMean(vel1(1), vel2(1));
        dv = ComputeJump(vel1(1), vel2(1));
        h_hat -= 0.5 * ComputeMean(std::pow(vel1(1), 2.0), std::pow(vel2(1), 2.0));
        h_hat += std::pow(v_mean, 2.0);
        // flux(2) = flux(0) * v_mean;
        if (fluxFunction.dim > 2)
        {
            w_mean = ComputeMean(vel1(2), vel2(2));
            dw = ComputeJump(vel1(2), vel2(2));
            h_hat -= 0.5 * ComputeMean(std::pow(vel1(2), 2.0), std::pow(vel2(2), 2.0));
            h_hat += std::pow(w_mean, 2.0);
            // flux(3) = flux(0) * w_mean;
        }
    }
    // flux(flux.Size() - 1) = flux(0) * h_hat;

    Vector momentum(fluxFunction.dim);
    momentum(0) = rho_ln * u_mean;
    if(fluxFunction.dim > 1)
    {
        momentum(1) = rho_ln * v_mean;
        if (fluxFunction.dim > 2)
        {
            momentum(2) = rho_ln * w_mean;
        }
    }
    flux(0) = momentum * unit_nor;
    const real_t velN = flux(0) / rho_ln;
    for (int d = 0; d < fluxFunction.dim; d++)
    {
        flux(1 + d) = velN * momentum(d) + p_hat * unit_nor(d);
    }
    flux(1 + fluxFunction.dim) = rho_ln * velN * h_hat;


    for (int i = 0; i < fluxFunction.num_equations - 1; i++)
    {
        flux(i) -= 0.5 * lambda_max * (state2(i) - state1(i));
    }
    diss = 0.5 * gammaM1Inverse / beta_ln + 0.5 * vel1(0) * vel2(0);
    if (fluxFunction.dim > 1)
    {
        diss += 0.5 * vel1(1) * vel2(1);
        if (fluxFunction.dim > 2)
        {
            diss += 0.5 * vel1(2) * vel2(2);
        }
    }
    diss *= drho;
    diss += rho_mean * u_mean * du;
    if (fluxFunction.dim > 1)
    {
        diss += rho_mean * v_mean * dv;
        if (fluxFunction.dim > 2)
        {
            diss += rho_mean * w_mean * dw;
        }
    }


    diss += 0.5 * rho_mean * gammaM1Inverse * (1.0 / beta2 - 1.0 / beta1);

    flux(flux.Size() - 1) -= 0.5 * lambda_max * diss;

    flux *= std::sqrt(nor * nor);

    // if (nor.Size() > 1)
    // {
    //     RotateBack(flux, unit_nor);
    // }

    return lambda_max;
}

}