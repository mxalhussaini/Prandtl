#include "ChandrashekarFlux.hpp"
#include "BasicOperations.hpp"


namespace Prandtl
{

ChandrashekarFlux::ChandrashekarFlux(const NavierStokesFlux &fluxFunction, real_t gamma_)
    : NumericalFlux(fluxFunction), metric(dim), mom(dim), gamma(gamma_), gammaM1(gamma - 1.0), gammaM1Inverse(1.0 / gammaM1) {}

real_t ChandrashekarFlux::ComputeVolumeFlux(const Vector &state1, const Vector &state2,
                                            const Vector &metric1, const Vector &metric2,
                                            Vector &F_tilde)
{
    ComputeMean(metric1, metric2, metric);

    rho1 = state1(0);
    rho2 = state2(0);
    rho_ln = ComputeLogMean(rho1, rho2);
    drho = rho2 - rho1;

    u1 = state1(1) / rho1;
    u2 = state2(1) / rho2;
    u_mean = 0.5 * (u1 + u2);
    du = u2 - u1;
    V_sq1 = u1 * u1;
    V_sq2 = u2 * u2;
    Vn = u_mean * metric(0);
    mom(0) = rho_ln * u_mean;
    h_hat = -0.25 * (u1 * u1 + u2 * u2) + u_mean * u_mean;

    if (dim > 1)
    {
        v1 = state1(2) / rho1;
        v2 = state2(2) / rho2;
        v_mean = 0.5 * (v1 + v2);
        dv = v2 - v1;
        V_sq1 += v1 * v1;
        V_sq2 += v2 * v2;
        Vn += v_mean * metric(1);
        mom(1) = rho_ln * v_mean;
        h_hat += -0.25 * (v1 * v1 + v2 * v2) + v_mean * v_mean;

        if (dim > 2)
        {
            w1 = state1(3) / rho1;
            w2 = state2(3) / rho2;
            w_mean = 0.5 * (w1 + w2);
            dw = w2 - w1;
            V_sq1 += w1 * w1;
            V_sq2 += w2 * w2;
            Vn += w_mean * metric(2);
            mom(2) = rho_ln * w_mean;
            h_hat += -0.25 * (w1 * w1 + w2 * w2) + w_mean * w_mean;
        }
    }

    p1 = gammaM1 * (state1(num_equations - 1) - 0.5 * rho1 * V_sq1);
    p2 = gammaM1 * (state2(num_equations - 1) - 0.5 * rho2 * V_sq2);

    V_mag1 = std::sqrt(V_sq1);
    V_mag2 = std::sqrt(V_sq2);

    lambda_max = V_mag1 + ComputeSoundSpeed(p1, rho1, gamma);
    lambda_max = std::max(lambda_max, V_mag2 + ComputeSoundSpeed(p2, rho2, gamma));

    beta1 = 0.5 * rho1 / p1;
    beta2 = 0.5 * rho2 / p2;
    beta_ln = ComputeLogMean(beta1, beta2);

    p_hat = 0.5 * (rho1 + rho2) / (beta1 + beta2);

    h_hat += 0.5 / beta_ln * gammaM1Inverse + p_hat / rho_ln;

    F_tilde(0) = rho_ln * Vn;
    for (int d = 0; d < dim; d++)
    {
        F_tilde(1 + d) = Vn * mom(d) + p_hat * metric(d);
    }
    F_tilde(1 + dim) = rho_ln * Vn * h_hat;

    return lambda_max;
}

real_t ChandrashekarFlux::ComputeFaceFlux(const Vector &state1, const Vector &state2,
                                          const Vector &nor, Vector &flux) const
{
    nor_mag = nor.Norml2();

    rho1 = state1(0);
    rho2 = state2(0);
    rho_mean = 0.5 * (rho1 + rho2);
    rho_ln = ComputeLogMean(rho1, rho2);
    drho = rho2 - rho1;

    u1 = state1(1) / rho1;
    u2 = state2(1) / rho2;
    u_mean = 0.5 * (u1 + u2);
    du = u2 - u1;
    V_sq1 = u1 * u1;
    V_sq2 = u2 * u2;
    Vn = u_mean * nor(0);
    mom(0) = rho_ln * u_mean;
    h_hat = -0.25 * (u1 * u1 + u2 * u2) + u_mean * u_mean;
    diss = 0.5 * drho * u1 * u2 + rho_mean * du * u_mean;

    if (dim > 1)
    {
        v1 = state1(2) / rho1;
        v2 = state2(2) / rho2;
        v_mean = 0.5 * (v1 + v2);
        dv = v2 - v1;
        V_sq1 += v1 * v1;
        V_sq2 += v2 * v2;
        Vn += v_mean * nor(1);
        mom(1) = rho_ln * v_mean;
        h_hat += -0.25 * (v1 * v1 + v2 * v2) + v_mean * v_mean;
        diss += 0.5 * drho * v1 * v2 + rho_mean * dv * v_mean;

        if (dim > 2)
        {
            w1 = state1(3) / rho1;
            w2 = state2(3) / rho2;
            w_mean = 0.5 * (w1 + w2);
            dw = w2 - w1;
            V_sq1 += w1 * w1;
            V_sq2 += w2 * w2;
            Vn += w_mean * nor(2);
            mom(2) = rho_ln * w_mean;
            h_hat += -0.25 * (w1 * w1 + w2 * w2) + w_mean * w_mean;
            diss += 0.5 * drho * w1 * w2 + rho_mean * dw + w_mean;
        }
    }

    p1 = gammaM1 * (state1(num_equations - 1) - 0.5 * rho1 * V_sq1);
    p2 = gammaM1 * (state2(num_equations - 1) - 0.5 * rho2 * V_sq2);

    V_mag1 = std::sqrt(V_sq1);
    V_mag2 = std::sqrt(V_sq2);

    lambda_max = V_mag1 + ComputeSoundSpeed(p1, rho1, gamma);
    lambda_max = std::max(lambda_max, V_mag2 + ComputeSoundSpeed(p2, rho2, gamma));

    beta1 = 0.5 * rho1 / p1;
    beta2 = 0.5 * rho2 / p2;
    beta_ln = ComputeLogMean(beta1, beta2);

    p_hat = 0.5 * (rho1 + rho2) / (beta1 + beta2);

    h_hat += 0.5 / beta_ln * gammaM1Inverse + p_hat / rho_ln;
    diss += 0.5 * drho * gammaM1Inverse / beta_ln + 0.5 * rho_mean * gammaM1Inverse * (1.0 / beta2 - 1.0 / beta1);

    flux(0) = rho_ln * Vn - 0.5 * lambda_max * (rho2 - rho1) * nor_mag;
    for (int d = 0; d < dim; d++)
    {
        flux(1 + d) = Vn * mom(d) + p_hat * nor(d) - 0.5 * lambda_max * (state2(1 + d) - state1(1 + d)) * nor_mag;
    }
    flux(1 + dim) = rho_ln * Vn * h_hat - 0.5 * lambda_max * diss * nor_mag;

    return lambda_max;
}

}