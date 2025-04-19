#pragma once

#include "mfem.hpp"
#include "NumericalFlux.hpp"

namespace Prandtl
{

class ChandrashekarFlux : public NumericalFlux
{
private:
    mutable real_t beta1, beta2;
    mutable real_t p_hat, h_hat, rho_ln, beta_ln;
    mutable real_t rho_mean, u_mean, v_mean, w_mean;
    mutable Vector metric, mom;
    mutable real_t lambda_max, diss;
    mutable real_t drho, du, dv, dw;
    mutable real_t rho1, rho2, u1, u2, v1, v2, w1, w2, p1, p2, V_sq1, V_sq2, V_mag1, V_mag2, Vn, nor_mag;
    const real_t gamma, gammaM1, gammaM1Inverse;
public:
    ChandrashekarFlux(const NavierStokesFlux &fluxFunction, real_t gamma);
    real_t ComputeFaceFlux(const Vector &state1, const Vector &state2, const Vector &nor, Vector &flux) const override;
    real_t ComputeVolumeFlux(const Vector &state1, const Vector &state2, const Vector &metric1, const Vector &metric2, Vector &F_tilde) override;

};

}