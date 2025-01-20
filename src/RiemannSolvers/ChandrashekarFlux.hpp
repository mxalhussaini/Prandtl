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
    mutable Vector metric, metric1, metric2;
    mutable Vector state1, state2;
    mutable real_t q;
    mutable Vector unit_nor;
    mutable real_t lambda_max;
    mutable real_t diss;
    mutable real_t drho, du, dv, dw;

public:
    ChandrashekarFlux(const NavierStokesFlux &fluxFunction);
    real_t ComputeFaceFlux(const Vector &state1, const Vector &state2, 
                           const Vector &nor, Vector &flux) const override;
    
    real_t ComputeVolumeFlux(const Vector &state1, const Vector &state2,
                             const Vector &metric1, const Vector &metric2,
                             Vector &F_tilde) override;

};

}