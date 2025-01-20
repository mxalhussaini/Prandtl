#pragma once

#include "mfem.hpp"

namespace Prandtl
{

using namespace mfem;

extern const real_t gamma;
extern const real_t gammaInverse;
extern const real_t gammaP1;
extern const real_t gammaM1;
extern const real_t gammaP1Inverse;
extern const real_t gammaM1Inverse;
extern const real_t gamma_gammaM1Inverse; // gamma * gammaM1Inverse;
extern const real_t gammaM1_gammaInverse; // gammaM1 * gammaInverse; 

extern const real_t Pr;
extern const real_t PrInverse;

extern const real_t R_gas;
extern const real_t cp;

const real_t mu_bulk = 2.0 / 3.0;
const real_t mu0 = 1.716e-5;
const real_t T0 = 273.15;
const real_t Ts = 110.4;
const real_t T0pTs = T0 + Ts;

inline real_t ComputeViscosity(const real_t &rho, const real_t &p)
{
    real_t T = p / (rho * R_gas);
    return mu0 * T0pTs / (T + Ts) * (T / T0) * std::sqrt(T / T0);
}

}