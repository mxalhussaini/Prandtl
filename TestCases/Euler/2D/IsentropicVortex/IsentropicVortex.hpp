#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Isentropic Vortex initial condition
std::function<void(const Vector&, Vector&)> IsentropicVortexIC(real_t radius, real_t Minf, real_t beta, real_t R_gas, real_t gamma)
{
    return [radius, Minf, beta, R_gas, gamma](const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "");

        const real_t xc = 0.0, yc = 0.0;
  
        // Nice units
        const real_t vel_inf = 1.;
        const real_t den_inf = 1.;
  
        // Derive remainder of background state from this and Minf
        const real_t pres_inf = (den_inf / gamma) * (vel_inf / Minf) * (vel_inf / Minf);
        const real_t temp_inf = pres_inf / (den_inf * R_gas);
  
        real_t r2rad = 0.0;
        r2rad += (x(0) - xc) * (x(0) - xc);
        r2rad += (x(1) - yc) * (x(1) - yc);
        r2rad /= (radius * radius);
  
        const real_t shrinv1 = 1.0 / (gamma - 1.0);
  
        const real_t velX = vel_inf * (1 - beta * (x(1) - yc) / radius * std::exp(-0.5 * r2rad));
        const real_t velY = vel_inf * beta * (x(0) - xc) / radius * std::exp(-0.5 * r2rad);
        const real_t vel2 = velX * velX + velY * velY;
  
        const real_t specific_heat = R_gas * gamma * shrinv1;
        const real_t temp = temp_inf - 0.5 * (vel_inf * beta) * (vel_inf * beta) / specific_heat * std::exp(-r2rad);
  
        const real_t den = den_inf * std::pow(temp / temp_inf, shrinv1);
        const real_t pres = den * R_gas * temp;
        const real_t energy = shrinv1 * pres / den + 0.5 * vel2;
  
        y(0) = den;
        y(1) = den * velX;
        y(2) = den * velY;
        y(3) = den * energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterIsentropicVortex
{
    RegisterIsentropicVortex()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition5("IsentropicVortexIC", IsentropicVortexIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterIsentropicVortex regIsentropicVortex;

}