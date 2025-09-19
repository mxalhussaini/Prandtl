#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

static real_t gamma = 1.4;
static real_t R_gas = 287.05;
static real_t Ma = 2.0;                                      // Mach number 
static real_t p0 = 303975;                                   // stagnation pressure (Pa)
static real_t a0 = 360.63;                                   // stagnation speed (m/s)
static real_t C = 1 + 0.5 * (gamma-1) * Ma * Ma;             // 1 + (gamma-1)/2 * M^2
static real_t a = a0/std::sqrt(C);                           // speed of sound
static real_t ua = Ma * a;                                   // freestream velocity
static real_t pa = p0 / std::pow(C, gamma/(gamma-1));        // freestream pressure
static real_t pi = pa;                                       // pressure for the inlet and initilization
static real_t Ta = (a * a)/(gamma * R_gas);                  // freestream temperature
static real_t rhoa = pa / (R_gas * Ta);                      // freestream density
static real_t E = pi/(gamma - 1.0) + 0.5 * rhoa * ua * ua;   // Specific total energy


// Supersonic Freestream initial condition
std::function<void(const Vector&, Vector&)> NagashimaIC()
{
    return [] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "Nagashima Ramjet is a 2D problem");

            y(0) = rhoa;
            y(1) = rhoa * ua;
            y(2) = 0.0;
            y(3) = E;
    
    };
}

// Supersonic Freestream boundary condition
const BC_Vector NagashimaInletBCVector({rhoa, rhoa * ua, 0.0, E});

// Inlet velocity ramping boundary condition
std::function<void(const Vector&, real_t, Vector&)>NagashimaRampingInletBC(real_t t_ramp)
{
    return [t_ramp] (const Vector &x, real_t t, Vector &y)
    {
        real_t s    = std::tanh((t - 0.5*t_ramp)/(0.1*t_ramp));
        real_t ramp = std::clamp(0.5*(1.0 + s), 0.0, 1.0);
        real_t u   = Ma * a * ramp;
        y(0) = rhoa;
        y(1) = rhoa*u;
        y(2) = 0.0;
        y(3) = pi/(gamma-1.0) + 0.5*rhoa*u*u;
    };
}

// Subsonic outflow
const real_t NagashimaOutletPressure = pa; 

// Registration helper that automatically registers these functions
struct RegisterNagashima
{
    RegisterNagashima()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition0("NagashimaIC", NagashimaIC);
        // Register boundary condition.
        ConditionFactory::Instance().RegisterVectorBoundaryCondition("NagashimaInletBCVector", NagashimaInletBCVector);
        ConditionFactory::Instance().RegisterVectorTDFunctionBoundaryCondition1("NagashimaRampingInletBC", NagashimaRampingInletBC);
        ConditionFactory::Instance().RegisterScalarBoundaryCondition("NagashimaOutletPressure", NagashimaOutletPressure);

    }
};
// Global static instance to ensure registration happens at startup.
static RegisterNagashima registerNagashima;

}
