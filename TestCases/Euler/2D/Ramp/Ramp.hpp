#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Supersonic Freestream initial condition
std::function<void(const Vector&, Vector&)> RampIC()
{
    return [] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "Ramp is a 2D problem");
        const real_t gamma = 1.4;
        const real_t Ma = 2.0;
        const real_t a = 1; //340.294;
        const real_t rho = 1; //0.225;
        const real_t u = Ma * a;
        const real_t p = rho*a*a/1.4;
        y(0) = rho;
        y(1) = rho * u;
        y(2) = 0.0;
        y(3) = (p / (gamma-1)) + 0.5 * rho * u * u;
    };
}


const BC_Vector RampInletBCVector({1.0, 2.0, 0.0, ((1.0/1.4)/0.4) + 0.5 * 1.0 * 2.0 * 2.0});
// std::function<void(const Vector&, real_t, Vector&)>RampingInletBC(real_t t_ramp)
// {
//     return [t_ramp] (const Vector &x, real_t t, Vector &y)
//     {
//         const real_t gamma = 1.4;
//         const real_t Ma = 2.0;
//         const real_t a = 1;
//         const real_t rho = 1.0;
//         real_t s    = std::tanh((t - 0.5*t_ramp)/(0.1*t_ramp));
//         real_t ramp = std::clamp(0.5*(1.0 + s), 0.0, 1.0);
//         real_t u   = Ma * a * ramp;
//         real_t p   = 1.0; 
//         y(0) = rho;
//         y(1) = rho*u;
//         y(2) = 0.;
//         y(3) = p/(1.4-1.0) + 0.5*rho*u*u;
//     };
// }


// Registration helper that automatically registers these functions
struct RegisterRamp
{
    RegisterRamp()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition0("RampIC", RampIC);
        // Register boundary condition.
        ConditionFactory::Instance().RegisterVectorBoundaryCondition("RampInletBCVector", RampInletBCVector);
        // ConditionFactory::Instance().RegisterVectorTDFunctionBoundaryCondition1("RampingInletBC", RampingInletBC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterRamp registerRamp;

}
