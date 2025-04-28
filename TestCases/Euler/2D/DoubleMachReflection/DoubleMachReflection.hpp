#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Double Mach Reflection initial condition
std::function<void(const Vector&, Vector&)> DoubleMachReflectionIC(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        if (x(0) < 1.0 / 6.0 + x(1) / std::sqrt(3))
        {
            y(0) = 8.0;
            y(1) = 8.0 * 7.144709581221619;
            y(2) = -8.0 * 4.125;
            y(3) = 116.5 * gammaM1Inverse  + 0.5 * y(0) * (7.144709581221619 * 7.144709581221619 + 4.125 * 4.125);
        }
        else
        {
            y(0) = 1.4;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 1.0 * gammaM1Inverse;
        }
    };    
}

// Double Mach reflection boundary condition for top boundary
std::function<void(const Vector&, real_t, Vector&)> DoubleMachReflectionTopBCFunction(const real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, real_t t, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        if (x(0) < 1.0 / 6.0 + (x(1) + 20.0 * t) / std::sqrt(3))
        {
            y(0) = 8.0;
            y(1) = 8.0 * 7.144709581221619;
            y(2) = -8.0 * 4.125;
            y(3) = 116.5 * gammaM1Inverse + 0.5 * y(0) * (7.144709581221619 * 7.144709581221619 + 4.125 * 4.125);
        }
        else
        {
            y(0) = 1.4;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 1.0 * gammaM1Inverse;
        }
    };
}

// Double Mach Reflection conservative state boundary condition vector for left and bottom boundaries
const BC_Vector DoubleMachReflectionLeftBottom1BCVector({8.0, 8.0 * 7.144709581221619, -8.0 * 4.125, 116.5 * 1.0 / (1.4 - 1.0) + 0.5 * 8.0 * (7.144709581221619 * 7.144709581221619 + 4.125 * 4.125)}); 

// Registration helper that automatically registers these functions
struct RegisterDoubleMachReflection
{
    RegisterDoubleMachReflection()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("DoubleMachReflectionIC", DoubleMachReflectionIC);
        
        // Register boundary conditions
        ConditionFactory::Instance().RegisterVectorBoundaryCondition("DoubleMachReflectionLeftBottom1BCVector", DoubleMachReflectionLeftBottom1BCVector);
        ConditionFactory::Instance().RegisterVectorTDFunctionBoundaryCondition1("DoubleMachReflectionTopBCFunction", DoubleMachReflectionTopBCFunction);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterDoubleMachReflection regDoubleMachReflection;

}
