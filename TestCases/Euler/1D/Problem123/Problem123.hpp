#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// 123 Problem initial condition
std::function<void(const Vector&, Vector&)> Problem123IC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 0.5)
        {
            density = 1.0;
            velocity_x = -2.0;
            pressure = 0.4;
        }
        else
        {
            density = 1.0;
            velocity_x = 2.0;
            pressure = 0.4;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterProblem123
{
    RegisterProblem123()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("Problem123IC", Problem123IC);
        
        // Register boundary condition for left and right boundaries.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("Problem123LeftBC", Problem123IC);
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("Problem123RightBC", Problem123IC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterProblem123 regProblem123;

}