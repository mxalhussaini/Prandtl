#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Shu-Osher Shock initial condition
std::function<void(const Vector&, Vector&)> ShuOsherShockIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 1.0)
        {
            density = 3.857;
            velocity_x = 2.629;
            pressure = 10.333;
        }
        else
        {
            density = 1.0 + 0.2 * std::sin(5.0 * (x(0) - 5.0));
            velocity_x = 0.0;
            pressure = 1.0;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterShuOsherShock
{
    RegisterShuOsherShock()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("ShuOsherShockIC", ShuOsherShockIC);

        // Register boundary condition functions.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("ShuOsherShockLeftBC", ShuOsherShockIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterShuOsherShock regShuOsherShock;

}