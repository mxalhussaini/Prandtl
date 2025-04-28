#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Lax Shock Tube initial condition
std::function<void(const Vector&, Vector&)> LaxShockTubeIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 0.5)
        {
            density = 0.445;
            velocity_x = 0.698;
            pressure = 3.528;
        }
        else
        {
            density = 0.5;
            velocity_x = 0.0;
            pressure = 0.571;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
// along with associated boundary marker arrays.
struct RegisterLaxShockTube
{
    RegisterLaxShockTube()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("LaxShockTubeIC", LaxShockTubeIC);
        
        // Register boundary condition for left and right boundaries.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("LaxShockTubeLeftBC", LaxShockTubeIC);
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("LaxShockTubeRightBC", LaxShockTubeIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterLaxShockTube regLaxShockTube;

}