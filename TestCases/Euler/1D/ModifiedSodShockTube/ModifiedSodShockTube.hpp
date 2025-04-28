#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Sod Shock Tube initial condition
std::function<void(const Vector&, Vector&)> ModifiedSodShockTubeIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 0.5)
        {
            density = 1.0;
            velocity_x = 0.75;
            pressure = 1.0;
        }
        else
        {
            density = 0.125;
            velocity_x = 0.0;
            pressure = 0.1;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterModifiedSodShockTube
{
    RegisterModifiedSodShockTube()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("ModifiedSodShockTubeIC", ModifiedSodShockTubeIC);
        
        // Register boundary condition for left and right boundaries.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("ModifiedSodShockTubeLeftBC", ModifiedSodShockTubeIC);
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("ModifiedSodShockTubeRightBC", ModifiedSodShockTubeIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterModifiedSodShockTube regModifiedSodShockTube;

}