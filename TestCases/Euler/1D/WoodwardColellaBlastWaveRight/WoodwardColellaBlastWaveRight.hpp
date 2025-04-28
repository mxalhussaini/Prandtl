#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Woodward-Colella Blast Wave Right Half initial condition
std::function<void(const Vector&, Vector&)> WoodwardColellaBlastWaveRightIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 0.5)
        {
            density = 1.0;
            velocity_x = 0.0;
            pressure = 0.01;
        }
        else
        {
            density = 1.0;
            velocity_x = 0.0;
            pressure = 100.0;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
// along with associated boundary marker arrays.
struct RegisterWoodwardColellaBlastWaveRight
{
    RegisterWoodwardColellaBlastWaveRight()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("WoodwardColellaBlastWaveRightIC", WoodwardColellaBlastWaveRightIC);
        
        // Register boundary condition for left and right boundaries.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("WoodwardColellaBlastWaveRightLeftBC", WoodwardColellaBlastWaveRightIC);
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("WoodwardColellaBlastWaveRightRightBC", WoodwardColellaBlastWaveRightIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterWoodwardColellaBlastWaveRight regWoodwardColellaBlastWaveRight;

}