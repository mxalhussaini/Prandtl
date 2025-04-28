#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Woodward-Colella Blast Wave Left Half initial condition
std::function<void(const Vector&, Vector&)> WoodwardColellaBlastWaveLeftIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 0.5)
        {
            density = 1.0;
            velocity_x = 0.0;
            pressure = 1000.0;
        }
        else
        {
            density = 1.0;
            velocity_x = 0.0;
            pressure = 0.01;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
// along with associated boundary marker arrays.
struct RegisterWoodwardColellaBlastWaveLeft
{
    RegisterWoodwardColellaBlastWaveLeft()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("WoodwardColellaBlastWaveLeftIC", WoodwardColellaBlastWaveLeftIC);
        
        // Register boundary condition for left and right boundaries.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("WoodwardColellaBlastWaveLeftLeftBC", WoodwardColellaBlastWaveLeftIC);
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("WoodwardColellaBlastWaveLeftRightBC", WoodwardColellaBlastWaveLeftIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterWoodwardColellaBlastWaveLeft regWoodwardColellaBlastWaveLeft;

}