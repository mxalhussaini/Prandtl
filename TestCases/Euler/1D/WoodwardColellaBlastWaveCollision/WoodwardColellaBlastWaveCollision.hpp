#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Woodward-Colella Blast Wave Collision initial condition
std::function<void(const Vector&, Vector&)> WoodwardColellaBlastWaveCollisionIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 0.5)
        {
            density = 5.99924;
            velocity_x = 19.5975;
            pressure = 460.894;
        }
        else
        {
            density = 5.99242;
            velocity_x = -6.19633 ;
            pressure = 46.0950;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterWoodwardColellaBlastWaveCollision
{
    RegisterWoodwardColellaBlastWaveCollision()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("WoodwardColellaBlastWaveCollisionIC", WoodwardColellaBlastWaveCollisionIC);
        
        // Register boundary condition for left and right boundaries.
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("WoodwardColellaBlastWaveCollisionLeftBC", WoodwardColellaBlastWaveCollisionIC);
        ConditionFactory::Instance().RegisterVectorFunctionBoundaryCondition1("WoodwardColellaBlastWaveCollisionRightBC", WoodwardColellaBlastWaveCollisionIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterWoodwardColellaBlastWaveCollision regWoodwardColellaBlastWaveCollision;

}