#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Woodward-Colella Blast Wave initial condition
std::function<void(const Vector&, Vector&)> WoodwardColellaBlastWaveIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        density = 1.0;
        velocity_x = 0.0;

        if (x(0) <= 0.1)
        {
            pressure = 1000.0;
        }
        else if (x(0) <= 0.9)
        {
            pressure = 0.01;
        }
        else
        {
            pressure = 100.0;
        }

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterWoodwardColellaBlastWave
{
    RegisterWoodwardColellaBlastWave()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("WoodwardColellaBlastWaveIC", WoodwardColellaBlastWaveIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterWoodwardColellaBlastWave regWoodwardColellaBlastWave;

}