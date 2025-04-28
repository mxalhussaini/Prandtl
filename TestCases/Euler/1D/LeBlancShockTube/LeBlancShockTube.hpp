#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// LeBlanc Shock Tube initial condition
std::function<void(const Vector&, Vector&)> LeBlancShockTubeIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        if (x(0) < 3.0)
        {
            y(0) = 1.0;
            y(1) = 0.0;
            y(2) = 0.1;
        }
        else
        {
            y(0) = 1e-3;
            y(1) = 0.0;
            y(2) = 1e-9;
        }
    };
}

// Registration helper that automatically registers these functions
// along with associated boundary marker arrays.
struct RegisterLeBlancShockTube
{
    RegisterLeBlancShockTube()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("LeBlancShockTubeIC", LeBlancShockTubeIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterLeBlancShockTube regLeBlancShockTube;

}