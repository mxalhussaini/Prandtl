#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Triple Point Shock Interaction initial condition
std::function<void(const Vector&, Vector&)> TriplePointShockInteractionIC(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "Triple Point Shock Interaction is a 2D problem");
        real_t density, velocity_x, velocity_y, pressure, energy;
        if (x(0) <= 1.0)
        {
            density = 1.0;
            pressure = 1.0;
        }
        else if (x(0) <= 7.0 && x(1) <= 1.5)
        {
            density = 0.125;
            velocity_x = 0.0;
            pressure = 0.1;
        }
        else
        {
            density = 1.0;
            velocity_x = 0.0;
            pressure = 0.1;
        }
        energy = pressure * gammaM1Inverse;

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = density * velocity_y;
        y(3) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterTriplePointShockInteraction
{
    RegisterTriplePointShockInteraction()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("TriplePointShockInteractionIC", TriplePointShockInteractionIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterTriplePointShockInteraction regTriplePointShockInteraction;

}
