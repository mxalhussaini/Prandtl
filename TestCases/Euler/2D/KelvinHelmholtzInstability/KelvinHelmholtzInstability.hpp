#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Kelvin Helmholtz initial condition
std::function<void(const Vector&, Vector&)> KelvinHelmholtzInstabilyIC(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "KHI is a 2D problem");

        real_t density, velocity_x, velocity_y, pressure, energy, B;

        B = std::tanh(15.0 * x(1) + 7.5) - std::tanh(15.0 * x(1) - 7.5);
        density = 0.5 + 0.75 * B;
        velocity_x = 0.5 * (B - 1.0);
        velocity_y = 0.1 * std::sin(2.0 * M_PI * x(0));
        pressure = 1.0;
        energy = pressure * gammaM1Inverse + density * 0.5 * (velocity_x * velocity_x + velocity_y * velocity_y);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = density * velocity_y;
        y(3) = energy;

        // if (x(1) < 0.5 && x(1) > -0.5)
        // {
        //     y(0) = 2.0;
        //     y(1) = -1.0;
        //     y(2) = 2.0 * 0.01 * std::sin(M_PI * x(0));
        //     y(3) = 2.5 * gammaM1Inverse  + 0.5 * (y(1) * y(1) + y(2) * y(2)) / y(0);
        // }
        // else
        // {
        //     y(0) = 1.0;
        //     y(1) = 0.5;
        //     y(2) = 0.01 * std::sin(M_PI * x(0));;
        //     y(3) = 2.5 * gammaM1Inverse + 0.5 * (y(1) * y(1) + y(2) * y(2)) / y(0);
        // }
    }; 
}

// Registration helper that automatically registers these functions
// along with associated boundary marker arrays.
struct RegisterKelvinHelmholtzInstability
{
    RegisterKelvinHelmholtzInstability()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("KelvinHelmholtzInstabilityIC", KelvinHelmholtzInstabilyIC);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterKelvinHelmholtzInstability regKelvinHelmholtzInstability;

}
