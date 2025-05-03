#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Backward Facing Step initial condition
std::function<void(const Vector&, Vector&)> BackwardFacingStepIC(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "Backward Facing Step is a 2D problem");
        real_t density, velocity_x, velocity_y, pressure, energy;
        if (x(0) < 0.5)
        {
            density = 5.9970;
            velocity_x = 98.5914;
            pressure = 11666.5;
        }
        else
        {
            density = 1.0;
            velocity_x = 0.0;
            pressure = 1.0;
        }
        velocity_y = 0.0;
        energy = pressure * gammaM1Inverse + density * 0.5 * (velocity_x * velocity_x + velocity_y * velocity_y);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = density * velocity_y;
        y(3) = energy;
    };
}

// Backward Facing Step boundary condition
const BC_Vector BackwardFacingStepLeftBCVector({5.9970, 5.9970 * 98.5914, 0.0, 11666.5 / (1.4 - 1.0) + 0.5 * 5.9970 * 98.5914 * 98.5914});

// Registration helper that automatically registers these functions
struct RegisterBackwardFacingStep
{
    RegisterBackwardFacingStep()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("BackwardFacingStepIC", BackwardFacingStepIC);
        // Register boundary condition.
        ConditionFactory::Instance().RegisterVectorBoundaryCondition("BackwardFacingStepLeftBCVector", BackwardFacingStepLeftBCVector);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterBackwardFacingStep regBackwardFacingStep;

}
