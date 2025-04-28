#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Forward Facing Step initial condition
std::function<void(const Vector&, Vector&)> ForwardFacingStepIC(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "Forward Facing Step is a 2D problem");

        y(0) = 1.4;
        y(1) = 1.4 * 3.0;
        y(2) = 0.0;
        y(3) = 1.0 * gammaM1Inverse + 0.5 * 1.4 * 3.0 * 3.0;
    };
}

// Forward Facing Step boundary condition
const BC_Vector ForwardFacingStepLeftBCVector({1.4, 1.4 * 3.0, 0.0, 1.0 / (1.4 - 1.0) + 0.5 * 1.4 * 3.0 * 3.0});

// Registration helper that automatically registers these functions
struct RegisterForwardFacingStep
{
    RegisterForwardFacingStep()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("ForwardFacingStepIC", ForwardFacingStepIC);
        // Register boundary condition.
        ConditionFactory::Instance().RegisterVectorBoundaryCondition("ForwardFacingStepLeftBCVector", ForwardFacingStepLeftBCVector);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterForwardFacingStep regForwardFacingStep;

}
