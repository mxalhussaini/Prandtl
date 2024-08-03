#include "FixedStateBdrFaceIntegrator.hpp"

namespace Prandtl
{

FixedStateBdrFaceIntegrator::FixedStateBdrFaceIntegrator(
    const RiemannSolver &rsolver,
    const Vector &fixed_state_,
    const int IntOrderOffset=0)
    : BdrFaceIntegrator(rsolver, IntOrderOffset),
      fixed_state(fixed_state_) {}

void FixedStateBdrFaceIntegrator::ComputeOuterState(
    const Vector &state1,
    Vector &state2,
    FaceElementTransformations &Tr, 
    const IntegrationPoint &ip)
{
    state2 = fixed_state;
}

}
