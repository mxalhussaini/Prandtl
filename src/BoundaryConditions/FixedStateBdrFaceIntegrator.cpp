#include "FixedStateBdrFaceIntegrator.hpp"

namespace Prandtl
{

FixedStateBdrFaceIntegrator::FixedStateBdrFaceIntegrator(const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir,
    const Vector &fixed_state_)
    : BdrFaceIntegrator(rsolver, bdr_face_ir),
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
