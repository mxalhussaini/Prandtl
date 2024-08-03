#include "OutletBdrFaceIntegrator.hpp"

namespace Prandtl
{

OutletBdrFaceIntegrator::OutletBdrFaceIntegrator(
    const RiemannSolver &rsolver,
    const int IntOrderOffset=0)
    : BdrFaceIntegrator(rsolver, IntOrderOffset) {}

void OutletBdrFaceIntegrator::ComputeOuterState(
    const Vector &state1,
    Vector &state2,
    FaceElementTransformations &Tr, 
    const IntegrationPoint &ip)
{
    state2 = state1;
}

}