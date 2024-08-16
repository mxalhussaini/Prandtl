#include "OutletBdrFaceIntegrator.hpp"

namespace Prandtl
{

OutletBdrFaceIntegrator::OutletBdrFaceIntegrator(
    const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir)
    : BdrFaceIntegrator(rsolver, bdr_face_ir) {}

void OutletBdrFaceIntegrator::ComputeOuterState(
    const Vector &state1,
    Vector &state2,
    FaceElementTransformations &Tr, 
    const IntegrationPoint &ip)
{
    state2 = state1;
}

}
