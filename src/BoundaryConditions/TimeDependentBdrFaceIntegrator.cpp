#include "TimeDependentBdrFaceIntegrator.hpp"

namespace Prandtl
{

TimeDependentBdrFaceIntegrator::TimeDependentBdrFaceIntegrator(
    const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir,
    VectorFunctionCoefficient &vfc_, const real_t &t)
    : BdrFaceIntegrator(rsolver, bdr_face_ir),
      time(t), vfc(vfc_) {}

void TimeDependentBdrFaceIntegrator::ComputeOuterState(
    const Vector &state1,
    Vector &state2,
    FaceElementTransformations &Tr, 
    const IntegrationPoint &ip)
{
    vfc.SetTime(time);
    vfc.Eval(state2, Tr, ip);
}

}