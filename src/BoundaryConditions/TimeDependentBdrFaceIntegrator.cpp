#include "TimeDependentBdrFaceIntegrator.hpp"

namespace Prandtl
{

TimeDependentBdrFaceIntegrator::TimeDependentBdrFaceIntegrator(
    const RiemannSolver &rsolver,
    VectorFunctionCoefficient &vfc_, const real_t &t,
    const int IntOrderOffset=0)
    : BdrFaceIntegrator(rsolver, IntOrderOffset),
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