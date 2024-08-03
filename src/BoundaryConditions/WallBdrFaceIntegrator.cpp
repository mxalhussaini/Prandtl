#include "WallBdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

WallBdrFaceIntegrator::WallBdrFaceIntegrator(
    const RiemannSolver &rsolver,
    const Vector &fixed_state_,
    const int IntOrderOffset=0)
    : BdrFaceIntegrator(rsolver, IntOrderOffset) {}

void WallBdrFaceIntegrator::ComputeOuterState(
    const Vector &state1,
    Vector &state2,
    FaceElementTransformations &Tr, 
    const IntegrationPoint &ip)
{
    Vector nor(state1.Size() - 2);
    if (nor.Size() == 1)  // if 1D, use 1 or -1.
    {
        nor(0) = (Tr.GetElement1IntPoint().x - 0.5) * 2.0;
    }
    else
    {
        CalcOrtho(Tr.Jacobian(), nor);
    }
    state2 = state1;
    RotateState(state2, nor);
    state2(1) *= -1;
    RotateBack(state2, nor);
}

}