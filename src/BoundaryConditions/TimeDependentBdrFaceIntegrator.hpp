#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class TimeDependentBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    const real_t &time;
    VectorFunctionCoefficient vfc;

    virtual void ComputeOuterState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
public:
    TimeDependentBdrFaceIntegrator(const RiemannSolver &rsolver,
        VectorFunctionCoefficient &vfc_, const real_t &t,  const int IntOrderOffset=0);
};


}