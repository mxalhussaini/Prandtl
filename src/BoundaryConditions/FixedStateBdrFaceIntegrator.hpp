#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class FixedStateBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    Vector fixed_state;
    virtual void ComputeOuterState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
public:
    FixedStateBdrFaceIntegrator(const RiemannSolver &rsolver, const Vector &fixed_state_, const int IntOrderOffset=0);
};


}