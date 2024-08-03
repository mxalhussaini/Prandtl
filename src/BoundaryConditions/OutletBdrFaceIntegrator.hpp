#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class OutletBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    virtual void ComputeOuterState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
public:
    OutletBdrFaceIntegrator(const RiemannSolver &rsolver, const int IntOrderOffset=0);
};


}