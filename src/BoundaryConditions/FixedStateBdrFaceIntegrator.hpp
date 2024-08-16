#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class FixedStateBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    Vector fixed_state;
public:
    FixedStateBdrFaceIntegrator(const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir, const Vector &fixed_state_);
    virtual void ComputeOuterState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};


}