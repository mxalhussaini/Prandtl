#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class OutletBdrFaceIntegrator : public BdrFaceIntegrator
{
public:
    OutletBdrFaceIntegrator(const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir);
    virtual void ComputeOuterState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};


}