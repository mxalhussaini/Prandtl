#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class WallBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    Vector unit_nor;
    Vector prim;
    real_t p_star, an;
public:
    WallBdrFaceIntegrator(const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir);
    virtual void ComputeOuterState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual real_t ComputeBdrFaceFlux(const Vector &state1, Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};


}