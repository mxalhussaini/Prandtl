#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class SlipWallBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    Vector unit_nor;
    Vector prim;
    real_t p_star, an;
protected:
    const real_t gammaP1, gammaP1Inverse;
public:
    SlipWallBdrFaceIntegrator(std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, int Np, const real_t &time, real_t gamma, bool constant = true, bool t_dependent = false);
    virtual real_t ComputeBdrFaceInviscidFlux(const Vector &state1, Vector &state2, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};


}