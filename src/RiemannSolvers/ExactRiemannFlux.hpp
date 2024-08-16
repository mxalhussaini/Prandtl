#pragma once

#include "mfem.hpp"

namespace Prandtl
{

using namespace mfem;

class ExactRiemannFlux : public RiemannSolver
{
public:
    ExactRiemannFlux(const FluxFunction &FluxFunction)
        : RiemannSolver(fluxFunction) {}

    real_t Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, FaceElementTransformations &Tr,
               Vector &flux) const override;

private:
    mutable real_t p_star, a;

};



}