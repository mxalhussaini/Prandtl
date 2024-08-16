#include "ExactRiemannFlux.hpp"
#include "BasicOperations.hpp"


const mfem::real_t gamma = 1.4;

namespace Prandtl
{
real_t ExactRiemannFlux::Eval(const Vector &state1, const Vector &state2,
               const Vector &nor, FaceElementTransformations &Tr,
               Vector &flux) const
{

}

}