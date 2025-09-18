#include "SymmetryBdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

SymmetryBdrFaceIntegrator::SymmetryBdrFaceIntegrator(std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, int Np, const real_t &time, real_t gamma, bool constant, bool t_dependent)
    : BdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, constant, t_dependent) { }


real_t SymmetryBdrFaceIntegrator::ComputeBdrFaceInviscidFlux(const Vector &state1, Vector &state2,
    Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    unit_nor = nor;
    if (unit_nor.Norml2() == 0.0)
    {
        fluxN = 0.0;
        MFEM_WARNING("Zero face normal encountered at symmetry BC, skipping flux.");
        return 0.0;
    }

    Normalize(unit_nor);

    state2 = state1;

    Vector mom(state2.GetData() + 1, unit_nor.Size());
    const real_t mn = mom * unit_nor;
    mom.Add(-2.0 * mn, unit_nor);

    return rsolver.ComputeFaceFlux(state1, state2, nor, fluxN);
}

}






