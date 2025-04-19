#include "SlipWallBdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

SlipWallBdrFaceIntegrator::SlipWallBdrFaceIntegrator(std::shared_ptr<LiftingScheme> liftingScheme, const NumericalFlux &rsolver, int Np, const real_t &time, real_t gamma, bool constant, bool t_dependent)
    : BdrFaceIntegrator(liftingScheme, rsolver, Np, time, gamma, constant, t_dependent), gammaP1(1.0 + gamma), gammaP1Inverse(1.0 / gammaP1)
{
    unit_nor.SetSize(rsolver.GetFluxFunction().dim);
    prim.SetSize(rsolver.GetFluxFunction().num_equations);
}


real_t SlipWallBdrFaceIntegrator::ComputeBdrFaceInviscidFlux(const Vector &state1, Vector &state2,
    Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    unit_nor = nor;
    Normalize(unit_nor);
    state2 = state1;
    if (nor.Size() > 1)
    {
        RotateState(state2, unit_nor);
    }
    else
    {
        state2(1) = state2(1) * nor(0);
    }

    Conserv2Prim(state2, prim, gammaM1);

    Vector vel(prim.GetData() + 1, nor.Size());
    an = ComputeSoundSpeed(prim(prim.Size() - 1), prim(0), gamma);
    if (prim(1) > 0.0)
    {
        p_star = prim(prim.Size() - 1) +
            0.25 * prim(1) * gammaP1 * prim(0) *
            (prim(1) + std::sqrt(std::pow(prim(1), 2) +
            8.0 * gammaP1Inverse / prim(0) * prim(prim.Size() - 1) * (gammaM1 * gammaP1Inverse + 1.0)));
    }
    else
    {
        p_star = prim(prim.Size() - 1) *
            std::pow(std::max(1.0 + 0.5 * gammaM1 * prim(1) /
            an, 0.0001), 2.0 * gamma * gammaM1Inverse); 
    }

    fluxN = 0.0;
    fluxN(1) = p_star * nor(0);
    if (nor.Size() > 1)
    {
        fluxN(2) = p_star * nor(1);
        if (nor.Size() > 2)
        {
            fluxN(3) = p_star * nor(2);
        }
    }

    return std::sqrt(vel * vel) + an;
}

}