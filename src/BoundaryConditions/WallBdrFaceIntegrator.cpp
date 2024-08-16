#include "WallBdrFaceIntegrator.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{

WallBdrFaceIntegrator::WallBdrFaceIntegrator(
    const RiemannSolver &rsolver, const IntegrationRule *bdr_face_ir)
    : BdrFaceIntegrator(rsolver, bdr_face_ir)
{
    unit_nor.SetSize(rsolver.GetFluxFunction().dim);
    prim.SetSize(rsolver.GetFluxFunction().dim);
}

void WallBdrFaceIntegrator::ComputeOuterState(
    const Vector &state1,
    Vector &state2,
    FaceElementTransformations &Tr, 
    const IntegrationPoint &ip)
{
    unit_nor = nor;
    Normalize(unit_nor);
    state2 = state1;
    RotateState(state2, unit_nor);
    state2(1) *= -1;
    RotateBack(state2, unit_nor);
}

real_t WallBdrFaceIntegrator::ComputeBdrFaceFlux(const Vector &state1,
    Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    unit_nor = nor;
    Normalize(unit_nor);
    state2 = state1;
    RotateState(state2, unit_nor);

    prim = Conserv2Prim(state2, gamma);

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