#include "SupersonicOutflowBdrFaceIntegrator.hpp"

namespace Prandtl
{

SupersonicOutflowBdrFaceIntegrator::SupersonicOutflowBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time)
    : BdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false) {}

void SupersonicOutflowBdrFaceIntegrator::ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    state2 = state1;
}

void SupersonicOutflowBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    fluxFunction.ComputeViscousFlux(state2, grad_mat1, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

}