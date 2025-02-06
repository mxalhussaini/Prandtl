#include "NoSlipIsothWallBdrFaceIntegrator.hpp"
#include "Physics.hpp"

namespace Prandtl
{

NoSlipIsothWallBdrFaceIntegrator::NoSlipIsothWallBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &T_wall_, VectorFunctionCoefficient &V_wall_,
    bool t_dependent)
    : SlipWallBdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, false, t_dependent),
    T_wall(T_wall_), V_wall(V_wall_) {}

NoSlipIsothWallBdrFaceIntegrator::NoSlipIsothWallBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, real_t &T, const Vector &V)
    : SlipWallBdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, true, false),
    T(T), V(V), T_wall(std::function<real_t(const Vector&)>()), V_wall(dim, std::function<void(const Vector&, Vector&)>())
{
    v = 1.0 / (R_gas * T);
}

void NoSlipIsothWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1,
        Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    fluxFunction.ComputeViscousFlux(state1, grad_mat1, flux_mat);
    flux_mat.Mult(nor, fluxN);
}

void NoSlipIsothWallBdrFaceIntegrator::ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (!constant)
    {
        T_wall.SetTime(time);
        V_wall.SetTime(time);

        T = T_wall.Eval(Tr, ip);
        V_wall.Eval(V, Tr, ip);

        v = 1.0 / (R_gas * T);
    }
    fluxN(0) = state1(0);
    fluxN(1) = V(0) * v;
    if (dim > 1)
    {
        fluxN(2) = V(1) * v;
        if (dim > 2)
        {
            fluxN(3) = V(2) * v;
        }
    }
    fluxN(num_equations - 1) = -v;
    fluxN -= state1;
}

}