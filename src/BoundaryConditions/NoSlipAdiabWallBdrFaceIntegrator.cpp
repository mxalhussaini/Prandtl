#include "NoSlipAdiabWallBdrFaceIntegrator.hpp"

namespace Prandtl
{

NoSlipAdiabWallBdrFaceIntegrator::NoSlipAdiabWallBdrFaceIntegrator(
    const NumericalFlux &rsolver, const int Np,
    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
    const real_t &time, FunctionCoefficient &qn_wall_, VectorFunctionCoefficient &V_wall_,
    bool constant, bool t_dependent)
    : SlipWallBdrFaceIntegrator(rsolver, Np, grad_x, grad_y, grad_z, vfes, time, constant, t_dependent),
    qn_wall(qn_wall_), V_wall(V_wall_)
{
    V.SetSize(dim);
    if (!t_dependent)
    {
        qn_wall.SetTime(time);
        V_wall.SetTime(time);
        ElementTransformation *Tr = vfes->GetElementTransformation(0);
        IntegrationPoint ip;
        ip.x = 0.0;
        if (dim > 1)
        {
            ip.y = 0.0;
            if (dim > 2)
            {
                ip.z = 0.0;
            }
        }
        qn = qn_wall.Eval(*Tr, ip);
        V_wall.Eval(V, *Tr, ip);
    }
}

void NoSlipAdiabWallBdrFaceIntegrator::ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1,
        Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (t_dependent)
    {
        qn_wall.SetTime(time);
        V_wall.SetTime(time);

        qn = qn_wall.Eval(Tr, ip);
        V_wall.Eval(V, Tr, ip);
    }

    qn *= std::sqrt(nor * nor);

    fluxFunction.ComputeViscousFlux(state1, grad_mat1, flux_mat);
    flux_mat.Mult(nor, fluxN);
    fluxN(num_equations - 1) = V(0) * fluxN(1) + qn;
    if (dim > 1)
    {
        fluxN(num_equations - 1) += V(1) * fluxN(2);
        if (dim > 2)
        {
            fluxN(num_equations - 1) += V(2) * fluxN(3);
        }
    }
}

void NoSlipAdiabWallBdrFaceIntegrator::ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip)
{
    if (t_dependent)
    {
        V_wall.SetTime(time);
        V_wall.Eval(V, Tr, ip);        
    }
    v = -state1(num_equations - 1);
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