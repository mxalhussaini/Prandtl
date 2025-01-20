#include "NavierStokesFlux.hpp"
#include "BasicOperations.hpp"
#include "Physics.hpp"

namespace Prandtl
{

real_t NavierStokesFlux::ComputeInviscidFlux(const Vector &state, ElementTransformation &Tr, DenseMatrix &flux) const
{
    return ComputeFlux(state, Tr, flux);
}

void NavierStokesFlux::ComputeViscousFlux(const Vector &state, const DenseMatrix &grad_mat, DenseMatrix &flux) const
{
    Vector prim;
    Conserv2Prim(state, prim);

    grad_mat.GetColumn(0, grad_x_state);
    real_t &drho_dx = grad_x_state(0);
    real_t &du_dx = grad_x_state(1);
    div = grad_x_state(1); 
    real_t &dp_dx = grad_x_state(num_equations - 1);
    cv_dT_dx = Prandtl::gammaM1Inverse / prim(0) * (dp_dx - prim(num_equations - 1) / prim(0) * drho_dx);

#ifdef SUTHERLAND
    mu = ComputeViscosity(prim(0), prim(num_equations - 1));
#endif
    lambda = mu * Prandtl::gamma * Prandtl::PrInverse;
    

    flux(0, 0) = 0.0;
    flux(1, 0) = mu * 2.0 * du_dx;
    flux(num_equations - 1, 0) = lambda * cv_dT_dx + prim(1) * flux(1, 0);
    
    if (dim > 1)
    {
        real_t &dv_dx = grad_x_state(2);
        grad_mat.GetColumn(1, grad_y_state);
        real_t &drho_dy = grad_y_state(0);
        real_t &du_dy = grad_y_state(1);
        real_t &dv_dy = grad_y_state(2);
        real_t &dp_dy = grad_y_state(num_equations - 1);
        div += dv_dy;
        cv_dT_dy = Prandtl::gammaM1Inverse / prim(0) * (dp_dy - prim(num_equations - 1) / prim(0) * drho_dy);

        flux(2, 0) = mu * (du_dy + dv_dx);
        flux(num_equations - 1, 0) += prim(2) * flux(2, 0);

        flux(0, 1) = 0.0;
        flux(1, 1) = flux(2, 0);
        flux(2, 1) = mu * 2.0 * dv_dy;
        flux(num_equations - 1, 1) = lambda * cv_dT_dy + prim(1) * flux(1, 1) + prim(2) * flux(2, 1);


        if (dim > 2)
        {
            real_t &dw_dx = grad_x_state(3);
            real_t &dw_dy = grad_y_state(3);
            grad_mat.GetColumn(2, grad_z_state);
            real_t &drho_dz = grad_z_state(0);
            real_t &du_dz = grad_z_state(1);
            real_t &dv_dz = grad_z_state(2);
            real_t &dw_dz = grad_z_state(3);
            real_t &dp_dz = grad_z_state(3);
            div += dw_dz;
            cv_dT_dz = Prandtl::gammaM1Inverse / prim(0) * (dp_dz - prim(num_equations - 1) / prim(0) * drho_dz);

            flux(3, 0) = mu * (du_dz + dw_dx);
            flux(4, 0) += prim(3) * flux(3, 0);

            flux(3, 1) = mu * (dv_dz + dw_dy);
            flux(4, 1) += prim(3) * flux(3, 1);

            flux(0, 2) = 0.0;
            flux(1, 2) = flux(3, 0);
            flux(2, 2) = flux(3, 1);
            flux(3, 2) = mu * (2.0 * dw_dz - mu_bulk * div);
            flux(4, 2) = prim(1) * flux(1, 2) + prim(2) * flux(2, 2) + prim(3) * flux(3, 2) + lambda * cv_dT_dz;
        }
    }

    flux(1, 0) -= mu * mu_bulk * div;
    flux(num_equations, 0) -= mu * mu_bulk * div * prim(1); 
    flux(2, 1) -= mu * mu_bulk * div;
    flux(num_equations - 1, 1) -= mu * mu_bulk * div * prim(2);
}

real_t NavierStokesFlux::ComputeInviscidFluxDotN(const Vector &x, const Vector &nor, FaceElementTransformations &Tr, Vector &fluxN) const
{
    return ComputeFluxDotN(x, nor, Tr, fluxN);
}


}