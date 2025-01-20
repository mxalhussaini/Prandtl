#pragma once

#include "mfem.hpp"
#include "Physics.hpp"

namespace Prandtl
{

using namespace mfem;

class NavierStokesFlux : public EulerFlux
{
private:
    mutable Vector grad_x_state, grad_y_state, grad_z_state;
    mutable real_t div, cv_dT_dx, cv_dT_dy, cv_dT_dz, lambda, mu;


public:
    NavierStokesFlux(const int dim, real_t mu = mu0)
        : EulerFlux(dim, gamma), mu(mu)
    {
        grad_x_state.SetSize(dim + 2);
        if (dim > 1)
        {
            grad_y_state.SetSize(dim + 2);
            if (dim > 2)
            {
                grad_z_state.SetSize(dim + 2);
            }
        }
    }
    real_t ComputeInviscidFlux(const Vector &state, ElementTransformation &Tr, DenseMatrix &flux) const;
    void ComputeViscousFlux(const Vector &state, const DenseMatrix &grad_mat, DenseMatrix &flux) const;
    real_t ComputeInviscidFluxDotN(const Vector &x, const Vector &nor, FaceElementTransformations &Tr, Vector &fluxN) const;
};

}