#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{

using namespace mfem;

class SubsonicInflowPtTtAngBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    real_t p;
    real_t V_sq;
    real_t p0, T0;
    real_t theta, phi;
    Vector V_comps, pT;
    FunctionCoefficient pt, Tt;
public:
    SubsonicInflowPtTtAngBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                           std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                           std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                           const real_t &time, FunctionCoefficient &pt, FunctionCoefficient &Tt,
                                           real_t theta = 0.0, real_t phi = 0.0, bool t_dependent = false);
    SubsonicInflowPtTtAngBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                           std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                           std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                           const real_t &time, real_t pt, real_t Tt,
                                           real_t theta = 0.0, real_t phi = 0.0);

    virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};

}