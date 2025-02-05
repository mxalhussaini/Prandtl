#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{

using namespace mfem;

class SubsonicOutflowPBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    real_t p;
    real_t V_sq;
    FunctionCoefficient p_fun;
public:
    SubsonicOutflowPBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                      std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                      std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                      const real_t &time, FunctionCoefficient &p_out, bool t_dependent = false);
    SubsonicOutflowPBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                      std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                      std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                      const real_t &time, real_t p_out);
    virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};

}