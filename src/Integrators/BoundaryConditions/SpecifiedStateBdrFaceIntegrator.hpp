#pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl

{

using namespace mfem;

class SpecifiedStateBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    VectorFunctionCoefficient conserv_state_fun;
    Vector conserv_state;
public:
    SpecifiedStateBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                    const real_t &time, VectorFunctionCoefficient &conserv_state_fun, bool t_dependent = false);
    SpecifiedStateBdrFaceIntegrator(const NumericalFlux &rsolver, int Np,
                                    std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                    std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                    const real_t &time, const Vector &conserv_state);
    virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};

}