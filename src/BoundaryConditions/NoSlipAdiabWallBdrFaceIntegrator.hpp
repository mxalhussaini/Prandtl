#pragma once

#include "SlipWallBdrFaceIntegrator.hpp"

namespace Prandtl
{

using namespace mfem;

class NoSlipAdiabWallBdrFaceIntegrator : public SlipWallBdrFaceIntegrator
{
private:
    real_t qn, v;
    Vector V;
    FunctionCoefficient qn_wall;
    VectorFunctionCoefficient V_wall;
public:
    NoSlipAdiabWallBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                     std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                     std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                     const real_t &time, FunctionCoefficient &qn_wall, VectorFunctionCoefficient &V_wall,
                                     bool constant = true, bool t_dependent = false);
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor,
                                   FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceLiftingFlux(const Vector &state1, Vector &state2, Vector &fluxN, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
};

}