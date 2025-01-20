# pragma once

#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{
using namespace mfem;

class SupersonicInflowBdrFaceIntegrator : public BdrFaceIntegrator
{
private:
    real_t rho, p;
    Vector V;
    FunctionCoefficient rho_in, p_in;
    VectorFunctionCoefficient V_in;
public:
    SupersonicInflowBdrFaceIntegrator(const NumericalFlux &rsolver, const int Np,
                                      std::shared_ptr<ParGridFunction> grad_x, std::shared_ptr<ParGridFunction> grad_y,
                                      std::shared_ptr<ParGridFunction> grad_z, std::shared_ptr<ParFiniteElementSpace> vfes,
                                      const real_t &time, FunctionCoefficient &rho_in, FunctionCoefficient &p_in,
                                      VectorFunctionCoefficient &V_in, bool constant = true, bool t_dependent = false);
    virtual void ComputeOuterInviscidState(const Vector &state1, Vector &state2, FaceElementTransformations &Tr, const IntegrationPoint &ip) override;
    virtual void ComputeBdrFaceViscousFlux(const Vector &state1, const Vector &state2, const DenseMatrix &grad_mat1, Vector &fluxN, const Vector &nor, FaceElementTransformations &Tr, const IntegrationPoint &ip) override; 
};

}