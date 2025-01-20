#pragma once

#include "DGSEMIntegrator.hpp"
#include "BdrFaceIntegrator.hpp"
#include "ModalBasis.hpp"
#include "Indicator.hpp"
#include "BasicOperations.hpp"

namespace Prandtl
{

using namespace mfem;

class DGSEMOperator : public TimeDependentOperator
{
private:
    std::shared_ptr<ParFiniteElementSpace> vfes;
    std::shared_ptr<ParFiniteElementSpace> fes0;
    std::shared_ptr<ParMesh> pmesh;
    std::shared_ptr<ParGridFunction> eta;
    std::shared_ptr<ParGridFunction> alpha;
    std::shared_ptr<ParGridFunction> grad_x, grad_y, grad_z;
    std::unique_ptr<DGSEMIntegrator> integrator;
    std::unique_ptr<Indicator> indicator;
    std::unique_ptr<ParNonlinearForm> nonlinearForm;

    mutable Array<int> vdof_indices;
    mutable Vector el_vdofs, ent_vdofs, grad_vdofs;

    const int Ndofs;

    mutable BlockVector x_blocks;
    std::unique_ptr<ParNonlinearForm> nonlinearForm_Lifting;
    mutable Vector global_entropy;
    mutable std::vector<DenseMatrix> grad_mats;
    mutable DenseMatrix vdof_mat, ent_mat;
    mutable DenseMatrix grad_mat1, grad_mat2;
    mutable DenseMatrix grad_vdof_mats;
    
    const int num_equations, dim, order, num_elements;
    const real_t sharpness_fac = 9.21024;
    const real_t modalThreshold;
    const real_t alpha_min;
    const real_t alpha_max;

    mutable real_t max_char_speed;
    
    std::vector<BdrFaceIntegrator*> bfnfi;
    std::vector<Array<int>> bdr_marker;
    mutable Array<int> ind_indx;
    mutable Vector ind_dof;
    mutable real_t alpha_dof;
    
    void ComputeGlobalEntropyVector(const Vector &x, Vector &y) const;
    void ComputeBlendingCoefficient(const Vector &x) const;

public:
    DGSEMOperator(std::shared_ptr<ParFiniteElementSpace> vfes,
                  std::shared_ptr<ParFiniteElementSpace> fes0,
                  std::shared_ptr<ParMesh> pmesh,
                  std::shared_ptr<ParGridFunction> eta,
                  std::shared_ptr<ParGridFunction> alpha,
                  std::shared_ptr<ParGridFunction> grad_x,
                  std::shared_ptr<ParGridFunction> grad_y,
                  std::shared_ptr<ParGridFunction> grad_z,
                  std::unique_ptr<DGSEMIntegrator> integrator,
                  std::unique_ptr<Indicator> indicator,
                  std::vector<BdrFaceIntegrator*> bfnfi,
                  std::vector<Array<int>> bdr_marker,
                  const real_t alpha_max = 0.5, const real_t alpha_min = 0.001);
    
    void AddBdrFaceIntegrator(BdrFaceIntegrator *bfi, Array<int> &bdr_marker);
    
    void Mult(const Vector &x, Vector &y) const override;
    inline real_t GetMaxCharSpeed()
    {
        return max_char_speed;
    }

    inline real_t& GetTimeRef()
    {
        return t;
    }

};

}