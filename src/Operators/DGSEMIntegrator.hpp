#pragma once

#include "mfem.hpp"
#include "NumericalFlux.hpp"

namespace Prandtl
{

using namespace mfem;

class DGSEMIntegrator : public BilinearFormIntegrator
{
public:
    enum Mode
    {
        GRADIENT,
        DIVERGENCE
    };
private:
    Mode currentMode;
    std::shared_ptr<ParMesh> pmesh;
    std::shared_ptr<ParFiniteElementSpace> vfes;
    std::shared_ptr<ParFiniteElementSpace> fes0;
    std::shared_ptr<ParGridFunction> alpha;
    std::shared_ptr<ParGridFunction> grad_x, grad_y, grad_z;
    NumericalFlux &rsolver;
    const NavierStokesFlux &fluxFunction;
    DenseMatrix D, D_mod, D2_mod;
    const int Np_x, Np_y, Np_z;
    const int num_equations, dim, num_elements;
    const L2_SegmentElement *fe;
    IntegrationRules GLIntRules;
    const IntegrationRule *ir, *ir_face, *ir_vol;

    real_t max_char_speed;
    real_t J, J1, J2;
    int id1, id2;

    Vector shape1, shape2;
    Vector state1, state2;
    Vector f, g, h;

    Vector flux_num, flux1, flux2;
    DenseMatrix flux_mat1, flux_mat2, flux_mat_ref;

    DenseMatrix adj1, adj2;
    Vector metric1, metric2;
    Vector nor, nor1, nor2, unit_nor;

    DenseTensor F_inviscid, G_inviscid, H_inviscid;
    // std::vector<DenseMatrix> F_tilde;

    Vector Dx;

    Vector dU_inviscid, dU_viscous, dU, dU_sub;

    DenseTensor SubcellMetricXi, SubcellMetricEta, SubcellMetricZeta;

    Array<int> alpha_indx, vdof_indices, grad_indices1, grad_indices2;
    Vector alpha_dof;

    Vector grad_x_vdofs1, grad_x_vdofs2, grad_y_vdofs1, grad_y_vdofs2, grad_z_vdofs1, grad_z_vdofs2;
    DenseMatrix grad_x_mat1, grad_x_mat2, grad_y_mat1, grad_y_mat2, grad_z_mat1, grad_z_mat2;
    DenseMatrix grad_mat1, grad_mat2;
    Vector grad_state;
    Vector prim;
    int dir;

    Vector grad_xi, grad_eta, grad_zeta;

    real_t lambda;
    real_t div, cv_dT_dx, cv_dT_dy, cv_dT_dz;
    DenseMatrix F_viscous, G_viscous, H_viscous;

    void ComputeSubcellMetrics();
    void ComputeFVFluxes(const DenseMatrix &vdof_mat, real_t alpha_val, ElementTransformation &Tr, DenseMatrix &elvect_mat);

public:
    void SetMode(Mode op)
    {
        currentMode = op;
    }
    inline real_t GetMaxCharSpeed()
    {
        return max_char_speed;
    }

    void ChooseDirection(int dir_)
    {
        dir = dir_;
    }

    DGSEMIntegrator(std::shared_ptr<ParMesh> pmesh,
                    std::shared_ptr<ParFiniteElementSpace> vfes,
                    std::shared_ptr<ParFiniteElementSpace> fes0,
                    std::shared_ptr<ParGridFunction> grad_x,
                    std::shared_ptr<ParGridFunction> grad_y,
                    std::shared_ptr<ParGridFunction> grad_z,
                    std::shared_ptr<ParGridFunction> alpha,
                    NumericalFlux &rsolver, int Np);

    void AssembleFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &elfun, Vector &elvect) override;
    void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &elfun, Vector &elvect) override;
    void AssembleElementMatrix2(const FiniteElement &trial_fe, const FiniteElement &test_fe, ElementTransformation &Trans,  DenseMatrix &elmat) override; 
    
    void AssembleFaceVector(const FiniteElement &el, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, const Vector &el_dudx, const Vector &el_dudy, const Vector &el_dudz, Vector &el_dudt);
    void AssembleFaceVector(const FiniteElement &el, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, const Vector &el_dudx, const Vector &el_dudy, Vector &el_dudt);
    void AssembleFaceVector(const FiniteElement &el, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, const Vector &el_dudx, Vector &el_dudt);

    void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, const Vector &el_dudx, const Vector &el_dudy, const Vector &el_dudz, Vector &el_dudt);
    void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, const Vector &el_dudx, const Vector &el_dudy, Vector &el_dudt);
    void AssembleElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, const Vector &el_dudx, Vector &el_dudt);

    void AssembleLiftingFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy, Vector &el_dudz);
    void AssembleLiftingFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy);
    void AssembleLiftingFaceVector(const FiniteElement &el1, const FiniteElement &el2, FaceElementTransformations &Tr, const Vector &el_u, Vector &el_dudx);
    
    void AssembleLiftingElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy, Vector &el_dudz);
    void AssembleLiftingElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, Vector &el_dudx, Vector &el_dudy);
    void AssembleLiftingElementVector(const FiniteElement &el, ElementTransformation &Tr, const Vector &el_u, Vector &el_dudx);
    // ~DGSEMIntegrator() = default;
};

}