#pragma once

#include "mfem.hpp"
#include "Filter.hpp"
#include "BdrFaceIntegrator.hpp"
#include "ModalBasis.hpp"

namespace Prandtl
{

using namespace mfem;

class EntropyFilter : public Filter
{
private:
    // mutable GridFunction zeta;

    std::shared_ptr<ParFiniteElementSpace> vfes;
    std::shared_ptr<ParFiniteElementSpace> fes0;
    std::shared_ptr<ParGridFunction> sol;
    std::shared_ptr<ParMesh> mesh;
    std::unique_ptr<ParGridFunction> entropy;

    std::unique_ptr<ModalBasis> modalBasis;
    const DenseMatrix &VDM;
    Array<int> max_order;
    DenseMatrix VDM_vol_ir;
    std::vector<DenseMatrix> VDM_face_irs;

    int dim, ndofs, num_equations, num_elements, order;

    Vector ffac;
    Vector nor;
    FaceElementTransformations *Tr_face;
    IntegrationPointTransformation Tr_ip;
    IntegrationRule mapped_face_ir;
    const IntegrationRule *vol_ir, *face_ir, *bdr_face_ir;


    const Table &element2element;
    const Table &element2bdrelement;
    
    Array<int> bdr_attr_marker;

    const FiniteElement *fe1, *fe2;
    const Element *element;
    const Element::Type etype, ftype;
    const Geometry::Type egeom, fgeom;

    Vector el_vdofs, ent_dof;
    Vector shape, vdof_col, modes_col;
    Vector state1, state2;
    DenseMatrix vdof_mat1, vdof_mat2, modes_mat;

    real_t el_min_rho, el_min_p, el_min_e, rho, p, e;
    
    real_t zeta, zeta1, zeta2, zeta3;
    DenseMatrix hierarch_states;

    const real_t rho_min = 1e-8;
    const real_t p_min = 1e-8;
    const real_t ent_tol = 1e-4;
    const real_t epsilon = 1e-10;
    const int N_iter = 20;
    const real_t zeta_tol = 1e-4;
    

    virtual void FilterModes(const DenseMatrix &hiearch_states, Vector &state, real_t zeta) override;
    void ComputeMinima(const Vector &state, const char *mode = "e");
    void ComputeElementMinima(const DenseMatrix &vdof_mat, const FiniteElement *fe, Vector &state, const char *mode = "e");
    void ComputeBdrStateMinima(const DenseMatrix &vdof_mat,
        const FiniteElement *fe, const Array<int> &el2bdrel, Vector &state1, Vector &state2, const char *mode = "e");
    void ComputeRoot(Vector &state, real_t &zeta);
    void ComputeHierarchStates(const DenseMatrix &VDM, DenseMatrix &hierarch_states, int row);
public:
    EntropyFilter(
        std::shared_ptr<ParFiniteElementSpace> vfes_,
        std::shared_ptr<ParFiniteElementSpace> fes0_,
        std::shared_ptr<ParGridFunction> sol_,
        std::shared_ptr<ParMesh> mesh_,
        std::unique_ptr<ParGridFunction> entropy_,
        std::unique_ptr<ModalBasis> modalBasis_,
        const Table &el2el, const Table &el2bdrel,
        const IntegrationRule *vol_ir_, const IntegrationRule *face_ir_,
        const IntegrationRule *bdr_face_ir_
        );
    void GetFilterConstraints(const Vector &x) override;
    void GetEntropyConstraints(const Vector &x);
    virtual void FilterSolution(Vector& sol) override;
};

}