#pragma once

#include "mfem.hpp"
#include "Filter.hpp"

namespace Prandtl
{

using namespace mfem;

class EntropyFilter : public Filter
{
private:
    // mutable GridFunction zeta;
    mutable GridFunction entropy;
    ParGridFunction &sol;
    const ParFiniteElementSpace &vfes;
    const FiniteElementSpace &fes0;
    const Table &element2element;
    const Table &element2bdrelement;
    const ParMesh &pmesh;
    const NonlinearForm &nonlinearForm;
    const Array<BdrFaceIntegrator*> &bfnfi;
    const Array<Array<int>*> &bfnfi_marker;

    real_t el_min_e;

    const FiniteElement *fe1, *fe2;
    ElementTransformation *Tr;
    IntegrationPointTransformation *Tr_ip;
    FaceElementTransformations *Tr_face;
    Vector phys_ip;
    Array<int> bdr_attr_marker;
    const Element *element;
    const Element::Type etype, ftype;
    const Geometry::Type egeom, fgeom;
    int IntOrderOffset = 3;
    IntegrationRule *vol_ir, *face_ir, *bdr_face_ir, *mapped_face_ir;
    Vector el_vdofs;
    DenseMatrix vdof_mat1, vdof_mat2;
    real_t rho, p, e, gamma;

    Vector state1, state2, shape;

    int dim, ndofs, num_equations, num_elements;

    Vector nor;

    virtual void FilterModes(const DenseMatrix &hiearch_states, Vector &state, real_t zeta) const override;
    void ComputeMinEntropy(const Vector &state, real_t el_min_e);
public:
    EntropyFilter();
    void GetEntropyConstraints(const Vector &x);
    virtual void FilterSolution(Vector& sol) override;
};

}