#pragma once

#include "mfem.hpp"
#include "BdrFaceIntegrator.hpp"

namespace Prandtl
{
    
using namespace mfem;

class Filter
{
protected:
    Array<BdrFaceIntegrator*> bfnfi;
    Array<Array<int>*> bfnfi_marker;

    virtual void FilterModes(const DenseMatrix &hiearch_states, Vector &state, real_t zeta) = 0;
public:
    Filter() = default;

    virtual void FilterSolution(Vector& sol) = 0;

    inline void AddBdrFaceIntegrator(BdrFaceIntegrator* bfi, Array<int>& bdr_marker)
    {
        bfnfi_marker.Append(&bdr_marker);
        bfnfi.Append(bfi);
    }
    virtual void GetFilterConstraints(const Vector &x) = 0;
    virtual ~Filter() = default;
};
}