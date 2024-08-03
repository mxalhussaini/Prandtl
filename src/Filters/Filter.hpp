#pragma once

#include "mfem.hpp"

namespace Prandtl
{
using namespace mfem;

class Filter
{
protected:
    virtual void FilterModes(const DenseMatrix &hiearch_states, Vector &state, real_t zeta) const = 0;
public:
    Filter() = default;

    virtual void FilterSolution(Vector& sol) = 0;
    
    virtual ~Filter() = default;
};
}