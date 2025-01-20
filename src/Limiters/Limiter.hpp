#pragma once

#include "mfem.hpp"

namespace Prandtl
{

using namespace mfem;

class Limiter
{
public:
    Limiter() = default;
    virtual void LimitSolution(Vector &x) = 0;
    virtual ~Limiter() = default;
};

}