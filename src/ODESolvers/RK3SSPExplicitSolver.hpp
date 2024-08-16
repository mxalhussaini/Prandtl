#pragma once

#include "mfem.hpp"
#include "DGODESolver.hpp"

namespace Prandtl
{

using namespace mfem;

/// Third-order, strong stability preserving (SSP) Runge-Kutta method
class RK3SSPExplicitSolver : public DGODESolver
{
private:
   Vector y, k;

public:
   void Init(DGOperator &f_) override;

   void Step(Vector &x, real_t &t, real_t &dt) override;
};

}