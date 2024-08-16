#include "DGODESolver.hpp"

namespace Prandtl
{

void DGODESolver::Init(DGOperator &f_)
{
   this->f = &f_;
   mem_type = GetMemoryType(f_.GetMemoryClass());
}

}