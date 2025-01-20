#include "RK3SSPExplicitSolver.hpp"

namespace Prandtl
{

void RK3SSPExplicitSolver::Init(DGOperator &f_)
{
   DGODESolver::Init(f_);
   int n = f->Width();
   y.SetSize(n, mem_type);
   k.SetSize(n, mem_type);
}

void RK3SSPExplicitSolver::Step(Vector &x, real_t &t, real_t &dt)
{
   // x0 = x, t0 = t, k0 = dt*f(t0, x0)
   f->SetTime(t);
#ifdef USE_FILTER
   f->filter->GetFilterConstraints(x);
#endif
   f->Mult(x, k);

   // x1 = x + k0, t1 = t + dt, k1 = dt*f(t1, x1)
   add(x, dt, k, y);
#ifdef USE_FILTER
   f->filter->FilterSolution(y);
   f->filter->GetFilterConstraints(y);
#endif
   f->SetTime(t + dt);
   f->Mult(y, k);


   // x2 = 3/4*x + 1/4*(x1 + k1), t2 = t + 1/2*dt, k2 = dt*f(t2, x2)
   y.Add(dt, k);
   add(3./4, x, 1./4, y, y);
#ifdef USE_FILTER
   f->filter->FilterSolution(y);
   f->filter->GetFilterConstraints(y);
#endif
   f->SetTime(t + dt/2);
   f->Mult(y, k);

   // x3 = 1/3*x + 2/3*(x2 + k2), t3 = t + dt
   y.Add(dt, k);
   add(1./3, x, 2./3, y, x);
#ifdef USE_FILTER
   f->filter->FilterSolution(x);
#endif
   t += dt;
}

}