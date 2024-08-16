#pragma once

#include "mfem.hpp"
#include "DGOperator.hpp"

namespace Prandtl
{

using namespace mfem;

/// Abstract class for solving systems of ODEs: dx/dt = f(x,t)
class DGODESolver
{
protected:
   /// Pointer to the associated DGOperator.
   DGOperator *f;  // f(.,t) : R^n --> R^n
   MemoryType mem_type;

public:
   DGODESolver() : f(NULL) { mem_type = Device::GetHostMemoryType(); }

   /// Associate a DGOperator with the ODE solver.
   /** This method has to be called:
       - Before the first call to Step().
       - When the dimensions of the associated DGOperator change.
       - When a time stepping sequence has to be restarted.
       - To change the associated DGOperator. */
   virtual void Init(DGOperator &f_);

   /** @brief Perform a time step from time @a t [in] to time @a t [out] based
       on the requested step size @a dt [in]. */
   /** @param[in,out] x   Approximate solution.
       @param[in,out] t   Time associated with the approximate solution @a x.
       @param[in,out] dt  Time step size.

       The following rules describe the common behavior of the method:
       - The input @a x [in] is the approximate solution for the input time
         @a t [in].
       - The input @a dt [in] is the desired time step size, defining the desired
         target time: t [target] = @a t [in] + @a dt [in].
       - The output @a x [out] is the approximate solution for the output time
         @a t [out].
       - The output @a dt [out] is the last time step taken by the method which
         may be smaller or larger than the input @a dt [in] value, e.g. because
         of time step control.
       - The method may perform more than one time step internally; in this case
         @a dt [out] is the last internal time step size.
       - The output value of @a t [out] may be smaller or larger than
         t [target], however, it is not smaller than @a t [in] + @a dt [out], if
         at least one internal time step was performed.
       - The value @a x [out] may be obtained by interpolation using internally
         stored data.
       - In some cases, the contents of @a x [in] may not be used, e.g. when
         @a x [out] from a previous Step() call was obtained by interpolation.
       - In consecutive calls to this method, the output @a t [out] of one
         Step() call has to be the same as the input @a t [in] to the next
         Step() call.
       - If the previous rule has to be broken, e.g. to restart a time stepping
         sequence, then the ODE solver must be re-initialized by calling Init()
         between the two Step() calls. */
   virtual void Step(Vector &x, real_t &t, real_t &dt) = 0;

   /// Perform time integration from time @a t [in] to time @a tf [in].
   /** @param[in,out] x   Approximate solution.
       @param[in,out] t   Time associated with the approximate solution @a x.
       @param[in,out] dt  Time step size.
       @param[in]     tf  Requested final time.

       The default implementation makes consecutive calls to Step() until
       reaching @a tf.
       The following rules describe the common behavior of the method:
       - The input @a x [in] is the approximate solution for the input time
         @a t [in].
       - The input @a dt [in] is the initial time step size.
       - The output @a dt [out] is the last time step taken by the method which
         may be smaller or larger than the input @a dt [in] value, e.g. because
         of time step control.
       - The output value of @a t [out] is not smaller than @a tf [in]. */
   virtual void Run(Vector &x, real_t &t, real_t &dt, real_t tf)
   {
      while (t < tf) { Step(x, t, dt); }
   }

   /// Function for getting and setting the state vectors
   virtual int   GetMaxStateSize() { return 0; }
   virtual int   GetStateSize() { return 0; }
   virtual const Vector &GetStateVector(int i)
   {
      mfem_error("DGODESolver has no state vectors");
      Vector *s = NULL; return *s; // Make some compiler happy
   }
   virtual void  GetStateVector(int i, Vector &state)
   {
      mfem_error("DGODESolver has no state vectors");
   }
   virtual void  SetStateVector(int i, Vector &state)
   {
      mfem_error("DGODESolver has no state vectors");
   }

   virtual ~DGODESolver() { }
};

}