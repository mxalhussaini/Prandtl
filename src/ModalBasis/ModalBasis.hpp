#pragma once

#include "mfem.hpp"

namespace Prandtl
{

using namespace mfem;

class ModalBasis
{
private:
   int dim;
   Array2D<int> ubdegs; // Array of modal basis degrees along each dimension
   Vector umc; // Vector of solution modal basis coefficients
   DenseMatrix V, V_inv; // (Inverse) Vandermonde matrix

   void ComputeUBDegs(Geometry::Type &gtype);
   void ComputeVDM(IntegrationRule &solpts);

public:
   int order, npts;

   ModalBasis(DG_FECollection &fec, Geometry::Type &gtype, int order, int dim);
   ~ModalBasis() = default;

   void SetSolution(Vector &u_elem);
   double Eval(Vector &x);
   void Eval(const Vector& x, Vector& vec);
   Vector EvalGrad(Vector &x);
   Vector GetModes();
   void SetModes(Vector& modes);
   void SetNodes(Vector &u_elem);
   void SetNodes(const Vector& modes, Vector& nodes);
   DenseMatrix GetVandermonde();
   Array2D<int> GetPolyDegs();
};

}