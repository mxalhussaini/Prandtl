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
   real_t *x, *L, *Li, *Di;

   void ComputeUBDegs(Geometry::Type &gtype);
   void ComputeVDM(IntegrationRule &solpts);

public:
   int order, npts;

   ModalBasis(DG_FECollection &fec, Geometry::Type &gtype, int order, int dim);
   ~ModalBasis();

   void ComputeModes(const Vector &nodes);
   void ComputeModes(const Vector &nodes, Vector &modes);
   void SetModes(const Vector &modes);
   void GetModes(Vector &modes);
   real_t Eval(Vector &x);
   DenseMatrix ComputeVDM(const IntegrationRule *ir);
   Vector EvalGrad(Vector &x);
   void ComputeNodes(Vector &nodes);
   void ComputeNodes(const Vector& modes, Vector& nodes);
   const DenseMatrix& GetVandermonde();
   Array2D<int> GetPolyDegs();
};

}