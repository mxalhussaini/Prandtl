#include "MinmodLimiter.hpp"

namespace Prandtl
{

MinmodLimiter::MinmodLimiter(std::shared_ptr<ParFiniteElementSpace> vfes_,
    std::shared_ptr<ParMesh> pmesh_, const Table &element2element_)
    : Limiter(), vfes(vfes_), pmesh(pmesh_), num_elements(pmesh->GetNE()),
      element2element(element2element_), geom(pmesh->GetElementBaseGeometry(0))
{
    Tr_inv.SetSolverType(InverseElementTransformation::Newton);
    center1.SetSize(pmesh->SpaceDimension());
    center2.SetSize(pmesh->SpaceDimension());
    meanL.SetSize(vfes->GetVDim());
    meanR.SetSize(vfes->GetVDim());
    meanB.SetSize(vfes->GetVDim());
    meanT.SetSize(vfes->GetVDim());
    nbr_mean.SetSize(vfes->GetVDim());
}


void MinmodLimiter::LimitSolution(Vector &x)
{
    for (int i = 0; i < num_elements; i++)
    {
        pmesh->GetElementTransformation(i, &Tr);
        pmesh->GetElementCenter(i, center1);
        Tr_inv.SetTransformation(Tr);
        element2element.GetRow(i, el2el);
        
        for (int j = 0; j < el2el.Size(); j++)
        {
            Tr_nbr = el2el[j] >= num_elements ?
            vfes->GetFaceNbrElementTransformation(el2el[j] - num_elements)
            : vfes->GetElementTransformation(el2el[j]);
            if (el2el[j] >= num_elements)
            {
                shifted_indx = el2el[j] - num_elements;
                Tr_nbr->Transform(Geometries.GetCenter(geom), center2);
                vfes0->GetFaceNbrElementVDofs(shifted_indx, mean_indices);
                sol_mean->FaceNbrData().GetSubVector(mean_indices, nbr_mean);
            }
            else
            {
                pmesh->GetElementCenter(el2el[i], center2);
                vfes0->GetElementVDofs(el2el[j], mean_indices);
                sol_mean->GetSubVector(mean_indices, nbr_mean);
            }
            Tr_inv.Transform(center2, extrap_center);
            if (std::abs(center1(0) - center2(0)) > std::abs(center1(1) - center2(1)))
            {
                if (center1(0) > center2(0))
                {
                    meanL = nbr_mean;
                }
                else
                {
                    meanR = nbr_mean;
                }
            }
            else
            {
                if (center1(1) > center2(1))
                {
                    meanB = nbr_mean;
                }
                else
                {
                    meanT = nbr_mean;
                }
            }
        }
        


    }
}

}
