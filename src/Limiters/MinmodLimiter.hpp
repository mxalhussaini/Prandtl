#include "Limiter.hpp"
#include "ModalBasis.hpp"

namespace Prandtl
{

class MinmodLimiter : public Limiter
{

private:
    IsoparametricTransformation Tr;
    InverseElementTransformation Tr_inv;
    ElementTransformation *Tr_nbr;

    std::shared_ptr<ParMesh> pmesh;
    std::shared_ptr<ParFiniteElementSpace> vfes;
    std::shared_ptr<ParFiniteElementSpace> vfes0;
    std::shared_ptr<ParGridFunction> sol_mean;

    int num_elements;
    int shifted_indx;

    const Table &element2element;
    Array<int> el2el;

    int geom;

    Vector center1, center2;
    IntegrationPoint extrap_center;

    Vector meanL, meanR, meanB, meanT, nbr_mean;
    Array<int> mean_indices;

    std::shared_ptr<ModalBasis> modalBasis;

public:
    MinmodLimiter(std::shared_ptr<ParFiniteElementSpace> vfes_,
        std::shared_ptr<ParMesh> pmesh_, const Table &element2element_);
    virtual void LimitSolution(Vector &x) override;
    void UpdateElementAverages(const Vector &x);

};

}