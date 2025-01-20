#include "Indicator.hpp"

namespace Prandtl
{

Indicator::Indicator(std::shared_ptr<ParFiniteElementSpace> vfes,
                     std::shared_ptr<ParFiniteElementSpace> fes0,
                     std::shared_ptr<ParGridFunction> eta)
                    : vfes(vfes), fes0(fes0), eta(eta),
                      num_equations(vfes->GetVDim()), ndofs(vfes->GetFE(0)->GetDof()), order(vfes->GetElementOrder(0))
{
    state.SetSize(num_equations);
}

}