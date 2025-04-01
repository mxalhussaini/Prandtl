#include "LiftingScheme.hpp"

namespace Prandtl
{

void LiftingScheme::SetLiftingParameters(const IntegrationRule *ir_, const IntegrationRule *ir_face_, const IntegrationRule *ir_vol_, const DenseMatrix &D_T_, int num_equations_, int Np_, int dim)
{
    ir = ir_;
    ir_face = ir_face_;
    ir_vol = ir_vol_;
    D_T = D_T_;
    dof = dof1 = dof2 = ir_vol->GetNPoints();
    num_equations = num_equations_;
    Np_x = Np_;
    Np_y = dim > 1 ? Np_x : 1;
    Np_z = dim > 2 ? Np_x : 1;
    
    grad_state.SetSize(num_equations);
    grad_mat1.SetSize(num_equations, dim);
    grad_mat2.SetSize(num_equations, dim);

    adj.SetSize(dim);
    nor.SetSize(dim);

    shape1.SetSize(ir_vol->GetNPoints());
    shape2.SetSize(ir_vol->GetNPoints());

    state1.SetSize(num_equations);
    state2.SetSize(num_equations);

    f.SetSize(num_equations);
    g.SetSize(num_equations);
    h.SetSize(num_equations);

    dU.SetSize(num_equations);

    el_dudxi.SetSize(num_equations);
    el_dudeta.SetSize(num_equations);
    el_dudzeta.SetSize(num_equations);

}

}