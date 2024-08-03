#include "BasicOperations.hpp"

namespace Prandtl
{

real_t ComputePressure(const Vector &state, const real_t gamma)
{
    const real_t &rho = state(0);
    const Vector rhoV(state.GetData() + 1, state.Size() - 2);
    const real_t &rhoE = state(state.Size() - 1);

    real_t ke = 0.0;
    for (int i = 0; i < state.Size() - 2; i++)
    {
        ke += rhoV(i) * rhoV(i);
    }
    ke /= rho;
    ke *= 0.5;

    return (gamma - 1.0) * (rhoE - ke);
}

real_t ComputeEntropy(const real_t rho, const real_t p, const real_t gamma)
{
#ifdef HYPERBOLIC
    return rho > 0.0 && p > 0.0 ? rho * log(p * pow(rho, -gamma)) : __DBL_MAX__;
#elif PARABOLIC
    return 0;
#endif
}

real_t ComputeInternalEnergy(const real_t p, const real_t rho, const real_t gamma, const real_t b)
{
    return p * (1.0 - b * rho) / rho / (gamma - 1.0);
}

real_t ComputeSoundSpeed(const real_t p, const real_t rho, const real_t gamma, const real_t b)
{
    return std::sqrt(gamma * p / (1.0 - b * rho) / rho);
}

real_t ComputeEnthalpy(const real_t p, const real_t rho, const real_t e)
{
    return e + p / rho;
}

real_t ComputeTotalEnthalpy(const Vector &state, const real_t gamma)
{
    const real_t &rho = state(0);
    const real_t &rhoE = state(state.Size() - 1);
    const real_t p = ComputePressure(state, gamma);

    return (rhoE + p) / rho;
}


Vector Conserv2Prim(const Vector &state, const real_t gamma)
{
    Vector prim_state(state.Size());
    prim_state(0) = state(0);
    prim_state(1) = state(1) / state(0);
    if (state.Size() > 3)
    {
        prim_state(2) = state(2) / state(0);
        if (state.Size() > 4)
        {
            prim_state(3) = state(3) / state(0);
        }
    }
    prim_state(state.Size() - 1) = ComputePressure(state, gamma);

    return prim_state;
}

Vector Prim2Conserv(const Vector &state, const real_t gamma)
{
    Vector conserv_state(state.Size());
    conserv_state(0) = state(0);
    conserv_state(1) = state(0) * state(1);
    if (state.Size() > 3)
    {
        conserv_state(2) = state(0) * state(2);
        if (state.Size() > 4)
        {
            conserv_state(3) = state(0) * state(3);
        }
    }
    double ke = 0.0;
    for (int i = 0; i < state.Size() - 2; i++)
    {
        ke += state(i + 1) * state(i + 1);
    }
    ke *= state(0);
    ke *= 0.5;
    conserv_state(state.Size() - 1) = state(state.Size() - 1) / (gamma - 1.0) + ke;

    return conserv_state;
}

void RotateState(Vector &state, const Vector &nor)
{
    MFEM_ASSERT(nor.Size() > 1, "Rotate only in 2D or 3D");
    MFEM_ASSERT(nor.Size() < 4, "Rotate only in 2D or 3D");

    Vector tan1(nor.Size());
    Vector rhoV(state.GetData() + 1, nor.Size());

    if (nor.Size() == 2)
    {
        tan1 = Normal2D(nor);
        real_t rho_u = rhoV * nor;
        real_t rho_v = rhoV * tan1; 
        rhoV(0) = rho_u;
        rhoV(1) = rho_v;
    }
    else
    {
        Vector tan2(nor.Size());
        Vector vec(3);
        if (std::abs(nor(0)) > 1e-12 || std::abs(nor(1)) > 1e-12)
        {
            vec(2) = 1.0;
        }
        else
        {
            vec(0) = 1.0;
        }
        tan1 = Cross(nor, vec);
        Normalize(tan1);
        tan2 = Cross(nor, tan1);
        Normalize(tan2);
        real_t rho_u = rhoV * nor;
        real_t rho_v = rhoV * tan1;
        real_t rho_w = rhoV * tan2;
        rhoV(0) = rho_u;
        rhoV(1) = rho_v;
        rhoV(2) = rho_w;
    }
}



Vector ComputeRoeAverage(const Vector &state1, const Vector &state2, const real_t gamma)
{
    Vector RoeAverage(state1.Size() + 1);
    real_t V2 = 0.0;

    const Vector prim1 = Conserv2Prim(state1, gamma);
    const Vector prim2 = Conserv2Prim(state2, gamma);

    const real_t &rho1 = prim1(0);
    const real_t &rho2 = prim2(0);
    RoeAverage(0) = std::sqrt(rho1 * rho2);

    const real_t &u1 = prim1(1);
    const real_t &u2 = prim2(1);
    RoeAverage(1) = (std::sqrt(rho1) * u1 + std::sqrt(rho2) * u2) / (std::sqrt(rho1) + std::sqrt(rho2));
    V2 += RoeAverage(1) * RoeAverage(1);


    if (prim1.Size() > 3)
    {
        const real_t &v1 = prim1(2);
        const real_t &v2 = prim2(2);
        RoeAverage(2) = (std::sqrt(rho1) * v1 + std::sqrt(rho2) * v2) / (std::sqrt(rho1) + std::sqrt(rho2));
        V2 += RoeAverage(2) * RoeAverage(2);
        if (prim1.Size() > 4)
        {
            const real_t &w1 = prim1(3);
            const real_t &w2 = prim2(3);
            RoeAverage(3) = (std::sqrt(rho1) * w1 + std::sqrt(rho2) * w2) / (std::sqrt(rho1) + std::sqrt(rho2));
            V2 += RoeAverage(3) * RoeAverage(3);
        }
    }


    const real_t H1 = ComputeTotalEnthalpy(state1, gamma);
    const real_t H2 = ComputeTotalEnthalpy(state2, gamma);
    RoeAverage(RoeAverage.Size() - 2) = (std::sqrt(rho1) * H1 + std::sqrt(rho2) * H2) / (std::sqrt(rho1) + std::sqrt(rho2));

    RoeAverage(RoeAverage.Size() - 1) = std::sqrt((gamma - 1.0) * (RoeAverage(RoeAverage.Size() - 2) -  0.5 * V2));

    return RoeAverage;
}

const Table& ElementIndextoBdrElementIndex(Mesh &mesh)
{
    Array<Connection> conn;
    conn.Reserve(mesh.GetNBE());
    int e, info;
    for (int i = 0; i < mesh.GetNBE(); i++)
    {
        mesh.GetBdrElementAdjacentElement(i, e, info);
        conn.Append(Connection(i, e));
    }
    conn.Sort();
    conn.Unique();

    Table *bdr_el2el = new Table(mesh.GetNBE(), conn);
    Table *el2bdr_el = new Table;
    Transpose(*bdr_el2el, *el2bdr_el, -1);

    bdr_el2el->~Table();
    return *el2bdr_el;
}

}

