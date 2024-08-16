#pragma once

#include "mfem.hpp"
// #define HYPERBOLIC

namespace Prandtl
{
    using namespace mfem;

    real_t ComputePressure(const Vector &state, const real_t gamma);
    real_t ComputeEntropy(const real_t rho, const real_t p, const real_t gamma);
    real_t ComputeInternalEnergy(const real_t p, const real_t rho, const real_t gamma, const real_t b = 0.0);
    real_t ComputeSoundSpeed(const real_t p, const real_t rho, const real_t gamma, const real_t b = 0.0);
    real_t ComputeEnthalpy(const real_t p, const real_t rho, const real_t e);
    real_t ComputeTotalEnthalpy(const Vector &state, const real_t gamma);

    Vector Conserv2Prim(const Vector &state, const real_t gamma);
    Vector Prim2Conserv(const Vector &state, const real_t gamma);

    inline void Normalize(Vector &vec)
    {
        vec /= vec.Norml2();
    }

    inline Vector Cross(const Vector &vec1, const Vector &vec2)
    {
        MFEM_ASSERT(vec1.Size() == 3 && vec2.Size() = 3, "Apply cross product only to 3D vectors");
        Vector cross(3);

        cross(0) = vec1(1) * vec2(2) - vec1(2) * vec2(1);
        cross(1) = vec1(2) * vec2(0) - vec1(0) * vec2(2);
        cross(2) = vec1(0) * vec2(1) - vec1(1) * vec2(0);

        return cross;
    }

    inline Vector Normal(const Vector &vec)
    {
        MFEM_ASSERT(vec.Size() == 2 || vec.Size() == 3, "Apply only to 2D/3D vectors");
        
        Vector nor(vec.Size());

        nor(0) = vec.Size() == 2? -vec(1) :
            (std::abs(vec(0)) > 1e-10 ? -(vec(1) / vec(0)) / std::sqrt(1.0 + std::pow(vec(1) / vec(0), 2)) : -1.0);
        nor(1) = vec.Size() == 2? vec(0) : 
            (std::abs(vec(0)) > 1e-10 ? 1.0 / std::sqrt(1.0 + std::pow(vec(1) / vec(0), 2)) : 0.0);
        if (vec.Size() == 3)
        {
            nor(2) = 0.0;
        }
        return nor;
    }


    void RotateState(Vector &state, const Vector &nor);
    void RotateBack(Vector &state, const Vector &nor);

    Vector ComputeRoeAverage(const Vector &state1, const Vector &state2, const real_t gamma);

    const Table& ElementIndextoBdrElementIndex(Mesh &mesh);

    
}