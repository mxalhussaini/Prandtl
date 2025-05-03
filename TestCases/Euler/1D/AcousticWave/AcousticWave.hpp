#pragma once

#include "ConditionFactory.hpp"

namespace Prandtl
{

// Acoustic Wave initial condition
std::function<void(const Vector&, Vector&)> AcousticWaveIC(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy, V;
        real_t a_inf, M_inf = 0.5, rho_inf = 1.225, p_inf = 1.0;
        a_inf = std::sqrt(gamma * p_inf / rho_inf);
        // MFEM_ASSERT(x.Size() == 1, "");


        // V = 1.0 / (2.0 * M_PI) * std::sin(2.0 * M_PI * x(0));

        // velocity_x = a_inf * (M_inf + 2.0 / ((gamma + 1.0) * a_inf) * V);
        // density = rho_inf * std::pow(1.0 + (gamma - 1.0) / 2.0 * (velocity_x / a_inf - M_inf), 2.0 / (gamma - 1.0));
        // pressure = p_inf * std::pow(1.0 + (gamma - 1.0) / 2.0 * (velocity_x / a_inf - M_inf), 2.0 * gamma / (gamma - 1.0));

        density = rho_inf + 0.2 * std::sin(2.0 * M_PI * x(0));
        velocity_x = a_inf * M_inf;
        pressure = p_inf;

        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

std::function<void(const Vector&, real_t, Vector&)> AcousticWaveExactSolution(real_t gamma)
{
    return [gamma](const Vector &x, real_t t, Vector &y)
    {
        real_t density, velocity_x, pressure, energy, V, v_;
        real_t a_inf, M_inf = 0.5, rho_inf = 1.225, p_inf = 1.0;
        a_inf = std::sqrt(gamma * p_inf / rho_inf);
        // MFEM_ASSERT(x.Size() == 1, "");

        // v_ = 0.0;
        // for (int i = 0; i < 200; i++)
        // {
        //     V = 1.0 / (2.0 * M_PI) * std::sin(2.0 * M_PI * (x(0) - v_ * t));
        //     if (std::abs(V - v_) < 1e-14)
        //     {
        //         break;
        //     }
        //     v_ = V;
        // }

        // velocity_x = a_inf * (M_inf + 2.0 / ((gamma + 1.0) * a_inf) * V);
        // density = rho_inf * std::pow(1.0 + (gamma - 1.0) / 2.0 * (velocity_x / a_inf - M_inf), 2.0 / (gamma - 1.0));
        // pressure = p_inf * std::pow(1.0 + (gamma - 1.0) / 2.0 * (velocity_x / a_inf - M_inf), 2.0 * gamma / (gamma - 1.0));

        density = rho_inf + 0.2 * std::sin(2.0 * M_PI * (x(0) - a_inf * M_inf * t));
        velocity_x = a_inf * M_inf;
        pressure = p_inf;
        
        energy = pressure / (gamma - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

// Registration helper that automatically registers these functions
struct RegisterAcousticWave
{
    RegisterAcousticWave()
    {
        // Register initial condition.
        ConditionFactory::Instance().RegisterInitialCondition1("AcousticWaveIC", AcousticWaveIC);

        // Register exact solution.
        ConditionFactory::Instance().RegisterVectorTDFunctionBoundaryCondition1("AcousticWaveExactSolution", AcousticWaveExactSolution);
    }
};
// Global static instance to ensure registration happens at startup.
static RegisterAcousticWave regAcousticWave;

}