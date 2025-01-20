#pragma once

#include "mfem.hpp"
#include "DGOperator.hpp"
#include "DGODESolver.hpp"

namespace Prandtl
{

using namespace mfem;

class Simulation
{
private:
    DGOperator *dg_operator;
    DGODESolver *ode_solver;
    
    Simulation();
    ~Simulation();

public:
    Simulation(const Simulation&) = delete;
    Simulation& operator = (const Simulation&) = delete;

    static Simulation& StartSimulation();
    void Run() const;
    void Setup();

};

}