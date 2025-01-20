#include "Simulation.hpp"

namespace Prandtl
{

Simulation::Simulation()
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "DG Simulation has started" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
}

Simulation::~Simulation()
{
    std::cout << "------------------------------------------------" << std::endl;
    std::cout << "DG Simulation has ended" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
}

void Simulation::Setup()
{
    /*
    Parse the JSON configuration file with the mesh,
    initial and boundary conditions, and simulation settings
    */
}

void Simulation::Run() const
{
    std::cout << "Running time steps" << std::endl;
}

Simulation& Simulation::StartSimulation()
{
    static Simulation sim;
    return sim;
}


}