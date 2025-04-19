#include "Simulation.hpp"
#include "ConditionFactory.hpp"

#include "LidDrivenCavity.hpp"
#include "DoubleMachReflection.hpp"
#include "KelvinHelmholtzInstability.hpp"

#include "json.hpp"

namespace Prandtl
{

Simulation& Simulation::SimulationCreate()
{
    static Simulation sim;
    return sim;
}

Simulation::Simulation()
{
    Mpi::Init();
    numProcs = Mpi::WorldSize();
    myRank = Mpi::WorldRank();
    Hypre::Init();


    if (Mpi::Root())
    {
        std::cout << "================================================" << std::endl;
        std::cout << "DGSEM Simulation Starting with " << numProcs << " Processors!" << std::endl;
        std::cout << "================================================" << std::endl;       
        
#ifdef PARABOLIC
        std::cout << "The Navier-Stokes Equations will be Solved!" << std::endl;
#else
        std::cout << "The Euler Equations will be Solved!" << std::endl;
#endif

#ifdef SUBCELL_FV_BLENDING
        std::cout << "Subcell Finite Volume Blending will be Activated!" << std::endl; 
#endif
    }
}

Simulation::~Simulation()
{
    if (Mpi::Root())
    {
        std::cout << "================================================" << std::endl;
        std::cout << "DGSEM Simulation Destroyed!" << std::endl;
        std::cout << "================================================" << std::endl;
    }
}

void Simulation::LoadConfig(const std::string &config_file_path)
{
    /*
    Load the configuration file and parse the settings
    */
    std::ifstream config_file(config_file_path);
    if (!config_file.is_open())
    {
        std::cerr << "Error Opening Configuration File at " << config_file_path << std::endl;
        return;
    }

    // Parse the JSON configuration file
    nlohmann::json config;
    config_file >> config;
    auto runtime = config["runTime"];

    order = runtime.value("order", 3);
    dim = runtime.value("dim", 2);
    num_equations = runtime.value("num_equations", 4);

    precision = runtime.value("precision", 15);
    std::cout.precision(precision);
    
    visualize = runtime["visualize"].get<bool>();
    if (visualize)
    {
        vis_steps = runtime.value("vis_steps", 100);
        paraview = runtime["paraview"].get<bool>();
        visit = runtime["visit"].get<bool>();
        output_file_path = runtime["output_file_path"].get<std::string>();
        if (!paraview && !visit)
        {
            std::cerr << "Error: Both ParaView and VisIt visualization options are disabled. Please choose at least one." << std::endl;
            return;
        }
    }

    nancheck = runtime["nancheck"].get<bool>();
    if (nancheck)
    {
        nancheck_steps = runtime.value("nancheck_steps", 1000);
    }

    clock_simulation = runtime["clock_simulation"].get<bool>();
    variable_dt = runtime["variable_dt"].get<bool>();
    dt = runtime.value("dt", 1e-4);
    t_final = runtime["final_time"].get<real_t>();

    physicsConstants = std::make_shared<PhysicsConstants>(
        runtime.value("gamma", 1.4),
        runtime.value("Pr", 0.72),
        runtime.value("R_gas", 287.05),
        runtime.value("mu", 0.02));


    if (runtime.contains("lifting_scheme"))
    {
        if (runtime["lifting_scheme"].get<std::string>() == "BR1")
        {
            liftingScheme = std::make_shared<LiftingBR1>();
        }
        else
        {
            std::cerr << "Error: Invalid lifting scheme specified." << std::endl;
            return;
        }
    }
    else
    {
        liftingScheme = nullptr;
    }

    std::string ode_solver_string = runtime["ode_solver"].get<std::string>();
    if (ode_solver_string == "ForwardEuler")
    {
        ode_solver = std::make_shared<ForwardEulerSolver>();
    }
    else if (ode_solver_string == "RK2")
    {
        ode_solver = std::make_shared<RK2Solver>();
    }
    else if (ode_solver_string == "RK3SSP")
    {
        ode_solver = std::make_shared<RK3SSPSolver>();
    }
    else if (ode_solver_string == "RK4")
    {
        ode_solver = std::make_shared<RK4Solver>();
    }
    else if (ode_solver_string == "RK6")
    {
        ode_solver = std::make_shared<RK6Solver>();
    }
    else if (ode_solver_string == "RK8")
    {
        ode_solver = std::make_shared<RK8Solver>();
    }
    else
    {
        std::cerr << "Error: Invalid ODE solver specified." << std::endl;
        return;
    }

    flux = std::make_shared<NavierStokesFlux>(dim,
        physicsConstants->gamma, physicsConstants->Pr, physicsConstants->mu, physicsConstants->mu0,
        physicsConstants->mu_bulk, physicsConstants->R_gas, physicsConstants->Ts, physicsConstants->T0);

    if (runtime["numerical_flux"].get<std::string>() == "Chandrashekar")
    {
        numericalFlux = std::make_shared<ChandrashekarFlux>(*flux, physicsConstants->gamma);
    }
    else
    {
        std::cerr << "Error: Invalid numerical flux specified." << std::endl;
        return;
    }

    signature = runtime["conditions"]["initial_conditions"].value("signature", 0);
    std::string IC_key = runtime["conditions"]["initial_conditions"].value("function", "LidDrivenCavityIC");

    if (signature == 0)
    {
        u0 = std::make_unique<VectorFunctionCoefficient>(num_equations, ConditionFactory::Instance().GetInitialCondition0(IC_key)());
    }
    else if (signature == 1)
    {
        real_t x1 = runtime["conditions"]["initial_conditions"]["params"].value("x1", 0.0);
        u0 = std::make_unique<VectorFunctionCoefficient>(num_equations, ConditionFactory::Instance().GetInitialCondition1(IC_key)(x1));
    }
    else if (signature == 2)
    {
        real_t x1 = runtime["conditions"]["initial_conditions"]["params"].value("x1", 0.0);
        real_t x2 = runtime["conditions"]["initial_conditions"]["params"].value("x2", 0.0);
        u0 = std::make_unique<VectorFunctionCoefficient>(num_equations, ConditionFactory::Instance().GetInitialCondition2(IC_key)(x1, x2));
    }
    else
    {
        std::cerr << "Error: Invalid initial condition signature." << std::endl;
        return;
    }

    Mesh *mesh;
    if (dim == 1)
    {
        
    }
    else
    {
        mesh = new Mesh(runtime["mesh_file"].get<std::string>());
    }

    if (runtime.contains("mesh_ordering"))
    {
        if (runtime["mesh_ordering"].get<std::string>() == "Hilbert")
        {
            mesh->GetHilbertElementOrdering(mesh_ordering);
        }
        else if (runtime["mesh_ordering"].get<std::string>() == "Gecko")
        {
            mesh->GetGeckoElementOrdering(mesh_ordering);
        }
        mesh->ReorderElements(mesh_ordering);
    }

    if (dim > 1)
    {
        mesh->EnsureNCMesh();
    }

    pmesh = std::make_shared<ParMesh>(MPI_COMM_WORLD, *mesh);
    mesh->Clear();

    if (runtime.contains("ref_levels"))
    {
        ref_levels = runtime.value("ref_levels", 0);
        for (int lev = 0; lev < ref_levels; lev++)
        {
            pmesh->UniformRefinement();
        }
    }

    fec = std::make_shared<DG_FECollection>(order, dim, btype);
    fec0 = std::make_shared<DG_FECollection>(0, dim);
    vfes = std::make_shared<ParFiniteElementSpace>(pmesh.get(), fec.get(), num_equations, ordering);
    fes0 = std::make_shared<ParFiniteElementSpace>(pmesh.get(), fec0.get());
    dfes = std::make_unique<ParFiniteElementSpace>(pmesh.get(), fec.get(), dim, ordering);
    fes = std::make_unique<ParFiniteElementSpace>(pmesh.get(), fec.get());

    sol = std::make_shared<ParGridFunction>(vfes.get());
    sol->ProjectCoefficient(*u0);

    eta = std::make_shared<ParGridFunction>(fes0.get());
    alpha = std::make_shared<ParGridFunction>(fes0.get());

    dudx = std::make_shared<ParGridFunction>(vfes.get());
    dudy = std::make_shared<ParGridFunction>(vfes.get());
    dudz = std::make_shared<ParGridFunction>(vfes.get());

    Geometry::Type gtype = vfes->GetFE(0)->GetGeomType();

    NS = std::make_unique<DGSEMOperator>(vfes, fes0, pmesh, eta, alpha, dudx, dudy, dudz,
        std::make_unique<Prandtl::DGSEMIntegrator>(pmesh, fes0, alpha, liftingScheme, *numericalFlux, order + 1),
        std::make_unique<Prandtl::PerssonPeraireIndicator>(vfes, fes0, eta, std::make_shared<Prandtl::ModalBasis>(*fec, gtype, order, dim), physicsConstants->gamma), physicsConstants->gamma);


    if (runtime["conditions"].contains("boundary_conditions"))
    {
        max_bdr_attr = pmesh->bdr_attributes.Max();  
        bdr_marker_vector.reserve(max_bdr_attr + 1);
        auto boundaries = runtime["conditions"]["boundary_conditions"];
        for (auto& boundary : boundaries.items())
        {
            std::string boundaryName = boundary.key();
            bdr_marker_vector.push_back(Array<int>(max_bdr_attr));
            set_marker = pmesh->bdr_attribute_sets.GetAttributeSetMarker(boundaryName);
            for (int b = 0; b < max_bdr_attr; b++)
            {
                if (set_marker[b])
                {
                    bdr_marker_vector.back()[b] = 1;
                }
            }

            auto bc_props = boundary.value();  // This is a JSON object.
            std::string type = bc_props["type"].get<std::string>();

            if (type == "slip")
            {
                NS->AddBdrFaceIntegrator(new SlipWallBdrFaceIntegrator(liftingScheme, *numericalFlux, order + 1, NS->GetTimeRef(), physicsConstants->gamma), bdr_marker_vector.back());
            }
            else if (type == "no-slip-adiabatic")
            {   
                if (bc_props["velocity"].contains("vector"))
                {
                    std::string velBC_key = bc_props["velocity"]["vector"].get<std::string>();
                    std::string heatBC_key = bc_props["heat"]["scalar"].get<std::string>();
                    NS->AddBdrFaceIntegrator(
                        new NoSlipAdiabWallBdrFaceIntegrator(liftingScheme, *numericalFlux, order + 1, NS->GetTimeRef(), physicsConstants->gamma,
                            ConditionFactory::Instance().GetScalarBoundaryCondition(heatBC_key),
                            ConditionFactory::Instance().GetVectorBoundaryCondition(velBC_key)), bdr_marker_vector.back());
                }
                else if (bc_props["velocity"].contains("function"))
                {
                    std::shared_ptr<VectorFunctionCoefficient> velBC;
                    std::shared_ptr<FunctionCoefficient> heatBC;

                    std::string velBC_key = bc_props["velocity"]["function"].get<std::string>();
                    std::string heatBC_key = bc_props["heat"]["function"].get<std::string>();

                    bool td;
                    if (bc_props["velocity"].contains("time_dependent"))
                    {
                        td = bc_props["velocity"]["time_dependent"].get<bool>();
                    }
                    else
                    {
                        td = false;
                    }

                    signature = bc_props["velocity"]["signature"].get<int>();
                    if (signature == 0)
                    {
                        velBC = std::make_shared<VectorFunctionCoefficient>(dim, ConditionFactory::Instance().GetVectorFunctionBoundaryCondition0(velBC_key)());
                    }
                    else if (signature == 1)
                    {
                        real_t x1 = bc_props["velocity"]["params"].value("x1", 0.0);
                        velBC = std::make_shared<VectorFunctionCoefficient>(dim, ConditionFactory::Instance().GetVectorFunctionBoundaryCondition1(velBC_key)(x1));
                    }
                    else if (signature == 2)
                    {
                        real_t x1 = bc_props["velocity"]["params"].value("x1", 0.0);
                        real_t x2 = bc_props["velocity"]["params"].value("x2", 0.0);
                        velBC = std::make_shared<VectorFunctionCoefficient>(dim, ConditionFactory::Instance().GetVectorFunctionBoundaryCondition2(velBC_key)(x1, x2));
                    }
                    else
                    {
                        std::cerr << "Error: Invalid boundary condition signature." << std::endl;
                        return;
                    }

                    signature = bc_props["heat"]["signature"].get<int>();
                    if (signature == 0)
                    {
                        heatBC = std::make_shared<FunctionCoefficient>(ConditionFactory::Instance().GetScalarFunctionBoundaryCondition0(heatBC_key)());
                    }
                    else if (signature == 1)
                    {
                        real_t x1 = bc_props["heat"]["params"].value("x1", 0.0);
                        heatBC = std::make_shared<FunctionCoefficient>(ConditionFactory::Instance().GetScalarFunctionBoundaryCondition1(heatBC_key)(x1));
                    }
                    else if (signature == 2)
                    {
                        real_t x1 = bc_props["heat"]["params"].value("x1", 0.0);
                        real_t x2 = bc_props["heat"]["params"].value("x2", 0.0);
                        heatBC = std::make_shared<FunctionCoefficient>(ConditionFactory::Instance().GetScalarFunctionBoundaryCondition2(heatBC_key)(x1, x2));
                    }
                    else
                    {
                        std::cerr << "Error: Invalid boundary condition signature." << std::endl;
                        return;
                    }

                    NS->AddBdrFaceIntegrator(
                        new NoSlipAdiabWallBdrFaceIntegrator(liftingScheme, *numericalFlux, order + 1, NS->GetTimeRef(), physicsConstants->gamma, *heatBC, *velBC, td), bdr_marker_vector.back());
                }
                else
                {
                    std::cerr << "Error: Invalid boundary condition type specified." << std::endl;
                    return;
                }
            }
            else if (type == "no-slip-isothermal")
            {

            }
            else if (type == "supersonic-outflow")
            {
                NS->AddBdrFaceIntegrator(new SupersonicOutflowBdrFaceIntegrator(liftingScheme, *numericalFlux, order + 1, NS->GetTimeRef(), physicsConstants->gamma), bdr_marker_vector.back());
            }
            else if (type == "specified-state")
            {
                if (bc_props.contains("vector"))
                {
                    std::string state_key = bc_props["vector"].get<std::string>();
                    NS->AddBdrFaceIntegrator(new SpecifiedStateBdrFaceIntegrator(liftingScheme, *numericalFlux, order + 1, NS->GetTimeRef(), physicsConstants->gamma,
                        ConditionFactory::Instance().GetVectorBoundaryCondition(state_key)), bdr_marker_vector.back());
                }
                else
                {
                    std::string state_key = bc_props["function"].get<std::string>();
                    signature = bc_props["signature"].get<int>();
                    std::shared_ptr<VectorFunctionCoefficient> stateBC;
                    bool td;
                    if (bc_props.contains("time_dependent"))
                    {
                        td = bc_props["time_dependent"].get<bool>();
                    }
                    else
                    {
                        td = false;
                    }

                    if (signature == 0)
                    {
                        stateBC = std::make_shared<VectorFunctionCoefficient>(num_equations, ConditionFactory::Instance().GetVectorTDFunctionBoundaryCondition0(state_key)());
                    }
                    else if (signature == 1)
                    {
                        real_t x1 = bc_props["params"].value("x1", 0.0);
                        stateBC = std::make_shared<VectorFunctionCoefficient>(num_equations, ConditionFactory::Instance().GetVectorTDFunctionBoundaryCondition1(state_key)(x1));
                    }
                    else if (signature == 2)
                    {
                        real_t x1 = bc_props["params"].value("x1", 0.0);
                        real_t x2 = bc_props["params"].value("x2", 0.0);
                        stateBC = std::make_shared<VectorFunctionCoefficient>(num_equations, ConditionFactory::Instance().GetVectorTDFunctionBoundaryCondition2(state_key)(x1, x2));
                    
                    }
                    else
                    {
                        std::cerr << "Error: Invalid boundary condition signature." << std::endl;
                        return;
                    }
                    NS->AddBdrFaceIntegrator(new SpecifiedStateBdrFaceIntegrator(liftingScheme, *numericalFlux, order + 1, NS->GetTimeRef(), physicsConstants->gamma, *stateBC, td), bdr_marker_vector.back());
                }
            }
            else
            {
                std::cerr << "Error: Invalid boundary condition type specified." << std::endl;
                return;
            }
        }
    }

    num_dofs_scalar = fes->GetNDofs();
    num_dofs_system = vfes->GetVSize();

    if (Mpi::Root())
    {
        std::cout << "The Number of Degrees of Freedom per Conservative Variable per Rank: " << num_dofs_scalar << std::endl;
        std::cout << "The Number of Degrees of Freedom (System) per Rank: " << num_dofs_system << std::endl;
    }

    NS->SetTime(t);
    ode_solver->Init(*NS);

    rho.MakeRef(fes.get(), *sol, 0);
    mom.MakeRef(dfes.get(), *sol, num_dofs_scalar);
    energy.MakeRef(fes.get(), *sol, (dim + 1) * num_dofs_scalar);

    u = std::make_unique<ParGridFunction>(fes.get());
    if (dim > 1)
    {
        v = std::make_unique<ParGridFunction>(fes.get()); 
        if (dim > 2)
        {
            w = std::make_unique<ParGridFunction>(fes.get());
        }
    }
    p = std::make_unique<ParGridFunction>(fes.get());
    
    if (visualize)
    {
        if (paraview)
        {
            pd = std::make_unique<ParaViewDataCollection>("ParaView", pmesh.get());
            pd->SetPrefixPath(output_file_path);
            pd->RegisterField("Density", &rho);
            pd->RegisterField("Horizontal V", u.get());
            if (dim > 1)
            {
                pd->RegisterField("Vertical V", v.get());
                if (dim > 2)
                {
                    pd->RegisterField("Normal V", w.get());
                }
            }
            pd->RegisterField("Pressure", p.get());
            pd->RegisterField("Blending Coeff", alpha.get());
            pd->SetLevelsOfDetail(order);
            pd->SetDataFormat(VTKFormat::BINARY);
            pd->SetHighOrderOutput(true);
        }
        else if (visit)
        {
            vd = std::make_unique<VisItDataCollection>("VisIt", pmesh.get());
            vd->SetPrefixPath(output_file_path);
            vd->SetPrecision(precision);
            vd->SetFormat(DataCollection::PARALLEL_FORMAT);

            vd->RegisterField("Density", &rho);
            vd->RegisterField("Horizontal V", u.get());
            if (dim > 1)
            {
                vd->RegisterField("Vertical V", v.get());
                if (dim > 2)
                {
                    vd->RegisterField("Normal V", w.get());
                }
            }
            vd->RegisterField("Pressure", p.get());
            vd->RegisterField("Blending Coeff", alpha.get());
        }
    }
}

void Simulation::Run()
{
    if (Mpi::Root())
    {
        std::cout << "================================================" << std::endl;
        std::cout << "DGSEM Simulation Running Now!" << std::endl;
        std::cout << "================================================" << std::endl;
    }

    // Get the minimum characteristic element size and compute the initial time step if time step is variable
    if (variable_dt && cfl > 0.0)
    {
        hmin = infinity();
        for (int i = 0; i < pmesh->GetNE(); i++)
        {
            hmin = std::min(pmesh->GetElementSize(i, 1), hmin);
        }
        MPI_Allreduce(MPI_IN_PLACE, &hmin, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, pmesh->GetComm());
        Vector z(sol->Size());
        NS->Mult(*sol, z);
        real_t max_char_speed = NS->GetMaxCharSpeed();
        MPI_Allreduce(MPI_IN_PLACE, &max_char_speed, 1,  MPITypeMap<real_t>::mpi_type, MPI_MAX, pmesh->GetComm());
        dt = cfl * hmin / (max_char_speed * std::pow(order + 1, 2));
    }

    // Clock the simulation?
    if (clock_simulation)
    {
        tic_toc.Clear();
        tic_toc.Start();
    }

    // Visualize the initial condition?
    if (visualize)
    {
        for (int i = 0; i < num_dofs_scalar; i++)
        {
            (*u)(i) = mom(i) / rho(i);
            V_sq = (*u)(i) * (*u)(i);
            if (dim > 1)
            {
                (*v)(i) = mom(i + num_dofs_scalar) / rho(i);
                V_sq += (*v)(i) * (*v)(i);
                if (dim > 2)
                {
                    (*w)(i) = mom(i + 2 * num_dofs_scalar) / rho(i);
                    V_sq += (*w)(i) * (*w)(i);
                }
            }
            (*p)(i) = physicsConstants->gammaM1 * (energy(i) - 0.5 * rho(i) * V_sq);
        }

        if (paraview)
        {
            pd->SetCycle(0);
            pd->SetTime(0.0);
            pd->Save();
        }
        else if (visit)
        {
            vd->SetCycle(0);
            vd->SetTime(0.0);
            vd->Save();
        }
    }

    // Start the time-stepping loop
    for (int ti = 0; !done;)
    {
        // Compute the time step size
        dt_real = std::min(dt, t_final - t);

        // Perform the time step
        ode_solver->Step(*sol, t, dt_real);

        // Update the time step size with CFL?
        if (variable_dt && cfl > 0)
        {
            real_t max_char_speed = NS->GetMaxCharSpeed();
            MPI_Allreduce(MPI_IN_PLACE, &max_char_speed, 1, MPITypeMap<real_t>::mpi_type, MPI_MAX, pmesh->GetComm());
            dt = cfl * hmin / (max_char_speed * std::pow(order + 1, 2));
        }
        ti++;

        // Check for completion
        done = (t >= t_final - 1e-8 * dt);

        // Check for NaN/Inf values?
        if (nancheck && ti % nancheck_steps == 0)
        {
            for (const real_t &val : rho)
            {
                if (std::isnan(val) || std::isinf(val))
                {
                    std::cerr << "NaN/Inf Detected at Time Step " << ti << " on Rank " << myRank << "!" << std::endl;
                    break;
                }
            }
        }

        // Visualize the solution?
        if (visualize && ti % vis_steps == 0)
        {
            for (int i = 0; i < num_dofs_scalar; i++)
            {
                (*u)(i) = mom(i) / rho(i);
                V_sq = (*u)(i) * (*u)(i);
                if (dim > 1)
                {
                    (*v)(i) = mom(i + num_dofs_scalar) / rho(i);
                    V_sq += (*v)(i) * (*v)(i);
                    if (dim > 2)
                    {
                        (*w)(i) = mom(i + 2 * num_dofs_scalar) / rho(i);
                        V_sq += (*w)(i) * (*w)(i);
                    }
                }
                (*p)(i) = physicsConstants->gammaM1 * (energy(i) - 0.5 * rho(i) * V_sq);
            }
    
            if (paraview)
            {
                pd->SetCycle(ti);
                pd->SetTime(t);
                pd->Save();
            }
            else if (visit)
            {
                vd->SetCycle(ti);
                vd->SetTime(t);
                vd->Save();
            }

            if (Mpi::Root())
            {
                std::cout << "time step: " << ti << ", time: " << t << "\n";
            }
        }
    }

    // Stop the clock if enabled
    if (clock_simulation)
    {
        tic_toc.Stop();
        if (Mpi::Root())
        {
            std::cout << "================================================" << std::endl;
            std::cout << "DGSEM Simulation Completed in " << tic_toc.RealTime() << " seconds!" << std::endl;
            std::cout << "================================================" << std::endl;
        }
    }
    else
    {
        if (Mpi::Root())
        {
            std::cout << "================================================" << std::endl;
            std::cout << "DGSEM Simulation Completed!" << std::endl;
            std::cout << "================================================" << std::endl;
        }
    }
}

}