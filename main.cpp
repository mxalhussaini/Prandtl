#include "Prandtl.hpp"

using namespace mfem;

std::function<void(const Vector&, Vector&)> Cylinder(real_t Ma, real_t Re, real_t gamma, real_t mu, real_t D);
std::function<void(const Vector&, Vector&)> LidDrivenCavityIC(real_t Ma, real_t gamma);
std::function<real_t(const Vector&)> LidDrivenCavityAdiaBC();
std::function<void(const Vector&, Vector &)> LidDrivenCavityVelBC();
std::function<void(const Vector&, Vector&)> LidDrivenCavityVelLidBC();
std::function<void(const Vector&, Vector&)> Sod1D(real_t gamma);
std::function<void(const Vector&, Vector&)> quiescent_if(int dim);
std::function<void(const Vector&, real_t, Vector&)> dmr_top_bf(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> dmr_if(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> triple_point_if(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> KelvinHelmholtzInstability_if(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> airfoil_if(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> ramp_if(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> GetMovingVortexInit(const real_t radius, const real_t Minf, const real_t beta,
   const real_t gas_constant, const real_t specific_heat_ratio);

const real_t Prandtl::gamma = 1.4;
const real_t Prandtl::gammaInverse = 1.0 / gamma;
const real_t Prandtl::gammaM1 = Prandtl::gamma - 1.0;
const real_t Prandtl::gammaP1 = Prandtl::gamma + 1.0;
const real_t Prandtl::gammaM1Inverse = 1.0 / Prandtl::gammaM1;
const real_t Prandtl::gammaP1Inverse = 1.0 / Prandtl::gammaP1;
const real_t Prandtl::gammaM1_gammaInverse = Prandtl::gammaM1 * Prandtl::gammaM1Inverse;
const real_t Prandtl::gamma_gammaM1Inverse = Prandtl::gamma * Prandtl::gammaM1Inverse;
const real_t Prandtl::R_gas = 287.05;
const real_t Prandtl::Pr = 0.72;
const real_t Prandtl::PrInverse = 1.0 / Pr;
const real_t Prandtl::cp = Prandtl::gamma_gammaM1Inverse * Prandtl::R_gas;



int main (int argc, char* argv[])
{
    Mpi::Init(argc, argv);
    const int numProcs = Mpi::WorldSize();
    const int myRank = Mpi::WorldRank();
    Hypre::Init();

    int order = 3;
    real_t t_final = 3.0;
    real_t dt = -0.01;
    real_t cfl = 1.0;
    real_t mu = 0.02;
    real_t Ma = 0.1;
    real_t Re = 20.0;
    bool visualization = true;
    int vis_steps = 100;
    int nancheck_steps = 100;
    int precision = 15;
    std::cout.precision(precision);

    
    // std::string filename = "../../periodic-segment.mesh";
    // std::string filename = "../../inline-segment.mesh";
    // std::string filename = "../../periodic-cube.msh";
    // std::string filename = "../../periodic-square.mesh";
    // std::string filename = "../../SupercriticalAirfoilQuad.msh";
    // std::string filename = "../../CompressionRamp.msh";
    std::string filename = "../../cylinder.msh";
    Mesh *mesh = new Mesh(filename);
    // Mesh *mesh = new Mesh();
    // *mesh = Mesh::MakeCartesian1D(400, 1.0);
    // *mesh = Mesh::MakeCartesian2D(16, 16, Element::QUADRILATERAL, false, 2, 2, true);
    // *mesh = Mesh::MakeCartesian2D(384, 96, Element::QUADRILATERAL, false, 4.0, 1.0, true);
    // *mesh = Mesh::MakeCartesian2D(512, 384, Element::QUADRILATERAL, false, 4.0, 3.0, true);

    
    int dim = mesh->Dimension();
    int num_equations = dim + 2;

    if (dim > 1) mesh->EnsureNCMesh();
    std::shared_ptr<ParMesh> pmesh = std::make_shared<ParMesh>(MPI_COMM_WORLD, *mesh);
    mesh->Clear();

    // if (Mpi::Root())
    // {
    //     pmesh->bdr_attributes.Print(std::cout);
    //     pmesh->bdr_attribute_sets.Print(std::cout);
    // }

    // for (int k = 0; k < pmesh->GetNBE(); k++)
    // {
    //     if (pmesh->GetBdrAttribute(k) == 134)
    //     // if (pmesh->GetBdrAttribute(k) == 7)
    //     {
    //         pmesh->SetBdrAttribute(k, 1);
    //     }
    //     else if (pmesh->GetBdrAttribute(k) == 135)
    //     // else if (pmesh->GetBdrAttribute(k) == 10 || pmesh->GetBdrAttribute(k) == 8)
    //     {
    //         pmesh->SetBdrAttribute(k, 2);
    //     }
    //     else
    //     {
    //         pmesh->SetBdrAttribute(k, 3);
    //     }
    // }
    // pmesh->SetAttributes();

    // pmesh->PrintBdrVTU("cylinder_pmesh_bdr");

    for (int lev = 0; lev < 0; lev++)
    {
        pmesh->UniformRefinement();
    }

    DG_FECollection fec(order, dim, BasisType::GaussLobatto);
    DG_FECollection fec0(0, dim);
    std::shared_ptr<ParFiniteElementSpace> vfes = std::make_shared<ParFiniteElementSpace>(pmesh.get(), &fec, num_equations, Ordering::byNODES);
    std::shared_ptr<ParFiniteElementSpace> fes0 = std::make_shared<ParFiniteElementSpace>(pmesh.get(), &fec0);
    ParFiniteElementSpace* dfes =  new ParFiniteElementSpace(pmesh.get(), &fec, dim, Ordering::byNODES);
    ParFiniteElementSpace* fes = new ParFiniteElementSpace(pmesh.get(), &fec);

    if (Mpi::Root())
    {
        std::cout << "Number of unknowns: " << vfes->GetVSize() << std::endl;
    }
    
    // VectorFunctionCoefficient u0(num_equations, LidDrivenCavityIC(Ma, Prandtl::gamma));
    VectorFunctionCoefficient u0(num_equations, Cylinder(Ma, Re, Prandtl::gamma, mu, 2.0));

    std::shared_ptr<ParGridFunction> sol = std::make_shared<ParGridFunction>(vfes.get());
    sol->ProjectCoefficient(u0);

    std::shared_ptr<ParGridFunction> eta = std::make_shared<ParGridFunction>(fes0.get());
    std::shared_ptr<ParGridFunction> alpha = std::make_shared<ParGridFunction>(fes0.get());
    std::shared_ptr<ParGridFunction> grad_x = std::make_shared<ParGridFunction>(vfes.get());
    std::shared_ptr<ParGridFunction> grad_y = std::make_shared<ParGridFunction>(vfes.get());
    std::shared_ptr<ParGridFunction> grad_z = std::make_shared<ParGridFunction>(vfes.get());

    Prandtl::NavierStokesFlux flux(dim, mu);
    Prandtl::ChandrashekarFlux numericalFlux(flux);

    IntegrationRules GLIntRules(0, Quadrature1D::GaussLobatto);
    Geometry::Type gtype = vfes->GetFE(0)->GetGeomType();

    std::vector<Prandtl::BdrFaceIntegrator*> bfnfi;
    std::vector<Array<int>> bdr_marker;

    Prandtl::DGSEMOperator NS(vfes, fes0, pmesh, eta, alpha,
            grad_x, grad_y, grad_z,
            std::make_unique<Prandtl::DGSEMIntegrator>(pmesh, vfes, fes0, grad_x, grad_y, grad_z, alpha, numericalFlux, order + 1),
            std::make_unique<Prandtl::PerssonPeraireIndicator>(vfes, fes0, eta, std::make_shared<Prandtl::ModalBasis>(fec, gtype, order, dim)),
            bfnfi, bdr_marker);

    // const IntegrationRule *ir_face = &GLIntRules.Get(Geometry::SEGMENT, 2 * order - 1);

    // VectorFunctionCoefficient vn(dim, LidDrivenCavityVelBC());
    // VectorFunctionCoefficient vlid(dim, LidDrivenCavityVelLidBC());
    // FunctionCoefficient qn(LidDrivenCavityAdiaBC());
    // Array<int> walls(pmesh->bdr_attributes.Size());
    // walls = 1; walls[2] = 0;
    // Array<int> lid(pmesh->bdr_attributes.Size());
    // lid = 0; lid[2] = 1;
    // NS.AddBdrFaceIntegrator(new Prandtl::NoSlipAdiabWallBdrFaceIntegrator(numericalFlux, order + 1,
    //     grad_x, grad_y, nullptr, vfes, NS.GetTimeRef(), qn, vn), walls);
    // NS.AddBdrFaceIntegrator(new Prandtl::NoSlipAdiabWallBdrFaceIntegrator(numericalFlux, order + 1,
    //     grad_x, grad_y, nullptr, vfes, NS.GetTimeRef(), qn, vlid), lid);

    VectorFunctionCoefficient vn(dim, LidDrivenCavityVelBC());
    FunctionCoefficient qn(LidDrivenCavityAdiaBC());
    Array<int> box(pmesh->bdr_attributes.Size());
    box = 1; box[pmesh->bdr_attributes.Size() - 1] = 0;
    Array<int> circle(pmesh->bdr_attributes.Size());
    circle = 0; circle[pmesh->bdr_attributes.Size() - 1] = 1;
    NS.AddBdrFaceIntegrator(new Prandtl::NoSlipAdiabWallBdrFaceIntegrator(numericalFlux, order + 1,
        grad_x, grad_y, nullptr, vfes, NS.GetTimeRef(), qn, vn), circle);
    // NS.AddBdrFaceIntegrator(new Prandtl::SpecifiedStateBdrFaceIntegrator(numericalFlux, order + 1,
    //     grad_x, grad_y, nullptr, vfes, NS.GetTimeRef(), u0), box);

    Vector state(num_equations);
    ElementTransformation *Tr = vfes->GetElementTransformation(0);
    IntegrationPoint ip;
    ip.x = 0.0;
    if (dim > 1)
    {
        ip.y = 0.0;
        if (dim > 2)
        {
            ip.z = 0.0;
        }
    }
    u0.Eval(state, *Tr, ip); 

    NS.AddBdrFaceIntegrator(new Prandtl::SpecifiedStateBdrFaceIntegrator(numericalFlux, order + 1,
        grad_x, grad_y, nullptr, vfes, NS.GetTimeRef(), state), box);
    

    ParGridFunction rho, mom, energy;
    rho.MakeRef(fes, *sol, 0);
    mom.MakeRef(dfes, *sol, fes->GetNDofs());
    energy.MakeRef(fes, *sol, (dim + 1) * fes->GetNDofs());

    ParGridFunction vel(dfes);
    ParGridFunction vel_x, vel_y(fes);
    ParGridFunction pressure(fes);
    ParGridFunction tmp(fes);

    vel = mom;
    vel_x.MakeRef(fes, vel, 0);
    vel_x /= rho;
    if (dim == 2)
    {
        vel_y.MakeRef(fes, vel, fes->GetNDofs());
        vel_y /= rho;
    }
    tmp = vel_x;
    tmp *= vel_x;
    pressure = 0.0;
    pressure += tmp;
    if (dim == 2)
    {
        tmp = vel_y;
        tmp *= vel_y;
        pressure += tmp;
    }
    pressure *= 0.5;
    pressure *= rho;
    pressure.Neg();
    pressure += energy;
    pressure *= Prandtl::gammaM1;

    ParaViewDataCollection *pd = new ParaViewDataCollection("Cylinder", pmesh.get());
    pd->SetPrefixPath("ParaView");
    pd->RegisterField("density", &rho);
    pd->RegisterField("u", &vel_x);
    if (dim > 1)
    {
        pd->RegisterField("v", &vel_y);
    }
    pd->RegisterField("pressure", &pressure);
    pd->RegisterField("blending coeff", alpha.get());
    // pd->RegisterField("energy", &energy);
    // pd->RegisterField("momentum", &mom);
    pd->SetLevelsOfDetail(order);
    pd->SetDataFormat(VTKFormat::BINARY);
    pd->SetHighOrderOutput(true);
    pd->SetCycle(0);
    pd->SetTime(0.0);
    pd->Save();

    socketstream sout;
    if (visualization)
    {
        char vishost[] = "localhost";
        int visport = 19916;

        sout.open(vishost, visport);
        if (!sout)
        {
            visualization = false;
            if (Mpi::Root())
            {
            std::cout << "Unable to connect to GLVis server at " << vishost << ':'
                    << visport << std::endl;
            std::cout << "GLVis visualization disabled.\n";
            }
        }
        else
        {
            sout.precision(precision);
            sout << "parallel " << numProcs << " " << myRank << "\n";
            sout << "solution\n" << *pmesh << pressure;
            sout << "window_title 'density, t = 0'\n";
            // sout << "view 0 0\n";  // view from top
            // sout << "keys jlm\n";  // turn off 
            sout << "pause\n";
            sout << std::flush;
            if (Mpi::Root())
            {
                std::cout << "GLVis visualization paused."
                        << " Press space (in the GLVis window) to resume it.\n";
            }
            MPI_Barrier(pmesh->GetComm());
        }
    }


    real_t hmin = infinity();
    // if (cfl > 0)
    // {
        // for (int i = 0; i < pmesh->GetNE(); i++)
        // {
        //     hmin = std::min(pmesh->GetElementSize(i, 1), hmin);
        // }
        // MPI_Allreduce(MPI_IN_PLACE, &hmin, 1, MPITypeMap<real_t>::mpi_type, MPI_MIN, pmesh->GetComm());

        // Vector z(sol->Size());
        // NS.Mult(*sol, z);
        // real_t max_char_speed = NS.GetMaxCharSpeed();
        // MPI_Allreduce(MPI_IN_PLACE, &max_char_speed, 1,  MPITypeMap<real_t>::mpi_type,
        //             MPI_MAX,
        //             pmesh->GetComm());
        // dt = cfl * hmin / (max_char_speed * std::pow(order + 1, 2));
    // }

    dt = 1e-4;

    tic_toc.Clear();
    tic_toc.Start();

    ODESolver *ode_solver = new RK3SSPSolver();
    real_t t = 0.0;
    NS.SetTime(t);
    ode_solver->Init(NS);

    bool done = false;
    for (int ti = 0; !done;)
    {
        real_t dt_real = std::min(dt, t_final - t);

        ode_solver->Step(*sol, t, dt_real);
        // if (cfl > 0) // update time step size with CFL
        // {
        //     real_t max_char_speed = NS.GetMaxCharSpeed();
        //     MPI_Allreduce(MPI_IN_PLACE, &max_char_speed, 1,  MPITypeMap<real_t>::mpi_type,
        //                 MPI_MAX,
        //                 pmesh->GetComm());
        //     dt = cfl * hmin / (max_char_speed * std::pow(order + 1, 2));
        // }
        ti++;

        done = (t >= t_final - 1e-8 * dt);
        if (done || ti % nancheck_steps == 0)
        {
            for (auto el: *sol)
            {
                if (std::isnan(el) || std::isinf(el))
                {
                    MFEM_ABORT("nan/inf encountered");
                }
            }
        }
        if (done || ti % vis_steps == 0)
        {
            vel = mom;
            vel_x.MakeRef(fes, vel, 0);
            vel_x /= rho;
            if (dim == 2)
            {
                vel_y.MakeRef(fes, vel, fes->GetNDofs());
                vel_y /= rho;
            }
            tmp = vel_x;
            tmp *= vel_x;
            pressure = 0.0;
            pressure += tmp;
            if (dim == 2)
            {
                tmp = vel_y;
                tmp *= vel_y;
                pressure += tmp;
            }
            pressure *= 0.5;
            pressure *= rho;
            pressure.Neg();
            pressure += energy;
            pressure *= Prandtl::gammaM1;
            pd->SetCycle(ti);
            pd->SetTime(t);
            pd->Save();

            if (visualization)
            {
                sout << "window_title 'density, t = " << t << "'\n";
                sout << "parallel " << numProcs << " " << myRank << "\n";
                sout << "solution\n" << *pmesh << pressure << std::flush;
            }

            if (Mpi::Root())
            {
                // std::cout << "time step: " << ti << ", time: " << t << ", cL: " << lift << std::endl;
                std::cout << "time step: " << ti << ", time: " << t << std::endl;
            }
        }
        // break;
    }

    tic_toc.Stop();

    if (Mpi::Root())
    {
        std::cout << " done, " << tic_toc.RealTime() << "s." << std::endl;
    }

    // const real_t error = sol->ComputeLpError(2, u0);
    // if (Mpi::Root())
    // {
    //     std::cout << "Solution error: " << error << std::endl;
    // }

    delete ode_solver;
    
    return 0;
}

std::function<void(const Vector&, Vector&)> Cylinder(real_t Ma, real_t Re, real_t gamma, real_t mu, real_t D)
{
    return [Ma, Re, gamma, mu, D] (const Vector&x, Vector&y)
    {
        real_t rho = 1.0;
        real_t u = Re * mu / (D * rho);
        real_t p = rho * u * u / (gamma * Ma * Ma);
        y(0) = rho;
        y(1) = rho * u;
        y(2) = 0.0;
        y(3) = p / (gamma - 1.0) + 0.5 * rho * u * u;
    };
}

std::function<void(const Vector&, Vector&)> Sod1D(real_t gamma)
{
    return [gamma](const Vector &x, Vector &y)
    {
        real_t density, velocity_x, pressure, energy;
        MFEM_ASSERT(x.Size() == 1, "");

        density = 1.0;
        velocity_x = 0.0;
        pressure = 1000.0;

        energy = pressure / (1.4 - 1.0) + density * 0.5 * (velocity_x * velocity_x);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = energy;
    };
}

std::function<void(const Vector&, Vector&)> quiescent_if(int dim)
{
    return [dim] (const Vector &x, Vector &y)
    {
        y(0) = 1.0;
        y(1) = 0.0;
        if (dim > 1)
        {
            y(2) = 0.0;
            if (dim > 2)
            {
                y(3) = 0.0;
            }
        }
        y(dim + 1) = 1.0 ;
    };
}

std::function<void(const Vector&, Vector&)> LidDrivenCavityIC(real_t Ma, real_t gamma)
{
    return [Ma, gamma] (const Vector &x, Vector &y)
    {
        real_t p = 1.0 / (Ma * Ma * gamma);
        y(0) = 1.0;
        y(1) = 0.0;
        y(2) = 0.0;
        y(3) = p / (gamma - 1.0);
    };
}

std::function<real_t(const Vector&)> LidDrivenCavityAdiaBC()
{
    return [] (const Vector &x)
    {
        return 0.0;
    };
}

std::function<void(const Vector&, Vector&)> LidDrivenCavityVelBC()
{
    return [] (const Vector &x, Vector &vel)
    {
        vel(0) = 0.0;
        vel(1) = 0.0;
    };
}

std::function<void(const Vector&, Vector&)> LidDrivenCavityVelLidBC()
{
    return [] (const Vector &x, Vector &vel)
    {

        vel(0) = 1.0;
        vel(1) = 0.0;
    };
}

std::function<void(const Vector&, Vector&)> GetMovingVortexInit(
   const real_t radius, const real_t Minf, const real_t beta,
   const real_t gas_constant, const real_t specific_heat_ratio)
{
    return [specific_heat_ratio,
            gas_constant, Minf, radius, beta](const Vector &x, Vector &y)
   {
        MFEM_ASSERT(x.Size() == 2, "");

        const real_t xc = 0.0, yc = 0.0;

        // Nice units
        const real_t vel_inf = 1.;
        const real_t den_inf = 1.;

        // Derive remainder of background state from this and Minf
        const real_t pres_inf = (den_inf / specific_heat_ratio) *
                                (vel_inf / Minf) * (vel_inf / Minf);
        const real_t temp_inf = pres_inf / (den_inf * gas_constant);

        real_t r2rad = 0.0;
        r2rad += (x(0) - xc) * (x(0) - xc);
        r2rad += (x(1) - yc) * (x(1) - yc);
        r2rad /= (radius * radius);

        const real_t shrinv1 = 1.0 / (specific_heat_ratio - 1.);

        const real_t velX =
            vel_inf * (1 - beta * (x(1) - yc) / radius * std::exp(-0.5 * r2rad));
        const real_t velY =
            vel_inf * beta * (x(0) - xc) / radius * std::exp(-0.5 * r2rad);
        const real_t vel2 = velX * velX + velY * velY;

        const real_t specific_heat =
            gas_constant * specific_heat_ratio * shrinv1;
        const real_t temp = temp_inf - 0.5 * (vel_inf * beta) *
                            (vel_inf * beta) / specific_heat *
                            std::exp(-r2rad);

        const real_t den = den_inf * std::pow(temp / temp_inf, shrinv1);
        const real_t pres = den * gas_constant * temp;
        const real_t energy = shrinv1 * pres / den + 0.5 * vel2;

        y(0) = den;
        y(1) = den * velX;
        y(2) = den * velY;
        y(3) = den * energy;
   };
}

std::function<void(const Vector&, real_t, Vector&)> dmr_top_bf(const real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, real_t t, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        if (x(0) < 1.0 / 6.0 + (x(1) + 20.0 * t) / std::sqrt(3))
        {
            y(0) = 8.0;
            y(1) = 8.0 * 7.144709581221619;
            y(2) = -8.0 * 4.125;
            y(3) = 116.5 * gammaM1Inverse + 0.5 * y(0) * (7.144709581221619 * 7.144709581221619 + 4.125 * 4.125);
        }
        else
        {
            y(0) = 1.4;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 1.0 * gammaM1Inverse;
        }
    };
}

std::function<void(const Vector&, Vector&)> dmr_if(const real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        if (x(0) < 1.0 / 6.0 + x(1) / std::sqrt(3))
        {
            y(0) = 8.0;
            y(1) = 8.0 * 7.144709581221619;
            y(2) = -8.0 * 4.125;
            y(3) = 116.5 * gammaM1Inverse  + 0.5 * y(0) * (7.144709581221619 * 7.144709581221619 + 4.125 * 4.125);
        }
        else
        {
            y(0) = 1.4;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 1.0 * gammaM1Inverse;
        }
    };    
}

std::function<void(const Vector&, Vector&)> triple_point_if(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        if (x(0) <= 1.0 && x(1) <= 3.0)
        {
            y(0) = 1.0;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 1.0 * gammaM1Inverse;
        }
        else if (x(0) > 1.0 && x(1) <= 1.5)
        {
            y(0) = 1.0;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 0.1 * gammaM1Inverse;
        }
        else
        {
            y(0) = 0.125;
            y(1) = 0.0;
            y(2) = 0.0;
            y(3) = 0.1 * gammaM1Inverse;
        }
    };  
}

std::function<void(const Vector&, Vector&)> KelvinHelmholtzInstability_if(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        real_t density, velocity_x, velocity_y, pressure, energy, B;

        B = std::tanh(15.0 * x(1) + 7.5) - std::tanh(15.0 * x(1) - 7.5);
        density = 0.5 + 0.75 * B;
        velocity_x = 0.5 * (B - 1.0);
        velocity_y = 0.1 * std::sin(2.0 * M_PI * x(0));
        pressure = 1.0;
        energy = pressure / (1.4 - 1.0) + density * 0.5 * (velocity_x * velocity_x + velocity_y * velocity_y);

        y(0) = density;
        y(1) = density * velocity_x;
        y(2) = density * velocity_y;
        y(3) = energy;

        // if (x(1) < 0.5 && x(1) > -0.5)
        // {
        //     y(0) = 2.0;
        //     y(1) = -1.0;
        //     y(2) = 2.0 * 0.01 * std::sin(M_PI * x(0));
        //     y(3) = 2.5 * gammaM1Inverse  + 0.5 * (y(1) * y(1) + y(2) * y(2)) / y(0);
        // }
        // else
        // {
        //     y(0) = 1.0;
        //     y(1) = 0.5;
        //     y(2) = 0.01 * std::sin(M_PI * x(0));;
        //     y(3) = 2.5 * gammaM1Inverse + 0.5 * (y(1) * y(1) + y(2) * y(2)) / y(0);
        // }
    }; 
}

std::function<void(const Vector&, Vector&)> airfoil_if(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        y(0) = 1.0;
        y(1) = 1.5; // * std::cos(0.0349066 * 5.0);
        y(2) = 0.0; //1.5 * std::sin(0.0349066 * 5.0);
        y(3) = 1.0 * gammaM1Inverse + 0.5 * y(0) * (std::pow(y(1) / y(0), 2) + std::pow(y(2) / y(0), 2));
    };
}

std::function<void(const Vector&, Vector&)> ramp_if(real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        y(0) = 1.0;
        y(1) = 3.0;
        y(2) = 0.0;
        y(3) = 1.0 * gammaM1Inverse + 0.5 * y(0) * (std::pow(y(1) / y(0), 2) + std::pow(y(2) / y(0), 2));
    };
}