#include "Prandtl.hpp"

using namespace mfem;

std::function<void(const Vector&, real_t, Vector&)> dmr_top_bf(real_t gammaM1Inverse);
std::function<void(const Vector&, Vector&)> dmr_if(real_t gammaM1Inverse);

const real_t Prandtl::gamma = 1.4;
const real_t Prandtl::gammaM1 = Prandtl::gamma - 1.0;
const real_t Prandtl::gammaP1 = Prandtl::gamma + 1.0;
const real_t Prandtl::gammaM1Inverse = 1.0 / Prandtl::gammaM1;
const real_t Prandtl::gammaP1Inverse = 1.0 / Prandtl::gammaP1;

int main (int argc, char* argv[])
{

    // std::string filename = "../../SupercriticalAirfoilQuad.msh";
    // Mesh mesh(filename);

    // if (Mpi::Root())
    // {
    //     mesh.bdr_attributes.Print(std::cout);
    // }

    Mpi::Init(argc, argv);
    const int numProcs = Mpi::WorldSize();
    const int myRank = Mpi::WorldRank();
    Hypre::Init();

    int order = 3;
    real_t t_final = 2.0;
    real_t dt = -0.01;
    real_t cfl = 0.3;

    int precision = 15;
    std::cout.precision(precision);

    Mesh *mesh = new Mesh();
    *mesh = Mesh::MakeCartesian2D(384, 96, Element::QUADRILATERAL, false, 4.0, 1.0, true);

    int dim = mesh->Dimension();
    int num_equations = dim + 2;

    if (dim > 1) mesh->EnsureNCMesh();
    std::shared_ptr<ParMesh> pmesh = std::make_shared<ParMesh>(MPI_COMM_WORLD, *mesh);
    mesh->Clear();

    Array<int> bdr_elem_vtx;
    for (int k = 0; k < pmesh->GetNBE(); k++)
    {
        pmesh->GetBdrElementVertices(k, bdr_elem_vtx);
        double* coord = pmesh->GetVertex(bdr_elem_vtx[1]);

        if (coord[0] > 1.0 / 6.0 && coord[1] < 1e-15)
        {
            pmesh->SetBdrAttribute(k, 5);
        }
    }
    pmesh->SetAttributes();

    Array<int> top(pmesh->bdr_attributes.Max());
    Array<int> bottom1(pmesh->bdr_attributes.Max());
    Array<int> bottom2(pmesh->bdr_attributes.Max());
    Array<int> left(pmesh->bdr_attributes.Max());
    Array<int> right(pmesh->bdr_attributes.Max());

    top = 0; top[2] = 1;
    bottom1 = 0; bottom1[0] = 1;
    bottom2 = 0; bottom2[4] = 1;
    left = 0; left[3] = 1;
    right = 0; right[1] = 1;

    Vector fixed_state(num_equations);
    fixed_state(0) = 8.0;
    fixed_state(1) = 8.0 * 7.144709581221619;
    fixed_state(2) = -8.0 * 4.125;
    fixed_state(3) = 116.5 / 0.4 + 0.5 * fixed_state(0) * (7.144709581221619 * 7.144709581221619 + 4.125 * 4.125);



    DG_FECollection fec(order, dim, BasisType::GaussLobatto);
    DG_FECollection fec0(0, dim);
    std::shared_ptr<ParFiniteElementSpace> vfes = std::make_shared<ParFiniteElementSpace>(pmesh.get(), &fec, num_equations, Ordering::byNODES);
    std::shared_ptr<ParFiniteElementSpace> fes0 = std::make_shared<ParFiniteElementSpace>(pmesh.get(), &fec0);
    
    VectorFunctionCoefficient u0(num_equations, dmr_if(Prandtl::gammaM1Inverse));
    VectorFunctionCoefficient u_top(num_equations, dmr_top_bf(Prandtl::gammaM1Inverse));

    std::shared_ptr<ParGridFunction> sol = std::make_shared<ParGridFunction>(vfes.get());
    sol->ProjectCoefficient(u0);

    EulerFlux flux(dim, Prandtl::gamma);
    RusanovFlux numericalFlux(flux);

    pmesh->ExchangeFaceNbrData();
    const Table &el2el = pmesh->ElementToElementTable();
    const Table &el2bdrel = Prandtl::ElementIndextoBdrElementIndex(*pmesh);

    // std::cout << el2el.Size() << "\n";
    // std::cout << el2bdrel.Size() << "\n";
    // std::cout << vfes->GetNE() << "\n";
    // std::cout << pmesh->Nonconforming() << "\n";
    // std::cout << "-------" << std::endl;

    // Array<int> row;

    if (Mpi::Root())
    {
        // for (int i = 0; i < el2bdrel.Size(); i++)
        // {
            // el2bdrel.GetRow(i, row);
            // row.Prepend(i);
            // row.Print(std::cout);
            // for (int j = 0; j < row.Size(); j++)
            // {
            //     if (row[j] > vfes->GetNE())
            //     {
            //         std::cout << row[j] << std::endl;
            //         row.Print(std::cout);
            //         std::cout << row.Size() << std::endl;
            //         std::cout << "---------" << std::endl;
            //     }
            // }
        // }
    }
    

    const int int_order = order * 2 + 3;
    IntegrationRules GLIntRules(0, Quadrature1D::GaussLobatto);
    const IntegrationRule *vol_ir = &GLIntRules.Get(Geometry::SQUARE, int_order);
    const IntegrationRule *face_ir = &GLIntRules.Get(Geometry::SEGMENT, int_order);
    const IntegrationRule *bdr_face_ir = &GLIntRules.Get(Geometry::SEGMENT, int_order);

    Geometry::Type gtype = vfes->GetFE(0)->GetGeomType();

    Prandtl::DGOperator euler(vfes, fes0, pmesh, sol,
        std::make_unique<Prandtl::DGFormIntegrator>(numericalFlux, vol_ir, face_ir),
        std::make_unique<Prandtl::EntropyFilter>(vfes, fes0, sol,
        pmesh, std::make_unique<ParGridFunction>(fes0.get()),
        std::make_unique<Prandtl::ModalBasis>(fec, gtype, order, dim),
        el2el, el2bdrel, vol_ir, face_ir, bdr_face_ir), true); 

    real_t t = 0.0;
    euler.SetTime(t);

    // Array<Prandtl::BdrFaceIntegrator*> bfnfi;
    // Array<Array<int>> bfnfi_marker;

    // bfnfi.Append(new Prandtl::WallBdrFaceIntegrator(numericalFlux, bdr_face_ir));
    // bfnfi_marker.Append(bottom2);

    // bfnfi.Append(new Prandtl::OutletBdrFaceIntegrator(numericalFlux, bdr_face_ir));
    // bfnfi_marker.Append(right);

    // bfnfi.Append(new Prandtl::FixedStateBdrFaceIntegrator(numericalFlux, bdr_face_ir, fixed_state));
    // bfnfi_marker.Append(left);

    // bfnfi.Append(new Prandtl::FixedStateBdrFaceIntegrator(numericalFlux, bdr_face_ir, fixed_state));
    // bfnfi_marker.Append(bottom1);

    // bfnfi.Append(new Prandtl::TimeDependentBdrFaceIntegrator(numericalFlux, bdr_face_ir, u_top, euler.GetTimeRef()));
    // bfnfi_marker.Append(top);

    euler.AddBdrFaceIntegrator(new Prandtl::WallBdrFaceIntegrator(numericalFlux, bdr_face_ir), bottom2);
    euler.AddBdrFaceIntegrator(new Prandtl::OutletBdrFaceIntegrator(numericalFlux, bdr_face_ir), right);
    euler.AddBdrFaceIntegrator(new Prandtl::FixedStateBdrFaceIntegrator(numericalFlux, bdr_face_ir, fixed_state), left);
    euler.AddBdrFaceIntegrator(new Prandtl::FixedStateBdrFaceIntegrator(numericalFlux, bdr_face_ir, fixed_state), bottom1);
    euler.AddBdrFaceIntegrator(new Prandtl::TimeDependentBdrFaceIntegrator(numericalFlux, bdr_face_ir, u_top, euler.GetTimeRef()), top);

    Prandtl::DGODESolver *ode_solver = new Prandtl::RK3SSPExplicitSolver();

    tic_toc.Clear();
    tic_toc.Start();

    ode_solver->Init(euler);

    dt = 0.00001;


    ode_solver->Step(*sol, t, dt);

    delete ode_solver;
    
    return 0;
}

std::function<void(const Vector&, real_t, Vector&)> dmr_top_bf(const real_t gammaM1Inverse)
{
    return [gammaM1Inverse] (const Vector &x, real_t t, Vector &y)
    {
        MFEM_ASSERT(x.Size() == 2, "DMR is a 2D problem");

        if (x(0) <= 1.0 / 6.0 + (x(1) + 20.0 * t) / std::sqrt(3))
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

        if (x(0) <= 1.0 / 6.0 + x(1) / std::sqrt(3))
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
