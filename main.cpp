#include <iostream>
#include "mfem.hpp"
#include "DGOperator.hpp"
#include "BasicOperations.hpp"

using namespace mfem;

int main (int argc, char* argv[])
{
    // mfem::Vector vec(3);
    // vec = 1.0;
    // vec.Print(std::cout);
    Mpi::Init(argc, argv);
    // const int numProcs = Mpi::WorldSize();
    // const int myRank = Mpi::WorldRank();
    Hypre::Init();
    // if (Mpi::Root()) {
    //   mfem_warning(
    //       "The number of processor is larger than the number of elements.\n"
    //       "Refine serial meshes until the number of elements is large enough");
    // }

    Mesh *mesh = new Mesh();
    *mesh = Mesh::MakeCartesian2D(384, 96, Element::QUADRILATERAL, false, 4.0, 1.0, true);

    const Table& el2bdr_el = Prandtl::ElementIndextoBdrElementIndex(*mesh);    

    Array<int> row;

    for (int i = 0; i < 500; i++)
    {
        el2bdr_el.GetRow(i, row);
        if (row)
        {
            Vector center1(2);
            mesh->GetElementCenter(i, center1);

            for (int j = 0; j < row.Size(); j++)
            {
                ElementTransformation* Tr = mesh->GetBdrElementTransformation(row[j]);
                Vector center2(2);
                int geom = mesh->GetElementBaseGeometry(8643);
                Tr->Transform(Geometries.GetCenter(geom), center2);
                if (Mpi::Root())
                {
                    center1.Print(std::cout);
                    center2.Print(std::cout);
                    std::cout << i << "---" << row[j] << "\n";
                    std::cout << "-------" << std::endl;
                }
            }
        }
    }

    // Vector nor(3);

    // nor(0) = 1.5;
    // nor(1) = 2.4;
    // nor(2) = 3.7;

    // Prandtl::Normalize(nor);
    // if (std::abs(nor(0)) > 1e-12 || std::abs(nor(1)) > 1e-12)
    // {
    //     vec(2) = 1.0;
    // }
    // else
    // {
    //     vec(0) = 1.0;
    // }
    // Vector tan1 = Prandtl::Cross(nor, vec);
    // Prandtl::Normalize(tan1);
    // Vector tan2 = Prandtl::Cross(nor, tan1);
    // Prandtl::Normalize(tan2);

    // if (Mpi::Root())
    // {
    //     nor.Print(std::cout);
    //     tan1.Print(std::cout);
    //     tan2.Print(std::cout);

    //     std::cout << nor.Norml2() << " " << tan1.Norml2() << " " << tan2.Norml2() << "\n";

    //     Vector xy = Prandtl::Cross(nor, tan1);
    //     xy.Print(std::cout);

    //     Vector yz = Prandtl::Cross(tan1, tan2);
    //     yz.Print(std::cout);

    //     Vector zx = Prandtl::Cross(tan2, nor);
    //     zx.Print(std::cout);
    // }

    // real_t rho = 1.4;
    // real_t p = 100;
    // real_t e = Prandtl::ComputeEntropy(rho, p, 1.4);

    // std::cout << e << "\n";
    
    return 0;
}
