cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_COMPILER ${MPICXX})
set(MPI_CXX_INCLUDE_PATH "")
set(MPI_CXX_LIBRARIES "")


project(Prandtl)

add_executable(${PROJECT_NAME} main.cpp)

target_sources(${PROJECT_NAME}
    PRIVATE src/Operators/DGOperator.cpp 
    PRIVATE src/RiemannSolvers/RoeFlux.cpp
    PRIVATE src/RiemannSolvers/LaxFriedrichsFlux.cpp
    PRIVATE src/RiemannSolvers/HLLFlux.cpp
    PRIVATE src/BasicOperations/BasicOperations.cpp
    PRIVATE src/Filters/EntropyFilter.cpp
    PRIVATE src/BoundaryConditions/BdrFaceIntegrator.cpp
    PRIVATE src/BoundaryConditions/FixedStateBdrFaceIntegrator.cpp
    PRIVATE src/BoundaryConditions/OutletBdrFaceIntegrator.cpp
    PRIVATE src/BoundaryConditions/TimeDependentBdrFaceIntegrator.cpp
    PRIVATE src/BoundaryConditions/WallBdrFaceIntegrator.cpp
    PRIVATE src/ModalBasis/ModalBasis.cpp
    )

# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/mfem-4.7)
# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/hypre-2.31.0/src)
# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/metis-5.1.0)



target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/mfem
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/hypre/src/hypre/include
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/metis-5.1.0/include
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Operators
    PUBLIC ${CMAKE_SOURCE_DIR}/src/RiemannSolvers
    PUBLIC ${CMAKE_SOURCE_DIR}/src/BasicOperations
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Filters
    PUBLIC ${CMAKE_SOURCE_DIR}/src/BoundaryConditions
    PUBLIC ${CMAKE_SOURCE_DIR}/src/ModalBasis
    )

target_link_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/mfem
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/hypre/src/lib
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/metis-5.1.0/lib
    )

target_link_libraries(${PROJECT_NAME} mfem HYPRE  metis)


