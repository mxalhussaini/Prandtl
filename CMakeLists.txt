cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set(CMAKE_CXX_COMPILER ${MPICXX})
# set(MPI_CXX_INCLUDE_PATH "")
# set(MPI_CXX_LIBRARIES "")

# add_definitions(-DHYPERBOLIC)
add_definitions(-DUSE_FILTER)
add_definitions(-DPARABOLIC)
add_definitions(-DSUBCELL_FV_BLENDING)

# # Add AddressSanitizer flags
# set(SANITIZER_FLAGS "-fsanitize=address -g -O1")

# # Apply sanitizer flags to compile and link phases
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SANITIZER_FLAGS}")
# set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} ${SANITIZER_FLAGS}")

project(Prandtl)

find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
set(MPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES})

message(STATUS "MPI CXX COMPILER: ${MPI_CXX_COMPILER}")
message(STATUS "MPI INCLUDE PATH: ${MPI_INCLUDE_PATH}")
message(STATUS "MPI LIBRARIES: ${MPI_CXX_LIBRARIES}")

add_executable(${PROJECT_NAME} main.cpp)

target_sources(${PROJECT_NAME}
    PRIVATE src/Operators/DGSEMOperator.cpp 
    PRIVATE src/Integrators/DGSEMIntegrator/DGSEMIntegrator.cpp
    PRIVATE src/Integrators/DGSEMIntegrator/Lifting/LiftingScheme.cpp
    PRIVATE src/Integrators/DGSEMIntegrator/Lifting/LiftingBR1/LiftingBR1.cpp
    PRIVATE src/RiemannSolvers/RoeFlux.cpp
    PRIVATE src/RiemannSolvers/LaxFriedrichsFlux.cpp
    PRIVATE src/RiemannSolvers/HLLFlux.cpp
    PRIVATE src/RiemannSolvers/ChandrashekarFlux.cpp
    PRIVATE src/BasicOperations/BasicOperations.cpp
    PRIVATE src/Filters/EntropyFilter.cpp
    PRIVATE src/Integrators/BoundaryConditions/BdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/NoSlipAdiabWallBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/NoSlipIsothWallBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SpecifiedStateBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SlipWallBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/RiemannInvariantBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SubsonicInflowPtTtAngBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SubsonicInflowRVBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SubsonicOutflowPBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SupersonicInflowBdrFaceIntegrator.cpp
    PRIVATE src/Integrators/BoundaryConditions/SupersonicOutflowBdrFaceIntegrator.cpp
    PRIVATE src/Fluxes/NavierStokesFlux.cpp
    PRIVATE src/ModalBasis/ModalBasis.cpp
    PRIVATE src/Limiters/MinmodLimiter.cpp
    PRIVATE src/Indicators/Indicator.cpp
    PRIVATE src/Indicators/PerssonPeraireIndicator.cpp
    PRIVATE src/Simulation/Simulation.cpp
    PRIVATE src/Operators/NonlinearForm/DGSEMNonlinearForm.cpp
    )

# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/mfem-4.7)
# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/hypre-2.31.0/src)
# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/metis-5.1.0)



target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/mfem
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/hypre/src/hypre/include
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/metis-5.1.0/include
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Operators/NonlinearForm
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Operators
    PUBLIC ${CMAKE_SOURCE_DIR}/src/RiemannSolvers
    PUBLIC ${CMAKE_SOURCE_DIR}/src/BasicOperations
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Filters
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Integrators/BoundaryConditions
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Integrators/DGSEMIntegrator/Lifting/LiftingBR1
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Integrators/DGSEMIntegrator/Lifting
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Integrators/DGSEMIntegrator
    PUBLIC ${CMAKE_SOURCE_DIR}/src/ModalBasis
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Physics
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Limiters
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Indicators
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Simulation
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Include
    PUBLIC ${CMAKE_SOURCE_DIR}/src/Fluxes
    )

target_link_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/mfem
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/hypre/src/lib
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/metis-5.1.0/lib
    )

target_link_libraries(${PROJECT_NAME} mfem HYPRE metis MPI::MPI_CXX)


