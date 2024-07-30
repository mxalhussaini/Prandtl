cmake_minimum_required(VERSION 3.15)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_COMPILER ${MPICXX})
set(MPI_CXX_INCLUDE_PATH "")
set(MPI_CXX_LIBRARIES "")

project(Prandtl)

add_executable(${PROJECT_NAME} main.cpp)

# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/mfem-4.7)
# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/hypre-2.31.0/src)
# add_subdirectory(${CMAKE_SOURCE_DIR}/libs/metis-5.1.0)


target_include_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/mfem-4.7
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/hypre-2.31.0/src/hypre/include
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/metis-5.1.0/include
    )

target_link_directories(${PROJECT_NAME}
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/mfem-4.7
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/hypre-2.31.0/src/lib
    PUBLIC ${CMAKE_SOURCE_DIR}/libs/metis-5.1.0/lib
    )

target_link_libraries(${PROJECT_NAME} mfem HYPRE  metis)


