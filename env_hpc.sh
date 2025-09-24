#!/bin/bash

# A script for setting the enviroment variables for the HPC
# Usage: source env_hpc.sh

set -e

module unload gcc openmpi

echo "--- Loading Intel Compiler and Intel MPI Environment ---"
module load intel/24.0 impi/21.11 cmake/3.31.5
export CC=mpicc
export CXX=mpicxx

