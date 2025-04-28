#! /bin/sh


cd out/build

if [ -z "$1" ]; then
    mpiexec -n 4 Prandtl
    #ibrun Prandtl
else
    mpiexec -n 4 Prandtl "$1"
    #ibrun Prandtl "$1"
fi