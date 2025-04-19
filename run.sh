#! /bin/sh


cd out/build

if [ -z "$1" ]; then
    mpiexec -n 5 Prandtl
else
    mpiexec -n 5 Prandtl "$1"
fi