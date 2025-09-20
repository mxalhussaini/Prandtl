#! /bin/sh

# exit if a command fails
set -e

BUILD_DIR="out/build"

# first argument is the number of processes, default is 4
N_PROCS=${1:-4}
# second argument is the config file path
CONFIG_FILE=$2
cd "$BUILD_DIR"

if [ -z "$CONFIG_FILE" ]; then
    mpiexec -n ${N_PROCS} ./Prandtl
    #ibrun Prandtl
else
    mpiexec -n ${N_PROCS} ./Prandtl -c "../../${CONFIG_FILE}"
    #ibrun Prandtl "${CONFIG_FILE}"
fi