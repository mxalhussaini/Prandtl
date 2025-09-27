#! /bin/sh

# exit if a command fails
set -e

IS_HPC=false
# --- Environment Detection ---
if command -v module &> /dev/null; then
    IS_HPC=true
fi

# --- Load HPC Environment ---
if [ "$IS_HPC" = true ]; then
    source ./env_hpc.sh
fi

# Accept an optional command line argument for the config file path.
# If none is provided, default to DEFAULT_CONFIG
DEFAULT_CONFIG="TestCases/Euler/1D/SodShockTube/config.json"

BUILD_DIR="out/build"

CONFIG_FILE=${1:-$DEFAULT_CONFIG}

# Clean previous build
rm -rf "$BUILD_DIR"

# (Re-)run CMake configuration
cmake -S . -B "$BUILD_DIR" -DCONFIG_FILE="${CONFIG_FILE}" -DCMAKE_BUILD_TYPE="Release" -DCMAKE_CXX_FLAGS_RELEASE="-O2 -DNDEBUG"
