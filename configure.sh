#! /bin/sh

# exit if a command fails
set -e

# Accept an optional command line argument for the config file path.
# If none is provided, default to "TestCases/Euler/1D/SodShockTube/config.json"
DEFAULT_CONFIG="TestCases/Euler/1D/SodShockTube/config.json"

BUILD_DIR="out/build"

CONFIG_FILE=${1:-$DEFAULT_CONFIG}

# Clean previous build
rm -rf "$BUILD_DIR"

# (Re-)run CMake configuration
cmake -S . -B "$BUILD_DIR" -DCONFIG_FILE="${CONFIG_FILE}" -DCMAKE_BUILD_TYPE="Release"
# cmake -S . -B "$BUILD_DIR" -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCONFIG_FILE="${CONFIG_FILE}"
