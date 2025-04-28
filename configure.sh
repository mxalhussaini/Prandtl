#! /bin/sh

# Accept an optional command line argument for the config file path.
# If none is provided, default to "TestCases/Euler/1D/SodShockTube/config.json"
CONFIG_FILE=${1:-"TestCases/Euler/1D/SodShockTube/config.json"}

# Clean previous build
rm -rf out/build

# (Re-)run CMake configuration
cmake -S . -B out/build -DCONFIG_FILE="${CONFIG_FILE}"
# cmake -S . -B out/build -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DCONFIG_FILE="${CONFIG_FILE}"
