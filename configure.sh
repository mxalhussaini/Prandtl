#! /bin/sh

# Accept an optional command line argument for the config file path.
# If none is provided, default to "TestCases/NavierStokes/2D/LidDrivenCavity/config.json"
CONFIG_FILE=${1:-"TestCases/NavierStokes/2D/LidDrivenCavity/config.json"}

# Clean previous build
rm -rf out/build

# (Re-)run CMake configuration
cmake -S . -B out/build -DCONFIG_FILE="${CONFIG_FILE}"
