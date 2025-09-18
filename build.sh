#! /bin/sh

# exit if a command fails
set -e

cd out/build ; cmake --build . -- "$@"