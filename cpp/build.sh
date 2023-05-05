#!/usr/bin/env bash
# Script to build the KineticGas module and run tests. Can be used with options
# --fullclean : Wipe the release directory before rebuilding
# --cleancache : Delete the cache in the release directory before rebuilding
# --Debug : Build a debug version
# NOTE: The binary in `pykingas/` will not be overwritten unless tests succeed.
set -e

export CC=/usr/bin/gcc
export CXX=/usr/bin/g++

echo "Building KineticGas Release"

[ ! -d "cpp/release" ] && mkdir cpp/release
cd cpp/release
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ../..
echo "Copying binary to ${PWD}/pykingas/KineticGas_r.so"
cp cpp/release/KineticGas_r.* pykingas/KineticGas_r.so