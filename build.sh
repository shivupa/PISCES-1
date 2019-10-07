#!/bin/bash

mkdir -p build
cd build

CC=mpicc \
CXX=mpic++ \
FC=mpifort \
CXXFLAGS="-Wall -Wextra -O3" \
cmake .. -DHAVE_ARPACK=1 -DHAVE_FFTW=1
make
