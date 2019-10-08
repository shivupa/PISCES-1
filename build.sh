#!/bin/bash

mkdir -p build
cd build

CC=mpicc \
CXX=mpic++ \
FC=mpifort \
CFLAGS="-Wall -Wextra" \
CXXFLAGS="-Wall -Wextra" \
FFLAGS="-Wall -Wextra" \
cmake .. -DHAVE_ARPACK=1 -DHAVE_FFTW=1
make
