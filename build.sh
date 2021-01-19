#!/bin/bash
#rm -rf build
mkdir -p build
cd build

CC=mpicc \
CXX=mpic++ \
FC=mpifort \
CFLAGS="-Wall -Wextra" \
CXXFLAGS="-Wall -Wextra" \
FFLAGS="-Wall -Wextra" \
cmake .. -DHAVE_ARPACK=1 -DHAVE_FFTW=1 -DPISCES_DOCS=1 -DCMAKE_BUILD_TYPE=Release
make
#make doc
make test
#open build/doc/html/index.html
