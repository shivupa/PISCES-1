#!/bin/bash

mkdir -p build
cd build

module load intel/2017.3.196 \
 intel-mpi/2017.3.196 \
 fftw/3.3.7 \
 arpack/1 \
 cmake/3.13.3 \
 doxygen/1.8.13


CC=mpiicc \
CXX=mpiicpc \
FC=mpiifort \
CXXFLAGS="-Wall -O3" \
cmake .. -DHAVE_ARPACK=1 -DHAVE_FFTW=1 -DFFTW_OMP_LIBRARY=$FFTWROOT -DPISCES_DOCS=1
make
make doc
