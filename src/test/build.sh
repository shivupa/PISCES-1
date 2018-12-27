#################################################################################
##
## To build all versions of the pi program source this script:
## $ . build-pi.sh
##
## Comment/Uncomment lines to select between different versions
##   to build.
##
#################################################################################
##
## Set environment
##   Note: You must use the same environment on the compute
##         nodes to run the applications
##
module purge
module load intel
module load intel-mpi
module load gcc/5.4.0 
module load fftw
module load arpack

export ARPACK_HOME=/ihome/kjordan/thc9/fermi/Development/ARPACK
export ARPACK_ROOT=/ihome/kjordan/thc9/fermi/Development/ARPACK
export LIBRARY_PATH=$LIBRARY_PATH:/ihome/kjordan/thc9/fermi/Development/ARPACK/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/ihome/kjordan/thc9/fermi/Development/ARPACK/lib
export CMAKE_LIBRARY_PATH=$CMAKE_LIBRARY_PATH:/ihome/kjordan/thc9/fermi/Development/ARPACK/lib
export fft_path=/ihome/crc/install/intel-2017.1.132/intel-mpi-2017.1.132/fftw/3.3.5

##
#################################################################################
##
## Compile serial c++ version
#################################################################################
##
## Compile mpi c++ version
mpiicc fft_test.cpp -qopenmp -lifcore libarpack.a libblas.so.3.4.2 -o fft_test.x
#mpiicc fft_test.cpp vtx_FFT.cpp -qopenmp -lifcore libarpack.a libblas.so.3.4.2 libfftw3.so.3.5.5 -o fft_test.x
##
## Compile mpi fortran version
##
#################################################################################
