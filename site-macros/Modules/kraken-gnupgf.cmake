#
# Toolchain file for CMake build on Bigben for
# CATAMOUNT nodes
#
# Use as `cmake -DCMAKE_TOOLCHAIN_FILE=path/to/this/file'
#

#
# Kraken is not really Catamount.
# But they share the feature not to allow shared libs, 
# and Catamount sets this fine. The other way to do this is
# define a Platform module, which is more headache. 
#
SET(CMAKE_SYSTEM_NAME Catamount)

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
  /opt/fftw/3.1.1/cnos/include
  /opt/mpt/3.1.0/xt/mpich2-gnu
  /opt/acml/4.1.0/gnu64
  /opt/xt-libsci/10.3.1/gnu/snos64
  /opt/xt-libsci/10.3.1/gnu/snos64/include/superlu
  /opt/mpt/3.1.0/xt/sma
  /opt/mpt/3.1.0/xt/pmi
  /opt/xt-pe/2.1.56HD
  /opt/fftw/default/cnos
  /opt/mpt/default/xt/mpich2-pgi
  /opt/acml/default/pgi64
  /opt/xt-libsci/default/pgi/snos64
  /opt/pgi/7.2.5/linux86-64/7.2-5
)


# compilers
set( CMAKE_C_COMPILER cc )
set( CMAKE_CXX_COMPILER CC )
set( CMAKE_Fortran_COMPILER /opt/pgi/7.2.5/linux86-64/7.2/bin/pgf90 )
 # ftn -target=linux doesn't work, need to set this below

# libraries
set( Boost_ADDITIONAL_VERSIONS 1.39 1.39.0 )
set( BLAS_FOUND TRUE )
set( LAPACK_FOUND TRUE )
add_definitions( "-DMPICH_IGNORE_CXX_SEEK" "-DBOOST_UBLAS_UNSUPPORTED_COMPILER=0" )

#
# CXX
#
#set( CMAKE_CXX_FLAGS_INIT "-target=linux" )
# "--diag_suppress 236,9,177,185,368,381,1396 --display_error_number" )
#set( CMAKE_CXX_FLAGS_DEBUG_INIT "-g -O0")
#set( CMAKE_CXX_FLAGS_RELEASE_INIT "-fastsse -s")
#set( CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O2 -s")
#set( CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -gopt")


#
# C
#
#set( CMAKE_C_FLAGS_INIT "-target=linux")
# set( CMAKE_C_FLAGS_DEBUG_INIT "-g -O0")
# set( CMAKE_C_FLAGS_RELEASE_INIT "-fastsse -s")
# set( CMAKE_C_FLAGS_MINSIZEREL_INIT "-O2 -s")
# set( CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -gopt")


#
# Fortran
#
# set(CMAKE_Fortran_MODDIR_FLAG "-module ")
#set(CMAKE_Fortran_FLAGS_INIT "-target=catamount")
# -Kieee -Mextend)
# set(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g -O0 -Mbounds")
# set(CMAKE_Fortran_FLAGS_MINSIZEREL_INIT "-O2 -s")
# set(CMAKE_Fortran_FLAGS_RELEASE_INIT "-fast -O3 -Mipa=fast")
# set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -gopt")
# set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")


#
# Libs
#
# set(BOOST_ROOT "$ENV{HOME}")
# set(MPI_INCLUDE_PATH "/opt/xt-mpt/default/mpich2-64/P2")
# set(MPI_LIBRARY "mpich")
# set(BLA_STATIC 1)
# set(BLA_VENDOR "ACML")
