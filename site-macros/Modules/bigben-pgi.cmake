#
# Toolchain file for CMake build on Bigben for
# CATAMOUNT nodes
#
# Use as `cmake -DCMAKE_TOOLCHAIN_FILE=path/to/this/file'
#

# 
# Notes for Bigben:
# . when building with craypat, do NOT load craypat before building all object. 
#   Load it for link time only. Not sure, but Fortran objects are causeing some 
#   sort of a problem. 
#


set( CMAKE_SYSTEM_NAME Catamount )
set( CMAKE_C_COMPILER cc -target=catamount )
set( CMAKE_CXX_COMPILER CC -target=catamount )
set( CMAKE_Fortran_COMPILER  ftn )
# can not append -target to ftn. CMakeFortranInformation.cmake gives error.
# adding to FLAGS

# set the search path for the environment coming with the compiler
# and a directory where you can install your own compiled software
set(CMAKE_FIND_ROOT_PATH
    /opt/xt-mpt/1.5.60/mpich2-64/P2
    /opt/acml/3.0/pgi64
    /opt/xt-libsci/10.0.0/pgi/cnos64
    /opt/xt-mpt/1.5.60/sma/P2
    /opt/xt-lustre-ss/1.5.60/
    /opt/xt-libc/1.5.60/amd64
    /opt/pgi/7.2.2/linux86-64/7.2.2
    /opt/xt-pe/default/lib/cnos64
    )

#
# CXX
#
set( CMAKE_CXX_FLAGS_INIT "--diag_suppress 236,9,177,185,368,381,1396 --display_error_number" )
set( CMAKE_CXX_FLAGS_DEBUG_INIT "-g -O0")
set( CMAKE_CXX_FLAGS_RELEASE_INIT "-fastsse -s")
set( CMAKE_CXX_FLAGS_MINSIZEREL_INIT "-O2 -s")
set( CMAKE_CXX_FLAGS_RELWITHDEBINFO_INIT "-O2 -gopt")


#
# C
#
#set( CMAKE_C_FLAGS_INIT )
set( CMAKE_C_FLAGS_DEBUG_INIT "-g -O0")
set( CMAKE_C_FLAGS_RELEASE_INIT "-fastsse -s")
set( CMAKE_C_FLAGS_MINSIZEREL_INIT "-O2 -s")
set( CMAKE_C_FLAGS_RELWITHDEBINFO_INIT "-O2 -gopt")


#
# Fortran
#
set(CMAKE_Fortran_MODDIR_FLAG "-module ")
set(CMAKE_Fortran_FLAGS_INIT "-target=catamount -Kieee -Mextend")
set(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g -O0 -Mbounds")
set(CMAKE_Fortran_FLAGS_MINSIZEREL_INIT "-O2 -s")
set(CMAKE_Fortran_FLAGS_RELEASE_INIT "-fast -O3 -Mipa=fast")
set(CMAKE_Fortran_FLAGS_RELWITHDEBINFO_INIT "-O2 -gopt")
set(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")


#
# Libs
#
#set(MPI_INCLUDE_PATH "/opt/xt-mpt/default/mpich2-64/P2")
set(MPI_LIBRARY "mpich")
set(BLA_STATIC 1)
set(BLA_VENDOR "ACML")
add_definitions( -DBOOST_UBLAS_UNSUPPORTED_COMPILER=0 )

#
# Special loads
#
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREAD OFF)
set(BOOST_INCLUDEDIR /usr/users/7/yilmaz/include/boost-1_34_1 )
set(BOOST_LIBRARYDIR /usr/users/7/yilmaz/lib )

#
if ( MPI_FOUND )
  add_definitions( -DMPICH_IGNORE_CXX_SEEK )
endif ( MPI_FOUND )
