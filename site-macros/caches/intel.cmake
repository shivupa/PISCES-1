#
# INTEL compilers cache file. Use as
#
# cmake -C intel.cmake
#

set( CMAKE_C_COMPILER mpiicc CACHE STRING "C compiler  wrapper" FORCE )
set( CMAKE_CXX_COMPILER mpiicc CACHE STRING "C++ compiler  wrapper" FORCE )
set( CMAKE_Fortran_COMPILER ifort CACHE STRING "Fortran compiler wrapper" FORCE )

#set( BLA_VENDOR Intel CACHE STRING "BLA Vendor" FORCE )
set( msg "Flags used by the compiler during Release build")
set( CMAKE_CXX_FLAGS_RELEASE "-restrict -qopenmp -xHOST -O3 -ip -mkl" CACHE STRING ${msg})
set( CMAKE_C_FLAGS_RELEASE "-restrict -qopenmp -xHOST -O3 -ip -mkl" CACHE STRING ${msg} )
set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG" CACHE STRING ${msg} )

option( HIDE_Remarks "Sets the compiler option to disable compiler remarks" ON )
if (HIDE_Remarks)
  add_definitions("-diag-disable vec -diag-disable 279")
endif (HIDE_Remarks)


# ==============================================================================
# experimental stuff
# ==============================================================================

set(h "Release configuration to generate profiling data. Use PROFUSE for later compilation with profiler data.")
set(f "-DNDEBUG -O3 -prof-gen")
set(CMAKE_CXX_FLAGS_PROFGEN      ${f} CACHE STRING  ${h} FORCE)
set(CMAKE_C_FLAGS_PROFGEN        ${f} CACHE STRING  ${h} FORCE)
set(CMAKE_Fortran_FLAGS_PROFGEN  ${f} CACHE STRING  ${h} FORCE)
mark_as_advanced( CMAKE_CXX_FLAGS_PROFGEN CMAKE_C_FLAGS_PROFGEN CMAKE_Fortran_FLAGS_PROFGEN )

set(h "Release configuration that uses generated profiling data. Use PROFGEN to generate the data.")
set(f "-DNDEBUG -O3 -prof-use -prof-dir=${PROJECT_BINARY_DIR}")
set(CMAKE_CXX_FLAGS_PROFUSE      ${f} CACHE STRING  ${h} FORCE)
set(CMAKE_C_FLAGS_PROFUSE        ${f} CACHE STRING  ${h} FORCE)
set(CMAKE_Fortran_FLAGS_PROFUSE  ${f} CACHE STRING  ${h} FORCE)
mark_as_advanced( CMAKE_CXX_FLAGS_PROFUSE CMAKE_C_FLAGS_PROFUSE CMAKE_Fortran_FLAGS_PROFUSE )
