#
# Toolchain file for CMake to build  with intel compilers
#
# Use as `cmake -DCMAKE_TOOLCHAIN_FILE=path/to/this/file'
#

set( CMAKE_C_COMPILER icc )
set( CMAKE_CXX_COMPILER icpc )
set( CMAKE_Fortran_COMPILER ifort )

set( BLA_VENDOR Intel CACHE STRING "Using intel for BLAS libaries.")

option( HIDE_Remarks "Sets the compiler option to disable compiler remarks" ON )
if (HIDE_Remarks)
  add_definitions("-diag-disable vec -diag-disable 279")
endif (HIDE_Remarks)



set(h "Release configuration to generate profiling data. Use PROFUSE for later compilation with profiler data.")
set(f "-DNDEBUG -O3 -prof-gen")
set(CMAKE_CXX_FLAGS_PROFGEN      ${f} CACHE STRING  ${h})
set(CMAKE_C_FLAGS_PROFGEN        ${f} CACHE STRING  ${h})
set(CMAKE_Fortran_FLAGS_PROFGEN  ${f} CACHE STRING  ${h})
mark_as_advanced( CMAKE_CXX_FLAGS_PROFGEN CMAKE_C_FLAGS_PROFGEN CMAKE_Fortran_FLAGS_PROFGEN )

set(h "Release configuration that uses generated profiling data. Use PROFGEN to generate the data.")
set(f "-DNDEBUG -O3 -prof-use -prof-dir=${PROJECT_BINARY_DIR}")
set(CMAKE_CXX_FLAGS_PROFUSE      ${f} CACHE STRING  ${h})
set(CMAKE_C_FLAGS_PROFUSE        ${f} CACHE STRING  ${h})
set(CMAKE_Fortran_FLAGS_PROFUSE  ${f} CACHE STRING  ${h})
mark_as_advanced( CMAKE_CXX_FLAGS_PROFUSE CMAKE_C_FLAGS_PROFUSE CMAKE_Fortran_FLAGS_PROFUSE )
