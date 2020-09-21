#
# Toolchain file for CMake to build  on Fermi
#
# cmake -C initialCache.cmake
#

set( CMAKE_C_COMPILER CC CACHE STRING "C compiler  wrapper" FORCE )
set( CMAKE_CXX_COMPILER CXX CACHE STRING "C++ compiler  wrapper" FORCE )
set( CMAKE_Fortran_COMPILER FC CACHE STRING "Fortran compiler wrapper" FORCE )

#set( BLA_VENDOR Intel )
set( msg "Flags used by the compiler during Release build")
set( CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -axSSE4.2" CACHE STRING ${msg})
set( CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -axSSE4.2" CACHE STRING ${msg} )
set( CMAKE_Fortran_FLAGS_RELEASE "-O3 -DNDEBUG -axSSE4.2" CACHE STRING ${msg} )

option( HIDE_Remarks "Sets the compiler option to disable compiler remarks" ON )
if (HIDE_Remarks)
  add_definitions("-diag-disable vec -diag-disable 279")
endif (HIDE_Remarks)

set( MPI_COMPILER CXX CACHE STRING "MPI Compiler wrapper for path detection" FORCE)
