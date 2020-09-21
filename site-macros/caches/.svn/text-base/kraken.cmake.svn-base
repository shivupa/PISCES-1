#
# Toolchain file for CMake to build  on Fermi
#
# cmake -C initialCache.cmake
#

set( CMAKE_C_COMPILER cc CACHE STRING "C compiler  wrapper" FORCE )
set( CMAKE_CXX_COMPILER CC CACHE STRING "C++ compiler  wrapper" FORCE )
set( CMAKE_Fortran_COMPILER ftn CACHE STRING "Fortran compiler wrapper" FORCE )

#set( BLA_VENDOR Intel )
set( msg "Flags used by the compiler during Release build")
set( CMAKE_CXX_FLAGS_RELEASE "-target=compute_node  -DNDEBUG -O3 -ip" CACHE STRING ${msg})
set( CMAKE_C_FLAGS_RELEASE "-target=compute_node  -DNDEBUG -O3 -ip" CACHE STRING ${msg} )
set( CMAKE_Fortran_FLAGS_RELEASE "-target=compute_node  -DNDEBUG -O3 -ip" CACHE STRING ${msg} )

#set( msg "Cross compilation flag used for all builds on Kraken")
#set( CMAKE_CXX_FLAGS "-target=compute_node" CACHE STRING ${msg})
#set( CMAKE_C_FLAGS "-target=compute_node" CACHE STRING ${msg} )
#set( CMAKE_Fortran_FLAGS "-target=compute_node" CACHE STRING ${msg} )
#set( CMAKE_EXE_LINKER_FLAGS "-target=compute_node" CACHE STRING ${msg} )



#option( HIDE_Remarks "Sets the compiler option to disable compiler remarks" ON )
#if (HIDE_Remarks)
#  add_definitions("-diag-disable vec -diag-disable 279")
#endif (HIDE_Remarks)

#set( MPI_COMPILER CXX CACHE STRING "MPI Compiler wrapper for path detection" FORCE)
