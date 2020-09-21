#
# Determine which compiler is runnning
# Sets,
#
#   UsingIntelCompiler
#   Using...Compiler
#

if(CMAKE_CXX_COMPILER MATCHES "icpc$")
  set(UsingIntelCompiler 1)
  message( STATUS "Using Intel compiler" )
elseif( CMAKE_CXX_COMPILER_ID MATCHES "PGI" )
  set(UsingPGICompiler 1)
  message( STATUS "Using PGI compiler" )
endif(CMAKE_CXX_COMPILER MATCHES "icpc$")
