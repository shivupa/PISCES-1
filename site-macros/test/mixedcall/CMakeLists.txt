# ==============================================================================
#   PROJECT TEST_MIXEDCALL 
# ==============================================================================
project( test_mixedcall C CXX Fortran )
cmake_minimum_required(VERSION 2.6)


# ==============================================================================
#   FIND PACKAGES
# ==============================================================================
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../../Modules)
find_package( FORTRANRTL REQUIRED )


# ==============================================================================
#   TARGETS and LINKS
# ==============================================================================


add_library( mixedcallf sub.f90 )
add_executable( test_mixedcall main.cxx fun.cxx  )

target_link_libraries( test_mixedcall 
   mixedcallf ${FORTRANRTL_LIBRARIES} )

set_property(
   TARGET test_mixedcall
   APPEND
   PROPERTY LINK_FLAGS ${FORTRANRTL_LINK_FLAGS} )










