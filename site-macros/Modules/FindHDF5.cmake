#
# this module look for HDF5 (http://hdf.ncsa.uiuc.edu) support
# it will define the following values
#
# HDF5_INCLUDE_DIR  = where hdf5.h can be found
# HDF5_LIBRARY      = the library to link against
# HDF5_FOUND        = set to true after finding the library
#
# find_package( HDF5 COMPONENTS CXX ) will seeek CXX version of the
#  library as well and will set HDF5_FOUND to false if not found
#  even if C vesion is found.
#

set( HDF5_DETAIL_C_HEADER H5public.h )
set( HDF5_DETAIL_CXX_HEADER H5Cpp.h )
set( HDF5_DETAIL_C_LIBRARY hdf5 )
set( HDF5_DETAIL_CXX_LIBRARY hdf5_cpp )


set( HDF5_LIBRARIES ) 
set( HDF5_INCLUDE_DIRS ) 


foreach( HDF5_x C ${HDF5_FIND_COMPONENTS} )

  string(TOUPPER ${HDF5_x} HDF5_x)
  find_library( HDF5_${HDF5_x}_LIBRARY ${HDF5_DETAIL_${HDF5_x}_LIBRARY} )
  find_path(HDF5_${HDF5_x}_INCLUDE_DIR ${HDF5_DETAIL_${HDF5_x}_HEADER} )

  list(APPEND HDF5_LIBRARIES ${HDF5_${HDF5_x}_LIBRARY} )
  list(APPEND HDF5_INCLUDE_DIRS ${HDF5_${HDF5_x}_INCLUDE_DIR} )
  mark_as_advanced( HDF5_${HDF5_x}_LIBRARY HDF5_${HDF5_x}_INCLUDE_DIR )
endforeach( HDF5_x C ${HDF5_FIND_COMPONENTS} )

#
# Process variables
#
include(LibFindMacros)
set(HDF5_PROCESS_LIBS HDF5_LIBRARIES )
set(HDF5_PROCESS_INCLUDES HDF5_INCLUDE_DIRS)
libfind_process(HDF5)
