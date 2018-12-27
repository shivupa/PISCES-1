#
# this module look for ARPACK (http://hdf.ncsa.uiuc.edu) support
#
# ARPACK_PATH         = directory to find libarpack.a
#
# it will define the following values
#
# ARPACK_LIBRARY      = the library to link against
# ARPACK_FOUND        = set to true after finding the library
#

find_library(ARPACK_LIBRARY NAMES arpack
  HINTS ${ARPACK_PATH})

set(ARPACK_FOUND FALSE)
if(ARPACK_LIBRARY)
  set(ARPACK_FOUND TRUE)
  mark_as_advanced(ARPACK_LIBRARY)
endif(ARPACK_LIBRARY)

if(ARPACK_FOUND)
  message(STATUS "Found ARPACK : ${ARPACK_LIBRARY}")
else(ARPACK_FOUND)
  message(FATAL_ERROR "Could not find Arpack")
endif(ARPACK_FOUND)
