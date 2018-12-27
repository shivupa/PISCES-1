#
# this module look for FFTW  www.fftw.org support
#
# FFTW_PATH         = directory to find libfftw3.a libfftw3_omp.a
#
# it will define the following values
#
# FFTW_LIBRARY      = the library to link against
# FFTW_FOUND        = set to true after finding the library
#

find_path(FFTW_INCLUDES NAMES fftw3.h
   HINTS ${FFTW_PATH})
find_library(FFTW_LIBRARY NAMES fftw3
  HINTS ${FFTW_PATH})
find_library(FFTW_OMP_LIBRARY NAMES fftw3_omp
  HINTS ${FFTW_PATH})

set(FFTW_FOUND FALSE)
if(FFTW_LIBRARY)
  set(FFTW_FOUND TRUE)
  mark_as_advanced(FFTW_LIBRARY)
endif(FFTW_LIBRARY)

if(FFTW_FOUND)
  message(STATUS "Found FFTW : ${FFTW_LIBRARY}")
else(FFTW_FOUND)
  message(FATAL_ERROR "Could not find fftw")
endif(FFTW_FOUND)
