# - Try to find hwloc
# Once done this will define
#  hwloc_FOUND - System has hwloc
#  hwloc_INCLUDE_DIRS - The hwloc include directories
#  hwloc_LIBRARIES - The libraries needed to use hwloc


find_path(hwloc_INCLUDE_DIR hwloc.h
          HINTS ENV HWLOC_ROOT
          PATH_SUFFIXES include 
		  )

find_library(hwloc_LIBRARY NAMES hwloc
             HINTS ENV HWLOC_ROOT
			 PATH_SUFFIXES lib 
			 )

set(hwloc_LIBRARIES ${hwloc_LIBRARY} )
set(hwloc_INCLUDE_DIRS ${hwloc_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set hwloc_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(hwloc  DEFAULT_MSG
                                  hwloc_LIBRARY hwloc_INCLUDE_DIR)

mark_as_advanced(hwloc_INCLUDE_DIR hwloc_LIBRARY )