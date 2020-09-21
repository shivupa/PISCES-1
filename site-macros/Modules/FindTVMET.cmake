#
# this module look for TVMET support
# it will define the following values
#
# TVMET_INCLUDE_DIRS  = where tvmet/tvmet.h can be found
# TVMET_FOUND        = set to true after finding the library
#

find_path( TVMET_INCLUDE_DIRS tvmet/tvmet.h HINTS ENV CPATH )

#
# Process variables
#
include(LibFindMacros)
set(TVMET_PROCESS_INCLUDES TVMET_INCLUDE_DIRS)
libfind_process(TVMET)
