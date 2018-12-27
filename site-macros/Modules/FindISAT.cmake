#
# this module look for ISAT
# it will define the following values
#
# ISAT_TAB_LIBRARY
# ISAT_CK_LIBRARY
# ISAT_LIBRARIES
# ISAT_FOUND
#


set(paths "$ENV{HOME}/lib" /usr/lib /usr/local/lib )

find_library( ISAT_TAB_LIBRARY isatab ${paths} )
find_library( ISAT_CK_LIBRARY isatck ${paths} )

set( ISAT_LIBRARIES ${ISAT_CK_LIBRARY} ${ISAT_TAB_LIBRARY} )

#
# Set Process variables
#
include(LibFindMacros)
set(ISAT_PROCESS_LIBS ISAT_LIBRARIES )
libfind_process(ISAT)
