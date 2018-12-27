#
# this module look for signal handling support
# it will define the following values
#
# SIGNALHANDLING_FOUND
# 
# Currently incomplete. Assuming no additional libraries/headers 
# required. 
#

if (UNIX)
   set(SIGNALHANDLING_FOUND 1)
endif (UNIX)
