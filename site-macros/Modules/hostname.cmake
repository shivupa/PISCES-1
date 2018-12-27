#
# Determine the host name (particular name in network, not CMAKE_HOST_SYSTEM_NAME)
# Sets variable HOSTNAME to the name of the system, or "unknown" if not possible. 
#

set(HOSTNAME unknown)
find_program( HOSTNAME_EXE hostname )
if( HOSTNAME_EXE )
  execute_process( COMMAND ${HOSTNAME_EXE} OUTPUT_VARIABLE HOSTNAME )
  string( STRIP ${HOSTNAME} HOSTNAME )
endif( HOSTNAME_EXE )

set(HOSTNAME ${HOSTNAME} CACHE STRING "Name of the current system")
