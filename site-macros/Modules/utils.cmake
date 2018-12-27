# ==============================================================================
#  Miscellaneous utilities 
# ==============================================================================

# ==============================================================================
#  list_stringize(strlist item1 item2 ...) 
#      Convert a list of items into list of stringified items
# ==============================================================================

macro( list_stringize strlist )
   set(args ${ARGV})
   list(REMOVE_AT args 0)
   foreach( item ${args} )
      list(APPEND ${strlist} "\"${item}\"")
   endforeach( item )
endmacro( list_stringize )


# ==============================================================================
#  SET_BUILD_TYPE( conf )
#    set the default CMAKE_BUILD_TYPE
# ==============================================================================
macro( set_build_type conf )
  if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE ${conf} CACHE STRING
      "Choose the type of build, options are: <None> ${CMAKE_CONFIGURATION_TYPES}"
      FORCE
      )
  endif (NOT CMAKE_BUILD_TYPE)
endmacro( set_build_type )
