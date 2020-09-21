#
# FORTRANRTL ===================================================================
#     Find Fortran Runtime Libraries
# 
# Supported compilers [todo: and platforms]:
#    * Intel  
#    * PGI
#    * [todo: G77]
#    * [todo: GFortran]
#    * [todo: G95]
#
# Once called, this module will define
#
#    FORTRANRTL_FOUND      
#         (not cached) successfully located required libraries
#
#    FORTRANRTL_LIBRARIES  
#         (not cached) list of libs to link
#
#    FORTRANRTL_LIBRARY_NAMES
#         (cached) This is a list of library names (not paths) that will be
#         found. For the supported compilers, it is automatically populated and
#         to be edited only if something goes wrong. If you are using an
#         UNSUPPORTED compiler, you must specify the correct names (again, not
#         paths!) separated by space or semi-colon.
#
#    FORTRANRTL_LIBRARY_PATHS
#         (cached) Hint paths to search for libraries in the list  
#         FORTRANRTL_LIBARARY_NAMES. Some "magic" locations are also searched. 
#
#    FORTRANRTL_LIBRARY_<lib>  
#         (cached) Location of the RTL with name <lib>.
#
# Following variables will affect the way this module works
#
#    FORTRANRTL_LINKCHECK  [todo]
#         Set it to TRUE in order compile-and-link check a small mixed language
#         app. Set to FALSE if fails, or by default.
#
#    FORTRANRTL_COMPILER
#         Set to the ID of the compiler to get the RTL of.  Defaults to
#         CMAKE_Fortran_COMPILER_ID. Currently supported are Intel, PGI, GNU,
#         G95.
#
# Note: In  very very rare cases (which  I personally never ran  into) you might
# need  RTL of  multiple  fortran compilers.  You  can NOT  simply redefine  the
# FORTRANRTL_COMPILER and  call this module,  since variables from a  first call
# will  be cached  and won't  be replaced.  Instead, you  should simply  set the
# FORTRANRTL_COMPILER variable  to a  dummy value (you  cannot not set  it), and
# edit the FORTRANRTL_LIBRARY_NAMES and FotranRTL_LIBRARY_PATHS yourself.
#
# Let know if something goes wrong, using the CMAKE user group of course, 
# not my email! ;) Flames not welcome >:|  
#
# Copyright (c) 2008-2009 S. Levent Yilmaz / leventyilmaz AT gmail DOT com
#
# THANKS: 
#
# * Tronic for LibFindMacros (see http://www.vtk.org/Wiki/CMake:How_To_Find_Libraries)
# 
#
# DETAILS:
#
# o Predefined Lists 
#
#    * Intel: ifport ifcode
#    * PGI  : pgf90  pgf90_rpm1  pgf902  pgf90rtl  pgftnrtl 
#
# o Usage Errors
#
#    * {todo:Describe not setting any compiler or not enabling fortran language}
#


# ==============================================================================
#  set module variates
# ==============================================================================

set( FORTRANRTL_DETAIL_LIBLIST_Intel_Windows libifcoremd )
set( FORTRANRTL_DETAIL_LIBLIST_Intel_Linux  ifcore ifport)
set( FORTRANRTL_DETAIL_LIBLIST_PGI_Linux pgf90 pgf90_rpm1 pgf902 pgf90rtl pgftnrtl pgc )
set( FORTRANRTL_DETAIL_LIBLIST_GNU gfortran )


# ==============================================================================
# Let the fun begin!  ========================================================== 
# ==============================================================================
function(showme x)
   message(STATUS "${x} is ${${x}}")
endfunction(showme)

# ==============================================================================
#   SET the target Fortran Compiler
# ==============================================================================
#
# First call, and caller did not set the compiler manually
# (it is expected that Fortran language is enabled)
#
if (NOT FORTRANRTL_COMPILER)
   # 
   macro(setcom val)
     set(FORTRANRTL_COMPILER ${val} CACHE STRING
        "Fortran compiler to get the runtime libraries of" FORCE )
   endmacro(setcom)
    
   
   # Put compiler setting to cache
   setcom(NOTFOUND)
   mark_as_advanced( FORTRANRTL_COMPILER )


   # Try detecting from preset ID
   if (CMAKE_Fortran_COMPILER_ID)
      setcom(${CMAKE_Fortran_COMPILER_ID})
   
   # Try detecting from compiler path
   elseif (CMAKE_Fortran_COMPILER)
      get_filename_component(x ${CMAKE_Fortran_COMPILER} NAME_WE)
      
      # try a few defaults
      if (x MATCHES "ifort")
         setcom(Intel)
      elseif (CMAKE_COMPILER_IS_GNUG77)
         setcom(G77)
      endif(x MATCHES "ifort")

   endif(CMAKE_Fortran_COMPILER_ID)
   
   # Failed to detect
   # It is an error, see [DETAILS/Usage Errors]
   if (NOT FORTRANRTL_COMPILER)
      message(SEND_ERROR "FORTRANRTL: Can not find a Fortran compiler, or a compiler not set.")
      return()
   endif (NOT FORTRANRTL_COMPILER)

   # Append platform
   setcom( ${FORTRANRTL_COMPILER}_${CMAKE_Fortran_PLATFORM_ID} )
   
endif(NOT FORTRANRTL_COMPILER)



message(STATUS "Fortran run-time using ${FORTRANRTL_COMPILER}")





# ==============================================================================
#   LIBRARY Names
# ==============================================================================

# Set LIBRARY_NAMES
#   put to cache
set(FORTRANRTL_LIBRARY_NAMES ${FORTRANRTL_DETAIL_LIBLIST_${FORTRANRTL_COMPILER}} 
   CACHE STRING "Fortran runtime library names, not paths, separated by space.")
mark_as_advanced( FORTRANRTL_LIBRARY_NAMES )

message(STATUS "Fortran run-time using libs ${FORTRANRTL_LIBRARY_NAMES}")


# ==============================================================================
#   Possible LIBRARY PATH HINTS
# ==============================================================================
if( NOT FORTRANRTL_LIBRARY_PATHS )
  
  file(TO_CMAKE_PATH 
    "$ENV{LIBRARY_PATH};$ENV{LD_LIBRARY_PATH}"
    paths
    )
  
  set(FORTRANRTL_LIBRARY_PATHS  ${paths}
    CACHE 
    FILENAME 
    "Fortran runtime paths as hints, separated by semi-colon."
    FORCE
    )
  mark_as_advanced( FORTRANRTL_LIBRARY_PATHS )
endif( NOT FORTRANRTL_LIBRARY_PATHS )
  

#
# Locate the libraries
# 

# clean up
foreach(x ${FORTRANRTL_DETAIL_libnames})
   unset(FORTRANRTL_LIBRARY_${x} CACHE)
endforeach(x ${FORTRANRTL_DETAIL_libnames})
set( FORTRANRTL_DETAIL_libnames ${FORTRANRTL_LIBRARY_NAMES} CACHE INTERNAL "" FORCE )

# find and populate FORTRANRTL_LIBRARIES
foreach( x ${FORTRANRTL_LIBRARY_NAMES} )
  #file(TO_CMAKE_PATH ${FORTRANRTL_LIBRARY_PATHS} paths)
  find_library( FORTRANRTL_LIBRARY_${x}
    NAMES ${x}
    HINTS ${FORTRANRTL_LIBRARY_PATHS}
    )
  list(APPEND FORTRANRTL_LIBRARIES ${FORTRANRTL_LIBRARY_${x}} )
  mark_as_advanced( FORTRANRTL_LIBRARY_${x} )
endforeach( x ${FORTRANRTL_LIBRARY_NAMES} )

#
# Quirks
#
if(FORTRANRTL_COMPILER STREQUAL "Intel_Windows")
   list(GET FORTRANRTL_LIBRARY_NAMES 0 lib)
   if( FORTRANRTL_LIBRARY_${lib} )
      get_filename_component(dir ${FORTRANRTL_LIBRARY_${lib}} PATH)
      link_directories(${dir})
#      list(APPEND FORTRANRTL_LINK_FLAGS_ ${CMAKE_LIBRARY_PATH_FLAG}${dir})
      #list(APPEND FORTRANRTL_LINK_FLAGS_DEBUG /NODEFAULTLIB:LIBCMT       
      
   endif( FORTRANRTL_LIBRARY_${lib} )   
endif(FORTRANRTL_COMPILER STREQUAL "Intel_Windows")


# store in cache to display
set( FORTRANRTL_LIBRARIES_CACHED ${FORTRANRTL_LIBRARIES} CACHE STRING 
   "FORTRANRTL_LIBRARIES variable to be referenced in your script. This is just a display, changes here will have no effect"
   FORCE)

#
# Set Process variables
#
include(LibFindMacros)
set(FORTRANRTL_PROCESS_LIBS FORTRANRTL_LIBRARIES )
libfind_process(FORTRANRTL)
