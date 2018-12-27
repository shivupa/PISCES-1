


set( CMAKE_Fortran_FLAGS "/W1 /nologo /fpp /libs:dll /threads /iface:mixed_str_len_arg" CACHE STRING "Intel Fortran flags" )
set( CMAKE_Fortran_FLAGS_DEBUG "/debug:full /dbglibs" CACHE STRING "Intel Fortran debug flags" )
set( CMAKE_Fortran_FLAGS_RELEASE "/O1 /D NDEBUG" CACHE STRING "Intel Fortran debug flags" )
