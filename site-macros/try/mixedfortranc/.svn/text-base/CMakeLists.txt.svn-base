
project( TestBLASLAPACK C Fortran )
cmake_minimum_required(VERSION 2.8)


add_library( clib showmsg.c )
add_library( flib showmsg.f90 )

add_executable( fcallc fcallc.f90 )
add_executable( ccallf ccallf.c   )

target_link_libraries( fcallc clib )
target_link_libraries( ccallf flib )


# TODO do not assume name mangling is froutine_

