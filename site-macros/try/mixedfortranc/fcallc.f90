!
! A simple fortran program that calls a C routine
!
program fcallc
  implicit none
  external showmsg
  
  
  call showmsg("fcallc success")
end program fcallc
