subroutine sub()
  print*, "this is fortran routine"
end subroutine sub


subroutine sub2fun()
  print *, "this is fortran routine calling c"
  call fun()
end subroutine sub2fun
