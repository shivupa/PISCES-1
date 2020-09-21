!
! Multiply a matrix by vector using  DGEMV
!
program simpleblas
implicit none
real(8) :: a(10,10),x(10),y(10)
real(8) :: alpha=1, beta=0
integer i

a = 0
x = 0
y = 0
do i = 1,10
   a(i,i) = 1
   x(i) = i
enddo

call dgemv('n', 10, 10, alpha, a, 10, x, 1, beta, y, 1)

if (all( x  == y )) then
   print *, "success"
else
   print *, "unexpected matrix-vector result"
   stop 1
endif
end program simpleblas
