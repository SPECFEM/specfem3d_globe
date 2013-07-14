
program test_shape_vectorization

!! DK DK July 2013: test if we can use (ijk,1,1) instead of (i,j,k) when the group is not located at the first index
!! DK DK July 2013:  (do NOT compile with range checking of course)

implicit none

integer, parameter :: NX = 8
integer, parameter :: NY = 8
integer, parameter :: NZ = 5
integer, parameter :: NL = 7
integer, parameter :: NM = 4

integer, dimension(NL,NX,NY,NZ,NM) :: a

integer :: i,j,k,l,m,ijk
integer :: a_from_first_order,difference

real :: r

call random_seed()

open(unit=10,file='first_result.dat',status='unknown')
do m = 1,NM
do k = 1,NZ
do j = 1,NY
do i = 1,NX
do l = 1,NL
  call random_number(r)
  a(l,i,j,k,m) = nint(r * 20000.)  ! create test values
  print *,a(l,i,j,k,m)
  write(10,*) a(l,i,j,k,m)
enddo
enddo
enddo
enddo
enddo
close(10)

print *
!! DK DK in practice it gives the exact same order, thus the trick works fine
print *,'now in the vectorized order order'
print *

open(unit=10,file='first_result.dat',status='old')
do m = 1,NM
do ijk = 1,NX*NY*NZ
do l = 1,NL
  read(10,*) a_from_first_order
  difference = abs(a(l,ijk,1,1,m) - a_from_first_order)
  print *,a(l,ijk,1,1,m),difference
  if(difference /= 0) stop 'error, difference between the two orders is not zero'
enddo
enddo
enddo
close(10)

print *
print *,'the test is successful, the two orders are 100% identical'
print *

end program test_shape_vectorization

