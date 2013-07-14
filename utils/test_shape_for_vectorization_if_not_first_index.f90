
program test

!! DK DK July 2013: test if we can use (ijk,1,1) instead of (i,j,k) when the group is not located at the first index
!! DK DK July 2013:  (do NOT compile with range checking of course)

implicit none

integer, parameter :: NX = 5
integer, parameter :: NY = 6
integer, parameter :: NZ = 7
integer, parameter :: NL = 4

integer, dimension(NL,NX,NY,NZ) :: a

integer :: i,j,k,l,ijk

do l = 1,NL
do k = 1,NZ
do j = 1,NY
do i = 1,NX
  a(l,i,j,k) = l + i + j + k    ! create test values
  print *,a(l,i,j,k)
enddo
enddo
enddo
enddo

print *
!! DK DK in practice it gives the exact same order, thus the trick works fine
print *,'now in different order'
print *

do l = 1,NL
do ijk = 1,NX*NY*NZ
  print *,a(l,ijk,1,1)
enddo
enddo

end

