
program convert_points

implicit none

integer, parameter :: NPOIN = 67

integer i
real a,b,c

  do i = 1,NPOIN
    read(*,*) a,b,c
    write(*,*) i,',',a,',',b,',',c
  enddo

end program convert_points

