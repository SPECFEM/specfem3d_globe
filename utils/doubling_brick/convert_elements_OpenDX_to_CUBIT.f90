
program convert_elements

implicit none

integer, parameter :: NELEM = 32

integer i,i1,i2,i3,i4,i5,i6,i7,i8

  do i = 1,NELEM
! point order is inverted in OpenDX and numbers start at 0, not 1
    read(*,*) i4,i1,i8,i5,i3,i2,i7,i6
    write(*,230) i,i1+1,i2+1,i3+1,i4+1,i5+1,i6+1,i7+1,i8+1
  enddo

 230 format(i4,',',i4,',',i4,',',i4,',',i4,',',i4,',',i4,',',i4,',',i4)

end program convert_elements

