!
!--- read model 1066a and output it in Fortran90 format
!

  program read_model_1066a

  implicit none

! number of layers in DATA/1066a/1066a.dat
  integer, parameter :: NR_1066A = 160

! these are the arrays we need
  double precision radius_1066a(NR_1066A),density_1066a(NR_1066A)
  double precision vp_1066a(NR_1066A),vs_1066a(NR_1066A)
  double precision Qkappa_1066a(NR_1066A),Qmu_1066a(NR_1066A)

  integer i

  character(len=100) texttowrite

! 1066a layercake model
  open(unit=10,file='1066a.dat',status='old',action='read')

  do i=1,NR_1066A

! depth: m
! density density_1066a: g/m^3
! compressional wave speed vp_1066a: m/s
! shear wave speed vs_1066a: m/s

    read(10,*) radius_1066a(i),density_1066a(i),vp_1066a(i),vs_1066a(i),Qkappa_1066a(i),Qmu_1066a(i)
! make kg/m^3
    density_1066a(i) = density_1066a(i)/1.0d03
! make m/s
    vp_1066a(i) = vp_1066a(i)/1.0d03
    vs_1066a(i) = vs_1066a(i)/1.0d03

  enddo

  close(10)

  print *

  do i=1,NR_1066A
    write(texttowrite,201) i
    print *,trim(texttowrite),radius_1066a(i)
  enddo
 201 format(' radius_1066a(',i3,') = ')
  print *

  do i=1,NR_1066A
    write(texttowrite,202) i
    print *,trim(texttowrite),density_1066a(i)
  enddo
 202 format(' density_1066a(',i3,') = ')
  print *

  do i=1,NR_1066A
    write(texttowrite,203) i
    print *,trim(texttowrite),vp_1066a(i)
  enddo
 203 format(' vp_1066a(',i3,') = ')
  print *

  do i=1,NR_1066A
    write(texttowrite,204) i
    print *,trim(texttowrite),vs_1066a(i)
  enddo
 204 format(' vs_1066a(',i3,') = ')
  print *

  do i=1,NR_1066A
    write(texttowrite,205) i
    print *,trim(texttowrite),Qkappa_1066a(i)
  enddo
 205 format(' Qkappa_1066a(',i3,') = ')
  print *

  do i=1,NR_1066A
    write(texttowrite,206) i
    print *,trim(texttowrite),Qmu_1066a(i)
  enddo
 206 format(' Qmu_1066a(',i3,') = ')
  print *

  end program read_model_1066a

