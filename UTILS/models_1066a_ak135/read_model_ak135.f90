!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program read_model_ak135

  implicit none

! number of layers in DATA/ak135/ak135.dat
  integer, parameter :: NR_AK135 = 144

! these are the values we need
  double precision radius_ak135(NR_AK135),density_ak135(NR_AK135)
  double precision vp_ak135(NR_AK135),vs_ak135(NR_AK135)
  double precision Qkappa_ak135(NR_AK135),Qmu_ak135(NR_AK135)

  integer i

  double precision depth

! R_EARTH is the radius of the bottom of the oceans (radius of Earth in m)
  double precision, parameter :: R_EARTH = 6371000.d0

  character(len=100) texttowrite

! ak135 layercake model
  open(unit=10,file='ak135.dat',status='old',action='read')

  do i=NR_AK135,1,-1

! depth: km
! density density_ak135: kg/m^3
! compressional wave speed vp_ak135: km/s
! shear wave speed vs_ak135: km/s

    read(10,*) depth,density_ak135(i),vp_ak135(i),vs_ak135(i),Qkappa_ak135(i),Qmu_ak135(i)

! get radius in meters
    radius_ak135(i)=R_EARTH-depth*1.0d03

  enddo
  close(10)

  do i=NR_AK135-4,NR_AK135
! strip the ocean and sediment layers and replace them by upper crust
    density_ak135(i) = density_ak135(NR_AK135-5)
    vp_ak135(i) = vp_ak135(NR_AK135-5)
    vs_ak135(i) = vs_ak135(NR_AK135-5)
    Qkappa_ak135(i) = Qkappa_ak135(NR_AK135-5)
    Qmu_ak135(i) = Qmu_ak135(NR_AK135-5)
  enddo

  print *

  do i=1,NR_AK135
    write(texttowrite,201) i
    print *,trim(texttowrite),radius_ak135(i)
  enddo
 201 format(' radius_ak135(',i3,') = ')
  print *

  do i=1,NR_AK135
    write(texttowrite,202) i
    print *,trim(texttowrite),density_ak135(i)
  enddo
 202 format(' density_ak135(',i3,') = ')
  print *

  do i=1,NR_AK135
    write(texttowrite,203) i
    print *,trim(texttowrite),vp_ak135(i)
  enddo
 203 format(' vp_ak135(',i3,') = ')
  print *

  do i=1,NR_AK135
    write(texttowrite,204) i
    print *,trim(texttowrite),vs_ak135(i)
  enddo
 204 format(' vs_ak135(',i3,') = ')
  print *

  do i=1,NR_AK135
    write(texttowrite,205) i
    print *,trim(texttowrite),Qkappa_ak135(i)
  enddo
 205 format(' Qkappa_ak135(',i3,') = ')
  print *

  do i=1,NR_AK135
    write(texttowrite,206) i
    print *,trim(texttowrite),Qmu_ak135(i)
  enddo
 206 format(' Qmu_ak135(',i3,') = ')
  print *

  end program read_model_ak135

