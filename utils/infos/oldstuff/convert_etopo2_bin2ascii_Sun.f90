!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  program convert_etopo2_bin2ascii
!
!---- convert etopo-2 raw binary Sun file to ASCII - use on Sun GPS machines at Caltech
!
  implicit none

!--- ETOPO2 2-minute model
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY = 10800,NY_BATHY = 5400

! use integer array to store values
  integer(kind=2) ibathy_topo(NX_BATHY,NY_BATHY)

  integer itopo_x,itopo_y

! read the topography file
  print *,'reading topo file in Sun-Sparc binary format'
  call read_topo_bathy_bin_etopo2(ibathy_topo)

  print *,'min and max of original topography = ',minval(ibathy_topo),maxval(ibathy_topo)

! save the ASCII file
!! DK DK beware, longitude blocks MUST be swapped
  stop 'beware, longitude blocks MUST be swapped, see file swap_topo_bathy_Sun.f90 for details'
  open(unit=13,file='topo_bathy_etopo2.dat',status='unknown')
  do itopo_y = 1,NY_BATHY
    do itopo_x = 1,NX_BATHY
      write(13,*) ibathy_topo(itopo_x,itopo_y)
    enddo
  enddo
  close(13)

  print *
  print *,'can filter white spaces from topo_bathy_etopo2.dat to save space'
  print *

  end program convert_etopo2_bin2ascii

! -------------------------------------------

  subroutine read_topo_bathy_bin_etopo2(ibathy_topo)
!
!---- read topography and bathymetry file once and for all
!
  implicit none

!--- ETOPO2 2-minute model
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY = 10800,NY_BATHY = 5400

! use integer array to store values
  integer(kind=2) ibathy_topo(NX_BATHY,NY_BATHY)

  integer iadd1,iel1,icurrent_rec

! this is the path on the Caltech GPS Sun network
  open(unit=13,file='/home/datalib/Topography/ETOPO-02/etopo2.raw',status='old',access='direct',recl=2)

  icurrent_rec = 1

  do iadd1 = 1,NY_BATHY
    do iel1 = 1,NX_BATHY
      read(13,rec=icurrent_rec) ibathy_topo(iel1,iadd1)
      icurrent_rec = icurrent_rec + 1
    enddo
  enddo

  close(13)

  end subroutine read_topo_bathy_bin_etopo2

