!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_topo_bathy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

  implicit none

  include "constants.h"

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

  double precision xlat,xlon,value

  integer iadd1,iel1
  double precision samples_per_degree_topo
  double precision xlo

  xlo = xlon
  if(xlon < 0.d0) xlo = xlo + 360.d0

! compute number of samples per degree
  samples_per_degree_topo = dble(RESOLUTION_TOPO_FILE) / 60.d0

! compute offset in data file and avoid edge effects
  iadd1 = 1 + int((90.d0-xlat)/samples_per_degree_topo)
  if(iadd1 < 1) iadd1 = 1
  if(iadd1 > NY_BATHY) iadd1 = NY_BATHY

  iel1 = int(xlo/samples_per_degree_topo)
  if(iel1 <= 0 .or. iel1 > NX_BATHY) iel1 = NX_BATHY

! convert integer value to double precision
  value = dble(ibathy_topo(iel1,iadd1))

  end subroutine get_topo_bathy

! -------------------------------------------

  subroutine read_topo_bathy_file(ibathy_topo)
!
!---- read topography and bathymetry file once and for all
!
  implicit none

  include "constants.h"

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

  integer itopo_x,itopo_y

! use either Etopo-5 or Etopo-2
! Harvard Etopo-5 is a smoothed model
  if(RESOLUTION_TOPO_FILE == 5) then
    open(unit=13,file='DATA/topo_bathy/topo_bathy_etopo5.dat',status='old')
! ETOPO-2 not implemented yet
!  else if(RESOLUTION_TOPO_FILE == 2) then
!    open(unit=13,file='DATA/topo_bathy/topo_bathy_etopo2.dat',status='old')
  else
    stop 'incorrect resolution selected for topography/bathymetry file'
  endif

  do itopo_y=1,NY_BATHY
    do itopo_x=1,NX_BATHY

      read(13,*) ibathy_topo(itopo_x,itopo_y)

! impose maximum depth of oceans, to suppress oscillations near deep trenches
  if (USE_MAXIMUM_DEPTH_OCEANS .and. ibathy_topo(itopo_x,itopo_y) < MAXIMUM_DEPTH_OCEANS) &
    ibathy_topo(itopo_x,itopo_y) = MAXIMUM_DEPTH_OCEANS

    enddo
  enddo

  close(13)

  end subroutine read_topo_bathy_file

