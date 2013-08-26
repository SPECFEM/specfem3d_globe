!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  program convert_epicenter_OpenDX
!
!--- convert the coordinates of the earthquake epicenter to OpenDX format
!

  implicit none

  include "constants_topo.h"

! small offset to add to the epicenter glyph (in meters)
  double precision, parameter :: SMALL_OFFSET = 125000.d0

! filter final surface using box filter
  logical, parameter :: SMOOTH_THE_MODEL = .true.
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 5

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)
  integer ibathy_topo_ori(NX_BATHY,NY_BATHY)

  integer ix,iy,minvalue,maxvalue
  integer ix_current,iy_current,ix_min,ix_max,iy_min,iy_max,ix_value,iy_value

  double precision value_sum,area_window
  double precision :: phi,theta,r,x,y,z,elevation

  double precision :: long,lat

!----

! latitude and longitude of the epicenter
  lat = 31.4400d0
  long = 104.1000d0

!----

! read the topography file
  print *,'NX_BATHY,NY_BATHY = ',NX_BATHY,NY_BATHY
  print *,'file used has a resolution of ',RESOLUTION_TOPO_FILE,' minutes'

  print *
  print *,'reading topo file'

  open(unit=13,file='topo_bathy_etopo4_from_etopo2_subsampled.dat',status='old')
  do iy=1,NY_BATHY
    do ix=1,NX_BATHY
      read(13,*) ibathy_topo_ori(ix,iy)
    enddo
  enddo
  close(13)

! compute min and max before smoothing
  minvalue = minval(ibathy_topo_ori)
  maxvalue = maxval(ibathy_topo_ori)
  print *,'min and max of topography before smoothing = ',minvalue,maxvalue

!----

! smooth topography/bathymetry model
  if(SMOOTH_THE_MODEL) then

  print *
  print *,'smoothing topo file'
  if(SIZE_FILTER_ONE_SIDE < 1) stop 'SIZE_FILTER_ONE_SIDE must be greater than 1 for filter'
  print *,'size of window filter is ',2*SIZE_FILTER_ONE_SIDE+1,' x ',2*SIZE_FILTER_ONE_SIDE+1
  area_window = dble((2*SIZE_FILTER_ONE_SIDE+1)**2)

  do iy_current = 1,NY_BATHY

   if(mod(iy_current,10) == 0) print *,'smoothing line ',iy_current,' out of ',NY_BATHY

    do ix_current = 1,NX_BATHY

      value_sum = 0.d0

! compute min and max of window
      ix_min = ix_current - SIZE_FILTER_ONE_SIDE
      ix_max = ix_current + SIZE_FILTER_ONE_SIDE

      iy_min = iy_current - SIZE_FILTER_ONE_SIDE
      iy_max = iy_current + SIZE_FILTER_ONE_SIDE

! loop on points in window to compute sum
      do iy = iy_min,iy_max
      do ix = ix_min,ix_max

! copy current value
        ix_value = ix
        iy_value = iy

! avoid edge effects, use periodic boundary in Xmin and Xmax
      if(ix_value < 1) ix_value = ix_value + NX_BATHY
      if(ix_value > NX_BATHY) ix_value = ix_value - NX_BATHY

! avoid edge effects, use rigid boundary in Ymin and Ymax
! *not* periodic, because South and North poles must not be merged
      if(iy_value < 1) iy_value = 1
      if(iy_value > NY_BATHY) iy_value = NY_BATHY

! compute sum
      value_sum = value_sum + dble(ibathy_topo_ori(ix_value,iy_value))

      enddo
      enddo

! assign mean value to filtered point
      ibathy_topo(ix_current,iy_current) = nint(value_sum / area_window)

    enddo
  enddo

  else

! no smoothing
    ibathy_topo = ibathy_topo_ori

  endif

!----

! compute min and max after smoothing
  minvalue = minval(ibathy_topo)
  maxvalue = maxval(ibathy_topo)
  print *,'min and max of topography after smoothing = ',minvalue,maxvalue

! store file in OpenDX format

! create the OpenDX file
  print *
  print *,'creating OpenDX file'

! creating the header
  open(unit=11,file='epicenter.dx',status='unknown')

! write point coordinates
  write(11,*) 'object 1 class array type float rank 1 shape 3 items 1 data follows'
!! DK DK in principle this should take the ellipticity of the Earth into account
  phi = long * PI / 180.d0
  theta = (90.d0 - lat) * PI / 180.d0

! use an Earth model of mean radius equal to 1, and add elevation to each point, amplified in order to see it
  call get_topo_bathy(lat,long,elevation,ibathy_topo)
  r = (R_EARTH + EXAGGERATION_FACTOR_TOPO * elevation + SMALL_OFFSET) / R_EARTH

  x = r*sin(theta)*cos(phi)
  y = r*sin(theta)*sin(phi)
  z = r*cos(theta)

  write(11,*) sngl(x),sngl(y),sngl(z)

! write elements
  write(11,*) 'object 2 class array type int rank 1 shape 2 items 1 data follows'
! point numbers start at 0 in OpenDX, not 1 as in AVS
  write(11,*) '0 0'

  write(11,*) 'attribute "element type" string "lines"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items 1 data follows'

! write data values (size of the glyph, therefore do NOT set it to 0)
  write(11,*) '2'

! define OpenDX field
  write(11,*) 'attribute "dep" string "positions"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

  close(11)

  end program convert_epicenter_OpenDX

!----------------------------

  subroutine get_topo_bathy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

  implicit none

  include "constants_topo.h"

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

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

