!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

  program smooth_topo_bathy_PPM_image
!
!--- read topography and bathymetry ASCII file, smooth it and create PNM image to check it
!

! this program uses a very simple window filter
! it works fine because the model has a very high resolution
! in principle it would be better to use a circular cap (i.e. a Gaussian filter) to take
! the notion of angular distance between points into account in the filter
! in practice though this simple filter works perfectly fine

  implicit none

!---  ETOPO4 4-minute model created by subsampling and smoothing etopo-2
! size of topography and bathymetry file
! integer, parameter :: NX_BATHY = 5400,NY_BATHY = 2700

!--- ETOPO2 2-minute model, not implemented yet
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY = 10800,NY_BATHY = 5400

!--- ETOPO1 1-minute model, implemented now, but data file must be created first
! size of topography and bathymetry file
! integer, parameter :: NX_BATHY = 21600,NY_BATHY = 10800

! filter final surface using box filter
  logical, parameter :: SMOOTH_THE_MODEL = .true.
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 3 !! 7

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)
  integer ibathy_topo_ori(NX_BATHY,NY_BATHY)

  integer ix,iy,minvalue,maxvalue
  integer ix_current,iy_current,ix_min,ix_max,iy_min,iy_max,ix_value,iy_value

  double precision value_sum,area_window

!----

! read the topography file
  print *,'NX_BATHY,NY_BATHY = ',NX_BATHY,NY_BATHY

  print *
  print *,'reading topo file'

  open(unit=13,file='topo_bathy_etopo2v2c_original_unmodified_unsmoothed.dat',status='old')
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
  if (SMOOTH_THE_MODEL) then

  print *
  print *,'smoothing topo file'
  if (SIZE_FILTER_ONE_SIDE < 1) stop 'SIZE_FILTER_ONE_SIDE must be greater than 1 for filter'
  print *,'size of window filter is ',2*SIZE_FILTER_ONE_SIDE+1,' x ',2*SIZE_FILTER_ONE_SIDE+1
  area_window = dble((2*SIZE_FILTER_ONE_SIDE+1)**2)

  do iy_current = 1,NY_BATHY

   if (mod(iy_current,10) == 0) print *,'smoothing line ',iy_current,' out of ',NY_BATHY

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
      if (ix_value < 1) ix_value = ix_value + NX_BATHY
      if (ix_value > NX_BATHY) ix_value = ix_value - NX_BATHY

! avoid edge effects, use rigid boundary in Ymin and Ymax
! *not* periodic, because South and North poles must not be merged
      if (iy_value < 1) iy_value = 1
      if (iy_value > NY_BATHY) iy_value = NY_BATHY

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

! save the smoothed model
  if (SMOOTH_THE_MODEL) then
    print *
    print *,'saving the smoothed model'
    open(unit=13,file='topo_bathy_etopo2v2c_smoothed_window_3.dat',status='unknown')
    do iy=1,NY_BATHY
      do ix=1,NX_BATHY
        write(13,*) ibathy_topo(ix,iy)
      enddo
    enddo
    close(13)
    print *,'can suppress white spaces in filtered model to save space if needed'
  endif

! create image with grey levels
  do iy = 1,NY_BATHY
    do ix = 1,NX_BATHY
      ibathy_topo(ix,iy) = 255 * (ibathy_topo(ix,iy) - minvalue) / (maxvalue - minvalue)
      if (ibathy_topo(ix,iy) < 1) ibathy_topo(ix,iy) = 1
      if (ibathy_topo(ix,iy) > 255) ibathy_topo(ix,iy) = 255
    enddo
  enddo

! store image in PNM format with grey levels

! create the PNM image
  print *
  print *,'creating PNM image'

! creating the header
  open(unit=27,file='image_topo_bathy.pnm',status='unknown')
  write(27,100)
  write(27,102) NX_BATHY,NY_BATHY
  write(27,103)

 100 format('P3')
 102 format(i6,' ',i6)
 103 format('255')
 104 format(i3)

  do iy = 1,NY_BATHY
  do ix = 1,NX_BATHY

! write data value (red = green = blue to produce grey levels)
      write(27,104) ibathy_topo(ix,iy)
      write(27,104) ibathy_topo(ix,iy)
      write(27,104) ibathy_topo(ix,iy)

  enddo
  enddo

  close(27)

  end program smooth_topo_bathy_PPM_image

