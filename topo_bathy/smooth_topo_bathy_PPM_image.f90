!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program smooth_topo_bathy_PPM_image
!
!--- read topography and bathymetry ASCII file, smooth it and create PNM image to check it
!

!! DK DK this program uses a very simple window filter
!! DK DK it works fine because the model has a very high resolution
!! DK DK in principle it would be better to use a circular cap to take
!! DK DK the notion of angular distance between points into account in the filter
!! DK DK in practice though this simple filter works perfectly fine

  implicit none

  include "../../constants.h"

! filter final surface using box filter
  logical, parameter :: SMOOTH_THE_MODEL = .true.
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 7

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)
  integer ibathy_topo_ori(NX_BATHY,NY_BATHY)

  integer ivalue,icurrent_rec,ix,iy
  integer minvalue,maxvalue
  integer ix_current,iy_current,ix_min,ix_max,iy_min,iy_max,ix_value,iy_value

  double precision value_sum,area_window

!----

! read the topography file
  print *,'NX_BATHY,NY_BATHY = ',NX_BATHY,NY_BATHY
  print *,'file used has a resolution of ',RESOLUTION_TOPO_FILE,' minutes'

  print *
  print *,'reading topo file'

  open(unit=13,file='topo_bathy_etopo4_ori_from_etopo2_subsampled.dat',status='old')
  do iy=1,NY_BATHY
    do ix=1,NX_BATHY
      read(13,*) ibathy_topo_ori(ix,iy)
    enddo
  enddo
  close(13)

!! DK DK compute min and max before smoothing
  minvalue = minval(ibathy_topo_ori)
  maxvalue = maxval(ibathy_topo_ori)
  print *,'min and max of topography before smoothing = ',minvalue,maxvalue

!----

!! DK DK smooth topography/bathymetry model

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

! avoid edge effects, use periodic boundary
      if(ix_value < 1) ix_value = ix_value + NX_BATHY
      if(ix_value > NX_BATHY) ix_value = ix_value - NX_BATHY

      if(iy_value < 1) iy_value = iy_value + NY_BATHY
      if(iy_value > NY_BATHY) iy_value = iy_value - NY_BATHY

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

!! DK DK compute min and max after smoothing
  minvalue = minval(ibathy_topo)
  maxvalue = maxval(ibathy_topo)
  print *,'min and max of topography after smoothing = ',minvalue,maxvalue

!! DK DK save the smoothed model
  if(SMOOTH_THE_MODEL) then
    print *
    print *,'saving the smoothed model'
    open(unit=13,file='topo_bathy_etopo4_smoothed_window7.dat',status='unknown')
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
      if(ibathy_topo(ix,iy) < 1) ibathy_topo(ix,iy) = 1
      if(ibathy_topo(ix,iy) > 255) ibathy_topo(ix,iy) = 255
    enddo
  enddo

! store image in PNM format with grey levels

! create the PNM image
  print *
  print *,'creating PNM image'

! creating the header
  open(unit=27,file='image_topo_bathy.pnm',status='unknown')
  write(27,100)
  write(27,101)
  write(27,102) NX_BATHY,NY_BATHY
  write(27,103)

 100 format('P3')
 101 format('# creator DK')
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

