!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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
 integer, parameter :: NX_BATHY_4 = 5400,NY_BATHY_4 = 2700
!--- ETOPO2 2-minute model, not implemented yet
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY_2 = 10800,NY_BATHY_2 = 5400
!--- ETOPO1 1-minute model, implemented now, but data file must be created first
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY_1 = 21600,NY_BATHY_1 = 10800

  integer :: NX_BATHY,NY_BATHY
  integer :: ETOPO_INTERVAL

! filter final surface using box filter
  logical :: SMOOTH_THE_MODEL
  integer :: SIZE_FILTER_ONE_SIDE ! e.g. = 3 , 7

! use integer array to store values
  integer,dimension(:,:),allocatable :: ibathy_topo
  integer,dimension(:,:),allocatable :: ibathy_topo_ori

  integer :: ix,iy,minvalue,maxvalue
  integer :: ix_current,iy_current,ix_min,ix_max,iy_min,iy_max,ix_value,iy_value

  double precision :: value_sum,area_window

  integer :: file_size
  integer :: i,ier,idummy,icount
  logical :: file_exits
  character(len=256) :: arg(2)
  character(len=256) :: filename,file_ending

!----

  ! gets input parameters
  do i = 1, 2
    call get_command_argument(i,arg(i))
    if (i <= 2 .and. trim(arg(i)) == '') then
      print *, ' '
      print *, ' Usage: ./xsmooth_topo_bathy filename filter_size'
      print *, '    with'
      print *, '        filename      - input file, e.g. topo_bathy_etopo2v2c_original_unmodified_unsmoothed.dat'
      print *, '        filter_size   - file size (one side), e.g. 3'
      print *, ' '
      stop
    endif
  enddo
  filename = trim(arg(1))
  read(arg(2),*) SIZE_FILTER_ONE_SIDE

  print *,'smooth topo/bathy file'
  print *
  print *,'file name           : ',trim(filename)
  print *,'size filter one side: ',SIZE_FILTER_ONE_SIDE
  print *

  inquire(file=trim(filename),exist=file_exits)
  if (.not. file_exits) stop 'Error file not found'

  inquire(file=trim(filename),size=file_size,iostat=ier)
  if (ier /= 0 .or. file_size < 0) stop 'Error getting file size'

  print *,'file size = ',file_size / 1024. / 1024., '(MB)'

  open(unit=13,file=trim(filename),status='old',iostat=ier)
  if (ier /= 0) stop 'Error opening file'
  icount = 0
  do while(.true.)
    read(13,*,iostat=ier) idummy
    if (ier /= 0) exit
    icount = icount + 1
  enddo
  close(13)

  ! selects topo size
  select case(icount)
  case (NX_BATHY_1 * NY_BATHY_1)
    print *,'matches etopo1'
    NX_BATHY = NX_BATHY_1
    NY_BATHY = NY_BATHY_1
    ETOPO_INTERVAL = 1
  case (NX_BATHY_2 * NY_BATHY_2)
    print *,'matches etopo2'
    NX_BATHY = NX_BATHY_2
    NY_BATHY = NY_BATHY_2
    ETOPO_INTERVAL = 2
  case (NX_BATHY_4 * NY_BATHY_4)
    print *,'matches etopo4'
    NX_BATHY = NX_BATHY_4
    NY_BATHY = NY_BATHY_4
    ETOPO_INTERVAL = 4
  case default
    print *,'no match'
    print *,'invalid file size, does not match etopo1, etopo2 or etopo4 size. please check file...'
    stop 'Invalid file size'
  end select

! read the topography file
  print *,'ETOPO',ETOPO_INTERVAL
  print *,'NX_BATHY,NY_BATHY = ',NX_BATHY,NY_BATHY

  allocate(ibathy_topo(NX_BATHY,NY_BATHY),ibathy_topo_ori(NX_BATHY,NY_BATHY),stat=ier)
  if (ier /= 0) stop 'Error allocating ibathy_topo arrays'

  print *
  print *,'reading topo file'

  open(unit=13,file=trim(filename),status='old')
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

  ! checks if smoothing
  if (SIZE_FILTER_ONE_SIDE <= 0) then
    print *,'no smoothing'
    SMOOTH_THE_MODEL = .false.
  else
    print *,'using smoothing'
    SMOOTH_THE_MODEL = .true.
  endif

! smooth topography/bathymetry model
  if (SMOOTH_THE_MODEL) then
    print *
    print *,'smoothing topo file'
    if (SIZE_FILTER_ONE_SIDE < 1) stop 'SIZE_FILTER_ONE_SIDE must be greater than 1 for filter'
    print *,'size of window filter is ',2*SIZE_FILTER_ONE_SIDE+1,' x ',2*SIZE_FILTER_ONE_SIDE+1
    area_window = dble((2*SIZE_FILTER_ONE_SIDE+1)**2)

    do iy_current = 1,NY_BATHY

      if (mod(iy_current,100) == 0) print *,'smoothing line ',iy_current,' out of ',NY_BATHY

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
    ibathy_topo(:,:) = ibathy_topo_ori(:,:)

  endif

!----

! compute min and max after smoothing
  minvalue = minval(ibathy_topo)
  maxvalue = maxval(ibathy_topo)
  print *
  print *,'min and max of topography after smoothing = ',minvalue,maxvalue
  print *

! save the smoothed model
  if (SMOOTH_THE_MODEL) then
    print *,'saving the smoothed model'

    if (SIZE_FILTER_ONE_SIDE < 10) then
      write(file_ending,"('_smoothed_window_',i1,'.dat')") SIZE_FILTER_ONE_SIDE
    else
      write(file_ending,"('_smoothed_window_',i2,'.dat')") SIZE_FILTER_ONE_SIDE
    endif
    write(filename,"('topo_bathy_etopo',i1,a)") ETOPO_INTERVAL,trim(file_ending)

    open(unit=13,file=trim(filename),status='unknown')
    do iy=1,NY_BATHY
      do ix=1,NX_BATHY
        write(13,*) ibathy_topo(ix,iy)
      enddo
    enddo
    close(13)

    print *,'written to: ',trim(filename)
    print *,'can suppress white spaces in filtered model to save space if needed'
    print *
  endif

! create the PNM image
  print *,'creating PNM image'

! create image with grey levels
  do iy = 1,NY_BATHY
    do ix = 1,NX_BATHY
      ibathy_topo(ix,iy) = 255 * (ibathy_topo(ix,iy) - minvalue) / (maxvalue - minvalue)
      if (ibathy_topo(ix,iy) < 1) ibathy_topo(ix,iy) = 1
      if (ibathy_topo(ix,iy) > 255) ibathy_topo(ix,iy) = 255
    enddo
  enddo

! store image in PNM format with grey levels
  filename = 'image_topo_bathy.pnm'

! creating the header
  open(unit=27,file=trim(filename),status='unknown')
  write(27,100)
  write(27,102) NX_BATHY,NY_BATHY
  write(27,103)

 100 format('P3')
 102 format(i6,' ',i6)
 103 format('255')
 104 format(i3,1x,i3,1x,i3)

  do iy = 1,NY_BATHY
    do ix = 1,NX_BATHY
! write data value (red = green = blue to produce grey levels)
      write(27,104) ibathy_topo(ix,iy),ibathy_topo(ix,iy),ibathy_topo(ix,iy)
    enddo
  enddo

  close(27)
  print *, 'image written to: ',trim(filename)
  print *
  print *, 'done'

  end program smooth_topo_bathy_PPM_image

