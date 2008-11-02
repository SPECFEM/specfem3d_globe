!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
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

!
!---  create a movie of radial component of surface displacement
!---  in AVS or OpenDX format
!

  program create_movie_AVS_DX

  implicit none

  include "constants_topo.h"

! only use the four corners of each element to subsample the movie
! and therefore create much smaller movie files
  logical, parameter :: SUBSAMPLE_MOVIE = .true.

! remove elements that do not carry any wave
  logical, parameter :: REMOVE_ZERO_DATA = .true.

! add topography to the display of the wave field
! add a small offset to make sure wave field is displayed on top of the topography map
  logical, parameter :: ADD_TOPOGRAPHY_TO_DISPLAY = .true.
  real(kind=CUSTOM_REAL), parameter :: SMALL_OFFSET_TOPO = 7000. ! in meters

! amplify the radius of the sphere by a constant factor
  real(kind=CUSTOM_REAL), parameter :: AMPLIFY_RADIUS = 1.01_CUSTOM_REAL

! threshold in percent of the maximum below which we cut the amplitude
! and flag to cut amplitude below a certain threshold
  logical, parameter :: APPLY_THRESHOLD = .true.
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD = 1._CUSTOM_REAL / 100._CUSTOM_REAL

! flag to apply non linear scaling to normalized norm of displacement
! and coefficient of power law used for non linear scaling
  logical, parameter :: NONLINEAR_SCALING = .true.
  real(kind=CUSTOM_REAL), parameter :: POWER_SCALING = 0.30_CUSTOM_REAL

! filter final surface using box filter
  logical, parameter :: SMOOTH_THE_MODEL = .true.
  integer, parameter :: SIZE_FILTER_ONE_SIDE = 5

! parameters read from parameter file
  integer NEX_XI,NEX_ETA
  integer NSTEP,NTSTEP_BETWEEN_FRAMES,NCHUNKS
  integer NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  logical MOVIE_SURFACE

  integer i,j,it,it1,it2,iformat,ilimit,jlimit,iadd,jadd
  integer nspectot_AVS_max,nspectot_AVS_max_nonzero
  integer ispec
  integer ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4
  integer nframes,iframe

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,displn
  real(kind=CUSTOM_REAL) xcoord,ycoord,zcoord,rval,thetaval,phival,lat,long
  real(kind=CUSTOM_REAL) displx,disply,displz
  real(kind=CUSTOM_REAL) normal_x,normal_y,normal_z

  double precision :: phi,theta,r,elevation,latitude,longitude
  double precision min_field_current,max_field_current,max_absol

  logical USE_OPENDX,USE_GMT,USE_AVS,already_done

  character(len=150) outputname

  integer iproc,ipoin

! for sorting routine
  integer npointot,ilocnum,nglob,ielm,ieoff,ispecloc
  integer, dimension(:), allocatable :: iglob,loc,ireorder
  logical, dimension(:), allocatable :: ifseg,mask_point
  double precision, dimension(:), allocatable :: xp,yp,zp,xp_save,yp_save,zp_save,field_display

! for dynamic memory allocation
  integer ierror

! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz

  character(len=150) OUTPUT_FILES

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)
  integer ibathy_topo_ori(NX_BATHY,NY_BATHY)

  integer ix,iy,minvalue,maxvalue
  integer ix_current,iy_current,ix_min,ix_max,iy_min,iy_max,ix_value,iy_value

  double precision value_sum,area_window

! ************** PROGRAM STARTS HERE **************

  call read_AVS_DX_parameters(NEX_XI,NEX_ETA, &
           NSTEP,NTSTEP_BETWEEN_FRAMES, &
           NCHUNKS,MOVIE_SURFACE, &
           NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  if(.not. MOVIE_SURFACE) stop 'movie frames were not saved by the solver'

  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format'
  print *,'3 = create files in GMT xyz ASCII long/lat/Uz format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if(iformat<1 .or. iformat>3) stop 'exiting...'

  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  print *,'enter first time step of movie (e.g. 1)'
  read(5,*) it1

  print *,'enter last time step of movie (e.g. ',NSTEP,') (-1 for all frames)'
  read(5,*) it2
  if(it2 <= 0) it2 = NSTEP

  if(iformat == 1) then
    USE_OPENDX = .true.
    USE_AVS = .false.
    USE_GMT = .false.
  else if(iformat == 2) then
    USE_OPENDX = .false.
    USE_AVS = .true.
    USE_GMT = .false.
  else if(iformat == 3) then
    USE_OPENDX = .false.
    USE_AVS = .false.
    USE_GMT = .true.
  else
    stop 'error: invalid format'
  endif

!----

! add topography to the display of the wave field
  if(ADD_TOPOGRAPHY_TO_DISPLAY) then

! read the topography file
  print *,'NX_BATHY,NY_BATHY = ',NX_BATHY,NY_BATHY
  print *,'file used has a resolution of ',RESOLUTION_TOPO_FILE,' minutes'

  print *
  print *,'reading topo file'

  open(unit=13,file='DATA/topo_bathy/topo_bathy_etopo4_from_etopo2_subsampled.dat',status='old')
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

! compute min and max after smoothing
  minvalue = minval(ibathy_topo)
  maxvalue = maxval(ibathy_topo)
  print *,'min and max of topography after smoothing = ',minvalue,maxvalue

  endif

!----

  print *
  print *,'Generating all the movie frames, one after the other'
  print *

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *

  ilocnum = NGLLX * NGLLY * NEX_PER_PROC_XI * NEX_PER_PROC_ETA

  print *
  print *,'Allocating arrays of size ',ilocnum*NPROCTOT
  print *

  allocate(store_val_x(ilocnum,0:NPROCTOT-1),stat=ierror)
  if(ierror /= 0) stop 'error while allocating store_val_x'

  allocate(store_val_y(ilocnum,0:NPROCTOT-1),stat=ierror)
  if(ierror /= 0) stop 'error while allocating store_val_y'

  allocate(store_val_z(ilocnum,0:NPROCTOT-1),stat=ierror)
  if(ierror /= 0) stop 'error while allocating store_val_z'

  allocate(store_val_ux(ilocnum,0:NPROCTOT-1),stat=ierror)
  if(ierror /= 0) stop 'error while allocating store_val_ux'

  allocate(store_val_uy(ilocnum,0:NPROCTOT-1),stat=ierror)
  if(ierror /= 0) stop 'error while allocating store_val_uy'

  allocate(store_val_uz(ilocnum,0:NPROCTOT-1),stat=ierror)
  if(ierror /= 0) stop 'error while allocating store_val_uz'

  allocate(x(NGLLX,NGLLY),stat=ierror)
  if(ierror /= 0) stop 'error while allocating x'

  allocate(y(NGLLX,NGLLY),stat=ierror)
  if(ierror /= 0) stop 'error while allocating y'

  allocate(z(NGLLX,NGLLY),stat=ierror)
  if(ierror /= 0) stop 'error while allocating z'

  allocate(displn(NGLLX,NGLLY),stat=ierror)
  if(ierror /= 0) stop 'error while allocating displn'

  print *
  print *,'looping from ',it1,' to ',it2,' every ',NTSTEP_BETWEEN_FRAMES,' time steps'

! count number of movie frames
  nframes = 0
  do it = it1,it2
    if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) nframes = nframes + 1
  enddo
  print *
  print *,'total number of frames will be ',nframes
  if(nframes == 0) stop 'null number of frames'

! Make OpenDX think that each "grid cell" between GLL points is actually
! a finite element with four corners. This means that inside each real
! spectral element one should have (NGLL-1)^2 OpenDX "elements"

! define the total number of OpenDX "elements" at the surface
  if(SUBSAMPLE_MOVIE) then
    nspectot_AVS_max = NCHUNKS * NEX_XI * NEX_ETA
  else
    nspectot_AVS_max = NCHUNKS * NEX_XI * NEX_ETA * (NGLLX-1) * (NGLLY-1)
  endif

  print *
  print *,'there are a total of ',nspectot_AVS_max,' OpenDX "elements" at the surface'
  print *

! maximum theoretical number of points at the surface
  npointot = NGNOD2D_AVS_DX * nspectot_AVS_max

  print *
  print *,'Allocating arrays of size ',npointot
  print *

! allocate arrays for sorting routine
  allocate(iglob(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating iglob'

  allocate(loc(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating loc'

  allocate(ifseg(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating ifseg'

  allocate(xp(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating xp'

  allocate(yp(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating yp'

  allocate(zp(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating zp'

  allocate(xp_save(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating xp_save'

  allocate(yp_save(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating yp_save'

  allocate(zp_save(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating zp_save'

  allocate(field_display(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating field_display'

  allocate(mask_point(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating mask_point'

  allocate(ireorder(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating ireorder'

!--- ****** read data saved by solver ******

  print *

  if(APPLY_THRESHOLD) print *,'Will apply a threshold to amplitude below ',100.*THRESHOLD,' %'

  if(NONLINEAR_SCALING) print *,'Will apply a non linear scaling with coef ',POWER_SCALING

! --------------------------------------

  iframe = 0

  already_done = .false.

! loop on all the time steps in the range entered
  do it = it1,it2

! check if time step corresponds to a movie frame
  if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

  iframe = iframe + 1

  print *
  print *,'reading snapshot time step ',it,' out of ',NSTEP
  print *

! read all the elements from the same file
!! DK DK changed that for now  write(outputname,"('/moviedata',i6.6)") it
!! DK DK changed that for now  open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='old',action='read',form='unformatted')
  write(outputname,"('/scratch/komatits/moviedata',i6.6)") it
  open(unit=IOUT,file=outputname,status='old',action='read',form='unformatted')
  read(IOUT) store_val_x
  read(IOUT) store_val_y
  read(IOUT) store_val_z
  read(IOUT) store_val_ux
  read(IOUT) store_val_uy
  read(IOUT) store_val_uz
  close(IOUT)

! clear number of elements kept
  ispec = 0

! read points for all the slices
  do iproc = 0,NPROCTOT-1

! reset point number
    ipoin = 0

    do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA

      do j = 1,NGLLY
        do i = 1,NGLLX

          ipoin = ipoin + 1

          if(SUBSAMPLE_MOVIE .and. .not. &
             ((i == 1 .and. j == 1) .or. &
              (i == 1 .and. j == NGLLY) .or. &
              (i == NGLLX .and. j == 1) .or. &
              (i == NGLLX .and. j == NGLLY))) cycle

          xcoord = store_val_x(ipoin,iproc)
          ycoord = store_val_y(ipoin,iproc)
          zcoord = store_val_z(ipoin,iproc)

          displx = store_val_ux(ipoin,iproc)
          disply = store_val_uy(ipoin,iproc)
          displz = store_val_uz(ipoin,iproc)

! coordinates actually contain r theta phi, therefore convert back to x y z
          rval = xcoord
          thetaval = ycoord
          phival = zcoord
          call rthetaphi_2_xyz(xcoord,ycoord,zcoord,rval,thetaval,phival)

! compute unit normal vector to the surface
          normal_x = xcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
          normal_y = ycoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
          normal_z = zcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)

! save the results for this element
          x(i,j) = xcoord
          y(i,j) = ycoord
          z(i,j) = zcoord
          displn(i,j) = displx*normal_x + disply*normal_y + displz*normal_z

        enddo
      enddo

! assign the values of the corners of the OpenDX "elements"
      ispec = ispec + 1

      if(SUBSAMPLE_MOVIE) then
        ielm = ispec-1
        ilimit = 1
        jlimit = 1
        iadd = NGLLX-1
        jadd = NGLLY-1
      else
        ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
        ilimit = NGLLX-1
        jlimit = NGLLY-1
        iadd = 1
        jadd = 1
      endif

      do j = 1,jlimit
        do i = 1,ilimit
          ieoff = NGNOD2D_AVS_DX*(ielm + (i-1) + (j-1)*(NGLLX-1))
          do ilocnum = 1,NGNOD2D_AVS_DX
            if(ilocnum == 1) then
              xp(ieoff+ilocnum) = dble(x(i,j))
              yp(ieoff+ilocnum) = dble(y(i,j))
              zp(ieoff+ilocnum) = dble(z(i,j))
              field_display(ieoff+ilocnum) = dble(displn(i,j))
            elseif(ilocnum == 2) then
              xp(ieoff+ilocnum) = dble(x(i+iadd,j))
              yp(ieoff+ilocnum) = dble(y(i+iadd,j))
              zp(ieoff+ilocnum) = dble(z(i+iadd,j))
              field_display(ieoff+ilocnum) = dble(displn(i+iadd,j))
            elseif(ilocnum == 3) then
              xp(ieoff+ilocnum) = dble(x(i+iadd,j+jadd))
              yp(ieoff+ilocnum) = dble(y(i+iadd,j+jadd))
              zp(ieoff+ilocnum) = dble(z(i+iadd,j+jadd))
              field_display(ieoff+ilocnum) = dble(displn(i+iadd,j+jadd))
            else
              xp(ieoff+ilocnum) = dble(x(i,j+jadd))
              yp(ieoff+ilocnum) = dble(y(i,j+jadd))
              zp(ieoff+ilocnum) = dble(z(i,j+jadd))
              field_display(ieoff+ilocnum) = dble(displn(i,j+jadd))
            endif
          enddo
        enddo
      enddo

    enddo

  enddo

! compute min and max of data value to normalize
  min_field_current = minval(field_display(:))
  max_field_current = maxval(field_display(:))

! make sure range is always symmetric and center is in zero
! this assumption works only for fields that can be negative
! would not work for norm of vector for instance
! (we would lose half of the color palette if no negative values)
  max_absol = max(abs(min_field_current),abs(max_field_current))
  min_field_current = - max_absol
  max_field_current = + max_absol

! print minimum and maximum amplitude in current snapshot
  print *
  print *,'minimum amplitude in current snapshot = ',min_field_current
  print *,'maximum amplitude in current snapshot = ',max_field_current
  print *

! normalize field to [0:1]
  field_display(:) = (field_display(:) - min_field_current) / (max_field_current - min_field_current)

! rescale to [-1,1]
  field_display(:) = 2.*field_display(:) - 1.

! apply threshold to normalized field
  if(APPLY_THRESHOLD) where(abs(field_display(:)) <= THRESHOLD) field_display = 0.

! apply non linear scaling to normalized field if needed
  if(NONLINEAR_SCALING) then
    where(field_display(:) >= 0.)
      field_display = field_display ** POWER_SCALING
    elsewhere
      field_display = - abs(field_display) ** POWER_SCALING
    endwhere
  endif

! copy coordinate arrays since the sorting routine does not preserve them
! do this only once for each movie because all movie frames use the same mesh
  if(.not. already_done) then
    already_done = .true.

! add topography to the display of the wave field
    if(ADD_TOPOGRAPHY_TO_DISPLAY) then

    do ipoin = 1,npointot

      call xyz_2_rthetaphi_dble(xp(ipoin),yp(ipoin),zp(ipoin),r,theta,phi)
!! DK DK in principle this should take ellipticity into account
      longitude = phi * 180.d0 / PI
      latitude = 90.d0 - theta * 180.d0 / PI

! use an Earth model of mean radius equal to 1, and add elevation to each point, amplified in order to see it
      call get_topo_bathy_local_copy(latitude,longitude,elevation,ibathy_topo)
      r = (R_EARTH + EXAGGERATION_FACTOR_TOPO*elevation + SMALL_OFFSET_TOPO) / R_EARTH

      xp_save(ipoin) = r*sin(theta)*cos(phi)
      yp_save(ipoin) = r*sin(theta)*sin(phi)
      zp_save(ipoin) = r*cos(theta)
    enddo

    else
      xp_save(:) = xp(:)
      yp_save(:) = yp(:)
      zp_save(:) = zp(:)
    endif

!--- sort the list based upon coordinates to get rid of multiples
    print *,'sorting list of points'
    call get_global_AVS(nspectot_AVS_max,xp,yp,zp,iglob,loc,ifseg,nglob,npointot)
  endif

!--- print total number of points found
  print *
  print *,'found a total of ',nglob,' points'
  print *,'initial number of points (with multiples) was ',npointot

!--- ****** create AVS file using sorted list ******

! create file name and open file
  if(USE_OPENDX) then
    write(outputname,"('/DX_movie_',i6.6,'.dx')") it
    open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
  else if(USE_AVS) then
    write(outputname,"('/AVS_movie_',i6.6,'.inp')") it
    open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
  else if(USE_GMT) then
    write(outputname,"('/gmt_movie_',i6.6,'.xyz')") it
    open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
  else
    stop 'wrong output format selected'
  endif

  if(USE_GMT) then

! output list of points
    mask_point = .false.
    do ispec = 1,nspectot_AVS_max
      ieoff = NGNOD2D_AVS_DX*(ispec-1)
! four points for each element
      do ilocnum = 1,NGNOD2D_AVS_DX
        ibool_number = iglob(ilocnum+ieoff)
        if(.not. mask_point(ibool_number)) then
          xcoord = sngl(xp_save(ilocnum+ieoff))
          ycoord = sngl(yp_save(ilocnum+ieoff))
          zcoord = sngl(zp_save(ilocnum+ieoff))
          call xyz_2_rthetaphi(xcoord,ycoord,zcoord,rval,thetaval,phival)
          lat = (PI/2.0-thetaval)*180.0/PI
          long = phival*180.0/PI
          if(long > 180.0) long = long-360.0
          write(11,*) long,lat,sngl(field_display(ilocnum+ieoff))
        endif
        mask_point(ibool_number) = .true.
      enddo
    enddo

  else

! output list of points
  mask_point = .false.
  ipoin = 0
  do ispec=1,nspectot_AVS_max
    ieoff = NGNOD2D_AVS_DX*(ispec-1)
! four points for each element
    do ilocnum = 1,NGNOD2D_AVS_DX
      ibool_number = iglob(ilocnum+ieoff)
      if(.not. mask_point(ibool_number)) then
        ipoin = ipoin + 1
        ireorder(ibool_number) = ipoin
        if(USE_OPENDX) then
          write(11,"(f10.7,1x,f10.7,1x,f10.7)") &
            xp_save(ilocnum+ieoff),yp_save(ilocnum+ieoff),zp_save(ilocnum+ieoff)
        else if(USE_AVS) then
          write(11,"(i6,1x,f10.7,1x,f10.7,1x,f10.7)") ireorder(ibool_number), &
            xp_save(ilocnum+ieoff),yp_save(ilocnum+ieoff),zp_save(ilocnum+ieoff)
        endif
      endif
      mask_point(ibool_number) = .true.
    enddo
  enddo

  if(USE_OPENDX) then
    if(REMOVE_ZERO_DATA) then
      nspectot_AVS_max_nonzero = 0
      do ispec=1,nspectot_AVS_max
        ieoff = NGNOD2D_AVS_DX*(ispec-1)
        if(REMOVE_ZERO_DATA .and. &
            field_display(ieoff + 1) == 0. .and. &
            field_display(ieoff + 2) == 0. .and. &
            field_display(ieoff + 3) == 0. .and. &
            field_display(ieoff + 4) == 0.) cycle
        nspectot_AVS_max_nonzero = nspectot_AVS_max_nonzero + 1
      enddo
    else
      nspectot_AVS_max_nonzero = nspectot_AVS_max
    endif
    write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspectot_AVS_max_nonzero,' data follows'
  endif

! output list of elements
  do ispec=1,nspectot_AVS_max
    ieoff = NGNOD2D_AVS_DX*(ispec-1)

    if(REMOVE_ZERO_DATA .and. &
        field_display(ieoff + 1) == 0. .and. &
        field_display(ieoff + 2) == 0. .and. &
        field_display(ieoff + 3) == 0. .and. &
        field_display(ieoff + 4) == 0.) cycle

! four points for each element
    ibool_number1 = iglob(ieoff + 1)
    ibool_number2 = iglob(ieoff + 2)
    ibool_number3 = iglob(ieoff + 3)
    ibool_number4 = iglob(ieoff + 4)
    if(USE_OPENDX) then
! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
      write(11,"(i10,1x,i10,1x,i10,1x,i10)") ireorder(ibool_number1)-1, &
        ireorder(ibool_number4)-1,ireorder(ibool_number2)-1,ireorder(ibool_number3)-1
    else
      write(11,"(i10,' 1 quad ',i10,1x,i10,1x,i10,1x,i10)") ispec,ireorder(ibool_number1), &
        ireorder(ibool_number2),ireorder(ibool_number3),ireorder(ibool_number4)
    endif
  enddo

  if(USE_OPENDX) then
    write(11,*) 'attribute "element type" string "quads"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',nglob,' data follows'
  else
! dummy text for labels
    write(11,*) '1 1'
    write(11,*) 'a, b'
  endif

! output data values
  mask_point = .false.

! output point data
  do ispec=1,nspectot_AVS_max
  ieoff = NGNOD2D_AVS_DX*(ispec-1)
! four points for each element
  do ilocnum = 1,NGNOD2D_AVS_DX
    ibool_number = iglob(ilocnum+ieoff)
    if(.not. mask_point(ibool_number)) then
      if(USE_OPENDX) then
        write(11,"(f7.2)") field_display(ilocnum+ieoff)
      else
        write(11,"(i10,1x,f7.2)") ireorder(ibool_number),field_display(ilocnum+ieoff)
      endif
    endif
    mask_point(ibool_number) = .true.
  enddo
  enddo

! define OpenDX field
  if(USE_OPENDX) then
    write(11,*) 'attribute "dep" string "positions"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'
  endif

! end of test for GMT format
  endif

  close(11)

! end of loop and test on all the time steps for all the movie images
  endif

  enddo

  print *
  print *,'done creating movie'
  print *

  if(USE_OPENDX) print *,'DX files are stored in ', trim(OUTPUT_FILES), '/DX_*.dx'
  if(USE_AVS) print *,'AVS files are stored in ', trim(OUTPUT_FILES), '/AVS_*.inp'
  if(USE_GMT) print *,'GMT files are stored in ', trim(OUTPUT_FILES), '/gmt_*.xyz'
  print *

  end program create_movie_AVS_DX

!
!=====================================================================
!

  subroutine read_AVS_DX_parameters(NEX_XI,NEX_ETA, &
          NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NCHUNKS,MOVIE_SURFACE, &
          NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  implicit none

  include "constants_topo.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,IT_LAST_VALUE_DUMPED,INTERVAL_DUMP_FILES,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
          ifirst_layer_aniso,ilast_layer_aniso

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_layer_has_a_doubling

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               NGLOB

  character(len=150) MODEL
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  print *
  print *,'reading parameter file'
  print *

! read the parameter file and compute additional parameters
    call read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,IT_LAST_VALUE_DUMPED,INTERVAL_DUMP_FILES,NCHUNKS,DT, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
          NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
          NSPEC, &
          NSPEC2D_XI, &
          NSPEC2D_ETA, &
          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
          NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
          NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
          ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_layer_has_a_doubling,rmins,rmaxs,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
          DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
          USE_BINARY_FOR_LARGE_FILE,ifirst_layer_aniso,ilast_layer_aniso,.false.)

  if(MOVIE_COARSE) stop 'create_movie_AVS_DX does not work with MOVIE_COARSE'

  end subroutine read_AVS_DX_parameters

! ------------------------------------------------------------------

  subroutine get_global_AVS(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  include "constants_topo.h"

  integer npointot
  integer iglob(npointot),loc(npointot)
  logical ifseg(npointot)
  double precision xp(npointot),yp(npointot),zp(npointot)
  integer nspec,nglob

  integer ispec,i,j
  integer ieoff,ilocnum,nseg,ioff,iseg,ig

! for dynamic memory allocation
  integer ierror

  integer, dimension(:), allocatable :: ind,ninseg,iwork
  double precision, dimension(:), allocatable :: work

  print *
  print *,'Allocating arrays of size ',npointot
  print *

! dynamically allocate arrays
  allocate(ind(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating ind'

  allocate(ninseg(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating ninseg'

  allocate(iwork(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating iwork'

  allocate(work(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating work'

! establish initial pointers
  do ispec=1,nspec
    ieoff=NGNOD2D_AVS_DX*(ispec-1)
    do ilocnum=1,NGNOD2D_AVS_DX
      loc(ieoff+ilocnum)=ieoff+ilocnum
    enddo
  enddo

  ifseg(:)=.false.

  nseg=1
  ifseg(1)=.true.
  ninseg(1)=npointot

  do j=1,NDIM

! sort within each segment
  ioff=1
  do iseg=1,nseg
    if(j == 1) then
      call rank(xp(ioff),ind,ninseg(iseg))
    else if(j == 2) then
      call rank(yp(ioff),ind,ninseg(iseg))
    else
      call rank(zp(ioff),ind,ninseg(iseg))
    endif
    call swap_all(loc(ioff),xp(ioff),yp(ioff),zp(ioff),iwork,work,ind,ninseg(iseg))
    ioff=ioff+ninseg(iseg)
  enddo

! check for jumps in current coordinate
! compare the coordinates of the points within a small tolerance
  if(j == 1) then
    do i=2,npointot
      if(dabs(xp(i)-xp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  else if(j == 2) then
    do i=2,npointot
      if(dabs(yp(i)-yp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  else
    do i=2,npointot
      if(dabs(zp(i)-zp(i-1)) > SMALLVALTOL) ifseg(i)=.true.
    enddo
  endif

! count up number of different segments
  nseg=0
  do i=1,npointot
    if(ifseg(i)) then
      nseg=nseg+1
      ninseg(nseg)=1
    else
      ninseg(nseg)=ninseg(nseg)+1
    endif
  enddo
  enddo

! assign global node numbers (now sorted lexicographically)
  ig=0
  do i=1,npointot
    if(ifseg(i)) ig=ig+1
    iglob(loc(i))=ig
  enddo

  nglob=ig

! deallocate arrays
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iwork)
  deallocate(work)

! -----------------------------------

! get_global_AVS internal procedures follow

! sorting routines put in same file to allow for inlining

  contains

! -----------------------------------

  subroutine rank(A,IND,N)
!
! Use Heap Sort (Numerical Recipes)
!
  implicit none

  integer n
  double precision A(n)
  integer IND(n)

  integer i,j,l,ir,indx
  double precision q

  do j=1,n
   IND(j)=j
  enddo

  if (n == 1) return

  L=n/2+1
  ir=n
  100 CONTINUE
   IF (l>1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   ENDIF
   i=l
   j=l+l
  200    CONTINUE
   IF (J <= IR) THEN
      IF (J<IR) THEN
         IF ( A(IND(j))<A(IND(j+1)) ) j=j+1
      ENDIF
      IF (q<A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      ENDIF
   goto 200
   ENDIF
   IND(I)=INDX
  goto 100
  end subroutine rank

! ------------------------------------------------------------------

  subroutine swap_all(IA,A,B,C,IW,W,ind,n)
!
! swap arrays IA, A, B and C according to addressing in array IND
!
  implicit none

  integer n

  integer IND(n)
  integer IA(n),IW(n)
  double precision A(n),B(n),C(n),W(n)

  integer i

  IW(:) = IA(:)
  W(:) = A(:)

  do i=1,n
    IA(i)=IW(ind(i))
    A(i)=W(ind(i))
  enddo

  W(:) = B(:)

  do i=1,n
    B(i)=W(ind(i))
  enddo

  W(:) = C(:)

  do i=1,n
    C(i)=W(ind(i))
  enddo

  end subroutine swap_all

! ------------------------------------------------------------------

  end subroutine get_global_AVS

!----------------------------

  subroutine get_topo_bathy_local_copy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

!! indentical code to that on the SVN server but with NX_BATHY set to real value
!! instead of a fictitious value of 1 to reduce memory size and turn topography off

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

  end subroutine get_topo_bathy_local_copy

