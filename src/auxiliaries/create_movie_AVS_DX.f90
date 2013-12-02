!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  program xcreate_movie_AVS_DX

  use constants
  use shared_parameters,only: &
    NEX_XI,NEX_ETA,NSTEP,NTSTEP_BETWEEN_FRAMES, &
    NCHUNKS,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
    OUTPUT_FILES,MOVIE_COARSE

  implicit none

!---------------------
! USER PARAMETER

! threshold and normalization parameters

! single parameter to turn off all thresholding
  logical, parameter :: ALL_THRESHOLD_OFF = .true.

! threshold in percent of the maximum below which we cut the amplitude
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD = 1._CUSTOM_REAL / 100._CUSTOM_REAL

! flag to apply non linear scaling to normalized norm of displacement
  logical, parameter :: NONLINEAR_SCALING = .false.

! uses fixed max_value to normalize instead of max of current wavefield
  logical, parameter :: FIX_SCALING = .false.
  real,parameter:: MAX_VALUE = 1.0

! coefficient of power law used for non linear scaling
  real(kind=CUSTOM_REAL), parameter :: POWER_SCALING = 0.30_CUSTOM_REAL

! flag to cut amplitude below a certain threshold
  logical, parameter :: APPLY_THRESHOLD = .false.

!---------------------

  ! user input
  integer :: it1,it2
  integer :: iformat
  integer :: USE_COMPONENT

  integer :: i,j,it

  integer :: nspectot_AVS_max
  integer :: ispec
  integer :: ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,displn
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord,rval,thetaval,phival,lat,long
  real(kind=CUSTOM_REAL) :: displx,disply,displz
  real(kind=CUSTOM_REAL) :: normal_x,normal_y,normal_z
  real(kind=CUSTOM_REAL) :: RRval,rhoval
  real(kind=CUSTOM_REAL) :: thetahat_x,thetahat_y,thetahat_z
  real(kind=CUSTOM_REAL) :: phihat_x,phihat_y

  double precision :: min_field_current,max_field_current,max_absol

  logical :: USE_OPENDX,UNIQUE_FILE,USE_GMT,USE_AVS

  integer :: nframes,iframe

  character(len=150) :: outputname

  integer :: iproc,ipoin

! for sorting routine
  integer :: npointot,ilocnum,nglob,ielm,ieoff,ispecloc
  integer :: NIT
  integer, dimension(:), allocatable :: iglob,loc,ireorder
  logical, dimension(:), allocatable :: ifseg,mask_point
  double precision, dimension(:), allocatable :: xp,yp,zp,xp_save,yp_save,zp_save,field_display

! for dynamic memory allocation
  integer :: ierror

! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz


! ************** PROGRAM STARTS HERE **************

  call read_AVS_DX_parameters()

  ! user input
  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format with individual files'
  print *,'3 = create files in AVS UCD format with one time-dependent file'
  print *,'4 = create files in GMT xyz Ascii long/lat/U format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if( iformat < 1 .or. iformat > 4 ) stop 'exiting...'

  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  print *,'enter first time step of movie (e.g. 1)'
  read(5,*) it1

  print *,'enter last time step of movie (e.g. ',NSTEP,'or -1 for all)'
  read(5,*) it2

  print *,'enter component (e.g. 1=Z, 2=N, 3=E)'
  read(5,*) USE_COMPONENT
  if( USE_COMPONENT < 1 .or. USE_COMPONENT > 3 ) stop 'component must be 1, 2 or 3'

  print *
  print *,'--------'
  print *

  ! checks options
  if( it2 == -1 ) it2 = NSTEP

! --------------------------------------

  ! runs the main program
  ! sets flags
  if(iformat == 1) then
    USE_OPENDX = .true.
    USE_AVS = .false.
    USE_GMT = .false.
    UNIQUE_FILE = .false.
  else if(iformat == 2) then
    USE_OPENDX = .false.
    USE_AVS = .true.
    USE_GMT = .false.
    UNIQUE_FILE = .false.
  else if(iformat == 3) then
    USE_OPENDX = .false.
    USE_AVS = .true.
    USE_GMT = .false.
    UNIQUE_FILE = .true.
  else if(iformat == 4) then
    USE_OPENDX = .false.
    USE_AVS = .false.
    USE_GMT = .true.
    UNIQUE_FILE = .false.
  else
    stop 'error: invalid format'
  endif

  print *
  print *,'Recombining all movie frames to create a movie'
  print *
  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *

  if(MOVIE_COARSE) then
    ! note:
    ! nex_per_proc_xi*nex_per_proc_eta = nex_xi*nex_eta/nproc = nspec2d_top(iregion_crust_mantle)
    ! used in specfem3D.f90
    ! and ilocnum = nmovie_points = 2 * 2 * NEX_XI * NEX_ETA / NPROC
    ilocnum = 2 * 2 * NEX_PER_PROC_XI*NEX_PER_PROC_ETA
    NIT = NGLLX-1
  else
    ilocnum = NGLLX*NGLLY*NEX_PER_PROC_XI*NEX_PER_PROC_ETA
    NIT = 1
  endif

  print *
  print *,'Allocating arrays for reading data of size ',ilocnum*NPROCTOT,'=', &
                            6*ilocnum*NPROCTOT*CUSTOM_REAL/1024./1024.,'MB'
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
  if(MOVIE_COARSE) then
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

  ! movie point locations
  outputname = "/moviedata_xyz.bin"
  open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname), &
       status='old',action='read',form='unformatted',iostat=ierror)
  if(ierror /= 0 ) then
    print*,'error opening file: ',trim(OUTPUT_FILES)//trim(outputname)
    stop 'error opening moviedata file'
  endif

  ! reads in point locations
  ! (given as r theta phi for geocentric coordinate system)
  read(IOUT) store_val_x
  read(IOUT) store_val_y
  read(IOUT) store_val_z
  close(IOUT)

! --------------------------------------

  iframe = 0

  ! loop on all the time steps in the range entered
  do it = it1,it2

    ! check if time step corresponds to a movie frame
    if(mod(it,NTSTEP_BETWEEN_FRAMES) /= 0) cycle

    iframe = iframe + 1

    print *
    print *,'reading snapshot time step ',it,' out of ',NSTEP
    print *

    ! read all the elements from the same file
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname, &
         status='old',action='read',form='unformatted',iostat=ierror)
    if(ierror /= 0 ) then
      print*,'error opening file: ',trim(OUTPUT_FILES)//outputname
      stop 'error opening moviedata file'
    endif

    ! reads in associated values (displacement or velocity..)
    read(IOUT) store_val_ux
    read(IOUT) store_val_uy
    read(IOUT) store_val_uz

    close(IOUT)

    print *,'  file ',trim(OUTPUT_FILES)//trim(outputname)
    print *
    !debug
    print *,"  data ux min/max: ",minval(store_val_ux),maxval(store_val_ux)
    print *,"  data uy min/max: ",minval(store_val_uy),maxval(store_val_uy)
    print *,"  data uz min/max: ",minval(store_val_uz),maxval(store_val_uz)
    print *

    ! clear number of elements kept
    ispec = 0

    ! read points for all the slices
    do iproc = 0,NPROCTOT-1

      ! reset point number
      ipoin = 0

      do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA

        ! gets element values
        do j = 1,NGLLY,NIT
          do i = 1,NGLLX,NIT
            ipoin = ipoin + 1

            ! coordinates actually contain r theta phi
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

            ! save the results for this element
            x(i,j) = xcoord
            y(i,j) = ycoord
            z(i,j) = zcoord

            ! saves the desired component
            if(USE_COMPONENT == 1) then
               ! compute unit normal vector to the surface
               RRval = sqrt(xcoord**2 + ycoord**2 + zcoord**2)
               if( RRval < 1.e-10 ) stop 'error unit normal vector'
               normal_x = xcoord / RRval
               normal_y = ycoord / RRval
               normal_z = zcoord / RRval

               displn(i,j) = displx*normal_x + disply*normal_y + displz*normal_z

            else if(USE_COMPONENT == 2) then

               ! compute unit tangent vector to the surface (N-S)
               RRval = sqrt(xcoord**2 + ycoord**2 + zcoord**2)
               if( RRval < 1.e-10 ) stop 'error unit normal vector'
               rhoval = sqrt(xcoord**2 + ycoord**2)
               if( rhoval < 1.e-10 ) then
                ! location at pole
                thetahat_x = 0.0
                thetahat_y = 0.0
               else
                thetahat_x = (zcoord*xcoord) / (rhoval*RRval)
                thetahat_y = (zcoord*ycoord) / (rhoval*RRval)
               endif
               thetahat_z = - rhoval/RRval

               displn(i,j) = - (displx*thetahat_x + disply*thetahat_y + displz*thetahat_z)

            else if(USE_COMPONENT == 3) then

               ! compute unit tangent to the surface (E-W)
               rhoval = sqrt(xcoord**2 + ycoord**2)
               if( rhoval < 1.e-10 ) then
                ! location at pole
                phihat_x = 0.0
                phihat_y = 0.0
               else
                phihat_x = -ycoord / rhoval
                phihat_y = xcoord / rhoval
               endif

               displn(i,j) = displx*phihat_x + disply*phihat_y
            endif

          enddo
        enddo

        ! assign the values of the corners of the OpenDX "elements"
        ispec = ispec + 1
        if(MOVIE_COARSE) then
          ielm = ispec-1
        else
          ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
        endif

        do j = 1,NGLLY-NIT
          do i = 1,NGLLX-NIT
            ! offset (starts at 0)
            if(MOVIE_COARSE) then
              ieoff = NGNOD2D_AVS_DX * ielm
            else
              ieoff = NGNOD2D_AVS_DX * (ielm+(i-1)+(j-1)*(NGLLX-1))
            endif

            ! stores 4 points as element corners
            do ilocnum = 1,NGNOD2D_AVS_DX
              ! for movie_coarse e.g. x(i,j) is defined at x(1,1), x(1,NGLLY), x(NGLLX,1) and x(NGLLX,NGLLY)
              ! be aware that for the cubed sphere, the mapping changes for different chunks,
              ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates
              if(MOVIE_COARSE) then
                if(ilocnum == 1) then
                  xp(ieoff+ilocnum) = dble(x(1,1))
                  yp(ieoff+ilocnum) = dble(y(1,1))
                  zp(ieoff+ilocnum) = dble(z(1,1))
                  field_display(ieoff+ilocnum) = dble(displn(1,1))
                else if(ilocnum == 2) then
                  xp(ieoff+ilocnum) = dble(x(NGLLX,1))
                  yp(ieoff+ilocnum) = dble(y(NGLLX,1))
                  zp(ieoff+ilocnum) = dble(z(NGLLX,1))
                  field_display(ieoff+ilocnum) = dble(displn(NGLLX,1))
                else if(ilocnum == 3) then
                  xp(ieoff+ilocnum) = dble(x(NGLLX,NGLLY))
                  yp(ieoff+ilocnum) = dble(y(NGLLX,NGLLY))
                  zp(ieoff+ilocnum) = dble(z(NGLLX,NGLLY))
                  field_display(ieoff+ilocnum) = dble(displn(NGLLX,NGLLY))
                else
                  xp(ieoff+ilocnum) = dble(x(1,NGLLY))
                  yp(ieoff+ilocnum) = dble(y(1,NGLLY))
                  zp(ieoff+ilocnum) = dble(z(1,NGLLY))
                  field_display(ieoff+ilocnum) = dble(displn(1,NGLLY))
                endif
              else
                ! movie saved on fine grid (all gll points)
                if(ilocnum == 1) then
                  xp(ieoff+ilocnum) = dble(x(i,j))
                  yp(ieoff+ilocnum) = dble(y(i,j))
                  zp(ieoff+ilocnum) = dble(z(i,j))
                  field_display(ieoff+ilocnum) = dble(displn(i,j))
                else if(ilocnum == 2) then
                  xp(ieoff+ilocnum) = dble(x(i+1,j))
                  yp(ieoff+ilocnum) = dble(y(i+1,j))
                  zp(ieoff+ilocnum) = dble(z(i+1,j))
                  field_display(ieoff+ilocnum) = dble(displn(i+1,j))
                else if(ilocnum == 3) then
                  xp(ieoff+ilocnum) = dble(x(i+1,j+1))
                  yp(ieoff+ilocnum) = dble(y(i+1,j+1))
                  zp(ieoff+ilocnum) = dble(z(i+1,j+1))
                  field_display(ieoff+ilocnum) = dble(displn(i+1,j+1))
                else
                  xp(ieoff+ilocnum) = dble(x(i,j+1))
                  yp(ieoff+ilocnum) = dble(y(i,j+1))
                  zp(ieoff+ilocnum) = dble(z(i,j+1))
                  field_display(ieoff+ilocnum) = dble(displn(i,j+1))
                endif
              endif ! MOVIE_COARSE
            enddo ! ilocnum

          enddo
        enddo

      enddo ! ispecloc

    enddo

    !---------------------------------
    ! apply normalization and thresholding, if desired

    if (.not. ALL_THRESHOLD_OFF) then

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
      if( FIX_SCALING ) then
        print *,'  to be normalized by : ',MAX_VALUE
        if( max_field_current > MAX_VALUE ) stop 'increase MAX_VALUE'
      endif
      print *

      ! normalize field to [0:1]
      print *,'normalizing... '
      if( FIX_SCALING ) then
        field_display(:) = (field_display(:) + MAX_VALUE) / (2.0*MAX_VALUE)
      else
        field_display(:) = (field_display(:) - min_field_current) / (max_field_current - min_field_current)
      endif
      ! rescale to [-1,1]
      field_display(:) = 2.*field_display(:) - 1.

      ! apply threshold to normalized field
      if(APPLY_THRESHOLD) then
        print *,'thresholding... '
        where(abs(field_display(:)) <= THRESHOLD) field_display = 0.
      endif

      ! apply non linear scaling to normalized field if needed
      if(NONLINEAR_SCALING) then
        print *,'nonlinear scaling... '
        where(field_display(:) >= 0.)
          field_display = field_display ** POWER_SCALING
        elsewhere
          field_display = - abs(field_display) ** POWER_SCALING
        endwhere
      endif

      print *,'color scaling... '
      ! map back to [0,1]
      field_display(:) = (field_display(:) + 1.) / 2.

      ! map field to [0:255] for AVS color scale
      field_display(:) = 255. * field_display(:)

    else

      ! compute min and max of data value
      min_field_current = minval(field_display(:))
      max_field_current = maxval(field_display(:))

      ! print minimum and maximum amplitude in current snapshot
      print *
      print *,'minimum amplitude in current snapshot = ',min_field_current
      print *,'maximum amplitude in current snapshot = ',max_field_current
      print *

    endif

    !---------------------------------

    ! copy coordinate arrays since the sorting routine does not preserve them
    print *,'sorting... '
    xp_save(:) = xp(:)
    yp_save(:) = yp(:)
    zp_save(:) = zp(:)

    !--- sort the list based upon coordinates to get rid of multiples
    print *,'sorting list of points'
    call get_global_AVS(nspectot_AVS_max,xp,yp,zp,iglob,loc,ifseg,nglob,npointot)

    !--- print total number of points found
    print *
    print *,'found a total of ',nglob,' points'
    print *,'initial number of points (with multiples) was ',npointot

    !--- ****** create AVS file using sorted list ******

    ! create file name and open file
    if(USE_OPENDX) then
      if( USE_COMPONENT == 1) then
        write(outputname,"('/DX_movie_',i6.6,'.Z.dx')") it
      else if( USE_COMPONENT == 2) then
        write(outputname,"('/DX_movie_',i6.6,'.N.dx')") it
      else if( USE_COMPONENT == 3) then
        write(outputname,"('/DX_movie_',i6.6,'.E.dx')") it
      endif
      open(unit=11,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown')
      write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
    else if(USE_AVS) then
      if(UNIQUE_FILE .and. iframe == 1) then
        if( USE_COMPONENT == 1) then
          outputname = '/AVS_movie_all.Z.inp'
        else if( USE_COMPONENT == 2) then
          outputname = '/AVS_movie_all.N.inp'
        else if( USE_COMPONENT == 3) then
          outputname = '/AVS_movie_all.E.inp'
        endif
        open(unit=11,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown')
        write(11,*) nframes
        write(11,*) 'data'
        write(11,"('step',i1,' image',i1)") 1,1
        write(11,*) nglob,' ',nspectot_AVS_max
      else if(.not. UNIQUE_FILE) then
        if( USE_COMPONENT == 1) then
          write(outputname,"('/AVS_movie_',i6.6,'.Z.inp')") it
        else if( USE_COMPONENT == 2) then
          write(outputname,"('/AVS_movie_',i6.6,'.N.inp')") it
        else if( USE_COMPONENT == 3) then
          write(outputname,"('/AVS_movie_',i6.6,'.E.inp')") it
        endif
        open(unit=11,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown')
        write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
      endif
    else if(USE_GMT) then
      if( USE_COMPONENT == 1) then
        write(outputname,"('/gmt_movie_',i6.6,'.Z.xyz')") it
      else if( USE_COMPONENT == 2) then
        write(outputname,"('/gmt_movie_',i6.6,'.N.xyz')") it
      else if( USE_COMPONENT == 3) then
        write(outputname,"('/gmt_movie_',i6.6,'.E.xyz')") it
      endif
      open(unit=11,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown')
    else
      stop 'wrong output format selected'
    endif

    if(USE_GMT) then

      ! output list of points
      mask_point = .false.
      do ispec=1,nspectot_AVS_max
        ieoff = NGNOD2D_AVS_DX*(ispec-1)
        ! four points for each element
        do ilocnum = 1,NGNOD2D_AVS_DX
          ibool_number = iglob(ilocnum+ieoff)
          if(.not. mask_point(ibool_number)) then
            ! gets cartesian coordinates
            xcoord = sngl(xp_save(ilocnum+ieoff))
            ycoord = sngl(yp_save(ilocnum+ieoff))
            zcoord = sngl(zp_save(ilocnum+ieoff))

            ! location latitude/longitude (with geocentric colatitude theta )
            call xyz_2_rthetaphi(xcoord,ycoord,zcoord,rval,thetaval,phival)

            ! note: converts the geocentric colatitude to a geographic colatitude
            call geocentric_2_geographic_cr(thetaval,thetaval)

            ! gets geographic latitude and longitude in degrees
            lat = (PI_OVER_TWO-thetaval)*RADIANS_TO_DEGREES
            long = phival*RADIANS_TO_DEGREES
            if(long > 180.0) long = long-360.0

            write(11,*) long,lat,sngl(field_display(ilocnum+ieoff))
          endif
          mask_point(ibool_number) = .true.
        enddo
      enddo

    else
      ! if unique file, output geometry only once
      if(.not. UNIQUE_FILE .or. iframe == 1) then

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
                write(11,"(i10,1x,f10.7,1x,f10.7,1x,f10.7)") ireorder(ibool_number), &
                  xp_save(ilocnum+ieoff),yp_save(ilocnum+ieoff),zp_save(ilocnum+ieoff)
              endif
            endif
            mask_point(ibool_number) = .true.
          enddo
        enddo

        if(USE_OPENDX) &
          write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspectot_AVS_max,' data follows'

        ! output list of elements
        do ispec=1,nspectot_AVS_max
          ieoff = NGNOD2D_AVS_DX*(ispec-1)
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

      endif

      if(USE_OPENDX) then
        write(11,*) 'attribute "element type" string "quads"'
        write(11,*) 'attribute "ref" string "positions"'
        write(11,*) 'object 3 class array type float rank 0 items ',nglob,' data follows'
      else
        if(UNIQUE_FILE) then
          ! step number for AVS multistep file
          if(iframe > 1) then
            if(iframe < 10) then
              write(11,"('step',i1,' image',i1)") iframe,iframe
            else if(iframe < 100) then
              write(11,"('step',i2,' image',i2)") iframe,iframe
            else if(iframe < 1000) then
              write(11,"('step',i3,' image',i3)") iframe,iframe
            else
              write(11,"('step',i4,' image',i4)") iframe,iframe
            endif
          endif
          write(11,*) '1 0'
        endif
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
            if( .not. ALL_THRESHOLD_OFF ) then
              if(USE_OPENDX) then
                write(11,"(f11.4)") field_display(ilocnum+ieoff)
              else
                write(11,"(i10,1x,f11.4)") ireorder(ibool_number),field_display(ilocnum+ieoff)
              endif
            else
              if(USE_OPENDX) then
                write(11,"(e11.4)") field_display(ilocnum+ieoff)
              else
                ! format spezifier has problems w/ very small values
                !write(11,"(i10,1x,e7.4)") ireorder(ibool_number),field_display(ilocnum+ieoff)
                write(11,*) ireorder(ibool_number),field_display(ilocnum+ieoff)
              endif
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

    if(.not. UNIQUE_FILE) close(11)

  ! end of loop and test on all the time steps for all the movie images
  enddo

  if(UNIQUE_FILE) close(11)

  print *
  print *,'done creating movie'
  print *
  if(USE_OPENDX) print *,'DX files are stored in ', trim(OUTPUT_FILES), '/DX_*.dx'
  if(USE_AVS) print *,'AVS files are stored in ', trim(OUTPUT_FILES), '/AVS_*.inp'
  if(USE_GMT) print *,'GMT files are stored in ', trim(OUTPUT_FILES), '/gmt_*.xyz'
  print *

  end program xcreate_movie_AVS_DX


!
!=====================================================================
!

  subroutine read_AVS_DX_parameters()

  use constants
  use shared_parameters

  implicit none

  print *
  print *,'reading parameter file'
  print *

  ! read the parameter file and compute additional parameters
  call read_compute_parameters()

  ! checks
  if(.not. MOVIE_SURFACE) stop 'movie frames were not saved by the solver'

!  if(MOVIE_COARSE) stop 'create_movie_AVS_DX does not work with MOVIE_COARSE'

  end subroutine read_AVS_DX_parameters

!
!=====================================================================
!


  subroutine get_global_AVS(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  use constants

  implicit none

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
