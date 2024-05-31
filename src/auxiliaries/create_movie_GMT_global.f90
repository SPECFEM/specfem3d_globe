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

!
!---  create a movie of radial component of surface displacement/velocity in GMT format
!

  program create_movie_GMT_global

! reads in files: OUTPUT_FILES/moviedata******
!
! and creates new files: ascii_movie_*** (ASCII option)
!                     or   bin_movie_*** (binary option)
!
! these files can then be visualized using GMT, the Generic Mapping Tools
! ( http://www.soest.hawaii.edu/GMT/ )
! or the shakemovie renderer (renderOnSphere)
!
! example scripts can be found in: UTILS/Visualization/GMT/

  use constants
  use shared_parameters

  implicit none

!---------------------
! USER PARAMETER

  ! to avoid flickering in movies, the displacement/velocity field will get normalized with an
  ! averaged maximum value over the past few, available snapshots
  logical :: USE_AVERAGED_MAXIMUM = .true.
  ! minimum number of frames to average maxima
  integer,parameter :: AVERAGE_MINIMUM = 5
  ! normalizes output values
  logical, parameter :: AVERAGE_NORMALIZE_VALUES = .true.

  ! uses an absolute value given from input to normalize wavefield
  logical :: USE_ABSOLUTE_VALUE_NORMALIZATION = .false.
  double precision :: ABSOLUTE_VALUE_NORM = 0.d0

  ! muting source region
  logical :: MUTE_SOURCE = .true.
  real(kind=CUSTOM_REAL) :: RADIUS_TO_MUTE = 0.5    ! start radius in degrees
  real(kind=CUSTOM_REAL) :: STARTTIME_TO_MUTE = 0.5 ! adds seconds to shift starttime
  real(kind=CUSTOM_REAL) :: MUTE_SOURCE_MINIMAL_DISTANCE = 2.0 ! minimum taper around the source (in case of displacement movies)
  real(kind=CUSTOM_REAL) :: SURFACE_WAVE_VELOCITY = 3.5     ! speed of surface waves (km/s)
!---------------------

  integer :: i,j,it
  integer :: it1,it2
  integer :: ispec

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,displn
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord,rval,thetaval,phival,latval,longval
  real(kind=CUSTOM_REAL) :: RRval,rhoval
  real(kind=CUSTOM_REAL) :: displx,disply,displz
  real(kind=CUSTOM_REAL) :: normal_x,normal_y,normal_z
  real(kind=CUSTOM_REAL) :: thetahat_x,thetahat_y,thetahat_z
  real(kind=CUSTOM_REAL) :: phihat_x,phihat_y

  ! to average maxima over past few steps
  double precision :: min_field_current,max_field_current,max_absol,max_average
  double precision,dimension(:),allocatable :: max_history
  integer :: nmax_history,imax

  real :: disp,lat,long
  integer :: nframes,iframe,USE_COMPONENT

  character(len=MAX_STRING_LEN) :: outputname

  integer :: iproc,ipoin

! for sorting routine
  integer :: npointot,ilocnum,ielm,ieoff,ispecloc
  integer :: NIT
  double precision, dimension(:), allocatable :: xp,yp,zp,field_display

! for dynamic memory allocation
  integer :: ier

! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz

  logical :: OUTPUT_BINARY,temp_l

  real(kind=CUSTOM_REAL) :: LAT_SOURCE,LON_SOURCE,DEP_SOURCE
  real(kind=CUSTOM_REAL) :: dist_lon,dist_lat,distance,mute_factor,val
  real(kind=CUSTOM_REAL) :: cmt_hdur,cmt_t_shift,t0,hdur
  real(kind=CUSTOM_REAL) :: xmesh,ymesh,zmesh
  integer :: istamp1,istamp2

  character(len=256) :: line

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all movie frames to create a movie'
  print *,'Run this program from the directory containing directories DATA and OUTPUT_FILES'

  print *
  print *,'reading parameter file'
  print *

  ! read the parameter file and compute additional parameters
  call read_compute_parameters()

  ! get the base pathname for output files
  OUTPUT_FILES = OUTPUT_FILES_BASE

  if (.not. MOVIE_SURFACE) stop 'movie frames were not saved by the solver'

  ! user input
  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  print *,'enter first time step of movie (e.g. 1)'
  read(5,*) it1

  print *,'enter last time step of movie (e.g. ',NSTEP,'or -1 for all)'
  read(5,*) it2

  print *,'enter component (e.g. 1=Z, 2=N, 3=E)'
  read(5,*) USE_COMPONENT
  if (USE_COMPONENT < 1 .or. USE_COMPONENT > 3 ) stop 'component must be 1, 2 or 3'

  print *,'enter output ASCII (F) or binary (T)'
  read(5,*) OUTPUT_BINARY

  print *,'(optional) mute source area (T) or not (F)'
  read(5,*,iostat=ier) temp_l
  if (ier == 0) MUTE_SOURCE = temp_l

  print *,'(optional) use moving average for normalization (T) or not (F)'
  read(5,*,iostat=ier) temp_l
  if (ier == 0) USE_AVERAGED_MAXIMUM = temp_l

  print *,'(optional) enter absolute value for normalization (e.g. 1.e-7)'
  read(5,*,iostat=ier) ABSOLUTE_VALUE_NORM
  if (ier == 0) USE_ABSOLUTE_VALUE_NORMALIZATION = .true.

  print *
  print *,'--------'
  print *

  ! checks options
  if (it1 < 1 ) stop 'the first time step must be >= 1'
  if (it2 == -1 ) it2 = NSTEP

  ! user info
  print *, 'Movie surface:'
  print *, '  moviedata*** files read in directory: ',trim(OUTPUT_FILES)
  if (MOVIE_VOLUME_TYPE == 5) then
    print *, '  movie output    : displacement'
  else
    print *, '  movie output    : velocity'
  endif
  print *, '  time steps every: ',NTSTEP_BETWEEN_FRAMES
  print *
  if (MUTE_SOURCE) &
    print *, '  using mute source region'
  if (USE_AVERAGED_MAXIMUM) then
    print *, '  using averaged maximum of wavefields for scaling:'
    print *, '    averaging history over ',AVERAGE_MINIMUM,' wavefields'
    print *, '    normalizes values: ',AVERAGE_NORMALIZE_VALUES
  endif
  if (USE_ABSOLUTE_VALUE_NORMALIZATION) then
    print *, '  using absolute value for normalization: norm = ',ABSOLUTE_VALUE_NORM
  endif
  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *

  if (MOVIE_COARSE) then
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

  ! allocates movie arrays
  allocate(store_val_x(ilocnum,0:NPROCTOT-1),stat=ier)
  if (ier /= 0) stop 'Error while allocating store_val_x'

  allocate(store_val_y(ilocnum,0:NPROCTOT-1),stat=ier)
  if (ier /= 0) stop 'Error while allocating store_val_y'

  allocate(store_val_z(ilocnum,0:NPROCTOT-1),stat=ier)
  if (ier /= 0) stop 'Error while allocating store_val_z'

  allocate(store_val_ux(ilocnum,0:NPROCTOT-1),stat=ier)
  if (ier /= 0) stop 'Error while allocating store_val_ux'

  allocate(store_val_uy(ilocnum,0:NPROCTOT-1),stat=ier)
  if (ier /= 0) stop 'Error while allocating store_val_uy'

  allocate(store_val_uz(ilocnum,0:NPROCTOT-1),stat=ier)
  if (ier /= 0) stop 'Error while allocating store_val_uz'

  allocate(x(NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'Error while allocating x'

  allocate(y(NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'Error while allocating y'

  allocate(z(NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'Error while allocating z'

  allocate(displn(NGLLX,NGLLY),stat=ier)
  if (ier /= 0) stop 'Error while allocating displn'


  print *
  print *,'looping from ',it1,' to ',it2,' every ',NTSTEP_BETWEEN_FRAMES,' time steps'

  ! counts number of movie frames
  nframes = 0
  do it = it1,it2
    if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) nframes = nframes + 1
  enddo
  print *
  print *,'total number of frames will be ',nframes
  if (nframes == 0) stop 'null number of frames'

  ! maximum theoretical number of points at the surface
  if (MOVIE_COARSE) then
    npointot = NCHUNKS * NEX_XI * NEX_ETA
  else
    npointot = NCHUNKS * NEX_XI * NEX_ETA * (NGLLX-1) * (NGLLY-1)
  endif

  print *
  print *,'there are a total of ',npointot,' points on the surface.'
  print *


  print *
  print *,'Allocating 4 outputdata arrays of size 4 * CUSTOM_REAL *',npointot,'=', &
                                      4*npointot*CUSTOM_REAL/1024.0/1024.0,' MB'
  print *

  allocate(xp(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating xp'

  allocate(yp(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating yp'

  allocate(zp(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating zp'

  allocate(field_display(npointot),stat=ier)
  if (ier /= 0) stop 'Error while allocating field_display'

  ! initializes maxima history
  if (USE_AVERAGED_MAXIMUM) then
    ! determines length of history
    nmax_history = AVERAGE_MINIMUM + int( HDUR_MOVIE / (DT*NTSTEP_BETWEEN_FRAMES) * 1.5 )

    ! allocates history array
    allocate(max_history(nmax_history))
    max_history(:) = 0.0d0

    print *,'averages wavefield maxima'
    print *
    print *,'Movie half-duration: ',HDUR_MOVIE,'(s)'
    print *,'DT per time step   : ',DT,'(s)'
    print *,'Frame step size    : ',DT*NTSTEP_BETWEEN_FRAMES,'(s)'
    print *,'Normalization by averaged maxima over ',nmax_history,'snapshots'
    print *
  endif

  ! initializes
  mute_factor = 0.0_CUSTOM_REAL
  t0 = 0.0_CUSTOM_REAL

  if (MUTE_SOURCE) then
    ! initializes
    LAT_SOURCE = -1000.0
    LON_SOURCE = -1000.0
    DEP_SOURCE = 0.0
    cmt_t_shift = 0.0
    cmt_hdur = 0.0

    ! reads in source lat/lon
    if (USE_FORCE_POINT_SOURCE) then
      open(IIN,file='DATA/FORCESOLUTION',status='old',action='read',iostat=ier )
      if (ier /= 0) stop 'DATA/FORCESOLUTION not found, necessary for muting source area'
      ! skip first line
      read(IIN,*,iostat=ier ) line ! FORCE 001 ..
      ! timeshift
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(12:len_trim(line)),*) cmt_t_shift
      ! halfduration
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(15:len_trim(line)),*) cmt_hdur
      ! latitude
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(10:len_trim(line)),*) LAT_SOURCE
      ! longitude
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(11:len_trim(line)),*) LON_SOURCE
      ! depth
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(7:len_trim(line)),*) DEP_SOURCE
      close(IIN)
    else
      open(IIN,file='DATA/CMTSOLUTION',status='old',action='read',iostat=ier )
      if (ier /= 0) stop 'DATA/CMTSOLUTION not found, necessary for muting source area'
      ! skip first line, event name,timeshift,half duration
      read(IIN,*,iostat=ier ) line ! PDE line
      read(IIN,*,iostat=ier ) line ! event name
      ! timeshift
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(12:len_trim(line)),*) cmt_t_shift
      ! halfduration
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(15:len_trim(line)),*) cmt_hdur
      ! latitude
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(10:len_trim(line)),*) LAT_SOURCE
      ! longitude
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(11:len_trim(line)),*) LON_SOURCE
      ! depth
      read(IIN,'(a256)',iostat=ier ) line
      if (ier == 0 ) read(line(7:len_trim(line)),*) DEP_SOURCE
      close(IIN)
    endif
    ! effective half duration in movie runs
    hdur = real(sqrt( cmt_hdur**2 + HDUR_MOVIE**2),kind=CUSTOM_REAL)
    ! start time of simulation
    t0 = real(-1.5d0*( cmt_t_shift - hdur ),kind=CUSTOM_REAL)

    ! becomes time (s) from hypocenter to reach surface (using average 8 km/s s-wave speed)
    ! note: especially for deep sources, this helps determine a better starttime to mute
    DEP_SOURCE = DEP_SOURCE / 8.0

    ! time when muting starts
    ! note: this starttime is supposed to be the time when displacements at the surface
    !          can be observed;
    !          it helps to mute out numerical noise before the source effects actually start showing up
    STARTTIME_TO_MUTE = STARTTIME_TO_MUTE + DEP_SOURCE
    if (STARTTIME_TO_MUTE < 0.0 ) STARTTIME_TO_MUTE = 0.0

    print *,'mutes source area'
    print *
    print *,'source lat/lon/dep: ',LAT_SOURCE,LON_SOURCE,DEP_SOURCE
    print *,'muting radius: ',RADIUS_TO_MUTE,'(degrees)'
    print *,'muting starttime: ',STARTTIME_TO_MUTE,'(s)'
    print *,'simulation starttime: ',-t0,'(s)'
    print *

    ! converts values into radians
    ! colatitude [0, PI]
    LAT_SOURCE = real((90.d0 - LAT_SOURCE)*DEGREES_TO_RADIANS,kind=CUSTOM_REAL)

    ! longitude [-PI, PI]
    if (LON_SOURCE < -180.0 ) LON_SOURCE = LON_SOURCE + 360.0
    if (LON_SOURCE > 180.0 ) LON_SOURCE = LON_SOURCE - 360.0
    LON_SOURCE = real(LON_SOURCE * DEGREES_TO_RADIANS,kind=CUSTOM_REAL)

    ! mute radius in rad
    RADIUS_TO_MUTE = real(RADIUS_TO_MUTE * DEGREES_TO_RADIANS,kind=CUSTOM_REAL)
  endif

  print *,'--------'

!--- ****** read data saved by solver ******

  ! movie point locations
  outputname = "/moviedata_xyz.bin"
  open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname), &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(OUTPUT_FILES)//trim(outputname)
    stop 'Error opening moviedata file'
  endif

  ! reads in point locations
  ! (given as r theta phi for geocentric coordinate system)
  read(IOUT) store_val_x
  read(IOUT) store_val_y
  read(IOUT) store_val_z
  close(IOUT)

! --------------------------------------

  istamp1 = 0
  istamp2 = 0
  iframe = 0

! loop on all the time steps in the range entered
  do it = it1,it2

    ! check if time step corresponds to a movie frame
    if (mod(it,NTSTEP_BETWEEN_FRAMES) /= 0) cycle

    iframe = iframe + 1

    print *
    print *,'reading snapshot time step ',it,' out of ',NSTEP
    print *

    ! read all the elements from the same file
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname), &
        status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(OUTPUT_FILES)//trim(outputname)
      stop 'Error opening moviedata file'
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

    ! mutes source region
    if (MUTE_SOURCE) then
      ! initialize factor
      mute_factor = 1.0

      print *,'simulation time: ',(it-1)*DT - t0,'(s)'

      ! muting radius grows/shrinks with time
      if ((it-1)*DT - t0 > STARTTIME_TO_MUTE) then

        ! approximate wavefront travel distance in degrees
        ! (~3.5 km/s wave speed for surface waves)
        distance = real(SURFACE_WAVE_VELOCITY * ((it-1)*DT-t0) / (R_PLANET/1000.d0) * RADIANS_TO_DEGREES,kind=CUSTOM_REAL)

        print *,'distance approximate: ',distance,'(degrees)'

        ! approximate distance to source (in degrees)
        ! (shrinks if waves travel back from antipode)
        !do while ( distance > 360. )
        !  distance = distance - 360.
        !enddo
        ! waves are back at origin, no source tapering anymore
        if (distance > 360.0 ) distance = 0.0

        ! shrinks when waves reached antipode
        !if (distance > 180. ) distance = 360. - distance
        ! shrinks when waves reached half-way to antipode
        if (distance > 90.0 ) distance = 90.0 - distance

        ! limit size around source (in degrees)
        if (distance < 0.0 ) distance = 0.0
        if (distance > 80.0 ) distance = 80.0

        ! minimal taper
        if (MOVIE_VOLUME_TYPE == 5) then
          ! in case movie shows displacement, there will be a static offset at the source which will tamper all waves
          if (distance <= TINYVAL) distance = MUTE_SOURCE_MINIMAL_DISTANCE
        endif

        print *,'muting radius: ',0.7 * distance,'(degrees)'

        ! new radius of mute area (in rad)
        RADIUS_TO_MUTE = real(0.7d0 * distance * DEGREES_TO_RADIANS,kind=CUSTOM_REAL)
      else
        ! mute_factor used at the beginning for scaling displacement values
        if (STARTTIME_TO_MUTE > TINYVAL) then
          ! mute factor 1: no masking out
          !                     0: masks out values (within mute radius)
          ! linear scaling between [0,1]:
          ! from 0 (simulation time equal to zero )
          ! to 1 (simulation time equals starttime_to_mute)
          mute_factor = real(1.d0 - ( STARTTIME_TO_MUTE - ((it-1)*DT-t0) ) / (STARTTIME_TO_MUTE+t0),kind=CUSTOM_REAL)
          ! threshold value for mute_factor
          if (mute_factor < TINYVAL ) mute_factor = real(TINYVAL,kind=CUSTOM_REAL)
          if (mute_factor > 1.0 ) mute_factor = 1.0
        endif
      endif
    endif

    ! clear number of elements kept
    ispec = 0

    print *,'Converting to geo-coordinates'

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
                if (USE_COMPONENT == 1) then
                   ! compute unit normal vector to the surface
                   RRval = sqrt(xcoord**2 + ycoord**2 + zcoord**2)
                   if (RRval < 1.e-10 ) stop 'Error unit normal vector'
                   normal_x = xcoord / RRval
                   normal_y = ycoord / RRval
                   normal_z = zcoord / RRval

                   displn(i,j) = displx*normal_x + disply*normal_y + displz*normal_z

                else if (USE_COMPONENT == 2) then

                   ! compute unit tangent vector to the surface (N-S)
                   RRval = sqrt(xcoord**2 + ycoord**2 + zcoord**2)
                   if (RRval < 1.e-10 ) stop 'Error unit normal vector'
                   rhoval = sqrt(xcoord**2 + ycoord**2)
                   if (rhoval < 1.e-10) then
                    ! location at pole
                    thetahat_x = 0.0
                    thetahat_y = 0.0
                   else
                    thetahat_x = (zcoord*xcoord) / (rhoval*RRval)
                    thetahat_y = (zcoord*ycoord) / (rhoval*RRval)
                   endif
                   thetahat_z = - rhoval/RRval

                   displn(i,j) = - (displx*thetahat_x + disply*thetahat_y + displz*thetahat_z)

                else if (USE_COMPONENT == 3) then

                   ! compute unit tangent to the surface (E-W)
                   rhoval = sqrt(xcoord**2 + ycoord**2)
                   if (rhoval < 1.e-10) then
                    ! location at pole
                    phihat_x = 0.0
                    phihat_y = 0.0
                   else
                    phihat_x = -ycoord / rhoval
                    phihat_y = xcoord / rhoval
                   endif

                   displn(i,j) = displx*phihat_x + disply*phihat_y
                endif

                ! mute values
                if (MUTE_SOURCE) then

                  ! distance in colatitude (in rad)
                  ! note: this mixes geocentric (point location) and geographic (source location) coordinates;
                  !          since we only need approximate distances here,
                  !          this should be fine for the muting region
                  dist_lat = thetaval - LAT_SOURCE

                  ! distance in longitude (in rad)
                  ! checks source longitude range
                  if (LON_SOURCE - RADIUS_TO_MUTE < -PI .or. LON_SOURCE + RADIUS_TO_MUTE > PI) then
                    ! source close to 180. longitudes, shifts range to [0, 2PI]
                    if (phival < 0.0 ) phival = phival + real(TWO_PI,kind=CUSTOM_REAL)
                    if (LON_SOURCE < 0.0) then
                      dist_lon = phival - real((LON_SOURCE + TWO_PI),kind=CUSTOM_REAL)
                    else
                      dist_lon = phival - LON_SOURCE
                    endif
                  else
                    ! source well between range to [-PI, PI]
                    ! shifts phival to be like LON_SOURCE between [-PI,PI]
                    if (phival > PI ) phival = phival - real(TWO_PI,kind=CUSTOM_REAL)
                    if (phival < -PI ) phival = phival + real(TWO_PI,kind=CUSTOM_REAL)

                    dist_lon = phival - LON_SOURCE
                  endif
                  ! distance of point to source (in rad)
                  distance = sqrt(dist_lat**2 + dist_lon**2)

                  ! mutes source region values
                  if (distance < RADIUS_TO_MUTE) then
                    ! muting takes account of the event time
                    if ((it-1)*DT-t0 > STARTTIME_TO_MUTE) then
                      ! wavefield will be tapered to mask out noise in source area
                      ! factor from 0 to 1
                      mute_factor = real( ( 0.d5*(1.d0 - cos(distance/RADIUS_TO_MUTE*PI)) )**6,kind=CUSTOM_REAL)
                      ! factor from 0.01 to 1
                      mute_factor = mute_factor * 0.99 + 0.01
                      displn(i,j) = displn(i,j) * mute_factor
                    else
                      ! wavefield will initially be scaled down to avoid noise being amplified at beginning
                      displn(i,j) = displn(i,j) * mute_factor
                    endif
                  endif

                endif

             enddo !i
          enddo  !j

          ispec = ispec + 1
          if (MOVIE_COARSE) then
            ielm = ispec-1
          else
            ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
          endif

          do j = 1,NGLLY-NIT
             do i = 1,NGLLX-NIT
                ! offset (starts at 1)
                if (MOVIE_COARSE) then
                  ieoff = ielm+1
                else
                  ieoff = (ielm+(i-1)+(j-1)*(NGLLX-1))+1
                endif

                ! for movie_coarse e.g. x(i,j) is defined at x(1,1), x(1,NGLLY), x(NGLLX,1) and x(NGLLX,NGLLY)
                ! be aware that for the cubed sphere, the mapping changes for different chunks,
                ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates
                if (MOVIE_COARSE) then
                  if (NCHUNKS == 6) then
                    ! chunks mapped such that element corners increase in long/lat
                    select case (iproc/NPROC+1)
                      case (CHUNK_AB)
                        xp(ieoff) = dble(x(1,NGLLY))
                        yp(ieoff) = dble(y(1,NGLLY))
                        zp(ieoff) = dble(z(1,NGLLY))
                        field_display(ieoff) = dble(displn(1,NGLLY))
                      case (CHUNK_AB_ANTIPODE)
                        xp(ieoff) = dble(x(1,1))
                        yp(ieoff) = dble(y(1,1))
                        zp(ieoff) = dble(z(1,1))
                        field_display(ieoff) = dble(displn(1,1))
                      case (CHUNK_AC)
                        xp(ieoff) = dble(x(1,NGLLY))
                        yp(ieoff) = dble(y(1,NGLLY))
                        zp(ieoff) = dble(z(1,NGLLY))
                        field_display(ieoff) = dble(displn(1,NGLLY))
                      case (CHUNK_AC_ANTIPODE)
                        xp(ieoff) = dble(x(1,1))
                        yp(ieoff) = dble(y(1,1))
                        zp(ieoff) = dble(z(1,1))
                        field_display(ieoff) = dble(displn(1,1))
                      case (CHUNK_BC)
                        xp(ieoff) = dble(x(1,NGLLY))
                        yp(ieoff) = dble(y(1,NGLLY))
                        zp(ieoff) = dble(z(1,NGLLY))
                        field_display(ieoff) = dble(displn(1,NGLLY))
                      case (CHUNK_BC_ANTIPODE)
                        xp(ieoff) = dble(x(NGLLX,NGLLY))
                        yp(ieoff) = dble(y(NGLLX,NGLLY))
                        zp(ieoff) = dble(z(NGLLX,NGLLY))
                        field_display(ieoff) = dble(displn(NGLLX,NGLLY))
                      case default
                        stop 'incorrect chunk number'
                    end select
                  else
                    ! takes lower-left point only
                    xp(ieoff) = dble(x(1,1))
                    yp(ieoff) = dble(y(1,1))
                    zp(ieoff) = dble(z(1,1))
                    field_display(ieoff) = dble(displn(1,1))
                  endif ! NCHUNKS
                else
                  xp(ieoff) = dble(x(i,j))
                  yp(ieoff) = dble(y(i,j))
                  zp(ieoff) = dble(z(i,j))
                  field_display(ieoff) = dble(displn(i,j))
                endif ! MOVIE_COARSE

                ! determines North / South pole index for stamping maximum values
                if (USE_AVERAGED_MAXIMUM .and. AVERAGE_NORMALIZE_VALUES) then
                  xmesh = real(xp(ieoff),kind=CUSTOM_REAL)
                  ymesh = real(yp(ieoff),kind=CUSTOM_REAL)
                  zmesh = real(zp(ieoff),kind=CUSTOM_REAL)
                  if (zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = - real(SMALL_VAL_ANGLE,kind=CUSTOM_REAL)
                  if (zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = real(SMALL_VAL_ANGLE,kind=CUSTOM_REAL)
                  thetaval = atan2(sqrt(xmesh*xmesh+ymesh*ymesh),zmesh)
                  ! thetaval between 0 and PI / 2
                  !print *,'thetaval:',thetaval * 180. / PI
                  ! close to north pole
                  if (thetaval >= 0.495 * PI ) istamp1 = ieoff
                  ! close to south pole
                  if (thetaval <= 0.01 ) istamp2 = ieoff
                endif

             enddo !i
          enddo  !j

       enddo !ispec

    enddo !nproc

    ! compute min and max of data value to normalize
    min_field_current = minval(field_display(:))
    max_field_current = maxval(field_display(:))

    ! print minimum and maximum amplitude in current snapshot
    print *
    print *,'minimum amplitude in current snapshot = ',min_field_current
    print *,'maximum amplitude in current snapshot = ',max_field_current
    print *

    ! normalizes with an absolute value
    if (USE_ABSOLUTE_VALUE_NORMALIZATION) then
      ! info
      print *,'normalizes with maximum absolute value = ',ABSOLUTE_VALUE_NORM

      ! thresholds positive & negative maximum values
      where( field_display(:) > ABSOLUTE_VALUE_NORM ) field_display = ABSOLUTE_VALUE_NORM
      where( field_display(:) < - ABSOLUTE_VALUE_NORM ) field_display = - ABSOLUTE_VALUE_NORM

      ! normalizes wavefield
      if (abs(ABSOLUTE_VALUE_NORM) > 0.d0) field_display = field_display / ABSOLUTE_VALUE_NORM

      ! re-computes min and max of data value
      min_field_current = minval(field_display(:))
      max_field_current = maxval(field_display(:))
      print *,'  new minimum amplitude in current snapshot = ',min_field_current
      print *,'  new maximum amplitude in current snapshot = ',max_field_current
      print *
    endif

    ! takes average over last few snapshots available and uses it
    ! to normalize field values
    if (USE_AVERAGED_MAXIMUM) then

      ! (average) maximum between positive and negative values
      max_absol = (abs(min_field_current)+abs(max_field_current))/2.0

      ! stores last few maxima
      ! index between 1 and nmax_history
      imax = mod(iframe-1,nmax_history) + 1
      max_history( imax ) = max_absol

      ! average over history
      max_average = sum( max_history )
      if (iframe < nmax_history) then
        ! history not filled yet, only average over available entries
        max_average = max_average / iframe
      else
        ! average over all history entries
        max_average = max_average / nmax_history
      endif

      print *,'maximum amplitude averaged in current snapshot = ',max_absol
      print *,'maximum amplitude over averaged last snapshots = ',max_average
      print *

      ! thresholds positive & negative maximum values
      where( field_display(:) > max_absol ) field_display = max_absol
      where( field_display(:) < - max_absol ) field_display = -max_absol

      ! sets new maxima for decaying wavefield
      ! this should avoid flickering when normalizing wavefields
      if (AVERAGE_NORMALIZE_VALUES) then
        ! checks stamp indices for maximum values
        if (istamp1 == 0 ) istamp1 = ieoff
        if (istamp2 == 0 ) istamp2 = ieoff-1
        !print *, 'stamp: ',istamp1,istamp2

        if (max_absol < max_average) then
          ! distance (in degree) of surface waves travelled
          distance = real(SURFACE_WAVE_VELOCITY * ((it-1)*DT-t0) / (R_PLANET/1000.d0) * RADIANS_TO_DEGREES,kind=CUSTOM_REAL)
          if (distance > 10.0 .and. distance <= 20.0) then
            ! smooth transition between 10 and 20 degrees
            ! sets positive and negative maximum
            field_display(istamp1) = + max_absol + (max_average-max_absol) * (distance - 10.d0)/10.d0
            field_display(istamp2) = - ( max_absol + (max_average-max_absol) * (distance - 10.d0)/10.d0 )
          else if (distance > 20.0) then
            ! sets positive and negative maximum
            field_display(istamp1) = + max_average
            field_display(istamp2) = - max_average
          endif
        else
          ! thresholds positive & negative maximum values
          where( field_display(:) > max_average ) field_display = max_average
          where( field_display(:) < - max_average ) field_display = -max_average
          ! sets positive and negative maximum
          field_display(istamp1) = + max_average
          field_display(istamp2) = - max_average
        endif
        ! updates current wavefield maxima
        min_field_current = minval(field_display(:))
        max_field_current = maxval(field_display(:))
        max_absol = (abs(min_field_current)+abs(max_field_current))/2.d0
      endif

      ! scales field values up to match average
      if (abs(max_absol) > TINYVAL) &
        field_display = field_display * max_average / max_absol

      ! thresholds after scaling positive & negative maximum values
      where( field_display(:) > max_average ) field_display = max_average
      where( field_display(:) < - max_average ) field_display = -max_average

      ! normalizes field values
      if (AVERAGE_NORMALIZE_VALUES) then
        if (MUTE_SOURCE) then
          ! checks if source wavefield kicked in
          if ((it-1)*DT - t0 > STARTTIME_TO_MUTE) then
            ! wavefield should be visible at surface now
            ! normalizes wavefield
            if (abs(max_average) > TINYVAL ) field_display = field_display / max_average
          else
            ! no wavefield yet assumed

            ! we set two single field values (last in array)
            ! to: +/- 100 * max_average
            ! to avoid further amplifying when
            ! a normalization routine is used for rendering images
            ! (which for example is the case for shakemovies)
            if (STARTTIME_TO_MUTE > TINYVAL) then
              ! with additional scale factor:
              ! linear scaling between [0,1]:
              ! from 0 (simulation time equal to -t0 )
              ! to 1 (simulation time equals starttime_to_mute)
              mute_factor = real(1.d0 - ( STARTTIME_TO_MUTE - ((it-1)*DT-t0) ) / (STARTTIME_TO_MUTE+t0),kind=CUSTOM_REAL)
              ! takes complement and shifts scale to (1,100)
              ! thus, mute factor is 100 at simulation start and 1.0 at starttime_to_mute
              mute_factor = abs(1.0 - mute_factor) * 99.0 + 1.0
              ! positive and negative maximum reach average when wavefield appears
              val = real(mute_factor * max_average,kind=CUSTOM_REAL)
            else
              ! uses a constant factor
              val = real(100.d0 * max_average,kind=CUSTOM_REAL)
            endif
            ! positive and negative maximum
            field_display(istamp1) = + val
            field_display(istamp2) = - val
            if (abs(max_average) > TINYVAL ) field_display = field_display / val
          endif
        else
          ! no source to mute
          ! normalizes wavefield
          if (abs(max_average) > TINYVAL ) field_display = field_display / max_average
        endif
      endif

      ! re-computes min and max of data value
      min_field_current = minval(field_display(:))
      max_field_current = maxval(field_display(:))
      print *,'  new minimum amplitude in current snapshot = ',min_field_current
      print *,'  new maximum amplitude in current snapshot = ',max_field_current
      print *
    endif

    print *,'initial number of points (with multiples) was ',npointot
    print *,'final number of points is                     ',ieoff

    !--- ****** create GMT file ******

    ! create file name and open file
    if (OUTPUT_BINARY) then
      if (USE_COMPONENT == 1) then
       write(outputname,"('/bin_movie_',i6.6,'.d')") it
      else if (USE_COMPONENT == 2) then
       write(outputname,"('/bin_movie_',i6.6,'.N')") it
      else if (USE_COMPONENT == 3) then
       write(outputname,"('/bin_movie_',i6.6,'.E')") it
      endif
      open(unit=11,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown', &
            form='unformatted',action='write',iostat=ier)
      if (ier /= 0) stop 'Error opening bin_movie file'

      if (iframe == 1) then
        open(unit=12,file=trim(OUTPUT_FILES)//'/bin_movie.xy',status='unknown', &
            form='unformatted',action='write',iostat=ier)
        if (ier /= 0) stop 'Error opening bin_movie.xy file'
      endif
    else
      if (USE_COMPONENT == 1) then
       write(outputname,"('/ascii_movie_',i6.6,'.d')") it
      else if (USE_COMPONENT == 2) then
       write(outputname,"('/ascii_movie_',i6.6,'.N')") it
      else if (USE_COMPONENT == 3) then
       write(outputname,"('/ascii_movie_',i6.6,'.E')") it
      endif
      open(unit=11,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown', &
            action='write',iostat=ier)
      if (ier /= 0) stop 'Error opening ascii_movie file'

      if (iframe == 1) then
        open(unit=12,file=trim(OUTPUT_FILES)//'/ascii_movie.xy',status='unknown', &
            action='write',iostat=ier)
        if (ier /= 0) stop 'Error opening ascii_movie.xy file'
      endif
    endif
    print *
    print *,'Writing output: ',trim(OUTPUT_FILES)//trim(outputname)
    print *

    ! clear number of elements kept
    ispec = 0

    ! read points for all the slices
    do iproc = 0,NPROCTOT-1
      do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA
        ispec = ispec + 1
        if (MOVIE_COARSE) then
          ielm = ispec - 1
        else
          ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
        endif

        do j = 1,NGLLY-NIT
          do i = 1,NGLLX-NIT
            if (MOVIE_COARSE) then
              ieoff = ielm + 1
            else
              ieoff = (ielm+(i-1)+(j-1)*(NGLLX-1))+1
            endif

            ! point position
            if (iframe == 1) then
              ! gets Cartesian coordinates
              xcoord = sngl(xp(ieoff))
              ycoord = sngl(yp(ieoff))
              zcoord = sngl(zp(ieoff))

              ! location latitude/longitude (with geographic latitude in degrees)
              call xyz_2_rlatlon_cr(xcoord,ycoord,zcoord,rval,latval,longval,ELLIPTICITY)

              ! converts to real
              lat = latval
              long = longval

              ! puts lon in range [-180,180]
              if (long > 180.0) long = long-360.0
            endif

            ! displacement
            disp = sngl(field_display(ieoff))

            ! writes displacement and longitude/latitude to corresponding files
            if (OUTPUT_BINARY) then
              ! binary output
              write(11) disp
              if (iframe == 1) write(12) long,lat
            else
              ! ascii output
              write(11,*) disp
              if (iframe == 1) write(12,'(2f18.6)') long,lat
            endif

          enddo !i
        enddo !j
      enddo !ispecloc
    enddo !iproc
    close(11)
    if (iframe == 1) close(12)

  ! end of loop and test on all the time steps for all the movie images
  enddo

  print *,'done creating movie'
  print *,'GMT ASCII files are stored in ascii_movie_*.{xy,d,E,N}'
  print *,'binary files are stored in bin_movie_*.{xy,d,E,N}'

  end program create_movie_GMT_global


