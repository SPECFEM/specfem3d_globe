!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
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
!---  create a movie of radial component of surface displacement in GMT format
!

  program create_movie_GMT_global

! reads in files: OUTPUT_FILES/moviedata******
!
! and creates new files: ascii_movie_*** (ascii option) /or/ bin_movie_*** (binary option)
!
! these files can then be visualized using GMT, the Generic Mapping Tools
! ( http://www.soest.hawaii.edu/GMT/ )
!
! example scripts can be found in: UTILS/Visualization/GMT/

  implicit none

  include "constants.h"

!---------------------
! USER PARAMETER

  ! to avoid flickering in movies, the displacement field will get normalized with an 
  ! averaged maximum value over the past few, available snapshots
  logical,parameter :: USE_AVERAGED_MAXIMUM = .true.
  
  ! minimum number of frames to average maxima
  integer,parameter :: AVERAGE_MINIMUM = 5

  ! muting source region
  logical, parameter :: MUTE_SOURCE = .true.
  real(kind=CUSTOM_REAL) :: RADIUS_TO_MUTE = 1.0    ! start radius in degrees  
  real(kind=CUSTOM_REAL) :: STARTTIME_TO_MUTE = 2.0 ! factor times hdur_movie

  ! normalizes output values
  logical, parameter :: NORMALIZE_VALUES = .true.
  
!---------------------

  integer i,j,it
  integer it1,it2
  integer ispec

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,displn
  real(kind=CUSTOM_REAL) xcoord,ycoord,zcoord,rval,thetaval,phival
  real(kind=CUSTOM_REAL) RRval,rhoval
  real(kind=CUSTOM_REAL) displx,disply,displz
  real(kind=CUSTOM_REAL) normal_x,normal_y,normal_z
  real(kind=CUSTOM_REAL) thetahat_x,thetahat_y,thetahat_z
  real(kind=CUSTOM_REAL) phihat_x,phihat_y
  
  ! to average maxima over past few steps
  double precision min_field_current,max_field_current,max_absol,max_average
  double precision,dimension(:),allocatable :: max_history
  integer :: nmax_history,imax

  real disp,lat,long
  integer nframes,iframe,USE_COMPONENT

  character(len=150) outputname

  integer iproc,ipoin

! for sorting routine
  integer npointot,ilocnum,ielm,ieoff,ispecloc,NIT
  double precision, dimension(:), allocatable :: xp,yp,zp,field_display

! for dynamic memory allocation
  integer ierror

! movie files stored by solver
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
         store_val_x,store_val_y,store_val_z, &
         store_val_ux,store_val_uy,store_val_uz

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,RMOHO_FICTITIOUS_IN_MESHER

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               NGLOB


  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs
  logical :: CASE_3D,OUTPUT_BINARY

  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  real(kind=CUSTOM_REAL) :: LAT_SOURCE,LON_SOURCE,DEP_SOURCE 
  real(kind=CUSTOM_REAL) :: dist_lon,dist_lat,mute_factor
  character(len=256) line

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all movie frames to create a movie'
  print *,'Run this program from the directory containing directories DATA and OUTPUT_FILES'

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
         NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
         ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
         CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
         RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
         R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
         MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
         TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
         ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
         MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
         PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
         ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
         NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
         DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
         WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY)

  if(.not. MOVIE_SURFACE) stop 'movie frames were not saved by the solver'

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *
  if(MOVIE_COARSE) then
    ! note:
    ! nex_per_proc_xi*nex_per_proc_eta = nex_xi*nex_eta/nproc = nspec2d_top(iregion_crust_mantle) used in specfem3D.f90
    ! and ilocnum = nmovie_points = 2 * 2 * NEX_XI * NEX_ETA / NPROC
    ilocnum = 2 * 2 * NEX_PER_PROC_XI*NEX_PER_PROC_ETA
    NIT =NGLLX-1
  else
    ilocnum = NGLLX*NGLLY*NEX_PER_PROC_XI*NEX_PER_PROC_ETA
    NIT = 1
  endif
  print *
  print *,'Allocating arrays for reading data of size ',ilocnum*NPROCTOT,'=',6*ilocnum*NPROCTOT*CUSTOM_REAL/1000000,'MB'
  print *

  ! allocates movie arrays
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

  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  ! user input
  print *,'--------'
  print *,'enter first time step of movie (e.g. 1)'
  read(5,*) it1

  print *,'enter last time step of movie (e.g. ',NSTEP,'or -1 for all)'
  read(5,*) it2
  
  print *,'enter component (e.g. 1=Z, 2=N, 3=E)'
  read(5,*) USE_COMPONENT

  print *,'enter output ascii (F) or binary (T)'
  read(5,*) OUTPUT_BINARY
  print *,'--------'

  ! checks options
  if( it2 == -1 ) it2 = NSTEP

  print *
  print *,'looping from ',it1,' to ',it2,' every ',NTSTEP_BETWEEN_FRAMES,' time steps'

  ! counts number of movie frames
  nframes = 0
  do it = it1,it2
    if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) nframes = nframes + 1
  enddo
  print *
  print *,'total number of frames will be ',nframes
  if(nframes == 0) stop 'null number of frames'

  ! maximum theoretical number of points at the surface
  if(MOVIE_COARSE) then
    npointot = NCHUNKS * NEX_XI * NEX_ETA
  else
    npointot = NCHUNKS * NEX_XI * NEX_ETA * (NGLLX-1) * (NGLLY-1)
  endif

  print *
  print *,'there are a total of ',npointot,' points on the surface.'
  print *


  print *
  print *,'Allocating 4 outputdata arrays of size 4*CUSTOM_REAL',npointot,'=',4*npointot*CUSTOM_REAL/1000000,' MB'
  print *

  allocate(xp(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating xp'

  allocate(yp(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating yp'

  allocate(zp(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating zp'

  allocate(field_display(npointot),stat=ierror)
  if(ierror /= 0) stop 'error while allocating field_display'


  ! initializes maxima history
  if( USE_AVERAGED_MAXIMUM ) then
    ! determines length of history
    nmax_history = AVERAGE_MINIMUM + int( HDUR_MOVIE / (DT*NTSTEP_BETWEEN_FRAMES) * 1.5 )
        
    ! allocates history array
    allocate(max_history(nmax_history))
    max_history(:) = 0.0d0

    print *
    print *,'Movie half-duration: ',HDUR_MOVIE,'(s)'
    print *,'Frame step size    : ',DT*NTSTEP_BETWEEN_FRAMES,'(s)'
    print *,'Normalization by averaged maxima over ',nmax_history,'snapshots'
    print *
    
    if( MUTE_SOURCE ) then
      ! initializes
      LAT_SOURCE = -1000.0
      LON_SOURCE = -1000.0
      
      ! reads in source lat/lon
      open(22,file="DATA/CMTSOLUTION",status='old',action='read',iostat=ierror )
      if( ierror == 0 ) then
        ! skip first line, event name,timeshift,half duration
        read(22,*,iostat=ierror ) line ! PDE line
        read(22,*,iostat=ierror ) line ! event name
        read(22,*,iostat=ierror ) line ! timeshift 
        read(22,*,iostat=ierror ) line ! halfduration 
        ! latitude
        read(22,'(a256)',iostat=ierror ) line
        if( ierror == 0 ) read(line(10:len_trim(line)),*) LAT_SOURCE
        ! longitude
        read(22,'(a256)',iostat=ierror ) line
        if( ierror == 0 ) read(line(11:len_trim(line)),*) LON_SOURCE
        ! depth
        read(22,'(a256)',iostat=ierror ) line
        if( ierror == 0 ) read(line(11:len_trim(line)),*) DEP_SOURCE
        close(22)
      endif
      
      print *,'muting source lat/lon/dep: ',LAT_SOURCE,LON_SOURCE,DEP_SOURCE
      
      ! becomes time (s) from hypocenter to reach surface (using average 8 km/s p-wave speed)
      DEP_SOURCE = DEP_SOURCE / 8.0
      
      ! time when muting starts 
      STARTTIME_TO_MUTE = STARTTIME_TO_MUTE * HDUR_MOVIE + DEP_SOURCE
      
      print *,'muting radius: ',RADIUS_TO_MUTE
      print *,'muting starttime: ',STARTTIME_TO_MUTE,'(s)'
      print *
      
      ! colatitude [0, PI]
      LAT_SOURCE = (90. - LAT_SOURCE)*PI/180.0
      
      ! longitude [-PI, PI]
      if( LON_SOURCE < -180.0 ) LON_SOURCE = LON_SOURCE + 360.0
      if( LON_SOURCE > 180.0 ) LON_SOURCE = LON_SOURCE - 360.0
      LON_SOURCE = LON_SOURCE *PI/180.0
      
      ! mute radius in rad
      RADIUS_TO_MUTE = RADIUS_TO_MUTE*PI/180.0
    endif
    
    
  endif
  print *,'--------'

!--- ****** read data saved by solver ******

! --------------------------------------

  iframe = 0

! loop on all the time steps in the range entered
  do it = it1,it2
     ! check if time step corresponds to a movie frame
     if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

        iframe = iframe + 1

        ! mutes source region
        if( MUTE_SOURCE ) then

          ! muting radius grows/shrinks with time
          if( (it-1)*DT > STARTTIME_TO_MUTE  ) then                          

            ! approximate wavefront travel distance in degrees (~3.5 km/s wave speed for surface waves)
            mute_factor = 3.5 * (it-1)*DT / 6371. * 180./PI
            
            ! approximate distance to source (in degrees)
            do while ( mute_factor > 360. )
              mute_factor = mute_factor - 360.
            enddo
            if( mute_factor > 180. ) mute_factor = 360. - mute_factor

            ! limit size around source (in degrees)
            !if( mute_factor < 10. ) then
            !  mute_factor = 0.0
            !endif
            if( mute_factor > 80. ) then
              mute_factor = 80.0
            endif
            
            print*,'muting radius: ',0.7 * mute_factor
                        
            RADIUS_TO_MUTE = 0.7 * mute_factor * PI/180.
            
          else
            ! mute_factor used at the beginning for scaling displacement values
            if( STARTTIME_TO_MUTE > TINYVAL ) then
              ! scales from 1 to 0
              mute_factor = ( STARTTIME_TO_MUTE - (it-1)*DT ) / STARTTIME_TO_MUTE  
              if( mute_factor < TINYVAL ) mute_factor = TINYVAL
            else
              mute_factor = 1.0
            endif          
          endif
          
        endif

        ! read all the elements from the same file
        write(outputname,"('OUTPUT_FILES/moviedata',i6.6)") it
        open(unit=IOUT,file=outputname,status='old',form='unformatted')

        print *
        print *,'reading snapshot time step ',it,' out of ',NSTEP,' file ',outputname
        !print *

        ! reads in point locations
        ! (given as r theta phi for geocentric coordinate system)
        read(IOUT) store_val_x
        read(IOUT) store_val_y
        read(IOUT) store_val_z
        
        ! reads in associated values (velocity..)
        read(IOUT) store_val_ux
        read(IOUT) store_val_uy
        read(IOUT) store_val_uz
        
        close(IOUT)
        !print *, 'finished reading ',outputname

        ! clear number of elements kept
        ispec = 0

        ! read points for all the slices
        print *,'Converting to geo-coordinates'
        do iproc = 0,NPROCTOT-1
           ! reset point number
           ipoin = 0
           do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA
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


                    if(USE_COMPONENT == 1) then
                       ! compute unit normal vector to the surface
                       RRval = sqrt(xcoord**2 + ycoord**2 + zcoord**2)
                       normal_x = xcoord / RRval
                       normal_y = ycoord / RRval
                       normal_z = zcoord / RRval

                       displn(i,j) = displx*normal_x   + disply*normal_y   + displz*normal_z

                    elseif(USE_COMPONENT == 2) then

                       ! compute unit tangent vector to the surface (N-S)
                       RRval = sqrt(xcoord**2 + ycoord**2 + zcoord**2)
                       rhoval = sqrt(xcoord**2 + ycoord**2)
                       thetahat_x = (zcoord*xcoord) / (rhoval*RRval)
                       thetahat_y = (zcoord*ycoord) / (rhoval*RRval)
                       thetahat_z = - rhoval/RRval

                       displn(i,j) = - (displx*thetahat_x + disply*thetahat_y + displz*thetahat_z)
                    elseif(USE_COMPONENT == 3) then

                       ! compute unit tangent to the surface (E-W)
                       rhoval = sqrt(xcoord**2 + ycoord**2)
                       phihat_x = -ycoord / rhoval
                       phihat_y = xcoord / rhoval

                       displn(i,j) = displx*phihat_x   + disply*phihat_y
                    endif
                    
                    
                    ! mute values
                    if( MUTE_SOURCE ) then
                      
                      ! distance in colatitude                                                
                      ! note: this mixes geocentric (point location) and geographic (source location) coordinates;
                      !          since we only need approximate distances here, this should be fine for the muting region  
                      dist_lat = thetaval - LAT_SOURCE
                      
                      ! distance in longitude
                      ! checks source longitude range
                      if( LON_SOURCE - RADIUS_TO_MUTE < -PI .or. LON_SOURCE + RADIUS_TO_MUTE > PI ) then
                        ! source close to 180. longitudes, shifts range to [0, 2PI]
                        if( phival < 0.0 ) phival = phival + 2.0*PI
                        if( LON_SOURCE < 0.0 ) then
                          dist_lon = phival - (LON_SOURCE + 2.0*PI)
                        else
                          dist_lon = phival - LON_SOURCE
                        endif
                      else
                        ! source well between range to [-PI, PI]
                        ! shifts phival to be like LON_SOURCE between [-PI,PI]
                        if( phival > PI ) phival = phival - 2.0*PI
                        if( phival < -PI ) phival = phival + 2.0*PI
                        
                        dist_lon = phival - LON_SOURCE
                      endif
                      
                      ! mutes source region values
                      if ( ( dist_lat**2 + dist_lon**2 ) < RADIUS_TO_MUTE**2 ) then                          
                        ! muting takes account of the event time
                        if( (it-1)*DT > STARTTIME_TO_MUTE  ) then                          
                          displn(i,j) = displn(i,j) * TINYVAL
                        else
                          displn(i,j) = displn(i,j) * mute_factor                          
                        endif
                      endif
                        
                    endif
                    
                    
                 enddo !i
              enddo  !j

              ispec = ispec + 1
              if(MOVIE_COARSE) then
                ielm = ispec-1
              else
                ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
              endif
              do j = 1,NGLLY-NIT
                 do i = 1,NGLLX-NIT
                    if(MOVIE_COARSE) then
                      ieoff = ielm+1
                    else
                      ieoff = (ielm+(i-1)+(j-1)*(NGLLX-1))+1
                    endif

! for movie_coarse e.g. x(i,j) is defined at x(1,1), x(1,NGLLY), x(NGLLX,1) and x(NGLLX,NGLLY)
! be aware that for the cubed sphere, the mapping changes for different chunks,
! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates
                    if(MOVIE_COARSE) then
                      if(NCHUNKS == 6) then
                        ! chunks mapped such that element corners increase in long/lat
                        select case (iproc/NPROC+1)
                          case(CHUNK_AB)
                            xp(ieoff) = dble(x(1,NGLLY))
                            yp(ieoff) = dble(y(1,NGLLY))
                            zp(ieoff) = dble(z(1,NGLLY))
                            field_display(ieoff) = dble(displn(1,NGLLY))
                          case(CHUNK_AB_ANTIPODE)
                            xp(ieoff) = dble(x(1,1))
                            yp(ieoff) = dble(y(1,1))
                            zp(ieoff) = dble(z(1,1))
                            field_display(ieoff) = dble(displn(1,1))
                          case(CHUNK_AC)
                            xp(ieoff) = dble(x(1,NGLLY))
                            yp(ieoff) = dble(y(1,NGLLY))
                            zp(ieoff) = dble(z(1,NGLLY))
                            field_display(ieoff) = dble(displn(1,NGLLY))
                          case(CHUNK_AC_ANTIPODE)
                            xp(ieoff) = dble(x(1,1))
                            yp(ieoff) = dble(y(1,1))
                            zp(ieoff) = dble(z(1,1))
                            field_display(ieoff) = dble(displn(1,1))
                          case(CHUNK_BC)
                            xp(ieoff) = dble(x(1,NGLLY))
                            yp(ieoff) = dble(y(1,NGLLY))
                            zp(ieoff) = dble(z(1,NGLLY))
                            field_display(ieoff) = dble(displn(1,NGLLY))
                          case(CHUNK_BC_ANTIPODE)
                            xp(ieoff) = dble(x(NGLLX,NGLLY))
                            yp(ieoff) = dble(y(NGLLX,NGLLY))
                            zp(ieoff) = dble(z(NGLLX,NGLLY))
                            field_display(ieoff) = dble(displn(NGLLX,NGLLY))
                          case default
                            stop 'incorrect chunk number'
                        end select
                      else
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

        ! takes average over last few snapshots available and uses it
        ! to normalize field values
        if( USE_AVERAGED_MAXIMUM ) then

          ! (average) maximum between positive and negative values
          max_absol = (abs(min_field_current)+abs(max_field_current))/2.0
        
          ! stores last few maxima
          ! index between 1 and nmax_history
          imax = mod(iframe-1,nmax_history) + 1                
          max_history( imax ) = max_absol
        
          ! average over history 
          max_average = sum( max_history )
          if( iframe < nmax_history ) then
            ! history not filled yet, only average over available entries
            max_average = max_average / iframe
          else
            ! average over all history entries
            max_average = max_average / nmax_history        
          endif

          print *,'maximum amplitude over averaged last snapshots = ',max_average

          ! scales field values up to match average 
          if( abs(max_absol) > TINYVAL) &
            field_display = field_display * max_average / max_absol 

          ! thresholds positive & negative maximum values          
          where( field_display(:) > max_average ) field_display = max_average          
          where( field_display(:) < - max_average ) field_display = -max_average
          
          ! normalizes field values
          if( NORMALIZE_VALUES ) then
            if( abs(max_average) > TINYVAL ) field_display = field_display / max_average
          endif
          
        endif

        print *
        print *,'initial number of points (with multiples) was ',npointot
        print *,'final number of points is                     ',ieoff

        !--- ****** create GMT file ******

        ! create file name and open file
        if(OUTPUT_BINARY) then
          if(USE_COMPONENT == 1) then
           write(outputname,"('bin_movie_',i6.6,'.d')") it
          elseif(USE_COMPONENT == 2) then
           write(outputname,"('bin_movie_',i6.6,'.N')") it
          elseif(USE_COMPONENT == 3) then
           write(outputname,"('bin_movie_',i6.6,'.E')") it
          endif
          open(unit=11,file='OUTPUT_FILES/'//trim(outputname),status='unknown',form='unformatted')
          if(iframe == 1) open(unit=12,file='OUTPUT_FILES/bin_movie.xy',status='unknown',form='unformatted')
        else
          if(USE_COMPONENT == 1) then
           write(outputname,"('ascii_movie_',i6.6,'.d')") it
          elseif(USE_COMPONENT == 2) then
           write(outputname,"('ascii_movie_',i6.6,'.N')") it
          elseif(USE_COMPONENT == 3) then
           write(outputname,"('ascii_movie_',i6.6,'.E')") it
          endif
          open(unit=11,file='OUTPUT_FILES/'//trim(outputname),status='unknown')
          if(iframe == 1) open(unit=12,file='OUTPUT_FILES/ascii_movie.xy',status='unknown')
        endif
        ! clear number of elements kept
        ispec = 0

        ! read points for all the slices
        print *,'Writing output',outputname
        do iproc = 0,NPROCTOT-1

          ! reset point number
          ipoin = 0

          do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA
            ispec = ispec + 1
            if(MOVIE_COARSE) then
              ielm = ispec - 1
            else
              ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
            endif
            
            do j = 1,NGLLY-NIT
              do i = 1,NGLLX-NIT
                if(MOVIE_COARSE) then
                  ieoff = ielm + 1
                else
                  ieoff = (ielm+(i-1)+(j-1)*(NGLLX-1))+1
                endif
                
                ! point position
                if(iframe == 1) then
                  ! gets cartesian coordinates
                  xcoord = sngl(xp(ieoff))
                  ycoord = sngl(yp(ieoff))
                  zcoord = sngl(zp(ieoff))
                
                  ! location latitude/longitude (with geocentric colatitude theta )
                  call xyz_2_rthetaphi(xcoord,ycoord,zcoord,rval,thetaval,phival)
                  
                  ! converts the geocentric colatitude to a geographic colatitude
                  if(.not. ASSUME_PERFECT_SPHERE) then
                    thetaval = PI/2.0d0 - &
                      datan(1.006760466d0*dcos(dble(thetaval))/dmax1(TINYVAL,dble(sin(thetaval))))
                  endif
                  
                  ! gets geographic latitude and longitude in degrees
                  lat = sngl(90.d0 - thetaval*180.0/PI)
                  long = sngl(phival*180.0/PI)
                  if(long > 180.0) long = long-360.0
                endif
                
                ! displacement
                disp = sngl(field_display(ieoff))
                
                ! writes displacement and latitude/longitude to corresponding files
                if(OUTPUT_BINARY) then
                  write(11) disp
                  if(iframe == 1) write(12) long,lat
                else
                  write(11,*) disp
                  if(iframe == 1) write(12,*) long,lat
                endif
                
              enddo !i
            enddo !j
          enddo !ispecloc
        enddo !iproc
        close(11)
        if(iframe == 1) close(12)


! end of loop and test on all the time steps for all the movie images
     endif
  enddo

  print *,'done creating movie'
  print *,'GMT ascii files are stored in ascii_movie_*.{xy,d,E,N}'
  print *,'binary files are stored in bin_movie_*.{xy,d,E,N}'
  
  end program create_movie_GMT_global


