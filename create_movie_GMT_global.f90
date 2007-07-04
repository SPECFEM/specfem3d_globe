!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
!---  create a movie of radial component of surface displacement in GMT format
!

  program create_movie_GMT_global

  implicit none

  include "constants.h"

  integer i,j,it
  integer it1,it2
  integer ispec

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,displn
  real(kind=CUSTOM_REAL) xcoord,ycoord,zcoord,rval,thetaval,phival,lat,long
  real(kind=CUSTOM_REAL) RRval,rhoval
  real(kind=CUSTOM_REAL) displx,disply,displz
  real(kind=CUSTOM_REAL) normal_x,normal_y,normal_z
  real(kind=CUSTOM_REAL) thetahat_x,thetahat_y,thetahat_z
  real(kind=CUSTOM_REAL) phihat_x,phihat_y
  double precision min_field_current,max_field_current,max_absol
  real disp
  integer nframes,iframe,USE_COMPONENT

  character(len=150) outputname

  integer iproc,ipoin

! for sorting routine
  integer npointot,ilocnum,ielm,ieoff,ispecloc
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
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE,REFERENCE_1D_MODEL

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

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
  logical :: CASE_3D

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Recombining all movie frames to create a movie'
  print *

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
         R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
         TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
         ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
         MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
         PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
         ATTENUATION,REFERENCE_1D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
         NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC, &
         NSPEC2D_XI, &
         NSPEC2D_ETA, &
         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
         NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
         NGLOB, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
	 OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
	 ROTATE_SEISMOGRAMS_RT)

  if(.not. MOVIE_SURFACE) stop 'movie frames were not saved by the solver'

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *

  ilocnum = NGLLX*NGLLY*NEX_PER_PROC_XI*NEX_PER_PROC_ETA

  print *
  print *,'Allocating arrays for reading data of size ',ilocnum*NPROCTOT,'=',6*ilocnum*NPROCTOT*CUSTOM_REAL/1000000,'MB'
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



  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  print *,'enter first time step of movie (e.g. 1)'
  read(5,*) it1

  print *,'enter last time step of movie (e.g. ',NSTEP,')'
  read(5,*) it2

  print *,'enter component (e.g. 1=Z, 2=N, 3=E)'
  read(5,*) USE_COMPONENT

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

! maximum theoretical number of points at the surface
  npointot = NCHUNKS * NEX_XI * NEX_ETA * (NGLLX-1) * (NGLLY-1)

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


!--- ****** read data saved by solver ******

! --------------------------------------

  iframe = 0

! loop on all the time steps in the range entered
  do it = it1,it2

     ! check if time step corresponds to a movie frame
     if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

        iframe = iframe + 1

        ! read all the elements from the same file
        write(outputname,"('moviedata',i6.6)") it
        open(unit=IOUT,file=outputname,status='old',form='unformatted')

        print *
        print *,'reading snapshot time step ',it,' out of ',NSTEP,' file ',outputname
        print *

        read(IOUT) store_val_x
        read(IOUT) store_val_y
        read(IOUT) store_val_z
        read(IOUT) store_val_ux
        read(IOUT) store_val_uy
        read(IOUT) store_val_uz
        close(IOUT)
        print *, 'finished reading ',outputname
        ! clear number of elements kept
        ispec = 0

        ! read points for all the slices
        print *,'Converting to geo-coordinates'
        do iproc = 0,NPROCTOT-1

           ! reset point number
           ipoin = 0

           do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA

              do j = 1,NGLLY
                 do i = 1,NGLLX

                    ipoin = ipoin + 1

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
                 enddo
              enddo

              ispec = ispec + 1
              ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
              do j = 1,NGLLY-1
                 do i = 1,NGLLX-1
                    ieoff = (ielm+(i-1)+(j-1)*(NGLLX-1))
                    xp(ieoff) = dble(x(i,j))
                    yp(ieoff) = dble(y(i,j))
                    zp(ieoff) = dble(z(i,j))
                    field_display(ieoff) = dble(displn(i,j))
                 enddo
              enddo

           enddo

        enddo

        ! compute min and max of data value to normalize
        min_field_current = minval(field_display(:))
        max_field_current = maxval(field_display(:))

        ! print minimum and maximum amplitude in current snapshot
        print *
        print *,'minimum amplitude in current snapshot = ',min_field_current
        print *,'maximum amplitude in current snapshot = ',max_field_current
        print *

        ! make sure range is always symmetric and center is in zero
        ! this assumption works only for fields that can be negative
        ! would not work for norm of vector for instance
        ! (we would lose half of the color palette if no negative values)
        max_absol = max(abs(min_field_current),abs(max_field_current))
        min_field_current = - max_absol
        max_field_current = + max_absol


        print *
        print *,'initial number of points (with multiples) was ',npointot

        !--- ****** create GMT file ******

        ! create file name and open file
        if(USE_COMPONENT == 1) then
           write(outputname,"('ascii_movie_',i6.6,'.d')") it
        elseif(USE_COMPONENT == 2) then
           write(outputname,"('ascii_movie_',i6.6,'.N')") it
        elseif(USE_COMPONENT == 3) then
           write(outputname,"('ascii_movie_',i6.6,'.E')") it
        endif

        ! clear number of elements kept
        ispec = 0

        ! read points for all the slices
        print *,'Writing output'
        do iproc = 0,NPROCTOT-1

           ! reset point number
           ipoin = 0

           do ispecloc = 1,NEX_PER_PROC_XI*NEX_PER_PROC_ETA
              ispec = ispec + 1
              ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
              do j = 1,NGLLY-1
                 do i = 1,NGLLX-1
                    ieoff = (ielm+(i-1)+(j-1)*(NGLLX-1))
                    xcoord = sngl(xp(ieoff))
                    ycoord = sngl(yp(ieoff))
                    zcoord = sngl(zp(ieoff))
                    call xyz_2_rthetaphi(xcoord,ycoord,zcoord,rval,thetaval,phival)
                    lat = sngl((PI/2.0-thetaval)*180.0/PI)
                    long = sngl(phival*180.0/PI)
                    disp = sngl(field_display(ieoff))
                   if(long > 180.0) long = long-360.0
                    write(11,*) long, lat, disp
                 enddo
              enddo
           enddo
        enddo


! end of loop and test on all the time steps for all the movie images
     endif
  enddo

  print *,'done creating movie'
  print *,'GMT files are stored in ascii_files_*.{xy,d,E,N}'

end program create_movie_GMT_global


