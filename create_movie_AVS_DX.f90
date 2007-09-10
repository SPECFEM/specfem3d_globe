!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, October 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
!---  create a movie of radial component of surface displacement
!---  in AVS or OpenDX format
!

  subroutine create_movie_AVS_DX(iformat,it1,it2,NEX_XI,NEX_ETA,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NCHUNKS,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  implicit none

  include "constants.h"

! threshold in percent of the maximum below which we cut the amplitude
  real(kind=CUSTOM_REAL), parameter :: THRESHOLD = 1._CUSTOM_REAL / 100._CUSTOM_REAL

! flag to apply non linear scaling to normalized norm of displacement
  logical, parameter :: NONLINEAR_SCALING = .true.

! coefficient of power law used for non linear scaling
  real(kind=CUSTOM_REAL), parameter :: POWER_SCALING = 0.30_CUSTOM_REAL

! flag to cut amplitude below a certain threshold
  logical, parameter :: APPLY_THRESHOLD = .false.

  integer i,j,it
  integer it1,it2
  integer nspectot_AVS_max
  integer ispec
  integer ibool_number,ibool_number1,ibool_number2,ibool_number3,ibool_number4
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: x,y,z,displn
  real(kind=CUSTOM_REAL) xcoord,ycoord,zcoord,rval,thetaval,phival,lat,long
  real(kind=CUSTOM_REAL) displx,disply,displz
  real(kind=CUSTOM_REAL) normal_x,normal_y,normal_z
  double precision min_field_current,max_field_current,max_absol
  logical USE_OPENDX,UNIQUE_FILE,USE_GMT,USE_AVS
  integer iformat,nframes,iframe

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

! parameters read from file or deduced from parameters read from file
  integer NEX_XI,NEX_ETA
  integer NSTEP,NTSTEP_BETWEEN_FRAMES,NCHUNKS
  integer NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  character(len=150) OUTPUT_FILES

! --------------------------------------

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
  nspectot_AVS_max = NCHUNKS * NEX_XI * NEX_ETA * (NGLLX-1) * (NGLLY-1)

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

! loop on all the time steps in the range entered
  do it = it1,it2

! check if time step corresponds to a movie frame
  if(mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

  iframe = iframe + 1

  print *
  print *,'reading snapshot time step ',it,' out of ',NSTEP
  print *

! read all the elements from the same file
  write(outputname,"('/moviedata',i6.6)") it
  open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='old',action='read',form='unformatted')
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
      ielm = (NGLLX-1)*(NGLLY-1)*(ispec-1)
      do j = 1,NGLLY-1
        do i = 1,NGLLX-1
          ieoff = NGNOD2D_AVS_DX*(ielm+(i-1)+(j-1)*(NGLLX-1))
          do ilocnum = 1,NGNOD2D_AVS_DX
            if(ilocnum == 1) then
              xp(ieoff+ilocnum) = dble(x(i,j))
              yp(ieoff+ilocnum) = dble(y(i,j))
              zp(ieoff+ilocnum) = dble(z(i,j))
              field_display(ieoff+ilocnum) = dble(displn(i,j))
            elseif(ilocnum == 2) then
              xp(ieoff+ilocnum) = dble(x(i+1,j))
              yp(ieoff+ilocnum) = dble(y(i+1,j))
              zp(ieoff+ilocnum) = dble(z(i+1,j))
              field_display(ieoff+ilocnum) = dble(displn(i+1,j))
            elseif(ilocnum == 3) then
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
  if(APPLY_THRESHOLD) &
    where(abs(field_display(:)) <= THRESHOLD) field_display = 0.

! apply non linear scaling to normalized field if needed
  if(NONLINEAR_SCALING) then
    where(field_display(:) >= 0.)
      field_display = field_display ** POWER_SCALING
    elsewhere
      field_display = - abs(field_display) ** POWER_SCALING
    endwhere
  endif

! map back to [0,1]
  field_display(:) = (field_display(:) + 1.) / 2.

! map field to [0:255] for AVS color scale
  field_display(:) = 255. * field_display(:)

! copy coordinate arrays since the sorting routine does not preserve them
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
    write(outputname,"('/DX_movie_',i6.6,'.dx')") it
    open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',nglob,' data follows'
  else if(USE_AVS) then
    if(UNIQUE_FILE .and. iframe == 1) then
      open(unit=11,file=trim(OUTPUT_FILES)//'/AVS_movie_all.inp',status='unknown')
      write(11,*) nframes
      write(11,*) 'data'
      write(11,"('step',i1,' image',i1)") 1,1
      write(11,*) nglob,' ',nspectot_AVS_max
    else if(.not. UNIQUE_FILE) then
      write(outputname,"('/AVS_movie_',i6.6,'.inp')") it
      open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
      write(11,*) nglob,' ',nspectot_AVS_max,' 1 0 0'
    endif
  else if(USE_GMT) then
    write(outputname,"('/gmt_movie_',i6.6,'.xyz')") it
    open(unit=11,file=trim(OUTPUT_FILES)//outputname,status='unknown')
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
          write(11,"(i6,1x,f10.7,1x,f10.7,1x,f10.7)") ireorder(ibool_number), &
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

  if(.not. UNIQUE_FILE) close(11)

! end of loop and test on all the time steps for all the movie images
  endif
  enddo

  if(UNIQUE_FILE) close(11)

  print *
  print *,'done creating movie'
  print *
  if(USE_OPENDX) print *,'DX files are stored in ', trim(OUTPUT_FILES), '/DX_*.dx'
  if(USE_AVS) print *,'AVS files are stored in ', trim(OUTPUT_FILES), '/AVS_*.inp'
  if(USE_GMT) print *,'GMT files are stored in ', trim(OUTPUT_FILES), '/gmt_*.xyz'
  print *

  end subroutine create_movie_AVS_DX

!
!=====================================================================
!

  subroutine read_params_and_create_movie

!
! This routine is called by the Pyrized version.
!

  implicit none

  integer it1,it2
  integer iformat

! parameters read from parameter file
  integer NEX_XI,NEX_ETA
  integer NSTEP,NTSTEP_BETWEEN_FRAMES,NCHUNKS
  integer NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  logical MOVIE_SURFACE

  integer, external :: err_occurred

  call read_AVS_DX_parameters(NEX_XI,NEX_ETA, &
           NSTEP,NTSTEP_BETWEEN_FRAMES, &
           NCHUNKS,MOVIE_SURFACE, &
           NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

! read additional parameters for making movies
  call read_value_integer(iformat, 'format')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(it1, 'beginning')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'
  call read_value_integer(it2, 'end')
  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! run the main program
  call create_movie_AVS_DX(iformat,it1,it2, &
           NEX_XI,NEX_ETA, &
           NSTEP,NTSTEP_BETWEEN_FRAMES, &
           NCHUNKS, &
           NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  end subroutine read_params_and_create_movie

! ------------------------------------------------------------------

  subroutine read_AVS_DX_parameters(NEX_XI,NEX_ETA, &
          NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NCHUNKS,MOVIE_SURFACE, &
          NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               NGLOB

  character(len=150) LOCAL_PATH,MODEL

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
         ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
         NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC, &
         NSPEC2D_XI, &
         NSPEC2D_ETA, &
         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
         NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
         NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube)

  end subroutine read_AVS_DX_parameters

! ------------------------------------------------------------------

  subroutine get_global_AVS(nspec,xp,yp,zp,iglob,loc,ifseg,nglob,npointot)

! this routine MUST be in double precision to avoid sensitivity
! to roundoff errors in the coordinates of the points

! leave sorting subroutines in same source file to allow for inlining

  implicit none

  include "constants.h"

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
