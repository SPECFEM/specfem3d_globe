!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! code to check that all the internal MPI buffers are okay
! inside any given chunk, along both xi and eta
! we compare the coordinates of the points in the buffers

  program check_buffers_2D

  implicit none

  include "constants.h"

  integer ithisproc,iotherproc

  integer ipoin

  integer npoin2d_xi_save,npoin2d_xi_mesher,npoin2d_xi
  integer npoin2d_eta_save,npoin2d_eta_mesher,npoin2d_eta

! for addressing of the slices
  integer ichunk,iproc_xi,iproc_eta,iproc
  integer iproc_read,iregion_code
  integer, dimension(:,:,:), allocatable :: addressing

  double precision diff

! 2-D addressing and buffers for summation between slices
  integer, dimension(:), allocatable :: iboolleft_xi,iboolright_xi, &
    iboolleft_eta,iboolright_eta

! coordinates of the points to compare
  double precision, dimension(:), allocatable :: xleft_xi,yleft_xi,zleft_xi, &
     xright_xi,yright_xi,zright_xi,xleft_eta,yleft_eta,zleft_eta, &
     xright_eta,yright_eta,zright_eta

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE,REFERENCE_1D_MODEL

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

! now this is for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
               nglob_AB,nglob_AC,nglob_BC

! processor identification
  character(len=150) prname,prname_other

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers along xi and eta inside each chunk'
  print *

! read the parameter file
  call read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, HDUR_MOVIE, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,REFERENCE_1D_MODEL,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD)

! compute other parameters based upon values read
  call compute_parameters(NER_CRUST,NER_220_MOHO,NER_400_220, &
      NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
      NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
      NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB, &
      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NSPEC_AB,NSPEC_AC,NSPEC_BC, &
      NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
      NGLOB_AB,NGLOB_AC,NGLOB_BC,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NCHUNKS,INCLUDE_CENTRAL_CUBE)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *,'There are ',NCHUNKS,' chunks'
  print *,'There are ',NPROC_XI,' slices along xi in each chunk'
  print *,'There are ',NPROC_ETA,' slices along eta in each chunk'
  print *

! dynamic memory allocation for arrays
  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))

! open file with global slice number addressing
  print *,'reading slice addressing'
  open(unit=34,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read')
  do iproc = 0,NPROCTOT-1
      read(34,*) iproc_read,ichunk,iproc_xi,iproc_eta
      if(iproc_read /= iproc) stop 'incorrect slice number read'
      addressing(ichunk,iproc_xi,iproc_eta) = iproc
  enddo
  close(34)

! loop over all the regions of the mesh
  do iregion_code = 1,MAX_NUM_REGIONS

  print *
  print *,' ********* checking region ',iregion_code,' *********'
  print *

! dynamic memory allocation for arrays
  allocate(iboolleft_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(iboolright_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(iboolleft_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(iboolright_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(xleft_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(yleft_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(zleft_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(xright_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(yright_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(zright_xi(NPOIN2DMAX_XMIN_XMAX(iregion_code)))
  allocate(xleft_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(yleft_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(zleft_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(xright_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(yright_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))
  allocate(zright_eta(NPOIN2DMAX_YMIN_YMAX(iregion_code)))

! loop on the chunks
  do ichunk = 1,NCHUNKS

  print *
  print *,'Checking xi in chunk ',ichunk
  print *

! double loop on NPROC_XI and NPROC_ETA
  do iproc_eta=0,NPROC_ETA-1

  print *,'checking row ',iproc_eta

  do iproc_xi=0,NPROC_XI-2

  print *,'checking slice ixi = ',iproc_xi,' in that row'

  ithisproc = addressing(ichunk,iproc_xi,iproc_eta)
  iotherproc = addressing(ichunk,iproc_xi+1,iproc_eta)

! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolright_xi of this slice
  write(*,*) 'reading MPI buffer iboolright_xi slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 360  continue
  read(34,*) iboolright_xi(npoin2D_xi), &
              xright_xi(npoin2D_xi),yright_xi(npoin2D_xi),zright_xi(npoin2D_xi)
  if(iboolright_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 360
  endif
  npoin2D_xi = npoin2D_xi - 1
  write(*,*) 'found ',npoin2D_xi,' points in iboolright_xi slice ',ithisproc
  read(34,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX(iregion_code) .or. npoin2D_xi /= npoin2D_xi_mesher) then
      stop 'incorrect iboolright_xi read'
  endif
  close(34)

! save to compare to other side
  npoin2D_xi_save = npoin2D_xi

! read iboolleft_xi of other slice
  write(*,*) 'reading MPI buffer iboolleft_xi slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'iboolleft_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 350  continue
  read(34,*) iboolleft_xi(npoin2D_xi), &
              xleft_xi(npoin2D_xi),yleft_xi(npoin2D_xi),zleft_xi(npoin2D_xi)
  if(iboolleft_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 350
  endif
  npoin2D_xi = npoin2D_xi - 1
  write(*,*) 'found ',npoin2D_xi,' points in iboolleft_xi slice ',iotherproc
  read(34,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX(iregion_code) .or. npoin2D_xi /= npoin2D_xi_mesher) then
      stop 'incorrect iboolleft_xi read'
  endif
  close(34)

  if(npoin2D_xi_save == npoin2D_xi) then
      print *,'okay, same size for both buffers'
  else
      stop 'wrong buffer size'
  endif

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin2D_xi
      diff = dmax1(dabs(xleft_xi(ipoin)-xright_xi(ipoin)), &
       dabs(yleft_xi(ipoin)-yright_xi(ipoin)),dabs(zleft_xi(ipoin)-zright_xi(ipoin)))
      if(diff > 0.0000001d0) then
            print *,'different: ',ipoin,iboolleft_xi(ipoin),iboolright_xi(ipoin),diff
            stop 'error: different'
      endif
  enddo

  enddo
  enddo
  enddo


! loop on the chunks
  do ichunk = 1,NCHUNKS

  print *
  print *,'Checking eta in chunk ',ichunk
  print *

! double loop on NPROC_XI and NPROC_ETA
  do iproc_xi=0,NPROC_XI-1

  print *,'checking row ',iproc_xi

  do iproc_eta=0,NPROC_ETA-2

  print *,'checking slice ieta = ',iproc_eta,' in that row'

  ithisproc = addressing(ichunk,iproc_xi,iproc_eta)
  iotherproc = addressing(ichunk,iproc_xi,iproc_eta+1)

! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolright_eta of this slice
  write(*,*) 'reading MPI buffer iboolright_eta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 460  continue
  read(34,*) iboolright_eta(npoin2D_eta), &
              xright_eta(npoin2D_eta),yright_eta(npoin2D_eta),zright_eta(npoin2D_eta)
  if(iboolright_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 460
  endif
  npoin2D_eta = npoin2D_eta - 1
  write(*,*) 'found ',npoin2D_eta,' points in iboolright_eta slice ',ithisproc
  read(34,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX(iregion_code) .or. npoin2D_eta /= npoin2D_eta_mesher) then
      stop 'incorrect iboolright_eta read'
  endif
  close(34)

! save to compare to other side
  npoin2D_eta_save = npoin2D_eta

! read iboolleft_eta of other slice
  write(*,*) 'reading MPI buffer iboolleft_eta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'iboolleft_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 450  continue
  read(34,*) iboolleft_eta(npoin2D_eta), &
              xleft_eta(npoin2D_eta),yleft_eta(npoin2D_eta),zleft_eta(npoin2D_eta)
  if(iboolleft_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 450
  endif
  npoin2D_eta = npoin2D_eta - 1
  write(*,*) 'found ',npoin2D_eta,' points in iboolleft_eta slice ',iotherproc
  read(34,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX(iregion_code) .or. npoin2D_eta /= npoin2D_eta_mesher) then
      stop 'incorrect iboolleft_eta read'
  endif
  close(34)

  if(npoin2D_eta_save == npoin2D_eta) then
      print *,'okay, same size for both buffers'
  else
      stop 'wrong buffer size'
  endif

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin2D_eta
      diff = dmax1(dabs(xleft_eta(ipoin)-xright_eta(ipoin)), &
       dabs(yleft_eta(ipoin)-yright_eta(ipoin)),dabs(zleft_eta(ipoin)-zright_eta(ipoin)))
      if(diff > 0.0000001d0) then
            print *,'different: ',ipoin,iboolleft_eta(ipoin),iboolright_eta(ipoin),diff
            stop 'error: different'
      endif
  enddo

  enddo
  enddo
  enddo

! deallocate arrays
  deallocate(iboolleft_xi)
  deallocate(iboolright_xi)
  deallocate(iboolleft_eta)
  deallocate(iboolright_eta)
  deallocate(xleft_xi)
  deallocate(yleft_xi)
  deallocate(zleft_xi)
  deallocate(xright_xi)
  deallocate(yright_xi)
  deallocate(zright_xi)
  deallocate(xleft_eta)
  deallocate(yleft_eta)
  deallocate(zleft_eta)
  deallocate(xright_eta)
  deallocate(yright_eta)
  deallocate(zright_eta)

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_2D

