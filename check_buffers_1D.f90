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

! code to check that all the internal MPI 1D buffers are okay
! inside any given chunk, along both xi and eta
! we compare the coordinates of the points in the buffers

  program check_buffers_1D

  implicit none

  include "constants.h"

  integer ithisproc,iotherproc
  integer ipoin

  double precision diff

  integer npoin1D_mesher,npoin1D

! for addressing of the slices
  integer ichunk,iproc_xi,iproc_eta,iproc,icorners,iregion_code
  integer iproc_read
  integer, dimension(:,:,:), allocatable :: addressing

! 1D addressing for copy of edges between slices
! we add one to the size of the array for the final flag
  integer, dimension(:), allocatable :: iboolleft,iboolright
  double precision, dimension(:), allocatable :: xleft,yleft,zleft,xright,yright,zright

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
          NPROC_ETA,NPROC_XI,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_AVS_DX_MESH_FILES,ATTENUATION,IASPEI, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

! this is for all the regions
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
          NER_TOP_CENTRAL_CUBE_ICB,NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
          NPROC_ETA,NPROC_XI,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_AVS_DX_MESH_FILES, &
          ATTENUATION,IASPEI,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL)

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
  open(unit=34,file='OUTPUT_FILES/addressing.txt',status='old')
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
  allocate(iboolleft(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(iboolright(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(xleft(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(yleft(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(zleft(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(xright(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(yright(NPOIN1D_RADIAL(iregion_code)+1))
  allocate(zright(NPOIN1D_RADIAL(iregion_code)+1))

! ********************************************************
! ***************  check along xi
! ********************************************************

! loop for both corners for 1D buffers
  do icorners=1,2

  print *
  print *,'Checking for xi in set of corners # ',icorners
  print *

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
  call create_serial_name_database(prname,ithisproc,iregion_code,LOCAL_PATH,NPROCTOT)
  call create_serial_name_database(prname_other,iotherproc,iregion_code,LOCAL_PATH,NPROCTOT)

! read 1D addressing buffers for copy between slices along xi with MPI

  if(icorners == 1) then
! read ibool1D_rightxi_lefteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_lefteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_rightxi_lefteta.txt',status='old')
  else if(icorners == 2) then
! read ibool1D_rightxi_righteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_righteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt',status='old')
  else
      stop 'incorrect corner number'
  endif

  npoin1D = 1
 360  continue
  read(34,*) iboolright(npoin1D),xright(npoin1D),yright(npoin1D),zright(npoin1D)
  if(iboolright(npoin1D) > 0) then
      npoin1D = npoin1D + 1
      goto 360
  endif
  npoin1D = npoin1D - 1
  write(*,*) 'found ',npoin1D,' points in iboolright slice ',ithisproc
  read(34,*) npoin1D_mesher
  if(npoin1D /= NPOIN1D_RADIAL(iregion_code)) stop 'incorrect iboolright read'
  close(34)

  if(icorners == 1) then
! read ibool1D_leftxi_lefteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_lefteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_leftxi_lefteta.txt',status='old')
  else if(icorners == 2) then
! read ibool1D_leftxi_righteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_righteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_leftxi_righteta.txt',status='old')
  else
      stop 'incorrect corner number'
  endif

  npoin1D = 1
 350  continue
  read(34,*) iboolleft(npoin1D),xleft(npoin1D),yleft(npoin1D),zleft(npoin1D)
  if(iboolleft(npoin1D) > 0) then
      npoin1D = npoin1D + 1
      goto 350
  endif
  npoin1D = npoin1D - 1
  write(*,*) 'found ',npoin1D,' points in iboolleft slice ',iotherproc
  read(34,*) npoin1D_mesher
  if(npoin1D /= NPOIN1D_RADIAL(iregion_code)) stop 'incorrect iboolleft read'
  close(34)

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin1D
      diff = dmax1(dabs(xleft(ipoin)-xright(ipoin)), &
       dabs(yleft(ipoin)-yright(ipoin)),dabs(zleft(ipoin)-zright(ipoin)))
      if(diff > 0.0000001d0) then
            print *,'different: ',ipoin,iboolleft(ipoin),iboolright(ipoin),diff
            stop 'error: different'
      endif
  enddo

  enddo
  enddo
  enddo

  enddo


! ********************************************************
! ***************  check along eta
! ********************************************************

! added loop for both corners for 1D buffers
  do icorners=1,2

  print *
  print *,'Checking for eta in set of corners # ',icorners
  print *

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
  call create_serial_name_database(prname,ithisproc,iregion_code,LOCAL_PATH,NPROCTOT)
  call create_serial_name_database(prname_other,iotherproc,iregion_code,LOCAL_PATH,NPROCTOT)

! read 1D addressing buffers for copy between slices along xi with MPI

  if(icorners == 1) then
! read ibool1D_leftxi_righteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_righteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_leftxi_righteta.txt',status='old')
  else if(icorners == 2) then
! read ibool1D_rightxi_righteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_righteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt',status='old')
  else
      stop 'incorrect corner number'
  endif

  npoin1D = 1
 460  continue
  read(34,*) iboolright(npoin1D),xright(npoin1D),yright(npoin1D),zright(npoin1D)
  if(iboolright(npoin1D) > 0) then
      npoin1D = npoin1D + 1
      goto 460
  endif
  npoin1D = npoin1D - 1
  write(*,*) 'found ',npoin1D,' points in iboolright slice ',ithisproc
  read(34,*) npoin1D_mesher
  if(npoin1D /= NPOIN1D_RADIAL(iregion_code)) stop 'incorrect iboolright read'
  close(34)

  if(icorners == 1) then
! read ibool1D_leftxi_lefteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_lefteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_leftxi_lefteta.txt',status='old')
  else if(icorners == 2) then
! read ibool1D_rightxi_lefteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_lefteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_rightxi_lefteta.txt',status='old')
  else
      stop 'incorrect corner number'
  endif

  npoin1D = 1
 450  continue
  read(34,*) iboolleft(npoin1D),xleft(npoin1D),yleft(npoin1D),zleft(npoin1D)
  if(iboolleft(npoin1D) > 0) then
      npoin1D = npoin1D + 1
      goto 450
  endif
  npoin1D = npoin1D - 1
  write(*,*) 'found ',npoin1D,' points in iboolleft slice ',iotherproc
  read(34,*) npoin1D_mesher
  if(npoin1D /= NPOIN1D_RADIAL(iregion_code)) stop 'incorrect iboolleft read'
  close(34)

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin1D
      diff = dmax1(dabs(xleft(ipoin)-xright(ipoin)), &
       dabs(yleft(ipoin)-yright(ipoin)),dabs(zleft(ipoin)-zright(ipoin)))
      if(diff > 0.0000001d0) then
            print *,'different: ',ipoin,iboolleft(ipoin),iboolright(ipoin),diff
            stop 'error: different'
      endif
  enddo

  enddo
  enddo
  enddo

  enddo

! deallocate arrays
  deallocate(iboolleft)
  deallocate(iboolright)
  deallocate(xleft)
  deallocate(yleft)
  deallocate(zleft)
  deallocate(xright)
  deallocate(yright)
  deallocate(zright)

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_1D

