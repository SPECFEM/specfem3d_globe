!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! code to check that all the 1D buffers between chunk corners are okay

  program check_buffers_corners_chunks

  implicit none

  include "constants.h"

  integer imsg
  integer ipoin1D
  integer iboolmaster,iboolslave1,iboolslave2
  integer npoin1D_master,npoin1D_slave1,npoin1D_slave2
  integer iregion_code,iproc

! number of corners between chunks
  integer NCORNERSCHUNKS

  double precision xmaster,ymaster,zmaster
  double precision xslave1,yslave1,zslave1
  double precision xslave2,yslave2,zslave2
  double precision diff1,diff2

! communication pattern for corners between chunks
  integer, dimension(:), allocatable :: iproc_master_corners,iproc_slave1_corners,iproc_slave2_corners

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
             NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
             NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
             NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
             NPROC_ETA,NPROC_XI,NSEIS,NSTEP

  double precision DT

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY, &
             GRAVITY,ONE_CRUST,ROTATION, &
             THREE_D,TOPOGRAPHY,ATTENUATION,OCEANS, &
             MOVIE_SURFACE,MOVIE_VOLUME
  integer NSOURCES,NMOVIE,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB
  double precision RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC,HDUR_MIN_MOVIES

  character(len=150) LOCAL_PATH

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

  character(len=150) filename,prname

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers between chunk corners'
  print *

! read the parameter file
  call read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
        NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
        NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB,NER_DOUBLING_OUTER_CORE, &
        NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP, &
        DT,TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,OCEANS,ELLIPTICITY, &
        GRAVITY,ONE_CRUST,ATTENUATION, &
        ROTATION,THREE_D,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        MOVIE_SURFACE,MOVIE_VOLUME,NMOVIE,HDUR_MIN_MOVIES, &
        NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC)

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
      NGLOB_AB,NGLOB_AC,NGLOB_BC,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB)

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *,'There are ',NCHUNKS,' chunks'
  print *,'There are ',NPROC_XI,' slices along xi in each chunk'
  print *,'There are ',NPROC_ETA,' slices along eta in each chunk'
  print *

  if(NCHUNKS == 1 .or. NCHUNKS == 3) then
    NCORNERSCHUNKS = 1
  else if(NCHUNKS == 6) then
    NCORNERSCHUNKS = 8
  else
    stop 'number of chunks must be either 1, 3 or 6'
  endif

  if(NCHUNKS == 1) stop 'only one chunk, nothing to check'

  print *,'There are ',NCORNERSCHUNKS,' messages to assemble all the corners'
  print *

! allocate array for messages for corners
  allocate(iproc_master_corners(NCORNERSCHUNKS))
  allocate(iproc_slave1_corners(NCORNERSCHUNKS))
  allocate(iproc_slave2_corners(NCORNERSCHUNKS))

! file with the list of processors for each message for corners
  open(unit=IIN,file='OUTPUT_FILES/list_messages_corners.txt',status='old')
  do imsg = 1,NCORNERSCHUNKS
  read(IIN,*) iproc_master_corners(imsg),iproc_slave1_corners(imsg), &
                          iproc_slave2_corners(imsg)
  if    (iproc_master_corners(imsg) < 0 &
    .or. iproc_slave1_corners(imsg) < 0 &
    .or. iproc_slave2_corners(imsg) < 0 &
    .or. iproc_master_corners(imsg) > NPROCTOT-1 &
    .or. iproc_slave1_corners(imsg) > NPROCTOT-1 &
    .or. iproc_slave2_corners(imsg) > NPROCTOT-1) &
      stop 'incorrect chunk corner numbering'
  enddo
  close(IIN)

! loop over all the regions of the mesh
  do iregion_code = 1,MAX_NUM_REGIONS

  print *
  print *,' ********* checking region ',iregion_code,' *********'
  print *

! loop on all the messages between corners
  do imsg = 1,NCORNERSCHUNKS

  print *
  print *,'Checking message ',imsg,' out of ',NCORNERSCHUNKS

! read 1-D buffers for the corners
  write(filename,200) imsg
  iproc = iproc_master_corners(imsg)
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)
  open(unit=34,file=prname(1:len_trim(prname))//filename,status='old')

  write(filename,210) imsg
  iproc = iproc_slave1_corners(imsg)
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)
  open(unit=35,file=prname(1:len_trim(prname))//filename,status='old')

  write(filename,220) imsg
  iproc = iproc_slave2_corners(imsg)
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)
  open(unit=36,file=prname(1:len_trim(prname))//filename,status='old')

  200 format('buffer_corners_chunks_master_msg',i4.4,'.txt')
  210 format('buffer_corners_chunks_slave1_msg',i4.4,'.txt')
  220 format('buffer_corners_chunks_slave2_msg',i4.4,'.txt')

  write(*,*) 'reading MPI 1D buffers for 3 procs corner'

  read(34,*) npoin1D_master
  read(35,*) npoin1D_slave1
  read(36,*) npoin1D_slave2

  if(npoin1D_master /= NPOIN1D_RADIAL(iregion_code) .or. &
     npoin1D_master /= NPOIN1D_RADIAL(iregion_code) .or. &
     npoin1D_master /= NPOIN1D_RADIAL(iregion_code)) then
              stop 'incorrect total number of points'
  else
    print *,'number of points is correct: ',NPOIN1D_RADIAL(iregion_code)
  endif

! check all the points based upon their coordinates
  do ipoin1D = 1, NPOIN1D_RADIAL(iregion_code)

  read(34,*) iboolmaster,xmaster,ymaster,zmaster
  read(35,*) iboolslave1,xslave1,yslave1,zslave1
  read(36,*) iboolslave2,xslave2,yslave2,zslave2

  diff1 = dmax1(dabs(xmaster-xslave1),dabs(ymaster-yslave1),dabs(zmaster-zslave1))
  diff2 = dmax1(dabs(xmaster-xslave2),dabs(ymaster-yslave2),dabs(zmaster-zslave2))
  if(diff1 > 0.0000001d0 .or. diff2 > 0.0000001d0) then
    print *,'different : ',ipoin1D,iboolmaster,iboolslave1,iboolslave2,diff1,diff2
    print *,'xmaster,xslave1 = ',xmaster,xslave1
    print *,'ymaster,yslave1 = ',ymaster,yslave1
    print *,'zmaster,zslave1 = ',zmaster,zslave1
    print *,'xmaster,xslave2 = ',xmaster,xslave2
    print *,'ymaster,yslave2 = ',ymaster,yslave2
    print *,'zmaster,zslave2 = ',zmaster,zslave2
    stop 'error: different'
  endif

  enddo

  close(34)
  close(35)
  close(36)

  enddo

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_corners_chunks

