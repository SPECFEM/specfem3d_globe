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

! code to check that all the 2D buffers between chunk faces are okay

  program check_buffers_faces_chunks

  implicit none

  include "constants.h"

  integer imsg

  integer npoin2D_sender,npoin2D_receiver
  integer iboolsend,iboolreceive,ipoin2D
  integer iregion_code,iproc

! number of faces between chunks
  integer NUM_FACES,NUMMSGS_FACES

! number of message types
  integer NUM_MSG_TYPES

  double precision xsend,ysend,zsend
  double precision xreceive,yreceive,zreceive
  double precision diff

  integer NPROC_ONE_DIRECTION

! communication pattern for faces between chunks
  integer, dimension(:), allocatable :: iprocfrom_faces,iprocto_faces,imsg_type

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
          NPROC_ETA,NPROC_XI,NSEIS,NSTEP,NSOURCES,NMOVIE,NER_ICB_BOTTOMDBL, &
          NER_TOPDBL_CMB,ITAFF_TIME_STEPS,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC,HDUR_MIN_MOVIES, &
          ANGULAR_WIDTH_XI_DEG,ANGULAR_WIDTH_ETA_DEG

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,THREE_D,TOPOGRAPHY, &
          ATTENUATION,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME, ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCT,SAVE_AVS_DX_MESH_FILES

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
  print *,'Check all MPI buffers between chunk faces'
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
        NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, ATTENUATION_3D, &
        RECEIVERS_CAN_BE_BURIED,ANGULAR_WIDTH_XI_DEG,ANGULAR_WIDTH_ETA_DEG, &
        SAVE_AVS_DX_MESH_FILES,ITAFF_TIME_STEPS,PRINT_SOURCE_TIME_FUNCT, &
        NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN)

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

! number of corners and faces shared between chunks and number of message types
  if(NCHUNKS == 1 .or. NCHUNKS == 2) then
    NUM_FACES = 1
    NUM_MSG_TYPES = 1
  else if(NCHUNKS == 3) then
    NUM_FACES = 1
    NUM_MSG_TYPES = 3
  else if(NCHUNKS == 6) then
    NUM_FACES = 4
    NUM_MSG_TYPES = 3
  else
    stop 'number of chunks must be either 1, 2, 3 or 6'
  endif

! if more than one chunk then same number of processors in each direction
  NPROC_ONE_DIRECTION = NPROC_XI

! total number of messages corresponding to these common faces
  NUMMSGS_FACES = NPROC_ONE_DIRECTION*NUM_FACES*NUM_MSG_TYPES

  if(NCHUNKS == 1) stop 'only one chunk, nothing to check'

  print *,'There are ',NUMMSGS_FACES,' messages to assemble all the faces'
  print *

! allocate array for messages for faces
  allocate(iprocfrom_faces(NUMMSGS_FACES))
  allocate(iprocto_faces(NUMMSGS_FACES))
  allocate(imsg_type(NUMMSGS_FACES))

! file with the list of processors for each message for faces
  open(unit=IIN,file='OUTPUT_FILES/list_messages_faces.txt',status='old')
  do imsg = 1,NUMMSGS_FACES
  read(IIN,*) imsg_type(imsg),iprocfrom_faces(imsg),iprocto_faces(imsg)
  if      (iprocfrom_faces(imsg) < 0 &
        .or. iprocto_faces(imsg) < 0 &
        .or. iprocfrom_faces(imsg) > NPROCTOT-1 &
        .or. iprocto_faces(imsg) > NPROCTOT-1) &
    stop 'incorrect chunk faces numbering'
  if (imsg_type(imsg) < 1 .or. imsg_type(imsg) > 3) &
    stop 'incorrect message type labeling'
  enddo
  close(IIN)

! loop over all the regions of the mesh
  do iregion_code = 1,MAX_NUM_REGIONS

  print *
  print *,' ********* checking region ',iregion_code,' *********'
  print *

! loop on all the messages between faces
  do imsg = 1,NUMMSGS_FACES

  print *
  print *,'Checking message ',imsg,' out of ',NUMMSGS_FACES

! read 2-D buffer for the sender and the receiver
  write(filename,200) imsg
  iproc = iprocfrom_faces(imsg)
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)
  open(unit=34,file=prname(1:len_trim(prname))//filename,status='old')

  write(filename,210) imsg
  iproc = iprocto_faces(imsg)
  call create_serial_name_database(prname,iproc,iregion_code,LOCAL_PATH,NPROCTOT)
  open(unit=35,file=prname(1:len_trim(prname))//filename,status='old')

  200 format('buffer_faces_chunks_sender_msg',i4.4,'.txt')
  210 format('buffer_faces_chunks_receiver_msg',i4.4,'.txt')

  write(*,*) 'reading MPI 2D buffer for sender'
  read(34,*) npoin2D_sender
  read(35,*) npoin2D_receiver

! check that number of points is the same in both buffers
  if(npoin2D_sender /= npoin2D_receiver) &
        stop 'different number of points in the two buffers'

  print *,'this message contains ',npoin2D_sender,' points'

! check all the points based upon their coordinates
  do ipoin2D = 1,npoin2D_sender
  read(34,*) iboolsend,xsend,ysend,zsend
  read(35,*) iboolreceive,xreceive,yreceive,zreceive

  diff = dmax1(dabs(xsend-xreceive),dabs(ysend-yreceive),dabs(zsend-zreceive))
  if(diff > 0.0000001d0) then
    print *,'different : ',ipoin2D,iboolsend,iboolreceive,diff
    print *,'xsend,xreceive = ',xsend,xreceive
    print *,'ysend,yreceive = ',ysend,yreceive
    print *,'zsend,zreceive = ',zsend,zreceive
    stop 'error: different'
  endif

  enddo

  enddo

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_faces_chunks

