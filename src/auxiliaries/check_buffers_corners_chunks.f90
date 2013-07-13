!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

! code to check that all the 1D buffers between chunk corners are okay

  program check_buffers_corners_chunks

  implicit none

  include "constants.h"

  integer imsg
  integer ipoin1D
  integer iboolmaster,iboolworker1,iboolworker2
  integer npoin1D_master,npoin1D_worker1,npoin1D_worker2
  integer iregion_code,iproc

! number of corners between chunks
  integer NCORNERSCHUNKS

  double precision xmaster,ymaster,zmaster
  double precision xworker1,yworker1,zworker1
  double precision xworker2,yworker2,zworker2
  double precision diff1,diff2

! communication pattern for corners between chunks
  integer, dimension(:), allocatable :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY, &
          ATT1,ATT2,ATT3,ATT4,ATT5

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,RMOHO_FICTITIOUS_IN_MESHER, &
          RATIO_BY_WHICH_TO_INCREASE_IT

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,SAVE_REGULAR_KL, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX, &
          EXACT_MASS_MATRIX_FOR_ROTATION

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

! this is for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               nglob

  character(len=150) filename,prname

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  logical :: PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION
  integer :: NT_DUMP_ATTENUATION

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers between chunk corners'
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
         INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE, &
         LOCAL_PATH,MODEL, &
         SIMULATION_TYPE,SAVE_FORWARD, &
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
         ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
         DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
         WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,.false.,NOISE_TOMOGRAPHY,&
         SAVE_REGULAR_KL,PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION,NT_DUMP_ATTENUATION, &
          USE_LDDRK,INCREASE_CFL_FOR_LDDRK,ANISOTROPIC_KL,SAVE_TRANSVERSE_KL_ONLY,APPROXIMATE_HESS_KL, &
          USE_FULL_TISO_MANTLE,SAVE_SOURCE_MASK,GPU_MODE,ADIOS_ENABLED,ADIOS_FOR_FORWARD_ARRAYS, &
          ADIOS_FOR_MPI_ARRAYS,ADIOS_FOR_ARRAYS_SOLVER,ADIOS_FOR_AVS_DX,RATIO_BY_WHICH_TO_INCREASE_IT, &
          ATT1,ATT2,ATT3,ATT4,ATT5,EXACT_MASS_MATRIX_FOR_ROTATION)

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *,'There are ',NCHUNKS,' chunks'
  print *,'There are ',NPROC_XI,' slices along xi in each chunk'
  print *,'There are ',NPROC_ETA,' slices along eta in each chunk'
  print *

! number of corners shared between chunks
  if(NCHUNKS == 1 .or. NCHUNKS == 2 .or. NCHUNKS == 3) then
    NCORNERSCHUNKS = 1
  else if(NCHUNKS == 6) then
    NCORNERSCHUNKS = 8
  else
    stop 'number of chunks must be either 1, 2, 3 or 6'
  endif

  if(NCHUNKS == 1) stop 'only one chunk, nothing to check'

  print *,'There are ',NCORNERSCHUNKS,' messages to assemble all the corners'
  print *

! allocate array for messages for corners
  allocate(iproc_master_corners(NCORNERSCHUNKS))
  allocate(iproc_worker1_corners(NCORNERSCHUNKS))
  allocate(iproc_worker2_corners(NCORNERSCHUNKS))

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! file with the list of processors for each message for corners
  open(unit=IIN,file=trim(OUTPUT_FILES)//'/list_messages_corners.txt',status='old',action='read')
  do imsg = 1,NCORNERSCHUNKS
  read(IIN,*) iproc_master_corners(imsg),iproc_worker1_corners(imsg), &
                          iproc_worker2_corners(imsg)
  if    (iproc_master_corners(imsg) < 0 &
    .or. iproc_worker1_corners(imsg) < 0 &
    .or. iproc_worker2_corners(imsg) < 0 &
    .or. iproc_master_corners(imsg) > NPROCTOT-1 &
    .or. iproc_worker1_corners(imsg) > NPROCTOT-1 &
    .or. iproc_worker2_corners(imsg) > NPROCTOT-1) &
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

! master
  write(filename,"('buffer_corners_chunks_master_msg',i6.6,'.txt')") imsg
  iproc = iproc_master_corners(imsg)
  call create_serial_name_database(prname,iproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  open(unit=34,file=prname(1:len_trim(prname))//filename,status='old',action='read')

! first worker
  write(filename,"('buffer_corners_chunks_worker1_msg',i6.6,'.txt')") imsg
  iproc = iproc_worker1_corners(imsg)
  call create_serial_name_database(prname,iproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  open(unit=35,file=prname(1:len_trim(prname))//filename,status='old',action='read')

! second worker
! if only two chunks then there is no second worker
  if(NCHUNKS /= 2) then
    write(filename,"('buffer_corners_chunks_worker2_msg',i6.6,'.txt')") imsg
    iproc = iproc_worker2_corners(imsg)
    call create_serial_name_database(prname,iproc,iregion_code, &
        LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
    open(unit=36,file=prname(1:len_trim(prname))//filename,status='old',action='read')
  endif

  write(*,*) 'reading MPI 1D buffers for 3 procs corner'

  read(34,*) npoin1D_master
  read(35,*) npoin1D_worker1
! if only two chunks then there is no second worker
  if(NCHUNKS /= 2) then
    read(36,*) npoin1D_worker2
  else
    npoin1D_worker2 = npoin1D_worker1
  endif

  if(npoin1D_master /= NGLOB1D_RADIAL(iregion_code) .or. &
     npoin1D_worker1 /= NGLOB1D_RADIAL(iregion_code) .or. &
     npoin1D_worker2 /= NGLOB1D_RADIAL(iregion_code)) then
              stop 'incorrect total number of points'
  else
    print *,'number of points is correct: ',NGLOB1D_RADIAL(iregion_code)
  endif

! check all the points based upon their coordinates
  do ipoin1D = 1, NGLOB1D_RADIAL(iregion_code)

  read(34,*) iboolmaster,xmaster,ymaster,zmaster
  read(35,*) iboolworker1,xworker1,yworker1,zworker1
! if only two chunks then there is no second worker
  if(NCHUNKS /= 2) read(36,*) iboolworker2,xworker2,yworker2,zworker2

  diff1 = dmax1(dabs(xmaster-xworker1),dabs(ymaster-yworker1),dabs(zmaster-zworker1))
  if(diff1 > 0.0000001d0) then
    print *,'different : ',ipoin1D,iboolmaster,iboolworker1,diff1
    print *,'xmaster,xworker1 = ',xmaster,xworker1
    print *,'ymaster,yworker1 = ',ymaster,yworker1
    print *,'zmaster,zworker1 = ',zmaster,zworker1
    stop 'error: different'
  endif

! if only two chunks then there is no second worker
  if(NCHUNKS /= 2) then
    diff2 = dmax1(dabs(xmaster-xworker2),dabs(ymaster-yworker2),dabs(zmaster-zworker2))
    if(diff2 > 0.0000001d0) then
      print *,'different : ',ipoin1D,iboolmaster,iboolworker2,diff2
      print *,'xmaster,xworker2 = ',xmaster,xworker2
      print *,'ymaster,yworker2 = ',ymaster,yworker2
      print *,'zmaster,zworker2 = ',zmaster,zworker2
      stop 'error: different'
    endif
  endif

  enddo

  close(34)
  close(35)
! if only two chunks then there is no second worker
  if(NCHUNKS /= 2) close(36)

  enddo

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_corners_chunks

