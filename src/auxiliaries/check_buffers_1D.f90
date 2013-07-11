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
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE,REFERENCE_1D_MODEL, &
          THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY, &
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
          RECOMPUTE_STRAIN_DO_NOT_STORE,EXACT_MASS_MATRIX_FOR_ROTATION

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

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! processor identification
  character(len=150) prname,prname_other

  integer :: NGLOB1D_RADIAL_MAX
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NGLOB1D_RADIAL_CORNER
  integer, dimension(NB_SQUARE_CORNERS) :: NGLOB1D_RADIAL_SPEC_THIS
  integer, dimension(NB_SQUARE_CORNERS) :: NGLOB1D_RADIAL_SPEC_OTHER

  logical :: PARTIAL_PHYS_DISPERSION_ONLY,UNDO_ATTENUATION
  integer :: NT_DUMP_ATTENUATION

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers along xi and eta inside each chunk'
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
          ATT1,ATT2,ATT3,ATT4,ATT5,RECOMPUTE_STRAIN_DO_NOT_STORE,EXACT_MASS_MATRIX_FOR_ROTATION)

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

  NGLOB1D_RADIAL_CORNER(iregion_code,:) = NGLOB1D_RADIAL(iregion_code)
  NGLOB1D_RADIAL_MAX = NGLOB1D_RADIAL(iregion_code)
  if (iregion_code == IREGION_OUTER_CORE .and. (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA)) then
    NGLOB1D_RADIAL_MAX = NGLOB1D_RADIAL_MAX + maxval(DIFF_NSPEC1D_RADIAL(:,:))*(NGLLZ-1)
  endif

! dynamic memory allocation for arrays
  allocate(iboolleft(NGLOB1D_RADIAL_MAX+1))
  allocate(iboolright(NGLOB1D_RADIAL_MAX+1))
  allocate(xleft(NGLOB1D_RADIAL_MAX+1))
  allocate(yleft(NGLOB1D_RADIAL_MAX+1))
  allocate(zleft(NGLOB1D_RADIAL_MAX+1))
  allocate(xright(NGLOB1D_RADIAL_MAX+1))
  allocate(yright(NGLOB1D_RADIAL_MAX+1))
  allocate(zright(NGLOB1D_RADIAL_MAX+1))

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

  NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_CORNER(iregion_code,:)
  if (iregion_code==IREGION_OUTER_CORE) then
    if (CUT_SUPERBRICK_XI) then
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_xi,2) == 0) then
          if (mod(iproc_eta,2) == 0) then
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
          endif
        else
          if (mod(iproc_eta,2) == 0) then
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
          endif
        endif
      else
        if (mod(iproc_xi,2) == 0) then
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    else
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_eta,2) == 0) then
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    endif
  endif
  NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_CORNER(iregion_code,:)
  if (iregion_code==IREGION_OUTER_CORE) then
    if (CUT_SUPERBRICK_XI) then
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_xi+1,2) == 0) then
          if (mod(iproc_eta,2) == 0) then
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
          endif
        else
          if (mod(iproc_eta,2) == 0) then
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
          endif
        endif
      else
        if (mod(iproc_xi+1,2) == 0) then
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    else
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_eta,2) == 0) then
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    endif
  endif
! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! read 1D addressing buffers for copy between slices along xi with MPI

  if(icorners == 1) then
! read ibool1D_rightxi_lefteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_lefteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_rightxi_lefteta.txt',status='old',action='read')
  else if(icorners == 2) then
! read ibool1D_rightxi_righteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_righteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt',status='old',action='read')
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
  if(icorners == 1) then
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_THIS(2)) stop 'incorrect iboolright read'
  else
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_THIS(3)) stop 'incorrect iboolright read'
  endif
  close(34)

  if(icorners == 1) then
! read ibool1D_leftxi_lefteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_lefteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_leftxi_lefteta.txt',status='old',action='read')
  else if(icorners == 2) then
! read ibool1D_leftxi_righteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_righteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_leftxi_righteta.txt',status='old',action='read')
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
  if(icorners == 1) then
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_OTHER(1)) stop 'incorrect iboolleft read'
  else
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_OTHER(4)) stop 'incorrect iboolleft read'
  endif
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

  NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_CORNER(iregion_code,:)
  if (iregion_code==IREGION_OUTER_CORE) then
    if (CUT_SUPERBRICK_XI) then
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_xi,2) == 0) then
          if (mod(iproc_eta,2) == 0) then
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
          endif
        else
          if (mod(iproc_eta,2) == 0) then
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
          endif
        endif
      else
        if (mod(iproc_xi,2) == 0) then
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    else
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_eta,2) == 0) then
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_THIS(:) = NGLOB1D_RADIAL_SPEC_THIS(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    endif
  endif
  NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_CORNER(iregion_code,:)
  if (iregion_code==IREGION_OUTER_CORE) then
    if (CUT_SUPERBRICK_XI) then
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_xi,2) == 0) then
          if (mod(iproc_eta+1,2) == 0) then
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
          endif
        else
          if (mod(iproc_eta+1,2) == 0) then
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
          else
            NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
          endif
        endif
      else
        if (mod(iproc_xi,2) == 0) then
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    else
      if (CUT_SUPERBRICK_ETA) then
        if (mod(iproc_eta+1,2) == 0) then
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NGLOB1D_RADIAL_SPEC_OTHER(:) = NGLOB1D_RADIAL_SPEC_OTHER(:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      endif
    endif
  endif
! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! read 1D addressing buffers for copy between slices along xi with MPI

  if(icorners == 1) then
! read ibool1D_leftxi_righteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_righteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_leftxi_righteta.txt',status='old',action='read')
  else if(icorners == 2) then
! read ibool1D_rightxi_righteta of this slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_righteta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt',status='old',action='read')
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

  if(icorners == 1) then
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_THIS(4)) stop 'incorrect iboolright read'
  else
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_THIS(3)) stop 'incorrect iboolright read'
  endif
  close(34)

  if(icorners == 1) then
! read ibool1D_leftxi_lefteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_leftxi_lefteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_leftxi_lefteta.txt',status='old',action='read')
  else if(icorners == 2) then
! read ibool1D_rightxi_lefteta of other slice
  write(*,*) 'reading MPI 1D buffer ibool1D_rightxi_lefteta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'ibool1D_rightxi_lefteta.txt',status='old',action='read')
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

  if(icorners == 1) then
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_OTHER(1)) stop 'incorrect iboolleft read'
  else
    if(npoin1D /= NGLOB1D_RADIAL_SPEC_OTHER(2)) stop 'incorrect iboolleft read'
  endif
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

