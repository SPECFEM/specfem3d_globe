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

  subroutine read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
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
                        MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,&
                        MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
                        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
                        ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
                        ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
                        MOVIE_VOLUME,MOVIE_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
                        PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
                        ATTENUATION,ATTENUATION_NEW,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
                        INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE, &
                        LOCAL_PATH,LOCAL_TMP_PATH,MODEL, &
                        SIMULATION_TYPE,SAVE_FORWARD, &
                        NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                        NSPEC,NSPEC2D_XI,NSPEC2D_ETA,NSPEC2DMAX_XMIN_XMAX, &
                        NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
                        ratio_sampling_array, ner, doubling_index,r_bottom,r_top, &
                        this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
                        OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
                        ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube, &
                        HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
                        DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
                        WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,&
                        USE_BINARY_FOR_LARGE_FILE,EMULATE_ONLY,NOISE_TOMOGRAPHY)

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
          NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read,NOISE_TOMOGRAPHY

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
          HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM, &
          MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,&
          MOVIE_SOUTH_DEG,RECORD_LENGTH_IN_MINUTES

  logical ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,&
         MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
         RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
         SAVE_MESH_FILES,ATTENUATION,ATTENUATION_NEW, &
         ABSORBING_CONDITIONS,SAVE_FORWARD, &
         OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
         ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
         SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

  character(len=150) OUTPUT_FILES,LOCAL_PATH,LOCAL_TMP_PATH,MODEL

! parameters to be computed based upon parameters above read from file
  integer NSTEP,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,REFERENCE_1D_MODEL,THREE_D_MODEL

  double precision DT,ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400, &
          R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
          RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER

  double precision MOVIE_TOP,MOVIE_BOTTOM,MOVIE_EAST,MOVIE_WEST,&
          MOVIE_NORTH,MOVIE_SOUTH

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ONE_CRUST,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
          ATTENUATION_3D,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE, &
          EMULATE_ONLY

  integer NEX_MAX

  double precision ELEMENT_WIDTH

  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
      NGLOB

  integer nblocks_xi,nblocks_eta

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                          DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval

! honor PREM Moho or not
! doing so drastically reduces the stability condition and therefore the time step
  logical :: HONOR_1D_SPHERICAL_MOHO,CASE_3D

  integer :: ifirst_region, ilast_region, iter_region, iter_layer, doubling, padding, tmp_sum, tmp_sum_xi, tmp_sum_eta
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs
  integer ::  NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
              nb_lay_sb, nspec_sb, nglob_vol, nglob_surf, nglob_edge

! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer :: last_doubling_layer, cut_doubling, nglob_int_surf_xi, nglob_int_surf_eta,nglob_ext_surf,&
              normal_doubling, nglob_center_edge, nglob_corner_edge, nglob_border_edge
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  integer :: tmp_sum_nglob2D_xi, tmp_sum_nglob2D_eta,divider,nglob_edges_h,nglob_edge_v,to_remove

  ! reads in Par_file values
  call read_parameter_file(OUTPUT_FILES, &
                          LOCAL_PATH,LOCAL_TMP_PATH,MODEL, &
                          NTSTEP_BETWEEN_OUTPUT_SEISMOS,NTSTEP_BETWEEN_READ_ADJSRC,NTSTEP_BETWEEN_FRAMES, &
                          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
                          NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
                          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP, &
                          NEX_XI_read,NEX_ETA_read,NPROC_XI_read,NPROC_ETA_read, &
                          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
                          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,&
                          HDUR_MOVIE,MOVIE_TOP_KM,MOVIE_BOTTOM_KM,RECORD_LENGTH_IN_MINUTES, &
                          MOVIE_EAST_DEG,MOVIE_WEST_DEG,MOVIE_NORTH_DEG,MOVIE_SOUTH_DEG,&
                          ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,&
                          MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
                          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
                          SAVE_MESH_FILES,ATTENUATION,ATTENUATION_NEW,ABSORBING_CONDITIONS,SAVE_FORWARD, &
                          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
                          ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER, &
                          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,NOISE_TOMOGRAPHY)

  ! converts values to radians
  MOVIE_EAST = MOVIE_EAST_DEG * DEGREES_TO_RADIANS
  MOVIE_WEST = MOVIE_WEST_DEG * DEGREES_TO_RADIANS
  MOVIE_NORTH = (90.0d0 - MOVIE_NORTH_DEG) * DEGREES_TO_RADIANS ! converting from latitude to colatitude
  MOVIE_SOUTH = (90.0d0 - MOVIE_SOUTH_DEG) * DEGREES_TO_RADIANS
  ! converts movie top/bottom depths to radii
  MOVIE_TOP = (R_EARTH_KM-MOVIE_TOP_KM)/R_EARTH_KM
  MOVIE_BOTTOM = (R_EARTH_KM-MOVIE_BOTTOM_KM)/R_EARTH_KM

  ! include central cube or not
  ! use regular cubed sphere instead of cube for large distances
  if(NCHUNKS == 6) then
    INCLUDE_CENTRAL_CUBE = .true.
    INFLATE_CENTRAL_CUBE = .false.
  else
    INCLUDE_CENTRAL_CUBE = .false.
    INFLATE_CENTRAL_CUBE = .true.
  endif

  if(.not. EMULATE_ONLY) then
    NEX_XI = NEX_XI_read
    NEX_ETA = NEX_ETA_read
    NPROC_XI = NPROC_XI_read
    NPROC_ETA = NPROC_ETA_read
  else
    ! this is used in UTILS/estimate_best_values_runs.f90 only, to estimate memory use
    NEX_ETA = NEX_XI
    NPROC_ETA = NPROC_XI
  endif

  ! turns on/off corresponding 1-D/3-D model flags
  ! and sets radius for each discontinuity and ocean density values
  call get_model_parameters(MODEL,REFERENCE_1D_MODEL,THREE_D_MODEL, &
                        ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
                        CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO, &
                        ISOTROPIC_3D_MANTLE,ONE_CRUST,TRANSVERSE_ISOTROPY, &
                        OCEANS,TOPOGRAPHY, &
                        ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400,R600,R670,R771, &
                        RTOPDDOUBLEPRIME,RCMB,RICB,RMOHO_FICTITIOUS_IN_MESHER, &
                        R80_FICTITIOUS_IN_MESHER,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS)

  ! sets time step size and number of layers
  ! right distribution is determined based upon maximum value of NEX
  NEX_MAX = max(NEX_XI,NEX_ETA)
  call get_timestep_and_layers(DT,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD, &
                          NER_CRUST,NER_80_MOHO,NER_220_80,NER_400_220,&
                          NER_600_400,NER_670_600,NER_771_670, &
                          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
                          NER_TOP_CENTRAL_CUBE_ICB,R_CENTRAL_CUBE, &
                          NEX_MAX,NCHUNKS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
                          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,&
                          ONE_CRUST,HONOR_1D_SPHERICAL_MOHO,CASE_3D,CRUSTAL, &
                          ANISOTROPIC_INNER_CORE)

  ! initial guess : compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)

!! DK DK make sure NSTEP is a multiple of NT_DUMP_ATTENUATION
  if(UNDO_ATTENUATION .and. mod(NSTEP,NT_DUMP_ATTENUATION) /= 0) NSTEP = (NSTEP/NT_DUMP_ATTENUATION + 1)*NT_DUMP_ATTENUATION

! if doing benchmark runs to measure scaling of the code for a limited number of time steps only
  if (DO_BENCHMARK_RUN_ONLY) NSTEP = NSTEP_FOR_BENCHMARK

!<YANGL
  if ( NOISE_TOMOGRAPHY /= 0 )   NSTEP = 2*NSTEP-1   ! time steps needs to be doubled, due to +/- branches
!>YANGL

  ! subsets used to save seismograms must not be larger than the whole time series,
  ! otherwise we waste memory
  if(NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! computes a default hdur_movie that creates nice looking movies.
  ! Sets HDUR_MOVIE as the minimum period the mesh can resolve
  if(HDUR_MOVIE <= TINYVAL) &
    HDUR_MOVIE = 1.2d0*max(240.d0/NEX_XI*18.d0*ANGULAR_WIDTH_XI_IN_DEGREES/90.d0, &
                           240.d0/NEX_ETA*18.d0*ANGULAR_WIDTH_ETA_IN_DEGREES/90.d0)

  ! checks parameters
  call rcp_check_parameters(NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
                        NCHUNKS,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                        ATTENUATION,ATTENUATION_NEW,ABSORBING_CONDITIONS, &
                        INCLUDE_CENTRAL_CUBE,OUTPUT_SEISMOS_SAC_ALPHANUM)

  ! check that mesh can be coarsened in depth three or four times
  CUT_SUPERBRICK_XI=.false.
  CUT_SUPERBRICK_ETA=.false.

  if (SUPPRESS_CRUSTAL_MESH .and. .not. ADD_4TH_DOUBLING) then
    if(mod(NEX_XI,8) /= 0) stop 'NEX_XI must be a multiple of 8'
    if(mod(NEX_ETA,8) /= 0) stop 'NEX_ETA must be a multiple of 8'
    if(mod(NEX_XI/4,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 4*NPROC_XI'
    if(mod(NEX_ETA/4,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 4*NPROC_ETA'
    if(mod(NEX_XI/8,NPROC_XI) /=0) CUT_SUPERBRICK_XI = .true.
    if(mod(NEX_ETA/8,NPROC_ETA) /=0) CUT_SUPERBRICK_ETA = .true.
  else if (SUPPRESS_CRUSTAL_MESH .or. .not. ADD_4TH_DOUBLING) then
    if(mod(NEX_XI,16) /= 0) stop 'NEX_XI must be a multiple of 16'
    if(mod(NEX_ETA,16) /= 0) stop 'NEX_ETA must be a multiple of 16'
    if(mod(NEX_XI/8,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 8*NPROC_XI'
    if(mod(NEX_ETA/8,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 8*NPROC_ETA'
    if(mod(NEX_XI/16,NPROC_XI) /=0) CUT_SUPERBRICK_XI = .true.
    if(mod(NEX_ETA/16,NPROC_ETA) /=0) CUT_SUPERBRICK_ETA = .true.
  else
    if(mod(NEX_XI,32) /= 0) stop 'NEX_XI must be a multiple of 32'
    if(mod(NEX_ETA,32) /= 0) stop 'NEX_ETA must be a multiple of 32'
    if(mod(NEX_XI/16,NPROC_XI) /= 0) stop 'NEX_XI must be a multiple of 16*NPROC_XI'
    if(mod(NEX_ETA/16,NPROC_ETA) /= 0) stop 'NEX_ETA must be a multiple of 16*NPROC_ETA'
    if(mod(NEX_XI/32,NPROC_XI) /=0) CUT_SUPERBRICK_XI = .true.
    if(mod(NEX_ETA/32,NPROC_ETA) /=0) CUT_SUPERBRICK_ETA = .true.
  endif

  ELEMENT_WIDTH = ANGULAR_WIDTH_XI_IN_DEGREES/dble(NEX_MAX) * DEGREES_TO_RADIANS

!
!--- compute additional parameters
!

  ! number of elements horizontally in each slice (i.e. per processor)
  ! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI / NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA / NPROC_ETA

  ! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

  ! total number of processors in the full Earth composed of the six chunks
  NPROCTOT = NCHUNKS * NPROC

  !  definition of general mesh parameters
  call define_all_layers(NER_CRUST,NER_80_MOHO,NER_220_80,&
                        NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
                        NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
                        NER_TOP_CENTRAL_CUBE_ICB,&
                        RMIDDLE_CRUST,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
                        R_CENTRAL_CUBE,RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER,&
                        ONE_CRUST,ner,ratio_sampling_array,&
                        NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer, &
                        r_bottom,r_top,this_region_has_a_doubling,&
                        ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,&
                        elem_doubling_bottom_outer_core,&
                        DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                        DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval,&
                        doubling_index,rmins,rmaxs)

  ! calculates number of elements (NSPEC)
  call count_elements(NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NPROC,&
                        NEX_PER_PROC_ETA,ratio_divide_central_cube,&
                        NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        NSPEC1D_RADIAL, &
                        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                        ner,ratio_sampling_array,this_region_has_a_doubling, &
                        ifirst_region,ilast_region,iter_region,iter_layer,&
                        doubling,tmp_sum,tmp_sum_xi,tmp_sum_eta, &
                        NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
                        nb_lay_sb, nspec_sb, nglob_surf, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, INCLUDE_CENTRAL_CUBE, &
                        last_doubling_layer, &
                        DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
                        tmp_sum_nglob2D_xi, tmp_sum_nglob2D_eta,divider,nglob_edges_h,&
                        nglob_edge_v,to_remove)

  ! calculates number of points (NGLOB)
  call count_points(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube,&
                        NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
                        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB,&
                        nblocks_xi,nblocks_eta,ner,ratio_sampling_array,&
                        this_region_has_a_doubling,&
                        ifirst_region, ilast_region, iter_region, iter_layer, &
                        doubling, padding, tmp_sum, &
                        INCLUDE_CENTRAL_CUBE,NER_TOP_CENTRAL_CUBE_ICB,NEX_XI, &
                        NUMBER_OF_MESH_LAYERS,layer_offset, &
                        nb_lay_sb, nglob_vol, nglob_surf, nglob_edge, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                        last_doubling_layer, cut_doubling, nglob_int_surf_xi, nglob_int_surf_eta,nglob_ext_surf,&
                        normal_doubling, nglob_center_edge, nglob_corner_edge, nglob_border_edge)

  end subroutine read_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rcp_check_parameters(NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
                        NCHUNKS,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
                        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
                        ATTENUATION,ATTENUATION_NEW,ABSORBING_CONDITIONS, &
                        INCLUDE_CENTRAL_CUBE,OUTPUT_SEISMOS_SAC_ALPHANUM)

  implicit none

  include "constants.h"

  integer  NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA,NCHUNKS,NTSTEP_BETWEEN_OUTPUT_SEISMOS

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES

  logical ATTENUATION,ATTENUATION_NEW,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,OUTPUT_SEISMOS_SAC_ALPHANUM

! checks parameters

  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    stop 'NCHUNKS must be either 1, 2, 3 or 6'

  ! this MUST be 90 degrees for two chunks or more to match geometrically
  if(NCHUNKS > 1 .and. abs(ANGULAR_WIDTH_XI_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_XI_IN_DEGREES must be 90 for more than one chunk'

  ! this can be any value in the case of two chunks
  if(NCHUNKS > 2 .and. abs(ANGULAR_WIDTH_ETA_IN_DEGREES - 90.d0) > 0.00000001d0) &
    stop 'ANGULAR_WIDTH_ETA_IN_DEGREES must be 90 for more than two chunks'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 6) &
    stop 'cannot have absorbing conditions in the full Earth'

  if(ABSORBING_CONDITIONS .and. NCHUNKS == 3) &
    stop 'absorbing conditions not supported for three chunks yet'

  if(ATTENUATION_NEW .and. .not. ATTENUATION) &
    stop 'need ATTENUATION to use ATTENUATION_NEW'

  if (OUTPUT_SEISMOS_SAC_ALPHANUM .and. (mod(NTSTEP_BETWEEN_OUTPUT_SEISMOS,5)/=0)) &
    stop 'if OUTPUT_SEISMOS_SAC_ALPHANUM = .true. then NTSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of 5, check the Par_file'

  ! check that reals are either 4 or 8 bytes
  if(CUSTOM_REAL /= SIZE_REAL .and. CUSTOM_REAL /= SIZE_DOUBLE) &
    stop 'wrong size of CUSTOM_REAL for reals'

  ! check that the parameter file is correct
  if(NGNOD /= 27) &
    stop 'number of control nodes must be 27'
  if(NGNOD == 27 .and. NGNOD2D /= 9) &
    stop 'elements with 27 points should have NGNOD2D = 9'

  ! for the number of standard linear solids for attenuation
  if(N_SLS /= 3) &
    stop 'number of SLS must be 3'

  ! check number of slices in each direction
  if(NCHUNKS < 1) &
    stop 'must have at least one chunk'
  if(NPROC_XI < 1) &
    stop 'NPROC_XI must be at least 1'
  if(NPROC_ETA < 1) &
    stop 'NPROC_ETA must be at least 1'

  ! check number of chunks
  if(NCHUNKS /= 1 .and. NCHUNKS /= 2 .and. NCHUNKS /= 3 .and. NCHUNKS /= 6) &
    stop 'only one, two, three or six chunks can be meshed'

  ! check that the central cube can be included
  if(INCLUDE_CENTRAL_CUBE .and. NCHUNKS /= 6) &
    stop 'need six chunks to include central cube'

  ! check that sphere can be cut into slices without getting negative Jacobian
  if(NCHUNKS == 6 ) then
    if(NEX_XI < 48) &
      stop 'NEX_XI must be greater than 48 to cut the sphere into slices with positive Jacobian'
    if(NEX_ETA < 48) &
      stop 'NEX_ETA must be greater than 48 to cut the sphere into slices with positive Jacobian'
  endif

  ! check that topology is correct if more than two chunks
  if(NCHUNKS > 2 .and. NEX_XI /= NEX_ETA) &
    stop 'must have NEX_XI = NEX_ETA for more than two chunks'

  if(NCHUNKS > 2 .and. NPROC_XI /= NPROC_ETA) &
    stop 'must have NPROC_XI = NPROC_ETA for more than two chunks'

  ! small meshes useful for testing, also for GPU version
  if(NCHUNKS > 1 .and. (NPROC_XI == 1 .or. NPROC_ETA == 1) ) then
    if( NUMFACES_SHARED < 4 ) &
      stop 'NPROC_XI,NPROC_ETA==1: please set in constants.h NUMFACES_SHARED and NUMCORNERS_SHARED equal to 4 and recompile'
    if( NUMCORNERS_SHARED < 4 ) &
      stop 'NPROC_XI,NPROC_ETA==1: please set in constants.h NUMFACES_SHARED and NUMCORNERS_SHARED equal to 4 and recompile'
  endif

  end subroutine rcp_check_parameters

