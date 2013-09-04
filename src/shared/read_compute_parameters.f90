!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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


  subroutine read_compute_parameters()

  use constants,only: &
    TINYVAL,R_EARTH_KM,DEGREES_TO_RADIANS, &
    SUPPRESS_CRUSTAL_MESH,ADD_4TH_DOUBLING, &
    DO_BENCHMARK_RUN_ONLY,NSTEP_FOR_BENCHMARK, &
    IREGION_CRUST_MANTLE,IREGION_INNER_CORE, &
    NGLLX,NGLLY,NGLLZ


  use shared_parameters

  implicit none

  ! local parameters
  integer :: NEX_MAX
  double precision :: ELEMENT_WIDTH
  integer :: nblocks_xi,nblocks_eta
  ! doubling layers
  integer :: ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                          DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval
  integer :: ifirst_region, ilast_region, iter_region, iter_layer, doubling, padding, tmp_sum, tmp_sum_xi, tmp_sum_eta
  ! layers
  integer ::  NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
              nb_lay_sb, nspec_sb, nglob_vol, nglob_surf, nglob_edge
  ! for the cut doublingbrick improvement
  integer :: last_doubling_layer, cut_doubling, nglob_int_surf_xi, nglob_int_surf_eta,nglob_ext_surf,&
              normal_doubling, nglob_center_edge, nglob_corner_edge, nglob_border_edge
  integer :: tmp_sum_nglob2D_xi, tmp_sum_nglob2D_eta,divider,nglob_edges_h,nglob_edge_v,to_remove

  ! reads in Par_file values
  call read_parameter_file()

  ! count the total number of sources in the CMTSOLUTION file
  call count_number_of_sources(NSOURCES)

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
  call get_timestep_and_layers(NEX_MAX)

  ! initial guess : compute total number of time steps, rounded to next multiple of 100
  NSTEP = 100 * (int(RECORD_LENGTH_IN_MINUTES * 60.d0 / (100.d0*DT)) + 1)

  ! if doing benchmark runs to measure scaling of the code for a limited number of time steps only
  if (DO_BENCHMARK_RUN_ONLY) NSTEP = NSTEP_FOR_BENCHMARK

  ! noise simulations
  if ( NOISE_TOMOGRAPHY /= 0 ) then
    ! time steps needs to be doubled, due to +/- branches
    NSTEP = 2*NSTEP-1
  endif

!! DK DK make sure NSTEP is a multiple of NT_DUMP_ATTENUATION
  if(UNDO_ATTENUATION .and. mod(NSTEP,NT_DUMP_ATTENUATION) /= 0) then
    NSTEP = (NSTEP/NT_DUMP_ATTENUATION + 1)*NT_DUMP_ATTENUATION
  endif

  ! subsets used to save seismograms must not be larger than the whole time series,
  ! otherwise we waste memory
  if(NTSTEP_BETWEEN_OUTPUT_SEISMOS > NSTEP) NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP

  ! computes a default hdur_movie that creates nice looking movies.
  ! Sets HDUR_MOVIE as the minimum period the mesh can resolve
  if(HDUR_MOVIE <= TINYVAL) then
    HDUR_MOVIE = 1.2d0*max(240.d0/NEX_XI*18.d0*ANGULAR_WIDTH_XI_IN_DEGREES/90.d0, &
                           240.d0/NEX_ETA*18.d0*ANGULAR_WIDTH_ETA_IN_DEGREES/90.d0)
  endif

  ! checks parameters
  call rcp_check_parameters()

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

  if(ATTENUATION) then
!! DK DK July 2013: to save a huge amount of memory, when 3D attenuation is off it is sufficient to save a single point
!! DK DK July 2013: per spectral element because the Q attenuation factor is then constant per layer of the geological model
    if(ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ATT1     = NGLLX
      ATT2     = NGLLY
      ATT3     = NGLLZ
    else
      ATT1     = 1
      ATT2     = 1
      ATT3     = 1
    endif
    ATT4     = NSPEC(IREGION_CRUST_MANTLE)
    ATT5     = NSPEC(IREGION_INNER_CORE)
  else
     ATT1 = 1
     ATT2 = 1
     ATT3 = 1
     ATT4 = 1
     ATT5 = 1
  endif

  end subroutine read_compute_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rcp_check_parameters()

  use constants,only: &
    CUSTOM_REAL,SIZE_REAL,SIZE_DOUBLE,NUMFACES_SHARED,NUMCORNERS_SHARED, &
    N_SLS,NGNOD,NGNOD2D

  use shared_parameters

  implicit none

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

  if(ATTENUATION_3D .and. .not. ATTENUATION) &
    stop 'need ATTENUATION to use ATTENUATION_3D'

  if(SAVE_TRANSVERSE_KL_ONLY .and. .not. ANISOTROPIC_KL) &
    stop 'need ANISOTROPIC_KL to use SAVE_TRANSVERSE_KL_ONLY'

  if(PARTIAL_PHYS_DISPERSION_ONLY .and. UNDO_ATTENUATION) &
    stop 'cannot have both PARTIAL_PHYS_DISPERSION_ONLY and UNDO_ATTENUATION, they are mutually exclusive'

  if(UNDO_ATTENUATION .and. MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE == 4 ) &
    stop 'MOVIE_VOLUME_TYPE == 4 is not implemented for UNDO_ATTENUATION in order to save memory'

  !! DK DK this should not be difficult to fix and test, but not done yet by lack of time
  if( UNDO_ATTENUATION .and. NUMBER_OF_RUNS /= 1) &
    stop 'NUMBER_OF_RUNS should be == 1 for now when using UNDO_ATTENUATION'

  if(UNDO_ATTENUATION .and. NUMBER_OF_THIS_RUN > 1) &
    stop 'we currently do not support NUMBER_OF_THIS_RUN > 1 in the case of UNDO_ATTENUATION'

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

