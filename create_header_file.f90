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

! create file OUTPUT_FILES/values_from_mesher.h based upon DATA/Par_file
! in order to compile the solver with the right array sizes

  subroutine create_header_file

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE,REFERENCE_1D_MODEL

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION,CASE_3D, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD

  character(len=150) LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               nglob
  integer :: ilayer, ner_without_doubling

  integer ispec_aniso,NUMBER_OF_MESH_LAYERS,static_size,dynamic_size
  character(len=150) HEADER_FILE

! variables added for memory size computation
  integer :: NSPECMAX_CRUST_MANTLE_ATTENUAT,NSPEC_INNER_CORE_ATTENUATION,NSPECMAX_ISO_MANTLE, &
            NSPECMAX_TISO_MANTLE,NSPECMAX_ANISO_MANTLE,NSPECMAX_ANISO_IC,NSPECMAX_OUTER_CORE_ROTATION

! ************** PROGRAM STARTS HERE **************

  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')
  print *
  print *,'creating file ', trim(HEADER_FILE), ' to compile solver with correct values'

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
         NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
         ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D)

! count the total number of sources in the CMTSOLUTION file
  call count_number_of_sources(NSOURCES)

! generate the elements in all the regions of the mesh
  ispec_aniso = 0

  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif

! count anisotropic elements
do ilayer = 1, NUMBER_OF_MESH_LAYERS
    if (doubling_index(ilayer) == IFLAG_220_80 .or. doubling_index(ilayer) == IFLAG_80_MOHO) then
        ner_without_doubling = ner(ilayer)
        if(this_region_has_a_doubling(ilayer)) then
            ner_without_doubling = ner_without_doubling - 2
            ispec_aniso = ispec_aniso + &
            (NSPEC_DOUBLING_SUPERBRICK*(NEX_PER_PROC_XI/ratio_sampling_array(ilayer)/2)* &
            (NEX_PER_PROC_ETA/ratio_sampling_array(ilayer)/2))
        endif
        ispec_aniso = ispec_aniso + &
        ((NEX_PER_PROC_XI/ratio_sampling_array(ilayer))*(NEX_PER_PROC_ETA/ratio_sampling_array(ilayer))*ner_without_doubling)
    endif
enddo

  print *
  print *,'found ',ispec_aniso,' tranversely isotropic elements in the mantle'

! create include file for the solver
  call save_header_file(NSPEC, &
        nglob,NEX_XI,NEX_ETA,ispec_aniso,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION,ATTENUATION_3D, &
        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
        INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP,&
        SIMULATION_TYPE,SAVE_FORWARD,MOVIE_VOLUME)

  print *
  print *,'edit file OUTPUT_FILES/values_from_mesher.h to see some statistics about the mesh'
  print *
  print *,'on NEC SX and Earth Simulator, make sure "loopcnt=" parameter'
! use fused loops on the ES
  print *,'in Makefile is greater than max vector length = ',nglob(IREGION_CRUST_MANTLE)*NDIM

! DM memory size evaluation
    static_size = 0
    dynamic_size = 0

  if(ATTENUATION) then
    if(ATTENUATION_3D) then
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
! factor_scale_crust_mantle,one_minus_sum_beta_crust_mantle,factor_common_crust_mantle
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*3*CUSTOM_REAL
! factor_scale_inner_core, one_minus_sum_beta_inner_core, factor_common_inner_core
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*3*SIZE_DOUBLE
! factor_scale_crust_mantle_dble, omsb_crust_mantle_dble, factor_common_crust_mantle_dble
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*3*SIZE_DOUBLE
! factor_scale_inner_core_dble, omsb_inner_core_dble, factor_common_inner_core_dble
    endif
    NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_ATTENUATION = NSPEC(IREGION_INNER_CORE)
  else
    NSPECMAX_CRUST_MANTLE_ATTENUAT = 1
    NSPEC_INNER_CORE_ATTENUATION = 1
  endif
  if(ANISOTROPIC_3D_MANTLE) then
    NSPECMAX_ISO_MANTLE = 1
    NSPECMAX_TISO_MANTLE = 1
    NSPECMAX_ANISO_MANTLE = NSPEC(IREGION_CRUST_MANTLE)
  else
    NSPECMAX_ISO_MANTLE = NSPEC(IREGION_CRUST_MANTLE)
    if(TRANSVERSE_ISOTROPY) then
      NSPECMAX_TISO_MANTLE = ispec_aniso
    else
      NSPECMAX_TISO_MANTLE = 1
    endif
    NSPECMAX_ANISO_MANTLE = 1
  endif
  if(ANISOTROPIC_INNER_CORE) then
    NSPECMAX_ANISO_IC = NSPEC(IREGION_INNER_CORE)
  else
    NSPECMAX_ANISO_IC = 1
  endif
  if(ROTATION) then
    NSPECMAX_OUTER_CORE_ROTATION = NSPEC(IREGION_OUTER_CORE)
  else
    NSPECMAX_OUTER_CORE_ROTATION = 1
  endif

! size of static arrays
! ---------------------

    static_size = static_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPECMAX_CRUST_MANTLE_ATTENUAT*CUSTOM_REAL   !R_memory_crust_mantle
    static_size = static_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION*CUSTOM_REAL     !R_memory_inner_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*12*CUSTOM_REAL
! rho_vp_crust_mantle,rho_vs_crust_mantle,xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle
! etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,
! gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,jacobian_crust_mantle
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*SIZE_REAL
! ibool_crust_mantle
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*13*CUSTOM_REAL
! vp_outer_core,xix_outer_core,xiy_outer_core,xiz_outer_core,
! etax_outer_core,etay_outer_core,etaz_outer_core,
! gammax_outer_core,gammay_outer_core,gammaz_outer_core,jacobian_outer_core
! rhostore_outer_core,kappavstore_outer_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*SIZE_REAL
! ibool_outer_core
    static_size = static_size + nglob(IREGION_CRUST_MANTLE)*2*SIZE_REAL
! updated_dof_ocean_load, idoubling_crust_mantle
    static_size = static_size + nglob(IREGION_CRUST_MANTLE)*5*CUSTOM_REAL
! xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
! rmass_crust_mantle,rmass_ocean_load
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_ISO_MANTLE*3*CUSTOM_REAL
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_TISO_MANTLE*3*CUSTOM_REAL
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_ANISO_MANTLE*21*CUSTOM_REAL
    static_size = static_size + NDIM*nglob(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
! displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
    static_size = static_size + nglob(IREGION_OUTER_CORE)*7*CUSTOM_REAL
! xstore_outer_core, ystore_outer_core, zstore_outer_core, rmass_outer_core, displ_outer_core, veloc_outer_core, accel_outer_core

    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*SIZE_REAL
! ibool_inner_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*13*CUSTOM_REAL
! xix_inner_core,xiy_inner_core,xiz_inner_core,
! etax_inner_core,etay_inner_core,etaz_inner_core,
! gammax_inner_core,gammay_inner_core,gammaz_inner_core,jacobian_inner_core,
! rhostore_inner_core, kappavstore_inner_core,muvstore_inner_core
    static_size = static_size + nglob(IREGION_INNER_CORE)*4*CUSTOM_REAL
! xstore_inner_core,ystore_inner_core,zstore_inner_core,rmass_inner_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_ANISO_IC*5*CUSTOM_REAL
! c11store_inner_core,c33store_inner_core,c12store_inner_core,
! c13store_inner_core,c44store_inner_core
    static_size = static_size + NDIM*nglob(IREGION_INNER_CORE)*3*CUSTOM_REAL
! displ_inner_core,veloc_inner_core,accel_inner_core

    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_OUTER_CORE_ROTATION*2*CUSTOM_REAL
! A_array_rotation,B_array_rotation

! size of dynamic arrays
! ----------------------

    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*2 !dummy_rho,dummy_vstore
  if (SIMULATION_TYPE > 1) then
    dynamic_size = dynamic_size + NDIM*nglob(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + nglob(IREGION_OUTER_CORE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NDIM*nglob(IREGION_INNER_CORE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*2*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*2*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_OUTER_CORE_ROTATION*2*CUSTOM_REAL
  endif
  if (ATTENUATION .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    dynamic_size = dynamic_size + 7*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*CUSTOM_REAL
    dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*CUSTOM_REAL
! epsilondev_crust_mantle(5), eps_trace_over_3_crust_mantle, eps_trace_over_3_inner_core
    if (SIMULATION_TYPE > 1) then
      dynamic_size = dynamic_size + 7*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*CUSTOM_REAL
      dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*CUSTOM_REAL
    endif
  else
    dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPECMAX_CRUST_MANTLE_ATTENUAT*CUSTOM_REAL
    dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION*CUSTOM_REAL
  endif
  if(ATTENUATION .and. SIMULATION_TYPE > 1) then
    dynamic_size = dynamic_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPECMAX_CRUST_MANTLE_ATTENUAT*CUSTOM_REAL
    dynamic_size = dynamic_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION*CUSTOM_REAL
  endif

    print *
    print *, 'memory needed for the solver : '
    print *
    print *, 'approximate size of static arrays per slice : ',static_size/(1024**2),' MB'
    print *, 'approximate size of static arrays for all slices : ',((static_size/(1024**2))*NPROCTOT)/1024.d0,' GB'
    print *
    print *, 'approximate size of dynamic arrays per slice : ',dynamic_size/(1024**2),' MB'
    print *, 'approximate size of dynamic arrays for all slices : ',((dynamic_size/(1024**2))*NPROCTOT)/1024.d0,' GB'
    print *
    print *, 'approximate total size of arrays per slice : ',(dynamic_size+static_size)/(1024**2),' MB'
    print *, 'approximate total size of arrays for all slices : ', &
                    (((dynamic_size+static_size)/(1024**2))*NPROCTOT)/1024.d0,' GB'
    print *
    print *,'done'
    print *

! DM memory size evaluation

  end subroutine create_header_file

