!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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
!
! United States and French Government Sponsorship Acknowledged.

module constants

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

end module constants

!=====================================================================

module specfem_par

! main parameter module for specfem simulations

  use constants

  implicit none

  !-----------------------------------------------------------------
  ! attenuation parameters
  !-----------------------------------------------------------------

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval, betaval, gammaval
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: b_alphaval, b_betaval, b_gammaval

  !-----------------------------------------------------------------
  ! topography/bathymetry & oceans
  !-----------------------------------------------------------------

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! additional mass matrix for ocean load
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load

  ! flag to mask ocean-bottom degrees of freedom for ocean load
  logical, dimension(NGLOB_CRUST_MANTLE_OCEANS) :: updated_dof_ocean_load

  !-----------------------------------------------------------------
  ! ellipticity
  !-----------------------------------------------------------------

  ! for ellipticity
  integer :: nspl
  double precision :: rspl(NR),espl(NR),espl2(NR)

  !-----------------------------------------------------------------
  ! rotation
  !-----------------------------------------------------------------

  ! non-dimensionalized rotation rate of the Earth times two
  real(kind=CUSTOM_REAL) :: two_omega_earth

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION) :: &
    A_array_rotation,B_array_rotation

  !ADJOINT
  real(kind=CUSTOM_REAL) b_two_omega_earth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT) :: &
    b_A_array_rotation,b_B_array_rotation

  !-----------------------------------------------------------------
  ! gravity
  !-----------------------------------------------------------------

  ! lookup table every km for gravity
  real(kind=CUSTOM_REAL) :: minus_g_cmb,minus_g_icb
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table, &
    minus_deriv_gravity_table,density_table,d_ln_density_dr_table,minus_rho_g_over_kappa_fluid

  !-----------------------------------------------------------------
  ! time scheme
  !-----------------------------------------------------------------

  integer :: it

  ! Newmark time scheme parameters and non-dimensionalization
  double precision :: scale_t,scale_t_inv,scale_displ,scale_veloc
  real(kind=CUSTOM_REAL) :: deltat,deltatover2,deltatsqover2
  ! ADJOINT
  real(kind=CUSTOM_REAL) :: b_deltat,b_deltatover2,b_deltatsqover2

#ifdef _HANDOPT
  integer :: imodulo_NGLOB_CRUST_MANTLE,imodulo_NGLOB_CRUST_MANTLE4, &
            imodulo_NGLOB_INNER_CORE,imodulo_NGLOB_OUTER_CORE
#endif

  !-----------------------------------------------------------------
  ! sources
  !-----------------------------------------------------------------

  ! parameters for the source
  integer :: NSOURCES,nsources_local
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sourcearrays
  double precision, dimension(:,:,:) ,allocatable:: nu_source

  double precision, dimension(:), allocatable :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(:), allocatable :: xi_source,eta_source,gamma_source
  double precision, dimension(:), allocatable :: tshift_cmt,hdur,hdur_gaussian
  double precision, dimension(:), allocatable :: theta_source,phi_source
  double precision, external :: comp_source_time_function
  double precision :: t0

  !-----------------------------------------------------------------
  ! receivers
  !-----------------------------------------------------------------

  ! receiver information
  integer :: nrec,nrec_local
  integer, dimension(:), allocatable :: islice_selected_rec,ispec_selected_rec
  integer, dimension(:), allocatable :: number_receiver_global
  double precision, dimension(:), allocatable :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(:,:,:), allocatable :: nu
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur
  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable  :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: network_name
  character(len=150) :: STATIONS,rec_filename

  ! Lagrange interpolators at receivers
  double precision, dimension(:,:), allocatable :: hxir_store,hetar_store,hgammar_store

  !ADJOINT
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: adj_sourcearrays
  integer :: nrec_simulation, nadj_rec_local
  integer :: NSTEP_SUB_ADJ  ! to read input in chunks
  integer, dimension(:,:), allocatable :: iadjsrc ! to read input in chunks
  integer, dimension(:), allocatable :: iadjsrc_len,iadj_vec
  ! source frechet derivatives
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: moment_der
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sloc_der
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: stshift_der, shdur_der
  double precision, dimension(:,:), allocatable :: hpxir_store,hpetar_store,hpgammar_store
  integer :: nadj_hprec_local

  !-----------------------------------------------------------------
  ! seismograms
  !-----------------------------------------------------------------

  ! seismograms
  integer :: it_begin,it_end,nit_written
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms
  integer :: seismo_offset, seismo_current

  ! for SAC headers for seismograms
  integer :: yr_SAC,jda_SAC,ho_SAC,mi_SAC
  real :: mb_SAC
  double precision :: t_cmt_SAC,t_shift_SAC,elat_SAC,elon_SAC,depth_SAC, &
    cmt_lat_SAC,cmt_lon_SAC,cmt_depth_SAC,cmt_hdur_SAC,sec_SAC
  character(len=20) :: event_name_SAC


  !-----------------------------------------------------------------
  ! GLL points & weights
  !-----------------------------------------------------------------

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  ! product of weights for gravity term
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz


  !-----------------------------------------------------------------
  ! file parameters
  !-----------------------------------------------------------------

  ! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NTSTEP_BETWEEN_OUTPUT_SEISMOS,&
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY

  double precision DT,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
          ANGULAR_WIDTH_XI_IN_DEGREES

  logical ONE_CRUST,TOPOGRAPHY,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL
  ! process/partition name
  character(len=150) :: prname


  !-----------------------------------------------------------------
  ! mesh
  !-----------------------------------------------------------------

  ! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
    NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NGLOB1D_RADIAL, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX

  ! computed in read_compute_parameters
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  ! temporary arrays for elements on slices or edges
  logical, dimension(:),allocatable :: is_on_a_slice_edge_crust_mantle, &
    is_on_a_slice_edge_inner_core,is_on_a_slice_edge_outer_core

  !-----------------------------------------------------------------
  ! MPI partitions
  !-----------------------------------------------------------------

  ! for addressing of the slices
  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing
  integer, dimension(0:NPROCTOT_VAL-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice

  ! proc numbers for MPI
  integer :: myrank
  integer :: ichunk,iproc_xi,iproc_eta

  ! time loop timing
  double precision :: time_start,tCPU

  !-----------------------------------------------------------------
  ! assembly
  !-----------------------------------------------------------------

  ! ---- arrays to assemble between chunks
  ! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: iprocfrom_faces,iprocto_faces,imsg_type
  ! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  ! indirect addressing for each message for faces and corners of the chunks
  ! a given slice can belong to at most one corner and at most two faces
  integer :: NGLOB2DMAX_XY

  ! this for non blocking MPI

  ! buffers for send and receive between faces of the slices and the chunks
  ! we use the same buffers to assemble scalars and vectors because vectors are
  ! always three times bigger and therefore scalars can use the first part
  ! of the vector buffer in memory even if it has an additional index here
  integer :: npoin2D_max_all_CM_IC
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_faces,buffer_received_faces, &
                                                           b_buffer_send_faces,b_buffer_received_faces

  ! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL_CM) :: &
    buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
    b_buffer_send_chunkcorn_scalar,b_buffer_recv_chunkcorn_scalar

  ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC) :: &
     buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector, &
     b_buffer_send_chunkcorn_vector,b_buffer_recv_chunkcorn_vector

  ! request ids for non-blocking MPI
  integer :: request_send,request_receive
  integer, dimension(NUMFACES_SHARED) :: request_send_array,request_receive_array
  integer :: request_send_cc,request_receive_cc
  integer, dimension(NPROC_XI_VAL+4) :: request_send_array_cc,request_receive_array_cc

  integer :: b_request_send,b_request_receive
  integer, dimension(NUMFACES_SHARED) :: b_request_send_array,b_request_receive_array
  integer :: b_request_send_cc,b_request_receive_cc
  integer, dimension(NPROC_XI_VAL+4) :: b_request_send_array_cc,b_request_receive_array_cc


  logical, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool

  ! number of faces between chunks
  integer :: NUMMSGS_FACES

  ! number of corners between chunks
  integer :: NCORNERSCHUNKS

  ! number of message types
  integer :: NUM_MSG_TYPES


  ! collected MPI interfaces
  ! MPI crust/mantle mesh
  integer :: num_interfaces_crust_mantle
  integer :: max_nibool_interfaces_crust_mantle
  integer, dimension(:), allocatable :: my_neighbours_crust_mantle,nibool_interfaces_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_crust_mantle,b_buffer_recv_vector_crust_mantle

  integer, dimension(:), allocatable :: request_send_vector_crust_mantle,request_recv_vector_crust_mantle
  integer, dimension(:), allocatable :: b_request_send_vector_crust_mantle,b_request_recv_vector_crust_mantle

  ! MPI inner core mesh
  integer :: num_interfaces_inner_core
  integer :: max_nibool_interfaces_inner_core
  integer, dimension(:), allocatable :: my_neighbours_inner_core,nibool_interfaces_inner_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_inner_core,buffer_recv_vector_inner_core
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core

  integer, dimension(:), allocatable :: request_send_vector_inner_core,request_recv_vector_inner_core
  integer, dimension(:), allocatable :: b_request_send_vector_inner_core,b_request_recv_vector_inner_core

  ! MPI outer core mesh
  integer :: num_interfaces_outer_core
  integer :: max_nibool_interfaces_outer_core
  integer, dimension(:), allocatable :: my_neighbours_outer_core,nibool_interfaces_outer_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core  

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core
  
  integer, dimension(:), allocatable :: request_send_scalar_outer_core,request_recv_scalar_outer_core
  integer, dimension(:), allocatable :: b_request_send_scalar_outer_core,b_request_recv_scalar_outer_core

  !-----------------------------------------------------------------
  ! gpu
  !-----------------------------------------------------------------

  ! CUDA mesh pointer<->integer wrapper
  integer(kind=8) :: Mesh_pointer
  logical :: GPU_MODE

end module specfem_par


!=====================================================================

module specfem_par_crustmantle

! parameter module for elastic solver in crust/mantle region

  use constants
  implicit none

  ! ----------------- crust, mantle and oceans ---------------------
  ! mesh parameters
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  ! arrays for isotropic elements stored only where needed to save space
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
    rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle

  ! arrays for anisotropic elements stored only where needed to save space
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
    kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle

  ! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
    c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
    c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
    c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
    c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
    c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
    c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
    c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle

  ! local to global mapping
  !  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling_crust_mantle
  logical, dimension(NSPEC_CRUST_MANTLE) :: ispec_is_tiso_crust_mantle

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmass_crust_mantle

  ! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(ATT1,ATT2,ATT3,ATT4) :: &
    one_minus_sum_beta_crust_mantle, factor_scale_crust_mantle
  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1,ATT2,ATT3,ATT4) :: &
    factor_common_crust_mantle
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: &
!    R_memory_crust_mantle
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: &
    R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle,R_xz_crust_mantle,R_yz_crust_mantle

!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) :: &
!    epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: &
    eps_trace_over_3_crust_mantle

  ! ADJOINT
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: &
!    b_R_memory_crust_mantle
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: &
    b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle,b_R_xz_crust_mantle,b_R_yz_crust_mantle

!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
!    b_epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
    b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    b_eps_trace_over_3_crust_mantle

  ! for crust/oceans coupling
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_CM) :: ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_CM) :: ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM) :: ibelm_bottom_crust_mantle
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_CM) :: &
    jacobian2D_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_CM) :: &
    jacobian2D_top_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    normal_xmin_crust_mantle,normal_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    normal_ymin_crust_mantle,normal_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_CM) :: &
    normal_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_CM) :: &
    normal_top_crust_mantle

  ! Stacey
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY) :: &
    rho_vp_crust_mantle,rho_vs_crust_mantle
  integer :: nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
            nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle
  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_CM) :: nimin_crust_mantle,nimax_crust_mantle,nkmin_eta_crust_mantle
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_CM) :: njmin_crust_mantle,njmax_crust_mantle,nkmin_xi_crust_mantle

  ! daniel: not sure why name ...5 and dimensions with one additional?
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: absorb_xmin_crust_mantle5, &
  !  absorb_xmax_crust_mantle5, absorb_ymin_crust_mantle5, absorb_ymax_crust_mantle5
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: absorb_xmin_crust_mantle, &
    absorb_xmax_crust_mantle, absorb_ymin_crust_mantle, absorb_ymax_crust_mantle

  integer :: nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm

  integer :: reclen_xmin_crust_mantle, reclen_xmax_crust_mantle, &
            reclen_ymin_crust_mantle,reclen_ymax_crust_mantle

  ! assembly
  integer :: npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle

  ! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle

  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_crust_mantle

  ! kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    rho_kl_crust_mantle,beta_kl_crust_mantle, alpha_kl_crust_mantle, Sigma_kl_crust_mantle
  ! For anisotropic kernels (see compute_kernels.f90 for a definition of the array)
  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    cijkl_kl_crust_mantle
  ! approximate hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: hess_kl_crust_mantle

  ! Boundary Mesh and Kernels
  integer :: k_top,k_bot,iregion_code
  integer, dimension(NSPEC2D_MOHO) :: ibelm_moho_top,ibelm_moho_bot
  integer, dimension(NSPEC2D_400) :: ibelm_400_top,ibelm_400_bot
  integer, dimension(NSPEC2D_670) :: ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO) :: normal_moho
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_400) :: normal_400
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_670) :: normal_670
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_MOHO) :: moho_kl, moho_kl_top, moho_kl_bot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_400) :: d400_kl, d400_kl_top, d400_kl_bot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_670) ::  d670_kl, d670_kl_top, d670_kl_bot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_CMB) :: cmb_kl, cmb_kl_top, cmb_kl_bot

  ! NOISE_TOMOGRAPHY
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_sourcearray
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    normal_x_noise,normal_y_noise,normal_z_noise, mask_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: noise_surface_movie
  integer :: irec_master_noise

  ! inner / outer elements crust/mantle region
  integer :: num_phase_ispec_crust_mantle
  integer :: nspec_inner_crust_mantle,nspec_outer_crust_mantle
  integer, dimension(:,:), allocatable :: phase_ispec_inner_crust_mantle

end module specfem_par_crustmantle

!=====================================================================

module specfem_par_innercore

! parameter module for elastic solver in inner core region

  use constants
  implicit none

  ! ----------------- inner core ---------------------
  ! mesh parameters
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    xix_inner_core,xiy_inner_core,xiz_inner_core,&
    etax_inner_core,etay_inner_core,etaz_inner_core, &
    gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
    rhostore_inner_core, kappavstore_inner_core,muvstore_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: &
    xstore_inner_core,ystore_inner_core,zstore_inner_core

  ! arrays for inner-core anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC) :: &
    c11store_inner_core,c33store_inner_core,c12store_inner_core, &
    c13store_inner_core,c44store_inner_core

  ! local to global mapping
  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core
  ! only needed for compute_boundary_kernel() routine
  !logical, dimension(NSPEC_INNER_CORE) :: ispec_is_tiso_inner_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: rmass_inner_core

  ! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: &
    displ_inner_core,veloc_inner_core,accel_inner_core
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(ATT1,ATT2,ATT3,ATT5) :: &
    one_minus_sum_beta_inner_core, factor_scale_inner_core
  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1,ATT2,ATT3,ATT5) :: &
    factor_common_inner_core
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: &
!    R_memory_inner_core
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: &
    R_xx_inner_core,R_yy_inner_core,R_xy_inner_core,R_xz_inner_core,R_yz_inner_core

!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT) :: &
!    epsilondev_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT) :: &
    epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
    epsilondev_xz_inner_core,epsilondev_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY) :: &
    eps_trace_over_3_inner_core
  ! ADJOINT
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_AND_ATT) :: &
!    b_R_memory_inner_core
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_AND_ATT) :: &
    b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core,b_R_xz_inner_core,b_R_yz_inner_core


!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
!    b_epsilondev_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
    b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core


  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    b_eps_trace_over_3_inner_core

  ! assembly
  ! for matching with central cube in inner core
  integer, dimension(:), allocatable :: sender_from_slices_to_cube
  integer, dimension(:,:), allocatable :: ibool_central_cube
  double precision, dimension(:,:), allocatable :: buffer_slices,b_buffer_slices,buffer_slices2
  double precision, dimension(:,:,:), allocatable :: buffer_all_cube_from_slices,b_buffer_all_cube_from_slices
  integer nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,receiver_cube_from_slices

  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
            nspec2D_ymin_inner_core,nspec2D_ymax_inner_core

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_IC) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_IC) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_IC) :: ibelm_bottom_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

  integer :: npoin2D_faces_inner_core(NUMFACES_SHARED)
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_inner_core,npoin2D_eta_inner_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED) :: iboolcorner_inner_core

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_inner_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    rho_kl_inner_core,beta_kl_inner_core, alpha_kl_inner_core

  ! Boundary Mesh and Kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_ICB) :: icb_kl, icb_kl_top, icb_kl_bot
  logical :: fluid_solid_boundary

  ! inner / outer elements inner core region
  integer :: num_phase_ispec_inner_core
  integer :: nspec_inner_inner_core,nspec_outer_inner_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_inner_core


end module specfem_par_innercore

!=====================================================================

module specfem_par_outercore

! parameter module for acoustic solver in outer core region

  use constants
  implicit none

  ! ----------------- outer core ---------------------
  ! mesh parameters
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
    xix_outer_core,xiy_outer_core,xiz_outer_core,&
    etax_outer_core,etay_outer_core,etaz_outer_core, &
    gammax_outer_core,gammay_outer_core,gammaz_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    xstore_outer_core,ystore_outer_core,zstore_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
    rhostore_outer_core,kappavstore_outer_core

  ! local to global mapping
  !integer, dimension(NSPEC_OUTER_CORE) :: idoubling_outer_core
  ! only needed for compute_boundary_kernel()
  !logical, dimension(NSPEC_OUTER_CORE) :: ispec_is_tiso_outer_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: rmass_outer_core

  ! velocity potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    displ_outer_core,veloc_outer_core,accel_outer_core
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core


  ! Stacey
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_STACEY) :: vp_outer_core
  integer :: nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
            nspec2D_ymin_outer_core,nspec2D_ymax_outer_core
  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_OC) :: nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_OC) :: njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: absorb_xmin_outer_core, &
     absorb_xmax_outer_core, absorb_ymin_outer_core, absorb_ymax_outer_core, &
     absorb_zmin_outer_core
  integer :: nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc,nabs_zmin_oc

  integer :: reclen_xmin_outer_core, reclen_xmax_outer_core, &
            reclen_ymin_outer_core, reclen_ymax_outer_core
  integer :: reclen_zmin

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_OUTER_CORE) :: vector_accel_outer_core,&
             vector_displ_outer_core, b_vector_displ_outer_core

  ! arrays to couple with the fluid regions by pointwise matching
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_OC) :: ibelm_xmin_outer_core,ibelm_xmax_outer_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_OC) :: ibelm_ymin_outer_core,ibelm_ymax_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core
  integer, dimension(NSPEC2D_TOP_OC) :: ibelm_top_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: normal_xmin_outer_core,normal_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: normal_ymin_outer_core,normal_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC) :: normal_top_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core

  ! assembly
  integer :: npoin2D_faces_outer_core(NUMFACES_SHARED)
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_outer_core,npoin2D_eta_outer_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL_OC,NUMCORNERS_SHARED) :: iboolcorner_outer_core

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_OC) :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_OC) :: iboolleft_eta_outer_core,iboolright_eta_outer_core

  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_outer_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: &
    rho_kl_outer_core,alpha_kl_outer_core

  ! kernel runs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: &
    div_displ_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: &
    b_div_displ_outer_core

  ! check for deviatoric kernel for outer core region
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: beta_kl_outer_core
  integer :: nspec_beta_kl_outer_core
  logical,parameter:: deviatoric_outercore = .false.

  ! inner / outer elements outer core region
  integer :: num_phase_ispec_outer_core
  integer :: nspec_inner_outer_core,nspec_outer_outer_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_outer_core


end module specfem_par_outercore


!=====================================================================

module specfem_par_movie

! parameter module for movies/shakemovies

  use constants

  implicit none

  ! to save movie frames
  integer :: nmovie_points,NIT
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
      store_val_x,store_val_y,store_val_z, &
      store_val_ux,store_val_uy,store_val_uz
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
      store_val_x_all,store_val_y_all,store_val_z_all, &
      store_val_ux_all,store_val_uy_all,store_val_uz_all

  ! to save movie volume
  integer :: npoints_3dmovie,nspecel_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE) :: num_ibool_3dmovie
  double precision :: scalingval
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) :: &
    muvstore_crust_mantle_3dmovie
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nu_3dmovie
  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: mask_3dmovie

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) :: &
    Iepsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: &
    Ieps_trace_over_3_crust_mantle

end module specfem_par_movie
