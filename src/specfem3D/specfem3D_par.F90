!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

module constants_solver

  use constants

  implicit none

! daniel debug: todo
#ifdef USE_STATIC_COMPILATION
  ! static compilation

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

#else
  ! "dynamic" compilation
  !
  ! we will read in mesh parameters stored in file DATABASES_MPI/mesh_parameters.bin at runtime
  !
  ! no static compilation, will use dynamic arrays:
  ! all these parameters must be stored in mesh files and read in by solver before allocating arrays...

  ! mesh
  integer :: NEX_XI_VAL
  integer :: NEX_ETA_VAL

  integer :: NSPEC_CRUST_MANTLE
  integer :: NSPEC_OUTER_CORE
  integer :: NSPEC_INNER_CORE
  integer :: NSPEC_TRINFINITE
  integer :: NSPEC_INFINITE

  integer :: NGLOB_CRUST_MANTLE
  integer :: NGLOB_OUTER_CORE
  integer :: NGLOB_INNER_CORE
  integer :: NGLOB_TRINFINITE
  integer :: NGLOB_INFINITE

  ! element types
  integer :: NSPECMAX_ANISO_IC
  integer :: NSPECMAX_ISO_MANTLE
  integer :: NSPECMAX_TISO_MANTLE
  integer :: NSPECMAX_ANISO_MANTLE

  ! attenuation & strain
  integer :: NSPEC_CRUST_MANTLE_ATTENUATION
  integer :: NSPEC_INNER_CORE_ATTENUATION

  integer :: NSPEC_CRUST_MANTLE_STR_OR_ATT
  integer :: NSPEC_INNER_CORE_STR_OR_ATT

  integer :: NSPEC_CRUST_MANTLE_STR_AND_ATT
  integer :: NSPEC_INNER_CORE_STR_AND_ATT

  integer :: NSPEC_CRUST_MANTLE_STRAIN_ONLY
  integer :: NSPEC_INNER_CORE_STRAIN_ONLY

  ! adjoint
  integer :: NSPEC_CRUST_MANTLE_ADJOINT
  integer :: NSPEC_OUTER_CORE_ADJOINT
  integer :: NSPEC_INNER_CORE_ADJOINT
  integer :: NSPEC_TRINFINITE_ADJOINT
  integer :: NSPEC_INFINITE_ADJOINT

  integer :: NGLOB_CRUST_MANTLE_ADJOINT
  integer :: NGLOB_OUTER_CORE_ADJOINT
  integer :: NGLOB_INNER_CORE_ADJOINT
  integer :: NGLOB_TRINFINITE_ADJOINT
  integer :: NGLOB_INFINITE_ADJOINT

  integer :: NSPEC_OUTER_CORE_ROT_ADJOINT

  ! absorbing boundary
  integer :: NSPEC_CRUST_MANTLE_STACEY
  integer :: NSPEC_OUTER_CORE_STACEY

  ! ocean/bathymetry
  integer :: NGLOB_CRUST_MANTLE_OCEANS

  integer :: NX_BATHY_VAL
  integer :: NY_BATHY_VAL

  ! rotation
  integer :: NSPEC_OUTER_CORE_ROTATION

  ! MPI partitions
  integer :: NPROC_XI_VAL
  integer :: NPROC_ETA_VAL
  integer :: NCHUNKS_VAL
  integer :: NPROCTOT_VAL

  ! attenuation
  integer :: ATT1_VAL
  integer :: ATT2_VAL
  integer :: ATT3_VAL
  integer :: ATT4_VAL
  integer :: ATT5_VAL

  ! boundaries
  integer :: NSPEC2DMAX_XMIN_XMAX_CM
  integer :: NSPEC2DMAX_YMIN_YMAX_CM
  integer :: NSPEC2D_BOTTOM_CM
  integer :: NSPEC2D_TOP_CM

  integer :: NSPEC2DMAX_XMIN_XMAX_IC
  integer :: NSPEC2DMAX_YMIN_YMAX_IC
  integer :: NSPEC2D_BOTTOM_IC
  integer :: NSPEC2D_TOP_IC

  integer :: NSPEC2DMAX_XMIN_XMAX_OC
  integer :: NSPEC2DMAX_YMIN_YMAX_OC
  integer :: NSPEC2D_BOTTOM_OC
  integer :: NSPEC2D_TOP_OC

  integer :: NSPEC2DMAX_XMIN_XMAX_TRINF
  integer :: NSPEC2DMAX_YMIN_YMAX_TRINF
  integer :: NSPEC2D_BOTTOM_TRINF
  integer :: NSPEC2D_TOP_TRINF

  integer :: NSPEC2DMAX_XMIN_XMAX_INF
  integer :: NSPEC2DMAX_YMIN_YMAX_INF
  integer :: NSPEC2D_BOTTOM_INF
  integer :: NSPEC2D_TOP_INF

  integer :: NSPEC2D_MOHO
  integer :: NSPEC2D_400
  integer :: NSPEC2D_670
  integer :: NSPEC2D_CMB
  integer :: NSPEC2D_ICB

  ! movies
  integer :: NSPEC_CRUST_MANTLE_3DMOVIE
  integer :: NGLOB_CRUST_MANTLE_3DMOVIE
  integer :: NSPEC_OUTER_CORE_3DMOVIE

  integer :: NGLOB_XY_CM
  integer :: NGLOB_XY_IC
  integer :: NT_DUMP_ATTENUATION_VAL

  ! model
  logical :: TRANSVERSE_ISOTROPY_VAL
  logical :: ANISOTROPIC_3D_MANTLE_VAL
  logical :: ANISOTROPIC_INNER_CORE_VAL
  logical :: ATTENUATION_VAL
  logical :: ATTENUATION_3D_VAL
  logical :: ELLIPTICITY_VAL
  logical :: GRAVITY_VAL
  logical :: OCEANS_VAL

  ! full gravity support
  logical :: FULL_GRAVITY_VAL

  logical :: ROTATION_VAL
  logical :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL
  logical :: PARTIAL_PHYS_DISPERSION_ONLY_VAL

  logical :: USE_DEVILLE_PRODUCTS_VAL
  logical :: ATTENUATION_1D_WITH_3D_STORAGE_VAL
  logical :: UNDO_ATTENUATION_VAL

  ! 1-chunk
  double precision :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL
  double precision :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL
  double precision :: CENTER_LATITUDE_IN_DEGREES_VAL
  double precision :: CENTER_LONGITUDE_IN_DEGREES_VAL
  double precision :: GAMMA_ROTATION_AZIMUTH_VAL

  ! we use this vectorization flag for solver routines in files **.f90
#ifdef FORCE_VECTORIZATION
  logical, parameter :: FORCE_VECTORIZATION_VAL = .true.
#else
  logical, parameter :: FORCE_VECTORIZATION_VAL = .false.
#endif

#endif

end module constants_solver

!=====================================================================

module specfem_par

! main parameter module for specfem simulations

  use constants_solver

  use shared_parameters

  implicit none

  !-----------------------------------------------------------------
  ! GLL points & weights
  !-----------------------------------------------------------------

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  ! product of weights for gravity term
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  ! arrays for Deville and force_vectorization
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  !-----------------------------------------------------------------
  ! attenuation parameters
  !-----------------------------------------------------------------

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval, betaval, gammaval
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: b_alphaval, b_betaval, b_gammaval

  ! attenuation: predictor
  double precision, dimension(N_SLS) :: tau_sigma_dble

  !-----------------------------------------------------------------
  ! topography/bathymetry & oceans
  !-----------------------------------------------------------------

  ! use integer array to store values
  integer, dimension(:,:),allocatable :: ibathy_topo

  ! additional mass matrix for ocean load
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: rmass_ocean_load

  integer :: npoin_oceans
  integer, dimension(:),allocatable :: ibool_ocean_load
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: normal_ocean_load
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: rmass_ocean_load_selected

  !-----------------------------------------------------------------
  ! ellipticity
  !-----------------------------------------------------------------

  ! for ellipticity
  integer :: nspl_ellip
  double precision,dimension(NR_DENSITY) :: rspl_ellip,ellipicity_spline,ellipicity_spline2

  !-----------------------------------------------------------------
  ! rotation
  !-----------------------------------------------------------------

  ! non-dimensionalized rotation rate of the Earth times two
  real(kind=CUSTOM_REAL) :: two_omega_earth
  !ADJOINT
  real(kind=CUSTOM_REAL) :: b_two_omega_earth

  !-----------------------------------------------------------------
  ! gravity
  !-----------------------------------------------------------------

  ! lookup table every km for gravity
  real(kind=CUSTOM_REAL) :: minus_g_cmb,minus_g_icb
  double precision, dimension(:),allocatable :: minus_gravity_table, &
    minus_deriv_gravity_table,density_table,d_ln_density_dr_table,minus_rho_g_over_kappa_fluid

  ! pre-computed vectors
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: gravity_pre_store_outer_core,gravity_pre_store_crust_mantle, &
    gravity_pre_store_inner_core
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: gravity_H_crust_mantle,gravity_H_inner_core

  !-----------------------------------------------------------------
  ! sources
  !-----------------------------------------------------------------

  ! parameters for the source
  integer :: nsources_local
  integer, dimension(:), allocatable :: islice_selected_source,ispec_selected_source
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sourcearrays
  double precision, dimension(:,:,:) ,allocatable:: nu_source

  double precision, dimension(:), allocatable :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(:), allocatable :: xi_source,eta_source,gamma_source
  double precision, dimension(:), allocatable :: tshift_src,hdur,hdur_Gaussian

  double precision :: Mrr,Mtt,Mpp,Mrt,Mrp,Mtp,Mw,M0
  double precision :: t0
  double precision :: min_tshift_src_original
  double precision :: source_final_distance_max
  double precision :: source_theta_ref,source_phi_ref   ! reference (geocentric) theta/phi position for source

  ! External source time function.
  double precision, dimension(:), allocatable :: user_source_time_function

  ! parameters for a force source located exactly at a grid point
  integer, dimension(:), allocatable :: force_stf
  double precision, dimension(:), allocatable :: factor_force_source
  double precision, dimension(:), allocatable :: comp_dir_vect_source_E
  double precision, dimension(:), allocatable :: comp_dir_vect_source_N
  double precision, dimension(:), allocatable :: comp_dir_vect_source_Z_UP

  !-----------------------------------------------------------------
  ! receivers
  !-----------------------------------------------------------------

  ! receiver information
  integer :: nrec,nrec_local
  integer, dimension(:), allocatable :: islice_selected_rec,ispec_selected_rec
  integer, dimension(:), allocatable :: islice_num_rec_local

  double precision, dimension(:), allocatable :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(:,:,:), allocatable :: nu_rec
  double precision, allocatable, dimension(:) :: stlat,stlon,stele,stbur
  double precision :: receiver_final_distance_max

  character(len=MAX_LENGTH_STATION_NAME), dimension(:), allocatable :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(:), allocatable :: network_name
  character(len=MAX_STRING_LEN) :: STATIONS_FILE

  ! Lagrange interpolators at receivers
  integer, dimension(:), allocatable, target :: number_receiver_global
  double precision, dimension(:,:), allocatable, target :: hxir_store,hetar_store,hgammar_store

  ! adjoint sources
  integer :: nadj_rec_local
  integer, dimension(:), pointer :: number_adjsources_global
  double precision, dimension(:,:), pointer :: hxir_adjstore,hetar_adjstore,hgammar_adjstore
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: source_adjoint
  ! asynchronous read buffer when IO_ASYNC_COPY is set
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_source_adjoint

  integer :: NSTEP_SUB_ADJ  ! to read input in chunks

  integer, dimension(:,:), allocatable :: iadjsrc ! to read input in chunks
  integer, dimension(:), allocatable :: iadjsrc_len,iadj_vec

  ! source Frechet derivatives
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: moment_der
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: sloc_der
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: stshift_der, shdur_der
  double precision, dimension(:,:), allocatable :: hpxir_store,hpetar_store,hpgammar_store


  !-----------------------------------------------------------------
  ! seismograms
  !-----------------------------------------------------------------

  ! seismograms
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms

  integer :: nlength_seismogram
  integer :: seismo_offset, seismo_current
  logical :: do_save_seismograms

  ! strain seismograms
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms_eps

  ! for SAC headers for seismograms
  integer :: yr_SAC,jda_SAC,mo_SAC,da_SAC,ho_SAC,mi_SAC
  double precision :: sec_SAC

  real :: mb_SAC ! body-wave magnitude
  real :: ms_SAC ! surface-wave magnitude (for ASDF QuakeML file)
  double precision :: t_cmt_SAC,t_shift_SAC
  double precision :: elat_SAC,elon_SAC,depth_SAC, &
    cmt_lat_SAC,cmt_lon_SAC,cmt_depth_SAC,cmt_hdur_SAC
  character(len=20) :: event_name_SAC

  ! strain flag
  logical :: COMPUTE_AND_STORE_STRAIN

  ! process/partition name
  character(len=MAX_STRING_LEN) :: prname

  !-----------------------------------------------------------------
  ! MPI partitions
  !-----------------------------------------------------------------

  integer :: ichunk ! needed for Stacey boundaries
  integer, dimension(:,:,:), allocatable :: addressing

  ! time loop timing
  double precision :: time_start

  !-----------------------------------------------------------------
  ! assembly
  !-----------------------------------------------------------------

  ! collected MPI interfaces
  ! MPI crust/mantle mesh
  integer :: num_interfaces_crust_mantle
  integer :: max_nibool_interfaces_cm
  integer, dimension(:), allocatable :: my_neighbors_crust_mantle,nibool_interfaces_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle

  ! note: MPI buffers are declared as pointers instead of allocatable arrays.
  !       this will allow for CUDA-aware MPI handling, where buffers have to be allocated on the GPU device.
  real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: b_buffer_send_vector_cm,b_buffer_recv_vector_cm

  integer, dimension(:), allocatable :: request_send_vector_cm,request_recv_vector_cm
  integer, dimension(:), allocatable :: b_request_send_vector_cm,b_request_recv_vector_cm

  ! MPI inner core mesh
  integer :: num_interfaces_inner_core
  integer :: max_nibool_interfaces_ic
  integer, dimension(:), allocatable :: my_neighbors_inner_core,nibool_interfaces_inner_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: buffer_send_vector_inner_core,buffer_recv_vector_inner_core
  real(kind=CUSTOM_REAL), dimension(:,:,:), pointer :: b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core

  integer, dimension(:), allocatable :: request_send_vector_ic,request_recv_vector_ic
  integer, dimension(:), allocatable :: b_request_send_vector_ic,b_request_recv_vector_ic

  ! MPI outer core mesh
  integer :: num_interfaces_outer_core
  integer :: max_nibool_interfaces_oc
  integer, dimension(:), allocatable :: my_neighbors_outer_core,nibool_interfaces_outer_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core

  integer, dimension(:), allocatable :: request_send_scalar_oc,request_recv_scalar_oc
  integer, dimension(:), allocatable :: b_request_send_scalar_oc,b_request_recv_scalar_oc

  !-----------------------------------------------------------------
  ! GPU
  !-----------------------------------------------------------------

  ! CUDA mesh pointer to integer wrapper
  integer(kind=8) :: Mesh_pointer

  ! flags for CUDA-aware MPI handling
  logical :: USE_CUDA_AWARE_MPI = .false.

  !-----------------------------------------------------------------
  ! ADIOS
  !-----------------------------------------------------------------
  ! adios file handle
  integer(kind=8) :: current_adios_handle

  !-----------------------------------------------------------------
  ! ASDF
  !-----------------------------------------------------------------
  ! asdf file handle
  !
  ! note: ASDF uses hdf5 file i/o routines. therefore, ASDF c routines define the file_id handle as hid_t.
  !       the datatype hid_t is defined by the hdf5 library, and for Fortran in file H5f90i_gen.h as:
  !          #define c_int_8 long long
  !          typedef c_int_8 hid_t_f
  !       in Fortran codes, one could use the hdf5 module for this
  !          use hdf5, only: HID_T
  !          integer(HID_T) :: file_id
  !       which will required the hdf5 library paths set for compilation and linking.
  !       instead here, the c_int_8 corresponds to long long, which in Fortran would be an 8-byte integer
  integer(kind=8) :: current_asdf_handle

  !-----------------------------------------------------------------
  ! time scheme
  !-----------------------------------------------------------------

  integer :: it
  integer :: it_begin,it_end

  ! non-dimensionalization
  double precision :: scale_t,scale_t_inv,scale_displ,scale_displ_inv,scale_veloc

  ! time scheme parameters
  real(kind=CUSTOM_REAL) :: deltat,deltatover2,deltatsqover2
  ! ADJOINT
  real(kind=CUSTOM_REAL) :: b_deltat,b_deltatover2,b_deltatsqover2

  ! LDDRK time scheme
  integer :: NSTAGE_TIME_SCHEME,istage
  real(kind=CUSTOM_REAL),dimension(N_SLS) :: tau_sigmainv_CUSTOM_REAL

  ! UNDO_ATTENUATION
  integer :: NT_DUMP_ATTENUATION,NSUBSET_ITERATIONS
  integer :: iteration_on_subset,it_of_this_subset
  integer :: it_subset_end

! if one wants to detect slow nodes compared to a reference time on normal nodes on a given cluster
! and exclude them from the runs in order not to slow down all the others when NUMBER_OF_SIMULTANEOUS_RUNS > 1.
  logical :: I_am_running_on_a_slow_node

  ! serial i/o mesh reading
#ifdef USE_SERIAL_CASCADE_FOR_IOs
  logical :: you_can_start_doing_IOs
#endif

  ! for EXACT_UNDOING_TO_DISK
  integer, dimension(:), allocatable :: integer_mask_ibool_exact_undo
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_for_disk
  integer :: it_exact_subset
  integer :: it_of_this_exact_subset, it_exact_subset_offset, it_exact_subset_end

  ! for optimized arrays
  logical :: use_inversed_arrays

  ! for saving/reading stacey boundary contributions
  logical :: SAVE_STACEY

  !-----------------------------------------------------------------
  ! point search
  !-----------------------------------------------------------------
  ! (i,j,k) indices of the control/anchor points of the element
  integer :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! coordinates of element midpoints
  double precision, allocatable, dimension(:,:) :: xyz_midpoints

  ! search range: single slice dimensions
  double precision :: lat_min,lat_max,lon_min,lon_max

  ! search margin in degrees
  double precision,parameter :: LAT_LON_MARGIN = 2.d0

  ! estimated typical element size (at surface)
  double precision :: element_size

  ! typical element search distance squared
  double precision :: typical_size_squared

  ! adjacency arrays
  integer,dimension(:),allocatable :: xadj   ! adjacency indexing
  integer,dimension(:),allocatable :: adjncy ! adjacency
  integer :: num_neighbors_all

end module specfem_par


!=====================================================================

module specfem_par_crustmantle

! parameter module for elastic solver in crust/mantle region

  use constants_solver
  implicit none

  ! ----------------- crust, mantle and oceans ---------------------
  ! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool_crust_mantle

  ! optimized arrays
  integer, dimension(:,:), allocatable :: ibool_inv_tbl_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_inv_st_crust_mantle
  integer, dimension(:,:), allocatable :: phase_iglob_crust_mantle
  integer, dimension(2) :: num_globs_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    deriv_mapping_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sum_terms_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: rstore_crust_mantle

  ! arrays for isotropic elements stored only where needed to save space
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle

  ! arrays for anisotropic elements stored only where needed to save space
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle

  ! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
    c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
    c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
    c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
    c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
    c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
    c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: mu0store_crust_mantle

  ! flag for transversely isotropic elements
  logical, dimension(:),allocatable :: ispec_is_tiso_crust_mantle

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  real(kind=CUSTOM_REAL), dimension(:), allocatable, target :: rmassz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:), pointer :: rmassx_crust_mantle,rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassx_crust_mantle,b_rmassy_crust_mantle

  ! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
     displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    one_minus_sum_beta_crust_mantle, factor_scale_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    factor_common_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle,R_xz_crust_mantle,R_yz_crust_mantle
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle,b_R_xz_crust_mantle,b_R_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, target :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, target :: &
    eps_trace_over_3_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), pointer :: &
    b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
    b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), pointer :: &
    b_eps_trace_over_3_crust_mantle

  ! coupling/boundary surfaces
  ! for coupling fluid-solid interfaces (at CMB and ICB)
  integer :: nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
             nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle
  ! for crust/oceans and CMB coupling
  integer, dimension(:), allocatable :: ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle
  integer, dimension(:), allocatable :: ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle
  integer, dimension(:), allocatable :: ibelm_bottom_crust_mantle
  integer, dimension(:), allocatable :: ibelm_top_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_top_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin_crust_mantle,normal_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_ymin_crust_mantle,normal_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_top_crust_mantle

  ! Stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp_crust_mantle,rho_vs_crust_mantle

  ! old version
  !integer, dimension(:,:), allocatable :: nimin_crust_mantle,nimax_crust_mantle,nkmin_eta_crust_mantle
  !integer, dimension(:,:), allocatable :: njmin_crust_mantle,njmax_crust_mantle,nkmin_xi_crust_mantle
  ! buffer arrays for file I/O
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: absorb_xmin_crust_mantle, &
  !  absorb_xmax_crust_mantle, absorb_ymin_crust_mantle, absorb_ymax_crust_mantle, absorb_zmin_crust_mantle
  !integer :: reclen_xmin_crust_mantle, reclen_xmax_crust_mantle, &
  !           reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
  !           reclen_zmin_crust_mantle

  ! absorbing boundary arrays
  integer :: num_abs_boundary_faces_crust_mantle
  integer, dimension(:), allocatable :: abs_boundary_ispec_crust_mantle
  integer, dimension(:), allocatable :: abs_boundary_npoin_crust_mantle
  integer, dimension(:,:,:), allocatable :: abs_boundary_ijk_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: abs_boundary_normal_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian2Dw_crust_mantle
  ! buffer array for file I/O
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: absorb_buffer_crust_mantle
  integer :: reclen_absorb_buffer_crust_mantle

  ! kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    rho_kl_crust_mantle,beta_kl_crust_mantle,alpha_kl_crust_mantle
  ! For anisotropic kernels (see compute_kernels.f90 for a definition of the array)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: cijkl_kl_crust_mantle
  ! approximate Hessian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_kl_crust_mantle
  ! approximate Hessian using rho, kappa and mu parametrization (diagonal part only)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_rho_kl_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_kappa_kl_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: hess_mu_kl_crust_mantle

  ! Boundary Mesh and Kernels
  integer, dimension(:), allocatable :: ibelm_moho_top,ibelm_moho_bot
  integer, dimension(:), allocatable :: ibelm_400_top,ibelm_400_bot
  integer, dimension(:), allocatable :: ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_moho
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_400
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_670
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: moho_kl, moho_kl_top, moho_kl_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: d400_kl, d400_kl_top, d400_kl_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: d670_kl, d670_kl_top, d670_kl_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: cmb_kl, cmb_kl_top, cmb_kl_bot

  ! regular kernels
  ! For saving kernels on a regular grid
  integer :: npoints_slice_reg
  integer, dimension(:), allocatable :: points_slice_reg
  integer, dimension(:), allocatable :: ispec_reg
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: hxir_reg
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: hetar_reg
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: hgammar_reg
  type kl_reg_grid_variables
    sequence
    real :: dlat
    real :: dlon
    real,dimension(:),allocatable :: rlayer
    integer :: nlayer
    integer :: npts_total ! order matters with sequence keyword, otherwise a misalignement warning appears
    integer,dimension(:),allocatable :: ndoubling
    integer,dimension(:),allocatable :: nlat
    integer,dimension(:),allocatable :: nlon
    integer,dimension(:),allocatable :: npts_before_layer
  end type kl_reg_grid_variables
  type(kl_reg_grid_variables) :: kl_reg_grid

  ! inner / outer elements crust/mantle region
  integer :: num_phase_ispec_crust_mantle
  integer :: nspec_inner_crust_mantle,nspec_outer_crust_mantle
  integer, dimension(:,:), allocatable :: phase_ispec_inner_crust_mantle

  ! number of surface elements
  integer :: NSPEC_TOP

  ! mesh coloring
  integer :: num_colors_outer_crust_mantle,num_colors_inner_crust_mantle
  integer,dimension(:),allocatable :: num_elem_colors_crust_mantle

  ! LDDRK time scheme
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: displ_crust_mantle_lddrk,veloc_crust_mantle_lddrk
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: &
    R_xx_crust_mantle_lddrk,R_yy_crust_mantle_lddrk,R_xy_crust_mantle_lddrk, &
    R_xz_crust_mantle_lddrk,R_yz_crust_mantle_lddrk
  ! adjoint
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_displ_crust_mantle_lddrk,b_veloc_crust_mantle_lddrk
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: &
    b_R_xx_crust_mantle_lddrk,b_R_yy_crust_mantle_lddrk,b_R_xy_crust_mantle_lddrk, &
    b_R_xz_crust_mantle_lddrk,b_R_yz_crust_mantle_lddrk

end module specfem_par_crustmantle

!=====================================================================

module specfem_par_innercore

! parameter module for elastic solver in inner core region

  use constants_solver
  implicit none

  ! ----------------- inner core ---------------------
  ! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool_inner_core

  ! optimized arrays
  integer, dimension(:,:), allocatable :: ibool_inv_tbl_inner_core
  integer, dimension(:,:), allocatable :: ibool_inv_st_inner_core
  integer, dimension(:,:), allocatable :: phase_iglob_inner_core
  integer, dimension(2) :: num_globs_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    xix_inner_core,xiy_inner_core,xiz_inner_core, &
    etax_inner_core,etay_inner_core,etaz_inner_core, &
    gammax_inner_core,gammay_inner_core,gammaz_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    deriv_mapping_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: sum_terms_inner_core

  ! material parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    xstore_inner_core,ystore_inner_core,zstore_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rstore_inner_core

  ! arrays for inner-core anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store_inner_core,c33store_inner_core,c12store_inner_core, &
    c13store_inner_core,c44store_inner_core

  ! local to global mapping
  integer, dimension(:), allocatable :: idoubling_inner_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable, target :: rmassz_inner_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassz_inner_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: rmassx_inner_core,rmassy_inner_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassx_inner_core,b_rmassy_inner_core

  ! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    displ_inner_core,veloc_inner_core,accel_inner_core
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    one_minus_sum_beta_inner_core, factor_scale_inner_core
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    factor_common_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    R_xx_inner_core,R_yy_inner_core,R_xy_inner_core,R_xz_inner_core,R_yz_inner_core
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core,b_R_xz_inner_core,b_R_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, target :: &
    epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
    epsilondev_xz_inner_core,epsilondev_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, target :: &
    eps_trace_over_3_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), pointer :: &
    b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
    b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), pointer :: &
    b_eps_trace_over_3_inner_core

  ! coupling/boundary surfaces
  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
             nspec2D_ymin_inner_core,nspec2D_ymax_inner_core
  integer, dimension(:), allocatable :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(:), allocatable :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(:), allocatable :: ibelm_bottom_inner_core
  integer, dimension(:), allocatable :: ibelm_top_inner_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    rho_kl_inner_core,beta_kl_inner_core, alpha_kl_inner_core

  ! Boundary Mesh and Kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: icb_kl, icb_kl_top, icb_kl_bot
  logical :: fluid_solid_boundary

  ! inner / outer elements inner core region
  integer :: num_phase_ispec_inner_core
  integer :: nspec_inner_inner_core,nspec_outer_inner_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_inner_core

  ! mesh coloring
  integer :: num_colors_outer_inner_core,num_colors_inner_inner_core
  integer,dimension(:),allocatable :: num_elem_colors_inner_core

  ! LDDRK time scheme
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: displ_inner_core_lddrk,veloc_inner_core_lddrk
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: &
    R_xx_inner_core_lddrk,R_yy_inner_core_lddrk,R_xy_inner_core_lddrk, &
    R_xz_inner_core_lddrk,R_yz_inner_core_lddrk
  ! adjoint
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: b_displ_inner_core_lddrk,b_veloc_inner_core_lddrk
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: &
    b_R_xx_inner_core_lddrk,b_R_yy_inner_core_lddrk,b_R_xy_inner_core_lddrk, &
    b_R_xz_inner_core_lddrk,b_R_yz_inner_core_lddrk

end module specfem_par_innercore

!=====================================================================

module specfem_par_outercore

! parameter module for acoustic solver in outer core region

  use constants_solver
  implicit none

  ! ----------------- outer core ---------------------
  ! mesh parameters
  integer, dimension(:,:,:,:), allocatable :: ibool_outer_core

  ! optimized arrays
  integer, dimension(:,:), allocatable :: ibool_inv_tbl_outer_core
  integer, dimension(:,:), allocatable :: ibool_inv_st_outer_core
  integer, dimension(:,:), allocatable :: phase_iglob_outer_core
  integer, dimension(2) :: num_globs_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    xix_outer_core,xiy_outer_core,xiz_outer_core, &
    etax_outer_core,etay_outer_core,etaz_outer_core, &
    gammax_outer_core,gammay_outer_core,gammaz_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    deriv_mapping_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sum_terms_outer_core

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    xstore_outer_core,ystore_outer_core,zstore_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: rstore_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore_outer_core,kappavstore_outer_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable, target :: rmass_outer_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmass_outer_core

  ! potentials
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    displ_outer_core,veloc_outer_core,accel_outer_core
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    A_array_rotation,B_array_rotation
  !ADJOINT
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    b_A_array_rotation,b_B_array_rotation

  ! Stacey
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: vp_outer_core

  ! old version...
  !integer, dimension(:,:), allocatable :: nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core
  !integer, dimension(:,:), allocatable :: njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core
  ! buffer arrays for file I/O
  !real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: absorb_xmin_outer_core, &
  !  absorb_xmax_outer_core, absorb_ymin_outer_core, absorb_ymax_outer_core, absorb_zmin_outer_core
  !integer :: reclen_xmin_outer_core, reclen_xmax_outer_core, &
  !           reclen_ymin_outer_core, reclen_ymax_outer_core, &
  !           reclen_zmin_outer_core

  ! absorbing boundary arrays
  integer :: num_abs_boundary_faces_outer_core
  integer, dimension(:), allocatable :: abs_boundary_ispec_outer_core
  integer, dimension(:), allocatable :: abs_boundary_npoin_outer_core
  integer, dimension(:,:,:), allocatable :: abs_boundary_ijk_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: abs_boundary_jacobian2Dw_outer_core
  ! buffer array for file I/O
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: absorb_buffer_outer_core
  integer :: reclen_absorb_buffer_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    vector_accel_outer_core,vector_displ_outer_core,b_vector_displ_outer_core

  ! coupling/boundary surfaces
  integer :: nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
             nspec2D_ymin_outer_core,nspec2D_ymax_outer_core
  ! arrays to couple with the fluid regions by pointwise matching
  integer, dimension(:), allocatable :: ibelm_xmin_outer_core,ibelm_xmax_outer_core
  integer, dimension(:), allocatable :: ibelm_ymin_outer_core,ibelm_ymax_outer_core
  integer, dimension(:), allocatable :: ibelm_bottom_outer_core
  integer, dimension(:), allocatable :: ibelm_top_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_xmin_outer_core,normal_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_ymin_outer_core,normal_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_top_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_kl_outer_core,alpha_kl_outer_core

  ! check for deviatoric kernel for outer core region
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: beta_kl_outer_core
  integer :: nspec_beta_kl_outer_core
  logical, parameter :: deviatoric_outercore = .false.

  ! inner / outer elements outer core region
  integer :: num_phase_ispec_outer_core
  integer :: nspec_inner_outer_core,nspec_outer_outer_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_outer_core

  ! mesh coloring
  integer :: num_colors_outer_outer_core,num_colors_inner_outer_core
  integer,dimension(:),allocatable :: num_elem_colors_outer_core

  ! LDDRK time scheme
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: displ_outer_core_lddrk,veloc_outer_core_lddrk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: A_array_rotation_lddrk,B_array_rotation_lddrk
  ! adjoint
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_displ_outer_core_lddrk,b_veloc_outer_core_lddrk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: b_A_array_rotation_lddrk,b_B_array_rotation_lddrk

end module specfem_par_outercore


!=====================================================================

module specfem_par_noise

! parameter module for noise simulations

  use constants_solver

  implicit none

  ! NOISE_TOMOGRAPHY
  integer :: reclen_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_sourcearray
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    normal_x_noise,normal_y_noise,normal_z_noise, mask_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: noise_surface_movie

  integer :: num_noise_surface_points
  integer :: irec_main_noise
  integer :: nsources_local_noise

  ! noise buffer for file i/o
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_buffer
  integer :: NT_DUMP_NOISE_BUFFER
  integer :: icounter_noise_buffer,nstep_subset_noise_buffer

  ! noise strength kernel
  ! crust/mantle region
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sigma_kl_crust_mantle

end module specfem_par_noise

!=====================================================================

module specfem_par_movie

! parameter module for movies/shakemovies

  use constants_solver

  implicit none

  ! to save movie frames
  integer :: nmovie_points,NIT

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_ux,store_val_uy,store_val_uz
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: store_val_ux_all,store_val_uy_all,store_val_uz_all

  ! to save movie volume
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    muvstore_crust_mantle_3dmovie

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    Iepsilondev_xx_crust_mantle,Iepsilondev_yy_crust_mantle,Iepsilondev_xy_crust_mantle, &
    Iepsilondev_xz_crust_mantle,Iepsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    Ieps_trace_over_3_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nu_3dmovie

  integer :: npoints_3dmovie,nspecel_3dmovie

  logical, dimension(:,:,:,:), allocatable :: mask_3dmovie

  ! outer core
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    div_displ_outer_core

  ! VTK run-time visualization
#ifdef HAVE_VTK
  ! VTK window
  logical, parameter :: VTK_MODE = .true.
#else
  logical, parameter :: VTK_MODE = .false.
#endif
  real,dimension(:),allocatable :: vtkdata
  logical,dimension(:),allocatable :: vtkmask
  ! multi-MPI processes, gather data arrays on main
  real,dimension(:),allocatable :: vtkdata_all
  integer,dimension(:),allocatable :: vtkdata_points_all
  integer,dimension(:),allocatable :: vtkdata_offset_all
  integer :: vtkdata_numpoints_all
  real :: vtkdata_source_x,vtkdata_source_y,vtkdata_source_z

end module specfem_par_movie

!=====================================================================

module specfem_par_full_gravity

  use constants_solver

  implicit none

  ! non-dimensionalization
  double precision :: scale_accel,scale_pgrav

  ! seismograms
  ! perturbed gravitational potential
  real(kind=CUSTOM_REAL),allocatable :: seismograms_phi(:,:,:)
  ! perturbed gravity
  real(kind=CUSTOM_REAL),allocatable :: seismograms_pgrav(:,:,:)
  ! perturbed gravity gradient
  real(kind=CUSTOM_REAL),allocatable :: seismograms_Hgrav(:,:,:,:)
  ! Due to background gravity. Free-air change in the gravity (vertical) or tilt of the ground surface(horizontal)
  real(kind=CUSTOM_REAL),allocatable :: seismograms_grav(:,:,:)
  ! Coriolis acceleration
  real(kind=CUSTOM_REAL),allocatable :: seismograms_corio(:,:,:)

  ! gradient of the background gravity at receiver location: \nabla g, where g = ||g||
  real(kind=CUSTOM_REAL),allocatable :: g_spec_rec(:,:,:),gradg_rec(:,:)

  ! Lagrange interpolators at receivers
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storederiv_rec

  ! full gravity solver
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: gravity_rho_crust_mantle,gravity_rho_inner_core
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: gravity_rho_g_over_kappa_outer_core

  ! MPI infinite mesh
  ! transition infinite
  integer :: num_interfaces_trinfinite
  integer :: max_nibool_interfaces_trinfinite
  integer, dimension(:), allocatable :: my_neighbors_trinfinite,nibool_interfaces_trinfinite
  integer, dimension(:,:), allocatable :: ibool_interfaces_trinfinite

  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: buffer_send_scalar_trinfinite,buffer_recv_scalar_trinfinite
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: b_buffer_send_scalar_trinfinite,b_buffer_recv_scalar_trinfinite

  integer, dimension(:), allocatable :: request_send_scalar_trinfinite,request_recv_scalar_trinfinite
  integer, dimension(:), allocatable :: b_request_send_scalar_trinfinite,b_request_recv_scalar_trinfinite

  ! infinite
  integer :: num_interfaces_infinite
  integer :: max_nibool_interfaces_infinite
  integer, dimension(:), allocatable :: my_neighbors_infinite,nibool_interfaces_infinite
  integer, dimension(:,:), allocatable :: ibool_interfaces_infinite

  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: buffer_send_scalar_infinite,buffer_recv_scalar_infinite
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: b_buffer_send_scalar_infinite,b_buffer_recv_scalar_infinite

  integer, dimension(:), allocatable :: request_send_scalar_infinite,request_recv_scalar_infinite
  integer, dimension(:), allocatable :: b_request_send_scalar_infinite,b_request_recv_scalar_infinite

  ! Level-1 solver
  ! MPI crust/mantle mesh
  integer :: num_interfaces_crust_mantle1
  integer :: max_nibool_interfaces_crust_mantle1
  integer, dimension(:), allocatable :: my_neighbors_crust_mantle1,nibool_interfaces_crust_mantle1
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle1

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_crust_mantle1,buffer_recv_vector_crust_mantle1
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_crust_mantle1,b_buffer_recv_vector_crust_mantle1

  integer, dimension(:), allocatable :: request_send_vector_crust_mantle1,request_recv_vector_crust_mantle1
  integer, dimension(:), allocatable :: b_request_send_vector_crust_mantle1,b_request_recv_vector_crust_mantle1

  ! MPI inner core mesh
  integer :: num_interfaces_inner_core1
  integer :: max_nibool_interfaces_inner_core1
  integer, dimension(:), allocatable :: my_neighbors_inner_core1,nibool_interfaces_inner_core1
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core1

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_inner_core1,buffer_recv_vector_inner_core1
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_inner_core1,b_buffer_recv_vector_inner_core1

  integer, dimension(:), allocatable :: request_send_vector_inner_core1,request_recv_vector_inner_core1
  integer, dimension(:), allocatable :: b_request_send_vector_inner_core1,b_request_recv_vector_inner_core1

  ! MPI outer core mesh
  integer :: num_interfaces_outer_core1
  integer :: max_nibool_interfaces_outer_core1
  integer, dimension(:), allocatable :: my_neighbors_outer_core1,nibool_interfaces_outer_core1
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core1

  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: buffer_send_scalar_outer_core1,buffer_recv_scalar_outer_core1
  real(kind=CUSTOM_REAL), dimension(:,:), pointer :: b_buffer_send_scalar_outer_core1,b_buffer_recv_scalar_outer_core1

  integer, dimension(:), allocatable :: request_send_scalar_outer_core1,request_recv_scalar_outer_core1
  integer, dimension(:), allocatable :: b_request_send_scalar_outer_core1,b_request_recv_scalar_outer_core1

  ! MPI infinite mesh
  ! transition infinite
  integer :: num_interfaces_trinfinite1
  integer :: max_nibool_interfaces_trinfinite1
  integer, dimension(:), allocatable :: my_neighbors_trinfinite1,nibool_interfaces_trinfinite1
  integer, dimension(:,:), allocatable :: ibool_interfaces_trinfinite1
  ! infinite
  integer :: num_interfaces_infinite1
  integer :: max_nibool_interfaces_infinite1
  integer, dimension(:), allocatable :: my_neighbors_infinite1,nibool_interfaces_infinite1
  integer, dimension(:,:), allocatable :: ibool_interfaces_infinite1

  ! parameters for Poisson's equation
  integer :: neq, nnode
  integer :: neq1, nnode1
  double precision,dimension(:,:),allocatable :: lagrange_gll
  double precision,dimension(:,:),allocatable :: lagrange_gll1

  ! Level-1 solver arrays
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: pgrav_ic1, pgrav_oc1, pgrav_cm1, pgrav_trinf1, pgrav_inf1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: pgrav1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon1, gravload1
  ! Level-2 solver arrays
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: pgrav_ic, pgrav_oc, pgrav_cm, pgrav_trinf, pgrav_inf
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: pgrav
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon, gravload

  ! Adjoint arrays for Poisson Equation
  ! Level-1 solver arrays
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: b_pgrav_ic1, b_pgrav_oc1, b_pgrav_cm1, b_pgrav_trinf1, b_pgrav_inf1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: b_pgrav1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: b_gravload1
  ! Level-2 solver arrays
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: b_pgrav_ic, b_pgrav_oc, b_pgrav_cm, b_pgrav_trinf, b_pgrav_inf
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: b_pgrav
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: b_gravload

  ! Level-2 solver
  ! number of global degrees of freedom
  integer :: ngdof,nsparse
  integer,dimension(:),allocatable :: l2gdof
  integer,dimension(:),allocatable :: krow_sparse, kcol_sparse, kgrow_sparse, kgcol_sparse

  ! Level-1 solver
  ! number of global degrees of freedom
  integer :: ngdof1,nsparse1
  integer,dimension(:),allocatable :: l2gdof1
  integer,dimension(:),allocatable :: krow_sparse1, kcol_sparse1, kgrow_sparse1, kgcol_sparse1

  integer :: nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1

  ! mirror nodes
  integer,dimension(:),allocatable :: nmir_ic, nmir_oc, nmir_cm, nmir_trinf, nmir_inf

  ! mapping
  integer,dimension(:,:),allocatable :: inode_map_ic,inode_map_oc,inode_map_cm,inode_map_trinf,inode_map_inf

  integer,dimension(:,:),allocatable :: inode_elmt_ic,inode_elmt_oc,inode_elmt_cm,inode_elmt_trinf,inode_elmt_inf
  integer,dimension(:,:),allocatable :: inode_elmt_ic1,inode_elmt_oc1,inode_elmt_cm1,inode_elmt_trinf1,inode_elmt_inf1

  ! active GLL points
  integer,dimension(NGLLCUBE_INF) :: igll_active_on
  logical,dimension(NGLLCUBE) :: is_active_gll

  ! parameters for ENSIGHT GOLD files
  !integer :: nnode4_ic,nnode4_oc,nnode4_cm
  !integer,allocatable :: inode4_ic(:),inode4_oc(:),inode4_cm(:)

  ! CG solver scaling
  logical, parameter :: CG_SCALING = .true.
  ! CG solver non-dim scaling
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: ndscale1, ndscale

  ! full gravity arrays for kernels
  ! crust/mantle
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: g_cm
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: gradg_cm
  ! inner core - kernels not implemented yet
  !real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: g_ic
  !real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: gradg_ic
  ! outer core - kernels not implemented yet
  !real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: g_oc

  ! parameters for Poisson's solver
  ! crust/mantle
  integer,dimension(:),allocatable :: gdof_cm, gdof_cm1
  integer,dimension(:,:),allocatable :: ggdof_cm, ggdof_cm1

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_crust_mantle
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_crust_mantle
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_crust_mantle1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_crust_mantle1

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: storederiv_cm
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: storerhojw_cm
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: storederiv_cm1
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: storerhojw_cm1,storejw_cm1

  ! outer core
  integer,dimension(:),allocatable :: gdof_oc, gdof_oc1
  integer,dimension(:,:),allocatable :: ggdof_oc, ggdof_oc1

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_outer_core
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_outer_core
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_outer_core1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_outer_core1

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: storederiv_oc
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: storerhojw_oc
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: storederiv_oc1
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: storerhojw_oc1

  ! inner core
  integer,dimension(:),allocatable :: gdof_ic, gdof_ic1
  integer,dimension(:,:),allocatable :: ggdof_ic, ggdof_ic1

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_inner_core
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_inner_core
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable ::  storekmat_inner_core1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_inner_core1

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: storederiv_ic
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: storerhojw_ic
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: storederiv_ic1
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: storerhojw_ic1

  ! transition-to-infinite
  integer,dimension(:),allocatable :: gdof_trinf, gdof_trinf1
  integer,dimension(:,:),allocatable :: ggdof_trinf, ggdof_trinf1

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_trinfinite
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_trinfinite
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_trinfinite1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_trinfinite1

  ! infinite
  integer,dimension(:),allocatable :: gdof_inf, gdof_inf1
  integer,dimension(:,:),allocatable :: ggdof_inf, ggdof_inf1

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_infinite
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_infinite
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: storekmat_infinite1
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: dprecon_infinite1

  ! Full gravity kernels
  ! density kernel 1 (related to gravity matrix \Pi)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: rho1siem_kl_crust_mantle
  ! density kernel 2
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: rho2siem_kl_crust_mantle
  ! 1. and 2. gravity kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable ::  gknl1, gknl2

  ! for debugging only: delete at some point
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: debug_rho_kl_cm

  ! for the Euler scheme for rotation in linear indexing
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: A_array_rotationL, B_array_rotationL

  ! for the Euler scheme for rotation in linear indexing for 3-GLLX points
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: A_array_rotationL3, B_array_rotationL3

  ! for the adjoint Euler scheme for rotation in linear indexing for 3-GLLX points
  real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: b_A_array_rotationL3, b_B_array_rotationL3

end module specfem_par_full_gravity

!=====================================================================

module specfem_par_trinfinite

! parameter module for acoustic solver in transition infinite region

  use constants_solver

  implicit none

  ! ----------------- trinfinite ---------------------
  ! mesh parameters
  integer,dimension(:,:,:,:),allocatable :: ibool_trinfinite

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: &
    xix_trinfinite, xiy_trinfinite, xiz_trinfinite, &
    etax_trinfinite, etay_trinfinite, etaz_trinfinite, &
    gammax_trinfinite, gammay_trinfinite, gammaz_trinfinite

  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore_trinfinite, ystore_trinfinite, zstore_trinfinite

  integer :: nspec2D_xmin_trinfinite,nspec2D_xmax_trinfinite, &
             nspec2D_ymin_trinfinite,nspec2D_ymax_trinfinite

  ! arrays to couple with the fluid regions by pointwise matching
  integer,dimension(:),allocatable :: ibelm_xmin_trinfinite,ibelm_xmax_trinfinite
  integer,dimension(:),allocatable :: ibelm_ymin_trinfinite,ibelm_ymax_trinfinite
  integer,dimension(:),allocatable :: ibelm_bottom_trinfinite
  integer,dimension(:),allocatable :: ibelm_top_trinfinite

  ! inner / outer elements outer core region
  integer :: num_phase_ispec_trinfinite
  integer :: nspec_inner_trinfinite,nspec_outer_trinfinite
  integer,dimension(:,:),allocatable :: phase_ispec_inner_trinfinite

  ! mesh coloring
  integer :: num_colors_outer_trinfinite,num_colors_inner_trinfinite
  integer,dimension(:),allocatable :: num_elem_colors_trinfinite

end module specfem_par_trinfinite

!=====================================================================

module specfem_par_infinite

! parameter module for acoustic solver in infinite region

  use constants_solver

  implicit none

  ! ----------------- infinite ---------------------
  ! mesh parameters
  integer,dimension(:,:,:,:),allocatable :: ibool_infinite

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: &
    xix_infinite, xiy_infinite, xiz_infinite, &
    etax_infinite, etay_infinite, etaz_infinite, &
    gammax_infinite, gammay_infinite, gammaz_infinite

  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore_infinite, ystore_infinite, zstore_infinite

  integer :: nspec2D_xmin_infinite,nspec2D_xmax_infinite, &
             nspec2D_ymin_infinite,nspec2D_ymax_infinite

  ! arrays to couple with the fluid regions by pointwise matching
  integer, dimension(:),allocatable :: ibelm_xmin_infinite,ibelm_xmax_infinite
  integer, dimension(:),allocatable :: ibelm_ymin_infinite,ibelm_ymax_infinite
  integer, dimension(:),allocatable :: ibelm_bottom_infinite
  integer, dimension(:),allocatable :: ibelm_top_infinite

  ! inner / outer elements outer core region
  integer :: num_phase_ispec_infinite
  integer :: nspec_inner_infinite,nspec_outer_infinite
  integer, dimension(:,:), allocatable :: phase_ispec_inner_infinite

  ! mesh coloring
  integer :: num_colors_outer_infinite,num_colors_inner_infinite
  integer,dimension(:),allocatable :: num_elem_colors_infinite

end module specfem_par_infinite

!=====================================================================

! not used yet...

!module specfem_par_freesurface
!
!  use constants_solver, only: CUSTOM_REAL,NGLLX,NGLLY
!
!  implicit none
!
!  integer :: nspec_fs,nnode_fs
!
!  integer,dimension(NGLLX,NGLLY) :: igll1d
!
!  integer,allocatable :: gnode_fs(:)
!  integer,allocatable :: nvalency_fs(:)
!
!  ! Connectivity renumberd to the free surface nodes.
!  integer,allocatable :: rgnum_fs(:,:)
!  real, dimension(:), allocatable :: div_fs,pgrav_fs
!  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: curl_fs
!
!end module specfem_par_freesurface

!=====================================================================

#ifdef USE_XSMM

module my_libxsmm

  ! for C_LOC(A)
  ! note: using C_LOC(A) however needs array A to be defined with "..,target :: A" attribute to become an interoperable pointer.
  !       we haven't had issues yet though with array pointers when passing to C-functions, thus we omit it for now.
  !use,intrinsic :: ISO_C_BINDING

  !use libxsmm, only: C_LOC,LIBXSMM_SMMfunction,libxsmm_dispatch,libxsmm_call

  ! function calls
  use libxsmm, only: libxsmm_init,libxsmm_finalize
  use libxsmm, only: libxsmm_dispatch,libxsmm_available,LIBXSMM_SMMfunction
  use libxsmm, only: libxsmm_smmcall_abc
  !use libxsmm, only: libxsmm_smm_25_5_5,libxsmm_smm_5_25_5,libxsmm_smm_5_5_5

  implicit none

  ! function pointers
  ! (note: defined for single precision, thus needs CUSTOM_REAL to be SIZE_REAL)
  type(LIBXSMM_SMMfunction) :: xmm1, xmm2, xmm3


end module my_libxsmm

#endif


