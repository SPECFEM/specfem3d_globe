!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

module constants_solver

  use constants

  implicit none

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

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
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

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
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load

  ! flag to mask ocean-bottom degrees of freedom for ocean load
  logical, dimension(NGLOB_CRUST_MANTLE_OCEANS) :: updated_dof_ocean_load

  integer :: npoin_oceans
  integer, dimension(:),allocatable :: ibool_ocean_load
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: normal_ocean_load
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: rmass_ocean_load_selected

  !-----------------------------------------------------------------
  ! ellipticity
  !-----------------------------------------------------------------

  ! for ellipticity
  integer :: nspl
  double precision,dimension(NR) :: rspl,espl,espl2

  !-----------------------------------------------------------------
  ! rotation
  !-----------------------------------------------------------------

  ! non-dimensionalized rotation rate of the Earth times two
  real(kind=CUSTOM_REAL) :: two_omega_earth
  !ADJOINT
  real(kind=CUSTOM_REAL) b_two_omega_earth

  !-----------------------------------------------------------------
  ! gravity
  !-----------------------------------------------------------------

  ! lookup table every km for gravity
  real(kind=CUSTOM_REAL) :: minus_g_cmb,minus_g_icb
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table, &
    minus_deriv_gravity_table,density_table,d_ln_density_dr_table,minus_rho_g_over_kappa_fluid

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
  double precision, dimension(:), allocatable :: tshift_cmt,hdur,hdur_gaussian
  double precision, dimension(:), allocatable :: theta_source,phi_source
  double precision :: t0

  ! External source time function.
  double precision, dimension(:), allocatable :: stfArray_external

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
  character(len=MAX_STRING_LEN) :: STATIONS,rec_filename

  ! Lagrange interpolators at receivers
  double precision, dimension(:,:), allocatable :: hxir_store,hetar_store,hgammar_store
  double precision, dimension(:,:,:,:), allocatable :: hlagrange_store

  ! ADJOINT sources
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: adj_sourcearrays
  ! asynchronous read buffer when IO_ASYNC_COPY is set
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:,:), allocatable :: buffer_sourcearrays

  integer :: nadj_rec_local
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
  integer :: it_begin,it_end
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: seismograms
  integer :: seismo_offset, seismo_current
  ! adjoint seismograms
  integer :: it_adj_written

  ! for SAC headers for seismograms
  integer :: yr_SAC,jda_SAC,ho_SAC,mi_SAC
  double precision :: sec_SAC
  real :: mb_SAC
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

  ! proc numbers for MPI
  integer :: myrank
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
  integer, dimension(:), allocatable :: my_neighbours_crust_mantle,nibool_interfaces_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_cm,b_buffer_recv_vector_cm

  integer, dimension(:), allocatable :: request_send_vector_cm,request_recv_vector_cm
  integer, dimension(:), allocatable :: b_request_send_vector_cm,b_request_recv_vector_cm

  ! MPI inner core mesh
  integer :: num_interfaces_inner_core
  integer :: max_nibool_interfaces_ic
  integer, dimension(:), allocatable :: my_neighbours_inner_core,nibool_interfaces_inner_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_inner_core,buffer_recv_vector_inner_core
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: b_buffer_send_vector_inner_core,b_buffer_recv_vector_inner_core

  integer, dimension(:), allocatable :: request_send_vector_ic,request_recv_vector_ic
  integer, dimension(:), allocatable :: b_request_send_vector_ic,b_request_recv_vector_ic

  ! MPI outer core mesh
  integer :: num_interfaces_outer_core
  integer :: max_nibool_interfaces_oc
  integer, dimension(:), allocatable :: my_neighbours_outer_core,nibool_interfaces_outer_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core

  integer, dimension(:), allocatable :: request_send_scalar_oc,request_recv_scalar_oc
  integer, dimension(:), allocatable :: b_request_send_scalar_oc,b_request_recv_scalar_oc

  !-----------------------------------------------------------------
  ! GPU
  !-----------------------------------------------------------------

  ! CUDA mesh pointer <-> integer wrapper
  integer(kind=8) :: Mesh_pointer

  !-----------------------------------------------------------------
  ! ADIOS
  !-----------------------------------------------------------------
  ! adios file handle
  integer(kind=8) :: current_adios_handle

  !-----------------------------------------------------------------
  ! time scheme
  !-----------------------------------------------------------------

  integer :: it

  ! non-dimensionalization
  double precision :: scale_t,scale_t_inv,scale_displ,scale_veloc

  ! time scheme parameters
  real(kind=CUSTOM_REAL) :: deltat,deltatover2,deltatsqover2
  ! ADJOINT
  real(kind=CUSTOM_REAL) :: b_deltat,b_deltatover2,b_deltatsqover2

  ! LDDRK time scheme
  integer :: NSTAGE_TIME_SCHEME,istage
  real(kind=CUSTOM_REAL),dimension(N_SLS) :: tau_sigma_CUSTOM_REAL

  ! undo_attenuation
  integer :: NSUBSET_ITERATIONS
  integer :: iteration_on_subset,it_of_this_subset

  ! serial i/o mesh reading
#ifdef USE_SERIAL_CASCADE_FOR_IOs
  logical :: you_can_start_doing_IOs
#endif

  ! for EXACT_UNDOING_TO_DISK
  integer, dimension(:), allocatable :: integer_mask_ibool_exact_undo
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: buffer_for_disk


end module specfem_par


!=====================================================================

module specfem_par_crustmantle

! parameter module for elastic solver in crust/mantle region

  use constants_solver
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

  ! flag for transversely isotropic elements
  logical, dimension(NSPEC_CRUST_MANTLE) :: ispec_is_tiso_crust_mantle

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE),target :: rmassz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:), pointer :: rmassx_crust_mantle,rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassx_crust_mantle,b_rmassy_crust_mantle

  ! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle

  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL) :: &
    one_minus_sum_beta_crust_mantle, factor_scale_crust_mantle
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,ATT4_VAL) :: &
    factor_common_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION) :: &
    R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle,R_xz_crust_mantle,R_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT),target :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY),target :: &
    eps_trace_over_3_crust_mantle

  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT) :: &
    b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle,b_R_xz_crust_mantle,b_R_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),pointer :: &
    b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
    b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),pointer :: &
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

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: absorb_xmin_crust_mantle, &
    absorb_xmax_crust_mantle, absorb_ymin_crust_mantle, absorb_ymax_crust_mantle

  integer :: reclen_xmin_crust_mantle, reclen_xmax_crust_mantle, &
            reclen_ymin_crust_mantle,reclen_ymax_crust_mantle

  ! kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    rho_kl_crust_mantle,beta_kl_crust_mantle,alpha_kl_crust_mantle

  ! noise strength kernel
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: sigma_kl_crust_mantle
  ! For anisotropic kernels (see compute_kernels.f90 for a definition of the array)
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: cijkl_kl_crust_mantle
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

  ! For saving kernels on a regular grid
  integer :: npoints_slice
  integer, dimension(NM_KL_REG_PTS_VAL) :: points_slice
  integer, dimension(NM_KL_REG_PTS_VAL) :: ispec_reg
  real(kind=CUSTOM_REAL), dimension(NGLLX, NM_KL_REG_PTS_VAL) :: hxir_reg
  real(kind=CUSTOM_REAL), dimension(NGLLY, NM_KL_REG_PTS_VAL) :: hetar_reg
  real(kind=CUSTOM_REAL), dimension(NGLLZ, NM_KL_REG_PTS_VAL) :: hgammar_reg

  ! NOISE_TOMOGRAPHY
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: noise_sourcearray
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    normal_x_noise,normal_y_noise,normal_z_noise, mask_noise
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: noise_surface_movie
  integer :: num_noise_surface_points
  integer :: irec_master_noise
  integer :: NSPEC_TOP

  ! inner / outer elements crust/mantle region
  integer :: num_phase_ispec_crust_mantle
  integer :: nspec_inner_crust_mantle,nspec_outer_crust_mantle
  integer, dimension(:,:), allocatable :: phase_ispec_inner_crust_mantle

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
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    xix_inner_core,xiy_inner_core,xiz_inner_core,&
    etax_inner_core,etay_inner_core,etaz_inner_core, &
    gammax_inner_core,gammay_inner_core,gammaz_inner_core

  ! material parameters
  ! (note: muvstore also needed for attenuation in case of anisotropic inner core)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
    rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: &
    xstore_inner_core,ystore_inner_core,zstore_inner_core

  ! arrays for inner-core anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC) :: &
    c11store_inner_core,c33store_inner_core,c12store_inner_core, &
    c13store_inner_core,c44store_inner_core

  ! local to global mapping
  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE), target :: rmassz_inner_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassz_inner_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: rmassx_inner_core,rmassy_inner_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmassx_inner_core,b_rmassy_inner_core

  ! displacement, velocity, acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: &
    displ_inner_core,veloc_inner_core,accel_inner_core
  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL) :: &
    one_minus_sum_beta_inner_core, factor_scale_inner_core
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,ATT5_VAL) :: &
    factor_common_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION) :: &
    R_xx_inner_core,R_yy_inner_core,R_xy_inner_core,R_xz_inner_core,R_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT),target :: &
    epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
    epsilondev_xz_inner_core,epsilondev_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY),target :: &
    eps_trace_over_3_inner_core

  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT) :: &
    b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core,b_R_xz_inner_core,b_R_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),pointer :: &
    b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
    b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),pointer :: &
    b_eps_trace_over_3_inner_core

  ! coupling/boundary surfaces
  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
            nspec2D_ymin_inner_core,nspec2D_ymax_inner_core
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_IC) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_IC) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_IC) :: ibelm_bottom_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

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
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
    xix_outer_core,xiy_outer_core,xiz_outer_core,&
    etax_outer_core,etay_outer_core,etaz_outer_core, &
    gammax_outer_core,gammay_outer_core,gammaz_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    xstore_outer_core,ystore_outer_core,zstore_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: rhostore_outer_core,kappavstore_outer_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE), target :: rmass_outer_core
  real(kind=CUSTOM_REAL), dimension(:), pointer :: b_rmass_outer_core

  ! velocity potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    displ_outer_core,veloc_outer_core,accel_outer_core

  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION) :: &
    A_array_rotation,B_array_rotation
  !ADJOINT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT) :: &
    b_A_array_rotation,b_B_array_rotation


  ! Stacey
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_STACEY) :: vp_outer_core
  integer :: nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
             nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
             nspec2D_zmin_outer_core
  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_OC) :: nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_OC) :: njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: absorb_xmin_outer_core, &
     absorb_xmax_outer_core, absorb_ymin_outer_core, absorb_ymax_outer_core, &
     absorb_zmin_outer_core

  integer :: reclen_xmin_outer_core, reclen_xmax_outer_core, &
             reclen_ymin_outer_core, reclen_ymax_outer_core
  integer :: reclen_zmin

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_OUTER_CORE_ADJOINT) :: &
    vector_accel_outer_core,vector_displ_outer_core,b_vector_displ_outer_core

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

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: rho_kl_outer_core,alpha_kl_outer_core

  ! check for deviatoric kernel for outer core region
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: beta_kl_outer_core
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

module specfem_par_movie

! parameter module for movies/shakemovies

  use constants_solver

  implicit none

  ! to save movie frames
  integer :: nmovie_points,NIT

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_ux,store_val_uy,store_val_uz
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: store_val_ux_all,store_val_uy_all,store_val_uz_all

  ! to save movie volume
  double precision :: scalingval

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: &
    muvstore_crust_mantle_3dmovie

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: &
    Iepsilondev_xx_crust_mantle,Iepsilondev_yy_crust_mantle,Iepsilondev_xy_crust_mantle, &
    Iepsilondev_xz_crust_mantle,Iepsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: &
    Ieps_trace_over_3_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: nu_3dmovie

  integer :: npoints_3dmovie,nspecel_3dmovie
  integer, dimension(NGLOB_CRUST_MANTLE_3DMOVIE) :: num_ibool_3dmovie

  logical, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: mask_3dmovie
  logical, dimension(NGLOB_CRUST_MANTLE_3DMOVIE) :: mask_ibool

  ! outer core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE) :: &
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
  ! multi-MPI processes, gather data arrays on master
  real,dimension(:),allocatable :: vtkdata_all
  integer,dimension(:),allocatable :: vtkdata_points_all
  integer,dimension(:),allocatable :: vtkdata_offset_all
  integer :: vtkdata_numpoints_all
  real :: vtkdata_source_x,vtkdata_source_y,vtkdata_source_z

end module specfem_par_movie
