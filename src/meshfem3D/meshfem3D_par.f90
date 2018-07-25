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

!
!-------------------------------------------------------------------------------------------------
!

  module meshfem3D_models_par

!---
!
! ADD YOUR MODEL HERE
!
!---

  use constants

  use shared_input_parameters, only: &
    ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS, &
    ATTENUATION,USE_FULL_TISO_MANTLE

  use shared_compute_parameters, only: &
    REFERENCE_1D_MODEL,REFERENCE_CRUSTAL_MODEL,THREE_D_MODEL, &
    HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY, &
    ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE, &
    ATTENUATION_3D, &
    ANISOTROPIC_INNER_CORE, &
    CEM_REQUEST, CEM_ACCEPT

  implicit none

  ! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision :: min_period, max_period
    double precision :: QT_c_source        ! Source Frequency
    double precision, dimension(:), allocatable   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), allocatable   :: QrDisc             ! Discontinuities Defined
    double precision, dimension(:), allocatable   :: Qr                 ! Radius
    double precision, dimension(:), allocatable   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), allocatable :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), allocatable   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), allocatable :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), allocatable   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), allocatable            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), allocatable            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), allocatable            :: interval_Q                 ! Steps
    integer :: Qn                 ! Number of points
    integer :: dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables
  type (model_attenuation_variables) :: AM_V
  ! model_attenuation_variables

  ! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), allocatable :: tau_e_storage
    double precision, dimension(:), allocatable   :: Qmu_storage
    integer :: Q_resolution
    integer :: Q_max
  end type model_attenuation_storage_var
  type (model_attenuation_storage_var) :: AM_S
  ! model_attenuation_storage_var

  ! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision :: Q  ! Q     = Desired Value of Attenuation or Q
    double precision :: iQ ! iQ    = 1/Q
    double precision, dimension(:), allocatable ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), allocatable :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer :: nf          ! nf    = Number of Frequencies
    integer :: nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables
  ! attenuation_simplex_variables

  ! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    logical :: MODEL_GLL
    logical,dimension(3) :: dummy_pad ! padding 3 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) :: MGLL_V

  ! bathymetry and topography: use integer array to store values
  integer, dimension(:,:),allocatable :: ibathy_topo

  ! for ellipticity
  double precision,dimension(NR) :: rspl,espl,espl2
  integer :: nspl

  end module meshfem3D_models_par


!
!-------------------------------------------------------------------------------------------------
!


  module meshfem3D_par

! main parameter module for specfem simulations

  use constants

  use shared_parameters

  implicit none

  ! correct number of spectral elements in each block depending on chunk type
  integer :: npointot

  ! check area and volume of the final mesh
  double precision :: volume_total

  ! check Earth mass and Earth center of mass computed in the final mesh
  double precision :: Earth_mass_total
  double precision :: Earth_center_of_mass_x_total,Earth_center_of_mass_y_total,Earth_center_of_mass_z_total
  double precision :: distance_to_center_in_km

  ! arrays containing the positions of the observation points in non-dimensionalized value for gravity integrals
  double precision, dimension(NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX) :: x_observation,y_observation,z_observation

  ! daniel: compilers don't won't like equivalence statements, as they hinder internal optimizations
  !         let's go without it and do the 1D-loops with FORCE_VECTORIZATION when no bounds-checking is done,
  !         i.e., when compiled with run flags (already done in other code parts)
  !
  ! left here for referencing original way:
  ! the 1D equivalenced versions are for the FORCE_VECTORIZATION version of the loops
  !double precision, dimension(NTOTAL_OBSERVATION) :: x_observation1D,y_observation1D,z_observation1D
  !equivalence(x_observation,x_observation1D)
  !equivalence(y_observation,y_observation1D)
  !equivalence(z_observation,z_observation1D)

  double precision, dimension(NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX) :: lon_observation,lat_observation

  ! arrays containing the computed fields for gravity integrals at the observation points
  double precision, dimension(NX_OBSERVATION,NY_OBSERVATION,NCHUNKS_MAX) :: g_x,g_y,g_z,G_xx,G_yy,G_zz,G_xy,G_xz,G_yz, &
                                                                            temporary_array_for_sum

  ! for loop on all the slices
  integer :: iregion_code
  integer :: iproc_xi,iproc_eta,ichunk

  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix
  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! for some statistics for the mesh
  integer :: numelem_crust_mantle,numelem_outer_core,numelem_inner_core
  integer :: numelem_total

  ! timer MPI
  double precision :: time_start

  ! addressing for all the slices
  integer, dimension(:), allocatable :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer, dimension(:,:,:), allocatable :: addressing

  ! memory size of all the static arrays
  double precision :: static_memory_size

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
             NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
             NSPEC_INNER_CORE_ATTENUATION, &
             NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
             NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
             NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
             NSPEC_CRUST_MANTLE_ADJOINT, &
             NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
             NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
             NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
             NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
             NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

  ! this for the different corners of the slice (which are different if the superbrick is cut)
  ! 1 : xi_min, eta_min
  ! 2 : xi_max, eta_min
  ! 3 : xi_max, eta_max
  ! 4 : xi_min, eta_max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER

  ! 1 -> min, 2 -> max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE

  ! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

  ! temporary arrays with mesh only on global points
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: xstore_glob,ystore_glob,zstore_glob

  ! parameters needed to store the radii of the grid points
  ! in the spherically symmetric Earth
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool

  ! this for non blocking MPI
  logical, dimension(:), allocatable :: is_on_a_slice_edge

  ! mesh global time step
  real(kind=CUSTOM_REAL) :: dt_max_glob
  ! mesh minimum period resolved
  real(kind=CUSTOM_REAL) :: pmax_glob

  ! number of spectral elements (in current region)
  integer :: nspec

  ! number of global GLL points (in current region)
  integer :: nglob

  end module meshfem3D_par

!
!-------------------------------------------------------------------------------------------------
!

  module regions_mesh_par

  use constants, only: NGLLX,NGLLY,NGLLZ,NGNOD,NGNOD2D,NDIM,NDIM2D

  implicit none

  ! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

  ! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

  ! 3D shape functions and their derivatives
  double precision, dimension(NGNOD,NGLLX,NGLLY,NGLLZ) :: shape3D
  double precision, dimension(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ) :: dershape3D

  ! 2D shape functions and their derivatives
  double precision, dimension(NGNOD2D,NGLLY,NGLLZ) :: shape2D_x
  double precision, dimension(NGNOD2D,NGLLX,NGLLZ) :: shape2D_y
  double precision, dimension(NGNOD2D,NGLLX,NGLLY) :: shape2D_bottom,shape2D_top
  double precision, dimension(NDIM2D,NGNOD2D,NGLLY,NGLLZ) :: dershape2D_x
  double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLZ) :: dershape2D_y
  double precision, dimension(NDIM2D,NGNOD2D,NGLLX,NGLLY) :: dershape2D_bottom,dershape2D_top

  end module regions_mesh_par

!
!-------------------------------------------------------------------------------------------------
!

  module regions_mesh_par2

  use constants, only: CUSTOM_REAL,N_SLS,MAX_STRING_LEN

  implicit none

  integer :: nspec_stacey,nspec_actually,nspec_att

  integer :: ifirst_region,ilast_region
  integer, dimension(:), allocatable :: perm_layer

  ! for model density and anisotropy
  integer :: nspec_ani
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,dvpstore, &
    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! boundary locator
  logical, dimension(:,:), allocatable :: iboun

  ! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassy,rmassz
  integer :: nglob_xy

  ! mass matrices for backward simulation when ROTATION is .true.
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_rmassx,b_rmassy

  ! mass matrix and bathymetry for ocean load
  integer :: nglob_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

  ! number of elements on the boundaries
  integer :: nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  ! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax, &
    ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

  ! 2-D Jacobians and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_xmin,jacobian2D_xmax, &
    jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom,jacobian2D_top

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top

  ! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

  ! Stacey, indices for Clayton-Engquist absorbing conditions
  integer, dimension(:,:), allocatable :: nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

  ! attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),   allocatable :: Qmu_store
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: tau_e_store
  double precision, dimension(N_SLS) :: tau_s
  double precision :: T_c_source

  ! element layers
  integer :: NUMBER_OF_MESH_LAYERS,layer_shift, &
    first_layer_aniso,last_layer_aniso
  logical :: USE_ONE_LAYER_SB

  ! layer stretching
  double precision, dimension(:,:), allocatable :: stretch_tab

  ! Boundary kernel Mesh
  integer :: NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho
  integer, dimension(:), allocatable :: ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_moho,normal_400,normal_670
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_moho,jacobian2D_400,jacobian2D_670

  integer :: ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot
  double precision :: r_moho,r_400,r_670

  ! flags for transverse isotropic elements
  logical, dimension(:), allocatable :: ispec_is_tiso

  ! name of the database file
  character(len=MAX_STRING_LEN) :: prname, prname_adios

  end module regions_mesh_par2

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_interfaces_par

  use constants, only: &
    CUSTOM_REAL,NDIM,IMAIN, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE, &
    NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  ! indirect addressing for each message for faces and corners of the chunks
  ! a given slice can belong to at most one corner and at most two faces

  ! number of faces between chunks
  integer :: NUMMSGS_FACES

  ! number of corners between chunks
  integer :: NCORNERSCHUNKS

  ! number of message types
  integer :: NUM_MSG_TYPES

  !-----------------------------------------------------------------
  ! assembly
  !-----------------------------------------------------------------

  ! ---- arrays to assemble between chunks
  ! communication pattern for faces between chunks
  integer, dimension(:),allocatable :: iprocfrom_faces,iprocto_faces,imsg_type
  ! communication pattern for corners between chunks
  integer, dimension(:),allocatable :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner

  ! chunk faces
  integer, dimension(:,:),allocatable :: iboolfaces
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces
  integer :: NGLOB2DMAX_XY

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi,iboolright_xi
  integer, dimension(:),allocatable :: iboolleft_eta,iboolright_eta
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_all,npoin2D_eta_all
  integer :: npoin2D_xi,npoin2D_eta

  ! 1-D addressing
  integer :: NGLOB1D_RADIAL_MAX
  integer,dimension(:),allocatable :: ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                                    ibool1D_leftxi_righteta,ibool1D_rightxi_righteta

  double precision,dimension(:,:),allocatable :: xyz1D_leftxi_lefteta,xyz1D_rightxi_lefteta, &
                                    xyz1D_leftxi_righteta,xyz1D_rightxi_righteta

  ! this for non blocking MPI

  ! buffers for send and receive between faces of the slices and the chunks
  ! we use the same buffers to assemble scalars and vectors because vectors are
  ! always three times bigger and therefore scalars can use the first part
  ! of the vector buffer in memory even if it has an additional index here
  integer :: npoin2D_max_all_CM_IC

  ! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar

  ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(:,:),allocatable :: &
     buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector

  end module MPI_interfaces_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_crust_mantle_par

  use constants, only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  ! collected MPI interfaces
  !--------------------------------------
  ! MPI crust/mantle mesh
  !--------------------------------------
  integer :: num_interfaces_crust_mantle
  integer :: max_nibool_interfaces_cm
  integer, dimension(:), allocatable :: my_neighbors_crust_mantle,nibool_interfaces_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle

  !--------------------------------------
  ! crust mantle
  !--------------------------------------
  integer :: NSPEC_CRUST_MANTLE
  integer :: NGLOB_CRUST_MANTLE

  integer :: NGLOB1D_RADIAL_CM
  integer :: NGLOB2DMAX_XMIN_XMAX_CM
  integer :: NGLOB2DMAX_YMIN_YMAX_CM
  integer :: NSPEC2DMAX_XMIN_XMAX_CM
  integer :: NSPEC2DMAX_YMIN_YMAX_CM
  integer :: NSPEC2D_BOTTOM_CM
  integer :: NSPEC2D_TOP_CM

  ! assembly
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_crust_mantle
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner_crust_mantle

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(:),allocatable :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle

  integer, dimension(:,:),allocatable :: iboolfaces_crust_mantle

  ! inner / outer elements crust/mantle region
  integer :: num_phase_ispec_crust_mantle
  integer :: nspec_inner_crust_mantle,nspec_outer_crust_mantle
  integer, dimension(:,:), allocatable :: phase_ispec_inner_crust_mantle

  ! mesh coloring
  integer :: num_colors_outer_crust_mantle,num_colors_inner_crust_mantle
  integer,dimension(:),allocatable :: num_elem_colors_crust_mantle

  end module MPI_crust_mantle_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_inner_core_par

  use constants, only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  !--------------------------------------
  ! MPI inner core mesh
  !--------------------------------------
  integer :: num_interfaces_inner_core
  integer :: max_nibool_interfaces_ic
  integer, dimension(:), allocatable :: my_neighbors_inner_core,nibool_interfaces_inner_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core

  !--------------------------------------
  ! inner core
  !--------------------------------------
  integer :: NSPEC_INNER_CORE
  integer :: NGLOB_INNER_CORE

  integer :: NGLOB1D_RADIAL_IC
  integer :: NGLOB2DMAX_XMIN_XMAX_IC
  integer :: NGLOB2DMAX_YMIN_YMAX_IC
  integer :: NSPEC2DMAX_XMIN_XMAX_IC
  integer :: NSPEC2DMAX_YMIN_YMAX_IC
  integer :: NSPEC2D_BOTTOM_IC
  integer :: NSPEC2D_TOP_IC

  ! for matching with central cube in inner core
  integer, dimension(:), allocatable :: sender_from_slices_to_cube
  integer, dimension(:,:), allocatable :: ibool_central_cube
  double precision, dimension(:,:), allocatable :: buffer_slices,buffer_slices2
  double precision, dimension(:,:,:), allocatable :: buffer_all_cube_from_slices
  integer :: nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube, &
    npoin2D_cube_from_slices,receiver_cube_from_slices

  ! bottom inner core / top central cube
  integer, dimension(:),allocatable :: ibelm_bottom_inner_core

  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_inner_core
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_inner_core,npoin2D_eta_inner_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner_inner_core

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(:),allocatable :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer, dimension(:,:),allocatable :: iboolfaces_inner_core

  ! inner / outer elements inner core region
  integer :: num_phase_ispec_inner_core
  integer :: nspec_inner_inner_core,nspec_outer_inner_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_inner_core

  ! mesh coloring
  integer :: num_colors_outer_inner_core,num_colors_inner_inner_core
  integer,dimension(:),allocatable :: num_elem_colors_inner_core

  end module MPI_inner_core_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_outer_core_par

  use constants, only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR

  implicit none

  !--------------------------------------
  ! MPI outer core mesh
  !--------------------------------------
  integer :: num_interfaces_outer_core
  integer :: max_nibool_interfaces_oc
  integer, dimension(:), allocatable :: my_neighbors_outer_core,nibool_interfaces_outer_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core

  !--------------------------------------
  ! outer core
  !--------------------------------------
  integer :: NSPEC_OUTER_CORE
  integer :: NGLOB_OUTER_CORE

  integer :: NGLOB1D_RADIAL_OC
  integer :: NGLOB2DMAX_XMIN_XMAX_OC
  integer :: NGLOB2DMAX_YMIN_YMAX_OC
  integer :: NSPEC2DMAX_XMIN_XMAX_OC
  integer :: NSPEC2DMAX_YMIN_YMAX_OC
  integer :: NSPEC2D_BOTTOM_OC
  integer :: NSPEC2D_TOP_OC

  ! assembly
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_outer_core
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_outer_core,npoin2D_eta_outer_core

  ! indirect addressing for each corner of the chunks
  integer, dimension(:,:),allocatable :: iboolcorner_outer_core

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(:),allocatable :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(:),allocatable :: iboolleft_eta_outer_core,iboolright_eta_outer_core

  integer, dimension(:,:),allocatable :: iboolfaces_outer_core

  ! inner / outer elements outer core region
  integer :: num_phase_ispec_outer_core
  integer :: nspec_inner_outer_core,nspec_outer_outer_core
  integer, dimension(:,:), allocatable :: phase_ispec_inner_outer_core

  ! mesh coloring
  integer :: num_colors_outer_outer_core,num_colors_inner_outer_core
  integer,dimension(:),allocatable :: num_elem_colors_outer_core

  end module MPI_outer_core_par
