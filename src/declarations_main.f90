
!! DK DK for merged version, all the arrays below are allocated statically instead

  integer npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer npoin2D_faces_outer_core(NUMFACES_SHARED)
  integer npoin2D_faces_inner_core(NUMFACES_SHARED)

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
             npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
             npoin2D_xi_inner_core,npoin2D_eta_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_CM) :: normal_top_crust_mantle
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle

! number of elements on the boundaries
  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core,nspec2D_ymin_inner_core,nspec2D_ymax_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmass_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: rmass_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: rmass_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: kappavstore_crust_mantle,muvstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
    kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore_inner_core,muvstore_inner_core

! 2-D jacobians and normals
  integer, dimension(NSPEC2D_BOTTOM_CM) :: ibelm_bottom_crust_mantle

  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core
  integer, dimension(NSPEC2D_TOP_OC) :: ibelm_top_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC) :: normal_top_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC) :: jacobian2D_top_outer_core

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_IC) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_IC) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_IC) :: ibelm_bottom_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

  integer, dimension(NGLOB2DMAX_XY_VAL_CM,NUMFACES_SHARED) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_VAL_OC,NUMFACES_SHARED) :: iboolfaces_outer_core
  integer, dimension(NGLOB2DMAX_XY_VAL_IC,NUMFACES_SHARED) :: iboolfaces_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling_crust_mantle

  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  integer, dimension(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_OC,NUMCORNERS_SHARED) :: iboolcorner_outer_core
  integer, dimension(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED) :: iboolcorner_inner_core

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_OC) :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_OC) :: iboolleft_eta_outer_core,iboolright_eta_outer_core
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer :: npoin2D_max_all,NDIM_smaller_buffers

! receiver information
  integer :: nrec,ios
  character(len=150) :: STATIONS,rec_filename,dummystring

!---- arrays to assemble between chunks

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: imsg_type,iprocfrom_faces,iprocto_faces

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(N_SLS)        :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: Qs                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

  type (attenuation_model_variables) AM_V
! attenuation_model_variables

