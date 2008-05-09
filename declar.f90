
!! DK DK added this temporarily
!!!!!!!!!!!!!  integer, dimension(NSPEC_CRUST_MANTLE) :: perm,invperm
!! DK DK suppressed for now because useless when CUTHILL-MCKEE is off
  integer, dimension(1) :: perm,invperm

!! DK DK added this for merged version
!! DK DK stored in single precision for merged version, check if it precise enough (probably yes)
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
          xelm_store_crust_mantle,yelm_store_crust_mantle,zelm_store_crust_mantle, &
          xelm_store_outer_core,yelm_store_outer_core,zelm_store_outer_core, &
          xelm_store_inner_core,yelm_store_inner_core,zelm_store_inner_core

  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:,:) :: R_memory_crust_mantle
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:,:) :: R_memory_inner_core
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:,:) :: epsilondev_inner_core
  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: rmass_ocean_load
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:) :: &
        xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
        etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
        gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: &
     displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:) :: &
        xix_outer_core,xiy_outer_core,xiz_outer_core,&
        etax_outer_core,etay_outer_core,etaz_outer_core, &
        gammax_outer_core,gammay_outer_core,gammaz_outer_core
  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: xstore_outer_core,ystore_outer_core,zstore_outer_core
  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: displ_outer_core,veloc_outer_core,accel_outer_core
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:,:,:) :: &
        xix_inner_core,xiy_inner_core,xiz_inner_core,&
        etax_inner_core,etay_inner_core,etaz_inner_core, &
        gammax_inner_core,gammay_inner_core,gammaz_inner_core

  real(kind=CUSTOM_REAL), allocatable, dimension(:) :: xstore_inner_core,ystore_inner_core,zstore_inner_core
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: displ_inner_core,veloc_inner_core,accel_inner_core

!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead
!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead
!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead
!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead
!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead
!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead
!!!!!!!!!!!!!!!! DK DK for merged version, all the arrays below are allocated statically instead

  integer npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer npoin2D_faces_outer_core(NUMFACES_SHARED)
  integer npoin2D_faces_inner_core(NUMFACES_SHARED)

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
             npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
             npoin2D_xi_inner_core,npoin2D_eta_inner_core

! number of elements on the boundaries
  integer :: nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle,nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle
  integer :: nspec2D_xmin_outer_core,nspec2D_xmax_outer_core,nspec2D_ymin_outer_core,nspec2D_ymax_outer_core
  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core,nspec2D_ymin_inner_core,nspec2D_ymax_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmass_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: rmass_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: rmass_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: kappavstore_crust_mantle,muvstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
    kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore_inner_core,muvstore_inner_core

!! DK DK added this for the merged version
!! DK DK these arrays are useless in the solver and will therefore be allocated with a dummy size of 1
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: kappahstore_inner_core,muhstore_inner_core,eta_anisostore_inner_core
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: kappavstore_outer_core,muvstore_outer_core
  real(kind=CUSTOM_REAL), dimension(1,1,1,1) :: kappahstore_outer_core,muhstore_outer_core,eta_anisostore_outer_core

! 2-D jacobians and normals
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_IC) :: jacobian2D_xmin_inner_core,jacobian2D_xmax_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_IC) :: jacobian2D_ymin_inner_core,jacobian2D_ymax_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_IC) :: jacobian2D_bottom_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_IC) :: jacobian2D_top_inner_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_IC) :: normal_xmin_inner_core,normal_xmax_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_IC) :: normal_ymin_inner_core,normal_ymax_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_IC) :: normal_bottom_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_IC) :: normal_top_inner_core

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_CM) :: ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_CM) :: ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM) :: ibelm_bottom_crust_mantle
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_CM) :: jacobian2D_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_CM) :: jacobian2D_top_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: normal_xmin_crust_mantle,normal_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2DMAX_YMIN_YMAX_CM) :: normal_ymin_crust_mantle,normal_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_CM) :: normal_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_CM) :: normal_top_crust_mantle

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_OC) :: ibelm_xmin_outer_core,ibelm_xmax_outer_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_OC) :: ibelm_ymin_outer_core,ibelm_ymax_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core
  integer, dimension(NSPEC2D_TOP_OC) :: ibelm_top_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: normal_xmin_outer_core,normal_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: normal_ymin_outer_core,normal_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC) :: normal_top_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC) :: jacobian2D_top_outer_core

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_IC) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_IC) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_IC) :: ibelm_bottom_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_crust_mantle,iboolfaces_outer_core,iboolfaces_inner_core

!!!!  real(kind=CUSTOM_REAL), dimension(NGLOB2DMAX_XY_VAL) :: buffer_send_faces_scalar,buffer_received_faces_scalar

!!!!  real(kind=CUSTOM_REAL), dimension(NDIM,2*NGLOB2DMAX_XY_VAL) :: buffer_send_faces_vector,buffer_received_faces_vector

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling_crust_mantle

!! DK DK this array is useless in the solver and is therefore allocated with a dummy size of 1
  integer, dimension(1) :: idoubling_outer_core

  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  integer, dimension(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_OC,NUMCORNERS_SHARED) :: iboolcorner_outer_core
  integer, dimension(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED) :: iboolcorner_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL_CM) :: buffer_send_chunkcorners_scalar,buffer_recv_chunkcorners_scalar

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC) :: &
    buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_OC) :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_OC) :: iboolleft_eta_outer_core,iboolright_eta_outer_core
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  double precision, dimension(ATT1,ATT2,ATT3,ATT4) :: omsb_crust_mantle_dble, factor_scale_crust_mantle_dble

  double precision, dimension(ATT1,ATT2,ATT3,ATT5) :: omsb_inner_core_dble, factor_scale_inner_core_dble

  real(kind=CUSTOM_REAL), dimension(ATT1,ATT2,ATT3,ATT4) :: one_minus_sum_beta_crust_mantle, factor_scale_crust_mantle

  real(kind=CUSTOM_REAL), dimension(ATT1,ATT2,ATT3,ATT5) :: one_minus_sum_beta_inner_core, factor_scale_inner_core

  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1,ATT2,ATT3,ATT4) :: factor_common_crust_mantle

  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1,ATT2,ATT3,ATT5) :: factor_common_inner_core

  double precision, dimension(N_SLS,ATT1,ATT2,ATT3,ATT4) :: factor_common_crust_mantle_dble

  double precision, dimension(N_SLS,ATT1,ATT2,ATT3,ATT5) :: factor_common_inner_core_dble

