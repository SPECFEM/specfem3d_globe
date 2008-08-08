
! allocate these automatic arrays in the memory stack to avoid memory fragmentation with "allocate()"

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll,wxgll
  double precision, dimension(NGLLY) :: yigll,wygll
  double precision, dimension(NGLLZ) :: zigll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(NGNOD,NGLLX,NGLLY,NGLLZ) :: shape3D
  double precision, dimension(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ) :: dershape3D

! 2D shape functions and their derivatives
  double precision shape2D_x(NGNOD2D,NGLLY,NGLLZ)
  double precision shape2D_y(NGNOD2D,NGLLX,NGLLZ)
  double precision shape2D_bottom(NGNOD2D,NGLLX,NGLLY)
  double precision shape2D_top(NGNOD2D,NGLLX,NGLLY)

  double precision dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ)
  double precision dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)
  double precision dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: rhostore_local,kappavstore_local

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

  logical, dimension(6,NSPEC_CRUST_MANTLE) :: iboun

! very large arrays used for the sorting routine
  integer, dimension(NSPEC_CRUST_MANTLE * NGLLX * NGLLY * NGLLZ) :: locval
  logical, dimension(NSPEC_CRUST_MANTLE * NGLLX * NGLLY * NGLLZ) :: ifseg
  double precision, dimension(NSPEC_CRUST_MANTLE * NGLLX * NGLLY * NGLLZ) :: xp,yp,zp

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load

  integer, dimension(NGLOB_CRUST_MANTLE) :: mask_ibool
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: copy_ibool_ori

  logical, dimension(2,NSPEC_CRUST_MANTLE) :: iMPIcut_xi,iMPIcut_eta

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY) :: rho_vp,rho_vs

  double precision, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: Qmu_store
  double precision, dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: tau_e_store

  integer ibelm_moho_top(NSPEC2D_MOHO)
  integer ibelm_moho_bot(NSPEC2D_MOHO)
  integer ibelm_400_top(NSPEC2D_400)
  integer ibelm_400_bot(NSPEC2D_400)
  integer ibelm_670_top(NSPEC2D_670)
  integer ibelm_670_bot(NSPEC2D_670)

  real(kind=CUSTOM_REAL) normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO)
  real(kind=CUSTOM_REAL) normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400)
  real(kind=CUSTOM_REAL) normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670)

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

  integer, dimension(NGLOB2DMAX_XY_VAL_CM,NUMFACES_SHARED) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_VAL_OC,NUMFACES_SHARED) :: iboolfaces_outer_core
  integer, dimension(NGLOB2DMAX_XY_VAL_IC,NUMFACES_SHARED) :: iboolfaces_inner_core

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

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_OC) :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_OC) :: iboolleft_eta_outer_core,iboolright_eta_outer_core
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

!! DK DK added this for the merged version
!---- arrays to assemble between chunks

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: imsg_type,iprocfrom_faces,iprocto_faces

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

