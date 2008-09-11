!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
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

  subroutine create_regions_mesh(iregion_code,ibool,idoubling,xstore,ystore,zstore,rmins,rmaxs,iproc_xi,iproc_eta,ichunk,nspec, &
  nspec_tiso,volume_local,area_local_bottom,area_local_top,nspl,rspl,espl,espl2,nglob_theor,npointot,NEX_XI,NEX_PER_PROC_XI, &
  NEX_PER_PROC_ETA,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP,ELLIPTICITY,TOPOGRAPHY, &
  TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,NPROC_XI,NPROC_ETA, &
  myrank,OCEANS,ibathy_topo,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,ATTENUATION,ATTENUATION_3D,NCHUNKS, &
  INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL,R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO, &
  RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN,ner,ratio_sampling_array,doubling_index,r_bottom,r_top, &
  this_layer_has_a_doubling,CASE_3D,AMM_V,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V,AM_S,AS_V, &
  numker,numhpa,numcof,ihpa,lmax,nylm,lmxhpa,itypehpa,ihpakern,numcoe,ivarkern,nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
  coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ratio_divide_central_cube,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
  offset_proc_xi,offset_proc_eta,nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax,ibelm_xmin,ibelm_xmax,ibelm_ymin, &
  ibelm_ymax,ibelm_bottom,ibelm_top,jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom, &
  jacobian2D_top,normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top,kappavstore,kappahstore,muvstore, &
  muhstore,eta_anisostore,rmass,xelm_store,yelm_store,zelm_store,xigll,wxgll,yigll,wygll,zigll,wzgll,shape3D,dershape3D, &
  shape2D_x,shape2D_y,shape2D_bottom,shape2D_top,dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top,rhostore_local, &
  kappavstore_local,c11store,c12store,c13store,c14store,c15store,c16store,c22store,c23store,c24store,c25store,c26store, &
  c33store,c34store,c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store,iboun,locval,ifseg, &
  rmass_ocean_load,mask_ibool,copy_ibool_ori,iMPIcut_xi,iMPIcut_eta, &
#ifdef USE_MPI
  NGLOB1D_RADIAL_MAX,NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE,NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,NGLOB2DMAX_XMIN_XMAX, &
  NGLOB2DMAX_YMIN_YMAX,npoin2D_xi,npoin2D_eta,iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta,ibool1D_leftxi_lefteta, &
  ibool1D_rightxi_lefteta,ibool1D_leftxi_righteta,ibool1D_rightxi_righteta,xread1D_leftxi_lefteta,xread1D_rightxi_lefteta, &
  xread1D_leftxi_righteta,xread1D_rightxi_righteta,yread1D_leftxi_lefteta,yread1D_rightxi_lefteta,yread1D_leftxi_righteta, &
  yread1D_rightxi_righteta,zread1D_leftxi_lefteta,zread1D_rightxi_lefteta,zread1D_leftxi_righteta,zread1D_rightxi_righteta, &
#endif
  rho_vp,rho_vs,Qmu_store,tau_e_store,ifirst_layer_aniso,ilast_layer_aniso,SAVE_MESH_FILES)

! create the different regions of the mesh

  implicit none

#ifdef USE_MPI
  include "mpif.h"
#endif
  include "constants.h"

!! DK DK for the merged version
! include values created by the mesher
  include "values_from_mesher.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! DK DK to debug the second sorting routine
  logical, parameter :: DEBUG = .true.

!! DK DK added this for merged version
#ifdef USE_MPI
  integer :: npoin2D_xi,npoin2D_eta
#endif

! mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob_theor) :: rmass

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

! the jacobian
  real(kind=CUSTOM_REAL) jacobianl

!! DK DK changed this for merged version: made it local
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xixstore,xiystore,xizstore, &
        etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

!! DK DK added this for merged version
  logical :: add_contrib_this_element

#ifdef USE_MPI
!! DK DK for merged version
  integer :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer, dimension(NGLOB2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NGLOB2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

!! DK DK added this for merged version
  integer :: NGLOB1D_RADIAL_MAX
  double precision, dimension(NGLOB1D_RADIAL_MAX) :: &
  xread1D_leftxi_lefteta, xread1D_rightxi_lefteta, xread1D_leftxi_righteta, xread1D_rightxi_righteta, &
  yread1D_leftxi_lefteta, yread1D_rightxi_lefteta, yread1D_leftxi_righteta, yread1D_rightxi_righteta, &
  zread1D_leftxi_lefteta, zread1D_rightxi_lefteta, zread1D_leftxi_righteta, zread1D_rightxi_righteta

  integer ibool1D_leftxi_lefteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D_rightxi_lefteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D_leftxi_righteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D_rightxi_righteta(NGLOB1D_RADIAL_MAX)

! this to cut the doubling brick
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE
#endif
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer :: step_mult,offset_proc_xi,offset_proc_eta
  integer :: case_xi,case_eta,subblock_num

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_layer_has_a_doubling

  integer :: ignod,ner_without_doubling,ispec_superbrick,ilayer,ilayer_loop,ix_elem,iy_elem,iz_elem, &
               ifirst_layer,ilast_layer,ratio_divide_central_cube

! allocate this automatic array in the memory stack to avoid memory fragmentation with "allocate()"
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: permutation_layer

! mesh doubling superbrick
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick

  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick

! aniso_mantle_model_variables
  type aniso_mantle_model_variables
    sequence
    double precision beta(AMM_V_DIM1,AMM_V_DIM2,AMM_V_DIM3,AMM_V_DIM4)
    double precision pro(47)
    integer npar1
  end type aniso_mantle_model_variables

  type (aniso_mantle_model_variables) AMM_V
! aniso_mantle_model_variables

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(N_SLS)        :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
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

! model_1066a_variables
  type model_1066a_variables
    sequence
      double precision, dimension(NR_1066A) :: radius_1066a
      double precision, dimension(NR_1066A) :: density_1066a
      double precision, dimension(NR_1066A) :: vp_1066a
      double precision, dimension(NR_1066A) :: vs_1066a
      double precision, dimension(NR_1066A) :: Qkappa_1066a
      double precision, dimension(NR_1066A) :: Qmu_1066a
  end type model_1066a_variables

  type (model_1066a_variables) M1066a_V
! model_1066a_variables

! model_ak135_variables
  type model_ak135_variables
    sequence
    double precision, dimension(NR_AK135) :: radius_ak135
    double precision, dimension(NR_AK135) :: density_ak135
    double precision, dimension(NR_AK135) :: vp_ak135
    double precision, dimension(NR_AK135) :: vs_ak135
    double precision, dimension(NR_AK135) :: Qkappa_ak135
    double precision, dimension(NR_AK135) :: Qmu_ak135
  end type model_ak135_variables

 type (model_ak135_variables) Mak135_V
! model_ak135_variables

! model_ref_variables
  type model_ref_variables
    sequence
     double precision, dimension(NR_REF) :: radius_ref
     double precision, dimension(NR_REF) :: density_ref
     double precision, dimension(NR_REF) :: vpv_ref
     double precision, dimension(NR_REF) :: vph_ref
     double precision, dimension(NR_REF) :: vsv_ref
     double precision, dimension(NR_REF) :: vsh_ref
     double precision, dimension(NR_REF) :: eta_ref
     double precision, dimension(NR_REF) :: Qkappa_ref
     double precision, dimension(NR_REF) :: Qmu_ref
  end type model_ref_variables

 type (model_ref_variables) Mref_V
! model_ref_variables

! sea1d_model_variables
  type sea1d_model_variables
    sequence
     double precision, dimension(NR_SEA1D) :: radius_sea1d
     double precision, dimension(NR_SEA1D) :: density_sea1d
     double precision, dimension(NR_SEA1D) :: vp_sea1d
     double precision, dimension(NR_SEA1D) :: vs_sea1d
     double precision, dimension(NR_SEA1D) :: Qkappa_sea1d
     double precision, dimension(NR_SEA1D) :: Qmu_sea1d
  end type sea1d_model_variables

  type (sea1d_model_variables) SEA1DM_V
! sea1d_model_variables

! three_d_mantle_model_variables
  type three_d_mantle_model_variables
    sequence
    double precision dvs_a(0:NK,0:NS,0:NS)
    double precision dvs_b(0:NK,0:NS,0:NS)
    double precision dvp_a(0:NK,0:NS,0:NS)
    double precision dvp_b(0:NK,0:NS,0:NS)
    double precision spknt(NK+1)
    double precision qq0(NK+1,NK+1)
    double precision qq(3,NK+1,NK+1)
  end type three_d_mantle_model_variables

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables

! jp3d_model_variables
  type jp3d_model_variables
    sequence
! vmod3d
  integer :: NPA
  integer :: NRA
  integer :: NHA
  integer :: NPB
  integer :: NRB
  integer :: NHB
  double precision :: PNA(MPA)
  double precision :: RNA(MRA)
  double precision :: HNA(MHA)
  double precision :: PNB(MPB)
  double precision :: RNB(MRB)
  double precision :: HNB(MHB)
  double precision :: VELAP(MPA,MRA,MHA)
  double precision :: VELBP(MPB,MRB,MHB)
! discon
  double precision :: PN(51)
  double precision :: RRN(63)
  double precision :: DEPA(51,63)
  double precision :: DEPB(51,63)
  double precision :: DEPC(51,63)
! locate
  integer :: IPLOCA(MKA)
  integer :: IRLOCA(MKA)
  integer :: IHLOCA(MKA)
  integer :: IPLOCB(MKB)
  integer :: IRLOCB(MKB)
  integer :: IHLOCB(MKB)
  double precision :: PLA
  double precision :: RLA
  double precision :: HLA
  double precision :: PLB
  double precision :: RLB
  double precision :: HLB
! weight
  integer :: IP
  integer :: JP
  integer :: KP
  integer :: IP1
  integer :: JP1
  integer :: KP1
  double precision :: WV(8)
! prhfd
  double precision :: P
  double precision :: R
  double precision :: H
  double precision :: PF
  double precision :: RF
  double precision :: HF
  double precision :: PF1
  double precision :: RF1
  double precision :: HF1
  double precision :: PD
  double precision :: RD
  double precision :: HD
! jpmodv
  double precision :: VP(29)
  double precision :: VS(29)
  double precision :: RA(29)
  double precision :: DEPJ(29)
  end type jp3d_model_variables

  type (jp3d_model_variables) JP3DM_V
! jp3d_model_variables

! sea99_s_model_variables
  type sea99_s_model_variables
    sequence
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    double precision :: sea99_vs(SEA99_VS_DIM1,SEA99_VS_DIM2,SEA99_VS_DIM3)
    double precision :: sea99_depth(100)
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
 end type sea99_s_model_variables

 type (sea99_s_model_variables) SEA99M_V
! sea99_s_model_variables

! crustal_model_variables
  type crustal_model_variables
    sequence
    real(kind=4) velocp(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) velocs(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) dens(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) thlr(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

! attenuation_model_storage
  type attenuation_model_storage
    sequence
    integer Q_resolution
    integer Q_max
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
  end type attenuation_model_storage

  type (attenuation_model_storage) AM_S
! attenuation_model_storage

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

! correct number of spectral elements in each block depending on chunk type
  integer nspec,nspec_tiso,nspec_stacey

  integer NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NCHUNKS,REFERENCE_1D_MODEL,THREE_D_MODEL

  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer NPROC_XI,NPROC_ETA

  integer npointot

  logical ELLIPTICITY,TOPOGRAPHY,SAVE_MESH_FILES
  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,OCEANS

  logical ATTENUATION,ATTENUATION_3D,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS

  double precision R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO, &
          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN

  character(len=150) errmsg,prname

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! meshing parameters
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! to define the central cube in the inner core
  integer nx_central_cube,ny_central_cube,nz_central_cube
  double precision radius_cube
  double precision :: xgrid_central_cube,ygrid_central_cube,zgrid_central_cube

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! auxiliary variables to generate the mesh
  integer ix,iy,iz

! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

! code for the four regions of the mesh
  integer iregion_code

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

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

!! DK DK added this for merged version
!! DK DK stored in single precision for merged version, check if it precise enough (probably yes)
  real(kind=CUSTOM_REAL), dimension(NGNOD,nspec) :: xelm_store,yelm_store,zelm_store

  double precision, dimension(NGNOD) :: xelm,yelm,zelm,offset_x,offset_y,offset_z

  integer idoubling(nspec)

! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision :: rmin,rmax

! for model density and anisotropy
  integer nspec_ani

!! DK DK changed this for the merged version
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: rhostore_local,kappavstore_local

!! DK DK added this for merged version
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! boundary locator
  logical, dimension(6,nspec) :: iboun

! proc numbers for MPI
  integer myrank

! check area and volume of the final mesh
  double precision weight
  double precision area_local_bottom,area_local_top
  double precision volume_local

! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(npointot) :: locval
  logical, dimension(npointot) :: ifseg

  integer :: nglob,nglob_theor,ier
#ifdef USE_MPI
  integer :: errorcode
#endif

! mass matrix and bathymetry for ocean load
  integer ix_oceans,iy_oceans,iz_oceans,ispec_oceans
  integer ispec2D_top_crust
  integer nglob_oceans
  double precision xval,yval,zval,rval,thetaval,phival
  double precision lat,lon,colat
  double precision elevation,height_oceans
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load

! mask to sort ibool
  integer :: inumber
  integer, dimension(nglob_theor) :: mask_ibool
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: copy_ibool_ori

! boundary parameters locator
  integer, dimension(NSPEC2DMAX_XMIN_XMAX) :: ibelm_xmin,ibelm_xmax
  integer, dimension(NSPEC2DMAX_YMIN_YMAX) :: ibelm_ymin,ibelm_ymax
  integer, dimension(NSPEC2D_BOTTOM) :: ibelm_bottom
  integer, dimension(NSPEC2D_TOP) :: ibelm_top

! MPI cut-planes parameters along xi and along eta
  logical, dimension(2,nspec) :: iMPIcut_xi,iMPIcut_eta

! Stacey indices for Clayton-Engquist absorbing conditions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY) :: rho_vp,rho_vs

! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer i,j,k,ia,ispec,iglobnum
  integer iproc_xi,iproc_eta,ichunk

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

! attenuation
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: Qmu_store
  double precision, dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: tau_e_store
  double precision, dimension(N_SLS) :: tau_s
  double precision  T_c_source

! **************
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
  logical :: USE_ONE_LAYER_SB,CASE_3D
  integer :: nspec_sb

  integer NUMBER_OF_MESH_LAYERS,layer_shift,ifirst_layer_aniso,ilast_layer_aniso
  double precision, dimension(:,:), allocatable :: stretch_tab

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)
  real(kind=4) vercof(maxker)
  real(kind=4) vercofd(maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=80) kerstr
  character(len=40) varstr(maxker)

! to perform two passes of the whole routine to be able to save memory
  integer :: ipass

! the height at which the central cube is cut
  integer :: nz_inf_limit

!! DK DK added this for the merged version
! 2-D jacobians and normals
  real(kind=CUSTOM_REAL) :: jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) :: jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) :: jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) :: jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) :: jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) :: jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP)

  real(kind=CUSTOM_REAL) :: normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) :: normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) :: normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) :: normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) :: normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) :: normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP)

! perform two passes of the whole routine to be able to save memory
  do ipass = 1,2

! attenuation
  if(ATTENUATION .and. ATTENUATION_3D) then
    T_c_source = AM_V%QT_c_source
    tau_s(:)   = AM_V%Qtau_s(:)
  else
    Qmu_store(1,1,1,1) = 0.0d0
    tau_e_store(:,1,1,1,1) = 0.0d0
  endif

! Stacey
  if(NCHUNKS /= 6) then
    nspec_stacey = nspec
  else
    nspec_stacey = 1
  endif

! anisotropy
  nspec_ani = 1
  if((ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) .or. &
     (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) nspec_ani = nspec

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY)

! define models 1066a and ak135 and ref
  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
    call define_model_1066a(CRUSTAL, M1066a_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
    call define_model_ak135(CRUSTAL, Mak135_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then
    call define_model_ref(Mref_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
    call define_model_sea1d(CRUSTAL, SEA1DM_V)
  endif

!------------------------------------------------------------------------

! create the shape of the corner nodes of a regular mesh element
  call hex_nodes(iaddx,iaddy,iaddz)

! reference element has size one here, not two
  iaddx(:) = iaddx(:) / 2
  iaddy(:) = iaddy(:) / 2
  iaddz(:) = iaddz(:) / 2

  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
    layer_shift = 0
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
    layer_shift = 1
  endif

  if (.not. ADD_4TH_DOUBLING) NUMBER_OF_MESH_LAYERS = NUMBER_OF_MESH_LAYERS - 1

! define the first and last layers that define this region
  if(iregion_code == IREGION_CRUST_MANTLE) then
    ifirst_layer = 1
    ilast_layer = 10 + layer_shift

  else if(iregion_code == IREGION_OUTER_CORE) then
    ifirst_layer = 11 + layer_shift
    ilast_layer = NUMBER_OF_MESH_LAYERS - 1

  else if(iregion_code == IREGION_INNER_CORE) then
    ifirst_layer = NUMBER_OF_MESH_LAYERS
    ilast_layer = NUMBER_OF_MESH_LAYERS

  else
    call exit_MPI(myrank,'incorrect region code detected')

  endif

  if(iregion_code == IREGION_CRUST_MANTLE) then

! create anisotropic (transversely isotropic) layers first to save memory when
! storing the anisotropic arrays
    ilayer = 0
    do ilayer_loop = ifirst_layer_aniso,ilast_layer_aniso
      ilayer = ilayer + 1
      permutation_layer(ilayer) = ilayer_loop
    enddo

! and then create all the isotropic layers
    do ilayer_loop = ifirst_layer,ilast_layer
      if(ilayer_loop < ifirst_layer_aniso .or. ilayer_loop > ilast_layer_aniso) then
        ilayer = ilayer + 1
        permutation_layer(ilayer) = ilayer_loop
      endif
    enddo

  else

! use identity permutation for regions that do not have transversely isotropic layer
    do ilayer_loop = ifirst_layer,ilast_layer
      permutation_layer(ilayer_loop) = ilayer_loop
    enddo

  endif

! initialize mesh arrays
!! DK DK merged version: we exclude the outer core because the doubling array is useless there and therefore not allocated
  if(iregion_code /= IREGION_OUTER_CORE) idoubling(:) = 0

  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  if(ipass == 1) ibool(:,:,:,:) = 0

! initialize boundary arrays
  iboun(:,:) = .false.
  iMPIcut_xi(:,:) = .false.
  iMPIcut_eta(:,:) = .false.

!! DK DK added this for merged version
! creating mass matrix in this slice (will be fully assembled in the solver)
  if(ipass == 2) rmass(:) = 0._CUSTOM_REAL

  if (.not. PATCH_FOR_GORDON_BELL .and. (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH)) then
    allocate(stretch_tab(2,ner(1)),STAT=ier )
    if (ier /= 0) then
      print *,"ABORTING can not allocate in create_regions_mesh ier=",ier
#ifdef USE_MPI
      call MPI_Abort(MPI_COMM_WORLD,errorcode,ier)
#else
      stop 'fatal error'
#endif
    endif

    call stretching_function(r_top(1),r_bottom(1),ner(1),stretch_tab)
  endif

! generate and count all the elements in this region of the mesh
  ispec = 0

! loop on all the layers in this region of the mesh
  do ilayer_loop = ifirst_layer,ilast_layer

    ilayer = permutation_layer(ilayer_loop)

! determine the radii that define the shell
  rmin = rmins(ilayer)
  rmax = rmaxs(ilayer)

    ner_without_doubling = ner(ilayer)

! if there is a doubling at the top of this region, we implement it in the last two layers of elements
! and therefore we suppress two layers of regular elements here
    USE_ONE_LAYER_SB = .false.
    if(this_layer_has_a_doubling(ilayer)) then
      if (ner(ilayer) == 1) then
        ner_without_doubling = ner_without_doubling - 1
        USE_ONE_LAYER_SB = .true.
      else
        ner_without_doubling = ner_without_doubling - 2
        USE_ONE_LAYER_SB = .false.
      endif
    endif

!----
!----   regular mesh elements
!----

! loop on all the elements
   do ix_elem = 1,NEX_PER_PROC_XI,ratio_sampling_array(ilayer)
   do iy_elem = 1,NEX_PER_PROC_ETA,ratio_sampling_array(ilayer)
   do iz_elem = 1,ner_without_doubling
! loop on all the corner nodes of this element
   do ignod = 1,NGNOD_EIGHT_CORNERS
! define topological coordinates of this mesh point
      offset_x(ignod) = (ix_elem - 1) + iaddx(ignod) * ratio_sampling_array(ilayer)
      offset_y(ignod) = (iy_elem - 1) + iaddy(ignod) * ratio_sampling_array(ilayer)
      if (.not. PATCH_FOR_GORDON_BELL .and. (ilayer == 1 .and. CASE_3D)) then
        offset_z(ignod) = iaddz(ignod)
      else
        offset_z(ignod) = (iz_elem - 1) + iaddz(ignod)
      endif
   enddo
     call add_missing_nodes(offset_x,offset_y,offset_z)

! compute the actual position of all the grid points of that element
  if (.not. PATCH_FOR_GORDON_BELL .and. (ilayer == 1 .and. CASE_3D .and. .not. SUPPRESS_CRUSTAL_MESH)) then
! crustal elements are stretched to be thinner in the upper crust than in lower crust in the 3D case
    call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               stretch_tab(1,ner_without_doubling-iz_elem+1),&
               stretch_tab(2,ner_without_doubling-iz_elem+1),1,ilayer,ichunk,rotation_matrix, &
               NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)
  else
     call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer,ichunk,rotation_matrix, &
               NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)
  endif
! add one spectral element to the list
     ispec = ispec + 1
     if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')

! new get_flag_boundaries
! xmin & xmax
  if (ix_elem == 1) then
      iMPIcut_xi(1,ispec) = .true.
      if (iproc_xi == 0) iboun(1,ispec)= .true.
  endif
  if (ix_elem == (NEX_PER_PROC_XI-ratio_sampling_array(ilayer)+1)) then
      iMPIcut_xi(2,ispec) = .true.
      if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= .true.
  endif
! ymin & ymax
  if (iy_elem == 1) then
      iMPIcut_eta(1,ispec) = .true.
      if (iproc_eta == 0) iboun(3,ispec)= .true.
  endif
  if (iy_elem == (NEX_PER_PROC_ETA-ratio_sampling_array(ilayer)+1)) then
      iMPIcut_eta(2,ispec) = .true.
      if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= .true.
  endif
! zmin & zmax
  if (iz_elem == ner(ilayer) .and. ilayer == ifirst_layer) then
      iboun(6,ispec)= .true.
  endif
  if (iz_elem == 1 .and. ilayer == ilast_layer) then    ! defined if no doubling in this layer
      iboun(5,ispec)= .true.
  endif

! define the doubling flag of this element
     if(iregion_code /= IREGION_OUTER_CORE) idoubling(ispec) = doubling_index(ilayer)

! compute several rheological and geometrical properties for this spectral element
     call compute_element_properties(ispec,iregion_code,idoubling, &
           xstore,ystore,zstore,nspec, &
           nspl,rspl,espl,espl2,ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           myrank,ibathy_topo,ATTENUATION,ATTENUATION_3D, &
           ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           RICB,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           xelm,yelm,zelm,shape3D,dershape3D,rmin,rmax,rhostore_local,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
!! DK DK added this for the merged version
           kappavstore_local, &
           xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           nspec_ani,nspec_stacey,Qmu_store,tau_e_store,tau_s,T_c_source,rho_vp,rho_vs,&
           AMM_V,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V,AM_S,AS_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ipass)

!! DK DK added this for merged version
    include "comp_mass_matrix_one_element.f90"
    include "store_xelm_yelm_zelm.f90"

! end of loop on all the regular elements
  enddo
  enddo
  enddo

!----
!----   mesh doubling elements
!----

! If there is a doubling at the top of this region, let us add these elements.
! The superbrick implements a symmetric four-to-two doubling and therefore replaces
! a basic regular block of 2 x 2 = 4 elements.
! We have imposed that NEX be a multiple of 16 therefore we know that we can always create
! these 2 x 2 blocks because NEX_PER_PROC_XI / ratio_sampling_array(ilayer) and
! NEX_PER_PROC_ETA / ratio_sampling_array(ilayer) are always divisible by 2.
    if(this_layer_has_a_doubling(ilayer)) then
      if (USE_ONE_LAYER_SB) then
        call define_superbrick_one_layer(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
        nspec_sb = NSPEC_SUPERBRICK_1L
        iz_elem = ner(ilayer)
        step_mult = 2
      else
        if(iregion_code==IREGION_OUTER_CORE .and. ilayer==ilast_layer .and. (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA)) then
          nspec_sb = NSPEC_DOUBLING_BASICBRICK
          step_mult = 1
        else
          call define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
          nspec_sb = NSPEC_DOUBLING_SUPERBRICK
          step_mult = 2
        endif
! the doubling is implemented in the last two radial elements
! therefore we start one element before the last one
        iz_elem = ner(ilayer) - 1
      endif

! loop on all the elements in the 2 x 2 blocks
      do ix_elem = 1,NEX_PER_PROC_XI,step_mult*ratio_sampling_array(ilayer)
        do iy_elem = 1,NEX_PER_PROC_ETA,step_mult*ratio_sampling_array(ilayer)

          if (step_mult == 1) then
! for xi direction
            if (.not. CUT_SUPERBRICK_XI) then
              if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                case_xi = 1
              else
                case_xi = 2
              endif
            else
              if (offset_proc_xi == 0) then
                if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                  case_xi = 1
                else
                  case_xi = 2
                endif
              else
                if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))/=0) then
                  case_xi = 1
                else
                  case_xi = 2
                endif
              endif
            endif
! for eta direction
            if (.not. CUT_SUPERBRICK_ETA) then
              if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                case_eta = 1
              else
                case_eta = 2
              endif
            else
              if (offset_proc_eta == 0) then
                if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))==0) then
                  case_eta = 1
                else
                  case_eta = 2
                endif
              else
                if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer)))/=0) then
                  case_eta = 1
                else
                  case_eta = 2
                endif
              endif
            endif
! determine the current sub-block
            if (case_xi == 1) then
              if (case_eta == 1) then
                subblock_num = 1
              else
                subblock_num = 2
              endif
            else
              if (case_eta == 1) then
                subblock_num = 3
              else
                subblock_num = 4
              endif
            endif
! then define the geometry for this sub-block
            call define_basic_doubling_brick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb,subblock_num)
          endif
! loop on all the elements in the mesh doubling superbrick
          do ispec_superbrick = 1,nspec_sb
! loop on all the corner nodes of this element
            do ignod = 1,NGNOD_EIGHT_CORNERS

! define topological coordinates of this mesh point
              offset_x(ignod) = (ix_elem - 1) + &
         x_superbrick(ibool_superbrick(ignod,ispec_superbrick)) * ratio_sampling_array(ilayer)
              offset_y(ignod) = (iy_elem - 1) + &
         y_superbrick(ibool_superbrick(ignod,ispec_superbrick)) * ratio_sampling_array(ilayer)
              offset_z(ignod) = (iz_elem - 1) + &
         z_superbrick(ibool_superbrick(ignod,ispec_superbrick))

            enddo
! the rest of the 27 nodes are missing, therefore add them
     call add_missing_nodes(offset_x,offset_y,offset_z)

! compute the actual position of all the grid points of that element
     call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
               ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
               NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
               r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer,ichunk,rotation_matrix, &
               NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)

! add one spectral element to the list
     ispec = ispec + 1
     if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')

! new get_flag_boundaries

! xmin and xmax
  if (ix_elem == 1) then
      iMPIcut_xi(1,ispec) = iboun_sb(ispec_superbrick,1)
      if (iproc_xi == 0) iboun(1,ispec)= iboun_sb(ispec_superbrick,1)
  endif
  if (ix_elem == (NEX_PER_PROC_XI-step_mult*ratio_sampling_array(ilayer)+1)) then
      iMPIcut_xi(2,ispec) = iboun_sb(ispec_superbrick,2)
      if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= iboun_sb(ispec_superbrick,2)
  endif

! ymin and ymax
  if (iy_elem == 1) then
      iMPIcut_eta(1,ispec) = iboun_sb(ispec_superbrick,3)
      if (iproc_eta == 0) iboun(3,ispec)= iboun_sb(ispec_superbrick,3)
  endif
  if (iy_elem == (NEX_PER_PROC_ETA-step_mult*ratio_sampling_array(ilayer)+1)) then
      iMPIcut_eta(2,ispec) = iboun_sb(ispec_superbrick,4)
      if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= iboun_sb(ispec_superbrick,4)
  endif

! zmax only
  if (ilayer==ifirst_layer) iboun(6,ispec)= iboun_sb(ispec_superbrick,6)
  if (ilayer==ilast_layer .and. iz_elem==1) iboun(5,ispec)= iboun_sb(ispec_superbrick,5)

! define the doubling flag of this element
     if(iregion_code /= IREGION_OUTER_CORE) idoubling(ispec) = doubling_index(ilayer)

! compute several rheological and geometrical properties for this spectral element
     call compute_element_properties(ispec,iregion_code,idoubling, &
           xstore,ystore,zstore,nspec, &
           nspl,rspl,espl,espl2,ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           myrank,ibathy_topo,ATTENUATION,ATTENUATION_3D, &
           ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           RICB,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           xelm,yelm,zelm,shape3D,dershape3D,rmin,rmax,rhostore_local,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
!! DK DK added this for the merged version
           kappavstore_local, &
           xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           nspec_ani,nspec_stacey,Qmu_store,tau_e_store,tau_s,T_c_source,rho_vp,rho_vs,&
           AMM_V,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V,AM_S,AS_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ipass)

!! DK DK added this for merged version
    include "comp_mass_matrix_one_element.f90"
    include "store_xelm_yelm_zelm.f90"

! end of loops on the mesh doubling elements
          enddo
        enddo
      enddo
    endif

! end of loop on all the layers of the mesh
  enddo

  if (.not. PATCH_FOR_GORDON_BELL .and. (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH)) then
    deallocate(stretch_tab,STAT=ier )
    if (ier /= 0) then
      print *,"ERROR can not deallocate stretch_tab in create_regions_mesh ier=",ier
      stop 'error in deallocate'
    endif
  endif

!---

! define central cube in inner core

  if(INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) then

! create the shape of a regular mesh element in the inner core
  call hex_nodes(iaddx,iaddy,iaddz)

! define vertical slice in central cube on current processor
! we can assume that NEX_XI = NEX_ETA, otherwise central cube cannot be defined
  nx_central_cube = NEX_PER_PROC_XI / ratio_divide_central_cube
  ny_central_cube = NEX_PER_PROC_ETA / ratio_divide_central_cube
  nz_central_cube = NEX_XI / ratio_divide_central_cube

! size of the cube along cartesian axes before rotation
  radius_cube = (R_CENTRAL_CUBE / R_EARTH) / sqrt(3.d0)

! define spectral elements in central cube
  do iz = 0,2*nz_central_cube-2,2
    do iy = 0,2*ny_central_cube-2,2
      do ix = 0,2*nx_central_cube-2,2

!       radii that define the shell, we know that we are in the central cube
        rmin = 0.d0
        rmax = R_CENTRAL_CUBE / R_EARTH

!       loop over the NGNOD nodes
        do ia=1,NGNOD

! flat cubed sphere with correct mapping
          call compute_coord_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia), &
                  xgrid_central_cube,ygrid_central_cube,zgrid_central_cube, &
                  iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,nx_central_cube,ny_central_cube,nz_central_cube,radius_cube)

          if(ichunk == CHUNK_AB) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = + xgrid_central_cube
            zelm(ia) = + zgrid_central_cube

          else if(ichunk == CHUNK_AB_ANTIPODE) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = - xgrid_central_cube
            zelm(ia) = - zgrid_central_cube

          else if(ichunk == CHUNK_AC) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = - zgrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else if(ichunk == CHUNK_AC_ANTIPODE) then
            xelm(ia) = - ygrid_central_cube
            yelm(ia) = + zgrid_central_cube
            zelm(ia) = - xgrid_central_cube

          else if(ichunk == CHUNK_BC) then
            xelm(ia) = - zgrid_central_cube
            yelm(ia) = + ygrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else if(ichunk == CHUNK_BC_ANTIPODE) then
            xelm(ia) = + zgrid_central_cube
            yelm(ia) = - ygrid_central_cube
            zelm(ia) = + xgrid_central_cube

          else
            call exit_MPI(myrank,'wrong chunk number in flat cubed sphere definition')
          endif

        enddo

! add one spectral element to the list
        ispec = ispec + 1
        if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in central cube creation')

! new get_flag_boundaries
! xmin & xmax
  if (ix == 0) then
      iMPIcut_xi(1,ispec) = .true.
      if (iproc_xi == 0) iboun(1,ispec)= .true.
  endif
  if (ix == 2*nx_central_cube-2) then
      iMPIcut_xi(2,ispec) = .true.
      if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= .true.
  endif
! ymin & ymax
  if (iy == 0) then
      iMPIcut_eta(1,ispec) = .true.
      if (iproc_eta == 0) iboun(3,ispec)= .true.
  endif
  if (iy == 2*ny_central_cube-2) then
      iMPIcut_eta(2,ispec) = .true.
      if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= .true.
  endif

! define the doubling flag of this element
! only two active central cubes, the four others are fictitious

! determine where we cut the central cube to share it between CHUNK_AB & CHUNK_AB_ANTIPODE
! in the case of mod(NPROC_XI,2)/=0, the cut is asymetric and the bigger part is for CHUNK_AB
  if (mod(NPROC_XI,2)/=0) then
    if (ichunk == CHUNK_AB) then
      nz_inf_limit = ((nz_central_cube*2)/NPROC_XI)*floor(NPROC_XI/2.d0)
    elseif (ichunk == CHUNK_AB_ANTIPODE) then
      nz_inf_limit = ((nz_central_cube*2)/NPROC_XI)*ceiling(NPROC_XI/2.d0)
    endif
  else
    nz_inf_limit = nz_central_cube
  endif

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    if(iz == nz_inf_limit) then
      idoubling(ispec) = IFLAG_BOTTOM_CENTRAL_CUBE
    else if(iz == 2*nz_central_cube-2) then
      idoubling(ispec) = IFLAG_TOP_CENTRAL_CUBE
    else if (iz > nz_inf_limit .and. iz < 2*nz_central_cube-2) then
      idoubling(ispec) = IFLAG_MIDDLE_CENTRAL_CUBE
    else
      idoubling(ispec) = IFLAG_IN_FICTITIOUS_CUBE
    endif
  else
    idoubling(ispec) = IFLAG_IN_FICTITIOUS_CUBE
  endif

! compute several rheological and geometrical properties for this spectral element
     call compute_element_properties(ispec,iregion_code,idoubling, &
           xstore,ystore,zstore,nspec, &
           nspl,rspl,espl,espl2,ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           myrank,ibathy_topo,ATTENUATION,ATTENUATION_3D, &
           ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
           RICB,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
           xelm,yelm,zelm,shape3D,dershape3D,rmin,rmax,rhostore_local,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
!! DK DK added this for the merged version
           kappavstore_local, &
           xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           nspec_ani,nspec_stacey,Qmu_store,tau_e_store,tau_s,T_c_source,rho_vp,rho_vs,&
           AMM_V,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V,AM_S,AS_V, &
           numker,numhpa,numcof,ihpa,lmax,nylm, &
           lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
           nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
           coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ipass)

!! DK DK added this for merged version
    include "comp_mass_matrix_one_element.f90"
    include "store_xelm_yelm_zelm.f90"

      enddo
    enddo
  enddo

  endif    ! end of definition of central cube in inner core

!---

! check total number of spectral elements created
  if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

! only create global addressing and the MPI buffers in the first pass
  if(ipass == 1) then

    locval = 0
    ifseg = .false.

! arrays xstore,ystore,zstore are destroyed by this sorting routine but
! they will be recreated in the second pass
    call get_global(nspec,xstore,ystore,zstore,ibool,locval,ifseg,nglob,npointot)

  ! check that number of points found equals theoretical value
    if(nglob /= nglob_theor) then
      write(errmsg,*) 'incorrect total number of points found: myrank,nglob,nglob_theor,ipass,iregion_code = ',&
        myrank,nglob,nglob_theor,ipass,iregion_code
      call exit_MPI(myrank,errmsg)
    endif

    if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering')

! create a new indirect addressing to reduce cache misses in memory access in the solver
! this is *critical* to improve performance in the solver
    mask_ibool(:) = -1
    copy_ibool_ori(:,:,:,:) = ibool(:,:,:,:)

    inumber = 0
    do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(DEBUG .and. copy_ibool_ori(i,j,k,ispec) < 1) then
              print *,'error: copy_ibool is < 1 in debug1 of second sorting: ',copy_ibool_ori(i,j,k,ispec)
              call exit_MPI(myrank,'error: copy_ibool is < 1 in debug1 of second sorting')
            endif
            if(DEBUG .and. copy_ibool_ori(i,j,k,ispec) > nglob_theor) then
              print *,'error: copy_ibool is too high in debug1 of second sorting: ',copy_ibool_ori(i,j,k,ispec),nglob_theor
              call exit_MPI(myrank,'error: copy_ibool is too high in debug1 of second sorting')
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          if(mask_ibool(copy_ibool_ori(i,j,k,ispec)) == -1) then
  ! create a new point
            inumber = inumber + 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(DEBUG .and. inumber < 1) then
              print *,'error: inumber is < 1 in debug2 of second sorting: ',inumber
              call exit_MPI(myrank,'error: inumber is < 1 in debug2 of second sorting')
            endif
            if(DEBUG .and. inumber > nglob_theor) then
              print *,'error: inumber is too high in debug2 of second sorting: ',inumber,nglob_theor
              call exit_MPI(myrank,'error: inumber is too high in debug2 of second sorting')
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ibool(i,j,k,ispec) = inumber
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(DEBUG .and. copy_ibool_ori(i,j,k,ispec) < 1) then
              print *,'error: copy_ibool is < 1 in debug3 of second sorting: ',copy_ibool_ori(i,j,k,ispec)
              call exit_MPI(myrank,'error: copy_ibool is < 1 in debug3 of second sorting')
            endif
            if(DEBUG .and. copy_ibool_ori(i,j,k,ispec) > nglob_theor) then
              print *,'error: copy_ibool is too high in debug3 of second sorting: ',copy_ibool_ori(i,j,k,ispec),nglob_theor
              call exit_MPI(myrank,'error: copy_ibool is too high in debug3 of second sorting')
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            mask_ibool(copy_ibool_ori(i,j,k,ispec)) = inumber
          else
  ! use an existing point created previously
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(DEBUG .and. copy_ibool_ori(i,j,k,ispec) < 1) then
              print *,'error: copy_ibool is < 1 in debug4 of second sorting: ',copy_ibool_ori(i,j,k,ispec)
              call exit_MPI(myrank,'error: copy_ibool is < 1 in debug4 of second sorting')
            endif
            if(DEBUG .and. copy_ibool_ori(i,j,k,ispec) > nglob_theor) then
              print *,'error: copy_ibool is too high in debug4 of second sorting: ',copy_ibool_ori(i,j,k,ispec),nglob_theor
              call exit_MPI(myrank,'error: copy_ibool is too high in debug4 of second sorting')
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ibool(i,j,k,ispec) = mask_ibool(copy_ibool_ori(i,j,k,ispec))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(DEBUG .and. ibool(i,j,k,ispec) < 1) then
              print *,'error: ibool is < 1 in debug5 of second sorting: ',ibool(i,j,k,ispec)
              call exit_MPI(myrank,'error: ibool is < 1 in debug5 of second sorting')
            endif
            if(DEBUG .and. ibool(i,j,k,ispec) > nglob_theor) then
              print *,'error: ibool is too high in debug5 of second sorting: ',ibool(i,j,k,ispec) ,nglob_theor
              call exit_MPI(myrank,'error: ibool is too high in debug5 of second sorting')
            endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          endif
        enddo
      enddo
    enddo
    enddo

    if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) then
      print *,'Incorrect global numbering after second sorting: mask_ibool ',minval(mask_ibool),maxval(mask_ibool)
      print *,'Incorrect global numbering after second sorting: ibool ',minval(ibool),maxval(ibool)
      print *,'Incorrect global numbering after second sorting: expected min/max = 1, ',nglob_theor
      call exit_MPI(myrank,'incorrect global numbering after second sorting')
    endif

! create MPI buffers and mass matrix in the second pass
  else if(ipass == 2) then

! copy the theoretical number of points for the second pass
  nglob = nglob_theor

! create MPI buffers
! arrays locval(npointot) and ifseg(npointot) used to save memory
#ifdef USE_MPI
    call get_MPI_cutplanes_xi(myrank,nspec,iMPIcut_xi,ibool, &
                  ifseg,npointot, &
                  NSPEC2D_ETA_FACE,iregion_code,nglob,iboolleft_xi,iboolright_xi,NGLOB2DMAX_XMIN_XMAX,npoin2D_xi)

    call get_MPI_cutplanes_eta(myrank,nspec,iMPIcut_eta,ibool, &
                  ifseg,npointot, &
                  NSPEC2D_XI_FACE,iregion_code,nglob,iboolleft_eta,iboolright_eta,NGLOB2DMAX_YMIN_YMAX,npoin2D_eta)

    call get_MPI_1D_buffers(myrank,nspec,iMPIcut_xi,iMPIcut_eta,ibool,idoubling, &
                  xstore,ystore,zstore,ifseg,npointot, &
                  NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,iregion_code,nglob, &
                  ibool1D_leftxi_lefteta,ibool1D_rightxi_lefteta, &
                  ibool1D_leftxi_righteta,ibool1D_rightxi_righteta,NGLOB1D_RADIAL_MAX, &
  xread1D_leftxi_lefteta, xread1D_rightxi_lefteta, xread1D_leftxi_righteta, xread1D_rightxi_righteta, &
  yread1D_leftxi_lefteta, yread1D_rightxi_lefteta, yread1D_leftxi_righteta, yread1D_rightxi_righteta, &
  zread1D_leftxi_lefteta, zread1D_rightxi_lefteta, zread1D_leftxi_righteta, zread1D_rightxi_righteta, &
  iregion_code)
#endif

! count number of anisotropic elements in current region
! should be zero in all the regions except in the mantle
  if(iregion_code == IREGION_CRUST_MANTLE)  then
    nspec_tiso = count(idoubling(1:nspec) == IFLAG_220_80) + count(idoubling(1:nspec) == IFLAG_80_MOHO)
  else
    nspec_tiso = 0
  endif

  call get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
      dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
      ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
      nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
              jacobian2D_xmin,jacobian2D_xmax, &
              jacobian2D_ymin,jacobian2D_ymax, &
              jacobian2D_bottom,jacobian2D_top, &
              normal_xmin,normal_xmax, &
              normal_ymin,normal_ymax, &
              normal_bottom,normal_top, &
              NSPEC2D_BOTTOM,NSPEC2D_TOP, &
              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

! save the binary files
! save ocean load mass matrix as well if oceans
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then

! adding ocean load mass matrix at the top of the crust for oceans
  nglob_oceans = nglob

! create ocean load mass matrix for degrees of freedom at ocean bottom
  rmass_ocean_load(:) = 0._CUSTOM_REAL

! add contribution of the oceans
! for surface elements exactly at the top of the crust (ocean bottom)
    do ispec2D_top_crust = 1,NSPEC2D_TOP

      ispec_oceans = ibelm_top(ispec2D_top_crust)

      iz_oceans = NGLLZ

      do ix_oceans = 1,NGLLX
        do iy_oceans = 1,NGLLY

        iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! if 3D Earth, compute local height of oceans
        if(CASE_3D) then

! get coordinates of current point
          xval = xstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          yval = ystore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          zval = zstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! map to latitude and longitude for bathymetry routine
          call xyz_2_rthetaphi_dble(xval,yval,zval,rval,thetaval,phival)
          call reduce(thetaval,phival)

! convert the geocentric colatitude to a geographic colatitude
          colat = PI/2.0d0 - datan(1.006760466d0*dcos(thetaval)/dmax1(TINYVAL,dsin(thetaval)))

! get geographic latitude and longitude in degrees
          lat = 90.0d0 - colat*180.0d0/PI
          lon = phival*180.0d0/PI
          elevation = 0.d0

! compute elevation at current point
          call get_topo_bathy(lat,lon,elevation,ibathy_topo)

! non-dimensionalize the elevation, which is in meters
! and suppress positive elevation, which means no oceans
          if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation) / R_EARTH
          endif

        else
! if 1D Earth, use oceans of constant thickness everywhere
          height_oceans = THICKNESS_OCEANS_PREM
        endif

! take into account inertia of water column
        weight = wxgll(ix_oceans)*wygll(iy_oceans)*dble(jacobian2D_top(ix_oceans,iy_oceans,ispec2D_top_crust)) &
                   * dble(RHO_OCEANS) * height_oceans

! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + sngl(weight)
        else
          rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + weight
        endif

        enddo
      enddo

    enddo

! add regular mass matrix to ocean load contribution
  rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

! allocate dummy array if no oceans
    nglob_oceans = 1

  endif

!! DK DK save Brian's attenuation files to a shared disk
!! DK DK obviously we should do this with MPI or with subroutine arguments
!! DK DK shared by the mesher and the solver subroutines at some point
  if(ATTENUATION_VAL) call attenuation_save_arrays(iregion_code, AM_V)

! create AVS or DX mesh data for the slices
  if(SAVE_MESH_FILES) then
    call create_name_database(prname,myrank,iregion_code)

    call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot,iregion_code)

    call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot,iregion_code)

    call write_AVS_DX_global_chunks_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot, &
!! DK DK changed for now because array rhostore is not available in v4.1 anymore
!! DK DK      rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
              nspl,rspl,espl,espl2, &
              ELLIPTICITY,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,REFERENCE_1D_MODEL, &
              RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
              RMIDDLE_CRUST,ROCEAN,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,iregion_code)

    call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot,iregion_code)
  endif

! compute volume, bottom and top area of that part of the slice
  volume_local = ZERO
  area_local_bottom = ZERO
  area_local_top = ZERO

  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

! compute the jacobian
!! DK DK in merged version the jacobian is not stored anymore and therefore not valid anymore
  goto 777
          xixl = xixstore(i,j,k)
          xiyl = xiystore(i,j,k)
          xizl = xizstore(i,j,k)
          etaxl = etaxstore(i,j,k)
          etayl = etaystore(i,j,k)
          etazl = etazstore(i,j,k)
          gammaxl = gammaxstore(i,j,k)
          gammayl = gammaystore(i,j,k)
          gammazl = gammazstore(i,j,k)

          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          volume_local = volume_local + dble(jacobianl)*weight
!! DK DK in merged version the jacobian is not stored anymore and therefore not valid anymore
  777 continue

        enddo
      enddo
    enddo
  enddo

  do ispec = 1,NSPEC2D_BOTTOM
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_bottom = area_local_bottom + dble(jacobian2D_bottom(i,j,ispec))*weight
      enddo
    enddo
  enddo

  do ispec = 1,NSPEC2D_TOP
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_top = area_local_top + dble(jacobian2D_top(i,j,ispec))*weight
      enddo
    enddo
  enddo

  else
    stop 'there cannot be more than two passes in mesh creation'
  endif  ! end of test if first or second pass

  enddo ! of loop on ipass = 1,2

  end subroutine create_regions_mesh

