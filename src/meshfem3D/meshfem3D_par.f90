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

  end module constants

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

  implicit none


! model_aniso_mantle_variables
  type model_aniso_mantle_variables
    sequence
    double precision beta(14,34,37,73)
    double precision pro(47)
    integer npar1
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_aniso_mantle_variables
  type (model_aniso_mantle_variables) AMM_V
! model_aniso_mantle_variables

! model_attenuation_variables
  type model_attenuation_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    integer                                   :: Qn                 ! Number of points
    integer dummy_pad ! padding 4 bytes to align the structure
  end type model_attenuation_variables
  type (model_attenuation_variables) AM_V
! model_attenuation_variables

! model_atten3D_QRFSI12_variables
  type model_atten3D_QRFSI12_variables
    sequence
    double precision dqmu(NKQ,NSQ)
    double precision spknt(NKQ)
    double precision refdepth(NDEPTHS_REFQ)
    double precision refqmu(NDEPTHS_REFQ)
  end type model_atten3D_QRFSI12_variables
  type (model_atten3D_QRFSI12_variables) QRFSI12_Q
! model_atten3D_QRFSI12_variables

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

! model_1dref_variables
  type model_1dref_variables
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
  end type model_1dref_variables
 type (model_1dref_variables) Mref_V
! model_1dref_variables

! model_sea1d_variables
  type model_sea1d_variables
    sequence
     double precision, dimension(NR_SEA1D) :: radius_sea1d
     double precision, dimension(NR_SEA1D) :: density_sea1d
     double precision, dimension(NR_SEA1D) :: vp_sea1d
     double precision, dimension(NR_SEA1D) :: vs_sea1d
     double precision, dimension(NR_SEA1D) :: Qkappa_sea1d
     double precision, dimension(NR_SEA1D) :: Qmu_sea1d
  end type model_sea1d_variables
  type (model_sea1d_variables) SEA1DM_V
! model_sea1d_variables

! model_s20rts_variables
  type model_s20rts_variables
    sequence
    double precision dvs_a(0:NK_20,0:NS_20,0:NS_20)   !a = positive m  (radial, theta, phi) --> (k,l,m) (maybe other way around??)
    double precision dvs_b(0:NK_20,0:NS_20,0:NS_20)   !b = negative m  (radial, theta, phi) --> (k,l,-m)
    double precision dvp_a(0:NK_20,0:NS_20,0:NS_20)
    double precision dvp_b(0:NK_20,0:NS_20,0:NS_20)
    double precision spknt(NK_20+1)
    double precision qq0(NK_20+1,NK_20+1)
    double precision qq(3,NK_20+1,NK_20+1)
  end type model_s20rts_variables
  type (model_s20rts_variables) S20RTS_V
! model_s20rts_variables

! model_s40rts_variables
  type model_s40rts_variables
    sequence
    double precision dvs_a(0:NK_20,0:NS_40,0:NS_40)
    double precision dvs_b(0:NK_20,0:NS_40,0:NS_40)
    double precision dvp_a(0:NK_20,0:NS_40,0:NS_40)
    double precision dvp_b(0:NK_20,0:NS_40,0:NS_40)
    double precision spknt(NK_20+1)
    double precision qq0(NK_20+1,NK_20+1)
    double precision qq(3,NK_20+1,NK_20+1)
  end type model_s40rts_variables
  type (model_s40rts_variables) S40RTS_V
! model_s40rts_variables

! model_heterogen_m_variables
  type model_heterogen_m_variables
    sequence
    double precision rho_in(N_R*N_THETA*N_PHI)
  end type model_heterogen_m_variables
  type (model_heterogen_m_variables) HMM
! model_heterogen_m_variables

! model_jp3d_variables
  type model_jp3d_variables
    sequence
    ! vmod3d
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
    double precision :: PLA
    double precision :: RLA
    double precision :: HLA
    double precision :: PLB
    double precision :: RLB
    double precision :: HLB
    ! weight
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
    ! locate integers
    integer :: IPLOCA(MKA)
    integer :: IRLOCA(MKA)
    integer :: IHLOCA(MKA)
    integer :: IPLOCB(MKB)
    integer :: IRLOCB(MKB)
    integer :: IHLOCB(MKB)
    ! vmod3D integers
    integer :: NPA
    integer :: NRA
    integer :: NHA
    integer :: NPB
    integer :: NRB
    integer :: NHB
    ! weight integers
    integer :: IP
    integer :: JP
    integer :: KP
    integer :: IP1
    integer :: JP1
    integer :: KP1
  end type model_jp3d_variables
  type (model_jp3d_variables) JP3DM_V
! model_jp3d_variables

! model_sea99_s_variables
  type model_sea99_s_variables
    sequence
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    integer :: dummy_pad ! padding 4 bytes to align the structure
 end type model_sea99_s_variables
 type (model_sea99_s_variables) SEA99M_V
! model_sea99_s_variables

! crust 2.0 model_crust_variables
  type model_crust_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
    character(len=2) dummy_pad ! padding 2 bytes to align the structure
  end type model_crust_variables
  type (model_crust_variables) CM_V
! model_crust_variables

! EUcrust
  type model_eucrust_variables
    sequence
    double precision, dimension(:),pointer :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
    integer :: num_eucrust
    integer :: dummy_pad ! padding 4 bytes to align the structure
  end type model_eucrust_variables
  type (model_eucrust_variables) EUCM_V

! type for EPCRUST 1.0
  type model_epcrust_variables
    sequence
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT):: lon_ep,lat_ep,topo_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: thickness_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vp_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vs_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: rho_ep
  end type model_epcrust_variables
  type (model_epcrust_variables) EPCRUST

! model_crustmaps_variables combined crustal maps
  type model_crustmaps_variables
    sequence
    double precision, dimension(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: thickness
    double precision, dimension(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: density
    double precision, dimension(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocp
    double precision, dimension(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP) :: velocs
    double precision thicknessnp(NLAYERS_CRUSTMAP)
    double precision densitynp(NLAYERS_CRUSTMAP)
    double precision velocpnp(NLAYERS_CRUSTMAP)
    double precision velocsnp(NLAYERS_CRUSTMAP)
    double precision thicknesssp(NLAYERS_CRUSTMAP)
    double precision densitysp(NLAYERS_CRUSTMAP)
    double precision velocpsp(NLAYERS_CRUSTMAP)
    double precision velocssp(NLAYERS_CRUSTMAP)
  end type model_crustmaps_variables
  type (model_crustmaps_variables) GC_V
!model_crustmaps_variables

! model_attenuation_storage_var
  type model_attenuation_storage_var
    sequence
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
    integer Q_resolution
    integer Q_max
  end type model_attenuation_storage_var
  type (model_attenuation_storage_var) AM_S
! model_attenuation_storage_var

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
  end type attenuation_simplex_variables
  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables

! point profile model_variables
  type model_ppm_variables
    sequence
    double precision,dimension(:),pointer :: dvs,lat,lon,depth
    double precision :: maxlat,maxlon,minlat,minlon,maxdepth,mindepth
    double precision :: dlat,dlon,ddepth,max_dvs,min_dvs
    integer :: num_v,num_latperlon,num_lonperdepth
    integer :: dummy_pad ! padding 4 bytes to align the structure
  end type model_ppm_variables
  type (model_ppm_variables) PPM_V

! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    logical :: MODEL_GLL
    logical,dimension(3) :: dummy_pad ! padding 3 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) MGLL_V

! bathymetry and topography: use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

! used for 3D Harvard models s362ani, s362wmani, s362ani_prem and s2.9ea
  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

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

  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)
  integer itpspl(maxcoe,maxhpa)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm

  character(len=80) kerstr
  character(len=80) refmdl
  character(len=40) varstr(maxker)
  character(len=80) hsplfl(maxhpa)
  character(len=40) dskker(maxker)


! for ellipticity
  double precision rspl(NR),espl(NR),espl2(NR)
  integer nspl

! model parameter and flags
  integer REFERENCE_1D_MODEL,THREE_D_MODEL

  logical ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS

  logical HONOR_1D_SPHERICAL_MOHO,CRUSTAL,ONE_CRUST,CASE_3D,TRANSVERSE_ISOTROPY

  logical ISOTROPIC_3D_MANTLE,ANISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE

  logical ATTENUATION,ATTENUATION_NEW,ATTENUATION_3D

  logical ANISOTROPIC_INNER_CORE

! to create a reference model based on 1D_REF but with 3D crust and 410/660 topography
  logical,parameter :: USE_1D_REFERENCE = .false.

  end module meshfem3D_models_par


!
!-------------------------------------------------------------------------------------------------
!


  module meshfem3D_par

! main parameter module for specfem simulations

  use constants

  implicit none

  ! correct number of spectral elements in each block depending on chunk type
  integer npointot

  ! proc numbers for MPI
  integer myrank,sizeprocs

  ! check area and volume of the final mesh
  double precision volume_total

  ! for loop on all the slices
  integer iregion_code
  integer iproc_xi,iproc_eta,ichunk

  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! for some statistics for the mesh
  integer numelem_crust_mantle,numelem_outer_core,numelem_inner_core
  integer numelem_total

  ! timer MPI
  double precision time_start,tCPU

  ! addressing for all the slices
  integer, dimension(:), allocatable :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer, dimension(:,:,:), allocatable :: addressing

  ! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN, &
          NCHUNKS,SIMULATION_TYPE, &
          MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP,NOISE_TOMOGRAPHY

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
          RMOHO_FICTITIOUS_IN_MESHER

  logical MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_COARSE, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

  character(len=150) :: OUTPUT_FILES,LOCAL_PATH,LOCAL_TMP_PATH,MODEL

  ! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

  ! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP,NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               NGLOB

  ! computed in read_compute_parameters
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  ! memory size of all the static arrays
  double precision :: static_memory_size

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
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

  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL_TEMP

  end module meshfem3D_par

!
!-------------------------------------------------------------------------------------------------
!

  module create_regions_mesh_par

  use constants,only: NGLLX,NGLLY,NGLLZ,NGNOD,NGNOD2D,NDIM,NDIM2D

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

  end module create_regions_mesh_par

!
!-------------------------------------------------------------------------------------------------
!

  module create_regions_mesh_par2

  use constants,only: CUSTOM_REAL,N_SLS

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

  ! mass matrix and bathymetry for ocean load
  integer nglob_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

  ! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax, &
    ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

  ! 2-D jacobians and normals
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


  ! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  ! attenuation
  double precision, dimension(:,:,:,:),   allocatable :: Qmu_store
  double precision, dimension(:,:,:,:,:), allocatable :: tau_e_store
  double precision, dimension(N_SLS)                  :: tau_s
  double precision  T_c_source

  logical :: USE_ONE_LAYER_SB

  integer NUMBER_OF_MESH_LAYERS,layer_shift,cpt, &
    first_layer_aniso,last_layer_aniso,FIRST_ELT_NON_ANISO

  double precision, dimension(:,:), allocatable :: stretch_tab

  integer :: nb_layer_above_aniso,FIRST_ELT_ABOVE_ANISO

  logical :: ACTUALLY_STORE_ARRAYS

  ! Boundary Mesh
  integer NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho
  integer, dimension(:), allocatable :: ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_moho,normal_400,normal_670
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_moho,jacobian2D_400,jacobian2D_670
  integer ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot
  double precision r_moho,r_400,r_670

  ! flags for transverse isotropic elements
  logical, dimension(:), allocatable :: ispec_is_tiso

  ! name of the database file
  character(len=150) :: prname

  end module create_regions_mesh_par2

!
!-------------------------------------------------------------------------------------------------
!

  module create_MPI_interfaces_par

  use constants,only: &
    CUSTOM_REAL,NDIM,IMAIN, &
    IREGION_CRUST_MANTLE,IREGION_OUTER_CORE,IREGION_INNER_CORE

  implicit none
  
  ! indirect addressing for each message for faces and corners of the chunks
  ! a given slice can belong to at most one corner and at most two faces
  integer :: NGLOB2DMAX_XY

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
  
  end module create_MPI_interfaces_par

!
!-------------------------------------------------------------------------------------------------
!

  module MPI_crust_mantle_par

  use constants,only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR
  
  implicit none
  
  ! collected MPI interfaces
  !--------------------------------------
  ! MPI crust/mantle mesh
  !--------------------------------------
  integer :: num_interfaces_crust_mantle
  integer :: max_nibool_interfaces_crust_mantle
  integer, dimension(:), allocatable :: my_neighbours_crust_mantle,nibool_interfaces_crust_mantle
  integer, dimension(:,:), allocatable :: ibool_interfaces_crust_mantle

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_crust_mantle,buffer_recv_vector_crust_mantle

  integer, dimension(:), allocatable :: request_send_vector_crust_mantle,request_recv_vector_crust_mantle

  ! temporary arrays for elements on slices or edges
  logical, dimension(:),allocatable :: is_on_a_slice_edge_crust_mantle

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

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  integer, dimension(:),allocatable :: idoubling_crust_mantle
  integer, dimension(:,:,:,:),allocatable :: ibool_crust_mantle

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

  use constants,only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR
  
  implicit none
  
  !--------------------------------------
  ! MPI inner core mesh
  !--------------------------------------
  integer :: num_interfaces_inner_core
  integer :: max_nibool_interfaces_inner_core
  integer, dimension(:), allocatable :: my_neighbours_inner_core,nibool_interfaces_inner_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_inner_core

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: buffer_send_vector_inner_core,buffer_recv_vector_inner_core

  integer, dimension(:), allocatable :: request_send_vector_inner_core,request_recv_vector_inner_core

  ! temporary arrays for elements on slices or edges
  logical, dimension(:),allocatable :: is_on_a_slice_edge_inner_core

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

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    xstore_inner_core,ystore_inner_core,zstore_inner_core
  integer, dimension(:),allocatable :: idoubling_inner_core
  integer, dimension(:,:,:,:),allocatable :: ibool_inner_core


  ! for matching with central cube in inner core
  integer, dimension(:), allocatable :: sender_from_slices_to_cube
  integer, dimension(:,:), allocatable :: ibool_central_cube
  double precision, dimension(:,:), allocatable :: buffer_slices,buffer_slices2
  double precision, dimension(:,:,:), allocatable :: buffer_all_cube_from_slices
  integer nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,receiver_cube_from_slices

  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
            nspec2D_ymin_inner_core,nspec2D_ymax_inner_core

  integer, dimension(:),allocatable :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(:),allocatable :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(:),allocatable :: ibelm_bottom_inner_core
  integer, dimension(:),allocatable :: ibelm_top_inner_core

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

  use constants,only: CUSTOM_REAL,NUMFACES_SHARED,NB_SQUARE_EDGES_ONEDIR
  
  implicit none
  
  !--------------------------------------
  ! MPI outer core mesh
  !--------------------------------------
  integer :: num_interfaces_outer_core
  integer :: max_nibool_interfaces_outer_core
  integer, dimension(:), allocatable :: my_neighbours_outer_core,nibool_interfaces_outer_core
  integer, dimension(:,:), allocatable :: ibool_interfaces_outer_core

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core

  integer, dimension(:), allocatable :: request_send_scalar_outer_core,request_recv_scalar_outer_core

  ! temporary arrays for elements on slices or edges
  logical, dimension(:),allocatable :: is_on_a_slice_edge_outer_core

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

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    xstore_outer_core,ystore_outer_core,zstore_outer_core
  integer, dimension(:),allocatable :: idoubling_outer_core
  integer, dimension(:,:,:,:),allocatable :: ibool_outer_core

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
