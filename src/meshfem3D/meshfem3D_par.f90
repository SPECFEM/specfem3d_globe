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
