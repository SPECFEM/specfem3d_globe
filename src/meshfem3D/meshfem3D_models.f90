!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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

  module meshfem3D_models_par

!---
!
! ADD YOUR MODEL HERE
!
!---

  implicit none

  include "constants.h"

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
    double precision, dimension(NR_AK135F_NO_MUD) :: radius_ak135
    double precision, dimension(NR_AK135F_NO_MUD) :: density_ak135
    double precision, dimension(NR_AK135F_NO_MUD) :: vp_ak135
    double precision, dimension(NR_AK135F_NO_MUD) :: vs_ak135
    double precision, dimension(NR_AK135F_NO_MUD) :: Qkappa_ak135
    double precision, dimension(NR_AK135F_NO_MUD) :: Qmu_ak135
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

  logical ATTENUATION,ATTENUATION_3D

  logical ANISOTROPIC_INNER_CORE

! to create a reference model based on 1D_REF but with 3D crust and 410/660 topography
  logical,parameter :: USE_1D_REFERENCE = .false.

  end module meshfem3D_models_par


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_broadcast(myrank,NSPEC, &
                        MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,&
                        R80,R220,R670,RCMB,RICB)

! preparing model parameter coefficients on all processes

  use meshfem3D_models_par
  use mpi

  implicit none

  integer myrank
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC

  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD

  double precision R80,R220,R670,RCMB,RICB

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! sets up spline coefficients for ellipticity
  if(ELLIPTICITY) &
    call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

  ! GLL model uses s29ea as reference 3D model
  if( THREE_D_MODEL == THREE_D_MODEL_GLL ) then
    MGLL_V%MODEL_GLL = .true.
    ! sets to initial reference model from which iterations started
    THREE_D_MODEL = GLL_REFERENCE_MODEL
  else
    MGLL_V%MODEL_GLL = .false.
  endif

  ! reads in 3D mantle models
  if(ISOTROPIC_3D_MANTLE) then

    select case( THREE_D_MODEL )

      case(THREE_D_MODEL_S20RTS)
        call model_s20rts_broadcast(myrank,S20RTS_V)

      case(THREE_D_MODEL_S40RTS)
        call model_s40rts_broadcast(myrank,S40RTS_V)

      case(THREE_D_MODEL_SEA99_JP3D)
        ! the variables read are declared and stored in structure SEA99M_V and JP3DM_V
        call model_sea99_s_broadcast(myrank,SEA99M_V)
        call model_jp3d_broadcast(myrank,JP3DM_V)

      case(THREE_D_MODEL_SEA99)
        ! the variables read are declared and stored in structure SEA99M_V
        call model_sea99_s_broadcast(myrank,SEA99M_V)

      case(THREE_D_MODEL_JP3D)
        ! the variables read are declared and stored in structure JP3DM_V
        call model_jp3d_broadcast(myrank,JP3DM_V)

      case(THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
           THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
        call model_s362ani_broadcast(myrank,THREE_D_MODEL,numker,numhpa,ihpa,&
                                lmxhpa,itypehpa,ihpakern,numcoe,ivarkern,itpspl, &
                                xlaspl,xlospl,radspl,coe,hsplfl,dskker,kerstr,varstr,refmdl)

      case(THREE_D_MODEL_PPM)
        ! Point Profile Models
        ! the variables read are declared and stored in structure PPM_V
        call model_ppm_broadcast(myrank,PPM_V)

        ! could use EUcrust07 Vp crustal structure
        !call model_eucrust_broadcast(myrank,EUCM_V)

      case(THREE_D_MODEL_GAPP2)
        ! GAP model
        call model_gapp2_broadcast(myrank)

      case default
        call exit_MPI(myrank,'3D model not defined')

    end select

  endif

  ! arbitrary mantle models
  if(HETEROGEN_3D_MANTLE) &
    call model_heterogen_mntl_broadcast(myrank,HMM)

  ! anisotropic mantle
  if(ANISOTROPIC_3D_MANTLE) &
    call model_aniso_mantle_broadcast(myrank,AMM_V)

  ! crustal model
  if(CRUSTAL) &
    call meshfem3D_crust_broadcast(myrank)

  ! GLL model
  if( MGLL_V%MODEL_GLL ) &
    call model_gll_broadcast(myrank,MGLL_V,NSPEC)

  ! attenuation
  if(ATTENUATION ) then
    call model_attenuation_broadcast(myrank,AM_V,MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

    ! 3D attenuation
    if( ATTENUATION_3D) then
      ! Colleen's model defined originally between 24.4km and 650km
      call model_atten3D_QRFSI12_broadcast(myrank,QRFSI12_Q)
    else
      ! sets up attenuation coefficients according to the chosen, "pure" 1D model
      ! (including their 1D-crustal profiles)
      call model_attenuation_setup(REFERENCE_1D_MODEL, RICB, RCMB, &
              R670, R220, R80,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,AM_S,AS_V)
    endif

  endif

  ! read topography and bathymetry file
  if(TOPOGRAPHY .or. OCEANS) &
    call model_topo_bathy_broadcast(myrank,ibathy_topo)

  ! re-defines/initializes models 1066a and ak135 and ref
  ! ( with possible external crustal model: if CRUSTAL is set to true
  !    it strips the 1-D crustal profile and replaces it with mantle properties)
  select case( REFERENCE_1D_MODEL )

    case(REFERENCE_MODEL_1066A)
      call model_1066a_broadcast(CRUSTAL,M1066a_V)

    case( REFERENCE_MODEL_AK135F_NO_MUD)
      call model_ak135_broadcast(CRUSTAL,Mak135_V)

    case(REFERENCE_MODEL_1DREF)
      call model_1dref_broadcast(CRUSTAL,Mref_V)

    case(REFERENCE_MODEL_SEA1D)
      call model_sea1d_broadcast(CRUSTAL,SEA1DM_V)

  end select

  end subroutine meshfem3D_models_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_crust_broadcast(myrank)

! preparing model parameter coefficients on all processes

  use meshfem3D_models_par
  use mpi

  implicit none

  integer myrank

!---
!
! ADD YOUR MODEL HERE
!
!---

  select case (ITYPE_CRUSTAL_MODEL )

    case (ICRUST_CRUST2)
      ! crust 2.0
      call model_crust_broadcast(myrank,CM_V)

    case (ICRUST_CRUSTMAPS)
      ! general crustmaps
      call model_crustmaps_broadcast(myrank,GC_V)

    case (ICRUST_EPCRUST)
      ! EPcrust
      call model_epcrust_broadcast(myrank,EPCRUST)

    case default
      stop 'crustal model type not defined'

  end select


  end subroutine meshfem3D_crust_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_get1D_val(myrank,iregion_code,idoubling, &
                              r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                              Qkappa,Qmu,RICB,RCMB, &
                              RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                              RMOHO,RMIDDLE_CRUST,ROCEAN)
! reference model values
!
! for a given location radius (r_prem, which is the point's radius with tolerance factor),
! this calculates density and velocities
!
! note: if CRUSTAL is set, it strips the 1-D crustal profile and mantle gets expanded
!          up to the surface.
!          only exception is JP1D...
!
! routine returns: rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu

  use meshfem3D_models_par

  implicit none

  integer myrank,iregion_code,idoubling
  double precision r_prem,rho
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision Qkappa,Qmu
  double precision RICB,RCMB,RTOPDDOUBLEPRIME,R80,R120,R220,R400, &
    R600,R670,R771,RMOHO,RMIDDLE_CRUST,ROCEAN

  ! local parameters
  double precision drhodr,vp,vs

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! gets 1-D reference model parameters
  select case ( REFERENCE_1D_MODEL )

    case(REFERENCE_MODEL_PREM)
      ! PREM (by Dziewonski & Anderson) - used also as background for 3D models
      if(TRANSVERSE_ISOTROPY) then
        ! get the anisotropic PREM parameters
        call model_prem_aniso(myrank,r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                  Qkappa,Qmu,idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
                  R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
      else
        ! isotropic model
        call model_prem_iso(myrank,r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
                  ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                  R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
      endif

    case(REFERENCE_MODEL_1DREF)
      ! 1D-REF also known as STW105 (by Kustowski et al.) - used also as background for 3D models
      call model_1dref(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL,Mref_V)
      if(.not. TRANSVERSE_ISOTROPY) then
        if(.not. ISOTROPIC_3D_MANTLE) then
          ! this case here is only executed for 1D_ref_iso
          ! calculates isotropic values
          vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                    + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
          vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                    + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
        endif
      endif

    case(REFERENCE_MODEL_1066A)
      ! 1066A (by Gilbert & Dziewonski) - pure isotropic model, used in 1D model mode only
      call model_1066a(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,M1066a_V)

    case(REFERENCE_MODEL_AK135F_NO_MUD)
      ! AK135 (by Kennett et al. ) - pure isotropic model, used in 1D model mode only
      call model_ak135(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,Mak135_V)

    case(REFERENCE_MODEL_IASP91)
      ! IASP91 (by Kennett & Engdahl) - pure isotropic model, used in 1D model mode only
      call model_iasp91(myrank,r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
                    ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                    R771,R670,R400,R220,R120,RMOHO,RMIDDLE_CRUST)

    case(REFERENCE_MODEL_JP1D)
      !JP1D (by Zhao et al.) - pure isotropic model, used also as background for 3D models
      call model_jp1d(myrank,r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
                      .true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                      R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST)

    case(REFERENCE_MODEL_SEA1D)
      ! SEA1D (by Lebedev & Nolet) - pure isotropic model, used also as background for 3D models
      call model_sea1d(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,SEA1DM_V)

    case default
      stop 'unknown 1D reference Earth model in meshfem3D_models_get1D_val()'

  end select

  ! needs to set vpv,vph,vsv,vsh and eta_aniso for isotropic models
  if( .not. TRANSVERSE_ISOTROPY ) then
     ! in the case of s362iso we want to save the anisotropic constants for the Voight average
     if(.not. (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF .and. ISOTROPIC_3D_MANTLE)) then
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0
     endif
  endif ! TRANSVERSE_ISOTROPY

  end subroutine meshfem3D_models_get1D_val


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho,dvp,&
                              vpv,vph,vsv,vsh,eta_aniso, &
                              RCMB,R670,RMOHO, &
                              xmesh,ymesh,zmesh,r, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use meshfem3D_models_par

  implicit none

  integer iregion_code
  double precision r_prem
  double precision rho,dvp
  double precision vpv,vph,vsv,vsh,eta_aniso

  double precision RCMB,R670,RMOHO
  double precision xmesh,ymesh,zmesh,r

  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: r_used,r_dummy,theta,phi
  double precision :: dvs,drho,vp,vs
  real(kind=4) :: xcolat,xlon,xrad,dvpv,dvph,dvsv,dvsh
  logical :: found_crust,suppress_mantle_extension

  ! initializes perturbation values
  dvs = ZERO
  dvp = ZERO
  drho = ZERO
  dvpv = 0.
  dvph = 0.
  dvsv = 0.
  dvsh = 0.
  r_used = ZERO
  suppress_mantle_extension = .false.

  ! gets point's theta/phi
  call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
  call reduce(theta,phi)

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! sets flag when mantle should not be extended to surface
  if(r_prem >= RMOHO/R_EARTH .and. .not. CRUSTAL) then
    suppress_mantle_extension = .true.
  endif

  ! gets parameters for isotropic 3D mantle model
  !
  ! note: there can be tranverse isotropy in the mantle, but only lamé parameters
  !           like kappav,kappah,muv,muh and eta_aniso are used for these simulations
  !
  ! note: in general, models here make use of perturbation values with respect to their
  !          corresponding 1-D reference models
  if( ISOTROPIC_3D_MANTLE .and. r_prem > RCMB/R_EARTH .and. .not. suppress_mantle_extension) then

    ! extend 3-D mantle model above the Moho to the surface before adding the crust
    if(r_prem > RCMB/R_EARTH .and. r_prem < RMOHO/R_EARTH) then
      ! GLL point is in mantle region, takes exact location
      r_used = r
    else ! else if(r_prem >= RMOHO/R_EARTH) then
      if( CRUSTAL ) then
        ! GLL point is above moho
        ! takes radius slightly below moho radius, this will then "extend the mantle up to the surface";
        ! crustal values will be superimposed later on
        r_used = 0.999999d0*RMOHO/R_EARTH
      endif
    endif

    ! gets model parameters
    select case( THREE_D_MODEL )

      case(THREE_D_MODEL_S20RTS)
        ! s20rts
        call mantle_s20rts(r_used,theta,phi,dvs,dvp,drho,S20RTS_V)
        vpv=vpv*(1.0d0+dvp)
        vph=vph*(1.0d0+dvp)
        vsv=vsv*(1.0d0+dvs)
        vsh=vsh*(1.0d0+dvs)
        rho=rho*(1.0d0+drho)

      case(THREE_D_MODEL_S40RTS)
        ! s40rts
        call mantle_s40rts(r_used,theta,phi,dvs,dvp,drho,S40RTS_V)
        vpv=vpv*(1.0d0+dvp)
        vph=vph*(1.0d0+dvp)
        vsv=vsv*(1.0d0+dvs)
        vsh=vsh*(1.0d0+dvs)
        rho=rho*(1.0d0+drho)

      case(THREE_D_MODEL_SEA99_JP3D)
        ! sea99 + jp3d1994
        call model_sea99_s(r_used,theta,phi,dvs,SEA99M_V)
        vsv=vsv*(1.0d0+dvs)
        vsh=vsh*(1.0d0+dvs)
        ! use Lebedev model sea99 as background and add vp & vs perturbation from Zhao 1994 model jp3d
        if(theta>=(PI/2.d0 - LAT_MAX*DEGREES_TO_RADIANS) .and. theta<=(PI/2.d0 - LAT_MIN*DEGREES_TO_RADIANS) &
            .and. phi>=LON_MIN*DEGREES_TO_RADIANS .and. phi<=LON_MAX*DEGREES_TO_RADIANS) then
          if(r_used > (R_EARTH - DEP_MAX*1000.d0)/R_EARTH) then
            call model_jp3d_iso_zhao(r_used,theta,phi,vp,vs,dvp,dvs,rho,found_crust,JP3DM_V)
            vpv=vpv*(1.0d0+dvp)
            vph=vph*(1.0d0+dvp)
            vsv=vsv*(1.0d0+dvs)
            vsh=vsh*(1.0d0+dvs)
          endif
        endif

      case(THREE_D_MODEL_SEA99)
        ! sea99 Vs-only
        call model_sea99_s(r_used,theta,phi,dvs,SEA99M_V)
        vsv=vsv*(1.0d0+dvs)
        vsh=vsh*(1.0d0+dvs)

      case(THREE_D_MODEL_JP3D)
        ! jp3d1994
        if(theta>=(PI/2.d0 - LAT_MAX*DEGREES_TO_RADIANS) .and. theta<=(PI/2.d0 - LAT_MIN*DEGREES_TO_RADIANS) &
            .and. phi>=LON_MIN*DEGREES_TO_RADIANS .and. phi<=LON_MAX*DEGREES_TO_RADIANS) then
          if(r_used > (R_EARTH - DEP_MAX*1000.d0)/R_EARTH) then
            call model_jp3d_iso_zhao(r_used,theta,phi,vp,vs,dvp,dvs,rho,found_crust,JP3DM_V)
            vpv=vpv*(1.0d0+dvp)
            vph=vph*(1.0d0+dvp)
            vsv=vsv*(1.0d0+dvs)
            vsh=vsh*(1.0d0+dvs)
          endif
        endif

      case(THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
           THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
        ! 3D Harvard models s362ani, s362wmani, s362ani_prem and s2.9ea
        xcolat = sngl(theta*180.0d0/PI)
        xlon = sngl(phi*180.0d0/PI)
        xrad = sngl(r_used*R_EARTH_KM)
        call model_s362ani_subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv, &
                    numker,numhpa,numcof,ihpa,lmax,nylm, &
                    lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
                    nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
                    coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr)

        ! to use speed values from the 1D reference model but with 3D mesh variations
        if( USE_1D_REFERENCE ) then
          ! sets all 3D variations in the mantle to zero
          dvpv = 0.d0
          dvph = 0.d0
          dvsv = 0.d0
          dvsh = 0.d0
        endif

        if(TRANSVERSE_ISOTROPY) then
          vpv=vpv*(1.0d0+dble(dvpv))
          vph=vph*(1.0d0+dble(dvph))
          vsv=vsv*(1.0d0+dble(dvsv))
          vsh=vsh*(1.0d0+dble(dvsh))
        else
          vpv=vpv+dvpv
          vph=vph+dvph
          vsv=vsv+dvsv
          vsh=vsh+dvsh
          vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                    + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
          vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                    + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
          vpv=vp
          vph=vp
          vsv=vs
          vsh=vs
          eta_aniso=1.0d0
        endif

      case(THREE_D_MODEL_PPM )
        ! point profile model
        call model_PPM(r_used,theta,phi,dvs,dvp,drho,PPM_V)
        vpv=vpv*(1.0d0+dvp)
        vph=vph*(1.0d0+dvp)
        vsv=vsv*(1.0d0+dvs)
        vsh=vsh*(1.0d0+dvs)
        rho=rho*(1.0d0+drho)

      case(THREE_D_MODEL_GAPP2 )
        ! 3D GAP model (Obayashi)
        call mantle_gapmodel(r_used,theta,phi,dvs,dvp,drho)
        vpv=vpv*(1.0d0+dvp)
        vph=vph*(1.0d0+dvp)
        vsv=vsv*(1.0d0+dvs)
        vsh=vsh*(1.0d0+dvs)
        rho=rho*(1.0d0+drho)

      case default
        stop 'unknown 3D Earth model in meshfem3D_models_get3Dmntl_val() '

    end select ! THREE_D_MODEL

  endif ! ISOTROPIC_3D_MANTLE

  ! heterogen model
  if( HETEROGEN_3D_MANTLE .and. .not. suppress_mantle_extension ) then
    call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_used,theta,phi)
    call reduce(theta,phi)
    call model_heterogen_mantle(r_used,theta,phi,dvs,dvp,drho,HMM)
    vpv=vpv*(1.0d0+dvp)
    vph=vpv*(1.0d0+dvp)
    vsv=vsv*(1.0d0+dvs)
    vsh=vsh*(1.0d0+dvs)
    rho=rho*(1.0d0+drho)
  endif ! HETEROGEN_3D_MANTLE

  if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) &
    call model_aniso_inner_core(r_prem,c11,c33,c12,c13,c44,REFERENCE_1D_MODEL, &
                                vpv,vph,vsv,vsh,rho,eta_aniso)

  if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

    ! anisotropic model between the Moho and 670 km (change to CMB if desired)
    if( r_prem > R670/R_EARTH .and. .not. suppress_mantle_extension ) then

      ! extend 3-D mantle model above the Moho to the surface before adding the crust
      if( r_prem < RMOHO/R_EARTH) then
        r_used = r_prem
      else
        if( CRUSTAL ) then
          ! fills 3-D mantle model above the Moho with the values at moho depth
          r_used = RMOHO/R_EARTH
        endif
      endif
      call model_aniso_mantle(r_used,theta,phi,rho,c11,c12,c13,c14,c15,c16, &
                        c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,AMM_V)

    else
      ! fills the rest of the mantle with the isotropic model
      c11 = rho*vpv*vpv
      c12 = rho*(vpv*vpv-2.*vsv*vsv)
      c13 = c12
      c14 = 0.d0
      c15 = 0.d0
      c16 = 0.d0
      c22 = c11
      c23 = c12
      c24 = 0.d0
      c25 = 0.d0
      c26 = 0.d0
      c33 = c11
      c34 = 0.d0
      c35 = 0.d0
      c36 = 0.d0
      c44 = rho*vsv*vsv
      c45 = 0.d0
      c46 = 0.d0
      c55 = c44
      c56 = 0.d0
      c66 = c44
    endif
  endif ! ANISOTROPIC_3D_MANTLE

!> Hejun
! Assign Attenuation after get 3-D crustal model
! This is here to identify how and where to include 3D attenuation
!       if(ATTENUATION .and. ATTENUATION_3D) then
!         call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
!         call reduce(theta,phi)
!         theta_degrees = theta / DEGREES_TO_RADIANS
!         phi_degrees = phi / DEGREES_TO_RADIANS
!         tau_e(:)   = 0.0d0
!         ! Get the value of Qmu (Attenuation) dependedent on
!         ! the radius (r_prem) and idoubling flag
!         !call model_attenuation_1D_PREM(r_prem, Qmu, idoubling)
!          call model_atten3D_QRFSI12(r_prem*R_EARTH_KM,theta_degrees,phi_degrees,Qmu,QRFSI12_Q,idoubling)
!          ! Get tau_e from tau_s and Qmu
!         call model_attenuation_getstored_tau(Qmu, T_c_source, tau_s, tau_e, AM_V, AM_S, AS_V)
!       endif

  end subroutine meshfem3D_models_get3Dmntl_val

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_get3Dcrust_val(iregion_code,xmesh,ymesh,zmesh,r, &
                              vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                              c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              elem_in_crust,moho)

! returns velocities and density for points in 3D crustal region

  use meshfem3D_models_par

  implicit none

  integer :: iregion_code
  ! note: r is the exact radius (and not r_prem with tolerance)
  double precision :: xmesh,ymesh,zmesh,r
  double precision :: vpv,vph,vsv,vsh,rho,eta_aniso,dvp

  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66

  logical :: elem_in_crust
  double precision :: moho

  ! local parameters
  double precision :: r_dummy,theta,phi
  double precision :: lat,lon
  double precision :: vpc,vsc,rhoc !,vpc_eu
  double precision :: dvs
  logical :: found_crust !,found_eucrust

  ! checks if anything to do, that is, there is nothing to do
  ! for point radius smaller than deepest possible crust radius (~80 km depth)
  if( r < R_DEEPEST_CRUST ) return

  ! gets point's position theta/phi, lat/lon
  call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
  call reduce(theta,phi)

  lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  lon = phi * RADIANS_TO_DEGREES
  if( lon > 180.0d0 ) lon = lon - 360.0d0

!---
!
! ADD YOUR MODEL HERE
!
!---
  found_crust = .false.

  ! crustal model can vary for different 3-D models
  select case (THREE_D_MODEL )

    case(THREE_D_MODEL_SEA99_JP3D,THREE_D_MODEL_JP3D)
      ! tries to use Zhao's model of the crust
      if(theta>=(PI/2.d0 - LAT_MAX*DEGREES_TO_RADIANS) .and. theta<=(PI/2.d0 - LAT_MIN*DEGREES_TO_RADIANS) &
        .and. phi>=LON_MIN*DEGREES_TO_RADIANS .and. phi<=LON_MAX*DEGREES_TO_RADIANS) then
        ! makes sure radius is fine
        if(r > (R_EARTH - DEP_MAX*1000.d0)/R_EARTH) then
          call model_jp3d_iso_zhao(r,theta,phi,vpc,vsc,dvp,dvs,rhoc,found_crust,JP3DM_V)
        endif
      else
        ! default crust
        call meshfem3D_model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)
      endif

    case ( THREE_D_MODEL_PPM )
      ! takes vs,rho from default crust
      call meshfem3D_model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)

      ! takes vp from eucrust07
      !call model_eucrust(lat,lon,r,vpc_eu,found_eucrust,EUCM_V)
      !if( found_eucrust) then
      !  vpc=vpc_eu
      !endif

    case default
      ! default crust
      call meshfem3D_model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)

  end select

  ! sets crustal values
  if( found_crust ) then
    vpv=vpc
    vph=vpc
    vsv=vsc
    vsh=vsc
    rho=rhoc
    eta_aniso=1.0d0

    ! sets anisotropy in crustal region as well
    if( ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      c11 = rho*vpv*vpv
      c12 = rho*(vpv*vpv-2.*vsv*vsv)
      c13 = c12
      c14 = 0.
      c15 = 0.
      c16 = 0.
      c22 = c11
      c23 = c12
      c24 = 0.
      c25 = 0.
      c26 = 0.
      c33 = c11
      c34 = 0.
      c35 = 0.
      c36 = 0.
      c44 = rho*vsv*vsv
      c45 = 0.
      c46 = 0.
      c55 = c44
      c56 = 0.
      c66 = c44
    endif
  endif

  end subroutine meshfem3D_models_get3Dcrust_val

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)

! returns velocity/density for default crust

  use meshfem3D_models_par
  use mpi

  implicit none

  double precision,intent(in) :: lat,lon,r
  double precision,intent(out) :: vpc,vsc,rhoc
  double precision,intent(out) :: moho
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust

  ! initializes
  vpc = 0.d0
  vsc = 0.d0
  rhoc = 0.d0
  moho = 0.d0
  found_crust = .false.

!---
!
! ADD YOUR MODEL HERE
!
!---

  select case (ITYPE_CRUSTAL_MODEL )

    case (ICRUST_CRUST2)
      ! crust 2.0
      call model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,CM_V,elem_in_crust)

    case (ICRUST_CRUSTMAPS)
      ! general crustmaps
      call model_crustmaps(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,GC_V,elem_in_crust)

    case (ICRUST_EPCRUST)
!      call model_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,CM_V,elem_in_crust)
      ! within EPCRUST region
!      if (lat >= EPCRUST_LAT_MIN .and. lat <= EPCRUST_LAT_MAX &
!          .and. lon >= EPCRUST_LON_MIN .and. lon<=EPCRUST_LON_MAX ) then
!          vpc=0.0d0
!          vsc=0.0d0
!          rhoc=0.0d0
!          moho=0.0d0
!          found_crust = .false.
          call model_epcrust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,EPCRUST,elem_in_crust)
!      endif

    case default
      stop 'crustal model type not defined'

  end select


  end subroutine meshfem3D_model_crust

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_getatten_val(idoubling,xmesh,ymesh,zmesh,r_prem, &
                              tau_e,tau_s,T_c_source, &
                              moho,Qmu,Qkappa,elem_in_crust)

! sets attenuation values tau_e and Qmu for a given point
!
! note:  only Qmu attenuation considered, Qkappa attenuation not used so far in solver...

  use meshfem3D_models_par

  implicit none

  integer idoubling

  double precision xmesh,ymesh,zmesh

  double precision r_prem
  double precision moho

  ! attenuation values
  double precision :: Qkappa,Qmu
  double precision, dimension(N_SLS) :: tau_s, tau_e
  double precision  :: T_c_source

  logical elem_in_crust

  ! local parameters
  double precision r_dummy,theta,phi,theta_degrees,phi_degrees
  double precision r_used
  double precision, parameter :: rmoho_prem = 6371.d0 - 24.4d0

  ! initializes
  tau_e(:)   = 0.0d0

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! Get the value of Qmu (Attenuation) dependent on
  ! the radius (r_prem) and idoubling flag
  if (ATTENUATION_3D) then
    ! used for models: s362ani_3DQ, s362iso_3DQ, 3D_attenuation

    ! gets spherical coordinates
    call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
    call reduce(theta,phi)
    theta_degrees = theta / DEGREES_TO_RADIANS
    phi_degrees = phi / DEGREES_TO_RADIANS

    ! in case models incorporate a 3D crust, attenuation values for mantle
    ! get expanded up to surface, and for the crustal points Qmu for PREM crust is imposed
    r_used = r_prem*R_EARTH_KM
    if( CRUSTAL ) then
      if ( r_prem > (ONE-moho) .or. elem_in_crust) then
        ! points in actual crust: puts point radius into prem crust
        r_used = rmoho_prem*1.0001
      else if( r_prem*R_EARTH_KM >= rmoho_prem ) then
        ! points below actual crust (e.g. oceanic crust case), but above prem moho:
        ! puts point slightly below prem moho to expand mantle values at that depth
        r_used = rmoho_prem*0.99999
      endif
    endif ! CRUSTAL

    ! gets value according to radius/theta/phi location and idoubling flag
    call model_atten3D_QRFSI12(r_used,theta_degrees,phi_degrees,Qmu,QRFSI12_Q,idoubling)

  else

    select case (REFERENCE_1D_MODEL)

      ! case(REFERENCE_MODEL_PREM)
      ! this case is probably not needed since Qmu is 600. between R80 and surface
      !   call model_attenuation_1D_PREM(r_prem, Qmu)

      case(REFERENCE_MODEL_1DREF)
        ! 1D Ref changes Qmu at moho depth of 24.4km
        ! we take the crustal value and assign it to points only inside actual crust,
        ! otherwise the mantle values is taken
        ! makes sense especially for points below thin oceanic and thick continental crust
        if ( CRUSTAL ) then
          ! takes crustal Q value only if point is in actual crust
          if ( r_prem > (ONE-moho) .or. elem_in_crust) then
            ! reference from 1D-REF aka STW105
            Qmu=300.0d0
            Qkappa=57822.5d0 !  not used so far...
          endif
        endif ! CRUSTAL

      case(REFERENCE_MODEL_SEA1D)
        ! SEA1D changes Qmu at 25km (moho) depth. we take the crustal value
        ! for points only inside actual crust
        if ( CRUSTAL ) then
          ! takes crustal Q value only if point is in actual crust
          if ( r_prem > (ONE-moho) .or. elem_in_crust) then
            ! reference from Sea1D
            Qmu = 300.0d0
            Qkappa = 57822.5d0  ! not used so far...
          endif
        endif

    end select

  endif

  ! Get tau_e from tau_s and Qmu
  call model_attenuation_getstored_tau(Qmu, T_c_source, tau_s, tau_e, AM_V, AM_S, AS_V)

  end subroutine meshfem3D_models_getatten_val


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_impose_val(vpv,vph,vsv,vsh,rho,dvp,eta_aniso,&
                                         myrank,iregion_code,ispec,i,j,k)

! overwrites values with updated model values (from iteration step) here, given at all GLL points

  use meshfem3D_models_par

  implicit none

  double precision :: vpv,vph,vsv,vsh,rho,dvp,eta_aniso
  integer :: myrank,iregion_code,ispec,i,j,k

  ! local parameters
  double precision :: vp,vs

  ! model GLL
  if( MGLL_V%MODEL_GLL .and. iregion_code == IREGION_CRUST_MANTLE ) then

    ! isotropic model
    if( .not. TRANSVERSE_ISOTROPY ) then

      !check
      if( ispec > size(MGLL_V%vp_new(1,1,1,:)) ) then
        call exit_MPI(myrank,'model gll: ispec too big')
      endif

      ! takes stored gll values from file
      ! ( note that these values are non-dimensionalized)
      if(CUSTOM_REAL == SIZE_REAL) then
        vp = dble( MGLL_V%vp_new(i,j,k,ispec) )
        vs = dble( MGLL_V%vs_new(i,j,k,ispec) )
        rho = dble( MGLL_V%rho_new(i,j,k,ispec) )
      else
        vp = MGLL_V%vp_new(i,j,k,ispec)
        vs = MGLL_V%vs_new(i,j,k,ispec)
        rho = MGLL_V%rho_new(i,j,k,ispec)
      endif
      ! isotropic model
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      rho = rho
      eta_aniso = 1.0d0

    ! transverse isotropic model
    else

      !check
      if( ispec > size(MGLL_V%vpv_new(1,1,1,:)) ) then
        call exit_MPI(myrank,'model gll: ispec too big')
      endif

      ! takes stored gll values from file
      if(CUSTOM_REAL == SIZE_REAL) then
        vph = dble( MGLL_V%vph_new(i,j,k,ispec) )
        vpv = dble( MGLL_V%vpv_new(i,j,k,ispec) )
        vsh = dble( MGLL_V%vsh_new(i,j,k,ispec) )
        vsv = dble( MGLL_V%vsv_new(i,j,k,ispec) )
        rho = dble( MGLL_V%rho_new(i,j,k,ispec) )
        eta_aniso = dble( MGLL_V%eta_new(i,j,k,ispec) )
      else
        vph = MGLL_V%vph_new(i,j,k,ispec)
        vpv = MGLL_V%vpv_new(i,j,k,ispec)
        vsh = MGLL_V%vsh_new(i,j,k,ispec)
        vsv = MGLL_V%vsv_new(i,j,k,ispec)
        rho = MGLL_V%rho_new(i,j,k,ispec)
        eta_aniso = MGLL_V%eta_new(i,j,k,ispec)
      endif
    endif
    ! no mantle vp perturbation
    dvp = 0.0d0

  endif ! MODEL_GLL

  end subroutine meshfem3D_models_impose_val


