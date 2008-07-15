!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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
! United States Government Sponsorship Acknowledged.

  program write_profile

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

! aniso_mantle_model_variables
  type aniso_mantle_model_variables
    sequence
    double precision beta(14,34,37,73)
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
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
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

! atten_model_QRFSI12_variables
  type atten_model_QRFSI12_variables
    sequence
    double precision dqmu(NKQ,NSQ)
    double precision spknt(NKQ)
    double precision refdepth(NDEPTHS_REFQ)
    double precision refqmu(NDEPTHS_REFQ)
  end type atten_model_QRFSI12_variables

  type (atten_model_QRFSI12_variables) QRFSI12_Q
! atten_model_QRFSI12_variables

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

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables
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
    integer :: sea99_ndep
    integer :: sea99_nlat
    integer :: sea99_nlon
    double precision :: sea99_ddeg
    double precision :: alatmin
    double precision :: alatmax
    double precision :: alonmin
    double precision :: alonmax
    double precision :: sea99_vs(100,100,100)
    double precision :: sea99_depth(100)
 end type sea99_s_model_variables
 
 type (sea99_s_model_variables) SEA99M_V
! sea99_s_model_variables

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
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


! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer idoubling


! proc numbers for MPI
  integer myrank,ier

! for loop on all the slices
  integer iregion_code

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE, &
          REFERENCE_1D_MODEL,THREE_D_MODEL,MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP

  double precision DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH


  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, HETEROGEN_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_VOLUME_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,EMULATE_ONLY

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

  integer, external :: err_occurred

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, &
               NSPEC2D_XI, &
               NSPEC2D_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
               nglob

! computed in read_compute_parameters
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

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
  integer itpspl(maxcoe,maxhpa)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)
  character(len=80) hsplfl(maxhpa)
  character(len=40) dskker(maxker)
  real(kind=4) vercof(maxker)
  real(kind=4) vercofd(maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=80) kerstr
  character(len=80) refmdl
  character(len=40) varstr(maxker)


!  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
!         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
!         NSPEC_INNER_CORE_ATTENUATION, &
!         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
!         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
!         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
!         NSPEC_CRUST_MANTLE_ADJOINT, &
!         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
!         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
!         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
!         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
!         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

! this for the different corners of the slice (which are different if the superbrick is cut)
! 1 : xi_min, eta_min
! 2 : xi_max, eta_min
! 3 : xi_max, eta_max
! 4 : xi_min, eta_max

! 1 -> min, 2 -> max

  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA


  integer i,j
  double precision rho,drhodr,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision dvp,dvs,drho
  real(kind=4) xcolat,xlon,xrad,dvpv,dvph,dvsv,dvsh
  double precision r,r_prem,r_moho,theta,phi,theta_deg,phi_deg
  double precision lat,lon,elevation
  double precision vpc,vsc,rhoc,moho
  integer NUMBER_OF_MESH_LAYERS

  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision rmin,rmax,rmax_last
! Attenuation values
  double precision, dimension(N_SLS)                     :: tau_s, tau_e
  double precision  T_c_source

  logical found_crust

  integer nit,ilayer,islice,iline,iline_icb,iline_cmb,iline_moho,iline_ocean
  integer ilayers_ocean,nlayers_ocean
  double precision delta,scaleval,r_ocean
  character(len=200) outfile

! ************** PROGRAM STARTS HERE **************
  call MPI_INIT(ier)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_profile.txt',status='unknown')
  write(IMAIN,*) 'reading parameter file..'
  print *,'reading par file'

! read the parameter file and compute additional parameters
    call read_compute_parameters(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,DT, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE,MOVIE_VOLUME_TYPE, &
          MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH,MOVIE_START,MOVIE_STOP, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,HETEROGEN_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,MOVIE_VOLUME_COARSE,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,REFERENCE_1D_MODEL,THREE_D_MODEL,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD, &
          NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
          NSPEC, &
          NSPEC2D_XI, &
          NSPEC2D_ETA, &
          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
          NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
          NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB, &
          ratio_sampling_array, ner, doubling_index,r_bottom,r_top,this_region_has_a_doubling,rmins,rmaxs,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO,CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,&
          DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,EMULATE_ONLY)

! count the total number of sources in the CMTSOLUTION file
    call count_number_of_sources(NSOURCES)


  write(IMAIN,*)
  if(ELLIPTICITY) then
    write(IMAIN,*) 'incorporating ellipticity'
  else
    write(IMAIN,*) 'no ellipticity'
  endif

  write(IMAIN,*)
  if(TOPOGRAPHY) then
    write(IMAIN,*) 'incorporating surface topography'
  else
    write(IMAIN,*) 'no surface topography'
  endif

  write(IMAIN,*)
  if(ISOTROPIC_3D_MANTLE) then
    write(IMAIN,*) 'incorporating 3-D lateral variations'
  else
    write(IMAIN,*) 'no 3-D lateral variations'
  endif

  write(IMAIN,*)
  if(CRUSTAL) then
    write(IMAIN,*) 'incorporating crustal variations'
  else
    write(IMAIN,*) 'no crustal variations'
  endif

  write(IMAIN,*)
  if(ONE_CRUST) then
    write(IMAIN,*) 'using one layer only in PREM crust'
  else
    write(IMAIN,*) 'using unmodified 1D crustal model with two layers'
  endif

  write(IMAIN,*)
  if(GRAVITY) then
    write(IMAIN,*) 'incorporating self-gravitation (Cowling approximation)'
  else
    write(IMAIN,*) 'no self-gravitation'
  endif

  write(IMAIN,*)
  if(ROTATION) then
    write(IMAIN,*) 'incorporating rotation'
  else
    write(IMAIN,*) 'no rotation'
  endif

  write(IMAIN,*)
  if(TRANSVERSE_ISOTROPY) then
    write(IMAIN,*) 'incorporating anisotropy'
  else
    write(IMAIN,*) 'no anisotropy'
  endif

  write(IMAIN,*)
  if(ATTENUATION) then
    write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
    if(ATTENUATION_3D) write(IMAIN,*)'using 3D attenuation'
  else
    write(IMAIN,*) 'no attenuation'
  endif

  write(IMAIN,*)
  if(OCEANS) then
    write(IMAIN,*) 'incorporating the oceans using equivalent load'
  else
    write(IMAIN,*) 'no oceans'
  endif

  write(IMAIN,*)

  if(ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)
  write(IMAIN,*) 'ellipticity done'

print *,'defining models'
! define models 1066a and ak135 and ref
  if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
    call define_model_1066a(CRUSTAL, M1066a_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
    call define_model_ak135(CRUSTAL, Mak135_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then
    call define_model_ref(Mref_V)
  elseif(REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
    call define_model_sea1d(CRUSTAL, SEA1DM_V)
  else
    print *,'not using ref model 1066a, ak135, ref or sea1d'
  endif
! done defining

  if(ISOTROPIC_3D_MANTLE) then
    if(THREE_D_MODEL /= 0) call read_smooth_moho
    if(THREE_D_MODEL == THREE_D_MODEL_S20RTS) then
       call read_mantle_model(D3MM_V)
    elseif(THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
           .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA) then
       call read_model_s362ani(THREE_D_MODEL,THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
                              THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA, &
                              numker,numhpa,ihpa,lmxhpa,itypehpa,ihpakern,numcoe,ivarkern,itpspl, &
                              xlaspl,xlospl,radspl,coe,hsplfl,dskker,kerstr,varstr,refmdl)
   elseif(THREE_D_MODEL == THREE_D_MODEL_SEA99_JP3D) then
       call read_sea99_s_model(SEA99M_V)
       call read_iso3d_dpzhao_model(JP3DM_V)
   elseif(THREE_D_MODEL == THREE_D_MODEL_SEA99) then
       call read_sea99_s_model(SEA99M_V)
   elseif(THREE_D_MODEL == THREE_D_MODEL_JP3D) then
       call read_iso3d_dpzhao_model(JP3DM_V)
   else
      call exit_MPI(myrank,'3D model not defined')
    endif
  endif

  if(ANISOTROPIC_3D_MANTLE) then
! the variables read are declared and stored in structure AMM_V
    call read_aniso_mantle_model(AMM_V)
  endif

  if(CRUSTAL) then
! the variables read are declared and stored in structure CM_V
    call read_crustal_model(CM_V)
  endif
  if(ANISOTROPIC_INNER_CORE) then
    call read_aniso_inner_core_model
!   one should add an MPI_BCAST here if one adds a read_aniso_inner_core_model subroutine
  endif

  if(ATTENUATION .and. ATTENUATION_3D) then
    call read_attenuation_model(MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD, AM_V)
    ! need something here to do the attenuation model setup!!
    !call attenuation_model_setup(REFERENCE_1D_MODEL, RICB, RCMB, R670, R220, R80,AM_V,M1066a_V,Mak135_V,Mref_V,AM_S,AS_V)
    call read_atten_model_3D_QRFSI12(QRFSI12_Q)
    T_c_source = AM_V%QT_c_source
    tau_s(:)   = AM_V%Qtau_s(:)
  endif

  if(ATTENUATION .and. .not. ATTENUATION_3D) then
    call read_attenuation_model(MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD, AM_V)
    call attenuation_model_setup(REFERENCE_1D_MODEL, RICB, RCMB, R670, R220, R80,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,AM_S,AS_V)
  endif
  
! read topography and bathymetry file
  if(TOPOGRAPHY .or. OCEANS) then
    call read_topo_bathy_file(ibathy_topo)
  endif

    write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
    write(IMAIN,*)
    write(IMAIN,*) 'Central cube is at a radius of ',R_CENTRAL_CUBE/1000.d0,' km'

! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * PI / 180.d0
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * PI / 180.d0
  if(NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

  if(ONE_CRUST) then
     NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS-1
  else
     NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif
!!!!!!
 do i=0,89
!  do i=0,1
  theta_deg = 1.0d0 + i*2.0d0
!   do j=0,1
  do j=0,179
   phi_deg = 1.0d0 + j*2.0d0
!  theta_deg = 90.0d0
!  phi_deg   = 270.0d0
!  theta_deg = 30.0d0
!  phi_deg   = 90.0d0

  write(*,'(a,i04.4,a,i04.4)') &
     'OUTPUT_FILES/CARDS_th',int(theta_deg),'_ph',int(phi_deg)

  write(outfile,'(a,i04.4,a,i04.4)') &
     'OUTPUT_FILES/CARDS_th',int(theta_deg),'_ph',int(phi_deg)
  open(unit=57,file=outfile,status='unknown')
!!!!!!
    rmax_last = 0.0d0
    theta = theta_deg*TWO_PI/360.0d0
    phi   = phi_deg  *TWO_PI/360.0d0
!    write(57,*) 'PREM_MODEL' 
!    write(57,*) '1 1. 1 1'
!    write(57,*) '661  124  351  647  658'
    iline = 0
    if(CRUSTAL) then
       call reduce(theta,phi)
       lat=(PI/2.0d0-theta)*180.0d0/PI
       lon=phi*180.0d0/PI
       if(lon>180.0d0) lon=lon-360.0d0
       call crustal_model(lat,lon,1.0d0,vpc,vsc,rhoc,moho,found_crust,CM_V)
       print *, 'moho depth [km]:',moho*R_EARTH_KM!, 'moho radius:',1.0d0 - moho, 'in km: ',(1.0d0-moho)*R_EARTH_KM
    endif
    if(TOPOGRAPHY .or. OCEANS) then
       call get_topo_bathy(lat,lon,elevation,ibathy_topo)
       print *, 'elevation [km]:',elevation/1000.0d0!, 'surface radius:',1.0d0 + elevation /R_EARTH
    endif
    do ilayer = 1,NUMBER_OF_MESH_LAYERS

       if(ilayer == 1) then
         rmin = 0.0d0
         rmax = rmins(NUMBER_OF_MESH_LAYERS-1)
         idoubling = IFLAG_INNER_CORE_NORMAL 
       else
         rmin = rmins(NUMBER_OF_MESH_LAYERS-ilayer+1)
         rmax = rmaxs(NUMBER_OF_MESH_LAYERS-ilayer+1)
         idoubling = doubling_index(NUMBER_OF_MESH_LAYERS-ilayer+1)
       endif
      if(CRUSTAL) then
        if(rmin == RMOHO_FICTITIOUS_IN_MESHER/R_EARTH) then
          rmin = 1.0d0 - moho
!          write(*,*) 'rmin == RMOHO',iline
        endif
        if(rmax == RMOHO_FICTITIOUS_IN_MESHER/R_EARTH) rmax = 1.0d0 - moho
       endif


      if(rmin == rmax_last) then
       if(rmin>(R_EARTH_KM-100.0d0)/R_EARTH_KM) then
         delta = 1.0d0/R_EARTH_KM
       else
         delta = 10.0d0/R_EARTH_KM
       endif
       if(TOPOGRAPHY .or. OCEANS) then
        if(rmax == 1.0d0) rmax = 1.0d0 + elevation /R_EARTH
       endif
       rmax_last = rmax
       nit = floor((rmax - rmin)/delta) + 1
       do islice = 1,nit+1
         r = rmin + (islice-1)*delta
         if(rmin == RICB/R_EARTH .and. islice == 1) iline_icb = iline
         if(rmin == RCMB/R_EARTH .and. islice == 1) iline_cmb = iline
         if(CRUSTAL) then
           if(rmin == (1.0d0 - moho) .and. islice == 1) then
              iline_moho = iline
           endif
         else
           if(rmin == RMOHO/R_EARTH .and. islice == 1) iline_moho = iline
         endif

!!start GET_MODEL
!       print *,'starting get model'
!      make sure we are within the right shell in PREM to honor discontinuities
!      use small geometrical tolerance
       r_prem = r
       if(r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
       if(r >= rmax*0.999999d0) r_prem = rmax*0.999999d0

!      get the anisotropic PREM parameters
       if(TRANSVERSE_ISOTROPY) then
         if(REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
           call prem_aniso(myrank,r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
           Qkappa,Qmu,idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
           R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then
!           print *,'calling model ref'
!           print *,r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL
!           print *,'Mref_V',Mref_V%Qkappa_ref(750),MREF_V%radius_ref(750),Mref_V%density_ref(750),Mref_V%vpv_ref(750),Mref_V%qMu_ref(750)
           call model_ref(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL,Mref_V)
!           print *,r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL
!           print *,'called model ref'
         else
           stop 'unknown 1D transversely isotropic reference Earth model in get_model'
         endif

       else
         if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
           call model_iasp91(myrank,r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
             ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME,R771,R670,R400,R220,R120,RMOHO,RMIDDLE_CRUST)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
             call prem_iso(myrank,r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
             ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
             R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
           call model_1066a(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,M1066a_V)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
           call model_ak135(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,Mak135_V)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then
           call model_ref(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL,Mref_V)
           if(.not. ISOTROPIC_3D_MANTLE) then
             vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
             vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
           endif
         else
           stop 'unknown 1D reference Earth model in get_model'
         endif

         ! in the case of s362iso we want to save the anisotropic constants for the Voight average
         if(.not. (REFERENCE_1D_MODEL == REFERENCE_MODEL_REF .and. ISOTROPIC_3D_MANTLE)) then
          vpv = vp
          vph = vp
          vsv = vs
          vsh = vs
          eta_aniso = 1.d0
         endif
       endif
!       print *, 'get 3D model'

!      get the 3-D model parameters
       if(ISOTROPIC_3D_MANTLE) then
         if(r_prem > RCMB/R_EARTH .and. r_prem < RMOHO/R_EARTH) then

           call reduce(theta,phi)
           if(THREE_D_MODEL == THREE_D_MODEL_S20RTS) then
! s20rts
             dvs = ZERO
             dvp = ZERO
             drho = ZERO
             call mantle_model(r,theta,phi,dvs,dvp,drho,D3MM_V)
             vpv=vpv*(1.0d0+dvp)
             vph=vph*(1.0d0+dvp)
             vsv=vsv*(1.0d0+dvs)
             vsh=vsh*(1.0d0+dvs)
             rho=rho*(1.0d0+drho)
           elseif(THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
                  .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA) then
! 3D Harvard models s362ani, s362wmani, s362ani_prem and s2.9ea
             dvpv = 0.
             dvph = 0.
             dvsv = 0.
             dvsh = 0.
             xcolat = sngl(theta*180.0d0/PI)
             xlon = sngl(phi*180.0d0/PI)
             xrad = sngl(r*R_EARTH_KM)
             call subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv, &
                          numker,numhpa,numcof,ihpa,lmax,nylm, &
                          lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
                          nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
                          coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr)
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
               vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
               vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
               vpv=vp
               vph=vp
               vsv=vs
               vsh=vs
               eta_aniso=1.0d0
             endif
           else
             stop 'unknown 3D Earth model in get_model'
           endif

! extend 3-D mantle model above the Moho to the surface before adding the crust
         else if(r_prem >= RMOHO/R_EARTH) then
          ! write(*,*) 'rmoho:',RMOHO,(R_EARTH-RMOHO)/1000.0d0
           call reduce(theta,phi)
           r_moho = 0.999999d0*RMOHO/R_EARTH
           if(THREE_D_MODEL == THREE_D_MODEL_S20RTS) then
! s20rts
             dvs = ZERO
             dvp = ZERO
             drho = ZERO
             call mantle_model(r_moho,theta,phi,dvs,dvp,drho,D3MM_V)
!          write(*,'(6F10.5)') r_moho,dvpv,vpv*R_EARTH*scaleval/1000.0d0
             vpv=vpv*(1.0d0+dvp)
             vph=vph*(1.0d0+dvp)
             vsv=vsv*(1.0d0+dvs)
             vsh=vsh*(1.0d0+dvs)
             rho=rho*(1.0d0+drho)
           elseif(THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
                  .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA) then
! 3D Harvard models s362ani, s362wmani, s362ani_prem and s2.9ea
             dvpv = 0.
             dvph = 0.
             dvsv = 0.
             dvsh = 0.
             xcolat = sngl(theta*180.0d0/PI)
             xlon = sngl(phi*180.0d0/PI)
             xrad = sngl(r_moho*R_EARTH_KM)
             call subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv, &
                          numker,numhpa,numcof,ihpa,lmax,nylm, &
                          lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
                          nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
                          coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr)
!          write(*,'(15F10.5)') dvpv,vpv*R_EARTH*scaleval/1000.0d0
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
               vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
               vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
               vpv=vp
               vph=vp
               vsv=vs
               vsh=vs
               eta_aniso=1.0d0
             endif
             
  else
             stop 'unknown 3D Earth model in get_model'
           endif
         endif
       endif
!       print *,'get aniso inner core'
       if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) &
           call aniso_inner_core_model(r_prem,c11,c33,c12,c13,c44,REFERENCE_1D_MODEL)
!       print *,'get aniso 3D mantle'
       if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

! anisotropic model between the Moho and 670 km (change to CMB if desired)
         if(r_prem < RMOHO/R_EARTH .and. r_prem > R670/R_EARTH) then

           call reduce(theta,phi)
           call aniso_mantle_model(r_prem,theta,phi,rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,AMM_V)
! extend 3-D mantle model above the Moho to the surface before adding the crust
         elseif(r_prem >= RMOHO/R_EARTH) then

           call reduce(theta,phi)
           r_moho = RMOHO/R_EARTH
           call aniso_mantle_model(r_moho,theta,phi,rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,AMM_V)
! fill the rest of the mantle with the isotropic model
         else
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

! This is here to identify how and where to include 3D attenuation
       if(ATTENUATION .and. ATTENUATION_3D) then
!         print *,'calling attenuation model'
         call reduce(theta,phi)
         tau_e(:)   = 0.0d0
         ! Get the value of Qmu (Attenuation) dependedent on
         ! the radius (r_prem) and idoubling flag
!         call attenuation_model_1D_PREM(r_prem, Qmu, idoubling)
!         print *,r_prem*R_EARTH_KM,theta,phi,Qmu,idoubling
         call attenuation_model_3D_QRFSI12(r_prem*R_EARTH_KM,theta,phi,Qmu,QRFSI12_Q,idoubling)
!         print *,'atten:',theta*180.0d0/PI,phi*180.0d0/PI,r_prem*R_EARTH_KM,Qmu
         ! Get tau_e from tau_s and Qmu
!         print *,'calling attenuation conversion',Qmu,T_c_source,tau_s,tau_e
        ! call attenuation_conversion(Qmu, T_c_source, tau_s, tau_e, AM_V, AM_S, AS_V)
!         print *,'done with attenuation conversion'
       endif

!      get the 3-D crustal model
       if(CRUSTAL) then
         if(r > R_DEEPEST_CRUST) then

           call reduce(theta,phi)

           lat=(PI/2.0d0-theta)*180.0d0/PI
           lon=phi*180.0d0/PI
           if(lon>180.0d0) lon=lon-360.0d0
           call crustal_model(lat,lon,r_prem,vpc,vsc,rhoc,moho,found_crust,CM_V)
!           write(*,'(a,5F10.4)') 'crust:',(1.0d0-r)*R_EARTH_KM,vpc,vsc,rhoc,moho*R_EARTH_KM
           if (found_crust) then
             vpv=vpc
             vph=vpc
             vsv=vsc
             vsh=vsc
             rho=rhoc
             eta_aniso=1.0d0
             if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
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
           endif !found_crust
         endif !r>R_DEEPEST_CRUST
       endif !CRUSTAL
 


!END GET_MODEL

      if(islice == 1) then
         r_prem = rmin
      elseif(islice == nit+1) then
         r_prem = rmax
      endif
      scaleval = dsqrt(PI*GRAV*RHOAV)
      rho = rho*RHOAV/1000.0d0
      vpv  = vpv*R_EARTH*scaleval/1000.0d0
      vph  = vph*R_EARTH*scaleval/1000.0d0
      vsv  = vsv*R_EARTH*scaleval/1000.0d0
      vsh  = vsh*R_EARTH*scaleval/1000.0d0

       iline = iline +1
!       write(57,'(i3,11F10.4)') &
!            iline,sngl(rmin*R_EARTH_KM),sngl(rmax*R_EARTH_KM),sngl(r_prem*R_EARTH_KM),sngl(r*R_EARTH_KM), &
!            sngl(vpv),sngl(vph),sngl(vsv),sngl(vsh),sngl(rho),sngl(eta_aniso),sngl(Qmu)
       write(57,'(F8.0,7F9.2,F9.5)') &
            sngl(r_prem*R_EARTH),sngl(rho*1000.d0),sngl(vpv*1000.d0),sngl(vsv*1000.d0), &
            sngl(Qkappa),sngl(Qmu),sngl(vph*1000.d0),sngl(vsh*1000.d0),sngl(eta_aniso)
    enddo !islice
   endif !rmin == rmax_last
  enddo !ilayer
  if(OCEANS .and. elevation < -500.0) then
     iline_ocean = iline
     nlayers_ocean = floor(-elevation/500.0d0)
     do ilayers_ocean=0,nlayers_ocean
        r_ocean = r_prem + ilayers_ocean*0.5d0/R_EARTH_KM
        write(57,'(F8.0,7F9.2,F9.5)') &
            sngl(r_ocean*R_EARTH),1.02,1450.,0.0,57822.5,0.0,1450.0,0.0,1.0
        iline = iline +1
     enddo
     write(57,'(F8.0,7F9.2,F9.5)') &
       sngl(1.0d0*R_EARTH),1.02,1450.,0.0,57822.5,0.0,1450.,0.0,1.0
     iline = iline+1
     write(57,*) iline,iline_icb,iline_cmb,iline_moho,iline_ocean
  else
     write(57,*) iline,iline_icb,iline_cmb,iline_moho
  endif
  enddo !sum over phi
  enddo !sum over theta
!!!!!!!!

  if(CUSTOM_REAL == SIZE_REAL) then
    write(IMAIN,*) 'using single precision for the calculations'
  else
    write(IMAIN,*) 'using double precision for the calculations'
  endif
  write(IMAIN,*)
  write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
  write(IMAIN,*)
  close(IMAIN)
  close(57)

  call MPI_FINALIZE(ier)

  end program write_profile

