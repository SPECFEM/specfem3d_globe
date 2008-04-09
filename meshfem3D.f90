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

  program xmeshfem3D

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

!=====================================================================!
!                                                                     !
!  meshfem3D produces a spectral element grid for the Earth.          !
!  This is accomplished based upon a mapping of the face of a cube    !
!  to a portion of the sphere (Ronchi et al., The Cubed Sphere).      !
!  Grid density is decreased by a factor of two                       !
!  three times in the radial direction.                               !
!                                                                     !
!=====================================================================!
!
! If you use this code for your own research, please cite some of these articles:
!
! @ARTICLE{KoRiTr02,
! author={D. Komatitsch and J. Ritsema and J. Tromp},
! year=2002,
! title={The Spectral-Element Method, {B}eowulf Computing, and Global Seismology},
! journal={Science},
! volume=298,
! number=5599,
! pages={1737-1742},
! doi={10.1126/science.1076024}}
!
! @ARTICLE{KoTr02a,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-I. V}alidation},
! journal={Geophys. J. Int.},
! volume=149,
! number=2,
! pages={390-412},
! doi={10.1046/j.1365-246X.2002.01653.x}}
!
! @ARTICLE{KoTr02b,
! author={D. Komatitsch and J. Tromp},
! year=2002,
! title={Spectral-Element Simulations of Global Seismic Wave Propagation{-II. 3-D} Models, Oceans, Rotation, and Self-Gravitation},
! journal={Geophys. J. Int.},
! volume=150,
! pages={303-318},
! number=1,
! doi={10.1046/j.1365-246X.2002.01716.x}}
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! year=1999,
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}
!
! If you use the kernel capabilities of the code, please cite
!
! @ARTICLE{LiTr06,
! author={Qinya Liu and Jeroen Tromp},
! title={Finite-frequency kernels based on adjoint methods},
! journal={Bull. Seismol. Soc. Am.},
! year=2006,
! volume=96,
! number=6,
! pages={2383-2397},
! doi={10.1785/0120060041}}
!
! If you use 3-D model S20RTS, please cite
!
! @ARTICLE{RiVa00,
! author={J. Ritsema and H. J. {Van Heijst}},
! year=2000,
! title={Seismic imaging of structural heterogeneity in {E}arth's mantle: Evidence for large-scale mantle flow},
! journal={Science Progress},
! volume=83,
! pages={243-259}}
!
! Reference frame - convention:
! ----------------------------
!
! The code uses the following convention for the reference frame:
!
!  - X axis is East
!  - Y axis is North
!  - Z axis is up
!
! Note that this convention is different from both the Aki-Richards convention
! and the Harvard CMT convention.
!
! Let us recall that the Aki-Richards convention is:
!
!  - X axis is North
!  - Y axis is East
!  - Z axis is down
!
! and that the Harvard CMT convention is:
!
!  - X axis is South
!  - Y axis is East
!  - Z axis is up
!
! To report bugs or suggest improvements to the code, please send an email
! to Jeroen Tromp <jtromp AT caltech.edu> and/or use our online
! bug tracking system at http://www.geodynamics.org/roundup .
!
! Evolution of the code:
! ---------------------
!
! v. 4.0 David Michea and Dimitri Komatitsch, University of Pau, France, February 2008:
!      new doubling brick in the mesh, new perfectly load-balanced mesh,
!      more flexible routines for mesh design, new inflated central cube
!      with optimized shape, far fewer mesh files saved by the mesher,
!      global arrays sorted to speed up the simulation, seismos can be
!      written by the master
! v. 3.6 Many people, many affiliations, September 2006:
!      adjoint and kernel calculations, fixed IASP91 model,
!      added AK135 and 1066a, fixed topography/bathymetry routine,
!      new attenuation routines, faster and better I/Os on very large
!      systems, many small improvements and bug fixes, new "configure"
!      script, new Pyre version, new user's manual etc.
! v. 3.5 Dimitri Komatitsch, Brian Savage and Jeroen Tromp, Caltech, July 2004:
!      any size of chunk, 3D attenuation, case of two chunks,
!      more precise topography/bathymetry model, new Par_file structure
! v. 3.4 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2003:
!      merged global and regional codes, no iterations in fluid, better movies
! v. 3.3 Dimitri Komatitsch, Caltech, September 2002:
!      flexible mesh doubling in outer core, inlined code, OpenDX support
! v. 3.2 Jeroen Tromp, Caltech, July 2002:
!      multiple sources and flexible PREM reading
! v. 3.1 Dimitri Komatitsch, Caltech, June 2002:
!      vectorized loops in solver and merged central cube
! v. 3.0 Dimitri Komatitsch and Jeroen Tromp, Caltech, May 2002:
!   ported to SGI and Compaq, double precision solver, more general anisotropy
! v. 2.3 Dimitri Komatitsch and Jeroen Tromp, Caltech, August 2001:
!                       gravity, rotation, oceans and 3-D models
! v. 2.2 Dimitri Komatitsch and Jeroen Tromp, Caltech, March 2001:
!                       final MPI package
! v. 2.0 Dimitri Komatitsch, Harvard, January 2000: MPI code for the globe
! v. 1.0 Dimitri Komatitsch, Mexico, June 1999: first MPI code for a chunk
! Jeroen Tromp, Harvard, July 1998: first chunk solver using OpenMP on Sun
! Dimitri Komatitsch, IPG Paris, December 1996: first 3-D solver for the CM5
!

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

! correct number of spectral elements in each block depending on chunk type

  integer nspec_aniso,npointot

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer, dimension(:), allocatable :: idoubling
  integer, dimension(:,:,:,:), allocatable :: ibool

! arrays with the mesh in double precision
  double precision, dimension(:,:,:,:), allocatable :: xstore,ystore,zstore

! proc numbers for MPI
  integer myrank,sizeprocs,ier

! check area and volume of the final mesh
  double precision area_local_bottom,area_total_bottom
  double precision area_local_top,area_total_top
  double precision volume_local,volume_total,volume_total_region

  integer iprocnum

! for loop on all the slices
  integer iregion_code,iregion
  integer iproc_xi,iproc_eta,ichunk

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

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
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,MOVIE_VOLUME_COARSE,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D, &
          OUTPUT_SEISMOS_ASCII_TEXT,OUTPUT_SEISMOS_SAC_ALPHANUM,OUTPUT_SEISMOS_SAC_BINARY, &
          ROTATE_SEISMOGRAMS_RT,HONOR_1D_SPHERICAL_MOHO,WRITE_SEISMOGRAMS_BY_MASTER,&
          SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE

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

! DK DK UGLY if running on MareNostrum in Barcelona
  integer :: sender, receiver, dummy1, dummy2
  integer msg_status(MPI_STATUS_SIZE)
  character(len=400) system_command

! computed in read_compute_parameters
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

! memory size of all the static arrays
  double precision :: static_memory_size

! arrays for BCAST
  integer, dimension(38) :: bcast_integer
  double precision, dimension(30) :: bcast_double_precision
  logical, dimension(26) :: bcast_logical

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

  integer :: ipass

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

! ************** PROGRAM STARTS HERE **************

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

! sizeprocs returns number of processes started (should be equal to NPROCTOT).
! myrank is the rank of each process, between 0 and NPROCTOT-1.
! as usual in MPI, process 0 is in charge of coordinating everything
! and also takes care of the main output
! do not create anything for the inner core here, will be done in solver
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! open main output file, only written to by process 0
  if(myrank == 0 .and. IMAIN /= ISTANDARD_OUTPUT) &
    open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesher.txt',status='unknown')

! get MPI starting time
  time_start = MPI_WTIME()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '****************************'
    write(IMAIN,*) '*** Specfem3D MPI Mesher ***'
    write(IMAIN,*) '****************************'
    write(IMAIN,*)
  endif

  if (myrank==0) then
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
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
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
          WRITE_SEISMOGRAMS_BY_MASTER,SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,.false.)

    if(err_occurred() /= 0) then
          call exit_MPI(myrank,'an error occurred while reading the parameter file')
    endif

! count the total number of sources in the CMTSOLUTION file
    call count_number_of_sources(NSOURCES)

    bcast_integer = (/MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
            NER_80_MOHO,NER_220_80,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
            NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
            NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,RMOHO_FICTITIOUS_IN_MESHER, &
            NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
            NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
            NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NCHUNKS,&
            SIMULATION_TYPE,REFERENCE_1D_MODEL,THREE_D_MODEL,NPROC,NPROCTOT, &
            NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube,&
            MOVIE_VOLUME_TYPE,MOVIE_START,MOVIE_STOP/)

    bcast_logical = (/TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
            CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
            TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
            RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
            SAVE_MESH_FILES,ATTENUATION, &
            ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD,CASE_3D,&
            CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,SAVE_ALL_SEISMOS_IN_ONE_FILE/)

    bcast_double_precision = (/DT,ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
            CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
            RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
            R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
            MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH/)

  endif

! broadcast the information read on the master to the nodes
    call MPI_BCAST(NSOURCES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(bcast_integer,38,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(bcast_double_precision,30,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(bcast_logical,25,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(LOCAL_PATH,150,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(MODEL,150,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(ner,MAX_NUMBER_OF_MESH_LAYERS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(ratio_sampling_array,MAX_NUMBER_OF_MESH_LAYERS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(doubling_index,MAX_NUMBER_OF_MESH_LAYERS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(r_bottom,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(r_top,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(rmins,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(rmaxs,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(rmaxs,MAX_NUMBER_OF_MESH_LAYERS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(this_region_has_a_doubling,MAX_NUMBER_OF_MESH_LAYERS,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(NSPEC,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC2D_XI,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC2D_ETA,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC2DMAX_XMIN_XMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC2DMAX_YMIN_YMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC2D_BOTTOM,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC2D_TOP,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NSPEC1D_RADIAL,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NGLOB1D_RADIAL,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NGLOB2DMAX_XMIN_XMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NGLOB2DMAX_YMIN_YMAX,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(NGLOB,MAX_NUM_REGIONS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

    call MPI_BCAST(DIFF_NSPEC1D_RADIAL,NB_SQUARE_CORNERS*NB_CUT_CASE,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(DIFF_NSPEC2D_ETA,NB_SQUARE_EDGES_ONEDIR*NB_CUT_CASE,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(DIFF_NSPEC2D_XI,NB_SQUARE_EDGES_ONEDIR*NB_CUT_CASE,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  if (myrank /=0) then

    MIN_ATTENUATION_PERIOD = bcast_integer(1)
    MAX_ATTENUATION_PERIOD = bcast_integer(2)
    NER_CRUST = bcast_integer(3)
    NER_80_MOHO = bcast_integer(4)
    NER_220_80 = bcast_integer(5)
    NER_400_220 = bcast_integer(6)
    NER_600_400 = bcast_integer(7)
    NER_670_600 = bcast_integer(8)
    NER_771_670 = bcast_integer(9)
    NER_TOPDDOUBLEPRIME_771 = bcast_integer(10)
    NER_CMB_TOPDDOUBLEPRIME = bcast_integer(11)
    NER_OUTER_CORE = bcast_integer(12)
    NER_TOP_CENTRAL_CUBE_ICB = bcast_integer(13)
    NEX_XI = bcast_integer(14)
    NEX_ETA = bcast_integer(15)
    RMOHO_FICTITIOUS_IN_MESHER = bcast_integer(16)
    NPROC_XI = bcast_integer(17)
    NPROC_ETA = bcast_integer(18)
    NTSTEP_BETWEEN_OUTPUT_SEISMOS = bcast_integer(19)
    NTSTEP_BETWEEN_READ_ADJSRC = bcast_integer(20)
    NSTEP = bcast_integer(21)
    NSOURCES = bcast_integer(22)
    NTSTEP_BETWEEN_FRAMES = bcast_integer(23)
    NTSTEP_BETWEEN_OUTPUT_INFO = bcast_integer(24)
    NUMBER_OF_RUNS = bcast_integer(25)
    NUMBER_OF_THIS_RUN = bcast_integer(26)
    NCHUNKS = bcast_integer(27)
    SIMULATION_TYPE = bcast_integer(28)
    REFERENCE_1D_MODEL = bcast_integer(29)
    THREE_D_MODEL = bcast_integer(30)
    NPROC = bcast_integer(31)
    NPROCTOT = bcast_integer(32)
    NEX_PER_PROC_XI = bcast_integer(33)
    NEX_PER_PROC_ETA = bcast_integer(34)
    ratio_divide_central_cube = bcast_integer(35)
    MOVIE_VOLUME_TYPE = bcast_integer(36)
    MOVIE_START = bcast_integer(37)
    MOVIE_STOP = bcast_integer(38)

    TRANSVERSE_ISOTROPY = bcast_logical(1)
    ANISOTROPIC_3D_MANTLE = bcast_logical(2)
    ANISOTROPIC_INNER_CORE = bcast_logical(3)
    CRUSTAL = bcast_logical(4)
    ELLIPTICITY = bcast_logical(5)
    GRAVITY = bcast_logical(6)
    ONE_CRUST = bcast_logical(7)
    ROTATION = bcast_logical(8)
    ISOTROPIC_3D_MANTLE = bcast_logical(9)
    TOPOGRAPHY = bcast_logical(10)
    OCEANS = bcast_logical(11)
    MOVIE_SURFACE = bcast_logical(12)
    MOVIE_VOLUME = bcast_logical(13)
    ATTENUATION_3D = bcast_logical(14)
    RECEIVERS_CAN_BE_BURIED = bcast_logical(15)
    PRINT_SOURCE_TIME_FUNCTION = bcast_logical(16)
    SAVE_MESH_FILES = bcast_logical(17)
    ATTENUATION = bcast_logical(18)
    ABSORBING_CONDITIONS = bcast_logical(19)
    INCLUDE_CENTRAL_CUBE = bcast_logical(20)
    INFLATE_CENTRAL_CUBE = bcast_logical(21)
    SAVE_FORWARD = bcast_logical(22)
    CASE_3D = bcast_logical(23)
    CUT_SUPERBRICK_XI = bcast_logical(24)
    CUT_SUPERBRICK_ETA = bcast_logical(25)
    SAVE_ALL_SEISMOS_IN_ONE_FILE = bcast_logical(26)

    DT = bcast_double_precision(1)
    ANGULAR_WIDTH_XI_IN_DEGREES = bcast_double_precision(2)
    ANGULAR_WIDTH_ETA_IN_DEGREES = bcast_double_precision(3)
    CENTER_LONGITUDE_IN_DEGREES = bcast_double_precision(4)
    CENTER_LATITUDE_IN_DEGREES = bcast_double_precision(5)
    GAMMA_ROTATION_AZIMUTH = bcast_double_precision(6)
    ROCEAN = bcast_double_precision(7)
    RMIDDLE_CRUST = bcast_double_precision(8)
    RMOHO = bcast_double_precision(9)
    R80 = bcast_double_precision(10)
    R120 = bcast_double_precision(11)
    R220 = bcast_double_precision(12)
    R400 = bcast_double_precision(13)
    R600 = bcast_double_precision(14)
    R670 = bcast_double_precision(15)
    R771 = bcast_double_precision(16)
    RTOPDDOUBLEPRIME = bcast_double_precision(17)
    RCMB = bcast_double_precision(18)
    RICB = bcast_double_precision(19)
    R_CENTRAL_CUBE = bcast_double_precision(20)
    RHO_TOP_OC = bcast_double_precision(21)
    RHO_BOTTOM_OC = bcast_double_precision(22)
    RHO_OCEANS = bcast_double_precision(23)
    HDUR_MOVIE = bcast_double_precision(24)
    MOVIE_TOP = bcast_double_precision(25)
    MOVIE_BOTTOM = bcast_double_precision(26)
    MOVIE_WEST = bcast_double_precision(27)
    MOVIE_EAST = bcast_double_precision(28)
    MOVIE_NORTH = bcast_double_precision(29)
    MOVIE_SOUTH = bcast_double_precision(30)

  endif

! DK DK UGLY if running on MareNostrum in Barcelona
  if(RUN_ON_MARENOSTRUM_BARCELONA) then

! check that we combine the seismograms in one large file to avoid GPFS overloading
    if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) call exit_MPI(myrank,'should use SAVE_ALL_SEISMOS_IN_ONE_FILE for GPFS in Barcelona')

! clean the local scratch space using a cascade (serial removal, one process after the other)
    if(myrank == 0) then

      receiver = myrank + 1
      call system('rm -f -r /scratch/komatits_new* > /dev/null')
      call MPI_SEND(dummy1,1,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD,ier)

    else

      sender = myrank - 1
      receiver = myrank + 1
      call MPI_RECV(dummy2,1,MPI_INTEGER,sender,itag,MPI_COMM_WORLD,msg_status,ier)
      call system('rm -f -r /scratch/komatits_new* > /dev/null')
      if(myrank < sizeprocs - 1) call MPI_SEND(dummy1,1,MPI_INTEGER,receiver,itag,MPI_COMM_WORLD,ier)

    endif

    call MPI_BARRIER(MPI_COMM_WORLD,ier)

! use the local scratch disk to save all the files, ignore the path that is given in the Par_file
    LOCAL_PATH = '/scratch/komatits_new'

! add processor name to local /scratch/komatits_new path
    write(system_command,"('_proc',i4.4)") myrank
    LOCAL_PATH = trim(LOCAL_PATH) // trim(system_command)

! create a local directory to store all the local files
    write(system_command,"('mkdir /scratch/komatits_new_proc',i4.4)") myrank
    call system(system_command)

  endif

! check that the code is running with the requested number of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes')

! dynamic allocation of mesh arrays
  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))
  allocate(ichunk_slice(0:NPROCTOT-1))
  allocate(iproc_xi_slice(0:NPROCTOT-1))
  allocate(iproc_eta_slice(0:NPROCTOT-1))

  addressing(:,:,:) = 0
  ichunk_slice(:) = 0
  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

! loop on all the chunks to create global slice addressing for solver
  if(myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/addressing.txt',status='unknown')
    write(IMAIN,*) 'creating global slice addressing'
    write(IMAIN,*)
  endif
  do ichunk = 1,NCHUNKS
    do iproc_eta=0,NPROC_ETA-1
      do iproc_xi=0,NPROC_XI-1
        iprocnum = (ichunk-1)*NPROC + iproc_eta * NPROC_XI + iproc_xi
        addressing(ichunk,iproc_xi,iproc_eta) = iprocnum
        ichunk_slice(iprocnum) = ichunk
        iproc_xi_slice(iprocnum) = iproc_xi
        iproc_eta_slice(iprocnum) = iproc_eta
        if(myrank == 0) write(IOUT,*) iprocnum,ichunk,iproc_xi,iproc_eta
      enddo
    enddo
  enddo
  if(myrank == 0) close(IOUT)

! this for the different counters (which are now different if the superbrick is cut in the outer core)
  do iregion=1,MAX_NUM_REGIONS
    NSPEC1D_RADIAL_CORNER(iregion,:) = NSPEC1D_RADIAL(iregion)
    NSPEC2D_XI_FACE(iregion,:) = NSPEC2D_XI(iregion)
    NSPEC2D_ETA_FACE(iregion,:) = NSPEC2D_ETA(iregion)
    NGLOB1D_RADIAL_CORNER(iregion,:) = NGLOB1D_RADIAL(iregion)
  enddo

  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      if (mod(iproc_xi_slice(myrank),2) == 0) then
        if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      else
        if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,3)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,3)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,3)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
        else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,4)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,4)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,4)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
        endif
      endif
    else
      if (mod(iproc_xi_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
      else
        NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
        NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
        NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
      endif
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
      else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
   NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
      endif
    endif
  endif

  if(myrank == 0) then
    write(IMAIN,*) 'This is process ',myrank
    write(IMAIN,*) 'There are ',sizeprocs,' MPI processes'
    write(IMAIN,*) 'Processes are numbered from 0 to ',sizeprocs-1
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NEX_XI,' elements along xi in each chunk'
    write(IMAIN,*) 'There are ',NEX_ETA,' elements along eta in each chunk'
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi in each chunk'
    write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta in each chunk'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices in each chunk'
    write(IMAIN,*) 'There are ',NCHUNKS,' chunks in the global mesh'
    write(IMAIN,*) 'There is a total of ',NPROCTOT,' slices in the global mesh'
    write(IMAIN,*)
    write(IMAIN,*) 'NGLLX = ',NGLLX
    write(IMAIN,*) 'NGLLY = ',NGLLY
    write(IMAIN,*) 'NGLLZ = ',NGLLZ

    write(IMAIN,*)
    write(IMAIN,*) 'Shape functions defined by NGNOD = ',NGNOD,' control nodes'
    write(IMAIN,*) 'Surface shape functions defined by NGNOD2D = ',NGNOD2D,' control nodes'
    write(IMAIN,*)
    write(IMAIN,*)
  endif

  if(myrank == 0) then

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

  endif
  if(ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

  if(ISOTROPIC_3D_MANTLE) then
    if(THREE_D_MODEL /= 0) call read_smooth_moho
    if(THREE_D_MODEL == THREE_D_MODEL_S20RTS) then
! the variables read are declared and stored in structure D3MM_V
      if(myrank == 0) call read_mantle_model(D3MM_V)
! broadcast the information read on the master to the nodes
      call MPI_BCAST(D3MM_V%dvs_a,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(D3MM_V%dvs_b,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(D3MM_V%dvp_a,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(D3MM_V%dvp_b,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(D3MM_V%spknt,NK+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(D3MM_V%qq0,(NK+1)*(NK+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(D3MM_V%qq,3*(NK+1)*(NK+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
   elseif(THREE_D_MODEL == THREE_D_MODEL_SEA99_JP3D) then
! the variables read are declared and stored in structure SEA99M_V and JP3DM_V
      if(myrank == 0) then
         call read_sea99_s_model(SEA99M_V)
         call read_iso3d_dpzhao_model(JP3DM_V)
      endif
! broadcast the information read on the master to the nodes
! SEA99M_V
      call MPI_BCAST(SEA99M_V%sea99_ndep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_nlat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_nlon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_ddeg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alatmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alatmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alonmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alonmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_vs,100*100*100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_depth,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
! JP3DM_V
      call MPI_BCAST(JP3DM_V%NPA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NRA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NHA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NPB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NRB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NHB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PNA,MPA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RNA,MRA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HNA,MHA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PNB,MPB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RNB,MRB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HNB,MHB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VELAP,MPA*MRA*MHA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VELBP,MPB*MRB*MHB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PN,51,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RRN,63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPA,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPB,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPC,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IPLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IRLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IHLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IPLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IRLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IHLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%JP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%KP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%JP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%KP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%WV,8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%P,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%R,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%H,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VP,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VS,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RA,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPJ,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
   elseif(THREE_D_MODEL == THREE_D_MODEL_SEA99) then
! the variables read are declared and stored in structure SEA99M_V
      if(myrank == 0) call read_sea99_s_model(SEA99M_V)
! broadcast the information read on the master to the nodes
! SEA99M_V
      call MPI_BCAST(SEA99M_V%sea99_ndep,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_nlat,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_nlon,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_ddeg,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alatmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alatmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alonmin,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%alonmax,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_vs,100*100*100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(SEA99M_V%sea99_depth,100,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
   elseif(THREE_D_MODEL == THREE_D_MODEL_JP3D) then
! the variables read are declared and stored in structure JP3DM_V
      if(myrank == 0) call read_iso3d_dpzhao_model(JP3DM_V)
! JP3DM_V
      call MPI_BCAST(JP3DM_V%NPA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NRA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NHA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NPB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NRB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%NHB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PNA,MPA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RNA,MRA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HNA,MHA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PNB,MPB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RNB,MRB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HNB,MHB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VELAP,MPA*MRA*MHA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VELBP,MPB*MRB*MHB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PN,51,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RRN,63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPA,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPB,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPC,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IPLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IRLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IHLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IPLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IRLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IHLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%JP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%KP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%IP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%JP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%KP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%WV,8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%P,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%R,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%H,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%PD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%HD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VP,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%VS,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%RA,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
      call MPI_BCAST(JP3DM_V%DEPJ,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
   elseif(THREE_D_MODEL == THREE_D_MODEL_S362ANI .or. THREE_D_MODEL == THREE_D_MODEL_S362WMANI &
           .or. THREE_D_MODEL == THREE_D_MODEL_S362ANI_PREM .or. THREE_D_MODEL == THREE_D_MODEL_S29EA) then
      if(myrank == 0) call read_model_s362ani(THREE_D_MODEL,THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
                              THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA, &
                              numker,numhpa,ihpa,lmxhpa,itypehpa,ihpakern,numcoe,ivarkern,itpspl, &
                              xlaspl,xlospl,radspl,coe,hsplfl,dskker,kerstr,varstr,refmdl)
  call MPI_BCAST(numker,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(numhpa,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ihpa,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lmxhpa,maxhpa,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(itypehpa,maxhpa,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ihpakern,maxker,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(numcoe,maxhpa,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ivarkern,maxker,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(itpspl,maxcoe*maxhpa,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(xlaspl,maxcoe*maxhpa,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(xlospl,maxcoe*maxhpa,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(radspl,maxcoe*maxhpa,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(coe,maxcoe*maxker,MPI_REAL,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(hsplfl,80*maxhpa,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(dskker,40*maxker,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(kerstr,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(refmdl,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(varstr,40*maxker,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)

    else
      call exit_MPI(myrank,'3D model not defined')
    endif
  endif

  if(ANISOTROPIC_3D_MANTLE) then
! the variables read are declared and stored in structure AMM_V
    if(myrank == 0) call read_aniso_mantle_model(AMM_V)
! broadcast the information read on the master to the nodes
    call MPI_BCAST(AMM_V%npar1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(AMM_V%beta,14*34*37*73,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(AMM_V%pro,47,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  endif

  if(CRUSTAL) then
! the variables read are declared and stored in structure CM_V
    if(myrank == 0) call read_crustal_model(CM_V)
! broadcast the information read on the master to the nodes
    call MPI_BCAST(CM_V%thlr,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(CM_V%velocp,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(CM_V%velocs,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(CM_V%dens,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(CM_V%abbreviation,NCAP_CRUST*NCAP_CRUST,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(CM_V%code,2*NKEYS_CRUST,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  endif

  if(ANISOTROPIC_INNER_CORE) then
    if(myrank == 0) call read_aniso_inner_core_model
!   one should add an MPI_BCAST here if one adds a read_aniso_inner_core_model subroutine
  endif

  if(ATTENUATION .and. ATTENUATION_3D) then
    if(myrank == 0) call read_attenuation_model(MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD, AM_V)

    if(myrank /= 0) allocate(AM_V%Qtau_s(N_SLS))
    call MPI_BCAST(AM_V%min_period,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%max_period,  1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%QT_c_source, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%Qtau_s(1),   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%Qtau_s(2),   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%Qtau_s(3),   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  endif

  if(ATTENUATION .and. .not. ATTENUATION_3D) then
     if(myrank == 0) call read_attenuation_model(MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD, AM_V)

    call MPI_BCAST(AM_V%min_period, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%max_period, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)

    call attenuation_model_setup(REFERENCE_1D_MODEL, RICB, RCMB, R670, R220, R80,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,AM_S,AS_V)
  endif

! read topography and bathymetry file
  if(TOPOGRAPHY .or. OCEANS) then
    if(myrank == 0) call read_topo_bathy_file(ibathy_topo)
! broadcast the information read on the master to the nodes
    call MPI_BCAST(ibathy_topo,NX_BATHY*NY_BATHY,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  endif

! get addressing for this process
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  if(myrank == 0) then
    write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
    write(IMAIN,*)
    write(IMAIN,*) 'Central cube is at a radius of ',R_CENTRAL_CUBE/1000.d0,' km'
  endif

! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * PI / 180.d0
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * PI / 180.d0
  if(NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

! volume of the slice
  volume_total = ZERO

! make sure everybody is synchronized
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

!----
!----  loop on all the regions of the mesh
!----

! number of regions in full Earth
  do iregion_code = 1,MAX_NUM_REGIONS

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '*******************************************'
    write(IMAIN,*) 'creating mesh in region ',iregion_code

    select case(iregion_code)

      case(IREGION_CRUST_MANTLE)
        write(IMAIN,*) 'this region is the crust and mantle'

      case(IREGION_OUTER_CORE)
        write(IMAIN,*) 'this region is the outer core'

      case(IREGION_INNER_CORE)
        write(IMAIN,*) 'this region is the inner core'

      case default
        call exit_MPI(myrank,'incorrect region code')
    end select
    write(IMAIN,*) '*******************************************'
    write(IMAIN,*)
  endif

! compute maximum number of points
  npointot = NSPEC(iregion_code) * NGLLX * NGLLY * NGLLZ

! use dynamic allocation to allocate memory for arrays
  allocate(idoubling(NSPEC(iregion_code)))
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))
  allocate(xstore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))
  allocate(ystore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))
  allocate(zstore(NGLLX,NGLLY,NGLLZ,NSPEC(iregion_code)))

! create all the regions of the mesh
! perform two passes in this part to be able to save memory
  do ipass = 1,2
    call create_regions_mesh(iregion_code,ibool,idoubling, &
         xstore,ystore,zstore,rmins,rmaxs, &
         iproc_xi,iproc_eta,ichunk,NSPEC(iregion_code),nspec_aniso, &
         volume_local,area_local_bottom,area_local_top, &
         nspl,rspl,espl,espl2, &
         nglob(iregion_code),npointot, &
         NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NSPEC2DMAX_XMIN_XMAX(iregion_code), &
         NSPEC2DMAX_YMIN_YMAX(iregion_code),NSPEC2D_BOTTOM(iregion_code),NSPEC2D_TOP(iregion_code), &
         ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
         NPROC_XI,NPROC_ETA,NSPEC2D_XI_FACE, &
         NSPEC2D_ETA_FACE,NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
         max(NGLOB2DMAX_XMIN_XMAX(iregion_code),NGLOB2DMAX_YMIN_YMAX(iregion_code)), &
         myrank,LOCAL_PATH,OCEANS,ibathy_topo, &
         rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
         ATTENUATION,ATTENUATION_3D,SAVE_MESH_FILES, &
         NCHUNKS,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
         R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
         ner,ratio_sampling_array,doubling_index,r_bottom, r_top,this_region_has_a_doubling,CASE_3D, &
         AMM_V,AM_V,M1066a_V,Mak135_V,Mref_V,SEA1DM_V,D3MM_V,JP3DM_V,SEA99M_V,CM_V,AM_S,AS_V, &
         numker,numhpa,numcof,ihpa,lmax,nylm, &
         lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
         nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
         coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr,ipass,ratio_divide_central_cube,HONOR_1D_SPHERICAL_MOHO, &
         CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,mod(iproc_xi_slice(myrank),2),mod(iproc_eta_slice(myrank),2))
  enddo

! store number of anisotropic elements found in the mantle
  if(nspec_aniso /= 0 .and. iregion_code /= IREGION_CRUST_MANTLE) &
    call exit_MPI(myrank,'found anisotropic elements outside of the mantle')

  if(iregion_code == IREGION_CRUST_MANTLE .and. nspec_aniso == 0) &
    call exit_MPI(myrank,'found no anisotropic elements in the mantle')

! use MPI reduction to compute total area and volume
  volume_total_region = ZERO
  area_total_bottom   = ZERO
  area_total_top   = ZERO
  call MPI_REDUCE(area_local_bottom,area_total_bottom,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ier)
  call MPI_REDUCE(area_local_top,area_total_top,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ier)
  call MPI_REDUCE(volume_local,volume_total_region,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                          MPI_COMM_WORLD,ier)

  if(myrank == 0) then
!   sum volume over all the regions
    volume_total = volume_total + volume_total_region

!   check volume of chunk, and bottom and top area

      write(IMAIN,*)
      write(IMAIN,*) '   calculated top area: ',area_total_top

! compare to exact theoretical value
    if(NCHUNKS == 6 .and. .not. TOPOGRAPHY) then
      select case(iregion_code)

        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*R_UNIT_SPHERE**2

        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2

        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2

        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
  endif

      write(IMAIN,*)
      write(IMAIN,*) 'calculated bottom area: ',area_total_bottom

! compare to exact theoretical value
  if(NCHUNKS == 6 .and. .not. TOPOGRAPHY) then

      select case(iregion_code)

        case(IREGION_CRUST_MANTLE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RCMB/R_EARTH)**2

        case(IREGION_OUTER_CORE)
          write(IMAIN,*) '            exact area: ',dble(NCHUNKS)*(4.0d0/6.0d0)*PI*(RICB/R_EARTH)**2

        case(IREGION_INNER_CORE)
          write(IMAIN,*) '            exact area (central cube): ',dble(NCHUNKS)*(2.*(R_CENTRAL_CUBE / R_EARTH)/sqrt(3.))**2

        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
  endif

  endif

  ! create chunk buffers if more than one chunk
  if(NCHUNKS > 1) then
    call create_chunk_buffers(iregion_code,NSPEC(iregion_code),ibool,idoubling,xstore,ystore,zstore, &
      nglob(iregion_code), &
      NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
      NPROC_XI,NPROC_ETA,NPROC,NPROCTOT,NGLOB1D_RADIAL_CORNER,maxval(NGLOB1D_RADIAL_CORNER(iregion_code,:)), &
      NGLOB2DMAX_XMIN_XMAX(iregion_code),NGLOB2DMAX_YMIN_YMAX(iregion_code), &
      myrank,LOCAL_PATH,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice,NCHUNKS)
  else
    if(myrank == 0) then

      write(IMAIN,*)
      write(IMAIN,*) 'only one chunk, no need to create chunk buffers'
      write(IMAIN,*)
    endif
  endif

! deallocate arrays used for that region
  deallocate(idoubling)
  deallocate(ibool)
  deallocate(xstore)
  deallocate(ystore)
  deallocate(zstore)

! make sure everybody is synchronized
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! end of loop on all the regions
  enddo

  if(myrank == 0) then
! check volume of chunk
      write(IMAIN,*)
      write(IMAIN,*) 'calculated volume: ',volume_total
    if(.not. TOPOGRAPHY) then
! take the central cube into account
! it is counted 6 times because of the fictitious elements
      if(INCLUDE_CENTRAL_CUBE) then
        write(IMAIN,*) '     exact volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)+5.*(2.*(R_CENTRAL_CUBE/R_EARTH)/sqrt(3.))**3)/6.d0
      else
        write(IMAIN,*) '     exact volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3)-(2.*(R_CENTRAL_CUBE/R_EARTH)/sqrt(3.))**3)/6.d0
      endif
    endif
  endif

!--- print number of points and elements in the mesh for each region

  if(myrank == 0) then

    numelem_crust_mantle = NSPEC(IREGION_CRUST_MANTLE)
    numelem_outer_core = NSPEC(IREGION_OUTER_CORE)
    numelem_inner_core = NSPEC(IREGION_INNER_CORE)

    numelem_total = numelem_crust_mantle + numelem_outer_core + numelem_inner_core

  write(IMAIN,*)
  write(IMAIN,*) 'Repartition of elements in regions:'
  write(IMAIN,*) '----------------------------------'
  write(IMAIN,*)
  write(IMAIN,*) 'total number of elements in each slice: ',numelem_total
  write(IMAIN,*)
  write(IMAIN,*) ' - crust and mantle: ',sngl(100.d0*dble(numelem_crust_mantle)/dble(numelem_total)),' %'
  write(IMAIN,*) ' - outer core: ',sngl(100.d0*dble(numelem_outer_core)/dble(numelem_total)),' %'
  write(IMAIN,*) ' - inner core: ',sngl(100.d0*dble(numelem_inner_core)/dble(numelem_total)),' %'
  write(IMAIN,*)
  write(IMAIN,*) 'for some mesh statistics, see comments in file OUTPUT_FILES/values_from_mesher.h'
  write(IMAIN,*)

! load balancing
  write(IMAIN,*) 'Load balancing = 100 % by definition'
  write(IMAIN,*)

  write(IMAIN,*)
  write(IMAIN,*) 'total number of time steps in the solver will be: ',NSTEP
  write(IMAIN,*)

  write(IMAIN,*)
  write(IMAIN,*) 'time-stepping of the solver will be: ',DT
  write(IMAIN,*)

! write information about precision used for floating-point operations
  if(CUSTOM_REAL == SIZE_REAL) then
    write(IMAIN,*) 'using single precision for the calculations'
  else
    write(IMAIN,*) 'using double precision for the calculations'
  endif
  write(IMAIN,*)
  write(IMAIN,*) 'smallest and largest possible floating-point numbers are: ',tiny(1._CUSTOM_REAL),huge(1._CUSTOM_REAL)
  write(IMAIN,*)

! evaluate the amount of static memory needed by the solver
  call memory_eval(OCEANS,ABSORBING_CONDITIONS,ATTENUATION,ANISOTROPIC_3D_MANTLE,&
                   TRANSVERSE_ISOTROPY,ANISOTROPIC_INNER_CORE,ROTATION,&
                   ONE_CRUST,doubling_index,this_region_has_a_doubling,&
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_sampling_array,&
                   NSPEC,nglob,SIMULATION_TYPE,MOVIE_VOLUME,SAVE_FORWARD, &
         NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
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
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,static_memory_size)

  NGLOB1D_RADIAL_TEMP(:) = &
  (/maxval(NGLOB1D_RADIAL_CORNER(1,:)),maxval(NGLOB1D_RADIAL_CORNER(2,:)),maxval(NGLOB1D_RADIAL_CORNER(3,:))/)

! create include file for the solver
  call save_header_file(NSPEC,nglob,NEX_XI,NEX_ETA,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION,ATTENUATION_3D, &
        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
        INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP, &
        static_memory_size,NGLOB1D_RADIAL_TEMP, &
        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NSPEC2D_TOP,NSPEC2D_BOTTOM, &
        NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX, &
        NPROC_XI,NPROC_ETA, &
         NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
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
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION)

  endif   ! end of section executed by main process only

! deallocate arrays used for mesh generation
  deallocate(addressing)
  deallocate(ichunk_slice)
  deallocate(iproc_xi_slice)
  deallocate(iproc_eta_slice)

! elapsed time since beginning of mesh generation
  if(myrank == 0) then
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for mesh generation and buffer creation in seconds = ',tCPU
    write(IMAIN,*) 'End of mesh generation'
    write(IMAIN,*)
! close main output file
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end program xmeshfem3D

