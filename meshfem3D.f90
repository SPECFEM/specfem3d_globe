!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 6
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!       (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================
!
! Copyright September 2006, by the California Institute of Technology.
! ALL RIGHTS RESERVED. United States Government Sponsorship Acknowledged.
!
! Any commercial use must be negotiated with the Office of Technology
! Transfer at the California Institute of Technology. This software may be
! subject to U.S. export control laws and regulations. By accepting
! this software, the user agrees to comply with all applicable U.S. export laws
! and regulations, including the International Traffic and Arms Regulations,
! 22 C.F.R. 120-130 and the Export Administration Regulations,
! 15 C.F.R. 730-744. User has the responsibility to obtain export licenses,
! or other export authority as may be required before exporting such
! information to foreign countries or providing access to foreign nationals.
! In no event shall the California Institute of Technology be liable to any
! party for direct, indirect, special, incidental or consequential damages,
! including lost profits, arising out of the use of this software and its
! documentation, even if the California Institute of Technology has been
! advised of the possibility of such damage.
!
! The California Institute of Technology specifically disclaims any
! warranties, including the implied warranties or merchantability and fitness
! for a particular purpose. The software and documentation provided hereunder
! is on an "as is" basis, and the California Institute of Technology has no
! obligations to provide maintenance, support, updates, enhancements or
! modifications.
!

  subroutine meshfem3D

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

!=====================================================================!
!                                                                     !
!  meshfem3D produces a spectral element grid for the Earth.          !
!  This is accomplished based upon a mapping of the face of a cube    !
!  to a portion of the sphere (Ronchi et al., The Cubed Sphere)       !
!  Grid density is decreased by a factor of two                       !
!  below the moho, below the 670 and at the bottom of the outer core  !
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

  external read_crustal_model,read_mantle_model,  &
       read_aniso_mantle_model,read_aniso_inner_core_model, &
       crustal_model, mantle_model,aniso_mantle_model, &
       aniso_inner_core_model,attenuation_model

! correct number of spectral elements in each block depending on chunk type

  integer nspec,nspec_aniso,nspec_aniso_mantle,nspec_aniso_mantle_all,npointot

! meshing parameters
  double precision, dimension(:), allocatable :: rns,rmins,rmaxs

! auxiliary variables to generate the mesh
  integer ix,iy,ir

  double precision xin,etan,xi,eta,rn
  double precision x,y,gamma,rgt,rgb,rg_icb
  double precision x_top,y_top,z_top
  double precision x_bot,y_bot,z_bot
  double precision x_icb,y_icb,z_icb
  double precision x_central_cube,y_central_cube,z_central_cube

  double precision, dimension(:,:,:), allocatable :: xgrid,ygrid,zgrid

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

  integer iprocnum,npx,npy

! for loop on all the slices
  integer iregion_code
  integer iproc_xi,iproc_eta,ichunk

! rotation matrix from Euler angles
  integer i,j
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

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
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NSOURCES,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,SIMULATION_TYPE,REFERENCE_1D_MODEL

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST,ROTATION,ISOTROPIC_3D_MANTLE, &
          TOPOGRAPHY,OCEANS,MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D, &
          RECEIVERS_CAN_BE_BURIED,PRINT_SOURCE_TIME_FUNCTION, &
          SAVE_MESH_FILES,ATTENUATION, &
          ABSORBING_CONDITIONS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,SAVE_FORWARD

  character(len=150) OUTPUT_FILES,LOCAL_PATH,MODEL

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

  integer, external :: err_occurred

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
               nglob_AB,nglob_AC,nglob_BC

!! DK DK UGLY if running on MareNostrum in Barcelona
  character(len=400) system_command

! ************** PROGRAM STARTS HERE **************

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

! read the parameter file
  call read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST, &
          NER_220_MOHO,NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
          NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB, &
          NER_TOP_CENTRAL_CUBE_ICB,NEX_XI,NEX_ETA,NER_DOUBLING_OUTER_CORE, &
          NPROC_XI,NPROC_ETA,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
          NTSTEP_BETWEEN_READ_ADJSRC,NSTEP,NTSTEP_BETWEEN_FRAMES, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NTSTEP_BETWEEN_OUTPUT_INFO,NUMBER_OF_RUNS, &
          NUMBER_OF_THIS_RUN,NCHUNKS,DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
          ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,CENTER_LONGITUDE_IN_DEGREES, &
          CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          R_CENTRAL_CUBE,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS,HDUR_MOVIE, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
          ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY,GRAVITY,ONE_CRUST, &
          ROTATION,ISOTROPIC_3D_MANTLE,TOPOGRAPHY,OCEANS,MOVIE_SURFACE, &
          MOVIE_VOLUME,ATTENUATION_3D,RECEIVERS_CAN_BE_BURIED, &
          PRINT_SOURCE_TIME_FUNCTION,SAVE_MESH_FILES, &
          ATTENUATION,REFERENCE_1D_MODEL,ABSORBING_CONDITIONS, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,LOCAL_PATH,MODEL,SIMULATION_TYPE,SAVE_FORWARD)

  if(err_occurred() /= 0) stop 'an error occurred while reading the parameter file'

! count the total number of sources in the CMTSOLUTION file
  call count_number_of_sources(NSOURCES)
! broadcast the information read on the master to the nodes
  call MPI_BCAST(NSOURCES,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

!! DK DK UGLY if running on MareNostrum in Barcelona
  if(RUN_ON_MARENOSTRUM_BARCELONA) then

! use the local scratch disk to save all the files, ignore the path that is given in the Par_file
    LOCAL_PATH = '/scratch/komatits'

! add processor name to local /scratch/komatits path
    write(system_command,"('_proc',i4.4)") myrank
    LOCAL_PATH = trim(LOCAL_PATH) // trim(system_command)

! create a local directory to store all the local files
    write(system_command,"('rm -r -f /scratch/komatits_proc',i4.4,' ; mkdir /scratch/komatits_proc',i4.4)") myrank,myrank
    call system(system_command)

  endif

! compute other parameters based upon values read
  call compute_parameters(NER_CRUST,NER_220_MOHO,NER_400_220, &
      NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
      NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
      NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB, &
      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NSPEC_AB,NSPEC_AC,NSPEC_BC, &
      NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
      NGLOB_AB,NGLOB_AC,NGLOB_BC,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,NCHUNKS,INCLUDE_CENTRAL_CUBE)

! check that the code is running with the requested nb of processes
  if(sizeprocs /= NPROCTOT) call exit_MPI(myrank,'wrong number of MPI processes')

! dynamic allocation of mesh arrays
  allocate(rns(0:2*NER))
  allocate(rmins(0:2*NER))
  allocate(rmaxs(0:2*NER))

  allocate(xgrid(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER))
  allocate(ygrid(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER))
  allocate(zgrid(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER))

  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))
  allocate(ichunk_slice(0:NPROCTOT-1))
  allocate(iproc_xi_slice(0:NPROCTOT-1))
  allocate(iproc_eta_slice(0:NPROCTOT-1))

! clear arrays
  xgrid(:,:,:) = 0.
  ygrid(:,:,:) = 0.
  zgrid(:,:,:) = 0.

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

! one doubling layer in outer core (block size multiple of 16)
    if(mod(NEX_PER_PROC_XI,16) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_XI must be a multiple of 16 for outer core doubling')
    if(mod(NEX_PER_PROC_ETA,16) /= 0) call exit_MPI(myrank,'NEX_PER_PROC_ETA must be a multiple of 16 for outer core doubling')

! check that number of elements per processor is the same in both directions
  if(NCHUNKS > 2 .and. NEX_PER_PROC_XI /= NEX_PER_PROC_ETA) &
    call exit_MPI(myrank,'must have the same number of elements per processor in both directions for more than two chunks')

! check that mesh doubling can be implemented
  if(NER_670_400/2<2) call exit_MPI(myrank,'NER_670_400 should be at least 2')
  if(NER_CMB_670<16) call exit_MPI(myrank,'NER_CMB_670 should be at least 16')
  if(NER_CENTRAL_CUBE_CMB<4) call exit_MPI(myrank,'NER_CENTRAL_CUBE_CMB should be at least 4')

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
    write(IMAIN,*) 'using real PREM crust with two layers'
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
    if(myrank == 0) call read_mantle_model(D3MM_V)
! broadcast the information read on the master to the nodes
    call MPI_BCAST(D3MM_V%dvs_a,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(D3MM_V%dvs_b,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(D3MM_V%dvp_a,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(D3MM_V%dvp_b,(NK+1)*(NS+1)*(NS+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(D3MM_V%spknt,NK+1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(D3MM_V%qq0,(NK+1)*(NK+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(D3MM_V%qq,3*(NK+1)*(NK+1),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  endif

  if(ANISOTROPIC_3D_MANTLE) then
! the variables read are declared and stored in module aniso_mantle_model_variables
    if(myrank == 0) call read_aniso_mantle_model(AMM_V)
! broadcast the information read on the master to the nodes
    call MPI_BCAST(AMM_V%npar1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(AMM_V%beta,14*34*37*73,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
    call MPI_BCAST(AMM_V%pro,47,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  endif

  if(CRUSTAL) then
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

    call MPI_BCAST(AM_V%min_period, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%max_period, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)

  endif
  if(ATTENUATION .AND. .NOT. ATTENUATION_3D) then
     if(myrank == 0) call read_attenuation_model(MIN_ATTENUATION_PERIOD, MAX_ATTENUATION_PERIOD, AM_V)

    call MPI_BCAST(AM_V%min_period, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
    call MPI_BCAST(AM_V%max_period, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)

    call attenuation_model_setup(REFERENCE_1D_MODEL, RICB, RCMB, R670, R220, R80,AM_V,M1066a_V,Mak135_V,AM_S,AS_V)
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
    write(IMAIN,*) 'Central cube is at a radius of ',R_CENTRAL_CUBE*R_EARTH/1000.d0,' km'
  endif

! create the radial mesh
  call mesh_radial(myrank,rns,rmins,rmaxs,NER,NER_TOP_CENTRAL_CUBE_ICB, &
      NER_CMB_TOPDDOUBLEPRIME,NER_TOPDDOUBLEPRIME_771, &
      NER_771_670,NER_670_600,NER_600_400,NER_400_220,NER_220_MOHO,NER_CRUST, &
      NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
      CRUSTAL,TOPOGRAPHY,ONE_CRUST,RMIDDLE_CRUST,R220,R400,R600,R670,R771, &
      RTOPDDOUBLEPRIME,RCMB,RICB,RMOHO,R_CENTRAL_CUBE,NCHUNKS)

! number of elements in each slice in a given chunk
  npx = 2*NEX_PER_PROC_XI
  npy = 2*NEX_PER_PROC_ETA

! compute rotation matrix from Euler angles
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * PI / 180.
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * PI / 180.
  if(NCHUNKS /= 6) call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

! fill the region between the central cube and the free surface
  do iy=0,npy
  do ix=0,npx

! full Earth (cubed sphere)

    xin=dble(ix)/dble(npx)
    xi= - (ANGULAR_WIDTH_XI_RAD/2.) + (dble(iproc_xi)+xin)*ANGULAR_WIDTH_XI_RAD/dble(NPROC_XI)
    x=dtan(xi)

    etan=dble(iy)/dble(npy)
    eta= - (ANGULAR_WIDTH_ETA_RAD/2.) + (dble(iproc_eta)+etan)*ANGULAR_WIDTH_ETA_RAD/dble(NPROC_ETA)
    y=dtan(eta)

    gamma=ONE/dsqrt(ONE+x*x+y*y)

    rgt=R_UNIT_SPHERE*gamma
    rgb=R_CENTRAL_CUBE*gamma
    rg_icb = (RICB/R_EARTH)*gamma

!   define the mesh points on the top and the bottom
!   in the six regions of the cubed shpere

    if(ichunk == CHUNK_AB) then
      x_top=-y*rgt
      y_top=x*rgt
      z_top=rgt

      x_bot=-y*rgb
      y_bot=x*rgb
      z_bot=rgb

      x_icb = -y*rg_icb
      y_icb = x*rg_icb
      z_icb = rg_icb

      if(INFLATE_CENTRAL_CUBE) then
        x_central_cube = x_bot
        y_central_cube = y_bot
        z_central_cube = z_bot
      else
        x_central_cube = -y
        y_central_cube = x
        z_central_cube = 1.
      endif

    else if(ichunk == CHUNK_AB_ANTIPODE) then
      x_top=-y*rgt
      y_top=-x*rgt
      z_top=-rgt

      x_bot=-y*rgb
      y_bot=-x*rgb
      z_bot=-rgb

      x_icb = -y*rg_icb
      y_icb = -x*rg_icb
      z_icb = -rg_icb

      if(INFLATE_CENTRAL_CUBE) then
        x_central_cube = x_bot
        y_central_cube = y_bot
        z_central_cube = z_bot
      else
        x_central_cube = -y
        y_central_cube = -x
        z_central_cube = -1.
      endif

    else if(ichunk == CHUNK_AC) then
      x_top=-y*rgt
      y_top=-rgt
      z_top=x*rgt

      x_bot=-y*rgb
      y_bot=-rgb
      z_bot=x*rgb

      x_icb = -y*rg_icb
      y_icb = -rg_icb
      z_icb = x*rg_icb

      if(INFLATE_CENTRAL_CUBE) then
        x_central_cube = x_bot
        y_central_cube = y_bot
        z_central_cube = z_bot
      else
        x_central_cube = -y
        y_central_cube = -1.
        z_central_cube = x
      endif

    else if(ichunk == CHUNK_AC_ANTIPODE) then
      x_top=-y*rgt
      y_top=rgt
      z_top=-x*rgt

      x_bot=-y*rgb
      y_bot=rgb
      z_bot=-x*rgb

      x_icb = -y*rg_icb
      y_icb = rg_icb
      z_icb = -x*rg_icb

      if(INFLATE_CENTRAL_CUBE) then
        x_central_cube = x_bot
        y_central_cube = y_bot
        z_central_cube = z_bot
      else
        x_central_cube = -y
        y_central_cube = 1.
        z_central_cube = -x
      endif

    else if(ichunk == CHUNK_BC) then
      x_top=-rgt
      y_top=y*rgt
      z_top=x*rgt

      x_bot=-rgb
      y_bot=y*rgb
      z_bot=x*rgb

      x_icb = -rg_icb
      y_icb = y*rg_icb
      z_icb = x*rg_icb

      if(INFLATE_CENTRAL_CUBE) then
        x_central_cube = x_bot
        y_central_cube = y_bot
        z_central_cube = z_bot
      else
        x_central_cube = -1.
        y_central_cube = y
        z_central_cube = x
      endif

    else if(ichunk == CHUNK_BC_ANTIPODE) then
      x_top=rgt
      y_top=-y*rgt
      z_top=x*rgt

      x_bot=rgb
      y_bot=-y*rgb
      z_bot=x*rgb

      x_icb = rg_icb
      y_icb = -y*rg_icb
      z_icb = x*rg_icb

      if(INFLATE_CENTRAL_CUBE) then
        x_central_cube = x_bot
        y_central_cube = y_bot
        z_central_cube = z_bot
      else
        x_central_cube = 1.
        y_central_cube = -y
        z_central_cube = x
      endif

    else
      call exit_MPI(myrank,'incorrect chunk numbering in meshfem3D')
    endif

  if(NCHUNKS /= 6) then

! rotate bottom
    vector_ori(1) = x_bot
    vector_ori(2) = y_bot
    vector_ori(3) = z_bot
    do i=1,3
      vector_rotated(i)=0.0d0
      do j=1,3
        vector_rotated(i)=vector_rotated(i)+rotation_matrix(i,j)*vector_ori(j)
      enddo
    enddo
    x_bot = vector_rotated(1)
    y_bot = vector_rotated(2)
    z_bot = vector_rotated(3)

! rotate top
    vector_ori(1) = x_top
    vector_ori(2) = y_top
    vector_ori(3) = z_top
    do i=1,3
      vector_rotated(i)=0.0d0
      do j=1,3
        vector_rotated(i)=vector_rotated(i)+rotation_matrix(i,j)*vector_ori(j)
      enddo
    enddo
    x_top = vector_rotated(1)
    y_top = vector_rotated(2)
    z_top = vector_rotated(3)

! rotate icb
    vector_ori(1) = x_icb
    vector_ori(2) = y_icb
    vector_ori(3) = z_icb
    do i=1,3
      vector_rotated(i)=0.0d0
      do j=1,3
        vector_rotated(i)=vector_rotated(i)+rotation_matrix(i,j)*vector_ori(j)
      enddo
    enddo
    x_icb = vector_rotated(1)
    y_icb = vector_rotated(2)
    z_icb = vector_rotated(3)

! rotate central cube
    vector_ori(1) = x_central_cube
    vector_ori(2) = y_central_cube
    vector_ori(3) = z_central_cube
    do i=1,3
      vector_rotated(i)=0.0d0
      do j=1,3
        vector_rotated(i)=vector_rotated(i)+rotation_matrix(i,j)*vector_ori(j)
      enddo
    enddo
    x_central_cube = vector_rotated(1)
    y_central_cube = vector_rotated(2)
    z_central_cube = vector_rotated(3)

  endif

! rescale central cube to match cubed sphere
  if(.not. INFLATE_CENTRAL_CUBE) then
    x_central_cube = x_central_cube * R_CENTRAL_CUBE / dsqrt(3.d0)
    y_central_cube = y_central_cube * R_CENTRAL_CUBE / dsqrt(3.d0)
    z_central_cube = z_central_cube * R_CENTRAL_CUBE / dsqrt(3.d0)
  endif

! fill the volume
  do ir=0,2*NER
    rn=rns(ir)
    xgrid(ix,iy,ir)=x_top*rn+x_bot*(ONE-rn)
    ygrid(ix,iy,ir)=y_top*rn+y_bot*(ONE-rn)
    zgrid(ix,iy,ir)=z_top*rn+z_bot*(ONE-rn)
  enddo

! modify in the inner core to match the central cube instead of a sphere
! this mesh works because even if ellipticity and/or topography are turned on
! at this stage the reference Earth is still purely spherical
! therefore it will always perfectly match the sphere defined above
    do ir=0,2*NER_TOP_CENTRAL_CUBE_ICB
      rn = dble(ir) / dble(2*NER_TOP_CENTRAL_CUBE_ICB)
      xgrid(ix,iy,ir) = x_icb * rn + x_central_cube * (ONE-rn)
      ygrid(ix,iy,ir) = y_icb * rn + y_central_cube * (ONE-rn)
      zgrid(ix,iy,ir) = z_icb * rn + z_central_cube * (ONE-rn)
    enddo

  enddo
  enddo

! volume of the slice
  volume_total = ZERO

! number of anisotropic elements in the mantle
  nspec_aniso_mantle = 0

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

! volume of bottom and top area of the slice
  volume_local = ZERO
  area_local_bottom = ZERO
  area_local_top = ZERO

! check chunk number and assign theoretical number of elements
  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    nspec = NSPEC_AB(iregion_code)
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    nspec = NSPEC_AC(iregion_code)
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    nspec = NSPEC_BC(iregion_code)
  else
    call exit_MPI(myrank,'incorrect chunk number')
  endif

! compute maximum number of points
  npointot = nspec * NGLLX * NGLLY * NGLLZ

! use dynamic allocation to allocate memory for arrays
  allocate(idoubling(nspec))
  allocate(ibool(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(ystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(zstore(NGLLX,NGLLY,NGLLZ,nspec))

! create all the regions of the mesh
  call create_regions_mesh(iregion_code,xgrid,ygrid,zgrid,ibool,idoubling, &
         xstore,ystore,zstore,npx,npy,rmins,rmaxs, &
         iproc_xi,iproc_eta,ichunk,nspec,nspec_aniso, &
         volume_local,area_local_bottom,area_local_top, &
         nspl,rspl,espl,espl2, &
         nglob_AB(iregion_code),nglob_AC(iregion_code),nglob_BC(iregion_code),npointot, &
         NER,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
         NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400,NER_400_220, &
         NER_220_MOHO,NER_CRUST,NER_DOUBLING_OUTER_CORE, &
         NSPEC2DMAX_XMIN_XMAX(iregion_code), &
         NSPEC2DMAX_YMIN_YMAX(iregion_code),NSPEC2D_BOTTOM(iregion_code),NSPEC2D_TOP(iregion_code), &
         ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE, &
         ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
         NPROC_XI,NPROC_ETA,NSPEC2D_A_XI(iregion_code),NSPEC2D_B_XI(iregion_code),NSPEC2D_C_XI(iregion_code), &
         NSPEC2D_A_ETA(iregion_code),NSPEC2D_B_ETA(iregion_code), &
         NSPEC2D_C_ETA(iregion_code),NSPEC1D_RADIAL(iregion_code),NPOIN1D_RADIAL(iregion_code), &
         myrank,LOCAL_PATH,OCEANS,ibathy_topo,NER_ICB_BOTTOMDBL, &
         crustal_model,mantle_model,aniso_mantle_model, &
         aniso_inner_core_model,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
         attenuation_model,ATTENUATION,ATTENUATION_3D,SAVE_MESH_FILES, &
         NCHUNKS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL, &
         R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN,&
         AMM_V, AM_V, M1066a_V, Mak135_V,D3MM_V,CM_V, AM_S,AS_V)

! store number of anisotropic elements found in the mantle
  if(nspec_aniso /= 0 .and. iregion_code /= IREGION_CRUST_MANTLE) &
    call exit_MPI(myrank,'found anisotropic elements outside of the mantle')

  if(iregion_code == IREGION_CRUST_MANTLE .and. nspec_aniso == 0) &
    call exit_MPI(myrank,'found no anisotropic elements in the mantle')

  if(iregion_code == IREGION_CRUST_MANTLE) nspec_aniso_mantle = nspec_aniso

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

! if only one or two chunks, create dummy values for array sizes
    if(NCHUNKS == 1) then
      nglob_AC(iregion_code) = nglob_AB(iregion_code)
      nglob_BC(iregion_code) = nglob_AB(iregion_code)
    else if(NCHUNKS == 2) then
      nglob_BC(iregion_code) = nglob_AC(iregion_code)
    endif

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
          write(IMAIN,*) '            exact area (central cube): ',dble(NCHUNKS)*(2.*R_CENTRAL_CUBE/sqrt(3.))**2

        case default
          call exit_MPI(myrank,'incorrect region code')
      end select
  endif

  endif

! create chunk buffers if more than one chunk
  if(NCHUNKS > 1) then
    call create_chunk_buffers(iregion_code,nspec,ibool,idoubling,xstore,ystore,zstore, &
      nglob_AB(iregion_code),nglob_AC(iregion_code),nglob_BC(iregion_code), &
      NSPEC2DMAX_XMIN_XMAX(iregion_code),NSPEC2DMAX_YMIN_YMAX(iregion_code), &
      NPROC_XI,NPROC_ETA,NPROC,NPROCTOT,NPOIN1D_RADIAL(iregion_code), &
      NPOIN2DMAX_XMIN_XMAX(iregion_code),NPOIN2DMAX_YMIN_YMAX(iregion_code), &
      NSPEC_AB(iregion_code),NSPEC_AC(iregion_code),NSPEC_BC(iregion_code), &
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

! compute the maximum number of anisotropic elements found in all the slices
  call MPI_ALLREDUCE(nspec_aniso_mantle,nspec_aniso_mantle_all,1,MPI_INTEGER, &
              MPI_MAX,MPI_COMM_WORLD,ier)

  if(myrank == 0) then
! check volume of chunk
      write(IMAIN,*)
      write(IMAIN,*) 'calculated volume: ',volume_total
    if(.not. TOPOGRAPHY) then
! take the central cube into account
! it is counted 6 times because of the fictitious elements
      if(INCLUDE_CENTRAL_CUBE) then
        write(IMAIN,*) '     exact volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3) + 5.*(2.*R_CENTRAL_CUBE/sqrt(3.))**3) / 6.d0
      else
        write(IMAIN,*) '     exact volume: ', &
          dble(NCHUNKS)*((4.0d0/3.0d0)*PI*(R_UNIT_SPHERE**3) - (2.*R_CENTRAL_CUBE/sqrt(3.))**3) / 6.d0
      endif
    endif
  endif

!--- print number of points and elements in the mesh for each region
!--- in chunks of type AB

  if(myrank == 0) then

    numelem_crust_mantle = NSPEC_AB(IREGION_CRUST_MANTLE)
    numelem_outer_core = NSPEC_AB(IREGION_OUTER_CORE)
    numelem_inner_core = NSPEC_AB(IREGION_INNER_CORE)

    numelem_total = numelem_crust_mantle + numelem_outer_core + numelem_inner_core

  write(IMAIN,*)
  write(IMAIN,*) 'Repartition of elements in regions:'
  write(IMAIN,*) '----------------------------------'
  write(IMAIN,*)
  write(IMAIN,*) 'total number of elements in each slice in type AB chunk: ',numelem_total
  write(IMAIN,*)
  write(IMAIN,*) ' - crust and mantle: ',sngl(100.d0*dble(numelem_crust_mantle)/dble(numelem_total)),' %'
  write(IMAIN,*) ' - outer core: ',sngl(100.d0*dble(numelem_outer_core)/dble(numelem_total)),' %'
  write(IMAIN,*) ' - inner core: ',sngl(100.d0*dble(numelem_inner_core)/dble(numelem_total)),' %'
  write(IMAIN,*)
  write(IMAIN,*) 'for some mesh statistics, see comments in file OUTPUT_FILES/values_from_mesher.h'
  write(IMAIN,*)

! load balancing with respect to chunk of type BC (largest number of elements)
  if(NCHUNKS > 1) then
    write(IMAIN,*) 'Load balancing AB blocks = ',100.d0*dble(sum(nspec_AB))/dble(sum(nspec_BC)),' %'
    write(IMAIN,*) 'Load balancing AC blocks = ',100.d0*dble(sum(nspec_AC))/dble(sum(nspec_BC)),' %'
    write(IMAIN,*) 'Load balancing BC blocks = 100 % by definition'
  else
    write(IMAIN,*) 'Load balancing = 100 % by definition for one chunk'
  endif
  write(IMAIN,*)

  write(IMAIN,*)
  write(IMAIN,*) 'total number of time steps in the solver will be: ',NSTEP
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

! create include file for the solver
  call save_header_file(NSPEC_AB,NSPEC_AC,NSPEC_BC,nglob_AB,nglob_AC,nglob_BC, &
        NEX_XI,NEX_ETA,nspec_aniso_mantle_all,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION,ATTENUATION_3D, &
        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
        INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP)

  endif   ! end of section executed by main process only

! deallocate arrays used for mesh generation
  deallocate(rns)
  deallocate(rmins)
  deallocate(rmaxs)

  deallocate(xgrid)
  deallocate(ygrid)
  deallocate(zgrid)

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
  endif

! close main output file
  if(myrank == 0) then
    write(IMAIN,*) 'done'
    write(IMAIN,*)
    close(IMAIN)
  endif

! synchronize all the processes to make sure everybody has finished
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  end subroutine meshfem3D

