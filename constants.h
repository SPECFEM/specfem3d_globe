!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
!--- user can modify parameters below
!

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
!  ALSO CHANGE FILE  precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8

! uncomment this to run in single precision
  integer, parameter :: CUSTOM_REAL = SIZE_REAL
! uncomment this to run in double precision (increases memory size by 2)
! integer, parameter :: CUSTOM_REAL = SIZE_DOUBLE

! to use restart files
  logical, parameter :: USE_RESTART_FILES = .false.
  logical, parameter :: FIRST_PART_OF_RUN = .true.

! on some processors (e.g. Pentiums) it is necessary to suppress underflows
! by using a small initial field instead of zero
  logical, parameter :: FIX_UNDERFLOW_PROBLEM = .true.

! if files on a local path on each node are also seen as global with same path
! set to .true. typically on a shared-memory machine with a common file system
! set to .false. typically on a cluster of nodes, e.g. on a Beowulf-type machine
! if running on a Beowulf-type machine, also customize global path to local
! files in create_serial_name_database.f90 ("20 format ...")
! Flag is used only when one checks the mesh with the serial codes
! ("xcheck_buffers_1D" etc.), ignore it if you do not plan to use them
  logical, parameter :: LOCAL_PATH_IS_ALSO_GLOBAL = .false.

! save AVS or OpenDX files in mesher or not
! do not use if you do not plan to use AVS or OpenDX to visualize the mesh
! because this option can create very large files
  logical, parameter :: SAVE_AVS_DX_MESH_FILES = .false.

! save a movie or not, and interval in time steps at which we save movie frames
  logical, parameter :: SAVE_AVS_DX_MOVIE = .false.
  integer, parameter :: NMOVIE = 200

! input, output and main MPI I/O files
  integer, parameter :: ISTANDARD_OUTPUT = 6
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen (slows down the code)
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT

! R_EARTH: radius of Earth (m)
! ROCEAN: radius of the ocean (m)
! RMIDDLE_CRUST: radius of the middle crust (m)
! RMOHO: radius of the Moho (m)
! R80: radius of 80 km discontinuity (m)
! R220: radius of 220 km discontinuity (m)
! R400: radius of 400 km discontinuity (m)
! R600: radius of 600 km 2nd order discontinuity (m)
! R670: radius of 670 km discontinuity (m)
! R771: radius of 771 km 2nd order discontinuity (m)
! RTOPDDOUBLEPRIME: radius of top of D" 2nd order discontinuity (m)
! RCMB: radius of CMB (m)
! RICB: radius of ICB (m)

! R_EARTH is the radius of the bottom of the oceans
  double precision, parameter :: R_EARTH = 6371000.d0
! uncomment line below for PREM with oceans
! double precision, parameter :: R_EARTH = 6368000.d0

! average density in the full Earth to normalize equation
  double precision, parameter :: RHOAV = 5514.3d0

! radii in PREM or IASPEI
! and normalized density at fluid-solid interface on fluid size for coupling

! ************************************
! ************** PREM ****************
! ************************************

  logical, parameter :: IASPEI = .false.
  double precision, parameter :: ROCEAN = 6368000.d0
  double precision, parameter :: RMIDDLE_CRUST = 6356000.d0
  double precision, parameter :: RMOHO = 6346600.d0
  double precision, parameter :: R80 = 6291000.d0
  double precision, parameter :: R220 = 6151000.d0
  double precision, parameter :: R400 = 5971000.d0
  double precision, parameter :: R600 = 5771000.d0
  double precision, parameter :: R670 = 5701000.d0
  double precision, parameter :: R771 = 5600000.d0
  double precision, parameter :: RTOPDDOUBLEPRIME = 3630000.d0
  double precision, parameter :: RCMB = 3480000.d0
  double precision, parameter :: RICB = 1221000.d0
  real(kind=CUSTOM_REAL), parameter :: RHO_TOP_OC = 9903.4384 / RHOAV
  real(kind=CUSTOM_REAL), parameter :: RHO_BOTTOM_OC = 12166.5885 / RHOAV
  real(kind=CUSTOM_REAL), parameter :: RHO_OCEANS = 1020.0 / RHOAV

! ************************************
! ************* IASPEI ***************
! ************************************

!  logical, parameter :: IASPEI = .true.
!  double precision, parameter :: RMOHO = 6341000.d0
!  double precision, parameter :: R220 = 6161000.d0
!  double precision, parameter :: R400 = 5961000.d0
!  double precision, parameter :: R600 = 5781000.d0
!  double precision, parameter :: R670 = 5711000.d0
!  double precision, parameter :: R771 = 5611000.d0
!  double precision, parameter :: RTOPDDOUBLEPRIME = 3631000.d0
!  double precision, parameter :: RCMB = 3482000.d0
!  double precision, parameter :: RICB = 1217000.d0
!  real(kind=CUSTOM_REAL), parameter :: RHO_TOP_OC = 9900.2379 / RHOAV
!  real(kind=CUSTOM_REAL), parameter :: RHO_BOTTOM_OC = 12168.6383 / RHOAV
!  real(kind=CUSTOM_REAL), parameter :: RHO_OCEANS = 1020.0 / RHOAV

! crustal model parameters for crust2.0
  integer, parameter :: NKEYS_CRUST = 359
  integer, parameter :: NLAYERS_CRUST = 8
  integer, parameter :: NCAP_CRUST = 180

! use sedimentary layers of crust 2.0
  logical, parameter :: INCLUDE_SEDIMENTS_CRUST = .true.

! minimum thickness in meters to include the effect of the oceans
! to avoid taking into account spurious oscillations in global model ETOPO
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 100.d0

! NK+1: number of radial eigenfunctions for mantle model
! NS+1: number of angular degrees for mantle model
! ND: number of discontinuities for mantle model
! uncomment only one of these two lines
! SKS12WM13
! integer, parameter :: NK = 13,NS = 12,ND = 1
! S20RTS
  integer, parameter :: NK = 20,NS = 20,ND = 1

! for topography/bathymetry model

!--- ETOPO5 5-minute model, smoothed Harvard version
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY = 4320,NY_BATHY = 2160
! resolution of topography file in minutes
  integer, parameter :: RESOLUTION_TOPO_FILE = 5

! interval at which we output time step info and max of norm of displacement
  integer, parameter :: ITAFF_TIME_STEPS = 200

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX

! number of points per surface element
  integer, parameter :: NGLLSQUARE = NGLLX * NGLLY

! number of points per spectral element
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! number of chunks (1, 3 or 6, full Earth -> six chunks)
  integer, parameter :: NCHUNKS = 6

! number of iterations in the time stepping scheme
  integer, parameter :: NITER = 2

! flag to display detailed information about location of stations
  logical, parameter :: DISPLAY_DETAILS_STATIONS = .false.

!
!--- debugging flags
!

! include central cube in the case of 6 chunks or not
! should always be set to true except when debugging code
  logical, parameter :: INCLUDE_CENTRAL_CUBE = .true.

! flags to actually assemble with MPI or not
! and to actually match fluid and solid regions of the Earth or not
! should always be set to true except when debugging code
  logical, parameter :: ACTUALLY_ASSEMBLE_MPI_SLICES = .true.
  logical, parameter :: ACTUALLY_ASSEMBLE_MPI_CHUNKS = .true.
  logical, parameter :: ACTUALLY_COUPLE_FLUID_CMB = .true.
  logical, parameter :: ACTUALLY_COUPLE_FLUID_ICB = .true.

!
!--- do NOT modify parameters below
!

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI
  double precision, parameter :: PI_OVER_TWO = PI / 2.d0,PI_OVER_FOUR = PI / 4.d0

  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0

! 3-D simulation
  integer, parameter :: NDIM = 3

! dimension of the boundaries of the slices
  integer, parameter :: NDIM2D = 2

! number of nodes for 2D and 3D shape functions for hexahedra with 27 nodes
  integer, parameter :: NGNOD = 27, NGNOD2D = 9

! gravitational constant
  double precision, parameter :: GRAV = 6.6723d-11

  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0,TWO = 2.d0,HALF = 0.5d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    TWO_THIRDS  = 2._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! declare real value independently of the machine
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_SNGL = 1.e+30_CUSTOM_REAL

  integer, parameter :: HUGEINT = 100000000

! normalized radius of free surface
  double precision, parameter :: R_UNIT_SPHERE = ONE

! same radius in km
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0

! fixed thickness of 3 km for PREM oceans
  double precision, parameter :: THICKNESS_OCEANS_PREM = 3000.d0 / R_EARTH

! non-dimensionalized size of central cube in the inner core
! This is where the central cube in the inner core and the rest of the mesh
! are matched (150 km below the ICB is optimal)
  double precision, parameter :: R_CENTRAL_CUBE = (RICB - 150000.d0) / R_EARTH

! shortest radius at which crust is implemented
  double precision, parameter :: R_DEEPEST_CRUST = (R_EARTH - 90000.d0) / R_EARTH

! define block type based upon chunk number (between 1 and 6)
! do not change this numbering, chunk AB must be number 1 for central cube
  integer, parameter :: CHUNK_AB = 1
  integer, parameter :: CHUNK_AC = 2
  integer, parameter :: CHUNK_BC = 3
  integer, parameter :: CHUNK_AC_ANTIPODE = 4
  integer, parameter :: CHUNK_BC_ANTIPODE = 5
  integer, parameter :: CHUNK_AB_ANTIPODE = 6

! maximum number of regions in the mesh
  integer, parameter :: MAX_NUM_REGIONS = 3

! define flag for regions of the global Earth mesh
  integer, parameter :: IREGION_CRUST_MANTLE = 1
  integer, parameter :: IREGION_OUTER_CORE = 2
  integer, parameter :: IREGION_INNER_CORE = 3

! define flag for elements
  integer, parameter :: IFLAG_CRUST = 1

  integer, parameter :: IFLAG_220_MOHO = 2
  integer, parameter :: IFLAG_670_220 = 3
  integer, parameter :: IFLAG_DOUBLING_670 = 4
  integer, parameter :: IFLAG_MANTLE_NORMAL = 5
  integer, parameter :: IFLAG_BOTTOM_MANTLE_LEV2 = 6
  integer, parameter :: IFLAG_BOTTOM_MANTLE = 7

  integer, parameter :: IFLAG_TOP_OUTER_CORE = 8
  integer, parameter :: IFLAG_TOP_OUTER_CORE_LEV2 = 9
  integer, parameter :: IFLAG_OUTER_CORE_NORMAL = 10
  integer, parameter :: IFLAG_BOTTOM_OUTER_CORE_LEV2 = 11
  integer, parameter :: IFLAG_BOTTOM_OUTER_CORE = 12

  integer, parameter :: IFLAG_TOP_INNER_CORE = 13
  integer, parameter :: IFLAG_TOP_INNER_CORE_LEV2 = 14
  integer, parameter :: IFLAG_INNER_CORE_NORMAL = 15

  integer, parameter :: IFLAG_IN_CENTRAL_CUBE = 16
  integer, parameter :: IFLAG_BOTTOM_CENTRAL_CUBE = 17
  integer, parameter :: IFLAG_TOP_CENTRAL_CUBE = 18
  integer, parameter :: IFLAG_IN_FICTITIOUS_CUBE = 19

! dummy flag
  integer, parameter :: IFLAG_DUMMY = 100

! define flag for regions of the global Earth for attenuation
  integer, parameter :: NUM_REGIONS_ATTENUATION = 5

  integer, parameter :: IREGION_ATTENUATION_INNER_CORE = 1
  integer, parameter :: IREGION_ATTENUATION_CMB_670 = 2
  integer, parameter :: IREGION_ATTENUATION_670_220 = 3
  integer, parameter :: IREGION_ATTENUATION_220_80 = 4
  integer, parameter :: IREGION_ATTENUATION_80_SURFACE = 5

! number of standard linear solids for attenuation
  integer, parameter :: N_SLS = 3

! define flag for buffer points
  integer, parameter :: IFLAG_NORMAL_POINT = 1
  integer, parameter :: IFLAG_MATCHING_POINT = 2

! flag for the four edges of each slice and for the bottom edge
  integer, parameter :: XI_MIN = 1
  integer, parameter :: XI_MAX = 2
  integer, parameter :: ETA_MIN = 3
  integer, parameter :: ETA_MAX = 4
  integer, parameter :: BOTTOM = 5

! number of points in each AVS or OpenDX quadrangular cell for movies
  integer, parameter :: NGNOD2D_AVS_DX = 4

! number of faces a given slice can share with other slices
! this is at most 2, except when there is only once slice per chunk
! in which case it is 4
  integer, parameter :: NUMFACES_SHARED = 4

! number of corners a given slice can share with other slices
! this is at most 1, except when there is only once slice per chunk
! in which case it is 4
  integer, parameter :: NUMCORNERS_SHARED = 4

! number of slaves per corner
  integer, parameter :: NUMSLAVES = 2

! number of elements used in radial direction in outer core doubling region
  integer, parameter :: NER_BOTTOMDBL_TOPDBL = 6*4

! for vectorization of loops
  integer, parameter :: NGLLSQUARE_NDIM = NGLLSQUARE * NDIM
  integer, parameter :: NGLLCUBE_NDIM = NGLLCUBE * NDIM

! number of layers in PREM
  integer, parameter :: NR = 640

! smallest real number on the Pentium and the SGI =  1.1754944E-38
! largest real number on the Pentium and the SGI  =  3.4028235E+38
! small negligible initial value to avoid very slow underflow trapping
! but not too small to avoid trapping on velocity and acceleration in Newmark
  real(kind=CUSTOM_REAL), parameter :: VERYSMALLVAL = 1.E-24_CUSTOM_REAL

! displacement threshold above which we consider that the code became unstable
  real(kind=CUSTOM_REAL), parameter :: STABILITY_THRESHOLD = 1.E+25_CUSTOM_REAL

! geometrical tolerance for boundary detection
  double precision, parameter :: SMALLVAL = 0.00001d0

! small tolerance for conversion from x y z to r theta phi
  double precision, parameter :: SMALL_VAL_ANGLE = 1.d-10

! geometry tolerance parameter to calculate number of independent grid points
! sensitive to actual size of model, assumes reference sphere of radius 1
! this is an absolute value for normalized coordinates in the Earth
  double precision, parameter :: SMALLVALTOL = 1.d-10

! do not use tags for MPI messages, use dummy tag instead
  integer, parameter :: itag = 0,itag2 = 0

! for the Gauss-Lobatto-Legendre points and weights
  double precision, parameter :: GAUSSALPHA = 0.d0,GAUSSBETA = 0.d0

! number of iterations to solve the non linear system for xi and eta
  integer, parameter :: NUM_ITER = 5

! number of hours per day for rotation rate of the Earth
  double precision, parameter :: HOURS_PER_DAY = 24.d0

! for lookup table for gravity every 100 m in radial direction of Earth model
  integer, parameter :: NRAD_GRAVITY = 70000

