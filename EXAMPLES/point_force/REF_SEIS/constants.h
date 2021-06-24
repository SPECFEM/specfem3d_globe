!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

! setup/constants.h.  Generated from constants.h.in by configure.

!
!--- user can modify parameters below
!

!
! solver in single or double precision depending on the machine (4 or 8 bytes)
!
!  ALSO CHANGE FILE precision.h ACCORDINGLY
!
  integer, parameter :: SIZE_REAL = 4, SIZE_DOUBLE = 8

! usually the size of integer and logical variables is the same as regular single-precision real variable
  integer, parameter :: SIZE_INTEGER = SIZE_REAL
  integer, parameter :: SIZE_LOGICAL = SIZE_REAL

! set to SIZE_REAL to run in single precision
! set to SIZE_DOUBLE to run in double precision (increases memory size by 2)
  integer, parameter :: CUSTOM_REAL = SIZE_REAL

!*********************************************************************************************************

!! DK DK added this temporarily here to make SPECFEM3D and SPECFEM3D_GLOBE much more similar
!! DK DK in terms of the structure of their main time iteration loop; these are future features
!! DK DK that are missing in this code but implemented in the other and that could thus be cut and pasted one day
  logical, parameter :: OUTPUT_ENERGY = .false.
  integer, parameter :: IOUT_ENERGY = 937  ! file number for the energy file
  logical, parameter :: GRAVITY_SIMULATION = .false.

! if files on a local path on each node are also seen as global with same path
! set to .true. typically on a machine with a common (shared) file system, e.g. LUSTRE, GPFS or NFS-mounted /home
! set to .false. typically on a cluster of nodes with local disks used to store the mesh (not very common these days)
! if running on a cluster of nodes with local disks, also customize global path
! to local files in create_serial_name_database.f90 ("20 format ...")
! Flag is used only when one checks the mesh with the serial codes
! ("xcheck_buffers_1D" etc.), ignore it if you do not plan to use them
  logical, parameter :: LOCAL_PATH_IS_ALSO_GLOBAL = .true.

! maximum length of strings used for paths, reading from files, etc.
  integer, parameter :: MAX_STRING_LEN = 512

! input, output and main MPI I/O files
! note: careful with these unit numbers, we mostly use units in the 40-50 range.
!       cray Fortran e.g. reserves 0,5,6 (standard error,input,output units) and 100,101,102 (input,output,error unit)
  integer, parameter :: ISTANDARD_OUTPUT = 6     ! or for cray: 101
! I/O unit for file input,output
  integer, parameter :: IIN = 40,IOUT = 41
! uncomment this to write messages to a text file
  integer, parameter :: IMAIN = 42
! uncomment this to write messages to the screen (slows down the code)
! integer, parameter :: IMAIN = ISTANDARD_OUTPUT

! I/O unit for noise data files
  integer, parameter :: IIN_NOISE = 43,IOUT_NOISE = 44
! I/O unit for adjoint source files
  integer, parameter :: IIN_ADJ = 45
! local file unit for output of buffers
  integer, parameter :: IOUT_BUFFERS = 46
! I/O unit for source and receiver vtk file
  integer, parameter :: IOUT_VTK = 47
! I/O unit for sac files
  integer, parameter :: IOUT_SAC = 48

!!-----------------------------------------------------------
!!
!! Earth
!!
!!-----------------------------------------------------------
! EARTH_R is the radius of the bottom of the oceans (radius of Earth in m)
  double precision, parameter :: EARTH_R = 6371000.d0
! uncomment line below for PREM with oceans
! double precision, parameter :: EARTH_R = 6368000.d0

! radius of the Earth in km
  double precision, parameter :: EARTH_R_KM = EARTH_R / 1000.d0

! average density in the full Earth to normalize equation
  double precision, parameter :: EARTH_RHOAV = 5514.3d0

!!-----------------------------------------------------------
!!
!! for topography/bathymetry model
!!
!!-----------------------------------------------------------
!--- ETOPO5 5-minute model, smoothed Harvard version (not provided in package, see DATA/topo_bathy/ for more infos)
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY_5 = 4320,NY_BATHY_5 = 2160
! resolution of topography file in minutes
  double precision, parameter :: RESOLUTION_TOPO_FILE_5 = 5.0
! pathname of the topography file
  character (len=*), parameter :: PATHNAME_TOPO_FILE_5 = &
    'DATA/topo_bathy/topo_bathy_etopo5_smoothed_Harvard.dat'

!---  ETOPO4 4-minute model (default)
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY_4 = 5400,NY_BATHY_4 = 2700
! resolution of topography file in minutes
  double precision, parameter :: RESOLUTION_TOPO_FILE_4 = 4.0
! pathname of the topography file
! character (len=*), parameter :: PATHNAME_TOPO_FILE_4 = &
!    'DATA/topo_bathy/topo_bathy_etopo4_from_etopo2_subsampled.bin'
  character (len=*), parameter :: PATHNAME_TOPO_FILE_4 = &
    'DATA/topo_bathy/topo_bathy_etopo4_smoothed_window_7.bin'

!--- ETOPO2 2-minute model (not provided in package, see DATA/topo_bathy/ for more infos)
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY_2 = 10800,NY_BATHY_2 = 5400
! resolution of topography file in minutes
  double precision, parameter :: RESOLUTION_TOPO_FILE_2 = 2.0
! pathname of the topography file
! character (len=*), parameter :: PATHNAME_TOPO_FILE_2 = &
!    'DATA/topo_bathy/topo_bathy_etopo1_ice_c_resampled_at_2minutes_original_unmodified_unsmoothed.bin'
! character (len=*), parameter :: PATHNAME_TOPO_FILE_2 = &
!    'DATA/topo_bathy/topo_bathy_etopo1_ice_c_resampled_at_2minutes_smoothed_window_3.bin'
! character (len=*), parameter :: PATHNAME_TOPO_FILE_2 = &
!    'DATA/topo_bathy/topo_bathy_etopo2v2c_original_unmodified_unsmoothed.bin'
  character (len=*), parameter :: PATHNAME_TOPO_FILE_2 = &
    'DATA/topo_bathy/topo_bathy_etopo2v2c_smoothed_window_3.bin'

!--- ETOPO1 1-minute model (not provided in package, see DATA/topo_bathy/ for more infos)
! size of topography and bathymetry file
  integer, parameter :: NX_BATHY_1 = 21600,NY_BATHY_1 = 10800
! resolution of topography file in minutes
  double precision, parameter :: RESOLUTION_TOPO_FILE_1 = 1.0
! pathname of the topography file
!  character (len=*), parameter :: PATHNAME_TOPO_FILE_1 = &
!    'DATA/topo_bathy/topo_bathy_etopo1_ice_c_original_unmodified_unsmoothed.bin'
  character (len=*), parameter :: PATHNAME_TOPO_FILE_1 = &
    'DATA/topo_bathy/topo_bathy_etopo1_ice_c_smoothed_window_3.bin'

!--- Default
! Topography defaults to ETOPO4
  integer, parameter :: EARTH_NX_BATHY = NX_BATHY_4
  integer, parameter :: EARTH_NY_BATHY = NY_BATHY_4
  double precision, parameter :: EARTH_RESOLUTION_TOPO_FILE = RESOLUTION_TOPO_FILE_4
  character (len=*), parameter :: EARTH_PATHNAME_TOPO_FILE = PATHNAME_TOPO_FILE_4

! For the reference ellipsoid to convert geographic latitudes to geocentric:
!
! From Dahlen and Tromp (1998): "Spherically-symmetric Earth models all have the same hydrostatic
! surface ellipticity 1/299.8. This is 0.5 percent smaller than observed flattening of best-fitting ellipsoid 1/298.3.
! The discrepancy is referred to as the "excess equatorial bulge of the Earth",
! an early discovery of artificial satellite geodesy."
!
! From Paul Melchior, IUGG General Assembly, Vienna, Austria, August 1991 Union lecture,
! available at http://www.agu.org/books/sp/v035/SP035p0047/SP035p0047.pdf :
! "It turns out that the spheroidal models constructed on the basis of the spherically-symmetric models (PREM, 1066A)
! by using the Clairaut differential equation to calculate the flattening in function of the radius vector imply hydrostaticity.
! These have surface ellipticity 1/299.8 and a corresponding dynamical flattening of .0033 (PREM).
! The actual ellipticty of the Earth for a best-fitting ellipsoid is 1/298.3 with a corresponding dynamical flattening of .0034."
!
! Thus, flattening f = 1/299.8 is what is used in SPECFEM3D_GLOBE, as it should.
! And thus eccentricity squared e^2 = 1 - (1-f)^2 = 1 - (1 - 1/299.8)^2 = 0.00665998813529,
! and the correction factor used in the code to convert geographic latitudes to geocentric
! is 1 - e^2 = (1-f)^2 = (1 - 1/299.8)^2 = 0.9933400118647.
!
! As a comparison, the classical World Geodetic System reference ellipsoid WGS 84
! (see e.g. http://en.wikipedia.org/wiki/World_Geodetic_System) has f = 1/298.2572236.
  double precision, parameter :: EARTH_FLATTENING_F = 1.d0 / 299.8d0
  double precision, parameter :: EARTH_ONE_MINUS_F_SQUARED = (1.d0 - EARTH_FLATTENING_F)**2

! Use GLL points to capture TOPOGRAPHY and ELLIPTICITY (experimental feature, currently does not work, DO NOT USE)
  logical, parameter :: USE_GLL = .false.

! the code will print an error message and stop if it finds that the topography input file
! contains values outside this range.
! we take a safety margin just in case of a smoothed or modified model, which can locally create slightly different values
  integer, parameter :: EARTH_TOPO_MINIMUM = - 11200 ! (max depth in m, Mariana trench)
  integer, parameter :: EARTH_TOPO_MAXIMUM = + 9000 ! (height in m, Mount Everest)

! minimum thickness in meters to include the effect of the oceans and topo
  double precision, parameter :: MINIMUM_THICKNESS_3D_OCEANS = 50.d0

! plots pnm-image (showing used topography)
! file can become fairly big for large topo-files, e.g. ETOPO1 creates a ~2.7 GB pnm-image
  logical,parameter :: PLOT_PNM_IMAGE_TOPO_BATHY = .false.


!!-----------------------------------------------------------
!!
!! for crustal model
!!
!!-----------------------------------------------------------
!! (uncomment desired model)

! crustal model cap smoothing - extension range in degree
! note: using a smaller range, e.g. 0.5 degrees, leads to undefined Jacobian error at different places.
!       this is probably due to stretching elements below sharp gradients, especially with deep moho values.
!       so far, the only thing that works is to smooth out values and take special care of the Andes...
! TODO: one could try to adapt this degree range to the simulation resolution in the future
  double precision, parameter :: CAP_SMOOTHING_DEGREE_DEFAULT = 1.0d0

! increase smoothing for critical regions (for instance high mountains in Europe and in the Andes) to increases mesh stability
  logical, parameter :: SMOOTH_CRUST_EVEN_MORE = .true.

! use sedimentary layers in crustal model
  logical, parameter :: INCLUDE_SEDIMENTS_IN_CRUST = .true.
  logical, parameter :: INCLUDE_ICE_IN_CRUST       = .false. ! always set this to false except for gravity integral calculations
  double precision, parameter :: MINIMUM_SEDIMENT_THICKNESS = 2.d0 ! minimim thickness in km


!!-----------------------------------------------------------
!!
!! Gauss-Lobatto-Legendre resolution
!!
!!-----------------------------------------------------------

! number of GLL points in each direction of an element (degree plus one)
  integer, parameter :: NGLLX = 5
  integer, parameter :: NGLLY = NGLLX
  integer, parameter :: NGLLZ = NGLLX


!!-----------------------------------------------------------
!!
!! source/receiver setup
!!
!!-----------------------------------------------------------

! flag to exclude elements that are too far from target in source detection
  logical, parameter :: USE_DISTANCE_CRITERION = .true.

! flag to display detailed information about location of stations
  logical, parameter :: DISPLAY_DETAILS_STATIONS = .false.

! maximum length of station and network name for receivers
  integer, parameter :: MAX_LENGTH_STATION_NAME = 32
  integer, parameter :: MAX_LENGTH_NETWORK_NAME = 8

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual. This source decay rate to mimic an equivalent triangle
! was found by trial and error
  double precision, parameter :: SOURCE_DECAY_MIMIC_TRIANGLE = 1.628d0

! maximum number of sources and receivers to locate simultaneously
  integer, parameter :: NSOURCES_SUBSET_MAX = 100
  integer, parameter :: NREC_SUBSET_MAX = 200

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision, parameter :: USER_T0 = 0.0d0

! distance threshold (in km) above which we consider that a receiver
! is located outside the mesh and therefore excluded from the station list
  double precision, parameter :: THRESHOLD_EXCLUDE_STATION = 5.d0

! This parameter flags whether or not we will use an external source
! time function. Set to false by default.
  logical, parameter :: EXTERNAL_SOURCE_TIME_FUNCTION = .false.

!!-----------------------------------------------------------
!!
!! for attenuation
!!
!!-----------------------------------------------------------

! in the case of 1D attenuation models, use NGLL^3 storage of the attenuation constants in each GLL spectral element anyway
! (always safe to leave that to true; in the case of 1D attenuation models, setting it to false can save some memory storage)
  logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE  = .true.

! reference frequency of anelastic model
! by default, we use PREM values at 1 Hz
  double precision, parameter :: ATTENUATION_f0_REFERENCE = 1.d0   ! in Hz

! saves velocity model files shifted to the center frequency of the simulation attenuation period band
  logical, parameter :: ATTENUATION_SAVE_MODEL_AT_SHIFTED_CENTER_FREQ = .false.

!!-----------------------------------------------------------
!!
!! noise simulations
!!
!!-----------------------------------------------------------

! noise buffer for file i/o
! maximum noise buffer size in MB (will be used to evaluate number of steps for buffer, when UNDO_ATTENUATION == .false.)
! good performance results obtained for buffer sizes ~ 150 MB
  double precision,parameter :: MAXIMUM_NOISE_BUFFER_SIZE_IN_MB = 150.d0


!!-----------------------------------------------------------
!!
!! mesh optimization
!!
!!-----------------------------------------------------------

! maximum number of chunks (full sphere)
  integer, parameter :: NCHUNKS_MAX = 6

! Courant number for time step suggestion
! (empirical choice for distorted elements to estimate time step)
  double precision,parameter :: COURANT_SUGGESTED = 0.55d0

! number of points per minimum wavelength for minimum period estimate
  integer,parameter :: NPTS_PER_WAVELENGTH = 5

! the first doubling is implemented right below the Moho
! it seems optimal to implement the three other doublings at these depths
! in the mantle
  double precision, parameter :: EARTH_DEPTH_SECOND_DOUBLING_OPTIMAL = 1650000.d0
! in the outer core
  double precision, parameter :: EARTH_DEPTH_THIRD_DOUBLING_OPTIMAL  = 3860000.d0
! in the outer core
  double precision, parameter :: EARTH_DEPTH_FOURTH_DOUBLING_OPTIMAL = 5000000.d0

! to suppress element stretching for 3D moho surface
  logical,parameter :: SUPPRESS_MOHO_STRETCHING = .false.

! to suppress element stretching at 410/660 internal topography
! (i.e. creates mesh without 410/660 topography for Harvard model (s362ani,..))
  logical,parameter :: SUPPRESS_INTERNAL_TOPOGRAPHY = .false.

! by default, SAVE_MESH_FILES will use VTK formats; this will also output additional AVS_DX format files
  logical,parameter :: SAVE_MESHFILES_AVS_DX_FORMAT = .false.

!!-----------------------------------------------------------
!!
!! GPU optimization
!!
!!-----------------------------------------------------------
! added these parameters for the GPU version of the solver

! asynchronuous memcopy between CPU and GPU
  logical, parameter :: GPU_ASYNC_COPY = .true.

! mesh coloring
! add mesh coloring for the GPU + MPI implementation
! this is needed on NVIDIA hardware up to FERMI boards, included.
! Starting on KEPLER boards you can leave it off because on KEPLER hardware
! or higher atomic reduction operations have become as fast as resorting to mesh coloring.
  logical, parameter :: USE_MESH_COLORING_GPU = .false.
  integer, parameter :: MAX_NUMBER_OF_COLORS = 1000
! enhanced coloring:
! using Droux algorithm
! try several times with one more color before giving up
  logical, parameter :: USE_DROUX_OPTIMIZATION = .false.
  integer, parameter :: MAX_NB_TRIES_OF_DROUX_1993 = 15
! using balancing algorithm
! postprocess the colors to balance them if Droux (1993) algorithm is not used
  logical, parameter :: BALANCE_COLORS_SIMPLE_ALGO = .true.


!!-----------------------------------------------------------
!!
!! ADIOS Related values
!!
!!-----------------------------------------------------------

  integer, parameter :: ADIOS_BUFFER_SIZE_IN_MB = 400
  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD = "MPI"

! ADIOS transport methods for undo save_frame** data files (see ADIOS manual for details)
!! MPI (default)
  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD_UNDO_ATT = "MPI"   ! or check with: "POSIX", "MPI_AGGREGATE",..
  character(len=*), parameter :: ADIOS_METHOD_PARAMS_UNDO_ATT =  ''
! or
!! POSIX
!  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD_UNDO_ATT = "POSIX"
!  character(len=*), parameter :: ADIOS_METHOD_PARAMS_UNDO_ATT =  ''
! or
!! MPI_AGGREGATE
!  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD_UNDO_ATT = "MPI_AGGREGATE"
!  character(len=*), parameter :: ADIOS_METHOD_PARAMS_UNDO_ATT =  "num_aggregators=64,num_ost=672"
! or
!! MPI_LUSTRE
!  character(len=*), parameter :: ADIOS_TRANSPORT_METHOD_UNDO_ATT = "MPI_LUSTRE"
!  character(len=*), parameter :: ADIOS_METHOD_PARAMS_UNDO_ATT =  "stripe_count=16,stripe_size=4194304,block_size=4194304"


!!-----------------------------------------------------------
!!
!! ADIOS2 Related values
!!
!!-----------------------------------------------------------

! ADIOS2 engines
!! note on native engine types:
!!  - "MPI" is not supported yet by adios2 version (current 2.4.0), check out in future.
!!  - "BPfile" doesn't support file appending yet, which is needed at the moment.
!!  - "BP3" would allow for backward compatibility to ADIOS 1.x, but doesn't support file appending yet.
!!  - "BP4" is the new adios2 format with enhanced capabilities.
!! we will use "BP4" by default.
!!
!! BP4
!! format details: https://adios2.readthedocs.io/en/latest/engines/engines.html#bp4
!!
!! note: parameter SubStreams=64 for larger runs with NPROCS > 64 creates problems when reading scalar values (in appended mode),
!!       try to avoid it for now as default parameter.
!!       for undo_att, it seems to work however and can be used in ADIOS2_ENGINE_PARAMS_UNDO_ATT setting.
!!
!!       in future adios2 versions, re-evalute if parameters could be "SubStreams=64,MaxBufferSize=800Mb" for larger runs
  character(len=*), parameter :: ADIOS2_ENGINE_DEFAULT = "BP4"
  character(len=*), parameter :: ADIOS2_ENGINE_PARAMS_DEFAULT = "" ! add "MaxBufferSize=800Mb" for larger runs

  character(len=*), parameter :: ADIOS2_ENGINE_UNDO_ATT = "BP4"
  character(len=*), parameter :: ADIOS2_ENGINE_PARAMS_UNDO_ATT = "" ! add "SubStreams=64,MaxBufferSize=800Mb" for larger runs


!!-----------------------------------------------------------
!!
!! ASDF parameters
!!
!!-----------------------------------------------------------

! keeps track of everything
  logical, parameter :: ASDF_OUTPUT_PROVENANCE = .false.

! ASDF string lengths
  integer, parameter :: ASDF_MAX_STRING_LENGTH = 1024
  integer, parameter :: ASDF_MAX_QUAKEML_LENGTH = 8096
  integer, parameter :: ASDF_MAX_STATIONXML_LENGTH = 16182
  integer, parameter :: ASDF_MAX_PARFILE_LENGTH = 20000
  integer, parameter :: ASDF_MAX_CONSTANTS_LENGTH = 45000
  integer, parameter :: ASDF_MAX_TIME_STRING_LENGTH = 22


!!-----------------------------------------------------------
!!
!! adjoint kernel outputs
!!
!!-----------------------------------------------------------

! asynchronuous reading of adjoint sources
  logical, parameter :: IO_ASYNC_COPY = .true.

! regular kernel parameters
  character(len=*), parameter :: PATHNAME_KL_REG = 'DATA/kl_reg_grid.txt'
  integer, parameter :: NM_KL_REG_LAYER = 100
  integer, parameter :: NM_KL_REG_PTS = 300000
  real, parameter :: KL_REG_MIN_LON = 0.0
  real, parameter :: KL_REG_MAX_LON = 360.0
  real, parameter :: KL_REG_MIN_LAT = -90.0
  real, parameter :: KL_REG_MAX_LAT = +90.0

! by default, we turn kernel computations in the outer and inner core and for topology kernels off
!
! kernels for outer core
  logical, parameter :: SAVE_KERNELS_OC = .false.

! kernels for inner core
  logical, parameter :: SAVE_KERNELS_IC = .false.

! kernels for MOHO, 400, 670, CMB and ICB discontinuities
  logical, parameter :: SAVE_KERNELS_BOUNDARY = .false.

! Boundary Mesh -- save Moho, 400, 670 km discontinuity topology files (in
! the mesher) and use them for the computation of boundary kernel (in the solver)
  logical, parameter :: SAVE_BOUNDARY_MESH = SAVE_KERNELS_BOUNDARY


!!-----------------------------------------------------------
!!
!! new version compatibility
!!
!!-----------------------------------------------------------
! old version 5.1.5 uses full 3d attenuation arrays (set to .true.), CUSTOM_REAL for attenuation arrays (Qmu_store, tau_e_store)
! new version uses optional full 3d attenuation
  logical, parameter :: USE_OLD_VERSION_5_1_5_FORMAT = .false.
! for comparison against older versions: in version 7.0 tiso was constrained to below moho, new version allows also in crust
  logical, parameter :: USE_OLD_VERSION_7_0_0_FORMAT = .false.


!!-----------------------------------------------------------
!!
!! time stamp information
!!
!!-----------------------------------------------------------

! print date and time estimate of end of run in another country,
! in addition to local time.
! For instance: the code runs at Caltech in California but the person
! running the code is connected remotely from France, which has 9 hours more.
! The time difference with that remote location can be positive or negative
  logical, parameter :: ADD_TIME_ESTIMATE_ELSEWHERE = .false.
  integer, parameter :: HOURS_TIME_DIFFERENCE = +9
  integer, parameter :: MINUTES_TIME_DIFFERENCE = +0


!!-----------------------------------------------------------
!!
!! directory structure
!!
!!-----------------------------------------------------------

! paths for inputs and outputs files
  character(len=*), parameter :: OUTPUT_FILES_BASE = './OUTPUT_FILES/'


!!-----------------------------------------------------------
!!
!! movie outputs
!!
!!-----------------------------------------------------------

! runs external bash script after storing movie files
! (useful for postprocessing during simulation time)
  logical, parameter :: RUN_EXTERNAL_MOVIE_SCRIPT = .false.
  character(len=*),parameter :: MOVIE_SCRIPT_NAME = "./tar_movie_files.sh"


!!-----------------------------------------------------------
!!
!! for gravity integrals
!!
!!-----------------------------------------------------------
! (for using gravity integral computations, please make sure you compile with double-precision --enable-double-precision)

  logical, parameter :: GRAVITY_INTEGRALS = .false.

! reuse an existing observation surface created in another run and stored to disk,
! so that we are sure that they are exactly the same (for instance when comparing results for a reference ellipsoidal Earth
! and results for a 3D Earth with topography)
  logical, parameter :: REUSE_EXISTING_OBSERVATION_SURF = .false.

! only compute the center of mass of the 3D model (very cheap), or also compute the gravity integrals (expensive)
  logical, parameter :: ONLY_COMPUTE_CENTER_OF_MASS = .false.

! we may want to shift the reference frame to a pre-computed center of mass
  logical, parameter :: SHIFT_TO_THIS_CENTER_OF_MASS = .true.

! position of the center of mass of the Earth for this density model and mesh,
! but first convert it to meters and then make it non-dimensional
  double precision, parameter :: x_shift =   0.606220633681674d0 * 1000.d0 / EARTH_R !    km
  double precision, parameter :: y_shift =   0.433103991863316d0 * 1000.d0 / EARTH_R !    km
  double precision, parameter :: z_shift =   0.520078637496872d0 * 1000.d0 / EARTH_R !    km
!    distance to center =   0.908605697566306       km

! altitude of the observation points in km
  double precision, parameter :: altitude_of_observation_points = 255.d0

! elevation ratio of the points at which we observe
  double precision, parameter :: observation_elevation_ratio = (EARTH_R_KM + altitude_of_observation_points) / EARTH_R_KM

! compute the contribution of the crust only (otherwise by default compute the contribution of the whole Earth)
  logical, parameter :: COMPUTE_CRUST_CONTRIB_ONLY = .false.

! check for negative Jacobians in the calculation of integrals or not
! (can safely be done once to check that the mesh is OK at a given resolution, and then permanently
!  turned off in future runs because the mesh does not change)
  logical, parameter :: CHECK_FOR_NEGATIVE_JACOBIANS = .true.

! number of points in each horizontal direction of the observation grid of each cubed-sphere chunk
! at the altitude of the observation point
!! DK DK 4 is a fictitious value used to save memory when the GRAVITY_INTEGRALS option is off
  integer, parameter :: NX_OBSERVATION = 4 ! 500
  integer, parameter :: NY_OBSERVATION = NX_OBSERVATION

! the code will display sample output values at this particular point as a check
  integer, parameter :: ixr = max(int(NX_OBSERVATION / 3.0), 1)
  integer, parameter :: iyr = max(int(NY_OBSERVATION / 4.0), 1)
  integer, parameter :: ichunkr = 3

! how often (every how many spectral elements computed) we print a timestamp to monitor the behavior of the code
  integer, parameter :: NSPEC_DISPLAY_INTERVAL = 100


!!-----------------------------------------------------------
!!
!! debugging flags
!!
!!-----------------------------------------------------------

! flags to actually assemble with MPI or not
! and to actually match fluid and solid regions of the Earth or not
! should always be set to true except when debugging code
  logical, parameter :: ACTUALLY_ASSEMBLE_MPI_SLICES = .true.
  logical, parameter :: ACTUALLY_ASSEMBLE_MPI_CHUNKS = .true.
  logical, parameter :: ACTUALLY_COUPLE_FLUID_CMB = .true.
  logical, parameter :: ACTUALLY_COUPLE_FLUID_ICB = .true.

! flag to turn off the conversion of geographic to geocentric coordinates for
! the seismic source and the stations; i.e. assume a perfect sphere, which
! can be useful for benchmarks of a spherical Earth with fictitious sources and stations
  logical, parameter :: ASSUME_PERFECT_SPHERE = .false.

! flags to do benchmark runs to measure scaling of the code
! for a limited number of time steps only, setting the initial field to 1
! to make sure gradual underflow trapping does not slow down the code
  logical, parameter :: DO_BENCHMARK_RUN_ONLY = .false.
  integer, parameter :: NSTEP_FOR_BENCHMARK = 300
  logical, parameter :: SET_INITIAL_FIELD_TO_1_IN_BENCH = .true.

! for validation of undoing of attenuation versus an exact solution saved to disk
! should never be needed any more, it was developed and used for Figure 3 of
! D. Komatitsch, Z. Xie, E. Bozdag, E. Sales de Andrade, D. Peter, Q. Liu and J. Tromp,
! Anelastic sensitivity kernels with parsimonious storage for adjoint tomography and full waveform inversion,
! Geophysical Journal International, vol. 206(3), p. 1467-1478, doi: 10.1093/gji/ggw224 (2016).
!
! Compute an alpha sensitivity kernel directly by dumping the whole forward
! run to disk instead of rebuilding it backwards from the final time step in a second stage;
! used for debugging and validation purposes only, in order to have an exact reference for the sensitivity kernels.
! use wisely, this requires several terabytes (yes, tera) of disk space.
! In the future that option should become standard at some point though.
! For now, in order to save disk space, it is implemented for the alpha kernel only
! (because it only requires saving a scalar: the trace of epsilon) and for surface wave
! kernels only, this way we can save the upper mantle only and ignore the rest.
! In the future it will be easy to generalize this though.
  logical, parameter :: EXACT_UNDOING_TO_DISK = .false.
  ! ID of the huge file in which we dump all the time steps of the simulation
  integer, parameter :: IFILE_FOR_EXACT_UNDOING = 244
  ! Number of time steps between dumping the forward run
  integer, parameter :: NSTEP_FOR_EXACT_UNDOING = 500

!------------------------------------------------------
!----------- do not modify anything below -------------
!------------------------------------------------------

! on some processors (e.g. some Intel chips) it is necessary to suppress underflows
! by using a small initial field instead of zero
!! DK DK August 2018: on modern processors this does not happen any more,
!! DK DK August 2018: and thus no need to purposely lose accuracy to avoid underflows; thus turning it off by default
  logical, parameter :: FIX_UNDERFLOW_PROBLEM = .false. ! .true.

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: TWO_PI = 2.d0 * PI
  double precision, parameter :: PI_OVER_FOUR = PI / 4.d0
  double precision, parameter :: PI_OVER_TWO = PI / 2.0d0

! to convert angles from degrees to radians
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0
! to convert angles from radians to degrees
  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI

! 3-D simulation
  integer, parameter :: NDIM = 3

! dimension of the boundaries of the slices
  integer, parameter :: NDIM2D = 2

! number of nodes for 2D and 3D shape functions for hexahedra with 27 nodes
  integer, parameter :: NGNOD = 27, NGNOD2D = 9

! Deville routines optimized for NGLLX = NGLLY = NGLLZ = 5
  integer, parameter :: m1 = NGLLX, m2 = NGLLX * NGLLY
  integer, parameter :: NGLLCUBE = NGLLX * NGLLY * NGLLZ

! mid-points inside a GLL element
  integer, parameter :: MIDX = (NGLLX+1)/2
  integer, parameter :: MIDY = (NGLLY+1)/2
  integer, parameter :: MIDZ = (NGLLZ+1)/2

! gravitational constant in S.I. units i.e. in m3 kg-1 s-2, or equivalently in N.(m/kg)^2
!! DK DK April 2014: switched to the 2010 Committee on Data for Science and Technology (CODATA) recommended value
!! DK DK see e.g. http://www.physics.nist.gov/cgi-bin/cuu/Value?bg
!! DK DK and http://en.wikipedia.org/wiki/Gravitational_constant
!! double precision, parameter :: GRAV = 6.6723d-11
!! double precision, parameter :: GRAV = 6.67430d-11  ! newer suggestion by CODATA 2018
  double precision, parameter :: GRAV = 6.67384d-11  ! CODATA 2010

! standard gravity at the surface of the Earth
  double precision, parameter :: EARTH_STANDARD_GRAVITY = 9.80665d0 ! in m.s-2

! a few useful constants
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0,TWO = 2.d0,HALF = 0.5d0
  double precision, parameter :: ONE_HALF = HALF, ONE_FOURTH = 0.25d0, ONE_EIGHTH = 0.125d0, ONE_SIXTEENTH = 0.0625d0

  real(kind=CUSTOM_REAL), parameter :: &
    ONE_THIRD   = 1._CUSTOM_REAL/3._CUSTOM_REAL, &
    TWO_THIRDS  = 2._CUSTOM_REAL/3._CUSTOM_REAL, &
    FOUR_THIRDS = 4._CUSTOM_REAL/3._CUSTOM_REAL

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30,TINYVAL = 1.d-9

! very large real value declared independently of the machine
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_SNGL = 1.e+30_CUSTOM_REAL

! very large integer value
  integer, parameter :: HUGEINT = 100000000

! normalized radius of free surface
  double precision, parameter :: R_UNIT_SPHERE = ONE

! fixed thickness of 3 km for PREM oceans
  double precision, parameter :: THICKNESS_OCEANS_PREM = 3000.d0 / EARTH_R

! shortest radius at which crust is implemented (80 km depth)
! to be constistent with the D80 discontinuity, we impose the crust only above it
  double precision, parameter :: EARTH_R_DEEPEST_CRUST = (EARTH_R - 80000.d0) / EARTH_R

! definition of an Eotvos compared to S.I. units.
! The unit of gravity gradient is the Eotvos, which is equivalent to 1e-9 s-2 (or 1e-4 mGal/m).
! A person walking at a distance of 2 meters provides a gravity gradient signal of approximately one Eotvos.
! Mountains can create signals of several hundred Eotvos.
  double precision, parameter :: SI_UNITS_TO_EOTVOS = 1.d+9

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

! define flag for planet
! strictly speaking, the moon is not a planet but a natural satellite. well, let's stay with IPLANET_** for now.
  integer, parameter :: IPLANET_EARTH = 1
  integer, parameter :: IPLANET_MARS = 2
  integer, parameter :: IPLANET_MOON = 3

! define flag for regions of the global Earth mesh
  integer, parameter :: IREGION_CRUST_MANTLE = 1
  integer, parameter :: IREGION_OUTER_CORE = 2
  integer, parameter :: IREGION_INNER_CORE = 3

! define flag for elements
  integer, parameter :: IFLAG_CRUST = 1

  integer, parameter :: IFLAG_80_MOHO = 2
  integer, parameter :: IFLAG_220_80 = 3
  integer, parameter :: IFLAG_670_220 = 4
  integer, parameter :: IFLAG_MANTLE_NORMAL = 5

  integer, parameter :: IFLAG_OUTER_CORE_NORMAL = 6

  integer, parameter :: IFLAG_INNER_CORE_NORMAL = 7
  integer, parameter :: IFLAG_MIDDLE_CENTRAL_CUBE = 8
  integer, parameter :: IFLAG_BOTTOM_CENTRAL_CUBE = 9
  integer, parameter :: IFLAG_TOP_CENTRAL_CUBE = 10
  integer, parameter :: IFLAG_IN_FICTITIOUS_CUBE = 11

  integer, parameter :: NSPEC2D_XI_SUPERBRICK = 8
  integer, parameter :: NSPEC2D_ETA_SUPERBRICK = 8
  integer, parameter :: NSPEC2D_XI_SUPERBRICK_1L = 6
  integer, parameter :: NSPEC2D_ETA_SUPERBRICK_1L = 6

! dummy flag used for mesh display purposes only
  integer, parameter :: IFLAG_DUMMY = 100

! max number of layers that are used in the radial direction to build the full mesh
  integer, parameter :: MAX_NUMBER_OF_MESH_LAYERS = 15

! define number of spectral elements and points in basic symmetric mesh doubling superbrick
  integer, parameter :: NSPEC_DOUBLING_SUPERBRICK = 32
  integer, parameter :: NGLOB_DOUBLING_SUPERBRICK = 67
  integer, parameter :: NSPEC_SUPERBRICK_1L = 28
  integer, parameter :: NGLOB_SUPERBRICK_1L = 58
  integer, parameter :: NGNOD_EIGHT_CORNERS = 8

! define flag for reference 1D Earth model
  integer, parameter :: REFERENCE_MODEL_PREM   = 1
  integer, parameter :: REFERENCE_MODEL_IASP91 = 2
  integer, parameter :: REFERENCE_MODEL_1066A  = 3
  integer, parameter :: REFERENCE_MODEL_AK135F_NO_MUD = 4
  integer, parameter :: REFERENCE_MODEL_1DREF = 5
  integer, parameter :: REFERENCE_MODEL_JP1D  = 6
  integer, parameter :: REFERENCE_MODEL_SEA1D = 7
  integer, parameter :: REFERENCE_MODEL_SOHL = 8        ! Mars
  integer, parameter :: REFERENCE_MODEL_SOHL_B = 9      ! Mars
  integer, parameter :: REFERENCE_MODEL_CASE65TAY = 10  ! Mars
  integer, parameter :: REFERENCE_MODEL_VPREMOON = 11   ! Moon
  integer, parameter :: REFERENCE_MODEL_MOON_MEENA = 12 ! Moon

! crustal model constants
  integer, parameter :: ICRUST_CRUST1    = 1    ! Crust1.0
  integer, parameter :: ICRUST_CRUST2    = 2    ! Crust2.0
  integer, parameter :: ICRUST_CRUSTMAPS = 3    ! Crustmaps
  integer, parameter :: ICRUST_EPCRUST   = 4    ! EPcrust
  integer, parameter :: ICRUST_CRUST_SH  = 5    ! Spherical-Harmonics Crust
  integer, parameter :: ICRUST_EUCRUST   = 6    ! EUcrust07

! define flag for 3D Earth model
  integer, parameter :: THREE_D_MODEL_S20RTS        = 101
  integer, parameter :: THREE_D_MODEL_S362ANI       = 102
  integer, parameter :: THREE_D_MODEL_S362WMANI     = 103
  integer, parameter :: THREE_D_MODEL_S362ANI_PREM  = 104
  integer, parameter :: THREE_D_MODEL_S29EA         = 105
  integer, parameter :: THREE_D_MODEL_SEA99_JP3D    = 106
  integer, parameter :: THREE_D_MODEL_SEA99         = 107
  integer, parameter :: THREE_D_MODEL_JP3D          = 108
  integer, parameter :: THREE_D_MODEL_PPM           = 109            ! format for point profile models
  integer, parameter :: THREE_D_MODEL_GLL           = 110            ! format for iterations with GLL mesh
  integer, parameter :: THREE_D_MODEL_S40RTS        = 111
  integer, parameter :: THREE_D_MODEL_GAPP2         = 112
  integer, parameter :: THREE_D_MODEL_MANTLE_SH     = 113
  integer, parameter :: THREE_D_MODEL_SGLOBE        = 114
  integer, parameter :: THREE_D_MODEL_SGLOBE_ISO    = 115
  integer, parameter :: THREE_D_MODEL_ANISO_MANTLE  = 116
  ! inner core model
  integer, parameter :: THREE_D_MODEL_INNER_CORE_ISHII = 201

!! attenuation
! number of standard linear solids for attenuation
  integer, parameter :: N_SLS = 3

! computation of standard linear solids in meshfem3D
! ATTENUATION_COMP_RESOLUTION: Number of Digits after decimal
! ATTENUATION_COMP_MAXIMUM:    Maximum Q Value
  integer, parameter :: ATTENUATION_COMP_RESOLUTION = 1
  integer, parameter :: ATTENUATION_COMP_MAXIMUM    = 9000

! for determination of the attenuation period range
! if this is set to .true. then the hardcoded values will be used
! otherwise they are computed automatically from the Number of elements
! This *may* be a useful parameter for Benchmarking against older versions
  logical, parameter :: ATTENUATION_RANGE_PREDEFINED = .false.

! flag for the four edges of each slice and for the bottom edge
  integer, parameter :: XI_MIN  = 1
  integer, parameter :: XI_MAX  = 2
  integer, parameter :: ETA_MIN = 3
  integer, parameter :: ETA_MAX = 4
  integer, parameter :: BOTTOM  = 5

! flags to select the right corner in each slice
  integer, parameter :: ILOWERLOWER = 1
  integer, parameter :: ILOWERUPPER = 2
  integer, parameter :: IUPPERLOWER = 3
  integer, parameter :: IUPPERUPPER = 4

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

! number of layers in PREM
  integer, parameter :: NR_DENSITY = 640

! smallest real number on many machines =  1.1754944E-38
! largest real number on many machines =  3.4028235E+38
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

! number of lines per source in CMTSOLUTION file
  integer, parameter :: NLINES_PER_CMTSOLUTION_SOURCE = 13
  integer, parameter :: NLINES_PER_FORCESOLUTION_SOURCE = 11

! number of iterations to solve the system for xi and eta
! setting it to 5 instead of 4 ensures that the result obtained is not compiler dependent
! (when using 4 only some small discrepancies were observed)
  integer, parameter :: NUM_ITER = 5

! number of hours per day for rotation rate of the Earth
  double precision, parameter :: EARTH_HOURS_PER_DAY = 24.d0
  double precision, parameter :: EARTH_SECONDS_PER_HOUR = 3600.d0

! for lookup table for gravity every 100 m in radial direction of Earth model
  integer, parameter :: NRAD_GRAVITY = 70000


!!-----------------------------------------------------------
!!
!! GLL model constants
!!
!!-----------------------------------------------------------

! parameters for GLL model (used for iterative model inversions)
  character(len=*), parameter :: PATHNAME_GLL_modeldir = 'DATA/GLL/'

! to create a reference model based on 1D_REF but with 3D crust and 410/660 topography
  logical,parameter :: USE_1D_REFERENCE = .false.

!!-- uncomment for using PREM as reference model (used in CEM inversion)
!  integer, parameter :: GLL_REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
!  integer, parameter :: GLL_REFERENCE_MODEL = REFERENCE_MODEL_PREM

!!-- uncomment for using S362ANI as reference model
  integer, parameter :: GLL_REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
  integer, parameter :: GLL_REFERENCE_MODEL = THREE_D_MODEL_S362ANI

!!-- uncomment for using S29EA as reference model
!  integer, parameter :: GLL_REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
!  integer, parameter :: GLL_REFERENCE_MODEL = THREE_D_MODEL_S29EA

!!-----------------------------------------------------------
!!
!! crustal stretching
!!
!!-----------------------------------------------------------

! for the stretching of crustal elements in the case of 3D models
! (values are chosen for 3D models to have RMOHO_FICTICIOUS at 35 km
!  and RMIDDLE_CRUST to become 15 km with stretching function stretch_tab)
  double precision, parameter :: EARTH_MAX_RATIO_CRUST_STRETCHING = 0.75d0
  double precision, parameter :: EARTH_RMOHO_STRETCH_ADJUSTMENT = 5000.d0 ! moho up to 35km
  double precision, parameter :: EARTH_R80_STRETCH_ADJUSTMENT = -40000.d0 ! r80 down to 120km

! adapted regional moho stretching
! 1 chunk simulations, 3-layer crust
  logical, parameter :: EARTH_REGIONAL_MOHO_MESH = .false.
  logical, parameter :: EARTH_REGIONAL_MOHO_MESH_EUROPE = .false. ! used only for fixing time step
  logical, parameter :: EARTH_REGIONAL_MOHO_MESH_ASIA = .false.   ! used only for fixing time step
  logical, parameter :: EARTH_HONOR_DEEP_MOHO = .false.
!!-- uncomment for e.g. Europe case, where deep moho is rare
!  double precision, parameter :: EARTH_RMOHO_STRETCH_ADJUSTMENT = -15000.d0  ! moho mesh boundary down to 55km
!!-- uncomment for deep moho cases, e.g. Asia case (Himalayan moho)
!  double precision, parameter :: EARTH_RMOHO_STRETCH_ADJUSTMENT = -20000.d0  ! moho mesh boundary down to 60km

!!-----------------------------------------------------------
!!
!! mesh tweaking
!!
!!-----------------------------------------------------------

! to suppress the crustal layers
! (replaced by an extension of the mantle: EARTH_R is not modified, but no more crustal doubling)
  logical, parameter :: SUPPRESS_CRUSTAL_MESH = .false.

! to inflate the central cube (set to 0.d0 for a non-inflated cube)
  double precision, parameter :: CENTRAL_CUBE_INFLATE_FACTOR = 0.41d0

! to add a fourth doubling at the bottom of the outer core
  logical, parameter :: ADD_4TH_DOUBLING = .false.

! parameters to cut the doubling brick

! this to cut the superbrick: 3 possibilities, 4 cases max / possibility
! three possibilities: (cut in xi and eta) or (cut in xi) or (cut in eta)
! case 1: (ximin and etamin) or ximin or etamin
! case 2: (ximin and etamax) or ximax or etamax
! case 3: ximax and etamin
! case 4: ximax and etamax
  integer, parameter :: NB_CUT_CASE = 4

! corner 1: ximin and etamin
! corner 2: ximax and etamin
! corner 3: ximax and etamax
! corner 4: ximin and etamax
  integer, parameter :: NB_SQUARE_CORNERS = 4

! two possibilities: xi or eta
! face 1: ximin or etamin
! face 2: ximax or etamax
  integer, parameter :: NB_SQUARE_EDGES_ONEDIR = 2

! this for the geometry of the basic doubling brick
  integer, parameter :: NSPEC_DOUBLING_BASICBRICK = 8
  integer, parameter :: NGLOB_DOUBLING_BASICBRICK = 27

!!-----------------------------------------------------------
!!
!! for LDDRK high-order time scheme
!!
!!-----------------------------------------------------------

! Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm
!
! reference:
! J. Berland, C. Bogey, and C. Bailly.
! Low-dissipation and low-dispersion fourth-order Runge-Kutta algorithm.
! Computers and Fluids, 35:1459-1463, 2006
!
! see: http://www.sciencedirect.com/science/article/pii/S0045793005000575?np=y

! number of stages
  integer, parameter :: NSTAGE = 6

! coefficients from Table 1, Berland et al. (2006)
  real(kind=CUSTOM_REAL), dimension(NSTAGE), parameter :: ALPHA_LDDRK = &
    (/0.0_CUSTOM_REAL,-0.737101392796_CUSTOM_REAL, -1.634740794341_CUSTOM_REAL, &
      -0.744739003780_CUSTOM_REAL,-1.469897351522_CUSTOM_REAL,-2.813971388035_CUSTOM_REAL/)

  real(kind=CUSTOM_REAL), dimension(NSTAGE), parameter :: BETA_LDDRK = &
    (/0.032918605146_CUSTOM_REAL,0.823256998200_CUSTOM_REAL,0.381530948900_CUSTOM_REAL, &
      0.200092213184_CUSTOM_REAL,1.718581042715_CUSTOM_REAL,0.27_CUSTOM_REAL/)

  real(kind=CUSTOM_REAL), dimension(NSTAGE), parameter :: C_LDDRK = &
    (/0.0_CUSTOM_REAL,0.032918605146_CUSTOM_REAL,0.249351723343_CUSTOM_REAL, &
      0.466911705055_CUSTOM_REAL,0.582030414044_CUSTOM_REAL,0.847252983783_CUSTOM_REAL/)


!!-----------------------------------------------------------
!!
!! Mars
!!
!!-----------------------------------------------------------
! Model MARS
  double precision, parameter :: MARS_R = 3390000.d0
  double precision, parameter :: R_MARS = MARS_R

! radius of the MARS in km
  double precision, parameter :: MARS_R_KM = MARS_R / 1000.d0
  double precision, parameter :: R_MARS_KM = MARS_R_KM

! average density in the full MARS to normalize equation
  double precision, parameter :: MARS_RHOAV = 3393.0d0

! Topography
!--- 4-minute model
! size of topography and bathymetry file
  integer, parameter :: MARS_NX_BATHY_4 = 5400,MARS_NY_BATHY_4 = 2700
! resolution of topography file in minutes
  double precision, parameter :: MARS_RESOLUTION_TOPO_FILE_4 = 4.0
! pathname of the topography file
  character (len=*), parameter :: MARS_PATHNAME_TOPO_FILE_4 = 'DATA/topo_bathy/topo_bathy_marstopo4_smoothed_window_3.dat.bin'
!--- Default
! Topography defaults to 4-minute
  integer, parameter :: MARS_NX_BATHY = MARS_NX_BATHY_4
  integer, parameter :: MARS_NY_BATHY = MARS_NY_BATHY_4
  double precision, parameter :: MARS_RESOLUTION_TOPO_FILE = MARS_RESOLUTION_TOPO_FILE_4
  character (len=*), parameter :: MARS_PATHNAME_TOPO_FILE = MARS_PATHNAME_TOPO_FILE_4

! Topography value min/max range
  integer, parameter :: MARS_TOPO_MINIMUM = - 8000 ! (max depth in m, Mars)
  integer, parameter :: MARS_TOPO_MAXIMUM = + 23200 ! (height in m, Olympus Mons)

! For the reference ellipsoid to convert geographic latitudes to geocentric:
! Mars flattening (https://en.wikipedia.org/wiki/Mars): 0.00589 +/- 0.00015 -> 1/f with f ~ 169.77
! Smith et al. 1999, The global topography of Mars and implications for surface evolution. Science 284(5419), 1495-1503
! uses f = 169.8
  double precision, parameter :: MARS_FLATTENING_F = 1.d0 / 169.8d0
  double precision, parameter :: MARS_ONE_MINUS_F_SQUARED = (1.d0 - MARS_FLATTENING_F)**2

! doubling layers
  double precision, parameter :: MARS_DEPTH_SECOND_DOUBLING_OPTIMAL = 1200000.d0
  double precision, parameter :: MARS_DEPTH_THIRD_DOUBLING_OPTIMAL  = 2090000.d0
  double precision, parameter :: MARS_DEPTH_FOURTH_DOUBLING_OPTIMAL = 2690000.d0

! standard gravity at the surface
  double precision, parameter :: MARS_STANDARD_GRAVITY = 3.71d0 ! in m.s-2

! shortest radius at which crust is implemented (150 km depth)
! we impose the crust only above it
  double precision, parameter :: MARS_R_DEEPEST_CRUST = (MARS_R - 150000.d0) / MARS_R

! number of hours per day for rotation rate of the planet Mars (24:39:35)
  double precision, parameter :: MARS_HOURS_PER_DAY = 24.658d0
  double precision, parameter :: MARS_SECONDS_PER_HOUR = 3600.d0

! for the stretching of crustal elements in the case of 3D models
  double precision, parameter :: MARS_MAX_RATIO_CRUST_STRETCHING = 0.9d0 ! choose ratio < 1 for stretching effect
  double precision, parameter :: MARS_RMOHO_STRETCH_ADJUSTMENT = 0.d0 ! moho at 110 km
  double precision, parameter :: MARS_R80_STRETCH_ADJUSTMENT =   0.d0 ! r80  at 334.5 km

!!double precision, parameter :: MARS_MAX_RATIO_CRUST_STRETCHING = 0.7d0 ! choose ratio < 1 for stretching effect
!!double precision, parameter :: MARS_RMOHO_STRETCH_ADJUSTMENT = -10000.d0 ! moho down to 120 km
!!double precision, parameter :: MARS_R80_STRETCH_ADJUSTMENT =   -20000.d0 ! r80 down to 120km

! adapted regional moho stretching (use only for special areas to optimize a local mesh)
! 1~6 chunk simulation, 5-layer crust
  logical, parameter :: MARS_REGIONAL_MOHO_MESH = .false.
  logical, parameter :: MARS_HONOR_DEEP_MOHO = .false.


!!-----------------------------------------------------------
!!
!! Moon
!!
!!-----------------------------------------------------------
! Model Moon
!
! mostly based on VPREMOON model:
! Garcia, R.F., J. Gagnepain-Beyneix, S. Chevrot and Ph. Lognonne, 2011.
! Very preliminary reference Moon model,
! Physics of the Earth and Planetary Interiors, 188, 96-113.
!
! NASA fact infos: https://nssdc.gsfc.nasa.gov/planetary/factsheet/moonfact.html
! some of these NASA values are different with respect to the seismologically preferred ones below.
!
! average Moon radius
  double precision, parameter :: MOON_R = 1737100.d0
  double precision, parameter :: R_MOON = MOON_R

! radius of the Moon in km
  double precision, parameter :: MOON_R_KM = MOON_R / 1000.d0
  double precision, parameter :: R_MOON_KM = MOON_R_KM

! average density in the full MOON to normalize equation
  double precision, parameter :: MOON_RHOAV = 3344.0d0    ! from NASA fact sheet

! Topography
!--- 4-minute model
  integer, parameter :: MOON_NX_BATHY_4 = 5400,MOON_NY_BATHY_4 = 2700
! resolution of topography file in minutes
  double precision, parameter :: MOON_RESOLUTION_TOPO_FILE_4 = 4.0
! pathname of the topography file
  character (len=*), parameter :: MOON_PATHNAME_TOPO_FILE_4 = 'DATA/topo_bathy/topo_bathy_moon4_smoothed_window_0.dat.bin'
!--- 2-minute
  integer, parameter :: MOON_NX_BATHY_2 = 10800,MOON_NY_BATHY_2 = 5400
! resolution of topography file in minutes
  double precision, parameter :: MOON_RESOLUTION_TOPO_FILE_2 = 2.0
! pathname of the topography file
  character (len=*), parameter :: MOON_PATHNAME_TOPO_FILE_2 = 'DATA/topo_bathy/topo_bathy_moon2_smoothed_window_0.dat.bin'
!--- 1-minute
  integer, parameter :: MOON_NX_BATHY_1 = 21600,MOON_NY_BATHY_1 = 10800
! resolution of topography file in minutes
  double precision, parameter :: MOON_RESOLUTION_TOPO_FILE_1 = 1.0
! pathname of the topography file
  character (len=*), parameter :: MOON_PATHNAME_TOPO_FILE_1 = 'DATA/topo_bathy/topo_bathy_moon1_smoothed_window_0.dat.bin'
!--- custom resolution < 1-minute model
! size of topography and bathymetry file
  integer, parameter :: MOON_NX_BATHY_C = 23040,MOON_NY_BATHY_C = 11520
! resolution of topography file in minutes
  double precision, parameter :: MOON_RESOLUTION_TOPO_FILE_C = 0.94d0
! pathname of the topography file
  character (len=*), parameter :: MOON_PATHNAME_TOPO_FILE_C = 'DATA/topo_bathy/interface.bin'
!--- Default
! Topography defaults to 4-minute
  integer, parameter :: MOON_NX_BATHY = MOON_NX_BATHY_4
  integer, parameter :: MOON_NY_BATHY = MOON_NY_BATHY_4
  double precision, parameter :: MOON_RESOLUTION_TOPO_FILE = MOON_RESOLUTION_TOPO_FILE_4
  character (len=*), parameter :: MOON_PATHNAME_TOPO_FILE = MOON_PATHNAME_TOPO_FILE_4

! Topography value min/max range
  integer, parameter :: MOON_TOPO_MINIMUM = -9130 ! (max depth in m)
  integer, parameter :: MOON_TOPO_MAXIMUM = 10780 ! (height in m)

! For the reference ellipsoid to convert geographic latitudes to geocentric:
! Moon flattening:
! NASA fact sheet has flattening = 0.0012 = 1/833.3
!
! Note: Steinberger et al. (2015, PEPI 245, 26-39) argue that due to non-equilibrum flattening of the Moon
!       there is no reference ellipsoid to refer to. Thus, they prefer a sphere with a "frozen"-in shape.
!
!       We could in principle do the same, as the Moon topography is also given in absolute values (radius).
!       However, for now we take topography as the +/- elevation value with respect to a reference surface radius.
!
! For flattening, here preferred: 1/901 = 0.00110987791
  double precision, parameter :: MOON_FLATTENING_F = 1.d0 / 901.0d0
  double precision, parameter :: MOON_ONE_MINUS_F_SQUARED = (1.d0 - MOON_FLATTENING_F)**2

! doubling layers
  double precision, parameter :: MOON_DEPTH_SECOND_DOUBLING_OPTIMAL = 771000.d0    ! between RTOPDDOUBLEPRIME and R771
  double precision, parameter :: MOON_DEPTH_THIRD_DOUBLING_OPTIMAL  = 1380000.d0    ! inside outer core (closer to top)
  double precision, parameter :: MOON_DEPTH_FOURTH_DOUBLING_OPTIMAL = 1420000.d0    ! inside outer core (closer to bottom)

! standard gravity at the surface
  double precision, parameter :: MOON_STANDARD_GRAVITY = 1.62d0 ! in m.s-2

! shortest radius at which crust is implemented (39 km depth)
! we impose the crust only above it
  double precision, parameter :: MOON_R_DEEPEST_CRUST = (MOON_R - 39000.d0) / MOON_R

! number of hours per day for rotation rate of the Moon
  double precision, parameter :: MOON_HOURS_PER_DAY = 27.322d0
  double precision, parameter :: MOON_SECONDS_PER_HOUR = 3600.d0

! for the stretching of crustal elements in the case of 3D models
  double precision, parameter :: MOON_MAX_RATIO_CRUST_STRETCHING = 0.7d0 ! 0.75
  double precision, parameter :: MOON_RMOHO_STRETCH_ADJUSTMENT = -10000.d0 ! moho down
  double precision, parameter :: MOON_R80_STRETCH_ADJUSTMENT =   -20000.d0 ! r80 down

! adapted regional moho stretching (use only for special areas to optimize a local mesh)
  logical, parameter :: MOON_REGIONAL_MOHO_MESH = .false.
  logical, parameter :: MOON_HONOR_DEEP_MOHO = .false.
