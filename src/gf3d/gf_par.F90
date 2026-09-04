!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

!----
!---- Green function extraction library: derived types, error codes and
!---- module-level state.
!----
!---- This module is deliberately free of `use hdf5`, `#ifdef USE_HDF5`
!---- and `use specfem_par`, so that it builds from a plain ./configure
!---- (see tests/gf3d/ and gf3d_kernels in src/gf3d/rules.mk).
!----
!---- Precision policy, stated once and held to throughout src/gf3d/:
!----   read float32 from the database, compute in double, return double.
!---- The database is written with CUSTOM_REAL = 4; every quantity kept
!---- in t_gfdb is double precision.
!----

  module gf_par

  use constants, only: MAX_STRING_LEN,MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME,NDIM

  implicit none

  ! version of the library and of the on-disk database layout it reads
  character(len=*), parameter :: GF3D_VERSION = '0.1.0'

  ! number of force components (N,E,Z) and displacement components (x,y,z)
  ! stored per element; mirrors GF_NCOMP_FORCE / GF_NCOMP_DISP in the writer
  integer, parameter :: GF_NCOMP = 3

  ! Morton codes are formatted as 16 hexadecimal digits by the writer,
  ! see green_function_morton.F90 -> write(gf_morton_hex(i),'(Z16.16)')
  integer, parameter :: GF_MORTON_HEXLEN = 16

  ! combined "NET.STA" identifier
  integer, parameter :: GF_STATION_ID_LEN = MAX_LENGTH_NETWORK_NAME + 1 + MAX_LENGTH_STATION_NAME

  !-----------------------------------------------------------------
  ! source location
  !-----------------------------------------------------------------

  ! Containment tolerance on the local coordinates.
  !
  ! This replicates the solver rather than asking for a mathematically
  ! containing element, and the value is not arbitrary:
  ! find_local_coordinates() (locate_point.f90:513-519) clamps each of
  ! xi/eta/gamma to +-1.10, and locate_point.f90:377 then treats anything
  ! above 1.099 as "this point probably belongs to a neighbour" and starts
  ! an adjacency walk. So 1.099 is exactly the solver's own test for "the
  ! element I have is the element I want", and accepting <= 1.1 instead
  ! would accept an iterate that merely hit the clamp, i.e. one that failed.
  !
  ! The shipped global example needs this: output_solver.txt reports
  ! gamma = 1.05304706 for its CMTSOLUTION, so the forward run itself used
  ! an element that does not strictly contain its source. Matching that
  ! choice matters because the forward source was distributed onto that
  ! element's GLL points with that basis, and reproducing the forward
  ! seismogram is what this library is judged on.
  double precision, parameter :: GF_XI_TOL = 1.099d0

  ! number of candidate elements tried, in centroid-distance order
  integer, parameter :: GF_NCAND = 10

  ! Tolerance of the 27-anchor consistency guard, non-dimensional.
  !
  ! With USE_GLL = .false. (setup/constants.h:USE_GLL) the mesher applies
  ! topography and ellipticity to the 27 anchors and re-interpolates the GLL
  ! points with the tri-quadratic shape functions, so the stored xyz(3,5,5,5)
  ! *is* a tri-quadratic sampled at GLL points and the anchors must reproduce
  ! it exactly -- in exact arithmetic.
  !
  ! In practice the floor is the storage precision, not the arithmetic. The
  ! solver holds xstore_crust_mantle as real(CUSTOM_REAL) and the writer
  ! copies it straight through (green_function_metadata.F90:82,119-121), so
  ! with CUSTOM_REAL = 4 both the anchors and the values they must reproduce
  ! carry a float32 rounding. Measured over all 82 elements x 125 points of
  ! the shipped global example the worst residual is 6.2e-8, which is
  ! float32 epsilon times the O(0.5) magnitude of the coordinates.
  !
  ! 1e-6 is therefore a comfortable margin over the float32 floor while still
  ! failing loudly on the thing this guard exists to catch: a USE_GLL = .true.
  ! database, where the map is not tri-quadratic at all and the residual is
  ! 1e-3 to 1e-2.
  double precision, parameter :: GF_ANCHOR_TOL = 1.d-6

  !-----------------------------------------------------------------
  ! error codes
  !
  ! Every public entry point of the library returns one of these in an
  ! `ierr` argument. Nothing in src/gf3d/ calls stop or exit_MPI: the
  ! library is meant to be loaded into a Python interpreter (Stage 9),
  ! where a stop would kill the caller.
  !-----------------------------------------------------------------

  integer, parameter :: GF_OK              =  0
  integer, parameter :: GF_ERR_NO_HDF5     =  1   ! built without --with-hdf5
  integer, parameter :: GF_ERR_NO_PATH     =  2   ! database directory missing
  integer, parameter :: GF_ERR_NO_FILE     =  3   ! an expected file is missing
  integer, parameter :: GF_ERR_HDF5        =  4   ! an HDF5 call failed
  integer, parameter :: GF_ERR_IO          =  5   ! a Fortran read/open failed
  integer, parameter :: GF_ERR_FORMAT      =  6   ! file present but malformed
  integer, parameter :: GF_ERR_MISMATCH    =  7   ! database incompatible with this build
  integer, parameter :: GF_ERR_INCOMPLETE  =  8   ! element/station data not fully written
  integer, parameter :: GF_ERR_ALLOC       =  9   ! allocation failed
  integer, parameter :: GF_ERR_ARG         = 10   ! invalid argument from the caller
  integer, parameter :: GF_ERR_NO_ELEMENT  = 11   ! no database element contains the point
  integer, parameter :: GF_ERR_GEOMETRY    = 12   ! degenerate element geometry

  !-----------------------------------------------------------------
  ! per-station metadata, read from {GFDB}/stations/{net}.{sta}.h5
  !-----------------------------------------------------------------

  type :: t_gf_station
    character(len=MAX_LENGTH_NETWORK_NAME) :: network = ''
    character(len=MAX_LENGTH_STATION_NAME) :: station = ''
    character(len=GF_STATION_ID_LEN) :: id = ''          ! 'NET.STA'

    ! reciprocal source location (the station itself)
    double precision :: latitude  = 0.d0
    double precision :: longitude = 0.d0
    double precision :: depth     = 0.d0                 ! in m, below the surface

    ! source time function used for the reciprocal runs
    double precision :: hdur                = 0.d0
    double precision :: f_cutoff            = 0.d0
    double precision :: factor_force_source = 0.d0
    double precision :: time_shift          = 0.d0

    ! stf(nstep), on the *unsubsampled* solver time axis
    double precision, dimension(:), allocatable :: stf
  end type t_gf_station

  !-----------------------------------------------------------------
  ! the database handle
  !
  ! Contract: open once, extract many. gf_open() reads all of the
  ! metadata (a few hundred kB) but none of the bulk arrays; the 58 MB
  ! ibathy_topo grid is loaded on demand by gf_load_topo().
  !-----------------------------------------------------------------

  type :: t_gfdb
    logical :: is_open = .false.
    character(len=MAX_STRING_LEN) :: path = ''

    !--- simulation parameters, from mesh_info.h5 ---
    double precision :: dt          = 0.d0    ! solver time step, s
    double precision :: t0          = 0.d0    ! time of the first sample, s before origin
    double precision :: scale_displ = 0.d0    ! non-dimensionalisation of the stored displacement
    double precision :: R_PLANET    = 0.d0    ! m
    double precision :: RHOAV       = 0.d0    ! kg/m^3, see note in gf_database.F90

    integer :: nstep          = 0             ! solver time steps
    integer :: nt_subsampled  = 0             ! stored time samples = nstep / subsample_step
    integer :: subsample_step = 0
    integer :: buffer_size    = 0
    integer :: neighbor_shells = 0

    integer :: ngllx = 0, nglly = 0, ngllz = 0

    logical :: topography  = .false.
    logical :: ellipticity = .false.
    logical :: rotation    = .false.
    logical :: attenuation = .false.
    logical :: gravity     = .false.

    !--- topography grid (loaded lazily by gf_load_topo) ---
    integer :: NX_BATHY = 0, NY_BATHY = 0
    double precision :: RESOLUTION_TOPO_FILE = 0.d0
    logical :: topo_loaded = .false.
    integer, dimension(:,:), allocatable :: ibathy_topo     ! (NX_BATHY,NY_BATHY)

    !--- ellipticity splines ---
    integer :: nspl = 0
    double precision, dimension(:), allocatable :: rspl
    double precision, dimension(:), allocatable :: ellipicity_spline
    double precision, dimension(:), allocatable :: ellipicity_spline2

    !--- element index ---
    integer :: nelem = 0
    integer(kind=8), dimension(:), allocatable :: morton              ! (nelem)
    character(len=GF_MORTON_HEXLEN), dimension(:), allocatable :: morton_hex  ! (nelem)
    double precision, dimension(:,:), allocatable :: centroid         ! (3,nelem), non-dimensional

    ! where the element index came from: 'centroids.bin' or 'manifest.csv'
    character(len=16) :: index_source = ''

    !--- stations ---
    integer :: nstations = 0
    type(t_gf_station), dimension(:), allocatable :: stations
  end type t_gfdb

  !-----------------------------------------------------------------
  ! a located source
  !
  ! Everything the extraction needs about where a source sits, produced
  ! once by gf_locate_source() and consumed by every later stage. The
  ! inverse Jacobian is carried here rather than recomputed because the
  ! strain (Stage 4) needs exactly the values the locator already had.
  !-----------------------------------------------------------------

  type :: t_gf_location
    !--- which element ---
    integer :: ielem = 0                                     ! 1..db%nelem
    character(len=GF_MORTON_HEXLEN) :: morton_hex = ''

    !--- where in it ---
    double precision :: xi = 0.d0, eta = 0.d0, gamma = 0.d0

    !--- the mapped position, non-dimensional Cartesian ---
    double precision, dimension(NDIM) :: xyz = 0.d0

    !--- the target position, before the element map ---
    double precision, dimension(NDIM) :: xyz_target = 0.d0

    !--- inverse Jacobian d(xi,eta,gamma)/d(x,y,z) at the source ---
    ! rows are xi/eta/gamma, columns x/y/z, i.e. jinv(1,2) is xiy
    double precision, dimension(NDIM,NDIM) :: jinv = 0.d0
    double precision :: jacobian = 0.d0

    !--- source orientation: rows are N, E, Z-up in Cartesian ---
    double precision, dimension(NDIM,NDIM) :: nu = 0.d0

    !--- geographic intermediates, kept for reporting and for Stage 8 ---
    double precision :: theta = 0.d0       ! geocentric colatitude, radians
    double precision :: phi   = 0.d0       ! longitude, radians
    double precision :: r_surface = 0.d0   ! surface radius above the source

    !--- diagnostics ---
    ! |mapped - target|, in km: the analogue of the solver's
    ! "Error in location of the source" (locate_sources.f90:615)
    double precision :: distance_km = 0.d0
    ! worst 27-anchor reconstruction residual on the accepted element
    double precision :: anchor_err = 0.d0
  end type t_gf_location

  !-----------------------------------------------------------------
  ! a seismic source
  !
  ! Filled by gf_read_force_source() / gf_read_cmt_source(), which wrap
  ! src/specfem3D/get_force.f90 and get_cmt.f90 rather than re-parsing the
  ! files here. That is deliberate: those readers also apply the
  ! non-dimensionalisation (scaleF, scaleM), so the scaling comes from the
  ! solver's own source and cannot drift away from it.
  !-----------------------------------------------------------------

  integer, parameter :: GF_SRC_FORCE = 1
  integer, parameter :: GF_SRC_CMT   = 2

  type :: t_gf_source
    integer :: source_type = 0                     ! GF_SRC_FORCE or GF_SRC_CMT
    character(len=MAX_STRING_LEN) :: filename = ''

    double precision :: latitude  = 0.d0
    double precision :: longitude = 0.d0
    double precision :: depth     = 0.d0           ! km, as the files give it

    ! get_cmt returns the *triangle* half duration; the conversion to a
    ! Gaussian width is hdur/SOURCE_DECAY_MIMIC_TRIANGLE and happens in
    ! setup_sources_receivers.f90:777, not in the reader. For a force
    ! source get_force stores the dominant frequency f0 here instead.
    double precision :: hdur = 0.d0

    ! get_cmt/get_force zero this when NSOURCES == 1 and return the original
    ! in min_tshift_src_original, so t = 0 is the centroid time. Those
    ! seconds are origin-time metadata (Stage 10's SAC headers), not part
    ! of the trace.
    double precision :: tshift_src = 0.d0
    double precision :: min_tshift_src_original = 0.d0

    !--- force sources ---
    integer :: force_stf = 0
    double precision :: factor_force_source = 0.d0     ! non-dimensional, /scaleF
    double precision :: comp_dir_vect_source_E    = 0.d0
    double precision :: comp_dir_vect_source_N    = 0.d0
    double precision :: comp_dir_vect_source_Z_UP = 0.d0

    !--- moment tensor sources (Stage 4) ---
    ! spherical (Mrr,Mtt,Mpp,Mrt,Mrp,Mtp), non-dimensional, /scaleM
    double precision, dimension(6) :: moment_tensor = 0.d0

    !--- PDE origin time, from the CMTSOLUTION header (Stage 10) ---
    integer :: yr = 0, jda = 0, mo = 0, da = 0, ho = 0, mi = 0
    double precision :: sec = 0.d0
  end type t_gf_source

  !-----------------------------------------------------------------
  ! last error message
  !
  ! Set by gf_set_error() alongside the returned code, so a caller that
  ! wants a human-readable reason has one without the library printing
  ! anything itself.
  !-----------------------------------------------------------------

  character(len=MAX_STRING_LEN) :: gf_errmsg = ''

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_set_error(ierr_out,code,msg)

! records an error and its message, and returns the code

  implicit none

  integer, intent(out) :: ierr_out
  integer, intent(in) :: code
  character(len=*), intent(in) :: msg

  ierr_out = code
  gf_errmsg = msg

  end subroutine gf_set_error

!
!-------------------------------------------------------------------------------------------------
!

  function gf_error_string(code) result(str)

! short name of an error code, for callers that report the code rather than the message

  implicit none

  integer, intent(in) :: code
  character(len=24) :: str

  select case (code)
  case (GF_OK)             ; str = 'ok'
  case (GF_ERR_NO_HDF5)    ; str = 'not built with HDF5'
  case (GF_ERR_NO_PATH)    ; str = 'database not found'
  case (GF_ERR_NO_FILE)    ; str = 'file not found'
  case (GF_ERR_HDF5)       ; str = 'HDF5 error'
  case (GF_ERR_IO)         ; str = 'I/O error'
  case (GF_ERR_FORMAT)     ; str = 'malformed database'
  case (GF_ERR_MISMATCH)   ; str = 'incompatible database'
  case (GF_ERR_INCOMPLETE) ; str = 'incomplete database'
  case (GF_ERR_ALLOC)      ; str = 'allocation failed'
  case (GF_ERR_ARG)        ; str = 'invalid argument'
  case (GF_ERR_NO_ELEMENT) ; str = 'point not in database'
  case (GF_ERR_GEOMETRY)   ; str = 'degenerate geometry'
  case default             ; str = 'unknown error'
  end select

  end function gf_error_string

  end module gf_par
