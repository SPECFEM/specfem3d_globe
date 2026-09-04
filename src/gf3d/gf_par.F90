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

  use constants, only: MAX_STRING_LEN,MAX_LENGTH_STATION_NAME,MAX_LENGTH_NETWORK_NAME

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
  case default             ; str = 'unknown error'
  end select

  end function gf_error_string

  end module gf_par
