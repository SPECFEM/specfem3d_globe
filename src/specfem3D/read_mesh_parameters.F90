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


  subroutine read_mesh_parameters()

! reads in file with all mesh parameters

#ifdef USE_STATIC_COMPILATION
  ! static compilation
  ! no need to read parameter file, array values are included in constants_solver through values_from_mesher.h
  return

#else
  ! for dynamic compilation/allocation of arrays

  use constants, only: MAX_STRING_LEN,IIN

  use shared_parameters, only: LOCAL_PATH,SIMULATION_TYPE,SAVE_FORWARD, &
    ABSORBING_CONDITIONS,ATTENUATION, &
    MOVIE_VOLUME,MOVIE_VOLUME_TYPE, &
    OCEANS,ROTATION,FULL_GRAVITY

  use constants_solver

  implicit none

  ! local parameters
  integer :: ier
  ! processor identification
  character(len=MAX_STRING_LEN) :: filename

  if (.not. I_should_read_the_database) return

  ! create full name with path
  filename = trim(LOCAL_PATH) // "/mesh_parameters.bin"

  open(unit=IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) then
    print *,'Error opening file: ',trim(filename)
    stop 'Error opening mesh_parameters.bin'
  endif

  ! same values from mesher as would have been saved by save_header_file.F90

  read(IIN) NEX_XI_VAL
  read(IIN) NEX_ETA_VAL

  read(IIN) NSPEC_CRUST_MANTLE   ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
  read(IIN) NSPEC_OUTER_CORE     ! NSPEC_REGIONS(IREGION_OUTER_CORE)
  read(IIN) NSPEC_INNER_CORE     ! NSPEC_REGIONS(IREGION_INNER_CORE)
  read(IIN) NSPEC_TRINFINITE     ! NSPEC_REGIONS(IREGION_TRINFINITE)
  read(IIN) NSPEC_INFINITE       ! NSPEC_REGIONS(IREGION_INFINITE)

  read(IIN) NGLOB_CRUST_MANTLE   ! NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  read(IIN) NGLOB_OUTER_CORE     ! NGLOB_REGIONS(IREGION_OUTER_CORE)
  read(IIN) NGLOB_INNER_CORE     ! NGLOB_REGIONS(IREGION_INNER_CORE)
  read(IIN) NGLOB_TRINFINITE     ! NGLOB_REGIONS(IREGION_TRINFINITE)
  read(IIN) NGLOB_INFINITE       ! NGLOB_REGIONS(IREGION_INFINITE)

  read(IIN) NSPECMAX_ANISO_IC

  read(IIN) NSPECMAX_ISO_MANTLE
  read(IIN) NSPECMAX_TISO_MANTLE
  read(IIN) NSPECMAX_ANISO_MANTLE

  ! if attenuation is off, set dummy size of arrays to one
  if (ATTENUATION) then
    NSPEC_CRUST_MANTLE_ATTENUATION = NSPEC_CRUST_MANTLE   ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_ATTENUATION = NSPEC_INNER_CORE       ! NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_ATTENUATION = 0
    NSPEC_INNER_CORE_ATTENUATION = 0
  endif

  if (ATTENUATION .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    NSPEC_CRUST_MANTLE_STR_OR_ATT = NSPEC_CRUST_MANTLE    ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STR_OR_ATT = NSPEC_INNER_CORE        ! NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STR_OR_ATT = 0
    NSPEC_INNER_CORE_STR_OR_ATT = 0
  endif

  if (ATTENUATION .and. &
    ( SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    NSPEC_CRUST_MANTLE_STR_AND_ATT = NSPEC_CRUST_MANTLE   ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STR_AND_ATT = NSPEC_INNER_CORE       ! NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STR_AND_ATT = 0
    NSPEC_INNER_CORE_STR_AND_ATT = 0
  endif

  if (SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = NSPEC_CRUST_MANTLE   ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STRAIN_ONLY = NSPEC_INNER_CORE       ! NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = 0
    NSPEC_INNER_CORE_STRAIN_ONLY = 0
  endif

  ! full gravity region sizes
  if (.not. FULL_GRAVITY) then
    NSPEC_TRINFINITE = 0
    NGLOB_TRINFINITE = 0

    NSPEC_INFINITE = 0
    NGLOB_INFINITE = 0
  endif

  ! adjoint sizes
  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    NSPEC_CRUST_MANTLE_ADJOINT = NSPEC_CRUST_MANTLE       ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_OUTER_CORE_ADJOINT = NSPEC_OUTER_CORE           ! NSPEC_REGIONS(IREGION_OUTER_CORE)
    NSPEC_INNER_CORE_ADJOINT = NSPEC_INNER_CORE           ! NSPEC_REGIONS(IREGION_INNER_CORE)
    NSPEC_TRINFINITE_ADJOINT = NSPEC_TRINFINITE
    NSPEC_INFINITE_ADJOINT   = NSPEC_INFINITE

    NGLOB_CRUST_MANTLE_ADJOINT = NGLOB_CRUST_MANTLE       ! NGLOB_REGIONS(IREGION_CRUST_MANTLE)
    NGLOB_OUTER_CORE_ADJOINT = NGLOB_OUTER_CORE           ! NGLOB_REGIONS(IREGION_OUTER_CORE)
    NGLOB_INNER_CORE_ADJOINT = NGLOB_INNER_CORE           ! NGLOB_REGIONS(IREGION_INNER_CORE)
    NGLOB_TRINFINITE_ADJOINT = NGLOB_TRINFINITE
    NGLOB_INFINITE_ADJOINT   = NGLOB_INFINITE

    if (ROTATION) then
      NSPEC_OUTER_CORE_ROT_ADJOINT = NSPEC_OUTER_CORE     ! NSPEC_REGIONS(IREGION_OUTER_CORE)
    else
      NSPEC_OUTER_CORE_ROT_ADJOINT = 0
    endif
  else
    NSPEC_CRUST_MANTLE_ADJOINT = 0
    NSPEC_OUTER_CORE_ADJOINT = 0
    NSPEC_INNER_CORE_ADJOINT = 0
    NSPEC_TRINFINITE_ADJOINT = 0
    NSPEC_INFINITE_ADJOINT   = 0

    NGLOB_CRUST_MANTLE_ADJOINT = 0
    NGLOB_OUTER_CORE_ADJOINT = 0
    NGLOB_INNER_CORE_ADJOINT = 0
    NGLOB_TRINFINITE_ADJOINT = 0
    NGLOB_INFINITE_ADJOINT   = 0

    NSPEC_OUTER_CORE_ROT_ADJOINT = 0
  endif

  ! if absorbing conditions are off, set dummy size of arrays to one
  if (ABSORBING_CONDITIONS) then
    NSPEC_CRUST_MANTLE_STACEY = NSPEC_CRUST_MANTLE   ! NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_OUTER_CORE_STACEY = NSPEC_OUTER_CORE       ! NSPEC_REGIONS(IREGION_OUTER_CORE)
  else
    NSPEC_CRUST_MANTLE_STACEY = 0
    NSPEC_OUTER_CORE_STACEY = 0
  endif

  ! if oceans are off, set dummy size of arrays to one
  if (OCEANS) then
    NGLOB_CRUST_MANTLE_OCEANS = NGLOB_CRUST_MANTLE   ! NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  else
    NGLOB_CRUST_MANTLE_OCEANS = 0
  endif

  ! this to allow for code elimination by the compiler in the solver for performance
  read(IIN) TRANSVERSE_ISOTROPY_VAL

  read(IIN) ANISOTROPIC_3D_MANTLE_VAL
  read(IIN) ANISOTROPIC_INNER_CORE_VAL

  read(IIN) ATTENUATION_VAL
  read(IIN) ATTENUATION_3D_VAL

  read(IIN) ELLIPTICITY_VAL
  read(IIN) GRAVITY_VAL
  read(IIN) OCEANS_VAL

  ! full gravity support
  read(IIN) FULL_GRAVITY_VAL

  read(IIN) NX_BATHY_VAL
  read(IIN) NY_BATHY_VAL

  read(IIN) ROTATION_VAL
  read(IIN) EXACT_MASS_MATRIX_FOR_ROTATION_VAL

  if (ROTATION) then
    NSPEC_OUTER_CORE_ROTATION = NSPEC_OUTER_CORE     ! NSPEC_REGIONS(IREGION_OUTER_CORE)
  else
    NSPEC_OUTER_CORE_ROTATION = 0
  endif

  read(IIN) PARTIAL_PHYS_DISPERSION_ONLY_VAL

  read(IIN) NPROC_XI_VAL
  read(IIN) NPROC_ETA_VAL
  read(IIN) NCHUNKS_VAL
  read(IIN) NPROCTOT_VAL

  read(IIN) ATT1_VAL
  read(IIN) ATT2_VAL
  read(IIN) ATT3_VAL
  read(IIN) ATT4_VAL
  read(IIN) ATT5_VAL

  read(IIN) NSPEC2DMAX_XMIN_XMAX_CM      ! NSPEC2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE)
  read(IIN) NSPEC2DMAX_YMIN_YMAX_CM      ! NSPEC2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE)
  read(IIN) NSPEC2D_BOTTOM_CM            ! NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)
  read(IIN) NSPEC2D_TOP_CM               ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)

  read(IIN) NSPEC2DMAX_XMIN_XMAX_IC      ! NSPEC2DMAX_XMIN_XMAX(IREGION_INNER_CORE)
  read(IIN) NSPEC2DMAX_YMIN_YMAX_IC      ! NSPEC2DMAX_YMIN_YMAX(IREGION_INNER_CORE)
  read(IIN) NSPEC2D_BOTTOM_IC            ! NSPEC2D_BOTTOM(IREGION_INNER_CORE)
  read(IIN) NSPEC2D_TOP_IC               ! NSPEC2D_TOP(IREGION_INNER_CORE)

  read(IIN) NSPEC2DMAX_XMIN_XMAX_OC      ! NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE)
  read(IIN) NSPEC2DMAX_YMIN_YMAX_OC      ! NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE)
  read(IIN) NSPEC2D_BOTTOM_OC            ! NSPEC2D_BOTTOM(IREGION_OUTER_CORE)
  read(IIN) NSPEC2D_TOP_OC               ! NSPEC2D_TOP(IREGION_OUTER_CORE)

  read(IIN) NSPEC2DMAX_XMIN_XMAX_TRINF   ! NSPEC2DMAX_XMIN_XMAX(IREGION_TRINFINITE)
  read(IIN) NSPEC2DMAX_YMIN_YMAX_TRINF   ! NSPEC2DMAX_YMIN_YMAX(IREGION_TRINFINITE)
  read(IIN) NSPEC2D_BOTTOM_TRINF         ! NSPEC2D_BOTTOM(IREGION_TRINFINITE)
  read(IIN) NSPEC2D_TOP_TRINF            ! NSPEC2D_TOP(IREGION_TRINFINITE)

  read(IIN) NSPEC2DMAX_XMIN_XMAX_INF     ! NSPEC2DMAX_XMIN_XMAX(IREGION_INFINITE)
  read(IIN) NSPEC2DMAX_YMIN_YMAX_INF     ! NSPEC2DMAX_YMIN_YMAX(IREGION_INFINITE)
  read(IIN) NSPEC2D_BOTTOM_INF           ! NSPEC2D_BOTTOM(IREGION_INFINITE)
  read(IIN) NSPEC2D_TOP_INF              ! NSPEC2D_TOP(IREGION_INFINITE)

  ! for boundary kernels
  read(IIN) NSPEC2D_MOHO                 ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
  read(IIN) NSPEC2D_400                  ! NSPEC2D_MOHO / 4
  read(IIN) NSPEC2D_670                  ! NSPEC2D_400
  read(IIN) NSPEC2D_CMB                  ! NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)
  read(IIN) NSPEC2D_ICB                  ! NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

  ! boundary kernels only needed for kernel simulations
  if (.not. (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3)) then
    NSPEC2D_MOHO = 0
    NSPEC2D_400 = 0
    NSPEC2D_670 = 0
    NSPEC2D_CMB = 0
    NSPEC2D_ICB = 0
  endif

  ! Deville routines only implemented for NGLLX = NGLLY = NGLLZ = 5
  if (NGLLX == 5 .and. NGLLY == 5 .and. NGLLZ == 5) then
    USE_DEVILLE_PRODUCTS_VAL = .true.
  else
    USE_DEVILLE_PRODUCTS_VAL = .false.
  endif

  if (MOVIE_VOLUME) then
    NSPEC_CRUST_MANTLE_3DMOVIE = NSPEC_CRUST_MANTLE
    NGLOB_CRUST_MANTLE_3DMOVIE = NGLOB_CRUST_MANTLE
  else
    NSPEC_CRUST_MANTLE_3DMOVIE = 0
    NGLOB_CRUST_MANTLE_3DMOVIE = 0
  endif

  if (MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE == 4) then
    NSPEC_OUTER_CORE_3DMOVIE = NSPEC_OUTER_CORE
  else
    NSPEC_OUTER_CORE_3DMOVIE = 0
  endif

  ! in the case of Stacey boundary conditions, add C*delta/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be fictitious / unused

  read(IIN) NGLOB_XY_CM
  read(IIN) NGLOB_XY_IC

  read(IIN) ATTENUATION_1D_WITH_3D_STORAGE_VAL

  ! for UNDO_ATTENUATION
  read(IIN) UNDO_ATTENUATION_VAL
  read(IIN) NT_DUMP_ATTENUATION_VAL    ! NT_DUMP_ATTENUATION_optimal

  ! mesh geometry (with format specifier to avoid writing double values on a newline)
  read(IIN) ANGULAR_WIDTH_ETA_IN_DEGREES_VAL
  read(IIN) ANGULAR_WIDTH_XI_IN_DEGREES_VAL
  read(IIN) CENTER_LATITUDE_IN_DEGREES_VAL
  read(IIN) CENTER_LONGITUDE_IN_DEGREES_VAL
  read(IIN) GAMMA_ROTATION_AZIMUTH_VAL

  close(IIN)

#endif

  end subroutine read_mesh_parameters

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_parameters()

#ifdef USE_STATIC_COMPILATION
  ! static compilation
  ! no need to broadcast, each process has array values included already during compilation
  return

#else
  ! for dynamic compilation/allocation of arrays
  ! (no need for values_from_mesher.h file to compile)

  use constants, only: myrank
  use constants_solver

  implicit none

  ! local parameters
  ! broadcast parameter arrays
  integer, parameter :: nparam_i = 81
  integer, dimension(nparam_i) :: bcast_integer

  integer, parameter :: nparam_l = 15
  logical, dimension(nparam_l) :: bcast_logical

  integer, parameter :: nparam_dp = 5
  double precision, dimension(nparam_dp) :: bcast_double_precision

  ! initializes containers
  bcast_integer(:) = 0
  bcast_logical(:) = .false.
  bcast_double_precision(:) = 0.d0

  ! main process prepares broadcasting arrays
  if (I_should_read_the_database) then
    if (myrank == 0) then
      ! simple way to pass parameters in arrays from main to all other processes
      ! rather than single values one by one to reduce MPI communication calls:
      ! sets up broadcasting array
      bcast_integer = (/ &
        NEX_XI_VAL,NEX_ETA_VAL, &
        NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE, &
        NGLOB_CRUST_MANTLE,NGLOB_OUTER_CORE,NGLOB_INNER_CORE, &
        NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE,NSPECMAX_ANISO_MANTLE, &
        NSPEC_CRUST_MANTLE_ATTENUATION,NSPEC_INNER_CORE_ATTENUATION, &
        NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
        NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
        NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
        NSPEC_CRUST_MANTLE_ADJOINT,NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
        NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT,NGLOB_INNER_CORE_ADJOINT, &
        NSPEC_OUTER_CORE_ROT_ADJOINT, &
        NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
        NGLOB_CRUST_MANTLE_OCEANS, &
        NX_BATHY_VAL,NY_BATHY_VAL, &
        NSPEC_OUTER_CORE_ROTATION, &
        NPROC_XI_VAL,NPROC_ETA_VAL,NCHUNKS_VAL,NPROCTOT_VAL, &
        ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL,ATT5_VAL, &
        NSPEC2DMAX_XMIN_XMAX_CM,NSPEC2DMAX_YMIN_YMAX_CM,NSPEC2D_BOTTOM_CM,NSPEC2D_TOP_CM, &
        NSPEC2DMAX_XMIN_XMAX_IC,NSPEC2DMAX_YMIN_YMAX_IC,NSPEC2D_BOTTOM_IC,NSPEC2D_TOP_IC, &
        NSPEC2DMAX_XMIN_XMAX_OC,NSPEC2DMAX_YMIN_YMAX_OC,NSPEC2D_BOTTOM_OC,NSPEC2D_TOP_OC, &
        NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,NSPEC2D_CMB,NSPEC2D_ICB, &
        NSPEC_CRUST_MANTLE_3DMOVIE,NGLOB_CRUST_MANTLE_3DMOVIE,NSPEC_OUTER_CORE_3DMOVIE, &
        NGLOB_XY_CM,NGLOB_XY_IC,NT_DUMP_ATTENUATION_VAL, &
        NSPEC_TRINFINITE,NSPEC_INFINITE, &
        NGLOB_TRINFINITE,NGLOB_INFINITE, &
        NSPEC_TRINFINITE_ADJOINT,NSPEC_INFINITE_ADJOINT, &
        NGLOB_TRINFINITE_ADJOINT,NGLOB_INFINITE_ADJOINT, &
        NSPEC2DMAX_XMIN_XMAX_TRINF,NSPEC2DMAX_YMIN_YMAX_TRINF,NSPEC2D_BOTTOM_TRINF,NSPEC2D_TOP_TRINF, &
        NSPEC2DMAX_XMIN_XMAX_INF,NSPEC2DMAX_YMIN_YMAX_INF,NSPEC2D_BOTTOM_INF,NSPEC2D_TOP_INF &
        /)

      bcast_logical = (/ &
        TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL,ANISOTROPIC_INNER_CORE_VAL, &
        ATTENUATION_VAL,ATTENUATION_3D_VAL, &
        ELLIPTICITY_VAL,GRAVITY_VAL,OCEANS_VAL, &
        ROTATION_VAL,EXACT_MASS_MATRIX_FOR_ROTATION_VAL, &
        PARTIAL_PHYS_DISPERSION_ONLY_VAL, &
        USE_DEVILLE_PRODUCTS_VAL, &
        ATTENUATION_1D_WITH_3D_STORAGE_VAL,UNDO_ATTENUATION_VAL, &
        FULL_GRAVITY_VAL &
        /)

      bcast_double_precision = (/ &
        ANGULAR_WIDTH_ETA_IN_DEGREES_VAL,ANGULAR_WIDTH_XI_IN_DEGREES_VAL, &
        CENTER_LATITUDE_IN_DEGREES_VAL,CENTER_LONGITUDE_IN_DEGREES_VAL, &
        GAMMA_ROTATION_AZIMUTH_VAL &
        /)
    endif

    ! broadcasts the information read on main rank to the nodes (within same group)
    call bcast_all_i(bcast_integer,nparam_i)
    call bcast_all_l(bcast_logical,nparam_l)
    call bcast_all_dp(bcast_double_precision,nparam_dp)
  endif

  ! broadcasts the information read on the main to the other run groups
  call bcast_all_i_for_database(bcast_integer(1),size(bcast_integer,kind=4))
  call bcast_all_l_for_database(bcast_logical(1),size(bcast_logical,kind=4))
  call bcast_all_dp_for_database(bcast_double_precision(1),size(bcast_double_precision,kind=4))

  ! non-main processes set their parameters
  if (myrank /= 0 .or. (.not. I_should_read_the_database)) then
    ! please, be careful with ordering and counting here
    ! integers
    NEX_XI_VAL = bcast_integer(1)
    NEX_ETA_VAL = bcast_integer(2)
    NSPEC_CRUST_MANTLE = bcast_integer(3)
    NSPEC_OUTER_CORE = bcast_integer(4)
    NSPEC_INNER_CORE = bcast_integer(5)
    NGLOB_CRUST_MANTLE = bcast_integer(6)
    NGLOB_OUTER_CORE = bcast_integer(7)
    NGLOB_INNER_CORE = bcast_integer(8)
    NSPECMAX_ANISO_IC = bcast_integer(9)
    NSPECMAX_ISO_MANTLE = bcast_integer(10)
    NSPECMAX_TISO_MANTLE = bcast_integer(11)
    NSPECMAX_ANISO_MANTLE = bcast_integer(12)
    NSPEC_CRUST_MANTLE_ATTENUATION = bcast_integer(13)
    NSPEC_INNER_CORE_ATTENUATION = bcast_integer(14)
    NSPEC_CRUST_MANTLE_STR_OR_ATT = bcast_integer(15)
    NSPEC_INNER_CORE_STR_OR_ATT = bcast_integer(16)
    NSPEC_CRUST_MANTLE_STR_AND_ATT = bcast_integer(17)
    NSPEC_INNER_CORE_STR_AND_ATT = bcast_integer(18)
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = bcast_integer(19)
    NSPEC_INNER_CORE_STRAIN_ONLY = bcast_integer(20)
    NSPEC_CRUST_MANTLE_ADJOINT = bcast_integer(21)
    NSPEC_OUTER_CORE_ADJOINT = bcast_integer(22)
    NSPEC_INNER_CORE_ADJOINT = bcast_integer(23)
    NGLOB_CRUST_MANTLE_ADJOINT = bcast_integer(24)
    NGLOB_OUTER_CORE_ADJOINT = bcast_integer(25)
    NGLOB_INNER_CORE_ADJOINT = bcast_integer(26)
    NSPEC_OUTER_CORE_ROT_ADJOINT = bcast_integer(27)
    NSPEC_CRUST_MANTLE_STACEY = bcast_integer(28)
    NSPEC_OUTER_CORE_STACEY = bcast_integer(29)
    NGLOB_CRUST_MANTLE_OCEANS = bcast_integer(30)
    NX_BATHY_VAL = bcast_integer(31)
    NY_BATHY_VAL = bcast_integer(32)
    NSPEC_OUTER_CORE_ROTATION = bcast_integer(33)
    NPROC_XI_VAL = bcast_integer(34)
    NPROC_ETA_VAL = bcast_integer(35)
    NCHUNKS_VAL = bcast_integer(36)
    NPROCTOT_VAL = bcast_integer(37)
    ATT1_VAL = bcast_integer(38)
    ATT2_VAL = bcast_integer(39)
    ATT3_VAL = bcast_integer(40)
    ATT4_VAL = bcast_integer(41)
    ATT5_VAL = bcast_integer(42)
    NSPEC2DMAX_XMIN_XMAX_CM = bcast_integer(43)
    NSPEC2DMAX_YMIN_YMAX_CM = bcast_integer(44)
    NSPEC2D_BOTTOM_CM = bcast_integer(45)
    NSPEC2D_TOP_CM = bcast_integer(46)
    NSPEC2DMAX_XMIN_XMAX_IC = bcast_integer(47)
    NSPEC2DMAX_YMIN_YMAX_IC = bcast_integer(48)
    NSPEC2D_BOTTOM_IC = bcast_integer(49)
    NSPEC2D_TOP_IC = bcast_integer(50)
    NSPEC2DMAX_XMIN_XMAX_OC = bcast_integer(51)
    NSPEC2DMAX_YMIN_YMAX_OC = bcast_integer(52)
    NSPEC2D_BOTTOM_OC = bcast_integer(53)
    NSPEC2D_TOP_OC = bcast_integer(54)
    NSPEC2D_MOHO = bcast_integer(55)
    NSPEC2D_400 = bcast_integer(56)
    NSPEC2D_670 = bcast_integer(57)
    NSPEC2D_CMB = bcast_integer(58)
    NSPEC2D_ICB = bcast_integer(59)
    NSPEC_CRUST_MANTLE_3DMOVIE = bcast_integer(60)
    NGLOB_CRUST_MANTLE_3DMOVIE = bcast_integer(61)
    NSPEC_OUTER_CORE_3DMOVIE = bcast_integer(62)
    NGLOB_XY_CM = bcast_integer(63)
    NGLOB_XY_IC = bcast_integer(64)
    NT_DUMP_ATTENUATION_VAL = bcast_integer(65)
    ! infinite element mesh
    NSPEC_TRINFINITE  = bcast_integer(66)
    NSPEC_INFINITE = bcast_integer(67)
    NGLOB_TRINFINITE = bcast_integer(68)
    NGLOB_INFINITE = bcast_integer(69)
    NSPEC_TRINFINITE_ADJOINT = bcast_integer(70)
    NSPEC_INFINITE_ADJOINT = bcast_integer(71)
    NGLOB_TRINFINITE_ADJOINT = bcast_integer(72)
    NGLOB_INFINITE_ADJOINT = bcast_integer(73)
    NSPEC2DMAX_XMIN_XMAX_TRINF = bcast_integer(74)
    NSPEC2DMAX_YMIN_YMAX_TRINF = bcast_integer(75)
    NSPEC2D_BOTTOM_TRINF = bcast_integer(76)
    NSPEC2D_TOP_TRINF = bcast_integer(77)
    NSPEC2DMAX_XMIN_XMAX_INF = bcast_integer(78)
    NSPEC2DMAX_YMIN_YMAX_INF = bcast_integer(79)
    NSPEC2D_BOTTOM_INF = bcast_integer(80)
    NSPEC2D_TOP_INF = bcast_integer(81)

    ! logicals
    TRANSVERSE_ISOTROPY_VAL = bcast_logical(1)
    ANISOTROPIC_3D_MANTLE_VAL = bcast_logical(2)
    ANISOTROPIC_INNER_CORE_VAL = bcast_logical(3)
    ATTENUATION_VAL = bcast_logical(4)
    ATTENUATION_3D_VAL = bcast_logical(5)
    ELLIPTICITY_VAL = bcast_logical(6)
    GRAVITY_VAL = bcast_logical(7)
    OCEANS_VAL = bcast_logical(8)
    ROTATION_VAL = bcast_logical(9)
    EXACT_MASS_MATRIX_FOR_ROTATION_VAL = bcast_logical(10)
    PARTIAL_PHYS_DISPERSION_ONLY_VAL = bcast_logical(11)
    USE_DEVILLE_PRODUCTS_VAL = bcast_logical(12)
    ATTENUATION_1D_WITH_3D_STORAGE_VAL = bcast_logical(13)
    UNDO_ATTENUATION_VAL = bcast_logical(14)
    FULL_GRAVITY_VAL = bcast_logical(15)

    ! double precisions
    ANGULAR_WIDTH_ETA_IN_DEGREES_VAL = bcast_double_precision(1)
    ANGULAR_WIDTH_XI_IN_DEGREES_VAL = bcast_double_precision(2)
    CENTER_LATITUDE_IN_DEGREES_VAL = bcast_double_precision(3)
    CENTER_LONGITUDE_IN_DEGREES_VAL = bcast_double_precision(4)
    GAMMA_ROTATION_AZIMUTH_VAL = bcast_double_precision(5)
  endif

#endif

  end subroutine bcast_mesh_parameters
