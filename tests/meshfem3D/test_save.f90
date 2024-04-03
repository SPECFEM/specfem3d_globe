program test_save

  use meshfem_par
  use manager_adios

  implicit none

  include 'version.fh'

  ! local parameters
  integer :: sizeprocs
  ! for some statistics for the mesh
  integer :: numelem_crust_mantle,numelem_outer_core,numelem_inner_core
  integer :: numelem_trinfinite,numelem_infinite
  integer :: numelem_total

  ! timing
  double precision, external :: wtime

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call init_mpi()

! from initialize_mesher.f90

  call world_size(sizeprocs)
  call world_rank(myrank)

  if (myrank == 0) print *,'program: test_save'

  ! set the base pathname for output files
  OUTPUT_FILES = OUTPUT_FILES_BASE

  ! open main output file, only written to by process 0
  if (myrank == 0) then
    if (IMAIN /= ISTANDARD_OUTPUT) &
      open(unit=IMAIN,file=trim(OUTPUT_FILES)//'/output_mesher.txt',status='unknown')

    write(IMAIN,*)
    write(IMAIN,*) '*** Specfem3D MPI Mesher ***'
    write(IMAIN,*) 'Version: ', git_package_version
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

  ! initialize mesher
  if (myrank == 0) then
    ! reads the parameter file and computes additional parameters
    call read_parameter_file()

    ! overimposes a fixed model setup for testing
    MODEL = "transversely_isotropic_prem_plus_3D_crust_1.0"

    NCHUNKS = 1

    ANGULAR_WIDTH_XI_IN_DEGREES  = 90.d0    !angular size of a chunk
    ANGULAR_WIDTH_ETA_IN_DEGREES = 90.d0
    CENTER_LATITUDE_IN_DEGREES   = 90.d0
    CENTER_LONGITUDE_IN_DEGREES  = 0.d0
    GAMMA_ROTATION_AZIMUTH       = 0.d0

    NPROC_XI_read  = 2
    NPROC_ETA_read = 2

    NEX_XI_read    = 48
    NEX_ETA_read   = 48

    OCEANS      = .true.
    TOPOGRAPHY  = .true.
    ELLIPTICITY = .true.
    ATTENUATION = .true.
    GRAVITY     = .true.
    ROTATION    = .true.

    ABSORBING_CONDITIONS = .true.

    ! all other parameters are from the default (given through the configuration step)

    ! user output
    print *,'MODEL: ',trim(MODEL)
    print *

    ! sets compute parameters accordingly
    call rcp_set_compute_parameters()

    ! sets mesh parameters
    call rcp_set_mesh_parameters()
  endif

  ! broadcast parameters read from main process to all processes
  call broadcast_computed_parameters()

  ! tests further calls from initialize_mesher.f90
  call im_initialize_system()

  ! setup addressing and models
  call setup_model()

  ! creates meshes for regions crust/mantle, outer core and inner core
  call create_meshes()

  ! outputs mesh info and saves new header file
  call finalize_mesher()

  ! results only on main available
  if (myrank == 0) then
    print *
    print *,'Checks:'
    print *

    ! volume
    ! without topo/ellip: 0.69581709728628893d0
    ! with    topo/ellip: 0.69301991059962509
    !
    ! accuracy:
    !    gcc version 9.3: 0.69301991059962509 (MacOS)
    !    gcc version 5.4: 0.69301991059962254
    !    gcc version 7.5: 0.69301991064577684 (IBM Power)
    !
    ! new PREM model with inner core at 1221.5km instead of 1221.0km
    !    calculated volume:   0.69301991575060418
    !
    print *,'volume_total = ',volume_total
    ! more sensitive threshold (works only for gfortran tests)
    !if (abs(volume_total - 0.6930199157d0) > 1.d-10) then
    ! less sensitive threshold (works also on intel ifort tests)
    if (abs(volume_total - 0.6930199157d0) > 1.d-5) then
      print *,'volume expected: ',0.6930199157d0,' difference: ',abs(volume_total - 0.6930199157d0)
      print *,'ERROR: volume_total value invalid'
      stop 1
    else
      print *,'  result is correct'
    endif

    ! number of elements
    numelem_crust_mantle = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    numelem_outer_core   = NSPEC_REGIONS(IREGION_OUTER_CORE)
    numelem_inner_core   = NSPEC_REGIONS(IREGION_INNER_CORE)
    numelem_trinfinite   = NSPEC_REGIONS(IREGION_TRINFINITE)
    numelem_infinite     = NSPEC_REGIONS(IREGION_INFINITE)

    numelem_total = numelem_crust_mantle + numelem_outer_core + numelem_inner_core &
                    + numelem_trinfinite + numelem_infinite

    print *,'numelem_total = ',numelem_total
    if (numelem_total /= 5265) then
      print *,'ERROR: numelem_total value invalid'
      stop 1
    else
      print *,'  result is correct'
    endif
  endif

  ! done
  if (myrank == 0) print *,'test_save done successfully'

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program test_save


