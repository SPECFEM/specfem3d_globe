program test_locate

  use specfem_par
  use specfem_par_movie
  use manager_adios

  implicit none

  integer :: iproc

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call init_mpi()
  call world_rank(myrank)

  if (myrank == 0) print *,'program: test_locate'

  ! initializes simulation parameters
  call initialize_simulation()

  ! sets up reference element GLL points/weights/derivatives
  call setup_GLL_points()

  ! starts reading the databases
  call read_mesh_databases()

  ! reads topography & bathymetry & ellipticity
  call read_topography_bathymetry()

  ! prepares sources and receivers
  call setup_sources_receivers()

  call synchronize_all()

  if (myrank == 0) then
    print *
    print *,'Checks:'
    print *
  endif

  do iproc = 0,NPROCTOT-1
    if (myrank == iproc) then
      print *,'rank ',iproc
      print *

      print *,'  sources:'
      print *,'  NSOURCES = ',NSOURCES
      if (NSOURCES /= 7) then
        print *,'ERROR: rank ',myrank,' - NSOURCES value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif

      ! total Mw = 7.7356437096222379  ! old
      !            7.7356439003571005  ! new
      print *,'  Mw = ',Mw
      if (abs(Mw - 7.7356439003571005d0) > 1.e-6) then
        print *,'ERROR: rank ',myrank,' - Mw value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif

      if (myrank == 0) then
        print *,'  final distance max = ',source_final_distance_max
        if (source_final_distance_max > 1.e-10) then
          print *,'Error: source final distance max too large'
          stop 1
        else
          print *,'  result is correct'
        endif
      endif
      print *

      print *,'  receivers:'
      print *,'  nrec = ',nrec
      if (nrec /= 10) then
        print *,'ERROR: rank ',myrank,' - nrec value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif

      if (myrank == 0) then
        print *,'  final distance max = ',receiver_final_distance_max
        if (receiver_final_distance_max > 1.e-10) then
          print *,'Error: receiver final distance max too large'
          stop 1
        else
          print *,'  result is correct'
        endif
      endif
      print *

    endif
    call synchronize_all()
  enddo

  call synchronize_all()

  ! done
  if (myrank == 0) print *,'test_locate done successfully'

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program test_locate

