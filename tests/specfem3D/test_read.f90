program test_read

  use specfem_par
  use specfem_par_movie
  use manager_adios

  implicit none

  ! local parameters
  integer :: iproc

  ! some test values from the default Par_file
  integer,parameter :: PAR_FILE_NPROC = 4
  double precision,parameter :: PAR_FILE_RECORD_LENGTH_IN_MINUTES = 2.5
  integer,parameter :: PAR_FILE_NTSTEP_BETWEEN_OUTPUT_INFO = 500
  logical,parameter :: PAR_FILE_RECEIVERS_CAN_BE_BURIED = .true.

  ! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call init_mpi()
  call world_rank(myrank)

  if (myrank == 0) print *,'program: test_read'

  ! initializes simulation parameters
  call initialize_simulation()

  call synchronize_all()

  if (myrank == 0) then
    print *
    print *,'Checks:'
    print *
  endif

  do iproc = 0,NPROCTOT-1
    if (myrank == iproc) then
      print *,'rank: ',iproc
      print *

      ! punctual check of values for given default Par_file in SPECFEM3D/DATA/ directory
      print *,'  NPROC = ',NPROCTOT
      if (NPROCTOT /= PAR_FILE_NPROC) then
        print *,'ERROR: NPROCTOT value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif

      print *,'  RECORD_LENGTH_IN_MINUTES = ',RECORD_LENGTH_IN_MINUTES
      if (abs(RECORD_LENGTH_IN_MINUTES - PAR_FILE_RECORD_LENGTH_IN_MINUTES) > 1.e-9) then
        print *,'ERROR: RECORD_LENGTH_IN_MINUTES value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif

      print *,'  NTSTEP_BETWEEN_OUTPUT_INFO = ',NTSTEP_BETWEEN_OUTPUT_INFO
      if (NTSTEP_BETWEEN_OUTPUT_INFO /= PAR_FILE_NTSTEP_BETWEEN_OUTPUT_INFO) then
        print *,'ERROR: NTSTEP_BETWEEN_OUTPUT_INFO value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif

      print *,'  RECEIVERS_CAN_BE_BURIED = ',RECEIVERS_CAN_BE_BURIED
      if (RECEIVERS_CAN_BE_BURIED .neqv. PAR_FILE_RECEIVERS_CAN_BE_BURIED) then
        print *,'ERROR: RECEIVERS_CAN_BE_BURIED value invalid'
        stop 1
      else
        print *,'  result is correct'
      endif
      print *
    endif
    call synchronize_all()
  enddo

  ! done
  if (myrank == 0) print *,'test_read done successfully'

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program test_read

