program test_read

  use meshfem_par
  use meshfem_models_par

!  use manager_adios

  implicit none

  ! some test values from the default Par_file
  integer,parameter :: PAR_FILE_NPROC = 4
  double precision,parameter :: PAR_FILE_RECORD_LENGTH_IN_MINUTES = 2.5
  integer,parameter :: PAR_FILE_NTSTEP_BETWEEN_OUTPUT_INFO = 500
  logical,parameter :: PAR_FILE_RECEIVERS_CAN_BE_BURIED = .true.

!  integer :: myrank

  ! reads DATA/Par_file
  myrank = 0

  ! reads the parameter file and computes additional parameters
  call read_compute_parameters()

  ! punctual check of values for given default Par_file in SPECFEM3D/DATA/ directory
  print *,'NPROC = ',NPROCTOT
  if (NPROCTOT /= PAR_FILE_NPROC) then
    print *,'ERROR: NPROCTOT value invalid'
    stop 1
  else
    print *,'  result is correct'
  endif

  print *,'RECORD_LENGTH_IN_MINUTES = ',RECORD_LENGTH_IN_MINUTES
  if (abs(RECORD_LENGTH_IN_MINUTES - PAR_FILE_RECORD_LENGTH_IN_MINUTES) > 1.e-9) then
    print *,'ERROR: RECORD_LENGTH_IN_MINUTES value invalid'
    stop 1
  else
    print *,'  result is correct'
  endif

  print *,'NTSTEP_BETWEEN_OUTPUT_INFO = ',NTSTEP_BETWEEN_OUTPUT_INFO
  if (NTSTEP_BETWEEN_OUTPUT_INFO /= PAR_FILE_NTSTEP_BETWEEN_OUTPUT_INFO) then
    print *,'ERROR: NTSTEP_BETWEEN_OUTPUT_INFO value invalid'
    stop 1
  else
    print *,'  result is correct'
  endif

  print *,'RECEIVERS_CAN_BE_BURIED = ',RECEIVERS_CAN_BE_BURIED
  if (RECEIVERS_CAN_BE_BURIED .neqv. PAR_FILE_RECEIVERS_CAN_BE_BURIED) then
    print *,'ERROR: RECEIVERS_CAN_BE_BURIED value invalid'
    stop 1
  else
    print *,'  result is correct'
  endif

  ! done
  print *,'test_read done successfully'

end program test_read

