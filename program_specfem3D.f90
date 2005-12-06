
  program xspecfem3D

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer ier

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

! run the main program
  call specfem3D

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

  end program xspecfem3D
