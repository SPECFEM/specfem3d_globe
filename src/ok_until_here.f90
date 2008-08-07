
  print *,'ok until here'
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
! stop all the MPI processes, and exit
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'end of the test'

