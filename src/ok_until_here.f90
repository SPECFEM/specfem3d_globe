
  print *,'ok until here'
  call flush(6) ! flush the screen buffer (supported at least in Intel ifort)
  call MPI_BARRIER(MPI_COMM_WORLD,ier)
! stop all the MPI processes, and exit
  call MPI_ABORT(MPI_COMM_WORLD,30,ier)
  stop 'end of the test'

