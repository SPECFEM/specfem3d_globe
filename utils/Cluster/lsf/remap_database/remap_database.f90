program remap_databases

implicit none

! standard include of the MPI library
  include 'mpif.h'

  integer, parameter :: MAX_PROCS = 1000
  integer ier,sizeprocs,myrank,ios,i
  character(len=150) old_machine_file,junk,junk2,junk3
  character(len=150) slice_to_old_machine(MAX_PROCS),mymachine, &
  old_local_data_base, new_local_data_base,scp_outfile, command_string

  integer num_slices, num_slices2,num, old_jobid, new_jobid
  logical use_jobid

! initialize the MPI communicator and start the NPROCTOT MPI processes.
  call MPI_INIT(ier)

! run the main program
  call MPI_COMM_SIZE(MPI_COMM_WORLD,sizeprocs,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

  call getarg(1,old_machine_file)
  call getarg(2,junk)
  if (trim(old_machine_file) == '' .or. trim(junk) == '') call exit_mpi(myrank,'Usage: remap old-mach num-slice [old-jobid new-jobid]')
  read(junk,*) num_slices

  call getarg(3,junk2)
  if (trim(junk2) == '') then
     use_jobid=.false.
  else
     call getarg(4,junk3)
     if (trim(junk3) == '') call exit_mpi(myrank,'Usage: remap old-mach num-slice [old-jobid new-jobid]')
     read(junk2,*) old_jobid
     read(junk3,*) new_jobid
     use_jobid=.true.
  endif
  if (num_slices /= sizeprocs) call exit_mpi(myrank,'number of slices does not match')

  num_slices2 = num_slices
  open(11,file=trim(old_machine_file),status='old',iostat=ios)
  if (ios /= 0) stop 'Error opening old machine file'
  do while (1 == 1)
    read(11,'(a)',iostat=ios) junk2
    if (ios /= 0) exit
    read(11,'(i)',iostat=ios) num
    if (ios /= 0) exit
    do i = 1, num
      slice_to_old_machine(num_slices2-i+1) = junk2
    enddo
    num_slices2 = num_slices2 - num
  enddo
  if (num_slices2 /= 0) stop 'Error counting number of slices'
  close(11)

  mymachine = slice_to_old_machine(myrank+1)

  if (use_jobid) then
    write(old_local_data_base,'(a,i0)') '/scratch/lqy/DATABASES_MPI.',old_jobid
    write(new_local_data_base,'(a,i0)') '/scratch/lqy/DATABASES_MPI.',new_jobid
  else
    old_local_data_base = '/scratch/lqy/DATABASES_MPI'
    new_local_data_base = '/scratch/lqy/DATABASES_MPI'
  endif

  write(scp_outfile,'(a,i4.4)') 'OUTPUT_FILES/scp_out.',myrank

  write(command_string,'(a,i4.4,a)') 'scp lqy@'//trim(mymachine)//':'//trim(old_local_data_base)//'/*', &
             myrank, '*  '//trim(new_local_data_base)

!  call system('echo '//trim(command_string)//' > '//trim(scp_outfile))

  call system(trim(command_string)) !//' >> '//trim(scp_outfile))

! stop all the MPI processes, and exit
  call MPI_FINALIZE(ier)

end program remap_databases

