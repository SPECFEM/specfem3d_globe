!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

! XCOMBINE_SEM
!
! USAGE
!   mpirun -np NPROC bin/xcombine_sem KERNEL_NAMES INPUT_FILE OUTPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   KERNEL_NAMES           - one or more material parameter names separated by commas
!   INPUT_FILE             - text file containing list of kernel directories
!   OUTPUT_PATH            - directory to which summed kernels are written
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, sums kernels from directories specified in
!   INPUT_FILE. Writes the resulting sum to OUTPUT_DIR.
!
!   INPUT_FILE is a text file containing a list of absolute or relative paths to
!   kernel directories, one directory per line.
!
!   KERNEL_NAMES is comma-delimited list of kernel names,
!   e.g.'alpha_kernel,beta_kernel,rho_kernel'.
!
!   This program's primary use case is to sum kernels. It can be used though on
!   any "reg1" array, i.e. any array of dimension
!   (NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE).
!
!   This is a parallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarrassingly-parallel
!   fashion.



program combine_sem_globe

  use postprocess_par

  implicit none

  integer, parameter :: NARGS = 3

  character(len=MAX_STRING_LEN) :: kernel_paths(MAX_KERNEL_PATHS),kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: sline,output_dir,input_file,kernel_names_comma_delimited, kernel_name
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  integer :: i,ier,iker,npath,nker


  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! check command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
      print *, 'USAGE: mpirun -np NPROC bin/xcombine_sem KERNEL_NAMES INPUT_FILE OUTPUT_DIR'
      stop ' Please check command line arguments'
    endif
  endif
  call synchronize_all()

  ! check number of MPI processes
  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print *,''
      print *,'Expected number of MPI processes: ', NPROCTOT_VAL
      print *,'Actual number of MPI processes: ', sizeprocs
      print *,''
    endif
    call synchronize_all()
    stop 'Error wrong number of MPI processes'
  endif
  call synchronize_all()

  if(myrank==0) then
    write(*,*) 'Running XCOMBINE_SEM'
    write(*,*)
  endif

  ! parse command line arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
  enddo
  read(arg(1),'(a)') kernel_names_comma_delimited
  read(arg(2),'(a)') input_file
  read(arg(3),'(a)') output_dir

 ! parse names from KERNEL_NAMES
  call parse_kernel_names(kernel_names_comma_delimited, kernel_names, nker)

  ! user output
  if (myrank == 0) then
    print *,'kernel names     : ',(trim(kernel_names(i))//' ',i=1,nker)
    print *,'input file       : ',trim(input_file)
    print *,'output directory : ',trim(output_dir)
    print *,''
  endif

  ! parse paths from INPUT_FILE
  npath=0
  open(unit = IIN, file = trim(input_file), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(input_file),myrank
     stop 1
  endif
  do while (.true.)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     npath = npath+1
     if (npath > MAX_KERNEL_PATHS) stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     kernel_paths(npath) = sline
  enddo
  close(IIN)
  ! user output
  if (myrank == 0) then
    print *,'number of events: ',npath
    print *,''
    print *,'summing kernels in:'
    print *,(trim(kernel_paths(i))//' ',i=1,npath)
    print *,''
  endif

  call synchronize_all()

  ! call kernel summation subroutine
  do iker=1,nker
      kernel_name = kernel_names(iker) ! e.g. alpha_kernel, beta_kernel, rho_kernel
      call sum_kernel(kernel_names(iker),kernel_paths,output_dir,npath)
  enddo

  ! user output
  if(myrank==0) then
    print *,'done writing all kernels, see directory ',trim(output_dir)
  endif

  ! stop all the processes, and exit
  call finalize_mpi()

end program combine_sem_globe

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel(kernel_name,kernel_paths,output_dir,npath)

  use postprocess_par

  implicit none

  character(len=*),parameter :: reg = '_reg1_'

  character(len=MAX_STRING_LEN) :: output_dir,kernel_name,kernel_paths(MAX_KERNEL_PATHS)
  integer :: ipath,npath

  ! local parameters
  character(len=MAX_STRING_LEN) :: k_file
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: kernel,total_kernel
  double precision :: norm,norm_sum
  integer :: ier

  ! initializes arrays
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           total_kernel(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  ! loops over all paths
  total_kernel = 0._CUSTOM_REAL
  do ipath = 1, npath
    if(myrank==0) then
    write(*,*) 'reading in event kernel for: ',trim(kernel_name)
    write(*,*) '    ',ipath, ' out of ', npath
    endif

    ! read kernel
    kernel = 0._CUSTOM_REAL
    write(k_file,'(a,i6.6,a)') trim(kernel_paths(ipath))//'/proc',myrank,reg//trim(kernel_name)//'.bin'
    open(IIN,file=trim(k_file),status='old',form='unformatted',action='read',iostat=ier)
    if (ier /= 0) then
     write(*,*) '  kernel not found: ',trim(k_file)
     stop 'Error kernel file not found'
    endif
    read(IIN) kernel
    close(IIN)

    norm = sum( kernel * kernel )
    call sum_all_dp(norm,norm_sum)
    if (myrank == 0) then
      print *,'  norm kernel: ',sqrt(norm_sum)
      print *
    endif

    ! keep track of sum
    total_kernel = total_kernel + kernel
  enddo


  ! write sum
  if(myrank==0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)
  write(k_file,'(a,i6.6,a)') trim(output_dir)//'/'//'proc',myrank,reg//trim(kernel_name)//'.bin'

  open(IOUT,file=trim(k_file),form='unformatted',status='unknown',action='write',iostat=ier)
  if (ier /= 0) then
    write(*,*) 'Error kernel not written: ',trim(k_file)
    stop 'Error kernel write'
  endif
  write(IOUT) total_kernel
  close(IOUT)

  if(myrank==0) write(*,*)
  deallocate(kernel,total_kernel)

end subroutine sum_kernel


