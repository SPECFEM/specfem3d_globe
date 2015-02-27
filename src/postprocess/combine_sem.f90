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
!   mpirun -np NPROC bin/xcombine_sem INPUT_FILE OUTPUT_DIR KERNELL_NAMES
!
!
! COMMAND LINE ARGUMENTS
!   INPUT_FILE             - text file containing list of kernel directories
!   OUTPUT_PATH            - directory to which summed kernels are written
!   KERNEL_NAMES           - one or more material parameter names separated by commas
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, sums kernels from directories specified in
!   INPUT_FILE. Writes the resulting sum to OUTPUT_DIR.
!
!   INPUT_FILE is a text file containing a list of absolute or relative paths to
!   kernel direcotires, one directoy per line.
!
!   KERNEL_NAMES is comma-delimited list of kernel names, 
!   e.g.'reg1_alpha_kernel,reg1_beta_kernel,reg1_rho_kernel'.
!
!   Currently only names beginning with 'reg1' are supported.
!
!   This program's primary use case is to sum kernels. It can be used though on
!   any scalar field of dimension (NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE). 
!
!   This is a parrallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.



program combine_sem

  use postprocess_par

  implicit none

  character(len=MAX_STRING_LEN) :: kernel_paths(MAX_KERNEL_PATHS),kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: sline,output_dir,input_file,kernel_names_comma_delimited, kernel_name
  character(len=MAX_STRING_LEN) :: arg(3)
  integer :: i,ier,iker,npath,nker


  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! parse command line arguments
  do i = 1, 3
    call get_command_argument(i,arg(i), status=ier)
    if (i <= 1 .and. trim(arg(i)) == '') then
      if (myrank == 0) then
      stop ' Reenter command line options'
      endif
    endif
  enddo
  read(arg(1),'(a)') input_file
  read(arg(2),'(a)') output_dir
  read(arg(3),'(a)') kernel_names_comma_delimited

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print*,''
      print*,'Error: run xcombine_sem with the same number of MPI processes '
      print*,'       as specified when slices were created'
      print*,''
      print*,'for example: mpirun -np ',NPROCTOT_VAL,' ./xcombine_sem ...'
      print*,''
    endif
    call synchronize_all()
    stop 'Error total number of slices'
  endif
  call synchronize_all()


  ! start execution
  if(myrank==0) then
    write(*,*) 'Running XCOMBINE_SEM'
    write(*,*)
  endif

 ! parse names from KERNEL_NAMES
  call parse_kernel_names(kernel_names_comma_delimited, kernel_names, nker)

  ! currently, only names beginning with 'reg1' supported
  if (myrank == 0) then
    do iker = 1, nker
      kernel_name = kernel_names(iker)
      if (kernel_name(1:4) /= 'reg1') stop 'Currently, only names beginning with reg1 supported.'
    enddo
  endif
  call synchronize_all()

  ! parse paths from INPUT_FILE
  npath=0
  open(unit = IIN, file = trim(input_file), status = 'old',iostat = ier)
  if (ier /= 0) then
     print *,'Error opening ',trim(input_file),myrank
     stop 1
  endif
  do while (1 == 1)
     read(IIN,'(a)',iostat=ier) sline
     if (ier /= 0) exit
     npath = npath+1
     if (npath > MAX_KERNEL_PATHS) stop 'Error number of paths exceeds MAX_KERNEL_PATHS'
     kernel_paths(npath) = sline
  enddo
  close(IIN)
  if (myrank == 0) then
    write(*,*) '  ',npath,' events'
    write(*,*)
  endif

  if(myrank == 0) then
    print*,'summing kernels in:'
    print*,kernel_paths(1:npath)
    print*
  endif

  call synchronize_all()

  ! call kernel summation subroutine
  do iker=1,nker
      kernel_name = kernel_names(iker) ! e.g. alpha_kernel, beta_kernel, rho_kernel
      call sum_kernel(kernel_names(iker),kernel_paths,output_dir,npath)
  enddo

  if(myrank==0) write(*,*) 'done writing all kernels, see directory ', output_dir

  ! stop all the processes, and exit
  call finalize_mpi()

end program combine_sem

!
!-------------------------------------------------------------------------------------------------
!

subroutine sum_kernel(kernel_name,kernel_paths,output_dir,npath)

  use postprocess_par

  implicit none

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
    write(k_file,'(a,i6.6,a)') trim(kernel_paths(ipath))//'/proc',myrank,'_'//trim(kernel_name)//'.bin'
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
      print*,'  norm kernel: ',sqrt(norm_sum)
      print*
    endif

    ! keep track of sum
    total_kernel = total_kernel + kernel
  enddo


  ! write sum
  if(myrank==0) write(*,*) 'writing out summed kernel for: ',trim(kernel_name)
  write(k_file,'(a,i6.6,a)') trim(output_dir)//'/'//'proc',myrank,'_'//trim(kernel_name)//'.bin'

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


