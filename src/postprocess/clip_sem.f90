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

! XCLIP_SEM
!
! USAGE
!   mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL KERNEL_NAMES INPUT_FILE OUTPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   MIN_VAL                - threshold below which array values are clipped
!   MAX_VAL                - threshold above which array values are clipped
!   KERNEL_NAMES           - one or more material parameter names separated by commas
!   INPUT_DIR              - directory from which arrays are read
!   OUTPUT_DIR             - directory to which clipped array are written
!
!
! DESCRIPTION
!   For each name in KERNEL_NAMES, reads kernels from INPUT_DIR, applies
!   thresholds, and writes the resulting clipped kernels to OUTPUT_DIR.
!
!   KERNEL_NAMES is comma-delimited list of material names,
!   e.g. 'reg1_alphav_kernel,reg1_alphah_kernel'
!
!   Files written to OUTPUT_DIR have the suffix 'clip' appended,
!   e.g. proc***alphav_kernel.bin becomes proc***alphav_kernel_clip.bin
!
!   This program's primary use case is to sum kernels. It can be used though on
!   any "reg1" array, i.e. any array of dimension
!   (NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE). The region suffix must be included
!   explicitly in all names supplied through the KERNEL_NAMES arugment.
!
!   This is a parrallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.



program clip_sem_globe

  use postprocess_par

  implicit none

  integer, parameter :: NARGS = 5

  character(len=MAX_STRING_LEN) :: input_dir,output_dir,kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: filename, kernel_name, kernel_names(MAX_KERNEL_NAMES)
  character(len=MAX_STRING_LEN) :: arg(NARGS)
  integer :: ier, iker,nker,i,j,k,ispec

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: sem_array

  double precision :: min_val, max_val

  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! parse command line arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i), status=ier)
    if (i <= 1 .and. trim(arg(i)) == '') then
      if (myrank == 0) then
      print *, 'USAGE: mpirun -np NPROC bin/xclip_sem MIN_VAL MAX_VAL KERNEL_NAMES INPUT_FILE OUTPUT_DIR'
      print *, ''
      stop 'Please check command line arguments'
      endif
    endif
  enddo

  read(arg(1),*) min_val
  read(arg(2),*) max_val
  read(arg(3),'(a)') kernel_names_comma_delimited
  read(arg(4),'(a)') input_dir
  read(arg(5),'(a)') output_dir

  ! checks if number of MPI process as specified
  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print*,''
      print*,'Expected number of MPI processes: ', NPROCTOT_VAL
      print*,'Actual number of MPI processes: ', sizeprocs
      print*,''
    endif
    call synchronize_all()
    stop 'Error total number of slices'
  endif
  call synchronize_all()


  ! start execution
  if(myrank==0) then
    write(*,*) 'Running XCLIP_SEM'
    write(*,*)
  endif
  call synchronize_all()

  ! parse kernel names
  call parse_kernel_names(kernel_names_comma_delimited, kernel_names, nker)

  ! initialize array
  allocate(sem_array(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE))

  ! clip kernels
  do iker=1,nker

      kernel_name = trim(kernel_names(iker))
      write(filename,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(kernel_name)//'.bin'

      ! read array
      open(IIN,file=trim(filename),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  file not found: ',trim(filename)
        stop 'File not found'
      endif
      read(IIN) sem_array
      close(IIN)

     if (myrank==0) write(*,*) 'clipping array: ',trim(kernel_names(iker))

     ! apply thresholds
      do ispec=1,NSPEC_CRUST_MANTLE
        do k=1,NGLLZ
          do j=1,NGLLY
            do i=1,NGLLX
              if (sem_array(i,j,k,ispec) < min_val) sem_array(i,j,k,ispec) = min_val
              if (sem_array(i,j,k,ispec) > max_val) sem_array(i,j,k,ispec) = max_val
            enddo
          enddo
        enddo
      enddo

      ! write clipped array
      kernel_name = trim(kernel_names(iker))//'_clip'
      write(filename,'(a,i6.6,a)') trim(input_dir)//'/proc',myrank,'_'//trim(kernel_name)//'.bin'

      open(IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0) then
        write(*,*) '  error opening file: ',trim(filename)
        stop 'Error opening file'
      endif
      write(IOUT) sem_array
      close(IOUT)

  enddo

  if(myrank==0) write(*,*) 'done writing all kernels, see directory ', output_dir

  ! stop all the processes, and exit
  call finalize_mpi()

end program clip_sem_globe

