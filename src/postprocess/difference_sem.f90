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

!------------------------------------------------------------------------------------------------
! difference_sem
!
! this runs in parallel, please submit as parallel job.
! takes the difference between proc***.bin from two different input directories
!
! usage: mpirun -np * ./xdifference_sem filename INPUT_dir_1/ INPUT_dir_2/ OUTPUT_DIR/ [region]
!
!------------------------------------------------------------------------------------------------

program difference_sem

  use postprocess_par,only: &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IIN,IOUT,MAX_STRING_LEN, &
    NPROCTOT_VAL,NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE

  implicit none

  character(len=MAX_STRING_LEN) :: arg(5)
  character(len=MAX_STRING_LEN) :: file1name,file2name
  character(len=MAX_STRING_LEN) :: input1dir,input2dir
  character(len=MAX_STRING_LEN) :: outputdir,kernel_name
  character(len=20) :: reg_name

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: sem_data,sem_data_2

  real(kind=CUSTOM_REAL) :: min,max,min_rel,max_rel
  real(kind=CUSTOM_REAL) :: min_all,max_all,min_rel_all,max_rel_all

  integer :: i,iproc,ier
  integer :: iregion,ir,irs,ire
  integer :: nspec

  ! mpi
  integer :: sizeprocs,myrank

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks compilation setup
  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print*, 'Error number of processors supposed to run on : ',NPROCTOT_VAL
      print*, 'Error number of MPI processors actually run on: ',sizeprocs
      print*
      print*, 'please rerun with: mpirun -np ',NPROCTOT_VAL,' bin/xdifference_sem .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! checks arguments
  do i = 1, 5
    call get_command_argument(i,arg(i))

    ! usage info
    if (i <= 4 .and. trim(arg(i)) == '') then
      if (myrank == 0) then
        print *, ' '
        print *, ' Usage: difference_sem kernel_name input1_dir/ input2_dir/ output_dir/ [region]'
        print *, ' '
        print *, ' with'
        print *, '   kernel_name   - takes files with this kernel name'
        print *, '                     e.g. "vsv" for proc***_reg1_vsv.bin'
        print *, '   input1_dir/   - input directory for first files'
        print *, '   input2_dir/   - input directory for second files'
        print *, '   output_dir/   - output directory for (first - second) file values'
        print *, '   [region]      - optional: if region (1/2/3) is not specified, all 3 regions will be taken,'
        print *, '                             otherwise, only takes region specified'
        print *, ' '
        print *, ' possible kernel_names are: '
        print *, '   "alpha_kernel", "beta_kernel", .., "vsv", "rho_vp", "kappastore", "mustore", etc.'
        print *
        print *, '   that are stored in the local directories input1_dir/ and input2_dir/ '
        print *, '   as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) in proc***_reg1_filename.bin'
        print *, ' '
      endif
      call synchronize_all()
      call exit_MPI(myrank,'Reenter command line options')
    endif
  enddo

  ! get region id
  if (len_trim(arg(5)) == 0) then
    iregion = 0
  else
    read(arg(5),*) iregion
  endif
  if (iregion > 3 .or. iregion < 0) stop 'Error region: region id must be either = 0,1,2,3'
  if (iregion == 0) then
    irs = 1
    ire = 3
  else
    irs = iregion
    ire = irs
  endif

  ! gets kernel and directory names from argument call
  kernel_name = trim(arg(1))
  input1dir = trim(arg(2))
  input2dir = trim(arg(3))
  outputdir = trim(arg(4))

  ! loops over slices
  min_all = 0.0
  max_all = 0.0
  min_rel_all = 0.0
  max_rel_all = 0.0


  ! user output
  if (myrank == 0) then
    write(*,*) 'differencing files: ',sizeprocs,' slices'
    write(*,*)
    write(*,*) 'kernel name: ',trim(kernel_name)
    write(*,*) 'input 1 directory: ',trim(input1dir)
    write(*,*) 'input 2 directory: ',trim(input2dir)
    write(*,*) 'output directory : ',trim(outputdir)
    write(*,*) 'regions: start =', irs, ' to end =', ire
    write(*,*)
  endif

  do ir = irs, ire
    if (myrank == 0) write(*,*) '----------- Region ', ir, '----------------'

    ! initializes
    sem_data(:,:,:,:) = 0.0_CUSTOM_REAL
    sem_data_2(:,:,:,:) = 0.0_CUSTOM_REAL

    ! sets number of elements
    select case (ir)
    case (1)
      ! crust/mantle
      nspec = NSPEC_CRUST_MANTLE
    case (2)
      ! outer core
      nspec = NSPEC_OUTER_CORE
    case (3)
      ! inner core
      nspec = NSPEC_INNER_CORE
    case default
      stop 'Error region id not recognized'
    end select

    ! prefix for (crust/mantle) region
    write(reg_name,'(a,i1,a)') 'reg',ir,'_'

    ! process file name
    iproc = myrank
    write(file1name,'(a,i6.6,a)') trim(input1dir)//'/proc',iproc,'_'//trim(reg_name)//trim(kernel_name)//'.bin'
    write(file2name,'(a,i6.6,a)') trim(input2dir)//'/proc',iproc,'_'//trim(reg_name)//trim(kernel_name)//'.bin'

    ! reads in file from first directory
    if(myrank==0) write(*,*) '  data_1: ',trim(file1name)
    open(IIN,file=trim(file1name),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(file1name)
      stop 'Error opening first data file'
    endif
    read(IIN) sem_data(:,:,:,1:nspec)
    close(IIN)

    ! reads in file from second directory
    if (myrank == 0) write(*,*) '  data_2: ',trim(file2name)
    open(IIN,file=trim(file2name),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(file2name)
      stop 'Error opening second data file'
    endif
    read(IIN) sem_data_2(:,:,:,1:nspec)
    close(IIN)

    ! user output
    if (myrank == 0) then
      write(*,*) '  min/max data_1 value: ',minval(sem_data(:,:,:,1:nspec)),maxval(sem_data(:,:,:,1:nspec))
      write(*,*) '  min/max data_2 value: ',minval(sem_data_2(:,:,:,1:nspec)),maxval(sem_data_2(:,:,:,1:nspec))
      write(*,*)
    endif

    ! stores difference between kernel files
    if (myrank == 0) write(*,*) '  difference: (data_1 - data_2)'

    ! absolute values
    write(file1name,'(a,i6.6,a)') trim(outputdir)//'/proc',iproc,'_'//trim(reg_name)//trim(kernel_name)//'_diff.bin'
    if (myrank == 0) write(*,*) '  file: ',trim(file1name)
    open(IOUT,file=trim(file1name),form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(file1name)
      stop 'Error opening output data file'
    endif

    ! takes the difference
    sem_data(:,:,:,1:nspec) = sem_data(:,:,:,1:nspec) - sem_data_2(:,:,:,1:nspec)

    write(IOUT) sem_data(:,:,:,1:nspec)
    close(IOUT)

    ! min/max
    min = minval(sem_data(:,:,:,1:nspec))
    max = maxval(sem_data(:,:,:,1:nspec))
    call min_all_cr(min,min_all)
    call max_all_cr(max,max_all)

    ! stores relative difference (k1 - k2)/ k2 with respect to second input file
    write(file1name,'(a,i6.6,a)') trim(outputdir)//'/proc',iproc,'_'//trim(reg_name)//trim(kernel_name)//'_diff_relative.bin'
    if (myrank == 0) write(*,*) '  file: ',trim(file1name)
    open(IOUT,file=trim(file1name),form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(file1name)
      stop 'Error opening output data file'
    endif

    ! relative difference (k1 - k2)/ k2 with respect to second input file
    where( sem_data_2(:,:,:,1:nspec) /= 0.0_CUSTOM_REAL)
      sem_data(:,:,:,1:nspec) = sem_data(:,:,:,1:nspec) / sem_data_2(:,:,:,1:nspec)
    elsewhere
      sem_data(:,:,:,1:nspec) = 0.0_CUSTOM_REAL
    endwhere

    write(IOUT) sem_data(:,:,:,1:nspec)
    close(IOUT)

    ! min/max of relative values
    min_rel = minval(sem_data(:,:,:,1:nspec))
    max_rel = maxval(sem_data(:,:,:,1:nspec))
    call min_all_cr(min_rel,min_rel_all)
    call max_all_cr(max_rel,max_rel_all)

    ! output
    if (myrank == 0) then
      write(*,*) '  min/max value         : ',min,max
      write(*,*) '  min/max relative value: ',min_rel,max_rel
      write(*,*)
    endif

  enddo ! region

  ! user output
  if (myrank == 0) then
    write(*,*)
    write(*,*) 'statistics:'
    write(*,*) '  total min/max         : ',min_all,max_all
    write(*,*) '  total relative min/max: ',min_rel_all,max_rel_all
    write(*,*)
    write(*,*) 'done writing all difference and relative difference files'
    write(*,*) 'see output directory: ',trim(outputdir)
    write(*,*)
  endif

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program difference_sem


