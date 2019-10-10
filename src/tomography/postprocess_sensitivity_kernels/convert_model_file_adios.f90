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
! the Free Software Foundation; either version 3 of the License, or
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

! Ebru Bozdag and Daniel Peter,
! Nice & Zurich, September 2014
!
! converts between adios and binary format for a model file 'model_gll.bp'

program convert_model_file_adios

  use constants, only: ADIOS_TRANSPORT_METHOD

  use postprocess_par, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IIN,IOUT, &
    MAX_STRING_LEN,NPROCTOT_VAL,NSPEC_CRUST_MANTLE

  use adios_read_mod
  use adios_write_mod
  use adios_helpers_mod
  use manager_adios

  implicit none

  !-------------------------------------------------------------------
  ! USER PARAMETERS
  ! ADIOS model file name
  character(len=64) :: model_adios_file = 'model_gll.bp'

  ! model type
  ! isotropic model parameters (vp,vs,rho) or
  logical :: HAS_ISOTROPY
  ! transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)
  logical :: HAS_TRANSVERSE_ISOTROPY
  ! azimuthal anisotropic (vpv,vph,vsv,vsh,eta,Gs,Gc,rho)
  logical :: HAS_AZIMUTHAL_ANISO
  ! shear attenuation (Qmu)
  logical :: HAS_ATTENUATION_Q

  !-------------------------------------------------------------------

  ! model files
  integer, parameter :: NSPEC = NSPEC_CRUST_MANTLE

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: &
    model_par

  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: model

  ! mpi
  integer :: myrank, sizeprocs

  ! model parameters
  character(len=MAX_STRING_LEN) :: m_file
  character(len=MAX_STRING_LEN) :: m_adios_file

  integer, parameter :: nparams = 12
  character(len=16) :: model_name(nparams),fname(nparams),found_fname(nparams),found_model_name(nparams)

  integer :: i,iker,icount
  logical :: exist,exist_all

  ! input arguments
  character(len=MAX_STRING_LEN) :: arg
  character(len=MAX_STRING_LEN) :: input_model_dir, output_model_dir
  integer :: convert_format

  ! adios
  character(len=MAX_STRING_LEN) :: group_name
  integer(kind=8) :: group, model_handle
  integer(kind=8) :: totalsize
  integer(kind=8) :: group_size_inc
  integer :: local_dim
  integer :: comm
  integer :: ier

  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks number of processes
  ! note: must run as with same number of process as file was created
  if (sizeprocs /= NPROCTOT_VAL) then
    ! usage info
    if (myrank == 0) then
      print *, "this program must be executed in parallel with NPROCTOT_VAL = ",NPROCTOT_VAL,"processes"
      print *, "Invalid number of processes used: ", sizeprocs, " procs"
      print *
      print *, "Please run: mpirun -np ",NPROCTOT_VAL," ./bin/xconvert_model_file_adios .."
    endif
    call abort_mpi()
  endif
  call synchronize_all()

  ! reads input arguments
  do i = 1, 3
    call get_command_argument(i,arg)

    ! usage info
    if (len_trim(arg) == 0) then
      if (myrank == 0) then
        print *, ' '
        print *, ' Usage: xconvert_model_file_adios type MODEL_INPUT_DIR/ MODEL_OUTPUT_DIR/ '
        print *, ' '
        print *, ' with'
        print *, '   type = 1 - convert from ADIOS to (old) binaries'
        print *, '        = 2 - convert from (old) binaries to ADIOS'
        print *, '   MODEL_INPUT_DIR/  - directory which holds input model file(s)'
        print *, '   MODEL_OUTPUT_DIR/ - directory which will hold output model file(s)'
        print *, ' '
      endif
      stop ' Reenter command line options'
    endif

    select case (i)
    case (1)
      read(arg(1:len_trim(arg)),*) convert_format
    case (2)
      input_model_dir = trim(arg)
    case (3)
      output_model_dir = trim(arg)
    end select
  enddo

  ! checks arguments
  if (convert_format /= 1 .and. convert_format /= 2) then
    print *,'Error: invalid format type',convert_format
    stop ' Reenter format type in command line options'
  endif
  call synchronize_all()

  ! sets adios model file name
  if ( convert_format == 1) then
    ! from adios to old binaries
    m_adios_file = trim(input_model_dir) // '/' // trim(model_adios_file)
  else
    ! from old binaries to adios
    m_adios_file = trim(output_model_dir) // '/' // trim(model_adios_file)
  endif

  ! defines model parameters
  ! note: adds space at endings to equal number of characters
  !       to avoid compiler error: "Different shape for array assignment.."
  model_name(1:3) = (/character(len=16) :: "reg1/vp ","reg1/vs ","reg1/rho"/)
  fname(1:3) = (/character(len=16) :: "vp ","vs ","rho"/)
  ! tiso
  model_name(4:8) = (/character(len=16) :: "reg1/vpv","reg1/vph","reg1/vsv","reg1/vsh","reg1/eta"/)
  fname(4:8) = (/character(len=16) :: "vpv","vph","vsv","vsh","eta"/)
  ! azi
  model_name(9:10) = (/character(len=16) :: "reg1/Gc_prime","reg1/Gs_prime"/)
  fname(9:10) = (/character(len=16) :: "Gc_prime","Gs_prime"/)
  model_name(11) = "reg1/mu0"
  fname(11) = "mu0"
  ! adds shear attenuation
  model_name(12) = "reg1/qmu"
  fname(12) = "qmu"

  ! user output
  if (myrank == 0) then
    print *, 'program convert_model_file_adios:'
    print *, ' '
    if ( convert_format == 1) then
      print *, 'conversion type : ',convert_format,' - converts from adios to binary files'
    else
      print *, 'conversion type : ',convert_format,' - converts from binary files to adios'
    endif
    print *, 'input  directory: ',trim(input_model_dir)
    print *, 'output directory: ',trim(output_model_dir)
    print *, ' '
    print *, 'crust/mantle region:'
    print *, '  number of spectral elements = ',NSPEC
    print *, ' '
  endif

  ! initializes ADIOS
  if (myrank == 0) then
    print *, 'initializing ADIOS...'
    print *, ' '
  endif
  call initialize_adios()

  ! initializes model values
  allocate(model(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating model array'
  model(:,:,:,:) = 0.0_CUSTOM_REAL

  HAS_ISOTROPY = .false.
  HAS_TRANSVERSE_ISOTROPY = .false.
  HAS_AZIMUTHAL_ANISO = .false.
  HAS_ATTENUATION_Q = .false.

  ! gets MPI communicator for adios calls
  call world_duplicate(comm)

!--------------------------------------------
  if (convert_format == 1) then ! from adios to old binaries
!--------------------------------------------

! READ THE MODEL IN ADIOS FORMAT

    ! user output
    if (myrank == 0) then
      print *, 'reading in ADIOS model file: ',trim(m_adios_file)
    endif

    ! opens adios file
    call open_file_adios_read(m_adios_file)

    ! debug
    !if (myrank == 0) call show_adios_file_variables(m_adios_file)

    icount = 0
    do iker = 1,nparams
      model(:,:,:,:) = 0.0_CUSTOM_REAL

      ! reads in associated model array
      call read_adios_array_gll_check(myrank,NSPEC,model_name(iker),model,ier)
      ! checks if read was successful
      if (ier /= 0) then
        exist = .false.
      else
        exist = .true.
      endif

      ! makes sure all processes have same flag
      call any_all_l(exist,exist_all)

      if (exist_all) then
        ! set name
        icount = icount + 1
        found_fname(icount) = fname(iker)

        ! checks if vp,vs,rho found
        if (iker == 3 .and. icount == 3) HAS_ISOTROPY = .true.
        ! checks if rho,vpv,vph,vsv,vsh,eta found
        if (iker == 8 .and. icount >= 6) HAS_TRANSVERSE_ISOTROPY = .true.
        ! checks if rho,vpv,vph,vsv,vsh,eta,Gc_prime,Gs_prime,mu0 found
        if (iker == 11 .and. icount >= 9) HAS_AZIMUTHAL_ANISO = .true.
        ! checks if qmu
        if (iker == nparams) HAS_ATTENUATION_Q = .true.

        ! writes out binary file
        write(m_file,'(a,i6.6,a)') trim(output_model_dir)//'/proc',myrank,'_reg1_'//trim(fname(iker))//'.bin'
        open(IOUT,file=trim(m_file),form='unformatted',action='write',iostat=ier)
        if (ier /= 0) then
          print *, 'Error opening binary parameter file: ',trim(m_file)
          stop 'Error opening binary parameter file'
        endif
        ! model values
        write(IOUT) model
        close(IOUT)
      endif

    enddo

    ! closes adios file
    call close_file_adios_read()

    ! check if file written
    if (icount == 0) then
      if (myrank == 0) then
        print *,'Error: no model parameters found in ADIOS file. please check file model_gll.bp in input directory...'
        stop 'No parameters found'
      endif
    endif
    call synchronize_all()

    ! user output
    if (myrank == 0) then
      print *, ' '
      print *, 'written out binary files...'
      print *, ' '
      print *, 'for model parameters:'
      if (HAS_TRANSVERSE_ISOTROPY) &
        print *, '  transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)'
      if (HAS_ISOTROPY) &
        print *, '  isotropic model parameters (vp,vs,rho)'
      if (HAS_AZIMUTHAL_ANISO) &
        print *, '  azimuthal anisotropic model parameters (Gc_prime,Gs_prime)'
      if (HAS_ATTENUATION_Q) &
        print *, '  attenuation model parameter (qmu)'
      print *, ' '
      do iker = 1,icount
        print *, '  found parameter: ',trim(found_fname(iker))
      enddo
      print *, ' '
      print *, 'done writing the model in binary format'
    endif

!--------------------------------------------
  else if (convert_format == 2) then ! from binaries to adios
!--------------------------------------------

! READ THE MODEL IN OLD BINARY FORMAT

    ! user output
    if (myrank == 0) then
      print *, 'reading in binary model files'
    endif

    icount = 0
    do iker = 1,nparams
      ! file name
      write(m_file,'(a,i6.6,a)') trim(input_model_dir)//'/proc',myrank,'_reg1_'//trim(fname(iker))//'.bin'

      inquire(file=trim(m_file),EXIST=exist)
      ! makes sure all processes have same flag
      call any_all_l(exist,exist_all)

      if (exist_all) then
        icount = icount + 1
        found_fname(icount) = fname(iker)
        ! stores adios name
        found_model_name(icount) = model_name(iker)

        ! user output
        if (myrank == 0) then
          print *, '  found parameter: ',trim(fname(iker))
        endif

        ! checks if vp,vs,rho found
        if (iker == 3 .and. icount == 3) HAS_ISOTROPY = .true.
        ! checks if rho,vpv,vph,vsv,vsh,eta found
        if (iker == 8 .and. icount >= 6) HAS_TRANSVERSE_ISOTROPY = .true.
        ! checks if rho,vpv,vph,vsv,vsh,eta,Gc_prime,Gs_prime,mu0 found
        if (iker == 11 .and. icount >= 9) HAS_AZIMUTHAL_ANISO = .true.
        ! checks if qmu
        if (iker == nparams) HAS_ATTENUATION_Q = .true.
      endif
    enddo
    ! checks if any found
    if (icount == 0) then
      print *,'Error: rank ',myrank,' found no binary parameter files. please check your input directory...'
      stop 'No binary files found'
    endif
    call synchronize_all()

    allocate(model_par(NGLLX,NGLLY,NGLLZ,NSPEC,icount),stat=ier)
    if (ier /= 0) stop 'Error allocating model_par array'

    ! reads in parameter values
    do iker = 1,icount
      ! file name
      write(m_file,'(a,i6.6,a)') trim(input_model_dir)//'/proc',myrank,'_reg1_'//trim(found_fname(iker))//'.bin'
      open(IIN,file=trim(m_file),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        print *, 'Error: rank ',myrank,' could not open binary parameter file: ',trim(m_file)
        stop 'Error opening parameter file'
      endif
      read(IIN) model
      close(IIN)

      ! stores parameter
      model_par(:,:,:,:,iker) = model(:,:,:,:)
    enddo

! WRITE OUT THE MODEL IN ADIOS FORMAT

    ! user output
    if (myrank == 0) then
      print *, ' '
      print *, 'writing out ADIOS model file: ',trim(m_adios_file)
      print *, ' '
      if (HAS_TRANSVERSE_ISOTROPY) &
        print *, '  transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)'
      if (HAS_ISOTROPY) &
        print *, '  isotropic model parameters (vp,vs,rho)'
      if (HAS_AZIMUTHAL_ANISO) &
        print *, '  azimuthal anisotropic model parameters (Gc_prime,Gs_prime)'
      if (HAS_ATTENUATION_Q) &
        print *, '  attenuation model parameter (qmu)'
      print *, ' '
    endif

    group_size_inc = 0
    group_name = "MODELS_GROUP"

    call adios_declare_group(group, group_name, '', 1, ier)
    call adios_select_method(group, ADIOS_TRANSPORT_METHOD, '', '', ier)
    call define_adios_scalar(group, group_size_inc, '', "NSPEC", nspec)

    ! Setup ADIOS for the current group
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC

    ! Define ADIOS Variables
    do iker=1,icount
      call define_adios_global_array1D(group, group_size_inc,local_dim,'',trim(found_model_name(iker)),model)
    enddo

    ! Open an handler to the ADIOS file and setup the group size
    call adios_open(model_handle, group_name, trim(m_adios_file), "w", comm, ier);
    if (ier /= 0) then
      print *, 'Error opening adios model file: ',trim(m_adios_file)
      stop 'Error opening adios model file'
    endif

    call adios_group_size (model_handle, group_size_inc, totalsize, ier)
    if (ier /= 0 ) stop 'Error calling adios_group_size() routine failed'

    call adios_write(model_handle, "NSPEC", nspec, ier)

    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC
    do iker=1,icount
      model(:,:,:,:) = model_par(:,:,:,:,iker)
      ! Write previously defined ADIOS variables
      call write_adios_global_1d_array(model_handle, myrank, sizeprocs, local_dim, &
                                       trim(found_model_name(iker)),model(:,:,:,:))
    enddo
    ! Perform the actual write to disk
    call adios_set_path(model_handle, '', ier)
    call adios_close(model_handle, ier)

    ! free memory
    deallocate(model_par)

    if (myrank == 0) print *, 'done writing the model in adios format'
  endif

  ! free memory
  deallocate(model)

  ! user output
  if (myrank == 0) then
    print *, ' '
    print *, 'see output file(s) in directory: ',trim(output_model_dir)
    print *, ' '
  endif

  call adios_finalize (myrank, ier)
  call finalize_mpi()

end program convert_model_file_adios
