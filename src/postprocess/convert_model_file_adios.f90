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

! Ebru & Daniel
! Nice & Zurich, September 2014
!
! converts between adios and binary format for a model file 'model_gll.bp'

program convert_model_file_adios

  use adios_read_mod
  use adios_write_mod
  use adios_helpers_mod

  use postprocess_par,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IIN,IOUT, &
    MAX_STRING_LEN,NPROCTOT_VAL,NSPEC_CRUST_MANTLE

  use constants,only: ADIOS_TRANSPORT_METHOD

  implicit none

  !-------------------------------------------------------------------
  ! USER PARAMETERS
  ! ADIOS model file name
  character(len=64) :: model_adios_file = 'model_gll.bp'

  ! model type
  ! isotropic model parameters (vp,vs,rho) or
  ! transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)
  ! defaults: TI models
  logical,parameter :: USE_TRANSVERSE_ISOTROPY = .true.

  ! shear attenuation
  logical,parameter :: USE_ATTENUATION_Q = .false.

  !-------------------------------------------------------------------

  ! model files
  integer, parameter :: NSPEC = NSPEC_CRUST_MANTLE

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    model_vpv,model_vph,model_vsv,model_vsh,model_eta,model_rho,model_qmu

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: model

  ! mpi
  integer :: myrank, sizeprocs

  ! model parameters
  character(len=MAX_STRING_LEN) :: m_file
  character(len=MAX_STRING_LEN) :: m_adios_file
  character(len=16) :: model_name(7),fname(7)
  integer :: nparams
  integer :: i,iker

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
  integer(kind=8) :: sel
  integer(kind=8) :: sel_start(1),count_ad(1)
  integer :: comm
  integer :: ier

  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)
  call world_get_comm(comm)

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
    print*,'Error: invalid format type',convert_format
    stop ' Reenter format type in command line options'
  endif

  ! sets adios model file name
  if ( convert_format == 1) then
    ! from adios to old binaries
    m_adios_file = trim(input_model_dir) // '/' // trim(model_adios_file)
  else
    ! from old binaries to adios
    m_adios_file = trim(output_model_dir) // '/' // trim(model_adios_file)
  endif

  ! defines model parameters
  if (USE_TRANSVERSE_ISOTROPY) then
    ! transversly isotropic (TI) model
    nparams = 6
    model_name(1:6) = (/character(len=16) :: "reg1/vpv","reg1/vph","reg1/vsv","reg1/vsh","reg1/eta","reg1/rho"/)
    fname(1:6) = (/character(len=16) :: "vpv","vph","vsv","vsh","eta","rho"/)
  else
    ! isotropic model
    nparams = 3
    ! note: adds space at endings to equal number of characters
    !       to avoid compiler error: "Different shape for array assignment.."
    model_name(1:3) = (/character(len=16) :: "reg1/vp ","reg1/vs ","reg1/rho"/)
    fname(1:3) = (/character(len=16) :: "vp ","vs ","rho"/)
  endif
  ! adds shear attenuation
  if (USE_ATTENUATION_Q) then
    nparams = nparams + 1
    model_name(nparams) = "reg1/qmu"
    fname(nparams) = "qmu"
  endif

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
    if (USE_TRANSVERSE_ISOTROPY) then
      print *, 'model parameters:',nparams,' - transversely isotropic model'
    else
      print *, 'model parameters:',nparams,' - isotropic model'
    endif
    if (USE_ATTENUATION_Q) then
      print *, '  includes qmu model parameter'
    endif
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
  call adios_setup()

  ! initializes model values
  model(:,:,:,:) = 0.0_CUSTOM_REAL
  model_vpv(:,:,:,:) = 0.0_CUSTOM_REAL
  model_vph(:,:,:,:) = 0.0_CUSTOM_REAL
  model_vsv(:,:,:,:) = 0.0_CUSTOM_REAL
  model_vsh(:,:,:,:) = 0.0_CUSTOM_REAL
  model_eta(:,:,:,:) = 0.0_CUSTOM_REAL
  model_rho(:,:,:,:) = 0.0_CUSTOM_REAL
  model_qmu(:,:,:,:) = 0.0_CUSTOM_REAL

!--------------------------------------------
  if (convert_format == 1) then ! from adios to old binaries
!--------------------------------------------

! READ THE MODEL IN ADIOS FORMAT

    ! user output
    if (myrank == 0) then
      print *, 'reading in ADIOS model file: ',trim(m_adios_file)
    endif

    call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=1", ier)
    if (ier /= 0 ) stop 'Error in adios_read_init_method()'

    call adios_read_open_file (model_handle, trim(m_adios_file), 0, comm, ier)
    if (ier /= 0) then
      print *, 'Error opening adios model file: ',trim(m_adios_file)
      stop 'Error opening adios model file'
    endif

    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC

    do iker = 1,nparams
      model(:,:,:,:) = 0.0_CUSTOM_REAL

      call adios_get_scalar(model_handle, trim(model_name(iker))//"/local_dim",local_dim, ier)
      if (ier /= 0 ) stop 'Error adios get scalar'

      sel_start(1) = local_dim * myrank
      count_ad(1) = NGLLX * NGLLY * NGLLZ * NSPEC

      call adios_selection_boundingbox(sel, 1, sel_start, count_ad)
      call adios_schedule_read(model_handle, sel,trim(model_name(iker))//"/array",0, 1, model, ier)
      if (ier /= 0 ) stop 'Error adios schedule read'

      call adios_perform_reads(model_handle, ier)
      if (ier /= 0 ) stop 'Error adios perform read'

      if (USE_TRANSVERSE_ISOTROPY) then
        ! TI model
        select case (iker)
        case (1)
           model_vpv(:,:,:,:) = model(:,:,:,:)
        case (2)
           model_vph(:,:,:,:) = model(:,:,:,:)
        case (3)
           model_vsv(:,:,:,:) = model(:,:,:,:)
        case (4)
           model_vsh(:,:,:,:) = model(:,:,:,:)
        case (5)
           model_eta(:,:,:,:) = model(:,:,:,:)
        case (6)
           model_rho(:,:,:,:) = model(:,:,:,:)
        end select
      else
        ! isotropic model
        select case (iker)
        case (1)
           model_vpv(:,:,:,:) = model(:,:,:,:)
        case (2)
           model_vsv(:,:,:,:) = model(:,:,:,:)
        case (3)
           model_rho(:,:,:,:) = model(:,:,:,:)
        end select
      endif
      ! adds shear attenuation
      if (USE_ATTENUATION_Q) then
        if (iker == nparams) then
           model_qmu(:,:,:,:) = model(:,:,:,:)
        endif
      endif
    enddo

    call adios_read_close(model_handle,ier)
    if (ier /= 0 ) stop 'Error adios read close'

    call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
    if (ier /= 0 ) stop 'Error adios read finalize'

    ! WRITE OUT THE MODEL IN OLD BINARIES

    ! user output
    if (myrank == 0) then
      print *, ' '
      print *, 'writing out binary files...'
    endif

    do iker = 1,nparams
      ! user output
      if (myrank == 0) then
        print *, '  for parameter: ',trim(fname(iker))
      endif

      ! file name
      write(m_file,'(a,i6.6,a)') trim(output_model_dir)//'/proc',myrank,'_reg1_'//trim(fname(iker))//'.bin'
      open(IOUT,file=trim(m_file),form='unformatted',action='write',iostat=ier)
      if (ier /= 0) then
        print *, 'Error opening binary parameter file: ',trim(m_file)
        stop 'Error opening binary parameter file'
      endif

      ! selects output
      if (USE_TRANSVERSE_ISOTROPY) then
        ! TI model
        select case (iker)
        case (1)
          ! vpv model
          write(IOUT) model_vpv
        case (2)
          ! vph model
          write(IOUT) model_vph
        case (3)
          ! vsv model
          write(IOUT) model_vsv
        case (4)
          ! vsh model
          write(IOUT) model_vsh
        case (5)
          ! eta model
          write(IOUT) model_eta
        case (6)
          ! rho model
          write(IOUT) model_rho
        end select
      else
        ! isotropic model
        select case (iker)
        case (1)
          ! vp model
          write(IOUT) model_vpv
        case (2)
          ! vs model
          write(IOUT) model_vsv
        case (3)
          ! rho model
          write(IOUT) model_rho
        end select
      endif
      ! adds shear attenuation
      if (USE_ATTENUATION_Q) then
        if (iker == nparams) then
          ! qmu model
          write(IOUT) model_qmu
        endif
      endif

      close(IOUT)
    enddo

    if (myrank==0) print *, 'done writing the model in binary format'

!--------------------------------------------
  else if (convert_format == 2) then ! from binaries to adios
!--------------------------------------------

! READ THE MODEL IN OLD BINARY FORMAT

    ! user output
    if (myrank == 0) then
      print *, 'reading in binary model files'
    endif

    do iker = 1,nparams
      ! user output
      if (myrank == 0) then
        print *, '  for parameter: ',trim(fname(iker))
      endif

      ! file name
      write(m_file,'(a,i6.6,a)') trim(input_model_dir)//'/proc',myrank,'_reg1_'//trim(fname(iker))//'.bin'
      open(IIN,file=trim(m_file),status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0) then
        print *, 'Error opening binary parameter file: ',trim(m_file)
        stop 'Error opening binary parameter file'
      endif

      ! selects output
      if (USE_TRANSVERSE_ISOTROPY) then
        ! TI model
        select case (iker)
        case (1)
          ! vpv model
          read(IIN) model_vpv
        case (2)
          ! vph model
          read(IIN) model_vph
        case (3)
          ! vsv model
          read(IIN) model_vsv
        case (4)
          ! vsh model
          read(IIN) model_vsh
        case (5)
          ! eta model
          read(IIN) model_eta
        case (6)
          ! rho model
          read(IIN) model_rho
        end select
      else
        ! isotropic model
        select case (iker)
        case (1)
          ! vp model
          read(IIN) model_vpv
        case (2)
          ! vs model
          read(IIN) model_vsv
        case (3)
          ! rho model
          read(IIN) model_rho
        end select
      endif
      ! adds shear attenuation
      if (USE_ATTENUATION_Q) then
        if (iker == nparams) then
          ! qmu model
          read(IIN) model_qmu
        endif
      endif

      close(IIN)
    enddo

! WRITE OUT THE MODEL IN ADIOS FORMAT

    ! user output
    if (myrank == 0) then
      print *, ' '
      print *, 'writing out ADIOS model file: ',trim(m_adios_file)
    endif

    group_size_inc = 0
    group_name = "MODELS_GROUP"

    call adios_declare_group(group, group_name, "", 1, ier)
    call adios_select_method(group, ADIOS_TRANSPORT_METHOD, "", "", ier)
    call define_adios_scalar(group, group_size_inc, "", "NSPEC", nspec)

    ! Setup ADIOS for the current group
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC

    ! Define ADIOS Variables
    do iker=1,nparams
      call define_adios_global_array1D(group, group_size_inc,local_dim,"",trim(model_name(iker)),model)
    enddo

    ! Open an handler to the ADIOS file and setup the group size
    call adios_open(model_handle, group_name, trim(m_adios_file), "w", comm, ier);
    if (ier /= 0) then
      print *, 'Error opening adios model file: ',trim(m_adios_file)
      stop 'Error opening adios model file'
    endif

    call adios_group_size (model_handle, group_size_inc, totalsize, ier)
    call adios_write(model_handle, "NSPEC", nspec, ier)

    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC
    do iker=1,nparams
      if (USE_TRANSVERSE_ISOTROPY) then
        ! TI model
        select case (iker)
        case (1)
          model(:,:,:,:) = model_vpv(:,:,:,:)
        case (2)
          model(:,:,:,:) = model_vph(:,:,:,:)
        case (3)
          model(:,:,:,:) = model_vsv(:,:,:,:)
        case (4)
          model(:,:,:,:) = model_vsh(:,:,:,:)
        case (5)
          model(:,:,:,:) = model_eta(:,:,:,:)
        case (6)
          model(:,:,:,:) = model_rho(:,:,:,:)
        end select
      else
        ! isotropic model
        select case (iker)
        case (1)
          model(:,:,:,:) = model_vpv(:,:,:,:)
        case (2)
          model(:,:,:,:) = model_vsv(:,:,:,:)
        case (3)
          model(:,:,:,:) = model_rho(:,:,:,:)
        end select
      endif
      ! adds shear attenuation
      if (USE_ATTENUATION_Q) then
        if (iker == nparams) then
          ! qmu model
          model(:,:,:,:) = model_qmu(:,:,:,:)
        endif
      endif

      ! Write previously defined ADIOS variables
      call write_adios_global_1d_array(model_handle, myrank, sizeprocs, local_dim, &
                                       trim(model_name(iker)),model(:,:,:,:))
    enddo
    ! Perform the actual write to disk
    call adios_set_path(model_handle, "", ier)
    call adios_close(model_handle, ier)

    if (myrank==0) print *, 'done writing the model in adios format'
  endif
  ! user output
  if (myrank==0) then
    print *, ' '
    print *, 'see output file(s) in directory: ',trim(output_model_dir)
    print *, ' '
  endif

  call adios_finalize (myrank, ier)
  call finalize_mpi()

end program convert_model_file_adios
