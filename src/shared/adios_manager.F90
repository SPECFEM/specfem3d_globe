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

!> Tools to setup and cleanup ADIOS,
!> contains wrapper subroutines for common adios calls
!
! note: adios library calls use a format like "adios_do_something***()"
!       our own wrapper functions thus will rather use something like "do_something_adios***()"
!       to better distinguish between library functions and wrappers.

module manager_adios

  implicit none

  private

  ! MPI copies of communicator and rank
  integer :: comm_adios
  integer :: myrank_adios

  integer,public :: sizeprocs_adios

  ! initialized flag
  logical :: is_adios_initialized

#ifdef ADIOS_INPUT
  ! adios error message
  character(len=1024) :: err_message
  character(len=*),parameter :: ADIOS_VERBOSITY = "verbose=1" ! lowest level: verbose=1
#endif

  ! debugging
  logical,parameter :: DEBUG = .false.

  ! public accessibility

  ! file handle for read/write
  integer(kind=8),public :: file_handle_adios

  ! public routines
  public :: initialize_adios
  public :: finalize_adios
  public :: close_file_adios

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT
  public :: close_file_adios_read
  public :: init_adios_group
  public :: set_adios_group_size
  ! reading
  public :: read_adios_array_gll
  public :: read_adios_array_gll_int
  public :: read_adios_array_1d
  public :: read_adios_array_1d_int
  public :: read_adios_scalar_int
  public :: read_adios_scalar_int_only_rank
  ! writing
  public :: write_adios_array_gll
  public :: write_adios_scalar_int
  ! file opening
  public :: open_file_adios_read
  public :: open_file_adios_read_only_rank
  public :: open_file_adios_write
  public :: open_file_adios_write_append
#endif

contains


!-------------------------------------------------------------------------------
!
! public ADIOS wrapper routines (also available without adios compilation support)
!
!-------------------------------------------------------------------------------

  subroutine initialize_adios()

!> Initialize ADIOS and setup the xml output file

  use constants, only: ADIOS_BUFFER_SIZE_IN_MB

#ifdef ADIOS_INPUT
  use adios_write_mod, only: adios_init_noxml,adios_set_max_buffer_size
#endif

  implicit none

  ! local parameters
#ifdef ADIOS_INPUT
  integer :: ier
#endif

  ! initializes
  is_adios_initialized = .false.
  comm_adios = 0
  myrank_adios = -1
  file_handle_adios = 0

#ifdef ADIOS_INPUT

  ! gets MPI communicator for adios calls
  call world_duplicate(comm_adios)

  ! gets rank from (duplicate) adios communicator
  call world_rank_comm(myrank_adios,comm_adios)

  ! number of MPI processes
  call world_size(sizeprocs_adios)

  call adios_init_noxml (comm_adios, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !       e.g., version 1.5.0 returns 1 here
  !print *,'adios init return: ',ier
  !if (ier /= 0 ) stop 'Error setting up ADIOS: calling adios_init_noxml() routine failed'

! ask/check at configuration step for adios version 1.10 or higher?
#ifdef ADIOS_VERSION_OLD
  ! ADIOS versions <= 1.9
  ! note: for newer versions ( >= 1.10), this will produce a warning might not be supported anymore
  call adios_allocate_buffer(ADIOS_BUFFER_SIZE_IN_MB, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !       e.g., version 1.5.0 returns 1 if called first time, 0 if already called
  !print *,'adios allocate buffer return: ',ier
  !call check_adios_err(myrank_adios,ier)
#else
  ! ADIOS versions >= 1.10
  call adios_set_max_buffer_size(ADIOS_BUFFER_SIZE_IN_MB)
#endif

  ! sets flag
  is_adios_initialized = .true.

#else

  ! gets rank
  call world_rank(myrank_adios)

  ! compilation without ADIOS support
  if (myrank_adios == 0) then
    print *, "Error: ADIOS enabled without ADIOS Support."
    print *, "To enable ADIOS support, reconfigure with --with-adios flag."
  endif
  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: intitialize called without compilation support")

#endif

  end subroutine initialize_adios

!
!-------------------------------------------------------------------------------
!

  subroutine finalize_adios()

!> Finalize ADIOS. Must be called once everything is written down.

#ifdef ADIOS_INPUT
  use adios_write_mod, only: adios_finalize
#endif

  implicit none

  ! local parameters
#ifdef ADIOS_INPUT
  integer :: ier
#endif

  ! synchronizes all first
  call synchronize_all_comm(comm_adios)

#ifdef ADIOS_INPUT
  ! finalize
  call adios_finalize(myrank_adios, ier)
  if (ier /= 0 ) stop 'Error cleaning up ADIOS: calling adios_finalize() routine failed'

  ! frees (duplicate) MPI communicator
  call world_comm_free(comm_adios)

#else
  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: finalize called without compilation support")
#endif

  end subroutine finalize_adios

!
!-------------------------------------------------------------------------------
!


  subroutine close_file_adios()

#ifdef ADIOS_INPUT
  use adios_write_mod, only: adios_close
  use adios_read_mod, only: adios_errmsg
#endif

  implicit none

  ! Variables
#ifdef ADIOS_INPUT
  integer :: ier

  ! closes file
  call adios_close(file_handle_adios, ier)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error closing adios file'
    print *,trim(err_message)
    stop 'Error closing ADIOS file: calling adios_close() routine failed'
  endif

  ! sets explicitly to zero
  file_handle_adios = 0

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

#else

  ! safety stop
  call exit_MPI(myrank_adios,"Error ADIOS manager: close file called without compilation support")

#endif

  end subroutine close_file_adios


!-------------------------------------------------------------------------------
!
! ADIOS wrapper routines (only available with adios compilation support)
!
!-------------------------------------------------------------------------------


! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine close_file_adios_read()

  use adios_read_mod, only: adios_read_close,adios_read_finalize_method,ADIOS_READ_METHOD_BP

  implicit none

  ! local parameters
  integer :: ier

  call adios_read_close(file_handle_adios,ier)
  if (ier /= 0 ) stop 'Error helper adios read close'

  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
  if (ier /= 0 ) stop 'Error helper adios read finalize'

  end subroutine close_file_adios_read

#endif

!
!---------------------------------------------------------------------------------
!


! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine init_adios_group(adios_group,group_name)

! useful to read in data from the same number of processors
! as the data was written from

  use constants, only: ADIOS_TRANSPORT_METHOD

  use adios_write_mod, only: adios_declare_group,adios_select_method

  implicit none

  integer(kind=8),intent(inout) :: adios_group
  character(len=*),intent(in) :: group_name

  ! local parameters
  integer :: ier

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios group init: please initialize adios first using intitialize_adios() routine'

  ! debug
  if (DEBUG) print *,'***debug ADIOS: rank ',myrank_adios,' init group ',trim(group_name),'****'

  ! initializes adios group
  call adios_declare_group(adios_group, group_name, '', 0, ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,ier)

  ! We set the transport method to 'MPI'. This seems to be the correct choice
  ! for now. We might want to move this to the constant.h file later on.
  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, '', '', ier)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,ier)

  end subroutine init_adios_group

#endif

!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine set_adios_group_size(groupsize)

  use adios_write_mod, only: adios_group_size

  implicit none

  integer(kind=8),intent(in) :: groupsize

  ! local parameters
  integer :: ier
  integer(kind=8) :: totalsize

  ! checks if file handle valid
  if (file_handle_adios == 0) stop 'Invalid ADIOS file handle in set_adios_group_size()'

  call adios_group_size(file_handle_adios, groupsize, totalsize, ier)
  if (ier /= 0 ) stop 'Error calling adios_group_size() routine failed'

  end subroutine set_adios_group_size

#endif


!
!---------------------------------------------------------------------------------
!


! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine open_file_adios_read(filename)

! useful to read in data from the same number of processors
! as the data was written from

  use adios_read_mod, only: adios_read_open_file,adios_read_init_method,ADIOS_READ_METHOD_BP,adios_errmsg

  implicit none

  character(len=*),intent(in) :: filename

  ! local parameters
  integer :: ier

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios read: please initialize adios first using intitialize_adios() routine'

  ! initializes read method
  call adios_read_init_method(ADIOS_READ_METHOD_BP, comm_adios, ADIOS_VERBOSITY, ier)
  if (ier /= 0 ) then
    call adios_errmsg(err_message)
    print *,'Error initializing read adios for file: ',trim(filename)
    print *,trim(err_message)
    stop 'Error initializing adios read method'
  endif

  ! debug
  if (DEBUG) print *,'***debug ADIOS: rank ',myrank_adios,' read open file ',trim(filename),'****'

  ! opens file
  call adios_read_open_file(file_handle_adios, trim(filename), 0, comm_adios, ier)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error opening adios file for reading: ',trim(filename)
    print *,trim(err_message)
    stop 'Error opening adios file'
  endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_read

#endif

!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine open_file_adios_read_only_rank(rank,filename)

! only single process is reading, useful for file inquiry

  use adios_read_mod, only: adios_read_open_file,adios_read_init_method,ADIOS_READ_METHOD_BP, &
    adios_inq_file,adios_inq_varnames,adios_inq_attrnames

  implicit none

  integer,intent(in) :: rank
  character(len=*),intent(in) :: filename

  ! local parameters
  integer :: ier
  integer :: comm_dummy
  character(len=128) :: name !uses a string copy, trying to prevent a memory corruption issue somewhere in adios...

  ! only specified rank proceeds
  if (myrank_adios /= rank) return

  ! gets MPI communicator for only a single process
  call world_get_comm_self(comm_dummy)

  ! copies name
  if (len_trim(filename) > 128) then
    stop 'Error filename provided too long'
  else
    name = trim(filename)
  endif

  ! initializes read method
  call adios_read_init_method(ADIOS_READ_METHOD_BP, comm_dummy, ADIOS_VERBOSITY, ier)
  if (ier /= 0 ) then
    call adios_errmsg(err_message)
    print *, 'Error initializing read adios by master for file: ',trim(name)
    print *,trim(err_message)
    stop 'Error initializing adios read by master method'
  endif

  ! debug
  if (DEBUG) print *,'***debug ADIOS: only rank ',myrank_adios,' read open file ',trim(filename),'****'

  ! opens file
  call adios_read_open_file(file_handle_adios, trim(name), 0, comm_dummy, ier)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *, 'Error opening adios file for reading: ',trim(name)
    print *,trim(err_message)
    stop 'Error opening adios file'
  endif

  ! shows file contents
  call show_adios_file_variables(name)

  end subroutine open_file_adios_read_only_rank

#endif

!
!-------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine open_file_adios_write(filename,group_name)

! opens adios file for writing

  use adios_write_mod, only: adios_open

  implicit none

  character(len=*),intent(in) :: filename
  character(len=*) :: group_name

  ! local parameters
  integer :: ier

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios write: please initialize adios first using intitialize_adios() routine'

  ! debug
  if (DEBUG) print *,'***debug ADIOS: rank ',myrank_adios,' open file ',trim(filename),' (for writing) ****'

  ! opens file
  call adios_open(file_handle_adios, group_name, trim(filename), "w", comm_adios, ier)

  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error opening adios file for writing: ',trim(filename)
    print *,'Error rank',myrank_adios,'- group name ',trim(group_name)
    print *,trim(err_message)
    stop 'Error opening adios file'
  endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_write

#endif


!
!-------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine open_file_adios_write_append(filename,group_name)

! open adios file for appending data

  use adios_write_mod, only: adios_open

  implicit none

  character(len=*),intent(in) :: filename
  character(len=*) :: group_name

  ! local parameters
  integer :: ier

  ! check
  if (.not. is_adios_initialized) &
    stop 'Error intializing for adios write: please initialize adios first using intitialize_adios() routine'

  ! debug
  if (DEBUG) print *,'***debug ADIOS: rank ',myrank_adios,' open file ',trim(filename),' (for appending) ****'

  ! opens file
  call adios_open(file_handle_adios, group_name, trim(filename), "a", comm_adios, ier)

  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error opening adios file for appending: ',trim(filename)
    print *,'Error rank',myrank_adios,'- group name ',trim(group_name)
    print *,trim(err_message)
    stop 'Error opening adios file'
  endif

  ! synchronizes all processes
  call synchronize_all_comm(comm_adios)

  end subroutine open_file_adios_write_append

#endif


!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine read_adios_scalar_int(rank,scalar_name,scalar)

! reads in a single integer value

  use adios_read_mod, only: adios_schedule_read,adios_perform_reads

  implicit none

  integer, intent(in) :: rank

  integer, intent(out) :: scalar
  character(len=*),intent(in) :: scalar_name

  ! local parameters
  integer :: ier
  integer(kind=8) :: sel

  ! selects data block
  call adios_selection_writeblock(sel, rank)

  ! reads array
  call adios_schedule_read(file_handle_adios, sel, trim(scalar_name), 0, 1, scalar, ier)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error adios: could not read parameter: ',trim(scalar_name)
    print *,trim(err_message)
    stop 'Error adios helper read scalar'
  endif

  call adios_perform_reads(file_handle_adios, ier)
  if (ier /= 0 ) stop 'Error helper adios read scalar failed'

  end subroutine read_adios_scalar_int

#endif

!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine read_adios_scalar_int_only_rank(rank,scalar_name,scalar)

  use adios_read_mod, only: adios_get_scalar

  implicit none

  integer, intent(in) :: rank

  integer, intent(out) :: scalar
  character(len=*),intent(in) :: scalar_name

  ! local parameters
  integer :: ier
  !integer(kind=8) :: sel
  character(len=128) :: name !uses a string copy, trying to prevent a memory corruption issue somewhere in adios...

  ! checks if anything to do
  if (myrank_adios /= rank) return

  ! copies name
  if (len_trim(scalar_name) > 128) then
    stop 'Error scalar name provided too long'
  else
    name = trim(scalar_name)
  endif

  ! quick read of metadata
  call adios_get_scalar(file_handle_adios, trim(name), scalar, ier)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error adios: could not read parameter: ',trim(name)
    print *,trim(err_message)
    stop 'Error adios helper read scalar'
  endif

  ! single read of an integer
  !sel = 0
  !call adios_schedule_read(file_handle_adios, sel, trim(scalar_name), 0, 1, scalar, ier)
  !if (ier /= 0) then
  !  call adios_errmsg(err_message)
  !  print * ,'Error adios: could not read parameter: ',trim(scalar_name)
  !  print *,trim(err_message)
  !  stop 'Error adios helper read scalar'
  !endif
  !call adios_perform_reads(file_handle_adios, ier)
  !if (ier /= 0 ) stop 'Error helper adios read scalar failed'

  end subroutine read_adios_scalar_int_only_rank

#endif

!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine show_adios_file_variables(filename)

! file inquiry showing all adios file variables

  use adios_read_mod, only: adios_read_open_file,adios_read_init_method,ADIOS_READ_METHOD_BP, &
    adios_inq_file,adios_inq_varnames,adios_inq_attrnames

  implicit none

  character(len=*),intent(in) :: filename

  ! local parameters
  integer :: i,ier
  integer :: variable_count, attribute_count, group_count
  integer :: timestep_first, timestep_last
  character (len=128), dimension(:), allocatable :: fnamelist

  ! file inquiry
  call adios_inq_file(file_handle_adios,variable_count,attribute_count,timestep_first,timestep_last,ier)
  if (ier /= 0) then
    ! show error message
    call adios_errmsg(err_message)
    print *,'Error inquiring adios file for reading: ',trim(filename)
    print *,trim(err_message)
    stop 'Error inquiring adios file'
  endif

  ! variables
  if (variable_count > 0) then
    allocate (fnamelist(variable_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets variable names
    call adios_inq_varnames(file_handle_adios, fnamelist, ier)

    ! user output
    print *,'variables: ',variable_count
    do i = 1,variable_count
      print *,'  ',trim(fnamelist(i))
    enddo
    deallocate(fnamelist)
  else
    print *,'  no variables'
  endif
  print *

  ! attributes
  if (attribute_count > 0) then
    allocate (fnamelist(variable_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets attribute names
    call adios_inq_attrnames(file_handle_adios, fnamelist, ier)

    ! user output
    print *,'attributes: ',attribute_count
    do i = 1,attribute_count
      print *,'  ',trim(fnamelist(i))
    enddo
    deallocate(fnamelist)
  else
    print *,'  no attributes'
  endif
  print *

  ! groups
  call adios_inq_ngroups(file_handle_adios, group_count, ier)

  if (group_count > 0) then
    allocate (fnamelist(group_count),stat=ier)
    if (ier /= 0) stop 'Error allocating namelist array'

    ! gets group names
    call adios_inq_groupnames(file_handle_adios, fnamelist, ier)

    ! user output
    print *,'groups: ',group_count
    do i = 1,group_count
      print *,'  ',trim(fnamelist(i))
    enddo
    deallocate (fnamelist)
  else
    print *,'  no groups'
  endif
  print *

  end subroutine show_adios_file_variables

#endif


!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine read_adios_array_gll(rank,nspec,array_name,array_gll)

  use adios_read_mod, only: adios_get_scalar,adios_selection_boundingbox,adios_schedule_read,adios_perform_reads

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: rank
  integer,intent(in) :: nspec

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(out) :: array_gll
  character(len=*),intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer :: local_dim
  integer(kind=8) :: sel
  integer(kind=8) :: sel_start(1),count_ad(1)

  array_gll(:,:,:,:) = 0.0_CUSTOM_REAL

  call adios_get_scalar(file_handle_adios, trim(array_name)//"/local_dim",local_dim, ier)
  if (ier /= 0 ) then
    call adios_errmsg(err_message)
    print *,'Error adios: reading array ',trim(array_name),' failed'
    print *,trim(err_message)
    stop 'Error adios helper read array'
  endif

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  sel_start(1) = local_dim * rank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox(sel, 1, sel_start, count_ad)

  ! reads selected array
  call adios_schedule_read(file_handle_adios, sel,trim(array_name)//"/array",0, 1, array_gll, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(array_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(file_handle_adios, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(array_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  end subroutine read_adios_array_gll

#endif

!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine read_adios_array_gll_int(rank,nspec,array_name,array_gll)

  use adios_read_mod, only: adios_get_scalar,adios_selection_boundingbox,adios_schedule_read,adios_perform_reads

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: rank
  integer,intent(in) :: nspec

  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(out) :: array_gll
  character(len=*),intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer :: local_dim
  integer(kind=8) :: sel
  integer(kind=8) :: sel_start(1),count_ad(1)

  array_gll(:,:,:,:) = 0

  call adios_get_scalar(file_handle_adios, trim(array_name)//"/local_dim",local_dim, ier)
  if (ier /= 0 ) then
    print *,'Error adios: reading array ',trim(array_name),' failed'
    stop 'Error adios helper read array'
  endif

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  sel_start(1) = local_dim * rank
  count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox(sel, 1, sel_start, count_ad)

  ! reads selected array
  call adios_schedule_read(file_handle_adios, sel,trim(array_name)//"/array",0, 1, array_gll, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(array_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(file_handle_adios, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(array_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  end subroutine read_adios_array_gll_int

#endif

!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine read_adios_array_1d(rank,nsize,array_name,array_1d)

  use adios_read_mod, only: adios_get_scalar,adios_selection_boundingbox,adios_schedule_read,adios_perform_reads

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: rank
  integer,intent(in) :: nsize

  real(kind=CUSTOM_REAL),dimension(nsize),intent(out) :: array_1d
  character(len=*),intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer :: local_dim
  integer(kind=8) :: sel
  integer(kind=8) :: sel_start(1),count_ad(1)
  !integer :: offset

  ! allow any rank to read from another rank-segment
  !! checks rank
  !!if (myrank_adios /= rank) &
  !!  stop 'Error invalid rank for reading adios GLL array'

  ! initializes
  array_1d(:) = 0

  ! gets local_dim metadata
  call adios_get_scalar(file_handle_adios, trim(array_name)//"/local_dim",local_dim, ier)
  if (ier /= 0 ) then
    print *,'Error adios: reading array ',trim(array_name),' failed'
    stop 'Error adios helper read array'
  endif

  sel_start(1) = local_dim * rank
  count_ad(1) = nsize

  ! gets offset (offset is the same as: local_dim * rank)
  !call adios_selection_writeblock(sel,rank)
  !call adios_schedule_read(file_handle_adios, sel, trim(array_name)//"/offset", 0, 1, offset, ier)
  !if (ier /= 0 ) stop 'Error adios: reading offset'
  !call adios_perform_reads(file_handle_adios, ier)
  !if (ier /= 0 ) stop 'Error adios: perform reading mesh file offsets failed'
  !sel_start(1) = offset
  !count_ad(1) = nsize

  call adios_selection_boundingbox(sel, 1, sel_start, count_ad)

  ! reads selected array
  call adios_schedule_read(file_handle_adios, sel,trim(array_name)//"/array",0, 1, array_1d, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(array_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(file_handle_adios, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(array_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  end subroutine read_adios_array_1d

#endif


!
!---------------------------------------------------------------------------------
!

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine read_adios_array_1d_int(rank,nsize,array_name,array_1d)

  use adios_read_mod, only: adios_get_scalar,adios_selection_boundingbox,adios_schedule_read,adios_perform_reads

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: rank
  integer,intent(in) :: nsize

  integer,dimension(nsize),intent(out) :: array_1d
  character(len=*),intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer :: local_dim
  integer(kind=8) :: sel
  integer(kind=8) :: sel_start(1),count_ad(1)
  !integer :: offset

  array_1d(:) = 0

  call adios_get_scalar(file_handle_adios, trim(array_name)//"/local_dim",local_dim, ier)
  if (ier /= 0 ) then
    print *,'Error adios: reading array ',trim(array_name),' failed'
    stop 'Error adios helper read array'
  endif

  sel_start(1) = local_dim * rank
  count_ad(1) = nsize

  call adios_selection_boundingbox(sel, 1, sel_start, count_ad)

  ! reads selected array
  call adios_schedule_read(file_handle_adios, sel,trim(array_name)//"/array", 0, 1, array_1d, ier)
  if (ier /= 0 ) then
    print *,'Error adios: scheduling read of array ',trim(array_name),' failed'
    stop 'Error adios helper schedule read array'
  endif

  call adios_perform_reads(file_handle_adios, ier)
  if (ier /= 0 ) then
    print *,'Error adios: performing read of array ',trim(array_name),' failed'
    stop 'Error adios helper perform read array'
  endif

  end subroutine read_adios_array_1d_int

#endif

!-------------------------------------------------------------------------------
! ADIOS writing
!-------------------------------------------------------------------------------

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine write_adios_scalar_int(scalar_name,scalar)

! writes a single scalar

  use adios_write_mod, only: adios_write

  implicit none

  integer, intent(in) :: scalar
  character(len=*),intent(in) :: scalar_name

  ! local parameters
  integer :: ier

  ! checks if file handle valid
  if (file_handle_adios == 0) stop 'Invalid ADIOS file handle in write_adios_scalar_int()'

  ! writes scalar (either buffered or directly to disk)
  call adios_write(file_handle_adios, trim(scalar_name), scalar, ier)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error adios: could not write parameter: ',trim(scalar_name)
    print *,trim(err_message)
    stop 'Error adios helper write scalar'
  endif

  end subroutine write_adios_scalar_int

#endif

!-------------------------------------------------------------------------------

! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools
#ifdef ADIOS_INPUT

  subroutine write_adios_array_gll(rank,nspec,array_name,array_gll)

! writes a GLL array

  use adios_helpers_writers_mod, only: write_adios_global_1d_array

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer,intent(in) :: rank
  integer,intent(in) :: nspec

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: array_gll
  character(len=*),intent(in) :: array_name

  ! local parameters
  integer :: ier
  integer :: local_dim

  ! checks if file handle valid
  if (file_handle_adios == 0) stop 'Invalid ADIOS file handle in write_adios_array_gll()'

  ! array size
  local_dim = NGLLX * NGLLY * NGLLZ * nspec

  ! writes array
  call write_adios_global_1d_array(file_handle_adios,rank,sizeprocs_adios,local_dim,trim(array_name),array_gll)
  if (ier /= 0) then
    call adios_errmsg(err_message)
    print *,'Error adios: could not write array parameter: ',trim(array_name)
    print *,trim(err_message)
    stop 'Error adios helper write array'
  endif

  end subroutine write_adios_array_gll

#endif

end module manager_adios


