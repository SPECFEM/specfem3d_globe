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

!> Tools to setup and cleanup ADIOS2,
!> contains wrapper subroutines for common adios calls
!
! note: adios library calls use a format like "adios_do_something***()"
!       our own wrapper functions thus will rather use something like "do_something_adios***()"
!       to better distinguish between library functions and wrappers.

module manager_adios2

#ifdef HAVE_ADIOS2
  use adios2, only: adios2_adios, adios2_engine
#endif
  implicit none

  private

  ! MPI copies of communicator and rank
  integer,public :: comm_adios2
  integer,public :: myrank_adios2
  integer,public :: sizeprocs_adios2

  ! initialized flag
  logical :: is_adios2_initialized

  ! ADIOS Real type based on CUSTOM_REAL
  integer, public :: adios2_custom_real

#ifdef HAVE_ADIOS2
  ! adios2 handlers
  type(adios2_adios), public:: adios2obj
  ! File handlers that needs to be closed at exit
  type(adios2_engine), public:: adios2_file_fwdatt
#endif

  ! debugging
  logical,parameter :: DEBUG = .false.

  ! public routines
  public :: initialize_adios2
  public :: finalize_adios2

  ! only available with ADIOS2 compilation support
  ! to clearly separate adios version and non-adios version of same tools
#ifdef HAVE_ADIOS2xxx
  public :: check_adios2_err
  public :: read_adios2_global_real_1d_array
  public :: read_adios2_global_double_1d_array
  public :: read_adios2_global_integer_1d_array
  public :: read_adios2_global_long_1d_array
  public :: read_adios2_global_logical_1d_array
  public :: read_adios2_global_string_1d_array
  public :: read_adios2_global_1d_array
#endif

contains


!-------------------------------------------------------------------------------
!
! public ADIOS2 wrapper routines (also available without adios compilation support)
!
!-------------------------------------------------------------------------------

  subroutine initialize_adios2()

!> Initialize ADIOS and setup the xml output file

#ifdef HAVE_ADIOS2
  use adios2, only: adios2_init, adios2_debug_mode_on, adios2_type_real4, adios2_type_real8
#endif
  use constants, only: CUSTOM_REAL, SIZE_REAL, SIZE_DOUBLE

  implicit none

  ! local parameters
#ifdef HAVE_ADIOS2
  integer :: ier
#endif

  ! initializes
  is_adios2_initialized = .false.
  call world_get_comm(comm_adios2)
  call world_rank(myrank_adios2)
  call world_size(sizeprocs_adios2)
  adios2_custom_real = 0;


#ifdef HAVE_ADIOS2

  ! Create adios handler passing the communicator, debug mode and error flag
  ! adios2 duplicates the communicator for its internal use
  call adios2_init(adios2obj, comm_adios2, adios2_debug_mode_on, ier)

  ! sets flag
  is_adios2_initialized = .true.

  if (CUSTOM_REAL == SIZE_REAL) then
    adios2_custom_real = adios2_type_real4
  else
    adios2_custom_real = adios2_type_real8
  endif

#else

  ! compilation without ADIOS support
  if (myrank_adios2 == 0) then
    print *, "Error: ADIOS2 enabled without ADIOS2 Support."
    print *, "To enable ADIOS2 support, reconfigure with --with-adios2 flag."
  endif
  ! safety stop
  call exit_MPI(myrank_adios2,"Error ADIOS2 manager: intitialize called without compilation support")

#endif

  end subroutine initialize_adios2

!
!-------------------------------------------------------------------------------
!

  subroutine finalize_adios2()

!> Finalize ADIOS. Must be called once everything is written down.

#ifdef HAVE_ADIOS2
  use adios2, only: adios2_close, adios2_finalize
#endif

  implicit none

  ! local parameters
#ifdef HAVE_ADIOS2
  integer :: ier
#endif

  ! synchronizes all first
  call synchronize_all_comm(comm_adios2)

#ifdef HAVE_ADIOS2
  ! close files that are growing until the end of run
  if (adios2_file_fwdatt%valid) then
     call adios2_close(adios2_file_fwdatt, ier)
  endif

  ! finalize
  call adios2_finalize(adios2obj, ier)
  if (ier /= 0 ) stop 'Error cleaning up ADIOS2: calling adios2_finalize() routine failed'

#else
  ! safety stop
  call exit_MPI(myrank_adios2,"Error ADIOS2 manager: finalize called without compilation support")
#endif

  end subroutine finalize_adios2


!-------------------------------------------------------------------------------
!
! ADIOS wrapper routines (only available with adios compilation support)
!
!-------------------------------------------------------------------------------
#ifdef HAVE_ADIOS2xxx



subroutine read_var_adios2(file, io, varname, ndims, start, count, data, src, func )

!> Read (a selection of) variable data from a file
!! (set up a selection and schedule for reading in data)
!! Note: This function calls check_adios2_err which aborts the program on error
!! \param file adios2_engine object from adios2_open
!! \param io adios2_io object from adios2_declare_io
!! \param varname Name of the variable
!! \param ndims Number of dimensions of the variable
!! \param start Offsets in global array for reading
!! \param count Local sizes for reading
!! \param data Pre-allocated array to receive the data from file
!! \param src Sourcefile string (for error print)
!! \param func Calling Function string (for error print)


  use adios2
  implicit none
  type(adios2_engine), intent(in)   :: file
  type(adios2_io), intent(in)       :: io
  character(len=:), intent(in)      :: varname
  integer(kind=4), intent(in)       :: ndims
  integer(kind=8), dimension(1), intent(in) :: start, count
  (kind), dimension(:), intent(out) :: data
  character(len=:), intent(in)      :: sourcefile, func, msg

  integer                 :: ier
  type(adios2_variable)   :: var

  call adios2_inquire_variable(var, io, varname, ier)
  call check_adios2_err(myrank_adios2, ier, src, unc, "Inquire variable "//trim(varname))
  call adios2_set_selection(var, ndims, start, count, ier)
  call adios2_get(file, var, data, adios2_mode_sync, ier)
  call check_adios2_err(myrank_adios2, ier, src, unc, "adios2_get("//trim(varname)//")")

end subroutine read_var_adios2



#endif /* HAVE_ADIOS2 */
end module manager_adios2


