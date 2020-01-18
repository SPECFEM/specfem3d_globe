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

#ifdef USE_ADIOS2
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
  integer, public :: adios2_CUSTOM_REAL

#ifdef USE_ADIOS2
  ! adios2 handlers
!  type(adios2_adios), public:: adios2obj

  ! File handlers that needs to be closed at exit
!  type(adios2_engine), public:: adios2_file_fwdatt

  ! file handle for read/write
!  type(adios2_engine),public :: adios2_file_handle
#endif

  ! debugging
  logical,parameter :: DEBUG = .false.

  ! public routines

  ! only available with ADIOS2 compilation support
  ! to clearly separate adios version and non-adios version of same tools
#ifdef USE_ADIOS2
  !public :: read_adios2_global_real_1d_array
  !public :: read_adios2_global_double_1d_array
  !public :: read_adios2_global_integer_1d_array
  !public :: read_adios2_global_long_1d_array
  !public :: read_adios2_global_logical_1d_array
  !public :: read_adios2_global_string_1d_array
  !public :: read_adios2_global_1d_array
#endif

contains


!-------------------------------------------------------------------------------
!
! ADIOS wrapper routines (only available with adios compilation support)
!
!-------------------------------------------------------------------------------
#ifdef USE_ADIOS2


! only available with ADIOS compilation support
! to clearly separate adios version and non-adios version of same tools

! not used yet...
!
!  subroutine read_var_adios2(file, io, varname, ndims, start, count, data, src, func )
!
!!> Read (a selection of) variable data from a file
!!! (set up a selection and schedule for reading in data)
!!! Note: This function calls check_adios2_err which aborts the program on error
!!! \param file adios2_engine object from adios2_open
!!! \param io adios2_io object from adios2_declare_io
!!! \param varname Name of the variable
!!! \param ndims Number of dimensions of the variable
!!! \param start Offsets in global array for reading
!!! \param count Local sizes for reading
!!! \param data Pre-allocated array to receive the data from file
!!! \param src Sourcefile string (for error print)
!!! \param func Calling function string (for error print)
!
!  use adios2
!
!  implicit none
!
!  type(adios2_engine), intent(in)   :: file
!  type(adios2_io), intent(in)       :: io
!  character(len=:), intent(in)      :: varname
!  integer(kind=4), intent(in)       :: ndims
!  integer(kind=8), dimension(1), intent(in) :: start, count
!  (kind), dimension(:), intent(out) :: data
!  character(len=:), intent(in)      :: sourcefile, func, msg
!
!  integer                 :: ier
!  type(adios2_variable)   :: var
!
!  call adios2_inquire_variable(var, io, varname, ier)
!  call check_adios2_err(myrank_adios2, ier, src, func, "Inquire variable "//trim(varname))
!
!  call adios2_set_selection(var, ndims, start, count, ier)
!
!  call adios2_get(file, var, data, adios2_mode_sync, ier)
!  call check_adios2_err(myrank_adios2, ier, src, func, "adios2_get("//trim(varname)//")")
!
!  end subroutine read_var_adios2

#endif

end module manager_adios2
