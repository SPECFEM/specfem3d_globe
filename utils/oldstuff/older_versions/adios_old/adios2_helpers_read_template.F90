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


!===============================================================================
!> Helpers to read data with ADIOS2
!!
!-------------------------------------------------------------------------------
module adios2_helpers_read_mod

  use adios2

  implicit none

  public :: check_adios2_err
!
!  public :: read_adios2_real_array
!  public :: read_adios2_double_array
!  public :: read_adios2_integer_array
!  public :: read_adios2_long_array
!  public :: read_adios2_logical_array
!  public :: read_adios2_array

  !! read_adios2_array(file, io, varname, ndims, start, count, data, src,func,step)
  !!
  !! Read (a selection of) variable data from a file
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
  !! \param func Calling function string (for error print)
  !! \param step (optional) step to read from a multi-step file

  interface read_adios2
    module procedure read_adios2_integer_1d
    module procedure read_adios2_integer_2d
    module procedure read_adios2_integer_3d
    module procedure read_adios2_integer_4d
    module procedure read_adios2_integer_5d

    module procedure read_adios2_long_1d
    module procedure read_adios2_long_2d
    module procedure read_adios2_long_3d
    module procedure read_adios2_long_4d
    module procedure read_adios2_long_5d

    module procedure read_adios2_real_1d
    module procedure read_adios2_real_2d
    module procedure read_adios2_real_3d
    module procedure read_adios2_real_4d
    module procedure read_adios2_real_5d

    module procedure read_adios2_double_1d
    module procedure read_adios2_double_2d
    module procedure read_adios2_double_3d
    module procedure read_adios2_double_4d
    module procedure read_adios2_double_5d
  end interface read_adios2

contains


!===============================================================================

!> Check ADIOS return code, print a message and abort on error
!! \param adios2_err The error code considered.
  subroutine check_adios2_err(myrank, adios2_err, sourcefile, func, msg)

  implicit none
  integer, intent(in) :: myrank, adios2_err
  character(len=*), intent(in) :: sourcefile, func, msg

  if (adios2_err /= adios2_error_none) then
    print *, "process ", myrank, "has ADIOS2 error ",adios2_err,' in ', &
    trim(sourcefile), ':', trim(func), ": ", trim(msg)
    stop 'ADIOS2 error'
  endif

  end subroutine check_adios2_err


!===============================================================================

#define CODE_ARGS \
  use adios2                                                    _NL_\
  implicit none                                                 _NL_\
  type(adios2_engine), intent(in)   :: file                     _NL_\
  type(adios2_io), intent(in)       :: io                       _NL_\
  character(len=*), intent(in)      :: varname                  _NL_\
  integer(kind=4), intent(in)       :: ndims                    _NL_\
  integer(kind=8), dimension(1), intent(in) :: start, count     _NL_\
  integer, intent(in)               :: rank                     _NL_\
  character(len=*), intent(in)      :: src, func                _NL_\
  integer(kind=8), intent(in), optional :: step


#define CODE_BODY \
  integer                 :: ier                                                               _NL_\
  type(adios2_variable)   :: var                                                               _NL_\
  call adios2_inquire_variable(var, io, varname, ier)                                          _NL_\
  call check_adios2_err(rank, ier, src, func, "Inquire variable "//trim(varname))              _NL_\
  call adios2_set_selection(var, ndims, start, count, ier)                                     _NL_\
  if (present(step)) then                                                                      _NL_\
    call adios2_set_step_selection(var, step, 1_8, ier)                                        _NL_\
    call check_adios2_err(rank, ier, src, func, "Set step variable("//trim(varname)//")")      _NL_\
  endif                                                                                        _NL_\
  call adios2_get(file, var, data, adios2_mode_sync, ier)                                      _NL_\
  call check_adios2_err(rank, ier, src, func, "Read variable("//trim(varname)//")")


!# SUB(TYPENAME,TYPEDEF,NDIMS,DIMDEF) \
!subroutine read_adios2_##TYPENAME##_##NDIMS##d(file, io, varname, ndims, start, count, data, myrank, src, func, step)  _NL_\
!##CODE _ARGS                                                                                     _NL_\
!  TYPEDEF, dimension(DIMDEF), intent(out) :: data                                               _NL_\
!##CODE _BODY                                                                                     _NL_\
!end subroutine read_adios2_##TYPENAME##_##NDIMS##d

! #  !define subroutineTYPE(TYPENAME,TYPEDEF) _NL_\
! # #SUB(TYPENAME,TYPEDEF,1,:)              _NL_\
! # #(TYPENAME,TYPEDEF,2,@@:,:@@)          _NL_\


!! subroutineTYPE(integer,integer(kind=4))

!===============================================================================
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
!! \param func Calling function string (for error print)
  subroutine read_adios2_real_1d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=4), dimension(:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_real_1d


!===============================================================================
  subroutine read_adios2_real_2d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=4), dimension(:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_real_2d


!===============================================================================
  subroutine read_adios2_real_3d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=4), dimension(:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_real_3d


!===============================================================================
  subroutine read_adios2_real_4d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=4), dimension(:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_real_4d


!===============================================================================
  subroutine read_adios2_real_5d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=4), dimension(:,:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_real_5d


!===============================================================================
  subroutine read_adios2_double_1d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=8), dimension(:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_double_1d


!===============================================================================
  subroutine read_adios2_double_2d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=8), dimension(:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_double_2d


!===============================================================================
  subroutine read_adios2_double_3d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=8), dimension(:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_double_3d


!===============================================================================
  subroutine read_adios2_double_4d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=8), dimension(:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_double_4d


!===============================================================================
  subroutine read_adios2_double_5d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  real(kind=8), dimension(:,:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_double_5d


!===============================================================================
  subroutine read_adios2_integer_1d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer, dimension(:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_integer_1d


!===============================================================================
  subroutine read_adios2_integer_2d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer, dimension(:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_integer_2d


!===============================================================================
  subroutine read_adios2_integer_3d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer, dimension(:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_integer_3d


!===============================================================================
  subroutine read_adios2_integer_4d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer, dimension(:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_integer_4d


!===============================================================================
  subroutine read_adios2_integer_5d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer, dimension(:,:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_integer_5d


!===============================================================================
  subroutine read_adios2_long_1d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer(kind=8), dimension(:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_long_1d


!===============================================================================
  subroutine read_adios2_long_2d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer(kind=8), dimension(:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_long_2d


!===============================================================================
  subroutine read_adios2_long_3d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer(kind=8), dimension(:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_long_3d


!===============================================================================
  subroutine read_adios2_long_4d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer(kind=8), dimension(:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_long_4d


!===============================================================================
  subroutine read_adios2_long_5d(file, io, varname, ndims, start, count, data, rank, src,func,step)
  CODE_ARGS
  integer(kind=8), dimension(:,:,:,:,:), intent(out) :: data
  CODE_BODY
  end subroutine read_adios2_long_5d


end module adios2_helpers_read_mod
