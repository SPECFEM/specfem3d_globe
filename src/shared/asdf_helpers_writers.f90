!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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


!===============================================================================
!> Helpers to set up adios features.
!! * Scalar definition
!! * Global arrays definition
!!
!! \note We do not define function to write scalars variables into adios
!!       since it is already a single function call.
!!
!! \author MPBL
!-------------------------------------------------------------------------------
module asdf_helpers_writers_mod
  implicit none

  private

  public :: write_adios_global_real_1d_array
  public :: write_adios_global_double_1d_array
  public :: write_adios_global_integer_1d_array
  public :: write_adios_global_long_1d_array
  public :: write_adios_global_logical_1d_array
  public :: write_adios_global_string_1d_array
  public :: write_adios_global_1d_array

  interface write_adios_global_real_1d_array
    module procedure write_adios_global_1d_real_1d
  end interface write_adios_global_real_1d_array

  interface write_adios_global_double_1d_array
    module procedure write_adios_global_1d_double_1d
  end interface write_adios_global_double_1d_array

  interface write_adios_global_integer_1d_array
    module procedure write_adios_global_1d_integer_1d
  end interface write_adios_global_integer_1d_array

  interface write_adios_global_long_1d_array
    module procedure write_adios_global_1d_long_1d
  end interface write_adios_global_long_1d_array

  interface write_adios_global_logical_1d_array
    module procedure write_adios_global_1d_logical_1d
  end interface write_adios_global_logical_1d_array

  interface write_adios_global_string_1d_array
    module procedure write_adios_global_1d_string_1d
  end interface write_adios_global_string_1d_array

  interface write_adios_global_1d_array
    module procedure write_adios_global_1d_integer_1d

    module procedure write_adios_global_1d_long_1d

    module procedure write_adios_global_1d_logical_1d

    module procedure write_adios_global_1d_real_1d

    module procedure write_adios_global_1d_double_1d
  end interface write_adios_global_1d_array

contains


!===============================================================================
subroutine write_1D_global_array_adios_dims(adios_handle, myrank, &
    local_dim, global_dim, offset, sizeprocs, path)
  use adios_write_mod

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs, local_dim, myrank
  integer, intent(in) :: global_dim, offset
  character(len=*), intent(in) :: path

  integer :: adios_err

  call adios_write(adios_handle, trim(path)// "/local_dim", &
                   local_dim, adios_err)
  call adios_write(adios_handle, trim(path)// "/global_dim", &
                   global_dim, adios_err)
  call adios_write(adios_handle, trim(path)// "/offset", &
                   offset, adios_err)
end subroutine write_1D_global_array_adios_dims


!===============================================================================
!> Schedule an ADIOS single precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
subroutine write_adios_global_1d_real_1d(adios_handle, myrank, sizeprocs, &
    local_dim, global_dim, offset, array_name, array)
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs, local_dim
  integer, intent(in) :: global_dim, offset
  character(len=*) :: array_name
  real, dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, global_dim, offset, sizeprocs, array_name)
  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
end subroutine write_adios_global_1d_real_1d


!===============================================================================
!> Schedule an ADIOS double precision global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
subroutine write_adios_global_1d_double_1d(adios_handle, myrank, sizeprocs, &
    local_dim, global_dim, offset, array_name, array)
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs, local_dim
  integer, intent(in) :: global_dim, offset
  character(len=*) :: array_name
  real(kind=8), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, global_dim, offset, sizeprocs, array_name)
  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
end subroutine write_adios_global_1d_double_1d


!===============================================================================
!> Schedule an ADIOS integer global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
subroutine write_adios_global_1d_integer_1d(adios_handle, myrank, sizeprocs, &
    local_dim, global_dim, offset, array_name, array)
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs, local_dim
  integer, intent(in) :: global_dim, offset
  character(len=*) :: array_name
  integer(kind=4), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, global_dim, offset, sizeprocs, array_name)
  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
end subroutine write_adios_global_1d_integer_1d


!===============================================================================
!> Schedule an ADIOS long global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
subroutine write_adios_global_1d_long_1d(adios_handle, myrank, sizeprocs, &
    local_dim, global_dim, offset, array_name, array)
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs, local_dim
  integer, intent(in) :: global_dim, offset
  character(len=*) :: array_name
  integer(kind=8), dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, global_dim, offset, sizeprocs, array_name)
  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
end subroutine write_adios_global_1d_long_1d


!===============================================================================
!> Schedule an ADIOS logical global 1D array for write
!! \param adios_handle The adios handle to the file to be written
!! \param myrank The rank of the MPI process involved
!! \param sizeprocs The number of MPI process in the communicator writing the
!!                  variable
!! \param local_dim The number of elements to be written by each process. Might
!!                  eventually be padded.
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param array_name The array name in the ADIOS file.
!! \param array The array to be written
subroutine write_adios_global_1d_logical_1d(adios_handle, myrank, sizeprocs, &
    local_dim, global_dim, offset, array_name, array)
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs, local_dim
  integer, intent(in) :: global_dim, offset
  character(len=*) :: array_name
  logical, dimension(:), intent(in) :: array
  ! Variables
  integer :: adios_err

  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, global_dim, offset, sizeprocs, array_name)
  call adios_write(adios_handle, trim(array_name)// "/array", array, adios_err)
end subroutine write_adios_global_1d_logical_1d


!===============================================================================
!string subroutine added
subroutine write_adios_global_1d_string_1d(adios_handle, myrank, sizeprocs, &
    local_dim, global_dim, offset, array_name, array)
  ! Parameters
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs, local_dim
  integer, intent(in) :: global_dim, offset
  character(len=*) :: array_name
  character(len=*), intent(in) :: array
  ! Variables
  integer :: adios_err

  print *,"tag2:",trim(array_name)
  print *,"tag2:",trim(array)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, global_dim, offset, sizeprocs, array_name)
  call adios_write(adios_handle, trim(array_name)// "/array", array(1:local_dim), adios_err)
end subroutine write_adios_global_1d_string_1d

end module asdf_helpers_writers_mod
