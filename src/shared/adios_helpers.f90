!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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


!-------------------------------------------------------------------------------
!> \file adios_helpers.f90
!! \brief Helpers to set up adios features.
!! \author MPBL
!-------------------------------------------------------------------------------

!===============================================================================
!> Get the ADIOS error message from an adios error number if there is an error.
!! \param adios_err The error code considered.
subroutine check_adios_err(myrank, adios_err)
  use adios_read_mod
  implicit none
  integer, intent(in) :: myrank, adios_err
  character(len=1024) :: msg

  if (adios_err /= 0) then
    call adios_errmsg(msg)
    print *, "process: ", myrank, ", error: ", msg
    stop
  endif
end subroutine check_adios_err


!===============================================================================
!> Define an ADIOS scalar double precision variable and autoincrement
!! the adios group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_double_scalar (adios_group, name, path, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 6 == real(kind=8)
  call adios_define_var (adios_group, name, path, 6,  "", "", "", varid)
  group_size_inc = group_size_inc + 8
end subroutine define_adios_double_scalar

!===============================================================================
!> Define an ADIOS scalar integer variable and autoincrement the adios
!! group size by (4).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_integer_scalar (adios_group, name, path, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 2 == integer(kind=4)
  call adios_define_var (adios_group, name, path, adios_integer,  "", "", "", varid)
  group_size_inc = group_size_inc + 4
end subroutine define_adios_integer_scalar

!===============================================================================
!> Define an ADIOS scalar byte variable and autoincrement the adios
!! group size by (1).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_byte_scalar (adios_group, name, path, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 0 == byte == any_data_type(kind=1)
  call adios_define_var (adios_group, name, path, 0,  "", "", "", varid)
  group_size_inc = group_size_inc + 1
end subroutine define_adios_byte_scalar

!===============================================================================
!> Define a local ADIOS array of integers and autoincrement the adios
!! group size by (4 * number of elements).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param dim The number of elements in the 1D array. Required to
!!            correctly increment adios group size.
!! \param dim_str The "stringified" version of dim. Needed by adios
!!                to define variables
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_integer_local_array1D (adios_group, name, path, dim, dim_str, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path, dim_str
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer, intent(in)              :: dim
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 2 == integer
  call adios_define_var (adios_group, name, path, 2,  dim_str, "", "", varid)
  group_size_inc = group_size_inc + 4*dim
end subroutine define_adios_integer_local_array1D

!===============================================================================
!> Define a local ADIOS array of doubles and autoincrement the adios
!! group size by (8 * number of elements).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param dim The number of elements in the 1D array. Required to
!!            correctly increment adios group size.
!! \param dim_str The "stringified" version of dim. Needed by adios
!!                to define variables
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
subroutine define_adios_double_local_array1D (adios_group, name, path, dim, dim_str, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path, dim_str
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer, intent(in)              :: dim
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 6 == real(kind=8)
  call adios_define_var (adios_group, name, path, 6, dim_str, "", "", varid)
  group_size_inc = group_size_inc + 8*dim
end subroutine define_adios_double_local_array1D

!===============================================================================
!> Define a local ADIOS string and autoincrement the adios
!! group size by (1 * string's length).
!! \param adios_group The adios group where the variables belongs
!! \param name The variable to be defined
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param len The length of the string(number of character. in Fortran
!!            it does not include a final '\0' -- null -- character)
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \note Adios string are scalar values counting for (1) byte. It is
!!       mandatory to increase the group size by the length of the
!!       string in order not to overlap 'data regions'.
subroutine define_adios_string (adios_group, name, path, length, group_size_inc)
  use adios_write_mod
  implicit none
  ! Arguments
  integer(kind=8),  intent(in)     :: adios_group
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer                          :: length
  ! Local Variables
  integer(kind=8)                  :: varid ! dummy variable, adios use var name

  ! adios: 9 == string
  call adios_define_var (adios_group, name, path, 9,  "", "", "", varid)
  group_size_inc = group_size_inc + 1*length
end subroutine define_adios_string

!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param len The local dimension of the array.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
subroutine define_adios_global_real_1d_array(adios_group, array_name, &
    local_dim, group_size_inc)
  use adios_write_mod
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_group
  character(len=*), intent(in) :: array_name
  integer, intent(in) :: local_dim
  integer(kind=8), intent(inout) :: group_size_inc
  ! Variables
  integer(kind=8) :: var_id

  call define_adios_integer_scalar (adios_group, "local_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "global_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "offset", array_name, &
      group_size_inc)
!  call adios_define_var(adios_group, "array", array_name, 6, &
!      "local_dim", "global_dim", "offset", var_id)
  call adios_define_var(adios_group, "array", array_name, 5, &
      array_name // "/local_dim", array_name // "/global_dim", &
      array_name // "/offset", var_id)
  group_size_inc = group_size_inc + local_dim*4
end subroutine define_adios_global_real_1d_array

!===============================================================================
!> Define a global ADIOS 1D integer array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param len The local dimension of the array.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
subroutine define_adios_global_integer_1d_array(adios_group, array_name, &
    local_dim, group_size_inc)
  use adios_write_mod
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_group
  character(len=*), intent(in) :: array_name
  integer, intent(in) :: local_dim
  integer(kind=8), intent(inout) :: group_size_inc
  ! Variables
  integer(kind=8) :: var_id

  call define_adios_integer_scalar (adios_group, "local_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "global_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "offset", array_name, &
      group_size_inc)
!  call adios_define_var(adios_group, "array", array_name, 6, &
!      "local_dim", "global_dim", "offset", var_id)
  call adios_define_var(adios_group, "array", array_name, 2, &
      array_name // "/local_dim", array_name // "/global_dim", &
      array_name // "/offset", var_id)
  group_size_inc = group_size_inc + local_dim*4
end subroutine define_adios_global_integer_1d_array

!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param len The local dimension of the array.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
subroutine define_adios_global_logical_1d_array(adios_group, array_name, &
    local_dim, group_size_inc)
  use adios_write_mod
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_group
  character(len=*), intent(in) :: array_name
  integer, intent(in) :: local_dim
  integer(kind=8), intent(inout) :: group_size_inc
  ! Variables
  integer(kind=8) :: var_id

  call define_adios_integer_scalar (adios_group, "local_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "global_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "offset", array_name, &
      group_size_inc)
!  call adios_define_var(adios_group, "array", array_name, 6, &
!      "local_dim", "global_dim", "offset", var_id)
  call adios_define_var(adios_group, "array", array_name, 1, &
      array_name // "/local_dim", array_name // "/global_dim", &
      array_name // "/offset", var_id)
  group_size_inc = group_size_inc + local_dim*1
end subroutine define_adios_global_logical_1d_array

!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param len The local dimension of the array.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
subroutine define_adios_global_double_1d_array(adios_group, array_name, &
    local_dim, group_size_inc)
  use adios_write_mod
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: adios_group
  character(len=*), intent(in) :: array_name
  integer, intent(in) :: local_dim
  integer(kind=8), intent(inout) :: group_size_inc
  ! Variables
  integer(kind=8) :: var_id

  call define_adios_integer_scalar (adios_group, "local_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "global_dim", array_name, &
      group_size_inc)
  call define_adios_integer_scalar (adios_group, "offset", array_name, &
      group_size_inc)
  call adios_define_var(adios_group, "array", array_name, 6, &
      array_name // "/local_dim", array_name // "/global_dim", &
      array_name // "/offset", var_id)
  group_size_inc = group_size_inc + local_dim*8
end subroutine define_adios_global_double_1d_array
