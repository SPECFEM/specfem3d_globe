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


!-------------------------------------------------------------------------------
!> Helpers to set up adios features.
!! * Scalar definition
!! * Global arrays definition
!!
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"


module adios_helpers_definitions_mod

  !use manager_adios
#if defined(USE_ADIOS)
  use adios_write_mod, only: adios_define_var, adios_long, adios_integer
#elif defined(USE_ADIOS2)
  use adios2
#endif

  use manager_adios, only: check_adios_err
#if defined(USE_ADIOS2)
  use manager_adios, only: sizeprocs_adios,myrank_adios
#endif

  implicit none

  private

  public :: define_adios_scalar
  public :: define_adios_global_real_1d_array
  public :: define_adios_global_double_1d_array
  public :: define_adios_global_integer_1d_array
  public :: define_adios_global_long_1d_array
  public :: define_adios_global_logical_1d_array
  public :: define_adios_global_string_1d_array
  public :: define_adios_local_string_1d_array
  public :: define_adios_global_array1D

  ! Generic interface to define scalar variables in ADIOS
  interface define_adios_scalar
    module procedure define_adios_scalar_double
    module procedure define_adios_scalar_float
    module procedure define_adios_scalar_integer
    module procedure define_adios_scalar_long
    module procedure define_adios_scalar_byte
  end interface define_adios_scalar

  interface define_adios_global_real_1d_array
    module procedure define_adios_global_1d_real_1d
    ! unused so far...
    !module procedure define_adios_global_1d_real_2d
    !module procedure define_adios_global_1d_real_3d
    !module procedure define_adios_global_1d_real_4d
    !module procedure define_adios_global_1d_real_5d
  end interface define_adios_global_real_1d_array

  interface define_adios_global_double_1d_array
    module procedure define_adios_global_1d_double_1d
    ! unused so far..
    !module procedure define_adios_global_1d_double_2d
    !module procedure define_adios_global_1d_double_3d
    !module procedure define_adios_global_1d_double_4d
    !module procedure define_adios_global_1d_double_5d
  end interface define_adios_global_double_1d_array

  interface define_adios_global_integer_1d_array
    module procedure define_adios_global_1d_int_1d
    ! unused so far..
    !module procedure define_adios_global_1d_int_2d
    !module procedure define_adios_global_1d_int_3d
    !module procedure define_adios_global_1d_int_4d
    !module procedure define_adios_global_1d_int_5d
  end interface define_adios_global_integer_1d_array

  interface define_adios_global_long_1d_array
    module procedure define_adios_global_1d_long_1d
    ! unused so far..
    !module procedure define_adios_global_1d_long_2d
    !module procedure define_adios_global_1d_long_3d
    !module procedure define_adios_global_1d_long_4d
    !module procedure define_adios_global_1d_long_5d
  end interface define_adios_global_long_1d_array

  interface define_adios_global_logical_1d_array
    module procedure define_adios_global_1d_logical_1d
    ! unused so far..
    !module procedure define_adios_global_1d_logical_2d
    !module procedure define_adios_global_1d_logical_3d
    !module procedure define_adios_global_1d_logical_4d
    !module procedure define_adios_global_1d_logical_5d
  end interface define_adios_global_logical_1d_array

  interface define_adios_global_string_1d_array
    module procedure define_adios_global_1d_string_1d
  end interface define_adios_global_string_1d_array

  interface define_adios_local_string_1d_array
    module procedure define_adios_local_1d_string_1d
  end interface define_adios_local_string_1d_array

  ! Cannot include an interface in another interface
  interface define_adios_global_array1D
    module procedure define_adios_global_1d_int_1d
    module procedure define_adios_global_1d_int_2d
    module procedure define_adios_global_1d_int_3d
    module procedure define_adios_global_1d_int_4d
    module procedure define_adios_global_1d_int_5d

    module procedure define_adios_global_1d_long_1d
    module procedure define_adios_global_1d_long_2d
    module procedure define_adios_global_1d_long_3d
    module procedure define_adios_global_1d_long_4d
    module procedure define_adios_global_1d_long_5d

    module procedure define_adios_global_1d_logical_1d
    module procedure define_adios_global_1d_logical_2d
    module procedure define_adios_global_1d_logical_3d
    module procedure define_adios_global_1d_logical_4d
    module procedure define_adios_global_1d_logical_5d

    module procedure define_adios_global_1d_real_1d
    module procedure define_adios_global_1d_real_2d
    module procedure define_adios_global_1d_real_3d
    module procedure define_adios_global_1d_real_4d
    module procedure define_adios_global_1d_real_5d

    module procedure define_adios_global_1d_double_1d
    module procedure define_adios_global_1d_double_2d
    module procedure define_adios_global_1d_double_3d
    module procedure define_adios_global_1d_double_4d
    module procedure define_adios_global_1d_double_5d

    module procedure define_adios_global_1d_string_1d

  end interface define_adios_global_array1D

contains


!===============================================================================
!
! scalars
!
!===============================================================================

!===============================================================================
!> Define an ADIOS scalar double precision variable and autoincrement
!! the adios group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!!            ignored.
!!
!! \note 'name' and 'var' are written as successive arguments on purpose.
!!       One should be able to define a macro such as:
!!       #define STRINGIFY_VAR(x) #x, x
!!       Calling define_adi os_double_scalar with such a macro will be done as:
!!       call define_adios_scalar_double(group, size, path, STRINGIFY_VAR(x))
!!       as STRINGIFY_VAR(x) expand as:
!!       "x", x
!!       x being the variable name inside the code.
subroutine define_adios_scalar_double (adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in)    :: name, path
  real(kind=8),     intent(in)    :: var
  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  real(kind=8) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_double: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_double()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 6 == real(kind=8)
  call adios_define_var (adios_group, trim(name), trim(path), 6,  '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines scalar as global variable (same for all processes, would be enough to be stored by master process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real8, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real8, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real8, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 8

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_double


!===============================================================================
!> Define an ADIOS scalar single precision variable and autoincrement
!! the adios group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_float(adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in)    :: name, path
  real(kind=4),     intent(in)    :: var
  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  real(kind=4) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_float: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_float()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 6 == real(kind=8), 5 == real(kind=4)
  call adios_define_var (adios_group, trim(name), trim(path), 5,  '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines scalar as global variable (same for all processes, would be enough to be stored by master process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real4, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real4, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_real4, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 4

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_float


!===============================================================================
!> Define an ADIOS scalar integer variable and autoincrement the adios
!! group size by (4).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_integer(adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)    :: adios_group
#endif
  character(len=*), intent(in)    :: name, path
  integer(kind=8),  intent(inout) :: group_size_inc
  integer(kind=4),  intent(in)    :: var
  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  integer(kind=4) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_integer',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_integer()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

  !debug
  !print *,'debug adios: ',myrank_adios,' define integer scalar: ',trim(full_name),' path: ',trim(path),' name: ',trim(name)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 2 ~ integer(kind=4)
  call adios_define_var (adios_group, trim(name), trim(path), adios_integer, '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note 1: we won't store the variable object, but get it back by adios2_inquire_** for writing
  !
  ! note 2: we will store a scalar as a 1-D array with a single entry instead.
  !         this is due to appending to a file will increase the step count for variables.
  !         retrieving local scalar variables would use adios2_set_block_selection() which then fails in such cases with an error:
  !          " ERROR: invalid blockID 0 from steps start 0 in variable reg2/nspec,
  !                   check argument to Variable < T>::SetBlockID, in call to Get "
  !
  !         however, using 1-D arrays will use adios2_set_selection() which succeeds also for appended variables.
  !         until adios2 fixes this, we will use the 1-D work-around.

  ! defines scalar as global variable (same for all processes, would be enough to be stored by master process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer4, ier)
  !
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer4, &
  !                                1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer4, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 4

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_integer


!===============================================================================
!> Define an ADIOS scalar long integer variable and autoincrement the adios
!! group size by (8).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_long(adios_group, group_size_inc, path, name, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  integer(kind=8),  intent(in)  :: var
  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  integer(kind=8) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_long: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_long()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 2 ~ integer(kind=4)
  call adios_define_var (adios_group, trim(name), trim(path), adios_long, '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines global variable
  ! defines scalar as global variable (same for all processes, would be enough to be stored by master process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer8, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer8, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer8, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 8

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_long

!===============================================================================
!> Define an ADIOS scalar byte variable and autoincrement the adios
!! group size by (1).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_scalar_double()
subroutine define_adios_scalar_byte (adios_group, group_size_inc, name, path, var)

  implicit none
  ! Arguments
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in)     :: name, path
  integer(kind=8),  intent(inout)  :: group_size_inc
  ! note: byte is non-standard gnu Fortran
  !byte,     intent(in)             :: var
  integer(kind=1),  intent(in)     :: var
  ! Local Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: varid ! dummy variable, adios use var name
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
  !integer(kind=8), parameter :: count_dims(1) = (/ 1 /)
  integer(kind=8) :: ldim(1),gdim(1),offs(1)
#endif
  integer(kind=1) :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_scalar_byte: ',trim(name))

  ! check
  if (len_trim(name) == 0) stop 'Error adios: invalid name in define_adios_scalar_byte()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(name)
  else
    full_name = trim(path) // '/' // trim(name)
  endif

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! adios: 0 == byte == any_data_type(kind=1)
  call adios_define_var (adios_group, trim(name), trim(path), 0,  '', '', '', varid)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines global variable
  ! defines scalar as global variable (same for all processes, would be enough to be stored by master process myrank==0)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer1, ier)
  ! defines scalar as local variable (value can vary between processes)
  ! > call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer1, &
  !                               1, adios2_null_dims, adios2_null_dims, count_dims, adios2_constant_dims, ier)
  !
  ! defines scalar as 1-d array with single entry
  ldim(1) = 1
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_integer1, &
                              1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + 1

  ! to avoid compiler warnings
  idummy = var

end subroutine define_adios_scalar_byte


!===============================================================================
!
! arrays
!
!===============================================================================

!===============================================================================
!> Define the dimensions that will be written along a global array in ADIOS.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
subroutine define_adios_global_dims_1d(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in) :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in) :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc

  TRACE_ADIOS_L2_ARG('define_adios_global_dims_1d: ',trim(array_name))

  ! array_name should be defined
  if (len_trim(array_name) == 0) stop 'Error adios: invalid array_name in define_adios_global_dims_1d()'

  ! uses local_dim as dummy variable
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "local_dim", local_dim)
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "global_dim", int(local_dim,kind=8)) ! long type
  call define_adios_scalar(adios_group, group_size_inc, trim(array_name), "offset", local_dim)

end subroutine define_adios_global_dims_1d


!===============================================================================
!> Define a real global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_real(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  integer :: ier
#endif

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_real: ',trim(array_name))

  ! Define the dimensions of the array.
  ! local_dim used as a dummy variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, "array", trim(array_name), 5, &
                        trim(array_name) // "/local_dim", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = local_dim
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_real4, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 4

  end subroutine define_adios_global_1d_generic_real

!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_1d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_2d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_3d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_4d


!===============================================================================
!> Define a global ADIOS 1D real array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_real_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real, dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE real
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_real_5d


!===============================================================================
!> Define a double global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_double(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  integer :: ier
#endif

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_double: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, "array", trim(array_name), 6, &
                        trim(array_name) // "/local_dim", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = local_dim
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_real8, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 8

  end subroutine define_adios_global_1d_generic_double

!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_1d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_2d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_3d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
subroutine define_adios_global_1d_double_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

end subroutine define_adios_global_1d_double_4d


!===============================================================================
!> Define a global ADIOS 1D double array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_double_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  real(kind=8), dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE double
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_double_5d


!===============================================================================
!> Define a integer global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_int(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  integer :: ier
#endif

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_int: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, "array", trim(array_name), 2, &
                        trim(array_name) // "/local_dim", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = local_dim
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  full_name = trim(array_name) // "/array"

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_integer4, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 4

  end subroutine define_adios_global_1d_generic_int

!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_1d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_2d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_3d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_4d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_int_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=4), dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE int
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_int_5d


!===============================================================================
!> Define a long integer global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
  subroutine define_adios_global_1d_generic_long(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  integer :: ier
#endif

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_long: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, "array", trim(array_name), adios_long, &
                        trim(array_name) // "/local_dim", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable will not be stored globally, instead will be retrieved by _inquire** function again
  ldim(1) = local_dim
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_integer8, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 8

  end subroutine define_adios_global_1d_generic_long

!===============================================================================
!> Define a global ADIOS 1D long array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_1d


!===============================================================================
!> Define a global ADIOS 1D long array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_2d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_3d


!===============================================================================
!> Define a global ADIOS 1D int array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_4d


!===============================================================================
!> Define a global ADIOS 1D long array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_long_5d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  integer(kind=8), dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE long
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_long_5d

!===============================================================================
!> Define a logical global array in ADIOS regardless of the array shape
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param local_dim The local dimension of the array.
subroutine define_adios_global_1d_generic_logical(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  integer :: ier
#endif

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_generic_logical: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  ! The Fortran standard does not specify how variables of LOGICAL type are
  ! represented, beyond requiring that LOGICAL variables of default kind
  ! have the same storage size as default INTEGER and REAL variables.
  ! Hence the 'adios_integer' (2) data type to store logical values
  call adios_define_var(adios_group, "array", trim(array_name), 2, &
                        trim(array_name) // "/local_dim", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  !
  !       also, adios2 has no adios2_get()/adios2_put() routine for logicals.
  !       we need to use integer array instead for storing/reading.
  ldim(1) = local_dim
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_integer4, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 4

end subroutine define_adios_global_1d_generic_logical

!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_1d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_2d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_2d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_3d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_3d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_4d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:,:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_4d


!===============================================================================
!> Define a global ADIOS 1D logical array and autoincrement the adios
!! group size.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param local_dim The local dimension of the array.
!! \param path The path where array name lie.
!! \param array_name The variable to be defined. This is actually the path for
!!                   ADIOS. The values are stored in array_name/array
!! \param var The variable to define. Used for type and shape inference.
!! \note This function define local, global and offset sizes as well as the
!!       array to store the values in.
  subroutine define_adios_global_1d_logical_5d(adios_group, group_size_inc,local_dim, path, array_name, var)

  implicit none
  ! Parameters
  logical, dimension(:,:,:,:,:), intent(in) :: var

#define VAR_TYPE logical
#include "adios_helpers_definitions_1d_generic.inc"
#undef VAR_TYPE

  end subroutine define_adios_global_1d_logical_5d

!===============================================================================

!string added
subroutine define_adios_global_1d_string_generic(adios_group, group_size_inc, array_name, local_dim)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  ! Variables
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array
  type(adios2_variable) :: v
  character(len=256) :: full_name
  integer :: ier
#endif

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_string_generic: ',trim(array_name))

  ! Define the dimensions of the array. local_dim used as a dummy
  ! variable to call the integer routine.
  call define_adios_global_dims_1d(adios_group, group_size_inc, trim(array_name), local_dim)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, "array", trim(array_name), 9, &
                        trim(array_name) // "/local_dim", &
                        trim(array_name) // "/global_dim", &
                        trim(array_name) // "/offset", var_id)

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: variable parameter v will not be stored globally, instead will be retrieved by _inquire** function again
  !
  !       we define the array as a global array (adios2_constant_dims = .true.) over all processes,
  !       but specify the corresponding sub-array dimension for each process using the dims

  ! local array dimensions
  ldim(1) = local_dim
  gdim(:) = sizeprocs_adios * ldim(:)
  offs(:) = myrank_adios * ldim(:)

  full_name = trim(array_name) // '/array'

  ! defines variable  (adios2_constant_dims = .true. argument for constant dimensions, won't change)
  call adios2_define_variable(v, adios_group, trim(full_name), &
                              adios2_type_string, 1, gdim, offs, ldim, adios2_constant_dims, ier)
  call check_adios_err(ier,"Error adios2 define variable "//trim(full_name)//" failed")

  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 1

end subroutine define_adios_global_1d_string_generic

!===============================================================================

subroutine define_adios_global_1d_string_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: path, array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in) :: var
  ! Local vars
  integer :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_global_1d_string_1d: ',trim(array_name))

  ! check
  if (len_trim(array_name) == 0) stop 'Error adios: invalid name in define_adios_global_1d_string_1d()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(array_name)
  else
    full_name = trim(path) // '/' // trim(array_name)
  endif

  ! debug
  !print *,"full name", trim(full_name),"local_dim:",local_dim

  call define_adios_global_1d_string_generic(adios_group, group_size_inc, full_name, local_dim)

  ! to avoid compiler warnings
  idummy = len(var)

end subroutine define_adios_global_1d_string_1d

!
!------------
!

subroutine  define_adios_local_1d_string_1d(adios_group, group_size_inc, local_dim, path, array_name, var)

  implicit none
  ! Parameters
#if defined(USE_ADIOS)
  integer(kind=8),  intent(in)    :: adios_group
#elif defined(USE_ADIOS2)
  type(adios2_io),  intent(in)     :: adios_group
#endif
  character(len=*), intent(in) :: path, array_name
  integer,          intent(in) :: local_dim
  integer(kind=8),  intent(inout) :: group_size_inc
  character(len=*), intent(in) :: var
  ! Local
#if defined(USE_ADIOS)
  integer(kind=8) :: var_id
#elif defined(USE_ADIOS2)
  type(adios2_variable) :: v
  integer :: ier
#endif
  integer :: idummy
  character(len=256) :: full_name

  TRACE_ADIOS_L2_ARG('define_adios_local_1d_string_1d: ',trim(array_name))

  ! check
  if (len_trim(array_name) == 0) stop 'Error adios: invalid name in define_adios_local_1d_string_1d()'

  ! sets full variable name
  if (len_trim(path) == 0) then
    full_name = trim(array_name)
  else
    full_name = trim(path) // '/' // trim(array_name)
  endif

  !debug
  !print *,"in define local: full_name:", trim(full_name)

#if defined(USE_ADIOS)
  ! ADIOS 1
  call adios_define_var(adios_group, trim(array_name), trim(path), 9, '', '', '', var_id )

#elif defined(USE_ADIOS2)
  ! ADIOS 2
  ! note: we won't store the variable object, but get it back by adios2_inquire_** for writing
  ! defines global variable
  call adios2_define_variable(v, adios_group, trim(full_name), adios2_type_string, ier)
  call check_adios_err(ier,"Error adios2 could not define parameter: "//trim(full_name))
  if (.not. v%valid) stop 'Error adios2 defined variable is invalid'

#endif

  group_size_inc = group_size_inc + local_dim * 1

  ! to avoid compiler warnings
  idummy = len(var)

end subroutine define_adios_local_1d_string_1d

end module adios_helpers_definitions_mod
