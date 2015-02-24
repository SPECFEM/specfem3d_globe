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

!-------------------------------------------------------------------------------
!> \file save_forward_arrays_adios.F90
!! \brief Save forward arrays with the help of the ADIOS library.
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"


!-------------------------------------------------------------------------------
!> \brief Write intermediate forward arrays in an ADIOS file.
!!
!! This subroutine is only used when NUMBER_OF_RUNS > 1 and
!! NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS.

  subroutine save_intermediate_forward_arrays_adios()

  ! External imports
  use adios_write_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err

  implicit none
  ! Local parameters
  integer :: comm
  character(len=MAX_STRING_LEN) :: outputname
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle
  integer(kind=8)         :: adios_totalsize

  outputname = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"

  call world_duplicate(comm)

  group_size_inc = 0

  call adios_declare_group(adios_group, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  ! Define ADIOS variables
  call define_common_forward_arrays_adios(adios_group, group_size_inc)
  call define_epsilon_forward_arrays_adios(adios_group, group_size_inc)
  call define_rotation_forward_arrays_adios(adios_group, group_size_inc)
  call define_attenuation_forward_arrays_adios(adios_group, group_size_inc)

  ! Open an ADIOS handler to the restart file.
  call adios_open (adios_handle, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", outputname, "w", comm, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed for SPECFEM3D_GLOBE_FORWARD_ARRAYS'

  call adios_group_size (adios_handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios(adios_handle)
  call write_epsilon_forward_arrays_adios(adios_handle)
  call write_rotation_forward_arrays_adios(adios_handle)
  call write_attenuation_forward_arrays_adios(adios_handle)
  ! Reset the path to its original value to avoid bugs.
  call adios_set_path (adios_handle, "", adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_close(adios_handle, adios_err)

  end subroutine save_intermediate_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_adios()

  ! External imports
  use adios_write_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err


  implicit none
  ! Local parameters
  integer :: comm
  character(len=MAX_STRING_LEN) :: outputname
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle
  integer(kind=8)         :: adios_totalsize

  outputname = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"

  call world_duplicate(comm)

  group_size_inc = 0

  call adios_declare_group(adios_group, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", "", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD, "", "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  ! Define ADIOS variables
  call define_common_forward_arrays_adios(adios_group, group_size_inc)
  call define_epsilon_forward_arrays_adios(adios_group, group_size_inc)
  ! TODO check following:
  ! conditional definition of vars seem to mess with the group size,
  ! even if the variables are conditionally written.
  !  if (ROTATION_VAL) then
    call define_rotation_forward_arrays_adios(adios_group, group_size_inc)
  !  endif
  !  if (ATTENUATION_VAL) then
    call define_attenuation_forward_arrays_adios(adios_group, group_size_inc)
  !  endif

  ! Open an ADIOS handler to the restart file.
  call adios_open (adios_handle, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", outputname, "w", comm, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed for SPECFEM3D_GLOBE_FORWARD_ARRAYS'

  call adios_group_size (adios_handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios(adios_handle)

  call write_epsilon_forward_arrays_adios(adios_handle)

  if (ROTATION_VAL) then
      call write_rotation_forward_arrays_adios(adios_handle)
  endif

  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios(adios_handle)
  endif

  ! Reset the path to its original value to avoid bugs.
  call adios_set_path (adios_handle, "", adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_close(adios_handle, adios_err)

  end subroutine save_forward_arrays_adios

!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are always dumped.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_common_forward_arrays_adios(adios_group, group_size_inc)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(displ_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(veloc_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(accel_crust_mantle))

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(displ_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(veloc_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(accel_inner_core))

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(displ_outer_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(veloc_outer_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(accel_outer_core))

  ! strains
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_xx_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_yy_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_xy_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_xz_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_yz_crust_mantle))


  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_xx_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_yy_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_xy_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_xz_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_yz_inner_core))

  end subroutine define_common_forward_arrays_adios


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are always dumped
!! except for undo attenuation.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_epsilon_forward_arrays_adios(adios_group, group_size_inc)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  ! strains
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_xx_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_yy_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_xy_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_xz_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", &
                                   STRINGIFY_VAR(epsilondev_yz_crust_mantle))


  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_xx_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_yy_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_xy_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_xz_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(epsilondev_yz_inner_core))

  end subroutine define_epsilon_forward_arrays_adios


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are dumped if ROTATION is true.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_rotation_forward_arrays_adios(adios_group, group_size_inc)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(A_array_rotation))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(B_array_rotation))

  end subroutine define_rotation_forward_arrays_adios

!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are dumped if ATTENUATION is true.
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_attenuation_forward_arrays_adios(adios_group, group_size_inc)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  ! attenuation memory variables
  ! crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_xx_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_yy_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_xy_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_xz_crust_mantle))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_yz_crust_mantle))

  ! inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_xx_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_yy_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_xy_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_xz_inner_core))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, "", STRINGIFY_VAR(R_yz_inner_core))

  end subroutine define_attenuation_forward_arrays_adios

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.
!! \param adios_handle The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writing

  subroutine write_common_forward_arrays_adios(adios_handle)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_handle

  integer :: local_dim !, adios_err

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(displ_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(veloc_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(accel_crust_mantle))

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(displ_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(veloc_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(accel_inner_core))

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(displ_outer_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(veloc_outer_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(accel_outer_core))

  end subroutine write_common_forward_arrays_adios


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.
!! \param adios_handle The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writing

  subroutine write_epsilon_forward_arrays_adios(adios_handle)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_handle

  integer :: local_dim !, adios_err

  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ! strains
  ! crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xx_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(epsilondev_yy_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xy_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xz_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, &
                                   STRINGIFY_VAR(epsilondev_yz_crust_mantle))

  ! inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(epsilondev_xx_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(epsilondev_yy_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(epsilondev_xy_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(epsilondev_xz_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(epsilondev_yz_inner_core))

  end subroutine write_epsilon_forward_arrays_adios


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ROTATION is true.
!! \param adios_handle The handle to the adios bp file

  subroutine write_rotation_forward_arrays_adios(adios_handle)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_handle

  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(A_array_rotation))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(B_array_rotation))

  end subroutine write_rotation_forward_arrays_adios

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ATTENUATION
!!  is true.
!! \param adios_handle The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writing

  subroutine write_attenuation_forward_arrays_adios(adios_handle)

  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod

  implicit none

  integer(kind=8), intent(in) :: adios_handle

  integer :: local_dim !, adios_err
  integer :: sizeprocs

  ! number of MPI processes
  call world_size(sizeprocs)

  ! attenuation memory variables
  ! crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_xx_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_yy_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_xy_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_xz_crust_mantle))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_yz_crust_mantle))

  ! inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_xx_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_yy_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_xy_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_xz_inner_core))
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, local_dim, STRINGIFY_VAR(R_yz_inner_core))

  end subroutine write_attenuation_forward_arrays_adios

!!-------------------------------------------------------------------------------
!!> Write local, global and offset dimensions to ADIOS
!!! \param adios_handle Handle to the adios file
!!! \param local_dim Number of elements to be written by one process
!!! \param sizeprocs Number of MPI processes
!subroutine write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  !use adios_write_mod
  !use specfem_par, only: myrank

  !use adios_helpers_mod, only: check_adios_err

  !implicit none

  !integer(kind=8), intent(in) :: adios_handle
  !integer, intent(in) :: sizeprocs, local_dim

  !integer :: adios_err

  !call adios_write(adios_handle, "local_dim", local_dim, adios_err)
  !call check_adios_err(myrank,adios_err)
  !call adios_write(adios_handle, "global_dim", local_dim*sizeprocs, adios_err)
  !call check_adios_err(myrank,adios_err)
  !call adios_write(adios_handle, "offset", local_dim*myrank, adios_err)
  !call check_adios_err(myrank,adios_err)
!end subroutine write_1D_global_array_adios_dims

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_undoatt_adios()

  ! External imports
  use adios_write_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err


  implicit none
  ! Local parameters
  integer :: comm, iteration_on_subset_tmp
  character(len=MAX_STRING_LEN) :: outputname
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=MAX_STRING_LEN) :: group_name
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle
  integer(kind=8)         :: adios_totalsize

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  write(outputname,'(a, a, i6.6, a)') trim(LOCAL_PATH), '/save_frame_at', iteration_on_subset_tmp,'.bp'
  write(group_name, '(a, i6)') "SPECFEM3D_GLOBE_FORWARD_ARRAYS", iteration_on_subset_tmp

  call world_duplicate(comm)

  group_size_inc = 0

  call adios_declare_group(adios_group, group_name, "iter", 1, adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  call adios_select_method(adios_group, ADIOS_TRANSPORT_METHOD_UNDO_ATT, &
                           ADIOS_METHOD_PARAMS_UNDO_ATT, "", adios_err)
  ! note: return codes for this function have been fixed for ADIOS versions >= 1.6
  !call check_adios_err(myrank,adios_err)

  ! Define ADIOS variables
  call define_common_forward_arrays_adios(adios_group, group_size_inc)
  ! TODO check following:
  ! conditional definition of vars seem to mess with the group size,
  ! even if the variables are conditionally written.
  !if (ROTATION_VAL) then
  call define_rotation_forward_arrays_adios(adios_group, group_size_inc)
  !endif
  !if (ATTENUATION_VAL) then
  call define_attenuation_forward_arrays_adios(adios_group, group_size_inc)
  !endif

  ! Open an ADIOS handler to the restart file.
  call adios_open (adios_handle, group_name, outputname, "w", comm, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_open() routine failed'

  call adios_group_size (adios_handle, group_size_inc, adios_totalsize, adios_err)
  if (adios_err /= 0 ) stop 'Error calling adios_group_size() routine failed'

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios(adios_handle)

  if (ROTATION_VAL) then
      call write_rotation_forward_arrays_adios(adios_handle)
  endif

  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios(adios_handle)
  endif

  ! Reset the path to its original value to avoid bugs.
  call adios_set_path (adios_handle, "", adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_close(adios_handle, adios_err)

end subroutine save_forward_arrays_undoatt_adios
