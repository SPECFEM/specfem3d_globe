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

!-------------------------------------------------------------------------------
!> \file save_forward_arrays_adios.F90
!! \brief Save forward arrays with the help of the ADIOS library.
!! \author MPBL
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> \brief Write intermediate forward arrays in an ADIOS file.
!!
!! This subroutine is only used when NUMBER_OF_RUNS > 1 and
!! NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS.
subroutine save_intermediate_forward_arrays_adios()
  ! External imports
  use mpi
  use adios_write_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none
  ! Local parameters
  integer :: sizeprocs, comm, ierr
  character(len=150) :: outputname
  integer(kind=8) :: group_size_inc
  integer :: local_dim, global_dim, offset
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  outputname = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"
  call world_size(sizeprocs) ! TODO keep it in parameters
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  group_size_inc = 0
  call adios_declare_group(adios_group, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", &
      "", 1, adios_err)
!  call check_adios_err(myrank,adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)
!  call check_adios_err(myrank,adios_err)

  ! Define ADIOS variables
  call define_common_forward_arrays_adios(adios_group, group_size_inc)
  call define_rotation_forward_arrays_adios(adios_group, group_size_inc)
  call define_attenuation_forward_arrays_adios(adios_group, group_size_inc)

  ! Open an ADIOS handler to the restart file.
  call adios_open (adios_handle, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", &
      outputname, "w", comm, adios_err);
!  call check_adios_err(myrank,adios_err)
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)
!  call check_adios_err(myrank,adios_err)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios(adios_handle, sizeprocs)
  call write_rotation_forward_arrays_adios(adios_handle, sizeprocs)
  call write_attenuation_forward_arrays_adios(adios_handle, sizeprocs)
  ! Reset the path to its original value to avoid bugs.
  call adios_set_path (adios_handle, "", adios_err)
!  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_close(adios_handle, adios_err)
!  call check_adios_err(myrank,adios_err)
end subroutine save_intermediate_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simualtions when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios() execpt than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.
subroutine save_forward_arrays_adios()
  ! External imports
  use mpi
  use adios_write_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none
  ! Local parameters
  integer :: sizeprocs, comm, ierr
  character(len=150) :: outputname
  integer(kind=8) :: group_size_inc
  integer :: local_dim, global_dim, offset
!  integer, parameter :: num_arrays = 9 ! TODO correct number
!  character(len=256), dimension(num_arrays) :: local_dims1, local_dims2, &
!      global_dims1, global_dims2, offsets1, offsets2, array_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  outputname = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"
  call world_size(sizeprocs)
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  group_size_inc = 0

  call adios_declare_group(adios_group, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", &
      "", 1, adios_err)
!  call check_adios_err(myrank,adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)
!  call check_adios_err(myrank,adios_err)

  ! Define ADIOS variables
  call define_common_forward_arrays_adios(adios_group, group_size_inc)
  ! conditional definition of vars seem to mess with the group size,
  ! even if the variables are conditionnaly written.
!  if (ROTATION_VAL) then
    call define_rotation_forward_arrays_adios(adios_group, group_size_inc)
!  endif
!  if (ATTENUATION_VAL) then
    call define_attenuation_forward_arrays_adios(adios_group, group_size_inc)
!  endif

  ! Open an ADIOS handler to the restart file.
  call adios_open (adios_handle, "SPECFEM3D_GLOBE_FORWARD_ARRAYS", &
      outputname, "w", comm, adios_err);
!  call check_adios_err(myrank,adios_err)
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)
!  call check_adios_err(myrank,adios_err)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios(adios_handle, sizeprocs)
  if (ROTATION_VAL) then
      call write_rotation_forward_arrays_adios(adios_handle, sizeprocs)
  endif
  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios(adios_handle, sizeprocs)
  endif
  ! Reset the path to its original value to avoid bugs.
  call adios_set_path (adios_handle, "", adios_err)
!  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_close(adios_handle, adios_err)
!  call check_adios_err(myrank,adios_err)
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

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  local_dim = NDIM * NGLOB_CRUST_MANTLE
  call define_adios_global_real_1d_array(adios_group, "displ_crust_mantle", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "veloc_crust_mantle", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "accel_crust_mantle", &
      local_dim, group_size_inc)

  local_dim = NDIM * NGLOB_INNER_CORE
  call define_adios_global_real_1d_array(adios_group, "displ_inner_core", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "veloc_inner_core", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "accel_inner_core", &
      local_dim, group_size_inc)

  local_dim = NGLOB_OUTER_CORE
  call define_adios_global_real_1d_array(adios_group, "displ_outer_core", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "veloc_outer_core", &
      local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, "accel_outer_core", &
      local_dim, group_size_inc)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_xx_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_yy_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_xy_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_xz_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_yz_crust_mantle", local_dim, group_size_inc)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_xx_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_yy_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_xy_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_xz_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "epsilondev_yz_inner_core", local_dim, group_size_inc)
end subroutine define_common_forward_arrays_adios

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

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  call define_adios_global_real_1d_array(adios_group, &
      "A_array_rotation", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "B_array_rotation", local_dim, group_size_inc)
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

  implicit none

  integer(kind=8), intent(in) :: adios_group
  integer(kind=8), intent(inout) :: group_size_inc

  integer :: local_dim

  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUAT
  call define_adios_global_real_1d_array(adios_group, &
      "R_xx_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_yy_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_xy_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_xz_crust_mantle", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_yz_crust_mantle", local_dim, group_size_inc)

  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  call define_adios_global_real_1d_array(adios_group, &
      "R_xx_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_yy_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_xy_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_xz_inner_core", local_dim, group_size_inc)
  call define_adios_global_real_1d_array(adios_group, &
      "R_yz_inner_core", local_dim, group_size_inc)
end subroutine define_attenuation_forward_arrays_adios

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.
!! \param adios_handle The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writting
subroutine write_common_forward_arrays_adios(adios_handle, sizeprocs)
  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs

  integer :: local_dim, adios_err

  local_dim = NDIM * NGLOB_CRUST_MANTLE
  call adios_set_path (adios_handle, "displ_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", displ_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "veloc_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", veloc_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "accel_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", accel_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM * NGLOB_INNER_CORE
  call adios_set_path (adios_handle, "displ_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", displ_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "veloc_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", veloc_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "accel_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", accel_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_set_path (adios_handle, "", adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLOB_OUTER_CORE
  call adios_set_path (adios_handle, "displ_outer_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", displ_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "veloc_outer_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", veloc_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "accel_outer_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", accel_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  call adios_set_path (adios_handle, "epsilondev_xx_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_yy_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_xy_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_xz_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_yz_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  call adios_set_path (adios_handle, "epsilondev_xx_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_yy_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_xy_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_xz_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "epsilondev_yz_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", epsilondev_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
end subroutine write_common_forward_arrays_adios

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ROTATION is true.
!! \param adios_handle The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writting
subroutine write_rotation_forward_arrays_adios(adios_handle, sizeprocs)
  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs

  integer :: local_dim, adios_err

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  call adios_set_path (adios_handle, "A_array_rotation", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", A_array_rotation, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "B_array_rotation", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", B_array_rotation, adios_err)
  call check_adios_err(myrank,adios_err)
end subroutine write_rotation_forward_arrays_adios

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ATTENUATION
!!  is true.
!! \param adios_handle The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writting
subroutine write_attenuation_forward_arrays_adios(adios_handle, sizeprocs)
  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs

  integer :: local_dim, adios_err

  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUAT
  call adios_set_path (adios_handle, "R_xx_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_yy_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_xy_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_xz_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_yz_crust_mantle", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  call adios_set_path (adios_handle, "R_xx_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_yy_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_xy_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_xz_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_set_path (adios_handle, "R_yz_inner_core", adios_err)
  call check_adios_err(myrank,adios_err)
  call write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  call adios_write(adios_handle, "array", R_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
end subroutine write_attenuation_forward_arrays_adios

!-------------------------------------------------------------------------------
!> Write local, global and offset dimensions to ADIOS
!! \param adios_handle Handle to the adios file
!! \param local_dim Number of elements to be written by one process
!! \param sizeprocs Number of MPI processes
subroutine write_1D_global_array_adios_dims(adios_handle, local_dim, sizeprocs)
  use adios_write_mod
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: sizeprocs, local_dim

  integer :: adios_err

  call adios_write(adios_handle, "local_dim", local_dim, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(adios_handle, "global_dim", local_dim*sizeprocs, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_write(adios_handle, "offset", local_dim*myrank, adios_err)
  call check_adios_err(myrank,adios_err)
end subroutine write_1D_global_array_adios_dims

