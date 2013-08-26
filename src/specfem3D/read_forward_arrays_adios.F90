!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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
!> \file read_forward_arrays_adios.F90
!! \brief Read saved forward arrays with the help of the ADIOS library.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios()
!!       and read_forward_arrays_adios() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays
subroutine read_intermediate_forward_arrays_adios()
  ! External imports
  use mpi
  use adios_read_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none
  ! Local parameters
  integer :: sizeprocs, comm, ierr
  character(len=150) :: file_name
  integer(kind=8) :: group_size_inc
  integer :: local_dim, global_dim, offset
!  integer, parameter :: num_arrays = 9 ! TODO correct number
!  character(len=256), dimension(num_arrays) :: local_dims1, local_dims2, &
!      global_dims1, global_dims2, offsets1, offsets2, array_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid, sel
  integer(kind=8)         :: adios_groupsize, adios_totalsize
  integer :: vars_count, attrs_count, current_step, last_step, vsteps
  character(len=128), dimension(:), allocatable :: adios_names
  integer(kind=8), dimension(1) :: start, count


  file_name = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"
  call world_size(sizeprocs)
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)


  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_crust_mantle/array", 0, 1, &
      displ_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_crust_mantle/array", 0, 1, &
      veloc_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_crust_mantle/array", 0, 1, &
      accel_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_inner_core/array", 0, 1, &
      displ_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_inner_core/array", 0, 1, &
      veloc_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_inner_core/array", 0, 1, &
      accel_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_outer_core/array", 0, 1, &
      displ_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_outer_core/array", 0, 1, &
      veloc_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_outer_core/array", 0, 1, &
      accel_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_crust_mantle/array",&
      0, 1, epsilondev_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_crust_mantle/array",&
      0, 1, epsilondev_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_crust_mantle/array",&
      0, 1, epsilondev_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_crust_mantle/array",&
      0, 1, epsilondev_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_crust_mantle/array",&
      0, 1, epsilondev_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_inner_core/array",&
      0, 1, epsilondev_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_inner_core/array",&
      0, 1, epsilondev_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_inner_core/array",&
      0, 1, epsilondev_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_inner_core/array",&
      0, 1, epsilondev_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_inner_core/array",&
      0, 1, epsilondev_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "A_array_rotation/array", 0, 1, &
      A_array_rotation, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "B_array_rotation/array", 0, 1, &
      B_array_rotation, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "R_xx_crust_mantle/array", 0, 1, &
      R_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yy_crust_mantle/array", 0, 1, &
      R_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xy_crust_mantle/array", 0, 1, &
      R_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xz_crust_mantle/array", 0, 1, &
      R_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yz_crust_mantle/array", 0, 1, &
      R_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "R_xx_inner_core/array", 0, 1, &
      R_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yy_inner_core/array", 0, 1, &
      R_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xy_inner_core/array", 0, 1, &
      R_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xz_inner_core/array", 0, 1, &
      R_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yz_inner_core/array", 0, 1, &
      R_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call sync_all()

end subroutine read_intermediate_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios()
!!       and read_forward_arrays_adios() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays
subroutine read_forward_arrays_adios()
  ! External imports
  use mpi
  use adios_read_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none
  ! Local parameters
  integer :: sizeprocs, comm, ierr
  character(len=150) :: file_name
  integer(kind=8) :: group_size_inc
  integer :: local_dim, global_dim, offset
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid, sel
  integer(kind=8)         :: adios_groupsize, adios_totalsize
  integer :: vars_count, attrs_count, current_step, last_step, vsteps
  character(len=128), dimension(:), allocatable :: adios_names
  integer(kind=8), dimension(1) :: start, count


  file_name = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"
  call world_size(sizeprocs)
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)


  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_crust_mantle/array", 0, 1, &
      b_displ_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_crust_mantle/array", 0, 1, &
      b_veloc_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_crust_mantle/array", 0, 1, &
      b_accel_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_inner_core/array", 0, 1, &
      b_displ_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_inner_core/array", 0, 1, &
      b_veloc_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_inner_core/array", 0, 1, &
      b_accel_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_outer_core/array", 0, 1, &
      b_displ_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_outer_core/array", 0, 1, &
      b_veloc_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_outer_core/array", 0, 1, &
      b_accel_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_crust_mantle/array",&
      0, 1, b_epsilondev_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_crust_mantle/array",&
      0, 1, b_epsilondev_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_crust_mantle/array",&
      0, 1, b_epsilondev_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_crust_mantle/array",&
      0, 1, b_epsilondev_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_crust_mantle/array",&
      0, 1, b_epsilondev_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_inner_core/array",&
      0, 1, b_epsilondev_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_inner_core/array",&
      0, 1, b_epsilondev_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_inner_core/array",&
      0, 1, b_epsilondev_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_inner_core/array",&
      0, 1, b_epsilondev_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_inner_core/array",&
      0, 1, b_epsilondev_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  if (ROTATION_VAL) then
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "A_array_rotation/array", 0, 1, &
        b_A_array_rotation, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "B_array_rotation/array", 0, 1, &
        b_B_array_rotation, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  if (ATTENUATION_VAL) then
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "R_xx_crust_mantle/array", 0, 1, &
        b_R_xx_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yy_crust_mantle/array", 0, 1, &
        b_R_yy_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xy_crust_mantle/array", 0, 1, &
        b_R_xy_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xz_crust_mantle/array", 0, 1, &
        b_R_xz_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yz_crust_mantle/array", 0, 1, &
        b_R_yz_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "R_xx_inner_core/array", 0, 1, &
        b_R_xx_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yy_inner_core/array", 0, 1, &
        b_R_yy_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xy_inner_core/array", 0, 1, &
        b_R_xy_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xz_inner_core/array", 0, 1, &
        b_R_xz_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yz_inner_core/array", 0, 1, &
        b_R_yz_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call sync_all()

end subroutine read_forward_arrays_adios
