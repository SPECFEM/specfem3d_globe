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
  use adios_read_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err

  implicit none
  ! Local parameters
  integer :: comm
  character(len=MAX_STRING_LEN) :: file_name, path_to_add
  integer :: local_dim
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_handle, sel
  integer(kind=8), dimension(1) :: start, count


  file_name = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    file_name = path_to_add(1:len_trim(path_to_add))//file_name(1:len_trim(file_name))
  endif

  call world_duplicate(comm)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_open_file (adios_handle, file_name, 0, comm, adios_err)
  if (adios_err /= 0) then
    print*,'Error rank ',myrank,' opening adios file: ',trim(file_name)
    call check_adios_err(myrank,adios_err)
  endif

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_crust_mantle/array", 0, 1, displ_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_crust_mantle/array", 0, 1, veloc_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_crust_mantle/array", 0, 1, accel_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_inner_core/array", 0, 1, displ_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_inner_core/array", 0, 1, veloc_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_inner_core/array", 0, 1, accel_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_outer_core/array", 0, 1, displ_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_outer_core/array", 0, 1, veloc_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_outer_core/array", 0, 1, accel_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! strains crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_crust_mantle/array", 0, 1, epsilondev_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_crust_mantle/array", 0, 1, epsilondev_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_crust_mantle/array", 0, 1, epsilondev_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_crust_mantle/array", 0, 1, epsilondev_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_crust_mantle/array", 0, 1, epsilondev_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! strains inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_inner_core/array", 0, 1, epsilondev_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_inner_core/array", 0, 1, epsilondev_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_inner_core/array", 0, 1, epsilondev_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_inner_core/array", 0, 1, epsilondev_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_inner_core/array", 0, 1, epsilondev_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! rotation
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "A_array_rotation/array", 0, 1, A_array_rotation, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "B_array_rotation/array", 0, 1, B_array_rotation, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! attenuation memory variables crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "R_xx_crust_mantle/array", 0, 1, R_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yy_crust_mantle/array", 0, 1, R_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xy_crust_mantle/array", 0, 1, R_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xz_crust_mantle/array", 0, 1, R_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yz_crust_mantle/array", 0, 1, R_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! attenuation memory variables inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "R_xx_inner_core/array", 0, 1, R_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yy_inner_core/array", 0, 1, R_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xy_inner_core/array", 0, 1, R_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_xz_inner_core/array", 0, 1, R_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "R_yz_inner_core/array", 0, 1, R_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)

  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call synchronize_all_comm(comm)

  end subroutine read_intermediate_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios() and read_forward_arrays_adios() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays

  subroutine read_forward_arrays_adios()

  ! External imports
  use adios_read_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err

  implicit none
  ! Local parameters
  integer :: comm
  character(len=MAX_STRING_LEN) :: file_name, path_to_add
  integer :: local_dim
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_handle, sel
  integer(kind=8), dimension(1) :: start, count

  file_name = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    file_name = path_to_add(1:len_trim(path_to_add))//file_name(1:len_trim(file_name))
  endif

  call world_duplicate(comm)

  ! opens adios file
  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_open_file (adios_handle, file_name, 0, comm, adios_err)
  if (adios_err /= 0) then
    print*,'Error rank ',myrank,' opening adios file: ',trim(file_name)
    call check_adios_err(myrank,adios_err)
  endif

  ! reads in arrays
  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_crust_mantle/array", 0, 1, b_displ_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_crust_mantle/array", 0, 1, b_veloc_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_crust_mantle/array", 0, 1, b_accel_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_inner_core/array", 0, 1, b_displ_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_inner_core/array", 0, 1, b_veloc_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_inner_core/array", 0, 1, b_accel_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_outer_core/array", 0, 1, b_displ_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_outer_core/array", 0, 1, b_veloc_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_outer_core/array", 0, 1, b_accel_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! strains crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_crust_mantle/array", 0, 1, b_epsilondev_xx_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_crust_mantle/array", 0, 1, b_epsilondev_yy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_crust_mantle/array", 0, 1, b_epsilondev_xy_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_crust_mantle/array", 0, 1, b_epsilondev_xz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_crust_mantle/array", 0, 1, b_epsilondev_yz_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! strains inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xx_inner_core/array", 0, 1, b_epsilondev_xx_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yy_inner_core/array", 0, 1, b_epsilondev_yy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xy_inner_core/array", 0, 1, b_epsilondev_xy_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_xz_inner_core/array", 0, 1, b_epsilondev_xz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "epsilondev_yz_inner_core/array", 0, 1, b_epsilondev_yz_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! rotation
  if (ROTATION_VAL) then
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "A_array_rotation/array", 0, 1, b_A_array_rotation, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "B_array_rotation/array", 0, 1, b_B_array_rotation, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! attenuation memory variables
  if (ATTENUATION_VAL) then
    ! crust/mantle
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "R_xx_crust_mantle/array", 0, 1, b_R_xx_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yy_crust_mantle/array", 0, 1, b_R_yy_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xy_crust_mantle/array", 0, 1, b_R_xy_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xz_crust_mantle/array", 0, 1, b_R_xz_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yz_crust_mantle/array", 0, 1, b_R_yz_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    ! inner core
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "R_xx_inner_core/array", 0, 1, b_R_xx_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yy_inner_core/array", 0, 1, b_R_yy_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xy_inner_core/array", 0, 1, b_R_xy_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xz_inner_core/array", 0, 1, b_R_xz_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yz_inner_core/array", 0, 1, b_R_yz_inner_core, adios_err)
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

  call synchronize_all_comm(comm)

  end subroutine read_forward_arrays_adios


!-------------------------------------------------------------------------------
!> \brief Read forward arrays for undo attenuation from an ADIOS file.

  subroutine read_forward_arrays_undoatt_adios(iteration_on_subset_tmp)

  ! External imports
  use adios_read_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err

  implicit none
  ! Arguments
  integer, intent(in) :: iteration_on_subset_tmp
  ! Local parameters
  integer :: comm
  character(len=MAX_STRING_LEN) :: file_name, path_to_add
  integer :: local_dim
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_handle, sel
  integer(kind=8), dimension(1) :: start, count


  write(file_name,'(a,a,i6.6,a)') trim(LOCAL_TMP_PATH), '/save_frame_at',iteration_on_subset_tmp,'.bp'
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    file_name = path_to_add(1:len_trim(path_to_add))//file_name(1:len_trim(file_name))
  endif

  call world_duplicate(comm)

  ! opens adios file
  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_open_file (adios_handle, file_name, 0, comm, adios_err)
  if (adios_err /= 0) then
    print*,'Error rank ',myrank,' opening adios file: ',trim(file_name)
    call check_adios_err(myrank,adios_err)
  endif

  ! reads in arrays
  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_crust_mantle/array", 0, 1, b_displ_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_crust_mantle/array", 0, 1, b_veloc_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_crust_mantle/array", 0, 1, b_accel_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_inner_core/array", 0, 1, b_displ_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_inner_core/array", 0, 1, b_veloc_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_inner_core/array", 0, 1, b_accel_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "displ_outer_core/array", 0, 1, b_displ_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "veloc_outer_core/array", 0, 1, b_veloc_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "accel_outer_core/array", 0, 1, b_accel_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! rotation
  if (ROTATION_VAL) then
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "A_array_rotation/array", 0, 1, b_A_array_rotation, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "B_array_rotation/array", 0, 1, b_B_array_rotation, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! attenuation memory variables
  if (ATTENUATION_VAL) then
    ! crust/mantle
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "R_xx_crust_mantle/array", 0, 1, b_R_xx_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yy_crust_mantle/array", 0, 1, b_R_yy_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xy_crust_mantle/array", 0, 1, b_R_xy_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xz_crust_mantle/array", 0, 1, b_R_xz_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yz_crust_mantle/array", 0, 1, b_R_yz_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    ! inner core
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "R_xx_inner_core/array", 0, 1, b_R_xx_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yy_inner_core/array", 0, 1, b_R_yy_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xy_inner_core/array", 0, 1, b_R_xy_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_xz_inner_core/array", 0, 1, b_R_xz_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "R_yz_inner_core/array", 0, 1, b_R_yz_inner_core, adios_err)
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

  call synchronize_all_comm(comm)

  end subroutine read_forward_arrays_undoatt_adios
