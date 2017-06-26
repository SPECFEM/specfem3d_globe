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
!> \brief Read adios boundary arrays created by the mesher
!!        (file: regX_boundary.bp)
subroutine read_mesh_databases_coupling_adios()

! to couple mantle with outer core

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use manager_adios

  implicit none

  ! local parameters
  integer :: njunk1,njunk2,njunk3
  character(len=MAX_STRING_LEN) :: file_name
  integer :: local_dim
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  character(len=128)      :: region_name

  file_name = trim(LOCAL_PATH) // "/boundary.bp"

  ! opens adios file
  call open_file_adios_read(file_name)

  ! crust and mantle
  write(region_name,"('reg',i1, '/')") IREGION_CRUST_MANTLE

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_xmin", 0, 1, &
     nspec2D_xmin_crust_mantle, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_xmax", 0, 1, &
     nspec2D_xmax_crust_mantle, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_ymin", 0, 1, &
     nspec2D_ymin_crust_mantle, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_ymax", 0, 1, &
     nspec2D_ymax_crust_mantle, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "NSPEC2D_BOTTOM", 0, 1, &
     njunk1, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "NSPEC2D_TOP", 0, 1, &
     njunk2, adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! boundary elements

  local_dim = NSPEC2DMAX_XMIN_XMAX_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_xmin/array", 0, 1, &
    ibelm_xmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_xmax/array", 0, 1, &
    ibelm_xmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2DMAX_YMIN_YMAX_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_ymin/array", 0, 1, &
    ibelm_ymin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_ymax/array", 0, 1, &
    ibelm_ymax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2D_BOTTOM_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_bottom/array", 0, 1, &
    ibelm_bottom_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2D_TOP_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_top/array", 0, 1, &
    ibelm_top_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_xmin/array", 0, 1, &
    normal_xmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_xmax/array", 0, 1, &
    normal_xmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_ymin/array", 0, 1, &
    normal_ymin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_ymax/array", 0, 1, &
    normal_ymax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_bottom/array", 0, 1, &
    normal_bottom_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_top/array", 0, 1, &
    normal_top_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_xmin/array", 0, 1, &
    jacobian2D_xmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_xmax/array", 0, 1, &
    jacobian2D_xmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_ymin/array", 0, 1, &
    jacobian2D_ymin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_ymax/array", 0, 1, &
    jacobian2D_ymax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_bottom/array", 0, 1, &
    jacobian2D_bottom_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX*NGLLY*NSPEC2D_TOP_CM
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_top/array", 0, 1, &
    jacobian2D_top_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)


  ! read parameters to couple fluid and solid regions
  !
  ! outer core

  write(region_name,"('reg',i1, '/')") IREGION_OUTER_CORE

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_xmin", 0, 1, &
     nspec2D_xmin_outer_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_xmax", 0, 1, &
     nspec2D_xmax_outer_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_ymin", 0, 1, &
     nspec2D_ymin_outer_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_ymax", 0, 1, &
     nspec2D_ymax_outer_core, adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! boundary elements

  local_dim = NSPEC2DMAX_XMIN_XMAX_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_xmin/array", 0, 1, &
    ibelm_xmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_xmax/array", 0, 1, &
    ibelm_xmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2DMAX_YMIN_YMAX_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_ymin/array", 0, 1, &
    ibelm_ymin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_ymax/array", 0, 1, &
    ibelm_ymax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2D_BOTTOM_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_bottom/array", 0, 1, &
    ibelm_bottom_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2D_TOP_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_top/array", 0, 1, &
    ibelm_top_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! normals

  local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_xmin/array", 0, 1, &
    normal_xmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_xmax/array", 0, 1, &
    normal_xmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_ymin/array", 0, 1, &
    normal_ymin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_ymax/array", 0, 1, &
    normal_ymax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_bottom/array", 0, 1, &
    normal_bottom_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_top/array", 0, 1, &
    normal_top_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Jacobians

  local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_xmin/array", 0, 1, &
    jacobian2D_xmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_xmax/array", 0, 1, &
    jacobian2D_xmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_ymin/array", 0, 1, &
    jacobian2D_ymin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_ymax/array", 0, 1, &
    jacobian2D_ymax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_bottom/array", 0, 1, &
    jacobian2D_bottom_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX*NGLLY*NSPEC2D_TOP_OC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "jacobian2D_top/array", 0, 1, &
    jacobian2D_top_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)


  ! inner core

  write(region_name,"('reg',i1, '/')") IREGION_INNER_CORE

  ! number of elements
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_xmin", 0, 1, &
     nspec2D_xmin_inner_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_xmax", 0, 1, &
     nspec2D_xmax_inner_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_ymin", 0, 1, &
     nspec2D_ymin_inner_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec2D_ymax", 0, 1, &
     nspec2D_ymax_inner_core, adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! boundary elements

  local_dim = NSPEC2DMAX_XMIN_XMAX_IC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_xmin/array", 0, 1, &
    ibelm_xmin_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_xmax/array", 0, 1, &
    ibelm_xmax_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2DMAX_YMIN_YMAX_IC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_ymin/array", 0, 1, &
    ibelm_ymin_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_ymax/array", 0, 1, &
    ibelm_ymax_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2D_BOTTOM_IC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_bottom/array", 0, 1, &
    ibelm_bottom_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NSPEC2D_TOP_IC
  start(1) = local_dim * myrank; count(1) = local_dim

  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_top/array", 0, 1, &
    ibelm_top_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)

  ! closes adios file
  call close_file_adios_read()

  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then

    file_name = trim(LOCAL_PATH) // "boundary_disc.bp"

    ! opens adios file
    call open_file_adios_read(file_name)

    ! number of elements
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "NSPEC2D_MOHO", 0, 1, &
       njunk1, adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "NSPEC2D_400", 0, 1, &
       njunk2, adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "NSPEC2D_670", 0, 1, &
       njunk3, adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! checks dimensions
    if (njunk1 /= NSPEC2D_MOHO .and. njunk2 /= NSPEC2D_400 .and. &
        njunk3 /= NSPEC2D_670) &
        call exit_mpi(myrank, 'Error reading boundary_disc.bp file')

    ! boundary elements

    ! moho
    local_dim = NSPEC2D_MOHO
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_moho_top/array", 0, 1, &
      ibelm_moho_bot, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_moho_bot/array", 0, 1, &
      ibelm_moho_top, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! 400
    local_dim = NSPEC2D_400
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_400_top/array", 0, 1, &
      ibelm_400_bot, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_400_bot/array", 0, 1, &
      ibelm_400_top, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! 670
    local_dim = NSPEC2D_670
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_670_top/array", 0, 1, &
      ibelm_670_bot, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "ibelm_670_bot/array", 0, 1, &
      ibelm_670_top, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! normals

    ! moho
    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_moho/array", 0, 1, &
      normal_moho, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! 400
    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_400/array", 0, 1, &
      normal_400, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! 670
    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "normal_670/array", 0, 1, &
      normal_670, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! Close ADIOS handler to the restart file.
    call adios_selection_delete(sel)

    ! closes adios file
    call close_file_adios_read()

  endif

end subroutine read_mesh_databases_coupling_adios

!===============================================================================

subroutine read_mesh_databases_addressing_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod, only: check_adios_err

  implicit none

  ! local parameters
  integer, dimension(0:NPROCTOT_VAL-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer :: ierr,iproc,iproc_read,iproc_xi,iproc_eta

  ! open file with global slice number addressing
  if (myrank == 0) then
    open(unit=IIN,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read',iostat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error opening addressing.txt')

    do iproc = 0,NPROCTOT_VAL-1
      read(IIN,*) iproc_read,ichunk,iproc_xi,iproc_eta

      if (iproc_read /= iproc) call exit_MPI(myrank,'incorrect slice number read')

      addressing(ichunk,iproc_xi,iproc_eta) = iproc
      ichunk_slice(iproc) = ichunk
      iproc_xi_slice(iproc) = iproc_xi
      iproc_eta_slice(iproc) = iproc_eta
    enddo
    close(IIN)
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_i(addressing,NCHUNKS_VAL*NPROC_XI_VAL*NPROC_ETA_VAL)
  call bcast_all_i(ichunk_slice,NPROCTOT_VAL)
  call bcast_all_i(iproc_xi_slice,NPROCTOT_VAL)
  call bcast_all_i(iproc_eta_slice,NPROCTOT_VAL)


  ! output a topology map of slices - fix 20x by nproc
  if (myrank == 0) then
    if (NCHUNKS_VAL == 6 .and. NPROCTOT_VAL < 1000) then
      write(IMAIN,*) 'Spatial distribution of the slices'
      do iproc_xi = NPROC_XI_VAL-1, 0, -1
        write(IMAIN,'(20x)',advance='no')
        do iproc_eta = NPROC_ETA_VAL -1, 0, -1
          ichunk = CHUNK_AB
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      do iproc_xi = NPROC_XI_VAL-1, 0, -1
        write(IMAIN,'(1x)',advance='no')
        do iproc_eta = NPROC_ETA_VAL -1, 0, -1
          ichunk = CHUNK_BC
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(3x)',advance='no')
        do iproc_eta = NPROC_ETA_VAL -1, 0, -1
          ichunk = CHUNK_AC
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(3x)',advance='no')
        do iproc_eta = NPROC_ETA_VAL -1, 0, -1
          ichunk = CHUNK_BC_ANTIPODE
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      do iproc_xi = NPROC_XI_VAL-1, 0, -1
        write(IMAIN,'(20x)',advance='no')
        do iproc_eta = NPROC_ETA_VAL -1, 0, -1
          ichunk = CHUNK_AB_ANTIPODE
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      do iproc_xi = NPROC_XI_VAL-1, 0, -1
        write(IMAIN,'(20x)',advance='no')
        do iproc_eta = NPROC_ETA_VAL -1, 0, -1
          ichunk = CHUNK_AC_ANTIPODE
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
    endif
  endif

  ! determine chunk number and local slice coordinates using addressing
  ! (needed for Stacey conditions)
  ichunk = ichunk_slice(myrank)

end subroutine read_mesh_databases_addressing_adios


!===============================================================================
!> \brief Read crust mantle MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_CM_adios()

  use specfem_par
  use specfem_par_crustmantle

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use manager_adios

  implicit none

  ! local parameters
  integer :: ierr
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  integer :: offset_my_neighbors, offset_nibool_interfaces, &
             offset_ibool_interfaces, offset_phase_ispec_inner, &
             offset_num_elem_colors

  character(len=128)      :: region_name

  write(region_name,"('reg',i1, '/')") IREGION_CRUST_MANTLE

  file_name = trim(LOCAL_PATH) // "/solver_data_mpi.bp"

  ! opens adios file
  call open_file_adios_read(file_name)

  ! MPI interfaces
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_interfaces", 0, 1, &
     num_interfaces_crust_mantle, adios_err)
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec_inner", &
    0, 1, nspec_inner_crust_mantle, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec_outer", &
    0, 1, nspec_outer_crust_mantle, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_phase_ispec", &
    0, 1, num_phase_ispec_crust_mantle, adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  !--------------------------------------.
  ! Get offsets to avoid buffer overflow |
  !--------------------------------------'
  call adios_selection_writeblock(sel, myrank)
  if (num_interfaces_crust_mantle > 0) then
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "my_neighbors/offset", &
                             0, 1,  offset_my_neighbors, adios_err)
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "nibool_interfaces/offset", &
                             0, 1, offset_nibool_interfaces, adios_err)
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "ibool_interfaces/offset", &
                             0, 1, offset_ibool_interfaces, adios_err)
  endif
  if (num_phase_ispec_crust_mantle > 0) then
    call adios_schedule_read(file_handle_adios, sel, &
                              trim(region_name) // "phase_ispec_inner/offset", &
                              0, 1, offset_phase_ispec_inner, adios_err)
  endif
  if (USE_MESH_COLORING_GPU) then
    call adios_schedule_read(file_handle_adios, sel, &
                              trim(region_name) // "num_elem_colors/offset", &
                              0, 1,  offset_num_elem_colors, adios_err)
  endif
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  allocate(my_neighbors_crust_mantle(num_interfaces_crust_mantle), &
          nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
          stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank, &
      'Error allocating array my_neighbors_crust_mantle etc.')

  if (num_interfaces_crust_mantle > 0) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "max_nibool_interfaces", 0, 1, &
       max_nibool_interfaces_cm, adios_err)
    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm, &
        num_interfaces_crust_mantle), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank, &
        'Error allocating array ibool_interfaces_crust_mantle')

    start(1) = offset_my_neighbors
    count(1) = num_interfaces_crust_mantle
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "my_neighbors/array", 0, 1, &
      my_neighbors_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nibool_interfaces/array", &
      0, 1, nibool_interfaces_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    start(1) = offset_ibool_interfaces
    count(1) = max_nibool_interfaces_cm * num_interfaces_crust_mantle
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "ibool_interfaces/array", 0, 1, &
      ibool_interfaces_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! dummy array
    max_nibool_interfaces_cm = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank, &
        'Error allocating array dummy ibool_interfaces_crust_mantle')
  endif

  ! inner / outer elements

  if (num_phase_ispec_crust_mantle < 0 ) &
      call exit_mpi(myrank,'Error num_phase_ispec_crust_mantle is < zero')

  allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2), &
          stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank, &
      'Error allocating array phase_ispec_inner_crust_mantle')

  if (num_phase_ispec_crust_mantle > 0) then
    start(1) = offset_phase_ispec_inner
    count(1) = num_phase_ispec_crust_mantle * 2
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "phase_ispec_inner/array", 0, 1, &
      phase_ispec_inner_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_colors_outer", &
      0, 1, num_colors_outer_crust_mantle, adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_colors_inner", &
      0, 1, num_colors_inner_crust_mantle, adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! colors

    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle +&
        num_colors_inner_crust_mantle), stat=ierr)
    if (ierr /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')

    start(1) = offset_num_elem_colors
    count(1)= num_colors_outer_crust_mantle + num_colors_inner_crust_mantle
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "num_elem_colors/array", 0, 1, &
      num_elem_colors_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + &
        num_colors_inner_crust_mantle), stat=ierr)
    if (ierr /= 0 ) &
      call exit_mpi(myrank, &
          'Error allocating num_elem_colors_crust_mantle array')
  endif
  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)

  ! closes adios file
  call close_file_adios_read()

end subroutine read_mesh_databases_MPI_CM_adios

!===============================================================================
!> \brief Read outer core MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_OC_adios()

  use specfem_par
  use specfem_par_outercore

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use manager_adios

  implicit none

  ! local parameters
  integer :: ierr
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  integer :: offset_my_neighbors, offset_nibool_interfaces, &
             offset_ibool_interfaces, offset_phase_ispec_inner, &
             offset_num_elem_colors

  character(len=128)      :: region_name

  write(region_name,"('reg',i1, '/')") IREGION_OUTER_CORE

  file_name = trim(LOCAL_PATH) // "/solver_data_mpi.bp"

  ! opens adios file
  call open_file_adios_read(file_name)

  ! MPI interfaces
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_interfaces", &
    0, 1, num_interfaces_outer_core, adios_err)
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! inner / outer elements
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec_inner", &
  0, 1, nspec_inner_outer_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec_outer", &
    0, 1, nspec_outer_outer_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_phase_ispec", &
    0, 1, num_phase_ispec_outer_core, adios_err)
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  allocate(my_neighbors_outer_core(num_interfaces_outer_core), &
          nibool_interfaces_outer_core(num_interfaces_outer_core), &
          stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank, &
      'Error allocating array my_neighbors_outer_coreetc.')

  !--------------------------------------.
  ! Get offsets to avoid buffer overflow |
  !--------------------------------------'
  call adios_selection_writeblock(sel, myrank)
  if (num_interfaces_outer_core > 0) then
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "my_neighbors/offset", &
                             0, 1,  offset_my_neighbors, adios_err)
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "nibool_interfaces/offset", &
                             0, 1, offset_nibool_interfaces, adios_err)
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "ibool_interfaces/offset", &
                             0, 1, offset_ibool_interfaces, adios_err)
  endif
  if (num_phase_ispec_outer_core > 0) then
    call adios_schedule_read(file_handle_adios, sel, &
                              trim(region_name) // "phase_ispec_inner/offset", &
                              0, 1, offset_phase_ispec_inner, adios_err)
  endif
  if (USE_MESH_COLORING_GPU) then
    call adios_schedule_read(file_handle_adios, sel, &
                              trim(region_name) // "num_elem_colors/offset", &
                              0, 1,  offset_num_elem_colors, adios_err)
  endif
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  if (num_interfaces_outer_core > 0) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "max_nibool_interfaces", &
      0, 1, max_nibool_interfaces_oc, adios_err)
    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc, &
        num_interfaces_outer_core), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank, &
        'Error allocating array ibool_interfaces_outer_core')

    start(1) = offset_my_neighbors
    count(1) = num_interfaces_outer_core
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "my_neighbors/array", 0, 1, &
      my_neighbors_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nibool_interfaces/array", &
      0, 1, nibool_interfaces_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    start(1) = offset_ibool_interfaces
    count(1) = max_nibool_interfaces_oc * num_interfaces_outer_core
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "ibool_interfaces/array", 0, 1, &
      ibool_interfaces_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! dummy array
    max_nibool_interfaces_oc = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank, &
        'Error allocating array dummy ibool_interfaces_outer_core')
  endif

  if (num_phase_ispec_outer_core < 0 ) &
      call exit_mpi(myrank,'Error num_phase_ispec_outer_core is < zero')

  allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2), &
          stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank, &
      'Error allocating array phase_ispec_inner_outer_core')

  if (num_phase_ispec_outer_core > 0) then
    start(1) = offset_phase_ispec_inner
    count(1) = num_phase_ispec_outer_core * 2
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "phase_ispec_inner/array", 0, 1, &
      phase_ispec_inner_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_colors_outer", &
      0, 1, num_colors_outer_outer_core, adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_colors_inner", &
      0, 1, num_colors_inner_outer_core, adios_err)
    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    ! colors

    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core+&
        num_colors_inner_outer_core), stat=ierr)
    if (ierr /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')

    start(1) = offset_num_elem_colors
    count(1)= num_colors_outer_outer_core + num_colors_inner_outer_core
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "num_elem_colors/array", 0, 1, &
      num_elem_colors_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core+ &
        num_colors_inner_outer_core), stat=ierr)
    if (ierr /= 0 ) &
      call exit_mpi(myrank, &
          'Error allocating num_elem_colors_outer_core array')
  endif
  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)

  ! closes adios file
  call close_file_adios_read()

end subroutine read_mesh_databases_MPI_OC_adios


!===============================================================================
!> \brief Read outer core MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_IC_adios()

  use specfem_par
  use specfem_par_innercore

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use manager_adios

  implicit none

  ! local parameters
  integer :: ierr
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  integer :: offset_my_neighbors, offset_nibool_interfaces, &
             offset_ibool_interfaces, offset_phase_ispec_inner, &
             offset_num_elem_colors

  character(len=128)      :: region_name

  write(region_name,"('reg',i1, '/')") IREGION_INNER_CORE

  file_name = trim(LOCAL_PATH) // "/solver_data_mpi.bp"

  ! opens adios file
  call open_file_adios_read(file_name)

  ! MPI interfaces
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_interfaces", &
    0, 1, num_interfaces_inner_core, adios_err)
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! inner / outer elements
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec_inner", &
    0, 1, nspec_inner_inner_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nspec_outer", &
    0, 1, nspec_outer_inner_core, adios_err)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_phase_ispec", &
    0, 1, num_phase_ispec_inner_core, adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  !--------------------------------------.
  ! Get offsets to avoid buffer overflow |
  !--------------------------------------'
  call adios_selection_writeblock(sel, myrank)
  if (num_interfaces_inner_core > 0) then
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "my_neighbors/offset", &
                             0, 1,  offset_my_neighbors, adios_err)
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "nibool_interfaces/offset", &
                             0, 1, offset_nibool_interfaces, adios_err)
    call adios_schedule_read(file_handle_adios, sel, &
                             trim(region_name) // "ibool_interfaces/offset", &
                             0, 1, offset_ibool_interfaces, adios_err)
  endif
  if (num_phase_ispec_inner_core > 0) then
    call adios_schedule_read(file_handle_adios, sel, &
                              trim(region_name) // "phase_ispec_inner/offset", &
                              0, 1, offset_phase_ispec_inner, adios_err)
  endif
  if (USE_MESH_COLORING_GPU) then
    call adios_schedule_read(file_handle_adios, sel, &
                              trim(region_name) // "num_elem_colors/offset", &
                              0, 1,  offset_num_elem_colors, adios_err)
  endif
  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  allocate(my_neighbors_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &
          stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank, &
      'Error allocating array my_neighbors_inner_core etc.')

  if (num_interfaces_inner_core > 0) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "max_nibool_interfaces", &
      0, 1, max_nibool_interfaces_ic, adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic, &
        num_interfaces_inner_core), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank, &
        'Error allocating array ibool_interfaces_inner_core')

    start(1) = offset_my_neighbors
    count(1) = num_interfaces_inner_core
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "my_neighbors/array", 0, 1, &
      my_neighbors_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nibool_interfaces/array", &
      0, 1, nibool_interfaces_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)

    start(1) = offset_ibool_interfaces
    count(1) = max_nibool_interfaces_ic * num_interfaces_inner_core
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "ibool_interfaces/array", 0, 1, &
      ibool_interfaces_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! dummy array
    max_nibool_interfaces_ic = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank, &
        'Error allocating array dummy ibool_interfaces_inner_core')
  endif

  if (num_phase_ispec_inner_core < 0 ) &
      call exit_mpi(myrank,'Error num_phase_ispec_inner_core is < zero')

  allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2), &
          stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank, &
      'Error allocating array phase_ispec_inner_inner_core')

  if (num_phase_ispec_inner_core > 0) then
    start(1) = offset_phase_ispec_inner
    count(1) = num_phase_ispec_inner_core * 2
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // &
      "phase_ispec_inner/array", 0, 1, &
      phase_ispec_inner_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_colors_outer", &
      0, 1, num_colors_outer_inner_core, adios_err)
    call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "num_colors_inner", &
      0, 1, num_colors_inner_inner_core, adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
    ! colors

    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core +&
        num_colors_inner_inner_core), stat=ierr)
    if (ierr /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')

    start(1) = offset_num_elem_colors
    count(1)= num_colors_outer_inner_core + num_colors_inner_inner_core
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(file_handle_adios, sel, &
      "num_elem_colors/array", 0, 1, &
      num_elem_colors_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(file_handle_adios, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + &
        num_colors_inner_inner_core), stat=ierr)
    if (ierr /= 0 ) &
      call exit_mpi(myrank, &
          'Error allocating num_elem_colors_inner_core array')
  endif
  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)

  ! closes adios file
  call close_file_adios_read()

end subroutine read_mesh_databases_MPI_IC_adios


!===============================================================================
!> \brief Read Stacey BC arrays from an ADIOS file.
subroutine read_mesh_databases_stacey_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use manager_adios

  implicit none

  ! local parameters
  integer :: local_dim
  ! processor identification
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  character(len=128)      :: region_name

  file_name = trim(LOCAL_PATH) // "/stacey.bp"

  ! crust and mantle
  write(region_name,"('reg',i1, '/')") IREGION_CRUST_MANTLE

  ! opens adios file
  call open_file_adios_read(file_name)

  ! read arrays for Stacey conditions
  local_dim = 2*NSPEC2DMAX_XMIN_XMAX_CM
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "njmin/array", 0, 1, &
      njmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "njmax/array", 0, 1, &
      njmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nkmin_xi/array", 0, 1, &
      nkmin_xi_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX_CM
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nimin/array", 0, 1, &
      nimin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nimax/array", 0, 1, &
      nimax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nkmin_eta/array", 0, 1, &
      nkmin_eta_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_selection_delete(sel)


  ! outer core

  write(region_name,"('reg',i1, '/')") IREGION_OUTER_CORE

  ! read arrays for Stacey conditions

  local_dim = 2*NSPEC2DMAX_XMIN_XMAX_OC
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "njmin/array", 0, 1, &
      njmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "njmax/array", 0, 1, &
      njmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nkmin_xi/array", 0, 1, &
      nkmin_xi_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX_OC
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nimin/array", 0, 1, &
      nimin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nimax/array", 0, 1, &
      nimax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "nkmin_eta/array", 0, 1, &
      nkmin_eta_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_selection_delete(sel)

  ! closes adios file
  call close_file_adios_read()

end subroutine read_mesh_databases_stacey_adios

