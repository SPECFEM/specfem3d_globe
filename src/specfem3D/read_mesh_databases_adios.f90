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


!-------------------------------------------------------------------------------------------------
!

subroutine read_mesh_databases_coupling_adios()
  use adios_read_mod

! to couple mantle with outer core

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  include 'mpif.h'

  ! local parameters
  integer :: njunk1,njunk2,njunk3
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

  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  ! crust and mantle
  ! create name of database
  call create_name_database_adios(prname, IREGION_CRUST_MANTLE, LOCAL_PATH)
  file_name= trim(prname) // "boundary.bp" 

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec2D_xmin", 0, 1, &
     nspec2D_xmin_crust_mantle, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_xmax", 0, 1, &
     nspec2D_xmax_crust_mantle, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_ymin", 0, 1, &
     nspec2D_ymin_crust_mantle, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_ymax", 0, 1, &
     nspec2D_ymax_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_xmin", &
      !nspec2D_xmin_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_xmax", &
      !nspec2D_xmax_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_ymin", &
      !nspec2D_ymin_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_ymax", &
      !nspec2D_ymax_crust_mantle, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2DMAX_XMIN_XMAX_CM
  local_dim = size (ibelm_xmin_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_xmin/array", 0, 1, &
    ibelm_xmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ibelm_xmax/array", 0, 1, &
    ibelm_xmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2DMAX_YMIN_YMAX_CM
  local_dim = size (ibelm_ymin_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_ymin/array", 0, 1, &
    ibelm_ymin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ibelm_ymax/array", 0, 1, &
    ibelm_ymax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2D_BOTTOM_CM
  local_dim = size (ibelm_bottom_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_bottom/array", 0, 1, &
    ibelm_bottom_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2D_TOP_CM
  local_dim = size (ibelm_top_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_top/array", 0, 1, &
    ibelm_top_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_CM
  local_dim = size (normal_xmin_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_xmin/array", 0, 1, &
    normal_xmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "normal_xmax/array", 0, 1, &
    normal_xmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_CM
  local_dim = size (normal_ymin_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_ymin/array", 0, 1, &
    normal_ymin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "normal_ymax/array", 0, 1, &
    normal_ymax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM_CM
  local_dim = size (ibelm_bottom_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_bottom/array", 0, 1, &
    normal_bottom_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP_CM
  local_dim = size (ibelm_top_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_top/array", 0, 1, &
    normal_top_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_CM
  local_dim = size (jacobian2D_xmin_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_xmin/array", 0, 1, &
    jacobian2D_xmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_xmax/array", 0, 1, &
    jacobian2D_xmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_CM
  local_dim = size (jacobian2D_ymin_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_ymin/array", 0, 1, &
    jacobian2D_ymin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_ymax/array", 0, 1, &
    jacobian2D_ymax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM_CM
  local_dim = size (jacobian2D_bottom_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_bottom/array", 0, 1, &
    jacobian2D_bottom_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLX*NGLLY*NSPEC2D_TOP_CM
  local_dim = size (jacobian2D_top_crust_mantle) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_top/array", 0, 1, &
    jacobian2D_top_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
! boundary parameters

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  ! read parameters to couple fluid and solid regions
  !
  ! outer core

  ! create name of database
  call create_name_database_adios(prname, IREGION_OUTER_CORE, LOCAL_PATH)
  file_name= trim(prname) // "boundary.bp" 

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec2D_xmin", 0, 1, &
     nspec2D_xmin_outer_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_xmax", 0, 1, &
     nspec2D_xmax_outer_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_ymin", 0, 1, &
     nspec2D_ymin_outer_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_ymax", 0, 1, &
     nspec2D_ymax_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_xmin", &
      !nspec2D_xmin_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_xmax", &
      !nspec2D_xmax_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_ymin", &
      !nspec2D_ymin_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_ymax", &
      !nspec2D_ymax_outer_core, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  nspec2D_zmin_outer_core = NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

  !local_dim = NSPEC2DMAX_XMIN_XMAX_OC
  local_dim = size (ibelm_xmin_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_xmin/array", 0, 1, &
    ibelm_xmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ibelm_xmax/array", 0, 1, &
    ibelm_xmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2DMAX_YMIN_YMAX_OC
  local_dim = size (ibelm_ymin_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_ymin/array", 0, 1, &
    ibelm_ymin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ibelm_ymax/array", 0, 1, &
    ibelm_ymax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2D_BOTTOM_OC
  local_dim = size (ibelm_bottom_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_bottom/array", 0, 1, &
    ibelm_bottom_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2D_TOP_OC
  local_dim = size (ibelm_top_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_top/array", 0, 1, &
    ibelm_top_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_OC
  local_dim = size (normal_xmin_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_xmin/array", 0, 1, &
    normal_xmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "normal_xmax/array", 0, 1, &
    normal_xmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_OC
  local_dim = size (normal_ymin_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_ymin/array", 0, 1, &
    normal_ymin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "normal_ymax/array", 0, 1, &
    normal_ymax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM_OC
  local_dim = size (normal_bottom_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_bottom/array", 0, 1, &
    normal_bottom_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP_OC
  local_dim = size (normal_top_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "normal_top/array", 0, 1, &
    normal_top_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_OC
  local_dim = size (jacobian2D_xmin_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_xmin/array", 0, 1, &
    jacobian2D_xmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_xmax/array", 0, 1, &
    jacobian2D_xmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_OC
  local_dim = size (jacobian2D_ymin_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_ymin/array", 0, 1, &
    jacobian2D_ymin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_ymax/array", 0, 1, &
    jacobian2D_ymax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)


  !local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM_OC
  local_dim = size (jacobian2D_bottom_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_bottom/array", 0, 1, &
    jacobian2D_bottom_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NGLLX*NGLLY*NSPEC2D_TOP_OC
  local_dim = size (jacobian2D_top_outer_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "jacobian2D_top/array", 0, 1, &
    jacobian2D_top_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  ! boundary parameters

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)


  ! inner core

  ! create name of database
  call create_name_database_adios(prname, IREGION_INNER_CORE, LOCAL_PATH)
  file_name= trim(prname) // "boundary.bp" 

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec2D_xmin", 0, 1, &
     nspec2D_xmin_inner_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_xmax", 0, 1, &
     nspec2D_xmax_inner_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_ymin", 0, 1, &
     nspec2D_ymin_inner_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec2D_ymax", 0, 1, &
     nspec2D_ymax_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_xmin", &
      !nspec2D_xmin_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_xmax", &
      !nspec2D_xmax_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_ymin", &
      !nspec2D_ymin_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec2D_ymax", &
      !nspec2D_ymax_inner_core, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2DMAX_XMIN_XMAX_IC
  local_dim = size (ibelm_xmin_inner_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_xmin/array", 0, 1, &
    ibelm_xmin_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ibelm_xmax/array", 0, 1, &
    ibelm_xmax_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2DMAX_YMIN_YMAX_IC
  local_dim = size (ibelm_ymin_inner_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_ymin/array", 0, 1, &
    ibelm_ymin_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ibelm_ymax/array", 0, 1, &
    ibelm_ymax_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2D_BOTTOM_IC
  local_dim = size (ibelm_bottom_inner_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_bottom/array", 0, 1, &
    ibelm_bottom_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  !local_dim = NSPEC2D_TOP_IC
  local_dim = size (ibelm_top_inner_core) 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibelm_top/array", 0, 1, &
    ibelm_top_inner_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
    file_name = LOCAL_PATH // "boundary_disc.bp" 
    call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
        "verbose=1", adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
    call check_adios_err(myrank,adios_err)

    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(adios_handle, sel, "NSPEC2D_MOHO", 0, 1, &
       njunk1, adios_err)
    call adios_schedule_read(adios_handle, sel, "NSPEC2D_400", 0, 1, &
       njunk2, adios_err)
    call adios_schedule_read(adios_handle, sel, "NSPEC2D_670", 0, 1, &
       njunk3, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
    if (njunk1 /= NSPEC2D_MOHO .and. njunk2 /= NSPEC2D_400 .and. &
        njunk3 /= NSPEC2D_670) &
        call exit_mpi(myrank, 'Error reading boundary_disc.bp file')

    local_dim = NSPEC2D_MOHO
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "ibelm_moho_top/array", 0, 1, &
      ibelm_moho_bot, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "ibelm_moho_bot/array", 0, 1, &
      ibelm_moho_top, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = NSPEC2D_400
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "ibelm_400_top/array", 0, 1, &
      ibelm_400_bot, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "ibelm_400_bot/array", 0, 1, &
      ibelm_400_top, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = NSPEC2D_670
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "ibelm_670_top/array", 0, 1, &
      ibelm_670_bot, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "ibelm_670_bot/array", 0, 1, &
      ibelm_670_top, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "normal_moho/array", 0, 1, &
      normal_moho, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "normal_400/array", 0, 1, &
      normal_400, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "normal_670/array", 0, 1, &
      normal_670, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    ! Close ADIOS handler to the restart file.
    call adios_selection_delete(sel)
    call adios_read_close(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
    call check_adios_err(myrank,adios_err)

    k_top = 1
    k_bot = NGLLZ

    ! initialization
    moho_kl = 0.; d400_kl = 0.; d670_kl = 0.; cmb_kl = 0.; icb_kl = 0.
  endif

  end subroutine read_mesh_databases_coupling_adios

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_addressing_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  include 'mpif.h'

  ! local parameters
  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing
  integer, dimension(0:NPROCTOT_VAL-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer :: ierr,iproc,iproc_read,iproc_xi,iproc_eta

  ! open file with global slice number addressing
  if(myrank == 0) then
    open(unit=IIN,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read',iostat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank,'error opening addressing.txt')

    do iproc = 0,NPROCTOT_VAL-1
      read(IIN,*) iproc_read,ichunk,iproc_xi,iproc_eta

      if(iproc_read /= iproc) call exit_MPI(myrank,'incorrect slice number read')

      addressing(ichunk,iproc_xi,iproc_eta) = iproc
      ichunk_slice(iproc) = ichunk
      iproc_xi_slice(iproc) = iproc_xi
      iproc_eta_slice(iproc) = iproc_eta
    enddo
    close(IIN)
  endif

  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(addressing,NCHUNKS_VAL*NPROC_XI_VAL*NPROC_ETA_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(ichunk_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iproc_xi_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(iproc_eta_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

  ! output a topology map of slices - fix 20x by nproc
  if (myrank == 0 ) then
    if( NCHUNKS_VAL == 6 .and. NPROCTOT_VAL < 1000 ) then
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
  ! (needed for stacey conditions)
  ichunk = ichunk_slice(myrank)

  end subroutine read_mesh_databases_addressing_adios


!-------------------------------------------------------------------------------
!> \brief Read crust mantle MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_CM_adios()
  ! External imports
  use mpi
  use adios_read_mod
  ! Internal imports
  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
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

  ! create the name for the database of the current slide and region
  call create_name_database_adios(prname, IREGION_CRUST_MANTLE, LOCAL_PATH)

  file_name= trim(prname) // "solver_data_mpi.bp" 
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  ! MPI interfaces
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "num_interfaces", 0, 1, &
     num_interfaces_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "num_interfaces", &
      !num_interfaces_crust_mantle, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  allocate(my_neighbours_crust_mantle(num_interfaces_crust_mantle), &
          nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
          stat=ierr)
  if( ierr /= 0 ) call exit_mpi(myrank, &
      'error allocating array my_neighbours_crust_mantle etc.')

  if( num_interfaces_crust_mantle > 0 ) then
    call adios_selection_writeblock(sel, myrank)
    !call adios_get_scalar(adios_handle, "max_nibool_interfaces", &
      !max_nibool_interfaces_cm, adios_err)
    call adios_schedule_read(adios_handle, sel, "max_nibool_interfaces", 0, 1, &
       max_nibool_interfaces_cm, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm, &
        num_interfaces_crust_mantle), stat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank, &
        'error allocating array ibool_interfaces_crust_mantle')

    local_dim = num_interfaces_crust_mantle
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "my_neighbours/array", 0, 1, &
      my_neighbours_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "nibool_interfaces/array", &
      0, 1, nibool_interfaces_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = max_nibool_interfaces_cm * num_interfaces_crust_mantle
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "ibool_interfaces/array", 0, 1, &
      ibool_interfaces_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! dummy array
    max_nibool_interfaces_cm = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank, &
        'error allocating array dummy ibool_interfaces_crust_mantle')
  endif

  ! inner / outer elements
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec_inner", &
    0, 1, nspec_inner_crust_mantle, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec_outer", &
    0, 1, nspec_outer_crust_mantle, adios_err)
  call adios_schedule_read(adios_handle, sel, "num_phase_ispec", &
    0, 1, num_phase_ispec_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "nspec_inner", &
      !nspec_inner_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "nspec_outer", &
      !nspec_outer_crust_mantle, adios_err)
  !call adios_get_scalar(adios_handle, "num_phase_ispec", &
      !num_phase_ispec_crust_mantle, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  if( num_phase_ispec_crust_mantle < 0 ) &
      call exit_mpi(myrank,'error num_phase_ispec_crust_mantle is < zero')

  allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),&
          stat=ierr)
  if( ierr /= 0 ) call exit_mpi(myrank, &
      'error allocating array phase_ispec_inner_crust_mantle')

  if(num_phase_ispec_crust_mantle > 0 ) then
    local_dim = num_phase_ispec_crust_mantle * 2
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "phase_ispec_inner/array", 0, 1, &
      phase_ispec_inner_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(adios_handle, sel, "num_colors_outer", &
      0, 1, num_colors_outer_crust_mantle, adios_err)
    call adios_schedule_read(adios_handle, sel, "num_colors_inner", &
      0, 1, num_colors_inner_crust_mantle, adios_err)
    !call adios_get_scalar(adios_handle, "num_colors_outer", &
        !num_colors_outer_crust_mantle, adios_err)
    !call adios_get_scalar(adios_handle, "num_colors_inner", &
        !num_colors_inner_crust_mantle, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
    ! colors

    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle +&
        num_colors_inner_crust_mantle), stat=ierr)
    if( ierr /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_crust_mantle array')

    local_dim = num_colors_outer_crust_mantle + num_colors_inner_crust_mantle 
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "num_elem_colors/array", 0, 1, &
      num_elem_colors_crust_mantle, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + &
        num_colors_inner_crust_mantle), stat=ierr)
    if( ierr /= 0 ) &
      call exit_mpi(myrank, &
          'error allocating num_elem_colors_crust_mantle array')
  endif
  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call MPI_Barrier(comm, ierr)

end subroutine read_mesh_databases_MPI_CM_adios

!-------------------------------------------------------------------------------
!> \brief Read outer core MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_OC_adios()
  use mpi
  use adios_read_mod
  use specfem_par
  use specfem_par_outercore
  implicit none

  ! local parameters
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

  ! create the name for the database of the current slide and region
  call create_name_database_adios(prname, IREGION_OUTER_CORE, LOCAL_PATH)

  file_name= trim(prname) // "solver_data_mpi.bp" 
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  ! MPI interfaces
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "num_interfaces", &
    0, 1, num_interfaces_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "num_interfaces", &
      !num_interfaces_outer_core, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  allocate(my_neighbours_outer_core(num_interfaces_outer_core), &
          nibool_interfaces_outer_core(num_interfaces_outer_core), &
          stat=ierr)
  if( ierr /= 0 ) call exit_mpi(myrank, &
      'error allocating array my_neighbours_outer_coreetc.')

  if( num_interfaces_outer_core> 0 ) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(adios_handle, sel, "max_nibool_interfaces", &
      0, 1, max_nibool_interfaces_oc, adios_err)
    !call adios_get_scalar(adios_handle, "max_nibool_interfaces", &
      !max_nibool_interfaces_oc, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc, &
        num_interfaces_outer_core), stat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank, &
        'error allocating array ibool_interfaces_outer_core')

    local_dim = num_interfaces_outer_core
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "my_neighbours/array", 0, 1, &
      my_neighbours_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "nibool_interfaces/array", &
      0, 1, nibool_interfaces_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = max_nibool_interfaces_oc * num_interfaces_outer_core
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "ibool_interfaces/array", 0, 1, &
      ibool_interfaces_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! dummy array
    max_nibool_interfaces_oc = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank, &
        'error allocating array dummy ibool_interfaces_outer_core')
  endif

  ! inner / outer elements
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec_inner", &
  0, 1, nspec_inner_outer_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec_outer", &
    0, 1, nspec_outer_outer_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "num_phase_ispec", &
    0, 1, num_phase_ispec_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec_inner", &
      !nspec_inner_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec_outer", &
      !nspec_outer_outer_core, adios_err)
  !call adios_get_scalar(adios_handle, "num_phase_ispec", &
      !num_phase_ispec_outer_core, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  if( num_phase_ispec_outer_core< 0 ) &
      call exit_mpi(myrank,'error num_phase_ispec_outer_core is < zero')

  allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),&
          stat=ierr)
  if( ierr /= 0 ) call exit_mpi(myrank, &
      'error allocating array phase_ispec_inner_outer_core')

  if(num_phase_ispec_outer_core> 0 ) then
    local_dim = num_phase_ispec_outer_core * 2
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "phase_ispec_inner/array", 0, 1, &
      phase_ispec_inner_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(adios_handle, sel, "num_colors_outer", &
      0, 1, num_colors_outer_outer_core, adios_err)
    call adios_schedule_read(adios_handle, sel, "num_colors_inner", &
      0, 1, num_colors_inner_outer_core, adios_err)
    !call adios_get_scalar(adios_handle, "num_colors_outer", &
        !num_colors_outer_outer_core, adios_err)
    !call adios_get_scalar(adios_handle, "num_colors_inner", &
        !num_colors_inner_outer_core, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
    ! colors

    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core+&
        num_colors_inner_outer_core), stat=ierr)
    if( ierr /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_outer_core array')

    local_dim = num_colors_outer_outer_core+ num_colors_inner_outer_core
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "num_elem_colors/array", 0, 1, &
      num_elem_colors_outer_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core+ &
        num_colors_inner_outer_core), stat=ierr)
    if( ierr /= 0 ) &
      call exit_mpi(myrank, &
          'error allocating num_elem_colors_outer_core array')
  endif
  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call MPI_Barrier(comm, ierr)

end subroutine read_mesh_databases_MPI_OC_adios


!-------------------------------------------------------------------------------
!> \brief Read outer core MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_IC_adios()
  use mpi
  use adios_read_mod

  use specfem_par
  use specfem_par_innercore
  implicit none

  ! local parameters
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

  ! create the name for the database of the current slide and region
  call create_name_database_adios(prname, IREGION_INNER_CORE, LOCAL_PATH)

  file_name= trim(prname) // "solver_data_mpi.bp" 
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  ! MPI interfaces
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "num_interfaces", &
    0, 1, num_interfaces_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "num_interfaces", &
      !num_interfaces_inner_core, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  allocate(my_neighbours_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &
          stat=ierr)
  if( ierr /= 0 ) call exit_mpi(myrank, &
      'error allocating array my_neighbours_inner_core etc.')

  if( num_interfaces_inner_core > 0 ) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(adios_handle, sel, "max_nibool_interfaces", &
      0, 1, max_nibool_interfaces_ic, adios_err)
    !call adios_get_scalar(adios_handle, "max_nibool_interfaces", &
      !max_nibool_interfaces_ic, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic, &
        num_interfaces_inner_core), stat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank, &
        'error allocating array ibool_interfaces_inner_core')

    local_dim = num_interfaces_inner_core
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "my_neighbours/array", 0, 1, &
      my_neighbours_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "nibool_interfaces/array", &
      0, 1, nibool_interfaces_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)

    local_dim = max_nibool_interfaces_ic * num_interfaces_inner_core
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "ibool_interfaces/array", 0, 1, &
      ibool_interfaces_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! dummy array
    max_nibool_interfaces_ic = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ierr)
    if( ierr /= 0 ) call exit_mpi(myrank, &
        'error allocating array dummy ibool_interfaces_inner_core')
  endif

  ! inner / outer elements
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec_inner", &
    0, 1, nspec_inner_inner_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "nspec_outer", &
    0, 1, nspec_outer_inner_core, adios_err)
  call adios_schedule_read(adios_handle, sel, "num_phase_ispec", &
    0, 1, num_phase_ispec_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec_inner", &
      !nspec_inner_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "nspec_outer", &
      !nspec_outer_inner_core, adios_err)
  !call adios_get_scalar(adios_handle, "num_phase_ispec", &
      !num_phase_ispec_inner_core, adios_err)
  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  if( num_phase_ispec_inner_core < 0 ) &
      call exit_mpi(myrank,'error num_phase_ispec_inner_core is < zero')

  allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),&
          stat=ierr)
  if( ierr /= 0 ) call exit_mpi(myrank, &
      'error allocating array phase_ispec_inner_inner_core')

  if(num_phase_ispec_inner_core > 0 ) then
    local_dim = num_phase_ispec_inner_core * 2
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "phase_ispec_inner/array", 0, 1, &
      phase_ispec_inner_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    call adios_selection_writeblock(sel, myrank)
    call adios_schedule_read(adios_handle, sel, "num_colors_outer", &
      0, 1, num_colors_outer_inner_core, adios_err)
    call adios_schedule_read(adios_handle, sel, "num_colors_inner", &
      0, 1, num_colors_inner_inner_core, adios_err)
    !call adios_get_scalar(adios_handle, "num_colors_outer", &
        !num_colors_outer_inner_core, adios_err)
    !call adios_get_scalar(adios_handle, "num_colors_inner", &
        !num_colors_inner_inner_core, adios_err)
    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
    ! colors

    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core +&
        num_colors_inner_inner_core), stat=ierr)
    if( ierr /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_inner_core array')

    local_dim = num_colors_outer_inner_core + num_colors_inner_inner_core 
    start(1) = local_dim*myrank; count(1) = local_dim
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, &
      "num_elem_colors/array", 0, 1, &
      num_elem_colors_inner_core, adios_err)
    call check_adios_err(myrank,adios_err)

    call adios_perform_reads(adios_handle, adios_err)
    call check_adios_err(myrank,adios_err)
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + &
        num_colors_inner_inner_core), stat=ierr)
    if( ierr /= 0 ) &
      call exit_mpi(myrank, &
          'error allocating num_elem_colors_inner_core array')
  endif
  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call MPI_Barrier(comm, ierr)

  end subroutine read_mesh_databases_MPI_IC_adios


!
!-------------------------------------------------------------------------------------------------
!

subroutine read_mesh_databases_stacey_adios()

  use mpi
  use adios_read_mod

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ierr, comm, lnspec, lnglob, local_dim
  ! processor identification
  character(len=150) :: reg_name, file_name
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid, sel
  integer(kind=8)         :: adios_groupsize, adios_totalsize
  integer :: vars_count, attrs_count, current_step, last_step, vsteps
  character(len=128), dimension(:), allocatable :: adios_names 
  integer(kind=8), dimension(1) :: start, count

  ! crust and mantle

  ! create name of database
  call create_name_database_adios(reg_name, IREGION_CRUST_MANTLE, LOCAL_PATH)

  file_name= trim(reg_name) // "stacey.bp" 
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)
  ! read arrays for Stacey conditions

  local_dim = 2*NSPEC2DMAX_XMIN_XMAX_CM 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "njmin/array", 0, 1, &
      njmin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "njmax/array", 0, 1, &
      njmax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nkmin_xi/array", 0, 1, &
      nkmin_xi_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX_CM 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nimin/array", 0, 1, &
      nimin_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nimax/array", 0, 1, &
      nimax_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nkmin_eta/array", 0, 1, &
      nkmin_eta_crust_mantle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  ! outer core

  ! create name of database
  call create_name_database_adios(reg_name, IREGION_OUTER_CORE, LOCAL_PATH)

  file_name= trim(reg_name) // "stacey.bp" 
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)
  ! read arrays for Stacey conditions

  local_dim = 2*NSPEC2DMAX_XMIN_XMAX_OC 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "njmin/array", 0, 1, &
      njmin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "njmax/array", 0, 1, &
      njmax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nkmin_xi/array", 0, 1, &
      nkmin_xi_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX_OC 
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nimin/array", 0, 1, &
      nimin_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nimax/array", 0, 1, &
      nimax_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "nkmin_eta/array", 0, 1, &
      nkmin_eta_outer_core, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)


end subroutine read_mesh_databases_stacey_adios

