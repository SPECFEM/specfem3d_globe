!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

  use specfem_par_trinfinite
  use specfem_par_infinite

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer :: njunk1,njunk2,njunk3
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer(kind=8) :: local_dim
  integer(kind=8) :: sel
  integer(kind=8), dimension(1) :: start, count

  character(len=128) :: region_name

  file_name = get_adios_filename(trim(LOCAL_PATH) // "/boundary")

  ! opens adios file
  call init_adios_group(myadios_group,"BoundaryReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  ! crust and mantle
  if (NSPEC_CRUST_MANTLE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_CRUST_MANTLE

    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmin",nspec2D_xmin_crust_mantle)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmax",nspec2D_xmax_crust_mantle)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymin",nspec2D_ymin_crust_mantle)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymax",nspec2D_ymax_crust_mantle)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "NSPEC2D_BOTTOM",njunk1)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "NSPEC2D_TOP",njunk2)

    ! boundary elements
    local_dim = NSPEC2DMAX_XMIN_XMAX_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmin/array", ibelm_xmin_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmax/array", ibelm_xmax_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2DMAX_YMIN_YMAX_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymin/array", ibelm_ymin_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymax/array", ibelm_ymax_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_BOTTOM_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_bottom/array", ibelm_bottom_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_TOP_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_top/array", ibelm_top_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_xmin/array", normal_xmin_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_xmax/array", normal_xmax_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_ymin/array", normal_ymin_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_ymax/array", normal_ymax_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_bottom/array", normal_bottom_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_top/array", normal_top_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_xmin/array", jacobian2D_xmin_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_xmax/array", jacobian2D_xmax_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_ymin/array", jacobian2D_ymin_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_ymax/array", jacobian2D_ymax_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_bottom/array", jacobian2D_bottom_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLX*NGLLY*NSPEC2D_TOP_CM
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_top/array", jacobian2D_top_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! read parameters to couple fluid and solid regions
  !
  ! outer core
  if (NSPEC_OUTER_CORE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_OUTER_CORE

    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmin",nspec2D_xmin_outer_core)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmax",nspec2D_xmax_outer_core)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymin",nspec2D_ymin_outer_core)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymax",nspec2D_ymax_outer_core)

    ! boundary elements
    local_dim = NSPEC2DMAX_XMIN_XMAX_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmin/array", ibelm_xmin_outer_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmax/array", ibelm_xmax_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2DMAX_YMIN_YMAX_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymin/array", ibelm_ymin_outer_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymax/array", ibelm_ymax_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_BOTTOM_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_bottom/array", ibelm_bottom_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_TOP_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_top/array", ibelm_top_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    ! normals
    local_dim = NDIM*NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_xmin/array", normal_xmin_outer_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_xmax/array", normal_xmax_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_ymin/array", normal_ymin_outer_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_ymax/array", normal_ymax_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_BOTTOM_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_bottom/array", normal_bottom_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_TOP_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "normal_top/array", normal_top_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    ! Jacobians
    local_dim = NGLLY*NGLLZ*NSPEC2DMAX_XMIN_XMAX_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_xmin/array",  jacobian2D_xmin_outer_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_xmax/array", jacobian2D_xmax_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLX*NGLLZ*NSPEC2DMAX_YMIN_YMAX_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_ymin/array", jacobian2D_ymin_outer_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_ymax/array", jacobian2D_ymax_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLX*NGLLY*NSPEC2D_BOTTOM_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_bottom/array", jacobian2D_bottom_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NGLLX*NGLLY*NSPEC2D_TOP_OC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "jacobian2D_top/array", jacobian2D_top_outer_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! inner core
  if (NSPEC_INNER_CORE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_INNER_CORE

    ! number of elements
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmin",nspec2D_xmin_inner_core)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmax",nspec2D_xmax_inner_core)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymin",nspec2D_ymin_inner_core)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymax",nspec2D_ymax_inner_core)


    ! boundary elements
    local_dim = NSPEC2DMAX_XMIN_XMAX_IC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmin/array", ibelm_xmin_inner_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmax/array", ibelm_xmax_inner_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2DMAX_YMIN_YMAX_IC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymin/array", ibelm_ymin_inner_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymax/array", ibelm_ymax_inner_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_BOTTOM_IC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_bottom/array", ibelm_bottom_inner_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_TOP_IC
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_top/array", ibelm_top_inner_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! transition-to-infinite
  if (NSPEC_TRINFINITE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_TRINFINITE

    ! number of elements
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmin",nspec2D_xmin_trinfinite)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmax",nspec2D_xmax_trinfinite)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymin",nspec2D_ymin_trinfinite)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymax",nspec2D_ymax_trinfinite)

    ! boundary elements
    local_dim = NSPEC2DMAX_XMIN_XMAX_TRINF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmin/array", ibelm_xmin_trinfinite)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmax/array", ibelm_xmax_trinfinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2DMAX_YMIN_YMAX_TRINF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymin/array", ibelm_ymin_trinfinite)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymax/array", ibelm_ymax_trinfinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_BOTTOM_TRINF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_bottom/array", ibelm_bottom_trinfinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_TOP_TRINF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_top/array", ibelm_top_trinfinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! infinite
  if (NSPEC_INFINITE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_INFINITE

    ! number of elements
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmin",nspec2D_xmin_infinite)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_xmax",nspec2D_xmax_infinite)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymin",nspec2D_ymin_infinite)
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec2D_ymax",nspec2D_ymax_infinite)

    ! boundary elements
    local_dim = NSPEC2DMAX_XMIN_XMAX_INF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmin/array", ibelm_xmin_infinite)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_xmax/array", ibelm_xmax_infinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2DMAX_YMIN_YMAX_INF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymin/array", ibelm_ymin_infinite)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_ymax/array", ibelm_ymax_infinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_BOTTOM_INF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_bottom/array", ibelm_bottom_infinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    local_dim = NSPEC2D_TOP_INF
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibelm_top/array", ibelm_top_infinite)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"BoundaryReader")

  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
    if (NSPEC_CRUST_MANTLE > 0) then
      file_name = get_adios_filename(trim(LOCAL_PATH) // "boundary_disc")

      ! opens adios file
      call init_adios_group(myadios_group,"BoundaryDiscReader")
      call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

      ! number of elements
      call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "NSPEC2D_MOHO",njunk1)
      call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "NSPEC2D_400",njunk2)
      call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "NSPEC2D_670",njunk3)

      ! checks dimensions
      if (njunk1 /= NSPEC2D_MOHO .or. njunk2 /= NSPEC2D_400 .or. njunk3 /= NSPEC2D_670) then
        print *,'Error: invalid NSPEC2D values read in for solver: ',njunk1,njunk2,njunk3,'(boundary_disc.bin)'
        print *,'       should be MOHO/400/670 : ',NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,'(mesh_parameters.bin)'
        call exit_mpi(myrank, 'Error reading adios boundary_disc file')
      endif

      ! boundary elements

      ! moho
      local_dim = NSPEC2D_MOHO
      start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "ibelm_moho_top/array",ibelm_moho_bot)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "ibelm_moho_bot/array",ibelm_moho_top)

      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! 400
      local_dim = NSPEC2D_400
      start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "ibelm_400_top/array",ibelm_400_bot)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "ibelm_400_bot/array",ibelm_400_top)

      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! 670
      local_dim = NSPEC2D_670
      start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "ibelm_670_top/array",ibelm_670_bot)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "ibelm_670_bot/array",ibelm_670_top)

      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! normals

      ! moho
      local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_MOHO
      start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "normal_moho/array",normal_moho)

      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! 400
      local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_400
      start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "normal_400/array",normal_400)

      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! 670
      local_dim = NDIM*NGLLX*NGLLY*NSPEC2D_670
      start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "normal_670/array",normal_670)

      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! closes adios file
      call close_file_adios_read_and_finalize_method(myadios_file)
      call delete_adios_group(myadios_group,"BoundaryDiscReader")
    endif
  endif

end subroutine read_mesh_databases_coupling_adios

!===============================================================================

subroutine read_mesh_databases_addressing_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

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

  ! broadcast the information read on the main to the nodes
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
!> \brief Read MPI arrays from an ADIOS file.
subroutine read_mesh_databases_MPI_adios(iregion_code)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore

  use specfem_par_trinfinite
  use specfem_par_infinite
  use specfem_par_full_gravity

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer, intent(in) :: iregion_code

  ! local parameters
  integer :: ierr,rank,i
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  integer :: num_interfaces,max_nibool_interfaces, &
             num_phase_ispec,num_colors_outer,num_colors_inner
  integer :: nspec_inner,nspec_outer

  integer(kind=8) :: offset_my_neighbors, offset_nibool_interfaces, &
                     offset_ibool_interfaces, offset_phase_ispec_inner, &
                     offset_num_elem_colors

  ! temporary read arrays
  integer, dimension(:),allocatable :: tmp_nibool_interfaces,tmp_my_neighbors
  integer, dimension(:,:),allocatable :: tmp_ibool_interfaces,tmp_phase_ispec_inner
  integer, dimension(:),allocatable :: tmp_num_elem_colors

  character(len=128)      :: region_name
  integer :: nglob_tmp,nspec_tmp

  file_name = get_adios_filename(trim(LOCAL_PATH) // "/solver_data_mpi")

  ! opens adios file
  call init_adios_group(myadios_group,"SolverMPIReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  write(region_name,"('reg',i1, '/')") iregion_code

  ! file read checking
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "myrank",rank)
  if (rank /= myrank) then
    print *,'Error: reading scalar values from adios file ',trim(file_name)
    print *,'region ',trim(region_name),' got invalid rank number ',rank,' instead of local rank ',myrank
    call exit_mpi(myrank,'Error reading adios file solver_data_mpi')
  endif

  ! MPI interfaces
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "num_interfaces",num_interfaces)
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "max_nibool_interfaces",max_nibool_interfaces)

  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "num_phase_ispec",num_phase_ispec)

  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec_inner",nspec_inner)
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec_outer",nspec_outer)

  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "num_colors_outer",num_colors_outer)
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "num_colors_inner",num_colors_inner)

  ! checks
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    nglob_tmp = NGLOB_CRUST_MANTLE
    nspec_tmp = NSPEC_CRUST_MANTLE
  case (IREGION_OUTER_CORE)
    nglob_tmp = NGLOB_OUTER_CORE
    nspec_tmp = NSPEC_OUTER_CORE
  case (IREGION_INNER_CORE)
    nglob_tmp = NGLOB_INNER_CORE
    nspec_tmp = NSPEC_INNER_CORE
  case (IREGION_TRINFINITE)
    nglob_tmp = NGLOB_TRINFINITE
    nspec_tmp = NSPEC_TRINFINITE
  case (IREGION_INFINITE)
    nglob_tmp = NGLOB_INFINITE
    nspec_tmp = NSPEC_INFINITE
  case default
    call exit_mpi(myrank,'Error invalid iregion_code in read_mesh_databases_MPI_adios() routine')
  end select

  if (num_interfaces < 0) then
    print *,'Error: adios rank ',myrank,' num_interfaces: ',num_interfaces,'should be positive'
    call exit_mpi(myrank,'Error invalid value reading num_interfaces')
  endif
  if (max_nibool_interfaces < 0 .or. max_nibool_interfaces > nglob_tmp) then
    print *,'Error: adios rank ',myrank,' max_nibool_interfaces: ',max_nibool_interfaces,'should be between 1 and ',nglob_tmp
    call exit_mpi(myrank,'Error invalid value reading max_nibool_interfaces')
  endif
  if (num_phase_ispec < 0 .or. num_phase_ispec > nspec_tmp) then
    print *,'Error: adios rank ',myrank,' max_nibool_interfaces: ',num_phase_ispec,'should be between 1 and ',nspec_tmp
    call exit_mpi(myrank,'Error invalid value reading num_phase_ispec')
  endif
  if ( (nspec_inner + nspec_outer) /= nspec_tmp) then
    print *,'Error: adios rank ',myrank,' nspec_inner + nspec_outer: ',nspec_inner,nspec_outer,'should be ',nspec_tmp
    call exit_mpi(myrank,'Error invalid value reading nspec_inner & nspec_outer')
  endif
  if (num_colors_outer < 0 .or. num_colors_outer > nspec_tmp) &
    call exit_mpi(myrank,'Error invalid value reading num_colors_outer')
  if (num_colors_inner < 0 .or. num_colors_inner > nspec_tmp) &
    call exit_mpi(myrank,'Error invalid value reading num_colors_inner')

  !--------------------------------------.
  ! Get offsets to avoid buffer overflow |
  !--------------------------------------'

  if (num_interfaces > 0) then
    call read_adios_scalar(myadios_file,myadios_group,myrank, &
                           trim(region_name) // "my_neighbors/offset",offset_my_neighbors)
    call read_adios_scalar(myadios_file,myadios_group,myrank, &
                           trim(region_name) // "nibool_interfaces/offset",offset_nibool_interfaces)
    call read_adios_scalar(myadios_file,myadios_group,myrank, &
                           trim(region_name) // "ibool_interfaces/offset",offset_ibool_interfaces)
  endif
  if (num_phase_ispec > 0) then
    call read_adios_scalar(myadios_file,myadios_group,myrank, &
                           trim(region_name) // "phase_ispec_inner/offset",offset_phase_ispec_inner)
  endif
  if (USE_MESH_COLORING_GPU) then
    call read_adios_scalar(myadios_file,myadios_group,myrank, &
                           trim(region_name) // "num_elem_colors/offset",offset_num_elem_colors)
  endif

  ! temporary arrays
  allocate(tmp_my_neighbors(num_interfaces), &
           tmp_nibool_interfaces(num_interfaces),stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_crust_mantle etc.')

  allocate(tmp_ibool_interfaces(max_nibool_interfaces,num_interfaces), stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')

  if (num_interfaces > 0) then
! note: we set offset values which usually are equal to local_dim * myrank.
!       this is more flexible than setting it directly as local_dim * myrank in case local_dim varies for different processes.
    start(1) = offset_my_neighbors
    count(1) = int(num_interfaces,kind=8)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "my_neighbors/array", tmp_my_neighbors)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "nibool_interfaces/array", tmp_nibool_interfaces)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    tmp_ibool_interfaces(:,:) = 0

    start(1) = offset_ibool_interfaces
    count(1) = int(max_nibool_interfaces,kind=8) * int(num_interfaces,kind=8)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "ibool_interfaces/array", tmp_ibool_interfaces)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! inner / outer elements
  if (num_phase_ispec < 0 ) call exit_mpi(myrank,'Error num_phase_ispec is < zero')

  allocate(tmp_phase_ispec_inner(num_phase_ispec,2),stat=ierr)
  if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner')

  if (num_phase_ispec > 0) then
    tmp_phase_ispec_inner(:,:) = 0

    start(1) = offset_phase_ispec_inner
    count(1) = int(num_phase_ispec,kind=8) * 2
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "phase_ispec_inner/array", tmp_phase_ispec_inner)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! mesh coloring for GPUs
  allocate(tmp_num_elem_colors(num_colors_outer + num_colors_inner), stat=ierr)
  if (ierr /= 0) call exit_mpi(myrank,'Error allocating num_elem_colors array')

  if (USE_MESH_COLORING_GPU) then
    ! colors
    start(1) = offset_num_elem_colors
    count(1) = int(num_colors_outer,kind=8) + int(num_colors_inner,kind=8)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "num_elem_colors/array", tmp_num_elem_colors)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"SolverMPIReader")

  ! sets region MPI parameters
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle
    num_interfaces_crust_mantle = num_interfaces
    max_nibool_interfaces_cm = max_nibool_interfaces
    num_phase_ispec_crust_mantle = num_phase_ispec
    nspec_inner_crust_mantle = nspec_inner
    nspec_outer_crust_mantle = nspec_outer
    num_colors_outer_crust_mantle = num_colors_outer
    num_colors_inner_crust_mantle = num_colors_inner

    ! MPI arrays
    allocate(my_neighbors_crust_mantle(num_interfaces), &
             nibool_interfaces_crust_mantle(num_interfaces),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_crust_mantle etc.')
    if (num_interfaces > 0) then
      my_neighbors_crust_mantle(1:num_interfaces) = tmp_my_neighbors(1:num_interfaces)
      nibool_interfaces_crust_mantle(1:num_interfaces) = tmp_nibool_interfaces(1:num_interfaces)
    endif

    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
    if (num_interfaces > 0) then
      ibool_interfaces_crust_mantle(:,:) = 0
      do i = 1,num_interfaces
        ibool_interfaces_crust_mantle(1:max_nibool_interfaces,i) = tmp_ibool_interfaces(1:max_nibool_interfaces,i)
      enddo
    endif

    allocate(phase_ispec_inner_crust_mantle(num_phase_ispec,2),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_crust_mantle(:,:) = 0
      ! fills actual values
      do i = 1,2
        phase_ispec_inner_crust_mantle(1:num_phase_ispec,i) = tmp_phase_ispec_inner(1:num_phase_ispec,i)
      enddo
    endif

    ! mesh coloring for GPUs
    allocate(num_elem_colors_crust_mantle(num_colors_outer + num_colors_inner), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
    if (USE_MESH_COLORING_GPU) then
      ! colors
      num_elem_colors_crust_mantle(1:(num_colors_outer + num_colors_inner)) = &
              tmp_num_elem_colors(1:(num_colors_outer + num_colors_inner))
    endif

  case (IREGION_OUTER_CORE)
    ! outer core
    num_interfaces_outer_core = num_interfaces
    max_nibool_interfaces_oc = max_nibool_interfaces
    num_phase_ispec_outer_core = num_phase_ispec
    nspec_inner_outer_core = nspec_inner
    nspec_outer_outer_core = nspec_outer
    num_colors_outer_outer_core = num_colors_outer
    num_colors_inner_outer_core = num_colors_inner

    ! MPI arrays
    allocate(my_neighbors_outer_core(num_interfaces), &
             nibool_interfaces_outer_core(num_interfaces),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_outer_core etc.')
    if (num_interfaces > 0) then
      my_neighbors_outer_core(1:num_interfaces) = tmp_my_neighbors(1:num_interfaces)
      nibool_interfaces_outer_core(1:num_interfaces) = tmp_nibool_interfaces(1:num_interfaces)
    endif

    allocate(ibool_interfaces_outer_core(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_outer_core')
    if (num_interfaces > 0) then
      ibool_interfaces_outer_core(:,:) = 0
      do i = 1,num_interfaces
        ibool_interfaces_outer_core(1:max_nibool_interfaces,i) = tmp_ibool_interfaces(1:max_nibool_interfaces,i)
      enddo
    endif

    allocate(phase_ispec_inner_outer_core(num_phase_ispec,2),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_outer_core(:,:) = 0
      ! fills actual values
      do i = 1,2
        phase_ispec_inner_outer_core(1:num_phase_ispec,i) = tmp_phase_ispec_inner(1:num_phase_ispec,i)
      enddo
    endif

    ! mesh coloring for GPUs
    allocate(num_elem_colors_outer_core(num_colors_outer + num_colors_inner), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
    if (USE_MESH_COLORING_GPU) then
      ! colors
      num_elem_colors_outer_core(1:(num_colors_outer + num_colors_inner)) = &
              tmp_num_elem_colors(1:(num_colors_outer + num_colors_inner))
    endif

  case (IREGION_INNER_CORE)
    ! inner core
    num_interfaces_inner_core = num_interfaces
    max_nibool_interfaces_ic = max_nibool_interfaces
    num_phase_ispec_inner_core = num_phase_ispec
    nspec_inner_inner_core = nspec_inner
    nspec_outer_inner_core = nspec_outer
    num_colors_outer_inner_core = num_colors_outer
    num_colors_inner_inner_core = num_colors_inner

    ! MPI arrays
    allocate(my_neighbors_inner_core(num_interfaces), &
             nibool_interfaces_inner_core(num_interfaces),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_inner_core etc.')
    if (num_interfaces > 0) then
      my_neighbors_inner_core(1:num_interfaces) = tmp_my_neighbors(1:num_interfaces)
      nibool_interfaces_inner_core(1:num_interfaces) = tmp_nibool_interfaces(1:num_interfaces)
    endif

    allocate(ibool_interfaces_inner_core(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_inner_core')
    if (num_interfaces > 0) then
      ibool_interfaces_inner_core(:,:) = 0
      do i = 1,num_interfaces
        ibool_interfaces_inner_core(1:max_nibool_interfaces,i) = tmp_ibool_interfaces(1:max_nibool_interfaces,i)
      enddo
    endif

    allocate(phase_ispec_inner_inner_core(num_phase_ispec,2),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_inner_core(:,:) = 0
      ! fills actual values
      do i = 1,2
        phase_ispec_inner_inner_core(1:num_phase_ispec,i) = tmp_phase_ispec_inner(1:num_phase_ispec,i)
      enddo
    endif

    ! mesh coloring for GPUs
    allocate(num_elem_colors_inner_core(num_colors_outer + num_colors_inner), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
    if (USE_MESH_COLORING_GPU) then
      ! colors
      num_elem_colors_inner_core(1:(num_colors_outer + num_colors_inner)) = &
              tmp_num_elem_colors(1:(num_colors_outer + num_colors_inner))
    endif

  case (IREGION_TRINFINITE)
    ! transition-to-infinite
    num_interfaces_trinfinite = num_interfaces
    max_nibool_interfaces_trinfinite = max_nibool_interfaces
    num_phase_ispec_trinfinite = num_phase_ispec
    nspec_inner_trinfinite = nspec_inner
    nspec_outer_trinfinite = nspec_outer
    num_colors_outer_trinfinite = num_colors_outer
    num_colors_inner_trinfinite = num_colors_inner

    ! MPI arrays
    allocate(my_neighbors_trinfinite(num_interfaces), &
             nibool_interfaces_trinfinite(num_interfaces),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_trinfinite etc.')
    if (num_interfaces > 0) then
      my_neighbors_trinfinite(1:num_interfaces) = tmp_my_neighbors(1:num_interfaces)
      nibool_interfaces_trinfinite(1:num_interfaces) = tmp_nibool_interfaces(1:num_interfaces)
    endif

    allocate(ibool_interfaces_trinfinite(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_trinfinite')
    if (num_interfaces > 0) then
      ibool_interfaces_trinfinite(:,:) = 0
      do i = 1,num_interfaces
        ibool_interfaces_trinfinite(1:max_nibool_interfaces,i) = tmp_ibool_interfaces(1:max_nibool_interfaces,i)
      enddo
    endif

    allocate(phase_ispec_inner_trinfinite(num_phase_ispec,2),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_trinfinite')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_trinfinite(:,:) = 0
      ! fills actual values
      do i = 1,2
        phase_ispec_inner_trinfinite(1:num_phase_ispec,i) = tmp_phase_ispec_inner(1:num_phase_ispec,i)
      enddo
    endif

    ! mesh coloring for GPUs
    allocate(num_elem_colors_trinfinite(num_colors_outer + num_colors_inner), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_trinfinite array')
    if (USE_MESH_COLORING_GPU) then
      ! colors
      num_elem_colors_trinfinite(1:(num_colors_outer + num_colors_inner)) = &
              tmp_num_elem_colors(1:(num_colors_outer + num_colors_inner))
    endif

  case (IREGION_INFINITE)
    ! infinite
    num_interfaces_infinite = num_interfaces
    max_nibool_interfaces_infinite = max_nibool_interfaces
    num_phase_ispec_infinite = num_phase_ispec
    nspec_inner_infinite = nspec_inner
    nspec_outer_infinite = nspec_outer
    num_colors_outer_infinite = num_colors_outer
    num_colors_inner_infinite = num_colors_inner

    ! MPI arrays
    allocate(my_neighbors_infinite(num_interfaces), &
             nibool_interfaces_infinite(num_interfaces),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_infinite etc.')
    if (num_interfaces > 0) then
      my_neighbors_infinite(1:num_interfaces) = tmp_my_neighbors(1:num_interfaces)
      nibool_interfaces_infinite(1:num_interfaces) = tmp_nibool_interfaces(1:num_interfaces)
    endif

    allocate(ibool_interfaces_infinite(max_nibool_interfaces,num_interfaces), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_infinite')
    if (num_interfaces > 0) then
      ibool_interfaces_infinite(:,:) = 0
      do i = 1,num_interfaces
        ibool_interfaces_infinite(1:max_nibool_interfaces,i) = tmp_ibool_interfaces(1:max_nibool_interfaces,i)
      enddo
    endif

    allocate(phase_ispec_inner_infinite(num_phase_ispec,2),stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_infinite')
    if (num_phase_ispec > 0) then
      phase_ispec_inner_infinite(:,:) = 0
      ! fills actual values
      do i = 1,2
        phase_ispec_inner_infinite(1:num_phase_ispec,i) = tmp_phase_ispec_inner(1:num_phase_ispec,i)
      enddo
    endif

    ! mesh coloring for GPUs
    allocate(num_elem_colors_infinite(num_colors_outer + num_colors_inner), stat=ierr)
    if (ierr /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_infinite array')
    if (USE_MESH_COLORING_GPU) then
      ! colors
      num_elem_colors_infinite(1:(num_colors_outer + num_colors_inner)) = &
              tmp_num_elem_colors(1:(num_colors_outer + num_colors_inner))
    endif

  case default
    stop 'Invalid region case in read_mesh_databases_MPI_arrays_adios'
  end select

  ! frees temporary arrays
  deallocate(tmp_my_neighbors,tmp_nibool_interfaces)
  deallocate(tmp_ibool_interfaces)

end subroutine read_mesh_databases_MPI_adios


!===============================================================================
!> \brief Read Stacey BC arrays from an ADIOS file.
subroutine read_mesh_databases_stacey_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer :: ier
  ! processor identification
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer(kind=8) :: sel
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: offset_ispec,offset_ijk,offset_jacobian,offset_normal

  character(len=128) :: region_name

  file_name = get_adios_filename(trim(LOCAL_PATH) // "/stacey")

  ! opens adios file
  call init_adios_group(myadios_group,"StaceyReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  ! crust and mantle
  if (NSPEC_CRUST_MANTLE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_CRUST_MANTLE

    ! read arrays for Stacey conditions
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "num_abs_boundary_faces", &
                           num_abs_boundary_faces_crust_mantle)

    ! reads in arrays
    if (num_abs_boundary_faces_crust_mantle > 0) then
      ! allocates absorbing boundary arrays
      allocate(abs_boundary_ispec_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_crust_mantle(3,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_crust_mantle(NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      allocate(abs_boundary_normal_crust_mantle(NDIM,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
      allocate(abs_boundary_npoin_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

      ! gets offsets
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_ispec/offset",offset_ispec)
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_ijk/offset",offset_ijk)
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_jacobian2Dw/offset",offset_jacobian)
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_normal/offset",offset_normal)

      ! we set offset values which usually are equal to local_dim * myrank.
      ! this is more flexible than setting it directly as local_dim * myrank in case local_dim varies for different processes.
      start(1) = offset_ispec
      count(1) = int(num_abs_boundary_faces_crust_mantle,kind=8)
      call set_selection_boundingbox(sel, start, count)
      ! ispec and npoin arrays have same dimensions
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_ispec/array", abs_boundary_ispec_crust_mantle)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_npoin/array", abs_boundary_npoin_crust_mantle)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! ijk array
      start(1) = offset_ijk
      count(1) = int(3 * NGLLSQUARE,kind=8) * int(num_abs_boundary_faces_crust_mantle,kind=8)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_ijk/array", abs_boundary_ijk_crust_mantle)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! jacobian array
      start(1) = offset_jacobian
      count(1) = int(NGLLSQUARE,kind=8) * int(num_abs_boundary_faces_crust_mantle,kind=8)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_jacobian2Dw/array", &
                                     abs_boundary_jacobian2Dw_crust_mantle)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! normal array
      start(1) = offset_normal
      count(1) = int(NDIM * NGLLSQUARE,kind=8) * int(num_abs_boundary_faces_crust_mantle,kind=8)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_normal/array", abs_boundary_normal_crust_mantle)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)
    else
      ! dummy arrays
      allocate(abs_boundary_ispec_crust_mantle(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_crust_mantle(1,1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_crust_mantle(1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      allocate(abs_boundary_normal_crust_mantle(1,1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
      allocate(abs_boundary_npoin_crust_mantle(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      abs_boundary_ispec_crust_mantle(:) = 0; abs_boundary_npoin_crust_mantle(:) = 0
      abs_boundary_ijk_crust_mantle(:,:,:) = 0
      abs_boundary_jacobian2Dw_crust_mantle(:,:) = 0.0; abs_boundary_normal_crust_mantle(:,:,:) = 0.0
    endif
  endif

  ! outer core
  if (NSPEC_OUTER_CORE > 0) then
    write(region_name,"('reg',i1, '/')") IREGION_OUTER_CORE

    ! read arrays for Stacey conditions
    call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "num_abs_boundary_faces", &
                           num_abs_boundary_faces_outer_core)

    ! reads in arrays
    if (num_abs_boundary_faces_outer_core > 0) then
      ! allocates absorbing boundary arrays
      allocate(abs_boundary_ispec_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_outer_core(3,NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_outer_core(NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      allocate(abs_boundary_npoin_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

      ! gets offsets
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_ispec/offset",offset_ispec)
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_ijk/offset",offset_ijk)
      call read_adios_scalar(myadios_file,myadios_group,myrank, &
                             trim(region_name) // "abs_boundary_jacobian2Dw/offset",offset_jacobian)

      ! we set offset values which usually are equal to local_dim * myrank.
      ! this is more flexible than setting it directly as local_dim * myrank in case local_dim varies for different processes.
      start(1) = offset_ispec
      count(1) = int(num_abs_boundary_faces_outer_core,kind=8)
      call set_selection_boundingbox(sel, start, count)
      ! ispec and npoin arrays have same dimensions
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_ispec/array", abs_boundary_ispec_outer_core)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_npoin/array", abs_boundary_npoin_outer_core)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! ijk array
      start(1) = offset_ijk
      count(1) = int(3 * NGLLSQUARE,kind=8) * int(num_abs_boundary_faces_outer_core,kind=8)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_ijk/array", abs_boundary_ijk_outer_core)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)

      ! jacobian array
      start(1) = offset_jacobian
      count(1) = int(NGLLSQUARE,kind=8) * int(num_abs_boundary_faces_outer_core,kind=8)
      call set_selection_boundingbox(sel, start, count)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "abs_boundary_jacobian2Dw/array", &
                                     abs_boundary_jacobian2Dw_outer_core)
      call read_adios_perform(myadios_file)
      call delete_adios_selection(sel)
    else
      ! dummy arrays
      allocate(abs_boundary_ispec_outer_core(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
      allocate(abs_boundary_ijk_outer_core(1,1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
      allocate(abs_boundary_jacobian2Dw_outer_core(1,1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
      allocate(abs_boundary_npoin_outer_core(1),stat=ier)
      if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      abs_boundary_ispec_outer_core(:) = 0; abs_boundary_npoin_outer_core(:) = 0
      abs_boundary_ijk_outer_core(:,:,:) = 0
      abs_boundary_jacobian2Dw_outer_core(:,:) = 0.0
    endif
  endif

  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"StaceyReader")

end subroutine read_mesh_databases_stacey_adios

