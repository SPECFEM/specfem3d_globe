!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

module combine_vol_data_adios_mod
  use adios_read_mod
  use mpi
  implicit none
contains

!=============================================================================
!> Print help message.
subroutine print_usage_adios()
  print *, 'Usage: '
  print *, '   xcombine_data slice_list varname var_file mesh_file ' // &
           'output_dir high/low-resolution region'
  print *
  print *, '* possible varnames are '
  print *, '   rho, rho, kappastore, mustore, alpha_kernel, etc'
  print *
  print *, '   that are stored in the local directory as ' // &
           'real(kind=CUSTOM_REAL) varname(NGLLX,NGLLY,NGLLZ,NSPEC)  '
  print *, '   in var_file.bp'
  print *
  print *, '* mesh_files are used to link variable to the topology'
  print *, '* output_dir indicates where var_name.vtk will be written'
  print *, '* give 0 for low resolution and 1 for high resolution'
  print *

  stop ' Reenter command line options'
end subroutine print_usage_adios

!=============================================================================
!> Interpret command line arguments
subroutine read_args_adios(arg, MAX_NUM_NODES, node_list, num_node,   &
                           var_name, value_file_name, mesh_file_name, &
                           outdir, ires, irs, ire)
  implicit none
  ! Arguments
  character(len=*), intent(in) :: arg(:)
  integer, intent(in) :: MAX_NUM_NODES
  integer, intent(out) :: node_list(:)
  integer, intent(out) :: num_node, ires, irs, ire
  character(len=*), intent(out) :: var_name, value_file_name, mesh_file_name, &
                                   outdir
  ! Variables
  character(len=256) :: sline
  integer :: ios, njunk, iregion

  if ((command_argument_count() == 6) &
     .or. (command_argument_count() == 7)) then
    num_node = 0
    open(unit = 20, file = trim(arg(1)), status = 'unknown',iostat = ios)
    if (ios /= 0) then
      print *,'Error opening slice file ',trim(arg(1))
      stop
    endif
    do while ( 1 == 1)
      read(20,'(a)',iostat=ios) sline
      if (ios /= 0) exit
      read(sline,*,iostat=ios) njunk
      if (ios /= 0) exit
      num_node = num_node + 1
      if( num_node > MAX_NUM_NODES ) &
          stop 'error number of slices exceeds MAX_NUM_NODES...'
      node_list(num_node) = njunk
    enddo
    close(20)
    var_name = arg(2)
    value_file_name = arg(3)
    mesh_file_name = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  else
    call print_usage_adios()
  endif

  iregion = 0
  if (command_argument_count() == 7) then
    read(arg(7),*) iregion
  endif
  if (iregion > 3 .or. iregion < 0) stop 'Iregion = 0,1,2,3'
  if (iregion == 0) then
    irs = 1
    ire = 3
  else
    irs = iregion
    ire = irs
  endif

end subroutine read_args_adios


!=============================================================================
!> Open ADIOS value and mesh files, read mode
subroutine init_adios(value_file_name, mesh_file_name, &
                      value_handle, mesh_handle)
  implicit none
  ! Parameters
  character(len=*), intent(in) :: value_file_name, mesh_file_name
  integer(kind=8), intent(out) :: value_handle, mesh_handle
  ! Variables
  integer :: ier

  call adios_read_init_method(ADIOS_READ_METHOD_BP, MPI_COMM_WORLD, &
                              "verbose=1", ier)
  call adios_read_open_file(mesh_handle, trim(mesh_file_name), 0, &
                            MPI_COMM_WORLD, ier)
  call adios_read_open_file(value_handle, trim(value_file_name), 0, &
                            MPI_COMM_WORLD, ier)
end subroutine init_adios


!=============================================================================
!> Open ADIOS value and mesh files, read mode
subroutine clean_adios(value_handle, mesh_handle)
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: value_handle, mesh_handle
  ! Variables
  integer :: ier

  call adios_read_close(mesh_handle,ier)
  call adios_read_close(value_handle,ier)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, ier)
end subroutine clean_adios


!=============================================================================
subroutine read_scalars_adios_mesh(mesh_handle, iproc, ir, nglob, nspec)
  implicit none
  ! Parameters
  integer(kind=8), intent(in) :: mesh_handle
  integer, intent(in) :: iproc, ir
  integer, intent(out) :: nglob, nspec
  ! Variables
  integer(kind=8) :: sel
  character(len=256) :: reg_name
  integer :: ier

  write(reg_name, '(a,i1)') trim("reg"), ir

  call adios_selection_writeblock(sel, iproc)
  call adios_schedule_read(mesh_handle, sel, trim(reg_name) // "/nglob", &
                           0, 1, nglob, ier)
  call adios_schedule_read(mesh_handle, sel, trim(reg_name) // "/nspec", &
                           0, 1, nspec, ier)
  call adios_perform_reads(mesh_handle, ier)
end subroutine read_scalars_adios_mesh


!=============================================================================
subroutine read_coordinates_adios_mesh(mesh_handle, iproc, ir, nglob, nspec, &
                                       xstore, ystore, zstore, ibool)
  use constants
  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: mesh_handle
  integer, intent(in) :: iproc, ir, nglob, nspec
  real(kind=CUSTOM_REAL),dimension(:), intent(inout) :: xstore, ystore, zstore
  integer, dimension(:,:,:,:), intent(inout) :: ibool
  ! Variables
  character(len=256) :: reg_name
  integer(kind=8), dimension(1) :: start, count_ad
  integer(kind=8) :: sel_coord, sel_ibool, sel_scalar
  integer :: offset_coord, offset_ibool, ier

  write(reg_name, '(a,i1, a)') trim("reg"), ir, "/"

  call adios_selection_writeblock(sel_scalar, iproc)
  call adios_schedule_read(mesh_handle, sel_scalar,          &
                           trim(reg_name) // "ibool/offset", &
                           0, 1, offset_ibool, ier)
  call adios_schedule_read(mesh_handle, sel_scalar,             &
                           trim(reg_name) // "x_global/offset", &
                           0, 1, offset_coord, ier)
  call adios_perform_reads(mesh_handle, ier)

  start(1) = offset_ibool
  count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox (sel_ibool , 1, start, count_ad)
  call adios_schedule_read(mesh_handle, sel_ibool, &
                           trim(reg_name) // "ibool/array", 0, 1, &
                           ibool, ier)

  start(1) = offset_coord
  count_ad(1) = nglob
  call adios_selection_boundingbox (sel_coord , 1, start, count_ad)
  call adios_schedule_read(mesh_handle, sel_coord, &
                           trim(reg_name) // "x_global/array", 0, 1, &
                           xstore, ier)
  call adios_schedule_read(mesh_handle, sel_coord, &
                           trim(reg_name) // "y_global/array", 0, 1, &
                           ystore, ier)
  call adios_schedule_read(mesh_handle, sel_coord, &
                           trim(reg_name) // "z_global/array", 0, 1, &
                           zstore, ier)
  call adios_perform_reads(mesh_handle, ier)
end subroutine read_coordinates_adios_mesh


!=============================================================================
subroutine read_values_adios(value_handle, var_name, iproc, ir, &
                             nspec, data)
  use constants
  implicit none

  ! Parameters
  integer(kind=8), intent(in) :: value_handle
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: iproc, ir, nspec
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), intent(inout) :: data
  ! Variables
  character(len=256) :: reg_name
  integer(kind=8), dimension(1) :: start, count_ad
  integer(kind=8) :: sel
  integer :: offset, ier

  write(reg_name, '(a,i1, a)') trim("reg"), ir, "/"

  call adios_selection_writeblock(sel, iproc)
  call adios_schedule_read(value_handle, sel,                             &
                           trim(reg_name) // trim(var_name) // "/offset", &
                           0, 1, offset, ier)
  call adios_perform_reads(value_handle, ier)

  start(1) = offset
  count_ad(1) = NGLLX * NGLLY * NGLLZ * nspec
  call adios_selection_boundingbox (sel , 1, start, count_ad)
  call adios_schedule_read(value_handle, sel, &
                           trim(reg_name) // trim(var_name) // "/array", 0, 1,&
                           data, ier)
  call adios_perform_reads(value_handle, ier)
end subroutine read_values_adios

end module combine_vol_data_adios_mod
