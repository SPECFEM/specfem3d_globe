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

subroutine get_absorb_adios(myrank, iregion, nimin, nimax, njmin, njmax, &
  nkmin_xi, nkmin_eta, NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  use mpi
  use adios_write_mod
  use meshfem3D_par, only: LOCAL_PATH

  ! Stacey, define flags for absorbing boundaries
  implicit none

  include "constants.h"

  integer :: myrank
  integer :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX

  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nimin,nimax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX) :: njmin,njmax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX) :: nkmin_xi
  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nkmin_eta

  character(len=150) :: reg_name, outputname, group_name
  integer :: sizeprocs, comm, local_dim, ierr, iregion
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize

  call create_name_database_adios(reg_name,iregion,LOCAL_PATH)

  outputname = trim(reg_name) // "stacey.bp" 

  ! save these temporary arrays for the solver for Stacey conditions
  write(group_name,"('SPECFEM3D_GLOBE_STACEY_reg',i1)") iregion
  call world_size(sizeprocs) ! TODO keep it in parameters
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)
  group_size_inc = 0
  call adios_declare_group(adios_group, group_name, &
      "", 0, adios_err)
  call adios_select_method(adios_group, "MPI", "", "", adios_err)

  local_dim = 2*NSPEC2DMAX_XMIN_XMAX 
  call define_adios_global_integer_1d_array(adios_group, "njmin", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "njmax", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "nkmin_xi", &
      local_dim, group_size_inc)
  local_dim = 2*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_integer_1d_array(adios_group, "nimin", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "nimax", &
      local_dim, group_size_inc)
  call define_adios_global_integer_1d_array(adios_group, "nkmin_eta", &
      local_dim, group_size_inc)

  ! Open an ADIOS handler to the restart file.
  call adios_open (adios_handle, group_name, &
      outputname, "w", comm, adios_err);
  call adios_group_size (adios_handle, group_size_inc, &
                         adios_totalsize, adios_err)

  local_dim = 2*NSPEC2DMAX_XMIN_XMAX 
  call adios_set_path (adios_handle, "njmin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", njmin, adios_err)

  call adios_set_path (adios_handle, "njmax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", njmax, adios_err)

  call adios_set_path (adios_handle, "nkmin_xi", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", nkmin_xi, adios_err)

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX 
  call adios_set_path (adios_handle, "nimin", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", nimin, adios_err)

  call adios_set_path (adios_handle, "nimax", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", nimax, adios_err)

  call adios_set_path (adios_handle, "nkmin_eta", adios_err)
  call write_1D_global_array_adios_dims(adios_handle, myrank, &
      local_dim, sizeprocs)
  call adios_write(adios_handle, "array", nkmin_eta, adios_err)

  call adios_set_path (adios_handle, "", adios_err)
  call adios_close(adios_handle, adios_err)

end subroutine get_absorb_adios

