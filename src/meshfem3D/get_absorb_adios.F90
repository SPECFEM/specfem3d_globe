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
!> \file get_absorb_adios.f90
!! \brief function to write Stacey boundary condition to disk with ADIOS.
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"

!===============================================================================
!> \brief Write Stacey boundary conditions to a single file using ADIOS
!!
!! \param iregion The region the absorbing condition is written for. Check
!!                constant.h files to see what these regions are.
!! \param nimin An array to be written
!! \param nimax An array to be written
!! \param njmin An array to be written
!! \param njmax An array to be written
!! \param nkmin_xi An array to be written
!! \param nkmin_eta An array to be written
!! \param NSPEC2DMAX_XMIN_XMAX Integer to compute the size of the arrays
!!                             in argument
!! \param NSPEC2DMAX_YMIN_YMAX Integer to compute the size of the arrays
!!                             in argument
!!
!! \note This routine only call adios to write the file to disk, Note that he
!!       necessary data preparation is done by the get_absorb() routine.
  subroutine get_absorb_adios(iregion, &
                              nimin, nimax, njmin, njmax, nkmin_xi, nkmin_eta, &
                              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

  use constants

  use meshfem3D_par, only: LOCAL_PATH

  use adios_helpers_mod, only: define_adios_global_array1D,write_adios_global_1d_array
  use manager_adios

  ! Stacey, define flags for absorbing boundaries
  implicit none

  integer :: iregion
  integer :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX

  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nimin,nimax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX) :: njmin,njmax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX) :: nkmin_xi
  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nkmin_eta

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  integer(kind=8)         :: adios_group
  character(len=128)      :: region_name

  integer, save :: num_regions_written = 0

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  !call create_name_database_adios(reg_name,iregion,LOCAL_PATH)

  write(region_name,"('reg',i1, '/')") iregion

  ! Append the actual file name.
  !outputname = trim(reg_name) // "stacey.bp"
  outputname = trim(LOCAL_PATH) // "/stacey.bp"

  ! save these temporary arrays for the solver for Stacey conditions
  write(group_name,"('SPECFEM3D_GLOBE_STACEY_reg',i1)") iregion

  ! set the adios group size to 0 before incremented by calls to helpers functions.
  group_size_inc = 0
  call init_adios_group(adios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  local_dim = 2*NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, &
                                   region_name, STRINGIFY_VAR(njmin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, &
                                   region_name, STRINGIFY_VAR(njmax))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, &
                                   region_name, STRINGIFY_VAR(nkmin_xi))

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, &
                                   region_name, STRINGIFY_VAR(nimin))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, &
                                   region_name, STRINGIFY_VAR(nimax))
  call define_adios_global_array1D(adios_group, group_size_inc, local_dim, &
                                   region_name, STRINGIFY_VAR(nkmin_eta))

  !--- Open an ADIOS handler to the restart file. ---------
  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(outputname,group_name)
  endif
  call set_adios_group_size(group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  local_dim = 2*NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(njmin))
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(njmax))
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nkmin_xi))

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nimin))
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nimax))
  call write_adios_global_1d_array(file_handle_adios, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nkmin_eta))

  !--- Reset the path to zero and perform the actual write to disk
  ! closes file
  call close_file_adios()

  num_regions_written = num_regions_written + 1

  end subroutine get_absorb_adios

