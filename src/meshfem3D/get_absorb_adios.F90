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

! old version: not used anymore, but left for reference...

  use constants

  use meshfem_par, only: LOCAL_PATH

  use adios_helpers_mod
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
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=128)      :: region_name

  integer, save :: num_regions_written = 0

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  !call create_name_database_adios(reg_name,iregion,LOCAL_PATH)

  write(region_name,"('reg',i1, '/')") iregion

  ! Append the actual file name.
  !outputname = trim(reg_name) // "stacey.bp"
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/stacey.old_format")

  ! save these temporary arrays for the solver for Stacey conditions
  write(group_name,"('SPECFEM3D_GLOBE_STACEY_reg',i1)") iregion

  ! set the adios group size to 0 before incremented by calls to helpers functions.
  group_size_inc = 0
  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  local_dim = 2*NSPEC2DMAX_XMIN_XMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(njmin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(njmax))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(nkmin_xi))

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(nimin))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(nimax))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, STRINGIFY_VAR(nkmin_eta))

  !--- Open an ADIOS handler to the restart file. ---------
  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
  endif

  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  local_dim = 2*NSPEC2DMAX_XMIN_XMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(njmin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(njmax))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nkmin_xi))

  local_dim = 2*NSPEC2DMAX_YMIN_YMAX
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nimin))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nimax))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(nkmin_eta))

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  num_regions_written = num_regions_written + 1

  end subroutine get_absorb_adios


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_absorb_stacey_boundary_adios(iregion, num_abs_boundary_faces, &
                                              abs_boundary_ispec,abs_boundary_npoin, &
                                              abs_boundary_ijk,abs_boundary_normal,abs_boundary_jacobian2Dw)


  use constants, only: NDIM,NGLLX,NGLLY,NGLLSQUARE,CUSTOM_REAL,MAX_STRING_LEN,myrank,IREGION_CRUST_MANTLE

  use meshfem_par, only: LOCAL_PATH

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer,intent(in) :: iregion

  ! absorbing boundary arrays
  integer,intent(in) :: num_abs_boundary_faces
  integer, dimension(num_abs_boundary_faces), intent(in) :: abs_boundary_ispec
  integer, dimension(num_abs_boundary_faces), intent(in) :: abs_boundary_npoin
  integer, dimension(3,NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_ijk
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_normal
  real(kind=CUSTOM_REAL), dimension(NGLLSQUARE,num_abs_boundary_faces), intent(in) :: abs_boundary_jacobian2Dw

  ! local parameters
  integer :: num_abs_boundary_faces_max

  ! ADIOS variables
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  character(len=128) :: region_name,region_name_scalar

  integer, save :: num_regions_written = 0

  ! since different slice will have different number of absorbing faces,
  ! for adios we determine the maximum size for the local_dim assignements of the arrays
  call max_all_all_i(num_abs_boundary_faces,num_abs_boundary_faces_max)

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  !call create_name_database_adios(reg_name,iregion,LOCAL_PATH)

  write(region_name,"('reg',i1, '/')") iregion
  write(region_name_scalar,"('reg',i1)") iregion

  ! Append the actual file name.
  !outputname = trim(reg_name) // "stacey.bp"
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/stacey")

  ! save these temporary arrays for the solver for Stacey conditions
  write(group_name,"('SPECFEM3D_GLOBE_STACEY_reg',i1)") iregion

  ! set the adios group size to 0 before incremented by calls to helpers functions.
  call init_adios_group(myadios_group,group_name)

  !--- Define ADIOS variables -----------------------------
  group_size_inc = 0
  call define_adios_scalar (myadios_group, group_size_inc, region_name_scalar, STRINGIFY_VAR(num_abs_boundary_faces))

  local_dim = num_abs_boundary_faces_max
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, &
                                   STRINGIFY_VAR(abs_boundary_ispec))
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, &
                                   STRINGIFY_VAR(abs_boundary_npoin))

  local_dim = 3 * NGLLSQUARE * num_abs_boundary_faces_max
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, &
                                   STRINGIFY_VAR(abs_boundary_ijk))

  local_dim = NGLLSQUARE * num_abs_boundary_faces_max
  call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, &
                                   STRINGIFY_VAR(abs_boundary_jacobian2Dw))

  ! normals only needed for elastic boundary conditions
  if (iregion == IREGION_CRUST_MANTLE) then
    local_dim = NDIM * NGLLSQUARE * num_abs_boundary_faces_max
    call define_adios_global_array1D(myadios_group, group_size_inc, local_dim, region_name, &
                                     STRINGIFY_VAR(abs_boundary_normal))
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  if (num_regions_written == 0) then
    ! opens file for writing
    call open_file_adios_write(myadios_file,myadios_group,outputname,group_name)
  else
    ! opens file for writing in append mode
    call open_file_adios_write_append(myadios_file,myadios_group,outputname,group_name)
  endif

  call set_adios_group_size(myadios_file,group_size_inc)

  !--- Schedule writes for the previously defined ADIOS variables
  call write_adios_scalar(myadios_file,myadios_group,trim(region_name) // STRINGIFY_VAR(num_abs_boundary_faces))

  ! note: using the *_wmax array sizes for local_dim is providing the same local_dim/global_dim/offset values
  !       in the adios file for all rank processes. this mimicks the same chunk sizes for all processes in ADIOS files.
  !       this helps when reading back arrays using offsets based on the local_dim value.
  !
  ! we thus use num_abs_boundary_faces_max here rather than num_abs_boundary_faces
  local_dim = num_abs_boundary_faces_max
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(abs_boundary_ispec))
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(abs_boundary_npoin))

  local_dim = 3 * NGLLSQUARE * num_abs_boundary_faces_max
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(abs_boundary_ijk))

  local_dim = NGLLSQUARE * num_abs_boundary_faces_max
  call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                   local_dim, trim(region_name) // STRINGIFY_VAR(abs_boundary_jacobian2Dw))

  ! normals only needed for elastic boundary conditions
  if (iregion == IREGION_CRUST_MANTLE) then
    local_dim = NDIM * NGLLSQUARE * num_abs_boundary_faces_max
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, &
                                     local_dim, trim(region_name) // STRINGIFY_VAR(abs_boundary_normal))
  endif

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_file)
  ! closes file
  call close_file_adios(myadios_file)

  num_regions_written = num_regions_written + 1

  end subroutine get_absorb_stacey_boundary_adios
