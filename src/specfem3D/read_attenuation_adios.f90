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
!> \brief Read adios attenuation arrays created by the mesher
!         (regX_attenuation.bp)
  subroutine read_attenuation_adios(iregion_code, factor_common, scale_factor, tau_s, vnspec, f_c_source)

  use constants_solver
  use specfem_par, only: ATTENUATION_VAL,LOCAL_PATH

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer,intent(in) :: iregion_code

  integer,intent(in) :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec),intent(inout) :: scale_factor
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec),intent(inout) :: factor_common
  double precision, dimension(N_SLS),intent(inout) :: tau_s
  double precision,intent(inout) :: f_c_source

  ! local parameters
  character(len=MAX_STRING_LEN) :: file_name
  ! ADIOS variables
  integer(kind=8) :: local_dim
  integer(kind=8) :: sel
  integer(kind=8), dimension(1) :: start, count

  character(len=128) :: region_name, region_name_scalar

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  ! checks if attenuation is on and anything to do
  if (.not. ATTENUATION_VAL) return

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read
  file_name = get_adios_filename(trim(LOCAL_PATH) // "/attenuation")

  ! opens adios file
  call init_adios_group(myadios_group,"AttenuationReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "f_c_source",f_c_source)

  local_dim = size (tau_s)
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "tau_s_store/array", tau_s)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  local_dim = size (factor_common)
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "tau_e_store/array", factor_common)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  local_dim = size (scale_factor)
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "Qmu_store/array", scale_factor)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! closes ADIOS handler to the restart file.
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"AttenuationReader")

  end subroutine read_attenuation_adios
