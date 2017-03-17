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


!===============================================================================
!> \brief Read adios attenuation arrays created by the mesher
!         (regX_attenuation.bp)
subroutine read_attenuation_adios(iregion_code, &
                                  factor_common, scale_factor, tau_s, vnspec, T_c_source)

  use constants_solver
  use specfem_par, only: ATTENUATION_VAL,LOCAL_PATH

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err
  use manager_adios

  implicit none

  integer :: vnspec
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: scale_factor
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,N_SLS,vnspec) :: factor_common
  double precision, dimension(N_SLS)                 :: tau_s

  integer :: iregion_code

  ! local parameters
  double precision :: T_c_source
  character(len=MAX_STRING_LEN) :: file_name
  integer :: local_dim
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: sel
  integer(kind=8), dimension(1) :: start, count

  character(len=128)      :: region_name, region_name_scalar

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  ! checks if attenuation is on and anything to do
  if (.not. ATTENUATION_VAL) return

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read
  file_name= trim(LOCAL_PATH) // "/attenuation.bp"

  ! opens adios file
  call open_file_adios_read(file_name)

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "t_c_source", 0, 1, &
     T_c_source, adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = size (tau_s)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "tau_s/array", 0, 1, &
    tau_s, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = size (factor_common)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "tau_e_store/array", 0, 1, &
    factor_common, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = size (scale_factor)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(file_handle_adios, sel, trim(region_name) // "Qmu_store/array", 0, 1, &
    scale_factor, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)

  ! closes adios file
  call close_file_adios_read()

end subroutine read_attenuation_adios
