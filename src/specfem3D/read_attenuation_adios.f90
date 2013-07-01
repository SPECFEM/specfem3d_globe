!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
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


!===============================================================================
!> \brief Read adios attenuation arrays created by the mesher
!         (regX_attenuation.bp)
subroutine read_attenuation_adios(myrank, prname, &
   factor_common, scale_factor, tau_s, vx, vy, vz, vnspec, T_c_source)

  use adios_read_mod
  use specfem_par,only: ATTENUATION_VAL

  implicit none

  include 'constants.h'
  include 'mpif.h'

  integer :: myrank

  integer :: vx,vy,vz,vnspec
  double precision, dimension(vx,vy,vz,vnspec)       :: scale_factor
  double precision, dimension(N_SLS,vx,vy,vz,vnspec) :: factor_common
  double precision, dimension(N_SLS)                 :: tau_s

  character(len=150) :: prname

  ! local parameters
  integer :: i,j,k,ispec,ier
  double precision, dimension(N_SLS) :: tau_e, fc
  double precision :: omsb, Q_mu, sf, T_c_source, scale_t
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

  ! checks if attenuation is on and anything to do
  if( .not. ATTENUATION_VAL) return

  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read
  file_name= trim(prname) // "attenuation.bp"

  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "T_c_source", 0, 1, &
     T_c_source, adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = size (tau_s)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "tau_s/array", 0, 1, &
    tau_s, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = size (factor_common)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "tau_e_store/array", 0, 1, &
    factor_common, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = size (scale_factor)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "Qmu_store/array", 0, 1, &
    scale_factor, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  ! Close ADIOS handler to the restart file.
  call adios_selection_delete(sel)
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

end subroutine read_attenuation_adios
