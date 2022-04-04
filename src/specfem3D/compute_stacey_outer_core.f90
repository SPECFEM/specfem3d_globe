!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  subroutine compute_stacey_oc_forward()

  use constants_solver

  use specfem_par, only: &
    SIMULATION_TYPE,SAVE_STACEY,it

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_outercore, only: &
    veloc_outer_core,accel_outer_core, &
    ibool_outer_core, &
    vp_outer_core, &
    num_abs_boundary_faces_outer_core,abs_boundary_ispec_outer_core,abs_boundary_npoin_outer_core, &
    abs_boundary_ijk_outer_core,abs_boundary_jacobian2Dw_outer_core, &
    absorb_buffer_outer_core,reclen_absorb_buffer_outer_core

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: sn,weight
  integer :: i,j,k,igll,iface,ispec,iglob

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE == 3) return

  ! outer core

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0 .or. num_abs_boundary_faces_outer_core == 0) return

  if (.not. GPU_MODE) then
    ! on CPU
! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces_outer_core > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,i,j,k,iglob,sn,weight)
!$OMP DO
    do iface = 1,num_abs_boundary_faces_outer_core

      ispec = abs_boundary_ispec_outer_core(iface)

      do igll = 1,abs_boundary_npoin_outer_core(iface)
        i = abs_boundary_ijk_outer_core(1,igll,iface)
        j = abs_boundary_ijk_outer_core(2,igll,iface)
        k = abs_boundary_ijk_outer_core(3,igll,iface)

        weight = abs_boundary_jacobian2Dw_outer_core(igll,iface)

        iglob = ibool_outer_core(i,j,k,ispec)

        sn = veloc_outer_core(iglob) / vp_outer_core(i,j,k,ispec)

        !$OMP ATOMIC
        accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

        if (SAVE_STACEY) then
          absorb_buffer_outer_core(igll,iface) = weight*sn
        endif

      enddo
    enddo
!$OMP ENDDO
!$OMP END PARALLEL

  else
    ! on GPU
    ! adds boundary contributions
    call compute_stacey_acoustic_gpu(Mesh_pointer,absorb_buffer_outer_core)
  endif

  ! writes absorbing boundary values to file
  if (SAVE_STACEY) call write_abs(4,absorb_buffer_outer_core,reclen_absorb_buffer_outer_core,it)

  end subroutine compute_stacey_oc_forward


!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_stacey_oc_backward()

  use constants_solver

  use specfem_par, only: &
    SIMULATION_TYPE,NSTEP,it

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_outercore, only: &
    b_accel_outer_core, &
    ibool_outer_core, &
    num_abs_boundary_faces_outer_core,abs_boundary_ispec_outer_core,abs_boundary_npoin_outer_core, &
    abs_boundary_ijk_outer_core, &
    absorb_buffer_outer_core,reclen_absorb_buffer_outer_core

  implicit none

  ! local parameters
  integer :: i,j,k,igll,iface,ispec,iglob

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! outer core

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0 .or. num_abs_boundary_faces_outer_core == 0) return

  ! reads absorbing boundary values
  ! note: backward/reconstructed wavefields are read in after the Newmark time scheme in the first time loop
  !          this leads to a corresponding boundary condition at time index NSTEP - (it-1) = NSTEP - it + 1
  call read_abs(4,absorb_buffer_outer_core,reclen_absorb_buffer_outer_core,NSTEP-it+1)

  ! adds boundary contributions
  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1,num_abs_boundary_faces_outer_core

      ispec = abs_boundary_ispec_outer_core(iface)

      do igll = 1,abs_boundary_npoin_outer_core(iface)
        i = abs_boundary_ijk_outer_core(1,igll,iface)
        j = abs_boundary_ijk_outer_core(2,igll,iface)
        k = abs_boundary_ijk_outer_core(3,igll,iface)

        iglob = ibool_outer_core(i,j,k,ispec)

        b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_buffer_outer_core(igll,iface)
      enddo
    enddo

  else
    ! on GPU
    call compute_stacey_acoustic_backward_gpu(Mesh_pointer,absorb_buffer_outer_core)
  endif

  end subroutine compute_stacey_oc_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_oc_backward_undoatt()

  use constants_solver

  use specfem_par, only: &
    SIMULATION_TYPE,SAVE_FORWARD,UNDO_ATTENUATION

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_outercore, only: &
    b_veloc_outer_core,b_accel_outer_core, &
    ibool_outer_core, &
    vp_outer_core, &
    num_abs_boundary_faces_outer_core,abs_boundary_ispec_outer_core,abs_boundary_npoin_outer_core, &
    abs_boundary_ijk_outer_core,abs_boundary_jacobian2Dw_outer_core

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: sn,weight
  integer :: i,j,k,igll,iface,ispec,iglob

  ! checks
  if (SIMULATION_TYPE /= 3) return
  if (SAVE_FORWARD) return
  if (.not. UNDO_ATTENUATION) stop 'Error invalid UNDO_ATTENUATION flag for compute_stacey_oc_backward_undoatt()'

  ! outer core

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0 .or. num_abs_boundary_faces_outer_core == 0) return

  ! adds boundary contribution
  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1,num_abs_boundary_faces_outer_core

      ispec = abs_boundary_ispec_outer_core(iface)

      do igll = 1,abs_boundary_npoin_outer_core(iface)
        i = abs_boundary_ijk_outer_core(1,igll,iface)
        j = abs_boundary_ijk_outer_core(2,igll,iface)
        k = abs_boundary_ijk_outer_core(3,igll,iface)

        weight = abs_boundary_jacobian2Dw_outer_core(igll,iface)

        iglob = ibool_outer_core(i,j,k,ispec)

        sn = b_veloc_outer_core(iglob) / vp_outer_core(i,j,k,ispec)

        b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - weight*sn

      enddo
    enddo
  else
    ! on GPU
    call compute_stacey_acoustic_undoatt_gpu(Mesh_pointer,4) ! <= xmin
  endif

  end subroutine compute_stacey_oc_backward_undoatt


