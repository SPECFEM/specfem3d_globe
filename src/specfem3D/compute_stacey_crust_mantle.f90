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

  subroutine compute_stacey_cm_forward()

! Stacey conditions for forward or adjoint wavefields (SIMULATION_TYPE == 1 or 2)

  use constants_solver

  use specfem_par, only: SIMULATION_TYPE,SAVE_STACEY,it

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    veloc_crust_mantle,accel_crust_mantle, &
    ibool_crust_mantle, &
    rho_vp_crust_mantle,rho_vs_crust_mantle, &
    num_abs_boundary_faces_crust_mantle,abs_boundary_ispec_crust_mantle,abs_boundary_npoin_crust_mantle, &
    abs_boundary_ijk_crust_mantle,abs_boundary_jacobian2Dw_crust_mantle,abs_boundary_normal_crust_mantle, &
    absorb_buffer_crust_mantle,reclen_absorb_buffer_crust_mantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: weight
  real(kind=CUSTOM_REAL) :: vn,vx,vy,vz,nx,ny,nz,tx,ty,tz
  integer :: i,j,k,ispec,iglob,iface,igll

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE == 3) return

  ! crust & mantle

  ! checks if anything to do
  if (num_abs_boundary_faces_crust_mantle == 0) return

  if (.not. GPU_MODE) then
    ! on CPU
! openmp solver
!$OMP PARALLEL if (num_abs_boundary_faces_crust_mantle > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(iface,ispec,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,weight)
!$OMP DO
    do iface = 1,num_abs_boundary_faces_crust_mantle

      ispec = abs_boundary_ispec_crust_mantle(iface)

      do igll = 1,abs_boundary_npoin_crust_mantle(iface)
        i = abs_boundary_ijk_crust_mantle(1,igll,iface)
        j = abs_boundary_ijk_crust_mantle(2,igll,iface)
        k = abs_boundary_ijk_crust_mantle(3,igll,iface)

        weight = abs_boundary_jacobian2Dw_crust_mantle(igll,iface)

        iglob = ibool_crust_mantle(i,j,k,ispec)

        vx = veloc_crust_mantle(1,iglob)
        vy = veloc_crust_mantle(2,iglob)
        vz = veloc_crust_mantle(3,iglob)

        nx = abs_boundary_normal_crust_mantle(1,igll,iface)
        ny = abs_boundary_normal_crust_mantle(2,igll,iface)
        nz = abs_boundary_normal_crust_mantle(3,igll,iface)

        vn = vx*nx+vy*ny+vz*nz

        tx = rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
        ty = rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
        tz = rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

!$OMP ATOMIC
        accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) - tx*weight
!$OMP ATOMIC
        accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) - ty*weight
!$OMP ATOMIC
        accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) - tz*weight

        if (SAVE_STACEY) then
          absorb_buffer_crust_mantle(1,igll,iface) = tx*weight
          absorb_buffer_crust_mantle(2,igll,iface) = ty*weight
          absorb_buffer_crust_mantle(3,igll,iface) = tz*weight
        endif

      enddo
    enddo
!$OMP ENDDO
!$OMP END PARALLEL

  else
    ! on GPU
    call compute_stacey_elastic_gpu(Mesh_pointer,absorb_buffer_crust_mantle)
  endif

  ! writes absorbing boundary values
  if (SAVE_STACEY) call write_abs(0,absorb_buffer_crust_mantle, reclen_absorb_buffer_crust_mantle,it)

  end subroutine compute_stacey_cm_forward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_cm_backward()

! Stacey for backward/reconstructed wavefield

  use constants_solver

  use specfem_par, only: SIMULATION_TYPE,NSTEP,it

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    b_accel_crust_mantle, &
    ibool_crust_mantle, &
    num_abs_boundary_faces_crust_mantle,abs_boundary_ispec_crust_mantle,abs_boundary_npoin_crust_mantle, &
    abs_boundary_ijk_crust_mantle, &
    absorb_buffer_crust_mantle,reclen_absorb_buffer_crust_mantle

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,iglob,iface,igll

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE /= 3) return

  ! crust & mantle

  ! checks if anything to do
  if (num_abs_boundary_faces_crust_mantle == 0) return

  ! reads absorbing boundary values
  ! note: backward/reconstructed wavefields are read in after the Newmark time scheme in the first time loop
  !          this leads to a corresponding boundary condition at time index NSTEP - (it-1) = NSTEP - it + 1
  call read_abs(0,absorb_buffer_crust_mantle,reclen_absorb_buffer_crust_mantle,NSTEP-it+1)

  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1,num_abs_boundary_faces_crust_mantle

      ispec = abs_boundary_ispec_crust_mantle(iface)

      do igll = 1,abs_boundary_npoin_crust_mantle(iface)
        i = abs_boundary_ijk_crust_mantle(1,igll,iface)
        j = abs_boundary_ijk_crust_mantle(2,igll,iface)
        k = abs_boundary_ijk_crust_mantle(3,igll,iface)

        iglob = ibool_crust_mantle(i,j,k,ispec)

        ! note:
        ! this leads to an internal cray compiler error (Cray Fortran Version 8.6.5):
        !b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) - absorb_xmin_crust_mantle(:,j,k,ispec2D)
        ! it helps to be more explicit:
        b_accel_crust_mantle(1,iglob) = b_accel_crust_mantle(1,iglob) - absorb_buffer_crust_mantle(1,igll,iface)
        b_accel_crust_mantle(2,iglob) = b_accel_crust_mantle(2,iglob) - absorb_buffer_crust_mantle(2,igll,iface)
        b_accel_crust_mantle(3,iglob) = b_accel_crust_mantle(3,iglob) - absorb_buffer_crust_mantle(3,igll,iface)

      enddo
    enddo

  else
    ! on GPU
    call compute_stacey_elastic_backward_gpu(Mesh_pointer,absorb_buffer_crust_mantle)
  endif

  end subroutine compute_stacey_cm_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_cm_backward_undoatt()

! Stacey conditions for backward/reconstructed wavefields in UNDO_ATTENUATION case

  use constants_solver

  use specfem_par, only: SIMULATION_TYPE,SAVE_FORWARD,UNDO_ATTENUATION

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    b_veloc_crust_mantle,b_accel_crust_mantle, &
    ibool_crust_mantle, &
    rho_vp_crust_mantle,rho_vs_crust_mantle, &
    num_abs_boundary_faces_crust_mantle,abs_boundary_ispec_crust_mantle,abs_boundary_npoin_crust_mantle, &
    abs_boundary_ijk_crust_mantle,abs_boundary_jacobian2Dw_crust_mantle,abs_boundary_normal_crust_mantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: weight
  real(kind=CUSTOM_REAL) :: vn,vx,vy,vz,nx,ny,nz,tx,ty,tz
  integer :: i,j,k,ispec,iglob,iface,igll

  ! note: we need no file i/o to read in contributions from previous forward run.
  !       the stacey contribution here is calculated again for the reconstruction of the forward wavefield
  !       based on the way we run the undo_attenuation time stepping and temporary storage of the forward wavefield to a buffer.

  ! checks
  if (SIMULATION_TYPE /= 3) return
  if (SAVE_FORWARD) return
  if (.not. UNDO_ATTENUATION) stop 'Error invalid UNDO_ATTENUATION flag for compute_stacey_cm_backward_undoatt()'

  ! crust & mantle

  ! checks if anything to do
  if (num_abs_boundary_faces_crust_mantle == 0) return


  if (.not. GPU_MODE) then
    ! on CPU
    do iface = 1,num_abs_boundary_faces_crust_mantle

      ispec = abs_boundary_ispec_crust_mantle(iface)

      do igll = 1,abs_boundary_npoin_crust_mantle(iface)
        i = abs_boundary_ijk_crust_mantle(1,igll,iface)
        j = abs_boundary_ijk_crust_mantle(2,igll,iface)
        k = abs_boundary_ijk_crust_mantle(3,igll,iface)

        weight = abs_boundary_jacobian2Dw_crust_mantle(igll,iface)

        iglob = ibool_crust_mantle(i,j,k,ispec)

        vx = b_veloc_crust_mantle(1,iglob)
        vy = b_veloc_crust_mantle(2,iglob)
        vz = b_veloc_crust_mantle(3,iglob)

        nx = abs_boundary_normal_crust_mantle(1,igll,iface)
        ny = abs_boundary_normal_crust_mantle(2,igll,iface)
        nz = abs_boundary_normal_crust_mantle(3,igll,iface)

        vn = vx*nx+vy*ny+vz*nz

        tx = rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
        ty = rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
        tz = rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

        b_accel_crust_mantle(1,iglob) = b_accel_crust_mantle(1,iglob) - tx*weight
        b_accel_crust_mantle(2,iglob) = b_accel_crust_mantle(2,iglob) - ty*weight
        b_accel_crust_mantle(3,iglob) = b_accel_crust_mantle(3,iglob) - tz*weight

      enddo
    enddo
  else
    ! on GPU
    call compute_stacey_elastic_undoatt_gpu(Mesh_pointer)
  endif

  end subroutine compute_stacey_cm_backward_undoatt

