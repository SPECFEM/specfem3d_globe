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

  subroutine compute_stacey_cm_forward()

! Stacey conditions for forward or adjoint wavefields (SIMULATION_TYPE == 1 or 2)

  use constants_solver

  use specfem_par, only: &
    ichunk,SIMULATION_TYPE,SAVE_STACEY,it, &
    wgllwgll_xz,wgllwgll_yz

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    veloc_crust_mantle,accel_crust_mantle, &
    ibool_crust_mantle, &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
    normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
    normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
    rho_vp_crust_mantle,rho_vs_crust_mantle, &
    ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
    ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
    nimin_crust_mantle,nimax_crust_mantle, &
    njmin_crust_mantle,njmax_crust_mantle, &
    nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
    nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
    reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
    reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
    absorb_xmin_crust_mantle,absorb_xmax_crust_mantle, &
    absorb_ymin_crust_mantle,absorb_ymax_crust_mantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: weight
  real(kind=CUSTOM_REAL) :: vn,vx,vy,vz,nx,ny,nz,tx,ty,tz
  integer :: i,j,k,ispec,iglob,ispec2D
  !integer :: reclen1,reclen2

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE == 3 ) return

  ! crust & mantle
  if (.not. GPU_MODE) then
    ! on CPU
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,i,j,k,iglob,vx,vy,vz,vn,nx,ny,nz,tx,ty,tz,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xz,wgllwgll_yz)

    !   xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
!$OMP DO
      do ispec2D = 1,nspec2D_xmin_crust_mantle

        ispec=ibelm_xmin_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_xi_crust_mantle(1,ispec2D) == 0 .or. njmin_crust_mantle(1,ispec2D) == 0) cycle

        i = 1
        do k = nkmin_xi_crust_mantle(1,ispec2D),NGLLZ
          do j = njmin_crust_mantle(1,ispec2D),njmax_crust_mantle(1,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            vx=veloc_crust_mantle(1,iglob)
            vy=veloc_crust_mantle(2,iglob)
            vz=veloc_crust_mantle(3,iglob)

            nx=normal_xmin_crust_mantle(1,j,k,ispec2D)
            ny=normal_xmin_crust_mantle(2,j,k,ispec2D)
            nz=normal_xmin_crust_mantle(3,j,k,ispec2D)

            vn=vx*nx+vy*ny+vz*nz

            tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
            ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
            tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

            weight=jacobian2D_xmin_crust_mantle(j,k,ispec2D)*wgllwgll_yz(j,k)

!$OMP ATOMIC
            accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
!$OMP ATOMIC
            accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
!$OMP ATOMIC
            accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

            if (SAVE_STACEY) then
              absorb_xmin_crust_mantle(1,j,k,ispec2D) = tx*weight
              absorb_xmin_crust_mantle(2,j,k,ispec2D) = ty*weight
              absorb_xmin_crust_mantle(3,j,k,ispec2D) = tz*weight
            endif
          enddo
        enddo
      enddo
!$OMP enddo NOWAIT

    endif

    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
!$OMP DO
      do ispec2D = 1,nspec2D_xmax_crust_mantle

        ispec=ibelm_xmax_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_xi_crust_mantle(2,ispec2D) == 0 .or. njmin_crust_mantle(2,ispec2D) == 0) cycle

        i=NGLLX
        do k=nkmin_xi_crust_mantle(2,ispec2D),NGLLZ
          do j=njmin_crust_mantle(2,ispec2D),njmax_crust_mantle(2,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            vx=veloc_crust_mantle(1,iglob)
            vy=veloc_crust_mantle(2,iglob)
            vz=veloc_crust_mantle(3,iglob)

            nx=normal_xmax_crust_mantle(1,j,k,ispec2D)
            ny=normal_xmax_crust_mantle(2,j,k,ispec2D)
            nz=normal_xmax_crust_mantle(3,j,k,ispec2D)

            vn=vx*nx+vy*ny+vz*nz

            tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
            ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
            tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

            weight=jacobian2D_xmax_crust_mantle(j,k,ispec2D)*wgllwgll_yz(j,k)

!$OMP ATOMIC
            accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
!$OMP ATOMIC
            accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
!$OMP ATOMIC
            accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

            if (SAVE_STACEY) then
              absorb_xmax_crust_mantle(1,j,k,ispec2D) = tx*weight
              absorb_xmax_crust_mantle(2,j,k,ispec2D) = ty*weight
              absorb_xmax_crust_mantle(3,j,k,ispec2D) = tz*weight
            endif

          enddo
        enddo
      enddo
!$OMP enddo NOWAIT

    endif

    !   ymin
!$OMP DO
    do ispec2D = 1,nspec2D_ymin_crust_mantle

      ispec=ibelm_ymin_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_eta_crust_mantle(1,ispec2D) == 0 .or. nimin_crust_mantle(1,ispec2D) == 0) cycle

      j = 1
      do k = nkmin_eta_crust_mantle(1,ispec2D),NGLLZ
        do i = nimin_crust_mantle(1,ispec2D),nimax_crust_mantle(1,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          vx=veloc_crust_mantle(1,iglob)
          vy=veloc_crust_mantle(2,iglob)
          vz=veloc_crust_mantle(3,iglob)

          nx=normal_ymin_crust_mantle(1,i,k,ispec2D)
          ny=normal_ymin_crust_mantle(2,i,k,ispec2D)
          nz=normal_ymin_crust_mantle(3,i,k,ispec2D)

          vn=vx*nx+vy*ny+vz*nz

          tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
          ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
          tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

          weight=jacobian2D_ymin_crust_mantle(i,k,ispec2D)*wgllwgll_xz(i,k)

!$OMP ATOMIC
          accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
!$OMP ATOMIC
          accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
!$OMP ATOMIC
          accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

          if (SAVE_STACEY) then
            absorb_ymin_crust_mantle(1,i,k,ispec2D) = tx*weight
            absorb_ymin_crust_mantle(2,i,k,ispec2D) = ty*weight
            absorb_ymin_crust_mantle(3,i,k,ispec2D) = tz*weight
          endif

        enddo
      enddo
    enddo
!$OMP enddo NOWAIT

    !   ymax
!$OMP DO
    do ispec2D = 1,nspec2D_ymax_crust_mantle

      ispec=ibelm_ymax_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_eta_crust_mantle(2,ispec2D) == 0 .or. nimin_crust_mantle(2,ispec2D) == 0) cycle

      j=NGLLY
      do k=nkmin_eta_crust_mantle(2,ispec2D),NGLLZ
        do i=nimin_crust_mantle(2,ispec2D),nimax_crust_mantle(2,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          vx=veloc_crust_mantle(1,iglob)
          vy=veloc_crust_mantle(2,iglob)
          vz=veloc_crust_mantle(3,iglob)

          nx=normal_ymax_crust_mantle(1,i,k,ispec2D)
          ny=normal_ymax_crust_mantle(2,i,k,ispec2D)
          nz=normal_ymax_crust_mantle(3,i,k,ispec2D)

          vn=vx*nx+vy*ny+vz*nz

          tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
          ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
          tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

          weight=jacobian2D_ymax_crust_mantle(i,k,ispec2D)*wgllwgll_xz(i,k)

!$OMP ATOMIC
          accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
!$OMP ATOMIC
          accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
!$OMP ATOMIC
          accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

          if (SAVE_STACEY) then
            absorb_ymax_crust_mantle(1,i,k,ispec2D) = tx*weight
            absorb_ymax_crust_mantle(2,i,k,ispec2D) = ty*weight
            absorb_ymax_crust_mantle(3,i,k,ispec2D) = tz*weight
          endif

        enddo
      enddo
    enddo
!$OMP enddo
!$OMP END PARALLEL

  else
    ! on GPU
    !   xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
      if (nspec2D_xmin_crust_mantle > 0 ) call compute_stacey_elastic_gpu(Mesh_pointer,absorb_xmin_crust_mantle,0) ! <= xmin
    endif
    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
      if (nspec2D_xmax_crust_mantle > 0 ) call compute_stacey_elastic_gpu(Mesh_pointer,absorb_xmax_crust_mantle,1) ! <= xmin
    endif
    ! ymin
    if (nspec2D_ymin_crust_mantle > 0 ) call compute_stacey_elastic_gpu(Mesh_pointer,absorb_ymin_crust_mantle,2) ! <= ymin
    !   ymax
    if (nspec2D_ymax_crust_mantle > 0 ) call compute_stacey_elastic_gpu(Mesh_pointer,absorb_ymax_crust_mantle,3) ! <= ymax
  endif

  ! writes absorbing boundary values
  if (SAVE_STACEY) then
    ! xmin
    if (nspec2D_xmin_crust_mantle > 0) call write_abs(0,absorb_xmin_crust_mantle, reclen_xmin_crust_mantle,it)
    ! xmax
    if (nspec2D_xmax_crust_mantle > 0) call write_abs(1,absorb_xmax_crust_mantle,reclen_xmax_crust_mantle,it)
    ! ymin
    if (nspec2D_ymin_crust_mantle > 0) call write_abs(2,absorb_ymin_crust_mantle,reclen_ymin_crust_mantle,it)
    ! ymax
    if (nspec2D_ymax_crust_mantle > 0) call write_abs(3,absorb_ymax_crust_mantle,reclen_ymax_crust_mantle,it)
  endif


  end subroutine compute_stacey_cm_forward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_cm_backward()

! Stacey for backward/reconstructed wavefield

  use constants_solver

  use specfem_par, only: &
    ichunk,SIMULATION_TYPE,NSTEP,it

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    b_accel_crust_mantle, &
    ibool_crust_mantle, &
    ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
    ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
    nimin_crust_mantle,nimax_crust_mantle, &
    njmin_crust_mantle,njmax_crust_mantle, &
    nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
    nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
    reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
    reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
    absorb_xmin_crust_mantle,absorb_xmax_crust_mantle, &
    absorb_ymin_crust_mantle,absorb_ymax_crust_mantle

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,iglob,ispec2D

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! crust & mantle

  !   xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then

    ! reads absorbing boundary values
    if (nspec2D_xmin_crust_mantle > 0) then
      ! note: backward/reconstructed wavefields are read in after the Newmark time scheme in the first time loop
      !          this leads to a corresponding boundary condition at time index NSTEP - (it-1) = NSTEP - it + 1
      call read_abs(0,absorb_xmin_crust_mantle,reclen_xmin_crust_mantle,NSTEP-it+1)
    endif

    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_xmin_crust_mantle

        ispec=ibelm_xmin_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_xi_crust_mantle(1,ispec2D) == 0 .or. njmin_crust_mantle(1,ispec2D) == 0) cycle

        i = 1
        do k = nkmin_xi_crust_mantle(1,ispec2D),NGLLZ
          do j = njmin_crust_mantle(1,ispec2D),njmax_crust_mantle(1,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_xmin_crust_mantle(:,j,k,ispec2D)
          enddo
        enddo
      enddo

    else
      ! on GPU
      if (nspec2D_xmin_crust_mantle > 0 ) call compute_stacey_elastic_backward_gpu(Mesh_pointer, &
                                                                                   absorb_xmin_crust_mantle, &
                                                                                   0) ! <= xmin
    endif

  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC

  !   xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then

    ! reads absorbing boundary values
    if (nspec2D_xmax_crust_mantle > 0) then
      call read_abs(1,absorb_xmax_crust_mantle,reclen_xmax_crust_mantle,NSTEP-it+1)
    endif

    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_xmax_crust_mantle

        ispec=ibelm_xmax_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_xi_crust_mantle(2,ispec2D) == 0 .or. njmin_crust_mantle(2,ispec2D) == 0) cycle

        i=NGLLX
        do k=nkmin_xi_crust_mantle(2,ispec2D),NGLLZ
          do j=njmin_crust_mantle(2,ispec2D),njmax_crust_mantle(2,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_xmax_crust_mantle(:,j,k,ispec2D)
          enddo
        enddo
      enddo

    else
      ! on GPU
      if (nspec2D_xmax_crust_mantle > 0 ) call compute_stacey_elastic_backward_gpu(Mesh_pointer, &
                                                                                   absorb_xmax_crust_mantle, &
                                                                                   1) ! <= xmin
    endif

  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB

  !   ymin

  ! reads absorbing boundary values
  if (nspec2D_ymin_crust_mantle > 0) then
    call read_abs(2,absorb_ymin_crust_mantle, reclen_ymin_crust_mantle,NSTEP-it+1)
  endif

  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_ymin_crust_mantle

      ispec=ibelm_ymin_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_eta_crust_mantle(1,ispec2D) == 0 .or. nimin_crust_mantle(1,ispec2D) == 0) cycle

      j = 1
      do k = nkmin_eta_crust_mantle(1,ispec2D),NGLLZ
        do i = nimin_crust_mantle(1,ispec2D),nimax_crust_mantle(1,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_ymin_crust_mantle(:,i,k,ispec2D)
        enddo
      enddo
    enddo

  else
    ! on GPU
    if (nspec2D_ymin_crust_mantle > 0 ) call compute_stacey_elastic_backward_gpu(Mesh_pointer, &
                                                                                 absorb_ymin_crust_mantle, &
                                                                                 2) ! <= ymin
  endif

  !   ymax

  ! reads absorbing boundary values
  if (nspec2D_ymax_crust_mantle > 0) then
    call read_abs(3,absorb_ymax_crust_mantle,reclen_ymax_crust_mantle,NSTEP-it+1)
  endif

  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_ymax_crust_mantle

      ispec=ibelm_ymax_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_eta_crust_mantle(2,ispec2D) == 0 .or. nimin_crust_mantle(2,ispec2D) == 0) cycle

      j=NGLLY
      do k=nkmin_eta_crust_mantle(2,ispec2D),NGLLZ
        do i=nimin_crust_mantle(2,ispec2D),nimax_crust_mantle(2,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_ymax_crust_mantle(:,i,k,ispec2D)
        enddo
      enddo
    enddo

  else
    ! on GPU
    if (nspec2D_ymax_crust_mantle > 0 ) call compute_stacey_elastic_backward_gpu(Mesh_pointer, &
                                                                                 absorb_ymax_crust_mantle, &
                                                                                 3) ! <= ymax
  endif

  end subroutine compute_stacey_cm_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_cm_backward_undoatt()

! Stacey conditions for backward/reconstructed wavefields in UNDO_ATTENUATION case

  use constants_solver

  use specfem_par, only: &
    ichunk,SIMULATION_TYPE,SAVE_FORWARD,UNDO_ATTENUATION, &
    wgllwgll_xz,wgllwgll_yz

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    b_veloc_crust_mantle,b_accel_crust_mantle, &
    ibool_crust_mantle, &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
    normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
    normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
    rho_vp_crust_mantle,rho_vs_crust_mantle, &
    ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
    ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
    nimin_crust_mantle,nimax_crust_mantle, &
    njmin_crust_mantle,njmax_crust_mantle, &
    nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
    nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: weight
  real(kind=CUSTOM_REAL) :: vn,vx,vy,vz,nx,ny,nz,tx,ty,tz
  integer :: i,j,k,ispec,iglob,ispec2D
  !integer :: reclen1,reclen2

  ! note: we need no file i/o to read in contributions from previous forward run.
  !       the stacey contribution here is calculated again for the reconstruction of the forward wavefield
  !       based on the way we run the undo_attenuation time stepping and temporary storage of the forward wavefield to a buffer.

  ! checks
  if (SIMULATION_TYPE /= 3) return
  if (SAVE_FORWARD) return
  if (.not. UNDO_ATTENUATION) stop 'Error invalid UNDO_ATTENUATION flag for compute_stacey_cm_backward_undoatt()'

  ! crust & mantle

  !   xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then

    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_xmin_crust_mantle

        ispec=ibelm_xmin_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_xi_crust_mantle(1,ispec2D) == 0 .or. njmin_crust_mantle(1,ispec2D) == 0) cycle

        i = 1
        do k = nkmin_xi_crust_mantle(1,ispec2D),NGLLZ
          do j = njmin_crust_mantle(1,ispec2D),njmax_crust_mantle(1,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            vx = b_veloc_crust_mantle(1,iglob)
            vy = b_veloc_crust_mantle(2,iglob)
            vz = b_veloc_crust_mantle(3,iglob)

            nx=normal_xmin_crust_mantle(1,j,k,ispec2D)
            ny=normal_xmin_crust_mantle(2,j,k,ispec2D)
            nz=normal_xmin_crust_mantle(3,j,k,ispec2D)

            vn=vx*nx+vy*ny+vz*nz

            tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
            ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
            tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

            weight=jacobian2D_xmin_crust_mantle(j,k,ispec2D)*wgllwgll_yz(j,k)

            b_accel_crust_mantle(1,iglob) = b_accel_crust_mantle(1,iglob) - tx*weight
            b_accel_crust_mantle(2,iglob) = b_accel_crust_mantle(2,iglob) - ty*weight
            b_accel_crust_mantle(3,iglob) = b_accel_crust_mantle(3,iglob) - tz*weight

          enddo
        enddo
      enddo

    else
      ! on GPU
      if (nspec2D_xmin_crust_mantle > 0 ) call compute_stacey_elastic_undoatt_gpu(Mesh_pointer,0) ! <= xmin
    endif

  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC

  !   xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then

    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_xmax_crust_mantle

        ispec=ibelm_xmax_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_xi_crust_mantle(2,ispec2D) == 0 .or. njmin_crust_mantle(2,ispec2D) == 0) cycle

        i=NGLLX
        do k=nkmin_xi_crust_mantle(2,ispec2D),NGLLZ
          do j=njmin_crust_mantle(2,ispec2D),njmax_crust_mantle(2,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            vx = b_veloc_crust_mantle(1,iglob)
            vy = b_veloc_crust_mantle(2,iglob)
            vz = b_veloc_crust_mantle(3,iglob)

            nx=normal_xmax_crust_mantle(1,j,k,ispec2D)
            ny=normal_xmax_crust_mantle(2,j,k,ispec2D)
            nz=normal_xmax_crust_mantle(3,j,k,ispec2D)

            vn=vx*nx+vy*ny+vz*nz

            tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
            ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
            tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

            weight=jacobian2D_xmax_crust_mantle(j,k,ispec2D)*wgllwgll_yz(j,k)

            b_accel_crust_mantle(1,iglob) = b_accel_crust_mantle(1,iglob) - tx*weight
            b_accel_crust_mantle(2,iglob) = b_accel_crust_mantle(2,iglob) - ty*weight
            b_accel_crust_mantle(3,iglob) = b_accel_crust_mantle(3,iglob) - tz*weight

          enddo
        enddo
      enddo

    else
      ! on GPU
      if (nspec2D_xmax_crust_mantle > 0 ) call compute_stacey_elastic_undoatt_gpu(Mesh_pointer,1) ! <= xmin
    endif

  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB

  !   ymin

  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_ymin_crust_mantle

      ispec=ibelm_ymin_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_eta_crust_mantle(1,ispec2D) == 0 .or. nimin_crust_mantle(1,ispec2D) == 0) cycle

      j = 1
      do k = nkmin_eta_crust_mantle(1,ispec2D),NGLLZ
        do i = nimin_crust_mantle(1,ispec2D),nimax_crust_mantle(1,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          vx = b_veloc_crust_mantle(1,iglob)
          vy = b_veloc_crust_mantle(2,iglob)
          vz = b_veloc_crust_mantle(3,iglob)

          nx=normal_ymin_crust_mantle(1,i,k,ispec2D)
          ny=normal_ymin_crust_mantle(2,i,k,ispec2D)
          nz=normal_ymin_crust_mantle(3,i,k,ispec2D)

          vn=vx*nx+vy*ny+vz*nz

          tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
          ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
          tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

          weight=jacobian2D_ymin_crust_mantle(i,k,ispec2D)*wgllwgll_xz(i,k)

          b_accel_crust_mantle(1,iglob) = b_accel_crust_mantle(1,iglob) - tx*weight
          b_accel_crust_mantle(2,iglob) = b_accel_crust_mantle(2,iglob) - ty*weight
          b_accel_crust_mantle(3,iglob) = b_accel_crust_mantle(3,iglob) - tz*weight

        enddo
      enddo
    enddo

  else
    ! on GPU
    if (nspec2D_ymin_crust_mantle > 0 ) call compute_stacey_elastic_undoatt_gpu(Mesh_pointer,2) ! <= ymin
  endif

  !   ymax

  if (.not. GPU_MODE) then
    ! on CPU
    do ispec2D = 1,nspec2D_ymax_crust_mantle

      ispec=ibelm_ymax_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_eta_crust_mantle(2,ispec2D) == 0 .or. nimin_crust_mantle(2,ispec2D) == 0) cycle

      j=NGLLY
      do k=nkmin_eta_crust_mantle(2,ispec2D),NGLLZ
        do i=nimin_crust_mantle(2,ispec2D),nimax_crust_mantle(2,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          vx = b_veloc_crust_mantle(1,iglob)
          vy = b_veloc_crust_mantle(2,iglob)
          vz = b_veloc_crust_mantle(3,iglob)

          nx=normal_ymax_crust_mantle(1,i,k,ispec2D)
          ny=normal_ymax_crust_mantle(2,i,k,ispec2D)
          nz=normal_ymax_crust_mantle(3,i,k,ispec2D)

          vn=vx*nx+vy*ny+vz*nz

          tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
          ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
          tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

          weight=jacobian2D_ymax_crust_mantle(i,k,ispec2D)*wgllwgll_xz(i,k)

          b_accel_crust_mantle(1,iglob) = b_accel_crust_mantle(1,iglob) - tx*weight
          b_accel_crust_mantle(2,iglob) = b_accel_crust_mantle(2,iglob) - ty*weight
          b_accel_crust_mantle(3,iglob) = b_accel_crust_mantle(3,iglob) - tz*weight

        enddo
      enddo
    enddo

  else
    ! on GPU
    if (nspec2D_ymax_crust_mantle > 0 ) call compute_stacey_elastic_undoatt_gpu(Mesh_pointer,3) ! <= ymax
  endif

  end subroutine compute_stacey_cm_backward_undoatt

