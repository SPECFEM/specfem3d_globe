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

  subroutine compute_stacey_oc_forward()

  use constants_solver

  use specfem_par, only: &
    ichunk,SIMULATION_TYPE,SAVE_STACEY,it, &
    wgllwgll_xz,wgllwgll_yz,wgllwgll_xy

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_outercore, only: &
    veloc_outer_core,accel_outer_core, &
    ibool_outer_core, &
    jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
    jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
    jacobian2D_bottom_outer_core, &
    nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core,nspec2D_zmin_outer_core, &
    vp_outer_core, &
    nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core, &
    njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core, &
    absorb_xmin_outer_core,absorb_xmax_outer_core, &
    absorb_ymin_outer_core,absorb_ymax_outer_core, &
    absorb_zmin_outer_core, &
    reclen_xmin_outer_core,reclen_xmax_outer_core, &
    reclen_ymin_outer_core,reclen_ymax_outer_core, &
    reclen_zmin, &
    ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
    ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
    ibelm_bottom_outer_core
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: sn,weight
  integer :: i,j,k,ispec2D,ispec,iglob

  ! note: we use C functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storage (every 8 step) and/or asynchronous
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! checks
  if (SIMULATION_TYPE == 3 ) return

  ! outer core

  if (.not. GPU_MODE) then
    ! on CPU
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,i,j,k,iglob,sn,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy,wgllwgll_xz,wgllwgll_yz)

    !   xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
      if (nspec2D_xmin_outer_core > 0) then
        ! adds boundary contributions
!$OMP DO
        do ispec2D = 1,nspec2D_xmin_outer_core

          ispec = ibelm_xmin_outer_core(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi_outer_core(1,ispec2D) == 0 .or. njmin_outer_core(1,ispec2D) == 0) cycle

          i = 1
          do k = nkmin_xi_outer_core(1,ispec2D),NGLLZ
            do j = njmin_outer_core(1,ispec2D),njmax_outer_core(1,ispec2D)
              iglob=ibool_outer_core(i,j,k,ispec)

              sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

              weight = jacobian2D_xmin_outer_core(j,k,ispec2D)*wgllwgll_yz(j,k)

!$OMP ATOMIC
              accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

              if (SAVE_STACEY) then
                absorb_xmin_outer_core(j,k,ispec2D) = weight*sn
              endif
            enddo
          enddo
        enddo
!$OMP ENDDO NOWAIT
      endif
    endif

    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
      if (nspec2D_xmax_outer_core > 0) then
        ! adds boundary contributions
!$OMP DO
        do ispec2D = 1,nspec2D_xmax_outer_core

          ispec = ibelm_xmax_outer_core(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi_outer_core(2,ispec2D) == 0 .or. njmin_outer_core(2,ispec2D) == 0) cycle

          i=NGLLX
          do k=nkmin_xi_outer_core(2,ispec2D),NGLLZ
            do j=njmin_outer_core(2,ispec2D),njmax_outer_core(2,ispec2D)
              iglob=ibool_outer_core(i,j,k,ispec)

              sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

              weight = jacobian2D_xmax_outer_core(j,k,ispec2D)*wgllwgll_yz(j,k)

!$OMP ATOMIC
              accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

              if (SAVE_STACEY) then
                absorb_xmax_outer_core(j,k,ispec2D) = weight*sn
              endif

            enddo
          enddo
        enddo
!$OMP ENDDO NOWAIT
      endif
    endif

    !   ymin
    if (nspec2D_ymin_outer_core > 0) then
      ! adds boundary contributions
!$OMP DO
      do ispec2D = 1,nspec2D_ymin_outer_core

        ispec = ibelm_ymin_outer_core(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_eta_outer_core(1,ispec2D) == 0 .or. nimin_outer_core(1,ispec2D) == 0) cycle

        j = 1
        do k = nkmin_eta_outer_core(1,ispec2D),NGLLZ
          do i = nimin_outer_core(1,ispec2D),nimax_outer_core(1,ispec2D)
            iglob=ibool_outer_core(i,j,k,ispec)

            sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

            weight = jacobian2D_ymin_outer_core(i,k,ispec2D)*wgllwgll_xz(i,k)

!$OMP ATOMIC
            accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

            if (SAVE_STACEY) then
              absorb_ymin_outer_core(i,k,ispec2D) = weight*sn
            endif

          enddo
        enddo
      enddo
!$OMP ENDDO NOWAIT
    endif

    !   ymax
    if (nspec2D_ymax_outer_core > 0) then
      ! adds boundary contributions
!$OMP DO
      do ispec2D = 1,nspec2D_ymax_outer_core

        ispec = ibelm_ymax_outer_core(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_eta_outer_core(2,ispec2D) == 0 .or. nimin_outer_core(2,ispec2D) == 0) cycle

        j=NGLLY
        do k=nkmin_eta_outer_core(2,ispec2D),NGLLZ
          do i=nimin_outer_core(2,ispec2D),nimax_outer_core(2,ispec2D)
            iglob=ibool_outer_core(i,j,k,ispec)

            sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

            weight = jacobian2D_ymax_outer_core(i,k,ispec2D)*wgllwgll_xz(i,k)

!$OMP ATOMIC
            accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

            if (SAVE_STACEY) then
              absorb_ymax_outer_core(i,k,ispec2D) = weight*sn
            endif

          enddo
        enddo
      enddo
!$OMP ENDDO NOWAIT
    endif

    ! zmin
    ! for surface elements exactly on the ICB
    if (nspec2D_zmin_outer_core > 0) then
      ! adds boundary contributions
!$OMP DO
      do ispec2D = 1,nspec2D_zmin_outer_core

        ispec = ibelm_bottom_outer_core(ispec2D)

        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)

            sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

            weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

!$OMP ATOMIC
            accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

            if (SAVE_STACEY) then
              absorb_zmin_outer_core(i,j,ispec2D) = weight*sn
            endif

          enddo
        enddo
      enddo
!$OMP ENDDO NOWAIT
    endif

!$OMP END PARALLEL

  else
    ! on GPU

    ! xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
      if (nspec2D_xmin_outer_core > 0) then
        ! adds boundary contributions
        call compute_stacey_acoustic_gpu(Mesh_pointer,absorb_xmin_outer_core,4) ! <= xmin
      endif
    endif

    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
      if (nspec2D_xmax_outer_core > 0) then
        call compute_stacey_acoustic_gpu(Mesh_pointer,absorb_xmax_outer_core,5) ! <= xmax
      endif
    endif

    !   ymin
    if (nspec2D_ymin_outer_core > 0) then
      ! adds boundary contributions
      call compute_stacey_acoustic_gpu(Mesh_pointer,absorb_ymin_outer_core,6) ! <= ymin
    endif

    !   ymax
    if (nspec2D_ymax_outer_core > 0) then
      ! adds boundary contributions
      call compute_stacey_acoustic_gpu(Mesh_pointer,absorb_ymax_outer_core,7) ! <= ymax
    endif

    ! zmin
    ! for surface elements exactly on the ICB
    if (nspec2D_zmin_outer_core > 0) then
      ! adds boundary contributions
      call compute_stacey_acoustic_gpu(Mesh_pointer,absorb_zmin_outer_core,8) ! <= zmin
    endif
  endif

  ! writes absorbing boundary values to file
  if (SAVE_STACEY) then
    ! xmin
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
      if (nspec2D_xmin_outer_core > 0) then
        call write_abs(4,absorb_xmin_outer_core,reclen_xmin_outer_core,it)
      endif
    endif

    !   xmax
    ! if two chunks exclude this face for one of them
    if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
      if (nspec2D_xmax_outer_core > 0) then
        call write_abs(5,absorb_xmax_outer_core,reclen_xmax_outer_core,it)
      endif
    endif

    !   ymin
    if (nspec2D_ymin_outer_core > 0) then
      call write_abs(6,absorb_ymin_outer_core,reclen_ymin_outer_core,it)
    endif

    !   ymax
    if (nspec2D_ymax_outer_core > 0) then
      call write_abs(7,absorb_ymax_outer_core,reclen_ymax_outer_core,it)
    endif

    ! zmin
    ! for surface elements exactly on the ICB
    if (nspec2D_zmin_outer_core > 0) then
      call write_abs(8,absorb_zmin_outer_core,reclen_zmin,it)
    endif
  endif ! SAVE_STACEY

  end subroutine compute_stacey_oc_forward


!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_stacey_oc_backward()

  use constants_solver

  use specfem_par, only: &
    ichunk,SIMULATION_TYPE,NSTEP,it

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_outercore, only: &
    b_accel_outer_core, &
    ibool_outer_core, &
    nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core,nspec2D_zmin_outer_core, &
    nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core, &
    njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core, &
    absorb_xmin_outer_core,absorb_xmax_outer_core, &
    absorb_ymin_outer_core,absorb_ymax_outer_core, &
    absorb_zmin_outer_core, &
    reclen_xmin_outer_core,reclen_xmax_outer_core, &
    reclen_ymin_outer_core,reclen_ymax_outer_core, &
    reclen_zmin, &
    ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
    ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
    ibelm_bottom_outer_core
  implicit none

  ! local parameters
  integer :: i,j,k,ispec2D,ispec,iglob

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

  !   xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
    ! reads absorbing boundary values
    if (nspec2D_xmin_outer_core > 0) then
      ! note: backward/reconstructed wavefields are read in after the Newmark time scheme in the first time loop
      !          this leads to a corresponding boundary condition at time index NSTEP - (it-1) = NSTEP - it + 1
      call read_abs(4,absorb_xmin_outer_core,reclen_xmin_outer_core,NSTEP-it+1)

      ! adds boundary contributions
      if (.not. GPU_MODE) then
        ! on CPU
        do ispec2D = 1,nspec2D_xmin_outer_core

          ispec=ibelm_xmin_outer_core(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi_outer_core(1,ispec2D) == 0 .or. njmin_outer_core(1,ispec2D) == 0) cycle

          i = 1
          do k = nkmin_xi_outer_core(1,ispec2D),NGLLZ
            do j = njmin_outer_core(1,ispec2D),njmax_outer_core(1,ispec2D)
              iglob=ibool_outer_core(i,j,k,ispec)

              b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_xmin_outer_core(j,k,ispec2D)
            enddo
          enddo
        enddo
      else
        ! on GPU
        call compute_stacey_acoustic_backward_gpu(Mesh_pointer,absorb_xmin_outer_core,4) ! <= xmin
      endif
    endif
  endif

  !   xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
    if (nspec2D_xmax_outer_core > 0) then
      ! reads in boundary contributions from file
      call read_abs(5,absorb_xmax_outer_core,reclen_xmax_outer_core,NSTEP-it+1)

      ! adds boundary contributions
      if (.not. GPU_MODE) then
        ! on CPU
        do ispec2D = 1,nspec2D_xmax_outer_core

          ispec=ibelm_xmax_outer_core(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi_outer_core(2,ispec2D) == 0 .or. njmin_outer_core(2,ispec2D) == 0) cycle

          i=NGLLX
          do k=nkmin_xi_outer_core(2,ispec2D),NGLLZ
            do j=njmin_outer_core(2,ispec2D),njmax_outer_core(2,ispec2D)
              iglob=ibool_outer_core(i,j,k,ispec)

              b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_xmax_outer_core(j,k,ispec2D)
            enddo
          enddo
        enddo
      else
        ! on GPU
        call compute_stacey_acoustic_backward_gpu(Mesh_pointer,absorb_xmax_outer_core,5) ! <= xmax
      endif
    endif
  endif

  !   ymin
  if (nspec2D_ymin_outer_core > 0) then
    ! reads in boundary contributions from file
    call read_abs(6,absorb_ymin_outer_core,reclen_ymin_outer_core,NSTEP-it+1)

    ! adds boundary contributions
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_ymin_outer_core

        ispec=ibelm_ymin_outer_core(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_eta_outer_core(1,ispec2D) == 0 .or. nimin_outer_core(1,ispec2D) == 0) cycle

        j = 1
        do k = nkmin_eta_outer_core(1,ispec2D),NGLLZ
          do i = nimin_outer_core(1,ispec2D),nimax_outer_core(1,ispec2D)
            iglob=ibool_outer_core(i,j,k,ispec)

            b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_ymin_outer_core(i,k,ispec2D)
          enddo
        enddo
      enddo
    else
      ! on GPU
      call compute_stacey_acoustic_backward_gpu(Mesh_pointer,absorb_ymin_outer_core,6) ! <= ymin
    endif
  endif

  !   ymax
  if ( nspec2D_ymax_outer_core > 0) then
    ! reads in boundary contributions from file
    call read_abs(7,absorb_ymax_outer_core,reclen_ymax_outer_core,NSTEP-it+1)

    ! adds boundary contributions
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_ymax_outer_core

        ispec=ibelm_ymax_outer_core(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_eta_outer_core(2,ispec2D) == 0 .or. nimin_outer_core(2,ispec2D) == 0) cycle

        j=NGLLY
        do k=nkmin_eta_outer_core(2,ispec2D),NGLLZ
          do i=nimin_outer_core(2,ispec2D),nimax_outer_core(2,ispec2D)
            iglob=ibool_outer_core(i,j,k,ispec)

            b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_ymax_outer_core(i,k,ispec2D)
          enddo
        enddo
      enddo
    else
      ! on GPU
      call compute_stacey_acoustic_backward_gpu(Mesh_pointer,absorb_ymax_outer_core,7) ! <= ymax
    endif
  endif

  ! zmin

  ! for surface elements exactly on the ICB
  if (nspec2D_zmin_outer_core > 0) then
    ! reads in boundary contributions from file
    call read_abs(8,absorb_zmin_outer_core,reclen_zmin,NSTEP-it+1)

    ! adds boundary contributions
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_zmin_outer_core

        ispec = ibelm_bottom_outer_core(ispec2D)

        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)

            b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_zmin_outer_core(i,j,ispec2D)
          enddo
        enddo
      enddo
    else
      ! on GPU
      call compute_stacey_acoustic_backward_gpu(Mesh_pointer,absorb_zmin_outer_core,8) ! <= zmin
    endif
  endif

  end subroutine compute_stacey_oc_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_stacey_oc_backward_undoatt()

  use constants_solver

  use specfem_par, only: &
    ichunk,SIMULATION_TYPE,SAVE_FORWARD,UNDO_ATTENUATION, &
    wgllwgll_xz,wgllwgll_yz,wgllwgll_xy

  use specfem_par, only: GPU_MODE,Mesh_pointer

  use specfem_par_outercore, only: &
    b_veloc_outer_core,b_accel_outer_core, &
    ibool_outer_core, &
    jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
    jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
    jacobian2D_bottom_outer_core, &
    nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core,nspec2D_zmin_outer_core, &
    vp_outer_core, &
    nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core, &
    njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core, &
    ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
    ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
    ibelm_bottom_outer_core
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: sn,weight
  integer :: i,j,k,ispec2D,ispec,iglob

  ! checks
  if (SIMULATION_TYPE /= 3 ) return
  if (SAVE_FORWARD) return
  if (.not. UNDO_ATTENUATION) stop 'Error invalid UNDO_ATTENUATION flag for compute_stacey_oc_backward_undoatt()'

  ! outer core

  !   xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then
    if (nspec2D_xmin_outer_core > 0) then
      ! adds boundary contribution
      if (.not. GPU_MODE) then
        ! on CPU
        do ispec2D = 1,nspec2D_xmin_outer_core

          ispec=ibelm_xmin_outer_core(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi_outer_core(1,ispec2D) == 0 .or. njmin_outer_core(1,ispec2D) == 0) cycle

          i = 1
          do k = nkmin_xi_outer_core(1,ispec2D),NGLLZ
            do j = njmin_outer_core(1,ispec2D),njmax_outer_core(1,ispec2D)
              iglob=ibool_outer_core(i,j,k,ispec)

              sn = b_veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

              weight = jacobian2D_xmin_outer_core(j,k,ispec2D)*wgllwgll_yz(j,k)

              b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - weight*sn

            enddo
          enddo
        enddo
      else
        ! on GPU
        call compute_stacey_acoustic_undoatt_gpu(Mesh_pointer,4) ! <= xmin
      endif
    endif
  endif

  !   xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then
    if (nspec2D_xmax_outer_core > 0) then
      ! adds boundary contribution
      if (.not. GPU_MODE) then
        ! on CPU
        do ispec2D = 1,nspec2D_xmax_outer_core

          ispec=ibelm_xmax_outer_core(ispec2D)

          ! exclude elements that are not on absorbing edges
          if (nkmin_xi_outer_core(2,ispec2D) == 0 .or. njmin_outer_core(2,ispec2D) == 0) cycle

          i=NGLLX
          do k=nkmin_xi_outer_core(2,ispec2D),NGLLZ
            do j=njmin_outer_core(2,ispec2D),njmax_outer_core(2,ispec2D)
              iglob=ibool_outer_core(i,j,k,ispec)

              sn = b_veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

              weight = jacobian2D_xmax_outer_core(j,k,ispec2D)*wgllwgll_yz(j,k)

              b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - weight*sn

            enddo
          enddo
        enddo
      else
        ! on GPU
        call compute_stacey_acoustic_undoatt_gpu(Mesh_pointer,5) ! <= xmax
      endif
    endif
  endif

  !   ymin
  if (nspec2D_ymin_outer_core > 0) then
    ! adds boundary contribution
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_ymin_outer_core

        ispec=ibelm_ymin_outer_core(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_eta_outer_core(1,ispec2D) == 0 .or. nimin_outer_core(1,ispec2D) == 0) cycle

        j = 1
        do k = nkmin_eta_outer_core(1,ispec2D),NGLLZ
          do i = nimin_outer_core(1,ispec2D),nimax_outer_core(1,ispec2D)
            iglob=ibool_outer_core(i,j,k,ispec)

            sn = b_veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

            weight=jacobian2D_ymin_outer_core(i,k,ispec2D)*wgllwgll_xz(i,k)

            b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - weight*sn

          enddo
        enddo
      enddo
    else
      ! on GPU
      call compute_stacey_acoustic_undoatt_gpu(Mesh_pointer,6) ! <= ymin
    endif
  endif

  !   ymax
  if (nspec2D_ymax_outer_core > 0) then
    ! adds boundary contribution
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_ymax_outer_core

        ispec=ibelm_ymax_outer_core(ispec2D)

        ! exclude elements that are not on absorbing edges
        if (nkmin_eta_outer_core(2,ispec2D) == 0 .or. nimin_outer_core(2,ispec2D) == 0) cycle

        j=NGLLY
        do k=nkmin_eta_outer_core(2,ispec2D),NGLLZ
          do i=nimin_outer_core(2,ispec2D),nimax_outer_core(2,ispec2D)
            iglob=ibool_outer_core(i,j,k,ispec)

            sn = b_veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

            weight=jacobian2D_ymax_outer_core(i,k,ispec2D)*wgllwgll_xz(i,k)

            b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - weight*sn

          enddo
        enddo
      enddo
    else
      ! on GPU
      call compute_stacey_acoustic_undoatt_gpu(Mesh_pointer,7) ! <= ymax
    endif
  endif

  ! zmin

  ! for surface elements exactly on the ICB
  if (nspec2D_zmin_outer_core > 0) then
    ! adds boundary contribution
    if (.not. GPU_MODE) then
      ! on CPU
      do ispec2D = 1,nspec2D_zmin_outer_core

        ispec = ibelm_bottom_outer_core(ispec2D)

        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool_outer_core(i,j,k,ispec)

            sn = b_veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

            weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

            b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - weight*sn

          enddo
        enddo
      enddo
    else
      ! on GPU
      call compute_stacey_acoustic_undoatt_gpu(Mesh_pointer,8) ! <= zmin
    endif
  endif

  end subroutine compute_stacey_oc_backward_undoatt


