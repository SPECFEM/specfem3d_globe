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

  subroutine compute_forces_acoustic()

! acoustic domains for forward or adjoint simulations (SIMULATION_TYPE == 1 or 2 )

  use specfem_par
  use specfem_par_outercore
  use specfem_par_movie, only: div_displ_outer_core
  use specfem_par_full_gravity, only: pgrav_oc

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: timeval
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements in the outer_core,
  !              iphase = 2 is for computing inner elements in the outer core (former icall parameter)
  integer :: iphase

  ! checks if anything to do
  ! for regional mesh cut-offs, there are no outer core elements
  if (NSPEC_OUTER_CORE == 0) return

  ! compute internal forces in the fluid region

  ! current simulated time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing
    !       for the source.
    timeval = real((dble(it-1-1)*DT + dble(C_LDDRK(istage))*DT - t0)*scale_t_inv, kind=CUSTOM_REAL)
  else
    timeval = real((dble(it-1)*DT - t0)*scale_t_inv, kind=CUSTOM_REAL)
  endif

  ! ****************************************************
  !   big loop over all spectral elements in the fluid
  ! ****************************************************

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! first, iphase == 1 for points on MPI interfaces (thus outer elements)
    ! second, iphase == 2 for points purely inside partition (thus inner elements)
    !
    ! compute all the outer elements first, then sends out non blocking MPI communication
    ! and continues computing inner elements (overlapping)

    if (.not. GPU_MODE) then
      ! on CPU
      call compute_forces_outer_core(timeval,deltat,two_omega_earth, &
                                     NSPEC_OUTER_CORE_ROTATION,NGLOB_OUTER_CORE, &
                                     A_array_rotation,B_array_rotation, &
                                     A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                     displ_outer_core,accel_outer_core, &
                                     div_displ_outer_core,iphase, &
                                     pgrav_oc)
    else
      ! on GPU
      ! includes FORWARD_OR_ADJOINT == 1
      call compute_forces_outer_core_gpu(Mesh_pointer,iphase,timeval,ALPHA_LDDRK(istage),BETA_LDDRK(istage),1)

      ! initiates asynchronous MPI transfer
      if (NPROCTOT_VAL > 1) then
        if (GPU_ASYNC_COPY .and. iphase == 2) then
          ! crust/mantle region
          ! wait for asynchronous copy to finish
          call sync_copy_from_device(Mesh_pointer,iphase,buffer_send_scalar_outer_core,IREGION_OUTER_CORE,1)
          ! sends MPI buffers
          call assemble_MPI_scalar_send_gpu(NPROCTOT_VAL, &
                                            buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
                                            num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                            nibool_interfaces_outer_core, &
                                            my_neighbors_outer_core, &
                                            request_send_scalar_oc,request_recv_scalar_oc)
        endif
      endif
    endif

    ! computes additional contributions to acceleration field
    if (iphase == 1) then

       ! Stacey absorbing boundaries
       if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) call compute_stacey_oc_forward()

       ! ****************************************************
       ! **********  add matching with solid part  **********
       ! ****************************************************
       call compute_coupling_fluid()

    endif ! iphase == 1

    ! assemble all the contributions between slices using MPI
    if (NPROCTOT_VAL > 1) then
      ! in outer core
      if (iphase == 1) then
        ! sends out MPI interface data (non-blocking)
        if (.not. GPU_MODE) then
          ! on CPU
          call assemble_MPI_scalar_s(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                     accel_outer_core, &
                                     buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
                                     num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                     nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                                     my_neighbors_outer_core, &
                                     request_send_scalar_oc,request_recv_scalar_oc)
        else
          ! on GPU
          ! sends accel values to corresponding MPI interface neighbors

          ! preparation of the contribution between partitions using MPI
          ! transfers MPI buffers to CPU
          ! note: for asynchronous transfers, this transfers boundary region to host asynchronously. The
          !       MPI-send is done after compute_forces_outer_core_gpu,
          !       once the inner element kernels are launched, and the memcpy has finished.
          call transfer_boun_pot_from_device(Mesh_pointer,buffer_send_scalar_outer_core,1)

          if (.not. GPU_ASYNC_COPY) then
            ! for synchronous transfers, sending over MPI can directly proceed
            ! outer core
            call assemble_MPI_scalar_send_gpu(NPROCTOT_VAL, &
                                              buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
                                              num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                              nibool_interfaces_outer_core, &
                                              my_neighbors_outer_core, &
                                              request_send_scalar_oc,request_recv_scalar_oc)
          endif
        endif
      else
        ! make sure the last communications are finished and processed
        ! waits for send/receive requests to be completed and assembles values
        if (.not. GPU_MODE) then
          ! on CPU
          call assemble_MPI_scalar_w(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                     accel_outer_core, &
                                     buffer_recv_scalar_outer_core,num_interfaces_outer_core, &
                                     max_nibool_interfaces_oc, &
                                     nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                                     my_neighbors_outer_core, &
                                     request_send_scalar_oc,request_recv_scalar_oc)
        else
          ! on GPU
          if (GPU_ASYNC_COPY) then
            ! while inner elements compute "Kernel_2", we wait for MPI to
            ! finish and transfer the boundary terms to the device asynchronously
            !
            ! transfers MPI buffers onto GPU
            call transfer_boundarypot_to_device(Mesh_pointer,NPROCTOT_VAL,buffer_recv_scalar_outer_core, &
                                                num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                                request_recv_scalar_oc, &
                                                IREGION_OUTER_CORE,1)
          endif

          ! waits for MPI send/receive requests to be completed and assembles values
          call assemble_MPI_scalar_write_gpu(Mesh_pointer,NPROCTOT_VAL, &
                                             buffer_recv_scalar_outer_core, &
                                             num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                             request_send_scalar_oc,request_recv_scalar_oc, &
                                             1) ! -- 1 == fwd accel
        endif
      endif ! iphase == 1
    endif

  enddo ! iphase

  ! multiply by the inverse of the mass matrix
  if (.not. GPU_MODE) then
    ! on CPU
    call multiply_accel_acoustic(NGLOB_OUTER_CORE,accel_outer_core,rmass_outer_core)
  else
    ! on GPU
    ! includes FORWARD_OR_ADJOINT == 1
    call multiply_accel_acoustic_gpu(Mesh_pointer,1)
  endif

  ! time schemes
  if (USE_LDDRK) then
    ! Runge-Kutta scheme
    call update_veloc_acoustic_lddrk()
  else
    ! Newmark time scheme
    call update_veloc_acoustic_newmark()
  endif

  end subroutine compute_forces_acoustic

!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_forces_acoustic_backward()

! backward/reconstructed wavefields only

  use specfem_par
  use specfem_par_outercore
  use specfem_par_movie, only: div_displ_outer_core
  use specfem_par_full_gravity, only: b_pgrav_oc

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: b_timeval
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements in the outer_core,
  !              iphase = 2 is for computing inner elements in the outer core (former icall parameter)
  integer :: iphase
  integer :: it_tmp

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

  ! compute internal forces in the fluid region

  ! note on backward/reconstructed wavefields:
  !       b_timeval for b_displ( it=1 ) corresponds to (NSTEP - 1)*DT - t0  (after Newmark scheme...)
  !       as we start with saved wavefields b_displ( 1 ), displ( NSTEP ) which correspond
  !       to a time (NSTEP - (it-1) - 1)*DT - t0
  !       for reconstructing the rotational contributions
  if (UNDO_ATTENUATION) then
    !valid for multiples only:
    !it_tmp = iteration_on_subset * NT_DUMP_ATTENUATION - it_of_this_subset + 1
    !
    ! example: NSTEP is **NOT** a multiple of NT_DUMP_ATTENUATION
    !          NT_DUMP_ATTENUATION = 301, NSTEP = 900, NSUBSET_ITERATIONS = 3, iteration_on_subset = 1 -> 3
    !              1. subset, it_temp goes from (900 - 602) = 298 down to 1
    !              2. subset, it_temp goes from (900 - 301) = 599 down to 299
    !              3. subset, it_temp goes from (900 - 0)   = 900 down to 600
    !works always:
    it_tmp = NSTEP - (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION - it_of_this_subset + 1
  else
    it_tmp = it
  endif

  ! current simulated time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing
    !       for the source.
    if (UNDO_ATTENUATION) then
      ! stepping moves forward from snapshot position
      b_timeval = real((dble(NSTEP-it_tmp-1)*DT + dble(C_LDDRK(istage))*DT - t0)*scale_t_inv, kind=CUSTOM_REAL)
    else
      ! stepping backwards
      b_timeval = real((dble(NSTEP-it_tmp-1)*DT - dble(C_LDDRK(istage))*DT - t0)*scale_t_inv, kind=CUSTOM_REAL)
    endif
  else
    b_timeval = real((dble(NSTEP-it_tmp)*DT - t0)*scale_t_inv, kind=CUSTOM_REAL)
  endif

  !debug
  !if (myrank == 0 ) print *,'compute_forces_acoustic_backward: it = ',it_tmp

  ! ****************************************************
  !   big loop over all spectral elements in the fluid
  ! ****************************************************

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase = 1,2

    ! first, iphase == 1 for points on MPI interfaces (thus outer elements)
    ! second, iphase == 2 for points purely inside partition (thus inner elements)
    !
    ! compute all the outer elements first, then sends out non blocking MPI communication
    ! and continues computing inner elements (overlapping)

    if (.not. GPU_MODE) then
      ! on CPU
      ! adjoint / kernel runs
      call compute_forces_outer_core(b_timeval,b_deltat,b_two_omega_earth, &
                                     NSPEC_OUTER_CORE_ROT_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                                     b_A_array_rotation,b_B_array_rotation, &
                                     b_A_array_rotation_lddrk,b_B_array_rotation_lddrk, &
                                     b_displ_outer_core,b_accel_outer_core, &
                                     div_displ_outer_core,iphase, &
                                     b_pgrav_oc)
    else
      ! on GPU
      ! includes FORWARD_OR_ADJOINT == 3
      call compute_forces_outer_core_gpu(Mesh_pointer,iphase,b_timeval,ALPHA_LDDRK(istage),BETA_LDDRK(istage),3)

      ! initiates asynchronous MPI transfer
      if (GPU_ASYNC_COPY .and. iphase == 2) then
        ! crust/mantle region
        ! wait for asynchronous copy to finish
        call sync_copy_from_device(Mesh_pointer,iphase,b_buffer_send_scalar_outer_core,IREGION_OUTER_CORE,3)
        ! sends MPI buffers
        call assemble_MPI_scalar_send_gpu(NPROCTOT_VAL, &
                              b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core, &
                              num_interfaces_outer_core,max_nibool_interfaces_oc, &
                              nibool_interfaces_outer_core, &
                              my_neighbors_outer_core, &
                              b_request_send_scalar_oc,b_request_recv_scalar_oc)
      endif

    endif


    ! computes additional contributions to acceleration field
    if (iphase == 1) then

      ! Stacey absorbing boundaries
      if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
        if (UNDO_ATTENUATION) then
          call compute_stacey_oc_backward_undoatt()
        else
          call compute_stacey_oc_backward()
        endif
      endif

      ! ****************************************************
      ! **********  add matching with solid part  **********
      ! ****************************************************
      call compute_coupling_fluid_backward()

    endif ! iphase == 1

    ! assemble all the contributions between slices using MPI
    ! in outer core
    if (iphase == 1) then
      ! sends out MPI interface data (non-blocking)
      ! adjoint simulations
      if (.not. GPU_MODE) then
        ! on CPU
        call assemble_MPI_scalar_s(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                              b_accel_outer_core, &
                              b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core, &
                              num_interfaces_outer_core,max_nibool_interfaces_oc, &
                              nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                              my_neighbors_outer_core, &
                              b_request_send_scalar_oc,b_request_recv_scalar_oc)
      else
        ! on GPU
        ! sends accel values to corresponding MPI interface neighbors

        ! preparation of the contribution between partitions using MPI
        ! transfers MPI buffers to CPU
        ! note: for asynchronous transfers, this transfers boundary region to host asynchronously. The
        !       MPI-send is done after compute_forces_outer_core_gpu,
        !       once the inner element kernels are launched, and the memcpy has finished.
        call transfer_boun_pot_from_device(Mesh_pointer,b_buffer_send_scalar_outer_core,3)

        if (.not. GPU_ASYNC_COPY) then
          ! for synchronous transfers, sending over MPI can directly proceed
          ! outer core
          call assemble_MPI_scalar_send_gpu(NPROCTOT_VAL, &
                                b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                nibool_interfaces_outer_core, &
                                my_neighbors_outer_core, &
                                b_request_send_scalar_oc,b_request_recv_scalar_oc)
        endif
      endif ! GPU
    else
      ! make sure the last communications are finished and processed
      ! waits for send/receive requests to be completed and assembles values
      ! adjoint simulations
      if (.not. GPU_MODE) then
        ! on CPU
        call assemble_MPI_scalar_w(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                   b_accel_outer_core, &
                                   b_buffer_recv_scalar_outer_core,num_interfaces_outer_core, &
                                   max_nibool_interfaces_oc, &
                                   nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                                   my_neighbors_outer_core, &
                                   b_request_send_scalar_oc,b_request_recv_scalar_oc)
      else
        ! on GPU
        if (GPU_ASYNC_COPY) then
          ! while inner elements compute "Kernel_2", we wait for MPI to
          ! finish and transfer the boundary terms to the device asynchronously
          !
          ! transfers MPI buffers onto GPU
          call transfer_boundarypot_to_device(Mesh_pointer,NPROCTOT_VAL,b_buffer_recv_scalar_outer_core, &
                                              num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                              b_request_recv_scalar_oc, &
                                              IREGION_OUTER_CORE,3)
        endif

        ! waits for MPI send/receive requests to be completed and assembles values
        call assemble_MPI_scalar_write_gpu(Mesh_pointer,NPROCTOT_VAL, &
                                            b_buffer_recv_scalar_outer_core, &
                                            num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                            b_request_send_scalar_oc,b_request_recv_scalar_oc, &
                                            3) ! -- 3 == adjoint b_accel
      endif
    endif ! iphase == 1

  enddo ! iphase

  ! multiply by the inverse of the mass matrix
  if (.not. GPU_MODE) then
    ! on CPU
    call multiply_accel_acoustic(NGLOB_OUTER_CORE_ADJOINT,b_accel_outer_core,b_rmass_outer_core)
  else
    ! on GPU
    ! includes FORWARD_OR_ADJOINT == 3
    call multiply_accel_acoustic_gpu(Mesh_pointer,3)
  endif

  ! time schemes
  if (USE_LDDRK) then
    ! Runge-Kutta scheme
    call update_veloc_acoustic_lddrk_backward()
  else
    ! Newmark time scheme
    call update_veloc_acoustic_newmark_backward()
  endif

  end subroutine compute_forces_acoustic_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_forces_outer_core(timeval,deltat,two_omega_earth, &
                                       NSPEC_ROT,NGLOB, &
                                       A_array_rotation,B_array_rotation, &
                                       A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                       displfluid,accelfluid, &
                                       div_displfluid,iphase, &
                                       pgrav_outer_core)

! wrapper function, decides about Deville optimization
!
! (left in this file to let compiler decide about inlining)

  use constants_solver, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE,USE_DEVILLE_PRODUCTS_VAL

  use specfem_par_outercore, only: sum_terms_outer_core

  implicit none

  integer,intent(in) :: NSPEC_ROT,NGLOB

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL),intent(in) :: timeval,deltat,two_omega_earth

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ROT),intent(inout) :: &
    A_array_rotation,B_array_rotation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ROT),intent(inout) :: &
    A_array_rotation_lddrk,B_array_rotation_lddrk

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: displfluid
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(inout) :: accelfluid

  ! divergence of displacement
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE),intent(out) :: div_displfluid

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: pgrav_outer_core

  if (USE_DEVILLE_PRODUCTS_VAL) then
    ! uses Deville et al. (2002) routine
    call compute_forces_outer_core_Dev(timeval,deltat,two_omega_earth, &
                                       NSPEC_ROT,NGLOB, &
                                       A_array_rotation,B_array_rotation, &
                                       A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                       displfluid,accelfluid, &
                                       div_displfluid,iphase,sum_terms_outer_core, &
                                       pgrav_outer_core)
  else
    ! div_displ_outer_core is initialized to zero in the following subroutine.
    call compute_forces_outer_core_noDev(timeval,deltat,two_omega_earth, &
                                         NSPEC_ROT,NGLOB, &
                                         A_array_rotation,B_array_rotation, &
                                         A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                         displfluid,accelfluid, &
                                         div_displfluid,iphase, &
                                         pgrav_outer_core)
  endif

  end subroutine compute_forces_outer_core

