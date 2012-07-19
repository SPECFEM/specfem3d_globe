!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine compute_forces_acoustic()

  use specfem_par
  use specfem_par_crustmantle,only: displ_crust_mantle,b_displ_crust_mantle, &
                                    ibool_crust_mantle,ibelm_bottom_crust_mantle
  use specfem_par_innercore,only: displ_inner_core,b_displ_inner_core, &
                                  ibool_inner_core,ibelm_top_inner_core
  use specfem_par_outercore
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: time,b_time
  ! non blocking MPI
  ! iphase: iphase = 1 is for computing outer elements in the outer_core,
  !              iphase = 2 is for computing inner elements in the outer core (former icall parameter)
  integer :: iphase
  logical :: phase_is_inner

  ! compute internal forces in the fluid region
  if(CUSTOM_REAL == SIZE_REAL) then
    time = sngl((dble(it-1)*DT-t0)*scale_t_inv)
  else
    time = (dble(it-1)*DT-t0)*scale_t_inv
  endif
  if (SIMULATION_TYPE == 3) then
    ! note on backward/reconstructed wavefields:
    !       b_time for b_displ( it=1 ) corresponds to (NSTEP - 1)*DT - t0  (after Newmark scheme...)
    !       as we start with saved wavefields b_displ( 1 ) <-> displ( NSTEP ) which correspond
    !       to a time (NSTEP - (it-1) - 1)*DT - t0
    !       for reconstructing the rotational contributions
    if(CUSTOM_REAL == SIZE_REAL) then
      b_time = sngl((dble(NSTEP-it)*DT-t0)*scale_t_inv)
    else
      b_time = (dble(NSTEP-it)*DT-t0)*scale_t_inv
    endif
  endif

  ! ****************************************************
  !   big loop over all spectral elements in the fluid
  ! ****************************************************

  ! distinguishes two runs: for elements on MPI interfaces, and elements within the partitions
  do iphase=1,2

    ! first, iphase == 1 for points on MPI interfaces (thus outer elements)
    ! second, iphase == 2 for points purely inside partition (thus inner elements)
    !
    ! compute all the outer elements first, then sends out non blocking MPI communication
    ! and continues computing inner elements (overlapping)
    if( iphase == 1 ) then
      phase_is_inner = .false.
    else
      phase_is_inner = .true.
    endif

    if( .not. GPU_MODE ) then
      ! on CPU
      if( USE_DEVILLE_PRODUCTS_VAL ) then
        ! uses Deville et al. (2002) routine
        call compute_forces_outer_core_Dev(time,deltat,two_omega_earth, &
                                      NSPEC_OUTER_CORE_ROTATION,NGLOB_OUTER_CORE, &
                                      A_array_rotation,B_array_rotation, &
                                      displ_outer_core,accel_outer_core, &
                                      div_displ_outer_core,phase_is_inner)
      else
        ! div_displ_outer_core is initialized to zero in the following subroutine.
        call compute_forces_outer_core(time,deltat,two_omega_earth, &
                                      NSPEC_OUTER_CORE_ROTATION,NGLOB_OUTER_CORE, &
                                      A_array_rotation,B_array_rotation, &
                                      displ_outer_core,accel_outer_core, &
                                      div_displ_outer_core,phase_is_inner)
      endif

      ! adjoint / kernel runs
      if (SIMULATION_TYPE == 3) then
        if( USE_DEVILLE_PRODUCTS_VAL ) then
          ! uses Deville et al. (2002) routine
          call compute_forces_outer_core_Dev(b_time,b_deltat,b_two_omega_earth, &
                                      NSPEC_OUTER_CORE_ROT_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                                      b_A_array_rotation,b_B_array_rotation, &
                                      b_displ_outer_core,b_accel_outer_core, &
                                      b_div_displ_outer_core,phase_is_inner)
        else
          call compute_forces_outer_core(b_time,b_deltat,b_two_omega_earth, &
                                      NSPEC_OUTER_CORE_ROT_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                                      b_A_array_rotation,b_B_array_rotation, &
                                      b_displ_outer_core,b_accel_outer_core, &
                                      b_div_displ_outer_core,phase_is_inner)
        endif
      endif

    else
      ! on GPU
      ! includes both forward and adjoint/kernel simulations
      call compute_forces_outer_core_cuda(Mesh_pointer,iphase,time,b_time)
    endif


    ! computes additional contributions to acceleration field
    if( iphase == 1 ) then

       ! Stacey absorbing boundaries
       if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) call compute_stacey_outer_core()

       ! ****************************************************
       ! **********  add matching with solid part  **********
       ! ****************************************************
       ! only for elements in first matching layer in the fluid
       if( .not. GPU_MODE ) then
          ! on CPU
          !---
          !--- couple with mantle at the top of the outer core
          !---
          if(ACTUALLY_COUPLE_FLUID_CMB) &
               call compute_coupling_fluid_CMB(displ_crust_mantle,b_displ_crust_mantle, &
               ibool_crust_mantle,ibelm_bottom_crust_mantle,  &
               accel_outer_core,b_accel_outer_core, &
               normal_top_outer_core,jacobian2D_top_outer_core, &
               wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
               SIMULATION_TYPE,NSPEC2D_TOP(IREGION_OUTER_CORE))

          !---
          !--- couple with inner core at the bottom of the outer core
          !---
          if(ACTUALLY_COUPLE_FLUID_ICB) &
               call compute_coupling_fluid_ICB(displ_inner_core,b_displ_inner_core, &
               ibool_inner_core,ibelm_top_inner_core,  &
               accel_outer_core,b_accel_outer_core, &
               normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
               wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
               SIMULATION_TYPE,NSPEC2D_BOTTOM(IREGION_OUTER_CORE))

       else
          ! on GPU
          !---
          !--- couple with mantle at the top of the outer core
          !---
          if( ACTUALLY_COUPLE_FLUID_CMB ) &
               call compute_coupling_fluid_cmb_cuda(Mesh_pointer)
          !---
          !--- couple with inner core at the bottom of the outer core
          !---
          if( ACTUALLY_COUPLE_FLUID_ICB ) &
               call compute_coupling_fluid_icb_cuda(Mesh_pointer)

       endif
    endif ! iphase == 1

    ! assemble all the contributions between slices using MPI
    ! in outer core
    if( iphase == 1 ) then
      ! sends out MPI interface data (non-blocking)

      if(.NOT. GPU_MODE) then
        ! on CPU
        call assemble_MPI_scalar_ext_mesh_s(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                accel_outer_core, &
                                buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,ibool_interfaces_outer_core,&
                                my_neighbours_outer_core, &
                                request_send_scalar_outer_core,request_recv_scalar_outer_core)
      else
        ! on GPU
        ! outer core
        call assemble_MPI_scalar_send_cuda(NPROCTOT_VAL, &
                                buffer_send_scalar_outer_core,buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,&
                                my_neighbours_outer_core, &
                                request_send_scalar_outer_core,request_recv_scalar_outer_core, &
                                1) ! <-- 1 == fwd accel
      endif

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        if(.NOT. GPU_MODE) then
          ! on CPU
          call assemble_MPI_scalar_ext_mesh_s(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                b_accel_outer_core, &
                                b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,ibool_interfaces_outer_core,&
                                my_neighbours_outer_core, &
                                b_request_send_scalar_outer_core,b_request_recv_scalar_outer_core)
        else
          ! on GPU
          ! outer core
          call assemble_MPI_scalar_send_cuda(NPROCTOT_VAL, &
                                b_buffer_send_scalar_outer_core,b_buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,&
                                my_neighbours_outer_core, &
                                b_request_send_scalar_outer_core,b_request_recv_scalar_outer_core, &
                                3) ! <-- 3 == adjoint b_accel
        endif ! GPU
      endif ! SIMULATION_TYPE == 3

    else
      ! make sure the last communications are finished and processed
      ! waits for send/receive requests to be completed and assembles values
      if(.NOT. GPU_MODE) then
        ! on CPU
        call assemble_MPI_scalar_ext_mesh_w(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                accel_outer_core, &
                                buffer_recv_scalar_outer_core,num_interfaces_outer_core,&
                                max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                                request_send_scalar_outer_core,request_recv_scalar_outer_core)
      else
        ! on GPU
        call assemble_MPI_scalar_write_cuda(NPROCTOT_VAL, &
                                buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                request_send_scalar_outer_core,request_recv_scalar_outer_core, &
                                1) ! <-- 1 == fwd accel
      endif

      ! adjoint simulations
      if( SIMULATION_TYPE == 3 ) then
        if(.NOT. GPU_MODE) then
          ! on CPU
          call assemble_MPI_scalar_ext_mesh_w(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                                b_accel_outer_core, &
                                b_buffer_recv_scalar_outer_core,num_interfaces_outer_core,&
                                max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                                b_request_send_scalar_outer_core,b_request_recv_scalar_outer_core)
        else
          ! on GPU
          call assemble_MPI_scalar_write_cuda(NPROCTOT_VAL, &
                                b_buffer_recv_scalar_outer_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                b_request_send_scalar_outer_core,b_request_recv_scalar_outer_core, &
                                3) ! <-- 3 == adjoint b_accel
        endif
      endif ! SIMULATION_TYPE == 3
    endif ! iphase == 1

  enddo ! iphase

  ! Newmark time scheme:
  ! corrector terms for fluid parts
  ! (multiply by the inverse of the mass matrix and update velocity)
  if(.NOT. GPU_MODE) then
    ! on CPU
    call compute_forces_ac_update_veloc(NGLOB_OUTER_CORE,veloc_outer_core,accel_outer_core, &
                                       deltatover2,rmass_outer_core)

   ! adjoint / kernel runs
    if (SIMULATION_TYPE == 3) &
      call compute_forces_ac_update_veloc(NGLOB_OUTER_CORE_ADJOINT,b_veloc_outer_core,b_accel_outer_core, &
                                         b_deltatover2,rmass_outer_core)

  else
    ! on GPU
    call kernel_3_outer_core_cuda(Mesh_pointer, &
                                 deltatover2,SIMULATION_TYPE,b_deltatover2)
  endif

  end subroutine compute_forces_acoustic

!=====================================================================

  subroutine compute_forces_ac_update_veloc(NGLOB,veloc_outer_core,accel_outer_core, &
                                           deltatover2,rmass_outer_core)

  use constants_solver,only: CUSTOM_REAL

#ifdef _HANDOPT
  use specfem_par,only: imodulo_NGLOB_OUTER_CORE
#endif

  implicit none

  integer :: NGLOB

  ! velocity potential
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: veloc_outer_core,accel_outer_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: rmass_outer_core

  real(kind=CUSTOM_REAL) :: deltatover2

  ! local parameters
  integer :: i

  ! Newmark time scheme
  ! multiply by the inverse of the mass matrix and update velocity

#ifdef _HANDOPT_NEWMARK
! way 2:
  ! outer core
  if(imodulo_NGLOB_OUTER_CORE >= 1) then
    do i=1,imodulo_NGLOB_OUTER_CORE
      accel_outer_core(i) = accel_outer_core(i)*rmass_outer_core(i)
      veloc_outer_core(i) = veloc_outer_core(i) + deltatover2*accel_outer_core(i)
    enddo
  endif
  do i=imodulo_NGLOB_OUTER_CORE+1,NGLOB,3
    accel_outer_core(i) = accel_outer_core(i)*rmass_outer_core(i)
    veloc_outer_core(i) = veloc_outer_core(i) + deltatover2*accel_outer_core(i)

    accel_outer_core(i+1) = accel_outer_core(i+1)*rmass_outer_core(i+1)
    veloc_outer_core(i+1) = veloc_outer_core(i+1) + deltatover2*accel_outer_core(i+1)

    accel_outer_core(i+2) = accel_outer_core(i+2)*rmass_outer_core(i+2)
    veloc_outer_core(i+2) = veloc_outer_core(i+2) + deltatover2*accel_outer_core(i+2)
  enddo
#else
! way 1:
  do i=1,NGLOB
    accel_outer_core(i) = accel_outer_core(i)*rmass_outer_core(i)
    veloc_outer_core(i) = veloc_outer_core(i) + deltatover2*accel_outer_core(i)
  enddo
#endif

  end subroutine compute_forces_ac_update_veloc

