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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


  subroutine prepare_wavefields()

! initializes arrays
!
! note: we first allocate arrays, and initialization will be done in extra routine.
!       This will enable us to call the initialization routine more than once in case we want to loop over different
!       earthquakes in the same simulation setup

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: ier

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing wavefields"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  allocating wavefields"
    call flush_IMAIN()
  endif

  ! allocates arrays
  allocate(displ_crust_mantle(NDIM,NGLOB_CRUST_MANTLE), &
           veloc_crust_mantle(NDIM,NGLOB_CRUST_MANTLE), &
           accel_crust_mantle(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating displ,veloc,accel in crust_mantle'

  allocate(displ_outer_core(NGLOB_OUTER_CORE), &
           veloc_outer_core(NGLOB_OUTER_CORE), &
           accel_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating displ,veloc,accel in outer_core'

  allocate(displ_inner_core(NDIM,NGLOB_INNER_CORE), &
           veloc_inner_core(NDIM,NGLOB_INNER_CORE), &
           accel_inner_core(NDIM,NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating displ,veloc,accel in inner_core'

  ! for strain/attenuation
  allocate(epsilondev_xx_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT), &
           epsilondev_yy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT), &
           epsilondev_xy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT), &
           epsilondev_xz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT), &
           epsilondev_yz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT), &
           eps_trace_over_3_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays epsilondev_xx_crust_mantle,..'

  allocate(epsilondev_xx_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT), &
           epsilondev_yy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT), &
           epsilondev_xy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT), &
           epsilondev_xz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT), &
           epsilondev_yz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT), &
           eps_trace_over_3_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays epsilondev_xx_inner_core,..'

  ! backward/reconstructed strain fields
  if (SIMULATION_TYPE == 3) then
    if (UNDO_ATTENUATION) then
      ! for undo_attenuation, whenever strain is needed it will be computed locally.
      ! pointers are using the allocated arrays for adjoint strain, however values stored in those arrays will be overwritten
      ! crust/mantle
      b_epsilondev_xx_crust_mantle => epsilondev_xx_crust_mantle
      b_epsilondev_yy_crust_mantle => epsilondev_yy_crust_mantle
      b_epsilondev_xy_crust_mantle => epsilondev_xy_crust_mantle
      b_epsilondev_xz_crust_mantle => epsilondev_xz_crust_mantle
      b_epsilondev_yz_crust_mantle => epsilondev_yz_crust_mantle
      b_eps_trace_over_3_crust_mantle => eps_trace_over_3_crust_mantle
      ! inner core
      b_epsilondev_xx_inner_core => epsilondev_xx_inner_core
      b_epsilondev_yy_inner_core => epsilondev_yy_inner_core
      b_epsilondev_xy_inner_core => epsilondev_xy_inner_core
      b_epsilondev_xz_inner_core => epsilondev_xz_inner_core
      b_epsilondev_yz_inner_core => epsilondev_yz_inner_core
      b_eps_trace_over_3_inner_core => eps_trace_over_3_inner_core
    else
      ! allocates actual arrays
      allocate(b_epsilondev_xx_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_yy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_xy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_xz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_yz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_eps_trace_over_3_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_epsilondev*** arrays for crust/mantle')

      allocate(b_epsilondev_xx_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_yy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_xy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_xz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_yz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_eps_trace_over_3_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_epsilondev*** arrays for inner core')
    endif
  else
    ! initializes pointers
    nullify(b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
            b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle)
    nullify(b_eps_trace_over_3_crust_mantle)
    nullify(b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
            b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core)
    nullify(b_eps_trace_over_3_inner_core)
  endif

  allocate(Iepsilondev_xx_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE), &
           Iepsilondev_yy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE), &
           Iepsilondev_xy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE), &
           Iepsilondev_xz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE), &
           Iepsilondev_yz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE), &
           Ieps_trace_over_3_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays Iepsilondev_xx_crust_mantle,..'

  allocate(R_xx_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
           R_yy_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
           R_xy_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
           R_xz_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
           R_yz_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays R_xx_crust_mantle,..'

  allocate(R_xx_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
           R_yy_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
           R_xy_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
           R_xz_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
           R_yz_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays R_xx_inner_core,..'

  if (ROTATION_VAL) then
    allocate(A_array_rotation(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION), &
             B_array_rotation(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays A_array_rotation,..'
  endif

  ! allocates backward/reconstructed arrays (dummy in case of forward simulation)
  allocate(b_displ_crust_mantle(NDIM,NGLOB_CRUST_MANTLE_ADJOINT), &
           b_veloc_crust_mantle(NDIM,NGLOB_CRUST_MANTLE_ADJOINT), &
           b_accel_crust_mantle(NDIM,NGLOB_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating b_displ,b_veloc,b_accel in crust_mantle'

  allocate(b_displ_outer_core(NGLOB_OUTER_CORE_ADJOINT), &
           b_veloc_outer_core(NGLOB_OUTER_CORE_ADJOINT), &
           b_accel_outer_core(NGLOB_OUTER_CORE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating b_displ,b_veloc,b_accel in outer_core'

  allocate(b_displ_inner_core(NDIM,NGLOB_INNER_CORE_ADJOINT), &
           b_veloc_inner_core(NDIM,NGLOB_INNER_CORE_ADJOINT), &
           b_accel_inner_core(NDIM,NGLOB_INNER_CORE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating b_displ,b_veloc,b_accel in inner_core'

  allocate(b_R_xx_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
           b_R_yy_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
           b_R_xy_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
           b_R_xz_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
           b_R_yz_crust_mantle(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays b_R_xx_crust_mantle,..'

  allocate(b_R_xx_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
           b_R_yy_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
           b_R_xy_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
           b_R_xz_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
           b_R_yz_inner_core(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays b_R_xx_inner_core,..'

  ! initializes backward/reconstructed arrays
  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL) then
      allocate(b_A_array_rotation(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT), &
               b_B_array_rotation(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error allocating arrays b_A_array_rotation,..'
    endif

  endif

  ! Runge-Kutta time scheme
  if (USE_LDDRK) then
    ! checks
    if (SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. NOISE_TOMOGRAPHY /= 0) &
        stop 'Error: LDDRK is not implemented for adjoint or noise tomography'

    ! number of stages for scheme
    NSTAGE_TIME_SCHEME = NSTAGE   ! 6 stages

    ! scheme wavefields
    allocate(displ_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array displ_crust_mantle_lddrk'
    allocate(veloc_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array veloc_crust_mantle_lddrk'
    allocate(displ_outer_core_lddrk(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array displ_outer_core_lddrk'
    allocate(veloc_outer_core_lddrk(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array veloc_outer_core_lddrk'
    allocate(displ_inner_core_lddrk(NDIM,NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array displ_inner_core_lddrk'
    allocate(veloc_inner_core_lddrk(NDIM,NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array veloc_inner_core_lddrk'

    if (SIMULATION_TYPE == 3) then
      ! scheme adjoint wavefields
      allocate(b_displ_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_displ_crust_mantle_lddrk'
      allocate(b_veloc_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_veloc_crust_mantle_lddrk'
      allocate(b_displ_outer_core_lddrk(NGLOB_OUTER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_displ_outer_core_lddrk'
      allocate(b_veloc_outer_core_lddrk(NGLOB_OUTER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_veloc_outer_core_lddrk'
      allocate(b_displ_inner_core_lddrk(NDIM,NGLOB_INNER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_displ_inner_core_lddrk'
      allocate(b_veloc_inner_core_lddrk(NDIM,NGLOB_INNER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_veloc_inner_core_lddrk'
    endif

    ! rotation in fluid outer core
    allocate(A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array A_array_rotation_lddrk'
    allocate(B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array B_array_rotation_lddrk'

    if (SIMULATION_TYPE == 3) then
      allocate(b_A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_A_array_rotation_lddrk'
      allocate(b_B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_B_array_rotation_lddrk'
    endif

    ! attenuation memory variables
    ! crust/mantle
    allocate(R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_crust_mantle_lddrk'
    ! inner core
    allocate(R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_inner_core_lddrk'

    if (SIMULATION_TYPE == 3) then
      ! crust/mantle
      allocate(b_R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_crust_mantle_lddrk'
      ! inner core
      allocate(b_R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_inner_core_lddrk'
    endif

  else
    ! default Newmark time scheme

    ! only 1 stage for Newmark time scheme
    NSTAGE_TIME_SCHEME = 1

    ! dummy arrays needed for passing as function arguments
    allocate(A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array A_array_rotation_lddrk'
    allocate(B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array B_array_rotation_lddrk'
    allocate(R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_crust_mantle_lddrk'
    allocate(R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_inner_core_lddrk'
    if (SIMULATION_TYPE == 3) then
      allocate(b_A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_A_array_rotation_lddrk'
      allocate(b_B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_B_array_rotation_lddrk'
      allocate(b_R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_R_memory_crust_mantle_lddrk'
      allocate(b_R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_R_memory_inner_core_lddrk'
    endif
  endif

  ! movies
  allocate(div_displ_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE),stat=ier)
  if (ier /= 0) stop 'Error allocating array div_displ_outer_core'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  initializing wavefields"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! initializes wavefield values
  call init_wavefields()

  ! sensitivity kernels
  ! allocates arrays
  allocate(rho_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           beta_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
           alpha_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rho_kl,.. in crust_mantle'

  allocate(rho_kl_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
           beta_kl_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
           alpha_kl_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rho_kl,.. in inner_core'

  allocate(rho_kl_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT), &
           alpha_kl_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rho_kl,.. in outer_core'

  allocate(vector_accel_outer_core(NDIM,NGLOB_OUTER_CORE_ADJOINT), &
           vector_displ_outer_core(NDIM,NGLOB_OUTER_CORE_ADJOINT), &
           b_vector_displ_outer_core(NDIM,NGLOB_OUTER_CORE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays vector_accel_outer_core,..'


  if (SIMULATION_TYPE == 3) then
    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      allocate( sigma_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise sigma kernel')
    endif

    ! approximate Hessian
    if (APPROXIMATE_HESS_KL) then
      allocate( hess_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating Hessian')
    endif

    ! For anisotropic kernels (in crust_mantle only)
    if (ANISOTROPIC_KL) then
      allocate( cijkl_kl_crust_mantle(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating full cijkl kernel in crust_mantle')
    endif

    ! deviatoric kernel check
    if (deviatoric_outercore) then
      nspec_beta_kl_outer_core = NSPEC_OUTER_CORE_ADJOINT
    else
      nspec_beta_kl_outer_core = 1
    endif
    allocate(beta_kl_outer_core(NGLLX,NGLLY,NGLLZ,nspec_beta_kl_outer_core),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating beta outercore')
  endif

  ! initializes kernel values
  if (SIMULATION_TYPE == 3) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  initializing sensitivity kernels"
      call flush_IMAIN()
    endif
    call synchronize_all()

    ! initializes kernel values
    call init_kernels()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_wavefields

!
!-------------------------------------------------------------------------------------------------
!

  subroutine init_wavefields()

! initializes arrays

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  use specfem_par_movie
  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: init_value
  ! openmp
#ifdef WAVEFIELD_INIT_WITH_OMP_PER_REGION
  integer :: iglob,iphase,ispec_p,ispec,num_elements
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif
#elif WAVEFIELD_INIT_WITH_OMP
  integer :: i
#endif

  ! put negligible initial value to avoid very slow underflow trapping
  if (FIX_UNDERFLOW_PROBLEM) then
    init_value = VERYSMALLVAL
  else
    init_value = 0._CUSTOM_REAL
  endif

#ifdef WAVEFIELD_INIT_WITH_OMP_PER_REGION
! note: after allocation, arrays have not been mapped to memory yet. this will be done with the first initialization here.
!       we initialize arrays the same way as we access them with OpenMP threads in compute_forces***() routines.
!       this ensures that memory blocks close to the thread location (on the corresponding CPU-core) will be mapped,
!       which should speedup (at least the OpenMP-) code.

  ! crust/mantle
  do iphase = 1,2
    if (iphase == 1) then
      ! outer elements (halo region)
      num_elements = nspec_outer_crust_mantle
    else
      ! inner elements
      num_elements = nspec_inner_crust_mantle
    endif
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(init_value,num_elements,iphase,phase_ispec_inner_crust_mantle,ibool_crust_mantle, &
!$OMP displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle) &
!$OMP PRIVATE( ispec,ispec_p, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP iglob)
!$OMP DO SCHEDULE(GUIDED)
    do ispec_p = 1,num_elements
      ! only compute elements which belong to current phase (inner or outer elements)
      ispec = phase_ispec_inner_crust_mantle(ispec_p,iphase)
      DO_LOOP_IJK
        iglob = ibool_crust_mantle(INDEX_IJK,ispec)
        ! initialize arrays to zero
        displ_crust_mantle(:,iglob) = init_value
        veloc_crust_mantle(:,iglob) = 0._CUSTOM_REAL
        accel_crust_mantle(:,iglob) = 0._CUSTOM_REAL
      ENDDO_LOOP_IJK
      ! if doing benchmark runs to measure scaling of the code,
      ! set the initial field to 1 to make sure gradual underflow trapping does not slow down the code
      if (DO_BENCHMARK_RUN_ONLY .and. SET_INITIAL_FIELD_TO_1_IN_BENCH) then
        DO_LOOP_IJK
          displ_crust_mantle(:,iglob) = 1._CUSTOM_REAL
          veloc_crust_mantle(:,iglob) = 1._CUSTOM_REAL
          accel_crust_mantle(:,iglob) = 1._CUSTOM_REAL
        ENDDO_LOOP_IJK
      endif
    enddo
!$OMP ENDDO
!$OMP END PARALLEL
  enddo !iphase

  ! outer core
  do iphase = 1,2
    if (iphase == 1) then
      ! outer elements (halo region)
      num_elements = nspec_outer_outer_core
    else
      ! inner elements
      num_elements = nspec_inner_outer_core
    endif
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(init_value,num_elements,iphase,phase_ispec_inner_outer_core,ibool_outer_core, &
!$OMP displ_outer_core,veloc_outer_core,accel_outer_core) &
!$OMP PRIVATE( ispec,ispec_p, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP iglob)
!$OMP DO SCHEDULE(GUIDED)
    do ispec_p = 1,num_elements
      ! only compute elements which belong to current phase (inner or outer elements)
      ispec = phase_ispec_inner_outer_core(ispec_p,iphase)
      DO_LOOP_IJK
        iglob = ibool_outer_core(INDEX_IJK,ispec)
        ! initialize arrays to zero
        displ_outer_core(iglob) = init_value
        veloc_outer_core(iglob) = 0._CUSTOM_REAL
        accel_outer_core(iglob) = 0._CUSTOM_REAL
      ENDDO_LOOP_IJK
      ! if doing benchmark runs to measure scaling of the code,
      ! set the initial field to 1 to make sure gradual underflow trapping does not slow down the code
      if (DO_BENCHMARK_RUN_ONLY .and. SET_INITIAL_FIELD_TO_1_IN_BENCH) then
        DO_LOOP_IJK
          displ_outer_core(iglob) = 1._CUSTOM_REAL
          veloc_outer_core(iglob) = 1._CUSTOM_REAL
          accel_outer_core(iglob) = 1._CUSTOM_REAL
        ENDDO_LOOP_IJK
      endif
    enddo
!$OMP ENDDO
!$OMP END PARALLEL
  enddo !iphase

  ! inner core
  do iphase = 1,2
    if (iphase == 1) then
      ! outer elements (halo region)
      num_elements = nspec_outer_inner_core
    else
      ! inner elements
      num_elements = nspec_inner_inner_core
    endif
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(init_value,num_elements,iphase,phase_ispec_inner_inner_core,ibool_inner_core, &
!$OMP displ_inner_core,veloc_inner_core,accel_inner_core) &
!$OMP PRIVATE( ispec,ispec_p, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
#else
!$OMP i,j,k, &
#endif
!$OMP iglob)
!$OMP DO SCHEDULE(GUIDED)
    do ispec_p = 1,num_elements
      ! only compute elements which belong to current phase (inner or outer elements)
      ispec = phase_ispec_inner_inner_core(ispec_p,iphase)
      DO_LOOP_IJK
        iglob = ibool_inner_core(INDEX_IJK,ispec)
        ! initialize arrays to zero
        displ_inner_core(:,iglob) = init_value
        veloc_inner_core(:,iglob) = 0._CUSTOM_REAL
        accel_inner_core(:,iglob) = 0._CUSTOM_REAL
      ENDDO_LOOP_IJK
      ! if doing benchmark runs to measure scaling of the code,
      ! set the initial field to 1 to make sure gradual underflow trapping does not slow down the code
      if (DO_BENCHMARK_RUN_ONLY .and. SET_INITIAL_FIELD_TO_1_IN_BENCH) then
        DO_LOOP_IJK
          displ_inner_core(:,iglob) = 1._CUSTOM_REAL
          veloc_inner_core(:,iglob) = 1._CUSTOM_REAL
          accel_inner_core(:,iglob) = 1._CUSTOM_REAL
        ENDDO_LOOP_IJK
      endif
    enddo
!$OMP ENDDO
!$OMP END PARALLEL
  enddo !iphase

#elif WAVEFIELD_INIT_WITH_OMP
  ! uses OpenMP to initialize wavefield arrays
  ! the initialization follows the Newmark update routines, where the looping is directly over global points

  if (FORCE_VECTORIZATION_VAL) then
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(init_value, &
!$OMP displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
!$OMP displ_outer_core,veloc_outer_core,accel_outer_core, &
!$OMP displ_inner_core,veloc_inner_core,accel_inner_core) &
!$OMP PRIVATE(i)
!$OMP DO
    do i = 1,NGLOB_CRUST_MANTLE * NDIM
      displ_crust_mantle(i,1) = init_value
      veloc_crust_mantle(i,1) = 0._CUSTOM_REAL
      accel_crust_mantle(i,1) = 0._CUSTOM_REAL
    enddo
!$OMP ENDDO
!$OMP DO
    do i = 1,NGLOB_OUTER_CORE
      displ_outer_core(i) = init_value
      veloc_outer_core(i) = 0._CUSTOM_REAL
      accel_outer_core(i) = 0._CUSTOM_REAL
    enddo
!$OMP ENDDO
!$OMP DO
    do i = 1,NGLOB_INNER_CORE * NDIM
      displ_inner_core(i,1) = init_value
      veloc_inner_core(i,1) = 0._CUSTOM_REAL
      accel_inner_core(i,1) = 0._CUSTOM_REAL
    enddo
!$OMP ENDDO
!$OMP END PARALLEL
  else
    do i = 1,NGLOB_CRUST_MANTLE
      displ_crust_mantle(:,i) = init_value
      veloc_crust_mantle(:,i) = 0._CUSTOM_REAL
      accel_crust_mantle(:,i) = 0._CUSTOM_REAL
    enddo
    do i = 1,NGLOB_OUTER_CORE
      displ_outer_core(i) = init_value
      veloc_outer_core(i) = 0._CUSTOM_REAL
      accel_outer_core(i) = 0._CUSTOM_REAL
    enddo
    do i = 1,NGLOB_INNER_CORE
      displ_inner_core(:,i) = init_value
      veloc_inner_core(:,i) = 0._CUSTOM_REAL
      accel_inner_core(:,i) = 0._CUSTOM_REAL
    enddo
  endif

#else
  ! initialization without OpenMP, most direct way
  !
  ! note: when using OpenMP for simulations, this initialization will be executed by the master thread (for each MPI) process
  !       and thus memory might be mapped closer to master thread than later on the additional OpenMP threads

  ! initialize arrays to zero
  displ_crust_mantle(:,:) = init_value
  veloc_crust_mantle(:,:) = 0._CUSTOM_REAL
  accel_crust_mantle(:,:) = 0._CUSTOM_REAL

  displ_outer_core(:) = init_value
  veloc_outer_core(:) = 0._CUSTOM_REAL
  accel_outer_core(:) = 0._CUSTOM_REAL

  displ_inner_core(:,:) = init_value
  veloc_inner_core(:,:) = 0._CUSTOM_REAL
  accel_inner_core(:,:) = 0._CUSTOM_REAL
#endif

  ! if doing benchmark runs to measure scaling of the code,
  ! set the initial field to 1 to make sure gradual underflow trapping does not slow down the code
  if (DO_BENCHMARK_RUN_ONLY .and. SET_INITIAL_FIELD_TO_1_IN_BENCH) then
    displ_crust_mantle(:,:) = 1._CUSTOM_REAL
    veloc_crust_mantle(:,:) = 1._CUSTOM_REAL
    accel_crust_mantle(:,:) = 1._CUSTOM_REAL

    displ_outer_core(:) = 1._CUSTOM_REAL
    veloc_outer_core(:) = 1._CUSTOM_REAL
    accel_outer_core(:) = 1._CUSTOM_REAL

    displ_inner_core(:,:) = 1._CUSTOM_REAL
    veloc_inner_core(:,:) = 1._CUSTOM_REAL
    accel_inner_core(:,:) = 1._CUSTOM_REAL
  endif

  ! initialize to be on the save side for adjoint runs SIMULATION_TYPE == 2
  ! crust/mantle
  eps_trace_over_3_crust_mantle(:,:,:,:) = init_value
  epsilondev_xx_crust_mantle(:,:,:,:) = init_value
  epsilondev_yy_crust_mantle(:,:,:,:) = init_value
  epsilondev_xy_crust_mantle(:,:,:,:) = init_value
  epsilondev_xz_crust_mantle(:,:,:,:) = init_value
  epsilondev_yz_crust_mantle(:,:,:,:) = init_value

  ! inner core
  eps_trace_over_3_inner_core(:,:,:,:) = init_value
  epsilondev_xx_inner_core(:,:,:,:) = init_value
  epsilondev_yy_inner_core(:,:,:,:) = init_value
  epsilondev_xy_inner_core(:,:,:,:) = init_value
  epsilondev_xz_inner_core(:,:,:,:) = init_value
  epsilondev_yz_inner_core(:,:,:,:) = init_value

  ! backward/reconstructed strain fields
  if (COMPUTE_AND_STORE_STRAIN) then
    if (MOVIE_VOLUME .and. (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3)) then
      Iepsilondev_xx_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_yy_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_xy_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_xz_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_yz_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Ieps_trace_over_3_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif

  ! clear memory variables if attenuation
  if (ATTENUATION_VAL) then
    R_xx_crust_mantle(:,:,:,:,:) = init_value
    R_yy_crust_mantle(:,:,:,:,:) = init_value
    R_xy_crust_mantle(:,:,:,:,:) = init_value
    R_xz_crust_mantle(:,:,:,:,:) = init_value
    R_yz_crust_mantle(:,:,:,:,:) = init_value

    R_xx_inner_core(:,:,:,:,:) = init_value
    R_yy_inner_core(:,:,:,:,:) = init_value
    R_xy_inner_core(:,:,:,:,:) = init_value
    R_xz_inner_core(:,:,:,:,:) = init_value
    R_yz_inner_core(:,:,:,:,:) = init_value
  endif

  if (ROTATION_VAL) then
    A_array_rotation(:,:,:,:) = 0._CUSTOM_REAL
    B_array_rotation(:,:,:,:) = 0._CUSTOM_REAL
  endif

  ! initializes backward/reconstructed arrays
  if (SIMULATION_TYPE == 3) then
    ! initializes wavefields
    b_displ_crust_mantle = 0._CUSTOM_REAL
    b_veloc_crust_mantle = 0._CUSTOM_REAL
    b_accel_crust_mantle = 0._CUSTOM_REAL

    b_displ_inner_core = 0._CUSTOM_REAL
    b_veloc_inner_core = 0._CUSTOM_REAL
    b_accel_inner_core = 0._CUSTOM_REAL

    b_displ_outer_core = 0._CUSTOM_REAL
    b_veloc_outer_core = 0._CUSTOM_REAL
    b_accel_outer_core = 0._CUSTOM_REAL

    b_epsilondev_xx_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_yy_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_xy_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_xz_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_yz_crust_mantle = 0._CUSTOM_REAL
    b_eps_trace_over_3_crust_mantle = init_value

    b_epsilondev_xx_inner_core = 0._CUSTOM_REAL
    b_epsilondev_yy_inner_core = 0._CUSTOM_REAL
    b_epsilondev_xy_inner_core = 0._CUSTOM_REAL
    b_epsilondev_xz_inner_core = 0._CUSTOM_REAL
    b_epsilondev_yz_inner_core = 0._CUSTOM_REAL
    b_eps_trace_over_3_inner_core = init_value

    if (ROTATION_VAL) then
      b_A_array_rotation(:,:,:,:) = 0._CUSTOM_REAL
      b_B_array_rotation(:,:,:,:) = 0._CUSTOM_REAL
    endif

    if (ATTENUATION_VAL) then
      b_R_xx_crust_mantle = 0._CUSTOM_REAL
      b_R_yy_crust_mantle = 0._CUSTOM_REAL
      b_R_xy_crust_mantle = 0._CUSTOM_REAL
      b_R_xz_crust_mantle = 0._CUSTOM_REAL
      b_R_yz_crust_mantle = 0._CUSTOM_REAL

      b_R_xx_inner_core = 0._CUSTOM_REAL
      b_R_yy_inner_core = 0._CUSTOM_REAL
      b_R_xy_inner_core = 0._CUSTOM_REAL
      b_R_xz_inner_core = 0._CUSTOM_REAL
      b_R_yz_inner_core = 0._CUSTOM_REAL
    endif
  endif

  ! Runge-Kutta time scheme
  if (USE_LDDRK) then
    displ_crust_mantle_lddrk(:,:) = init_value
    veloc_crust_mantle_lddrk(:,:) = 0._CUSTOM_REAL

    displ_outer_core_lddrk(:) = init_value
    veloc_outer_core_lddrk(:) = 0._CUSTOM_REAL

    displ_inner_core_lddrk(:,:) = init_value
    veloc_inner_core_lddrk(:,:) = 0._CUSTOM_REAL

    if (SIMULATION_TYPE == 3) then
      b_displ_crust_mantle_lddrk(:,:) = init_value
      b_veloc_crust_mantle_lddrk(:,:) = 0._CUSTOM_REAL
      b_displ_outer_core_lddrk(:) = init_value
      b_veloc_outer_core_lddrk(:) = 0._CUSTOM_REAL
      b_displ_inner_core_lddrk(:,:) = init_value
      b_veloc_inner_core_lddrk(:,:) = 0._CUSTOM_REAL
    endif

    ! rotation in fluid outer core
    if (ROTATION_VAL) then
      A_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
      B_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
      if (SIMULATION_TYPE == 3) then
        b_A_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
        b_B_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
      endif
    endif

    ! attenuation memory variables
    if (ATTENUATION_VAL) then
      R_xx_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_yy_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_xy_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_xz_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_yz_crust_mantle_lddrk(:,:,:,:,:) = init_value

      R_xx_inner_core_lddrk(:,:,:,:,:) = init_value
      R_yy_inner_core_lddrk(:,:,:,:,:) = init_value
      R_xy_inner_core_lddrk(:,:,:,:,:) = init_value
      R_xz_inner_core_lddrk(:,:,:,:,:) = init_value
      R_yz_inner_core_lddrk(:,:,:,:,:) = init_value
      if (SIMULATION_TYPE == 3) then
        b_R_xx_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_yy_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_xy_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_xz_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_yz_crust_mantle_lddrk(:,:,:,:,:) = init_value

        b_R_xx_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_yy_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_xy_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_xz_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_yz_inner_core_lddrk(:,:,:,:,:) = init_value
      endif
    endif
  endif

  ! movies
  div_displ_outer_core(:,:,:,:) = 0._CUSTOM_REAL

  end subroutine init_wavefields


!
!-------------------------------------------------------------------------------------------------
!

  subroutine init_kernels()

! initializes kernel arrays

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  implicit none

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  ! sensitivity kernels
  rho_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  beta_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  alpha_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL

  ! noise strength kernel
  if (NOISE_TOMOGRAPHY == 3) then
    sigma_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  endif

  ! approximate Hessian
  if (APPROXIMATE_HESS_KL) then
    hess_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  endif

  ! For anisotropic kernels (in crust_mantle only)
  if (ANISOTROPIC_KL) then
    cijkl_kl_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
  endif

  rho_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL
  alpha_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL

  rho_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL
  beta_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL
  alpha_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL

  ! deviatoric kernel check
  beta_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL

  end subroutine init_kernels

