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

  subroutine update_displ_lddrk()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! mantle
  accel_crust_mantle(:,:) = 0._CUSTOM_REAL
  ! outer core
  accel_outer_core(:) = 0._CUSTOM_REAL
  ! inner core
  accel_inner_core(:,:) = 0._CUSTOM_REAL

  end subroutine update_displ_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_displ_lddrk_backward()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! checks
  if (SIMULATION_TYPE /= 3 ) return

  ! mantle
  b_accel_crust_mantle(:,:) = 0._CUSTOM_REAL
  ! outer core
  b_accel_outer_core(:) = 0._CUSTOM_REAL
  ! inner core
  b_accel_inner_core(:,:) = 0._CUSTOM_REAL

  end subroutine update_displ_lddrk_backward


!
!-------------------------------------------------------------------------------------------------
!
! acoustic domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_acoustic_lddrk()

! updates acceleration, velocity and displacement in acoustic region (outer core)

  use specfem_par
  use specfem_par_outercore

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_LDDRK(istage)
  beta = BETA_LDDRK(istage)

  ! forward wavefields
  call update_acoustic_lddrk(NGLOB_OUTER_CORE,displ_outer_core,veloc_outer_core,accel_outer_core, &
                             displ_outer_core_lddrk,veloc_outer_core_lddrk, &
                             deltat,alpha,beta)

  end subroutine update_veloc_acoustic_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_acoustic_lddrk_backward()

! updates acceleration, velocity and displacement in acoustic region (outer core)

  use specfem_par
  use specfem_par_outercore

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_LDDRK(istage)
  beta = BETA_LDDRK(istage)

  ! backward/reconstructed wavefields
  call update_acoustic_lddrk(NGLOB_OUTER_CORE_ADJOINT,b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                             b_displ_outer_core_lddrk,b_veloc_outer_core_lddrk, &
                             b_deltat,alpha,beta)

  end subroutine update_veloc_acoustic_lddrk_backward


!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_acoustic_lddrk(NGLOB,displ,veloc,accel,displ_lddrk,veloc_lddrk,deltat,alpha,beta)

! updates acceleration and velocity in outer core

  use constants,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: NGLOB
  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(inout) :: displ_lddrk,veloc_lddrk

  real(kind=CUSTOM_REAL),intent(in) :: deltat

  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta

  ! local parameters
  integer :: i

  ! Runge-Kutta scheme update

  ! note: splitting the do-loops seems to be slightly more effective

  ! low-memory Runge-Kutta: intermediate storage wavefields
  do i = 1,NGLOB
    veloc_lddrk(i) =  alpha * veloc_lddrk(i) + deltat * accel(i)
    displ_lddrk(i) =  alpha * displ_lddrk(i) + deltat * veloc(i)
  enddo
  ! updates wavefields
  do i = 1,NGLOB
    veloc(i) = veloc(i) + beta * veloc_lddrk(i)
    displ(i) = displ(i) + beta * displ_lddrk(i)
  enddo

  end subroutine update_acoustic_lddrk


!
!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_lddrk()

! updates acceleration,velocity and displacement in elastic regions (crust/mantle,inner core)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_LDDRK(istage)
  beta = BETA_LDDRK(istage)

  ! forward wavefields
  ! crust/mantle
  call update_elastic_lddrk(NGLOB_CRUST_MANTLE,displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                            displ_crust_mantle_lddrk,veloc_crust_mantle_lddrk, &
                            deltat,alpha,beta)

  ! inner core
  call update_elastic_lddrk(NGLOB_INNER_CORE,displ_inner_core,veloc_inner_core,accel_inner_core, &
                            displ_inner_core_lddrk,veloc_inner_core_lddrk, &
                            deltat,alpha,beta)


  end subroutine update_veloc_elastic_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic_lddrk_backward()

! updates acceleration,velocity and displacement in elastic regions (crust/mantle,inner core)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none
  ! local parameters
  real(kind=CUSTOM_REAL) :: alpha,beta

  ! current Runge-Kutta coefficients
  alpha = ALPHA_LDDRK(istage)
  beta = BETA_LDDRK(istage)

  ! backward/reconstructed wavefields
  ! crust/mantle
  call update_elastic_lddrk(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                            b_displ_crust_mantle_lddrk,b_veloc_crust_mantle_lddrk, &
                            b_deltat,alpha,beta)

  ! inner core
  call update_elastic_lddrk(NGLOB_INNER_CORE_ADJOINT,b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                            b_displ_inner_core_lddrk,b_veloc_inner_core_lddrk, &
                            b_deltat,alpha,beta)

  end subroutine update_veloc_elastic_lddrk_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_elastic_lddrk(NGLOB,displ,veloc,accel, &
                                  displ_lddrk,veloc_lddrk, &
                                  deltat,alpha,beta)


  use constants_solver,only: CUSTOM_REAL,NDIM,FORCE_VECTORIZATION_VAL

  implicit none

  integer,intent(in) :: NGLOB

  ! wavefields
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(inout) :: displ_lddrk,veloc_lddrk

  real(kind=CUSTOM_REAL),intent(in) :: deltat
  ! Runge-Kutta coefficients
  real(kind=CUSTOM_REAL),intent(in) :: alpha,beta

  ! local parameters
  integer :: i

  ! low-memory Runge-Kutta scheme

  if (FORCE_VECTORIZATION_VAL) then

    ! note: splitting the do-loops seems to be slightly more effective

    ! low-memory Runge-Kutta: intermediate storage wavefields
    do i = 1,NGLOB * NDIM
      veloc_lddrk(i,1) = alpha * veloc_lddrk(i,1) + deltat * accel(i,1)
      displ_lddrk(i,1) = alpha * displ_lddrk(i,1) + deltat * veloc(i,1)
    enddo
    ! updates wavefields
    do i = 1,NGLOB * NDIM
      veloc(i,1) = veloc(i,1) + beta * veloc_lddrk(i,1)
      displ(i,1) = displ(i,1) + beta * displ_lddrk(i,1)
    enddo

  else

    ! non-vectorized loops

    do i = 1,NGLOB
      ! low-memory Runge-Kutta: intermediate storage wavefields
      veloc_lddrk(:,i) = alpha * veloc_lddrk(:,i) + deltat * accel(:,i)
      displ_lddrk(:,i) = alpha * displ_lddrk(:,i) + deltat * veloc(:,i)
      ! updates wavefields
      veloc(:,i) = veloc(:,i) + beta * veloc_lddrk(:,i)
      displ(:,i) = displ(:,i) + beta * displ_lddrk(:,i)
    enddo

  endif

  end subroutine update_elastic_lddrk


