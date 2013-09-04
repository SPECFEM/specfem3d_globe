!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  subroutine update_displacement_lddrk()

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

  end subroutine update_displacement_lddrk

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_displacement_lddrk_backward()

! low-memory Runge-Kutta time scheme

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! checks
  if( SIMULATION_TYPE /= 3 ) return

  ! mantle
  b_accel_crust_mantle(:,:) = 0._CUSTOM_REAL
  ! outer core
  b_accel_outer_core(:) = 0._CUSTOM_REAL
  ! inner core
  b_accel_inner_core(:,:) = 0._CUSTOM_REAL

  end subroutine update_displacement_lddrk_backward


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
  integer :: i

  do i=1,NGLOB_OUTER_CORE
    veloc_outer_core_lddrk(i) =  ALPHA_LDDRK(istage) * veloc_outer_core_lddrk(i) + deltat * accel_outer_core(i)

    displ_outer_core_lddrk(i) =  ALPHA_LDDRK(istage) * displ_outer_core_lddrk(i) + deltat * veloc_outer_core(i)

    veloc_outer_core(i) = veloc_outer_core(i) + BETA_LDDRK(istage) * veloc_outer_core_lddrk(i)

    displ_outer_core(i) = displ_outer_core(i) + BETA_LDDRK(istage) * displ_outer_core_lddrk(i)
  enddo

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
  integer :: i

  do i=1,NGLOB_OUTER_CORE
    b_veloc_outer_core_lddrk(i) =  ALPHA_LDDRK(istage) * b_veloc_outer_core_lddrk(i) + b_deltat * b_accel_outer_core(i)

    b_displ_outer_core_lddrk(i) =  ALPHA_LDDRK(istage) * b_displ_outer_core_lddrk(i) + b_deltat * b_veloc_outer_core(i)

    b_veloc_outer_core(i) = b_veloc_outer_core(i) + BETA_LDDRK(istage) * b_veloc_outer_core_lddrk(i)

    b_displ_outer_core(i) = b_displ_outer_core(i) + BETA_LDDRK(istage) * b_displ_outer_core_lddrk(i)
  enddo

  end subroutine update_veloc_acoustic_lddrk_backward

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
  integer :: i

  ! crust/mantle
  do i=1,NGLOB_CRUST_MANTLE
    veloc_crust_mantle_lddrk(:,i) = ALPHA_LDDRK(istage) * veloc_crust_mantle_lddrk(:,i) &
                                    + deltat * accel_crust_mantle(:,i)

    displ_crust_mantle_lddrk(:,i) = ALPHA_LDDRK(istage) * displ_crust_mantle_lddrk(:,i) &
                                    + deltat * veloc_crust_mantle(:,i)

    veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) + BETA_LDDRK(istage) * veloc_crust_mantle_lddrk(:,i)

    displ_crust_mantle(:,i) = displ_crust_mantle(:,i) + BETA_LDDRK(istage) * displ_crust_mantle_lddrk(:,i)
  enddo

  ! inner core
  do i=1,NGLOB_INNER_CORE
    veloc_inner_core_lddrk(:,i) = ALPHA_LDDRK(istage) * veloc_inner_core_lddrk(:,i) &
                                    + deltat * accel_inner_core(:,i)

    displ_inner_core_lddrk(:,i) = ALPHA_LDDRK(istage) * displ_inner_core_lddrk(:,i) &
                                    + deltat * veloc_inner_core(:,i)

    veloc_inner_core(:,i) = veloc_inner_core(:,i) + BETA_LDDRK(istage) * veloc_inner_core_lddrk(:,i)

    displ_inner_core(:,i) = displ_inner_core(:,i) + BETA_LDDRK(istage) * displ_inner_core_lddrk(:,i)
  enddo

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
  integer :: i

  ! crust/mantle
  do i=1,NGLOB_CRUST_MANTLE
    b_veloc_crust_mantle_lddrk(:,i) = ALPHA_LDDRK(istage) * b_veloc_crust_mantle_lddrk(:,i) &
                                    + b_deltat * b_accel_crust_mantle(:,i)

    b_displ_crust_mantle_lddrk(:,i) = ALPHA_LDDRK(istage) * b_displ_crust_mantle_lddrk(:,i) &
                                    + b_deltat * b_veloc_crust_mantle(:,i)

    b_veloc_crust_mantle(:,i) = b_veloc_crust_mantle(:,i) + BETA_LDDRK(istage) * b_veloc_crust_mantle_lddrk(:,i)

    b_displ_crust_mantle(:,i) = b_displ_crust_mantle(:,i) + BETA_LDDRK(istage) * b_displ_crust_mantle_lddrk(:,i)
  enddo

  ! inner core
  do i=1,NGLOB_INNER_CORE
    b_veloc_inner_core_lddrk(:,i) = ALPHA_LDDRK(istage) * b_veloc_inner_core_lddrk(:,i) &
                                    + b_deltat * b_accel_inner_core(:,i)

    b_displ_inner_core_lddrk(:,i) = ALPHA_LDDRK(istage) * b_displ_inner_core_lddrk(:,i) &
                                    + b_deltat * b_veloc_inner_core(:,i)

    b_veloc_inner_core(:,i) = b_veloc_inner_core(:,i) + BETA_LDDRK(istage) * b_veloc_inner_core_lddrk(:,i)

    b_displ_inner_core(:,i) = b_displ_inner_core(:,i) + BETA_LDDRK(istage) * b_displ_inner_core_lddrk(:,i)
  enddo

  end subroutine update_veloc_elastic_lddrk_backward

