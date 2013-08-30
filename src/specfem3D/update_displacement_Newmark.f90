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


  subroutine update_displacement_Newmark()

! explicit Newmark time scheme with acoustic & elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 delta_t chi_dot_dot(t+delta_t)
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!
! note that this stage calculates the predictor terms
!
!   for
!   potential chi_dot(t+delta) requires + 1/2 delta_t chi_dot_dot(t+delta_t)
!                                   at a later stage (corrector) once where chi_dot_dot(t+delta) is calculated
!   and similar,
!   velocity v(t+delta_t) requires  + 1/2 delta_t a(t+delta_t)
!                                   at a later stage once where a(t+delta) is calculated
! also:
!   boundary term B_elastic requires chi_dot_dot(t+delta)
!                                   thus chi_dot_dot has to be updated first before the elastic boundary term is considered

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! updates wavefields
  if( .not. GPU_MODE) then
    ! on CPU

    ! Newmark time scheme update
    ! mantle
    call update_displ_elastic(NGLOB_CRUST_MANTLE,displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                                deltat,deltatover2,deltatsqover2)
    ! outer core
    call update_displ_acoustic(NGLOB_OUTER_CORE,displ_outer_core,veloc_outer_core,accel_outer_core, &
                                deltat,deltatover2,deltatsqover2)
    ! inner core
    call update_displ_elastic(NGLOB_INNER_CORE,displ_inner_core,veloc_inner_core,accel_inner_core, &
                                deltat,deltatover2,deltatsqover2)
  else
    ! on GPU
    ! Includes FORWARD_OR_ADJOINT == 1
    ! outer core region
    call it_update_displacement_oc_cuda(Mesh_pointer,deltat,deltatsqover2,deltatover2,1)
    ! inner core region
    call it_update_displacement_ic_cuda(Mesh_pointer,deltat,deltatsqover2,deltatover2,1)
    ! crust/mantle region
    call it_update_displacement_cm_cuda(Mesh_pointer,deltat,deltatsqover2,deltatover2,1)
  endif

  end subroutine update_displacement_Newmark

!
!-------------------------------------------------------------------------------------------------
!


  subroutine update_displacement_Newmark_backward()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! checks
  if( SIMULATION_TYPE /= 3 ) return

  ! updates wavefields
  if( .not. GPU_MODE) then
    ! on CPU
    ! Newmark time scheme update for backward/reconstructed fields
    ! mantle
    call update_displ_elastic(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                              b_deltat,b_deltatover2,b_deltatsqover2)
    ! outer core
    call update_displ_acoustic(NGLOB_OUTER_CORE_ADJOINT,b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                              b_deltat,b_deltatover2,b_deltatsqover2)
    ! inner core
    call update_displ_elastic(NGLOB_INNER_CORE_ADJOINT,b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                              b_deltat,b_deltatover2,b_deltatsqover2)
  else
    ! on GPU
    ! Includes FORWARD_OR_ADJOINT == 3
    ! outer core region
    call it_update_displacement_oc_cuda(Mesh_pointer,b_deltat,b_deltatsqover2,b_deltatover2,3)
    ! inner core region
    call it_update_displacement_ic_cuda(Mesh_pointer,b_deltat,b_deltatsqover2,b_deltatover2,3)
    ! crust/mantle region
    call it_update_displacement_cm_cuda(Mesh_pointer,b_deltat,b_deltatsqover2,b_deltatover2,3)
  endif

  end subroutine update_displacement_Newmark_backward

!
!-------------------------------------------------------------------------------------------------
!


  subroutine update_displ_elastic(NGLOB,displ,veloc,accel, &
                                  deltat,deltatover2,deltatsqover2)

  use constants_solver,only: CUSTOM_REAL,NDIM,FORCE_VECTORIZATION_VAL

  implicit none

  integer,intent(in) :: NGLOB
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL),intent(in) :: deltat,deltatover2,deltatsqover2

  ! local parameters
  integer :: i

  ! Newmark time scheme update
  if(FORCE_VECTORIZATION_VAL) then
    do i=1,NGLOB * NDIM
      displ(i,1) = displ(i,1) + deltat * veloc(i,1) + deltatsqover2 * accel(i,1)
      veloc(i,1) = veloc(i,1) + deltatover2 * accel(i,1)
      accel(i,1) = 0._CUSTOM_REAL
    enddo
  else
    do i=1,NGLOB
      displ(:,i) = displ(:,i) + deltat * veloc(:,i) + deltatsqover2 * accel(:,i)
      veloc(:,i) = veloc(:,i) + deltatover2 * accel(:,i)
      accel(:,i) = 0._CUSTOM_REAL
    enddo
  endif

  end subroutine update_displ_elastic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_displ_acoustic(NGLOB,displ,veloc,accel, &
                                   deltat,deltatover2,deltatsqover2)

  use constants,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: NGLOB
  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(inout) :: displ,veloc,accel
  real(kind=CUSTOM_REAL),intent(in) :: deltat,deltatover2,deltatsqover2

  ! local parameters
  integer :: i

  ! Newmark time scheme update
  do i=1,NGLOB
    displ(i) = displ(i) + deltat * veloc(i) + deltatsqover2 * accel(i)
    veloc(i) = veloc(i) + deltatover2 * accel(i)
    accel(i) = 0._CUSTOM_REAL
  enddo

  end subroutine update_displ_acoustic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_acoustic(NGLOB,veloc_outer_core,accel_outer_core, &
                                  deltatover2,rmass_outer_core)

! updates acceleration and velocity in outer core

  use constants_solver,only: CUSTOM_REAL

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

  do i=1,NGLOB
    accel_outer_core(i) = accel_outer_core(i)*rmass_outer_core(i)
    veloc_outer_core(i) = veloc_outer_core(i) + deltatover2*accel_outer_core(i)
  enddo

  end subroutine update_veloc_acoustic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_accel_elastic(NGLOB,NGLOB_XY,veloc_crust_mantle,accel_crust_mantle, &
                                  two_omega_earth, &
                                  rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle)

! updates acceleration in crust/mantle region

  use constants_solver,only: CUSTOM_REAL,NDIM,NCHUNKS_VAL

  use specfem_par,only: ABSORBING_CONDITIONS

  implicit none

  integer :: NGLOB,NGLOB_XY

  ! velocity & acceleration
  ! crust/mantle region
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: veloc_crust_mantle,accel_crust_mantle

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY) :: rmassx_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY) :: rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB)    :: rmassz_crust_mantle

  real(kind=CUSTOM_REAL) :: two_omega_earth

  ! local parameters
  integer :: i

  ! updates acceleration w/ rotation in crust/mantle region only

  if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then

     do i=1,NGLOB
        accel_crust_mantle(1,i) = accel_crust_mantle(1,i)*rmassx_crust_mantle(i) &
             + two_omega_earth*veloc_crust_mantle(2,i)
        accel_crust_mantle(2,i) = accel_crust_mantle(2,i)*rmassy_crust_mantle(i) &
             - two_omega_earth*veloc_crust_mantle(1,i)
        accel_crust_mantle(3,i) = accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
     enddo

  else

     do i=1,NGLOB
        accel_crust_mantle(1,i) = accel_crust_mantle(1,i)*rmassz_crust_mantle(i) &
             + two_omega_earth*veloc_crust_mantle(2,i)
        accel_crust_mantle(2,i) = accel_crust_mantle(2,i)*rmassz_crust_mantle(i) &
             - two_omega_earth*veloc_crust_mantle(1,i)
        accel_crust_mantle(3,i) = accel_crust_mantle(3,i)*rmassz_crust_mantle(i)
     enddo

  endif

  end subroutine update_accel_elastic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine update_veloc_elastic(NGLOB_CM,veloc_crust_mantle,accel_crust_mantle, &
                                  NGLOB_IC,veloc_inner_core,accel_inner_core, &
                                  deltatover2,two_omega_earth,rmass_inner_core)

! updates velocity in crust/mantle region, and acceleration and velocity in inner core

  use constants_solver,only: CUSTOM_REAL,NDIM,FORCE_VECTORIZATION_VAL

  implicit none

  integer :: NGLOB_CM,NGLOB_IC

  ! acceleration & velocity
  ! crust/mantle region
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM) :: veloc_crust_mantle,accel_crust_mantle
  ! inner core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC) :: veloc_inner_core,accel_inner_core

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_IC) :: rmass_inner_core

  real(kind=CUSTOM_REAL) :: deltatover2,two_omega_earth

  ! local parameters
  integer :: i

  ! Newmark time scheme:
  !
  ! note:
  !   - crust/mantle region
  !         needs only velocity corrector terms
  !         (acceleration already updated before)
  !   - inner core region
  !         needs both, acceleration update & velocity corrector terms

  ! mantle
  if(FORCE_VECTORIZATION_VAL) then
    do i=1,NGLOB_CM * NDIM
      veloc_crust_mantle(i,i) = veloc_crust_mantle(i,i) + deltatover2*accel_crust_mantle(i,i)
    enddo
  else
    do i=1,NGLOB_CM
      veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) + deltatover2*accel_crust_mantle(:,i)
    enddo
  endif

  ! inner core
  if(FORCE_VECTORIZATION_VAL) then
    do i=1,NGLOB_IC
      accel_inner_core(1,i) = accel_inner_core(1,i)*rmass_inner_core(i) &
             + two_omega_earth*veloc_inner_core(2,i)
      accel_inner_core(2,i) = accel_inner_core(2,i)*rmass_inner_core(i) &
             - two_omega_earth*veloc_inner_core(1,i)
      accel_inner_core(3,i) = accel_inner_core(3,i)*rmass_inner_core(i)
    enddo
    do i=1,NGLOB_IC * NDIM
      veloc_inner_core(i,1) = veloc_inner_core(i,1) + deltatover2*accel_inner_core(i,1)
    enddo
  else
    do i=1,NGLOB_IC
      accel_inner_core(1,i) = accel_inner_core(1,i)*rmass_inner_core(i) &
             + two_omega_earth*veloc_inner_core(2,i)
      accel_inner_core(2,i) = accel_inner_core(2,i)*rmass_inner_core(i) &
             - two_omega_earth*veloc_inner_core(1,i)
      accel_inner_core(3,i) = accel_inner_core(3,i)*rmass_inner_core(i)

      veloc_inner_core(:,i) = veloc_inner_core(:,i) + deltatover2*accel_inner_core(:,i)
    enddo
  endif

  end subroutine update_veloc_elastic
