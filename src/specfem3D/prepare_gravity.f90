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


  subroutine prepare_gravity()

! precomputes gravity factors

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore

  implicit none

  ! local parameters
  double precision :: rspl_gravity(NR),gspl(NR),gspl2(NR)
  double precision :: radius,radius_km,g,dg
  double precision :: g_cmb_dble,g_icb_dble
  double precision :: rho,drhodr,vp,vs,Qkappa,Qmu
  double precision :: theta,phi,fac

  double precision :: minus_g,minus_dg
  double precision :: minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
  double precision :: cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision :: gxl,gyl,gzl
  double precision :: Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  integer :: int_radius,idoubling,nspl_gravity,iglob,ier

  ! minimum radius in inner core (to avoid zero radius)
  double precision, parameter :: MINIMUM_RADIUS_INNER_CORE = 100.d0 / R_EARTH

  if (myrank == 0) then
    write(IMAIN,*) "preparing gravity arrays"
    call flush_IMAIN()
  endif

  ! outer core
  ! to hold pre-computed vectors (different between gravity or no gravity case)
  allocate(gravity_pre_store_outer_core(NDIM,NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating gravity_grad_ln_density_dr array'
  gravity_pre_store_outer_core(:,:) = 0._CUSTOM_REAL

  ! store g, rho and dg/dr=dg using normalized radius in lookup table every 100 m
  ! get density and velocity from PREM model using dummy doubling flag
  ! this assumes that the gravity perturbations are small and smooth
  ! and that we can neglect the 3D model and use PREM every 100 m in all cases
  ! this is probably a rather reasonable assumption
  if (GRAVITY_VAL) then
    ! gravity term
    call make_gravity(nspl_gravity,rspl_gravity,gspl,gspl2,ONE_CRUST)

    do int_radius = 1,NRAD_GRAVITY
      radius = dble(int_radius) / (R_EARTH_KM * 10.d0)
      call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g)

      ! use PREM density profile to calculate gravity (fine for other 1D models)
      idoubling = 0
      call model_prem_iso(radius,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,.false., &
                          ONE_CRUST,.false.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                          R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

      dg = 4.0d0 * rho - 2.0d0 * g / radius

      minus_gravity_table(int_radius) = - g
      minus_deriv_gravity_table(int_radius) = - dg
      density_table(int_radius) = rho
      minus_rho_g_over_kappa_fluid(int_radius) = - g / vp**2
    enddo

    ! make sure fluid array is only assigned in outer core between 1222 and 3478 km
    ! lookup table is defined every 100 m
    do int_radius = 1,NRAD_GRAVITY
      radius_km = dble(int_radius) / 10.d0
      if (radius_km > RCMB/1000.d0 - 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RCMB/1000.d0 - 3.d0)*10.d0))
      if (radius_km < RICB/1000.d0 + 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RICB/1000.d0 + 3.d0)*10.d0))
    enddo

    ! compute gravity value at CMB and ICB once and for all
    radius = RCMB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_cmb_dble)

    radius = RICB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_icb_dble)

    ! distinguish between single and double precision for reals
    minus_g_cmb = real(- g_cmb_dble, kind=CUSTOM_REAL)
    minus_g_icb = real(- g_icb_dble, kind=CUSTOM_REAL)

    ! pre-computes gradient
    do iglob = 1,NGLOB_OUTER_CORE
      ! pre-computes factors
      ! use mesh coordinates to get theta and phi
      ! x y z contain r theta phi
      radius = dble(rstore_outer_core(1,iglob))
      theta = dble(rstore_outer_core(2,iglob))
      phi = dble(rstore_outer_core(3,iglob))

      ! integrated and multiply by rho / Kappa
      int_radius = nint(10.d0 * radius * R_EARTH_KM)
      fac = minus_rho_g_over_kappa_fluid(int_radius)

      ! gravitational acceleration (integrated and multiply by rho / Kappa) in Cartesian coordinates
      gxl = dsin(theta) * dcos(phi) * fac
      gyl = dsin(theta) * dsin(phi) * fac
      gzl = dcos(theta) * fac
      gravity_pre_store_outer_core(1,iglob) = real(gxl,kind=CUSTOM_REAL)
      gravity_pre_store_outer_core(2,iglob) = real(gyl,kind=CUSTOM_REAL)
      gravity_pre_store_outer_core(3,iglob) = real(gzl,kind=CUSTOM_REAL)
    enddo

  else
    ! no gravity
    ! tabulate d ln(rho)/dr needed for the no gravity fluid potential
    do int_radius = 1,NRAD_GRAVITY
       radius = dble(int_radius) / (R_EARTH_KM * 10.d0)
       idoubling = 0
       call model_prem_iso(radius,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,.false., &
                           ONE_CRUST,.false.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                           R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

       d_ln_density_dr_table(int_radius) = drhodr/rho
    enddo

    ! pre-computes gradient
    do iglob = 1,NGLOB_OUTER_CORE
      ! pre-computes factors
      ! use mesh coordinates to get theta and phi
      ! x y z contain r theta phi
      radius = rstore_outer_core(1,iglob)
      theta = rstore_outer_core(2,iglob)
      phi = rstore_outer_core(3,iglob)

      ! d ln(rho) / dr
      int_radius = nint(10.d0 * radius * R_EARTH_KM)
      fac = d_ln_density_dr_table(int_radius)

      ! gradient of d ln(rho)/dr in Cartesian coordinates
      gxl = dsin(theta) * dcos(phi) * fac
      gyl = dsin(theta) * dsin(phi) * fac
      gzl = dcos(theta) * fac
      gravity_pre_store_outer_core(1,iglob) = real(gxl,kind=CUSTOM_REAL)
      gravity_pre_store_outer_core(2,iglob) = real(gyl,kind=CUSTOM_REAL)
      gravity_pre_store_outer_core(3,iglob) = real(gzl,kind=CUSTOM_REAL)
    enddo
  endif

  ! crust/mantle
  if (GRAVITY_VAL) then
    ! to hold pre-computed vectors (different between gravity or no gravity case)
    allocate(gravity_pre_store_crust_mantle(NDIM,NGLOB_CRUST_MANTLE), &
             gravity_H_crust_mantle(6,NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error allocating gravity arrays for crust/mantle'
    gravity_pre_store_crust_mantle(:,:) = 0._CUSTOM_REAL
    gravity_H_crust_mantle(:,:) = 0._CUSTOM_REAL

    do iglob = 1,NGLOB_CRUST_MANTLE
      ! use mesh coordinates to get theta and phi
      ! x y and z contain r theta and phi
      radius = dble(rstore_crust_mantle(1,iglob))
      theta = dble(rstore_crust_mantle(2,iglob))
      phi = dble(rstore_crust_mantle(3,iglob))

      cos_theta = dcos(theta)
      sin_theta = dsin(theta)
      cos_phi = dcos(phi)
      sin_phi = dsin(phi)

      cos_theta_sq = cos_theta*cos_theta
      sin_theta_sq = sin_theta*sin_theta
      cos_phi_sq = cos_phi*cos_phi
      sin_phi_sq = sin_phi*sin_phi

      ! get g, rho and dg/dr=dg
      ! spherical components of the gravitational acceleration

      ! for efficiency replace with lookup table every 100 m in radial direction
      int_radius = nint(10.d0 * radius * R_EARTH_KM )

      rho = density_table(int_radius)
      minus_g = minus_gravity_table(int_radius)
      minus_dg = minus_deriv_gravity_table(int_radius)

      ! Cartesian components of the gravitational acceleration
      ! multiplied by common factor rho
      gxl = minus_g * sin_theta * cos_phi * rho
      gyl = minus_g * sin_theta * sin_phi * rho
      gzl = minus_g * cos_theta * rho
      gravity_pre_store_crust_mantle(1,iglob) = real(gxl,kind=CUSTOM_REAL)
      gravity_pre_store_crust_mantle(2,iglob) = real(gyl,kind=CUSTOM_REAL)
      gravity_pre_store_crust_mantle(3,iglob) = real(gzl,kind=CUSTOM_REAL)

      ! H-matrix
      ! Cartesian components of gradient of gravitational acceleration
      ! obtained from spherical components
      ! multiplied by common factor rho
      minus_g_over_radius = minus_g / radius
      minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius

      Hxxl = (minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq) * rho
      Hyyl = (minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq) * rho
      Hzzl = (cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq) * rho
      Hxyl = (cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq) * rho
      Hxzl = (cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta) * rho
      Hyzl = (cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta) * rho

      gravity_H_crust_mantle(1,iglob) = real(Hxxl,kind=CUSTOM_REAL)
      gravity_H_crust_mantle(2,iglob) = real(Hyyl,kind=CUSTOM_REAL)
      gravity_H_crust_mantle(3,iglob) = real(Hzzl,kind=CUSTOM_REAL)
      gravity_H_crust_mantle(4,iglob) = real(Hxyl,kind=CUSTOM_REAL)
      gravity_H_crust_mantle(5,iglob) = real(Hxzl,kind=CUSTOM_REAL)
      gravity_H_crust_mantle(6,iglob) = real(Hyzl,kind=CUSTOM_REAL)
    enddo
  else
    ! dummy allocation
    allocate(gravity_pre_store_crust_mantle(NDIM,1), &
             gravity_H_crust_mantle(6,1),stat=ier)
    if (ier /= 0) stop 'Error allocating gravity arrays for crust/mantle'
  endif


  ! inner core
  if (GRAVITY_VAL) then
    ! to hold pre-computed vectors (different between gravity or no gravity case)
    allocate(gravity_pre_store_inner_core(NDIM,NGLOB_INNER_CORE), &
             gravity_H_inner_core(6,NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0) stop 'Error allocating gravity arrays for inner core'
    gravity_pre_store_inner_core(:,:) = 0._CUSTOM_REAL
    gravity_H_inner_core(:,:) = 0._CUSTOM_REAL

    do iglob = 1,NGLOB_INNER_CORE
      ! use mesh coordinates to get theta and phi
      ! x y and z contain r theta and phi
      radius = dble(rstore_inner_core(1,iglob))
      theta = dble(rstore_inner_core(2,iglob))
      phi = dble(rstore_inner_core(3,iglob))

      ! make sure radius is never zero even for points at center of cube
      ! because we later divide by radius
      if (radius < MINIMUM_RADIUS_INNER_CORE) radius = MINIMUM_RADIUS_INNER_CORE

      cos_theta = dcos(theta)
      sin_theta = dsin(theta)
      cos_phi = dcos(phi)
      sin_phi = dsin(phi)

      cos_theta_sq = cos_theta*cos_theta
      sin_theta_sq = sin_theta*sin_theta
      cos_phi_sq = cos_phi*cos_phi
      sin_phi_sq = sin_phi*sin_phi

      ! get g, rho and dg/dr=dg
      ! spherical components of the gravitational acceleration

      ! for efficiency replace with lookup table every 100 m in radial direction
      ! make sure we never use zero for point exactly at the center of the Earth
      int_radius = max(1,nint(10.d0 * radius * R_EARTH_KM))

      rho = density_table(int_radius)
      minus_g = minus_gravity_table(int_radius)
      minus_dg = minus_deriv_gravity_table(int_radius)

      ! Cartesian components of the gravitational acceleration
      ! multiplied by common factor rho
      gxl = minus_g * sin_theta * cos_phi * rho
      gyl = minus_g * sin_theta * sin_phi * rho
      gzl = minus_g * cos_theta * rho

      gravity_pre_store_inner_core(1,iglob) = real(gxl,kind=CUSTOM_REAL)
      gravity_pre_store_inner_core(2,iglob) = real(gyl,kind=CUSTOM_REAL)
      gravity_pre_store_inner_core(3,iglob) = real(gzl,kind=CUSTOM_REAL)

      ! H-matrix
      ! Cartesian components of gradient of gravitational acceleration
      ! obtained from spherical components
      ! multiplied by common factor rho
      minus_g_over_radius = minus_g / radius
      minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius

      Hxxl = (minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq) * rho
      Hyyl = (minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq) * rho
      Hzzl = (cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq) * rho
      Hxyl = (cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq) * rho
      Hxzl = (cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta) * rho
      Hyzl = (cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta) * rho

      gravity_H_inner_core(1,iglob) = real(Hxxl,kind=CUSTOM_REAL)
      gravity_H_inner_core(2,iglob) = real(Hyyl,kind=CUSTOM_REAL)
      gravity_H_inner_core(3,iglob) = real(Hzzl,kind=CUSTOM_REAL)
      gravity_H_inner_core(4,iglob) = real(Hxyl,kind=CUSTOM_REAL)
      gravity_H_inner_core(5,iglob) = real(Hxzl,kind=CUSTOM_REAL)
      gravity_H_inner_core(6,iglob) = real(Hyzl,kind=CUSTOM_REAL)
    enddo
  else
    ! dummy allocation
    allocate(gravity_pre_store_inner_core(NDIM,1), &
             gravity_H_inner_core(6,1),stat=ier)
    if (ier /= 0) stop 'Error allocating gravity arrays for inner core'
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_gravity
