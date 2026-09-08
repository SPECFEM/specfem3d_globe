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

!----
!---- The derivative of the geographic map: d(x,y,z)/d(lat,lon,depth).
!----
!---- gf_geometry's gf_geographic_to_cartesian is the solver's chain
!---- (locate_sources.f90:220-354): geographic latitude to geocentric
!---- colatitude, the surface radius from the topography, ellipticity on
!---- that surface radius, and only then the depth. This module
!---- differentiates each of those steps, in the same order, with the same
!---- module variables, so that the analytic centroid partials of Stage 8
!---- and a finite difference through the forward map agree by
!---- construction and not by luck:
!----
!----   theta = pi/2 - atan(k tan(lat)),  k = ONE_MINUS_F_SQUARED   (rthetaphi_xyz.f90:214)
!----   dtheta/dlat = -k sec^2(lat) / (1 + k^2 tan^2(lat)) * pi/180
!----   dphi/dlon   = pi/180
!----   r0  = 1 + elev(lat,lon)/R                                    (gf_geometry.F90:155-158)
!----   r_s = r0 (1 - (2/3) ell(r0) P2(cos theta))                   (make_ellipticity.f90:503-516)
!----   r   = r_s - depth/R
!----   x   = r e_r(theta,phi)
!----
!---- so, with e_r, e_theta (south) and e_phi (east) the spherical basis,
!----
!----   dx/ds = dr/ds e_r + r dtheta/ds e_theta + r sin(theta) dphi/ds e_phi
!----
!---- and dx/ddepth = -(1000/R) e_r exactly, because the ellipticity is
!---- applied to the surface radius before the depth is removed. The
!---- topography enters through delev/dlat and delev/dlon, supplied by the
!---- caller (gf_database's gf_topo_gradient) for the same reason
!---- gf_geometry takes the elevation as a number: this stays a kernel
!---- module with no get_topo_bathy behind it.
!----
!---- The ellipticity's radial derivative ell'(r0) comes from the same cubic
!---- spline table the forward map evaluates (spline_routines.f90); the
!---- spline's derivative is closed form from its coefficients, below.
!----
!---- No `use hdf5`, no `use specfem_par`: a kernel module.
!----

  module gf_geo_chain

  use gf_par, only: gf_set_error,GF_OK,GF_ERR_ARG,GF_ERR_MISMATCH

  use gf_geometry, only: gf_geographic_to_cartesian

  implicit none

  private

  public :: gf_spline_derivative
  public :: gf_geographic_jacobian

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_spline_derivative(xpoint,ypoint,spline_coefficients,npoint,x,dy,ierr)

! the derivative of the cubic spline spline_evaluation() evaluates
!
! spline_evaluation (spline_routines.f90:82-137) returns, on the interval
! [x_lo, x_hi] of width h found by bisection,
!
!   y = c1 y_lo + c2 y_hi + [(c1^3 - c1) s_lo + (c2^3 - c2) s_hi] h^2/6
!
! with c1 = (x_hi - x)/h, c2 = (x - x_lo)/h and s the second-derivative
! coefficients. With dc1/dx = -1/h and dc2/dx = 1/h,
!
!   dy/dx = (y_hi - y_lo)/h + [(3 c2^2 - 1) s_hi - (3 c1^2 - 1) s_lo] h/6
!
! The bisection is repeated here rather than shared, because the original
! stops on a zero-width interval and this library returns ierr instead.

  implicit none

  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint,spline_coefficients
  double precision, intent(in) :: x
  double precision, intent(out) :: dy
  integer, intent(out) :: ierr

  ! local parameters
  integer :: index_loop,index_lower,index_higher
  double precision :: h,c1,c2

  dy = 0.d0

  if (npoint < 2) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_spline_derivative: a spline needs at least two points')
    return
  endif

  ! the interval, by dichotomy, exactly as spline_evaluation finds it
  index_lower = 1
  index_higher = npoint
  do while (index_higher - index_lower > 1)
    index_loop = (index_higher + index_lower) / 2
    if (index_loop < 1) index_loop = 1
    if (xpoint(index_loop) > x) then
      index_higher = index_loop
    else
      index_lower = index_loop
    endif
  enddo

  h = xpoint(index_higher) - xpoint(index_lower)
  if (h == 0.d0) then
    call gf_set_error(ierr,GF_ERR_MISMATCH,'gf_spline_derivative: a zero-width interval in the spline table')
    return
  endif

  c1 = (xpoint(index_higher) - x) / h
  c2 = (x - xpoint(index_lower)) / h

  dy = (ypoint(index_higher) - ypoint(index_lower)) / h &
       + ((3.d0*c2**2 - 1.d0)*spline_coefficients(index_higher) &
          - (3.d0*c1**2 - 1.d0)*spline_coefficients(index_lower)) * h / 6.d0

  ierr = GF_OK

  end subroutine gf_spline_derivative

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_geographic_jacobian(lat,lon,depth_km,is_elliptical,elevation,delev_dlat,delev_dlon, &
                                    nspl,rspl,ell1,ell2,R_PLANET,dxds,dtheta_dlat,dphi_dlon,ierr)

! d(x,y,z)/d(lat,lon,depth) at a geographic position
!
! `dxds(NDIM,3)`: columns lat (per degree), lon (per degree), depth (per
! km). `elevation` is in metres and `delev_dlat`, `delev_dlon` in metres
! per degree, as gf_topo_gradient returns them (zero without topography).
! The spline arguments are the database's, as for
! gf_geographic_to_cartesian. `dtheta_dlat` and `dphi_dlon` (radians per
! degree) are returned as well, because the moment tensor's rotation is
! differentiated with respect to theta and phi and needs the same chain.
!
! The forward map is called first, so every intermediate -- the reduced
! angles, the surface radius -- is the one the located source used.

  use constants, only: NDIM,NR_DENSITY,DEGREES_TO_RADIANS,R_UNIT_SPHERE,ASSUME_PERFECT_SPHERE

  use shared_parameters, only: ONE_MINUS_F_SQUARED

  implicit none

  double precision, intent(in) :: lat,lon,depth_km
  logical, intent(in) :: is_elliptical
  double precision, intent(in) :: elevation,delev_dlat,delev_dlon
  integer, intent(in) :: nspl
  double precision, dimension(:), intent(in) :: rspl,ell1,ell2
  double precision, intent(in) :: R_PLANET
  double precision, dimension(NDIM,3), intent(out) :: dxds
  double precision, intent(out) :: dtheta_dlat,dphi_dlon
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(NDIM) :: xyz,e_r,e_t,e_p
  double precision, dimension(NR_DENSITY) :: rspl_p,ell1_p,ell2_p
  double precision :: theta,phi,r_surface,r
  double precision :: l,k,tl
  double precision :: r0,dr0_dlat,dr0_dlon
  double precision :: ell,dell,cost,sint,p20,dp20_dtheta,factor
  double precision :: drs_dlat,drs_dlon,dr_ddep
  double precision :: sinp,cosp

  dxds(:,:) = 0.d0
  dtheta_dlat = 0.d0
  dphi_dlon = 0.d0

  ! the forward map, for the reduced angles and the surface radius
  call gf_geographic_to_cartesian(lat,lon,depth_km,is_elliptical,elevation, &
                                  nspl,rspl,ell1,ell2,R_PLANET,xyz,theta,phi,r_surface,ierr)
  if (ierr /= GF_OK) return

  !--- the angles ---------------------------------------------------------
  !
  ! lat_2_geocentric_colat_dble (rthetaphi_xyz.f90:191-229): the
  ! ellipticity branch only when the mesh is elliptical and the constant
  ! ASSUME_PERFECT_SPHERE allows it, with the flattening from the module
  ! variable the forward map used (never a literal here)
  l = lat*DEGREES_TO_RADIANS
  if (.not. ASSUME_PERFECT_SPHERE .and. is_elliptical) then
    k = ONE_MINUS_F_SQUARED
    tl = tan(l)
    dtheta_dlat = -k*(1.d0 + tl**2)/(1.d0 + (k*tl)**2) * DEGREES_TO_RADIANS
  else
    dtheta_dlat = -DEGREES_TO_RADIANS
  endif
  ! the longitude wrap to [0,360] is a constant shift
  dphi_dlon = DEGREES_TO_RADIANS

  !--- the surface radius --------------------------------------------------
  !
  ! r0 = 1 + elev/R before the ellipticity (gf_geometry.F90:155-158)
  r0 = R_UNIT_SPHERE + elevation/R_PLANET
  dr0_dlat = delev_dlat/R_PLANET
  dr0_dlon = delev_dlon/R_PLANET

  sint = sin(theta)
  cost = cos(theta)

  if (is_elliptical) then
    if (nspl <= 0 .or. nspl > NR_DENSITY) then
      call gf_set_error(ierr,GF_ERR_MISMATCH,'gf_geographic_jacobian: unusable spline table length')
      return
    endif
    rspl_p(:) = 0.d0
    ell1_p(:) = 0.d0
    ell2_p(:) = 0.d0
    rspl_p(1:nspl) = rspl(1:nspl)
    ell1_p(1:nspl) = ell1(1:nspl)
    ell2_p(1:nspl) = ell2(1:nspl)

    ! ell(r0) as add_ellipticity_rtheta evaluates it, and its derivative
    call spline_evaluation(rspl_p,ell1_p,ell2_p,nspl,r0,ell)
    call gf_spline_derivative(rspl_p,ell1_p,ell2_p,nspl,r0,dell,ierr)
    if (ierr /= GF_OK) return

    ! r_s = r0 (1 - (2/3) ell P2(cos theta)), P2 = (3 cos^2 - 1)/2
    p20 = 0.5d0*(3.d0*cost*cost - 1.d0)
    dp20_dtheta = -3.d0*cost*sint
    factor = 1.d0 - (2.d0/3.d0)*ell*p20

    drs_dlat = dr0_dlat*factor - r0*(2.d0/3.d0)*(dell*dr0_dlat*p20 + ell*dp20_dtheta*dtheta_dlat)
    drs_dlon = dr0_dlon*factor - r0*(2.d0/3.d0)*(dell*dr0_dlon*p20)
  else
    drs_dlat = dr0_dlat
    drs_dlon = dr0_dlon
  endif

  !--- the depth, subtracted last ------------------------------------------

  r = r_surface - depth_km*1000.d0/R_PLANET
  dr_ddep = -1000.d0/R_PLANET

  !--- the spherical basis at the point ------------------------------------

  sinp = sin(phi)
  cosp = cos(phi)
  e_r(1) = sint*cosp ; e_r(2) = sint*sinp ; e_r(3) = cost
  e_t(1) = cost*cosp ; e_t(2) = cost*sinp ; e_t(3) = -sint
  e_p(1) = -sinp     ; e_p(2) = cosp      ; e_p(3) = 0.d0

  dxds(:,1) = drs_dlat*e_r(:) + r*dtheta_dlat*e_t(:)
  dxds(:,2) = drs_dlon*e_r(:) + r*sint*dphi_dlon*e_p(:)
  dxds(:,3) = dr_ddep*e_r(:)

  ierr = GF_OK

  end subroutine gf_geographic_jacobian

  end module gf_geo_chain
