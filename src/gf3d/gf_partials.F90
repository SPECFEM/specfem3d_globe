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
!---- Partial derivatives of a moment-tensor seismogram with respect to the
!---- source parameters, from quantities the seismogram itself already
!---- produced.
!----
!---- The parameters, their order and their units follow GF3DF's get_sdp and
!---- the GF3D Python package, so a downstream inversion changes nothing:
!----
!----    1..6   Mrr Mtt Mpp Mrt Mrp Mtp    m per dyne-cm
!----    7..9   lat lon dep                m per degree, per degree, per km
!----   10      tim                        m per second (centroid time shift)
!----
!---- `kind` 1 returns the first six, 2 all ten; there is no half-duration
!---- partial. Slots 1..9 are weight vectors (gf_weights_mt, gf_weights_loc:
!---- the strain gradient from gf_strain, the rotation's derivative from
!---- gf_moment and the geographic map's from gf_geo_chain, folded into one
!---- vector each), contracted with the element block and converted by
!---- gf_seis exactly as the seismogram is. What stays here is the slot
!---- table and the centroid-time partial, which is not a weight but a
!---- shift of the converted trace.
!----
!---- The moment-tensor partials, and why they are exact
!---- --------------------------------------------------
!---- The seismogram is linear in the moment tensor:
!----
!----   u_a(t) = STF[ scale * SUM_pq M_pq eps^a_pq(t) ]
!----
!---- with M in Cartesian coordinates, i.e. M_cart = R(theta,phi) M_sph R^T.
!---- The rotation is linear too, so the partial with respect to the v-th
!---- *spherical* component is the same block contracted with the weights
!---- of the rotated unit tensor R e_v R^T, pushed through the same
!---- conversion. No finite difference and no second element read. Per
!---- dyne-cm means the non-dimensional
!---- scale get_cmt applied is divided out again (src%scale_moment), so that
!----
!----   SUM_v M_v(CMTSOLUTION, dyne-cm) * dp(v) == seismogram
!----
!---- holds with the file's own numbers -- the identity tests/gf3d/
!---- test_gf_partials.f90 asserts at 1e-12 and test_gf_partials_db.f90
!---- repeats on a real database.
!----
!---- The centroid-time partial
!---- -------------------------
!---- The forward run evaluates its source time function at t - tshift, so
!---- the partial with respect to the shift is minus the time derivative of
!---- the seismogram. That derivative is analytic here: the seismogram is
!---- x * H_h, the pre-conversion trace convolved with the quasi-Heaviside of
!---- width h = hdur_corr, and dH_h/dt is the unit Gaussian g_h, so
!----
!----   dp(10) = -(x * g_h)
!----
!---- with g_h sampled on the stored grid and normalised to unit sum
!---- (gf_stf_kernel_gauss_unit, which says why). This replaces GF3DF's central
!---- difference of the convolved trace, whose (w dt)^2/6 error is two
!---- percent at 60 s on the 3.4 s grid -- and its `gradient`, which was off
!---- by one. In the guard case h = 0 the kernel is the identity and dp(10)
!---- is minus the pre-conversion trace: the derivative of an integral is
!---- its integrand.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module, driven
!---- by gf_seismograms with arrays read from the database and by the tests
!---- with manufactured ones.
!----

  module gf_partials

  use gf_par, only: t_gf_stf,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_ALLOC,GF_STF_HEAVI

  use gf_stf, only: gf_stf_kernel_gauss_unit,gf_conv_sym

  implicit none

  private

  ! how many partials each kind returns
  integer, parameter, public :: GF_NDP_MT  = 6
  integer, parameter, public :: GF_NDP_LOC = 10

  ! the slots, named so that no caller has to count
  integer, parameter, public :: GF_DP_MRR = 1, GF_DP_MTT = 2, GF_DP_MPP = 3
  integer, parameter, public :: GF_DP_MRT = 4, GF_DP_MRP = 5, GF_DP_MTP = 6
  integer, parameter, public :: GF_DP_LAT = 7, GF_DP_LON = 8, GF_DP_DEP = 9
  integer, parameter, public :: GF_DP_TIM = 10

  character(len=3), dimension(GF_NDP_LOC), parameter, public :: GF_DP_NAME = &
    (/ 'Mrr','Mtt','Mpp','Mrt','Mrp','Mtp','lat','lon','dep','tim' /)

  character(len=9), dimension(GF_NDP_LOC), parameter, public :: GF_DP_UNIT = &
    (/ 'm/dyne-cm','m/dyne-cm','m/dyne-cm','m/dyne-cm','m/dyne-cm','m/dyne-cm', &
       'm/deg    ','m/deg    ','m/km     ','m/s      ' /)

  public :: gf_partials_ndp
  public :: gf_partials_time

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_partials_ndp(kind,ndp,ierr)

! the number of partials a kernel type returns: 0, 6 or 10

  implicit none

  integer, intent(in) :: kind
  integer, intent(out) :: ndp
  integer, intent(out) :: ierr

  ! local parameters
  character(len=16) :: tmp

  ierr = GF_OK
  select case (kind)
  case (0)
    ndp = 0
  case (1)
    ndp = GF_NDP_MT
  case (2)
    ndp = GF_NDP_LOC
  case default
    ndp = 0
    write(tmp,'(i0)') kind
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_ndp: kind must be 0, 1 or 2, not '//trim(tmp))
  end select

  end subroutine gf_partials_ndp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_partials_time(xpad,nt,dt_sub,stf,dp10,wsum_raw,ierr)

! the centroid-time partial, from the padded pre-conversion trace
!
! `xpad(nt)` is the trace the Heaviside conversion was applied to -- the
! contracted, scaled trace on the extended axis, i.e. gf_seis_station's `xpad`
! -- and `dp10(nt)` comes back as -(xpad * g_h) with the normalised sampled
! Gaussian of the plan's width and half length. `wsum_raw` is the sum of the
! sampled Gaussian before normalisation, the aliasing measure, for the
! trace header.

  implicit none

  integer, intent(in) :: nt
  double precision, dimension(nt), intent(in) :: xpad
  double precision, intent(in) :: dt_sub
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(nt), intent(out) :: dp10
  double precision, intent(out) :: wsum_raw
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(:), allocatable :: wg,y
  integer :: it,ier

  dp10(:) = 0.d0
  wsum_raw = 0.d0

  if (stf%kind_stf /= GF_STF_HEAVI) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_time: the plan is not a Heaviside conversion')
    return
  endif
  if (dt_sub <= 0.d0) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_time: dt_sub must be positive')
    return
  endif

  allocate(wg(-stf%khalf:stf%khalf),y(nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'gf_partials_time: could not allocate the work arrays')
    return
  endif

  call gf_stf_kernel_gauss_unit(stf%hdur_corr,dt_sub,stf%khalf,wg,wsum_raw)
  call gf_conv_sym(xpad,nt,stf%khalf,wg,y)

  do it = 1,nt
    dp10(it) = -y(it)
  enddo

  deallocate(wg,y)

  ierr = GF_OK

  end subroutine gf_partials_time


  end module gf_partials
