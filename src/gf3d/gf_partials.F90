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
!---- `itypsokern` 1 returns the first six, 2 all ten; GF3DF's 3 (the
!---- half-duration partial) is dropped. Slots 1..6 and 10 are Stage 6's,
!---- 7..9 Stage 8's (gf_partials_loc, with the strain gradient from
!---- gf_strain, the rotation's derivative from gf_moment and the
!---- geographic map's from gf_geo_chain).
!----
!---- The moment-tensor partials, and why they are exact
!---- --------------------------------------------------
!---- The seismogram is linear in the moment tensor:
!----
!----   u_a(t) = STF[ scale * SUM_pq M_pq eps^a_pq(t) ]
!----
!---- with M in Cartesian coordinates, i.e. M_cart = R(theta,phi) M_sph R^T.
!---- The rotation is linear too, so the partial with respect to the v-th
!---- *spherical* component is the same strain contracted with the rotated
!---- unit tensor R e_v R^T, pushed through the same conversion. No finite
!---- difference and no second element read: the strain trace is the one
!---- the seismogram was made from. Per dyne-cm means the non-dimensional
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

  use gf_par, only: t_gf_stf,t_gf_taxis,gf_set_error, &
                    GF_OK,GF_ERR_ARG,GF_ERR_ALLOC,GF_NCOMP,GF_STF_HEAVI

  use gf_strain, only: GF_VOIGT

  use gf_moment, only: gf_rotate_moment_tensor,gf_moment_contract

  use gf_stf, only: gf_pad_left,gf_cumsum,gf_stf_apply,gf_stf_kernel_gauss_unit,gf_conv_sym

  implicit none

  private

  ! how many partials each itypsokern returns
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
  public :: gf_partials_mt
  public :: gf_partials_time
  public :: gf_partials_loc

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_partials_ndp(itypsokern,ndp,ierr)

! the number of partials a kernel type returns: 0, 6 or 10

  implicit none

  integer, intent(in) :: itypsokern
  integer, intent(out) :: ndp
  integer, intent(out) :: ierr

  ! local parameters
  character(len=16) :: tmp

  ierr = GF_OK
  select case (itypsokern)
  case (0)
    ndp = 0
  case (1)
    ndp = GF_NDP_MT
  case (2)
    ndp = GF_NDP_LOC
  case default
    ndp = 0
    write(tmp,'(i0)') itypsokern
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_ndp: itypsokern must be 0, 1 or 2, not '//trim(tmp))
  end select

  end subroutine gf_partials_ndp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_partials_mt(eps,nt_db,theta,phi,scale,tax,stf,w,dp,ierr)

! the six moment-tensor partials from a strain trace
!
! `eps(6,3,nt_db)` is the Voigt strain of the reciprocal field at the
! source on the database axis, as gf_strain_trace returns it; `theta,phi`
! the source's geocentric colatitude and longitude, about which the
! spherical unit tensors are rotated exactly as the moment tensor itself
! is; `scale` the amplitude factor per dyne-cm, i.e. the seismogram's
! 1/factor_force_source divided by the source's scale_moment; `tax`, `stf`
! and `w` the planned axis, conversion and Heaviside kernel the seismogram
! used. `dp(6,3,tax%nt)` comes back on the output axis.
!
! Everything after the contraction is the seismogram's own statement
! sequence (gf_seismograms.F90, gf_seis_cmt), so a partial computed for a
! unit tensor is bitwise the seismogram of that unit tensor.

  implicit none

  integer, intent(in) :: nt_db
  double precision, dimension(GF_VOIGT,GF_NCOMP,nt_db), intent(in) :: eps
  double precision, intent(in) :: theta,phi,scale
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(-stf%khalf:stf%khalf), intent(in) :: w
  double precision, dimension(GF_NDP_MT,GF_NCOMP,tax%nt), intent(out) :: dp
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(:), allocatable :: trace,xpad,p,y
  double precision, dimension(6) :: e_sph
  double precision, dimension(3,3) :: m_unit
  integer :: v,icomp,it,nt,ier

  dp(:,:,:) = 0.d0

  if (tax%nt_db /= nt_db .or. tax%nt < nt_db) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_mt: the time axis was not planned for this trace')
    return
  endif
  if (stf%kind_stf /= GF_STF_HEAVI) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_mt: the plan is not a Heaviside conversion')
    return
  endif

  nt = tax%nt

  allocate(trace(nt_db),xpad(nt),p(0:nt),y(nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'gf_partials_mt: could not allocate the work arrays')
    return
  endif

  do v = 1,GF_NDP_MT

    ! the v-th spherical unit tensor, rotated the way the moment tensor is
    e_sph(:) = 0.d0
    e_sph(v) = 1.d0
    call gf_rotate_moment_tensor(theta,phi,e_sph,m_unit)

    do icomp = 1,GF_NCOMP

      do it = 1,nt_db
        call gf_moment_contract(m_unit,eps(:,icomp,it),trace(it))
        trace(it) = scale * trace(it)
      enddo

      ! extend, then convert: the seismogram's own three steps
      call gf_pad_left(trace,nt_db,tax%npad,xpad)
      call gf_cumsum(xpad,nt,p)
      call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)

      do it = 1,nt
        dp(v,icomp,it) = y(it)
      enddo

    enddo

  enddo

  deallocate(trace,xpad,p,y)

  ierr = GF_OK

  end subroutine gf_partials_mt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_partials_time(xpad,nt,dt_sub,stf,dp10,wsum_raw,ierr)

! the centroid-time partial, from the padded pre-conversion trace
!
! `xpad(nt)` is the trace the Heaviside conversion was applied to -- the
! contracted, scaled strain on the extended axis, i.e. gf_seis_cmt's `xpad`
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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_partials_loc(eps,deps,nt_db,m_cart,dm_dtheta,dm_dphi,dtheta_dlat,dphi_dlon, &
                             jinv,dxds,scale,tax,stf,w,dp,ierr)

! the three centroid-position partials (Stage 8): lat, lon, depth
!
! With the seismogram u = STF[ scale SUM_pq M_pq eps_pq ] and the source
! position s = (lat, lon, depth),
!
!   du/ds_a = STF[ scale ( SUM_pq (dM_pq/ds_a) eps_pq
!                        + SUM_pq M_pq SUM_m (d eps_pq/dx_m) (dx_m/ds_a) ) ]
!
! `deps(6,3,nt_db,NDIM)` holds d eps/d xi_b -- gf_strain's kernel run with
! the differentiated weight table -- so d eps/dx_m = SUM_b jinv(b,m)
! d eps/d xi_b; `dxds(NDIM,3)` is d(x,y,z)/d(lat,lon,depth) from
! gf_geo_chain, `dm_dtheta`/`dm_dphi` the rotation's derivative from
! gf_moment, and dtheta/dlat, dphi/dlon the chain into them (the moment
! tensor does not depend on depth). `scale` is the seismogram's own
! 1/factor_force_source. Per sample the Cartesian gradient traces G_m and
! the rotation traces R_a are formed first, then combined; the conversion
! is applied once per partial, because it is linear.
!
! Units: m per degree, per degree, per km, i.e. the units of dxds.

  use constants, only: NDIM

  implicit none

  integer, intent(in) :: nt_db
  double precision, dimension(GF_VOIGT,GF_NCOMP,nt_db), intent(in) :: eps
  double precision, dimension(GF_VOIGT,GF_NCOMP,nt_db,NDIM), intent(in) :: deps
  double precision, dimension(NDIM,NDIM), intent(in) :: m_cart,dm_dtheta,dm_dphi,jinv
  double precision, intent(in) :: dtheta_dlat,dphi_dlon
  double precision, dimension(NDIM,3), intent(in) :: dxds
  double precision, intent(in) :: scale
  type(t_gf_taxis), intent(in) :: tax
  type(t_gf_stf), intent(in) :: stf
  double precision, dimension(-stf%khalf:stf%khalf), intent(in) :: w
  double precision, dimension(3,GF_NCOMP,tax%nt), intent(out) :: dp
  integer, intent(out) :: ierr

  ! local parameters
  double precision, dimension(:), allocatable :: trace,xpad,p,y
  double precision, dimension(NDIM) :: g,gx
  double precision :: r_lat,r_lon
  integer :: ia,icomp,it,b,m,nt,ier

  dp(:,:,:) = 0.d0

  if (tax%nt_db /= nt_db .or. tax%nt < nt_db) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_loc: the time axis was not planned for this trace')
    return
  endif
  if (stf%kind_stf /= GF_STF_HEAVI) then
    call gf_set_error(ierr,GF_ERR_ARG,'gf_partials_loc: the plan is not a Heaviside conversion')
    return
  endif

  nt = tax%nt

  allocate(trace(nt_db),xpad(nt),p(0:nt),y(nt),stat=ier)
  if (ier /= 0) then
    call gf_set_error(ierr,GF_ERR_ALLOC,'gf_partials_loc: could not allocate the work arrays')
    return
  endif

  do ia = 1,3
    do icomp = 1,GF_NCOMP

      do it = 1,nt_db
        ! the reference-coordinate gradient of the contracted strain,
        ! g(b) = SUM_pq M_pq d eps_pq / d xi_b, then physical, gx(m)
        do b = 1,NDIM
          call gf_moment_contract(m_cart,deps(:,icomp,it,b),g(b))
        enddo
        do m = 1,NDIM
          gx(m) = jinv(1,m)*g(1) + jinv(2,m)*g(2) + jinv(3,m)*g(3)
        enddo

        ! the position term, and the rotation term for lat and lon
        trace(it) = gx(1)*dxds(1,ia) + gx(2)*dxds(2,ia) + gx(3)*dxds(3,ia)
        if (ia == 1) then
          call gf_moment_contract(dm_dtheta,eps(:,icomp,it),r_lat)
          trace(it) = trace(it) + r_lat*dtheta_dlat
        else if (ia == 2) then
          call gf_moment_contract(dm_dphi,eps(:,icomp,it),r_lon)
          trace(it) = trace(it) + r_lon*dphi_dlon
        endif
        trace(it) = scale * trace(it)
      enddo

      call gf_pad_left(trace,nt_db,tax%npad,xpad)
      call gf_cumsum(xpad,nt,p)
      call gf_stf_apply(stf,tax%dt_sub,w,p,xpad,nt,y)

      do it = 1,nt
        dp(ia,icomp,it) = y(it)
      enddo

    enddo
  enddo

  deallocate(trace,xpad,p,y)

  ierr = GF_OK

  end subroutine gf_partials_loc

  end module gf_partials
