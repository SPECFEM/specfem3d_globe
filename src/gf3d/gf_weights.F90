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
!---- Weights: the source side of an extraction as one vector per trace.
!----
!---- Everything an extraction does to an element block before the source
!---- time function is linear in the block, so for each force component a
!---- and time t it is a dot product over the block's 375 values
!---- u(a,p,i,j,k,t):
!----
!----   trace(a,t) = SUM_{p,i,j,k} u(a,p,i,j,k,t) * w(p,i,j,k)
!----
!---- A weight vector is w(GF_NCOMP,NGLLX,NGLLY,NGLLZ); seen as w(375) its
!---- index is m = p + 3*(i-1) + 15*(j-1) + 75*(k-1), the order of the
!---- writer's chunks. This module builds the vectors; it reads no file and
!---- knows no database, so it is a kernel.
!----
!---- Moment tensor. With eps the symmetric gradient and M symmetric,
!----
!----   SUM_pq M_pq eps_pq = SUM_pq M_pq du_p/dx_q
!----                      = SUM_{p,ijk} u(p,ijk) SUM_q M_pq dw(ijk,q)
!----
!---- so w(p,ijk) = SUM_q M_pq dw(ijk,q), dw the physical derivative weights
!---- of gf_strain_dweights. gf_moment_contract reads only the upper
!---- triangle of M; these weights do the same, mirroring it, so that a
!---- tensor that is symmetric only to rounding (dM/dtheta, dM/dphi) is
!---- read the same way by both.
!----
!---- Force. gf_interp_trace interpolates every displacement component with
!---- hxi(i) heta(j) hgam(k), and the force direction fhat contracts p:
!---- w(p,i,j,k) = fhat(p) hxi(i) heta(j) hgam(k).
!----
!---- The partials are weights too: one per spherical unit moment tensor,
!---- and one per centroid coordinate (lat, lon, depth), the latter the
!---- weight form of the strain-gradient partials of Stage 8 (the strain
!---- gradient, the rotation's derivative and the geographic map's, see
!---- gf_partials). None carries an amplitude scale: the
!---- caller applies the station's and, for the moment-tensor columns,
!---- 1/scale_moment, as gf_seis does.
!----
!---- Summation order differs from the strain route, so a trace contracted
!---- here agrees with gf_strain_trace + gf_moment_contract to rounding,
!---- not bit for bit. The traces are differences of terms far larger than
!---- themselves (the derivative weights sum to zero over the element), so
!---- "rounding" is eps times SUM |u w|, not eps times the trace.
!----
!---- gf_weights_contract is the extraction's contraction: every weight
!---- vector of an extraction against one station's block, sample by
!---- sample, with the force component kept. It reads the block in its
!---- stored order, the 1125 values of one time sample as one run, against
!---- each weight vector repeated three times (gf_weights_replicate), so
!---- that value k of the run meets weight (k-1)/3 + 1 and belongs to force
!---- component mod(k-1,3) + 1. Accumulating in GF_NLANE lanes, a multiple
!---- of three, keeps each lane on one component, so the loop vectorises
!---- without reassociating anything; the lanes are summed at the end.
!----

  module gf_weights

  use gf_par, only: GF_NCOMP

  implicit none

  private

  ! columns of gf_weights_mt and gf_weights_loc
  integer, parameter, public :: GF_NW_MT = 6
  integer, parameter, public :: GF_NW_LOC = 3

  ! the most weight vectors one extraction uses: the seismogram's, six
  ! moment-tensor and three position partials
  integer, parameter, public :: GF_NW_MAX = 1 + GF_NW_MT + GF_NW_LOC

  ! accumulator lanes of gf_weights_contract, a multiple of GF_NCOMP
  integer, parameter :: GF_NLANE = 48

  public :: gf_weights_moment
  public :: gf_weights_force
  public :: gf_weights_mt
  public :: gf_weights_loc
  public :: gf_weights_replicate
  public :: gf_weights_contract

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_weights_moment(dw,m,w)

! w(p,i,j,k) = SUM_q M_pq dw(i,j,k,q), with M read from its upper triangle

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, dimension(NDIM,NDIM), intent(in) :: m
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(out) :: w

  ! local parameters
  double precision, dimension(NDIM,NDIM) :: ms
  integer :: i,j,k,p,q

  ! the symmetric tensor gf_moment_contract sees
  do q = 1,NDIM
    do p = 1,NDIM
      ms(p,q) = m(min(p,q),max(p,q))
    enddo
  enddo

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        do p = 1,GF_NCOMP
          w(p,i,j,k) = ms(p,1)*dw(i,j,k,1) + ms(p,2)*dw(i,j,k,2) + ms(p,3)*dw(i,j,k,3)
        enddo
      enddo
    enddo
  enddo

  end subroutine gf_weights_moment

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_weights_force(hxi,heta,hgam,fhat,w)

! w(p,i,j,k) = fhat(p) hxi(i) heta(j) hgam(k)

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(NGLLX), intent(in) :: hxi
  double precision, dimension(NGLLY), intent(in) :: heta
  double precision, dimension(NGLLZ), intent(in) :: hgam
  double precision, dimension(NDIM), intent(in) :: fhat
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ), intent(out) :: w

  ! local parameters
  integer :: i,j,k,p

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        do p = 1,GF_NCOMP
          w(p,i,j,k) = fhat(p) * (hxi(i)*heta(j)*hgam(k))
        enddo
      enddo
    enddo
  enddo

  end subroutine gf_weights_force

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_weights_mt(dw,theta,phi,wmt)

! the weights of the six spherical unit moment tensors (Mrr, Mtt, Mpp, Mrt,
! Mrp, Mtp, get_cmt's order), each rotated about (theta,phi) exactly as the
! moment tensor itself is: the moment-tensor partials' weights, before
! their 1/scale_moment

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM
  use gf_moment, only: gf_rotate_moment_tensor

  implicit none

  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, intent(in) :: theta,phi
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ,GF_NW_MT), intent(out) :: wmt

  ! local parameters
  double precision, dimension(GF_NW_MT) :: e_sph
  double precision, dimension(NDIM,NDIM) :: m_unit
  integer :: v

  do v = 1,GF_NW_MT
    e_sph(:) = 0.d0
    e_sph(v) = 1.d0
    call gf_rotate_moment_tensor(theta,phi,e_sph,m_unit)
    call gf_weights_moment(dw,m_unit,wmt(:,:,:,:,v))
  enddo

  end subroutine gf_weights_mt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_weights_loc(dw,ddw,m_cart,dm_dtheta,dm_dphi,dtheta_dlat,dphi_dlon,jinv,dxds,wloc)

! the weights of the three centroid-position partials, lat, lon, depth:
! the per-sample chain-rule sum of Stage 8 with the strain replaced by its
! weights,
!
!   wloc(ia) = SUM_m dxds(m,ia) SUM_b jinv(b,m) W(ddw(b), M)
!            + [ia = lat] dtheta/dlat W(dw, dM/dtheta)
!            + [ia = lon] dphi/dlon   W(dw, dM/dphi)
!
! W(d, M) being gf_weights_moment and ddw(:,:,:,:,b) the derivative weights
! differentiated along reference direction b (gf_strain_ddweights). Units
! those of dxds: per degree, per degree, per km.

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM

  implicit none

  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM), intent(in) :: dw
  double precision, dimension(NGLLX,NGLLY,NGLLZ,NDIM,NDIM), intent(in) :: ddw
  double precision, dimension(NDIM,NDIM), intent(in) :: m_cart,dm_dtheta,dm_dphi,jinv
  double precision, intent(in) :: dtheta_dlat,dphi_dlon
  double precision, dimension(NDIM,GF_NW_LOC), intent(in) :: dxds
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ,GF_NW_LOC), intent(out) :: wloc

  ! local parameters
  ! g(:,b): the weights of d(M:eps)/d xi_b; gx(:,m): of d(M:eps)/dx_m
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ,NDIM) :: g,gx
  double precision, dimension(GF_NCOMP,NGLLX,NGLLY,NGLLZ) :: r
  integer :: b,mm,ia

  do b = 1,NDIM
    call gf_weights_moment(ddw(:,:,:,:,b),m_cart,g(:,:,:,:,b))
  enddo
  do mm = 1,NDIM
    gx(:,:,:,:,mm) = jinv(1,mm)*g(:,:,:,:,1) + jinv(2,mm)*g(:,:,:,:,2) + jinv(3,mm)*g(:,:,:,:,3)
  enddo

  do ia = 1,GF_NW_LOC
    wloc(:,:,:,:,ia) = gx(:,:,:,:,1)*dxds(1,ia) + gx(:,:,:,:,2)*dxds(2,ia) + gx(:,:,:,:,3)*dxds(3,ia)
  enddo

  ! the rotation of the moment tensor with the source: lat and lon only
  call gf_weights_moment(dw,dm_dtheta,r)
  wloc(:,:,:,:,1) = wloc(:,:,:,:,1) + r(:,:,:,:)*dtheta_dlat
  call gf_weights_moment(dw,dm_dphi,r)
  wloc(:,:,:,:,2) = wloc(:,:,:,:,2) + r(:,:,:,:)*dphi_dlon

  end subroutine gf_weights_loc

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_weights_replicate(w,wr)

! a weight vector laid out against one time sample of a block: wr(k) is
! w(m) for k = 3*(m-1) + a, every force component a

  use constants, only: NGLLX,NGLLY,NGLLZ

  implicit none

  double precision, dimension(GF_NCOMP*NGLLX*NGLLY*NGLLZ), intent(in) :: w
  double precision, dimension(GF_NCOMP,GF_NCOMP*NGLLX*NGLLY*NGLLZ), intent(out) :: wr

  ! local parameters
  integer :: m

  do m = 1,GF_NCOMP*NGLLX*NGLLY*NGLLZ
    wr(:,m) = w(m)
  enddo

  end subroutine gf_weights_replicate

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_weights_contract(u,nt,wr,nw,x)

! x(t,a,c) = SUM_m u(a,m,t) w_c(m), for every sample t, force component a and
! weight vector c of wr (each laid out by gf_weights_replicate)
!
! u is one station's block, displ(a,p,i,j,k,t) as stored, seen as
! u(1125,nt). The run of one sample is widened to double once and then met
! by each weight vector in turn; see the module header for the lanes.

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  integer, parameter :: NRUN = GF_NCOMP*GF_NCOMP*NGLLX*NGLLY*NGLLZ
  integer, parameter :: NFULL = (NRUN/GF_NLANE)*GF_NLANE

  integer, intent(in) :: nt,nw
  real(kind=CUSTOM_REAL), dimension(NRUN,nt), intent(in) :: u
  double precision, dimension(NRUN,nw), intent(in) :: wr
  double precision, dimension(nt,GF_NCOMP,nw), intent(out) :: x

  ! local parameters
  double precision, dimension(NRUN) :: ud
  double precision, dimension(GF_NLANE) :: acc
  integer :: it,c,k0,l,a

  do it = 1,nt
    ud(:) = dble(u(:,it))
    do c = 1,nw
      acc(:) = 0.d0
      do k0 = 0,NFULL-GF_NLANE,GF_NLANE
        do l = 1,GF_NLANE
          acc(l) = acc(l) + wr(k0+l,c)*ud(k0+l)
        enddo
      enddo
      ! the remainder starts on a multiple of GF_NLANE, hence of three, so
      ! its lanes keep their components
      do l = 1,NRUN-NFULL
        acc(l) = acc(l) + wr(NFULL+l,c)*ud(NFULL+l)
      enddo
      do a = 1,GF_NCOMP
        x(it,a,c) = sum(acc(a:GF_NLANE:GF_NCOMP))
      enddo
    enddo
  enddo

  end subroutine gf_weights_contract

  end module gf_weights
