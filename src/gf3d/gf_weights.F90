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
!---- weight form of gf_partials_loc. None carries an amplitude scale: the
!---- caller applies the station's and, for the moment-tensor columns,
!---- 1/scale_moment, as gf_seis does.
!----
!---- Summation order differs from the strain route, so a trace contracted
!---- here agrees with gf_strain_trace + gf_moment_contract to rounding,
!---- not bit for bit.
!----

  module gf_weights

  use gf_par, only: GF_NCOMP

  implicit none

  private

  ! columns of gf_weights_mt and gf_weights_loc
  integer, parameter, public :: GF_NW_MT = 6
  integer, parameter, public :: GF_NW_LOC = 3

  public :: gf_weights_moment
  public :: gf_weights_force
  public :: gf_weights_mt
  public :: gf_weights_loc

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
! gf_partials_loc's per-sample sum with the strain replaced by its weights,
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

  end module gf_weights
