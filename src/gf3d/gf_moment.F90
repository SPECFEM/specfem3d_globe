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
!---- The moment tensor: spherical to Cartesian, and its contraction with
!---- the strain.
!----
!---- CMTSOLUTION gives the moment tensor in the spherical basis
!---- (r, theta, phi) as (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp), and get_cmt.f90
!---- returns it in that order, already non-dimensionalised by
!----
!----     scaleM = 1.d7 * RHOAV * R_PLANET**5 * PI*GRAV*RHOAV
!----
!---- (get_cmt.f90:426). The rotation to Cartesian below is transcribed from
!---- locate_sources.f90:257-273, which writes out M_cart = R M_sph R^T for
!---- the (r,theta,phi) -> (x,y,z) rotation term by term.
!----
!---- Transcribed rather than re-derived on purpose. The expanded form is
!---- long and every term is a place to lose a sign, and a sign error in an
!---- off-diagonal is not a crash -- it is a seismogram with the wrong
!---- radiation pattern, which looks entirely reasonable. Keeping it
!---- character-for-character identical to the solver's makes the two
!---- diffable.
!----
!---- The contraction and its factor of two
!---- -------------------------------------
!---- The seismogram is linear in M:
!----
!----   u_a(t) = SUM_pq M_pq eps^a_pq
!----
!---- and with eps stored in Voigt form [xx,yy,zz,xy,xz,yz] -- where eps(4)
!---- is eps_xy itself, not 2*eps_xy -- the six-term sum needs the
!---- off-diagonal terms counted twice, because each stands for two equal
!---- entries of the symmetric pair:
!----
!----   u_a = Mxx exx + Myy eyy + Mzz ezz + 2(Mxy exy + Mxz exz + Myz eyz)
!----
!---- That factor is the single most common way to get a moment-tensor
!---- contraction wrong, so gf_moment_contract_full() exists beside it: it
!---- performs the plain 3x3 double sum with no Voigt packing at all, and
!---- tests/gf3d/test_gf_strain.f90 asserts the two agree. One of them is
!---- obviously right; the other is the one production uses.
!----
!---- Linearity is also what makes Stage 6 nearly free: the partials with
!---- respect to the six moment-tensor components are the same strain
!---- contracted with six unit tensors, which is why nothing here recomputes
!---- a strain.
!----
!---- No `use hdf5`, no `use specfem_par`: this is a kernel module.
!----

  module gf_moment

  use gf_par, only: GF_NCOMP

  use gf_strain, only: GF_VOIGT,GF_XX,GF_YY,GF_ZZ,GF_XY,GF_XZ,GF_YZ

  implicit none

  private

  public :: gf_rotate_moment_tensor
  public :: gf_moment_contract
  public :: gf_moment_contract_full
  public :: gf_moment_unit_tensor

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_rotate_moment_tensor(theta,phi,m_sph,m_cart)

! spherical (r,theta,phi) moment tensor to Cartesian
!
! `m_sph` is (Mrr, Mtt, Mpp, Mrt, Mrp, Mtp) -- get_cmt's moment_tensor(1:6)
! order -- and `m_cart` comes back as the full symmetric 3x3.
!
! Transcribed from locate_sources.f90:257-273.

  use constants, only: NDIM

  implicit none

  double precision, intent(in) :: theta,phi
  double precision, dimension(6), intent(in) :: m_sph
  double precision, dimension(NDIM,NDIM), intent(out) :: m_cart

  ! local parameters
  double precision :: sint,cost,sinp,cosp
  double precision :: Mrr,Mtt,Mpp,Mrt,Mrp,Mtp
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  sint = sin(theta)
  cost = cos(theta)
  sinp = sin(phi)
  cosp = cos(phi)

  Mrr = m_sph(1)
  Mtt = m_sph(2)
  Mpp = m_sph(3)
  Mrt = m_sph(4)
  Mrp = m_sph(5)
  Mtp = m_sph(6)

  ! convert from a spherical to a Cartesian representation of the moment tensor
  Mxx = sint*sint*cosp*cosp*Mrr + cost*cost*cosp*cosp*Mtt + sinp*sinp*Mpp &
      + 2.0d0*sint*cost*cosp*cosp*Mrt - 2.0d0*sint*sinp*cosp*Mrp - 2.0d0*cost*sinp*cosp*Mtp

  Myy = sint*sint*sinp*sinp*Mrr + cost*cost*sinp*sinp*Mtt + cosp*cosp*Mpp &
      + 2.0d0*sint*cost*sinp*sinp*Mrt + 2.0d0*sint*sinp*cosp*Mrp + 2.0d0*cost*sinp*cosp*Mtp

  Mzz = cost*cost*Mrr + sint*sint*Mtt - 2.0d0*sint*cost*Mrt

  Mxy = sint*sint*sinp*cosp*Mrr + cost*cost*sinp*cosp*Mtt - sinp*cosp*Mpp &
      + 2.0d0*sint*cost*sinp*cosp*Mrt + sint*(cosp*cosp-sinp*sinp)*Mrp &
      + cost*(cosp*cosp-sinp*sinp)*Mtp

  Mxz = sint*cost*cosp*Mrr - sint*cost*cosp*Mtt &
      + (cost*cost-sint*sint)*cosp*Mrt - cost*sinp*Mrp + sint*sinp*Mtp

  Myz = sint*cost*sinp*Mrr - sint*cost*sinp*Mtt &
      + (cost*cost-sint*sint)*sinp*Mrt + cost*cosp*Mrp - sint*cosp*Mtp

  m_cart(1,1) = Mxx
  m_cart(2,2) = Myy
  m_cart(3,3) = Mzz
  m_cart(1,2) = Mxy
  m_cart(2,1) = Mxy
  m_cart(1,3) = Mxz
  m_cart(3,1) = Mxz
  m_cart(2,3) = Myz
  m_cart(3,2) = Myz

  end subroutine gf_rotate_moment_tensor

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_moment_contract(m_cart,eps,val)

! SUM_pq M_pq eps_pq, from the Voigt-packed strain
!
! The off-diagonal terms are counted twice because each Voigt slot stands
! for two equal entries of the symmetric pair. See the module header.

  use constants, only: NDIM

  implicit none

  double precision, dimension(NDIM,NDIM), intent(in) :: m_cart
  double precision, dimension(GF_VOIGT), intent(in) :: eps
  double precision, intent(out) :: val

  val = m_cart(1,1)*eps(GF_XX) &
      + m_cart(2,2)*eps(GF_YY) &
      + m_cart(3,3)*eps(GF_ZZ) &
      + 2.d0*m_cart(1,2)*eps(GF_XY) &
      + 2.d0*m_cart(1,3)*eps(GF_XZ) &
      + 2.d0*m_cart(2,3)*eps(GF_YZ)

  end subroutine gf_moment_contract

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_moment_contract_full(m_cart,eps_full,val)

! the same contraction, written as the plain 3x3 double sum
!
! Deliberately naive, and deliberately not used in production: it is the
! oracle for the Voigt form above, which carries the factor-2 weights that
! are the easiest thing here to get wrong. tests/gf3d/test_gf_strain.f90
! asserts the two agree.

  use constants, only: NDIM

  implicit none

  double precision, dimension(NDIM,NDIM), intent(in) :: m_cart
  double precision, dimension(NDIM,NDIM), intent(in) :: eps_full
  double precision, intent(out) :: val

  ! local parameters
  integer :: p,q

  val = 0.d0
  do q = 1,NDIM
    do p = 1,NDIM
      val = val + m_cart(p,q)*eps_full(p,q)
    enddo
  enddo

  end subroutine gf_moment_contract_full

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_moment_unit_tensor(ivoigt,m_cart)

! the Cartesian moment tensor of a unit Voigt component
!
! Slot 4, 5 or 6 sets *both* off-diagonal entries, because a moment tensor
! is symmetric and the Voigt component names the pair. That is the
! convention under which Stage 6's six partials sum back to the full
! seismogram, and it is what makes that sum an exact identity rather than
! an approximation.

  use constants, only: NDIM

  implicit none

  integer, intent(in) :: ivoigt
  double precision, dimension(NDIM,NDIM), intent(out) :: m_cart

  m_cart(:,:) = 0.d0

  select case (ivoigt)
  case (GF_XX) ; m_cart(1,1) = 1.d0
  case (GF_YY) ; m_cart(2,2) = 1.d0
  case (GF_ZZ) ; m_cart(3,3) = 1.d0
  case (GF_XY) ; m_cart(1,2) = 1.d0 ; m_cart(2,1) = 1.d0
  case (GF_XZ) ; m_cart(1,3) = 1.d0 ; m_cart(3,1) = 1.d0
  case (GF_YZ) ; m_cart(2,3) = 1.d0 ; m_cart(3,2) = 1.d0
  end select

  end subroutine gf_moment_unit_tensor

  end module gf_moment
