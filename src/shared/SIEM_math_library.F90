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


! Spectral-Infinite Element Method (SIEM)
! helper module functions
!
! note: using lower-case naming for module names ("siem_math_*" instead of "SIEM_math_*") to avoid problems
!       with different compilers (like Cray) when generating the *.mod files.

module siem_math_library

  use constants, only: CUSTOM_REAL

  implicit none

  private

  public :: determinant
  public :: invert

  public :: i_uniinv
  public :: compute_g_gradg
  public :: compute_g_gradg_elliptical

  !public :: dotmat

contains

  function determinant(jac) result(det)

  ! this function returns the determinant of a 3x3 jacobian matrix.

  implicit none
  real(kind=CUSTOM_REAL),intent(in) :: jac(3,3)
  real(kind=CUSTOM_REAL) :: det

  det = jac(1,1) * (jac(2,2)*jac(3,3) - jac(3,2)*jac(2,3))
  det = det - jac(1,2) * (jac(2,1)*jac(3,3) - jac(3,1)*jac(2,3))
  det = det + jac(1,3) * (jac(2,1)*jac(3,2) - jac(3,1)*jac(2,2))

  end function determinant

  !=======================================================

! not used yet...

!  function determinant_generic(jac) result(det)
!
!  ! this generic function returns the determinant of a 1x1, 2x2 or 3x3 jacobian matrix.
!
!  ! the routine is following the idea from
!  ! Smith and Griffiths (2004): Programming the finite element method
!  ! see: https://github.com/ParaFEM/ParaFEM
!
!  implicit none
!  real(kind=CUSTOM_REAL),intent(in) :: jac(:,:)
!  real(kind=CUSTOM_REAL) :: det
!  integer :: it
!
!  it = ubound(jac,1)
!
!  select case(it)
!  case(1)
!    det = 1.0_CUSTOM_REAL
!  case(2)
!    det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
!  case(3)
!    det = jac(1,1)*(jac(2,2)*jac(3,3) - jac(3,2)*jac(2,3))
!    det = det - jac(1,2)*(jac(2,1)*jac(3,3) - jac(3,1)*jac(2,3))
!    det = det + jac(1,3)*(jac(2,1)*jac(3,2) - jac(3,1)*jac(2,2))
!  case default
!    stop 'ERROR: wrong dimension for jacobian matrix!'
!  end select
!
!  return
!
!  end function determinant_generic

  !=======================================================

  subroutine invert(matrix)

  ! this subroutine inverts a small square 3x3 matrix onto itself.

  implicit none
  real(kind=CUSTOM_REAL),intent(inout) :: matrix(3,3)
  real(kind=CUSTOM_REAL) :: det,detInv
  real(kind=CUSTOM_REAL) :: j11,j12,j13,j21,j22,j23,j31,j32,j33

  det = matrix(1,1) * (matrix(2,2)*matrix(3,3) - matrix(3,2)*matrix(2,3))
  det = det - matrix(1,2) * (matrix(2,1)*matrix(3,3) - matrix(3,1)*matrix(2,3))
  det = det + matrix(1,3) * (matrix(2,1)*matrix(3,2) - matrix(3,1)*matrix(2,2))

  if (det == 0.0_CUSTOM_REAL) then
    matrix(:,:) = 0.0_CUSTOM_REAL
    return
  endif

  detInv = 1.0_CUSTOM_REAL / det

  j11 = matrix(2,2)*matrix(3,3) - matrix(3,2)*matrix(2,3)
  j21 = -matrix(2,1)*matrix(3,3) + matrix(3,1)*matrix(2,3)
  j31 = matrix(2,1)*matrix(3,2) - matrix(3,1)*matrix(2,2)
  j12 = -matrix(1,2)*matrix(3,3) + matrix(3,2)*matrix(1,3)
  j22 = matrix(1,1)*matrix(3,3) - matrix(3,1)*matrix(1,3)
  j32 = -matrix(1,1)*matrix(3,2) + matrix(3,1)*matrix(1,2)
  j13 = matrix(1,2)*matrix(2,3) - matrix(2,2)*matrix(1,3)
  j23 = -matrix(1,1)*matrix(2,3) + matrix(2,1)*matrix(1,3)
  j33 = matrix(1,1)*matrix(2,2) - matrix(2,1)*matrix(1,2)

  matrix(1,1) = j11
  matrix(1,2) = j12
  matrix(1,3) = j13
  matrix(2,1) = j21
  matrix(2,2) = j22
  matrix(2,3) = j23
  matrix(3,1) = j31
  matrix(3,2) = j32
  matrix(3,3) = j33

  matrix(:,:) = matrix(:,:) * detInv

  end subroutine invert

  !=======================================================

! not used yet...

!  subroutine invert_generic(matrix)
!
!  ! this subroutine inverts a small square matrix onto itself.
!
!  ! the routine is following the idea from
!  ! Smith and Griffiths (2004): Programming the finite element method
!  ! see: https://github.com/ParaFEM/ParaFEM
!
!  implicit none
!  real(kind=CUSTOM_REAL),intent(inout) :: matrix(:,:)
!  real(kind=CUSTOM_REAL) :: det,j11,j12,j13,j21,j22,j23,j31,j32,j33,con
!  integer :: ndim,i,k
!
!  ndim = ubound(matrix,1)
!
!  if (ndim == 2) then
!    det = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
!    j11 = matrix(1,1)
!    matrix(1,1) = matrix(2,2)
!    matrix(2,2) = j11
!    matrix(1,2) = -matrix(1,2)
!    matrix(2,1) = -matrix(2,1)
!    matrix = matrix/det
!
!  else if (ndim == 3) then
!    det = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
!    det = det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
!    det = det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
!
!    j11 = matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
!    j21 = -matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
!    j31 = matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
!    j12 = -matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
!    j22 = matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
!    j32 = -matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
!    j13 = matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
!    j23 = -matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
!    j33 = matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
!
!    matrix(1,1) = j11
!    matrix(1,2) = j12
!    matrix(1,3) = j13
!    matrix(2,1) = j21
!    matrix(2,2) = j22
!    matrix(2,3) = j23
!    matrix(3,1) = j31
!    matrix(3,2) = j32
!    matrix(3,3) = j33
!    matrix = matrix/det
!
!  else
!    do k = 1,ndim
!      con = matrix(k,k)
!      matrix(k,k) = 1.0_CUSTOM_REAL
!      matrix(k,:) = matrix(k,:)/con
!      do i = 1,ndim
!        if (i /= k) then
!          con = matrix(i,k)
!          matrix(i,k) = 0.0_CUSTOM_REAL
!          matrix(i,:) = matrix(i,:) - matrix(k,:)*con
!        endif
!      enddo
!    enddo
!  endif
!
!  end subroutine invert_generic
!
  !=======================================================


! not used ...

!  function dotmat(m,n,x1,x2) result(dotm)
!
!  implicit none
!  integer,intent(in) :: m,n
!  real(kind=CUSTOM_REAL),intent(in) :: x1(m,n),x2(m,n)
!  real(kind=CUSTOM_REAL) :: dotm(n)
!  integer :: i
!
!  do i = 1,n
!    dotm(i) = dot_product(x1(:,i),x2(:,i))
!  enddo
!
!  return
!
!  end function dotmat

  !=======================================================

  ! Author: Michel Olagnon
  ! orderpack 2.0
  ! source: http://www.Fortran-2000.com/rank/

  subroutine i_uniinv (XDONT, IGOEST)

  ! UNIINV = Merge-sort inverse ranking of an array, with removal of duplicate entries.

  ! this routine is similar to pure merge-sort ranking, but on
  ! the last pass, it sets indices in IGOEST to the rank
  ! of the value in the ordered set with duplicates removed.
  ! for performance reasons, the first 2 passes are taken
  ! out of the standard loop, and use dedicated coding.
  implicit none
  integer,intent(in)  :: XDONT(:)
  integer,intent(out) :: IGOEST(:)

  integer :: XTST, XDONA, XDONB
  integer, dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
  integer :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
  integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

  NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
  select case (NVAL)
  case (:0)
    return
  case (1)
    IGOEST (1) = 1
    return
  case default
    continue
  end select

  ! fill-in the index array, creating ordered couples
  do IIND = 2, NVAL, 2
    if (XDONT(IIND-1) < XDONT(IIND)) then
      IRNGT (IIND-1) = IIND - 1
      IRNGT (IIND) = IIND
    else
      IRNGT (IIND-1) = IIND
      IRNGT (IIND) = IIND - 1
    endif
  enddo
  if (modulo(NVAL,2) /= 0) then
    IRNGT (NVAL) = NVAL
  endif

  ! we will now have ordered subsets A - B - A - B - ...
  ! and merge A and B couples into     C   -   C   - ...
  LMTNA = 2
  LMTNC = 4

  ! first iteration. The length of the ordered subsets goes from 2 to 4
  do
    if (NVAL <= 4) Exit
    ! loop on merges of A and B into C
    do IWRKD = 0, NVAL - 1, 4
      if ((IWRKD+4) > NVAL) then
        if ((IWRKD+2) >= NVAL) Exit
        !   1 2 3
        if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
        !   1 3 2
        if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
          IRNG2 = IRNGT (IWRKD+2)
          IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
          IRNGT (IWRKD+3) = IRNG2
          !   3 1 2
        else
          IRNG1 = IRNGT (IWRKD+1)
          IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
          IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
          IRNGT (IWRKD+2) = IRNG1
        endif
        exit
      endif
      !   1 2 3 4
      if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
      !   1 3 x x
      if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
        IRNG2 = IRNGT (IWRKD+2)
        IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
        if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
          !   1 3 2 4
          IRNGT (IWRKD+3) = IRNG2
        else
          !   1 3 4 2
          IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
          IRNGT (IWRKD+4) = IRNG2
        endif
        !   3 x x x
      else
        IRNG1 = IRNGT (IWRKD+1)
        IRNG2 = IRNGT (IWRKD+2)
        IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
        if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
          IRNGT (IWRKD+2) = IRNG1
          if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
            !   3 1 2 4
            IRNGT (IWRKD+3) = IRNG2
          else
            !   3 1 4 2
            IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
            IRNGT (IWRKD+4) = IRNG2
          endif
        else
          !   3 4 1 2
          IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
          IRNGT (IWRKD+3) = IRNG1
          IRNGT (IWRKD+4) = IRNG2
        endif
      endif
    enddo

  ! the Cs become As and Bs
    LMTNA = 4
    Exit
  enddo

  ! iteration loop. Each time, the length of the ordered subsets
  ! is doubled.
  do
    if (2*LMTNA >= NVAL) Exit
    IWRKF = 0
    LMTNC = 2 * LMTNC

    ! loop on merges of A and B into C
    do
      IWRK = IWRKF
      IWRKD = IWRKF + 1
      JINDA = IWRKF + LMTNA
      IWRKF = IWRKF + LMTNC
      if (IWRKF >= NVAL) then
        if (JINDA >= NVAL) Exit
        IWRKF = NVAL
      endif
      IINDA = 1
      IINDB = JINDA + 1

      ! one steps in the C subset, that we create in the final rank array
      ! make a copy of the rank array for the iteration
      JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
      XDONA = XDONT (JWRKT(IINDA))
      XDONB = XDONT (IRNGT(IINDB))
      do
        IWRK = IWRK + 1
        ! we still have unprocessed values in both A and B
        if (XDONA > XDONB) then
          IRNGT (IWRK) = IRNGT (IINDB)
          IINDB = IINDB + 1
          if (IINDB > IWRKF) then
            ! only A still with unprocessed values
            IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
            Exit
          endif
          XDONB = XDONT (IRNGT(IINDB))
        else
          IRNGT (IWRK) = JWRKT (IINDA)
          IINDA = IINDA + 1
          if (IINDA > LMTNA) Exit! Only B still with unprocessed values
          XDONA = XDONT (JWRKT(IINDA))
        endif

      enddo
    enddo

    ! the Cs become As and Bs
    LMTNA = 2 * LMTNA
  enddo

  ! last merge of A and B into C, with removal of duplicates.
  IINDA = 1
  IINDB = LMTNA + 1
  NUNI = 0

  ! one steps in the C subset, that we create in the final rank array
  JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
  if (IINDB <= NVAL) then
    XTST = i_nearless(Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
  else
    XTST = i_nearless(XDONT(JWRKT(1)))
  endif

  do IWRK = 1, NVAL
    ! we still have unprocessed values in both A and B
    if (IINDA <= LMTNA) then
      if (IINDB <= NVAL) then
        if (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) then
          IRNG = IRNGT (IINDB)
          IINDB = IINDB + 1
        else
          IRNG = JWRKT (IINDA)
          IINDA = IINDA + 1
        endif
      else
        ! only A still with unprocessed values
        IRNG = JWRKT (IINDA)
        IINDA = IINDA + 1
      endif
    else
      ! only B still with unprocessed values
      IRNG = IRNGT (IWRK)
    endif
    if (XDONT(IRNG) > XTST) then
      XTST = XDONT (IRNG)
      NUNI = NUNI + 1
    endif
    IGOEST (IRNG) = NUNI
  enddo

  end subroutine i_uniinv

  !=======================================================

  function i_nearless (XVAL) result (I_nl)

  ! nearest value less than given value

  implicit none
  integer,intent(in) :: XVAL
  integer :: I_nl

  I_nl = XVAL - 1

  end function i_nearless

  !=======================================================

  subroutine dspherical_unitvect(theta,phi,unitr,unittheta,unitphi)

  ! this function computes unit vectors in spherical coordinates at a point
  ! defined by (r,theta,phi)

  implicit none
  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: unitr(3),unittheta(3),unitphi(3)
  double precision :: ctheta,cphi,stheta,sphi
  double precision,parameter :: zero = 0.d0

  ctheta = cos(theta); stheta = sin(theta)
  cphi = cos(phi);     sphi = sin(phi)

  unitr(1) = stheta*cphi;   unitr(2) = stheta*sphi;   unitr(3) = ctheta
  unittheta(1) = ctheta*cphi; unittheta(2) = ctheta*sphi; unittheta(3) = -stheta
  unitphi(1) = -sphi;   unitphi(2) = cphi;    unitphi(3) = zero

  end subroutine dspherical_unitvect

  !=======================================================

! not used yet...

!  subroutine spherical_unitvect(theta,phi,unitr,unittheta,unitphi)
!
!  ! this function computes unit vectors in spherical coordinates at a point
!  ! defined by (r,theta,phi)
!
!  implicit none
!  real(kind=CUSTOM_REAL),intent(in) :: theta,phi
!  real(kind=CUSTOM_REAL),intent(out) :: unitr(3),unittheta(3),unitphi(3)
!  real(kind=CUSTOM_REAL) :: ctheta,cphi,stheta,sphi
!  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL
!
!  ctheta = cos(theta); stheta = sin(theta)
!  cphi = cos(phi);     sphi = sin(phi)
!
!  unitr(1) = stheta*cphi;   unitr(2) = stheta*sphi;   unitr(3) = ctheta
!  unittheta(1) = ctheta*cphi; unittheta(2) = ctheta*sphi; unittheta(3) = -stheta
!  unitphi(1) = -sphi;   unitphi(2) = cphi;    unitphi(3) = zero
!
!  end subroutine spherical_unitvect

  !=======================================================

  subroutine dlegendreP2_costheta(theta,P2,dP2,d2P2)

  implicit none
  double precision,intent(in) :: theta ! radian
  double precision,intent(out) :: P2,dP2,d2P2
  double precision :: ctheta
  double precision,parameter :: half = 0.5d0,one = 1.d0,two = 2.d0,three = 3.d0

  ctheta = cos(theta)

  P2 = half*(three*ctheta*ctheta-one)
  dP2 = -three*sin(theta)*ctheta
  d2P2 = -3*cos(two*theta)

  end subroutine dlegendreP2_costheta

  !=======================================================

! not used yet...

!  subroutine legendreP2_costheta(theta,P2,dP2,d2P2)
!
!  implicit none
!  real(kind=CUSTOM_REAL),intent(in) :: theta ! radian
!  real(kind=CUSTOM_REAL),intent(out) :: P2,dP2,d2P2
!  real(kind=CUSTOM_REAL) :: ctheta
!  real(kind=CUSTOM_REAL),parameter :: half=0.5_CUSTOM_REAL,one=1.0_CUSTOM_REAL,two=2.0_CUSTOM_REAL,three=3.0_CUSTOM_REAL
!
!  ctheta = cos(theta)
!
!  P2 = half*(three*ctheta*ctheta-one)
!  dP2 = -three*sin(theta)*ctheta
!  d2P2 = -3*cos(two*theta)
!
!  end subroutine legendreP2_costheta

  !=======================================================

  subroutine compute_g_gradg_elliptical(ndim,r,theta,phi,rho,dotrho,g0,eps,eta,twothirdOmega2,g,gradg)

  implicit none
  integer,intent(in) :: ndim
  double precision,intent(in) :: r,theta,phi,rho,dotrho,g0,eps,eta,twothirdOmega2
  double precision,intent(out) :: g(ndim)
  double precision,optional,intent(out) :: gradg(6)

  double precision,parameter :: two = 2.d0,four = 4.d0,two_third = two/3.d0

  double precision :: hmat(ndim,ndim)
  double precision :: lfac,rinv
  double precision :: cotthetaXdP2,doteps,ddoteps,doteta,dotg
  double precision :: facP2,P2,dP2,d2P2
  double precision :: unitr(ndim,1),unittheta(ndim,1),unitphi(ndim,1)
  double precision :: unitrT(1,ndim),unitthetaT(1,ndim),unitphiT(1,ndim)

  ! compute unit vectors
  call dspherical_unitvect(theta,phi,unitr,unittheta,unitphi)

  unitrT(1,:) = unitr(:,1);
  unitthetaT(1,:) = unittheta(:,1)
  unitphiT(1,:) = unitphi(:,1)

  call dlegendreP2_costheta(theta,P2,dP2,d2P2)

  rinv = 1.d0/r

  dotg = four*rho-two*rinv*g0 !dotg = four_pi_G*rho-two*rinv*g_ss
  doteps = eta*eps*rinv !eta*eps/r
  ddoteps = 6.d0*rinv*rinv*eps-8.0d0*rho*(doteps+rinv*eps)/g0 !two*four

  !doteta=doteps/eps-0.5_CUSTOM_REAL*r0*rinv*rinv*doteps*doteps+r0*ddoteps/eps
  !doteta=doteps/eps-r*doteps*doteps/(eps*eps)+r*ddoteps/eps
  doteta = doteps/eps+r*(eps*ddoteps-doteps*doteps)/(eps*eps)

  !cottheta=cos(theta)/sin(theta)
  cotthetaXdp2 = -3.d0*cos(theta)*cos(theta)

  ! lfac=g_ss+two_third*(eta*eps*g_ss+four_pi_G*r0*rho*eps-eps*g_ss)*P2-two_third*Omega2*r0
  lfac = g0+two_third*(eta*eps*g0+four*r*rho*eps-eps*g0)*P2-twothirdOmega2*r

  ! compute g
  g = -lfac*unitr(:,1)-two_third*eps*g0*dP2*unittheta(:,1)

  if (.not.present(gradg)) return

  ! compute grad g
  facP2 = (doteta*eps*g0+eta*doteps*g0+eta*eps*dotg+four*(rho*eps+r*dotrho*eps+r*rho*doteps)-doteps*g0-eps*dotg)

  !hmat=-(dotg+two_third*facP2*P2-two_third*Omega2)*matmul(unitr,unitrT)     &
  !     -two_third*(doteps*g0+eps+dotg)*dP2*(matmul(unitr,unitthetaT)+matmul(unittheta,unitrT))     &
  !     -rinv*(lfac-two_third*Omega2*r+two_third*rinv*eps*g0*d2p2)*matmul(unittheta,unitthetaT)  &
  !     -rinv*(lfac-two_third*Omega2*r+two_third*rinv*eps*g0*cotthetaXdp2)*matmul(unitphi,unitphiT)
  hmat = -(dotg+two_third*facP2*P2-twothirdOmega2)*matmul(unitr,unitrT)     &
        -two_third*(doteps*g0+eps+dotg)*dP2*(matmul(unitr,unitthetaT)+matmul(unittheta,unitrT))     &
        -rinv*(lfac+two_third*rinv*eps*g0*d2p2)*matmul(unittheta,unitthetaT)  &
        -rinv*(lfac+two_third*rinv*eps*g0*cotthetaXdp2)*matmul(unitphi,unitphiT)

  gradg = (/ hmat(1,1),hmat(2,2),hmat(3,3),hmat(1,2),hmat(1,3),hmat(2,3) /)

  end subroutine compute_g_gradg_elliptical

  !=======================================================

  subroutine compute_g_gradg(ndim,r,theta,phi,rho,g0,g,gradg)

  implicit none
  integer,intent(in) :: ndim
  double precision,intent(in) :: r,theta,phi,rho,g0
  double precision,intent(out) :: g(ndim)
  double precision,optional,intent(out) :: gradg(6)

  double precision,parameter :: two = 2.d0,four = 4.d0
  double precision :: hmat(ndim,ndim)
  double precision :: lfac,rinv
  double precision :: dotg
  double precision :: ctheta,P2 !,dP2,d2P2
  double precision :: unitr(ndim,1),unittheta(ndim,1),unitphi(ndim,1)
  double precision :: unitrT(1,ndim),unitthetaT(1,ndim),unitphiT(1,ndim)

  ! compute unit vectors
  call dspherical_unitvect(theta,phi,unitr,unittheta,unitphi)

  unitrT(1,:) = unitr(:,1);
  unitthetaT(1,:) = unittheta(:,1)
  unitphiT(1,:) = unitphi(:,1)

  !call dlegendreP2_costheta(theta,P2,dP2,d2P2)

  ctheta = cos(theta)
  P2 = 0.5d0*(3.0d0*ctheta*ctheta-1.0d0)

  rinv = 1.d0/r
  dotg = four*rho-two*rinv*g0 !dotg = four_pi_G*rho-two*rinv*g_ss

  lfac = g0

  ! compute g
  g = -lfac*unitr(:,1)

  if (.not.present(gradg)) return

  ! compute grad g
  hmat = -dotg*matmul(unitr,unitrT)-rinv*lfac*matmul(unittheta,unitthetaT) &
        -rinv*lfac*matmul(unitphi,unitphiT)

  gradg = (/ hmat(1,1),hmat(2,2),hmat(3,3),hmat(1,2),hmat(1,3),hmat(2,3) /)

  end subroutine compute_g_gradg

end module siem_math_library


!
!-----------------------------------------------------------------------------------
!


! MPI math library
module siem_math_library_mpi

  use constants, only: CUSTOM_REAL

  implicit none

  private

  public :: maxscal
  public :: minvec
  public :: maxvec

  public :: dot_product_all_proc


  ! global sum of a scalar in all processors
  !interface sumscal
  !  module procedure isumscal
  !  module procedure fsumscal
  !end interface

  ! global maximum of a scalar in all processors
  !interface minscal
  !  module procedure iminscal
  !  module procedure fminscal
  !end interface

  ! global maximum of a scalar in all processors
  interface maxscal
    module procedure imaxscal
    module procedure fmaxscal
  end interface

  ! global maximum of a vector in all processors
  interface maxvec
    module procedure imaxvec
    module procedure fmaxvec
  end interface

  ! global minimum of a scalar in all processors
  interface minvec
    module procedure iminvec
    module procedure fminvec
  end interface

contains

! not used yet ...

!  function iminscal(scal) result(gmin)
!  !
!  ! this finds a global minimum of a scalar across the processors
!  !
!  implicit none
!  integer,intent(in)::scal
!  integer :: gmin
!
!  call min_all_all_i(scal,gmin)
!
!  return
!  end function iminscal

  !=======================================================

! not used yet ...

!  function fminscal(scal) result(gmin)
!  !
!  ! this finds a global minimum of a scalar across the processors
!  !
!  implicit none
!  real(kind=CUSTOM_REAL),intent(in)::scal
!  real(kind=CUSTOM_REAL) :: gmin
!
!  call min_all_all_cr(scal,gmin)
!
!  return
!  end function fminscal

  !=======================================================

  function imaxscal(scal) result(gmax)
  !
  ! this finds a global maximum of a scalar across the processors
  !
  implicit none
  integer,intent(in)::scal
  integer :: gmax

  call max_all_all_i(scal,gmax)

  return
  end function imaxscal

  !=======================================================

  function fmaxscal(scal) result(gmax)
  !
  ! this finds a global maximum of a scalar across the processors
  !
  implicit none
  real(kind=CUSTOM_REAL),intent(in)::scal
  real(kind=CUSTOM_REAL) :: gmax

  call max_all_all_cr(scal,gmax)

  return
  end function fmaxscal

  !=======================================================

  function imaxvec(vec) result(gmax)
  implicit none
  integer,intent(in)::vec(:)
  integer :: lmax,gmax ! local and global

  lmax = maxval(vec)

  call max_all_all_i(lmax,gmax)

  return
  end function imaxvec

  !=======================================================

  function fmaxvec(vec) result(gmax)
  implicit none
  real(kind=CUSTOM_REAL),intent(in)::vec(:)
  real(kind=CUSTOM_REAL) :: lmax,gmax ! local and global

  lmax = maxval(vec)

  call max_all_all_cr(lmax,gmax)

  return
  end function fmaxvec

  !=======================================================

  function iminvec(vec) result(gmin)
  implicit none
  integer,intent(in)::vec(:)
  integer :: lmin,gmin ! local and global

  lmin = minval(vec)

  call min_all_all_i(lmin,gmin)

  return
  end function iminvec

  !=======================================================

  function fminvec(vec) result(gmin)
  implicit none
  real(kind=CUSTOM_REAL),intent(in)::vec(:)
  real(kind=CUSTOM_REAL) :: lmin,gmin ! local and global

  lmin = minval(vec)

  call min_all_all_cr(lmin,gmin)

  return
  end function fminvec

  !=======================================================

! not used yet...

!  function isumscal(scal) result(gsum)
!  !
!  ! this finds a global summation of a scalar across the processors
!  !
!  implicit none
!  integer,intent(in)::scal
!  integer :: gsum
!
!  call sum_all_all_i(scal,gsum)
!
!  return
!  end function isumscal

  !=======================================================

! not used yet...

!  function fsumscal(scal) result(gsum)
!  !
!  ! this finds a global summation of a scalar across the processors
!  !
!  implicit none
!  real(kind=CUSTOM_REAL),intent(in)::scal
!  real(kind=CUSTOM_REAL) :: gsum
!
!  call sum_all_all_cr(scal,gsum)
!
!  return
!  end function fsumscal

  !=======================================================

  function dot_product_all_proc(vec1,vec2) result(gdot)
  !
  ! this finds global dot product of two vectors across the processors
  !
  implicit none
  real(kind=CUSTOM_REAL),intent(in)::vec1(:),vec2(:)
  real(kind=CUSTOM_REAL) :: ldot,gdot

  ! find local dot
  ldot = dot_product(vec1,vec2)

  call sum_all_all_cr(ldot,gdot)

  return
  end function dot_product_all_proc

end module siem_math_library_mpi


!
!-----------------------------------------------------------------------------------
!

! NOTES:
!  - gll_points not needed can be removed
!  - can be make faster using orthoginality of shape functions
! this module contains routines to compute Gauss-Legendre-Lobatto quadrature

module siem_gll_library

  implicit none

  private

  ! double precision
  integer,parameter :: kdble = selected_real_kind(15)

  ! number of shape function for 8-noded Hex
  integer,parameter :: NGNOD_INF = 8

  public :: kdble
  public :: NGNOD_INF

  public :: dshape_function_hex8
  public :: dshape_function_hex8_point
  public :: gll_quadrature
  public :: gll_quadrature3inNGLL
  public :: gll_lagrange3d_point
  public :: lagrange1dGLLAS
  public :: lagrange1dGENAS
  public :: zwgljd

contains

! this subroutines computes derivatives of the shape fucntions at gll
! points. the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention

  subroutine dshape_function_hex8(ndim,ngnod,ngllx,nglly,ngllz,ngll,xigll,etagll, &
                                  zetagll,dshape_hex8)

  implicit none
  integer,intent(in) :: ndim,ngnod,ngllx,nglly,ngllz,ngll

  ! gauss-lobatto-legendre points of integration
  real(kind=kdble),intent(in) :: xigll(ngllx),etagll(nglly),zetagll(ngllz)

  ! derivatives of the 3d shape functions
  real(kind=kdble),intent(out) :: dshape_hex8(ndim,ngnod,ngll)

  integer :: i,j,k,i_gnod,igll

  ! location of the nodes of the 3d quadrilateral elements
  !real(kind=kdble) :: xi,eta,gamma
  real(kind=kdble) :: xip,xim,etap,etam,zetap,zetam

  ! for checking the 3d shape functions
  real(kind=kdble) :: sum_dshapexi,sum_dshapeeta,sum_dshapezeta

  real(kind=kdble),parameter :: one=1.0_kdble,one_eighth = 0.125_kdble, &
  zero = 0.0_kdble,zerotol = 1.0e-12_kdble !1.0e-12_kdble WARNING: please correct immediately

  ! check that the parameter file is correct
  if (ngnod /= NGNOD_INF) then
    stop 'ERROR: elements must have 8 geometrical nodes!'
  endif

  !print *,'after:',ndim,ngnod,ngll
  ! compute the derivatives of 3d shape functions
  !dshape_hex8=zero
  igll = 0
  do k = 1,ngllz
    zetap = one + zetagll(k)
    zetam = one - zetagll(k)
    do j = 1,nglly
      etap = one + etagll(j)
      etam = one - etagll(j)
      do i = 1,ngllx
        igll = igll+1

        !xi = xigll(i)
        !eta = etagll(j)
        !gamma = zetagll(k)

        xip = one + xigll(i)
        xim = one - xigll(i)

        dshape_hex8(1,1,igll) = - one_eighth*etam*zetam
        dshape_hex8(1,2,igll) = one_eighth*etam*zetam
        dshape_hex8(1,3,igll) = one_eighth*etap*zetam
        dshape_hex8(1,4,igll) = - one_eighth*etap*zetam
        dshape_hex8(1,5,igll) = - one_eighth*etam*zetap
        dshape_hex8(1,6,igll) = one_eighth*etam*zetap
        dshape_hex8(1,7,igll) = one_eighth*etap*zetap
        dshape_hex8(1,8,igll) = - one_eighth*etap*zetap

        dshape_hex8(2,1,igll) = - one_eighth*xim*zetam
        dshape_hex8(2,2,igll) = - one_eighth*xip*zetam
        dshape_hex8(2,3,igll) = one_eighth*xip*zetam
        dshape_hex8(2,4,igll) = one_eighth*xim*zetam
        dshape_hex8(2,5,igll) = - one_eighth*xim*zetap
        dshape_hex8(2,6,igll) = - one_eighth*xip*zetap
        dshape_hex8(2,7,igll) = one_eighth*xip*zetap
        dshape_hex8(2,8,igll) = one_eighth*xim*zetap

        dshape_hex8(3,1,igll) = - one_eighth*xim*etam
        dshape_hex8(3,2,igll) = - one_eighth*xip*etam
        dshape_hex8(3,3,igll) = - one_eighth*xip*etap
        dshape_hex8(3,4,igll) = - one_eighth*xim*etap
        dshape_hex8(3,5,igll) = one_eighth*xim*etam
        dshape_hex8(3,6,igll) = one_eighth*xip*etam
        dshape_hex8(3,7,igll) = one_eighth*xip*etap
        dshape_hex8(3,8,igll) = one_eighth*xim*etap

      enddo
    enddo
  enddo

  ! check the shape functions and their derivatives
  do i = 1,ngll
    sum_dshapexi = zero
    sum_dshapeeta = zero
    sum_dshapezeta = zero

    do i_gnod = 1,ngnod
      sum_dshapexi = sum_dshapexi + dshape_hex8(1,i_gnod,i)
      sum_dshapeeta = sum_dshapeeta + dshape_hex8(2,i_gnod,i)
      sum_dshapezeta = sum_dshapezeta + dshape_hex8(3,i_gnod,i)
      !print *,sum_dshapexi,sum_dshapeeta,sum_dshapezeta
    enddo

    ! sum of derivative of shape functions should be zero
    if (abs(sum_dshapexi) > zerotol) then
      stop 'ERROR: derivative xi shape functions!'
    endif

    if (abs(sum_dshapeeta) > zerotol) then
      stop 'ERROR: derivative eta shape functions!'
    endif

    if (abs(sum_dshapezeta) > zerotol) then
      print *,'ERROR: derivative gamma shape functions!'
      print *,ngllx,nglly,ngllz,sum_dshapexi,sum_dshapeeta,sum_dshapezeta,zerotol
      print *,xigll
      print *,etagll
      print *,zetagll
      if (all(xigll == etagll)) print *,'yes0'
      if (all(xigll == zetagll)) print *,'yes1'
      stop
    endif
  enddo

  end subroutine dshape_function_hex8

!
!=======================================================
!

! This subroutines computes derivatives of the shape fucntions at given point.
! the 8-noded hexahedra is conformed to the exodus/cubit numbering
! convention

!NOTE: dimension of dshape_hex8 is (ngnod,3) NOT (3,ngnode)

  subroutine dshape_function_hex8_point(ngnod,xi,eta,zeta,dshape_hex8)

  implicit none
  integer,intent(in) :: ngnod

  ! given point
  double precision,intent(in) :: xi,eta,zeta

  ! derivatives of the 3d shape functions
  double precision :: dshape_hex8(ngnod,3)

  integer :: i_gnod

  double precision :: xip,xim,etap,etam,zetap,zetam

  ! for checking the 3d shape functions
  double precision :: sum_dshapexi,sum_dshapeeta,sum_dshapezeta

  real(kind=kdble),parameter :: one=1.0_kdble,one_eighth = 0.125_kdble, &
                                zero = 0.0_kdble,zerotol = 1.0e-12_kdble !1.0e-12_kdble WARNING: please correct immediately

  ! check that the parameter file is correct
  if (ngnod /= NGNOD_INF) then
    stop 'ERROR: elements must have 8 geometrical nodes!'
  endif

  ! compute the derivatives of 3d shape functions
  zetap = one + zeta
  zetam = one - zeta
  etap =  one + eta
  etam =  one - eta
  xip =   one + xi
  xim =   one - xi

  dshape_hex8 = zero

  dshape_hex8(1,1) = - one_eighth*etam*zetam
  dshape_hex8(2,1) = one_eighth*etam*zetam
  dshape_hex8(3,1) = one_eighth*etap*zetam
  dshape_hex8(4,1) = - one_eighth*etap*zetam
  dshape_hex8(5,1) = - one_eighth*etam*zetap
  dshape_hex8(6,1) = one_eighth*etam*zetap
  dshape_hex8(7,1) = one_eighth*etap*zetap
  dshape_hex8(8,1) = - one_eighth*etap*zetap

  dshape_hex8(1,2) = - one_eighth*xim*zetam
  dshape_hex8(2,2) = - one_eighth*xip*zetam
  dshape_hex8(3,2) = one_eighth*xip*zetam
  dshape_hex8(4,2) = one_eighth*xim*zetam
  dshape_hex8(5,2) = - one_eighth*xim*zetap
  dshape_hex8(6,2) = - one_eighth*xip*zetap
  dshape_hex8(7,2) = one_eighth*xip*zetap
  dshape_hex8(8,2) = one_eighth*xim*zetap

  dshape_hex8(1,3) = - one_eighth*xim*etam
  dshape_hex8(2,3) = - one_eighth*xip*etam
  dshape_hex8(3,3) = - one_eighth*xip*etap
  dshape_hex8(4,3) = - one_eighth*xim*etap
  dshape_hex8(5,3) = one_eighth*xim*etam
  dshape_hex8(6,3) = one_eighth*xip*etam
  dshape_hex8(7,3) = one_eighth*xip*etap
  dshape_hex8(8,3) = one_eighth*xim*etap

  ! check the shape functions and their derivatives
  sum_dshapexi = zero
  sum_dshapeeta = zero
  sum_dshapezeta = zero

  do i_gnod = 1,ngnod
    sum_dshapexi = sum_dshapexi + dshape_hex8(i_gnod,1)
    sum_dshapeeta = sum_dshapeeta + dshape_hex8(i_gnod,2)
    sum_dshapezeta = sum_dshapezeta + dshape_hex8(i_gnod,3)
  enddo

  ! sum of derivative of shape functions should be zero
  if (abs(sum_dshapexi) > zerotol) then
    stop 'ERROR: derivative xi shape functions!'
  endif
  if (abs(sum_dshapeeta) > zerotol) then
    stop 'ERROR: derivative eta shape functions!'
  endif
  if (abs(sum_dshapezeta) > zerotol) then
    stop 'ERROR: derivative gamma shape functions!'
  endif

  end subroutine dshape_function_hex8_point

!
!===============================================================================
!

! this subroutine computes GLL quadrature points and weights for 3D

  subroutine gll_quadrature(ndim,ngllx,nglly,ngllz,ngll,gll_points,gll_weights, &
                            lagrange_gll,dlagrange_gll)

  implicit none
  integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll
  real(kind=kdble),dimension(ngll),intent(out) :: gll_weights
  real(kind=kdble),dimension(ndim,ngll),intent(out) :: gll_points
  real(kind=kdble),dimension(ngll,ngll),intent(out) :: lagrange_gll
  real(kind=kdble),dimension(ndim,ngll,ngll),intent(out) :: dlagrange_gll

  !real(kind=kdble),dimension(ngll,ngll) :: flagrange_gll
  !real(kind=kdble),dimension(ndim,ngll,ngll) ::fdlagrange_gll
  real(kind=kdble),parameter :: zero=0.0_kdble,one=1.0_kdble
  real(kind=kdble),parameter :: jacobi_alpha=zero,jacobi_beta=zero
  integer :: i,j,k,n
  integer :: ip,ipx,ipy,ipz ! integration points
  real(kind=kdble) :: xi,eta,zeta
  real(kind=kdble),dimension(ngllx) :: gllpx,gllwx ! GLL points and weights
  real(kind=kdble),dimension(nglly) :: gllpy,gllwy ! GLL points and weights
  real(kind=kdble),dimension(ngllz) :: gllpz,gllwz ! GLL points and weights
  real(kind=kdble),dimension(ngllx) :: lagrange_x,lagrange_dx
  real(kind=kdble),dimension(nglly) :: lagrange_y,lagrange_dy
  real(kind=kdble),dimension(ngllz) :: lagrange_z,lagrange_dz

  ! compute everything in indexed order

  ! get GLL points
  ! for alpha=beta=0, jacobi polynomial is legendre polynomial
  ! for ngllx=nglly=ngllz=ngll, need to call only once
  call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
  call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
  call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

  ip = 0
  do ipz = 1,ngllz
    do ipy = 1,nglly
      do ipx = 1,ngllx
        ip = ip+1
        ! integration points
        gll_points(1,ip)=gllpx(ipx)
        gll_points(2,ip)=gllpy(ipy)
        gll_points(3,ip)=gllpz(ipz)

        ! integration weights
        gll_weights(ip)=gllwx(ipx)*gllwy(ipy)*gllwz(ipz)
      enddo
    enddo
  enddo

  !! faster approach
  !ngllxy=ngllx*nglly
  !flagrange_gll=zero
  !fdlagrange_gll=zero
  !ip=0
  !! do ii=1,ngll ! ngllx*nglly*ngllz
  !do ipz=1,ngllz
  !  do ipy=1,nglly
  !    do ipx=1,ngllx
  !      ip=ip+1
  !      ! integration point
  !      xi=gll_points(1,ip)
  !      eta=gll_points(2,ip)
  !      zeta=gll_points(3,ip)

  !      ! interpolation function is orthogonal
  !      flagrange_gll(ip,ip)=one

  !      ! compute 1d lagrange polynomials on GLL points
  !      ! this can also be computed in a simple manner due to the orthogonality
  !      call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
  !      call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
  !      call lagrange1dGLL(ngllz,gllpz,zeta,lagrange_z,lagrange_dz)

  !      ! derivatives are nonzeros only on the lines intersecting the integration
  !      ! point

  !      ! x line
  !      inpz=ngllxy*(ipz-1); inpy=ngllx*(ipy-1)
  !      do npx=1,ngllx
  !        np=inpz+inpy+npx
  !        fdlagrange_gll(1,ip,np)=lagrange_dx(npx)
  !      enddo

  !      ! y line
  !      !inpz=ngllxy*(ipz-1) ! this was computed just above
  !      do npy=1,nglly
  !        np=inpz+ngllx*(npy-1)+ipx
  !        fdlagrange_gll(2,ip,np)=lagrange_dy(npy)
  !      enddo

  !      ! z line
  !      inpy=ngllx*(ipy-1)
  !      do npz=1,ngllz
  !        np=ngllxy*(npz-1)+inpy+ipx
  !        fdlagrange_gll(3,ip,np)=lagrange_dz(npz)
  !      enddo
  !    enddo
  !  enddo
  !enddo

  ! easier and general approach
  do ip = 1,ngll ! ngllx*nglly*ngllz
    xi = gll_points(1,ip)
    eta = gll_points(2,ip)
    zeta = gll_points(3,ip)

    ! compute 1d lagrange polynomials on GLL points
    ! this can also be computed in a simple manner due to the orthogonality as
    ! above
    call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
    call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
    call lagrange1dGLL(ngllz,gllpz,zeta,lagrange_z,lagrange_dz)
    n = 0
    do k = 1,ngllz
      do j = 1,nglly
        do i = 1,ngllx
          n = n+1
          lagrange_gll(ip,n) = lagrange_x(i)*lagrange_y(j)*lagrange_z(k)
          dlagrange_gll(1,ip,n) = lagrange_dx(i)*lagrange_y(j)*lagrange_z(k)
          dlagrange_gll(2,ip,n) = lagrange_x(i)*lagrange_dy(j)*lagrange_z(k)
          dlagrange_gll(3,ip,n) = lagrange_x(i)*lagrange_y(j)*lagrange_dz(k)
        enddo
      enddo
    enddo
  enddo

  !print *,maxval(abs(lagrange_gll-flagrange_gll))
  !print *,maxval(abs(dlagrange_gll-fdlagrange_gll))

  end subroutine gll_quadrature

!
!===========================================
!

! this subroutine computes lagrange polynomials and their derivatives defined on
! GLL points at an arbitrary point

  subroutine gll_lagrange3d_point(ndim,ngllx,nglly,ngllz,ngll,xi,eta,zeta, &
                                  lagrange_gll,dlagrange_gll)

  implicit none
  integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll
  double precision,intent(in) :: xi,eta,zeta
  real(kind=kdble),dimension(ngll),intent(out) :: lagrange_gll
  real(kind=kdble),dimension(ndim,ngll),intent(out) :: dlagrange_gll

  integer :: i,j,k,n
  real(kind=kdble),dimension(ngllx) :: lagrange_x,lagrange_dx
  real(kind=kdble),dimension(nglly) :: lagrange_y,lagrange_dy
  real(kind=kdble),dimension(ngllz) :: lagrange_z,lagrange_dz

  ! compute 1d lagrange polynomials
  call lagrange1d(ngllx,xi,lagrange_x,lagrange_dx)
  call lagrange1d(nglly,eta,lagrange_y,lagrange_dy)
  call lagrange1d(ngllz,zeta,lagrange_z,lagrange_dz)

  n = 0
  do k = 1,ngllz
    do j = 1,nglly
      do i = 1,ngllx
        n = n+1
        lagrange_gll(n) = lagrange_x(i)*lagrange_y(j)*lagrange_z(k)
        dlagrange_gll(1,n) = lagrange_dx(i)*lagrange_y(j)*lagrange_z(k)
        dlagrange_gll(2,n) = lagrange_x(i)*lagrange_dy(j)*lagrange_z(k)
        dlagrange_gll(3,n) = lagrange_x(i)*lagrange_y(j)*lagrange_dz(k)
      enddo
    enddo
  enddo

  end subroutine gll_lagrange3d_point

!
!===========================================
!

! subroutine below is only applicable ngllx=nglly=ngllz and NGLLX_INF=NGLLY_INF=NGLLZ_INF=3
! this subroutine computes GLL quadrature points and weights for 3D

  subroutine gll_quadrature3inNGLL(ndim,ngllx,ngll,gllpx,gllpx1,lagrange_gll)

  implicit none
  integer,intent(in) :: ndim,ngllx,ngll
  real(kind=kdble),intent(in) :: gllpx(ngllx),gllpx1(3)
  real(kind=kdble),dimension(ngll,27),intent(out) :: lagrange_gll

  real(kind=kdble),parameter :: zero=0.0_kdble,one=1.0_kdble
  integer :: i,j,k
  integer :: ip,ipx,ipy,ipz ! integration points
  integer :: np ! interpolation function points
  real(kind=kdble) :: xi,eta,zeta
  real(kind=kdble),dimension(ngllx) :: lagrange_x,lagrange_dx
  real(kind=kdble),dimension(ngllx) :: lagrange_y,lagrange_dy
  real(kind=kdble),dimension(ngllx) :: lagrange_z,lagrange_dz

  real(kind=kdble),dimension(ndim,ngll) :: gll_points

  ! compute everything in indexed order
  ip = 0
  do ipz = 1,ngllx
    do ipy = 1,ngllx
      do ipx = 1,ngllx
        ip = ip+1;
        ! integration points
        gll_points(1,ip)=gllpx(ipx)
        gll_points(2,ip)=gllpx(ipy)
        gll_points(3,ip)=gllpx(ipz)
      enddo
    enddo
  enddo

  ! easier and general approach
  do ip = 1,ngll ! ngllx*nglly*ngllz
    xi=gll_points(1,ip)
    eta=gll_points(2,ip)
    zeta=gll_points(3,ip)

    ! compute 1d lagrange polynomials on GLL points
    ! this can also be computed in a simple manner due to the orthogonality as
    ! above
    call lagrange1dGLL(3,gllpx1,xi,lagrange_x,lagrange_dx)
    call lagrange1dGLL(3,gllpx1,eta,lagrange_y,lagrange_dy)
    call lagrange1dGLL(3,gllpx1,zeta,lagrange_z,lagrange_dz)
    np = 0
    do k = 1,3
      do j = 1,3
        do i = 1,3
          np = np+1
          lagrange_gll(ip,np)=lagrange_x(i)*lagrange_y(j)*lagrange_z(k)
        enddo
      enddo
    enddo
  enddo

  end subroutine gll_quadrature3inNGLL

!
!===========================================
!

! this subroutine computes GLL quadrature points and weights for 2D

! not used yet ...

!  subroutine gll_quadrature2d(ndim,ngllx,nglly,ngll,gll_points2d,gll_weights2d, &
!                              lagrange_gll2d,dlagrange_gll2d)
!  implicit none
!  integer,intent(in) :: ndim,ngllx,nglly,ngll
!  real(kind=kdble),dimension(ngll),intent(out) :: gll_weights2d
!  real(kind=kdble),dimension(2,ngll),intent(out) :: gll_points2d
!  real(kind=kdble),dimension(ngll,ngll),intent(out) :: lagrange_gll2d
!  real(kind=kdble),dimension(ndim,ngll,ngll),intent(out) :: dlagrange_gll2d
!
!  real(kind=kdble),parameter :: jacobi_alpha=0.0_kdble,jacobi_beta=0.0_kdble
!  integer :: i,ii,j,n
!  real(kind=kdble) :: xi,eta !,zeta
!  real(kind=kdble),dimension(ngllx) :: gllpx,gllwx ! GLL points and weights
!  real(kind=kdble),dimension(nglly) :: gllpy,gllwy ! GLL points and weights
!  real(kind=kdble),dimension(ngllx) :: lagrange_x,lagrange_dx
!  real(kind=kdble),dimension(nglly) :: lagrange_y,lagrange_dy
!
!  ! compute everything in indexed order
!
!  ! GLL points and weights (source: http://mathworld.wolfram.com/lobattoquadrature.html)
!  !gllp(1)=-1.0_kdble   ; gllw(1)=1.0_kdble/3.0_kdble
!  !gllp(2)= 0.0_kdble   ; gllw(2)=4.0_kdble/3.0_kdble
!  !gllp(3)= 1.0_kdble   ; gllw(3)=gllw(1)
!
!  ! get GLL points
!  ! for alpha=beta=0, jacobi polynomial is legendre polynomial
!  ! for ngllx=nglly=ngllz=ngll, need to call only once
!  call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
!  call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
!
!  n = 0
!  do j = 1,nglly
!    do i = 1,ngllx
!      n = n+1
!      ! integration points
!      gll_points2d(1,n)=gllpx(i)
!      gll_points2d(2,n)=gllpy(j)
!
!      ! integration weights
!      gll_weights2d(n)=gllwx(i)*gllwy(j)
!    enddo
!  enddo
!
!  do ii = 1,ngll ! ngllx*nglly
!    xi=gll_points2d(1,ii)
!    eta=gll_points2d(2,ii)
!
!    ! compute 1d lagrange polynomials
!    ! this can also be computed in a simple manner due to the orthogonality
!    call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
!    call lagrange1dGLL(nglly,gllpy,eta,lagrange_y,lagrange_dy)
!    n = 0
!    do j = 1,nglly
!      do i = 1,ngllx
!        n = n+1
!        lagrange_gll2d(ii,n)=lagrange_x(i)*lagrange_y(j)
!        dlagrange_gll2d(1,ii,n)=lagrange_dx(i)*lagrange_y(j)
!        dlagrange_gll2d(2,ii,n)=lagrange_x(i)*lagrange_dy(j)
!      enddo
!    enddo
!  enddo
!
!  end subroutine gll_quadrature2d

!
!===========================================
!

! this subroutine computes GLL quadrature points and weights for 1D

! not used yet...

!  subroutine gll_quadrature1d(ndim,ngllx,ngll,gll_points1d,gll_weights1d, &
!                              lagrange_gll1d,dlagrange_gll1d)
!  implicit none
!  integer,intent(in) :: ndim,ngllx,ngll
!  real(kind=kdble),dimension(ngll),intent(out) :: gll_weights1d
!  real(kind=kdble),dimension(ndim,ngll),intent(out) :: gll_points1d
!  real(kind=kdble),dimension(ngll,ngll),intent(out) :: lagrange_gll1d
!  real(kind=kdble),dimension(ndim,ngll,ngll),intent(out) :: dlagrange_gll1d
!
!  real(kind=kdble),parameter :: jacobi_alpha=0.0_kdble,jacobi_beta=0.0_kdble
!  integer :: i,ii,n
!  real(kind=kdble) :: xi
!  real(kind=kdble),dimension(ngllx) :: gllpx,gllwx ! GLL points and weights
!  real(kind=kdble),dimension(ngllx) :: lagrange_x,lagrange_dx
!
!  ! compute everything in indexed order
!
!  ! GLL points and weights (source: http://mathworld.wolfram.com/lobattoquadrature.html)
!  !gllp(1)=-1.0_kdble   ; gllw(1)=1.0_kdble/3.0_kdble
!  !gllp(2)= 0.0_kdble   ; gllw(2)=4.0_kdble/3.0_kdble
!  !gllp(3)= 1.0_kdble   ; gllw(3)=gllw(1)
!
!  ! get GLL points
!  ! for alpha=beta=0, jacobi polynomial is legendre polynomial
!  ! for ngllx=nglly=ngllz=ngll, need to call only once
!  call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
!
!  n = 0
!  do i = 1,ngllx
!    n = n+1
!    ! integration points
!    gll_points1d(1,n)=gllpx(i)
!
!    ! integration weights
!    gll_weights1d(n)=gllwx(i)
!  enddo
!
!  do ii = 1,ngll ! ngllx
!    xi=gll_points1d(1,ii)
!
!    ! compute 1d lagrange polynomials
!    ! this can also be computed in a simple manner due to the orthogonality
!    call lagrange1dGLL(ngllx,gllpx,xi,lagrange_x,lagrange_dx)
!
!    n = 0
!    do i = 1,ngllx
!      n = n+1
!      lagrange_gll1d(ii,n)=lagrange_x(i)
!      dlagrange_gll1d(1,ii,n)=lagrange_dx(i)
!    enddo
!  enddo
!
!  end subroutine gll_quadrature1d

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

  subroutine lagrange1d(nenode,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenode ! number of nodes in an 1d element
  integer :: i,j,k
  real(kind=kdble),intent(in) :: xi ! point where to calculate lagrange function and
  !its derivative
  real(kind=kdble),dimension(nenode),intent(out) :: phi,dphi_dxi
  real(kind=kdble),dimension(nenode) :: xii,term,dterm,sum_term
  real(kind=kdble) :: dx

  ! compute natural coordnates
  dx = 2.0_kdble/real((nenode-1),kdble)! length = 2.0 as xi is taken -1 to +1
  do i = 1,nenode
    ! coordinates when origin is in the left
    xii(i) = real((i-1),kdble)*dx
  enddo

  ! origin is tranformed to mid point
  xii = xii-1.0_kdble

  do i = 1,nenode
    k = 0
    phi(i) = 1.0_kdble
    do j = 1,nenode
      if (j /= i) then
        k = k+1
        term(k) = (xi-xii(j))/(xii(i)-xii(j))
        dterm(k) = 1.0_kdble/(xii(i)-xii(j)) ! derivative of the term wrt xi

        phi(i) = phi(i)*(xi-xii(j))/(xii(i)-xii(j))
      endif
    enddo

    sum_term = 1.0_kdble
    do j = 1,nenode-1
      do k = 1,nenode-1
        if (k == j) then
          sum_term(j) = sum_term(j)*dterm(k)
        else
          sum_term(j) = sum_term(j)*term(k)
        endif
      enddo
    enddo
    dphi_dxi(i) = 0.0_kdble
    do j = 1,nenode-1
      dphi_dxi(i) = dphi_dxi(i) + sum_term(j)
    enddo
  enddo

  end subroutine lagrange1d

!
!===============================================================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

! not used ...

!  subroutine lagrange1dGEN(nenod,xi,phi,dphi_dxi)
!
!  implicit none
!  integer,intent(in) :: nenod ! number of nodes in an 1d element
!  integer :: i,j,k
!  real(kind=kdble),intent(in) :: xi ! point where to calculate lagrange function and
!  !its derivative
!  real(kind=kdble),dimension(nenod),intent(out) :: phi,dphi_dxi
!  real(kind=kdble),dimension(nenod) :: xii
!  real(kind=kdble),dimension(nenod-1) :: term,dterm,sum_term
!  real(kind=kdble) :: dx
!  real(kind=kdble),parameter :: one=1.0_kdble
!
!  ! compute natural coordnates
!  dx = 2.0_kdble/real((nenod-1),kdble)! length = 2.0 as xi is taken -1 to +1
!  do i = 1,nenod
!    ! coordinates when origin is in the left
!    xii(i)=real((i-1),kdble)*dx
!  enddo
!
!  ! origin is tranformed to mid point
!  xii = xii-one
!
!  do i = 1,nenod
!    k = 0
!    phi(i) = one
!    do j = 1,nenod
!      if (j /= i) then
!        k = k+1
!        term(k) = (xi-xii(j))/(xii(i)-xii(j))
!        dterm(k) = one/(xii(i)-xii(j)) ! derivative of the term wrt xi
!
!        phi(i) = phi(i)*term(k) !(xi-xii(j))/(xii(i)-xii(j))
!      endif
!    enddo
!
!    ! derivative of the polynomial: product rule
!    sum_term = one
!    do j = 1,nenod-1
!      do k = 1,nenod-1
!        if (k == j) then
!          sum_term(j) = sum_term(j)*dterm(k)
!        else
!          sum_term(j) = sum_term(j)*term(k)
!        endif
!      enddo
!    enddo
!    dphi_dxi(i) = sum(sum_term)
!    !dphi_dxi(i) = 0.0_kdble
!    !do j=1,nenod-1
!    !  dphi_dxi(i) = dphi_dxi(i)+sum_term(j)
!    !enddo
!  enddo
!
!  end subroutine lagrange1dGEN

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
! Assumed Shape array: pass pointer, subarray or allocatable array

  subroutine lagrange1dGENAS(nenod,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  integer :: i,j,k
  real(kind=kdble),intent(in) :: xi ! point where to calculate lagrange function and
  !its derivative
  real(kind=kdble),dimension(:),intent(out) :: phi,dphi_dxi !,dimension(nenod)
  real(kind=kdble),dimension(nenod) :: xii
  real(kind=kdble),dimension(nenod-1) :: term,dterm,sum_term
  real(kind=kdble) :: dx
  real(kind=kdble),parameter :: one=1.0_kdble

  ! compute natural coordnates
  dx = 2.0_kdble/real((nenod-1),kdble)! length = 2.0 as xi is taken -1 to +1
  do i = 1,nenod
    ! coordinates when origin is in the left
    xii(i)=real((i-1),kdble)*dx
  enddo

  ! origin is tranformed to mid point
  xii = xii-one

  do i = 1,nenod
    k = 0
    phi(i) = one
    do j = 1,nenod
      if (j /= i) then
        k = k+1
        term(k) = (xi-xii(j))/(xii(i)-xii(j))
        dterm(k) = one/(xii(i)-xii(j)) ! derivative of the term wrt xi

        phi(i) = phi(i)*term(k) !(xi-xii(j))/(xii(i)-xii(j))
      endif
    enddo

    ! derivative of the polynomial: product rule
    sum_term = one
    do j = 1,nenod-1
      do k = 1,nenod-1
        if (k == j) then
          sum_term(j) = sum_term(j)*dterm(k)
        else
          sum_term(j) = sum_term(j)*term(k)
        endif
      enddo
    enddo
    dphi_dxi(i) = sum(sum_term)
    !dphi_dxi(i) = 0.0_kdble
    !do j=1,nenod-1
    !  dphi_dxi(i) = dphi_dxi(i)+sum_term(j)
    !enddo
  enddo

  end subroutine lagrange1dGENAS

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

  subroutine lagrange1dGLL(nenod,xii,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  real(kind=kdble),dimension(nenod),intent(in) :: xii
  real(kind=kdble),intent(in) :: xi ! point where to calculate lagrange function and
  !its derivative
  real(kind=kdble),dimension(nenod),intent(out) :: phi,dphi_dxi

  integer :: i,j,k
  real(kind=kdble),dimension(nenod-1) :: term,dterm,sum_term
  real(kind=kdble),parameter :: one=1.0_kdble

  do i = 1,nenod
    k = 0
    phi(i) = one
    do j = 1,nenod
      if (j /= i) then
        k = k+1
        term(k) = (xi-xii(j))/(xii(i)-xii(j))
        dterm(k) = one/(xii(i)-xii(j)) ! derivative of the term wrt xi

        phi(i) = phi(i)*term(k)!(xi-xii(j))/(xii(i)-xii(j))
      endif
    enddo

    ! derivative of the polynomial: product rule
    sum_term = one
    do j = 1,nenod-1
      do k = 1,nenod-1
        if (k == j) then
          sum_term(j) = sum_term(j)*dterm(k)
        else
          sum_term(j) = sum_term(j)*term(k)
        endif
      enddo
    enddo
    dphi_dxi(i) = sum(sum_term)
    !dphi_dxi(i) = 0.0_kdble
    !do j=1,nenod-1
    !  dphi_dxi(i) = dphi_dxi(i)+sum_term(j)
    !enddo
  enddo

  end subroutine lagrange1dGLL

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
! Assumed Shape array: pass pointer, subarray or allocatable array

  subroutine lagrange1dGLLAS(nenod,xii,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  real(kind=kdble),dimension(nenod),intent(in) :: xii
  real(kind=kdble),intent(in) :: xi ! point where to calculate lagrange function and
  !its derivative
  real(kind=kdble),dimension(:),intent(out) :: phi,dphi_dxi !,dimension(nenod)

  integer :: i,j,k
  real(kind=kdble),dimension(nenod-1) :: term,dterm,sum_term
  real(kind=kdble),parameter :: one=1.0_kdble

  do i = 1,nenod
    k = 0
    phi(i) = one
    do j = 1,nenod
      if (j /= i) then
        k = k+1
        term(k) = (xi-xii(j))/(xii(i)-xii(j))
        dterm(k) = one/(xii(i)-xii(j)) ! derivative of the term wrt xi

        phi(i) = phi(i)*term(k)!(xi-xii(j))/(xii(i)-xii(j))
      endif
    enddo

    ! derivative of the polynomial: product rule
    sum_term = one
    do j = 1,nenod-1
      do k = 1,nenod-1
        if (k == j) then
          sum_term(j) = sum_term(j)*dterm(k)
        else
          sum_term(j) = sum_term(j)*term(k)
        endif
      enddo
    enddo
    dphi_dxi(i) = sum(sum_term)
    !dphi_dxi(i) = 0.0_kdble
    !do j=1,nenod-1
    !  dphi_dxi(i) = dphi_dxi(i)+sum_term(j)
    !enddo
  enddo

  end subroutine lagrange1dGLLAS

!===========================================
!
!  Library to compute the Gauss-Lobatto-Legendre points and weights
!  Based on Gauss-Lobatto routines from M.I.T.
!  Department of Mechanical Engineering
!
!===========================================

  real(kind=kdble) function endw1(n,alpha,beta) !double precision

  implicit none

  integer n
  real(kind=kdble) alpha,beta !double precision

  real(kind=kdble), parameter :: zero=0._kdble,one=1._kdble,two=2._kdble, &
  three = 3._kdble,four = 4._kdble !double precision
  real(kind=kdble) apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3 !double precision
  !double precision, external :: gammaf
  integer i

  f3 = zero
  apb   = alpha+beta
  if (n == 0) then
    endw1 = zero
    return
  endif
  f1   = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
    endw1 = f1
    return
  endif
  fint1 = gammaf(alpha+two)*gammaf(beta+one)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (-two*(beta+two)*fint1 + (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
    endw1 = f2
    return
  endif
  do i = 3,n
    di   = dble(i-1)
    abn  = alpha+beta+di
    abnn = abn+di
    a1   = -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
    a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
    a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
    f3   =  -(a2*f2+a1*f1)/a3
    f1   = f2
    f2   = f3
  enddo
  endw1  = f3

  end function endw1
  !=======================================================================

  real(kind=kdble) function endw2(n,alpha,beta) !double precision

  implicit none

  integer n
  real(kind=kdble) alpha,beta !double precision

  real(kind=kdble), parameter :: zero=0._kdble,one=1._kdble,two=2._kdble, &
  three = 3._kdble,four = 4._kdble !double precision
  real(kind=kdble) apb,f1,fint1,fint2,f2,di,abn,abnn,a1,a2,a3,f3 !double precision
  !real(kind=kdble), external :: gammaf
  integer i

  apb   = alpha+beta
  f3 = zero
  if (n == 0) then
    endw2 = zero
    return
  endif
  f1   = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  f1   = f1*(apb+two)*two**(apb+two)/two
  if (n == 1) then
    endw2 = f1
    return
  endif
  fint1 = gammaf(alpha+one)*gammaf(beta+two)/gammaf(apb+three)
  fint1 = fint1*two**(apb+two)
  fint2 = gammaf(alpha+two)*gammaf(beta+two)/gammaf(apb+four)
  fint2 = fint2*two**(apb+three)
  f2    = (two*(alpha+two)*fint1 - (apb+four)*fint2) * (apb+three)/four
  if (n == 2) then
    endw2 = f2
    return
  endif
  do i = 3,n
    di   = dble(i-1)
    abn  = alpha+beta+di
    abnn = abn+di
    a1   =  -(two*(di+alpha)*(di+beta))/(abn*abnn*(abnn+one))
    a2   =  (two*(alpha-beta))/(abnn*(abnn+two))
    a3   =  (two*(abn+one))/((abnn+two)*(abnn+one))
    f3   =  -(a2*f2+a1*f1)/a3
    f1   = f2
    f2   = f3
  enddo
  endw2  = f3

  end function endw2

  !
  !=======================================================================
  !

  real(kind=kdble) function gammaf (x) !double precision

  implicit none

  real(kind=kdble), parameter :: pi = 3.141592653589793_kdble !double precision

  real(kind=kdble) x !double precision

  real(kind=kdble), parameter :: half=0.5_kdble,one=1._kdble,two=2._kdble !double precision

  gammaf = one

  if (x == -half) gammaf = -two*sqrt(pi)
  if (x == half) gammaf =  sqrt(pi)
  if (x == one ) gammaf =  one
  if (x == two ) gammaf =  one
  if (x == 1.5_kdble) gammaf =  sqrt(pi)/2._kdble
  if (x == 2.5_kdble) gammaf =  1.5_kdble*sqrt(pi)/2._kdble
  if (x == 3.5_kdble) gammaf =  2.5_kdble*1.5_kdble*sqrt(pi)/2._kdble
  if (x == 3._kdble ) gammaf =  2._kdble
  if (x == 4._kdble ) gammaf = 6._kdble
  if (x == 5._kdble ) gammaf = 24._kdble
  if (x == 6._kdble ) gammaf = 120._kdble

  end function gammaf

  !
  !=====================================================================
  !

  subroutine jacg (xjac,np,alpha,beta)

  !=======================================================================
  !
  ! computes np Gauss points, which are the zeros of the
  ! Jacobi polynomial with parameters alpha and beta
  !
  !                  .alpha = beta =  0.0  ->  Legendre points
  !                  .alpha = beta = -0.5  ->  Chebyshev points
  !
  !=======================================================================

  implicit none

  integer np
  real(kind=kdble) alpha,beta !double precision
  real(kind=kdble) xjac(np) !double precision

  integer k,j,i,jmin,jm,n
  real(kind=kdble) xlast,dth,x,x1,x2,recsum,delx,xmin,swap !double precision
  real(kind=kdble) p,pd,pm1,pdm1,pm2,pdm2 !double precision

  integer, parameter :: K_MAX_ITER = 10
  real(kind=kdble), parameter :: zero = 0._kdble, eps = 1.0e-12_kdble !double precision

  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  xlast = 0._kdble
  n   = np-1
  dth = 4._kdble*atan(1._kdble)/(2._kdble*dble(n)+2._kdble)
  p = 0._kdble
  pd = 0._kdble
  jmin = 0
  do j = 1,np
    if (j == 1) then
      x = cos((2._kdble*(dble(j)-1._kdble)+1._kdble)*dth)
    else
      x1 = cos((2._kdble*(dble(j)-1._kdble)+1._kdble)*dth)
      x2 = xlast
      x  = (x1+x2)/2._kdble
    endif
    do k = 1,K_MAX_ITER
      call jacobf (p,pd,pm1,pdm1,pm2,pdm2,np,alpha,beta,x)
      recsum = 0._kdble
      jm = j-1
      do i = 1,jm
          recsum = recsum+1._kdble/(x-xjac(np-i+1))
      enddo
      delx = -p/(pd-recsum*p)
      x    = x+delx
      if (abs(delx) < eps) goto 31
    enddo
  31  continue
    xjac(np-j+1) = x
    xlast        = x
  enddo
  do i = 1,np
    xmin = 2._kdble
    do j = i,np
      if (xjac(j) < xmin) then
          xmin = xjac(j)
          jmin = j
      endif
    enddo
    if (jmin /= i) then
      swap = xjac(i)
      xjac(i) = xjac(jmin)
      xjac(jmin) = swap
    endif
  enddo

  end subroutine jacg

  !
  !=====================================================================
  !

  subroutine jacobf (poly,pder,polym1,pderm1,polym2,pderm2,n,alp,bet,x)

  !=======================================================================
  !
  ! Computes the Jacobi polynomial of degree n and its derivative at x
  !
  !=======================================================================

  implicit none

  real(kind=kdble) poly,pder,polym1,pderm1,polym2,pderm2,alp,bet,x !double precision
  integer n

  real(kind=kdble) apb,polyl,pderl,dk,a1,a2,b3,a3,a4,polyn,pdern,psave,pdsave !double precision
  integer k

  apb  = alp+bet
  poly = 1._kdble
  pder = 0._kdble
  psave = 0._kdble
  pdsave = 0._kdble

  if (n == 0) return

  polyl = poly
  pderl = pder
  poly  = (alp-bet+(apb+2._kdble)*x)/2._kdble
  pder  = (apb+2._kdble)/2._kdble
  if (n == 1) return

  do k = 2,n
    dk = dble(k)
    a1 = 2._kdble*dk*(dk+apb)*(2._kdble*dk+apb-2._kdble)
    a2 = (2._kdble*dk+apb-1._kdble)*(alp**2-bet**2)
    b3 = (2._kdble*dk+apb-2._kdble)
    a3 = b3*(b3+1._kdble)*(b3+2._kdble)
    a4 = 2._kdble*(dk+alp-1._kdble)*(dk+bet-1._kdble)*(2._kdble*dk+apb)
    polyn  = ((a2+a3*x)*poly-a4*polyl)/a1
    pdern  = ((a2+a3*x)*pder-a4*pderl+a3*poly)/a1
    psave  = polyl
    pdsave = pderl
    polyl  = poly
    poly   = polyn
    pderl  = pder
    pder   = pdern
  enddo

  polym1 = polyl
  pderm1 = pderl
  polym2 = psave
  pderm2 = pdsave

  end subroutine jacobf

  !
  !------------------------------------------------------------------------
  !

  real(kind=kdble) function PNDLEG (Z,N) !double precision

  !------------------------------------------------------------------------
  !
  !     Compute the derivative of the Nth order Legendre polynomial at Z.
  !     Based on the recursion formula for the Legendre polynomials.
  !
  !------------------------------------------------------------------------
  implicit none

  real(kind=kdble) z !double precision
  integer n

  real(kind=kdble) P1,P2,P1D,P2D,P3D,FK,P3 !double precision
  integer k

  P1   = 1._kdble
  P2   = Z
  P1D  = 0._kdble
  P2D  = 1._kdble
  P3D  = 1._kdble

  do K = 1, N-1
    FK  = dble(K)
    P3  = ((2._kdble*FK+1._kdble)*Z*P2-FK*P1)/(FK+1._kdble)
    P3D = ((2._kdble*FK+1._kdble)*P2+(2._kdble*FK+1._kdble)*Z*P2D-FK*P1D)/(FK+1._kdble)
    P1  = P2
    P2  = P3
    P1D = P2D
    P2D = P3D
  enddo

  PNDLEG = P3D

  end function pndleg

  !
  !------------------------------------------------------------------------
  !

  real(kind=kdble) function PNLEG (Z,N) !double precision

  !------------------------------------------------------------------------
  !
  !     Compute the value of the Nth order Legendre polynomial at Z.
  !     Based on the recursion formula for the Legendre polynomials.
  !
  !------------------------------------------------------------------------
  implicit none

  real(kind=kdble) z !double precision
  integer n

  real(kind=kdble) P1,P2,P3,FK !double precision
  integer k

  P1   = 1._kdble
  P2   = Z
  P3   = P2

  do K = 1, N-1
    FK  = dble(K)
    P3  = ((2._kdble*FK+1._kdble)*Z*P2 - FK*P1)/(FK+1._kdble)
    P1  = P2
    P2  = P3
  enddo

  PNLEG = P3

  end function pnleg

  !
  !------------------------------------------------------------------------
  !

  real(kind=kdble) function pnormj (n,alpha,beta) !double precision

  implicit none

  real(kind=kdble) alpha,beta !double precision
  integer n

  real(kind=kdble) one,two,dn,const,prod,dindx,frac !double precision
  !real(kind=kdble), external :: gammaf
  integer i

  one   = 1._kdble
  two   = 2._kdble
  dn    = dble(n)
  const = alpha+beta+one

  if (n <= 1) then
    prod   = gammaf(dn+alpha)*gammaf(dn+beta)
    prod   = prod/(gammaf(dn)*gammaf(dn+alpha+beta))
    pnormj = prod * two**const/(two*dn+const)
    return
  endif

  prod  = gammaf(alpha+one)*gammaf(beta+one)
  prod  = prod/(two*(one+const)*gammaf(const+one))
  prod  = prod*(one+alpha)*(two+alpha)
  prod  = prod*(one+beta)*(two+beta)

  do i = 3,n
    dindx = dble(i)
    frac  = (dindx+alpha)*(dindx+beta)/(dindx*(dindx+alpha+beta))
    prod  = prod*frac
  enddo

  pnormj = prod * two**const/(two*dn+const)

  end function pnormj

  !
  !------------------------------------------------------------------------
  !

  subroutine zwgjd(z,w,np,alpha,beta)

  !=======================================================================
  !
  !     Z w g j d : Generate np Gauss-Jacobi points and weights
  !                 associated with Jacobi polynomial of degree n = np-1
  !
  !     Note : Coefficients alpha and beta must be greater than -1.
  !     ----
  !=======================================================================

  implicit none

  real(kind=kdble), parameter :: zero=0._kdble,one=1._kdble,two=2._kdble !double precision

  integer np
  real(kind=kdble) z(np),w(np) !double precision
  real(kind=kdble) alpha,beta !double precision

  integer n,np1,np2,i
  real(kind=kdble) p,pd,pm1,pdm1,pm2,pdm2 !double precision
  real(kind=kdble) apb,dnp1,dnp2,fac1,fac2,fac3,fnorm,rcoef !double precision
  !real(kind=kdble), external :: gammaf,pnormj

  pd = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n    = np-1
  apb  = alpha+beta
  p    = zero
  pdm1 = zero

  if (np <= 0) then
    stop 'ERROR: number of Gauss points < 1!'
  endif

  if ((alpha <= -one) .or. (beta <= -one)) then
    stop 'ERROR: alpha and beta must be greater than -1!'
  endif

  if (np == 1) then
    z(1) = (beta-alpha)/(apb+two)
    w(1) = gammaf(alpha+one)*gammaf(beta+one)/gammaf(apb+two) * two**(apb+one)
    return
  endif

  call jacg(z,np,alpha,beta)

  np1   = n+1
  np2   = n+2
  dnp1  = dble(np1)
  dnp2  = dble(np2)
  fac1  = dnp1+alpha+beta+one
  fac2  = fac1+dnp1
  fac3  = fac2+one
  fnorm = pnormj(np1,alpha,beta)
  rcoef = (fnorm*fac2*fac3)/(two*fac1*dnp2)
  do i = 1,np
    call jacobf(p,pd,pm1,pdm1,pm2,pdm2,np2,alpha,beta,z(i))
    w(i) = -rcoef/(p*pdm1)
  enddo

  end subroutine zwgjd

  !
  !------------------------------------------------------------------------
  !

  subroutine zwgljd(z,w,np,alpha,beta)

  !=======================================================================
  !
  !     Z w g l j d : Generate np Gauss-Lobatto-Jacobi points and the
  !     -----------   weights associated with Jacobi polynomials of degree
  !                   n = np-1.
  !
  !     Note : alpha and beta coefficients must be greater than -1.
  !            Legendre polynomials are special case of Jacobi polynomials
  !            just by setting alpha and beta to 0.
  !
  !=======================================================================

  implicit none

  real(kind=kdble), parameter :: zero=0._kdble,one=1._kdble,two=2._kdble !double precision

  integer np
  real(kind=kdble) alpha,beta !double precision
  real(kind=kdble) z(np), w(np) !double precision

  integer n,nm1,i
  real(kind=kdble) p,pd,pm1,pdm1,pm2,pdm2 !double precision
  real(kind=kdble) alpg,betg !double precision
  !real(kind=kdble), external :: endw1,endw2

  p = zero
  pm1 = zero
  pm2 = zero
  pdm1 = zero
  pdm2 = zero

  n   = np-1
  nm1 = n-1
  pd  = zero

  if (np <= 1) then
    stop 'ERROR: number of Gauss-Lobatto points < 2!'
  endif

  ! with spectral elements, use at least 3 points
  if (np < 3) then
    stop 'WARNING: number of Gauss-Lobatto points < 3!'
  endif
  !if (np <= 2) stop 'minimum number of Gauss-Lobatto points for the SEM is 3'

  if ((alpha <= -one) .or. (beta <= -one)) then
    stop 'ERROR: alpha and beta must be greater than -1!'
  endif

  if (nm1 > 0) then
    alpg  = alpha+one
    betg  = beta+one
    call zwgjd(z(2),w(2),nm1,alpg,betg)
  endif

  z(1)  = - one
  z(np) =  one

  do i = 2,np-1
    w(i) = w(i)/(one-z(i)**2)
  enddo

  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(1))
  w(1)  = endw1(n,alpha,beta)/(two*pd)
  call jacobf(p,pd,pm1,pdm1,pm2,pdm2,n,alpha,beta,z(np))
  w(np) = endw2(n,alpha,beta)/(two*pd)

  end subroutine zwgljd


end module siem_gll_library

