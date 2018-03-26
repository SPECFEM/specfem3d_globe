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


  subroutine ylm(xlat,xlon,LMAX,Y)

! spherical harmonics Ylm evaluated at (lat,lon) up to degree LMAX using Legendre polynomials

  implicit none

  real(kind=4),intent(in) :: xlat,xlon

  integer,intent(in) :: LMAX
  real(kind=4), dimension((LMAX+1)**2),intent(inout) :: Y

  ! local parameters
!
!     WK1 SHOULD BE DIMENSIONED AT LEAST (LMAX+1)*4
!
  real(kind=4), dimension(LMAX+1) :: wk1
  integer :: IM,IL1,IND,LM1,L

  real(kind=4) :: theta,phi
  complex :: TEMP,FAC,DFAC

  real(kind=4), parameter :: RADIAN = 57.29577951308232 ! 180.0/PI

  ! converts to lat/lon (in degrees) to colat/lon (in rad)
  theta = (90.0 - xlat)/RADIAN
  phi = xlon/RADIAN

  IND = 0
  LM1 = LMAX + 1

  do IL1 = 1,LM1
    ! index L goes from 0 to LMAX
    L = IL1 - 1

    ! legendre polynomials
    call legndr(theta,L,L,wk1)

    FAC = (1.0,0.0)
    DFAC = CEXP(CMPLX(0.0,phi))

    ! loops over M
    do IM = 1,IL1
      ! index IM goes maximum from 1 to LMAX+1
      TEMP = FAC * CMPLX(wk1(IM),0.0)
      IND = IND + 1
      Y(IND) = REAL(TEMP)
      if (IM /= 1) then
        IND = IND + 1
        Y(IND) = AIMAG(TEMP)
      endif
      FAC = FAC * DFAC
    enddo
  enddo

  end subroutine ylm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine legndr(theta,L,M,X)

! Legendre polynomials and associated Legendre functions Plm(THETA) for degree l and m

  implicit none

  real(kind=4),intent(in) :: theta

  integer, intent(in) :: L,M
  real(kind=4),intent(inout) :: X(M+1)

  ! local parameters
  real(kind=4) :: XP(M+1),XCOSEC(M+1)
  real(kind=4) :: DSFL3,COSEC,SFL3
  double precision :: CT,ST,FCT,COT,X1,X2,X3,F1,F2,XM,TH
  double precision :: SMALL,sumval,COMPAR
  integer :: i,k,MP1,LP1

  double precision, parameter :: FPI = 12.56637062d0

!!!!!! illegal statement, removed by Dimitri Komatitsch   DFLOAT(I)=FLOAT(I)

  TH = theta
  CT = DCOS(TH)
  ST = DSIN(TH)
  LP1 = L+1
  MP1 = M+1

  FCT = DSQRT(dble(2*L+1)/FPI)
  SFL3 = SQRT(FLOAT(L*(L+1)))
  DSFL3 = SFL3

  COMPAR = dble(2*L+1)/FPI
  SMALL = 1.D-16 * COMPAR

  ! initializes work arrays
  do i = 1,MP1
    X(i) = 0.0
    XCOSEC(i) = 0.0
    XP(i) = 0.0
  enddo

  if (L > 1 .and. ABS(THETA) > 1.E-5) goto 3

  X(1) = FCT
  if (L == 0) return

  X(1) = CT * FCT
  X(2) = -ST * FCT/DSFL3
  XP(1) = -ST * FCT
  XP(2) = -0.5d0 * CT * FCT * DSFL3
  if (ABS(THETA) < 1.E-5) XCOSEC(2) = XP(2)
  if (ABS(THETA) >= 1.E-5) XCOSEC(2) = X(2)/ST

  return

 3 continue

  X1 = 1.d0
  X2 = CT

  do i = 2,L
    X3 = (dble(2*i-1)*CT*X2-dble(i-1)*X1)/dble(i)
    X1 = X2
    X2 = X3
  enddo

  COT = CT/ST
  COSEC = 1.0/ST
  X3 = X2*FCT
  X2 = dble(L)*(X1-CT*X2)*FCT/ST
  X(1) = X3
  X(2) = X2
  sumval = X3*X3

  XP(1) = -X2
  XP(2) = dble(L*(L+1))*X3-COT*X2
  X(2) = -X(2)/SFL3
  XCOSEC(2) = X(2)*COSEC
  XP(2) = -XP(2)/SFL3

  sumval = sumval + 2.d0*X(2)*X(2)
  if (sumval - COMPAR > SMALL) return

  X1 = X3
  X2 = -X2/DSQRT(dble(L*(L+1)))

  do i = 3,MP1
    K = i - 1
    F1 = DSQRT(dble((L+i-1)*(L-i+2)))
    F2 = DSQRT(dble((L+i-2)*(L-i+3)))
    XM = K
    X3 = -(2.d0*COT*(XM-1.D0)*X2+F2*X1)/F1

    sumval = sumval + 2.d0*X3*X3
    if (sumval - COMPAR > SMALL .and. i /= LP1) return

    X(i) = X3
    XCOSEC(i) = X(i)*COSEC
    X1 = X2
    XP(i) = -(F1*X2+XM*COT*X3)
    X2 = X3
  enddo

  end subroutine legndr

