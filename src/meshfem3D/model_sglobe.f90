!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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


!--------------------------------------------------------------------------------------------------
! SGLOBE-rani
!
! 3D global shear-wave isotropic and radially anisotropic model from a joint inversion for multiple data sets
!
! SGLOBE-rani is a global radially anisotropic shear wave speed model with radial anisotropy allowed in the whole mantle.
! It is based on a seismic data set of over 43M seismic surface wave (fundamental and overtones) and body wave measurements.
! It models simultaneously crustal thickness and mantle structure, and its reference model is PREM.
!
! NOTE: Kustowski et al., 2008 and Chang et al., 2014 showed that the anisotropic structure in the lowermost mantle
!       retrieved from global tomographic inversions can be strongly affected by leakage effects, so we discourage
!       interpreting SGLOBE-rani's anisotropic structure below ~1500 km depth.
!
! reference:
!   Chang, S.-J., A.M.G. Ferreira, J. Ritsema, H.J. van Heijst, and J.H. Woodhouse (2015),
!   Joint inversion for global isotropic and radially anisotropic mantle structure including crustal thickness perturbations,
!   J. Geophys. Res., 120, 4278-4300, https://doi.org/10.1002/2014JB011824.
!
! implementation:
! Elodie Kendall, 2018 - spherical harmonics model, up to degree 35
!                        (routines based on model_s40rts.f90 implementation)
!
!                        P-wave velocity perturbations (dvp) taken from P12 of S20RTS/S40RTS by default;
!                        density perturbations (drho) scaled from Vsv perturbations (dvsv);
!                        mantle model defined between Moho and CMB;
!                        uses PREM as 1D reference (also for attenuation & eta-parameter)
!
!--------------------------------------------------------------------------------------------------

  module model_sglobe_par

  ! three_d_mantle_model_constants
  integer, parameter :: NK_20 = 20  ! radial basis
  integer, parameter :: NS_35 = 35  ! horizontal basis

  ! model_sglobe_variables
  !a = positive m  (radial, theta, phi) --> (k,l,m) (maybe other way around??)
  !b = negative m  (radial, theta, phi) --> (k,l,-m)
  double precision,dimension(:,:,:),allocatable :: &
    SGLOBE_V_dvsv_a,SGLOBE_V_dvsv_b,SGLOBE_V_dvp_a,SGLOBE_V_dvp_b, &
    SGLOBE_V_dvsh_a,SGLOBE_V_dvsh_b

  ! splines
  double precision,dimension(:),allocatable :: SGLOBE_V_spknt
  double precision,dimension(:,:),allocatable :: SGLOBE_V_qq0
  double precision,dimension(:,:,:),allocatable :: SGLOBE_V_qq

  end module model_sglobe_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_sglobe_broadcast()

! standard routine to setup model
  use constants
  use model_sglobe_par

  implicit none

  ! local parameters
  integer :: ier

 ! allocates memory
  allocate(SGLOBE_V_dvsv_a(0:NK_20,0:NS_35,0:NS_35), &
           SGLOBE_V_dvsv_b(0:NK_20,0:NS_35,0:NS_35), &
           SGLOBE_V_dvsh_a(0:NK_20,0:NS_35,0:NS_35), &
           SGLOBE_V_dvsh_b(0:NK_20,0:NS_35,0:NS_35), &
           SGLOBE_V_dvp_a(0:NK_20,0:NS_35,0:NS_35), &
           SGLOBE_V_dvp_b(0:NK_20,0:NS_35,0:NS_35), &
           SGLOBE_V_spknt(NK_20+1), &
           SGLOBE_V_qq0(NK_20+1,NK_20+1), &
           SGLOBE_V_qq(3,NK_20+1,NK_20+1), &
           stat=ier)
  if ( ier /= 0 ) call exit_MPI(myrank,'error allocating SGLOBE_V arrays')

  ! the variables read are declared and stored in structure SGLOBE_V
  if (myrank == 0) call read_model_SGLOBE()

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(SGLOBE_V_dvsv_a,(NK_20+1)*(NS_35+1)*(NS_35+1))
  call bcast_all_dp(SGLOBE_V_dvsv_b,(NK_20+1)*(NS_35+1)*(NS_35+1))
  call bcast_all_dp(SGLOBE_V_dvsh_a,(NK_20+1)*(NS_35+1)*(NS_35+1))
  call bcast_all_dp(SGLOBE_V_dvsh_b,(NK_20+1)*(NS_35+1)*(NS_35+1))
  call bcast_all_dp(SGLOBE_V_dvp_a,(NK_20+1)*(NS_35+1)*(NS_35+1))
  call bcast_all_dp(SGLOBE_V_dvp_b,(NK_20+1)*(NS_35+1)*(NS_35+1))
  call bcast_all_dp(SGLOBE_V_spknt,NK_20+1)
  call bcast_all_dp(SGLOBE_V_qq0,(NK_20+1)*(NK_20+1))
  call bcast_all_dp(SGLOBE_V_qq,3*(NK_20+1)*(NK_20+1))

  end subroutine model_sglobe_broadcast

!-------------------------------------------------------------------------------------------------
!

  subroutine read_model_sglobe()

  use constants
  use model_sglobe_par

  implicit none

  ! local parameters
  integer :: k,l,m,ier
  character(len=*), parameter :: SGLOBEv = 'DATA/sglobe/dvsv.dat'
  character(len=*), parameter :: SGLOBEh = 'DATA/sglobe/dvsh.dat'
  character(len=*), parameter :: P12 = 'DATA/s20rts/P12.dat'

  ! SGLOBE degree 35 S model from Chang at al.
  ! dvsv
  open(unit=10,file=SGLOBEv,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) call exit_MPI(0,'error opening file SGLOBE.dat')

  SGLOBE_V_dvsv_a(:,:,:) = 0.d0
  SGLOBE_V_dvsv_b(:,:,:) = 0.d0
  do k = 0,NK_20
    do l = 0,NS_35
      read(10,*) SGLOBE_V_dvsv_a(k,l,0),(SGLOBE_V_dvsv_a(k,l,m),SGLOBE_V_dvsv_b(k,l,m),m=1,l)
    enddo
  enddo
  close(10)

  ! dvsh
  open(unit=10,file=SGLOBEh,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) call exit_MPI(0,'error opening file SGLOBE.dat')

  SGLOBE_V_dvsh_a(:,:,:) = 0.d0
  SGLOBE_V_dvsh_b(:,:,:) = 0.d0
  do k = 0,NK_20
    do l = 0,NS_35
      read(10,*) SGLOBE_V_dvsh_a(k,l,0),(SGLOBE_V_dvsh_a(k,l,m),SGLOBE_V_dvsh_b(k,l,m),m=1,l)
    enddo
  enddo
  close(10)

  ! P12 degree 12 P model from Ritsema
  open(unit=10,file=P12,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) call exit_MPI(0,'error opening file P12.dat')

  SGLOBE_V_dvp_a(:,:,:) = 0.d0
  SGLOBE_V_dvp_b(:,:,:) = 0.d0
  do k = 0,NK_20
    do l = 0,12
      read(10,*) SGLOBE_V_dvp_a(k,l,0),(SGLOBE_V_dvp_a(k,l,m),SGLOBE_V_dvp_b(k,l,m),m=1,l)
    enddo
    do l = 13,20
      SGLOBE_V_dvp_a(k,l,0) = 0.0d0
      do m = 1,l
        SGLOBE_V_dvp_a(k,l,m) = 0.0d0
        SGLOBE_V_dvp_b(k,l,m) = 0.0d0
      enddo
    enddo
  enddo
  close(10)

  ! set up the splines used as radial basis functions by Ritsema
  call sglobe_splhsetup()

  end subroutine read_model_sglobe

!---------------------------

  subroutine mantle_sglobe(radius,theta,phi,dvsv,dvsh,dvp,drho)

  use constants
  use model_sglobe_par

  implicit none

  double precision, intent(in) :: radius,theta,phi
  double precision, intent(out) :: dvsv,dvsh,dvp,drho

  ! local parameters

  !-----------------------------------------------------------------------------
  ! Density scaling:
  ! factor to convert perturbations in shear speed to perturbations in density
  double precision, parameter :: SCALE_RHO = 0.40d0

  ! Vp model scaling:
  ! factor to convert perturbations in vs speed to perturbations in vp
  logical :: USE_VP_SCALING = .false.             ! for .false., the P12 model will be taken for vp perturbations to PREM
  double precision, parameter :: SCALE_VP = 0.5d0 ! scaling factor used in the inversion by Chang et al.

  !-----------------------------------------------------------------------------

  double precision, parameter :: RMOHO_ = 6346600.d0
  double precision, parameter :: RCMB_ = 3480000.d0
  double precision, parameter :: R_EARTH_ = 6371000.d0

  integer :: l,m,k
  double precision :: r_moho,r_cmb,xr
  double precision :: dvsv_alm,dvsv_blm,dvsh_alm,dvsh_blm
  double precision :: dvp_alm,dvp_blm
  double precision :: sglobe_rsple,radial_basis(0:NK_20)
  double precision :: sint,cost,x(2*NS_35+1),dx(2*NS_35+1)
  double precision :: dvs

  dvsv = 0.d0
  dvsh = 0.d0
  dvp = 0.d0
  drho = 0.d0

  ! model defined between Moho and CMB
  r_moho = RMOHO_ / R_EARTH_
  r_cmb = RCMB_ / R_EARTH_
  if (radius >= r_moho .or. radius <= r_cmb) return

  xr = -1.0d0+2.0d0*(radius-r_cmb)/(r_moho-r_cmb)
  do k = 0,NK_20
    radial_basis(k) = sglobe_rsple(1,NK_20+1,SGLOBE_V_spknt(1),SGLOBE_V_qq0(1,NK_20+1-k),SGLOBE_V_qq(1,1,NK_20+1-k),xr)
  enddo

  do l = 0,NS_35
    sint = dsin(theta)
    cost = dcos(theta)
    call lgndr(l,cost,sint,x,dx)

    dvsv_alm = 0.0d0
    dvsh_alm = 0.0d0
    dvp_alm = 0.0d0
    do k = 0,NK_20
      dvsv_alm = dvsv_alm+radial_basis(k)*SGLOBE_V_dvsv_a(k,l,0)
      dvsh_alm = dvsh_alm+radial_basis(k)*SGLOBE_V_dvsh_a(k,l,0)
      dvp_alm = dvp_alm+radial_basis(k)*SGLOBE_V_dvp_a(k,l,0)
    enddo
    dvsv = dvsv+dvsv_alm*x(1)
    dvsh = dvsh+dvsh_alm*x(1)
    dvp = dvp+dvp_alm*x(1)

    do m = 1,l
      dvsv_alm = 0.0d0
      dvsh_alm = 0.0d0
      dvp_alm  = 0.0d0
      dvsv_blm = 0.0d0
      dvsh_blm = 0.0d0
      dvp_blm  = 0.0d0
      do k = 0,NK_20
        dvsv_alm = dvsv_alm + radial_basis(k)*SGLOBE_V_dvsv_a(k,l,m)
        dvsh_alm = dvsh_alm + radial_basis(k)*SGLOBE_V_dvsh_a(k,l,m)
        dvp_alm  = dvp_alm  + radial_basis(k)*SGLOBE_V_dvp_a(k,l,m)
        dvsv_blm = dvsv_blm + radial_basis(k)*SGLOBE_V_dvsv_b(k,l,m)
        dvsh_blm = dvsh_blm + radial_basis(k)*SGLOBE_V_dvsh_b(k,l,m)
        dvp_blm  = dvp_blm  + radial_basis(k)*SGLOBE_V_dvp_b(k,l,m)
      enddo
      dvsv = dvsv+(dvsv_alm*dcos(dble(m)*phi)+dvsv_blm*dsin(dble(m)*phi))*x(m+1)
      dvsh = dvsh+(dvsh_alm*dcos(dble(m)*phi)+dvsh_blm*dsin(dble(m)*phi))*x(m+1)
      dvp = dvp+(dvp_alm*dcos(dble(m)*phi)+dvp_blm*dsin(dble(m)*phi))*x(m+1)
    enddo
  enddo

  ! scales density perturbation from Vsv
  drho = SCALE_RHO*dvsv

  ! alternative Vp perturbation:
  ! instead of taking dvp from P12 (like S20RTS), uses a scaling from Vs as mentioned in Chang et al. (2015)
  if (USE_VP_SCALING) then
    dvs = sqrt( SCALE_VP * (dvsv**2 + dvsh**2) )
    dvp = 0.5 * dvs
  endif

  end subroutine mantle_sglobe

!----------------------------------

  subroutine sglobe_splhsetup()

  use constants
  use model_sglobe_par

  implicit none

  ! local parameters
  integer :: i,j
  double precision :: qqwk(3,NK_20+1)

  SGLOBE_V_spknt(1) = -1.00000d0
  SGLOBE_V_spknt(2) = -0.78631d0
  SGLOBE_V_spknt(3) = -0.59207d0
  SGLOBE_V_spknt(4) = -0.41550d0
  SGLOBE_V_spknt(5) = -0.25499d0
  SGLOBE_V_spknt(6) = -0.10909d0
  SGLOBE_V_spknt(7) = 0.02353d0
  SGLOBE_V_spknt(8) = 0.14409d0
  SGLOBE_V_spknt(9) = 0.25367d0
  SGLOBE_V_spknt(10) = 0.35329d0
  SGLOBE_V_spknt(11) = 0.44384d0
  SGLOBE_V_spknt(12) = 0.52615d0
  SGLOBE_V_spknt(13) = 0.60097d0
  SGLOBE_V_spknt(14) = 0.66899d0
  SGLOBE_V_spknt(15) = 0.73081d0
  SGLOBE_V_spknt(16) = 0.78701d0
  SGLOBE_V_spknt(17) = 0.83810d0
  SGLOBE_V_spknt(18) = 0.88454d0
  SGLOBE_V_spknt(19) = 0.92675d0
  SGLOBE_V_spknt(20) = 0.96512d0
  SGLOBE_V_spknt(21) = 1.00000d0

  do i = 1,NK_20+1
    do j = 1,NK_20+1
      if (i == j) then
        SGLOBE_V_qq0(j,i)=1.0d0
      else
        SGLOBE_V_qq0(j,i)=0.0d0
      endif
    enddo
  enddo
  do i = 1,NK_20+1
    call sglobe_rspln(1,NK_20+1,SGLOBE_V_spknt(1),SGLOBE_V_qq0(1,i),SGLOBE_V_qq(1,1,i),qqwk(1,1))
  enddo

  end subroutine sglobe_splhsetup

!----------------------------------

! changed the obsolescent f77 features in the two routines below
! now still awful Fortran, but at least conforms to f90 standard

  double precision function sglobe_rsple(I1,I2,X,Y,Q,S)

  implicit none

! rsple returns the value of the function y(x) evaluated at point S
! using the cubic spline coefficients computed by rspln and saved in Q.
! If S is outside the interval (x(i1),x(i2)) rsple extrapolates
! using the first or last interpolation polynomial. The arrays must
! be dimensioned at least - x(i2), y(i2), and q(3,i2).

  integer :: i1,i2
  double precision :: X(*),Y(*),Q(3,*),s

  ! local parameters
  integer :: i,ii
  double precision :: h

  i = 1
  II = I2-1

!   GUARANTEE I WITHIN BOUNDS.
  I = MAX0(I,I1)
  I = MIN0(I,II)

!   SEE IF X IS INCREASING OR DECREASING.
  if (X(I2)-X(I1) < 0) goto 1
  if (X(I2)-X(I1) >= 0) goto 2

!   X IS DECREASING.  CHANGE I AS NECESSARY.
1 continue
  if (S-X(I) <= 0) goto 3
  if (S-X(I) > 0) goto 4

4 continue
  I = I-1

  if (I-I1 < 0) goto 11
  if (I-I1 == 0) goto 6
  if (I-I1 > 0) goto 1

3 continue
  if (S-X(I+1) < 0) goto 5
  if (S-X(I+1) >= 0) goto 6

5 continue
  I = I+1

  if (I-II < 0) goto 3
  if (I-II == 0) goto 6
  if (I-II > 0) goto 7

!   X IS INCREASING.  CHANGE I AS NECESSARY.
2 continue
  if (S-X(I+1) <= 0) goto 8
  if (S-X(I+1) > 0) goto 9

9 continue
  I = I+1

  if (I-II < 0) goto 2
  if (I-II == 0) goto 6
  if (I-II > 0) goto 7

8 continue
  if (S-X(I) < 0) goto 10
  if (S-X(I) >= 0) goto 6

10 continue
  I = I-1
  if (I-I1 < 0) goto 11
  if (I-I1 == 0) goto 6
  if (I-I1 > 0) goto 8

7 continue
  I = II
  goto 6

11 continue
  I = I1

!   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
6 continue
  H = S-X(I)
  SGLOBE_RSPLE = Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))

  end function sglobe_rsple

!----------------------------------

  subroutine sglobe_rspln(I1,I2,X,Y,Q,F)

  implicit none

! The subroutine rspln computes cubic spline interpolation coefficients
! for y(x) between grid points i1 and i2 saving them in q.The
! interpolation is continuous with continuous first and second
! derivatives. It agrees exactly with y at grid points and with the
! three point first derivatives at both end points (i1 and i2).
! X must be monotonic but if two successive values of x are equal
! a discontinuity is assumed and separate interpolation is done on
! each strictly monotonic segment. The arrays must be dimensioned at
! least - x(i2), y(i2), q(3,i2), and f(3,i2).
! F is working storage for rspln.

  integer :: i1,i2
  double precision :: X(*),Y(*),Q(3,*),F(3,*)

  ! local parameters
  integer :: i,j,k,j1,j2
  double precision :: y0,a0,b0,b1,h,h2,ha,h2a,h3a,h2b
  double precision :: YY(3),SMALL

  equivalence (YY(1),Y0)
  data SMALL/1.0d-08/,YY/0.0d0,0.0d0,0.0d0/

  J1 = I1+1
  Y0 = 0.0d0

!   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL
  if (I2-I1 < 0) return
  if (I2-I1 == 0) goto 17
  if (I2-I1 > 0) goto 8

8 continue
  A0 = X(J1-1)

!   SEARCH FOR DISCONTINUITIES.
  do I=J1,I2
    B0 = A0
    A0 = X(I)
    if (DABS((A0-B0)/DMAX1(A0,B0)) < SMALL) goto 4
  enddo

17 continue
  J1 = J1-1
  J2 = I2-2
  goto 5

4 continue
  J1 = J1-1
  J2 = I-3

!   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
5 continue
  if (J2+1-J1 < 0) goto 9
  if (J2+1-J1 == 0) goto 10
  if (J2+1-J1 > 0) goto 11

!   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
10 continue
  J2 = J2+2
  Y0 =(Y(J2)-Y(J1))/(X(J2)-X(J1))
  do J=1,3
    Q(J,J1)=YY(J)
    Q(J,J2)=YY(J)
  enddo
  goto 12

!   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
11  continue
  A0 = 0.
  H = X(J1+1)-X(J1)
  H2 = X(J1+2)-X(J1)
  Y0 = H*H2*(H2-H)
  H = H*H
  H2 = H2*H2

!   CALCULATE DERIVITIVE AT NEAR END.
  B0 = (Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
  B1 = B0

!   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
  DO I=J1,J2
    H = X(I+1)-X(I)
    Y0 = Y(I+1)-Y(I)
    H2 = H*H
    HA = H-A0
    H2A = H-2.0d0*A0
    H3A = 2.0d0*H-3.0d0*A0
    H2B = H2*B0
    Q(1,I) = H2/HA
    Q(2,I) = -HA/(H2A*H2)
    Q(3,I) = -H*H2A/H3A
    F(1,I) = (Y0-H*B0)/(H*HA)
    F(2,I) = (H2B-Y0*(2.0d0*H-A0))/(H*H2*H2A)
    F(3,I) = -(H2B-3.0d0*Y0*HA)/(H*H3A)
    A0 = Q(3,I)
    B0 = F(3,I)
  enddo

!   TAKE CARE OF LAST TWO ROWS.
  I = J2+1
  H = X(I+1)-X(I)
  Y0 = Y(I+1)-Y(I)
  H2 = H*H
  HA = H-A0
  H2A = H*HA
  H2B = H2*B0-Y0*(2.0d0*H-A0)
  Q(1,I) = H2/HA
  F(1,I) = (Y0-H*B0)/H2A
  HA = X(J2)-X(I+1)
  Y0 = -H*HA*(HA+H)
  HA = HA*HA

!   CALCULATE DERIVATIVE AT FAR END.
  Y0 = (Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
  Q(3,I) = (Y0*H2A+H2B)/(H*H2*(H-2.0d0*A0))
  Q(2,I) = F(1,I)-Q(1,I)*Q(3,I)

!   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
  DO J=J1,J2
    K=I-1
    Q(1,I)=F(3,K)-Q(3,K)*Q(2,I)
    Q(3,K)=F(2,K)-Q(2,K)*Q(1,I)
    Q(2,K)=F(1,K)-Q(1,K)*Q(3,K)
    I=K
  enddo
  Q(1,I)=B1

!   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
9 continue
  J2 = J2+2
  DO J=1,3
    Q(J,J2)=YY(J)
  enddo

!   SEE IF THIS DISCONTINUITY IS THE LAST.
12 continue
  if (J2-I2 < 0) then
    goto 6
  else
    return
  endif

!   NO.  GO BACK FOR MORE.
6 continue
  J1 = J2+2
  if (J1-I2 <= 0) goto 8
  if (J1-I2 > 0) goto 7

!   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
7 continue
  do J=1,3
    Q(J,I2)=YY(J)
  enddo

  end subroutine sglobe_rspln

