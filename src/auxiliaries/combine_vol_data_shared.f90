!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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
! the Free Software Foundation; either version 2 of the License, or
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

!
! ------------------------------------------------------------------------------------------------
!
  subroutine reverse_ellipticity(x,y,z,nspl,rspl,espl,espl2)

  use constants

  implicit none

  real(kind=CUSTOM_REAL) :: x,y,z
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)
  double precision x1,y1,z1

  double precision ell
  double precision r,theta,phi,factor
  double precision cost,p20

  ! gets spherical coordinates
  x1 = x
  y1 = y
  z1 = z
  call xyz_2_rthetaphi_dble(x1,y1,z1,r,theta,phi)

  cost=dcos(theta)
  p20=0.5d0*(3.0d0*cost*cost-1.0d0)

  ! get ellipticity using spline evaluation
  call spline_evaluation(rspl,espl,espl2,nspl,r,ell)

  factor=ONE-(TWO/3.0d0)*ell*p20

  ! removes ellipticity factor
  x = x / factor
  y = y / factor
  z = z / factor

  end subroutine reverse_ellipticity

!
! ------------------------------------------------------------------------------------------------
!

! copy from make_ellipticity.f90 to avoid compiling issues

  subroutine make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

! creates a spline for the ellipticity profile in PREM
! radius and density are non-dimensional

  use constants

  implicit none

  integer nspl

  logical ONE_CRUST

! radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_ELLIPTICITY = 6371000.d0
! radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_ELLIPTICITY = 6368000.d0

  double precision rspl(NR),espl(NR),espl2(NR)

  integer i
  double precision ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                   R771,RTOPDDOUBLEPRIME,RCMB,RICB
  double precision r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision r(NR),rho(NR),epsilonval(NR),eta(NR)
  double precision radau(NR),z,k(NR),g_a,bom,exponentval,i_rho,i_radau
  double precision s1(NR),s2(NR),s3(NR)
  double precision yp1,ypn

! PREM
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0
  R80  = 6291000.d0
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

! non-dimensionalize
  r_icb = RICB/R_EARTH_ELLIPTICITY
  r_cmb = RCMB/R_EARTH_ELLIPTICITY
  r_topddoubleprime = RTOPDDOUBLEPRIME/R_EARTH_ELLIPTICITY
  r_771 = R771/R_EARTH_ELLIPTICITY
  r_670 = R670/R_EARTH_ELLIPTICITY
  r_600 = R600/R_EARTH_ELLIPTICITY
  r_400 = R400/R_EARTH_ELLIPTICITY
  r_220 = R220/R_EARTH_ELLIPTICITY
  r_80 = R80/R_EARTH_ELLIPTICITY
  r_moho = RMOHO/R_EARTH_ELLIPTICITY
  r_middle_crust = RMIDDLE_CRUST/R_EARTH_ELLIPTICITY
  r_ocean = ROCEAN_ELLIPTICITY/R_EARTH_ELLIPTICITY
  r_0 = 1.d0

  do i=1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  do i=164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  do i=324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  do i=337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  do i=518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  do i=531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  do i=541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  do i=566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  do i=591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  do i=610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  do i=620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  do i=627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  do i=634,NR
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
  do i=1,NR
    call prem_density(r(i),rho(i),ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
    radau(i)=rho(i)*r(i)*r(i)
  enddo

  eta(1)=0.0d0

  k(1)=0.0d0

  do i=2,NR
    call intgrl(i_rho,r,1,i,rho,s1,s2,s3)
    call intgrl(i_radau,r,1,i,radau,s1,s2,s3)
    z=(2.0d0/3.0d0)*i_radau/(i_rho*r(i)*r(i))
    eta(i)=(25.0d0/4.0d0)*((1.0d0-(3.0d0/2.0d0)*z)**2.0d0)-1.0d0
    k(i)=eta(i)/(r(i)**3.0d0)
  enddo

  g_a=4.0D0*i_rho

  bom=TWO_PI/(HOURS_PER_DAY*SECONDS_PER_HOUR)

! non-dimensionalized version
  bom=bom/sqrt(PI*GRAV*RHOAV)

!! DK DK I think 24.d0 below stands for HOURS_PER_DAY and thus I replace it here
!! DK DK in order to be consistent if someone uses the code one day for other planets
!! DK DK with a different rotation rate, or for the Earth in the past of in the future i.e. with a different rate as well.
!! DK DK Please do not hesitate to fix it back if my assumption below was wrong.
!! DK DK  epsilonval(NR)=15.0d0*(bom**2.0d0)/(24.0d0*i_rho*(eta(NR)+2.0d0))
  epsilonval(NR)=15.0d0*(bom**2.0d0)/(HOURS_PER_DAY*i_rho*(eta(NR)+2.0d0))

  do i=1,NR-1
    call intgrl(exponentval,r,i,NR,k,s1,s2,s3)
    epsilonval(i)=epsilonval(NR)*exp(-exponentval)
  enddo

! get ready to spline epsilonval
  nspl=1
  rspl(1)=r(1)
  espl(1)=epsilonval(1)
  do i=2,NR
    if(r(i) /= r(i-1)) then
      nspl=nspl+1
      rspl(nspl)=r(i)
      espl(nspl)=epsilonval(i)
    endif
  enddo

! spline epsilonval
  yp1=0.0d0
  ypn=(5.0d0/2.0d0)*(bom**2)/g_a-2.0d0*epsilonval(NR)
  call spline_construction(rspl,espl,nspl,yp1,ypn,espl2)

  end subroutine make_ellipticity

!
! ------------------------------------------------------------------------------------------------
!

! copy from model_prem.f90 to avoid compiling issues

  subroutine prem_density(x,rho,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

  use constants

  implicit none

  double precision x,rho,RICB,RCMB,RTOPDDOUBLEPRIME, &
      R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN

  logical ONE_CRUST

  double precision r

  ! compute real physical radius in meters
  r = x * R_EARTH

  ! calculates density according to radius
  if(r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
  else if(r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
  else if(r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
  else if(r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
  else if(r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
  else if(r > R220 .and. r <= R80) then
    rho=2.6910d0+0.6924d0*x
  else
    if(r > R80 .and. r <= RMOHO) then
      rho=2.6910d0+0.6924d0*x
    else if(r > RMOHO .and. r <= RMIDDLE_CRUST) then
      if(ONE_CRUST) then
        rho=2.6d0
      else
        rho=2.9d0
      endif
    else if(r > RMIDDLE_CRUST .and. r <= ROCEAN) then
      rho=2.6d0
    else if(r > ROCEAN) then
      rho=2.6d0
    endif
  endif

  rho=rho*1000.0d0/RHOAV

  end subroutine prem_density

!
! ------------------------------------------------------------------------------------------------
!

! copy from intgrl.f90 to avoid compiling issues


 subroutine intgrl(sumval,r,nir,ner,f,s1,s2,s3)

! Computes the integral of f[i]*r[i]*r[i] from i=nir to i=ner for
! radii values as in model PREM_an640

  implicit none

! Argument variables
  integer ner,nir
  double precision f(640),r(640),s1(640),s2(640)
  double precision s3(640),sumval

! Local variables
  double precision, parameter :: third = 1.0d0/3.0d0
  double precision, parameter :: fifth = 1.0d0/5.0d0
  double precision, parameter :: sixth = 1.0d0/6.0d0

  double precision rji,yprime(640)
  double precision s1l,s2l,s3l

  integer i,j,n,kdis(28)
  integer ndis,nir1

  data kdis/163,323,336,517,530,540,565,590,609,619,626,633,16*0/

  ndis = 12
  n = 640

  call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)
  nir1 = nir + 1
  sumval = 0.0d0
  do i=nir1,ner
    j = i-1
    rji = r(i) - r(j)
    s1l = s1(j)
    s2l = s2(j)
    s3l = s3(j)
    sumval = sumval + r(j)*r(j)*rji*(f(j) &
              + rji*(0.5d0*s1l + rji*(third*s2l + rji*0.25d0*s3l))) &
              + 2.0d0*r(j)*rji*rji*(0.5d0*f(j) + rji*(third*s1l + rji*(0.25d0*s2l + rji*fifth*s3l))) &
              + rji*rji*rji*(third*f(j) + rji*(0.25d0*s1l + rji*(fifth*s2l + rji*sixth*s3l)))
  enddo

  end subroutine intgrl

! -------------------------------

  subroutine deriv(y,yprime,n,r,ndis,kdis,s1,s2,s3)

  implicit none

! Argument variables
  integer kdis(28),n,ndis
  double precision r(n),s1(n),s2(n),s3(n)
  double precision y(n),yprime(n)

! Local variables
  integer i,j,j1,j2
  integer k,nd,ndp
  double precision a0,b0,b1
  double precision f(3,1000),h,h2,h2a
  double precision h2b,h3a,ha,s13
  double precision s21,s32,yy(3)

  yy(1) = 0.d0
  yy(2) = 0.d0
  yy(3) = 0.d0

  ndp=ndis+1
  do 3 nd=1,ndp
  if(nd == 1) goto 4
  if(nd == ndp) goto 5
  j1=kdis(nd-1)+1
  j2=kdis(nd)-2
  goto 6
    4 j1=1
  j2=kdis(1)-2
  goto 6
    5 j1=kdis(ndis)+1
  j2=n-2
    6 if((j2+1-j1)>0) goto 11
  j2=j2+2
  yy(1)=(y(j2)-y(j1))/(r(j2)-r(j1))
  s1(j1)=yy(1)
  s1(j2)=yy(1)
  s2(j1)=yy(2)
  s2(j2)=yy(2)
  s3(j1)=yy(3)
  s3(j2)=yy(3)
  goto 3
   11 a0=0.0d0
  if(j1 == 1) goto 7
  h=r(j1+1)-r(j1)
  h2=r(j1+2)-r(j1)
  yy(1)=h*h2*(h2-h)
  h=h*h
  h2=h2*h2
  b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/yy(1)
  goto 8
 7 b0=0.0d0
 8 b1=b0

  if(j2 > 1000) stop 'error in subroutine deriv for j2'

  do i=j1,j2
    h=r(i+1)-r(i)
    yy(1)=y(i+1)-y(i)
    h2=h*h
    ha=h-a0
    h2a=h-2.0d0*a0
    h3a=2.0d0*h-3.0d0*a0
    h2b=h2*b0
    s1(i)=h2/ha
    s2(i)=-ha/(h2a*h2)
    s3(i)=-h*h2a/h3a
    f(1,i)=(yy(1)-h*b0)/(h*ha)
    f(2,i)=(h2b-yy(1)*(2.0d0*h-a0))/(h*h2*h2a)
    f(3,i)=-(h2b-3.0d0*yy(1)*ha)/(h*h3a)
    a0=s3(i)
    b0=f(3,i)
  enddo

  i=j2+1
  h=r(i+1)-r(i)
  yy(1)=y(i+1)-y(i)
  h2=h*h
  ha=h-a0
  h2a=h*ha
  h2b=h2*b0-yy(1)*(2.d0*h-a0)
  s1(i)=h2/ha
  f(1,i)=(yy(1)-h*b0)/h2a
  ha=r(j2)-r(i+1)
  yy(1)=-h*ha*(ha+h)
  ha=ha*ha
  yy(1)=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/yy(1)
  s3(i)=(yy(1)*h2a+h2b)/(h*h2*(h-2.0d0*a0))
  s13=s1(i)*s3(i)
  s2(i)=f(1,i)-s13

  do j=j1,j2
    k=i-1
    s32=s3(k)*s2(i)
    s1(i)=f(3,k)-s32
    s21=s2(k)*s1(i)
    s3(k)=f(2,k)-s21
    s13=s1(k)*s3(k)
    s2(k)=f(1,k)-s13
    i=k
  enddo

  s1(i)=b1
  j2=j2+2
  s1(j2)=yy(1)
  s2(j2)=yy(2)
  s3(j2)=yy(3)
 3 continue

  do i=1,n
    yprime(i)=s1(i)
  enddo

  end subroutine deriv

!
! ------------------------------------------------------------------------------------------------
!

! copy from spline_routines.f90 to avoid compiling issues

! compute spline coefficients

  subroutine spline_construction(xpoint,ypoint,npoint,tangent_first_point,tangent_last_point,spline_coefficients)

  implicit none

! tangent to the spline imposed at the first and last points
  double precision, intent(in) :: tangent_first_point,tangent_last_point

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients output by the routine
  double precision, dimension(npoint), intent(out) :: spline_coefficients

  integer :: i

  double precision, dimension(:), allocatable :: temporary_array

  allocate(temporary_array(npoint))

  spline_coefficients(1) = - 1.d0 / 2.d0

  temporary_array(1) = (3.d0/(xpoint(2)-xpoint(1)))*((ypoint(2)-ypoint(1))/(xpoint(2)-xpoint(1))-tangent_first_point)

  do i = 2,npoint-1

    spline_coefficients(i) = ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))-1.d0) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

    temporary_array(i) = (6.d0*((ypoint(i+1)-ypoint(i))/(xpoint(i+1)-xpoint(i)) &
       - (ypoint(i)-ypoint(i-1))/(xpoint(i)-xpoint(i-1)))/(xpoint(i+1)-xpoint(i-1)) &
       - (xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*temporary_array(i-1)) &
       / ((xpoint(i)-xpoint(i-1))/(xpoint(i+1)-xpoint(i-1))*spline_coefficients(i-1)+2.d0)

  enddo

  spline_coefficients(npoint) = ((3.d0/(xpoint(npoint)-xpoint(npoint-1))) &
      * (tangent_last_point-(ypoint(npoint)-ypoint(npoint-1))/(xpoint(npoint)-xpoint(npoint-1))) &
      - 1.d0/2.d0*temporary_array(npoint-1))/(1.d0/2.d0*spline_coefficients(npoint-1)+1.d0)

  do i = npoint-1,1,-1
    spline_coefficients(i) = spline_coefficients(i)*spline_coefficients(i+1) + temporary_array(i)
  enddo

  deallocate(temporary_array)

  end subroutine spline_construction

! --------------

! evaluate a spline

  subroutine spline_evaluation(xpoint,ypoint,spline_coefficients,npoint,x_evaluate_spline,y_spline_obtained)

  implicit none

! number of input points and coordinates of the input points
  integer, intent(in) :: npoint
  double precision, dimension(npoint), intent(in) :: xpoint,ypoint

! spline coefficients to use
  double precision, dimension(npoint), intent(in) :: spline_coefficients

! abscissa at which we need to evaluate the value of the spline
  double precision, intent(in):: x_evaluate_spline

! ordinate evaluated by the routine for the spline at this abscissa
  double precision, intent(out):: y_spline_obtained

  integer :: index_loop,index_lower,index_higher

  double precision :: coef1,coef2

! initialize to the whole interval
  index_lower = 1
  index_higher = npoint

! determine the right interval to use, by dichotomy
  do while (index_higher - index_lower > 1)
! compute the middle of the interval
    index_loop = (index_higher + index_lower) / 2
    if(xpoint(index_loop) > x_evaluate_spline) then
      index_higher = index_loop
    else
      index_lower = index_loop
    endif
  enddo

! test that the interval obtained does not have a size of zero
! (this could happen for instance in the case of duplicates in the input list of points)
  if(xpoint(index_higher) == xpoint(index_lower)) stop 'incorrect interval found in spline evaluation'

  coef1 = (xpoint(index_higher) - x_evaluate_spline) / (xpoint(index_higher) - xpoint(index_lower))
  coef2 = (x_evaluate_spline - xpoint(index_lower)) / (xpoint(index_higher) - xpoint(index_lower))

  y_spline_obtained = coef1*ypoint(index_lower) + coef2*ypoint(index_higher) + &
        ((coef1**3 - coef1)*spline_coefficients(index_lower) + &
         (coef2**3 - coef2)*spline_coefficients(index_higher))*((xpoint(index_higher) - xpoint(index_lower))**2)/6.d0

  end subroutine spline_evaluation

