program attenuation
!    Complile:
!      f90 -o attenuation_test atteuation.f90 attenuation_simplex.f90
!
!    Compiling this code with pgf90 truncates the relaxation times on output
!        to single precision (real(4)), whereas using ifort (ifc)  does not
!        For sanity, use ifort to compile and use the run_attenuation2.sh
!        script  to obtain the values for prem
!
! Brian Savage 19/01/05
!  This code should not produce the exact values that attenuation_prem.c
!    It has been updated to use a more robust inversion for the
!    the stress and strain relaxation times
!    A simplex inversion is used
!
! Brian Savage 22/03/04
!  This code should produce the exact values that attenuation_prem.c output
!  It is similar and the basis of the subroutines in attenuation_model.f90

  real(8) t1,t2,Q
  integer i, n

  integer write_central_period

  real(8) omega_not
  real(8), dimension(3) :: tau_e,tau_s

  write_central_period = 0

!  write(*,*)'longest period (seconds): '
  read(5,*, end=13)t1
!  write(*,*)'shortest period (seconds): '
  read(5,*, end=13)t2
  write(*,*)
  write(*,'(A16F11.6A4F11.6A2)')'! period range: ', t2, ' -- ', t1, ' s'
  write(*,*)
!  write(*,*)'number of mechanisms: '
  read(5,*, end=13)n
42 continue
!  write(*,*)'Q: '
  read(5,*, end=13)Q

  tau_e(:)  = 0.0d0
  tau_s(:)  = 0.0d0
  omega_not = 0.0d0

  !  call attenuation_liu(t1, t2, n, Q, omega_not, tau_s, tau_e)
  call attenuation_simplex(t1, t2, n, Q, omega_not, tau_s, tau_e)
  if(write_central_period == 0) then
     write(*,*)"! define central period of source in seconds using values from Jeroen's code"
     write(*,'(A27,F20.15,A2)')'  T_c_source = 1000.d0 / ', omega_not ,'d0'
     write(*,*)
     write(*,*)'! tau sigma evenly spaced in log frequency, does not depend on value of Q'
     do i = 1,n
        write(*,'(A13,I1,A4,F30.20,A2)')'  tau_sigma(',i,') = ', tau_s(i), 'd0'
     enddo
     write(*,*)
     write(*,*)"! check in which region we are based upon doubling flag"
     write(*,*)
     write(*,*)'  select case(iregion_attenuation)'
     write(*,*)
     write_central_period = 1

  endif
  if(Q == 84.6d0) then
     write(*,*)'!--- inner core, target Q_mu: 84.60'
     write(*,*)
     write(*,*)'  case(IREGION_ATTENUATION_INNER_CORE)'
     write(*,*)
  endif
  if(Q == 312.0d0) then
     write(*,*)'!--- CMB -> d670 (no attenuation in fluid outer core), target Q_mu = 312.'
     write(*,*)
     write(*,*)'  case(IREGION_ATTENUATION_CMB_670)'
     write(*,*)
  endif
  if(Q == 143.0d0) then
     write(*,*)'!--- d670 -> d220, target Q_mu: 143.'
     write(*,*)
     write(*,*)'  case(IREGION_ATTENUATION_670_220)'
     write(*,*)
  endif
  if(Q == 80.0d0) then
     write(*,*)'!--- d220 -> depth of 80 km, target Q_mu:  80.'
     write(*,*)
     write(*,*)'  case(IREGION_ATTENUATION_220_80)'
     write(*,*)
  endif
  if(Q == 600.0d0) then
     write(*,*)'!--- depth of 80 km -> surface, target Q_mu: 600.'
     write(*,*)
     write(*,*)'  case(IREGION_ATTENUATION_80_SURFACE)'
     write(*,*)
  endif

  do i = 1,n
     write(*,'(A12,I1,A4,F30.20,A2)')'     tau_mu(',i,') = ', tau_e(i), 'd0'
  enddo
  write(*,'(A17,F20.10,A2)')'       Q_mu = ', Q, 'd0'
  write(*,*)

  goto 42


13 continue
  write(*,*)'!--- do nothing for fluid outer core (no attenuation there)'
  write(*,*)
  write(*,*)'  case default'
  write(*,*)
  write(*,*)"    call exit_MPI(myrank,'wrong attenuation flag in mesh')"
  write(*,*)
  write(*,*)'  end select'
  write(*,*)

end program attenuation


subroutine attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)
  implicit none
  integer, parameter :: N_SLS = 3
  real(8), dimension(N_SLS) :: tau_s, alphaval, betaval,gammaval
  real(8)  deltat

  real(8), dimension(N_SLS) :: tauinv

  tauinv(:) = - 1.0 / tau_s(:)

  alphaval(:)  = 1 + deltat*tauinv(:) + deltat**2*tauinv(:)**2 / 2. + &
       deltat**3*tauinv(:)**3 / 6. + deltat**4*tauinv(:)**4 / 24.
  betaval(:)   = deltat / 2. + deltat**2*tauinv(:) / 3. + deltat**3*tauinv(:)**2 / 8. + deltat**4*tauinv(:)**3 / 24.
  gammaval(:)  = deltat / 2. + deltat**2*tauinv(:) / 6. + deltat**3*tauinv(:)**2 / 24.0

end subroutine attenuation_memory_values

subroutine attenuation_scale_factor(myrank, T_c_source, tau_mu, tau_sigma, Q_mu, scale_factor)
  implicit none
  integer, parameter :: N_SLS = 3

  integer myrank
  real(8) scale_factor, Q_mu, T_c_source
  real(8), dimension(N_SLS) :: tau_mu, tau_sigma

  real(8) scale_t, T_c_source_nondim
  real(8) f_c_source, w_c_source, f_0_prem
  real(8) factor_scale_mu0, factor_scale_mu
  real(8) a_val, b_val
  real(8) big_omega
  integer i

  real(8) PI, ZERO, ONE, TWO, GRAV, RHOAV, TWO_PI

  PI     = 3.14159265358979d0
  ZERO   = 0.0d0
  ONE    = 1.0d0
  TWO    = 2.0d0
  GRAV   = 6.6723d-11
  RHOAV  = 5514.3d0
  TWO_PI = 6.28318530717958

  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

  T_c_source_nondim = T_c_source / scale_t

!--- compute central angular frequency of source (non dimensionalized)
  f_c_source = ONE / T_c_source_nondim
  w_c_source = TWO_PI * f_c_source

!--- non dimensionalize PREM reference of 1 second
  f_0_prem = ONE / ( ONE / scale_t)

!--- quantity by which to scale mu_0 to get mu
  factor_scale_mu0 = ONE + TWO * log(f_c_source / f_0_prem) / (PI * Q_mu)

!--- compute a, b and Omega parameters, also compute one minus sum of betas
  a_val = ONE
  b_val = ZERO

  do i = 1,N_SLS
    a_val = a_val - w_c_source * w_c_source * tau_mu(i) * &
      (tau_mu(i) - tau_sigma(i)) / (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
    b_val = b_val + w_c_source * (tau_mu(i) - tau_sigma(i)) / &
      (1.d0 + w_c_source * w_c_source * tau_mu(i) * tau_mu(i))
  enddo

  big_omega = a_val*(sqrt(1.d0 + b_val*b_val/(a_val*a_val))-1.d0);

!--- quantity by which to scale mu to get mu_relaxed
  factor_scale_mu = b_val * b_val / (TWO * big_omega)

!--- total factor by which to scale mu0
  scale_factor = factor_scale_mu * factor_scale_mu0

!--- check that the correction factor is close to one
  if(scale_factor < 0.9 .or. scale_factor > 1.1) &
       call exit(-1)

end subroutine attenuation_scale_factor

subroutine attenuation_property_values(myrank, tau_s, tau_e, factor_common, one_minus_sum_beta)
  implicit none
  integer, parameter :: N_SLS = 3
  integer myrank
  real(8), dimension(N_SLS) :: tau_s, tau_e, beta, factor_common
  real(8)  one_minus_sum_beta

  real(8), dimension(N_SLS) :: tauinv
  integer i

  tauinv(:) = -1.0 / tau_s(:)

  beta(:) = 1.0d0 - tau_e(:) / tau_s(:)
  one_minus_sum_beta = 1.0d0

  do i = 1,N_SLS
     one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  factor_common(:) = 2.0 * beta(:) * tauinv(:)

end subroutine attenuation_property_values



!!!!!!
! Inversion by SVD

subroutine attenuation_invert_SVD(t1,t2,n,Q_real,omega_not,tau_s,tau_e)
  implicit none
  real(8)  t1, t2
  integer  n
  real(8)  Q_real
  real(8)  omega_not
  real(8), dimension(n)   :: tau_s, tau_e

  integer i, j, k

  real(8), dimension(n)   :: gradient, dadp, dbdp, dqdp, x1, x2
  real(8), dimension(n,n) :: hessian

  real(8) a, b, demon, Q_omega, Q_ratio, F_ratio, PI2, q
  real(8) f1, f2, exp1, exp2, expo, dexp, omega, df, d2qdp2

  gradient(:)  = 0.0d0
  hessian(:,:) = 0.0d0
  tau_e(:)     = 0.0d0
  tau_s(:)     = 0.0d0

  PI2 = 6.28318530717958d0
  write(*,*)'pi2: ', PI2

  f1 = 1.0d0/t1
  f2 = 1.0d0/t2

  if(f2 < f1 .OR. Q_real < 0.0d0 .OR. n < 1) then
     write(*,*)'bad parameters'
     call exit(-1)
  endif

  omega_not =  1.0e+03 * 10.0d0**(0.5 * (log10(f1) + log10(f2)))

  exp1 = log10(f1)
  exp2 = log10(f2)

  dexp = (exp2 - exp1) / real(n - 1.0d0)
  write(*,*)dexp
  q = 1.0d0 / (real(n - 1.0d0) * Q_real )
  write(*,*)'Q: ',q
  do i = 1,n
     expo     = exp1 + real(i-1) * dexp
     omega    = PI2 * 10.d0**expo
     write(*,*)omega,expo,PI2
     tau_s(i) = 1.0d0 / omega
     tau_e(i) = tau_s(i) * (1.0d0 + q) / (1.0d0 - q)
  enddo

  x1(:) = tau_e(:) - tau_s(:)
  x2(:) = tau_s(:)
  write(*,*)'e: ', tau_e
  write(*,*)'x1: ', x1
  write(*,*)'x2: ', x2

  exp1 = log10(f1);
  exp2 = log10(f2);
  dexp = (exp2 - exp1) / 100.0d0
  write(*,*)'exp1, exp2: ', exp1, exp2, dexp
  expo = exp1 - dexp
  do i=1,100
     expo = expo + dexp
     df       = 10.0d0**(expo+dexp) - 10.0d0**(expo)
     omega    = PI2 * 10.0d0**(expo)
     write(*,*) 'df,expo,omega: ',df,expo,omega
     a = real(1.0d0 - n)
     b = 0.0d0
     do j = 1,n
        tau_e(j) = x1(j) + x2(j)
        tau_s(j) = x2(j)
        demon   = 1.0d0 + omega**2 * tau_s(j)**2
        a       = a + (1.0d0 + omega**2.0d0 * tau_e(j) * tau_s(j)) / demon
        b       = b + ( omega * ( tau_e(j) - tau_s(j) ) ) / demon
        dadp(j) = omega**2.0d0 * tau_s(j) / demon
        dbdp(j) = omega / demon
     enddo
!     write(*,*)'ab: ', a, b,demon,expo
     Q_omega = a / b
     Q_ratio = 1.0d0 / Q_omega - 1.0d0 / Q_real
     F_ratio = df / (f2 - f1)
     do j = 1,n
        dqdp(j)     = (dbdp(j) - ( b / a ) * dadp(j)) / a
        gradient(j) = gradient(j) + 2.0d0 * (Q_ratio) * dqdp(j) * F_ratio
        do k = 1,j
           d2qdp2   = -(dadp(j) * dbdp(k) + dbdp(j) * dadp(k) - 2.0d0 * (b / a) * dadp(j) * dadp(k)) / (a * a)
           hessian(j,k) = hessian(j,k) + (2.0d0 * dqdp(j) * dqdp(k) + 2.0d0 * Q_ratio * d2qdp2) * F_ratio
           hessian(k,j) = hessian(j,k)
        enddo
     enddo
  enddo
  write(*,*)
  write(*,*)'blah: ', x1
  write(*,*)
  write(*,*)gradient
  write(*,*)
  write(*,*)hessian
  write(*,*)
  call invert(x1, gradient, hessian, n)

  tau_e(:) = x1(:) + x2(:)
  tau_s(:) = x2(:)

end subroutine attenuation_invert_SVD


subroutine invert(x,b,A,n)

  implicit none

  integer n
  real(8), dimension(n)   :: x, b
  real(8), dimension(n,n) :: A

  integer i, j, k
  real(8), dimension(n)   :: W, xp
  real(8), dimension(n,n) :: V
  real(8), dimension(n,n) :: A_inverse


  call svdcmp_dp(A,W,V,n)
  do i = 1,n
     do j = 1,n
        V(i,j) = (1.0d0 / W(i)) * A(j,i)
     enddo
  enddo

  do i = 1,n
     do j = 1,n
        A_inverse(i,j) = 0.0d0
        do k = 1,n
           A_inverse(i,j) = A_inverse(i,j) + A(i,k) * V(k,j)
        enddo
     enddo
  enddo

  do i = 1,n
     xp(i) = x(i)
     do j = 1, n
        xp(i) = xp(i) - A_inverse(i,j) * b(j)
     enddo
     x(i) = xp(i)
  enddo

end subroutine invert




FUNCTION pythag_dp(a,b)
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), INTENT(IN) :: a,b
  REAL(DP) :: pythag_dp
  REAL(DP) :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if (absa > absb) then
     pythag_dp=absa*sqrt(1.0d0+(absb/absa)**2)
  else
     if (absb == 0.0d0) then
        pythag_dp=0.0d0
     else
        pythag_dp=absb*sqrt(1.0d0+(absa/absb)**2)
     endif
  endif
END FUNCTION pythag_dp

SUBROUTINE svdcmp_dp(a,w,v,p)
!  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerprod
!  USE nr, ONLY : pythag
  IMPLICIT NONE
  integer p
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  REAL(DP), DIMENSION(p,p), INTENT(INOUT) :: a
  REAL(DP), DIMENSION(p), INTENT(OUT) :: w
  REAL(DP), DIMENSION(p,p), INTENT(OUT) :: v
  INTEGER(4) :: i,its,j,k,l,m,n,nm
  REAL(DP) :: anorm,c,f,g,h,s,scale,x,y,z
  REAL(DP), DIMENSION(size(a,1)) :: tempm
  REAL(DP), DIMENSION(size(a,2)) :: rv1,tempn
  REAL(DP)  PYTHAG_DP
  REAL(DP) OUTERPROD_D
  m=size(a,1)
!  n=assert_eq(size(a,2),size(v,1),size(v,2),size(w),'svdcmp_dp')
  n = size(a,2)
!  write(*,*)'Inside svdcmp_dp'
  g=0.0d0
  scale=0.0d0
  do i=1,n
     l=i+1
     rv1(i)=scale*g
     g=0.0d0
     scale=0.0d0
!     write(*,*)i,n,m,l
     if (i <= m) then
        scale=sum(abs(a(i:m,i)))
        if (scale /= 0.0d0) then
           a(i:m,i)=a(i:m,i)/scale
           s=dot_product(a(i:m,i),a(i:m,i))
           f=a(i,i)
           g=-sign(sqrt(s),f)
           h=f*g-s
           a(i,i)=f-g
           tempn(l:n)=matmul(a(i:m,i),a(i:m,l:n))/h
!           a(i:m,l:n)=a(i:m,l:n)+outerprod_d(a(i:m,i),m-1+1,tempn(l:n),n-l+1)
           a(i:m,l:n)=a(i:m,l:n)+spread(a(i:m,i),dim=2,ncopies=size(tempn(l:n))) * &
                spread(tempn(l:n),dim=1,ncopies=size(a(i:m,i)))
           a(i:m,i)=scale*a(i:m,i)
        endif
     endif
     w(i)=scale*g
     g=0.0d0
     scale=0.0d0
     if ((i <= m) .and. (i /= n)) then
        scale=sum(abs(a(i,l:n)))
        if (scale /= 0.0d0) then
           a(i,l:n)=a(i,l:n)/scale
           s=dot_product(a(i,l:n),a(i,l:n))
           f=a(i,l)
           g=-sign(sqrt(s),f)
           h=f*g-s
           a(i,l)=f-g
           rv1(l:n)=a(i,l:n)/h
           tempm(l:m)=matmul(a(l:m,l:n),a(i,l:n))
!           a(l:m,l:n)=a(l:m,l:n)+outerprod_d(tempm(l:m),m-l+1,rv1(l:n),n-l+1)
           a(l:m,l:n)=a(l:m,l:n)+spread(tempm(l:m),dim=2,ncopies=size(rv1(l:n))) * &
                spread(rv1(l:n),dim=1,ncopies=size(tempm(l:m)))
           a(i,l:n)=scale*a(i,l:n)
        endif
     endif
  enddo
  anorm=maxval(abs(w)+abs(rv1))
!  write(*,*)W
  do i=n,1,-1
     if (i < n) then
        if (g /= 0.0d0) then
           v(l:n,i)=(a(i,l:n)/a(i,l))/g
           tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
!           v(l:n,l:n)=v(l:n,l:n)+outerprod_d(v(l:n,i),n-1+1,tempn(l:n),n-l+1)
           v(l:n,l:n)=v(l:n,l:n)+spread(v(l:n,i),dim=2,ncopies=size(tempn(l:n))) * &
                spread(tempn(l:n), dim=1, ncopies=size(v(l:n,i)))
        endif
        v(i,l:n)=0.0d0
        v(l:n,i)=0.0d0
     endif
     v(i,i)=1.0d0
     g=rv1(i)
     l=i
  enddo
  do i=min(m,n),1,-1
     l=i+1
     g=w(i)
     a(i,l:n)=0.0d0
     if (g /= 0.0d0) then
        g=1.0d0/g
        tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
!        a(i:m,l:n)=a(i:m,l:n)+outerprod_d(a(i:m,i),m-i+1,tempn(l:n),n-l)
        a(i:m,l:n)=a(i:m,l:n)+spread(a(i:m,i),dim=2,ncopies=size(tempn(l:n))) * &
             spread(tempn(l:n),dim=1,ncopies=size(a(i:m,i)))
        a(i:m,i)=a(i:m,i)*g
     else
        a(i:m,i)=0.0d0
     endif
     a(i,i)=a(i,i)+1.0d0
  enddo
  do k=n,1,-1
     do its=1,30
        do l=k,1,-1
           nm=l-1
           if ((abs(rv1(l))+anorm) == anorm) exit
           if ((abs(w(nm))+anorm) == anorm) then
              c=0.0d0
              s=1.0d0
              do i=l,k
                 f=s*rv1(i)
                 rv1(i)=c*rv1(i)
                 if ((abs(f)+anorm) == anorm) exit
                 g=w(i)
                 h=pythag_dp(f,g)
                 w(i)=h
                 h=1.0d0/h
                 c= (g*h)
                 s=-(f*h)
                 tempm(1:m)=a(1:m,nm)
                 a(1:m,nm)=a(1:m,nm)*c+a(1:m,i)*s
                 a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
              enddo
              exit
           endif
        enddo
        z=w(k)
        if (l == k) then
           if (z < 0.0d0) then
              w(k)=-z
              v(1:n,k)=-v(1:n,k)
           endif
           exit
        endif
        if (its == 30) then
           write(*,*) 'svdcmp_dp: no convergence in svdcmp'
           call exit(-1)
        endif
        x=w(l)
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k)
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
        g=pythag_dp(f,1.0d0)
        f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
        c=1.0d0
        s=1.0d0
        do j=l,nm
           i=j+1
           g=rv1(i)
           y=w(i)
           h=s*g
           g=c*g
           z=pythag_dp(f,h)
           rv1(j)=z
           c=f/z
           s=h/z
           f= (x*c)+(g*s)
           g=-(x*s)+(g*c)
           h=y*s
           y=y*c
           tempn(1:n)=v(1:n,j)
           v(1:n,j)=v(1:n,j)*c+v(1:n,i)*s
           v(1:n,i)=-tempn(1:n)*s+v(1:n,i)*c
           z=pythag_dp(f,h)
           w(j)=z
           if (z /= 0.0d0) then
              z=1.0d0/z
              c=f*z
              s=h*z
           endif
           f= (c*g)+(s*y)
           x=-(s*g)+(c*y)
           tempm(1:m)=a(1:m,j)
           a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
           a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
        enddo
        rv1(l)=0.0d0
        rv1(k)=f
        w(k)=x
     enddo
  enddo
END SUBROUTINE svdcmp_dp

