
module attenuation_model_constants
  implicit none
  include 'constants.h'

end module attenuation_model_constants

module attenuation_model_variables
  use attenuation_model_constants
  implicit none
  double precision min_period, max_period

end module attenuation_model_variables

subroutine attenuation_model(myrank, xlat, xlon, x, Qmu, tau_s, tau_e, T_c_source)

!! DK DK
!! DK DK xlat, xlon currently not used in this routine (which uses PREM).
!! DK DK The user needs to modify this routine if he wants to use
!! DK DK a particular 3D attenuation model. The current version is 1D.
!! DK DK

  use attenuation_model_variables
  implicit none

  integer myrank
  double precision xlat, xlon, r, x, Qmu
  double precision Qkappa, T_c_source
  double precision, dimension(N_SLS) :: tau_s, tau_e

  r = x * R_EARTH

! PREM
!
!--- inner core
!
  if(r >= 0.d0 .and. r <= RICB) then
     Qmu=84.6d0
     Qkappa=1327.7d0
!
!--- outer core
!
  else if(r > RICB .and. r <= RCMB) then
     Qmu=0.0d0
     Qkappa=57827.0d0
     if(RCMB - r .LT. r - RICB) then
        Qmu = 312.0d0  ! CMB
     else
        Qmu = 84.6d0   ! ICB
     endif
!
!--- D" at the base of the mantle
!
  else if(r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
     Qmu=312.0d0
     Qkappa=57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if(r > RTOPDDOUBLEPRIME .and. r <= R771) then
     Qmu=312.0d0
     Qkappa=57827.0d0
  else if(r > R771 .and. r <= R670) then
     Qmu=312.0d0
     Qkappa=57827.0d0
!
!--- mantle: above d670
!
  else if(r > R670 .and. r <= R600) then
     Qmu=143.0d0
     Qkappa=57827.0d0
  else if(r > R600 .and. r <= R400) then
     Qmu=143.0d0
     Qkappa=57827.0d0
  else if(r > R400 .and. r <= R220) then
     Qmu=143.0d0
     Qkappa=57827.0d0
  else if(r > R220 .and. r <= R80) then
     Qmu=80.0d0
     Qkappa=57827.0d0
  else if(r > R80) then
     Qmu=600.0d0
     Qkappa=57827.0d0
  endif

  call attenuation_liu(myrank, min_period, max_period, N_SLS, Qmu, T_c_source, tau_s, tau_e)

end subroutine attenuation_model


subroutine attenuation_conversion(myrank, Qmu, T_c_source, tau_s, tau_e)
  use attenuation_model_variables ! includes min_period, max_period, and N_SLS
  implicit none
  integer myrank
  double precision Qmu, T_c_source
  double precision, dimension(N_SLS) :: tau_s, tau_e

  call attenuation_liu(myrank, min_period, max_period, N_SLS, Qmu, T_c_source, tau_s, tau_e)

end subroutine attenuation_conversion

subroutine read_attenuation_model(min, max)

  use attenuation_model_variables
  implicit none

  integer min, max

  min_period = min * 1.0d0
  max_period = max * 1.0d0

end subroutine read_attenuation_model

subroutine attenuation_memory_values(tau_s, deltat, alphaval,betaval,gammaval)
  use attenuation_model_variables
  implicit none

  double precision, dimension(N_SLS) :: tau_s, alphaval, betaval,gammaval
  real(kind=CUSTOM_REAL) deltat

  double precision, dimension(N_SLS) :: tauinv

  tauinv(:) = - 1.0 / tau_s(:)

  alphaval(:)  = 1 + deltat*tauinv(:) + deltat**2*tauinv(:)**2 / 2. + &
       deltat**3*tauinv(:)**3 / 6. + deltat**4*tauinv(:)**4 / 24.
  betaval(:)   = deltat / 2. + deltat**2*tauinv(:) / 3. + deltat**3*tauinv(:)**2 / 8. + deltat**4*tauinv(:)**3 / 24.
  gammaval(:)  = deltat / 2. + deltat**2*tauinv(:) / 6. + deltat**3*tauinv(:)**2 / 24.0

end subroutine attenuation_memory_values

subroutine attenuation_scale_factor(myrank, T_c_source, tau_mu, tau_sigma, Q_mu, scale_factor)
  use attenuation_model_variables
  implicit none

  integer myrank
  double precision scale_factor, Q_mu, T_c_source
  double precision, dimension(N_SLS) :: tau_mu, tau_sigma

  double precision scale_t
  double precision f_c_source, w_c_source, f_0_prem
  double precision factor_scale_mu0, factor_scale_mu
  double precision a_val, b_val
  double precision big_omega
  integer i

  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

!--- compute central angular frequency of source (non dimensionalized)
  f_c_source = ONE / T_c_source
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
  if(scale_factor < 0.9 .or. scale_factor > 1.1) call exit_MPI(myrank,'incorrect correction factor in attenuation model')

end subroutine attenuation_scale_factor

!----

subroutine attenuation_property_values(tau_s, tau_e, factor_common, one_minus_sum_beta)

  use attenuation_model_variables
  implicit none

  double precision, dimension(N_SLS) :: tau_s, tau_e, beta, factor_common
  double precision  one_minus_sum_beta

  double precision, dimension(N_SLS) :: tauinv
  integer i

  tauinv(:) = -1.0d0 / tau_s(:)

  beta(:) = 1.0d0 - tau_e(:) / tau_s(:)
  one_minus_sum_beta = 1.0d0

  do i = 1,N_SLS
     one_minus_sum_beta = one_minus_sum_beta - beta(i)
  enddo

  factor_common(:) = 2.0d0 * beta(:) * tauinv(:)

end subroutine attenuation_property_values

subroutine get_attenuation_model_3D(myrank, prname, one_minus_sum_beta, factor_common, scale_factor, tau_s, vnspec)
  use attenuation_model_variables
  implicit none

  integer myrank, vnspec
  character(len=150) prname
  double precision, dimension(NGLLX,NGLLY,NGLLZ,vnspec)       :: one_minus_sum_beta, scale_factor
  double precision, dimension(N_SLS,NGLLX,NGLLY,NGLLZ,vnspec) :: factor_common
  double precision, dimension(N_SLS)                          :: tau_s

  integer i,j,k,ispec

  double precision, dimension(N_SLS) :: tau_e, fc
  double precision  omsb, Q_mu, sf, T_c_source, scale_t

  ! All of the following reads use the output parameters as their temporary arrays
  ! use the filename to determine the actual contents of the read

  open(unit=27, file=prname(1:len_trim(prname))//'tau_s.bin',status='old',form='unformatted')
  read(27) tau_s
  close(27)

  open(unit=27, file=prname(1:len_trim(prname))//'T_c_source.bin',status='old',form='unformatted')
  read(27) T_c_source
  close(27);

  open(unit=27, file=prname(1:len_trim(prname))//'Q.bin',status='old',form='unformatted')
  read(27) scale_factor
  close(27)

  open(unit=27, file=prname(1:len_trim(prname))//'tau_e.bin',status='old',form='unformatted')
  read(27) factor_common
  close(27)


  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

  factor_common(:,:,:,:,:) = factor_common(:,:,:,:,:) / scale_t
  tau_s(:)                 = tau_s(:) / scale_t
  T_c_source               = 1000.0d0 / T_c_source
  T_c_source               = T_c_source / scale_t

  do ispec = 1, vnspec
     do k = 1, NGLLZ
        do j = 1, NGLLY
           do i = 1, NGLLX
              tau_e(:) = factor_common(:,i,j,k,ispec)
              Q_mu     = scale_factor(i,j,k,ispec)

              ! Determine the factor_common and one_minus_sum_beta from tau_s and tau_e
              call attenuation_property_values(tau_s, tau_e, fc, omsb)

              factor_common(:,i,j,k,ispec)    = fc(:)
              one_minus_sum_beta(i,j,k,ispec) = omsb

              ! Deteremine the "scale_factor" from tau_s, tau_e, central source frequency, and Q
              call attenuation_scale_factor(myrank, T_c_source, tau_e, tau_s, Q_mu, sf)
              scale_factor(i,j,k,ispec) = sf
           enddo
        enddo
     enddo
  enddo
end subroutine get_attenuation_model_3D



subroutine attenuation_liu(myrank, t2, t1, n,Q_real,omega_not,tau_s,tau_e)
  use attenuation_model_variables
  implicit none
  double precision  t1, t2
  integer  myrank, n
  double precision  Q_real
  double precision  omega_not
  double precision, dimension(n)   :: tau_s, tau_e

  integer i, j, k

  double precision, dimension(n)   :: gradient, dadp, dbdp, dqdp, x1, x2
  double precision, dimension(n,n) :: hessian

  double precision a, b, demon, Q_omega, Q_ratio, F_ratio, PI2, q
  double precision f1, f2, exp1, exp2, expo, dexp, omega, df, d2qdp2

  gradient(:)  = 0.0d0
  hessian(:,:) = 0.0d0
  tau_e(:)     = 0.0d0
  tau_s(:)     = 0.0d0

  PI2 = 6.28318530717958d0

  f1 = 1.0d0/t1
  f2 = 1.0d0/t2

  if(f2 < f1) call exit_MPI(myrank, 'max frequency is less than min frequency')

  if(Q_real < 0.0) call exit_MPI(myrank, 'Attenuation is less than zero')

  if(n < 1) call exit_MPI(myrank, 'Number of standard linear solids is less than one')

  omega_not = 1.0e+3 * 10.0**(0.5 * (log10(f1) + log10(f2)))

  exp1 = log10(f1)
  exp2 = log10(f2)

  dexp = (exp2 - exp1) / real(n - 1.0)

  q = 1.0 / (real(n - 1.0) * Q_real )

  do i = 1,n
     expo     = exp1 + (i-1) * dexp
     omega    = PI2 * 10.0**expo
     tau_s(i) = 1.0 / omega
     tau_e(i) = tau_s(i) * (1.0 + q) / (1.0 - q)
  enddo

  x1(:) = tau_e(:) - tau_s(:)
  x2(:) = tau_s(:)

  exp1 = log10(f1);
  exp2 = log10(f2);
  dexp = (exp2 - exp1) / 100.0

  expo = exp1 - dexp
  do i = 1,100
     expo = expo + dexp
     df       = 10.0**(expo+dexp) - 10.0**(expo)
     omega    = PI2 * 10.0**(expo)
     a = real(1.0 - n)
     b = 0.0
     do j = 1,n
        tau_e(j) = x1(j) + x2(j)
        tau_s(j) = x2(j)
        demon   = 1.0 + omega**2.0 * tau_s(j)**2.0
        a       = a + (1.0 + omega**2.0 * tau_e(j) * tau_s(j)) / demon
        b       = b + ( omega * ( tau_e(j) - tau_s(j) ) ) / demon
        dadp(j) = omega**2.0 * tau_s(j) / demon
        dbdp(j) = omega / demon
     enddo

     Q_omega = a / b
     Q_ratio = 1.0 / Q_omega - 1.0 / Q_real
     F_ratio = df / (f2 - f1)
     do j = 1,n
        dqdp(j)     = (dbdp(j) - ( b / a ) * dadp(j)) / a
        gradient(j) = gradient(j) + 2.0 * (Q_ratio) * dqdp(j) * F_ratio
        do k = 1,j
           d2qdp2   = -(dadp(j) * dbdp(k) + dbdp(j) * dadp(k) - 2.0 * (b / a) * dadp(j) * dadp(k)) / (a * a)
           hessian(j,k) = hessian(j,k) + (2.0 * dqdp(j) * dqdp(k) + 2.0 * Q_ratio * d2qdp2) * F_ratio
           hessian(k,j) = hessian(j,k)
        enddo
     enddo
  enddo

  call invert(x1, gradient, hessian, n)
  tau_e(:) = x1(:) + x2(:)
  tau_s(:) = x2(:)

end subroutine attenuation_liu

!----

subroutine invert(x,b,A,n)

  use attenuation_model_variables
  implicit none

  integer n

  double precision, dimension(n)   :: x, b
  double precision, dimension(n,n) :: A

  integer i, j, k
  double precision, dimension(n)   :: W, xp
  double precision, dimension(n,n) :: V
  double precision, dimension(n,n) :: A_inverse

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
  use attenuation_model_variables
  IMPLICIT NONE

  double precision, INTENT(IN) :: a,b
  double precision :: pythag_dp
  double precision :: absa,absb
  absa=abs(a)
  absb=abs(b)
  if (absa > absb) then
     pythag_dp=absa*sqrt(1.0d0+(absb/absa)**2)
  else
     if (absb == 0.0d0) then
        pythag_dp=0.0d0
     else
        pythag_dp=absb*sqrt(1.0+(absa/absb)**2)
     end if
  end if
END FUNCTION pythag_dp

SUBROUTINE svdcmp_dp(a,w,v,p)
!  USE nrtype; USE nrutil, ONLY : assert_eq,nrerror,outerprod
!  USE nr, ONLY : pythag
  use attenuation_model_variables
  IMPLICIT NONE
  integer p
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  double precision, DIMENSION(p,p), INTENT(INOUT) :: a
  double precision, DIMENSION(p), INTENT(OUT) :: w
  double precision, DIMENSION(p,p), INTENT(OUT) :: v
  INTEGER(4) :: i,its,j,k,l,m,n,nm
  double precision :: anorm,c,f,g,h,s,scale,x,y,z
  double precision, DIMENSION(size(a,1)) :: tempm
  double precision, DIMENSION(size(a,2)) :: rv1,tempn
  double precision PYTHAG_DP

  m=size(a,1)
  n = size(a,2)

  g=0.0d0
  scale=0.0d0
  do i=1,n
     l=i+1
     rv1(i)=scale*g
     g=0.0d0
     scale=0.0d0

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

           a(i:m,l:n)=a(i:m,l:n)+spread(a(i:m,i),dim=2,ncopies=size(tempn(l:n))) * &
                spread(tempn(l:n),dim=1,ncopies=size(a(i:m,i)))
           a(i:m,i)=scale*a(i:m,i)
        end if
     end if
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

           a(l:m,l:n)=a(l:m,l:n)+spread(tempm(l:m),dim=2,ncopies=size(rv1(l:n))) * &
                spread(rv1(l:n),dim=1,ncopies=size(tempm(l:m)))
           a(i,l:n)=scale*a(i,l:n)
        end if
     end if
  end do
  anorm=maxval(abs(w)+abs(rv1))

  do i=n,1,-1
     if (i < n) then
        if (g /= 0.0d0) then
           v(l:n,i)=(a(i,l:n)/a(i,l))/g
           tempn(l:n)=matmul(a(i,l:n),v(l:n,l:n))
           v(l:n,l:n)=v(l:n,l:n)+spread(v(l:n,i),dim=2,ncopies=size(tempn(l:n))) * &
                spread(tempn(l:n), dim=1, ncopies=size(v(l:n,i)))
        end if
        v(i,l:n)=0.0d0
        v(l:n,i)=0.0d0
     end if
     v(i,i)=1.0d0
     g=rv1(i)
     l=i
  end do
  do i=min(m,n),1,-1
     l=i+1
     g=w(i)
     a(i,l:n)=0.0d0
     if (g /= 0.0d0) then
        g=1.0d0/g
        tempn(l:n)=(matmul(a(l:m,i),a(l:m,l:n))/a(i,i))*g
        a(i:m,l:n)=a(i:m,l:n)+spread(a(i:m,i),dim=2,ncopies=size(tempn(l:n))) * &
             spread(tempn(l:n),dim=1,ncopies=size(a(i:m,i)))
        a(i:m,i)=a(i:m,i)*g
     else
        a(i:m,i)=0.0d0
     end if
     a(i,i)=a(i,i)+1.0d0
  end do
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
              end do
              exit
           end if
        end do
        z=w(k)
        if (l == k) then
           if (z < 0.0d0) then
              w(k)=-z
              v(1:n,k)=-v(1:n,k)
           end if
           exit
        end if

        if (its == 30) stop 'svdcmp_dp: no convergence in svdcmp'

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
           end if
           f= (c*g)+(s*y)
           x=-(s*g)+(c*y)
           tempm(1:m)=a(1:m,j)
           a(1:m,j)=a(1:m,j)*c+a(1:m,i)*s
           a(1:m,i)=-tempm(1:m)*s+a(1:m,i)*c
        end do
        rv1(l)=0.0d0
        rv1(k)=f
        w(k)=x
     end do
  end do

END SUBROUTINE svdcmp_dp

