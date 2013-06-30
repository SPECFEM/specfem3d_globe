

module attenuation_simplex_variables
  implicit none

  ! nf    = Number of Frequencies
  ! nsls  = Number of Standard Linear Solids
  integer  nf, nsls

  ! f     = Frequencies at which to evaluate the solution
  real(8), allocatable, dimension(:) ::  f

  ! Q     = Desired Value of Attenuation or Q
  ! iQ    = 1/Q
  real(8)  Q, iQ

  ! tau_s = Tau_sigma defined by the frequency range and
  !             number of standard linear solids
  real(8), allocatable, dimension(:) :: tau_s

end module attenuation_simplex_variables

subroutine attenuation_simplex(t1, t2, n, Q_real, omega_not, tau_s, tau_e)
  implicit none
  ! Input / Output
  real(8)  t1, t2
  real(8)  Q_real
  real(8)  omega_not
  integer  n
  real(8), dimension(n)   :: tau_s, tau_e

  ! Internal
  integer i, iterations, err,prnt
  real(8) f1, f2, exp1,exp2,dexp, min_value
  real(8), allocatable, dimension(:) :: f
  real(8), parameter :: PI = 3.14159265358979d0
  integer, parameter :: nf = 100
  real(8) attenuation_eval
  EXTERNAL attenuation_eval

  ! Values to be passed into the simplex minimization routine
  iterations = -1
  min_value  = -1.0e-4
  err        = 0
  prnt       = 0

  allocate(f(nf))

  ! Determine the min and max frequencies
  f1 = 1.0d0 / t1
  f2 = 1.0d0 / t2

  ! Determine the exponents of the frequencies
  exp1 = log10(f1);
  exp2 = log10(f2);

  if(f2 < f1 .OR. Q_real < 0.0d0 .OR. n < 1) then
     write(*,*)'bad parameters'
     call exit(-1)
  endif

  ! Determine the Source frequency
  omega_not =  1.0e+03 * 10.0d0**(0.5 * (log10(f1) + log10(f2)))

  ! Determine the Frequencies at which to compare solutions
  !   The frequencies should be equally spaced in log10 frequency
  do i = 1,nf
     f(i) = exp1 + ((i-1)*1.0d0 * (exp2-exp1) / ((nf-1)*1.0d0))
  enddo

  ! Set the Tau_sigma (tau_s) to be equally spaced in log10 frequency
  dexp = (exp2-exp1) / ((n*1.0d0) - 1)
  do i = 1,n
     tau_s(i) = 1.0 / (PI * 2.0d0 * 10**(exp1 + (i - 1)* 1.0d0 *dexp))
  enddo

  ! Shove the paramters into the module
  call attenuation_simplex_setup(nf,n,f,Q_real,tau_s)

  ! Set the Tau_epsilon (tau_e) to an initial value
  ! at omega*tau = 1; tan_delta = 1/Q = (tau_e - tau_s)/(2 * sqrt(tau e*tau_s))
  !    if we assume tau_e =~ tau_s
  !    we get the equation below
  do i = 1,n
     tau_e(i) = tau_s(i) + (tau_s(i) * 2.0d0/Q_real)
  enddo

  ! Run a simplex search to determine the optimum values of tau_e
  call fminsearch(attenuation_eval, tau_e, n, iterations, min_value, prnt, err)
  if(err > 0) then
     write(*,*)'Search did not converge for an attenuation of ', Q_real
     write(*,*)'    Iterations: ', iterations
     write(*,*)'    Min Value:  ', min_value
     write(*,*)'    Aborting program'
     call exit(-1)
  endif

  deallocate(f)
  call attenuation_simplex_finish()

end subroutine attenuation_simplex

subroutine attenuation_simplex_finish()
  use attenuation_simplex_variables
  implicit none

  deallocate(f)
  deallocate(tau_s)

end subroutine attenuation_simplex_finish

!!!!!!!
! subroutine simplex_setup
!   - Inserts necessary parameters into the module attenuation_simplex_variables
!   - See module for explaination
subroutine attenuation_simplex_setup(nf_in,nsls_in,f_in,Q_in,tau_s_in)
  use attenuation_simplex_variables
  implicit none

  integer nf_in, nsls_in
  real(8) Q_in
  real(8), dimension(nf_in)   :: f_in
  real(8), dimension(nsls_in) :: tau_s_in

  allocate(f(nf_in))
  allocate(tau_s(nsls_in))

  nf    = nf_in
  nsls  = nsls_in
  f     = f_in
  Q     = Q_in
  iQ    = 1.0d0/Q
  tau_s = tau_s_in

end subroutine attenuation_simplex_setup

!!!!!!!!
! subroutine attenuation_maxwell
!   - Computes the Moduli (Maxwell Solid) for a series of
!         Standard Linear Solids
!   - Computes M1 and M2 parameters after Dahlen and Tromp pp.203
!         here called B and A after Liu et al. 1976
!   - Another formulation uses Kelvin-Voigt Solids and computes
!         Compliences J1 and J2 after Dahlen and Tromp pp.203
!
!   Input
!     nf    = Number of Frequencies
!     nsls  = Number of Standard Linear Solids
!     f     = Frequencies (in log10 of frequencies)
!                dimension(nf)
!     tau_s = Tau_sigma  Stress relaxation time (see References)
!                dimension(nsls)
!     tau_e = Tau_epislon Strain relaxation time (see References)
!                dimension(nsls)!
!   Output
!     B     = Real Moduli      ( M2 Dahlen and Tromp pp.203 )
!                dimension(nf)
!     A     = Imaginary Moduli ( M1 Dahlen and Tromp pp.203 )
!                dimension(nf)
!
!   Dahlen and Tromp, 1998
!      Theoritical Global Seismology
!
!   Liu et al. 1976
!      Velocity dispersion due to anelasticity; implications for seismology and mantle composition
!      Geophys, J. R. asts. Soc, Vol 47, pp. 41-58
subroutine attenuation_maxwell(nf,nsls,f,tau_s,tau_e,B,A)
  implicit none

  ! Input
  integer nf, nsls
  real(8), dimension(nf)   :: f
  real(8), dimension(nsls) :: tau_s, tau_e
  ! Output
  real(8), dimension(nf)   :: A,B

  integer i,j
  real(8) w, pi, demon

  PI = 3.14159265358979d0

  A(:) = 1.0d0 -  nsls*1.0d0
  B(:) = 0.0d0
  do i = 1,nf
     w = 2.0d0 * PI * 10**f(i)
     do j = 1,nsls
!        write(*,*)j,tau_s(j),tau_e(j)
        demon = 1.0d0 + w**2 * tau_s(j)**2
        A(i) = A(i) + ((1.0d0 + (w**2 * tau_e(j) * tau_s(j)))/ demon)
        B(i) = B(i) + ((w * (tau_e(j) - tau_s(j))) / demon)
     enddo
!     write(*,*)A(i),B(i),10**f(i)
  enddo

end subroutine attenuation_maxwell

!!!!!!!!
! subroutine attenuation_eval
!    - Computes the misfit from a set of relaxation paramters
!          given a set of frequencies and target attenuation
!    - Evaluates only at the given frequencies
!    - Evaluation is done with an L2 norm
!
!    Input
!      Xin = Tau_epsilon, Strain Relaxation Time
!                Note: Tau_sigma the Stress Relaxation Time is loaded
!                      with attenuation_simplex_setup and stored in
!                      attenuation_simplex_variables
!
!    Xi = Sum_i^N sqrt [ (1/Qc_i - 1/Qt_i)^2 / 1/Qt_i^2 ]
!
!     where Qc_i is the computed attenuation at a specific frequency
!           Qt_i is the desired attenuaiton at that frequency
!
!    Uses attenuation_simplex_variables to store constant values
!
!    See atteunation_simplex_setup
!
real(8) function attenuation_eval(Xin)
  use attenuation_simplex_variables
  implicit none
   ! Input
  real(8), dimension(nsls) :: Xin
  real(8), dimension(nsls) :: tau_e

  real(8), dimension(nf)   :: A, B, tan_delta

  integer i
  real(8) xi, iQ2

  tau_e = Xin

  call attenuation_maxwell(nf,nsls,f,tau_s,tau_e,B,A)

  tan_delta = B / A

  attenuation_eval = 0.0d0
  iQ2 = iQ**2
  do i = 1,nf
     xi = sqrt(( ( (tan_delta(i) - iQ) ** 2 ) / iQ2 ))
     attenuation_eval = attenuation_eval + xi
  enddo

end function attenuation_eval


!!!!!!!!!!!!!1
! subroutine fminsearch
!   - Computes the minimization of funk(x(n)) using the simplex method
!   - This subroutine is copied from Matlab fminsearch.m
!         and modified to suit my nefarious needs
!   Input
!     funk = real(8) function with one input parameter
!                real(8) function the_funk(x)
!     x    = Input/Output
!               variables to be minimized
!               dimension(n)
!            Input:  Initial Value
!            Output: Mimimized Value
!     n    = number of variables
!     itercount = Input/Output
!                 Input:  maximum number of iterations
!                         if < 0 default is used (200 * n)
!                 Output: total number of iterations on output
!     tolf      = Input/Output
!                 Input:  minimium tolerance of the function funk(x)
!                 Output: minimium value of funk(x)(i.e. "a" solution)
!     prnt      = Input
!                 3 => report every iteration
!                 4 => report every iteration, total simplex
!     err       = Output
!                 0 => Normal exeecution, converged within desired range
!                 1 => Function Evaluation exceeded limit
!                 2 => Iterations exceeded limit
!
!     See Matlab fminsearch
!
subroutine fminsearch(funk, x, n, itercount, tolf, prnt, err)
  implicit none

  ! Input
  real(8) funk
  EXTERNAL funk

  integer n
  real(8) x(n) ! Also Output
  integer itercount, prnt, err
  real(8) tolf

  !Internal
  integer i,j, how
  integer, parameter :: none             = 0
  integer, parameter :: initial          = 1
  integer, parameter :: expand           = 2
  integer, parameter :: reflect          = 3
  integer, parameter :: contract_outside = 4
  integer, parameter :: contract_inside  = 5
  integer, parameter :: shrink           = 6

  integer maxiter, maxfun
  integer func_evals
  real(8) tolx

  real(8) rho, chi, psi, sigma
  real(8) xin(n), y(n), v(n,n+1), fv(n+1)
  real(8) vtmp(n,n+1)
  real(8) usual_delta, zero_term_delta
  real(8) xbar(n), xr(n), fxr, xe(n), fxe, xc(n), fxc, fxcc, xcc(n)
  integer place(n+1)

  real(8) max_size_simplex, max_value

  rho   = 1.0d0
  chi   = 2.0d0
  psi   = 0.5d0
  sigma = 0.5d0


  if(itercount > 0) then
     maxiter = itercount
  else
     maxiter = 200 * n
  endif
  itercount = 0
  maxfun  = 200 * n

  if(tolf > 0.0d0) then
     tolx = 1.0e-4
  else
     tolx = 1.0e-4
     tolf = 1.0e-4
  endif

  err = 0

  xin    = x
  v(:,:) = 0.0d0
  fv(:)  = 0.0d0

  v(:,1) = xin
  x      = xin

  fv(1) = funk(xin)

  usual_delta = 0.05
  zero_term_delta = 0.00025

  do j = 1,n
     y = xin
     if(y(j) /= 0.0d0) then
        y(j) = (1.0d0 + usual_delta) * y(j)
     else
        y(j) = zero_term_delta
     endif
     v(:,j+1) = y
     x(:) = y
     fv(j+1) = funk(x)
  enddo

  call bubble_sort(fv,n+1,place)

  do i = 1,n+1
     vtmp(:,i) = v(:,place(i))
  enddo
  v = vtmp

  how = initial
  itercount = 1
  func_evals = n+1
  if(prnt == 3) then
     write(*,*)'Iterations   Funk Evals   Value How'
     write(*,*)itercount, func_evals, fv(1), how
  endif
  if(prnt == 4) then
     write(*,*)'How: ',how
     write(*,*)'V: ', v
     write(*,*)'fv: ',fv
     write(*,*)'evals: ',func_evals
  endif

  do while (func_evals < maxfun .AND. itercount < maxiter)

!     if(max(max(abs(v(:,2:n+1) - v(:,1)))) <= tolx .AND. &
!          max(abs(fv(1) - fv(2:n+1))) <= tolf) then

     if(max_size_simplex(v,n) <= tolx .AND. &
          max_value(fv,n+1) <= tolf) then
        goto 666
     endif
     how = none

     ! xbar = average of the n (NOT n+1) best points
     !     xbar = sum(v(:,1:n), 2)/n
     xbar(:) = 0.0d0
     do i = 1,n
        do j = 1,n
           xbar(i) = xbar(i) + v(i,j)
        enddo
        xbar(i) = xbar(i) / (n*1.0d0)
     enddo
     xr = (1 + rho)*xbar - rho*v(:,n+1)
     x(:) = xr
     fxr = funk(x)
     func_evals = func_evals + 1
     if (fxr < fv(1)) then
        ! Calculate the expansion point
        xe = (1 + rho*chi)*xbar - rho*chi*v(:,n+1)
        x = xe
        fxe = funk(x)
        func_evals = func_evals+1
        if (fxe < fxr) then
           v(:,n+1) = xe
           fv(n+1) = fxe
           how = expand
        else
           v(:,n+1) = xr
           fv(n+1) = fxr
           how = reflect
        endif
     else ! fv(:,1) <= fxr
        if (fxr < fv(n)) then
           v(:,n+1) = xr
           fv(n+1) = fxr
           how = reflect
        else ! fxr >= fv(:,n)
           ! Perform contraction
           if (fxr < fv(n+1)) then
              ! Perform an outside contraction
              xc = (1 + psi*rho)*xbar - psi*rho*v(:,n+1)
              x(:) = xc
              fxc = funk(x)
              func_evals = func_evals+1

              if (fxc <= fxr) then
                 v(:,n+1) = xc
                 fv(n+1) = fxc
                 how = contract_outside
              else
                 ! perform a shrink
                 how = shrink
              endif
           else
              ! Perform an inside contraction
              xcc = (1-psi)*xbar + psi*v(:,n+1)
              x(:) = xcc
              fxcc = funk(x)
              func_evals = func_evals+1

              if (fxcc < fv(n+1)) then
                 v(:,n+1) = xcc
                 fv(n+1) = fxcc
                 how = contract_inside
              else
                 ! perform a shrink
                 how = shrink
              endif
           endif
           if (how == shrink) then
              do j=2,n+1
                 v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1))
                 x(:) = v(:,j)
                 fv(j) = funk(x)
              enddo
              func_evals = func_evals + n
           endif
        endif
     endif

     call bubble_sort(fv,n+1,place)
     do i = 1,n+1
        vtmp(:,i) = v(:,place(i))
     enddo
     v = vtmp

     itercount = itercount + 1
     if (prnt == 3) then
        write(*,*)itercount, func_evals, fv(1), how
     else if (prnt == 4) then
        write(*,*)
        write(*,*)'How: ',how
        write(*,*)'v: ',v
        write(*,*)'fv: ',fv
        write(*,*)'evals: ',func_evals
     endif
  enddo

  if(func_evals > maxfun) then
     write(*,*)'function evaluations exceeded prescribed limit', maxfun
     err = 1
  endif
  if(itercount > maxiter) then
     write(*,*)'iterations exceeded prescribed limit', maxiter
     err = 2
  endif

666 continue
  x = v(:,1)
  tolf = fv(1)

end subroutine fminsearch


!!!!!!!
! real(8) function max_value
!    - Finds the maximim value of the difference of between the first
!          value and the remaining values of a vector
!    Input
!      fv = Input
!             Vector
!             dimension(n)
!      n  = Input
!             Length of fv
!
!      Returns:
!         Xi = max( || fv(1)- fv(i) || ); i=2:n
!
real(8) function max_value(fv,n)
  implicit none
  integer n
  real(8) fv(n)

  integer i
  real(8) m, z

  m = 0.0d0
  do i = 2,n
     z = abs(fv(1) - fv(i))
     if(z > m) then
        m = z
     endif
  enddo

  max_value = m

end function max_value

!!!!!!!!
! function max_size_simplex
!   - Determines the maximum distance between two point in a simplex
!   Input
!     v  = Input
!            Simplex Verticies
!            dimension(n, n+1)
!     n  = Pseudo Length of n
!
!     Returns:
!       Xi = max( max( || v(:,1) - v(:,i) || ) ) ; i=2:n+1
!
real(8) function max_size_simplex(v,n)
  implicit none
  integer n
  real(8) v(n,n+1)

  integer i,j
  real(8) m, z

  m = 0.0d0
  do i = 1,n
     do j = 2,n+1
        z = abs(v(i,j) - v(i,1))
        if(z > m) then
           m = z
        endif
     enddo
  enddo

  max_size_simplex = m

end function max_size_simplex


!!!!!!!
!    - Implementation of a Bubble Sort Routine
!    Input
!      X = Input/Output
!         Vector to be sorted
!         dimension(n)
!      n = Input
!         Length of X
!      I = Output
!         Sorted Indicies of vecotr X
!
!      Example:
!         X = [ 4 3 1 2 ] on Input
!         I = [ 1 2 3 4 ] Computed Internally (in order)
!
!         X = [ 1 2 3 4 ] on Output
!         I = [ 3 4 2 1 ] on Output
!
subroutine bubble_sort(X,n,I)
  implicit none
  integer n
  real(8) X(n)
  integer I(n)

  integer j,k
  real(8) rtmp
  integer itmp

  do j = 1,n
     I(j) = j
  enddo

  do j = 1,n
     do k = 1,n-j
        if(X(k+1) < X(k)) then
           rtmp   = X(k)
           X(k)   = X(k+1)
           X(k+1) = rtmp

           itmp   = I(k)
           I(k)   = I(k+1)
           I(k+1) = itmp
        endif
     enddo
  enddo

end subroutine bubble_sort



