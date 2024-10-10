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

!---------------------------------------------------------------------------------
!
! Source time function - UCB Heaviside
!
!---------------------------------------------------------------------------------

module ucb_heaviside

  use shared_input_parameters, only: UCB_SOURCE_T1, UCB_SOURCE_T2, UCB_SOURCE_T3, UCB_SOURCE_T4, UCB_TAU

  implicit none

  private

  ! expose only initialization and
  public :: comp_source_time_function_ucb_stf, init_ucb_heaviside

  ! hard-wired source passband, origin time, and (optional) amplitude are set
  ! here:
  !double precision, parameter :: &
  !        f1h = 1.d0 / 400.d0, &
  !        f2h = 1.d0 / 250.d0, &
  !        f3h = 1.d0 /  53.d0, &
  !        f4h = 1.d0 /  40.d0, &
  !        t0  = 400.d0, &
  !        amp = 1.d7

  ! source amplitude (optional)
  double precision, parameter ::  amp = 1.d0
  double precision :: f1h, f2h, f3h, f4h, t0

  ! the STF type (displacement vs. velocity vs. acceleration) is set below in init_ucb_heaviside

  logical :: initialized = .false.

  integer :: ucb_nstep, ucb_itsource_max

  double precision :: ucb_dt
  double precision, dimension(:), allocatable :: g

  integer, parameter :: CUSTOM_CMPLX = 8   ! for double precision complex

contains

  double precision function comp_source_time_function_ucb_stf(t)

  implicit none
  double precision, intent(in) :: t
  ! local parameters
  integer :: it

  if (.not. initialized) stop "Error [ucb_heaviside]: call to ucb_stf _before_ STF init"

  ! infer the time-step index, always starting from zero time
  it = 1 + NINT(max(0.d0,t) / ucb_dt)

  comp_source_time_function_ucb_stf = g(min(it,ucb_nstep))

  end function comp_source_time_function_ucb_stf

  !---------------------------------------------------------------------------------

  subroutine init_ucb_heaviside(nstep_,dt_)

  use constants, only: myrank

  implicit none

  doubleprecision, parameter :: zero = 0.d0

  integer, intent(in) :: nstep_
  double precision, intent(in) :: dt_

  ! local parameters
  integer :: it0,j,nstep2,i,ic,istat
  double precision :: wt,freq,t1,t2,tmax,seuil

  complex(kind=CUSTOM_CMPLX) :: dphi
  complex(kind=CUSTOM_CMPLX), dimension(:), allocatable :: spectre

  double precision :: pi = 4.d0 * datan(1.d0)

  ! frequencies
  f1h = 1.d0 / UCB_SOURCE_T1
  f2h = 1.d0 / UCB_SOURCE_T2
  f3h = 1.d0 / UCB_SOURCE_T3
  f4h = 1.d0 / UCB_SOURCE_T4

  t0  = UCB_TAU

  ! user info
  if (myrank == 0) then
    write(*,*) "*****************************************"
    write(*,*) "Heaviside Source T1:", UCB_SOURCE_T1
    write(*,*) "Heaviside Source T2:", UCB_SOURCE_T2
    write(*,*) "Heaviside Source T3:", UCB_SOURCE_T3
    write(*,*) "Heaviside Source T4:", UCB_SOURCE_T4
    write(*,*) "Source time-shift:", t0
    write(*,*) "*****************************************"
  endif

  ucb_nstep = nstep_
  ucb_dt = dt_

  allocate(g(ucb_nstep))

  g(:) = zero

  nstep2 = int(2.d0**(int(log(dble(ucb_nstep))/log(2.d0))+1))

  allocate(spectre(nstep2))
  spectre(:) = cmplx(0.d0,0.d0,kind=CUSTOM_CMPLX)
  do j = 1,nstep2
    if (j <= nstep2/2) then
      freq = (j-1)/(ucb_dt * nstep2)
    else if (j == nstep2/2+1) then
      freq = 1/(2.d0 * ucb_dt)
    else
      freq = -(nstep2-j+1)/(ucb_dt * nstep2)
    endif
    dphi = exp(-2.d0 * pi * freq * t0 * cmplx(0.d0,1.d0,kind=CUSTOM_CMPLX))

    call wtcoef(abs(freq),f1h,f2h,f3h,f4h,wt)

    ! >>> STF TYPE IS SET HERE <  <  < 
    ! displacement
    ! if (j /= 1) spectre(j) = wt * dphi / cmplx(0.d0, 2.d0 * pi * freq)
    ! velocity
    if (j /= 1) spectre(j) = wt * dphi
    ! acceleration
    !if (j /= 1) spectre(j) = wt * dphi * cmplx(0.d0, 2.d0 * pi * freq)
  enddo

  call dfour1(spectre,nstep2,1)

  g(:) = amp * real(spectre(1:ucb_nstep))/nstep2/ucb_dt

  !on met les premier pas de temps a zero:
  tmax = nstep2 * ucb_dt
  t1 = 0.d0
  t2 = t0/5.d0
  do i = 1,ucb_nstep
    call wtcoef((i-1) * ucb_dt,t1,t2,tmax,tmax,wt)
    g(i) = g(i) * wt
  enddo

  deallocate(spectre,stat=istat)
  if (istat /= 0) stop 'time_function deallocate error'

  ! calcul du temps de fin de la source
  seuil = 1.d-3 * maxval(abs(g(:)))
  ucb_itsource_max = -1
  ic = 50 !50 pas de temps
  it0 = int(t0/ucb_dt) + 1
  do i = it0,ucb_nstep-ic
    if (ucb_itsource_max < 0) then
      if (maxval(abs(g(i:i+ic))) <= seuil) ucb_itsource_max = i
    endif
  enddo
  if (ucb_itsource_max <= 0) ucb_itsource_max = ucb_nstep-ic

  initialized = .true.

  contains

    subroutine wtcoef(f,f1,f2,f3,f4,wt_out)

    implicit none
    doubleprecision, intent(in) ::  f,f1,f2,f3,f4
    doubleprecision, intent(out)::  wt_out

    if (f3 > f4) stop 'wtcoef: f3 > f4 '
    if (f1 > f2) stop 'wtcoef: f1 > f2 '

    if (f <= f3.and.f >= f2) then
       wt_out = 1.d0
    else if (f > f4.or.f < f1 ) then
       wt_out = 0.d0
    else if (f > f3.and.f <= f4) then
       wt_out = 0.5d0*(1.d0+cos(pi*(f-f3)/(f4-f3)))
    else if (f >= f1.and.f < f2) then
       wt_out = 0.5d0*(1.d0+cos(pi*(f-f2)/(f2-f1)))
    endif

    end subroutine wtcoef

  end subroutine init_ucb_heaviside

  !---------------------------------------------------------------------------------

  subroutine dfour1(u,n,isign)

  implicit none
  integer :: n,isign
  complex(kind=CUSTOM_CMPLX) :: u(n)
  ! local parameters
  double precision :: a(n),b(n),trigs(2*n)
  integer :: inc,jump,lot,i

  inc = 1
  jump = 1
  lot = 1

  call setgpfa(trigs,n)

  do i = 1,n
    a(i) = real(u(i))
    b(i) = aimag(u(i))
  enddo

  call gpfa(a,b,trigs,inc,jump,n,lot,isign)

  do i = 1,n
    u(i) = cmplx( a(i),b(i), kind=CUSTOM_CMPLX)
  enddo

  return
  end

  !---------------------------------------------------------------------------------
  ! subroutine 'SETGPFA'
  ! SETUP ROUTINE FOR SELF-SORTING IN-PLACE
  !     GENERALIZED PRIME FACTOR (COMPLEX) FFT [GPFA]
  !
  ! CALL SETGPFA(TRIGS,N)
  !
  ! INPUT :
  ! -----
  ! N IS THE LENGTH OF THE TRANSFORMS. N MUST BE OF THE FORM:
  !   -----------------------------------
  !     N = (2**IP) * (3**IQ) * (5**IR)
  !   -----------------------------------
  !
  ! OUTPUT:
  ! ------
  ! TRIGS IS A TABLE OF TWIDDLE FACTORS,
  !   OF LENGTH 2*IPQR (REAL) WORDS, WHERE:
  !   --------------------------------------
  !     IPQR = (2**IP) + (3**IQ) + (5**IR)
  !   --------------------------------------
  !
  ! WRITTEN BY CLIVE TEMPERTON 1990
  !
  !----------------------------------------------------------------------

  subroutine SETGPFA(TRIGS,N)

  implicit none
  double precision :: TRIGS(*)
  integer :: N

  double precision :: TWOPI,ANGLE,DEL
  integer, DIMENSION(3) :: NJ

  integer :: I,K,NN,KK,LL,IFAC,IP,IQ,IR,IROT,KINK,NI

!     DECOMPOSE N INTO FACTORS 2,3,5
!     ------------------------------
  NN = N
  IFAC = 2
!
  DO 30 LL = 1 , 3
    KK = 0
10 continue
    if (MOD(NN,IFAC) /= 0) goto 20
    KK = KK + 1
    NN = NN / IFAC
    goto 10
20 continue
    NJ(LL) = KK
    IFAC = IFAC + LL
30 continue

  if (NN /= 1) then
    write(*,40) N
40    FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
    return
  endif

  IP = NJ(1)
  IQ = NJ(2)
  IR = NJ(3)

!     COMPUTE LIST OF ROTATED TWIDDLE FACTORS
!     ---------------------------------------
  NJ(1) = 2**IP
  NJ(2) = 3**IQ
  NJ(3) = 5**IR

  TWOPI = 4.0 * ASIN(1.0)
  I = 1

  DO 60 LL = 1 , 3
    NI = NJ(LL)
    if (NI == 1) goto 60

    DEL = TWOPI / FLOAT(NI)
    IROT = N / NI
    KINK = MOD(IROT,NI)
    KK = 0

    DO 50 K = 1 , NI
      ANGLE = FLOAT(KK) * DEL
      TRIGS(I) = COS(ANGLE)
      TRIGS(I+1) = SIN(ANGLE)
      I = I + 2
      KK = KK + KINK
      if (KK > NI) KK = KK - NI
50 continue
60 continue

  return
  END

  !---------------------------------------------------------------------------------
  ! subroutine 'GPFA'
  ! SELF-SORTING IN-PLACE GENERALIZED PRIME FACTOR (COMPLEX) FFT
  !
  ! *** THIS IS THE ALL-FORTRAN VERSION ***
  !     -------------------------------
  !
  ! CALL GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)
  !
  ! A IS FIRST REAL INPUT/OUTPUT VECTOR
  ! B IS FIRST IMAGINARY INPUT/OUTPUT VECTOR
  ! TRIGS IS A TABLE OF TWIDDLE FACTORS, PRECALCULATED
  !BY CALLING subroutine 'SETGPFA'
  ! INC IS THE INCREMENT WITHIN EACH DATA VECTOR
  ! JUMP IS THE INCREMENT BETWEEN DATA VECTORS
  ! N IS THE LENGTH OF THE TRANSFORMS:
  !   -----------------------------------
  !     N = (2**IP) * (3**IQ) * (5**IR)
  !   -----------------------------------
  ! LOT IS THE NUMBER OF TRANSFORMS
  ! ISIGN = +1 FOR FORWARD TRANSFORM
  != -1 FOR INVERSE TRANSFORM
  !
  ! WRITTEN BY CLIVE TEMPERTON
  ! RECHERCHE EN PREVISION NUMERIQUE
  ! ATMOSPHERIC ENVIRONMENT SERVICE, CANADA
  !
  !----------------------------------------------------------------------
  !
  ! DEFINITION OF TRANSFORM
  ! -----------------------
  !
  ! X(J) = SUM(K=0,...,N-1)(C(K)*EXP(ISIGN*2*I*J*K*PI/N))
  !
  !---------------------------------------------------------------------
  !
  ! FOR A MATHEMATICAL DEVELOPMENT OF THE ALGORITHM USED,
  ! SEE:
  !
  ! C TEMPERTON : "A GENERALIZED PRIME FACTOR FFT ALGORITHM
  !   FOR ANY N = (2**P)(3**Q)(5**R)",
  !   SIAM J. SCI. STAT. COMP., MAY 1992.
  !
  !----------------------------------------------------------------------

  subroutine GPFA(A,B,TRIGS,INC,JUMP,N,LOT,ISIGN)

  implicit none
  double precision :: A(*), B(*), TRIGS(*)
  integer :: INC,JUMP,N,LOT,ISIGN
  integer,DIMENSION(3) :: NJ

  integer :: I,NN,KK,LL,IFAC,IP,IQ,IR

!     DECOMPOSE N INTO FACTORS 2,3,5
!     ------------------------------
  NN = N
  IFAC = 2

  DO 30 LL = 1 , 3
    KK = 0
10 continue
    if (MOD(NN,IFAC) /= 0) goto 20
    KK = KK + 1
    NN = NN / IFAC
    goto 10
20 continue
    NJ(LL) = KK
    IFAC = IFAC + LL
30 continue

  if (NN /= 1) then
    write(*,40) N
40    FORMAT(' *** WARNING!!!',I10,' IS NOT A LEGAL VALUE OF N ***')
    return
  endif

  IP = NJ(1)
  IQ = NJ(2)
  IR = NJ(3)

!     COMPUTE THE TRANSFORM
!     ---------------------
  I = 1
  if (IP > 0) then
    CALL GPFA2F(A,B,TRIGS,INC,JUMP,N,IP,LOT,ISIGN)
    I = I + 2 * ( 2**IP)
  endif
  if (IQ > 0) then
    CALL GPFA3F(A,B,TRIGS(I),INC,JUMP,N,IQ,LOT,ISIGN)
    I = I + 2 * (3**IQ)
  endif
  if (IR > 0) then
    CALL GPFA5F(A,B,TRIGS(I),INC,JUMP,N,IR,LOT,ISIGN)
  endif

  return
  END

  !---------------------------------------------------------------------------------
  !     Fortran version of *gpfa2* -
  !     radix-2 section of self-sorting, in-place, generalized pfa
  !     central radix-2 and radix-8 passes included
  !      so that transform length can be any power of 2
  !
  !---------------------------------------------------------------------------------

  subroutine gpfa2f(a,b,trigs,inc,jump,n,mm,lot,isign)

  implicit none
  double precision :: a(*), b(*), trigs(*)
  integer :: inc,jump,n,mm,lot,isign

  double precision :: s
  double precision :: aja,ajc,t0,t2,ajb,ajd,t1,t3
  double precision :: bja,bjc,u0,u2,bjb,bjd,u1,u3
  double precision :: co1,si1,co2,si2,co3,si3
  double precision :: c1,c2,c3
  double precision :: aje,ajg,ajf,ajh
  double precision :: bje,bjg,bjf,bjh
  double precision :: co4,si4,co5,si5,co6,si6,co7,si7
  double precision :: aji,ajk,ajl,ajm,ajn,ajo,ajp
  double precision :: bji,bjk,bjl,bjm,bjn,bjo,bjp

  integer,parameter :: lvr = 128

  double precision :: ajj
  double precision :: bjj
  double precision :: ss

  integer :: ink,inq,ipass,istart
  integer :: j,ja,jb,jc,jd,je,jf,jg,jh,ji,jj,jjj,jk,jl,jm,jn,jo,jp
  integer :: jstep,jstepl,jstepx
  integer :: k,kk
  integer :: l,la,laincl,left,ll
  integer :: m,m2,m8,mh,mu
  integer :: n2,nb,nblox,ninc,nu,nvex

!     ***************************************************************
!     *     *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *     *
!     ***************************************************************

  n2 = 2**mm
  inq = n/n2
  jstepx = (n2-n) * inc
  ninc = n * inc
  ink = inc * inq

  m = 0
  m2 = 0
  m8 = 0
  if (mod(mm,2) == 0) then
    m = mm/2
  else if (mod(mm,4) == 1) then
    m = (mm-1)/2
    m2 = 1
  else if (mod(mm,4) == 3) then
    m = (mm-3)/2
    m8 = 1
  endif
  mh = (m+1)/2

  nblox = 1 + (lot-1)/lvr
  left = lot
  s = float(isign)
  istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
  do 500 nb = 1 , nblox

    if (left <= lvr) then
      nvex = left
    else if (left < (2*lvr)) then
      nvex = left/2
      nvex = nvex + mod(nvex,2)
    else
      nvex = lvr
    endif
    left = left - nvex

    la = 1
!
!  loop on type I radix-4 passes
!  -----------------------------
    mu = mod(inq,4)
    if (isign == -1) mu = 4 - mu
    ss = 1.0
    if (mu == 3) ss = -1.0

    if (mh == 0) goto 200

    do 160 ipass = 1 , mh
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 4*jstep
        ja = istart + jjj
!
!     "transverse" loop
!     -----------------
        do 115 nu = 1 , inq
          jb = ja + jstepl
          if (jb < istart) jb = jb + ninc
          jc = jb + jstepl
          if (jc < istart) jc = jc + ninc
          jd = jc + jstepl
          if (jd < istart) jd = jd + ninc
          j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
          do 110 l = 1 , nvex
            aja = a(ja+j)
            ajc = a(jc+j)
            t0 = aja + ajc
            t2 = aja - ajc
            ajb = a(jb+j)
            ajd = a(jd+j)
            t1 = ajb + ajd
            t3 = ss * ( ajb - ajd )
            bja = b(ja+j)
            bjc = b(jc+j)
            u0 = bja + bjc
            u2 = bja - bjc
            bjb = b(jb+j)
            bjd = b(jd+j)
            u1 = bjb + bjd
            u3 = ss * ( bjb - bjd )
            a(ja+j) = t0 + t1
            a(jc+j) = t0 - t1
            b(ja+j) = u0 + u1
            b(jc+j) = u0 - u1
            a(jb+j) = t2 - u3
            a(jd+j) = t2 + u3
            b(jb+j) = u2 + t3
            b(jd+j) = u2 - t3
            j = j + jump
110 continue
          ja = ja + jstepx
          if (ja < istart) ja = ja + ninc
115 continue
120 continue
!
!  finished if n2 = 4
!  ------------------
      if (n2 == 4) goto 490
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 150 k = ink , jstep-ink , ink
        co1 = trigs(kk+1)
        si1 = s*trigs(kk+2)
        co2 = trigs(2*kk+1)
        si2 = s*trigs(2*kk+2)
        co3 = trigs(3*kk+1)
        si3 = s*trigs(3*kk+2)
!
!  loop along transform
!  --------------------
        do 140 jjj = k , (n-1)*inc , 4*jstep
          ja = istart + jjj
!
!     "transverse" loop
!     -----------------
          do 135 nu = 1 , inq
            jb = ja + jstepl
            if (jb < istart) jb = jb + ninc
            jc = jb + jstepl
            if (jc < istart) jc = jc + ninc
            jd = jc + jstepl
            if (jd < istart) jd = jd + ninc
            j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep,shortloop
            do 130 l = 1 , nvex
              aja = a(ja+j)
              ajc = a(jc+j)
              t0 = aja + ajc
              t2 = aja - ajc
              ajb = a(jb+j)
              ajd = a(jd+j)
              t1 = ajb + ajd
              t3 = ss * ( ajb - ajd )
              bja = b(ja+j)
              bjc = b(jc+j)
              u0 = bja + bjc
              u2 = bja - bjc
              bjb = b(jb+j)
              bjd = b(jd+j)
              u1 = bjb + bjd
              u3 = ss * ( bjb - bjd )
              a(ja+j) = t0 + t1
              b(ja+j) = u0 + u1
              a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
              b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
              a(jc+j) = co2*(t0-t1) - si2*(u0-u1)
              b(jc+j) = si2*(t0-t1) + co2*(u0-u1)
              a(jd+j) = co3*(t2+u3) - si3*(u2-t3)
              b(jd+j) = si3*(t2+u3) + co3*(u2-t3)
              j = j + jump
130 continue
!-----( end of loop across transforms )
            ja = ja + jstepx
            if (ja < istart) ja = ja + ninc
135 continue
140 continue
!-----( end of loop along transforms )
        kk = kk + 2*la
150 continue
!-----( end of loop on nonzero k )
      la = 4*la
160 continue
!-----( end of loop on type I radix-4 passes)
!
!  central radix-2 pass
!  --------------------
200 continue
    if (m2 == 0) goto 300

    jstep = (n*inc) / (2*la)
    jstepl = jstep - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
    do 220 jjj = 0 , (n-1)*inc , 2*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 215 nu = 1 , inq
        jb = ja + jstepl
        if (jb < istart) jb = jb + ninc
        j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
        do 210 l = 1 , nvex
          aja = a(ja+j)
          ajb = a(jb+j)
          t0 = aja - ajb
          a(ja+j) = aja + ajb
          a(jb+j) = t0
          bja = b(ja+j)
          bjb = b(jb+j)
          u0 = bja - bjb
          b(ja+j) = bja + bjb
          b(jb+j) = u0
          j = j + jump
210 continue
!-----(end of loop across transforms)
        ja = ja + jstepx
        if (ja < istart) ja = ja + ninc
215 continue
220 continue

!  finished if n2=2
!  ----------------
    if (n2 == 2) goto 490

    kk = 2 * la

!  loop on nonzero k
!  -----------------
    do 260 k = ink , jstep - ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)

!  loop along transforms
!  ---------------------
      do 250 jjj = k , (n-1)*inc , 2*jstep
        ja = istart + jjj

!     "transverse" loop
!     -----------------
        do 245 nu = 1 , inq
          jb = ja + jstepl
          if (jb < istart) jb = jb + ninc
          j = 0
!
!  loop across transforms
!  ----------------------
          if (kk == n2/2) then
!cdiri$ ivdep, shortloop
            do 230 l = 1 , nvex
              aja = a(ja+j)
              ajb = a(jb+j)
              t0 = ss * ( aja - ajb )
              a(ja+j) = aja + ajb
              bjb = b(jb+j)
              bja = b(ja+j)
              a(jb+j) = ss * ( bjb - bja )
              b(ja+j) = bja + bjb
              b(jb+j) = t0
              j = j + jump
230 continue

          else

!cdiri$ ivdep, shortloop
            do 240 l = 1 , nvex
              aja = a(ja+j)
              ajb = a(jb+j)
              t0 = aja - ajb
              a(ja+j) = aja + ajb
              bja = b(ja+j)
              bjb = b(jb+j)
              u0 = bja - bjb
              b(ja+j) = bja + bjb
              a(jb+j) = co1*t0 - si1*u0
              b(jb+j) = si1*t0 + co1*u0
              j = j + jump
240 continue

          endif

!-----(end of loop across transforms)
          ja = ja + jstepx
          if (ja < istart) ja = ja + ninc
245 continue
250 continue
!-----(end of loop along transforms)
    kk = kk + 2 * la
260 continue
!-----(end of loop on nonzero k)
!-----(end of radix-2 pass)
!
    la = 2 * la
    goto 400
!
!  central radix-8 pass
!  --------------------
300 continue
    if (m8 == 0) goto 400
    jstep = (n*inc) / (8*la)
    jstepl = jstep - ninc
    mu = mod(inq,8)
    if (isign == -1) mu = 8 - mu
    c1 = 1.0
    if (mu == 3.or.mu == 7) c1 = -1.0
    c2 = sqrt(0.5)
    if (mu == 3.or.mu == 5) c2 = -c2
    c3 = c1 * c2
!
!  stage 1
!  -------
    do 320 k = 0 , jstep - ink , ink
      do 315 jjj = k , (n-1)*inc , 8*jstep
        ja = istart + jjj
!
!     "transverse" loop
!     -----------------
        do 312 nu = 1 , inq
          jb = ja + jstepl
          if (jb < istart) jb = jb + ninc
          jc = jb + jstepl
          if (jc < istart) jc = jc + ninc
          jd = jc + jstepl
          if (jd < istart) jd = jd + ninc
          je = jd + jstepl
          if (je < istart) je = je + ninc
          jf = je + jstepl
          if (jf < istart) jf = jf + ninc
          jg = jf + jstepl
          if (jg < istart) jg = jg + ninc
          jh = jg + jstepl
          if (jh < istart) jh = jh + ninc
          j = 0
!cdiri$ ivdep, shortloop
          do 310 l = 1 , nvex
            aja = a(ja+j)
            aje = a(je+j)
            t0 = aja - aje
            a(ja+j) = aja + aje
            ajc = a(jc+j)
            ajg = a(jg+j)
            t1 = c1 * ( ajc - ajg )
            a(je+j) = ajc + ajg
            ajb = a(jb+j)
            ajf = a(jf+j)
            t2 = ajb - ajf
            a(jc+j) = ajb + ajf
            ajd = a(jd+j)
            ajh = a(jh+j)
            t3 = ajd - ajh
            a(jg+j) = ajd + ajh
            a(jb+j) = t0
            a(jf+j) = t1
            a(jd+j) = c2 * ( t2 - t3 )
            a(jh+j) = c3 * ( t2 + t3 )
            bja = b(ja+j)
            bje = b(je+j)
            u0 = bja - bje
            b(ja+j) = bja + bje
            bjc = b(jc+j)
            bjg = b(jg+j)
            u1 = c1 * ( bjc - bjg )
            b(je+j) = bjc + bjg
            bjb = b(jb+j)
            bjf = b(jf+j)
            u2 = bjb - bjf
            b(jc+j) = bjb + bjf
            bjd = b(jd+j)
            bjh = b(jh+j)
            u3 = bjd - bjh
            b(jg+j) = bjd + bjh
            b(jb+j) = u0
            b(jf+j) = u1
            b(jd+j) = c2 * ( u2 - u3 )
            b(jh+j) = c3 * ( u2 + u3 )
            j = j + jump
310 continue
          ja = ja + jstepx
          if (ja < istart) ja = ja + ninc
312 continue
315 continue
320 continue
!
!  stage 2
!  -------
!
!  k=0 (no twiddle factors)
!  ------------------------
    do 330 jjj = 0 , (n-1)*inc , 8*jstep
      ja = istart + jjj
!
!     "transverse" loop
!     -----------------
      do 328 nu = 1 , inq
        jb = ja + jstepl
        if (jb < istart) jb = jb + ninc
        jc = jb + jstepl
        if (jc < istart) jc = jc + ninc
        jd = jc + jstepl
        if (jd < istart) jd = jd + ninc
        je = jd + jstepl
        if (je < istart) je = je + ninc
        jf = je + jstepl
        if (jf < istart) jf = jf + ninc
        jg = jf + jstepl
        if (jg < istart) jg = jg + ninc
        jh = jg + jstepl
        if (jh < istart) jh = jh + ninc
        j = 0
!cdiri$ ivdep, shortloop
        do 325 l = 1 , nvex
          aja = a(ja+j)
          aje = a(je+j)
          t0 = aja + aje
          t2 = aja - aje
          ajc = a(jc+j)
          ajg = a(jg+j)
          t1 = ajc + ajg
          t3 = c1 * ( ajc - ajg )
          bja = b(ja+j)
          bje = b(je+j)
          u0 = bja + bje
          u2 = bja - bje
          bjc = b(jc+j)
          bjg = b(jg+j)
          u1 = bjc + bjg
          u3 = c1 * ( bjc - bjg )
          a(ja+j) = t0 + t1
          a(je+j) = t0 - t1
          b(ja+j) = u0 + u1
          b(je+j) = u0 - u1
          a(jc+j) = t2 - u3
          a(jg+j) = t2 + u3
          b(jc+j) = u2 + t3
          b(jg+j) = u2 - t3
          ajb = a(jb+j)
          ajd = a(jd+j)
          t0 = ajb + ajd
          t2 = ajb - ajd
          ajf = a(jf+j)
          ajh = a(jh+j)
          t1 = ajf - ajh
          t3 = ajf + ajh
          bjb = b(jb+j)
          bjd = b(jd+j)
          u0 = bjb + bjd
          u2 = bjb - bjd
          bjf = b(jf+j)
          bjh = b(jh+j)
          u1 = bjf - bjh
          u3 = bjf + bjh
          a(jb+j) = t0 - u3
          a(jh+j) = t0 + u3
          b(jb+j) = u0 + t3
          b(jh+j) = u0 - t3
          a(jd+j) = t2 + u1
          a(jf+j) = t2 - u1
          b(jd+j) = u2 - t1
          b(jf+j) = u2 + t1
          j = j + jump
325 continue
        ja = ja + jstepx
        if (ja < istart) ja = ja + ninc
328 continue
330 continue

    if (n2 == 8) goto 490

!  loop on nonzero k
!  -----------------
    kk = 2 * la

    do 350 k = ink , jstep - ink , ink

      co1 = trigs(kk+1)
      si1 = s * trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s * trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s * trigs(3*kk+2)
      co4 = trigs(4*kk+1)
      si4 = s * trigs(4*kk+2)
      co5 = trigs(5*kk+1)
      si5 = s * trigs(5*kk+2)
      co6 = trigs(6*kk+1)
      si6 = s * trigs(6*kk+2)
      co7 = trigs(7*kk+1)
      si7 = s * trigs(7*kk+2)

      do 345 jjj = k , (n-1)*inc , 8*jstep
        ja = istart + jjj

!     "transverse" loop
!     -----------------
        do 342 nu = 1 , inq
          jb = ja + jstepl
          if (jb < istart) jb = jb + ninc
          jc = jb + jstepl
          if (jc < istart) jc = jc + ninc
          jd = jc + jstepl
          if (jd < istart) jd = jd + ninc
          je = jd + jstepl
          if (je < istart) je = je + ninc
          jf = je + jstepl
          if (jf < istart) jf = jf + ninc
          jg = jf + jstepl
          if (jg < istart) jg = jg + ninc
          jh = jg + jstepl
          if (jh < istart) jh = jh + ninc
          j = 0
!cdiri$ ivdep, shortloop
          do 340 l = 1 , nvex
            aja = a(ja+j)
            aje = a(je+j)
            t0 = aja + aje
            t2 = aja - aje
            ajc = a(jc+j)
            ajg = a(jg+j)
            t1 = ajc + ajg
            t3 = c1 * ( ajc - ajg )
            bja = b(ja+j)
            bje = b(je+j)
            u0 = bja + bje
            u2 = bja - bje
            bjc = b(jc+j)
            bjg = b(jg+j)
            u1 = bjc + bjg
            u3 = c1 * ( bjc - bjg )
            a(ja+j) = t0 + t1
            b(ja+j) = u0 + u1
            a(je+j) = co4*(t0-t1) - si4*(u0-u1)
            b(je+j) = si4*(t0-t1) + co4*(u0-u1)
            a(jc+j) = co2*(t2-u3) - si2*(u2+t3)
            b(jc+j) = si2*(t2-u3) + co2*(u2+t3)
            a(jg+j) = co6*(t2+u3) - si6*(u2-t3)
            b(jg+j) = si6*(t2+u3) + co6*(u2-t3)
            ajb = a(jb+j)
            ajd = a(jd+j)
            t0 = ajb + ajd
            t2 = ajb - ajd
            ajf = a(jf+j)
            ajh = a(jh+j)
            t1 = ajf - ajh
            t3 = ajf + ajh
            bjb = b(jb+j)
            bjd = b(jd+j)
            u0 = bjb + bjd
            u2 = bjb - bjd
            bjf = b(jf+j)
            bjh = b(jh+j)
            u1 = bjf - bjh
            u3 = bjf + bjh
            a(jb+j) = co1*(t0-u3) - si1*(u0+t3)
            b(jb+j) = si1*(t0-u3) + co1*(u0+t3)
            a(jh+j) = co7*(t0+u3) - si7*(u0-t3)
            b(jh+j) = si7*(t0+u3) + co7*(u0-t3)
            a(jd+j) = co3*(t2+u1) - si3*(u2-t1)
            b(jd+j) = si3*(t2+u1) + co3*(u2-t1)
            a(jf+j) = co5*(t2-u1) - si5*(u2+t1)
            b(jf+j) = si5*(t2-u1) + co5*(u2+t1)
            j = j + jump
340 continue
          ja = ja + jstepx
          if (ja < istart) ja = ja + ninc
342 continue
345 continue
      kk = kk + 2 * la
350 continue
!
      la = 8 * la
!
!  loop on type II radix-4 passes
!  ------------------------------
400 continue
    mu = mod(inq,4)
    if (isign == -1) mu = 4 - mu
    ss = 1.0
    if (mu == 3) ss = -1.0

    do 480 ipass = mh+1 , m
      jstep = (n*inc) / (4*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc

!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 4*jstep

        do 420 jjj = ll , (n-1)*inc , 4*la*ink
          ja = istart + jjj

!     "transverse" loop
!     -----------------
          do 415 nu = 1 , inq
            jb = ja + jstepl
            if (jb < istart) jb = jb + ninc
            jc = jb + jstepl
            if (jc < istart) jc = jc + ninc
            jd = jc + jstepl
            if (jd < istart) jd = jd + ninc
            je = ja + laincl
            if (je < istart) je = je + ninc
            jf = je + jstepl
            if (jf < istart) jf = jf + ninc
            jg = jf + jstepl
            if (jg < istart) jg = jg + ninc
            jh = jg + jstepl
            if (jh < istart) jh = jh + ninc
            ji = je + laincl
            if (ji < istart) ji = ji + ninc
            jj = ji + jstepl
            if (jj < istart) jj = jj + ninc
            jk = jj + jstepl
            if (jk < istart) jk = jk + ninc
            jl = jk + jstepl
            if (jl < istart) jl = jl + ninc
            jm = ji + laincl
            if (jm < istart) jm = jm + ninc
            jn = jm + jstepl
            if (jn < istart) jn = jn + ninc
            jo = jn + jstepl
            if (jo < istart) jo = jo + ninc
            jp = jo + jstepl
            if (jp < istart) jp = jp + ninc
            j = 0

!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
            do 410 l = 1 , nvex
              aja = a(ja+j)
              ajc = a(jc+j)
              t0 = aja + ajc
              t2 = aja - ajc
              ajb = a(jb+j)
              ajd = a(jd+j)
              t1 = ajb + ajd
              t3 = ss * ( ajb - ajd )
              aji = a(ji+j)
              ajc =  aji
              bja = b(ja+j)
              bjc = b(jc+j)
              u0 = bja + bjc
              u2 = bja - bjc
              bjb = b(jb+j)
              bjd = b(jd+j)
              u1 = bjb + bjd
              u3 = ss * ( bjb - bjd )
              aje = a(je+j)
              ajb =  aje
              a(ja+j) = t0 + t1
              a(ji+j) = t0 - t1
              b(ja+j) = u0 + u1
              bjc =  u0 - u1
              bjm = b(jm+j)
              bjd =  bjm
              a(je+j) = t2 - u3
              ajd =  t2 + u3
              bjb =  u2 + t3
              b(jm+j) = u2 - t3
              !----------------------
              ajg = a(jg+j)
              t0 = ajb + ajg
              t2 = ajb - ajg
              ajf = a(jf+j)
              ajh = a(jh+j)
              t1 = ajf + ajh
              t3 = ss * ( ajf - ajh )
              ajj = a(jj+j)
              ajg =  ajj
              bje = b(je+j)
              bjg = b(jg+j)
              u0 = bje + bjg
              u2 = bje - bjg
              bjf = b(jf+j)
              bjh = b(jh+j)
              u1 = bjf + bjh
              u3 = ss * ( bjf - bjh )
              b(je+j) = bjb
              a(jb+j) = t0 + t1
              a(jj+j) = t0 - t1
              bjj = b(jj+j)
              bjg =  bjj
              b(jb+j) = u0 + u1
              b(jj+j) = u0 - u1
              a(jf+j) = t2 - u3
              ajh =  t2 + u3
              b(jf+j) = u2 + t3
              bjh =  u2 - t3
              !----------------------
              ajk = a(jk+j)
              t0 = ajc + ajk
              t2 = ajc - ajk
              ajl = a(jl+j)
              t1 = ajg + ajl
              t3 = ss * ( ajg - ajl )
              bji = b(ji+j)
              bjk = b(jk+j)
              u0 = bji + bjk
              u2 = bji - bjk
              ajo = a(jo+j)
              ajl =  ajo
              bjl = b(jl+j)
              u1 = bjg + bjl
              u3 = ss * ( bjg - bjl )
              b(ji+j) = bjc
              a(jc+j) = t0 + t1
              a(jk+j) = t0 - t1
              bjo = b(jo+j)
              bjl =  bjo
              b(jc+j) = u0 + u1
              b(jk+j) = u0 - u1
              a(jg+j) = t2 - u3
              a(jo+j) = t2 + u3
              b(jg+j) = u2 + t3
              b(jo+j) = u2 - t3
              !----------------------
              ajm = a(jm+j)
              t0 = ajm + ajl
              t2 = ajm - ajl
              ajn = a(jn+j)
              ajp = a(jp+j)
              t1 = ajn + ajp
              t3 = ss * ( ajn - ajp )
              a(jm+j) = ajd
              u0 = bjd + bjl
              u2 = bjd - bjl
              bjn = b(jn+j)
              bjp = b(jp+j)
              u1 = bjn + bjp
              u3 = ss * ( bjn - bjp )
              a(jn+j) = ajh
              a(jd+j) = t0 + t1
              a(jl+j) = t0 - t1
              b(jd+j) = u0 + u1
              b(jl+j) = u0 - u1
              b(jn+j) = bjh
              a(jh+j) = t2 - u3
              a(jp+j) = t2 + u3
              b(jh+j) = u2 + t3
              b(jp+j) = u2 - t3
              j = j + jump
410 continue
!-----( end of loop across transforms )
          ja = ja + jstepx
          if (ja < istart) ja = ja + ninc
415 continue
420 continue
430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
    if (ipass == m) goto 490

    kk = 2*la

!     loop on nonzero k
!     -----------------
    do 470 k = ink , jstep-ink , ink
      co1 = trigs(kk+1)
      si1 = s*trigs(kk+2)
      co2 = trigs(2*kk+1)
      si2 = s*trigs(2*kk+2)
      co3 = trigs(3*kk+1)
      si3 = s*trigs(3*kk+2)
!
!  double loop along first transform in block
!  ------------------------------------------
      do 460 ll = k , (la-1)*ink , 4*jstep

        do 450 jjj = ll , (n-1)*inc , 4*la*ink
          ja = istart + jjj

!     "transverse" loop
!     -----------------
          do 445 nu = 1 , inq
            jb = ja + jstepl
            if (jb < istart) jb = jb + ninc
            jc = jb + jstepl
            if (jc < istart) jc = jc + ninc
            jd = jc + jstepl
            if (jd < istart) jd = jd + ninc
            je = ja + laincl
            if (je < istart) je = je + ninc
            jf = je + jstepl
            if (jf < istart) jf = jf + ninc
            jg = jf + jstepl
            if (jg < istart) jg = jg + ninc
            jh = jg + jstepl
            if (jh < istart) jh = jh + ninc
            ji = je + laincl
            if (ji < istart) ji = ji + ninc
            jj = ji + jstepl
            if (jj < istart) jj = jj + ninc
            jk = jj + jstepl
            if (jk < istart) jk = jk + ninc
            jl = jk + jstepl
            if (jl < istart) jl = jl + ninc
            jm = ji + laincl
            if (jm < istart) jm = jm + ninc
            jn = jm + jstepl
            if (jn < istart) jn = jn + ninc
            jo = jn + jstepl
            if (jo < istart) jo = jo + ninc
            jp = jo + jstepl
            if (jp < istart) jp = jp + ninc
            j = 0

!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
            do 440 l = 1 , nvex
              aja = a(ja+j)
              ajc = a(jc+j)
              t0 = aja + ajc
              t2 = aja - ajc
              ajb = a(jb+j)
              ajd = a(jd+j)
              t1 = ajb + ajd
              t3 = ss * ( ajb - ajd )
              aji = a(ji+j)
              ajc =  aji
              bja = b(ja+j)
              bjc = b(jc+j)
              u0 = bja + bjc
              u2 = bja - bjc
              bjb = b(jb+j)
              bjd = b(jd+j)
              u1 = bjb + bjd
              u3 = ss * ( bjb - bjd )
              aje = a(je+j)
              ajb =  aje
              a(ja+j) = t0 + t1
              b(ja+j) = u0 + u1
              a(je+j) = co1*(t2-u3) - si1*(u2+t3)
              bjb =  si1*(t2-u3) + co1*(u2+t3)
              bjm = b(jm+j)
              bjd =  bjm
              a(ji+j) = co2*(t0-t1) - si2*(u0-u1)
              bjc =  si2*(t0-t1) + co2*(u0-u1)
              ajd =  co3*(t2+u3) - si3*(u2-t3)
              b(jm+j) = si3*(t2+u3) + co3*(u2-t3)
              !----------------------------------------
              ajg = a(jg+j)
              t0 = ajb + ajg
              t2 = ajb - ajg
              ajf = a(jf+j)
              ajh = a(jh+j)
              t1 = ajf + ajh
              t3 = ss * ( ajf - ajh )
              ajj = a(jj+j)
              ajg =  ajj
              bje = b(je+j)
              bjg = b(jg+j)
              u0 = bje + bjg
              u2 = bje - bjg
              bjf = b(jf+j)
              bjh = b(jh+j)
              u1 = bjf + bjh
              u3 = ss * ( bjf - bjh )
              b(je+j) = bjb
              a(jb+j) = t0 + t1
              b(jb+j) = u0 + u1
              bjj = b(jj+j)
              bjg =  bjj
              a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
              b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
              a(jj+j) = co2*(t0-t1) - si2*(u0-u1)
              b(jj+j) = si2*(t0-t1) + co2*(u0-u1)
              ajh =  co3*(t2+u3) - si3*(u2-t3)
              bjh =  si3*(t2+u3) + co3*(u2-t3)
              !----------------------------------------
              ajk = a(jk+j)
              t0 = ajc + ajk
              t2 = ajc - ajk
              ajl = a(jl+j)
              t1 = ajg + ajl
              t3 = ss * ( ajg - ajl )
              bji = b(ji+j)
              bjk = b(jk+j)
              u0 = bji + bjk
              u2 = bji - bjk
              ajo = a(jo+j)
              ajl =  ajo
              bjl = b(jl+j)
              u1 = bjg + bjl
              u3 = ss * ( bjg - bjl )
              b(ji+j) = bjc
              a(jc+j) = t0 + t1
              b(jc+j) = u0 + u1
              bjo = b(jo+j)
              bjl =  bjo
              a(jg+j) = co1*(t2-u3) - si1*(u2+t3)
              b(jg+j) = si1*(t2-u3) + co1*(u2+t3)
              a(jk+j) = co2*(t0-t1) - si2*(u0-u1)
              b(jk+j) = si2*(t0-t1) + co2*(u0-u1)
              a(jo+j) = co3*(t2+u3) - si3*(u2-t3)
              b(jo+j) = si3*(t2+u3) + co3*(u2-t3)
              !----------------------------------------
              ajm = a(jm+j)
              t0 = ajm + ajl
              t2 = ajm - ajl
              ajn = a(jn+j)
              ajp = a(jp+j)
              t1 = ajn + ajp
              t3 = ss * ( ajn - ajp )
              a(jm+j) = ajd
              u0 = bjd + bjl
              u2 = bjd - bjl
              a(jn+j) = ajh
              bjn = b(jn+j)
              bjp = b(jp+j)
              u1 = bjn + bjp
              u3 = ss * ( bjn - bjp )
              b(jn+j) = bjh
              a(jd+j) = t0 + t1
              b(jd+j) = u0 + u1
              a(jh+j) = co1*(t2-u3) - si1*(u2+t3)
              b(jh+j) = si1*(t2-u3) + co1*(u2+t3)
              a(jl+j) = co2*(t0-t1) - si2*(u0-u1)
              b(jl+j) = si2*(t0-t1) + co2*(u0-u1)
              a(jp+j) = co3*(t2+u3) - si3*(u2-t3)
              b(jp+j) = si3*(t2+u3) + co3*(u2-t3)
              j = j + jump
440 continue
!-----(end of loop across transforms)
            ja = ja + jstepx
            if (ja < istart) ja = ja + ninc
445 continue
450 continue
460 continue
!-----( end of double loop for this k )
      kk = kk + 2*la
470 continue
!-----( end of loop over values of k )
    la = 4*la
480 continue
!-----( end of loop on type II radix-4 passes )
!-----( nvex transforms completed)
490 continue
    istart = istart + nvex * jump
500 continue
!-----( end of loop on blocks of transforms )
  return
  end

  !---------------------------------------------------------------------------------
  !     Fortran version of *gpfa3* -
  !     radix-3 section of self-sorting, in-place
  ! generalized PFA
  !
  !---------------------------------------------------------------------------------

  subroutine gpfa3f(a,b,trigs,inc,jump,n,mm,lot,isign)

  implicit none
  double precision :: a(*), b(*), trigs(*)
  integer :: inc,jump,n,mm,lot,isign

  double precision, parameter :: sin60 = 0.866025403784437

  double precision :: s
  double precision :: aja,ajc,t2,ajb,ajd,t1,t3
  double precision :: bja,bjc,u2,bjb,bjd,u1,u3
  double precision :: co1,si1,co2,si2
  double precision :: c1
  double precision :: aje,ajg,ajf,ajh
  double precision :: bje,bjg,bjf,bjh
  double precision :: bji

  integer, parameter :: lvr = 128

  double precision :: aji

  integer :: ink,inq,ipass,istart
  integer :: j,ja,jb,jc,jd,je,jf,jg,jh,ji,jjj
  integer :: jstep,jstepl,jstepx
  integer :: k,kk
  integer :: l,la,laincl,left,ll
  integer :: m,mh,mu
  integer :: n3,nb,nblox,ninc,nu,nvex

!     ***************************************************************
!     *     *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *     *
!     ***************************************************************

  n3 = 3**mm
  inq = n/n3
  jstepx = (n3-n) * inc
  ninc = n * inc
  ink = inc * inq
  mu = mod(inq,3)
  if (isign == -1) mu = 3-mu
  m = mm
  mh = (m+1)/2
  s = float(isign)
  c1 = sin60
  if (mu == 2) c1 = -c1

  nblox = 1 + (lot-1)/lvr
  left = lot
  s = float(isign)
  istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
  do 500 nb = 1 , nblox

    if (left <= lvr) then
      nvex = left
    else if (left < (2*lvr)) then
      nvex = left/2
      nvex = nvex + mod(nvex,2)
    else
      nvex = lvr
    endif
    left = left - nvex

    la = 1

!  loop on type I radix-3 passes
!  -----------------------------
    do 160 ipass = 1 , mh
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
!
!  k = 0 loop (no twiddle factors)
!  -------------------------------
      do 120 jjj = 0 , (n-1)*inc , 3*jstep
        ja = istart + jjj
!
!  "transverse" loop
!  -----------------
        do 115 nu = 1 , inq
          jb = ja + jstepl
          if (jb < istart) jb = jb + ninc
          jc = jb + jstepl
          if (jc < istart) jc = jc + ninc
          j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
          do 110 l = 1 , nvex
            ajb = a(jb+j)
            ajc = a(jc+j)
            t1 = ajb + ajc
            aja = a(ja+j)
            t2 = aja - 0.5 * t1
            t3 = c1 * ( ajb - ajc )
            bjb = b(jb+j)
            bjc = b(jc+j)
            u1 = bjb + bjc
            bja = b(ja+j)
            u2 = bja - 0.5 * u1
            u3 = c1 * ( bjb - bjc )
            a(ja+j) = aja + t1
            b(ja+j) = bja + u1
            a(jb+j) = t2 - u3
            b(jb+j) = u2 + t3
            a(jc+j) = t2 + u3
            b(jc+j) = u2 - t3
            j = j + jump
  110 continue
          ja = ja + jstepx
          if (ja < istart) ja = ja + ninc
  115 continue
  120 continue
!
!  finished if n3 = 3
!  ------------------
      if (n3 == 3) goto 490
      kk = 2 * la
!
!  loop on nonzero k
!  -----------------
      do 150 k = ink , jstep-ink , ink
        co1 = trigs(kk+1)
        si1 = s*trigs(kk+2)
        co2 = trigs(2*kk+1)
        si2 = s*trigs(2*kk+2)
!
!  loop along transform
!  --------------------
        do 140 jjj = k , (n-1)*inc , 3*jstep
          ja = istart + jjj
!
!  "transverse" loop
!  -----------------
          do 135 nu = 1 , inq
            jb = ja + jstepl
            if (jb < istart) jb = jb + ninc
            jc = jb + jstepl
            if (jc < istart) jc = jc + ninc
            j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep,shortloop
            do 130 l = 1 , nvex
              ajb = a(jb+j)
              ajc = a(jc+j)
              t1 = ajb + ajc
              aja = a(ja+j)
              t2 = aja - 0.5 * t1
              t3 = c1 * ( ajb - ajc )
              bjb = b(jb+j)
              bjc = b(jc+j)
              u1 = bjb + bjc
              bja = b(ja+j)
              u2 = bja - 0.5 * u1
              u3 = c1 * ( bjb - bjc )
              a(ja+j) = aja + t1
              b(ja+j) = bja + u1
              a(jb+j) = co1*(t2-u3) - si1*(u2+t3)
              b(jb+j) = si1*(t2-u3) + co1*(u2+t3)
              a(jc+j) = co2*(t2+u3) - si2*(u2-t3)
              b(jc+j) = si2*(t2+u3) + co2*(u2-t3)
              j = j + jump
  130 continue
!-----( end of loop across transforms )
            ja = ja + jstepx
            if (ja < istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
        kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 3*la
  160 continue
!-----( end of loop on type I radix-3 passes)
!
!  loop on type II radix-3 passes
!  ------------------------------
!  400 continue

    do 480 ipass = mh+1 , m
      jstep = (n*inc) / (3*la)
      jstepl = jstep - ninc
      laincl = la*ink - ninc
!
!  k=0 loop (no twiddle factors)
!  -----------------------------
      do 430 ll = 0 , (la-1)*ink , 3*jstep

        do 420 jjj = ll , (n-1)*inc , 3*la*ink
          ja = istart + jjj
!
!  "transverse" loop
!  -----------------
          do 415 nu = 1 , inq
            jb = ja + jstepl
            if (jb < istart) jb = jb + ninc
            jc = jb + jstepl
            if (jc < istart) jc = jc + ninc
            jd = ja + laincl
            if (jd < istart) jd = jd + ninc
            je = jd + jstepl
            if (je < istart) je = je + ninc
            jf = je + jstepl
            if (jf < istart) jf = jf + ninc
            jg = jd + laincl
            if (jg < istart) jg = jg + ninc
            jh = jg + jstepl
            if (jh < istart) jh = jh + ninc
            ji = jh + jstepl
            if (ji < istart) ji = ji + ninc
            j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
            do 410 l = 1 , nvex
              ajb = a(jb+j)
              ajc = a(jc+j)
              t1 = ajb + ajc
              aja = a(ja+j)
              t2 = aja - 0.5 * t1
              t3 = c1 * ( ajb - ajc )
              ajd = a(jd+j)
              ajb =  ajd
              bjb = b(jb+j)
              bjc = b(jc+j)
              u1 = bjb + bjc
              bja = b(ja+j)
              u2 = bja - 0.5 * u1
              u3 = c1 * ( bjb - bjc )
              bjd = b(jd+j)
              bjb =  bjd
              a(ja+j) = aja + t1
              b(ja+j) = bja + u1
              a(jd+j) = t2 - u3
              b(jd+j) = u2 + t3
              ajc =  t2 + u3
              bjc =  u2 - t3
              !----------------------
              aje = a(je+j)
              ajf = a(jf+j)
              t1 = aje + ajf
              t2 = ajb - 0.5 * t1
              t3 = c1 * ( aje - ajf )
              ajh = a(jh+j)
              ajf =  ajh
              bje = b(je+j)
              bjf = b(jf+j)
              u1 = bje + bjf
              u2 = bjb - 0.5 * u1
              u3 = c1 * ( bje - bjf )
              bjh = b(jh+j)
              bjf =  bjh
              a(jb+j) = ajb + t1
              b(jb+j) = bjb + u1
              a(je+j) = t2 - u3
              b(je+j) = u2 + t3
              a(jh+j) = t2 + u3
              b(jh+j) = u2 - t3
              !----------------------
              aji = a(ji+j)
              t1 = ajf + aji
              ajg = a(jg+j)
              t2 = ajg - 0.5 * t1
              t3 = c1 * ( ajf - aji )
              t1 = ajg + t1
              a(jg+j) = ajc
              bji = b(ji+j)
              u1 = bjf + bji
              bjg = b(jg+j)
              u2 = bjg - 0.5 * u1
              u3 = c1 * ( bjf - bji )
              u1 = bjg + u1
              b(jg+j) = bjc
              a(jc+j) = t1
              b(jc+j) = u1
              a(jf+j) = t2 - u3
              b(jf+j) = u2 + t3
              a(ji+j) = t2 + u3
              b(ji+j) = u2 - t3
              j = j + jump
  410 continue
!-----( end of loop across transforms )
            ja = ja + jstepx
            if (ja < istart) ja = ja + ninc
  415 continue
  420 continue
  430 continue
!-----( end of double loop for k=0 )
!
!  finished if last pass
!  ---------------------
      if (ipass == m) goto 490

      kk = 2*la

!     loop on nonzero k
!     -----------------
      do 470 k = ink , jstep-ink , ink
        co1 = trigs(kk+1)
        si1 = s*trigs(kk+2)
        co2 = trigs(2*kk+1)
        si2 = s*trigs(2*kk+2)

!  double loop along first transform in block
!  ------------------------------------------
        do 460 ll = k , (la-1)*ink , 3*jstep

          do 450 jjj = ll , (n-1)*inc , 3*la*ink
            ja = istart + jjj

!  "transverse" loop
!  -----------------
            do 445 nu = 1 , inq
              jb = ja + jstepl
              if (jb < istart) jb = jb + ninc
              jc = jb + jstepl
              if (jc < istart) jc = jc + ninc
              jd = ja + laincl
              if (jd < istart) jd = jd + ninc
              je = jd + jstepl
              if (je < istart) je = je + ninc
              jf = je + jstepl
              if (jf < istart) jf = jf + ninc
              jg = jd + laincl
              if (jg < istart) jg = jg + ninc
              jh = jg + jstepl
              if (jh < istart) jh = jh + ninc
              ji = jh + jstepl
              if (ji < istart) ji = ji + ninc
              j = 0
!
!  loop across transforms
!  ----------------------
!cdiri$ ivdep, shortloop
              do 440 l = 1 , nvex
                ajb = a(jb+j)
                ajc = a(jc+j)
                t1 = ajb + ajc
                aja = a(ja+j)
                t2 = aja - 0.5 * t1
                t3 = c1 * ( ajb - ajc )
                ajd = a(jd+j)
                ajb =  ajd
                bjb = b(jb+j)
                bjc = b(jc+j)
                u1 = bjb + bjc
                bja = b(ja+j)
                u2 = bja - 0.5 * u1
                u3 = c1 * ( bjb - bjc )
                bjd = b(jd+j)
                bjb =  bjd
                a(ja+j) = aja + t1
                b(ja+j) = bja + u1
                a(jd+j) = co1*(t2-u3) - si1*(u2+t3)
                b(jd+j) = si1*(t2-u3) + co1*(u2+t3)
                ajc =  co2*(t2+u3) - si2*(u2-t3)
                bjc =  si2*(t2+u3) + co2*(u2-t3)
                !----------------------
                aje = a(je+j)
                ajf = a(jf+j)
                t1 = aje + ajf
                t2 = ajb - 0.5 * t1
                t3 = c1 * ( aje - ajf )
                ajh = a(jh+j)
                ajf =  ajh
                bje = b(je+j)
                bjf = b(jf+j)
                u1 = bje + bjf
                u2 = bjb - 0.5 * u1
                u3 = c1 * ( bje - bjf )
                bjh = b(jh+j)
                bjf =  bjh
                a(jb+j) = ajb + t1
                b(jb+j) = bjb + u1
                a(je+j) = co1*(t2-u3) - si1*(u2+t3)
                b(je+j) = si1*(t2-u3) + co1*(u2+t3)
                a(jh+j) = co2*(t2+u3) - si2*(u2-t3)
                b(jh+j) = si2*(t2+u3) + co2*(u2-t3)
                !----------------------
                aji = a(ji+j)
                t1 = ajf + aji
                ajg = a(jg+j)
                t2 = ajg - 0.5 * t1
                t3 = c1 * ( ajf - aji )
                t1 = ajg + t1
                a(jg+j) = ajc
                bji = b(ji+j)
                u1 = bjf + bji
                bjg = b(jg+j)
                u2 = bjg - 0.5 * u1
                u3 = c1 * ( bjf - bji )
                u1 = bjg + u1
                b(jg+j) = bjc
                a(jc+j) = t1
                b(jc+j) = u1
                a(jf+j) = co1*(t2-u3) - si1*(u2+t3)
                b(jf+j) = si1*(t2-u3) + co1*(u2+t3)
                a(ji+j) = co2*(t2+u3) - si2*(u2-t3)
                b(ji+j) = si2*(t2+u3) + co2*(u2-t3)
                j = j + jump
  440 continue
!-----(end of loop across transforms)
              ja = ja + jstepx
              if (ja < istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
        kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 3*la
  480 continue
!-----( end of loop on type II radix-3 passes )
!-----( nvex transforms completed)
  490 continue
    istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
  return
  end

  !---------------------------------------------------------------------------------
  !     Fortran version of *gpfa5* -
  !     radix-5 section of self-sorting, in-place,
  ! generalized pfa
  !
  !-------------------------------------------------------------------

  subroutine gpfa5f(a,b,trigs,inc,jump,n,mm,lot,isign)

  implicit none
  double precision :: a(*), b(*), trigs(*)
  integer :: inc,jump,n,mm,lot,isign

  double precision, parameter :: sin36 = 0.587785252292473
  double precision, parameter :: sin72 = 0.951056516295154
  double precision, parameter :: qrt5 = 0.559016994374947

  double precision :: s
  double precision :: t1,t2,t3,t4,t5,t6,t8,t9,t10,t11
  double precision :: u1,u2,u3,u4,u5,u6,u8,u9,u10,u11
  double precision :: aja,ajc,ajb,ajd
  double precision :: bja,bjc,bjb,bjd
  double precision :: co1,si1,co2,si2,co3,si3
  double precision :: c1,c2,c3
  double precision :: aje,ajg,ajf,ajh
  double precision :: bje,bjg,bjf,bjh
  double precision :: co4,si4
  double precision :: aji,ajk,ajl,ajm,ajn,ajo,ajp
  double precision :: bji,bjk,bjl,bjm,bjn,bjo,bjp

  integer, parameter :: lvr = 128

  double precision :: ajj,ajq,ajr,ajs,ajt,aju,ajv,ajw,ajx,ajy,ax
  double precision :: bjj,bjq,bjr,bjs,bjt,bju,bjv,bjw,bjx,bjy,bx
  double precision :: t7,u7

  integer :: ink,inq,ipass,istart
  integer :: j,ja,jb,jc,jd,je,jf,jg,jh,ji,jj,jjj,jk,jl,jm,jn,jo,jp,jq,jr,js,jt,ju,jv,jw,jx,jy
  integer :: jstep,jstepl,jstepx
  integer :: k,kk
  integer :: l,la,laincl,left,ll
  integer :: m,mh,mu
  integer :: n5,nb,nblox,ninc,nu,nvex

!     ***************************************************************
!     *     *
!     *  N.B. LVR = LENGTH OF VECTOR REGISTERS, SET TO 128 FOR C90. *
!     *  RESET TO 64 FOR OTHER CRAY MACHINES, OR TO ANY LARGE VALUE *
!     *  (GREATER THAN OR EQUAL TO LOT) FOR A SCALAR COMPUTER.      *
!     *     *
!     ***************************************************************

  n5 = 5 ** mm
  inq = n / n5
  jstepx = (n5-n) * inc
  ninc = n * inc
  ink = inc * inq
  mu = mod(inq,5)
  if (isign == -1) mu = 5 - mu

  m = mm
  mh = (m+1)/2
  s = float(isign)
  c1 = qrt5
  c2 = sin72
  c3 = sin36
  if (mu == 2.or.mu == 3) then
    c1 = -c1
    c2 = sin36
    c3 = sin72
  endif
  if (mu == 3.or.mu == 4) c2 = -c2
  if (mu == 2.or.mu == 4) c3 = -c3

  nblox = 1 + (lot-1)/lvr
  left = lot
  s = float(isign)
  istart = 1
!
!  loop on blocks of lvr transforms
!  --------------------------------
  do 500 nb = 1 , nblox

    if (left <= lvr) then
      nvex = left
    else if (left < (2*lvr)) then
      nvex = left/2
      nvex = nvex + mod(nvex,2)
    else
      nvex = lvr
    endif
    left = left - nvex

    la = 1

!  loop on type I radix-5 passes
!  -----------------------------
    do 160 ipass = 1 , mh
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      kk = 0

!  loop on k
!  ---------
      do 150 k = 0 , jstep-ink , ink

        if (k > 0) then
          co1 = trigs(kk+1)
          si1 = s*trigs(kk+2)
          co2 = trigs(2*kk+1)
          si2 = s*trigs(2*kk+2)
          co3 = trigs(3*kk+1)
          si3 = s*trigs(3*kk+2)
          co4 = trigs(4*kk+1)
          si4 = s*trigs(4*kk+2)
        endif

!  loop along transform
!  --------------------
        do 140 jjj = k , (n-1)*inc , 5*jstep
          ja = istart + jjj
!
!     "transverse" loop
!     -----------------
          do 135 nu = 1 , inq
            jb = ja + jstepl
            if (jb < istart) jb = jb + ninc
            jc = jb + jstepl
            if (jc < istart) jc = jc + ninc
            jd = jc + jstepl
            if (jd < istart) jd = jd + ninc
            je = jd + jstepl
            if (je < istart) je = je + ninc
            j = 0
!
!  loop across transforms
!  ----------------------
            if (k == 0) then
!
!cdiri$ ivdep, shortloop
              do 110 l = 1 , nvex
                ajb = a(jb+j)
                aje = a(je+j)
                t1 = ajb + aje
                ajc = a(jc+j)
                ajd = a(jd+j)
                t2 = ajc + ajd
                t3 = ajb - aje
                t4 = ajc - ajd
                t5 = t1 + t2
                t6 = c1 * ( t1 - t2 )
                aja = a(ja+j)
                t7 = aja - 0.25 * t5
                a(ja+j) = aja + t5
                t8 = t7 + t6
                t9 = t7 - t6
                t10 = c3 * t3 - c2 * t4
                t11 = c2 * t3 + c3 * t4
                bjb = b(jb+j)
                bje = b(je+j)
                u1 = bjb + bje
                bjc = b(jc+j)
                bjd = b(jd+j)
                u2 = bjc + bjd
                u3 = bjb - bje
                u4 = bjc - bjd
                u5 = u1 + u2
                u6 = c1 * ( u1 - u2 )
                bja = b(ja+j)
                u7 = bja - 0.25 * u5
                b(ja+j) = bja + u5
                u8 = u7 + u6
                u9 = u7 - u6
                u10 = c3 * u3 - c2 * u4
                u11 = c2 * u3 + c3 * u4
                a(jb+j) = t8 - u11
                b(jb+j) = u8 + t11
                a(je+j) = t8 + u11
                b(je+j) = u8 - t11
                a(jc+j) = t9 - u10
                b(jc+j) = u9 + t10
                a(jd+j) = t9 + u10
                b(jd+j) = u9 - t10
                j = j + jump
  110 continue

            else

!cdiri$ ivdep,shortloop
              do 130 l = 1 , nvex
                ajb = a(jb+j)
                aje = a(je+j)
                t1 = ajb + aje
                ajc = a(jc+j)
                ajd = a(jd+j)
                t2 = ajc + ajd
                t3 = ajb - aje
                t4 = ajc - ajd
                t5 = t1 + t2
                t6 = c1 * ( t1 - t2 )
                aja = a(ja+j)
                t7 = aja - 0.25 * t5
                a(ja+j) = aja + t5
                t8 = t7 + t6
                t9 = t7 - t6
                t10 = c3 * t3 - c2 * t4
                t11 = c2 * t3 + c3 * t4
                bjb = b(jb+j)
                bje = b(je+j)
                u1 = bjb + bje
                bjc = b(jc+j)
                bjd = b(jd+j)
                u2 = bjc + bjd
                u3 = bjb - bje
                u4 = bjc - bjd
                u5 = u1 + u2
                u6 = c1 * ( u1 - u2 )
                bja = b(ja+j)
                u7 = bja - 0.25 * u5
                b(ja+j) = bja + u5
                u8 = u7 + u6
                u9 = u7 - u6
                u10 = c3 * u3 - c2 * u4
                u11 = c2 * u3 + c3 * u4
                a(jb+j) = co1*(t8-u11) - si1*(u8+t11)
                b(jb+j) = si1*(t8-u11) + co1*(u8+t11)
                a(je+j) = co4*(t8+u11) - si4*(u8-t11)
                b(je+j) = si4*(t8+u11) + co4*(u8-t11)
                a(jc+j) = co2*(t9-u10) - si2*(u9+t10)
                b(jc+j) = si2*(t9-u10) + co2*(u9+t10)
                a(jd+j) = co3*(t9+u10) - si3*(u9-t10)
                b(jd+j) = si3*(t9+u10) + co3*(u9-t10)
                j = j + jump
  130 continue

            endif

!-----( end of loop across transforms )

            ja = ja + jstepx
            if (ja < istart) ja = ja + ninc
  135 continue
  140 continue
!-----( end of loop along transforms )
        kk = kk + 2*la
  150 continue
!-----( end of loop on nonzero k )
      la = 5*la
  160 continue
!-----( end of loop on type I radix-5 passes)

    if (n == 5) goto 490
!
!  loop on type II radix-5 passes
!  ------------------------------
!  400 continue

    do 480 ipass = mh+1 , m
      jstep = (n*inc) / (5*la)
      jstepl = jstep - ninc
      laincl = la * ink - ninc
      kk = 0

!     loop on k
!     ---------
      do 470 k = 0 , jstep-ink , ink

        if (k > 0) then
          co1 = trigs(kk+1)
          si1 = s*trigs(kk+2)
          co2 = trigs(2*kk+1)
          si2 = s*trigs(2*kk+2)
          co3 = trigs(3*kk+1)
          si3 = s*trigs(3*kk+2)
          co4 = trigs(4*kk+1)
          si4 = s*trigs(4*kk+2)
        endif

!  double loop along first transform in block
!  ------------------------------------------
        do 460 ll = k , (la-1)*ink , 5*jstep

          do 450 jjj = ll , (n-1)*inc , 5*la*ink
            ja = istart + jjj

!     "transverse" loop
!     -----------------
            do 445 nu = 1 , inq
              jb = ja + jstepl
              if (jb < istart) jb = jb + ninc
              jc = jb + jstepl
              if (jc < istart) jc = jc + ninc
              jd = jc + jstepl
              if (jd < istart) jd = jd + ninc
              je = jd + jstepl
              if (je < istart) je = je + ninc
              jf = ja + laincl
              if (jf < istart) jf = jf + ninc
              jg = jf + jstepl
              if (jg < istart) jg = jg + ninc
              jh = jg + jstepl
              if (jh < istart) jh = jh + ninc
              ji = jh + jstepl
              if (ji < istart) ji = ji + ninc
              jj = ji + jstepl
              if (jj < istart) jj = jj + ninc
              jk = jf + laincl
              if (jk < istart) jk = jk + ninc
              jl = jk + jstepl
              if (jl < istart) jl = jl + ninc
              jm = jl + jstepl
              if (jm < istart) jm = jm + ninc
              jn = jm + jstepl
              if (jn < istart) jn = jn + ninc
              jo = jn + jstepl
              if (jo < istart) jo = jo + ninc
              jp = jk + laincl
              if (jp < istart) jp = jp + ninc
              jq = jp + jstepl
              if (jq < istart) jq = jq + ninc
              jr = jq + jstepl
              if (jr < istart) jr = jr + ninc
              js = jr + jstepl
              if (js < istart) js = js + ninc
              jt = js + jstepl
              if (jt < istart) jt = jt + ninc
              ju = jp + laincl
              if (ju < istart) ju = ju + ninc
              jv = ju + jstepl
              if (jv < istart) jv = jv + ninc
              jw = jv + jstepl
              if (jw < istart) jw = jw + ninc
              jx = jw + jstepl
              if (jx < istart) jx = jx + ninc
              jy = jx + jstepl
              if (jy < istart) jy = jy + ninc
              j = 0

!  loop across transforms
!  ----------------------
              if (k == 0) then

!cdiri$ ivdep, shortloop
                do 410 l = 1 , nvex
                  ajb = a(jb+j)
                  aje = a(je+j)
                  t1 = ajb + aje
                  ajc = a(jc+j)
                  ajd = a(jd+j)
                  t2 = ajc + ajd
                  t3 = ajb - aje
                  t4 = ajc - ajd
                  ajf = a(jf+j)
                  ajb =  ajf
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  aja = a(ja+j)
                  t7 = aja - 0.25 * t5
                  a(ja+j) = aja + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  ajk = a(jk+j)
                  ajc =  ajk
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  bjb = b(jb+j)
                  bje = b(je+j)
                  u1 = bjb + bje
                  bjc = b(jc+j)
                  bjd = b(jd+j)
                  u2 = bjc + bjd
                  u3 = bjb - bje
                  u4 = bjc - bjd
                  bjf = b(jf+j)
                  bjb =  bjf
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  bja = b(ja+j)
                  u7 = bja - 0.25 * u5
                  b(ja+j) = bja + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  bjk = b(jk+j)
                  bjc =  bjk
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  a(jf+j) = t8 - u11
                  b(jf+j) = u8 + t11
                  aje =  t8 + u11
                  bje =  u8 - t11
                  a(jk+j) = t9 - u10
                  b(jk+j) = u9 + t10
                  ajd =  t9 + u10
                  bjd =  u9 - t10
                  !----------------------
                  ajg = a(jg+j)
                  ajj = a(jj+j)
                  t1 = ajg + ajj
                  ajh = a(jh+j)
                  aji = a(ji+j)
                  t2 = ajh + aji
                  t3 = ajg - ajj
                  t4 = ajh - aji
                  ajl = a(jl+j)
                  ajh =  ajl
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  t7 = ajb - 0.25 * t5
                  a(jb+j) = ajb + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  ajq = a(jq+j)
                  aji =  ajq
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  bjg = b(jg+j)
                  bjj = b(jj+j)
                  u1 = bjg + bjj
                  bjh = b(jh+j)
                  bji = b(ji+j)
                  u2 = bjh + bji
                  u3 = bjg - bjj
                  u4 = bjh - bji
                  bjl = b(jl+j)
                  bjh =  bjl
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  u7 = bjb - 0.25 * u5
                  b(jb+j) = bjb + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  bjq = b(jq+j)
                  bji =  bjq
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  a(jg+j) = t8 - u11
                  b(jg+j) = u8 + t11
                  ajj =  t8 + u11
                  bjj =  u8 - t11
                  a(jl+j) = t9 - u10
                  b(jl+j) = u9 + t10
                  a(jq+j) = t9 + u10
                  b(jq+j) = u9 - t10
                  !----------------------
                  ajo = a(jo+j)
                  t1 = ajh + ajo
                  ajm = a(jm+j)
                  ajn = a(jn+j)
                  t2 = ajm + ajn
                  t3 = ajh - ajo
                  t4 = ajm - ajn
                  ajr = a(jr+j)
                  ajn =  ajr
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  t7 = ajc - 0.25 * t5
                  a(jc+j) = ajc + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  ajw = a(jw+j)
                  ajo =  ajw
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  bjo = b(jo+j)
                  u1 = bjh + bjo
                  bjm = b(jm+j)
                  bjn = b(jn+j)
                  u2 = bjm + bjn
                  u3 = bjh - bjo
                  u4 = bjm - bjn
                  bjr = b(jr+j)
                  bjn =  bjr
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  u7 = bjc - 0.25 * u5
                  b(jc+j) = bjc + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  bjw = b(jw+j)
                  bjo =  bjw
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  a(jh+j) = t8 - u11
                  b(jh+j) = u8 + t11
                  a(jw+j) = t8 + u11
                  b(jw+j) = u8 - t11
                  a(jm+j) = t9 - u10
                  b(jm+j) = u9 + t10
                  a(jr+j) = t9 + u10
                  b(jr+j) = u9 - t10
                  !----------------------
                  ajt = a(jt+j)
                  t1 = aji + ajt
                  ajs = a(js+j)
                  t2 = ajn + ajs
                  t3 = aji - ajt
                  t4 = ajn - ajs
                  ajx = a(jx+j)
                  ajt =  ajx
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  ajp = a(jp+j)
                  t7 = ajp - 0.25 * t5
                  ax = ajp + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  a(jp+j) = ajd
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  a(jd+j) = ax
                  bjt = b(jt+j)
                  u1 = bji + bjt
                  bjs = b(js+j)
                  u2 = bjn + bjs
                  u3 = bji - bjt
                  u4 = bjn - bjs
                  bjx = b(jx+j)
                  bjt =  bjx
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  bjp = b(jp+j)
                  u7 = bjp - 0.25 * u5
                  bx = bjp + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  b(jp+j) = bjd
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  b(jd+j) = bx
                  a(ji+j) = t8 - u11
                  b(ji+j) = u8 + t11
                  a(jx+j) = t8 + u11
                  b(jx+j) = u8 - t11
                  a(jn+j) = t9 - u10
                  b(jn+j) = u9 + t10
                  a(js+j) = t9 + u10
                  b(js+j) = u9 - t10
                  !----------------------
                  ajv = a(jv+j)
                  ajy = a(jy+j)
                  t1 = ajv + ajy
                  t2 = ajo + ajt
                  t3 = ajv - ajy
                  t4 = ajo - ajt
                  a(jv+j) = ajj
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  aju = a(ju+j)
                  t7 = aju - 0.25 * t5
                  ax = aju + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  a(ju+j) = aje
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  a(je+j) = ax
                  bjv = b(jv+j)
                  bjy = b(jy+j)
                  u1 = bjv + bjy
                  u2 = bjo + bjt
                  u3 = bjv - bjy
                  u4 = bjo - bjt
                  b(jv+j) = bjj
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  bju = b(ju+j)
                  u7 = bju - 0.25 * u5
                  bx = bju + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  b(ju+j) = bje
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  b(je+j) = bx
                  a(jj+j) = t8 - u11
                  b(jj+j) = u8 + t11
                  a(jy+j) = t8 + u11
                  b(jy+j) = u8 - t11
                  a(jo+j) = t9 - u10
                  b(jo+j) = u9 + t10
                  a(jt+j) = t9 + u10
                  b(jt+j) = u9 - t10
                  j = j + jump
  410 continue

              else

!cdiri$ ivdep, shortloop
                do 440 l = 1 , nvex
                  ajb = a(jb+j)
                  aje = a(je+j)
                  t1 = ajb + aje
                  ajc = a(jc+j)
                  ajd = a(jd+j)
                  t2 = ajc + ajd
                  t3 = ajb - aje
                  t4 = ajc - ajd
                  ajf = a(jf+j)
                  ajb =  ajf
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  aja = a(ja+j)
                  t7 = aja - 0.25 * t5
                  a(ja+j) = aja + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  ajk = a(jk+j)
                  ajc =  ajk
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  bjb = b(jb+j)
                  bje = b(je+j)
                  u1 = bjb + bje
                  bjc = b(jc+j)
                  bjd = b(jd+j)
                  u2 = bjc + bjd
                  u3 = bjb - bje
                  u4 = bjc - bjd
                  bjf = b(jf+j)
                  bjb =  bjf
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  bja = b(ja+j)
                  u7 = bja - 0.25 * u5
                  b(ja+j) = bja + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  bjk = b(jk+j)
                  bjc =  bjk
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  a(jf+j) = co1*(t8-u11) - si1*(u8+t11)
                  b(jf+j) = si1*(t8-u11) + co1*(u8+t11)
                  aje =  co4*(t8+u11) - si4*(u8-t11)
                  bje =  si4*(t8+u11) + co4*(u8-t11)
                  a(jk+j) = co2*(t9-u10) - si2*(u9+t10)
                  b(jk+j) = si2*(t9-u10) + co2*(u9+t10)
                  ajd =  co3*(t9+u10) - si3*(u9-t10)
                  bjd =  si3*(t9+u10) + co3*(u9-t10)
                  !----------------------
                  ajg = a(jg+j)
                  ajj = a(jj+j)
                  t1 = ajg + ajj
                  ajh = a(jh+j)
                  aji = a(ji+j)
                  t2 = ajh + aji
                  t3 = ajg - ajj
                  t4 = ajh - aji
                  ajl = a(jl+j)
                  ajh =  ajl
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  t7 = ajb - 0.25 * t5
                  a(jb+j) = ajb + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  ajq = a(jq+j)
                  aji =  ajq
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  bjg = b(jg+j)
                  bjj = b(jj+j)
                  u1 = bjg + bjj
                  bjh = b(jh+j)
                  bji = b(ji+j)
                  u2 = bjh + bji
                  u3 = bjg - bjj
                  u4 = bjh - bji
                  bjl = b(jl+j)
                  bjh =  bjl
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  u7 = bjb - 0.25 * u5
                  b(jb+j) = bjb + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  bjq = b(jq+j)
                  bji =  bjq
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  a(jg+j) = co1*(t8-u11) - si1*(u8+t11)
                  b(jg+j) = si1*(t8-u11) + co1*(u8+t11)
                  ajj =  co4*(t8+u11) - si4*(u8-t11)
                  bjj =  si4*(t8+u11) + co4*(u8-t11)
                  a(jl+j) = co2*(t9-u10) - si2*(u9+t10)
                  b(jl+j) = si2*(t9-u10) + co2*(u9+t10)
                  a(jq+j) = co3*(t9+u10) - si3*(u9-t10)
                  b(jq+j) = si3*(t9+u10) + co3*(u9-t10)
                  !----------------------
                  ajo = a(jo+j)
                  t1 = ajh + ajo
                  ajm = a(jm+j)
                  ajn = a(jn+j)
                  t2 = ajm + ajn
                  t3 = ajh - ajo
                  t4 = ajm - ajn
                  ajr = a(jr+j)
                  ajn =  ajr
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  t7 = ajc - 0.25 * t5
                  a(jc+j) = ajc + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  ajw = a(jw+j)
                  ajo =  ajw
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  bjo = b(jo+j)
                  u1 = bjh + bjo
                  bjm = b(jm+j)
                  bjn = b(jn+j)
                  u2 = bjm + bjn
                  u3 = bjh - bjo
                  u4 = bjm - bjn
                  bjr = b(jr+j)
                  bjn =  bjr
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  u7 = bjc - 0.25 * u5
                  b(jc+j) = bjc + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  bjw = b(jw+j)
                  bjo =  bjw
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  a(jh+j) = co1*(t8-u11) - si1*(u8+t11)
                  b(jh+j) = si1*(t8-u11) + co1*(u8+t11)
                  a(jw+j) = co4*(t8+u11) - si4*(u8-t11)
                  b(jw+j) = si4*(t8+u11) + co4*(u8-t11)
                  a(jm+j) = co2*(t9-u10) - si2*(u9+t10)
                  b(jm+j) = si2*(t9-u10) + co2*(u9+t10)
                  a(jr+j) = co3*(t9+u10) - si3*(u9-t10)
                  b(jr+j) = si3*(t9+u10) + co3*(u9-t10)
                  !----------------------
                  ajt = a(jt+j)
                  t1 = aji + ajt
                  ajs = a(js+j)
                  t2 = ajn + ajs
                  t3 = aji - ajt
                  t4 = ajn - ajs
                  ajx = a(jx+j)
                  ajt =  ajx
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  ajp = a(jp+j)
                  t7 = ajp - 0.25 * t5
                  ax = ajp + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  a(jp+j) = ajd
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  a(jd+j) = ax
                  bjt = b(jt+j)
                  u1 = bji + bjt
                  bjs = b(js+j)
                  u2 = bjn + bjs
                  u3 = bji - bjt
                  u4 = bjn - bjs
                  bjx = b(jx+j)
                  bjt =  bjx
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  bjp = b(jp+j)
                  u7 = bjp - 0.25 * u5
                  bx = bjp + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  b(jp+j) = bjd
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  b(jd+j) = bx
                  a(ji+j) = co1*(t8-u11) - si1*(u8+t11)
                  b(ji+j) = si1*(t8-u11) + co1*(u8+t11)
                  a(jx+j) = co4*(t8+u11) - si4*(u8-t11)
                  b(jx+j) = si4*(t8+u11) + co4*(u8-t11)
                  a(jn+j) = co2*(t9-u10) - si2*(u9+t10)
                  b(jn+j) = si2*(t9-u10) + co2*(u9+t10)
                  a(js+j) = co3*(t9+u10) - si3*(u9-t10)
                  b(js+j) = si3*(t9+u10) + co3*(u9-t10)
                  !----------------------
                  ajv = a(jv+j)
                  ajy = a(jy+j)
                  t1 = ajv + ajy
                  t2 = ajo + ajt
                  t3 = ajv - ajy
                  t4 = ajo - ajt
                  a(jv+j) = ajj
                  t5 = t1 + t2
                  t6 = c1 * ( t1 - t2 )
                  aju = a(ju+j)
                  t7 = aju - 0.25 * t5
                  ax = aju + t5
                  t8 = t7 + t6
                  t9 = t7 - t6
                  a(ju+j) = aje
                  t10 = c3 * t3 - c2 * t4
                  t11 = c2 * t3 + c3 * t4
                  a(je+j) = ax
                  bjv = b(jv+j)
                  bjy = b(jy+j)
                  u1 = bjv + bjy
                  u2 = bjo + bjt
                  u3 = bjv - bjy
                  u4 = bjo - bjt
                  b(jv+j) = bjj
                  u5 = u1 + u2
                  u6 = c1 * ( u1 - u2 )
                  bju = b(ju+j)
                  u7 = bju - 0.25 * u5
                  bx = bju + u5
                  u8 = u7 + u6
                  u9 = u7 - u6
                  b(ju+j) = bje
                  u10 = c3 * u3 - c2 * u4
                  u11 = c2 * u3 + c3 * u4
                  b(je+j) = bx
                  a(jj+j) = co1*(t8-u11) - si1*(u8+t11)
                  b(jj+j) = si1*(t8-u11) + co1*(u8+t11)
                  a(jy+j) = co4*(t8+u11) - si4*(u8-t11)
                  b(jy+j) = si4*(t8+u11) + co4*(u8-t11)
                  a(jo+j) = co2*(t9-u10) - si2*(u9+t10)
                  b(jo+j) = si2*(t9-u10) + co2*(u9+t10)
                  a(jt+j) = co3*(t9+u10) - si3*(u9-t10)
                  b(jt+j) = si3*(t9+u10) + co3*(u9-t10)
                  j = j + jump
  440 continue

              endif
!
!-----(end of loop across transforms)
!
              ja = ja + jstepx
              if (ja < istart) ja = ja + ninc
  445 continue
  450 continue
  460 continue
!-----( end of double loop for this k )
        kk = kk + 2*la
  470 continue
!-----( end of loop over values of k )
      la = 5*la
  480 continue
!-----( end of loop on type II radix-5 passes )
!-----( nvex transforms completed)
  490 continue
    istart = istart + nvex * jump
  500 continue
!-----( end of loop on blocks of transforms )
!
  return
  end

end module ucb_heaviside

