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

  subroutine FFT(npow,xi,zign,dtt,mpow)

! Fourier transform
! This inputs AND outputs a complex function.
! The convention is FFT --> e^(-iwt)
! numerical factor for Plancherel theorem: planch_fac = dble(NPT * dt * dt)

  use constants, only: CUSTOM_REAL

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer,intent(in) :: npow

  complex(kind=CUSTOM_CMPLX),dimension(*) :: xi
  real(kind=CUSTOM_REAL),intent(in) :: zign,dtt

! here is the hardwired maximum size of the array
! Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
! Aug 2016: this should be moved to the calling program and precomputed once and for all
  real(kind=CUSTOM_REAL),intent(in) :: mpow(30)

  ! local parameters
  integer :: lblock,k,FK,jh,ii,istart
  integer :: l,iblock,nblock,i,lbhalf,j,lx
  complex(kind=CUSTOM_CMPLX) :: wk, hold, q
  real(kind=CUSTOM_REAL) :: flx,inv_of_flx,v

  real(kind=CUSTOM_REAL), parameter :: PI = acos(-1.0)

  ! added this sanity check
  if (npow > 30) stop 'Error: FTT routine has an hardwired maximum of 30 levels'
  ! in any case the line below imposes a maximum of 31, otherwise the integer 2**n will overflow

  lx = 2**npow
  flx = lx

  do l = 1,npow

    nblock = 2**(l-1)
    lblock = lx/nblock
    lbhalf = lblock/2

    k = 0
    ii = 0

    do iblock = 1,nblock

      FK = k
      flx = lx

      v = zign * 2.0_CUSTOM_REAL * PI * FK / flx         ! Fourier convention

      ! - sign: MATLAB convention: forward e^{-i om t}
      ! + sign: engineering convention: forward e^{i om t}
      wk = cmplx(cos(v),-sin(v))   ! sign change to -sin(v) or sin(v)
      istart = lblock*(iblock-1)

      do i = 1,lbhalf
        j  = istart+i
        jh = j+lbhalf
        q = xi(jh)*wk
        xi(jh) = xi(j)-q
        xi(j)  = xi(j)+q
      enddo

      do i = 2,npow
        ii = i
        if (k < int(mpow(i))) goto 4
        k = k - int(mpow(i))
      enddo

    4 k = k + int(mpow(ii))

    enddo
  enddo

  k = 0
  ii = 0

  do j = 1,lx
    if (k < j) goto 5

    hold = xi(j)
    xi(j) = xi(k+1)
    xi(k+1) = hold

5   do i = 1,npow
      ii = i
      if (k < int(mpow(i))) goto 7
      k = k - int(mpow(i))
    enddo

7   k = k + int(mpow(ii))
  enddo

  ! final steps deal with dt factors
  if (zign > 0.0_CUSTOM_REAL) then
    ! FORWARD FFT
    xi(1:lx) = xi(1:lx) * dtt    ! multiplication by dt
  else
    ! REVERSE FFT
    flx = flx*dtt
    inv_of_flx = 1._CUSTOM_REAL / flx
! Aug 2016: changed to multiplication by the precomputed inverse to make the routine faster
!       xi(1:lx) = xi(1:lx) / flx         ! division by dt
    xi(1:lx) = xi(1:lx) * inv_of_flx  ! division by dt
  endif

  end subroutine FFT

!
!-------------------------------------------------------------------------------------------------
!

  subroutine FFTinv(npow,s,zign,dtt,r,mpow)

! inverse Fourier transform -- calls FFT

  use constants

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer,intent(in) :: npow
  real(kind=CUSTOM_REAL),intent(in)  :: dtt,zign

  complex(kind=CUSTOM_CMPLX) :: s(*)
  real(kind=CUSTOM_REAL) :: r(*)       ! note that this is real, not double precision

! here is the hardwired maximum size of the array
! Aug 2016: if this routine is called many times (for different mesh points at which the SEM is coupled with FK)
! Aug 2016: this should be moved to the calling program and precomputed once and for all
  real(kind=CUSTOM_REAL) :: mpow(30)

  ! local parameters
  integer :: nsmp, nhalf

  nsmp = 2**npow
  nhalf = nsmp/2

  call rspec(s,nhalf)               ! restructuring

  call FFT(npow,s,zign,dtt,mpow)    ! Fourier transform

  r(1:nsmp) = real(s(1:nsmp),kind=CUSTOM_REAL)       ! take the real part

  end subroutine FFTinv

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rspec(s,np2)

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  complex(kind=CUSTOM_CMPLX) :: s(*)
  integer, intent(in) :: np2

  ! local parameters
  integer :: n,n1,i

  n = 2*np2
  n1 = np2 + 1

  s(n1) = 0.0
  s(1)  = cmplx(real(s(1)),0.0,kind=CUSTOM_CMPLX)

  do i = 1,np2
    s(np2+i) = conjg(s(np2+2-i))
  enddo

  end subroutine rspec

!
!-------------------------------------------------------------------------------------------------
!

  subroutine fft_apply_3D_symmetry(array,N)

! applies symmetries and restructuring of complex 3D array(N,N,N) for inverse FFT

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer,intent(in) :: N
  complex(kind=CUSTOM_CMPLX),dimension(N,N,N),intent(inout) :: array

  ! local parameters
  integer :: i,j,k

  ! define symmetry conditions and nugget effect
  array(1,1,(/1,N/2+1/)) = 0.0
  array(N/2+1,1,(/1,N/2+1/)) = abs(array(N/2+1,1,(/1,N/2+1/)))
  array(1,N/2+1,(/1,N/2+1/)) = abs(array(1,N/2+1,(/1,N/2+1/)))
  array(N/2+1,N/2+1,(/1,N/2+1/)) = abs(array(N/2+1,N/2+1,(/1,N/2+1/)))

! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(array) &
!$OMP PRIVATE(i,j,k) &
!$OMP FIRSTPRIVATE(N)

!$OMP DO
  do j = 2,N/2
    array(1,N-j+2,(/1,N/2+1/)) = conjg(array(1,j,(/1,N/2+1/)))
  enddo
!$OMP ENDDO

!$OMP DO
  do i = 2,N/2
    array(N-i+2,1,(/1,N/2+1/)) = conjg(array(i,1,(/1,N/2+1/)))
  enddo
!$OMP ENDDO

!$OMP DO
  do j = 2,N/2
    array(N/2+1,N-j+2,(/1,N/2+1/)) = conjg(array(N/2+1,j,(/1,N/2+1/)))
  enddo
!$OMP ENDDO

!$OMP DO
  do i = 2,N/2
    array(N-i+2,N/2+1,(/1,N/2+1/)) = conjg(array(i,N/2+1,(/1,N/2+1/)))
  enddo
!$OMP ENDDO

!$OMP DO
  do j = 2,N/2
    do i = 2,N/2
      array(N-i+2,N-j+2,(/1,N/2+1/)) = conjg(array(i,j,(/1,N/2+1/)))
      array(N-i+2,N/2-j+2,(/1,N/2+1/)) = conjg(array(i,N/2+j,(/1,N/2+1/)))
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 2,N/2
    array(1,1,N-k+2) = conjg(array(1,1,k))
    array(1,N/2+1,N-k+2) = conjg(array(1,N/2+1,k))
    array(N/2+1,1,N-k+2) = conjg(array(N/2+1,1,k))
    array(N/2+1,N/2+1,N-k+2) = conjg(array(N/2+1,N/2+1,k))
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 2,N/2
    do j = 2,N/2
      array(1,N-j+2,N-k+2) = conjg(array(1,j,k))
      array(1,j,N-k+2) = conjg(array(1,N-j+2,k))
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 2,N/2
    do i = 2,N/2
      array(N-i+2,1,N-k+2) = conjg(array(i,1,k))
      array(i,1,N-k+2) = conjg(array(N-i+2,1,k))
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 2,N/2
    do j = 2,N/2
      array(N/2+1,N-j+2,N-k+2) = conjg(array(N/2+1,j,k))
      array(N/2+1,j,N-k+2) = conjg(array(N/2+1,N-j+2,k))
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 2,N/2
    do i = 2,N/2
      array(N-i+2,N/2+1,N-k+2) = conjg(array(i,N/2+1,k))
      array(i,N/2+1,N-k+2) = conjg(array(N-i+2,N/2+1,k))
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 2,N/2
    do j = 2,N/2
      do i = 2,N/2
        array(N-i+2,N-j+2,N-k+2) = conjg(array(i,j,k))
        array(N-i+2,N/2-j+2,N-k+2) = conjg(array(i,N/2+j,k))
        array(i,j,N-k+2) = conjg(array(N-i+2,N-j+2,k))
        array(i,N/2+j,N-k+2) = conjg(array(N-i+2,N/2-j+2,k))
      enddo
    enddo
  enddo
!$OMP ENDDO

!$OMP END PARALLEL

  end subroutine fft_apply_3D_symmetry

!
!-------------------------------------------------------------------------------------------------
!


  subroutine FFT_3D(N,npower_of_2,array,zign,dk,mpow)

  use constants, only: CUSTOM_REAL

  implicit none

  integer, parameter :: CUSTOM_CMPLX = 8

  integer,intent(in) :: N,npower_of_2
  complex(kind=CUSTOM_CMPLX),dimension(N,N,N),intent(inout) :: array

  real(kind=CUSTOM_REAL),intent(in) :: zign,dk

  ! pre-computed helper array
  real(kind=CUSTOM_REAL),intent(in) :: mpow(30)

  ! local parameters
  integer :: i,j,k

  ! checks
  if (npower_of_2 > 30) stop 'Error: the FTTinv_3D routine has an hardwired maximum of 30 levels'

! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(array) &
!$OMP PRIVATE(i,j,k) &
!$OMP FIRSTPRIVATE(N,npower_of_2,zign,dk,mpow)

!$OMP DO
  do k = 1,N
    do j = 1,N
      call FFT(npower_of_2, array(:,j,k), zign, dk, mpow)
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do k = 1,N
    do i = 1,N
      call FFT(npower_of_2, array(i,:,k), zign, dk, mpow)
    enddo
  enddo
!$OMP ENDDO

!$OMP DO
  do j = 1,N
    do i = 1,N
      call FFT(npower_of_2, array(i,j,:), zign, dk, mpow)
    enddo
  enddo
!$OMP ENDDO

!$OMP END PARALLEL

  end subroutine FFT_3D
