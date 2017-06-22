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


! we put these multiplications in a separate routine AND IN A SEPARATE FILE because otherwise
! some compilers try to unroll the six loops above and take forever to compile

! we leave this in a separate file otherwise many compilers perform subroutine inlining when
! two subroutines are in the same file and one calls the other


!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------


  subroutine multiply_accel_elastic(two_omega_earth, &
                                    NGLOB_CM,veloc_cm,accel_cm, &
                                    rmassx_cm,rmassy_cm,rmassz_cm, &
                                    NGLOB_IC,veloc_ic,accel_ic, &
                                    rmassx_ic,rmassy_ic,rmassz_ic)

! multiplies acceleration with inverse of mass matrices in crust/mantle,solid inner core region

  use constants_solver, only: CUSTOM_REAL,NDIM,ROTATION_VAL

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: two_omega_earth

  ! crust/mantle
  integer,intent(in) :: NGLOB_CM
  ! velocity & acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM),intent(in) :: veloc_cm
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM),intent(inout) :: accel_cm
  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(NGLOB_CM),intent(in) :: rmassx_cm,rmassy_cm,rmassz_cm

  ! inner core
  integer,intent(in) :: NGLOB_IC
  ! velocity & acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC),intent(in) :: veloc_ic
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC),intent(inout) :: accel_ic
  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(NGLOB_IC),intent(in) :: rmassx_ic,rmassy_ic,rmassz_ic

  ! local parameters
  integer :: i

  ! note: mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete

  ! updates acceleration w/ rotation in elastic region

  ! see input call, differs for corrected mass matrices for rmassx,rmassy,rmassz

  ! divides by mass matrix
  ! (rmassx holds inverted mass matrix; numerically multiplication is faster than division)
  if (ROTATION_VAL) then
    ! adds contributions due to rotation
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(two_omega_earth, &
!$OMP NGLOB_CM, accel_cm, veloc_cm, rmassx_cm, rmassy_cm, rmassz_cm, &
!$OMP NGLOB_IC, accel_ic, veloc_ic, rmassx_ic, rmassy_ic, rmassz_ic) &
!$OMP PRIVATE(i)

    ! crust/mantle
!$OMP DO
    do i = 1,NGLOB_CM
      accel_cm(1,i) = accel_cm(1,i)*rmassx_cm(i) + two_omega_earth * veloc_cm(2,i)
      accel_cm(2,i) = accel_cm(2,i)*rmassy_cm(i) - two_omega_earth * veloc_cm(1,i)
      accel_cm(3,i) = accel_cm(3,i)*rmassz_cm(i)
    enddo
!$OMP enddo NOWAIT

    ! inner core
!$OMP DO
    do i = 1,NGLOB_IC
      accel_ic(1,i) = accel_ic(1,i)*rmassx_ic(i) + two_omega_earth * veloc_ic(2,i)
      accel_ic(2,i) = accel_ic(2,i)*rmassy_ic(i) - two_omega_earth * veloc_ic(1,i)
      accel_ic(3,i) = accel_ic(3,i)*rmassz_ic(i)
    enddo
!$OMP enddo
!$OMP END PARALLEL

  else
    ! no rotation
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(two_omega_earth, &
!$OMP NGLOB_CM, accel_cm, veloc_cm, rmassx_cm, rmassy_cm, rmassz_cm, &
!$OMP NGLOB_IC, accel_ic, veloc_ic, rmassx_ic, rmassy_ic, rmassz_ic) &
!$OMP PRIVATE(i)

    ! crust/mantle
!$OMP DO
    do i = 1,NGLOB_CM
      accel_cm(1,i) = accel_cm(1,i)*rmassx_cm(i)
      accel_cm(2,i) = accel_cm(2,i)*rmassy_cm(i)
      accel_cm(3,i) = accel_cm(3,i)*rmassz_cm(i)
    enddo
!$OMP enddo NOWAIT

    ! inner core
!$OMP DO
    do i = 1,NGLOB_IC
      accel_ic(1,i) = accel_ic(1,i)*rmassx_ic(i)
      accel_ic(2,i) = accel_ic(2,i)*rmassy_ic(i)
      accel_ic(3,i) = accel_ic(3,i)*rmassz_ic(i)
    enddo
!$OMP enddo
!$OMP END PARALLEL

  endif

  end subroutine multiply_accel_elastic



!-------------------------------------------------------------------------------------------------
!
! acoustic/fluid domains
!
!-------------------------------------------------------------------------------------------------


  subroutine multiply_accel_acoustic(NGLOB,accel,rmass)

! multiplies acceleration with inverse of mass matrix in outer core region

  use constants_solver, only: CUSTOM_REAL

  implicit none

  integer :: NGLOB

  ! velocity & acceleration
  ! crust/mantle region
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: accel

  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: rmass

  ! local parameters
  integer :: i

  ! note: mass matrices for fluid region has no Stacey or rotation correction
  !       it is also the same for forward and backward/reconstructed wavefields

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(NGLOB, accel, rmass) &
!$OMP PRIVATE(i)
!$OMP DO
  do i = 1,NGLOB
    accel(i) = accel(i)*rmass(i)
  enddo
!$OMP enddo
!$OMP END PARALLEL

  end subroutine multiply_accel_acoustic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine multiply_arrays_adjoint(sourcearrayd,hxir,hetar,hgammar,adj_src_ud)

  use constants

  implicit none

  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX) :: hxir
  double precision, dimension(NGLLY) :: hetar
  double precision, dimension(NGLLZ) :: hgammar
  double precision, dimension(NDIM) :: adj_src_ud

  integer :: i,j,k

  ! adds interpolated source contribution to all GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        sourcearrayd(:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src_ud(:)
      enddo
    enddo
  enddo

  end subroutine multiply_arrays_adjoint

