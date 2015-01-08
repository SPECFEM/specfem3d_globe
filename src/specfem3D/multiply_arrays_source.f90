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


! we put these multiplications in a separate routine AND IN A SEPARATE FILE because otherwise
! some compilers try to unroll the six loops above and take forever to compile

! we leave this in a separate file otherwise many compilers perform subroutine inlining when
! two subroutines are in the same file and one calls the other


!-------------------------------------------------------------------------------------------------
!
! elastic domains
!
!-------------------------------------------------------------------------------------------------


  subroutine multiply_accel_elastic(NGLOB,veloc,accel, &
                                    two_omega_earth, &
                                    rmassx,rmassy,rmassz)

! multiplies acceleration with inverse of mass matrices in crust/mantle,solid inner core region

  use constants_solver,only: CUSTOM_REAL,NDIM

  implicit none

  integer :: NGLOB

  ! velocity & acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: veloc,accel

  real(kind=CUSTOM_REAL) :: two_omega_earth

  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: rmassx,rmassy,rmassz

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
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(NGLOB, accel, rmassx, rmassy, rmassz, two_omega_earth, veloc) &
!$OMP PRIVATE(i)
!$OMP DO SCHEDULE(GUIDED)
  do i = 1,NGLOB
    accel(1,i) = accel(1,i)*rmassx(i) + two_omega_earth*veloc(2,i)
    accel(2,i) = accel(2,i)*rmassy(i) - two_omega_earth*veloc(1,i)
    accel(3,i) = accel(3,i)*rmassz(i)
  enddo
!$OMP enddo
!$OMP END PARALLEL

  end subroutine multiply_accel_elastic



!-------------------------------------------------------------------------------------------------
!
! acoustic/fluid domains
!
!-------------------------------------------------------------------------------------------------


  subroutine multiply_accel_acoustic(NGLOB,accel,rmass)

! multiplies acceleration with inverse of mass matrix in outer core region

  use constants_solver,only: CUSTOM_REAL

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

  do i = 1,NGLOB
    accel(i) = accel(i)*rmass(i)
  enddo

  end subroutine multiply_accel_acoustic



!-------------------------------------------------------------------------------------------------
!
! interpolated source arrays
!
!-------------------------------------------------------------------------------------------------

  subroutine multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)

  use constants

  implicit none

  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer :: k,l,m

  ! local parameters
  integer :: ir,it,iv

  ! initializes
  sourcearrayd(:,k,l,m) = ZERO

  do iv = 1,NGLLZ
    do it = 1,NGLLY
      do ir = 1,NGLLX

        sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                           *(G11(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                           +G12(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                           +G13(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

        sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                           *(G21(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                           +G22(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                           +G23(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

        sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                           *(G31(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                           +G32(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                           +G33(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

      enddo
    enddo
  enddo

  end subroutine multiply_arrays_source

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

