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


  subroutine rotate_tensor_Love_to_global(theta,phi, &
                                          A,C,N,L,F, &
                                          c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                          c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! rotates from local radial symmetry given by Love parameterization
! to global global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(in) :: A,C,N,L,F

  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  ! Anderson & Dziewonski, 1982: "Upper mantle anisotropy: evidence from free oscillations", GJR
  ! A = rho * vph**2
  ! C = rho * vpv**2
  ! N = rho * vsh**2
  ! L = rho * vsv**2
  ! F = eta * (A - 2*L)
  !
  ! and therefore (assuming radial axis symmetry)
  ! C11 = A = rho * vph**2
  ! C33 = C = rho * vpv**2
  ! C44 = L = rho * vsv**2
  ! C13 = F = eta * (A - 2*L)
  ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
  ! C22 = C11
  ! C23 = C13
  ! C55 = C44
  ! C66 = N = rho * vsh**2 = (C11-C12)/2

  ! sets up local (radial) coordinate system
  d11 = A
  d12 = A - 2.d0*N
  d13 = F
  d14 = 0.d0
  d15 = 0.d0
  d16 = 0.d0
  d22 = A
  d23 = F
  d24 = 0.d0
  d25 = 0.d0
  d26 = 0.d0
  d33 = C
  d34 = 0.d0
  d35 = 0.d0
  d36 = 0.d0
  d44 = L
  d45 = 0.d0
  d46 = 0.d0
  d55 = L
  d56 = 0.d0
  d66 = N

  ! rotates to global reference system
  call rotate_tensor_radial_to_global(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine rotate_tensor_Love_to_global

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_tensor_azimuthal_to_global(theta,phi, &
                                               A,C,N,L,F, &
                                               Gc,Gs, &
                                               c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                               c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! rotates from local azimuthal symmetry given by Love & Gc,Gs parameterization
! to global global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(in) :: A,C,N,L,F
  double precision,intent(in) :: Gc,Gs

  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  ! only A,C,L,N,F and Gc,Gs are non-zero
  d11 = A
  d12 = A - 2.d0*N
  d13 = F
  d14 = 0.d0
  d15 = 0.d0
  d16 = 0.d0
  d22 = A
  d23 = F
  d24 = 0.d0
  d25 = 0.d0
  d26 = 0.d0
  d33 = C
  d34 = 0.d0
  d35 = 0.d0
  d36 = 0.d0
  d44 = L - Gc
  d45 = -Gs
  d46 = 0.d0
  d55 = L + Gc
  d56 = 0.d0
  d66 = N

  ! rotates to global reference system
  call rotate_tensor_radial_to_global(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine rotate_tensor_azimuthal_to_global


!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_tensor_aniso_to_global(theta,phi, &
                                               A,C,N,L,F, &
                                               Gc,Gs, &
                                               Jc,Js,Kc,Ks,Mc,Ms,Bc,Bs,Hc,Hs,Dc,Ds,Ec,Es, &
                                               c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                               c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! rotates from local symmetry given by aniso (Montagner) parameterization
! to global global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(in) :: A,C,N,L,F
  double precision,intent(in) :: Gc,Gs
  double precision,intent(in) :: Jc,Js,Kc,Ks,Mc,Ms,Bc,Bs,Hc,Hs,Dc,Ds,Ec,Es

  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  ! general case (see also model_aniso_mantle.f90)
  d11 = A + Ec + Bc
  d12 = A - 2.d0 * N - Ec
  d13 = F + Hc
  d14 = Ds + 2.d0 * Js + 2.d0 * Ms
  d15 = 2.d0 * Jc + Dc
  d16 = -0.5d0 * Bs - Es
  d22 = A + Ec - Bc
  d23 = F - Hc
  d24 = 2.d0 * Js - Ds
  d25 = 2.d0 * Jc - 2.d0 * Mc - Dc
  d26 = -Bs/2.d0 + Es

  d33 = C
  d34 = 2.d0 * (Js - Ks)
  d35 = 2.d0 * (Jc - Kc)
  d36 = -Hs

  d44 = L - Gc
  d45 = -Gs
  d46 = Mc - Dc

  d55 = L + Gc
  d56 = Ds - Ms

  d66 = N - Ec

  ! rotates to global reference system
  call rotate_tensor_radial_to_global(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine rotate_tensor_aniso_to_global

!
!-------------------------------------------------------------------------------------------------
!

 subroutine rotate_tensor_radial_to_global(theta,phi, &
                                           d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                           d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! rotates from local (d_ij) to global (c_ij) anisotropic parameters.
! The c_ij are the coefficients in the global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(in) :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                 d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  double precision,intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision,dimension(3,3) :: rotmat
  double precision,dimension(6,6) :: cij,dij
  integer :: i,j
  ! original
  double precision :: costheta,sintheta,cosphi,sinphi
  double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
  double precision :: costwotheta,sintwotheta,costwophi,sintwophi
  double precision :: cosfourtheta,sinfourtheta
  double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
  double precision :: sintwophisq,sintwothetasq

  ! by default, we use Bond's matrix approach which is also used for the backward rotation
  logical, parameter :: USE_BOND_MATRIX_ROTATION = .true.

  ! tensor rotation
  if (USE_BOND_MATRIX_ROTATION) then
    ! rotation matrix
    ! rotates from local (radial) to pole (Cartesian) reference
    ! First column
    rotmat(1,1) =  cos(phi) * cos(theta)
    rotmat(2,1) =  sin(phi) * cos(theta)
    rotmat(3,1) = -sin(theta)

    ! Second column
    rotmat(1,2) = -sin(phi)
    rotmat(2,2) =  cos(phi)
    rotmat(3,2) =  0.d0

    ! Third column
    rotmat(1,3) =  cos(phi)*sin(theta)
    rotmat(2,3) =  sin(phi)*sin(theta)
    rotmat(3,3) =  cos(theta)

    ! Cij Voigt notation
    dij(1,1) = d11
    dij(1,2) = d12
    dij(1,3) = d13
    dij(1,4) = d14
    dij(1,5) = d15
    dij(1,6) = d16
    dij(2,2) = d22
    dij(2,3) = d23
    dij(2,4) = d24
    dij(2,5) = d25
    dij(2,6) = d26
    dij(3,3) = d33
    dij(3,4) = d34
    dij(3,5) = d35
    dij(3,6) = d36
    dij(4,4) = d44
    dij(4,5) = d45
    dij(4,6) = d46
    dij(5,5) = d55
    dij(5,6) = d56
    dij(6,6) = d66
    ! fills lower-triangle, for example: C(2,1) = C(1,2) <->  c21 = c12
    !                                    C(3,1) = C(1,3) <->  c31 = c13
    !                                    C(3,2) = C(2,3) <->  c32 = c23
    do j = 2,6
      do i = 1,j - 1
        dij(j,i) = dij(i,j)
      enddo
    enddo

    call rotate_tensor(dij,rotmat,cij)

    ! returns global cij (SPECFEM x/y/z reference)
    c11 = cij(1,1)
    c12 = cij(1,2)
    c13 = cij(1,3)
    c14 = cij(1,4)
    c15 = cij(1,5)
    c16 = cij(1,6)
    c22 = cij(2,2)
    c23 = cij(2,3)
    c24 = cij(2,4)
    c25 = cij(2,5)
    c26 = cij(2,6)
    c33 = cij(3,3)
    c34 = cij(3,4)
    c35 = cij(3,5)
    c36 = cij(3,6)
    c44 = cij(4,4)
    c45 = cij(4,5)
    c46 = cij(4,6)
    c55 = cij(5,5)
    c56 = cij(5,6)
    c66 = cij(6,6)

  else
    ! original routine... kept for reference
    !
    ! based on original routine rotate_aniso_tensor() from model_aniso_mantle.f90
    costheta = dcos(theta)
    sintheta = dsin(theta)
    cosphi = dcos(phi)
    sinphi = dsin(phi)

    costhetasq = costheta * costheta
    sinthetasq = sintheta * sintheta
    cosphisq = cosphi * cosphi
    sinphisq = sinphi * sinphi

    costhetafour = costhetasq * costhetasq
    sinthetafour = sinthetasq * sinthetasq
    cosphifour = cosphisq * cosphisq
    sinphifour = sinphisq * sinphisq

    costwotheta = dcos(2.d0 * theta)
    sintwotheta = dsin(2.d0 * theta)
    costwophi = dcos(2.d0 * phi)
    sintwophi = dsin(2.d0 * phi)

    cosfourtheta = dcos(4.d0 * theta)
    sinfourtheta = dsin(4.d0 * theta)
    sintwothetasq = sintwotheta * sintwotheta
    sintwophisq = sintwophi * sintwophi

    ! recompute 21 anisotropic coefficients for full anisotropic model using Mathematica

    c11 = d22 * sinphifour - 2.d0 * sintwophi*sinphisq*(d26 * costheta + d24 * sintheta) - &
          2.d0 * cosphisq*sintwophi * (d16 * costhetasq * costheta + &
          (d14 + 2.d0 * d56) * costhetasq*sintheta + &
          (d36 + 2.d0 * d45) * costheta*sinthetasq + d34 * sintheta*sinthetasq) + &
          cosphifour*(d11 * costhetafour + 2.d0 * d15 * costhetasq*sintwotheta + &
          (d13 + 2.d0 * d55)*sintwothetasq * 0.5d0 + &
          2.d0 * d35 * sintwotheta*sinthetasq + d33 * sinthetafour) + &
          (sintwophisq * 0.25d0)*(d12 + d23 + 2.d0 * (d44 + d66) + &
          (d12 - d23 - 2.d0 * d44 + 2.d0 * d66)*costwotheta + &
          2.d0 * (d25 + 2.d0 * d46)*sintwotheta)

    c12 = -((sintwophi * 0.5d0)*sinphisq*((3.d0 * d16 - 4.d0 * d26 + d36 + 2.d0 * d45)*costheta + &
          (d16 - d36 - 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (d14 - 2.d0 * d24 + d34 + 2.d0 * d56 + &
          (d14 - d34 + 2.d0 * d56)*costwotheta)*sintheta)) * 0.5d0 + &
          cosphisq*sintwophi*(d16*costhetasq*costheta - d24*sintheta + &
          (d14 + 2.d0 * d56)*costhetasq*sintheta + d34*sintheta*sinthetasq + &
          costheta*(-d26 + (d36 + 2.d0 * d45)*sinthetasq)) + &
          (sintwophisq * 0.25d0)*(d22 + d11*costhetafour + &
          2.d0 * d15*costhetasq*sintwotheta - 4.d0 * d44*sinthetasq + &
          d33*sinthetafour + costhetasq*(-4.d0 * d66 + &
          2.d0 * (d13 + 2.d0 * d55)*sinthetasq) + &
          costheta*(-8.d0 * d46*sintheta + 4.d0 * d35*sintheta*sinthetasq)) + &
          (cosphifour + sinphifour)*(d12*costhetasq + &
          d23*sinthetasq + d25*sintwotheta)

    c13 = sinphisq*(d23*costhetasq - d25*sintwotheta + d12*sinthetasq) - &
          sintwophi*(d36*costhetasq*costheta + &
          (d34 - 2.d0 * d56)*costhetasq*sintheta + &
          (d16 - 2.d0 * d45)*costheta*sinthetasq + d14*sintheta*sinthetasq) + &
          (cosphisq*(d11 + 6.*d13 + d33 - 4.d0 * d55 - &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*cosfourtheta + &
          4.d0 * (-d15 + d35)*sinfourtheta)) * 0.125d0

    c14 = (-4.d0 * cosphi*sinphisq*((-d14 - 2.d0 * d24 + d34 + 2.d0 * d56)*costheta + &
          (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (-d16 + d26 + d36 + (-d16 + d36 + 2.d0 * d45)*costwotheta)*sintheta) + &
          8.d0 * cosphisq*cosphi*(d14*costhetasq*costheta - &
          (d16 - 2.d0 * d45)*costhetasq*sintheta + &
          (d34 - 2.d0 * d56)*costheta*sinthetasq - d36*sintheta*sinthetasq) + &
          4.d0 * sinphi*sinphisq*(2.d0 * d25*costwotheta + (-d12 + d23)*sintwotheta) + &
          cosphisq*sinphi*(4.d0 * (d15 + d35 - 4*d46)*costwotheta + &
          4.d0 * (d15 - d35)*cosfourtheta - &
          2.d0 * (d11 - d33 + 4.d0 * d44 - 4.d0 * d66 + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*costwotheta)*sintwotheta)) * 0.125d0

    c15 = (8.d0 * sinphi*sinphisq*(-(d24*costheta) + d26*sintheta) + &
          4.d0 * cosphi*sinphisq*(2.d0 * (d25 + 2.d0 * d46)*costwotheta + &
          (-d12 + d23 + 2.d0 * d44 - 2.d0 * d66)*sintwotheta) + &
          cosphisq*cosphi*(4.d0 * (d15 + d35)*costwotheta + &
          4.d0 * (d15 - d35)*cosfourtheta - 2.d0 * (d11 - d33 + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*costwotheta)*sintwotheta) - &
          2.d0 * cosphisq*sinphi*((d14 + 3.d0 * d34 + 2.d0 * d56)*costheta + &
          3.d0 * (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) - &
          (3.d0 * d16 + d36 + 2.d0 * d45)*sintheta + &
          3.d0 * (-d16 + d36 + 2.d0 * d45)*(-4.d0 * sinthetasq*sintheta + 3.d0 * sintheta))) * 0.125d0

    c16 = -(sinphifour*(d26*costheta + d24*sintheta)) - &
          (3.d0 * (sintwophisq * 0.25d0)*((3.d0 * d16 - 4.d0 * d26 + d36 + 2.d0 * d45)*costheta + &
          (d16 - d36 - 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (d14 - 2.d0 * d24 + d34 + 2.d0 * d56 + &
          (d14 - d34 + 2.d0 * d56)*costwotheta)*sintheta)) * 0.25d0 + &
          cosphifour*(d16*costhetasq*costheta + &
          (d14 + 2.d0 * d56)*costhetasq*sintheta + &
          (d36 + 2.d0 * d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
          (sintwophi * 0.5d0)*sinphisq*(-d22 + (d12 + 2.d0 * d66)*costhetasq + &
          2.d0 * d46*sintwotheta + (d23 + 2.d0 * d44)*sinthetasq + d25*sintwotheta) + &
          cosphisq*(sintwophi * 0.5d0)*(d11*costhetafour + &
          2.d0 * d15*costhetasq*sintwotheta - (d23 + 2.d0 * d44)*sinthetasq + &
          d33*sinthetafour - costhetasq*(d12 + &
          2.d0 * d66 - 2.d0 * (d13 + 2.d0 * d55)*sinthetasq) - &
          (d25 - d35 + 2.d0 * d46 + d35*costwotheta)*sintwotheta)

    c22 = d22*cosphifour + 2.d0 * cosphisq*sintwophi*(d26*costheta + d24*sintheta) + &
          2.d0 * sintwophi*sinphisq*(d16*costhetasq*costheta + &
          (d14 + 2.d0 * d56)*costhetasq*sintheta + &
          (d36 + 2.d0 * d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
          sinphifour*(d11*costhetafour + 2.d0 * d15*costhetasq*sintwotheta + &
          (d13 + 2.d0 * d55)*sintwothetasq * 0.5d0 + &
          2.d0 * d35*sintwotheta*sinthetasq + d33*sinthetafour) + &
          (sintwophisq * 0.25d0)*(d12 + d23 + 2.d0 * (d44 + d66) + &
          (d12 - d23 - 2.d0 * d44 + 2.d0 * d66)*costwotheta + &
          2.d0 * (d25 + 2.d0 * d46)*sintwotheta)

    c23 = d13*costhetafour*sinphisq + &
          sintheta*sinthetasq*(d14*sintwophi + d13*sinphisq*sintheta) + &
          costheta*sinthetasq*((d16 - 2.d0 * d45)*sintwophi + &
          2.d0 * (d15 - d35)*sinphisq*sintheta) + &
          costhetasq*costheta*(d36*sintwophi + &
          2.d0 * (-d15 + d35)*sinphisq*sintheta) + &
          costhetasq*sintheta*((d34 - 2.d0 * d56)*sintwophi + &
          (d11 + d33 - 4.d0 * d55)*sinphisq*sintheta) + &
          cosphisq*(d23*costhetasq - d25*sintwotheta + d12*sinthetasq)

    c24 = (8.d0 * cosphisq*cosphi*(d24*costheta - d26*sintheta) + &
          4.d0 * cosphisq*sinphi*(2.d0 * (d25 + 2.d0 * d46)*costwotheta + &
          (-d12 + d23 + 2.d0 * d44 - 2.d0 * d66)*sintwotheta) + &
          sinphi*sinphisq*(4.d0 * (d15 + d35)*costwotheta + &
          4.d0 * (d15 - d35)*cosfourtheta - &
          2.d0 * (d11 - d33 + (d11 - 2.d0 * d13 + &
          d33 - 4.d0 * d55)*costwotheta)*sintwotheta) + &
          2.d0 * cosphi*sinphisq*((d14 + 3.d0 * d34 + 2.d0 * d56)*costheta + &
          3.d0 * (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) - &
          (3.d0 * d16 + d36 + 2.d0 * d45)*sintheta + &
          3.d0 * (-d16 + d36 + 2.d0 * d45)*(-4.d0 * sinthetasq*sintheta + 3.d0 * sintheta))) * 0.125d0

    c25 = (4.d0 * cosphisq*sinphi*((-d14 - 2.d0 * d24 + d34 + 2.d0 * d56)*costheta + &
          (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (-d16 + d26 + d36 + (-d16 + d36 + 2.d0 * d45)*costwotheta)*sintheta) - &
          8.d0 * sinphi*sinphisq*(d14*costhetasq*costheta - &
          (d16 - 2.d0 * d45)*costhetasq*sintheta + &
          (d34 - 2.d0 * d56)*costheta*sinthetasq - d36*sintheta*sinthetasq) + &
          4.d0 * cosphisq*cosphi*(2.d0 * d25*costwotheta + (-d12 + d23)*sintwotheta) + &
          cosphi*sinphisq*(4.d0 * (d15 + d35 - 4.d0 * d46)*costwotheta + &
          4.d0 * (d15 - d35)*cosfourtheta - 2.d0 * (d11 - d33 + 4.d0 * d44 - 4.d0 * d66 + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*costwotheta)*sintwotheta)) * 0.125d0

    c26 = cosphifour*(d26*costheta + d24*sintheta) + &
          (3.d0 * (sintwophisq * 0.25d0)*((3.d0 * d16 - 4.d0 * d26 + d36 + 2.d0 * d45)*costheta + &
          (d16 - d36 - 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (d14 - 2.d0 * d24 + d34 + 2.d0 * d56 + &
          (d14 - d34 + 2.d0 * d56)*costwotheta)*sintheta)) * 0.25d0 - &
          sinphifour*(d16*costhetasq*costheta + &
          (d14 + 2.d0 * d56)*costhetasq*sintheta + &
          (d36 + 2.d0 * d45)*costheta*sinthetasq + d34*sintheta*sinthetasq) + &
          cosphisq*(sintwophi * 0.5d0)*(-d22 + (d12 + 2.d0 * d66)*costhetasq + &
          2.d0 * d46*sintwotheta + (d23 + 2.d0 * d44)*sinthetasq + &
          d25*sintwotheta) + (sintwophi * 0.5d0)*sinphisq*(d11*costhetafour + &
          2.d0 * d15*costhetasq*sintwotheta - (d23 + 2.d0 * d44)*sinthetasq + &
          d33*sinthetafour - costhetasq*(d12 + &
          2.d0 * d66 - 2.d0 * (d13 + 2.d0 * d55)*sinthetasq) - &
          (d25 - d35 + 2.d0 * d46 + d35*costwotheta)*sintwotheta)

    c33 = d33*costhetafour - 2.d0 * d35*costhetasq*sintwotheta + &
          (d13 + 2.d0 * d55)*sintwothetasq * 0.5d0 - &
          2.d0 * d15*sintwotheta*sinthetasq + d11*sinthetafour

    c34 = cosphi*(d34*costhetasq*costheta - (d36 + 2.d0 * d45)*costhetasq*sintheta + &
          (d14 + 2.d0 * d56)*costheta*sinthetasq - d16*sintheta*sinthetasq) + &
          (sinphi*(4.d0 * (d15 + d35)*costwotheta + 4.d0 * (-d15 + d35)*cosfourtheta + &
          2.d0 * (-d11 + d33)*sintwotheta + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*sinfourtheta)) * 0.125d0

    c35 = sinphi*(-(d34*costhetasq*costheta) + &
          (d36 + 2.d0 * d45)*costhetasq*sintheta - &
          (d14 + 2.d0 * d56)*costheta*sinthetasq + d16*sintheta*sinthetasq) + &
          (cosphi*(4.d0 * (d15 + d35)*costwotheta + 4.d0 * (-d15 + d35)*cosfourtheta + &
          2.d0 * (-d11 + d33)*sintwotheta + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*sinfourtheta)) * 0.125d0

    c36 = (4.d0 * costwophi*((d16 + 3.d0 * d36 - 2.d0 * d45)*costheta + &
          (-d16 + d36 + 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          (3.d0 * d14 + d34 - 2.d0 * d56)*sintheta + &
          (-d14 + d34 - 2.d0 * d56)*(-4.d0 * sinthetasq*sintheta + 3.d0 * sintheta)) + &
          sintwophi*(d11 - 4.d0 * d12 + 6.*d13 - 4.d0 * d23 + d33 - 4.d0 * d55 + &
          4.d0 * (d12 - d23)*costwotheta - &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*cosfourtheta + &
          8.d0 * d25*sintwotheta + 4.d0 * (-d15 + d35)*sinfourtheta)) * 0.0625d0

    c44 = (d11 - 2.d0 * d13 + d33 + 4.d0 * (d44 + d55 + d66) - &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * (d44 - d55 + d66))*costwophi + &
          4.d0 * sintwophi*((d16 - d36 + 2.d0 * d45)*costheta + &
          (-d16 + d36 + 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) - &
          2.d0 * (d14 - d34 + (d14 - d34 + 2.d0 * d56)*costwotheta)*sintheta) + &
          8.d0 * cosphisq*((d44 - d66)*costwotheta - 2.d0 * d46*sintwotheta) + &
          2.d0 * sinphisq*(-((d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*cosfourtheta) + &
          4.d0 * (-d15 + d35)*sinfourtheta)) * 0.0625d0

    c45 = (4.d0 * costwophi*((d16 - d36 + 2.d0 * d45)*costheta + &
          (-d16 + d36 + 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) - &
          2.d0 * (d14 - d34 + (d14 - d34 + 2.d0 * d56)*costwotheta)*sintheta) + &
          sintwophi*(d11 - 2.d0 * d13 + d33 - 4.d0 * (d44 - d55 + d66) + &
          4.d0 * (-d44 + d66)*costwotheta - &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*cosfourtheta + 8.d0 * d46*sintwotheta + &
          4.d0 * (-d15 + d35)*sinfourtheta)) * 0.0625d0

    c46 = (-2.d0 * sinphi*sinphisq*((-d14 + d34 + 2.d0 * d56)*costheta + &
          (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (-d16 + d36 + (-d16 + d36 + 2.d0 * d45)*costwotheta)*sintheta) + &
          4.d0 * cosphisq*cosphi*(2.d0 * d46*costwotheta + (d44 - d66)*sintwotheta) + &
          cosphi*sinphisq*(4.d0 * (d15 - 2.d0 * d25 + d35 - 2.d0 * d46)*costwotheta + &
          4.d0 * (d15 - d35)*cosfourtheta - &
          2.d0 * (d11 - 2.d0 * d12 + 2.d0 * d23 - d33 + 2.d0 * d44 - 2.d0 * d66 + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*costwotheta)*sintwotheta) + &
          4.d0 * cosphisq*sinphi*((d14 - 2.d0 * d24 + d34)*costheta + &
          (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) - &
          (d16 - 2.d0 * d26 + d36)*sintheta + &
          (-d16 + d36 + 2.d0 * d45)*(-4.d0 * sinthetasq*sintheta + 3.d0 * sintheta))) * 0.125d0

    c55 = d66 * sinphisq*sinthetasq + (sintwotheta * 0.5d0)*(-2.d0 * d46*sinphisq + &
          (d36 + d45)*sintwophi*sintheta) + &
          costhetasq*(d44 * sinphisq + (d14 + d56)*sintwophi*sintheta) - &
          sintwophi*(d45 * costhetasq*costheta + d34*costhetasq*sintheta + &
          d16 * costheta*sinthetasq + d56*sintheta*sinthetasq) + &
          (cosphisq*(d11 - 2.d0 * d13 + d33 + 4.d0 * d55 - &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*cosfourtheta + &
          4.d0 * (-d15 + d35)*sinfourtheta)) * 0.125d0

    c56 = (8.d0 * cosphisq*cosphi*(d56*costhetasq*costheta - &
          (d16 - d36 - d45)*costhetasq*sintheta - &
          (d14 - d34 + d56)*costheta*sinthetasq - d45*sintheta*sinthetasq) + &
          4.d0 * sinphi*sinphisq*(2.d0 * d46*costwotheta + (d44 - d66)*sintwotheta) + &
          cosphisq*sinphi*(4.d0 * (d15 - 2.d0 * d25 + d35 - 2.d0 * d46)*costwotheta + &
          4.d0 * (d15 - d35)*cosfourtheta - &
          2.d0 * (d11 - 2.d0 * d12 + 2.d0 * d23 - d33 + 2.d0 * d44 - 2.d0 * d66 + &
          (d11 - 2.d0 * d13 + d33 - 4.d0 * d55)*costwotheta)*sintwotheta) - &
          4.d0 * cosphi*sinphisq*((d14 - 2.d0 * d24 + d34)*costheta + &
          (d14 - d34 + 2.d0 * d56)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) - &
          (d16 - 2.d0 * d26 + d36)*sintheta + &
          (-d16 + d36 + 2.d0 * d45)*(-4.d0 * sinthetasq*sintheta + 3.d0 * sintheta))) * 0.125d0

    c66 = -((sintwophi * 0.5d0)*sinphisq*((3.d0 * d16 - 4.d0 * d26 + d36 + 2.d0 * d45)*costheta + &
          (d16 - d36 - 2.d0 * d45)*(4.d0 * costhetasq*costheta - 3.d0 * costheta) + &
          2.d0 * (d14 - 2.d0 * d24 + d34 + 2.d0 * d56 + &
          (d14 - d34 + 2.d0 * d56)*costwotheta)*sintheta)) * 0.5d0 + &
          cosphisq*sintwophi*(d16*costhetasq*costheta - d24*sintheta + &
          (d14 + 2.d0 * d56)*costhetasq*sintheta + d34*sintheta*sinthetasq + &
          costheta*(-d26 + (d36 + 2.d0 * d45)*sinthetasq)) + &
          (sintwophisq * 0.25d0)*(d22 + d11*costhetafour + &
          2.d0 * d15*costhetasq*sintwotheta - 2.d0 * (d23 + d44)*sinthetasq + &
          d33*sinthetafour - 2.d0 * sintwotheta*(d25 + d46 - d35*sinthetasq) - &
          2.d0 * costhetasq*(d12 + d66 - (d13 + 2.d0 * d55)*sinthetasq)) + &
          (cosphifour + sinphifour)*(d66*costhetasq + &
          d44*sinthetasq + d46*sintwotheta)

  endif ! USE_BOND_MATRIX_ROTATION

  end subroutine rotate_tensor_radial_to_global


!
!-------------------------------------------------------------------------------------------------
!


  subroutine rotate_tensor_tiso_to_cij(theta_in,phi_in, &
                                       rhovphsq,rhovpvsq,rhovsvsq,rhovshsq,eta_aniso, &
                                       c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                       c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use constants, only: CUSTOM_REAL
  implicit none

! note: CUSTOM_REAL parameters, not double precision like other routines

  real(kind=CUSTOM_REAL),intent(in) :: theta_in,phi_in
  real(kind=CUSTOM_REAL),intent(in) :: rhovphsq,rhovpvsq,rhovsvsq,rhovshsq,eta_aniso

  ! aniso element
  real(kind=CUSTOM_REAL),intent(out) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                        c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: theta,phi
  double precision :: A,C,N,L,F
  double precision :: t11,t12,t13,t14,t15,t16,t22,t23,t24,t25,t26, &
                      t33,t34,t35,t36,t44,t45,t46,t55,t56,t66

  ! Love parameterization
  A = rhovphsq
  C = rhovpvsq
  L = rhovsvsq
  N = rhovshsq
  F = eta_aniso * (A - 2.d0 * L)

  theta = theta_in
  phi = phi_in

  call rotate_tensor_Love_to_global(theta,phi, &
                                    A,C,N,L,F, &
                                    t11,t12,t13,t14,t15,t16,t22,t23,t24,t25,t26, &
                                    t33,t34,t35,t36,t44,t45,t46,t55,t56,t66)

  ! converts and returns back in custom real
  c11 = real(t11,kind=CUSTOM_REAL)
  c12 = real(t12,kind=CUSTOM_REAL)
  c13 = real(t13,kind=CUSTOM_REAL)
  c14 = real(t14,kind=CUSTOM_REAL)
  c15 = real(t15,kind=CUSTOM_REAL)
  c16 = real(t16,kind=CUSTOM_REAL)
  c22 = real(t22,kind=CUSTOM_REAL)
  c23 = real(t23,kind=CUSTOM_REAL)
  c24 = real(t24,kind=CUSTOM_REAL)
  c25 = real(t25,kind=CUSTOM_REAL)
  c26 = real(t26,kind=CUSTOM_REAL)
  c33 = real(t33,kind=CUSTOM_REAL)
  c34 = real(t34,kind=CUSTOM_REAL)
  c35 = real(t35,kind=CUSTOM_REAL)
  c36 = real(t36,kind=CUSTOM_REAL)
  c44 = real(t44,kind=CUSTOM_REAL)
  c45 = real(t45,kind=CUSTOM_REAL)
  c46 = real(t46,kind=CUSTOM_REAL)
  c55 = real(t55,kind=CUSTOM_REAL)
  c56 = real(t56,kind=CUSTOM_REAL)
  c66 = real(t66,kind=CUSTOM_REAL)


! original routine... kept for reference
!
! note: the routine has somewhat lower precision, we thus use the (full) rotate Love routine
!       since this will be called in a preparation step and shouldn't increase the time loop costs
!
!  ! tiso elements
!  double precision :: sinphifour,cosphisq,sinphisq,costhetasq,sinthetasq, &
!        cosphifour,costhetafour,sinthetafour,cosfourphi, &
!        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
!        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq
!
!  double precision :: two_L,two_N
!  double precision :: four_L,four_N
!
!  double precision :: twoetaminone,etaminone
!  double precision :: two_eta_aniso,four_eta_aniso,six_eta_aniso
!
!  double precision :: templ1,templ1_cos,templ2,templ2_cos,templ3,templ3_two,templ3_cos
!
!  ! precompute some products to reduce the CPU time
!
!  costheta = cos(theta)
!  sintheta = sin(theta)
!  cosphi = cos(phi)
!  sinphi = sin(phi)
!
!  costhetasq = costheta * costheta
!  sinthetasq = sintheta * sintheta
!  cosphisq = cosphi * cosphi
!  sinphisq = sinphi * sinphi
!
!  costhetafour = costhetasq * costhetasq
!  sinthetafour = sinthetasq * sinthetasq
!  cosphifour = cosphisq * cosphisq
!  sinphifour = sinphisq * sinphisq
!
!  costwotheta = cos(2.d0 * theta)
!  sintwotheta = sin(2.d0 * theta)
!  costwophi = cos(2.d0 * phi)
!  sintwophi = sin(2.d0 * phi)
!
!  cosfourtheta = cos(4.d0 * theta)
!  cosfourphi = cos(4.d0 * phi)
!
!  costwothetasq = costwotheta * costwotheta
!
!  costwophisq = costwophi * costwophi
!  sintwophisq = sintwophi * sintwophi
!
!  etaminone = eta_aniso - 1.d0
!  twoetaminone = 2.d0 * eta_aniso - 1.d0
!
!  ! precompute some products to reduce the CPU time
!  two_eta_aniso = 2.d0 * eta_aniso
!  four_eta_aniso = 4.d0 * eta_aniso
!  six_eta_aniso = 6.d0 * eta_aniso
!
!  two_L = 2.d0 * L
!  two_N = 2.d0 * N
!
!  four_L = 4.d0 * L
!  four_N = 4.d0 * N
!
!  ! pre-compute temporary values
!  templ1 = four_L - C + twoetaminone*A - four_eta_aniso*L
!  templ1_cos = A - C + costwotheta*templ1
!  templ2 = four_L - C - A + two_eta_aniso*A - four_eta_aniso*L
!  templ2_cos = C - A + costwotheta*templ2
!  templ3 = A + C - two_eta_aniso*A + four_eta_aniso*L
!  templ3_two = templ3 - two_N - two_L
!  templ3_cos = templ3_two + costwotheta*templ2
!
!  ! reordering operations to facilitate compilation, avoiding divisions, using locality for temporary values
!  c11 = A*sinphifour &
!        + 2.d0*cosphisq*sinphisq* &
!        ( A*costhetasq + sinthetasq*(eta_aniso*A + two_L - two_eta_aniso*L) ) &
!        + cosphifour*(A*costhetafour &
!          + 2.d0*costhetasq*sinthetasq*(eta_aniso*A + two_L - two_eta_aniso*L) &
!          + C*sinthetafour)
!
!  c12 = 0.25d0*costhetasq &
!        *(A - two_N)*(3.d0 + cosfourphi) &
!        - four_N*cosphisq*costhetasq*sinphisq &
!        + 0.03125d0*A*sintwophisq*(11.d0 + cosfourtheta + 4.0*costwotheta) &
!        + eta_aniso*sinthetasq*(A - two_L) &
!                   *(cosphifour + sinphifour + 2.d0*cosphisq*costhetasq*sinphisq) &
!        + C*cosphisq*sinphisq*sinthetafour &
!        - L*sintwophisq*sinthetafour
!
!  c13 = 0.125d0*cosphisq &
!        *(A + six_eta_aniso*A + C - four_L &
!              - 12.d0*eta_aniso*L + cosfourtheta*templ1) &
!        + sinphisq*(eta_aniso*costhetasq*(A - two_L) + sinthetasq*(A - two_N))
!
!  ! uses temporary templ1 from c13
!  c15 = cosphi*costheta*sintheta* &
!        ( 0.5d0*cosphisq* (C - A + costwotheta*templ1) &
!          + etaminone*sinphisq*(A - two_L))
!
!  c14 = costheta*sinphi*sintheta* &
!        ( 0.5d0*cosphisq*(templ2_cos + four_N - four_L) &
!          + sinphisq*(etaminone*A + 2.d0*(N - eta_aniso*L)) )
!
!  ! uses temporary templ2_cos from c14
!  c16 = 0.5d0*cosphi*sinphi*sinthetasq* &
!        ( cosphisq*templ2_cos &
!          + 2.d0*etaminone*sinphisq*(A - two_L) )
!
!  c22 = A*cosphifour + 2.d0*cosphisq*sinphisq* &
!        (A*costhetasq + sinthetasq*(eta_aniso*A + two_L - two_eta_aniso*L)) &
!        + sinphifour* &
!        (A*costhetafour + 2.d0*costhetasq*sinthetasq*(eta_aniso*A &
!              + two_L - two_eta_aniso*L) + C*sinthetafour)
!
!  ! uses temporary templ1 from c13
!  c23 = 0.125d0*sinphisq*(A + six_eta_aniso*A &
!          + C - four_L - 12.d0*eta_aniso*L + cosfourtheta*templ1) &
!        + cosphisq*(eta_aniso*costhetasq*(A - two_L) + sinthetasq*(A - two_N))
!
!  ! uses temporary templ1 from c13
!  c24 = costheta*sinphi*sintheta* &
!        ( etaminone*cosphisq*(A - two_L) &
!          + 0.5d0*sinphisq*(C - A + costwotheta*templ1) )
!
!  ! uses temporary templ2_cos from c14
!  c25 = cosphi*costheta*sintheta* &
!        ( cosphisq*(etaminone*A + 2.d0*(N - eta_aniso*L)) &
!          + 0.5d0*sinphisq*(templ2_cos + four_N - four_L) )
!
!  ! uses temporary templ2_cos from c14
!  c26 = 0.5d0*cosphi*sinphi*sinthetasq* &
!        ( 2.d0*etaminone*cosphisq*(A - two_L) &
!          + sinphisq*templ2_cos )
!
!  c33 = C*costhetafour &
!        + 2.d0*costhetasq*sinthetasq*(two_L + eta_aniso*(A - two_L)) &
!        + A*sinthetafour
!
!  ! uses temporary templ1_cos from c13
!  c34 = - 0.25d0*sinphi*sintwotheta*templ1_cos
!
!  ! uses temporary templ1_cos from c34
!  c35 = - 0.25d0*cosphi*sintwotheta*templ1_cos
!
!  ! uses temporary templ1_cos from c34
!  c36 = - 0.25d0*sintwophi*sinthetasq &
!        *(templ1_cos - four_N + four_L)
!
!  c44 = cosphisq*(L*costhetasq + N*sinthetasq) &
!        + sinphisq*(L*costwothetasq + costhetasq*sinthetasq*templ3)
!
!  ! uses temporary templ3 from c44
!  c46 = - cosphi*costheta*sintheta* &
!          ( cosphisq*(N - L) - 0.5d0*sinphisq*templ3_cos  )
!
!  ! uses templ3 from c46
!  c45 = 0.25d0*sintwophi*sinthetasq* &
!        (templ3_two + costwotheta*(A + C - two_eta_aniso*A + 4.d0*etaminone*L))
!
!  c55 = sinphisq*(L*costhetasq + N*sinthetasq) &
!        + cosphisq*(L*costwothetasq &
!            + costhetasq*sinthetasq*(A - two_eta_aniso*A + C + four_eta_aniso*L) )
!
!  ! uses temporary templ3_cos from c46
!  c56 = costheta*sinphi*sintheta* &
!        ( 0.5d0*cosphisq*templ3_cos + sinphisq*(L - N) )
!
!  c66 = N*costwophisq*costhetasq &
!        - 2.d0*cosphisq*costhetasq*sinphisq*(A - two_N) &
!        + 0.03125d0*A*sintwophisq*(11.d0 + 4.d0*costwotheta + cosfourtheta) &
!        - 0.125d0*L*sinthetasq* &
!        ( -6.d0 - 2.d0*costwotheta - 2.d0*cosfourphi &
!                  + cos(4.d0*phi - 2.d0*theta) &
!                  + cos(2.d0*(2.d0*phi + theta)) ) &
!        + C*cosphisq*sinphisq*sinthetafour &
!        - 0.5d0*eta_aniso*sintwophisq*sinthetafour*(A - two_L)

  end subroutine rotate_tensor_tiso_to_cij



!
!-------------------------------------------------------------------------------------------------
!
! opposite direction (SPECFEM reference system to local radial system)
!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_tensor_global_to_Love(theta,phi, &
                                          A,C,N,L,F, &
                                          c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                          c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! rotates from local radial symmetry given by Love parameterization
! to global global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: A,C,N,L,F

  double precision,intent(in) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  ! Anderson & Dziewonski, 1982: "Upper mantle anisotropy: evidence from free oscillations", GJR
  ! A = rho * vph**2
  ! C = rho * vpv**2
  ! N = rho * vsh**2
  ! L = rho * vsv**2
  ! F = eta * (A - 2*L)
  !
  ! and therefore (assuming radial axis symmetry)
  ! C11 = A = rho * vph**2
  ! C33 = C = rho * vpv**2
  ! C44 = L = rho * vsv**2
  ! C13 = F = eta * (A - 2*L)
  ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
  ! C22 = C11
  ! C23 = C13
  ! C55 = C44
  ! C66 = N = rho * vsh**2 = (C11-C12)/2

  ! rotates from global to local (radial) reference system
  call rotate_tensor_global_to_radial(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! sets up Love parameters from local (radial) system
  ! for pure tiso tensors:
  ! A = d22
  ! F = d23
  ! C = d33
  ! L = d44
  ! N = d66
  !
  ! for more general tensors:
  !
  ! Sieminski, 2007:
  ! A = 1/8 (3 C11 + 3 C22 + 2 C12 + 4 C66)
  ! C = C33
  ! N = 1/8 (C11 + C22 - 2 C12 + 4 C66)
  ! L = 1/2 (C44 + C55)
  ! F = 1/2 (C13 + C23)
  ! and derived eta = F / (A - 2 L)
  !
  A = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)
  C = d33
  N = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
  L = 0.5d0 * (d44 + d55)
  F = 0.5d0 * (d13 + d23)

  end subroutine rotate_tensor_global_to_Love


!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_tensor_global_to_azi(theta,phi, &
                                          A,C,N,L,F, &
                                          Gc,Gs, &
                                          c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                          c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! rotates from local radial symmetry given by Love parameterization
! to global global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: A,C,N,L,F,Gc,Gs

  double precision,intent(in) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  ! Anderson & Dziewonski, 1982: "Upper mantle anisotropy: evidence from free oscillations", GJR
  ! A = rho * vph**2
  ! C = rho * vpv**2
  ! N = rho * vsh**2
  ! L = rho * vsv**2
  ! F = eta * (A - 2*L)
  !
  ! and therefore (assuming radial axis symmetry)
  ! C11 = A = rho * vph**2
  ! C33 = C = rho * vpv**2
  ! C44 = L = rho * vsv**2
  ! C13 = F = eta * (A - 2*L)
  ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
  ! C22 = C11
  ! C23 = C13
  ! C55 = C44
  ! C66 = N = rho * vsh**2 = (C11-C12)/2

  ! rotates from global to local (radial) reference system
  call rotate_tensor_global_to_radial(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! sets up Love parameters from local (radial) system
  ! for pure tiso tensors:
  ! A = d22
  ! F = d23
  ! C = d33
  ! L = d44
  ! N = d66
  !
  ! for more general tensors:
  !
  ! Sieminski, 2007:
  ! A = 1/8 (3 C11 + 3 C22 + 2 C12 + 4 C66)
  ! C = C33
  ! N = 1/8 (C11 + C22 - 2 C12 + 4 C66)
  ! L = 1/2 (C44 + C55)
  ! F = 1/2 (C13 + C23)
  ! and derived eta = F / (A - 2 L)
  !
  A = 0.125d0 * (3.d0 * d11 + 3.d0 * d22 + 2.d0 * d12 + 4.d0 * d66)
  C = d33
  N = 0.125d0 * (d11 + d22 - 2.d0 * d12 + 4.d0 * d66)
  L = 0.5d0 * (d44 + d55)
  F = 0.5d0 * (d13 + d23)

  ! see rotate for aniso rotations:
  ! in local (radial) symmetry
  !   d44 = L - Gc
  !   d45 = -Gs
  !   d55 = L + Gc
  !
  ! therefor: Gs = - d45
  !           Gc = 1/2 ( d55 - d44 )
  Gs = - d45
  Gc = 0.5d0 * (d55 - d44)

  end subroutine rotate_tensor_global_to_azi


!
!-------------------------------------------------------------------------------------------------
!


  subroutine rotate_tensor_global_to_radial_vector(cij_kl,cij_kl_spherical,theta_in,phi_in)

! Purpose : compute the kernels in r,theta,phi (cij_kl_spherical)
! from the kernels in x,y,z (cij_kl) (x,y,z to r,theta,phi)
! At r,theta,phi fixed
! theta and phi are in radians

! Coeff from Min's routine rotate_anisotropic_tensor
! with the help of Collect[Expand[cij],{dij}] in Mathematica

! Definition of the output array cij_kl_spherical :
! cij_kl_spherical(1) = C11 ; cij_kl_spherical(2) = C12 ; cij_kl_spherical(3) = C13
! cij_kl_spherical(4) = C14 ; cij_kl_spherical(5) = C15 ; cij_kl_spherical(6) = C16
! cij_kl_spherical(7) = C22 ; cij_kl_spherical(8) = C23 ; cij_kl_spherical(9) = C24
! cij_kl_spherical(10) = C25 ; cij_kl_spherical(11) = C26 ; cij_kl_spherical(12) = C33
! cij_kl_spherical(13) = C34 ; cij_kl_spherical(14) = C35 ; cij_kl_spherical(15) = C36
! cij_kl_spherical(16) = C44 ; cij_kl_spherical(17) = C45 ; cij_kl_spherical(18) = C46
! cij_kl_spherical(19) = C55 ; cij_kl_spherical(20) = C56 ; cij_kl_spherical(21) = C66
!
! where the Cij (Voigt's notation) are defined as function of
! the components of the elastic tensor in spherical coordinates
! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(21), intent(in) :: cij_kl
  real(kind=CUSTOM_REAL), dimension(21), intent(out) :: cij_kl_spherical

  real(kind=CUSTOM_REAL), intent(in) :: theta_in,phi_in

  ! local parameters
  double precision :: theta,phi
  double precision :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66


  theta = dble(theta_in)
  phi = dble(phi_in)

  ! input SPECFEM reference
  c11 = cij_kl(1)
  c12 = cij_kl(2)
  c13 = cij_kl(3)
  c14 = cij_kl(4)
  c15 = cij_kl(5)
  c16 = cij_kl(6)
  c22 = cij_kl(7)
  c23 = cij_kl(8)
  c24 = cij_kl(9)
  c25 = cij_kl(10)
  c26 = cij_kl(11)
  c33 = cij_kl(12)
  c34 = cij_kl(13)
  c35 = cij_kl(14)
  c36 = cij_kl(15)
  c44 = cij_kl(16)
  c45 = cij_kl(17)
  c46 = cij_kl(18)
  c55 = cij_kl(19)
  c56 = cij_kl(20)
  c66 = cij_kl(21)

  call rotate_tensor_global_to_radial(theta,phi, &
                                      d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                      d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                      c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! return radial reference
  cij_kl_spherical(1) = real(d11,kind=CUSTOM_REAL)
  cij_kl_spherical(2) = real(d12,kind=CUSTOM_REAL)
  cij_kl_spherical(3) = real(d13,kind=CUSTOM_REAL)
  cij_kl_spherical(4) = real(d14,kind=CUSTOM_REAL)
  cij_kl_spherical(5) = real(d15,kind=CUSTOM_REAL)
  cij_kl_spherical(6) = real(d16,kind=CUSTOM_REAL)
  cij_kl_spherical(7) = real(d22,kind=CUSTOM_REAL)
  cij_kl_spherical(8) = real(d23,kind=CUSTOM_REAL)
  cij_kl_spherical(9) = real(d24,kind=CUSTOM_REAL)
  cij_kl_spherical(10) = real(d25,kind=CUSTOM_REAL)
  cij_kl_spherical(11) = real(d26,kind=CUSTOM_REAL)
  cij_kl_spherical(12) = real(d33,kind=CUSTOM_REAL)
  cij_kl_spherical(13) = real(d34,kind=CUSTOM_REAL)
  cij_kl_spherical(14) = real(d35,kind=CUSTOM_REAL)
  cij_kl_spherical(15) = real(d36,kind=CUSTOM_REAL)
  cij_kl_spherical(16) = real(d44,kind=CUSTOM_REAL)
  cij_kl_spherical(17) = real(d45,kind=CUSTOM_REAL)
  cij_kl_spherical(18) = real(d46,kind=CUSTOM_REAL)
  cij_kl_spherical(19) = real(d55,kind=CUSTOM_REAL)
  cij_kl_spherical(20) = real(d56,kind=CUSTOM_REAL)
  cij_kl_spherical(21) = real(d66,kind=CUSTOM_REAL)


! for reference, original routine...
!
!  double precision :: costheta,sintheta,cosphi,sinphi
!  double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
!  double precision :: costwotheta,sintwotheta,costwophi,sintwophi
!  double precision :: cosfourtheta,sinfourtheta,cosfourphi,sinfourphi
!  double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
!  double precision :: sintwophisq,sintwothetasq
!  double precision :: costhreetheta,sinthreetheta,costhreephi,sinthreephi
!
!  costheta = dcos(theta)
!  sintheta = dsin(theta)
!  cosphi = dcos(phi)
!  sinphi = dsin(phi)
!
!  costhetasq = costheta * costheta
!  sinthetasq = sintheta * sintheta
!  cosphisq = cosphi * cosphi
!  sinphisq = sinphi * sinphi
!
!  costhetafour = costhetasq * costhetasq
!  sinthetafour = sinthetasq * sinthetasq
!  cosphifour = cosphisq * cosphisq
!  sinphifour = sinphisq * sinphisq
!
!  costwotheta = dcos(2.d0*theta)
!  sintwotheta = dsin(2.d0*theta)
!  costwophi = dcos(2.d0*phi)
!  sintwophi = dsin(2.d0*phi)
!
!  costhreetheta=dcos(3.d0*theta)
!  sinthreetheta=dsin(3.d0*theta)
!  costhreephi=dcos(3.d0*phi)
!  sinthreephi=dsin(3.d0*phi)
!
!  cosfourtheta = dcos(4.d0*theta)
!  sinfourtheta = dsin(4.d0*theta)
!  cosfourphi = dcos(4.d0*phi)
!  sinfourphi = dsin(4.d0*phi)
!  sintwothetasq = sintwotheta * sintwotheta
!  sintwophisq = sintwophi * sintwophi
!
! cij_kl_spherical(1) = ONE_SIXTEENTH * (cij_kl(16) - cij_kl(16)* costwophi + &
!     16.d0* cosphi*cosphisq* costhetafour* (cij_kl(1)* cosphi + cij_kl(6)* sinphi) + &
!     2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq - &
!     2.d0* (cij_kl(16)* cosfourtheta* sinphisq + &
!     2.d0* costhetafour* (-4* cij_kl(7)* sinphifour - &
!     (cij_kl(2) + cij_kl(21))* sintwophisq) + &
!     8.d0* cij_kl(5)* cosphi*cosphisq* costheta*costhetasq* sintheta - &
!     8.d0* cij_kl(8)* costhetasq* sinphisq* sinthetasq - &
!     8.d0* cij_kl(12)* sinthetafour + &
!     8.d0* cosphisq* costhetasq* sintheta* ((cij_kl(4) + &
!     cij_kl(20))* costheta* sinphi - &
!     (cij_kl(3) + cij_kl(19))*sintheta) + &
!     8.d0* cosphi* costheta* (-cij_kl(11)* costheta*costhetasq* &
!     sinphi*sinphisq + (cij_kl(10) + cij_kl(18))* costhetasq* sinphisq* sintheta + &
!     cij_kl(14)* sintheta*sinthetasq) + 2.d0* sinphi* (cij_kl(13) + &
!     cij_kl(9)* sinphisq)* sintwotheta + &
!     sinphi* (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta))
!
! cij_kl_spherical(2) = ONE_FOURTH * (costhetasq* (cij_kl(1) + 3.d0* cij_kl(2) + cij_kl(7) - &
!      cij_kl(21) + (-cij_kl(1) + cij_kl(2) - cij_kl(7) + &
!      cij_kl(21))* cosfourphi + (-cij_kl(6) + cij_kl(11))* sinfourphi) + &
!      4.d0* (cij_kl(8)* cosphisq - cij_kl(15)* cosphi* sinphi + &
!      cij_kl(3)* sinphisq)* sinthetasq - &
!      2.d0* (cij_kl(10)* cosphisq*cosphi + &
!      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
!      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
!      cij_kl(4)* sinphisq*sinphi)* sintwotheta)
!
! cij_kl_spherical(3) = ONE_EIGHTH * (sintwophi* (3.d0* cij_kl(15) - cij_kl(17) + &
!     4.d0* (cij_kl(2) + cij_kl(21))* costhetasq* sintwophi* sinthetasq) + &
!     4.d0* cij_kl(12)* sintwothetasq + 4.d0* cij_kl(1)* cosphifour* sintwothetasq + &
!     2.d0* cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
!     cij_kl(5)* sinfourtheta) + 2.d0* cosphisq* (3.d0* cij_kl(3) -  cij_kl(19) + &
!     (cij_kl(3) + cij_kl(19))* cosfourtheta + &
!     (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
!     2.d0* sinphi* (sinphi* (3.d0* cij_kl(8) - &
!     cij_kl(16) + (cij_kl(8) + cij_kl(16))* cosfourtheta + &
!     2.d0* cij_kl(7)* sinphisq* sintwothetasq)+ &
!     (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta)+ &
!     2.d0* cosphi* ((cij_kl(15) + cij_kl(17))* cosfourtheta* sinphi + &
!     8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
!     (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)*sinfourtheta))
!
! cij_kl_spherical(4) = ONE_EIGHTH * (cosphi* costheta *(5.d0* cij_kl(4) - &
!     cij_kl(9) + 4.d0* cij_kl(13) - &
!     3.d0* cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
!     4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
!     ONE_HALF* (cij_kl(4) - cij_kl(9) + &
!     cij_kl(20))* costhreephi * (costheta + 3.d0* costhreetheta) - &
!     costheta* (-cij_kl(5) + 5.d0* cij_kl(10) + &
!     4.d0* cij_kl(14) - 3.d0* cij_kl(18) + &
!     (3.d0* cij_kl(5) + cij_kl(10) - &
!     4.d0* cij_kl(14) + cij_kl(18))* costwotheta)* sinphi - &
!     ONE_HALF* (cij_kl(5) - cij_kl(10) - cij_kl(18))* (costheta + &
!     3.d0* costhreetheta)* sinthreephi + &
!     4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* costhetasq* sintheta - &
!     4.d0* (cij_kl(1) + cij_kl(3) - cij_kl(7) - cij_kl(8) + cij_kl(16) - cij_kl(19) + &
!     (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + &
!     cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi* sintheta - &
!     4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
!     cij_kl(21))* costhetasq* sinfourphi* sintheta + &
!     costwophi* ((cij_kl(6) + cij_kl(11) + 6.d0* cij_kl(15) - &
!     2.d0* cij_kl(17))* sintheta + &
!     (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))
!
! cij_kl_spherical(5) = ONE_FOURTH * (2.d0* (cij_kl(4) + &
!     cij_kl(20))* cosphisq* (costwotheta + cosfourtheta)* sinphi + &
!     2.d0* cij_kl(9)* (costwotheta + cosfourtheta)* sinphi*sinphisq + &
!     16.d0* cij_kl(1)* cosphifour* costheta*costhetasq* sintheta + &
!     4.d0* costheta*costhetasq* (-2.d0* cij_kl(8)* sinphisq + &
!     4.d0* cij_kl(7)* sinphifour + &
!     (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta + &
!     4.d0* cij_kl(13)* (1.d0 + 2.d0* costwotheta)* sinphi* sinthetasq + &
!     8.d0* costheta* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta*sinthetasq + &
!     2.d0* cosphi*cosphisq* (cij_kl(5)* (costwotheta + cosfourtheta) + &
!     8.d0* cij_kl(6)* costheta*costhetasq* sinphi* sintheta) + &
!     2.d0* cosphi* (cosfourtheta* (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
!     costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
!     8.d0* cij_kl(11)* costheta*costhetasq* sinphi*sinphisq* sintheta) - &
!     (cij_kl(3) + cij_kl(16) + cij_kl(19) + &
!     (cij_kl(3) - cij_kl(16) + cij_kl(19))* costwophi + &
!     (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)
!
! cij_kl_spherical(6) = ONE_HALF * costheta*costhetasq* ((cij_kl(6) + cij_kl(11))* costwophi + &
!      (cij_kl(6) - cij_kl(11))* cosfourphi + 2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
!      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi) + &
!      ONE_FOURTH* costhetasq* (-(cij_kl(4) + 3* cij_kl(9) + cij_kl(20))* cosphi - &
!      3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
!      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
!      3.d0* (cij_kl(5) - cij_kl(10) - cij_kl(18))* sinthreephi)* sintheta + &
!      costheta* ((cij_kl(15) + cij_kl(17))* costwophi + &
!      (-cij_kl(3) + cij_kl(8) + cij_kl(16) - cij_kl(19))* sintwophi)* sinthetasq + &
!      (-cij_kl(13)* cosphi + cij_kl(14)* sinphi)* sintheta*sinthetasq
!
! cij_kl_spherical(7) = cij_kl(7) * cosphifour - cij_kl(11)* cosphi*cosphisq* sinphi + &
!      (cij_kl(2) + cij_kl(21))* cosphisq* sinphisq - &
!      cij_kl(6)* cosphi* sinphi*sinphisq + &
!      cij_kl(1)* sinphifour
!
! cij_kl_spherical(8) = ONE_HALF * (2.d0* costhetasq* sinphi* (-cij_kl(15)* cosphi + &
!      cij_kl(3)* sinphi) + 2.d0* cij_kl(2)* cosphifour* sinthetasq + &
!      (2.d0* cij_kl(2)* sinphifour + &
!      (cij_kl(1) + cij_kl(7) - cij_kl(21))* sintwophisq)* sinthetasq + &
!      cij_kl(4)* sinphi*sinphisq* sintwotheta + &
!      cosphi*cosphisq* (2.d0* (-cij_kl(6) + cij_kl(11))* sinphi* sinthetasq + &
!      cij_kl(10)* sintwotheta) + cosphi* sinphisq* (2.d0* (cij_kl(6) - &
!      cij_kl(11))* sinphi* sinthetasq + &
!      (cij_kl(5) - cij_kl(18))* sintwotheta) + &
!      cosphisq* (2.d0* cij_kl(8)* costhetasq + &
!      (cij_kl(9) - cij_kl(20))* sinphi* sintwotheta))
!
! cij_kl_spherical(9) = cij_kl(11) * cosphifour* sintheta - sinphi*sinphisq* (cij_kl(5)* costheta + &
!      cij_kl(6)* sinphi* sintheta) +  cosphisq* sinphi* (-(cij_kl(10) + &
!      cij_kl(18))* costheta + &
!      3.d0* (cij_kl(6) - cij_kl(11))* sinphi* sintheta) + &
!      cosphi* sinphisq* ((cij_kl(4) + cij_kl(20))* costheta + &
!      2.d0* (-2.d0* cij_kl(1) + cij_kl(2) + cij_kl(21))* sinphi* sintheta) + &
!      cosphi*cosphisq* (cij_kl(9)* costheta - 2.d0* (cij_kl(2) - 2.d0* cij_kl(7) + &
!      cij_kl(21))* sinphi* sintheta)
!
! cij_kl_spherical(10) = ONE_FOURTH * (4.d0* costwotheta* (cij_kl(10)* cosphi*cosphisq + &
!      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
!      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
!      cij_kl(4)* sinphi*sinphisq) + (cij_kl(1) + 3.d0* cij_kl(2) - &
!      2.d0* cij_kl(3) + cij_kl(7) - &
!      2.d0* cij_kl(8) - cij_kl(21) + 2.d0* (cij_kl(3) - cij_kl(8))* costwophi + &
!      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
!      2.d0* cij_kl(15)* sintwophi + &
!      (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)
!
! cij_kl_spherical(11) = ONE_FOURTH * (2.d0* costheta* ((cij_kl(6) + cij_kl(11))* costwophi + &
!      (-cij_kl(6) + cij_kl(11))* cosfourphi + &
!      2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
!      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(21))* sinfourphi) + &
!      (-(cij_kl(4) + 3.d0* cij_kl(9) + cij_kl(20))* cosphi + &
!      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
!      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
!      (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sintheta)
!
! cij_kl_spherical(12) = ONE_SIXTEENTH * (cij_kl(16) - 2.d0* cij_kl(16)* cosfourtheta* sinphisq + &
!      costwophi* (-cij_kl(16) + 8.d0* costheta* sinthetasq* ((cij_kl(3) - &
!      cij_kl(8) + cij_kl(19))* costheta + &
!      (cij_kl(5) - cij_kl(10) - cij_kl(18))* cosphi* sintheta)) + &
!      2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq + &
!      2.d0* (8.d0* cij_kl(12)* costhetafour + &
!      8.d0* cij_kl(14)* cosphi* costheta*costhetasq* sintheta + &
!      4.d0* cosphi* costheta* (cij_kl(5) + cij_kl(10) + cij_kl(18) + &
!      (cij_kl(4) + cij_kl(20))* sintwophi)* &
!      sintheta*sinthetasq + 8.d0* cij_kl(1)* cosphifour* sinthetafour + &
!      8.d0* cij_kl(6)* cosphi*cosphisq* sinphi* sinthetafour + &
!      8.d0* cij_kl(11)* cosphi* sinphi*sinphisq* sinthetafour + &
!      8.d0* cij_kl(7)* sinphifour* sinthetafour + &
!      2.d0* cij_kl(2)* sintwophisq* sinthetafour + &
!      2.d0* cij_kl(21)* sintwophisq* sinthetafour + &
!      2.d0* cij_kl(13)* sinphi* sintwotheta + &
!      2.d0* cij_kl(9)* sinphi*sinphisq* sintwotheta + &
!      cij_kl(3)* sintwothetasq + cij_kl(8)* sintwothetasq + &
!      cij_kl(19)* sintwothetasq + cij_kl(13)* sinphi* sinfourtheta - &
!      cij_kl(9)* sinphi*sinphisq* sinfourtheta))
!
! cij_kl_spherical(13) = ONE_EIGHTH * (cosphi* costheta* (cij_kl(4) + 3.d0* cij_kl(9) + &
!      4.d0* cij_kl(13) + cij_kl(20) - (cij_kl(4) + 3.d0* cij_kl(9) - &
!      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + 4.d0* (-cij_kl(1) - &
!      cij_kl(3) + cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19) + &
!      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
!      cij_kl(19))* costwotheta)* sintwophi* sintheta + &
!      4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* sinthetasq*sintheta - &
!      4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
!      cij_kl(21))* sinfourphi* sinthetasq*sintheta + &
!      costheta* ((-3.d0* cij_kl(5) - cij_kl(10) - 4.d0* cij_kl(14) - &
!      cij_kl(18) + (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + &
!      cij_kl(18))* costwotheta)* sinphi + 6.d0* ((cij_kl(4) - cij_kl(9) + &
!      cij_kl(20))* costhreephi + (-cij_kl(5) + cij_kl(10) + &
!      cij_kl(18))* sinthreephi)* sinthetasq) + costwophi* ((3* cij_kl(6) + &
!      3.d0* cij_kl(11) + 2.d0* (cij_kl(15) + cij_kl(17)))* sintheta - &
!      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
!      cij_kl(17)))* sinthreetheta))
!
! cij_kl_spherical(14) = ONE_FOURTH * (2.d0* cij_kl(13)* (costwotheta + cosfourtheta)* sinphi + &
!      8.d0* costheta*costhetasq* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta + &
!      4.d0* (cij_kl(4) + cij_kl(20))* cosphisq* (1.d0 + &
!      2.d0* costwotheta)* sinphi* sinthetasq + &
!      4.d0* cij_kl(9)* (1.d0 + 2.d0* costwotheta)* sinphi*sinphisq* sinthetasq + &
!      16.d0* cij_kl(1)* cosphifour* costheta* sintheta*sinthetasq + &
!      4.d0* costheta* (-2.d0* cij_kl(8)* sinphisq + 4.d0* cij_kl(7)* sinphifour + &
!      (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta*sinthetasq + &
!      4.d0* cosphi*cosphisq* sinthetasq* (cij_kl(5) + 2.d0* cij_kl(5)* costwotheta + &
!      4.d0* cij_kl(6)* costheta* sinphi* sintheta) + &
!      2.d0* cosphi* (cosfourtheta* (cij_kl(14) - (cij_kl(10) + cij_kl(18))* sinphisq) + &
!      costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
!      8.d0* cij_kl(11)* costheta* sinphi*sinphisq* sintheta*sinthetasq) + &
!      (cij_kl(3) + cij_kl(16) + cij_kl(19) + (cij_kl(3) - cij_kl(16) + &
!      cij_kl(19))* costwophi + (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)
!
! cij_kl_spherical(15) = costwophi * costheta* (-cij_kl(17) + (cij_kl(15) + cij_kl(17))* costhetasq) + &
!       ONE_SIXTEENTH* (-((11.d0* cij_kl(4) + cij_kl(9) + 4.d0* cij_kl(13) - &
!       5.d0* cij_kl(20))* cosphi + (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
!       (cij_kl(5) + 11.d0* cij_kl(10) + 4.d0* cij_kl(14) - &
!       5.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
!       cij_kl(18))* sinthreephi)* sintheta + &
!       8.d0* costheta* ((-cij_kl(1) - cij_kl(3) + cij_kl(7) + cij_kl(8) - cij_kl(16) +&
!       cij_kl(19) + (cij_kl(1) - cij_kl(3) - &
!       cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi +&
!       ((cij_kl(6) + cij_kl(11))* costwophi + &
!       (cij_kl(6) - cij_kl(11))* cosfourphi + (-cij_kl(1) + cij_kl(2) - cij_kl(7) +&
!       cij_kl(21))* sinfourphi)* sinthetasq) +&
!       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
!       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
!       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
!       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)
!
! cij_kl_spherical(16) = ONE_FOURTH *(cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
!       cij_kl(19) + cij_kl(21) + 2.d0*(cij_kl(16) - cij_kl(19))*costwophi* costhetasq + &
!       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(16) + &
!       cij_kl(19) - cij_kl(21))*costwotheta - 2.d0* cij_kl(17)* costhetasq* sintwophi + &
!       2.d0* ((-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
!       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sinthetasq + ((cij_kl(5) - cij_kl(10) +&
!       cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) + cij_kl(18))* costhreephi +&
!       (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - &
!       (cij_kl(4) - cij_kl(9) + cij_kl(20))* sinthreephi)* sintwotheta)
!
! cij_kl_spherical(17) = ONE_EIGHTH * (4.d0* costwophi* costheta* (cij_kl(6) + cij_kl(11) - &
!       2.d0* cij_kl(15) - (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
!       cij_kl(17)))* costwotheta) - (2.d0* cosphi* (-3.d0* cij_kl(4) +&
!       cij_kl(9) + 2.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) - cij_kl(9) + &
!       cij_kl(20))* costwophi) - (cij_kl(5) - 5.d0* cij_kl(10) + &
!       4.d0* cij_kl(14) + 3.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
!       cij_kl(18))* sinthreephi)* sintheta + &
!       8.d0* costheta* ((-cij_kl(1) + cij_kl(3) + cij_kl(7) - cij_kl(8) + &
!       (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
!       cij_kl(19))* costwotheta)* sintwophi + ((cij_kl(6) - cij_kl(11))* cosfourphi + &
!       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi)* sinthetasq) +&
!       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
!       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
!       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
!       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)
!
! cij_kl_spherical(18) = ONE_HALF * ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi* costwotheta - &
!       (cij_kl(5) - cij_kl(10) - cij_kl(18))* costhreephi* costwotheta - &
!       2.d0* (cij_kl(4) - cij_kl(9) + &
!       (cij_kl(4) - cij_kl(9) + cij_kl(20))* costwophi)* costwotheta* sinphi + &
!       (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + cij_kl(21) + &
!       (-cij_kl(16) + cij_kl(19))* costwophi + &
!       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
!       cij_kl(17)* sintwophi + &
!       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)
!
! cij_kl_spherical(19) = ONE_FOURTH * (cij_kl(16) - cij_kl(16)* costwophi + &
!      (-cij_kl(15) + cij_kl(17))* sintwophi + &
!      4.d0* cij_kl(12)* sintwothetasq + &
!      2.d0* (2.d0* cij_kl(1)* cosphifour* sintwothetasq + &
!      cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
!      cij_kl(5)* sinfourtheta) + cosphisq* (-cij_kl(3) + cij_kl(19) + (cij_kl(3) +&
!      cij_kl(19))* cosfourtheta + (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
!      sinphi* (cosfourtheta* ((cij_kl(15) + cij_kl(17))* cosphi + &
!      cij_kl(16)* sinphi) + (cij_kl(2) + cij_kl(7) - 2.d0* cij_kl(8) + cij_kl(21) + &
!      (cij_kl(2) - cij_kl(7) + cij_kl(21))* costwophi)* sinphi* sintwothetasq + &
!      (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta) + &
!      cosphi* (8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
!      (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)* sinfourtheta)))
!
! cij_kl_spherical(20) = ONE_EIGHTH * (2.d0* cosphi* costheta* (-3.d0* cij_kl(4) - cij_kl(9) + &
!      4.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
!      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
!      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi* (costheta + &
!      3.d0* costhreetheta) - &
!      2.d0* costheta* (-cij_kl(5) - 3.d0* cij_kl(10) + 4.d0* cij_kl(14) + &
!      cij_kl(18) + (3.d0* cij_kl(5) + &
!      cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))*costwotheta)* sinphi - &
!      (cij_kl(5) - cij_kl(10) - cij_kl(18))* &
!      (costheta + 3.d0* costhreetheta)* sinthreephi + 8.d0* (cij_kl(6) - &
!      cij_kl(11))* cosfourphi* costhetasq* sintheta - 8.d0* (cij_kl(1) - &
!      cij_kl(3) - cij_kl(7) + cij_kl(8) + &
!      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
!      cij_kl(19))* costwotheta)* sintwophi* sintheta - &
!      8.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
!      cij_kl(21))* costhetasq* sinfourphi* sintheta + &
!      2.d0* costwophi* ((cij_kl(6) + cij_kl(11) - 2.d0* cij_kl(15) + &
!      2.d0* cij_kl(17))* sintheta + &
!      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))
!
! cij_kl_spherical(21) = ONE_FOURTH * (cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
!      cij_kl(19) + cij_kl(21) - 2.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
!      cij_kl(21))* cosfourphi* costhetasq + &
!      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + &
!      cij_kl(21))* costwotheta + &
!      2.d0* (-cij_kl(6) + cij_kl(11))* costhetasq* sinfourphi - &
!      2.d0* ((-cij_kl(16) + cij_kl(19))* costwophi + cij_kl(17)* sintwophi)* sinthetasq - &
!      ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) +&
!      cij_kl(18))* costhreephi + &
!      (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - (cij_kl(4) - cij_kl(9) + &
!      cij_kl(20))* sinthreephi)* sintwotheta)

  end subroutine rotate_tensor_global_to_radial_vector

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_tensor_global_to_radial(theta,phi, &
                                            d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                            d33,d34,d35,d36,d44,d45,d46,d55,d56,d66, &
                                            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! rotates from global (c_ij) to local (d_ij) anisotropic parameters.
  ! The c_ij are the coefficients in the global reference frame used in SPECFEM3D

  implicit none

  double precision,intent(in) :: theta,phi
  double precision,intent(out) :: d11,d12,d13,d14,d15,d16,d22,d23,d24,d25,d26, &
                                 d33,d34,d35,d36,d44,d45,d46,d55,d56,d66

  double precision,intent(in) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                 c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

! Purpose : compute the kernels in r,theta,phi (cij_kl_spherical)
! from the kernels in x,y,z (cij_kl) (x,y,z to r,theta,phi)
! At r,theta,phi fixed
! theta and phi are in radians

  ! local parameters
  double precision,dimension(3,3) :: rotmat
  double precision,dimension(6,6) :: cij,dij
  integer :: i,j

  ! rotation matrix
  ! rotates pole (Cartesian) to spherical (radial) position
  ! First column
  rotmat(1,1) =  cos(phi) * cos(theta)
  rotmat(2,1) = -sin(phi)
  rotmat(3,1) =  cos(phi) * sin(theta)

  ! Second column
  rotmat(1,2) =  sin(phi) * cos(theta)
  rotmat(2,2) =  cos(phi)
  rotmat(3,2) =  sin(phi) * sin(theta)

  ! Third column
  rotmat(1,3) = -sin(theta)
  rotmat(2,3) =  0.d0
  rotmat(3,3) =  cos(theta)

  ! Cij Voigt notation
  cij(1,1) = c11
  cij(1,2) = c12
  cij(1,3) = c13
  cij(1,4) = c14
  cij(1,5) = c15
  cij(1,6) = c16
  cij(2,2) = c22
  cij(2,3) = c23
  cij(2,4) = c24
  cij(2,5) = c25
  cij(2,6) = c26
  cij(3,3) = c33
  cij(3,4) = c34
  cij(3,5) = c35
  cij(3,6) = c36
  cij(4,4) = c44
  cij(4,5) = c45
  cij(4,6) = c46
  cij(5,5) = c55
  cij(5,6) = c56
  cij(6,6) = c66
  ! fills lower-triangle, for example: C(2,1) = C(1,2) <->  c21 = c12
  !                                    C(3,1) = C(1,3) <->  c31 = c13
  !                                    C(3,2) = C(2,3) <->  c32 = c23
  do j = 2,6
    do i = 1,j - 1
      cij(j,i) = cij(i,j)
    enddo
  enddo

  call rotate_tensor(cij,rotmat,dij)

  ! returns local dij
  d11 = dij(1,1)
  d12 = dij(1,2)
  d13 = dij(1,3)
  d14 = dij(1,4)
  d15 = dij(1,5)
  d16 = dij(1,6)
  d22 = dij(2,2)
  d23 = dij(2,3)
  d24 = dij(2,4)
  d25 = dij(2,5)
  d26 = dij(2,6)
  d33 = dij(3,3)
  d34 = dij(3,4)
  d35 = dij(3,5)
  d36 = dij(3,6)
  d44 = dij(4,4)
  d45 = dij(4,5)
  d46 = dij(4,6)
  d55 = dij(5,5)
  d56 = dij(5,6)
  d66 = dij(6,6)

! original routines
!
! note: these calculations were not consistent with the rotation from local to global above.
!       we will replace them here with the rotation by "simpler" multiplication with a Bond matrix.
!
!       one can double check the resulting coefficients with sympy formulations created by running the script:
!       ./utils/scripts/rotate_elastic_tensor_symbolic.py 2
!
!
! Coeff from Min's routine rotate_anisotropic_tensor
! with the help of Collect[Expand[cij],{dij}] in Mathematica

! Definition of the output array cij_kl_spherical :
! cij_kl_spherical(1) = C11 ; cij_kl_spherical(2) = C12 ; cij_kl_spherical(3) = C13
! cij_kl_spherical(4) = C14 ; cij_kl_spherical(5) = C15 ; cij_kl_spherical(6) = C16
! cij_kl_spherical(7) = C22 ; cij_kl_spherical(8) = C23 ; cij_kl_spherical(9) = C24
! cij_kl_spherical(10) = C25 ; cij_kl_spherical(11) = C26 ; cij_kl_spherical(12) = C33
! cij_kl_spherical(13) = C34 ; cij_kl_spherical(14) = C35 ; cij_kl_spherical(15) = C36
! cij_kl_spherical(16) = C44 ; cij_kl_spherical(17) = C45 ; cij_kl_spherical(18) = C46
! cij_kl_spherical(19) = C55 ; cij_kl_spherical(20) = C56 ; cij_kl_spherical(21) = C66
!
! where the Cij (Voigt's notation) are defined as function of
! the components of the elastic tensor in spherical coordinates
! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)
!
! based on original routine rotate_kernels_dble() in save_kernel.F90
!
!  use constants, only: ONE_SIXTEENTH,ONE_FOURTH,ONE_EIGHTH,ONE_HALF
!  implicit none
!
!  double precision :: costheta,sintheta,cosphi,sinphi
!  double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
!  double precision :: costwotheta,sintwotheta,costwophi,sintwophi
!  double precision :: cosfourtheta,sinfourtheta,cosfourphi,sinfourphi
!  double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
!  double precision :: sintwophisq,sintwothetasq
!  double precision :: costhreetheta,sinthreetheta,costhreephi,sinthreephi
!
!  costheta = cos(theta)
!  sintheta = sin(theta)
!  cosphi = cos(phi)
!  sinphi = sin(phi)
!
!  costhetasq = costheta * costheta
!  sinthetasq = sintheta * sintheta
!  cosphisq = cosphi * cosphi
!  sinphisq = sinphi * sinphi
!
!  costhetafour = costhetasq * costhetasq
!  sinthetafour = sinthetasq * sinthetasq
!  cosphifour = cosphisq * cosphisq
!  sinphifour = sinphisq * sinphisq
!
!  costwotheta = cos(2.d0 * theta)
!  sintwotheta = sin(2.d0 * theta)
!  costwophi = cos(2.d0 * phi)
!  sintwophi = sin(2.d0 * phi)
!
!  costhreetheta = cos(3.d0 * theta)
!  sinthreetheta = sin(3.d0 * theta)
!  costhreephi = cos(3.d0 * phi)
!  sinthreephi = sin(3.d0 * phi)
!
!  cosfourtheta = cos(4.d0 * theta)
!  sinfourtheta = sin(4.d0 * theta)
!  cosfourphi = cos(4.d0 * phi)
!  sinfourphi = sin(4.d0 * phi)
!  sintwothetasq = sintwotheta * sintwotheta
!  sintwophisq = sintwophi * sintwophi
!
!  d11 = ONE_SIXTEENTH * (c44 - c44 * costwophi + &
!     16.d0 * cosphi*cosphisq * costhetafour * (c11 * cosphi + c16 * sinphi) + &
!     2.d0 * (c36 + c45) * sintwophi * sintwothetasq - &
!     2.d0 * (c44 * cosfourtheta * sinphisq + &
!     2.d0 * costhetafour * (-4.d0 * c22 * sinphifour - &
!     (c12 + c66) * sintwophisq) + &
!     8.d0 * c15 * cosphi*cosphisq * costheta*costhetasq * sintheta - &
!     8.d0 * c23 * costhetasq * sinphisq * sinthetasq - &
!     8.d0 * c33 * sinthetafour + &
!     8.d0 * cosphisq * costhetasq * sintheta * ((c14 + &
!     c56) * costheta * sinphi - &
!     (c13 + c55)*sintheta) + &
!     8.d0 * cosphi * costheta * (-c26 * costheta*costhetasq * &
!     sinphi*sinphisq + (c25 + c46) * costhetasq * sinphisq * sintheta + &
!     c35 * sintheta*sinthetasq) + 2.d0 * sinphi * (c34 + &
!     c24 * sinphisq) * sintwotheta + &
!     sinphi * (-c34 + c24 * sinphisq) * sinfourtheta))
!
!  d12 = ONE_FOURTH * (costhetasq * (c11 + 3.d0 * c12 + c22 - &
!      c66 - (-c11 + c12 - c22 + &
!      c66) * cosfourphi + (-c16 + c26) * sinfourphi) + &
!      4.d0 * (c23 * cosphisq - c36 * cosphi* sinphi + &
!      c13 * sinphisq) * sinthetasq - &
!      2.d0 * (c25 * cosphisq*cosphi + &
!      (c24 - c56) * cosphisq * sinphi + &
!      (c15 - c46) * cosphi* sinphisq + &
!      c14 * sinphisq*sinphi) * sintwotheta)
!
!  d13 = ONE_EIGHTH * (sintwophi* (3.d0 * c36 - c45 + &
!     4.d0* (c12 + c66) * costhetasq * sintwophi* sinthetasq) + &
!     4.d0* c33 * sintwothetasq + 4.d0 * c11 * cosphifour* sintwothetasq + &
!     2.d0* cosphi*cosphisq * (8.d0 * c16 * costhetasq * sinphi* sinthetasq + &
!     c15* sinfourtheta) + 2.d0* cosphisq * (3.d0* c13 -  c55 + &
!     (c13 + c55) * cosfourtheta + &
!     (c14 + c56) * sinphi* sinfourtheta) + &
!     2.d0* sinphi* (sinphi* (3.d0* c23 - &
!     c44 + (c23 + c44) * cosfourtheta + &
!     2.d0* c22* sinphisq * sintwothetasq)+ &
!     (-c34 + c24* sinphisq) * sinfourtheta)+ &
!     2.d0* cosphi* ((c36 + c45) * cosfourtheta* sinphi + &
!     8.d0* c26* costhetasq * sinphi*sinphisq * sinthetasq + &
!     (-c35 + (c25 + c46) * sinphisq)*sinfourtheta))
!
!  d14 = ONE_EIGHTH * (cosphi* costheta *(5.d0* c14 - &
!     c24 + 4.d0* c34 - &
!     3.d0* c56 + (c14 + 3.d0* c24 - &
!     4.d0* c34 + c56) * costwotheta) + &
!     ONE_HALF* (c14 - c24 + &
!     c56) * costhreephi * (costheta + 3.d0* costhreetheta) - &
!     costheta* (-c15 + 5.d0* c25 + &
!     4.d0* c35 - 3.d0* c46 + &
!     (3.d0* c15 + c25 - &
!     4.d0* c35 + c46) * costwotheta) * sinphi - &
!     ONE_HALF* (c15 - c25 - c46) * (costheta + &
!     3.d0* costhreetheta) * sinthreephi + &
!     4.d0* (c16 - c26) * cosfourphi* costhetasq * sintheta - &
!     4.d0* (c11 + c13 - c22 - c23 + c44 - c55 + &
!     (c11 - c13 - c22 + c23 + &
!     c44 - c55) * costwotheta) * sintwophi* sintheta - &
!     4.d0* (c11 - c12 + c22 - &
!     c66) * costhetasq * sinfourphi* sintheta + &
!     costwophi* ((c16 + c26 + 6.d0* c36 - &
!     2.d0* c45) * sintheta + &
!     (c16 + c26 - 2.d0* (c36 + c45)) * sinthreetheta))
!
!  d15 = ONE_FOURTH * (2.d0* (c14 + &
!     c56) * cosphisq * (costwotheta + cosfourtheta) * sinphi + &
!     2.d0* c24* (costwotheta + cosfourtheta) * sinphi*sinphisq + &
!     16.d0* c11* cosphifour* costheta*costhetasq * sintheta + &
!     4.d0* costheta*costhetasq * (-2.d0* c23* sinphisq + &
!     4.d0* c22* sinphifour + &
!     (c12 + c66) * sintwophisq) * sintheta + &
!     4.d0* c34* (1.d0 + 2.d0* costwotheta) * sinphi* sinthetasq + &
!     8.d0* costheta* (-2.d0* c33 + c23* sinphisq) * sintheta*sinthetasq + &
!     2.d0* cosphi*cosphisq * (c15* (costwotheta + cosfourtheta) + &
!     8.d0* c16* costheta*costhetasq * sinphi* sintheta) + &
!     2.d0* cosphi* (cosfourtheta* (-c35 + (c25 + c46) * sinphisq) + &
!     costwotheta* (c35 + (c25 + c46) * sinphisq) + &
!     8.d0* c26* costheta*costhetasq * sinphi*sinphisq * sintheta) - &
!     (c13 + c44 + c55 + &
!     (c13 - c44 + c55) * costwophi + &
!     (c36 + c45) * sintwophi) * sinfourtheta)
!
!  d16 = ONE_HALF * costheta*costhetasq * ((c16 + c26) * costwophi + &
!      (c16 - c26) * cosfourphi + 2.d0* (-c11 + c22) * sintwophi + &
!      (-c11 + c12 - c22 + c66)) * sinfourphi) + &
!      ONE_FOURTH* costhetasq * (-(c14 + 3.d0* c24 + c56) * cosphi - &
!      3.d0* (c14 - c24 + c56) * costhreephi + &
!      (3.d0* c15 + c25 + c46) * sinphi + &
!      3.d0* (c15 - c25 - c46) * sinthreephi) * sintheta + &
!      costheta* ((c36 + c45) * costwophi + &
!      (-c13 + c23 + c44 - c55) * sintwophi) * sinthetasq + &
!      (-c34* cosphi + c35* sinphi) * sintheta*sinthetasq
!
!  d22 = c22 * cosphifour - c26* cosphi*cosphisq * sinphi + &
!      (c12 + c66) * cosphisq * sinphisq - &
!      c16* cosphi* sinphi*sinphisq + &
!      c11* sinphifour
!
!  d23 = ONE_HALF * (2.d0* costhetasq * sinphi* (-c36* cosphi + &
!      c13* sinphi) + 2.d0* c12* cosphifour* sinthetasq + &
!      (2.d0* c12* sinphifour + &
!      (c11 + c22 - c66) * sintwophisq) * sinthetasq + &
!      c14* sinphi*sinphisq * sintwotheta + &
!      cosphi*cosphisq * (2.d0* (-c16 + c26) * sinphi* sinthetasq + &
!      c25* sintwotheta) + cosphi* sinphisq * (2.d0* (c16 - &
!      c26) * sinphi* sinthetasq + &
!      (c15 - c46) * sintwotheta) + &
!      cosphisq * (2.d0* c23* costhetasq + &
!      (c24 - c56) * sinphi* sintwotheta))
!
!  d24 = c26 * cosphifour* sintheta - sinphi*sinphisq * (c15* costheta + &
!      c16* sinphi* sintheta) +  cosphisq * sinphi* (-(c25 + &
!      c46) * costheta + &
!      3.d0* (c16 - c26) * sinphi* sintheta) + &
!      cosphi* sinphisq * ((c14 + c56) * costheta + &
!      2.d0* (-2.d0* c11 + c12 + c66) * sinphi* sintheta) + &
!      cosphi*cosphisq * (c24* costheta - 2.d0* (c12 - 2.d0* c22 + &
!      c66) * sinphi* sintheta)
!
!  d25 = ONE_FOURTH * (4.d0* costwotheta* (c25* cosphi*cosphisq + &
!      (c24 - c56) * cosphisq * sinphi + &
!      (c15 - c46) * cosphi* sinphisq + &
!      c14* sinphi*sinphisq) + (c11 + 3.d0* c12 - &
!      2.d0* c13 + c22 - &
!      2.d0* c23 - c66 + 2.d0* (c13 - c23) * costwophi + &
!      (-c11 + c12 - c22 + c66) * cosfourphi + &
!      2.d0* c36* sintwophi + &
!      (-c16 + c26) * sinfourphi) * sintwotheta)
!
!  d26 = ONE_FOURTH * (2.d0* costheta* ((c16 + c26) * costwophi + &
!      (-c16 + c26) * cosfourphi + &
!      2.d0* (-c11 + c22) * sintwophi + &
!      (c11 - c12 + c22 - c66) * sinfourphi) + &
!      (-(c14 + 3.d0* c24 + c56) * cosphi + &
!      (c14 - c24 + c56) * costhreephi + &
!      (3.d0* c15 + c25 + c46) * sinphi + &
!      (-c15 + c25 + c46) * sinthreephi) * sintheta)
!
!  d33 = ONE_SIXTEENTH * (c44 - 2.d0* c44* cosfourtheta* sinphisq + &
!      costwophi* (-c44 + 8.d0* costheta* sinthetasq * ((c13 - &
!      c23 + c55) * costheta + &
!      (c15 - c25 - c46) * cosphi* sintheta)) + &
!      2.d0* (c36 + c45) * sintwophi* sintwothetasq + &
!      2.d0* (8.d0* c33* costhetafour + &
!      8.d0* c35* cosphi* costheta*costhetasq * sintheta + &
!      4.d0* cosphi* costheta* (c15 + c25 + c46 + &
!      (c14 + c56) * sintwophi) * &
!      sintheta*sinthetasq + 8.d0* c11* cosphifour* sinthetafour + &
!      8.d0* c16* cosphi*cosphisq * sinphi* sinthetafour + &
!      8.d0* c26* cosphi* sinphi*sinphisq * sinthetafour + &
!      8.d0* c22* sinphifour* sinthetafour + &
!      2.d0* c12* sintwophisq * sinthetafour + &
!      2.d0* c66* sintwophisq * sinthetafour + &
!      2.d0* c34* sinphi* sintwotheta + &
!      2.d0* c24* sinphi*sinphisq * sintwotheta + &
!      c13* sintwothetasq + c23* sintwothetasq + &
!      c55* sintwothetasq + c34* sinphi* sinfourtheta - &
!      c24* sinphi*sinphisq * sinfourtheta))
!
!  d34 = ONE_EIGHTH * (cosphi* costheta* (c14 + 3.d0* c24 + &
!      4.d0* c34 + c56 - (c14 + 3.d0* c24 - &
!      4.d0* c34 + c56) * costwotheta) + 4.d0* (-c11 - &
!      c13 + c22 + c23 + c44 - c55 + &
!      (c11 - c13 - c22 + c23 + c44 - &
!      c55) * costwotheta) * sintwophi* sintheta + &
!      4.d0* (c16 - c26) * cosfourphi* sinthetasq*sintheta - &
!      4.d0* (c11 - c12 + c22 - &
!      c66) * sinfourphi* sinthetasq*sintheta + &
!      costheta* ((-3.d0* c15 - c25 - 4.d0* c35 - &
!      c46 + (3.d0* c15 + c25 - 4.d0* c35 + &
!      c46) * costwotheta) * sinphi + 6.d0* ((c14 - c24 + &
!      c56) * costhreephi + (-c15 + c25 + &
!      c46) * sinthreephi) * sinthetasq) + costwophi* ((3.d0 * c16 + &
!      3.d0* c26 + 2.d0* (c36 + c45)) * sintheta - &
!      (c16 + c26 - 2.d0* (c36 + &
!      c45)) * sinthreetheta))
!
!  d35 = ONE_FOURTH * (2.d0* c34* (costwotheta + cosfourtheta) * sinphi + &
!      8.d0* costheta*costhetasq * (-2.d0* c33 + c23* sinphisq) * sintheta + &
!      4.d0* (c14 + c56) * cosphisq * (1.d0 + &
!      2.d0* costwotheta) * sinphi* sinthetasq + &
!      4.d0* c24* (1.d0 + 2.d0* costwotheta) * sinphi*sinphisq * sinthetasq + &
!      16.d0* c11* cosphifour* costheta* sintheta*sinthetasq + &
!      4.d0* costheta* (-2.d0* c23* sinphisq + 4.d0* c22* sinphifour + &
!      (c12 + c66) * sintwophisq) * sintheta*sinthetasq + &
!      4.d0* cosphi*cosphisq * sinthetasq * (c15 + 2.d0* c15* costwotheta + &
!      4.d0* c16* costheta* sinphi* sintheta) + &
!      2.d0* cosphi* (cosfourtheta* (c35 - (c25 + c46) * sinphisq) + &
!      costwotheta* (c35 + (c25 + c46) * sinphisq) + &
!      8.d0* c26* costheta* sinphi*sinphisq * sintheta*sinthetasq) + &
!      (c13 + c44 + c55 + (c13 - c44 + &
!      c55) * costwophi + (c36 + c45) * sintwophi) * sinfourtheta)
!
!  d36 = costwophi * costheta* (-c45 + (c36 + c45) * costhetasq) + &
!       ONE_SIXTEENTH* (-((11.d0* c14 + c24 + 4.d0* c34 - &
!       5.d0* c56) * cosphi + (c14 - c24 + c56) * costhreephi - &
!       (c15 + 11.d0* c25 + 4.d0* c35 - &
!       5.d0* c46) * sinphi + (-c15 + c25 + &
!       c46) * sinthreephi) * sintheta + &
!       8.d0* costheta* ((-c11 - c13 + c22 + c23 - c44 +&
!       c55 + (c11 - c13 - &
!       c22 + c23 + c44 - c55) * costwotheta) * sintwophi +&
!       ((c16 + c26) * costwophi + &
!       (c16 - c26) * cosfourphi + (-c11 + c12 - c22 +&
!       c66) * sinfourphi) * sinthetasq) +&
!       ((c14 + 3.d0* c24 - 4.d0* c34 + c56) * cosphi + &
!       3.d0* (c14 - c24 + c56) * costhreephi - &
!       (3.d0* c15 + c25 - 4.d0* c35 + c46) * sinphi + &
!       3.d0* (-c15 + c25 + c46) * sinthreephi) * sinthreetheta)
!
!  d44 = ONE_FOURTH *(c11 - c12 + c22 + c44 + &
!       c55 + c66 + 2.d0*(c44 - c55)*costwophi* costhetasq + &
!       (-c11 + c12 - c22 + c44 + &
!       c55 - c66)*costwotheta - 2.d0* c45* costhetasq * sintwophi + &
!       2.d0* ((-c11 + c12 - c22 + c66) * cosfourphi + &
!       (-c16 + c26) * sinfourphi) * sinthetasq + ((c15 - c25 +&
!       c46) * cosphi + (-c15 + c25 + c46) * costhreephi +&
!       (-c14 + c24 + c56) * sinphi - &
!       (c14 - c24 + c56) * sinthreephi) * sintwotheta)
!
!  d45 = ONE_EIGHTH * (4.d0* costwophi* costheta* (c16 + c26 - &
!       2.d0* c36 - (c16 + c26 - 2.d0* (c36 + &
!       c45)) * costwotheta) - (2.d0* cosphi* (-3.d0* c14 +&
!       c24 + 2.d0* c34 + c56 + (c14 - c24 + &
!       c56) * costwophi) - (c15 - 5.d0* c25 + &
!       4.d0* c35 + 3.d0* c46) * sinphi + (-c15 + c25 + &
!       c46) * sinthreephi) * sintheta + &
!       8.d0* costheta* ((-c11 + c13 + c22 - c23 + &
!       (c11 - c13 - c22 + c23 + c44 - &
!       c55) * costwotheta) * sintwophi + ((c16 - c26) * cosfourphi + &
!       (-c11 + c12 - c22 + c66) * sinfourphi) * sinthetasq) +&
!       ((c14 + 3.d0* c24 - 4.d0* c34 + c56) * cosphi + &
!       3.d0* (c14 - c24 + c56) * costhreephi - &
!       (3.d0* c15 + c25 - 4.d0* c35 + c46) * sinphi + &
!       3.d0* (-c15 + c25 + c46) * sinthreephi) * sinthreetheta)
!
!  d46 = ONE_HALF * ((c15 - c25 + c46) * cosphi* costwotheta - &
!       (c15 - c25 - c46) * costhreephi* costwotheta - &
!       2.d0* (c14 - c24 + &
!       (c14 - c24 + c56) * costwophi) * costwotheta* sinphi + &
!       (c11 - c12 + c22 - c44 - c55 + c66 + &
!       (-c44 + c55) * costwophi + &
!       (-c11 + c12 - c22 + c66) * cosfourphi + &
!       c45* sintwophi + &
!       (-c16 + c26) * sinfourphi) * sintwotheta)
!
!  d55 = ONE_FOURTH * (c44 - c44* costwophi + &
!      (-c36 + c45) * sintwophi + &
!      4.d0* c33* sintwothetasq + &
!      2.d0* (2.d0* c11* cosphifour* sintwothetasq + &
!      cosphi*cosphisq * (8.d0* c16* costhetasq * sinphi* sinthetasq + &
!      c15* sinfourtheta) + cosphisq * (-c13 + c55 + (c13 +&
!      c55) * cosfourtheta + (c14 + c56) * sinphi* sinfourtheta) + &
!      sinphi* (cosfourtheta* ((c36 + c45) * cosphi + &
!      c44* sinphi) + (c12 + c22 - 2.d0* c23 + c66 + &
!      (c12 - c22 + c66) * costwophi) * sinphi* sintwothetasq + &
!      (-c34 + c24* sinphisq) * sinfourtheta) + &
!      cosphi* (8.d0* c26* costhetasq * sinphi*sinphisq * sinthetasq + &
!      (-c35 + (c25 + c46) * sinphisq) * sinfourtheta)))
!
!  d56 = ONE_EIGHTH * (2.d0* cosphi* costheta* (-3.d0* c14 - c24 + &
!      4.d0* c34 + c56 + (c14 + 3.d0* c24 - &
!      4.d0* c34 + c56) * costwotheta) + &
!      (c14 - c24 + c56) * costhreephi* (costheta + &
!      3.d0* costhreetheta) - &
!      2.d0* costheta* (-c15 - 3.d0* c25 + 4.d0* c35 + &
!      c46 + (3.d0* c15 + &
!      c25 - 4.d0* c35 + c46)*costwotheta) * sinphi - &
!      (c15 - c25 - c46) * &
!      (costheta + 3.d0* costhreetheta) * sinthreephi + 8.d0* (c16 - &
!      c26) * cosfourphi* costhetasq * sintheta - 8.d0* (c11 - &
!      c13 - c22 + c23 + &
!      (c11 - c13 - c22 + c23 + c44 - &
!      c55) * costwotheta) * sintwophi* sintheta - &
!      8.d0* (c11 - c12 + c22 - &
!      c66) * costhetasq * sinfourphi* sintheta + &
!      2.d0* costwophi* ((c16 + c26 - 2.d0* c36 + &
!      2.d0* c45) * sintheta + &
!      (c16 + c26 - 2.d0* (c36 + c45)) * sinthreetheta))
!
!  d66 = ONE_FOURTH * (c11 - c12 + c22 + c44 + &
!      c55 + c66 + 2.d0* (c11 - c12 + c22 - &
!      c66) * cosfourphi* costhetasq + &
!      (c11 - c12 + c22 - c44 - c55 + &
!      c66) * costwotheta + &
!      2.d0* (-c16 + c26) * costhetasq * sinfourphi - &
!      2.d0* ((-c44 + c55) * costwophi + c45* sintwophi) * sinthetasq - &
!      ((c15 - c25 + c46) * cosphi + (-c15 + c25 + &
!      c46) * costhreephi + &
!      (-c14 + c24 + c56) * sinphi - (c14 - c24 + &
!      c56) * sinthreephi) * sintwotheta)

  end subroutine rotate_tensor_global_to_radial


!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_tensor(cij_in,rotmat,cij_out)

  ! rotates from (6,6)-tensor cij_in to cij_out using the (3,3)-rotation matrix rotmat
  !
  ! cij (Voigt) tensors only need upper triangle set

  implicit none

  double precision,dimension(3,3),intent(in) :: rotmat
  double precision,dimension(6,6),intent(in) :: cij_in
  double precision,dimension(6,6),intent(out) :: cij_out

  ! local parameters
  double precision,dimension(6,6) :: bond,bond_t
  double precision,dimension(6,6) :: tensor_tmp
  integer :: i,j,k

  ! creates Bond matrix (see e.g. Auld, 1973)
  ! First column
  bond(1,1) = rotmat(1,1)*rotmat(1,1)
  bond(2,1) = rotmat(2,1)*rotmat(2,1)
  bond(3,1) = rotmat(3,1)*rotmat(3,1)
  bond(4,1) = rotmat(2,1)*rotmat(3,1)
  bond(5,1) = rotmat(3,1)*rotmat(1,1)
  bond(6,1) = rotmat(1,1)*rotmat(2,1)

  ! Second column
  bond(1,2) = rotmat(1,2)*rotmat(1,2)
  bond(2,2) = rotmat(2,2)*rotmat(2,2)
  bond(3,2) = rotmat(3,2)*rotmat(3,2)
  bond(4,2) = rotmat(2,2)*rotmat(3,2)
  bond(5,2) = rotmat(3,2)*rotmat(1,2)
  bond(6,2) = rotmat(1,2)*rotmat(2,2)

  ! Third column
  bond(1,3) = rotmat(1,3)*rotmat(1,3)
  bond(2,3) = rotmat(2,3)*rotmat(2,3)
  bond(3,3) = rotmat(3,3)*rotmat(3,3)
  bond(4,3) = rotmat(2,3)*rotmat(3,3)
  bond(5,3) = rotmat(3,3)*rotmat(1,3)
  bond(6,3) = rotmat(1,3)*rotmat(2,3)

  ! Fourth column
  bond(1,4) = 2.d0*rotmat(1,2)*rotmat(1,3)
  bond(2,4) = 2.d0*rotmat(2,2)*rotmat(2,3)
  bond(3,4) = 2.d0*rotmat(3,2)*rotmat(3,3)
  bond(4,4) = rotmat(2,2)*rotmat(3,3) + rotmat(2,3)*rotmat(3,2)
  bond(5,4) = rotmat(1,2)*rotmat(3,3) + rotmat(1,3)*rotmat(3,2)
  bond(6,4) = rotmat(1,2)*rotmat(2,3) + rotmat(1,3)*rotmat(2,2)

  ! Fifth column
  bond(1,5) = 2.d0*rotmat(1,3)*rotmat(1,1)
  bond(2,5) = 2.d0*rotmat(2,3)*rotmat(2,1)
  bond(3,5) = 2.d0*rotmat(3,3)*rotmat(3,1)
  bond(4,5) = rotmat(2,1)*rotmat(3,3) + rotmat(2,3)*rotmat(3,1)
  bond(5,5) = rotmat(1,3)*rotmat(3,1) + rotmat(1,1)*rotmat(3,3)
  bond(6,5) = rotmat(1,3)*rotmat(2,1) + rotmat(1,1)*rotmat(2,3)

  ! Sixth column
  bond(1,6) = 2.d0*rotmat(1,1)*rotmat(1,2)
  bond(2,6) = 2.d0*rotmat(2,1)*rotmat(2,2)
  bond(3,6) = 2.d0*rotmat(3,1)*rotmat(3,2)
  bond(4,6) = rotmat(2,2)*rotmat(3,1) + rotmat(2,1)*rotmat(3,2)
  bond(5,6) = rotmat(1,1)*rotmat(3,2) + rotmat(1,2)*rotmat(3,1)
  bond(6,6) = rotmat(1,1)*rotmat(2,2) + rotmat(1,2)*rotmat(2,1)

  bond_t = transpose(bond)

  ! rotates Cij
  ! First compute C M^t
  tensor_tmp(:,:) = 0.d0
  do j = 1,6
    do k = 1,6
      do i = 1,6
        tensor_tmp(i,j) = tensor_tmp(i,j) + cij_in(i,k) * bond_t(k,j)
      enddo
    enddo
  enddo
  ! Second compute M * (C M^t)
  cij_out(:,:) = 0.d0
  do j = 1,6
    do k = 1,6
      do i = 1,j ! half only
        cij_out(i,j) = cij_out(i,j) + bond(i,k) * tensor_tmp(k,j)
      enddo
    enddo
  enddo

  end subroutine rotate_tensor
