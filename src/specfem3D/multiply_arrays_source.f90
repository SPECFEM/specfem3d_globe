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

! we put these multiplications in a separate routine because otherwise
! some compilers try to unroll the six loops above and take forever to compile
  subroutine multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)

  implicit none

  include "constants.h"

! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,l,m

  integer ir,it,iv

  sourcearrayd(:,k,l,m) = ZERO

  do iv=1,NGLLZ
    do it=1,NGLLY
      do ir=1,NGLLX

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

!================================================================

  subroutine multiply_arrays_adjoint(sourcearrayd,hxir,hetar,hgammar,adj_src_ud)

  implicit none

  include "constants.h"

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
