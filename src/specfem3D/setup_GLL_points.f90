!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  subroutine setup_GLL_points()

  use specfem_par
  implicit none

  ! local parameters
  integer :: i,j

  ! set up GLL points, weights and derivation matrices
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                                 hprime_xx,hprime_yy,hprime_zz, &
                                 hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                                 wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube)

  if( USE_DEVILLE_PRODUCTS_VAL ) then

  ! check that optimized routines from Deville et al. (2002) can be used
    if(NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) &
      stop 'Deville et al. (2002) routines can only be used if NGLLX = NGLLY = NGLLZ = 5'

    ! define transpose of derivation matrix
    do j = 1,NGLLY
      do i = 1,NGLLX
        hprime_xxT(j,i) = hprime_xx(i,j)
        hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
      enddo
    enddo
  endif

  end subroutine setup_GLL_points
