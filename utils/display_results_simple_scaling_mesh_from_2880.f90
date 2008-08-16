!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
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

  program display_results_simple_scaling

  implicit none

  integer :: NEX_MAX, NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
          NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
          NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB

! good values for NEX_XI = 2880 determined by trial and error using OpenDX on one slice
 NER_CRUST =            5 !!!!!!!!!!!!!!!!!! 6
 NER_80_MOHO =            8 !!!!!!!!!!!!! 9
 NER_220_80 =           18
 NER_400_220 =           24 !!!!!!!!!! 26
 NER_600_400 =           26 !!!!!!!!!! 28
 NER_670_600 =           9 !!!!!!!!!!! 10
 NER_771_670 =           13
 NER_TOPDDOUBLEPRIME_771 =          216
 NER_CMB_TOPDDOUBLEPRIME =           18
 NER_OUTER_CORE =          268
 NER_TOP_CENTRAL_CUBE_ICB =           33

  print *,'enter NEX_MAX to use:'
  read(*,*) NEX_MAX

  if(mod(NEX_MAX,32) /= 0) stop 'not a multiple of 32, exiting'

    print *
    print *,'! scaling from 2880 for NEX_MAX = ',NEX_MAX
    print *
    print *,'NER_CRUST = ',nint(NER_CRUST * dble(NEX_MAX) / 2880.d0)
    print *,'NER_80_MOHO = ',nint(NER_80_MOHO * dble(NEX_MAX) / 2880.d0)
    print *,'NER_220_80 = ',nint(NER_220_80 * dble(NEX_MAX) / 2880.d0)
    print *,'NER_400_220 = ',nint(NER_400_220 * dble(NEX_MAX) / 2880.d0)
    print *,'NER_600_400 = ',nint(NER_600_400 * dble(NEX_MAX) / 2880.d0)
    print *,'NER_670_600 = ',nint(NER_670_600 * dble(NEX_MAX) / 2880.d0)
    print *,'NER_771_670 = ',nint(NER_771_670 * dble(NEX_MAX) / 2880.d0)
    print *,'NER_TOPDDOUBLEPRIME_771 = ',nint(NER_TOPDDOUBLEPRIME_771 * dble(NEX_MAX) / 2880.d0)
    print *,'NER_CMB_TOPDDOUBLEPRIME = ',nint(NER_CMB_TOPDDOUBLEPRIME * dble(NEX_MAX) / 2880.d0)
    print *,'NER_OUTER_CORE = ',nint(NER_OUTER_CORE * dble(NEX_MAX) / 2880.d0)
    print *,'NER_TOP_CENTRAL_CUBE_ICB = ',nint(NER_TOP_CENTRAL_CUBE_ICB * dble(NEX_MAX) / 2880.d0)
    print *

  end program display_results_simple_scaling

