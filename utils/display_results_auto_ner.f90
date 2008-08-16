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

  program display_results_auto_ner

  implicit none

  include "constants.h"

  integer, parameter :: NCASES = 11

  double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES = 90.d0

  logical, parameter :: CASE_3D = .false.

  integer :: NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
          NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
          NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB

  double precision :: R_CENTRAL_CUBE,DT

  integer, dimension(NCASES) :: NEX_MAX = (/ 480, 512, 640, 672, 864, 1056, 1152, 1248, 2368, 2880, 3264 /)

  integer :: icase

  do icase = 1,NCASES

     call auto_ner(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX(icase), &
          NER_CRUST, NER_80_MOHO, NER_220_80, NER_400_220, NER_600_400, &
          NER_670_600, NER_771_670, NER_TOPDDOUBLEPRIME_771, &
          NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB, &
          R_CENTRAL_CUBE, CASE_3D)

     call auto_time_stepping(ANGULAR_WIDTH_XI_IN_DEGREES, NEX_MAX(icase), DT)
!! DK DK I realize that auto_NER DT values are too small by a factor 1.5, therefore fix it
     DT = DT * 1.5d0

goto 777
    print *
    print *,'Case NEX_XI = ',NEX_MAX(icase)
    print *,'NER_CRUST = ',NER_CRUST
    print *,'NER_80_MOHO = ',NER_80_MOHO
    print *,'NER_220_80 = ',NER_220_80
    print *,'NER_400_220 = ',NER_400_220
    print *,'NER_600_400 = ',NER_600_400
    print *,'NER_670_600 = ',NER_670_600
    print *,'NER_771_670 = ',NER_771_670
    print *,'NER_TOPDDOUBLEPRIME_771 = ',NER_TOPDDOUBLEPRIME_771
    print *,'NER_CMB_TOPDDOUBLEPRIME = ',NER_CMB_TOPDDOUBLEPRIME
    print *,'NER_OUTER_CORE = ',NER_OUTER_CORE
    print *,'NER_TOP_CENTRAL_CUBE_ICB = ',NER_TOP_CENTRAL_CUBE_ICB
    print *,'R_CENTRAL_CUBE = ',R_CENTRAL_CUBE
    print *,'DT = ',DT
    print *
777 continue

!! print the evolution of the time step for Gnuplot
    print *,NEX_MAX(icase),DT

  enddo

  end program display_results_auto_ner

!-----------

  include "auto_ner.f90"

