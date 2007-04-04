!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 7
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!          (c) California Institute of Technology January 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program compute_thickness_crust

  implicit none

! number of different cases to study
  integer, parameter :: NUM_CASES = 9
  integer, dimension(NUM_CASES) :: NEX_val
  integer icase

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: R_EARTH_KM = 6371.d0

! case of multiples of 32
  NEX_val(1) = 160
  NEX_val(2) = 256
  NEX_val(3) = 320
  NEX_val(4) = 480
  NEX_val(5) = 512
  NEX_val(6) = 640
  NEX_val(7) = 864
  NEX_val(8) = 1152
  NEX_val(9) = 1248

! loop on all the cases to study
  do icase = 1,NUM_CASES
    print *,'NEX = ',NEX_val(icase)
    print *,'horizontal size = thickness element = ',2*PI*R_EARTH_KM/(4*NEX_val(icase)),' km'
    print *
  enddo

  end program compute_thickness_crust

