!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_codes_buffers(idoubling,iz,icode)

! assign codes to buffers for general MPI assembling routine

  implicit none

  include "constants.h"

  integer idoubling
  integer iz,icode

! by default the point is not a matching point
  icode = IFLAG_NORMAL_POINT

! full coupling layer (first layer in contact with a fluid/solid interface)
  if(idoubling == IFLAG_BOTTOM_MANTLE .or. &
     idoubling == IFLAG_BOTTOM_OUTER_CORE .or. &
     idoubling == IFLAG_TOP_OUTER_CORE .or. &
     idoubling == IFLAG_TOP_INNER_CORE) &
        icode = IFLAG_MATCHING_POINT

! one-point coupling layer (second layer in contact with a fluid/solid
! interface located at the bottom in iz = 1)
  if((idoubling == IFLAG_BOTTOM_MANTLE_LEV2 .or. &
      idoubling == IFLAG_BOTTOM_OUTER_CORE_LEV2) .and. iz == 1) &
        icode = IFLAG_MATCHING_POINT

! one-point coupling layer (second layer in contact with a fluid/solid
! interface located at the top in iz = NGLLZ)
  if((idoubling == IFLAG_TOP_OUTER_CORE_LEV2 .or. &
      idoubling == IFLAG_TOP_INNER_CORE_LEV2) .and. iz == NGLLZ) &
        icode = IFLAG_MATCHING_POINT

  end subroutine get_codes_buffers

