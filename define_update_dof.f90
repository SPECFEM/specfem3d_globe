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

  subroutine define_update_dof(update_dof,ibool,idoubling,nspec,nglob, &
    IFLAG_BOTTOM,IFLAG_BOTTOM_LEV2,IFLAG_TOP,IFLAG_TOP_LEV2)

!
! create flag for DOFs that need to be updated in the doubling regions
!

  implicit none

  include "constants.h"

  integer nspec,nglob
  integer IFLAG_BOTTOM,IFLAG_BOTTOM_LEV2,IFLAG_TOP,IFLAG_TOP_LEV2
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)
  integer idoubling(nspec)
  logical update_dof(nglob)

  integer i,j,k,ispec

! clear mask
  update_dof(:) = .false.

! define mask
  do ispec = 1,nspec

! update all DOFs for layers in contact with other regions (at CMB and ICB)
  if(idoubling(ispec) == IFLAG_TOP .or. idoubling(ispec) == IFLAG_BOTTOM) then
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          update_dof(ibool(i,j,k,ispec)) = .true.
        enddo
      enddo
    enddo
  endif

! update top layer of DOFs for second matching layer at top of current region
  if(idoubling(ispec) == IFLAG_TOP_LEV2) then
    do j=1,NGLLY
      do i=1,NGLLX
        update_dof(ibool(i,j,NGLLZ,ispec)) = .true.
      enddo
    enddo
  endif

! update bottom layer of DOFs for second matching layer at bottom of current region
  if(idoubling(ispec) == IFLAG_BOTTOM_LEV2) then
    do j=1,NGLLY
      do i=1,NGLLX
        update_dof(ibool(i,j,1,ispec)) = .true.
      enddo
    enddo
  endif

  enddo

  end subroutine define_update_dof

