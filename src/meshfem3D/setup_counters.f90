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


  subroutine setup_counters(NSPEC1D_RADIAL,NSPEC2D_XI,NSPEC2D_ETA,NGLOB1D_RADIAL, &
                        DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                        NPROCTOT,iproc_xi_slice,iproc_eta_slice, &
                        NSPEC1D_RADIAL_CORNER,NSPEC2D_XI_FACE, &
                        NSPEC2D_ETA_FACE,NGLOB1D_RADIAL_CORNER)

! returns: NSPEC1D_RADIAL_CORNER,NSPEC2D_XI_FACE,
!              NSPEC2D_ETA_FACE,NGLOB1D_RADIAL_CORNER

  use constants

  implicit none

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_XI,NSPEC2D_ETA, &
                                         NSPEC1D_RADIAL,NGLOB1D_RADIAL

  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  ! addressing for all the slices
  integer :: NPROCTOT
  integer, dimension(0:NPROCTOT-1) :: iproc_xi_slice,iproc_eta_slice

  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA

! this for the different corners of the slice (which are different if the superbrick is cut)
! 1 : xi_min, eta_min
! 2 : xi_max, eta_min
! 3 : xi_max, eta_max
! 4 : xi_min, eta_max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: &
    NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER
! 1 -> min, 2 -> max
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE


  ! local parameters
  integer :: iregion

  do iregion = 1,MAX_NUM_REGIONS
    NSPEC1D_RADIAL_CORNER(iregion,:) = NSPEC1D_RADIAL(iregion)
    NSPEC2D_XI_FACE(iregion,:) = NSPEC2D_XI(iregion)
    NSPEC2D_ETA_FACE(iregion,:) = NSPEC2D_ETA(iregion)
    NGLOB1D_RADIAL_CORNER(iregion,:) = NGLOB1D_RADIAL(iregion)
  enddo

  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      if (mod(iproc_xi_slice(myrank),2) == 0) then
        if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
        else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
        endif
      else
        if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,3)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,3)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,3)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,3)*(NGLLZ-1))
        else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,4)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,4)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,4)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,4)*(NGLLZ-1))
        endif
      endif
    else
      if (mod(iproc_xi_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
      else
        NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
        NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
        NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
        NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                      + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
      endif
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      if (mod(iproc_eta_slice(myrank),2) == 0) then
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,1)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,1)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,1)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,1)*(NGLLZ-1))
      else
          NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NSPEC1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) + DIFF_NSPEC1D_RADIAL(:,2)
          NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_XI_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_XI(:,2)
          NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) = NSPEC2D_ETA_FACE(IREGION_OUTER_CORE,:) + DIFF_NSPEC2D_ETA(:,2)
          NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) = NGLOB1D_RADIAL_CORNER(IREGION_OUTER_CORE,:) &
                                                        + (DIFF_NSPEC1D_RADIAL(:,2)*(NGLLZ-1))
      endif
    endif
  endif

  end subroutine setup_counters

