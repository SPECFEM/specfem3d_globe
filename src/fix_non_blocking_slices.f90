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

!! DK DK fix the non-blocking arrays to assemble inside the chunks: elements
!! DK DK in contact with the MPI faces by an edge or a corner only but not
!! DK DK a full face are missing, therefore let us add them
  subroutine fix_non_blocking_slices(is_on_a_slice_edge,iboolright_xi, &
         iboolleft_xi,iboolright_eta,iboolleft_eta, &
         npoin2D_xi,npoin2D_eta,ibool, &
         mask_ibool,nspec,nglob,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX)

  implicit none

  include "constants.h"

  integer :: nspec,nglob,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi,npoin2D_eta

  logical, dimension(nspec) :: is_on_a_slice_edge

  integer, dimension(NGLOB2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NGLOB2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

! this mask is declared as integer in the calling program because it is used elsewhere
! to store integers, and it is reused here as a logical to save memory
  logical, dimension(nglob) :: mask_ibool

  integer :: ipoin,ispec,i,j,k

! clean the mask
  mask_ibool(:) = .false.

! mark all the points that are in the MPI buffers to assemble inside each chunk
  do ipoin = 1,npoin2D_xi(1)
    mask_ibool(iboolleft_xi(ipoin)) = .true.
  enddo

  do ipoin = 1,npoin2D_eta(1)
    mask_ibool(iboolleft_eta(ipoin)) = .true.
  enddo

  do ipoin = 1,npoin2D_xi(2)
    mask_ibool(iboolright_xi(ipoin)) = .true.
  enddo

  do ipoin = 1,npoin2D_eta(2)
    mask_ibool(iboolright_eta(ipoin)) = .true.
  enddo

! now label all the elements that have at least one corner belonging
! to any of these buffers as elements that must contribute to the
! first step of the calculations (performed on the edges before starting
! the non-blocking communications); there is no need to examine the inside
! of the elements, checking their eight corners is sufficient
  do ispec = 1,nspec
    do k = 1,NGLLZ,NGLLZ-1
      do j  = 1,NGLLY,NGLLY-1
        do i = 1,NGLLX,NGLLX-1
          if(mask_ibool(ibool(i,j,k,ispec))) then
            is_on_a_slice_edge(ispec) = .true.
            goto 888
          endif
        enddo
      enddo
    enddo
  888 continue
  enddo

  end subroutine fix_non_blocking_slices

