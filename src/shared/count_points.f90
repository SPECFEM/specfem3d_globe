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


  subroutine count_points(NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube,&
                        NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
                        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB,&
                        nblocks_xi,nblocks_eta,ner,ratio_sampling_array,&
                        this_region_has_a_doubling,&
                        ifirst_region, ilast_region, iter_region, iter_layer, &
                        doubling, padding, tmp_sum, &
                        INCLUDE_CENTRAL_CUBE,NER_TOP_CENTRAL_CUBE_ICB,NEX_XI, &
                        NUMBER_OF_MESH_LAYERS,layer_offset, &
                        nb_lay_sb, nglob_vol, nglob_surf, nglob_edge, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, &
                        last_doubling_layer, cut_doubling, nglob_int_surf_xi, nglob_int_surf_eta,nglob_ext_surf,&
                        normal_doubling, nglob_center_edge, nglob_corner_edge, nglob_border_edge)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  calculation of number of points (NGLOB) below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants

  implicit none

! parameters read from parameter file

! parameters to be computed based upon parameters above read from file
  integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

  integer, dimension(MAX_NUM_REGIONS) :: &
      NSPEC1D_RADIAL,NGLOB1D_RADIAL, &
      NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
      NGLOB

  integer NER_TOP_CENTRAL_CUBE_ICB,NEX_XI
  integer nblocks_xi,nblocks_eta

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ifirst_region, ilast_region, iter_region, iter_layer, doubling, padding, tmp_sum
  integer ::  NUMBER_OF_MESH_LAYERS,layer_offset, &
              nb_lay_sb, nglob_vol, nglob_surf, nglob_edge

! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,INCLUDE_CENTRAL_CUBE
  integer :: last_doubling_layer, cut_doubling, nglob_int_surf_xi, nglob_int_surf_eta,nglob_ext_surf,&
              normal_doubling, nglob_center_edge, nglob_corner_edge, nglob_border_edge



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  1D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! theoretical number of Gauss-Lobatto points in radial direction
  NGLOB1D_RADIAL(:) = NSPEC1D_RADIAL(:)*(NGLLZ-1)+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  2D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! 2-D addressing and buffers for summation between slices
! we add one to number of points because of the flag after the last point
  NGLOB2DMAX_XMIN_XMAX(:) = NGLOB2DMAX_XMIN_XMAX(:) + 1
  NGLOB2DMAX_YMIN_YMAX(:) = NGLOB2DMAX_YMIN_YMAX(:) + 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  3D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! exact number of global points in each region

! initialize array
  NGLOB(:) = 0

! in the inner core (no doubling region + eventually central cube)
  if (INCLUDE_CENTRAL_CUBE) then
    NGLOB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/ratio_divide_central_cube) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/ratio_divide_central_cube) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB + NEX_XI / ratio_divide_central_cube)*(NGLLZ-1)+1)
  else
    NGLOB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/ratio_divide_central_cube) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/ratio_divide_central_cube) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB)*(NGLLZ-1)+1)
  endif

! in the crust-mantle and outercore
  do iter_region = IREGION_CRUST_MANTLE,IREGION_OUTER_CORE
      if (iter_region == IREGION_CRUST_MANTLE) then
            ifirst_region = 1
            ilast_region = 10 + layer_offset
      else if (iter_region == IREGION_OUTER_CORE) then
            ifirst_region = 11 + layer_offset
            ilast_region = NUMBER_OF_MESH_LAYERS - 1
      else
            stop 'incorrect region code detected'
      endif
      tmp_sum = 0;
      do iter_layer = ifirst_region, ilast_region
        nglob_int_surf_eta = 0
        nglob_int_surf_xi = 0
        nglob_ext_surf = 0
        nglob_center_edge = 0
        nglob_corner_edge = 0
        nglob_border_edge = 0
        if (this_region_has_a_doubling(iter_layer)) then
            if (iter_region == IREGION_OUTER_CORE .and. iter_layer == last_doubling_layer .and. &
               (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA)) then
              doubling = 1
              normal_doubling = 0
              cut_doubling = 1
              nb_lay_sb = 2
              nglob_edge = 0
              nglob_surf = 0
              nglob_vol = 8*NGLLX**3 - 12*NGLLX**2 + 6*NGLLX - 1
              nglob_int_surf_eta = 6*NGLLX**2 - 7*NGLLX + 2
              nglob_int_surf_xi = 5*NGLLX**2 - 5*NGLLX + 1
              nglob_ext_surf = 4*NGLLX**2-4*NGLLX+1
              nglob_center_edge = 4*(NGLLX-1)+1
              nglob_corner_edge = 2*(NGLLX-1)+1
              nglob_border_edge = 3*(NGLLX-1)+1
            else
              if (ner(iter_layer) == 1) then
                nb_lay_sb = 1
                nglob_vol = 28*NGLLX**3 - 62*NGLLX**2 + 47*NGLLX - 12
                nglob_surf = 6*NGLLX**2-8*NGLLX+3
                nglob_edge = NGLLX
              else
                nb_lay_sb = 2
                nglob_vol = 32*NGLLX**3 - 70*NGLLX**2 + 52*NGLLX - 13
                nglob_surf = 8*NGLLX**2-11*NGLLX+4
                nglob_edge = 2*NGLLX-1
              endif
              doubling = 1
              normal_doubling = 1
              cut_doubling = 0
            endif
            padding = -1
        else
            doubling = 0
            normal_doubling = 0
            cut_doubling = 0
            padding = 0
            nb_lay_sb = 0
            nglob_vol = 0
            nglob_surf = 0
            nglob_edge = 0
        endif
        if (iter_layer == ilast_region) padding = padding +1
        nblocks_xi = NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)
        nblocks_eta = NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)

        tmp_sum = tmp_sum + &
        ((nblocks_xi)*(NGLLX-1)+1) * ((nblocks_eta)*(NGLLX-1)+1) * ((ner(iter_layer) - doubling*nb_lay_sb)*(NGLLX-1)+padding)+&
        normal_doubling * ((((nblocks_xi*nblocks_eta)/4)*nglob_vol) - &
        (((nblocks_eta/2-1)*nblocks_xi/2+(nblocks_xi/2-1)*nblocks_eta/2)*nglob_surf) + &
        ((nblocks_eta/2-1)*(nblocks_xi/2-1)*nglob_edge)) + &
        cut_doubling*(nglob_vol*(nblocks_xi*nblocks_eta) - &
            ( nblocks_eta*(int(nblocks_xi/2)*nglob_int_surf_xi + int((nblocks_xi-1)/2)*nglob_ext_surf) + &
              nblocks_xi*(int(nblocks_eta/2)*nglob_int_surf_eta + int((nblocks_eta-1)/2)*nglob_ext_surf)&
            ) + &
            ( int(nblocks_xi/2)*int(nblocks_eta/2)*nglob_center_edge + &
              int((nblocks_xi-1)/2)*int((nblocks_eta-1)/2)*nglob_corner_edge + &
              ((int(nblocks_eta/2)*int((nblocks_xi-1)/2))+(int((nblocks_eta-1)/2)*int(nblocks_xi/2)))*nglob_border_edge&
            ))
      enddo
      NGLOB(iter_region) = tmp_sum
  enddo

!!! example :
!!!                        nblocks_xi/2=5
!!!                  ____________________________________
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!! nblocks_eta/2=3  I______+______+______+______+______I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I______+______+______+______+______I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I      I      I      I      I      I
!!!                  I______I______I______I______I______I
!!!
!!! NGLOB for this doubling layer = 3*5*Volume - ((3-1)*5+(5-1)*3)*Surface + (3-1)*(5-1)*Edge
!!!
!!! 32*NGLLX**3 - 70*NGLLX**2 + 52*NGLLX - 13 -> nb GLL points in a superbrick (Volume)
!!! 8*NGLLX**2-11*NGLLX+4 -> nb GLL points on a superbrick side (Surface)
!!! 2*NGLLX-1 -> nb GLL points on a corner edge of a superbrick (Edge)

!!! for the one layer superbrick :
!!! NGLOB = 28.NGLL^3 - 62.NGLL^2 + 47.NGLL - 12 (Volume)
!!! NGLOB = 6.NGLL^2 - 8.NGLL + 3 (Surface)
!!! NGLOB = NGLL (Edge)
!!!
!!! those results were obtained by using the script UTILS/doubling_brick/count_nglob_analytical.pl
!!! with an OpenDX file of the superbrick's geometry

!!! for the basic doubling bricks (two layers)
!!! NGLOB = 8.NGLL^3 - 12.NGLL^2 + 6.NGLL - 1 (VOLUME)
!!! NGLOB = 5.NGLL^2 - 5.NGLL + 1 (SURFACE 1)
!!! NGLOB = 6.NGLL^2 - 7.NGLL + 2 (SURFACE 2)

  end subroutine count_points

