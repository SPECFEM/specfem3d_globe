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

  subroutine count_elements(NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NPROC,&
                        NEX_PER_PROC_ETA,ratio_divide_central_cube,&
                        NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                        NSPEC1D_RADIAL, &
                        NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                        ner,ratio_sampling_array,this_region_has_a_doubling, &
                        ifirst_region,ilast_region,iter_region,iter_layer, &
                        doubling,tmp_sum,tmp_sum_xi,tmp_sum_eta, &
                        NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
                        nb_lay_sb, nspec_sb, nglob_surf, &
                        CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA, INCLUDE_CENTRAL_CUBE, &
                        last_doubling_layer, &
                        DIFF_NSPEC1D_RADIAL,DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA,&
                        tmp_sum_nglob2D_xi, tmp_sum_nglob2D_eta,divider,nglob_edges_h,&
                        nglob_edge_v,to_remove)

  use constants

  implicit none

  ! parameters to be computed based upon parameters above read from file
  integer NPROC,NEX_XI,NEX_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,ratio_divide_central_cube

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC,NSPEC2D_XI,NSPEC2D_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX


  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array


  integer :: ifirst_region, ilast_region, iter_region, iter_layer, doubling, tmp_sum, tmp_sum_xi, tmp_sum_eta
  integer ::  NUMBER_OF_MESH_LAYERS,layer_offset,nspec2D_xi_sb,nspec2D_eta_sb, &
              nb_lay_sb, nspec_sb, nglob_surf


  ! for the cut doublingbrick improvement
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  logical :: INCLUDE_CENTRAL_CUBE
  integer :: last_doubling_layer
  integer, dimension(NB_SQUARE_CORNERS,NB_CUT_CASE) :: DIFF_NSPEC1D_RADIAL
  integer, dimension(NB_SQUARE_EDGES_ONEDIR,NB_CUT_CASE) :: DIFF_NSPEC2D_XI,DIFF_NSPEC2D_ETA

  integer :: tmp_sum_nglob2D_xi, tmp_sum_nglob2D_eta,divider,nglob_edges_h,nglob_edge_v,to_remove

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  calculation of number of elements (NSPEC) below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ratio_divide_central_cube = maxval(ratio_sampling_array(1:NUMBER_OF_MESH_LAYERS))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  1D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! theoretical number of spectral elements in radial direction
  do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if (iter_region == IREGION_CRUST_MANTLE) then
      ifirst_region = 1
      ilast_region = 10 + layer_offset
    else if (iter_region == IREGION_OUTER_CORE) then
      ifirst_region = 11 + layer_offset
      ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if (iter_region == IREGION_INNER_CORE) then
      ifirst_region = NUMBER_OF_MESH_LAYERS
      ilast_region = NUMBER_OF_MESH_LAYERS
    else
      stop 'incorrect region code detected'
    endif
    NSPEC1D_RADIAL(iter_region) = sum(ner(ifirst_region:ilast_region))
  enddo

  ! difference of radial number of element for outer core if the superbrick is cut
  DIFF_NSPEC1D_RADIAL(:,:) = 0
  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC1D_RADIAL(2,1) = 1
      DIFF_NSPEC1D_RADIAL(3,1) = 2
      DIFF_NSPEC1D_RADIAL(4,1) = 1

      DIFF_NSPEC1D_RADIAL(1,2) = 1
      DIFF_NSPEC1D_RADIAL(2,2) = 2
      DIFF_NSPEC1D_RADIAL(3,2) = 1

      DIFF_NSPEC1D_RADIAL(1,3) = 1
      DIFF_NSPEC1D_RADIAL(3,3) = 1
      DIFF_NSPEC1D_RADIAL(4,3) = 2

      DIFF_NSPEC1D_RADIAL(1,4) = 2
      DIFF_NSPEC1D_RADIAL(2,4) = 1
      DIFF_NSPEC1D_RADIAL(4,4) = 1
    else
      DIFF_NSPEC1D_RADIAL(2,1) = 1
      DIFF_NSPEC1D_RADIAL(3,1) = 1

      DIFF_NSPEC1D_RADIAL(1,2) = 1
      DIFF_NSPEC1D_RADIAL(4,2) = 1
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC1D_RADIAL(3,1) = 1
      DIFF_NSPEC1D_RADIAL(4,1) = 1

      DIFF_NSPEC1D_RADIAL(1,2) = 1
      DIFF_NSPEC1D_RADIAL(2,2) = 1
    endif
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  2D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! exact number of surface elements for faces along XI and ETA

  do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if (iter_region == IREGION_CRUST_MANTLE) then
      ifirst_region = 1
      ilast_region = 10 + layer_offset
    else if (iter_region == IREGION_OUTER_CORE) then
      ifirst_region = 11 + layer_offset
      ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if (iter_region == IREGION_INNER_CORE) then
      ifirst_region = NUMBER_OF_MESH_LAYERS
      ilast_region = NUMBER_OF_MESH_LAYERS
    else
      stop 'incorrect region code detected'
    endif
    tmp_sum_xi = 0
    tmp_sum_eta = 0
    tmp_sum_nglob2D_xi = 0
    tmp_sum_nglob2D_eta = 0
    do iter_layer = ifirst_region, ilast_region
      if (this_region_has_a_doubling(iter_layer)) then
        if (iter_region == IREGION_OUTER_CORE .and. iter_layer == last_doubling_layer) then
          ! simple brick
          divider = 1
          nglob_surf = 6*NGLLX**2 - 7*NGLLX + 2
          nglob_edges_h = 2*(NGLLX-1)+1 + NGLLX
          ! minimum value to be safe
          nglob_edge_v = NGLLX-2
          nb_lay_sb = 2
          nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK
          nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK
        else
          ! double brick
          divider = 2
          if (ner(iter_layer) == 1) then
            nglob_surf = 6*NGLLX**2 - 8*NGLLX + 3
            nglob_edges_h = 4*(NGLLX-1)+1 + 2*(NGLLX-1)+1
            nglob_edge_v = NGLLX-2
            nb_lay_sb = 1
            nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK_1L
            nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK_1L
          else
            nglob_surf = 8*NGLLX**2 - 11*NGLLX + 4
            nglob_edges_h = 4*(NGLLX-1)+1 + 2*(NGLLX-1)+1
            nglob_edge_v = 2*(NGLLX-1)+1 -2
            nb_lay_sb = 2
            nspec2D_xi_sb = NSPEC2D_XI_SUPERBRICK
            nspec2D_eta_sb = NSPEC2D_ETA_SUPERBRICK
            divider = 2
          endif
        endif
        doubling = 1
        to_remove = 1
      else
        if (iter_layer /= ifirst_region) then
          to_remove = 0
        else
          to_remove = 1
        endif
        ! dummy values to avoid a warning
        nglob_surf = 0
        nglob_edges_h = 0
        nglob_edge_v = 0
        divider = 1
        doubling = 0
        nb_lay_sb = 0
        nspec2D_xi_sb = 0
        nspec2D_eta_sb = 0
      endif

      tmp_sum_xi = tmp_sum_xi + ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * (nspec2D_xi_sb/2))

      tmp_sum_eta = tmp_sum_eta + ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * (nspec2D_eta_sb/2))

      tmp_sum_nglob2D_xi = tmp_sum_nglob2D_xi + (((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb))*NGLLX*NGLLX) - &
                ((((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - doubling*nb_lay_sb)) + &
                ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))*(ner(iter_layer) - to_remove - doubling*nb_lay_sb))*NGLLX) + &
                (((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - to_remove - doubling*nb_lay_sb)) + &
                doubling * (((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))/divider) * (nglob_surf-nglob_edges_h) - &
                ((NEX_PER_PROC_XI / ratio_sampling_array(iter_layer))/divider -1) * nglob_edge_v)

      tmp_sum_nglob2D_eta = tmp_sum_nglob2D_eta + (((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb))*NGLLX*NGLLX) - &
                ((((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - doubling*nb_lay_sb)) + &
                ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))* &
                   (ner(iter_layer) - to_remove - doubling*nb_lay_sb))*NGLLX) + &
                (((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))-1)*(ner(iter_layer) - to_remove - doubling*nb_lay_sb)) + &
                doubling * (((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))/divider) * (nglob_surf-nglob_edges_h) - &
                ((NEX_PER_PROC_ETA / ratio_sampling_array(iter_layer))/divider -1) * nglob_edge_v)

    enddo ! iter_layer

    NSPEC2D_XI(iter_region) = tmp_sum_xi
    NSPEC2D_ETA(iter_region) = tmp_sum_eta

    NGLOB2DMAX_YMIN_YMAX(iter_region) = tmp_sum_nglob2D_xi
    NGLOB2DMAX_XMIN_XMAX(iter_region) = tmp_sum_nglob2D_eta

    if (iter_region == IREGION_INNER_CORE .and. INCLUDE_CENTRAL_CUBE) then
      NSPEC2D_XI(iter_region) = NSPEC2D_XI(iter_region) + &
          ((NEX_PER_PROC_XI / ratio_divide_central_cube)*(NEX_XI / ratio_divide_central_cube))
      NSPEC2D_ETA(iter_region) = NSPEC2D_ETA(iter_region) + &
          ((NEX_PER_PROC_ETA / ratio_divide_central_cube)*(NEX_XI / ratio_divide_central_cube))

      NGLOB2DMAX_YMIN_YMAX(iter_region) = NGLOB2DMAX_YMIN_YMAX(iter_region) + &
          (((NEX_PER_PROC_XI / ratio_divide_central_cube)*(NGLLX-1)+1)*((NEX_XI / ratio_divide_central_cube)*(NGLLX-1)+1))

      NGLOB2DMAX_XMIN_XMAX(iter_region) = NGLOB2DMAX_XMIN_XMAX(iter_region) + &
          (((NEX_PER_PROC_ETA / ratio_divide_central_cube)*(NGLLX-1)+1)*((NEX_XI / ratio_divide_central_cube)*(NGLLX-1)+1))
    endif
  enddo ! iter_region

  ! difference of number of surface elements along xi or eta for outer core if the superbrick is cut
  DIFF_NSPEC2D_XI(:,:) = 0
  DIFF_NSPEC2D_ETA(:,:) = 0
  if (CUT_SUPERBRICK_XI) then
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC2D_XI(2,1) = 2
      DIFF_NSPEC2D_XI(1,2) = 2
      DIFF_NSPEC2D_XI(2,3) = 2
      DIFF_NSPEC2D_XI(1,4) = 2

      DIFF_NSPEC2D_ETA(2,1) = 1
      DIFF_NSPEC2D_ETA(2,2) = 1
      DIFF_NSPEC2D_ETA(1,3) = 1
      DIFF_NSPEC2D_ETA(1,4) = 1
    else
      DIFF_NSPEC2D_ETA(2,1) = 1
      DIFF_NSPEC2D_ETA(1,2) = 1
    endif
  else
    if (CUT_SUPERBRICK_ETA) then
      DIFF_NSPEC2D_XI(2,1) = 2
      DIFF_NSPEC2D_XI(1,2) = 2
    endif
  endif
  DIFF_NSPEC2D_XI(:,:) = DIFF_NSPEC2D_XI(:,:) * (NEX_PER_PROC_XI / ratio_divide_central_cube)
  DIFF_NSPEC2D_ETA(:,:) = DIFF_NSPEC2D_ETA(:,:) * (NEX_PER_PROC_ETA / ratio_divide_central_cube)

! exact number of surface elements on the bottom and top boundaries

  ! in the crust and mantle
  NSPEC2D_TOP(IREGION_CRUST_MANTLE) = (NEX_XI/ratio_sampling_array(1))*(NEX_ETA/ratio_sampling_array(1))/NPROC
  NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE) = (NEX_XI/ratio_sampling_array(10+layer_offset))*&
                                         (NEX_ETA/ratio_sampling_array(10+layer_offset))/NPROC

  ! in the outer core with mesh doubling
  if (ADD_4TH_DOUBLING) then
    NSPEC2D_TOP(IREGION_OUTER_CORE) = (NEX_XI/(ratio_divide_central_cube/4))*(NEX_ETA/(ratio_divide_central_cube/4))/NPROC
    NSPEC2D_BOTTOM(IREGION_OUTER_CORE) = (NEX_XI/ratio_divide_central_cube)*(NEX_ETA/ratio_divide_central_cube)/NPROC
  else
    NSPEC2D_TOP(IREGION_OUTER_CORE) = (NEX_XI/(ratio_divide_central_cube/2))*(NEX_ETA/(ratio_divide_central_cube/2))/NPROC
    NSPEC2D_BOTTOM(IREGION_OUTER_CORE) = (NEX_XI/ratio_divide_central_cube)*(NEX_ETA/ratio_divide_central_cube)/NPROC
  endif

  ! in the top of the inner core
  NSPEC2D_TOP(IREGION_INNER_CORE) = (NEX_XI/ratio_divide_central_cube)*(NEX_ETA/ratio_divide_central_cube)/NPROC
  NSPEC2D_BOTTOM(IREGION_INNER_CORE) = NSPEC2D_TOP(IREGION_INNER_CORE)

  ! maximum number of surface elements on vertical boundaries of the slices
  NSPEC2DMAX_XMIN_XMAX(:) = NSPEC2D_ETA(:)
  NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE) = NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE) + maxval(DIFF_NSPEC2D_ETA(:,:))
  NSPEC2DMAX_YMIN_YMAX(:) = NSPEC2D_XI(:)
  NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE) = NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE) + maxval(DIFF_NSPEC2D_XI(:,:))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  3D case
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! exact number of spectral elements in each region

  do iter_region = IREGION_CRUST_MANTLE,IREGION_INNER_CORE
    if (iter_region == IREGION_CRUST_MANTLE) then
        ifirst_region = 1
        ilast_region = 10 + layer_offset
    else if (iter_region == IREGION_OUTER_CORE) then
        ifirst_region = 11 + layer_offset
        ilast_region = NUMBER_OF_MESH_LAYERS - 1
    else if (iter_region == IREGION_INNER_CORE) then
        ifirst_region = NUMBER_OF_MESH_LAYERS
        ilast_region = NUMBER_OF_MESH_LAYERS
    else
        stop 'incorrect region code detected'
    endif
    tmp_sum = 0;
    do iter_layer = ifirst_region, ilast_region
      if (this_region_has_a_doubling(iter_layer)) then
        if (ner(iter_layer) == 1) then
          nb_lay_sb = 1
          nspec_sb = NSPEC_SUPERBRICK_1L
        else
          nb_lay_sb = 2
          nspec_sb = NSPEC_DOUBLING_SUPERBRICK
        endif
        doubling = 1
      else
        doubling = 0
        nb_lay_sb = 0
        nspec_sb = 0
      endif
      tmp_sum = tmp_sum + (((NEX_XI / ratio_sampling_array(iter_layer)) * (NEX_ETA / ratio_sampling_array(iter_layer)) * &
                (ner(iter_layer) - doubling*nb_lay_sb)) + &
                doubling * ((NEX_XI / ratio_sampling_array(iter_layer)) * (NEX_ETA / ratio_sampling_array(iter_layer)) * &
                (nspec_sb/4))) / NPROC
    enddo
    NSPEC(iter_region) = tmp_sum
  enddo

  if (INCLUDE_CENTRAL_CUBE) NSPEC(IREGION_INNER_CORE) = NSPEC(IREGION_INNER_CORE) + &
         (NEX_PER_PROC_XI / ratio_divide_central_cube) * &
         (NEX_PER_PROC_ETA / ratio_divide_central_cube) * &
         (NEX_XI / ratio_divide_central_cube)

  if (minval(NSPEC) <= 0) stop 'negative NSPEC, there is a problem somewhere, try to recompile :) '

  end subroutine count_elements
