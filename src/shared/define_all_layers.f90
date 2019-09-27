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
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine define_all_layers(NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer, &
                               ielem,elem_doubling_mantle,elem_doubling_middle_outer_core, &
                               elem_doubling_bottom_outer_core, &
                               DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                               DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,zval, &
                               rmins,rmaxs)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!
!!!!!!  definition of general mesh parameters below
!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use constants, only: ADD_4TH_DOUBLING,MAX_NUMBER_OF_MESH_LAYERS,SUPPRESS_CRUSTAL_MESH, &
    DEPTH_FOURTH_DOUBLING_OPTIMAL,DEPTH_SECOND_DOUBLING_OPTIMAL,DEPTH_THIRD_DOUBLING_OPTIMAL, &
    HUGEVAL,ONE,R_EARTH, &
    IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80,IFLAG_670_220,IFLAG_MANTLE_NORMAL, &
    IFLAG_OUTER_CORE_NORMAL,IFLAG_INNER_CORE_NORMAL

  use shared_parameters, only: ner_mesh_layers, &
    ratio_sampling_array,this_region_has_a_doubling,doubling_index,r_bottom,r_top

  use shared_parameters, only: &
    NER_CRUST,NER_80_MOHO,NER_220_80, &
    NER_400_220,NER_600_400,NER_670_600,NER_771_670, &
    NER_TOPDDOUBLEPRIME_771,NER_CMB_TOPDDOUBLEPRIME,NER_OUTER_CORE, &
    NER_TOP_CENTRAL_CUBE_ICB, &
    RMIDDLE_CRUST,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
    R_CENTRAL_CUBE,RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER, &
    ONE_CRUST


  implicit none

  ! layers
  integer :: NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer

  ! doubling elements
  integer :: ielem,elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                      DEPTH_FOURTH_DOUBLING_REAL
  double precision :: distance,distance_min,zval

  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  ! note: all these parameters must be set to .true. for now,
  !       otherwise mesher will fail since it assumes to have at least 3 doubling layers
  logical, parameter :: ADD_1ST_DOUBLING = .true.
  logical, parameter :: ADD_2ND_DOUBLING = .true.
  logical, parameter :: ADD_3RD_DOUBLING = .true.

! find element below top of which we should implement the second doubling in the mantle
  DEPTH_SECOND_DOUBLING_REAL = 0.d0
  if (ADD_2ND_DOUBLING) then
    ! locate element closest to optimal value
    elem_doubling_mantle = -1
    distance_min = HUGEVAL
    do ielem = 2,NER_TOPDDOUBLEPRIME_771
      zval = RTOPDDOUBLEPRIME + ielem * (R771 - RTOPDDOUBLEPRIME) / dble(NER_TOPDDOUBLEPRIME_771)
      distance = abs(zval - (R_EARTH - DEPTH_SECOND_DOUBLING_OPTIMAL))
      if (distance < distance_min) then
        elem_doubling_mantle = ielem
        distance_min = distance
        DEPTH_SECOND_DOUBLING_REAL = R_EARTH - zval
      endif
    enddo
    if (elem_doubling_mantle == -1) stop 'Unable to determine second doubling element'
  endif

! find element below top of which we should implement the third doubling in the middle of the outer core
  DEPTH_THIRD_DOUBLING_REAL = 0.d0
  if (ADD_3RD_DOUBLING) then
    ! locate element closest to optimal value
    elem_doubling_middle_outer_core = -1
    distance_min = HUGEVAL
    ! start at element number 4 because we need at least two elements below for the fourth doubling
    ! implemented at the bottom of the outer core
    do ielem = 4,NER_OUTER_CORE
      zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
      distance = abs(zval - (R_EARTH - DEPTH_THIRD_DOUBLING_OPTIMAL))
      if (distance < distance_min) then
        elem_doubling_middle_outer_core = ielem
        distance_min = distance
        DEPTH_THIRD_DOUBLING_REAL = R_EARTH - zval
      endif
    enddo
    if (elem_doubling_middle_outer_core == -1) stop 'Unable to determine third doubling element'
  endif

  if (ADD_4TH_DOUBLING) then
! find element below top of which we should implement the fourth doubling in the middle of the outer core
! locate element closest to optimal value
    elem_doubling_bottom_outer_core = -1
    DEPTH_FOURTH_DOUBLING_REAL = 0.d0
    distance_min = HUGEVAL
! end two elements before the top because we need at least two elements above for the third doubling
! implemented in the middle of the outer core
    do ielem = 2,NER_OUTER_CORE-2
      zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
      distance = abs(zval - (R_EARTH - DEPTH_FOURTH_DOUBLING_OPTIMAL))
      if (distance < distance_min) then
        elem_doubling_bottom_outer_core = ielem
        distance_min = distance
        DEPTH_FOURTH_DOUBLING_REAL = R_EARTH - zval
      endif
    enddo
    if (elem_doubling_bottom_outer_core == -1) stop 'Unable to determine fourth doubling element'
! make sure that the two doublings in the outer core are found in the right order
    if (elem_doubling_bottom_outer_core >= elem_doubling_middle_outer_core) &
                    stop 'Error in location of the two doublings in the outer core'
  endif

  ratio_sampling_array(15) = 0
  last_doubling_layer = 0

! define all the layers of the mesh
  if (.not. ADD_4TH_DOUBLING) then

    ! default case:
    !     no fourth doubling at the bottom of the outer core

    if (SUPPRESS_CRUSTAL_MESH) then

      ! suppress the crustal layers
      ! will be replaced by an extension of the mantle: R_EARTH is not modified,
      ! but no more crustal doubling

      NUMBER_OF_MESH_LAYERS = 14
      layer_offset = 1

      ! now only one region
      ner_mesh_layers( 1) = NER_CRUST + NER_80_MOHO
      ner_mesh_layers( 2) = 0
      ner_mesh_layers( 3) = 0

      ner_mesh_layers( 4) = NER_220_80
      ner_mesh_layers( 5) = NER_400_220
      ner_mesh_layers( 6) = NER_600_400
      ner_mesh_layers( 7) = NER_670_600
      ner_mesh_layers( 8) = NER_771_670
      ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner_mesh_layers(10) = elem_doubling_mantle
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core
      ner_mesh_layers(14) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:9) = 1
      if (ADD_1ST_DOUBLING) then
        ratio_sampling_array(10:12) = 2
      else
        ratio_sampling_array(10:12) = ratio_sampling_array(9)
      endif
      if (ADD_2ND_DOUBLING) then
        ratio_sampling_array(13:14) = 4 ! 2 * ratio_sampling_array(12)
      else
        ratio_sampling_array(13:14) = ratio_sampling_array(12)
      endif

      ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:3) = IFLAG_CRUST !!!!! IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:13) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(14) = IFLAG_INNER_CORE_NORMAL

      ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      if (ADD_1ST_DOUBLING) then
        this_region_has_a_doubling(10) = .true.
        last_doubling_layer = 10
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth
      r_top(1) = R_EARTH
      r_bottom(1) = R80_FICTITIOUS_IN_MESHER

      r_top(2) = RMIDDLE_CRUST    !!!! now fictitious
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = R80_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH    !!!! now fictitious
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:13) = RCMB / R_EARTH
      rmins(12:13) = RICB / R_EARTH

      rmaxs(14) = RICB / R_EARTH
      rmins(14) = R_CENTRAL_CUBE / R_EARTH

    else if (ONE_CRUST) then

      ! 1D models:
      ! in order to increase CFL stability of the time scheme and therefore to allow cheaper
      ! simulations (larger time step), 1D models can be run with just one average crustal
      ! layer instead of two.

      !daniel debug
      !print *,'one_crust case in define_all_layers'

      NUMBER_OF_MESH_LAYERS = 13
      layer_offset = 0

      ner_mesh_layers( 1) = NER_CRUST
      ner_mesh_layers( 2) = NER_80_MOHO
      ner_mesh_layers( 3) = NER_220_80
      ner_mesh_layers( 4) = NER_400_220
      ner_mesh_layers( 5) = NER_600_400
      ner_mesh_layers( 6) = NER_670_600
      ner_mesh_layers( 7) = NER_771_670
      ner_mesh_layers( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner_mesh_layers( 9) = elem_doubling_mantle
      ner_mesh_layers(10) = NER_CMB_TOPDDOUBLEPRIME
      ner_mesh_layers(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(12) = elem_doubling_middle_outer_core
      ner_mesh_layers(13) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1) = 1
      if (ADD_1ST_DOUBLING) then
        ratio_sampling_array(2:8) = 2
      else
        ratio_sampling_array(2:8) = ratio_sampling_array(1)
      endif
      if (ADD_2ND_DOUBLING) then
        ratio_sampling_array(9:11) = 4 ! 2 * ratio_sampling_array(8)
      else
        ratio_sampling_array(9:11) = ratio_sampling_array(8)
      endif
      if (ADD_3RD_DOUBLING) then
        ratio_sampling_array(12:13) = 8 ! 2 * ratio_sampling_array(11)
      else
        ratio_sampling_array(12:13) = ratio_sampling_array(11)
      endif

      ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1) = IFLAG_CRUST
      doubling_index(2) = IFLAG_80_MOHO
      doubling_index(3) = IFLAG_220_80
      doubling_index(4:6) = IFLAG_670_220
      doubling_index(7:10) = IFLAG_MANTLE_NORMAL
      doubling_index(11:12) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(13) = IFLAG_INNER_CORE_NORMAL

      ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      if (ADD_1ST_DOUBLING) then
        this_region_has_a_doubling(2)  = .true.
        last_doubling_layer = 2
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(9)  = .true.
        last_doubling_layer = 9
      endif
      if (ADD_3RD_DOUBLING) then
        this_region_has_a_doubling(12) = .true.
        last_doubling_layer = 12
      endif

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth

      !!!!!!!!!!! DK DK: beware, is there a bug when 3D crust crosses anisotropy in the mantle?
      !!!!!!!!!!! DK DK: i.e. if there is no thick crust there, some elements above the Moho
      !!!!!!!!!!! DK DK: should be anisotropic but anisotropy is currently only
      !!!!!!!!!!! DK DK: stored between d220 and MOHO to save memory? Clarify this one day.
      !!!!!!!!!!! DK DK: The Moho stretching and squishing that Jeroen added to V4.0
      !!!!!!!!!!! DK DK: should partly deal with this problem.

      r_top(1) = R_EARTH
      r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(2) = R80_FICTITIOUS_IN_MESHER

      r_top(3) = R80_FICTITIOUS_IN_MESHER
      r_bottom(3) = R220

      r_top(4) = R220
      r_bottom(4) = R400

      r_top(5) = R400
      r_bottom(5) = R600

      r_top(6) = R600
      r_bottom(6) = R670

      r_top(7) = R670
      r_bottom(7) = R771

      r_top(8) = R771
      r_bottom(8) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(9) = RTOPDDOUBLEPRIME

      r_top(10) = RTOPDDOUBLEPRIME
      r_bottom(10) = RCMB

      r_top(11) = RCMB
      r_bottom(11) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(12) = RICB

      r_top(13) = RICB
      r_bottom(13) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(2) = R80_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(3) = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(3) = R220 / R_EARTH

      rmaxs(4) = R220 / R_EARTH
      rmins(4) = R400 / R_EARTH

      rmaxs(5) = R400 / R_EARTH
      rmins(5) = R600 / R_EARTH

      rmaxs(6) = R600 / R_EARTH
      rmins(6) = R670 / R_EARTH

      rmaxs(7) = R670 / R_EARTH
      rmins(7) = R771 / R_EARTH

      rmaxs(8:9) = R771 / R_EARTH
      rmins(8:9) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(10) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(10) = RCMB / R_EARTH

      rmaxs(11:12) = RCMB / R_EARTH
      rmins(11:12) = RICB / R_EARTH

      rmaxs(13) = RICB / R_EARTH
      rmins(13) = R_CENTRAL_CUBE / R_EARTH

    else

      ! default case for 3D models:
      !   contains the crustal layers
      !   doubling at the base of the crust

      !daniel debug
      !print *,'default case in define_all_layers'

      NUMBER_OF_MESH_LAYERS = 14
      layer_offset = 1
      if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER) < (R_EARTH-RMIDDLE_CRUST)) then
        ner_mesh_layers( 1) = ceiling (NER_CRUST / 2.d0)
        ner_mesh_layers( 2) = floor (NER_CRUST / 2.d0)
      else
        ner_mesh_layers( 1) = floor (NER_CRUST / 2.d0)      ! regional mesh: ner_mesh_layers(1) = 1 since NER_CRUST=3
        ner_mesh_layers( 2) = ceiling (NER_CRUST / 2.d0)    !                ner_mesh_layers(2) = 2
      endif
      ner_mesh_layers( 3) = NER_80_MOHO
      ner_mesh_layers( 4) = NER_220_80
      ner_mesh_layers( 5) = NER_400_220
      ner_mesh_layers( 6) = NER_600_400
      ner_mesh_layers( 7) = NER_670_600
      ner_mesh_layers( 8) = NER_771_670
      ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner_mesh_layers(10) = elem_doubling_mantle
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core
      ner_mesh_layers(14) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:2) = 1
      if (ADD_1ST_DOUBLING) then
        ratio_sampling_array(3:9) = 2
      else
        ratio_sampling_array(3:9) = ratio_sampling_array(2)
      endif
      if (ADD_2ND_DOUBLING) then
        ratio_sampling_array(10:12) = 4 ! 2 * ratio_sampling_array(9)
      else
        ratio_sampling_array(10:12) = ratio_sampling_array(9)
      endif
      if (ADD_3RD_DOUBLING) then
        ratio_sampling_array(13:14) = 8 ! 2 * ratio_sampling_array(12)
      else
        ratio_sampling_array(13:14) = ratio_sampling_array(12)
      endif

      ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:2) = IFLAG_CRUST
      doubling_index(3) = IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:13) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(14) = IFLAG_INNER_CORE_NORMAL

      ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      if (ADD_1ST_DOUBLING) then
        this_region_has_a_doubling(3)  = .true.
        last_doubling_layer = 3
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(10) = .true.
        last_doubling_layer = 10
      endif
      if (ADD_3RD_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif
      this_region_has_a_doubling(14) = .false.

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = RMIDDLE_CRUST

      r_top(2) = RMIDDLE_CRUST
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMIDDLE_CRUST / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:13) = RCMB / R_EARTH
      rmins(12:13) = RICB / R_EARTH

      rmaxs(14) = RICB / R_EARTH
      rmins(14) = R_CENTRAL_CUBE / R_EARTH

    endif
  else

    ! 4th doubling case:
    !     includes fourth doubling at the bottom of the outer core

    if (SUPPRESS_CRUSTAL_MESH) then

      ! suppress the crustal layers
      ! will be replaced by an extension of the mantle: R_EARTH is not modified,
      ! but no more crustal doubling

      NUMBER_OF_MESH_LAYERS = 15
      layer_offset = 1

      ! now only one region
      ner_mesh_layers( 1) = NER_CRUST + NER_80_MOHO
      ner_mesh_layers( 2) = 0
      ner_mesh_layers( 3) = 0

      ner_mesh_layers( 4) = NER_220_80
      ner_mesh_layers( 5) = NER_400_220
      ner_mesh_layers( 6) = NER_600_400
      ner_mesh_layers( 7) = NER_670_600
      ner_mesh_layers( 8) = NER_771_670
      ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner_mesh_layers(10) = elem_doubling_mantle
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner_mesh_layers(14) = elem_doubling_bottom_outer_core
      ner_mesh_layers(15) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:9) = 1
      if (ADD_1ST_DOUBLING) then
        ratio_sampling_array(10:12) = 2
      else
        ratio_sampling_array(10:12) = ratio_sampling_array(9)
      endif
      if (ADD_2ND_DOUBLING) then
        ratio_sampling_array(13) = 4 ! 2 * ratio_sampling_array(12)
      else
        ratio_sampling_array(13) = ratio_sampling_array(12)
      endif
      if (ADD_3RD_DOUBLING) then
        ratio_sampling_array(14:15) = 8 ! 2 * ratio_sampling_array(13)
      else
        ratio_sampling_array(14:15) = ratio_sampling_array(13)
      endif

      ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:3) = IFLAG_CRUST !!!!! IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:14) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(15) = IFLAG_INNER_CORE_NORMAL

      ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      if (ADD_1ST_DOUBLING) then
        this_region_has_a_doubling(10) = .true.
        last_doubling_layer = 10
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif
      if (ADD_3RD_DOUBLING) then
        this_region_has_a_doubling(14) = .true.
        last_doubling_layer = 14
      endif

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = R80_FICTITIOUS_IN_MESHER

      r_top(2) = RMIDDLE_CRUST    !!!! now fictitious
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

      r_top(14) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(14) = RICB

      r_top(15) = RICB
      r_bottom(15) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = R80_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH    !!!! now fictitious
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_EARTH    !!!! now fictitious

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:14) = RCMB / R_EARTH
      rmins(12:14) = RICB / R_EARTH

      rmaxs(15) = RICB / R_EARTH
      rmins(15) = R_CENTRAL_CUBE / R_EARTH

    else if (ONE_CRUST) then

      ! 1D models:
      ! in order to increase CFL stability of the time scheme and therefore to allow cheaper
      ! simulations (larger time step), 1D models can be run with just one average crustal
      ! layer instead of two.

      NUMBER_OF_MESH_LAYERS = 14
      layer_offset = 0

      ner_mesh_layers( 1) = NER_CRUST
      ner_mesh_layers( 2) = NER_80_MOHO
      ner_mesh_layers( 3) = NER_220_80
      ner_mesh_layers( 4) = NER_400_220
      ner_mesh_layers( 5) = NER_600_400
      ner_mesh_layers( 6) = NER_670_600
      ner_mesh_layers( 7) = NER_771_670
      ner_mesh_layers( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner_mesh_layers( 9) = elem_doubling_mantle
      ner_mesh_layers(10) = NER_CMB_TOPDDOUBLEPRIME
      ner_mesh_layers(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(12) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner_mesh_layers(13) = elem_doubling_bottom_outer_core
      ner_mesh_layers(14) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1) = 1
      if (ADD_1ST_DOUBLING) then
        ratio_sampling_array(2:8) = 2
      else
        ratio_sampling_array(2:8) = ratio_sampling_array(1)
      endif
      if (ADD_2ND_DOUBLING) then
        ratio_sampling_array(9:11) = 4 ! 2 * ratio_sampling_array(8)
      else
        ratio_sampling_array(9:11) = ratio_sampling_array(8)
      endif
      if (ADD_3RD_DOUBLING) then
        ratio_sampling_array(12) = 8 ! 2 * ratio_sampling_array(11)
      else
        ratio_sampling_array(12) = ratio_sampling_array(11)
      endif
      if (ADD_4TH_DOUBLING) then
        ratio_sampling_array(13:14) = 16 ! 2 * ratio_sampling_array(12)
      else
        ratio_sampling_array(13:14) = ratio_sampling_array(12)
      endif

      ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1) = IFLAG_CRUST
      doubling_index(2) = IFLAG_80_MOHO
      doubling_index(3) = IFLAG_220_80
      doubling_index(4:6) = IFLAG_670_220
      doubling_index(7:10) = IFLAG_MANTLE_NORMAL
      doubling_index(11:13) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(14) = IFLAG_INNER_CORE_NORMAL

      ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      if (ADD_1ST_DOUBLING) then
        this_region_has_a_doubling(2)  = .true.
        last_doubling_layer = 2
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(9)  = .true.
        last_doubling_layer = 9
      endif
      if (ADD_3RD_DOUBLING) then
        this_region_has_a_doubling(12) = .true.
        last_doubling_layer = 12
      endif
      if (ADD_4TH_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth

      !!!!!!!!!!! DK DK: beware, is there a bug when 3D crust crosses anisotropy in the mantle?
      !!!!!!!!!!! DK DK: i.e. if there is no thick crust there, some elements above the Moho
      !!!!!!!!!!! DK DK: should be anisotropic but anisotropy is currently only
      !!!!!!!!!!! DK DK: stored between d220 and MOHO to save memory? Clarify this one day.
      !!!!!!!!!!! DK DK: The Moho stretching and squishing that Jeroen added to V4.0
      !!!!!!!!!!! DK DK: should partly deal with this problem.

      r_top(1) = R_EARTH
      r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(2) = R80_FICTITIOUS_IN_MESHER

      r_top(3) = R80_FICTITIOUS_IN_MESHER
      r_bottom(3) = R220

      r_top(4) = R220
      r_bottom(4) = R400

      r_top(5) = R400
      r_bottom(5) = R600

      r_top(6) = R600
      r_bottom(6) = R670

      r_top(7) = R670
      r_bottom(7) = R771

      r_top(8) = R771
      r_bottom(8) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(9) = RTOPDDOUBLEPRIME

      r_top(10) = RTOPDDOUBLEPRIME
      r_bottom(10) = RCMB

      r_top(11) = RCMB
      r_bottom(11) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(12) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(2) = R80_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(3) = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(3) = R220 / R_EARTH

      rmaxs(4) = R220 / R_EARTH
      rmins(4) = R400 / R_EARTH

      rmaxs(5) = R400 / R_EARTH
      rmins(5) = R600 / R_EARTH

      rmaxs(6) = R600 / R_EARTH
      rmins(6) = R670 / R_EARTH

      rmaxs(7) = R670 / R_EARTH
      rmins(7) = R771 / R_EARTH

      rmaxs(8:9) = R771 / R_EARTH
      rmins(8:9) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(10) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(10) = RCMB / R_EARTH

      rmaxs(11:13) = RCMB / R_EARTH
      rmins(11:13) = RICB / R_EARTH

      rmaxs(14) = RICB / R_EARTH
      rmins(14) = R_CENTRAL_CUBE / R_EARTH

    else

      ! for 3D models:
      !   contains the crustal layers
      !   doubling at the base of the crust

      NUMBER_OF_MESH_LAYERS = 15
      layer_offset = 1
      if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER) < (R_EARTH-RMIDDLE_CRUST)) then
        ner_mesh_layers( 1) = ceiling (NER_CRUST / 2.d0)
        ner_mesh_layers( 2) = floor (NER_CRUST / 2.d0)
      else
        ner_mesh_layers( 1) = floor (NER_CRUST / 2.d0)
        ner_mesh_layers( 2) = ceiling (NER_CRUST / 2.d0)
      endif
      ner_mesh_layers( 3) = NER_80_MOHO
      ner_mesh_layers( 4) = NER_220_80
      ner_mesh_layers( 5) = NER_400_220
      ner_mesh_layers( 6) = NER_600_400
      ner_mesh_layers( 7) = NER_670_600
      ner_mesh_layers( 8) = NER_771_670
      ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
      ner_mesh_layers(10) = elem_doubling_mantle
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner_mesh_layers(14) = elem_doubling_bottom_outer_core
      ner_mesh_layers(15) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:2) = 1
      if (ADD_1ST_DOUBLING) then
        ratio_sampling_array(3:9) = 2
      else
        ratio_sampling_array(3:9) = ratio_sampling_array(2)
      endif
      if (ADD_2ND_DOUBLING) then
        ratio_sampling_array(10:12) = 4 ! 2 * ratio_sampling_array(9)
      else
        ratio_sampling_array(10:12) = ratio_sampling_array(9)
      endif
      if (ADD_3RD_DOUBLING) then
        ratio_sampling_array(13) = 8 ! 2 * ratio_sampling_array(12)
      else
        ratio_sampling_array(13) = ratio_sampling_array(12)
      endif
      if (ADD_4TH_DOUBLING) then
        ratio_sampling_array(14:15) = 16 ! 2 * ratio_sampling_array(13)
      else
        ratio_sampling_array(14:15) = ratio_sampling_array(13)
      endif

      ! value of the doubling index flag in each radial region of the mesh
      doubling_index(1:2) = IFLAG_CRUST
      doubling_index(3) = IFLAG_80_MOHO
      doubling_index(4) = IFLAG_220_80
      doubling_index(5:7) = IFLAG_670_220
      doubling_index(8:11) = IFLAG_MANTLE_NORMAL
      doubling_index(12:14) = IFLAG_OUTER_CORE_NORMAL
      doubling_index(15) = IFLAG_INNER_CORE_NORMAL

      ! define the three regions in which we implement a mesh doubling at the top of that region
      this_region_has_a_doubling(:)  = .false.
      if (ADD_1ST_DOUBLING) then
        this_region_has_a_doubling(3)  = .true.
        last_doubling_layer = 3
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(10) = .true.
        last_doubling_layer = 10
      endif
      if (ADD_3RD_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif
      if (ADD_4TH_DOUBLING) then
        this_region_has_a_doubling(14) = .true.
        last_doubling_layer = 14
      endif

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_EARTH
      r_bottom(1) = RMIDDLE_CRUST

      r_top(2) = RMIDDLE_CRUST
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      r_top(5) = R220
      r_bottom(5) = R400

      r_top(6) = R400
      r_bottom(6) = R600

      r_top(7) = R600
      r_bottom(7) = R670

      r_top(8) = R670
      r_bottom(8) = R771

      r_top(9) = R771
      r_bottom(9) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

      r_top(10) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
      r_bottom(10) = RTOPDDOUBLEPRIME

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

      r_top(14) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(14) = RICB

      r_top(15) = RICB
      r_bottom(15) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = ONE
      rmins(1) = RMIDDLE_CRUST / R_EARTH

      rmaxs(2) = RMIDDLE_CRUST / R_EARTH
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_EARTH

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_EARTH
      rmins(4) = R220 / R_EARTH

      rmaxs(5) = R220 / R_EARTH
      rmins(5) = R400 / R_EARTH

      rmaxs(6) = R400 / R_EARTH
      rmins(6) = R600 / R_EARTH

      rmaxs(7) = R600 / R_EARTH
      rmins(7) = R670 / R_EARTH

      rmaxs(8) = R670 / R_EARTH
      rmins(8) = R771 / R_EARTH

      rmaxs(9:10) = R771 / R_EARTH
      rmins(9:10) = RTOPDDOUBLEPRIME / R_EARTH

      rmaxs(11) = RTOPDDOUBLEPRIME / R_EARTH
      rmins(11) = RCMB / R_EARTH

      rmaxs(12:14) = RCMB / R_EARTH
      rmins(12:14) = RICB / R_EARTH

      rmaxs(15) = RICB / R_EARTH
      rmins(15) = R_CENTRAL_CUBE / R_EARTH
    endif
  endif

  end subroutine define_all_layers

