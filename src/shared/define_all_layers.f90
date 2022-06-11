!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  use constants, only: myrank,ADD_4TH_DOUBLING,MAX_NUMBER_OF_MESH_LAYERS,SUPPRESS_CRUSTAL_MESH,HUGEVAL,R_UNIT_SPHERE, &
    IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80,IFLAG_670_220,IFLAG_MANTLE_NORMAL, &
    IFLAG_OUTER_CORE_NORMAL,IFLAG_INNER_CORE_NORMAL

  use shared_parameters, only: PLANET_TYPE,IPLANET_MARS,IPLANET_EARTH,IPLANET_MOON, &
    R_PLANET

  ! Earth
  use constants, only: &
    EARTH_DEPTH_SECOND_DOUBLING_OPTIMAL,EARTH_DEPTH_THIRD_DOUBLING_OPTIMAL,EARTH_DEPTH_FOURTH_DOUBLING_OPTIMAL
  ! Mars
  use constants, only: &
    MARS_DEPTH_SECOND_DOUBLING_OPTIMAL,MARS_DEPTH_THIRD_DOUBLING_OPTIMAL,MARS_DEPTH_FOURTH_DOUBLING_OPTIMAL
  ! Moon
  use constants, only: &
    MOON_DEPTH_SECOND_DOUBLING_OPTIMAL,MOON_DEPTH_THIRD_DOUBLING_OPTIMAL,MOON_DEPTH_FOURTH_DOUBLING_OPTIMAL

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

  use shared_parameters, only: REGIONAL_MESH_CUTOFF,REGIONAL_MESH_CUTOFF_DEPTH,REGIONAL_MESH_ADD_2ND_DOUBLING, &
    USE_LOCAL_MESH

  implicit none

  ! layers
  integer :: NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer

  ! doubling elements
  integer :: ielem
  integer :: elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  integer :: ner_start,ner_end,ner_layer

  ! doubling elements positions
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                      DEPTH_FOURTH_DOUBLING_REAL
  double precision :: distance,distance_min,zval
  double precision :: r_layer_top,r_layer_bottom
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  double precision :: DEPTH_SECOND_DOUBLING_OPTIMAL,DEPTH_THIRD_DOUBLING_OPTIMAL,DEPTH_FOURTH_DOUBLING_OPTIMAL

  ! doubling layer selection
  logical :: ADD_1ST_DOUBLING
  logical :: ADD_2ND_DOUBLING
  logical :: ADD_3RD_DOUBLING

  ! debugging
  logical, parameter :: DEBUG = .false.

  ! initializes
  NUMBER_OF_MESH_LAYERS = 0
  layer_offset = 0

  ner_mesh_layers(:) = 0

  this_region_has_a_doubling(:)  = .false.
  last_doubling_layer = 0

  ratio_sampling_array(:) = 0
  doubling_index(:) = 0

  r_top(:) = 0.d0
  r_bottom(:) = 0.d0

  rmins(:) = 0.d0
  rmaxs(:) = 0.d0

  elem_doubling_mantle = 0
  elem_doubling_middle_outer_core = 0
  elem_doubling_bottom_outer_core = 0

  DEPTH_SECOND_DOUBLING_REAL = 0.d0
  DEPTH_THIRD_DOUBLING_REAL = 0.d0
  DEPTH_FOURTH_DOUBLING_REAL = 0.d0

  if (REGIONAL_MESH_CUTOFF) then
    ! regional mesh cut-off
    ! skips doubling layers
    ADD_1ST_DOUBLING = .false.
    ADD_2ND_DOUBLING = .false.
    ADD_3RD_DOUBLING = .false.
    DEPTH_SECOND_DOUBLING_REAL = R_PLANET - R771  ! for placing radius of auxiliary layers to correct dpeths
    DEPTH_THIRD_DOUBLING_REAL = R_PLANET - RICB
    DEPTH_FOURTH_DOUBLING_REAL = R_PLANET - RICB

    ! first-doubling below Moho in default cases (unless SUPPRESS_CRUSTAL_MESH was set in constant.h)
    if (REGIONAL_MESH_CUTOFF_DEPTH > 80.d0) ADD_1ST_DOUBLING = .true.
    ! second-doubling will be moved to 220km depth (default would be below 771km)
    if (REGIONAL_MESH_CUTOFF_DEPTH > 220.d0 .and. REGIONAL_MESH_ADD_2ND_DOUBLING) ADD_2ND_DOUBLING = .true.
  else
    ! default mesh
    ! note: all these parameters must be set to .true. for now,
    !       otherwise mesher will fail since it assumes to have at least 3 doubling layers
    ADD_1ST_DOUBLING = .true.
    ADD_2ND_DOUBLING = .true.
    ADD_3RD_DOUBLING = .true.
  endif

  ! desired doubling depths
  select case(PLANET_TYPE)
  case (IPLANET_EARTH)
    ! Earth
    DEPTH_SECOND_DOUBLING_OPTIMAL = EARTH_DEPTH_SECOND_DOUBLING_OPTIMAL
    DEPTH_THIRD_DOUBLING_OPTIMAL = EARTH_DEPTH_THIRD_DOUBLING_OPTIMAL
    DEPTH_FOURTH_DOUBLING_OPTIMAL = EARTH_DEPTH_FOURTH_DOUBLING_OPTIMAL

  case (IPLANET_MARS)
    ! Mars
    DEPTH_SECOND_DOUBLING_OPTIMAL = MARS_DEPTH_SECOND_DOUBLING_OPTIMAL
    DEPTH_THIRD_DOUBLING_OPTIMAL = MARS_DEPTH_THIRD_DOUBLING_OPTIMAL
    DEPTH_FOURTH_DOUBLING_OPTIMAL = MARS_DEPTH_FOURTH_DOUBLING_OPTIMAL

  case (IPLANET_MOON)
    ! Moon
    DEPTH_SECOND_DOUBLING_OPTIMAL = MOON_DEPTH_SECOND_DOUBLING_OPTIMAL
    DEPTH_THIRD_DOUBLING_OPTIMAL = MOON_DEPTH_THIRD_DOUBLING_OPTIMAL
    DEPTH_FOURTH_DOUBLING_OPTIMAL = MOON_DEPTH_FOURTH_DOUBLING_OPTIMAL

  case default
    ! avoiding exit_MPI(), since we also call this routine in create_header_file
    ! which can be compiled without MPI - using stop instead
    !call exit_MPI(myrank,'Invalid planet, defining all layers not implemented yet')
    print *,'Invalid planet, defining all layers not implemented yet'
    stop 'Invalid planet, defining all layers not implemented yet'
  end select

  ! find element below top of which we should implement the second doubling in the mantle
  if (ADD_2ND_DOUBLING) then
    ! locate element closest to optimal value
    elem_doubling_mantle = -1
    distance_min = HUGEVAL
    ! doubling between D'' and 771
    ner_layer = NER_TOPDDOUBLEPRIME_771
    r_layer_top = R771
    r_layer_bottom = RTOPDDOUBLEPRIME

    ! regional mesh cutoff
    if (REGIONAL_MESH_CUTOFF) then
      ! sets doubling layer around 220km
      DEPTH_SECOND_DOUBLING_OPTIMAL = 220.d0
      ner_layer = NER_400_220
      r_layer_top = R220
      r_layer_bottom = R400
    endif

    ! finds best layer
    do ielem = 2,ner_layer
      zval = r_layer_bottom + ielem * (r_layer_top - r_layer_bottom) / dble(ner_layer)
      distance = abs(zval - (R_PLANET - DEPTH_SECOND_DOUBLING_OPTIMAL))

      ! debug
      if (DEBUG .and. myrank == 0) &
        print *,'debug: 2nd doubling',ielem,ner_layer,'dist/zval',distance,distance_min,zval

      ! checks if closer and sets as new depth
      if (distance < distance_min) then
        elem_doubling_mantle = ielem
        distance_min = distance
        DEPTH_SECOND_DOUBLING_REAL = R_PLANET - zval
      endif
    enddo

    !debug
    if (DEBUG .and. myrank == 0) &
      print *,'debug: 2nd doubling index = ',elem_doubling_mantle,DEPTH_SECOND_DOUBLING_REAL,'(in mantle D" - 771)'

    ! check if layer found
    if (REGIONAL_MESH_CUTOFF) then
      if (elem_doubling_mantle == -1) then
        ! skip doubling layer
        ADD_2ND_DOUBLING = .false.
        elem_doubling_mantle = 0
        DEPTH_SECOND_DOUBLING_REAL = (R_PLANET - R771)
      endif
    else
      if (elem_doubling_mantle == -1) stop 'Unable to determine second doubling element'
    endif
  endif

  ! find element below top of which we should implement the third doubling in the middle of the outer core
  if (ADD_3RD_DOUBLING) then
    ! locate element closest to optimal value
    elem_doubling_middle_outer_core = -1
    distance_min = HUGEVAL
    ! start at element number 4 because we need at least two elements below for the fourth doubling
    ! implemented at the bottom of the outer core
    ! also, end at least one element layer before the CMB, otherwise the MPI interface points might cause problems.
    if (ADD_4TH_DOUBLING) then
      ner_start = 4
    else
      ner_start = 2
    endif
    ner_end = NER_OUTER_CORE-1
    do ielem = ner_start,ner_end
      zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
      distance = abs(zval - (R_PLANET - DEPTH_THIRD_DOUBLING_OPTIMAL))

      ! debug
      if (DEBUG .and. myrank == 0) &
        print *,'debug: 3nd doubling',ielem,ner_end,'dist/zval',distance,distance_min,zval

      ! checks if closer and sets as new depth
      if (distance < distance_min) then
        elem_doubling_middle_outer_core = ielem
        distance_min = distance
        DEPTH_THIRD_DOUBLING_REAL = R_PLANET - zval
      endif
    enddo

    ! debug
    if (DEBUG .and. myrank == 0) &
      print *,'debug: 3rd doubling index = ',elem_doubling_middle_outer_core,DEPTH_THIRD_DOUBLING_REAL,'(in outer core)'
    ! check
    if (elem_doubling_middle_outer_core == -1) stop 'Unable to determine third doubling element'
  endif

  if (ADD_4TH_DOUBLING) then
    ! find element below top of which we should implement the fourth doubling in the middle of the outer core
    ! locate element closest to optimal value
    elem_doubling_bottom_outer_core = -1
    distance_min = HUGEVAL
    ! end two elements before the top because we need at least two elements above for the third doubling
    ! implemented in the middle of the outer core
    do ielem = 2,NER_OUTER_CORE-2
      zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
      distance = abs(zval - (R_PLANET - DEPTH_FOURTH_DOUBLING_OPTIMAL))

      ! debug
      if (DEBUG .and. myrank == 0) &
        print *,'debug: 2nd doubling',ielem,NER_OUTER_CORE-2,'dist/zval',distance,distance_min,zval

      ! checks if closer and sets as new depth
      if (distance < distance_min) then
        elem_doubling_bottom_outer_core = ielem
        distance_min = distance
        DEPTH_FOURTH_DOUBLING_REAL = R_PLANET - zval
      endif
    enddo

    !debug
    if (DEBUG .and. myrank == 0) &
      print *,'debug: 4th doubling index = ',elem_doubling_bottom_outer_core,DEPTH_FOURTH_DOUBLING_REAL,'(in outer core)'
    ! check
    if (elem_doubling_bottom_outer_core == -1) stop 'Unable to determine fourth doubling element'
    ! make sure that the two doublings in the outer core are found in the right order
    if (elem_doubling_bottom_outer_core >= elem_doubling_middle_outer_core) then
      print *,'Error: invalid 4th doubling layer with bottom doubling layer index = ',elem_doubling_bottom_outer_core
      print *,'       should be below 3rd doubling layer at ',elem_doubling_middle_outer_core
      stop 'Error in location of the two doublings in the outer core, 4th should be below 3rd doubling layer'
    endif
  endif

  ! sets total number of layers
  !   NUMBER_OF_MESH_LAYERS  -  total number of mesh layers
  !   layer_offset           -  crust/mantle region has 10 + layer_offset layers
  call define_all_layers_number_and_offset(NUMBER_OF_MESH_LAYERS,layer_offset)

  ! define all the layers of the mesh
  if (.not. ADD_4TH_DOUBLING) then

    ! default case:
    !     no fourth doubling at the bottom of the outer core

    if (SUPPRESS_CRUSTAL_MESH) then

      ! suppress the crustal layers
      ! will be replaced by an extension of the mantle: R_PLANET is not modified,
      ! but no more crustal doubling

      ! check with define_all_layers_number_and_offset()
      ! NUMBER_OF_MESH_LAYERS = 14
      ! layer_offset = 1

      ! crust/mantle
      ! now only one region
      ner_mesh_layers( 1) = NER_CRUST + NER_80_MOHO
      ner_mesh_layers( 2) = 0
      ner_mesh_layers( 3) = 0
      ner_mesh_layers( 4) = NER_220_80

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ner_mesh_layers( 5) = NER_400_220 - elem_doubling_mantle
        ner_mesh_layers( 6) = elem_doubling_mantle
        ner_mesh_layers( 7) = NER_600_400
        ner_mesh_layers( 8) = NER_670_600
        ner_mesh_layers( 9) = NER_771_670
        ner_mesh_layers(10) = NER_TOPDDOUBLEPRIME_771
      else
        ! 2nd doubling after 771
        ner_mesh_layers( 5) = NER_400_220
        ner_mesh_layers( 6) = NER_600_400
        ner_mesh_layers( 7) = NER_670_600
        ner_mesh_layers( 8) = NER_771_670
        ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
        ner_mesh_layers(10) = elem_doubling_mantle
      endif
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME

      ! outer core
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core

      ! inner core
      ner_mesh_layers(14) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ratio_sampling_array(1:5) = 1
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(6:12) = 2
        else
          ratio_sampling_array(6:12) = ratio_sampling_array(5)
        endif
      else
        ! 2nd doubling after 771
        ratio_sampling_array(1:9) = 1
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(10:12) = 2
        else
          ratio_sampling_array(10:12) = ratio_sampling_array(9)
        endif
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
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_1ST_DOUBLING) then
          this_region_has_a_doubling(6) = .true.
          last_doubling_layer = 6
        endif
      else
        ! 2nd doubling after 771
        if (ADD_1ST_DOUBLING) then
          this_region_has_a_doubling(10) = .true.
          last_doubling_layer = 10
        endif
      endif
      if (ADD_2ND_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth
      r_top(1) = R_PLANET
      r_bottom(1) = R80_FICTITIOUS_IN_MESHER

      r_top(2) = RMIDDLE_CRUST    !!!! now fictitious
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        r_top(5) = R220
        r_bottom(5) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(6) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(6) = R400

        r_top(7) = R400
        r_bottom(7) = R600

        r_top(8) = R600
        r_bottom(8) = R670

        r_top(9) = R670
        r_bottom(9) = R771

        r_top(10) = R771
        r_bottom(10) = RTOPDDOUBLEPRIME
      else
        ! 2nd doubling after 771
        r_top(5) = R220
        r_bottom(5) = R400

        r_top(6) = R400
        r_bottom(6) = R600

        r_top(7) = R600
        r_bottom(7) = R670

        r_top(8) = R670
        r_bottom(8) = R771

        r_top(9) = R771
        r_bottom(9) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(10) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(10) = RTOPDDOUBLEPRIME
      endif

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = R_UNIT_SPHERE
      rmins(1) = R80_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(2) = RMIDDLE_CRUST / R_PLANET    !!!! now fictitious
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET    !!!! now fictitious

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET    !!!! now fictitious
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET    !!!! now fictitious

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(4) = R220 / R_PLANET

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        rmaxs(5:6) = R220 / R_PLANET
        rmins(5:6) = R400 / R_PLANET

        rmaxs(7) = R400 / R_PLANET
        rmins(7) = R600 / R_PLANET

        rmaxs(8) = R600 / R_PLANET
        rmins(8) = R670 / R_PLANET

        rmaxs(9) = R670 / R_PLANET
        rmins(9) = R771 / R_PLANET

        rmaxs(10) = R771 / R_PLANET
        rmins(10) = RTOPDDOUBLEPRIME / R_PLANET
      else
        ! 2nd doubling after 771
        rmaxs(5) = R220 / R_PLANET
        rmins(5) = R400 / R_PLANET

        rmaxs(6) = R400 / R_PLANET
        rmins(6) = R600 / R_PLANET

        rmaxs(7) = R600 / R_PLANET
        rmins(7) = R670 / R_PLANET

        rmaxs(8) = R670 / R_PLANET
        rmins(8) = R771 / R_PLANET

        rmaxs(9:10) = R771 / R_PLANET
        rmins(9:10) = RTOPDDOUBLEPRIME / R_PLANET
      endif

      rmaxs(11) = RTOPDDOUBLEPRIME / R_PLANET
      rmins(11) = RCMB / R_PLANET

      rmaxs(12:13) = RCMB / R_PLANET
      rmins(12:13) = RICB / R_PLANET

      rmaxs(14) = RICB / R_PLANET
      rmins(14) = R_CENTRAL_CUBE / R_PLANET

    else if (ONE_CRUST) then

      ! 1D models:
      ! in order to increase CFL stability of the time scheme and therefore to allow cheaper
      ! simulations (larger time step), 1D models can be run with just one average crustal
      ! layer instead of two.

      !daniel debug
      !print *,'one_crust case in define_all_layers'

      ! check with define_all_layers_number_and_offset()
      ! NUMBER_OF_MESH_LAYERS = 13
      ! layer_offset = 0

      ! crust/mantle
      ner_mesh_layers( 1) = NER_CRUST
      ner_mesh_layers( 2) = NER_80_MOHO
      ner_mesh_layers( 3) = NER_220_80

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ner_mesh_layers( 4) = NER_400_220 - elem_doubling_mantle
        ner_mesh_layers( 5) = elem_doubling_mantle
        ner_mesh_layers( 6) = NER_600_400
        ner_mesh_layers( 7) = NER_670_600
        ner_mesh_layers( 8) = NER_771_670
        ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771
      else
        ! 2nd doubling after 771
        ner_mesh_layers( 4) = NER_400_220
        ner_mesh_layers( 5) = NER_600_400
        ner_mesh_layers( 6) = NER_670_600
        ner_mesh_layers( 7) = NER_771_670
        ner_mesh_layers( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
        ner_mesh_layers( 9) = elem_doubling_mantle
      endif
      ner_mesh_layers(10) = NER_CMB_TOPDDOUBLEPRIME

      ! outer core
      ner_mesh_layers(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(12) = elem_doubling_middle_outer_core

      ! inner core
      ner_mesh_layers(13) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1) = 1
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(2:4) = 2
        else
          ratio_sampling_array(2:4) = ratio_sampling_array(1)
        endif
        if (ADD_2ND_DOUBLING) then
          ratio_sampling_array(5:11) = 4 ! 2 * ratio_sampling_array(8)
        else
          ratio_sampling_array(5:11) = ratio_sampling_array(4)
        endif
      else
        ! 2nd doubling after 771
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
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(5)  = .true.
          last_doubling_layer = 5
        endif
      else
        ! 2nd doubling after 771
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(9)  = .true.
          last_doubling_layer = 9
        endif
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

      r_top(1) = R_PLANET
      r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(2) = R80_FICTITIOUS_IN_MESHER

      r_top(3) = R80_FICTITIOUS_IN_MESHER
      r_bottom(3) = R220

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        r_top(4) = R220
        r_bottom(4) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(5) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(5) = R400

        r_top(6) = R400
        r_bottom(6) = R600

        r_top(7) = R600
        r_bottom(7) = R670

        r_top(8) = R670
        r_bottom(8) = R771

        r_top(9) = R771
        r_bottom(9) = RTOPDDOUBLEPRIME
      else
        ! 2nd doubling after 771
        r_top(4) = R220
        r_bottom(4) = R400

        r_top(5) = R400
        r_bottom(5) = R600

        r_top(6) = R600
        r_bottom(6) = R670

        r_top(7) = R670
        r_bottom(7) = R771

        r_top(8) = R771
        r_bottom(8) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(9) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(9) = RTOPDDOUBLEPRIME
      endif
      r_top(10) = RTOPDDOUBLEPRIME
      r_bottom(10) = RCMB

      r_top(11) = RCMB
      r_bottom(11) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL

      r_top(12) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(12) = RICB

      r_top(13) = RICB
      r_bottom(13) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = R_UNIT_SPHERE
      rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(2) = R80_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(3) = R220 / R_PLANET

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        rmaxs(4:5) = R220 / R_PLANET
        rmins(4:5) = R400 / R_PLANET

        rmaxs(6) = R400 / R_PLANET
        rmins(6) = R600 / R_PLANET

        rmaxs(7) = R600 / R_PLANET
        rmins(7) = R670 / R_PLANET

        rmaxs(8) = R670 / R_PLANET
        rmins(8) = R771 / R_PLANET

        rmaxs(9) = R771 / R_PLANET
        rmins(9) = RTOPDDOUBLEPRIME / R_PLANET
      else
        ! 2nd doubling after 771
        rmaxs(4) = R220 / R_PLANET
        rmins(4) = R400 / R_PLANET

        rmaxs(5) = R400 / R_PLANET
        rmins(5) = R600 / R_PLANET

        rmaxs(6) = R600 / R_PLANET
        rmins(6) = R670 / R_PLANET

        rmaxs(7) = R670 / R_PLANET
        rmins(7) = R771 / R_PLANET

        rmaxs(8:9) = R771 / R_PLANET
        rmins(8:9) = RTOPDDOUBLEPRIME / R_PLANET
      endif
      rmaxs(10) = RTOPDDOUBLEPRIME / R_PLANET
      rmins(10) = RCMB / R_PLANET

      rmaxs(11:12) = RCMB / R_PLANET
      rmins(11:12) = RICB / R_PLANET

      rmaxs(13) = RICB / R_PLANET
      rmins(13) = R_CENTRAL_CUBE / R_PLANET

    else

      ! default case for 3D models:
      !   contains the crustal layers
      !   doubling at the base of the crust

      !daniel debug
      !print *,'default case in define_all_layers'

      ! check with define_all_layers_number_and_offset()
      ! NUMBER_OF_MESH_LAYERS = 14
      ! layer_offset = 1

      ! crust/mantle
      if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER) < (R_PLANET-RMIDDLE_CRUST)) then
        ner_mesh_layers( 1) = ceiling (NER_CRUST / 2.d0)
        ner_mesh_layers( 2) = floor (NER_CRUST / 2.d0)
      else
        ner_mesh_layers( 1) = floor (NER_CRUST / 2.d0)      ! regional mesh: ner_mesh_layers(1) = 1 since NER_CRUST=3
        ner_mesh_layers( 2) = ceiling (NER_CRUST / 2.d0)    !                ner_mesh_layers(2) = 2
      endif
      ner_mesh_layers( 3) = NER_80_MOHO
      ner_mesh_layers( 4) = NER_220_80

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ner_mesh_layers( 5) = NER_400_220 - elem_doubling_mantle
        ner_mesh_layers( 6) = elem_doubling_mantle
        ner_mesh_layers( 7) = NER_600_400
        ner_mesh_layers( 8) = NER_670_600
        ner_mesh_layers( 9) = NER_771_670
        ner_mesh_layers(10) = NER_TOPDDOUBLEPRIME_771
      else
        ! 2nd doubling after 771
        ner_mesh_layers( 5) = NER_400_220
        ner_mesh_layers( 6) = NER_600_400
        ner_mesh_layers( 7) = NER_670_600
        ner_mesh_layers( 8) = NER_771_670
        ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
        ner_mesh_layers(10) = elem_doubling_mantle
      endif
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME

      ! outer core
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core

      ! inner core
      ner_mesh_layers(14) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:2) = 1
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(3:5) = 2
        else
          ratio_sampling_array(3:5) = ratio_sampling_array(2)
        endif
        if (ADD_2ND_DOUBLING) then
          ratio_sampling_array(6:12) = 4 ! 2 * ratio_sampling_array(9)
        else
          ratio_sampling_array(6:12) = ratio_sampling_array(9)
        endif
      else
        ! 2nd doubling after 771
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
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(6) = .true.
          last_doubling_layer = 6
        endif
      else
        ! 2nd doubling after 771
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(10) = .true.
          last_doubling_layer = 10
        endif
      endif
      if (ADD_3RD_DOUBLING) then
        this_region_has_a_doubling(13) = .true.
        last_doubling_layer = 13
      endif
      this_region_has_a_doubling(14) = .false.

      ! define the top and bottom radii of all the regions of the mesh in the radial direction
      ! the first region is the crust at the surface of the Earth
      ! the last region is in the inner core near the center of the Earth

      r_top(1) = R_PLANET
      r_bottom(1) = RMIDDLE_CRUST

      r_top(2) = RMIDDLE_CRUST
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        r_top(5) = R220
        r_bottom(5) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(6) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(6) = R400

        r_top(7) = R400
        r_bottom(7) = R600

        r_top(8) = R600
        r_bottom(8) = R670

        r_top(9) = R670
        r_bottom(9) = R771

        r_top(10) = R771
        r_bottom(10) = RTOPDDOUBLEPRIME
      else
        ! 2nd doubling after 771
        r_top(5) = R220
        r_bottom(5) = R400

        r_top(6) = R400
        r_bottom(6) = R600

        r_top(7) = R600
        r_bottom(7) = R670

        r_top(8) = R670
        r_bottom(8) = R771

        r_top(9) = R771
        r_bottom(9) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(10) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(10) = RTOPDDOUBLEPRIME
      endif

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = R_UNIT_SPHERE
      rmins(1) = RMIDDLE_CRUST / R_PLANET

      rmaxs(2) = RMIDDLE_CRUST / R_PLANET
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(4) = R220 / R_PLANET

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        rmaxs(5:6) = R220 / R_PLANET
        rmins(5:6) = R400 / R_PLANET

        rmaxs(7) = R400 / R_PLANET
        rmins(7) = R600 / R_PLANET

        rmaxs(8) = R600 / R_PLANET
        rmins(8) = R670 / R_PLANET

        rmaxs(9) = R670 / R_PLANET
        rmins(9) = R771 / R_PLANET

        rmaxs(10) = R771 / R_PLANET
        rmins(10) = RTOPDDOUBLEPRIME / R_PLANET
      else
        ! 2nd doubling after 771
        rmaxs(5) = R220 / R_PLANET
        rmins(5) = R400 / R_PLANET

        rmaxs(6) = R400 / R_PLANET
        rmins(6) = R600 / R_PLANET

        rmaxs(7) = R600 / R_PLANET
        rmins(7) = R670 / R_PLANET

        rmaxs(8) = R670 / R_PLANET
        rmins(8) = R771 / R_PLANET

        rmaxs(9:10) = R771 / R_PLANET
        rmins(9:10) = RTOPDDOUBLEPRIME / R_PLANET
      endif

      rmaxs(11) = RTOPDDOUBLEPRIME / R_PLANET
      rmins(11) = RCMB / R_PLANET

      rmaxs(12:13) = RCMB / R_PLANET
      rmins(12:13) = RICB / R_PLANET

      rmaxs(14) = RICB / R_PLANET
      rmins(14) = R_CENTRAL_CUBE / R_PLANET

    endif
  else

    ! 4th doubling case:
    !     includes fourth doubling at the bottom of the outer core

    if (SUPPRESS_CRUSTAL_MESH) then

      ! suppress the crustal layers
      ! will be replaced by an extension of the mantle: R_PLANET is not modified,
      ! but no more crustal doubling

      ! check with define_all_layers_number_and_offset()
      ! NUMBER_OF_MESH_LAYERS = 15
      ! layer_offset = 1

      ! crust/mantle
      ! now only one region
      ner_mesh_layers( 1) = NER_CRUST + NER_80_MOHO
      ner_mesh_layers( 2) = 0
      ner_mesh_layers( 3) = 0
      ner_mesh_layers( 4) = NER_220_80

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ner_mesh_layers( 5) = NER_400_220 - elem_doubling_mantle
        ner_mesh_layers( 6) = elem_doubling_mantle
        ner_mesh_layers( 7) = NER_600_400
        ner_mesh_layers( 8) = NER_670_600
        ner_mesh_layers( 9) = NER_771_670
        ner_mesh_layers(10) = NER_TOPDDOUBLEPRIME_771
      else
        ! 2nd doubling after 771
        ner_mesh_layers( 5) = NER_400_220
        ner_mesh_layers( 6) = NER_600_400
        ner_mesh_layers( 7) = NER_670_600
        ner_mesh_layers( 8) = NER_771_670
        ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
        ner_mesh_layers(10) = elem_doubling_mantle
      endif
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME

      ! outer core
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner_mesh_layers(14) = elem_doubling_bottom_outer_core

      ! inner core
      ner_mesh_layers(15) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ratio_sampling_array(1:5) = 1
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(6:12) = 2
        else
          ratio_sampling_array(6:12) = ratio_sampling_array(9)
        endif
      else
        ! 2nd doubling after 771
        ratio_sampling_array(1:9) = 1
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(10:12) = 2
        else
          ratio_sampling_array(10:12) = ratio_sampling_array(9)
        endif
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
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_1ST_DOUBLING) then
          this_region_has_a_doubling(6) = .true.
          last_doubling_layer = 6
        endif
      else
        ! 2nd doubling after 771
        if (ADD_1ST_DOUBLING) then
          this_region_has_a_doubling(10) = .true.
          last_doubling_layer = 10
        endif
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

      r_top(1) = R_PLANET
      r_bottom(1) = R80_FICTITIOUS_IN_MESHER

      r_top(2) = RMIDDLE_CRUST    !!!! now fictitious
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER    !!!! now fictitious
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER    !!!! now fictitious

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        r_top(5) = R220
        r_bottom(5) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(6) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(6) = R400

        r_top(7) = R400
        r_bottom(7) = R600

        r_top(8) = R600
        r_bottom(8) = R670

        r_top(9) = R670
        r_bottom(9) = R771

        r_top(10) = R771
        r_bottom(10) = RTOPDDOUBLEPRIME
      else
        ! 2nd doubling after 771
        r_top(5) = R220
        r_bottom(5) = R400

        r_top(6) = R400
        r_bottom(6) = R600

        r_top(7) = R600
        r_bottom(7) = R670

        r_top(8) = R670
        r_bottom(8) = R771

        r_top(9) = R771
        r_bottom(9) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(10) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(10) = RTOPDDOUBLEPRIME
      endif

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = R_PLANET - DEPTH_FOURTH_DOUBLING_REAL

      r_top(14) = R_PLANET - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(14) = RICB

      r_top(15) = RICB
      r_bottom(15) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = R_UNIT_SPHERE
      rmins(1) = R80_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(2) = RMIDDLE_CRUST / R_PLANET    !!!! now fictitious
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET    !!!! now fictitious

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET    !!!! now fictitious
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET    !!!! now fictitious

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(4) = R220 / R_PLANET

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        rmaxs(5:6) = R220 / R_PLANET
        rmins(5:6) = R400 / R_PLANET

        rmaxs(7) = R400 / R_PLANET
        rmins(7) = R600 / R_PLANET

        rmaxs(8) = R600 / R_PLANET
        rmins(8) = R670 / R_PLANET

        rmaxs(9) = R670 / R_PLANET
        rmins(9) = R771 / R_PLANET

        rmaxs(10) = R771 / R_PLANET
        rmins(10) = RTOPDDOUBLEPRIME / R_PLANET
      else
        ! 2nd doubling after 771
        rmaxs(5) = R220 / R_PLANET
        rmins(5) = R400 / R_PLANET

        rmaxs(6) = R400 / R_PLANET
        rmins(6) = R600 / R_PLANET

        rmaxs(7) = R600 / R_PLANET
        rmins(7) = R670 / R_PLANET

        rmaxs(8) = R670 / R_PLANET
        rmins(8) = R771 / R_PLANET

        rmaxs(9:10) = R771 / R_PLANET
        rmins(9:10) = RTOPDDOUBLEPRIME / R_PLANET
      endif

      rmaxs(11) = RTOPDDOUBLEPRIME / R_PLANET
      rmins(11) = RCMB / R_PLANET

      rmaxs(12:14) = RCMB / R_PLANET
      rmins(12:14) = RICB / R_PLANET

      rmaxs(15) = RICB / R_PLANET
      rmins(15) = R_CENTRAL_CUBE / R_PLANET

    else if (ONE_CRUST) then

      ! 1D models:
      ! in order to increase CFL stability of the time scheme and therefore to allow cheaper
      ! simulations (larger time step), 1D models can be run with just one average crustal
      ! layer instead of two.

      ! check with define_all_layers_number_and_offset()
      ! NUMBER_OF_MESH_LAYERS = 14
      ! layer_offset = 0

      ! crust/mantle
      ner_mesh_layers( 1) = NER_CRUST
      ner_mesh_layers( 2) = NER_80_MOHO
      ner_mesh_layers( 3) = NER_220_80

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ner_mesh_layers( 4) = NER_400_220 - elem_doubling_mantle
        ner_mesh_layers( 5) = elem_doubling_mantle
        ner_mesh_layers( 6) = NER_600_400
        ner_mesh_layers( 7) = NER_670_600
        ner_mesh_layers( 8) = NER_771_670
        ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771
      else
        ! 2nd doubling after 771
        ner_mesh_layers( 4) = NER_400_220
        ner_mesh_layers( 5) = NER_600_400
        ner_mesh_layers( 6) = NER_670_600
        ner_mesh_layers( 7) = NER_771_670
        ner_mesh_layers( 8) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
        ner_mesh_layers( 9) = elem_doubling_mantle
      endif
      ner_mesh_layers(10) = NER_CMB_TOPDDOUBLEPRIME

      ! outer core
      ner_mesh_layers(11) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(12) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner_mesh_layers(13) = elem_doubling_bottom_outer_core

      ! inner core
      ner_mesh_layers(14) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1) = 1
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(2:4) = 2
        else
          ratio_sampling_array(2:4) = ratio_sampling_array(1)
        endif
        if (ADD_2ND_DOUBLING) then
          ratio_sampling_array(5:11) = 4 ! 2 * ratio_sampling_array(8)
        else
          ratio_sampling_array(5:11) = ratio_sampling_array(4)
        endif
      else
        ! 2nd doubling after 771
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
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(5)  = .true.
          last_doubling_layer = 5
        endif
      else
        ! 2nd doubling after 771
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(9)  = .true.
          last_doubling_layer = 9
        endif
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

      r_top(1) = R_PLANET
      r_bottom(1) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(2) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(2) = R80_FICTITIOUS_IN_MESHER

      r_top(3) = R80_FICTITIOUS_IN_MESHER
      r_bottom(3) = R220

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        r_top(4) = R220
        r_bottom(4) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(5) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(5) = R400

        r_top(6) = R400
        r_bottom(6) = R600

        r_top(7) = R600
        r_bottom(7) = R670

        r_top(8) = R670
        r_bottom(8) = R771

        r_top(9) = R771
        r_bottom(9) = RTOPDDOUBLEPRIME
      else
        ! 2nd doubling after 771
        r_top(4) = R220
        r_bottom(4) = R400

        r_top(5) = R400
        r_bottom(5) = R600

        r_top(6) = R600
        r_bottom(6) = R670

        r_top(7) = R670
        r_bottom(7) = R771

        r_top(8) = R771
        r_bottom(8) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(9) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(9) = RTOPDDOUBLEPRIME
      endif

      r_top(10) = RTOPDDOUBLEPRIME
      r_bottom(10) = RCMB

      r_top(11) = RCMB
      r_bottom(11) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL

      r_top(12) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(12) = R_PLANET - DEPTH_FOURTH_DOUBLING_REAL

      r_top(13) = R_PLANET - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(13) = RICB

      r_top(14) = RICB
      r_bottom(14) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = R_UNIT_SPHERE
      rmins(1) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(2) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(2) = R80_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(3) = R220 / R_PLANET

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        rmaxs(4:5) = R220 / R_PLANET
        rmins(4:5) = R400 / R_PLANET

        rmaxs(6) = R400 / R_PLANET
        rmins(6) = R600 / R_PLANET

        rmaxs(7) = R600 / R_PLANET
        rmins(7) = R670 / R_PLANET

        rmaxs(8) = R670 / R_PLANET
        rmins(8) = R771 / R_PLANET

        rmaxs(9) = R771 / R_PLANET
        rmins(9) = RTOPDDOUBLEPRIME / R_PLANET
      else
        ! 2nd doubling after 771
        rmaxs(4) = R220 / R_PLANET
        rmins(4) = R400 / R_PLANET

        rmaxs(5) = R400 / R_PLANET
        rmins(5) = R600 / R_PLANET

        rmaxs(6) = R600 / R_PLANET
        rmins(6) = R670 / R_PLANET

        rmaxs(7) = R670 / R_PLANET
        rmins(7) = R771 / R_PLANET

        rmaxs(8:9) = R771 / R_PLANET
        rmins(8:9) = RTOPDDOUBLEPRIME / R_PLANET
      endif

      rmaxs(10) = RTOPDDOUBLEPRIME / R_PLANET
      rmins(10) = RCMB / R_PLANET

      rmaxs(11:13) = RCMB / R_PLANET
      rmins(11:13) = RICB / R_PLANET

      rmaxs(14) = RICB / R_PLANET
      rmins(14) = R_CENTRAL_CUBE / R_PLANET

    else

      ! for 3D models:
      !   contains the crustal layers
      !   doubling at the base of the crust

      ! check with define_all_layers_number_and_offset()
      ! NUMBER_OF_MESH_LAYERS = 15
      ! layer_offset = 1

      ! crust/mantle
      if ((RMIDDLE_CRUST-RMOHO_FICTITIOUS_IN_MESHER) < (R_PLANET-RMIDDLE_CRUST)) then
        ner_mesh_layers( 1) = ceiling (NER_CRUST / 2.d0)
        ner_mesh_layers( 2) = floor (NER_CRUST / 2.d0)
      else
        ner_mesh_layers( 1) = floor (NER_CRUST / 2.d0)
        ner_mesh_layers( 2) = ceiling (NER_CRUST / 2.d0)
      endif
      ner_mesh_layers( 3) = NER_80_MOHO
      ner_mesh_layers( 4) = NER_220_80

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        ner_mesh_layers( 5) = NER_400_220 - elem_doubling_mantle
        ner_mesh_layers( 6) = elem_doubling_mantle
        ner_mesh_layers( 7) = NER_600_400
        ner_mesh_layers( 8) = NER_670_600
        ner_mesh_layers( 9) = NER_771_670
        ner_mesh_layers(10) = NER_TOPDDOUBLEPRIME_771
      else
        ! 2nd doubling after 771
        ner_mesh_layers( 5) = NER_400_220
        ner_mesh_layers( 6) = NER_600_400
        ner_mesh_layers( 7) = NER_670_600
        ner_mesh_layers( 8) = NER_771_670
        ner_mesh_layers( 9) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
        ner_mesh_layers(10) = elem_doubling_mantle
      endif
      ner_mesh_layers(11) = NER_CMB_TOPDDOUBLEPRIME

      ! outer core
      ner_mesh_layers(12) = NER_OUTER_CORE - elem_doubling_middle_outer_core
      ner_mesh_layers(13) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
      ner_mesh_layers(14) = elem_doubling_bottom_outer_core

      ! inner core
      ner_mesh_layers(15) = NER_TOP_CENTRAL_CUBE_ICB

      ! value of the doubling ratio in each radial region of the mesh
      ratio_sampling_array(1:2) = 1
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_1ST_DOUBLING) then
          ratio_sampling_array(3:5) = 2
        else
          ratio_sampling_array(3:5) = ratio_sampling_array(2)
        endif
        if (ADD_2ND_DOUBLING) then
          ratio_sampling_array(6:12) = 4 ! 2 * ratio_sampling_array(9)
        else
          ratio_sampling_array(6:12) = ratio_sampling_array(5)
        endif
      else
        ! 2nd doubling after 771
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
      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(6) = .true.
          last_doubling_layer = 6
        endif
      else
        ! 2nd doubling after 771
        if (ADD_2ND_DOUBLING) then
          this_region_has_a_doubling(10) = .true.
          last_doubling_layer = 10
        endif
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

      r_top(1) = R_PLANET
      r_bottom(1) = RMIDDLE_CRUST

      r_top(2) = RMIDDLE_CRUST
      r_bottom(2) = RMOHO_FICTITIOUS_IN_MESHER

      r_top(3) = RMOHO_FICTITIOUS_IN_MESHER
      r_bottom(3) = R80_FICTITIOUS_IN_MESHER

      r_top(4) = R80_FICTITIOUS_IN_MESHER
      r_bottom(4) = R220

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        r_top(5) = R220
        r_bottom(5) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(6) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(6) = R400

        r_top(7) = R400
        r_bottom(7) = R600

        r_top(8) = R600
        r_bottom(8) = R670

        r_top(9) = R670
        r_bottom(9) = R771

        r_top(10) = R771
        r_bottom(10) = RTOPDDOUBLEPRIME
      else
        ! 2nd doubling after 771
        r_top(5) = R220
        r_bottom(5) = R400

        r_top(6) = R400
        r_bottom(6) = R600

        r_top(7) = R600
        r_bottom(7) = R670

        r_top(8) = R670
        r_bottom(8) = R771

        r_top(9) = R771
        r_bottom(9) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL

        r_top(10) = R_PLANET - DEPTH_SECOND_DOUBLING_REAL
        r_bottom(10) = RTOPDDOUBLEPRIME
      endif

      r_top(11) = RTOPDDOUBLEPRIME
      r_bottom(11) = RCMB

      r_top(12) = RCMB
      r_bottom(12) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL

      r_top(13) = R_PLANET - DEPTH_THIRD_DOUBLING_REAL
      r_bottom(13) = R_PLANET - DEPTH_FOURTH_DOUBLING_REAL

      r_top(14) = R_PLANET - DEPTH_FOURTH_DOUBLING_REAL
      r_bottom(14) = RICB

      r_top(15) = RICB
      r_bottom(15) = R_CENTRAL_CUBE

      ! new definition of rmins & rmaxs
      rmaxs(1) = R_UNIT_SPHERE
      rmins(1) = RMIDDLE_CRUST / R_PLANET

      rmaxs(2) = RMIDDLE_CRUST / R_PLANET
      rmins(2) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(3) = RMOHO_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET

      rmaxs(4) = R80_FICTITIOUS_IN_MESHER / R_PLANET
      rmins(4) = R220 / R_PLANET

      if (REGIONAL_MESH_CUTOFF .and. ADD_2ND_DOUBLING) then
        ! 2nd doubling after 220
        rmaxs(5:6) = R220 / R_PLANET
        rmins(5:6) = R400 / R_PLANET

        rmaxs(7) = R400 / R_PLANET
        rmins(7) = R600 / R_PLANET

        rmaxs(8) = R600 / R_PLANET
        rmins(8) = R670 / R_PLANET

        rmaxs(9) = R670 / R_PLANET
        rmins(9) = R771 / R_PLANET

        rmaxs(10) = R771 / R_PLANET
        rmins(10) = RTOPDDOUBLEPRIME / R_PLANET
      else
        ! 2nd doubling after 771
        rmaxs(5) = R220 / R_PLANET
        rmins(5) = R400 / R_PLANET

        rmaxs(6) = R400 / R_PLANET
        rmins(6) = R600 / R_PLANET

        rmaxs(7) = R600 / R_PLANET
        rmins(7) = R670 / R_PLANET

        rmaxs(8) = R670 / R_PLANET
        rmins(8) = R771 / R_PLANET

        rmaxs(9:10) = R771 / R_PLANET
        rmins(9:10) = RTOPDDOUBLEPRIME / R_PLANET
      endif

      rmaxs(11) = RTOPDDOUBLEPRIME / R_PLANET
      rmins(11) = RCMB / R_PLANET

      rmaxs(12:14) = RCMB / R_PLANET
      rmins(12:14) = RICB / R_PLANET

      rmaxs(15) = RICB / R_PLANET
      rmins(15) = R_CENTRAL_CUBE / R_PLANET
    endif
  endif

  ! simpler local mesh feature
  if (REGIONAL_MESH_CUTOFF .and. USE_LOCAL_MESH) then
    ! re-defines a local mesh (if selected) with more number of doubling layers
    call define_all_layers_for_local_mesh(NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer,rmins,rmaxs)
  endif

  ! debug
  if (DEBUG .and. myrank == 0) then
    print *,'debug: define_all_layers:',NUMBER_OF_MESH_LAYERS
    do ielem = 1,NUMBER_OF_MESH_LAYERS
      print *,'debug:  layer ',ielem,': top/bottom ',sngl(r_top(ielem)),sngl(r_bottom(ielem)), &
              'rmin/rmax = ',sngl(rmins(ielem)),sngl(rmaxs(ielem)),'ner',ner_mesh_layers(ielem), &
              'doubling',this_region_has_a_doubling(ielem),ratio_sampling_array(ielem)
    enddo
  endif

  ! checks arrays
  if (NUMBER_OF_MESH_LAYERS <= 0) &
    stop 'Error invalid number of mesh layers in define_all_layers()'
  if (minval(ratio_sampling_array(1:NUMBER_OF_MESH_LAYERS)) <= 0) &
    stop 'Error invalid ratio array in define_all_layers()'
  if (minval(doubling_index(1:NUMBER_OF_MESH_LAYERS)) <= 0) &
    stop 'Error invalid doubling array in define_all_layers()'
  if (minval(r_top(1:NUMBER_OF_MESH_LAYERS)) <= 0.d0) &
    stop 'Error invalid r_top array in define_all_layers()'
  if (minval(r_bottom(1:NUMBER_OF_MESH_LAYERS)) <= 0.d0) &
    stop 'Error invalid r_bottom array in define_all_layers()'
  if (minval(rmins(1:NUMBER_OF_MESH_LAYERS)) <= 0.d0) &
    stop 'Error invalid rmins array in define_all_layers()'
  if (minval(rmaxs(1:NUMBER_OF_MESH_LAYERS)) <= 0.d0) &
    stop 'Error invalid rmaxs array in define_all_layers()'

  ! checks layering
  do ielem = 1,NUMBER_OF_MESH_LAYERS-1
    if (r_top(ielem) - r_top(ielem+1) < 0.d0) then
      print *,'Error: define_all_layers rank ',myrank,'found invalid layer ',ielem, &
              ' with negative radius separation: r_top ',r_top(ielem),r_top(ielem+1)
      stop 'Error invalid layering in define all layers'
    endif
    if (r_bottom(ielem) - r_bottom(ielem+1) < 0.d0) then
      print *,'Error: define_all_layers rank ',myrank,'found invalid layer ',ielem, &
              ' with negative radius separation: r_bottom ',r_bottom(ielem),r_bottom(ielem+1)
      stop 'Error invalid layering in define all layers'
    endif
  enddo

  end subroutine define_all_layers


!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_all_layers_number_and_offset(NUMBER_OF_MESH_LAYERS,layer_offset)

  use constants, only: MAX_NUMBER_OF_MESH_LAYERS,ADD_4TH_DOUBLING
  use shared_parameters, only: ONE_CRUST

  implicit none

  integer, intent(inout) :: NUMBER_OF_MESH_LAYERS,layer_offset

  ! sets total number of layers
  !   NUMBER_OF_MESH_LAYERS  -  total number of mesh layers
  !   layer_offset           -  crust/mantle region has 10 + layer_offset layers
  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
    layer_offset = 0
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
    layer_offset = 1
  endif

  if (.not. ADD_4TH_DOUBLING) NUMBER_OF_MESH_LAYERS = NUMBER_OF_MESH_LAYERS - 1

  end subroutine define_all_layers_number_and_offset

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_all_layers_for_local_mesh(NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer,rmins,rmaxs)

  use constants, only: myrank,MAX_NUMBER_OF_MESH_LAYERS,R_UNIT_SPHERE, &
    IFLAG_CRUST,IFLAG_MANTLE_NORMAL,IFLAG_OUTER_CORE_NORMAL,IFLAG_INNER_CORE_NORMAL

  use shared_parameters, only: R_PLANET,RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER, &
    RTOPDDOUBLEPRIME,RCMB,RICB,R_CENTRAL_CUBE

  use shared_parameters, only: ner_mesh_layers, &
    ratio_sampling_array,this_region_has_a_doubling,doubling_index,r_bottom,r_top

  use shared_parameters, only: REGIONAL_MESH_CUTOFF,REGIONAL_MESH_CUTOFF_DEPTH, &
    NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  use shared_parameters, only: USE_LOCAL_MESH,LOCAL_MESH_NUMBER_OF_LAYERS_CRUST,LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE, &
    NDOUBLINGS,NZ_DOUBLING_1,NZ_DOUBLING_2,NZ_DOUBLING_3,NZ_DOUBLING_4,NZ_DOUBLING_5

  implicit none

  integer,intent(inout) :: NUMBER_OF_MESH_LAYERS,layer_offset,last_doubling_layer
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS),intent(inout) :: rmins,rmaxs

  ! local parameters
  ! mesh layering
  integer :: ilayer,ilayer_top,ilayer_bottom,ilayer_last,idoubling,ilocal
  integer :: LOCAL_MESH_NUMBER_OF_LAYERS
  ! element layer numbers of doubling zone
  integer, dimension(5) :: ner_doublings

  integer :: num_layers_top,N
  double precision :: layer_thickness,rmesh_bottom,rmesh_moho
  double precision :: ratio,th_layer,F
  logical :: has_doubling

  !---------------------------------------
  ! MESHING PARAMETERS

  ! scaling factor of vertical layer width with depth (with 1.0 == no stretching)
  double precision, parameter :: LAYER_SCALING_FACTOR = 1.8d0

  ! debugging
  logical, parameter :: DEBUG = .false.

  !---------------------------------------

  ! checks if anything to do
  if (.not. REGIONAL_MESH_CUTOFF) return
  if (.not. USE_LOCAL_MESH) return

  ! initializes
  ! fills element layer numbers
  ner_doublings(:) = 0
  ner_doublings(1) = NZ_DOUBLING_1
  ner_doublings(2) = NZ_DOUBLING_2
  ner_doublings(3) = NZ_DOUBLING_3
  ner_doublings(4) = NZ_DOUBLING_4
  ner_doublings(5) = NZ_DOUBLING_5

  ! user info
  if (myrank == 0) then
    print *,'using local mesh layout:'
    print *,'  number of layers in crust  = ',LOCAL_MESH_NUMBER_OF_LAYERS_CRUST
    print *,'  number of layers in mantle = ',LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE
    print *
    print *,'  number of doubling layers  = ',NDOUBLINGS
    do idoubling = 1,NDOUBLINGS
      print *,'           doubling layer at = ',ner_doublings(idoubling)
    enddo
    print *,'  fictitious moho at depth   : ',sngl((R_PLANET - RMOHO_FICTITIOUS_IN_MESHER)/1000.d0),'(km)'
    print *
  endif

  !debug
  if (DEBUG) then
    print *,'debug: local mesh - NUMBER_OF_MESH_LAYERS = ',NUMBER_OF_MESH_LAYERS,' / layer_offset = ',layer_offset
    print *,'debug: local mesh - ner ',ner_mesh_layers(:)
  endif

  ! default mesh layering:
  ! NUMBER_OF_MESH_LAYERS - total number of ner_mesh_layers()
  !                         crust/mantle: 1 to 10 + layer_offset, i.e., ner_mesh_layers(1 : 10+layer_offset)
  !                         outer core: entries between, i.e., ner_mesh_layers(10+layer_offset + 1 : NUMBER_OF_MESH_LAYERS-1)
  !                         inner core: last entry, i.e., ner_mesh_layers(NUMBER_OF_MESH_LAYERS)
  !
  ! for REGIONAL_MESH_CUTOFF, we assigned ner values of 0 to layers ner_mesh_layers(i) below the cut-off, say 400.
  ! here, we will re-assign ner_mesh_layers(i) and rmins/rmax values to create a new local mesh, with a cut-off below the desired depth.
  !
  ! NDOUBLINGS determines how many different ner_mesh_layer(i) entries we will need.
  ! REGIONAL_MESH_CUTOFF_DEPTH tells where to set rmins/rmax

  ! re-initializes for local mesh
  ner_mesh_layers(:) = 0      ! number of element layers

  this_region_has_a_doubling(:) = .false.
  last_doubling_layer = 0

  ratio_sampling_array(:) = 1 ! doubling ratio

  r_top(:) = 0.d0
  r_bottom(:) = 0.d0

  rmins(:) = 0.d0
  rmaxs(:) = 0.d0

  ! global setup
  ! inner core mesh
  ! assigns last layer to inner core
  doubling_index(NUMBER_OF_MESH_LAYERS) = IFLAG_INNER_CORE_NORMAL

  r_top(NUMBER_OF_MESH_LAYERS) = RICB
  r_bottom(NUMBER_OF_MESH_LAYERS) = R_CENTRAL_CUBE

  rmaxs(NUMBER_OF_MESH_LAYERS) = RICB / R_PLANET
  rmins(NUMBER_OF_MESH_LAYERS) = R_CENTRAL_CUBE / R_PLANET

  ! outer core mesh
  ! next layer(s) to outer core
  ilayer_top = 10 + layer_offset + 1
  ilayer_bottom = NUMBER_OF_MESH_LAYERS - 1

  doubling_index(ilayer_top:ilayer_bottom) = IFLAG_OUTER_CORE_NORMAL

  layer_thickness = (RICB - RCMB) / (ilayer_bottom - ilayer_top + 1)
  do ilayer = ilayer_top,ilayer_bottom
    r_top(ilayer) = RCMB + (ilayer - ilayer_top) * layer_thickness
    r_bottom(ilayer) = RCMB + (ilayer - ilayer_top + 1) * layer_thickness
  enddo
  ! makes sure bottom has RICB limit, in case above incremental contribution has some numerical round-off error
  r_bottom(ilayer_bottom) = RICB

  rmaxs(ilayer_top:ilayer_bottom) = RCMB / R_PLANET
  rmins(ilayer_top:ilayer_bottom) = RICB / R_PLANET

  ! crust/mantle mesh
  ! we'll first use a bottom layer of mantle from CMB to D'' with a zero ner entry
  ! to make sure that counting elements/points later will work (might not be needed anymore though, just to be on the safe side)
  doubling_index(10+layer_offset) = IFLAG_MANTLE_NORMAL

  r_top(10+layer_offset) = RTOPDDOUBLEPRIME
  r_bottom(10+layer_offset) = RCMB

  rmaxs(10+layer_offset) = RTOPDDOUBLEPRIME / R_PLANET
  rmins(10+layer_offset) = RCMB / R_PLANET

  ! mantle and crust zones
  ! we have a total of (LOCAL_MESH_NUMBER_OF_LAYERS_CRUST + LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE) element layers
  ! in case there is a doubling in such a zone, we need to create a new ner_mesh_layer(i) entry for it.

  ! at least, we will need (NDOUBLINGS + 1) entries in ner_mesh_layers to have different ner_mesh_layers(i) entries for each zone;
  !! we will also need +1 entry to separate crust & mantle zones.
  ! also, we can have a maximum of (10+layer_offset) - 1 for these crust/mantle entries.
  ! this limits the NDOUBLINGS value
  ilayer_last = NDOUBLINGS + 1
  if (ilayer_last > 10 + layer_offset - 1) stop 'Invalid NDOUBLINGS values, must be less than 10+layer_offset-1'

  ! sets remaining layers between mesh cut-off depth and last layer in mantle which goes to till top of D''
  ilayer_top = ilayer_last + 1
  ilayer_bottom = 10 + layer_offset - 1

  doubling_index(ilayer_top:ilayer_bottom) = IFLAG_MANTLE_NORMAL  ! will assign normal mantle flag to remaining layers

  ! radius of bottom depth (in m)
  rmesh_bottom = R_PLANET - REGIONAL_MESH_CUTOFF_DEPTH * 1000.d0
  rmesh_moho = RMOHO_FICTITIOUS_IN_MESHER

  ! debug
  if (DEBUG) then
    print *,'debug: local mesh - depth  moho (km) = ',sngl((R_PLANET-rmesh_moho)/1000.d0), &
                                    ' bottom (km) = ',sngl(REGIONAL_MESH_CUTOFF_DEPTH/1000.d0)
    print *,'debug: local mesh - radius moho (km) = ',sngl(rmesh_moho/1000.d0), &
                                    ' bottom (km) = ',sngl(rmesh_bottom/1000.d0)
  endif

  ! sets top/bottom radius of layers
  layer_thickness = (RTOPDDOUBLEPRIME - rmesh_bottom) / (ilayer_bottom - ilayer_top + 1)
  do ilayer = ilayer_top,ilayer_bottom
    r_top(ilayer) = rmesh_bottom + (ilayer - ilayer_top) * layer_thickness
    r_bottom(ilayer) = rmesh_bottom + (ilayer - ilayer_top + 1) * layer_thickness
  enddo
  ! makes sure last layer is on D''
  r_bottom(ilayer_bottom) = RTOPDDOUBLEPRIME

  rmaxs(ilayer_top:ilayer_bottom) = rmesh_bottom / R_PLANET
  rmins(ilayer_top:ilayer_bottom) = RTOPDDOUBLEPRIME / R_PLANET

  ! assigns layering for crust/mantle
  if (NDOUBLINGS == 0) then
    ! no doubling layers
    ! entry in ner_mesh_layers for crust
    ner_mesh_layers(1) = LOCAL_MESH_NUMBER_OF_LAYERS_CRUST
    doubling_index(1) = IFLAG_CRUST  ! will assign crust flag
    r_top(1) = R_PLANET
    r_bottom(1) = rmesh_moho
    rmaxs(1) = R_UNIT_SPHERE
    rmins(1) = rmesh_moho / R_PLANET

    ! entry in ner_mesh_layers for mantle
    if (LAYER_SCALING_FACTOR <= 1.d0 .or. rmesh_bottom >= R80_FICTITIOUS_IN_MESHER) then
      ! regular layer thickness
      ner_mesh_layers(2) = LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE
      doubling_index(2) = IFLAG_CRUST  ! will assign crust flag (to include topography stretching, in case, down to 220km)
      r_top(2) = rmesh_moho
      r_bottom(2) = rmesh_bottom
      rmaxs(2) = rmesh_moho / R_PLANET
      rmins(2) = rmesh_bottom / R_PLANET
    else
      ! mesh zone between moho - R80
      ! determines number of element layers between moho - R80
      ratio = (rmesh_moho - R80_FICTITIOUS_IN_MESHER) / (rmesh_moho - rmesh_bottom)
      num_layers_top = ceiling(LAYER_SCALING_FACTOR * ratio * LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE)

      ! puts at least 1 layer between moho
      if (num_layers_top <= 0) num_layers_top = 1
      ! in case rmesh_bottom is close to R80, let's at least 1 element layer for bottom zone below R80
      if (num_layers_top >= LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE) num_layers_top = LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE - 1

      ! debugging
      if (DEBUG) then
        print *,'debug: num_layers_top ',num_layers_top,LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE,' - R80',sngl(R80_FICTITIOUS_IN_MESHER)
      endif

      if (num_layers_top == 0) then
        ! only single zone below moho
        ner_mesh_layers(2) = LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE
        doubling_index(2) = IFLAG_CRUST  ! will assign crust flag (to include topography stretching, in case, down to 220km)
        r_top(2) = rmesh_moho
        r_bottom(2) = rmesh_bottom
        rmaxs(2) = rmesh_moho / R_PLANET
        rmins(2) = rmesh_bottom / R_PLANET
      else
        ! splits zones at R80, with a refined layering for top moho - R80 zone
        ner_mesh_layers(2) = num_layers_top
        doubling_index(2) = IFLAG_CRUST
        r_top(2) = rmesh_moho
        r_bottom(2) = R80_FICTITIOUS_IN_MESHER
        rmaxs(2) = rmesh_moho / R_PLANET
        rmins(2) = R80_FICTITIOUS_IN_MESHER / R_PLANET

        ! mesh zone between R80 - mesh bottom
        ner_mesh_layers(3) = LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE - num_layers_top
        doubling_index(3) = IFLAG_CRUST
        r_top(3) = R80_FICTITIOUS_IN_MESHER
        r_bottom(3) = rmesh_bottom
        rmaxs(3) = R80_FICTITIOUS_IN_MESHER / R_PLANET
        rmins(3) = rmesh_bottom / R_PLANET
      endif
    endif

  else
    ! with doubling layers
    ! checks
    if (any(ner_doublings(:) == 1)) stop 'Cannot have doubling in first element layer'

    ! initial layer
    ilayer = 1
    r_top(1) = R_PLANET
    r_bottom(1) = R_PLANET ! initial bottom radius

    rmaxs(1) = R_UNIT_SPHERE
    rmins(1) = R_UNIT_SPHERE  ! initial minimum

    doubling_index(1) = IFLAG_CRUST

    ! total number of mesh layers
    LOCAL_MESH_NUMBER_OF_LAYERS = LOCAL_MESH_NUMBER_OF_LAYERS_CRUST + LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE

    ! assigns number of element layers, top down
    ilayer = 1
    do ilocal = 1,LOCAL_MESH_NUMBER_OF_LAYERS
      ! increases layer if layer has a doubling
      has_doubling = .false.
      do idoubling = 1,NDOUBLINGS
        ! layer number
        if (ner_doublings(idoubling) == ilocal) then
          has_doubling = .true.
          exit
        endif
      enddo

      if (has_doubling .or. ilocal == LOCAL_MESH_NUMBER_OF_LAYERS_CRUST+1) then
        ! increase ner_mesh_layer
        ilayer = ilayer + 1

        ! sets doubling flag
        if (has_doubling) then
          this_region_has_a_doubling(ilayer) = .true.

          ! increases sampling ratio for doubling (and all following layers below will have same new ratio)
          ratio_sampling_array(ilayer:MAX_NUMBER_OF_MESH_LAYERS) = 2 * ratio_sampling_array(ilayer)

          ! checks ratio: result of NEX / ratio must have at least a value of 2, otherwise element count is off
          if (NEX_PER_PROC_XI / ratio_sampling_array(ilayer) < 2 .or. NEX_PER_PROC_ETA / ratio_sampling_array(ilayer) < 2) then
            print *,'Error invalid ratio_sampling_array value: layer ',ilayer
            print *,'  NEX_PER_PROC XI/ETA  = ',NEX_PER_PROC_XI,NEX_PER_PROC_ETA
            print *,'  ratio_sampling_array = ',ratio_sampling_array(ilayer)
            print *,'Please increase NEX_XI / NEX_ETA value in Par_file, or decrease the number of doubling layers'
            stop 'Invalid number of doubling layers for NEX'
          endif

          ! sets last doubling layer index
          last_doubling_layer = ilocal
        endif
        doubling_index(ilayer) = IFLAG_CRUST  ! will assign crust flag (to include stretching for topography)

        ! top/bottom layer
        r_top(ilayer) = r_bottom(ilayer - 1) ! from previous layer
        r_bottom(ilayer) = r_top(ilayer)     ! initial bottom radius

        ! rmin/rmax
        rmaxs(ilayer) = r_top(ilayer) / R_PLANET
        rmins(ilayer) = rmaxs(ilayer)        ! initial max radius
      endif

      ! increase count of layer in this ner_mesh_layer
      ner_mesh_layers(ilayer) = ner_mesh_layers(ilayer) + 1

      ! layering, thickness per layer (note that thickness has positive value)
      if (ilocal <= LOCAL_MESH_NUMBER_OF_LAYERS_CRUST) then
        ! crust layers
        layer_thickness = (R_PLANET - rmesh_moho) / LOCAL_MESH_NUMBER_OF_LAYERS_CRUST
      else
        ! mantle layers
        if (LAYER_SCALING_FACTOR <= 1.d0 .or. REGIONAL_MESH_CUTOFF_DEPTH <= 80.d0) then
          ! regular layer thickness
          layer_thickness = (rmesh_moho - rmesh_bottom) / LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE
        else
          ! progressive layering, with refined layers from moho down
          ! factors
          N = LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE
          if (N > 1) then
            ! to determine total thickness:
            ! th_layer + (th_layer + (F-1) * th_layer * 1/(N-1)) + (th_layer + (F-1) * th_layer * 2/(N-1)) + ..
            ! = N * th_layer + (F-1) * th_layer * 1/(N-1) + (F-1) * th_layer * 2/(N-1) + .. + F * th_layer * (N-1)/(N-1)
            ! = N * th_layer + (F-1) * th_layer * (1 + 2 + .. + (N-1))/(N-1)
            ! = th_layer * ( N + (F-1) * (1 + 2 + .. (N-1))/(N-1) )
            ! = th_layer * ( N + (F-1) * (N-1) * N / 2 / (N-1) )
            ! should be equal to (Rmoho - Rbottom)
            F = LAYER_SCALING_FACTOR
            th_layer = (rmesh_moho - rmesh_bottom) / dble(N + (F - 1.d0) * (N-1) * N / 2 / (N-1))
            ratio = dble(ilocal - LOCAL_MESH_NUMBER_OF_LAYERS_CRUST - 1) / (N - 1)

            layer_thickness = th_layer + ratio * (F - 1.d0) * th_layer

            ! debugging
            if (DEBUG) then
              print *,'debug: layer thickness ',layer_thickness,' - th_layer',th_layer,ilocal,ratio,F,N
            endif
          else
            ! regular layer thickness (single element layer for mantle)
            layer_thickness = (rmesh_moho - rmesh_bottom) / LOCAL_MESH_NUMBER_OF_LAYERS_MANTLE
          endif
        endif
      endif

      ! moves bottom layer down
      r_bottom(ilayer) = r_bottom(ilayer) - layer_thickness

      ! minimum radius (non-dimensionalized)
      rmins(ilayer) = r_bottom(ilayer) / R_PLANET

      ! debug
      if (DEBUG) then
        print *,'debug: local_mesh ilocal ',ilocal, &
                ' - ilayer',ilayer,'thickness',sngl(layer_thickness), &
                'top/bottom',sngl(r_top(ilayer)),sngl(r_bottom(ilayer)),sngl(rmesh_bottom)
      endif
    enddo
  endif

  ! alternative: uses same cutoff mesh as provided with one additional doubling between upper/lower crust
  ! mesh doubling layers in crust
  ! crust - moho - 80 - 220
  ! default
  !ner_mesh_layers( 1) = ner_mesh_layers( 1) !+ 1      ! upper crust
  !ner_mesh_layers( 2) = ner_mesh_layers( 2) !+ 1      ! lower crust

  ! value of the doubling ratio in each radial region of the mesh
  !ratio_sampling_array( 1) = 1
  !ratio_sampling_array( 2) = 2    ! doubling at lower crust
  !do ilayer = 3,NUMBER_OF_MESH_LAYERS
  !  ratio_sampling_array(ilayer) = 2 * ratio_sampling_array(ilayer)
  !enddo

  ! define the three regions in which we implement a mesh doubling at the top of that region
  !this_region_has_a_doubling(2)  = .true.
  !if (last_doubling_layer < 1) last_doubling_layer = 2

  !debug
  if (DEBUG) then
    print *
    do ilayer = 1,NUMBER_OF_MESH_LAYERS
      print *,'debug: local_mesh ilayer ',ilayer,'out of ',NUMBER_OF_MESH_LAYERS
      print *,'debug:   ner / this_region_has_a_doubling / ratio_sampling_array ', &
              ner_mesh_layers(ilayer),this_region_has_a_doubling(ilayer),ratio_sampling_array(ilayer)
      print *,'debug:   r_top/r_bottom                                          ',r_top(ilayer),r_bottom(ilayer)
      print *,'debug:   rmins/rmaxs                                             ',rmins(ilayer),rmaxs(ilayer)
    enddo
    print *
  endif

  end subroutine define_all_layers_for_local_mesh

