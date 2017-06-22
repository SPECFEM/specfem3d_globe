!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 7
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!          (c) California Institute of Technology January 2007
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

  program display_prem_sampling_doubling

  implicit none

  include "constants_modified_2D.h"

! honor PREM Moho or not
! doing so drastically reduces the stability condition and therefore the time step
  logical, parameter :: HONOR_1D_SPHERICAL_MOHO = .false.

! resolution target for minimum number of points per S wavelength in the whole mesh
  double precision, parameter :: RESOLUTION_TARGET = 4.d0

! for outer core, in which there is no S wave velocity, use P velocity potential with a margin to be safe
  double precision, parameter :: MARGIN_SAFE = 1.25d0

  integer :: NEX,NEX_XI

  integer :: NER_220_CRUST, NER_400_220, NER_600_400, NER_670_600, NER_771_670, &
    NER_TOPDDOUBLEPRIME_771, NER_CMB_TOPDDOUBLEPRIME, NER_OUTER_CORE, NER_TOP_CENTRAL_CUBE_ICB

  double precision :: DT

  integer, dimension(NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling
  integer number_basic_elems_horizontal,ner_without_doubling

! small offset to avoid being exactly on an interface
  double precision, parameter :: SMALL_OFFSET = 200.d0

  double precision, parameter :: R_CENTRAL_CUBE = RICB - 150000.d0

  integer ilayer,ipoin,ielem

  double precision :: r,r_next,rho,vp,vs,ratio_sampling,spectral_element_size, &
               average_mesh_point_size,num_points_per_lambda_S,stability,period_min,gamma, &
               ratio_sampling_new,vp_new,vs_new,r_gnuplot

!!!!!! DK DK exclude center of inner core to avoid division by zero
!!!!!! DK DK we have: 1000 * 1200/6371 = 188.35347
!!!!!! DK DK therefore the central cube ends around i = 188
!!!!!! DK DK  do i=1,1000
  integer, parameter :: NINIT = 140

! mesh chunk has a size of 90 x 90 degrees
  double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES = 90.d0
  double precision, parameter :: ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * PI / 180.d0

  integer ix,ir
  integer ix_elem,ir_elem

! number of nodes of an OpenDX element
  integer, parameter :: NGNOD_OPENDX = 4

! topology of the elements
  integer, dimension(NGNOD_OPENDX) :: iaddx,iaddr

  integer :: ispec,ispec_final,nspec,npoin,ioffset,ignod,ignod2
  double precision, dimension(NGNOD_OPENDX) :: xelm,zelm
  double precision :: xval,zval,rval

! list of corners defining the edges
  integer, parameter :: NEDGES = 4
  integer, dimension(NEDGES,2) :: list_corners_edge
  integer iedge

  integer :: elem_doubling_mantle,elem_doubling_middle_outer_core,elem_doubling_bottom_outer_core
  double precision :: DEPTH_SECOND_DOUBLING_REAL,DEPTH_THIRD_DOUBLING_REAL, &
                          DEPTH_FOURTH_DOUBLING_REAL,distance,distance_min,distance_max

! mesh doubling superbrick
  integer, dimension(NGNOD_DOUBLING_SUPERBRICK,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,z_superbrick

! for the stability condition
! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLLX_MAX_STABILITY = 15
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

! allocate arrays to store Courant stability value and mesh resolution for S waves
  integer, parameter :: NSPEC_WORST_ELEMENTS = 50
  integer, dimension(NSPEC_WORST_ELEMENTS) :: ispec_worst_elements_stability,ispec_worst_elements_sampling
  integer, dimension(:), allocatable :: color
  double precision, dimension(:), allocatable :: courant_stability_number, &
      copy_courant_stability_number,number_points_S_wavelength,copy_number_points_S_wavelength

! number of different cases to study
  integer, parameter :: NUM_CASES = 9
  integer, dimension(NUM_CASES) :: NEX_val
  integer icase

! variable file name for output files
  character(len=150) filename

! remove old files
  call system('rm -f DX_fullmesh*.dx DX_worst*.dx prem_horizontal_sampling_*.dat prem_radial_sampling*.dat')

! case of multiples of 16
  NEX_val(1) = 160
  NEX_val(2) = 128
  NEX_val(3) = 320
  NEX_val(4) = 480
  NEX_val(5) = 512
  NEX_val(6) = 640
  NEX_val(7) = 864
  NEX_val(8) = 1152
  NEX_val(9) = 1248

! loop on all the cases to study
!!!!!!!!!!!! DK DK only one case    do icase = 1,NUM_CASES
  do icase = 2,2

! define value of NEX for this case
  NEX = NEX_val(icase)

  NEX_XI = NEX

! check that mesh can be coarsened in depth three times (block size must be a multiple of 16)
  if (mod(NEX_XI,16) /= 0) stop 'NEX_XI must be a multiple of 16'

  print *,'-----------------------------------------'
  print *
  print *,'NEX used = ',NEX
  print *

! period min for    4.000000 points per lambda S min in horizontal direction =    17.11689
! element width =   0.3515625      degrees =    39.09196      km
  if (NEX == 128) then
    DT                       = 0.20d0 * 0.30d0 / 0.2565
! stability max = approximately 0.2565
    NER_220_CRUST            = 4
    NER_400_220              = 3
    NER_600_400              = 3
    NER_670_600              = 1
    NER_771_670              = 1
    NER_TOPDDOUBLEPRIME_771  = 22
    NER_CMB_TOPDDOUBLEPRIME  = 2
    NER_OUTER_CORE           = 24
    NER_TOP_CENTRAL_CUBE_ICB = 2

! period min for 4.000000 points per lambda S min in horizontal direction =    13.69351
! element width =   0.2812500      degrees =    31.27357      km
  else if (NEX == 320) then
    DT                       = 0.125d0 * 0.30d0 / 0.2429
! stability max = approximately 0.2429
    NER_220_CRUST            = 2
    NER_400_220              = 4
    NER_600_400              = 4
    NER_670_600              = 1
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 29
    NER_CMB_TOPDDOUBLEPRIME  = 2
    NER_OUTER_CORE           = 32
    NER_TOP_CENTRAL_CUBE_ICB = 2

! period min for 4.000000 points per lambda S min in horizontal direction =    9.129005
! element width =   0.1875000      degrees =    20.84905      km
  else if (NEX == 480) then
    DT                       = 0.125d0 * 0.30d0 / 0.3235
! stability max = approximately 0.3235
    NER_220_CRUST            = 3
    NER_400_220              = 5
    NER_600_400              = 6
    NER_670_600              = 2
    NER_771_670              = 2
    NER_TOPDDOUBLEPRIME_771  = 44
    NER_CMB_TOPDDOUBLEPRIME  = 3
    NER_OUTER_CORE           = 48
    NER_TOP_CENTRAL_CUBE_ICB = 3

! period min for 4.000000 points per lambda S min in horizontal direction =    8.558443
! element width =   0.1757812      degrees =    19.54598      km
  else if (NEX == 512) then
    DT                       = 0.125d0 * 0.30d0 / 0.3235
! stability max = approximately 0.3235
    NER_220_CRUST            = 4
    NER_400_220              = 6
    NER_600_400              = 6
    NER_670_600              = 2
    NER_771_670              = 3
    NER_TOPDDOUBLEPRIME_771  = 47
    NER_CMB_TOPDDOUBLEPRIME  = 3
    NER_OUTER_CORE           = 51
    NER_TOP_CENTRAL_CUBE_ICB = 3

! period min for 4.000000 points per lambda S min in horizontal direction =    6.846754
! element width =   0.1406250      degrees =    15.63679      km
  else if (NEX == 640) then
    DT                       = 0.125d0 * 0.30d0 / 0.4067
! stability max = approximately 0.4067
    NER_220_CRUST            = 4
    NER_400_220              = 7
    NER_600_400              = 8
    NER_670_600              = 3
    NER_771_670              = 3
    NER_TOPDDOUBLEPRIME_771  = 59
    NER_CMB_TOPDDOUBLEPRIME  = 4
    NER_OUTER_CORE           = 64
    NER_TOP_CENTRAL_CUBE_ICB = 4

! period min for 4.000000 points per lambda S min in horizontal direction =    5.071670
! element width =   0.1041667      degrees =    11.58280      km
  else if (NEX == 864) then
    DT                       = 0.0555555555d0 * 0.30d0 / 0.2565
! stability max = approximately 0.2565
    NER_220_CRUST            = 6
    NER_400_220              = 10
    NER_600_400              = 10
    NER_670_600              = 3
    NER_771_670              = 4
    NER_TOPDDOUBLEPRIME_771  = 79
    NER_CMB_TOPDDOUBLEPRIME  = 5
    NER_OUTER_CORE           = 86
    NER_TOP_CENTRAL_CUBE_ICB = 6

! period min for 4.000000 points per lambda S min in horizontal direction =    3.803752
! element width =   7.8125000E-02  degrees =    8.687103      km
  else if (NEX == 1152) then
    DT                       = 0.0555555555d0 * 0.30d0 / 0.3504
! stability max = approximately 0.3504
    NER_220_CRUST            = 8
    NER_400_220              = 13
    NER_600_400              = 13
    NER_670_600              = 4
    NER_771_670              = 6
    NER_TOPDDOUBLEPRIME_771  = 106
    NER_CMB_TOPDDOUBLEPRIME  = 7
    NER_OUTER_CORE           = 116
    NER_TOP_CENTRAL_CUBE_ICB = 6

! period min for 4.000000 points per lambda S min in horizontal direction =    3.511156
! element width =   7.2115384E-02  degrees =    8.018865      km
  else if (NEX == 1248) then
    DT                       = 0.05d0 * 0.30d0 / 0.3318
! stability max = approximately 0.3318
    NER_220_CRUST            = 8
    NER_400_220              = 14
    NER_600_400              = 14
    NER_670_600              = 5
    NER_771_670              = 6
    NER_TOPDDOUBLEPRIME_771  = 114
    NER_CMB_TOPDDOUBLEPRIME  = 8
    NER_OUTER_CORE           = 124
    NER_TOP_CENTRAL_CUBE_ICB = 7

  else
    stop 'incorrect value of NEX, should use an updated version of auto_NER'
  endif

!! DK DK UGLY
  print *,'DT computed for Courant number of 0.30 = ',DT
  print *

!
!--- display the PREM model alone
!
  open(unit=27,file='prem_model.dat',status='unknown')
  do ipoin = 2,1000
    r = R_EARTH * dble(ipoin)/1000.d0
    call prem_iso(r,rho,vp,vs,ratio_sampling)
    write(27,*) sngl((R_EARTH-r)/1000.d0),sngl(vp),sngl(vs)
  enddo
  close(27)

!---------------------------------------------------
!--- mesh resolution in the horizontal direction
!---------------------------------------------------

! estimate of mesh resolution: minimum resolution is for S waves right below the second doubling
  r = R_EARTH - DEPTH_SECOND_DOUBLING_OPTIMAL - SMALL_OFFSET ! in meters
  call prem_iso(r,rho,vp,vs,ratio_sampling)
  period_min = RESOLUTION_TARGET / (vs * (NGLLX - 1) * (4 * NEX / ratio_sampling) / (2.d0 * PI * r))

  print *
  print *,'! period min for ',sngl(RESOLUTION_TARGET),' points per lambda S min in horizontal direction = ',sngl(period_min)
  print *
  print *,'! element width = ',sngl(360.d0 / (4*NEX)),' degrees = ',sngl(2 * PI * R_EARTH / (4*NEX*1000.d0)),' km'
  print *

! define percentage of smallest distance between GLL points for NGLLX points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

  if (NGLLX > NGLLX_MAX_STABILITY) stop 'cannot estimate the stability condition for that degree'

! list of corners defining the edges
! the edge number is sorted according to the numbering convention defined in file hex_nodes.f90
! as well as in DATA/util/YYYYYYYYYYYYYYYYYYYYYYYYYYY DK DK UGLY YYYYYYYYYYYYYYYYYYY

  list_corners_edge( 1,1) = 1
  list_corners_edge( 1,2) = 2

  list_corners_edge( 2,1) = 2
  list_corners_edge( 2,2) = 3

  list_corners_edge( 3,1) = 3
  list_corners_edge( 3,2) = 4

  list_corners_edge( 4,1) = 4
  list_corners_edge( 4,2) = 1

  write(filename,"('prem_horizontal_sampling_',i4.4,'.dat')") NEX
  open(unit=27,file=filename,status='unknown')

  do ipoin = NINIT,1000

  r = R_EARTH * dble(ipoin)/1000.d0

  call prem_iso(r,rho,vp,vs,ratio_sampling)
! for outer core, in which there is no S wave velocity, use P velocity potential with a margin to be safe
  if (vs < 0.001d0) vs = vp / MARGIN_SAFE

! compute mesh size in the horizontal direction at that depth
  spectral_element_size = (2.d0 * PI * r) / (4 * NEX / ratio_sampling)
  average_mesh_point_size = spectral_element_size / (NGLLX - 1)

! compute number of elements per wavelength at that depth
  num_points_per_lambda_S = period_min * vs / average_mesh_point_size

! compute stability condition (Courant number)
  stability = vp * DT / (spectral_element_size * percent_GLL(NGLLX))

  r_gnuplot = (R_EARTH-r)/1000.d0
  write(27,*) sngl(r_gnuplot),sngl(num_points_per_lambda_S * vp/vs),sngl(num_points_per_lambda_S),sngl(stability)

  enddo

  close(27)

!------------------------------------------------
!--- mesh resolution in the radial direction
!------------------------------------------------

  ner(1) = NER_220_CRUST
  ner(2) = NER_400_220
  ner(3) = NER_600_400
  ner(4) = NER_670_600
  ner(5) = NER_771_670
  ner(6) = NER_TOPDDOUBLEPRIME_771
  ner(7) = NER_CMB_TOPDDOUBLEPRIME
  ner(8) = NER_OUTER_CORE
  ner(9) = NER_TOP_CENTRAL_CUBE_ICB

! define the top and bottom radii of all the regions of the mesh in the radial direction
! the first region is the crust at the surface of the Earth
! the last region is in the inner core near the center of the Earth
  r_top(1) = R_EARTH
  r_bottom(1) = R220

  r_top(2) = R220
  r_bottom(2) = R400

  r_top(3) = R400
  r_bottom(3) = R600

  r_top(4) = R600
  r_bottom(4) = R670

  r_top(5) = R670
  r_bottom(5) = R771

  r_top(6) = R771
  r_bottom(6) = RTOPDDOUBLEPRIME

  r_top(7) = RTOPDDOUBLEPRIME
  r_bottom(7) = RCMB

  r_top(8) = RCMB
  r_bottom(8) = RICB

  r_top(9) = RICB
  r_bottom(9) = R_CENTRAL_CUBE

  write(filename,"('prem_radial_sampling_',i4.4,'.dat')") NEX
  open(unit=27,file=filename,status='unknown')

! loop on the layers, start from the last one to go from Earth center to surface
  do ilayer = NB_LAYERS_SAMPLING_STUDY,1,-1

! loop on all the elements in this layer
  do ielem = 1,ner(ilayer)

! compute at the bottom of the layer
  ipoin = ielem - 1
  gamma =  dble(ipoin) / dble(ner(ilayer))
  r = r_bottom(ilayer) * (ONE - gamma) + r_top(ilayer) * gamma

! define next point at the top of the layer
  ipoin = ielem - 1 + 1
  gamma =  dble(ipoin) / dble(ner(ilayer))
  r_next = r_bottom(ilayer) * (ONE - gamma) + r_top(ilayer) * gamma

! add security margin to avoid being exactly on an interface
  r = r + SMALL_OFFSET
  r_next = r_next - SMALL_OFFSET

! compute mesh size in the horizontal direction at that depth
  spectral_element_size = r_next - r
  average_mesh_point_size = spectral_element_size / (NGLLX - 1)

  call prem_iso(r,rho,vp,vs,ratio_sampling)
! for outer core, in which there is no S wave velocity, use P velocity potential with a margin to be safe
  if (vs < 0.001d0) vs = vp / MARGIN_SAFE

! compute number of elements per wavelength at that depth
  num_points_per_lambda_S = period_min * vs / average_mesh_point_size

! compute stability condition (Courant number)
  stability = vp * DT / (spectral_element_size * percent_GLL(NGLLX))

!!!!!! DK DK quick hack to detect doublings and show problems with small doubling brick
  call prem_iso(r_next,rho,vp_new,vs_new,ratio_sampling_new)
  if (ratio_sampling_new < ratio_sampling) stability = stability * 3.d0 / 2.d0
!!!!!! DK DK quick hack to detect doublings and show problems with small doubling brick

  r_gnuplot = (R_EARTH-r)/1000.d0
  write(27,*) sngl(r_gnuplot),sngl(num_points_per_lambda_S * vp/vs),sngl(num_points_per_lambda_S),sngl(stability)

! add one more point at the top of the last layer
  if (ielem == ner(ilayer)) then
    call prem_iso(r_next,rho,vp,vs,ratio_sampling)
! for outer core, in which there is no S wave velocity, use P velocity potential with a margin to be safe
    if (vs < 0.001d0) vs = vp / MARGIN_SAFE

! compute number of elements per wavelength at that depth
    num_points_per_lambda_S = period_min * vs / average_mesh_point_size

! compute stability condition (Courant number)
    stability = vp * DT / (spectral_element_size * percent_GLL(NGLLX))

    r_gnuplot = (R_EARTH-r_next)/1000.d0
    write(27,*) sngl(r_gnuplot),sngl(num_points_per_lambda_S * vp/vs),sngl(num_points_per_lambda_S),sngl(stability)

  endif

  enddo

  enddo

  close(27)

!---------------------------------------------------------------

! find element below top of which we should implement the second doubling in the mantle
! locate element closest to optimal value
  distance_min = HUGEVAL
  do ielem = 1,NER_TOPDDOUBLEPRIME_771
    zval = RTOPDDOUBLEPRIME + ielem * (R771 - RTOPDDOUBLEPRIME) / dble(NER_TOPDDOUBLEPRIME_771)
    distance = abs(zval - (R_EARTH - DEPTH_SECOND_DOUBLING_OPTIMAL))
    if (distance < distance_min) then
      elem_doubling_mantle = ielem
      distance_min = distance
      DEPTH_SECOND_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo
  print *,'for second doubling, selected element ',elem_doubling_mantle,' at depth ', &
              DEPTH_SECOND_DOUBLING_REAL,' and distance ',distance_min,' from optimal depth'
  print *

! find element below top of which we should implement the third doubling in the middle of the outer core
! locate element closest to optimal value
  distance_min = HUGEVAL
! start at element number 2 because we need at least one element below for the fourth doubling
! implemented at the bottom of the outer core
  do ielem = 2,NER_OUTER_CORE
    zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
    distance = abs(zval - (R_EARTH - DEPTH_THIRD_DOUBLING_OPTIMAL))
    if (distance < distance_min) then
      elem_doubling_middle_outer_core = ielem
      distance_min = distance
      DEPTH_THIRD_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo
  print *,'for third doubling, selected element ',elem_doubling_middle_outer_core,' at depth ', &
              DEPTH_THIRD_DOUBLING_REAL,' and distance ',distance_min,' from optimal depth'
  print *

! find element below top of which we should implement the fourth doubling in the middle of the outer core
! locate element closest to optimal value
  distance_min = HUGEVAL
! stop one element before the top because we need at least one element above for the third doubling
! implemented in the middle of the outer core
  do ielem = 1,NER_OUTER_CORE-1
    zval = RICB + ielem * (RCMB - RICB) / dble(NER_OUTER_CORE)
    distance = abs(zval - (R_EARTH - DEPTH_FOURTH_DOUBLING_OPTIMAL))
    if (distance < distance_min) then
      elem_doubling_bottom_outer_core = ielem
      distance_min = distance
      DEPTH_FOURTH_DOUBLING_REAL = R_EARTH - zval
    endif
  enddo
  print *,'for fourth doubling, selected element ',elem_doubling_bottom_outer_core,' at depth ', &
              DEPTH_FOURTH_DOUBLING_REAL,' and distance ',distance_min,' from optimal depth'
  print *

! make sure that the two doublings in the outer core are found in the right order
  if (elem_doubling_bottom_outer_core >= elem_doubling_middle_outer_core) &
                  stop 'error in location of the two doublings in the outer core'

! define all the layers for the mesh
  ner( 1) = NER_220_CRUST
  ner( 2) = NER_400_220
  ner( 3) = NER_600_400
  ner( 4) = NER_670_600
  ner( 5) = NER_771_670
  ner( 6) = NER_TOPDDOUBLEPRIME_771 - elem_doubling_mantle
  ner( 7) = elem_doubling_mantle
  ner( 8) = NER_CMB_TOPDDOUBLEPRIME
  ner( 9) = NER_OUTER_CORE - elem_doubling_middle_outer_core
  ner(10) = elem_doubling_middle_outer_core - elem_doubling_bottom_outer_core
  ner(11) = elem_doubling_bottom_outer_core
  ner(12) = NER_TOP_CENTRAL_CUBE_ICB

! value of the doubling ratio in each radial region of the mesh
  ratio_sampling_array(1:6) = 1
  ratio_sampling_array(7:9) = 2
  ratio_sampling_array(10) = 4
  ratio_sampling_array(11:12) = 8

! define the three regions in which we implement a mesh doubling at the top of that region
  this_region_has_a_doubling(:)  = .false.
  this_region_has_a_doubling(7)  = .true.
  this_region_has_a_doubling(10) = .true.
  this_region_has_a_doubling(11) = .true.

! define the top and bottom radii of all the regions of the mesh in the radial direction
! the first region is the crust at the surface of the Earth
! the last region is in the inner core near the center of the Earth
  r_top(1) = R_EARTH
  r_bottom(1) = R220

  r_top(2) = R220
  r_bottom(2) = R400

  r_top(3) = R400
  r_bottom(3) = R600

  r_top(4) = R600
  r_bottom(4) = R670

  r_top(5) = R670
  r_bottom(5) = R771

  r_top(6) = R771
  r_bottom(6) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL

  r_top(7) = R_EARTH - DEPTH_SECOND_DOUBLING_REAL
  r_bottom(7) = RTOPDDOUBLEPRIME

  r_top(8) = RTOPDDOUBLEPRIME
  r_bottom(8) = RCMB

  r_top(9) = RCMB
  r_bottom(9) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL

  r_top(10) = R_EARTH - DEPTH_THIRD_DOUBLING_REAL
  r_bottom(10) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL

  r_top(11) = R_EARTH - DEPTH_FOURTH_DOUBLING_REAL
  r_bottom(11) = RICB

  r_top(12) = RICB
  r_bottom(12) = R_CENTRAL_CUBE

! create the mesh doubling superbrick
  call define_superbrick_2D(x_superbrick,z_superbrick,ibool_superbrick)

! compute total number of spectral elements in the mesh
  nspec = 0

  do ilayer = 1,NUMBER_OF_MESH_LAYERS

    ner_without_doubling = ner(ilayer)

! if there is a doubling at the top of this region, we implement it in the last layer of elements
! and therefore we suppress one layer of regular elements here
    if (this_region_has_a_doubling(ilayer)) ner_without_doubling = ner_without_doubling - 1

    number_basic_elems_horizontal = NEX_XI/ratio_sampling_array(ilayer)

    nspec = nspec + ner_without_doubling * number_basic_elems_horizontal

! If there is a doubling at the top of this region, let us add these elements.
! The superbrick implements a symmetric four-to-two doubling and therefore replaces
! a basic regular block of 2 elements.
! We have imposed that NEX be a multiple of 16 therefore we know that we can always create
! these 2 blocks because NEX_XI / ratio_sampling_array(ilayer) is always divisible by 2.
    if (this_region_has_a_doubling(ilayer)) nspec = nspec + NSPEC_DOUBLING_SUPERBRICK * number_basic_elems_horizontal / 2

  enddo

! compute total number of OpenDX points in the mesh (many are identical, but simpler to duplicate them here)
  npoin = nspec * NGNOD_OPENDX

! corner nodes
  iaddx(1) = 0
  iaddr(1) = 0

  iaddx(2) = 1
  iaddr(2) = 0

  iaddx(3) = 1
  iaddr(3) = 1

  iaddx(4) = 0
  iaddr(4) = 1

!---
!--- create an OpenDX file with the whole mesh and
!--- another OpenDX file with the worst elements for the stability condition (Courant number) only
!---

! allocate arrays to store Courant stability value and mesh resolution for S waves
  allocate(color(nspec))
  allocate(courant_stability_number(nspec))
  allocate(copy_courant_stability_number(nspec))
  allocate(number_points_S_wavelength(nspec))
  allocate(copy_number_points_S_wavelength(nspec))

! write OpenDX header with element data
  write(filename,"('DX_fullmesh_',i4.4,'.dx')") NEX
  open(unit=11,file=filename,status='unknown')
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

  write(filename,"('DX_worst_stability_',i4.4,'.dx')") NEX
  open(unit=12,file=filename,status='unknown')
  write(12,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

  write(filename,"('DX_worst_sampling_',i4.4,'.dx')") NEX
  open(unit=13,file=filename,status='unknown')
  write(13,*) 'object 1 class array type float rank 1 shape 3 items ',npoin,' data follows'

! final spectral element number
  ispec_final = 0

! output global OpenDX points
  do ilayer = 1,NUMBER_OF_MESH_LAYERS

! loop on all the elements
   do ix_elem = 1,NEX_XI,ratio_sampling_array(ilayer)

    ner_without_doubling = ner(ilayer)

! if there is a doubling at the top of this region, we implement it in the last layer of elements
! and therefore we suppress one layer of regular elements here
    if (this_region_has_a_doubling(ilayer)) ner_without_doubling = ner_without_doubling - 1

   do ir_elem = 1,ner_without_doubling

! final spectral element number
   ispec_final = ispec_final + 1

! loop on all the nodes of this element
   do ignod = 1,NGNOD_OPENDX

! define topological coordinates of this mesh point
    ix = (ix_elem - 1) + iaddx(ignod) * ratio_sampling_array(ilayer)
    ir = (ir_elem - 1) + iaddr(ignod)

! compute the actual position of that grid point
    call compute_value_grid_main_mesh(dble(ix),dble(ir),xelm(ignod),zelm(ignod), &
               ANGULAR_WIDTH_XI_RAD,NEX_XI,R_CENTRAL_CUBE, &
               r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer)

    write(11,"(f10.7,1x,f10.7,' 0')") xelm(ignod),zelm(ignod)
    write(12,"(f10.7,1x,f10.7,' 0')") xelm(ignod),zelm(ignod)
    write(13,"(f10.7,1x,f10.7,' 0')") xelm(ignod),zelm(ignod)

! end of loop on all the nodes of this element
  enddo

! scale arrays from unit sphere to real Earth
   xelm(:) = xelm(:) * R_EARTH
   zelm(:) = zelm(:) * R_EARTH

! compute minimum and maximum distance using the 8 corners of the element
   distance_min = + HUGEVAL
   distance_max = - HUGEVAL
   do iedge = 1,NEDGES
     ignod = list_corners_edge(iedge,1)
     ignod2 = list_corners_edge(iedge,2)
     distance = sqrt((xelm(ignod2) - xelm(ignod))**2 + (zelm(ignod2) - zelm(ignod))**2)
     distance_min = min(distance_min,distance)
     distance_max = max(distance_max,distance)
   enddo

! determine P and S velocity at the barycenter of the element
   xval = sum(xelm(:)) / NGNOD_OPENDX
   zval = sum(zelm(:)) / NGNOD_OPENDX
   r = sqrt(xval**2 + zval**2)
   call prem_iso(r,rho,vp,vs,ratio_sampling)
! for outer core, in which there is no S wave velocity, use P velocity potential with a margin to be safe
  if (vs < 0.001d0) vs = vp / MARGIN_SAFE
! if we are in the last region (the inner core), ignore S-wave velocity because we purposely
! make sampling bad there since Poisson's ratio is too high
  if (ilayer == NUMBER_OF_MESH_LAYERS) vs = HUGEVAL

! store stability condition (Courant number) and number of points per S wavelength for that element
   courant_stability_number(ispec_final) = vp * DT / (distance_min * percent_GLL(NGLLX))
   number_points_S_wavelength(ispec_final) = period_min * vs / (distance_max / (NGLLX - 1))

! store region number to color that element

! we are in the crust or mantle
   if (ilayer <= 8) then
     color(ispec_final) = 1

! we are in the inner core
   else if (ilayer == NUMBER_OF_MESH_LAYERS) then
     color(ispec_final) = 3

! we are in the outer core
   else
     color(ispec_final) = 2
   endif

! end of loop on all the regular elements
  enddo
  enddo

! If there is a doubling at the top of this region, let us add these elements.
! The superbrick implements a symmetric four-to-two doubling and therefore replaces
! a basic regular block of 2 elements.
! We have imposed that NEX be a multiple of 16 therefore we know that we can always create
! these 2 blocks because NEX_XI / ratio_sampling_array(ilayer) is always divisible by 2.
    if (this_region_has_a_doubling(ilayer)) then

! the doubling is implemented in the last radial element
      ir_elem = ner(ilayer)

! loop on all the elements in the 2 x 2 blocks
      do ix_elem = 1,NEX_XI,2*ratio_sampling_array(ilayer)

! loop on all the elements in the mesh doubling superbrick
          do ispec = 1,NSPEC_DOUBLING_SUPERBRICK

! final spectral element number
            ispec_final = ispec_final + 1

! loop on all the nodes of this element
            do ignod = 1,NGNOD_OPENDX

! define topological coordinates of this mesh point
              xval = (ix_elem - 1) + x_superbrick(ibool_superbrick(ignod,ispec)) * ratio_sampling_array(ilayer)
              rval = (ir_elem - 1) + z_superbrick(ibool_superbrick(ignod,ispec))

! compute the actual position of that grid point
              call compute_value_grid_main_mesh(xval,rval,xelm(ignod),zelm(ignod), &
                       ANGULAR_WIDTH_XI_RAD,NEX_XI,R_CENTRAL_CUBE, &
                       r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer)

              write(11,"(f10.7,1x,f10.7,' 0')") xelm(ignod),zelm(ignod)
              write(12,"(f10.7,1x,f10.7,' 0')") xelm(ignod),zelm(ignod)
              write(13,"(f10.7,1x,f10.7,' 0')") xelm(ignod),zelm(ignod)

! end of loop on all the nodes of this element
            enddo

! scale arrays from unit sphere to real Earth
   xelm(:) = xelm(:) * R_EARTH
   zelm(:) = zelm(:) * R_EARTH

! compute minimum and maximum distance using the 8 corners of the element
   distance_min = + HUGEVAL
   distance_max = - HUGEVAL
   do iedge = 1,NEDGES
     ignod = list_corners_edge(iedge,1)
     ignod2 = list_corners_edge(iedge,2)
     distance = sqrt((xelm(ignod2) - xelm(ignod))**2 + (zelm(ignod2) - zelm(ignod))**2)
     distance_min = min(distance_min,distance)
     distance_max = max(distance_max,distance)
   enddo

! determine P and S velocity at the barycenter of the element
   xval = sum(xelm(:)) / NGNOD_OPENDX
   zval = sum(zelm(:)) / NGNOD_OPENDX
   r = sqrt(xval**2 + zval**2)
   call prem_iso(r,rho,vp,vs,ratio_sampling)
! for outer core, in which there is no S wave velocity, use P velocity potential with a margin to be safe
  if (vs < 0.001d0) vs = vp / MARGIN_SAFE
! if we are in the last region (the inner core), ignore S-wave velocity because we purposely
! make sampling bad there since Poisson's ratio is too high, therefore compute for Vp instead
  if (ilayer == NUMBER_OF_MESH_LAYERS) vs = HUGEVAL

! store stability condition (Courant number) for that element
   courant_stability_number(ispec_final) = vp * DT / (distance_min * percent_GLL(NGLLX))
   number_points_S_wavelength(ispec_final) = period_min * vs / (distance_max / (NGLLX - 1))

! store region number to color that element

! we are in the crust or mantle
   if (ilayer <= 8) then
     color(ispec_final) = 1

! we are in the inner core
   else if (ilayer == NUMBER_OF_MESH_LAYERS) then
     color(ispec_final) = 3

! we are in the outer core
   else
     color(ispec_final) = 2
   endif

! end of loops on the mesh doubling elements
        enddo
      enddo

    endif

  enddo

!---------------------------------------------

! determine the NSPEC_WORST_ELEMENTS worst elements and store their element number
  copy_courant_stability_number(:) = courant_stability_number(:)
  copy_number_points_S_wavelength(:) = number_points_S_wavelength(:)

  do ispec = 1,NSPEC_WORST_ELEMENTS
    ispec_worst_elements_stability(ispec) = maxloc(copy_courant_stability_number(:),dim=1)
    ispec_worst_elements_sampling(ispec) = minloc(copy_number_points_S_wavelength(:),dim=1)
! set it to a fictitious value to make sure we do not detect it a second time
    copy_courant_stability_number(ispec_worst_elements_stability(ispec)) = - HUGEVAL
    copy_number_points_S_wavelength(ispec_worst_elements_sampling(ispec)) = + HUGEVAL
  enddo

!---------------------------------------------

! write element header
  write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspec,' data follows'

! output global OpenDX elements
  ioffset = 0

  do ispec = 1,nspec

! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
! in the case of OpenDX, node numbers start at zero
    write(11,"(i6,1x,i6,1x,i6,1x,i6)") ioffset+0,ioffset+3,ioffset+1,ioffset+2

    ioffset = ioffset + NGNOD_OPENDX

  enddo

! write element header
! represent the NSPEC_WORST_ELEMENTS worst elements
  write(12,*) 'object 2 class array type int rank 1 shape 4 items ',NSPEC_WORST_ELEMENTS,' data follows'

  do ispec = 1,NSPEC_WORST_ELEMENTS

    ioffset = (ispec_worst_elements_stability(ispec) - 1) * NGNOD_OPENDX

! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
! in the case of OpenDX, node numbers start at zero
    write(12,"(i6,1x,i6,1x,i6,1x,i6)") ioffset+0,ioffset+3,ioffset+1,ioffset+2

  enddo

! write element header
! represent the NSPEC_WORST_ELEMENTS worst elements
  write(13,*) 'object 2 class array type int rank 1 shape 4 items ',NSPEC_WORST_ELEMENTS,' data follows'

  do ispec = 1,NSPEC_WORST_ELEMENTS

    ioffset = (ispec_worst_elements_sampling(ispec) - 1) * NGNOD_OPENDX

! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
! in the case of OpenDX, node numbers start at zero
    write(13,"(i6,1x,i6,1x,i6,1x,i6)") ioffset+0,ioffset+3,ioffset+1,ioffset+2

  enddo

! output OpenDX header for data
! label for quadrangles in OpenDX is "quads"
  write(11,*) 'attribute "element type" string "quads"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

! write element data
! color elements according to the region they belong to in order to represent the mesh
  do ispec = 1,nspec
    write(11,*) color(ispec)
  enddo

! represent the NSPEC_WORST_ELEMENTS worst elements
  write(12,*) 'attribute "element type" string "quads"'
  write(12,*) 'attribute "ref" string "positions"'
  write(12,*) 'object 3 class array type float rank 0 items ',NSPEC_WORST_ELEMENTS,' data follows'

! write element data
! color elements according to the Courant stability number
  do ispec = 1,NSPEC_WORST_ELEMENTS
    write(12,*) courant_stability_number(ispec_worst_elements_stability(ispec))
  enddo

! represent the NSPEC_WORST_ELEMENTS worst elements
  write(13,*) 'attribute "element type" string "quads"'
  write(13,*) 'attribute "ref" string "positions"'
  write(13,*) 'object 3 class array type float rank 0 items ',NSPEC_WORST_ELEMENTS,' data follows'

! write element data
! color elements according to the number of points per S wavelength
  do ispec = 1,NSPEC_WORST_ELEMENTS
    write(13,*) number_points_S_wavelength(ispec_worst_elements_sampling(ispec))
  enddo

! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

  write(12,*) 'attribute "dep" string "connections"'
  write(12,*) 'object "irregular positions irregular connections" class field'
  write(12,*) 'component "positions" value 1'
  write(12,*) 'component "connections" value 2'
  write(12,*) 'component "data" value 3'
  write(12,*) 'end'

  write(13,*) 'attribute "dep" string "connections"'
  write(13,*) 'object "irregular positions irregular connections" class field'
  write(13,*) 'component "positions" value 1'
  write(13,*) 'component "connections" value 2'
  write(13,*) 'component "data" value 3'
  write(13,*) 'end'

  close(11)
  close(12)
  close(13)

! print information about the maximum Courant number and minimum number of points per wavelength
! detected in the mesh
  print *,'Minimum number of points per S wavelength detected in the mesh = ',minval(number_points_S_wavelength), &
              ' in element ',minloc(number_points_S_wavelength,dim=1)
  print *
  print *,'Maximum stability condition (Courant number) detected in the mesh = ',maxval(courant_stability_number), &
              ' in element ',maxloc(courant_stability_number,dim=1)
  print *

! deallocate arrays
  deallocate(color)
  deallocate(courant_stability_number)
  deallocate(copy_courant_stability_number)
  deallocate(number_points_S_wavelength)
  deallocate(copy_number_points_S_wavelength)

! end of loop on all the cases to study
  enddo

! create a gnuplot script to display all the results
  open(unit=27,file='plot_mesh_resolution.gnu',status='unknown')
  write(27,*) 'set term x11'
  write(27,*) '######set xrange [0:6371]'
  write(27,*) '## this to be able to see values in the crust better on a larger interval'
  write(27,*) 'set xrange [0:6371]'
  write(27,*) 'set xlabel "Depth (km)"'
  write(27,*) 'set grid'
  write(27,*)

! display the PREM model alone
  write(27,"(a140)") 'plot "prem_model.dat" us 1:2 t "Vp PREM" w l 1, "prem_model.dat" us 1:3 t "Vs PREM" w l 3'
  write(27,*) '########pause -1 "hit key"'
  write(27,*)

! loop on all the cases to study
!!!!!!!!!!!! DK DK only one case    do icase = 1,NUM_CASES
  do icase = 2,2

! define value of NEX for this case
  NEX = NEX_val(icase)

  write(27,*) 'set title "NEX =',NEX,'"'
  write(27,*)

  write(27,*) 'set ylabel "Number of points per S wavelength"'
  write(27,*) 'set yrange [0:14]'
  write(27,*) 'set ytics 0,1,14'
  write(27,100) NEX,NEX
  write(27,*) 'pause -1 "hit key"'
  write(27,*)

  write(27,*) 'set ylabel "Number of points per P wavelength"'
  write(27,*) 'set yrange [0:20]'
  write(27,*) 'set ytics 0,1,20'
  write(27,200) NEX,NEX
  write(27,*) 'pause -1 "hit key"'
  write(27,*)

  write(27,*) 'set ylabel "Stability condition (Courant number)"'
  write(27,*) 'set yrange [0:1]'
  write(27,*) 'set ytics 0,0.1,1'
  write(27,300) NEX,NEX
  write(27,*) 'pause -1 "hit key"'
  write(27,*)

! end of loop on all the cases to study
  enddo

! formats for output
 100 format('plot "prem_radial_sampling_',i4.4,'.dat" us 1:3 t "Vs sampling radial" w linesp 1, &
    & "prem_horizontal_sampling_',i4.4,'.dat" us 1:3 t "Vs sampling horizontal" w l 3, &
    & "cmb.dat" t '''' w l 6, "icb.dat" t '''' w l 6')

 200 format('plot "prem_radial_sampling_',i4.4,'.dat" us 1:2 t "Vp sampling radial" w linesp 1, &
    & "prem_horizontal_sampling_',i4.4,'.dat" us 1:2 t "Vp sampling horizontal" w l 3, &
    & "cmb.dat" t '''' w l 6, "icb.dat" t '''' w l 6')

 300 format('plot "prem_radial_sampling_',i4.4,'.dat" us 1:4 t "Stability radial" w linesp 1, &
    & "prem_horizontal_sampling_',i4.4,'.dat" us 1:4 t "Stability horizontal" w l 3, &
    & "cmb.dat" t '''' w l 6, "icb.dat" t '''' w l 6')

  close(27)

  end program display_prem_sampling_doubling

!---------------------------------------------------------------------------

! include code to define the mesh doubling superbrick

  include "define_superbrick_2D.f90"

!---------------------------------------------------------------------------

  subroutine compute_value_grid_main_mesh(xval,rval,xgrid,zgrid, &
               ANGULAR_WIDTH_XI_RAD,NEX_XI,R_CENTRAL_CUBE,r_top,r_bottom,ner,ilayer)

  implicit none

  include "constants_modified_2D.h"

  integer :: NEX_XI,ner,ilayer

  double precision :: xval,rval,xgrid,zgrid,ANGULAR_WIDTH_XI_RAD,R_CENTRAL_CUBE,r_top,r_bottom

! local variables
  double precision :: xi,gamma,x,rgb,rgt,rn
  double precision :: x_bot,z_bot
  double precision :: x_top,z_top

! full Earth (cubed sphere)
  xi = - ANGULAR_WIDTH_XI_RAD/2.d0 + (xval/dble(NEX_XI))*ANGULAR_WIDTH_XI_RAD
  x = tan(xi)

  gamma = ONE / sqrt(ONE + x*x)

  rgt = (r_top / R_EARTH)*gamma
  rgb = (r_bottom / R_EARTH)*gamma

! define the mesh points on the top and the bottom in the cubed shpere
  x_top = -x*rgt
  z_top = rgt

  x_bot = -x*rgb
  z_bot = rgb

! modify in the inner core to match the central cube instead of a sphere
! this mesh works because even if ellipticity and/or topography are turned on
! at this stage the reference Earth is still purely spherical
! therefore it will always perfectly match the sphere defined above
  if (ilayer == NUMBER_OF_MESH_LAYERS) then ! if we are in the last region (the inner core)

    x_bot = -x
    z_bot = ONE

! rescale central cube to match cubed sphere
    x_bot = x_bot * (R_CENTRAL_CUBE/R_EARTH) / sqrt(3.d0)
    z_bot = z_bot * (R_CENTRAL_CUBE/R_EARTH) / sqrt(3.d0)

  endif

! compute the position of the point
  rn = rval / dble(ner)
  xgrid = x_top*rn + x_bot*(ONE-rn)
  zgrid = z_top*rn + z_bot*(ONE-rn)

  end subroutine compute_value_grid_main_mesh

!---------------------------------------------------------------------------

  subroutine prem_iso(r,rho,vp,vs,ratio_sampling)

  implicit none

  include "constants_modified_2D.h"

  double precision :: r,x,rho,vp,vs,ratio_sampling

  x = r / R_EARTH

!
! PREM
!
!--- inner core
!
  if (r <= RICB) then
    rho=13.0885d0-8.8381d0*x*x
    vp=11.2622d0-6.3640d0*x*x
    vs=3.6678d0-4.4475d0*x*x
    if (IMPLEMENT_FOURTH_DOUBLING) then
      ratio_sampling = 8.d0
    else
      ratio_sampling = 4.d0
    endif
!
!--- outer core
!
  else if (r > RICB .and. r <= RCMB) then
    rho=12.5815d0-1.2638d0*x-3.6426d0*x*x-5.5281d0*x*x*x
    vp=11.0487d0-4.0362d0*x+4.8023d0*x*x-13.5732d0*x*x*x
    vs=0.0d0

    if (IMPLEMENT_FOURTH_DOUBLING) then
! third and fourth doublings
      if (r > R_EARTH - DEPTH_THIRD_DOUBLING_OPTIMAL) then
        ratio_sampling = 2.d0
      else if (r > R_EARTH - DEPTH_FOURTH_DOUBLING_OPTIMAL) then
        ratio_sampling = 4.d0
      else
        ratio_sampling = 8.d0
      endif
    else
! third doubling
      if (r > R_EARTH - DEPTH_THIRD_DOUBLING_OPTIMAL) then
        ratio_sampling = 2.d0
      else
        ratio_sampling = 4.d0
      endif
    endif

!
!--- D" at the base of the mantle
!
  else if (r > RCMB .and. r <= RTOPDDOUBLEPRIME) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=15.3891d0-5.3181d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=6.9254d0+1.4672d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    ratio_sampling = 2.d0
!
!--- mantle: from top of D" to d670
!
  else if (r > RTOPDDOUBLEPRIME .and. r <= R771) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=24.9520d0-40.4673d0*x+51.4832d0*x*x-26.6419d0*x*x*x
    vs=11.1671d0-13.7818d0*x+17.4575d0*x*x-9.2777d0*x*x*x

! second doubling
    if (r > R_EARTH - DEPTH_SECOND_DOUBLING_OPTIMAL) then
      ratio_sampling = 1.d0
    else
      ratio_sampling = 2.d0
    endif

  else if (r > R771 .and. r <= R670) then
    rho=7.9565d0-6.4761d0*x+5.5283d0*x*x-3.0807d0*x*x*x
    vp=29.2766d0-23.6027d0*x+5.5242d0*x*x-2.5514d0*x*x*x
    vs=22.3459d0-17.2473d0*x-2.0834d0*x*x+0.9783d0*x*x*x
    ratio_sampling = 1.d0
!
!--- mantle: above d670
!
  else if (r > R670 .and. r <= R600) then
    rho=5.3197d0-1.4836d0*x
    vp=19.0957d0-9.8672d0*x
    vs=9.9839d0-4.9324d0*x
    ratio_sampling = 1.d0
  else if (r > R600 .and. r <= R400) then
    rho=11.2494d0-8.0298d0*x
    vp=39.7027d0-32.6166d0*x
    vs=22.3512d0-18.5856d0*x
    ratio_sampling = 1.d0
  else if (r > R400 .and. r <= R220) then
    rho=7.1089d0-3.8045d0*x
    vp=20.3926d0-12.2569d0*x
    vs=8.9496d0-4.4597d0*x
    ratio_sampling = 1.d0
  else if (r > R220) then
!! DK DK completely suppressed the crust for PKP study
    rho=2.6910d0+0.6924d0*x
    vp=4.1875d0+3.9382d0*x
    vs=2.1519d0+2.3481d0*x
    ratio_sampling = 1.d0
  endif

  rho=rho*1000.0d0
  vp=vp*1000.0d0
  vs=vs*1000.0d0

  end subroutine prem_iso

