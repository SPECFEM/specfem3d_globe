!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

!--------------------------------------------------------------------------------------------------
! SGLOBE-rani crustal model: modified CRUST2.0 by Chang et al. (2015)
!
! Sung-Joon Chang, Ana M. G. Ferreira, Jeroen Ritsema, Hendrik J. van Heijst, and John H. Woodhouse.
! Joint inversion for global isotropic and radially anisotropic mantle structure including crustal thickness perturbations.
! JGR Solid Earth, 10.1002/2014JB011824.
!
!
! From original Crust2.0, the 7 crustal layers:
! ====================
! 1) ice
! 2) water
! 3) soft sediments
! 4) hard sediments
! 5) upper crust
! 6) middle crust
! 7) lower crust
! + Parameters VP, VS and rho are given explicitly for these 7 layers as well as the mantle below the Moho.
!
! the SGLOBE-rani crustal model modifies crustal thickness compared to CRUST2.0, the rest follows the default
! setup as for the crust2.0 crustal model.
!--------------------------------------------------------------------------------------------------

  module model_sglobecrust_par

  ! crustal_model_constants
  ! modified crustal model parameters for crust2.0
  integer, parameter :: CRUST_NP = 8
  integer, parameter :: CRUST_NLO = 16200
  integer, parameter :: CRUST_NLA = 180

  ! model_crust_variables
  ! Vp, Vs and density
  double precision, dimension(:,:), allocatable :: crust_vp,crust_vs,crust_rho
  character(len=5),dimension(:,:),allocatable :: abbreviation
  character(len=5),dimension(:),allocatable :: code

  ! layer thickness
  double precision, dimension(:,:), allocatable :: crust_thickness

  ! hash table
  integer, dimension(:), allocatable :: crustalhash_to_key

  end module model_sglobecrust_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_sglobecrust_broadcast()

! standard routine to setup model

  use constants
  use model_sglobecrust_par

  implicit none

  integer :: ier,i,ihash

  ! allocate crustal arrays
  allocate(crust_thickness(CRUST_NP,CRUST_NLO), &
           crust_vp(CRUST_NP,CRUST_NLO), &
           crust_vs(CRUST_NP,CRUST_NLO), &
           crust_rho(CRUST_NP,CRUST_NLO), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating crustal arrays')

  ! initializes
  crust_vp(:,:) = ZERO
  crust_vs(:,:) = ZERO
  crust_rho(:,:) = ZERO
  crust_thickness(:,:) = ZERO

  ! allocates arrays
  allocate(abbreviation(CRUST_NLA/2,CRUST_NLA), &
           code(CRUST_NLO),stat=ier)
  if (ier /= 0) stop 'Error allocating abbrev.. arrays'
  abbreviation(:,:) = ''
  code(:) = ''

  ! the variables read are declared and stored in structure model_sglobecrust_par
  if (myrank == 0) call read_sglobecrust_model()

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(crust_thickness,CRUST_NLO*CRUST_NP)
  call bcast_all_dp(crust_vp,CRUST_NLO*CRUST_NP)
  call bcast_all_dp(crust_vs,CRUST_NLO*CRUST_NP)
  call bcast_all_dp(crust_rho,CRUST_NLO*CRUST_NP)

  call bcast_all_ch_array2(abbreviation,CRUST_NLA/2,CRUST_NLA,5)
  call bcast_all_ch_array(code,CRUST_NLO,5)

  ! fill in the hash table
  allocate(crustalhash_to_key(32401),stat=ier)
  if (ier /= 0) stop 'Error allocating crustalhash table'
  crustalhash_to_key(:) = -1
  do i = 1,CRUST_NLO
    call hash_sglobecrust_type(code(i), ihash)
    if (crustalhash_to_key(ihash) /= -1) stop 'Error in sglobecrust_CAPsmoothed: hash table collision'
    crustalhash_to_key(ihash) = i
  enddo

  end subroutine model_sglobecrust_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_sglobecrust(lat,lon,x,vp,vs,rho,moho,sediment,found_crust,elem_in_crust,moho_only)

  use constants
  use shared_parameters, only: R_PLANET_KM,RHOAV

  use model_sglobecrust_par

  implicit none

  double precision,intent(in) :: lat,lon,x
  double precision,intent(out) :: vp,vs,rho,moho,sediment
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust,moho_only

  ! local parameters
  double precision :: thicks_1
  double precision :: h_sed,h_uc
  double precision :: x1,x3,x4,x5,x6,x7
  double precision :: scaleval
  double precision,dimension(CRUST_NP):: vps,vss,rhos,thicks

  ! initializes
  vp = ZERO
  vs = ZERO
  rho = ZERO
  moho = ZERO
  sediment = ZERO
  found_crust = .true.

  ! gets smoothed structure
  call sglobecrust_CAPsmoothed(lat,lon,vps,vss,rhos,thicks,abbreviation, &
                               crust_thickness,crust_vp,crust_vs,crust_rho)

  ! note: for seismic wave propagation in general we ignore the water and ice sheets (oceans are re-added later as an ocean load)
  ! note: but for gravity integral calculations we include the ice
  if (INCLUDE_ICE_IN_CRUST) then
    thicks_1 = thicks(1)
  else
    thicks_1 = ZERO
  endif

  ! whole sediment thickness (with ice if included)
  h_sed = thicks_1 + thicks(3) + thicks(4)

  ! upper crust thickness (including sediments above, and also ice if included)
  h_uc = h_sed + thicks(5)

  ! non-dimensionalization factor
  scaleval = ONE / R_PLANET_KM

  ! non-dimensionalize thicknesses (given in km)
  ! ice
  x1 = ONE - thicks_1 * scaleval
  ! upper sediment
  x3 = ONE - (thicks_1 + thicks(3)) * scaleval
  ! all sediments
  x4 = ONE - h_sed * scaleval
  ! upper crust
  x5 = ONE - h_uc * scaleval
  ! middle crust
  x6 = ONE - (h_uc + thicks(6)) * scaleval
  ! lower crust
  x7 = ONE - (h_uc + thicks(6) + thicks(7)) * scaleval

  ! no matter if found_crust is true or false, compute moho thickness
  moho = (h_uc + thicks(6) + thicks(7)) * scaleval

  ! checks if anything further to do
  if (moho_only) return

  ! sediment thickness
  if (INCLUDE_SEDIMENTS_IN_CRUST) then
    sediment = h_sed * scaleval
  endif

  ! gets corresponding crustal velocities and density
  if (x > x1 .and. INCLUDE_ICE_IN_CRUST) then
    vp = vps(1)
    vs = vss(1)
    rho = rhos(1)
  else if (x > x3 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
    vp = vps(3)
    vs = vss(3)
    rho = rhos(3)
  else if (x > x4 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
    vp = vps(4)
    vs = vss(4)
    rho = rhos(4)
  else if (x > x5) then
    vp = vps(5)
    vs = vss(5)
    rho = rhos(5)
  else if (x > x6) then
    vp = vps(6)
    vs = vss(6)
    rho = rhos(6)
  else if (x > x7 .or. elem_in_crust) then
    ! takes lower crustal values only if x is slightly above moho depth or
    ! if elem_in_crust is set
    !
    ! note: it looks like this does distinguish between GLL points at the exact moho boundary,
    !          where the point is on the interface between both,
    !          oceanic elements and mantle elements below
    vp = vps(7)
    vs = vss(7)
    rho = rhos(7)
  else
    ! note: if x is exactly the moho depth this will return false
    found_crust = .false.
  endif

  ! non-dimensionalize
  if (found_crust) then
    scaleval = ONE / ( R_PLANET_KM * dsqrt(PI*GRAV*RHOAV) )
    vp = vp * scaleval
    vs = vs * scaleval
    rho = rho * 1000.0d0 / RHOAV
  endif

  end subroutine model_sglobecrust

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_sglobecrust_model()

  use constants
  use model_sglobecrust_par

  implicit none

  ! local variables
  integer :: i,ila,icolat,ikey,ier

  double precision :: moho
  double precision :: h_moho_min,h_moho_max
  double precision :: moho_total,h_moho_min_total,h_moho_max_total

  character(len=*), parameter :: CNtype2 = 'DATA/sglobe/CNtype2.txt'
  character(len=*), parameter :: CNtype2_key_modif = 'DATA/sglobe/CNtype2_key_modif.txt'

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: sglobecrust'
  write(IMAIN,*)
  call flush_IMAIN()

  open(unit = IIN,file=CNtype2,status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(CNtype2), '": ', ier
    call flush_IMAIN()
    ! stop
    call exit_MPI(0,'Error model sglobecrust')
  endif

  do ila = 1,CRUST_NLA/2
    read(IIN,*) icolat,(abbreviation(ila,i),i = 1,CRUST_NLA)
  enddo
  close(IIN)

  open(unit = IIN,file=CNtype2_key_modif,status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(CNtype2_key_modif), '": ', ier
    call exit_MPI(0,'Error model sglobecrust')
  endif

  h_moho_min = HUGEVAL
  h_moho_max = -HUGEVAL
  h_moho_min_total = HUGEVAL
  h_moho_max_total = -HUGEVAL

  do ikey = 1,CRUST_NLO
    read(IIN,"(a5)") code(ikey)
    read(IIN,*) (crust_vp(i,ikey),i = 1,CRUST_NP)
    read(IIN,*) (crust_vs(i,ikey),i = 1,CRUST_NP)
    read(IIN,*) (crust_rho(i,ikey),i = 1,CRUST_NP)
    read(IIN,*) (crust_thickness(i,ikey),i = 1,CRUST_NP-1),crust_thickness(CRUST_NP,ikey)

    ! thickness = ice (layer index 1) + sediment (index 3+4) + crystalline crust (index 5+6+7)
    ! crustal thickness without ice
    ! note: etopo1 has topography including ice ("ice surface" version) and at base of ice sheets ("bedrock" version)
    !       see: http://www.ngdc.noaa.gov/mgg/global/global.html
    moho = crust_thickness(3,ikey) + crust_thickness(4,ikey) &
           + crust_thickness(5,ikey) + crust_thickness(6,ikey) + crust_thickness(7,ikey)

    ! limit moho thickness
    if (moho > h_moho_max) h_moho_max = moho
    if (moho < h_moho_min) h_moho_min = moho

    ! last entry interpreted as total crustal thickness
    moho_total = crust_thickness(CRUST_NP,ikey)
    if (moho_total > h_moho_max_total) h_moho_max_total = moho_total
    if (moho_total < h_moho_min_total) h_moho_min_total = moho_total

  enddo
  close(IIN)

  ! user output
  write(IMAIN,*) '  Moho crustal thickness (without ice) min/max = ',sngl(h_moho_min),sngl(h_moho_max),' km'
  write(IMAIN,*) '               total crustal thickness min/max = ',sngl(h_moho_min_total),sngl(h_moho_max_total),' km'
  write(IMAIN,*)
  call flush_IMAIN()

  ! checks
  if (h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) stop 'incorrect moho depths in read_sglobecrust_model'

  end subroutine read_sglobecrust_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sglobecrust_CAPsmoothed(lat,lon,velp,vels,rho,thick,abbreviation,crust_thickness,crust_vp,crust_vs,crust_rho)

! crustal vp and vs in km/s, layer thickness in km
!
! sglobecrust based on crust2.0 gets smoothed with a cap of size CAP using NTHETA points
! in the theta direction and NPHI in the phi direction.
! The cap is first rotated to the North Pole for easier implementation.

  use constants
  use model_sglobecrust_par, only: CRUST_NP,CRUST_NLO,CRUST_NLA,crustalhash_to_key

  implicit none

  ! sampling rate for CAP points
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 20

  ! argument variables
  double precision :: lat,lon
  double precision,dimension(CRUST_NP) :: rho,thick,velp,vels
  double precision,dimension(CRUST_NP,CRUST_NLO) :: crust_thickness,crust_vp,crust_vs,crust_rho

  character(len=5) :: abbreviation(CRUST_NLA/2,CRUST_NLA)

  !-------------------------------
  ! work-around to avoid Jacobian problems when stretching mesh elements;
  ! one could also try to slightly change the shape of the doubling element bricks (which cause the problem)...
  !
  ! defines a "critical" region around the andes to have at least a 2-degree smoothing;
  ! critical region can lead to negative Jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_ANDES = -20.0d0
  double precision,parameter :: LON_CRITICAL_ANDES = -70.0d0
  double precision,parameter :: CRITICAL_RANGE = 70.0d0
  !-------------------------------

  ! local variables
  double precision :: xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)
  double precision :: rhol(CRUST_NP),thickl(CRUST_NP),velpl(CRUST_NP),velsl(CRUST_NP)

  double precision :: weightl,cap_degree
  double precision :: dist
  double precision :: h_sed
  integer :: i,icolat,ilon
  character(len=5) :: crustaltype

  ! small hash table to convert crustal types to key
  integer :: ihash, crustalkey

  ! checks latitude/longitude
  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
    print *,'Error in lat/lon:',lat,lon
    stop 'Error in latitude/longitude range in sglobecrust'
  endif

  ! makes sure lat/lon are within sglobecrust range
  if (lat == 90.0d0) lat=89.9999d0
  if (lat == -90.0d0) lat=-89.9999d0
  if (lon == 180.0d0) lon=179.9999d0
  if (lon == -180.0d0) lon=-179.9999d0

  ! sets up smoothing points based on cap smoothing
  cap_degree = CAP_SMOOTHING_DEGREE_DEFAULT

  ! checks if inside/outside of critical region for mesh stretching
  if (SMOOTH_CRUST_EVEN_MORE) then
    dist = dsqrt( (lon-LON_CRITICAL_ANDES)**2 + (lat-LAT_CRITICAL_ANDES )**2 )
    if (dist < CRITICAL_RANGE) then
      ! increases cap smoothing degree
      ! scales between -1 at center and 0 at border
      dist = dist / CRITICAL_RANGE - ONE
      ! shifts value to 1 at center and 0 to the border with exponential decay
      dist = ONE - exp( - dist*dist*10.0d0 )
      ! increases smoothing degree inside of critical region to 2 degree
      cap_degree = cap_degree + dist
    endif
  endif

  ! gets smoothing points and weights
  call smooth_weights_CAP_vardegree(lon,lat,xlon,xlat,weight,cap_degree,NTHETA,NPHI)

  ! initializes
  velp(:) = ZERO
  vels(:) = ZERO
  rho(:) = ZERO
  thick(:) = ZERO

  ! loops over weight points
  do i = 1,NTHETA*NPHI
    ! gets lat/lon indices
    call sglobecrust_icolat_ilon(xlat(i),xlon(i),icolat,ilon)

    crustaltype = abbreviation(icolat,ilon)

    call hash_sglobecrust_type(crustaltype, ihash)
    crustalkey = crustalhash_to_key(ihash)
    if (crustalkey == -1) stop 'Error in retrieving crust type key'

    ! gets crust values
    call get_sglobecrust_structure(crustalkey,velpl,velsl,rhol,thickl, &
                                   crust_thickness,crust_vp,crust_vs,crust_rho)

    ! sediment thickness
    h_sed = thickl(3) + thickl(4)

    ! takes upper crust value if sediment too thin
    if (h_sed < MINIMUM_SEDIMENT_THICKNESS) then
      velpl(3) = velpl(5)
      velpl(4) = velpl(5)

      velsl(3) = velsl(5)
      velsl(4) = velsl(5)

      rhol(3) = rhol(5)
      rhol(4) = rhol(5)
    endif

    ! weighting value
    weightl = weight(i)

    ! total, smoothed values
    rho(:) = rho(:) + weightl*rhol(:)
    thick(:) = thick(:) + weightl*thickl(:)
    velp(:) = velp(:) + weightl*velpl(:)
    vels(:) = vels(:) + weightl*velsl(:)
  enddo

  end subroutine sglobecrust_CAPsmoothed

!
!-------------------------------------------------------------------------------------------------
!

! hash table to define the crustal type using an integer instead of characters
! because working with integers in the rest of the routines results in much faster code

  subroutine hash_sglobecrust_type(crustaltype, ihash)

  implicit none

  character(len=5), intent(in) :: crustaltype
  integer, intent(out) :: ihash
  ! local parameters
  integer :: ier

  ! original crust2.0 uses crustal types (CNtype2_key_modif.txt) with character codes, like D0, D1, .., YN, YP
  !ihash = iachar(crustaltype(1:1)) + 128*iachar(crustaltype(2:2)) + 1

  ! the new modified SGLOBE crust uses numbers as types, like 16201, 16202, .., 32399, 32400
  read(crustaltype,*,iostat=ier) ihash
  if (ier /= 0) stop 'Error reading ihash from crustaltype in model sglobecrust'

  end subroutine hash_sglobecrust_type

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_sglobecrust_structure(ikey,vptyp,vstyp,rhtyp,thtp,crust_thickness,crust_vp,crust_vs,crust_rho)

  use model_sglobecrust_par, only: CRUST_NP,CRUST_NLO

  implicit none

  ! argument variables
  integer, intent(in) :: ikey
  double precision, intent(out) :: rhtyp(CRUST_NP),thtp(CRUST_NP)
  double precision, intent(out) :: vptyp(CRUST_NP),vstyp(CRUST_NP)
  double precision :: crust_thickness(CRUST_NP,CRUST_NLO),crust_vp(CRUST_NP,CRUST_NLO)
  double precision :: crust_vs(CRUST_NP,CRUST_NLO),crust_rho(CRUST_NP,CRUST_NLO)

  ! local variables
  integer :: i

  ! set vp,vs and rho for all layers
  do i = 1,CRUST_NP
    vptyp(i) = crust_vp(i,ikey)
    vstyp(i) = crust_vs(i,ikey)
    rhtyp(i) = crust_rho(i,ikey)
    thtp(i)  = crust_thickness(i,ikey)
  enddo

  ! get distance to Moho from the bottom of the ocean or the ice
  ! value could be used for checking, but is unused so far...
  thtp(CRUST_NP) = thtp(CRUST_NP) - thtp(1) - thtp(2)

  end subroutine get_sglobecrust_structure

!
!-------------------------------------------------------------------------------------------------
!

  subroutine sglobecrust_icolat_ilon(xlat,xlon,icolat,ilon)

  implicit none

! argument variables
  double precision :: xlat,xlon
  integer :: icolat,ilon

  if (xlat > 90.0d0 .or. xlat < -90.0d0 .or. xlon > 180.0d0 .or. xlon < -180.0d0) &
    stop 'Error in latitude/longitude range in sglobecrust_icolat_ilon'

  icolat = int(1+( (90.d0-xlat)*0.5d0 ))
  if (icolat == 91) icolat = 90

  ilon = int(1+( (180.d0+xlon)*0.5d0 ))
  if (ilon == 181) ilon = 1

  if (icolat > 90 .or. icolat < 1) stop 'Error in routine sglobecrust_icolat_ilon'
  if (ilon < 1 .or. ilon > 180) stop 'Error in routine sglobecrust_icolat_ilon'

  end subroutine sglobecrust_icolat_ilon
