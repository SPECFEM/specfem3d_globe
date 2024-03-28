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
! VPREMOON
!
!   Garcia et al. (2011),
!   Very preliminary reference Moon model, PEPI, 188, 96 - 113.
!   https://doi.org/10.1016/j.pepi.2011.06.015
!
! modified seismic model, starting from table 6 (left side, not geodesic model)
!
! added velocity values for outer core and solid inner core based on:
!
!   Weber et al. (2011),
!   Seismic Detection of the Lunar Core, Science, 331.
!   DOI: 10.1126/science.1199375
!
!--------------------------------------------------------------------------------------------------

  module model_vpremoon_par

  implicit none

  ! Very preliminary reference moon model
  ! radius at physical surface
  double precision, parameter :: VPREMOON_RSURFACE = 1737100.d0

  ! VPREMOON radii
  ! no ocean, same as surface radius (Physical surface)
  double precision, parameter :: VPREMOON_ROCEAN = VPREMOON_RSURFACE ! no ocean
  double precision, parameter :: VPREMOON_RMIDDLE_CRUST = 1725100.d0 ! depth = 12 km
  ! Crust-mantle boundary: at 28km depth corresponds to seismic model, not geodesic model
  double precision, parameter :: VPREMOON_RMOHO = 1709100.d0         ! depth = 28 km (Garcia et al. 2011)
  ! Rheological lithosphere
  double precision, parameter :: VPREMOON_R80  =  1607100.d0         ! depth = 130 km
  ! Thermal lithosphere
  double precision, parameter :: VPREMOON_R220 = 1517100.d0          ! depth = 220 km
  double precision, parameter :: VPREMOON_R400 = 1337100.d0          ! depth = 400 km
  double precision, parameter :: VPREMOON_R600 = 1137100.d0          ! depth = 600 km
  ! alpha-olivine-beta-spinel transition likely not reached in the mantle for Moon (due to smaller pressure)
  double precision, parameter :: VPREMOON_R670 = 1067100.d0          ! depth = 670 km
  ! beta-spinel-gamma-spinel transition likely not reached in the mantle for Moon (due to smaller pressure)
  double precision, parameter :: VPREMOON_R771 =  966100.d0           ! depth = 771 km
  ! lower thermal boundary layer: 480km suggested as partial melt boundary (PMB) Weber et al. (2011)
  double precision, parameter :: VPREMOON_RTOPDDOUBLEPRIME = 480000.d0
  ! Core-mantle boundary
  double precision, parameter :: VPREMOON_RCMB = 380000.d0            ! depth = 1357.1 km (Garcia et al. 2011)
  ! inner core radius based on Weber et al. 2011 suggestion
  ! original VPREMOON model from Garcia would have no constraints on solid inner core.
  ! (current meshfem3D code also needs a "small" inner core)
  double precision, parameter :: VPREMOON_RICB = 240000.d0            ! depth = 1497.1 km

  ! note: Based on Weber et al 2011, we add an inner core with a radius of 240 km and density of 8.0 g/cm^3.
  !       This will increase the average density for the whole core (inner & outer) which should be around 5.2 g/cm^3.
  !       Thus, this moon model actually becomes too heavy. Still, it should allow for seismic simulations.

  ! densities [g/cm^3]
  double precision, parameter :: VPREMOON_RHO_OCEANS    = 1020.d0  ! will not be used
  double precision, parameter :: VPREMOON_RHO_TOP_OC    = 5171.d0  ! outer core density fluid outer core (from VPREMOON)
  double precision, parameter :: VPREMOON_RHO_BOTTOM_OC = 5171.d0  ! outer core density (constant value)

  ! number of layers
  integer, parameter :: NR_VPREMOON_layers = 78  ! check with vpremoon.dat in DATA/moon/ folder

  ! model arrays
  double precision, dimension(:),allocatable :: VPREMOON_radius,VPREMOON_density, &
    VPREMOON_vp,VPREMOON_vs,VPREMOON_Qkappa,VPREMOON_Qmu

  ! storage of initial Qmu profile (used for attenuation setup)
  double precision, dimension(:),allocatable :: VPREMOON_Qmu_original

  end module model_vpremoon_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_vpremoon_broadcast(CRUSTAL)

! standard routine to setup model

  use constants, only: SUPPRESS_CRUSTAL_MESH
  use shared_parameters, only: ONE_CRUST

  use model_vpremoon_par

  implicit none

  logical,intent(in) :: CRUSTAL

  ! local parameters
  double precision :: vp_tmp,vs_tmp,rho_tmp,Qmu_tmp,Qkappa_tmp
  integer :: i

  ! check if already done
  if (allocated(VPREMOON_radius)) return

  ! all processes will define same parameters
  call define_model_vpremoon()

  ! ONE_CRUST uses a single layer value for the crust
  if (ONE_CRUST) then
    ! upper crust is in (3:4)
    rho_tmp = VPREMOON_density(NR_VPREMOON_layers-3)
    vp_tmp = VPREMOON_vp(NR_VPREMOON_layers-3)
    vs_tmp = VPREMOON_vs(NR_VPREMOON_layers-3)
    Qkappa_tmp = VPREMOON_Qkappa(NR_VPREMOON_layers-3)
    Qmu_tmp = VPREMOON_Qmu(NR_VPREMOON_layers-3)
    ! assign all crust (1:6) values to upper crust
    VPREMOON_density(NR_VPREMOON_layers-5:NR_VPREMOON_layers) = rho_tmp
    VPREMOON_vp(NR_VPREMOON_layers-5:NR_VPREMOON_layers) = vp_tmp
    VPREMOON_vs(NR_VPREMOON_layers-5:NR_VPREMOON_layers) = vs_tmp
    VPREMOON_Qkappa(NR_VPREMOON_layers-5:NR_VPREMOON_layers) = Qkappa_tmp
    VPREMOON_Qmu(NR_VPREMOON_layers-5:NR_VPREMOON_layers) = Qmu_tmp
  endif

  ! in case an external crustal model will be superimposed on top, we extend mantle values to the surface
  !
  ! strip the crust and replace it by mantle if we use an external crustal model
  ! note: assumes that the crust is given by 6 layers
  !       see in vpremoon.dat: regolith layer (1:2), upper crust (3:4), lower crust (5:6), the mantle
  if (SUPPRESS_CRUSTAL_MESH .or. CRUSTAL) then
    do i = NR_VPREMOON_layers-5,NR_VPREMOON_layers
      VPREMOON_density(i) = VPREMOON_density(NR_VPREMOON_layers-6)
      VPREMOON_vp(i) = VPREMOON_vp(NR_VPREMOON_layers-6)
      VPREMOON_vs(i) = VPREMOON_vs(NR_VPREMOON_layers-6)
      VPREMOON_Qkappa(i) = VPREMOON_Qkappa(NR_VPREMOON_layers-6)
      VPREMOON_Qmu(i) = VPREMOON_Qmu(NR_VPREMOON_layers-6)
    enddo
  endif

  end subroutine model_vpremoon_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_vpremoon(x,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
                            check_doubling_flag,iregion_code)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV

  use model_vpremoon_par

  implicit none

! input:
! radius r: meters

! output:
! density rho: kg/m^3
! compressional wave speed vp: km/s
! shear wave speed vs: km/s

! note: although we provide Qkappa, it is not further used in the code.
!       SPECFEM3D_GLOBE only implements shear attenuation at the moment.

  double precision,intent(in) :: x
  double precision,intent(inout) :: rho,drhodr
  double precision,intent(inout) :: vp,vs,Qmu,Qkappa
  integer,intent(in) :: idoubling
  integer,intent(in) :: iregion_code
  logical,intent(in) :: CRUSTAL,check_doubling_flag

  ! local parameters
  double precision :: r,frac,scaleval,dr,drho
  integer :: i

  ! first allocates and defines model arrays (will be done only once)
  if (.not. allocated(VPREMOON_radius)) then
    call model_vpremoon_broadcast(CRUSTAL)
  endif

  ! compute real physical radius in meters
  r = x * R_PLANET

  ! check flags to make sure we correctly honor the discontinuities
  ! we use strict inequalities since r has been slightly changed in mesher

  ! note: using stop statements, not exit_mpi() calls to avoid the need for MPI libraries when linking xcreate_header_file

  if (check_doubling_flag) then
    !
    !--- inner core
    !
    if (r >= 0.d0 .and. r < VPREMOON_RICB) then
      ! debug
      !if (myrank == 0) print *,'debug: model Sohl:',idoubling,IFLAG_INNER_CORE_NORMAL,IFLAG_MIDDLE_CENTRAL_CUBE, &
      !                            IFLAG_BOTTOM_CENTRAL_CUBE,IFLAG_TOP_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE

      ! checks with inner core flags
      if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
          idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
           stop 'wrong doubling flag for inner core point in model_vpremoon()'
    !
    !--- outer core
    !
    else if (r > VPREMOON_RICB .and. r < VPREMOON_RCMB) then
      if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
        stop 'wrong doubling flag for outer core point in model_vpremoon()'
    !
    !--- D" at the base of the mantle
    !
    else if (r > VPREMOON_RCMB .and. r < VPREMOON_RTOPDDOUBLEPRIME) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for D" point in model_vpremoon()'
    !
    !--- mantle: from top of D" to d670
    !
    else if (r > VPREMOON_RTOPDDOUBLEPRIME .and. r < VPREMOON_R670) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for top D" - > d670 point in model_vpremoon()'
    !
    !--- mantle: from d670 to d220
    !
    else if (r > VPREMOON_R670 .and. r < VPREMOON_R220) then
      if (idoubling /= IFLAG_670_220) &
        stop 'wrong doubling flag for d670 - > d220 point in model_vpremoon()'
    !
    !--- mantle and crust: from d220 to MOHO and then to surface
    !
    else if (r > VPREMOON_R220) then
      if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
        stop 'wrong doubling flag for d220 - > Moho - > surface point in model_vpremoon()'
    endif
  endif

  ! searches radius index
  i = 1
  do while(r >= VPREMOON_radius(i) .and. i /= NR_VPREMOON_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_VPREMOON_layers) i = NR_VPREMOON_layers

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (iregion_code /= 0) then
    ! check
    if (iregion_code /= IREGION_INNER_CORE .and. &
        iregion_code /= IREGION_OUTER_CORE .and. &
        iregion_code /= IREGION_CRUST_MANTLE) stop 'Invalid iregion_code in model_vpremoon() routine'

    ! inner core range [1:2]
    if (iregion_code == IREGION_INNER_CORE .and. i > 2) i = 2

    ! outer core range [3:4]
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_OUTER_CORE .and. i <= 3) i = 3 + 1
    if (iregion_code == IREGION_OUTER_CORE .and. i > 4) i = 4

    ! crust/mantle range [5:78] NR_VPREMOON_layers == 78
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_CRUST_MANTLE .and. i <= 5) i = 5 + 1
    if (iregion_code == IREGION_CRUST_MANTLE .and. i > NR_VPREMOON_layers) i = NR_VPREMOON_layers
  endif

  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    drhodr = 0.d0
    rho = VPREMOON_density(i)
    vp = VPREMOON_vp(i)
    vs = VPREMOON_vs(i)
    Qmu = VPREMOON_Qmu(i)
    Qkappa = VPREMOON_Qkappa(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = VPREMOON_radius(i) - VPREMOON_radius(i-1)
    if (abs(dr) > 1.d-9) then
      drho = (VPREMOON_density(i)-VPREMOON_density(i-1))
      drhodr = drho/dr
      frac = (r - VPREMOON_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      drhodr = 0.d0
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = VPREMOON_density(i-1) + frac * (VPREMOON_density(i)-VPREMOON_density(i-1))
    vp = VPREMOON_vp(i-1) + frac * (VPREMOON_vp(i)-VPREMOON_vp(i-1))
    vs = VPREMOON_vs(i-1) + frac * (VPREMOON_vs(i)-VPREMOON_vs(i-1))
    Qmu = VPREMOON_Qmu(i-1) + frac * (VPREMOON_Qmu(i)-VPREMOON_Qmu(i-1))
    Qkappa = VPREMOON_Qkappa(i-1) + frac * (VPREMOON_Qkappa(i)-VPREMOON_Qkappa(i-1))
  endif

! make sure Vs is zero in the outer core even if roundoff errors on depth
! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if (iregion_code == IREGION_OUTER_CORE) then
    vs = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  rho = rho*1000.0d0/RHOAV
  vp = vp*1000.0d0/(R_PLANET*scaleval)
  vs = vs*1000.0d0/(R_PLANET*scaleval)

  end subroutine model_vpremoon

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_model_vpremoon()

  use constants, only: myrank,SUPPRESS_CRUSTAL_MESH,IIN
  use shared_parameters, only: RICB,RCMB,ROCEAN
  use model_vpremoon_par

  implicit none

  ! local parameters
  integer :: nlines,ir,ier
  double precision :: r,vp,vs,rho
  double precision :: Qp,Qs,Qkappa,Qmu
  double precision :: Qkappa_inv,Qp_inv,dr
  character(len=256) :: datafile,line

  double precision, parameter :: TOL = 1.d-9
  ! debugging
  logical, parameter :: DEBUG = .false.

  ! file name
  datafile = 'DATA/moon/vpremoon.dat'

  ! reads in data file
  open(IIN,file=trim(datafile),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(datafile)
    stop 'Error opening vpremoon.dat file'
  endif

  ! checks number of line
  nlines = 0
  do while(ier == 0)
    read(IIN,"(a256)",iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    line = adjustl(line)        ! suppress leading white spaces, if any
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    nlines = nlines + 1
  enddo
  rewind(IIN)

  ! safety check
  if (nlines /= NR_VPREMOON_layers) then
    print *,'Error: number of data entries ',nlines,' should match number of layers ',NR_VPREMOON_layers,' for vpremoon model'
    print *,'Please check file: ',trim(datafile)
    stop 'Invalid number of data line entries in vpremoon.dat'
  endif

  ! allocates model arrays
  allocate(VPREMOON_radius(NR_VPREMOON_layers), &
           VPREMOON_density(NR_VPREMOON_layers), &
           VPREMOON_vp(NR_VPREMOON_layers), &
           VPREMOON_vs(NR_VPREMOON_layers), &
           VPREMOON_Qkappa(NR_VPREMOON_layers), &
           VPREMOON_Qmu(NR_VPREMOON_layers), &
           VPREMOON_Qmu_original(NR_VPREMOON_layers), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating VPREMOON arrays'
  VPREMOON_radius(:) = 0.d0; VPREMOON_density(:) = 0.d0
  VPREMOON_vp(:) = 0.d0; VPREMOON_vs(:) = 0.d0
  VPREMOON_Qkappa(:) = 0.d0; VPREMOON_Qmu(:) = 0.d0; VPREMOON_Qmu_original(:) = 0.d0

  ! read in layering
  ! note: table 6 starts from top surface to inner core
  !       for the model arrays here, we choose an opposite order, with index 1 being the center of Moon
  !       and then increasing radius.
  nlines = 0
  do while(ier == 0)
    ! read line
    read(IIN,"(a256)",iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    line = adjustl(line)        ! suppress leading white spaces, if any
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    ! read values from line
    ! format: # radius(km) #vp(km/s) #vs(km/s) #density(kg/cm^3) #Qp  #Qs
    read(line,*,iostat=ier) r,vp,vs,rho,Qp,Qs

    if (ier /= 0) then
      print *,'Error: reading vpremoon.dat data, on line ',nlines
      print *,'       data line: ',trim(line)
      print *,'Please check data file vpremoon.dat'
      stop 'Error reading data in vpremoon data'
    endif

    if (r < 0.d0) stop 'Error found negative r value'
    if (vs < 0.d0) stop 'Error found negative vs value'
    if (vp <= 0.d0) stop 'Error found negative or zero vp value'
    if (Qp <= 0.d0) stop 'Error found negative or zero Qp value'
    if (Qs <= 0.d0) stop 'Error found negative or zero Qs value'

    ! converts Qp,Qs to Qkappa,Qmu
    ! see for example Dahlen and Tromp 1998, formulas 9.59 and 9.60
    Qmu = Qs

    if (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2) <= 0.d0) then
      print *,'Error: invalid vp,vs = ',vp,vs,' leads to factor ',(1.d0 - 4.d0/3.d0*(vs**2)/(vp**2))
      print *,'       data line: ',trim(line)
      stop 'Error invalid vp,vs values'
    endif

    ! Qkappa converted from (Qp,Qs)
    ! note: in case Qp == Qs then it follows that Qkappa = Qs
    Qkappa_inv = (1.d0/Qp - 4.d0/3.d0*(vs**2)/(vp**2)/Qs) / (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2))
    if (Qkappa_inv <= 0.d0) stop 'Error got negative or zero 1/Qkappa value'

    Qkappa = 1.d0/Qkappa_inv

    ! check back with Qp
    if (DEBUG) then
      Qp_inv = (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2)) * 1.d0/Qkappa + 4.d0/3.d0*(vs**2)/(vp**2) * 1.d0/Qmu
      ! debug output
      if (myrank == 0) &
        print *,'debug: vpremoon values ',r,vp,vs,rho,Qp,Qs
        print *,'debug: Qkappa,Qmu = ',Qkappa,Qmu,' check Qp',Qp,1.d0/Qp_inv
      ! check
      if (abs(Qp - 1.d0/Qp_inv) > 1.d-9) stop 'Error invalid Qp check'
    endif

    ! increases layer number count
    nlines = nlines + 1

    ! define all the values in the model
    ! change ordering from core to top
    ir = NR_VPREMOON_layers - (nlines-1)
    if (ir < 1) stop 'Error reading in vpremoon data lines'

    ! we need radius in [m], density in [kg/cm^3], vp,vs in [km/s] and Qkappa and Qmu
    VPREMOON_radius(ir) = r * 1000.d0
    VPREMOON_density(ir) = rho
    VPREMOON_vp(ir) = vp
    VPREMOON_vs(ir) = vs
    VPREMOON_Qkappa(ir) = Qkappa
    VPREMOON_Qmu(ir) = Qmu
  enddo
  close(IIN)

  ! debug
  if (DEBUG) then
    if (myrank == 0) then
      print *,'debug: model VPREMOON'
      print *,'format: #line #r #vp #vs #rho #Qkappa #Qmu'
      do ir = 1,NR_VPREMOON_layers
        print *,ir,sngl(VPREMOON_radius(ir)),sngl(VPREMOON_vp(ir)),sngl(VPREMOON_vs(ir)),sngl(VPREMOON_density(ir)), &
                   sngl(VPREMOON_Qkappa(ir)),sngl(VPREMOON_Qmu(ir))
      enddo
      print *
    endif
  endif

  ! checks radius is increasing
  do ir = 2,NR_VPREMOON_layers
    dr = VPREMOON_radius(ir) - VPREMOON_radius(ir-1)
    ! slightly negative steps should not occur
    if (dr < - TOL) then
      print *,'Error layering inversed: radius should decrease in file vpremoon.dat and start from top down to center.'
      print *,'  wrong layering: ',ir,VPREMOON_radius(ir),'and ',ir-1,VPREMOON_radius(ir-1)
      stop 'Error layering in vpremoon.dat has wrong ordering'
    endif
  enddo

  ! checks data table
  ! ICB at radius(33)
  if (abs(RICB - VPREMOON_radius(2)) > TOL) &
    stop 'Error: vpremoon radius RICB and model radius index not matching'
  if (abs(RCMB - VPREMOON_radius(4)) > TOL) &
    stop 'Error: vpremoon radius RCMB and model radius index not matching'
  if (abs(ROCEAN - VPREMOON_radius(NR_VPREMOON_layers)) > TOL) &
    stop 'Error: vpremoon radius ROCEAN and model radius index not matching'

  ! stores Qmu original values (without CRUSTAL modifications) for further use in attenuation routine
  VPREMOON_Qmu_original(:) = VPREMOON_Qmu(:)

  end subroutine define_model_vpremoon

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_vpremoon_density(x,rho)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,CRUSTAL

  use model_vpremoon_par

  implicit none

  double precision,intent(in) :: x
  double precision,intent(out) :: rho

  ! local parameters
  integer :: i
  double precision :: r,dr,frac

  ! first allocates and defines model arrays (will be done only once)
  if (.not. allocated(VPREMOON_radius)) then
    call model_vpremoon_broadcast(CRUSTAL)
  endif

  ! initializes
  rho = 0.d0

  ! compute real physical radius in meters
  r = x * R_PLANET

  i = 1
  do while(r >= VPREMOON_radius(i) .and. i /= NR_VPREMOON_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_VPREMOON_layers) i = NR_VPREMOON_layers

  ! note: this routine helps to sample the density profile of the 1D model
  !       we might assign the "upper" value for points on interfaces like ICB, CMB etc.
  !       this effect should not matter much though since the sampling rate is quite high along the profile.

  ! calculates density according to radius
  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    rho = VPREMOON_density(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = VPREMOON_radius(i) - VPREMOON_radius(i-1)
    ! fraction for linear interpolation
    if (abs(dr) > 1.d-9) then
      frac = (r - VPREMOON_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = VPREMOON_density(i-1) + frac * (VPREMOON_density(i)-VPREMOON_density(i-1))
  endif

  ! non-dimensionalize
  rho = rho*1000.0d0/RHOAV

  end subroutine model_vpremoon_density
