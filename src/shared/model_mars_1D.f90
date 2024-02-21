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

!=====================================================================
!
! Mars 1D model - generic model
!
!=====================================================================

! Mars 1D model defined by table in file DATA/mars/mars_1D.dat
!
! similar implementation and file format as VPREMOON model (in DATA/moon/vpremoon.dat)

  module model_mars_1d_par

  implicit none

  ! radius at physical surface
  double precision, parameter :: MARS_1D_RSURFACE = 3389500.d0

  ! radii
  ! no ocean, same as surface radius (Physical surface)
  double precision, parameter :: MARS_1D_ROCEAN = MARS_1D_RSURFACE ! no ocean
  ! Crust
  double precision, parameter :: MARS_1D_RMIDDLE_CRUST = 3365500.d0  ! depth = 24 km
  ! Crust-mantle boundary
  double precision, parameter :: MARS_1D_RMOHO = 3341500.d0      ! depth = 48 km average crustal thickness
  ! Rheological lithosphere
  double precision, parameter :: MARS_1D_R80  = 3239500.d0       ! depth = 150 km, if smaller will cause negative Jacobian err
  ! Thermal lithosphere
  double precision, parameter :: MARS_1D_R220 = 3139000.d0       ! depth = ~250 km
  double precision, parameter :: MARS_1D_R400 = 2885000.d0       ! depth = ~500 km (added for meshing)
  double precision, parameter :: MARS_1D_R600 = 2635000.d0       ! depth = 750 km (added for meshing)
  ! alpha-olivine-beta-spinel transition
  double precision, parameter :: MARS_1D_R670 = 2485000.d0       ! depth = 900 km
  ! beta-spinel-gamma-spinel transition
  double precision, parameter :: MARS_1D_R771 = 2133000.d0       ! depth = 1357 km, below which the second doubling implemented
  ! lower thermal boundary layer
  double precision, parameter :: MARS_1D_RTOPDDOUBLEPRIME = 1915380.d0 ! ~60km layer thickness to CMB
  ! Core-mantle boundary
  double precision, parameter :: MARS_1D_RCMB = 1855380.d0
  ! inner core radius
  double precision, parameter :: MARS_1D_RICB = 515000.d0          ! good for both stability and efficiency

  ! densities [g/cm^3]
  double precision, parameter :: MARS_1D_RHO_OCEANS    = 1020.d0  ! will not be used
  double precision, parameter :: MARS_1D_RHO_TOP_OC    = 5480.d0  ! density outer core (top)
  double precision, parameter :: MARS_1D_RHO_BOTTOM_OC = 6020.d0  ! density outer core (bottom)

  ! number of layers
  integer :: NR_mars_1D_layers
  integer :: interface_index_moho,interface_index_cmb,interface_index_icb

  ! model arrays
  double precision, dimension(:),allocatable :: mars_1D_radius,mars_1D_density, &
    mars_1D_vp,mars_1D_vs,mars_1D_Qkappa,mars_1D_Qmu

  ! storage of initial Qmu profile (used for attenuation setup)
  double precision, dimension(:),allocatable :: MARS_1D_Qmu_original

  end module model_mars_1d_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_mars_1D(x,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
                            check_doubling_flag,iregion_code)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV

  use model_mars_1d_par

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
  if (.not. allocated(mars_1D_radius)) then
    ! all processes will define same parameters
    call define_model_mars_1D(CRUSTAL)
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
    if (r >= 0.d0 .and. r < MARS_1D_RICB) then
      ! debug
      !if (myrank == 0) print *,'debug: model Sohl:',idoubling,IFLAG_INNER_CORE_NORMAL,IFLAG_MIDDLE_CENTRAL_CUBE, &
      !                            IFLAG_BOTTOM_CENTRAL_CUBE,IFLAG_TOP_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE

      ! checks with inner core flags
      if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
          idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
           stop 'wrong doubling flag for inner core point in model_mars_1D()'
    !
    !--- outer core
    !
    else if (r > MARS_1D_RICB .and. r < MARS_1D_RCMB) then
      if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
        stop 'wrong doubling flag for outer core point in model_mars_1D()'
    !
    !--- D" at the base of the mantle
    !
    else if (r > MARS_1D_RCMB .and. r < MARS_1D_RTOPDDOUBLEPRIME) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for D" point in model_mars_1D()'
    !
    !--- mantle: from top of D" to d670
    !
    else if (r > MARS_1D_RTOPDDOUBLEPRIME .and. r < MARS_1D_R670) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for top D" - > d670 point in model_mars_1D()'
    !
    !--- mantle: from d670 to d220
    !
    else if (r > MARS_1D_R670 .and. r < MARS_1D_R220) then
      if (idoubling /= IFLAG_670_220) &
        stop 'wrong doubling flag for d670 - > d220 point in model_mars_1D()'
    !
    !--- mantle and crust: from d220 to MOHO and then to surface
    !
    else if (r > MARS_1D_R220) then
      if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
        stop 'wrong doubling flag for d220 - > Moho - > surface point in model_mars_1D()'
    endif
  endif

  ! searches radius index
  i = 1
  do while(r >= MARS_1D_radius(i) .and. i /= NR_mars_1D_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_mars_1D_layers) i = NR_mars_1D_layers

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (iregion_code /= 0) then
    ! check
    if (iregion_code /= IREGION_INNER_CORE .and. &
        iregion_code /= IREGION_OUTER_CORE .and. &
        iregion_code /= IREGION_CRUST_MANTLE) stop 'Invalid iregion_code in model_mars_1D() routine'

    ! inner core range [1:interface_index_icb]
    if (iregion_code == IREGION_INNER_CORE .and. i > interface_index_icb) i = interface_index_icb

    ! outer core range [interface_icb+1:interface_cmb]
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_OUTER_CORE .and. i <= (interface_index_icb+1)) i = (interface_index_icb+1) + 1
    if (iregion_code == IREGION_OUTER_CORE .and. i > interface_index_cmb) i = interface_index_cmb

    ! crust/mantle range [interface_cmb+1:NR]
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_CRUST_MANTLE .and. i <= (interface_index_cmb+1)) i = (interface_index_cmb+1) + 1
    if (iregion_code == IREGION_CRUST_MANTLE .and. i > NR_mars_1D_layers) i = NR_mars_1D_layers
  endif

  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    drhodr = 0.d0
    rho = mars_1D_density(i)
    vp = mars_1D_vp(i)
    vs = mars_1D_vs(i)
    Qmu = mars_1D_Qmu(i)
    Qkappa = mars_1D_Qkappa(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = mars_1D_radius(i) - mars_1D_radius(i-1)
    if (abs(dr) > 1.d-9) then
      drho = (mars_1D_density(i)-mars_1D_density(i-1))
      drhodr = drho/dr
      frac = (r - mars_1D_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      drhodr = 0.d0
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = mars_1D_density(i-1) + frac * (mars_1D_density(i)-mars_1D_density(i-1))
    vp = mars_1D_vp(i-1) + frac * (mars_1D_vp(i)-mars_1D_vp(i-1))
    vs = mars_1D_vs(i-1) + frac * (mars_1D_vs(i)-mars_1D_vs(i-1))
    Qmu = mars_1D_Qmu(i-1) + frac * (mars_1D_Qmu(i)-mars_1D_Qmu(i-1))
    Qkappa = mars_1D_Qkappa(i-1) + frac * (mars_1D_Qkappa(i)-mars_1D_Qkappa(i-1))
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

  end subroutine model_mars_1D

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_mars_1D_broadcast(CRUSTAL)

! standard routine to setup model

  implicit none

  logical,intent(in) :: CRUSTAL

  ! all processes will define same parameters
  call define_model_mars_1D(CRUSTAL)

  end subroutine model_mars_1D_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_model_mars_1D(USE_EXTERNAL_CRUSTAL_MODEL)

  use constants, only: myrank,SUPPRESS_CRUSTAL_MESH,IIN
  use shared_parameters, only: R_PLANET,RICB,RCMB,ROCEAN,ONE_CRUST
  use model_mars_1d_par

  implicit none

  logical, intent(in) :: USE_EXTERNAL_CRUSTAL_MODEL

  ! local parameters
  integer :: nlines,ir,i,ier
  double precision :: r,vp,vs,rho
  double precision :: Qp,Qs,Qkappa,Qmu
  double precision :: Qkappa_inv,Qp_inv,dr
  double precision :: vp_tmp,vs_tmp,rho_tmp,Qmu_tmp,Qkappa_tmp
  double precision :: r_top,r_interface1,r_interface2
  character(len=256) :: datafile,line,tmp_line

  ! discontinuity layer - tolerance (discontinuities repeat same depth for upper/lower values)
  double precision, parameter :: TOL = 1.d-9

  ! debugging
  logical, parameter :: DEBUG = .false.

  ! checks if anything to do, only needs to read in model once
  if (allocated(mars_1D_radius)) return

  ! file name
  datafile = 'DATA/mars/mars_1D.dat'

  ! reads in data file
  open(IIN,file=trim(datafile),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(datafile)
    stop 'Error opening mars_1D.dat file'
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

  ! sets number of data entries
  NR_mars_1D_layers = nlines

  ! allocates model arrays
  allocate(mars_1D_radius(NR_mars_1D_layers), &
           mars_1D_density(NR_mars_1D_layers), &
           mars_1D_vp(NR_mars_1D_layers), &
           mars_1D_vs(NR_mars_1D_layers), &
           mars_1D_Qkappa(NR_mars_1D_layers), &
           mars_1D_Qmu(NR_mars_1D_layers), &
           mars_1D_Qmu_original(NR_mars_1D_layers), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating mars_1D arrays'
  mars_1D_radius(:) = 0.d0; mars_1D_density(:) = 0.d0
  mars_1D_vp(:) = 0.d0; mars_1D_vs(:) = 0.d0
  mars_1D_Qkappa(:) = 0.d0; mars_1D_Qmu(:) = 0.d0; mars_1D_Qmu_original(:) = 0.d0

  interface_index_moho = 0
  interface_index_cmb = 0
  interface_index_icb = 0

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
      print *,'Error: reading mars_1D.dat data, on line ',nlines
      print *,'       data line: ',trim(line)
      print *,'Please check data file mars_1D.dat'
      stop 'Error reading data in mars_1D data'
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

    !debug
    !print *,'debug: ',nlines,'line ',trim(line),'leads to Qkappa = ',1.d0/Qkappa_inv
    !print *,'debug: ',nlines,'line Qkappa=9999.9 for Qp ', &
    !  1.d0/((1.d0 - 4.d0/3.d0*(vs**2)/(vp**2)) * 1.d0/9999.9 + 4.d0/3.d0*(vs**2)/(vp**2) * 1.d0/Qmu)

    ! checks for positive value
    if (Qkappa_inv <= 0.d0) then
      print *,'Error: invalid Qp, Qs values ',Qp,Qs,' leads to negative or zero 1/Qkappa value ',Qkappa_inv
      print *,'       please check model in mars_1D.dat, on line ',nlines
      print *,'       data line: ',trim(line)
      stop 'Error got negative or zero 1/Qkappa value'
    endif
    Qkappa = 1.d0/Qkappa_inv

    ! check back with Qp
    if (DEBUG) then
      Qp_inv = (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2)) * 1.d0/Qkappa + 4.d0/3.d0*(vs**2)/(vp**2) * 1.d0/Qmu
      ! debug output
      if (myrank == 0) &
        print *,'debug: mars_1D values ',r,vp,vs,rho,Qp,Qs
        print *,'debug: Qkappa,Qmu = ',Qkappa,Qmu,' check Qp',Qp,1.d0/Qp_inv
      ! check
      if (abs(Qp - 1.d0/Qp_inv) > 1.d-9) stop 'Error invalid Qp check'
    endif

    ! increases layer number count
    nlines = nlines + 1

    ! define all the values in the model
    ! change ordering from core to top
    ir = NR_mars_1D_layers - (nlines-1)
    if (ir < 1) stop 'Error reading in mars_1D data lines'

    ! we need radius in [m], density in [kg/cm^3], vp,vs in [km/s] and Qkappa and Qmu
    mars_1D_radius(ir) = r * 1000.d0
    mars_1D_density(ir) = rho
    mars_1D_vp(ir) = vp
    mars_1D_vs(ir) = vs
    mars_1D_Qkappa(ir) = Qkappa
    mars_1D_Qmu(ir) = Qmu

    ! checks if line is defining a moho/CMB/ICB interface by a comment at the end like:
    !  .. # moho ..
    !  .. # CMB ..
    !  .. # ICB ..
    ! converts all string characters to lowercase (to make user input case-insensitive)
    call convert_to_lowercase(line,tmp_line)
    line = tmp_line
    ! checks index for moho/cmb/icb
    ! note that the line comment .. # moho
    ! could be added on both sides of the interface, e.g., .. # moho crust and .. # moho mantle
    ! this will recognize the lower line as well, that is with a higher nline number .. # moho mantle
    ! and therefore take the smaller index  ir = NR - (nlines-1) as index_moho
    if (index(line,'moho') > 0) interface_index_moho = ir
    if (index(line,'cmb') > 0) interface_index_cmb = ir
    if (index(line,'icb') > 0) interface_index_icb = ir

    !debug
    !print *,"debug: line ",line,interface_index_moho,interface_index_cmb,interface_index_icb
  enddo
  close(IIN)

  ! makes sure that the interface index is for the lower layer, starting from the inner core up to the surface
  ! that is, interface_index_icb should be the uppermost layer in the inner core
  !          interface_index_cmb should be the uppermost layer in the outer core
  !          interface_index_moho should be the uppermost layer in the mantle
  if (interface_index_moho == 0) &
    stop 'Error: Moho interface not recognized in mars_1D.dat, please add a comment .. # moho at the end of the data line'
  if (interface_index_cmb == 0) &
    stop 'Error: CMB interface not recognized in mars_1D.dat, please add a comment .. # CMB at the end of the data line'
  if (interface_index_icb == 0) &
    stop 'Error: ICB interface not recognized in mars_1D.dat, please add a comment .. # ICB at the end of the data line'

  ! moho interface
  ! checks if layer below has the same radius, that is the .. # moho keyword was added at the top layer in the crust
  r_interface1 = mars_1D_radius(interface_index_moho)
  r_interface2 = mars_1D_radius(interface_index_moho-1)
  if (abs(r_interface1 - r_interface2) < TOL) interface_index_moho = interface_index_moho - 1

  ! CMB interface
  ! checks if layer below has the same radius, that is the .. # CMB keyword was added at the top layer in the mantle
  r_interface1 = mars_1D_radius(interface_index_cmb)
  r_interface2 = mars_1D_radius(interface_index_cmb-1)
  if (abs(r_interface1 - r_interface2) < TOL) interface_index_cmb = interface_index_cmb - 1

  ! ICB interface
  ! checks if layer below at the same radius, that is the .. # ICB keyword was added at the top layer in the outer core
  r_interface1 = mars_1D_radius(interface_index_moho)
  r_interface2 = mars_1D_radius(interface_index_moho-1)
  if (abs(r_interface1 - r_interface2) < TOL) interface_index_icb = interface_index_icb - 1

  !debug
  if (DEBUG) then
    if (myrank == 0) then
      print *,'debug: moho interface at index ',interface_index_moho, &
               mars_1D_radius(interface_index_moho),mars_1D_radius(interface_index_moho+1)
      print *,'debug: cmb  interface at index ',interface_index_cmb, &
               mars_1D_radius(interface_index_cmb),mars_1D_radius(interface_index_cmb+1)
      print *,'debug: icb  interface at index ',interface_index_icb, &
              mars_1D_radius(interface_index_icb),mars_1D_radius(interface_index_icb+1)
      print *
      print *,'debug: model mars_1D'
      print *,'format: #line #r #vp #vs #rho #Qkappa #Qmu'
      do ir = 1,NR_mars_1D_layers
        print *,ir,sngl(mars_1D_radius(ir)),sngl(mars_1D_vp(ir)),sngl(mars_1D_vs(ir)),sngl(mars_1D_density(ir)), &
                   sngl(mars_1D_Qkappa(ir)),sngl(mars_1D_Qmu(ir))
      enddo
      print *
    endif
  endif

  ! checks radius is increasing
  ! note that table entries have reversed order compared to file mars_1D.dat
  do ir = 2,NR_mars_1D_layers
    dr = mars_1D_radius(ir) - mars_1D_radius(ir-1)
    ! slightly negative steps should not occur
    if (dr < - TOL) then
      print *,'Error layering inversed: radius should decrease in file mars_1D.dat and start from top down to center.'
      print *,'  wrong layering: ',ir,mars_1D_radius(ir),'and ',ir-1,mars_1D_radius(ir-1)
      stop 'Error layering in mars_1D.dat has wrong ordering'
    endif
  enddo

  ! checks data table
  ! currently, interface depths in file mars_1D.dat must match with the parameters defined in the module model_mars_1d_par above
  ! as those values are used for the meshing and definition of layering in get_model_parameters.F90, e.g.,
  ! in routine get_model_parameters_radii()
  !
  ! ICB at radius(33)
  if (abs(RICB - mars_1D_radius(interface_index_icb)) > TOL) &
    stop 'Error: mars_1D radius RICB and model radius index not matching'
  if (abs(RCMB - mars_1D_radius(interface_index_cmb)) > TOL) &
    stop 'Error: mars_1D radius RCMB and model radius index not matching'
  if (abs(ROCEAN - mars_1D_radius(NR_mars_1D_layers)) > TOL) &
    stop 'Error: mars_1D radius ROCEAN and model radius index not matching'

  ! safety checks if surface is consistent with R_PLANET setting and mars 1D parameter from module
  r_top = mars_1D_radius(nlines)      ! arrays starts from bottom to top
  if (abs(MARS_1D_RSURFACE - r_top) > TOL) then
    print *,'Model error: mars_1D.dat model has a different surface radius than defined in module mars_1D_par:'
    print *,'             table: ',r_top/1000.d0,' versus module value: ',MARS_1D_RSURFACE/1000.d0
    stop 'Error: mars_1D.dat has invalid module surface radius'
  endif
  if (abs(R_PLANET - r_top) > TOL) then
    print *,'Model error: mars_1D.dat model has a different surface radius than used for mars planet:'
    print *,'             table: ',r_top/1000.d0,' versus mars radius value: ',R_PLANET/1000.d0
    stop 'Error: mars_1D.dat has invalid mars surface radius'
  endif

  ! stores Qmu original values (without CRUSTAL modifications) for further use in attenuation routine
  mars_1D_Qmu_original(:) = mars_1D_Qmu(:)

  ! ONE_CRUST uses a single layer value for the crust
  if (ONE_CRUST) then
    ! mid crust is in (2:3)
    rho_tmp = mars_1D_density(NR_mars_1D_layers-3)
    vp_tmp = mars_1D_vp(NR_mars_1D_layers-3)
    vs_tmp = mars_1D_vs(NR_mars_1D_layers-3)
    Qkappa_tmp = mars_1D_Qkappa(NR_mars_1D_layers-3)
    Qmu_tmp = mars_1D_Qmu(NR_mars_1D_layers-3)
    ! assign all crust (1:6) values to upper crust
    mars_1D_density(NR_mars_1D_layers-5:NR_mars_1D_layers) = rho_tmp
    mars_1D_vp(NR_mars_1D_layers-5:NR_mars_1D_layers) = vp_tmp
    mars_1D_vs(NR_mars_1D_layers-5:NR_mars_1D_layers) = vs_tmp
    mars_1D_Qkappa(NR_mars_1D_layers-5:NR_mars_1D_layers) = Qkappa_tmp
    mars_1D_Qmu(NR_mars_1D_layers-5:NR_mars_1D_layers) = Qmu_tmp
  endif

  ! in case an external crustal model will be superimposed on top, we extend mantle values to the surface
  !
  ! strip the crust and replace it by mantle if we use an external crustal model
  ! note: assumes that the crust is given by 6 layers
  !       see in mars_1D.dat: regolith layer (1:2), upper crust (3:4), lower crust (5:6), the mantle
  if (SUPPRESS_CRUSTAL_MESH .or. USE_EXTERNAL_CRUSTAL_MODEL) then
    do i = NR_mars_1D_layers-5,NR_mars_1D_layers
      mars_1D_density(i) = mars_1D_density(NR_mars_1D_layers-6)
      mars_1D_vp(i) = mars_1D_vp(NR_mars_1D_layers-6)
      mars_1D_vs(i) = mars_1D_vs(NR_mars_1D_layers-6)
      mars_1D_Qkappa(i) = mars_1D_Qkappa(NR_mars_1D_layers-6)
      mars_1D_Qmu(i) = mars_1D_Qmu(NR_mars_1D_layers-6)
    enddo
  endif

  end subroutine define_model_mars_1D

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_mars_1D_density(x,rho)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,CRUSTAL

  use model_mars_1d_par

  implicit none

  double precision,intent(in) :: x
  double precision,intent(out) :: rho

  ! local parameters
  integer :: i
  double precision :: r,dr,frac

  ! first allocates and defines model arrays (will be done only once)
  if (.not. allocated(mars_1D_radius)) then
    ! all processes will define same parameters
    call define_model_mars_1D(CRUSTAL)
  endif

  ! initializes
  rho = 0.d0

  ! compute real physical radius in meters
  r = x * R_PLANET

  i = 1
  do while(r >= mars_1D_radius(i) .and. i /= NR_mars_1D_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_mars_1D_layers) i = NR_mars_1D_layers

  ! note: this routine helps to sample the density profile of the 1D model
  !       we might assign the "upper" value for points on interfaces like ICB, CMB etc.
  !       this effect should not matter much though since the sampling rate is quite high along the profile.

  ! calculates density according to radius
  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    rho = mars_1D_density(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = mars_1D_radius(i) - mars_1D_radius(i-1)
    ! fraction for linear interpolation
    if (abs(dr) > 1.d-9) then
      frac = (r - mars_1D_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = mars_1D_density(i-1) + frac * (mars_1D_density(i)-mars_1D_density(i-1))
  endif

  ! non-dimensionalize
  rho = rho*1000.0d0/RHOAV

  end subroutine model_mars_1D_density

