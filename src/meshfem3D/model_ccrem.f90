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
! CCREM
!
! 1D isotropic Earth model by:
!    Ma, X. and H. Tkalcic (2021),
!    CCREM: New Reference Earth Model From the Global Coda-Correlation Wavefield,
!    JGR Solid Earth, 126, e2021JB022515.
!    https://doi.org/10.1029/2021JB022515
!
!  Model file is in DATA/CCREM:
!    file ccrem.dat contains the model values provided by table S3 in the supporting material.
!
!--------------------------------------------------------------------------------------------------

  module model_ccrem_par

  implicit none

  ! radius at physical surface
  double precision, parameter :: CCREM_RSURFACE = 6371000.d0

  ! CCREM radii
  ! no ocean, same as surface radius (Physical surface)
  double precision, parameter :: CCREM_ROCEAN = CCREM_RSURFACE ! no ocean
  double precision, parameter :: CCREM_RMIDDLE_CRUST = CCREM_RSURFACE - 20000.d0 ! depth = 20 km
  ! Crust-mantle boundary
  double precision, parameter :: CCREM_RMOHO = CCREM_RSURFACE - 35000.d0         ! depth = 35 km
  ! Rheological lithosphere
  double precision, parameter :: CCREM_R80  =  CCREM_RSURFACE - 80000.d00         ! depth = 80 km
  ! Thermal lithosphere
  double precision, parameter :: CCREM_R220 = CCREM_RSURFACE - 220000.d0          ! depth = 220 km
  double precision, parameter :: CCREM_R400 = CCREM_RSURFACE - 410000.d0          ! depth = 410 km - CCREM depth
  double precision, parameter :: CCREM_R600 = CCREM_RSURFACE - 600000.d0          ! depth = 600 km
  ! alpha-olivine-beta-spinel transition
  double precision, parameter :: CCREM_R670 = CCREM_RSURFACE - 660000.d0          ! depth = 660 km - CCREM depth
  ! beta-spinel-gamma-spinel transition
  double precision, parameter :: CCREM_R771 = CCREM_RSURFACE - 771000.d0          ! depth = 771 km (PREM)
  ! lower thermal boundary layer
  double precision, parameter :: CCREM_RTOPDDOUBLEPRIME = CCREM_RSURFACE - 2741000.d0 ! depth = 2741 km (PREM)
  ! Core-mantle boundary
  double precision, parameter :: CCREM_RCMB = CCREM_RSURFACE - 2891000.d0         ! depth = 2891 km (PREM)
  ! inner core radius
  double precision, parameter :: CCREM_RICB = CCREM_RSURFACE - 5153500.d0         ! depth = 5153.5 km (ak135)

  ! densities in [kg/m^3]: rho [kg/m^3] = rho * 1000 [g/cm^3]
  double precision, parameter :: CCREM_RHO_OCEANS    = 1020.d0              ! will not be used
  double precision, parameter :: CCREM_RHO_TOP_OC    = 9.9131d0 * 1000.d0   ! outer core density fluid outer core (from PREM)
  double precision, parameter :: CCREM_RHO_BOTTOM_OC = 12.1478d0 * 1000.d0  ! outer core density

  ! number of layers
  integer, parameter :: NR_CCREM_layers = 222  ! check with ccrem.dat in DATA/CCREM/ folder

  ! model arrays
  double precision, dimension(:),allocatable :: CCREM_radius,CCREM_density, &
    CCREM_vp,CCREM_vs,CCREM_Qkappa,CCREM_Qmu

  ! storage of initial Qmu profile (used for attenuation setup)
  double precision, dimension(:),allocatable :: CCREM_Qmu_original

  end module model_ccrem_par


!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_ccrem_broadcast(CRUSTAL)

! standard routine to setup model

  use constants, only: myrank
  use model_ccrem_par

  implicit none

  logical,intent(in) :: CRUSTAL

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(CCREM_radius(NR_CCREM_layers), &
           CCREM_density(NR_CCREM_layers), &
           CCREM_vp(NR_CCREM_layers), &
           CCREM_vs(NR_CCREM_layers), &
           CCREM_Qkappa(NR_CCREM_layers), &
           CCREM_Qmu(NR_CCREM_layers), &
           CCREM_Qmu_original(NR_CCREM_layers), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating CCREM arrays'
  CCREM_radius(:) = 0.d0; CCREM_density(:) = 0.d0
  CCREM_vp(:) = 0.d0; CCREM_vs(:) = 0.d0
  CCREM_Qkappa(:) = 0.d0; CCREM_Qmu(:) = 0.d0; CCREM_Qmu_original(:) = 0.d0

  ! main process will read in parameters
  if (myrank == 0) call define_model_ccrem(CRUSTAL)

  ! broadcast the information read on the main to the nodes
  call bcast_all_dp(CCREM_radius,NR_CCREM_layers)
  call bcast_all_dp(CCREM_density,NR_CCREM_layers)
  call bcast_all_dp(CCREM_vp,NR_CCREM_layers)
  call bcast_all_dp(CCREM_vs,NR_CCREM_layers)
  call bcast_all_dp(CCREM_Qkappa,NR_CCREM_layers)
  call bcast_all_dp(CCREM_Qmu,NR_CCREM_layers)
  call bcast_all_dp(CCREM_Qmu_original,NR_CCREM_layers)

  end subroutine model_ccrem_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_ccrem(x,rho,vp,vs,Qkappa,Qmu,idoubling,iregion_code)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV

  use model_ccrem_par

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
  double precision,intent(inout) :: rho
  double precision,intent(inout) :: vp,vs,Qmu,Qkappa
  integer,intent(in) :: idoubling
  integer,intent(in) :: iregion_code

  ! local parameters
  double precision :: r,frac,scaleval,dr
  integer :: i

  ! checks if model arrays are defined
  if (.not. allocated(CCREM_radius)) then
    stop 'Error in model_ccrem(), arrays not defined yet'
  endif

  ! compute real physical radius in meters
  r = x * R_PLANET

  ! check flags to make sure we correctly honor the discontinuities
  ! we use strict inequalities since r has been slightly changed in mesher

  ! note: using stop statements, not exit_mpi() calls to avoid the need for MPI libraries when linking xcreate_header_file

  ! checks doubling flag
  !
  !--- inner core
  !
  if (r >= 0.d0 .and. r < CCREM_RICB) then
    ! checks with inner core flags
    if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
        idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
        idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
        idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
        idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
         stop 'wrong doubling flag for inner core point in model_ccrem()'
  !
  !--- outer core
  !
  else if (r > CCREM_RICB .and. r < CCREM_RCMB) then
    if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
      stop 'wrong doubling flag for outer core point in model_ccrem()'
  !
  !--- D" at the base of the mantle
  !
  else if (r > CCREM_RCMB .and. r < CCREM_RTOPDDOUBLEPRIME) then
    if (idoubling /= IFLAG_MANTLE_NORMAL) &
      stop 'wrong doubling flag for D" point in model_ccrem()'
  !
  !--- mantle: from top of D" to d670
  !
  else if (r > CCREM_RTOPDDOUBLEPRIME .and. r < CCREM_R670) then
    if (idoubling /= IFLAG_MANTLE_NORMAL) &
      stop 'wrong doubling flag for top D" - > d670 point in model_ccrem()'
  !
  !--- mantle: from d670 to d220
  !
  else if (r > CCREM_R670 .and. r < CCREM_R220) then
    if (idoubling /= IFLAG_670_220) &
      stop 'wrong doubling flag for d670 - > d220 point in model_ccrem()'
  !
  !--- mantle and crust: from d220 to MOHO and then to surface
  !
  else if (r > CCREM_R220) then
    if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
      stop 'wrong doubling flag for d220 - > Moho - > surface point in model_ccrem()'
  endif

  ! searches radius index
  i = 1
  do while(r >= CCREM_radius(i) .and. i /= NR_CCREM_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_CCREM_layers) i = NR_CCREM_layers

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (iregion_code /= 0) then
    ! check
    if (iregion_code /= IREGION_INNER_CORE .and. &
        iregion_code /= IREGION_OUTER_CORE .and. &
        iregion_code /= IREGION_CRUST_MANTLE) stop 'Invalid iregion_code in model_ccrem() routine'

    ! see table in ccrem.dat with index reversed, i.e., starting with index 1 at the center (radius == 0)
    ! inner core range [1:41]
    if (iregion_code == IREGION_INNER_CORE .and. i > 41) i = 41

    ! outer core range [42:118]
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_OUTER_CORE .and. i <= 42) i = 42 + 1
    if (iregion_code == IREGION_OUTER_CORE .and. i > 118) i = 118

    ! crust/mantle range [119:222] NR_CCREM_layers == 222
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_CRUST_MANTLE .and. i <= 119) i = 119 + 1
    if (iregion_code == IREGION_CRUST_MANTLE .and. i > NR_CCREM_layers) i = NR_CCREM_layers
  endif

  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    rho = CCREM_density(i)
    vp = CCREM_vp(i)
    vs = CCREM_vs(i)
    Qmu = CCREM_Qmu(i)
    Qkappa = CCREM_Qkappa(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = CCREM_radius(i) - CCREM_radius(i-1)
    if (abs(dr) > 1.d-9) then
      frac = (r - CCREM_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = CCREM_density(i-1) + frac * (CCREM_density(i)-CCREM_density(i-1))
    vp = CCREM_vp(i-1) + frac * (CCREM_vp(i)-CCREM_vp(i-1))
    vs = CCREM_vs(i-1) + frac * (CCREM_vs(i)-CCREM_vs(i-1))
    Qmu = CCREM_Qmu(i-1) + frac * (CCREM_Qmu(i)-CCREM_Qmu(i-1))
    Qkappa = CCREM_Qkappa(i-1) + frac * (CCREM_Qkappa(i)-CCREM_Qkappa(i-1))
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
  rho = rho*1000.0d0/RHOAV              ! [g/cm^3] - > [kg/m^3]: rho [kg/m^3] = rho * 1000 [g/cm^3]
  vp = vp*1000.0d0/(R_PLANET*scaleval)  ! [km/s] - > [m/s]: v [m/s] = v * 1000 [km/s]
  vs = vs*1000.0d0/(R_PLANET*scaleval)

  end subroutine model_ccrem

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_model_ccrem(CRUSTAL)

  use constants, only: myrank,IIN,SUPPRESS_CRUSTAL_MESH,ATTENUATION_COMP_MAXIMUM
  use shared_parameters, only: RICB,RCMB,ROCEAN,ONE_CRUST

  use model_ccrem_par

  implicit none

  logical,intent(in) :: CRUSTAL

  ! local parameters
  integer :: nlines,ir,ier
  double precision :: dep,r,vp,vs,rho
  double precision :: dr,rmin,rmax,r_prem,Qmu_tmp,Qkappa_tmp,rho_tmp,vp_tmp,vs_tmp
  double precision :: Qkappa,Qmu
  character(len=256) :: datafile,line

  ! tolerance for checking values
  double precision, parameter :: TOL = 1.d-9

  ! debugging output
  logical, parameter :: DEBUG = .false.

  ! note: this routine is only called by the main process (myrank == 0)

  ! file name
  datafile = 'DATA/CCREM/ccrem.dat'

  ! reads in data file
  open(IIN,file=trim(datafile),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(datafile)
    stop 'Error opening ccrem.dat file'
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
  if (nlines /= NR_CCREM_layers) then
    print *,'Error: number of data entries ',nlines,' should match number of layers ',NR_CCREM_layers,' for ccrem model'
    print *,'Please check file: ',trim(datafile)
    stop 'Invalid number of data line entries in ccrem.dat'
  endif

  ! read in layering
  ! note: table starts from top surface to inner core
  !       for the model arrays here, we choose an opposite order, with index 1 being the center of Earth
  !       and then increasing radius.
  ier = 0
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
    ! format: #depth(km) # radius(km) #vp(km/s) #vs(km/s) #density(g/cm^3)
    read(line,*,iostat=ier) dep,r,vp,vs,rho

    if (ier /= 0) then
      print *,'Error: reading ccrem.dat data, on line ',nlines
      print *,'       data line: ',trim(line)
      print *,'Please check data file ccrem.dat'
      stop 'Error reading data in ccrem data'
    endif

    ! checks values
    if (r < 0.d0) stop 'Error found negative r value'
    if (vs < 0.d0) stop 'Error found negative vs value'
    if (vp <= 0.d0) stop 'Error found negative or zero vp value'

    if (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2) <= 0.d0) then
      print *,'Error: invalid vp,vs = ',vp,vs,' leads to a rock type with negative factor ',(1.d0 - 4.d0/3.d0*(vs**2)/(vp**2))
      print *,'       data line: ',trim(line)
      stop 'Error invalid vp,vs values'
    endif

    ! increases layer number count
    nlines = nlines + 1

    ! define all the values in the model
    ! change ordering from core to top
    ir = NR_CCREM_layers - (nlines-1)
    if (ir < 1) stop 'Error reading in ccrem data lines'

    ! we need radius in [m], density in [g/cm^3], vp,vs in [km/s] (and Qkappa and Qmu from PREM)
    CCREM_radius(ir) = r * 1000.d0
    CCREM_density(ir) = rho
    CCREM_vp(ir) = vp
    CCREM_vs(ir) = vs

    ! debug
    !print *,'debug: CCREM: ',ir,CCREM_radius(ir),CCREM_vp(ir),CCREM_vs(ir),CCREM_density(ir)
  enddo
  close(IIN)

  ! checks radius is increasing
  do ir = 2,NR_CCREM_layers
    dr = CCREM_radius(ir) - CCREM_radius(ir-1)
    ! slightly negative steps should not occur
    if (dr < - TOL) then
      print *,'Error layering inversed: radius should decrease in file ccrem.dat and start from top down to center.'
      print *,'  wrong layering: ',ir,CCREM_radius(ir),'and ',ir-1,CCREM_radius(ir-1)
      stop 'Error layering in ccrem.dat has wrong ordering'
    endif
  enddo

  ! checks data table
  if (abs(RICB - CCREM_radius(42)) > TOL) &
    stop 'Error: CCREM radius RICB and model radius index not matching'
  if (abs(RCMB - CCREM_radius(119)) > TOL) &
    stop 'Error: CCREM radius RCMB and model radius index not matching'
  if (abs(ROCEAN - CCREM_radius(NR_CCREM_layers)) > TOL) &
    stop 'Error: CCREM radius ROCEAN and model radius index not matching'

  ! sets attenuation model values taken from PREM
  do ir = 1,NR_CCREM_layers
    r = CCREM_radius(ir)
    ! make sure we are within the right shell in CCREM to honor discontinuities
    if (ir <= 42) then
      ! inner core
      rmin = 0.d0
      rmax = CCREM_RICB
    else if (ir <= 119) then
      ! outer core
      rmin = CCREM_RICB
      rmax = CCREM_RCMB
    else if (ir <= 195) then
      ! lower mantle
      rmin = CCREM_RCMB
      rmax = CCREM_R670
    else if (ir <= 205) then
      ! 660km - 410km
      rmin = CCREM_R670
      rmax = CCREM_R400
    else if (ir <= 219) then
      ! 410km - 35km (no discontinuities between 410 and 35km depth in CCREM)
      rmin = CCREM_R400
      rmax = CCREM_RMOHO
    else
      ! 35km - top
      rmin = CCREM_RMOHO
      rmax = CCREM_RSURFACE
    endif

    ! use small geometrical tolerance
    r_prem = r
    if (r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
    if (r >= rmax*0.999999d0) r_prem = rmax*0.999999d0
    r = r_prem

    ! determines attenuation Q factors from corresponding PREM shell values
    Qmu = ATTENUATION_COMP_MAXIMUM
    Qkappa = ATTENUATION_COMP_MAXIMUM

    ! inner core
    if (r >= 0.d0 .and. r <= CCREM_RICB) then
      Qmu = 84.6d0
      Qkappa = 1327.7d0
    ! outer core
    else if (r > CCREM_RICB .and. r <= CCREM_RCMB) then
      Qmu = 0.0d0
      Qkappa = 57827.0d0
    ! D" at the base of the mantle
    else if (r > CCREM_RCMB .and. r <= CCREM_RTOPDDOUBLEPRIME) then
      Qmu = 312.0d0
      Qkappa = 57827.0d0
    ! mantle: from top of D" to d670
    else if (r > CCREM_RTOPDDOUBLEPRIME .and. r <= CCREM_R771) then
      Qmu = 312.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_R771 .and. r <= CCREM_R670) then
      Qmu = 312.0d0
      Qkappa = 57827.0d0
    ! mantle: above d670
    else if (r > CCREM_R670 .and. r <= CCREM_R600) then
      Qmu = 143.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_R600 .and. r <= CCREM_R400) then
      Qmu = 143.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_R400 .and. r <= CCREM_R220) then
      Qmu = 143.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_R220 .and. r <= CCREM_R80) then
      Qmu = 80.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_R80 .and. r <= CCREM_RMOHO) then
      Qmu = 600.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_RMOHO .and. r <= CCREM_RMIDDLE_CRUST) then
      Qmu = 600.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_RMIDDLE_CRUST .and. r <= CCREM_ROCEAN) then
      Qmu = 600.0d0
      Qkappa = 57827.0d0
    else if (r > CCREM_ROCEAN) then
      Qmu = 600.0d0
      Qkappa = 57827.0d0
    endif

    ! stores
    CCREM_Qkappa(ir) = Qkappa
    CCREM_Qmu(ir) = Qmu
  enddo

  ! stores Qmu original values (without CRUSTAL modifications) for further use in attenuation routine
  CCREM_Qmu_original(:) = CCREM_Qmu(:)

  ! ONE_CRUST uses a single layer value for the crust
  if (ONE_CRUST) then
    ! crust is in (1:3)
    ! takes upper most crustal value (surface values)
    rho_tmp = CCREM_density(NR_CCREM_layers)
    vp_tmp = CCREM_vp(NR_CCREM_layers)
    vs_tmp = CCREM_vs(NR_CCREM_layers)
    Qkappa_tmp = CCREM_Qkappa(NR_CCREM_layers)
    Qmu_tmp = CCREM_Qmu(NR_CCREM_layers)
    ! assign all crust (1:3) values to upper crust
    CCREM_density(NR_CCREM_layers-2:NR_CCREM_layers) = rho_tmp
    CCREM_vp(NR_CCREM_layers-2:NR_CCREM_layers) = vp_tmp
    CCREM_vs(NR_CCREM_layers-2:NR_CCREM_layers) = vs_tmp
    CCREM_Qkappa(NR_CCREM_layers-2:NR_CCREM_layers) = Qkappa_tmp
    CCREM_Qmu(NR_CCREM_layers-2:NR_CCREM_layers) = Qmu_tmp
  endif

  ! in case an external crustal model will be superimposed on top, we extend mantle values to the surface
  !
  ! strip the crust and replace it by mantle if we use an external crustal model
  ! note: assumes that the crust is given by 3 layers
  !       see in ccrem.dat: upper crust (1:1), middle crust (2:2), lower crust (3:3), then mantle
  if (SUPPRESS_CRUSTAL_MESH .or. CRUSTAL) then
    do ir = NR_CCREM_layers-2,NR_CCREM_layers
      CCREM_density(ir) = CCREM_density(NR_CCREM_layers-3)
      CCREM_vp(ir) = CCREM_vp(NR_CCREM_layers-3)
      CCREM_vs(ir) = CCREM_vs(NR_CCREM_layers-3)
      CCREM_Qkappa(ir) = CCREM_Qkappa(NR_CCREM_layers-3)
      CCREM_Qmu(ir) = CCREM_Qmu(NR_CCREM_layers-3)
    enddo
  endif

  ! debug
  if (DEBUG) then
    if (myrank == 0) then
      print *,'debug: model CCREM'
      print *,'format: #line #r #vp #vs #rho #Qkappa #Qmu'
      do ir = 1,NR_CCREM_layers
        print *,ir,sngl(CCREM_radius(ir)),sngl(CCREM_vp(ir)),sngl(CCREM_vs(ir)),sngl(CCREM_density(ir)), &
                   sngl(CCREM_Qkappa(ir)),sngl(CCREM_Qmu(ir))
      enddo
      print *
    endif
  endif

  end subroutine define_model_ccrem

!
!-------------------------------------------------------------------------------------------------
!

! not needed yet...
!
! would be in case we want to setup gravity and ellipticity with the density profile of this reference model;
! for now, gravity & ellipticity are based on PREM's density profile.
!
!  subroutine model_ccrem_density(x,rho,ONE_CRUST)
!
!  use constants
!  use shared_parameters, only: R_PLANET,RHOAV
!
!  use model_ccrem_par
!
!  implicit none
!
!  double precision,intent(in) :: x
!  double precision,intent(out) :: rho
!  logical,intent(in) :: ONE_CRUST
!
!  ! local parameters
!  integer :: i
!  double precision :: r,dr,frac
!  double precision :: vp_tmp,vs_tmp,rho_tmp,Qmu_tmp,Qkappa_tmp
!
!  ! makes sure arrays are allocated
!  if (.not. allocated(CCREM_radius)) then
!    stop 'Error in model_ccrem_density(), density & radius arrays not defined yet'
!  endif
!
!  ! initializes
!  rho = 0.d0
!
!  ! compute real physical radius in meters
!  r = x * R_PLANET
!
!  i = 1
!  do while(r >= CCREM_radius(i) .and. i /= NR_CCREM_layers)
!    i = i + 1
!  enddo
!
!  ! just double check
!  if (i > NR_CCREM_layers) i = NR_CCREM_layers
!
!  ! note: this routine helps to sample the density profile of the 1D model
!  !       we might assign the "upper" value for points on interfaces like ICB, CMB etc.
!  !       this effect should not matter much though since the sampling rate is quite high along the profile.
!
!  ! calculates density according to radius
!  ! interpolates model values
!  if (i == 1) then
!    ! no interpolation with previous layer possible
!    rho = CCREM_density(i)
!  else
!    ! interpolate from radius(i-1) to r using the values at i-1 and i
!    dr = CCREM_radius(i) - CCREM_radius(i-1)
!    ! fraction for linear interpolation
!    if (abs(dr) > 1.d-9) then
!      frac = (r - CCREM_radius(i-1))/dr
!    else
!      ! both current and lower layer have same radius
!      frac = 1.d0 ! will take values from current layer i
!    endif
!    rho = CCREM_density(i-1) + frac * (CCREM_density(i)-CCREM_density(i-1))
!  endif
!
!  ! non-dimensionalize
!  rho = rho*1000.0d0/RHOAV
!
!  end subroutine model_ccrem_density
!
