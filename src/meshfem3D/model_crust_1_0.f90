!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

!--------------------------------------------------------------------------------------------------
! CRUST 1.0 model
! see: http://igppweb.ucsd.edu/~gabi/crust1.html
!
! Initial release:
! ===============
! 15 July 2013: this is the initial release of the essential files. As
! described on the website, http://igppweb.ucsd.edu/~gabi/crust1.html,
! the structure in the crystalline crust is defined using a crustal type
! assignment. The crustal type file can be obtained upon request.
!
! The 8 crustal layers:
! ====================
! 1) water
! 2) ice
! 3) upper sediments   (VP, VS, rho not defined in all cells)
! 4) middle sediments  "
! 5) lower sediments   "
! 6) upper crystalline crust
! 7) middle crystalline crust
! 8) lower crystalline crust
! + a ninth layer gives V_Pn, V_Sn and rho below the Moho. The values
!   are associated with LLNL model G3Cv3 on continents and a thermal
!   model in the oceans.
!
! reads and smooths crust1.0 model
!--------------------------------------------------------------------------------------------------

  module model_crust_1_0_par

  ! crustal_model_constants
  ! crustal model parameters for crust1.0
  integer, parameter :: CRUST1_NP  = 9
  integer, parameter :: CRUST1_NLO = 360
  integer, parameter :: CRUST1_NLA = 180

  ! cap smoothing
  integer, parameter :: NCAP_CRUST1 = 180

  ! model_crust_variables
  ! Vp, Vs and density
  double precision, dimension(:,:,:),allocatable :: crust1_vp,crust1_vs,crust1_rho

  ! layer thickness
  double precision, dimension(:,:,:),allocatable :: crust1_thk

  end module model_crust_1_0_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_crust_1_0_broadcast(myrank)

! standard routine to setup model

  use constants
  use model_crust_1_0_par

  implicit none

  integer :: myrank
  integer :: ier

  ! allocate crustal arrays
  allocate(crust1_vp(CRUST1_NP,CRUST1_NLA,CRUST1_NLO), &
           crust1_vs(CRUST1_NP,CRUST1_NLA,CRUST1_NLO), &
           crust1_rho(CRUST1_NP,CRUST1_NLA,CRUST1_NLO), &
           crust1_thk(CRUST1_NP,CRUST1_NLA,CRUST1_NLO), &
           stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating crustal arrays')

  ! initializes
  crust1_vp(:,:,:) = ZERO
  crust1_vs(:,:,:) = ZERO
  crust1_rho(:,:,:) = ZERO
  crust1_thk(:,:,:) = ZERO

  ! the variables read are declared and stored in structure model_crust_1_0_par
  if(myrank == 0) call read_crust_1_0_model()

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(crust1_vp,CRUST1_NP*CRUST1_NLA*CRUST1_NLO)
  call bcast_all_dp(crust1_vs,CRUST1_NP*CRUST1_NLA*CRUST1_NLO)
  call bcast_all_dp(crust1_rho,CRUST1_NP*CRUST1_NLA*CRUST1_NLO)
  call bcast_all_dp(crust1_thk,CRUST1_NP*CRUST1_NLA*CRUST1_NLO)

  end subroutine model_crust_1_0_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_crust_1_0(lat,lon,x,vp,vs,rho,moho,found_crust,elem_in_crust)

  use constants
  use model_crust_1_0_par

  implicit none

  double precision,intent(in) :: lat,lon,x
  double precision,intent(out) :: vp,vs,rho,moho
  logical,intent(out) :: found_crust
  logical,intent(in):: elem_in_crust

  ! local parameters
  double precision :: h_sed,h_uc
  double precision :: x3,x4,x5,x6,x7,x8
  double precision :: scaleval
  double precision,dimension(CRUST1_NP):: vps,vss,rhos,thicks

  ! initializes
  vp = ZERO
  vs = ZERO
  rho = ZERO
  moho = ZERO

  ! gets smoothed structure
  call crust_1_0_CAPsmoothed(lat,lon,vps,vss,rhos,thicks)

  ! note: we ignore water & ice sheets
  ! (only elastic layers are considered)

  ! whole sediment thickness
  h_sed = thicks(3) + thicks(4) + thicks(5)

  ! upper crust thickness (including sediments above)
  h_uc = h_sed + thicks(6)

  ! non-dimensionalization factor
  scaleval = ONE / R_EARTH_KM

  ! non-dimensionalizes thickness (given in km)
  ! upper sediment
  x3 = ONE - thicks(3) * scaleval
  ! middle sediment
  x4 = ONE - (thicks(3) + thicks(4)) * scaleval
  ! all sediments
  x5 = ONE - h_sed * scaleval
  ! upper crust
  x6 = ONE - h_uc * scaleval
  ! middle crust
  x7 = ONE - (h_uc+thicks(7)) * scaleval
  ! lower crust
  x8 = ONE - (h_uc+thicks(7)+thicks(8)) * scaleval

  ! checks moho value
  !moho = h_uc + thicks(6) + thicks(7)
  !if( moho /= thicks(NLAYERS_CRUST) ) then
  ! print*,'moho:',moho,thicks(NLAYERS_CRUST)
  ! print*,'  lat/lon/x:',lat,lon,x
  !endif

  ! No matter found_crust true or false, output moho thickness
  moho = (h_uc+thicks(7)+thicks(8)) * scaleval

  ! initializes
  found_crust = .true.

  ! gets corresponding crustal velocities and density
  if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST ) then
    vp = vps(3)
    vs = vss(3)
    rho = rhos(3)
  else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST ) then
    vp = vps(4)
    vs = vss(4)
    rho = rhos(4)
  else if(x > x5 .and. INCLUDE_SEDIMENTS_CRUST ) then
    vp = vps(5)
    vs = vss(5)
    rho = rhos(5)
  else if(x > x6) then
    vp = vps(6)
    vs = vss(6)
    rho = rhos(6)
  else if(x > x7) then
    vp = vps(7)
    vs = vss(7)
    rho = rhos(7)
  else if(x > x8 .or. elem_in_crust) then
    ! takes lower crustal values only if x is slightly above moho depth or
    ! if elem_in_crust is set
    !
    ! note: it looks like this does distinguish between GLL points at the exact moho boundary,
    !          where the point is on the interface between both,
    !          oceanic elements and mantle elements below
    vp = vps(8)
    vs = vss(8)
    rho = rhos(8)
  else
    ! note: if x is exactly the moho depth this will return false
    found_crust = .false.
  endif

  ! non-dimensionalize
  if (found_crust) then
    scaleval = ONE / ( R_EARTH_KM * dsqrt(PI*GRAV*RHOAV) )
    vp = vp * scaleval
    vs = vs * scaleval
    rho = rho * 1000.0d0 / RHOAV
 endif

 end subroutine model_crust_1_0

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_crust_1_0_model()

  use constants
  use model_crust_1_0_par

  implicit none

  ! local variables
  integer :: ier
  integer :: i,j,k
  ! boundaries
  double precision, dimension(:,:,:),allocatable :: bnd
  ! crustal / sediment thickness
  double precision, dimension(:,:),allocatable :: thc,ths
  double precision :: lat,lon,x
  double precision :: vp,vs,rho,moho
  double precision :: h_moho_min,h_moho_max
  logical :: found_crust

  !-------------------------------------------------------------
  ! debugging user parameter
  ! outputs files for inspection
  logical,parameter :: DEBUG_FILE_OUTPUT = .false.
  !-------------------------------------------------------------

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: CRUST1.0'
  write(IMAIN,*)

  ! allocates temporary array
  allocate(bnd(CRUST1_NP,CRUST1_NLA,CRUST1_NLO), &
           stat=ier)
  if( ier /= 0 ) call exit_MPI(0,'error allocating crustal arrays in read routine')

  ! initializes
  bnd(:,:,:) = ZERO

  ! opens crust1.0 data files
  open(51,file='DATA/crust1.0/crust1.vp',action='read',status='old',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "DATA/crust1.0/crust1.vp": ', ier
    call exit_MPI(0,'error model crust1.0: file not found DATA/crust1.0/crust1.vp')
  endif

  open(52,file='DATA/crust1.0/crust1.vs',action='read',status='old',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "DATA/crust1.0/crust1.vs": ', ier
    call exit_MPI(0,'error model crust1.0: file not found DATA/crust1.0/crust1.vs')
  endif

  open(53,file='DATA/crust1.0/crust1.rho',action='read',status='old',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "DATA/crust1.0/crust1.rho": ', ier
    call exit_MPI(0,'error model crust1.0: file not found DATA/crust1.0/crust1.rho')
  endif

  open(54,file='DATA/crust1.0/crust1.bnds',action='read',status='old',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "DATA/crust1.0/crust1.bnds": ', ier
    call exit_MPI(0,'error model crust1.0: file not found DATA/crust1.0/crust1.bnds')
  endif

  ! reads in data values
  do j = 1,CRUST1_NLA
    do i = 1,CRUST1_NLO
      read(51,*)(crust1_vp(k,j,i),k=1,CRUST1_NP)
      read(52,*)(crust1_vs(k,j,i),k=1,CRUST1_NP)
      read(53,*)(crust1_rho(k,j,i),k=1,CRUST1_NP)
      read(54,*)(bnd(k,j,i),k=1,CRUST1_NP)
    enddo
  enddo

  ! closes files
  close(51)
  close(52)
  close(53)
  close(54)

  ! determines layer thickness
  do j = 1,CRUST1_NLA
    do i = 1,CRUST1_NLO
      do k = 1,CRUST1_NP - 1
        crust1_thk(k,j,i) = - (bnd(k+1,j,i) - bnd(k,j,i))
      enddo
    enddo
  enddo

  ! frees memory
  deallocate(bnd)

  ! additional info
  if( DEBUG_FILE_OUTPUT ) then
    ! allocates temporary arrays
    allocate(thc(CRUST1_NLA,CRUST1_NLO), &
             ths(CRUST1_NLA,CRUST1_NLO), &
             stat=ier)
    if( ier /= 0 ) call exit_MPI(0,'error allocating crustal arrays in read routine')

    thc(:,:) = ZERO
    ths(:,:) = ZERO

    ! debug: file output for original data
    open(77,file='tmp-crust1.0.dat',status='unknown')
    write(77,*)'#crustal thickness: #lat (degree) #lon (degree) #crust (km) (including ice) #crust (w/out ice) #sediment #ice'

    h_moho_min = HUGEVAL
    h_moho_max = -HUGEVAL

    ! crustal thickness
    ! thickness = ice (layer index 2) + sediment (index 3+4+5) + crystalline crust (index 6+7+8)
    do j = 1,CRUST1_NLA
      do i = 1,CRUST1_NLO
        ! crustal thickness with ice
        !thc(j,i) = crust1_thk(2,j,i) &
        !         + crust1_thk(3,j,i) + crust1_thk(4,j,i) + crust1_thk(5,j,i)  &
        !         + crust1_thk(6,j,i) + crust1_thk(7,j,i) + crust1_thk(8,j,i)

        ! sediment thickness
        ths(j,i) = crust1_thk(3,j,i) + crust1_thk(4,j,i) + crust1_thk(5,j,i)


        ! crustal thickness without ice
        ! note: etopo1 has topography including ice ("ice surface" version) and at base of ice sheets ("bedrock" version)
        !       see: http://www.ngdc.noaa.gov/mgg/global/global.html
        thc(j,i) = crust1_thk(3,j,i) + crust1_thk(4,j,i) + crust1_thk(5,j,i)  &
                 + crust1_thk(6,j,i) + crust1_thk(7,j,i) + crust1_thk(8,j,i)

        ! limit moho thickness
        if(thc(j,i) > h_moho_max) h_moho_max = thc(j,i)
        if(thc(j,i) < h_moho_min) h_moho_min = thc(j,i)

        write(77,*)(90.0-j+0.5),(-180.0+i-0.5),thc(j,i)+crust1_thk(2,j,i),thc(j,i),ths(j,i),crust1_thk(2,j,i)
      enddo
    enddo
    close(77)

    ! checks min/max
    if(h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) &
      stop 'incorrect moho depths in read_crust_1_0_model'

    ! debug: file output for smoothed data
    open(77,file='tmp-crust1.0-smooth.dat',status='unknown')
    write(77,*)'#crustal thickness: #lat (degree) #lon (degree) #moho (km) (w/out ice) #vp (at surface) #vs (at surface)'

    h_moho_min = HUGEVAL
    h_moho_max = -HUGEVAL

    ! smoothed version
    do j = 1,CRUST1_NLA
      lat = 90.d0 - j + 0.5
      do i = 1,CRUST1_NLO
        lon = -180.d0 + i - 0.5
        x = 1.0d0
        call model_crust_1_0(lat,lon,x,vp,vs,rho,moho,found_crust,.false.)

        ! limit moho thickness
        if(moho > h_moho_max) h_moho_max = moho
        if(moho < h_moho_min) h_moho_min = moho

        write(77,*)lat,lon,moho*R_EARTH_KM, &
         vp*(R_EARTH_KM*dsqrt(PI*GRAV*RHOAV)),vs*(R_EARTH_KM*dsqrt(PI*GRAV*RHOAV)),rho*(RHOAV/1000.0d0)
      enddo
    enddo
    close(77)

    ! checks min/max
    if(h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) &
      stop 'incorrect moho depths in read_crust_1_0_model'

    ! frees memory
    deallocate(ths,thc)
  endif

  end subroutine read_crust_1_0_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crust_1_0_CAPsmoothed(lat,lon,velp,vels,rho,thick)

! crustal vp and vs in km/s, layer thickness in km
!
! crust1.0 gets smoothed with a cap of size CAP using NTHETA points
! in the theta direction and NPHI in the phi direction.
! The cap is rotated to the North Pole.

  use constants
  use model_crust_1_0_par

  implicit none

  ! sampling rate for CAP points
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 20

  ! argument variables
  double precision :: lat,lon
  double precision,dimension(CRUST1_NP) :: rho,thick,velp,vels


  !-------------------------------
  ! CAP smoothing - extension range in degree
  ! note: using a smaller range, e.g. 0.5 degrees, leads to undefined Jacobian error at different places.
  !       this is probably due to stretching elements below sharp gradients, especially with deep moho values.
  !       so far, the only thing that works is to smooth out values and take special care of the Andes...
  ! TODO: one could try to adapt this degree range to the simulation resolution in future
  double precision,parameter :: CAP_DEFAULT_DEGREE = 1.0d0

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
  double precision :: rhol(CRUST1_NP),thickl(CRUST1_NP),velpl(CRUST1_NP),velsl(CRUST1_NP)

  double precision :: weightl,cap_degree
  double precision :: dist
  double precision :: h_sed
  integer :: i,icolat,ilon

  ! checks latitude/longitude
  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) then
    print*,'error in lat/lon:',lat,lon
    stop 'error in latitude/longitude range in crust1.0'
  endif

  ! makes sure lat/lon are within range
  if(lat==90.0d0) lat=89.9999d0
  if(lat==-90.0d0) lat=-89.9999d0
  if(lon==180.0d0) lon=179.9999d0
  if(lon==-180.0d0) lon=-179.9999d0

  ! sets up smoothing points
  ! by default uses CAP smoothing with 1 degree
  cap_degree = CAP_DEFAULT_DEGREE

  ! checks if inside/outside of critical region for mesh stretching
  if( SMOOTH_CRUST ) then
    dist = dsqrt( (lon-LON_CRITICAL_ANDES)**2 + (lat-LAT_CRITICAL_ANDES )**2 )
    if( dist < CRITICAL_RANGE ) then
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
  ! (see routine in model_crust.f90)
  call CAP_vardegree(lon,lat,xlon,xlat,weight,cap_degree,NTHETA,NPHI)

  ! initializes
  velp(:) = ZERO
  vels(:) = ZERO
  rho(:) = ZERO
  thick(:) = ZERO

  ! loops over weight points
  do i=1,NTHETA*NPHI
    ! gets lat/lon indices

    ! checks latitude/longitude value
    if(xlat(i) > 90.0d0 .or. xlat(i) < -90.0d0 .or. xlon(i) > 180.0d0 .or. xlon(i) < -180.0d0) then
      print*,'error in lat/lon range:',xlat(i),xlon(i)
      stop 'error in latitude/longitude range in crust1.0'
    endif

    icolat = int(1 + (90.d0-xlat(i)) )
    if(icolat == 181) icolat = 180
    ! checks
    if(icolat>180 .or. icolat<1) then
      print*,'error in lat/lon range: icolat = ',icolat
      stop 'error in routine icolat/ilon crust1.0'
    endif

    ilon = int(1 + (180.d0+xlon(i)) )
    if(ilon == 361) ilon = 1
    ! checks
    if(ilon<1 .or. ilon>360) then
      print*,'error in lat/lon range: ilon = ',ilon
      stop 'error in routine icolat/ilon crust1.0'
    endif

    ! gets crust values
    call get_crust_1_0_structure(icolat,ilon,velpl,velsl,rhol,thickl)

    ! sediment thickness
    h_sed = thickl(3) + thickl(4) + thickl(5)

    ! takes upper crust value if sediment too thin
    if( h_sed < MINIMUM_SEDIMENT_THICKNESS ) then
      velpl(3) = velpl(6)
      velpl(4) = velpl(6)
      velpl(5) = velpl(6)

      velsl(3) = velsl(6)
      velsl(4) = velsl(6)
      velsl(5) = velsl(6)

      rhol(3) = rhol(6)
      rhol(4) = rhol(6)
      rhol(5) = rhol(6)
    endif

    ! weighting value
    weightl = weight(i)

    ! total, smoothed values
    rho(:) = rho(:) + weightl*rhol(:)
    thick(:) = thick(:) + weightl*thickl(:)
    velp(:) = velp(:) + weightl*velpl(:)
    vels(:) = vels(:) + weightl*velsl(:)
  enddo

  end subroutine crust_1_0_CAPsmoothed


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_crust_1_0_structure(icolat,ilon,vptyp,vstyp,rhtyp,thtp)

  use model_crust_1_0_par

  implicit none

  ! argument variables
  integer,intent(in) :: icolat,ilon
  double precision,intent(out) :: rhtyp(CRUST1_NP),thtp(CRUST1_NP)
  double precision,intent(out) :: vptyp(CRUST1_NP),vstyp(CRUST1_NP)

  ! sets vp,vs and rho for all layers
  vptyp(:) = crust1_vp(:,icolat,ilon)
  vstyp(:) = crust1_vs(:,icolat,ilon)
  rhtyp(:) = crust1_rho(:,icolat,ilon)
  thtp(:) = crust1_thk(:,icolat,ilon)

  ! get distance to Moho from the bottom of the ocean or the ice
  ! value could be used for checking, but is unused so far...
  thtp(CRUST1_NP) = thtp(CRUST1_NP) - thtp(1) - thtp(2)

  end subroutine get_crust_1_0_structure


