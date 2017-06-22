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

!--------------------------------------------------------------------------------------------------
! General Crustmaps
!
! combines Crust2.0 and EUcrust07 for moho depths; the crustal maps
! interpolate the crustal velocities from Crust2.0 onto the more detailed EUcrust
! crustal depths where ever they are defined.

! current crustmaps (cmaps) take sediment thickness
! and moho depths from EUcrust07 if possible and interpolate corresponding
! velocity/densities given from Crust2.0.
!
! main author: Matthias Meschede (meschede AT princeton DOT edu)
!--------------------------------------------------------------------------------------------------

  module model_crustmaps_par

  ! General Crustmaps parameters
  integer, parameter :: CRUSTMAP_RESOLUTION = 4 !means 1/4 degrees
  integer, parameter :: NLAYERS_CRUSTMAP = 5

  ! model_crustmaps_variables combined crustal maps
  double precision, dimension(:,:,:),allocatable :: thickness,density,velocp,velocs

  double precision,dimension(:),allocatable :: thicknessnp,densitynp, &
    velocpnp,velocsnp,thicknesssp,densitysp,velocpsp,velocssp

  end module model_crustmaps_par

!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_crustmaps_broadcast()

! standard routine to setup model

  use constants
  use model_crustmaps_par

  implicit none

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(thickness(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP), &
           density(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP), &
           velocp(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP), &
           velocs(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION,NLAYERS_CRUSTMAP), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating crustmaps arrays')

  allocate(thicknessnp(NLAYERS_CRUSTMAP), &
           densitynp(NLAYERS_CRUSTMAP), &
           velocpnp(NLAYERS_CRUSTMAP), &
           velocsnp(NLAYERS_CRUSTMAP), &
           thicknesssp(NLAYERS_CRUSTMAP), &
           densitysp(NLAYERS_CRUSTMAP), &
           velocpsp(NLAYERS_CRUSTMAP), &
           velocssp(NLAYERS_CRUSTMAP), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating crustmaps np/sp arrays')

  ! master reads in crust maps
  if (myrank == 0) call read_general_crustmap()

  ! broadcasts values to all processes
  call bcast_all_dp(thickness,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP)
  call bcast_all_dp(velocp,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP)
  call bcast_all_dp(velocs,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP)
  call bcast_all_dp(density,180*360*CRUSTMAP_RESOLUTION*CRUSTMAP_RESOLUTION*NLAYERS_CRUSTMAP)

  ! north pole
  call bcast_all_dp(thicknessnp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(densitynp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(velocpnp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(velocsnp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(densitynp,NLAYERS_CRUSTMAP)

  ! south pole
  call bcast_all_dp(thicknesssp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(densitysp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(velocpsp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(velocssp,NLAYERS_CRUSTMAP)
  call bcast_all_dp(densitysp,NLAYERS_CRUSTMAP)

  end subroutine model_crustmaps_broadcast

!
!-------------------------------------------------------------------------------------------------
!

! read general crustmap by Matthias Meschede

  subroutine read_general_crustmap()

  use constants
  use model_crustmaps_par

  implicit none

  ! local parameters
  integer :: i,k,l
  double precision :: moho,moho_min,moho_max

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: crustMap'
  write(IMAIN,*)
  call flush_IMAIN()

  do l = 1,NLAYERS_CRUSTMAP
    ! file index: from 3 to 7
    i = l + 2

    call read_general_crustmap_layer(thickness(:,:,l), 't', i)
    call read_general_crustmap_layer(density(:,:,l),   'r', i)
    call read_general_crustmap_layer(velocp(:,:,l),    'p', i)
    call read_general_crustmap_layer(velocs(:,:,l),    's', i)
  enddo

  ! North pole / South pole
  thicknessnp(:) = ZERO
  thicknesssp(:) = ZERO
  densitynp(:) = ZERO
  densitysp(:) = ZERO
  velocpnp(:) = ZERO
  velocpsp(:) = ZERO
  velocsnp(:) = ZERO
  velocssp(:) = ZERO

  ! compute average values for North and South pole
  do l = 1,NLAYERS_CRUSTMAP
    do i=1,360*CRUSTMAP_RESOLUTION
      thicknessnp(l) =  thicknessnp(l)+thickness(1,i,l)
      thicknesssp(l) = thicknesssp(l)+thickness(180*CRUSTMAP_RESOLUTION,i,l)
      densitynp(l) = densitynp(l)+density(1,i,l)
      densitysp(l) = densitysp(l)+density(180*CRUSTMAP_RESOLUTION,i,l)
      velocpnp(l) = velocpnp(l)+velocp(1,i,l)
      velocpsp(l) = velocpsp(l)+velocp(180*CRUSTMAP_RESOLUTION,i,l)
      velocsnp(l) = velocsnp(l)+velocs(1,i,l)
      velocssp(l) = velocssp(l)+velocs(180*CRUSTMAP_RESOLUTION,i,l)
    enddo
    thicknessnp(l) = thicknessnp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    thicknesssp(l) = thicknesssp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    densitynp(l) = densitynp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    densitysp(l) = densitysp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    velocpnp(l) = velocpnp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    velocpsp(l) = velocpsp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    velocsnp(l) = velocsnp(l)/360.0/dble(CRUSTMAP_RESOLUTION)
    velocssp(l) = velocssp(l)/360.0/dble(CRUSTMAP_RESOLUTION)

!   print *,'thicknessnp(',l,')',thicknessnp(l)
  enddo

  ! moho thickness = sediment (index 1+2) + upper crust (index 3) + middle crust (index 4) + lower crust (index 5)
  moho_min = HUGEVAL
  moho_max = -HUGEVAL
  do k = 1,180*CRUSTMAP_RESOLUTION
    do i = 1,360*CRUSTMAP_RESOLUTION
      moho = sum(thickness(k,i,:)) !thickness(:,:,1) + thickness(:,:,2)+ thickness(:,:,3) + thickness(:,:,4) + thickness(:,:,5)
      ! limit moho thickness
      if (moho > moho_max) moho_max = moho
      if (moho < moho_min) moho_min = moho
    enddo
  enddo

  ! user output
  write(IMAIN,*) '  Moho crustal thickness (without ice) min/max = ',sngl(moho_min),sngl(moho_max),' km'
  write(IMAIN,*)
  call flush_IMAIN()

  end subroutine read_general_crustmap

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_general_crustmap_layer(var,var_letter,ind)

  use constants
  use model_crustmaps_par, only: CRUSTMAP_RESOLUTION

  implicit none

  double precision, intent(out), dimension(180*CRUSTMAP_RESOLUTION,360*CRUSTMAP_RESOLUTION) :: var
  character(len=1), intent(in) :: var_letter
  integer, intent(in) :: ind

  ! local variables
  character(len=MAX_STRING_LEN) :: eucrust
  integer :: ier, ila, iln

  write(eucrust,'(a,a1,i1,a5)') 'DATA/crustmap/eucrust', var_letter, ind,'.cmap'

  open(unit = IIN,file=trim(eucrust),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(eucrust), '": ', ier
    call flush_IMAIN()
    ! stop
    call exit_MPI(0, 'Error model crustmap')
  endif

  do ila=1,180*CRUSTMAP_RESOLUTION
    read(IIN,*) (var(ila,iln),iln=1,360*CRUSTMAP_RESOLUTION)
  enddo
  close(IIN)

  end subroutine read_general_crustmap_layer


!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_crustmaps(lat,lon,x,vp,vs,rho,moho,found_crust,elem_in_crust)

! Matthias Meschede
! read smooth crust2.0 model (0.25 degree resolution) with eucrust
! based on software routines provided with the crust2.0 model by Bassin et al.

  use constants
  use model_crustmaps_par

  implicit none

  double precision,intent(in) :: lat,lon,x
  double precision,intent(out) :: vp,vs,rho,moho
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust

  double precision :: h_sed,h_uc
  double precision :: x3,x4,x5,x6,x7,scaleval
  double precision,dimension(NLAYERS_CRUSTMAP) :: vps,vss,rhos,thicks

  call read_crustmaps(lat,lon,vps,vss,rhos,thicks)

  x3 = (R_EARTH-thicks(1)*1000.0d0)/R_EARTH
  h_sed = thicks(1) + thicks(2)
  x4 = (R_EARTH-h_sed*1000.0d0)/R_EARTH
  h_uc = h_sed + thicks(3)
  x5 = (R_EARTH-h_uc*1000.0d0)/R_EARTH
  x6 = (R_EARTH-(h_uc+thicks(4))*1000.0d0)/R_EARTH
  x7 = (R_EARTH-(h_uc+thicks(4)+thicks(5))*1000.0d0)/R_EARTH

  found_crust = .true.
! if (x > x3 .and. INCLUDE_SEDIMENTS_IN_CRUST .and. h_sed > MINIMUM_SEDIMENT_THICKNESS) then
  if (x > x3 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
   vp = vps(1)
   vs = vss(1)
   rho = rhos(1)
! else if (x > x4 .and. INCLUDE_SEDIMENTS_IN_CRUST .and. h_sed > MINIMUM_SEDIMENT_THICKNESS) then
  else if (x > x4 .and. INCLUDE_SEDIMENTS_IN_CRUST) then
   vp = vps(2)
   vs = vss(2)
   rho = rhos(2)
  else if (x > x5) then
   vp = vps(3)
   vs = vss(3)
   rho = rhos(3)
  else if (x > x6) then
   vp = vps(4)
   vs = vss(4)
   rho = rhos(4)
  else if (x > x7 .or. elem_in_crust) then
   vp = vps(5)
   vs = vss(5)
   rho = rhos(5)
  else
   found_crust = .false.
  endif

  !   non-dimensionalize
  scaleval = dsqrt(PI*GRAV*RHOAV)

  if (found_crust) then
    vp = vp*1000.0d0/(R_EARTH*scaleval)
    vs = vs*1000.0d0/(R_EARTH*scaleval)
    rho = rho*1000.0d0/RHOAV
  else
    ! takes ficticious values
    vp = 20.0*1000.0d0/(R_EARTH*scaleval)
    vs = 20.0*1000.0d0/(R_EARTH*scaleval)
    rho = 20.0*1000.0d0/RHOAV
  endif

  moho = (h_uc+thicks(4)+thicks(5))*1000.0d0/R_EARTH

  end subroutine model_crustmaps

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_crustmaps(lat,lon,velp,vels,rhos,thicks)

! crustal vp and vs in km/s, layer thickness in km

  use constants
  use model_crustmaps_par

  implicit none

  ! argument variables
  double precision,intent(in) :: lat,lon
  double precision,intent(out) :: rhos(5),thicks(5),velp(5),vels(5)

  !-------------------------------
  ! work-around to avoid Jacobian problems when stretching mesh elements;
  ! one could also try to slightly change the shape of the doubling element bricks (which cause the problem)...
  !
  ! defines a "critical" region to have at least a 1-degree smoothing;
  ! critical region can lead to negative Jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_EUROPE = 50.0d0
  double precision,parameter :: LON_CRITICAL_EUROPE = 22.0d0
  double precision,parameter :: CRITICAL_RANGE_EUROPE = 50.0d0

  ! defines a "critical" region around the andes to have at least a 1-degree smoothing;
  ! critical region can lead to negative Jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_ANDES = -20.0d0
  double precision,parameter :: LON_CRITICAL_ANDES = -70.0d0
  double precision,parameter :: CRITICAL_RANGE_ANDES = 70.0d0

  ! sampling rate for CAP points
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 20
  !-------------------------------

  ! local variables
  double precision,dimension(NTHETA*NPHI) :: xlon,xlat,weight
  double precision,dimension(NLAYERS_CRUSTMAP) :: rhol,thickl,velpl,velsl

  double precision :: lat_used,lon_used
  double precision :: weightup,weightleft,weightul,weightur,weightll,weightlr
  double precision :: weightl,cap_degree,dist
  double precision :: h_sed
  integer :: num_points
  integer :: i,ipoin,iupcolat,ileftlng,irightlng

! get integer colatitude and longitude of crustal cap
! -90 < lat < 90 -180 < lon < 180
  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
    write(*,*) lat,' ',lon, ' error in latitude/longitude range in crust'

  ! avoid rounding at poles
  lat_used = lat
  lon_used = lon
  if (lat == 90.0d0) lat_used = 89.9999d0
  if (lat == -90.0d0) lat_used = -89.9999d0
  if (lon == 180.0d0) lon_used = 179.9999d0
  if (lon == -180.0d0) lon_used = -179.9999d0

  ! by defaults uses only 1 point location
  num_points = 1

  ! initializes cap arrays (to avoid compiler warning)
  xlat(:) = 0.d0
  xlon(:) = 0.d0
  weight(:) = 0.d0

  ! checks if inside/outside of critical region for mesh stretching
  if (SMOOTH_CRUST_EVEN_MORE) then

    dist = dsqrt( (lon_used - LON_CRITICAL_EUROPE)**2 + (lat_used - LAT_CRITICAL_EUROPE )**2 )
    if (dist < CRITICAL_RANGE_EUROPE) then
      ! sets up smoothing points
      ! by default uses CAP smoothing with crustmap resolution, e.g. 1/4 degree
      cap_degree = 1.d0 / CRUSTMAP_RESOLUTION

      ! increases cap smoothing degree
      ! scales between -1 at center and 0 at border
      dist = dist / CRITICAL_RANGE_EUROPE - 1.0d0
      ! shifts value to 1 at center and 0 to the border with exponential decay
      dist = 1.0d0 - exp( - dist*dist*10.0d0 )
      ! increases smoothing degree inside of critical region
      cap_degree = cap_degree + dist

      ! gets smoothing points and weights
      call smooth_weights_CAP_vardegree(lon_used,lat_used,xlon,xlat,weight,cap_degree,NTHETA,NPHI)
      num_points = NTHETA*NPHI
    endif

    dist = dsqrt( (lon_used - LON_CRITICAL_ANDES)**2 + (lat_used - LAT_CRITICAL_ANDES )**2 )
    if (dist < CRITICAL_RANGE_ANDES) then
      ! sets up smoothing points
      ! by default uses CAP smoothing with crustmap resolution, e.g. 1/4 degree
      cap_degree = 1.d0 / CRUSTMAP_RESOLUTION

      ! increases cap smoothing degree
      ! scales between -1 at center and 0 at border
      dist = dist / CRITICAL_RANGE_ANDES - 1.0d0
      ! shifts value to 1 at center and 0 to the border with exponential decay
      dist = 1.0d0 - exp( - dist*dist*10.0d0 )
      ! increases smoothing degree inside of critical region
      cap_degree = cap_degree + dist

      ! gets smoothing points and weights
      call smooth_weights_CAP_vardegree(lon_used,lat_used,xlon,xlat,weight,cap_degree,NTHETA,NPHI)
      num_points = NTHETA*NPHI
    endif

  endif

  ! initializes
  velp(:) = 0.0d0
  vels(:) = 0.0d0
  rhos(:) = 0.0d0
  thicks(:) = 0.0d0

  ! loops over weight points
  do ipoin = 1,num_points
    ! checks if more than one weighting points are taken
    if (num_points > 1) then
      lat_used = xlat(ipoin)
      lon_used = xlon(ipoin)
      ! weighting value
      weightl = weight(ipoin)
    else
      weightl = 1.0d0
    endif

    ! gets crust value indices
    call ibilinearmap(lat_used,lon_used,iupcolat,ileftlng,weightup,weightleft)

    ! interpolates location and crust values
    if (iupcolat == 0) then
       weightup=weightup*2
    else if (iupcolat == 180*CRUSTMAP_RESOLUTION) then
       weightup=2*weightup-1
    endif

    if (ileftlng == 360*CRUSTMAP_RESOLUTION) then
      irightlng = 1
    else
      irightlng = ileftlng+1
    endif

    weightul=weightup*weightleft
    weightur=weightup*(1.0-weightleft)
    weightll=(1.0-weightup)*weightleft
    weightlr=(1.0-weightup)*(1.0-weightleft)

    if (iupcolat == 0) then
      ! North pole
      do i = 1,NLAYERS_CRUSTMAP
       thickl(i)=weightul*thicknessnp(i)+weightur*thicknessnp(i)+&
                 weightll*thickness(1,ileftlng,i)+weightlr*thickness(1,irightlng,i)

       rhol(i)=weightul*densitynp(i)+weightur*densitynp(i)+&
               weightll*density(1,ileftlng,i)+weightlr*density(1,irightlng,i)
       velpl(i)=weightul*velocpnp(i)+weightur*velocpnp(i)+&
               weightll*velocp(1,ileftlng,i)+weightlr*velocp(1,irightlng,i)
       velsl(i)=weightul*velocsnp(i)+weightur*velocsnp(i)+&
               weightll*velocs(1,ileftlng,i)+weightlr*velocs(1,irightlng,i)
      enddo
    else if (iupcolat == 180*CRUSTMAP_RESOLUTION) then
      ! South pole
      do i = 1,NLAYERS_CRUSTMAP
       thickl(i)=weightul*thickness(iupcolat,ileftlng,i)+weightur*thickness(iupcolat,irightlng,i)+&
                 weightll*thicknesssp(i)+weightlr*thicknesssp(i)
       rhol(i)=weightul*density(iupcolat,ileftlng,i)+weightur*density(iupcolat,irightlng,i)+&
               weightll*densitysp(i)+weightlr*densitysp(i)
       velpl(i)=weightul*velocp(iupcolat,ileftlng,i)+weightur*velocp(iupcolat,irightlng,i)+&
               weightll*velocpsp(i)+weightlr*velocpsp(i)
       velsl(i)=weightul*velocs(iupcolat,ileftlng,i)+weightur*velocs(iupcolat,irightlng,i)+&
               weightll*velocssp(i)+weightlr*velocssp(i)
      enddo
    else
      do i = 1,NLAYERS_CRUSTMAP
       thickl(i)=weightul*thickness(iupcolat,ileftlng,i)+weightur*thickness(iupcolat,irightlng,i)+&
                 weightll*thickness(iupcolat+1,ileftlng,i)+weightlr*thickness(iupcolat+1,irightlng,i)
       rhol(i)=weightul*density(iupcolat,ileftlng,i)+weightur*density(iupcolat,irightlng,i)+&
               weightll*density(iupcolat+1,ileftlng,i)+weightlr*density(iupcolat+1,irightlng,i)
       velpl(i)=weightul*velocp(iupcolat,ileftlng,i)+weightur*velocp(iupcolat,irightlng,i)+&
               weightll*velocp(iupcolat+1,ileftlng,i)+weightlr*velocp(iupcolat+1,irightlng,i)
       velsl(i)=weightul*velocs(iupcolat,ileftlng,i)+weightur*velocs(iupcolat,irightlng,i)+&
               weightll*velocs(iupcolat+1,ileftlng,i)+weightlr*velocs(iupcolat+1,irightlng,i)
    !  thicks(i)=1.0
    !  rhos(i)=1.0
    !  velp(i)=1.0
    !  vels(i)=1.0i
      enddo
    endif

    ! sediment thickness
    h_sed = thickl(1) + thickl(2)

    ! takes upper crust value if sediment too thin
    if (h_sed < MINIMUM_SEDIMENT_THICKNESS) then
      velpl(1) = velpl(3)
      velpl(2) = velpl(3)
      velsl(1) = velsl(3)
      velsl(2) = velsl(3)
      rhol(1) = rhol(3)
      rhol(2) = rhol(3)
    endif

    ! total, smoothed values
    rhos(:) = rhos(:) + weightl*rhol(:)
    thicks(:) = thicks(:) + weightl*thickl(:)
    velp(:) = velp(:) + weightl*velpl(:)
    vels(:) = vels(:) + weightl*velsl(:)
  enddo

  end subroutine read_crustmaps

!--------------------------------------------------------------------------------------------

  subroutine ibilinearmap(lat,lng,iupcolat,ileftlng,weightup,weightleft)

  use constants
  use model_crustmaps_par, only: CRUSTMAP_RESOLUTION

  implicit none

  ! argument variables
  double precision weightup,weightleft
  double precision lat,lng, xlng
  double precision buffer
  integer iupcolat
  integer ileftlng

  if (lat > 90.0d0 .or. lat < -90.0d0 .or. lng > 180.0d0 .or. lng < -180.0d0) &
    stop 'Error in latitude/longitude range in ibilinearmap'

! map longitudes to [0,360]
  if (lng < 0) then
    xlng=lng+360.0
  else
    xlng=lng
  endif

  buffer=0.5+((90.0-lat)*CRUSTMAP_RESOLUTION)
  iupcolat=int(buffer)
  weightup=1.0-(buffer-dble(iupcolat))

  if (iupcolat < 0) iupcolat = 0
  if (iupcolat > 180*CRUSTMAP_RESOLUTION)  iupcolat=180*CRUSTMAP_RESOLUTION


  buffer=0.5+(xlng*CRUSTMAP_RESOLUTION)
  ileftlng=int(buffer)
  weightleft=1.0-(buffer-dble(ileftlng))

  if (ileftlng < 1) ileftlng=360*CRUSTMAP_RESOLUTION
  if (ileftlng > 360*CRUSTMAP_RESOLUTION) ileftlng=1

  end subroutine ibilinearmap

