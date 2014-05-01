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
! CRUST 2.0 model by Bassin et al. (2000)
!
! C. Bassin, G. Laske, and G. Masters.
! The current limits of resolution for surface wave tomography in North America.
! EOS, 81: F897, 2000.
!
! The 7 crustal layers:
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
! reads and smooths crust2.0 model
!--------------------------------------------------------------------------------------------------

  module model_crust_par

  ! crustal_model_constants
  ! crustal model parameters for crust2.0
  integer, parameter :: NKEYS_CRUST = 359
  integer, parameter :: NLAYERS_CRUST = 8
  integer, parameter :: NCAP_CRUST = 180

  ! model_crust_variables
  double precision, dimension(:,:),allocatable :: thlr,velocp,velocs,dens
  character(len=2) :: abbreviation(NCAP_CRUST/2,NCAP_CRUST)
  character(len=2) :: code(NKEYS_CRUST)

  end module model_crust_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_crust_broadcast(myrank)

! standard routine to setup model

  use constants
  use model_crust_par

  implicit none

  integer :: myrank
  integer :: ier

  ! allocate crustal arrays
  allocate(thlr(NLAYERS_CRUST,NKEYS_CRUST), &
           velocp(NLAYERS_CRUST,NKEYS_CRUST), &
           velocs(NLAYERS_CRUST,NKEYS_CRUST), &
           dens(NLAYERS_CRUST,NKEYS_CRUST), &
           stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating crustal arrays')

  ! the variables read are declared and stored in structure model_crust_par
  if(myrank == 0) call read_crust_model()

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(thlr,NKEYS_CRUST*NLAYERS_CRUST)
  call bcast_all_dp(velocp,NKEYS_CRUST*NLAYERS_CRUST)
  call bcast_all_dp(velocs,NKEYS_CRUST*NLAYERS_CRUST)
  call bcast_all_dp(dens,NKEYS_CRUST*NLAYERS_CRUST)

  call bcast_all_ch_array2(abbreviation,NCAP_CRUST/2,NCAP_CRUST,2)
  call bcast_all_ch_array(code,NKEYS_CRUST,2)

  end subroutine model_crust_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_crust(lat,lon,x,vp,vs,rho,moho,found_crust,elem_in_crust)

  use constants
  use model_crust_par

  implicit none

  double precision,intent(in) :: lat,lon,x
  double precision,intent(out) :: vp,vs,rho,moho
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust

  ! local parameters
  double precision :: h_sed,h_uc
  double precision :: x3,x4,x5,x6,x7,scaleval
  double precision,dimension(NLAYERS_CRUST):: vps,vss,rhos,thicks

  ! initializes
  vp = ZERO
  vs = ZERO
  rho = ZERO

  ! gets smoothed crust2.0 structure
  call crust_CAPsmoothed(lat,lon,vps,vss,rhos,thicks,abbreviation, &
                        code,thlr,velocp,velocs,dens)

  scaleval = ONE / R_EARTH_KM

  ! non-dimensionalizes thickness (given in km)
  x3 = ONE - thicks(3) * scaleval
  h_sed = thicks(3) + thicks(4)
  x4 = ONE - h_sed * scaleval
  h_uc = h_sed + thicks(5)
  x5 = ONE - h_uc * scaleval
  x6 = ONE - (h_uc+thicks(6)) * scaleval
  x7 = ONE - (h_uc+thicks(6)+thicks(7)) * scaleval

  ! checks moho value
  !moho = h_uc + thicks(6) + thicks(7)
  !if( moho /= thicks(NLAYERS_CRUST) ) then
  ! print*,'moho:',moho,thicks(NLAYERS_CRUST)
  ! print*,'  lat/lon/x:',lat,lon,x
  !endif

  ! No matter found_crust true or false, output moho thickness
  moho = (h_uc+thicks(6)+thicks(7)) * scaleval

  ! gets corresponding crustal velocities and density
  found_crust = .true.
!  if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST &
!   .and. h_sed >= MINIMUM_SEDIMENT_THICKNESS) then
  if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST ) then
    vp = vps(3)
    vs = vss(3)
    rho = rhos(3)
!  else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST &
!   .and. h_sed >= MINIMUM_SEDIMENT_THICKNESS) then
  else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST ) then
    vp = vps(4)
    vs = vss(4)
    rho = rhos(4)
  else if(x > x5) then
    vp = vps(5)
    vs = vss(5)
    rho = rhos(5)
  else if(x > x6) then
    vp = vps(6)
    vs = vss(6)
    rho = rhos(6)
  else if(x > x7 .or. elem_in_crust) then
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
    scaleval = ONE / ( R_EARTH_KM * dsqrt(PI*GRAV*RHOAV) )
    vp = vp * scaleval
    vs = vs * scaleval
    rho = rho * 1000.0d0 / RHOAV
 endif

 end subroutine model_crust

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_crust_model()

  use constants
  use model_crust_par

  implicit none

  ! local variables
  integer :: i,ila,icolat,ikey,ier

  double precision :: h_moho_min,h_moho_max

  character(len=150) :: CNtype2, CNtype2_key_modif

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: CRUST2.0'
  write(IMAIN,*)

  call get_value_string(CNtype2, 'model.CNtype2', 'DATA/crust2.0/CNtype2.txt')
  open(unit=1,file=CNtype2,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "', trim(CNtype2), '": ', ier
    call flush_IMAIN()
    ! stop
    call exit_MPI(0, 'error model crust2.0')
  endif

  do ila=1,NCAP_CRUST/2
    read(1,*) icolat,(abbreviation(ila,i),i=1,NCAP_CRUST)
  enddo
  close(1)

  call get_value_string(CNtype2_key_modif, 'model.CNtype2_key_modif', 'DATA/crust2.0/CNtype2_key_modif.txt')
  open(unit=1,file=CNtype2_key_modif,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "', trim(CNtype2_key_modif), '": ', ier
    call exit_MPI(0, 'error model crust2.0')
  endif

  h_moho_min = HUGEVAL
  h_moho_max = -HUGEVAL

  do ikey=1,NKEYS_CRUST
    read (1,"(a2)") code(ikey)
    read (1,*) (velocp(i,ikey),i=1,NLAYERS_CRUST)
    read (1,*) (velocs(i,ikey),i=1,NLAYERS_CRUST)
    read (1,*) (dens(i,ikey),i=1,NLAYERS_CRUST)
    read (1,*) (thlr(i,ikey),i=1,NLAYERS_CRUST-1),thlr(NLAYERS_CRUST,ikey)

    ! limit moho thickness
    if(thlr(NLAYERS_CRUST,ikey) > h_moho_max) h_moho_max = thlr(NLAYERS_CRUST,ikey)
    if(thlr(NLAYERS_CRUST,ikey) < h_moho_min) h_moho_min = thlr(NLAYERS_CRUST,ikey)
  enddo
  close(1)

  if(h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) &
    stop 'incorrect moho depths in read_crust_model'

  end subroutine read_crust_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine crust_CAPsmoothed(lat,lon,velp,vels,rho,thick,abbreviation,&
                              code,thlr,velocp,velocs,dens)

! crustal vp and vs in km/s, layer thickness in km
!
! crust2.0 is smoothed with a cap of size CAP using NTHETA points
! in the theta direction and NPHI in the phi direction.
! The cap is rotated to the North Pole.

  use constants
  use model_crust_par,only: NLAYERS_CRUST,NKEYS_CRUST,NCAP_CRUST

  implicit none

  ! sampling rate for CAP points
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 20

  ! argument variables
  double precision :: lat,lon
  double precision,dimension(NLAYERS_CRUST) :: rho,thick,velp,vels
  double precision,dimension(NLAYERS_CRUST,NKEYS_CRUST) :: thlr,velocp,velocs,dens

  character(len=2) :: code(NKEYS_CRUST)
  character(len=2) :: abbreviation(NCAP_CRUST/2,NCAP_CRUST)

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
  double precision :: rhol(NLAYERS_CRUST),thickl(NLAYERS_CRUST),velpl(NLAYERS_CRUST),velsl(NLAYERS_CRUST)

  double precision :: weightl,cap_degree,dist
  double precision :: h_sed
  integer :: i,icolat,ilon
  character(len=2) :: crustaltype

  ! small hash table to convert crustal types to key
  integer, dimension(128*128) :: crustalhash_to_key
  integer :: ihash, crustalkey

  ! fill in the hash table
  crustalhash_to_key = -1
  do i=1,NKEYS_CRUST
    call hash_crustal_type(code(i), ihash)
    if (crustalhash_to_key(ihash) /= -1) stop 'error in crust_CAPsmoothed: hash table collision'
    crustalhash_to_key(ihash) = i
  enddo

  ! checks latitude/longitude
  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
    stop 'error in latitude/longitude range in crust'

  ! makes sure lat/lon are within crust2.0 range
  if(lat==90.0d0) lat=89.9999d0
  if(lat==-90.0d0) lat=-89.9999d0
  if(lon==180.0d0) lon=179.9999d0
  if(lon==-180.0d0) lon=-179.9999d0

  ! sets up smoothing points
  ! by default uses CAP smoothing with 1 degree
  cap_degree = 1.0d0

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
  call CAP_vardegree(lon,lat,xlon,xlat,weight,cap_degree,NTHETA,NPHI)

  ! initializes
  velp(:) = ZERO
  vels(:) = ZERO
  rho(:) = ZERO
  thick(:) = ZERO

  ! loops over weight points
  do i=1,NTHETA*NPHI
    ! gets crust values
    call icolat_ilon(xlat(i),xlon(i),icolat,ilon)

    crustaltype = abbreviation(icolat,ilon)

    call hash_crustal_type(crustaltype, ihash)
    crustalkey = crustalhash_to_key(ihash)
    if(crustalkey == -1) stop 'error in retrieving crust type key'

    call get_crust_structure(crustalkey,velpl,velsl,rhol,thickl, &
                            thlr,velocp,velocs,dens)

    ! sediment thickness
    h_sed = thickl(3) + thickl(4)

    ! takes upper crust value if sediment too thin
    if( h_sed < MINIMUM_SEDIMENT_THICKNESS ) then
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

  end subroutine crust_CAPsmoothed

!
!-------------------------------------------------------------------------------------------------
!

! hash table to define the crustal type using an integer instead of characters
! because working with integers in the rest of the routines results in much faster code
  subroutine hash_crustal_type(crustaltype, ihash)

    implicit none

    character(len=2), intent(in) :: crustaltype
    integer, intent(out) :: ihash

    ihash = iachar(crustaltype(1:1)) + 128*iachar(crustaltype(2:2)) + 1

  end subroutine hash_crustal_type

!
!-------------------------------------------------------------------------------------------------
!

  subroutine icolat_ilon(xlat,xlon,icolat,ilon)

  implicit none

! argument variables
  double precision :: xlat,xlon
  integer :: icolat,ilon

  if(xlat > 90.0d0 .or. xlat < -90.0d0 .or. xlon > 180.0d0 .or. xlon < -180.0d0) &
    stop 'error in latitude/longitude range in icolat_ilon'

  icolat = int(1+( (90.d0-xlat)*0.5d0 ))
  if(icolat == 91) icolat = 90

  ilon = int(1+( (180.d0+xlon)*0.5d0 ))
  if(ilon == 181) ilon = 1

  if(icolat>90 .or. icolat<1) stop 'error in routine icolat_ilon'
  if(ilon<1 .or. ilon>180) stop 'error in routine icolat_ilon'

  end subroutine icolat_ilon

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_crust_structure(ikey,vptyp,vstyp,rhtyp,thtp,thlr,velocp,velocs,dens)

  use model_crust_par,only: NLAYERS_CRUST,NKEYS_CRUST

  implicit none

  ! argument variables
  double precision :: rhtyp(NLAYERS_CRUST),thtp(NLAYERS_CRUST)
  double precision :: vptyp(NLAYERS_CRUST),vstyp(NLAYERS_CRUST)
  double precision :: thlr(NLAYERS_CRUST,NKEYS_CRUST),velocp(NLAYERS_CRUST,NKEYS_CRUST)
  double precision :: velocs(NLAYERS_CRUST,NKEYS_CRUST),dens(NLAYERS_CRUST,NKEYS_CRUST)
  integer :: ikey

  ! local variables
  integer :: i

  do i=1,NLAYERS_CRUST
    vptyp(i)=velocp(i,ikey)
    vstyp(i)=velocs(i,ikey)
    rhtyp(i)=dens(i,ikey)
  enddo

  do i=1,NLAYERS_CRUST-1
    thtp(i)=thlr(i,ikey)
  enddo

  !   get distance to Moho from the bottom of the ocean or the ice
  thtp(NLAYERS_CRUST) = thlr(NLAYERS_CRUST,ikey) - thtp(1) - thtp(2)

  end subroutine get_crust_structure


!
!-------------------------------------------------------------------------------------------------
!

  subroutine CAP_vardegree(lon,lat,xlon,xlat,weight,CAP_DEGREE,NTHETA,NPHI)

! calculates weighting points to smooth around lon/lat location with
! a smoothing range of CAP_DEGREE
!
! The cap is rotated to the North Pole.
!
! returns: xlon,xlat,weight

  use constants

  implicit none

  ! sampling rate
  integer :: NTHETA
  integer :: NPHI

  ! smoothing size (in degrees)
  double precision :: CAP_DEGREE

  ! argument variables
  double precision :: lat,lon
  double precision,dimension(NTHETA*NPHI) :: xlon,xlat,weight

  ! local variables
  double precision :: CAP
  double precision :: theta,phi,sint,cost,sinp,cosp,wght,total
  double precision :: r_rot,theta_rot,phi_rot
  double precision :: rotation_matrix(3,3),x(3),xc(3)
  double precision :: dtheta,dphi,cap_area,dweight,pi_over_nphi
  integer :: i,j,k,itheta,iphi

  ! initializes
  xlon(:) = ZERO
  xlat(:) = ZERO
  weight(:) = ZERO

  ! checks cap degree size
  if( CAP_DEGREE < TINYVAL ) then
    ! no cap smoothing
    print*,'error cap:',CAP_DEGREE
    print*,'  lat/lon:',lat,lon
    stop 'error cap_degree too small'
  endif

  ! pre-compute parameters
  CAP = CAP_DEGREE * DEGREES_TO_RADIANS
  dtheta = 0.5d0 * CAP / dble(NTHETA)
  dphi = TWO_PI / dble(NPHI)

  cap_area = TWO_PI * ( ONE - dcos(CAP) )
  dweight = CAP / dble(NTHETA) * dphi / cap_area
  pi_over_nphi = PI/dble(NPHI)

  ! colatitude/longitude in radian
  theta = (90.0d0 - lat ) * DEGREES_TO_RADIANS
  phi = lon * DEGREES_TO_RADIANS

  sint = dsin(theta)
  cost = dcos(theta)
  sinp = dsin(phi)
  cosp = dcos(phi)

  ! set up rotation matrix to go from cap at North pole
  ! to cap around point of interest
  rotation_matrix(1,1) = cosp*cost
  rotation_matrix(1,2) = -sinp
  rotation_matrix(1,3) = cosp*sint
  rotation_matrix(2,1) = sinp*cost
  rotation_matrix(2,2) = cosp
  rotation_matrix(2,3) = sinp*sint
  rotation_matrix(3,1) = -sint
  rotation_matrix(3,2) = ZERO
  rotation_matrix(3,3) = cost

  ! calculates points over a cap at the North pole and rotates them to specified lat/lon point
  i = 0
  total = ZERO
  do itheta = 1,NTHETA

    theta = dble(2*itheta-1)*dtheta
    cost = dcos(theta)
    sint = dsin(theta)
    wght = sint*dweight

    do iphi = 1,NPHI

      i = i+1

      !  get the weight associated with this integration point (same for all phi)
      weight(i) = wght

      total = total + weight(i)
      phi = dble(2*iphi-1)*pi_over_nphi
      cosp = dcos(phi)
      sinp = dsin(phi)

      ! x,y,z coordinates of integration point in cap at North pole
      xc(1) = sint*cosp
      xc(2) = sint*sinp
      xc(3) = cost

      ! get x,y,z coordinates in cap around point of interest
      do j=1,3
        x(j) = ZERO
        do k=1,3
          x(j) = x(j)+rotation_matrix(j,k)*xc(k)
        enddo
      enddo

      ! get latitude and longitude (degrees) of integration point
      call xyz_2_rthetaphi_dble(x(1),x(2),x(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      xlat(i) = (PI_OVER_TWO - theta_rot) * RADIANS_TO_DEGREES
      xlon(i) = phi_rot * RADIANS_TO_DEGREES
      if(xlon(i) > 180.0d0) xlon(i) = xlon(i) - 360.0d0

    enddo

  enddo
  if(abs(total - ONE) > 0.001d0) then
    print*,'error cap:',total,CAP_DEGREE
    stop 'error in cap integration for variable degree'
  endif

  end subroutine CAP_vardegree

