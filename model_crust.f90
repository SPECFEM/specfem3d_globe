!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
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
! reads and smooths crust2.0 model
!--------------------------------------------------------------------------------------------------


  subroutine model_crust_broadcast(myrank,CM_V)

! standard routine to setup model

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  ! model_crust_variables
  type model_crust_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
    character(len=2) dummy_pad ! padding 2 bytes to align the structure
  end type model_crust_variables

  type (model_crust_variables) CM_V
  ! model_crust_variables

  integer :: myrank
  integer :: ier

  ! the variables read are declared and stored in structure CM_V
  if(myrank == 0) call read_crust_model(CM_V)

  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(CM_V%thlr,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(CM_V%velocp,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(CM_V%velocs,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(CM_V%dens,NKEYS_CRUST*NLAYERS_CRUST,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(CM_V%abbreviation,NCAP_CRUST*NCAP_CRUST,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(CM_V%code,2*NKEYS_CRUST,MPI_CHARACTER,0,MPI_COMM_WORLD,ier)


  end subroutine model_crust_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_crust(lat,lon,x,vp,vs,rho,moho,found_crust,CM_V,elem_in_crust)

  implicit none
  include "constants.h"

! model_crust_variables
  type model_crust_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
    character(len=2) dummy_pad ! padding 2 bytes to align the structure
  end type model_crust_variables

  type (model_crust_variables) CM_V
! model_crust_variables

  double precision lat,lon,x,vp,vs,rho,moho
  logical found_crust,elem_in_crust

  ! local parameters
  double precision h_sed,h_uc
  double precision x3,x4,x5,x6,x7,scaleval
  double precision vps(NLAYERS_CRUST),vss(NLAYERS_CRUST),rhos(NLAYERS_CRUST),thicks(NLAYERS_CRUST)

  ! initializes
  vp = 0.d0
  vs = 0.d0
  rho = 0.d0

  ! gets smoothed crust2.0 structure
  call crust_CAPsmoothed(lat,lon,vps,vss,rhos,thicks,CM_V%abbreviation, &
                        CM_V%code,CM_V%thlr,CM_V%velocp,CM_V%velocs,CM_V%dens)

  x3 = (R_EARTH-thicks(3)*1000.0d0)/R_EARTH
  h_sed = thicks(3) + thicks(4)
  x4 = (R_EARTH-h_sed*1000.0d0)/R_EARTH
  h_uc = h_sed + thicks(5)
  x5 = (R_EARTH-h_uc*1000.0d0)/R_EARTH
  x6 = (R_EARTH-(h_uc+thicks(6))*1000.0d0)/R_EARTH
  x7 = (R_EARTH-(h_uc+thicks(6)+thicks(7))*1000.0d0)/R_EARTH

  found_crust = .true.

  if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST &
   .and. h_sed >= MINIMUM_SEDIMENT_THICKNESS) then
    vp = vps(3)
    vs = vss(3)
    rho = rhos(3)
  else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST &
   .and. h_sed >= MINIMUM_SEDIMENT_THICKNESS) then
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
    scaleval = dsqrt(PI*GRAV*RHOAV)
    vp = vp*1000.0d0/(R_EARTH*scaleval)
    vs = vs*1000.0d0/(R_EARTH*scaleval)
    rho = rho*1000.0d0/RHOAV
 endif

 ! checks moho value
 !moho = h_uc + thicks(6) + thicks(7)
 !if( moho /= thicks(NLAYERS_CRUST) ) then
 ! print*,'moho:',moho,thicks(NLAYERS_CRUST)
 ! print*,'  lat/lon/x:',lat,lon,x
 !endif

 ! No matter found_crust true or false, output moho thickness
 moho = (h_uc+thicks(6)+thicks(7))*1000.0d0/R_EARTH

 end subroutine model_crust

!---------------------------

  subroutine read_crust_model(CM_V)

  implicit none
  include "constants.h"

! model_crust_variables
  type model_crust_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
    character(len=2) dummy_pad ! padding 2 bytes to align the structure
  end type model_crust_variables

  type (model_crust_variables) CM_V
! model_crust_variables

! local variables
  integer i
  integer ila,icolat
  integer ikey

  double precision h_moho_min,h_moho_max

  character(len=150) CNtype2, CNtype2_key_modif

  call get_value_string(CNtype2, 'model.CNtype2', 'DATA/crust2.0/CNtype2.txt')
  call get_value_string(CNtype2_key_modif, 'model.CNtype2_key_modif', 'DATA/crust2.0/CNtype2_key_modif.txt')

  open(unit=1,file=CNtype2,status='old',action='read')
  do ila=1,NCAP_CRUST/2
    read(1,*) icolat,(CM_V%abbreviation(ila,i),i=1,NCAP_CRUST)
  enddo
  close(1)

  open(unit=1,file=CNtype2_key_modif,status='old',action='read')
  h_moho_min=HUGEVAL
  h_moho_max=-HUGEVAL
  do ikey=1,NKEYS_CRUST
    read (1,"(a2)") CM_V%code(ikey)
    read (1,*) (CM_V%velocp(ikey,i),i=1,NLAYERS_CRUST)
    read (1,*) (CM_V%velocs(ikey,i),i=1,NLAYERS_CRUST)
    read (1,*) (CM_V%dens(ikey,i),i=1,NLAYERS_CRUST)
    read (1,*) (CM_V%thlr(ikey,i),i=1,NLAYERS_CRUST-1),CM_V%thlr(ikey,NLAYERS_CRUST)
    if(CM_V%thlr(ikey,NLAYERS_CRUST) > h_moho_max) h_moho_max=CM_V%thlr(ikey,NLAYERS_CRUST)
    if(CM_V%thlr(ikey,NLAYERS_CRUST) < h_moho_min) h_moho_min=CM_V%thlr(ikey,NLAYERS_CRUST)
  enddo
  close(1)

  if(h_moho_min == HUGEVAL .or. h_moho_max == -HUGEVAL) &
    stop 'incorrect moho depths in read_crust_model'

  end subroutine read_crust_model

!---------------------------

  subroutine crust_CAPsmoothed(lat,lon,velp,vels,rho,thick,abbreviation,&
                              code,thlr,velocp,velocs,dens)

! crustal vp and vs in km/s, layer thickness in km
!
! crust2.0 is smoothed with a cap of size CAP using NTHETA points
! in the theta direction and NPHI in the phi direction.
! The cap is rotated to the North Pole.

  implicit none
  include "constants.h"

  ! sampling rate for CAP points
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 20

  ! argument variables
  double precision lat,lon
  double precision rho(NLAYERS_CRUST),thick(NLAYERS_CRUST),velp(NLAYERS_CRUST),vels(NLAYERS_CRUST)
  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)
  character(len=2) code(NKEYS_CRUST),abbreviation(NCAP_CRUST/2,NCAP_CRUST)

  !-------------------------------
  ! work-around to avoid jacobian problems when stretching mesh elements;
  ! one could also try to slightly change the shape of the doulbing element bricks (which cause the problem)...
  !
  ! defines a "critical" region around the andes to have at least a 2-degree smoothing;
  ! critical region can lead to negative jacobians for mesh stretching when CAP smoothing is too small
  double precision,parameter :: LAT_CRITICAL_ANDES = -20.0d0
  double precision,parameter :: LON_CRITICAL_ANDES = -70.0d0
  double precision,parameter :: CRITICAL_RANGE = 70.0d0
  !-------------------------------

  ! local variables
  double precision xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)
  double precision rhol(NLAYERS_CRUST),thickl(NLAYERS_CRUST),velpl(NLAYERS_CRUST),velsl(NLAYERS_CRUST)
  double precision weightl,cap_degree,dist
  integer i,icolat,ilon,ierr
  character(len=2) crustaltype

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
      dist = dist / CRITICAL_RANGE - 1.0d0
      ! shifts value to 1 at center and 0 to the border with exponential decay
      dist = 1.0d0 - exp( - dist*dist*10.0d0 )
      ! increases smoothing degree inside of critical region to 2 degree
      cap_degree = cap_degree + dist
    endif
  endif

  ! gets smoothing points and weights
  call CAP_vardegree(lon,lat,xlon,xlat,weight,cap_degree,NTHETA,NPHI)

  ! initializes
  velp(:) = 0.0d0
  vels(:) = 0.0d0
  rho(:) = 0.0d0
  thick(:) = 0.0d0

  ! loops over weight points
  do i=1,NTHETA*NPHI
    ! gets crust values
    call icolat_ilon(xlat(i),xlon(i),icolat,ilon)
    crustaltype = abbreviation(icolat,ilon)
    call get_crust_structure(crustaltype,velpl,velsl,rhol,thickl, &
                            code,thlr,velocp,velocs,dens,ierr)
    if(ierr /= 0) stop 'error in routine get_crust_structure'

    ! weighting value
    weightl = weight(i)

    ! total, smoothed values
    rho(:) = rho(:) + weightl*rhol(:)
    thick(:) = thick(:) + weightl*thickl(:)
    velp(:) = velp(:) + weightl*velpl(:)
    vels(:) = vels(:) + weightl*velsl(:)
  enddo

  end subroutine crust_CAPsmoothed


!------------------------------------------------------

  subroutine icolat_ilon(xlat,xlon,icolat,ilon)

  implicit none


! argument variables
  double precision xlat,xlon
  integer icolat,ilon

  if(xlat > 90.0d0 .or. xlat < -90.0d0 .or. xlon > 180.0d0 .or. xlon < -180.0d0) &
    stop 'error in latitude/longitude range in icolat_ilon'
  icolat=int(1+((90.d0-xlat)/2.d0))
  if(icolat == 91) icolat=90
  ilon=int(1+((180.d0+xlon)/2.d0))
  if(ilon == 181) ilon=1

  if(icolat>90 .or. icolat<1) stop 'error in routine icolat_ilon'
  if(ilon<1 .or. ilon>180) stop 'error in routine icolat_ilon'

  end subroutine icolat_ilon

!---------------------------------------------------------------------

  subroutine get_crust_structure(type,vptyp,vstyp,rhtyp,thtp, &
               code,thlr,velocp,velocs,dens,ierr)

  implicit none
  include "constants.h"


! argument variables
  integer ierr
  double precision rhtyp(NLAYERS_CRUST),thtp(NLAYERS_CRUST)
  double precision vptyp(NLAYERS_CRUST),vstyp(NLAYERS_CRUST)
  character(len=2) type,code(NKEYS_CRUST)
  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)

! local variables
  integer i,ikey

  ierr=1
  do ikey=1,NKEYS_CRUST
    if (code(ikey) == type) then
      do i=1,NLAYERS_CRUST
        vptyp(i)=velocp(ikey,i)
        vstyp(i)=velocs(ikey,i)
        rhtyp(i)=dens(ikey,i)
      enddo
      do i=1,NLAYERS_CRUST-1
        thtp(i)=thlr(ikey,i)
      enddo
      !   get distance to Moho from the bottom of the ocean or the ice
      thtp(NLAYERS_CRUST)=thlr(ikey,NLAYERS_CRUST)-thtp(1)-thtp(2)
      ierr=0
    endif
  enddo

  end subroutine get_crust_structure


!---------------------------

  subroutine CAP_vardegree(lon,lat,xlon,xlat,weight,CAP_DEGREE,NTHETA,NPHI)

! calculates weighting points to smooth around lon/lat location with
! a smoothing range of CAP_DEGREE
!
! The cap is rotated to the North Pole.
!
! returns: xlon,xlat,weight

  implicit none
  include "constants.h"

  ! sampling rate
  integer :: NTHETA
  integer :: NPHI
  ! smoothing size (in degrees)
  double precision :: CAP_DEGREE

  ! argument variables
  double precision lat,lon
  double precision xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)

  ! local variables
  double precision CAP
  double precision theta,phi,sint,cost,sinp,cosp,wght,total
  double precision r_rot,theta_rot,phi_rot
  double precision rotation_matrix(3,3),x(3),xc(3)
  double precision dtheta,dphi,cap_area,dweight,pi_over_nphi
  integer i,j,k
  integer itheta,iphi

  double precision, parameter :: RADIANS_TO_DEGREES = 180.d0 / PI
  double precision, parameter :: PI_OVER_TWO = PI / 2.0d0

  ! initializes
  xlon(:) = 0.d0
  xlat(:) = 0.d0
  weight(:) = 0.d0

  ! checks cap degree size
  if( CAP_DEGREE < TINYVAL ) then
    ! no cap smoothing
    print*,'error cap:',CAP_DEGREE
    print*,'  lat/lon:',lat,lon
    stop 'error cap_degree too small'
  endif

  ! pre-compute parameters
  CAP = CAP_DEGREE * PI/180.0d0
  dtheta = 0.5d0 * CAP / dble(NTHETA)
  dphi = TWO_PI / dble(NPHI)
  cap_area = TWO_PI * (1.0d0 - dcos(CAP))
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
  rotation_matrix(3,2) = 0.0d0
  rotation_matrix(3,3) = cost

  ! calculates points over a cap at the North pole and rotates them to specified lat/lon point
  i = 0
  total = 0.0d0
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
        x(j) = 0.0d0
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
  if(abs(total-1.0d0) > 0.001d0) then
    print*,'error cap:',total,CAP_DEGREE
    stop 'error in cap integration for variable degree'
  endif

  end subroutine


!---------------------------
! unused routines...
!
!  subroutine crust_singlevalue(lat,lon,velp,vels,rho,thick,abbreviation,&
!                              code,thlr,velocp,velocs,dens)
!
!! crustal vp and vs in km/s, layer thickness in km
!
!!  uses crust2.0 as is, without smoothing
!
!  implicit none
!  include "constants.h"
!
!! argument variables
!  double precision lat,lon
!  double precision rho(NLAYERS_CRUST),thick(NLAYERS_CRUST),velp(NLAYERS_CRUST),vels(NLAYERS_CRUST)
!  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
!  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)
!  character(len=2) code(NKEYS_CRUST),abbreviation(NCAP_CRUST/2,NCAP_CRUST)
!
!! local variables
!  integer icolat,ilon,ierr
!  character(len=2) crustaltype
!
!
!! get integer colatitude and longitude of crustal cap
!! -90<lat<90 -180<lon<180
!  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
!    stop 'error in latitude/longitude range in crust'
!  if(lat==90.0d0) lat=89.9999d0
!  if(lat==-90.0d0) lat=-89.9999d0
!  if(lon==180.0d0) lon=179.9999d0
!  if(lon==-180.0d0) lon=-179.9999d0
!
!  call icolat_ilon(lat,lon,icolat,ilon)
!  crustaltype = abbreviation(icolat,ilon)
!  call get_crust_structure(crustaltype,velp,vels,rho,thick, &
!                          code,thlr,velocp,velocs,dens,ierr)
!  if( ierr /= 0 ) stop 'error in routine get_crust_structure'
!
!  end subroutine crust_singlevalue
!
!---------------------------
!
!
!  subroutine crust_org(lat,lon,velp,vels,rho,thick,abbreviation,code,thlr,velocp,velocs,dens)
!
!! crustal vp and vs in km/s, layer thickness in km
!! crust2.0 is smoothed with a cap of size CAP using NTHETA points
!! in the theta direction and NPHI in the phi direction.
!! The cap is rotated to the North Pole.
!
!  implicit none
!  include "constants.h"
!! Change the CAP function to smooth crustal model
!  integer, parameter :: NTHETA = 4         !2
!  integer, parameter :: NPHI = 20          !10
!  double precision, parameter :: CAP = 1.0d0*PI/180.0d0   ! 2.0d0*PI/180.0d0
!
!! argument variables
!  double precision lat,lon
!  double precision rho(NLAYERS_CRUST),thick(NLAYERS_CRUST),velp(NLAYERS_CRUST),vels(NLAYERS_CRUST)
!  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
!  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)
!  character(len=2) code(NKEYS_CRUST),abbreviation(NCAP_CRUST/2,NCAP_CRUST)
!
!! local variables
!  integer i,j,k,icolat,ilon,ierr
!  integer itheta,iphi,npoints
!  double precision theta,phi,sint,cost,sinp,cosp,dtheta,dphi,cap_area,wght,total
!  double precision r_rot,theta_rot,phi_rot
!  double precision rotation_matrix(3,3),x(3),xc(3)
!  double precision xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)
!  double precision rhol(NLAYERS_CRUST),thickl(NLAYERS_CRUST),velpl(NLAYERS_CRUST),velsl(NLAYERS_CRUST)
!  character(len=2) crustaltype
!
!! get integer colatitude and longitude of crustal cap
!! -90<lat<90 -180<lon<180
!  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
!    stop 'error in latitude/longitude range in crust'
!  if(lat==90.0d0) lat=89.9999d0
!  if(lat==-90.0d0) lat=-89.9999d0
!  if(lon==180.0d0) lon=179.9999d0
!  if(lon==-180.0d0) lon=-179.9999d0
!
!  call icolat_ilon(lat,lon,icolat,ilon)
!  crustaltype=abbreviation(icolat,ilon)
!  call get_crust_structure(crustaltype,velp,vels,rho,thick, &
!                    code,thlr,velocp,velocs,dens,ierr)
!
!!  uncomment the following line to use crust2.0 as is, without smoothing
!!
!!  return
!
!  theta = (90.0-lat)*PI/180.0
!  phi = lon*PI/180.0
!
!  sint = sin(theta)
!  cost = cos(theta)
!  sinp = sin(phi)
!  cosp = cos(phi)
!
!! set up rotation matrix to go from cap at North pole
!! to cap around point of interest
!  rotation_matrix(1,1) = cosp*cost
!  rotation_matrix(1,2) = -sinp
!  rotation_matrix(1,3) = cosp*sint
!  rotation_matrix(2,1) = sinp*cost
!  rotation_matrix(2,2) = cosp
!  rotation_matrix(2,3) = sinp*sint
!  rotation_matrix(3,1) = -sint
!  rotation_matrix(3,2) = 0.0
!  rotation_matrix(3,3) = cost
!
!  dtheta = CAP/dble(NTHETA)
!  dphi = 2.0*PI/dble(NPHI)
!  cap_area = 2.0*PI*(1.0-cos(CAP))
!
!! integrate over a cap at the North pole
!  i = 0
!  total = 0.0
!  do itheta = 1,NTHETA
!
!    theta = 0.5*dble(2*itheta-1)*CAP/dble(NTHETA)
!    cost = cos(theta)
!    sint = sin(theta)
!    wght = sint*dtheta*dphi/cap_area
!
!    do iphi = 1,NPHI
!
!      i = i+1
!!     get the weight associated with this integration point (same for all phi)
!      weight(i) = wght
!      total = total + weight(i)
!      phi = dble(2*iphi-1)*PI/dble(NPHI)
!      cosp = cos(phi)
!      sinp = sin(phi)
!!     x,y,z coordinates of integration point in cap at North pole
!      xc(1) = sint*cosp
!      xc(2) = sint*sinp
!      xc(3) = cost
!!     get x,y,z coordinates in cap around point of interest
!      do j=1,3
!        x(j) = 0.0
!        do k=1,3
!          x(j) = x(j)+rotation_matrix(j,k)*xc(k)
!        enddo
!      enddo
!!     get latitude and longitude (degrees) of integration point
!      call xyz_2_rthetaphi_dble(x(1),x(2),x(3),r_rot,theta_rot,phi_rot)
!      call reduce(theta_rot,phi_rot)
!      xlat(i) = (PI/2.0-theta_rot)*180.0/PI
!      xlon(i) = phi_rot*180.0/PI
!      if(xlon(i) > 180.0) xlon(i) = xlon(i)-360.0
!
!    enddo
!
!  enddo
!
!  if(abs(total-1.0) > 0.001) stop 'error in cap integration for crust2.0'
!
!  npoints = i
!
!  do j=1,NLAYERS_CRUST
!    rho(j)=0.0d0
!    thick(j)=0.0d0
!    velp(j)=0.0d0
!    vels(j)=0.0d0
!  enddo
!
!  do i=1,npoints
!    call icolat_ilon(xlat(i),xlon(i),icolat,ilon)
!    crustaltype=abbreviation(icolat,ilon)
!    call get_crust_structure(crustaltype,velpl,velsl,rhol,thickl, &
!                    code,thlr,velocp,velocs,dens,ierr)
!    if(ierr /= 0) stop 'error in routine get_crust_structure'
!    do j=1,NLAYERS_CRUST
!      rho(j)=rho(j)+weight(i)*rhol(j)
!      thick(j)=thick(j)+weight(i)*thickl(j)
!      velp(j)=velp(j)+weight(i)*velpl(j)
!      vels(j)=vels(j)+weight(i)*velsl(j)
!    enddo
!  enddo
!
!  end subroutine crust_org

