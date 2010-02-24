!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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
! EUCRUST-07
!
! Tesauro, M., M. K. Kaban and S. A. P. L. Cloetingh, 2008. 
! Eucrust-07: A New Reference Model for the European Crust, 
! Geophysical Research Letters, 35: p. L05313.208 
!--------------------------------------------------------------------------------------------------

  subroutine model_eucrust_broadcast(myrank,EUCM_V)

! standard routine to setup model 

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  ! EUcrust
  type model_eucrust_variables
    double precision, dimension(:),pointer :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
    integer :: num_eucrust
  end type model_eucrust_variables
  type (model_eucrust_variables) EUCM_V

  integer :: myrank
  integer :: ier
  
  ! EUcrust07 Vp crustal structure
  if( myrank == 0 ) call read_EuCrust(EUCM_V)
  
  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(EUCM_V%num_eucrust,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)    
  
  if( myrank /= 0 ) then
    allocate(EUCM_V%eucrust_vp_uppercrust(EUCM_V%num_eucrust),EUCM_V%eucrust_vp_lowercrust(EUCM_V%num_eucrust),&
            EUCM_V%eucrust_mohodepth(EUCM_V%num_eucrust),EUCM_V%eucrust_basement(EUCM_V%num_eucrust),&
            EUCM_V%eucrust_ucdepth(EUCM_V%num_eucrust), EUCM_V%eucrust_lon(EUCM_V%num_eucrust),&
            EUCM_V%eucrust_lat(EUCM_V%num_eucrust))
  endif      
  
  call MPI_BCAST(EUCM_V%eucrust_lat(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(EUCM_V%eucrust_lon(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(EUCM_V%eucrust_vp_uppercrust(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(EUCM_V%eucrust_vp_lowercrust(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(EUCM_V%eucrust_mohodepth(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(EUCM_V%eucrust_basement(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(EUCM_V%eucrust_ucdepth(1:EUCM_V%num_eucrust),EUCM_V%num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine model_eucrust_broadcast
  
!----------------------------------------------------------------------------------------------------

  subroutine read_EuCrust(EUCM_V)

  implicit none

  include "constants.h"

  type model_eucrust_variables
    double precision, dimension(:),pointer :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
    integer :: num_eucrust
  end type model_eucrust_variables
  type (model_eucrust_variables) EUCM_V

    
  ! local variables
  character(len=80):: line
  character(len=150):: filename
  integer:: i,ierror
  double precision:: vp_uppercrust,vp_lowercrust,vp_avg,topo,basement
  double precision:: upper_lower_depth,moho_depth,lat,lon
    
  ! original file size entries    
  EUCM_V%num_eucrust = 36058
  
  allocate(EUCM_V%eucrust_vp_uppercrust(EUCM_V%num_eucrust),EUCM_V%eucrust_vp_lowercrust(EUCM_V%num_eucrust),&
        EUCM_V%eucrust_mohodepth(EUCM_V%num_eucrust),EUCM_V%eucrust_basement(EUCM_V%num_eucrust),&
        EUCM_V%eucrust_ucdepth(EUCM_V%num_eucrust), EUCM_V%eucrust_lon(EUCM_V%num_eucrust),&
        EUCM_V%eucrust_lat(EUCM_V%num_eucrust))
        
  EUCM_V%eucrust_vp_uppercrust(:) = 0.0
  EUCM_V%eucrust_vp_lowercrust(:) = 0.0
  EUCM_V%eucrust_mohodepth(:) = 0.0
  EUCM_V%eucrust_basement(:) = 0.0
  EUCM_V%eucrust_ucdepth(:) = 0.0
  
  ! opens data file
  call get_value_string(filename, 'model.eu', 'DATA/eucrust-07/ds01.txt')
  open(unit=11,file=filename,status='old',action='read')
  
  ! skip first line      
  read(11,*) 
    
  ! data
  do i=1,36058

    read(11,'(a80)',iostat=ierror) line
    if(ierror .ne. 0 ) stop

    read(line,*)lon,lat,vp_uppercrust,vp_lowercrust,vp_avg,topo,basement,upper_lower_depth,moho_depth
 
    ! stores moho values
    EUCM_V%eucrust_lon(i) = lon
    EUCM_V%eucrust_lat(i) = lat
    EUCM_V%eucrust_vp_uppercrust(i) = vp_uppercrust
    EUCM_V%eucrust_vp_lowercrust(i) = vp_lowercrust
    EUCM_V%eucrust_mohodepth(i) = moho_depth     
    EUCM_V%eucrust_basement(i) = basement     
    EUCM_V%eucrust_ucdepth(i) = upper_lower_depth
    
  enddo
  close(11)

  end subroutine read_EuCrust

!  
!--------------------------------------------------------------------------------------------------
!

  subroutine model_eucrust(lat,lon,x,vp,found_crust,EUCM_V)

  implicit none

  type model_eucrust_variables
    double precision, dimension(:),pointer :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
    integer :: num_eucrust
  end type model_eucrust_variables
  type (model_eucrust_variables) EUCM_V

  double precision :: lat,lon,x,vp 
  logical :: found_crust  
  double precision :: lon_min,lon_max,lat_min,lat_max  
  double precision, external:: crust_eu
  
  ! initializes
  vp = 0.d0
  
  ! eucrust boundary region
  lon_min = -24.875
  lon_max = 35.375
      
  lat_min = 34.375
  lat_max = 71.375

  found_crust = .false.
  if( lon < lon_min .or. lon > lon_max ) return
  if( lat < lat_min .or. lat > lat_max ) return

  ! smoothing over 1.0 degrees
  call eu_cap_smoothing(lat,lon,x,vp,found_crust,EUCM_V)

  ! without smoothing
  !vp = crust_eu(lat,lon,x,vp,found_crust,EUCM_V)

  end subroutine model_eucrust

!  
!--------------------------------------------------------------------------------------------------
!

  double precision function crust_eu(lat,lon,x,vp,found_crust,EUCM_V)

! returns Vp at the specific location lat/lon

  implicit none

  include "constants.h"

  type model_eucrust_variables
    double precision, dimension(:),pointer :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
    integer :: num_eucrust
  end type model_eucrust_variables
  type (model_eucrust_variables) EUCM_V

  double precision :: lat,lon,x,vp !,vs,rho,moho
  logical :: found_crust
  
  double precision :: longitude_min,longitude_max,latitude_min,latitude_max
  double precision :: h_basement,h_uc,h_moho,x3,x4,x5
  double precision :: scaleval
  
  integer :: i,j
  integer,parameter :: ilons = 242  ! number of different longitudes
  integer,parameter :: ilats = 149  ! number of different latitudes
  
  ! eucrust boundary region
  longitude_min = -24.875
  longitude_max = 35.375
      
  latitude_min = 34.375
  latitude_max = 71.375

  found_crust = .false.
  crust_eu = 0.0
  if( lon < longitude_min .or. lon > longitude_max ) return
  if( lat < latitude_min .or. lat > latitude_max ) return

  ! search
  do i=1,ilons-1   
    if( lon >= EUCM_V%eucrust_lon(i) .and. lon < EUCM_V%eucrust_lon(i+1) ) then
          do j=0,ilats-1
            if(lat>=EUCM_V%eucrust_lat(i+j*ilons) .and. lat<EUCM_V%eucrust_lat(i+(j+1)*ilons)) then
              
              h_basement = EUCM_V%eucrust_basement(i+j*ilons)
              h_uc = EUCM_V%eucrust_ucdepth(i+j*ilons)
              h_moho = EUCM_V%eucrust_mohodepth(i+j*ilons)
              
              x3=(R_EARTH - h_basement*1000.0d0)/R_EARTH
              x4=(R_EARTH - h_uc*1000.0d0)/R_EARTH
              x5=(R_EARTH - h_moho*1000.0d0)/R_EARTH
              
              scaleval = dsqrt(PI*GRAV*RHOAV)
              
              if( x > x3 ) then
                return
              else if( x > x4 ) then              
                found_crust = .true.
                vp = EUCM_V%eucrust_vp_uppercrust(i+j*ilons) *1000.0d0/(R_EARTH*scaleval)
                crust_eu = vp
                return
              else if( x > x5 ) then
                found_crust = .true.
                vp = EUCM_V%eucrust_vp_lowercrust(i+j*ilons) *1000.0d0/(R_EARTH*scaleval)
                crust_eu = vp
                return
              endif
              return
            endif
          enddo
        endif
      enddo

  end function crust_eu

!  
!--------------------------------------------------------------------------------------------------
!
  subroutine eu_cap_smoothing(lat,lon,radius,value,found,EUCM_V)

! smooths with a cap of size CAP (in degrees)
! using NTHETA points in the theta direction (latitudal) 
! and NPHI in the phi direction (longitudal).
! The cap is rotated to the North Pole.

  implicit none
  include "constants.h"

  ! argument variables
  double precision lat,lon,radius
  double precision :: value
  logical :: found

  type model_eucrust_variables
    double precision, dimension(:),pointer :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
    integer :: num_eucrust
  end type model_eucrust_variables
  type (model_eucrust_variables) EUCM_V
  
  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 10
  double precision, parameter :: CAP = 1.0d0*PI/180.0d0   ! 1 degree smoothing

  double precision,external :: crust_eu
  
  ! local variables
  integer i,j,k !,icolat,ilon,ierr
  integer itheta,iphi,npoints
  double precision theta,phi,sint,cost,sinp,cosp,dtheta,dphi,cap_area,wght,total,valuel
  double precision r_rot,theta_rot,phi_rot
  double precision rotation_matrix(3,3),x(3),xc(3)
  double precision xlon(NTHETA*NPHI),xlat(NTHETA*NPHI),weight(NTHETA*NPHI)

  ! get integer colatitude and longitude of crustal cap
  ! -90<lat<90 -180<lon<180
  if(lat > 90.0d0 .or. lat < -90.0d0 .or. lon > 180.0d0 .or. lon < -180.0d0) &
    stop 'error in latitude/longitude range in crust'    
  if(lat==90.0d0) lat=89.9999d0
  if(lat==-90.0d0) lat=-89.9999d0
  if(lon==180.0d0) lon=179.9999d0
  if(lon==-180.0d0) lon=-179.9999d0

  !call icolat_ilon(lat,lon,icolat,ilon)
  !crustaltype=abbreviation(icolat,ilon)
  !call get_crust_structure(crustaltype,velp,vels,rho,thick,code,thlr,velocp,velocs,dens,ierr)

  !  uncomment the following line to use as is, without smoothing
  !  value = func(lat,lon,x,value,found,EUCM_V)
  !  return

  theta = (90.0-lat)*PI/180.0
  phi = lon*PI/180.0

  sint = sin(theta)
  cost = cos(theta)
  sinp = sin(phi)
  cosp = cos(phi)

  ! set up rotation matrix to go from cap at North pole
  ! to cap around point of interest
  rotation_matrix(1,1) = cosp*cost
  rotation_matrix(1,2) = -sinp
  rotation_matrix(1,3) = cosp*sint
  rotation_matrix(2,1) = sinp*cost
  rotation_matrix(2,2) = cosp
  rotation_matrix(2,3) = sinp*sint
  rotation_matrix(3,1) = -sint
  rotation_matrix(3,2) = 0.0
  rotation_matrix(3,3) = cost

  dtheta = CAP/dble(NTHETA)
  dphi = 2.0*PI/dble(NPHI)
  cap_area = 2.0*PI*(1.0-cos(CAP))

  ! integrate over a cap at the North pole
  i = 0
  total = 0.0
  do itheta = 1,NTHETA

    theta = 0.5*dble(2*itheta-1)*CAP/dble(NTHETA)
    cost = cos(theta)
    sint = sin(theta)
    wght = sint*dtheta*dphi/cap_area

    do iphi = 1,NPHI

      i = i+1
      !get the weight associated with this integration point (same for all phi)
      weight(i) = wght
      total = total + weight(i)
      phi = dble(2*iphi-1)*PI/dble(NPHI)
      cosp = cos(phi)
      sinp = sin(phi)
      !     x,y,z coordinates of integration point in cap at North pole
      xc(1) = sint*cosp
      xc(2) = sint*sinp
      xc(3) = cost
      !     get x,y,z coordinates in cap around point of interest
      do j=1,3
        x(j) = 0.0
        do k=1,3
          x(j) = x(j)+rotation_matrix(j,k)*xc(k)
        enddo
      enddo
      !     get latitude and longitude (degrees) of integration point
      call xyz_2_rthetaphi_dble(x(1),x(2),x(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      xlat(i) = (PI/2.0-theta_rot)*180.0/PI
      xlon(i) = phi_rot*180.0/PI
      if(xlon(i) > 180.0) xlon(i) = xlon(i)-360.0

    enddo

  enddo

  if(abs(total-1.0) > 0.001) stop 'error in cap integration for crust2.0'

  npoints = i

  ! at this point: 
  !
  ! xlat(i),xlon(i) are point locations to be used for interpolation
  ! with weights weight(i)

  ! integrates value
  value = 0.0d0  
  do i=1,npoints        
    valuel = crust_eu(xlat(i),xlon(i),radius,value,found,EUCM_V)        
    value = value + weight(i)*valuel    
  enddo

  if( abs(value) < TINYVAL) found = .false.

  end subroutine eu_cap_smoothing
  
  