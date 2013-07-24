!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

  module model_eucrust_par

  ! EUcrust
  double precision, dimension(:),allocatable :: eucrust_lat,eucrust_lon,&
      eucrust_vp_uppercrust,eucrust_vp_lowercrust,eucrust_mohodepth,&
      eucrust_basement,eucrust_ucdepth
  integer :: num_eucrust

  end module model_eucrust_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_eucrust_broadcast(myrank)

! standard routine to setup model

  use mpi
  use model_eucrust_par

  implicit none

  include "constants.h"

  integer :: myrank
  integer :: ier

  ! EUcrust07 Vp crustal structure
  if( myrank == 0 ) call read_EuCrust()

  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(num_eucrust,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! allocates on all other processes
  if( myrank /= 0 ) then
    allocate(eucrust_vp_uppercrust(num_eucrust), &
            eucrust_vp_lowercrust(num_eucrust),&
            eucrust_mohodepth(num_eucrust), &
            eucrust_basement(num_eucrust),&
            eucrust_ucdepth(num_eucrust), &
            eucrust_lon(num_eucrust),&
            eucrust_lat(num_eucrust), &
            stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating EUcrust arrays')
  endif

  call MPI_BCAST(eucrust_lat(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eucrust_lon(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eucrust_vp_uppercrust(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eucrust_vp_lowercrust(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eucrust_mohodepth(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eucrust_basement(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(eucrust_ucdepth(1:num_eucrust),num_eucrust,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)

  end subroutine model_eucrust_broadcast

!----------------------------------------------------------------------------------------------------

  subroutine read_EuCrust()

  use model_eucrust_par

  implicit none

  include "constants.h"

  ! local variables
  character(len=80):: line
  character(len=150):: filename
  integer:: i,ier
  double precision:: vp_uppercrust,vp_lowercrust,vp_avg,topo,basement
  double precision:: upper_lower_depth,moho_depth,lat,lon

  ! original file size entries
  num_eucrust = 36058

  ! only on master we allocate these arrays here
  allocate(eucrust_vp_uppercrust(num_eucrust), &
          eucrust_vp_lowercrust(num_eucrust),&
          eucrust_mohodepth(num_eucrust), &
          eucrust_basement(num_eucrust),&
          eucrust_ucdepth(num_eucrust), &
          eucrust_lon(num_eucrust),&
          eucrust_lat(num_eucrust), &
          stat=ier)
  if( ier /= 0 ) call exit_MPI(0,'error allocating EUcrust arrays on master')

  eucrust_vp_uppercrust(:) = ZERO
  eucrust_vp_lowercrust(:) = ZERO
  eucrust_mohodepth(:) = ZERO
  eucrust_basement(:) = ZERO
  eucrust_ucdepth(:) = ZERO

  ! opens data file
  call get_value_string(filename, 'model.eu', 'DATA/eucrust-07/ds01.txt')
  open(unit=11,file=filename,status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(0,'error opening EUcrust file')

  ! skip first line
  read(11,*)

  ! data
  do i=1,36058

    read(11,'(a80)',iostat=ier) line
    if( ier /= 0 ) stop 'error reading EUcrust file'

    read(line,*)lon,lat,vp_uppercrust,vp_lowercrust,vp_avg,topo,basement,upper_lower_depth,moho_depth

    ! stores moho values
    eucrust_lon(i) = lon
    eucrust_lat(i) = lat
    eucrust_vp_uppercrust(i) = vp_uppercrust
    eucrust_vp_lowercrust(i) = vp_lowercrust
    eucrust_mohodepth(i) = moho_depth
    eucrust_basement(i) = basement
    eucrust_ucdepth(i) = upper_lower_depth

  enddo
  close(11)

  end subroutine read_EuCrust

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_eucrust(lat,lon,x,vp,found_crust)

  use model_eucrust_par

  implicit none

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
  call eu_cap_smoothing(lat,lon,x,vp,found_crust)

  ! without smoothing
  !vp = crust_eu(lat,lon,x,vp,found_crust)

  end subroutine model_eucrust

!
!--------------------------------------------------------------------------------------------------
!

  double precision function crust_eu(lat,lon,x,vp,found_crust)

! returns Vp at the specific location lat/lon

  use model_eucrust_par

  implicit none

  include "constants.h"

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
    if( lon >= eucrust_lon(i) .and. lon < eucrust_lon(i+1) ) then
          do j=0,ilats-1
            if(lat>=eucrust_lat(i+j*ilons) .and. lat<eucrust_lat(i+(j+1)*ilons)) then

              h_basement = eucrust_basement(i+j*ilons)
              h_uc = eucrust_ucdepth(i+j*ilons)
              h_moho = eucrust_mohodepth(i+j*ilons)

              x3=(R_EARTH - h_basement*1000.0d0)/R_EARTH
              x4=(R_EARTH - h_uc*1000.0d0)/R_EARTH
              x5=(R_EARTH - h_moho*1000.0d0)/R_EARTH

              scaleval = dsqrt(PI*GRAV*RHOAV)

              if( x > x3 .and. INCLUDE_SEDIMENTS_CRUST &
                .and. h_basement > MINIMUM_SEDIMENT_THICKNESS) then
                ! above sediment basement, returns average upper crust value
                ! since no special sediment values are given
                found_crust = .true.
                vp = eucrust_vp_uppercrust(i+j*ilons) *1000.0d0/(R_EARTH*scaleval)
                crust_eu = vp
                return
              else if( x > x4 ) then
                found_crust = .true.
                vp = eucrust_vp_uppercrust(i+j*ilons) *1000.0d0/(R_EARTH*scaleval)
                crust_eu = vp
                return
              else if( x > x5 ) then
                found_crust = .true.
                vp = eucrust_vp_lowercrust(i+j*ilons) *1000.0d0/(R_EARTH*scaleval)
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
  subroutine eu_cap_smoothing(lat,lon,radius,value,found)

! smooths with a cap of size CAP (in degrees)
! using NTHETA points in the theta direction (latitudal)
! and NPHI in the phi direction (longitudal).
! The cap is rotated to the North Pole.

  use model_eucrust_par

  implicit none

  include "constants.h"

  ! argument variables
  double precision lat,lon,radius
  double precision :: value
  logical :: found

  integer, parameter :: NTHETA = 4
  integer, parameter :: NPHI = 10
  double precision, parameter :: CAP = 1.0d0 * DEGREES_TO_RADIANS  ! 1 degree smoothing

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
  !  value = func(lat,lon,x,value,found)
  !  return

  theta = (90.0-lat)*DEGREES_TO_RADIANS
  phi = lon*DEGREES_TO_RADIANS

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
  dphi = TWO_PI/dble(NPHI)
  cap_area = TWO_PI*(1.0-cos(CAP))

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
      xlat(i) = (PI_OVER_TWO-theta_rot)*RADIANS_TO_DEGREES
      xlon(i) = phi_rot*RADIANS_TO_DEGREES
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
    valuel = crust_eu(xlat(i),xlon(i),radius,value,found)
    value = value + weight(i)*valuel
  enddo

  if( abs(value) < TINYVAL) found = .false.

  end subroutine eu_cap_smoothing


