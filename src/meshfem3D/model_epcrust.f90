!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2011
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
! EPCRUST 1.0
! I. Molinari and A. Morelli, 2011.
! EPcrust: A reference crustal model for the European plate
! GJI, 185 (1), pages 352-364
!--------------------------------------------------------------------------------------------------


  subroutine model_epcrust_broadcast(myrank,EPCRUST)

  implicit none

  include "constants.h"
  include 'mpif.h'

  type model_epcrust_variables
    sequence
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT):: lon_ep,lat_ep,topo_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: thickness_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vp_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vs_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: rho_ep
  end type model_epcrust_variables
  type (model_epcrust_variables) EPCRUST
  integer:: myrank,ierr

  ! read EPCRUST model on master
  if(myrank == 0) call read_epcrust_model(EPCRUST)

  ! broadcast EPCRUST model
  call MPI_BCAST(EPCRUST%lon_ep,EPCRUST_NLON*EPCRUST_NLAT, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(EPCRUST%lat_ep,EPCRUST_NLON*EPCRUST_NLAT, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(EPCRUST%topo_ep,EPCRUST_NLON*EPCRUST_NLAT, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(EPCRUST%thickness_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(EPCRUST%vp_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(EPCRUST%vs_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_BCAST(EPCRUST%rho_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER, &
                MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

  end subroutine model_epcrust_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_epcrust_model(EPCRUST)

  implicit none
  include "constants.h"

  type model_epcrust_variables
    sequence
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT):: lon_ep,lat_ep,topo_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: thickness_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vp_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vs_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: rho_ep
  end type model_epcrust_variables
  type (model_epcrust_variables) EPCRUST

  character(len=150) EPCRUST_FNM
  character(len=150),dimension(15) :: header
  double precision,dimension(15) :: tmp
  integer:: ilon, jlat

  call get_value_string(EPCRUST_FNM,'model.EPCRUST_FNM',PATHNAME_EPCRUST)
  open(unit=1001,file=EPCRUST_FNM,status='old',action='read')
  read(1001,*) header

  do jlat = 1,EPCRUST_NLAT
    do ilon=1,EPCRUST_NLON
      read(1001,*) tmp

      EPCRUST%lon_ep(ilon,jlat) = tmp(1)
      EPCRUST%lat_ep(ilon,jlat) = tmp(2)
      EPCRUST%topo_ep(ilon,jlat) = tmp(3)
      EPCRUST%thickness_ep(ilon,jlat,1:3) = tmp(4:6)
      EPCRUST%vp_ep(ilon,jlat,1:3) = tmp(7:9)
      EPCRUST%vs_ep(ilon,jlat,1:3) = tmp(10:12)
      EPCRUST%rho_ep(ilon,jlat,1:3) = tmp(13:15)
    end do
  end do
  close(1001)

  end subroutine read_epcrust_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_epcrust(lat,lon,dep,vp,vs,rho,moho,found_crust,EPCRUST,elem_in_crust)
  implicit none
  include "constants.h"

  ! INPUT & OUTPUT
  type model_epcrust_variables
    sequence
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT):: lon_ep,lat_ep,topo_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: thickness_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vp_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: vs_ep
    double precision,dimension(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER):: rho_ep
  end type model_epcrust_variables
  type (model_epcrust_variables) EPCRUST

  double precision:: lat, lon, dep, vp, vs, rho, moho
  logical :: found_crust, elem_in_crust

  ! INTERIOR
  integer:: ilon, jlat, k
  double precision :: z0 , topo, basement, conrad, moho_top, scaleval
  double precision,dimension(3):: zsmooth, vpsmooth, vssmooth, rhosmooth
  double precision,dimension(NTHETA_EP*NPHI_EP) :: x1,y1,weight
  double precision:: weightl
  double precision:: cut, min_sed

  !if ( lat < EPCRUST_LAT_MIN .or. lat > EPCRUST_LAT_MAX &
  !        .or. lon < EPCRUST_LON_MIN .or. lon > EPCRUST_LON_MAX ) then
  !        stop 'incorrect enter EPCRUST model, check lat and lon'
  !end if

  vp=0.0d0
  vs=0.0d0
  rho=0.0d0

  if ( .not. flag_smooth_epcrust) then
    call ilon_jlat(lon,lat,ilon,jlat)
    z0=EPCRUST%topo_ep(ilon,jlat)
    zsmooth(:)=EPCRUST%thickness_ep(ilon,jlat,:)
    vpsmooth(:)=EPCRUST%vp_ep(ilon,jlat,:)
    vssmooth(:)=EPCRUST%vs_ep(ilon,jlat,:)
    rhosmooth(:)=EPCRUST%rho_ep(ilon,jlat,:)
  else
    call epcrust_smooth_base(lon,lat,x1,y1,weight)
    z0=0.d0
    zsmooth(:)=0.0d0
    vpsmooth(:)=0.0d0
    vssmooth(:)=0.0d0
    rhosmooth(:)=0.0d0

    do k = 1,NTHETA_EP*NPHI_EP
      call ilon_jlat(x1(k),y1(k),ilon,jlat)
      weightl=weight(k)

      z0=z0+weightl*EPCRUST%topo_ep(ilon,jlat)
      zsmooth(:)=zsmooth(:)+weightl*EPCRUST%thickness_ep(ilon,jlat,:)
      vpsmooth(:)=vpsmooth(:)+weightl*EPCRUST%vp_ep(ilon,jlat,:)
      vssmooth(:)=vssmooth(:)+weightl*EPCRUST%vs_ep(ilon,jlat,:)
      rhosmooth(:)=rhosmooth(:)+weightl*EPCRUST%rho_ep(ilon,jlat,:)
    end do
  end if

  !topo=(R_EARTH_KM+z0)/R_EARTH_KM
  !basement=(R_EARTH_KM+z0-zsmooth(1))/R_EARTH_KM
  !conrad=(R_EARTH_KM+z0-zsmooth(1)-zsmooth(2))/R_EARTH_KM
  !moho_top=(R_EARTH_KM+z0-zsmooth(1)-zsmooth(2)-zsmooth(3))/R_EARTH_KM

  topo=(R_EARTH_KM+z0)/R_EARTH_KM
  basement=(R_EARTH_KM-zsmooth(1))/R_EARTH_KM
  conrad=(R_EARTH_KM-zsmooth(1)-zsmooth(2))/R_EARTH_KM
  moho_top=(R_EARTH_KM-zsmooth(1)-zsmooth(2)-zsmooth(3))/R_EARTH_KM
  min_sed=1.0 - MINIMUM_SEDIMENT_THICKNESS/R_EARTH_KM

  found_crust=.true.

  if ( dep > basement .and. INCLUDE_SEDIMENTS_CRUST &
          .and. zsmooth(1) >= MINIMUM_SEDIMENT_THICKNESS ) then ! Hejun Zhu add minimum sediment thickness
    vp=vpsmooth(1)
    vs=vssmooth(1)
    rho=rhosmooth(1)
  else if ( dep > conrad ) then
    vp=vpsmooth(2)
    vs=vssmooth(2)
    rho=rhosmooth(2)
  else if ( dep > moho_top .or. elem_in_crust ) then
    vp=vpsmooth(3)
    vs=vssmooth(3)
    rho=rhosmooth(3)
  else
    found_crust=.false.
  end if

  if (found_crust ) then
    scaleval=dsqrt(PI*GRAV*RHOAV)
    vp=vp*1000.d0/(R_EARTH*scaleval)
    vs=vs*1000.d0/(R_EARTH*scaleval)
    rho=rho*1000.d0/RHOAV
  end if

  !moho = -(z0-zsmooth(1)-zsmooth(2)-zsmooth(3))/R_EARTH_KM
  moho = (zsmooth(1)+zsmooth(2)+zsmooth(3))/R_EARTH_KM

  ! Hejun Zhu, delete moho thickness less than 7 km
  cut=7.0/R_EARTH_KM
  if ( moho < cut ) then
    moho = cut
  end if

  end subroutine model_epcrust

!
!-------------------------------------------------------------------------------------------------
!


  subroutine epcrust_smooth_base(x,y,x1,y1,weight)

  implicit none
  include "constants.h"

  ! INPUT & OUTPUT
  double precision:: x, y
  double precision,dimension(NTHETA_EP*NPHI_EP):: x1,y1,weight

  ! INTERIOR
  double precision:: CAP,dtheta,dphi,cap_area,dweight,pi_over_nphi,total,wght
  double precision:: theta,phi,sint,cost,sinp,cosp
  double precision:: r_rot,theta_rot,phi_rot
  double precision,dimension(3,3):: rotation_matrix
  double precision,dimension(3):: xx,xc
  integer:: i,j,k,itheta,iphi
  double precision:: RADIANS_TO_DEGREES = 180.d0/PI
  double precision:: PI_OVER_TWO = PI/2.0d0

  x1(:)=0.0d0
  y1(:)=0.0d0
  weight(:)=0.0d0

  if (cap_degree_EP < TINYVAL ) then
          print*, 'error cap:', cap_degree_EP
          print*, 'lat/lon:', x,y
          stop 'error cap_degree too small'
  end if

  CAP=cap_degree_EP*PI/180.0d0
  dtheta=0.5d0*CAP/dble(NTHETA_EP)
  dphi=TWO_PI/dble(NPHI_EP)
  cap_area=TWO_PI*(1.0d0-dcos(CAP))
  dweight=CAP/dble(NTHETA_EP)*dphi/cap_area
  pi_over_nphi=PI/dble(NPHI_EP)

  phi=x*DEGREES_TO_RADIANS
  theta=(90.0d0-y)*DEGREES_TO_RADIANS

  sint=dsin(theta)
  cost=dcos(theta)
  sinp=dsin(phi)
  cosp=dcos(phi)

  rotation_matrix(1,1)=cosp*cost
  rotation_matrix(1,2)=-sinp
  rotation_matrix(1,3)=cosp*sint
  rotation_matrix(2,1)=sinp*cost
  rotation_matrix(2,2)=cosp
  rotation_matrix(2,3)=sinp*sint
  rotation_matrix(3,1)=-sint
  rotation_matrix(3,2)=0.0d0
  rotation_matrix(3,3)=cost

  i=0
  total=0.0d0
  do itheta=1,NTHETA_EP
    theta=dble(2*itheta-1)*dtheta
    cost=dcos(theta)
    sint=dsin(theta)
    wght=sint*dweight
    do iphi = 1,NPHI_EP
      i=i+1
      weight(i)=wght

      total=total+weight(i)
      phi=dble(2*iphi-1)*pi_over_nphi
      cosp=dcos(phi)
      sinp=dsin(phi)

      xc(1)=sint*cosp
      xc(2)=sint*sinp
      xc(3)=cost
      do j =1 ,3
              xx(j)=0.0d0
              do k = 1,3
                      xx(j)=xx(j)+rotation_matrix(j,k)*xc(k)
              end do
      end do
      call xyz_2_rthetaphi_dble(xx(1),xx(2),xx(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      x1(i)=phi_rot*RADIANS_TO_DEGREES
      y1(i)=(PI_OVER_TWO-theta_rot)*RADIANS_TO_DEGREES
      if (x1(i) > 180.d0) x1(i)=x1(i)-360.d0
    end do
  end do

  if (abs(total-1.0d0) > 0.001d0) then
    print*,'error cap:',total,cap_degree_EP
    stop
  end if

  end subroutine epcrust_smooth_base

!
!-------------------------------------------------------------------------------------------------
!

  subroutine ilon_jlat(lon,lat,ilon,jlat)

  implicit none
  include "constants.h"

  double precision:: lon,lat
  integer:: ilon,jlat

  ilon=int((lon-EPCRUST_LON_MIN)/EPCRUST_SAMPLE)+1
  jlat=int((EPCRUST_LAT_MAX-lat)/EPCRUST_SAMPLE)+1

  if (ilon < 1) ilon = 1
  if (ilon > EPCRUST_NLON) ilon = EPCRUST_NLON
  if (jlat < 1) jlat = 1
  if (jlat > EPCRUST_NLAT) jlat = EPCRUST_NLAT

  end subroutine ilon_jlat
