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

  module model_epcrust_par

  ! parameters for EPCRUST , from Molinari & Morelli model(2011)
  !       latitude :  9.0N - 89.5N
  !       longitude:  56.0W - 70.0E
  character(len=*), parameter :: PATHNAME_EPCRUST = 'DATA/epcrust/EPcrust_0_5.txt'

  integer, parameter :: EPCRUST_NLON = 253, EPCRUST_NLAT = 162, EPCRUST_NLAYER = 3
  double precision, parameter :: EPCRUST_LON_MIN = -56.0d0
  double precision, parameter :: EPCRUST_LON_MAX =  70.0d0
  double precision, parameter :: EPCRUST_LAT_MIN =   9.0d0
  double precision, parameter :: EPCRUST_LAT_MAX =  89.5d0
  double precision, parameter :: EPCRUST_SAMPLE = 0.5d0
  logical, parameter :: flag_smooth_epcrust = .true.
  integer, parameter :: NTHETA_EP = 4, NPHI_EP = 20
  double precision, parameter :: cap_degree_EP = 0.2d0

  ! arrays for EPCRUST 1.0
  double precision,dimension(:,:),allocatable :: lon_ep,lat_ep,topo_ep
  double precision,dimension(:,:,:),allocatable :: thickness_ep,vp_ep,vs_ep,rho_ep

  end module model_epcrust_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_epcrust_broadcast()

  use constants
  use model_epcrust_par

  implicit none

  integer :: ier

  ! allocates arrays for model
  allocate(lon_ep(EPCRUST_NLON,EPCRUST_NLAT), &
          lat_ep(EPCRUST_NLON,EPCRUST_NLAT), &
          topo_ep(EPCRUST_NLON,EPCRUST_NLAT), &
          thickness_ep(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER), &
          vp_ep(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER), &
          vs_ep(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER), &
          rho_ep(EPCRUST_NLON,EPCRUST_NLAT,EPCRUST_NLAYER), &
          stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating EPcrust arrays')

  ! read EPCRUST model on master
  if (myrank == 0) call read_epcrust_model()

  ! broadcast EPCRUST model
  call bcast_all_dp(lon_ep,EPCRUST_NLON*EPCRUST_NLAT)
  call bcast_all_dp(lat_ep,EPCRUST_NLON*EPCRUST_NLAT)
  call bcast_all_dp(topo_ep,EPCRUST_NLON*EPCRUST_NLAT)
  call bcast_all_dp(thickness_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER)
  call bcast_all_dp(vp_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER)
  call bcast_all_dp(vs_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER)
  call bcast_all_dp(rho_ep,EPCRUST_NLON*EPCRUST_NLAT*EPCRUST_NLAYER)

  end subroutine model_epcrust_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_epcrust_model()

  use constants
  use model_epcrust_par

  implicit none

  character(len=MAX_STRING_LEN),dimension(15) :: header
  double precision,dimension(15) :: tmp
  integer:: ilon,jlat,ier

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: EPcrust 1.0'
  write(IMAIN,*)

  open(unit=IIN,file=trim(PATHNAME_EPCRUST),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(PATHNAME_EPCRUST), '": ', ier
    call flush_IMAIN()
    ! stop
    call exit_MPI(0, 'Error model epcrust')
  endif

  read(IIN,*) header

  do jlat = 1,EPCRUST_NLAT
    do ilon = 1,EPCRUST_NLON
      read(IIN,*) tmp

      lon_ep(ilon,jlat) = tmp(1)
      lat_ep(ilon,jlat) = tmp(2)
      topo_ep(ilon,jlat) = tmp(3)
      thickness_ep(ilon,jlat,1:3) = tmp(4:6)
      vp_ep(ilon,jlat,1:3) = tmp(7:9)
      vs_ep(ilon,jlat,1:3) = tmp(10:12)
      rho_ep(ilon,jlat,1:3) = tmp(13:15)
    enddo
  enddo
  close(IIN)

  end subroutine read_epcrust_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_epcrust(lat,lon,dep,vp,vs,rho,moho,found_crust,elem_in_crust)

  use constants
  use model_epcrust_par

  implicit none

  ! INPUT & OUTPUT
  double precision,intent(in):: lat, lon, dep
  double precision,intent(out):: vp, vs, rho, moho
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust

  ! INTERIOR
  integer:: ilon, jlat, k
  double precision :: z0 , topo, basement, conrad, moho_top, scaleval
  double precision,dimension(3):: zsmooth, vpsmooth, vssmooth, rhosmooth
  double precision,dimension(NTHETA_EP*NPHI_EP) :: x1,y1,weight
  double precision:: weightl
  double precision:: cut, min_sed

  !if (lat < EPCRUST_LAT_MIN .or. lat > EPCRUST_LAT_MAX &
  ! .or. lon < EPCRUST_LON_MIN .or. lon > EPCRUST_LON_MAX) then
  !        stop 'incorrect enter EPCRUST model, check lat and lon'
  !endif

  vp = ZERO
  vs = ZERO
  rho = ZERO

  if (.not. flag_smooth_epcrust) then
    call ilon_jlat(lon,lat,ilon,jlat)
    z0 = topo_ep(ilon,jlat)
    zsmooth(:) = thickness_ep(ilon,jlat,:)
    vpsmooth(:) = vp_ep(ilon,jlat,:)
    vssmooth(:) = vs_ep(ilon,jlat,:)
    rhosmooth(:) = rho_ep(ilon,jlat,:)
  else
    call epcrust_smooth_base(lon,lat,x1,y1,weight)
    z0 = ZERO
    zsmooth(:) = ZERO
    vpsmooth(:) = ZERO
    vssmooth(:) = ZERO
    rhosmooth(:) = ZERO

    do k = 1,NTHETA_EP*NPHI_EP
      call ilon_jlat(x1(k),y1(k),ilon,jlat)
      weightl=weight(k)

      z0=z0+weightl*topo_ep(ilon,jlat)
      zsmooth(:)=zsmooth(:)+weightl*thickness_ep(ilon,jlat,:)
      vpsmooth(:)=vpsmooth(:)+weightl*vp_ep(ilon,jlat,:)
      vssmooth(:)=vssmooth(:)+weightl*vs_ep(ilon,jlat,:)
      rhosmooth(:)=rhosmooth(:)+weightl*rho_ep(ilon,jlat,:)
    enddo
  endif

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

  if (dep > basement .and. INCLUDE_SEDIMENTS_IN_CRUST &
          .and. zsmooth(1) >= MINIMUM_SEDIMENT_THICKNESS) then ! Hejun Zhu add minimum sediment thickness
    vp=vpsmooth(1)
    vs=vssmooth(1)
    rho=rhosmooth(1)
  else if (dep > conrad) then
    vp=vpsmooth(2)
    vs=vssmooth(2)
    rho=rhosmooth(2)
  else if (dep > moho_top .or. elem_in_crust) then
    vp=vpsmooth(3)
    vs=vssmooth(3)
    rho=rhosmooth(3)
  else
    found_crust=.false.
  endif

  if (found_crust) then
    scaleval=dsqrt(PI*GRAV*RHOAV)
    vp=vp*1000.d0/(R_EARTH*scaleval)
    vs=vs*1000.d0/(R_EARTH*scaleval)
    rho=rho*1000.d0/RHOAV
  endif

  !moho = -(z0-zsmooth(1)-zsmooth(2)-zsmooth(3))/R_EARTH_KM
  moho = (zsmooth(1)+zsmooth(2)+zsmooth(3))/R_EARTH_KM

  ! Hejun Zhu, delete moho thickness less than 7 km
  cut=7.0/R_EARTH_KM
  if (moho < cut) then
    moho = cut
  endif

  end subroutine model_epcrust

!
!-------------------------------------------------------------------------------------------------
!

  subroutine epcrust_smooth_base(x,y,x1,y1,weight)

  use constants
  use model_epcrust_par, only: NTHETA_EP,NPHI_EP,cap_degree_EP

  implicit none

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

  x1(:)=ZERO
  y1(:)=ZERO
  weight(:)=ZERO

  if (cap_degree_EP < TINYVAL) then
          print *, 'Error cap:', cap_degree_EP
          print *, 'lat/lon:', x,y
          stop 'Error cap_degree too small'
  endif

  CAP=cap_degree_EP * DEGREES_TO_RADIANS
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
  rotation_matrix(3,2)=ZERO
  rotation_matrix(3,3)=cost

  i = 0
  total=0.0d0
  do itheta = 1,NTHETA_EP
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
      do j = 1,3
              xx(j)=0.0d0
              do k = 1,3
                      xx(j)=xx(j)+rotation_matrix(j,k)*xc(k)
              enddo
      enddo
      call xyz_2_rthetaphi_dble(xx(1),xx(2),xx(3),r_rot,theta_rot,phi_rot)
      call reduce(theta_rot,phi_rot)
      x1(i)=phi_rot*RADIANS_TO_DEGREES
      y1(i)=(PI_OVER_TWO-theta_rot)*RADIANS_TO_DEGREES
      if (x1(i) > 180.d0) x1(i)=x1(i)-360.d0
    enddo
  enddo

  if (abs(total-1.0d0) > 0.001d0) then
    print *,'Error cap:',total,cap_degree_EP
    stop
  endif

  end subroutine epcrust_smooth_base

!
!-------------------------------------------------------------------------------------------------
!

  subroutine ilon_jlat(lon,lat,ilon,jlat)

  use constants
  use model_epcrust_par, only: &
    EPCRUST_LON_MIN,EPCRUST_LAT_MAX,EPCRUST_SAMPLE, &
    EPCRUST_NLON,EPCRUST_NLAT

  implicit none

  double precision:: lon,lat
  integer:: ilon,jlat

  ilon=int((lon-EPCRUST_LON_MIN)/EPCRUST_SAMPLE)+1
  jlat=int((EPCRUST_LAT_MAX-lat)/EPCRUST_SAMPLE)+1

  if (ilon < 1) ilon = 1
  if (ilon > EPCRUST_NLON) ilon = EPCRUST_NLON
  if (jlat < 1) jlat = 1
  if (jlat > EPCRUST_NLAT) jlat = EPCRUST_NLAT

  end subroutine ilon_jlat
