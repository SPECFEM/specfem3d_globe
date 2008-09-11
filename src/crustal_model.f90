!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
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

!
! read or evaluate smoothed crust2.0 model
!

  subroutine get_smoothed_crust(lat,lon,x,vp,vs,rho,moho,found_crust,CM_V)

  implicit none

  include "constants.h"

  logical found_crust
  double precision lat,lon,x,vp,vs,rho,moho

  double precision h_sed,h_uc
  double precision x3,x4,x5,x6,x7,scaleval
  double precision vps(3:7),vss(3:7),rhos(3:7),thicks(3:7)

  double precision vps1(3:7),vss1(3:7),rhos1(3:7),thicks1(3:7)
  double precision vps2(3:7),vss2(3:7),rhos2(3:7),thicks2(3:7)
  double precision vps3(3:7),vss3(3:7),rhos3(3:7),thicks3(3:7)
  double precision vps4(3:7),vss4(3:7),rhos4(3:7),thicks4(3:7)

! crustal_model_variables
  type crustal_model_variables
    sequence
    real(kind=4) velocp(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) velocs(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) dens(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) thlr(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

  integer ilon_scaled,ilat_scaled
  double precision lon_scaled,lat_scaled
  double precision gamma_lon,gamma_lat

! scale interval and make all values positive
  lon_scaled = lon * NFACTOR_CRUST + NLON_CRUST
  lat_scaled = lat * NFACTOR_CRUST + NLAT_CRUST

! select the right cell to use in smoothed file
  ilon_scaled = int(lon_scaled)
  ilat_scaled = int(lat_scaled)

! compute coordinates in interpolation cell
  gamma_lon = lon_scaled - dble(ilon_scaled)
  gamma_lat = lat_scaled - dble(ilat_scaled)

! prevent edge effects
  if(ilon_scaled < 0) then
    ilon_scaled = 0
    gamma_lon = 0.d0
  endif

  if(ilat_scaled < 0) then
    ilat_scaled = 0
    gamma_lat = 0.d0
  endif

  if(ilon_scaled >= 2*NLON_CRUST) then
    ilon_scaled = 2*NLON_CRUST - 1
    gamma_lon = 1.d0
  endif

  if(ilat_scaled >= 2*NLAT_CRUST) then
    ilat_scaled = 2*NLAT_CRUST - 1
    gamma_lat = 1.d0
  endif

! copy four corners of interpolation cell
  vps1(:) = dble(CM_V%velocp(ilon_scaled,ilat_scaled,:))
  vss1(:) = dble(CM_V%velocs(ilon_scaled,ilat_scaled,:))
  rhos1(:) = dble(CM_V%dens(ilon_scaled,ilat_scaled,:))
  thicks1(:) = dble(CM_V%thlr(ilon_scaled,ilat_scaled,:))

  vps2(:) = dble(CM_V%velocp(ilon_scaled+1,ilat_scaled,:))
  vss2(:) = dble(CM_V%velocs(ilon_scaled+1,ilat_scaled,:))
  rhos2(:) = dble(CM_V%dens(ilon_scaled+1,ilat_scaled,:))
  thicks2(:) = dble(CM_V%thlr(ilon_scaled+1,ilat_scaled,:))

  vps3(:) = dble(CM_V%velocp(ilon_scaled+1,ilat_scaled+1,:))
  vss3(:) = dble(CM_V%velocs(ilon_scaled+1,ilat_scaled+1,:))
  rhos3(:) = dble(CM_V%dens(ilon_scaled+1,ilat_scaled+1,:))
  thicks3(:) = dble(CM_V%thlr(ilon_scaled+1,ilat_scaled+1,:))

  vps4(:) = dble(CM_V%velocp(ilon_scaled,ilat_scaled+1,:))
  vss4(:) = dble(CM_V%velocs(ilon_scaled,ilat_scaled+1,:))
  rhos4(:) = dble(CM_V%dens(ilon_scaled,ilat_scaled+1,:))
  thicks4(:) = dble(CM_V%thlr(ilon_scaled,ilat_scaled+1,:))

! perform interpolation from the four corners
  vps(:) = (1.d0-gamma_lon)*(1.d0-gamma_lat)*vps1(:) + &
           gamma_lon*(1.d0-gamma_lat)*vps2(:) + &
           gamma_lon*gamma_lat*vps3(:) + &
           (1.d0-gamma_lon)*gamma_lat*vps4(:)

  vss(:) = (1.d0-gamma_lon)*(1.d0-gamma_lat)*vss1(:) + &
           gamma_lon*(1.d0-gamma_lat)*vss2(:) + &
           gamma_lon*gamma_lat*vss3(:) + &
           (1.d0-gamma_lon)*gamma_lat*vss4(:)

  rhos(:) = (1.d0-gamma_lon)*(1.d0-gamma_lat)*rhos1(:) + &
           gamma_lon*(1.d0-gamma_lat)*rhos2(:) + &
           gamma_lon*gamma_lat*rhos3(:) + &
           (1.d0-gamma_lon)*gamma_lat*rhos4(:)

  thicks(:) = (1.d0-gamma_lon)*(1.d0-gamma_lat)*thicks1(:) + &
           gamma_lon*(1.d0-gamma_lat)*thicks2(:) + &
           gamma_lon*gamma_lat*thicks3(:) + &
           (1.d0-gamma_lon)*gamma_lat*thicks4(:)

 x3 = (R_EARTH-thicks(3)*1000.0d0)/R_EARTH
 h_sed = thicks(3) + thicks(4)
 x4 = (R_EARTH-h_sed*1000.0d0)/R_EARTH
 h_uc = h_sed + thicks(5)
 x5 = (R_EARTH-h_uc*1000.0d0)/R_EARTH
 x6 = (R_EARTH-(h_uc+thicks(6))*1000.0d0)/R_EARTH
 x7 = (R_EARTH-(h_uc+thicks(6)+thicks(7))*1000.0d0)/R_EARTH

 found_crust = .true.

 if(x > x3 .and. INCLUDE_SEDIMENTS_CRUST) then
   vp = vps(3)
   vs = vss(3)
   rho = rhos(3)
 else if(x > x4 .and. INCLUDE_SEDIMENTS_CRUST) then
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
 else if(x > x7) then
   vp = vps(7)
   vs = vss(7)
   rho = rhos(7)
 else
   found_crust = .false.
 endif

  if(found_crust) then

!   non-dimensionalize
    scaleval = dsqrt(PI*GRAV*RHOAV)
    vp = vp*1000.0d0/(R_EARTH*scaleval)
    vs = vs*1000.0d0/(R_EARTH*scaleval)
    rho = rho*1000.0d0/RHOAV
    moho = (h_uc+thicks(6)+thicks(7))*1000.0d0/R_EARTH

  endif

 end subroutine get_smoothed_crust

!---------------------------

  subroutine read_smoothed_3D_crustal_model(CM_V)

  implicit none

  include "constants.h"

! crustal_model_variables
  type crustal_model_variables
    sequence
    real(kind=4) velocp(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) velocs(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) dens(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    real(kind=4) thlr(0:2*NLON_CRUST,0:2*NLAT_CRUST,3:7)
    end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

  integer ilon,ilat

  open(unit=1,file='DATA/crust2.0/smoothed_crust2.0.dat',status='old',action='read')

! loop on latitude and longitude
  do ilat = 0,2*NLAT_CRUST
  do ilon = 0,2*NLON_CRUST

! we only saved layers 3 to 7 because 1, 2 and 8 are never used in the mesher

  read(1,*) CM_V%velocp(ilon,ilat,3)
  read(1,*) CM_V%velocp(ilon,ilat,4)
  read(1,*) CM_V%velocp(ilon,ilat,5)
  read(1,*) CM_V%velocp(ilon,ilat,6)
  read(1,*) CM_V%velocp(ilon,ilat,7)

  read(1,*) CM_V%velocs(ilon,ilat,3)
  read(1,*) CM_V%velocs(ilon,ilat,4)
  read(1,*) CM_V%velocs(ilon,ilat,5)
  read(1,*) CM_V%velocs(ilon,ilat,6)
  read(1,*) CM_V%velocs(ilon,ilat,7)

  read(1,*) CM_V%dens(ilon,ilat,3)
  read(1,*) CM_V%dens(ilon,ilat,4)
  read(1,*) CM_V%dens(ilon,ilat,5)
  read(1,*) CM_V%dens(ilon,ilat,6)
  read(1,*) CM_V%dens(ilon,ilat,7)

  read(1,*) CM_V%thlr(ilon,ilat,3)
  read(1,*) CM_V%thlr(ilon,ilat,4)
  read(1,*) CM_V%thlr(ilon,ilat,5)
  read(1,*) CM_V%thlr(ilon,ilat,6)
  read(1,*) CM_V%thlr(ilon,ilat,7)

  enddo
  enddo

  close(1)

  end subroutine read_smoothed_3D_crustal_model

