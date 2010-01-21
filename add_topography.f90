!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
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

  subroutine add_topography(myrank,xelm,yelm,zelm,ibathy_topo,R220)

  implicit none

  include "constants.h"

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  integer myrank

! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  integer ia

  double precision lat,lon,elevation,R220
  double precision r,theta,phi,colat
  double precision gamma

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
! slightly move points to avoid roundoff problem when exactly on the polar axis
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    theta = theta + 0.0000001d0
    phi = phi + 0.0000001d0
    call reduce(theta,phi)

! convert the geocentric colatitude to a geographic colatitude
  colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))

! get geographic latitude and longitude in degrees
  lat = 90.0d0 - colat*180.0d0/PI
  lon = phi*180.0d0/PI
  elevation = 0.d0

! compute elevation at current point
  call get_topo_bathy(lat,lon,elevation,ibathy_topo)

! non-dimensionalize the elevation, which is in meters
  elevation = elevation / R_EARTH

! stretching topography between d220 and the surface
  gamma = (r - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)

! add elevation to all the points of that element
! also make sure gamma makes sense
  if(gamma < -0.02 .or. gamma > 1.02) call exit_MPI(myrank,'incorrect value of gamma for topography')

  xelm(ia) = xelm(ia)*(ONE + gamma * elevation / r)
  yelm(ia) = yelm(ia)*(ONE + gamma * elevation / r)
  zelm(ia) = zelm(ia)*(ONE + gamma * elevation / r)

  enddo

  end subroutine add_topography

  !> Hejun
  ! This subroutine uses GLL points to capture topography variation rather
  ! than using control nodes
  ! Hejun Zhu, OCT16, 2009

  ! input parameters: myrank,
  !                   xstore,ystore,zstore,
  !                   ispec,nspec,
  !                   ibathy_topo
  !                   R220

  subroutine add_topography_gll(myrank,xstore,ystore,zstore,ispec,nspec,&
                                ibathy_topo,R220)

  implicit none

  include "constants.h"

  ! input parameters
  integer:: myrank
  integer:: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec):: xstore,ystore,zstore
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo
  double precision:: R220

  ! local parameters used in this subroutine                                      
  integer:: i,j,k
  double precision:: r,theta,phi,colat
  double precision:: lat,lon,elevation,gamma

  do k = 1,NGLLZ
     do j = 1,NGLLY
        do i = 1,NGLLX

           ! convert to r theta phi
           ! slightly move points to avoid roundoff problem when exactly on the polar axis
           call xyz_2_rthetaphi_dble(xstore(i,j,k,ispec),ystore(i,j,k,ispec),zstore(i,j,k,ispec),&
                                          r,theta,phi)
           theta = theta + 0.0000001d0
           phi = phi + 0.0000001d0
           call reduce(theta,phi)


           ! convert the geocentric colatitude to a geographic colatitude
           colat = PI/2.0d0 - datan(1.006760466d0*dcos(theta)/dmax1(TINYVAL,dsin(theta)))

           ! get geographic latitude and longitude in degrees
           lat = 90.0d0 - colat*180.0d0/PI
           lon = phi*180.0d0/PI
           elevation = 0.d0

           ! compute elevation at current point
           call get_topo_bathy(lat,lon,elevation,ibathy_topo)
           ! non-dimensionalize the elevation, which is in meters

           elevation = elevation / R_EARTH

           ! stretching topography between d220 and the surface
           gamma = (r - R220/R_EARTH) / (R_UNIT_SPHERE - R220/R_EARTH)
           !

           ! add elevation to all the points of that element
           ! also make sure factor makes sense
           if(gamma < -0.02 .or. gamma > 1.02) then
                call exit_MPI(myrank,'incorrect value of factor for topography gll points')
           end if 
           !

           ! since not all GLL points are exactlly at R220, use a small
           ! tolerance for R220 detection
           if (abs(gamma) < SMALLVAL) then
               gamma = 0.0
           end if 
           xstore(i,j,k,ispec) = xstore(i,j,k,ispec)*(ONE + gamma * elevation / r)
           ystore(i,j,k,ispec) = ystore(i,j,k,ispec)*(ONE + gamma * elevation / r)
           zstore(i,j,k,ispec) = zstore(i,j,k,ispec)*(ONE + gamma * elevation / r)

        end do 
     end do 
  end do
  end subroutine add_topography_gll
  !< Hejun
