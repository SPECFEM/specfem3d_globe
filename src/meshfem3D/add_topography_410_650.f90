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

  subroutine add_topography_410_650(myrank,xelm,yelm,zelm,R220,R400,R670,R771, &
    numker,numhpa,numcof,ihpa,lmax,nylm, &
    lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
    nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
    coe,ylmcof,wk1,wk2,wk3,varstr)

  implicit none

  include "constants.h"

  integer myrank

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision R220,R400,R670,R771

  integer ia

  real(kind=4) xcolat,xlon
  real(kind=4) topo410out,topo650out
  double precision topo410,topo650

  double precision r,theta,phi
  double precision gamma

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=40) varstr(maxker)

! we loop on all the points of the element
  do ia = 1,NGNOD

! convert to r theta phi
    call xyz_2_rthetaphi_dble(xelm(ia),yelm(ia),zelm(ia),r,theta,phi)
    call reduce(theta,phi)

! get colatitude and longitude in degrees
    xcolat = sngl(theta*180.0d0/PI)
    xlon = sngl(phi*180.0d0/PI)

! compute topography on 410 and 650 at current point
    call subtopo(xcolat,xlon,topo410out,topo650out, &
                 numker,numhpa,numcof,ihpa,lmax,nylm, &
                 lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
                 nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
                 coe,ylmcof,wk1,wk2,wk3,varstr)

! non-dimensionalize the topography, which is in km
! positive for a depression, so change the sign for a perturbation in radius
    topo410 = -dble(topo410out) / R_EARTH_KM
    topo650 = -dble(topo650out) / R_EARTH_KM

    gamma = 0.d0
    if(r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
! stretching between R220 and R400
      gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
      xelm(ia) = xelm(ia)*(ONE + gamma * topo410 / r)
      yelm(ia) = yelm(ia)*(ONE + gamma * topo410 / r)
      zelm(ia) = zelm(ia)*(ONE + gamma * topo410 / r)
    else if(r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
! stretching between R771 and R670
      gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
      xelm(ia) = xelm(ia)*(ONE + gamma * topo650 / r)
      yelm(ia) = yelm(ia)*(ONE + gamma * topo650 / r)
      zelm(ia) = zelm(ia)*(ONE + gamma * topo650 / r)
    else if(r > R670/R_EARTH .and. r < R400/R_EARTH) then
! stretching between R670 and R400
      gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
      xelm(ia) = xelm(ia)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      yelm(ia) = yelm(ia)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      zelm(ia) = zelm(ia)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
    endif
    if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

  enddo

  end subroutine add_topography_410_650

!
!-------------------------------------------------------------------------------------------------
!

  !> Hejun
  ! use GLL points to capture 410_650 topography
  ! JAN08, 2010
  subroutine add_topography_410_650_gll(myrank,xstore,ystore,zstore,ispec,nspec,R220,R400,R670,R771, &
        numker,numhpa,numcof,ihpa,lmax,nylm, &
        lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
        nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
        coe,ylmcof,wk1,wk2,wk3,varstr)

  implicit none

  include "constants.h"

  integer myrank
  integer:: ispec,nspec
  double precision,dimension(NGLLX,NGLLY,NGLLZ,nspec):: xstore,ystore,zstore

  double precision R220,R400,R670,R771

  integer i,j,k

  real(kind=4) xcolat,xlon
  real(kind=4) topo410out,topo650out
  double precision topo410,topo650

  double precision r,theta,phi
  double precision gamma

  integer, parameter :: maxker=200
  integer, parameter :: maxl=72
  integer, parameter :: maxcoe=2000
  integer, parameter :: maxver=1000
  integer, parameter :: maxhpa=2

  integer numker
  integer numhpa,numcof
  integer ihpa,lmax,nylm
  integer lmxhpa(maxhpa)
  integer itypehpa(maxhpa)
  integer ihpakern(maxker)
  integer numcoe(maxhpa)
  integer ivarkern(maxker)

  integer nconpt(maxhpa),iver
  integer iconpt(maxver,maxhpa)
  real(kind=4) conpt(maxver,maxhpa)

  real(kind=4) xlaspl(maxcoe,maxhpa)
  real(kind=4) xlospl(maxcoe,maxhpa)
  real(kind=4) radspl(maxcoe,maxhpa)
  real(kind=4) coe(maxcoe,maxker)

  real(kind=4) ylmcof((maxl+1)**2,maxhpa)
  real(kind=4) wk1(maxl+1)
  real(kind=4) wk2(maxl+1)
  real(kind=4) wk3(maxl+1)

  character(len=40) varstr(maxker)

  ! we loop on all GLL points of the element
  do k = 1,NGLLZ
     do j = 1,NGLLY
        do i = 1,NGLLX

        ! convert to r theta phi
        call xyz_2_rthetaphi_dble(xstore(i,j,k,ispec),ystore(i,j,k,ispec),zstore(i,j,k,ispec),r,theta,phi)
        call reduce(theta,phi)

        ! get colatitude and longitude in degrees
        xcolat = sngl(theta*180.0d0/PI)
        xlon = sngl(phi*180.0d0/PI)

        ! compute topography on 410 and 650 at current point
        call subtopo(xcolat,xlon,topo410out,topo650out, &
                 numker,numhpa,numcof,ihpa,lmax,nylm, &
                 lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
                 nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
                 coe,ylmcof,wk1,wk2,wk3,varstr)

        ! non-dimensionalize the topography, which is in km
        ! positive for a depression, so change the sign for a perturbation in radius
        topo410 = -dble(topo410out) / R_EARTH_KM
        topo650 = -dble(topo650out) / R_EARTH_KM

        gamma = 0.d0
        if(r >= R400/R_EARTH .and. r <= R220/R_EARTH) then
        ! stretching between R220 and R400
                gamma = (R220/R_EARTH - r) / (R220/R_EARTH - R400/R_EARTH)
                xstore(i,j,k,ispec) = xstore(i,j,k,ispec)*(ONE + gamma * topo410 / r)
                ystore(i,j,k,ispec) = ystore(i,j,k,ispec)*(ONE + gamma * topo410 / r)
                zstore(i,j,k,ispec) = zstore(i,j,k,ispec)*(ONE + gamma * topo410 / r)
        else if(r>= R771/R_EARTH .and. r <= R670/R_EARTH) then
        ! stretching between R771 and R670
                gamma = (r - R771/R_EARTH) / (R670/R_EARTH - R771/R_EARTH)
                xstore(i,j,k,ispec) = xstore(i,j,k,ispec)*(ONE + gamma * topo650 / r)
                ystore(i,j,k,ispec) = ystore(i,j,k,ispec)*(ONE + gamma * topo650 / r)
                zstore(i,j,k,ispec) = zstore(i,j,k,ispec)*(ONE + gamma * topo650 / r)
        else if(r > R670/R_EARTH .and. r < R400/R_EARTH) then
        ! stretching between R670 and R400
                gamma = (R400/R_EARTH - r) / (R400/R_EARTH - R670/R_EARTH)
                xstore(i,j,k,ispec) = xstore(i,j,k,ispec)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
                ystore(i,j,k,ispec) = ystore(i,j,k,ispec)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
                zstore(i,j,k,ispec) = zstore(i,j,k,ispec)*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
        endif
        if(gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

        enddo
     enddo
  enddo

  end subroutine add_topography_410_650_gll
