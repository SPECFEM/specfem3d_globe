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
!
! GAP P2 model - Global automatic parameterization model
!
! 3D Vp mantle model (version P2) from Masayuki Obayashi
!
!--------------------------------------------------------------------------------------------------


  module gapp2_mantle_model_constants
    ! data file resolution
    integer, parameter :: ma=228,mo=576,mr=32,mr1=64
    integer no,na,nnr,nr1
    real dela,delo
    ! allocatable model arrays
    real,dimension(:),allocatable :: dep,dep1,vp1
    real,dimension(:,:,:),allocatable :: vp3
  end module gapp2_mantle_model_constants

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_gapp2_broadcast(myrank)

! standard routine to setup model

  use gapp2_mantle_model_constants

  implicit none
  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'
  integer :: myrank
  integer :: ier

  ! allocates arrays only when called and needed
  allocate(dep(0:mr),dep1(0:mr1),vp1(0:mr1),vp3(ma,mo,mr), &
          stat=ier)
  if( ier /= 0 ) then
    call exit_mpi(myrank,'error allocation GAP model')
  endif

  ! the variables read are declared in the module
  if(myrank == 0) call read_mantle_gapmodel()

  ! master process broadcasts data to all processes
  call MPI_BCAST( dep,mr+1,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(dep1,mr1+1,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST( vp1,mr1+1,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST( vp3,ma*mo*mr,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST( nnr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST( nr1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(  no,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(  na,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST( dela,1,MPI_REAL,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST( delo,1,MPI_REAL,0,MPI_COMM_WORLD,ier)

  end subroutine model_gapp2_broadcast

!
!--------------------------------------------------------------------------------------------------
!

  subroutine read_mantle_gapmodel()

  use gapp2_mantle_model_constants

  implicit none
  include "constants.h"
  integer i,ir,ia,io
  character(len=150) GAPP2

!...........................................input data

  ! default model: 3dvpGAP_P2
  call get_value_string(GAPP2, 'model.GAPP2', 'DATA/3dvpGAP_P2')

  ! reads in GAP-P2 model from Obayashi
  open(unit=10,file=GAPP2,status='old',action='read')

  read(10,'(3i4,2f10.6)') no,na,nnr,dela,delo
  read(10,'(34f8.2)') (dep(i),i=0,nnr)
  read(10,*) nr1
  read(10,'(67f8.2)') (dep1(i),i=0,nr1)
  read(10,'(67f8.3)') (vp1(i),i=0,nr1)
  do ir=1,nnr
    do ia=1,na
      read(10,'(256f7.3)') (vp3(ia,io,ir),io=1,no)
    enddo
  enddo
  write(6,*) vp3(1,1,1),vp3(na,no,nnr)
  close(10)

  end subroutine read_mantle_gapmodel

!
!--------------------------------------------------------------------------------------------------
!

  subroutine mantle_gapmodel(radius,theta,phi,dvs,dvp,drho)

    use gapp2_mantle_model_constants

    implicit none
    include "constants.h"
    integer id,ia,io,icon
    real d,dtheta,dphi

    double precision radius,theta,phi,dvs,dvp,drho

! factor to convert perturbations in shear speed to perturbations in density
    double precision, parameter :: SCALE_VS =  1.40d0
    double precision, parameter :: SCALE_RHO = 0.0d0

    double precision, parameter :: R_EARTH_ = 6371.d0
    double precision, parameter :: ZERO_ = 0.d0

!.....................................

    dvs = ZERO_
    dvp = ZERO_
    drho = ZERO_

    ! increments in latitude/longitude (in rad)
    dtheta = dela * PI / 180.0
    dphi = delo * PI / 180.0

    ! depth given in km
    d=R_EARTH_-radius*R_EARTH_

    call d2id(d,nnr,dep,id,icon)
    if(icon.ne.0) then
       write(6,*)icon
       write(6,*) radius,theta,phi,dvp,dvs,drho
    endif

    ! latitude
    if(theta.ge.PI) then
       ia = na
    else
       ia = theta / dtheta + 1
    endif
    ! longitude
    if(phi .lt. 0.0d0) phi = phi + 2.*PI
    io=phi / dphi + 1
    if(io.gt.no) io=io-no

    ! velocity and density perturbations
    dvp = vp3(ia,io,id)/100.d0
    dvs = SCALE_VS*dvp
    drho = SCALE_RHO*dvs

  end subroutine mantle_gapmodel

!
!--------------------------------------------------------------------------------------------------
!

  subroutine d2id(d,mr,di,id,icon)
!.................................................................
!     radial section index for a given depth d
!.................................................................
!   d     i   depth(km)
!   mr    i   number of radial division
!   di    i   depth table
!   id    o   depth section index for d
!              shallow .... di(id-1) <= d < di(id) .... deep
!   icon  o   condition code
!              0:normal, -99:above the surface, 99:below the cmb
!.................................................................
    integer i, mr, id, icon
    real d,dmax,dmin
    real di(0:mr)
    icon=0
    dmax=di(mr)
    dmin=di(0)
    if(d.gt.dmax) then
       icon=99
    else if(d.lt.dmin) then
       icon=-99
    else if(d.eq.dmax) then
       id=mr+1
    else
       do i = 0, mr
          if(d.lt.di(i)) then
             id=i
             goto 900
          endif
       enddo
    end if
900 continue

!..................................................................
    return

  end subroutine d2id
