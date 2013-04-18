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
!
! GAP P2 model - Global automatic parameterization model
!
! 3D Vp mantle model (version P2) from Masayuki Obayashi
!
!--------------------------------------------------------------------------------------------------


  module gapp2_mantle_model_constants
    ! data file resolution
    integer, parameter :: ma=288,mo=576,mr=32,mr1=64
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
  integer i,ir,ia,io,ier
  character(len=150) GAPP2

!...........................................input data

  ! default model: 3dvpGAP_P2
  call get_value_string(GAPP2, 'model.GAPP2', 'DATA/3dvpGAP_P2')

  ! reads in GAP-P2 model from Obayashi
  open(unit=10,file=GAPP2,status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(0,'error opening file for GAPP2 model')

  read(10,*) no,na,nnr,dela,delo

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) "model GAPP2: "
  write(IMAIN,*) "  dimensions no = ",no
  write(IMAIN,*) "             na,nnr = ",na,nnr
  write(IMAIN,*) "             dela,delon = ",dela,delo

  ! checks bounds
  if( nnr /= mr .or. no /= mo .or. na /= ma ) then
    print*,'error GAPP2 model bounds: '
    print*,'  file dimensions: nnr,no,na = ',nnr,no,na
    print*,'  module dimensions: mr,mo,ma = ',mr,mo,ma
    close(10)
    call exit_MPI(0,'please check GAPP2 model dimensions, and update model_gapp2.f90')
  endif

  read(10,*) (dep(i),i=0,nnr)
  read(10,*) nr1

  ! checks bounds
  write(IMAIN,*) "             nr1 = ",nr1
  if( nr1 /= mr1 ) then
    print*,'error GAPP2 model bounds: '
    print*,'  file dimensions: nr1 = ',nr1
    print*,'  module dimensions: mr1 = ',mr1
    close(10)
    call exit_MPI(0,'please check GAPP2 model dimensions, and update model_gapp2.f90')
  endif

  read(10,*) (dep1(i),i=0,nr1)
  read(10,*) (vp1(i),i=0,nr1)

  ! reads vp
  do ir=1,nnr
    do ia=1,na
      ! reads file 2 lines for all no values
      read(10,'(288f7.3)') (vp3(ia,io,ir),io=1,no)

      !read(10,*,iostat=ier) (vp3(ia,io,ir),io=1,no/2)
      !read(10,*,iostat=ier) (vp3(ia,io,ir),io=no/2,no)

      if( ier /= 0 ) then
        print*,'error GAPP2 read: ia,ir = ',ia,ir
        call exit_MPI(0,'error GAPP2 read')
      endif
    enddo
  enddo
  close(10)

  ! user output
  write(IMAIN,*) '  check vp3: ',vp3(1,1,1),vp3(na,no,nnr)
  write(IMAIN,*) '  check vp3: min/max = ',minval(vp3),maxval(vp3)
  write(IMAIN,*)

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
    if(icon/=0) then
       write(6,*)icon
       write(6,*) radius,theta,phi,dvp,dvs,drho
    endif

    ! latitude
    if(theta>=PI) then
       ia = na
    else
       ia = theta / dtheta + 1
    endif
    ! longitude
    if(phi < 0.0d0) phi = phi + 2.*PI
    io=phi / dphi + 1
    if(io>no) io=io-no

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
    if(d>dmax) then
       icon=99
    else if(d<dmin) then
       icon=-99
    else if(d==dmax) then
       id=mr+1
    else
       do i = 0, mr
          if(d<di(i)) then
             id=i
             goto 900
          endif
       enddo
    endif
900 continue

!..................................................................
    return

  end subroutine d2id
