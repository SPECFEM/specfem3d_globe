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
!   This file contains subroutines to read in and get values for
!   3-D attenuation model QRFSI12 (Dalton, Ekstrom, & Dziewonski, 2008)
!
! C.A. Dalton, G. Ekstr\"om and A.M. Dziewonski, 2008.
! The global attenuation structure of the upper mantle,
! J. Geophys. Res., 113, B05317,10.1029/2006JB004394
!
!   Last edit: Colleen Dalton, March 25, 2008
!
! Q1: what are theta and phi?
! A1: input theta is colatitude in degrees, phi is longitude in degrees
!
! Q2: units for radius?
! A2: radius is given in km
!
! Q3: what to do about core?
! A3: more research :)
!
!--------------------------------------------------------------------------------------------------


  subroutine model_atten3D_QRFSI12_broadcast(myrank,QRFSI12_Q)

! standard routine to setup model

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  ! model_atten3D_QRFSI12_variables
  type model_atten3D_QRFSI12_variables
    sequence
    double precision dqmu(NKQ,NSQ)
    double precision spknt(NKQ)
    double precision refdepth(NDEPTHS_REFQ)
    double precision refqmu(NDEPTHS_REFQ)
  end type model_atten3D_QRFSI12_variables

  type (model_atten3D_QRFSI12_variables) QRFSI12_Q
  ! model_atten3D_QRFSI12_variables

  integer :: myrank

  ! local parameters
  integer :: ier

  if(myrank == 0) call read_atten_model_3D_QRFSI12(QRFSI12_Q)

  call MPI_BCAST(QRFSI12_Q%dqmu,          NKQ*NSQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QRFSI12_Q%spknt,             NKQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QRFSI12_Q%refdepth, NDEPTHS_REFQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QRFSI12_Q%refqmu,   NDEPTHS_REFQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)

  if(myrank == 0) write(IMAIN,*) 'read 3D attenuation model'


  end subroutine

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_atten_model_3D_QRFSI12(QRFSI12_Q)

  implicit none

  include "constants.h"

! three_d_model_atten3D_QRFSI12_variables
  type model_atten3D_QRFSI12_variables
    sequence
    double precision dqmu(NKQ,NSQ)
    double precision spknt(NKQ)
    double precision refdepth(NDEPTHS_REFQ)
    double precision refqmu(NDEPTHS_REFQ)
  end type model_atten3D_QRFSI12_variables

  type (model_atten3D_QRFSI12_variables) QRFSI12_Q
! three_d_model_atten3D_QRFSI12_variables

  integer j,k,l,m
  integer index,ll,mm,ier
  double precision v1,v2

  character(len=150) QRFSI12,QRFSI12_ref

! read in QRFSI12
! hard-wire for now
  QRFSI12='DATA/QRFSI12/QRFSI12.dat'
  QRFSI12_ref='DATA/QRFSI12/ref_QRFSI12'

! get the dq model coefficients
  open(unit=10,file=QRFSI12,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "', trim(QRFSI12), '": ', ier
    call exit_MPI(0, 'error model QRFSI12')
  endif

  do k=1,NKQ
    read(10,*)index
    j=0
    do l=0,MAXL_Q
      do m=0,l
        if(m==0)then
          j=j+1
          read(10,*)ll,mm,v1
          QRFSI12_Q%dqmu(k,j)=v1
        else
          j=j+2
          read(10,*)ll,mm,v1,v2
  !        write(*,*) 'k,l,m,ll,mm:',k,l,m,ll,mm,v1
          QRFSI12_Q%dqmu(k,j-1)=2.*v1
          QRFSI12_Q%dqmu(k,j)=-2.*v2
        endif
      enddo
    enddo
  enddo
  close(10)

! get the depths (km) of the spline knots
  QRFSI12_Q%spknt(1) = 24.4
  QRFSI12_Q%spknt(2) = 75.0
  QRFSI12_Q%spknt(3) = 150.0
  QRFSI12_Q%spknt(4) = 225.0
  QRFSI12_Q%spknt(5) = 300.0
  QRFSI12_Q%spknt(6) = 410.0
  QRFSI12_Q%spknt(7) = 530.0
  QRFSI12_Q%spknt(8) = 650.0

! get the depths and 1/Q values of the reference model
  open(11,file=QRFSI12_ref,status='old',action='read',iostat=ier)
  if ( ier /= 0 ) then
    write(IMAIN,*) 'error opening "', trim(QRFSI12_ref), '": ', ier
    call exit_MPI(0, 'error model QRFSI12')
  endif

  do j=1,NDEPTHS_REFQ
    read(11,*)QRFSI12_Q%refdepth(j),QRFSI12_Q%refqmu(j)
  enddo
  close(11)


  end subroutine read_atten_model_3D_QRFSI12

!
!----------------------------------
!

  subroutine model_atten3D_QRFSI12(radius,theta,phi,Qmu,QRFSI12_Q,idoubling)

  implicit none

  include "constants.h"

! model_atten3D_QRFSI12_variables
  type model_atten3D_QRFSI12_variables
    sequence
    double precision dqmu(NKQ,NSQ)
    double precision spknt(NKQ)
    double precision refdepth(NDEPTHS_REFQ)
    double precision refqmu(NDEPTHS_REFQ)
  end type model_atten3D_QRFSI12_variables

  type (model_atten3D_QRFSI12_variables) QRFSI12_Q
! model_atten3D_QRFSI12_variables

  double precision :: radius,theta,phi,Qmu
  integer :: idoubling

  ! local parameters
  integer :: i,j,k,n
  integer :: ifnd
  double precision :: smallq,dqmu,smallq_ref
  real(kind=4) :: splpts(NKQ),splcon(NKQ),splcond(NKQ)
  real(kind=4) :: depth,ylat,xlon
  real(kind=4) :: shdep(NSQ)
  real(kind=4) :: xlmvec(NSQ),wk1(NSQ),wk2(NSQ),wk3(NSQ)
  double precision, parameter :: rmoho_prem = 6371.0-24.4
  double precision, parameter :: rcmb = 3480.0

  ! in Colleen's original code theta refers to the latitude.  Here we have redefined theta to be colatitude
  ! to agree with the rest of specfem

  ! debug
  !  print *,'entering QRFSI12 subroutine'

  ylat=90.0d0-theta
  xlon=phi

  ! only checks radius for crust, idoubling is missleading for oceanic crust
  ! when we want to expand mantle up to surface...

!  !if(idoubling == IFLAG_CRUST .or. radius >= rmoho) then
  if( radius >= rmoho_prem ) then
  !   print *,'QRFSI12: we are in the crust'
     Qmu = 600.0d0
  else if(idoubling == IFLAG_INNER_CORE_NORMAL .or. idoubling == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
       idoubling == IFLAG_BOTTOM_CENTRAL_CUBE .or. idoubling == IFLAG_TOP_CENTRAL_CUBE .or. &
       idoubling == IFLAG_IN_FICTITIOUS_CUBE) then
    ! we are in the inner core

    !debug
    !   print *,'QRFSI12: we are in the inner core'

    Qmu = 84.6d0

  else if(idoubling == IFLAG_OUTER_CORE_NORMAL) then
  !   print *,'QRFSI12: we are in the outer core'
     Qmu = 0.0d0
  else !we are in the mantle
    depth = 6371.-radius
!   print *,'QRFSI12: we are in the mantle at depth',depth
    ifnd=0
    do i=2,NDEPTHS_REFQ
      if(depth >= QRFSI12_Q%refdepth(i-1) .and. depth < QRFSI12_Q%refdepth(i))then
        ifnd=i
      endif
    enddo
    if(ifnd == 0)then
      write(6,"('problem finding reference Q value at depth: ',f8.3)") depth
      stop
    endif
    smallq_ref=QRFSI12_Q%refqmu(ifnd)
    smallq = smallq_ref

    if(depth < 650.d0) then !Colleen's model is only defined between depths of 24.4 and 650km
      do j=1,NSQ
        shdep(j)=0.
      enddo
      do n=1,NKQ
        splpts(n)=QRFSI12_Q%spknt(n)
      enddo
      call vbspl(depth,NKQ,splpts,splcon,splcond)
      do n=1,NKQ
        do j=1,NSQ
          shdep(j)=shdep(j)+(splcon(n)*QRFSI12_Q%dqmu(n,j))
        enddo
      enddo
      call ylm(ylat,xlon,MAXL_Q,xlmvec,wk1,wk2,wk3)
      dqmu=0.
      do k=1,NSQ
        dqmu=dqmu+xlmvec(k)*shdep(k)
      enddo
      smallq = smallq_ref + dqmu
    endif

    ! if smallq is small and negative (due to numerical error), Qmu is very large:
    if(smallq < 0.0d0) smallq = 1.0d0/ATTENUATION_COMP_MAXIMUM

    Qmu = 1/smallq

    ! if Qmu is larger than MAX_ATTENUATION_VALUE, set it to ATTENUATION_COMP_MAXIMUM.
    ! This assumes that this value is high enough that at this point there is almost no attenuation at all.
    if(Qmu >= ATTENUATION_COMP_MAXIMUM) Qmu = 0.99d0*ATTENUATION_COMP_MAXIMUM

  endif

  end subroutine model_atten3D_QRFSI12
