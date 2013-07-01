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

  module model_atten3D_QRFSI12_par

  ! QRFSI12 constants
  integer,parameter :: NKQ=8,MAXL_Q=12
  integer,parameter :: NSQ=(MAXL_Q+1)**2,NDEPTHS_REFQ=913

  ! model_atten3D_QRFSI12_variables
  double precision,dimension(:,:),allocatable :: QRFSI12_Q_dqmu
  double precision,dimension(:),allocatable :: QRFSI12_Q_spknt
  double precision,dimension(:),allocatable :: QRFSI12_Q_refdepth,QRFSI12_Q_refqmu

  end module model_atten3D_QRFSI12_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_atten3D_QRFSI12_broadcast(myrank)

! standard routine to setup model

  use model_atten3D_QRFSI12_par

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  integer :: myrank

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(QRFSI12_Q_dqmu(NKQ,NSQ), &
          QRFSI12_Q_spknt(NKQ), &
          QRFSI12_Q_refdepth(NDEPTHS_REFQ), &
          QRFSI12_Q_refqmu(NDEPTHS_REFQ), &
          stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating QRFSI12_Q model arrays')

  ! master process reads in file values
  if(myrank == 0) call read_atten_model_3D_QRFSI12()

  ! broadcasts to all processes
  call MPI_BCAST(QRFSI12_Q_dqmu,          NKQ*NSQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QRFSI12_Q_spknt,             NKQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QRFSI12_Q_refdepth, NDEPTHS_REFQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)
  call MPI_BCAST(QRFSI12_Q_refqmu,   NDEPTHS_REFQ,MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ier)

  ! user output
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'read 3D attenuation model'
    write(IMAIN,*) '  model: QRFSI12'
    write(IMAIN,*)
  endif

  end subroutine

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_atten_model_3D_QRFSI12()

  use model_atten3D_QRFSI12_par

  implicit none

  include "constants.h"

  ! local parameters
  integer :: j,k,l,m,ier
  integer :: index,ll,mm
  double precision :: v1,v2

  character(len=150) :: QRFSI12,QRFSI12_ref

! read in QRFSI12
! hard-wire for now
  QRFSI12='DATA/QRFSI12/QRFSI12.dat'
  QRFSI12_ref='DATA/QRFSI12/ref_QRFSI12'

! get the dq model coefficients
  open(unit=10,file=QRFSI12,status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(0,'error opening model file QRFSI12.dat')

  do k=1,NKQ
    read(10,*)index
    j=0
    do l=0,MAXL_Q
      do m=0,l
        if(m==0)then
          j=j+1
          read(10,*)ll,mm,v1
          QRFSI12_Q_dqmu(k,j)=v1
        else
          j=j+2
          read(10,*)ll,mm,v1,v2
  !        write(*,*) 'k,l,m,ll,mm:',k,l,m,ll,mm,v1
          QRFSI12_Q_dqmu(k,j-1)=2.*v1
          QRFSI12_Q_dqmu(k,j)=-2.*v2
        endif
      enddo
    enddo
  enddo
  close(10)

! get the depths (km) of the spline knots
  QRFSI12_Q_spknt(1) = 24.4d0
  QRFSI12_Q_spknt(2) = 75.d0
  QRFSI12_Q_spknt(3) = 150.d0
  QRFSI12_Q_spknt(4) = 225.d0
  QRFSI12_Q_spknt(5) = 300.d0
  QRFSI12_Q_spknt(6) = 410.d0
  QRFSI12_Q_spknt(7) = 530.d0
  QRFSI12_Q_spknt(8) = 650.d0

! get the depths and 1/Q values of the reference model
  open(11,file=QRFSI12_ref,status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(0,'error opening model file ref_QRFSI12')

  do j=1,NDEPTHS_REFQ
    read(11,*)QRFSI12_Q_refdepth(j),QRFSI12_Q_refqmu(j)
  enddo
  close(11)


  end subroutine read_atten_model_3D_QRFSI12

!
!----------------------------------
!

  subroutine model_atten3D_QRFSI12(radius,theta,phi,Qmu,idoubling)

  use model_atten3D_QRFSI12_par

  implicit none

  include "constants.h"

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

    ! we are in the outer core

    !debug
    !   print *,'QRFSI12: we are in the outer core'

    Qmu = 0.0d0

  else

    ! we are in the mantle

    depth = 6371.-radius

    !debug
    !   print *,'QRFSI12: we are in the mantle at depth',depth

    ifnd=0
    do i=2,NDEPTHS_REFQ
      if(depth >= QRFSI12_Q_refdepth(i-1) .and. depth < QRFSI12_Q_refdepth(i))then
        ifnd=i
      endif
    enddo
    if(ifnd == 0)then
      write(6,"('problem finding reference Q value at depth: ',f8.3)") depth
      stop
    endif
    smallq_ref=QRFSI12_Q_refqmu(ifnd)
    smallq = smallq_ref

    if(depth < 650.d0) then !Colleen's model is only defined between depths of 24.4 and 650km
      do j=1,NSQ
        shdep(j)=0.
      enddo
      do n=1,NKQ
        splpts(n)=QRFSI12_Q_spknt(n)
      enddo
      call vbspl(depth,NKQ,splpts,splcon,splcond)
      do n=1,NKQ
        do j=1,NSQ
          shdep(j)=shdep(j)+(splcon(n)*QRFSI12_Q_dqmu(n,j))
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

!
!----------------------------------
!

!!$  subroutine vbspl(x,np,xarr,splcon,splcond)
!!$!
!!$!---- this subroutine returns the spline contributions at a particular value of x
!!$!
!!$  implicit none
!!$
!!$  integer :: np
!!$
!!$  real(kind=4) :: xarr(np),x
!!$  real(kind=4) :: splcon(np)
!!$  real(kind=4) :: splcond(np)
!!$
!!$  real(kind=4) :: r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12,r13
!!$  real(kind=4) :: r1d,r2d,r3d,r4d,r5d,r6d,r7d,r8d,r9d,r10d,r11d,r12d,r13d,val,vald
!!$
!!$  real(kind=4) :: rr1,rr2,rr3,rr4,rr5,rr6,rr7,rr8,rr9,rr10,rr11,rr12
!!$  real(kind=4) :: rr1d,rr2d,rr3d,rr4d,rr5d,rr6d,rr7d,rr8d,rr9d,rr10d,rr11d,rr12d
!!$
!!$  integer :: iflag,interval,ik,ib
!!$
!!$!
!!$!---- iflag=1 ==>> second derivative is 0 at end points
!!$!---- iflag=0 ==>> first derivative is 0 at end points
!!$!
!!$  iflag=1
!!$!
!!$!---- first, find out within which interval x falls
!!$!
!!$  interval=0
!!$  ik=1
!!$  do while(interval == 0.and.ik < np)
!!$  ik=ik+1
!!$  if(x >= xarr(ik-1).and.x <= xarr(ik)) interval=ik-1
!!$  enddo
!!$  if(x > xarr(np)) then
!!$  interval=np
!!$  endif
!!$
!!$  if(interval == 0) then
!!$!        write(6,"('low value:',2f10.3)") x,xarr(1)
!!$  else if(interval > 0.and.interval < np) then
!!$!        write(6,"('bracket:',i5,3f10.3)") interval,xarr(interval),x,xarr(interval+1)
!!$  else
!!$!        write(6,"('high value:',2f10.3)") xarr(np),x
!!$  endif
!!$
!!$  do ib=1,np
!!$  val=0.
!!$  vald=0.
!!$  if(ib == 1) then
!!$
!!$    r1=(x-xarr(1))/(xarr(2)-xarr(1))
!!$    r2=(xarr(3)-x)/(xarr(3)-xarr(1))
!!$    r4=(xarr(2)-x)/(xarr(2)-xarr(1))
!!$    r5=(x-xarr(1))/(xarr(2)-xarr(1))
!!$    r6=(xarr(3)-x)/(xarr(3)-xarr(1))
!!$   r10=(xarr(2)-x)/(xarr(2)-xarr(1))
!!$   r11=(x-xarr(1))  /(xarr(2)-xarr(1))
!!$   r12=(xarr(3)-x)/(xarr(3)-xarr(2))
!!$   r13=(xarr(2)-x)/(xarr(2)-xarr(1))
!!$
!!$    r1d=1./(xarr(2)-xarr(1))
!!$    r2d=-1./(xarr(3)-xarr(1))
!!$    r4d=-1./(xarr(2)-xarr(1))
!!$    r5d=1./(xarr(2)-xarr(1))
!!$    r6d=-1./(xarr(3)-xarr(1))
!!$   r10d=-1./(xarr(2)-xarr(1))
!!$   r11d=1./(xarr(2)-xarr(1))
!!$   r12d=-1./(xarr(3)-xarr(2))
!!$   r13d=-1./(xarr(2)-xarr(1))
!!$
!!$    if(interval == ib.or.interval == 0) then
!!$         if(iflag == 0) then
!!$           val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11 +r13**3
!!$           vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
!!$           vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
!!$           vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
!!$           vald=vald+3.*r13d*r13**2
!!$         else if(iflag == 1) then
!!$           val=0.6667*(r1*r4*r10 + r2*r5*r10 + r2*r6*r11 &
!!$                    + 1.5*r13**3)
!!$           vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
!!$           vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
!!$           vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
!!$           vald=vald+4.5*r13d*r13**2
!!$           vald=0.6667*vald
!!$         endif
!!$    else if(interval == ib+1) then
!!$         if(iflag == 0) then
!!$           val=r2*r6*r12
!!$           vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
!!$         else if(iflag == 1) then
!!$           val=0.6667*r2*r6*r12
!!$           vald=0.6667*(r2d*r6*r12+r2*r6d*r12+r2*r6*r12d)
!!$         endif
!!$    else
!!$      val=0.
!!$    endif
!!$
!!$  else if(ib == 2) then
!!$
!!$    rr1=(x-xarr(1))/(xarr(2)-xarr(1))
!!$    rr2=(xarr(3)-x)/(xarr(3)-xarr(1))
!!$    rr4=(xarr(2)-x)/(xarr(2)-xarr(1))
!!$    rr5=(x-xarr(1))/(xarr(2)-xarr(1))
!!$    rr6=(xarr(3)-x)/(xarr(3)-xarr(1))
!!$   rr10=(xarr(2)-x)/(xarr(2)-xarr(1))
!!$   rr11=(x-xarr(1))  /(xarr(2)-xarr(1))
!!$   rr12=(xarr(3)-x)/(xarr(3)-xarr(2))
!!$
!!$    rr1d=1./(xarr(2)-xarr(1))
!!$    rr2d=-1./(xarr(3)-xarr(1))
!!$    rr4d=-1./(xarr(2)-xarr(1))
!!$    rr5d=1./(xarr(2)-xarr(1))
!!$    rr6d=-1./(xarr(3)-xarr(1))
!!$   rr10d=-1./(xarr(2)-xarr(1))
!!$   rr11d=1./(xarr(2)-xarr(1))
!!$   rr12d=-1./(xarr(3)-xarr(2))
!!$
!!$    r1=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
!!$    r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
!!$    r3=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
!!$    r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
!!$    r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
!!$    r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
!!$    r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
!!$    r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
!!$   r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
!!$   r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
!!$   r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))
!!$
!!$    r1d=1./(xarr(ib+1)-xarr(ib-1))
!!$    r2d=-1./(xarr(ib+2)-xarr(ib-1))
!!$    r3d=1./(xarr(ib)-xarr(ib-1))
!!$    r4d=-1./(xarr(ib+1)-xarr(ib-1))
!!$    r5d=1./(xarr(ib+1)-xarr(ib-1))
!!$    r6d=-1./(xarr(ib+2)-xarr(ib))
!!$    r8d=-1./  (xarr(ib)-xarr(ib-1))
!!$    r9d=1./(xarr(ib)-xarr(ib-1))
!!$   r10d=-1./(xarr(ib+1)-xarr(ib))
!!$   r11d=1./(xarr(ib+1)-xarr(ib))
!!$   r12d=-1./(xarr(ib+2)-xarr(ib+1))
!!$
!!$    if(interval == ib-1.or.interval == 0) then
!!$         val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
!!$         vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
!!$         vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
!!$         vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
!!$         if(iflag == 1) then
!!$           val=val+0.3333*(rr1*rr4*rr10 + rr2*rr5*rr10 + &
!!$                     rr2*rr6*rr11)
!!$           vald=vald+0.3333*(rr1d*rr4*rr10+rr1*rr4d*rr10+ &
!!$                    rr1*rr4*rr10d)
!!$           vald=vald+0.3333*(rr2d*rr5*rr10+rr2*rr5d*rr10+ &
!!$                    rr2*rr5*rr10d)
!!$           vald=vald+0.3333*(rr2d*rr6*rr11+rr2*rr6d*rr11+ &
!!$                    rr2*rr6*rr11d)
!!$         endif
!!$    else if(interval == ib) then
!!$         val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
!!$         vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
!!$         vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
!!$         vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
!!$         if(iflag == 1) then
!!$           val=val+0.3333*rr2*rr6*rr12
!!$           vald=vald+0.3333*(rr2d*rr6*rr12+rr2*rr6d*rr12+ &
!!$                    rr2*rr6*rr12d)
!!$         endif
!!$    else if(interval == ib+1) then
!!$         val=r2*r6*r12
!!$         vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
!!$    else
!!$         val=0.
!!$    endif
!!$  else if(ib == np-1) then
!!$
!!$    rr1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
!!$    rr2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
!!$    rr3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
!!$    rr4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
!!$    rr5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
!!$    rr7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
!!$    rr8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
!!$    rr9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
!!$
!!$    rr1d=1./(xarr(np)-xarr(np-2))
!!$    rr2d=-1./(xarr(np)-xarr(np-1))
!!$    rr3d=1./(xarr(np)-xarr(np-2))
!!$    rr4d=-1./(xarr(np)-xarr(np-1))
!!$    rr5d=1./(xarr(np)-xarr(np-1))
!!$    rr7d=1./(xarr(np-1)-xarr(np-2))
!!$    rr8d=-1./  (xarr(np)-xarr(np-1))
!!$    rr9d=1./(xarr(np)-xarr(np-1))
!!$
!!$    r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
!!$    r2=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
!!$    r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
!!$    r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
!!$    r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
!!$    r6=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
!!$    r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
!!$    r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
!!$    r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
!!$   r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
!!$   r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
!!$
!!$    r1d=1./(xarr(ib+1)-xarr(ib-2))
!!$    r2d=-1./(xarr(ib+1)-xarr(ib-1))
!!$    r3d=1./(xarr(ib)-xarr(ib-2))
!!$    r4d=-1./(xarr(ib+1)-xarr(ib-1))
!!$    r5d=1./(xarr(ib+1)-xarr(ib-1))
!!$    r6d=-1./(xarr(ib+1)-xarr(ib))
!!$    r7d=1./(xarr(ib-1)-xarr(ib-2))
!!$    r8d=-1./(xarr(ib)-xarr(ib-1))
!!$    r9d=1./(xarr(ib)-xarr(ib-1))
!!$   r10d=-1./(xarr(ib+1)-xarr(ib))
!!$   r11d=1./(xarr(ib+1)-xarr(ib))
!!$
!!$    if(interval == ib-2) then
!!$         val=r1*r3*r7
!!$         vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
!!$    else if(interval == ib-1) then
!!$         val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
!!$         vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
!!$         vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
!!$         vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
!!$         if(iflag == 1) then
!!$           val=val+0.3333*rr1*rr3*rr7
!!$           vald=vald+0.3333*(rr1d*rr3*rr7+rr1*rr3d*rr7+ &
!!$                    rr1*rr3*rr7d)
!!$         endif
!!$    else if(interval == ib.or.interval == np) then
!!$         val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
!!$         vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
!!$         vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
!!$         vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
!!$         if(iflag == 1) then
!!$           val=val+0.3333*(rr1*rr3*rr8 + rr1*rr4*rr9 + &
!!$                     rr2*rr5*rr9)
!!$           vald=vald+0.3333*(rr1d*rr3*rr8+rr1*rr3d*rr8+ &
!!$                    rr1*rr3*rr8d)
!!$           vald=vald+0.3333*(rr1d*rr4*rr9+rr1*rr4d*rr9+ &
!!$                    rr1*rr4*rr9d)
!!$           vald=vald+0.3333*(rr2d*rr5*rr9+rr2*rr5d*rr9+ &
!!$                    rr2*rr5*rr9d)
!!$         endif
!!$    else
!!$      val=0.
!!$    endif
!!$  else if(ib == np) then
!!$
!!$    r1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
!!$    r2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
!!$    r3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
!!$    r4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
!!$    r5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
!!$    r7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
!!$    r8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
!!$    r9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
!!$    r13=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
!!$
!!$    r1d=1./(xarr(np)-xarr(np-2))
!!$    r2d=-1./(xarr(np)-xarr(np-1))
!!$    r3d=1./(xarr(np)-xarr(np-2))
!!$    r4d=-1./(xarr(np)-xarr(np-1))
!!$    r5d=1./(xarr(np)-xarr(np-1))
!!$    r7d=1./(xarr(np-1)-xarr(np-2))
!!$    r8d=-1./  (xarr(np)-xarr(np-1))
!!$    r9d=1./(xarr(np)-xarr(np-1))
!!$    r13d=1./(xarr(np)-xarr(np-1))
!!$
!!$    if(interval == np-2) then
!!$         if(iflag == 0) then
!!$           val=r1*r3*r7
!!$           vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
!!$         else if(iflag == 1) then
!!$           val=0.6667*r1*r3*r7
!!$           vald=0.6667*(r1d*r3*r7+r1*r3d*r7+r1*r3*r7d)
!!$         endif
!!$    else if(interval == np-1.or.interval == np) then
!!$         if(iflag == 0) then
!!$           val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + r13**3
!!$           vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
!!$           vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
!!$           vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
!!$           vald=vald+3.*r13d*r13**2
!!$         else if(iflag == 1) then
!!$           val=0.6667*(r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + &
!!$                     1.5*r13**3)
!!$           vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
!!$           vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
!!$           vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
!!$           vald=vald+4.5*r13d*r13**2
!!$           vald=0.6667*vald
!!$         endif
!!$    else
!!$      val=0.
!!$    endif
!!$  else
!!$
!!$    r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
!!$    r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
!!$    r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
!!$    r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
!!$    r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
!!$    r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
!!$    r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
!!$    r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
!!$    r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
!!$   r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
!!$   r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
!!$   r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))
!!$
!!$    r1d=1./(xarr(ib+1)-xarr(ib-2))
!!$    r2d=-1./(xarr(ib+2)-xarr(ib-1))
!!$    r3d=1./(xarr(ib)-xarr(ib-2))
!!$    r4d=-1./(xarr(ib+1)-xarr(ib-1))
!!$    r5d=1./(xarr(ib+1)-xarr(ib-1))
!!$    r6d=-1./(xarr(ib+2)-xarr(ib))
!!$    r7d=1./(xarr(ib-1)-xarr(ib-2))
!!$    r8d=-1./  (xarr(ib)-xarr(ib-1))
!!$    r9d=1./(xarr(ib)-xarr(ib-1))
!!$   r10d=-1./(xarr(ib+1)-xarr(ib))
!!$   r11d=1./(xarr(ib+1)-xarr(ib))
!!$   r12d=-1./(xarr(ib+2)-xarr(ib+1))
!!$
!!$    if(interval == ib-2) then
!!$         val=r1*r3*r7
!!$         vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
!!$    else if(interval == ib-1) then
!!$         val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
!!$         vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
!!$         vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
!!$         vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
!!$    else if(interval == ib) then
!!$         val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
!!$         vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
!!$         vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
!!$         vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
!!$    else if(interval == ib+1) then
!!$         val=r2*r6*r12
!!$         vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
!!$    else
!!$      val=0.
!!$    endif
!!$  endif
!!$  splcon(ib)=val
!!$  splcond(ib)=vald
!!$  enddo
!!$
!!$  end subroutine vbspl

!----------------------------------
!----------------------------------

!!$  subroutine ylm(XLAT,XLON,LMAX,Y,WK1,WK2,WK3)
!!$
!!$  implicit none
!!$
!!$  complex TEMP,FAC,DFAC
!!$
!!$  real(kind=4) WK1(1),WK2(1),WK3(1),Y(1),XLAT,XLON
!!$
!!$  integer :: LMAX
!!$
!!$!
!!$!     WK1,WK2,WK3 SHOULD BE DIMENSIONED AT LEAST (LMAX+1)*4
!!$!
!!$  real(kind=4), parameter :: RADIAN = 57.2957795
!!$
!!$  integer :: IM,IL1,IND,LM1,L
!!$
!!$  real(kind=4) :: THETA,PHI
!!$
!!$  THETA=(90.-XLAT)/RADIAN
!!$  PHI=XLON/RADIAN
!!$
!!$  IND=0
!!$  LM1=LMAX+1
!!$
!!$  DO IL1=1,LM1
!!$
!!$  L=IL1-1
!!$  CALL legndr(THETA,L,L,WK1,WK2,WK3)
!!$
!!$  FAC=(1.,0.)
!!$  DFAC=CEXP(CMPLX(0.,PHI))
!!$
!!$  do IM=1,IL1
!!$    TEMP=FAC*CMPLX(WK1(IM),0.)
!!$    IND=IND+1
!!$    Y(IND)=REAL(TEMP)
!!$    IF(IM == 1) GOTO 20
!!$    IND=IND+1
!!$    Y(IND)=AIMAG(TEMP)
!!$ 20 FAC=FAC*DFAC
!!$  enddo
!!$
!!$  enddo
!!$
!!$  end subroutine ylm

!!$      subroutine legndr(THETA,L,M,X,XP,XCOSEC)
!!$      implicit none
!!$
!!$      integer :: L,M,i,k,LP1,MP1
!!$      real(kind=4) :: THETA,X,XP,XCOSEC,SFL3
!!$
!!$      DIMENSION X(2),XP(2),XCOSEC(2)
!!$      DOUBLE PRECISION SMALL,SUM,COMPAR,CT,ST,FCT,COT,FPI,X1,X2,X3,F1,F2,XM,TH,DSFL3,COSEC
!!$      DATA FPI/12.56637062D0/
!!$!      DFLOAT(I)=FLOAT(I)
!!$      SUM=0.D0
!!$      LP1=L+1
!!$      TH=THETA
!!$      CT=DCOS(TH)
!!$      ST=DSIN(TH)
!!$      MP1=M+1
!!$      FCT=DSQRT(dble(FLOAT(2*L+1))/FPI)
!!$      SFL3=SQRT(FLOAT(L*(L+1)))
!!$      COMPAR=dble(FLOAT(2*L+1))/FPI
!!$      DSFL3=SFL3
!!$      SMALL=1.D-16*COMPAR
!!$      do I=1,MP1
!!$       X(I)=0.
!!$       XCOSEC(I)=0.
!!$       XP(I)=0.
!!$      enddo
!!$      IF(L>1.AND.ABS(THETA)>1.E-5) GO TO 3
!!$      X(1)=FCT
!!$      IF(L==0) RETURN
!!$      X(1)=CT*FCT
!!$      X(2)=-ST*FCT/DSFL3
!!$      XP(1)=-ST*FCT
!!$      XP(2)=-.5D0*CT*FCT*DSFL3
!!$      IF(ABS(THETA)<1.E-5) XCOSEC(2)=XP(2)
!!$      IF(ABS(THETA)>=1.E-5) XCOSEC(2)=X(2)/ST
!!$      RETURN
!!$    3 X1=1.D0
!!$      X2=CT
!!$      DO  I=2,L
!!$       X3=(dble(FLOAT(2*I-1))*CT*X2-dble(FLOAT(I-1))*X1)/dble(FLOAT(I))
!!$       X1=X2
!!$       X2=X3
!!$      enddo
!!$      COT=CT/ST
!!$      COSEC=1./ST
!!$      X3=X2*FCT
!!$      X2=dble(FLOAT(L))*(X1-CT*X2)*FCT/ST
!!$      X(1)=X3
!!$      X(2)=X2
!!$      SUM=X3*X3
!!$      XP(1)=-X2
!!$      XP(2)=dble(FLOAT(L*(L+1)))*X3-COT*X2
!!$      X(2)=-X(2)/SFL3
!!$      XCOSEC(2)=X(2)*COSEC
!!$      XP(2)=-XP(2)/SFL3
!!$      SUM=SUM+2.D0*X(2)*X(2)
!!$      IF(SUM-COMPAR>SMALL) RETURN
!!$      X1=X3
!!$      X2=-X2/DSQRT(dble(FLOAT(L*(L+1))))
!!$      DO  I=3,MP1
!!$       K=I-1
!!$       F1=DSQRT(dble(FLOAT((L+I-1)*(L-I+2))))
!!$       F2=DSQRT(dble(FLOAT((L+I-2)*(L-I+3))))
!!$       XM=K
!!$       X3=-(2.D0*COT*(XM-1.D0)*X2+F2*X1)/F1
!!$       SUM=SUM+2.D0*X3*X3
!!$       IF(SUM-COMPAR>SMALL.AND.I/=LP1) RETURN
!!$       X(I)=X3
!!$       XCOSEC(I)=X(I)*COSEC
!!$       X1=X2
!!$       XP(I)=-(F1*X2+XM*COT*X3)
!!$       X2=X3
!!$      enddo
!!$      RETURN
!!$      end subroutine legndr

