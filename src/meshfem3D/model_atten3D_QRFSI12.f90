!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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
!   3-D attenuation model QRFSI12 (Dalton, Ekstr"om, & Dziewonski, 2008)
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
  integer,parameter :: NKQ = 8,MAXL_Q = 12
  integer,parameter :: NSQ=(MAXL_Q+1)**2,NDEPTHS_REFQ = 913

  ! model_atten3D_QRFSI12_variables
  double precision,dimension(:,:),allocatable :: QRFSI12_Q_dqmu
  double precision,dimension(:),allocatable :: QRFSI12_Q_spknt
  double precision,dimension(:),allocatable :: QRFSI12_Q_refdepth,QRFSI12_Q_refqmu

  ! helper array
  real(kind=4),dimension(:),allocatable :: xlmvec

  end module model_atten3D_QRFSI12_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_atten3D_QRFSI12_broadcast()

! standard routine to setup model

  use constants
  use model_atten3D_QRFSI12_par

  implicit none

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(QRFSI12_Q_dqmu(NKQ,NSQ), &
           QRFSI12_Q_spknt(NKQ), &
           QRFSI12_Q_refdepth(NDEPTHS_REFQ), &
           QRFSI12_Q_refqmu(NDEPTHS_REFQ), &
           stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating QRFSI12_Q model arrays')

  ! main process reads in file values
  if (myrank == 0) call read_atten_model_3D_QRFSI12()

  ! broadcasts to all processes
  call bcast_all_dp(QRFSI12_Q_dqmu,          NKQ*NSQ)
  call bcast_all_dp(QRFSI12_Q_spknt,             NKQ)
  call bcast_all_dp(QRFSI12_Q_refdepth, NDEPTHS_REFQ)
  call bcast_all_dp(QRFSI12_Q_refqmu,   NDEPTHS_REFQ)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'read 3D attenuation model'
    write(IMAIN,*) '  model: QRFSI12'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! helper array
  allocate(xlmvec(NSQ**2),stat=ier)
  if (ier /= 0) stop 'Error allocating helper array'
  xlmvec(:) = 0.0

  end subroutine

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_atten_model_3D_QRFSI12()

  use constants
  use model_atten3D_QRFSI12_par

  implicit none

  ! local parameters
  integer :: j,k,l,m,ier
  integer :: indexval,ll,mm
  double precision :: v1,v2

  character(len=MAX_STRING_LEN) :: QRFSI12,QRFSI12_ref

! read in QRFSI12
! hard-wire for now
  QRFSI12='DATA/QRFSI12/QRFSI12.dat'
  QRFSI12_ref='DATA/QRFSI12/ref_QRFSI12'

! get the dq model coefficients
  open(unit=IIN,file=QRFSI12,status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(QRFSI12), '": ', ier
    call exit_MPI(0, 'Error model QRFSI12')
  endif

  do k = 1,NKQ
    read(IIN,*) indexval
    j = 0
    do l = 0,MAXL_Q
      do m = 0,l
        if (m == 0) then
          j = j+1
          read(IIN,*)ll,mm,v1
          QRFSI12_Q_dqmu(k,j)=v1
        else
          j = j+2
          read(IIN,*)ll,mm,v1,v2
  !       write(*,*) 'k,l,m,ll,mm:',k,l,m,ll,mm,v1
          QRFSI12_Q_dqmu(k,j-1)=2.*v1
          QRFSI12_Q_dqmu(k,j)=-2.*v2
        endif
      enddo
    enddo
  enddo
  close(IIN)

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
  open(IIN,file=QRFSI12_ref,status='old',action='read',iostat=ier)
  if (ier /= 0) then
    write(IMAIN,*) 'Error opening "', trim(QRFSI12_ref), '": ', ier
    call exit_MPI(0, 'Error model QRFSI12')
  endif

  do j = 1,NDEPTHS_REFQ
    read(IIN,*)QRFSI12_Q_refdepth(j),QRFSI12_Q_refqmu(j)
  enddo
  close(IIN)


  end subroutine read_atten_model_3D_QRFSI12

!
!----------------------------------
!

  subroutine model_atten3D_QRFSI12(radius,theta,phi,Qmu,idoubling)

  use constants
  use shared_parameters, only: R_PLANET_KM

  use model_atten3D_QRFSI12_par

  use model_prem_par, only: PREM_RMOHO

  implicit none

  double precision :: radius,theta,phi,Qmu
  integer :: idoubling

  ! local parameters
  integer :: i,j,k,n
  integer :: ifnd
  double precision :: smallq,dqmu,smallq_ref
  real(kind=4) :: splpts(NKQ),splcon(NKQ),splcond(NKQ)
  real(kind=4) :: depth,ylat,xlon
  real(kind=4) :: shdep(NSQ)


  ! in Colleen's original code theta refers to the latitude.  Here we have redefined theta to be colatitude
  ! to agree with the rest of specfem

  ! debug
  !  print *,'entering QRFSI12 subroutine'

  ylat = real(90.0d0 - theta,kind=4)
  xlon = real(phi,kind=4)

  ! only checks radius for crust, idoubling is misleading for oceanic crust
  ! when we want to expand mantle up to surface...

  ! radius given in km

!  !if (idoubling == IFLAG_CRUST .or. radius >= rmoho) then
  if (radius >= PREM_RMOHO/1000.d0) then
  !   print *,'QRFSI12: we are in the crust'
     Qmu = 600.0d0
  else if (idoubling == IFLAG_INNER_CORE_NORMAL .or. idoubling == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
       idoubling == IFLAG_BOTTOM_CENTRAL_CUBE .or. idoubling == IFLAG_TOP_CENTRAL_CUBE .or. &
       idoubling == IFLAG_IN_FICTITIOUS_CUBE) then
    ! we are in the inner core

    !debug
    !   print *,'QRFSI12: we are in the inner core'

    Qmu = 84.6d0

  else if (idoubling == IFLAG_OUTER_CORE_NORMAL) then

    ! we are in the outer core

    !debug
    !   print *,'QRFSI12: we are in the outer core'

    Qmu = 0.0d0

  else

    ! we are in the mantle

    depth = real(R_PLANET_KM - radius,kind=4)

    !debug
    !   print *,'QRFSI12: we are in the mantle at depth',depth

    ifnd = 0
    do i = 2,NDEPTHS_REFQ
      if (depth >= QRFSI12_Q_refdepth(i-1) .and. depth < QRFSI12_Q_refdepth(i)) then
        ifnd = i
      endif
    enddo
    if (ifnd == 0) then
      write(*,"('problem finding reference Q value at depth: ',f8.3)") depth
      stop
    endif
    smallq_ref = QRFSI12_Q_refqmu(ifnd)
    smallq = smallq_ref

! Colleen's model is only defined between depths of 24.4 and 650km
    if (depth < 650.d0) then
      do j = 1,NSQ
        shdep(j) = 0.0
      enddo
      do n = 1,NKQ
        splpts(n) = real(QRFSI12_Q_spknt(n),kind=4)
      enddo

      call vbspl(depth,NKQ,splpts,splcon,splcond)

      do n = 1,NKQ
        do j = 1,NSQ
          shdep(j) = shdep(j) + real(splcon(n) * QRFSI12_Q_dqmu(n,j),kind=4)
        enddo
      enddo

      call ylm(ylat,xlon,MAXL_Q,xlmvec)

      dqmu = 0.0
      do k = 1,NSQ
        dqmu = dqmu+xlmvec(k)*shdep(k)
      enddo
      smallq = smallq_ref + dqmu
    endif

    ! if smallq is small and negative (due to numerical error), Qmu is very large:
    if (smallq < 0.0d0) smallq = 1.0d0/ATTENUATION_COMP_MAXIMUM

    Qmu = 1/smallq

    ! if Qmu is larger than MAX_ATTENUATION_VALUE, set it to ATTENUATION_COMP_MAXIMUM.
    ! This assumes that this value is high enough that at this point there is almost no attenuation at all.
    if (Qmu >= ATTENUATION_COMP_MAXIMUM) Qmu = 0.99d0*ATTENUATION_COMP_MAXIMUM

  endif

  end subroutine model_atten3D_QRFSI12
