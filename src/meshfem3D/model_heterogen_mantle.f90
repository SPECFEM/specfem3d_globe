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
! HMM
!
! generic heterogeneous mantle model
!--------------------------------------------------------------------------------------------------

  module model_heterogen_mantle_par

  ! heterogen_mantle_model_constants
  integer, parameter :: N_R = 256,N_THETA = 256,N_PHI = 256

  ! model array
  double precision,dimension(:),allocatable :: HMM_rho_in

  end module model_heterogen_mantle_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_heterogen_mntl_broadcast(myrank)

! standard routine to setup model

  use constants
  use model_heterogen_mantle_par

  implicit none

  integer :: myrank

  ! local parameters
  integer :: ier

  ! allocates model array
  allocate(HMM_rho_in(N_R*N_THETA*N_PHI), &
          stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating HMM array')

  ! master process reads in model
  if (myrank == 0) then
     write(IMAIN,*) 'Reading in model_heterogen_mantle.'
     call flush_IMAIN()

     call read_heterogen_mantle_model()

     write(IMAIN,*) 'model_heterogen_mantle is read in.'
     call flush_IMAIN()
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_dp(HMM_rho_in,N_R*N_THETA*N_PHI)

  if (myrank == 0) then
    write(IMAIN,*) 'model_heterogen_mantle is broadcast.'
    write(IMAIN,*) 'First value in HMM:',HMM_rho_in(1)
    write(IMAIN,*) 'Last value in HMM:',HMM_rho_in(N_R*N_THETA*N_PHI)
    call flush_IMAIN()
  endif

  end subroutine model_heterogen_mntl_broadcast

!
!-------------------------------------------------------------------------------------------------
!


!
! NOTE: CURRENTLY THIS ROUTINE ONLY WORKS FOR N_R=N_THETA=N_PHI !!!!!
!

  subroutine read_heterogen_mantle_model()

  use constants
  use model_heterogen_mantle_par

  implicit none

  ! local parameters
  integer :: i,j,ier

  ! open heterogen.dat
  open(unit=IIN,file='./DATA/heterogen/heterogen.dat',access='direct',&
       form='formatted',recl=20,status='old',action='read',iostat=ier)
  if (ier /= 0 ) call exit_MPI(0,'Error opening model file heterogen.dat')

  j = N_R*N_THETA*N_PHI

  do i = 1,j
    read(IIN,rec=i,fmt='(F20.15)') HMM_rho_in(i)
  enddo

  close(IIN)

  end subroutine read_heterogen_mantle_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_heterogen_mantle(radius,theta,phi,dvs,dvp,drho)

  use constants
  use model_heterogen_mantle_par

  implicit none

  ! variable declaration
  double precision :: radius,theta,phi            ! input coordinates
  double precision :: drho,dvp,dvs                ! output anomaly values

  ! local parameters
  double precision :: x,y,z                       ! input converted to Cartesian
  double precision :: x_low,x_high                ! x values used to interpolate
  double precision :: y_low,y_high                ! y values used to interpolate
  double precision :: z_low,z_high                ! z values used to interpolate
  double precision :: delta,delta2                ! weights in record# and in interpolation
  double precision :: rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8 ! rho values at the interpolation points
  double precision :: r_inner,r_outer             ! lower and upper domain bounds for r
  integer :: rec_read                             ! nr of record to be read from heterogen.dat (direct access file)
  double precision :: a,b,c                       ! substitutions in interpolation algorithm (weights)

  radius = radius*R_EARTH
  r_inner = 3.500d6  !lower bound for heterogeneity zone
! NOTE: r_outer NEEDS TO BE (just) SMALLER THAN R_EARTH!!!!!!!!
  r_outer = R_EARTH-1.0d1  !6.300d6  !upper bound for heterogeneity zone (lower mantle: e.g. 4.500d6)

  delta = 2.*R_EARTH/(real(N_R-1))
  delta2 = 2.*R_EARTH/(real(N_R-2))
  !delta2 = 2.*R_EARTH/(real(N_R))

  if ((radius >= r_inner) .and. (radius <= r_outer)) then
    ! convert spherical point to Cartesian point, move origin to corner
    x = R_EARTH + radius*sin(theta)*cos(phi)
    y = R_EARTH + radius*sin(theta)*sin(phi)
    z = R_EARTH + radius*cos(theta)

    ! determine which points to search for in heterogen.dat
    ! find x_low,y_low,z_low etc.
    x_low = floor(x/delta2) + 1
    x_high = x_low + 1
    y_low = floor(y/delta2) + 1
    y_high = y_low + 1
    z_low = floor(z/delta2) + 1
    z_high = z_low + 1

    ! rho1 at: x_low y_low z_low
    rec_read = int(1+(x_low*N_R*N_R)+(y_low*N_R)+z_low)
    rho1 = HMM_rho_in(rec_read)

    ! rho2 at: x_low y_high z_low
    rec_read = int(1+(x_low*N_R*N_R)+(y_high*N_R)+z_low)
    rho2 = HMM_rho_in(rec_read)

    ! rho3 at: x_high y_low z_low
    rec_read = int(1+(x_high*N_R*N_R)+(y_low*N_R)+z_low)
    rho3 = HMM_rho_in(rec_read)

    ! rho4 at: x_high y_high z_low
    rec_read = int(1+(x_high*N_R*N_R)+(y_high*N_R)+z_low)
    rho4 = HMM_rho_in(rec_read)

    ! rho5 at: x_low y_low z_high
    rec_read = int(1+(x_low*N_R*N_R)+(y_low*N_R)+z_high)
    rho5 = HMM_rho_in(rec_read)

    ! rho6 at: x_low y_high z_high
    rec_read = int(1+(x_low*N_R*N_R)+(y_high*N_R)+z_high)
    rho6 = HMM_rho_in(rec_read)

    ! rho7 at: x_high y_low z_high
    rec_read = int(1+(x_high*N_R*N_R)+(y_low*N_R)+z_high)
    rho7 = HMM_rho_in(rec_read)

    ! rho8 at: x_high y_high z_high
    rec_read = int(1+(x_high*N_R*N_R)+(y_high*N_R)+z_high)
    rho8 = HMM_rho_in(rec_read)

    ! perform linear interpolation between the 8 points
    a = (x-x_low*delta)/delta       ! weight for x
    b = (y-y_low*delta)/delta       ! weight for y
    c = (z-z_low*delta)/delta       ! weight for z

    drho = rho1*(1.-a)*(1.-b)*(1.-c) + rho2*(1.-a)*b*(1.-c) + &
           rho3*a*(1.-b)*(1.-c) + rho4*a*b*(1.-c) + rho5*(1.-a)*(1.-b)*c + &
           rho6*(1.-a)*b*c + rho7*a*(1.-b)*c + rho8*a*b*c

    ! calculate delta vp,vs from the interpolated delta rho
    dvp = (0.55/0.30)*drho
    dvs = (1.00/0.30)*drho

  else !outside of heterogeneity domain
    drho = 0.
    dvp = 0.
    dvs = 0.
  endif

  end subroutine model_heterogen_mantle
