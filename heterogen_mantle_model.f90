!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                 and University of Pau / CNRS, France
! (c) California Institute of Technology and University of Pau / CNRS, November 2007
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
! NOTE: CURRENTLY THIS ROUTINE ONLY WORKS FOR N_R=N_THETA=N_PHI !!!!!
!

  subroutine read_heterogen_mantle_model(HMM)

  implicit none

  include "constants.h"

  integer i,j

! heterogen_mod_variables
  type heterogen_mod_variables
    sequence
    double precision rho_in(N_R*N_THETA*N_PHI)
  end type heterogen_mod_variables

  type (heterogen_mod_variables) HMM
! heterogen_mod_variables


! open heterogen.dat
  open(unit=10,file='./DATA/heterogen/heterogen.dat',access='direct',&
       form='formatted',recl=20,status='old',action='read')

  j = N_R*N_THETA*N_PHI

  do i = 1,j
    read(10,rec=i,fmt='(F20.15)') HMM%rho_in(i)
  end do

  close(10)

  end subroutine read_heterogen_mantle_model

!====================================================================

  subroutine heterogen_mantle_model(radius,theta,phi,dvs,dvp,drho,HMM)

  implicit none

  include "constants.h"

  ! variable declaration
  double precision radius,theta,phi            ! input coordinates
  double precision x,y,z                       ! input converted to cartesian
  double precision drho,dvp,dvs                ! output anomaly values
  double precision x_low,x_high                ! x values used to interpolate
  double precision y_low,y_high                ! y values used to interpolate
  double precision z_low,z_high                ! z values used to interpolate
  double precision delta,delta2                ! weigts in record# and in interpolation
  double precision rho1,rho2,rho3,rho4,rho5,rho6,rho7,rho8 ! rho values at the interpolation points
  double precision r_inner,r_outer             ! lower and upper domain bounds for r
  integer rec_read                             ! nr of record to be read from heterogen.dat (direct access file)
  double precision a,b,c                       ! substitutions in interpolation algorithm (weights)


! heterogen_mod_variables
  type heterogen_mod_variables
    sequence
    double precision rho_in(N_R*N_THETA*N_PHI)
  end type heterogen_mod_variables

  type (heterogen_mod_variables) HMM
! heterogen_mod_variables

  radius = radius*R_EARTH
  r_inner = 3.500d6  !lower bound for heterogeneity zone
! NOTE: r_outer NEEDS TO BE (just) SMALLER THAN R_EARTH!!!!!!!!
  r_outer = R_EARTH-1.0d1  !6.300d6  !upper bound for heterogeneity zone (lower mantle: e.g. 4.500d6)

  delta = 2.*R_EARTH/(real(N_R-1))
  delta2 = 2.*R_EARTH/(real(N_R-2))
  !delta2 = 2.*R_EARTH/(real(N_R))

  if ((radius >= r_inner) .and. (radius <= r_outer)) then
    ! convert spherical point to cartesian point, move origin to corner
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
    rec_read = 1+(x_low*N_R*N_R)+(y_low*N_R)+z_low
    rho1 = HMM%rho_in(rec_read)

    ! rho2 at: x_low y_high z_low
    rec_read = 1+(x_low*N_R*N_R)+(y_high*N_R)+z_low
    rho2 = HMM%rho_in(rec_read)

    ! rho3 at: x_high y_low z_low
    rec_read = 1+(x_high*N_R*N_R)+(y_low*N_R)+z_low
    rho3 = HMM%rho_in(rec_read)

    ! rho4 at: x_high y_high z_low
    rec_read = 1+(x_high*N_R*N_R)+(y_high*N_R)+z_low
    rho4 = HMM%rho_in(rec_read)

    ! rho5 at: x_low y_low z_high
    rec_read = 1+(x_low*N_R*N_R)+(y_low*N_R)+z_high
    rho5 = HMM%rho_in(rec_read)

    ! rho6 at: x_low y_high z_high
    rec_read = 1+(x_low*N_R*N_R)+(y_high*N_R)+z_high
    rho6 = HMM%rho_in(rec_read)

    ! rho7 at: x_high y_low z_high
    rec_read = 1+(x_high*N_R*N_R)+(y_low*N_R)+z_high
    rho7 = HMM%rho_in(rec_read)

    ! rho8 at: x_high y_high z_high
    rec_read = 1+(x_high*N_R*N_R)+(y_high*N_R)+z_high
    rho8 = HMM%rho_in(rec_read)

    ! perform linear interpolation between the 8 points
    a = (x-x_low*delta)/delta       ! weight for x
    b = (y-y_low*delta)/delta       ! weight for y
    c = (z-z_low*delta)/delta       ! weight for z

    drho = rho1*(1.-a)*(1.-b)*(1.-c) + rho2*(1.-a)*b*(1.-c) + &
     & rho3*a*(1.-b)*(1.-c) + rho4*a*b*(1.-c) + rho5*(1.-a)*(1.-b)*c + &
     & rho6*(1.-a)*b*c + rho7*a*(1.-b)*c + rho8*a*b*c

    ! calculate delta vp,vs from the interpolated delta rho
    dvp = (0.55/0.30)*drho
    dvs = (1.00/0.30)*drho

  else !outside of heterogeneity domain
    drho = 0.
    dvp = 0.
    dvs = 0.
  end if

  end subroutine heterogen_mantle_model
