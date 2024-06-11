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

! get backazimuth baz from event and station coordinates the, phe, ths and phs
  subroutine get_backazimuth(the,phe,ths,phs,Baz)

  use constants, only: DEGREES_TO_RADIANS,RADIANS_TO_DEGREES

  implicit none

  ! event location
  double precision,intent(in) :: the, phe
  ! station location
  double precision,intent(in) :: ths, phs
  ! back-azimuth
  double precision,intent(inout) :: Baz

  double precision :: az,xdeg
  double precision :: a, a1, b, b1, c, c1
  double precision :: d, d1, e, e1
  double precision :: ec2, f, f1, g, g1, h, h1, onemec2, pherad
  double precision :: phsrad, sc, sd, ss
  double precision :: temp, therad, thg, thsrad

  !double precision, parameter :: rad = 6378.160d0
  double precision, parameter :: fl = 0.00335293d0
  double precision, parameter :: twopideg = 360.d0
  !double precision, parameter :: c00 = 1.d0
  !double precision, parameter :: c01 = 0.25d0
  !double precision, parameter :: c02 = -4.6875d-02
  !double precision, parameter :: c03 = 1.953125d-02
  !double precision, parameter :: c21 = -0.125d0
  !double precision, parameter :: c22 = 3.125d-02
  !double precision, parameter :: c23 = -1.46484375d-02
  !double precision, parameter :: c42 = -3.90625d-03
  !double precision, parameter :: c43 = 2.9296875d-03
  !double precision, parameter :: degtokm = 111.3199d0
  double precision, parameter :: TORAD = DEGREES_TO_RADIANS
  double precision, parameter :: TODEG = RADIANS_TO_DEGREES

  !=====================================================================
  ! PURPOSE:  To compute the distance and azimuth between locations.
  !=====================================================================
  ! INPUT ARGUMENTS:
  !    THE:     Event latitude in decimal degrees, North positive. [r]
  !    PHE:     Event longitude, East positive. [r]
  !    THS:     Array of station latitudes. [r]
  !    PHS:     Array of station longitudes. [r]
  !    NS:      Length of THS and PHS. [i]
  !=====================================================================
  ! OUTPUT ARGUMENTS:
  !    DIST:    Array of epicentral distances in km. [r]
  !    AZ:      Array of azimuths in degrees. [r]
  !    BAZ:     Array of back azimuths. [r]
  !    XDEG:    Array of great circle arc lengths. [r]
  !    NERR:    Error flag:
  !             =    0   No error.
  !             = 0904   Calculation failed internal consistency checks.
  !=====================================================================
  ! MODULE/LEVEL:  DFM/4
  !=====================================================================
  ! GLOBAL INPUT:
  !    MACH:
  !=====================================================================
  ! subroutines called:
  !    SACLIB:  SETMSG, APCMSG
  !=====================================================================
  ! LOCAL VARIABLES:
  !=====================================================================
  ! KNOWN ERRORS:
  ! - Problem with equation for distance. See discussion below.
  !=====================================================================
  ! PROCEDURE:
  ! - Calculations are based upon the reference spheroid of 1968 and
  !   are defined by the major radius (RAD) and the flattening (FL).
  ! - Initialize.
  !nerr = 0

  ec2 = 2.d0*fl - fl*fl
  onemec2 = 1.d0 - ec2

  ! - Convert event location to radians.
  !   (Equations are unstable for latitudes of exactly 0 degrees.)

  temp = the
  if (temp == 0.d0) temp = 1.0d-08
  therad = TORAD*temp
  pherad = TORAD*phe

  ! - Must convert from geographic to geocentric coordinates in order
  !   to use the spherical trig equations.  This requires a latitude
  !   correction given by: 1-EC2=1-2*FL+FL*FL

  if (the == 90.d0 .or. the == -90.d0) then         ! special attention at the poles
    thg = the*TORAD                   ! ... to avoid division by zero.
  else
    thg = atan( onemec2*tan( therad ) )
  endif

  d = sin( pherad )
  e = -cos( pherad )
  f = -cos( thg )
  c = sin( thg )
  a = f*e
  b = -f*d
  g = -c*e
  h = c*d

  ! -- Convert to radians.
  temp = Ths
  if (temp == 0.d0 ) temp = 1.0d-08
  thsrad = TORAD*temp
  phsrad = TORAD*Phs

  ! -- Calculate some trig constants.
  if (Ths == 90.d0 .or. Ths == -90.d0) then
    thg = Ths * TORAD
  else
    thg = atan( onemec2*tan( thsrad ) )
  endif

  d1 = sin( phsrad )
  e1 = -cos( phsrad )
  f1 = -cos( thg )
  c1 = sin( thg )
  a1 = f1*e1
  b1 = -f1*d1
  g1 = -c1*e1
  h1 = c1*d1
  sc = a*a1 + b*b1 + c*c1

  ! - Spherical trig relationships used to compute angles.

  sd = 0.5d0 * sqrt( ((a - a1)**2 + (b - b1)**2 + (c - c1)**2)*((a + a1)**2 + (b + b1)**2 + (c + c1)**2) )
  Xdeg = atan2( sd, sc )*TODEG
  if (Xdeg < 0.d0 ) Xdeg = Xdeg + twopideg

  ss = (a1 - d)**2 + (b1 - e)**2 + (c1)**2 - 2.
  sc = (a1 - g)**2 + (b1 - h)**2 + (c1 - f)**2 - 2.
  Az = atan2( ss, sc )*TODEG
  if (Az < 0.d0 ) Az = Az + twopideg

  ss = (a - d1)**2 + (b - e1)**2 + (c)**2 - 2.
  sc = (a - g1)**2 + (b - h1)**2 + (c - f1)**2 - 2.
  Baz = atan2( ss, sc )*TODEG
  if (Baz < 0.d0) Baz = Baz + twopideg

  end subroutine get_backazimuth
