!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
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

  program convert_palette

! convert GMT palette for the globe found at
! http://sview01.wiredworkplace.net/pub/cpt-city/gmt/tn/GMT_globe.png.index.html
! to proprietary OpenDX *.cm format

  implicit none

  integer, parameter :: NUM_LINES = 38

! opacity: 1=totally opaque  2=continents 80% opaque and oceans 30% opaque
  integer, parameter :: OPACITY = 2

  integer, dimension(NUM_LINES) :: r,g,b
  real, dimension(NUM_LINES) :: z

  integer :: i,icomp

  real :: rval,gval,bval,hval,sval,vval

! read the input GMT palette "GMT_globe_palette_cropped.dat" or "GMT_globe_palette.dat"
  do i = 1,NUM_LINES
    read(*,*) z(i),r(i),g(i),b(i)
  enddo

! for OpenDX *.cm format, output min and max of z
  print *,minval(z)
  print *,maxval(z)

! output the palette converted from RGB to HSV

! output H, then S, then V
  do icomp = 1,3

! for OpenDX *.cm format, output number of values
  print *,NUM_LINES
  do i = 1,NUM_LINES

    rval = r(i) / 255.d0
    gval = g(i) / 255.d0
    bval = b(i) / 255.d0

      call rgb_to_hsv ( rval, gval, bval, hval, sval, vval )

! for OpenDX *.cm format, output value of H,S or V and then three zeroes
! hue is an angle between 0 and 360 degrees but it is stored as a real value between 0 and 1
      if(icomp == 1) print *,(z(i)-z(1))/(z(NUM_LINES)-z(1)),hval/360.,' 0 0 0'
      if(icomp == 2) print *,(z(i)-z(1))/(z(NUM_LINES)-z(1)),sval,' 0 0 0'
      if(icomp == 3) print *,(z(i)-z(1))/(z(NUM_LINES)-z(1)),vval,' 0 0 0'

  enddo

  enddo

! output opacity, set to 1 (i.e. totally opaque)
  if(OPACITY == 1) then
    print *,'2'
    print *,'0.00000000000000000000   1.00000000000000000000  0  0  0'
    print *,'1.00000000000000000000   1.00000000000000000000  0  0  0'
  else if(OPACITY == 2) then
    print *,'4'
    print *,'0.0000   0.00  0  0  0'
    print *,'0.5934   0.00  0  0  0'
    print *,'0.5934   0.44  0  0  0'
    print *,'1.0000   0.44  0  0  0'
  else
    stop 'wrong opacity code'
  endif

  end program convert_palette

!---- Fortran conversion from RGB to HSV below
!---- taken from http://orion.math.iastate.edu/burkardt/f_src/colors/colors.f90

subroutine rgb_to_hsv ( r, g, b, h, s, v )
!
!*******************************************************************************
!
!! RGB_TO_HSV converts RGB to HSV color coordinates.
!
!
!  Definition:
!
!    The RGB color system describes a color based on the amounts of the
!    base colors red, green, and blue.  Thus, a particular color
!    has three coordinates, (R,G,B).  Each coordinate must be between
!    0 and 1.
!
!    The HSV color system describes a color based on the three qualities
!    of hue, saturation, and value.  A given color will be represented
!    by three numbers, (H,S,V).  H, the value of hue, is an angle
!    between 0 and 360 degrees, with 0 representing red.  S is the
!    saturation, and is between 0 and 1.  Finally, V is the "value",
!    a measure of brightness, which goes from 0 for black, increasing
!    to a maximum of 1 for the brightest colors.  The HSV color system
!    is sometimes also called HSB, where the B stands for brightness.
!
!  Reference:
!
!    Foley, van Dam, Feiner, and Hughes,
!    Computer Graphics, Principles and Practice,
!    Addison Wesley, Second Edition, 1990.
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real R, G, B, the RGB color coordinates to be converted.
!
!    Output, real H, S, V, the corresponding HSV color coordinates.

  implicit none

  real b
  real bc
  real g
  real gc
  real h
  real r
  real rc
  real rgbmax
  real rgbmin
  real r_modp
  real s
  real v

  rgbmax = max ( r, g, b )
  rgbmin = min ( r, g, b )

  v = rgbmax
!
!  Compute the saturation.
!
  if ( rgbmax /= 0.0E+00 ) then
    s = ( rgbmax - rgbmin ) / rgbmax
  else
    s = 0.0E+00
  end if
!
!  Compute the hue.
!
  if ( s == 0.0E+00 ) then

    h = 0.0E+00

  else

    rc = ( rgbmax - r ) / ( rgbmax - rgbmin )
    gc = ( rgbmax - g ) / ( rgbmax - rgbmin )
    bc = ( rgbmax - b ) / ( rgbmax - rgbmin )

    if ( r == rgbmax ) then
      h = bc - gc
    else if ( g == rgbmax ) then
      h = 2.0E+00 + rc - bc
    else
      h = 4.0E+00 + gc - rc
    end if

    h = h * 60.0E+00
!
!  Make sure H lies between 0 and 360.0E+00
!
    h = r_modp ( h, 360.0E+00 )

  end if

end

!-------------------

function r_modp ( x, y )
!
!*******************************************************************************
!
!! R_MODP returns the nonnegative remainder of real division.
!
!
!  Formula:
!
!    If
!      REM = R_MODP ( X, Y )
!      RMULT = ( X - REM ) / Y
!    then
!      X = Y * RMULT + REM
!    where REM is always nonnegative.
!
!  Comments:
!
!    The MOD function computes a result with the same sign as the
!    quantity being divided.  Thus, suppose you had an angle A,
!    and you wanted to ensure that it was between 0 and 360.
!    Then mod(A,360.0) would do, if A was positive, but if A
!    was negative, your result would be between -360 and 0.
!
!    On the other hand, R_MODP(A,360.0) is between 0 and 360, always.
!
!  Examples:
!
!        I         J     MOD   R_MODP   R_MODP Factorization
!
!      107        50       7       7    107 =  2 *  50 + 7
!      107       -50       7       7    107 = -2 * -50 + 7
!     -107        50      -7      43   -107 = -3 *  50 + 43
!     -107       -50      -7      43   -107 =  3 * -50 + 43
!
!  Modified:
!
!    29 August 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real X, the number to be divided.
!
!    Input, real Y, the number that divides X.
!
!    Output, real R_MODP, the nonnegative remainder when X is divided by Y.
!
  implicit none

  real r_modp
  real x
  real y

  if ( y == 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R_MODP - Fatal error!'
    write ( *, '(a,g14.6)' ) '  R_MODP ( X, Y ) called with Y = ', y
    stop
  end if

  r_modp = mod ( x, y )

  if ( r_modp < 0.0E+00 ) then
    r_modp = r_modp + abs ( y )
  end if

end

