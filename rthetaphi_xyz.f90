!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine xyz_2_rthetaphi(x,y,z,r,theta,phi)

! convert x y z to r theta phi, single precision call

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) x,y,z,r,theta,phi
  double precision xmesh,ymesh,zmesh

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

    xmesh = dble(x)
    ymesh = dble(y)
    zmesh = dble(z)

    if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
    if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE
    theta = sngl(datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh))
    if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
    if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE
    phi = sngl(datan2(ymesh,xmesh))

    r = sngl(dsqrt(xmesh**2 + ymesh**2 + zmesh**2))

  else

    xmesh = x
    ymesh = y
    zmesh = z

    if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
    if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE
    theta = datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh)
    if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
    if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE
    phi = datan2(ymesh,xmesh)

    r = dsqrt(xmesh**2 + ymesh**2 + zmesh**2)

  endif

  end subroutine xyz_2_rthetaphi

!-------------------------------------------------------------

  subroutine xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

! convert x y z to r theta phi, double precision call

  implicit none

  include "constants.h"

  double precision x,y,z,r,theta,phi
  double precision xmesh,ymesh,zmesh

  xmesh = x
  ymesh = y
  zmesh = z

  if(zmesh > -SMALL_VAL_ANGLE .and. zmesh <= ZERO) zmesh = -SMALL_VAL_ANGLE
  if(zmesh < SMALL_VAL_ANGLE .and. zmesh >= ZERO) zmesh = SMALL_VAL_ANGLE
  theta = datan2(dsqrt(xmesh*xmesh+ymesh*ymesh),zmesh)
  if(xmesh > -SMALL_VAL_ANGLE .and. xmesh <= ZERO) xmesh = -SMALL_VAL_ANGLE
  if(xmesh < SMALL_VAL_ANGLE .and. xmesh >= ZERO) xmesh = SMALL_VAL_ANGLE
  phi = datan2(ymesh,xmesh)

  r = dsqrt(xmesh**2 + ymesh**2 + zmesh**2)

  end subroutine xyz_2_rthetaphi_dble

!-------------------------------------------------------------

  subroutine rthetaphi_2_xyz(x,y,z,r,theta,phi)

! convert r theta phi to x y z

  implicit none

  include "constants.h"

  real(kind=CUSTOM_REAL) x,y,z,r,theta,phi

  x = r * sin(theta) * cos(phi)
  y = r * sin(theta) * sin(phi)
  z = r * cos(theta)

  end subroutine rthetaphi_2_xyz

