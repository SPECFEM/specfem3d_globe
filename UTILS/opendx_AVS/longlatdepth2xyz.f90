!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!----
!----  converts lat long depth to x y z
!----

  program convert

  implicit none

! R_EARTH is the radius of the 1-D Earth
  double precision, parameter :: R_EARTH = 6371000.d0

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0

  double precision theta,phi,latitude,longitude,depth,r0
  double precision x_target,y_target,z_target,r_target

! **************

! latitude and longitude in degrees
! depth of hypocenter in kilometers
  latitude = -16.0800
  longitude = 168.3100
  depth = 15.0000

! convert geographic latitude latitude (degrees)
! to geocentric colatitude theta (radians)
  theta=PI/2.0d0-atan(0.99329534d0*dtan(latitude*PI/180.0d0))
  phi=longitude*PI/180.0d0

! normalized Earth radius
  r0 = 1.d0

! compute the Cartesian position of the source
  r_target = r0 - depth*1000.0d0/R_EARTH
  x_target = r_target*dsin(theta)*dcos(phi)
  y_target = r_target*dsin(theta)*dsin(phi)
  z_target = r_target*dcos(theta)

! print result
  print *
  print *,'long = ',longitude
  print *,'lat = ',latitude
  print *,'depth (km) = ',depth
  print *
  print *,'x = ',x_target
  print *,'y = ',y_target
  print *,'z = ',z_target
  print *
  print *,'radius = ',sqrt(x_target**2 + y_target**2 + z_target**2)
  print *

  end program convert

