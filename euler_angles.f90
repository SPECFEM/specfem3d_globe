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

! compute the Euler angles and the associated rotation matrix

  subroutine euler_angles(rotation_matrix,ANGULAR_SIZE_CHUNK_RAD)

  implicit none

  include "constants.h"

  double precision rotation_matrix(3,3)
  double precision ANGULAR_SIZE_CHUNK_RAD

  double precision alpha,beta,gamma
  double precision sina,cosa,sinb,cosb,sing,cosg
  double precision center_colatitude,middle_face_colatitude,middle_face_longitude
  double precision center_longitude_loc_rad

! compute colatitute of the center and the corner of the first chunk
  center_colatitude = 90.0d0 - CENTER_LATITUDE_DEG
  middle_face_colatitude = 90.0d0 - MIDDLE_FACE_LATITUDE_DEG

! convert colatitudes and longitude to radians
  center_colatitude = center_colatitude * PI/180.0d0
  center_longitude_loc_rad  = CENTER_LONGITUDE_DEG  * PI/180.0d0
  middle_face_colatitude = middle_face_colatitude * PI/180.0d0

! deduce middle_face_longitude from the three other parameters
! also check particular case of North or South pole
!! DK DK UGLY I think there is a bug below: "PI/4.0D0" assumes that
!! DK DK UGLY the size of the chunk is 90 degrees, but it can now be different
!! DK DK UGLY since the size is defined in constants.h
!! DK DK UGLY therefore we should write something more general here one day
  if(center_colatitude > 0.01d0 .and. center_colatitude < 179.9d0 .and. &
    middle_face_colatitude > 0.01d0 .and. middle_face_colatitude < 179.9d0) then
      middle_face_longitude = center_longitude_loc_rad + dacos((dcos(PI/4.0D0) - &
         dcos(center_colatitude)*dcos(middle_face_colatitude)) / &
         (dsin(center_colatitude)*dsin(middle_face_colatitude)))
  else
    middle_face_longitude = 90.0d0*PI/180.0d0
  endif

  alpha = center_longitude_loc_rad
  beta = center_colatitude

  sing = ((dsin(center_colatitude)*dcos(middle_face_colatitude) - &
           dcos(center_colatitude)*dsin(middle_face_colatitude) * &
           dcos(middle_face_longitude-center_longitude_loc_rad)) &
           / dsin(ANGULAR_SIZE_CHUNK_RAD))

  cosg = dsin(middle_face_colatitude)*dsin(middle_face_longitude-center_longitude_loc_rad)/dsin(ANGULAR_SIZE_CHUNK_RAD)

  if(cosg > -SMALL_VAL_ANGLE .and. cosg <= ZERO) cosg=-SMALL_VAL_ANGLE
  if(cosg < SMALL_VAL_ANGLE .and. cosg >= ZERO) cosg=SMALL_VAL_ANGLE
  gamma=datan2(sing,cosg)
  if(gamma < ZERO) gamma=2.0d0*PI+gamma

  sina = dsin(alpha)
  cosa = dcos(alpha)
  sinb = dsin(beta)
  cosb = dcos(beta)
  sing = dsin(gamma)
  cosg = dcos(gamma)

! define rotation matrix
  rotation_matrix(1,1) = cosg*cosb*cosa-sing*sina
  rotation_matrix(1,2) = -sing*cosb*cosa-cosg*sina
  rotation_matrix(1,3) = sinb*cosa
  rotation_matrix(2,1) = cosg*cosb*sina+sing*cosa
  rotation_matrix(2,2) = -sing*cosb*sina+cosg*cosa
  rotation_matrix(2,3) = sinb*sina
  rotation_matrix(3,1) = -cosg*sinb
  rotation_matrix(3,2) = sing*sinb
  rotation_matrix(3,3) = cosb

  end subroutine euler_angles

