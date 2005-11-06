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

! save some statistics about the mesh

  subroutine save_mesh_statistics(NSPEC_AB,NSPEC_AC,NSPEC_BC, &
        nglob_AB,nglob_AC,nglob_BC,NEX_XI,NPROC,NPROCTOT, &
        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
        INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP)

  implicit none

  include "constants.h"

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               nglob_AB,nglob_AC,nglob_BC

  integer NEX_XI,NPROC,NPROCTOT,NCHUNKS,NSOURCES,NSTEP

  logical INCLUDE_CENTRAL_CUBE

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  integer subtract_central_cube_elems
  double precision subtract_central_cube_points

! for regional code
  double precision x,y,gamma,rgt,xi,eta
  double precision x_top,y_top,z_top
  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! rotation matrix from Euler angles
  integer i,j,ix,iy,icorner
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision r_corner,theta_corner,phi_corner,lat,long,colat_corner

  open(unit=IOUT,file='OUTPUT_FILES/mesh_statistics.txt',status='unknown')

  write(IOUT,*)
  write(IOUT,*) 'mesh statistics:'
  write(IOUT,*) '---------------'
  write(IOUT,*)
  write(IOUT,*) 'number of chunks = ',NCHUNKS
  write(IOUT,*)

! the central cube is counted 6 times, therefore remove 5 times
  if(INCLUDE_CENTRAL_CUBE) then
    write(IOUT,*) 'these statistics include the central cube'
    subtract_central_cube_elems = 5 * (NEX_XI/8)**3
    subtract_central_cube_points = 5.d0 * (dble(NEX_XI/8)*dble(NGLLX-1)+1.d0)**3
  else
    write(IOUT,*) 'these statistics do not include the central cube'
    subtract_central_cube_elems = 0
    subtract_central_cube_points = 0.d0
  endif

  write(IOUT,*)
  write(IOUT,*) 'number of processors = ',NPROCTOT
  write(IOUT,*)
  write(IOUT,*) 'number of ES nodes = ',real(NPROCTOT)/8.
  write(IOUT,*) 'percentage of total 640 ES nodes = ',100.*(real(NPROCTOT)/8.)/640.,' %'
  write(IOUT,*) 'total memory available on these ES nodes (Gb) = ',16.*real(NPROCTOT)/8.
  write(IOUT,*)
  write(IOUT,*) 'max points in largest region = ',nglob_BC(IREGION_CRUST_MANTLE)
  write(IOUT,*)
  write(IOUT,*) 'total elements per AB slice = ',sum(NSPEC_AB)
  write(IOUT,*) 'total points per AB slice = ',sum(nglob_AB)
  write(IOUT,*)
  write(IOUT,*) 'total elements per AC slice = ',sum(NSPEC_AC)
  write(IOUT,*) 'total points per AC slice = ',sum(nglob_AC)
  write(IOUT,*)
  write(IOUT,*) 'total elements per BC slice = ',sum(NSPEC_BC)
  write(IOUT,*) 'total points per BC slice = ',sum(nglob_BC)
  write(IOUT,*)
  write(IOUT,*) 'load balancing AB/BC for points = ',100.*real(sum(nglob_AB))/real(sum(nglob_BC)),' %'
  write(IOUT,*) 'load balancing AB/BC for elements = ',100.*real(sum(NSPEC_AB))/real(sum(NSPEC_BC)),' %'
  write(IOUT,*)
  write(IOUT,*) 'load balancing AC/BC for points = ',100.*real(sum(nglob_AC))/real(sum(nglob_BC)),' %'
  write(IOUT,*) 'load balancing AC/BC for elements = ',100.*real(sum(NSPEC_AC))/real(sum(NSPEC_BC)),' %'
  write(IOUT,*)
  write(IOUT,*) 'total for full 6-chunk mesh:'
  write(IOUT,*) '---------------------------'
  write(IOUT,*)
  write(IOUT,*) 'exact total number of spectral elements in entire mesh = '
  write(IOUT,*) '',2*NPROC*(sum(NSPEC_AB) + sum(NSPEC_AC) + sum(NSPEC_BC)) - subtract_central_cube_elems
  write(IOUT,*) 'approximate total number of points in entire mesh = '
  write(IOUT,*) '',2.d0*dble(NPROC)*(dble(sum(nglob_AB)) + dble(sum(nglob_AC)) + dble(sum(nglob_BC))) &
    - subtract_central_cube_points
! there are 3 DOFs in solid regions, but only 1 in fluid outer core
  write(IOUT,*) 'approximate total number of degrees of freedom in entire mesh = '
  write(IOUT,*) '',2.d0*dble(NPROC)*(3.d0*(dble(sum(nglob_AB)) + dble(sum(nglob_AC)) + dble(sum(nglob_BC))) &
    - 2.d0*dble(nglob_AB(IREGION_OUTER_CORE) + nglob_AC(IREGION_OUTER_CORE) + nglob_BC(IREGION_OUTER_CORE))) &
    - 3.d0*subtract_central_cube_points
  write(IOUT,*)

! display location of chunk if regional run
  if(NCHUNKS /= 6) then

  write(IOUT,*) 'position of the mesh chunk at the surface:'
  write(IOUT,*) '-----------------------------------------'
  write(IOUT,*)
  write(IOUT,*) 'angular size in first direction in degrees = ',sngl(ANGULAR_WIDTH_XI_IN_DEGREES)
  write(IOUT,*) 'angular size in second direction in degrees = ',sngl(ANGULAR_WIDTH_ETA_IN_DEGREES)
  write(IOUT,*)
  write(IOUT,*) 'longitude of center in degrees = ',sngl(CENTER_LONGITUDE_IN_DEGREES)
  write(IOUT,*) 'latitude of center in degrees = ',sngl(CENTER_LATITUDE_IN_DEGREES)
  write(IOUT,*)
  write(IOUT,*) 'angle of rotation of the first chunk = ',sngl(GAMMA_ROTATION_AZIMUTH)

! convert width to radians
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS

! compute rotation matrix from Euler angles
  call euler_angles(rotation_matrix,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH)

! loop on the four corners of the chunk to display their coordinates
  icorner = 0
  do iy = 0,1
    do ix = 0,1

    icorner = icorner + 1

    xi= - ANGULAR_WIDTH_XI_RAD/2. + dble(ix)*ANGULAR_WIDTH_XI_RAD
    eta= - ANGULAR_WIDTH_ETA_RAD/2. + dble(iy)*ANGULAR_WIDTH_ETA_RAD

    x=dtan(xi)
    y=dtan(eta)

    gamma=ONE/dsqrt(ONE+x*x+y*y)
    rgt=R_UNIT_SPHERE*gamma

! define the mesh points at the top surface
    x_top=-y*rgt
    y_top=x*rgt
    z_top=rgt

! rotate top
    vector_ori(1) = x_top
    vector_ori(2) = y_top
    vector_ori(3) = z_top
    do i=1,3
      vector_rotated(i)=0.0d0
      do j=1,3
        vector_rotated(i)=vector_rotated(i)+rotation_matrix(i,j)*vector_ori(j)
      enddo
    enddo
    x_top = vector_rotated(1)
    y_top = vector_rotated(2)
    z_top = vector_rotated(3)

! convert to latitude and longitude
    call xyz_2_rthetaphi_dble(x_top,y_top,z_top,r_corner,theta_corner,phi_corner)
    call reduce(theta_corner,phi_corner)

! convert geocentric to geographic colatitude
    colat_corner=PI/2.0d0-datan(1.006760466d0*dcos(theta_corner)/dmax1(TINYVAL,dsin(theta_corner)))
    if(phi_corner>PI) phi_corner=phi_corner-TWO_PI

! compute real position of the source
    lat = (PI/2.0d0-colat_corner)*180.0d0/PI
    long = phi_corner*180.0d0/PI

    write(IOUT,*)
    write(IOUT,*) 'corner ',icorner
    write(IOUT,*) 'longitude in degrees = ',long
    write(IOUT,*) 'latitude in degrees = ',lat

    enddo
  enddo

  write(IOUT,*)

  endif  ! regional chunk

  write(IOUT,*) 'resolution of the mesh at the surface:'
  write(IOUT,*) '-------------------------------------'
  write(IOUT,*)
  write(IOUT,*) 'spectral elements along a great circle = ',4*NEX_XI
  write(IOUT,*) 'GLL points along a great circle = ',4*NEX_XI*(NGLLX-1)
  write(IOUT,*) 'average distance between points in degrees = ',360./real(4*NEX_XI*(NGLLX-1))
  write(IOUT,*) 'average distance between points in km = ',real(TWO_PI*R_EARTH/1000.d0)/real(4*NEX_XI*(NGLLX-1))
  write(IOUT,*)
  write(IOUT,*) 'number of time steps = ',NSTEP
  write(IOUT,*)
  write(IOUT,*) 'number of seismic sources = ',NSOURCES
  write(IOUT,*)

  close(IOUT)

  end subroutine save_mesh_statistics

