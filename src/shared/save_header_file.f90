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

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC,nglob,NEX_XI,NEX_ETA,NPROC,NPROCTOT, &
                        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
                        ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY, &
                        OCEANS,ATTENUATION,ATTENUATION_3D, &
                        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
                        INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES, &
                        CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP,&
                        static_memory_size, &
                        NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                        NSPEC2D_TOP,NSPEC2D_BOTTOM, &
                        NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX, &
                        NPROC_XI,NPROC_ETA, &
                        NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                        NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
                        NSPEC_INNER_CORE_ATTENUATION, &
                        NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                        NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                        NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                        NSPEC_CRUST_MANTLE_ADJOINT, &
                        NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                        NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                        NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
                        NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                        NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION, &
                        SIMULATION_TYPE,SAVE_FORWARD,MOVIE_VOLUME,SAVE_REGULAR_KL,NOISE_TOMOGRAPHY, &
                        ATT1,ATT2,ATT3,ATT4,ATT5,APPROXIMATE_HESS_KL,ANISOTROPIC_KL,PARTIAL_PHYS_DISPERSION_ONLY)

  implicit none

  include "constants.h"

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, nglob

  integer NEX_XI,NEX_ETA,NPROC,NPROCTOT,NCHUNKS,NSOURCES,NSTEP,NOISE_TOMOGRAPHY

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          ELLIPTICITY,GRAVITY,ROTATION,TOPOGRAPHY,OCEANS,ATTENUATION,ATTENUATION_3D,INCLUDE_CENTRAL_CUBE, &
          SAVE_REGULAR_KL,APPROXIMATE_HESS_KL,ANISOTROPIC_KL,PARTIAL_PHYS_DISPERSION_ONLY

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  ! static memory size needed by the solver
  double precision :: static_memory_size

  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
        NSPEC2D_TOP,NSPEC2D_BOTTOM,NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX
  integer :: NPROC_XI,NPROC_ETA
  integer :: NCORNERSCHUNKS,NUM_FACES,NUM_MSG_TYPES

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUAT, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION, &
         NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670, NSPEC2D_CMB, NSPEC2D_ICB

  integer :: SIMULATION_TYPE
  logical :: SAVE_FORWARD,MOVIE_VOLUME

  ! local parameters
  double precision :: subtract_central_cube_elems,subtract_central_cube_points
  ! for regional code
  double precision x,y,gamma,rgt,xi,eta
  double precision x_top,y_top,z_top
  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD
  ! rotation matrix from Euler angles
  integer i,j,ix,iy,icorner
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision r_corner,theta_corner,phi_corner,lat,long,colat_corner
  integer :: ATT1,ATT2,ATT3,ATT4,ATT5
  integer :: ier
  character(len=150) HEADER_FILE

  ! copy number of elements and points in an include file for the solver
  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')
  open(unit=IOUT,file=HEADER_FILE,status='unknown',iostat=ier)
  if( ier /= 0 ) stop 'error opening OUTPUT_FILES/values_from_mesher.h'

  write(IOUT,*)

  write(IOUT,*) '!'
  write(IOUT,*) '! this is the parameter file for static compilation of the solver'
  write(IOUT,*) '!'
  write(IOUT,*) '! mesh statistics:'
  write(IOUT,*) '! ---------------'
  write(IOUT,*) '!'
  write(IOUT,*) '!'
  write(IOUT,*) '! number of chunks = ',NCHUNKS
  write(IOUT,*) '!'

! the central cube is counted 6 times, therefore remove 5 times
  if(INCLUDE_CENTRAL_CUBE) then
    write(IOUT,*) '! these statistics include the central cube'
    subtract_central_cube_elems = 5.d0 * dble((NEX_XI/8))**3
    subtract_central_cube_points = 5.d0 * (dble(NEX_XI/8)*dble(NGLLX-1)+1.d0)**3
  else
    write(IOUT,*) '! these statistics do not include the central cube'
    subtract_central_cube_elems = 0.d0
    subtract_central_cube_points = 0.d0
  endif

  write(IOUT,*) '!'
  write(IOUT,*) '! number of processors = ',NPROCTOT ! should be = NPROC
  write(IOUT,*) '!'
  write(IOUT,*) '! maximum number of points per region = ',nglob(IREGION_CRUST_MANTLE)
  write(IOUT,*) '!'
! use fused loops on NEC SX
  write(IOUT,*) '! on NEC SX, make sure "loopcnt=" parameter'
  write(IOUT,*) '! in Makefile is greater than max vector length = ',nglob(IREGION_CRUST_MANTLE)*NDIM
  write(IOUT,*) '!'

  write(IOUT,*) '! total elements per slice = ',sum(NSPEC)
  write(IOUT,*) '! total points per slice = ',sum(nglob)
  write(IOUT,*) '!'

  write(IOUT,'(1x,a,i1,a)') '! total for full ',NCHUNKS,'-chunk mesh:'
  write(IOUT,*) '! ---------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! exact total number of spectral elements in entire mesh = '
  write(IOUT,*) '! ',dble(NCHUNKS)*dble(NPROC)*dble(sum(NSPEC)) - subtract_central_cube_elems
  write(IOUT,*) '! approximate total number of points in entire mesh = '
  write(IOUT,*) '! ',dble(NCHUNKS)*dble(NPROC)*dble(sum(nglob)) - subtract_central_cube_points
! there are 3 DOFs in solid regions, but only 1 in fluid outer core
  write(IOUT,*) '! approximate total number of degrees of freedom in entire mesh = '
  write(IOUT,*) '! ',dble(NCHUNKS)*dble(NPROC)*(3.d0*(dble(sum(nglob))) &
    - 2.d0*dble(nglob(IREGION_OUTER_CORE))) &
    - 3.d0*subtract_central_cube_points
  write(IOUT,*) '!'

! display location of chunk if regional run
  if(NCHUNKS /= 6) then

  write(IOUT,*) '! position of the mesh chunk at the surface:'
  write(IOUT,*) '! -----------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! angular size in first direction in degrees = ',sngl(ANGULAR_WIDTH_XI_IN_DEGREES)
  write(IOUT,*) '! angular size in second direction in degrees = ',sngl(ANGULAR_WIDTH_ETA_IN_DEGREES)
  write(IOUT,*) '!'
  write(IOUT,*) '! longitude of center in degrees = ',sngl(CENTER_LONGITUDE_IN_DEGREES)
  write(IOUT,*) '! latitude of center in degrees = ',sngl(CENTER_LATITUDE_IN_DEGREES)
  write(IOUT,*) '!'
  write(IOUT,*) '! angle of rotation of the first chunk = ',sngl(GAMMA_ROTATION_AZIMUTH)

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
    colat_corner=PI_OVER_TWO-datan(1.006760466d0*dcos(theta_corner)/dmax1(TINYVAL,dsin(theta_corner)))
    if(phi_corner>PI) phi_corner=phi_corner-TWO_PI

! compute real position of the source
    lat = (PI_OVER_TWO-colat_corner)*RADIANS_TO_DEGREES
    long = phi_corner*RADIANS_TO_DEGREES

    write(IOUT,*) '!'
    write(IOUT,*) '! corner ',icorner
    write(IOUT,*) '! longitude in degrees = ',long
    write(IOUT,*) '! latitude in degrees = ',lat

    enddo
  enddo

  write(IOUT,*) '!'

  endif  ! regional chunk

  write(IOUT,*) '! resolution of the mesh at the surface:'
  write(IOUT,*) '! -------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! spectral elements along a great circle = ',4*NEX_XI
  write(IOUT,*) '! GLL points along a great circle = ',4*NEX_XI*(NGLLX-1)
  write(IOUT,*) '! average distance between points in degrees = ',360./real(4*NEX_XI*(NGLLX-1))
  write(IOUT,*) '! average distance between points in km = ',real(TWO_PI*R_EARTH/1000.d0)/real(4*NEX_XI*(NGLLX-1))
  write(IOUT,*) '! average size of a spectral element in km = ',real(TWO_PI*R_EARTH/1000.d0)/real(4*NEX_XI)
  write(IOUT,*) '!'
  write(IOUT,*) '! number of time steps = ',NSTEP
  write(IOUT,*) '!'
  write(IOUT,*) '! number of seismic sources = ',NSOURCES
  write(IOUT,*) '!'
  write(IOUT,*)

  write(IOUT,*) '! approximate static memory needed by the solver:'
  write(IOUT,*) '! ----------------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! (lower bound, usually the real amount used is 5% to 10% higher)'
  write(IOUT,*) '!'
  write(IOUT,*) '! (you can get a more precise estimate of the size used per MPI process'
  write(IOUT,*) '!  by typing "size -d bin/xspecfem3D"'
  write(IOUT,*) '!  after compiling the code with the DATA/Par_file you plan to use)'
  write(IOUT,*) '!'
  write(IOUT,*) '! size of static arrays per slice = ',static_memory_size/1.d6,' MB'
  write(IOUT,*) '!                                 = ',static_memory_size/1048576.d0,' MiB'
  write(IOUT,*) '!                                 = ',static_memory_size/1.d9,' GB'
  write(IOUT,*) '!                                 = ',static_memory_size/1073741824.d0,' GiB'
  write(IOUT,*) '!'

  ! note: using less memory becomes an issue only if the strong scaling of the code is poor.
  !          Some users will run simulations with an executable using far less than 80% RAM per core
  !          if they prefer having a faster computational time (and use a higher number of cores).

  write(IOUT,*) '! (should be below to 80% or 90% of the memory installed per core)'
  write(IOUT,*) '! (if significantly more, the job will not run by lack of memory )'
  write(IOUT,*) '! (note that if significantly less, you waste a significant amount'
  write(IOUT,*) '!  of memory per processor core)'
  write(IOUT,*) '! (but that can be perfectly acceptable if you can afford it and'
  write(IOUT,*) '!  want faster results by using more cores)'
  write(IOUT,*) '!'
  if(static_memory_size*dble(NPROCTOT)/1.d6 < 10000.d0) then
    write(IOUT,*) '! size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d6,' MB'
    write(IOUT,*) '!                                      = ',static_memory_size*dble(NPROCTOT)/1048576.d0,' MiB'
    write(IOUT,*) '!                                      = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
  else
    write(IOUT,*) '! size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
  endif
  write(IOUT,*) '!                                      = ',static_memory_size*dble(NPROCTOT)/1073741824.d0,' GiB'
  write(IOUT,*) '!                                      = ',static_memory_size*dble(NPROCTOT)/1.d12,' TB'
  write(IOUT,*) '!                                      = ',static_memory_size*dble(NPROCTOT)/1099511627776.d0,' TiB'
  write(IOUT,*) '!'

  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NEX_XI_VAL = ',NEX_XI
  write(IOUT,*) 'integer, parameter :: NEX_ETA_VAL = ',NEX_ETA
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE = ',NSPEC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE = ',NSPEC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE = ',NSPEC(IREGION_INNER_CORE)
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE = ',nglob(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE = ',nglob(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_INNER_CORE = ',nglob(IREGION_INNER_CORE)
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_IC = ',NSPECMAX_ANISO_IC
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPECMAX_ISO_MANTLE = ',NSPECMAX_ISO_MANTLE
  write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',NSPECMAX_TISO_MANTLE
  write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_MANTLE = ',NSPECMAX_ANISO_MANTLE
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUAT = ',NSPEC_CRUST_MANTLE_ATTENUAT
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = ',NSPEC_INNER_CORE_ATTENUATION
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_STR_OR_ATT = ',NSPEC_CRUST_MANTLE_STR_OR_ATT
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_STR_OR_ATT = ',NSPEC_INNER_CORE_STR_OR_ATT
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_STR_AND_ATT = ',NSPEC_CRUST_MANTLE_STR_AND_ATT
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_STR_AND_ATT = ',NSPEC_INNER_CORE_STR_AND_ATT
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_STRAIN_ONLY = ',NSPEC_CRUST_MANTLE_STRAIN_ONLY
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_STRAIN_ONLY = ',NSPEC_INNER_CORE_STRAIN_ONLY
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT = ',NSPEC_CRUST_MANTLE_ADJOINT
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_ADJOINT = ',NSPEC_OUTER_CORE_ADJOINT
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_ADJOINT = ',NSPEC_INNER_CORE_ADJOINT

  if(ANISOTROPIC_KL) then
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL = ',NSPEC_CRUST_MANTLE_ADJOINT
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL = ',1
  endif

  if(APPROXIMATE_HESS_KL) then
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_HESS = ',NSPEC_CRUST_MANTLE_ADJOINT
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_HESS = ',1
  endif

  if(NOISE_TOMOGRAPHY > 0) then
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_NOISE = ',NSPEC_CRUST_MANTLE_ADJOINT
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_NOISE = ',1
  endif

  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_ADJOINT = ',NGLOB_CRUST_MANTLE_ADJOINT
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_ADJOINT = ',NGLOB_OUTER_CORE_ADJOINT
  write(IOUT,*) 'integer, parameter :: NGLOB_INNER_CORE_ADJOINT = ',NGLOB_INNER_CORE_ADJOINT

  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_ROT_ADJOINT = ',NSPEC_OUTER_CORE_ROT_ADJOINT
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_STACEY = ',NSPEC_CRUST_MANTLE_STACEY
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_STACEY = ',NSPEC_OUTER_CORE_STACEY
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_OCEANS = ',NGLOB_CRUST_MANTLE_OCEANS
  write(IOUT,*)

! this to allow for code elimination by compiler in solver for performance

  if(TRANSVERSE_ISOTROPY) then
    write(IOUT,*) 'logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.'
  endif
  write(IOUT,*)

  if(ANISOTROPIC_3D_MANTLE) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.'
  endif
  write(IOUT,*)

  if(ANISOTROPIC_INNER_CORE) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.'
  endif
  write(IOUT,*)

  if(ATTENUATION) then
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .false.'
  endif
  write(IOUT,*)

  if(ATTENUATION_3D) then
    write(IOUT,*) 'logical, parameter :: ATTENUATION_3D_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_3D_VAL = .false.'
  endif
  write(IOUT,*)

  if(ELLIPTICITY) then
    write(IOUT,*) 'logical, parameter :: ELLIPTICITY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ELLIPTICITY_VAL = .false.'
  endif
  write(IOUT,*)

  if(GRAVITY) then
    write(IOUT,*) 'logical, parameter :: GRAVITY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: GRAVITY_VAL = .false.'
  endif
  write(IOUT,*)

  if(OCEANS) then
    write(IOUT,*) 'logical, parameter :: OCEANS_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: OCEANS_VAL = .false.'
  endif
  write(IOUT,*)

  if(TOPOGRAPHY .or. OCEANS) then
    write(IOUT,*) 'integer, parameter :: NX_BATHY_VAL = NX_BATHY'
    write(IOUT,*) 'integer, parameter :: NY_BATHY_VAL = NY_BATHY'
  else
    write(IOUT,*) 'integer, parameter :: NX_BATHY_VAL = 1'
    write(IOUT,*) 'integer, parameter :: NY_BATHY_VAL = 1'
  endif
  write(IOUT,*)

  if(ROTATION) then
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .false.'
  endif
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_ROTATION = ',NSPEC_OUTER_CORE_ROTATION
  write(IOUT,*)

  if(PARTIAL_PHYS_DISPERSION_ONLY) then
    write(IOUT,*) 'logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.'
  endif

  write(IOUT,*) 'integer, parameter :: NGLOB1D_RADIAL_CM = ',NGLOB1D_RADIAL(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB1D_RADIAL_OC = ',NGLOB1D_RADIAL(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB1D_RADIAL_IC = ',NGLOB1D_RADIAL(IREGION_INNER_CORE)

  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_XMIN_XMAX_CM = ',NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_XMIN_XMAX_OC = ',NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_XMIN_XMAX_IC = ',NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE)

  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_YMIN_YMAX_CM = ',NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_YMIN_YMAX_OC = ',NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_YMIN_YMAX_IC = ',NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE)

  write(IOUT,*) 'integer, parameter :: NPROC_XI_VAL = ',NPROC_XI
  write(IOUT,*) 'integer, parameter :: NPROC_ETA_VAL = ',NPROC_ETA
  write(IOUT,*) 'integer, parameter :: NCHUNKS_VAL = ',NCHUNKS
  write(IOUT,*) 'integer, parameter :: NPROCTOT_VAL = ',NPROCTOT
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_XY_CM_VAL = ', &
            max(NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE))
  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_XY_OC_VAL = ', &
            max(NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE))
  write(IOUT,*) 'integer, parameter :: NGLOB2DMAX_XY_IC_VAL = ', &
            max(NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE))

  if(NCHUNKS == 1 .or. NCHUNKS == 2) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 1
  else if(NCHUNKS == 3) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 3
  else if(NCHUNKS == 6) then
    NCORNERSCHUNKS = 8
    NUM_FACES = 4
    NUM_MSG_TYPES = 3
  else
    stop 'error nchunks in save_header_file()'
  endif

  write(IOUT,*) 'integer, parameter :: NUMMSGS_FACES_VAL = ',NPROC_XI*NUM_FACES*NUM_MSG_TYPES
  write(IOUT,*) 'integer, parameter :: NCORNERSCHUNKS_VAL = ',NCORNERSCHUNKS

  write(IOUT,*) 'integer, parameter :: ATT1_VAL = ',ATT1
  write(IOUT,*) 'integer, parameter :: ATT2_VAL = ',ATT2
  write(IOUT,*) 'integer, parameter :: ATT3_VAL = ',ATT3
  write(IOUT,*) 'integer, parameter :: ATT4_VAL = ',ATT4
  write(IOUT,*) 'integer, parameter :: ATT5_VAL = ',ATT5
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_XMIN_XMAX_CM = ',NSPEC2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_YMIN_YMAX_CM = ',NSPEC2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC2D_BOTTOM_CM = ',NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC2D_TOP_CM = ',NSPEC2D_TOP(IREGION_CRUST_MANTLE)

  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_XMIN_XMAX_IC = ',NSPEC2DMAX_XMIN_XMAX(IREGION_INNER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_YMIN_YMAX_IC = ',NSPEC2DMAX_YMIN_YMAX(IREGION_INNER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC2D_BOTTOM_IC = ',NSPEC2D_BOTTOM(IREGION_INNER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC2D_TOP_IC = ',NSPEC2D_TOP(IREGION_INNER_CORE)

  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_XMIN_XMAX_OC = ',NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC2DMAX_YMIN_YMAX_OC = ',NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC2D_BOTTOM_OC = ',NSPEC2D_BOTTOM(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC2D_TOP_OC = ',NSPEC2D_TOP(IREGION_OUTER_CORE)

  ! for boundary kernels

  if (SAVE_BOUNDARY_MESH) then
    NSPEC2D_MOHO = NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    NSPEC2D_400 = NSPEC2D_MOHO / 4
    NSPEC2D_670 = NSPEC2D_400
    NSPEC2D_CMB = NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)
    NSPEC2D_ICB = NSPEC2D_BOTTOM(IREGION_OUTER_CORE)
  else
    NSPEC2D_MOHO = 1
    NSPEC2D_400 = 1
    NSPEC2D_670 = 1
    NSPEC2D_CMB = 1
    NSPEC2D_ICB = 1
  endif

  write(IOUT,*) 'integer, parameter :: NSPEC2D_MOHO = ',NSPEC2D_MOHO
  write(IOUT,*) 'integer, parameter :: NSPEC2D_400 = ',NSPEC2D_400
  write(IOUT,*) 'integer, parameter :: NSPEC2D_670 = ',NSPEC2D_670
  write(IOUT,*) 'integer, parameter :: NSPEC2D_CMB = ',NSPEC2D_CMB
  write(IOUT,*) 'integer, parameter :: NSPEC2D_ICB = ',NSPEC2D_ICB
  write(IOUT,*)

  ! deville routines only implemented for NGLLX = NGLLY = NGLLZ = 5
  if( NGLLX == 5 .and. NGLLY == 5 .and. NGLLZ == 5 ) then
    write(IOUT,*) 'logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .false.'
  endif

  ! attenuation and/or adjoint simulations
  if (ATTENUATION .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD &
    .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    write(IOUT,*) 'logical, parameter :: COMPUTE_AND_STORE_STRAIN_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: COMPUTE_AND_STORE_STRAIN_VAL = .false.'
  endif
  write(IOUT,*)

  if (MOVIE_VOLUME) then
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = NSPEC_CRUST_MANTLE'
    write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = NGLOB_CRUST_MANTLE'
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1'
    write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1'
  endif
  write(IOUT,*)

  if (SAVE_REGULAR_KL) then
    write(IOUT,*) 'integer, parameter :: NM_KL_REG_PTS_VAL = NM_KL_REG_PTS'
  else
    write(IOUT,*) 'integer, parameter :: NM_KL_REG_PTS_VAL = 1'
  endif

  close(IOUT)

  end subroutine save_header_file

