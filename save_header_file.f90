!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC, &
        nglob,NEX_XI,NEX_ETA, &
        nspec_aniso_mantle,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION,ATTENUATION_3D, &
        ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
        INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH,NSOURCES,NSTEP,&
        SIMULATION_TYPE,SAVE_FORWARD,MOVIE_VOLUME)

  implicit none

  include "constants.h"

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC, nglob

  integer NEX_XI,NEX_ETA,NPROC,NPROCTOT,NCHUNKS,NSOURCES,NSTEP
  integer nspec_aniso_mantle,SIMULATION_TYPE

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION,ATTENUATION_3D, &
          INCLUDE_CENTRAL_CUBE,SAVE_FORWARD,MOVIE_VOLUME

  double precision ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES, &
          CENTER_LONGITUDE_IN_DEGREES,CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH

  integer subtract_central_cube_elems
  double precision subtract_central_cube_points

  character(len=150) HEADER_FILE

! for regional code
  double precision x,y,gamma,rgt,xi,eta
  double precision x_top,y_top,z_top
  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! rotation matrix from Euler angles
  integer i,j,ix,iy,icorner
  double precision rotation_matrix(3,3)
  double precision vector_ori(3),vector_rotated(3)
  double precision r_corner,theta_corner,phi_corner,lat,long,colat_corner

! solver's arrays size
  integer :: static_size, dynamic_size

! variables added for memory size computation
  integer :: NSPECMAX_CRUST_MANTLE_ATTENUAT,NSPEC_INNER_CORE_ATTENUATION,NSPECMAX_ISO_MANTLE, &
            NSPECMAX_TISO_MANTLE,NSPECMAX_ANISO_MANTLE,NSPECMAX_ANISO_IC,NSPECMAX_OUTER_CORE_ROTATION

! copy number of elements and points in an include file for the solver
  call get_value_string(HEADER_FILE, 'solver.HEADER_FILE', 'OUTPUT_FILES/values_from_mesher.h')
  open(unit=IOUT,file=HEADER_FILE,status='unknown')
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
    subtract_central_cube_elems = 5 * (NEX_XI/8)**3
    subtract_central_cube_points = 5.d0 * (dble(NEX_XI/8)*dble(NGLLX-1)+1.d0)**3
  else
    write(IOUT,*) '! these statistics do not include the central cube'
    subtract_central_cube_elems = 0
    subtract_central_cube_points = 0.d0
  endif

  write(IOUT,*) '!'
  write(IOUT,*) '! number of processors = ',NPROCTOT
  write(IOUT,*) '!'
  write(IOUT,*) '! number of ES nodes = ',real(NPROCTOT)/8.
  write(IOUT,*) '! percentage of total 640 ES nodes = ',100.*(real(NPROCTOT)/8.)/640.,' %'
  write(IOUT,*) '! total memory available on these ES nodes (Gb) = ',16.*real(NPROCTOT)/8.

  write(IOUT,*) '!'
  write(IOUT,*) '! total points per region = ',nglob(IREGION_CRUST_MANTLE)
! use fused loops on the ES
  write(IOUT,*) '! max vector length = ',nglob(IREGION_CRUST_MANTLE)*NDIM
  write(IOUT,*) '!'
  write(IOUT,*) '! on ES and SX-5, make sure "loopcnt=" parameter'
! use fused loops on the ES
  write(IOUT,*) '! in Makefile is greater than ',nglob(IREGION_CRUST_MANTLE)*NDIM
  write(IOUT,*) '!'

  write(IOUT,*) '! total elements per slice = ',sum(NSPEC)
  write(IOUT,*) '! total points per slice = ',sum(nglob)
  write(IOUT,*) '!'

  write(IOUT,*) '! total for full 6-chunk mesh:'
  write(IOUT,*) '! ---------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! exact total number of spectral elements in entire mesh = '
  write(IOUT,*) '! ',6*NPROC*(sum(NSPEC)) - subtract_central_cube_elems
  write(IOUT,*) '! approximate total number of points in entire mesh = '
  write(IOUT,*) '! ',2.d0*dble(NPROC)*(3.d0*dble(sum(nglob))) - subtract_central_cube_points
! there are 3 DOFs in solid regions, but only 1 in fluid outer core
  write(IOUT,*) '! approximate total number of degrees of freedom in entire mesh = '
  write(IOUT,*) '! ',6.d0*dble(NPROC)*(3.d0*(dble(sum(nglob))) &
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
    colat_corner=PI/2.0d0-datan(1.006760466d0*dcos(theta_corner)/dmax1(TINYVAL,dsin(theta_corner)))
    if(phi_corner>PI) phi_corner=phi_corner-TWO_PI

! compute real position of the source
    lat = (PI/2.0d0-colat_corner)*180.0d0/PI
    long = phi_corner*180.0d0/PI

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
  write(IOUT,*) '!'
  write(IOUT,*) '! number of time steps = ',NSTEP
  write(IOUT,*) '!'
  write(IOUT,*) '! number of seismic sources = ',NSOURCES
  write(IOUT,*) '!'
  write(IOUT,*)

! DM memory size evaluation
    static_size = 0
    dynamic_size = 0

  if(ATTENUATION) then
    if(ATTENUATION_3D) then
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
! factor_scale_crust_mantle,one_minus_sum_beta_crust_mantle,factor_common_crust_mantle
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*3*CUSTOM_REAL
! factor_scale_inner_core, one_minus_sum_beta_inner_core, factor_common_inner_core
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*3*SIZE_DOUBLE
! factor_scale_crust_mantle_dble, omsb_crust_mantle_dble, factor_common_crust_mantle_dble
      dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*3*SIZE_DOUBLE
! factor_scale_inner_core_dble, omsb_inner_core_dble, factor_common_inner_core_dble
    endif
    NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPEC(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_ATTENUATION = NSPEC(IREGION_INNER_CORE)
  else
    NSPECMAX_CRUST_MANTLE_ATTENUAT = 1
    NSPEC_INNER_CORE_ATTENUATION = 1
  endif
  if(ANISOTROPIC_3D_MANTLE) then
    NSPECMAX_ISO_MANTLE = 1
    NSPECMAX_TISO_MANTLE = 1
    NSPECMAX_ANISO_MANTLE = NSPEC(IREGION_CRUST_MANTLE)
  else
    NSPECMAX_ISO_MANTLE = NSPEC(IREGION_CRUST_MANTLE)
    if(TRANSVERSE_ISOTROPY) then
      NSPECMAX_TISO_MANTLE = nspec_aniso_mantle
    else
      NSPECMAX_TISO_MANTLE = 1
    endif
    NSPECMAX_ANISO_MANTLE = 1
  endif
  if(ANISOTROPIC_INNER_CORE) then
    NSPECMAX_ANISO_IC = NSPEC(IREGION_INNER_CORE)
  else
    NSPECMAX_ANISO_IC = 1
  endif
  if(ROTATION) then
    NSPECMAX_OUTER_CORE_ROTATION = NSPEC(IREGION_OUTER_CORE)
  else
    NSPECMAX_OUTER_CORE_ROTATION = 1
  endif

! size of static arrays
! ---------------------

    static_size = static_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPECMAX_CRUST_MANTLE_ATTENUAT*CUSTOM_REAL   !R_memory_crust_mantle
    static_size = static_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION*CUSTOM_REAL     !R_memory_inner_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*12*CUSTOM_REAL
! rho_vp_crust_mantle,rho_vs_crust_mantle,xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle
! etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,
! gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,jacobian_crust_mantle
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*SIZE_REAL
! ibool_crust_mantle
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*13*CUSTOM_REAL
! vp_outer_core,xix_outer_core,xiy_outer_core,xiz_outer_core,
! etax_outer_core,etay_outer_core,etaz_outer_core,
! gammax_outer_core,gammay_outer_core,gammaz_outer_core,jacobian_outer_core
! rhostore_outer_core,kappavstore_outer_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*SIZE_REAL
! ibool_outer_core
    static_size = static_size + nglob(IREGION_CRUST_MANTLE)*2*SIZE_REAL
! updated_dof_ocean_load, idoubling_crust_mantle
    static_size = static_size + nglob(IREGION_CRUST_MANTLE)*5*CUSTOM_REAL
! xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
! rmass_crust_mantle,rmass_ocean_load
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_ISO_MANTLE*3*CUSTOM_REAL
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_TISO_MANTLE*3*CUSTOM_REAL
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_ANISO_MANTLE*21*CUSTOM_REAL
    static_size = static_size + NDIM*nglob(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
! displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
    static_size = static_size + nglob(IREGION_OUTER_CORE)*7*CUSTOM_REAL
! xstore_outer_core, ystore_outer_core, zstore_outer_core, rmass_outer_core, displ_outer_core, veloc_outer_core, accel_outer_core

    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*SIZE_REAL
! ibool_inner_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*13*CUSTOM_REAL
! xix_inner_core,xiy_inner_core,xiz_inner_core,
! etax_inner_core,etay_inner_core,etaz_inner_core,
! gammax_inner_core,gammay_inner_core,gammaz_inner_core,jacobian_inner_core,
! rhostore_inner_core, kappavstore_inner_core,muvstore_inner_core
    static_size = static_size + nglob(IREGION_INNER_CORE)*4*CUSTOM_REAL
! xstore_inner_core,ystore_inner_core,zstore_inner_core,rmass_inner_core
    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_ANISO_IC*5*CUSTOM_REAL
! c11store_inner_core,c33store_inner_core,c12store_inner_core,
! c13store_inner_core,c44store_inner_core
    static_size = static_size + NDIM*nglob(IREGION_INNER_CORE)*3*CUSTOM_REAL
! displ_inner_core,veloc_inner_core,accel_inner_core

    static_size = static_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_OUTER_CORE_ROTATION*2*CUSTOM_REAL
! A_array_rotation,B_array_rotation

! size of dynamic arrays
! ----------------------

    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*2 !dummy_rho,dummy_vstore
  if (SIMULATION_TYPE > 1) then
    dynamic_size = dynamic_size + NDIM*nglob(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + nglob(IREGION_OUTER_CORE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NDIM*nglob(IREGION_INNER_CORE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*2*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_OUTER_CORE)*2*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*3*CUSTOM_REAL
    dynamic_size = dynamic_size + NGLLX*NGLLY*NGLLZ*NSPECMAX_OUTER_CORE_ROTATION*2*CUSTOM_REAL
  endif
  if (ATTENUATION .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    dynamic_size = dynamic_size + 7*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*CUSTOM_REAL
    dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*CUSTOM_REAL
    if (SIMULATION_TYPE > 1) then
      dynamic_size = dynamic_size + 7*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_CRUST_MANTLE)*CUSTOM_REAL
      dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPEC(IREGION_INNER_CORE)*CUSTOM_REAL
    endif
  else
    dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPECMAX_CRUST_MANTLE_ATTENUAT*CUSTOM_REAL
    dynamic_size = dynamic_size + 5*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION*CUSTOM_REAL
  endif
  if(ATTENUATION .and. SIMULATION_TYPE > 1) then
    dynamic_size = dynamic_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPECMAX_CRUST_MANTLE_ATTENUAT*CUSTOM_REAL
    dynamic_size = dynamic_size + 5*N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION*CUSTOM_REAL
  endif

    write(IOUT,*) '! approximate memory needed for the solver : '
    write(IOUT,*) '!'
    write(IOUT,*) '! size of static arrays per slice : ',static_size/(1024**2),' MB'
    write(IOUT,*) '! size of static arrays for all slices : ',((static_size/(1024**2))*NPROCTOT)/1024.d0,' GB'
    write(IOUT,*) '!'
    write(IOUT,*) '! size of dynamic arrays per slice : ',dynamic_size/(1024**2),' MB'
    write(IOUT,*) '! size of dynamic arrays for all slices : ',((dynamic_size/(1024**2))*NPROCTOT)/1024.d0,' GB'
    write(IOUT,*) '!'
    write(IOUT,*) '! total size of arrays per slice : ',(dynamic_size+static_size)/(1024**2),' MB'
    write(IOUT,*) '! total size of arrays for all slices : ', &
                    (((dynamic_size+static_size)/(1024**2))*NPROCTOT)/1024.d0,' GB'
    write(IOUT,*)

! DM memory size evaluation

  if(NCHUNKS == 1) write(IOUT,*) '! values for AC and BC below undefined for one chunk'
  if(NCHUNKS == 2) write(IOUT,*) '! values for BC below undefined for two chunks'

  write(IOUT,*) 'integer, parameter :: NEX_XI_VAL = ',NEX_XI
  write(IOUT,*) 'integer, parameter :: NEX_ETA_VAL = ',NEX_ETA
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_AB = ',NSPEC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_AC = ',NSPEC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_BC = ',NSPEC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_AB = ',NSPEC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_AC = ',NSPEC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_BC = ',NSPEC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE = ',NSPEC(IREGION_INNER_CORE)
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_AB = ',nglob(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_AC = ',nglob(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_BC = ',nglob(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_AB = ',nglob(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_AC = ',nglob(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_BC = ',nglob(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_INNER_CORE = ',nglob(IREGION_INNER_CORE)
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC'
  write(IOUT,*) 'integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC'
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC'
  write(IOUT,*) 'integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC'
  write(IOUT,*)

  if(ANISOTROPIC_INNER_CORE) then
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_IC = ',NSPEC(IREGION_INNER_CORE)
  else
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_IC = 1'
  endif

  if(ANISOTROPIC_3D_MANTLE) then
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ISO_MANTLE = ',1
    write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',1
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_MANTLE = NSPECMAX_CRUST_MANTLE'
  else

    write(IOUT,*) 'integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE'
    if(TRANSVERSE_ISOTROPY) then
      write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',nspec_aniso_mantle
    else
      write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',1
    endif

    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_MANTLE = 1'
  endif

  write(IOUT,*)

! if attenuation is off, set dummy size of arrays to one
  if(ATTENUATION) then
    write(IOUT,*) 'integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPECMAX_CRUST_MANTLE'
    write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = NSPEC_INNER_CORE'
  else
    write(IOUT,*) 'integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = 1'
    write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = 1'
  endif

! this to allow for code elimination by compiler in solver for performance
  write(IOUT,*)

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
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL_3D = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL_3D = .false.'
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

  if(ROTATION) then
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .true.'
    write(IOUT,*) 'integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = NSPECMAX_OUTER_CORE'
  else
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .false.'
    write(IOUT,*) 'integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = 1'
  endif
  write(IOUT,*)

  close(IOUT)

  end subroutine save_header_file

