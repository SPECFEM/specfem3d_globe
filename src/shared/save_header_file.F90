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

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC_REGIONS,NGLOB_REGIONS,NPROC,NPROCTOT, &
                              static_memory_size, &
                              NSPEC2D_TOP,NSPEC2D_BOTTOM, &
                              NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX, &
                              NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                              NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
                              NSPEC_INNER_CORE_ATTENUATION, &
                              NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                              NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                              NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                              NSPEC_CRUST_MANTLE_ADJOINT, &
                              NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                              NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
                              NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
                              NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                              NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION,NT_DUMP_ATTENUATION_optimal, &
                              PRINT_INFO_TO_SCREEN)

! Daniel note: the comment below is wrong, since e.g. NSPEC_CRUST_MANTLE_ADJOINT is either 1 (dummy value)
!              for SIMULATION_TYPE == 1 or equal to NSPEC_CRUST_MANTLE for SIMULATION_TYPE == 3 or SAVE_FORWARD set to .true.
!              the value is determined in routine memory_eval() and passed here as argument.
!
!              thus, in order to be able to run a kernel simulation, users still have to re-compile the binaries with the right
!              parameters set in Par_file otherwise the kernel runs will crash.
!
!              there is hardly another way of this usage when trying to take advantage of static compilation for the solver.
!              please ignore the following remark:
! ****************************************************************************************************
! IMPORTANT: this routine must *NOT* use flag SIMULATION_TYPE (nor SAVE_FORWARD), i.e. none of the parameters it computes
! should depend on SIMULATION_TYPE, because most users do not recompile the code nor rerun the mesher
! when switching from SIMULATION_TYPE == 1 to SIMULATION_TYPE == 3 and thus the header file created
! by this routine would become wrong in the case of a run with SIMULATION_TYPE == 3 if the code
! was compiled with SIMULATION_TYPE == 1
! ****************************************************************************************************

  use constants

  use shared_parameters, only: TOPOGRAPHY, &
    TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
    ELLIPTICITY,GRAVITY,ROTATION, &
    OCEANS,ATTENUATION,ATTENUATION_3D, &
    ANGULAR_WIDTH_XI_IN_DEGREES,ANGULAR_WIDTH_ETA_IN_DEGREES,NCHUNKS, &
    INCLUDE_CENTRAL_CUBE,CENTER_LONGITUDE_IN_DEGREES, &
    CENTER_LATITUDE_IN_DEGREES,GAMMA_ROTATION_AZIMUTH, &
    DT,NEX_XI,NEX_ETA, &
    NPROC_XI,NPROC_ETA, &
    PARTIAL_PHYS_DISPERSION_ONLY, &
    ABSORBING_CONDITIONS,EXACT_MASS_MATRIX_FOR_ROTATION, &
    ATT1,ATT2,ATT3,ATT4,ATT5, &
    MOVIE_VOLUME,MOVIE_VOLUME_TYPE,NTSTEP_BETWEEN_FRAMES,SIMULATION_TYPE,MOVIE_SURFACE, &
    UNDO_ATTENUATION,MEMORY_INSTALLED_PER_CORE_IN_GB,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
    this_region_has_a_doubling,doubling_index,ner,ratio_sampling_array, &
    RECORD_LENGTH_IN_MINUTES,NSTEP

  implicit none

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_REGIONS,NGLOB_REGIONS

  integer :: NPROC,NPROCTOT,NT_DUMP_ATTENUATION_optimal

  ! static memory size needed by the solver
  double precision :: static_memory_size

  integer, dimension(MAX_NUM_REGIONS) :: &
        NSPEC2D_TOP,NSPEC2D_BOTTOM,NSPEC2DMAX_YMIN_YMAX,NSPEC2DMAX_XMIN_XMAX

  integer :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
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

  integer :: NGLOB_XY_CM, NGLOB_XY_IC

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
  integer :: ier

  integer :: num_elem_gc,num_gll_gc
  double precision :: avg_dist_deg,avg_dist_km,avg_element_size

!! DK DK for UNDO_ATTENUATION
  integer :: saved_SIMULATION_TYPE
  integer :: number_of_dumpings_to_do
  double precision :: static_memory_size_GB,size_to_store_at_each_time_step,disk_size_of_each_dumping

  logical :: PRINT_INFO_TO_SCREEN

! evaluate the amount of static memory needed by the solver
  call memory_eval(doubling_index,this_region_has_a_doubling, &
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                   ratio_sampling_array,NPROCTOT,NSPEC_REGIONS,NGLOB_REGIONS, &
                   NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                   NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
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
                   NSPEC2D_BOTTOM,NSPEC2D_TOP,static_memory_size)

  if (PRINT_INFO_TO_SCREEN) then

    print *
    print *,'edit file OUTPUT_FILES/values_from_mesher.h to see'
    print *,'some statistics about the mesh'
    print *

    print *,'number of processors = ',NPROCTOT
    print *
    print *,'maximum number of points per region = ',NGLOB_REGIONS(IREGION_CRUST_MANTLE)
    print *
    print *,'total elements per slice = ',sum(NSPEC_REGIONS)
    print *,'total points per slice = ',sum(NGLOB_REGIONS)
    print *
    print *,'the time step of the solver will be DT = ',sngl(DT),' (s)'
    print *,'the (approximate) minimum period resolved will be = ', &
            sngl(max(ANGULAR_WIDTH_ETA_IN_DEGREES,ANGULAR_WIDTH_XI_IN_DEGREES)/90.0 * 256.0/min(NEX_ETA,NEX_XI) * 17.0),' (s)'
    print *
    print *,'current record length is = ',sngl(RECORD_LENGTH_IN_MINUTES),'min'
    print *,'current minimum number of time steps will be = ',NSTEP
    print *
    if (MOVIE_SURFACE .or. MOVIE_VOLUME) then
      print *,'MOVIE_VOLUME :',MOVIE_VOLUME
      print *,'MOVIE_SURFACE:',MOVIE_SURFACE
      print *,'Saving movie frames every',NTSTEP_BETWEEN_FRAMES
      print *
    endif
    print *,'on NEC SX, make sure "loopcnt=" parameter'
! use fused loops on NEC SX
    print *,'in Makefile is greater than max vector length = ',NGLOB_REGIONS(IREGION_CRUST_MANTLE) * NDIM
    print *

    print *,'approximate static memory needed by the solver:'
    print *,'----------------------------------------------'
    print *
    print *,'(lower bound, usually the real amount used is 5% to 10% higher)'
    print *
    print *,'(you can get a more precise estimate of the size used per MPI process'
    print *,' by typing "size -d bin/xspecfem3D"'
    print *,' after compiling the code with the DATA/Par_file you plan to use)'
    print *
    print *,'size of static arrays per slice = ',static_memory_size/1.d6,' MB'
    print *,'                                = ',static_memory_size/1048576.d0,' MiB'
    print *,'                                = ',static_memory_size/1.d9,' GB'
    print *,'                                = ',static_memory_size/1073741824.d0,' GiB'
    print *

! note: using less memory becomes an issue only if the strong scaling of the code is poor.
!          Some users will run simulations with an executable using far less than 80% RAM per core
!          if they prefer having a faster computational time (and use a higher number of cores).

    print *,'   (should be below 80% or 90% of the memory installed per core)'
    print *,'   (if significantly more, the job will not run by lack of memory)'
    print *,'   (note that if significantly less, you waste a significant amount'
    print *,'    of memory per processor core)'
    print *,'   (but that can be perfectly acceptable if you can afford it and'
    print *,'    want faster results by using more cores)'
    print *
    if (static_memory_size*dble(NPROCTOT)/1.d6 < 10000.d0) then
      print *,'size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d6,' MB'
      print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1048576.d0,' MiB'
      print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
    else
      print *,'size of static arrays for all slices = ',static_memory_size*dble(NPROCTOT)/1.d9,' GB'
    endif
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1073741824.d0,' GiB'
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1.d12,' TB'
    print *,'                                     = ',static_memory_size*dble(NPROCTOT)/1099511627776.d0,' TiB'
    print *

  endif ! of if (PRINT_INFO_TO_SCREEN)

  ! note: we will always calculate the value for NT_DUMP_ATTENUATION_optimal even if it will not be used
  !       this will avoid the need to recompile the solver if one wants to switch between simulations
  !       with UNDO_ATTENUATION set to .true. or .false.

  ! optimal dumping interval calculation can only be done when SIMULATION_TYPE == 3 in the Par_file,
  ! thus set it to that value here in this serial code even if it has a different value in the Par_file
  saved_SIMULATION_TYPE = SIMULATION_TYPE
  SIMULATION_TYPE = 3

  ! evaluate the amount of static memory needed by the solver, but imposing that SIMULATION_TYPE = 3
  ! because that is by far the most expensive setup for runs in terms of memory usage, thus that is
  ! the type of run for which we need to make sure that everything fits in memory
  call memory_eval(doubling_index,this_region_has_a_doubling, &
                   ner,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                   ratio_sampling_array,NPROCTOT,NSPEC_REGIONS,NGLOB_REGIONS, &
                   NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                   NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
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
                   NSPEC2D_BOTTOM,NSPEC2D_TOP,static_memory_size)

  call compute_optimized_dumping(static_memory_size,NT_DUMP_ATTENUATION_optimal,number_of_dumpings_to_do, &
                   static_memory_size_GB,size_to_store_at_each_time_step,disk_size_of_each_dumping)

  ! restore the simulation type that we have temporarily erased
  SIMULATION_TYPE = saved_SIMULATION_TYPE

  if (UNDO_ATTENUATION) then
    if (PRINT_INFO_TO_SCREEN) then
      print *,'*******************************************************************************'
      print *,'Estimating optimal disk dumping interval for UNDO_ATTENUATION:'
      print *,'*******************************************************************************'
      print *
      print *,'without undoing of attenuation you are using ',static_memory_size_GB,' GB per core'
      print *,'  i.e. ',sngl(100.d0 * static_memory_size_GB / MEMORY_INSTALLED_PER_CORE_IN_GB),'% of the installed memory'

      print *
      print *,'each time step to store in memory to undo attenuation will require storing ', &
                    size_to_store_at_each_time_step,' GB per core'
      print *
      print *,'*******************************************************************************'
      print *,'the optimal value is thus NT_DUMP_ATTENUATION = ',NT_DUMP_ATTENUATION_optimal
      print *,'*******************************************************************************'

      print *
      print *,'we will need to save a total of ',number_of_dumpings_to_do,' dumpings (restart files) to disk'

      print *
      print *,'each dumping on the disk to undo attenuation will require storing ',disk_size_of_each_dumping,' GB per core'

      print *
      print *,'each dumping on the disk will require storing ',disk_size_of_each_dumping*NPROCTOT,' GB for all cores'

      print *
      print *,'ALL dumpings on the disk will require storing ',disk_size_of_each_dumping*number_of_dumpings_to_do,' GB per core'

      print *
      print *,'*******************************************************************************'
      print *,'ALL dumpings on the disk will require storing ', &
                     disk_size_of_each_dumping*number_of_dumpings_to_do*NPROCTOT,' GB for all cores'
      print *,'  i.e. ',disk_size_of_each_dumping*number_of_dumpings_to_do*NPROCTOT/1000.d0,' TB'
      print *,'*******************************************************************************'
      print *

    endif
  endif

  ! copy number of elements and points in an include file for the solver
  open(unit=IOUT,file='OUTPUT_FILES/values_from_mesher.h',status='unknown',iostat=ier)
  if (ier /= 0 ) stop 'Error opening OUTPUT_FILES/values_from_mesher.h'

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
  if (INCLUDE_CENTRAL_CUBE) then
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
  write(IOUT,*) '! maximum number of points per region = ',NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  write(IOUT,*) '!'
! use fused loops on NEC SX
  write(IOUT,*) '! on NEC SX, make sure "loopcnt=" parameter'
  write(IOUT,*) '! in Makefile is greater than max vector length = ',NGLOB_REGIONS(IREGION_CRUST_MANTLE) * NDIM
  write(IOUT,*) '!'

  write(IOUT,*) '! total elements per slice = ',sum(NSPEC_REGIONS)
  write(IOUT,*) '! total points per slice = ',sum(NGLOB_REGIONS)
  write(IOUT,*) '!'
  write(IOUT,*) '! the time step of the solver will be DT = ',sngl(DT),' (s)'
  write(IOUT,*) '! the (approximate) minimum period resolved will be = ', &
            sngl(max(ANGULAR_WIDTH_ETA_IN_DEGREES,ANGULAR_WIDTH_XI_IN_DEGREES)/90.0 * 256.0/min(NEX_ETA,NEX_XI) * 17.0),' (s)'
  write(IOUT,*) '!'

  write(IOUT,'(1x,a,i1,a)') '! total for full ',NCHUNKS,'-chunk mesh:'
  write(IOUT,*) '! ---------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! exact total number of spectral elements in entire mesh = '
  write(IOUT,*) '! ',dble(NCHUNKS)*dble(NPROC)*dble(sum(NSPEC_REGIONS)) - subtract_central_cube_elems
  write(IOUT,*) '! approximate total number of points in entire mesh = '
  write(IOUT,*) '! ',dble(NCHUNKS)*dble(NPROC)*dble(sum(NGLOB_REGIONS)) - subtract_central_cube_points
! there are 3 DOFs in solid regions, but only 1 in fluid outer core
  write(IOUT,*) '! approximate total number of degrees of freedom in entire mesh = '
  write(IOUT,*) '! ',dble(NCHUNKS)*dble(NPROC)*(3.d0*(dble(sum(NGLOB_REGIONS))) &
    - 2.d0*dble(NGLOB_REGIONS(IREGION_OUTER_CORE))) &
    - 3.d0*subtract_central_cube_points
  write(IOUT,*) '!'

! convert width to radians
  ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES * DEGREES_TO_RADIANS
  ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES * DEGREES_TO_RADIANS

! display location of chunk if regional run
  if (NCHUNKS /= 6) then

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
      call geocentric_2_geographic_dble(theta_corner,colat_corner)

      if (phi_corner > PI) phi_corner=phi_corner-TWO_PI

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

  ! mesh averages
  if (NCHUNKS /= 6) then
    ! regional mesh, variable chunk sizes
    num_elem_gc = int( 90.d0 / ANGULAR_WIDTH_XI_IN_DEGREES * 4 * NEX_XI )
    num_gll_gc = int( 90.d0 / ANGULAR_WIDTH_XI_IN_DEGREES * 4 * NEX_XI *(NGLLX-1) )
    avg_dist_deg = max( ANGULAR_WIDTH_XI_IN_DEGREES/NEX_XI,ANGULAR_WIDTH_ETA_IN_DEGREES/NEX_ETA ) / dble(NGLLX-1)
    avg_dist_km = max( ANGULAR_WIDTH_XI_RAD/NEX_XI,ANGULAR_WIDTH_ETA_RAD/NEX_ETA ) * R_EARTH_KM / dble(NGLLX-1)
    avg_element_size = max( ANGULAR_WIDTH_XI_RAD/NEX_XI,ANGULAR_WIDTH_ETA_RAD/NEX_ETA ) * R_EARTH_KM
  else
    ! global mesh, chunks of 90 degrees
    num_elem_gc = 4 * NEX_XI
    num_gll_gc = 4*NEX_XI*(NGLLX-1)
    avg_dist_deg = 360.d0 / dble(4) / dble(NEX_XI*(NGLLX-1))
    avg_dist_km = TWO_PI / dble(4) * R_EARTH_KM / dble(NEX_XI*(NGLLX-1))
    avg_element_size = TWO_PI / dble(4) * R_EARTH_KM / dble(NEX_XI)
  endif

  write(IOUT,*) '! resolution of the mesh at the surface:'
  write(IOUT,*) '! -------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! spectral elements along a great circle = ',num_elem_gc
  write(IOUT,*) '! GLL points along a great circle = ',num_gll_gc
  write(IOUT,*) '! average distance between points in degrees = ',real(avg_dist_deg)
  write(IOUT,*) '! average distance between points in km = ',real(avg_dist_km)
  write(IOUT,*) '! average size of a spectral element in km = ',real(avg_element_size)
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
  if (static_memory_size*dble(NPROCTOT)/1.d6 < 10000.d0) then
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
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE = ',NSPEC_REGIONS(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE = ',NSPEC_REGIONS(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE = ',NSPEC_REGIONS(IREGION_INNER_CORE)
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE = ',NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE = ',NGLOB_REGIONS(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_INNER_CORE = ',NGLOB_REGIONS(IREGION_INNER_CORE)
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_IC = ',NSPECMAX_ANISO_IC
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPECMAX_ISO_MANTLE = ',NSPECMAX_ISO_MANTLE
  write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',NSPECMAX_TISO_MANTLE
  write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_MANTLE = ',NSPECMAX_ANISO_MANTLE
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ATTENUATION = ',NSPEC_CRUST_MANTLE_ATTENUATION
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

  ! unused... (dynamic allocation used)
  !if (ANISOTROPIC_KL) then
  !  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL = ',NSPEC_CRUST_MANTLE_ADJOINT
  !else
  !  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL = ',1
  !endif

  ! unused... (dynamic allocation used)
  !if (APPROXIMATE_HESS_KL) then
  !  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_HESS = ',NSPEC_CRUST_MANTLE_ADJOINT
  !else
  !  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_HESS = ',1
  !endif

  ! unused... (dynamic allocation used)
  !if (NOISE_TOMOGRAPHY > 0) then
  !  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_NOISE = ',NSPEC_CRUST_MANTLE_ADJOINT
  !else
  !  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_ADJOINT_NOISE = ',1
  !endif

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

! this to allow for code elimination by the compiler in the solver for performance

  if (TRANSVERSE_ISOTROPY) then
    write(IOUT,*) 'logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.'
  endif
  write(IOUT,*)

  if (ANISOTROPIC_3D_MANTLE) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_3D_MANTLE_VAL = .false.'
  endif
  write(IOUT,*)

  if (ANISOTROPIC_INNER_CORE) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.'
  endif
  write(IOUT,*)

  if (ATTENUATION) then
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .false.'
  endif
  write(IOUT,*)

  if (ATTENUATION_3D) then
    write(IOUT,*) 'logical, parameter :: ATTENUATION_3D_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_3D_VAL = .false.'
  endif
  write(IOUT,*)

  if (ELLIPTICITY) then
    write(IOUT,*) 'logical, parameter :: ELLIPTICITY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ELLIPTICITY_VAL = .false.'
  endif
  write(IOUT,*)

  if (GRAVITY) then
    write(IOUT,*) 'logical, parameter :: GRAVITY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: GRAVITY_VAL = .false.'
  endif
  write(IOUT,*)

  if (OCEANS) then
    write(IOUT,*) 'logical, parameter :: OCEANS_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: OCEANS_VAL = .false.'
  endif
  write(IOUT,*)

  if (TOPOGRAPHY .or. OCEANS) then
    write(IOUT,*) 'integer, parameter :: NX_BATHY_VAL = ',NX_BATHY
    write(IOUT,*) 'integer, parameter :: NY_BATHY_VAL = ',NY_BATHY
  else
    write(IOUT,*) 'integer, parameter :: NX_BATHY_VAL = 1'
    write(IOUT,*) 'integer, parameter :: NY_BATHY_VAL = 1'
  endif
  write(IOUT,*)

  if (ROTATION) then
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .false.'
  endif
  if (EXACT_MASS_MATRIX_FOR_ROTATION) then
    write(IOUT,*) 'logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: EXACT_MASS_MATRIX_FOR_ROTATION_VAL = .false.'
  endif
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_ROTATION = ',NSPEC_OUTER_CORE_ROTATION
  write(IOUT,*)

  if (PARTIAL_PHYS_DISPERSION_ONLY) then
    write(IOUT,*) 'logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: PARTIAL_PHYS_DISPERSION_ONLY_VAL = .false.'
  endif
  write(IOUT,*)

  write(IOUT,*) 'integer, parameter :: NPROC_XI_VAL = ',NPROC_XI
  write(IOUT,*) 'integer, parameter :: NPROC_ETA_VAL = ',NPROC_ETA
  write(IOUT,*) 'integer, parameter :: NCHUNKS_VAL = ',NCHUNKS
  write(IOUT,*) 'integer, parameter :: NPROCTOT_VAL = ',NPROCTOT
  write(IOUT,*)

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

  ! Deville routines only implemented for NGLLX = NGLLY = NGLLZ = 5
  if (NGLLX == 5 .and. NGLLY == 5 .and. NGLLZ == 5) then
    write(IOUT,*) 'logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: USE_DEVILLE_PRODUCTS_VAL = .false.'
  endif

  if (MOVIE_VOLUME) then
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = NSPEC_CRUST_MANTLE'
    write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = NGLOB_CRUST_MANTLE'
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_3DMOVIE = 1'
    write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_3DMOVIE = 1'
  endif
  write(IOUT,*)

  if (MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE == 4) then
    write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = NSPEC_OUTER_CORE'
  else
    write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_3DMOVIE = 1'
  endif

  ! in the case of Stacey boundary conditions, add C*delta/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be fictitious / unused

  if (NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then
     NGLOB_XY_CM = NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  else
     NGLOB_XY_CM = 1
  endif
  NGLOB_XY_IC = 1

  if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    NGLOB_XY_CM = NGLOB_REGIONS(IREGION_CRUST_MANTLE)
    NGLOB_XY_IC = NGLOB_REGIONS(IREGION_INNER_CORE)
  endif
  write(IOUT,*) 'integer, parameter :: NGLOB_XY_CM = ',NGLOB_XY_CM
  write(IOUT,*) 'integer, parameter :: NGLOB_XY_IC = ',NGLOB_XY_IC
  write(IOUT,*)

  if (ATTENUATION_1D_WITH_3D_STORAGE) then
    write(IOUT,*) 'logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_1D_WITH_3D_STORAGE_VAL = .false.'
  endif
  write(IOUT,*)

  ! we use this vectorization flag for solver routines in files **.f90
#ifdef FORCE_VECTORIZATION
  write(IOUT,*) 'logical, parameter :: FORCE_VECTORIZATION_VAL = .true.'
#else
  write(IOUT,*) 'logical, parameter :: FORCE_VECTORIZATION_VAL = .false.'
#endif
  write(IOUT,*)

  ! for UNDO_ATTENUATION
  if (UNDO_ATTENUATION) then
    write(IOUT,*) 'logical, parameter :: UNDO_ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: UNDO_ATTENUATION_VAL = .false.'
  endif
  write(IOUT,*) 'integer, parameter :: NT_DUMP_ATTENUATION_VAL = ',NT_DUMP_ATTENUATION_optimal
  write(IOUT,*)

  ! mesh geometry (with format specifier to avoid writing double values on a newline)
  write(IOUT,'(1x,a,f12.6)') 'double precision, parameter :: ANGULAR_WIDTH_ETA_IN_DEGREES_VAL = ',ANGULAR_WIDTH_ETA_IN_DEGREES
  write(IOUT,'(1x,a,f12.6)') 'double precision, parameter :: ANGULAR_WIDTH_XI_IN_DEGREES_VAL = ',ANGULAR_WIDTH_XI_IN_DEGREES
  write(IOUT,'(1x,a,f12.6)') 'double precision, parameter :: CENTER_LATITUDE_IN_DEGREES_VAL = ',CENTER_LATITUDE_IN_DEGREES
  write(IOUT,'(1x,a,f12.6)') 'double precision, parameter :: CENTER_LONGITUDE_IN_DEGREES_VAL = ',CENTER_LONGITUDE_IN_DEGREES
  write(IOUT,'(1x,a,f12.6)') 'double precision, parameter :: GAMMA_ROTATION_AZIMUTH_VAL = ',GAMMA_ROTATION_AZIMUTH
  write(IOUT,*)

  close(IOUT)

  end subroutine save_header_file

!
!-------------------------------------------------------------------------------------------------
!

! compute the optimal interval at which to dump restart files to disk to undo attenuation in an exact way

! Dimitri Komatitsch and Zhinan Xie, CNRS Marseille, France, June 2013.

  subroutine compute_optimized_dumping(static_memory_size,NT_DUMP_ATTENUATION_optimal,number_of_dumpings_to_do, &
                                       static_memory_size_GB,size_to_store_at_each_time_step,disk_size_of_each_dumping)

  use shared_parameters, only: NGLOB_REGIONS,NSPEC_REGIONS,NSTEP, &
    ROTATION,ATTENUATION,GPU_MODE, &
    MEMORY_INSTALLED_PER_CORE_IN_GB,PERCENT_OF_MEM_TO_USE_PER_CORE,NOISE_TOMOGRAPHY, &
    NSPEC2D_TOP

  use constants, only: NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,CUSTOM_REAL, &
    IREGION_CRUST_MANTLE,IREGION_INNER_CORE,IREGION_OUTER_CORE

  implicit none

  double precision, intent(in) :: static_memory_size
  integer, intent(out) :: NT_DUMP_ATTENUATION_optimal,number_of_dumpings_to_do
  double precision, intent(out) :: static_memory_size_GB,size_to_store_at_each_time_step,disk_size_of_each_dumping

  double precision :: what_we_can_use_in_GB

  if (MEMORY_INSTALLED_PER_CORE_IN_GB < 0.1d0) &
       stop 'less than 100 MB per core for MEMORY_INSTALLED_PER_CORE_IN_GB does not seem realistic; exiting...'
!! DK DK the value below will probably need to be increased one day, on future machines
  if (MEMORY_INSTALLED_PER_CORE_IN_GB > 512.d0) &
       stop 'more than 512 GB per core for MEMORY_INSTALLED_PER_CORE_IN_GB does not seem realistic; exiting...'

  if (PERCENT_OF_MEM_TO_USE_PER_CORE < 50.d0) &
       stop 'less than 50% for PERCENT_OF_MEM_TO_USE_PER_CORE does not seem realistic; exiting...'
  if (PERCENT_OF_MEM_TO_USE_PER_CORE > 100.d0) &
       stop 'more than 100% for PERCENT_OF_MEM_TO_USE_PER_CORE makes no sense; exiting...'
!! DK DK will need to remove the .and. .not. GPU_MODE test here
!! DK DK if the undo_attenuation buffers are stored on the GPU instead of on the host
  if (PERCENT_OF_MEM_TO_USE_PER_CORE > 92.d0 .and. .not. GPU_MODE) &
       stop 'more than 92% for PERCENT_OF_MEM_TO_USE_PER_CORE when not using GPUs is risky; exiting...'

  what_we_can_use_in_GB = MEMORY_INSTALLED_PER_CORE_IN_GB * PERCENT_OF_MEM_TO_USE_PER_CORE / 100.d0

! convert static memory size to GB
  static_memory_size_GB = static_memory_size / 1.d9

!! DK DK June 2014: TODO  this comment is true but the statement is commented out for now
!! DK DK June 2014: TODO  because there is no GPU support for UNDO_ATTENUATION yet
!! DK DK June 2014:
! in the case of GPUs, the buffers remain on the host i.e. on the CPU, thus static_memory_size_GB could be set to zero here
! because the solver uses almost no permanent host memory, since all calculations are performed and stored on the device;
! however we prefer not to do that here because we probably have some temporary copies of all the arrays created on the host first,
! and it is not clear if they are then suppressed when the time loop of the solver starts because static memory allocation
! is used for big arrays on the host rather than dynamic, thus there is no way of freeing it dynamically.
! Thus for now we prefer not to set static_memory_size_GB to zero here.
!
! if (GPU_MODE) static_memory_size_GB = 0.d0

  if (static_memory_size_GB >= MEMORY_INSTALLED_PER_CORE_IN_GB) &
    stop 'you are using more memory than what you told us is installed!!! there is an error'

  if (static_memory_size_GB >= what_we_can_use_in_GB) &
    stop 'you are using more memory than what you allowed us to use!!! there is an error'

! compute the size to store in memory at each time step
  size_to_store_at_each_time_step = 0

! displ_crust_mantle
  size_to_store_at_each_time_step = size_to_store_at_each_time_step &
    + dble(NDIM)*NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

! displ_inner_core
  size_to_store_at_each_time_step = size_to_store_at_each_time_step &
    + dble(NDIM)*NGLOB_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

! displ_outer_core and accel_outer_core (both being scalar arrays)
  size_to_store_at_each_time_step = size_to_store_at_each_time_step &
    + 2.d0*NGLOB_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

! noise_surface_movie
  if (NOISE_TOMOGRAPHY == 3) then
    size_to_store_at_each_time_step = size_to_store_at_each_time_step &
      + dble(NDIM)*dble(NGLLX*NGLLY)*dble(NSPEC2D_TOP(IREGION_CRUST_MANTLE))*dble(CUSTOM_REAL)
  endif

! convert to GB
  size_to_store_at_each_time_step = size_to_store_at_each_time_step / 1.d9

  NT_DUMP_ATTENUATION_optimal = int((what_we_can_use_in_GB - static_memory_size_GB) / size_to_store_at_each_time_step)

! compute the size of files to dump to disk
  disk_size_of_each_dumping = 0

! displ_crust_mantle, veloc_crust_mantle, accel_crust_mantle
  disk_size_of_each_dumping = disk_size_of_each_dumping + 3.d0*dble(NDIM)*NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

! displ_inner_core, veloc_inner_core, accel_inner_core
  disk_size_of_each_dumping = disk_size_of_each_dumping + 3.d0*dble(NDIM)*NGLOB_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

! displ_outer_core, veloc_outer_core, accel_outer_core (all scalar arrays)
  disk_size_of_each_dumping = disk_size_of_each_dumping + 3.d0*NGLOB_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

! A_array_rotation,B_array_rotation
  if (ROTATION) disk_size_of_each_dumping = disk_size_of_each_dumping + &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_OUTER_CORE)*2.d0*dble(CUSTOM_REAL)

  if (ATTENUATION) then
! R_memory_crust_mantle
    disk_size_of_each_dumping = disk_size_of_each_dumping + 5.d0*dble(N_SLS)*dble(NGLLX)* &
      dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

! R_memory_inner_core
    disk_size_of_each_dumping = disk_size_of_each_dumping + 5.d0*dble(N_SLS)*dble(NGLLX)* &
      dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)
  endif

! convert to GB
  disk_size_of_each_dumping = disk_size_of_each_dumping / 1.d9

!! DK DK this formula could be made more precise; currently in some cases it can probably be off by +1 or -1; does not matter much
  number_of_dumpings_to_do = ceiling( dble(NSTEP)/dble(NT_DUMP_ATTENUATION_optimal) )

  end subroutine compute_optimized_dumping

