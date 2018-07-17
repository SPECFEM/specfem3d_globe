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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_movie

  implicit none

  ! local parameters
  ! timing
  double precision :: tCPU
  double precision, external :: wtime

  ! get MPI starting time
  time_start = wtime()

  ! user output info
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! convert x/y/z into r/theta/phi spherical coordinates
  call prepare_timerun_convert_coord()

  ! sets up time increments and rotation constants
  call prepare_timerun_constants()

  ! sets up movie arrays
  call prepare_movie()

  ! precomputes gravity factors
  call prepare_gravity()

  ! precomputes attenuation factors
  call prepare_attenuation()

  ! precomputes iso/tiso/aniso elastic element factors
  call prepare_elastic_elements()

  ! allocates & initializes arrays
  call prepare_wavefields()

  ! reads files back from local disk or MT tape system if restart file
  ! note: for SIMULATION_TYPE 3 simulations, the stored wavefields
  !          will be read in the time loop after the Newmark time scheme update.
  !          this makes indexing and timing easier to match with adjoint wavefields indexing.
  call read_forward_arrays_startrun()

  ! prepares Stacey boundary arrays for re-construction of wavefields
  call prepare_stacey()

  ! prepares noise simulations
  call prepare_noise()

  ! prepares oceans
  call prepare_oceans()

  ! prepares GPU arrays
  call prepare_GPU()

  ! prepares VTK window visualization
  call prepare_vtk_window()

  ! optimizes array memory layout for better performance
  call prepare_optimized_arrays()

  ! synchronize all the processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for preparing timerun in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*) 'time loop:'
    write(IMAIN,*)
    if (USE_LDDRK) then
      write(IMAIN,*) '              scheme:         LDDRK with',NSTAGE_TIME_SCHEME,'stages'
    else
      write(IMAIN,*) '              scheme:         Newmark'
    endif
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(((NSTEP-1)*DT-t0)/60.d0),' minutes'
    write(IMAIN,*) 'start time          :',sngl(-t0),' seconds'
    write(IMAIN,*)
    call flush_IMAIN()

    ! daniel: total time estimation
    !  average time per element per time step:
    !     global simulations:
    !       isotropic models              ~ dt = 3.1857107305545455e-05
    !       transverse isotropic models   ~ dt = 3.7492335549202518e-05
    !                                            4.2252082718299598e-05 w/ attenuation (Intel Xeon @ 2.67GHz)
    !                                            1.0069202789238270e-04 w/ simulation_type == 3 (Intel Xeon @ 2.67GHz)
    !     regional simulations:
    !       transverse isotropic models   ~ dt = 4.3039939919998860e-05 w/ attenuation (Intel Xeon @ 2.67GHz)
    !                                            7.6099242919530619e-05 (Intel Xeon @ 2.27GHz)
    !
    !  total time per time step:
    !     T_total = dt * total_number_of_elements
    !
    !     (total_number_of_elements are listed in values_from_mesher.h)
    !
    !  total time using nproc processes (slices) for NSTEP time steps:
    !     T_simulation = T_total * NSTEP / nproc
  endif

  end subroutine prepare_timerun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_user_output()

  use specfem_par
  implicit none

  ! user output
  if (myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
    write(IMAIN,*)

    write(IMAIN,*)
    if (OCEANS_VAL) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) 'no oceans'
    endif

    write(IMAIN,*)
    if (ELLIPTICITY_VAL) then
      write(IMAIN,*) 'incorporating ellipticity'
    else
      write(IMAIN,*) 'no ellipticity'
    endif

    write(IMAIN,*)
    if (TOPOGRAPHY) then
      write(IMAIN,*) 'incorporating surface topography'
    else
      write(IMAIN,*) 'no surface topography'
    endif

    write(IMAIN,*)
    if (GRAVITY_VAL) then
      write(IMAIN,*) 'incorporating self-gravitation (Cowling approximation)'
    else
      write(IMAIN,*) 'no self-gravitation'
    endif

    write(IMAIN,*)
    if (ROTATION_VAL) then
      write(IMAIN,*) 'incorporating rotation'
      if (EXACT_MASS_MATRIX_FOR_ROTATION_VAL ) &
        write(IMAIN,*) '  using exact mass matrix for rotation'
    else
      write(IMAIN,*) 'no rotation'
    endif

    write(IMAIN,*)
    if (ATTENUATION_VAL) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if (ATTENUATION_3D_VAL) &
        write(IMAIN,*) '  using 3D attenuation model'
      if (PARTIAL_PHYS_DISPERSION_ONLY_VAL ) &
        write(IMAIN,*) '  mimicking effects on velocity only'
      if (UNDO_ATTENUATION ) &
        write(IMAIN,*) '  using undo_attenuation scheme'
    else
      write(IMAIN,*) 'no attenuation'
    endif

    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_user_output

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_mass_matrices()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  if (myrank == 0) then
    write(IMAIN,*) "preparing mass matrices"
    call flush_IMAIN()
  endif

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete


  ! mass matrices need to be assembled with MPI here once and for all
  call prepare_timerun_rmass_assembly()

  ! checks that all the mass matrices are positive
  ! ocean load
  if (OCEANS_VAL) then
    if (minval(rmass_ocean_load) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the oceans')
  endif

  ! checks mass matrices

  ! crust/mantle
  if (minval(rmassx_crust_mantle) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the crust_mantle rmassx')
  if (minval(rmassy_crust_mantle) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the crust_mantle rmassy')
  if (minval(rmassz_crust_mantle) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the crust_mantle rmassz')
  ! kernel simulations
  if (SIMULATION_TYPE == 3) then
    if (minval(b_rmassx_crust_mantle) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the b_crust_mantle b_rmassx')
    if (minval(b_rmassy_crust_mantle) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the b_crust_mantle b_rmassy')
    if (minval(b_rmassz_crust_mantle) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the b_crust_mantle b_rmassz')
  endif

  ! inner core
  ! checks mass matrices for rotation
  if (minval(rmassx_inner_core) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the inner core rmassx')
  if (minval(rmassy_inner_core) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the inner core rmassy')
  if (minval(rmassz_inner_core) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the inner core rmassz')
  ! kernel simulations
  if (SIMULATION_TYPE == 3) then
    if (minval(b_rmassx_inner_core) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the b_rmassx_inner_core')
    if (minval(b_rmassy_inner_core) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the b_rmassy_inner_core')
    if (minval(b_rmassz_inner_core) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the b_rmassz_inner_core')
  endif

  ! outer core
  if (minval(rmass_outer_core) <= 0._CUSTOM_REAL) &
    call exit_MPI(myrank,'negative mass matrix term for the outer core')
  if (SIMULATION_TYPE == 3) then
    if (minval(b_rmass_outer_core) <= 0._CUSTOM_REAL) &
      call exit_MPI(myrank,'negative mass matrix term for the outer core b_rmass')
  endif

  ! mass matrix inversions
  ! for efficiency, invert final mass matrix once and for all on each slice
  ! ocean load
  if (OCEANS_VAL) rmass_ocean_load = 1._CUSTOM_REAL / rmass_ocean_load

  ! mass matrices on Stacey edges
  ! crust/mantle
  if (((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
       (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL))) then
     rmassx_crust_mantle = 1._CUSTOM_REAL / rmassx_crust_mantle
     rmassy_crust_mantle = 1._CUSTOM_REAL / rmassy_crust_mantle
  endif
  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      b_rmassx_crust_mantle = 1._CUSTOM_REAL / b_rmassx_crust_mantle
      b_rmassy_crust_mantle = 1._CUSTOM_REAL / b_rmassy_crust_mantle
    endif
  endif
  rmassz_crust_mantle = 1._CUSTOM_REAL / rmassz_crust_mantle
  ! inner core
  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
     rmassx_inner_core = 1._CUSTOM_REAL / rmassx_inner_core
     rmassy_inner_core = 1._CUSTOM_REAL / rmassy_inner_core
     if (SIMULATION_TYPE == 3) then
       b_rmassx_inner_core = 1._CUSTOM_REAL / b_rmassx_inner_core
       b_rmassy_inner_core = 1._CUSTOM_REAL / b_rmassy_inner_core
     endif
  endif
  rmassz_inner_core = 1._CUSTOM_REAL / rmassz_inner_core
  ! outer core
  rmass_outer_core = 1._CUSTOM_REAL / rmass_outer_core

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_mass_matrices

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_rmass_assembly()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  ! the mass matrix needs to be assembled with MPI here once and for all

  ! ocean load
  if (OCEANS_VAL) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                        rmass_ocean_load, &
                        num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                        nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                        my_neighbors_crust_mantle)
  endif

  ! crust and mantle
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                           rmassz_crust_mantle, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                           my_neighbors_crust_mantle)

  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                             rmassx_crust_mantle, &
                             num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                             nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                             my_neighbors_crust_mantle)

    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                             rmassy_crust_mantle, &
                             num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                             nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                             my_neighbors_crust_mantle)
  endif

  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_CM, &
                               b_rmassx_crust_mantle, &
                               num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                               nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                               my_neighbors_crust_mantle)

      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_CM, &
                               b_rmassy_crust_mantle, &
                               num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                               nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                               my_neighbors_crust_mantle)
    endif
  endif



  ! outer core
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                           rmass_outer_core, &
                           num_interfaces_outer_core,max_nibool_interfaces_oc, &
                           nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                           my_neighbors_outer_core)

  ! inner core
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                           rmassz_inner_core, &
                           num_interfaces_inner_core,max_nibool_interfaces_ic, &
                           nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                           my_neighbors_inner_core)

  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                             rmassx_inner_core, &
                             num_interfaces_inner_core,max_nibool_interfaces_ic, &
                             nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                             my_neighbors_inner_core)

    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                             rmassy_inner_core, &
                             num_interfaces_inner_core,max_nibool_interfaces_ic, &
                             nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                             my_neighbors_inner_core)

    if (SIMULATION_TYPE == 3) then
      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                               b_rmassx_inner_core, &
                               num_interfaces_inner_core,max_nibool_interfaces_ic, &
                               nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                               my_neighbors_inner_core)

      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                               b_rmassy_inner_core, &
                               num_interfaces_inner_core,max_nibool_interfaces_ic, &
                               nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                               my_neighbors_inner_core)
    endif
  endif

  ! mass matrix including central cube
  if (INCLUDE_CENTRAL_CUBE) then
    ! suppress fictitious mass matrix elements in central cube
    ! because the slices do not compute all their spectral elements in the cube
    where(rmassz_inner_core(:) <= 0.0_CUSTOM_REAL) rmassz_inner_core = 1.0_CUSTOM_REAL

    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      where(rmassx_inner_core(:) <= 0.0_CUSTOM_REAL) rmassx_inner_core = 1.0_CUSTOM_REAL
      where(rmassy_inner_core(:) <= 0.0_CUSTOM_REAL) rmassy_inner_core = 1.0_CUSTOM_REAL
      if (SIMULATION_TYPE == 3) then
        where(b_rmassx_inner_core(:) <= 0.0_CUSTOM_REAL) b_rmassx_inner_core = 1.0_CUSTOM_REAL
        where(b_rmassy_inner_core(:) <= 0.0_CUSTOM_REAL) b_rmassy_inner_core = 1.0_CUSTOM_REAL
      endif
    endif
  endif

  call synchronize_all()

  end subroutine prepare_timerun_rmass_assembly

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_convert_coord()

! converts x/y/z into r/theta/phi spherical coordinates

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  ! local parameters
  integer :: i,ier
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival

  ! change x, y, z to r, theta and phi once and for all
  !
  ! note: we merged all 3 components into a single array rstore(:,:).
  !       this makes it more suitable for performance improvements by compilers (better prefetch, more efficient memory access)

  ! convert in the crust and mantle
  allocate(rstore_crust_mantle(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating rstore for crust/mantle'

#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,rstore_crust_mantle) &
!$OMP PRIVATE(i,rval,thetaval,phival)
!$OMP DO
#endif
  do i = 1,NGLOB_CRUST_MANTLE
    call xyz_2_rthetaphi(xstore_crust_mantle(i),ystore_crust_mantle(i),zstore_crust_mantle(i),rval,thetaval,phival)
    rstore_crust_mantle(1,i) = rval
    rstore_crust_mantle(2,i) = thetaval
    rstore_crust_mantle(3,i) = phival
  enddo
#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP ENDDO
!$OMP END PARALLEL
#endif

  ! convert in the outer core
  allocate(rstore_outer_core(NDIM,NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating rstore for outer core'

#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xstore_outer_core,ystore_outer_core,zstore_outer_core,rstore_outer_core) &
!$OMP PRIVATE(i,rval,thetaval,phival)
!$OMP DO
#endif
  do i = 1,NGLOB_OUTER_CORE
    call xyz_2_rthetaphi(xstore_outer_core(i),ystore_outer_core(i),zstore_outer_core(i),rval,thetaval,phival)
    rstore_outer_core(1,i) = rval
    rstore_outer_core(2,i) = thetaval
    rstore_outer_core(3,i) = phival
  enddo
#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP ENDDO
!$OMP END PARALLEL
#endif

  ! convert in the inner core
  allocate(rstore_inner_core(NDIM,NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating rstore for inner core'

#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xstore_inner_core,ystore_inner_core,zstore_inner_core,rstore_inner_core) &
!$OMP PRIVATE(i,rval,thetaval,phival)
!$OMP DO
#endif
  do i = 1,NGLOB_INNER_CORE
    call xyz_2_rthetaphi(xstore_inner_core(i),ystore_inner_core(i),zstore_inner_core(i),rval,thetaval,phival)
    rstore_inner_core(1,i) = rval
    rstore_inner_core(2,i) = thetaval
    rstore_inner_core(3,i) = phival
  enddo
#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP ENDDO
!$OMP END PARALLEL
#endif

  ! old x/y/z array not needed anymore
  deallocate(xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle)
  deallocate(xstore_outer_core,ystore_outer_core,zstore_outer_core)
  deallocate(xstore_inner_core,ystore_inner_core,zstore_inner_core)

  end subroutine prepare_timerun_convert_coord

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_constants()

! precomputes constants for time integration

  use specfem_par
  implicit none

  if (myrank == 0) then
    write(IMAIN,*) "preparing constants"
    call flush_IMAIN()
  endif

  ! define constants for the time integration
  ! scaling to make displacement in meters and velocity in meters per second
  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)
  scale_t_inv = dsqrt(PI*GRAV*RHOAV)

  scale_displ = R_EARTH
  scale_displ_inv = ONE / scale_displ

  scale_veloc = scale_displ * scale_t_inv

  ! distinguish between single and double precision for reals
  deltat = real(DT*scale_t_inv, kind=CUSTOM_REAL)
  deltatover2 = 0.5d0*deltat
  deltatsqover2 = 0.5d0*deltat*deltat

  if (SIMULATION_TYPE == 3) then
    if (UNDO_ATTENUATION) then
      ! moves forward
      b_deltat = deltat
      b_deltatover2 = deltatover2
      b_deltatsqover2 = deltatsqover2
    else
      ! reconstructed wavefield moves backward in time from last snapshot
      b_deltat = - real(DT*scale_t_inv, kind=CUSTOM_REAL)
      b_deltatover2 = 0.5d0*b_deltat
      b_deltatsqover2 = 0.5d0*b_deltat*b_deltat
    endif
  else
    ! will not be used, but initialized
    b_deltat = 0._CUSTOM_REAL
    b_deltatover2 = 0._CUSTOM_REAL
    b_deltatsqover2 = 0._CUSTOM_REAL
  endif

  ! non-dimensionalized rotation rate of the Earth times two
  if (ROTATION_VAL) then
    ! distinguish between single and double precision for reals
    if (SIMULATION_TYPE == 1) then
      ! spinning forward
      two_omega_earth = real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv), kind=CUSTOM_REAL)
    else
      ! adjoint wavefield (time-reversed) spins backward
      two_omega_earth = - real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv), kind=CUSTOM_REAL)
    endif

    if (SIMULATION_TYPE == 3) then
      ! reconstructed wavefield together with +/- b_deltat will spin backward/forward
      b_two_omega_earth = real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv), kind=CUSTOM_REAL)
    endif
  else
    ! will still be used (e.g. in GPU calculations), so initializes to zero
    two_omega_earth = 0._CUSTOM_REAL
    if (SIMULATION_TYPE == 3) b_two_omega_earth = 0._CUSTOM_REAL
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_constants

