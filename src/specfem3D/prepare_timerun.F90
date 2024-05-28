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

  ! full gravity preparation
  if (FULL_GRAVITY) call SIEM_prepare_solver()

  ! precomputes attenuation factors
  call prepare_attenuation()

  ! precomputes iso/tiso/aniso elastic element factors
  ! (careful with the order, prepare_attenuation() should be called before this one)
  call prepare_elastic_elements()

  ! allocates & initializes arrays
  call prepare_wavefields()

  ! sets up restarting/checkpointing
  call prepare_timerun_restarting()

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

  ! free up memory
  call prepare_deallocate_unused_arrays()

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
    write(IMAIN,*) '  current time steps: ',it_begin,' to ',it_end
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
    write(IMAIN,*) 'Reference radius of the globe used is ',R_PLANET_KM,' km'
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
  if (num_interfaces_outer_core > 0) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                             rmass_outer_core, &
                             num_interfaces_outer_core,max_nibool_interfaces_oc, &
                             nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                             my_neighbors_outer_core)
  endif

  ! inner core
  if (num_interfaces_inner_core > 0) then
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
! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,rstore_crust_mantle) &
!$OMP PRIVATE(i,rval,thetaval,phival)
!$OMP DO
#endif
  do i = 1,NGLOB_CRUST_MANTLE
    ! converts x/y/z to geocentric r/theta/phi
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
! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xstore_outer_core,ystore_outer_core,zstore_outer_core,rstore_outer_core) &
!$OMP PRIVATE(i,rval,thetaval,phival)
!$OMP DO
#endif
  do i = 1,NGLOB_OUTER_CORE
    ! converts x/y/z to geocentric r/theta/phi
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
! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(xstore_inner_core,ystore_inner_core,zstore_inner_core,rstore_inner_core) &
!$OMP PRIVATE(i,rval,thetaval,phival)
!$OMP DO
#endif
  do i = 1,NGLOB_INNER_CORE
    ! converts x/y/z to geocentric r/theta/phi
    call xyz_2_rthetaphi(xstore_inner_core(i),ystore_inner_core(i),zstore_inner_core(i),rval,thetaval,phival)

    rstore_inner_core(1,i) = rval
    rstore_inner_core(2,i) = thetaval
    rstore_inner_core(3,i) = phival
  enddo
#ifdef DANIEL_TEST_OMP_RSTORE
!$OMP ENDDO
!$OMP END PARALLEL
#endif

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
  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)        ! [s]
  scale_t_inv = dsqrt(PI*GRAV*RHOAV)        ! [1/s]

  scale_displ = R_PLANET                    ! [m]
  scale_displ_inv = ONE / scale_displ       ! [1/m]

  scale_veloc = scale_displ * scale_t_inv   ! [m/s]

  ! distinguish between single and double precision for reals
  deltat = real(DT*scale_t_inv, kind=CUSTOM_REAL)
  deltatover2 = real(0.5d0*deltat, kind=CUSTOM_REAL)
  deltatsqover2 = real(0.5d0*deltat*deltat, kind=CUSTOM_REAL)

  if (SIMULATION_TYPE == 3) then
    if (UNDO_ATTENUATION) then
      ! moves forward
      b_deltat = deltat
      b_deltatover2 = deltatover2
      b_deltatsqover2 = deltatsqover2
    else
      ! reconstructed wavefield moves backward in time from last snapshot
      b_deltat = - real(DT*scale_t_inv, kind=CUSTOM_REAL)
      b_deltatover2 = real(0.5d0*b_deltat, kind=CUSTOM_REAL)
      b_deltatsqover2 = real(0.5d0*b_deltat*b_deltat, kind=CUSTOM_REAL)
    endif
  else
    ! will not be used, but initialized
    b_deltat = 0._CUSTOM_REAL
    b_deltatover2 = 0._CUSTOM_REAL
    b_deltatsqover2 = 0._CUSTOM_REAL
  endif

  ! non-dimensionalized rotation rate of the Earth times two
  two_omega_earth = 0._CUSTOM_REAL
  b_two_omega_earth = 0._CUSTOM_REAL

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
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_constants

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_restarting()

! sets up restarting/checkpointing

  use specfem_par
  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing number of runs"
    write(IMAIN,*) "  number of runs    : ",NUMBER_OF_RUNS
    write(IMAIN,*) "  number of this run: ",NUMBER_OF_THIS_RUN
    call flush_IMAIN()
  endif

  ! safety checks for run/checkpoint number
  if (NUMBER_OF_RUNS < 1 .or. NUMBER_OF_RUNS > NSTEP) &
    stop 'number of restart runs can not be less than 1 or greater than NSTEP'

  if (NUMBER_OF_THIS_RUN < 1 .or. NUMBER_OF_THIS_RUN > NUMBER_OF_RUNS) &
    stop 'incorrect run number'

  if (SIMULATION_TYPE /= 1 .and. NUMBER_OF_RUNS /= 1) &
    stop 'Only 1 run for SIMULATION_TYPE = 2/3'

  if (USE_LDDRK .and. NUMBER_OF_RUNS > 1) &
     stop 'USE_LDDRK not supported yet for restarting simulations with NUMBER_OF_RUNS > 1'

  ! define correct time steps if restart files
  ! set start/end steps for time iteration loop
  it_begin = (NUMBER_OF_THIS_RUN - 1) * (NSTEP / NUMBER_OF_RUNS) + 1
  if (NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS) then
    it_end = NUMBER_OF_THIS_RUN * (NSTEP / NUMBER_OF_RUNS)
  else
    ! Last run may be a bit larger
    it_end = NSTEP
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  time stepping     : begin/end = ",it_begin,"/",it_end
    call flush_IMAIN()
  endif

  ! reads files back from local disk or MT tape system if restart file
  !
  ! note: for SIMULATION_TYPE 3 simulations, the stored wavefields
  !          will be read in the time loop after the Newmark time scheme update.
  !          this makes indexing and timing easier to match with adjoint wavefields indexing.
  call read_forward_arrays_startrun()

  end subroutine prepare_timerun_restarting

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_simultaneous_event_execution_shift_undoatt()

! for simultaneous events, the snapshot file I/O can lead to high peak bandwidth if done for all events at the same time.
! here, we shift the execution of events by a small, estimated time shift.

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  integer :: millisec_shift
  integer(kind=8) :: local_dim
  double precision :: snapshot_size_in_GB,estimated_io_time_in_millisec
  ! current time
  integer,dimension(8) :: tval

  ! overrides filesystem bandwidth (GB/s) for I/O
  !FILESYSTEM_IO_BANDWIDTH = 0.1d0  ! 0.1 GB/s
  ! (Frontera scratch file system: https://frontera-portal.tacc.utexas.edu/user-guide/files/)
  !FILESYSTEM_IO_BANDWIDTH = 60.d0  ! 60 GB/s

  ! only for multiple simultaneous events
  if (NUMBER_OF_SIMULTANEOUS_RUNS < 2) return
  ! only in case snapshot files are used
  if (.not. ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. (SIMULATION_TYPE == 3)) ) return

  ! shifts execution of (MPI) event subgroups
  if (SHIFT_SIMULTANEOUS_RUNS .and. FILESYSTEM_IO_BANDWIDTH > 0.d0) then
    ! estimate array sizes to store for wavefield snapshots (in GB)
    snapshot_size_in_GB = 0.d0

    ! wavefield arrays
    ! displ/veloc/accel crust-mantle
    local_dim = NDIM * NGLOB_CRUST_MANTLE
    snapshot_size_in_GB = snapshot_size_in_GB + 3.d0 * local_dim * dble(CUSTOM_REAL)

    ! displ/veloc/accel inner core
    local_dim = NDIM * NGLOB_INNER_CORE
    snapshot_size_in_GB = snapshot_size_in_GB + 3.d0 * local_dim * dble(CUSTOM_REAL)

    ! displ/veloc/accel outer core
    local_dim = NGLOB_OUTER_CORE
    snapshot_size_in_GB = snapshot_size_in_GB + 3.d0 * local_dim * dble(CUSTOM_REAL)

    ! rotation
    if (ROTATION_VAL) then
      ! A_array_rotation,..
      local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
      snapshot_size_in_GB = snapshot_size_in_GB + 2.d0 * local_dim * dble(CUSTOM_REAL)
    endif

    if (ATTENUATION_VAL) then
      ! memory variables crust-mantle R_xx_crust_mantle,..
      local_dim = N_SLS * NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_ATTENUATION
      snapshot_size_in_GB = snapshot_size_in_GB + 5.d0 * local_dim * dble(CUSTOM_REAL)

      ! memory variables crust-mantle R_xx_inner_core,..
      local_dim = N_SLS * NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_ATTENUATION
      snapshot_size_in_GB = snapshot_size_in_GB + 5.d0 * local_dim * dble(CUSTOM_REAL)
    endif

    ! converts size to GB
    snapshot_size_in_GB = snapshot_size_in_GB / 1024.d0 / 1024.d0 / 1024.d0

    ! file compression
    if (ADIOS_FOR_UNDO_ATTENUATION .and. ADIOS_COMPRESSION_ALGORITHM /= 0) then
      ! assumes a compression factor of 3x
      snapshot_size_in_GB = snapshot_size_in_GB / 3.d0
    endif

    ! total snapshot size for all processes in this event group
    snapshot_size_in_GB = snapshot_size_in_GB * NPROCTOT_VAL

    ! estimates file i/o speed (adding 10 percent for overheads)
    estimated_io_time_in_millisec = snapshot_size_in_GB / FILESYSTEM_IO_BANDWIDTH * 1000.d0 * (1.d0 + 0.01d0)

    ! determines time shift (in millisec) depending on group number
    if (estimated_io_time_in_millisec > 5000) then
      ! limits shifts to 5s
      millisec_shift = int(5000.d0 * mygroup)
    else
      millisec_shift = int(estimated_io_time_in_millisec * mygroup)
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "simultaneous events execution shift: "
      if (ADIOS_FOR_UNDO_ATTENUATION .and. ADIOS_COMPRESSION_ALGORITHM /= 0) then
        write(IMAIN,*) "  estimated snapshot size for event group (ADIOS compressed) = ",sngl(snapshot_size_in_GB * 1024.d0),"MB"
        write(IMAIN,*) "                                                             = ",sngl(snapshot_size_in_GB),"GB"
      else
        write(IMAIN,*) "  estimated snapshot size for event group        = ",sngl(snapshot_size_in_GB * 1024.d0),"MB"
        write(IMAIN,*) "                                                 = ",sngl(snapshot_size_in_GB),"GB"
      endif
      write(IMAIN,*)
      write(IMAIN,*) "  selected filesystem I/O bandwidth              = ",sngl(FILESYSTEM_IO_BANDWIDTH),"GB/s"
      write(IMAIN,*) "  estimated I/O time in millisec                 = ",sngl(estimated_io_time_in_millisec),"ms"
      write(IMAIN,*)
      write(IMAIN,*) "  (MPI) event group number                       = ",mygroup
      write(IMAIN,*) "  execution shift in millisec for this event     = ",millisec_shift,"ms"
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! let this group sleep
    call sleep_for_msec(millisec_shift)
  endif

  ! synchronize processes (for this MPI group)
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    ! outputs current time on system
    call date_and_time(VALUES=tval)

    ! user output
    write(IMAIN,'(a,i2.2,a,i2.2,a,i2.2,a,i3.3,a)') &
      "   Event starts at current clock time: ",tval(5),"h ",tval(6),"min ",tval(7),"sec ", tval(8),"msec"
    write(IMAIN,*)
    call flush_IMAIN()

    ! debug
    !print *,'debug: event group',mygroup,' shifted by ',millisec_shift, &
    !        'starts at:',tval(5),"h ",tval(6),"min ",tval(7),"sec ", tval(8),"msec"
  endif

  end subroutine prepare_simultaneous_event_execution_shift_undoatt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_deallocate_unused_arrays()

  ! free up memory by deallocating arrays that are no more needed

  use specfem_par, only: FULL_GRAVITY

  use specfem_par_crustmantle
  use specfem_par_outercore
  use specfem_par_innercore

  implicit none

  ! full gravity still needs xstore,.. arrays
  if (FULL_GRAVITY) return

  ! old x/y/z array not needed anymore
  deallocate(xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle)
  deallocate(xstore_outer_core,ystore_outer_core,zstore_outer_core)
  deallocate(xstore_inner_core,ystore_inner_core,zstore_inner_core)

  end subroutine prepare_deallocate_unused_arrays
