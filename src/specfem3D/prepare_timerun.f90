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

  ! allocate files to save surface movies
  call prepare_timerun_movie_surface()

  ! output point and element information for 3D movies
  call prepare_timerun_movie_volume()

  ! sets up time increments and rotation constants
  call prepare_timerun_constants()

  ! precomputes gravity factors
  call prepare_timerun_gravity()

  ! precomputes attenuation factors
  call prepare_timerun_attenuation()

  ! initializes arrays
  call prepare_timerun_init_wavefield()

  ! reads files back from local disk or MT tape system if restart file
  ! note: for SIMULATION_TYPE 3 simulations, the stored wavefields
  !          will be read in the time loop after the Newmark time scheme update.
  !          this makes indexing and timing easier to match with adjoint wavefields indexing.
  call read_forward_arrays_startrun()

  ! prepares Stacey boundary arrays for re-construction of wavefields
  call prepare_timerun_stacey()

  ! prepares noise simulations
  call prepare_timerun_noise()

  ! prepares GPU arrays
  call prepare_GPU()

  ! prepares VTK window visualization
  call prepare_vtk_window()

  ! optimizes array memory layout for better performance
  call prepare_optimized_arrays()

  ! output info for possible OpenMP
  call prepare_openmp()

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
        write(IMAIN,*) 'using exact mass matrix for rotation'
    else
      write(IMAIN,*) 'no rotation'
    endif

    write(IMAIN,*)
    if (ATTENUATION_VAL) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'
      if (ATTENUATION_3D_VAL) &
        write(IMAIN,*) 'using 3D attenuation model'
      if (PARTIAL_PHYS_DISPERSION_ONLY_VAL ) &
        write(IMAIN,*) 'mimicking effects on velocity only'
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
                        nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                        my_neighbours_crust_mantle)
  endif

  ! crust and mantle
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                           rmassz_crust_mantle, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                           my_neighbours_crust_mantle)

  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                           rmassx_crust_mantle, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                           my_neighbours_crust_mantle)

    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE, &
                           rmassy_crust_mantle, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                           my_neighbours_crust_mantle)
  endif

  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_CM, &
                           b_rmassx_crust_mantle, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                           my_neighbours_crust_mantle)

      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_CM, &
                           b_rmassy_crust_mantle, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle,&
                           my_neighbours_crust_mantle)
    endif
  endif



  ! outer core
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_OUTER_CORE, &
                           rmass_outer_core, &
                           num_interfaces_outer_core,max_nibool_interfaces_oc, &
                           nibool_interfaces_outer_core,ibool_interfaces_outer_core,&
                           my_neighbours_outer_core)

  ! inner core
  call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INNER_CORE, &
                           rmassz_inner_core, &
                           num_interfaces_inner_core,max_nibool_interfaces_ic, &
                           nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                           my_neighbours_inner_core)

  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                             rmassx_inner_core, &
                             num_interfaces_inner_core,max_nibool_interfaces_ic, &
                             nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                             my_neighbours_inner_core)

    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                             rmassy_inner_core, &
                             num_interfaces_inner_core,max_nibool_interfaces_ic, &
                             nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                             my_neighbours_inner_core)

    if (SIMULATION_TYPE == 3) then
      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                               b_rmassx_inner_core, &
                               num_interfaces_inner_core,max_nibool_interfaces_ic, &
                               nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                               my_neighbours_inner_core)

      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_XY_IC, &
                               b_rmassy_inner_core, &
                               num_interfaces_inner_core,max_nibool_interfaces_ic, &
                               nibool_interfaces_inner_core,ibool_interfaces_inner_core,&
                               my_neighbours_inner_core)
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

  do i = 1,NGLOB_CRUST_MANTLE
    call xyz_2_rthetaphi(xstore_crust_mantle(i),ystore_crust_mantle(i),zstore_crust_mantle(i),rval,thetaval,phival)
    rstore_crust_mantle(1,i) = rval
    rstore_crust_mantle(2,i) = thetaval
    rstore_crust_mantle(3,i) = phival
  enddo

  ! convert in the outer core
  allocate(rstore_outer_core(NDIM,NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating rstore for outer core'

  do i = 1,NGLOB_OUTER_CORE
    call xyz_2_rthetaphi(xstore_outer_core(i),ystore_outer_core(i),zstore_outer_core(i),rval,thetaval,phival)
    rstore_outer_core(1,i) = rval
    rstore_outer_core(2,i) = thetaval
    rstore_outer_core(3,i) = phival
  enddo

  ! convert in the inner core
  allocate(rstore_inner_core(NDIM,NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating rstore for inner core'

  do i = 1,NGLOB_INNER_CORE
    call xyz_2_rthetaphi(xstore_inner_core(i),ystore_inner_core(i),zstore_inner_core(i),rval,thetaval,phival)
    rstore_inner_core(1,i) = rval
    rstore_inner_core(2,i) = thetaval
    rstore_inner_core(3,i) = phival
  enddo

  ! old x/y/z array not needed anymore
  deallocate(xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle)
  deallocate(xstore_outer_core,ystore_outer_core,zstore_outer_core)
  deallocate(xstore_inner_core,ystore_inner_core,zstore_inner_core)

  end subroutine prepare_timerun_convert_coord

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_movie_surface()

  use specfem_par
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: ier

  ! checks if anything to do
  if (.not. MOVIE_SURFACE) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing movie surface"
    call flush_IMAIN()
  endif

  ! only output corners
  ! note: for noise tomography, must NOT be coarse (have to be saved on all GLL points)
  if (MOVIE_COARSE) then
    ! checks setup
    if (NGLLX /= NGLLY) &
      call exit_MPI(myrank,'MOVIE_COARSE together with MOVIE_SURFACE requires NGLLX=NGLLY')
    ! number of points
    nmovie_points = 2 * 2 * NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    NIT = NGLLX - 1
  else
    ! number of points
    nmovie_points = NGLLX * NGLLY * NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    NIT = 1
  endif

  ! checks exact number of points nmovie_points
  call movie_surface_count_points()

  ! those arrays are not necessary for noise tomography, so only allocate them in MOVIE_SURFACE case
  ! writes out movie point locations to file
  call write_movie_surface_mesh()

  ! allocates movie surface arrays for wavefield values
  allocate(store_val_ux(nmovie_points), &
           store_val_uy(nmovie_points), &
           store_val_uz(nmovie_points),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface arrays')

  ! allocates arrays for gathering wavefield values
  if (myrank == 0) then
    ! only master needs full arrays
    allocate(store_val_ux_all(nmovie_points,0:NPROCTOT_VAL-1), &
             store_val_uy_all(nmovie_points,0:NPROCTOT_VAL-1), &
             store_val_uz_all(nmovie_points,0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface all arrays')
  else
    ! slave processes only need dummy arrays
    allocate(store_val_ux_all(1,1), &
             store_val_uy_all(1,1), &
             store_val_uz_all(1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface all arrays')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Movie surface:'
    write(IMAIN,*) '  Writing to moviedata*** files in output directory'
    if (MOVIE_VOLUME_TYPE == 5) then
      write(IMAIN,*) '  movie output: displacement'
    else
      write(IMAIN,*) '  movie output: velocity'
    endif
    write(IMAIN,*) '  time steps every: ',NTSTEP_BETWEEN_FRAMES
    call flush_IMAIN()
  endif

  call synchronize_all()

  end subroutine prepare_timerun_movie_surface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_movie_volume()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: ier

  ! checks if anything to do
  if (.not. MOVIE_VOLUME) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing movie volume"
    call flush_IMAIN()
  endif

  ! checks
  ! the following has to be true for the the array dimensions of eps to match with those of rstore etc..
  ! note that epsilondev and eps_trace_over_3 don't have the same dimensions.. could cause trouble
  if (NSPEC_CRUST_MANTLE_STR_OR_ATT /= NSPEC_CRUST_MANTLE) &
    stop 'NSPEC_CRUST_MANTLE_STRAINS_ATT /= NSPEC_CRUST_MANTLE'
  if (NSPEC_CRUST_MANTLE_STRAIN_ONLY /= NSPEC_CRUST_MANTLE) &
    stop 'NSPEC_CRUST_MANTLE_STRAIN_ONLY /= NSPEC_CRUST_MANTLE'
  ! checks movie type
  if (MOVIE_VOLUME_TYPE < 1 .or. MOVIE_VOLUME_TYPE > 9) &
    stop 'MOVIE_VOLUME_TYPE has to be in range from 1 to 9'

  ! counts total number of points for movie file output
  call movie_volume_count_points()

  allocate(nu_3dmovie(3,3,npoints_3dmovie),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating nu for 3D movie')

  call write_movie_volume_mesh(nu_3dmovie,num_ibool_3dmovie,mask_3dmovie,mask_ibool, &
                                          muvstore_crust_mantle_3dmovie,npoints_3dmovie)

  if (myrank == 0) then
    write(IMAIN,*) '  Movie volume:'
    write(IMAIN,*) '  Writing to movie3D*** files on local disk databases directory'
    select case (MOVIE_VOLUME_TYPE)
    case (1)
      write(IMAIN,*) '  movie output: strains'
    case (2)
      write(IMAIN,*) '  movie output: time integral of strains'
    case (3)
      write(IMAIN,*) '  movie output: potency or integral of strain'
    case (4)
      write(IMAIN,*) '  movie output: divergence and curl'
    case (5)
      write(IMAIN,*) '  movie output: displacement'
    case (6)
      write(IMAIN,*) '  movie output: velocity'
    case (7)
      write(IMAIN,*) '  movie output: norm of displacement'
    case (8)
      write(IMAIN,*) '  movie output: norm of velocity'
    case (9)
      write(IMAIN,*) '  movie output: norm of acceleration'
    case default
      call exit_MPI(myrank, 'MOVIE_VOLUME_TYPE has to be in range from 1 to 9')
    end select
    write(IMAIN,*)
    write(IMAIN,*) '  depth(T,B):',MOVIE_TOP,MOVIE_BOTTOM
    write(IMAIN,*) '  lon(W,E)  :',MOVIE_WEST,MOVIE_EAST
    write(IMAIN,*) '  lat(S,N)  :',MOVIE_SOUTH,MOVIE_NORTH
    write(IMAIN,*) '  Starting at time step:',MOVIE_START, 'ending at:',MOVIE_STOP,'every: ',NTSTEP_BETWEEN_FRAMES
    call flush_IMAIN()
  endif

  call synchronize_all()

  end subroutine prepare_timerun_movie_volume

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

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_gravity()

! precomputes gravity factors

  use specfem_par
  implicit none

  ! local parameters
  double precision :: rspl_gravity(NR),gspl(NR),gspl2(NR)
  double precision :: radius,radius_km,g,dg
  double precision :: g_cmb_dble,g_icb_dble
  double precision :: rho,drhodr,vp,vs,Qkappa,Qmu
  integer :: int_radius,idoubling,nspl_gravity

  if (myrank == 0) then
    write(IMAIN,*) "preparing gravity arrays"
    call flush_IMAIN()
  endif

  ! store g, rho and dg/dr=dg using normalized radius in lookup table every 100 m
  ! get density and velocity from PREM model using dummy doubling flag
  ! this assumes that the gravity perturbations are small and smooth
  ! and that we can neglect the 3D model and use PREM every 100 m in all cases
  ! this is probably a rather reasonable assumption
  if (GRAVITY_VAL) then
    call make_gravity(nspl_gravity,rspl_gravity,gspl,gspl2,ONE_CRUST)
    do int_radius = 1,NRAD_GRAVITY
      radius = dble(int_radius) / (R_EARTH_KM * 10.d0)
      call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g)

      ! use PREM density profile to calculate gravity (fine for other 1D models)
      idoubling = 0
      call model_prem_iso(myrank,radius,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,.false., &
          ONE_CRUST,.false.,RICB,RCMB,RTOPDDOUBLEPRIME, &
          R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

      dg = 4.0d0*rho - 2.0d0*g/radius

      minus_gravity_table(int_radius) = - g
      minus_deriv_gravity_table(int_radius) = - dg
      density_table(int_radius) = rho
      minus_rho_g_over_kappa_fluid(int_radius) = - g / vp**2
    enddo

    ! make sure fluid array is only assigned in outer core between 1222 and 3478 km
    ! lookup table is defined every 100 m
    do int_radius = 1,NRAD_GRAVITY
      radius_km = dble(int_radius) / 10.d0
      if (radius_km > RCMB/1000.d0 - 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RCMB/1000.d0 - 3.d0)*10.d0))
      if (radius_km < RICB/1000.d0 + 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RICB/1000.d0 + 3.d0)*10.d0))
    enddo

    ! compute gravity value at CMB and ICB once and for all
    radius = RCMB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_cmb_dble)

    radius = RICB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_icb_dble)

    ! distinguish between single and double precision for reals
    minus_g_cmb = real(- g_cmb_dble, kind=CUSTOM_REAL)
    minus_g_icb = real(- g_icb_dble, kind=CUSTOM_REAL)

  else

    ! tabulate d ln(rho)/dr needed for the no gravity fluid potential
    do int_radius = 1,NRAD_GRAVITY
       radius = dble(int_radius) / (R_EARTH_KM * 10.d0)
       idoubling = 0
       call model_prem_iso(myrank,radius,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,.false., &
           ONE_CRUST,.false.,RICB,RCMB,RTOPDDOUBLEPRIME, &
           R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

       d_ln_density_dr_table(int_radius) = drhodr/rho

    enddo

  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_gravity

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_attenuation()

  ! precomputes attenuation factors

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_movie,only: muvstore_crust_mantle_3dmovie

  implicit none

  ! local parameters
  double precision, dimension(N_SLS) :: alphaval_dble, betaval_dble, gammaval_dble

  double precision :: scale_factor,scale_factor_minus_one
  real(kind=CUSTOM_REAL) :: mul
  integer :: ispec,i,j,k

  ! checks if attenuation is on and anything to do
  if (.not. ATTENUATION_VAL ) return

  ! get and store PREM attenuation model
  if (myrank == 0) then
    write(IMAIN,*) "preparing attenuation"
    call flush_IMAIN()
  endif

  ! reads in attenuation values
  ! CRUST_MANTLE ATTENUATION
  call get_attenuation_model_3D(myrank,IREGION_CRUST_MANTLE, &
                                one_minus_sum_beta_crust_mantle, &
                                factor_common_crust_mantle, &
                                factor_scale_crust_mantle,tau_sigma_dble, &
                                NSPEC_CRUST_MANTLE)
  call bcast_attenuation_model_3D(one_minus_sum_beta_crust_mantle, &
                                  factor_common_crust_mantle, &
                                  factor_scale_crust_mantle, &
                                  tau_sigma_dble, &
                                  NSPEC_CRUST_MANTLE)

  ! INNER_CORE ATTENUATION
  call get_attenuation_model_3D(myrank,IREGION_INNER_CORE, &
                                one_minus_sum_beta_inner_core, &
                                factor_common_inner_core, &
                                factor_scale_inner_core,tau_sigma_dble, &
                                NSPEC_INNER_CORE)
  call bcast_attenuation_model_3D(one_minus_sum_beta_inner_core, &
                                  factor_common_inner_core, &
                                  factor_scale_inner_core, &
                                  tau_sigma_dble, &
                                  NSPEC_INNER_CORE)

  ! if attenuation is on, shift PREM to right frequency
  ! rescale mu in PREM to average frequency for attenuation
  ! the formulas to implement the scaling can be found for instance in
  ! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
  ! anelasticity: implications for seismology and mantle composition,
  ! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
  ! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
  ! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170.
  ! Beware that in the book of Aki and Richards eq. (5.81) is given for velocities
  ! while we need an equation for "mu" and thus we have an additional factor of 2
  ! in the scaling factor below and in equation (49) of Komatitsch and Tromp, Geophys. J. Int. (2002) 149, 390-412,
  ! because "mu" is related to the square of velocity.

  ! rescale in crust and mantle

  do ispec = 1,NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            scale_factor = factor_scale_crust_mantle(i,j,k,ispec)
          else
            scale_factor = factor_scale_crust_mantle(1,1,1,ispec)
          endif

          if (ANISOTROPIC_3D_MANTLE_VAL) then
            scale_factor_minus_one = scale_factor - 1.d0

            mul = c44store_crust_mantle(i,j,k,ispec)
            c11store_crust_mantle(i,j,k,ispec) = c11store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c12store_crust_mantle(i,j,k,ispec) = c12store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c13store_crust_mantle(i,j,k,ispec) = c13store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c22store_crust_mantle(i,j,k,ispec) = c22store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c23store_crust_mantle(i,j,k,ispec) = c23store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c33store_crust_mantle(i,j,k,ispec) = c33store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c44store_crust_mantle(i,j,k,ispec) = c44store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
            c55store_crust_mantle(i,j,k,ispec) = c55store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
            c66store_crust_mantle(i,j,k,ispec) = c66store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
          else
            if (MOVIE_VOLUME .and. SIMULATION_TYPE == 3) then
              ! store the original value of \mu to compute \mu*\eps
              muvstore_crust_mantle_3dmovie(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec)
            endif

            muvstore_crust_mantle(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec) * scale_factor

            ! scales transverse isotropic values for mu_h
            if (ispec_is_tiso_crust_mantle(ispec)) then
              muhstore_crust_mantle(i,j,k,ispec) = muhstore_crust_mantle(i,j,k,ispec) * scale_factor
            endif
          endif

        enddo
      enddo
    enddo
  enddo ! enddo CRUST MANTLE

  ! rescale in inner core

  do ispec = 1,NSPEC_INNER_CORE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
            scale_factor = factor_scale_inner_core(i,j,k,ispec)
          else
            scale_factor = factor_scale_inner_core(1,1,1,ispec)
          endif

          if (ANISOTROPIC_INNER_CORE_VAL) then
            scale_factor_minus_one = scale_factor - 1.d0

            mul = muvstore_inner_core(i,j,k,ispec)
            c11store_inner_core(i,j,k,ispec) = c11store_inner_core(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c12store_inner_core(i,j,k,ispec) = c12store_inner_core(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c13store_inner_core(i,j,k,ispec) = c13store_inner_core(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c33store_inner_core(i,j,k,ispec) = c33store_inner_core(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c44store_inner_core(i,j,k,ispec) = c44store_inner_core(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
          endif

          muvstore_inner_core(i,j,k,ispec) = muvstore_inner_core(i,j,k,ispec) * scale_factor

        enddo
      enddo
    enddo
  enddo ! enddo INNER CORE

  ! precompute Runge-Kutta coefficients
  call get_attenuation_memory_values(tau_sigma_dble, deltat, alphaval_dble, betaval_dble, gammaval_dble)
  alphaval = real(alphaval_dble, kind=CUSTOM_REAL)
  betaval  = real(betaval_dble, kind=CUSTOM_REAL)
  gammaval = real(gammaval_dble, kind=CUSTOM_REAL)

  if (SIMULATION_TYPE == 3) then
   call get_attenuation_memory_values(tau_sigma_dble, b_deltat, alphaval_dble, betaval_dble, gammaval_dble)
   b_alphaval = real(alphaval_dble, kind=CUSTOM_REAL)
   b_betaval  = real(betaval_dble, kind=CUSTOM_REAL)
   b_gammaval = real(gammaval_dble, kind=CUSTOM_REAL)
  endif

  if (USE_LDDRK) then
    tau_sigma_CUSTOM_REAL(:) = real(tau_sigma_dble(:), kind=CUSTOM_REAL)
  endif

  if (UNDO_ATTENUATION) then
   b_alphaval = alphaval
   b_betaval = betaval
   b_gammaval = gammaval
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_attenuation

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_init_wavefield()

! initializes arrays

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_noise
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: ier
  real(kind=CUSTOM_REAL) :: init_value

  if (myrank == 0) then
    write(IMAIN,*) "preparing wavefields"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! put negligible initial value to avoid very slow underflow trapping
  if (FIX_UNDERFLOW_PROBLEM) then
    init_value = VERYSMALLVAL
  else
    init_value = 0._CUSTOM_REAL
  endif

  ! initialize arrays to zero
  displ_crust_mantle(:,:) = init_value
  veloc_crust_mantle(:,:) = 0._CUSTOM_REAL
  accel_crust_mantle(:,:) = 0._CUSTOM_REAL

  displ_outer_core(:) = init_value
  veloc_outer_core(:) = 0._CUSTOM_REAL
  accel_outer_core(:) = 0._CUSTOM_REAL

  displ_inner_core(:,:) = init_value
  veloc_inner_core(:,:) = 0._CUSTOM_REAL
  accel_inner_core(:,:) = 0._CUSTOM_REAL

  ! if doing benchmark runs to measure scaling of the code,
  ! set the initial field to 1 to make sure gradual underflow trapping does not slow down the code
  if (DO_BENCHMARK_RUN_ONLY .and. SET_INITIAL_FIELD_TO_1_IN_BENCH) then
    displ_crust_mantle(:,:) = 1._CUSTOM_REAL
    veloc_crust_mantle(:,:) = 1._CUSTOM_REAL
    accel_crust_mantle(:,:) = 1._CUSTOM_REAL

    displ_outer_core(:) = 1._CUSTOM_REAL
    veloc_outer_core(:) = 1._CUSTOM_REAL
    accel_outer_core(:) = 1._CUSTOM_REAL

    displ_inner_core(:,:) = 1._CUSTOM_REAL
    veloc_inner_core(:,:) = 1._CUSTOM_REAL
    accel_inner_core(:,:) = 1._CUSTOM_REAL
  endif

  ! sensitivity kernels
  if (SIMULATION_TYPE == 3) then
    rho_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    beta_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    alpha_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL

    ! noise strength kernel
    if (NOISE_TOMOGRAPHY == 3) then
      allocate( sigma_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise sigma kernel')
      sigma_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    endif

    ! approximate hessian
    if (APPROXIMATE_HESS_KL) then
      allocate( hess_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating hessian')
      hess_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    endif

    ! For anisotropic kernels (in crust_mantle only)
    if (ANISOTROPIC_KL) then
      allocate( cijkl_kl_crust_mantle(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating full cijkl kernel in crust_mantle')
      cijkl_kl_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
    endif

    rho_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL
    alpha_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL

    rho_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL
    beta_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL
    alpha_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL

    div_displ_outer_core(:,:,:,:) = 0._CUSTOM_REAL

    ! deviatoric kernel check
    if (deviatoric_outercore) then
      nspec_beta_kl_outer_core = NSPEC_OUTER_CORE_ADJOINT
    else
      nspec_beta_kl_outer_core = 1
    endif
    allocate(beta_kl_outer_core(NGLLX,NGLLY,NGLLZ,nspec_beta_kl_outer_core),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating beta outercore')
    beta_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL
  endif

  ! initialize to be on the save side for adjoint runs SIMULATION_TYPE == 2
  ! crust/mantle
  eps_trace_over_3_crust_mantle(:,:,:,:) = init_value
  epsilondev_xx_crust_mantle(:,:,:,:) = init_value
  epsilondev_yy_crust_mantle(:,:,:,:) = init_value
  epsilondev_xy_crust_mantle(:,:,:,:) = init_value
  epsilondev_xz_crust_mantle(:,:,:,:) = init_value
  epsilondev_yz_crust_mantle(:,:,:,:) = init_value

  ! backward/reconstructed strain fields
  if (SIMULATION_TYPE == 3) then
    if (UNDO_ATTENUATION) then
      ! for undo_attenuation, whenever strain is needed it will be computed locally.
      ! pointers are using the allocated arrays for adjoint strain, however values stored in those arrays will be overwritten
      ! crust/mantle
      b_epsilondev_xx_crust_mantle => epsilondev_xx_crust_mantle
      b_epsilondev_yy_crust_mantle => epsilondev_yy_crust_mantle
      b_epsilondev_xy_crust_mantle => epsilondev_xy_crust_mantle
      b_epsilondev_xz_crust_mantle => epsilondev_xz_crust_mantle
      b_epsilondev_yz_crust_mantle => epsilondev_yz_crust_mantle
      b_eps_trace_over_3_crust_mantle => eps_trace_over_3_crust_mantle
      ! inner core
      b_epsilondev_xx_inner_core => epsilondev_xx_inner_core
      b_epsilondev_yy_inner_core => epsilondev_yy_inner_core
      b_epsilondev_xy_inner_core => epsilondev_xy_inner_core
      b_epsilondev_xz_inner_core => epsilondev_xz_inner_core
      b_epsilondev_yz_inner_core => epsilondev_yz_inner_core
      b_eps_trace_over_3_inner_core => eps_trace_over_3_inner_core
    else
      ! allocates actual arrays
      allocate(b_epsilondev_xx_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_yy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_xy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_xz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_epsilondev_yz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
               b_eps_trace_over_3_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_epsilondev*** arrays for crust/mantle')

      allocate(b_epsilondev_xx_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_yy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_xy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_xz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_epsilondev_yz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT), &
               b_eps_trace_over_3_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT),stat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error allocating b_epsilondev*** arrays for inner core')
    endif
  else
    ! initializes pointers
    nullify(b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
            b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle)
    nullify(b_eps_trace_over_3_crust_mantle)
    nullify(b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
            b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core)
    nullify(b_eps_trace_over_3_inner_core)
  endif

  ! inner core
  eps_trace_over_3_inner_core(:,:,:,:) = init_value
  epsilondev_xx_inner_core(:,:,:,:) = init_value
  epsilondev_yy_inner_core(:,:,:,:) = init_value
  epsilondev_xy_inner_core(:,:,:,:) = init_value
  epsilondev_xz_inner_core(:,:,:,:) = init_value
  epsilondev_yz_inner_core(:,:,:,:) = init_value

  if (COMPUTE_AND_STORE_STRAIN) then
    if (MOVIE_VOLUME .and. (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3)) then
      Iepsilondev_xx_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_yy_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_xy_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_xz_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Iepsilondev_yz_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
      Ieps_trace_over_3_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    endif
  endif

  ! clear memory variables if attenuation
  if (ATTENUATION_VAL) then
    R_xx_crust_mantle(:,:,:,:,:) = init_value
    R_yy_crust_mantle(:,:,:,:,:) = init_value
    R_xy_crust_mantle(:,:,:,:,:) = init_value
    R_xz_crust_mantle(:,:,:,:,:) = init_value
    R_yz_crust_mantle(:,:,:,:,:) = init_value

    R_xx_inner_core(:,:,:,:,:) = init_value
    R_yy_inner_core(:,:,:,:,:) = init_value
    R_xy_inner_core(:,:,:,:,:) = init_value
    R_xz_inner_core(:,:,:,:,:) = init_value
    R_yz_inner_core(:,:,:,:,:) = init_value
  endif

  if (ROTATION_VAL) then
    A_array_rotation = 0._CUSTOM_REAL
    B_array_rotation = 0._CUSTOM_REAL
  endif

  ! initializes backward/reconstructed arrays
  if (SIMULATION_TYPE == 3) then
    ! initializes wavefields
    b_displ_crust_mantle = 0._CUSTOM_REAL
    b_veloc_crust_mantle = 0._CUSTOM_REAL
    b_accel_crust_mantle = 0._CUSTOM_REAL

    b_displ_inner_core = 0._CUSTOM_REAL
    b_veloc_inner_core = 0._CUSTOM_REAL
    b_accel_inner_core = 0._CUSTOM_REAL

    b_displ_outer_core = 0._CUSTOM_REAL
    b_veloc_outer_core = 0._CUSTOM_REAL
    b_accel_outer_core = 0._CUSTOM_REAL

    b_epsilondev_xx_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_yy_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_xy_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_xz_crust_mantle = 0._CUSTOM_REAL
    b_epsilondev_yz_crust_mantle = 0._CUSTOM_REAL
    b_eps_trace_over_3_crust_mantle = init_value

    b_epsilondev_xx_inner_core = 0._CUSTOM_REAL
    b_epsilondev_yy_inner_core = 0._CUSTOM_REAL
    b_epsilondev_xy_inner_core = 0._CUSTOM_REAL
    b_epsilondev_xz_inner_core = 0._CUSTOM_REAL
    b_epsilondev_yz_inner_core = 0._CUSTOM_REAL
    b_eps_trace_over_3_inner_core = init_value

    if (ROTATION_VAL) then
      b_A_array_rotation = 0._CUSTOM_REAL
      b_B_array_rotation = 0._CUSTOM_REAL
    endif

    if (ATTENUATION_VAL) then
      b_R_xx_crust_mantle = 0._CUSTOM_REAL
      b_R_yy_crust_mantle = 0._CUSTOM_REAL
      b_R_xy_crust_mantle = 0._CUSTOM_REAL
      b_R_xz_crust_mantle = 0._CUSTOM_REAL
      b_R_yz_crust_mantle = 0._CUSTOM_REAL

      b_R_xx_inner_core = 0._CUSTOM_REAL
      b_R_yy_inner_core = 0._CUSTOM_REAL
      b_R_xy_inner_core = 0._CUSTOM_REAL
      b_R_xz_inner_core = 0._CUSTOM_REAL
      b_R_yz_inner_core = 0._CUSTOM_REAL
    endif
  endif

  ! Runge-Kutta time scheme
  if (USE_LDDRK) then

    ! checks
    if (SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. NOISE_TOMOGRAPHY /= 0) &
        stop 'Error: LDDRK is not implemented for adjoint or noise tomography'

    ! number of stages for scheme
    NSTAGE_TIME_SCHEME = NSTAGE   ! 6 stages

    ! scheme wavefields
    allocate(displ_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array displ_crust_mantle_lddrk'
    allocate(veloc_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array veloc_crust_mantle_lddrk'
    allocate(displ_outer_core_lddrk(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array displ_outer_core_lddrk'
    allocate(veloc_outer_core_lddrk(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array veloc_outer_core_lddrk'
    allocate(displ_inner_core_lddrk(NDIM,NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array displ_inner_core_lddrk'
    allocate(veloc_inner_core_lddrk(NDIM,NGLOB_INNER_CORE),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array veloc_inner_core_lddrk'

    displ_crust_mantle_lddrk(:,:) = init_value
    veloc_crust_mantle_lddrk(:,:) = 0._CUSTOM_REAL

    displ_outer_core_lddrk(:) = init_value
    veloc_outer_core_lddrk(:) = 0._CUSTOM_REAL

    displ_inner_core_lddrk(:,:) = init_value
    veloc_inner_core_lddrk(:,:) = 0._CUSTOM_REAL

    if (SIMULATION_TYPE == 3) then
      ! scheme adjoint wavefields
      allocate(b_displ_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_displ_crust_mantle_lddrk'
      allocate(b_veloc_crust_mantle_lddrk(NDIM,NGLOB_CRUST_MANTLE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_veloc_crust_mantle_lddrk'
      allocate(b_displ_outer_core_lddrk(NGLOB_OUTER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_displ_outer_core_lddrk'
      allocate(b_veloc_outer_core_lddrk(NGLOB_OUTER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_veloc_outer_core_lddrk'
      allocate(b_displ_inner_core_lddrk(NDIM,NGLOB_INNER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_displ_inner_core_lddrk'
      allocate(b_veloc_inner_core_lddrk(NDIM,NGLOB_INNER_CORE_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_veloc_inner_core_lddrk'
      b_displ_crust_mantle_lddrk(:,:) = init_value
      b_veloc_crust_mantle_lddrk(:,:) = 0._CUSTOM_REAL
      b_displ_outer_core_lddrk(:) = init_value
      b_veloc_outer_core_lddrk(:) = 0._CUSTOM_REAL
      b_displ_inner_core_lddrk(:,:) = init_value
      b_veloc_inner_core_lddrk(:,:) = 0._CUSTOM_REAL
    endif

    ! rotation in fluid outer core
    allocate(A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array A_array_rotation_lddrk'
    allocate(B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array B_array_rotation_lddrk'
    if (ROTATION_VAL) then
      A_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
      B_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
    endif
    if (SIMULATION_TYPE == 3) then
      allocate(b_A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_A_array_rotation_lddrk'
      allocate(b_B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROT_ADJOINT),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_B_array_rotation_lddrk'
      if (ROTATION_VAL) then
        b_A_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
        b_B_array_rotation_lddrk(:,:,:,:) = 0._CUSTOM_REAL
      endif
    endif

    ! attenuation memory variables
    ! crust/mantle
    allocate(R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_ATTENUATION), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_crust_mantle_lddrk'
    ! inner core
    allocate(R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_ATTENUATION), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_inner_core_lddrk'
    if (ATTENUATION_VAL) then
      R_xx_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_yy_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_xy_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_xz_crust_mantle_lddrk(:,:,:,:,:) = init_value
      R_yz_crust_mantle_lddrk(:,:,:,:,:) = init_value

      R_xx_inner_core_lddrk(:,:,:,:,:) = init_value
      R_yy_inner_core_lddrk(:,:,:,:,:) = init_value
      R_xy_inner_core_lddrk(:,:,:,:,:) = init_value
      R_xz_inner_core_lddrk(:,:,:,:,:) = init_value
      R_yz_inner_core_lddrk(:,:,:,:,:) = init_value
    endif
    if (SIMULATION_TYPE == 3) then
      ! crust/mantle
      allocate(b_R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               b_R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_STR_AND_ATT), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_crust_mantle_lddrk'
      ! inner core
      allocate(b_R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               b_R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_INNER_CORE_STR_AND_ATT), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_inner_core_lddrk'
      if (ATTENUATION_VAL) then
        b_R_xx_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_yy_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_xy_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_xz_crust_mantle_lddrk(:,:,:,:,:) = init_value
        b_R_yz_crust_mantle_lddrk(:,:,:,:,:) = init_value

        b_R_xx_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_yy_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_xy_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_xz_inner_core_lddrk(:,:,:,:,:) = init_value
        b_R_yz_inner_core_lddrk(:,:,:,:,:) = init_value
      endif
    endif

  else

    ! default Newmark time scheme

    ! only 1 stage for Newmark time scheme
    NSTAGE_TIME_SCHEME = 1

    ! dummy arrays needed for passing as function arguments
    allocate(A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array A_array_rotation_lddrk'
    allocate(B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array B_array_rotation_lddrk'
    allocate(R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_crust_mantle_lddrk'
    allocate(R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
             stat=ier)
    if (ier /= 0) stop 'Error: not enough memory to allocate array R_memory_inner_core_lddrk'
    if (SIMULATION_TYPE == 3) then
      allocate(b_A_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_A_array_rotation_lddrk'
      allocate(b_B_array_rotation_lddrk(NGLLX,NGLLY,NGLLZ,1),stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_B_array_rotation_lddrk'
      allocate(b_R_xx_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xy_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yz_crust_mantle_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_R_memory_crust_mantle_lddrk'
      allocate(b_R_xx_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xy_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_xz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               b_R_yz_inner_core_lddrk(NGLLX,NGLLY,NGLLZ,N_SLS,1), &
               stat=ier)
      if (ier /= 0) stop 'Error: not enough memory to allocate array b_R_memory_inner_core_lddrk'
    endif

  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_init_wavefield


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_stacey()

! sets up arrays for Stacey conditions

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_outercore

  implicit none
  ! local parameters
  integer :: ier
  integer :: nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm
  integer :: nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc,nabs_zmin_oc
  integer(kind=8) :: filesize

  ! checks if anything to do
  if (.not. ABSORBING_CONDITIONS ) return

  ! sets up absorbing boundary buffer arrays
  if (myrank == 0) then
    write(IMAIN,*) "preparing absorbing boundaries"
    call flush_IMAIN()
  endif

  ! crust_mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! allocates buffers
  if (nspec2D_xmin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmin_cm = nspec2D_xmin_crust_mantle
  else
    nabs_xmin_cm = 1
  endif
  allocate(absorb_xmin_crust_mantle(NDIM,NGLLY,NGLLZ,nabs_xmin_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmin')

  if (nspec2D_xmax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmax_cm = nspec2D_xmax_crust_mantle
  else
    nabs_xmax_cm = 1
  endif
  allocate(absorb_xmax_crust_mantle(NDIM,NGLLY,NGLLZ,nabs_xmax_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmax')

  if (nspec2D_ymin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymin_cm = nspec2D_ymin_crust_mantle
  else
    nabs_ymin_cm = 1
  endif
  allocate(absorb_ymin_crust_mantle(NDIM,NGLLX,NGLLZ,nabs_ymin_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymin')

  if (nspec2D_ymax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymax_cm = nspec2D_ymax_crust_mantle
  else
    nabs_ymax_cm = 1
  endif
  allocate(absorb_ymax_crust_mantle(NDIM,NGLLX,NGLLZ,nabs_ymax_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymax')

  ! file I/O for re-construction of wavefields
  if (nspec2D_xmin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmin_crust_mantle)

    ! total file size
    filesize = reclen_xmin_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    endif
  endif
  if (nspec2D_xmax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmax_crust_mantle)

    ! total file size
    filesize = reclen_xmax_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    endif
  endif
  if (nspec2D_ymin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymin_crust_mantle)

    ! total file size
    filesize = reclen_ymin_crust_mantle
    filesize = filesize * NSTEP


    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    endif
  endif
  if (nspec2D_ymax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymax_crust_mantle)

    ! total file size
    filesize = reclen_ymax_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif


  ! outer_core
  ! create name of database
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  ! allocates buffers
  ! xmin
  if (nspec2D_xmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmin_oc = nspec2D_xmin_outer_core
  else
    nabs_xmin_oc = 1
  endif
  allocate(absorb_xmin_outer_core(NGLLY,NGLLZ,nabs_xmin_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmin')

  ! xmax
  if (nspec2D_xmax_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmax_oc = nspec2D_xmax_outer_core
  else
    nabs_xmax_oc = 1
  endif
  allocate(absorb_xmax_outer_core(NGLLY,NGLLZ,nabs_xmax_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmax')

  ! ymin
  if (nspec2D_ymin_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymin_oc = nspec2D_ymin_outer_core
  else
    nabs_ymin_oc = 1
  endif
  allocate(absorb_ymin_outer_core(NGLLX,NGLLZ,nabs_ymin_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymin')

  ! ymax
  if (nspec2D_ymax_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymax_oc = nspec2D_ymax_outer_core
  else
    nabs_ymax_oc = 1
  endif
  allocate(absorb_ymax_outer_core(NGLLX,NGLLZ,nabs_ymax_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymax')

  ! zmin
  if (nspec2D_zmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_zmin_oc = nspec2D_zmin_outer_core
  else
    nabs_zmin_oc = 1
  endif
  allocate(absorb_zmin_outer_core(NGLLX,NGLLY,nabs_zmin_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb zmin')

  ! file I/O for re-construction of wavefields
  ! xmin
  if (nspec2D_xmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmin_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmin_outer_core)

    ! total file size
    filesize = reclen_xmin_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif
  ! xmax
  if (nspec2D_xmax_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmax_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmax_outer_core)

    ! total file size
    filesize = reclen_xmax_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
   endif
  endif
  ! ymin
  if (nspec2D_ymin_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymin_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymin_outer_core)

    ! total file size
    filesize = reclen_ymin_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    endif
  endif
  ! ymanx
  if (nspec2D_ymax_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymax_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymax_outer_core)

    ! total file size
    filesize = reclen_ymax_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif
  ! zmin
  if (nspec2D_zmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_zmin = CUSTOM_REAL * (NGLLX * NGLLY * nspec2D_zmin_outer_core)

    ! total file size
    filesize = reclen_zmin
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    endif
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_noise()

  use specfem_par
  use specfem_par_crustmantle,only: NSPEC_TOP
  use specfem_par_noise

  implicit none
  ! local parameters
  integer :: ier
  double precision :: sizeval

  ! NOISE TOMOGRAPHY
  ! checks if anything to do
  if (NOISE_TOMOGRAPHY == 0) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "preparing noise arrays"
    write(IMAIN,*) "  NOISE_TOMOGRAPHY = ",NOISE_TOMOGRAPHY
    call flush_IMAIN()
  endif

  ! checks noise setup
  call check_parameters_noise()

  ! determines file i/o buffer size for surface movies
  ! (needed for better performance on clusters, otherwise i/o will become a serious bottleneck)
  ! size of a single noise movie snapshot at surface (in MB)
  sizeval = dble(CUSTOM_REAL) * dble(NDIM) * dble(NGLLX) * dble(NGLLY) * dble(NSPEC_TOP) / 1024.d0 / 1024.d0
  ! sets file i/o buffer size
  if (NOISE_TOMOGRAPHY == 3 .and. UNDO_ATTENUATION) then
    ! needs to align with attenuation buffer size, otherwise things will get very complicated
    NT_DUMP_NOISE_BUFFER = NT_DUMP_ATTENUATION
  else
    ! sets a user specified maximum size (given in MB)
    NT_DUMP_NOISE_BUFFER = int(MAXIMUM_NOISE_BUFFER_SIZE_IN_MB / sizeval)
    ! limits size
    if (NT_DUMP_NOISE_BUFFER > NSTEP) NT_DUMP_NOISE_BUFFER = NSTEP
  endif

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "  timing:"
    write(IMAIN,*) "    start time           = ",sngl(-t0)," seconds"
    write(IMAIN,*) "    time step            = ",sngl(DT)," s"
    write(IMAIN,*) "    number of time steps = ",NSTEP
    ! noise surface movie array size
    ! (holds displacement at surface for a single time step)
    write(IMAIN,*) "  arrays:"
    write(IMAIN,*) "    size of noise surface movie array = ",sngl(sizeval),"MB"
    write(IMAIN,*) "                                      = ",sngl(sizeval / 1024.d0),"GB"
    ! buffer size for file i/o
    write(IMAIN,*) "  noise buffer: "
    write(IMAIN,*) "    number of buffered time steps = ",NT_DUMP_NOISE_BUFFER
    write(IMAIN,*) "    size of noise buffer array for each slice = ",sngl(sizeval * dble(NT_DUMP_NOISE_BUFFER)),"MB"
    write(IMAIN,*) "                                              = ",sngl(sizeval * dble(NT_DUMP_NOISE_BUFFER) / 1024.d0),"GB"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  ! for noise tomography, number of surface (movie) points needed for 'surface movie';
  ! surface output must NOT be coarse (have to be saved on all GLL points)
  ! number of points
  num_noise_surface_points = NGLLX * NGLLY * NSPEC_TOP

  ! allocates noise arrays
  if (NOISE_TOMOGRAPHY == 1) then
    ! master noise source (only needed for 1. step)
    allocate(noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise source array')
  else
    ! dummy
    allocate(noise_sourcearray(1,1,1,1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise source array')
  endif
  ! initializes
  noise_sourcearray(:,:,:,:,:) = 0._CUSTOM_REAL

  ! ensemble surface noise
  allocate(normal_x_noise(num_noise_surface_points), &
           normal_y_noise(num_noise_surface_points), &
           normal_z_noise(num_noise_surface_points), &
           mask_noise(num_noise_surface_points), &
           noise_surface_movie(NDIM,NGLLX,NGLLY,NSPEC_TOP),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise arrays')

  ! initializes
  normal_x_noise(:)            = 0._CUSTOM_REAL
  normal_y_noise(:)            = 0._CUSTOM_REAL
  normal_z_noise(:)            = 0._CUSTOM_REAL
  mask_noise(:)                = 0._CUSTOM_REAL
  noise_surface_movie(:,:,:,:) = 0._CUSTOM_REAL

  ! file i/o buffer
  ! checks integer size limit: size of buffer must fit onto an 4-byte integer (<2 GB)
  if (NSPEC_TOP > 2147483646 / (CUSTOM_REAL * NGLLX * NGLLY * NDIM * NT_DUMP_NOISE_BUFFER)) then
    print *,'buffer of noise surface_movie needed exceeds integer 4-byte limit: ',dble(reclen_noise) * dble(NT_DUMP_NOISE_BUFFER)
    print *,'  ',CUSTOM_REAL, NDIM, NGLLX * NGLLY, NSPEC_TOP,NT_DUMP_NOISE_BUFFER
    print *,'bit size fortran: ',bit_size(NSPEC_TOP)
    print *,'NT_DUMP_NOISE_BUFFER: ',NT_DUMP_NOISE_BUFFER
    print *,'Please reduce size of noise buffer for file i/o ...'
    call exit_MPI(myrank,"Error NT_DUMP_NOISE_BUFFER leads to buffer length exceeding integer limit (2 GB)")
  endif

  ! allocates buffer memory
  allocate(noise_buffer(NDIM,NGLLX,NGLLY,NSPEC_TOP,NT_DUMP_NOISE_BUFFER),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise buffer array')

  ! initializes buffer and counters
  noise_buffer(:,:,:,:,:) = 0._CUSTOM_REAL
  icounter_noise_buffer = 0
  nstep_subset_noise_buffer = 0

  ! gets noise parameters and sets up arrays
  call read_parameters_noise()

  ! user output
  if(myrank == 0) then
    ! noise simulations ignore the CMTSOLUTIONS sources but employ a noise-spectrum source S_squared instead
    write(IMAIN,*) "  ignoring CMT sources"
    select case (NOISE_TOMOGRAPHY)
    case (1)
      write(IMAIN,*) "  noise source uses master record id = ",irec_master_noise
      write(IMAIN,*) "  noise master station: ",trim(network_name(irec_master_noise))//'.'//trim(station_name(irec_master_noise))
    case (2)
      write(IMAIN,*) "  noise source uses ensemble forward source"
    case (3)
      write(IMAIN,*) "  reconstructing ensemble forward wavefield"
      write(IMAIN,*) "  noise source uses ensemble adjoint sources"
    end select
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! user output of distances to master station
  if (NOISE_TOMOGRAPHY == 1) call print_master_distances_noise()

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_noise

