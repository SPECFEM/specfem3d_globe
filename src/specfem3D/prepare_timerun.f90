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
  call prepare_timerun_GPU()

  ! prepares VTK window visualization
  call prepare_vtk_window()

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
      if (EXACT_MASS_MATRIX_FOR_ROTATION ) &
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
       (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION))) then
     rmassx_crust_mantle = 1._CUSTOM_REAL / rmassx_crust_mantle
     rmassy_crust_mantle = 1._CUSTOM_REAL / rmassy_crust_mantle
  endif
  if (SIMULATION_TYPE == 3) then
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
      b_rmassx_crust_mantle = 1._CUSTOM_REAL / b_rmassx_crust_mantle
      b_rmassy_crust_mantle = 1._CUSTOM_REAL / b_rmassy_crust_mantle
    endif
  endif
  rmassz_crust_mantle = 1._CUSTOM_REAL / rmassz_crust_mantle
  ! inner core
  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
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
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION)) then
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
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
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

  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
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

    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
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
  integer :: i
  real(kind=CUSTOM_REAL) :: rval,thetaval,phival

  ! change x, y, z to r, theta and phi once and for all
  ! IMPROVE dangerous: old name kept (xstore ystore zstore) for new values

  ! convert in the crust and mantle
  do i = 1,NGLOB_CRUST_MANTLE
    call xyz_2_rthetaphi(xstore_crust_mantle(i),ystore_crust_mantle(i),zstore_crust_mantle(i),rval,thetaval,phival)
    xstore_crust_mantle(i) = rval
    ystore_crust_mantle(i) = thetaval
    zstore_crust_mantle(i) = phival
  enddo

  ! convert in the outer core
  do i = 1,NGLOB_OUTER_CORE
    call xyz_2_rthetaphi(xstore_outer_core(i),ystore_outer_core(i),zstore_outer_core(i),rval,thetaval,phival)
    xstore_outer_core(i) = rval
    ystore_outer_core(i) = thetaval
    zstore_outer_core(i) = phival
  enddo

  ! convert in the inner core
  do i = 1,NGLOB_INNER_CORE
    call xyz_2_rthetaphi(xstore_inner_core(i),ystore_inner_core(i),zstore_inner_core(i),rval,thetaval,phival)
    xstore_inner_core(i) = rval
    ystore_inner_core(i) = thetaval
    zstore_inner_core(i) = phival
  enddo

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
  ! the following has to be true for the the array dimensions of eps to match with those of xstore etc..
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

  scale_veloc = scale_displ * scale_t_inv

  ! distinguish between single and double precision for reals
  deltat = real(DT*scale_t_inv, kind=CUSTOM_REAL)
  deltatover2 = 0.5d0*deltat
  deltatsqover2 = 0.5d0*deltat*deltat

  if (SIMULATION_TYPE == 3) then
    b_deltat = - real(DT*scale_t_inv, kind=CUSTOM_REAL)
    b_deltatover2 = 0.5d0*b_deltat
    b_deltatsqover2 = 0.5d0*b_deltat*b_deltat
  endif

  ! non-dimensionalized rotation rate of the Earth times two
  if (ROTATION_VAL) then
    ! distinguish between single and double precision for reals
    if (SIMULATION_TYPE == 1) then
      two_omega_earth = real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv), kind=CUSTOM_REAL)
    else
      two_omega_earth = - real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv), kind=CUSTOM_REAL)
    endif

    if (SIMULATION_TYPE == 3) then
      b_two_omega_earth = real(2.d0 * TWO_PI / (HOURS_PER_DAY * SECONDS_PER_HOUR * scale_t_inv), kind=CUSTOM_REAL)
    endif
  else
    two_omega_earth = 0._CUSTOM_REAL
    if (SIMULATION_TYPE == 3) b_two_omega_earth = 0._CUSTOM_REAL
  endif

  if (UNDO_ATTENUATION) then
   b_deltat = deltat
   b_deltatover2 = deltatover2
   b_deltatsqover2 = deltatsqover2
   b_two_omega_earth = two_omega_earth
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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP


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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP

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
    filesize = filesize*NSTEP

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
  use specfem_par_crustmantle

  implicit none
  ! local parameters
  integer :: ier

  ! NOISE TOMOGRAPHY
  ! checks if anything to do
  if (NOISE_TOMOGRAPHY == 0) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "preparing noise arrays"
    write(IMAIN,*) "  NOISE_TOMOGRAPHY = ",NOISE_TOMOGRAPHY
    write(IMAIN,*) "  timing:"
    write(IMAIN,*) "    start time           = ",sngl(-t0)," seconds"
    write(IMAIN,*) "    time step            = ",sngl(DT)," s"
    write(IMAIN,*) "    number of time steps = ",NSTEP
    ! noise simulations ignore the CMTSOLUTIONS sources but employ a noise-spectrum source S_squared instead
    write(IMAIN,*) "  ignoring CMT sources"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  ! checks that number of spectral elements at surface is set (from read_mesh_databases_CM() routine)
  if (NSPEC_TOP /= NSPEC2D_TOP(IREGION_CRUST_MANTLE)) &
    call exit_MPI(myrank,'Error invalid number of NSPEC_TOP for noise simulation')

  ! for noise tomography, number of surface (movie) points needed for 'surface movie';
  ! surface output must NOT be coarse (have to be saved on all GLL points)
  ! number of points
  num_noise_surface_points = NGLLX * NGLLY * NSPEC_TOP

  ! allocates noise arrays
  allocate(noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP), &
           normal_x_noise(num_noise_surface_points), &
           normal_y_noise(num_noise_surface_points), &
           normal_z_noise(num_noise_surface_points), &
           mask_noise(num_noise_surface_points), &
           noise_surface_movie(NDIM,NGLLX,NGLLY,NSPEC_TOP),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating noise arrays')

  ! initializes
  noise_sourcearray(:,:,:,:,:) = 0._CUSTOM_REAL
  normal_x_noise(:)            = 0._CUSTOM_REAL
  normal_y_noise(:)            = 0._CUSTOM_REAL
  normal_z_noise(:)            = 0._CUSTOM_REAL
  mask_noise(:)                = 0._CUSTOM_REAL
  noise_surface_movie(:,:,:,:) = 0._CUSTOM_REAL

  ! gets noise parameters
  call read_parameters_noise()

  ! checks noise setup
  call check_parameters_noise()

  ! user output
  if(myrank == 0) then
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

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_timerun_noise

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_GPU()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: ier
  integer :: i,j,k,ispec,ispec2D,ipoin,iglob
  real :: free_mb,used_mb,total_mb
  ! dummy custom_real variables to convert from double precision
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable:: cr_wgll_cube
  real(kind=CUSTOM_REAL),dimension(:),allocatable:: &
    cr_d_ln_density_dr_table,cr_minus_rho_g_over_kappa_fluid, &
    cr_minus_gravity_table,cr_minus_deriv_gravity_table, &
    cr_density_table
  logical :: USE_3D_ATTENUATION_ARRAYS
  real(kind=CUSTOM_REAL) :: dummy

  ! checks if anything to do
  if (.not. GPU_MODE) return

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) "preparing fields and constants on GPU devices"
    call flush_IMAIN()
  endif
  call synchronize_all()

  if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
    USE_3D_ATTENUATION_ARRAYS = .true.
  else
    USE_3D_ATTENUATION_ARRAYS = .false.
  endif

  ! evaluates memory required
  call memory_eval_gpu()

  ! prepares general fields on GPU
  call prepare_constants_device(Mesh_pointer,myrank,NGLLX, &
                                hprime_xx,hprimewgll_xx, &
                                wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
                                NSOURCES, nsources_local, &
                                sourcearrays, &
                                islice_selected_source,ispec_selected_source, &
                                nrec, nrec_local, nadj_rec_local, &
                                number_receiver_global, &
                                islice_selected_rec,ispec_selected_rec, &
                                NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                                NSPEC_CRUST_MANTLE_STRAIN_ONLY, &
                                NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                                NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                NSPEC_INNER_CORE_STRAIN_ONLY, &
                                SIMULATION_TYPE,NOISE_TOMOGRAPHY, &
                                SAVE_FORWARD,ABSORBING_CONDITIONS, &
                                OCEANS_VAL,GRAVITY_VAL, &
                                ROTATION_VAL,EXACT_MASS_MATRIX_FOR_ROTATION, &
                                ATTENUATION_VAL,UNDO_ATTENUATION, &
                                PARTIAL_PHYS_DISPERSION_ONLY,USE_3D_ATTENUATION_ARRAYS, &
                                COMPUTE_AND_STORE_STRAIN, &
                                ANISOTROPIC_3D_MANTLE_VAL,ANISOTROPIC_INNER_CORE_VAL, &
                                SAVE_BOUNDARY_MESH, &
                                USE_MESH_COLORING_GPU, &
                                ANISOTROPIC_KL,APPROXIMATE_HESS_KL, &
                                deltat,b_deltat, &
                                GPU_ASYNC_COPY)
  call synchronize_all()

  ! prepares rotation arrays
  if (ROTATION_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading rotation arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (SIMULATION_TYPE == 3) then
      call prepare_fields_rotation_device(Mesh_pointer, &
                                          two_omega_earth, &
                                          A_array_rotation,B_array_rotation, &
                                          b_two_omega_earth, &
                                          b_A_array_rotation,b_B_array_rotation, &
                                          NSPEC_OUTER_CORE_ROTATION)
    else
      call prepare_fields_rotation_device(Mesh_pointer, &
                                          two_omega_earth, &
                                          A_array_rotation,B_array_rotation, &
                                          dummy, &
                                          dummy,dummy, &
                                          NSPEC_OUTER_CORE_ROTATION)
    endif
  endif
  call synchronize_all()

  ! prepares arrays related to gravity
  ! note: GPU will use only single-precision (or double precision) for all calculations
  !          we convert to wgll_cube to custom real (by default single-precision),
  !          using implicit conversion
  if (myrank == 0) then
    write(IMAIN,*) "  loading non-gravity/gravity arrays"
    call flush_IMAIN()
  endif
  call synchronize_all()

  allocate(cr_d_ln_density_dr_table(NRAD_GRAVITY), &
           cr_minus_rho_g_over_kappa_fluid(NRAD_GRAVITY), &
           cr_minus_gravity_table(NRAD_GRAVITY), &
           cr_minus_deriv_gravity_table(NRAD_GRAVITY), &
           cr_density_table(NRAD_GRAVITY), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating cr_minus_rho_g_over_kappa_fluid, etc...'

  allocate(cr_wgll_cube(NGLLX,NGLLY,NGLLZ),stat=ier)
  if (ier /= 0 ) stop 'Error allocating cr_wgll_cube'

  ! d_ln_density_dr_table needed for no gravity case
  cr_d_ln_density_dr_table(:) = real(d_ln_density_dr_table(:), kind=CUSTOM_REAL)
  ! these are needed for gravity cases only
  cr_minus_rho_g_over_kappa_fluid(:) = real(minus_rho_g_over_kappa_fluid(:), kind=CUSTOM_REAL)
  cr_minus_gravity_table(:) = real(minus_gravity_table(:), kind=CUSTOM_REAL)
  cr_minus_deriv_gravity_table(:) = real(minus_deriv_gravity_table(:), kind=CUSTOM_REAL)
  cr_density_table(:) = real(density_table(:), kind=CUSTOM_REAL)
  cr_wgll_cube(:,:,:) = real(wgll_cube(:,:,:), kind=CUSTOM_REAL)

  ! prepares on GPU
  call prepare_fields_gravity_device(Mesh_pointer, &
                                     cr_d_ln_density_dr_table, &
                                     cr_minus_rho_g_over_kappa_fluid, &
                                     cr_minus_gravity_table, &
                                     cr_minus_deriv_gravity_table, &
                                     cr_density_table, &
                                     cr_wgll_cube, &
                                     NRAD_GRAVITY, &
                                     minus_g_icb,minus_g_cmb, &
                                     RHO_BOTTOM_OC,RHO_TOP_OC)

  deallocate(cr_d_ln_density_dr_table,cr_minus_rho_g_over_kappa_fluid, &
            cr_minus_gravity_table,cr_minus_deriv_gravity_table, &
            cr_density_table)
  deallocate(cr_wgll_cube)
  call synchronize_all()

  ! prepares attenuation arrays
  if (ATTENUATION_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading attenuation"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (SIMULATION_TYPE == 3) then
      call prepare_fields_attenuat_device(Mesh_pointer, &
                                          R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                          R_xz_crust_mantle,R_yz_crust_mantle, &
                                          b_R_xx_crust_mantle,b_R_yy_crust_mantle,b_R_xy_crust_mantle, &
                                          b_R_xz_crust_mantle,b_R_yz_crust_mantle, &
                                          factor_common_crust_mantle, &
                                          one_minus_sum_beta_crust_mantle, &
                                          R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                          R_xz_inner_core,R_yz_inner_core, &
                                          b_R_xx_inner_core,b_R_yy_inner_core,b_R_xy_inner_core, &
                                          b_R_xz_inner_core,b_R_yz_inner_core, &
                                          factor_common_inner_core, &
                                          one_minus_sum_beta_inner_core, &
                                          alphaval,betaval,gammaval, &
                                          b_alphaval,b_betaval,b_gammaval)
    else
      call prepare_fields_attenuat_device(Mesh_pointer, &
                                          R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                          R_xz_crust_mantle,R_yz_crust_mantle, &
                                          dummy,dummy,dummy, &
                                          dummy,dummy, &
                                          factor_common_crust_mantle, &
                                          one_minus_sum_beta_crust_mantle, &
                                          R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                          R_xz_inner_core,R_yz_inner_core, &
                                          dummy,dummy,dummy, &
                                          dummy,dummy, &
                                          factor_common_inner_core, &
                                          one_minus_sum_beta_inner_core, &
                                          alphaval,betaval,gammaval, &
                                          dummy,dummy,dummy)

    endif
  endif
  call synchronize_all()


  ! prepares attenuation arrays
  if (COMPUTE_AND_STORE_STRAIN) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading strain"
      call flush_IMAIN()
    endif
    call synchronize_all()

    if (SIMULATION_TYPE == 3) then
      call prepare_fields_strain_device(Mesh_pointer, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                        b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle,b_epsilondev_xy_crust_mantle, &
                                        b_epsilondev_xz_crust_mantle,b_epsilondev_yz_crust_mantle, &
                                        eps_trace_over_3_crust_mantle, &
                                        b_eps_trace_over_3_crust_mantle, &
                                        epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                        epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                                        b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core,b_epsilondev_xy_inner_core, &
                                        b_epsilondev_xz_inner_core,b_epsilondev_yz_inner_core, &
                                        eps_trace_over_3_inner_core, &
                                        b_eps_trace_over_3_inner_core)
    else
      call prepare_fields_strain_device(Mesh_pointer, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                        dummy,dummy,dummy, &
                                        dummy,dummy, &
                                        eps_trace_over_3_crust_mantle, &
                                        dummy, &
                                        epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                                        epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                                        dummy,dummy,dummy, &
                                        dummy,dummy, &
                                        eps_trace_over_3_inner_core, &
                                        dummy)

    endif
  endif
  call synchronize_all()

  ! prepares absorbing arrays
  if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading absorbing boundaries"
      call flush_IMAIN()
    endif
    call synchronize_all()

    call prepare_fields_absorb_device(Mesh_pointer, &
                                      nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
                                      nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
                                      NSPEC2DMAX_XMIN_XMAX_CM,NSPEC2DMAX_YMIN_YMAX_CM, &
                                      nimin_crust_mantle,nimax_crust_mantle, &
                                      njmin_crust_mantle,njmax_crust_mantle, &
                                      nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
                                      ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
                                      ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
                                      normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
                                      normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
                                      jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
                                      jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
                                      rho_vp_crust_mantle,rho_vs_crust_mantle, &
                                      nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                                      nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                                      nspec2D_zmin_outer_core, &
                                      NSPEC2DMAX_XMIN_XMAX_OC,NSPEC2DMAX_YMIN_YMAX_OC, &
                                      nimin_outer_core,nimax_outer_core, &
                                      njmin_outer_core,njmax_outer_core, &
                                      nkmin_xi_outer_core,nkmin_eta_outer_core, &
                                      ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
                                      ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
                                      jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
                                      jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
                                      vp_outer_core)

  endif
  call synchronize_all()

  ! prepares MPI interfaces
  if (myrank == 0) then
    write(IMAIN,*) "  loading MPI interfaces"
    call flush_IMAIN()
  endif
  call synchronize_all()

  call prepare_mpi_buffers_device(Mesh_pointer, &
                                  num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                                  nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                                  num_interfaces_inner_core,max_nibool_interfaces_ic, &
                                  nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                                  num_interfaces_outer_core,max_nibool_interfaces_oc, &
                                  nibool_interfaces_outer_core,ibool_interfaces_outer_core)
  call synchronize_all()

  ! prepares fields on GPU for noise simulations
  if (NOISE_TOMOGRAPHY > 0) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading noise arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    call prepare_fields_noise_device(Mesh_pointer,NSPEC_TOP, &
                                     NSTEP, &
                                     ibelm_top_crust_mantle, &
                                     noise_sourcearray, &
                                     normal_x_noise,normal_y_noise,normal_z_noise, &
                                     mask_noise,jacobian2D_top_crust_mantle)

  endif
  call synchronize_all()

  ! prepares oceans arrays
  if (OCEANS_VAL) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading oceans arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    ! prepares GPU arrays for coupling with oceans
    !
    ! note: handling of coupling on GPU is slightly different than in CPU routine to avoid using a mutex
    !          to update acceleration; tests so far have shown, that with a simple mutex implementation
    !          the results differ between successive runs (probably still due to some race conditions?)
    !          here we now totally avoid mutex usage and still update each global point only once

    ! counts global points on surface to oceans
    updated_dof_ocean_load(:) = .false.
    ipoin = 0
    do ispec2D = 1,NSPEC_TOP
      ispec = ibelm_top_crust_mantle(ispec2D)
      k = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! get global point number
          iglob = ibool_crust_mantle(i,j,k,ispec)
          if (.not. updated_dof_ocean_load(iglob)) then
            ipoin = ipoin + 1
            updated_dof_ocean_load(iglob) = .true.
          endif
        enddo
      enddo
    enddo

    ! allocates arrays with all global points on ocean surface
    npoin_oceans = ipoin
    allocate(ibool_ocean_load(npoin_oceans), &
             normal_ocean_load(NDIM,npoin_oceans), &
             rmass_ocean_load_selected(npoin_oceans), &
             stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating oceans arrays')

    ! fills arrays for coupling surface at oceans
    updated_dof_ocean_load(:) = .false.
    ipoin = 0
    do ispec2D = 1,NSPEC_TOP
      ispec = ibelm_top_crust_mantle(ispec2D)
      k = NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! get global point number
          iglob = ibool_crust_mantle(i,j,k,ispec)
          if (.not. updated_dof_ocean_load(iglob)) then
            ipoin = ipoin + 1
            ! fills arrays
            ibool_ocean_load(ipoin) = iglob
            rmass_ocean_load_selected(ipoin) = rmass_ocean_load(iglob)
            normal_ocean_load(:,ipoin) = normal_top_crust_mantle(:,i,j,ispec2D)
            ! masks this global point
            updated_dof_ocean_load(iglob) = .true.
          endif
        enddo
      enddo
    enddo

    ! prepares arrays on GPU
    call prepare_oceans_device(Mesh_pointer,npoin_oceans, &
                              ibool_ocean_load, &
                              rmass_ocean_load_selected, &
                              normal_ocean_load)

    ! frees memory
    deallocate(ibool_ocean_load,rmass_ocean_load_selected,normal_ocean_load)

  endif
  call synchronize_all()

  ! prepares LDDRK arrays
  if (USE_LDDRK) then
    if (myrank == 0) then
      write(IMAIN,*) "  loading LDDRK arrays"
      call flush_IMAIN()
    endif
    call synchronize_all()

    stop 'prepare_lddrk_device not implemented yet'
    !call prepare_lddrk_device(Mesh_pointer)
  endif

  ! crust/mantle region
  if (myrank == 0) then
    write(IMAIN,*) "  loading crust/mantle region"
    call flush_IMAIN()
  endif
  call synchronize_all()

  if (SIMULATION_TYPE == 3) then
    call prepare_crust_mantle_device(Mesh_pointer, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     rhostore_crust_mantle, &
                                     kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle, &
                                     eta_anisostore_crust_mantle, &
                                     rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                     b_rmassx_crust_mantle,b_rmassy_crust_mantle, &
                                     ibool_crust_mantle, &
                                     xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                                     ispec_is_tiso_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                                     c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
                                     c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
                                     c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
                                     c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
                                     c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                                     nspec_outer_crust_mantle,nspec_inner_crust_mantle, &
                                     NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                     ibelm_bottom_crust_mantle, &
                                     NCHUNKS_VAL, &
                                     num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                                     num_elem_colors_crust_mantle)
  else
    call prepare_crust_mantle_device(Mesh_pointer, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     rhostore_crust_mantle, &
                                     kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle, &
                                     eta_anisostore_crust_mantle, &
                                     rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                     dummy,dummy, &
                                     ibool_crust_mantle, &
                                     xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                                     ispec_is_tiso_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                                     c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
                                     c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
                                     c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
                                     c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
                                     c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     num_phase_ispec_crust_mantle,phase_ispec_inner_crust_mantle, &
                                     nspec_outer_crust_mantle,nspec_inner_crust_mantle, &
                                     NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                     ibelm_bottom_crust_mantle, &
                                     NCHUNKS_VAL, &
                                     num_colors_outer_crust_mantle,num_colors_inner_crust_mantle, &
                                     num_elem_colors_crust_mantle)
  endif
  call synchronize_all()

  ! outer core region
  if (myrank == 0) then
    write(IMAIN,*) "  loading outer core region"
    call flush_IMAIN()
  endif
  call synchronize_all()

  call prepare_outer_core_device(Mesh_pointer, &
                                 xix_outer_core,xiy_outer_core,xiz_outer_core, &
                                 etax_outer_core,etay_outer_core,etaz_outer_core, &
                                 gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                                 rhostore_outer_core,kappavstore_outer_core, &
                                 rmass_outer_core, &
                                 ibool_outer_core, &
                                 xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                                 num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                                 nspec_outer_outer_core,nspec_inner_outer_core, &
                                 NSPEC2D_TOP(IREGION_OUTER_CORE), &
                                 NSPEC2D_BOTTOM(IREGION_OUTER_CORE), &
                                 normal_top_outer_core, &
                                 normal_bottom_outer_core, &
                                 jacobian2D_top_outer_core, &
                                 jacobian2D_bottom_outer_core, &
                                 ibelm_top_outer_core, &
                                 ibelm_bottom_outer_core, &
                                 num_colors_outer_outer_core,num_colors_inner_outer_core, &
                                 num_elem_colors_outer_core)
  call synchronize_all()

  ! inner core region
  if (myrank == 0) then
    write(IMAIN,*) "  loading inner core region"
    call flush_IMAIN()
  endif
  call synchronize_all()

  if (SIMULATION_TYPE == 3) then
    call prepare_inner_core_device(Mesh_pointer, &
                                   xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                   etax_inner_core,etay_inner_core,etaz_inner_core, &
                                   gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                   rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                   rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                                   b_rmassx_inner_core,b_rmassy_inner_core, &
                                   ibool_inner_core, &
                                   xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                                   c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                                   c33store_inner_core,c44store_inner_core, &
                                   idoubling_inner_core, &
                                   num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                                   nspec_outer_inner_core,nspec_inner_inner_core, &
                                   NSPEC2D_TOP(IREGION_INNER_CORE), &
                                   ibelm_top_inner_core, &
                                   num_colors_outer_inner_core,num_colors_inner_inner_core, &
                                   num_elem_colors_inner_core)
  else
    call prepare_inner_core_device(Mesh_pointer, &
                                   xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                   etax_inner_core,etay_inner_core,etaz_inner_core, &
                                   gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                   rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                   rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                                   dummy,dummy, &
                                   ibool_inner_core, &
                                   xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                                   c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                                   c33store_inner_core,c44store_inner_core, &
                                   idoubling_inner_core, &
                                   num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                                   nspec_outer_inner_core,nspec_inner_inner_core, &
                                   NSPEC2D_TOP(IREGION_INNER_CORE), &
                                   ibelm_top_inner_core, &
                                   num_colors_outer_inner_core,num_colors_inner_inner_core, &
                                   num_elem_colors_inner_core)
  endif
  call synchronize_all()

  ! transfer forward and backward fields to device with initial values
  if (myrank == 0) then
    write(IMAIN,*) "  transferring initial wavefield"
    call flush_IMAIN()
  endif
  call synchronize_all()

  call transfer_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                                    Mesh_pointer)

  call transfer_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,displ_inner_core,veloc_inner_core,accel_inner_core, &
                                    Mesh_pointer)

  call transfer_fields_oc_to_device(NGLOB_OUTER_CORE,displ_outer_core,veloc_outer_core,accel_outer_core, &
                                    Mesh_pointer)

  if (SIMULATION_TYPE == 3) then
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                        Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                        Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE,b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                        Mesh_pointer)
  endif

  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if (myrank == 0) then
    ! gets memory usage for main process
    call get_free_device_memory(free_mb,used_mb,total_mb)

    ! outputs info
    if (total_mb /= 0) then
       write(IMAIN,*)
       write(IMAIN,*)"  GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
       write(IMAIN,*)"             used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
       write(IMAIN,*)"             total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
       write(IMAIN,*)
    else
       write(IMAIN,*)
       write(IMAIN,*)"  GPU usage: not available."
    endif
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  contains

    subroutine memory_eval_gpu()

    implicit none

    ! local parameters
    double precision :: memory_size
    integer,parameter :: NGLL2 = 25
    integer,parameter :: NGLL3 = 125
    integer,parameter :: NGLL3_PADDED = 128
    integer :: NSPEC_AB,NGLOB_AB
    integer :: NGLL_ATT,NSPEC_2

    memory_size = 0.d0

    ! crust/mantle + inner core
    NSPEC_AB = NSPEC_CRUST_MANTLE + NSPEC_INNER_CORE
    NGLOB_AB = NGLOB_CRUST_MANTLE + NGLOB_INNER_CORE

    ! add size of each set of arrays multiplied by the number of such arrays
    ! d_hprime_xx,d_hprimewgll_xx
    memory_size = memory_size + 2.d0 * NGLL2 * dble(CUSTOM_REAL)

    ! sources
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      ! d_sourcearrays
      memory_size = memory_size + NGLL3 * NSOURCES * NDIM * dble(CUSTOM_REAL)
    endif
    ! d_islice_selected_source,d_ispec_selected_source
    memory_size = memory_size + 2.0 * NSOURCES * dble(SIZE_INTEGER)

    ! receivers
    !d_number_receiver_global
    memory_size = memory_size + nrec_local * dble(SIZE_INTEGER)
    ! d_station_seismo_field
    memory_size = memory_size + NDIM * NGLL3 * nrec_local * dble(CUSTOM_REAL)
    ! d_station_strain_field
    if (SIMULATION_TYPE == 2) then
      memory_size = memory_size + NGLL3 * nrec_local * dble(SIZE_INTEGER)
    endif
    ! d_ispec_selected_rec
    memory_size = memory_size + nrec * dble(SIZE_INTEGER)

    ! d_adj_sourcearrays
    memory_size = memory_size + NDIM * NGLL3 * nadj_rec_local * dble(CUSTOM_REAL)

    ! rotation
    if (ROTATION_VAL) then
      ! d_A_array_rotation,..
      memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
    endif

    ! gravity
    if (GRAVITY_VAL) then
      ! d_minus_rho_g_over_kappa_fluid,..
      memory_size = memory_size + 4.d0 * NRAD_GRAVITY * dble(CUSTOM_REAL)
    else
      ! d_d_ln_density_dr_table
      memory_size = memory_size + NRAD_GRAVITY * dble(CUSTOM_REAL)
    endif

    ! attenuation
    if (ATTENUATION_VAL) then
      if (USE_3D_ATTENUATION_ARRAYS) then
        NGLL_ATT = NGLL3
      else
        NGLL_ATT = 1
      endif
      ! d_one_minus_sum_beta_crust_mantle,..
      memory_size = memory_size + NGLL_ATT * NSPEC_AB * dble(CUSTOM_REAL)

      if (.not. PARTIAL_PHYS_DISPERSION_ONLY) then
        ! d_factor_common_crust_mantle,..
        memory_size = memory_size + N_SLS * NGLL_ATT * NSPEC_AB * dble(CUSTOM_REAL)
        ! d_R_xx_crust_mantle,..
        memory_size = memory_size + 5.d0 * N_SLS * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      endif
      ! alphaval,..
      memory_size = memory_size + 3.d0 * N_SLS * dble(CUSTOM_REAL)
    endif

    ! strains
    if (COMPUTE_AND_STORE_STRAIN) then
      ! d_epsilondev_xx_crust_mantle,..
      memory_size = memory_size + 5.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      ! d_eps_trace_over_3_crust_mantle,..
      NSPEC_2 = NSPEC_CRUST_MANTLE_STRAIN_ONLY + NSPEC_INNER_CORE_STRAIN_ONLY
      memory_size = memory_size + NGLL3 * NSPEC_2 * dble(CUSTOM_REAL)
    endif

    ! absorbing boundaries
    if (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
      ! d_rho_vp_crust_mantle,..
      memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_CRUST_MANTLE * dble(CUSTOM_REAL)
      ! d_nkmin_xi_crust_mantle,..,d_nkmin_eta_crust_mantle,..
      NSPEC_2 = NSPEC2DMAX_XMIN_XMAX_CM + NSPEC2DMAX_YMIN_YMAX_CM
      memory_size = memory_size + 6.d0 * NSPEC_2 * dble(SIZE_INTEGER)
      ! d_ibelm_xmin_crust_mantle,..
      NSPEC_2 = nspec2D_xmin_crust_mantle + nspec2D_xmax_crust_mantle
      memory_size = memory_size + NSPEC_2 * dble(SIZE_INTEGER)
      ! d_normal_xmax_crust_mantle,..
      memory_size = memory_size + NDIM * NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_jacobian2D_xmax_crust_mantle,..
      memory_size = memory_size + NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_ibelm_ymin_crust_mantle
      NSPEC_2 = nspec2D_ymin_crust_mantle + nspec2D_ymax_crust_mantle
      memory_size = memory_size + NSPEC_2 * dble(SIZE_INTEGER)
      ! d_normal_ymin_crust_mantle,..
      memory_size = memory_size + NDIM * NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_jacobian2D_ymin_crust_mantle,..
      memory_size = memory_size + NGLL2 * NSPEC_2 * dble(CUSTOM_REAL)
      ! d_vp_outer_core
      memory_size = memory_size + NGLL3 * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
      ! d_nkmin_xi_outer_core,..,d_nkmin_eta_outer_core,..
      NSPEC_2 = NSPEC2DMAX_XMIN_XMAX_OC + NSPEC2DMAX_YMIN_YMAX_OC
      memory_size = memory_size + 6.d0 * NSPEC_2 * dble(SIZE_INTEGER)
    endif

    ! mpi buffers
    ! d_ibool_interfaces_crust_mantle
    memory_size = memory_size + num_interfaces_crust_mantle * max_nibool_interfaces_cm * dble(SIZE_INTEGER)
    ! d_send_accel_buffer_crust_mantle
    memory_size = memory_size + NDIM * num_interfaces_crust_mantle * max_nibool_interfaces_cm * dble(CUSTOM_REAL)
    ! d_ibool_interfaces_inner_core
    memory_size = memory_size + num_interfaces_inner_core * max_nibool_interfaces_ic * dble(SIZE_INTEGER)
    ! d_send_accel_buffer_inner_core
    memory_size = memory_size + NDIM * num_interfaces_inner_core * max_nibool_interfaces_ic * dble(CUSTOM_REAL)
    ! d_ibool_interfaces_outer_core
    memory_size = memory_size + num_interfaces_outer_core * max_nibool_interfaces_oc * dble(SIZE_INTEGER)
    ! d_send_accel_buffer_outer_core
    memory_size = memory_size + num_interfaces_outer_core * max_nibool_interfaces_oc * dble(CUSTOM_REAL)

    ! noise
    if (NOISE_TOMOGRAPHY > 0) then
      ! d_noise_surface_movie
      memory_size = memory_size + NDIM * NGLL2 * NSPEC_TOP * dble(CUSTOM_REAL)
      ! d_noise_sourcearray
      if (NOISE_TOMOGRAPHY == 1) then
        memory_size = memory_size + NDIM * NGLL3 * NSTEP * dble(CUSTOM_REAL)
      endif
      ! d_normal_x_noise,..
      if (NOISE_TOMOGRAPHY == 1) then
        memory_size = memory_size + 4.d0 * NGLL2 * NSPEC_TOP * dble(CUSTOM_REAL)
      endif
    endif

    ! oceans
    if (OCEANS_VAL) then
      ! d_ibool_ocean_load
      memory_size = memory_size + npoin_oceans * dble(SIZE_INTEGER)
      ! d_rmass_ocean_load,..
      memory_size = memory_size + 4.d0 * npoin_oceans * dble(CUSTOM_REAL)
    endif

    ! crust/mantle + inner core
    ! padded xix,..gammaz
    memory_size = memory_size + 9.d0 * NGLL3_PADDED * NSPEC_AB * dble(CUSTOM_REAL)
    ! ibool
    memory_size = memory_size + NGLL3 * NSPEC_AB * dble(SIZE_INTEGER)
    ! crust/mantle
    if (.not. ANISOTROPIC_3D_MANTLE_VAL) then
      ! padded kappav,kappah,..
      memory_size = memory_size + 5.d0 * NGLL3_PADDED * NSPEC_CRUST_MANTLE * dble(CUSTOM_REAL)
    else
      ! padded c11,..
      memory_size = memory_size + 21.d0 * NGLL3_PADDED * NSPEC_CRUST_MANTLE * dble(CUSTOM_REAL)
    endif
    ! ystore,zstore
    memory_size = memory_size + 2.d0 * NGLOB_CRUST_MANTLE * dble(CUSTOM_REAL)
    ! xstore
    if (GRAVITY_VAL) then
      memory_size = memory_size + NGLOB_CRUST_MANTLE * dble(CUSTOM_REAL)
    endif
    ! inner core
    ! padded muv
    memory_size = memory_size + NGLL3_PADDED * NSPEC_INNER_CORE * dble(CUSTOM_REAL)
    if (.not. ANISOTROPIC_INNER_CORE_VAL) then
      ! padded kappav
      memory_size = memory_size + NGLL3_PADDED * NSPEC_INNER_CORE * dble(CUSTOM_REAL)
    else
      ! padded c11,..
      memory_size = memory_size + 5.d0 * NGLL3_PADDED * NSPEC_INNER_CORE * dble(CUSTOM_REAL)
    endif
    ! xstore,ystore,zstore
    if (GRAVITY_VAL) then
      memory_size = memory_size + 3.d0 * NGLOB_INNER_CORE * dble(CUSTOM_REAL)
    endif
    ! d_phase_ispec_inner_crust_mantle
    memory_size = memory_size + 2 * num_phase_ispec_crust_mantle * dble(SIZE_INTEGER)
    ! d_displ,..
    memory_size = memory_size + 3.d0 * NDIM * NGLOB_AB * dble(CUSTOM_REAL)
    ! crust/mantle
    ! d_rmassz
    memory_size = memory_size + NGLOB_AB * dble(CUSTOM_REAL)
    ! d_rmassx,..
    if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION)) then
      memory_size = memory_size + 2.d0 * NGLOB_CRUST_MANTLE * dble(CUSTOM_REAL)
    endif
    ! inner core
    if ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION)) then
      ! d_rmassx,..
      memory_size = memory_size + 2.d0 * NGLOB_INNER_CORE * dble(CUSTOM_REAL)
    endif

    ! outer core
    ! padded d_xix_outer_core,..
    memory_size = memory_size + 9.d0 * NGLL3_PADDED * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
    ! padded d_kappav
    memory_size = memory_size + NGLL3_PADDED * NSPEC_OUTER_CORE * dble(CUSTOM_REAL)
    ! ibool
    memory_size = memory_size + NGLL3 * NSPEC_OUTER_CORE * dble(SIZE_INTEGER)
    ! d_xstore_outer_core,..
    memory_size = memory_size + 3.d0 * NGLOB_OUTER_CORE * dble(CUSTOM_REAL)
    ! d_phase_ispec_inner_outer_core
    memory_size = memory_size + 2.d0 * num_phase_ispec_outer_core * dble(SIZE_INTEGER)
    ! d_displ_outer,..
    memory_size = memory_size + 3.d0 * NGLOB_OUTER_CORE * dble(CUSTOM_REAL)
    ! d_rmass_outer_core
    memory_size = memory_size + NGLOB_OUTER_CORE * dble(CUSTOM_REAL)

    ! poor estimate for kernel simulations...
    if (SIMULATION_TYPE == 3) memory_size = 2.d0 * memory_size

    ! user output
    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '  minimum memory requested     : ',memory_size / 1024. / 1024.,'MB per process'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    end subroutine memory_eval_gpu

  end subroutine prepare_timerun_GPU

!
!-------------------------------------------------------------------------------------------------
!

! VTK visualization

  subroutine prepare_vtk_window()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  ! local parameters
  integer :: i,j,k,iglob,ispec,inum,ier
  integer :: id1,id2,id3,id4,id5,id6,id7,id8
  integer :: ispec2D,NIT_res

  ! free surface points
  integer :: free_np,free_nspec
  real, dimension(:),allocatable :: free_x,free_y,free_z
  integer, dimension(:,:),allocatable :: free_conn
  integer, dimension(:),allocatable :: free_perm
  ! gather arrays for multi-MPI simulations
  real, dimension(:),allocatable :: free_x_all,free_y_all,free_z_all
  integer, dimension(:,:),allocatable :: free_conn_all
  integer, dimension(:),allocatable :: free_conn_offset_all,free_conn_nspec_all
  integer, dimension(:),allocatable :: free_points_all,free_offset_all
  integer :: free_np_all,free_nspec_all

  ! volume points
  integer :: vol_np,vol_nspec
  real, dimension(:),allocatable :: vol_x,vol_y,vol_z
  integer, dimension(:,:),allocatable :: vol_conn
  integer, dimension(:),allocatable :: vol_perm
  ! gather arrays for multi-MPI simulations
  real, dimension(:),allocatable :: vol_x_all,vol_y_all,vol_z_all
  integer, dimension(:,:),allocatable :: vol_conn_all
  integer, dimension(:),allocatable :: vol_conn_offset_all,vol_conn_nspec_all
  integer :: vol_nspec_all,ispec_start,ispec_end
  real,dimension(1) :: dummy
  integer,dimension(1) :: dummy_i

  real(kind=CUSTOM_REAL) :: x,y,z

  !-----------------------------------------------------------------------
  ! user parameter
  logical, parameter :: VTK_USE_HIRES         = .false.
  logical, parameter :: VTK_SHOW_FREESURFACE  = .true.
  logical, parameter :: VTK_SHOW_VOLUME       = .true.
  !-----------------------------------------------------------------------

  ! checks if anything to do
  if (.not. VTK_MODE) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing VTK runtime visualization"
    call flush_IMAIN()
  endif
  call synchronize_all()

  ! to avoid compiler warnings
  !NPROC = NPROCTOT_VAL

  ! adds source
  if (myrank == 0) then
    ! user output
    write(IMAIN,*) "  VTK source sphere:"
    call prepare_vtksource(vtkdata_source_x,vtkdata_source_y,vtkdata_source_z)
  endif
  call synchronize_all()

  ! mask
  allocate(vtkmask(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating arrays'

  if (VTK_USE_HIRES) then
    NIT_res = 1
  else
    NIT_res = NGLLX - 1
  endif

  ! free surface
  if (VTK_SHOW_FREESURFACE) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  VTK free surface:"
      write(IMAIN,*) "    free surface elements    : ",NSPEC_TOP
    endif

    ! counts global free surface points
    vtkmask(:) = .false.

    ! determines number of global points on surface
    do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
      ispec = ibelm_top_crust_mantle(ispec2D)
      ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
      ! the cubed sphere, the mapping changes for different chunks,
      ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
      ! for future consideration, like in create_movie_GMT_global.f90 ...
      k = NGLLZ
      ! loop on all the points inside the element
      do j = 1,NGLLY,NIT_res
        do i = 1,NGLLX,NIT_res
          iglob = ibool_crust_mantle(i,j,k,ispec)
          vtkmask(iglob) = .true.
        enddo
      enddo
    enddo

    ! loads free surface into data
    free_np = count(vtkmask(:))

    ! user output
    if (myrank == 0 ) write(IMAIN,*) "    loading surface points: ",free_np

    allocate(free_x(free_np),free_y(free_np),free_z(free_np),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    ! permutation array
    allocate(free_perm(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    free_perm(:) = 0
    inum = 0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! note: xstore/ystore/zstore have changed coordinates to r/theta/phi,
        !       converts back to x/y/z
        call rthetaphi_2_xyz(x,y,z,xstore_crust_mantle(iglob), &
                             ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))

        free_x(inum) = x
        free_y(inum) = y
        free_z(inum) = z
        ! stores permutation
        free_perm(iglob) = inum
      endif
    enddo
    if (inum /= free_np) stop 'Error free_np count in loading free surface points'

    ! hi/low resolution
    if (VTK_USE_HIRES) then
      ! point connectivity
      free_nspec = NSPEC_TOP*(NGLLX-1)*(NGLLY-1)

      allocate(free_conn(4,free_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      inum = 0
      free_conn(:,:) = -1
      do ispec2D = 1,NSPEC_TOP
        ispec = ibelm_top_crust_mantle(ispec2D)
        k = NGLLZ
        do j = 1, NGLLY-1
          do i = 1, NGLLX-1
            ! indices of corner points
            id1 = free_perm(ibool_crust_mantle(i,j,k,ispec))
            id2 = free_perm(ibool_crust_mantle(i+1,j,k,ispec))
            id3 = free_perm(ibool_crust_mantle(i+1,j+1,k,ispec))
            id4 = free_perm(ibool_crust_mantle(i,j+1,k,ispec))
            ! note: indices for VTK start at 0
            inum = inum+1
            free_conn(1,inum) = id1 - 1
            free_conn(2,inum) = id2 - 1
            free_conn(3,inum) = id3 - 1
            free_conn(4,inum) = id4 - 1
          enddo
        enddo
      enddo
    else
      ! point connectivity
      free_nspec = NSPEC_TOP

      allocate(free_conn(4,free_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      inum = 0
      free_conn(:,:) = -1
      do ispec2D = 1,NSPEC_TOP
        ispec = ibelm_top_crust_mantle(ispec2D)
        ! indices of corner points
        id1 = free_perm(ibool_crust_mantle(1,1,NGLLZ,ispec))
        id2 = free_perm(ibool_crust_mantle(NGLLX,1,NGLLZ,ispec))
        id3 = free_perm(ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,ispec))
        id4 = free_perm(ibool_crust_mantle(1,NGLLY,NGLLZ,ispec))
        ! note: indices for VTK start at 0
        inum = inum + 1
        free_conn(1,inum) = id1 - 1
        free_conn(2,inum) = id2 - 1
        free_conn(3,inum) = id3 - 1
        free_conn(4,inum) = id4 - 1
      enddo
    endif
    if (minval(free_conn(:,:)) < 0) stop 'Error VTK free surface point connectivity'

    ! gathers data from all MPI processes
    if (NPROC > 1) then
      ! multiple MPI processes

      ! user output
      !if (myrank == 0 ) print*,"    gathering all MPI info... "

      ! number of volume points for all partitions together
      call sum_all_i(free_np,free_np_all)
      if (myrank == 0 ) write(IMAIN,*) "    all freesurface points: ",free_np_all

      ! gathers point info
      allocate(free_points_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_points_all(:) = 0
      call gather_all_singlei(free_np,free_points_all,NPROC)

      ! array offsets
      allocate(free_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_offset_all(1) = 0
      do i = 2, NPROC
        free_offset_all(i) = sum(free_points_all(1:i-1))
      enddo

      ! number of volume elements
      call sum_all_i(free_nspec,free_nspec_all)
      if (myrank == 0 ) write(IMAIN,*) "    all freesurface elements: ",free_nspec_all

      ! freesurface elements
      allocate(free_conn_nspec_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_conn_nspec_all(:) = 0
      call gather_all_singlei(4*free_nspec,free_conn_nspec_all,NPROC)

      ! array offsets
      allocate(free_conn_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      free_conn_offset_all(1) = 0
      do i = 2, NPROC
        free_conn_offset_all(i) = sum(free_conn_nspec_all(1:i-1))
      enddo

      ! global data arrays (only needed on master process)
      if (myrank == 0) then
        ! gather locations
        allocate(free_x_all(free_np_all), &
                 free_y_all(free_np_all), &
                 free_z_all(free_np_all),stat=ier )
        if (ier /= 0 ) stop 'Error allocating free_x_all,... arrays'

        free_x_all(:) = 0.0
        free_y_all(:) = 0.0
        free_z_all(:) = 0.0

        ! connectivity
        allocate(free_conn_all(4,free_nspec_all),stat=ier)
        if (ier /= 0 ) stop 'Error allocating free_conn_all array'
        free_conn_all(:,:) = 0
      endif

      if (myrank == 0) then
        ! locations
        !if (myrank == 0 ) print*,"    locations..."
        call gatherv_all_r(free_x,free_np, &
                            free_x_all,free_points_all,free_offset_all, &
                            free_np_all,NPROC)
        call gatherv_all_r(free_y,free_np, &
                            free_y_all,free_points_all,free_offset_all, &
                            free_np_all,NPROC)
        call gatherv_all_r(free_z,free_np, &
                            free_z_all,free_points_all,free_offset_all, &
                            free_np_all,NPROC)

        ! connectivity
        !if (myrank == 0 ) print*,"    connectivity..."
        call gatherv_all_i(free_conn,4*free_nspec, &
                           free_conn_all,free_conn_nspec_all,free_conn_offset_all, &
                           free_nspec_all,NPROC)

        ! shifts connectivity ids for all additional slices
        do i = 2, NPROC
          ! divides by 4 to get nspec numbers
          ispec_start = free_conn_offset_all(i)/4 + 1
          ispec_end = free_conn_offset_all(i)/4 + free_conn_nspec_all(i)/4
          do ispec = ispec_start,ispec_end
            free_conn_all(:,ispec) = free_conn_all(:,ispec) + free_offset_all(i)
          enddo
        enddo

        !if (myrank == 0 ) print*,"    preparing VTK field..."

        ! adds free surface to VTK window
        call prepare_vtkfreesurface(free_np_all,free_x_all,free_y_all,free_z_all, &
                                    free_nspec_all,free_conn_all)

      else
        ! all other process just send data locations
        call gatherv_all_r(free_x,free_np, &
                            dummy,free_points_all,free_offset_all, &
                            1,NPROC)
        call gatherv_all_r(free_y,free_np, &
                            dummy,free_points_all,free_offset_all, &
                            1,NPROC)
        call gatherv_all_r(free_z,free_np, &
                            dummy,free_points_all,free_offset_all, &
                            1,NPROC)
        ! connectivity
        call gatherv_all_i(free_conn,4*free_nspec, &
                            dummy_i,free_conn_nspec_all,free_conn_offset_all, &
                            1,NPROC)

      endif
    else
      ! serial run
      ! creates VTK freesurface actor
      call prepare_vtkfreesurface(free_np,free_x,free_y,free_z, &
                                  free_nspec,free_conn)

    endif

    ! frees memory
    deallocate(free_x,free_y,free_z)
    deallocate(free_conn,free_perm)
    if (NPROC > 1) then
      deallocate(free_conn_nspec_all,free_conn_offset_all)
      deallocate(free_points_all,free_offset_all)
      if (myrank == 0 ) deallocate(free_x_all,free_y_all,free_z_all,free_conn_all)
    endif
  endif ! VTK_SHOW_FREESURFACE
  call synchronize_all()

  ! volume data
  if (VTK_SHOW_VOLUME) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  VTK volume:"
      write(IMAIN,*) "    spectral elements    : ",NSPEC_CRUST_MANTLE
    endif

    ! sets new point mask
    vtkmask(:) = .false.
    do ispec = 1,NSPEC_CRUST_MANTLE
      ! hi/low resolution
      ! loops only over points
      do k = 1,NGLLZ,NIT_res
        do j = 1,NGLLY,NIT_res
          do i = 1,NGLLX,NIT_res
            iglob = ibool_crust_mantle(i,j,k,ispec)
            ! sets mask
            vtkmask(iglob) = .true.
          enddo
        enddo
      enddo
    enddo
    vol_np = count(vtkmask(:))

    ! loads volume data arrays
    if (myrank == 0 ) write(IMAIN,*) "    loading volume points: ",vol_np

    allocate(vol_x(vol_np),vol_y(vol_np),vol_z(vol_np),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    ! permutation array
    allocate(vol_perm(NGLOB_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    vol_perm(:) = 0
    inum = 0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if (vtkmask(iglob) .eqv. .true.) then
        inum = inum + 1
        ! note: xstore/ystore/zstore have changed coordinates to r/theta/phi,
        !       converts back to x/y/z
        call rthetaphi_2_xyz(x,y,z,xstore_crust_mantle(iglob), &
                             ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        vol_x(inum) = x
        vol_y(inum) = y
        vol_z(inum) = z
        ! stores permutation
        vol_perm(iglob) = inum
      endif
    enddo
    if (inum /= vol_np) stop 'Error vol_np count in loading volume points'

    ! hi/low resolution
    if (VTK_USE_HIRES) then
      ! point connectivity
      vol_nspec = NSPEC_CRUST_MANTLE*(NGLLX-1)*(NGLLY-1)*(NGLLZ-1)

      allocate(vol_conn(8,vol_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      inum = 0
      vol_conn(:,:) = -1
      do ispec = 1,NSPEC_CRUST_MANTLE
        do k = 1, NGLLZ-1
          do j = 1, NGLLY-1
            do i = 1, NGLLX-1
              ! indices of corner points
              id1 = vol_perm(ibool_crust_mantle(i,j,k,ispec))
              id2 = vol_perm(ibool_crust_mantle(i+1,j,k,ispec))
              id3 = vol_perm(ibool_crust_mantle(i+1,j+1,k,ispec))
              id4 = vol_perm(ibool_crust_mantle(i,j+1,k,ispec))

              id5 = vol_perm(ibool_crust_mantle(i,j,k+1,ispec))
              id6 = vol_perm(ibool_crust_mantle(i+1,j,k+1,ispec))
              id7 = vol_perm(ibool_crust_mantle(i+1,j+1,k+1,ispec))
              id8 = vol_perm(ibool_crust_mantle(i,j+1,k+1,ispec))

              ! note: indices for VTK start at 0
              inum = inum+1
              vol_conn(1,inum) = id1 - 1
              vol_conn(2,inum) = id2 - 1
              vol_conn(3,inum) = id3 - 1
              vol_conn(4,inum) = id4 - 1
              vol_conn(5,inum) = id5 - 1
              vol_conn(6,inum) = id6 - 1
              vol_conn(7,inum) = id7 - 1
              vol_conn(8,inum) = id8 - 1
            enddo
          enddo
        enddo
      enddo
    else
      ! point connectivity
      vol_nspec = NSPEC_CRUST_MANTLE

      allocate(vol_conn(8,vol_nspec),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vol_conn(:,:) = -1
      do ispec = 1,NSPEC_CRUST_MANTLE
        ! indices of corner points
        id1 = vol_perm(ibool_crust_mantle(1,1,1,ispec))
        id2 = vol_perm(ibool_crust_mantle(NGLLX,1,1,ispec))
        id3 = vol_perm(ibool_crust_mantle(NGLLX,NGLLY,1,ispec))
        id4 = vol_perm(ibool_crust_mantle(1,NGLLY,1,ispec))

        id5 = vol_perm(ibool_crust_mantle(1,1,NGLLZ,ispec))
        id6 = vol_perm(ibool_crust_mantle(NGLLX,1,NGLLZ,ispec))
        id7 = vol_perm(ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,ispec))
        id8 = vol_perm(ibool_crust_mantle(1,NGLLY,NGLLZ,ispec))

        ! note: indices for VTK start at 0
        vol_conn(1,ispec) = id1 - 1
        vol_conn(2,ispec) = id2 - 1
        vol_conn(3,ispec) = id3 - 1
        vol_conn(4,ispec) = id4 - 1
        vol_conn(5,ispec) = id5 - 1
        vol_conn(6,ispec) = id6 - 1
        vol_conn(7,ispec) = id7 - 1
        vol_conn(8,ispec) = id8 - 1
      enddo
    endif
    if (minval(vol_conn(:,:)) < 0) stop 'Error VTK volume point connectivity'

    ! allocates local data array
    allocate(vtkdata(vol_np),stat=ier)
    if (ier /= 0 ) stop 'Error allocating arrays'

    vtkdata(:) = 0.0

    ! gathers data from all MPI processes
    if (NPROC > 1) then
      ! multiple MPI processes

      ! user output
      !if (myrank == 0 ) print*,"    gathering all MPI info... "

      ! number of volume points for all partitions together
      call sum_all_i(vol_np,vtkdata_numpoints_all)
      if (myrank == 0 ) write(IMAIN,*) "    all volume points: ",vtkdata_numpoints_all

      ! gathers point info
      allocate(vtkdata_points_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vtkdata_points_all(:) = 0
      call gather_all_singlei(vol_np,vtkdata_points_all,NPROC)

      ! array offsets
      allocate(vtkdata_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vtkdata_offset_all(1) = 0
      do i = 2, NPROC
        vtkdata_offset_all(i) = sum(vtkdata_points_all(1:i-1))
      enddo

      ! number of volume elements
      call sum_all_i(vol_nspec,vol_nspec_all)
      if (myrank == 0 ) write(IMAIN,*) "    all volume elements: ",vol_nspec_all

      ! volume elements
      allocate(vol_conn_nspec_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vol_conn_nspec_all(:) = 0
      call gather_all_singlei(8*vol_nspec,vol_conn_nspec_all,NPROC)

      ! array offsets
      allocate(vol_conn_offset_all(NPROC),stat=ier)
      if (ier /= 0 ) stop 'Error allocating arrays'

      vol_conn_offset_all(1) = 0
      do i = 2, NPROC
        vol_conn_offset_all(i) = sum(vol_conn_nspec_all(1:i-1))
      enddo

      ! global data arrays (only needed on master process)
      if (myrank == 0) then
        ! point data
        allocate(vtkdata_all(vtkdata_numpoints_all),stat=ier)
        if (ier /= 0 ) stop 'Error allocating vtkdata_all array'

        vtkdata_all(:) = 0.0

        ! gather locations
        allocate(vol_x_all(vtkdata_numpoints_all), &
                 vol_y_all(vtkdata_numpoints_all), &
                 vol_z_all(vtkdata_numpoints_all),stat=ier )
        if (ier /= 0 ) stop 'Error allocating vol_x_all,... arrays'

        vol_x_all(:) = 0.0
        vol_y_all(:) = 0.0
        vol_z_all(:) = 0.0

        ! connectivity
        allocate(vol_conn_all(8,vol_nspec_all),stat=ier)
        if (ier /= 0 ) stop 'Error allocating vol_conn_all array'

        vol_conn_all(:,:) = 0

      endif

      if (myrank == 0) then
        ! locations
        !if (myrank == 0 ) print*,"    locations..."
        call gatherv_all_r(vol_x,vol_np, &
                            vol_x_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)
        call gatherv_all_r(vol_y,vol_np, &
                            vol_y_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)
        call gatherv_all_r(vol_z,vol_np, &
                            vol_z_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROC)

        ! connectivity
        !if (myrank == 0 ) print*,"    connectivity..."
        call gatherv_all_i(vol_conn,8*vol_nspec, &
                           vol_conn_all,vol_conn_nspec_all,vol_conn_offset_all, &
                           vol_nspec_all,NPROC)

        ! shifts connectivity ids for all additional slices
        do i = 2, NPROC
          ! divides by 8 to get nspec numbers
          ispec_start = vol_conn_offset_all(i)/8 + 1
          ispec_end = vol_conn_offset_all(i)/8 + vol_conn_nspec_all(i)/8
          do ispec = ispec_start,ispec_end
            vol_conn_all(:,ispec) = vol_conn_all(:,ispec) + vtkdata_offset_all(i)
          enddo
        enddo

        !if (myrank == 0 ) print*,"    preparing VTK field..."

        ! adds total volume wavefield to VTK window
        call prepare_vtkfield(vtkdata_numpoints_all,vol_x_all,vol_y_all,vol_z_all, &
                              vol_nspec_all,vol_conn_all)

      else
        ! all other process just send data
        ! locations
        call gatherv_all_r(vol_x,vol_np, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
        call gatherv_all_r(vol_y,vol_np, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
        call gatherv_all_r(vol_z,vol_np, &
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROC)
        ! connectivity
        call gatherv_all_i(vol_conn,8*vol_nspec, &
                            dummy_i,vol_conn_nspec_all,vol_conn_offset_all, &
                            1,NPROC)
      endif

    else
      ! serial run
      !if (myrank == 0 ) print*,"    preparing VTK field..."

      ! adds volume wavefield to VTK window
      call prepare_vtkfield(vol_np,vol_x,vol_y,vol_z,vol_nspec,vol_conn)
    endif

    ! frees memory
    deallocate(vol_x,vol_y,vol_z)
    deallocate(vol_conn,vol_perm)
    if (NPROC > 1) then
      deallocate(vol_conn_nspec_all,vol_conn_offset_all)
      if (myrank == 0 ) deallocate(vol_x_all,vol_y_all,vol_z_all,vol_conn_all)
    endif
  endif ! VTK_SHOW_VOLUME

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  VTK visualization preparation done"
    call flush_IMAIN()
  endif
  call synchronize_all()

  end subroutine prepare_vtk_window

