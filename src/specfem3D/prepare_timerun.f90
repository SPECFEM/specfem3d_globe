!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine prepare_timerun()

  use specfem_par
  use specfem_par_movie
  implicit none

  include 'mpif.h'

  ! get MPI starting time
  time_start = MPI_WTIME()

  ! user output infos
  call prepare_timerun_user_output()

  ! sets up mass matrices
  call prepare_timerun_mass_matrices()

  ! convert x/y/z into r/theta/phi spherical coordinates
  call prepare_timerun_convert_coord()

  ! allocate files to save movies
  ! for noise tomography, store_val_x/y/z/ux/uy/uz needed for 'surface movie'
  if(MOVIE_SURFACE .or. NOISE_TOMOGRAPHY /= 0 ) then
    call prepare_timerun_movie_surface()
  endif

  ! output point and element information for 3D movies
  if(MOVIE_VOLUME) call prepare_timerun_movie_volume()

  ! sets up time increments and rotation constants
  call prepare_timerun_constants()

  ! precomputes gravity factors
  call prepare_timerun_gravity()

  ! precomputes attenuation factors
  if(ATTENUATION_VAL) call prepare_timerun_attenuation()

  ! initializes arrays
  call prepare_timerun_init_wavefield()

  ! reads files back from local disk or MT tape system if restart file
  ! note: for SIMULATION_TYPE 3 simulations, the stored wavefields
  !          will be read in the time loop after the Newmark time scheme update.
  !          this makes indexing and timing easier to match with adjoint wavefields indexing.
  call read_forward_arrays_startrun()

  ! prepares noise simulations
  call prepare_timerun_noise()

  ! prepares GPU arrays
  if(GPU_MODE) call prepare_timerun_GPU()

  ! user output
  if( myrank == 0 ) then
    ! elapsed time since beginning of mesh generation
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for preparing timerun in seconds = ',sngl(tCPU)
    write(IMAIN,*)
  endif

  end subroutine prepare_timerun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_user_output()

  use specfem_par
  implicit none

  ! user output
  if(myrank == 0) then

    write(IMAIN,*)
    write(IMAIN,*) 'Reference radius of the Earth used is ',R_EARTH_KM,' km'
    write(IMAIN,*)

    write(IMAIN,*)
    if(OCEANS_VAL) then
      write(IMAIN,*) 'incorporating the oceans using equivalent load'
    else
      write(IMAIN,*) 'no oceans'
    endif

    write(IMAIN,*)
    if(ELLIPTICITY_VAL) then
      write(IMAIN,*) 'incorporating ellipticity'
    else
      write(IMAIN,*) 'no ellipticity'
    endif

    write(IMAIN,*)
    if(TOPOGRAPHY) then
      write(IMAIN,*) 'incorporating surface topography'
    else
      write(IMAIN,*) 'no surface topography'
    endif

    write(IMAIN,*)
    if(GRAVITY_VAL) then
      write(IMAIN,*) 'incorporating self-gravitation (Cowling approximation)'
    else
      write(IMAIN,*) 'no self-gravitation'
    endif

    write(IMAIN,*)
    if(ROTATION_VAL) then
      write(IMAIN,*) 'incorporating rotation'
    else
      write(IMAIN,*) 'no rotation'
    endif

    write(IMAIN,*)
    if(ATTENUATION_VAL) then
      write(IMAIN,*) 'incorporating attenuation using ',N_SLS,' standard linear solids'

      if(ATTENUATION_3D_VAL) write(IMAIN,*) 'using 3D attenuation'

      if(USE_ATTENUATION_MIMIC ) write(IMAIN,*) 'mimicking effects on velocity only'
    else
      write(IMAIN,*) 'no attenuation'
    endif

    write(IMAIN,*)
    write(IMAIN,*)
    write(IMAIN,*)

  endif

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

  ! mass matrices need to be assembled with MPI here once and for all
  call prepare_timerun_rmass_assembly()

  ! check that all the mass matrices are positive
  if(OCEANS_VAL) then
    if(minval(rmass_ocean_load) <= 0.) call exit_MPI(myrank,'negative mass matrix term for the oceans')
  endif
  if(minval(rmass_crust_mantle) <= 0.) call exit_MPI(myrank,'negative mass matrix term for the crust_mantle')
  if(minval(rmass_inner_core) <= 0.) call exit_MPI(myrank,'negative mass matrix term for the inner core')
  if(minval(rmass_outer_core) <= 0.) call exit_MPI(myrank,'negative mass matrix term for the outer core')

  ! for efficiency, invert final mass matrix once and for all on each slice
  if(OCEANS_VAL) rmass_ocean_load = 1._CUSTOM_REAL / rmass_ocean_load

  rmass_crust_mantle = 1._CUSTOM_REAL / rmass_crust_mantle
  rmass_outer_core = 1._CUSTOM_REAL / rmass_outer_core
  rmass_inner_core = 1._CUSTOM_REAL / rmass_inner_core

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

  ! local parameters
  integer :: ndim_assemble

  ! temporary buffers for send and receive between faces of the slices and the chunks
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC) ::  &
    buffer_send_faces_scalar,buffer_received_faces_scalar

  ! synchronize all the processes before assembling the mass matrix
  ! to make sure all the nodes have finished to read their databases
  call sync_all()

  ! the mass matrix needs to be assembled with MPI here once and for all

  ! ocean load
  if (OCEANS_VAL) then
    call assemble_MPI_scalar_block(myrank,rmass_ocean_load,NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_XY,NCHUNKS_VAL)
  endif

  ! crust and mantle
  call assemble_MPI_scalar_block(myrank,rmass_crust_mantle,NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_XY,NCHUNKS_VAL)

  ! outer core
  call assemble_MPI_scalar_block(myrank,rmass_outer_core,NGLOB_OUTER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
            npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
            iboolfaces_outer_core,iboolcorner_outer_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE),NGLOB2DMAX_XY,NCHUNKS_VAL)

  ! inner core
  call assemble_MPI_scalar_block(myrank,rmass_inner_core,NGLOB_INNER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE),NGLOB2DMAX_XY,NCHUNKS_VAL)


  ! mass matrix including central cube
  if(INCLUDE_CENTRAL_CUBE) then
    ! the mass matrix to assemble is a scalar, not a vector
    ndim_assemble = 1

    ! use central cube buffers to assemble the inner core mass matrix with the central cube
    call assemble_MPI_central_cube_block(ichunk,nb_msgs_theor_in_cube, sender_from_slices_to_cube, &
                 npoin2D_cube_from_slices, buffer_all_cube_from_slices, &
                 buffer_slices, buffer_slices2, ibool_central_cube, &
                 receiver_cube_from_slices, ibool_inner_core, &
                 idoubling_inner_core, NSPEC_INNER_CORE, &
                 ibelm_bottom_inner_core, NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                 NGLOB_INNER_CORE, &
                 rmass_inner_core,ndim_assemble)

    ! suppress fictitious mass matrix elements in central cube
    ! because the slices do not compute all their spectral elements in the cube
    where(rmass_inner_core(:) <= 0.) rmass_inner_core = 1.
  endif

  if(myrank == 0) write(IMAIN,*) 'end assembling MPI mass matrix'

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
    call xyz_2_rthetaphi(xstore_crust_mantle(i), &
                        ystore_crust_mantle(i), &
                        zstore_crust_mantle(i),rval,thetaval,phival)
    xstore_crust_mantle(i) = rval
    ystore_crust_mantle(i) = thetaval
    zstore_crust_mantle(i) = phival
  enddo

  ! convert in the outer core
  do i = 1,NGLOB_OUTER_CORE
    call xyz_2_rthetaphi(xstore_outer_core(i), &
                        ystore_outer_core(i), &
                        zstore_outer_core(i),rval,thetaval,phival)
    xstore_outer_core(i) = rval
    ystore_outer_core(i) = thetaval
    zstore_outer_core(i) = phival
  enddo

  ! convert in the inner core
  do i = 1,NGLOB_INNER_CORE
    call xyz_2_rthetaphi(xstore_inner_core(i), &
                        ystore_inner_core(i), &
                        zstore_inner_core(i),rval,thetaval,phival)
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

  if(MOVIE_COARSE .and. NOISE_TOMOGRAPHY ==0) then  ! only output corners !for noise tomography, must NOT be coarse
     nmovie_points = 2 * 2 * NSPEC2D_TOP(IREGION_CRUST_MANTLE)
     if(NGLLX /= NGLLY) &
      call exit_MPI(myrank,'MOVIE_COARSE together with MOVIE_SURFACE requires NGLLX=NGLLY')
     NIT = NGLLX - 1
  else
     nmovie_points = NGLLX * NGLLY * NSPEC2D_TOP(IREGION_CRUST_MANTLE)
     NIT = 1
  endif
  allocate(store_val_x(nmovie_points), &
          store_val_y(nmovie_points), &
          store_val_z(nmovie_points), &
          store_val_ux(nmovie_points), &
          store_val_uy(nmovie_points), &
          store_val_uz(nmovie_points),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating movie surface arrays')

  if (MOVIE_SURFACE) then  ! those arrays are not neccessary for noise tomography, so only allocate them in MOVIE_SURFACE case
     allocate(store_val_x_all(nmovie_points,0:NPROCTOT_VAL-1), &
            store_val_y_all(nmovie_points,0:NPROCTOT_VAL-1), &
            store_val_z_all(nmovie_points,0:NPROCTOT_VAL-1), &
            store_val_ux_all(nmovie_points,0:NPROCTOT_VAL-1), &
            store_val_uy_all(nmovie_points,0:NPROCTOT_VAL-1), &
            store_val_uz_all(nmovie_points,0:NPROCTOT_VAL-1),stat=ier)
     if( ier /= 0 ) call exit_MPI(myrank,'error allocating movie surface all arrays')
  endif
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Movie surface:'
    write(IMAIN,*) '  Writing to moviedata*** files in output directory'
    if(MOVIE_VOLUME_TYPE == 5) then
      write(IMAIN,*) '  movie output: displacement'
    else
      write(IMAIN,*) '  movie output: velocity'
    endif
    write(IMAIN,*) '  time steps every: ',NTSTEP_BETWEEN_FRAMES
  endif

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

  ! the following has to be true for the the array dimensions of eps to match with those of xstore etc..
  ! note that epsilondev and eps_trace_over_3 don't have the same dimensions.. could cause trouble
  if (NSPEC_CRUST_MANTLE_STR_OR_ATT /= NSPEC_CRUST_MANTLE) &
    stop 'NSPEC_CRUST_MANTLE_STRAINS_ATT /= NSPEC_CRUST_MANTLE'
  if (NSPEC_CRUST_MANTLE_STRAIN_ONLY /= NSPEC_CRUST_MANTLE) &
    stop 'NSPEC_CRUST_MANTLE_STRAIN_ONLY /= NSPEC_CRUST_MANTLE'

  write(prname,'(a,i6.6,a)') trim(LOCAL_PATH)//'/'//'proc',myrank,'_'
  call count_points_movie_volume(prname,ibool_crust_mantle, xstore_crust_mantle,ystore_crust_mantle, &
              zstore_crust_mantle,MOVIE_TOP,MOVIE_BOTTOM,MOVIE_WEST,MOVIE_EAST,MOVIE_NORTH,MOVIE_SOUTH, &
              MOVIE_COARSE,npoints_3dmovie,nspecel_3dmovie,num_ibool_3dmovie,mask_ibool,mask_3dmovie)

  allocate(nu_3dmovie(3,3,npoints_3dmovie),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating nu for 3d movie')

  call write_movie_volume_mesh(npoints_3dmovie,prname,ibool_crust_mantle,xstore_crust_mantle, &
                         ystore_crust_mantle,zstore_crust_mantle, muvstore_crust_mantle_3dmovie, &
                         mask_3dmovie,mask_ibool,num_ibool_3dmovie,nu_3dmovie,MOVIE_COARSE)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Movie volume:'
    write(IMAIN,*) '  Writing to movie3D*** files on local disk databases directory'
    if(MOVIE_VOLUME_TYPE == 1) then
      write(IMAIN,*) '  movie output: strain'
    else if(MOVIE_VOLUME_TYPE == 2) then
      write(IMAIN,*) '  movie output: time integral of strain'
    else if(MOVIE_VOLUME_TYPE == 3) then
      write(IMAIN,*) '  movie output: potency or integral of strain'
    else if(MOVIE_VOLUME_TYPE == 4) then
      write(IMAIN,*) '  movie output: divergence and curl'
    else if(MOVIE_VOLUME_TYPE == 5) then
      write(IMAIN,*) '  movie output: displacement'
    else if(MOVIE_VOLUME_TYPE == 6) then
      write(IMAIN,*) '  movie output: velocity'
    endif
    write(IMAIN,*) '  depth(T,B):',MOVIE_TOP,MOVIE_BOTTOM
    write(IMAIN,*) '  lon(W,E)  :',MOVIE_WEST,MOVIE_EAST
    write(IMAIN,*) '  lat(S,N)  :',MOVIE_SOUTH,MOVIE_NORTH
    write(IMAIN,*) '  Starting at time step:',MOVIE_START, 'ending at:',MOVIE_STOP,'every: ',NTSTEP_BETWEEN_FRAMES
  endif

  if( MOVIE_VOLUME_TYPE < 1 .or. MOVIE_VOLUME_TYPE > 6) &
      call exit_MPI(myrank, 'MOVIE_VOLUME_TYPE has to be 1,2,3,4,5 or 6')

  end subroutine prepare_timerun_movie_volume

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_constants()

! precomputes constants for time integration

  use specfem_par
  implicit none

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(((NSTEP-1)*DT-t0)/60.d0),' minutes'
    write(IMAIN,*) 'start time:',sngl(-t0),' seconds'
    write(IMAIN,*)
  endif

  ! define constants for the time integration
  ! scaling to make displacement in meters and velocity in meters per second
  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)
  scale_t_inv = dsqrt(PI*GRAV*RHOAV)

  scale_displ = R_EARTH

  scale_veloc = scale_displ * scale_t_inv

  ! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    deltat = sngl(DT*scale_t_inv)
  else
    deltat = DT*scale_t_inv
  endif
  deltatover2 = 0.5d0*deltat
  deltatsqover2 = 0.5d0*deltat*deltat

  if (SIMULATION_TYPE == 3) then
    if(CUSTOM_REAL == SIZE_REAL) then
      b_deltat = - sngl(DT*scale_t_inv)
    else
      b_deltat = - DT*scale_t_inv
    endif
    b_deltatover2 = 0.5d0*b_deltat
    b_deltatsqover2 = 0.5d0*b_deltat*b_deltat
  endif

  ! non-dimensionalized rotation rate of the Earth times two
  if(ROTATION_VAL) then
    ! distinguish between single and double precision for reals
    if (SIMULATION_TYPE == 1) then
      if(CUSTOM_REAL == SIZE_REAL) then
        two_omega_earth = sngl(2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv))
      else
        two_omega_earth = 2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv)
      endif
    else
      if(CUSTOM_REAL == SIZE_REAL) then
        two_omega_earth = - sngl(2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv))
      else
        two_omega_earth = - 2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv)
      endif
    endif

    A_array_rotation = 0._CUSTOM_REAL
    B_array_rotation = 0._CUSTOM_REAL

    if (SIMULATION_TYPE == 3) then
      if(CUSTOM_REAL == SIZE_REAL) then
        b_two_omega_earth = sngl(2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv))
      else
        b_two_omega_earth = 2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv)
      endif
    endif
  else
    two_omega_earth = 0._CUSTOM_REAL
    if (SIMULATION_TYPE == 3) b_two_omega_earth = 0._CUSTOM_REAL
  endif


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

  ! store g, rho and dg/dr=dg using normalized radius in lookup table every 100 m
  ! get density and velocity from PREM model using dummy doubling flag
  ! this assumes that the gravity perturbations are small and smooth
  ! and that we can neglect the 3D model and use PREM every 100 m in all cases
  ! this is probably a rather reasonable assumption
  if(GRAVITY_VAL) then
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
      if(radius_km > RCMB/1000.d0 - 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RCMB/1000.d0 - 3.d0)*10.d0))
      if(radius_km < RICB/1000.d0 + 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RICB/1000.d0 + 3.d0)*10.d0))
    enddo

    ! compute gravity value at CMB and ICB once and for all
    radius = RCMB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_cmb_dble)

    radius = RICB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_icb_dble)

    ! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      minus_g_cmb = sngl(- g_cmb_dble)
      minus_g_icb = sngl(- g_icb_dble)
    else
      minus_g_cmb = - g_cmb_dble
      minus_g_icb = - g_icb_dble
    endif

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

  end subroutine prepare_timerun_gravity


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_attenuation()

  ! precomputes attenuation factors

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_movie
  implicit none

  ! local parameters
  double precision, dimension(ATT1,ATT2,ATT3,ATT4) :: omsb_crust_mantle_dble, factor_scale_crust_mantle_dble
  double precision, dimension(ATT1,ATT2,ATT3,ATT5) :: omsb_inner_core_dble, factor_scale_inner_core_dble
  double precision, dimension(N_SLS,ATT1,ATT2,ATT3,ATT4) :: factor_common_crust_mantle_dble
  double precision, dimension(N_SLS,ATT1,ATT2,ATT3,ATT5) :: factor_common_inner_core_dble
  double precision, dimension(N_SLS) :: alphaval_dble, betaval_dble, gammaval_dble
  double precision, dimension(N_SLS) :: tau_sigma_dble

  double precision :: scale_factor,scale_factor_minus_one
  real(kind=CUSTOM_REAL) :: mul
  integer :: ispec,i,j,k
  character(len=150) :: prnamel

  ! get and store PREM attenuation model

  ! CRUST_MANTLE ATTENUATION
  call create_name_database(prnamel, myrank, IREGION_CRUST_MANTLE, LOCAL_PATH)
  call get_attenuation_model_3D(myrank, prnamel, omsb_crust_mantle_dble, &
           factor_common_crust_mantle_dble,factor_scale_crust_mantle_dble,tau_sigma_dble,NSPEC_CRUST_MANTLE)

  ! INNER_CORE ATTENUATION
  call create_name_database(prnamel, myrank, IREGION_INNER_CORE, LOCAL_PATH)
  call get_attenuation_model_3D(myrank, prnamel, omsb_inner_core_dble, &
           factor_common_inner_core_dble,factor_scale_inner_core_dble,tau_sigma_dble,NSPEC_INNER_CORE)

  if(CUSTOM_REAL == SIZE_REAL) then
    factor_scale_crust_mantle       = sngl(factor_scale_crust_mantle_dble)
    one_minus_sum_beta_crust_mantle = sngl(omsb_crust_mantle_dble)
    factor_common_crust_mantle      = sngl(factor_common_crust_mantle_dble)

    factor_scale_inner_core         = sngl(factor_scale_inner_core_dble)
    one_minus_sum_beta_inner_core   = sngl(omsb_inner_core_dble)
    factor_common_inner_core        = sngl(factor_common_inner_core_dble)
  else
    factor_scale_crust_mantle       = factor_scale_crust_mantle_dble
    one_minus_sum_beta_crust_mantle = omsb_crust_mantle_dble
    factor_common_crust_mantle      = factor_common_crust_mantle_dble

    factor_scale_inner_core         = factor_scale_inner_core_dble
    one_minus_sum_beta_inner_core   = omsb_inner_core_dble
    factor_common_inner_core        = factor_common_inner_core_dble
  endif

  ! if attenuation is on, shift PREM to right frequency
  ! rescale mu in PREM to average frequency for attenuation
  ! the formulas to implement the scaling can be found for instance in
  ! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
  ! anelasticity: implications for seismology and mantle composition,
  ! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
  ! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
  ! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170

  ! rescale in crust and mantle

  do ispec = 1,NSPEC_CRUST_MANTLE
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          scale_factor = factor_scale_crust_mantle(i,j,k,ispec)

          if(ANISOTROPIC_3D_MANTLE_VAL) then
            scale_factor_minus_one = scale_factor - 1.
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
            if(MOVIE_VOLUME .and. SIMULATION_TYPE==3) then
              ! store the original value of \mu to comput \mu*\eps
              muvstore_crust_mantle_3dmovie(i,j,k,ispec)=muvstore_crust_mantle(i,j,k,ispec)
            endif
            muvstore_crust_mantle(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec) * scale_factor

            ! scales transverse isotropic values for mu_h
            !if(TRANSVERSE_ISOTROPY_VAL .and. (idoubling_crust_mantle(ispec) == IFLAG_220_80 &
            !    .or. idoubling_crust_mantle(ispec) == IFLAG_80_MOHO)) &
            if( ispec_is_tiso_crust_mantle(ispec) ) then
              muhstore_crust_mantle(i,j,k,ispec) = muhstore_crust_mantle(i,j,k,ispec) * scale_factor
            endif
          endif

        enddo
      enddo
    enddo
  enddo ! END DO CRUST MANTLE

  ! rescale in inner core

  do ispec = 1,NSPEC_INNER_CORE
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          scale_factor_minus_one = factor_scale_inner_core(i,j,k,ispec) - 1.0

          if(ANISOTROPIC_INNER_CORE_VAL) then
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

          muvstore_inner_core(i,j,k,ispec) = muvstore_inner_core(i,j,k,ispec) * factor_scale_inner_core(i,j,k,ispec)

        enddo
      enddo
    enddo
  enddo ! END DO INNER CORE

  ! precompute Runge-Kutta coefficients
  call get_attenuation_memory_values(tau_sigma_dble, deltat, alphaval_dble, betaval_dble, gammaval_dble)
  if(CUSTOM_REAL == SIZE_REAL) then
    alphaval = sngl(alphaval_dble)
    betaval  = sngl(betaval_dble)
    gammaval = sngl(gammaval_dble)
  else
    alphaval = alphaval_dble
    betaval  = betaval_dble
    gammaval = gammaval_dble
  endif

  if (SIMULATION_TYPE == 3) then
   call get_attenuation_memory_values(tau_sigma_dble, b_deltat, alphaval_dble, betaval_dble, gammaval_dble)
   if(CUSTOM_REAL == SIZE_REAL) then
     b_alphaval = sngl(alphaval_dble)
     b_betaval  = sngl(betaval_dble)
     b_gammaval = sngl(gammaval_dble)
   else
     b_alphaval = alphaval_dble
     b_betaval  = betaval_dble
     b_gammaval = gammaval_dble
   endif
  endif

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

  ! initialize arrays to zero
  displ_crust_mantle(:,:) = 0._CUSTOM_REAL
  veloc_crust_mantle(:,:) = 0._CUSTOM_REAL
  accel_crust_mantle(:,:) = 0._CUSTOM_REAL

  displ_outer_core(:) = 0._CUSTOM_REAL
  veloc_outer_core(:) = 0._CUSTOM_REAL
  accel_outer_core(:) = 0._CUSTOM_REAL

  displ_inner_core(:,:) = 0._CUSTOM_REAL
  veloc_inner_core(:,:) = 0._CUSTOM_REAL
  accel_inner_core(:,:) = 0._CUSTOM_REAL

  ! put negligible initial value to avoid very slow underflow trapping
  if(FIX_UNDERFLOW_PROBLEM) then
    displ_crust_mantle(:,:) = VERYSMALLVAL
    displ_outer_core(:) = VERYSMALLVAL
    displ_inner_core(:,:) = VERYSMALLVAL
  endif

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

  if (SIMULATION_TYPE == 3) then
    rho_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    beta_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    alpha_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    if (NOISE_TOMOGRAPHY == 3) Sigma_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL

    ! approximate hessian
    if( APPROXIMATE_HESS_KL ) then
      allocate( hess_kl_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
      if( ier /= 0 ) call exit_MPI(myrank,'error allocating hessian')
      hess_kl_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
    endif

    ! For anisotropic kernels (in crust_mantle only)
    cijkl_kl_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL

    rho_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL
    alpha_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL

    rho_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL
    beta_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL
    alpha_kl_inner_core(:,:,:,:) = 0._CUSTOM_REAL

    div_displ_outer_core(:,:,:,:) = 0._CUSTOM_REAL
    b_div_displ_outer_core(:,:,:,:) = 0._CUSTOM_REAL

    ! deviatoric kernel check
    if( deviatoric_outercore) then
      nspec_beta_kl_outer_core = NSPEC_OUTER_CORE_ADJOINT
    else
      nspec_beta_kl_outer_core = 1
    endif
    allocate(beta_kl_outer_core(NGLLX,NGLLY,NGLLZ,nspec_beta_kl_outer_core),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating beta outercore')
    beta_kl_outer_core(:,:,:,:) = 0._CUSTOM_REAL
  endif

  ! initialize to be on the save side for adjoint runs SIMULATION_TYPE==2
  eps_trace_over_3_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL

  epsilondev_xx_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_yy_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_xy_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_xz_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_yz_crust_mantle(:,:,:,:) = 0._CUSTOM_REAL

  eps_trace_over_3_inner_core(:,:,:,:) = 0._CUSTOM_REAL

  epsilondev_xx_inner_core(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_yy_inner_core(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_xy_inner_core(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_xz_inner_core(:,:,:,:) = 0._CUSTOM_REAL
  epsilondev_yz_inner_core(:,:,:,:) = 0._CUSTOM_REAL

  if(FIX_UNDERFLOW_PROBLEM) then
    eps_trace_over_3_crust_mantle(:,:,:,:) = VERYSMALLVAL

    epsilondev_xx_crust_mantle(:,:,:,:) = VERYSMALLVAL
    epsilondev_yy_crust_mantle(:,:,:,:) = VERYSMALLVAL
    epsilondev_xy_crust_mantle(:,:,:,:) = VERYSMALLVAL
    epsilondev_xz_crust_mantle(:,:,:,:) = VERYSMALLVAL
    epsilondev_yz_crust_mantle(:,:,:,:) = VERYSMALLVAL

    eps_trace_over_3_inner_core(:,:,:,:) = VERYSMALLVAL

    epsilondev_xx_inner_core(:,:,:,:) = VERYSMALLVAL
    epsilondev_yy_inner_core(:,:,:,:) = VERYSMALLVAL
    epsilondev_xy_inner_core(:,:,:,:) = VERYSMALLVAL
    epsilondev_xz_inner_core(:,:,:,:) = VERYSMALLVAL
    epsilondev_yz_inner_core(:,:,:,:) = VERYSMALLVAL

  endif

  if (COMPUTE_AND_STORE_STRAIN) then
    if(MOVIE_VOLUME .and. (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3)) then
      Iepsilondev_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
      Ieps_trace_over_3_crust_mantle(:,:,:,:)=0._CUSTOM_REAL
    endif
  endif

  ! clear memory variables if attenuation
  if(ATTENUATION_VAL) then
    R_xx_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yy_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xy_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xz_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yz_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL

    R_xx_inner_core(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yy_inner_core(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xy_inner_core(:,:,:,:,:) = 0._CUSTOM_REAL
    R_xz_inner_core(:,:,:,:,:) = 0._CUSTOM_REAL
    R_yz_inner_core(:,:,:,:,:) = 0._CUSTOM_REAL

    if(FIX_UNDERFLOW_PROBLEM) then
      R_xx_crust_mantle(:,:,:,:,:) = VERYSMALLVAL
      R_yy_crust_mantle(:,:,:,:,:) = VERYSMALLVAL
      R_xy_crust_mantle(:,:,:,:,:) = VERYSMALLVAL
      R_xz_crust_mantle(:,:,:,:,:) = VERYSMALLVAL
      R_yz_crust_mantle(:,:,:,:,:) = VERYSMALLVAL

      R_xx_inner_core(:,:,:,:,:) = VERYSMALLVAL
      R_yy_inner_core(:,:,:,:,:) = VERYSMALLVAL
      R_xy_inner_core(:,:,:,:,:) = VERYSMALLVAL
      R_xz_inner_core(:,:,:,:,:) = VERYSMALLVAL
      R_yz_inner_core(:,:,:,:,:) = VERYSMALLVAL
    endif
  endif

  end subroutine prepare_timerun_init_wavefield


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_noise()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  implicit none
  ! local parameters
  integer :: ier

  ! NOISE TOMOGRAPHY
  if ( NOISE_TOMOGRAPHY /= 0 ) then
    nspec_top = NSPEC2D_TOP(IREGION_CRUST_MANTLE)

    allocate(noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP), &
            normal_x_noise(nmovie_points), &
            normal_y_noise(nmovie_points), &
            normal_z_noise(nmovie_points), &
            mask_noise(nmovie_points), &
            noise_surface_movie(NDIM,NGLLX,NGLLY,nspec_top),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating noise arrays')

    noise_sourcearray(:,:,:,:,:) = 0._CUSTOM_REAL
    normal_x_noise(:)            = 0._CUSTOM_REAL
    normal_y_noise(:)            = 0._CUSTOM_REAL
    normal_z_noise(:)            = 0._CUSTOM_REAL
    mask_noise(:)                = 0._CUSTOM_REAL
    noise_surface_movie(:,:,:,:) = 0._CUSTOM_REAL

    call read_parameters_noise(myrank,nrec,NSTEP,nmovie_points, &
                              islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu, &
                              noise_sourcearray,xigll,yigll,zigll,nspec_top, &
                              NIT, ibool_crust_mantle, ibelm_top_crust_mantle, &
                              xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                              irec_master_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise)

    call check_parameters_noise(myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                              NUMBER_OF_RUNS, NUMBER_OF_THIS_RUN,ROTATE_SEISMOGRAMS_RT, &
                              SAVE_ALL_SEISMOS_IN_ONE_FILE, USE_BINARY_FOR_LARGE_FILE, &
                              MOVIE_COARSE,LOCAL_PATH,nspec_top,NSTEP)
  endif

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
  include 'mpif.h'

  ! local parameters
  integer :: ier
  real :: free_mb,used_mb,total_mb
  integer :: ncuda_devices,ncuda_devices_min,ncuda_devices_max
  ! dummy custom_real variables to convert from double precision
  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable:: cr_wgll_cube
  real(kind=CUSTOM_REAL),dimension(:),allocatable:: &
    cr_d_ln_density_dr_table,cr_minus_rho_g_over_kappa_fluid, &
    cr_minus_gravity_table,cr_minus_deriv_gravity_table, &
    cr_density_table

  ! GPU_MODE now defined in Par_file
  if(myrank == 0 ) then
    write(IMAIN,*)
    write(IMAIN,*) "GPU_MODE Active. Preparing Fields and Constants on Device."
    write(IMAIN,*)
  endif

  ! initializes GPU and outputs info to files for all processes
  call prepare_cuda_device(myrank,ncuda_devices)

  ! collects min/max of local devices found for statistics
  call MPI_REDUCE(ncuda_devices,ncuda_devices_min,1,MPI_INTEGER,MPI_MIN,0,MPI_COMM_WORLD,ier)
  call MPI_REDUCE(ncuda_devices,ncuda_devices_max,1,MPI_INTEGER,MPI_MAX,0,MPI_COMM_WORLD,ier)

  ! prepares general fields on GPU
  call prepare_constants_device(Mesh_pointer,myrank,NGLLX, &
                                  hprime_xx, hprime_yy, hprime_zz, &
                                  hprimewgll_xx, hprimewgll_yy, hprimewgll_zz, &
                                  wgllwgll_xy, wgllwgll_xz, wgllwgll_yz, &
                                  NSOURCES, nsources_local, &
                                  sourcearrays,islice_selected_source,ispec_selected_source, &
                                  number_receiver_global,islice_selected_rec,ispec_selected_rec, &
                                  nrec, nrec_local, nadj_rec_local, &
                                  NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                                  NGLOB_CRUST_MANTLE_OCEANS, &
                                  NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                                  NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                  SIMULATION_TYPE,NOISE_TOMOGRAPHY, &
                                  SAVE_FORWARD,ABSORBING_CONDITIONS, &
                                  OCEANS_VAL,GRAVITY_VAL,ROTATION_VAL, &
                                  ATTENUATION_VAL,USE_ATTENUATION_MIMIC, &
                                  COMPUTE_AND_STORE_STRAIN, &
                                  ANISOTROPIC_3D_MANTLE_VAL,ANISOTROPIC_INNER_CORE_VAL, &
                                  SAVE_BOUNDARY_MESH, &
                                  USE_MESH_COLORING_GPU, &
                                  ANISOTROPIC_KL,APPROXIMATE_HESS_KL)
  call sync_all()

  ! prepares rotation arrays
  if( ROTATION_VAL ) then
    if(myrank == 0 ) write(IMAIN,*) "  loading rotation arrays"

    call prepare_fields_rotation_device(Mesh_pointer, &
                                  two_omega_earth,deltat, &
                                  A_array_rotation,B_array_rotation, &
                                  b_two_omega_earth,b_deltat, &
                                  b_A_array_rotation,b_B_array_rotation, &
                                  NSPEC_OUTER_CORE_ROTATION)
  endif
  call sync_all()

  ! prepares arrays related to gravity
  ! note: GPU will use only single-precision (or double precision) for all calculations
  !          we convert to wgll_cube to custom real (by default single-precision),
  !          using implicit conversion
  if(myrank == 0 ) write(IMAIN,*) "  loading non-gravity/gravity arrays"

  allocate(cr_d_ln_density_dr_table(NRAD_GRAVITY), &
          cr_minus_rho_g_over_kappa_fluid(NRAD_GRAVITY), &
          cr_minus_gravity_table(NRAD_GRAVITY), &
          cr_minus_deriv_gravity_table(NRAD_GRAVITY), &
          cr_density_table(NRAD_GRAVITY), &
          stat=ier)
  if( ier /= 0 ) stop 'error allocating cr_minus_rho_g_over_kappa_fluid, etc...'
  ! d_ln_density_dr_table needed for no gravity case
  cr_d_ln_density_dr_table(:) = d_ln_density_dr_table(:)
  ! these are needed for gravity cases only
  cr_minus_rho_g_over_kappa_fluid(:) = minus_rho_g_over_kappa_fluid(:)
  cr_minus_gravity_table(:) = minus_gravity_table(:)
  cr_minus_deriv_gravity_table(:) = minus_deriv_gravity_table(:)
  cr_density_table(:) = density_table(:)

  allocate(cr_wgll_cube(NGLLX,NGLLY,NGLLZ),stat=ier)
  if( ier /= 0 ) stop 'error allocating cr_wgll_cube'
  cr_wgll_cube(:,:,:) = wgll_cube(:,:,:)

  ! prepares on GPU
  call prepare_fields_gravity_device(Mesh_pointer, &
                                    cr_d_ln_density_dr_table, &
                                    cr_minus_rho_g_over_kappa_fluid, &
                                    cr_minus_gravity_table, &
                                    cr_minus_deriv_gravity_table, &
                                    cr_density_table, &
                                    cr_wgll_cube, &
                                    NRAD_GRAVITY)
  deallocate(cr_d_ln_density_dr_table,cr_minus_rho_g_over_kappa_fluid, &
            cr_minus_gravity_table,cr_minus_deriv_gravity_table, &
            cr_density_table)
  deallocate(cr_wgll_cube)
  call sync_all()

  ! prepares attenuation arrays
  if( ATTENUATION_VAL ) then
    if(myrank == 0 ) write(IMAIN,*) "  loading attenuation"

    call prepare_fields_attenuat_device(Mesh_pointer, &
                                        R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                        R_xz_crust_mantle,R_yz_crust_mantle, &
                                        factor_common_crust_mantle, &
                                        one_minus_sum_beta_crust_mantle, &
                                        R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                        R_xz_inner_core,R_yz_inner_core, &
                                        factor_common_inner_core, &
                                        one_minus_sum_beta_inner_core, &
                                        alphaval,betaval,gammaval, &
                                        b_alphaval,b_betaval,b_gammaval)
  endif
  call sync_all()


  ! prepares attenuation arrays
  if( COMPUTE_AND_STORE_STRAIN ) then
    if(myrank == 0 ) write(IMAIN,*) "  loading strain"

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
  endif
  call sync_all()

  ! prepares absorbing arrays
  if(NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) then
    if(myrank == 0 ) write(IMAIN,*) "  loading absorbing boundaries"

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
                                    rho_vp_crust_mantle,rho_vs_crust_mantle,  &
                                    nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                                    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                                    nspec2D_zmin_outer_core, &
                                    NSPEC2DMAX_XMIN_XMAX_OC,NSPEC2DMAX_YMIN_YMAX_OC, &
                                    nimin_outer_core,nimax_outer_core, &
                                    njmin_outer_core,njmax_outer_core, &
                                    nkmin_xi_outer_core,nkmin_eta_outer_core, &
                                    ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
                                    ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
                                    ibelm_bottom_outer_core, &
                                    jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
                                    jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
                                    jacobian2D_bottom_outer_core, &
                                    vp_outer_core)

  endif
  call sync_all()

  ! prepares MPI interfaces
  if(myrank == 0 ) write(IMAIN,*) "  loading mpi interfaces"

  call prepare_mpi_buffers_device(Mesh_pointer, &
                                num_interfaces_crust_mantle,max_nibool_interfaces_crust_mantle, &
                                nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                                num_interfaces_inner_core,max_nibool_interfaces_inner_core, &
                                nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                                num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                                nibool_interfaces_outer_core,ibool_interfaces_outer_core)

  ! prepares fields on GPU for noise simulations
  if ( NOISE_TOMOGRAPHY > 0 ) then
    if(myrank == 0 ) write(IMAIN,*) "  loading noise arrays"

    call prepare_fields_noise_device(Mesh_pointer,nspec_top,ibelm_top_crust_mantle, &
                                    NSTEP,noise_sourcearray, &
                                    normal_x_noise,normal_y_noise,normal_z_noise, &
                                    mask_noise,jacobian2D_top_crust_mantle)

  endif

  ! prepares oceans arrays
  if ( OCEANS_VAL ) then
     if(myrank == 0 ) write(IMAIN,*) "  loading oceans arrays"

     call prepare_oceans_device(Mesh_pointer,rmass_ocean_load)

  endif

  ! crust/mantle region
  if(myrank == 0 ) write(IMAIN,*) "  loading crust/mantle region"
  call prepare_crust_mantle_device(Mesh_pointer, &
                                  xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                  etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                  gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                  rhostore_crust_mantle, &
                                  kappavstore_crust_mantle,muvstore_crust_mantle, &
                                  kappahstore_crust_mantle,muhstore_crust_mantle, &
                                  eta_anisostore_crust_mantle, &
                                  rmass_crust_mantle, &
                                  normal_top_crust_mantle, &
                                  ibelm_top_crust_mantle, &
                                  ibelm_bottom_crust_mantle, &
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
                                  NSPEC2D_TOP(IREGION_CRUST_MANTLE), &
                                  NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE))
  call sync_all()


  ! outer core region
  if(myrank == 0 ) write(IMAIN,*) "  loading outer core region"
  call prepare_outer_core_device(Mesh_pointer, &
                                  xix_outer_core,xiy_outer_core,xiz_outer_core, &
                                  etax_outer_core,etay_outer_core,etaz_outer_core, &
                                  gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                                  rhostore_outer_core,kappavstore_outer_core, &
                                  rmass_outer_core, &
                                  normal_top_outer_core, &
                                  normal_bottom_outer_core, &
                                  jacobian2D_top_outer_core, &
                                  jacobian2D_bottom_outer_core, &
                                  ibelm_top_outer_core, &
                                  ibelm_bottom_outer_core, &
                                  ibool_outer_core, &
                                  xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                                  num_phase_ispec_outer_core,phase_ispec_inner_outer_core, &
                                  nspec_outer_outer_core,nspec_inner_outer_core, &
                                  NSPEC2D_TOP(IREGION_OUTER_CORE), &
                                  NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
  call sync_all()


  ! inner core region
  if(myrank == 0 ) write(IMAIN,*) "  loading inner core region"
  call prepare_inner_core_device(Mesh_pointer, &
                                  xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                  etax_inner_core,etay_inner_core,etaz_inner_core, &
                                  gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                  rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                  rmass_inner_core, &
                                  ibelm_top_inner_core, &
                                  ibool_inner_core, &
                                  xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                                  c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                                  c33store_inner_core,c44store_inner_core, &
                                  idoubling_inner_core, &
                                  num_phase_ispec_inner_core,phase_ispec_inner_inner_core, &
                                  nspec_outer_inner_core,nspec_inner_inner_core, &
                                  NSPEC2D_TOP(IREGION_INNER_CORE))
  call sync_all()

  ! transfer forward and backward fields to device with initial values
  if(myrank == 0 ) write(IMAIN,*) "  transfering initial wavefield"
  call transfer_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                                   Mesh_pointer)

  call transfer_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,displ_inner_core,veloc_inner_core,accel_inner_core, &
                                   Mesh_pointer)

  call transfer_fields_oc_to_device(NGLOB_OUTER_CORE,displ_outer_core,veloc_outer_core,accel_outer_core, &
                                   Mesh_pointer)

  if(SIMULATION_TYPE == 3) then
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE, &
                                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                    Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE, &
                                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                    Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE, &
                                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                    Mesh_pointer)
  endif


  ! outputs GPU usage to files for all processes
  call output_free_device_memory(myrank)

  ! outputs usage for main process
  if( myrank == 0 ) then
    write(IMAIN,*)"  GPU number of devices per node: min =",ncuda_devices_min
    write(IMAIN,*)"                                  max =",ncuda_devices_max
    write(IMAIN,*)

    call get_free_device_memory(free_mb,used_mb,total_mb)
    write(IMAIN,*)"  GPU usage: free  =",free_mb," MB",nint(free_mb/total_mb*100.0),"%"
    write(IMAIN,*)"             used  =",used_mb," MB",nint(used_mb/total_mb*100.0),"%"
    write(IMAIN,*)"             total =",total_mb," MB",nint(total_mb/total_mb*100.0),"%"
    write(IMAIN,*)
  endif

  end subroutine prepare_timerun_GPU
