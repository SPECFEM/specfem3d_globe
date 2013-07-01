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

  subroutine iterate_time()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  implicit none

  ! timing
  double precision, external :: wtime

!
!   s t a r t   t i m e   i t e r a t i o n s
!

  ! synchronize all processes to make sure everybody is ready to start time loop
  call sync_all()

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Starting time iteration loop...'
    write(IMAIN,*)
  endif

  ! create an empty file to monitor the start of the simulation
  if(myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/starttimeloop.txt',status='unknown',action='write')
    write(IOUT,*) 'hello, starting time loop'
    close(IOUT)
  endif

  ! initialize variables for writing seismograms
  seismo_offset = it_begin-1
  seismo_current = 0

#ifdef _HANDOPT
  imodulo_NGLOB_CRUST_MANTLE = mod(NGLOB_CRUST_MANTLE,3)
  imodulo_NGLOB_CRUST_MANTLE4 = mod(NGLOB_CRUST_MANTLE,4)
  imodulo_NGLOB_INNER_CORE = mod(NGLOB_INNER_CORE,3)
  imodulo_NGLOB_OUTER_CORE = mod(NGLOB_OUTER_CORE,3)
#endif

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  do it = it_begin,it_end

    ! simulation status output and stability check
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call it_check_stability()
    endif

    ! update displacement using Newmark time scheme
    call it_update_displacement_scheme()

    ! acoustic solver for outer core
    ! (needs to be done first, before elastic one)
    call compute_forces_acoustic()

    ! elastic solver for crust/mantle and inner core
    call compute_forces_elastic()

    ! restores last time snapshot saved for backward/reconstruction of wavefields
    ! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
    !          and adjoint sources will become more complicated
    !          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
    if( SIMULATION_TYPE == 3 .and. it == 1 ) then
      call read_forward_arrays()
    endif

    ! write the seismograms with time shift
    if( nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) ) then
      call write_seismograms()
    endif

    ! adjoint simulations: kernels
    if( SIMULATION_TYPE == 3 ) then
      call compute_kernels()
    endif

    ! outputs movie files
    call write_movie_output()

    ! first step of noise tomography, i.e., save a surface movie at every time step
    ! modified from the subroutine 'write_movie_surface'
    if( NOISE_TOMOGRAPHY == 1 ) then
      call noise_save_surface_movie()
    endif

    ! updates vtk window
    if( VTK_MODE ) then
      call it_update_vtkwindow()
    endif

  enddo   ! end of main time loop

  !
  !---- end of time iteration loop
  !

  call it_print_elapsed_time()

  ! Transfer fields from GPU card to host for further analysis
  if(GPU_MODE) call it_transfer_from_GPU()

  end subroutine iterate_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_update_displacement_scheme()

! explicit Newmark time scheme with acoustic & elastic domains:
! (see e.g. Hughes, 1987; Chaljub et al., 2003)
!
! chi(t+delta_t) = chi(t) + delta_t chi_dot(t) + 1/2 delta_t**2 chi_dot_dot(t)
! chi_dot(t+delta_t) = chi_dot(t) + 1/2 delta_t chi_dot_dot(t) + 1/2 delta_t chi_dot_dot(t+delta_t)
! chi_dot_dot(t+delta_t) = 1/M_acoustic( -K_acoustic chi(t+delta) + B_acoustic u(t+delta_t) + f(t+delta_t) )
!
! u(t+delta_t) = u(t) + delta_t  v(t) + 1/2  delta_t**2 a(t)
! v(t+delta_t) = v(t) + 1/2 delta_t a(t) + 1/2 delta_t a(t+delta_t)
! a(t+delta_t) = 1/M_elastic ( -K_elastic u(t+delta) + B_elastic chi_dot_dot(t+delta_t) + f( t+delta_t) )
!
! where
!   chi, chi_dot, chi_dot_dot are acoustic (fluid) potentials ( dotted with respect to time)
!   u, v, a are displacement,velocity & acceleration
!   M is mass matrix, K stiffness matrix and B boundary term for acoustic/elastic domains
!   f denotes a source term (acoustic/elastic)
!
! note that this stage calculates the predictor terms
!
!   for
!   potential chi_dot(t+delta) requires + 1/2 delta_t chi_dot_dot(t+delta_t)
!                                   at a later stage (corrector) once where chi_dot_dot(t+delta) is calculated
!   and similar,
!   velocity v(t+delta_t) requires  + 1/2 delta_t a(t+delta_t)
!                                   at a later stage once where a(t+delta) is calculated
! also:
!   boundary term B_elastic requires chi_dot_dot(t+delta)
!                                   thus chi_dot_dot has to be updated first before the elastic boundary term is considered

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  ! local parameters
  integer :: i

  ! updates wavefields
  if( .not. GPU_MODE) then
  ! on CPU

    ! Newmark time scheme update
#ifdef _HANDOPT_NEWMARK
! way 2:
! One common technique in computational science to help enhance pipelining is loop unrolling
!
! we're accessing NDIM=3 components at each line,
! that is, for an iteration, the register must contain
! NDIM * displ_ + NDIM * veloc_ + NDIM * accel + deltat + deltatsq..
! in most cases a real (CUSTOM_REAL) value will have 4 bytes,
! assuming a default cache size of about 128 bytes, we unroll here in steps of 3, thus 29 reals or 118 bytes,
! rather than with steps of 4
    ! mantle
    if(imodulo_NGLOB_CRUST_MANTLE >= 1) then
      do i = 1,imodulo_NGLOB_CRUST_MANTLE
        displ_crust_mantle(:,i) = displ_crust_mantle(:,i) &
          + deltat*veloc_crust_mantle(:,i) + deltatsqover2*accel_crust_mantle(:,i)

        veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) &
          + deltatover2*accel_crust_mantle(:,i)

        accel_crust_mantle(:,i) = 0._CUSTOM_REAL
      enddo
    endif
    do i = imodulo_NGLOB_CRUST_MANTLE+1,NGLOB_CRUST_MANTLE, 3 ! in steps of 3
      displ_crust_mantle(:,i) = displ_crust_mantle(:,i) &
        + deltat*veloc_crust_mantle(:,i) + deltatsqover2*accel_crust_mantle(:,i)
      displ_crust_mantle(:,i+1) = displ_crust_mantle(:,i+1) &
        + deltat*veloc_crust_mantle(:,i+1) + deltatsqover2*accel_crust_mantle(:,i+1)
      displ_crust_mantle(:,i+2) = displ_crust_mantle(:,i+2) &
        + deltat*veloc_crust_mantle(:,i+2) + deltatsqover2*accel_crust_mantle(:,i+2)


      veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) &
        + deltatover2*accel_crust_mantle(:,i)
      veloc_crust_mantle(:,i+1) = veloc_crust_mantle(:,i+1) &
        + deltatover2*accel_crust_mantle(:,i+1)
      veloc_crust_mantle(:,i+2) = veloc_crust_mantle(:,i+2) &
        + deltatover2*accel_crust_mantle(:,i+2)

      ! set acceleration to zero
      ! note: we do initialize acceleration in this loop since it is read already into the cache,
      !           otherwise it would have to be read in again for this explicitly,
      !           which would make this step more expensive
      accel_crust_mantle(:,i) = 0._CUSTOM_REAL
      accel_crust_mantle(:,i+1) = 0._CUSTOM_REAL
      accel_crust_mantle(:,i+2) = 0._CUSTOM_REAL
    enddo

    ! outer core
    do i=1,NGLOB_OUTER_CORE
      displ_outer_core(i) = displ_outer_core(i) &
        + deltat*veloc_outer_core(i) + deltatsqover2*accel_outer_core(i)

      veloc_outer_core(i) = veloc_outer_core(i) &
        + deltatover2*accel_outer_core(i)

      accel_outer_core(i) = 0._CUSTOM_REAL
    enddo

    ! inner core
    if(imodulo_NGLOB_INNER_CORE >= 1) then
      do i = 1,imodulo_NGLOB_INNER_CORE
        displ_inner_core(:,i) = displ_inner_core(:,i) &
          + deltat*veloc_inner_core(:,i) + deltatsqover2*accel_inner_core(:,i)

        veloc_inner_core(:,i) = veloc_inner_core(:,i) &
          + deltatover2*accel_inner_core(:,i)

        accel_inner_core(:,i) = 0._CUSTOM_REAL
      enddo
    endif
    do i = imodulo_NGLOB_INNER_CORE+1,NGLOB_INNER_CORE, 3 ! in steps of 3
      displ_inner_core(:,i) = displ_inner_core(:,i) &
        + deltat*veloc_inner_core(:,i) + deltatsqover2*accel_inner_core(:,i)
      displ_inner_core(:,i+1) = displ_inner_core(:,i+1) &
        + deltat*veloc_inner_core(:,i+1) + deltatsqover2*accel_inner_core(:,i+1)
      displ_inner_core(:,i+2) = displ_inner_core(:,i+2) &
        + deltat*veloc_inner_core(:,i+2) + deltatsqover2*accel_inner_core(:,i+2)


      veloc_inner_core(:,i) = veloc_inner_core(:,i) &
        + deltatover2*accel_inner_core(:,i)
      veloc_inner_core(:,i+1) = veloc_inner_core(:,i+1) &
        + deltatover2*accel_inner_core(:,i+1)
      veloc_inner_core(:,i+2) = veloc_inner_core(:,i+2) &
        + deltatover2*accel_inner_core(:,i+2)

      accel_inner_core(:,i) = 0._CUSTOM_REAL
      accel_inner_core(:,i+1) = 0._CUSTOM_REAL
      accel_inner_core(:,i+2) = 0._CUSTOM_REAL
    enddo
#else
! way 1:
    ! mantle
    do i=1,NGLOB_CRUST_MANTLE
      displ_crust_mantle(:,i) = displ_crust_mantle(:,i) &
        + deltat*veloc_crust_mantle(:,i) + deltatsqover2*accel_crust_mantle(:,i)
      veloc_crust_mantle(:,i) = veloc_crust_mantle(:,i) &
        + deltatover2*accel_crust_mantle(:,i)
      accel_crust_mantle(:,i) = 0._CUSTOM_REAL
    enddo
    ! outer core
    do i=1,NGLOB_OUTER_CORE
      displ_outer_core(i) = displ_outer_core(i) &
        + deltat*veloc_outer_core(i) + deltatsqover2*accel_outer_core(i)
      veloc_outer_core(i) = veloc_outer_core(i) &
        + deltatover2*accel_outer_core(i)
      accel_outer_core(i) = 0._CUSTOM_REAL
    enddo
    ! inner core
    do i=1,NGLOB_INNER_CORE
      displ_inner_core(:,i) = displ_inner_core(:,i) &
        + deltat*veloc_inner_core(:,i) + deltatsqover2*accel_inner_core(:,i)
      veloc_inner_core(:,i) = veloc_inner_core(:,i) &
        + deltatover2*accel_inner_core(:,i)
      accel_inner_core(:,i) = 0._CUSTOM_REAL
    enddo
#endif

    ! backward field
    if (SIMULATION_TYPE == 3) then

#ifdef _HANDOPT_NEWMARK
! way 2:
      ! mantle
      if(imodulo_NGLOB_CRUST_MANTLE >= 1) then
        do i=1,imodulo_NGLOB_CRUST_MANTLE
          b_displ_crust_mantle(:,i) = b_displ_crust_mantle(:,i) &
            + b_deltat*b_veloc_crust_mantle(:,i) + b_deltatsqover2*b_accel_crust_mantle(:,i)
          b_veloc_crust_mantle(:,i) = b_veloc_crust_mantle(:,i) &
            + b_deltatover2*b_accel_crust_mantle(:,i)
          b_accel_crust_mantle(:,i) = 0._CUSTOM_REAL
        enddo
      endif
      do i=imodulo_NGLOB_CRUST_MANTLE+1,NGLOB_CRUST_MANTLE,3
        b_displ_crust_mantle(:,i) = b_displ_crust_mantle(:,i) &
          + b_deltat*b_veloc_crust_mantle(:,i) + b_deltatsqover2*b_accel_crust_mantle(:,i)
        b_displ_crust_mantle(:,i+1) = b_displ_crust_mantle(:,i+1) &
          + b_deltat*b_veloc_crust_mantle(:,i+1) + b_deltatsqover2*b_accel_crust_mantle(:,i+1)
        b_displ_crust_mantle(:,i+2) = b_displ_crust_mantle(:,i+2) &
          + b_deltat*b_veloc_crust_mantle(:,i+2) + b_deltatsqover2*b_accel_crust_mantle(:,i+2)


        b_veloc_crust_mantle(:,i) = b_veloc_crust_mantle(:,i) &
          + b_deltatover2*b_accel_crust_mantle(:,i)
        b_veloc_crust_mantle(:,i+1) = b_veloc_crust_mantle(:,i+1) &
          + b_deltatover2*b_accel_crust_mantle(:,i+1)
        b_veloc_crust_mantle(:,i+2) = b_veloc_crust_mantle(:,i+2) &
          + b_deltatover2*b_accel_crust_mantle(:,i+2)

        b_accel_crust_mantle(:,i) = 0._CUSTOM_REAL
        b_accel_crust_mantle(:,i+1) = 0._CUSTOM_REAL
        b_accel_crust_mantle(:,i+2) = 0._CUSTOM_REAL
      enddo

      ! outer core
      do i=1,NGLOB_OUTER_CORE
        b_displ_outer_core(i) = b_displ_outer_core(i) &
          + b_deltat*b_veloc_outer_core(i) + b_deltatsqover2*b_accel_outer_core(i)
        b_veloc_outer_core(i) = b_veloc_outer_core(i) &
          + b_deltatover2*b_accel_outer_core(i)
        b_accel_outer_core(i) = 0._CUSTOM_REAL
      enddo

      ! inner core
      if(imodulo_NGLOB_INNER_CORE >= 1) then
        do i=1,imodulo_NGLOB_INNER_CORE
          b_displ_inner_core(:,i) = b_displ_inner_core(:,i) &
            + b_deltat*b_veloc_inner_core(:,i) + b_deltatsqover2*b_accel_inner_core(:,i)
          b_veloc_inner_core(:,i) = b_veloc_inner_core(:,i) &
            + b_deltatover2*b_accel_inner_core(:,i)
          b_accel_inner_core(:,i) = 0._CUSTOM_REAL
        enddo
      endif
      do i=imodulo_NGLOB_INNER_CORE+1,NGLOB_INNER_CORE,3
        b_displ_inner_core(:,i) = b_displ_inner_core(:,i) &
          + b_deltat*b_veloc_inner_core(:,i) + b_deltatsqover2*b_accel_inner_core(:,i)
        b_displ_inner_core(:,i+1) = b_displ_inner_core(:,i+1) &
          + b_deltat*b_veloc_inner_core(:,i+1) + b_deltatsqover2*b_accel_inner_core(:,i+1)
        b_displ_inner_core(:,i+2) = b_displ_inner_core(:,i+2) &
          + b_deltat*b_veloc_inner_core(:,i+2) + b_deltatsqover2*b_accel_inner_core(:,i+2)

        b_veloc_inner_core(:,i) = b_veloc_inner_core(:,i) &
          + b_deltatover2*b_accel_inner_core(:,i)
        b_veloc_inner_core(:,i+1) = b_veloc_inner_core(:,i+1) &
          + b_deltatover2*b_accel_inner_core(:,i+1)
        b_veloc_inner_core(:,i+2) = b_veloc_inner_core(:,i+2) &
          + b_deltatover2*b_accel_inner_core(:,i+2)

        b_accel_inner_core(:,i) = 0._CUSTOM_REAL
        b_accel_inner_core(:,i+1) = 0._CUSTOM_REAL
        b_accel_inner_core(:,i+2) = 0._CUSTOM_REAL
      enddo
#else
! way 1:
      ! mantle
      do i=1,NGLOB_CRUST_MANTLE
        b_displ_crust_mantle(:,i) = b_displ_crust_mantle(:,i) &
          + b_deltat*b_veloc_crust_mantle(:,i) + b_deltatsqover2*b_accel_crust_mantle(:,i)
        b_veloc_crust_mantle(:,i) = b_veloc_crust_mantle(:,i) &
          + b_deltatover2*b_accel_crust_mantle(:,i)
        b_accel_crust_mantle(:,i) = 0._CUSTOM_REAL
      enddo
      ! outer core
      do i=1,NGLOB_OUTER_CORE
        b_displ_outer_core(i) = b_displ_outer_core(i) &
          + b_deltat*b_veloc_outer_core(i) + b_deltatsqover2*b_accel_outer_core(i)
        b_veloc_outer_core(i) = b_veloc_outer_core(i) &
          + b_deltatover2*b_accel_outer_core(i)
        b_accel_outer_core(i) = 0._CUSTOM_REAL
      enddo
      ! inner core
      do i=1,NGLOB_INNER_CORE
        b_displ_inner_core(:,i) = b_displ_inner_core(:,i) &
          + b_deltat*b_veloc_inner_core(:,i) + b_deltatsqover2*b_accel_inner_core(:,i)
        b_veloc_inner_core(:,i) = b_veloc_inner_core(:,i) &
          + b_deltatover2*b_accel_inner_core(:,i)
        b_accel_inner_core(:,i) = 0._CUSTOM_REAL
      enddo
#endif
    endif ! SIMULATION_TYPE == 3
  else
    ! on GPU
    ! Includes SIM_TYPE 1 & 3

    ! outer core region
    call it_update_displacement_oc_cuda(Mesh_pointer, &
                                       deltat, deltatsqover2, deltatover2, &
                                       b_deltat, b_deltatsqover2, b_deltatover2)
    ! inner core region
    call it_update_displacement_ic_cuda(Mesh_pointer, &
                                       deltat, deltatsqover2, deltatover2, &
                                       b_deltat, b_deltatsqover2, b_deltatover2)

    ! crust/mantle region
    call it_update_displacement_cm_cuda(Mesh_pointer, &
                                       deltat, deltatsqover2, deltatover2, &
                                       b_deltat, b_deltatsqover2, b_deltatover2)
  endif

  end subroutine it_update_displacement_scheme

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_print_elapsed_time()

  use specfem_par
  implicit none

  ! local parameters
  integer :: ihours,iminutes,iseconds,int_tCPU
  ! timing
  double precision, external :: wtime

  if(myrank == 0) then
     ! elapsed time since beginning of the simulation
     tCPU = wtime() - time_start

     int_tCPU = int(tCPU)
     ihours = int_tCPU / 3600
     iminutes = (int_tCPU - 3600*ihours) / 60
     iseconds = int_tCPU - 3600*ihours - 60*iminutes
     write(IMAIN,*) 'Time-Loop Complete. Timing info:'
     write(IMAIN,*) 'Total elapsed time in seconds = ',tCPU
     write(IMAIN,"(' Total elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
  endif
  end subroutine it_print_elapsed_time

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_check_stability()

! computes the maximum of the norm of the displacement
! in all the slices using an MPI reduction
! and output timestamp file to check that simulation is running fine

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  implicit none

  ! compute the maximum of the norm of the displacement
  ! in all the slices using an MPI reduction
  ! and output timestamp file to check that simulation is running fine
  call check_simulation_stability(it,displ_crust_mantle,displ_inner_core,displ_outer_core, &
                        b_displ_crust_mantle,b_displ_inner_core,b_displ_outer_core, &
                        eps_trace_over_3_crust_mantle, &
                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                        SIMULATION_TYPE,OUTPUT_FILES,time_start,DT,t0,NSTEP, &
                        myrank)

  ! debug output
  !if( maxval(displ_crust_mantle(1,:)**2 + &
  !                displ_crust_mantle(2,:)**2 + displ_crust_mantle(3,:)**2) > 1.e4 ) then
  !  print*,'slice',myrank
  !  print*,'  crust_mantle displ:', maxval(displ_crust_mantle(1,:)), &
  !           maxval(displ_crust_mantle(2,:)),maxval(displ_crust_mantle(3,:))
  !  print*,'  indxs: ',maxloc( displ_crust_mantle(1,:)),maxloc( displ_crust_mantle(2,:)),maxloc( displ_crust_mantle(3,:))
  !  indx = maxloc( displ_crust_mantle(3,:) )
  !  rval = xstore_crust_mantle(indx(1))
  !  thetaval = ystore_crust_mantle(indx(1))
  !  phival = zstore_crust_mantle(indx(1))
  !  !thetaval = PI_OVER_TWO-datan(1.006760466d0*dcos(dble(thetaval))/dmax1(TINYVAL,dsin(dble(thetaval))))
  !  print*,'r/lat/lon:',rval*R_EARTH_KM,90.0-thetaval*180./PI,phival*180./PI
  !  call rthetaphi_2_xyz(rval,thetaval,phival,xstore_crust_mantle(indx(1)),&
  !                     ystore_crust_mantle(indx(1)),zstore_crust_mantle(indx(1)))
  !  print*,'x/y/z:',rval,thetaval,phival
  !  call exit_MPI(myrank,'error stability')
  !endif

  end subroutine it_check_stability

!
!-------------------------------------------------------------------------------------------------
!

  subroutine it_transfer_from_GPU()

! transfers fields on GPU back onto CPU

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  ! to store forward wave fields
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then

    call transfer_fields_cm_from_device(NDIM*NGLOB_CRUST_MANTLE, &
                                    displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle, &
                                    Mesh_pointer)
    call transfer_fields_ic_from_device(NDIM*NGLOB_INNER_CORE, &
                                    displ_inner_core,veloc_inner_core,accel_inner_core, &
                                    Mesh_pointer)
    call transfer_fields_oc_from_device(NGLOB_OUTER_CORE, &
                                    displ_outer_core,veloc_outer_core,accel_outer_core, &
                                    Mesh_pointer)

    call transfer_strain_cm_from_device(Mesh_pointer,eps_trace_over_3_crust_mantle, &
                                    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                    epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                    epsilondev_yz_crust_mantle)
    call transfer_strain_ic_from_device(Mesh_pointer,eps_trace_over_3_inner_core, &
                                    epsilondev_xx_inner_core,epsilondev_yy_inner_core, &
                                    epsilondev_xy_inner_core,epsilondev_xz_inner_core, &
                                    epsilondev_yz_inner_core)

    if (ROTATION_VAL) then
      call transfer_rotation_from_device(Mesh_pointer,A_array_rotation,B_array_rotation)
    endif

    ! note: for kernel simulations (SIMULATION_TYPE == 3), attenuation is
    !          only mimicking effects on phase shifts, but not on amplitudes.
    !          flag USE_ATTENUATION_MIMIC will have to be set to true in this case.
    !
    ! arrays b_R_xx, ... are not used when USE_ATTENUATION_MIMIC is set,
    ! therefore no need to transfer arrays from GPU to CPU
    !if (ATTENUATION) then
    !endif

  else if (SIMULATION_TYPE == 3) then
    ! to store kernels
    call transfer_kernels_oc_to_host(Mesh_pointer, &
                                     rho_kl_outer_core,&
                                     alpha_kl_outer_core,NSPEC_OUTER_CORE)
    call transfer_kernels_cm_to_host(Mesh_pointer, &
                                     rho_kl_crust_mantle, &
                                     alpha_kl_crust_mantle, &
                                     beta_kl_crust_mantle, &
                                     cijkl_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    call transfer_kernels_ic_to_host(Mesh_pointer, &
                                     rho_kl_inner_core, &
                                     alpha_kl_inner_core, &
                                     beta_kl_inner_core,NSPEC_INNER_CORE)

    ! specific noise strength kernel
    if( NOISE_TOMOGRAPHY == 3 ) then
      call transfer_kernels_noise_to_host(Mesh_pointer,Sigma_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif

    ! approximative hessian for preconditioning kernels
    if ( APPROXIMATE_HESS_KL ) then
      call transfer_kernels_hess_cm_tohost(Mesh_pointer,hess_kl_crust_mantle,NSPEC_CRUST_MANTLE)
    endif
  endif

  ! frees allocated memory on GPU
  call prepare_cleanup_device(Mesh_pointer,NCHUNKS_VAL)

  end subroutine it_transfer_from_GPU

!=====================================================================


  subroutine it_update_vtkwindow()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  real :: currenttime
  integer :: iglob,inum
  real(kind=CUSTOM_REAL),dimension(1):: dummy

  ! vtk rendering at frame interval
  if( mod(it,NTSTEP_BETWEEN_FRAMES) == 0 ) then

    ! user output
    !if( myrank == 0 ) print*,"  vtk rendering..."

    ! updates time
    currenttime = sngl((it-1)*DT-t0)

    ! transfers fields from GPU to host
    if( GPU_MODE ) then
      !if( myrank == 0 ) print*,"  vtk: transfering velocity from gpu"
      call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
    endif

    ! updates wavefield
    !if( myrank == 0 ) print*,"  vtk: it = ",it," out of ",it_end," - norm of velocity field"
    inum = 0
    vtkdata(:) = 0.0
    do iglob = 1,NGLOB_CRUST_MANTLE
      if( vtkmask(iglob) .eqv. .true. ) then
        inum = inum + 1
        ! stores norm of velocity vector
        vtkdata(inum) = sqrt(veloc_crust_mantle(1,iglob)**2 &
                           + veloc_crust_mantle(2,iglob)**2 &
                           + veloc_crust_mantle(3,iglob)**2)
      endif
    enddo

    ! updates for multiple mpi process
    if( NPROCTOT_VAL > 1 ) then
      if( myrank == 0 ) then
        ! gather data
        call gatherv_all_cr(vtkdata,size(vtkdata),&
                            vtkdata_all,vtkdata_points_all,vtkdata_offset_all, &
                            vtkdata_numpoints_all,NPROCTOT_VAL)
        ! updates vtk window
        call visualize_vtkdata(it,currenttime,vtkdata_all)
      else
        ! all other process just send data
        call gatherv_all_cr(vtkdata,size(vtkdata),&
                            dummy,vtkdata_points_all,vtkdata_offset_all, &
                            1,NPROCTOT_VAL)
      endif
    else
      ! serial run
      ! updates vtk window
      call visualize_vtkdata(it,currenttime,vtkdata)
    endif

  endif

  end subroutine it_update_vtkwindow

