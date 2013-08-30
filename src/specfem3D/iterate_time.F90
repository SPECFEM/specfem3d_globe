!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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
    call flush_IMAIN()
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

  ! get MPI starting time
  time_start = wtime()

  ! *********************************************************
  ! ************* MAIN LOOP OVER THE TIME STEPS *************
  ! *********************************************************

  do it = it_begin,it_end

    ! simulation status output and stability check
    if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then
      call check_stability()
    endif

    ! update displacement using Newmark time scheme
    call update_displacement_Newmark()

    ! acoustic solver for outer core
    ! (needs to be done first, before elastic one)
    call compute_forces_acoustic()

    ! elastic solver for crust/mantle and inner core
    call compute_forces_viscoelastic()

    ! kernel simulations (forward and adjoint wavefields)
    if( SIMULATION_TYPE == 3 ) then
      ! reconstructs forward wavefields based on last store wavefield data

      ! update displacement using Newmark time scheme
      call update_displacement_Newmark_backward()

      ! acoustic solver for outer core
      ! (needs to be done first, before elastic one)
      call compute_forces_acoustic_backward()

      ! elastic solver for crust/mantle and inner core
      call compute_forces_viscoelastic_backward()

      ! restores last time snapshot saved for backward/reconstruction of wavefields
      ! note: this is done here after the Newmark time scheme, otherwise the indexing for sources
      !          and adjoint sources will become more complicated
      !          that is, index it for adjoint sources will match index NSTEP - 1 for backward/reconstructed wavefields
      if( it == 1 ) then
        call read_forward_arrays()
      endif

      ! adjoint simulations: kernels
      call compute_kernels()
    endif

    ! write the seismograms with time shift
    if( nrec_local > 0 .or. ( WRITE_SEISMOGRAMS_BY_MASTER .and. myrank == 0 ) ) then
      call write_seismograms()
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
    call flush_IMAIN()
  endif

  end subroutine it_print_elapsed_time


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
    !          flag PARTIAL_PHYS_DISPERSION_ONLY will have to be set to true in this case.
    !
    ! arrays b_R_xx, ... are not used when PARTIAL_PHYS_DISPERSION_ONLY is set,
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

