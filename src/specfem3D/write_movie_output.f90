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

  subroutine write_movie_output()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  implicit none

  ! local parameters
  ! debugging
  character(len=MAX_STRING_LEN) :: filename
  integer,dimension(:),allocatable :: dummy_i

  !-----------------------------------------------------------------------------
  ! user parameters
  ! outputs volume snapshot VTK files of displacement in crust/mantle region for debugging
  logical, parameter :: DEBUG_SNAPSHOT = .false.
  !-----------------------------------------------------------------------------

  character(len=MAX_STRING_LEN) :: command

  ! save movie on surface
  if (MOVIE_SURFACE) then
    if (mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then

      ! gets resulting array values onto CPU
      if (GPU_MODE) then
        ! transfers whole fields
        if (MOVIE_VOLUME_TYPE == 5) then
          call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,Mesh_pointer)
        else
          call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
        endif
      endif

      ! save velocity here to avoid static offset on displacement for movies
      call write_movie_surface()

      ! executes an external script on the node
      if (RUN_EXTERNAL_MOVIE_SCRIPT) then
        ! synchronizes outputs
        call synchronize_all()
        ! calls shell external command
        if (myrank == 0) then
          write(command,"(a,1x,i6.6,' >& out.',i6.6,'.log &')") trim(MOVIE_SCRIPT_NAME),it,it
          !print*,trim(command)
          call system_command(command)
        endif
      endif

    endif
  endif


  ! save movie in full 3D mesh
  if (MOVIE_VOLUME) then

    ! updates integral of strain for adjoint movie volume
    if (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3) then
      ! transfers strain arrays onto CPU
      if (GPU_MODE) then
        call transfer_strain_cm_from_device(Mesh_pointer,eps_trace_over_3_crust_mantle, &
                                         epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                         epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                         epsilondev_yz_crust_mantle)
      endif

      ! integrates strain
      call movie_volume_integrate_strain(deltat,NSPEC_CRUST_MANTLE_3DMOVIE, &
                                        eps_trace_over_3_crust_mantle, &
                                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                        epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                        epsilondev_yz_crust_mantle, &
                                        Ieps_trace_over_3_crust_mantle, &
                                        Iepsilondev_xx_crust_mantle,Iepsilondev_yy_crust_mantle, &
                                        Iepsilondev_xy_crust_mantle,Iepsilondev_xz_crust_mantle, &
                                        Iepsilondev_yz_crust_mantle)
    endif

    ! file output
    if (mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0 &
      .and. it >= MOVIE_START .and. it <= MOVIE_STOP) then

      select case (MOVIE_VOLUME_TYPE)
      case (1)
        ! output strains

        ! gets resulting array values onto CPU
        if (GPU_MODE) then
          call transfer_strain_cm_from_device(Mesh_pointer, &
                                eps_trace_over_3_crust_mantle, &
                                epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                epsilondev_yz_crust_mantle)
        endif

        call  write_movie_volume_strains(myrank,npoints_3dmovie, &
                    LOCAL_TMP_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,eps_trace_over_3_crust_mantle, &
                    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                    muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)

      case (2, 3)
        ! output the Time Integral of Strain, or \mu*TIS
        call  write_movie_volume_strains(myrank,npoints_3dmovie, &
                    LOCAL_TMP_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,Ieps_trace_over_3_crust_mantle, &
                    Iepsilondev_xx_crust_mantle,Iepsilondev_yy_crust_mantle,Iepsilondev_xy_crust_mantle, &
                    Iepsilondev_xz_crust_mantle,Iepsilondev_yz_crust_mantle, &
                    muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)

      case (4)
        ! output divergence and curl in whole volume

        ! gets resulting array values onto CPU
        if (GPU_MODE) then
          ! strains
          call transfer_strain_cm_from_device(Mesh_pointer, &
                                eps_trace_over_3_crust_mantle, &
                                epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle, &
                                epsilondev_xy_crust_mantle,epsilondev_xz_crust_mantle, &
                                epsilondev_yz_crust_mantle)
          call transfer_strain_ic_from_device(Mesh_pointer, &
                                eps_trace_over_3_inner_core, &
                                epsilondev_xx_inner_core,epsilondev_yy_inner_core, &
                                epsilondev_xy_inner_core,epsilondev_xz_inner_core, &
                                epsilondev_yz_inner_core)
          ! wavefields
          call transfer_fields_oc_from_device(NGLOB_OUTER_CORE, &
                                displ_outer_core,veloc_outer_core,accel_outer_core,Mesh_pointer)
        endif

        call write_movie_volume_divcurl(myrank,it,eps_trace_over_3_crust_mantle,&
                        div_displ_outer_core, &
                        accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                        eps_trace_over_3_inner_core, &
                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                        epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                        epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                        LOCAL_TMP_PATH)

      case (5)
        !output displacement
        if (GPU_MODE) then
          call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,Mesh_pointer)
        endif

        scalingval = scale_displ
        call write_movie_volume_vector(myrank,it,npoints_3dmovie, &
                                       LOCAL_TMP_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE,ibool_crust_mantle, &
                                       displ_crust_mantle, &
                                       scalingval,mask_3dmovie,nu_3dmovie)

      case (6)
        !output velocity
        if (GPU_MODE) then
          call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
        endif

        scalingval = scale_veloc
        call write_movie_volume_vector(myrank,it,npoints_3dmovie, &
                                       LOCAL_TMP_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE,ibool_crust_mantle, &
                                       veloc_crust_mantle, &
                                       scalingval,mask_3dmovie,nu_3dmovie)

      case (7)
        ! output norm of displacement

        ! gets resulting array values onto CPU
        if (GPU_MODE) then
          ! displacement wavefields
          call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,Mesh_pointer)
          call transfer_displ_ic_from_device(NDIM*NGLOB_INNER_CORE,displ_inner_core,Mesh_pointer)
          call transfer_displ_oc_from_device(NGLOB_OUTER_CORE,displ_outer_core,Mesh_pointer)
        endif

        call write_movie_volume_displnorm(myrank,it,LOCAL_TMP_PATH, &
                        displ_crust_mantle,displ_inner_core,displ_outer_core, &
                        ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

      case (8)
        ! output norm of velocity

        ! gets resulting array values onto CPU
        if (GPU_MODE) then
          ! velocity wavefields
          call transfer_veloc_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,veloc_crust_mantle,Mesh_pointer)
          call transfer_veloc_ic_from_device(NDIM*NGLOB_INNER_CORE,veloc_inner_core,Mesh_pointer)
          call transfer_veloc_oc_from_device(NGLOB_OUTER_CORE,veloc_outer_core,Mesh_pointer)
        endif

        call write_movie_volume_velnorm(myrank,it,LOCAL_TMP_PATH, &
                        veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                        ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

      case (9)
        ! output norm of acceleration

        ! gets resulting array values onto CPU
        if (GPU_MODE) then
          ! acceleration wavefields
          call transfer_accel_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,accel_crust_mantle,Mesh_pointer)
          call transfer_accel_ic_from_device(NDIM*NGLOB_INNER_CORE,accel_inner_core,Mesh_pointer)
          call transfer_accel_oc_from_device(NGLOB_OUTER_CORE,accel_outer_core,Mesh_pointer)
        endif

        call write_movie_volume_accelnorm(myrank,it,LOCAL_TMP_PATH, &
                        accel_crust_mantle,accel_inner_core,accel_outer_core, &
                        ibool_crust_mantle,ibool_inner_core,ibool_outer_core)

      case default
        call exit_MPI(myrank, 'MOVIE_VOLUME_TYPE has to be in range from 1 to 9')

      end select ! MOVIE_VOLUME_TYPE

      ! executes an external script on the node
      if (RUN_EXTERNAL_MOVIE_SCRIPT) then
        ! synchronizes outputs
        call synchronize_all()
        ! calls shell external command
        if (myrank == 0) then
          write(command,"(a,1x,i6.6,' >& out.',i6.6,'.log &')") trim(MOVIE_SCRIPT_NAME),it,it
          !print*,trim(command)
          call system_command(command)
        endif
      endif

    endif
  endif ! MOVIE_VOLUME

  ! debugging
  if (DEBUG_SNAPSHOT) then
    if (mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0 &
      .and. it >= MOVIE_START .and. it <= MOVIE_STOP) then

      !output displacement
      if (GPU_MODE) then
        call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,Mesh_pointer)
        call transfer_displ_ic_from_device(NDIM*NGLOB_INNER_CORE,displ_inner_core,Mesh_pointer)
      endif

      write(prname,'(a,i6.6,a)') 'OUTPUT_FILES/snapshot_proc',myrank,'_'

      ! VTK file output
      ! displacement values

      ! crust mantle
      allocate(dummy_i(NSPEC_CRUST_MANTLE))
      dummy_i(:) = IFLAG_CRUST

      ! one file per process
      write(filename,'(a,a,i6.6)') prname(1:len_trim(prname)),'reg_1_displ_',it
      call write_VTK_data_cr(dummy_i,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle, &
                          displ_crust_mantle,filename)

      ! backward/reconstructed field
      if (SIMULATION_TYPE == 3) then
        write(filename,'(a,a,i6.6)') prname(1:len_trim(prname)),'reg_1_b_displ_',it
        call write_VTK_data_cr(dummy_i,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                            xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle, &
                            b_displ_crust_mantle,filename)
      endif

      deallocate(dummy_i)

    endif
  endif

  end subroutine write_movie_output
