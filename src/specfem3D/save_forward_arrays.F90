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

!> Save the values of the arrays used for forward simulations.
!!
!! Two different options allow to save arrays used for forward simulations.
!!   - setting SAVE_FORWARD to true. That save the last frame.
!!   - setting NUMBERS_OF_RUNS to 2 or 3. This save intermediate frame to
!!     restart the remaining of the simulation later.
!! Moreover, one can choose to save these arrays with the help of the ADIOS
!! library. Set ADIOS_FOR_FORWARD_ARRAYS to true in the Par_file for this
!! feature.
  subroutine save_forward_arrays()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) outputname

  ! checks if anything to do
  if (UNDO_ATTENUATION) return

  ! save files to local disk or tape system if restart file
  if (NUMBER_OF_RUNS > 1 .and. NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  saving restart/checkpoint file for run ',NUMBER_OF_THIS_RUN
      call flush_IMAIN()
    endif

    ! note: we do not need to store the seismograms(:,:,:) array.
    !       the seismograms will be outputted at the end of each run. consecutive runs will append to previous seismograms.
    !       thus, it is similar to having NTSTEP_BETWEEN_OUTPUT_SEISMOS set to the number of timesteps of each run.
    !
    ! saves checkpoint
    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call save_intermediate_forward_arrays_adios()
    else
      write(outputname,"('dump_all_arrays',i6.6)") myrank
      open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname), &
           status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file dump_all_arrays*** for writing')

      ! wavefield
      write(IOUT) displ_crust_mantle
      write(IOUT) veloc_crust_mantle
      write(IOUT) accel_crust_mantle

      write(IOUT) displ_inner_core
      write(IOUT) veloc_inner_core
      write(IOUT) accel_inner_core

      write(IOUT) displ_outer_core
      write(IOUT) veloc_outer_core
      write(IOUT) accel_outer_core

      ! strain
      write(IOUT) epsilondev_xx_crust_mantle
      write(IOUT) epsilondev_yy_crust_mantle
      write(IOUT) epsilondev_xy_crust_mantle
      write(IOUT) epsilondev_xz_crust_mantle
      write(IOUT) epsilondev_yz_crust_mantle

      write(IOUT) epsilondev_xx_inner_core
      write(IOUT) epsilondev_yy_inner_core
      write(IOUT) epsilondev_xy_inner_core
      write(IOUT) epsilondev_xz_inner_core
      write(IOUT) epsilondev_yz_inner_core

      ! rotation
      if (ROTATION_VAL) then
        write(IOUT) A_array_rotation
        write(IOUT) B_array_rotation
      endif

      ! attenuation memory variables
      if (ATTENUATION_VAL) then
        write(IOUT) R_xx_crust_mantle
        write(IOUT) R_yy_crust_mantle
        write(IOUT) R_xy_crust_mantle
        write(IOUT) R_xz_crust_mantle
        write(IOUT) R_yz_crust_mantle

        write(IOUT) R_xx_inner_core
        write(IOUT) R_yy_inner_core
        write(IOUT) R_xy_inner_core
        write(IOUT) R_xz_inner_core
        write(IOUT) R_yz_inner_core
      endif

      ! full gravity
      if (FULL_GRAVITY_VAL) then
        write(IOUT) neq
        write(IOUT) neq1
        ! Level-1 solver array
        write(IOUT) pgrav1
      endif

      close(IOUT)
    endif
  endif

  ! save last frame of the forward simulation
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  saving forward arrays'
      call flush_IMAIN()
    endif

    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call save_forward_arrays_adios()
    else
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
      outputname = trim(LOCAL_TMP_PATH)//'/'//trim(outputname)

      open(unit=IOUT,file=trim(outputname),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_forward_arrays** for writing')

      write(IOUT) displ_crust_mantle
      write(IOUT) veloc_crust_mantle
      write(IOUT) accel_crust_mantle

      write(IOUT) displ_inner_core
      write(IOUT) veloc_inner_core
      write(IOUT) accel_inner_core

      write(IOUT) displ_outer_core
      write(IOUT) veloc_outer_core
      write(IOUT) accel_outer_core

      write(IOUT) epsilondev_xx_crust_mantle
      write(IOUT) epsilondev_yy_crust_mantle
      write(IOUT) epsilondev_xy_crust_mantle
      write(IOUT) epsilondev_xz_crust_mantle
      write(IOUT) epsilondev_yz_crust_mantle

      write(IOUT) epsilondev_xx_inner_core
      write(IOUT) epsilondev_yy_inner_core
      write(IOUT) epsilondev_xy_inner_core
      write(IOUT) epsilondev_xz_inner_core
      write(IOUT) epsilondev_yz_inner_core

      if (ROTATION_VAL) then
        write(IOUT) A_array_rotation
        write(IOUT) B_array_rotation
      endif

      if (ATTENUATION_VAL) then
         write(IOUT) R_xx_crust_mantle
         write(IOUT) R_yy_crust_mantle
         write(IOUT) R_xy_crust_mantle
         write(IOUT) R_xz_crust_mantle
         write(IOUT) R_yz_crust_mantle

         write(IOUT) R_xx_inner_core
         write(IOUT) R_yy_inner_core
         write(IOUT) R_xy_inner_core
         write(IOUT) R_xz_inner_core
         write(IOUT) R_yz_inner_core
      endif

      ! full gravity
      if (FULL_GRAVITY_VAL) then
        write(IOUT) neq
        write(IOUT) neq1
        ! Level-1 solver array
        write(IOUT) pgrav1
      endif

      close(IOUT)
    endif
  endif

  end subroutine save_forward_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_forward_arrays_undoatt()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname

  ! transfers wavefields from GPU device to CPU host
  if (GPU_MODE) then
    call transfer_fields_cm_from_device(NDIM*NGLOB_CRUST_MANTLE, &
                                        displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle,Mesh_pointer)
    call transfer_fields_ic_from_device(NDIM*NGLOB_INNER_CORE, &
                                        displ_inner_core,veloc_inner_core,accel_inner_core,Mesh_pointer)
    call transfer_fields_oc_from_device(NGLOB_OUTER_CORE, &
                                        displ_outer_core,veloc_outer_core,accel_outer_core,Mesh_pointer)
    if (ROTATION_VAL) then
      call transfer_rotation_from_device(Mesh_pointer,A_array_rotation,B_array_rotation)
    endif

    if (ATTENUATION_VAL) then
      call transfer_rmemory_cm_from_device(Mesh_pointer,R_xx_crust_mantle,R_yy_crust_mantle,R_xy_crust_mantle, &
                                           R_xz_crust_mantle,R_yz_crust_mantle)
      call transfer_rmemory_ic_from_device(Mesh_pointer,R_xx_inner_core,R_yy_inner_core,R_xy_inner_core, &
                                           R_xz_inner_core,R_yz_inner_core)
    endif
  endif

  if (ADIOS_FOR_UNDO_ATTENUATION) then
    call save_forward_arrays_undoatt_adios()
  else
    ! current subset iteration
    iteration_on_subset_tmp = iteration_on_subset

    ! saves frame of the forward simulation

    write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
    outputname = trim(LOCAL_PATH)//'/'//trim(outputname)

    ! debug
    !if (myrank == 0 ) print *,'saving in: ',trim(LOCAL_PATH)//'/'//trim(outputname), iteration_on_subset_tmp,it

    open(unit=IOUT,file=trim(outputname),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for writing')

    write(IOUT) displ_crust_mantle
    write(IOUT) veloc_crust_mantle
    write(IOUT) accel_crust_mantle

    write(IOUT) displ_inner_core
    write(IOUT) veloc_inner_core
    write(IOUT) accel_inner_core

    write(IOUT) displ_outer_core
    write(IOUT) veloc_outer_core
    write(IOUT) accel_outer_core

    if (ROTATION_VAL) then
      write(IOUT) A_array_rotation
      write(IOUT) B_array_rotation
    endif

    if (ATTENUATION_VAL) then
      write(IOUT) R_xx_crust_mantle
      write(IOUT) R_yy_crust_mantle
      write(IOUT) R_xy_crust_mantle
      write(IOUT) R_xz_crust_mantle
      write(IOUT) R_yz_crust_mantle

      write(IOUT) R_xx_inner_core
      write(IOUT) R_yy_inner_core
      write(IOUT) R_xy_inner_core
      write(IOUT) R_xz_inner_core
      write(IOUT) R_yz_inner_core
    endif

    ! full gravity
    if (FULL_GRAVITY_VAL) then
      write(IOUT) neq
      write(IOUT) neq1
      ! Level-1 solver array
      write(IOUT) pgrav1
    endif

    close(IOUT)
  endif

  end subroutine save_forward_arrays_undoatt

!
!-------------------------------------------------------------------------------------------------
!


  subroutine save_forward_model_at_shifted_frequency(factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core)

! outputs model files in binary format

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,PI,GRAV,FOUR_THIRDS
  use shared_parameters, only: R_PLANET,RHOAV,LOCAL_PATH,TRANSVERSE_ISOTROPY

  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  real(kind=CUSTOM_REAL),dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL) :: factor_scale_relaxed_crust_mantle
  real(kind=CUSTOM_REAL),dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL) :: factor_scale_relaxed_inner_core

  ! local parameters
  real(kind=CUSTOM_REAL) :: scaleval1,scale_factor_r
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: muv_shifted,muh_shifted
  integer :: i,j,k,ispec,ier
  character(len=MAX_STRING_LEN) :: filename, prname

  ! debug
  logical, parameter :: OUTPUT_RELAXED_MODEL = .false.

  ! scaling factor to re-dimensionalize units
  scaleval1 = real( sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0), kind=CUSTOM_REAL)  ! velocities km/s

  ! note: since we only use shear attenuation, the shift occurs in muv values.
  !       thus, we output here only vpv, vsv or vp,vs for crust/mantle and inner core regions
  !       which are affected by the attenuation shift. all other model arrays stay the same.

  ! crust/mantle region
  if (NSPEC_CRUST_MANTLE > 0) then
    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

    ! uses temporary array
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
             muv_shifted(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
             muh_shifted(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), stat=ier)
    if (ier /= 0) stop 'Error allocating temp_store array'
    temp_store(:,:,:,:) = 0._CUSTOM_REAL

    ! safety check
    if (ANISOTROPIC_3D_MANTLE_VAL) &
      call exit_mpi(myrank,'ANISOTROPIC_3D_MANTLE not supported yet for shifted model file output')

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  shifted model files in directory: ',trim(LOCAL_PATH)
    endif

    ! user output
    if (myrank == 0) write(IMAIN,*) '  crust/mantle:'

    ! moduli (muv,muh) are at relaxed values (only Qmu implemented),
    ! scales back to have values at center frequency
    muv_shifted(:,:,:,:) = muvstore_crust_mantle(:,:,:,:)
    muh_shifted(:,:,:,:) = muhstore_crust_mantle(:,:,:,:)
    do ispec = 1,NSPEC_CRUST_MANTLE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              scale_factor_r = factor_scale_relaxed_crust_mantle(i,j,k,ispec)
            else
              scale_factor_r = factor_scale_relaxed_crust_mantle(1,1,1,ispec)
            endif
            ! scaling back from relaxed to values at shifted frequency
            ! (see in prepare_attenuation.f90 for how muv,muh are scaled to become relaxed moduli)
            ! muv
            muv_shifted(i,j,k,ispec) = muv_shifted(i,j,k,ispec) / scale_factor_r
            ! muh
            if (ispec_is_tiso_crust_mantle(ispec)) then
              muh_shifted(i,j,k,ispec) = muh_shifted(i,j,k,ispec) / scale_factor_r
            endif
          enddo
        enddo
      enddo
    enddo

    ! transverse isotropic model
    if (TRANSVERSE_ISOTROPY) then
      ! vpv (at relaxed values)
      temp_store(:,:,:,:) = sqrt((kappavstore_crust_mantle(:,:,:,:) &
                            + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                            * scaleval1

      filename = prname(1:len_trim(prname))//'vpv_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,temp_store,"shifted vpv")

      ! vph
      temp_store(:,:,:,:) = sqrt((kappahstore_crust_mantle(:,:,:,:) &
                            + FOUR_THIRDS * muh_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                            * scaleval1

      filename = prname(1:len_trim(prname))//'vph_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,temp_store,"shifted vph")

      ! vsv
      temp_store(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1

      filename = prname(1:len_trim(prname))//'vsv_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,temp_store,"shifted vsv")

      ! vsh
      temp_store(:,:,:,:) = sqrt( muh_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1

      filename = prname(1:len_trim(prname))//'vsh_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,temp_store,"shifted vsh")

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  written:'
        write(IMAIN,*) '    proc*_reg1_vpv_shifted.bin'
        write(IMAIN,*) '    proc*_reg1_vph_shifted.bin'
        write(IMAIN,*) '    proc*_reg1_vsv_shifted.bin'
        write(IMAIN,*) '    proc*_reg1_vsh_shifted.bin'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

    else
      ! isotropic model
      ! vp
      temp_store(:,:,:,:) = sqrt((kappavstore_crust_mantle(:,:,:,:) &
                            + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                            * scaleval1

      filename = prname(1:len_trim(prname))//'vp_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,temp_store,"shifted vp")

      ! vs
      temp_store(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1

      filename = prname(1:len_trim(prname))//'vs_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,temp_store,"shifted vs")

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  written:'
        write(IMAIN,*) '    proc*_reg1_vp_shifted.bin'
        write(IMAIN,*) '    proc*_reg1_vs_shifted.bin'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

    endif ! TRANSVERSE_ISOTROPY

    ! frees temporary array
    deallocate(temp_store,muv_shifted,muh_shifted)
  endif

  ! inner core region
  if (NSPEC_INNER_CORE > 0) then
    call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_PATH)

    ! uses temporary array
    allocate(temp_store(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
             muv_shifted(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), stat=ier)
    if (ier /= 0) stop 'Error allocating temp_store array'
    temp_store(:,:,:,:) = 0._CUSTOM_REAL

    ! user output
    if (myrank == 0) write(IMAIN,*) '  inner core:'

    ! moduli (muv,muh) are at relaxed values, scale back to have shifted values at center frequency
    muv_shifted(:,:,:,:) = muvstore_inner_core(:,:,:,:)
    do ispec = 1,NSPEC_INNER_CORE
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
              scale_factor_r = factor_scale_relaxed_inner_core(i,j,k,ispec)
            else
              scale_factor_r = factor_scale_relaxed_inner_core(1,1,1,ispec)
            endif

            ! inverts to scale relaxed back to shifted factor
            ! scaling back from relaxed to values at shifted frequency
            ! (see in prepare_attenuation.f90 for how muv,muh are scaled to become relaxed moduli)
            ! muv
            muv_shifted(i,j,k,ispec) = muv_shifted(i,j,k,ispec) / scale_factor_r
          enddo
        enddo
      enddo
    enddo

    if (ANISOTROPIC_INNER_CORE_VAL) then
      call exit_mpi(myrank,'ANISOTROPIC_INNER_CORE not supported yet for shifted model file output')
    else
      ! isotropic model
      ! vp
      temp_store(:,:,:,:) = sqrt((kappavstore_inner_core(:,:,:,:) &
                            + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_inner_core(:,:,:,:)) &
                            * scaleval1

      filename = prname(1:len_trim(prname))//'vp_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_INNER_CORE,temp_store,"shifted vp")

      ! vs
      temp_store(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_inner_core(:,:,:,:) )*scaleval1

      filename = prname(1:len_trim(prname))//'vs_shifted.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) temp_store
      close(IOUT)
      call print_gll_min_max_all(NSPEC_INNER_CORE,temp_store,"shifted vs")

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  written:'
        write(IMAIN,*) '    proc*_reg3_vp_shifted.bin'
        write(IMAIN,*) '    proc*_reg3_vs_shifted.bin'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

    endif ! TRANSVERSE_ISOTROPY

    ! frees temporary array
    deallocate(temp_store,muv_shifted)
  endif

  ! outputs relaxed model values
  if (OUTPUT_RELAXED_MODEL) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  outputting relaxed model:'
      call flush_IMAIN()
    endif
    ! checks
    if (.not. TRANSVERSE_ISOTROPY) stop 'Outputting relaxed model requires TRANSVERSE_ISOTROPY'

    ! scaling factor to re-dimensionalize units
    ! the scale of GPa--[g/cm^3][(km/s)^2]
    scaleval1 = real( ((sqrt(PI*GRAV*RHOAV)*R_PLANET/1000.d0)**2)*(RHOAV/1000.d0), kind=CUSTOM_REAL) ! moduli GPa

    ! crust/mantle region
    if (NSPEC_CRUST_MANTLE > 0) then
      call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)
      ! muv
      filename = prname(1:len_trim(prname))//'muv_relaxed.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) muvstore_crust_mantle * scaleval1
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,muvstore_crust_mantle * scaleval1,"relaxed muv")
      ! muh
      filename = prname(1:len_trim(prname))//'muh_relaxed.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) muhstore_crust_mantle * scaleval1
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,muhstore_crust_mantle * scaleval1,"relaxed muh")
      ! kappav
      filename = prname(1:len_trim(prname))//'kappav_relaxed.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) kappavstore_crust_mantle * scaleval1
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,kappavstore_crust_mantle * scaleval1,"relaxed kappav")
      ! kappah
      filename = prname(1:len_trim(prname))//'kappah_relaxed.bin'
      open(unit=IOUT,file=trim(filename),status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening file '// trim(filename))
      write(IOUT) kappahstore_crust_mantle * scaleval1
      close(IOUT)
      call print_gll_min_max_all(NSPEC_CRUST_MANTLE,kappahstore_crust_mantle * scaleval1,"relaxed kappah")

      ! user output
      if (myrank == 0) then
        write(IMAIN,*) '  written:'
        write(IMAIN,*) '    proc*_reg1_muv_relaxed.bin, ..'
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    endif
  endif

  end subroutine save_forward_model_at_shifted_frequency

