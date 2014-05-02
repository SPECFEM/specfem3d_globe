!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

  subroutine read_forward_arrays_startrun()

! reads in saved wavefields

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) outputname

  ! checks run/checkpoint number
  if(NUMBER_OF_RUNS < 1 .or. NUMBER_OF_RUNS > NSTEP) &
    stop 'number of restart runs can not be less than 1 or greater than NSTEP'
  if(NUMBER_OF_THIS_RUN < 1 .or. NUMBER_OF_THIS_RUN > NUMBER_OF_RUNS) &
    stop 'incorrect run number'
  if (SIMULATION_TYPE /= 1 .and. NUMBER_OF_RUNS /= 1) &
    stop 'Only 1 run for SIMULATION_TYPE = 2/3'

  ! define correct time steps if restart files
  ! set start/end steps for time iteration loop
  it_begin = (NUMBER_OF_THIS_RUN - 1) * (NSTEP / NUMBER_OF_RUNS) + 1
  if (NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS) then
    it_end = NUMBER_OF_THIS_RUN * (NSTEP / NUMBER_OF_RUNS)
  else
    ! Last run may be a bit larger
    it_end = NSTEP
  endif

  ! checks if anything to do
  ! undoing attenuation doesn't support the following checkpointing
  if( UNDO_ATTENUATION ) return

  ! read files back from local disk or MT tape system if restart file
  if(NUMBER_OF_THIS_RUN > 1) then
    if( ADIOS_ENABLED .and. ADIOS_FOR_FORWARD_ARRAYS ) then
      call read_intermediate_forward_arrays_adios()
    else
      write(outputname,"('dump_all_arrays',i6.6)") myrank
      open(unit=IIN,file=trim(LOCAL_TMP_PATH)//'/'//outputname,status='old',action='read',form='unformatted',iostat=ier)
      if( ier /= 0 ) call exit_MPI(myrank,'error opening file dump_all_arrays*** for reading')

      read(IIN) displ_crust_mantle
      read(IIN) veloc_crust_mantle
      read(IIN) accel_crust_mantle

      read(IIN) displ_inner_core
      read(IIN) veloc_inner_core
      read(IIN) accel_inner_core

      read(IIN) displ_outer_core
      read(IIN) veloc_outer_core
      read(IIN) accel_outer_core

      read(IIN) epsilondev_xx_crust_mantle
      read(IIN) epsilondev_yy_crust_mantle
      read(IIN) epsilondev_xy_crust_mantle
      read(IIN) epsilondev_xz_crust_mantle
      read(IIN) epsilondev_yz_crust_mantle

      read(IIN) epsilondev_xx_inner_core
      read(IIN) epsilondev_yy_inner_core
      read(IIN) epsilondev_xy_inner_core
      read(IIN) epsilondev_xz_inner_core
      read(IIN) epsilondev_yz_inner_core

      read(IIN) A_array_rotation
      read(IIN) B_array_rotation

      read(IIN) R_xx_crust_mantle
      read(IIN) R_yy_crust_mantle
      read(IIN) R_xy_crust_mantle
      read(IIN) R_xz_crust_mantle
      read(IIN) R_yz_crust_mantle

      read(IIN) R_xx_inner_core
      read(IIN) R_yy_inner_core
      read(IIN) R_xy_inner_core
      read(IIN) R_xz_inner_core
      read(IIN) R_yz_inner_core

      close(IIN)
    endif
  endif

  end subroutine read_forward_arrays_startrun

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays()

! reads in saved wavefields to reconstruct/backward wavefield

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) outputname

  ! checks if anything to do
  if( UNDO_ATTENUATION ) return

  ! reads in file data
  if( ADIOS_ENABLED .and. ADIOS_FOR_FORWARD_ARRAYS ) then
    call read_forward_arrays_adios()
  else
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
    open(unit=IIN,file=trim(LOCAL_TMP_PATH)//'/'//outputname, &
          status='old',action='read',form='unformatted',iostat=ier)
    if( ier /= 0 ) then
      print*,'error: opening proc_****_save_forward_arrays.bin'
      print*,'path: ',trim(LOCAL_TMP_PATH)//'/'//outputname
      call exit_mpi(myrank,'error open file save_forward_arrays.bin')
    endif

    read(IIN) b_displ_crust_mantle
    read(IIN) b_veloc_crust_mantle
    read(IIN) b_accel_crust_mantle

    read(IIN) b_displ_inner_core
    read(IIN) b_veloc_inner_core
    read(IIN) b_accel_inner_core

    read(IIN) b_displ_outer_core
    read(IIN) b_veloc_outer_core
    read(IIN) b_accel_outer_core

    read(IIN) b_epsilondev_xx_crust_mantle
    read(IIN) b_epsilondev_yy_crust_mantle
    read(IIN) b_epsilondev_xy_crust_mantle
    read(IIN) b_epsilondev_xz_crust_mantle
    read(IIN) b_epsilondev_yz_crust_mantle

    read(IIN) b_epsilondev_xx_inner_core
    read(IIN) b_epsilondev_yy_inner_core
    read(IIN) b_epsilondev_xy_inner_core
    read(IIN) b_epsilondev_xz_inner_core
    read(IIN) b_epsilondev_yz_inner_core

    if (ROTATION_VAL) then
      read(IIN) b_A_array_rotation
      read(IIN) b_B_array_rotation
    endif

    if (ATTENUATION_VAL) then
       read(IIN) b_R_xx_crust_mantle
       read(IIN) b_R_yy_crust_mantle
       read(IIN) b_R_xy_crust_mantle
       read(IIN) b_R_xz_crust_mantle
       read(IIN) b_R_yz_crust_mantle

       read(IIN) b_R_xx_inner_core
       read(IIN) b_R_yy_inner_core
       read(IIN) b_R_xy_inner_core
       read(IIN) b_R_xz_inner_core
       read(IIN) b_R_yz_inner_core
    endif
    close(IIN)
  endif ! ADIOS_FOR_FORWARD_ARRAYS

  ! transfers fields onto GPU
  if( GPU_MODE ) then
    ! transfers initialized wavefields to GPU device
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE, &
                                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                    Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE, &
                                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                    Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE, &
                                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                    Mesh_pointer)
    ! strain
    if( .not. UNDO_ATTENUATION ) then
      call transfer_b_strain_cm_to_device(Mesh_pointer, &
                                    b_epsilondev_xx_crust_mantle,b_epsilondev_yy_crust_mantle, &
                                    b_epsilondev_xy_crust_mantle,b_epsilondev_xz_crust_mantle, &
                                    b_epsilondev_yz_crust_mantle)

      call transfer_b_strain_ic_to_device(Mesh_pointer, &
                                    b_epsilondev_xx_inner_core,b_epsilondev_yy_inner_core, &
                                    b_epsilondev_xy_inner_core,b_epsilondev_xz_inner_core, &
                                    b_epsilondev_yz_inner_core)
    endif
    ! rotation
    if (ROTATION_VAL) then
      call transfer_b_rotation_to_device(Mesh_pointer,b_A_array_rotation,b_B_array_rotation)
    endif
    ! attenuation memory variables
    if (ATTENUATION_VAL) then
      call transfer_b_rmemory_cm_to_device(Mesh_pointer, &
                                           b_R_xx_crust_mantle,b_R_yy_crust_mantle, &
                                           b_R_xy_crust_mantle,b_R_xz_crust_mantle, &
                                           b_R_yz_crust_mantle)
      call transfer_b_rmemory_ic_to_device(Mesh_pointer, &
                                           b_R_xx_inner_core,b_R_yy_inner_core, &
                                           b_R_xy_inner_core,b_R_xz_inner_core, &
                                           b_R_yz_inner_core)
    endif
  endif

  end subroutine read_forward_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_forward_arrays_undoatt()

! reads in saved wavefields

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=150) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = NSTEP/NT_DUMP_ATTENUATION - iteration_on_subset + 1

  ! reads in saved wavefield
  write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'

  ! debug
  !if(myrank == 0 ) print*,'reading in: ',trim(LOCAL_PATH)//'/'//outputname, NSTEP/NT_DUMP_ATTENUATION,iteration_on_subset

  ! opens corresponding snapshot file for reading
  open(unit=IIN,file=trim(LOCAL_PATH)//'/'//outputname, &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error opening file proc***_save_frame_at** for reading')

  read(IIN) b_displ_crust_mantle
  read(IIN) b_veloc_crust_mantle
  read(IIN) b_accel_crust_mantle

  read(IIN) b_displ_inner_core
  read(IIN) b_veloc_inner_core
  read(IIN) b_accel_inner_core

  read(IIN) b_displ_outer_core
  read(IIN) b_veloc_outer_core
  read(IIN) b_accel_outer_core

  if (ROTATION_VAL) then
    read(IIN) b_A_array_rotation
    read(IIN) b_B_array_rotation
  endif

  if (ATTENUATION_VAL) then
    read(IIN) b_R_xx_crust_mantle
    read(IIN) b_R_yy_crust_mantle
    read(IIN) b_R_xy_crust_mantle
    read(IIN) b_R_xz_crust_mantle
    read(IIN) b_R_yz_crust_mantle

    read(IIN) b_R_xx_inner_core
    read(IIN) b_R_yy_inner_core
    read(IIN) b_R_xy_inner_core
    read(IIN) b_R_xz_inner_core
    read(IIN) b_R_yz_inner_core
  endif

  close(IIN)

  ! transfers fields onto GPU
  if( GPU_MODE ) then
    ! transfers initialized wavefields to GPU device
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE, &
                                    b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                    Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE, &
                                    b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                    Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE, &
                                    b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                    Mesh_pointer)
    ! rotation
    if (ROTATION_VAL) then
      call transfer_b_rotation_to_device(Mesh_pointer,b_A_array_rotation,b_B_array_rotation)
    endif

    ! attenuation memory variables
    if (ATTENUATION_VAL) then
      call transfer_b_rmemory_cm_to_device(Mesh_pointer, &
                                           b_R_xx_crust_mantle,b_R_yy_crust_mantle, &
                                           b_R_xy_crust_mantle,b_R_xz_crust_mantle, &
                                           b_R_yz_crust_mantle)
      call transfer_b_rmemory_ic_to_device(Mesh_pointer, &
                                           b_R_xx_inner_core,b_R_yy_inner_core, &
                                           b_R_xy_inner_core,b_R_xz_inner_core, &
                                           b_R_yz_inner_core)
    endif
  endif

  end subroutine read_forward_arrays_undoatt

