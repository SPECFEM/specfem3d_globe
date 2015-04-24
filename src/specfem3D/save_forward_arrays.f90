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

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) outputname, path_to_add

  ! checks if anything to do
  if (UNDO_ATTENUATION ) return

  ! save files to local disk or tape system if restart file
  if (NUMBER_OF_RUNS > 1 .and. NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS) then
    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call save_intermediate_forward_arrays_adios()
    else
      write(outputname,"('dump_all_arrays',i6.6)") myrank
      open(unit=IOUT,file=trim(LOCAL_TMP_PATH)//'/'//trim(outputname), &
          status='unknown',form='unformatted',action='write',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file dump_all_arrays*** for writing')

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

      write(IOUT) A_array_rotation
      write(IOUT) B_array_rotation

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

      close(IOUT)
    endif
  endif

  ! save last frame of the forward simulation
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call save_forward_arrays_adios()
    else
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
      outputname = trim(LOCAL_TMP_PATH)//'/'//trim(outputname)

      if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
        write(path_to_add,"('run',i4.4,'/')") mygroup + 1
        outputname = path_to_add(1:len_trim(path_to_add))//outputname(1:len_trim(outputname))
      endif

      open(unit=IOUT,file=trim(outputname),status='unknown', &
          form='unformatted',action='write',iostat=ier)
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

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  character(len=MAX_STRING_LEN) :: outputname, path_to_add

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

    if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
      write(path_to_add,"('run',i4.4,'/')") mygroup + 1
      outputname = path_to_add(1:len_trim(path_to_add))//outputname(1:len_trim(outputname))
    endif

    ! debug
    !if (myrank == 0 ) print*,'saving in: ',trim(LOCAL_PATH)//'/'//trim(outputname), iteration_on_subset_tmp,it

    open(unit=IOUT,file=trim(outputname), &
         status='unknown',form='unformatted',action='write',iostat=ier)
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

    close(IOUT)
  endif

  end subroutine save_forward_arrays_undoatt


