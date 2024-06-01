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

  subroutine read_forward_arrays_startrun()

! reads in saved wavefields

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_full_gravity

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) outputname
  ! full gravity
  integer :: neq_read,neq1_read

  ! checks if anything to do
  ! undoing attenuation doesn't support the following checkpointing
  if (UNDO_ATTENUATION) return

  ! read files back from local disk or MT tape system if restart file
  if (NUMBER_OF_THIS_RUN > 1) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  reading startup file for run ',NUMBER_OF_THIS_RUN
      call flush_IMAIN()
    endif

    if (ADIOS_FOR_FORWARD_ARRAYS) then
      call read_intermediate_forward_arrays_adios()
    else
      write(outputname,"('dump_all_arrays',i6.6)") myrank
      outputname = trim(LOCAL_TMP_PATH) // '/' // outputname(1:len_trim(outputname))

      open(unit=IIN,file=trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening file dump_all_arrays*** for reading')

      ! wavefield
      read(IIN) displ_crust_mantle
      read(IIN) veloc_crust_mantle
      read(IIN) accel_crust_mantle

      read(IIN) displ_inner_core
      read(IIN) veloc_inner_core
      read(IIN) accel_inner_core

      read(IIN) displ_outer_core
      read(IIN) veloc_outer_core
      read(IIN) accel_outer_core

      ! strain
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

      ! rotation
      if (ROTATION_VAL) then
        read(IIN) A_array_rotation
        read(IIN) B_array_rotation
      endif

      ! attenuation memory variables
      if (ATTENUATION_VAL) then
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
      endif

      ! full gravity
      if (FULL_GRAVITY_VAL) then
        read(IIN) neq_read
        read(IIN) neq1_read
        ! check if array sizes match
        if (neq_read /= neq) then
          print *,'Error reading forward array for startrun: rank ',myrank,'has read neq =',neq_read,' - shoud be ',neq
          call exit_MPI(myrank,'Invalid forward array neq for startrun')
        endif
        if (neq1_read /= neq1) then
          print *,'Error reading forward array for startrun: rank ',myrank,'has read neq1 =',neq1_read,' - shoud be ',neq1
          call exit_MPI(myrank,'Invalid forward array neq1 for startrun')
        endif
        read(IIN) pgrav1
      endif

      close(IIN)

      ! full gravity
      if (FULL_GRAVITY) then
        ! need to interpolate the gravity values
        call SIEM_interpolate_gravity()
      endif
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
  use specfem_par_full_gravity

  implicit none

  ! local parameters
  integer :: ier
  integer :: b_neq_read,b_neq1_read
  character(len=MAX_STRING_LEN) :: outputname

  ! checks if anything to do
  if (UNDO_ATTENUATION) return

  ! reads in file data
  if (ADIOS_FOR_FORWARD_ARRAYS) then
    call read_forward_arrays_adios()
  else
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_save_forward_arrays.bin'
    outputname = trim(LOCAL_TMP_PATH) // '/' // outputname(1:len_trim(outputname))

    open(unit=IIN,file=trim(outputname), &
          status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error: opening proc_****_save_forward_arrays.bin'
      print *,'path: ',outputname
      call exit_mpi(myrank,'Error open file save_forward_arrays.bin')
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

    ! full gravity
    if (FULL_GRAVITY_VAL) then
      read(IIN) b_neq_read
      read(IIN) b_neq1_read
      ! check if array sizes match
      if (b_neq_read /= neq) then
        print *,'Error reading forward array: rank ',myrank,'has read neq =',b_neq_read,' - shoud be ',neq
        call exit_MPI(myrank,'Invalid forward array neq')
      endif
      if (b_neq1_read /= neq1) then
        print *,'Error reading forward array: rank ',myrank,'has read neq1 =',b_neq1_read,' - shoud be ',neq1
        call exit_MPI(myrank,'Invalid forward array neq1')
      endif
      ! Level-1 solver array
      read(IIN) b_pgrav1
    endif

    close(IIN)
  endif ! ADIOS_FOR_FORWARD_ARRAYS

  ! full gravity
  if (FULL_GRAVITY) then
    ! need to interpolate the gravity values
    call SIEM_interpolate_backward_gravity()
  endif

  ! transfers fields onto GPU
  if (GPU_MODE) then
    ! transfers initialized wavefields to GPU device
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                        Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                        Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE,b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
                                        Mesh_pointer)
    ! strain
    if (.not. UNDO_ATTENUATION) then
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
  use specfem_par_full_gravity

  implicit none

  ! local parameters
  integer :: iteration_on_subset_tmp
  integer :: ier
  integer :: b_neq_read,b_neq1_read
  character(len=MAX_STRING_LEN) :: outputname

  ! current subset iteration
  iteration_on_subset_tmp = NSUBSET_ITERATIONS - iteration_on_subset + 1

  if (ADIOS_FOR_UNDO_ATTENUATION) then
    call read_forward_arrays_undoatt_adios(iteration_on_subset_tmp)
  else
    ! reads in saved wavefield
    write(outputname,'(a,i6.6,a,i6.6,a)') 'proc',myrank,'_save_frame_at',iteration_on_subset_tmp,'.bin'
    outputname = trim(LOCAL_PATH) // '/' // outputname(1:len_trim(outputname))

    ! debug
    !if (myrank == 0 ) print *,'reading in: ',trim(LOCAL_PATH)//'/'//trim(outputname),iteration_on_subset_tmp,iteration_on_subset,it

    ! opens corresponding snapshot file for reading
    open(unit=IIN,file=trim(outputname),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file proc***_save_frame_at** for reading')

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

    ! full gravity
    if (FULL_GRAVITY_VAL) then
      read(IIN) b_neq_read
      read(IIN) b_neq1_read
      ! check if array sizes match
      if (b_neq_read /= neq) then
        print *,'Error reading forward array: rank ',myrank,'has read neq =',b_neq_read,' - shoud be ',neq
        call exit_MPI(myrank,'Invalid forward array neq')
      endif
      if (b_neq1_read /= neq1) then
        print *,'Error reading forward array: rank ',myrank,'has read neq1 =',b_neq1_read,' - shoud be ',neq1
        call exit_MPI(myrank,'Invalid forward array neq1')
      endif
      ! Level-1 solver array
      read(IIN) b_pgrav1
    endif

    close(IIN)
  endif

  ! full gravity
  if (FULL_GRAVITY) then
    ! need to interpolate the gravity values
    call SIEM_interpolate_backward_gravity()
  endif

  ! transfers fields onto GPU
  if (GPU_MODE) then
    ! transfers initialized wavefields to GPU device
    call transfer_b_fields_cm_to_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle, &
                                        Mesh_pointer)

    call transfer_b_fields_ic_to_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core, &
                                        Mesh_pointer)

    call transfer_b_fields_oc_to_device(NGLOB_OUTER_CORE,b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core, &
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

