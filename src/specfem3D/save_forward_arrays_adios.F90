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

!-------------------------------------------------------------------------------
!> \file save_forward_arrays_adios.F90
!! \brief Save forward arrays with the help of the ADIOS library.
!! \author MPBL
!-------------------------------------------------------------------------------

#include "config.fh"


!-------------------------------------------------------------------------------
!> \brief Write intermediate forward arrays in an ADIOS file.
!!
!! This subroutine is only used when NUMBER_OF_RUNS > 1 and
!! NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS.

  subroutine save_intermediate_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name,group_name
  integer(kind=8) :: group_size_inc

  file_name = get_adios_filename(trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios")

  group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS_RESTART"
  call init_adios_group(myadios_fwd_group,group_name)

  ! defines ADIOS variables
  group_size_inc = 0
  call define_common_forward_arrays_adios(group_size_inc)
  call define_epsilon_forward_arrays_adios(group_size_inc)
  call define_rotation_forward_arrays_adios(group_size_inc)
  call define_attenuation_forward_arrays_adios(group_size_inc)

  ! Open an ADIOS handler to the restart file.
  call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
  call set_adios_group_size(myadios_fwd_file,group_size_inc)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios()
  call write_epsilon_forward_arrays_adios()
  call write_rotation_forward_arrays_adios()
  call write_attenuation_forward_arrays_adios()

  ! full gravity
  if (FULL_GRAVITY_VAL) then
    !TODO: implement full gravity adios save forward
    stop 'FULL_GRAVITY save_intermediate_forward_arrays_adios() not fully implemented yet'
    !write(IOUT) neq1
    !write(IOUT pgrav1
  endif

  ! Reset the path to its original value to avoid bugs.
  call write_adios_perform(myadios_fwd_file)
  ! flushes all engines (makes sure i/o is all written out)
  call flush_adios_group_all(myadios_fwd_group)
  ! Close ADIOS handler to the restart file.
  call close_file_adios(myadios_fwd_file)

  end subroutine save_intermediate_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name,group_name
  integer(kind=8) :: group_size_inc

  file_name = get_adios_filename(trim(LOCAL_TMP_PATH) // "/save_forward_arrays")

  group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"
  call init_adios_group(myadios_fwd_group,group_name)

  ! defines ADIOS variables
  group_size_inc = 0
  call define_common_forward_arrays_adios(group_size_inc)
  call define_epsilon_forward_arrays_adios(group_size_inc)
  if (ROTATION_VAL) call define_rotation_forward_arrays_adios(group_size_inc)
  if (ATTENUATION_VAL) call define_attenuation_forward_arrays_adios(group_size_inc)

  ! Open an ADIOS handler to the restart file.
  call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
  call set_adios_group_size(myadios_fwd_file,group_size_inc)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios()
  call write_epsilon_forward_arrays_adios()
  if (ROTATION_VAL) call write_rotation_forward_arrays_adios()
  if (ATTENUATION_VAL) call write_attenuation_forward_arrays_adios()

  ! full gravity
  if (FULL_GRAVITY_VAL) then
    !TODO: implement full gravity adios save forward
    stop 'FULL_GRAVITY save_forward_arrays_adios() not fully implemented yet'
    !write(IOUT) neq
    !write(IOUT) neq1
    !write(IOUT pgrav1
  endif

  ! Reset the path to its original value to avoid bugs.
  call write_adios_perform(myadios_fwd_file)
  ! flushes all engines (makes sure i/o is all written out)
  call flush_adios_group_all(myadios_fwd_group)
  ! Close ADIOS handler to the restart file.
  call close_file_adios(myadios_fwd_file)

  end subroutine save_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_undoatt_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! Local parameters
  integer :: iteration_on_subset_tmp
  character(len=MAX_STRING_LEN) :: file_name
  integer(kind=8),save :: group_size_inc
  ! ADIOS variables
  character(len=MAX_STRING_LEN) :: group_name
  ! multiple/single file for storage of snapshots
  logical :: do_open_file,do_close_file,do_init_group

  ! current subset iteration
  iteration_on_subset_tmp = iteration_on_subset

  ! file handling
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
    ! single file for all steps
    do_open_file = .false.
    do_close_file = .false.
    do_init_group = .false.

    ! single file
    file_name = get_adios_filename(trim(LOCAL_TMP_PATH) // "/save_forward_arrays_undoatt",ADIOS2_ENGINE_UNDO_ATT)

    group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT"

    ! open file at first call of this routine
    if (is_adios_version1) then
      ! adds steps by appending to file
      do_open_file = .true.
      do_close_file = .true.
      if (iteration_on_subset_tmp == 1) do_init_group = .true. ! only needs to initialize group once
    else
      ! adds steps by commands
      if (.not. is_initialized_fwd_group) then
        do_open_file = .true.
        do_init_group = .true.
      endif
    endif
  else
    ! for each step a single file
    do_open_file = .true.
    do_close_file = .true.
    do_init_group = .true.

    ! files for each iteration step
    write(file_name,'(a, a, i6.6)') trim(LOCAL_TMP_PATH), '/save_frame_at', iteration_on_subset_tmp
    file_name = get_adios_filename(trim(file_name))

    write(group_name, '(a, i6)') "SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT", iteration_on_subset_tmp
  endif

  ! debug
  !if (myrank == 0) print *,'debug: undoatt adios: save forward step iteration_on_subset_tmp = ',iteration_on_subset_tmp, &
  !                         'open/close file',do_open_file,do_close_file,do_init_group

  ! opens file for writing
  if (do_open_file) then
    ! prepares group & metadata
    !
    ! note, see adios manual:
    ! "These routines prepare ADIOS metadata construction,
    ! for example, setting up groups, variables, attributes and IO transport method,
    ! and hence must be called before any other ADIOS I/O operations,
    ! i.e., adios_open, adios_group_size, adios_write, adios_close."
    if (do_init_group) then
      call init_adios_group_undo_att(myadios_fwd_group,group_name)

      ! adds wavefield compression
      if (ADIOS_COMPRESSION_ALGORITHM /= 0) then
        ! sets adios flag to add compression operation for the following define_adios_** function calls
        call define_adios_compression()
      endif

      ! defines ADIOS variables
      group_size_inc = 0
      ! iteration number
      call define_adios_scalar(myadios_fwd_group, group_size_inc, '', "iteration", iteration_on_subset_tmp)

      ! wavefields (displ/veloc/accel) for all regions
      call define_common_forward_arrays_adios(group_size_inc)
      ! rotation arrays
      if (ROTATION_VAL) call define_rotation_forward_arrays_adios(group_size_inc)
      ! attenuation memory variables
      if (ATTENUATION_VAL) call define_attenuation_forward_arrays_adios(group_size_inc)

      ! re-sets compression flag (in case other routines will call the define_adios_** function calls)
      if (ADIOS_COMPRESSION_ALGORITHM /= 0) use_adios_compression = .false.
    endif

    ! Open an ADIOS handler to the restart file.
    if (is_adios_version1) then
      ! checks if we open for first time or append
      if (iteration_on_subset_tmp == 1) then
        ! creates new file
        call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
      else
        ! append to existing file
        call open_file_adios_write_append(myadios_fwd_file,myadios_fwd_group,file_name,group_name)

        ! debug: note, do not call as the inquiry on the appended file handle will seg-fault
        ! call show_adios_file_variables(myadios_fwd_file,myadios_fwd_group,file_name)
      endif

      ! debug
      !if (myrank == 0) print *,'debug: undoatt adios: save forward step = ',iteration_on_subset_tmp,' handle ',myadios_fwd_file

    else
      ! version 2, only opens once at beginning
      call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
    endif

    call set_adios_group_size(myadios_fwd_file,group_size_inc)

    ! sets flag
    is_initialized_fwd_group = .true.
  endif

  ! indicate new step section
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) call write_adios_begin_step(myadios_fwd_file)

  ! iteration number
  call write_adios_scalar(myadios_fwd_file, myadios_fwd_group, "iteration", iteration_on_subset_tmp)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios()
  if (ROTATION_VAL) call write_rotation_forward_arrays_adios()
  if (ATTENUATION_VAL) call write_attenuation_forward_arrays_adios()

  ! full gravity
  if (FULL_GRAVITY_VAL) then
    !TODO: implement full gravity adios save forward
    stop 'FULL_GRAVITY save_forward_arrays_undoatt_adios() not fully implemented yet'
    !write(IOUT) neq
    !write(IOUT) neq1
    !write(IOUT pgrav1
  endif

  ! perform writing
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
    ! end step to indicate output is completed. ADIOS2 can do I/O
    call write_adios_end_step(myadios_fwd_file)
  else
    ! Reset the path to its original value to avoid bugs and write out arrays.
    call write_adios_perform(myadios_fwd_file)
  endif

  ! Close ADIOS handler to the restart file.
  if (do_close_file) then
    ! flushes all engines (makes sure i/o is all written out)
    call flush_adios_group_all(myadios_fwd_group)
    ! closes file
    call close_file_adios(myadios_fwd_file)
    ! re-sets flag
    is_initialized_fwd_group = .false.

    ! debug
    !if (myrank == 0) print *,'debug: undoatt adios: close file save forward step = ',iteration_on_subset_tmp, &
    !                         ' handle ',myadios_fwd_file
  endif

  end subroutine save_forward_arrays_undoatt_adios


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are always dumped.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_common_forward_arrays_adios(group_size_inc)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer(kind=8), intent(inout) :: group_size_inc

  ! local parameters
  integer(kind=8) :: local_dim

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(displ_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(veloc_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(accel_crust_mantle))

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(displ_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(veloc_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(accel_inner_core))

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(displ_outer_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(veloc_outer_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(accel_outer_core))

  end subroutine define_common_forward_arrays_adios


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are always dumped
!! except for undo attenuation.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_epsilon_forward_arrays_adios(group_size_inc)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer(kind=8), intent(inout) :: group_size_inc
  ! local parameters
  integer(kind=8) :: local_dim

  ! strains
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_xx_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_yy_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_xy_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_xz_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_yz_crust_mantle))


  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_xx_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_yy_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_xy_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_xz_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', &
                                   STRINGIFY_VAR(epsilondev_yz_inner_core))

  end subroutine define_epsilon_forward_arrays_adios


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are dumped if ROTATION is true.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_rotation_forward_arrays_adios(group_size_inc)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer(kind=8), intent(inout) :: group_size_inc

  ! local parameters
  integer(kind=8) :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(A_array_rotation))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(B_array_rotation))

  end subroutine define_rotation_forward_arrays_adios

!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are dumped if ATTENUATION is true.
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_attenuation_forward_arrays_adios(group_size_inc)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer(kind=8), intent(inout) :: group_size_inc

  ! local parameters
  integer(kind=8) :: local_dim

  ! attenuation memory variables
  ! crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_xx_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_yy_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_xy_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_xz_crust_mantle))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_yz_crust_mantle))

  ! inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_xx_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_yy_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_xy_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_xz_inner_core))
  call define_adios_global_array1D(myadios_fwd_group, group_size_inc, local_dim, '', STRINGIFY_VAR(R_yz_inner_core))

  end subroutine define_attenuation_forward_arrays_adios


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.

  subroutine write_common_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer(kind=8) :: local_dim

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(displ_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(veloc_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(accel_crust_mantle))

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(displ_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(veloc_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(accel_inner_core))

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(displ_outer_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(veloc_outer_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(accel_outer_core))

  end subroutine write_common_forward_arrays_adios


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.

  subroutine write_epsilon_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer(kind=8) :: local_dim

  ! strains
  ! crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT

  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xx_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_yy_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xy_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xz_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_yz_crust_mantle))

  ! inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT

  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xx_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_yy_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xy_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_xz_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(epsilondev_yz_inner_core))

  end subroutine write_epsilon_forward_arrays_adios


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ROTATION is true.

  subroutine write_rotation_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer(kind=8) :: local_dim

  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(A_array_rotation))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(B_array_rotation))

  end subroutine write_rotation_forward_arrays_adios

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ATTENUATION
!!  is true.

  subroutine write_attenuation_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  ! local parameters
  integer(kind=8) :: local_dim

  ! attenuation memory variables
  ! crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION

  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_xx_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_yy_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_xy_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_xz_crust_mantle))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_yz_crust_mantle))

  ! inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION

  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_xx_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_yy_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_xy_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_xz_inner_core))
  call write_adios_global_1d_array(myadios_fwd_file, myadios_fwd_group, myrank, sizeprocs_adios, local_dim, &
                                   STRINGIFY_VAR(R_yz_inner_core))

  end subroutine write_attenuation_forward_arrays_adios


!
!-------------------------------------------------------------------------------------------------
!


  subroutine save_forward_model_at_shifted_frequency_adios(factor_scale_relaxed_crust_mantle,factor_scale_relaxed_inner_core)

! outputs model files in binary format

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,LOCAL_PATH,TRANSVERSE_ISOTROPY

  use specfem_par_crustmantle
  use specfem_par_innercore

  use adios_helpers_mod
  use manager_adios

  implicit none

  real(kind=CUSTOM_REAL),dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL) :: factor_scale_relaxed_crust_mantle
  real(kind=CUSTOM_REAL),dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL) :: factor_scale_relaxed_inner_core

  ! local parameters
  integer :: ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2,scale_factor_r
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_vpv,temp_store_vph,temp_store_vsv,temp_store_vsh
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_vp_ic,temp_store_vs_ic
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store_rho,temp_store_rho_ic
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: muv_shifted,muh_shifted
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: dummy_ijke_crust_mantle
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: dummy_ijke_inner_core
  integer :: iregion_code,nspec,nglob,i,j,k,ispec
  character(len=MAX_STRING_LEN) :: outputname, group_name
  integer(kind=8) :: local_dim
  integer(kind=8) :: group_size_inc
  ! ADIOS variables
  character(len=128) :: region_name, region_name_scalar

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '    shifted model    in ADIOS 1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '    shifted model    in ADIOS 2 file format'
#endif
    call flush_IMAIN()
  endif

  ! scaling factors to re-dimensionalize units
  scaleval1 = real( sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0), kind=CUSTOM_REAL)  ! velocities
  scaleval2 = real( RHOAV/1000.d0, kind=CUSTOM_REAL)                            ! densities

  ! dummy for definitions
  allocate(dummy_ijke_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           dummy_ijke_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating dummy arrays'
  dummy_ijke_crust_mantle(:,:,:,:) = 0.0
  dummy_ijke_inner_core(:,:,:,:) = 0.0

  ! note: since we only use shear attenuation, the shift occurs in muv values.
  !       thus, we output here only vpv, vsv or vp,vs for crust/mantle and inner core regions
  !       which are affected by the attenuation shift. all other model arrays stay the same.

  ! initializes i/o group
  group_size_inc = 0
  group_name = 'SPECFEM3D_GLOBE_MODEL_SHIFTED'
  call init_adios_group(myadios_val_group,group_name)

  ! crust/mantle region
  if (NSPEC_CRUST_MANTLE > 0) then
    iregion_code = IREGION_CRUST_MANTLE
    nspec = NSPEC_CRUST_MANTLE
    nglob = NGLOB_CRUST_MANTLE

    ! region name
    write(region_name,"('reg',i1, '/')") iregion_code
    write(region_name_scalar,"('reg',i1)") iregion_code

    ! save nspec and nglob, to be used in combine_paraview_data
    call define_adios_scalar (myadios_val_group, group_size_inc, &
                              region_name_scalar, STRINGIFY_VAR(nspec))
    call define_adios_scalar (myadios_val_group, group_size_inc, &
                              region_name_scalar, STRINGIFY_VAR(nglob))

    ! array sizes
    local_dim = NGLLX * NGLLY * NGLLZ * nspec

    ! checks size
    if (size(kappavstore_crust_mantle) /= local_dim) then
      print *,'Error: size kappavstore ',size(kappavstore_crust_mantle), ' should be ',local_dim
      call exit_mpi(myrank,'Error size kappavstore_crust_mantle for storing meshfiles')
    endif

    ! safety check
    if (ANISOTROPIC_3D_MANTLE_VAL) &
      call exit_mpi(myrank,'ANISOTROPIC_3D_MANTLE not supported yet for shifted model file output')

    !--- Define ADIOS variables -----------------------------
    if (TRANSVERSE_ISOTROPY) then
      ! transverse isotropic model
      ! vpv
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vpv", dummy_ijke_crust_mantle)
      ! vph
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vph", dummy_ijke_crust_mantle)
      ! vsv
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vsv", dummy_ijke_crust_mantle)
      ! vsv
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vsh", dummy_ijke_crust_mantle)
    else
      ! isotropic model
      ! vp
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vp", dummy_ijke_crust_mantle)
      ! vs (will store it even for the outer core, although it should just be zero there)
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vs", dummy_ijke_crust_mantle)
    endif
    ! rho
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "rho", dummy_ijke_crust_mantle)
  endif

  ! inner core
  if (NSPEC_INNER_CORE > 0) then
    iregion_code = IREGION_INNER_CORE
    nspec = NSPEC_INNER_CORE
    nglob = NGLOB_INNER_CORE

    ! region name
    write(region_name,"('reg',i1, '/')") iregion_code
    write(region_name_scalar,"('reg',i1)") iregion_code

    ! save nspec and nglob, to be used in combine_paraview_data
    call define_adios_scalar (myadios_val_group, group_size_inc, &
                              region_name_scalar, STRINGIFY_VAR(nspec))
    call define_adios_scalar (myadios_val_group, group_size_inc, &
                              region_name_scalar, STRINGIFY_VAR(nglob))

    ! array sizes
    local_dim = NGLLX * NGLLY * NGLLZ * nspec

    ! checks size
    if (size(kappavstore_inner_core) /= local_dim) then
      print *,'Error: size kappavstore ',size(kappavstore_inner_core), ' should be ',local_dim
      call exit_mpi(myrank,'Error size kappavstore_inner_core for storing meshfiles')
    endif

    !--- Define ADIOS variables -----------------------------
    if (ANISOTROPIC_INNER_CORE_VAL) then
      call exit_mpi(myrank,'ANISOTROPIC_INNER_CORE not supported yet for shifted model file output')
    else
      ! isotropic model
      ! vp
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vp", dummy_ijke_inner_core)
      ! vs
      call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "vs", dummy_ijke_inner_core)
    endif
    ! rho
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, region_name, "rho", dummy_ijke_inner_core)

    deallocate(dummy_ijke_crust_mantle,dummy_ijke_inner_core)
  endif

  !--- Open an ADIOS handler to the restart file. ---------
  outputname = get_adios_filename(trim(LOCAL_PATH) // "/model_gll_shifted")

  ! user output
  if (myrank == 0) write(IMAIN,*) '    saving shifted model arrays in ADIOS file: ',trim(outputname)

  ! opens file for writing
  call open_file_adios_write(myadios_val_file,myadios_val_group,outputname,group_name)

  call set_adios_group_size(myadios_val_file,group_size_inc)

  ! crust/mantle region
  if (NSPEC_CRUST_MANTLE > 0) then
    iregion_code = IREGION_CRUST_MANTLE
    nspec = NSPEC_CRUST_MANTLE
    nglob = NGLOB_CRUST_MANTLE

    ! region name
    write(region_name,"('reg',i1, '/')") iregion_code
    write(region_name_scalar,"('reg',i1)") iregion_code

    ! save nspec and nglob, to be used in combine_paraview_data
    call write_adios_scalar(myadios_val_file,myadios_val_group,trim(region_name) // "nspec",nspec)
    call write_adios_scalar(myadios_val_file,myadios_val_group,trim(region_name) // "nglob",nglob)

  ! note: the following uses temporary arrays for array expressions like sqrt( (kappavstore+..)).
  !       since the write_adios_** calls might be in deferred mode, these temporary arrays should be valid
  !       until a perform/close/end_step call is done.
  !
  !       as a work-around, we will explicitly allocate temporary arrays and deallocate them after the file close.
    allocate(temp_store_vpv(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vph(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vsv(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vsh(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_rho(NGLLX,NGLLY,NGLLZ,nspec), &
             muv_shifted(NGLLX,NGLLY,NGLLZ,nspec), &
             muh_shifted(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating temp vp,.. arrays'
    temp_store_vpv(:,:,:,:) = 0.0_CUSTOM_REAL
    temp_store_vph(:,:,:,:) = 0.0_CUSTOM_REAL
    temp_store_vsv(:,:,:,:) = 0.0_CUSTOM_REAL
    temp_store_vsh(:,:,:,:) = 0.0_CUSTOM_REAL
    temp_store_rho(:,:,:,:) = 0.0_CUSTOM_REAL

    ! moduli (muv,muh) are at relaxed values (only Qmu implemented),
    ! scales back to have values at center frequency
    muv_shifted(:,:,:,:) = muvstore_crust_mantle(:,:,:,:)
    muh_shifted(:,:,:,:) = muhstore_crust_mantle(:,:,:,:)
    do ispec = 1,nspec
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
      ! vpv
      temp_store_vpv(:,:,:,:) = sqrt((kappavstore_crust_mantle(:,:,:,:) &
                            + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                            * scaleval1

      call print_gll_min_max_all(nspec,temp_store_vpv,"shifted vpv")

      ! vph
      temp_store_vph(:,:,:,:) = sqrt((kappahstore_crust_mantle(:,:,:,:) &
                            + FOUR_THIRDS * muh_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                            * scaleval1

      call print_gll_min_max_all(nspec,temp_store_vph,"shifted vph")

      ! vsv
      temp_store_vsv(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1

      call print_gll_min_max_all(nspec,temp_store_vsv,"shifted vsv")

      ! vsh
      temp_store_vsh(:,:,:,:) = sqrt( muh_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1

      call print_gll_min_max_all(nspec,temp_store_vsh,"shifted vsh")

      !--- Schedule writes for the previously defined ADIOS variables
      ! vpv
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vpv", temp_store_vpv)
      ! vph
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vph", temp_store_vph)
      ! vsv
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vsv", temp_store_vsv)
      ! vsh
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vsh", temp_store_vsh)
    else
      ! isotropic model
      ! vp
      temp_store_vpv(:,:,:,:) = sqrt((kappavstore_crust_mantle(:,:,:,:) &
                            + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_crust_mantle(:,:,:,:)) &
                            * scaleval1

      call print_gll_min_max_all(nspec,temp_store_vpv,"shifted vp")

      ! vs
      temp_store_vsv(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_crust_mantle(:,:,:,:) )*scaleval1

      call print_gll_min_max_all(nspec,temp_store_vsv,"shifted vs")

      !--- Schedule writes for the previously defined ADIOS variables
      ! vp
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vp", temp_store_vpv)
      ! vs
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vs", temp_store_vsv)
    endif ! TRANSVERSE_ISOTROPY
    ! rho
    temp_store_rho(:,:,:,:) = rhostore_crust_mantle(:,:,:,:) *scaleval2
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(region_name) // "rho", temp_store_rho)

    deallocate(muv_shifted,muh_shifted)
  endif

  ! inner core
  if (NSPEC_INNER_CORE > 0) then
    iregion_code = IREGION_INNER_CORE
    nspec = NSPEC_INNER_CORE
    nglob = NGLOB_INNER_CORE

    ! region name
    write(region_name,"('reg',i1, '/')") iregion_code
    write(region_name_scalar,"('reg',i1)") iregion_code

    ! uses temporary array
    allocate(temp_store_vp_ic(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_vs_ic(NGLLX,NGLLY,NGLLZ,nspec), &
             temp_store_rho_ic(NGLLX,NGLLY,NGLLZ,nspec), &
             muv_shifted(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating temp_store array'
    temp_store_vp_ic(:,:,:,:) = 0.0_CUSTOM_REAL
    temp_store_vs_ic(:,:,:,:) = 0.0_CUSTOM_REAL
    temp_store_rho_ic(:,:,:,:) = 0.0_CUSTOM_REAL

    ! moduli (muv,muh) are at relaxed values, scale back to have shifted values at center frequency
    muv_shifted(:,:,:,:) = muvstore_inner_core(:,:,:,:)
    do ispec = 1,nspec
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

    ! isotropic model
    if (ANISOTROPIC_INNER_CORE_VAL) then
      call exit_mpi(myrank,'ANISOTROPIC_INNER_CORE not supported yet for shifted model file output')
    else
      ! isotropic model
      ! vp
      temp_store_vp_ic(:,:,:,:) = sqrt((kappavstore_inner_core(:,:,:,:) &
                            + FOUR_THIRDS * muv_shifted(:,:,:,:))/rhostore_inner_core(:,:,:,:)) &
                            * scaleval1

      call print_gll_min_max_all(nspec,temp_store_vp_ic,"shifted vp")

      ! vs
      temp_store_vs_ic(:,:,:,:) = sqrt( muv_shifted(:,:,:,:)/rhostore_inner_core(:,:,:,:) )*scaleval1

      call print_gll_min_max_all(nspec,temp_store_vs_ic,"shifted vs")

      !--- Schedule writes for the previously defined ADIOS variables
      ! vp
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vp", temp_store_vp_ic)
      ! vs
      call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                       trim(region_name) // "vs", temp_store_vs_ic)
    endif ! TRANSVERSE_ISOTROPY
    ! rho
    temp_store_rho_ic(:,:,:,:) = rhostore_inner_core(:,:,:,:) *scaleval2
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                   trim(region_name) // "rho", temp_store_rho_ic)
  endif

  !--- Reset the path to zero and perform the actual write to disk
  call write_adios_perform(myadios_val_file)
  ! flushes all engines related to this io group (not really required, but used to make sure i/o has all written out)
  call flush_adios_group_all(myadios_val_group)
  ! closes file
  call close_file_adios(myadios_val_file)

  ! adios should be done with writing memory out.
  call synchronize_all()

  ! frees temporary array
  if (allocated(temp_store_vpv)) deallocate(temp_store_vpv,temp_store_vph,temp_store_vsv,temp_store_vsh)
  if (allocated(temp_store_vp_ic)) deallocate(temp_store_vp_ic,temp_store_vs_ic)
  if (allocated(muv_shifted)) deallocate(muv_shifted)

  end subroutine save_forward_model_at_shifted_frequency_adios
