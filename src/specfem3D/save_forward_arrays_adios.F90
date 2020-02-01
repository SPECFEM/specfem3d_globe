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

  file_name = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"

  group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"
  call init_adios_group(myadios_fwd_group,group_name)

  ! Define ADIOS variables
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

  ! Reset the path to its original value to avoid bugs.
  call write_adios_perform(myadios_fwd_file)
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

  file_name = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"

  group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"
  call init_adios_group(myadios_fwd_group,group_name)

  ! Define ADIOS variables
  group_size_inc = 0
  call define_common_forward_arrays_adios(group_size_inc)
  call define_epsilon_forward_arrays_adios(group_size_inc)
  ! TODO check following:
  ! conditional definition of vars seem to mess with the group size,
  ! even if the variables are conditionally written.
  !  if (ROTATION_VAL) then
    call define_rotation_forward_arrays_adios(group_size_inc)
  !  endif
  !  if (ATTENUATION_VAL) then
    call define_attenuation_forward_arrays_adios(group_size_inc)
  !  endif

  ! Open an ADIOS handler to the restart file.
  call open_file_adios_write(myadios_fwd_file,myadios_fwd_group,file_name,group_name)
  call set_adios_group_size(myadios_fwd_file,group_size_inc)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios()
  call write_epsilon_forward_arrays_adios()

  if (ROTATION_VAL) then
    call write_rotation_forward_arrays_adios()
  endif

  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios()
  endif

  ! Reset the path to its original value to avoid bugs.
  call write_adios_perform(myadios_fwd_file)
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
    file_name = trim(LOCAL_TMP_PATH) // "/save_forward_arrays_undoatt.bp"
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
    write(file_name,'(a, a, i6.6, a)') trim(LOCAL_TMP_PATH), '/save_frame_at', iteration_on_subset_tmp,'.bp'
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

      ! Define ADIOS variables
      group_size_inc = 0

      ! iteration number
      call define_adios_scalar(myadios_fwd_group, group_size_inc, '', "iteration", iteration_on_subset_tmp)

      call define_common_forward_arrays_adios(group_size_inc)
      ! TODO check following:
      ! conditional definition of vars seem to mess with the group size,
      ! even if the variables are conditionally written.
      if (ROTATION_VAL) then
        call define_rotation_forward_arrays_adios(group_size_inc)
      endif
      if (ATTENUATION_VAL) then
        call define_attenuation_forward_arrays_adios(group_size_inc)
      endif
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

  if (ROTATION_VAL) then
    call write_rotation_forward_arrays_adios()
  endif

  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios()
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
    ! > call flush_adios_group_all(myadios_fwd_group)
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

