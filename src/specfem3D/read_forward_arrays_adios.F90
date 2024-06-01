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
!> \file read_forward_arrays_adios.F90
!! \brief Read saved forward arrays with the help of the ADIOS library.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios()
!!       and read_forward_arrays_adios() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays

  subroutine read_intermediate_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none
  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name,group_name
  integer(kind=8) :: local_dim
  ! ADIOS variables
  integer(kind=8) :: sel
  integer(kind=8), dimension(1) :: start, count

  file_name = get_adios_filename(trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios")

  group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS_RESTART"
  call init_adios_group(myadios_group,group_name)

  ! opens adios file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

! note: we only use one sel pointer/variable here.
!       calling adios_set_selection will internally allocate each time a selection structure which at the end of reading
!       should be deleted by adios_delete_selection.
!       however, here we re-assign the sel variable and probably produce some garbage when calling set_selection again
!       with the same variable. furthermore, scheduling the reads and then re-assigning the sel variable leads to
!       segmentation faults, if the read hasn't been performed yet. thus, we call adios_perform_read before re-assigning
!       the set_selection_boundingbox() routine.
!
!       in future, a cleaner way might be found. we could store all different sel structures in an array and then
!       free the pointers/structures at the very end. this could also allow to call perform_read at the end.

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ_crust_mantle/array", displ_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc_crust_mantle/array", veloc_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel_crust_mantle/array", accel_crust_mantle)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ_inner_core/array", displ_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc_inner_core/array", veloc_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel_inner_core/array", accel_inner_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ_outer_core/array", displ_outer_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc_outer_core/array", veloc_outer_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel_outer_core/array", accel_outer_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! strains crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xx_crust_mantle/array", epsilondev_xx_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yy_crust_mantle/array", epsilondev_yy_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xy_crust_mantle/array", epsilondev_xy_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xz_crust_mantle/array", epsilondev_xz_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yz_crust_mantle/array", epsilondev_yz_crust_mantle)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! strains inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xx_inner_core/array", epsilondev_xx_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yy_inner_core/array", epsilondev_yy_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xy_inner_core/array", epsilondev_xy_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xz_inner_core/array", epsilondev_xz_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yz_inner_core/array", epsilondev_yz_inner_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! rotation
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "A_array_rotation/array", A_array_rotation)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "B_array_rotation/array", B_array_rotation)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! attenuation memory variables crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xx_crust_mantle/array", R_xx_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yy_crust_mantle/array", R_yy_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xy_crust_mantle/array", R_xy_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xz_crust_mantle/array", R_xz_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yz_crust_mantle/array", R_yz_crust_mantle)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! attenuation memory variables inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xx_inner_core/array", R_xx_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yy_inner_core/array", R_yy_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xy_inner_core/array", R_xy_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xz_inner_core/array", R_xz_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yz_inner_core/array", R_yz_inner_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! full gravity
  if (FULL_GRAVITY_VAL) then
    !TODO: implement full gravity adios read forward
    stop 'FULL_GRAVITY read_intermediate_forward_arrays_adios() not fully implemented yet'
    !read(IIN) neq1
    !allocate(pgrav1_oldrun(0:neq1))
    !read(IIN) pgrav1_oldrun
  endif

  ! closes ADIOS handler to the restart file.
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,group_name)

  end subroutine read_intermediate_forward_arrays_adios

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios() and read_forward_arrays_adios() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays

  subroutine read_forward_arrays_adios()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none
  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name,group_name
  ! ADIOS variables
  integer(kind=8) :: local_dim
  integer(kind=8) :: sel
  integer(kind=8), dimension(1) :: start, count

  file_name = get_adios_filename(trim(LOCAL_TMP_PATH) // "/save_forward_arrays")

  group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"
  call init_adios_group(myadios_group,group_name)

  ! opens adios file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  ! reads in arrays
  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ_crust_mantle/array", b_displ_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc_crust_mantle/array", b_veloc_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel_crust_mantle/array", b_accel_crust_mantle)

  ! NOTE: perform reads before changing selection, otherwise it will segfault
  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ_inner_core/array", b_displ_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc_inner_core/array", b_veloc_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel_inner_core/array", b_accel_inner_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "displ_outer_core/array", b_displ_outer_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "veloc_outer_core/array", b_veloc_outer_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "accel_outer_core/array", b_accel_outer_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! strains crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xx_crust_mantle/array",b_epsilondev_xx_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yy_crust_mantle/array",b_epsilondev_yy_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xy_crust_mantle/array",b_epsilondev_xy_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xz_crust_mantle/array",b_epsilondev_xz_crust_mantle)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yz_crust_mantle/array",b_epsilondev_yz_crust_mantle)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! strains inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xx_inner_core/array", b_epsilondev_xx_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yy_inner_core/array", b_epsilondev_yy_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xy_inner_core/array", b_epsilondev_xy_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_xz_inner_core/array", b_epsilondev_xz_inner_core)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "epsilondev_yz_inner_core/array", b_epsilondev_yz_inner_core)

  call read_adios_perform(myadios_file)
  call delete_adios_selection(sel)

  ! rotation
  if (ROTATION_VAL) then
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "A_array_rotation/array", b_A_array_rotation)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "B_array_rotation/array", b_B_array_rotation)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  ! attenuation memory variables
  if (ATTENUATION_VAL) then
    ! crust/mantle
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xx_crust_mantle/array", b_R_xx_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yy_crust_mantle/array", b_R_yy_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xy_crust_mantle/array", b_R_xy_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xz_crust_mantle/array", b_R_xz_crust_mantle)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yz_crust_mantle/array", b_R_yz_crust_mantle)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)

    ! inner core
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xx_inner_core/array", b_R_xx_inner_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yy_inner_core/array", b_R_yy_inner_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xy_inner_core/array", b_R_xy_inner_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_xz_inner_core/array", b_R_xz_inner_core)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, "R_yz_inner_core/array", b_R_yz_inner_core)

    call read_adios_perform(myadios_file)
    call delete_adios_selection(sel)
  endif

  if (FULL_GRAVITY_VAL) then
    ! Read in the gravity
    !TODO: implement full gravity adios read forward
    stop 'FULL_GRAVITY read_forward_arrays_adios() not fully implemented yet'
    !read(IIN) b_neq_read
    !read(IIN) b_neq1_read
    !read(IIN) b_pgrav1
  endif

  ! closes ADIOS handler to the restart file.
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,group_name)

  end subroutine read_forward_arrays_adios


!-------------------------------------------------------------------------------
!> \brief Read forward arrays for undo attenuation from an ADIOS file.

  subroutine read_forward_arrays_undoatt_adios(iteration_on_subset_tmp)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios_helpers_mod
  use manager_adios

  implicit none
  ! Arguments
  integer, intent(in) :: iteration_on_subset_tmp

  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name,group_name
  ! ADIOS variables
  integer(kind=8) :: local_dim
  integer(kind=8), dimension(1) :: start, count
  ! shorten the name of iteration variable and make it integer*8
  integer(kind=8) :: step
  integer :: t_tmp
  ! selection
  !integer(kind=8)         :: sel
  ! selections array
  integer(kind=8), dimension(8),target :: selections
  integer :: sel_num, i
  integer(kind=8), pointer :: sel => null()
  ! multiple/single file for storage of snapshots
  logical :: do_open_file,do_close_file

  ! selections array
  sel_num = 0

  ! file handling
  if (ADIOS_SAVE_ALL_SNAPSHOTS_IN_ONE_FILE) then
    ! iterations here go down from N to 1, but ADIOS files has steps 0..N-1
    step = iteration_on_subset_tmp - 1

    ! single file for all steps
    do_open_file = .false.
    do_close_file = .false.

    ! single file
    file_name = get_adios_filename(trim(LOCAL_TMP_PATH) // "/save_forward_arrays_undoatt",ADIOS2_ENGINE_UNDO_ATT)

    group_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT"

    ! open file at first call of this routine
    if (.not. is_initialized_fwd_group) do_open_file = .true.
    ! close at last step (step counting down from N-1 to 0)
    if (step == 0) do_close_file = .true.

  else
    ! single ADIOS files for each subset with step 0 entry
    step = 0

    ! for each step a single file
    do_open_file = .true.
    do_close_file = .true.

    ! files for each iteration step
    write(file_name,'(a,a,i6.6)') trim(LOCAL_TMP_PATH), '/save_frame_at',iteration_on_subset_tmp
    file_name = get_adios_filename(trim(file_name))

    write(group_name, '(a, i6)') "SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT", iteration_on_subset_tmp
  endif

  ! debug
  !if (myrank == 0) print *,'debug: undoatt adios: read forward step = ',step,'open/close file',do_open_file,do_close_file

  ! opens file for reading
  if (do_open_file) then
    ! creates adios group
    call init_adios_group_undo_att(myadios_fwd_group,group_name)

    ! opens adios file
    call open_file_adios_read(myadios_fwd_file,myadios_fwd_group,file_name)

    ! debug
    !if (myrank == 0) call show_adios_file_variables(myadios_fwd_file,myadios_fwd_group,file_name)

    ! sets flag
    is_initialized_fwd_group = .true.
  endif

  ! daniel todo: see if calling read_adios_perform() only once at the end would increase performance.
  !              this will need to have sel variables stored in an array.
  !
  !              note also that reading scalars is always blocking/synchronized.
  !
  ! checks if correct snapshot number of wavefields
  call read_adios_scalar(myadios_fwd_file,myadios_fwd_group,myrank,"iteration",t_tmp,step)
  if (t_tmp /= iteration_on_subset_tmp) then
    print *,'Error: invalid iteration step found in reading undoatt arrays: found ',t_tmp,' instead of ',iteration_on_subset_tmp
    call exit_mpi(myrank,'Invalid iteration step read in read_forward_arrays_undoatt_adios() routine')
  endif

  ! reads in arrays
  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel,start,count)

  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "displ_crust_mantle/array", b_displ_crust_mantle,step)
  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "veloc_crust_mantle/array", b_veloc_crust_mantle,step)
  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "accel_crust_mantle/array", b_accel_crust_mantle,step)
  ! NOTE: perform reads before changing selection, otherwise it will segfault if reusing the sel variable
  !call read_adios_perform(myadios_fwd_file)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "displ_inner_core/array", b_displ_inner_core,step)
  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "veloc_inner_core/array", b_veloc_inner_core,step)
  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "accel_inner_core/array", b_accel_inner_core,step)
  !call read_adios_perform(myadios_fwd_file)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "displ_outer_core/array", b_displ_outer_core,step)
  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "veloc_outer_core/array", b_veloc_outer_core,step)
  call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                 "accel_outer_core/array", b_accel_outer_core,step)
  !call read_adios_perform(myadios_fwd_file)

  ! rotation
  if (ROTATION_VAL) then
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "A_array_rotation/array", b_A_array_rotation,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "B_array_rotation/array", b_B_array_rotation,step)
    !call read_adios_perform(myadios_fwd_file)
  endif

  ! attenuation memory variables
  if (ATTENUATION_VAL) then
    ! crust/mantle
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_xx_crust_mantle/array", b_R_xx_crust_mantle,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_yy_crust_mantle/array", b_R_yy_crust_mantle,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_xy_crust_mantle/array", b_R_xy_crust_mantle,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_xz_crust_mantle/array", b_R_xz_crust_mantle,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_yz_crust_mantle/array", b_R_yz_crust_mantle,step)
    !call read_adios_perform(myadios_fwd_file)

    ! inner core
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
    start(1) = local_dim * int(myrank,kind=8); count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_xx_inner_core/array", b_R_xx_inner_core,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_yy_inner_core/array", b_R_yy_inner_core,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_xy_inner_core/array", b_R_xy_inner_core,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_xz_inner_core/array", b_R_xz_inner_core,step)
    call read_adios_schedule_array(myadios_fwd_file, myadios_fwd_group, sel, start, count, &
                                   "R_yz_inner_core/array", b_R_yz_inner_core,step)
    !call read_adios_perform(myadios_fwd_file)
  endif

  ! perform actual reading
  call read_adios_perform(myadios_fwd_file)

  ! free selection structures
  do i = 1, sel_num
    sel => selections(i)
    call delete_adios_selection(sel)
  enddo

  if (FULL_GRAVITY_VAL) then
    ! Read in the gravity
    !TODO: implement full gravity adios read forward
    stop 'FULL_GRAVITY read_forward_arrays_undoatt_adios() not fully implemented yet'
    !read(IIN) b_neq_read
    !read(IIN) b_neq1_read
    !read(IIN) b_pgrav1
  endif

  ! Close ADIOS handler to the restart file.
  ! note: for single file, only close at the very end
  if (do_close_file) then
    call close_file_adios_read_and_finalize_method(myadios_fwd_file)
    call delete_adios_group(myadios_fwd_group,group_name)

    is_initialized_fwd_group = .false.

    ! debug
    !if (myrank == 0) print *,'debug: undoatt adios: close file read forward step = ',step,' handle ',myadios_fwd_file
  endif

  end subroutine read_forward_arrays_undoatt_adios
