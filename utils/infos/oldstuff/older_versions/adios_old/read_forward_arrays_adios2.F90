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
!> \file read_forward_arrays_adios.F90
!! \brief Read saved forward arrays with the help of the ADIOS library.
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios2()
!!       and read_forward_arrays_adios2() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays

module forward_adios2_read
  use adios2, only: adios2_engine, adios2_io, adios2_variable, adios2_attribute
  use adios2_helpers_read_mod
  type(adios2_io), public :: io_fwdatt
  type(adios2_variable), public :: v_x
  character(len=*), parameter :: src="read_forward_arrays_adios2"

end module forward_adios2_read

  subroutine read_intermediate_forward_arrays_adios2()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_read

  implicit none
  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name, io_name
  character(len=*), parameter :: func="read_intermediate_forward_arrays_adios2()"
  integer :: local_dim
  ! ADIOS variables
  type(adios2_io)         :: io
  type(adios2_engine)     :: file
  integer                 :: ier
  integer(kind=8), dimension(1) :: start, count

  file_name = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"
  io_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"

  ! Create the ADIOS IO group which will contain all variables and attributes
  call adios2_declare_io(io, adios2obj, io_name, ier)
  call check_adios2_err(myrank,ier,src,func, "Declare IO group")

  ! Set engine and parameters
  call adios2_set_engine(io, ADIOS2_ENGINE_DEFAULT, ier)

  ! Open the handle to file containing all the ADIOS variables
  call adios2_open(file, io, file_name, adios2_mode_read, ier)
  call check_adios2_err(myrank,ier,src,func, "Open file "//file_name)

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim

  call read_adios2(file, io, "displ_crust_mantle", 1, start, count, displ_crust_mantle, myrank, src, func)
  ! the above call is equivalent to:
  ! call adios2_inquire_variable(v, io, "displ_crust_mantle", ier)
  ! call check_adios2_err(myrank,ier,src,func, "Inquire variable displ_crust_mantle")
  ! call adios2_set_selection(v, 1, start, count, ier)
  ! call adios2_get(file, v, displ_crust_mantle, ier)
  call read_adios2(file, io, "veloc_crust_mantle", 1, start, count, veloc_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "accel_crust_mantle", 1, start, count, accel_crust_mantle, myrank, src, func)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "displ_inner_core", 1, start, count, displ_inner_core, myrank, src, func)
  call read_adios2(file, io, "veloc_inner_core", 1, start, count, veloc_inner_core, myrank, src, func)
  call read_adios2(file, io, "accel_inner_core", 1, start, count, accel_inner_core, myrank, src, func)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "displ_outer_core", 1, start, count, displ_outer_core, myrank, src, func)
  call read_adios2(file, io, "veloc_outer_core", 1, start, count, veloc_outer_core, myrank, src, func)
  call read_adios2(file, io, "accel_outer_core", 1, start, count, accel_outer_core, myrank, src, func)

  ! strains crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "epsilondev_xx_crust_mantle", 1, start, count, epsilondev_xx_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yy_crust_mantle", 1, start, count, epsilondev_yy_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xy_crust_mantle", 1, start, count, epsilondev_xy_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xz_crust_mantle", 1, start, count, epsilondev_xz_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yz_crust_mantle", 1, start, count, epsilondev_yz_crust_mantle, myrank, src, func)

  ! strains inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "epsilondev_xx_inner_core", 1, start, count, epsilondev_xx_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yy_inner_core", 1, start, count, epsilondev_yy_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xy_inner_core", 1, start, count, epsilondev_xy_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xz_inner_core", 1, start, count, epsilondev_xz_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yz_inner_core", 1, start, count, epsilondev_yz_inner_core, myrank, src, func)

  ! rotation
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "A_array_rotation", 1, start, count, A_array_rotation, myrank, src, func)
  call read_adios2(file, io, "B_array_rotation", 1, start, count, B_array_rotation, myrank, src, func)

  ! attenuation memory variables crust/mantle
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "R_xx_crust_mantle", 1, start, count, R_xx_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "R_yy_crust_mantle", 1, start, count, R_yy_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "R_xy_crust_mantle", 1, start, count, R_xy_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "R_xz_crust_mantle", 1, start, count, R_xz_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "R_yz_crust_mantle", 1, start, count, R_yz_crust_mantle, myrank, src, func)

  ! attenuation memory variables inner core
  local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "R_xx_inner_core", 1, start, count, R_xx_inner_core, myrank, src, func)
  call read_adios2(file, io, "R_yy_inner_core", 1, start, count, R_yy_inner_core, myrank, src, func)
  call read_adios2(file, io, "R_xy_inner_core", 1, start, count, R_xy_inner_core, myrank, src, func)
  call read_adios2(file, io, "R_xz_inner_core", 1, start, count, R_xz_inner_core, myrank, src, func)
  call read_adios2(file, io, "R_yz_inner_core", 1, start, count, R_yz_inner_core, myrank, src, func)

  ! closes adios file
  call adios2_close(file, ier)
  call check_adios2_err(myrank,ier,src,func, "Close file "//file_name)

  end subroutine read_intermediate_forward_arrays_adios2

!-------------------------------------------------------------------------------
!> \brief Read forward arrays from an ADIOS file.
!> \note read_intermediate_forward_arrays_adios2() and read_forward_arrays_adios2() are not factorized, because
!>       the latest read the bp file in "b_" prefixed arrays

  subroutine read_forward_arrays_adios2()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_read

  implicit none
  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name, io_name
  character(len=*), parameter :: func="read_forward_arrays_adios2()"
  integer :: local_dim
  ! ADIOS variables
  type(adios2_io)         :: io
  type(adios2_engine)     :: file
  integer                 :: ier
  integer(kind=8), dimension(1) :: start, count

  file_name = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"
  io_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"

  ! Create the ADIOS IO group which will contain all variables and attributes
  call adios2_declare_io(io, adios2obj, io_name, ier)
  call check_adios2_err(myrank,ier,src,func, "Declare IO group")

  ! Set engine and parameters
  call adios2_set_engine(io, ADIOS2_ENGINE_DEFAULT, ier)

  ! Open the handle to file containing all the ADIOS variables
  call adios2_open(file, io, file_name, adios2_mode_read, ier)
  call check_adios2_err(myrank,ier,src,func, "Open file "//file_name)

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim

  call read_adios2(file, io, "displ_crust_mantle", 1, start, count, displ_crust_mantle, myrank, src, func)
  ! the above call is equivalent to:
  ! call adios2_inquire_variable(v, io, "displ_crust_mantle", ier)
  ! call check_adios2_err(myrank,ier,src,func, "Inquire variable displ_crust_mantle")
  ! call adios2_set_selection(v, 1, start, count, ier)
  ! call adios2_get(file, v, displ_crust_mantle, ier)
  call read_adios2(file, io, "veloc_crust_mantle", 1, start, count, veloc_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "accel_crust_mantle", 1, start, count, accel_crust_mantle, myrank, src, func)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "displ_inner_core", 1, start, count, displ_inner_core, myrank, src, func)
  call read_adios2(file, io, "veloc_inner_core", 1, start, count, veloc_inner_core, myrank, src, func)
  call read_adios2(file, io, "accel_inner_core", 1, start, count, accel_inner_core, myrank, src, func)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "displ_outer_core", 1, start, count, displ_outer_core, myrank, src, func)
  call read_adios2(file, io, "veloc_outer_core", 1, start, count, veloc_outer_core, myrank, src, func)
  call read_adios2(file, io, "accel_outer_core", 1, start, count, accel_outer_core, myrank, src, func)

  ! strains crust/mantle
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "epsilondev_xx_crust_mantle", 1, start, count, epsilondev_xx_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yy_crust_mantle", 1, start, count, epsilondev_yy_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xy_crust_mantle", 1, start, count, epsilondev_xy_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xz_crust_mantle", 1, start, count, epsilondev_xz_crust_mantle, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yz_crust_mantle", 1, start, count, epsilondev_yz_crust_mantle, myrank, src, func)

  ! strains inner core
  local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io, "epsilondev_xx_inner_core", 1, start, count, epsilondev_xx_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yy_inner_core", 1, start, count, epsilondev_yy_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xy_inner_core", 1, start, count, epsilondev_xy_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_xz_inner_core", 1, start, count, epsilondev_xz_inner_core, myrank, src, func)
  call read_adios2(file, io, "epsilondev_yz_inner_core", 1, start, count, epsilondev_yz_inner_core, myrank, src, func)

  ! rotation
  if (ROTATION_VAL) then
     local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
     start(1) = local_dim*myrank; count(1) = local_dim
     call read_adios2(file, io, "A_array_rotation", 1, start, count, A_array_rotation, myrank, src, func)
     call read_adios2(file, io, "B_array_rotation", 1, start, count, B_array_rotation, myrank, src, func)
  endif


  if (ATTENUATION_VAL) then
     ! attenuation memory variables crust/mantle
     local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
     start(1) = local_dim*myrank; count(1) = local_dim
     call read_adios2(file, io, "R_xx_crust_mantle", 1, start, count, R_xx_crust_mantle, myrank, src, func)
     call read_adios2(file, io, "R_yy_crust_mantle", 1, start, count, R_yy_crust_mantle, myrank, src, func)
     call read_adios2(file, io, "R_xy_crust_mantle", 1, start, count, R_xy_crust_mantle, myrank, src, func)
     call read_adios2(file, io, "R_xz_crust_mantle", 1, start, count, R_xz_crust_mantle, myrank, src, func)
     call read_adios2(file, io, "R_yz_crust_mantle", 1, start, count, R_yz_crust_mantle, myrank, src, func)

     ! attenuation memory variables inner core
     local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
     start(1) = local_dim*myrank; count(1) = local_dim
     call read_adios2(file, io, "R_xx_inner_core", 1, start, count, R_xx_inner_core, myrank, src, func)
     call read_adios2(file, io, "R_yy_inner_core", 1, start, count, R_yy_inner_core, myrank, src, func)
     call read_adios2(file, io, "R_xy_inner_core", 1, start, count, R_xy_inner_core, myrank, src, func)
     call read_adios2(file, io, "R_xz_inner_core", 1, start, count, R_xz_inner_core, myrank, src, func)
     call read_adios2(file, io, "R_yz_inner_core", 1, start, count, R_yz_inner_core, myrank, src, func)
  endif

  ! closes adios file
  call adios2_close(file, ier)
  call check_adios2_err(myrank,ier,src,func, "Close file "//file_name)

  end subroutine read_forward_arrays_adios2


!-------------------------------------------------------------------------------
!> \brief Read forward arrays for undo attenuation from an ADIOS file.

  subroutine read_forward_arrays_undoatt_adios2(iteration_on_subset_tmp)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_read

  implicit none
  ! Arguments
  integer, intent(in) :: iteration_on_subset_tmp
  ! Local parameters
  character(len=MAX_STRING_LEN) :: file_name, io_name
  character(len=*), parameter :: func="read_forward_arrays_undoatt_adios2()"
  integer :: local_dim
  ! ADIOS variables
  type(adios2_engine)     :: file
  integer                 :: ier
  integer(kind=8), dimension(1) :: start, count
  ! shorten the name of iteration variable and make it integer*8
  integer(kind=8) :: t

  !! iterations here go down from N to 1
  !! but ADIOS files has steps 0..N-1
  t = iteration_on_subset_tmp - 1

  file_name = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"
  io_name = "SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT"

  if (.not. io_fwdatt%valid) then
    ! Create the ADIOS IO group which will contain all variables and attributes
    call adios2_declare_io(io_fwdatt, adios2obj, io_name, ier)
    call check_adios2_err(myrank,ier,src,func, "Declare IO group")

    ! Set engine and parameters
    call adios2_set_engine(io_fwdatt, ADIOS2_ENGINE_DEFAULT, ier)
  endif

  if (.not. adios2_file_fwdatt%valid) then
    ! Open the handle to file containing all the ADIOS variables
    call adios2_open(adios2_file_fwdatt, io_fwdatt, file_name, adios2_mode_read, ier)
    call check_adios2_err(myrank,ier,src,func, "Open file "//file_name)
  endif

  file = adios2_file_fwdatt
  !
  !call adios2_set_step_selection(variable, iteration_on_subset_tmp, 1, ier)
  !
  call check_adios2_err(myrank,ier,src,func, "Seeking step in "//file_name)

  ! crust/mantle
  local_dim = NDIM * NGLOB_CRUST_MANTLE
  start(1) = local_dim*myrank; count(1) = local_dim

  call read_adios2(file, io_fwdatt, "displ_crust_mantle", 1, start, count, displ_crust_mantle, myrank, src,func, t)
  ! the above call is equivalent to:
  ! call adios2_inquire_variable(v, io, "displ_crust_mantle", ier)
  ! call check_adios2_err(myrank,ier,src,func, "Inquire variable displ_crust_mantle")
  ! call adios2_set_selection(v, 1, start, count, ier)
  ! call adios2_get(file, v, displ_crust_mantle, ier)
  call read_adios2(file, io_fwdatt, "veloc_crust_mantle", 1, start, count, veloc_crust_mantle, myrank, src,func, t)
  call read_adios2(file, io_fwdatt, "accel_crust_mantle", 1, start, count, accel_crust_mantle, myrank, src,func, t)

  ! inner core
  local_dim = NDIM * NGLOB_INNER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io_fwdatt, "displ_inner_core", 1, start, count, displ_inner_core, myrank, src,func, t)
  call read_adios2(file, io_fwdatt, "veloc_inner_core", 1, start, count, veloc_inner_core, myrank, src,func, t)
  call read_adios2(file, io_fwdatt, "accel_inner_core", 1, start, count, accel_inner_core, myrank, src,func, t)

  ! outer core
  local_dim = NGLOB_OUTER_CORE
  start(1) = local_dim*myrank; count(1) = local_dim
  call read_adios2(file, io_fwdatt, "displ_outer_core", 1, start, count, displ_outer_core, myrank, src,func, t)
  call read_adios2(file, io_fwdatt, "veloc_outer_core", 1, start, count, veloc_outer_core, myrank, src,func, t)
  call read_adios2(file, io_fwdatt, "accel_outer_core", 1, start, count, accel_outer_core, myrank, src,func, t)

  ! rotation
  if (ROTATION_VAL) then
    local_dim = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call read_adios2(file, io_fwdatt, "A_array_rotation", 1, start, count, A_array_rotation, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "B_array_rotation", 1, start, count, B_array_rotation, myrank, src,func, t)
  endif


  if (ATTENUATION_VAL) then
    ! attenuation memory variables crust/mantle
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call read_adios2(file, io_fwdatt, "R_xx_crust_mantle", 1, start, count, R_xx_crust_mantle, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_yy_crust_mantle", 1, start, count, R_yy_crust_mantle, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_xy_crust_mantle", 1, start, count, R_xy_crust_mantle, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_xz_crust_mantle", 1, start, count, R_xz_crust_mantle, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_yz_crust_mantle", 1, start, count, R_yz_crust_mantle, myrank, src,func, t)

    ! attenuation memory variables inner core
    local_dim = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
    start(1) = local_dim*myrank; count(1) = local_dim
    call read_adios2(file, io_fwdatt, "R_xx_inner_core", 1, start, count, R_xx_inner_core, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_yy_inner_core", 1, start, count, R_yy_inner_core, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_xy_inner_core", 1, start, count, R_xy_inner_core, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_xz_inner_core", 1, start, count, R_xz_inner_core, myrank, src,func, t)
    call read_adios2(file, io_fwdatt, "R_yz_inner_core", 1, start, count, R_yz_inner_core, myrank, src,func, t)
  endif

  ! Never close the file here. It will be closed in adios2_manager::finalize_adios2()

  end subroutine read_forward_arrays_undoatt_adios2
