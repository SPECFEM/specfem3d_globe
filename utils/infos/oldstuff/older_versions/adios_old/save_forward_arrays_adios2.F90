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


module forward_adios2_write
  use adios2, only: adios2_engine, adios2_io, adios2_variable, adios2_attribute
  type(adios2_engine), public :: file_fwd, file_fwdim
  type(adios2_io), public :: io_fwd, io_fwdim, io_fwdatt

  type(adios2_variable), public :: v_x

end module forward_adios2_write
!-------------------------------------------------------------------------------
!> \brief Write intermediate forward arrays in an ADIOS file.
!!
!! This subroutine is only used when NUMBER_OF_RUNS > 1 and
!! NUMBER_OF_THIS_RUN < NUMBER_OF_RUNS.

  subroutine save_intermediate_forward_arrays_adios2()

  use constants, only: MAX_STRING_LEN, ADIOS2_ENGINE_DEFAULT
  use shared_parameters, only: LOCAL_TMP_PATH

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  ! Local parameters
  character(len=MAX_STRING_LEN) :: outputname, ioname
  integer :: ier

  outputname = trim(LOCAL_TMP_PATH) // "/dump_all_arrays_adios.bp"
  ioname = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"

  ! Create the ADIOS IO group which will contain all variables and attributes
  call adios2_declare_io(io_fwdim, adios2obj, ioname, ier)
  if (ier /= 0) then
    print *, &
    'Error declaring an ADIOS2 IO group for itermediate forward arrays output in save_intermediate_forward_arrays_adios2()'
    stop 'Error declaring an ADIOS2 IO group: calling save_intermediate_forward_arrays_adios2() routine failed'
  endif

  ! Set engine and parameters
  call adios2_set_engine(io_fwdim, ADIOS2_ENGINE_DEFAULT, ier)

  ! Open the handle to file containing all the ADIOS variables
  call adios2_open(file_fwdim, io_fwdim, outputname, adios2_mode_write, ier)
  if (ier /= 0) then
    print *,'Error opening ADIOS2 file for itermediate forward arrays in save_intermediate_forward_arrays_adios2()'
    stop 'Error opening ADIOS2 file: calling save_intermediate_forward_arrays_adios2() routine failed'
  endif

  ! Define ADIOS variables
  call define_common_forward_arrays_adios2(io_fwdim)
  call define_epsilon_forward_arrays_adios2(io_fwdim)
  call define_rotation_forward_arrays_adios2(io_fwdim)
  call define_attenuation_forward_arrays_adios2(io_fwdim)

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios2(file_fwdim)
  call write_epsilon_forward_arrays_adios2(file_fwdim)
  call write_rotation_forward_arrays_adios2(file_fwdim)
  call write_attenuation_forward_arrays_adios2(file_fwdim)

  ! Close ADIOS handler to the restart file.
  call adios2_close(file_fwdim, ier)

  end subroutine save_intermediate_forward_arrays_adios2

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios2() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_adios2()

  use constants, only: MAX_STRING_LEN, ADIOS2_ENGINE_DEFAULT
  use constants_solver, only: ATTENUATION_VAL, ROTATION_VAL
  use shared_parameters, only: LOCAL_TMP_PATH

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  ! Local parameters
  character(len=MAX_STRING_LEN) :: outputname, ioname
  integer :: ier

  outputname = trim(LOCAL_TMP_PATH) // "/save_forward_arrays.bp"
  ioname = "SPECFEM3D_GLOBE_FORWARD_ARRAYS"

  ! Create the ADIOS IO group which will contain all variables and attributes
  call adios2_declare_io(io_fwd, adios2obj, ioname, ier)
  if (ier /= 0) then
    print *,'Error declaring an ADIOS2 IO group for forward arrays output in save_forward_arrays_adios2()'
    stop 'Error declaring an ADIOS2 IO group: calling save_forward_arrays_adios2() routine failed'
  endif

  ! Set engine and parameters
  call adios2_set_engine(io_fwd, ADIOS2_ENGINE_DEFAULT, ier)

  ! Open the handle to file containing all the ADIOS variables
  call adios2_open(file_fwd, io_fwd, outputname, adios2_mode_write, ier)
  if (ier /= 0) then
    print *,'Error opening ADIOS2 file for forward arrays in save_forward_arrays_adios2()'
    stop 'Error opening ADIOS2 file: calling save_forward_arrays_adios2() routine failed'
  endif

  ! Define ADIOS variables
  call define_common_forward_arrays_adios2(io_fwd)
  call define_epsilon_forward_arrays_adios2(io_fwd)
  if (ROTATION_VAL) then
    call define_rotation_forward_arrays_adios2(io_fwd)
  endif
  if (ATTENUATION_VAL) then
    call define_attenuation_forward_arrays_adios2(io_fwd)
  endif

  ! Issue the order to write the previously defined variable to the ADIOS file
  call write_common_forward_arrays_adios2(file_fwd)

  call write_epsilon_forward_arrays_adios2(file_fwd)

  if (ROTATION_VAL) then
      call write_rotation_forward_arrays_adios2(file_fwd)
  endif

  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios2(file_fwd)
  endif

  ! Close ADIOS handler to the restart file.
  call adios2_close(file_fwd, ier)

  end subroutine save_forward_arrays_adios2

!-------------------------------------------------------------------------------
!> \brief Write selected forward arrays in an ADIOS file.
!!
!! This subroutine is only used for forward simulations when
!! SAVE_FORWARD is set to .true. It dumps the same arrays than
!! save_intermediate_forward_arrays_adios2() except than some arrays
!! are only dumped if ROTATION and ATTENUATION are set to .true.

  subroutine save_forward_arrays_undoatt_adios2()

  use constants, only: MAX_STRING_LEN, ADIOS2_ENGINE_UNDO_ATT, ADIOS2_ENGINE_PARAMS_UNDO_ATT
  use constants_solver, only: ATTENUATION_VAL, ROTATION_VAL
  use shared_parameters, only: LOCAL_PATH
  use specfem_par, only: iteration_on_subset

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  ! Local parameters
  character(len=MAX_STRING_LEN) :: outputname, ioname
  integer :: ier
  type(adios2_variable), save :: v_iter

  ! current subset iteration
  !write(outputname,'(a, a, i6.6, a)') trim(LOCAL_PATH), '/save_frame_at', iteration_on_subset,'.bp'
  outputname = trim(LOCAL_PATH) // "/save_forward_arrays.bp"
  ioname = "SPECFEM3D_GLOBE_FORWARD_ARRAYS_UNDOATT"

  if (.not. io_fwdatt%valid) then
     ! Create the ADIOS IO group which will contain all variables and attributes
     call adios2_declare_io(io_fwdatt, adios2obj, ioname, ier)
     if (ier /= 0) then
       print *,'Error declaring an ADIOS2 IO group for forward arrays output in save_forward_arrays_undoatt_adios2()'
       stop 'Error declaring an ADIOS2 IO group: calling save_forward_arrays_undoatt_adios2() routine failed'
     endif

     ! Set engine and parameters
     call adios2_set_engine(io_fwdatt, ADIOS2_ENGINE_UNDO_ATT, ier)
     ! Set parameters to ADIOS2_ENGINE_PARAMS_UNDO_ATT
     call adios2_set_parameters(io_fwdatt, ADIOS2_ENGINE_PARAMS_UNDO_ATT, ier)

     ! Define ADIOS variables
     call adios2_define_variable(v_iter, io_fwdatt, "iteration", adios2_type_integer4, ier)
     call define_common_forward_arrays_adios2(io_fwdatt)
     if (ROTATION_VAL) then
        call define_rotation_forward_arrays_adios2(io_fwdatt)
     endif
     if (ATTENUATION_VAL) then
        call define_attenuation_forward_arrays_adios2(io_fwdatt)
     endif
  endif

  if (.not. adios2_file_fwdatt%valid) then
     ! Open the handle to file containing all the ADIOS variables
     call adios2_open(adios2_file_fwdatt, io_fwdatt, outputname, adios2_mode_write, ier)
     if (ier /= 0) then
       print *,'Error opening ADIOS2 file for undoatt forward arrays in save_forward_arrays_undoatt_adios2()'
       stop 'Error opening ADIOS2 file: calling save_forward_arrays_undoatt_adios2() routine failed'
     endif
  endif

  call adios2_begin_step(adios2_file_fwdatt, ier)
  ! Issue the order to write the previously defined variable to the ADIOS file
  call adios2_put(adios2_file_fwdatt, v_iter, iteration_on_subset, ier)
  call write_common_forward_arrays_adios2(adios2_file_fwdatt)

  if (ROTATION_VAL) then
    call write_rotation_forward_arrays_adios2(adios2_file_fwdatt)
  endif

  if (ATTENUATION_VAL) then
    call write_attenuation_forward_arrays_adios2(adios2_file_fwdatt)
  endif

  ! end step to indicate output is completed. ADIOS2 can do I/O
  call adios2_end_step(adios2_file_fwdatt, ier)
  ! Never close the file here. It will be closed in adios2_manager::finalize_adios2()

  end subroutine save_forward_arrays_undoatt_adios2


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are always dumped.
!! \param io The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_common_forward_arrays_adios2(io)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_io), intent(in) :: io
  type(adios2_variable) :: v

  integer :: ier
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array

  ! crust/mantle
  ldim(1) = NDIM * NGLOB_CRUST_MANTLE
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "displ_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "veloc_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "accel_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  ! inner core
  ldim(1) = NDIM * NGLOB_INNER_CORE
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "displ_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "veloc_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "accel_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  ! outer core
  ldim(1) = NGLOB_OUTER_CORE
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "displ_outer_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "veloc_outer_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "accel_outer_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  end subroutine define_common_forward_arrays_adios2


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are always dumped
!! except for undo attenuation.
!! \param io The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_epsilon_forward_arrays_adios2(io)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_io), intent(in) :: io
  type(adios2_variable) :: v

  integer :: ier
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array

  ! strains
  ldim(1) = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE_STR_OR_ATT
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "epsilondev_xx_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_yy_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_xy_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_xz_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_yz_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  ldim(1) = NGLLX * NGLLY * NGLLZ * NSPEC_INNER_CORE_STR_OR_ATT
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "epsilondev_xx_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_yy_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_xy_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_xz_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "epsilondev_yz_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  end subroutine define_epsilon_forward_arrays_adios2


!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are dumped if ROTATION is true.
!! \param io The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_rotation_forward_arrays_adios2(io)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_io), intent(in) :: io
  type(adios2_variable) :: v

  integer :: ier
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array

  ldim(1) = NGLLX * NGLLY * NGLLZ * NSPEC_OUTER_CORE_ROTATION
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "A_array_rotation", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "B_array_rotation", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  end subroutine define_rotation_forward_arrays_adios2

!-------------------------------------------------------------------------------
!> Define ADIOS forward arrays that are dumped if ATTENUATION is true.
!! \param io The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable

  subroutine define_attenuation_forward_arrays_adios2(io)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_io), intent(in) :: io
  type(adios2_variable) :: v

  integer :: ier
  integer(kind=8), dimension(1) :: gdim  ! Global shape of array
  integer(kind=8), dimension(1) :: ldim  ! Local size of array
  integer(kind=8), dimension(1) :: offs  ! Starting offset in global array

  ! attenuation memory variables
  ! crust/mantle
  ldim(1) = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE_ATTENUATION
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;

  call adios2_define_variable(v, io, "R_xx_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_yy_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_xy_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_xz_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_yz_crust_mantle", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  ! inner core
  ldim(1) = N_SLS*NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE_ATTENUATION
  gdim = sizeprocs_adios2 * ldim;
  offs = myrank_adios2 * ldim;
  call adios2_define_variable(v, io, "R_xx_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_yy_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_xy_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_xz_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)
  call adios2_define_variable(v, io, "R_yz_inner_core", adios2_CUSTOM_REAL, 1, gdim, offs, ldim, .true., ier)

  end subroutine define_attenuation_forward_arrays_adios2

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.
!! \param file The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writing

  subroutine write_common_forward_arrays_adios2(file)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_engine), intent(in) :: file
  integer :: ier

  ! crust/mantle
  call adios2_put(file, STRINGIFY_VAR(displ_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(veloc_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(accel_crust_mantle), ier)

  ! inner core
  call adios2_put(file, STRINGIFY_VAR(displ_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(veloc_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(accel_inner_core), ier)

  ! outer core
  call adios2_put(file, STRINGIFY_VAR(displ_outer_core), ier)
  call adios2_put(file, STRINGIFY_VAR(veloc_outer_core), ier)
  call adios2_put(file, STRINGIFY_VAR(accel_outer_core), ier)

  end subroutine write_common_forward_arrays_adios2


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are always dumped.
!! \param file The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writing

  subroutine write_epsilon_forward_arrays_adios2(file)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_engine), intent(in) :: file
  integer :: ier

  ! strains
  ! crust/mantle
  call adios2_put(file, STRINGIFY_VAR(epsilondev_xx_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_yy_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_xy_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_xz_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_yz_crust_mantle), ier)

  ! inner core
  call adios2_put(file, STRINGIFY_VAR(epsilondev_xx_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_yy_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_xy_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_xz_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(epsilondev_yz_inner_core), ier)

  end subroutine write_epsilon_forward_arrays_adios2


!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ROTATION is true.
!! \param file The handle to the adios bp file

  subroutine write_rotation_forward_arrays_adios2(file)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_engine), intent(in) :: file
  integer :: ier
  call adios2_put(file, STRINGIFY_VAR(A_array_rotation), ier)
  call adios2_put(file, STRINGIFY_VAR(B_array_rotation), ier)

  end subroutine write_rotation_forward_arrays_adios2

!-------------------------------------------------------------------------------
!>  Schedule writes of ADIOS forward arrays that are dumped if ATTENUATION
!!  is true.
!! \param file The handle to the adios bp file
!! \param group_size_inc The number of MPI processes involved in the writing

  subroutine write_attenuation_forward_arrays_adios2(file)

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use adios2
  use manager_adios2
  use forward_adios2_write

  implicit none

  type(adios2_engine), intent(in) :: file
  integer :: ier

  ! attenuation memory variables
  ! crust/mantle
  call adios2_put(file, STRINGIFY_VAR(R_xx_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(R_yy_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(R_xy_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(R_xz_crust_mantle), ier)
  call adios2_put(file, STRINGIFY_VAR(R_yz_crust_mantle), ier)

  ! inner core
  call adios2_put(file, STRINGIFY_VAR(R_xx_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(R_yy_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(R_xy_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(R_xz_inner_core), ier)
  call adios2_put(file, STRINGIFY_VAR(R_yz_inner_core), ier)

  end subroutine write_attenuation_forward_arrays_adios2

