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



  subroutine prepare_stacey()

! sets up arrays for Stacey conditions

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_outercore

  implicit none
  ! local parameters
  integer :: ier
  integer :: nabs_cm
  integer :: nabs_oc
  integer(kind=8) :: filesize

  ! checks if anything to do
  if (.not. ABSORBING_CONDITIONS ) return

  ! sets up absorbing boundary buffer arrays
  if (myrank == 0) then
    write(IMAIN,*) "preparing absorbing boundaries"
    call flush_IMAIN()
  endif

  ! sets flag to check if we need to save the stacey contributions to file
  if (UNDO_ATTENUATION) then
    ! not needed for undo_attenuation scheme
    SAVE_STACEY = .false.
  else
    ! used for simulation type 1 and 3
    if (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) then
      SAVE_STACEY = .true.
    else
      SAVE_STACEY = .false.
    endif
  endif

  ! crust_mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! allocates buffers
  if (num_abs_boundary_faces_crust_mantle > 0 .and. SAVE_STACEY) then
    nabs_cm = num_abs_boundary_faces_crust_mantle
  else
    nabs_cm = 1
  endif
  allocate(absorb_buffer_crust_mantle(NDIM,NGLLSQUARE,nabs_cm),stat=ier) ! assumes NGLLX == NGLLY == NGLLZ
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb buffer')
  absorb_buffer_crust_mantle(:,:,:) = 0.0_CUSTOM_REAL

  ! file I/O for re-construction of wavefields
  if (num_abs_boundary_faces_crust_mantle > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_absorb_buffer_crust_mantle = CUSTOM_REAL * (NDIM * NGLLSQUARE * num_abs_boundary_faces_crust_mantle)

    ! total file size
    filesize = reclen_absorb_buffer_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(0,trim(prname)//'absorb_buffer.bin',len_trim(trim(prname)//'absorb_buffer.bin'),filesize)
    else
      call open_file_abs_w(0,trim(prname)//'absorb_buffer.bin',len_trim(trim(prname)//'absorb_buffer.bin'),filesize)
    endif
  endif

  ! outer_core
  if (NSPEC_OUTER_CORE > 0) then
    ! create name of database
    call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

    ! allocates buffers
    if (num_abs_boundary_faces_outer_core > 0 .and. SAVE_STACEY) then
      nabs_oc = num_abs_boundary_faces_outer_core
    else
      nabs_oc = 1
    endif
    allocate(absorb_buffer_outer_core(NGLLSQUARE,nabs_oc),stat=ier)  ! assumes NGLLX == NGLLY == NGLLZ
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb buffer')
    absorb_buffer_outer_core(:,:) = 0.0_CUSTOM_REAL

    ! file I/O for re-construction of wavefields
    if (num_abs_boundary_faces_outer_core > 0 .and. SAVE_STACEY) then

      ! size of single record
      reclen_absorb_buffer_outer_core = CUSTOM_REAL * (NGLLSQUARE * num_abs_boundary_faces_outer_core)

      ! total file size
      filesize = reclen_absorb_buffer_outer_core
      filesize = filesize * NSTEP

      if (SIMULATION_TYPE == 3) then
        call open_file_abs_r(4,trim(prname)//'absorb_buffer.bin',len_trim(trim(prname)//'absorb_buffer.bin'),filesize)
      else
        call open_file_abs_w(4,trim(prname)//'absorb_buffer.bin',len_trim(trim(prname)//'absorb_buffer.bin'),filesize)
      endif
    endif
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_stacey

