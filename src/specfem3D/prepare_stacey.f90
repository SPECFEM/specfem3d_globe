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



  subroutine prepare_stacey()

! sets up arrays for Stacey conditions

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_outercore

  implicit none
  ! local parameters
  integer :: ier
  integer :: nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm
  integer :: nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc,nabs_zmin_oc
  integer(kind=8) :: filesize

  ! checks if anything to do
  if (.not. ABSORBING_CONDITIONS ) return

  ! sets up absorbing boundary buffer arrays
  if (myrank == 0) then
    write(IMAIN,*) "preparing absorbing boundaries"
    call flush_IMAIN()
  endif

  ! crust_mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

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

  ! allocates buffers
  if (nspec2D_xmin_crust_mantle > 0 .and. SAVE_STACEY) then
    nabs_xmin_cm = nspec2D_xmin_crust_mantle
  else
    nabs_xmin_cm = 1
  endif
  allocate(absorb_xmin_crust_mantle(NDIM,NGLLY,NGLLZ,nabs_xmin_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmin')

  if (nspec2D_xmax_crust_mantle > 0 .and. SAVE_STACEY) then
    nabs_xmax_cm = nspec2D_xmax_crust_mantle
  else
    nabs_xmax_cm = 1
  endif
  allocate(absorb_xmax_crust_mantle(NDIM,NGLLY,NGLLZ,nabs_xmax_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmax')

  if (nspec2D_ymin_crust_mantle > 0 .and. SAVE_STACEY) then
    nabs_ymin_cm = nspec2D_ymin_crust_mantle
  else
    nabs_ymin_cm = 1
  endif
  allocate(absorb_ymin_crust_mantle(NDIM,NGLLX,NGLLZ,nabs_ymin_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymin')

  if (nspec2D_ymax_crust_mantle > 0 .and. SAVE_STACEY) then
    nabs_ymax_cm = nspec2D_ymax_crust_mantle
  else
    nabs_ymax_cm = 1
  endif
  allocate(absorb_ymax_crust_mantle(NDIM,NGLLX,NGLLZ,nabs_ymax_cm),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymax')

  ! file I/O for re-construction of wavefields
  if (nspec2D_xmin_crust_mantle > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_xmin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmin_crust_mantle)

    ! total file size
    filesize = reclen_xmin_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    endif
  endif
  if (nspec2D_xmax_crust_mantle > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_xmax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmax_crust_mantle)

    ! total file size
    filesize = reclen_xmax_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    endif
  endif
  if (nspec2D_ymin_crust_mantle > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_ymin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymin_crust_mantle)

    ! total file size
    filesize = reclen_ymin_crust_mantle
    filesize = filesize * NSTEP


    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    endif
  endif
  if (nspec2D_ymax_crust_mantle > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_ymax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymax_crust_mantle)

    ! total file size
    filesize = reclen_ymax_crust_mantle
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif


  ! outer_core
  ! create name of database
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  ! allocates buffers
  ! xmin
  if (nspec2D_xmin_outer_core > 0 .and. SAVE_STACEY) then
    nabs_xmin_oc = nspec2D_xmin_outer_core
  else
    nabs_xmin_oc = 1
  endif
  allocate(absorb_xmin_outer_core(NGLLY,NGLLZ,nabs_xmin_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmin')

  ! xmax
  if (nspec2D_xmax_outer_core > 0 .and. SAVE_STACEY) then
    nabs_xmax_oc = nspec2D_xmax_outer_core
  else
    nabs_xmax_oc = 1
  endif
  allocate(absorb_xmax_outer_core(NGLLY,NGLLZ,nabs_xmax_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb xmax')

  ! ymin
  if (nspec2D_ymin_outer_core > 0 .and. SAVE_STACEY) then
    nabs_ymin_oc = nspec2D_ymin_outer_core
  else
    nabs_ymin_oc = 1
  endif
  allocate(absorb_ymin_outer_core(NGLLX,NGLLZ,nabs_ymin_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymin')

  ! ymax
  if (nspec2D_ymax_outer_core > 0 .and. SAVE_STACEY) then
    nabs_ymax_oc = nspec2D_ymax_outer_core
  else
    nabs_ymax_oc = 1
  endif
  allocate(absorb_ymax_outer_core(NGLLX,NGLLZ,nabs_ymax_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb ymax')

  ! zmin
  if (nspec2D_zmin_outer_core > 0 .and. SAVE_STACEY) then
    nabs_zmin_oc = nspec2D_zmin_outer_core
  else
    nabs_zmin_oc = 1
  endif
  allocate(absorb_zmin_outer_core(NGLLX,NGLLY,nabs_zmin_oc),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating absorb zmin')

  ! file I/O for re-construction of wavefields
  ! xmin
  if (nspec2D_xmin_outer_core > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_xmin_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmin_outer_core)

    ! total file size
    filesize = reclen_xmin_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif
  ! xmax
  if (nspec2D_xmax_outer_core > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_xmax_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmax_outer_core)

    ! total file size
    filesize = reclen_xmax_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
   endif
  endif
  ! ymin
  if (nspec2D_ymin_outer_core > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_ymin_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymin_outer_core)

    ! total file size
    filesize = reclen_ymin_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    endif
  endif
  ! ymanx
  if (nspec2D_ymax_outer_core > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_ymax_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymax_outer_core)

    ! total file size
    filesize = reclen_ymax_outer_core
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif
  ! zmin
  if (nspec2D_zmin_outer_core > 0 .and. SAVE_STACEY) then

    ! size of single record
    reclen_zmin = CUSTOM_REAL * (NGLLX * NGLLY * nspec2D_zmin_outer_core)

    ! total file size
    filesize = reclen_zmin
    filesize = filesize * NSTEP

    if (SIMULATION_TYPE == 3) then
      call open_file_abs_r(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    endif
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_stacey

