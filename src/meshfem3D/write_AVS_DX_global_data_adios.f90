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

!-------------------------------------------------------------------------------
!> \file write_AVS_DX_global_adios.f90
!! \brief Define a module to hold global AVS/DX data (points and elements) and
!!        provides function to deal with them.
!! \author MPBL
!-------------------------------------------------------------------------------

!===============================================================================
!> AVS_DX_global_mod module. Hold and write to ADIOS file global data (points
!! and elements).
module AVS_DX_global_mod

  implicit none

  ! ADIOS Arrays to write down
  type avs_dx_global_t
    integer(kind=4) :: npoin, nspec
    real(kind=4), dimension(:), allocatable :: x_adios, y_adios, z_adios
    integer(kind=4), dimension(:), allocatable :: idoubling, iglob1, iglob2, &
        iglob3, iglob4, iglob5, iglob6, iglob7, iglob8
  endtype

contains

!===============================================================================
!> Allocate the structure that hold data to be written; initialize adios vars.
!! \param adios_group ADIOS group where the variables belong
!! \param group_size_inc The size of the ADIOS group to increment
!! \param avs_dx_adios The structure holding the data to be allocated
subroutine define_AVS_DX_global_data_adios(adios_group, myrank, nspec, ibool, &
                                           npointot, mask_ibool, group_size_inc, avs_dx_adios)

  use constants
  use adios_write_mod
  use adios_helpers_mod

  implicit none

  !--- Arguments -------------------------------------------
  integer(kind=8), intent(in) :: adios_group
  integer(kind=4), intent(in) :: nspec, npointot, myrank
  integer(kind=4), intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nspec)
  logical, intent(inout) :: mask_ibool(npointot)
  integer(kind=8), intent(inout) :: group_size_inc
  type(avs_dx_global_t), intent(inout) :: avs_dx_adios
  !--- Variables -------------------------------------------
  integer ispec, npoin, ierr
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8

  ! Dummy arrays for type inference inside adios helpers
  real(kind=4), dimension(1) :: dummy_real1d
  integer(kind=4), dimension(1) :: dummy_int1d

  mask_ibool(:) = .false.

  ! mark global AVS or DX points
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

  ! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

  avs_dx_adios%npoin = npoin
  avs_dx_adios%nspec = nspec
  ! Allocate temporary arrays for AVS/DX points
  allocate(avs_dx_adios%x_adios(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating x_adios.")
  allocate(avs_dx_adios%y_adios(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating y_adios.")
  allocate(avs_dx_adios%z_adios(npoin), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating z_adios.")

  ! Allocate temporary arrays for AVS/DX elements.
  allocate(avs_dx_adios%idoubling(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating idoubling.")
  allocate(avs_dx_adios%iglob1(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob1.")
  allocate(avs_dx_adios%iglob2(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob2.")
  allocate(avs_dx_adios%iglob3(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob3.")
  allocate(avs_dx_adios%iglob4(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob4.")
  allocate(avs_dx_adios%iglob5(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob5.")
  allocate(avs_dx_adios%iglob6(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob6.")
  allocate(avs_dx_adios%iglob7(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob7.")
  allocate(avs_dx_adios%iglob8(nspec), stat=ierr)
  if (ierr /= 0) call exit_MPI(myrank, "Error allocating iglob8.")

  !--- Variables for '...AVS_DXpoints.txt'
  call define_adios_global_array1D(adios_group, group_size_inc, npoin, &
                                   "", "points/x_value", dummy_real1d)
  call define_adios_global_array1D(adios_group, group_size_inc, npoin, &
                                   "", "points/y_value", dummy_real1d)
  call define_adios_global_array1D(adios_group, group_size_inc, npoin, &
                                   "", "points/z_value", dummy_real1d)

  !--- Variables for AVS_DXelements.txt
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/idoubling", dummy_int1d)

  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob1", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob2", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob3", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob4", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob5", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob6", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob7", dummy_int1d)
  call define_adios_global_array1D(adios_group, group_size_inc, nspec, &
                                   "", "elements/num_ibool_AVS_DX_iglob8", dummy_int1d)

end subroutine define_AVS_DX_global_data_adios


!===============================================================================
!> Prepare the global AVS/DX data to be written; fill the structure.
!! \param adios_handle The handle to the ADIOS file to be written.
!! \param myrank The MPI rank of the current process.
!! \param avs_dx_adios The structure to be filled.
!!
!! Create AVS or DX 3D data for the slice, to be recombined in postprocessing.
  subroutine prepare_AVS_DX_global_data_adios(myrank, &
                                              nspec, ibool, idoubling, xstore, ystore, zstore, num_ibool_AVS_DX, &
                                              mask_ibool, npointot, avs_dx_adios)

  use constants
  use adios_write_mod

  implicit none

  integer nspec,myrank
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  ! logical mask used to output global points only once
  integer npointot
  logical mask_ibool(npointot)

  ! numbering of global AVS or DX points
  integer num_ibool_AVS_DX(npointot)

  integer ispec
  integer iglob1,iglob2,iglob3,iglob4,iglob5,iglob6,iglob7,iglob8
  integer npoin,numpoin

  type(avs_dx_global_t), intent(inout) :: avs_dx_adios

! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

! mark global AVS or DX points
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

  ! count global number of AVS or DX points
  npoin = count(mask_ibool(:))

  ! erase the logical mask used to mark points already found
  mask_ibool(:) = .false.

  ! fill the structure with global AVS or DX points
  numpoin = 0
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)
    if (.not. mask_ibool(iglob1)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob1) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,1,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,1,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,1,ispec))
    endif
    if (.not. mask_ibool(iglob2)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob2) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,1,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,1,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,1,ispec))
    endif
    if (.not. mask_ibool(iglob3)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob3) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,1,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,1,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,1,ispec))
    endif
    if (.not. mask_ibool(iglob4)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob4) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,1,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,1,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,1,ispec))
    endif
    if (.not. mask_ibool(iglob5)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob5) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,1,NGLLZ,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,1,NGLLZ,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,1,NGLLZ,ispec))
    endif
    if (.not. mask_ibool(iglob6)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob6) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,1,NGLLZ,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,1,NGLLZ,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,1,NGLLZ,ispec))
    endif
    if (.not. mask_ibool(iglob7)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob7) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(NGLLX,NGLLY,NGLLZ,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(NGLLX,NGLLY,NGLLZ,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(NGLLX,NGLLY,NGLLZ,ispec))
    endif
    if (.not. mask_ibool(iglob8)) then
      numpoin = numpoin + 1
      num_ibool_AVS_DX(iglob8) = numpoin
      avs_dx_adios%x_adios(numpoin) = sngl(xstore(1,NGLLY,NGLLZ,ispec))
      avs_dx_adios%y_adios(numpoin) = sngl(ystore(1,NGLLY,NGLLZ,ispec))
      avs_dx_adios%z_adios(numpoin) = sngl(zstore(1,NGLLY,NGLLZ,ispec))
    endif
    mask_ibool(iglob1) = .true.
    mask_ibool(iglob2) = .true.
    mask_ibool(iglob3) = .true.
    mask_ibool(iglob4) = .true.
    mask_ibool(iglob5) = .true.
    mask_ibool(iglob6) = .true.
    mask_ibool(iglob7) = .true.
    mask_ibool(iglob8) = .true.
  enddo

  ! check that number of global points output is okay
  if (numpoin /= npoin) &
    call exit_MPI(myrank, &
        'incorrect number of global points in AVS or DX file creation')

  ! AVS or DX elements
  do ispec = 1,nspec
    iglob1=ibool(1,1,1,ispec)
    iglob2=ibool(NGLLX,1,1,ispec)
    iglob3=ibool(NGLLX,NGLLY,1,ispec)
    iglob4=ibool(1,NGLLY,1,ispec)
    iglob5=ibool(1,1,NGLLZ,ispec)
    iglob6=ibool(NGLLX,1,NGLLZ,ispec)
    iglob7=ibool(NGLLX,NGLLY,NGLLZ,ispec)
    iglob8=ibool(1,NGLLY,NGLLZ,ispec)

    avs_dx_adios%iglob1 = num_ibool_AVS_DX(iglob1)
    avs_dx_adios%iglob2 = num_ibool_AVS_DX(iglob2)
    avs_dx_adios%iglob3 = num_ibool_AVS_DX(iglob3)
    avs_dx_adios%iglob4 = num_ibool_AVS_DX(iglob4)
    avs_dx_adios%iglob5 = num_ibool_AVS_DX(iglob5)
    avs_dx_adios%iglob6 = num_ibool_AVS_DX(iglob6)
    avs_dx_adios%iglob7 = num_ibool_AVS_DX(iglob7)
    avs_dx_adios%iglob8 = num_ibool_AVS_DX(iglob8)
  enddo
  avs_dx_adios%idoubling = idoubling

  end subroutine prepare_AVS_DX_global_data_adios

!===============================================================================
!> Schedule write to ADIOS file for global AVS/DX data
!! \param adios_handle The handle to the ADIOS file we want to write into
!! \param nspec Number of spectral elements
!! \avs_dx_adios Structure with the data that have to be written
subroutine write_AVS_DX_global_data_adios(adios_handle, myrank, sizeprocs, avs_dx_adios)

  use adios_write_mod
  use adios_helpers_mod

  implicit none
  !--- Arguments
  integer(kind=8), intent(in) :: adios_handle
  integer, intent(in) :: myrank, sizeprocs
  type(avs_dx_global_t), intent(inout) :: avs_dx_adios ! out for adios_write
  !--- Variables
  integer :: npoin, nspec

  npoin = avs_dx_adios%npoin
  nspec = avs_dx_adios%nspec

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, npoin, &
                                   "points/x_value", avs_dx_adios%x_adios)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, npoin, &
                                   "points/y_value", avs_dx_adios%y_adios)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, npoin, &
                                   "points/z_value", avs_dx_adios%z_adios)


  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/idoubling", avs_dx_adios%idoubling)

  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob1", avs_dx_adios%iglob1)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob2", avs_dx_adios%iglob2)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob3", avs_dx_adios%iglob3)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob4", avs_dx_adios%iglob4)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob5", avs_dx_adios%iglob5)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob6", avs_dx_adios%iglob6)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob7", avs_dx_adios%iglob7)
  call write_adios_global_1d_array(adios_handle, myrank, sizeprocs, nspec, &
                                   "elements/num_ibool_AVS_DX_iglob8", avs_dx_adios%iglob8)

end subroutine write_AVS_DX_global_data_adios

!===============================================================================
!> Free temporary structure filled to write AVS/DX global variable to file.
!! \param myrank The MPI rank of the process
!! \param avs_dx_adios The structure holding AVS/DX information
subroutine free_AVS_DX_global_data_adios(avs_dx_adios)

  implicit none
  !--- Arguments
  type(avs_dx_global_t), intent(inout) :: avs_dx_adios

  deallocate(avs_dx_adios%x_adios)
  deallocate(avs_dx_adios%y_adios)
  deallocate(avs_dx_adios%z_adios)

  deallocate(avs_dx_adios%idoubling)
  deallocate(avs_dx_adios%iglob1)
  deallocate(avs_dx_adios%iglob2)
  deallocate(avs_dx_adios%iglob3)
  deallocate(avs_dx_adios%iglob4)
  deallocate(avs_dx_adios%iglob5)
  deallocate(avs_dx_adios%iglob6)
  deallocate(avs_dx_adios%iglob7)
  deallocate(avs_dx_adios%iglob8)

  avs_dx_adios%npoin = 0
  avs_dx_adios%nspec = 0

end subroutine free_AVS_DX_global_data_adios

end module AVS_DX_global_mod
