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

!--------------------------------------------------------------------------------------------------
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

subroutine read_gll_model_adios(myrank,MGLL_V,NSPEC)

  use constants
  use adios_read_mod
  use adios_helpers_mod
  use meshfem3D_models_par,only: TRANSVERSE_ISOTROPY

  implicit none

  ! GLL model_variables
  type model_gll_variables
    sequence
    ! tomographic iteration model on GLL points
    double precision :: scale_velocity,scale_density
    ! isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vs_new,vp_new,rho_new
    ! transverse isotropic model
    real(kind=CUSTOM_REAL),dimension(:,:,:,:),pointer :: vsv_new,vpv_new, &
      vsh_new,vph_new,eta_new
    logical :: MODEL_GLL
    logical,dimension(3) :: dummy_pad ! padding 3 bytes to align the structure
  end type model_gll_variables
  type (model_gll_variables) MGLL_V

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC
  integer :: myrank

  ! local parameters
  integer :: local_dim
  character(len=MAX_STRING_LEN) :: file_name
  integer :: comm

  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_handle ! adios_group
  integer(kind=8), dimension(1) :: start, count

  integer(kind=8) :: sel

  ! only crust and mantle
  write(file_name,'(a)') trim(PATHNAME_GLL_modeldir) // 'model_gll.bp'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*)'reading in model from ',trim(PATHNAME_GLL_modeldir)
    write(IMAIN,*)'ADIOS file: ',trim(file_name)
    call flush_IMAIN()
  endif

  call world_get_comm(comm)

  ! Setup the ADIOS library to read the file
  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_open_file (adios_handle, trim(file_name), 0, comm, adios_err)
  if (adios_err /= 0) then
    print*,'Error rank ',myrank,' opening adios file: ',trim(file_name)
    call check_adios_err(myrank,adios_err)
  endif

  local_dim = NGLLX * NGLLY * NGLLZ * nspec(IREGION_CRUST_MANTLE)
  start(1) = local_dim*myrank; count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)

  ! reads in model for each partition
  if (.not. TRANSVERSE_ISOTROPY) then
    ! isotropic model
    ! vp mesh
    call adios_schedule_read(adios_handle, sel, "reg1/vp/array", 0, 1, &
        MGLL_V%vp_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)
    call adios_schedule_read(adios_handle, sel, "reg1/vs/array", 0, 1, &
        MGLL_V%vs_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)
  else
    ! transverse isotropic model
    ! WARNING previously wronly name 'vps' in the adios files
    ! vp mesh
    call adios_schedule_read(adios_handle, sel, "reg1/vpv/array", 0, 1, &
        MGLL_V%vpv_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)
    call adios_schedule_read(adios_handle, sel, "reg1/vph/array", 0, 1, &
        MGLL_V%vph_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)

    ! vs mesh
    call adios_schedule_read(adios_handle, sel, "reg1/vsv/array", 0, 1, &
        MGLL_V%vsv_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)
    call adios_schedule_read(adios_handle, sel, "reg1/vsh/array", 0, 1, &
        MGLL_V%vsh_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)

    ! eta mesh
    call adios_schedule_read(adios_handle, sel, "reg1/eta/array", 0, 1, &
        MGLL_V%eta_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)
  endif

  ! rho mesh
  call adios_schedule_read(adios_handle, sel, "reg1/rho/array", 0, 1, &
      MGLL_V%rho_new(:,:,:,1:nspec(IREGION_CRUST_MANTLE)), adios_err)

  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)

end subroutine read_gll_model_adios
