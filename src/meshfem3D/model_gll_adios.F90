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

!--------------------------------------------------------------------------------------------------
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

  subroutine read_gll_model_adios(rank)

  use constants, only: MAX_STRING_LEN,IMAIN,NGLLX,NGLLY,NGLLZ,PATHNAME_GLL_modeldir, &
    myrank

  use adios_read_mod
  use adios_helpers_mod
  use manager_adios, only: open_file_adios_read,file_handle_adios

  use model_gll_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer :: local_dim
  character(len=MAX_STRING_LEN) :: file_name

  ! ADIOS variables
  integer :: adios_err
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel

  ! only crust and mantle
  write(file_name,'(a)') trim(PATHNAME_GLL_modeldir) // 'model_gll.bp'

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)'reading in model from ',trim(PATHNAME_GLL_modeldir)
    write(IMAIN,*)'  ADIOS file: ',trim(file_name)
    if (rank /= myrank) write(IMAIN,*) '  mesh slice for rank: ',rank
    call flush_IMAIN()
  endif

  ! Setup the ADIOS library to read the file
  call open_file_adios_read(file_name)

  local_dim = NGLLX * NGLLY * NGLLZ * MGLL_V%nspec
  start(1) = local_dim * rank
  count(1) = local_dim
  call adios_selection_boundingbox (sel , 1, start, count)

  ! reads in model for each partition
  select case (MGLL_TYPE)
  case (1)
    ! isotropic model
    if (rank == 0) then
      write(IMAIN,*)'  ADIOS reads isotropic model values: vp,vs,rho'
      call flush_IMAIN()
    endif
    ! vp mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/vp/array", 0, 1, &
        MGLL_V%vp_new(:,:,:,1:MGLL_V%nspec), adios_err)
    call adios_schedule_read(file_handle_adios, sel, "reg1/vs/array", 0, 1, &
        MGLL_V%vs_new(:,:,:,1:MGLL_V%nspec), adios_err)
  case (2)
    ! transverse isotropic model
    if (rank == 0) then
      write(IMAIN,*)'  ADIOS reads transversely isotropic model values: vpv,vph,vsv,vsh,eta,rho'
      call flush_IMAIN()
    endif

    ! WARNING previously wronly name 'vps' in the adios files
    ! vpv/vph mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/vpv/array", 0, 1, &
        MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec), adios_err)
    call adios_schedule_read(file_handle_adios, sel, "reg1/vph/array", 0, 1, &
        MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec), adios_err)

    ! vsv/vsh mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/vsv/array", 0, 1, &
        MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec), adios_err)
    call adios_schedule_read(file_handle_adios, sel, "reg1/vsh/array", 0, 1, &
        MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec), adios_err)

    ! eta mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/eta/array", 0, 1, &
        MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec), adios_err)

  case (3)
    ! azimuthal model
    if (rank == 0) then
      write(IMAIN,*)'  reads azimuthal anisotropic model values: vpv,vph,vsv,vsh,eta,Gc_prime,Gs_prime,rho'
      call flush_IMAIN()
    endif

    ! vpv/vph mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/vpv/array", 0, 1, &
        MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec), adios_err)
    call adios_schedule_read(file_handle_adios, sel, "reg1/vph/array", 0, 1, &
        MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec), adios_err)

    ! vsv/vsh mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/vsv/array", 0, 1, &
        MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec), adios_err)
    call adios_schedule_read(file_handle_adios, sel, "reg1/vsh/array", 0, 1, &
        MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec), adios_err)

    ! eta mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/eta/array", 0, 1, &
        MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec), adios_err)

    ! Gc/Gs mesh
    call adios_schedule_read(file_handle_adios, sel, "reg1/Gc_prime/array", 0, 1, &
        MGLL_V%Gc_prime_new(:,:,:,1:MGLL_V%nspec), adios_err)
    call adios_schedule_read(file_handle_adios, sel, "reg1/Gs_prime/array", 0, 1, &
        MGLL_V%Gs_prime_new(:,:,:,1:MGLL_V%nspec), adios_err)

  case default
    stop 'Invalid MGLL_TYPE for reading ADIOS GLL model, type not implemented yet'
  end select

  ! rho mesh
  call adios_schedule_read(file_handle_adios, sel, "reg1/rho/array", 0, 1, &
      MGLL_V%rho_new(:,:,:,1:MGLL_V%nspec), adios_err)

  call adios_perform_reads(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  ! closes adios file
  call adios_read_close(file_handle_adios, adios_err)
  call check_adios_err(myrank,adios_err)

  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*)'  ADIOS reading done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_gll_model_adios
