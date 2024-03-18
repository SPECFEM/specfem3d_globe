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

!--------------------------------------------------------------------------------------------------
! GLL
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

  subroutine read_gll_model_adios(rank)

  use constants, only: MAX_STRING_LEN,IMAIN,NGLLX,NGLLY,NGLLZ,PATHNAME_GLL_modeldir,myrank, &
    NGLLX,NGLLY,NGLLZ,IREGION_INNER_CORE

  use shared_parameters, only: NSPEC_REGIONS

  use adios_helpers_mod
  use manager_adios

  use model_gll_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer(kind=8) :: local_dim
  character(len=MAX_STRING_LEN) :: file_name,model_name
  integer :: nspec,ier,iexist
  logical :: has_innercore,has_innercore_all
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_rho_ic

  ! ADIOS variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel

  ! only crust and mantle
  file_name = get_adios_filename(trim(PATHNAME_GLL_modeldir) // 'model_gll')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'reading in model from ',trim(PATHNAME_GLL_modeldir)
    write(IMAIN,*) '  ADIOS file: ',trim(file_name)
    if (rank /= myrank) write(IMAIN,*) '  mesh slice for rank: ',rank
    call flush_IMAIN()
  endif

  ! setup the ADIOS library to read the file
  call init_adios_group(myadios_group,"GLLReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  local_dim = NGLLX * NGLLY * NGLLZ * MGLL_V%nspec
  start(1) = local_dim * int(rank,kind=8)
  count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  ! reads in model for each partition
  select case (MGLL_TYPE)
  case (1)
    ! isotropic model
    if (rank == 0) then
      write(IMAIN,*) '  ADIOS reads isotropic model values: (rho,vp,vs)'
      call flush_IMAIN()
    endif
    ! vp mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vp/array", MGLL_V%vp_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vs/array", MGLL_V%vs_new(:,:,:,1:MGLL_V%nspec))

  case (2)
    ! transverse isotropic model
    if (rank == 0) then
      write(IMAIN,*) '  ADIOS reads transversely isotropic model values: (rho,vpv,vph,vsv,vsh,eta)'
      call flush_IMAIN()
    endif

    ! WARNING previously wronly name 'vps' in the adios files
    ! vpv/vph mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vpv/array", MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vph/array", MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec))

    ! vsv/vsh mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vsv/array", MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vsh/array", MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec))

    ! eta mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/eta/array", MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec))

  case (3)
    ! azimuthal model
    if (rank == 0) then
      write(IMAIN,*) '  reads azimuthal anisotropic model values: (rho,vpv,vph,vsv,vsh,eta,Gc_prime,Gs_prime,mu0)'
      call flush_IMAIN()
    endif

    ! vpv/vph mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vpv/array", MGLL_V%vpv_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vph/array", MGLL_V%vph_new(:,:,:,1:MGLL_V%nspec))

    ! vsv/vsh mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vsv/array", MGLL_V%vsv_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/vsh/array", MGLL_V%vsh_new(:,:,:,1:MGLL_V%nspec))

    ! eta mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/eta/array", MGLL_V%eta_new(:,:,:,1:MGLL_V%nspec))

    ! Gc_prime/Gs_prime/mu0 mesh
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/Gc_prime/array", MGLL_V%Gc_prime_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/Gs_prime/array", MGLL_V%Gs_prime_new(:,:,:,1:MGLL_V%nspec))
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg1/mu0/array", MGLL_V%mu0_new(:,:,:,1:MGLL_V%nspec))

  case default
    stop 'Invalid MGLL_TYPE for reading ADIOS GLL model, type not implemented yet'
  end select

  ! rho mesh
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "reg1/rho/array", MGLL_V%rho_new(:,:,:,1:MGLL_V%nspec))

  ! performs reads
  call read_adios_perform(myadios_file)
  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"GLLReader")

  ! inner core
  ! checks with rho if data available
  nspec = NSPEC_REGIONS(IREGION_INNER_CORE)

  allocate(temp_rho_ic(NGLLX,NGLLY,NGLLZ,nspec), stat=ier)
  if (ier /= 0 ) call exit_MPI(rank,'Error allocating inner core temp array')
  temp_rho_ic(:,:,:,:) = 0.0_CUSTOM_REAL

  ! uses density array to check if inner core arrays are available
  model_name = "reg3/rho"

  ! setup the ADIOS library to read the file
  call init_adios_group(myadios_group,"GLLReaderInnerCore")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  ! checks read of associated model array
  call read_adios_array_gll_check(myadios_file,myadios_group,myrank,nspec,model_name,temp_rho_ic,iexist)

  ! checks if read was successful
  if (iexist == 0) then
    has_innercore = .false.
  else
    has_innercore = .true.
  endif

  ! makes sure all processes have same flag
  call any_all_l(has_innercore,has_innercore_all)

  if (has_innercore_all) then
    if (myrank == 0) then
      write(IMAIN,*) '  ADIOS reads inner core model values: (vp,vs,rho)'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! only isotropic model file supported so far for inner core
    MGLL_V_IC%nspec = NSPEC_REGIONS(IREGION_INNER_CORE)

    ! allocates arrays
    allocate( MGLL_V_IC%vp_new(NGLLX,NGLLY,NGLLZ,MGLL_V_IC%nspec), &
              MGLL_V_IC%vs_new(NGLLX,NGLLY,NGLLZ,MGLL_V_IC%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(rank,'Error allocating inner core vp_new,.. arrays')
    MGLL_V_IC%vp_new(:,:,:,:) = 0.0_CUSTOM_REAL
    MGLL_V_IC%vs_new(:,:,:,:) = 0.0_CUSTOM_REAL

    allocate( MGLL_V_IC%rho_new(NGLLX,NGLLY,NGLLZ,MGLL_V_IC%nspec), stat=ier)
    if (ier /= 0 ) call exit_MPI(rank,'Error allocating inner core rho_new,.. arrays')
    MGLL_V_IC%rho_new(:,:,:,:) = temp_rho_ic(:,:,:,:)  ! already read in from GLL check routine above

    ! for vp,vs
    local_dim = NGLLX * NGLLY * NGLLZ * MGLL_V_IC%nspec
    start(1) = local_dim * int(rank,kind=8)
    count(1) = local_dim
    call set_selection_boundingbox(sel, start, count)
    ! vp
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg3/vp/array", MGLL_V_IC%vp_new(:,:,:,1:MGLL_V_IC%nspec))
    ! vs
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   "reg3/vs/array", MGLL_V_IC%vs_new(:,:,:,1:MGLL_V_IC%nspec))

    ! performs reads
    call read_adios_perform(myadios_file)
  endif

  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"GLLReaderInnerCore")

  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*) '  ADIOS reading done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! frees temporary array
  deallocate(temp_rho_ic)

  end subroutine read_gll_model_adios
