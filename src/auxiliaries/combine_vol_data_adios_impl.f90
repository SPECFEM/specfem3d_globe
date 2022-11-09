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

module combine_vol_data_adios_mod

  use adios_helpers_mod
  use manager_adios

  implicit none

contains

!=============================================================================
!> Print help message.
subroutine print_usage_adios()

  implicit none
  print *, ' Usage: '
  print *, '   xcombine_data slice_list varname var_file mesh_file output_dir high/low-resolution region'
  print *
  print *, ' with'
  print *, '   slice_list   - text file containing slice numbers to combine (or use name "all" for all slices)'
  print *, '   varname      - possible varnames are: '
  print *, '                    rho, vp, vs, kappastore, mustore, alpha_kl, beta_kl, etc.'
  print *, '   var_file     - datafile that holds array, as real(kind=CUSTOM_REAL):: varname(NGLLX,NGLLY,NGLLZ,NSPEC),'
  print *, '                  (e.g. OUTPUT_FILES/kernels.bp)'
  print *, '   mesh_file    - are used to link variable to the topology (e.g. DATABASES_MPI/solver_data.bp)'
  print *, '   output_dir   - indicates where var_name.vtk will be written'
  print *, '   high/low res - give 0 for low resolution and 1 for high resolution'
  print *, '   region       - (optional) region number, only use 1 == crust/mantle, 2 == outer core, 3 == inner core'
  print *

  stop ' Reenter command line options'

end subroutine print_usage_adios

!=============================================================================
!> Interpret command line arguments

subroutine read_args_adios(arg, var_name, value_file_name, mesh_file_name, slice_list_name, &
                           outdir, ires, iregion)

  use constants, only: IIN,MAX_STRING_LEN

  implicit none
  ! Arguments
  character(len=*), intent(in) :: arg(:)
  integer, intent(out) :: ires, iregion
  character(len=*), intent(out) :: var_name, value_file_name, mesh_file_name, &
                                   outdir, slice_list_name

  ! initializes
  iregion = 0

  ! gets arguments
  if ((command_argument_count() == 6) .or. (command_argument_count() == 7)) then
    slice_list_name = arg(1)
    var_name = arg(2)
    value_file_name = arg(3)
    mesh_file_name = arg(4)
    outdir = arg(5)
    read(arg(6),*) ires
  else
    call print_usage_adios()
  endif
  if (command_argument_count() == 7) then
    read(arg(7),*) iregion
  endif

  !debug
  !print *,'debug: read adios: arguments: ',trim(slice_list_name),"|",trim(var_name),"|", &
  !        trim(value_file_name),"|",trim(mesh_file_name),"|",trim(outdir),"|",ires,"|",iregion

end subroutine read_args_adios


!=============================================================================
!> Open ADIOS value and mesh files, read mode

subroutine init_adios(value_file_name, mesh_file_name)

  implicit none
  ! Parameters
  character(len=*), intent(in) :: value_file_name, mesh_file_name

  ! debug
  logical, parameter :: DEBUG = .false.
  integer :: nglob,nspec

  ! initializes adios
  call initialize_adios()

  ! initializes read method and opens mesh file (using default handle)
  call init_adios_group(myadios_group,"MeshReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,mesh_file_name)

  ! opens second adios file for reading data values
  call init_adios_group(myadios_val_group,"ValReader")
  call open_file_adios_read(myadios_val_file,myadios_val_group,value_file_name)

  ! debug output list variables and attributs in mesh file
  if (DEBUG) then
    call show_adios_file_variables(myadios_file,myadios_group,mesh_file_name)
    call show_adios_file_variables(myadios_val_file,myadios_val_group,value_file_name)
    ! checks reading nglob/nspec values
    print *,'mesh file: ',trim(mesh_file_name)
    call read_adios_scalar(myadios_file,myadios_group,0,"reg1/nglob",nglob)
    call read_adios_scalar(myadios_file,myadios_group,0,"reg1/nspec",nspec)
    print *,'  nglob/nspec = ',nglob,"/",nspec
  endif

end subroutine init_adios


!=============================================================================
!> Open ADIOS value and mesh files, read mode

subroutine clean_adios()

  implicit none

  ! closes file with data values
  call close_file_adios_read(myadios_val_file)

  ! closes default file and finalizes read method
  call close_file_adios_read_and_finalize_method(myadios_file)

  ! finalizes
  call finalize_adios()

end subroutine clean_adios


!=============================================================================

subroutine read_scalars_adios_mesh(iproc, ir, nglob, nspec)

  implicit none
  ! Parameters
  integer, intent(in) :: iproc, ir
  integer, intent(out) :: nglob, nspec
  ! Variables
  !integer(kind=8) :: sel
  character(len=80) :: reg_name
  ! debug
  logical, parameter :: DEBUG = .false.
  !character(len=1024) :: err_message
  !integer :: ier

  ! region name
  write(reg_name,"('reg',i1,'/')") ir

  ! debug output
  if (DEBUG) then
    call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nglob",nglob)
    ! note: adios_get_scalar here retrieves the same nglob for everyone (from writer rank 0)
    !call adios_get_scalar(myadios_file, trim(reg_name)//"nglob", nglob, ier)
    !if (ier /= 0) then
    !  call adios_errmsg(err_message)
    !  print *,'Error adios: could not read parameter: ',trim(reg_name)//"nglob"
    !  print *,trim(err_message)
    !  stop 'Error adios helper read scalar'
    !endif

    call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nspec",nspec)
    !call adios_get_scalar(myadios_file, trim(reg_name)//"nspec", nspec, ier)
    !if (ier /= 0) then
    !  call adios_errmsg(err_message)
    !  print *,'Error adios: could not read parameter: ',trim(reg_name)//"nspec"
    !  print *,trim(err_message)
    !  stop 'Error adios helper read scalar'
    !endif

    print *,'read nglob/nspec: ',iproc,ir,nglob,nspec
  endif

  ! reads nglob & nspec
  ! the following adios calls allow to have different nglob/nspec values for different processes (iproc)
  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nglob",nglob)
  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "nspec",nspec)

end subroutine read_scalars_adios_mesh


!=============================================================================

subroutine read_coordinates_adios_mesh(iproc, ir, nglob, nspec, &
                                       xstore, ystore, zstore, ibool)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none
  ! Parameters
  integer, intent(in) :: iproc, ir, nglob, nspec
  real(kind=CUSTOM_REAL),dimension(:), intent(inout) :: xstore, ystore, zstore
  integer, dimension(:,:,:,:), intent(inout) :: ibool
  ! Variables
  character(len=80) :: reg_name
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel_coord, sel_ibool
  integer(kind=8) :: offset_coord, offset_ibool

  write(reg_name,"('reg',i1,'/')") ir

  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "ibool/offset",offset_ibool)
  call read_adios_scalar(myadios_file,myadios_group,iproc,trim(reg_name) // "x_global/offset",offset_coord)

  start(1) = offset_ibool
  count(1) = NGLLX * NGLLY * NGLLZ * nspec
  call set_selection_boundingbox(sel_ibool , start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel_ibool, start, count, trim(reg_name) // "ibool/array", ibool)

  start(1) = offset_coord
  count(1) = nglob
  call set_selection_boundingbox(sel_coord, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel_coord, start, count, trim(reg_name) // "x_global/array", xstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel_coord, start, count, trim(reg_name) // "y_global/array", ystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel_coord, start, count, trim(reg_name) // "z_global/array", zstore)

  call read_adios_perform(myadios_file)

  call delete_adios_selection(sel_ibool)
  call delete_adios_selection(sel_coord)

end subroutine read_coordinates_adios_mesh


!=============================================================================
!> reads in data from ADIOS value file

subroutine read_values_adios(var_name, iproc, ir, nspec, data)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,IREGION_CRUST_MANTLE,IREGION_INNER_CORE,IREGION_OUTER_CORE

  implicit none
  ! Parameters
  character(len=*), intent(in) :: var_name
  integer, intent(in) :: iproc, ir, nspec
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec), intent(inout) :: data
  ! Variables
  integer(kind=8), dimension(1) :: start, count
  integer(kind=8) :: sel
  integer(kind=8) :: offset
  character(len=128) :: data_name
  character(len=8) :: reg_name
  logical :: is_kernel

  ! note: we can either visualize
  !         wavespeed arrays (in DATABASES_MPI/model_gll.bp)
  !                          (e.g. rho,vp,vs,vph,vpv,vsh,..)
  !       or
  !         sensitivity kernels (in OUTPUT_FILES/kernels.bp)
  !                             (e.g. rho_kl,alpha_kl,beta_kl,alphah_kl,alphav_kl,betah_kl,..)
  !
  !       unfortunately, they have different naming conventions for different Earth regions:
  !        - wavespeed arrays: reg1/vp/, reg2/vp/, reg3/vp/, ..
  !        - kernels: alpha_kl_crust_mantle/, alpha_kl_outer_core/, alpha_kl_inner_core/,..
  ! here we try to estimate the type by the ending of the variable name given by the user,
  ! i.e. if the ending is "_kl" we assume it is a kernel name
  !
  is_kernel = .false.
  ! example: alpha_kl checks ending '_kl'
  if (len_trim(var_name) > 3) then
    if (var_name(len_trim(var_name)-2:len_trim(var_name)) == '_kl') then
      is_kernel = .true.
    endif
  endif
  ! example: alpha_kl_crust_mantle checks ending '_kl_crust_mantle'
  if (len_trim(var_name) > 16) then
    if (var_name(len_trim(var_name)-15:len_trim(var_name)) == '_kl_crust_mantle') then
      is_kernel = .true.
    endif
  endif

  ! determines full data array name
  if (is_kernel) then
    ! for kernel name: alpha_kl, .. , **_kl only:
    ! adds region name to get kernel_name
    ! for example: var_name = "alpha_kl" -> alpha_kl_crust_mantle
    !
    ! note: this must match the naming convention used
    !       in file save_kernels_adios.F90
    if (var_name(len_trim(var_name)-2:len_trim(var_name)) == '_kl') then
      select case (ir)
      case (IREGION_CRUST_MANTLE)
        data_name = trim(var_name) // "_crust_mantle"
      case (IREGION_OUTER_CORE)
        data_name = trim(var_name) // "_outer_core"
      case (IREGION_INNER_CORE)
        data_name = trim(var_name) // "_inner_core"
      case default
        stop 'Error wrong region code in read_values_adios() routine'
      end select
    else
      ! example: alpha_kl_crust_mantle
      data_name = trim(var_name)
    endif
    print *,'  kernel data name: ',trim(data_name)
  else
    ! for wavespeed name: rho,vp,..
    ! adds region name: var_name = "rho" -> reg1/rho
    write(reg_name,"('reg',i1,'/')") ir
    data_name = trim(reg_name) // trim(var_name)
    print *,'  data name: ',trim(data_name)
  endif

  ! gets data values
  if (.true.) then
    ! default
    ! assumes GLL type array size (NGLLX,NGLLY,NGLLZ,nspec)
    call read_adios_array(myadios_val_file,myadios_val_group,iproc,nspec,trim(data_name),data)
  else
    ! reads in data offset
    call read_adios_scalar(myadios_val_file,myadios_val_group,iproc,trim(data_name) // "/offset",offset)

    ! reads in data array
    start(1) = offset
    count(1) = NGLLX * NGLLY * NGLLZ * nspec
    call set_selection_boundingbox(sel , start, count)

    call read_adios_schedule_array(myadios_val_file, myadios_val_group, sel, start, count, trim(data_name) // "/array", data)
    call read_adios_perform(myadios_val_file)

    call delete_adios_selection(sel)
  endif

end subroutine read_values_adios

end module combine_vol_data_adios_mod
