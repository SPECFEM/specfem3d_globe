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
! GLL QMU
!
! based on modified GLL mesh output from mesher
!
! used for iterative inversion procedures
!--------------------------------------------------------------------------------------------------

  subroutine read_gll_qmu_model_adios(rank)

  use constants, only: MAX_STRING_LEN,IMAIN,NGLLX,NGLLY,NGLLZ,PATHNAME_GLL_modeldir,myrank

  use adios_helpers_mod
  use manager_adios

  use model_gll_qmu_par

  implicit none

  integer,intent(in) :: rank

  ! local parameters
  integer(kind=8) :: local_dim
  character(len=MAX_STRING_LEN) :: file_name

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

  ! Setup the ADIOS library to read the file
  call init_adios_group(myadios_group,"GLLReaderQmu")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  local_dim = NGLLX * NGLLY * NGLLZ * MGLL_QMU_V%nspec
  start(1) = local_dim * int(rank,kind=8)
  count(1) = local_dim
  call set_selection_boundingbox(sel, start, count)

  ! reads in model for each partition
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 "reg1/qmu/array", MGLL_QMU_V%qmu_new(:,:,:,1:MGLL_QMU_V%nspec))
  call read_adios_perform(myadios_file)

  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"GLLReaderQmu")

  call synchronize_all()

  if (myrank == 0) then
    write(IMAIN,*) '  ADIOS reading done'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine read_gll_qmu_model_adios
