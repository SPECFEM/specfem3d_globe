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


subroutine write_coordinates_netcdf(fName)

  ! This writes the corresponding kernel coordinates.

  use netcdf
  use specfem_par
  use specfem_par_crustmantle

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: fName

  integer, parameter :: NDIMS = 2

  integer :: ncid
  integer :: sizeprocs
  integer :: procDimId, kernDimId, varIdRadius, varIdTheta, varIdPhi
  integer :: mpicomm,mpiinfo_null

  integer, dimension(NDIMS) :: start, count, ids

  ! Get the total number of processors.
  call world_size(sizeprocs)
  call world_get_comm(mpicomm)
  call world_get_info_null(mpiinfo_null)

  ! Create parallel NetCDF file.
  call checkNC (nf90_create (fName, IOR(NF90_NETCDF4, NF90_MPIIO), ncid, &
                             comm = mpicomm, info = mpiinfo_null))

  ! Define the processor array.
  call checkNC (nf90_def_dim (ncid, 'glob', NGLOB_CRUST_MANTLE, kernDimId))
  call checkNC (nf90_def_dim (ncid, 'proc', sizeprocs, procDimID))

  ! Sort ids into array.
  ids = (/ kernDimID, procDimID /)

  ! Define the kernel variable.
  call checkNC (nf90_def_var (ncid, 'radius', NF90_float, ids, varIdRadius))
  call checkNC (nf90_def_var (ncid, 'theta',  NF90_float, ids, varIdTheta))
  call checkNC (nf90_def_var (ncid, 'phi',    NF90_float, ids, varIdPhi))

  ! End definitions.
  call checkNC (nf90_enddef (ncid))

  ! Each processor writes one row.
  start = (/ 1, myRank + 1 /)
  count = (/ NGLOB_CRUST_MANTLE, 1 /)
  call checkNC (nf90_put_var (ncid, varIdRadius, xstore_crust_mantle, start = start, count = count))
  call checkNC (nf90_put_var (ncid, varIdTheta,  ystore_crust_mantle, start = start, count = count))
  call checkNC (nf90_put_var (ncid, varIdPhi,    zstore_crust_mantle, start = start, count = count))

  ! Close the netcdf file.
  call checkNC (nf90_close (ncid))

end subroutine write_coordinates_netcdf

!
!-------------------------------------------------------------------------------------------------
!

subroutine write_kernel_netcdf(fName, param)

  ! This guy converts the kernels to global numbering, and writes to a NetCDF file in parallel.

  use netcdf
  use specfem_par
  use specfem_par_crustmantle

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: fName

  real (kind=CUSTOM_REAL), dimension (NGLLX, NGLLY, NGLLZ, NSPEC_CRUST_MANTLE),intent(in) :: param

  ! local parameters
  real (kind=CUSTOM_REAL), dimension (:), allocatable :: paramUnpacked

  integer, parameter :: NDIMS = 2

  integer :: ncid
  integer :: i, j, k, ispec, iglob
  integer :: sizeprocs
  integer :: procDimId, kernDimId, varId
  integer :: mpicomm,mpiinfo_null

  integer, dimension (NDIMS) :: start, count, ids

  ! Get the total number of processors.
  call world_size(sizeprocs)
  call world_get_comm(mpicomm)
  call world_get_info_null(mpiinfo_null)

  ! Allocate the kernel dump array.
  allocate (paramUnpacked (NGLOB_CRUST_MANTLE))

  ! Convert kernels to global numbering.
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          iglob                 = ibool_crust_mantle (i, j, k, ispec)
          paramUnpacked (iglob) = param (i, j, k, ispec)

        enddo
      enddo
    enddo
  enddo

  ! Create the parallel NetCDF file.
  call checkNC (nf90_create (fName, IOR(NF90_NETCDF4, NF90_MPIIO), ncid, &
    comm = mpicomm, info = mpiinfo_null))

  ! Define the processor array.
  call checkNC (nf90_def_dim (ncid, 'glob', NGLOB_CRUST_MANTLE, kernDimId))
  call checkNC (nf90_def_dim (ncid, 'proc', sizeprocs, procDimID))

  ! Sort ids into array.
  ids = (/ kernDimID, procDimID /)

  ! Define the kernel variable.
  call checkNC (nf90_def_var (ncid, 'rawKernel', NF90_float, ids, varId))

  ! Deflate (compress) the variable. Doesn't work in parallel.
  ! call checkNC (nf90_def_var_deflate (ncid, varId, shuOn, comOn, comLvl))

  ! End definition.
  call checkNC (nf90_enddef (ncid))

  ! Each processor writes one row.
  start = (/ 1, myRank + 1 /)
  count = (/ NGLOB_CRUST_MANTLE, 1 /)
  call checkNC (nf90_put_var (ncid, varId, paramUnpacked, start = start, count = count))

  ! Close the netcdf file.
  call checkNC (nf90_close (ncid))

  ! Deallocate the kernel dump array.
  deallocate (paramUnpacked)

end subroutine write_kernel_netcdf

!
!-------------------------------------------------------------------------------------------------
!

subroutine checkNC (status)

  ! This little guy just checks for an error from the NetCDF libraries and throws a tantrum if
  ! one's found.

  use netcdf

  implicit none

  integer, intent (in) :: status

  if (status /= nf90_noerr) then
    print *,'Error: ', trim(nf90_strerror(status))
    stop 'Error Netcdf command failed'
  endif

end subroutine checkNC
