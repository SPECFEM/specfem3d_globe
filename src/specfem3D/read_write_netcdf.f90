subroutine write_kernel_netcdf (fName, param)
  
  use :: mpi
  use :: netcdf
  use :: specfem_par
  use :: specfem_par_crustmantle
  
  implicit none
  
  intent (in) :: param, fName
  
  character (1024) :: fName
  
  real (kind=CUSTOM_REAL), dimension (:), allocatable :: paramUnpacked
  real (kind=CUSTOM_REAL), dimension (NGLLX, NGLLY, NGLLZ, NSPEC_CRUST_MANTLE) :: param
  
  integer, parameter :: shuOn=1
  integer, parameter :: comLvl=9
  integer, parameter :: comOn=1
  integer            :: ncid, stat
  integer            :: i, j, k, ispec, iglob
  
  allocate (paramUnpacked (NGLOB_CRUST_MANTLE))
  
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          
          iglob                 = ibool_crust_mantle (i, j, k, ispec)
          paramUnpacked (iglob) = param (i, j, k, ispec)
          
        end do
      end do
    end do
  end do
  
  call checkNC ( nf90_create (fName, IOR(NF90_NETCDF4, NF90_MPIIO), ncid, &
    comm = MPI_COMM_WORLD, info = MPI_INFO_NULL))
  
end subroutine write_kernel_netcdf

subroutine checkNC (status)
  
  use :: netcdf
  
  implicit none
  
  integer, intent (in) :: status
  
  if (status /= nf90_noerr) then
    print *, trim (nf90_strerror (status))
    stop 2
  end if
  
end subroutine checkNC