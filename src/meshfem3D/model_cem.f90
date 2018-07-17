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

module cem_par

  use constants, only: PI, GRAV, RHOAV, R_EARTH, MAX_STRING_LEN

  implicit none

  double precision :: scaleval
  double precision :: scale_GPa

  integer, dimension (:), allocatable :: regCode
  integer, parameter                  :: shuOn=1
  integer, parameter                  :: comLvl=9
  integer, parameter                  :: comOn=1
  integer, parameter                  :: NDIMS=3
  integer                             :: nnodes_cem
  integer                             :: nelem_cem
  integer                             :: rank

  real, dimension (:,:), allocatable :: xyzOut

  type par

    double precision, dimension (:), allocatable :: vsv, vsh, vpv, vph, rho
    double precision, dimension (:), allocatable :: c11, c12, c13, c14
    double precision, dimension (:), allocatable :: c15, c16, c26, c33
    double precision, dimension (:), allocatable :: c22, c23, c24, c25
    double precision, dimension (:), allocatable :: c34, c35, c36, c44
    double precision, dimension (:), allocatable :: c55, c56, c66, c45
    double precision, dimension (:), allocatable :: c46

  end type par

  type (par) :: reg1Bc, reg2Bc, reg3Bc

end module cem_par

!
!-------------------------------------------------------------------------------------------
!

  subroutine model_cem_broadcast()

  use constants, only: myrank
  use cem_par
  use netcdf
  use meshfem3D_models_par, only: CEM_ACCEPT

  integer              :: wSize

  ! initializes
  rank = myrank
  call world_size (wSize)

  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_GPa = (RHOAV / 1000.d0) * ((R_EARTH * scaleval / 1000.d0) ** 2)

  if (CEM_ACCEPT) then

    call return_populated_arrays (reg1Bc, "vsv", 1)
    call return_populated_arrays (reg2Bc, "vsv", 2)
    call return_populated_arrays (reg3Bc, "vsv", 3)

    call return_populated_arrays (reg1Bc, "rho", 1)
    call return_populated_arrays (reg2Bc, "rho", 2)
    call return_populated_arrays (reg3Bc, "rho", 3)

    call return_populated_arrays (reg1Bc, "vsh", 1)
    call return_populated_arrays (reg2Bc, "vsh", 2)
    call return_populated_arrays (reg3Bc, "vsh", 3)

    call return_populated_arrays (reg1Bc, "vpv", 1)
    call return_populated_arrays (reg2Bc, "vpv", 2)
    call return_populated_arrays (reg3Bc, "vpv", 3)

    call return_populated_arrays (reg1Bc, "vph", 1)
    call return_populated_arrays (reg2Bc, "vph", 2)
    call return_populated_arrays (reg3Bc, "vph", 3)

    call synchronize_all ()

  endif

  end subroutine model_cem_broadcast

!
!-------------------------------------------------------------------------------------------
!

  subroutine request_cem (vsh, vsv, vph, vpv, rho, iregion_code, ispec, i, j, k)

  use cem_par
  use constants

  use meshfem3D_par, only: ibool

  implicit none

  double precision, intent (out)   :: vsh, vsv, vph, vpv, rho

  integer, intent (in) :: iregion_code, ispec, i, j, k
  integer              :: iglob

  if (iregion_code == IREGION_CRUST_MANTLE) then

    iglob = ibool(i,j,k,ispec)

    rho = reg1Bc%rho(iglob) * 1000.0d0 / (RHOAV)
    vpv = reg1Bc%vpv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsv = reg1Bc%vsv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsh = reg1Bc%vsh(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vph = reg1Bc%vph(iglob) * 1000.0d0 / (R_EARTH * scaleval)

  else if (iregion_code == IREGION_OUTER_CORE) then

    iglob = ibool(i,j,k,ispec)

    rho = reg2Bc%rho(iglob) * 1000.0d0 / (RHOAV)
    vpv = reg2Bc%vpv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsv = reg2Bc%vsv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsh = reg2Bc%vsh(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vph = reg2Bc%vph(iglob) * 1000.0d0 / (R_EARTH * scaleval)

  else if (iregion_code == IREGION_INNER_CORE) then

    iglob = ibool(i,j,k,ispec)

    rho = reg3Bc%rho(iglob) * 1000.0d0 / (RHOAV)
    vpv = reg3Bc%vpv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsv = reg3Bc%vsv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsh = reg3Bc%vsh(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vph = reg3Bc%vph(iglob) * 1000.0d0 / (R_EARTH * scaleval)

  endif

  end subroutine request_cem

!
!-------------------------------------------------------------------------------------------
!

  subroutine return_populated_arrays (structure_p, param, reg)

  use cem_par
  use netcdf

  implicit none

  integer, intent (in) :: reg
  integer :: ncid, dimid, nPar
  integer :: status, varid

  character(3), intent(in) :: param

  type (par) :: structure_p

  character (len=MAX_STRING_LEN) :: fileName
  character (len=MAX_STRING_LEN) :: fileNameTrim
  character (len=MAX_STRING_LEN) :: formatString

  formatString = "(A,A3,A,I0.2,A,I0.6,A)"
  write (fileName,formatString) "./DATA/cemRequest/", param, "_reg", reg, ".proc", rank, ".nc"
  fileNameTrim = trim(fileName)

  status = nf90_open              (fileNameTrim, NF90_NOWRITE, ncid)
  status = nf90_inq_dimid         (ncid, "param", dimid)
  status = nf90_inquire_dimension (ncid, dimid, len = nPar)

  if ( param == "vsv" ) allocate(structure_p%vsv(nPar))
  if ( param == "rho" ) allocate(structure_p%rho(nPar))
  if ( param == "vpv" ) allocate(structure_p%vpv(nPar))
  if ( param == "vph" ) allocate(structure_p%vph(nPar))
  if ( param == "vsh" ) allocate(structure_p%vsh(nPar))

  status = nf90_inq_varid (ncid, "data", varid)

  if ( param == "vsv" ) status = nf90_get_var (ncid, varid, structure_p%vsv)
  if ( param == "rho" ) status = nf90_get_var (ncid, varid, structure_p%rho)
  if ( param == "vpv" ) status = nf90_get_var (ncid, varid, structure_p%vpv)
  if ( param == "vph" ) status = nf90_get_var (ncid, varid, structure_p%vph)
  if ( param == "vsh" ) status = nf90_get_var (ncid, varid, structure_p%vsh)

  end subroutine return_populated_arrays

!
!-------------------------------------------------------------------------------------------
!

  subroutine write_cem_request (iregion_code)

  use netcdf
  use constants
  use CEM_par

  implicit none

  integer, parameter :: NDIMS_WRITE=2
  integer, dimension (NDIMS_WRITE) :: start, count, ids

  integer :: ncid, paramDimID, procDimID, varidX, varidY, varidZ
  integer :: varidR, iregion_code, commWorldSize, Comm, info
  character (len = MAX_STRING_LEN) :: fileName, fileNameTrim, formatString

  ! Get the total number of processors.
  call world_size(commWorldSize)
  call world_duplicate(Comm)
  call world_get_info_null(info)

  ! Define filename.
  formatString = "(A,I0.2,A)"
  write (fileName, formatString) "DATA/cemRequest/xyz_reg", iregion_code, ".nc"
  fileNameTrim = trim (fileName)

  ! Create parallel NetCDF file.
  call checkNC (nf90_create (fileNameTrim, IOR(NF90_NETCDF4, NF90_MPIIO), ncid, comm = Comm, info = info))

  ! Define the processor array.
  call checkNC (nf90_def_dim (ncid, 'glob', size (xyzOut(:,1)), paramDimID))
  call checkNC (nf90_def_dim (ncid, 'proc', commWorldSize,      procDimID))

  ! Sort ids into array.
  ids = (/ paramDimID, procDimID /)

  ! Define the kernel variable.
  call checkNC (nf90_def_var (ncid, 'x', NF90_float, ids, varidX))
  call checkNC (nf90_def_var (ncid, 'y', NF90_float, ids, varidY))
  call checkNC (nf90_def_var (ncid, 'z', NF90_float, ids, varidZ))
  call checkNC (nf90_def_var (ncid, 'r', NF90_SHORT, ids, varidR))

  ! End definitions.
  call checkNC (nf90_enddef (ncid))

  ! Each processor writes one row.
  start = (/ 1, myRank + 1 /)
  count = (/ size (xyzOut(:,1)), 1 /)
  call checkNC (nf90_put_var (ncid, varidX, xyzOut(:,1), start = start, count = count))
  call checkNC (nf90_put_var (ncid, varidY, xyzOut(:,2), start = start, count = count))
  call checkNC (nf90_put_var (ncid, varidZ, xyzOut(:,3), start = start, count = count))
  call checkNC (nf90_put_var (ncid, varidR, regCode,     start = start, count = count))

  ! Close the file.
  call checkNC (nf90_close (ncid))

  deallocate(xyzOut)
  deallocate(regCode)

  end subroutine write_cem_request

!
!-------------------------------------------------------------------------------------------
!


  subroutine build_global_coordinates (iregion_code)

  use constants
  use cem_par

  use meshfem3D_par, only: &
    nspec, nglob, &
    ibool,xstore,ystore,zstore

  implicit none

  integer, intent (in) :: iregion_code
  integer              :: i, j, k, iglob, ispec, region

  double precision, parameter :: R_020_KM=6351.0d0, R_052_KM=6319.0d0
  double precision, parameter :: R_100_KM=6271.0d0, R_400_KM=5971.0d0
  double precision, parameter :: R_670_KM=5701.0d0, R_CMB_KM=3480.0d0
  double precision, parameter :: R_ICB_KM=1221.0d0, R_THO_KM=5371.0d0
  double precision            :: x, y, z, rad

  allocate(xyzOut(nglob,NDIMS))
  allocate(regCode(nglob))

  ! Global to local x co-ordinate numbering.
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          x               = sngl(xstore(i,j,k,ispec)) * R_EARTH_KM
          iglob           = ibool(i,j,k,ispec)
          xyzOut(iglob,1) = sngl(x)

        enddo
      enddo
    enddo
  enddo

  ! Global to local y co-ordinate numbering.
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          y               = sngl(ystore(i,j,k,ispec)) * R_EARTH_KM
          iglob           = ibool(i,j,k,ispec)
          xyzOut(iglob,2) = sngl(y)

        enddo
      enddo
    enddo
  enddo

  ! Global to local z co-ordinate numbering.
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          z               = sngl(zstore(i,j,k,ispec)) * R_EARTH_KM
          iglob           = ibool(i,j,k,ispec)
          xyzOut(iglob,3) = sngl(z)

        enddo
      enddo
    enddo
  enddo

  if (iregion_code == 1) then

    do ispec = 1,nspec
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            x               = sngl(xstore(i,j,k,ispec)) * R_EARTH_KM
            iglob           = ibool(i,j,k,ispec)
            xyzOut(iglob,1) = sngl(x)

            y               = sngl(ystore(i,j,k,ispec)) * R_EARTH_KM
            iglob           = ibool(i,j,k,ispec)
            xyzOut(iglob,2) = sngl(y)

            z               = sngl(zstore(i,j,k,ispec)) * R_EARTH_KM
            iglob           = ibool(i,j,k,ispec)
            xyzOut(iglob,3) = sngl(z)

            rad = dsqrt ( x * x + y * y + z * z )

            if (rad < R_670_KM .and. rad >= R_CMB_KM) then
              region = 3
            else if (rad < R_400_KM .and. rad >= R_670_KM) then
              region = 2
            else if (rad >= R_400_KM) then
              region = 1
            endif

            iglob          = ibool(i,j,k,ispec)
            regCode(iglob) = region

          enddo
        enddo
      enddo
    enddo

  else if (iregion_code == 2) then

    regCode(:) = 4

  else if (iregion_code == 3) then

    regCode(:) = 5

  endif

  end subroutine build_global_coordinates

!
!-------------------------------------------------------------------------------------------
!

  subroutine checkNC (status)

  ! This little guy just checks for an error from the NetCDF libraries and throws a tantrum if
  ! one's found.

  use netcdf

  implicit none

  integer, intent (in) :: status

  if (status /= nf90_noerr) then
    print *, trim (nf90_strerror (status))
    stop 'Netcdf error.'
  endif

  end subroutine checkNC
