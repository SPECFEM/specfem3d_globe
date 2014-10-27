!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
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

    double precision, dimension (:), pointer :: vsv, vsh, vpv, vph, rho
    double precision, dimension (:), pointer :: c11, c12, c13, c14
    double precision, dimension (:), pointer :: c15, c16, c26, c33
    double precision, dimension (:), pointer :: c22, c23, c24, c25
    double precision, dimension (:), pointer :: c34, c35, c36, c44
    double precision, dimension (:), pointer :: c55, c56, c66, c45
    double precision, dimension (:), pointer :: c46

  end type par

  type (par) :: reg1Bc, reg2Bc, reg3Bc

end module cem_par

!
!--------------------------------------------------------------------------------------------------
!

subroutine model_cem_broadcast ( myrank )

  use cem_par
  use netcdf
  use meshfem3D_models_par, only: CEM_ACCEPT

  integer, intent (in) :: myrank
  integer              :: wSize

  ! initializes
  rank = myrank
  call world_size (wSize)

  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_GPa = (RHOAV / 1000.d0) * ((R_EARTH * scaleval / 1000.d0) ** 2)

  if ( CEM_ACCEPT ) then

    call return_populated_arrays (reg1Bc, "vsv", 1)
    call return_populated_arrays (reg2Bc, "vsv", 2)
    call return_populated_arrays (reg3Bc, "vsv", 3)

    call return_populated_arrays (reg1Bc, "rho", 1)
    call return_populated_arrays (reg2Bc, "rho", 2)
    call return_populated_arrays (reg3Bc, "rho", 3)

    call return_populated_arrays (reg1Bc, "vsh", 1)
    call return_populated_arrays (reg2Bc, "vsh", 2)
    call return_populated_arrays (reg3Bc, "vsh", 3)

    call return_populated_arrays (reg1Bc, "vpp", 1)
    call return_populated_arrays (reg2Bc, "vpp", 2)
    call return_populated_arrays (reg3Bc, "vpp", 3)

    call synchronize_all ()
  endif

end subroutine model_cem_broadcast

!
!--------------------------------------------------------------------------------------------------
!

subroutine request_cem ( vsh, vsv, vph, vpv, rho, iregion_code, ispec, i, j, k )

  use cem_par
  use constants

  use meshfem3D_par, only: ibool

  implicit none

  double precision, intent (out)   :: vsh, vsv, vph, vpv, rho

  integer, intent (in) :: iregion_code, ispec, i, j, k
  integer              :: iglob

  if      ( iregion_code == IREGION_CRUST_MANTLE ) then

    iglob = ibool(i,j,k,ispec)

    rho = reg1Bc%rho(iglob) * 1000.0d0 / (RHOAV)
    vpv = reg1Bc%vpv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsv = reg1Bc%vsv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsh = reg1Bc%vsh(iglob) * 1000.0d0 / (R_EARTH * scaleval)

  else if ( iregion_code == IREGION_OUTER_CORE ) then

    iglob = ibool(i,j,k,ispec)

    rho = reg2Bc%rho(iglob) * 1000.0d0 / (RHOAV)
    vpv = reg2Bc%vpv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsv = reg2Bc%vsv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsh = reg2Bc%vsh(iglob) * 1000.0d0 / (R_EARTH * scaleval)

  else if ( iregion_code == IREGION_INNER_CORE ) then

    iglob = ibool(i,j,k,ispec)

    rho = reg3Bc%rho(iglob) * 1000.0d0 / (RHOAV)
    vpv = reg3Bc%vpv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsv = reg3Bc%vsv(iglob) * 1000.0d0 / (R_EARTH * scaleval)
    vsh = reg3Bc%vsh(iglob) * 1000.0d0 / (R_EARTH * scaleval)

  endif

  vph = vpv

end subroutine request_cem

!
!--------------------------------------------------------------------------------------------------
!

subroutine return_populated_arrays ( structure, param, reg )

  use cem_par
  use netcdf

  implicit none

  integer, intent (in) :: reg
  integer :: ncid, dimid, nPar
  integer :: status, varid


  character(3), intent(in) :: param

  type (par) :: structure

  character (len=MAX_STRING_LEN) :: fileName
  character (len=MAX_STRING_LEN) :: fileNameTrim
  character (len=MAX_STRING_LEN) :: formatString


  formatString = "(A,I0.2,A,I0.4,A1,A3,A3)"
  write (fileName,formatString) "./DATA/cemRequest/xyz_reg", reg, "_proc", rank, ".", param, ".nc"

  fileNameTrim = trim(fileName)

  status = nf90_open              (fileNameTrim, NF90_NOWRITE, ncid)
  status = nf90_inq_dimid         (ncid, "param", dimid)
  status = nf90_inquire_dimension (ncid, dimid, len = nPar)

  if ( param == "vsv" ) allocate(structure%vsv(nPar))
  if ( param == "rho" ) allocate(structure%rho(nPar))
  if ( param == "vpp" ) allocate(structure%vpv(nPar))
  if ( param == "vsh" ) allocate(structure%vsh(nPar))

  status = nf90_inq_varid (ncid, "data", varid)

  if ( param == "vsv" ) status = nf90_get_var (ncid, varid, structure%vsv)
  if ( param == "rho" ) status = nf90_get_var (ncid, varid, structure%rho)
  if ( param == "vpp" ) status = nf90_get_var (ncid, varid, structure%vpv)
  if ( param == "vsh" ) status = nf90_get_var (ncid, varid, structure%vsh)


end subroutine return_populated_arrays

!
!--------------------------------------------------------------------------------------------------
!

subroutine write_cem_request ( iregion_code )

  use netcdf
  use constants
  use CEM_par

  implicit none

  integer :: ncid, status, x_dimind, y_dimind, z_dimind, varidX, varidY, varidZ
  integer :: varidR, r_dimind, iregion_code
  character (len = MAX_STRING_LEN) :: fileName, fileNameTrim, formatString

!  This line does not seem portable.
!  call execute_command_line ('mkdir -p cemRequest/')

  formatString = "(A,I0.2,A,I0.4)"
  write (fileName,formatString) "DATA/cemRequest/xyz_reg", iregion_code, "_proc", rank
  fileNameTrim = trim(fileName)

  status = nf90_create  (path = fileNameTrim, cmode = NF90_CLOBBER, ncid = ncid)

  status = nf90_def_dim (ncid, "x", size (xyzOut(:,1)), x_dimind )
  status = nf90_def_dim (ncid, "y", size (xyzOut(:,2)), y_dimind )
  status = nf90_def_dim (ncid, "z", size (xyzOut(:,3)), z_dimind )
  status = nf90_def_dim (ncid, "r", size (regCode),     r_dimind )

  status = nf90_def_var (ncid, "dataX", NF90_FLOAT, x_dimind, varidX)
  status = nf90_def_var (ncid, "dataY", NF90_FLOAT, y_dimind, varidY)
  status = nf90_def_var (ncid, "dataZ", NF90_FLOAT, z_dimind, varidZ)
  status = nf90_def_var (ncid, "regC_", NF90_SHORT, r_dimind, varidR)

  status = nf90_def_var_deflate (ncid, varidX, shuOn, comOn, comLvl)
  status = nf90_def_var_deflate (ncid, varidY, shuOn, comOn, comLvl)
  status = nf90_def_var_deflate (ncid, varidZ, shuOn, comOn, comLvl)
  status = nf90_def_var_deflate (ncid, varidR, shuOn, comOn, comLvl)

  status = nf90_enddef  (ncid)

  status = nf90_put_var (ncid, varidX, xyzOut(:,1))
  status = nf90_put_var (ncid, varidY, xyzOut(:,2))
  status = nf90_put_var (ncid, varidZ, xyzOut(:,3))
  status = nf90_put_var (ncid, varidR, regCode)
  status = nf90_close   (ncid)

  deallocate(xyzOut)
  deallocate(regCode)

end subroutine write_cem_request

!
!--------------------------------------------------------------------------------------------------
!

subroutine build_global_coordinates ( nspec, nglob, iregion_code )

  use constants
  use cem_par

  use meshfem3D_par, only: ibool,xstore,ystore,zstore

  implicit none

  integer, intent (in) :: nspec, nglob, iregion_code
  integer              :: i, j, k, iglob, ispec, region

  double precision, parameter :: R_020_KM=6351.0d0, R_052_KM=6319.0d0
  double precision, parameter :: R_100_KM=6271.0d0, R_400_KM=5971.0d0
  double precision, parameter :: R_670_KM=5701.0d0, R_OCR_KM=3480.0d0
  double precision, parameter :: R_ICR_KM=1221.0d0, R_THO_KM=5371.0d0
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
          xyzOut(iglob,1) = x

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
          xyzOut(iglob,2) = y

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
          xyzOut(iglob,3) = z

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
            xyzOut(iglob,1) = x

            y               = sngl(ystore(i,j,k,ispec)) * R_EARTH_KM
            iglob           = ibool(i,j,k,ispec)
            xyzOut(iglob,2) = y

            z               = sngl(zstore(i,j,k,ispec)) * R_EARTH_KM
            iglob           = ibool(i,j,k,ispec)
            xyzOut(iglob,3) = z

            rad = dsqrt( x * x + y * y + z * z )

            if      ( rad < R_THO_KM .and. rad >= R_ICR_KM ) then
              region = 7
            else if ( rad < R_670_KM .and. rad >= R_THO_KM ) then
              region = 6
            else if ( rad < R_400_KM .and. rad >= R_670_KM ) then
              region = 5
            else if ( rad < R_100_KM .and. rad >= R_400_KM ) then
              region = 4
            else if ( rad < R_052_KM .and. rad >= R_100_KM ) then
              region = 3
            else if ( rad < R_020_KM .and. rad >= R_052_KM ) then
              region = 2
            else if ( rad >= R_020_KM ) then
              region = 1
            endif

            iglob          = ibool(i,j,k,ispec)
            regCode(iglob) = region

          enddo
        enddo
      enddo
    enddo

  else if (iregion_code == 2) then

    regCode(:) = 8

  else if (iregion_code == 3) then

    regCode(:) = 9

  endif

end subroutine build_global_coordinates
