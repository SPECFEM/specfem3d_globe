!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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

! include the LIB_VTK_IO library by Stefano Zaghi
! taken from https://github.com/szaghi/Lib_VTK_IO

! this library handles this automatically in particular:
! from http://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python#Datasets_format :
! If you want to write your data in binary format in your "legacy" VTK file (XML files are not concerned),
! you have to write it in BigEndian, and not in LittleEndian, even if you use an Intel or AMD box
! which are known to be in LittleEndian. So, you have to convert your data in BigEndian
! before writing them in your VTK file. This is explicitly mentionned in the VTK file formats documentation

! the "legacy" VTK format has only one file extension for all kind of datasets format: ".vtk",
! which is not the case for the "XML" as you will see later. Writing "legacy" VTK files
! is well documented in the VTK file formats documentation and is quite straightforward.
! You have just to set a few parameters, such as:
!   * file format (ASCII or binary);
!   * data format (char, int, float, double, etc...)
!   * specific information related to the kind of the dataset you want to use.

include "LIB_VTK_IO.f90"

  program read_interpolated_kernel

USE LIB_VTK_IO

  implicit none

!! DK DK with the Intel ifort compiler, compile with the option below for this code to work fine:

! ifort -o xread -assume byterecl read_interpolated_kernel_for_one_chunk_with_binary_kernel_stored_only_and_VTK_output.f90

!  (you will also maybe need   -mcmodel=medium -shared-intel  )

! to debug or test, one can also use:
! -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments
! -warn ignore_loc -warn usage -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv

  integer, parameter :: NX = 512, NY = 512, NZ = 400

  integer, parameter :: subsampling_of_VTK_display = 1 !!!! 4

! number of the chunk to read and display in ParaView / VTK
  integer, parameter :: ichunk = 2

! use single precision reals
  integer, parameter :: SIZE_OF_REAL = 4

! define block type based upon chunk number (between 1 and 6)
! do not change this numbering, chunk AB must be number 1 for central cube
  integer, parameter :: CHUNK_AB = 1
  integer, parameter :: CHUNK_AC = 2
  integer, parameter :: CHUNK_BC = 3
  integer, parameter :: CHUNK_AC_ANTIPODE = 4
  integer, parameter :: CHUNK_BC_ANTIPODE = 5
  integer, parameter :: CHUNK_AB_ANTIPODE = 6

! R_EARTH is the radius of the bottom of the oceans (radius of Earth in m)
  double precision, parameter :: R_EARTH = 6371000.d0

! some useful constants
  double precision, parameter :: PI = 3.141592653589793d0

! RCMB: radius of CMB (m)
  ! default: PREM
  double precision, parameter :: RCMB = 3480000.d0

! local variables
  integer :: ix,iy,iz

  real(kind=4), dimension(NX,NY,NZ) :: alpha_kl

  double precision :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,r_top,r_bottom

  double precision :: gamma,x,y,rgb,rgt,rn
  double precision :: x_bot,y_bot,z_bot
  double precision :: x_top,y_top,z_top

  double precision :: ratio_xi, ratio_eta

  double precision :: xvalue,yvalue,zvalue

  character(len=200) :: filename

!!!!!!---- DK DK for ParaView / VTK output
   integer, parameter :: NX_VTK = 512/subsampling_of_VTK_display, &
                         NY_VTK = 512/subsampling_of_VTK_display, &
                         NZ_VTK = 400/subsampling_of_VTK_display

   integer(I4P),   parameter:: NN = Nx_VTK*Ny_VTK*Nz_VTK       ! total number of nodes

   real(R4P),   dimension(NN) :: XVAL,YVAL,ZVAL  ! x, y and z coordinates
   real(R4P),   dimension(Nx_VTK,Ny_VTK,Nz_VTK) :: X2,Y2,Z2  ! x, y and z coordinates
   equivalence(XVAL,X2)
   equivalence(YVAL,Y2)
   equivalence(ZVAL,Z2)

   real(R4P),   dimension(NN):: var_str_point
   real(R4P),   dimension(Nx_VTK,Ny_VTK,Nz_VTK):: var_str_point2
   equivalence(var_str_point,var_str_point2)

   integer(I4P):: E_IO
   integer :: ix_subsamp,iy_subsamp,iz_subsamp
!!!!!!---- DK DK for ParaView / VTK output

  if(mod(NX,subsampling_of_VTK_display) /= 0) stop 'error: NX must be a multiple of subsampling_of_VTK_display'
  if(mod(NY,subsampling_of_VTK_display) /= 0) stop 'error: NY must be a multiple of subsampling_of_VTK_display'
  if(mod(NZ,subsampling_of_VTK_display) /= 0) stop 'error: NZ must be a multiple of subsampling_of_VTK_display'

  print *
  print *,'processing chunk ',ichunk
  print *

! DK DK: below taken from http://www.atmos.washington.edu/~salathe/osx_unix/endian.html
!     open direct access unformatted file with records sized for the whole array
!     and write the whole array as one record
      write(filename,'(a,i1,a)') 'chunk',ichunk,'_interpolated_data_for_alpha_kernel.bin'
      open (unit=28,file=filename,form='unformatted',access='direct',status='unknown',action='read',recl=NX*NY*NZ*SIZE_OF_REAL)
      read (28,rec=1) alpha_kl
      close(28)

  print *,'minval maxval of alpha_kl = ',minval(alpha_kl),maxval(alpha_kl)

!-------------------------------------------------------------------------------

! the chunks have a size of 90 degrees = PI/2 radians
  ANGULAR_WIDTH_XI_RAD  = PI/2.d0
  ANGULAR_WIDTH_ETA_RAD = PI/2.d0

      r_top = R_EARTH
      r_bottom = RCMB

 iz_subsamp = 0
 do iz = 1,NZ,subsampling_of_VTK_display
      iz_subsamp = iz_subsamp + 1

      rn = dble(iz-1) / dble(NZ-1)

  iy_subsamp = 0
  do iy = 1,NY,subsampling_of_VTK_display
      iy_subsamp = iy_subsamp + 1

      ratio_eta = dble(iy-1) / dble(NY-1)
      y = 2.d0*ratio_eta-1
      y = tan((ANGULAR_WIDTH_ETA_RAD/2.d0) * y)

  ix_subsamp = 0
  do ix = 1,NX,subsampling_of_VTK_display
      ix_subsamp = ix_subsamp + 1

      ratio_xi = dble(ix-1) / dble(NX-1)
      x = 2.d0*ratio_xi-1
      x = tan((ANGULAR_WIDTH_XI_RAD/2.d0) * x)

      gamma = 1.d0 / sqrt(1.d0 + x*x + y*y)

      rgt = (r_top / R_EARTH)*gamma
      rgb = (r_bottom / R_EARTH)*gamma

    ! define the mesh points on the top and the bottom in the six regions of the cubed shpere
      select case(ichunk)

        case(CHUNK_AB)

          x_top = -y*rgt
          y_top = x*rgt
          z_top = rgt

          x_bot = -y*rgb
          y_bot = x*rgb
          z_bot = rgb

        case(CHUNK_AB_ANTIPODE)

          x_top = -y*rgt
          y_top = -x*rgt
          z_top = -rgt

          x_bot = -y*rgb
          y_bot = -x*rgb
          z_bot = -rgb

        case(CHUNK_AC)

          x_top = -y*rgt
          y_top = -rgt
          z_top = x*rgt

          x_bot = -y*rgb
          y_bot = -rgb
          z_bot = x*rgb

        case(CHUNK_AC_ANTIPODE)

          x_top = -y*rgt
          y_top = rgt
          z_top = -x*rgt

          x_bot = -y*rgb
          y_bot = rgb
          z_bot = -x*rgb

        case(CHUNK_BC)

          x_top = -rgt
          y_top = y*rgt
          z_top = x*rgt

          x_bot = -rgb
          y_bot = y*rgb
          z_bot = x*rgb

        case(CHUNK_BC_ANTIPODE)

          x_top = rgt
          y_top = -y*rgt
          z_top = x*rgt

          x_bot = rgb
          y_bot = -y*rgb
          z_bot = x*rgb

        case default
          stop 'incorrect chunk number in compute_coord_main_mesh'

      end select

    ! compute the position of the point in the cubed sphere
      xvalue = x_top*rn + x_bot*(1.d0-rn)
      yvalue = y_top*rn + y_bot*(1.d0-rn)
      zvalue = z_top*rn + z_bot*(1.d0-rn)

  x2(ix_subsamp,iy_subsamp,iz_subsamp) = sngl(xvalue)
  y2(ix_subsamp,iy_subsamp,iz_subsamp) = sngl(yvalue)
  z2(ix_subsamp,iy_subsamp,iz_subsamp) = sngl(zvalue)

  var_str_point2(ix_subsamp,iy_subsamp,iz_subsamp) = alpha_kl(ix,iy,iz)

  enddo
  enddo
  enddo

  print *,'minval maxval of x = ',minval(x2),maxval(x2)
  print *,'minval maxval of y = ',minval(y2),maxval(y2)
  print *,'minval maxval of z = ',minval(z2),maxval(z2)

! example below taken from https://github.com/szaghi/Lib_VTK_IO
! by Stefano Zaghi

! DK DK the VTK/ParaView legacy structured grid format is used

! create and open the VTK file in STRUCTURED_GRID (BINARY) format
 E_IO = VTK_INI(output_format = 'BINARY',                     &
                filename      = 'reg_1_alpha_kernel.vtk',      &
                title         = 'Structured Grid Alpha Kernel', &
                mesh_topology = 'STRUCTURED_GRID')

! save the geometry (location of the grid points)
 E_IO = VTK_GEO(Nx=Nx_VTK,Ny=Ny_VTK,Nz=Nz_VTK,NN=NN,X=XVAL,Y=YVAL,Z=ZVAL)

! create the header to save data values at the points
 E_IO = VTK_DAT(NC_NN = NN, var_location = 'node')

! save the data defined at the nodes
 E_IO = VTK_VAR(NC_NN = NN, varname = 'Alpha_kl', var = var_str_point)

! close the file
 E_IO = VTK_END()

  end program read_interpolated_kernel

