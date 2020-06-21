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


module postprocess_par

  use constants, only: CUSTOM_REAL,MAX_STRING_LEN, &
    NGLLX,NGLLY,NGLLZ,IIN,IOUT, &
    FOUR_THIRDS,GAUSSALPHA,GAUSSBETA, &
    myrank

  use shared_parameters, only: LOCAL_PATH,R_PLANET_KM

  implicit none

  ! array dimensions for static compilation
  include "OUTPUT_FILES/values_from_mesher.h"

  ! maximum number of kernel names (comma-separated e.g. vsv,vsh,vpv,vph,eta,rho -> 6 kernel names)
  integer,parameter :: MAX_KERNEL_NAMES = 24

  ! maximum number of kernel directory paths (e.g. for summing kernels from different events)
  integer,parameter :: MAX_KERNEL_PATHS = 65535

  ! mesh size
  integer :: NSPEC, NGLOB

  ! volume
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x, y, z
  integer, dimension(:,:,:,:),allocatable :: ibool
  logical, dimension(:),allocatable :: ispec_is_tiso

end module postprocess_par

