!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  subroutine initialize_layers(NEX_PER_PROC_ETA,nex_eta_moho,RMOHO,R400,R670,r_moho,r_400,r_670, &
                               NUMBER_OF_MESH_LAYERS, &
                               iregion_code,ifirst_region,ilast_region, &
                               first_layer_aniso,last_layer_aniso)

! create the different regions of the mesh

  use constants, only: myrank,IREGION_CRUST_MANTLE,IREGION_INNER_CORE,IREGION_OUTER_CORE

  use shared_parameters, only: R_PLANET,ONE_CRUST

  implicit none

  integer,intent(in) :: NEX_PER_PROC_ETA
  integer,intent(out) :: nex_eta_moho
  double precision,intent(in) :: RMOHO,R400,R670
  double precision,intent(out) :: r_moho,r_400,r_670

  integer,intent(out) :: NUMBER_OF_MESH_LAYERS

  ! code for the four regions of the mesh
  integer,intent(in) :: iregion_code
  integer,intent(out) :: ifirst_region,ilast_region
  integer,intent(out) :: first_layer_aniso,last_layer_aniso

  ! local parameters
  integer :: layer_offset

  ! sets number of layers
  call define_all_layers_number_and_offset(NUMBER_OF_MESH_LAYERS,layer_offset)

  ! define the first and last layers that define this region
  if (iregion_code == IREGION_CRUST_MANTLE) then
    ifirst_region = 1
    ilast_region = 10 + layer_offset

  else if (iregion_code == IREGION_OUTER_CORE) then
    ifirst_region = 11 + layer_offset
    ilast_region = NUMBER_OF_MESH_LAYERS - 1

  else if (iregion_code == IREGION_INNER_CORE) then
    ifirst_region = NUMBER_OF_MESH_LAYERS
    ilast_region = NUMBER_OF_MESH_LAYERS

  else
    call exit_MPI(myrank,'incorrect region code detected')
  endif

  ! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  ! note: in older versions, we assumed to have anisotropic elements at the beginning of ibool(..),etc. arrays
  !       to save memory and directly address only a portion of it for anisotropic stress computation.
  !       in newer versions, this mesh addressing has become more general and doesn't rely anymore of having
  !       anisotropic elements at the beginning of these arrays.
  !       still, it doesn't hurt to keep this and have it in case for backward compatibility.
  if (ONE_CRUST) then
    first_layer_aniso = 2         ! layer 80-MOHO   (see define_all_layers.f90)
    last_layer_aniso = 3          ! layer 220-80
  else
    first_layer_aniso = 3         ! layer 80-MOHO
    last_layer_aniso = 4          ! layer 220-80
  endif

  nex_eta_moho = NEX_PER_PROC_ETA

  r_moho = RMOHO / R_PLANET
  r_400 = R400 / R_PLANET
  r_670 = R670 / R_PLANET

  end subroutine initialize_layers
