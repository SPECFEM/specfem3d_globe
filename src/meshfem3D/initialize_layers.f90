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

  subroutine initialize_layers(myrank,ipass,xigll,yigll,zigll,wxgll,wygll,wzgll, &
                        shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        iaddx,iaddy,iaddz,nspec,xstore,ystore,zstore,ibool,idoubling, &
                        iboun,iMPIcut_xi,iMPIcut_eta,ispec2D_moho_top,ispec2D_moho_bot, &
                        ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                        NEX_PER_PROC_ETA,nex_eta_moho,RMOHO,R400,R670,r_moho,r_400,r_670, &
                        ONE_CRUST,NUMBER_OF_MESH_LAYERS,layer_shift, &
                        iregion_code,ifirst_region,ilast_region, &
                        first_layer_aniso,last_layer_aniso,is_on_a_slice_edge)

! create the different regions of the mesh

  use constants

  implicit none

  integer :: myrank,ipass

  double precision xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ)
  double precision wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ),dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ)

  double precision shape2D_x(NGNOD2D,NGLLY,NGLLZ),shape2D_y(NGNOD2D,NGLLX,NGLLZ)
  double precision shape2D_bottom(NGNOD2D,NGLLX,NGLLY),shape2D_top(NGNOD2D,NGLLX,NGLLY)
  double precision dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ),dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ)
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY),dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

  integer nspec
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)
  integer idoubling(nspec)

  logical iboun(6,nspec)
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  integer ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot
  integer NEX_PER_PROC_ETA,nex_eta_moho
  double precision RMOHO,R400,R670
  double precision r_moho,r_400,r_670

  logical ONE_CRUST
  integer NUMBER_OF_MESH_LAYERS,layer_shift

  ! code for the four regions of the mesh
  integer iregion_code,ifirst_region,ilast_region
  integer first_layer_aniso,last_layer_aniso

! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! get the 3-D shape functions
  call get_shape3D(myrank,shape3D,dershape3D,xigll,yigll,zigll)

! get the 2-D shape functions
  call get_shape2D(myrank,shape2D_x,dershape2D_x,yigll,zigll,NGLLY,NGLLZ)
  call get_shape2D(myrank,shape2D_y,dershape2D_y,xigll,zigll,NGLLX,NGLLZ)
  call get_shape2D(myrank,shape2D_bottom,dershape2D_bottom,xigll,yigll,NGLLX,NGLLY)
  call get_shape2D(myrank,shape2D_top,dershape2D_top,xigll,yigll,NGLLX,NGLLY)

! create the shape of the corner nodes of a regular mesh element
  call hex_nodes(iaddx,iaddy,iaddz)

! reference element has size one here, not two
  iaddx(:) = iaddx(:) / 2
  iaddy(:) = iaddy(:) / 2
  iaddz(:) = iaddz(:) / 2

! sets number of layers
  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
    layer_shift = 0
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
    layer_shift = 1
  endif

  if (.not. ADD_4TH_DOUBLING) NUMBER_OF_MESH_LAYERS = NUMBER_OF_MESH_LAYERS - 1

! define the first and last layers that define this region
  if (iregion_code == IREGION_CRUST_MANTLE) then
    ifirst_region = 1
    ilast_region = 10 + layer_shift

  else if (iregion_code == IREGION_OUTER_CORE) then
    ifirst_region = 11 + layer_shift
    ilast_region = NUMBER_OF_MESH_LAYERS - 1

  else if (iregion_code == IREGION_INNER_CORE) then
    ifirst_region = NUMBER_OF_MESH_LAYERS
    ilast_region = NUMBER_OF_MESH_LAYERS

  else
    call exit_MPI(myrank,'incorrect region code detected')
  endif

! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  if (ONE_CRUST) then
    first_layer_aniso = 2
    last_layer_aniso = 3
  else
    first_layer_aniso = 3
    last_layer_aniso = 4
  endif

! initialize mesh arrays
  idoubling(:) = 0

  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  if (ipass == 1) ibool(:,:,:,:) = 0

  ! initialize boundary arrays
  iboun(:,:) = .false.
  iMPIcut_xi(:,:) = .false.
  iMPIcut_eta(:,:) = .false.
  is_on_a_slice_edge(:) = .false.

  ! boundary mesh
  ispec2D_moho_top = 0; ispec2D_moho_bot = 0
  ispec2D_400_top = 0; ispec2D_400_bot = 0
  ispec2D_670_top = 0; ispec2D_670_bot = 0

  nex_eta_moho = NEX_PER_PROC_ETA

  r_moho = RMOHO/R_EARTH; r_400 = R400 / R_EARTH; r_670 = R670/R_EARTH

  end subroutine initialize_layers
