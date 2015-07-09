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

  subroutine create_regular_elements(myrank,ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner_without_doubling,ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore, &
                    iaddx,iaddy,iaddz,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                    gammaxstore,gammaystore,gammazstore,&
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB, &
                    stretch_tab, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                    ispec_is_tiso)


! adds a regular spectral element to the different regions of the mesh

  use meshfem3D_models_par

  implicit none

  integer :: myrank,ilayer,ichunk,ispec,ipass,ifirst_region,ilast_region
  ! code for the four regions of the mesh
  integer iregion_code
  ! correct number of spectral elements in each block depending on chunk type
  integer nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS
  integer NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer :: ner_without_doubling
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

! Gauss-Lobatto-Legendre points and weights of integration
  double precision xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)

! 2D shape functions and their derivatives
  double precision dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  logical INCLUDE_CENTRAL_CUBE

! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision rmin,rmax
  double precision r_moho,r_400,r_670

! for model density and anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  integer nspec_ani
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! arrays with mesh parameters
  integer nspec_actually
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_actually) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

! Stacey, indices for Clayton-Engquist absorbing conditions
  integer nspec_stacey
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_stacey) :: rho_vp,rho_vs

! boundary locator
  logical iboun(6,nspec)

! MPI cut-planes parameters along xi and along eta
  logical, dimension(2,nspec) :: iMPIcut_xi,iMPIcut_eta

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD
  integer iproc_xi,iproc_eta

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  integer idoubling(nspec)
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical :: USE_ONE_LAYER_SB

  double precision, dimension(2,ner(1)) :: stretch_tab

! Boundary Mesh
  integer NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho
  integer ibelm_moho_top(NSPEC2D_MOHO),ibelm_moho_bot(NSPEC2D_MOHO)
  integer ibelm_400_top(NSPEC2D_400),ibelm_400_bot(NSPEC2D_400)
  integer ibelm_670_top(NSPEC2D_670),ibelm_670_bot(NSPEC2D_670)
  real(kind=CUSTOM_REAL) normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO)
  real(kind=CUSTOM_REAL) normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400)
  real(kind=CUSTOM_REAL) normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670)
  real(kind=CUSTOM_REAL) jacobian2D_moho(NGLLX,NGLLY,NSPEC2D_MOHO)
  real(kind=CUSTOM_REAL) jacobian2D_400(NGLLX,NGLLY,NSPEC2D_400)
  real(kind=CUSTOM_REAL) jacobian2D_670(NGLLX,NGLLY,NSPEC2D_670)

  integer ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot

  logical, dimension(nspec) :: ispec_is_tiso

  ! local parameters
  double precision, dimension(NGNOD) :: offset_x,offset_y,offset_z
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision :: r1,r2,r3,r4,r5,r6,r7,r8
  integer :: ix_elem,iy_elem,iz_elem,ignod,ispec_superbrick
  logical :: is_superbrick

  ! loop on all the elements
  do ix_elem = 1,NEX_PER_PROC_XI,ratio_sampling_array(ilayer)
    do iy_elem = 1,NEX_PER_PROC_ETA,ratio_sampling_array(ilayer)
      do iz_elem = 1,ner_without_doubling
        ! loop on all the corner nodes of this element
        do ignod = 1,NGNOD_EIGHT_CORNERS
          ! define topological coordinates of this mesh point
          offset_x(ignod) = (ix_elem - 1) + iaddx(ignod) * ratio_sampling_array(ilayer)
          offset_y(ignod) = (iy_elem - 1) + iaddy(ignod) * ratio_sampling_array(ilayer)
          if (ilayer == 1 .and. CASE_3D) then
            offset_z(ignod) = iaddz(ignod)
          else
            offset_z(ignod) = (iz_elem - 1) + iaddz(ignod)
          endif
        enddo
        call add_missing_nodes(offset_x,offset_y,offset_z)

        ! compute the actual position of all the grid points of that element
        if (ilayer == 1 .and. CASE_3D .and. .not. SUPPRESS_CRUSTAL_MESH) then
          ! crustal elements are stretched to be thinner in the upper crust than in lower crust in the 3D case
          ! max ratio between size of upper crust elements and
          ! lower crust elements is given by the param MAX_RATIO_STRETCHING
          ! to avoid stretching, set MAX_RATIO_STRETCHING = 1.0d  in constants.h
          call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
             ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
             NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
             stretch_tab(1,ner_without_doubling-iz_elem+1),&
             stretch_tab(2,ner_without_doubling-iz_elem+1),1,ilayer,ichunk,rotation_matrix, &
             NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)
        else
          call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
             ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
             NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
             r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer,ichunk,rotation_matrix, &
             NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)
        endif
        ! add one spectral element to the list
        ispec = ispec + 1
        if (ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')

        ! new get_flag_boundaries
        ! xmin & xmax
        if (ix_elem == 1) then
          iMPIcut_xi(1,ispec) = .true.
          if (iproc_xi == 0) iboun(1,ispec)= .true.
        endif
        if (ix_elem == (NEX_PER_PROC_XI-ratio_sampling_array(ilayer)+1)) then
          iMPIcut_xi(2,ispec) = .true.
          if (iproc_xi == NPROC_XI-1) iboun(2,ispec)= .true.
        endif
        ! ymin & ymax
        if (iy_elem == 1) then
          iMPIcut_eta(1,ispec) = .true.
          if (iproc_eta == 0) iboun(3,ispec)= .true.
        endif
        if (iy_elem == (NEX_PER_PROC_ETA-ratio_sampling_array(ilayer)+1)) then
          iMPIcut_eta(2,ispec) = .true.
          if (iproc_eta == NPROC_ETA-1) iboun(4,ispec)= .true.
        endif
        ! zmin & zmax
        if (iz_elem == ner(ilayer) .and. ilayer == ifirst_region) then
          iboun(6,ispec)= .true.
        endif
        if (iz_elem == 1 .and. ilayer == ilast_region) then    ! defined if no doubling in this layer
          iboun(5,ispec)= .true.
        endif

        ! define the doubling flag of this element
        idoubling(ispec) = doubling_index(ilayer)

        ! save the radii of the nodes before modified through compute_element_properties()
        if (ipass == 2 .and. SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
          r1=sqrt(xelm(1)*xelm(1)+yelm(1)**2+zelm(1)**2)
          r2=sqrt(xelm(2)*xelm(2)+yelm(2)**2+zelm(2)**2)
          r3=sqrt(xelm(3)*xelm(3)+yelm(3)**2+zelm(3)**2)
          r4=sqrt(xelm(4)*xelm(4)+yelm(4)**2+zelm(4)**2)
          r5=sqrt(xelm(5)*xelm(5)+yelm(5)**2+zelm(5)**2)
          r6=sqrt(xelm(6)*xelm(6)+yelm(6)**2+zelm(6)**2)
          r7=sqrt(xelm(7)*xelm(7)+yelm(7)**2+zelm(7)**2)
          r8=sqrt(xelm(8)*xelm(8)+yelm(8)**2+zelm(8)**2)
        endif

        ! compute several rheological and geometrical properties for this spectral element
        call compute_element_properties(ispec,iregion_code,idoubling,ipass, &
                         xstore,ystore,zstore,nspec,myrank, &
                         xelm,yelm,zelm,shape3D,rmin,rmax,rhostore,dvpstore, &
                         kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                         xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                         gammaxstore,gammaystore,gammazstore,nspec_actually, &
                         c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                         c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                         c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                         nspec_ani,nspec_stacey, &
                         rho_vp,rho_vs, &
                         xigll,yigll,zigll,ispec_is_tiso)

        ! boundary mesh
        if (ipass == 2 .and. SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
          is_superbrick=.false.
          ispec_superbrick = 0
          call get_jacobian_discontinuities(myrank,ispec,ix_elem,iy_elem,rmin,rmax, &
                   r1,r2,r3,r4,r5,r6,r7,r8, &
                   xstore(:,:,:,ispec),ystore(:,:,:,ispec),zstore(:,:,:,ispec),dershape2D_bottom, &
                   ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                   normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                   ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,&
                   ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                   NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,r_moho,r_400,r_670, &
                   is_superbrick,USE_ONE_LAYER_SB,ispec_superbrick,nex_eta_moho,HONOR_1D_SPHERICAL_MOHO)
        endif

      ! end of loop on all the regular elements
      enddo
    enddo
  enddo

  end subroutine create_regular_elements
