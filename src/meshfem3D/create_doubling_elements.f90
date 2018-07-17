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

  subroutine create_doubling_elements(ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore, &
                    gammaxstore,gammaystore,gammazstore, &
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                    CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta, &
                    ispec_is_tiso)


! adds doubling elements to the different regions of the mesh

  use constants
  use meshfem3D_models_par, only: myrank,HONOR_1D_SPHERICAL_MOHO

  implicit none

  integer,intent(in) :: ilayer,ichunk,ipass,ifirst_region,ilast_region
  integer,intent(inout) :: ispec
  ! code for the four regions of the mesh
  integer,intent(in) :: iregion_code
  ! correct number of spectral elements in each block depending on chunk type
  integer,intent(in) :: nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS
  integer,intent(in) :: NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS),intent(in) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS),intent(in) :: r_bottom,r_top

! arrays with the mesh in double precision
  double precision,intent(inout) :: xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision,intent(inout) :: ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision,intent(inout) :: zstore(NGLLX,NGLLY,NGLLZ,nspec)

! Gauss-Lobatto-Legendre points and weights of integration
  double precision,intent(in) :: xigll(NGLLX),yigll(NGLLY),zigll(NGLLZ)

! 3D shape functions and their derivatives
  double precision,intent(in) :: shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)

! 2D shape functions and their derivatives
  double precision,intent(in) :: dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY)

  logical,intent(in) :: INCLUDE_CENTRAL_CUBE

! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision,intent(in) :: rmin,rmax
  double precision,intent(in) :: r_moho,r_400,r_670

! for model density and anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  integer,intent(in) :: nspec_ani
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani),intent(inout) :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! arrays with mesh parameters
  integer,intent(in) :: nspec_actually
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_actually),intent(inout) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

! Stacey, indices for Clayton-Engquist absorbing conditions
  integer,intent(in) :: nspec_stacey
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_stacey),intent(inout) :: rho_vp,rho_vs

! boundary locator
  logical,intent(inout) :: iboun(6,nspec)

! MPI cut-planes parameters along xi and along eta
  logical, dimension(2,nspec),intent(inout) :: iMPIcut_xi,iMPIcut_eta

  double precision,intent(in) :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD
  integer,intent(in) :: iproc_xi,iproc_eta

! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM),intent(in) :: rotation_matrix

  integer,intent(inout) :: idoubling(nspec)
  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS),intent(in) :: doubling_index
  logical,intent(in) :: USE_ONE_LAYER_SB

! Boundary Mesh
  integer :: NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho
  integer :: ibelm_moho_top(NSPEC2D_MOHO),ibelm_moho_bot(NSPEC2D_MOHO)
  integer :: ibelm_400_top(NSPEC2D_400),ibelm_400_bot(NSPEC2D_400)
  integer :: ibelm_670_top(NSPEC2D_670),ibelm_670_bot(NSPEC2D_670)
  real(kind=CUSTOM_REAL) :: normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO)
  real(kind=CUSTOM_REAL) :: normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400)
  real(kind=CUSTOM_REAL) :: normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670)
  real(kind=CUSTOM_REAL) :: jacobian2D_moho(NGLLX,NGLLY,NSPEC2D_MOHO)
  real(kind=CUSTOM_REAL) :: jacobian2D_400(NGLLX,NGLLY,NSPEC2D_400)
  real(kind=CUSTOM_REAL) :: jacobian2D_670(NGLLX,NGLLY,NSPEC2D_670)

  integer :: ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot

  integer,intent(in) :: offset_proc_xi,offset_proc_eta
  logical,intent(in) :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA

  logical, dimension(nspec),intent(inout) :: ispec_is_tiso

  ! local parameters
  double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick
  double precision, dimension(NGNOD) :: offset_x,offset_y,offset_z
  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision :: r1,r2,r3,r4,r5,r6,r7,r8
  ! mesh doubling superbrick
  integer, dimension(NGNOD_EIGHT_CORNERS,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick
  integer :: ix_elem,iy_elem,iz_elem,ignod,ispec_superbrick,case_xi,case_eta
  integer :: step_mult,subblock_num
  integer :: nspec_sb
  logical, dimension(NSPEC_DOUBLING_SUPERBRICK,6) :: iboun_sb
  logical :: is_superbrick

  integer, dimension(:), allocatable :: map_ispec
  integer :: nelements,ispec0,ielem,ier,ispec_loc

! If there is a doubling at the top of this region, let us add these elements.
! The superbrick implements a symmetric four-to-two doubling and therefore replaces
! a basic regular block of 2 x 2 = 4 elements.
! We have imposed that NEX be a multiple of 16 therefore we know that we can always create
! these 2 x 2 blocks because NEX_PER_PROC_XI / ratio_sampling_array(ilayer) and
! NEX_PER_PROC_ETA / ratio_sampling_array(ilayer) are always divisible by 2.

  if (USE_ONE_LAYER_SB) then
    call define_superbrick_one_layer(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
    nspec_sb = NSPEC_SUPERBRICK_1L
    iz_elem = ner(ilayer)
    step_mult = 2
  else
    if (iregion_code == IREGION_OUTER_CORE .and. ilayer == ilast_region &
      .and. (CUT_SUPERBRICK_XI .or. CUT_SUPERBRICK_ETA)) then
      nspec_sb = NSPEC_DOUBLING_BASICBRICK
      step_mult = 1
    else
      call define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick,iboun_sb)
      nspec_sb = NSPEC_DOUBLING_SUPERBRICK
      step_mult = 2
    endif
    ! the doubling is implemented in the last two radial elements
    ! therefore we start one element before the last one
    iz_elem = ner(ilayer) - 1
  endif

  ! stores original value
  ispec0 = ispec

  ! counts number of elements for this layer
  nelements = NEX_PER_PROC_XI/(step_mult*ratio_sampling_array(ilayer)) &
            * NEX_PER_PROC_ETA/(step_mult*ratio_sampling_array(ilayer)) &
            * nspec_sb

  ! fill mapping to be able to parallelize loops below
  allocate(map_ispec(nelements),stat=ier)
  if (ier /= 0) stop 'Error allocating map_ispec'
  map_ispec(:) = 0

  do ix_elem = 1,NEX_PER_PROC_XI,step_mult*ratio_sampling_array(ilayer)
    do iy_elem = 1,NEX_PER_PROC_ETA,step_mult*ratio_sampling_array(ilayer)
      do ispec_superbrick = 1,nspec_sb
        ! counts in increasing order 1,2,3,..
        ielem = ispec_superbrick + nspec_sb * ( (iy_elem-1)/(step_mult*ratio_sampling_array(ilayer)) &
                                            +   (ix_elem-1)/(step_mult*ratio_sampling_array(ilayer)) &
                                                * NEX_PER_PROC_ETA/(step_mult*ratio_sampling_array(ilayer)))
        map_ispec(ielem) = ispec0 + ielem
        ! check
        if (map_ispec(ielem) > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
      enddo
    enddo
  enddo

  ! loop on all the elements in the 2 x 2 blocks

!daniel: still debugging...
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(ix_elem,iy_elem,ispec_superbrick,ignod, &
!$OMP offset_x,offset_y,offset_z,xelm,yelm,zelm, &
!$OMP case_xi,case_eta,subblock_num, &
!$OMP r1,r2,r3,r4,r5,r6,r7,r8, &
!$OMP is_superbrick,ispec_loc,ielem)
!$OMP DO
  do ix_elem = 1,NEX_PER_PROC_XI,step_mult*ratio_sampling_array(ilayer)
    do iy_elem = 1,NEX_PER_PROC_ETA,step_mult*ratio_sampling_array(ilayer)

      if (step_mult == 1) then
        ! for xi direction
        if (.not. CUT_SUPERBRICK_XI) then
          if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer))) == 0) then
            case_xi = 1
          else
            case_xi = 2
          endif
        else
          if (offset_proc_xi == 0) then
            if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer))) == 0) then
              case_xi = 1
            else
              case_xi = 2
            endif
          else
            if (mod((ix_elem-1),(2*step_mult*ratio_sampling_array(ilayer))) /= 0) then
              case_xi = 1
            else
              case_xi = 2
            endif
          endif
        endif
        ! for eta direction
        if (.not. CUT_SUPERBRICK_ETA) then
          if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer))) == 0) then
            case_eta = 1
          else
            case_eta = 2
          endif
        else
          if (offset_proc_eta == 0) then
            if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer))) == 0) then
              case_eta = 1
            else
              case_eta = 2
            endif
          else
            if (mod((iy_elem-1),(2*step_mult*ratio_sampling_array(ilayer))) /= 0) then
              case_eta = 1
            else
              case_eta = 2
            endif
          endif
        endif
        ! determine the current sub-block
        if (case_xi == 1) then
          if (case_eta == 1) then
            subblock_num = 1
          else
            subblock_num = 2
          endif
        else
          if (case_eta == 1) then
            subblock_num = 3
          else
            subblock_num = 4
          endif
        endif
        ! then define the geometry for this sub-block
        call define_basic_doubling_brick(x_superbrick,y_superbrick, &
                        z_superbrick,ibool_superbrick,iboun_sb,subblock_num)
      endif

      ! loop on all the elements in the mesh doubling superbrick
      do ispec_superbrick = 1,nspec_sb
        ! loop on all the corner nodes of this element
        do ignod = 1,NGNOD_EIGHT_CORNERS
          ! define topological coordinates of this mesh point
          offset_x(ignod) = (ix_elem - 1) + &
            x_superbrick(ibool_superbrick(ignod,ispec_superbrick)) * ratio_sampling_array(ilayer)
          offset_y(ignod) = (iy_elem - 1) + &
            y_superbrick(ibool_superbrick(ignod,ispec_superbrick)) * ratio_sampling_array(ilayer)
          offset_z(ignod) = (iz_elem - 1) + &
            z_superbrick(ibool_superbrick(ignod,ispec_superbrick))
        enddo
        ! the rest of the 27 nodes are missing, therefore add them
        call add_missing_nodes(offset_x,offset_y,offset_z)

        ! compute the actual position of all the grid points of that element
        call compute_coord_main_mesh(offset_x,offset_y,offset_z,xelm,yelm,zelm, &
           ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
           NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
           r_top(ilayer),r_bottom(ilayer),ner(ilayer),ilayer,ichunk,rotation_matrix, &
           NCHUNKS,INCLUDE_CENTRAL_CUBE,NUMBER_OF_MESH_LAYERS)

        ! add one spectral element to the list
        !ispec = ispec + 1

        ! counts in increasing order 1,2,3,..
        ielem = ispec_superbrick + nspec_sb * ( (iy_elem-1)/(step_mult*ratio_sampling_array(ilayer)) &
                                            +   (ix_elem-1)/(step_mult*ratio_sampling_array(ilayer)) &
                                                * NEX_PER_PROC_ETA/(step_mult*ratio_sampling_array(ilayer)))
        ispec_loc = map_ispec(ielem)
        if (ispec_loc > nspec .or. ispec_loc < 1) call exit_MPI(myrank,'invalid ispec_loc in mesh creation')

        ! new get_flag_boundaries
        ! xmin & xmax
        if (ix_elem == 1) then
            iMPIcut_xi(1,ispec_loc) = iboun_sb(ispec_superbrick,1)
            if (iproc_xi == 0) iboun(1,ispec_loc)= iboun_sb(ispec_superbrick,1)
        endif
        if (ix_elem == (NEX_PER_PROC_XI-step_mult*ratio_sampling_array(ilayer)+1)) then
            iMPIcut_xi(2,ispec_loc) = iboun_sb(ispec_superbrick,2)
            if (iproc_xi == NPROC_XI-1) iboun(2,ispec_loc)= iboun_sb(ispec_superbrick,2)
        endif
        !! ymin & ymax
        if (iy_elem == 1) then
            iMPIcut_eta(1,ispec_loc) = iboun_sb(ispec_superbrick,3)
            if (iproc_eta == 0) iboun(3,ispec_loc)= iboun_sb(ispec_superbrick,3)
        endif
        if (iy_elem == (NEX_PER_PROC_ETA-step_mult*ratio_sampling_array(ilayer)+1)) then
            iMPIcut_eta(2,ispec_loc) = iboun_sb(ispec_superbrick,4)
            if (iproc_eta == NPROC_ETA-1) iboun(4,ispec_loc)= iboun_sb(ispec_superbrick,4)
        endif
        ! zmax only
        if (ilayer == ifirst_region) then
          iboun(6,ispec_loc)= iboun_sb(ispec_superbrick,6)
        endif
        if (ilayer == ilast_region .and. iz_elem == 1) then
          iboun(5,ispec_loc)= iboun_sb(ispec_superbrick,5)
        endif

        ! define the doubling flag of this element
        idoubling(ispec_loc) = doubling_index(ilayer)

        ! save the radii of the nodes before modified through compute_element_properties()
        if (ipass == 2 .and. SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
          r1 = sqrt(xelm(1)*xelm(1) + yelm(1)*yelm(1) + zelm(1)*zelm(1))
          r2 = sqrt(xelm(2)*xelm(2) + yelm(2)*yelm(2) + zelm(2)*zelm(2))
          r3 = sqrt(xelm(3)*xelm(3) + yelm(3)*yelm(3) + zelm(3)*zelm(3))
          r4 = sqrt(xelm(4)*xelm(4) + yelm(4)*yelm(4) + zelm(4)*zelm(4))
          r5 = sqrt(xelm(5)*xelm(5) + yelm(5)*yelm(5) + zelm(5)*zelm(5))
          r6 = sqrt(xelm(6)*xelm(6) + yelm(6)*yelm(6) + zelm(6)*zelm(6))
          r7 = sqrt(xelm(7)*xelm(7) + yelm(7)*yelm(7) + zelm(7)*zelm(7))
          r8 = sqrt(xelm(8)*xelm(8) + yelm(8)*yelm(8) + zelm(8)*zelm(8))
        endif

        ! compute several rheological and geometrical properties for this spectral element
        call compute_element_properties(ispec_loc,iregion_code,idoubling,ipass, &
                         xstore,ystore,zstore,nspec, &
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
          is_superbrick = .true.
          call get_jacobian_discontinuities(ispec_loc,ix_elem,iy_elem,rmin,rmax, &
              r1,r2,r3,r4,r5,r6,r7,r8, &
              xstore(:,:,:,ispec_loc),ystore(:,:,:,ispec_loc),zstore(:,:,:,ispec_loc),dershape2D_bottom, &
              ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
              normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
              ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top, &
              ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
              NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,r_moho,r_400,r_670, &
              is_superbrick,USE_ONE_LAYER_SB,ispec_superbrick,nex_eta_moho,HONOR_1D_SPHERICAL_MOHO)
        endif

      ! end of loops on the mesh doubling elements
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  ! end index
  ispec = ispec0 + nelements

  ! free array
  deallocate(map_ispec)

  end subroutine create_doubling_elements
