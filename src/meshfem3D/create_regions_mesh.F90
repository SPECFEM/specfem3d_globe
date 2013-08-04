!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

  subroutine create_regions_mesh(iregion_code,ibool,idoubling,is_on_a_slice_edge, &
                          xstore,ystore,zstore,rmins,rmaxs, &
                          iproc_xi,iproc_eta,ichunk,nspec,nspec_tiso, &
                          volume_local,area_local_bottom,area_local_top, &
                          nglob_theor,npointot, &
                          NSTEP,DT, &
                          NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                          NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                          NPROC_XI,NPROC_ETA,NSPEC2D_XI_FACE, &
                          NSPEC2D_ETA_FACE,NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER, &
                          myrank,LOCAL_PATH,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,&
                          SAVE_MESH_FILES,NCHUNKS,INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                          R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,&
                          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                          ner,ratio_sampling_array,doubling_index,r_bottom,r_top, &
                          this_region_has_a_doubling,ipass,ratio_divide_central_cube, &
                          CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta,USE_FULL_TISO_MANTLE, &
                          ATT1,ATT2,ATT3,SIMULATION_TYPE,USE_LDDRK,EXACT_MASS_MATRIX_FOR_ROTATION,ATTENUATION_1D_WITH_3D_STORAGE)

! creates the different regions of the mesh

  use meshfem3D_models_par

  use mpi

  implicit none

  ! this to cut the doubling brick
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER
  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_EDGES_ONEDIR) :: NSPEC2D_XI_FACE,NSPEC2D_ETA_FACE
  logical :: CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA
  integer :: offset_proc_xi,offset_proc_eta

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: ner,ratio_sampling_array
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: r_bottom,r_top
  logical, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: this_region_has_a_doubling

  integer :: ner_without_doubling,ilayer,ilayer_loop, &
               ifirst_region,ilast_region,ratio_divide_central_cube
  integer, dimension(:), allocatable :: perm_layer

  ! correct number of spectral elements in each block depending on chunk type
  integer nspec,nspec_tiso,nspec_stacey,nspec_actually,nspec_att,ATT1,ATT2,ATT3

  integer NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NCHUNKS

  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer NPROC_XI,NPROC_ETA

  integer npointot

  logical SAVE_MESH_FILES

  logical INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS,USE_FULL_TISO_MANTLE,ATTENUATION_1D_WITH_3D_STORAGE

  double precision R_CENTRAL_CUBE,RICB,RCMB,R670,RMOHO, &
          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
          RMOHO_FICTITIOUS_IN_MESHER

  double precision RHO_OCEANS

  character(len=150) LOCAL_PATH,errmsg

  ! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  ! meshing parameters
  double precision, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: rmins,rmaxs

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  ! topology of the elements
  integer, dimension(NGNOD) :: iaddx,iaddy,iaddz

  ! code for the four regions of the mesh
  integer iregion_code

  ! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

  ! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

  ! 2D shape functions and their derivatives
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y, &
    shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y, &
    dershape2D_bottom,dershape2D_top

  integer, dimension(nspec) :: idoubling

! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

  ! parameters needed to store the radii of the grid points in the spherically symmetric Earth
  double precision rmin,rmax

  ! for model density and anisotropy
  integer nspec_ani
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,dvpstore, &
    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! boundary locator
  logical, dimension(:,:), allocatable :: iboun

  ! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  ! proc numbers for MPI
  integer myrank

  ! check area and volume of the final mesh
  double precision area_local_bottom,area_local_top
  double precision volume_local

  ! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer nglob,nglob_theor,ieoff,ilocnum,ier

  ! mass matrices
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmassx,rmassy,rmassz
  integer :: nglob_xy

  ! mass matrices for backward simulation when ROTATION is .true.
  integer :: SIMULATION_TYPE
  logical :: EXACT_MASS_MATRIX_FOR_ROTATION,USE_LDDRK
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: b_rmassx,b_rmassy

  ! mass matrix and bathymetry for ocean load
  integer nglob_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

  ! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax, &
    ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

  ! 2-D jacobians and normals
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    jacobian2D_xmin,jacobian2D_xmax, &
    jacobian2D_ymin,jacobian2D_ymax,jacobian2D_bottom,jacobian2D_top

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top

  ! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi,iMPIcut_eta

  ! Stacey, indices for Clayton-Engquist absorbing conditions
  integer, dimension(:,:), allocatable :: nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rho_vp,rho_vs

  ! name of the database file
  character(len=150) prname

  ! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer i,j,k,ispec
  integer iproc_xi,iproc_eta,ichunk

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! rotation matrix from Euler angles
  double precision, dimension(NDIM,NDIM) :: rotation_matrix

  ! attenuation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: Qmu_store
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: tau_e_store
  double precision, dimension(N_SLS) :: tau_s
  double precision  T_c_source

  integer, dimension(MAX_NUMBER_OF_MESH_LAYERS) :: doubling_index
  logical :: USE_ONE_LAYER_SB

  integer NUMBER_OF_MESH_LAYERS,layer_shift,cpt, &
    first_layer_aniso,last_layer_aniso,FIRST_ELT_NON_ANISO

  double precision, dimension(:,:), allocatable :: stretch_tab

  integer :: nb_layer_above_aniso,FIRST_ELT_ABOVE_ANISO

  ! now perform two passes in this part to be able to save memory
  integer :: ipass

  logical :: ACTUALLY_STORE_ARRAYS

! added for color permutation
  integer :: nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer
  integer, dimension(:), allocatable :: perm
  integer, dimension(:), allocatable :: first_elem_number_in_this_color
  integer, dimension(:), allocatable :: num_of_elems_in_this_color

  integer :: icolor,ispec_counter
  integer :: nspec_outer_min_global,nspec_outer_max_global

  integer :: NSTEP
  integer, save :: npoin2D_xi,npoin2D_eta
  double precision :: DT

  ! boundary mesh
  integer NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho
  integer, dimension(:), allocatable :: ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
    ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: normal_moho,normal_400,normal_670
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: jacobian2D_moho,jacobian2D_400,jacobian2D_670
  integer ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,ispec2D_400_bot, &
    ispec2D_670_top,ispec2D_670_bot
  double precision r_moho,r_400,r_670

  ! flags for transverse isotropic elements
  logical, dimension(:), allocatable :: ispec_is_tiso

#ifdef USE_SERIAL_CASCADE_FOR_IOs
  logical :: you_can_start_doing_IOs
  integer msg_status(MPI_STATUS_SIZE)
#endif

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  ! new Attenuation definition on all GLL points
  ! attenuation
  if (ATTENUATION) then
    T_c_source = AM_V%QT_c_source
    tau_s(:)   = AM_V%Qtau_s(:)
    nspec_att = nspec
  else
    nspec_att = 1
  endif
  allocate(Qmu_store(ATT1,ATT2,ATT3,nspec_att), &
          tau_e_store(N_SLS,ATT1,ATT2,ATT3,nspec_att),stat=ier)
  if(ier /= 0) stop 'error in allocate 1'

  ! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX), &
          yigll(NGLLY), &
          zigll(NGLLZ),stat=ier)
  if(ier /= 0) stop 'error in allocate 2'

  ! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX), &
          wygll(NGLLY), &
          wzgll(NGLLZ),stat=ier)
  if(ier /= 0) stop 'error in allocate 3'

  ! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ), &
          dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ),stat=ier)
  if(ier /= 0) stop 'error in allocate 4'

  ! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ), &
          shape2D_y(NGNOD2D,NGLLX,NGLLZ), &
          shape2D_bottom(NGNOD2D,NGLLX,NGLLY), &
          shape2D_top(NGNOD2D,NGLLX,NGLLY), &
          dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ), &
          dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ), &
          dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY), &
          dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY),stat=ier)
  if(ier /= 0) stop 'error in allocate 5'

  ! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec), &
          dvpstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 6'

  ! for anisotropy
  allocate(kappavstore(NGLLX,NGLLY,NGLLZ,nspec), &
          muvstore(NGLLX,NGLLY,NGLLZ,nspec), &
          kappahstore(NGLLX,NGLLY,NGLLZ,nspec), &
          muhstore(NGLLX,NGLLY,NGLLZ,nspec), &
          eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec), &
          ispec_is_tiso(nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 7'

  ! initializes flags for transverse isotropic elements
  ispec_is_tiso(:) = .false.

  ! Stacey absorbing boundaries
  if(NCHUNKS /= 6) then
    nspec_stacey = nspec
  else
    nspec_stacey = 1
  endif
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey), &
          rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey),stat=ier)
  if(ier /= 0) stop 'error in allocate 8'

  ! anisotropy
  if((ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) .or. &
     (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) then
    nspec_ani = nspec
  else
    nspec_ani = 1
  endif
  allocate(c11store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c12store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c13store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c14store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c15store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c16store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c22store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c23store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c24store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c25store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c26store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c33store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c34store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c35store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c36store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c44store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c45store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c46store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c55store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c56store(NGLLX,NGLLY,NGLLZ,nspec_ani), &
          c66store(NGLLX,NGLLY,NGLLZ,nspec_ani),stat=ier)
  if(ier /= 0) stop 'error in allocate 9'

  ! boundary locator
  allocate(iboun(6,nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 10'

  ! boundary parameters locator
  allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX), &
          ibelm_xmax(NSPEC2DMAX_XMIN_XMAX), &
          ibelm_ymin(NSPEC2DMAX_YMIN_YMAX), &
          ibelm_ymax(NSPEC2DMAX_YMIN_YMAX), &
          ibelm_bottom(NSPEC2D_BOTTOM), &
          ibelm_top(NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error in allocate 11'

  ! 2-D jacobians and normals
  allocate(jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM), &
          jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error in allocate 12'

  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX), &
          normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX), &
          normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM), &
          normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP),stat=ier)
  if(ier /= 0) stop 'error in allocate 13'

  ! Stacey
  allocate(nimin(2,NSPEC2DMAX_YMIN_YMAX), &
          nimax(2,NSPEC2DMAX_YMIN_YMAX), &
          njmin(2,NSPEC2DMAX_XMIN_XMAX), &
          njmax(2,NSPEC2DMAX_XMIN_XMAX), &
          nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX), &
          nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX),stat=ier)
  if(ier /= 0) stop 'error in allocate 14'

  ! MPI cut-planes parameters along xi and along eta
  allocate(iMPIcut_xi(2,nspec), &
          iMPIcut_eta(2,nspec),stat=ier)
  if(ier /= 0) stop 'error in allocate 15'

  ! store and save the final arrays only in the second pass
  ! therefore in the first pass some arrays can be allocated with a dummy size
  if(ipass == 1) then
    ACTUALLY_STORE_ARRAYS = .false.
    nspec_actually = 1
  else
    ACTUALLY_STORE_ARRAYS = .true.
    nspec_actually = nspec
  endif
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          xiystore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          xizstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          etaxstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          etaystore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          etazstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          gammaxstore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          gammaystore(NGLLX,NGLLY,NGLLZ,nspec_actually), &
          gammazstore(NGLLX,NGLLY,NGLLZ,nspec_actually),stat=ier)
  if(ier /= 0) stop 'error in allocate 16'

  ! boundary mesh
  if (ipass == 2 .and. SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
    NSPEC2D_MOHO = NSPEC2D_TOP
    NSPEC2D_400 = NSPEC2D_MOHO / 4
    NSPEC2D_670 = NSPEC2D_400
  else
    NSPEC2D_MOHO = 1
    NSPEC2D_400 = 1
    NSPEC2D_670 = 1
  endif
  allocate(ibelm_moho_top(NSPEC2D_MOHO),ibelm_moho_bot(NSPEC2D_MOHO), &
          ibelm_400_top(NSPEC2D_400),ibelm_400_bot(NSPEC2D_400), &
          ibelm_670_top(NSPEC2D_670),ibelm_670_bot(NSPEC2D_670), &
          normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO), &
          normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400), &
          normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670), &
          jacobian2D_moho(NGLLX,NGLLY,NSPEC2D_MOHO), &
          jacobian2D_400(NGLLX,NGLLY,NSPEC2D_400), &
          jacobian2D_670(NGLLX,NGLLY,NSPEC2D_670),stat=ier)
  if(ier /= 0) stop 'error in allocate 17'

  ! initialize number of layers
  call crm_initialize_layers(myrank,ipass,xigll,yigll,zigll,wxgll,wygll,wzgll, &
                        shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        iaddx,iaddy,iaddz,nspec,xstore,ystore,zstore,ibool,idoubling, &
                        iboun,iMPIcut_xi,iMPIcut_eta,ispec2D_moho_top,ispec2D_moho_bot, &
                        ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                        NEX_PER_PROC_ETA,nex_eta_moho,RMOHO,R400,R670,r_moho,r_400,r_670, &
                        ONE_CRUST,NUMBER_OF_MESH_LAYERS,layer_shift, &
                        iregion_code,ifirst_region,ilast_region, &
                        first_layer_aniso,last_layer_aniso,nb_layer_above_aniso,is_on_a_slice_edge)

  ! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  allocate (perm_layer(ifirst_region:ilast_region),stat=ier)
  if(ier /= 0) stop 'error in allocate 18'
  perm_layer = (/ (i, i=ilast_region,ifirst_region,-1) /)

  if(iregion_code == IREGION_CRUST_MANTLE) then
    cpt=3
    perm_layer(1)=first_layer_aniso
    perm_layer(2)=last_layer_aniso
    do i = ilast_region,ifirst_region,-1
      if (i/=first_layer_aniso .and. i/=last_layer_aniso) then
        perm_layer(cpt) = i
        cpt=cpt+1
      endif
    enddo
  endif

  ! crustal layer stretching: element layer's top and bottom radii will get stretched when in crust
  ! (number of element layers in crust can vary for different resolutions and 1chunk simulations)
  allocate(stretch_tab(2,ner(1)),stat=ier)
  if(ier /= 0) stop 'error in allocate 19'
  if (CASE_3D .and. iregion_code == IREGION_CRUST_MANTLE .and. .not. SUPPRESS_CRUSTAL_MESH) then
    ! stretching function determines top and bottom of each element layer in the
    ! crust region (between r_top(1) and r_bottom(1)), where ner(1) is the
    ! number of element layers in this crust region

    ! differentiate between regional meshes or global meshes
    if( REGIONAL_MOHO_MESH ) then
      call stretching_function_regional(r_top(1),r_bottom(1),ner(1),stretch_tab)
    else
      call stretching_function(r_top(1),r_bottom(1),ner(1),stretch_tab)
    endif

    ! RMIDDLE_CRUST so far is only used for 1D - models with two layers in the crust
    ! (i.e. ONE_CRUST is set to .false.), those models do not use CASE_3D

    ! all 3D models use this stretching function to honor a 3D crustal model
    ! for those models, we set RMIDDLE_CRUST to the bottom of the first element layer
    ! this value will be used in moho_stretching.f90 to decide whether or not elements
    ! have to be stretched under oceanic crust.
    !
    ! note: stretch_tab uses (dimensionalized) radii from r_top and r_bottom
    !(with stretch_tab( index_radius(1=top,2=bottom), index_layer( 1=first layer, 2=second layer, 3= ...) )
    RMIDDLE_CRUST = stretch_tab(2,1)
  endif

!----
!----  creates mesh elements
!----

  ! loop on all the layers in this region of the mesh
  ispec = 0 ! counts all the elements in this region of the mesh
  do ilayer_loop = ifirst_region,ilast_region

    ilayer = perm_layer(ilayer_loop)

    ! determine the radii that define the shell
    rmin = rmins(ilayer)
    rmax = rmaxs(ilayer)

    if(iregion_code == IREGION_CRUST_MANTLE .and. ilayer_loop==3) then
      FIRST_ELT_NON_ANISO = ispec+1
    endif
    if(iregion_code == IREGION_CRUST_MANTLE &
      .and. ilayer_loop==(ilast_region-nb_layer_above_aniso+1)) then
      FIRST_ELT_ABOVE_ANISO = ispec+1
    endif

    ner_without_doubling = ner(ilayer)

    ! if there is a doubling at the top of this region, we implement it in the last two layers of elements
    ! and therefore we suppress two layers of regular elements here
    USE_ONE_LAYER_SB = .false.
    if(this_region_has_a_doubling(ilayer)) then
      if (ner(ilayer) == 1) then
        ner_without_doubling = ner_without_doubling - 1
        USE_ONE_LAYER_SB = .true.
      else
        ner_without_doubling = ner_without_doubling - 2
        USE_ONE_LAYER_SB = .false.
      endif
    endif

    ! regular mesh elements
    call create_regular_elements(myrank,ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner_without_doubling,ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore, &
                    iaddx,iaddy,iaddz,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                    RICB,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,RTOPDDOUBLEPRIME, &
                    R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                    gammaxstore,gammaystore,gammazstore,&
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    nspec_att,Qmu_store,tau_e_store,tau_s,T_c_source, &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB, &
                    stretch_tab,ACTUALLY_STORE_ARRAYS, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,&
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot,&
                    ispec_is_tiso,USE_FULL_TISO_MANTLE,ATT1,ATT2,ATT3,ATTENUATION_1D_WITH_3D_STORAGE)

    ! mesh doubling elements
    if( this_region_has_a_doubling(ilayer) ) &
      call create_doubling_elements(myrank,ilayer,ichunk,ispec,ipass, &
                    ifirst_region,ilast_region,iregion_code, &
                    nspec,NCHUNKS,NUMBER_OF_MESH_LAYERS, &
                    NPROC_XI,NPROC_ETA,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                    ner,ratio_sampling_array,r_top,r_bottom, &
                    xstore,ystore,zstore,xigll,yigll,zigll, &
                    shape3D,dershape2D_bottom, &
                    INCLUDE_CENTRAL_CUBE,ABSORBING_CONDITIONS, &
                    RICB,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,RTOPDDOUBLEPRIME, &
                    R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                    rmin,rmax,r_moho,r_400,r_670, &
                    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    nspec_actually,xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                    gammaxstore,gammaystore,gammazstore,&
                    nspec_stacey,rho_vp,rho_vs,iboun,iMPIcut_xi,iMPIcut_eta, &
                    ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,iproc_xi,iproc_eta, &
                    nspec_att,Qmu_store,tau_e_store,tau_s,T_c_source, &
                    rotation_matrix,idoubling,doubling_index,USE_ONE_LAYER_SB,ACTUALLY_STORE_ARRAYS, &
                    NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,nex_eta_moho, &
                    ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot,ibelm_670_top,ibelm_670_bot, &
                    normal_moho,normal_400,normal_670,jacobian2D_moho,jacobian2D_400,jacobian2D_670, &
                    ispec2D_moho_top,ispec2D_moho_bot,ispec2D_400_top,&
                    ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                    CUT_SUPERBRICK_XI,CUT_SUPERBRICK_ETA,offset_proc_xi,offset_proc_eta, &
                    ispec_is_tiso,USE_FULL_TISO_MANTLE,ATT1,ATT2,ATT3,ATTENUATION_1D_WITH_3D_STORAGE)

  enddo ! of ilayer_loop

  ! define central cube in inner core
  if(INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) &
    call create_central_cube(myrank,ichunk,ispec,iaddx,iaddy,iaddz, &
                        nspec,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,R_CENTRAL_CUBE, &
                        iproc_xi,iproc_eta,NPROC_XI,NPROC_ETA,ratio_divide_central_cube, &
                        iMPIcut_xi,iMPIcut_eta,iboun, &
                        idoubling,iregion_code,xstore,ystore,zstore, &
                        RICB,RCMB,R670,RMOHO,RMOHO_FICTITIOUS_IN_MESHER,RTOPDDOUBLEPRIME,&
                        R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
                        shape3D,rmin,rmax,rhostore,dvpstore,&
                        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                        xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
                        gammaxstore,gammaystore,gammazstore,nspec_actually, &
                        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                        c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                        nspec_ani,nspec_stacey,nspec_att,Qmu_store,tau_e_store,tau_s,T_c_source,&
                        rho_vp,rho_vs,ABSORBING_CONDITIONS,ACTUALLY_STORE_ARRAYS,xigll,yigll,zigll, &
                        ispec_is_tiso,USE_FULL_TISO_MANTLE,ATT1,ATT2,ATT3,ATTENUATION_1D_WITH_3D_STORAGE)

  ! check total number of spectral elements created
  if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

  ! if any of these flags is true, the element is on a communication edge
  ! this is not enough because it can also be in contact by an edge or a corner but not a full face
  ! therefore we will have to fix array "is_on_a_slice_edge" later in the solver to take this into account
  is_on_a_slice_edge(:) = &
      iMPIcut_xi(1,:) .or. iMPIcut_xi(2,:) .or. &
      iMPIcut_eta(1,:) .or. iMPIcut_eta(2,:) .or. &
      iboun(1,:) .or. iboun(2,:) .or. &
      iboun(3,:) .or. iboun(4,:) .or. &
      iboun(5,:) .or. iboun(6,:)

  ! no need to count fictitious elements on the edges
  ! for which communications cannot be overlapped with calculations
  where(idoubling == IFLAG_IN_FICTITIOUS_CUBE) is_on_a_slice_edge = .false.

  ! only create global addressing and the MPI buffers in the first pass
  if(ipass == 1) then

    !uncomment: adds model smoothing for point profile models
    !    if( THREE_D_MODEL == THREE_D_MODEL_PPM ) then
    !     call smooth_model(myrank, nproc_xi,nproc_eta,&
    !        rho_vp,rho_vs,nspec_stacey, &
    !        iregion_code,xixstore,xiystore,xizstore, &
    !        etaxstore,etaystore,etazstore, &
    !        gammaxstore,gammaystore,gammazstore, &
    !        xstore,ystore,zstore,rhostore,dvpstore, &
    !        kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,&
    !        nspec,HETEROGEN_3D_MANTLE, &
    !        NEX_XI,NCHUNKS,ABSORBING_CONDITIONS,PPM_V )

    ! allocate memory for arrays
    allocate(locval(npointot), &
            ifseg(npointot), &
            xp(npointot), &
            yp(npointot), &
            zp(npointot),stat=ier)
    if(ier /= 0) stop 'error in allocate 20'

    locval = 0
    ifseg = .false.
    xp = 0.d0
    yp = 0.d0
    zp = 0.d0

    ! we need to create a copy of the x, y and z arrays because sorting in get_global will swap
    ! these arrays and therefore destroy them
    do ispec=1,nspec
      ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
      ilocnum = 0
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            ilocnum = ilocnum + 1
            xp(ilocnum+ieoff) = xstore(i,j,k,ispec)
            yp(ilocnum+ieoff) = ystore(i,j,k,ispec)
            zp(ilocnum+ieoff) = zstore(i,j,k,ispec)
          enddo
        enddo
      enddo
    enddo

    call get_global(nspec,xp,yp,zp,ibool,locval,ifseg,nglob,npointot)

    deallocate(xp,stat=ier); if(ier /= 0) stop 'error in deallocate'
    deallocate(yp,stat=ier); if(ier /= 0) stop 'error in deallocate'
    deallocate(zp,stat=ier); if(ier /= 0) stop 'error in deallocate'

    ! check that number of points found equals theoretical value
    if(nglob /= nglob_theor) then
      write(errmsg,*) 'incorrect total number of points found: myrank,nglob,nglob_theor,ipass,iregion_code = ',&
        myrank,nglob,nglob_theor,ipass,iregion_code
      call exit_MPI(myrank,errmsg)
    endif
    if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering')

    ! creates a new indirect addressing to reduce cache misses in memory access in the solver
    ! this is *critical* to improve performance in the solver
    call get_global_indirect_addressing(nspec,nglob,ibool)

    ! checks again
    if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering after sorting')

    ! create MPI buffers
    ! arrays locval(npointot) and ifseg(npointot) used to save memory
    call get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                    xstore,ystore,zstore,ifseg,npointot, &
                    NSPEC2D_ETA_FACE,iregion_code,npoin2D_xi)

    call get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                    xstore,ystore,zstore,ifseg,npointot, &
                    NSPEC2D_XI_FACE,iregion_code,npoin2D_eta)

    call get_MPI_1D_buffers(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool,idoubling, &
                    xstore,ystore,zstore,ifseg,npointot, &
                    NSPEC1D_RADIAL_CORNER,NGLOB1D_RADIAL_CORNER,iregion_code)

    ! Stacey
    if(NCHUNKS /= 6) &
         call get_absorb(myrank,prname,iboun,nspec,nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
                         NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

    ! create AVS or DX mesh data for the slices
    if(SAVE_MESH_FILES) then
      call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot)

      call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot, &
              rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
              ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
              RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
              RMIDDLE_CRUST,ROCEAN,iregion_code)

      call write_AVS_DX_global_chunks_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot, &
              rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
              ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
              RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
              RMIDDLE_CRUST,ROCEAN,iregion_code)

      call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot, &
              rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
              ELLIPTICITY,ISOTROPIC_3D_MANTLE, &
              RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R120,R80,RMOHO, &
              RMIDDLE_CRUST,ROCEAN,iregion_code)

    endif

    deallocate(locval,stat=ier); if(ier /= 0) stop 'error in deallocate'
    deallocate(ifseg,stat=ier); if(ier /= 0) stop 'error in deallocate'

! only create mass matrix and save all the final arrays in the second pass
  else if(ipass == 2) then

    ! copy the theoretical number of points for the second pass
    nglob = nglob_theor

    ! count number of anisotropic elements in current region
    ! should be zero in all the regions except in the mantle
    ! (used only for checks in meshfem3D() routine)
    nspec_tiso = count(ispec_is_tiso(:))

  if(SORT_MESH_INNER_OUTER) then

! detection of the edges, coloring and permutation separately
  allocate(perm(nspec))

! implement mesh coloring for GPUs if needed, to create subsets of disconnected elements
! to remove dependencies and the need for atomic operations in the sum of elemental contributions in the solver
  if(USE_MESH_COLORING_GPU) then

    allocate(first_elem_number_in_this_color(MAX_NUMBER_OF_COLORS + 1))

    call get_perm_color_faster(is_on_a_slice_edge,ibool,perm,nspec,nglob, &
      nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer,first_elem_number_in_this_color,myrank)

! for the last color, the next color is fictitious and its first (fictitious) element number is nspec + 1
    first_elem_number_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements + 1) = nspec + 1

    allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements))

! save mesh coloring
    open(unit=99,file=prname(1:len_trim(prname))//'num_of_elems_in_this_color.dat',status='unknown')

! number of colors for outer elements
    write(99,*) nb_colors_outer_elements

! number of colors for inner elements
    write(99,*) nb_colors_inner_elements

! number of elements in each color
    do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
      num_of_elems_in_this_color(icolor) = first_elem_number_in_this_color(icolor+1) - first_elem_number_in_this_color(icolor)
      write(99,*) num_of_elems_in_this_color(icolor)
    enddo
    close(99)

! check that the sum of all the numbers of elements found in each color is equal
! to the total number of elements in the mesh
    if(sum(num_of_elems_in_this_color) /= nspec) then
      print *,'nspec = ',nspec
      print *,'total number of elements in all the colors of the mesh = ',sum(num_of_elems_in_this_color)
      stop 'incorrect total number of elements in all the colors of the mesh'
    endif

! check that the sum of all the numbers of elements found in each color for the outer elements is equal
! to the total number of outer elements found in the mesh
    if(sum(num_of_elems_in_this_color(1:nb_colors_outer_elements)) /= nspec_outer) then
      print *,'nspec_outer = ',nspec_outer
      print *,'total number of elements in all the colors of the mesh for outer elements = ',sum(num_of_elems_in_this_color)
      stop 'incorrect total number of elements in all the colors of the mesh for outer elements'
    endif

    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

    deallocate(first_elem_number_in_this_color)
    deallocate(num_of_elems_in_this_color)

  else

!! DK DK for regular C + MPI version for CPUs: do not use colors but nonetheless put all the outer elements
!! DK DK first in order to be able to overlap non-blocking MPI communications with calculations

!! DK DK nov 2010, for Rosa Badia / StarSs:
!! no need for mesh coloring, but need to implement inner/outer subsets for non blocking MPI for StarSs
    ispec_counter = 0
    perm(:) = 0

! first generate all the outer elements
    do ispec = 1,nspec
      if(is_on_a_slice_edge(ispec)) then
        ispec_counter = ispec_counter + 1
        perm(ispec) = ispec_counter
      endif
    enddo

! make sure we have detected some outer elements
    if(ispec_counter <= 0) stop 'fatal error: no outer elements detected!'

! store total number of outer elements
    nspec_outer = ispec_counter

! then generate all the inner elements
    do ispec = 1,nspec
      if(.not. is_on_a_slice_edge(ispec)) then
        ispec_counter = ispec_counter + 1
        perm(ispec) = ispec_counter
      endif
    enddo

! test that all the elements have been used once and only once
  if(ispec_counter /= nspec) stop 'fatal error: ispec_counter not equal to nspec'

! do basic checks
  if(minval(perm) /= 1) stop 'minval(perm) should be 1'
  if(maxval(perm) /= nspec) stop 'maxval(perm) should be nspec'

    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  endif

!! DK DK and Manh Ha, Nov 2011: added this to use the new mesher in the CUDA or C / StarSs test codes

  if (myrank == 0 .and. iregion_code == IREGION_CRUST_MANTLE) then

! write a header file for the Fortran version of the solver
    open(unit=99,file=prname(1:len_trim(prname))//'values_from_mesher_f90.h',status='unknown')
    write(99,*) 'integer, parameter :: NSPEC = ',nspec
    write(99,*) 'integer, parameter :: NGLOB = ',nglob
!!! DK DK use 1000 time steps only for the scaling tests
    write(99,*) 'integer, parameter :: NSTEP = 1000 !!!!!!!!!!! ',nstep
    write(99,*) 'real(kind=4), parameter :: deltat = ',DT
    write(99,*)
    write(99,*) 'integer, parameter ::  NGLOB2DMAX_XMIN_XMAX = ',npoin2D_xi
    write(99,*) 'integer, parameter ::  NGLOB2DMAX_YMIN_YMAX = ',npoin2D_eta
    write(99,*) 'integer, parameter ::  NGLOB2DMAX_ALL = ',max(npoin2D_xi,npoin2D_eta)
    write(99,*) 'integer, parameter ::  NPROC_XI = ',NPROC_XI
    write(99,*) 'integer, parameter ::  NPROC_ETA = ',NPROC_ETA
    write(99,*)
    write(99,*) '! element number of the source and of the station'
    write(99,*) '! after permutation of the elements by mesh coloring'
    write(99,*) '! and inner/outer set splitting in the mesher'
    write(99,*) 'integer, parameter :: NSPEC_SOURCE = ',perm(NSPEC/3)
    write(99,*) 'integer, parameter :: RANK_SOURCE = 0'
    write(99,*)
    write(99,*) 'integer, parameter :: RANK_STATION = (NPROC_XI*NPROC_ETA - 1)'
    write(99,*) 'integer, parameter :: NSPEC_STATION = ',perm(2*NSPEC/3)

! save coordinates of the seismic source
!   write(99,*) xstore(2,2,2,10);
!   write(99,*) ystore(2,2,2,10);
!   write(99,*) zstore(2,2,2,10);

! save coordinates of the seismic station
!   write(99,*) xstore(2,2,2,nspec-10);
!   write(99,*) ystore(2,2,2,nspec-10);
!   write(99,*) zstore(2,2,2,nspec-10);
    close(99)

!! write a header file for the C version of the solver
    open(unit=99,file=prname(1:len_trim(prname))//'values_from_mesher_C.h',status='unknown')
    write(99,*) '#define NSPEC ',nspec
    write(99,*) '#define NGLOB ',nglob
!!    write(99,*) '#define NSTEP ',nstep
!!! DK DK use 1000 time steps only for the scaling tests
    write(99,*) '// #define NSTEP ',nstep
    write(99,*) '#define NSTEP 1000'
! put an "f" at the end to force single precision
    write(99,"('#define deltat ',e18.10,'f')") DT
    write(99,*) '#define NGLOB2DMAX_XMIN_XMAX ',npoin2D_xi
    write(99,*) '#define NGLOB2DMAX_YMIN_YMAX ',npoin2D_eta
    write(99,*) '#define NGLOB2DMAX_ALL ',max(npoin2D_xi,npoin2D_eta)
    write(99,*) '#define NPROC_XI ',NPROC_XI
    write(99,*) '#define NPROC_ETA ',NPROC_ETA
    write(99,*)
    write(99,*) '// element and MPI slice number of the source and the station'
    write(99,*) '// after permutation of the elements by mesh coloring'
    write(99,*) '// and inner/outer set splitting in the mesher'
    write(99,*) '#define RANK_SOURCE 0'
    write(99,*) '#define NSPEC_SOURCE ',perm(NSPEC/3)
    write(99,*)
    write(99,*) '#define RANK_STATION (NPROC_XI*NPROC_ETA - 1)'
    write(99,*) '#define NSPEC_STATION ',perm(2*NSPEC/3)
    close(99)

    open(unit=99,file=prname(1:len_trim(prname))//'values_from_mesher_nspec_outer.h',status='unknown')
    write(99,*) '#define NSPEC_OUTER ',nspec_outer_max_global
    write(99,*) '// NSPEC_OUTER_min = ',nspec_outer_min_global
    write(99,*) '// NSPEC_OUTER_max = ',nspec_outer_max_global
    close(99)

  endif

!! DK DK and Manh Ha, Nov 2011: added this to use the new mesher in the CUDA or C / StarSs test codes

  deallocate(perm)

  else
!
    print *,'SORT_MESH_INNER_OUTER must always been set to .true. even for the regular C version for CPUs'
    print *,'in order to be able to use non blocking MPI to overlap communications'
!   print *,'generating identity permutation'
!   do ispec = 1,nspec
!     perm(ispec) = ispec
!   enddo
    stop 'please set SORT_MESH_INNER_OUTER to .true. and recompile the whole code'

  endif

!!!! David Michea: end of mesh coloring

!****************************************************************************************************

    ! precomputes jacobian for 2d absorbing boundary surfaces
    call get_jacobian_boundaries(myrank,iboun,nspec,xstore,ystore,zstore, &
              dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
              ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
              nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
              jacobian2D_xmin,jacobian2D_xmax, &
              jacobian2D_ymin,jacobian2D_ymax, &
              jacobian2D_bottom,jacobian2D_top, &
              normal_xmin,normal_xmax, &
              normal_ymin,normal_ymax, &
              normal_bottom,normal_top, &
              NSPEC2D_BOTTOM,NSPEC2D_TOP, &
              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,&
              xigll,yigll,zigll)

    ! allocates mass matrices in this slice (will be fully assembled in the solver)
    !
    ! in the case of Stacey boundary conditions, add C*delta/2 contribution to the mass matrix
    ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
    ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
    !
    ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
    ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be fictitious / unused
    nglob_xy = 1

    if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. (.not. USE_LDDRK)) then
       select case(iregion_code)
       case( IREGION_CRUST_MANTLE )
          nglob_xy = nglob
       case( IREGION_INNER_CORE, IREGION_OUTER_CORE )
          nglob_xy = 1
       endselect
    else
       nglob_xy = 1
    endif

    if(.not. USE_LDDRK)then
      if(ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION)then
        select case(iregion_code)
        case( IREGION_CRUST_MANTLE,IREGION_INNER_CORE )
           nglob_xy = nglob
        case( IREGION_OUTER_CORE )
           nglob_xy = 1
        endselect
      endif
    endif

    allocate(rmassx(nglob_xy),stat=ier)
    if(ier /= 0) stop 'error in allocate 21'
    allocate(rmassy(nglob_xy),stat=ier)
    if(ier /= 0) stop 'error in allocate 21'
    allocate(rmassz(nglob),stat=ier)
    if(ier /= 0) stop 'error in allocate 21'

    allocate(b_rmassx(nglob_xy),stat=ier)
    if(ier /= 0) stop 'error in allocate b_21'
    allocate(b_rmassy(nglob_xy),stat=ier)
    if(ier /= 0) stop 'error in allocate b_21'

    ! allocates ocean load mass matrix as well if oceans
    if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
      nglob_oceans = nglob
    else
      ! allocate dummy array if no oceans
      nglob_oceans = 1
    endif
    allocate(rmass_ocean_load(nglob_oceans),stat=ier)
    if(ier /= 0) stop 'error in allocate 22'

    ! creating mass matrices in this slice (will be fully assembled in the solver)
    call create_mass_matrices(myrank,nspec,idoubling,wxgll,wygll,wzgll,ibool, &
                          nspec_actually,xixstore,xiystore,xizstore, &
                          etaxstore,etaystore,etazstore, &
                          gammaxstore,gammaystore,gammazstore, &
                          iregion_code,rhostore,kappavstore, &
                          nglob_xy,nglob,prname, &
                          rmassx,rmassy,rmassz,DT,ichunk,NCHUNKS,ABSORBING_CONDITIONS, &
                          nglob_oceans,rmass_ocean_load, &
                          xstore,ystore,zstore,RHO_OCEANS, &
                          NSPEC2D_TOP,NSPEC2D_BOTTOM,NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                          ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                          nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                          normal_xmin,normal_xmax,normal_ymin,normal_ymax, &
                          rho_vp,rho_vs,nspec_stacey, &
                          jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
                          jacobian2D_bottom,jacobian2D_top,&
                          SIMULATION_TYPE,EXACT_MASS_MATRIX_FOR_ROTATION,USE_LDDRK, &
                          b_rmassx,b_rmassy)

    ! save the binary files
#ifdef USE_SERIAL_CASCADE_FOR_IOs
    you_can_start_doing_IOs = .false.
    if (myrank > 0) call MPI_RECV(you_can_start_doing_IOs, 1, MPI_LOGICAL, myrank-1, itag, MPI_COMM_WORLD, msg_status,ier)
!!!!    print *,'starting doing serialized I/Os on rank ',myrank
    print *,'starting doing serialized I/Os on rank ',myrank
#endif
    call save_arrays_solver(rho_vp,rho_vs,nspec_stacey, &
                  prname,iregion_code,xixstore,xiystore,xizstore, &
                  etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                  xstore,ystore,zstore,rhostore,dvpstore, &
                  kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                  nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                  c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                  c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                  ibool,idoubling,is_on_a_slice_edge,nglob_xy,nglob, &
                  rmassx,rmassy,rmassz,rmass_ocean_load,nglob_oceans, &
                  ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                  nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                  normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
                  jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
                  jacobian2D_bottom,jacobian2D_top,nspec, &
                  NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                  TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
                  ANISOTROPIC_INNER_CORE,OCEANS, &
                  tau_s,tau_e_store,Qmu_store,T_c_source,ATTENUATION, &
                  ATT1,ATT2,ATT3,size(tau_e_store,5),&
                  NCHUNKS,ABSORBING_CONDITIONS,SAVE_MESH_FILES,ispec_is_tiso,myrank,&
                  ROTATION,EXACT_MASS_MATRIX_FOR_ROTATION,USE_LDDRK,&
                  b_rmassx,b_rmassy)
#ifdef USE_SERIAL_CASCADE_FOR_IOs
    you_can_start_doing_IOs = .true.
    if (myrank < NPROC_XI*NPROC_ETA-1) call MPI_SEND(you_can_start_doing_IOs, 1, MPI_LOGICAL, myrank+1, itag, MPI_COMM_WORLD, ier)
#endif

    deallocate(rmassx,stat=ier); if(ier /= 0) stop 'error in deallocate rmassx'
    deallocate(rmassy,stat=ier); if(ier /= 0) stop 'error in deallocate rmassy'
    deallocate(rmassz,stat=ier); if(ier /= 0) stop 'error in deallocate rmassz'
    deallocate(b_rmassx,stat=ier); if(ier /= 0) stop 'error in deallocate b_rmassx'
    deallocate(b_rmassy,stat=ier); if(ier /= 0) stop 'error in deallocate b_rmassy'
    deallocate(rmass_ocean_load,stat=ier); if(ier /= 0) stop 'error in deallocate rmass_ocean_load'

    ! boundary mesh
    if (SAVE_BOUNDARY_MESH .and. iregion_code == IREGION_CRUST_MANTLE) then
      ! first check the number of surface elements are the same for Moho, 400, 670
      if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
        if (ispec2D_moho_top /= NSPEC2D_MOHO .or. ispec2D_moho_bot /= NSPEC2D_MOHO) &
               call exit_mpi(myrank, 'Not the same number of Moho surface elements')
      endif
      if (ispec2D_400_top /= NSPEC2D_400 .or. ispec2D_400_bot /= NSPEC2D_400) &
               call exit_mpi(myrank,'Not the same number of 400 surface elements')
      if (ispec2D_670_top /= NSPEC2D_670 .or. ispec2D_670_bot /= NSPEC2D_670) &
               call exit_mpi(myrank,'Not the same number of 670 surface elements')

      ! writing surface topology databases
      open(unit=27,file=prname(1:len_trim(prname))//'boundary_disc.bin',status='unknown',form='unformatted')
      write(27) NSPEC2D_MOHO, NSPEC2D_400, NSPEC2D_670
      write(27) ibelm_moho_top
      write(27) ibelm_moho_bot
      write(27) ibelm_400_top
      write(27) ibelm_400_bot
      write(27) ibelm_670_top
      write(27) ibelm_670_bot
      write(27) normal_moho
      write(27) normal_400
      write(27) normal_670
      close(27)
    endif

    ! compute volume, bottom and top area of that part of the slice
    call crm_compute_volumes(volume_local,area_local_bottom,area_local_top, &
                            nspec,wxgll,wygll,wzgll,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                            NSPEC2D_BOTTOM,jacobian2D_bottom,NSPEC2D_TOP,jacobian2D_top)


  else
    stop 'there cannot be more than two passes in mesh creation'

  endif  ! end of test if first or second pass

  deallocate(stretch_tab)
  deallocate(perm_layer)

  ! deallocate these arrays after each pass because they have a different size in each pass to save memory
  deallocate(xixstore,xiystore,xizstore,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(etaxstore,etaystore,etazstore,stat=ier); if(ier /= 0) stop 'error in deallocate'
  deallocate(gammaxstore,gammaystore,gammazstore,stat=ier); if(ier /= 0) stop 'error in deallocate'

  ! deallocate arrays
  deallocate(rhostore,dvpstore,kappavstore,kappahstore)
  deallocate(muvstore,muhstore)
  deallocate(eta_anisostore)
  deallocate(ispec_is_tiso)
  deallocate(c11store)
  deallocate(c12store)
  deallocate(c13store)
  deallocate(c14store)
  deallocate(c15store)
  deallocate(c16store)
  deallocate(c22store)
  deallocate(c23store)
  deallocate(c24store)
  deallocate(c25store)
  deallocate(c26store)
  deallocate(c33store)
  deallocate(c34store)
  deallocate(c35store)
  deallocate(c36store)
  deallocate(c44store)
  deallocate(c45store)
  deallocate(c46store)
  deallocate(c55store)
  deallocate(c56store)
  deallocate(c66store)
  deallocate(iboun)
  deallocate(xigll,yigll,zigll)
  deallocate(wxgll,wygll,wzgll)
  deallocate(shape3D,dershape3D)
  deallocate(shape2D_x,shape2D_y,shape2D_bottom,shape2D_top)
  deallocate(dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top)
  deallocate(ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax)
  deallocate(ibelm_bottom,ibelm_top)
  deallocate(jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax)
  deallocate(jacobian2D_bottom,jacobian2D_top)
  deallocate(normal_xmin,normal_xmax,normal_ymin,normal_ymax)
  deallocate(normal_bottom,normal_top)
  deallocate(iMPIcut_xi,iMPIcut_eta)
  deallocate(nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta)
  deallocate(rho_vp,rho_vs)
  deallocate(Qmu_store)
  deallocate(tau_e_store)
  deallocate(ibelm_moho_top,ibelm_moho_bot)
  deallocate(ibelm_400_top,ibelm_400_bot)
  deallocate(ibelm_670_top,ibelm_670_bot)
  deallocate(normal_moho,normal_400,normal_670)
  deallocate(jacobian2D_moho,jacobian2D_400,jacobian2D_670)

  end subroutine create_regions_mesh

!
!-------------------------------------------------------------------------------
!


  subroutine crm_initialize_layers(myrank,ipass,xigll,yigll,zigll,wxgll,wygll,wzgll, &
                        shape3D,dershape3D,shape2D_x,shape2D_y,shape2D_bottom,shape2D_top, &
                        dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top, &
                        iaddx,iaddy,iaddz,nspec,xstore,ystore,zstore,ibool,idoubling, &
                        iboun,iMPIcut_xi,iMPIcut_eta,ispec2D_moho_top,ispec2D_moho_bot, &
                        ispec2D_400_top,ispec2D_400_bot,ispec2D_670_top,ispec2D_670_bot, &
                        NEX_PER_PROC_ETA,nex_eta_moho,RMOHO,R400,R670,r_moho,r_400,r_670, &
                        ONE_CRUST,NUMBER_OF_MESH_LAYERS,layer_shift, &
                        iregion_code,ifirst_region,ilast_region, &
                        first_layer_aniso,last_layer_aniso,nb_layer_above_aniso,is_on_a_slice_edge)

! create the different regions of the mesh

  implicit none

  include "constants.h"

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
  integer first_layer_aniso,last_layer_aniso,nb_layer_above_aniso

! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

! if number of points is odd, the middle abscissa is exactly zero
  if(mod(NGLLX,2) /= 0) xigll((NGLLX-1)/2+1) = ZERO
  if(mod(NGLLY,2) /= 0) yigll((NGLLY-1)/2+1) = ZERO
  if(mod(NGLLZ,2) /= 0) zigll((NGLLZ-1)/2+1) = ZERO

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
  if(iregion_code == IREGION_CRUST_MANTLE) then
    ifirst_region = 1
    ilast_region = 10 + layer_shift

  else if(iregion_code == IREGION_OUTER_CORE) then
    ifirst_region = 11 + layer_shift
    ilast_region = NUMBER_OF_MESH_LAYERS - 1

  else if(iregion_code == IREGION_INNER_CORE) then
    ifirst_region = NUMBER_OF_MESH_LAYERS
    ilast_region = NUMBER_OF_MESH_LAYERS

  else
    call exit_MPI(myrank,'incorrect region code detected')
  endif

! to consider anisotropic elements first and to build the mesh from the bottom to the top of the region
  if (ONE_CRUST) then
    first_layer_aniso=2
    last_layer_aniso=3
    nb_layer_above_aniso = 1
  else
    first_layer_aniso=3
    last_layer_aniso=4
    nb_layer_above_aniso = 2
  endif

! initialize mesh arrays
  idoubling(:) = 0

  xstore(:,:,:,:) = 0.d0
  ystore(:,:,:,:) = 0.d0
  zstore(:,:,:,:) = 0.d0

  if(ipass == 1) ibool(:,:,:,:) = 0

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

  end subroutine crm_initialize_layers


!
!-------------------------------------------------------------------------------------------------
!

  subroutine crm_compute_volumes(volume_local,area_local_bottom,area_local_top, &
                            nspec,wxgll,wygll,wzgll,xixstore,xiystore,xizstore, &
                            etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                            NSPEC2D_BOTTOM,jacobian2D_bottom,NSPEC2D_TOP,jacobian2D_top)

  implicit none

  include "constants.h"

  double precision :: volume_local,area_local_bottom,area_local_top

  integer :: nspec
  double precision :: wxgll(NGLLX),wygll(NGLLY),wzgll(NGLLZ)

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  integer :: NSPEC2D_BOTTOM,NSPEC2D_TOP
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM) :: jacobian2D_bottom
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP) :: jacobian2D_top

  ! local parameters
  double precision :: weight
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  integer :: i,j,k,ispec

  ! initializes
  volume_local = ZERO
  area_local_bottom = ZERO
  area_local_top = ZERO

  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          weight = wxgll(i)*wygll(j)*wzgll(k)

          ! compute the jacobian
          xixl = xixstore(i,j,k,ispec)
          xiyl = xiystore(i,j,k,ispec)
          xizl = xizstore(i,j,k,ispec)
          etaxl = etaxstore(i,j,k,ispec)
          etayl = etaystore(i,j,k,ispec)
          etazl = etazstore(i,j,k,ispec)
          gammaxl = gammaxstore(i,j,k,ispec)
          gammayl = gammaystore(i,j,k,ispec)
          gammazl = gammazstore(i,j,k,ispec)

          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          volume_local = volume_local + dble(jacobianl)*weight

        enddo
      enddo
    enddo
  enddo

  do ispec = 1,NSPEC2D_BOTTOM
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_bottom = area_local_bottom + dble(jacobian2D_bottom(i,j,ispec))*weight
      enddo
    enddo
  enddo

  do ispec = 1,NSPEC2D_TOP
    do i=1,NGLLX
      do j=1,NGLLY
        weight=wxgll(i)*wygll(j)
        area_local_top = area_local_top + dble(jacobian2D_top(i,j,ispec))*weight
      enddo
    enddo
  enddo


  end subroutine crm_compute_volumes

