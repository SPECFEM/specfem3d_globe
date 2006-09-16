!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 6
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!       (c) California Institute of Technology September 2006
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine create_regions_mesh(iregion_code,xgrid,ygrid,zgrid,ibool,idoubling, &
           xstore,ystore,zstore,npx,npy,rmins,rmaxs, &
           iproc_xi,iproc_eta,ichunk,nspec,nspec_tiso, &
           volume_local,area_local_bottom,area_local_top, &
           nspl,rspl,espl,espl2, &
           nglob_AB,nglob_AC,nglob_BC,npointot, &
           NER,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
           NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400,NER_400_220, &
           NER_220_MOHO,NER_CRUST,NER_DOUBLING_OUTER_CORE, &
           NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
           ELLIPTICITY,TOPOGRAPHY,TRANSVERSE_ISOTROPY, &
           ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
           NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
           NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA,NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
           myrank,LOCAL_PATH,OCEANS,ibathy_topo,NER_ICB_BOTTOMDBL, &
           crustal_model,mantle_model,aniso_mantle_model, &
           aniso_inner_core_model,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,&
           attenuation_model,ATTENUATION,ATTENUATION_3D,SAVE_MESH_FILES, &
           NCHUNKS,INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL, &
           R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN)

! create the different regions of the mesh

  implicit none

  include "constants.h"

  external mantle_model,crustal_model,aniso_mantle_model, &
       aniso_inner_core_model,attenuation_model

! correct number of spectral elements in each block depending on chunk type
  integer nspec,nspec_tiso,nspec_stacey

  integer NER,NEX_XI,NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NCHUNKS,REFERENCE_1D_MODEL
  integer NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_670_400,NER_400_220,NER_220_MOHO
  integer NER_CRUST,NER_CMB_670,NER_DOUBLING_OUTER_CORE,NER_ICB_BOTTOMDBL

  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

  integer NPROC_XI,NPROC_ETA,NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI
  integer NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA,NSPEC1D_RADIAL,NPOIN1D_RADIAL

  integer npx,npy
  integer npointot

  logical ELLIPTICITY,TOPOGRAPHY
  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,OCEANS

  logical ATTENUATION,ATTENUATION_3D,SAVE_MESH_FILES, &
          INCLUDE_CENTRAL_CUBE,INFLATE_CENTRAL_CUBE,ABSORBING_CONDITIONS

  double precision R_CENTRAL_CUBE,RICB,RHO_OCEANS,RCMB,R670,RMOHO, &
          RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN

  character(len=150) LOCAL_PATH

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

! arrays with the mesh in double precision
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  double precision, dimension(0:2*NEX_PER_PROC_XI,0:2*NEX_PER_PROC_ETA,0:2*NER) :: xgrid,ygrid,zgrid

! to define the central cube in the inner core
  integer nx_central_cube,ny_central_cube,nz_central_cube
  double precision ratio_x,ratio_y,ratio_z,radius_cube
  double precision xi_central_cube,eta_central_cube,gamma_central_cube
  double precision, dimension(:,:,:), allocatable :: xgrid_central_cube,ygrid_central_cube,zgrid_central_cube

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! meshing parameters
  double precision rmins(0:2*NER),rmaxs(0:2*NER)

! auxiliary variables to generate the mesh
  integer ix,iy,iz,ir,ir1,ir2,dir
  integer ix1,ix2,dix,iy1,iy2,diy
  integer iax,iay,iar
  integer isubregion,nsubregions,doubling_index

! code for the four regions of the mesh
  integer iregion_code

! Gauss-Lobatto-Legendre points and weights of integration
  double precision, dimension(:), allocatable :: xigll,yigll,zigll,wxgll,wygll,wzgll

! 3D shape functions and their derivatives
  double precision, dimension(:,:,:,:), allocatable :: shape3D
  double precision, dimension(:,:,:,:,:), allocatable :: dershape3D

! 2D shape functions and their derivatives
  double precision, dimension(:,:,:), allocatable :: shape2D_x,shape2D_y,shape2D_bottom,shape2D_top
  double precision, dimension(:,:,:,:), allocatable :: dershape2D_x,dershape2D_y,dershape2D_bottom,dershape2D_top

! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

! topology of the elements
  integer iaddx(NGNOD)
  integer iaddy(NGNOD)
  integer iaddz(NGNOD)

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

! parameters needed to store the radii of the grid points
! in the spherically symmetric Earth
  integer idoubling(nspec)
  double precision rmin,rmax

! for model density and anisotropy
  integer nspec_ani
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store

! the jacobian
  real(kind=CUSTOM_REAL) jacobianl

! boundary locator
  logical, dimension(:,:), allocatable :: iboun

! arrays with mesh parameters
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore

! proc numbers for MPI
  integer myrank

! check area and volume of the final mesh
  double precision weight
  double precision area_local_bottom,area_local_top
  double precision volume_local

! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer, dimension(:), allocatable :: iglob,locval
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  integer nglob,nglob_theor,nglob_AB,nglob_AC,nglob_BC

  integer ieoff,ilocnum
  integer iboolmin,iboolmax

! mass matrix
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass

! mass matrix and bathymetry for ocean load
  integer ix_oceans,iy_oceans,iz_oceans,ispec_oceans
  integer ispec2D_top_crust
  integer nglob_oceans
  double precision xval,yval,zval,rval,thetaval,phival
  double precision lat,lon,colat
  double precision elevation,height_oceans
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_ocean_load

! boundary parameters locator
  integer, dimension(:), allocatable :: ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top

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

  integer i,j,k,ia,ispec,iglobnum
  integer iproc_xi,iproc_eta,ichunk

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! rotation matrix from Euler angles
  double precision rotation_matrix(3,3)

! attenuation
  double precision, dimension(:,:,:,:),   allocatable :: Qmu_store
  double precision, dimension(:,:,:,:,:), allocatable :: tau_e_store
  double precision, dimension(N_SLS)                  :: tau_s
  double precision  T_c_source

! **************

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

! Attenuation
  if(ATTENUATION .and. ATTENUATION_3D) then
    allocate(Qmu_store(NGLLX,NGLLY,NGLLZ,nspec))
    allocate(tau_e_store(N_SLS,NGLLX,NGLLY,NGLLZ,nspec))
  else
    allocate(Qmu_store(1,1,1,1))
    allocate(tau_e_store(N_SLS,1,1,1,1))
    Qmu_store(1,1,1,1) = 0.0d0
    tau_e_store(:,1,1,1,1) = 0.0d0
  endif

  T_c_source = 0.0
  tau_s(:) = 0.0

! Gauss-Lobatto-Legendre points of integration
  allocate(xigll(NGLLX))
  allocate(yigll(NGLLY))
  allocate(zigll(NGLLZ))

! Gauss-Lobatto-Legendre weights of integration
  allocate(wxgll(NGLLX))
  allocate(wygll(NGLLY))
  allocate(wzgll(NGLLZ))

! 3D shape functions and their derivatives
  allocate(shape3D(NGNOD,NGLLX,NGLLY,NGLLZ))
  allocate(dershape3D(NDIM,NGNOD,NGLLX,NGLLY,NGLLZ))

! 2D shape functions and their derivatives
  allocate(shape2D_x(NGNOD2D,NGLLY,NGLLZ))
  allocate(shape2D_y(NGNOD2D,NGLLX,NGLLZ))
  allocate(shape2D_bottom(NGNOD2D,NGLLX,NGLLY))
  allocate(shape2D_top(NGNOD2D,NGLLX,NGLLY))
  allocate(dershape2D_x(NDIM2D,NGNOD2D,NGLLY,NGLLZ))
  allocate(dershape2D_y(NDIM2D,NGNOD2D,NGLLX,NGLLZ))
  allocate(dershape2D_bottom(NDIM2D,NGNOD2D,NGLLX,NGLLY))
  allocate(dershape2D_top(NDIM2D,NGNOD2D,NGLLX,NGLLY))

! array with model density
  allocate(rhostore(NGLLX,NGLLY,NGLLZ,nspec))

! for anisotropy
  allocate(kappavstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(muvstore(NGLLX,NGLLY,NGLLZ,nspec))

  allocate(kappahstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(muhstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec))

! Stacey
  if(NCHUNKS /= 6) then
    nspec_stacey = nspec
  else
    nspec_stacey = 1
  endif
  allocate(rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey))
  allocate(rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey))

  nspec_ani = 1
  if((ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) .or. &
     (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) nspec_ani = nspec

  allocate(c11store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c12store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c13store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c14store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c15store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c16store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c22store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c23store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c24store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c25store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c26store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c33store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c34store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c35store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c36store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c44store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c45store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c46store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c55store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c56store(NGLLX,NGLLY,NGLLZ,nspec_ani))
  allocate(c66store(NGLLX,NGLLY,NGLLZ,nspec_ani))

! boundary locator
  allocate(iboun(6,nspec))

! arrays with mesh parameters
  allocate(xixstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xiystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(xizstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(etazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaxstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammaystore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(gammazstore(NGLLX,NGLLY,NGLLZ,nspec))
  allocate(jacobianstore(NGLLX,NGLLY,NGLLZ,nspec))

! boundary parameters locator
  allocate(ibelm_xmin(NSPEC2DMAX_XMIN_XMAX))
  allocate(ibelm_xmax(NSPEC2DMAX_XMIN_XMAX))
  allocate(ibelm_ymin(NSPEC2DMAX_YMIN_YMAX))
  allocate(ibelm_ymax(NSPEC2DMAX_YMIN_YMAX))
  allocate(ibelm_bottom(NSPEC2D_BOTTOM))
  allocate(ibelm_top(NSPEC2D_TOP))

! 2-D jacobians and normals
  allocate(jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM))
  allocate(jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP))

  allocate(normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX))
  allocate(normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX))
  allocate(normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM))
  allocate(normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP))

! Stacey
  allocate(nimin(2,NSPEC2DMAX_YMIN_YMAX))
  allocate(nimax(2,NSPEC2DMAX_YMIN_YMAX))
  allocate(njmin(2,NSPEC2DMAX_XMIN_XMAX))
  allocate(njmax(2,NSPEC2DMAX_XMIN_XMAX))
  allocate(nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX))
  allocate(nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX))

! MPI cut-planes parameters along xi and along eta
  allocate(iMPIcut_xi(2,nspec))
  allocate(iMPIcut_eta(2,nspec))

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

! allocate memory for arrays
  allocate(iglob(npointot))
  allocate(locval(npointot))
  allocate(ifseg(npointot))
  allocate(xp(npointot))
  allocate(yp(npointot))
  allocate(zp(npointot))

! define models 1066a and ak135
  call define_model_1066a(CRUSTAL)
  call define_model_ak135(CRUSTAL)

! generate the elements in all the regions of the mesh
  ispec = 0

! define number of subregions in current region
  select case(iregion_code)
    case(IREGION_CRUST_MANTLE)
      nsubregions = 32
    case(IREGION_OUTER_CORE)
      nsubregions = 19
    case(IREGION_INNER_CORE)
      nsubregions = 3
    case default
      call exit_MPI(myrank,'incorrect region code')
  end select

  do isubregion = 1,nsubregions

! define shape of elements in current subregion
  select case(iregion_code)

    case(IREGION_CRUST_MANTLE)
      call define_subregions_crust_mantle(isubregion,ichunk,iaddx,iaddy,iaddz, &
        ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
        doubling_index,npx,npy,NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400, &
        NER_400_220,NER_220_MOHO,NER_CRUST)

    case(IREGION_OUTER_CORE)
      call define_subregions_outer_core(isubregion,ichunk,iaddx,iaddy,iaddz, &
        ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
        doubling_index,npx,npy,NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_DOUBLING_OUTER_CORE,NER_ICB_BOTTOMDBL)

    case(IREGION_INNER_CORE)
      call define_subregions_inner_core(isubregion,iaddx,iaddy,iaddz, &
        ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
        doubling_index,npx,npy,NER_TOP_CENTRAL_CUBE_ICB)

    case default
      call exit_MPI(myrank,'incorrect region code')

  end select

  do ir = ir1,ir2,dir
    do iy = iy1,iy2,diy
      do ix = ix1,ix2,dix

!       determine the radii that define the shell
        rmin = rmins(ir)
        rmax = rmaxs(ir)

!       loop over the NGNOD nodes
        do ia=1,NGNOD
          xelm(ia) = xgrid(ix+iax*iaddx(ia),iy+iay*iaddy(ia),ir+iar*iaddz(ia))
          yelm(ia) = ygrid(ix+iax*iaddx(ia),iy+iay*iaddy(ia),ir+iar*iaddz(ia))
          zelm(ia) = zgrid(ix+iax*iaddx(ia),iy+iay*iaddy(ia),ir+iar*iaddz(ia))
        enddo

! add one spectral element to the list and store its material number
        ispec = ispec + 1
        if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in mesh creation')
        idoubling(ispec) = doubling_index

! compute values for the Earth model
        call get_model(myrank,iregion_code,nspec,iproc_xi,iproc_eta, &
          kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,rhostore, &
          nspec_ani, &
          c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
          c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
          c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
          xelm,yelm,zelm,shape3D,ispec,iboun,iMPIcut_xi,iMPIcut_eta, &
          rmin,rmax,ichunk,doubling_index, &
          rho_vp,rho_vs,nspec_stacey, &
          NPROC_XI,NPROC_ETA, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
          crustal_model,mantle_model,aniso_mantle_model, &
          aniso_inner_core_model,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
          attenuation_model, ATTENUATION, ATTENUATION_3D, tau_s, tau_e_store, Qmu_store, T_c_source, &
          size(tau_e_store,2), size(tau_e_store,3), size(tau_e_store,4), size(tau_e_store,5), &
          NCHUNKS,INFLATE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL, &
          R_CENTRAL_CUBE,RCMB,RICB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN)

! add topography without the crustal model
        if(TOPOGRAPHY .and. (idoubling(ispec) == IFLAG_CRUST &
                        .or. idoubling(ispec) == IFLAG_220_MOHO)) &
          call add_topography(myrank,xelm,yelm,zelm,ibathy_topo,R220)

! make the Earth elliptical
        if(ELLIPTICITY) call get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)

! recompute coordinates and jacobian for real 3-D model
        call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
               etaxstore,etaystore,etazstore, &
               gammaxstore,gammaystore,gammazstore,jacobianstore, &
               xstore,ystore,zstore, &
               xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)
      enddo
    enddo
  enddo

! end of loop on all the subregions of the current region the mesh
  enddo

!---

! define central cube in inner core

  if(INCLUDE_CENTRAL_CUBE .and. iregion_code == IREGION_INNER_CORE) then

! elements are always regular in inner core
  call usual_hex_nodes(iaddx,iaddy,iaddz)

! define vertical slice in central cube on current processor
! we can assume that NEX_XI = NEX_ETA, otherwise central cube cannot be defined
  nx_central_cube = NEX_PER_PROC_XI / 8
  ny_central_cube = NEX_PER_PROC_ETA / 8
  nz_central_cube = NEX_XI / 8

! allocate grid of points in central cube
  allocate(xgrid_central_cube(0:2*nx_central_cube,0:2*ny_central_cube,0:2*nz_central_cube))
  allocate(ygrid_central_cube(0:2*nx_central_cube,0:2*ny_central_cube,0:2*nz_central_cube))
  allocate(zgrid_central_cube(0:2*nx_central_cube,0:2*ny_central_cube,0:2*nz_central_cube))

! size of the cube along cartesian axes before rotation
  radius_cube = R_CENTRAL_CUBE / dsqrt(3.d0)

  do iz = 0,2*nz_central_cube
    do iy = 0,2*ny_central_cube
      do ix = 0,2*nx_central_cube

! the slice extends to the entire cube along Z
! but only to current block along X and Y
      ratio_x = (dble(iproc_xi) + dble(ix)/dble(2*nx_central_cube)) / dble(NPROC_XI)
      ratio_y = (dble(iproc_eta) + dble(iy)/dble(2*ny_central_cube)) / dble(NPROC_ETA)
      ratio_z = dble(iz)/dble(2*nz_central_cube)

      if(dabs(ratio_x) > 1.001d0 .or. dabs(ratio_y) > 1.001d0 .or. dabs(ratio_z) > 1.001d0) &
        call exit_MPI(myrank,'wrong ratio in central cube')

! use a "flat" cubed sphere to create the central cube

! define angle in [-pi/4,pi/4]
      xi_central_cube = - PI_OVER_FOUR*(1.d0 - ratio_x) + PI_OVER_FOUR*ratio_x
      eta_central_cube = - PI_OVER_FOUR*(1.d0 - ratio_y) + PI_OVER_FOUR*ratio_y
      gamma_central_cube = - PI_OVER_FOUR*(1.d0 - ratio_z) + PI_OVER_FOUR*ratio_z

! take tangent of angle to map to [-1,1] and then map to real radius
      xgrid_central_cube(ix,iy,iz) = radius_cube * tan(xi_central_cube)
      ygrid_central_cube(ix,iy,iz) = radius_cube * tan(eta_central_cube)
      zgrid_central_cube(ix,iy,iz) = radius_cube * tan(gamma_central_cube)

      enddo
    enddo
  enddo

! define spectral elements in central cube
  do iz = 0,2*nz_central_cube-2,2
    do iy = 0,2*ny_central_cube-2,2
      do ix = 0,2*nx_central_cube-2,2

!       radii that define the shell, we know that we are in the central cube
        rmin = 0.d0
        rmax = RICB / R_EARTH

!       loop over the NGNOD nodes
        do ia=1,NGNOD

! flat cubed sphere with correct mapping
        if(ichunk == CHUNK_AB) then
          xelm(ia) = -ygrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          yelm(ia) = xgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          zelm(ia) = zgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))

        else if(ichunk == CHUNK_AB_ANTIPODE) then
          xelm(ia) = -ygrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          yelm(ia) = -xgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          zelm(ia) = -zgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))

        else if(ichunk == CHUNK_AC) then
          xelm(ia) = -ygrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          yelm(ia) = -zgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          zelm(ia) = xgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))

        else if(ichunk == CHUNK_AC_ANTIPODE) then
          xelm(ia) = -ygrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          yelm(ia) = zgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          zelm(ia) = -xgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))

        else if(ichunk == CHUNK_BC) then
          xelm(ia) = -zgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          yelm(ia) = ygrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          zelm(ia) = xgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))

        else if(ichunk == CHUNK_BC_ANTIPODE) then
          xelm(ia) = zgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          yelm(ia) = -ygrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))
          zelm(ia) = xgrid_central_cube(ix+iaddx(ia),iy+iaddy(ia),iz+iaddz(ia))

        else
          call exit_MPI(myrank,'wrong chunk number in flat cubed sphere definition')
        endif

        enddo

! add one spectral element to the list and store its material number
        ispec = ispec + 1
        if(ispec > nspec) call exit_MPI(myrank,'ispec greater than nspec in central cube creation')

! only one active central cube, the 5 others are fictitious
        if(ichunk == CHUNK_AB) then
          if(iz == 0) then
            doubling_index = IFLAG_BOTTOM_CENTRAL_CUBE
          else if(iz == 2*nz_central_cube-2) then
            doubling_index = IFLAG_TOP_CENTRAL_CUBE
          else
            doubling_index = IFLAG_IN_CENTRAL_CUBE
          endif
        else
          doubling_index = IFLAG_IN_FICTITIOUS_CUBE
        endif
        idoubling(ispec) = doubling_index

! compute values for the Earth model
        call get_model(myrank,iregion_code,nspec,iproc_xi,iproc_eta, &
          kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,rhostore, &
          nspec_ani, &
          c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
          c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
          c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
          xelm,yelm,zelm,shape3D,ispec,iboun,iMPIcut_xi,iMPIcut_eta, &
          rmin,rmax,ichunk,doubling_index, &
          rho_vp,rho_vs,nspec_stacey, &
          NPROC_XI,NPROC_ETA, &
          TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE, &
          ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST, &
          crustal_model,mantle_model,aniso_mantle_model, &
          aniso_inner_core_model,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
          attenuation_model, ATTENUATION, ATTENUATION_3D, tau_s, tau_e_store, Qmu_store, T_c_source, &
          size(tau_e_store,2), size(tau_e_store,3), size(tau_e_store,4), size(tau_e_store,5), &
          NCHUNKS,INFLATE_CENTRAL_CUBE,ABSORBING_CONDITIONS,REFERENCE_1D_MODEL, &
          R_CENTRAL_CUBE,RCMB,RICB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN)

! make the Earth elliptical
        if(ELLIPTICITY) call get_ellipticity(xelm,yelm,zelm,nspl,rspl,espl,espl2)

! recompute coordinates and jacobian for real 3-D model
        call calc_jacobian(myrank,xixstore,xiystore,xizstore, &
               etaxstore,etaystore,etazstore, &
               gammaxstore,gammaystore,gammazstore,jacobianstore, &
               xstore,ystore,zstore, &
               xelm,yelm,zelm,shape3D,dershape3D,ispec,nspec)
      enddo
    enddo
  enddo

! deallocate grid of points in central cube
  deallocate(xgrid_central_cube)
  deallocate(ygrid_central_cube)
  deallocate(zgrid_central_cube)

  endif    ! end of definition of central cube in inner core

!---

! check total number of spectral elements created
  if(ispec /= nspec) call exit_MPI(myrank,'ispec should equal nspec')

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

  call get_global(nspec,xp,yp,zp,iglob,locval,ifseg,nglob,npointot)

! theoretical total number of points we should find
  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    nglob_theor = nglob_AB
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    nglob_theor = nglob_AC
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    nglob_theor = nglob_BC
  else
    call exit_MPI(myrank,'wrong chunk type')
  endif

! check that number of points found equals theoretical value
  if(nglob /= nglob_theor) call exit_MPI(myrank,'incorrect total number of points found')

! put in classical format
  do ispec=1,nspec
  ieoff = NGLLX * NGLLY * NGLLZ * (ispec-1)
  ilocnum = 0
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        ilocnum = ilocnum + 1
        ibool(i,j,k,ispec) = iglob(ilocnum+ieoff)
      enddo
    enddo
  enddo
  enddo

  iboolmin = minval(ibool(:,:,:,1:nspec))
  iboolmax = maxval(ibool(:,:,:,1:nspec))

  if(iboolmin /= 1 .or. iboolmax /= nglob_theor) call exit_MPI(myrank,'incorrect global numbering')

! count number of anisotropic elements in current region
! should be zero in all the regions except in the mantle
  nspec_tiso = count(idoubling(1:nspec) == IFLAG_220_MOHO)

! creating mass matrix in this slice (will be fully assembled in the solver)
  allocate(rmass(nglob))
  rmass(:) = 0._CUSTOM_REAL

  do ispec=1,nspec

! suppress fictitious elements in central cube
  if(idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        weight=wxgll(i)*wygll(j)*wzgll(k)
        iglobnum=ibool(i,j,k,ispec)

        jacobianl=jacobianstore(i,j,k,ispec)

! definition depends if region is fluid or solid
  if(iregion_code == IREGION_CRUST_MANTLE .or. iregion_code == IREGION_INNER_CORE) then

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(rhostore(i,j,k,ispec)) * dble(jacobianl) * weight)
    else
      rmass(iglobnum) = rmass(iglobnum) + rhostore(i,j,k,ispec) * jacobianl * weight
    endif

! fluid in outer core
  else if(iregion_code == IREGION_OUTER_CORE) then

! no anisotropy in the fluid, use kappav

! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      rmass(iglobnum) = rmass(iglobnum) + &
             sngl(dble(jacobianl) * weight * dble(rhostore(i,j,k,ispec)) / dble(kappavstore(i,j,k,ispec)))
    else
      rmass(iglobnum) = rmass(iglobnum) + &
             jacobianl * weight * rhostore(i,j,k,ispec) / kappavstore(i,j,k,ispec)
    endif

  else
    call exit_MPI(myrank,'wrong region code')
  endif

      enddo
    enddo
  enddo
  enddo

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
              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)

! create MPI buffers
! arrays locval(npointot) and ifseg(npointot) used to save memory
  call get_MPI_cutplanes_xi(myrank,prname,nspec,iMPIcut_xi,ibool, &
                  xstore,ystore,zstore,ifseg,ichunk,npointot, &
                  NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA)
  call get_MPI_cutplanes_eta(myrank,prname,nspec,iMPIcut_eta,ibool, &
                  xstore,ystore,zstore,ifseg,ichunk,npointot, &
                  NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI)
  call get_MPI_1D_buffers(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool,idoubling, &
                  xstore,ystore,zstore,ifseg,npointot, &
                  NSPEC1D_RADIAL,NPOIN1D_RADIAL)

! Stacey
  if(NCHUNKS /= 6) &
       call get_absorb(myrank,prname,iboun,nspec,nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
                       NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! create AVS or DX mesh data for the slices
  if(SAVE_MESH_FILES) then
    call write_AVS_DX_global_data(myrank,prname,nspec,ibool,idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
    call write_AVS_DX_mesh_quality_data(prname,nspec,xstore,ystore,zstore)
    call write_AVS_DX_global_faces_data(myrank,prname,nspec,iMPIcut_xi,iMPIcut_eta,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
    call write_AVS_DX_global_chunks_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot, &
              rhostore,kappavstore,muvstore,nspl,rspl,espl,espl2, &
              ELLIPTICITY,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST,REFERENCE_1D_MODEL, &
              RICB,RCMB,RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R80,RMOHO, &
              RMIDDLE_CRUST,ROCEAN)
    call write_AVS_DX_surface_data(myrank,prname,nspec,iboun,ibool, &
              idoubling,xstore,ystore,zstore,locval,ifseg,npointot)
  endif

! save the binary files
! save ocean load mass matrix as well if oceans
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then

! adding ocean load mass matrix at the top of the crust for oceans
  nglob_oceans = nglob
  allocate(rmass_ocean_load(nglob_oceans))

! create ocean load mass matrix for degrees of freedom at ocean bottom
  rmass_ocean_load(:) = 0._CUSTOM_REAL

! add contribution of the oceans
! for surface elements exactly at the top of the crust (ocean bottom)
    do ispec2D_top_crust = 1,NSPEC2D_TOP

      ispec_oceans = ibelm_top(ispec2D_top_crust)

      iz_oceans = NGLLZ

      do ix_oceans = 1,NGLLX
        do iy_oceans = 1,NGLLY

        iglobnum=ibool(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! compute local height of oceans
        if(ISOTROPIC_3D_MANTLE) then

! get coordinates of current point
          xval = xstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          yval = ystore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)
          zval = zstore(ix_oceans,iy_oceans,iz_oceans,ispec_oceans)

! map to latitude and longitude for bathymetry routine
          call xyz_2_rthetaphi_dble(xval,yval,zval,rval,thetaval,phival)
          call reduce(thetaval,phival)

! convert the geocentric colatitude to a geographic colatitude
          colat = PI/2.0d0 - datan(1.006760466d0*dcos(thetaval)/dmax1(TINYVAL,dsin(thetaval)))

! get geographic latitude and longitude in degrees
          lat = 90.0d0 - colat*180.0d0/PI
          lon = phival*180.0d0/PI
          elevation = 0.d0

! compute elevation at current point
          call get_topo_bathy(lat,lon,elevation,ibathy_topo)

! non-dimensionalize the elevation, which is in meters
! and suppress positive elevation, which means no oceans
          if(elevation >= - MINIMUM_THICKNESS_3D_OCEANS) then
            height_oceans = 0.d0
          else
            height_oceans = dabs(elevation) / R_EARTH
          endif

        else
          height_oceans = THICKNESS_OCEANS_PREM
        endif

! take into account inertia of water column
        weight = wxgll(ix_oceans)*wygll(iy_oceans)*dble(jacobian2D_top(ix_oceans,iy_oceans,ispec2D_top_crust)) &
                   * dble(RHO_OCEANS) * height_oceans

! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + sngl(weight)
        else
          rmass_ocean_load(iglobnum) = rmass_ocean_load(iglobnum) + weight
        endif

        enddo
      enddo

    enddo

! add regular mass matrix to ocean load contribution
  rmass_ocean_load(:) = rmass_ocean_load(:) + rmass(:)

  else

! allocate dummy array if no oceans
    nglob_oceans = 1
    allocate(rmass_ocean_load(nglob_oceans))

  endif

! save the binary files
    call save_arrays_solver(rho_vp,rho_vs,nspec_stacey, &
            prname,iregion_code,xixstore,xiystore,xizstore, &
            etaxstore,etaystore,etazstore, &
            gammaxstore,gammaystore,gammazstore,jacobianstore, &
            xstore,ystore,zstore, rhostore, &
            kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
            nspec_ani, &
            c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
            c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
            c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
            ibool,idoubling,rmass,rmass_ocean_load,nglob_oceans, &
            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
            jacobian2D_xmin,jacobian2D_xmax, &
            jacobian2D_ymin,jacobian2D_ymax, &
            jacobian2D_bottom,jacobian2D_top, &
            iMPIcut_xi,iMPIcut_eta,nspec,nglob, &
            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
            TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS, &
            tau_s,tau_e_store,Qmu_store,T_c_source, &
            ATTENUATION,ATTENUATION_3D, &
            size(tau_e_store,2),size(tau_e_store,3),size(tau_e_store,4),size(tau_e_store,5),&
            NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NEX_XI,ichunk,NCHUNKS)

  do ispec=1,nspec
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          weight=wxgll(i)*wygll(j)*wzgll(k)

          jacobianl=jacobianstore(i,j,k,ispec)

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

! deallocate arrays

  deallocate(rmass)
  deallocate(rmass_ocean_load)
  deallocate(rhostore,kappavstore,kappahstore)
  deallocate(muvstore,muhstore)
  deallocate(eta_anisostore)

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

  deallocate(xixstore,xiystore,xizstore)
  deallocate(etaxstore,etaystore,etazstore)
  deallocate(gammaxstore,gammaystore,gammazstore,jacobianstore)
  deallocate(iboun)
  deallocate(iglob)
  deallocate(locval)
  deallocate(ifseg)
  deallocate(xp,yp,zp)
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

! Stacey
  deallocate(nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta)
  deallocate(rho_vp,rho_vs)

! attenuation
  deallocate(Qmu_store)
  deallocate(tau_e_store)

  end subroutine create_regions_mesh

