!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

! spectral infinite-element method (SIEM) - meshing routines
!
! SIEM reference:
!   Gharti, H.N., J. Tromp, 2017.
!   A spectral-infinite-element solution of Poisson's equation: an application to self gravity
!   preprint, arXiv: 1706.00855
!   https://arxiv.org/abs/1706.00855
!
!
!-------------------------------------------------------------------------------------------------
!

  module SIEM_meshfem_par

  use constants, only: NDIM,NGLLZ

  implicit none

  ! size factor for infinite mesh with reference to hmin
  double precision, parameter :: hfac = 1.4d0

  ! minimum vertical element size for all processes
  double precision :: hmin_glob

  ! pole for transition and infinite meshing
  double precision, dimension(NDIM) :: xpole

  ! node location ratios (between bottom and top nodes)
  double precision :: ratio(NGLLZ-2),invratio(NGLLZ-2)

  ! reference (top) layer
  integer :: nspec0 = 0
  integer, dimension(:,:,:,:), allocatable :: ibool0
  logical, dimension(:,:), allocatable :: iboun0, iMPIcut0_xi, iMPIcut0_eta
  double precision, dimension(:,:,:,:), allocatable :: xstore0,ystore0,zstore0

  ! infinite element mesh ibool
  integer, dimension(:,:,:,:), allocatable :: ibool_inf
  integer :: ib_counter

  ! temporary layer arrays
  integer, dimension(:,:,:,:), allocatable :: iboolt
  logical, dimension(:,:), allocatable :: ibount,iMPIcutt_xi,iMPIcutt_eta
  double precision, dimension(:,:,:,:), allocatable :: xstoret,ystoret,zstoret

  end module SIEM_meshfem_par


!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_mesh_setup_layers(ipass)

  use constants
  use meshfem_par, only: iregion_code,NSPEC_REGIONS
  use regions_mesh_par2, only: iboun

  use SIEM_meshfem_par

  implicit none

  integer,intent(in) :: ipass

  ! local parameters
  integer :: ip,ier,nspec_all
  double precision :: left,right

  ! for Legendre polynomials
  double precision, parameter :: JACALPHA = 0.0d0, JACBETA = 0.0d0
  ! Legendre point position & weights
  double precision, dimension(NGLLZ) :: zgll, wgll

  ! check if anything to do
  if (iregion_code == IREGION_TRINFINITE .and. (.not. ADD_TRINF)) return

  ! only needs to be done for 1. pass
  if (ipass == 1) then
    ! note: xigll/yigll/zigll are already setup for Gauss-Lobatto-Legendre points
    !
    ! compute ratio for division formula (for Legendre polynomial)
    call zwgljd(zgll,wgll,NGLLZ,JACALPHA,JACBETA)

    do ip = 2,NGLLZ-1
      ! equidistant points
      !left  = real((ip-1),kind=8)
      !right = real(NGLLZ-1-(ip-1),kind=8)

      ! GLL points
      left  = 1.d0 + zgll(ip)
      right = 1.d0 - zgll(ip)

      ratio(ip-1) = left/right
      invratio(ip-1) = right/(left+right) ! 1/(ratio+1)
    enddo

    ! set pole for transition and infinite meshing
    xpole(:) = 0.d0    ! center of the Earth

    if (iregion_code == IREGION_INFINITE) then
      ! pole for infinite element region (see constants.h)
      xpole(:) = POLE_INF(:)
    endif

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  SIEM mesh setup: xpole = ',xpole(:)
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! infinite element mesh ibool
  if (allocated(ibool_inf)) deallocate(ibool_inf)

  ! initializes ibool_inf for this region
  nspec_all = NSPEC_REGIONS(iregion_code)

  allocate(ibool_inf(NGLLX,NGLLY,NGLLZ,nspec_all),stat=ier)
  if (ier /= 0) stop 'Error allocating ibool_inf array'
  ibool_inf(:,:,:,:) = -1

  ! boundary node counter
  ib_counter = 0

  ! initializes boundary element flags
  iboun(:,:) = .false.

  end subroutine SIEM_mesh_setup_layers

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_mesh_set_reference_arrays(iregion_code,ipass)

  use constants
  use meshfem_par, only: ibool,xstore,ystore,zstore,NSPEC_REGIONS
  use regions_mesh_par2, only: iboun,iMPIcut_xi,iMPIcut_eta

  use SIEM_meshfem_par

  implicit none

  integer,intent(in) :: iregion_code,ipass

  ! local parameters
  integer :: ispec,ielmt,ier

  ! checks to do only once after all passes are done
  if (ipass == 1) return

  ! only store after crust_mantle and transition-to-infinite regions are done
  if (iregion_code == IREGION_CRUST_MANTLE .or. &
      (ADD_TRINF .and. iregion_code == IREGION_TRINFINITE)) then
    ! sets reference values for transition-to-infinite and infinite layer
    ! only surface elements
    nspec0 = count(iboun(6,:))

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'setting surface reference for SIEM:'
      write(IMAIN,*) '  number of surface elements = ',nspec0
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! re-allocate reference arrays
    if (allocated(ibool0)) deallocate(ibool0,iboun0,xstore0,ystore0,zstore0,iMPIcut0_xi,iMPIcut0_eta)

    allocate(ibool0(NGLLX,NGLLY,NGLLZ,nspec0), &
             iboun0(6,nspec0), &
             xstore0(NGLLX,NGLLY,NGLLZ,nspec0), &
             ystore0(NGLLX,NGLLY,NGLLZ,nspec0), &
             zstore0(NGLLX,NGLLY,NGLLZ,nspec0), &
             iMPIcut0_xi(2,nspec0), &
             iMPIcut0_eta(2,nspec0),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool0 arrays'
    ibool0(:,:,:,:) = 0
    iboun0(:,:) = .false.
    xstore0(:,:,:,:) = 0.d0
    ystore0(:,:,:,:) = 0.d0
    zstore0(:,:,:,:) = 0.d0
    iMPIcut0_xi(:,:) = .false.
    iMPIcut0_eta(:,:) = .false.

    ! stores top layer elements
    ielmt = 0
    do ispec = 1,NSPEC_REGIONS(iregion_code)
      ! only if at the surface (top plane)
      if (iboun(6,ispec)) then
        ielmt = ielmt+1
        ibool0(:,:,:,ielmt) = ibool(:,:,:,ispec)
        iboun0(:,ielmt) = iboun(:,ispec)

        xstore0(:,:,:,ielmt) = xstore(:,:,:,ispec)
        ystore0(:,:,:,ielmt) = ystore(:,:,:,ispec)
        zstore0(:,:,:,ielmt) = zstore(:,:,:,ispec)

        iMPIcut0_xi(:,ielmt) = iMPIcut_xi(:,ispec)
        iMPIcut0_eta(:,ielmt) = iMPIcut_eta(:,ispec)
      endif
    enddo

    ! double-check
    if (ielmt /= nspec0) stop 'Error invalid number of surface elements set'
  endif

  end subroutine SIEM_mesh_set_reference_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_mesh_create_elements(ilayer,ispec_count,ipass,iregion_code)

! from original: create_regions_meshINF(..)

  use constants
  use shared_parameters, only: RINF

  use siem_math_library, only: i_uniinv

  use meshfem_par, only: xstore,ystore,zstore, &
    NSPEC_REGIONS,NGLOB_REGIONS

  use regions_mesh_par2, only: iboun,iMPIcut_xi,iMPIcut_eta

  ! for jacobian
  use regions_mesh_par2, only: &
    xixstore,xiystore,xizstore, &
    etaxstore,etaystore,etazstore, &
    gammaxstore,gammaystore,gammazstore

  use SIEM_meshfem_par

  !debugging
  use shared_parameters, only: NPROCTOT

  implicit none

  integer,intent(in) :: ilayer,ipass,iregion_code
  integer,intent(inout) :: ispec_count

  ! local parameters
  integer :: nlayer
  integer :: i,j,k,iglob,iglob_new
  integer :: ispec_layer,ispec
  integer :: nspec_all,nglob_all
  integer :: ibool_shift,ip,ier

  ! temporary arrays
  integer, dimension(:), allocatable :: nodelist,inode_order,ibnew
  integer, dimension(:), allocatable :: ispecnew
  integer, dimension(:,:), allocatable :: index_gll
  integer, dimension(:,:,:,:), allocatable :: iboolold
  logical, dimension(:), allocatable :: isnode

  double precision,allocatable :: xs(:,:,:)
  double precision :: gaminf,invgaminf
  double precision :: hmin
  double precision :: po(NDIM),px(NDIM),py(NDIM),pz(NDIM)

  integer :: igllx,iglly,inum,nglob_inf,nsnode,nsnode_all

  double precision :: r1

  ! for visualization
  !real, dimension(:,:,:,:), allocatable :: rxstore,rystore,rzstore

  ! initializes
  nlayer = 0
  select case(iregion_code)
  case (IREGION_TRINFINITE)
    nlayer = NLAYER_TRINF
  case (IREGION_INFINITE)
    nlayer = 1
  case default
    call exit_mpi(myrank,'Invalid region code for SIEM meshing')
  end select

  nspec_all = NSPEC_REGIONS(iregion_code)
  nglob_all = NGLOB_REGIONS(iregion_code)

  ! checks nspec_all and reference nspec0 values
  if (nspec_all /= nlayer * nspec0) stop 'ERROR: number of infinite elements mismatch!'

  ! we will need to determine the following arrays:
  ! - iMPIcut_xi(:,:) and iboun(:,:) for boundary elements
  ! - xstore/ystore/zstore with the element GLL point locations
  ! - xixstore/xiystore/xizstore,etaxstore/../../,gammaxstore/../.. for Jacobian values (in ipass == 2)
  ! - rhostore,kappavstore,.. for material properties

  ! allocates temporary arrays
  if (ilayer == 1) then
    allocate(iboolt(NGLLX,NGLLY,NGLLZ,nspec0), &
             ibount(6,nspec0), &
             iMPIcutt_xi(2,nspec0), &
             iMPIcutt_eta(2,nspec0), &
             xstoret(NGLLX,NGLLY,NGLLZ,nspec0), &
             ystoret(NGLLX,NGLLY,NGLLZ,nspec0), &
             zstoret(NGLLX,NGLLY,NGLLZ,nspec0),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolt arrays'

    ! sets surface reference elements
    iboolt(:,:,:,:) = ibool0(:,:,:,:)
    ibount(:,:) = iboun0(:,:)
    iMPIcutt_xi(:,:) = iMPIcut0_xi(:,:)
    iMPIcutt_eta(:,:) = iMPIcut0_eta(:,:)
    xstoret(:,:,:,:) = xstore0(:,:,:,:)
    ystoret(:,:,:,:) = ystore0(:,:,:,:)
    zstoret(:,:,:,:) = zstore0(:,:,:,:)
  endif

  ! for visualization
  ! extract surface mesh from which the infinite elements have to be created
  !allocate(rxstore(NGLLX,NGLLY,NGLLZ,nspec_all), &
  !         rystore(NGLLX,NGLLY,NGLLZ,nspec_all), &
  !         rzstore(NGLLX,NGLLY,NGLLZ,nspec_all),stat=ier)
  !if (ier /= 0) stop 'Error allocating rxstore arrays'

  ! surface nodes
  nsnode_all = NGLLX * NGLLY * nspec0
  allocate(index_gll(3,nsnode_all), &
           nodelist(nsnode_all), &
           inode_order(nsnode_all),stat=ier)
  if (ier /= 0) stop 'Error allocating index_gll arrays'
  index_gll(:,:) = 0
  nodelist(:) = 0
  inode_order(:) = 0

  allocate(ispecnew(nspec0),stat=ier)
  if (ier /= 0) stop 'Error allocating ispecnew array'
  ispecnew(:) = 0

  ! sets ispecnew
  ! (from a layer element numbering [1,nspec0] to a global element numbering [1,nspec_all])
  do ispec_layer = 1,nspec0
    ispecnew(ispec_layer) = (ilayer - 1) * nspec0 + ispec_layer
  enddo

  ! set boundary flags
  do ispec_layer = 1,nspec0
    ispec = ispecnew(ispec_layer)
    iMPIcut_xi(:,ispec) = iMPIcutt_xi(:,ispec_layer)
    iMPIcut_eta(:,ispec) = iMPIcutt_eta(:,ispec_layer)
    iboun(1:4,ispec) = ibount(1:4,ispec_layer)
  enddo

  ! bottom boundary
  if (ilayer == 1) then
    do ispec_layer = 1,nspec0
      ispec = ispecnew(ispec_layer)
      iboun(5,ispec) = .true.
    enddo
  endif

  ! top boundary
  if (ilayer == nlayer) then
    do ispec_layer = 1,nspec0
      ispec = ispecnew(ispec_layer)
      iboun(6,ispec) = .true.
    enddo
  endif

  ! surface node list
  inum = 0
  do ispec_layer = 1,nspec0
    do j = 1,NGLLY
      do i = 1,NGLLX
        inum = inum+1
        nodelist(inum) = iboolt(i,j,NGLLZ,ispec_layer)
        index_gll(1,inum) = i
        index_gll(2,inum) = j
        index_gll(3,inum) = ispec_layer
      enddo
    enddo
  enddo

  if (nsnode_all /= inum) then
    print *,'ERROR: rank',myrank,' total number of surface nodes mismatch!'
    stop 'Invalid number of surface nodes'
  endif

  ! checks that iglob values in nodelist are valid
  if (minval(nodelist) < 1) then
    print *,'ERROR: rank',myrank,' invalid nodelist entries'
    stop 'Invalid nodelist entries'
  endif

  ! sorts nodelist
  call i_uniinv(nodelist,inode_order)

  ! checks that ordering is within nodelist bounds [1,nsnode_all]
  if (minval(inode_order) < 1 .or. maxval(inode_order) > nsnode_all) then
    print *,'ERROR: rank',myrank,' invalid inode_order entries: min/max = ',minval(inode_order),maxval(inode_order)
    print *,'  node all: ',nsnode_all
    stop 'Invalid inode_order entries'
  endif

  ! in case all iglob values would be unique in nodelist(..), then the maxval(inode_order) would be equal to
  ! nsnode_all==NGLLX*NGLLY*nspec0. however, since elements share global nodes at corners, edges and surfaces,
  ! the maximum order number is likely smaller.
  nsnode = maxval(inode_order)

  allocate(isnode(nsnode),xs(ndim,nsnode,NGLLZ),stat=ier)
  if (ier /= 0) stop 'Error allocating isnode arrays'
  isnode(:) = .false.
  xs(:,:,:) = 0.d0

  ! assign surface nodes: xs
  do i = 1,nsnode_all
    if (.not. isnode(inode_order(i))) then
      !xs(:,inode_order(i))=g_coord(:,nodelist(i))
      igllx = index_gll(1,i)
      iglly = index_gll(2,i)
      ispec_layer = index_gll(3,i)
      xs(1,inode_order(i),1) = xstoret(igllx,iglly,NGLLZ,ispec_layer)    ! from NGLLZ - > top nodes
      xs(2,inode_order(i),1) = ystoret(igllx,iglly,NGLLZ,ispec_layer)
      xs(3,inode_order(i),1) = zstoret(igllx,iglly,NGLLZ,ispec_layer)
      isnode(inode_order(i)) = .true.
    endif
  enddo
  deallocate(isnode)

  ! determine minimum vertical element size for surface element
  if (ilayer == 1) then
    ! find minimum size hmin
    po = (/ xstoret(1,1,1,1),ystoret(1,1,1,1),zstoret(1,1,1,1) /)
    px = (/ xstoret(NGLLX,1,1,1),ystoret(NGLLX,1,1,1),zstoret(NGLLX,1,1,1) /)
    py = (/ xstoret(1,NGLLY,1,1),ystoret(1,NGLLY,1,1),zstoret(1,NGLLY,1,1) /)
    pz = (/ xstoret(1,1,NGLLZ,1),ystoret(1,1,NGLLZ,1),zstoret(1,1,NGLLZ,1) /)

    !hmin=minval((/ distance(po,px,NDIM),distance(po,py,NDIM),distance(po,pz,NDIM) /))
    hmin = distance_dble(po,pz,NDIM) ! only vertical distance

    ! find the common minimum across all processors
    ! local minimum gives the error, particularly in more infinite layers during
    ! MPI interface mismatching
    call min_all_all_dp(hmin,hmin_glob)
  endif

  ! apply size factor to vertical element size
  hmin_glob = hfac * hmin_glob

  ! checks
  if (abs(hmin_glob) < TINYVAL) then
    print *,'ERROR: rank',myrank,' hmin_glob is tiny: ',hmin_glob
    stop 'Invalid hmin_glob'
  endif

  ! compute mirror nodes
  do i = 1,nsnode
    r1 = distance_dble(xpole,xs(:,i,1),NDIM)

    RINF = r1 + hmin_glob ! infinite surface radius
    if (RINF <= r1) then
      print *,'ERROR: rank',myrank,' reference infinite radius smaller than the model!'
      print *,'  node          : ',i,xs(:,i,1)
      print *,'  r1/RINF       : ',r1,RINF
      print *,'  hmin_glob/hfac: ',hmin_glob,hfac
      stop 'Invalid infinite radius'
    endif

    ! checks
    if (abs(r1) < TINYVAL) then
      print *,'ERROR: rank',myrank,' distance r1 is tiny'
      stop 'Invalid r1 distance'
    endif

    gaminf = r1/(RINF-r1)
    invgaminf = (RINF-r1)/r1    !r1/(RINF-r1) use inverse instead for multiplication

    ! double check
    if (gaminf == 0.d0) then
      print *,'ERROR: rank',myrank,' division by zero!'
      stop 'Invalid gaminf zero division'
    endif

    ! division formula
    xs(:,i,NGLLZ) = ((gaminf + 1.d0) * xs(:,i,1) - xpole(:)) * invgaminf

    ! first and last points are known
    do ip = 2,NGLLZ-1
      xs(:,i,ip) = (ratio(ip-1) * xs(:,i,1) + xs(:,i,NGLLZ)) * invratio(ip-1)
    enddo
  enddo

  ! allocate global node - and element-arrays
  allocate(iboolold(NGLLX,NGLLY,NGLLZ,nspec0),stat=ier)
  if (ier /= 0) stop 'Error allocating iboolold array'
  iboolold(:,:,:,:) = -1

  ibool_shift = 0
  do k = 1,NGLLZ
    inum = 0
    do ispec_layer = 1,nspec0
      ispec = ispecnew(ispec_layer)
      do j = 1,NGLLY
        do i = 1,NGLLX
          inum = inum+1
          ip = inode_order(inum)
          iboolold(i,j,k,ispec_layer) = ip + ibool_shift

          xstore(i,j,k,ispec) = xs(1,ip,k)
          ystore(i,j,k,ispec) = xs(2,ip,k)
          zstore(i,j,k,ispec) = xs(3,ip,k)
        enddo
      enddo
    enddo
    ibool_shift = ibool_shift + nsnode
  enddo
  deallocate(xs)

  nglob_inf = NGLLZ * nsnode

  ! rearrange ibool so as to make consistent with the convention followed by other regions
  allocate(ibnew(nglob_inf),isnode(nglob_inf),stat=ier)
  if (ier /= 0) stop 'Error allocating ibnew arrays'

  ibnew(:) = -1
  isnode(:) = .false.

  ! bottom surface of layer 2 has same ibool as the top surface of the layer 1
  if (ilayer > 1) then
    do ispec_layer = 1,nspec0
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = iboolold(i,j,1,ispec_layer)

          ! checks iglob
          if (iglob < 1 .or. iglob > nglob_inf) then
            print *,'Error: rank',myrank,' invalid iglob ',iglob,'; should be between 1 to ',nglob_inf
            print *,'  ispec_layer: ',ispec_layer
            call exit_mpi(myrank,'Error invalid iglob for ibnew')
          endif

          ! stores top surface node index from previous layer
          isnode(iglob) = .true.
          ibnew(iglob) = iboolt(i,j,NGLLZ,ispec_layer)
        enddo
      enddo
    enddo
  endif

  ! sets new ibool indexing
  do ispec_layer = 1,nspec0
    ispec = ispecnew(ispec_layer)

    ! checks ispec
    if (ispec < 1 .or. ispec > nspec_all) then
      print *,'Error: rank',myrank,' invalid ispec ',ispec,'; should be between 1 to ',nspec_all
      print *,'  ispec_layer: ',ispec_layer
      call exit_mpi(myrank,'Error invalid ispec for ibool_inf')
    endif

    ! updates ibool mapping
    do k = 1,NGLLZ ! 1 was previously set for ilayer > 1 but not other
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = iboolold(i,j,k,ispec_layer)

          ! checks iglob
          if (iglob < 1 .or. iglob > nglob_inf) then
            print *,'Error: rank',myrank,' invalid iglob ',iglob,'; should be between 1 to ',nglob_inf
            print *,'  ispec_layer: ',ispec_layer
            call exit_mpi(myrank,'Error invalid iglob for ibnew/isnode')
          endif

          ! adds new index if not set yet
          if (.not. isnode(iglob)) then
            ib_counter = ib_counter + 1
            ibnew(iglob) = ib_counter
            isnode(iglob) = .true.
          endif

          ! new ibool index
          iglob_new = ibnew(iglob)

          ! checks iglob_new
          if (iglob_new < 1 .or. iglob_new > nglob_all) then
            print *,'Error: rank',myrank,' invalid iglob_new ',iglob_new,'; should be between 1 to ',nglob_all
            print *,'  ispec_layer/ispec/ib_counter: ',ispec_layer,ispec,ib_counter
            print *,'  iglob/nglob_inf             : ',iglob,nglob_inf
            call exit_mpi(myrank,'Error invalid iglob_new for ibool_inf')
          endif

          ! sets new index
          ibool_inf(i,j,k,ispec) = iglob_new
        enddo
      enddo
    enddo
  enddo

  deallocate(ibnew,iboolold,isnode)

  ! updates Jacobian
  ! (only needed for second meshing phase)
  if (ipass == 2) then
    do ispec_layer = 1,nspec0
      ispec = ispecnew(ispec_layer)
      ! sets an artifical jacobian matrix (J==unit matrix) for SIEM mesh
      ! (since we don't store the jacobian matrix itself, we set the corresponding mapping function derivatives xix/..)
      xixstore(:,:,:,ispec) = 1.0_CUSTOM_REAL
      xiystore(:,:,:,ispec) = 0.0_CUSTOM_REAL
      xizstore(:,:,:,ispec) = 0.0_CUSTOM_REAL
      etaxstore(:,:,:,ispec) = 0.0_CUSTOM_REAL
      etaystore(:,:,:,ispec) = 1.0_CUSTOM_REAL
      etazstore(:,:,:,ispec) = 0.0_CUSTOM_REAL
      gammaxstore(:,:,:,ispec) = 0.0_CUSTOM_REAL
      gammaystore(:,:,:,ispec) = 0.0_CUSTOM_REAL
      gammazstore(:,:,:,ispec) = 1.0_CUSTOM_REAL
    enddo
  endif

  ! re-set top arrays
  if (ilayer < nlayer) then
    iboolt(:,:,:,:) = ibool_inf(:,:,:,ispecnew(:))
    ibount(:,:) = iboun(:,ispecnew(:))
    iMPIcutt_xi(:,:) = iMPIcut_xi(:,ispecnew(:))
    iMPIcutt_eta(:,:) = iMPIcut_eta(:,ispecnew(:))
    xstoret(:,:,:,:) = xstore(:,:,:,ispecnew(:))
    ystoret(:,:,:,:) = ystore(:,:,:,ispecnew(:))
    zstoret(:,:,:,:) = zstore(:,:,:,ispecnew(:))
  endif

  ! debugging
  if (.false.) then
    if (ilayer == 1) then
      ! print out node locations (one process at a time)
      do i = 0,NPROCTOT-1
        if (myrank == i) then
          print *,'debug: rank',i,' ispecnew',ispecnew(1:10)
          do k = 1,NGLLZ
            print *,'debug: xstore ',k,xstore(:,:,k,ispecnew(1))
            call flush_stdout()
          enddo
          print *
          call flush_stdout()
        endif
        call synchronize_all()
      enddo
    endif
  endif

  ! number of elements created in this layer
  ispec_count = ispec_count + nspec0

  ! free temporary arrays
  if (ilayer == nlayer) then
    deallocate(iboolt,ibount,iMPIcutt_xi,iMPIcutt_eta,xstoret,ystoret,zstoret)
  endif
  deallocate(index_gll,nodelist,inode_order)
  deallocate(ispecnew)

contains

  function distance_dble(x1,x2,n) result (r)
  ! distance function in double precision
    implicit none
    integer,intent(in) :: n
    double precision,intent(in) :: x1(n),x2(n)
    double precision :: dx(n),r

    dx(:) = x1(:) - x2(:)
    r = sqrt(sum(dx*dx))
    return
  end function distance_dble

  end subroutine SIEM_mesh_create_elements
