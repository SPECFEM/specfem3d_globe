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

  use constants,only: CUSTOM_REAL,NDIM,NGLLZ

  implicit none

  ! for Legendre polynomials
  double precision, parameter :: JACALPHA = 0.0d0, JACBETA = 0.0d0

  ! size factor for infinite mesh with reference to hmin
  real(kind=CUSTOM_REAL), parameter :: hfac = 1.4_CUSTOM_REAL

  ! pole for transition and infinite meshing
  real(kind=CUSTOM_REAL), dimension(NDIM) :: xpole

  ! Legendre point position & weights
  double precision, dimension(NGLLZ) :: zgll, wgll

  ! reference (top) layer
  integer :: nspec0 = 0
  logical, dimension(:,:), allocatable :: iboun0, iMPIcut0_xi, iMPIcut0_eta
  integer, dimension(:,:,:,:), allocatable :: ibool0
  double precision, dimension(:,:,:,:), allocatable :: xstore0,ystore0,zstore0

  end module SIEM_meshfem_par


!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_mesh_setup_layers(ipass)

  use constants
  use meshfem_par, only: iregion_code

  use SIEM_meshfem_par

  implicit none

  integer,intent(in) :: ipass

  ! local parameters
  integer :: i
  double precision :: left,right !,r1
  double precision :: ratio(NGLLZ-2),invratio(NGLLZ-2)

  ! check
  if (iregion_code == IREGION_TRINFINITE .and. (.not. ADD_TRINF)) then
    call exit_mpi(myrank,'Error: SIEM meshing called for TRINFINITE for invalid ADD_TRINF')
  endif

  ! only needs to be done for 1. pass
  if (ipass == 1) then
    ! compute ratio for division formula (Legendre polynomial)
    call zwgljd(zgll,wgll,NGLLZ,JACALPHA,JACBETA)

    do i = 1,NGLLZ-2
      ! equidistant points
      ! left=real(i,kreal); right=real(NGLLZ-1-i,kreal)
      ! GLL points
      left  = zgll(i+1) + 1.d0
      right = 1.d0 - zgll(i+1)
      ratio(i) = left/right
      invratio(i) = right/(left+right) ! 1/(ratio+1)
    enddo

    ! set pole for transition and infinite meshing
    xpole(:) = 0.0_CUSTOM_REAL    ! center of the Earth
    if (iregion_code == IREGION_INFINITE) then
      ! center of the source or disturbance
      xpole = (/ -0.6334289, 0.4764568, 0.6045561 /)
    endif
  endif

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

    ! re-allocate reference arrays
    if (allocated(ibool0)) deallocate(ibool0,iboun0,xstore0,ystore0,zstore0,iMPIcut0_xi,iMPIcut0_eta)

    allocate(ibool0(ngllx,nglly,ngllz,nspec0), &
             iboun0(6,nspec0), &
             xstore0(ngllx,nglly,ngllz,nspec0), &
             ystore0(ngllx,nglly,ngllz,nspec0), &
             zstore0(ngllx,nglly,ngllz,nspec0), &
             iMPIcut0_xi(2,nspec0), &
             iMPIcut0_eta(2,nspec0),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool0 arrays'

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
  endif

  end subroutine SIEM_mesh_set_reference_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_mesh_create_elements(iregion_code,ilayer_loop)

! from original: create_regions_meshINF(..)

  use constants
  use shared_parameters, only: RINF
  use meshfem_par, only: xstore,ystore,zstore,NSPEC_REGIONS
  use regions_mesh_par2, only: iboun,iMPIcut_xi,iMPIcut_eta

  use SIEM_meshfem_par

  implicit none

  integer,intent(in) :: iregion_code,ilayer_loop

  ! local parameters
  integer :: ilayer,nlayer
  integer :: i,j,k,iglob,ispec,ispec_n,ier
  integer :: nspec
  integer :: ib,ibool_shift,ip

  ! temporary arrays
  integer, dimension(:), allocatable :: nodelist,inode_order,ibnew
  integer, dimension(:), allocatable :: ispecnew
  integer, dimension(:,:), allocatable :: index_gll
  integer, dimension(:,:,:,:), allocatable :: ibool_inf,iboolold
  logical, dimension(:), allocatable :: isnode

  integer, dimension(:,:,:,:), allocatable :: iboolt
  logical, dimension(:,:), allocatable :: ibount,iMPIcutt_xi,iMPIcutt_eta
  double precision, dimension(:,:,:,:), allocatable :: xstoret,ystoret,zstoret

  real(kind=CUSTOM_REAL),allocatable :: xp(:,:),xs(:,:,:)
  real(kind=CUSTOM_REAL) :: gaminf,invgaminf,hmin,hming
  real(kind=CUSTOM_REAL) :: po(NDIM),px(NDIM),py(NDIM),pz(NDIM)

  integer :: igllx,iglly,inum,nglob_inf,nsnode,nsnode_all

  double precision :: r1
  double precision :: ratio(NGLLZ-2),invratio(NGLLZ-2)

  ! for visualization
  !real, dimension(:,:,:,:), allocatable :: rxstore,rystore,rzstore

  ! initializes
  ilayer = ilayer_loop     ! no need for layer permutations
  nlayer = 0
  select case(iregion_code)
  case (IREGION_TRINFINITE)
    nlayer = NLAYER_TRINF
  case (IREGION_INFINITE)
    nlayer = 1
  case default
    call exit_mpi(myrank,'Invalid region code for SIEM meshing')
  end select

  nspec = NSPEC_REGIONS(iregion_code)

  ! checks nspec and reference nspec0 values
  if (nspec /= nlayer * nspec0) stop 'ERROR: number of infinite elements mismatch!'

  ! we will need to determine the following arrays:
  ! - iMPIcut_xi(:,:) and iboun(:,:) for boundary elements
  ! - xstore/ystore/zstore with the element GLL point locations
  ! - xixstore/xiystore/xizstore,etaxstore/../../,gammaxstore/../.. for Jacobian values (in ipass == 2)
  ! - rhostore,kappavstore,.. for material properties

  ! allocates temporary arrays
  allocate(iboolt(NGLLX,NGLLY,NGLLZ,nspec0), &
           ibount(6,nspec0),                     &
           iMPIcutt_xi(2,nspec0), &
           iMPIcutt_eta(2,nspec0), &
           xstoret(NGLLX,NGLLY,NGLLZ,nspec0), &
           ystoret(NGLLX,NGLLY,NGLLZ,nspec0), &
           zstoret(NGLLX,NGLLY,NGLLZ,nspec0),stat=ier)
  if (ier /= 0) stop 'Error allocating iboolt arrays'

  iboolt(:,:,:,:) = ibool0(:,:,:,:)
  ibount(:,:) = iboun0(:,:)
  iMPIcutt_xi(:,:) = iMPIcut0_xi(:,:)
  iMPIcutt_eta(:,:) = iMPIcut0_eta(:,:)
  xstoret(:,:,:,:) = xstore0(:,:,:,:)
  ystoret(:,:,:,:) = ystore0(:,:,:,:)
  zstoret(:,:,:,:) = zstore0(:,:,:,:)

  ! for visualization
  ! extract surface mesh from which the infinite elements have to be created
  !allocate(rxstore(NGLLX,NGLLY,NGLLZ,nspec), &
  !         rystore(NGLLX,NGLLY,NGLLZ,nspec), &
  !         rzstore(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
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

  ! inifinite element mesh ibool
  allocate(ibool_inf(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating ibool_inf array'
  ibool_inf(:,:,:,:) = -1

  allocate(ispecnew(nspec0),stat=ier)
  if (ier /= 0) stop 'Error allocating ispecnew array'
  ispecnew(:) = 0

  ib = 0
  iboun(:,:) = .false.

  ! ispecnew
  do ispec = 1,nspec0
    ispecnew(ispec) = (ilayer - 1) * nspec0 + ispec
  enddo

  ! set boundary flags
  do ispec = 1,nspec0
    ispec_n = ispecnew(ispec)
    iMPIcut_xi(:,ispec_n) = iMPIcutt_xi(:,ispec)
    iMPIcut_eta(:,ispec_n) = iMPIcutt_eta(:,ispec)
    iboun(1:4,ispec_n) = ibount(1:4,ispec)
  enddo

  ! bottom boundary
  if (ilayer == 1) then
    do ispec = 1,nspec0
      ispec_n = ispecnew(ispec)
      iboun(5,ispec_n) = .true.
    enddo
  endif

  ! top boundary
  if (ilayer == nlayer) then
    do ispec = 1,nspec0
      ispec_n = ispecnew(ispec)
      iboun(6,ispec_n) = .true.
    enddo
  endif

  ! surface node list
  inum = 0
  do ispec = 1,nspec0
    do j = 1,NGLLY
      do i = 1,NGLLX
        inum = inum+1
        nodelist(inum) = iboolt(i,j,NGLLZ,ispec)
        index_gll(1,inum) = i
        index_gll(2,inum) = j
        index_gll(3,inum) = ispec
      enddo
    enddo
  enddo

  if (nsnode_all /= inum) then
    write(*,*)'ERROR: total number of surface nodes mismatch!'
    stop 'Invalid number of surface nodes'
  endif

  call i_uniinv(nodelist,inode_order)

  nsnode = maxval(inode_order)

  allocate(isnode(nsnode),xp(ndim,nsnode),xs(ndim,nsnode,NGLLZ))
  isnode(:) = .false.

  ! assign surface nodes: xs
  do i = 1,nsnode_all
    if (.not. isnode(inode_order(i))) then
      !xs(:,inode_order(i))=g_coord(:,nodelist(i))
      igllx = index_gll(1,i)
      iglly = index_gll(2,i)
      ispec = index_gll(3,i)
      xs(1,inode_order(i),1) = xstoret(igllx,iglly,NGLLZ,ispec)
      xs(2,inode_order(i),1) = ystoret(igllx,iglly,NGLLZ,ispec)
      xs(3,inode_order(i),1) = zstoret(igllx,iglly,NGLLZ,ispec)
      isnode(inode_order(i)) = .true.
    endif
  enddo
  deallocate(isnode)

  ! pole specific to the spherical body which has the center at (0,0,0)
  !xpole=(/ -0.6334289, 0.4764568, 0.6045561 /) !0.0_kreal ! center of the source or disturbance

  if (ilayer == 1) then
    ! find minimum size hmin
    po = (/ xstoret(1,1,1,1),ystoret(1,1,1,1),zstoret(1,1,1,1) /)
    px = (/ xstoret(NGLLX,1,1,1),ystoret(NGLLX,1,1,1),zstoret(NGLLX,1,1,1) /)
    py = (/ xstoret(1,NGLLY,1,1),ystoret(1,NGLLY,1,1),zstoret(1,NGLLY,1,1) /)
    pz = (/ xstoret(1,1,NGLLZ,1),ystoret(1,1,NGLLZ,1),zstoret(1,1,NGLLZ,1) /)

    !hmin=minval((/ distance(po,px,NDIM),distance(po,py,NDIM),distance(po,pz,NDIM) /))
    hmin = distance(po,pz,NDIM) ! only vertical distance

    ! find the common minimum across all processors
    ! local minimum gives the error, particularly in more infinite layers during
    ! MPI interface mismatching
    call min_all_all_cr(hmin,hming)
  endif

  hming = hfac * hming

  ! compute mirror nodes
  do i = 1,nsnode
    r1 = distance(xpole,xs(:,i,1),NDIM)

    RINF = r1 + hming ! infinite surface radius
    if (RINF <= r1) then
      print *,i,xs(:,i,1),r1
      print *,'ERROR: reference infinite radius is smaller than the model!'
      stop 'Invalid infinite radius'
    endif

    gaminf = r1/(RINF-r1)
    invgaminf = (RINF-r1)/r1 !r1/(RINF-r1) use inverse instead for multiplication

    if(gaminf == 0.0_CUSTOM_REAL)then
      print *,'ERROR: zero division!'
      stop 'Invalid gaminf zero division'
    endif

    ! division formula
    xs(:,i,NGLLZ) = ((gaminf+one) * xs(:,i,1) - xpole) * invgaminf
    !mirxs2(:,i)=((gaminf+one)*xs(:,i)-xpole)*gaminf

    do ip = 2,NGLLZ-1 ! first and last points are known
       xs(:,i,ip) = (ratio(ip-1) * xs(:,i,1) + xs(:,i,NGLLZ)) * invratio(ip-1)
    enddo
    !mirxs1(:,i)=0.5_kreal*(mirxs2(:,i)+xs(:,i)) ! midpoint

    !g_numinf(i)=nnode+i
  enddo

  ! allocate global node - and element-arrays
  allocate(iboolold(NGLLX,NGLLY,NGLLZ,nspec0))

  iboolold(:,:,:,:) = -1
  ibool_shift = 0
  do k = 1,NGLLZ
    inum = 0
    xp = xs(:,:,k)
    !if(k.eq.1)xp=xs
    !if(k.eq.2)xp=mirxs1
    !if(k.eq.3)xp=mirxs2
    do ispec = 1,nspec0
      ispec_n = ispecnew(ispec)
      do j = 1,NGLLY
        do i = 1,NGLLX
          inum = inum+1
          ip = inode_order(inum)
          iboolold(i,j,k,ispec) = ip + ibool_shift

          xstore(i,j,k,ispec_n) = xp(1,ip)
          ystore(i,j,k,ispec_n) = xp(2,ip)
          zstore(i,j,k,ispec_n) = xp(3,ip)
        enddo
      enddo
    enddo
    ibool_shift = ibool_shift + nsnode
  enddo
  deallocate(xp,xs)  !,mirxs1,mirxs2)

  nglob_inf = NGLLZ * nsnode

  ! rearrange ibool so as to make consistent with the convention followed by other regions
  allocate(ibnew(nglob_inf),isnode(nglob_inf))

  ibnew(:) = -1
  isnode(:) = .false.

  ! bottom surface of layer 2 has same ibool as the top surface of the layer 1
  if (ilayer > 1) then
    !ibool_inf(:,:,1,ispecnew)=ibool0(:,:,NGLLZ,:)
    do ispec = 1,nspec0
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = iboolold(i,j,1,ispec)
          isnode(iglob) = .true.
          ibnew(iglob) = iboolt(i,j,NGLLZ,ispec)
        enddo
      enddo
    enddo
    !ib=ib+maxval(ibnew) ! next counting should start from this value
  endif

  do ispec = 1,nspec0
    ispec_n = ispecnew(ispec)
    do k = 1,NGLLZ ! 1 was previously set for ilayer>1 but not other
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = iboolold(i,j,k,ispec)
          if (.not. isnode(iglob)) then
            ib = ib+1
            isnode(iglob) = .true.
            ibnew(iglob) = ib
          endif
          ibool_inf(i,j,k,ispec_n) = ibnew(iglob)
        enddo
      enddo
    enddo
  enddo

  !    do ispec = 1,nspec0
  !      ispec_n = ispecnew(ispec)
  !      do k = 1,NGLLZ
  !        do j = 1,NGLLY
  !          do i = 1,NGLLX
  !            iglob = iboolold(i,j,k,ispec)
  !            if(.not.isnode(iglob))then
  !              ib = ib+1
  !              isnode(iglob) = .true.
  !              ibnew(iglob) = ib
  !            endif
  !            ibool_inf(i,j,k,ispec_n) = ibnew(iglob)
  !          enddo
  !        enddo
  !      enddo
  !    enddo

  deallocate(ibnew,iboolold,isnode)

  ! re-set top arrays
  if (ilayer < nlayer) then
    iboolt(:,:,:,:) = ibool_inf(:,:,:,ispecnew)
    ibount(:,:) = iboun(:,ispecnew)
    iMPIcutt_xi(:,:) = iMPIcut_xi(:,ispecnew)
    iMPIcutt_eta(:,:) = iMPIcut_eta(:,ispecnew)
    xstoret(:,:,:,:) = xstore(:,:,:,ispecnew)
    ystoret(:,:,:,:) = ystore(:,:,:,ispecnew)
    zstoret(:,:,:,:) = zstore(:,:,:,ispecnew)
    !hming=hfac*hming
  endif

  deallocate(iboolt,ibount,iMPIcutt_xi,iMPIcutt_eta,xstoret,ystoret,zstoret)
  deallocate(index_gll,nodelist,inode_order)
  deallocate(ispecnew)
  deallocate(xp,xs)

contains

  function distance(x1,x2,n) result(r)
    implicit none
    integer,intent(in) :: n
    real(kind=CUSTOM_REAL),intent(in) :: x1(n),x2(n)
    real(kind=CUSTOM_REAL) :: dx(n),r

    dx = x1-x2
    r = sqrt(sum(dx*dx))
    return
  end function distance

  !=======================================================

  ! Author: Michel Olagnon
  ! orderpack 2.0
  ! source: http://www.fortran-2000.com/rank/

  subroutine i_uniinv (XDONT, IGOEST)
    ! UNIINV = Merge-sort inverse ranking of an array, with removal of
    ! duplicate entries.
    ! this routine is similar to pure merge-sort ranking, but on
    ! the last pass, it sets indices in IGOEST to the rank
    ! of the value in the ordered set with duplicates removed.
    ! for performance reasons, the first 2 passes are taken
    ! out of the standard loop, and use dedicated coding.
    implicit none
    integer,intent(in)  :: XDONT(:)
    integer,intent(out) :: IGOEST(:)

    integer :: XTST, XDONA, XDONB
    integer, dimension (SIZE(IGOEST)) :: JWRKT, IRNGT
    integer :: LMTNA, LMTNC, IRNG, IRNG1, IRNG2, NUNI
    integer :: NVAL, IIND, IWRKD, IWRK, IWRKF, JINDA, IINDA, IINDB

    NVAL = Min (SIZE(XDONT), SIZE(IGOEST))
    select case (NVAL)
    case (:0)
      return
    case (1)
      IGOEST (1) = 1
      return
    case default
      continue
    end select

    ! fill-in the index array, creating ordered couples
    do IIND = 2, NVAL, 2
      if(XDONT(IIND-1) < XDONT(IIND)) then
        IRNGT (IIND-1) = IIND - 1
        IRNGT (IIND) = IIND
      else
        IRNGT (IIND-1) = IIND
        IRNGT (IIND) = IIND - 1
      endif
    enddo
    if(modulo(NVAL,2) /= 0) then
      IRNGT (NVAL) = NVAL
    endif

    ! we will now have ordered subsets A - B - A - B - ...
    ! and merge A and B couples into     C   -   C   - ...
    LMTNA = 2
    LMTNC = 4

    ! first iteration. The length of the ordered subsets goes from 2 to 4
    do
      if (NVAL <= 4) Exit
      ! loop on merges of A and B into C
      do IWRKD = 0, NVAL - 1, 4
        if ((IWRKD+4) > NVAL) then
          if ((IWRKD+2) >= NVAL) Exit
          !   1 2 3
          if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Exit
          !   1 3 2
          if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
            IRNG2 = IRNGT (IWRKD+2)
            IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
            IRNGT (IWRKD+3) = IRNG2
            !   3 1 2
          else
            IRNG1 = IRNGT (IWRKD+1)
            IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
            IRNGT (IWRKD+3) = IRNGT (IWRKD+2)
            IRNGT (IWRKD+2) = IRNG1
          endif
          exit
        endif
        !   1 2 3 4
        if (XDONT(IRNGT(IWRKD+2)) <= XDONT(IRNGT(IWRKD+3))) Cycle
        !   1 3 x x
        if (XDONT(IRNGT(IWRKD+1)) <= XDONT(IRNGT(IWRKD+3))) then
          IRNG2 = IRNGT (IWRKD+2)
          IRNGT (IWRKD+2) = IRNGT (IWRKD+3)
          if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
            !   1 3 2 4
            IRNGT (IWRKD+3) = IRNG2
          else
            !   1 3 4 2
            IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
            IRNGT (IWRKD+4) = IRNG2
          endif
          !   3 x x x
        else
          IRNG1 = IRNGT (IWRKD+1)
          IRNG2 = IRNGT (IWRKD+2)
          IRNGT (IWRKD+1) = IRNGT (IWRKD+3)
          if (XDONT(IRNG1) <= XDONT(IRNGT(IWRKD+4))) then
            IRNGT (IWRKD+2) = IRNG1
            if (XDONT(IRNG2) <= XDONT(IRNGT(IWRKD+4))) then
              !   3 1 2 4
              IRNGT (IWRKD+3) = IRNG2
            else
              !   3 1 4 2
              IRNGT (IWRKD+3) = IRNGT (IWRKD+4)
              IRNGT (IWRKD+4) = IRNG2
            endif
          else
            !   3 4 1 2
            IRNGT (IWRKD+2) = IRNGT (IWRKD+4)
            IRNGT (IWRKD+3) = IRNG1
            IRNGT (IWRKD+4) = IRNG2
          endif
        endif
      enddo

    ! the Cs become As and Bs
      LMTNA = 4
      Exit
    enddo

    ! iteration loop. Each time, the length of the ordered subsets
    ! is doubled.
    do
      if (2*LMTNA >= NVAL) Exit
      IWRKF = 0
      LMTNC = 2 * LMTNC

      ! loop on merges of A and B into C
      do
        IWRK = IWRKF
        IWRKD = IWRKF + 1
        JINDA = IWRKF + LMTNA
        IWRKF = IWRKF + LMTNC
        if (IWRKF >= NVAL) then
          if (JINDA >= NVAL) Exit
          IWRKF = NVAL
        endif
        IINDA = 1
        IINDB = JINDA + 1

        ! one steps in the C subset, that we create in the final rank array
        ! make a copy of the rank array for the iteration
        JWRKT (1:LMTNA) = IRNGT (IWRKD:JINDA)
        XDONA = XDONT (JWRKT(IINDA))
        XDONB = XDONT (IRNGT(IINDB))
        do
          IWRK = IWRK + 1
          ! we still have unprocessed values in both A and B
          if (XDONA > XDONB) then
            IRNGT (IWRK) = IRNGT (IINDB)
            IINDB = IINDB + 1
            if (IINDB > IWRKF) then
              ! only A still with unprocessed values
              IRNGT (IWRK+1:IWRKF) = JWRKT (IINDA:LMTNA)
              Exit
            endif
            XDONB = XDONT (IRNGT(IINDB))
          else
            IRNGT (IWRK) = JWRKT (IINDA)
            IINDA = IINDA + 1
            if (IINDA > LMTNA) Exit! Only B still with unprocessed values
            XDONA = XDONT (JWRKT(IINDA))
          endif

        enddo
      enddo

      ! the Cs become As and Bs
      LMTNA = 2 * LMTNA
    enddo

    ! last merge of A and B into C, with removal of duplicates.
    IINDA = 1
    IINDB = LMTNA + 1
    NUNI = 0

    ! one steps in the C subset, that we create in the final rank array
    JWRKT (1:LMTNA) = IRNGT (1:LMTNA)
    if (IINDB <= NVAL) then
      XTST = i_nearless(Min(XDONT(JWRKT(1)), XDONT(IRNGT(IINDB))))
    else
      XTST = i_nearless(XDONT(JWRKT(1)))
    endif

    do IWRK = 1, NVAL
      ! we still have unprocessed values in both A and B
      if (IINDA <= LMTNA) then
        if (IINDB <= NVAL) then
          if (XDONT(JWRKT(IINDA)) > XDONT(IRNGT(IINDB))) then
            IRNG = IRNGT (IINDB)
            IINDB = IINDB + 1
          else
            IRNG = JWRKT (IINDA)
            IINDA = IINDA + 1
          endif
        else
          ! only A still with unprocessed values
          IRNG = JWRKT (IINDA)
          IINDA = IINDA + 1
        endif
      else
        ! only B still with unprocessed values
        IRNG = IRNGT (IWRK)
      endif
      if (XDONT(IRNG) > XTST) then
        XTST = XDONT (IRNG)
        NUNI = NUNI + 1
      endif
      IGOEST (IRNG) = NUNI
    enddo
  end subroutine i_uniinv

  !=======================================================

  function i_nearless (XVAL) result (I_nl)
    ! nearest value less than given value
    implicit none
    integer,intent(in) :: XVAL
    integer :: I_nl
    I_nl = XVAL - 1
    return
  end function i_nearless

  end subroutine SIEM_mesh_create_elements
