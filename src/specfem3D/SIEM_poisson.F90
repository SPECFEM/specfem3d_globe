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


module siem_poisson

  use constants, only: CUSTOM_REAL

contains

!TODO: replace ibool with inode_elmt

  subroutine poisson_stiffness(iregion,nelmt,nnode,ibool,xstore,ystore,zstore, &
                               storekmat,dprecon)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_innercore, only: idoubling_inner_core
  !use specfem_par_crustmantle, only: rmassz_crust_mantle !TODO: remove this

  use siem_gll_library
  use siem_math_library, only: determinant,invert

  implicit none
  integer,intent(in) :: iregion,nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(out) :: storekmat(NGLLCUBE,NGLLCUBE,nelmt),dprecon(nnode)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,k,i_elmt
  integer :: dnx,dny,dnz,egdof(NGLLCUBE),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ)
  real(kind=kdble) :: dshape_hex8(ndim,ngnod,NGLLCUBE),gll_weights(NGLLCUBE), &
                      gll_points(ndim,NGLLCUBE)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),deriv(ndim,NGLLCUBE),jac(ndim,ndim), &
                            kmat(NGLLCUBE,NGLLCUBE)

  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gll
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gll

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE), &
           dlagrange_gll(ndim,NGLLCUBE,NGLLCUBE))
  lagrange_gll(:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  call zwgljd(xigll,wxgll,NGLLX,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  dnx = NGLLX-1; dny = NGLLY-1; dnz = NGLLZ-1
  storekmat = zero; dprecon = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong!!!!
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);               ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);       ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! top corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    kmat = zero
    do i = 1,NGLLCUBE
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL) ! use der for gll
      kmat = kmat + real(matmul(transpose(deriv),deriv)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    storekmat(:,:,i_elmt) = kmat
    do k = 1,NGLLCUBE
      dprecon(egdof(k)) = dprecon(egdof(k))+kmat(k,k)
    enddo
  enddo
  end subroutine poisson_stiffness

!
!===========================================
!

  subroutine poisson_stiffnessINF(nelmt,nnode,ibool,xstore,ystore,zstore, &
                                  storekmat,dprecon)

  use constants_solver, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use siem_infinite_element
  use siem_math_library, only: determinant,invert

  implicit none
  integer,intent(in) :: nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(out) :: storekmat(NGLLCUBE,NGLLCUBE,nelmt),dprecon(nnode)

  integer,parameter :: iface = 6,ndim = 3,ngnod = 8,nginf = 8
  ! GLL-Radau quadrature
  integer,parameter :: nipinf = NGLLCUBE,nipx = NGLLX !NGLLX = NGLLY = NGLLZ

  !! Gauss quadrature
  !integer,parameter :: nipinf=8,nipx=8

  integer :: i,k,i_elmt
  integer :: egdof(NGLLCUBE),ignod(ngnod) !dnx,dny,dnz,
  real(kind=CUSTOM_REAL),parameter :: one=1.0_CUSTOM_REAL,zero=0.0_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: detjac
  real(kind=CUSTOM_REAL) :: coordinf(nginf,ndim),deriv(ndim,NGLLCUBE)
  real(kind=CUSTOM_REAL) :: jac(ndim,ndim),kmat(NGLLCUBE,NGLLCUBE),x0(ndim)

  real(kind=kdble),parameter :: done=1.0_kdble
  real(kind=kdble) :: gaminf,GLw(nipinf)
  real(kind=kdble) :: shape_infinite(nipinf,nginf),dshape_infinite(ndim,nipinf,nginf)
  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gl
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gl

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gl(nipinf,NGLLCUBE), &
           dlagrange_gl(ndim,nipinf,NGLLCUBE))
  lagrange_gl(:,:) = 0.0_kdble
  dlagrange_gl(:,:,:) = 0.0_kdble

  ! ainf is irrevelant for the time being
  ! nd,ainf,gaminf can be removed from argument list
  gaminf = 1.0002_kdble !1.99_kdble

  x0 = (/ -0.6334289, 0.4764568, 0.6045561 /)!zero ! center of the Earth

  ! GLL-Radau quadrature
  call shape_function_infiniteGLHEX8ZW_GLLR(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,nipinf, &
                                            iface,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
  !! Gauss quadrature
  !call shape_function_infiniteGLHEX8ZW_GQ(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,nipx,nipinf, &
  !                                        iface,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-2
  storekmat = zero; dprecon = zero
  do i_elmt = 1,nelmt
    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ-1:dnz,i_elmt),(/ngnod/))
    ! indicial order NOT EXODUS order!!!
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);         ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(1,NGLLY,1,i_elmt);     ignod(4) = ibool(NGLLX,NGLLY,1,i_elmt)
    ! top corner nodes - second last
    ignod(5) = ibool(1,1,NGLLZ-1,i_elmt);     ignod(6) = ibool(NGLLX,1,NGLLZ-1,i_elmt)
    ignod(7) = ibool(1,NGLLY,NGLLZ-1,i_elmt); ignod(8) = ibool(NGLLX,NGLLY,NGLLZ-1,i_elmt);

    coordinf(:,1) = xstore(ignod)
    coordinf(:,2) = ystore(ignod)
    coordinf(:,3) = zstore(ignod)

    ! Zienkiewicz infinite coordinates
    coordinf(5,:) = x0 + real(gaminf*(coordinf(1,:)-x0),kind=CUSTOM_REAL)
    coordinf(6,:) = x0 + real(gaminf*(coordinf(2,:)-x0),kind=CUSTOM_REAL)
    coordinf(7,:) = x0 + real(gaminf*(coordinf(3,:)-x0),kind=CUSTOM_REAL)
    coordinf(8,:) = x0 + real(gaminf*(coordinf(4,:)-x0),kind=CUSTOM_REAL)

    ! Point X0 (Pole)
    coordinf(1,:) = x0; coordinf(2,:) = x0
    coordinf(3,:) = x0; coordinf(4,:) = x0

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    kmat = zero
    do i = 1,nipinf
      jac = real(matmul(dshape_infinite(:,i,:),coordinf),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gl(:,i,:)),kind=CUSTOM_REAL)
      kmat = kmat + real(matmul(transpose(deriv),deriv)*detjac*GLw(i),kind=CUSTOM_REAL)
    enddo
    storekmat(:,:,i_elmt) = kmat
    do k = 1,NGLLCUBE
      dprecon(egdof(k)) = dprecon(egdof(k))+kmat(k,k)
    enddo
  enddo

  end subroutine poisson_stiffnessINF

!
!===========================================
!

  subroutine poisson_stiffness3(iregion,nelmt,nnode,ibool,xstore,ystore,zstore, &
                                nnode1,ibool1,storekmat,dprecon)

  use constants_solver, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NGLLX,NGLLY,NGLLZ,NGLLCUBE, &
                              NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_innercore, only: idoubling_inner_core
  !use specfem_par_crustmantle, only: rmassz_crust_mantle !TODO: remove this

  use siem_gll_library
  use siem_math_library, only: determinant,invert

  implicit none
  integer,intent(in) :: iregion,nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(out) :: storekmat(NGLLCUBE_INF,NGLLCUBE_INF,nelmt),dprecon(nnode1)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,k,i_elmt
  integer :: dnx,dny,dnz,egdof(NGLLCUBE_INF),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll1(NGLLX_INF),wxgll1(NGLLX_INF),etagll1(NGLLY_INF), &
  wygll1(NGLLY_INF),zetagll1(NGLLZ_INF),wzgll1(NGLLZ_INF)
  real(kind=kdble) :: dshape_hex8(ndim,ngnod,NGLLCUBE_INF),gll_weights(NGLLCUBE_INF), &
                      gll_points(ndim,NGLLCUBE_INF), &
                      lagrange_gll(NGLLCUBE_INF,NGLLCUBE_INF),dlagrange_gll(ndim,NGLLCUBE_INF,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),deriv(ndim,NGLLCUBE_INF),jac(ndim,ndim), &
                      kmat(NGLLCUBE_INF,NGLLCUBE_INF)

  call zwgljd(xigll1,wxgll1,NGLLX_INF,jalpha,jbeta)
  call zwgljd(etagll1,wygll1,NGLLY_INF,jalpha,jbeta)
  call zwgljd(zetagll1,wzgll1,NGLLZ_INF,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,xigll1,etagll1, &
                            zetagll1,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  dnx = NGLLX-1; dny = NGLLY-1; dnz = NGLLZ-1
  storekmat = zero; dprecon = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong!!!!
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! top corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = ibool1(:,i_elmt)!reshape(ibool1(:,:,:,i_elmt),(/NGLLCUBE_INF/))
    kmat = zero
    do i = 1,NGLLCUBE_INF
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL) ! use der for gll
      kmat = kmat + real(matmul(transpose(deriv),deriv)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    storekmat(:,:,i_elmt) = kmat
    do k = 1,NGLLCUBE_INF
      dprecon(egdof(k))=dprecon(egdof(k))+kmat(k,k)
    enddo
  enddo

  end subroutine poisson_stiffness3

!
!===========================================
!

  subroutine poisson_stiffnessINF3(nelmt,nnode,ibool,xstore,ystore,zstore,nnode1, &
                                   ibool1,storekmat,dprecon)

  use constants_solver, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use siem_infinite_element
  use siem_math_library, only: determinant,invert

  implicit none
  integer,intent(in) :: nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(out) :: storekmat(NGLLCUBE_INF,NGLLCUBE_INF,nelmt),dprecon(nnode1)

  integer,parameter :: iface = 6,ndim = 3,ngnod = 8,nginf = 8
  ! GLL-Radau quadrature
  integer,parameter :: nipinf = NGLLCUBE_INF,nipx = NGLLX_INF !NGLLX_INF = NGLLY_INF = NGLLZ_INF

  !! Gauss quadrature
  !integer,parameter :: nipinf=8,nipx=8

  integer :: i,k,i_elmt
  integer :: egdof(NGLLCUBE_INF),ignod(ngnod) !dnx,dny,dnz,
  real(kind=CUSTOM_REAL),parameter :: one=1.0_CUSTOM_REAL,zero=0.0_CUSTOM_REAL
  real(kind=CUSTOM_REAL) :: detjac
  real(kind=CUSTOM_REAL) :: coordinf(nginf,ndim),deriv(ndim,NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: jac(ndim,ndim),kmat(NGLLCUBE_INF,NGLLCUBE_INF),x0(ndim)

  real(kind=kdble),parameter :: done=1.0_kdble
  real(kind=kdble) :: gaminf,GLw(nipinf)
  real(kind=kdble) :: shape_infinite(nipinf,nginf),dshape_infinite(ndim,nipinf,nginf)
  real(kind=kdble) :: lagrange_gl(nipinf,NGLLCUBE_INF),dlagrange_gl(ndim,nipinf,NGLLCUBE_INF)

  ! ainf is irrevelant for the time being
  ! nd,ainf,gaminf can be removed from argument list
  gaminf = 1.0002_kdble !1.99_kdble

  x0 = (/ -0.6334289, 0.4764568, 0.6045561 /)!zero ! center of the Earth

  ! GLL-Radau quadrature
  call shape_function_infiniteGLHEX8ZW_GLLR(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,nipinf, &
                                            iface,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
  !! Gauss quadrature
  !call shape_function_infiniteGLHEX8ZW_GQ(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,nipx,nipinf, &
  !                                        iface,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-2
  storekmat = zero; dprecon = zero
  do i_elmt = 1,nelmt
    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ-1:dnz,i_elmt),(/ngnod/))
    ! indicial order NOT EXODUS order!!!
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);         ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(1,NGLLY,1,i_elmt);     ignod(4) = ibool(NGLLX,NGLLY,1,i_elmt)
    ! top corner nodes - second last
    ignod(5) = ibool(1,1,NGLLZ-1,i_elmt);     ignod(6) = ibool(NGLLX,1,NGLLZ-1,i_elmt)
    ignod(7) = ibool(1,NGLLY,NGLLZ-1,i_elmt); ignod(8) = ibool(NGLLX,NGLLY,NGLLZ-1,i_elmt);

    coordinf(:,1) = xstore(ignod)
    coordinf(:,2) = ystore(ignod)
    coordinf(:,3) = zstore(ignod)

    ! Zienkiewicz infinite coordinates
    coordinf(5,:) = x0 + real(gaminf*(coordinf(1,:)-x0),kind=CUSTOM_REAL)
    coordinf(6,:) = x0 + real(gaminf*(coordinf(2,:)-x0),kind=CUSTOM_REAL)
    coordinf(7,:) = x0 + real(gaminf*(coordinf(3,:)-x0),kind=CUSTOM_REAL)
    coordinf(8,:) = x0 + real(gaminf*(coordinf(4,:)-x0),kind=CUSTOM_REAL)

    ! Point X0 (Pole)
    coordinf(1,:) = x0; coordinf(2,:) = x0
    coordinf(3,:) = x0; coordinf(4,:) = x0

    egdof = ibool1(:,i_elmt) !reshape(ibool1(:,:,:,i_elmt),(/NGLLCUBE_INF/))
    kmat = zero
    do i = 1,nipinf
      jac = real(matmul(dshape_infinite(:,i,:),coordinf),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gl(:,i,:)),kind=CUSTOM_REAL)
      kmat = kmat + real(matmul(transpose(deriv),deriv)*detjac*GLw(i),kind=CUSTOM_REAL)
    enddo
    storekmat(:,:,i_elmt) = kmat
    do k = 1,NGLLCUBE_INF
      dprecon(egdof(k)) = dprecon(egdof(k))+kmat(k,k)
    enddo
  enddo

  end subroutine poisson_stiffnessINF3

!
!===========================================
!

  subroutine compute_poisson_rhoload()

  use constants, only: IREGION_INNER_CORE,IREGION_OUTER_CORE,IREGION_CRUST_MANTLE

  use constants_solver, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NGLOB_INNER_CORE,NGLOB_OUTER_CORE, &
    NGLOB_CRUST_MANTLE

  use specfem_par_crustmantle, only: ibool_crust_mantle
  use specfem_par_outercore, only: ibool_outer_core
  use specfem_par_innercore, only: ibool_inner_core

  use specfem_par_full_gravity, only: load, &
    storerhojw_cm,gdof_cm, &
    storerhojw_oc,gdof_oc, &
    storerhojw_ic,gdof_ic

  implicit none
  real(kind=CUSTOM_REAL) :: load_ic(NGLOB_INNER_CORE),load_oc(NGLOB_OUTER_CORE), &
                            load_cm(NGLOB_CRUST_MANTLE)

  !note: NGLLCUBE = NGLLX*NGLLY*NGLLZ

  load = 0.0_CUSTOM_REAL

  ! crust-mantle
  call poisson_load_onlyrhoFAST(IREGION_CRUST_MANTLE,nspec_crust_mantle, &
                                nglob_crust_mantle,ibool_crust_mantle,storerhojw_cm,load_cm)

  ! outer core
  call poisson_load_onlyrhoFAST(IREGION_OUTER_CORE,nspec_outer_core, &
                                nglob_outer_core,ibool_outer_core,storerhojw_oc,load_oc)

  ! inner core
  call poisson_load_onlyrhoFAST(IREGION_INNER_CORE,nspec_inner_core, &
                                nglob_inner_core,ibool_inner_core,storerhojw_ic,load_ic)

  ! infinite
  ! this region has no contrbution

  ! assemble across the regions but not MPI
  ! crust_mantle
  load(gdof_cm) = load(gdof_cm) + load_cm

  ! outer core
  load(gdof_oc) = load(gdof_oc) + load_oc

  ! inner core
  load(gdof_ic) = load(gdof_ic) + load_ic

  ! infinite
  ! there no contribution from infinite region

  load(0) = 0.0_CUSTOM_REAL

  end subroutine compute_poisson_rhoload

!
!===========================================
!

  subroutine compute_poisson_rhoload3()

  use constants, only: IREGION_INNER_CORE,IREGION_OUTER_CORE, &
    IREGION_CRUST_MANTLE,NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use constants_solver, only: NSPEC_INNER_CORE,NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE

  use specfem_par_full_gravity, only: load1,nnode_ic1,nnode_oc1,nnode_cm1, &
    storerhojw_cm1,gdof_cm1,inode_elmt_cm1, &
    storerhojw_oc1,gdof_oc1,inode_elmt_oc1, &
    storerhojw_ic1,gdof_ic1,inode_elmt_ic1

  implicit none
  real(kind=CUSTOM_REAL) :: load_ic(nnode_ic1),load_oc(nnode_oc1), &
                            load_cm(nnode_cm1)

  !note: NGLLCUBE = NGLLX*NGLLY*NGLLZ

  load1 = 0.0_CUSTOM_REAL

  !! crust-mantle
  !call poisson_load_onlyrho3(1,nspec_crust_mantle, &
  !                           nglob_crust_mantle,ibool_crust_mantle,xstore0_crust_mantle,ystore0_crust_mantle, &
  !                           zstore0_crust_mantle,rhostore_crust_mantle,nnode_cm1,inode_elmt_cm1,load_cm)
  !
  !! outer core
  !call poisson_load_onlyrho3(2,nspec_outer_core, &
  !                           nglob_outer_core,ibool_outer_core,xstore0_outer_core,ystore0_outer_core, &
  !                           zstore0_outer_core,rhostore_outer_core,nnode_oc1,inode_elmt_oc1,load_oc)
  !
  !! inner core
  !call poisson_load_onlyrho3(3,nspec_inner_core, &
  !                           nglob_inner_core,ibool_inner_core,xstore0_inner_core,ystore0_inner_core, &
  !                           zstore0_inner_core,rhostore_inner_core,nnode_ic1,inode_elmt_ic1,load_ic)

  ! crust-mantle
  call poisson_load_onlyrhoFAST3(IREGION_CRUST_MANTLE,nspec_crust_mantle, &
                                 storerhojw_cm1,nnode_cm1,inode_elmt_cm1, &
                                 load_cm)

  ! outer core
  call poisson_load_onlyrhoFAST3(IREGION_OUTER_CORE,nspec_outer_core, &
                                 storerhojw_oc1,nnode_oc1,inode_elmt_oc1, &
                                 load_oc)

  ! inner core
  call poisson_load_onlyrhoFAST3(IREGION_INNER_CORE,nspec_inner_core, &
                                 storerhojw_ic1,nnode_ic1,inode_elmt_ic1, &
                                 load_ic)

  ! infinite
  ! this region has no contribution

  ! assemble across the regions but not MPI
  ! crust_mantle
  load1(gdof_cm1) = load1(gdof_cm1) + load_cm

  ! outer core
  load1(gdof_oc1) = load1(gdof_oc1) + load_oc

  ! inner core
  load1(gdof_ic1) = load1(gdof_ic1) + load_ic

  ! infinite
  ! there is no contribution from infinite region

  load1(0) = 0.0_CUSTOM_REAL

  end subroutine compute_poisson_rhoload3

!
!===========================================
!

  subroutine compute_grav_kl1_load(component)

  ! Computes load for the first gravity kernel that must be solved using SIEM

  use constants, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF
  use constants_solver, only: NSPEC_CRUST_MANTLE

  use specfem_par_full_gravity, only: load1, nnode_cm1, &
    storejw_cm1,gdof_cm1,inode_elmt_cm1, &
    rho1siem_kl_crust_mantle

  use siem_gll_library

  implicit none

  ! IO:
  integer :: component ! the component of the kernel to be calculated

  integer,parameter :: ndim = 3,ngnod = 8
  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  real(kind=CUSTOM_REAL) :: load_cm(nnode_cm1), evalue(NGLLCUBE_INF), eload(NGLLCUBE_INF)
  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF)

  ! Based off of poisson_load_solid3FAST
  load1 = 0.0_CUSTOM_REAL
  load_cm = zero
  do i_elmt = 1,nspec_crust_mantle
    evalue = reshape(rho1siem_kl_crust_mantle(component,1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    eload = zero
    do i = 1,NGLLCUBE_INF
      eload = eload + storejw_cm1(i,i_elmt)*evalue(i)
    enddo
    egdof = inode_elmt_cm1(:,i_elmt)
    load_cm(egdof) = load_cm(egdof)+eload
  enddo

  ! assemble across the regions but not MPI
  ! crust_mantle -  multiply by 4*PI*G! or scaled
  load1(gdof_cm1) = load1(gdof_cm1) + 4.0_CUSTOM_REAL*load_cm
  load1(0) = 0.0_CUSTOM_REAL

  end subroutine compute_grav_kl1_load

!
!===========================================
!

  subroutine compute_grav_kl2_load(icomp,jcomp)

  ! Computes load for the first gravity kernel that must be solved using SIEM

  use constants, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use constants_solver, only: NSPEC_CRUST_MANTLE

  use specfem_par_full_gravity, only: load1, nnode_cm1, &
    storejw_cm1,gdof_cm1,inode_elmt_cm1, &
    rho2siem_kl_crust_mantle

  use siem_gll_library

  implicit none
  ! IO variables
  integer :: icomp,jcomp, i_elmt, i
  !Local
  integer,parameter :: ndim = 3,ngnod = 8
  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL
  integer :: egdof(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: load_cm(nnode_cm1), evalue(NGLLCUBE_INF), eload(NGLLCUBE_INF)

  load1 = 0.0_CUSTOM_REAL

  ! Based off of poisson_load_solid3FAST
  load_cm = zero
  do i_elmt = 1,nspec_crust_mantle
    evalue = reshape(rho2siem_kl_crust_mantle(icomp,jcomp,1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    eload = zero
    do i = 1,NGLLCUBE_INF
      eload = eload + storejw_cm1(i,i_elmt)*evalue(i)
    enddo
    egdof = inode_elmt_cm1(:,i_elmt)
    load_cm(egdof)=load_cm(egdof)+eload
  enddo
  ! assemble across the regions but not MPI
  ! crust_mantle -  multiply by 4*PI*G! or scaled
  load1(gdof_cm1) = load1(gdof_cm1) + 4.0_CUSTOM_REAL*load_cm
  load1(0) = 0.0_CUSTOM_REAL

  end subroutine compute_grav_kl2_load

!
!===========================================
!

  subroutine compute_poisson_load()

  use constants, only: IREGION_INNER_CORE,IREGION_CRUST_MANTLE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use constants_solver, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NGLOB_INNER_CORE,NGLOB_OUTER_CORE, &
    NGLOB_CRUST_MANTLE

  use specfem_par, only: deltat,it,DT,t0,scale_t_inv,two_omega_earth

  use specfem_par_crustmantle, only: displ_crust_mantle,ibool_crust_mantle
  use specfem_par_outercore, only: displ_outer_core,ibool_outer_core
  use specfem_par_innercore, only: displ_inner_core,ibool_inner_core

  use specfem_par_full_gravity, only: load, &
    storederiv_cm,storerhojw_cm,gdof_cm, &
    storederiv_oc,storerhojw_oc,gdof_oc, &
    storederiv_ic,storerhojw_ic,gdof_ic, &
    A_array_rotationL,B_array_rotationL

  implicit none
  real(kind=CUSTOM_REAL) :: time
  real(kind=CUSTOM_REAL) :: load_ic(NGLOB_INNER_CORE),load_oc(NGLOB_OUTER_CORE), &
  load_cm(NGLOB_CRUST_MANTLE)

  !NGLLCUBE = NGLLX*NGLLY*NGLLZ

  load = 0.0_CUSTOM_REAL

  ! inner core
  !call poisson_load_solid(IREGION_INNER_CORE,nspec_inner_core, &
  !                        nglob_inner_core,ibool_inner_core,xstore0_inner_core,ystore0_inner_core, &
  !                        zstore0_inner_core,rhostore_inner_core,displ_inner_core,load_ic)

  call poisson_load_solidFAST(IREGION_INNER_CORE,nspec_inner_core, &
                              nglob_inner_core,ibool_inner_core,storederiv_ic,storerhojw_ic,displ_inner_core, &
                              load_ic)

  ! crust-mantle
  !call poisson_load_solid(IREGION_CRUST_MANTLE,nspec_crust_mantle, &
  !                        nglob_crust_mantle,ibool_crust_mantle,xstore0_crust_mantle,ystore0_crust_mantle, &
  !                        zstore0_crust_mantle,rhostore_crust_mantle,displ_crust_mantle,load_cm)

  call poisson_load_solidFAST(IREGION_CRUST_MANTLE,nspec_crust_mantle, &
                              nglob_crust_mantle,ibool_crust_mantle,storederiv_cm,storerhojw_cm, &
                              displ_crust_mantle,load_cm)

  ! outer core
  ! time
  time = real((dble(it-1)*DT-t0)*scale_t_inv,kind=CUSTOM_REAL)

  !call poisson_load_fluid(time,deltat,two_omega_earth,NSPEC_OUTER_CORE, &
  !                        NGLOB_OUTER_CORE,A_array_rotation,B_array_rotation,rhostore_outer_core, &
  !                        displ_outer_core,load_oc)

  !call poisson_load_fluidNEW(nspec_outer_core,nglob_outer_core,ibool_outer_core, &
  !                           xstore0_outer_core,ystore0_outer_core,zstore0_outer_core,rhostore_outer_core, &
  !                           time,deltat,two_omega_earth,A_array_rotationL,B_array_rotationL, &
  !                           displ_outer_core,load_oc)

  call poisson_load_fluidNEWFAST(nspec_outer_core,nglob_outer_core,ibool_outer_core, &
                                 storederiv_oc,storerhojw_oc,time,deltat,two_omega_earth,A_array_rotationL, &
                                 B_array_rotationL,displ_outer_core,load_oc)

  ! infinite
  ! this region has no contrbution

  ! assemble across the regions but not MPI
  ! crust_mantle
  load(gdof_cm) = load(gdof_cm) + load_cm

  ! outer core
  load(gdof_oc) = load(gdof_oc) + load_oc

  ! inner core
  load(gdof_ic) = load(gdof_ic) + load_ic

  ! infinite
  ! there no contribution from infinite region

  load(0) = 0.0_CUSTOM_REAL

  end subroutine compute_poisson_load

!
!===========================================
!

  subroutine compute_poisson_load3()

  use constants, only: IREGION_INNER_CORE,IREGION_CRUST_MANTLE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use constants_solver, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NGLOB_INNER_CORE,NGLOB_OUTER_CORE, &
    NGLOB_CRUST_MANTLE

  use specfem_par, only: deltat,it,DT,t0,scale_t_inv,two_omega_earth

  use specfem_par_crustmantle, only: displ_crust_mantle,ibool_crust_mantle
  use specfem_par_outercore, only: displ_outer_core,ibool_outer_core
  use specfem_par_innercore, only: displ_inner_core,ibool_inner_core

  use specfem_par_full_gravity, only: load1,nnode_ic1,nnode_oc1,nnode_cm1, &
    storederiv_cm1,storerhojw_cm1,gdof_cm1,inode_elmt_cm1, &
    storederiv_oc1,storerhojw_oc1,gdof_oc1,inode_elmt_oc1, &
    storederiv_ic1,storerhojw_ic1,gdof_ic1,inode_elmt_ic1, &
    A_array_rotationL3,B_array_rotationL3

  implicit none
  real(kind=CUSTOM_REAL) :: time
  real(kind=CUSTOM_REAL) :: load_ic(nnode_ic1),load_oc(nnode_oc1), &
  load_cm(nnode_cm1)

  !note: NGLLCUBE = NGLLX*NGLLY*NGLLZ

  load1 = 0.0_CUSTOM_REAL
  ! inner core
  call poisson_load_solid3FAST(IREGION_INNER_CORE,nspec_inner_core, &
                               nglob_inner_core,ibool_inner_core,storederiv_ic1,storerhojw_ic1, &
                               displ_inner_core,nnode_ic1,inode_elmt_ic1,load_ic)

  ! crust-mantle
  call poisson_load_solid3FAST(IREGION_CRUST_MANTLE,nspec_crust_mantle, &
                               nglob_crust_mantle,ibool_crust_mantle,storederiv_cm1,storerhojw_cm1, &
                               displ_crust_mantle,nnode_cm1,inode_elmt_cm1,load_cm)

  ! outer core
  ! time
  time = real((dble(it-1)*DT-t0)*scale_t_inv,kind=CUSTOM_REAL)

  call poisson_load_fluidNEW3FAST(nspec_outer_core,nglob_outer_core, &
                                  ibool_outer_core,storederiv_oc1,storerhojw_oc1,time,deltat,two_omega_earth, &
                                  A_array_rotationL3,B_array_rotationL3,displ_outer_core,nnode_oc1,inode_elmt_oc1, &
                                  load_oc)

  ! infinite
  ! this region has no contrbution

  ! assemble across the regions but not MPI
  ! crust_mantle
  load1(gdof_cm1) = load1(gdof_cm1) + load_cm

  ! outer core
  load1(gdof_oc1) = load1(gdof_oc1) + load_oc

  ! inner core
  load1(gdof_ic1) = load1(gdof_ic1) + load_ic

  load1(0) = 0.0_CUSTOM_REAL

  end subroutine compute_poisson_load3

!
!===========================================
!

  subroutine compute_backward_poisson_load3()

  use constants, only: IREGION_INNER_CORE,IREGION_CRUST_MANTLE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use constants_solver, only: NSPEC_INNER_CORE, &
    NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE,NGLOB_INNER_CORE,NGLOB_OUTER_CORE, &
    NGLOB_CRUST_MANTLE

  use specfem_par, only: b_deltat,it,DT,t0,scale_t_inv,two_omega_earth

  use specfem_par_crustmantle, only: b_displ_crust_mantle,ibool_crust_mantle
  use specfem_par_outercore, only: b_displ_outer_core,ibool_outer_core
  use specfem_par_innercore, only: b_displ_inner_core,ibool_inner_core

  use specfem_par_full_gravity, only: b_load1,nnode_ic1,nnode_oc1,nnode_cm1, &
    storederiv_cm1,storerhojw_cm1,gdof_cm1,inode_elmt_cm1, &
    storederiv_oc1,storerhojw_oc1,gdof_oc1,inode_elmt_oc1, &
    storederiv_ic1,storerhojw_ic1,gdof_ic1,inode_elmt_ic1, &
    b_A_array_rotationL3,b_B_array_rotationL3

  implicit none
  real(kind=CUSTOM_REAL) :: time
  ! local loads
  real(kind=CUSTOM_REAL) :: b_load_ic(nnode_ic1),b_load_oc(nnode_oc1), &
  b_load_cm(nnode_cm1)

  !note: NGLLCUBE = NGLLX*NGLLY*NGLLZ

  b_load1 = 0.0_CUSTOM_REAL

  ! inner core
  call poisson_load_solid3FAST(IREGION_INNER_CORE,nspec_inner_core, &
                               nglob_inner_core,ibool_inner_core,storederiv_ic1,storerhojw_ic1, &
                               b_displ_inner_core,nnode_ic1,inode_elmt_ic1,b_load_ic)

  ! crust-mantle
  call poisson_load_solid3FAST(IREGION_CRUST_MANTLE,nspec_crust_mantle, &
                               nglob_crust_mantle,ibool_crust_mantle,storederiv_cm1,storerhojw_cm1, &
                               b_displ_crust_mantle,nnode_cm1,inode_elmt_cm1,b_load_cm)

  ! outer core
  ! time
  time = real((dble(it-1)*DT-t0)*scale_t_inv,kind=CUSTOM_REAL)

  ! WE - unsure on whether to use b_deltat or deltat
  ! Note that two_omega_earth is already reversed
  call poisson_load_fluidNEW3FAST(nspec_outer_core,nglob_outer_core, &
                                  ibool_outer_core,storederiv_oc1,storerhojw_oc1,time,b_deltat,two_omega_earth, &
                                  b_A_array_rotationL3, b_B_array_rotationL3,b_displ_outer_core,nnode_oc1,inode_elmt_oc1, &
                                  b_load_oc)

  ! infinite
  ! this region has no contrbution


  ! assemble across the regions but not MPI
  ! crust_mantle
  b_load1(gdof_cm1) = b_load1(gdof_cm1) + b_load_cm

  ! outer core
  b_load1(gdof_oc1) = b_load1(gdof_oc1) + b_load_oc

  ! inner core
  b_load1(gdof_ic1) = b_load1(gdof_ic1) + b_load_ic

  b_load1(0) = 0.0_CUSTOM_REAL

  end subroutine compute_backward_poisson_load3

!
!===========================================
!

  subroutine poisson_load_solid(iregion,nelmt,nnode,ibool, &
                                xstore,ystore,zstore,rhostore,disp,load)

  use constants_solver, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NDIM, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(in) :: rhostore(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: disp(NDIM,nnode)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode)

  integer,parameter :: ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE),ignod(ngnod) !dnx,dny,dnz,
  real(kind=CUSTOM_REAL) :: detjac,rho(NGLLCUBE)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ)

  real(kind=kdble) :: gll_weights(NGLLCUBE),gll_points(ndim,NGLLCUBE), &
                      dshape_hex8(ndim,ngnod,NGLLCUBE)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),deriv(ndim,NGLLCUBE),jac(ndim,ndim)

  real(kind=CUSTOM_REAL) :: edisp(NDIM,NGLLCUBE),eload(NGLLCUBE)

  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gll
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gll

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE), &
           dlagrange_gll(ndim,NGLLCUBE,NGLLCUBE))
  lagrange_gll(:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  call zwgljd(xigll,wxgll,NGLLX,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = 0.0_CUSTOM_REAL
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    rho = reshape(rhostore(:,:,:,i_elmt),(/NGLLCUBE/))
    edisp = disp(:,egdof)

    eload = zero
    do i = 1,NGLLCUBE
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      eload = eload + real(rho(i)*dotmat(ndim,NGLLCUBE,deriv,edisp)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eload
  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_solid

!
!===========================================
!

  subroutine poisson_load_solid3(iregion,nelmt,nnode,ibool, &
                                 xstore,ystore,zstore,rhostore,disp,nnode1,ibool1,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NDIM, &
    NGLLX,NGLLY,NGLLZ, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(in) :: rhostore(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: disp(NDIM,nnode)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF),egdof1(NGLLCUBE_INF),ignod(ngnod) !dnx,dny,dnz,
  real(kind=CUSTOM_REAL) :: detjac,rho(NGLLCUBE_INF)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX_INF),wxgll(NGLLX_INF),etagll(NGLLY_INF),wygll(NGLLY_INF), &
  zetagll(NGLLZ_INF),wzgll(NGLLZ_INF)

  real(kind=kdble) :: gll_weights(NGLLCUBE_INF),gll_points(ndim,NGLLCUBE_INF),lagrange_gll(NGLLCUBE_INF,NGLLCUBE_INF), &
                      dlagrange_gll(ndim,NGLLCUBE_INF,NGLLCUBE_INF),dshape_hex8(ndim,ngnod,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),deriv(ndim,NGLLCUBE_INF),jac(ndim,ndim)

  real(kind=CUSTOM_REAL) :: edisp(NDIM,NGLLCUBE_INF),eload(NGLLCUBE_INF)

  call zwgljd(xigll,wxgll,NGLLX_INF,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY_INF,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ_INF,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-1

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = 0.0_CUSTOM_REAL
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof1 = reshape(ibool(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    egdof = ibool1(:,i_elmt)
    rho = reshape(rhostore(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    edisp = disp(:,egdof1)

    eload = zero
    do i = 1,NGLLCUBE_INF
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      eload = eload + real(rho(i)*dotmat(ndim,NGLLCUBE_INF,deriv,edisp)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eload
  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_solid3

!
!===========================================
!

!TODO: transpose can be avoided here doing so in preintegrate

  subroutine poisson_load_solidFAST(iregion,nelmt,nnode,ibool, &
                                    storederiv,storerhojw,disp,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NDIM, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_math_library, only: dotmat
  use siem_gll_library

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  storederiv(NDIM,NGLLCUBE,NGLLCUBE,nelmt),storerhojw(NGLLCUBE,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: disp(NDIM,nnode)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode)

  integer,parameter :: ngnod = 8
  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE)
  real(kind=CUSTOM_REAL) :: deriv(ndim,NGLLCUBE),edisp(NDIM,NGLLCUBE),eload(NGLLCUBE)
  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    edisp = disp(:,egdof)

    eload = zero
    do i = 1,NGLLCUBE
      deriv = storederiv(:,:,i,i_elmt)
      eload = eload + storerhojw(i,i_elmt)*matmul(transpose(deriv),edisp(:,i))
    enddo
    load(egdof) = load(egdof) + eload
  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_solidFAST

!
!===========================================
!

  subroutine poisson_load_solid3FAST1(iregion,nelmt,nnode, &
                                      ibool,storederiv,storerhojw,disp,nnode1,ibool1,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NDIM, &
    NGLLX,NGLLY,NGLLZ, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_innercore, only: idoubling_inner_core

  use specfem_par_full_gravity, only: lagrange_gll1

  use siem_math_library, only: dotmat
  use siem_gll_library

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  storederiv(NDIM,NGLLCUBE_INF,NGLLCUBE_INF,nelmt), &
                                        storerhojw(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: disp(NDIM,nnode)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ngnod = 8
  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF),egdof1(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: divs,deriv(ndim,NGLLCUBE_INF),edisp(NDIM,NGLLCUBE_INF),eload(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    egdof1 = reshape(ibool(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    egdof = ibool1(:,i_elmt)
    edisp = disp(:,egdof1)

    eload = zero
    do i = 1,NGLLCUBE_INF
      deriv = storederiv(:,:,i,i_elmt)
      divs = dot_product(deriv(1,:),edisp(1,:)) + dot_product(deriv(2,:),edisp(2,:)) + &
             dot_product(deriv(3,:),edisp(3,:)) ! rho should be included here
      eload = eload + real(storerhojw(i,i_elmt)*divs*lagrange_gll1(i,:),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eload
  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_solid3FAST1

!
!===========================================
!

  subroutine poisson_load_solid3FAST(iregion,nelmt,nnode,ibool, &
                                     storederiv,storerhojw,disp,nnode1,ibool1,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NDIM, &
    NGLLX,NGLLY,NGLLZ, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_math_library, only: dotmat
  use siem_gll_library

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  storederiv(NDIM,NGLLCUBE_INF,NGLLCUBE_INF,nelmt),storerhojw(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: disp(NDIM,nnode)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ngnod = 8
  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF),egdof1(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: deriv(ndim,NGLLCUBE_INF),edisp(NDIM,NGLLCUBE_INF),eload(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    egdof1 = reshape(ibool(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    egdof = ibool1(:,i_elmt)
    edisp = disp(:,egdof1)

    eload = zero
    do i = 1,NGLLCUBE_INF
      deriv = storederiv(:,:,i,i_elmt)
      eload = eload + storerhojw(i,i_elmt)*matmul(transpose(deriv),edisp(:,i))
    enddo
    load(egdof) = load(egdof) + eload
  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_solid3FAST

!
!===========================================
!

  subroutine poisson_load_fluid(time,deltat,two_omega_earth,NSPEC,NGLOB, &
                                A_array_rotation,B_array_rotation,rhostore,displfluid,load)

  use constants, only: NGLLX,NGLLY,NGLLZ
  use constants_solver, only: ROTATION_VAL

  use specfem_par, only: hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy, &
    hprimewgll_zz,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

  use specfem_par_outercore, only: &
    xix => xix_outer_core, &
    xiy => xiy_outer_core, &
    xiz => xiz_outer_core, &
    etax => etax_outer_core, &
    etay => etay_outer_core, &
    etaz => etaz_outer_core, &
    gammax => gammax_outer_core, &
    gammay => gammay_outer_core, &
    gammaz => gammaz_outer_core, &
    ibool => ibool_outer_core

  implicit none

  integer :: NSPEC,NGLOB

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) :: time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: A_array_rotation,B_array_rotation,rhostore

  ! displacement and acceleration
  real(kind=CUSTOM_REAL),dimension(NGLOB) :: displfluid,load

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t, &
    A_rotation,B_rotation,ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer :: ispec,i,j,k,l

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l,sum_terms

  load = 0.0_CUSTOM_REAL

  do ispec = 1,NSPEC
    ! only compute element which belong to current phase (inner or outer elements)
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l = 1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + displfluid(ibool(l,j,k,ispec)) * hprime_xx(i,l)
            tempx2l = tempx2l + displfluid(ibool(i,l,k,ispec)) * hprime_yy(j,l)
            tempx3l = tempx3l + displfluid(ibool(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          ! get derivatives of velocity potential with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! compute the jacobian
          jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl)       &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl)                          &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))

          dpotentialdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          dpotentialdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          dpotentialdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          ! compute contribution of rotation and add to gradient of potential
          ! this term has no Z component
          if (ROTATION_VAL) then

            ! store the source for the Euler scheme for A_rotation and B_rotation
            two_omega_deltat = deltat * two_omega_earth

            cos_two_omega_t = cos(two_omega_earth*time)
            sin_two_omega_t = sin(two_omega_earth*time)

            ! time step deltat of Euler scheme is included in the source
            source_euler_A(i,j,k) = two_omega_deltat                              &
                  *(cos_two_omega_t*dpotentialdyl+sin_two_omega_t*dpotentialdxl)
            source_euler_B(i,j,k) = two_omega_deltat                              &
                  *(sin_two_omega_t*dpotentialdyl-cos_two_omega_t*dpotentialdxl)

            A_rotation = A_array_rotation(i,j,k,ispec)
            B_rotation = B_array_rotation(i,j,k,ispec)

            ux_rotation =   A_rotation*cos_two_omega_t+B_rotation*sin_two_omega_t
            uy_rotation = - A_rotation*sin_two_omega_t+B_rotation*cos_two_omega_t

            dpotentialdx_with_rot = dpotentialdxl + ux_rotation
            dpotentialdy_with_rot = dpotentialdyl + uy_rotation

          else
            dpotentialdx_with_rot = dpotentialdxl
            dpotentialdy_with_rot = dpotentialdyl

          endif  ! end of section with rotation

          tempx1(i,j,k) = jacobianl*(xixl*dpotentialdx_with_rot+xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
          tempx2(i,j,k) = jacobianl*(etaxl*dpotentialdx_with_rot+etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
          tempx3(i,j,k) = jacobianl*(gammaxl*dpotentialdx_with_rot+gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

        enddo
      enddo
    enddo

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0.0_CUSTOM_REAL
          tempx2l = 0.0_CUSTOM_REAL
          tempx3l = 0.0_CUSTOM_REAL

          do l = 1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + tempx1(l,j,k) * hprimewgll_xx(l,i)
            tempx2l = tempx2l + tempx2(i,l,k) * hprimewgll_yy(l,j)
            tempx3l = tempx3l + tempx3(i,j,l) * hprimewgll_zz(l,k)
          enddo

          ! sum contributions from each element to the global mesh and add gravity term
          sum_terms = -(wgllwgll_yz(j,k)*tempx1l+wgllwgll_xz(i,k)*tempx2l+wgllwgll_xy(i,j)*tempx3l)

          load(ibool(i,j,k,ispec)) = load(ibool(i,j,k,ispec)) + rhostore(i,j,k,ispec)*sum_terms

        enddo
      enddo
    enddo

  enddo ! ispec = 1,NSPEC spectral element loop

  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_fluid

!
!===========================================
!

  subroutine poisson_load_fluid3(time,deltat,two_omega_earth,NSPEC,NGLOB, &
                                 A_array_rotation,B_array_rotation,rhostore,displfluid,load)

  use constants, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE_INF
  use constants_solver, only: ROTATION_VAL

  use specfem_par, only: hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx, &
    hprimewgll_yy,hprimewgll_zz,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz

  use specfem_par_outercore, only: &
    xix => xix_outer_core, &
    xiy => xiy_outer_core, &
    xiz => xiz_outer_core, &
    etax => etax_outer_core, &
    etay => etay_outer_core, &
    etaz => etaz_outer_core, &
    gammax => gammax_outer_core, &
    gammay => gammay_outer_core, &
    gammaz => gammaz_outer_core, &
    ibool => ibool_outer_core

  implicit none

  integer,intent(in) :: NSPEC,NGLOB
  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) :: time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: A_array_rotation,B_array_rotation,rhostore

  ! displacement and acceleration
  real(kind=CUSTOM_REAL),dimension(NGLOB) :: displfluid,load

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) two_omega_deltat,cos_two_omega_t,sin_two_omega_t, &
    A_rotation,B_rotation,ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer :: ispec,i,j,k,l

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l,sum_terms

  load = 0.0_CUSTOM_REAL

  do ispec = 1,NSPEC
    ! only compute element which belong to current phase (inner or outer elements)
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l = 1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + displfluid(ibool(l,j,k,ispec)) * hprime_xx(i,l)
            tempx2l = tempx2l + displfluid(ibool(i,l,k,ispec)) * hprime_yy(j,l)
            tempx3l = tempx3l + displfluid(ibool(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          ! get derivatives of velocity potential with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! compute the jacobian
          jacobianl = 1.0_CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl)       &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl)                          &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))

          dpotentialdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          dpotentialdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          dpotentialdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          ! compute contribution of rotation and add to gradient of potential
          ! this term has no Z component
          if (ROTATION_VAL) then

            ! store the source for the Euler scheme for A_rotation and B_rotation
            two_omega_deltat = deltat * two_omega_earth

            cos_two_omega_t = cos(two_omega_earth*time)
            sin_two_omega_t = sin(two_omega_earth*time)

            ! time step deltat of Euler scheme is included in the source
            source_euler_A(i,j,k) = two_omega_deltat &
                  * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
            source_euler_B(i,j,k) = two_omega_deltat &
                  * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

            A_rotation = A_array_rotation(i,j,k,ispec)
            B_rotation = B_array_rotation(i,j,k,ispec)

            ux_rotation =   A_rotation*cos_two_omega_t+B_rotation*sin_two_omega_t
            uy_rotation = - A_rotation*sin_two_omega_t+B_rotation*cos_two_omega_t

            dpotentialdx_with_rot = dpotentialdxl + ux_rotation
            dpotentialdy_with_rot = dpotentialdyl + uy_rotation

          else
            dpotentialdx_with_rot = dpotentialdxl
            dpotentialdy_with_rot = dpotentialdyl

          endif  ! end of section with rotation

          tempx1(i,j,k) = jacobianl*(xixl*dpotentialdx_with_rot+xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
          tempx2(i,j,k) = jacobianl*(etaxl*dpotentialdx_with_rot+etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
          tempx3(i,j,k) = jacobianl*(gammaxl*dpotentialdx_with_rot+gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

        enddo
      enddo
    enddo

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0.0_CUSTOM_REAL
          tempx2l = 0.0_CUSTOM_REAL
          tempx3l = 0.0_CUSTOM_REAL

          do l = 1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + tempx1(l,j,k) * hprimewgll_xx(l,i)
            tempx2l = tempx2l + tempx2(i,l,k) * hprimewgll_yy(l,j)
            tempx3l = tempx3l + tempx3(i,j,l) * hprimewgll_zz(l,k)
          enddo

          ! sum contributions from each element to the global mesh and add gravity term
          sum_terms = - (wgllwgll_yz(j,k)*tempx1l + wgllwgll_xz(i,k)*tempx2l + wgllwgll_xy(i,j)*tempx3l)

          load(ibool(i,j,k,ispec)) = load(ibool(i,j,k,ispec)) + rhostore(i,j,k,ispec)*sum_terms

        enddo
      enddo
    enddo

  enddo ! ispec = 1,NSPEC spectral element loop

  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_fluid3

!
!===========================================
!

  subroutine poisson_load_fluidNEW(nelmt,nnode,ibool,xstore,ystore,zstore, &
                                   rhostore,time,deltat,two_omega_earth,A_array_rot,B_array_rot,dispf,load)

  use constants, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE
  use constants_solver, only: ROTATION_VAL

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: xstore(nnode),ystore(nnode),zstore(nnode), &
                                       rhostore(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE,nelmt),intent(inout) :: A_array_rot,B_array_rot
  real(kind=CUSTOM_REAL),intent(in) :: dispf(1,nnode) !\chi
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac(NGLLCUBE),rho(NGLLCUBE)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ)

  real(kind=kdble) :: gll_weights(NGLLCUBE),gll_points(ndim,NGLLCUBE), &
                      dshape_hex8(ndim,ngnod,NGLLCUBE)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),jac(ndim,ndim)

  real(kind=CUSTOM_REAL) :: echi(NGLLCUBE,1),edisp(ndim,NGLLCUBE),eload(NGLLCUBE),gradchi(ndim,1)
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rot,B_rot,ux_rot,uy_rot
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE) :: source_euler_A,source_euler_B

  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gll
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gll

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: deriv

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE), &
           dlagrange_gll(ndim,NGLLCUBE,NGLLCUBE))
  lagrange_gll(:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  allocate(deriv(ndim,NGLLCUBE,NGLLCUBE))
  deriv(:,:,:) = 0.0_CUSTOM_REAL

  call zwgljd(xigll,wxgll,NGLLX,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-1

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = zero
  do i_elmt = 1,nelmt
    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    rho = reshape(rhostore(:,:,:,i_elmt),(/NGLLCUBE/))
    echi(:,1) = dispf(1,egdof)

    ! compute diaplacement
    do i = 1,NGLLCUBE
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac(i) = determinant(jac)
      call invert(jac)
      deriv(:,:,i) = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      gradchi = matmul(deriv(:,:,i),echi)

      edisp(:,i) = gradchi(:,1)
      ! compute contribution of rotation and add to gradient of potential
      ! this term has no Z component
      if (ROTATION_VAL) then
        ! store the source for the Euler scheme for A_rotation and B_rotation
        two_omega_deltat = deltat*two_omega_earth

        cos_two_omega_t = cos(two_omega_earth*time)
        sin_two_omega_t = sin(two_omega_earth*time)

        ! time step deltat of Euler scheme is included in the source
        source_euler_A(i) = two_omega_deltat*(cos_two_omega_t*gradchi(2,1)+         &
                              sin_two_omega_t*gradchi(1,1))
        source_euler_B(i) = two_omega_deltat*(sin_two_omega_t*gradchi(2,1)-         &
                              cos_two_omega_t*gradchi(1,1))

        A_rot = A_array_rot(i,i_elmt)
        B_rot = B_array_rot(i,i_elmt)

        ux_rot =  A_rot*cos_two_omega_t+B_rot*sin_two_omega_t
        uy_rot = -A_rot*sin_two_omega_t+B_rot*cos_two_omega_t

        edisp(1,i) = edisp(1,i) + ux_rot
        edisp(2,i) = edisp(2,i) + uy_rot
      endif
    enddo

    ! integration
    eload = zero
    do i = 1,NGLLCUBE
      eload = eload + real(rho(i)*dotmat(ndim,NGLLCUBE,deriv(:,:,i),edisp)*detjac(i)*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eload

    ! update rotation term with Euler scheme
    if (ROTATION_VAL) then
      ! use the source saved above
      A_array_rot(:,i_elmt) = A_array_rot(:,i_elmt) + source_euler_A
      B_array_rot(:,i_elmt) = B_array_rot(:,i_elmt) + source_euler_B
    endif

  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_fluidNEW

!
!===========================================
!

  subroutine poisson_load_fluidNEW3(nelmt,nnode,ibool,xstore,ystore,zstore, &
                                    rhostore,time,deltat,two_omega_earth,A_array_rot,B_array_rot,dispf,nnode1, &
                                    ibool1,load)

  use constants, only: NGLLX,NGLLY,NGLLZ,NGLLCUBE,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF
  use constants_solver, only: ROTATION_VAL

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: xstore(nnode),ystore(nnode),zstore(nnode), &
                                       rhostore(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE_INF,nelmt),intent(inout) :: A_array_rot,B_array_rot
  real(kind=CUSTOM_REAL),intent(in) :: dispf(1,nnode) !\chi
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF),egdof1(NGLLCUBE_INF),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac(NGLLCUBE_INF),rho(NGLLCUBE_INF)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX_INF),wxgll(NGLLX_INF),etagll(NGLLY_INF),wygll(NGLLY_INF), &
                      zetagll(NGLLZ_INF),wzgll(NGLLZ_INF)

  real(kind=kdble) :: gll_weights(NGLLCUBE_INF),gll_points(ndim,NGLLCUBE_INF), &
                      lagrange_gll(NGLLCUBE_INF,NGLLCUBE_INF),dlagrange_gll(ndim,NGLLCUBE_INF,NGLLCUBE_INF), &
                      dshape_hex8(ndim,ngnod,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),deriv(ndim,NGLLCUBE_INF,NGLLCUBE_INF),jac(ndim,ndim)

  real(kind=CUSTOM_REAL) :: echi(NGLLCUBE_INF,1),edisp(ndim,NGLLCUBE_INF),eload(NGLLCUBE_INF),gradchi(ndim,1)
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rot,B_rot,ux_rot,uy_rot
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE_INF) :: source_euler_A,source_euler_B

  call zwgljd(xigll,wxgll,NGLLX_INF,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY_INF,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ_INF,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-1

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = zero
  do i_elmt = 1,nelmt
    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/))
    ! this is wrong
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof1 = reshape(ibool(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    egdof = ibool1(:,i_elmt)
    rho = reshape(rhostore(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    echi(:,1) = dispf(1,egdof1)

    ! compute diaplacement
    do i = 1,NGLLCUBE_INF
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac(i) = determinant(jac)
      call invert(jac)
      deriv(:,:,i) = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      gradchi = matmul(deriv(:,:,i),echi)

      edisp(:,i) = gradchi(:,1)
      ! compute contribution of rotation and add to gradient of potential
      ! this term has no Z component
      if (ROTATION_VAL) then
        ! store the source for the Euler scheme for A_rotation and B_rotation
        two_omega_deltat = deltat*two_omega_earth

        cos_two_omega_t = cos(two_omega_earth*time)
        sin_two_omega_t = sin(two_omega_earth*time)

        ! time step deltat of Euler scheme is included in the source
        source_euler_A(i) = two_omega_deltat*(cos_two_omega_t*gradchi(2,1)+         &
                              sin_two_omega_t*gradchi(1,1))
        source_euler_B(i) = two_omega_deltat*(sin_two_omega_t*gradchi(2,1)-         &
                              cos_two_omega_t*gradchi(1,1))

        A_rot = A_array_rot(i,i_elmt)
        B_rot = B_array_rot(i,i_elmt)

        ux_rot =  A_rot*cos_two_omega_t+B_rot*sin_two_omega_t
        uy_rot = -A_rot*sin_two_omega_t+B_rot*cos_two_omega_t

        edisp(1,i) = edisp(1,i) + ux_rot
        edisp(2,i) = edisp(2,i) + uy_rot
      endif
    enddo

    ! integration
    eload = zero
    do i = 1,NGLLCUBE_INF
      eload = eload + real(rho(i)*dotmat(ndim,NGLLCUBE_INF,deriv(:,:,i),edisp)*detjac(i)*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eload

    ! update rotation term with Euler scheme
    if (ROTATION_VAL) then
      ! use the source saved above
      A_array_rot(:,i_elmt) = A_array_rot(:,i_elmt) + source_euler_A
      B_array_rot(:,i_elmt) = B_array_rot(:,i_elmt) + source_euler_B
    endif

  enddo

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_fluidNEW3

!
!===========================================
!

  subroutine poisson_load_fluidNEWFAST(nelmt,nnode,ibool,storederiv,storerhojw, &
                                       time,deltat,two_omega_earth,A_array_rot,B_array_rot,dispf,load)

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE
  use constants_solver, only: ROTATION_VAL

  use siem_math_library, only: dotmat
  use siem_gll_library

  implicit none

  integer,intent(in) :: nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: storederiv(NDIM,NGLLCUBE,NGLLCUBE,nelmt),storerhojw(NGLLCUBE,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE,nelmt),intent(inout) :: A_array_rot,B_array_rot
  real(kind=CUSTOM_REAL),intent(in) :: dispf(1,nnode) !\chi
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode)

  integer,parameter :: ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE)
  real(kind=CUSTOM_REAL) :: echi(NGLLCUBE,1),edisp(ndim,NGLLCUBE), &
                            eload(NGLLCUBE),gradchi(ndim,1)
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rot,B_rot,ux_rot,uy_rot
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE) :: source_euler_A,source_euler_B

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  real(kind=CUSTOM_REAL),dimension(:,:,:),allocatable :: deriv

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(deriv(ndim,NGLLCUBE,NGLLCUBE))
  deriv(:,:,:) = 0.0_CUSTOM_REAL

  load = zero
  do i_elmt = 1,nelmt
    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    echi(:,1) = dispf(1,egdof)

    ! compute diaplacement
    do i = 1,NGLLCUBE
      deriv(:,:,i) = storederiv(:,:,i,i_elmt)
      gradchi = matmul(deriv(:,:,i),echi)

      edisp(:,i) = gradchi(:,1)
      ! compute contribution of rotation and add to gradient of potential
      ! this term has no Z component
      if (ROTATION_VAL) then
        ! store the source for the Euler scheme for A_rotation and B_rotation
        two_omega_deltat = deltat*two_omega_earth

        cos_two_omega_t = cos(two_omega_earth*time)
        sin_two_omega_t = sin(two_omega_earth*time)

        ! time step deltat of Euler scheme is included in the source
        source_euler_A(i) = two_omega_deltat*(cos_two_omega_t*gradchi(2,1)+         &
                              sin_two_omega_t*gradchi(1,1))
        source_euler_B(i) = two_omega_deltat*(sin_two_omega_t*gradchi(2,1)-         &
                              cos_two_omega_t*gradchi(1,1))

        A_rot = A_array_rot(i,i_elmt)
        B_rot = B_array_rot(i,i_elmt)

        ux_rot =  A_rot*cos_two_omega_t+B_rot*sin_two_omega_t
        uy_rot = -A_rot*sin_two_omega_t+B_rot*cos_two_omega_t

        edisp(1,i) = edisp(1,i) + ux_rot
        edisp(2,i) = edisp(2,i) + uy_rot
      endif
    enddo

    ! integration
    eload = zero
    do i = 1,NGLLCUBE
      eload = eload + storerhojw(i,i_elmt)*matmul(transpose(deriv(:,:,i)),edisp(:,i))
    enddo
    load(egdof) = load(egdof) + eload

    ! update rotation term with Euler scheme
    if (ROTATION_VAL) then
      ! use the source saved above
      A_array_rot(:,i_elmt) = A_array_rot(:,i_elmt) + source_euler_A
      B_array_rot(:,i_elmt) = B_array_rot(:,i_elmt) + source_euler_B
    endif

  enddo !i_elmt = 1,nelmt

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_fluidNEWFAST

!
!===========================================
!

  subroutine poisson_load_fluidNEW3FAST(nelmt,nnode,ibool,storederiv,storerhojw, &
                                        time,deltat,two_omega_earth,A_array_rot,B_array_rot,dispf,nnode1,ibool1,load)

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF
  use constants_solver, only: ROTATION_VAL

  use siem_math_library, only: dotmat
  use siem_gll_library

  implicit none

  integer,intent(in) :: nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: storederiv(NDIM,NGLLCUBE_INF,NGLLCUBE_INF,nelmt),storerhojw(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE_INF,nelmt),intent(inout) :: A_array_rot,B_array_rot
  real(kind=CUSTOM_REAL),intent(in) :: dispf(1,nnode) !\chi
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF),egdof1(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: deriv(ndim,NGLLCUBE_INF,NGLLCUBE_INF),echi(NGLLCUBE_INF,1),edisp(ndim,NGLLCUBE_INF), &
                            eload(NGLLCUBE_INF),gradchi(ndim,1)
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rot,B_rot,ux_rot,uy_rot
  real(kind=CUSTOM_REAL),dimension(NGLLCUBE_INF) :: source_euler_A,source_euler_B

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  !! compute gauss-lobatto-legendre quadrature information
  !call gll_quadrature(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points,gll_weights, &
  !                    lagrange_gll,dlagrange_gll)

  load = zero
  do i_elmt = 1,nelmt
    egdof1 = reshape(ibool(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))
    egdof = ibool1(:,i_elmt)
    echi(:,1) = dispf(1,egdof1)

    ! compute diaplacement
    do i = 1,NGLLCUBE_INF
      deriv(:,:,i) = storederiv(:,:,i,i_elmt)
      gradchi = matmul(deriv(:,:,i),echi)

      edisp(:,i)=gradchi(:,1)
      ! compute contribution of rotation and add to gradient of potential
      ! this term has no Z component
      if (ROTATION_VAL) then
        ! store the source for the Euler scheme for A_rotation and B_rotation
        two_omega_deltat = deltat*two_omega_earth

        cos_two_omega_t = cos(two_omega_earth*time)
        sin_two_omega_t = sin(two_omega_earth*time)

        ! time step deltat of Euler scheme is included in the source
        source_euler_A(i) = two_omega_deltat*(cos_two_omega_t*gradchi(2,1)+         &
                              sin_two_omega_t*gradchi(1,1))
        source_euler_B(i) = two_omega_deltat*(sin_two_omega_t*gradchi(2,1)-         &
                              cos_two_omega_t*gradchi(1,1))

        A_rot = A_array_rot(i,i_elmt)
        B_rot = B_array_rot(i,i_elmt)

        ux_rot =  A_rot*cos_two_omega_t+B_rot*sin_two_omega_t
        uy_rot = -A_rot*sin_two_omega_t+B_rot*cos_two_omega_t

        edisp(1,i) = edisp(1,i) + ux_rot
        edisp(2,i) = edisp(2,i) + uy_rot
      endif
    enddo

    ! integration
    eload = zero
    do i = 1,NGLLCUBE_INF
      eload = eload+storerhojw(i,i_elmt)*matmul(transpose(deriv(:,:,i)),edisp(:,i))
    enddo
    load(egdof) = load(egdof) + eload

    ! update rotation term with Euler scheme
    if (ROTATION_VAL) then
      ! use the source saved above
      A_array_rot(:,i_elmt) = A_array_rot(:,i_elmt) + source_euler_A
      B_array_rot(:,i_elmt) = B_array_rot(:,i_elmt) + source_euler_B
    endif

  enddo !i_elmt = 1,nelmt

  ! multiply by 4*PI*G! or scaled
  load = -4.0_CUSTOM_REAL*load

  end subroutine poisson_load_fluidNEW3FAST

!
!===========================================
!

  subroutine poisson_load_onlyrho(iregion,nelmt,nnode, &
                                  ibool,xstore,ystore,zstore,rhostore,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(in) :: rhostore(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac,eld(NGLLCUBE),rho(NGLLCUBE)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ)

  real(kind=kdble) :: gll_weights(NGLLCUBE),gll_points(ndim,NGLLCUBE),dshape_hex8(ndim,ngnod,NGLLCUBE)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),jac(ndim,ndim)

  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gll
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gll

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE), &
           dlagrange_gll(ndim,NGLLCUBE,NGLLCUBE))
  lagrange_gll(:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  call zwgljd(xigll,wxgll,ngllx,jalpha,jbeta)
  call zwgljd(etagll,wygll,nglly,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,ngllz,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-1

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    rho = reshape(rhostore(:,:,:,i_elmt),(/NGLLCUBE/))

    eld = zero
    do i = 1,NGLLCUBE
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      eld = eld + real(rho(i)*lagrange_gll(i,:)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eld
  enddo

  ! multiply by 4*PI*G! or scaled
  load = 4.0_CUSTOM_REAL*load

  end subroutine poisson_load_onlyrho

!
!===========================================
!

  subroutine poisson_load_onlyrho3(iregion,nelmt,nnode, &
                                   ibool,xstore,ystore,zstore,rhostore,nnode1,ibool1,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode,nnode1
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt),ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(in) :: rhostore(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac,eld(NGLLCUBE_INF),rho(NGLLCUBE_INF)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX_INF),wxgll(NGLLX_INF),etagll(NGLLY_INF),wygll(NGLLY_INF), &
                      zetagll(NGLLZ_INF),wzgll(NGLLZ_INF)

  real(kind=kdble) :: gll_weights(NGLLCUBE_INF),gll_points(ndim,NGLLCUBE_INF),lagrange_gll(NGLLCUBE_INF,NGLLCUBE_INF), &
                      dlagrange_gll(ndim,NGLLCUBE_INF,NGLLCUBE_INF),dshape_hex8(ndim,ngnod,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),jac(ndim,ndim)

  call zwgljd(xigll,wxgll,NGLLX_INF,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY_INF,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ_INF,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !dnx=NGLLX-1; dny=NGLLY-1; dnz=NGLLZ-1

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);             ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);     ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = ibool1(:,i_elmt) !reshape(ibool1(:,:,:,i_elmt),(/NGLLCUBE_INF/))
    rho = reshape(rhostore(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/)) !WARNING: ONLY for 5 to 3 GLL extraction

    eld = zero
    do i = 1,NGLLCUBE_INF
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      eld = eld + real(rho(i)*lagrange_gll(i,:)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eld
  enddo

  ! multiply by 4*PI*G! or scaled
  load = 4.0_CUSTOM_REAL*load

  end subroutine poisson_load_onlyrho3

!
!===========================================
!

  subroutine poisson_load_onlyrhoFAST(iregion,nelmt,nnode, &
                                      ibool,storerhojw,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) :: storerhojw(NGLLCUBE,nelmt)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE)
  real(kind=CUSTOM_REAL) :: eld(NGLLCUBE)

  real(kind=kdble) :: gll_weights(NGLLCUBE),gll_points(ndim,NGLLCUBE)

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gll
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gll

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE), &
           dlagrange_gll(ndim,NGLLCUBE,NGLLCUBE))
  lagrange_gll(:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up
  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
    endif

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    eld = zero
    do i = 1,NGLLCUBE
      eld = eld + real(storerhojw(i,i_elmt)*lagrange_gll(i,:),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eld
  enddo

  ! multiply by 4*PI*G! or scaled
  load = 4.0_CUSTOM_REAL*load

  end subroutine poisson_load_onlyrhoFAST

!
!===========================================
!

  subroutine poisson_load_onlyrhoFAST3(iregion,nelmt, &
                                       storerhojw,nnode1,ibool1,load)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE, &
    NGLLX,NGLLY,NGLLZ,NGLLCUBE, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_innercore, only: idoubling_inner_core

  use siem_gll_library
  use siem_math_library, only: determinant,dotmat,invert

  implicit none

  integer,intent(in) :: iregion,nelmt,nnode1
  integer,intent(in) :: ibool1(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  storerhojw(NGLLCUBE_INF,nelmt)
  real(kind=CUSTOM_REAL),intent(out) :: load(nnode1)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: egdof(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: eld(NGLLCUBE_INF)

  real(kind=kdble) :: gll_weights(NGLLCUBE_INF),gll_points(ndim,NGLLCUBE_INF),lagrange_gll(NGLLCUBE_INF,NGLLCUBE_INF), &
                      dlagrange_gll(ndim,NGLLCUBE_INF,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL),parameter :: zero=0.0_CUSTOM_REAL

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  !TODO: can store deriv, and detjac*gll_weights(i) for speeding up

  load = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    egdof = ibool1(:,i_elmt) !reshape(ibool1(:,:,:,i_elmt),(/NGLLCUBE_INF/))
    eld = zero
    do i = 1,NGLLCUBE_INF
      eld = eld + real(storerhojw(i,i_elmt)*lagrange_gll(i,:),kind=CUSTOM_REAL)
    enddo
    load(egdof) = load(egdof) + eld
  enddo

  ! multiply by 4*PI*G! or scaled
  load = 4.0_CUSTOM_REAL*load

  end subroutine poisson_load_onlyrhoFAST3

!
!===========================================
!

  subroutine poisson_gravity(iregion,nelmt,nnode,ibool,xstore,ystore,zstore,pgrav, &
                             gradphi)

  use constants, only: IFLAG_IN_FICTITIOUS_CUBE,IREGION_INNER_CORE,NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_innercore, only: idoubling_inner_core
  !use specfem_par_crustmantle, only: rmassz_crust_mantle !TODO: remove this

  use siem_gll_library
  use siem_math_library, only: determinant,invert

  implicit none
  integer,intent(in) :: iregion,nelmt,nnode
  integer,intent(in) :: ibool(NGLLX,NGLLY,NGLLZ,nelmt)
  real(kind=CUSTOM_REAL),intent(in) ::  xstore(nnode),ystore(nnode),zstore(nnode)
  real(kind=CUSTOM_REAL),intent(in) :: pgrav(nnode)
  real(kind=CUSTOM_REAL),intent(out) :: gradphi(3,NGLLCUBE,nelmt)

  integer,parameter :: ndim = 3,ngnod = 8

  integer :: i,i_elmt
  integer :: dnx,dny,dnz,egdof(NGLLCUBE),ignod(ngnod)
  real(kind=CUSTOM_REAL) :: detjac

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble
  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ)
  real(kind=kdble) :: dshape_hex8(ndim,ngnod,NGLLCUBE),gll_weights(NGLLCUBE), &
                      gll_points(ndim,NGLLCUBE)

  real(kind=CUSTOM_REAL) :: coord(ngnod,ndim),deriv(ndim,NGLLCUBE),jac(ndim,ndim)

  real(kind=kdble),dimension(:,:), allocatable :: lagrange_gll
  real(kind=kdble),dimension(:,:,:), allocatable :: dlagrange_gll

  ! allocates arrays to avoid error about exceeding stack size limit
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE), &
           dlagrange_gll(ndim,NGLLCUBE,NGLLCUBE))
  lagrange_gll(:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  call zwgljd(xigll,wxgll,NGLLX,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(ndim,ngnod,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(ndim,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  dnx = NGLLX-1; dny = NGLLY-1; dnz = NGLLZ-1
  gradphi = zero
  do i_elmt = 1,nelmt
    ! suppress fictitious elements in central cube
    if (iregion == IREGION_INNER_CORE) then
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
    endif

    !ignod=reshape(ibool(1:NGLLX:dnx,1:NGLLY:dny,1:NGLLZ:dnz,i_elmt),(/ngnod/)) ! this is wrong!!!!
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool(1,1,1,i_elmt);               ignod(2) = ibool(NGLLX,1,1,i_elmt)
    ignod(3) = ibool(NGLLX,NGLLY,1,i_elmt);       ignod(4) = ibool(1,NGLLY,1,i_elmt)
    ! top corner nodes
    ignod(5) = ibool(1,1,NGLLZ,i_elmt);         ignod(6) = ibool(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool(NGLLX,NGLLY,NGLLZ,i_elmt); ignod(8) = ibool(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore(ignod)
    coord(:,2) = ystore(ignod)
    coord(:,3) = zstore(ignod)

    egdof = reshape(ibool(:,:,:,i_elmt),(/NGLLCUBE/))
    do i = 1,NGLLCUBE
      jac = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL) !jac = matmul(der,coord)
      detjac = determinant(jac)
      call invert(jac)
      deriv = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL) ! use der for gll
      gradphi(:,i,i_elmt) = matmul(deriv,pgrav(egdof))
    enddo
  enddo

  end subroutine poisson_gravity


end module siem_poisson

