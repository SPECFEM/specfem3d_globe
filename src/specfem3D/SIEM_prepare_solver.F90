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

  subroutine SIEM_prepare_solver()

  use specfem_par
  implicit none

  ! check if anything to do
  if (.not. FULL_GRAVITY) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing full gravity solver"
    call flush_IMAIN()
  endif

  ! safety stop
  stop 'FULL_GRAVITY not fully implemented yet'

! TODO: full gravity is not working yet, needs to fully implement solver...
#ifdef USE_PETSC_NOT_WORKING_YET

  ! compute and store integration coefficients
  call SIEM_prepare_solver_preintegrate3()
  call SIEM_prepare_solver_preintegrate()

  ! sets the stiffness matrices for Poisson's solver
  ! calculate dprecon
  ! allocates load, regional pgravs (e.g. pgrav_cm1)
  call SIEM_prepare_solver_poisson()

  ! create sparse matrix
  if (SOLVER == PETSC) call SIEM_prepare_solver_sparse()

#endif

  end subroutine SIEM_prepare_solver

!
!-------------------------------------------------------------------------------
!

! TODO: full gravity is not working yet, needs to fully implement solver...
#ifdef USE_PETSC_NOT_WORKING_YET

  subroutine prepare_solver_preintegrate()

  use specfem_par, only: CUSTOM_REAL,myrank,NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE

  use specfem_par_crustmantle, only: ibool_crust_mantle,NSPEC_CRUST_MANTLE, &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
    rhostore_crust_mantle

  use specfem_par_innercore, only: idoubling_inner_core,IFLAG_IN_FICTITIOUS_CUBE, &
    ibool_inner_core,NSPEC_INNER_CORE,xstore_inner_core,ystore_inner_core, &
    zstore_inner_core,rhostore_inner_core

  use specfem_par_outercore, only: ibool_outer_core,NSPEC_OUTER_CORE, &
    xstore_outer_core,ystore_outer_core,zstore_outer_core,rhostore_outer_core

  !use specfem_par_infinite
  use constants_solver, only: IFLAG_IN_FICTITIOUS_CUBE

  use gll_library1, only: kdble,zwgljd,dshape_function_hex8,gll_quadrature
  use math_library, only: determinant,invert

  use specfem_par_innercore, only: idoubling_inner_core

  use specfem_par_full_gravity, only: lagrange_gll, &
    storederiv_cm,storerhojw_cm, storedetjac_cm, SAVE_JACOBIAN_ENSIGHT, &
    storederiv_ic,storerhojw_ic, &
    storederiv_oc,storerhojw_oc

  implicit none

  integer,parameter :: ngnod = 8

  integer :: i,i_elmt, j, k
  integer :: ignod(ngnod) !dnx,dny,dnz,
  real(kind=CUSTOM_REAL) :: detjac,rho(NGLLCUBE)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble

  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ), detjac_cm_tmp(NGLLCUBE), element_detjac(NGLLX,NGLLY,NGLLZ)

  real(kind=kdble) :: gll_weights(NGLLCUBE),gll_points(NDIM,NGLLCUBE), &
                      dlagrange_gll(NDIM,NGLLCUBE,NGLLCUBE),dshape_hex8(NDIM,ngnod,NGLLCUBE)

  real(kind=CUSTOM_REAL) :: coord(ngnod,NDIM),deriv(NDIM,NGLLCUBE),jac(NDIM,NDIM)

  call zwgljd(xigll,wxgll,NGLLX,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(NDIM,ngnod,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll, &
                            zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights, &
                      lagrange_gll,dlagrange_gll)

  ! inner core
  storederiv_ic = 0.0_CUSTOM_REAL
  storerhojw_ic = 0.0_CUSTOM_REAL
  do i_elmt = 1,NSPEC_INNER_CORE
    ! suppress fictitious elements in central cube
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool_inner_core(1,1,1,i_elmt)
    ignod(2) = ibool_inner_core(NGLLX,1,1,i_elmt)
    ignod(3) = ibool_inner_core(NGLLX,NGLLY,1,i_elmt)
    ignod(4) = ibool_inner_core(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool_inner_core(1,1,NGLLZ,i_elmt)
    ignod(6) = ibool_inner_core(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool_inner_core(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8) = ibool_inner_core(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1)=xstore_inner_core(ignod)
    coord(:,2)=ystore_inner_core(ignod)
    coord(:,3)=zstore_inner_core(ignod)
    rho = reshape(rhostore_inner_core(:,:,:,i_elmt),(/NGLLCUBE/))

    do i = 1,NGLLCUBE
      jac = matmul(dshape_hex8(:,:,i),coord)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_ic(:,:,i,i_elmt) = matmul(jac,dlagrange_gll(:,i,:))
      storerhojw_ic(i,i_elmt) = rho(i)*detjac*gll_weights(i)
      enddo
  enddo

  ! outer core
  storederiv_oc = 0.0_CUSTOM_REAL
  storerhojw_oc = 0.0_CUSTOM_REAL
  do i_elmt = 1,NSPEC_OUTER_CORE
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool_outer_core(1,1,1,i_elmt)
    ignod(2) = ibool_outer_core(NGLLX,1,1,i_elmt)
    ignod(3) = ibool_outer_core(NGLLX,NGLLY,1,i_elmt)
    ignod(4) = ibool_outer_core(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool_outer_core(1,1,NGLLZ,i_elmt)
    ignod(6) = ibool_outer_core(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool_outer_core(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8) = ibool_outer_core(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore_outer_core(ignod)
    coord(:,2) = ystore_outer_core(ignod)
    coord(:,3) = zstore_outer_core(ignod)
    rho = reshape(rhostore_outer_core(:,:,:,i_elmt),(/NGLLCUBE/))

    do i = 1,NGLLCUBE
      jac = matmul(dshape_hex8(:,:,i),coord)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_oc(:,:,i,i_elmt) = matmul(jac,dlagrange_gll(:,i,:))
      storerhojw_oc(i,i_elmt) = rho(i)*detjac*gll_weights(i)
    enddo
  enddo

  ! crust mantle
  storederiv_cm = 0.0_CUSTOM_REAL
  storerhojw_cm = 0.0_CUSTOM_REAL
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1)=ibool_crust_mantle(1,1,1,i_elmt)
    ignod(2)=ibool_crust_mantle(NGLLX,1,1,i_elmt)
    ignod(3)=ibool_crust_mantle(NGLLX,NGLLY,1,i_elmt)
    ignod(4)=ibool_crust_mantle(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5)=ibool_crust_mantle(1,1,NGLLZ,i_elmt)
    ignod(6)=ibool_crust_mantle(NGLLX,1,NGLLZ,i_elmt)
    ignod(7)=ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8)=ibool_crust_mantle(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1)=xstore_crust_mantle(ignod)
    coord(:,2)=ystore_crust_mantle(ignod)
    coord(:,3)=zstore_crust_mantle(ignod)
    rho = reshape(rhostore_crust_mantle(:,:,:,i_elmt),(/NGLLCUBE/))

    do i = 1,NGLLCUBE
      jac = matmul(dshape_hex8(:,:,i),coord)
      detjac=determinant(jac)
      call invert(jac)
      storederiv_cm(:,:,i,i_elmt)=matmul(jac,dlagrange_gll(:,i,:))
      storerhojw_cm(i,i_elmt)=rho(i)*detjac*gll_weights(i)
      detjac_cm_tmp(i) = detjac
    enddo

    if (SAVE_JACOBIAN_ENSIGHT) then
      element_detjac(:,:,:) = reshape(detjac_cm_tmp(:), (/NGLLX,NGLLY,NGLLZ/) )
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            storedetjac_cm(ibool_crust_mantle(i,j,k,i_elmt)) = element_detjac(i,j,k)
          enddo
        enddo
      enddo
    endif

  enddo

  end subroutine prepare_solver_preintegrate

!
!-------------------------------------------------------------------------------
!

  subroutine prepare_solver_preintegrate3()

  use specfem_par, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF

  use specfem_par_crustmantle, only: ibool_crust_mantle,NSPEC_CRUST_MANTLE, &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
    rhostore_crust_mantle,storederiv_cm1,storerhojw_cm1,storejw_cm1

  use specfem_par_innercore, only: idoubling_inner_core,IFLAG_IN_FICTITIOUS_CUBE, &
    ibool_inner_core,NSPEC_INNER_CORE,xstore_inner_core,ystore_inner_core, &
    zstore_inner_core,rhostore_inner_core,storederiv_ic1,storerhojw_ic1

  use specfem_par_outercore, only: ibool_outer_core,NSPEC_OUTER_CORE, &
    xstore_outer_core,ystore_outer_core,zstore_outer_core,rhostore_outer_core, &
    storederiv_oc1,storerhojw_oc1

  !use specfem_par_infinite
  use constants_solver, only: IFLAG_IN_FICTITIOUS_CUBE

  use gll_library1, only: kdble,zwgljd,dshape_function_hex8,gll_quadrature
  use math_library, only: determinant,invert

  use specfem_par_innercore, only: idoubling_inner_core

  use specfem_par_full_gravity, only: lagrange_gll1

  implicit none

  integer,parameter :: ngnod = 8

  integer :: i,i_elmt
  integer :: ignod(ngnod) !dnx,dny,dnz,
  real(kind=CUSTOM_REAL) :: detjac,rho(NGLLCUBE_INF)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble

  real(kind=kdble) :: xigll(NGLLX_INF),wxgll(NGLLX_INF),etagll(NGLLY_INF),wygll(NGLLY_INF), &
                      zetagll(NGLLZ_INF),wzgll(NGLLZ_INF)

  real(kind=kdble) :: gll_weights1(NGLLCUBE_INF),gll_points1(NDIM,NGLLCUBE_INF), &
                      dlagrange_gll1(NDIM,NGLLCUBE_INF,NGLLCUBE_INF),dshape_hex8(NDIM,ngnod,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL) :: coord(ngnod,NDIM),deriv(NDIM,NGLLCUBE_INF),jac(NDIM,NDIM)

  call zwgljd(xigll,wxgll,NGLLX_INF,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY_INF,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ_INF,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(NDIM,ngnod,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,xigll,etagll, &
  zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(NDIM,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points1,gll_weights1, &
  lagrange_gll1,dlagrange_gll1)

  ! inner core
  storederiv_ic1 = 0.0_CUSTOM_REAL
  storerhojw_ic1 = 0.0_CUSTOM_REAL
  do i_elmt = 1,NSPEC_INNER_CORE
    ! suppress fictitious elements in central cube
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle

    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1)=ibool_inner_core(1,1,1,i_elmt)
    ignod(2)=ibool_inner_core(NGLLX,1,1,i_elmt)
    ignod(3)=ibool_inner_core(NGLLX,NGLLY,1,i_elmt)
    ignod(4)=ibool_inner_core(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5)=ibool_inner_core(1,1,NGLLZ,i_elmt)
    ignod(6)=ibool_inner_core(NGLLX,1,NGLLZ,i_elmt)
    ignod(7)=ibool_inner_core(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8)=ibool_inner_core(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1)=xstore_inner_core(ignod)
    coord(:,2)=ystore_inner_core(ignod)
    coord(:,3)=zstore_inner_core(ignod)
    rho = reshape(rhostore_inner_core(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))

    do i = 1,NGLLCUBE_INF
      jac = matmul(dshape_hex8(:,:,i),coord)
      detjac=determinant(jac)
      call invert(jac)
      storederiv_ic1(:,:,i,i_elmt)=matmul(jac,dlagrange_gll1(:,i,:))
      storerhojw_ic1(i,i_elmt)=rho(i)*detjac*gll_weights1(i)
    enddo
  enddo

  ! outer core
  storederiv_oc1 = 0.0_CUSTOM_REAL
  storerhojw_oc1 = 0.0_CUSTOM_REAL
  do i_elmt = 1,NSPEC_OUTER_CORE
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1)=ibool_outer_core(1,1,1,i_elmt)
    ignod(2)=ibool_outer_core(NGLLX,1,1,i_elmt)
    ignod(3)=ibool_outer_core(NGLLX,NGLLY,1,i_elmt)
    ignod(4)=ibool_outer_core(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5)=ibool_outer_core(1,1,NGLLZ,i_elmt)
    ignod(6)=ibool_outer_core(NGLLX,1,NGLLZ,i_elmt)
    ignod(7)=ibool_outer_core(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8)=ibool_outer_core(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1)=xstore_outer_core(ignod)
    coord(:,2)=ystore_outer_core(ignod)
    coord(:,3)=zstore_outer_core(ignod)
    rho = reshape(rhostore_outer_core(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt), (/NGLLCUBE_INF/))

    do i = 1,NGLLCUBE_INF
      jac = matmul(dshape_hex8(:,:,i),coord)
      detjac=determinant(jac)
      call invert(jac)
      storederiv_oc1(:,:,i,i_elmt)=matmul(jac,dlagrange_gll1(:,i,:))
      storerhojw_oc1(i,i_elmt)=rho(i)*detjac*gll_weights1(i)
    enddo
  enddo

  ! crust mantle
  storederiv_cm1 = 0.0_CUSTOM_REAL
  storerhojw_cm1 = 0.0_CUSTOM_REAL
  storejw_cm1 = 0.0_CUSTOM_REAL
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1)=ibool_crust_mantle(1,1,1,i_elmt)
    ignod(2)=ibool_crust_mantle(NGLLX,1,1,i_elmt)
    ignod(3)=ibool_crust_mantle(NGLLX,NGLLY,1,i_elmt)
    ignod(4)=ibool_crust_mantle(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5)=ibool_crust_mantle(1,1,NGLLZ,i_elmt)
    ignod(6)=ibool_crust_mantle(NGLLX,1,NGLLZ,i_elmt)
    ignod(7)=ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8)=ibool_crust_mantle(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1)=xstore_crust_mantle(ignod)
    coord(:,2)=ystore_crust_mantle(ignod)
    coord(:,3)=zstore_crust_mantle(ignod)
    rho = reshape(rhostore_crust_mantle(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt), (/NGLLCUBE_INF/))

    do i = 1,NGLLCUBE_INF
      jac = matmul(dshape_hex8(:,:,i),coord)
      detjac=determinant(jac)
      call invert(jac)
      storederiv_cm1(:,:,i,i_elmt)=matmul(jac,dlagrange_gll1(:,i,:))
      storerhojw_cm1(i,i_elmt)=rho(i)*detjac*gll_weights1(i)
      storejw_cm1(i,i_elmt)=detjac*gll_weights1(i)
    enddo
  enddo

  end subroutine prepare_solver_preintegrate3

!
!-------------------------------------------------------------------------------
!

! TODO: check this why is it so slow poisson_stiffness is very slow due to loop
! of jacobian and etc. If we store jacobian before hand it should be faster!

  subroutine prepare_solver_poisson()

  use math_library_mpi, only: maxscal,minscal
  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_trinfinite
  use specfem_par_infinite
  use index_region
  use poisson, only: poisson_stiffness,poisson_stiffnessINF,poisson_stiffness3, &
  poisson_stiffnessINF3

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  allocating poisson level-1 solver arrays"
    call flush_IMAIN()
  endif

  ! indexify regions
  call get_index_region()

  ! allocate inode arrays
  allocate(inode_elmt_cm(NGLLCUBE,NSPEC_CRUST_MANTLE))
  allocate(inode_elmt_cm1(NGLLCUBE_INF,NSPEC_CRUST_MANTLE))
  allocate(inode_elmt_ic(NGLLCUBE,NSPEC_INNER_CORE))
  allocate(inode_elmt_ic1(NGLLCUBE_INF,NSPEC_INNER_CORE))
  allocate(inode_elmt_oc(NGLLCUBE,NSPEC_OUTER_CORE))
  allocate(inode_elmt_oc1(NGLLCUBE_INF,NSPEC_OUTER_CORE))
  ! trinfinite arrays
  if (ADD_TRINF) then
    allocate(inode_elmt_trinf(NGLLCUBE,NSPEC_TRINFINITE))
    allocate(inode_elmt_trinf1(NGLLCUBE_INF,NSPEC_TRINFINITE))
  endif
  ! infinite arrays
  allocate(inode_elmt_inf(NGLLCUBE,NSPEC_INFINITE))
  allocate(inode_elmt_inf1(NGLLCUBE_INF,NSPEC_INFINITE))

  ! Level-1 solver-------------------
  allocate(storekmat_crust_mantle1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_CRUST_MANTLE), &
           dprecon_crust_mantle1(nnode_cm1))
  allocate(storekmat_outer_core1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_OUTER_CORE), &
           dprecon_outer_core1(nnode_oc1))
  allocate(storekmat_inner_core1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_INNER_CORE), &
           dprecon_inner_core1(nnode_ic1))
  if (ADD_TRINF) then
    allocate(storekmat_trinfinite1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_TRINFINITE), &
             dprecon_trinfinite1(nnode_trinf1))
  endif
  allocate(storekmat_infinite1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_INFINITE), &
           dprecon_infinite1(nnode_inf1))

  allocate(dprecon1(0:neq1),load1(0:neq1),pgrav_ic1(nnode_ic1), &
           pgrav_oc1(nnode_oc1),pgrav_cm1(nnode_cm1),pgrav_trinf1(nnode_trinf1), &
           pgrav_inf1(nnode_inf1))

  if (SIMULATION_TYPE == 3) then
    allocate(b_load1(0:neq1), b_pgrav_ic1(nnode_ic1), b_pgrav_oc1(nnode_oc1), &
             b_pgrav_cm1(nnode_cm1), b_pgrav_trinf1(nnode_trinf1), b_pgrav_inf1(nnode_inf1))
  endif

  ! crust mantle
  call poisson_stiffness3(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                          NGLOB_CRUST_MANTLE,ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle, &
                          zstore_crust_mantle,nnode_cm1,inode_elmt_cm1,storekmat_crust_mantle1, &
                          dprecon_crust_mantle1)

  ! outer core
  call poisson_stiffness3(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                          ibool_outer_core,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                          nnode_oc1,inode_elmt_oc1,storekmat_outer_core1,dprecon_outer_core1)

  ! inner core
  call poisson_stiffness3(IREGION_INNER_CORE,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                          ibool_inner_core,xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                          nnode_ic1,inode_elmt_ic1,storekmat_inner_core1,dprecon_inner_core1)

  if (ADD_TRINF) then
    ! transition infinite
    call poisson_stiffness3(IREGION_TRINFINITE,NSPEC_TRINFINITE,NGLOB_TRINFINITE, &
                            ibool_trinfinite,xstore_trinfinite,ystore_trinfinite,zstore_trinfinite, &
                            nnode_trinf1,inode_elmt_trinf1,storekmat_trinfinite1,dprecon_trinfinite1)
  endif

  ! infinite layer
  call poisson_stiffnessINF3(NSPEC_INFINITE,NGLOB_INFINITE,ibool_infinite, &
                             xstore_infinite,ystore_infinite,zstore_infinite,nnode_inf1,inode_elmt_inf1, &
                             storekmat_infinite1,dprecon_infinite1)

  call sync_all()

  ! assemble stiffness matrices
  ! assemble across the MPI processes in a region
  ! crust_mantle
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_cm1,dprecon_crust_mantle1, &
                           num_interfaces_crust_mantle1,max_nibool_interfaces_crust_mantle1, &
                           nibool_interfaces_crust_mantle1,ibool_interfaces_crust_mantle1, &
                           my_neighbors_crust_mantle1)

  ! outer core
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_oc1,dprecon_outer_core1, &
                           num_interfaces_outer_core1,max_nibool_interfaces_outer_core1, &
                           nibool_interfaces_outer_core1,ibool_interfaces_outer_core1, &
                           my_neighbors_outer_core1)

  ! inner core
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_ic1,dprecon_inner_core1, &
                           num_interfaces_inner_core1,max_nibool_interfaces_inner_core1, &
                           nibool_interfaces_inner_core1,ibool_interfaces_inner_core1, &
                           my_neighbors_inner_core1)

  ! transition infinite
  if (ADD_TRINF) then
    call assemble_MPI_scalar(NPROCTOT_VAL,nnode_trinf1,dprecon_trinfinite1, &
                             num_interfaces_trinfinite1,max_nibool_interfaces_trinfinite1, &
                             nibool_interfaces_trinfinite1,ibool_interfaces_trinfinite1, &
                             my_neighbors_trinfinite1)
  endif

  ! infinite
  call assemble_MPI_scalar(NPROCTOT_VAL,nnode_inf1,dprecon_infinite1, &
                           num_interfaces_infinite1,max_nibool_interfaces_infinite1, &
                           nibool_interfaces_infinite1,ibool_interfaces_infinite1, &
                           my_neighbors_infinite1)

  call sync_all()

  ! assemble across the different regions in a process
  dprecon1 = zero
  ! crust_mantle
  dprecon1(gdof_cm1) = dprecon1(gdof_cm1)+dprecon_crust_mantle1

  ! outer core
  dprecon1(gdof_oc1) = dprecon1(gdof_oc1)+dprecon_outer_core1

  ! inner core
  dprecon1(gdof_ic1) = dprecon1(gdof_ic1)+dprecon_inner_core1

  ! transition infinite
  if (ADD_TRINF) then
    dprecon1(gdof_trinf1) = dprecon1(gdof_trinf1)+dprecon_trinfinite1
  endif

  ! infinite
  dprecon1(gdof_inf1) = dprecon1(gdof_inf1)+dprecon_infinite1

  dprecon1(0) = 0.0_CUSTOM_REAL

  call sync_all()

  ! invert preconditioner
  !dprecon1(1:)=1.0_CUSTOM_REAL/dprecon1(1:)
  !-----------------Level-1 solver


  ! Level-2 solver------------------
  if (SOLVER_5GLL) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  allocating poisson level-2 solver arrays"
      call flush_IMAIN()
    endif

    allocate(storekmat_crust_mantle(NGLLCUBE,NGLLCUBE,NSPEC_CRUST_MANTLE), &
             dprecon_crust_mantle(NGLOB_CRUST_MANTLE))
    allocate(storekmat_outer_core(NGLLCUBE,NGLLCUBE,NSPEC_OUTER_CORE), &
             dprecon_outer_core(NGLOB_OUTER_CORE))
    allocate(storekmat_inner_core(NGLLCUBE,NGLLCUBE,NSPEC_INNER_CORE), &
             dprecon_inner_core(NGLOB_INNER_CORE))
    if (ADD_TRINF) then
      allocate(storekmat_trinfinite(NGLLCUBE,NGLLCUBE,NSPEC_TRINFINITE), &
      dprecon_trinfinite(NGLOB_TRINFINITE))
    endif
    allocate(storekmat_infinite(NGLLCUBE,NGLLCUBE,NSPEC_INFINITE), &
             dprecon_infinite(NGLOB_INFINITE))
    allocate(dprecon(0:neq),load(0:neq))

    ! better to make dprecon_* local rather than global

    ! crust mantle
    call poisson_stiffness(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                           NGLOB_CRUST_MANTLE,ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle, &
                           zstore_crust_mantle,storekmat_crust_mantle,dprecon_crust_mantle)
    ! outer core
    call poisson_stiffness(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                           ibool_outer_core,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                           storekmat_outer_core,dprecon_outer_core)
    ! inner core
    call poisson_stiffness(IREGION_INNER_CORE,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                           ibool_inner_core,xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                           storekmat_inner_core,dprecon_inner_core)

    ! transition infinite
    if (ADD_TRINF) then
      call poisson_stiffness(IREGION_TRINFINITE,NSPEC_TRINFINITE,NGLOB_TRINFINITE, &
                             ibool_trinfinite,xstore_trinfinite,ystore_trinfinite,zstore_trinfinite, &
                             storekmat_trinfinite,dprecon_trinfinite)
    endif

    ! infinite layer
    call poisson_stiffnessINF(NSPEC_INFINITE,NGLOB_INFINITE, &
                              ibool_infinite,xstore_infinite,ystore_infinite,zstore_infinite, &
                              storekmat_infinite,dprecon_infinite)

    call sync_all()

    ! assemble stiffness matrices
    ! assemble across the MPI processes in a region
    ! crust_mantle
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE,dprecon_crust_mantle, &
                             num_interfaces_crust_mantle,max_nibool_interfaces_crust_mantle, &
                             nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                             my_neighbors_crust_mantle)

    ! outer core
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_OUTER_CORE,dprecon_outer_core, &
                             num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                             nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                             my_neighbors_outer_core)

    ! inner core
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INNER_CORE,dprecon_inner_core, &
                             num_interfaces_inner_core,max_nibool_interfaces_inner_core, &
                             nibool_interfaces_inner_core,ibool_interfaces_inner_core, &
                             my_neighbors_inner_core)

    ! transition infinite

    if (ADD_TRINF) then
      call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_TRINFINITE,dprecon_trinfinite, &
                               num_interfaces_trinfinite,max_nibool_interfaces_trinfinite, &
                               nibool_interfaces_trinfinite,ibool_interfaces_trinfinite, &
                               my_neighbors_trinfinite)
    endif

    ! infinite
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INFINITE,dprecon_infinite, &
                             num_interfaces_infinite,max_nibool_interfaces_infinite, &
                             nibool_interfaces_infinite,ibool_interfaces_infinite, &
                             my_neighbors_infinite)

    call sync_all()

    ! assemble across the different regions in a process
    dprecon = zero
    ! crust_mantle
    dprecon(gdof_cm) = dprecon(gdof_cm)+dprecon_crust_mantle

    ! outer core
    dprecon(gdof_oc) = dprecon(gdof_oc)+dprecon_outer_core

    ! inner core
    dprecon(gdof_ic) = dprecon(gdof_ic)+dprecon_inner_core

    ! transition infinite

    if (ADD_TRINF) then
      dprecon(gdof_trinf) = dprecon(gdof_trinf)+dprecon_trinfinite
    endif

    ! infinite
    dprecon(gdof_inf) = dprecon(gdof_inf)+dprecon_infinite

    dprecon(0) = 0.0_CUSTOM_REAL

    call sync_all()
    !--------------------Level-2 solver
  endif ! if (SOLVER_5GLL) then

  return

  end subroutine prepare_solver_poisson

!
!-------------------------------------------------------------------------------
!

! TODO: check this why is it so slow poisson_stiffness is very slow due to loop
! of jacobian and etc. If we store jacobian before hand it may be faster!
! DEVELOPER
!   Hom N Gharti
! HISTORY
!   Sep 30,2013

  subroutine prepare_solver_sparse()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_trinfinite
  use specfem_par_infinite
  use index_region
  use math_library, only: i_uniinv
  use math_library_mpi, only: maxscal,minvec,maxvec
  implicit none
  logical :: ismpi
  integer :: errcode
  character(len=250) :: errtag

  ! sparse stage 0
  integer :: i,j,i_elmt,i_count,n,ncount
  integer :: igdof,jgdof
  integer :: nmax !,nsparse
  integer :: nedof_ic,nedof_oc,nedof_cm,nedof_trinf,nedof_inf
  integer :: gdof_elmt(NEDOF),ggdof_elmt(NEDOF)
  integer :: nmax1 !,nsparse1
  integer :: nedof_ic1,nedof_oc1,nedof_cm1,nedof_trinf1,nedof_inf1
  integer :: gdof_elmt1(NEDOF1),ggdof_elmt1(NEDOF1)
  integer,allocatable :: imap_ic(:),imap_oc(:),imap_cm(:),imap_trinf(:), &
  imap_inf(:)
  integer,allocatable :: ind0(:),iorder(:),row0(:),col0(:),grow0(:),gcol0(:)
  real(kind=CUSTOM_REAL),allocatable :: kmat0(:)
  integer,allocatable :: ind1(:),row1(:),col1(:)
  real(kind=CUSTOM_REAL),allocatable :: kmat1(:)

  integer :: nglob_ic,nglob_oc,nglob_cm,nglob_trinf,nglob_inf
  integer :: nglob_ic1,nglob_oc1,nglob_cm1,nglob_trinf1,nglob_inf1

  character(len=12) :: spm
  character(len=60) :: fname

  integer :: nx,ny,nz,nbyte,off0,gmin,gmax
  integer :: i_bool,ibool,i_s,i0,i1,ier,j_proc
  logical :: isbig

  ! counting nonzero elements in offdiagonal portion
  integer :: grow,ig0,ig1,ind,neq_part1
  integer,allocatable :: nzero_rowoff1(:)
  integer,allocatable :: igorder(:)

  ismpi = .true.

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  preparing sparse matrix solver"
    call flush_IMAIN()
  endif

  !===============================================================================
  ! Level-1 solver
  !===============================================================================
  ! Number of DOFs per element in each region
  nedof_ic1 = NEDOFU1+NEDOFPHI1
  nedof_oc1 = NEDOFCHI1+NEDOFP1+NEDOFPHI1
  nedof_cm1 = NEDOFU1+NEDOFPHI1
  nedof_trinf1 = NEDOFPHI1
  nedof_inf1 = NEDOFPHI1

  ! Maximum DOF in array - number of elements * Element_dof^2
  nmax1 = NSPEC_INNER_CORE*(nedof_ic1*nedof_ic1)+                                &
          NSPEC_OUTER_CORE*(nedof_oc1*nedof_oc1)+                                 &
          NSPEC_CRUST_MANTLE*(nedof_cm1*nedof_cm1)+                               &
          NSPEC_TRINFINITE*(nedof_trinf1*nedof_trinf1)+                           &
          NSPEC_INFINITE*(nedof_inf1*nedof_inf1)

  allocate(col0(nmax1),row0(nmax1),gcol0(nmax1),grow0(nmax1),kmat0(nmax1))

  !debug
  if (myrank == 0) then
    print *,' -- Elemental DOFs for IC : ', nedof_ic1
    print *,' -- Maximum DOFs (nmax1)  : ', nmax1
  endif

  ! Allocate map for each region
  allocate(imap_ic(nedof_ic1),imap_oc(nedof_oc1),imap_cm(nedof_cm1), &
           imap_trinf(nedof_trinf1),imap_inf(nedof_inf1))

  ! I THINK THIS SYNTAX MEANS CREATE A RANGE?
  imap_ic = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapu1, imapphi1 /)
  imap_oc = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapchi1, imapp1, imapphi1 /)
  imap_cm = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapu1, imapphi1 /)
  imap_trinf = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapphi1 /)
  imap_inf = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapphi1 /)

  ! read global degrees of freedoms from DATABASE files
  write(spm,*)myrank
  fname='DATABASES_MPI/gdof1_proc'//trim(adjustl(spm))
  open(10,file=fname,action='read',status='old')
  ! inner core
  read(10,*)nglob_ic1                   ! Global DOF in inner core
  allocate(ggdof_ic1(NNDOF,nglob_ic1))
  read(10,*)ggdof_ic1
  ! outer core
  read(10,*)nglob_oc1                   ! Global DOF in outer core
  allocate(ggdof_oc1(NNDOF,nglob_oc1))
  read(10,*)ggdof_oc1
  ! crust mantle
  read(10,*)nglob_cm1                   ! Global DOF in crust mantle
  allocate(ggdof_cm1(NNDOF,nglob_cm1))
  read(10,*)ggdof_cm1
  ! transition
  read(10,*)nglob_trinf1                ! Global DOF in transition
  allocate(ggdof_trinf1(NNDOF,nglob_trinf1))
  read(10,*)ggdof_trinf1
  ! infinite elements
  read(10,*)nglob_inf1                  ! Global DOF in infinite
  allocate(ggdof_inf1(NNDOF,nglob_inf1))
  read(10,*)ggdof_inf1
  close(10)

  ! Find maximum ID (dof value) for any of the regions
  ngdof1 = maxscal(maxval( (/ maxval(ggdof_ic1),maxval(ggdof_oc1), &
           maxval(ggdof_cm1),maxval(ggdof_trinf1),maxval(ggdof_inf1) /) ))

  !debug
  if (myrank == 0) write(*,'(a,i12)')' -- Total global degrees of freedom1: ',ngdof1

  ! stage 0: store all elements
  ncount = 0
  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    ! Skip fictitious inner core cube
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle

    ! Note: gdof_ic1 defined in specfem_par_innercore
    ! Fetch gdof for IC element only on NGLLCUBE_INF points
    gdof_elmt1 = reshape(gdof_ic1(inode_elmt_ic1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1 = reshape(ggdof_ic1(:,inode_elmt_ic1(:,i_elmt)),(/NEDOF1/))
    !if (myrank==0) print *,'ICkmat zeros1:',count(storekmat_inner_core1(:,:,i_elmt)==0.0_CUSTOM_REAL)
    !if (myrank==0.and.i_elmt==1)print *,'kmat1:',storekmat_inner_core1(1,:,i_elmt)
    do i = 1,nedof_ic1
      do j = 1,nedof_ic1
        igdof = gdof_elmt1(imap_ic(i))
        jgdof = gdof_elmt1(imap_ic(j))
        ! If a degree of freedom and a non-zero kmat:
        if (igdof > 0.and.jgdof > 0.and.storekmat_inner_core1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          ncount = ncount+1
          ! Local (MPI?) map?
          row0(ncount) = igdof
          col0(ncount) = jgdof
          ! Global map?
          grow0(ncount) = ggdof_elmt1(imap_ic(i))
          gcol0(ncount) = ggdof_elmt1(imap_ic(j))
        endif
      enddo
    enddo
  enddo
  call sync_all()

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    gdof_elmt1 = reshape(gdof_oc1(inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1 = reshape(ggdof_oc1(:,inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
    do i = 1,nedof_oc1
      do j = 1,nedof_oc1
        igdof = gdof_elmt1(imap_oc(i))
        jgdof = gdof_elmt1(imap_oc(j))
        if (igdof > 0.and.jgdof > 0.and.storekmat_outer_core1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          !if (myrank==0) write(1111,*) igdof,jgdof,storekmat_outer_core1(i,j,i_elmt)
          ncount = ncount+1
          row0(ncount) = igdof
          col0(ncount) = jgdof
          grow0(ncount) = ggdof_elmt1(imap_oc(i))
          gcol0(ncount) = ggdof_elmt1(imap_oc(j))
        endif
      enddo
    enddo
  enddo


  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    gdof_elmt1 = reshape(gdof_cm1(inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1 = reshape(ggdof_cm1(:,inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
    !if (myrank==0) print *,'CMkmat zeros1:',count(storekmat_crust_mantle1(:,:,i_elmt)==0.0_CUSTOM_REAL)
    do i = 1,nedof_cm1
      do j = 1,nedof_cm1
        igdof = gdof_elmt1(imap_cm(i))
        jgdof = gdof_elmt1(imap_cm(j))
        if (igdof > 0.and.jgdof > 0.and.storekmat_crust_mantle1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          ncount = ncount+1
          row0(ncount) = igdof
          col0(ncount) = jgdof
          grow0(ncount) = ggdof_elmt1(imap_cm(i))
          gcol0(ncount) = ggdof_elmt1(imap_cm(j))
        endif
      enddo
    enddo
  enddo

  ! transition infinite
  do i_elmt = 1,NSPEC_TRINFINITE
    gdof_elmt1 = reshape(gdof_trinf1(inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1 = reshape(ggdof_trinf1(:,inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
    !if (myrank==0) print *,'TRINFkmat zeros1:',count(storekmat_trinfinite1(:,:,i_elmt)==0.0_CUSTOM_REAL)
    do i = 1,nedof_trinf1
      do j = 1,nedof_trinf1
        igdof = gdof_elmt1(imap_trinf(i))
        jgdof = gdof_elmt1(imap_trinf(j))
        if (igdof > 0.and.jgdof > 0.and.storekmat_trinfinite1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
          ncount = ncount+1
          row0(ncount) = igdof
          col0(ncount) = jgdof
          grow0(ncount) = ggdof_elmt1(imap_trinf(i))
          gcol0(ncount) = ggdof_elmt1(imap_trinf(j))
        endif
      enddo
    enddo
  enddo

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    gdof_elmt1 = reshape(gdof_inf1(inode_elmt_inf1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1 = reshape(ggdof_inf1(:,inode_elmt_inf1(:,i_elmt)),(/NEDOF1/))
    !if (myrank==0) print *,'INFkmat zeros1:',count(storekmat_infinite1(:,:,i_elmt)==0.0_CUSTOM_REAL)
    do i = 1,nedof_inf1
      do j = 1,nedof_inf1
        igdof = gdof_elmt1(imap_inf(i))
        jgdof = gdof_elmt1(imap_inf(j))
        if (igdof > 0.and.jgdof > 0) then
          ncount = ncount+1
          row0(ncount) = igdof
          col0(ncount) = jgdof
          grow0(ncount) = ggdof_elmt1(imap_inf(i))
          gcol0(ncount) = ggdof_elmt1(imap_inf(j))
        endif
      enddo
    enddo
  enddo


  ! stage 1: assemble duplicates
  ! sort global indices
  allocate(ind0(ncount),iorder(ncount))
  ind0 = neq1*(row0(1:ncount)-1)+col0(1:ncount)
  call i_uniinv(ind0,iorder)
  nsparse1 = maxval(iorder)

  !debug
  if (myrank == 0) write(*,'(a,1x,i0,1x,a,1x,i0)')'  neq1:',neq1,' Nsparse1:',nsparse1
  call sync_all()

  allocate(krow_sparse1(nsparse1),kcol_sparse1(nsparse1))
  allocate(kgrow_sparse1(nsparse1),kgcol_sparse1(nsparse1))

  !kmat_sparse1=0.0_CUSTOM_REAL
  krow_sparse1 = -1
  kcol_sparse1 = -1
  kgrow_sparse1 = -1
  kgcol_sparse1 = -1
  do i_count = 1,ncount!nmax
    krow_sparse1(iorder(i_count)) = row0(i_count)
    kcol_sparse1(iorder(i_count)) = col0(i_count)
    kgrow_sparse1(iorder(i_count)) = grow0(i_count)
    kgcol_sparse1(iorder(i_count)) = gcol0(i_count)
  enddo
  if (minval(krow_sparse1) < 1.or.minval(kcol_sparse1) < 1.or.                  &
      minval(kgrow_sparse1) < 1.or.minval(kgcol_sparse1) < 1) then
    write(*,*) 'ERROR: local and global indices are less than 1!'
    stop 'Error local and global indices are less than 1'
  endif

  deallocate(row0,col0,grow0,gcol0,kmat0,ind0,iorder)
  deallocate(imap_ic,imap_oc,imap_cm,imap_trinf,imap_inf)

  ! stage 2: assemble across processors

  ! local DOF to global DOF mapping
  allocate(l2gdof1(0:neq1))
  l2gdof1 = -1
  l2gdof1(gdof_ic1) = ggdof_ic1(1,:)
  l2gdof1(gdof_oc1) = ggdof_oc1(1,:)
  l2gdof1(gdof_cm1) = ggdof_cm1(1,:)
  l2gdof1(gdof_trinf1) = ggdof_trinf1(1,:)
  l2gdof1(gdof_inf1) = ggdof_inf1(1,:)

  do i = 1,nsparse1
    if (kgrow_sparse1(i) /= l2gdof1(krow_sparse1(i)).or.kgcol_sparse1(i) /= l2gdof1(kcol_sparse1(i))) then
      print *,'VERY STRANGE!!!!!'
      stop 'Error very strange sparse dof numbers should not occur'
    endif
  enddo

  l2gdof1 = l2gdof1-1 ! PETSC uses 0 indexing
  gmin = minvec(l2gdof1(1:))
  gmax = maxvec(l2gdof1(1:))

  !debug
  if (myrank == 0) write(*,'(a,1x,i0,1x,i0)')'  l2gdof1 range:',gmin,gmax
  call sync_all()

  if (minval(l2gdof1(1:)) < 0) then
    write(*,*) 'ERROR: local-to-global indices are less than 1!'
    stop 'Error local-to-global indices are less than 1'
  endif

  !===============================================================================
  ! Level-2 solver
  !===============================================================================
  if (SOLVER_5GLL) then
    nedof_ic = NEDOFU+NEDOFPHI
    nedof_oc = NEDOFCHI+NEDOFP+NEDOFPHI
    nedof_cm = NEDOFU+NEDOFPHI
    nedof_trinf = NEDOFPHI
    nedof_inf = NEDOFPHI

    nmax = NSPEC_INNER_CORE*(nedof_ic*nedof_ic)+                           &
           NSPEC_OUTER_CORE*(nedof_oc*nedof_oc)+                            &
           NSPEC_CRUST_MANTLE*(nedof_cm*nedof_cm)+                          &
           NSPEC_TRINFINITE*(nedof_trinf*nedof_trinf)+                      &
           NSPEC_INFINITE*(nedof_inf*nedof_inf)
    allocate(col0(nmax),row0(nmax),gcol0(nmax),grow0(nmax),kmat0(nmax))
    !allocate(col0(nmax),row0(nmax),gcol0(nmax),grow0(nmax))

    !debug
    if (myrank == 0) print *,'nedof_ic = ',nedof_ic,nmax

    allocate(imap_ic(nedof_ic),imap_oc(nedof_oc),imap_cm(nedof_cm), &
             imap_trinf(nedof_trinf),imap_inf(nedof_inf))

    imap_ic = (/ (i,i = 1,NGLLCUBE) /) !(/ imapu, imapphi /)
    imap_oc = (/ (i,i = 1,NGLLCUBE) /) !(/ imapchi, imapp, imapphi /)
    imap_cm = (/ (i,i = 1,NGLLCUBE) /) !(/ imapu, imapphi /)
    imap_trinf = (/ (i,i = 1,NGLLCUBE) /) !(/ imapphi /)
    imap_inf = (/ (i,i = 1,NGLLCUBE) /) !(/ imapphi /)

    ! read global degrees of freedoms from DATABASE files
    ! inner core
    write(spm,*)myrank
    fname='DATABASES_MPI/gdof_proc'//trim(adjustl(spm))
    open(10,file=fname,action='read',status='old')
    read(10,*)nglob_ic !NGLOB_INNER_CORE
    allocate(ggdof_ic(NNDOF,nglob_ic))
    read(10,*)ggdof_ic
    read(10,*)nglob_oc !NGLOB_OUTER_CORE
    allocate(ggdof_oc(NNDOF,nglob_oc))
    read(10,*)ggdof_oc
    read(10,*)nglob_cm !NGLOB_CRUST_MANTLE
    allocate(ggdof_cm(NNDOF,nglob_cm))
    read(10,*)ggdof_cm
    read(10,*)nglob_trinf !NGLOB_TRINFINITE
    allocate(ggdof_trinf(NNDOF,nglob_trinf))
    read(10,*)ggdof_trinf
    read(10,*)nglob_inf !NGLOB_INFINITE
    allocate(ggdof_inf(NNDOF,nglob_inf))
    read(10,*)ggdof_inf
    close(10)

    ngdof = maxscal(maxval( (/ maxval(ggdof_ic),maxval(ggdof_oc),maxval(ggdof_cm), &
                               maxval(ggdof_trinf),maxval(ggdof_inf) /) ))

    !debug
    if (myrank == 0) write(*,'(a,i12)')'  Total global degrees of freedom:',ngdof

    ! stage 0: store all elements
    ncount = 0
    ! inner core
    do i_elmt = 1,NSPEC_INNER_CORE
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE)cycle
      gdof_elmt = reshape(gdof_ic(inode_elmt_ic(:,i_elmt)),(/NEDOF/))
      ggdof_elmt = reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
      !if (myrank==0) print *,'ICkmat zeros:',count(storekmat_inner_core(:,:,i_elmt)==0.0_CUSTOM_REAL)
      !if (myrank==0.and.i_elmt==1)print *,'kmat:',storekmat_inner_core(1,:,i_elmt)
      do i = 1,nedof_ic
        do j = 1,nedof_ic
          igdof = gdof_elmt(imap_ic(i))
          jgdof = gdof_elmt(imap_ic(j))
          if (igdof > 0.and.jgdof > 0.and.storekmat_inner_core(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt(imap_ic(i))
            gcol0(ncount) = ggdof_elmt(imap_ic(j))
          endif
        enddo
      enddo
    enddo
    call sync_all
    ! outer core
    do i_elmt = 1,NSPEC_OUTER_CORE
      gdof_elmt = reshape(gdof_oc(inode_elmt_oc(:,i_elmt)),(/NEDOF/))
      ggdof_elmt = reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
      do i = 1,nedof_oc
        do j = 1,nedof_oc
          igdof = gdof_elmt(imap_oc(i))
          jgdof = gdof_elmt(imap_oc(j))
          if (igdof > 0.and.jgdof > 0.and.storekmat_outer_core(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt(imap_oc(i))
            gcol0(ncount) = ggdof_elmt(imap_oc(j))
          endif
        enddo
      enddo
    enddo
    ! crust mantle
    do i_elmt = 1,NSPEC_CRUST_MANTLE
      gdof_elmt = reshape(gdof_cm(inode_elmt_cm(:,i_elmt)),(/NEDOF/))
      ggdof_elmt = reshape(ggdof_cm(:,inode_elmt_cm(:,i_elmt)),(/NEDOF/))
      !if (myrank==0) print *,'CMkmat zeros:',count(storekmat_crust_mantle(:,:,i_elmt)==0.0_CUSTOM_REAL)
      do i = 1,nedof_cm
        do j = 1,nedof_cm
          igdof = gdof_elmt(imap_cm(i))
          jgdof = gdof_elmt(imap_cm(j))
          if (igdof > 0.and.jgdof > 0.and.storekmat_crust_mantle(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt(imap_cm(i))
            gcol0(ncount) = ggdof_elmt(imap_cm(j))
          endif
        enddo
      enddo
    enddo
    ! transition infinite
    do i_elmt = 1,NSPEC_TRINFINITE
      gdof_elmt = reshape(gdof_trinf(inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
      ggdof_elmt = reshape(ggdof_trinf(:,inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
      !if (myrank==0) print *,'TRINFkmat zeros:',count(storekmat_trinfinite(:,:,i_elmt)==0.0_CUSTOM_REAL)
      do i = 1,nedof_trinf
        do j = 1,nedof_trinf
          igdof = gdof_elmt(imap_trinf(i))
          jgdof = gdof_elmt(imap_trinf(j))
          if (igdof > 0.and.jgdof > 0.and.storekmat_trinfinite(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt(imap_trinf(i))
            gcol0(ncount) = ggdof_elmt(imap_trinf(j))
          endif
        enddo
      enddo
    enddo
    ! infinite
    do i_elmt = 1,NSPEC_INFINITE
      gdof_elmt = reshape(gdof_inf(inode_elmt_inf(:,i_elmt)),(/NEDOF/))
      ggdof_elmt = reshape(ggdof_inf(:,inode_elmt_inf(:,i_elmt)),(/NEDOF/))
      !if (myrank==0) print *,'INFkmat zeros:',count(storekmat_infinite(:,:,i_elmt)==0.0_CUSTOM_REAL)
      do i = 1,nedof_inf
        do j = 1,nedof_inf
          igdof = gdof_elmt(imap_inf(i))
          jgdof = gdof_elmt(imap_inf(j))
          if (igdof > 0.and.jgdof > 0) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt(imap_inf(i))
            gcol0(ncount) = ggdof_elmt(imap_inf(j))
          endif
        enddo
      enddo
    enddo

    ! stage 1: assemble duplicates
    ! sort global indices
    allocate(ind0(ncount),iorder(ncount))
    ind0 = neq*(row0(1:ncount)-1)+col0(1:ncount)
    call i_uniinv(ind0,iorder)
    nsparse = maxval(iorder)

    !debug
    if (myrank == 0) write(*,'(a,1x,i0,1x,a,1x,i0)')'  neq:',neq,' Nsparse:',nsparse

    allocate(krow_sparse(nsparse),kcol_sparse(nsparse))
    allocate(kgrow_sparse(nsparse),kgcol_sparse(nsparse))

    krow_sparse = -1
    kcol_sparse = -1
    kgrow_sparse = -1
    kgcol_sparse = -1
    do i_count = 1,ncount
      krow_sparse(iorder(i_count)) = row0(i_count)
      kcol_sparse(iorder(i_count)) = col0(i_count)
      kgrow_sparse(iorder(i_count)) = grow0(i_count)
      kgcol_sparse(iorder(i_count)) = gcol0(i_count)
    enddo
    if (minval(krow_sparse) < 1.or.minval(kcol_sparse) < 1.or.                  &
        minval(kgrow_sparse) < 1.or.minval(kgcol_sparse) < 1) then
      write(*,*) 'ERROR: local and global indices are less than 1!'
      stop 'Error local and global indices are less than 1'
    endif

    deallocate(row0,col0,grow0,gcol0,kmat0,ind0,iorder)
    deallocate(imap_ic,imap_oc,imap_cm,imap_trinf,imap_inf)

    ! stage 2: assemble across processors

    ! local DOF to global DOF mapping
    allocate(l2gdof(0:neq))
    l2gdof = -1
    l2gdof(gdof_ic) = ggdof_ic(1,:)
    l2gdof(gdof_oc) = ggdof_oc(1,:)
    l2gdof(gdof_cm) = ggdof_cm(1,:)
    l2gdof(gdof_trinf) = ggdof_trinf(1,:)
    l2gdof(gdof_inf) = ggdof_inf(1,:)

    l2gdof = l2gdof-1 ! PETSC uses 0 indexing

    !debug
    if (myrank == 0) write(*,'(a,1x,i0,1x,i0)')'  l2gdof range:',minval(l2gdof(1:)),maxval(l2gdof(1:))
    call sync_all()

    if (minval(l2gdof(1:)) < 1) then
      write(*,*) 'ERROR: local-to-global indices are less than 1!'
      stop 'Error local-to-global indices are less than 1'
    endif
  endif !if (SOLVER_5GLL) then

  !debug
  if (myrank == 0) write(*,'(a)')'--------------------------------------------------'

  return

! not used yet...
!contains
!
!  ! subroutine within the prepare_solver_sparse subroutine
!  subroutine is_symmetric(myrank,mat)
!
!  implicit none
!  integer,parameter :: kreal = selected_real_kind(15)
!  integer,intent(in) :: myrank
!  real(kind=kreal),intent(in) :: mat(:,:)
!  real(kind=kreal) :: dx
!  integer :: i,j,n,nc
!
!  n=ubound(mat,1)
!  nc=ubound(mat,2)
!  if (nc /= n) then
!    write(*,*) 'ERROR: non-square matrix!',n,nc
!    stop
!  endif
!
!  do i = 1,n-1
!    do j = i+1,n
!      dx = mat(i,j)-mat(j,i)
!      if (abs(dx) > 1.0e-20_kreal) then
!        if (myrank == 0) write(*,*) 'ERROR: non-symmetric matrix!',mat(i,j),mat(j,i),dx
!        stop
!      endif
!    enddo
!  enddo
!  end subroutine is_symmetric

  end subroutine prepare_solver_sparse

#endif
