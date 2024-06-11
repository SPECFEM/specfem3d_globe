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
  use specfem_par_full_gravity

  implicit none

  ! timing
  double precision :: tstart,tstart0,tCPU
  double precision, external :: wtime

  ! check if anything to do
  if (.not. FULL_GRAVITY) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "preparing full gravity solver"
    if (POISSON_SOLVER == ISOLVER_BUILTIN) then
      write(IMAIN,*) "  Poisson solver: BUILTIN"
      if (CG_SCALING) then
        write(IMAIN,*) "  using CG w/ scaling"
      endif
    else if (POISSON_SOLVER == ISOLVER_PETSC) then
      write(IMAIN,*) "  Poisson solver: PETSc"
    else
      write(IMAIN,*) "  Poisson solver: unknown"
      call exit_MPI(myrank,'Error Poisson solver invalid, POISSON_SOLVER must be 0 or 1')
    endif
    if (USE_POISSON_SOLVER_5GLL) then
      write(IMAIN,*) "  using Level-1 and Level-2 solver"
    else
      write(IMAIN,*) "  using Level-1 solver"
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check PETSc compilation support
  if (POISSON_SOLVER == ISOLVER_PETSC) then
#ifdef USE_PETSC
    ! has compilation support
    continue
#else
    ! compilation without PETSc support
    if (myrank == 0) then
      print *, "Error: PETSc solver enabled without PETSc Support."
      print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
      print *
    endif
    call synchronize_all()
    ! safety stop
    call exit_MPI(myrank,"Error PETSc solver: PETSc solver setup called without compilation support")
#endif
  endif

  ! get MPI starting time
  tstart = wtime()
  tstart0 = tstart

  ! compute and store integration coefficients
  call SIEM_prepare_solver_preintegrate3()
  call SIEM_prepare_solver_preintegrate()

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) "  Elapsed time for preintegration: ",sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  tstart = wtime()

  ! indexify regions (determine inode_elmt_*,inode_map_*,nmir_* arrays)
  call SIEM_get_index_region()

  ! sets the stiffness matrices for Poisson's solver
  ! calculate dprecon
  ! allocates gravload, regional pgravs (e.g. pgrav_cm1)
  call SIEM_prepare_solver_poisson()

  ! create sparse matrix
  call SIEM_prepare_solver_sparse_petsc()

  ! setup seismograms
  call SIEM_prepare_seismos()

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) "  Elapsed time for solver preparation: ",sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  tstart = wtime()

  ! full gravity solver setup for time loop iteration
  call SIEM_prepare_iteration()

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) "  Elapsed time for preparing solver iteration: ",sngl(tCPU),'(s)'
    tCPU = wtime() - tstart0
    write(IMAIN,*) "  Total elapsed time for preparing full gravity solver: ",sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine SIEM_prepare_solver

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_prepare_solver_preintegrate()

  use constants, only: myrank,CUSTOM_REAL,IMAIN,NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,IFLAG_IN_FICTITIOUS_CUBE, &
    USE_POISSON_SOLVER_5GLL

  use specfem_par_crustmantle, only: ibool_crust_mantle,NSPEC_CRUST_MANTLE, &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
    rhostore_crust_mantle

  use specfem_par_innercore, only: idoubling_inner_core,ibool_inner_core,NSPEC_INNER_CORE, &
    xstore_inner_core,ystore_inner_core,zstore_inner_core, &
    rhostore_inner_core

  use specfem_par_outercore, only: ibool_outer_core,NSPEC_OUTER_CORE, &
    xstore_outer_core,ystore_outer_core,zstore_outer_core, &
    rhostore_outer_core

  use specfem_par_full_gravity, only: lagrange_gll, &
    storederiv_cm,storerhojw_cm, &
    storederiv_ic,storerhojw_ic, &
    storederiv_oc,storerhojw_oc

  !not used yet...
  !use specfem_par_full_gravity, only: SAVE_JACOBIAN_ENSIGHT,storedetjac_cm

  use siem_gll_library, only: kdble,NGNOD_INF,zwgljd,dshape_function_hex8,gll_quadrature
  use siem_math_library, only: determinant,invert

  implicit none

  integer :: i,i_elmt,ier
  integer :: ignod(NGNOD_INF)
  real(kind=CUSTOM_REAL) :: detjac,rho(NGLLCUBE)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble

  real(kind=kdble) :: xigll(NGLLX),wxgll(NGLLX),etagll(NGLLY),wygll(NGLLY), &
                      zetagll(NGLLZ),wzgll(NGLLZ), detjac_cm_tmp(NGLLCUBE) !, element_detjac(NGLLX,NGLLY,NGLLZ)

  real(kind=kdble) :: gll_weights(NGLLCUBE),gll_points(NDIM,NGLLCUBE)

  real(kind=kdble),dimension(:,:,:), allocatable :: dshape_hex8,dlagrange_gll

  real(kind=CUSTOM_REAL) :: coord(NGNOD_INF,NDIM)
  real(kind=CUSTOM_REAL) :: jac(NDIM,NDIM)

  double precision :: sizeval

  ! checks if anything to do
  if (.not. USE_POISSON_SOLVER_5GLL) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  preintegrating Level-2 solver arrays"
    call flush_IMAIN()
  endif

  ! estimated memory size required in MB
  ! storederiv_cm,storerhojw_cm
  sizeval = dble(NDIM)*dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
  ! storederiv_ic,storerhojw_ic
  sizeval = sizeval + dble(NDIM)*dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_INNER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_INNER_CORE)*dble(CUSTOM_REAL)
  ! storederiv_oc,storerhojw_oc
  sizeval = sizeval + dble(NDIM)*dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_OUTER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(NGLLCUBE)*dble(NSPEC_OUTER_CORE)*dble(CUSTOM_REAL)
  ! in MB
  sizeval = sizeval / 1024.d0 / 1024.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    size of Level-2 store arrays  = ',sngl(sizeval),'MB'
    write(IMAIN,*) '                                  = ',sngl(sizeval / 1024.d0),'GB'
    call flush_IMAIN()
  endif

  ! array allocations
  allocate(storederiv_cm(NDIM,NGLLCUBE,NGLLCUBE,NSPEC_CRUST_MANTLE), &
           storerhojw_cm(NGLLCUBE,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_cm,.. arrays'
  storederiv_cm(:,:,:,:) = 0.0_CUSTOM_REAL
  storerhojw_cm(:,:) = 0.0_CUSTOM_REAL

  allocate(storederiv_ic(NDIM,NGLLCUBE,NGLLCUBE,NSPEC_INNER_CORE), &
           storerhojw_ic(NGLLCUBE,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_ic,.. arrays'
  storederiv_ic(:,:,:,:) = 0.0_CUSTOM_REAL
  storerhojw_ic(:,:) = 0.0_CUSTOM_REAL

  allocate(storederiv_oc(NDIM,NGLLCUBE,NGLLCUBE,NSPEC_OUTER_CORE), &
           storerhojw_oc(NGLLCUBE,NSPEC_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_oc,.. arrays'
  storederiv_oc(:,:,:,:) = 0.0_CUSTOM_REAL
  storerhojw_oc(:,:) = 0.0_CUSTOM_REAL

  ! allocates local arrays to avoid error about exceeding stack size limit
  allocate(dshape_hex8(NDIM,NGNOD_INF,NGLLCUBE), &
           dlagrange_gll(NDIM,NGLLCUBE,NGLLCUBE))
  dshape_hex8(:,:,:) = 0.0_kdble
  dlagrange_gll(:,:,:) = 0.0_kdble

  ! GLL quadrature
  allocate(lagrange_gll(NGLLCUBE,NGLLCUBE),stat=ier)
  if (ier /= 0) stop 'Error allocating lagrange_gll array'
  lagrange_gll(:,:) = 0.d0

  ! Legendre polynomials
  call zwgljd(xigll,wxgll,NGLLX,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(NDIM,NGNOD_INF,NGLLX,NGLLY,NGLLZ,NGLLCUBE,xigll,etagll,zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,gll_points,gll_weights,lagrange_gll,dlagrange_gll)

  ! inner core
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

    coord(:,1) = xstore_inner_core(ignod(:))
    coord(:,2) = ystore_inner_core(ignod(:))
    coord(:,3) = zstore_inner_core(ignod(:))

    rho(:) = reshape(rhostore_inner_core(:,:,:,i_elmt),(/NGLLCUBE/))

    do i = 1,NGLLCUBE
      jac(:,:) = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_ic(:,:,i,i_elmt) = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      storerhojw_ic(i,i_elmt) = real(rho(i)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
  enddo

  ! outer core
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

    coord(:,1) = xstore_outer_core(ignod(:))
    coord(:,2) = ystore_outer_core(ignod(:))
    coord(:,3) = zstore_outer_core(ignod(:))

    rho(:) = reshape(rhostore_outer_core(:,:,:,i_elmt),(/NGLLCUBE/))

    do i = 1,NGLLCUBE
      jac(:,:) = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_oc(:,:,i,i_elmt) = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      storerhojw_oc(i,i_elmt) = real(rho(i)*detjac*gll_weights(i),kind=CUSTOM_REAL)
    enddo
  enddo

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool_crust_mantle(1,1,1,i_elmt)
    ignod(2) = ibool_crust_mantle(NGLLX,1,1,i_elmt)
    ignod(3) = ibool_crust_mantle(NGLLX,NGLLY,1,i_elmt)
    ignod(4) = ibool_crust_mantle(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool_crust_mantle(1,1,NGLLZ,i_elmt)
    ignod(6) = ibool_crust_mantle(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8) = ibool_crust_mantle(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore_crust_mantle(ignod(:))
    coord(:,2) = ystore_crust_mantle(ignod(:))
    coord(:,3) = zstore_crust_mantle(ignod(:))

    rho(:) = reshape(rhostore_crust_mantle(:,:,:,i_elmt),(/NGLLCUBE/))

    do i = 1,NGLLCUBE
      jac(:,:) = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_cm(:,:,i,i_elmt) = real(matmul(jac,dlagrange_gll(:,i,:)),kind=CUSTOM_REAL)
      storerhojw_cm(i,i_elmt) = real(rho(i)*detjac*gll_weights(i),kind=CUSTOM_REAL)
      detjac_cm_tmp(i) = detjac
    enddo

    !not used yet...
    !if (SAVE_JACOBIAN_ENSIGHT) then
    !  element_detjac(:,:,:) = reshape(detjac_cm_tmp(:), (/NGLLX,NGLLY,NGLLZ/) )
    !  do k = 1, NGLLZ
    !    do j = 1, NGLLY
    !      do i = 1, NGLLX
    !        storedetjac_cm(ibool_crust_mantle(i,j,k,i_elmt)) = element_detjac(i,j,k)
    !      enddo
    !    enddo
    !  enddo
    !endif
  enddo

  ! free temporary arrays
  deallocate(dshape_hex8,dlagrange_gll)

  end subroutine SIEM_prepare_solver_preintegrate

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_prepare_solver_preintegrate3()

  use constants, only: myrank,CUSTOM_REAL,IMAIN,NDIM,NGLLX,NGLLY,NGLLZ, &
    NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF, &
    IFLAG_IN_FICTITIOUS_CUBE

  use specfem_par_crustmantle, only: ibool_crust_mantle,NSPEC_CRUST_MANTLE, &
    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
    rhostore_crust_mantle

  use specfem_par_innercore, only: idoubling_inner_core,ibool_inner_core,NSPEC_INNER_CORE, &
    xstore_inner_core,ystore_inner_core,zstore_inner_core, &
    rhostore_inner_core

  use specfem_par_outercore, only: ibool_outer_core,NSPEC_OUTER_CORE, &
    xstore_outer_core,ystore_outer_core,zstore_outer_core, &
    rhostore_outer_core

  use specfem_par_full_gravity, only: lagrange_gll1, &
    storederiv_cm1,storerhojw_cm1,storejw_cm1, &
    storederiv_ic1,storerhojw_ic1, &
    storederiv_oc1,storerhojw_oc1

  use siem_gll_library, only: kdble,NGNOD_INF,zwgljd,dshape_function_hex8,gll_quadrature
  use siem_math_library, only: determinant,invert

  implicit none

  integer :: i,i_elmt,ier
  integer :: ignod(NGNOD_INF)
  real(kind=CUSTOM_REAL) :: detjac,rho(NGLLCUBE_INF)

  real(kind=kdble),parameter :: jalpha=0.0_kdble,jbeta=0.0_kdble,zero=0.0_kdble

  real(kind=kdble) :: xigll(NGLLX_INF),wxgll(NGLLX_INF),etagll(NGLLY_INF),wygll(NGLLY_INF), &
                      zetagll(NGLLZ_INF),wzgll(NGLLZ_INF)

  real(kind=kdble) :: gll_weights1(NGLLCUBE_INF),gll_points1(NDIM,NGLLCUBE_INF), &
                      dlagrange_gll1(NDIM,NGLLCUBE_INF,NGLLCUBE_INF),dshape_hex8(NDIM,NGNOD_INF,NGLLCUBE_INF)

  real(kind=CUSTOM_REAL) :: coord(NGNOD_INF,NDIM)
  real(kind=CUSTOM_REAL) :: jac(NDIM,NDIM)

  double precision :: sizeval

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  preintegrating Level-1 solver arrays"
    call flush_IMAIN()
  endif

  ! estimated memory size required in MB
  ! storederiv_cm1,storerhojw_cm1,storejw_cm1
  sizeval = dble(NDIM)*dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
  sizeval = sizeval + 2.d0 * dble(NGLLCUBE_INF)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
  ! storederiv_ic1,storerhojw_ic1
  sizeval = sizeval + dble(NDIM)*dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_INNER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_INNER_CORE)*dble(CUSTOM_REAL)
  ! storederiv_oc1,storerhojw_oc1
  sizeval = sizeval + dble(NDIM)*dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_OUTER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NSPEC_OUTER_CORE)*dble(CUSTOM_REAL)
  ! in MB
  sizeval = sizeval / 1024.d0 / 1024.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    size of Level-1 store arrays  = ',sngl(sizeval),'MB'
    write(IMAIN,*) '                                  = ',sngl(sizeval / 1024.d0),'GB'
    call flush_IMAIN()
  endif

  ! array allocations
  allocate(storederiv_cm1(NDIM,NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_CRUST_MANTLE), &
           storerhojw_cm1(NGLLCUBE_INF,NSPEC_CRUST_MANTLE), &
           storejw_cm1(NGLLCUBE_INF,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_cm1,.. arrays'
  storederiv_cm1(:,:,:,:) = 0.0_CUSTOM_REAL
  storerhojw_cm1(:,:) = 0.0_CUSTOM_REAL
  storejw_cm1(:,:) = 0.0_CUSTOM_REAL

  allocate(storederiv_ic1(NDIM,NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_INNER_CORE), &
           storerhojw_ic1(NGLLCUBE_INF,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_ic1,.. arrays'
  storederiv_ic1(:,:,:,:) = 0.0_CUSTOM_REAL
  storerhojw_ic1(:,:) = 0.0_CUSTOM_REAL

  allocate(storederiv_oc1(NDIM,NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_OUTER_CORE), &
           storerhojw_oc1(NGLLCUBE_INF,NSPEC_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_oc1,.. arrays'
  storederiv_oc1(:,:,:,:) = 0.0_CUSTOM_REAL
  storerhojw_oc1(:,:) = 0.0_CUSTOM_REAL

  ! GLL quadrature
  allocate(lagrange_gll1(NGLLCUBE_INF,NGLLCUBE_INF),stat=ier)
  if (ier /= 0) stop 'Error allocating lagrange_gll1 array'
  lagrange_gll1(:,:) = 0.d0

  ! Legendre polynomials
  call zwgljd(xigll,wxgll,NGLLX_INF,jalpha,jbeta)
  call zwgljd(etagll,wygll,NGLLY_INF,jalpha,jbeta)
  call zwgljd(zetagll,wzgll,NGLLZ_INF,jalpha,jbeta)

  ! get derivatives of shape functions for 8-noded hex
  call dshape_function_hex8(NDIM,NGNOD_INF,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,xigll,etagll,zetagll,dshape_hex8)

  ! compute gauss-lobatto-legendre quadrature information
  call gll_quadrature(NDIM,NGLLX_INF,NGLLY_INF,NGLLZ_INF,NGLLCUBE_INF,gll_points1,gll_weights1,lagrange_gll1,dlagrange_gll1)

  ! inner core
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

    coord(:,1) = xstore_inner_core(ignod(:))
    coord(:,2) = ystore_inner_core(ignod(:))
    coord(:,3) = zstore_inner_core(ignod(:))

    rho(:) = reshape(rhostore_inner_core(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt),(/NGLLCUBE_INF/))

    do i = 1,NGLLCUBE_INF
      jac(:,:) = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_ic1(:,:,i,i_elmt) = real(matmul(jac,dlagrange_gll1(:,i,:)),kind=CUSTOM_REAL)
      storerhojw_ic1(i,i_elmt) = real(rho(i)*detjac*gll_weights1(i),kind=CUSTOM_REAL)
    enddo
  enddo

  ! outer core
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

    coord(:,1) = xstore_outer_core(ignod(:))
    coord(:,2) = ystore_outer_core(ignod(:))
    coord(:,3) = zstore_outer_core(ignod(:))

    rho(:) = reshape(rhostore_outer_core(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt), (/NGLLCUBE_INF/))

    do i = 1,NGLLCUBE_INF
      jac(:,:) = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_oc1(:,:,i,i_elmt) = real(matmul(jac,dlagrange_gll1(:,i,:)),kind=CUSTOM_REAL)
      storerhojw_oc1(i,i_elmt) = real(rho(i)*detjac*gll_weights1(i),kind=CUSTOM_REAL)
    enddo
  enddo

  ! crust mantle
  do i_elmt = 1,NSPEC_CRUST_MANTLE
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool_crust_mantle(1,1,1,i_elmt)
    ignod(2) = ibool_crust_mantle(NGLLX,1,1,i_elmt)
    ignod(3) = ibool_crust_mantle(NGLLX,NGLLY,1,i_elmt)
    ignod(4) = ibool_crust_mantle(1,NGLLY,1,i_elmt)
    ! second-last corner nodes
    ignod(5) = ibool_crust_mantle(1,1,NGLLZ,i_elmt)
    ignod(6) = ibool_crust_mantle(NGLLX,1,NGLLZ,i_elmt)
    ignod(7) = ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,i_elmt)
    ignod(8) = ibool_crust_mantle(1,NGLLY,NGLLZ,i_elmt)

    coord(:,1) = xstore_crust_mantle(ignod(:))
    coord(:,2) = ystore_crust_mantle(ignod(:))
    coord(:,3) = zstore_crust_mantle(ignod(:))

    rho(:) = reshape(rhostore_crust_mantle(1:NGLLX:2,1:NGLLY:2,1:NGLLZ:2,i_elmt), (/NGLLCUBE_INF/))

    do i = 1,NGLLCUBE_INF
      jac(:,:) = real(matmul(dshape_hex8(:,:,i),coord),kind=CUSTOM_REAL)
      detjac = determinant(jac)
      call invert(jac)
      storederiv_cm1(:,:,i,i_elmt) = real(matmul(jac,dlagrange_gll1(:,i,:)),kind=CUSTOM_REAL)
      storerhojw_cm1(i,i_elmt) = real(rho(i)*detjac*gll_weights1(i),kind=CUSTOM_REAL)
      storejw_cm1(i,i_elmt) = real(detjac*gll_weights1(i),kind=CUSTOM_REAL)
    enddo
  enddo

  end subroutine SIEM_prepare_solver_preintegrate3

!
!-------------------------------------------------------------------------------
!

! TODO: check this why is it so slow poisson_stiffness is very slow due to loop
! of jacobian and etc. If we store jacobian before hand it should be faster!

  subroutine SIEM_prepare_solver_poisson()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_trinfinite
  use specfem_par_infinite

  use specfem_par_full_gravity

  use siem_poisson, only: poisson_stiffness,poisson_stiffnessINF,poisson_stiffness3, &
                          poisson_stiffnessINF3

  implicit none
  integer :: ier
  double precision :: sizeval,factor_b

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  allocating Poisson Level-1 solver arrays"
    call flush_IMAIN()
  endif

  ! estimated memory size required in MB
  ! storekmat_crust_mantle1,dprecon_crust_mantle1
  sizeval = dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(nnode_cm1)*dble(CUSTOM_REAL)
  ! storekmat_inner_core1,dprecon_inner_core1
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_INNER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(nnode_ic1)*dble(CUSTOM_REAL)
  ! storekmat_outer_core1,dprecon_outer_core1
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_OUTER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(nnode_oc1)*dble(CUSTOM_REAL)
  ! storekmat_trinfinite1,dprecon_trinfinite1
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_TRINFINITE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(nnode_trinf1)*dble(CUSTOM_REAL)
  ! storekmat_infinite1,dprecon_infinite1
  sizeval = sizeval + dble(NGLLCUBE_INF)*dble(NGLLCUBE_INF)*dble(NSPEC_INFINITE)*dble(CUSTOM_REAL)
  sizeval = sizeval + dble(nnode_inf1)*dble(CUSTOM_REAL)
  ! to account for array allocations like b_pgrav_cm,.., b_pgrav_cm1,..
  if (SIMULATION_TYPE == 3) then
    factor_b = 2.d0
  else
    factor_b = 1.d0
  endif
  ! pgrav_cm,..
  sizeval = sizeval + factor_b*dble(NGLOB_CRUST_MANTLE)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(NGLOB_INNER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(NGLOB_OUTER_CORE)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(NGLOB_TRINFINITE)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(NGLOB_INFINITE)*dble(CUSTOM_REAL)
  ! pgrav_cm1,..
  sizeval = sizeval + factor_b*dble(nnode_cm1)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(nnode_ic1)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(nnode_oc1)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(nnode_trinf1)*dble(CUSTOM_REAL)
  sizeval = sizeval + factor_b*dble(nnode_inf1)*dble(CUSTOM_REAL)
  ! dprecon1(0:neq1),gravload1(0:neq1))
  sizeval = sizeval + 2.d0*dble(neq1+1)*dble(CUSTOM_REAL)
  if (SIMULATION_TYPE == 3) then
    ! rho1siem_kl_crust_mantle
    sizeval = sizeval + 3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
    ! rho2siem_kl_crust_mantle
    sizeval = sizeval + 3.d0*3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
  endif
  ! in MB
  sizeval = sizeval / 1024.d0 / 1024.d0

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    size of Level-1 solver arrays = ',sngl(sizeval),'MB'
    write(IMAIN,*) '                                  = ',sngl(sizeval / 1024.d0),'GB'
    call flush_IMAIN()
  endif

  ! Level-1 solver-------------------
  allocate(storekmat_crust_mantle1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_CRUST_MANTLE), &
           dprecon_crust_mantle1(nnode_cm1))
  storekmat_crust_mantle1(:,:,:) = 0.0_CUSTOM_REAL; dprecon_crust_mantle1(:) = 0.0_CUSTOM_REAL

  allocate(storekmat_outer_core1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_OUTER_CORE), &
           dprecon_outer_core1(nnode_oc1))
  storekmat_outer_core1(:,:,:) = 0.0_CUSTOM_REAL; dprecon_outer_core1(:) = 0.0_CUSTOM_REAL

  allocate(storekmat_inner_core1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_INNER_CORE), &
           dprecon_inner_core1(nnode_ic1))
  storekmat_inner_core1(:,:,:) = 0.0_CUSTOM_REAL; dprecon_inner_core1(:) = 0.0_CUSTOM_REAL

  if (ADD_TRINF) then
    allocate(storekmat_trinfinite1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_TRINFINITE), &
             dprecon_trinfinite1(nnode_trinf1))
    storekmat_trinfinite1(:,:,:) = 0.0_CUSTOM_REAL; dprecon_trinfinite1(:) = 0.0_CUSTOM_REAL
  endif
  allocate(storekmat_infinite1(NGLLCUBE_INF,NGLLCUBE_INF,NSPEC_INFINITE), &
           dprecon_infinite1(nnode_inf1))
  storekmat_infinite1(:,:,:) = 0.0_CUSTOM_REAL; dprecon_infinite1(:) = 0.0_CUSTOM_REAL

  ! Allocate gravity perturbation arrays
  allocate(pgrav_ic(NGLOB_INNER_CORE), &
           pgrav_oc(NGLOB_OUTER_CORE), &
           pgrav_cm(NGLOB_CRUST_MANTLE), &
           pgrav_trinf(NGLOB_TRINFINITE), &
           pgrav_inf(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating pgrav_ic,.. arrays'
  pgrav_ic(:) = 0.0_CUSTOM_REAL; pgrav_oc(:) = 0.0_CUSTOM_REAL; pgrav_cm(:) = 0.0_CUSTOM_REAL;
  pgrav_trinf(:) = 0.0_CUSTOM_REAL; pgrav_inf(:) = 0.0_CUSTOM_REAL

  allocate(pgrav_ic1(nnode_ic1), &
           pgrav_oc1(nnode_oc1), &
           pgrav_cm1(nnode_cm1), &
           pgrav_trinf1(nnode_trinf1), &
           pgrav_inf1(nnode_inf1))
  pgrav_ic1(:) = 0.0_CUSTOM_REAL; pgrav_oc1(:) = 0.0_CUSTOM_REAL; pgrav_cm1(:) = 0.0_CUSTOM_REAL;
  pgrav_trinf1(:) = 0.0_CUSTOM_REAL; pgrav_inf1(:) = 0.0_CUSTOM_REAL

  ! Level-1 solver RHS
  allocate(gravload1(0:neq1))
  gravload1(:) = 0.0_CUSTOM_REAL

  ! preconditioner for PCG solver
  allocate(dprecon1(0:neq1))
  dprecon1(:) = 0.0_CUSTOM_REAL

  if (SIMULATION_TYPE == 3) then
    allocate(b_pgrav_ic(NGLOB_INNER_CORE_ADJOINT), &
             b_pgrav_oc(NGLOB_OUTER_CORE_ADJOINT), &
             b_pgrav_cm(NGLOB_CRUST_MANTLE_ADJOINT), &
             b_pgrav_trinf(NGLOB_TRINFINITE_ADJOINT), &
             b_pgrav_inf(NGLOB_INFINITE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating b_pgrav_ic,.. arrays'
    b_pgrav_ic(:) = 0.0_CUSTOM_REAL; b_pgrav_oc(:) = 0.0_CUSTOM_REAL; b_pgrav_cm(:) = 0.0_CUSTOM_REAL;
    b_pgrav_trinf(:) = 0.0_CUSTOM_REAL; b_pgrav_inf(:) = 0.0_CUSTOM_REAL

    allocate(b_pgrav_ic1(nnode_ic1), &
             b_pgrav_oc1(nnode_oc1), &
             b_pgrav_cm1(nnode_cm1), &
             b_pgrav_trinf1(nnode_trinf1), &
             b_pgrav_inf1(nnode_inf1),stat=ier)
    if (ier /= 0) stop 'Error allocating b_pgrav_ic,.. arrays'
    b_pgrav_ic1(:) = 0.0_CUSTOM_REAL; b_pgrav_oc1(:) = 0.0_CUSTOM_REAL; b_pgrav_cm1(:) = 0.0_CUSTOM_REAL;
    b_pgrav_trinf1(:) = 0.0_CUSTOM_REAL; b_pgrav_inf1(:) = 0.0_CUSTOM_REAL

    allocate(b_gravload1(0:neq1),stat=ier)
    if (ier /= 0) stop 'Error allocating b_gravload1 array'
    b_gravload1(:) = 0.0_CUSTOM_REAL
  endif

  ! rotation arrays
  allocate(A_array_rotationL(NGLLCUBE,NSPEC_OUTER_CORE_ROTATION), &
           B_array_rotationL(NGLLCUBE,NSPEC_OUTER_CORE_ROTATION))
  A_array_rotationL(:,:) = 0.0_CUSTOM_REAL; B_array_rotationL(:,:) = 0.0_CUSTOM_REAL

  allocate(A_array_rotationL3(NGLLCUBE_INF,NSPEC_OUTER_CORE_ROTATION), &
           B_array_rotationL3(NGLLCUBE_INF,NSPEC_OUTER_CORE_ROTATION))
  A_array_rotationL3(:,:) = 0.0_CUSTOM_REAL; B_array_rotationL3(:,:) = 0.0_CUSTOM_REAL

  if (SIMULATION_TYPE == 3) then
    allocate(b_A_array_rotationL3(NGLLCUBE_INF,NSPEC_OUTER_CORE_ROTATION), &
             b_B_array_rotationL3(NGLLCUBE_INF,NSPEC_OUTER_CORE_ROTATION))
    b_A_array_rotationL3(:,:) = 0.0_CUSTOM_REAL; b_B_array_rotationL3(:,:) = 0.0_CUSTOM_REAL
  endif

  ! kernels
  if (SIMULATION_TYPE == 3) then
    ! Full gravity density kernels - must be integrated using SIEM after.
    allocate(rho1siem_kl_crust_mantle(3,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT), &
             rho2siem_kl_crust_mantle(3,3,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT))
    rho1siem_kl_crust_mantle(:,:,:,:,:)   = 0._CUSTOM_REAL
    rho2siem_kl_crust_mantle(:,:,:,:,:,:) = 0._CUSTOM_REAL

    ! for debugging
    allocate(debug_rho_kl_cm(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT))
    debug_rho_kl_cm(:,:,:,:,:) = 0.0_CUSTOM_REAL
  endif

  ! setup stiffness matrix (storekmat_*) and preconditioner (dprecon_*) arrays for Level-1 solver
  ! crust mantle
  call poisson_stiffness3(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                          ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                          nnode_cm1,inode_elmt_cm1,storekmat_crust_mantle1,dprecon_crust_mantle1)

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
  call poisson_stiffnessINF3(NSPEC_INFINITE,NGLOB_INFINITE, &
                             ibool_infinite,xstore_infinite,ystore_infinite,zstore_infinite, &
                             nnode_inf1,inode_elmt_inf1,storekmat_infinite1,dprecon_infinite1)

  call synchronize_all()

  ! assemble preconditioner matrices
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

  call synchronize_all()

  ! assemble across the different regions in a process
  dprecon1(:) = zero

  ! crust_mantle
  dprecon1(gdof_cm1(:)) = dprecon1(gdof_cm1(:)) + dprecon_crust_mantle1(:)

  ! outer core
  dprecon1(gdof_oc1(:)) = dprecon1(gdof_oc1(:)) + dprecon_outer_core1(:)

  ! inner core
  dprecon1(gdof_ic1(:)) = dprecon1(gdof_ic1(:)) + dprecon_inner_core1(:)

  ! transition infinite
  if (ADD_TRINF) then
    dprecon1(gdof_trinf1(:)) = dprecon1(gdof_trinf1(:)) + dprecon_trinfinite1(:)
  endif

  ! infinite
  dprecon1(gdof_inf1(:)) = dprecon1(gdof_inf1(:)) + dprecon_infinite1(:)

  dprecon1(0) = 0.0_CUSTOM_REAL

  call synchronize_all()

  ! invert preconditioner
  !dprecon1(1:)=1.0_CUSTOM_REAL/dprecon1(1:)

  !-----------------Level-1 solver

  ! Level-2 solver------------------
  if (USE_POISSON_SOLVER_5GLL) then
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  allocating Poisson Level-2 solver arrays"
      call flush_IMAIN()
    endif

    ! estimated memory size required in MB
    ! storekmat_crust_mantle,dprecon_crust_mantle
    sizeval = dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_CRUST_MANTLE)*dble(CUSTOM_REAL)
    sizeval = sizeval + dble(NGLOB_CRUST_MANTLE)*dble(CUSTOM_REAL)
    ! storekmat_inner_core,dprecon_inner_core
    sizeval = sizeval + dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_INNER_CORE)*dble(CUSTOM_REAL)
    sizeval = sizeval + dble(NGLOB_INNER_CORE)*dble(CUSTOM_REAL)
    ! storekmat_outer_core,dprecon_outer_core
    sizeval = sizeval + dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_OUTER_CORE)*dble(CUSTOM_REAL)
    sizeval = sizeval + dble(NGLOB_OUTER_CORE)*dble(CUSTOM_REAL)
    ! storekmat_trinfinite,dprecon_trinfinite
    sizeval = sizeval + dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_TRINFINITE)*dble(CUSTOM_REAL)
    sizeval = sizeval + dble(NGLOB_TRINFINITE)*dble(CUSTOM_REAL)
    ! storekmat_infinite,dprecon_infinite
    sizeval = sizeval + dble(NGLLCUBE)*dble(NGLLCUBE)*dble(NSPEC_INFINITE)*dble(CUSTOM_REAL)
    sizeval = sizeval + dble(NGLOB_INFINITE)*dble(CUSTOM_REAL)
    ! dprecon(0:neq),gravload(0:neq))
    sizeval = sizeval + 2.d0*dble(neq+1)*dble(CUSTOM_REAL)
    ! in MB
    sizeval = sizeval / 1024.d0 / 1024.d0

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '    size of Level-2 solver arrays = ',sngl(sizeval),'MB'
      write(IMAIN,*) '                                  = ',sngl(sizeval / 1024.d0),'GB'
      call flush_IMAIN()
    endif

    ! Level-2 solver arrays
    allocate(storekmat_crust_mantle(NGLLCUBE,NGLLCUBE,NSPEC_CRUST_MANTLE), &
             dprecon_crust_mantle(NGLOB_CRUST_MANTLE))
    storekmat_crust_mantle(:,:,:) = 0.0_CUSTOM_REAL; dprecon_crust_mantle(:) = 0.0_CUSTOM_REAL

    allocate(storekmat_outer_core(NGLLCUBE,NGLLCUBE,NSPEC_OUTER_CORE), &
             dprecon_outer_core(NGLOB_OUTER_CORE))
    storekmat_outer_core(:,:,:) = 0.0_CUSTOM_REAL; dprecon_outer_core(:) = 0.0_CUSTOM_REAL

    allocate(storekmat_inner_core(NGLLCUBE,NGLLCUBE,NSPEC_INNER_CORE), &
             dprecon_inner_core(NGLOB_INNER_CORE))
    storekmat_inner_core(:,:,:) = 0.0_CUSTOM_REAL; dprecon_inner_core(:) = 0.0_CUSTOM_REAL

    if (ADD_TRINF) then
      allocate(storekmat_trinfinite(NGLLCUBE,NGLLCUBE,NSPEC_TRINFINITE), &
               dprecon_trinfinite(NGLOB_TRINFINITE))
      storekmat_trinfinite(:,:,:) = 0.0_CUSTOM_REAL; dprecon_trinfinite(:) = 0.0_CUSTOM_REAL
    endif
    allocate(storekmat_infinite(NGLLCUBE,NGLLCUBE,NSPEC_INFINITE), &
             dprecon_infinite(NGLOB_INFINITE))
    storekmat_infinite1(:,:,:) = 0.0_CUSTOM_REAL; dprecon_infinite1(:) = 0.0_CUSTOM_REAL

    allocate(dprecon(0:neq),gravload(0:neq))
    dprecon(:) = 0.0_CUSTOM_REAL; gravload(:) = 0.0_CUSTOM_REAL

    ! better to make dprecon_* local rather than global

    ! crust mantle
    call poisson_stiffness(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                           ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                           storekmat_crust_mantle,dprecon_crust_mantle)
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

    call synchronize_all()

    ! assemble stiffness matrices
    ! assemble across the MPI processes in a region
    ! crust_mantle
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_CRUST_MANTLE,dprecon_crust_mantle, &
                             num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                             nibool_interfaces_crust_mantle,ibool_interfaces_crust_mantle, &
                             my_neighbors_crust_mantle)

    ! outer core
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_OUTER_CORE,dprecon_outer_core, &
                             num_interfaces_outer_core,max_nibool_interfaces_oc, &
                             nibool_interfaces_outer_core,ibool_interfaces_outer_core, &
                             my_neighbors_outer_core)

    ! inner core
    call assemble_MPI_scalar(NPROCTOT_VAL,NGLOB_INNER_CORE,dprecon_inner_core, &
                             num_interfaces_inner_core,max_nibool_interfaces_ic, &
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

    call synchronize_all()

    ! assemble across the different regions in a process
    dprecon(:) = zero

    ! crust_mantle
    dprecon(gdof_cm(:)) = dprecon(gdof_cm(:)) + dprecon_crust_mantle(:)

    ! outer core
    dprecon(gdof_oc(:)) = dprecon(gdof_oc(:)) + dprecon_outer_core(:)

    ! inner core
    dprecon(gdof_ic(:)) = dprecon(gdof_ic(:)) + dprecon_inner_core(:)

    ! transition infinite
    if (ADD_TRINF) then
      dprecon(gdof_trinf(:)) = dprecon(gdof_trinf(:)) + dprecon_trinfinite(:)
    endif

    ! infinite
    dprecon(gdof_inf(:)) = dprecon(gdof_inf(:)) + dprecon_infinite(:)

    dprecon(0) = 0.0_CUSTOM_REAL

    call synchronize_all()
    !--------------------Level-2 solver
  endif

  end subroutine SIEM_prepare_solver_poisson

!
!-------------------------------------------------------------------------------
!

! TODO: check this why is it so slow poisson_stiffness is very slow due to loop
! of jacobian and etc. If we store jacobian before hand it may be faster!
! DEVELOPER
!   Hom N Gharti
! HISTORY
!   Sep 30,2013

  subroutine SIEM_prepare_solver_sparse_petsc()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_trinfinite
  use specfem_par_infinite

  use specfem_par_full_gravity

  use siem_math_library, only: i_uniinv
  use siem_math_library_mpi, only: maxscal,minvec,maxvec

  implicit none

  ! sparse stage 0
  integer :: i,j,i_elmt,i_count,ncount
  integer :: igdof,jgdof
  integer :: nmax
  integer :: nedof_ic,nedof_oc,nedof_cm,nedof_trinf,nedof_inf
  integer :: gdof_elmt(NEDOF),ggdof_elmt(NEDOF)
  integer :: nmax1
  integer :: nedof_ic1,nedof_oc1,nedof_cm1,nedof_trinf1,nedof_inf1
  integer :: gdof_elmt1(NEDOF1),ggdof_elmt1(NEDOF1)
  integer,allocatable :: imap_ic(:),imap_oc(:),imap_cm(:),imap_trinf(:),imap_inf(:)
  integer,allocatable :: ind0(:),iorder(:),row0(:),col0(:),grow0(:),gcol0(:)

  integer :: nglob_ic,nglob_oc,nglob_cm,nglob_trinf,nglob_inf
  integer :: nglob_ic1,nglob_oc1,nglob_cm1,nglob_trinf1,nglob_inf1

  character(len=12) :: spm
  character(len=60) :: fname
  integer :: ier
  integer :: gmin,gmax

  ! checks if anything to do
  if (POISSON_SOLVER /= ISOLVER_PETSC) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "  preparing PETSc sparse matrix solver"
    call flush_IMAIN()
  endif

  ! check PETSc compilation support
#ifdef USE_PETSC
  ! has compilation support
  continue
#else
  ! compilation without PETSc support
  if (myrank == 0) then
    print *, "Error: PETSc solver enabled without PETSc Support."
    print *, "       To enable PETSc support, reconfigure with --with-petsc flag."
    print *
  endif
  call synchronize_all()
  ! safety stop
  call exit_MPI(myrank,"Error PETSc solver: PETSc solver setup called without compilation support")
#endif

  !debug
  if (myrank == 0) print *
  if (myrank == 0) print *,'PETSc solver: --------------------------------------------------'

  !===============================================================================
  ! Level-1 solver
  !===============================================================================
  ! Number of DOFs per element in each region
  nedof_ic1 = NEDOFU1 + NEDOFPHI1
  nedof_oc1 = NEDOFCHI1 + NEDOFP1 + NEDOFPHI1
  nedof_cm1 = NEDOFU1 + NEDOFPHI1
  nedof_trinf1 = NEDOFPHI1
  nedof_inf1 = NEDOFPHI1

  ! Maximum DOF in array - number of elements * Element_dof^2
  nmax1 = NSPEC_INNER_CORE*(nedof_ic1*nedof_ic1) + &
          NSPEC_OUTER_CORE*(nedof_oc1*nedof_oc1) + &
          NSPEC_CRUST_MANTLE*(nedof_cm1*nedof_cm1) + &
          NSPEC_TRINFINITE*(nedof_trinf1*nedof_trinf1) + &
          NSPEC_INFINITE*(nedof_inf1*nedof_inf1)

  allocate(col0(nmax1),row0(nmax1),gcol0(nmax1),grow0(nmax1))
  col0(:) = 0; row0(:) = 0
  gcol0(:) = 0; grow0(:) = 0

  !debug
  if (myrank == 0) then
    print *,'PETSc solver: -- Elemental DOFs for IC : ', nedof_ic1
    print *,'PETSc solver: -- Maximum DOFs (nmax1)  : ', nmax1
  endif

  ! Allocate map for each region
  allocate(imap_ic(nedof_ic1),imap_oc(nedof_oc1),imap_cm(nedof_cm1), &
           imap_trinf(nedof_trinf1),imap_inf(nedof_inf1))

  ! initializes arrays with (1,2,3,..,27)
  imap_ic = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapu1, imapphi1 /)
  imap_oc = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapchi1, imapp1, imapphi1 /)
  imap_cm = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapu1, imapphi1 /)
  imap_trinf = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapphi1 /)
  imap_inf = (/ (i,i = 1,NGLLCUBE_INF) /) !(/ imapphi1 /)

  !TODO: check if these imap_* arrays are still valid if multiple degrees of freedom are chosen
  !      current setting in constants.h:
  !        NNDOFU   = 0 ! displacement components
  !        NNDOFCHI = 0 ! displacement potential
  !        NNDOFP   = 0 ! pressure
  !        NNDOFPHI = 1 ! gravitational potential
  !      and with
  !        NEDOF1    = NGLLCUBE_INF * NNDOF
  !      it follows that currently NEDOF1 == NGLLCUBE_INF == nedof_ic1 == nedof_oc1 == .. etc.
  !
  !      for multiple degrees of freedom, NEDOF1 could be different to nedof_ic1 etc.
  !      and these imap_* arrays might not work anymore.
  !
  !      they might need to wrap-around NGLLCUBE_INF, for example: 1,2,..,27,1,2,..,27,1,2,..,27 (?)
  !      this could be done by:
  !        imap_ic = (/ (mod(i-1,NGLLCUBE_INF)+1, i = 1,nedof_ic1) /)
  !      however, this depends on how the gdof_ic1(:),.. arrays are ordered by the SIEM_get_index_region() routine.
  !
  !      for now, this might work only with gravitational potential as only degree of freedom.

  ! safety check
  if (nedof_ic1 /= NGLLCUBE_INF .or. nedof_oc1 /= NGLLCUBE_INF .or. nedof_cm1 /= NGLLCUBE_INF .or. &
      nedof_trinf1 /= NGLLCUBE_INF .or. nedof_inf1 /= NGLLCUBE_INF) then
    print *,'Error: invalid degrees of freedom for Level-1 solver setup!'
    print *,'       nedof_ic1, nedof_oc1, nedof_cm1, nedof_trinf1, nedof_inf1 = ', &
                    nedof_ic1, nedof_oc1, nedof_cm1, nedof_trinf1, nedof_inf1
    print *,'       for now, all should be equal to NGLLCUBE_INF = ',NGLLCUBE_INF
    call exit_MPI(myrank,'Error invalid degrees of freedom for Level-1 solver setup')
  endif

  ! read global degrees of freedoms from DATABASE files (created by running xgindex3D)
  write(spm,*) myrank
  fname='DATABASES_MPI/gdof1_proc'//trim(adjustl(spm))

  open(IIN,file=trim(fname),action='read',status='old',iostat=ier)
  if (ier /= 0 ) then
    print *,'Error: could not open file: ',trim(fname)
    print *,'       Please check that xgindex3D was run prior to this simulation.'
    call exit_MPI(myrank,'Error opening file gdof1_proc_*** for reading')
  endif
  ! inner core
  read(IIN,*) nglob_ic1                   ! Global DOF in inner core
  allocate(ggdof_ic1(NNDOF,nglob_ic1))
  read(IIN,*) ggdof_ic1
  ! outer core
  read(IIN,*) nglob_oc1                   ! Global DOF in outer core
  allocate(ggdof_oc1(NNDOF,nglob_oc1))
  read(IIN,*) ggdof_oc1
  ! crust mantle
  read(IIN,*) nglob_cm1                   ! Global DOF in crust mantle
  allocate(ggdof_cm1(NNDOF,nglob_cm1))
  read(IIN,*) ggdof_cm1
  ! transition
  read(IIN,*) nglob_trinf1                ! Global DOF in transition
  allocate(ggdof_trinf1(NNDOF,nglob_trinf1))
  read(IIN,*) ggdof_trinf1
  ! infinite elements
  read(IIN,*) nglob_inf1                  ! Global DOF in infinite
  allocate(ggdof_inf1(NNDOF,nglob_inf1))
  read(IIN,*) ggdof_inf1
  close(IIN)

  ! Find maximum ID (dof value) for any of the regions
  ngdof1 = maxscal(maxval( (/ maxval(ggdof_ic1),maxval(ggdof_oc1), &
                              maxval(ggdof_cm1),maxval(ggdof_trinf1), &
                              maxval(ggdof_inf1) /) ))

  !debug
  if (myrank == 0) print *,'PETSc solver: -- Total global degrees of freedom1: ',ngdof1

  ! stage 0: store all elements
  ncount = 0

  ! inner core
  do i_elmt = 1,NSPEC_INNER_CORE
    ! suppress fictitious elements in central cube
    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle

    ! Note: gdof_ic1 defined in specfem_par_innercore
    ! Fetch gdof for IC element only on NGLLCUBE_INF points
    gdof_elmt1(:) = reshape(gdof_ic1(inode_elmt_ic1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1(:) = reshape(ggdof_ic1(:,inode_elmt_ic1(:,i_elmt)),(/NEDOF1/))
    do i = 1,nedof_ic1
      do j = 1,nedof_ic1
        igdof = gdof_elmt1(imap_ic(i))
        jgdof = gdof_elmt1(imap_ic(j))
        ! If a degree of freedom and a non-zero kmat:
        if (igdof > 0 .and. jgdof > 0 .and. storekmat_inner_core1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
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
  call synchronize_all()

  ! outer core
  do i_elmt = 1,NSPEC_OUTER_CORE
    gdof_elmt1(:) = reshape(gdof_oc1(inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1(:) = reshape(ggdof_oc1(:,inode_elmt_oc1(:,i_elmt)),(/NEDOF1/))
    do i = 1,nedof_oc1
      do j = 1,nedof_oc1
        igdof = gdof_elmt1(imap_oc(i))
        jgdof = gdof_elmt1(imap_oc(j))
        if (igdof > 0 .and. jgdof > 0 .and. storekmat_outer_core1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
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
    gdof_elmt1(:) = reshape(gdof_cm1(inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1(:) = reshape(ggdof_cm1(:,inode_elmt_cm1(:,i_elmt)),(/NEDOF1/))
    do i = 1,nedof_cm1
      do j = 1,nedof_cm1
        igdof = gdof_elmt1(imap_cm(i))
        jgdof = gdof_elmt1(imap_cm(j))
        if (igdof > 0 .and. jgdof > 0 .and. storekmat_crust_mantle1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
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
  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC_TRINFINITE
      gdof_elmt1(:) = reshape(gdof_trinf1(inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
      ggdof_elmt1(:) = reshape(ggdof_trinf1(:,inode_elmt_trinf1(:,i_elmt)),(/NEDOF1/))
      do i = 1,nedof_trinf1
        do j = 1,nedof_trinf1
          igdof = gdof_elmt1(imap_trinf(i))
          jgdof = gdof_elmt1(imap_trinf(j))
          if (igdof > 0 .and. jgdof > 0 .and. storekmat_trinfinite1(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt1(imap_trinf(i))
            gcol0(ncount) = ggdof_elmt1(imap_trinf(j))
          endif
        enddo
      enddo
    enddo
  endif

  ! infinite
  do i_elmt = 1,NSPEC_INFINITE
    gdof_elmt1(:) = reshape(gdof_inf1(inode_elmt_inf1(:,i_elmt)),(/NEDOF1/))
    ggdof_elmt1(:) = reshape(ggdof_inf1(:,inode_elmt_inf1(:,i_elmt)),(/NEDOF1/))
    do i = 1,nedof_inf1
      do j = 1,nedof_inf1
        igdof = gdof_elmt1(imap_inf(i))
        jgdof = gdof_elmt1(imap_inf(j))
        if (igdof > 0 .and. jgdof > 0) then
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
  ind0(:) = neq1*(row0(1:ncount)-1) + col0(1:ncount)
  call i_uniinv(ind0,iorder)
  nsparse1 = maxval(iorder)

  !debug
  if (myrank == 0) print *,'PETSc solver: -- neq1:',neq1,' Nsparse1:',nsparse1
  call synchronize_all()

  allocate(krow_sparse1(nsparse1),kcol_sparse1(nsparse1))
  krow_sparse1(:) = -1
  kcol_sparse1(:) = -1

  allocate(kgrow_sparse1(nsparse1),kgcol_sparse1(nsparse1))
  kgrow_sparse1(:) = -1
  kgcol_sparse1(:) = -1

  do i_count = 1,ncount!nmax
    krow_sparse1(iorder(i_count)) = row0(i_count)
    kcol_sparse1(iorder(i_count)) = col0(i_count)
    kgrow_sparse1(iorder(i_count)) = grow0(i_count)
    kgcol_sparse1(iorder(i_count)) = gcol0(i_count)
  enddo

  ! check
  if (minval(krow_sparse1) < 1 .or. minval(kcol_sparse1) < 1 .or. &
      minval(kgrow_sparse1) < 1 .or. minval(kgcol_sparse1) < 1) then
    print *,'ERROR: PETSc sparse matrix solver: local and global indices are less than 1!'
    stop 'Error local and global indices are less than 1'
  endif

  deallocate(row0,col0,grow0,gcol0)
  deallocate(ind0,iorder)
  deallocate(imap_ic,imap_oc,imap_cm,imap_trinf,imap_inf)

  ! stage 2: assemble across processors

  ! local DOF to global DOF mapping
  allocate(l2gdof1(0:neq1))
  l2gdof1(:) = -1
  l2gdof1(gdof_ic1(:)) = ggdof_ic1(1,:)
  l2gdof1(gdof_oc1(:)) = ggdof_oc1(1,:)
  l2gdof1(gdof_cm1(:)) = ggdof_cm1(1,:)
  if (ADD_TRINF) l2gdof1(gdof_trinf1(:)) = ggdof_trinf1(1,:)
  l2gdof1(gdof_inf1(:)) = ggdof_inf1(1,:)

  do i = 1,nsparse1
    if (kgrow_sparse1(i) /= l2gdof1(krow_sparse1(i)) .or. kgcol_sparse1(i) /= l2gdof1(kcol_sparse1(i))) then
      print *,'Error: PETSc sparse matrix solver: VERY STRANGE!!!!!'
      stop 'Error very strange sparse dof numbers should not occur'
    endif
  enddo

  l2gdof1(:) = l2gdof1(:) - 1 ! PETSc uses 0 indexing

  gmin = minvec(l2gdof1(1:))
  gmax = maxvec(l2gdof1(1:))

  !debug
  if (myrank == 0) print *,'PETSc solver: -- l2gdof1 range:',gmin,gmax
  call synchronize_all()

  if (minval(l2gdof1(1:)) < 0) then
    print *,'ERROR: PETSc Level-1 solver sparse matrix: local-to-global indices are less than 1!'
    stop 'Error Level-1 solver local-to-global indices are less than 1'
  endif

  !===============================================================================
  ! Level-2 solver
  !===============================================================================
  if (USE_POISSON_SOLVER_5GLL) then
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

    allocate(col0(nmax),row0(nmax),gcol0(nmax),grow0(nmax))
    col0(:) = 0; row0(:) = 0; gcol0(:) = 0; grow0(:) = 0

    !debug
    if (myrank == 0) print *,'PETSc solver: -- nedof_ic = ',nedof_ic,nmax

    allocate(imap_ic(nedof_ic),imap_oc(nedof_oc),imap_cm(nedof_cm), &
             imap_trinf(nedof_trinf),imap_inf(nedof_inf))

    imap_ic = (/ (i,i = 1,NGLLCUBE) /) !(/ imapu, imapphi /)
    imap_oc = (/ (i,i = 1,NGLLCUBE) /) !(/ imapchi, imapp, imapphi /)
    imap_cm = (/ (i,i = 1,NGLLCUBE) /) !(/ imapu, imapphi /)
    imap_trinf = (/ (i,i = 1,NGLLCUBE) /) !(/ imapphi /)
    imap_inf = (/ (i,i = 1,NGLLCUBE) /) !(/ imapphi /)

    ! read global degrees of freedoms from DATABASE files
    ! inner core
    write(spm,*) myrank
    fname='DATABASES_MPI/gdof_proc'//trim(adjustl(spm))

    open(IIN,file=trim(fname),action='read',status='old',iostat=ier)
    if (ier /= 0 ) then
      print *,'Error: could not open file: ',trim(fname)
      print *,'       Please check that xgindex3D was run prior to this simulation.'
      call exit_MPI(myrank,'Error opening file gdof_proc_*** for reading')
    endif
    read(IIN,*)nglob_ic !NGLOB_INNER_CORE
    allocate(ggdof_ic(NNDOF,nglob_ic))
    read(IIN,*)ggdof_ic
    read(IIN,*)nglob_oc !NGLOB_OUTER_CORE
    allocate(ggdof_oc(NNDOF,nglob_oc))
    read(IIN,*)ggdof_oc
    read(IIN,*)nglob_cm !NGLOB_CRUST_MANTLE
    allocate(ggdof_cm(NNDOF,nglob_cm))
    read(IIN,*)ggdof_cm
    read(IIN,*)nglob_trinf !NGLOB_TRINFINITE
    allocate(ggdof_trinf(NNDOF,nglob_trinf))
    read(IIN,*)ggdof_trinf
    read(IIN,*)nglob_inf !NGLOB_INFINITE
    allocate(ggdof_inf(NNDOF,nglob_inf))
    read(IIN,*)ggdof_inf
    close(IIN)

    ngdof = maxscal(maxval( (/ maxval(ggdof_ic),maxval(ggdof_oc),maxval(ggdof_cm), &
                               maxval(ggdof_trinf),maxval(ggdof_inf) /) ))

    !debug
    if (myrank == 0) print *,'PETSc solver: -- Total global degrees of freedom:',ngdof

    ! stage 0: store all elements
    ncount = 0

    ! inner core
    do i_elmt = 1,NSPEC_INNER_CORE
      ! suppress fictitious elements in central cube
      if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle

      gdof_elmt(:) = reshape(gdof_ic(inode_elmt_ic(:,i_elmt)),(/NEDOF/))
      ggdof_elmt(:) = reshape(ggdof_ic(:,inode_elmt_ic(:,i_elmt)),(/NEDOF/))
      do i = 1,nedof_ic
        do j = 1,nedof_ic
          igdof = gdof_elmt(imap_ic(i))
          jgdof = gdof_elmt(imap_ic(j))
          if (igdof > 0 .and. jgdof > 0 .and. storekmat_inner_core(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
            ncount = ncount+1
            row0(ncount) = igdof
            col0(ncount) = jgdof
            grow0(ncount) = ggdof_elmt(imap_ic(i))
            gcol0(ncount) = ggdof_elmt(imap_ic(j))
          endif
        enddo
      enddo
    enddo

    ! outer core
    do i_elmt = 1,NSPEC_OUTER_CORE
      gdof_elmt(:) = reshape(gdof_oc(inode_elmt_oc(:,i_elmt)),(/NEDOF/))
      ggdof_elmt(:) = reshape(ggdof_oc(:,inode_elmt_oc(:,i_elmt)),(/NEDOF/))
      do i = 1,nedof_oc
        do j = 1,nedof_oc
          igdof = gdof_elmt(imap_oc(i))
          jgdof = gdof_elmt(imap_oc(j))
          if (igdof > 0 .and. jgdof > 0 .and. storekmat_outer_core(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
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
      gdof_elmt(:) = reshape(gdof_cm(inode_elmt_cm(:,i_elmt)),(/NEDOF/))
      ggdof_elmt(:) = reshape(ggdof_cm(:,inode_elmt_cm(:,i_elmt)),(/NEDOF/))
      do i = 1,nedof_cm
        do j = 1,nedof_cm
          igdof = gdof_elmt(imap_cm(i))
          jgdof = gdof_elmt(imap_cm(j))
          if (igdof > 0 .and. jgdof > 0 .and. storekmat_crust_mantle(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
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
      gdof_elmt(:) = reshape(gdof_trinf(inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
      ggdof_elmt(:) = reshape(ggdof_trinf(:,inode_elmt_trinf(:,i_elmt)),(/NEDOF/))
      do i = 1,nedof_trinf
        do j = 1,nedof_trinf
          igdof = gdof_elmt(imap_trinf(i))
          jgdof = gdof_elmt(imap_trinf(j))
          if (igdof > 0 .and. jgdof > 0 .and. storekmat_trinfinite(i,j,i_elmt) /= 0.0_CUSTOM_REAL) then
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
      gdof_elmt(:) = reshape(gdof_inf(inode_elmt_inf(:,i_elmt)),(/NEDOF/))
      ggdof_elmt(:) = reshape(ggdof_inf(:,inode_elmt_inf(:,i_elmt)),(/NEDOF/))
      do i = 1,nedof_inf
        do j = 1,nedof_inf
          igdof = gdof_elmt(imap_inf(i))
          jgdof = gdof_elmt(imap_inf(j))
          if (igdof > 0 .and. jgdof > 0) then
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
    ind0(:) = neq*(row0(1:ncount)-1) + col0(1:ncount)
    call i_uniinv(ind0,iorder)
    nsparse = maxval(iorder)

    !debug
    if (myrank == 0) print *,'PETSc solver: -- neq:',neq,' Nsparse:',nsparse

    allocate(krow_sparse(nsparse),kcol_sparse(nsparse))
    krow_sparse(:) = -1
    kcol_sparse(:) = -1

    allocate(kgrow_sparse(nsparse),kgcol_sparse(nsparse))
    kgrow_sparse(:) = -1
    kgcol_sparse(:) = -1

    do i_count = 1,ncount
      krow_sparse(iorder(i_count)) = row0(i_count)
      kcol_sparse(iorder(i_count)) = col0(i_count)
      kgrow_sparse(iorder(i_count)) = grow0(i_count)
      kgcol_sparse(iorder(i_count)) = gcol0(i_count)
    enddo

    ! check
    if (minval(krow_sparse) < 1 .or. minval(kcol_sparse) < 1 .or. &
        minval(kgrow_sparse) < 1 .or. minval(kgcol_sparse) < 1) then
      print *,'ERROR: PETSc sparse matrix solver: local and global indices are less than 1!'
      stop 'Error local and global indices are less than 1'
    endif

    deallocate(row0,col0,grow0,gcol0)
    deallocate(ind0,iorder)
    deallocate(imap_ic,imap_oc,imap_cm,imap_trinf,imap_inf)

    ! stage 2: assemble across processors

    ! local DOF to global DOF mapping
    allocate(l2gdof(0:neq))
    l2gdof(:) = -1
    l2gdof(gdof_ic(:)) = ggdof_ic(1,:)
    l2gdof(gdof_oc(:)) = ggdof_oc(1,:)
    l2gdof(gdof_cm(:)) = ggdof_cm(1,:)
    l2gdof(gdof_trinf(:)) = ggdof_trinf(1,:)
    l2gdof(gdof_inf(:)) = ggdof_inf(1,:)

    l2gdof(:) = l2gdof(:) - 1 ! PETSc uses 0 indexing

    !debug
    if (myrank == 0) print *,'PETSc solver: -- l2gdof range:',minval(l2gdof(1:)),maxval(l2gdof(1:))
    call synchronize_all()

    if (minval(l2gdof(1:)) < 0) then
      print *,'ERROR: PETSc Level-2 solver sparse matrix: local-to-global indices are less than 1!'
      stop 'Error Level-2 solver local-to-global indices are less than 1'
    endif
  endif

  !debug
  if (myrank == 0) print *,'PETSc solver: --------------------------------------------------'
  if (myrank == 0) print *

  end subroutine SIEM_prepare_solver_sparse_petsc

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_prepare_seismos()

  use constants, only: NDIM,NGLLCUBE,CUSTOM_REAL

  use specfem_par, only: SIMULATION_TYPE, FULL_GRAVITY, SAVE_SEISMOGRAMS_STRAIN, ROTATION_VAL, &
    nrec_local,nlength_seismogram, &
    scale_t_inv, scale_displ, scale_veloc

  use specfem_par_full_gravity

  implicit none
  integer :: ier

  ! safety check
  if (.not. FULL_GRAVITY) return

  ! additional scaling factors for gravity seismograms
  scale_accel = scale_veloc * scale_t_inv    ! [m / s^2]
  scale_pgrav = scale_displ**2 * scale_t_inv**2  ! [m^2 / s^2]  ONE NEED TO BE CHECKED!!!

  ! allocate seismogram array
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    ! forward/kernel run
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) "  preparing full gravity seismogram arrays"
      call flush_IMAIN()
    endif

    ! (perturbed) gravitational potential
    allocate(seismograms_phi(1,nrec_local,nlength_seismogram),stat=ier)
    if (ier /= 0) stop 'Error allocating seismograms_phi'
    seismograms_phi = 0.0_CUSTOM_REAL

    ! perturbed gravity
    allocate(seismograms_pgrav(NDIM,nrec_local,nlength_seismogram),stat=ier)
    if (ier /= 0) stop 'Error allocating seismograms_pgrav'
    seismograms_pgrav = 0.0_CUSTOM_REAL

    ! perturbed gravity gradient
    if (SAVE_SEISMOGRAMS_STRAIN) then
      allocate(seismograms_Hgrav(NDIM,NDIM,nrec_local,nlength_seismogram),stat=ier)
      if (ier /= 0) stop 'Error allocating seismograms_Hgrav'
      seismograms_Hgrav = 0.0_CUSTOM_REAL
    endif

    ! Due to background gravity. Free-air change in the gravity (vertical) or tilt of the ground surface(horizontal)
    allocate(seismograms_grav(NDIM,nrec_local,nlength_seismogram),stat=ier)
    if (ier /= 0) stop 'Error allocating seismograms_grav'
    seismograms_grav = 0.0_CUSTOM_REAL

    allocate(g_spec_rec(NDIM,NGLLCUBE,nrec_local),stat=ier)
    if (ier /= 0) stop 'Error allocating g_spec_rec'
    g_spec_rec = 0.0_CUSTOM_REAL

    allocate(gradg_rec(NDIM,nrec_local),stat=ier)
    if (ier /= 0) stop 'Error allocating gradg_rec'
    gradg_rec = 0.0_CUSTOM_REAL

    ! Coriolis acceleration
    if (ROTATION_VAL) then
      allocate(seismograms_corio(NDIM,nrec_local,nlength_seismogram),stat=ier)
      if (ier /= 0) stop 'Error allocating seismograms_corio'
      seismograms_corio = 0.0_CUSTOM_REAL
    endif

    ! pre-computes receiver derivatives
    call SIEM_setup_receivers_precompute_dintp()

    ! compute background g at receivers
    call SIEM_compute_background_gravity_receiver()
  endif

  end subroutine SIEM_prepare_seismos

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_setup_receivers_precompute_dintp()

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,CUSTOM_REAL
  use specfem_par, only: nrec_local,number_receiver_global,ispec_selected_rec, &
    xi_receiver,eta_receiver,gamma_receiver

  use siem_gll_library, only: kdble,NGNOD_INF,dshape_function_hex8_point,gll_lagrange3d_point
  use siem_math_library, only: determinant,invert

  use specfem_par_crustmantle, only: ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  use specfem_par_full_gravity, only: storederiv_rec

  implicit none

  ! local parameters
  integer :: ier
  integer :: irec,irec_local,ispec
  integer :: ignod(NGNOD_INF)
  real(kind=kdble) :: lagrange_gllp(NGLLCUBE),dlagrange_gllp(NDIM,NGLLCUBE),dshape_hex8p(NDIM,NGNOD_INF)
  real(kind=CUSTOM_REAL) :: detjac
  real(kind=CUSTOM_REAL) :: coord(NGNOD_INF,NDIM),jac(NDIM,NDIM),deriv(NDIM,NGLLCUBE)

  allocate(storederiv_rec(NDIM,NGLLCUBE,nrec_local),stat=ier)
  if (ier /= 0) stop 'Error allocating storederiv_rec'
  storederiv_rec = 0.0_CUSTOM_REAL

  do irec_local = 1,nrec_local
    irec = number_receiver_global(irec_local)
    ispec = ispec_selected_rec(irec)

    ! get derivatives of shape functions for 8-noded hex
    call dshape_function_hex8_point(NGNOD_INF,xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec),dshape_hex8p)

    ! compute gauss-lobatto-legendre quadrature information
    call gll_lagrange3d_point(NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE, &
                              xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                              lagrange_gllp,dlagrange_gllp)
    ! crust mantle
    ! EXODUS order NOT indicial order
    ! bottom corner nodes
    ignod(1) = ibool_crust_mantle(1,1,1,ispec)
    ignod(2) = ibool_crust_mantle(NGLLX,1,1,ispec)
    ignod(3) = ibool_crust_mantle(NGLLX,NGLLY,1,ispec)
    ignod(4) = ibool_crust_mantle(1,NGLLY,1,ispec)
    ! second-last corner nodes
    ignod(5) = ibool_crust_mantle(1,1,NGLLZ,ispec)
    ignod(6) = ibool_crust_mantle(NGLLX,1,NGLLZ,ispec)
    ignod(7) = ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,ispec)
    ignod(8) = ibool_crust_mantle(1,NGLLY,NGLLZ,ispec)

    coord(:,1) = xstore_crust_mantle(ignod(:))
    coord(:,2) = ystore_crust_mantle(ignod(:))
    coord(:,3) = zstore_crust_mantle(ignod(:))

    jac = real(matmul(dshape_hex8p(:,:),coord),kind=CUSTOM_REAL)
    detjac = determinant(jac)
    call invert(jac)
    deriv = real(matmul(jac,dlagrange_gllp(:,:)),kind=CUSTOM_REAL)

    ! store for seismos
    storederiv_rec(:,:,irec_local) = deriv(:,:)
  enddo

  end subroutine SIEM_setup_receivers_precompute_dintp

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_compute_background_gravity_receiver()

! This subroutine should be called only once before the time loop to compute the
! background gravity at receiver location for the free air and tilt corrections.

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,CUSTOM_REAL,NR_DENSITY

  !use specfem_par, only: minus_gravity_table
  use specfem_par, only: ELLIPTICITY,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2
  use specfem_par, only: nrec_local,ispec_selected_rec,number_receiver_global

  use specfem_par_crustmantle, only: ibool_crust_mantle,rstore_crust_mantle
    !xstore => xstore_crust_mantle, &
    !ystore => ystore_crust_mantle, &
    !zstore => zstore_crust_mantle

  use specfem_par_full_gravity, only: storederiv_rec,g_spec_rec,gradg_rec

  implicit none

  ! local parameters
  ! phir = perturbed gravitational potential at receiver
  double precision :: g_elm(NGLLCUBE)
  double precision :: radius,theta,phi
  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
  double precision :: g_val,minus_g
  integer :: i,j,k,igll,ispec,iglob,irec_local,irec,ier
  real(kind=CUSTOM_REAL) :: deriv(NDIM,NGLLCUBE),gradg(NDIM)

  ! gravity
  integer :: nspl_gravity
  double precision,dimension(:),allocatable :: rspl_gravity,gravity_spline,gravity_spline2

  ! checks if anything to do
  if (nrec_local == 0) return

  ! helper arrays
  allocate(rspl_gravity(NR_DENSITY), &
           gravity_spline(NR_DENSITY), &
           gravity_spline2(NR_DENSITY),stat=ier)
  if (ier /= 0) stop 'Error allocating gravity helper arrays'

  ! initializes spline coefficients
  rspl_gravity(:) = 0.d0
  gravity_spline(:) = 0.d0
  gravity_spline2(:) = 0.d0

  ! gravity term
  call make_gravity(nspl_gravity,rspl_gravity,gravity_spline,gravity_spline2)

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)
    ispec = ispec_selected_rec(irec)

    ! first determine the nodal ||g||
    igll = 0
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          igll = igll+1
          iglob = ibool_crust_mantle(i,j,k,ispec)

          !radius = dble(xstore(iglob))
          !theta = ystore(iglob)
          !phi = zstore(iglob)

          radius = dble(rstore_crust_mantle(1,iglob))
          theta = dble(rstore_crust_mantle(2,iglob))
          phi = dble(rstore_crust_mantle(3,iglob))

          cos_theta = cos(theta)
          sin_theta = sin(theta)
          cos_phi = cos(phi)
          sin_phi = sin(phi)

          ! puts point locations back into a perfectly spherical shape by removing the ellipticity factor
          ! note: the gravity spline evaluation is done for a perfectly spherical model, thus we remove the ellipicity in case.
          !       however, due to topograpy the radius r might still be > 1.0
          if (ELLIPTICITY) &
            call revert_ellipticity_rtheta(radius,theta,nspl_ellip,rspl_ellip,ellipicity_spline,ellipicity_spline2)

          if (radius > 1.d0) radius = 1.d0

          ! get g
          call spline_evaluation(rspl_gravity,gravity_spline,gravity_spline2,nspl_gravity,radius,g_val)

          !debug
          !print *,'debug: gravity: receiver = ',irec_local,'r = ',radius,'g = ',g_val

          ! spherical components of the gravitational acceleration
          ! for efficiency replace with lookup table every 100 m in radial direction
          !int_radius = nint(radius*R_EARTH_KM*10.d0)
          !minus_g = minus_gravity_table(int_radius)
          !g_elm(igll) = abs(minus_g)

          minus_g = - g_val
          g_elm(igll) = abs(minus_g)

          ! store for seismos
          ! Cartesian components of the gravitational acceleration
          g_spec_rec(1,igll,irec_local) = real(minus_g * sin_theta * cos_phi,kind=CUSTOM_REAL) !gxl
          g_spec_rec(2,igll,irec_local) = real(minus_g * sin_theta * sin_phi,kind=CUSTOM_REAL) !gyl
          g_spec_rec(3,igll,irec_local) = real(minus_g * cos_theta,kind=CUSTOM_REAL)           !gzl
        enddo
      enddo
    enddo

    deriv = storederiv_rec(:,:,irec_local)
    gradg = real(matmul(deriv,g_elm),kind=CUSTOM_REAL)

    ! store for seismos
    gradg_rec(:,irec_local) = gradg(:)
  enddo

  ! free temporary arrays
  deallocate(rspl_gravity,gravity_spline,gravity_spline2)

  end subroutine SIEM_compute_background_gravity_receiver
