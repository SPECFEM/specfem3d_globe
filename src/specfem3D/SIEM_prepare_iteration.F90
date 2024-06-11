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


  subroutine SIEM_prepare_iteration()

! preparation for full gravity simulation

  use specfem_par
  use specfem_par_full_gravity

  use siem_solver_petsc, only: petsc_initialize1,petsc_initialize, &
                               petsc_set_matrix1,petsc_set_matrix, &
                               petsc_zero_initialguess1,petsc_zero_backwards_initialguess1

  implicit none
  integer :: ier

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  preparing solver iteration:'
    call flush_IMAIN()
  endif

  ! Allocate gravity perturbation arrays
  allocate(pgrav1(0:neq1),stat=ier)
  if (ier /= 0) stop 'Error allocating pgrav1 array'
  pgrav1(:) = 0.0_CUSTOM_REAL

  allocate(pgrav(0:neq),stat=ier)
  if (ier /= 0) stop 'Error allocating pgrav array'
  pgrav(:) = 0.0_CUSTOM_REAL

  if (SIMULATION_TYPE == 3) then
    allocate(b_pgrav1(0:neq1)); b_pgrav1(:) = 0.0_CUSTOM_REAL
    allocate(b_pgrav(0:neq));   b_pgrav(:) = 0.0_CUSTOM_REAL
  endif

  ! INITIALISE SOLVER
  if (POISSON_SOLVER == ISOLVER_BUILTIN .and. CG_SCALING) then
    ! builtin solver
    call initialise_inbuilt_solver()
  endif

  if (POISSON_SOLVER == ISOLVER_PETSC) then
    ! petsc
    call petsc_initialize1()
    call petsc_set_matrix1()

    if (USE_POISSON_SOLVER_5GLL) then
      ! Level-2 solver
      call petsc_initialize()
      call petsc_set_matrix()
    endif
  endif

  ! compute background gravity by using spectral-infinite-element method
  ! Note that the background gravity is already computed by using semi-analytical
  ! method in make_gravity.f90
  ! compute background gravitational field only once
  ! WARNING: We have not currently using this for the time marching.
  if (NUMBER_OF_THIS_RUN == 1) call SIEM_compute_background_gravity()

  ! initialize for dynamic solver
  pgrav1(:) = 0.0_CUSTOM_REAL
  pgrav(:) = 0.0_CUSTOM_REAL

  ! reset initial guess to zero otherwise due to the large values for the
  ! background gravity we will get very absurd values in timestepping solutions
  if (POISSON_SOLVER == ISOLVER_PETSC) then
    call petsc_zero_initialguess1()

    if (SIMULATION_TYPE == 3) then
      call petsc_zero_backwards_initialguess1()
    endif
  endif

  call synchronize_all()

  end subroutine SIEM_prepare_iteration

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_interpolate_gravity()

! interpolates Level-1 forward pgrav1 values onto Level-2 pgrav

  use constants, only: myrank,ADD_TRINF

  use constants_solver, only: NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE, &
    NSPEC_TRINFINITE,NSPEC_INFINITE, &
    NGLOB_CRUST_MANTLE,NGLOB_OUTER_CORE,NGLOB_INNER_CORE,NGLOB_TRINFINITE,NGLOB_INFINITE

  use specfem_par_full_gravity, only: is_active_gll,igll_active_on, &
    inode_elmt_ic,inode_elmt_oc,inode_elmt_cm,inode_elmt_trinf,inode_elmt_inf, &
    inode_map_ic,inode_map_oc,inode_map_cm,inode_map_trinf,inode_map_inf, &
    nmir_ic,nmir_oc,nmir_cm,nmir_trinf,nmir_inf, &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    gdof_ic1,gdof_oc1,gdof_cm1,gdof_trinf1,gdof_inf1

  use specfem_par_full_gravity, only: &
    pgrav1,pgrav_ic1,pgrav_oc1,pgrav_cm1,pgrav_trinf1,pgrav_inf1, &
    pgrav_ic,pgrav_oc,pgrav_cm,pgrav_trinf,pgrav_inf

  use siem_solver_mpi, only: interpolate3to5

  implicit none

  ! interpolate the gravity values
  ! transfer to array for each region (cm, ic, oc)
  pgrav_ic1(:) = pgrav1(gdof_ic1(:))
  pgrav_oc1(:) = pgrav1(gdof_oc1(:))
  pgrav_cm1(:) = pgrav1(gdof_cm1(:))
  if (ADD_TRINF) pgrav_trinf1(:) = pgrav1(gdof_trinf1(:))
  pgrav_inf1(:) = pgrav1(gdof_inf1(:))

  ! interpolate values for each region
  call interpolate3to5(NSPEC_INNER_CORE,NGLOB_INNER_CORE,nnode_ic1, &
                       inode_elmt_ic,nmir_ic,inode_map_ic,is_active_gll,igll_active_on,pgrav_ic1,pgrav_ic)

  call interpolate3to5(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,nnode_oc1, &
                       inode_elmt_oc,nmir_oc,inode_map_oc,is_active_gll,igll_active_on,pgrav_oc1,pgrav_oc)

  call interpolate3to5(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,nnode_cm1, &
                       inode_elmt_cm,nmir_cm,inode_map_cm,is_active_gll,igll_active_on,pgrav_cm1,pgrav_cm)

  if (ADD_TRINF) then
    call interpolate3to5(NSPEC_TRINFINITE,NGLOB_TRINFINITE,nnode_trinf1, &
                         inode_elmt_trinf,nmir_trinf,inode_map_trinf,is_active_gll,igll_active_on,pgrav_trinf1,pgrav_trinf)
  endif

  call interpolate3to5(NSPEC_INFINITE,NGLOB_INFINITE,nnode_inf1, &
                       inode_elmt_inf,nmir_inf,inode_map_inf,is_active_gll,igll_active_on,pgrav_inf1,pgrav_inf)

  !debug
  if (myrank == 0) print *,' Finished gravity interpolation of loaded (saved) forward wavefield'

  end subroutine SIEM_interpolate_gravity

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_interpolate_backward_gravity()

! interpolates Level-1 backward b_pgrav1 values onto Level-2 b_pgrav

  use constants, only: myrank,ADD_TRINF

  use constants_solver, only: NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE, &
    NSPEC_TRINFINITE,NSPEC_INFINITE, &
    NGLOB_CRUST_MANTLE,NGLOB_OUTER_CORE,NGLOB_INNER_CORE,NGLOB_TRINFINITE,NGLOB_INFINITE

  use specfem_par_full_gravity, only: is_active_gll,igll_active_on, &
    inode_elmt_ic,inode_elmt_oc,inode_elmt_cm,inode_elmt_trinf,inode_elmt_inf, &
    inode_map_ic,inode_map_oc,inode_map_cm,inode_map_trinf,inode_map_inf, &
    nmir_ic,nmir_oc,nmir_cm,nmir_trinf,nmir_inf, &
    nnode_ic1,nnode_oc1,nnode_cm1,nnode_trinf1,nnode_inf1, &
    gdof_ic1,gdof_oc1,gdof_cm1,gdof_trinf1,gdof_inf1

  use specfem_par_full_gravity, only: &
    b_pgrav1,b_pgrav_ic1,b_pgrav_oc1,b_pgrav_cm1,b_pgrav_trinf1,b_pgrav_inf1, &
    b_pgrav_ic,b_pgrav_oc,b_pgrav_cm,b_pgrav_trinf,b_pgrav_inf

  use siem_solver_mpi, only: interpolate3to5

  implicit none

  ! interpolate the gravity values
  ! transfer to array for each region (cm, ic, oc)
  b_pgrav_ic1(:) = b_pgrav1(gdof_ic1(:))
  b_pgrav_oc1(:) = b_pgrav1(gdof_oc1(:))
  b_pgrav_cm1(:) = b_pgrav1(gdof_cm1(:))
  if (ADD_TRINF) b_pgrav_trinf1(:) = b_pgrav1(gdof_trinf1(:))
  b_pgrav_inf1(:) = b_pgrav1(gdof_inf1(:))

  ! interpolate values for each region
  call interpolate3to5(NSPEC_INNER_CORE,NGLOB_INNER_CORE,nnode_ic1, &
                       inode_elmt_ic,nmir_ic,inode_map_ic,is_active_gll,igll_active_on,b_pgrav_ic1,b_pgrav_ic)

  call interpolate3to5(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,nnode_oc1, &
                       inode_elmt_oc,nmir_oc,inode_map_oc,is_active_gll,igll_active_on,b_pgrav_oc1,b_pgrav_oc)

  call interpolate3to5(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,nnode_cm1, &
                       inode_elmt_cm,nmir_cm,inode_map_cm,is_active_gll,igll_active_on,b_pgrav_cm1,b_pgrav_cm)

  if (ADD_TRINF) then
    call interpolate3to5(NSPEC_TRINFINITE,NGLOB_TRINFINITE,nnode_trinf1, &
                         inode_elmt_trinf,nmir_trinf,inode_map_trinf,is_active_gll,igll_active_on,b_pgrav_trinf1,b_pgrav_trinf)
  endif

  call interpolate3to5(NSPEC_INFINITE,NGLOB_INFINITE,nnode_inf1, &
                       inode_elmt_inf,nmir_inf,inode_map_inf,is_active_gll,igll_active_on,b_pgrav_inf1,b_pgrav_inf)

  !debug
  if (myrank == 0) print *,' Finished gravity interpolation of loaded (saved) backward wavefield'

  end subroutine SIEM_interpolate_backward_gravity

!
!-------------------------------------------------------------------------------
!

  subroutine initialise_inbuilt_solver()

  use constants, only: myrank,IMAIN,NGLLCUBE_INF,CUSTOM_REAL, &
    ADD_TRINF,USE_POISSON_SOLVER_5GLL

  use specfem_par, only: NSPEC_INNER_CORE,NSPEC_OUTER_CORE,NSPEC_CRUST_MANTLE, &
                         NSPEC_TRINFINITE,NSPEC_INFINITE

  use specfem_par_full_gravity, only: dprecon1, dprecon, neq, neq1, ndscale, ndscale1, &
    storekmat_crust_mantle1, inode_elmt_cm1, gdof_cm1, &
    storekmat_inner_core1, inode_elmt_ic1, gdof_ic1, &
    storekmat_outer_core1, inode_elmt_oc1, gdof_oc1, &
    storekmat_trinfinite1, inode_elmt_trinf1, gdof_trinf1, &
    storekmat_infinite1, inode_elmt_inf1, gdof_inf1

  implicit none

  ! IO/Global variables to be accounted for:
  integer :: igdof1(NGLLCUBE_INF)

  ! Local variables:
  integer :: istat, i_elmt, i,j
  integer :: inodes1(NGLLCUBE_INF)
  real(kind=CUSTOM_REAL) :: val

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '    initializing BUILTIN solver'
    call flush_IMAIN()
  endif

  !debug
  if (myrank == 0) print *,'builtin solver: Zeros:',count(dprecon1(1:) == 0.0_CUSTOM_REAL)

  ! TODO: create separate function for this

  ! Level-1 solver scaling
  allocate(ndscale1(0:neq1),stat=istat)
  if (istat /= 0) stop 'Error allocating ndscale1 array'
  ndscale1(:) = 1.0_CUSTOM_REAL

  ! scaling based on inverted preconditioner values
  do i = 1,neq1
    val = sqrt(abs(dprecon1(i)))
    if (val /= 0.0_CUSTOM_REAL) ndscale1(i) = 1.0_CUSTOM_REAL / val
  enddo

  ! scales stiffness matrices (storekmat_*)
  do i_elmt = 1,NSPEC_INNER_CORE
    inodes1 = inode_elmt_ic1(:,i_elmt)
    igdof1  = gdof_ic1(inodes1)
    do i = 1,NGLLCUBE_INF
      do j = 1,NGLLCUBE_INF
        storekmat_inner_core1(i,j,i_elmt) = ndscale1(igdof1(i)) * &
                                            storekmat_inner_core1(i,j,i_elmt) * ndscale1(igdof1(j))
      enddo
    enddo
  enddo

  do i_elmt = 1,NSPEC_OUTER_CORE
    inodes1 = inode_elmt_oc1(:,i_elmt)
    igdof1 = gdof_oc1(inodes1)
    do i = 1,NGLLCUBE_INF
      do j = 1,NGLLCUBE_INF
        storekmat_outer_core1(i,j,i_elmt) = ndscale1(igdof1(i)) * &
                                            storekmat_outer_core1(i,j,i_elmt)*ndscale1(igdof1(j))
      enddo
    enddo
  enddo

  do i_elmt = 1,NSPEC_CRUST_MANTLE
    inodes1 = inode_elmt_cm1(:,i_elmt)
    igdof1 = gdof_cm1(inodes1)
    do i = 1,NGLLCUBE_INF
      do j = 1,NGLLCUBE_INF
        storekmat_crust_mantle1(i,j,i_elmt) = ndscale1(igdof1(i)) * &
                                              storekmat_crust_mantle1(i,j,i_elmt)*ndscale1(igdof1(j))
      enddo
    enddo
  enddo

  if (ADD_TRINF) then
    do i_elmt = 1,NSPEC_TRINFINITE
      inodes1 = inode_elmt_trinf1(:,i_elmt)
      igdof1 = gdof_trinf1(inodes1)
      do i = 1,NGLLCUBE_INF
        do j = 1,NGLLCUBE_INF
          storekmat_trinfinite1(i,j,i_elmt) = ndscale1(igdof1(i)) * &
                                              storekmat_trinfinite1(i,j,i_elmt)*ndscale1(igdof1(j))
        enddo
      enddo
    enddo
  endif

  do i_elmt = 1,NSPEC_INFINITE
    inodes1 = inode_elmt_inf1(:,i_elmt)
    igdof1 = gdof_inf1(inodes1)
    do i = 1,NGLLCUBE_INF
      do j = 1,NGLLCUBE_INF
        storekmat_infinite1(i,j,i_elmt) = ndscale1(igdof1(i)) * &
                                          storekmat_infinite1(i,j,i_elmt)*ndscale1(igdof1(j))
      enddo
    enddo
  enddo

  !debug
  if (myrank == 0) print *,'builtin solver: scaling min/max = ',minval(ndscale1),'/',maxval(ndscale1)
  if (myrank == 0) print *
  call synchronize_all()

  ! Level-2 solver scaling
  if (USE_POISSON_SOLVER_5GLL) then
    allocate(ndscale(0:neq),stat=istat)
    if (istat /= 0) stop 'Error allocating ndscale array'
    ndscale(:) = 1.0_CUSTOM_REAL

    do i = 1,neq
      val = sqrt(abs(dprecon(i)))
      if (val /= 0.0_CUSTOM_REAL) ndscale(i) = 1.0_CUSTOM_REAL / val
    enddo

    !TODO: check scaling for storekmat_* arrays?
  endif

  end subroutine initialise_inbuilt_solver

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_compute_background_gravity()

  use specfem_par

  !ensight vis
  !use specfem_par_crustmantle, only: ibool_crust_mantle, &
  !  xstore_crust_mantle, ystore_crust_mantle, zstore_crust_mantle
  !use specfem_par_innercore, only: ibool_inner_core,idoubling_inner_core, &
  !  xstore_inner_core, ystore_inner_core, zstore_inner_core
  !use specfem_par_outercore, only: ibool_outer_core, &
  !  xstore_outer_core, ystore_outer_core, zstore_outer_core

  use specfem_par_full_gravity, only: gdof_cm, gdof_cm1, inode_elmt_cm, inode_map_cm, nmir_cm, &
    gdof_ic, gdof_ic1, inode_elmt_ic, inode_map_ic, nmir_ic, &
    gdof_oc, gdof_oc1, inode_elmt_oc, inode_map_oc, nmir_oc, &
    gdof_trinf, gdof_trinf1, inode_elmt_trinf, inode_map_trinf, nmir_trinf, &
    gdof_inf, gdof_inf1, inode_elmt_inf, inode_map_inf, nmir_inf, &
    nnode_cm1,nnode_ic1,nnode_oc1,nnode_trinf1,nnode_inf1

  use specfem_par_full_gravity, only: ndscale1,ndscale,dprecon1,dprecon, &
    gravload1,pgrav1,pgrav_ic1,pgrav_oc1,pgrav_cm1,pgrav_trinf1,pgrav_inf1, &
    gravload,pgrav,pgrav_ic,pgrav_oc,pgrav_cm,pgrav_trinf,pgrav_inf, &
    neq1,neq,CG_SCALING, &
    is_active_gll,igll_active_on

  use siem_solver_petsc, only: petsc_zero_initialguess1,petsc_set_vector1,petsc_solve1, &
    petsc_set_vector,petsc_solve

  use siem_poisson, only: compute_poisson_rhoload,compute_poisson_rhoload3,poisson_gravity
  use siem_solver_mpi, only: cg_solver,cg_solver3,diagpcg_solver,diagpcg_solver3, interpolate3to5

  implicit none
  ! Local variables
  integer :: cg_iter
  ! timing
  double precision :: tstart,tCPU
  double precision, external :: wtime

  ! ensight vis
  !integer :: i_elmt, i_node, ier
  !integer :: npart
  !integer :: num(NGLLCUBE)
  !real, allocatable :: varvar(:)
  !integer :: varvardim
  !real(kind=CUSTOM_REAL) :: id_ic(NSPEC_INNER_CORE)
  !real(kind=CUSTOM_REAL) :: grav0_ic(NGLOB_INNER_CORE), grav0_oc(NGLOB_OUTER_CORE), grav0_cm(NGLOB_CRUST_MANTLE)
  !real(kind=CUSTOM_REAL),dimension(:,:,:), allocatable :: gradphi_cm,gradphi_oc,gradphi_ic
  !real(kind=CUSTOM_REAL),dimension(:,:), allocatable :: ngradphi_cm,ngradphi_oc,ngradphi_ic
  !integer,dimension(:),allocatable :: nvalency_ic,nvalency_oc,nvalency_cm
  !
  !character(len=80) :: file_head
  !character(len=250)::  grav_file
  !character(len=80),parameter :: ensight_etype='hexa8'
  !character(len=80) :: destag ! this must be 80 characters long
  !character(len=20) :: ptail

  ! interfaces for visualization routine
  interface
    subroutine write_ensight_pernodeSCALAS(out_fname,destag,npart,n,var)
    implicit none
    character(len=250),intent(in) :: out_fname
    character(len=80),intent(in) :: destag
    integer,intent(in) :: npart,n
    real,dimension(:),intent(in) :: var
    end subroutine write_ensight_pernodeSCALAS
    !============================================
    subroutine write_ensight_perelementAS(out_fname,etype,destag,npart,var)
    implicit none
    character(len=250),intent(in) :: out_fname
    character(len=20),intent(in) :: etype
    character(len=80),intent(in) :: destag
    integer,intent(in) :: npart
    real,dimension(:),intent(in) :: var
    end subroutine write_ensight_perelementAS
    !============================================
  end interface

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '    computing SIEM background gravity'
    call flush_IMAIN()
  endif

  ! timing
  tstart = wtime()

  ! Level-1 solver
  if (myrank == 0) print *,'----- SOLVE BACKGROUND GRAVITY -----'
  if (myrank == 0) print *,'Level 1 solver starts...'

  ! Calculate the RHS (gravload1)
  call compute_poisson_rhoload3()

  if (POISSON_SOLVER == ISOLVER_BUILTIN) then
    ! built-in solver
    if (CG_SCALING) then
      gravload1(:) = ndscale1(:) * gravload1(:)

      call cg_solver3(myrank,neq1,pgrav1,gravload1,cg_iter)
      if (myrank == 0) print *,'CG solver: L1: solve done'

      pgrav1(:) = ndscale1(:) * pgrav1(:)
    else
      call diagpcg_solver3(myrank,neq1,pgrav1,gravload1,dprecon1,cg_iter)
      if (myrank == 0) print *,'PCG solver: L1: solve done'
    endif
  else
    ! petsc solver
    call petsc_set_vector1(gravload1)
    if (myrank == 0) print *,'PETSc solver: L1: Set RHS with PETSC done'

    call petsc_zero_initialguess1()

    call petsc_solve1(pgrav1(1:),cg_iter)
    if (myrank == 0) print *,'PETSc solver: L1: Solve with PETSC done'

    call synchronize_all()
  endif

  !debug
  if (myrank == 0) print *,'---- Level 1 solver complete! ','iterations: ',cg_iter
  if (myrank == 0) print *

  ! transfer to array for each region (cm, ic, oc)
  pgrav_ic1(:) = pgrav1(gdof_ic1(:))
  pgrav_oc1(:) = pgrav1(gdof_oc1(:))
  pgrav_cm1(:) = pgrav1(gdof_cm1(:))
  if (ADD_TRINF) pgrav_trinf1(:) = pgrav1(gdof_trinf1(:))
  pgrav_inf1(:) = pgrav1(gdof_inf1(:))

  ! interpolate values for each region
  call interpolate3to5(NSPEC_INNER_CORE,NGLOB_INNER_CORE,nnode_ic1, &
                       inode_elmt_ic,nmir_ic,inode_map_ic,is_active_gll,igll_active_on,pgrav_ic1,pgrav_ic)

  call interpolate3to5(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,nnode_oc1, &
                       inode_elmt_oc,nmir_oc,inode_map_oc,is_active_gll,igll_active_on,pgrav_oc1,pgrav_oc)

  call interpolate3to5(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,nnode_cm1, &
                       inode_elmt_cm,nmir_cm,inode_map_cm,is_active_gll,igll_active_on,pgrav_cm1,pgrav_cm)

  if (ADD_TRINF) then
    call interpolate3to5(NSPEC_TRINFINITE,NGLOB_TRINFINITE,nnode_trinf1, &
                         inode_elmt_trinf,nmir_trinf,inode_map_trinf,is_active_gll,igll_active_on,pgrav_trinf1,pgrav_trinf)
  endif

  call interpolate3to5(NSPEC_INFINITE,NGLOB_INFINITE,nnode_inf1, &
                       inode_elmt_inf,nmir_inf,inode_map_inf,is_active_gll,igll_active_on,pgrav_inf1,pgrav_inf)

  ! --> Finished interpolation

  ! ensight vis
  ! allocates arrays
  !allocate(gradphi_cm(NDIM,NGLLCUBE,NSPEC_CRUST_MANTLE), &
  !         gradphi_oc(NDIM,NGLLCUBE,NSPEC_OUTER_CORE), &
  !         gradphi_ic(NDIM,NGLLCUBE,NSPEC_INNER_CORE),stat=ier)
  !if (ier /= 0) stop 'Error allocating gradphi_cm,.. arrays'
  !gradphi_cm(:,:,:) = 0.0_CUSTOM_REAL; gradphi_oc(:,:,:) = 0.0_CUSTOM_REAL; gradphi_ic(:,:,:) = 0.0_CUSTOM_REAL
  !
  !! compute gradphi only for fluid/solid regions
  !call poisson_gravity(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
  !                     ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                     pgrav_cm,gradphi_cm)
  !
  !call poisson_gravity(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
  !                     ibool_outer_core,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                     pgrav_oc,gradphi_oc)
  !
  !call poisson_gravity(IREGION_INNER_CORE,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
  !                     ibool_inner_core,xstore_inner_core,ystore_inner_core,zstore_inner_core, &
  !                     pgrav_ic,gradphi_ic)
  !
  ! create mapping for 4 corner nodes required for plotting
  !call create_full_gravity_mapping()
  !
  ! compute global gravity
  !allocate(ngradphi_cm(NDIM,NGLOB_CRUST_MANTLE), &
  !         ngradphi_oc(NDIM,NGLOB_OUTER_CORE), &
  !         ngradphi_ic(NDIM,NGLOB_INNER_CORE),stat=ier)
  !if (ier /= 0) stop 'Error allocating ngradphi_cm,.. arrays'
  !ngradphi_cm(:,:) = 0.0_CUSTOM_REAL; ngradphi_oc(:,:) = 0.0_CUSTOM_REAL; ngradphi_ic(:,:) = 0.0_CUSTOM_REAL
  !! Inner core
  !ngradphi_ic(:,:) = zero
  !do i_elmt = 1,NSPEC_INNER_CORE
  !  if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
  !  num(:) = inode_elmt_ic(:,i_elmt)
  !  ngradphi_ic(:,num(:)) = gradphi_ic(:,:,i_elmt)
  !enddo
  !! Outer core
  !ngradphi_oc(:,:) = zero
  !do i_elmt = 1,NSPEC_OUTER_CORE
  !  num(:) = inode_elmt_oc(:,i_elmt)
  !  ngradphi_oc(:,num(:)) = gradphi_oc(:,:,i_elmt)
  !enddo
  !! Crust mantle
  !ngradphi_cm(:,:) = zero
  !do i_elmt = 1,NSPEC_CRUST_MANTLE
  !  num(:) = inode_elmt_cm(:,i_elmt)
  !  ngradphi_cm(:,num(:)) = gradphi_cm(:,:,i_elmt)
  !enddo
  !
  ! compute global gravity (magnitude?) - scalar g?
  !grav0_ic(:) = zero
  !do i_node = 1,NGLOB_INNER_CORE
  !  grav0_ic(i_node) = sqrt(sum(ngradphi_ic(:,i_node)*ngradphi_ic(:,i_node)))
  !enddo
  !grav0_oc(:) = zero
  !do i_node = 1,NGLOB_OUTER_CORE
  !  grav0_oc(i_node) = sqrt(sum(ngradphi_oc(:,i_node)*ngradphi_oc(:,i_node)))
  !enddo
  !grav0_cm(:) = zero
  !do i_node = 1,NGLOB_CRUST_MANTLE
  !  grav0_cm(i_node) = sqrt(sum(ngradphi_cm(:,i_node)*ngradphi_cm(:,i_node)))
  !enddo
  !deallocate(ngradphi_cm,ngradphi_oc,ngradphi_ic)
  !
  !destag = 'unstructured meshes'
  !npart  = 1
  ! Set parameters for saving to Ensight
  !write(ptail,*) myrank; ptail = adjustl(ptail)
  !
  ! (De)Active ids for inner core
  !id_ic(:) = 1.0_CUSTOM_REAL ! active
  !do i_elmt = 1,NSPEC_INNER_CORE
  !  if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) then
  !    id_ic(i_elmt) = 0.0_CUSTOM_REAL ! deactive
  !  endif
  !enddo
  !
  ! Write arrays to Ensight - potential
  ! crust mantle
  !write(file_head,'(a,i1,a)')'reg',1,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.gpot'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(pgrav_cm(inode4_cm)))
  ! outer core
  !write(file_head,'(a,i1,a)')'reg',2,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.gpot'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_oc,real(pgrav_oc(inode4_oc)))
  ! inner core
  !write(file_head,'(a,i1,a)')'reg',3,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.gpot'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_ic,real(pgrav_ic(inode4_ic)))
  !
  ! Write arrays to Ensight - gravity
  ! crust mantle
  !write(file_head,'(a,i1,a)')'reg',1,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(grav0_cm(inode4_cm)))
  ! outer core
  !write(file_head,'(a,i1,a)')'reg',2,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_oc,real(grav0_oc(inode4_oc)))
  ! inner core to .id
  !write(file_head,'(a,i1,a)')'reg',3,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.id'
  !call write_ensight_perelementAS(grav_file,ensight_etype,destag,npart,real(id_ic))
  ! inner core to .grav?
  !write(file_head,'(a,i1,a)')'reg',3,'_proc'
  !grav_file = trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_ic,real(grav0_ic(inode4_ic)))

  ! map to gdof array
  pgrav(:) = zero

  ! Assemble global pgrav array values from regional
  pgrav(gdof_ic(:)) = pgrav_ic(:)
  pgrav(gdof_oc(:)) = pgrav_oc(:)
  pgrav(gdof_cm(:)) = pgrav_cm(:)
  if (ADD_TRINF) pgrav(gdof_trinf(:)) = pgrav_trinf(:)
  pgrav(gdof_inf(:)) = pgrav_inf(:)

  ! Level-2 solver
  if (USE_POISSON_SOLVER_5GLL) then
    !debug
    if (myrank == 0) print *,' Level 2 solver starts...'

    ! Using solvers/rhoload etc for 5 GLL points

    ! Calculate the RHS (gravload)
    call compute_poisson_rhoload()

    if (POISSON_SOLVER == ISOLVER_BUILTIN) then
      ! built-in solver
      if (CG_SCALING) then
        gravload(:) = ndscale(:) * gravload(:)
        pgrav(1:) = pgrav(1:) / ndscale(1:) ! use scaling

        call cg_solver(myrank,neq,pgrav,gravload,cg_iter)
        if (myrank == 0) print *,'CG solver: L2: solve done'

        pgrav(:) = ndscale(:) * pgrav(:)
      else
        call diagpcg_solver(myrank,neq,pgrav,gravload,dprecon,cg_iter)
        if (myrank == 0) print *,'PCG solver: L2: solve done'
      endif
    else
      ! petsc solver
      call petsc_set_vector()
      if (myrank == 0) print *,'PETSc solver: L2: Set RHS with PETSC'

      call petsc_solve(pgrav(1:),cg_iter)
      if (myrank == 0) print *,'PETSc solver: L2: Solved with PETSC'
    endif
    call synchronize_all()

    !debug
    if (myrank == 0) print *,'---- Level 2 solver complete! ','iterations: ',cg_iter
    if (myrank == 0) print *

    ! transfer to array for each region (cm, ic, oc)
    ! note in Level-1 uses pgrav_ic1, here uses pgrav_ic
    pgrav_ic(:) = pgrav(gdof_ic(:))
    pgrav_oc(:) = pgrav(gdof_oc(:))
    pgrav_cm(:) = pgrav(gdof_cm(:))
    if (ADD_TRINF) pgrav_trinf(:) = pgrav(gdof_trinf(:))
    pgrav_inf(:) = pgrav(gdof_inf(:))
  endif !(SOLVER_5GLL) then

  ! ensight vis
  ! compute gradphi only for fluid/solid regions
  !call poisson_gravity(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
  !                     ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                     pgrav_cm,gradphi_cm)
  !
  !call poisson_gravity(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
  !                     ibool_outer_core,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                     pgrav_oc,gradphi_oc)
  !
  !call poisson_gravity(IREGION_INNER_CORE,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
  !                     ibool_inner_core,xstore_inner_core,ystore_inner_core,zstore_inner_core, &
  !                     pgrav_ic,gradphi_ic)
  !
  ! compute global gravity
  !
  !! valency
  !allocate(nvalency_ic(NGLOB_INNER_CORE), &
  !         nvalency_oc(NGLOB_OUTER_CORE), &
  !         nvalency_cm(NGLOB_CRUST_MANTLE),stat=ier)
  !if (ier /= 0) stop 'Error allocating nvalency_ic,.. arrays'
  !nvalency_ic(:) = 0; nvalency_oc(:) = 0; nvalency_cm(:) = 0
  !! compute node valency only once
  !call compute_full_gravity_node_valency(nvalency_ic, nvalency_oc, nvalency_cm)
  !
  !! Inner core
  !grav0_ic = zero
  !do i_elmt = 1,NSPEC_INNER_CORE
  !  if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
  !  num(:) = inode_elmt_ic(:,i_elmt)
  !  grav0_ic(num(:)) = grav0_ic(num(:)) + &
  !                     sqrt(gradphi_ic(1,:,i_elmt)*gradphi_ic(1,:,i_elmt)+ &
  !                          gradphi_ic(2,:,i_elmt)*gradphi_ic(2,:,i_elmt)+gradphi_ic(3,:,i_elmt)*gradphi_ic(3,:,i_elmt))
  !enddo
  !! compute average at sharing nodes
  !do i_node = 1,NGLOB_INNER_CORE
  !  grav0_ic(i_node) = grav0_ic(i_node)/real(nvalency_ic(i_node),CUSTOM_REAL)
  !enddo
  !! Outer core
  !grav0_oc = zero
  !do i_elmt = 1,NSPEC_OUTER_CORE
  !  num(:) = inode_elmt_oc(:,i_elmt)
  !  grav0_oc(num(:)) = grav0_oc(num(:)) + &
  !                     sqrt(gradphi_oc(1,:,i_elmt)*gradphi_oc(1,:,i_elmt)+ &
  !                          gradphi_oc(2,:,i_elmt)*gradphi_oc(2,:,i_elmt)+gradphi_oc(3,:,i_elmt)*gradphi_oc(3,:,i_elmt))
  !enddo
  !! compute average at sharing nodes
  !do i_node = 1,NGLOB_OUTER_CORE
  !  grav0_oc(i_node) = grav0_oc(i_node)/real(nvalency_oc(i_node),CUSTOM_REAL)
  !enddo
  !! Crust mantle
  !grav0_cm = zero
  !do i_elmt = 1,NSPEC_CRUST_MANTLE
  !  num(:) = inode_elmt_cm(:,i_elmt)
  !  grav0_cm(num(:)) = grav0_cm(num(:)) + &
  !                     sqrt(gradphi_cm(1,:,i_elmt)*gradphi_cm(1,:,i_elmt)+ &
  !                          gradphi_cm(2,:,i_elmt)*gradphi_cm(2,:,i_elmt)+gradphi_cm(3,:,i_elmt)*gradphi_cm(3,:,i_elmt))
  !enddo
  !! compute average at sharing nodes
  !do i_node = 1,NGLOB_INNER_CORE
  !  grav0_ic(i_node) = grav0_ic(i_node)/real(nvalency_ic(i_node),CUSTOM_REAL)
  !enddo
  !
  !deallocate(nvalency_cm,nvalency_ic,nvalency_oc)
  !
  ! Write to Ensight - note here that L2 solver doesnt calculate scalar grav0
  !write(file_head,'(a,i1,a)')'reg',1,'_proc'
  !grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav5'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(grav0_cm(inode4_cm)))
  !
  !write(file_head,'(a,i1,a)')'reg',2,'_proc'
  !grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav5'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_oc,real(grav0_oc(inode4_oc)))
  !
  !write(file_head,'(a,i1,a)')'reg',3,'_proc'
  !grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav5'
  !call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_ic,real(grav0_ic(inode4_ic)))
  !
  ! free temporary arrays
  !deallocate(gradphi_cm,gradphi_oc,gradphi_ic)

  call synchronize_all()

  ! user output
  if (myrank == 0) then
    tCPU = wtime() - tstart
    write(IMAIN,*) '    Elapsed time for background gravity solve: ',sngl(tCPU),'(s)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine SIEM_compute_background_gravity

!
!-------------------------------------------------------------------------------
!

! not used yet...

!  subroutine create_full_gravity_mapping()
!
!  ! this routine determines nnode4_* and inode4_* arrays for ensight visualizations
!
!  use constants, only: myrank, NGLLX, NGLLY, NGLLZ, IMAIN
!  use specfem_par, only: NGLOB_INNER_CORE, NGLOB_OUTER_CORE, NGLOB_CRUST_MANTLE, &
!                         NSPEC_INNER_CORE, NSPEC_OUTER_CORE, NSPEC_CRUST_MANTLE
!
!  use specfem_par_full_gravity, only: inode_elmt_cm, inode_elmt_ic, inode_elmt_oc
!
!  ! ensight
!  use specfem_par_full_gravity, only: nnode4_ic, nnode4_oc, nnode4_cm, inode4_ic, inode4_oc, inode4_cm
!
!  implicit none
!
!  ! Local variables:
!  integer :: dnx, dny, dnz, inum, i, j, k, igll4(8), i_elmt, i_node
!  logical,dimension(:),allocatable :: isnode_ic,isnode_oc,isnode_cm
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) '  creating full gravity mapping'
!    call flush_IMAIN()
!  endif
!
!  allocate(isnode_ic(NGLOB_INNER_CORE),isnode_oc(NGLOB_OUTER_CORE),isnode_cm(NGLOB_CRUST_MANTLE))
!
!  dnx = NGLLX-1; dny = NGLLY-1; dnz = NGLLZ-1
!
!  inum = 0
!  do k = 1,NGLLZ,dnz
!    do j = 1,NGLLY,dny
!      do i = 1,NGLLX,dnx
!        inum = inum+1
!        igll4(inum) = NGLLY*NGLLX*(k-1)+NGLLX*(j-1)+i
!      enddo
!    enddo
!  enddo
!
!  ! WE - setting bool for corner elements in inner core?
!  isnode_ic(:) = .false.
!  do i_elmt = 1,NSPEC_INNER_CORE
!    isnode_ic(inode_elmt_ic(igll4,i_elmt)) = .true.
!  enddo
!  ! WE - number of corner nodes in IC?
!  nnode4_ic = count(isnode_ic)
!  allocate(inode4_ic(nnode4_ic))
!  inum = 0
!  ! WE - array of indices to corner nodes in IC?
!  do i_node = 1,NGLOB_INNER_CORE
!    if (isnode_ic(i_node)) then
!      inum = inum+1
!      inode4_ic(inum) = i_node
!    endif
!  enddo
!
!  ! WE - setting bool for corner elements in outer core?
!  isnode_oc(:) = .false.
!  do i_elmt = 1,NSPEC_OUTER_CORE
!    isnode_oc(inode_elmt_oc(igll4,i_elmt)) = .true.
!  enddo
!  ! WE - number of corner nodes in OC?
!  nnode4_oc = count(isnode_oc)
!  allocate(inode4_oc(nnode4_oc))
!  inum = 0
!  ! WE - array of indices to corner nodes in OC?
!  do i_node = 1,NGLOB_OUTER_CORE
!    if (isnode_oc(i_node)) then
!      inum = inum+1
!      inode4_oc(inum) = i_node
!    endif
!  enddo
!
!  ! WE - Same for Crust mantle?
!  isnode_cm(:) = .false.
!  do i_elmt = 1,NSPEC_CRUST_MANTLE
!    isnode_cm(inode_elmt_cm(igll4,i_elmt)) = .true.
!  enddo
!  nnode4_cm = count(isnode_cm)
!  allocate(inode4_cm(nnode4_cm))
!  inum = 0
!  do i_node = 1,NGLOB_CRUST_MANTLE
!    if (isnode_cm(i_node)) then
!      inum = inum+1
!      inode4_cm(inum) = i_node
!    endif
!  enddo
!
!  ! free temporary arrays
!  deallocate(isnode_cm,isnode_ic,isnode_oc)
!
!  end subroutine create_full_gravity_mapping

!
!-------------------------------------------------------------------------------
!

! not used yet...

!  subroutine compute_full_gravity_node_valency(nvalency_ic, nvalency_oc, nvalency_cm)
!
!  use constants, only: myrank,IMAIN,NGLLCUBE, IFLAG_IN_FICTITIOUS_CUBE
!
!  use specfem_par, only: NSPEC_INNER_CORE, NGLOB_CRUST_MANTLE, NGLOB_OUTER_CORE, NGLOB_INNER_CORE, &
!                         NSPEC_OUTER_CORE, NSPEC_CRUST_MANTLE
!  use specfem_par_innercore, only: idoubling_inner_core
!
!  use specfem_par_full_gravity, only: inode_elmt_ic,inode_elmt_oc,inode_elmt_cm
!
!  implicit none
!  ! IO variables
!  integer,intent(inout) :: nvalency_ic(NGLOB_INNER_CORE),nvalency_oc(NGLOB_OUTER_CORE), nvalency_cm(NGLOB_CRUST_MANTLE)
!  ! local parameters
!  integer :: num(NGLLCUBE)
!  integer :: i_elmt
!
!  ! user output
!  if (myrank == 0) then
!    write(IMAIN,*) '    computing node valency'
!    call flush_IMAIN()
!  endif
!
!  nvalency_ic(:) = 0
!  do i_elmt = 1,NSPEC_INNER_CORE
!    ! do not consider fictitious elements
!    if (idoubling_inner_core(i_elmt) == IFLAG_IN_FICTITIOUS_CUBE) cycle
!    num(:) = inode_elmt_ic(:,i_elmt)
!    nvalency_ic(num(:)) = nvalency_ic(num(:))+1
!  enddo
!
!  nvalency_oc(:) = 0
!  do i_elmt = 1,NSPEC_OUTER_CORE
!    num(:) = inode_elmt_oc(:,i_elmt)
!    nvalency_oc(num(:)) = nvalency_oc(num(:))+1
!  enddo
!
!  nvalency_cm(:) = 0
!  do i_elmt = 1,NSPEC_CRUST_MANTLE
!    num(:) = inode_elmt_cm(:,i_elmt)
!    nvalency_cm(num(:)) = nvalency_cm(num(:))+1
!  enddo
!
!  ! Make sure nvalency is NOT equal to 0 to avoid the 0 division!
!  where(nvalency_ic == 0) nvalency_ic = 1
!  where(nvalency_oc == 0) nvalency_oc = 1
!  where(nvalency_cm == 0) nvalency_cm = 1
!
!  end subroutine compute_full_gravity_node_valency

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_finalize()

  use specfem_par
  use specfem_par_full_gravity

  use siem_solver_petsc, only: petsc_finalize1,petsc_finalize

  implicit none

  ! safety check
  if (.not. FULL_GRAVITY) return

  if (POISSON_SOLVER == ISOLVER_PETSC) then
    ! petsc solver
    call petsc_finalize1()
    if (USE_POISSON_SOLVER_5GLL) then
      ! Level-2 solver
      call petsc_finalize()
    endif
  endif

  ! free arrays
  ! Level-1 solver
  if (allocated(storederiv_cm1)) deallocate(storederiv_cm1,storerhojw_cm1)
  if (allocated(storederiv_ic1)) deallocate(storederiv_ic1,storerhojw_ic1)
  if (allocated(storederiv_oc1)) deallocate(storederiv_oc1,storerhojw_oc1)
  ! Level-2 solver
  if (allocated(storederiv_cm)) deallocate(storederiv_cm,storerhojw_cm)
  if (allocated(storederiv_ic)) deallocate(storederiv_ic,storerhojw_ic)
  if (allocated(storederiv_oc)) deallocate(storederiv_oc,storerhojw_oc)

  if (allocated(inode_elmt_cm1)) deallocate(inode_elmt_cm1,inode_elmt_ic1,inode_elmt_oc1, &
                                          inode_elmt_trinf1,inode_elmt_inf1)
  if (allocated(inode_elmt_cm)) deallocate(inode_elmt_cm,inode_elmt_ic,inode_elmt_oc, &
                                          inode_elmt_trinf,inode_elmt_inf)
  if (allocated(inode_map_ic)) deallocate(inode_map_ic,inode_map_oc,inode_map_cm, &
                                         inode_map_trinf,inode_map_inf)

  ! builtin solver scaling
  if (allocated(ndscale1)) deallocate(ndscale1)
  if (allocated(ndscale)) deallocate(ndscale)

  ! gravity perturbation arrays
  ! Level-1 solver
  if (allocated(pgrav1)) deallocate(pgrav1)
  if (allocated(pgrav_ic1)) deallocate(pgrav_ic1,pgrav_oc1,pgrav_cm1,pgrav_trinf1,pgrav_inf1)
  if (allocated(dprecon1)) deallocate(dprecon1,gravload1)
  if (allocated(b_gravload1)) deallocate(b_gravload1)
  ! Level-2 solver
  if (allocated(pgrav)) deallocate(pgrav)
  if (allocated(pgrav_ic)) deallocate(pgrav_ic,pgrav_oc,pgrav_cm,pgrav_trinf,pgrav_inf)
  if (allocated(dprecon)) deallocate(dprecon,gravload)
  ! sparse petsc solver arrays
  ! Level-1 solver
  if (allocated(l2gdof1)) deallocate(l2gdof1)
  if (allocated(krow_sparse1)) deallocate(krow_sparse1,kcol_sparse1)
  if (allocated(kgrow_sparse1)) deallocate(kgrow_sparse1,kgcol_sparse1)
  ! Level-2 solver
  if (allocated(l2gdof)) deallocate(l2gdof)
  if (allocated(krow_sparse)) deallocate(krow_sparse,kcol_sparse)
  if (allocated(kgrow_sparse)) deallocate(kgrow_sparse,kgcol_sparse)

  ! seismos
  if (allocated(seismograms_phi)) deallocate(seismograms_phi)
  if (allocated(seismograms_pgrav)) deallocate(seismograms_pgrav)
  if (allocated(seismograms_Hgrav)) deallocate(seismograms_Hgrav)
  if (allocated(seismograms_corio)) deallocate(seismograms_corio)
  if (allocated(seismograms_grav)) deallocate(seismograms_grav)
  if (allocated(g_spec_rec)) deallocate(g_spec_rec)
  if (allocated(gradg_rec)) deallocate(gradg_rec)
  if (allocated(storederiv_rec)) deallocate(storederiv_rec)

  end subroutine SIEM_finalize
