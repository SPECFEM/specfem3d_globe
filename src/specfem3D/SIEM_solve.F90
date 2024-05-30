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


  subroutine SIEM_solve_poisson()

  use specfem_par, only: SIMULATION_TYPE

  implicit none

  ! forward wavefield
  call solve_poisson_equation()

  ! backward/reconstructed wavefield
  if (SIMULATION_TYPE == 3) then
    call solve_poisson_equation_backward()
  endif

  end subroutine SIEM_solve_poisson

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_solve_poisson_backward()

  implicit none

  ! backward/reconstructed wavefield
  call solve_poisson_equation_backward()

  end subroutine SIEM_solve_poisson_backward

!
!-------------------------------------------------------------------------------
!

  subroutine solve_poisson_equation()

  use specfem_par

  use specfem_par_full_gravity, only: ndscale1,ndscale, &
    gdof_cm, gdof_cm1, inode_elmt_cm, inode_map_cm, nmir_cm, &
    gdof_ic, gdof_ic1, inode_elmt_ic, inode_map_ic, nmir_ic, &
    gdof_oc, gdof_oc1, inode_elmt_oc, inode_map_oc, nmir_oc, &
    gdof_trinf, gdof_trinf1, inode_elmt_trinf, inode_map_trinf, nmir_trinf, &
    gdof_inf, gdof_inf1, inode_elmt_inf, inode_map_inf, nmir_inf, &
    nnode_cm1,nnode_ic1,nnode_oc1,nnode_trinf1,nnode_inf1, &
    neq1,neq, &
    is_active_gll,igll_active_on, &
    CG_SCALING

  use specfem_par_full_gravity, only: gravload1,gravload, &
    pgrav1,pgrav_ic1,pgrav_oc1,pgrav_cm1,pgrav_trinf1,pgrav_inf1, &
    pgrav,pgrav_ic,pgrav_oc,pgrav_cm,pgrav_trinf,pgrav_inf

  use siem_math_library_mpi, only: maxvec

  use siem_solver_petsc, only: petsc_set_vector1,petsc_solve1, &
                               petsc_set_vector,petsc_solve

  use siem_poisson, only: compute_poisson_load,compute_poisson_load3 !,poisson_gravity

  use siem_solver_mpi, only: cg_solver,cg_solver3,diagpcg_solver,diagpcg_solver3, interpolate3to5

  implicit none

  ! Local variables
  real(kind=CUSTOM_REAL) :: maxf,upscale
  real(kind=CUSTOM_REAL),parameter :: rone = 1.0_CUSTOM_REAL

  ! load for Level-1 solver
  ! Compute the RHS (gravload1)
  call compute_poisson_load3()

  !debug
  !maxf = maxvec(abs(gravload1))
  !print *,'debug solve: L1: rank',myrank,' load:',minval(gravload1),maxval(gravload1),' max:',maxf

  ! Level-1 solver
  if (POISSON_SOLVER == ISOLVER_BUILTIN) then
    ! builtin solver
    maxf = maxvec(abs(gravload1))
    upscale = rone
    if (maxf > zero) upscale = rone / maxf

    if (CG_SCALING) then
      gravload1(:) = gravload1(:) * ndscale1(:) * upscale
      pgrav1(1:) = upscale * pgrav1(1:) / ndscale1(1:)
    endif

    call cg_solver3(myrank,neq1,pgrav1,gravload1)

    if (CG_SCALING) then
      pgrav1(:) = ndscale1(:) * pgrav1(:) / upscale
    endif
  else
    ! PETSc solver
    call petsc_set_vector1(gravload1)
    call petsc_solve1(pgrav1(1:))
  endif

  ! interpolate GLL3 results to GLL5
  pgrav_ic1(:) = pgrav1(gdof_ic1(:))
  pgrav_oc1(:) = pgrav1(gdof_oc1(:))
  pgrav_cm1(:) = pgrav1(gdof_cm1(:))
  if (ADD_TRINF) pgrav_trinf1(:) = pgrav1(gdof_trinf1(:))
  pgrav_inf1(:) = pgrav1(gdof_inf1(:))

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

  ! map to gdof array
  pgrav(:) = zero

  ! assemble across the regions
  pgrav(gdof_ic(:)) = pgrav_ic(:)
  pgrav(gdof_oc(:)) = pgrav_oc(:)
  pgrav(gdof_cm(:)) = pgrav_cm(:)
  if (ADD_TRINF) pgrav(gdof_trinf(:)) = pgrav_trinf(:)
  pgrav(gdof_inf(:)) = pgrav_inf(:)

  ! Level-2 solver
  if (USE_POISSON_SOLVER_5GLL) then
    ! Compute the RHS (gravload)
    call compute_poisson_load()

    if (POISSON_SOLVER == ISOLVER_BUILTIN) then
      ! builtin solver
      if (CG_SCALING) then
        gravload(:) = ndscale(:) * gravload(:)
        pgrav(1:) = pgrav(1:) / ndscale(1:)
      endif

      call cg_solver(myrank,neq,pgrav,gravload)

      if (CG_SCALING) then
        pgrav(:) = ndscale(:) * pgrav(:)
      endif
    else
      ! PETSc solver
      call petsc_set_vector()
      call petsc_solve(pgrav(1:))
    endif

    ! set to regional array
    pgrav_ic(:) = pgrav(gdof_ic(:))
    pgrav_oc(:) = pgrav(gdof_oc(:))
    pgrav_cm(:) = pgrav(gdof_cm(:))
    if (ADD_TRINF) pgrav_trinf(:) = pgrav(gdof_trinf(:))
    pgrav_inf(:) = pgrav(gdof_inf(:))
  endif

  ! TODO: \grad\phi computed here
  !
  ! compute gradphi only for solid regions
  ! call poisson_gravity(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
  !                      ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                      pgrav_cm,gradphi_cm)
  !
  ! call poisson_gravity(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
  !                      ibool_outer_core,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                      pgrav_oc,gradphi_oc)
  !
  ! call poisson_gravity(IREGION_INNER_CORE,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
  !                      ibool_inner_core,xstore_inner_core,ystore_inner_core,zstore_inner_core, &
  !                      pgrav_ic,gradphi_ic)

  end subroutine solve_poisson_equation

!
!-------------------------------------------------------------------------------
!

  subroutine solve_poisson_equation_backward()

  use specfem_par

  use specfem_par_full_gravity, only: &
    gdof_cm, gdof_cm1, inode_elmt_cm, inode_map_cm, nmir_cm, &
    gdof_ic, gdof_ic1, inode_elmt_ic, inode_map_ic, nmir_ic, &
    gdof_oc, gdof_oc1, inode_elmt_oc, inode_map_oc, nmir_oc, &
    gdof_trinf, gdof_trinf1, inode_elmt_trinf, inode_map_trinf, nmir_trinf, &
    gdof_inf, gdof_inf1, inode_elmt_inf, inode_map_inf, nmir_inf, &
    nnode_cm1,nnode_ic1,nnode_oc1,nnode_trinf1,nnode_inf1, &
    is_active_gll,igll_active_on, &
    CG_SCALING

  use specfem_par_full_gravity, only: b_gravload1, &
    b_pgrav1,b_pgrav_ic1,b_pgrav_oc1,b_pgrav_cm1,b_pgrav_trinf1,b_pgrav_inf1, &
    b_pgrav,b_pgrav_ic,b_pgrav_oc,b_pgrav_cm,b_pgrav_trinf,b_pgrav_inf

  use siem_solver_petsc, only: petsc_set_vector1,petsc_solve1
  use siem_poisson, only: compute_backward_poisson_load3 !,poisson_gravity
  use siem_solver_mpi, only: cg_solver,cg_solver3,diagpcg_solver,diagpcg_solver3, interpolate3to5

  implicit none

  ! safety check
  if (POISSON_SOLVER == ISOLVER_BUILTIN) then
    print *,'Error. Adjoint FG only implemented for petsc solver'
    stop 'FULL_GRAVITY adjoint FG not fully implemented yet of builtin solver'
  endif

  ! load for Level-1 solver

  ! Compute the RHS (b_gravload1)
  call compute_backward_poisson_load3()

  if (POISSON_SOLVER == ISOLVER_PETSC) then
    ! PETSc solver
    call petsc_set_vector1(b_gravload1)
    call petsc_solve1(b_pgrav1(1:))
  endif

  ! interpolate GLL3 results to GLL5
  b_pgrav_ic1(:) = b_pgrav1(gdof_ic1(:))
  b_pgrav_oc1(:) = b_pgrav1(gdof_oc1(:))
  b_pgrav_cm1(:) = b_pgrav1(gdof_cm1(:))
  if (ADD_TRINF) b_pgrav_trinf1(:) = b_pgrav1(gdof_trinf1(:))
  b_pgrav_inf1(:) = b_pgrav1(gdof_inf1(:))

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
                       inode_elmt_inf,nmir_inf,inode_map_inf,is_active_gll,igll_active_on,b_pgrav_inf1, b_pgrav_inf)

  ! map to gdof array
  b_pgrav(:) = zero

  ! assemble across the regions
  b_pgrav(gdof_ic(:)) = b_pgrav_ic(:)
  b_pgrav(gdof_oc(:)) = b_pgrav_oc(:)
  b_pgrav(gdof_cm(:)) = b_pgrav_cm(:)
  if (ADD_TRINF) b_pgrav(gdof_trinf(:)) = b_pgrav_trinf(:)
  b_pgrav(gdof_inf(:)) = b_pgrav_inf(:)

  ! TODO: \grad\phi computed here
  !
  ! compute gradphi only for solid regions
  ! call poisson_gravity(IREGION_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
  !                      ibool_crust_mantle,xstore_crust_mantle,ystore_crust_mantle, &
  !                      zstore_crust_mantle,pgrav_cm,gradphi_cm)
  !
  ! call poisson_gravity(IREGION_OUTER_CORE,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
  !                      ibool_outer_core,xstore_outer_core,ystore_outer_core, &
  !                      zstore_outer_core,pgrav_oc,gradphi_oc)
  !
  ! call poisson_gravity(IREGION_INNER_CORE,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
  !                      ibool_inner_core,xstore_inner_core,ystore_inner_core, &
  !                      zstore_inner_core,pgrav_ic,gradphi_ic)

  end subroutine solve_poisson_equation_backward

