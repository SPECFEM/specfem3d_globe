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


  subroutine compute_coupling_fluid()

  use specfem_par
  use specfem_par_crustmantle, only: displ_crust_mantle,ibool_crust_mantle,ibelm_bottom_crust_mantle
  use specfem_par_innercore, only: displ_inner_core,ibool_inner_core,ibelm_top_inner_core
  use specfem_par_outercore

  implicit none

  ! checks if anything to do
  if (.not. ACTUALLY_COUPLE_FLUID_CMB .and. .not. ACTUALLY_COUPLE_FLUID_ICB) return

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

  ! only for elements in first matching layer in the fluid
  if (.not. GPU_MODE) then
    ! on CPU
    if (ACTUALLY_COUPLE_FLUID_CMB .and. ACTUALLY_COUPLE_FLUID_ICB) then
      ! couple both mantle and inner core
      call compute_coupling_fluid_CMB_ICB(NGLOB_OUTER_CORE,accel_outer_core,wgllwgll_xy,ibool_outer_core, &
                                          NGLOB_CRUST_MANTLE,displ_crust_mantle, &
                                          ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                          NGLOB_INNER_CORE,displ_inner_core, &
                                          ibool_inner_core,ibelm_top_inner_core, &
                                          normal_top_outer_core,jacobian2D_top_outer_core, &
                                          normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                          ibelm_top_outer_core,NSPEC2D_TOP(IREGION_OUTER_CORE), &
                                          ibelm_bottom_outer_core,NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
    else
      ! separate only one coupling interface
      !---
      !--- couple with mantle at the top of the outer core
      !---
      if (ACTUALLY_COUPLE_FLUID_CMB) &
        call compute_coupling_fluid_CMB(NGLOB_CRUST_MANTLE,displ_crust_mantle, &
                                           ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                           NGLOB_OUTER_CORE,accel_outer_core, &
                                           normal_top_outer_core,jacobian2D_top_outer_core, &
                                           wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                           NSPEC2D_TOP(IREGION_OUTER_CORE))
      !---
      !--- couple with inner core at the bottom of the outer core
      !---
      if (ACTUALLY_COUPLE_FLUID_ICB) &
        call compute_coupling_fluid_ICB(NGLOB_INNER_CORE,displ_inner_core, &
                                           ibool_inner_core,ibelm_top_inner_core, &
                                           NGLOB_OUTER_CORE,accel_outer_core, &
                                           normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                           wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                           NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
    endif

  else
    ! on GPU
    !---
    !--- couple with mantle at the top of the outer core
    !---
    if (ACTUALLY_COUPLE_FLUID_CMB ) &
         call compute_coupling_fluid_cmb_gpu(Mesh_pointer,1)
    !---
    !--- couple with inner core at the bottom of the outer core
    !---
    if (ACTUALLY_COUPLE_FLUID_ICB ) &
         call compute_coupling_fluid_icb_gpu(Mesh_pointer,1)

  endif

  end subroutine compute_coupling_fluid

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_fluid_backward()

  use specfem_par
  use specfem_par_crustmantle, only: b_displ_crust_mantle,ibool_crust_mantle,ibelm_bottom_crust_mantle
  use specfem_par_innercore, only: b_displ_inner_core,ibool_inner_core,ibelm_top_inner_core
  use specfem_par_outercore

  implicit none

  ! checks if anything to do
  if (.not. ACTUALLY_COUPLE_FLUID_CMB .and. .not. ACTUALLY_COUPLE_FLUID_ICB) return

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

  ! only for elements in first matching layer in the fluid
  if (.not. GPU_MODE) then
    ! on CPU
    if (ACTUALLY_COUPLE_FLUID_CMB .and. ACTUALLY_COUPLE_FLUID_ICB) then
      ! couple both mantle and inner core
      call compute_coupling_fluid_CMB_ICB(NGLOB_OUTER_CORE,b_accel_outer_core,wgllwgll_xy,ibool_outer_core, &
                                          NGLOB_CRUST_MANTLE,b_displ_crust_mantle, &
                                          ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                          NGLOB_INNER_CORE,b_displ_inner_core, &
                                          ibool_inner_core,ibelm_top_inner_core, &
                                          normal_top_outer_core,jacobian2D_top_outer_core, &
                                          normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                          ibelm_top_outer_core,NSPEC2D_TOP(IREGION_OUTER_CORE), &
                                          ibelm_bottom_outer_core,NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
    else
      ! separate only one coupling interface
      !---
      !--- couple with mantle at the top of the outer core
      !---
      if (ACTUALLY_COUPLE_FLUID_CMB) &
        call compute_coupling_fluid_CMB(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle, &
                                        ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                        NGLOB_OUTER_CORE_ADJOINT,b_accel_outer_core, &
                                        normal_top_outer_core,jacobian2D_top_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                        NSPEC2D_TOP(IREGION_OUTER_CORE))

      !---
      !--- couple with inner core at the bottom of the outer core
      !---
      if (ACTUALLY_COUPLE_FLUID_ICB) &
        call compute_coupling_fluid_ICB(NGLOB_INNER_CORE_ADJOINT,b_displ_inner_core, &
                                        ibool_inner_core,ibelm_top_inner_core, &
                                        NGLOB_OUTER_CORE_ADJOINT,b_accel_outer_core, &
                                        normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                        NSPEC2D_BOTTOM(IREGION_OUTER_CORE))
    endif
  else
    ! on GPU
    !---
    !--- couple with mantle at the top of the outer core
    !---
    if (ACTUALLY_COUPLE_FLUID_CMB ) &
      call compute_coupling_fluid_cmb_gpu(Mesh_pointer,3)
    !---
    !--- couple with inner core at the bottom of the outer core
    !---
    if (ACTUALLY_COUPLE_FLUID_ICB ) &
      call compute_coupling_fluid_icb_gpu(Mesh_pointer,3)

  endif

  end subroutine compute_coupling_fluid_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_fluid_CMB_ICB(NGLOB_OC,accel_outer_core,wgllwgll_xy,ibool_outer_core, &
                                            NGLOB_CM,displ_crust_mantle, &
                                            ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                            NGLOB_IC,displ_inner_core, &
                                            ibool_inner_core,ibelm_top_inner_core, &
                                            normal_top_outer_core,jacobian2D_top_outer_core, &
                                            normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                            ibelm_top_outer_core,nspec2D_top, &
                                            ibelm_bottom_outer_core,nspec_bottom)

! this routine combines both coupling with mantle and inner core
! this will optimize OpenMP performance

  use constants_solver

  implicit none

  integer,intent(in) :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC),intent(inout) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY),intent(in) :: wgllwgll_xy
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),intent(in) :: ibool_outer_core

  integer,intent(in) :: NGLOB_CM
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM),intent(in) :: displ_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(in) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM),intent(in) :: ibelm_bottom_crust_mantle

  integer,intent(in) :: NGLOB_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC),intent(in) :: displ_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),intent(in) :: ibool_inner_core
  integer, dimension(NSPEC2D_TOP_IC),intent(in) :: ibelm_top_inner_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC),intent(in) :: normal_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC),intent(in) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC),intent(in) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC),intent(in) :: jacobian2D_bottom_outer_core

  integer, dimension(NSPEC2D_TOP_OC),intent(in) :: ibelm_top_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC),intent(in) :: ibelm_bottom_outer_core
  integer,intent(in) :: nspec2D_top
  integer,intent(in) :: nspec_bottom

  ! local parameters
  real(kind=CUSTOM_REAL) :: displ_x,displ_y,displ_z,displ_n,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob_cm,iglob_oc,iglob_ic,ispec_selected

  ! for surface elements exactly on the CMB

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

! openmp solver
!$OMP PARALLEL if (nspec2D_top > 500) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,ispec_selected,i,j,k,k_corresp,iglob_cm,iglob_ic,iglob_oc, &
!$OMP displ_x,displ_y,displ_z,nx,ny,nz,displ_n,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy)
!$OMP DO
  do ispec2D = 1,nspec2D_top !NSPEC2D_TOP(IREGION_OUTER_CORE)

    ispec = ibelm_top_outer_core(ispec2D)
    ispec_selected = ibelm_bottom_crust_mantle(ispec2D)

    ! only for DOFs exactly on the CMB (top of these elements)
    k = NGLLZ

    ! get displacement on the solid side using pointwise matching
    k_corresp = 1

    do j = 1,NGLLY
      do i = 1,NGLLX
        ! corresponding points are located at the bottom of the mantle
        iglob_cm = ibool_crust_mantle(i,j,k_corresp,ispec_selected)

        displ_x = displ_crust_mantle(1,iglob_cm)
        displ_y = displ_crust_mantle(2,iglob_cm)
        displ_z = displ_crust_mantle(3,iglob_cm)

        ! get normal on the CMB
        nx = normal_top_outer_core(1,i,j,ispec2D)
        ny = normal_top_outer_core(2,i,j,ispec2D)
        nz = normal_top_outer_core(3,i,j,ispec2D)

        ! compute dot product
        displ_n = displ_x*nx + displ_y*ny + displ_z*nz

        ! formulation with generalized potential
        weight = jacobian2D_top_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        ! get global point number
        iglob_oc = ibool_outer_core(i,j,k,ispec)

        ! update fluid acceleration/pressure
!$OMP ATOMIC
        accel_outer_core(iglob_oc) = accel_outer_core(iglob_oc) + weight*displ_n
      enddo
    enddo
  enddo
!$OMP ENDDO NOWAIT

  ! for surface elements exactly on the ICB
!$OMP DO
  do ispec2D = 1, nspec_bottom ! NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

    ispec = ibelm_bottom_outer_core(ispec2D)
    ispec_selected = ibelm_top_inner_core(ispec2D)

    ! only for DOFs exactly on the ICB (bottom of these elements)
    k = 1
    ! get displacement on the solid side using pointwise matching
    k_corresp = NGLLZ

    do j = 1,NGLLY
      do i = 1,NGLLX

        ! corresponding points are located at the bottom of the mantle
        iglob_ic = ibool_inner_core(i,j,k_corresp,ispec_selected)

        displ_x = displ_inner_core(1,iglob_ic)
        displ_y = displ_inner_core(2,iglob_ic)
        displ_z = displ_inner_core(3,iglob_ic)

        ! get normal on the ICB
        nx = normal_bottom_outer_core(1,i,j,ispec2D)
        ny = normal_bottom_outer_core(2,i,j,ispec2D)
        nz = normal_bottom_outer_core(3,i,j,ispec2D)

        ! compute dot product
        displ_n = displ_x*nx + displ_y*ny + displ_z*nz

        ! formulation with generalized potential
        weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        ! get global point number
        iglob_oc = ibool_outer_core(i,j,k,ispec)

        ! update fluid acceleration/pressure
!$OMP ATOMIC
        accel_outer_core(iglob_oc) = accel_outer_core(iglob_oc) - weight*displ_n
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_fluid_CMB_ICB

!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_coupling_fluid_CMB(NGLOB_CM,displ_crust_mantle, &
                                        ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                        NGLOB_OC,accel_outer_core, &
                                        normal_top_outer_core,jacobian2D_top_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                        nspec2D_top)

  use constants_solver

  implicit none

  integer :: NGLOB_CM
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM) :: displ_crust_mantle

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM) :: ibelm_bottom_crust_mantle

  integer :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC) :: normal_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core
  integer, dimension(NSPEC2D_TOP_OC) :: ibelm_top_outer_core

  integer :: nspec2D_top

  ! local parameters
  real(kind=CUSTOM_REAL) :: displ_x,displ_y,displ_z,displ_n,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob_cm,iglob_oc,ispec_selected

  ! for surface elements exactly on the CMB

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

! openmp solver
!$OMP PARALLEL if (nspec2D_top > 500) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,ispec_selected,i,j,k,k_corresp,iglob_cm,iglob_oc, &
!$OMP displ_x,displ_y,displ_z,nx,ny,nz,displ_n,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy)
!$OMP DO
  do ispec2D = 1,nspec2D_top !NSPEC2D_TOP(IREGION_OUTER_CORE)

    ispec = ibelm_top_outer_core(ispec2D)
    ispec_selected = ibelm_bottom_crust_mantle(ispec2D)

    ! only for DOFs exactly on the CMB (top of these elements)
    k = NGLLZ

    ! get displacement on the solid side using pointwise matching
    k_corresp = 1

    do j = 1,NGLLY
      do i = 1,NGLLX
        ! corresponding points are located at the bottom of the mantle
        iglob_cm = ibool_crust_mantle(i,j,k_corresp,ispec_selected)

        displ_x = displ_crust_mantle(1,iglob_cm)
        displ_y = displ_crust_mantle(2,iglob_cm)
        displ_z = displ_crust_mantle(3,iglob_cm)

        ! get normal on the CMB
        nx = normal_top_outer_core(1,i,j,ispec2D)
        ny = normal_top_outer_core(2,i,j,ispec2D)
        nz = normal_top_outer_core(3,i,j,ispec2D)

        ! compute dot product
        displ_n = displ_x*nx + displ_y*ny + displ_z*nz

        ! formulation with generalized potential
        weight = jacobian2D_top_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        ! get global point number
        iglob_oc = ibool_outer_core(i,j,k,ispec)

        ! update fluid acceleration/pressure
!$OMP ATOMIC
        accel_outer_core(iglob_oc) = accel_outer_core(iglob_oc) + weight*displ_n
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_fluid_CMB

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_fluid_ICB(NGLOB_IC,displ_inner_core, &
                                        ibool_inner_core,ibelm_top_inner_core, &
                                        NGLOB_OC,accel_outer_core, &
                                        normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                        nspec_bottom)

  use constants_solver

  implicit none

  integer :: NGLOB_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC) :: &
    displ_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

  integer :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core

  integer :: nspec_bottom

  ! local parameters
  real(kind=CUSTOM_REAL) :: displ_x,displ_y,displ_z,displ_n,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob_oc,iglob_ic,ispec_selected

  ! for surface elements exactly on the ICB

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

! openmp solver
!$OMP PARALLEL if (nspec_bottom > 500) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,ispec_selected,i,j,k,k_corresp,iglob_ic,iglob_oc, &
!$OMP displ_x,displ_y,displ_z,nx,ny,nz,displ_n,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy)
!$OMP DO
  do ispec2D = 1, nspec_bottom ! NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

    ispec = ibelm_bottom_outer_core(ispec2D)
    ispec_selected = ibelm_top_inner_core(ispec2D)

    ! only for DOFs exactly on the ICB (bottom of these elements)
    k = 1
    ! get displacement on the solid side using pointwise matching
    k_corresp = NGLLZ

    do j = 1,NGLLY
      do i = 1,NGLLX

        ! corresponding points are located at the bottom of the mantle
        iglob_ic = ibool_inner_core(i,j,k_corresp,ispec_selected)

        displ_x = displ_inner_core(1,iglob_ic)
        displ_y = displ_inner_core(2,iglob_ic)
        displ_z = displ_inner_core(3,iglob_ic)

        ! get normal on the ICB
        nx = normal_bottom_outer_core(1,i,j,ispec2D)
        ny = normal_bottom_outer_core(2,i,j,ispec2D)
        nz = normal_bottom_outer_core(3,i,j,ispec2D)

        ! compute dot product
        displ_n = displ_x*nx + displ_y*ny + displ_z*nz

        ! formulation with generalized potential
        weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        ! get global point number
        iglob_oc = ibool_outer_core(i,j,k,ispec)

        ! update fluid acceleration/pressure
!$OMP ATOMIC
        accel_outer_core(iglob_oc) = accel_outer_core(iglob_oc) - weight*displ_n
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_fluid_ICB




!
!-------------------------------------------------------------------------------------------------
!
! coupling solid to fluid region
!
!-------------------------------------------------------------------------------------------------
!


  subroutine compute_coupling_solid()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore, only: accel_outer_core, &
                                  normal_top_outer_core,jacobian2D_top_outer_core, &
                                  normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                  ibelm_top_outer_core,ibelm_bottom_outer_core, &
                                  ibool_outer_core
  use specfem_par_full_gravity, only: pgrav_oc
  implicit none

  ! checks if anything to do
  if (.not. ACTUALLY_COUPLE_FLUID_CMB .and. .not. ACTUALLY_COUPLE_FLUID_ICB) return

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

  ! only for elements in first matching layer in the solid
  if (.not. GPU_MODE) then
    ! on CPU
    if (ACTUALLY_COUPLE_FLUID_CMB .and. ACTUALLY_COUPLE_FLUID_ICB) then
      ! couple both mantle and inner core with outer core
      call compute_coupling_CMB_ICB_fluid(NGLOB_OUTER_CORE,accel_outer_core,wgllwgll_xy,ibool_outer_core, &
                                          NGLOB_CRUST_MANTLE,displ_crust_mantle,accel_crust_mantle, &
                                          ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                          NGLOB_INNER_CORE,displ_inner_core,accel_inner_core, &
                                          ibool_inner_core,ibelm_top_inner_core, &
                                          RHO_TOP_OC,minus_g_cmb, &
                                          RHO_BOTTOM_OC,minus_g_icb, &
                                          normal_top_outer_core,jacobian2D_top_outer_core, &
                                          normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                          ibelm_top_outer_core,NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                          ibelm_bottom_outer_core,NSPEC2D_TOP(IREGION_INNER_CORE), &
                                          pgrav_oc)
    else
      ! separate only one coupling interface
      !---
      !--- couple with outer core at the bottom of the mantle
      !---
      if (ACTUALLY_COUPLE_FLUID_CMB) &
        call compute_coupling_CMB_fluid(NGLOB_CRUST_MANTLE,displ_crust_mantle,accel_crust_mantle, &
                                        ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                        NGLOB_OUTER_CORE,accel_outer_core, &
                                        normal_top_outer_core,jacobian2D_top_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                        RHO_TOP_OC,minus_g_cmb, &
                                        NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                        pgrav_oc)

      !---
      !--- couple with outer core at the top of the inner core
      !---
      if (ACTUALLY_COUPLE_FLUID_ICB) &
        call compute_coupling_ICB_fluid(NGLOB_INNER_CORE,displ_inner_core,accel_inner_core, &
                                        ibool_inner_core,ibelm_top_inner_core, &
                                        NGLOB_OUTER_CORE,accel_outer_core, &
                                        normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                        RHO_BOTTOM_OC,minus_g_icb, &
                                        NSPEC2D_TOP(IREGION_INNER_CORE), &
                                        pgrav_oc)
    endif
  else
    ! on GPU
    !---
    !--- couple with outer core at the bottom of the mantle
    !---
    if (ACTUALLY_COUPLE_FLUID_CMB ) &
         call compute_coupling_cmb_fluid_gpu(Mesh_pointer,1)
    !---
    !--- couple with outer core at the top of the inner core
    !---
    if (ACTUALLY_COUPLE_FLUID_ICB ) &
         call compute_coupling_icb_fluid_gpu(Mesh_pointer,1)

  endif

  end subroutine compute_coupling_solid

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_solid_backward()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore, only: b_accel_outer_core, &
                                  normal_top_outer_core,jacobian2D_top_outer_core, &
                                  normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                  ibelm_top_outer_core,ibelm_bottom_outer_core, &
                                  ibool_outer_core
  use specfem_par_full_gravity, only: b_pgrav_oc
  implicit none

  ! checks if anything to do
  if (.not. ACTUALLY_COUPLE_FLUID_CMB .and. .not. ACTUALLY_COUPLE_FLUID_ICB) return

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

  ! only for elements in first matching layer in the solid
  if (.not. GPU_MODE) then
    ! on CPU
    if (ACTUALLY_COUPLE_FLUID_CMB .and. ACTUALLY_COUPLE_FLUID_ICB) then
      ! couple both mantle and inner core with outer core
      call compute_coupling_CMB_ICB_fluid(NGLOB_OUTER_CORE_ADJOINT,b_accel_outer_core,wgllwgll_xy,ibool_outer_core, &
                                          NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle,b_accel_crust_mantle, &
                                          ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                          NGLOB_INNER_CORE_ADJOINT,b_displ_inner_core,b_accel_inner_core, &
                                          ibool_inner_core,ibelm_top_inner_core, &
                                          RHO_TOP_OC,minus_g_cmb, &
                                          RHO_BOTTOM_OC,minus_g_icb, &
                                          normal_top_outer_core,jacobian2D_top_outer_core, &
                                          normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                          ibelm_top_outer_core,NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                          ibelm_bottom_outer_core,NSPEC2D_TOP(IREGION_INNER_CORE), &
                                          b_pgrav_oc)
    else
      ! separate only one coupling interface
      !---
      !--- couple with outer core at the bottom of the mantle
      !---
      if (ACTUALLY_COUPLE_FLUID_CMB) &
        call compute_coupling_CMB_fluid(NGLOB_CRUST_MANTLE_ADJOINT,b_displ_crust_mantle,b_accel_crust_mantle, &
                                        ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                        NGLOB_OUTER_CORE_ADJOINT,b_accel_outer_core, &
                                        normal_top_outer_core,jacobian2D_top_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                        RHO_TOP_OC,minus_g_cmb, &
                                        NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE), &
                                        b_pgrav_oc)

      !---
      !--- couple with outer core at the top of the inner core
      !---
      if (ACTUALLY_COUPLE_FLUID_ICB) &
        call compute_coupling_ICB_fluid(NGLOB_INNER_CORE_ADJOINT,b_displ_inner_core,b_accel_inner_core, &
                                        ibool_inner_core,ibelm_top_inner_core, &
                                        NGLOB_OUTER_CORE_ADJOINT,b_accel_outer_core, &
                                        normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                        RHO_BOTTOM_OC,minus_g_icb, &
                                        NSPEC2D_TOP(IREGION_INNER_CORE), &
                                        b_pgrav_oc)

    endif
  else
    ! on GPU
    !---
    !--- couple with outer core at the bottom of the mantle
    !---
    if (ACTUALLY_COUPLE_FLUID_CMB ) &
         call compute_coupling_cmb_fluid_gpu(Mesh_pointer,3)
    !---
    !--- couple with outer core at the top of the inner core
    !---
    if (ACTUALLY_COUPLE_FLUID_ICB ) &
         call compute_coupling_icb_fluid_gpu(Mesh_pointer,3)

  endif

  end subroutine compute_coupling_solid_backward

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_CMB_ICB_fluid(NGLOB_OC,accel_outer_core,wgllwgll_xy,ibool_outer_core, &
                                            NGLOB_CM,displ_crust_mantle,accel_crust_mantle, &
                                            ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                            NGLOB_IC,displ_inner_core,accel_inner_core, &
                                            ibool_inner_core,ibelm_top_inner_core, &
                                            RHO_TOP_OC,minus_g_cmb, &
                                            RHO_BOTTOM_OC,minus_g_icb, &
                                            normal_top_outer_core,jacobian2D_top_outer_core, &
                                            normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                            ibelm_top_outer_core,nspec_bottom, &
                                            ibelm_bottom_outer_core,nspec2D_top, &
                                            pgrav_oc)

  use constants_solver

  implicit none

  integer, intent(in) :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC), intent(in) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY), intent(in) :: wgllwgll_xy
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), intent(in) :: ibool_outer_core

  ! crust/mantel
  integer, intent(in) :: NGLOB_CM
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM), intent(in) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM), intent(inout) :: accel_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), intent(in) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM), intent(in) :: ibelm_bottom_crust_mantle

  ! inner core
  integer, intent(in) :: NGLOB_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC), intent(in) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC), intent(inout) :: accel_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), intent(in) :: ibool_inner_core
  integer, dimension(NSPEC2D_TOP_IC), intent(in) :: ibelm_top_inner_core

  ! outer core
  double precision, intent(in) :: RHO_TOP_OC
  real(kind=CUSTOM_REAL), intent(in) :: minus_g_cmb
  double precision, intent(in) :: RHO_BOTTOM_OC
  real(kind=CUSTOM_REAL), intent(in) :: minus_g_icb

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC), intent(in) :: normal_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC), intent(in) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC), intent(in) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC), intent(in) :: jacobian2D_bottom_outer_core

  integer, dimension(NSPEC2D_TOP_OC), intent(in) :: ibelm_top_outer_core
  integer, intent(in) :: nspec_bottom
  integer, dimension(NSPEC2D_BOTTOM_OC), intent(in) :: ibelm_bottom_outer_core
  integer, intent(in) :: nspec2D_top

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC), intent(in) :: pgrav_oc(NGLOB_OC)

  ! local parameters
  real(kind=CUSTOM_REAL) :: pressure,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob,iglob_mantle,iglob_inner_core,ispec_selected

  ! for surface elements exactly on the CMB

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

! openmp solver
!$OMP PARALLEL if (nspec_bottom > 500) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,ispec_selected,i,j,k,k_corresp,iglob,iglob_mantle,iglob_inner_core, &
!$OMP nx,ny,nz,pressure,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy)
!$OMP DO
  do ispec2D = 1,nspec_bottom ! NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)

    ispec = ibelm_bottom_crust_mantle(ispec2D)
    ispec_selected = ibelm_top_outer_core(ispec2D)

    ! only for DOFs exactly on the CMB (bottom of these elements)
    k = 1
    ! get potential on the fluid side using pointwise matching
    k_corresp = NGLLZ

    do j = 1,NGLLY
      do i = 1,NGLLX
        ! get normal at the CMB
        nx = normal_top_outer_core(1,i,j,ispec2D)
        ny = normal_top_outer_core(2,i,j,ispec2D)
        nz = normal_top_outer_core(3,i,j,ispec2D)

        ! get global point number
        ! corresponding points are located at the top of the outer core
        iglob = ibool_outer_core(i,j,k_corresp,ispec_selected)
        iglob_mantle = ibool_crust_mantle(i,j,k,ispec)

        ! compute pressure, taking gravity into account
        if (GRAVITY_VAL) then
          if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
            pressure = real(RHO_TOP_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_cmb *(displ_crust_mantle(1,iglob_mantle)*nx &
                             + displ_crust_mantle(2,iglob_mantle)*ny &
                             + displ_crust_mantle(3,iglob_mantle)*nz) - pgrav_oc(iglob))
          else
            pressure = real(RHO_TOP_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_cmb *(displ_crust_mantle(1,iglob_mantle)*nx &
                             + displ_crust_mantle(2,iglob_mantle)*ny &
                             + displ_crust_mantle(3,iglob_mantle)*nz))
          endif
        else
          pressure = - real(RHO_TOP_OC,kind=CUSTOM_REAL) * accel_outer_core(iglob)
        endif

        ! formulation with generalized potential
        weight = jacobian2D_top_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

!$OMP ATOMIC
        accel_crust_mantle(1,iglob_mantle) = accel_crust_mantle(1,iglob_mantle) + weight*nx*pressure
!$OMP ATOMIC
        accel_crust_mantle(2,iglob_mantle) = accel_crust_mantle(2,iglob_mantle) + weight*ny*pressure
!$OMP ATOMIC
        accel_crust_mantle(3,iglob_mantle) = accel_crust_mantle(3,iglob_mantle) + weight*nz*pressure
      enddo
    enddo
  enddo
!$OMP ENDDO NOWAIT

!$OMP DO
  do ispec2D = 1,nspec2D_top ! NSPEC2D_TOP(IREGION_INNER_CORE)

    ispec = ibelm_top_inner_core(ispec2D)
    ispec_selected = ibelm_bottom_outer_core(ispec2D)

    ! only for DOFs exactly on the ICB (top of these elements)
    k = NGLLZ
    ! get potential on the fluid side using pointwise matching
    k_corresp = 1

    do j = 1,NGLLY
      do i = 1,NGLLX
        ! get normal at the ICB
        nx = normal_bottom_outer_core(1,i,j,ispec2D)
        ny = normal_bottom_outer_core(2,i,j,ispec2D)
        nz = normal_bottom_outer_core(3,i,j,ispec2D)

        ! get global point number
        ! corresponding points are located at the bottom of the outer core
        iglob = ibool_outer_core(i,j,k_corresp,ispec_selected)
        iglob_inner_core = ibool_inner_core(i,j,k,ispec)

        ! compute pressure, taking gravity into account
        if (GRAVITY_VAL) then
          if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
            pressure = real(RHO_BOTTOM_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_icb *(displ_inner_core(1,iglob_inner_core)*nx &
                             + displ_inner_core(2,iglob_inner_core)*ny &
                             + displ_inner_core(3,iglob_inner_core)*nz) - pgrav_oc(iglob))
          else
            pressure = real(RHO_BOTTOM_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_icb *(displ_inner_core(1,iglob_inner_core)*nx &
                             + displ_inner_core(2,iglob_inner_core)*ny &
                             + displ_inner_core(3,iglob_inner_core)*nz))
          endif
        else
          pressure = - real(RHO_BOTTOM_OC,kind=CUSTOM_REAL) * accel_outer_core(iglob)
        endif

        ! formulation with generalized potential
        weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

!$OMP ATOMIC
        accel_inner_core(1,iglob_inner_core) = accel_inner_core(1,iglob_inner_core) - weight*nx*pressure
!$OMP ATOMIC
        accel_inner_core(2,iglob_inner_core) = accel_inner_core(2,iglob_inner_core) - weight*ny*pressure
!$OMP ATOMIC
        accel_inner_core(3,iglob_inner_core) = accel_inner_core(3,iglob_inner_core) - weight*nz*pressure

      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_CMB_ICB_fluid


!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_CMB_fluid(NGLOB_CM,displ_crust_mantle,accel_crust_mantle, &
                                        ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                        NGLOB_OC,accel_outer_core, &
                                        normal_top_outer_core,jacobian2D_top_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                        RHO_TOP_OC,minus_g_cmb, &
                                        nspec_bottom, &
                                        pgrav_oc)

  use constants_solver

  implicit none

  integer, intent(in) :: NGLOB_CM
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM), intent(in) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM), intent(inout) :: accel_crust_mantle

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), intent(in) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM), intent(in) :: ibelm_bottom_crust_mantle

  integer, intent(in) :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC), intent(in) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC), intent(in) :: normal_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC), intent(in) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY), intent(in) :: wgllwgll_xy

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), intent(in) :: ibool_outer_core
  integer, dimension(NSPEC2D_TOP_OC), intent(in) :: ibelm_top_outer_core

  double precision, intent(in) :: RHO_TOP_OC
  real(kind=CUSTOM_REAL), intent(in) :: minus_g_cmb

  integer, intent(in) :: nspec_bottom

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC), intent(in) :: pgrav_oc(NGLOB_OC)

  ! local parameters
  real(kind=CUSTOM_REAL) :: pressure,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob,iglob_mantle,ispec_selected

  ! for surface elements exactly on the CMB

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

! openmp solver
!$OMP PARALLEL if (nspec_bottom > 500) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,ispec_selected,i,j,k,k_corresp,iglob,iglob_mantle, &
!$OMP nx,ny,nz,pressure,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy)
!$OMP DO
  do ispec2D = 1,nspec_bottom ! NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)

    ispec = ibelm_bottom_crust_mantle(ispec2D)
    ispec_selected = ibelm_top_outer_core(ispec2D)

    ! only for DOFs exactly on the CMB (bottom of these elements)
    k = 1
    ! get potential on the fluid side using pointwise matching
    k_corresp = NGLLZ

    do j = 1,NGLLY
      do i = 1,NGLLX
        ! get normal at the CMB
        nx = normal_top_outer_core(1,i,j,ispec2D)
        ny = normal_top_outer_core(2,i,j,ispec2D)
        nz = normal_top_outer_core(3,i,j,ispec2D)

        ! get global point number
        ! corresponding points are located at the top of the outer core
        iglob = ibool_outer_core(i,j,k_corresp,ispec_selected)
        iglob_mantle = ibool_crust_mantle(i,j,k,ispec)

        ! compute pressure, taking gravity into account
        if (GRAVITY_VAL) then
          if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
            pressure = real(RHO_TOP_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_cmb *(displ_crust_mantle(1,iglob_mantle)*nx &
                             + displ_crust_mantle(2,iglob_mantle)*ny &
                             + displ_crust_mantle(3,iglob_mantle)*nz) - pgrav_oc(iglob))
          else
            pressure = real(RHO_TOP_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_cmb *(displ_crust_mantle(1,iglob_mantle)*nx &
                             + displ_crust_mantle(2,iglob_mantle)*ny &
                             + displ_crust_mantle(3,iglob_mantle)*nz))
          endif
        else
          pressure = - real(RHO_TOP_OC,kind=CUSTOM_REAL) * accel_outer_core(iglob)
        endif

        ! formulation with generalized potential
        weight = jacobian2D_top_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

!$OMP ATOMIC
        accel_crust_mantle(1,iglob_mantle) = accel_crust_mantle(1,iglob_mantle) + weight*nx*pressure
!$OMP ATOMIC
        accel_crust_mantle(2,iglob_mantle) = accel_crust_mantle(2,iglob_mantle) + weight*ny*pressure
!$OMP ATOMIC
        accel_crust_mantle(3,iglob_mantle) = accel_crust_mantle(3,iglob_mantle) + weight*nz*pressure
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_CMB_fluid


!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_ICB_fluid(NGLOB_IC,displ_inner_core,accel_inner_core, &
                                        ibool_inner_core,ibelm_top_inner_core, &
                                        NGLOB_OC,accel_outer_core, &
                                        normal_bottom_outer_core,jacobian2D_bottom_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_bottom_outer_core, &
                                        RHO_BOTTOM_OC,minus_g_icb, &
                                        nspec2D_top, &
                                        pgrav_oc)

  use constants_solver

  implicit none

  integer, intent(in) :: NGLOB_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC), intent(in) :: displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC), intent(inout) :: accel_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), intent(in) :: ibool_inner_core
  integer, dimension(NSPEC2D_TOP_IC), intent(in) :: ibelm_top_inner_core

  integer :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC), intent(in) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC), intent(in) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC), intent(in) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY), intent(in) :: wgllwgll_xy

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), intent(in) :: ibool_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC), intent(in) :: ibelm_bottom_outer_core

  double precision, intent(in) :: RHO_BOTTOM_OC
  real(kind=CUSTOM_REAL), intent(in) :: minus_g_icb

  integer, intent(in) :: nspec2D_top

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC), intent(in) :: pgrav_oc(NGLOB_OC)

  ! local parameters
  real(kind=CUSTOM_REAL) :: pressure,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob,iglob_inner_core,ispec_selected

  ! for surface elements exactly on the ICB

  ! checks if anything to do
  if (NSPEC_OUTER_CORE == 0) return

! openmp solver
!$OMP PARALLEL if (nspec2D_top > 500) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ispec2D,ispec,ispec_selected,i,j,k,k_corresp,iglob,iglob_inner_core, &
!$OMP nx,ny,nz,pressure,weight) &
!$OMP FIRSTPRIVATE(wgllwgll_xy)
!$OMP DO
  do ispec2D = 1,nspec2D_top ! NSPEC2D_TOP(IREGION_INNER_CORE)

    ispec = ibelm_top_inner_core(ispec2D)
    ispec_selected = ibelm_bottom_outer_core(ispec2D)

    ! only for DOFs exactly on the ICB (top of these elements)
    k = NGLLZ
    ! get potential on the fluid side using pointwise matching
    k_corresp = 1

    do j = 1,NGLLY
      do i = 1,NGLLX
        ! get normal at the ICB
        nx = normal_bottom_outer_core(1,i,j,ispec2D)
        ny = normal_bottom_outer_core(2,i,j,ispec2D)
        nz = normal_bottom_outer_core(3,i,j,ispec2D)

        ! get global point number
        ! corresponding points are located at the bottom of the outer core
        iglob = ibool_outer_core(i,j,k_corresp,ispec_selected)
        iglob_inner_core = ibool_inner_core(i,j,k,ispec)

        ! compute pressure, taking gravity into account
        if (GRAVITY_VAL) then
          if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
            pressure = real(RHO_BOTTOM_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_icb *(displ_inner_core(1,iglob_inner_core)*nx &
                             + displ_inner_core(2,iglob_inner_core)*ny &
                             + displ_inner_core(3,iglob_inner_core)*nz) - pgrav_oc(iglob))
          else
            pressure = real(RHO_BOTTOM_OC,kind=CUSTOM_REAL) * (- accel_outer_core(iglob) &
               + minus_g_icb *(displ_inner_core(1,iglob_inner_core)*nx &
                             + displ_inner_core(2,iglob_inner_core)*ny &
                             + displ_inner_core(3,iglob_inner_core)*nz))
          endif
        else
          pressure = - real(RHO_BOTTOM_OC,kind=CUSTOM_REAL) * accel_outer_core(iglob)
        endif

        ! formulation with generalized potential
        weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

!$OMP ATOMIC
        accel_inner_core(1,iglob_inner_core) = accel_inner_core(1,iglob_inner_core) - weight*nx*pressure
!$OMP ATOMIC
        accel_inner_core(2,iglob_inner_core) = accel_inner_core(2,iglob_inner_core) - weight*ny*pressure
!$OMP ATOMIC
        accel_inner_core(3,iglob_inner_core) = accel_inner_core(3,iglob_inner_core) - weight*nz*pressure
      enddo
    enddo
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_ICB_fluid

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_ocean(NGLOB,accel_crust_mantle, &
                                    rmassx_crust_mantle, rmassy_crust_mantle, rmassz_crust_mantle, &
                                    npoin_oceans,ibool_ocean_load,rmass_ocean_load_selected,normal_ocean_load)

  use constants_solver

  implicit none

  integer, intent(in) :: NGLOB
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB), intent(inout) :: accel_crust_mantle

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be pointers to it
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE), intent(in) :: rmassx_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE), intent(in) :: rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE), intent(in) :: rmassz_crust_mantle

  ! oceans arrays
  integer, intent(in) :: npoin_oceans
  real(kind=CUSTOM_REAL), dimension(npoin_oceans), intent(in) :: rmass_ocean_load_selected
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin_oceans), intent(in) :: normal_ocean_load
  integer, dimension(npoin_oceans), intent(in) :: ibool_ocean_load

  ! local parameters
  real(kind=CUSTOM_REAL) :: force_normal_comp,rmass
  real(kind=CUSTOM_REAL) :: additional_term_x,additional_term_y,additional_term_z
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  integer :: ipoin,iglob

  ! for surface elements exactly at the top of the crust (ocean bottom)

  ! checks if anything to do
  if (npoin_oceans == 0) return

! openmp solver
!$OMP PARALLEL if (npoin_oceans > 1000) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(ipoin,iglob, &
!$OMP nx,ny,nz,force_normal_comp,rmass,additional_term_x,additional_term_y,additional_term_z)
!$OMP DO
  do ipoin = 1,npoin_oceans !NSPEC2D_TOP(IREGION_CRUST_MANTLE)

    ! get global point number
    iglob = ibool_ocean_load(ipoin)

    ! get normal
    nx = normal_ocean_load(1,ipoin)
    ny = normal_ocean_load(2,ipoin)
    nz = normal_ocean_load(3,ipoin)

    ! make updated component of right-hand side
    ! we divide by rmass_crust_mantle() which is 1 / M
    ! we use the total force which includes the Coriolis term above
    force_normal_comp = accel_crust_mantle(1,iglob)*nx / rmassx_crust_mantle(iglob) &
                      + accel_crust_mantle(2,iglob)*ny / rmassy_crust_mantle(iglob) &
                      + accel_crust_mantle(3,iglob)*nz / rmassz_crust_mantle(iglob)

    rmass = rmass_ocean_load_selected(ipoin)

    additional_term_x = (rmass - rmassx_crust_mantle(iglob)) * force_normal_comp
    additional_term_y = (rmass - rmassy_crust_mantle(iglob)) * force_normal_comp
    additional_term_z = (rmass - rmassz_crust_mantle(iglob)) * force_normal_comp

!$OMP ATOMIC
    accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + additional_term_x * nx
!$OMP ATOMIC
    accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + additional_term_y * ny
!$OMP ATOMIC
    accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + additional_term_z * nz
  enddo
!$OMP ENDDO
!$OMP END PARALLEL

  end subroutine compute_coupling_ocean

