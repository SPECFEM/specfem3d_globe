!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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
        accel_outer_core(iglob_oc) = accel_outer_core(iglob_oc) + weight*displ_n
      enddo
    enddo
  enddo

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
        accel_outer_core(iglob_oc) = accel_outer_core(iglob_oc) - weight*displ_n
      enddo
    enddo
  enddo

  end subroutine compute_coupling_fluid_ICB

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_CMB_fluid(NGLOB_CM,displ_crust_mantle,accel_crust_mantle, &
                                        ibool_crust_mantle,ibelm_bottom_crust_mantle, &
                                        NGLOB_OC,accel_outer_core, &
                                        normal_top_outer_core,jacobian2D_top_outer_core, &
                                        wgllwgll_xy,ibool_outer_core,ibelm_top_outer_core, &
                                        RHO_TOP_OC,minus_g_cmb, &
                                        nspec_bottom)

  use constants_solver

  implicit none

  integer :: NGLOB_CM
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CM) :: displ_crust_mantle,accel_crust_mantle

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM) :: ibelm_bottom_crust_mantle

  integer :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC) :: normal_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core
  integer, dimension(NSPEC2D_TOP_OC) :: ibelm_top_outer_core

  double precision :: RHO_TOP_OC
  real(kind=CUSTOM_REAL) :: minus_g_cmb

  integer :: nspec_bottom

  ! local parameters
  real(kind=CUSTOM_REAL) :: pressure,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob_oc,iglob_mantle,ispec_selected

  ! for surface elements exactly on the CMB
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
        iglob_oc = ibool_outer_core(i,j,k_corresp,ispec_selected)
        iglob_mantle = ibool_crust_mantle(i,j,k,ispec)

        ! compute pressure, taking gravity into account
        if (GRAVITY_VAL) then
          pressure = RHO_TOP_OC * (- accel_outer_core(iglob_oc) &
             + minus_g_cmb *(displ_crust_mantle(1,iglob_mantle)*nx &
                            + displ_crust_mantle(2,iglob_mantle)*ny &
                            + displ_crust_mantle(3,iglob_mantle)*nz))
        else
          pressure = - RHO_TOP_OC * accel_outer_core(iglob_oc)
        endif

        ! formulation with generalized potential
        weight = jacobian2D_top_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        accel_crust_mantle(1,iglob_mantle) = accel_crust_mantle(1,iglob_mantle) + weight*nx*pressure
        accel_crust_mantle(2,iglob_mantle) = accel_crust_mantle(2,iglob_mantle) + weight*ny*pressure
        accel_crust_mantle(3,iglob_mantle) = accel_crust_mantle(3,iglob_mantle) + weight*nz*pressure
      enddo
    enddo
  enddo

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
                                        nspec2D_top)

  use constants_solver

  implicit none

  integer :: NGLOB_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_IC) :: displ_inner_core,accel_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

  integer :: NGLOB_OC
  real(kind=CUSTOM_REAL), dimension(NGLOB_OC) :: accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core

  double precision :: RHO_BOTTOM_OC
  real(kind=CUSTOM_REAL) :: minus_g_icb

  integer :: nspec2D_top

  ! local parameters
  real(kind=CUSTOM_REAL) :: pressure,nx,ny,nz,weight
  integer :: i,j,k,k_corresp,ispec,ispec2D,iglob,iglob_inner_core,ispec_selected

  ! for surface elements exactly on the ICB
  do ispec2D = 1,nspec2D_top ! NSPEC2D_TOP(IREGION_INNER_CORE)

    ispec = ibelm_top_inner_core(ispec2D)
    ispec_selected = ibelm_bottom_outer_core(ispec2D)

    ! only for DOFs exactly on the ICB (top of these elements)
    k = NGLLZ
    ! get velocity potential on the fluid side using pointwise matching
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
          pressure = RHO_BOTTOM_OC * (- accel_outer_core(iglob) &
             + minus_g_icb *(displ_inner_core(1,iglob_inner_core)*nx &
                             + displ_inner_core(2,iglob_inner_core)*ny &
                             + displ_inner_core(3,iglob_inner_core)*nz))
        else
          pressure = - RHO_BOTTOM_OC * accel_outer_core(iglob)
        endif

        ! formulation with generalized potential
        weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        accel_inner_core(1,iglob_inner_core) = accel_inner_core(1,iglob_inner_core) - weight*nx*pressure
        accel_inner_core(2,iglob_inner_core) = accel_inner_core(2,iglob_inner_core) - weight*ny*pressure
        accel_inner_core(3,iglob_inner_core) = accel_inner_core(3,iglob_inner_core) - weight*nz*pressure

      enddo
    enddo
  enddo

  end subroutine compute_coupling_ICB_fluid

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_coupling_ocean(NGLOB,accel_crust_mantle, &
                                    rmassx_crust_mantle, rmassy_crust_mantle, rmassz_crust_mantle, &
                                    rmass_ocean_load,normal_top_crust_mantle, &
                                    ibool_crust_mantle,ibelm_top_crust_mantle, &
                                    updated_dof_ocean_load, &
                                    nspec_top)

  use constants_solver

  implicit none

  integer :: NGLOB
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_crust_mantle

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be pointers to it
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmassx_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmassz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_CM) :: normal_top_crust_mantle

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle

  logical, dimension(NGLOB_CRUST_MANTLE_OCEANS) :: updated_dof_ocean_load

  integer :: nspec_top

  ! local parameters
  real(kind=CUSTOM_REAL) :: force_normal_comp
  real(kind=CUSTOM_REAL) :: additional_term_x,additional_term_y,additional_term_z
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  integer :: i,j,k,ispec,ispec2D,iglob

  !   initialize the updates
  updated_dof_ocean_load(:) = .false.

  ! for surface elements exactly at the top of the crust (ocean bottom)
  do ispec2D = 1,nspec_top !NSPEC2D_TOP(IREGION_CRUST_MANTLE)

    ispec = ibelm_top_crust_mantle(ispec2D)

    ! only for DOFs exactly at the top of the crust (ocean bottom)
    k = NGLLZ

    do j = 1,NGLLY
      do i = 1,NGLLX

        ! get global point number
        iglob = ibool_crust_mantle(i,j,k,ispec)

        ! only update once
        if (.not. updated_dof_ocean_load(iglob)) then

          ! get normal
          nx = normal_top_crust_mantle(1,i,j,ispec2D)
          ny = normal_top_crust_mantle(2,i,j,ispec2D)
          nz = normal_top_crust_mantle(3,i,j,ispec2D)

          ! make updated component of right-hand side
          ! we divide by rmass_crust_mantle() which is 1 / M
          ! we use the total force which includes the Coriolis term above
          force_normal_comp = accel_crust_mantle(1,iglob)*nx / rmassx_crust_mantle(iglob) &
                            + accel_crust_mantle(2,iglob)*ny / rmassy_crust_mantle(iglob) &
                            + accel_crust_mantle(3,iglob)*nz / rmassz_crust_mantle(iglob)

          additional_term_x = (rmass_ocean_load(iglob) - rmassx_crust_mantle(iglob)) * force_normal_comp
          additional_term_y = (rmass_ocean_load(iglob) - rmassy_crust_mantle(iglob)) * force_normal_comp
          additional_term_z = (rmass_ocean_load(iglob) - rmassz_crust_mantle(iglob)) * force_normal_comp

          accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + additional_term_x * nx
          accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + additional_term_y * ny
          accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + additional_term_z * nz

          ! done with this point
          updated_dof_ocean_load(iglob) = .true.

        endif

      enddo
    enddo
  enddo

  end subroutine compute_coupling_ocean

