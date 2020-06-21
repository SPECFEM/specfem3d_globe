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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



!--------------------------------------------------------------------------------------------
!
! isotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_iso(ispec, &
                                 gravity_pre_store,gravity_H, &
                                 deriv, &
                                 wgll_cube, &
                                 kappavstore,muvstore, &
                                 ibool, &
                                 R_xx,R_yy,R_xy,R_xz,R_yz, &
                                 epsilon_trace_over_3, &
                                 tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 epsilondev_loc,rho_s_H)

! isotropic element in crust/mantle region

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,FOUR_THIRDS

  use constants_solver, only: &
    NSPEC => NSPEC_CRUST_MANTLE, &
    NGLOB => NGLOB_CRUST_MANTLE, &
    NSPEC_ATTENUATION => NSPEC_CRUST_MANTLE_ATTENUATION, &
    NSPEC_STRAIN_ONLY => NSPEC_CRUST_MANTLE_STRAIN_ONLY, &
    ATTENUATION_VAL, &
    PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL

  use specfem_par, only: COMPUTE_AND_STORE_STRAIN

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: deriv

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: kappavstore,muvstore

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(in) :: gravity_pre_store
  real(kind=CUSTOM_REAL),dimension(6,NGLOB),intent(in) :: gravity_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(inout) :: epsilondev_loc

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl, duydyl, duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif
! note: profiling shows that this routine takes about 60% of the total time, another 30% is spend in the tiso routine below..


  ! isotropic element

  ! precomputes factors
  call compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                          deriv(:,:,:,:,ispec),jacobianl, &
                                          duxdxl,duydyl,duzdzl, &
                                          duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                                          duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)

  ! compute deviatoric strain
  if (COMPUTE_AND_STORE_STRAIN) then
    call compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
                                           duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
                                           ispec,NSPEC_STRAIN_ONLY, &
                                           epsilon_trace_over_3,epsilondev_loc)
  endif

  !
  ! compute  isotropic  elements
  !
  DO_LOOP_IJK
    ! layer with no transverse isotropy, use kappav and muv
    kappal = kappavstore(INDEX_IJK,ispec)
    mul = muvstore(INDEX_IJK,ispec)

    lambdalplus2mul = kappal + FOUR_THIRDS * mul
    lambdal = lambdalplus2mul - 2.0_CUSTOM_REAL*mul

    ! compute stress sigma
    sigma_xx(INDEX_IJK) = lambdalplus2mul*duxdxl(INDEX_IJK) + lambdal*duydyl_plus_duzdzl(INDEX_IJK)
    sigma_yy(INDEX_IJK) = lambdalplus2mul*duydyl(INDEX_IJK) + lambdal*duxdxl_plus_duzdzl(INDEX_IJK)
    sigma_zz(INDEX_IJK) = lambdalplus2mul*duzdzl(INDEX_IJK) + lambdal*duxdxl_plus_duydyl(INDEX_IJK)

    sigma_xy(INDEX_IJK) = mul*duxdyl_plus_duydxl(INDEX_IJK)
    sigma_xz(INDEX_IJK) = mul*duzdxl_plus_duxdzl(INDEX_IJK)
    sigma_yz(INDEX_IJK) = mul*duzdyl_plus_duydzl(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! attenuation contribution to stress
  ! subtract memory variables if attenuation
  if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
    call compute_element_stress_attenuation_contrib(R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                                    R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
  endif

  ! define symmetric components of sigma (to be general in case of gravity)
  DO_LOOP_IJK
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! compute non-symmetric terms for gravity
  if (GRAVITY_VAL) then
    call compute_element_gravity(ispec,NSPEC,NGLOB,ibool,jacobianl,wgll_cube, &
                                 gravity_pre_store,gravity_H, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product of stress tensor with test vector, non-symmetric form
  call compute_element_dot_product_stress(deriv(:,:,:,:,ispec),jacobianl, &
                                                sigma_xx,sigma_yy,sigma_zz, &
                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)

  end subroutine compute_element_iso


!--------------------------------------------------------------------------------------------

  subroutine compute_element_iso_ic(ispec, &
                                    gravity_pre_store,gravity_H, &
                                    deriv, &
                                    wgll_cube, &
                                    kappavstore,muvstore, &
                                    ibool, &
                                    R_xx,R_yy,R_xy,R_xz,R_yz, &
                                    epsilon_trace_over_3, &
                                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                    dummyx_loc,dummyy_loc,dummyz_loc, &
                                    epsilondev_loc,rho_s_H)

! isotropic element in inner core

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,FOUR_THIRDS

  use constants_solver, only: &
    NSPEC => NSPEC_INNER_CORE, &
    NGLOB => NGLOB_INNER_CORE, &
    NSPEC_ATTENUATION => NSPEC_INNER_CORE_ATTENUATION, &
    NSPEC_STRAIN_ONLY => NSPEC_INNER_CORE_STRAIN_ONLY, &
    ATTENUATION_VAL, &
    PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL

  use specfem_par, only: COMPUTE_AND_STORE_STRAIN

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: deriv

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: kappavstore,muvstore

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(in) :: gravity_pre_store
  real(kind=CUSTOM_REAL),dimension(6,NGLOB),intent(in) :: gravity_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(inout) :: epsilondev_loc

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl, duydyl, duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL) :: lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) :: kappal

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif
! note: profiling shows that this routine takes about 60% of the total time, another 30% is spend in the tiso routine below..


  ! isotropic element

  ! precomputes factors
  call compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                          deriv(:,:,:,:,ispec),jacobianl, &
                                          duxdxl,duydyl,duzdzl, &
                                          duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                                          duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)

  ! compute deviatoric strain
  if (COMPUTE_AND_STORE_STRAIN) then
    call compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
                                           duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
                                           ispec,NSPEC_STRAIN_ONLY, &
                                           epsilon_trace_over_3,epsilondev_loc)
  endif

  !
  ! compute  isotropic  elements
  !
  DO_LOOP_IJK
    ! inner core with no anisotropy, use kappav and muv for instance
    ! layer with no anisotropy, use kappav and muv for instance
    kappal = kappavstore(INDEX_IJK,ispec)
    mul = muvstore(INDEX_IJK,ispec)

    lambdalplus2mul = kappal + FOUR_THIRDS * mul
    lambdal = lambdalplus2mul - 2.0_CUSTOM_REAL*mul

    ! compute stress sigma
    sigma_xx(INDEX_IJK) = lambdalplus2mul*duxdxl(INDEX_IJK) + lambdal*duydyl_plus_duzdzl(INDEX_IJK)
    sigma_yy(INDEX_IJK) = lambdalplus2mul*duydyl(INDEX_IJK) + lambdal*duxdxl_plus_duzdzl(INDEX_IJK)
    sigma_zz(INDEX_IJK) = lambdalplus2mul*duzdzl(INDEX_IJK) + lambdal*duxdxl_plus_duydyl(INDEX_IJK)

    sigma_xy(INDEX_IJK) = mul*duxdyl_plus_duydxl(INDEX_IJK)
    sigma_xz(INDEX_IJK) = mul*duzdxl_plus_duxdzl(INDEX_IJK)
    sigma_yz(INDEX_IJK) = mul*duzdyl_plus_duydzl(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! attenuation contribution to stress
  ! subtract memory variables if attenuation
  if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
    call compute_element_stress_attenuation_contrib(R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                                    R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
  endif

  ! define symmetric components of sigma (to be general in case of gravity)
  DO_LOOP_IJK
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! compute non-symmetric terms for gravity
  if (GRAVITY_VAL) then
    call compute_element_gravity(ispec,NSPEC,NGLOB,ibool,jacobianl,wgll_cube, &
                                 gravity_pre_store,gravity_H, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product of stress tensor with test vector, non-symmetric form
  call compute_element_dot_product_stress(deriv(:,:,:,:,ispec),jacobianl, &
                                                sigma_xx,sigma_yy,sigma_zz, &
                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)

  end subroutine compute_element_iso_ic


!--------------------------------------------------------------------------------------------
!
! transversely isotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_tiso(ispec, &
                                  gravity_pre_store,gravity_H, &
                                  deriv, &
                                  wgll_cube, &
                                  c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                  c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                  c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                  ibool, &
                                  R_xx,R_yy,R_xy,R_xz,R_yz, &
                                  epsilon_trace_over_3, &
                                  tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                  dummyx_loc,dummyy_loc,dummyz_loc, &
                                  epsilondev_loc,rho_s_H)

! tiso element in crust/mantle

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS

  use constants_solver, only: &
    NSPEC => NSPEC_CRUST_MANTLE, &
    NGLOB => NGLOB_CRUST_MANTLE, &
    NSPECMAX_TISO => NSPECMAX_TISO_MANTLE, &
    NSPEC_ATTENUATION => NSPEC_CRUST_MANTLE_ATTENUATION, &
    NSPEC_STRAIN_ONLY => NSPEC_CRUST_MANTLE_STRAIN_ONLY, &
    ATTENUATION_VAL, &
    PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL

  use specfem_par, only: COMPUTE_AND_STORE_STRAIN

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: deriv

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO),intent(in) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(in) :: gravity_pre_store
  real(kind=CUSTOM_REAL),dimension(6,NGLOB),intent(in) :: gravity_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(inout) :: epsilondev_loc

  ! local parameters
  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) :: c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  ! local element arrays
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl, duydyl, duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif
! note: profiling shows that this routine takes about 30% of the total time, another 60% is spend in the iso routine above..

  ! transverse isotropic element

  ! precomputes factors
  call compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                          deriv(:,:,:,:,ispec),jacobianl, &
                                          duxdxl,duydyl,duzdzl, &
                                          duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                                          duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)

  ! compute deviatoric strain
  if (COMPUTE_AND_STORE_STRAIN) then
    call compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
                                           duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
                                           ispec,NSPEC_STRAIN_ONLY, &
                                           epsilon_trace_over_3,epsilondev_loc)
  endif

  !
  ! compute either transversely isotropic elements
  !
  DO_LOOP_IJK
! note: the mesh is built such that anisotropic elements are created first in anisotropic layers,
!           thus they are listed first ( see in create_regions_mesh.f90: perm_layer() ordering )
!           this is therefore still in bounds of 1:NSPECMAX_TISO_MANTLE even if NSPECMAX_TISO is less than NSPEC
    c11 = c11store(INDEX_IJK,ispec)
    c12 = c12store(INDEX_IJK,ispec)
    c13 = c13store(INDEX_IJK,ispec)
    c14 = c14store(INDEX_IJK,ispec)
    c15 = c15store(INDEX_IJK,ispec)
    c16 = c16store(INDEX_IJK,ispec)
    c22 = c22store(INDEX_IJK,ispec)
    c23 = c23store(INDEX_IJK,ispec)
    c24 = c24store(INDEX_IJK,ispec)
    c25 = c25store(INDEX_IJK,ispec)
    c26 = c26store(INDEX_IJK,ispec)
    c33 = c33store(INDEX_IJK,ispec)
    c34 = c34store(INDEX_IJK,ispec)
    c35 = c35store(INDEX_IJK,ispec)
    c36 = c36store(INDEX_IJK,ispec)
    c44 = c44store(INDEX_IJK,ispec)
    c45 = c45store(INDEX_IJK,ispec)
    c46 = c46store(INDEX_IJK,ispec)
    c55 = c55store(INDEX_IJK,ispec)
    c56 = c56store(INDEX_IJK,ispec)
    c66 = c66store(INDEX_IJK,ispec)

    ! general expression of stress tensor for full Cijkl with 21 coefficients
    sigma_xx(INDEX_IJK) = c11*duxdxl(INDEX_IJK) + c16*duxdyl_plus_duydxl(INDEX_IJK) + c12*duydyl(INDEX_IJK) + &
             c15*duzdxl_plus_duxdzl(INDEX_IJK) + c14*duzdyl_plus_duydzl(INDEX_IJK) + c13*duzdzl(INDEX_IJK)

    sigma_yy(INDEX_IJK) = c12*duxdxl(INDEX_IJK) + c26*duxdyl_plus_duydxl(INDEX_IJK) + c22*duydyl(INDEX_IJK) + &
             c25*duzdxl_plus_duxdzl(INDEX_IJK) + c24*duzdyl_plus_duydzl(INDEX_IJK) + c23*duzdzl(INDEX_IJK)

    sigma_zz(INDEX_IJK) = c13*duxdxl(INDEX_IJK) + c36*duxdyl_plus_duydxl(INDEX_IJK) + c23*duydyl(INDEX_IJK) + &
             c35*duzdxl_plus_duxdzl(INDEX_IJK) + c34*duzdyl_plus_duydzl(INDEX_IJK) + c33*duzdzl(INDEX_IJK)

    sigma_xy(INDEX_IJK) = c16*duxdxl(INDEX_IJK) + c66*duxdyl_plus_duydxl(INDEX_IJK) + c26*duydyl(INDEX_IJK) + &
             c56*duzdxl_plus_duxdzl(INDEX_IJK) + c46*duzdyl_plus_duydzl(INDEX_IJK) + c36*duzdzl(INDEX_IJK)

    sigma_xz(INDEX_IJK) = c15*duxdxl(INDEX_IJK) + c56*duxdyl_plus_duydxl(INDEX_IJK) + c25*duydyl(INDEX_IJK) + &
             c55*duzdxl_plus_duxdzl(INDEX_IJK) + c45*duzdyl_plus_duydzl(INDEX_IJK) + c35*duzdzl(INDEX_IJK)

    sigma_yz(INDEX_IJK) = c14*duxdxl(INDEX_IJK) + c46*duxdyl_plus_duydxl(INDEX_IJK) + c24*duydyl(INDEX_IJK) + &
             c45*duzdxl_plus_duxdzl(INDEX_IJK) + c44*duzdyl_plus_duydzl(INDEX_IJK) + c34*duzdzl(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! attenuation contribution to stress
  ! subtract memory variables if attenuation
  if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
    call compute_element_stress_attenuation_contrib(R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                                    R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
  endif

  ! define symmetric components of sigma (to be general in case of gravity)
  DO_LOOP_IJK
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! compute non-symmetric terms for gravity
  if (GRAVITY_VAL) then
    call compute_element_gravity(ispec,NSPEC,NGLOB,ibool,jacobianl,wgll_cube, &
                                 gravity_pre_store,gravity_H, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product of stress tensor with test vector, non-symmetric form
  call compute_element_dot_product_stress(deriv(:,:,:,:,ispec),jacobianl, &
                                                sigma_xx,sigma_yy,sigma_zz, &
                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)

  end subroutine compute_element_tiso

! left for reference: original routine...
!
!  subroutine compute_element_tiso(ispec, &
!                                  minus_gravity_table,density_table,minus_deriv_gravity_table, &
!                                  rstore, &
!                                  deriv, &
!                                  wgll_cube, &
!                                  kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
!                                  ibool, &
!                                  R_xx,R_yy,R_xy,R_xz,R_yz, &
!                                  epsilon_trace_over_3, &
!                                  one_minus_sum_beta,vnspec, &
!                                  tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!                                  dummyx_loc,dummyy_loc,dummyz_loc, &
!                                  epsilondev_loc,rho_s_H)
!
!! tiso element in crust/mantle
!
!  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS,NRAD_GRAVITY,FOUR_THIRDS,TWO_THIRDS
!
!  use constants_solver, only: &
!    NSPEC => NSPEC_CRUST_MANTLE, &
!    NGLOB => NGLOB_CRUST_MANTLE, &
!    NSPECMAX_TISO => NSPECMAX_TISO_MANTLE, &
!    NSPEC_ATTENUATION => NSPEC_CRUST_MANTLE_ATTENUATION, &
!    NSPEC_STRAIN_ONLY => NSPEC_CRUST_MANTLE_STRAIN_ONLY, &
!    ATT1_VAL,ATT2_VAL,ATT3_VAL, &
!    ATTENUATION_VAL,ATTENUATION_3D_VAL,ATTENUATION_1D_WITH_3D_STORAGE_VAL, &
!    PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL
!
!  use specfem_par, only: COMPUTE_AND_STORE_STRAIN
!
!#ifdef FORCE_VECTORIZATION
!  use constants, only: NGLLCUBE
!#endif
!
!  implicit none
!
!  ! element id
!  integer,intent(in) :: ispec
!
!  ! arrays with mesh parameters per slice
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
!
!  ! x y and z contain r theta and phi
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB),intent(in) :: rstore
!
!  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: deriv
!
!  ! array with derivatives of Lagrange polynomials and precalculated products
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube
!
!  ! store anisotropic properties only where needed to save memory
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO),intent(in) :: &
!        kappahstore,muhstore,eta_anisostore
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: kappavstore,muvstore
!
!  ! attenuation
!  ! memory variables for attenuation
!  ! memory variables R_ij are stored at the local rather than global level
!  ! to allow for optimization of cache access by compiler
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: &
!    R_xx,R_yy,R_xy,R_xz,R_yz
!
!  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(out) :: epsilon_trace_over_3
!
!  ! variable sized array variables
!  integer,intent(in) :: vnspec
!
!  ! [alpha,beta,gamma]val reduced to N_SLS  to N_SLS*NUM_NODES
!  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec),intent(in) :: one_minus_sum_beta
!
!  ! gravity
!  double precision, dimension(NRAD_GRAVITY),intent(in) :: minus_gravity_table,density_table,minus_deriv_gravity_table
!
!  ! element info
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
!    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc
!
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: rho_s_H
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(out) :: epsilondev_loc
!
!! local parameters
!  real(kind=CUSTOM_REAL) :: one_minus_sum_beta_use
!  ! the 21 coefficients for an anisotropic medium in reduced notation
!  real(kind=CUSTOM_REAL) :: c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56
!
!  real(kind=CUSTOM_REAL) :: rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
!        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
!        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
!        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta
!
!  real(kind=CUSTOM_REAL) :: two_rhovsvsq,two_rhovshsq ! two_rhovpvsq,two_rhovphsq
!  real(kind=CUSTOM_REAL) :: four_rhovsvsq,four_rhovshsq ! four_rhovpvsq,four_rhovphsq
!
!  real(kind=CUSTOM_REAL) :: twoetaminone,etaminone,eta_aniso
!  real(kind=CUSTOM_REAL) :: two_eta_aniso,four_eta_aniso,six_eta_aniso
!
!  ! local element arrays
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl, duydyl, duzdzl
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
!
!  real(kind=CUSTOM_REAL) :: templ1,templ1_cos,templ2,templ2_cos,templ3,templ3_two,templ3_cos
!  real(kind=CUSTOM_REAL) :: kappavl,kappahl,muvl,muhl
!
!  integer :: iglob
!
!#ifdef FORCE_VECTORIZATION
!! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
!! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
!! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
!  integer :: ijk
!#else
!  integer :: i,j,k
!#endif
!! note: profiling shows that this routine takes about 30% of the total time, another 60% is spend in the iso routine above..
!
!  ! transverse isotropic element
!
!  ! precomputes factors
!  call compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!                                          deriv(:,:,:,:,ispec),jacobianl, &
!                                          duxdxl,duydyl,duzdzl, &
!                                          duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
!                                          duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)
!
!  ! compute deviatoric strain
!  if (COMPUTE_AND_STORE_STRAIN) then
!    call compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
!                                           duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
!                                           ispec,NSPEC_STRAIN_ONLY, &
!                                           epsilon_trace_over_3,epsilondev_loc)
!  endif
!
!
!  !
!  ! compute either transversely isotropic elements
!  !
!  DO_LOOP_IJK
!
!! note: the mesh is built such that anisotropic elements are created first in anisotropic layers,
!!           thus they are listed first ( see in create_regions_mesh.f90: perm_layer() ordering )
!!           this is therefore still in bounds of 1:NSPECMAX_TISO_MANTLE even if NSPECMAX_TISO is less than NSPEC
!
!    ! use kappa and mu from transversely isotropic model
!    kappavl = kappavstore(INDEX_IJK,ispec)
!    muvl = muvstore(INDEX_IJK,ispec)
!
!    kappahl = kappahstore(INDEX_IJK,ispec)
!    muhl = muhstore(INDEX_IJK,ispec)
!
!    ! use unrelaxed parameters if attenuation
!    ! eta does not need to be shifted since it is a ratio
!    if (ATTENUATION_VAL) then
!      ! precompute terms for attenuation if needed
!      if (ATTENUATION_3D_VAL .or. ATTENUATION_1D_WITH_3D_STORAGE_VAL) then
!        one_minus_sum_beta_use = one_minus_sum_beta(INDEX_IJK,ispec)
!      else
!        one_minus_sum_beta_use = one_minus_sum_beta(1,1,1,ispec)
!      endif
!      muvl = muvl * one_minus_sum_beta_use
!      muhl = muhl * one_minus_sum_beta_use
!    endif
!
!    rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
!    rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A
!
!    rhovsvsq = muvl  !!! that is L
!    rhovshsq = muhl  !!! that is N
!
!    eta_aniso = eta_anisostore(INDEX_IJK,ispec)  !!! that is  F / (A - 2 L)
!
!    ! use mesh coordinates to get theta and phi
!    ! rstore contains theta and phi
!    iglob = ibool(INDEX_IJK,ispec)
!
!    theta = rstore(2,iglob)
!    phi = rstore(3,iglob)
!
!     ! precompute some products to reduce the CPU time
!
!    costheta = cos(theta)
!    sintheta = sin(theta)
!    cosphi = cos(phi)
!    sinphi = sin(phi)
!
!    costhetasq = costheta * costheta
!    sinthetasq = sintheta * sintheta
!    cosphisq = cosphi * cosphi
!    sinphisq = sinphi * sinphi
!
!    costhetafour = costhetasq * costhetasq
!    sinthetafour = sinthetasq * sinthetasq
!    cosphifour = cosphisq * cosphisq
!    sinphifour = sinphisq * sinphisq
!
!    costwotheta = cos(2.0_CUSTOM_REAL*theta)
!    sintwotheta = sin(2.0_CUSTOM_REAL*theta)
!    costwophi = cos(2.0_CUSTOM_REAL*phi)
!    sintwophi = sin(2.0_CUSTOM_REAL*phi)
!
!    cosfourtheta = cos(4.0_CUSTOM_REAL*theta)
!    cosfourphi = cos(4.0_CUSTOM_REAL*phi)
!
!    costwothetasq = costwotheta * costwotheta
!
!    costwophisq = costwophi * costwophi
!    sintwophisq = sintwophi * sintwophi
!
!    etaminone = eta_aniso - 1.0_CUSTOM_REAL
!    twoetaminone = 2.0_CUSTOM_REAL * eta_aniso - 1.0_CUSTOM_REAL
!
!    ! precompute some products to reduce the CPU time
!    two_eta_aniso = 2.0_CUSTOM_REAL*eta_aniso
!    four_eta_aniso = 4.0_CUSTOM_REAL*eta_aniso
!    six_eta_aniso = 6.0_CUSTOM_REAL*eta_aniso
!
!    two_rhovsvsq = 2.0_CUSTOM_REAL*rhovsvsq
!    two_rhovshsq = 2.0_CUSTOM_REAL*rhovshsq
!    four_rhovsvsq = 4.0_CUSTOM_REAL*rhovsvsq
!    four_rhovshsq = 4.0_CUSTOM_REAL*rhovshsq
!
!    ! pre-compute temporary values
!    templ1 = four_rhovsvsq - rhovpvsq + twoetaminone*rhovphsq - four_eta_aniso*rhovsvsq
!    templ1_cos = rhovphsq - rhovpvsq + costwotheta*templ1
!    templ2 = four_rhovsvsq - rhovpvsq - rhovphsq + two_eta_aniso*rhovphsq - four_eta_aniso*rhovsvsq
!    templ2_cos = rhovpvsq - rhovphsq + costwotheta*templ2
!    templ3 = rhovphsq + rhovpvsq - two_eta_aniso*rhovphsq + four_eta_aniso*rhovsvsq
!    templ3_two = templ3 - two_rhovshsq - two_rhovsvsq
!    templ3_cos = templ3_two + costwotheta*templ2
!
!    ! reordering operations to facilitate compilation, avoiding divisions, using locality for temporary values
!    c11 = rhovphsq*sinphifour &
!          + 2.0_CUSTOM_REAL*cosphisq*sinphisq* &
!          ( rhovphsq*costhetasq + sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq) ) &
!          + cosphifour*(rhovphsq*costhetafour &
!            + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq) &
!            + rhovpvsq*sinthetafour)
!
!    c12 = 0.25_CUSTOM_REAL*costhetasq*(rhovphsq - two_rhovshsq)*(3.0_CUSTOM_REAL + cosfourphi) &
!          - four_rhovshsq*cosphisq*costhetasq*sinphisq &
!          + 0.03125_CUSTOM_REAL*rhovphsq*sintwophisq*(11.0_CUSTOM_REAL + cosfourtheta + 4.0*costwotheta) &
!          + eta_aniso*sinthetasq*(rhovphsq - two_rhovsvsq) &
!                     *(cosphifour + sinphifour + 2.0_CUSTOM_REAL*cosphisq*costhetasq*sinphisq) &
!          + rhovpvsq*cosphisq*sinphisq*sinthetafour &
!          - rhovsvsq*sintwophisq*sinthetafour
!
!    c13 = 0.125_CUSTOM_REAL*cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq &
!                - 12.0_CUSTOM_REAL*eta_aniso*rhovsvsq + cosfourtheta*templ1) &
!          + sinphisq*(eta_aniso*costhetasq*(rhovphsq - two_rhovsvsq) + sinthetasq*(rhovphsq - two_rhovshsq))
!
!    ! uses temporary templ1 from c13
!    c15 = cosphi*costheta*sintheta* &
!          ( 0.5_CUSTOM_REAL*cosphisq* (rhovpvsq - rhovphsq + costwotheta*templ1) &
!            + etaminone*sinphisq*(rhovphsq - two_rhovsvsq))
!
!    c14 = costheta*sinphi*sintheta* &
!          ( 0.5_CUSTOM_REAL*cosphisq*(templ2_cos + four_rhovshsq - four_rhovsvsq) &
!            + sinphisq*(etaminone*rhovphsq + 2.0_CUSTOM_REAL*(rhovshsq - eta_aniso*rhovsvsq)) )
!
!    ! uses temporary templ2_cos from c14
!    c16 = 0.5_CUSTOM_REAL*cosphi*sinphi*sinthetasq* &
!          ( cosphisq*templ2_cos &
!            + 2.0_CUSTOM_REAL*etaminone*sinphisq*(rhovphsq - two_rhovsvsq) )
!
!    c22 = rhovphsq*cosphifour + 2.0_CUSTOM_REAL*cosphisq*sinphisq* &
!          (rhovphsq*costhetasq + sinthetasq*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)) &
!          + sinphifour* &
!          (rhovphsq*costhetafour + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(eta_aniso*rhovphsq &
!                + two_rhovsvsq - two_eta_aniso*rhovsvsq) + rhovpvsq*sinthetafour)
!
!    ! uses temporary templ1 from c13
!    c23 = 0.125_CUSTOM_REAL*sinphisq*(rhovphsq + six_eta_aniso*rhovphsq &
!            + rhovpvsq - four_rhovsvsq - 12.0_CUSTOM_REAL*eta_aniso*rhovsvsq + cosfourtheta*templ1) &
!          + cosphisq*(eta_aniso*costhetasq*(rhovphsq - two_rhovsvsq) + sinthetasq*(rhovphsq - two_rhovshsq))
!
!    ! uses temporary templ1 from c13
!    c24 = costheta*sinphi*sintheta* &
!          ( etaminone*cosphisq*(rhovphsq - two_rhovsvsq) &
!            + 0.5_CUSTOM_REAL*sinphisq*(rhovpvsq - rhovphsq + costwotheta*templ1) )
!
!    ! uses temporary templ2_cos from c14
!    c25 = cosphi*costheta*sintheta* &
!          ( cosphisq*(etaminone*rhovphsq + 2.0_CUSTOM_REAL*(rhovshsq - eta_aniso*rhovsvsq)) &
!            + 0.5_CUSTOM_REAL*sinphisq*(templ2_cos + four_rhovshsq - four_rhovsvsq) )
!
!    ! uses temporary templ2_cos from c14
!    c26 = 0.5_CUSTOM_REAL*cosphi*sinphi*sinthetasq* &
!          ( 2.0_CUSTOM_REAL*etaminone*cosphisq*(rhovphsq - two_rhovsvsq) &
!            + sinphisq*templ2_cos )
!
!    c33 = rhovpvsq*costhetafour &
!          + 2.0_CUSTOM_REAL*costhetasq*sinthetasq*(two_rhovsvsq + eta_aniso*(rhovphsq - two_rhovsvsq)) &
!          + rhovphsq*sinthetafour
!
!    ! uses temporary templ1_cos from c13
!    c34 = - 0.25_CUSTOM_REAL*sinphi*sintwotheta*templ1_cos
!
!    ! uses temporary templ1_cos from c34
!    c35 = - 0.25_CUSTOM_REAL*cosphi*sintwotheta*templ1_cos
!
!    ! uses temporary templ1_cos from c34
!    c36 = - 0.25_CUSTOM_REAL*sintwophi*sinthetasq*(templ1_cos - four_rhovshsq + four_rhovsvsq)
!
!    c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) &
!          + sinphisq*(rhovsvsq*costwothetasq + costhetasq*sinthetasq*templ3)
!
!    ! uses temporary templ3 from c44
!    c46 = - cosphi*costheta*sintheta* &
!            ( cosphisq*(rhovshsq - rhovsvsq) - 0.5_CUSTOM_REAL*sinphisq*templ3_cos  )
!
!    ! uses templ3 from c46
!    c45 = 0.25_CUSTOM_REAL*sintwophi*sinthetasq* &
!          (templ3_two + costwotheta*(rhovphsq + rhovpvsq - two_eta_aniso*rhovphsq + 4.0_CUSTOM_REAL*etaminone*rhovsvsq))
!
!    c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) &
!          + cosphisq*(rhovsvsq*costwothetasq &
!              + costhetasq*sinthetasq*(rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq) )
!
!    ! uses temporary templ3_cos from c46
!    c56 = costheta*sinphi*sintheta* &
!          ( 0.5_CUSTOM_REAL*cosphisq*templ3_cos + sinphisq*(rhovsvsq - rhovshsq) )
!
!    c66 = rhovshsq*costwophisq*costhetasq &
!          - 2.0_CUSTOM_REAL*cosphisq*costhetasq*sinphisq*(rhovphsq - two_rhovshsq) &
!          + 0.03125_CUSTOM_REAL*rhovphsq*sintwophisq*(11.0_CUSTOM_REAL + 4.0_CUSTOM_REAL*costwotheta + cosfourtheta) &
!          - 0.125_CUSTOM_REAL*rhovsvsq*sinthetasq* &
!          ( -6.0_CUSTOM_REAL - 2.0_CUSTOM_REAL*costwotheta - 2.0_CUSTOM_REAL*cosfourphi &
!                    + cos(4.0_CUSTOM_REAL*phi - 2.0_CUSTOM_REAL*theta) &
!                    + cos(2.0_CUSTOM_REAL*(2.0_CUSTOM_REAL*phi + theta)) ) &
!          + rhovpvsq*cosphisq*sinphisq*sinthetafour &
!          - 0.5_CUSTOM_REAL*eta_aniso*sintwophisq*sinthetafour*(rhovphsq - two_rhovsvsq)
!
!    ! general expression of stress tensor for full Cijkl with 21 coefficients
!    sigma_xx(INDEX_IJK) = c11*duxdxl(INDEX_IJK) + c16*duxdyl_plus_duydxl(INDEX_IJK) + c12*duydyl(INDEX_IJK) + &
!             c15*duzdxl_plus_duxdzl(INDEX_IJK) + c14*duzdyl_plus_duydzl(INDEX_IJK) + c13*duzdzl(INDEX_IJK)
!
!    sigma_yy(INDEX_IJK) = c12*duxdxl(INDEX_IJK) + c26*duxdyl_plus_duydxl(INDEX_IJK) + c22*duydyl(INDEX_IJK) + &
!             c25*duzdxl_plus_duxdzl(INDEX_IJK) + c24*duzdyl_plus_duydzl(INDEX_IJK) + c23*duzdzl(INDEX_IJK)
!
!    sigma_zz(INDEX_IJK) = c13*duxdxl(INDEX_IJK) + c36*duxdyl_plus_duydxl(INDEX_IJK) + c23*duydyl(INDEX_IJK) + &
!             c35*duzdxl_plus_duxdzl(INDEX_IJK) + c34*duzdyl_plus_duydzl(INDEX_IJK) + c33*duzdzl(INDEX_IJK)
!
!    sigma_xy(INDEX_IJK) = c16*duxdxl(INDEX_IJK) + c66*duxdyl_plus_duydxl(INDEX_IJK) + c26*duydyl(INDEX_IJK) + &
!             c56*duzdxl_plus_duxdzl(INDEX_IJK) + c46*duzdyl_plus_duydzl(INDEX_IJK) + c36*duzdzl(INDEX_IJK)
!
!    sigma_xz(INDEX_IJK) = c15*duxdxl(INDEX_IJK) + c56*duxdyl_plus_duydxl(INDEX_IJK) + c25*duydyl(INDEX_IJK) + &
!             c55*duzdxl_plus_duxdzl(INDEX_IJK) + c45*duzdyl_plus_duydzl(INDEX_IJK) + c35*duzdzl(INDEX_IJK)
!
!    sigma_yz(INDEX_IJK) = c14*duxdxl(INDEX_IJK) + c46*duxdyl_plus_duydxl(INDEX_IJK) + c24*duydyl(INDEX_IJK) + &
!             c45*duzdxl_plus_duxdzl(INDEX_IJK) + c44*duzdyl_plus_duydzl(INDEX_IJK) + c34*duzdzl(INDEX_IJK)
!
!  ENDDO_LOOP_IJK
!
!  ! attenuation contribution to stress
!  ! subtract memory variables if attenuation
!  if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
!    call compute_element_stress_attenuation_contrib(R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
!                                                    R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
!                                                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
!  endif
!
!  ! define symmetric components of sigma (to be general in case of gravity)
!  DO_LOOP_IJK
!
!    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
!    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
!    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
!
!  ENDDO_LOOP_IJK
!
!  ! compute non-symmetric terms for gravity
!  if (GRAVITY_VAL) then
!    call compute_element_gravity(ispec,NSPEC,NGLOB,ibool,rstore,jacobianl,wgll_cube, &
!                                 minus_gravity_table,minus_deriv_gravity_table,density_table, &
!                                 dummyx_loc,dummyy_loc,dummyz_loc, &
!                                 sigma_xx,sigma_yy,sigma_zz, &
!                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
!                                 rho_s_H)
!  endif
!
!  ! dot product of stress tensor with test vector, non-symmetric form
!  call compute_element_dot_product_stress(deriv(:,:,:,:,ispec),jacobianl, &
!                                                sigma_xx,sigma_yy,sigma_zz, &
!                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
!                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)
!
!  end subroutine compute_element_tiso
!


!--------------------------------------------------------------------------------------------
!
! anisotropic element
!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_aniso(ispec, &
                                   gravity_pre_store,gravity_H, &
                                   deriv, &
                                   wgll_cube, &
                                   c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                   c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                   c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                   ibool, &
                                   R_xx,R_yy,R_xy,R_xz,R_yz, &
                                   epsilon_trace_over_3, &
                                   tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                   dummyx_loc,dummyy_loc,dummyz_loc, &
                                   epsilondev_loc,rho_s_H)

! fully anisotropic element in crust/mantle region

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS

  use constants_solver, only: &
    NSPEC => NSPEC_CRUST_MANTLE, &
    NGLOB => NGLOB_CRUST_MANTLE, &
    NSPECMAX_ANISO => NSPECMAX_ANISO_MANTLE, &
    NSPEC_ATTENUATION => NSPEC_CRUST_MANTLE_ATTENUATION, &
    NSPEC_STRAIN_ONLY => NSPEC_CRUST_MANTLE_STRAIN_ONLY, &
    ATTENUATION_VAL, &
    PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL

  use specfem_par, only: COMPUTE_AND_STORE_STRAIN

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: deriv

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO),intent(in) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(in) :: gravity_pre_store
  real(kind=CUSTOM_REAL),dimension(6,NGLOB),intent(in) :: gravity_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(inout) :: epsilondev_loc

  ! local parameters
  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) :: c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  ! local element arrays
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl, duydyl, duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  !  anisotropic elements

  ! precomputes factors
  call compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                          deriv(:,:,:,:,ispec),jacobianl, &
                                          duxdxl,duydyl,duzdzl, &
                                          duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                                          duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)

  ! compute deviatoric strain
  if (COMPUTE_AND_STORE_STRAIN) then
    call compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
                                           duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
                                           ispec,NSPEC_STRAIN_ONLY, &
                                           epsilon_trace_over_3,epsilondev_loc)
  endif

  !
  ! compute anisotropic elements
  !
  DO_LOOP_IJK
    c11 = c11store(INDEX_IJK,ispec)
    c12 = c12store(INDEX_IJK,ispec)
    c13 = c13store(INDEX_IJK,ispec)
    c14 = c14store(INDEX_IJK,ispec)
    c15 = c15store(INDEX_IJK,ispec)
    c16 = c16store(INDEX_IJK,ispec)
    c22 = c22store(INDEX_IJK,ispec)
    c23 = c23store(INDEX_IJK,ispec)
    c24 = c24store(INDEX_IJK,ispec)
    c25 = c25store(INDEX_IJK,ispec)
    c26 = c26store(INDEX_IJK,ispec)
    c33 = c33store(INDEX_IJK,ispec)
    c34 = c34store(INDEX_IJK,ispec)
    c35 = c35store(INDEX_IJK,ispec)
    c36 = c36store(INDEX_IJK,ispec)
    c44 = c44store(INDEX_IJK,ispec)
    c45 = c45store(INDEX_IJK,ispec)
    c46 = c46store(INDEX_IJK,ispec)
    c55 = c55store(INDEX_IJK,ispec)
    c56 = c56store(INDEX_IJK,ispec)
    c66 = c66store(INDEX_IJK,ispec)

    sigma_xx(INDEX_IJK) = c11*duxdxl(INDEX_IJK) + c16*duxdyl_plus_duydxl(INDEX_IJK) + c12*duydyl(INDEX_IJK) + &
             c15*duzdxl_plus_duxdzl(INDEX_IJK) + c14*duzdyl_plus_duydzl(INDEX_IJK) + c13*duzdzl(INDEX_IJK)

    sigma_yy(INDEX_IJK) = c12*duxdxl(INDEX_IJK) + c26*duxdyl_plus_duydxl(INDEX_IJK) + c22*duydyl(INDEX_IJK) + &
             c25*duzdxl_plus_duxdzl(INDEX_IJK) + c24*duzdyl_plus_duydzl(INDEX_IJK) + c23*duzdzl(INDEX_IJK)

    sigma_zz(INDEX_IJK) = c13*duxdxl(INDEX_IJK) + c36*duxdyl_plus_duydxl(INDEX_IJK) + c23*duydyl(INDEX_IJK) + &
             c35*duzdxl_plus_duxdzl(INDEX_IJK) + c34*duzdyl_plus_duydzl(INDEX_IJK) + c33*duzdzl(INDEX_IJK)

    sigma_xy(INDEX_IJK) = c16*duxdxl(INDEX_IJK) + c66*duxdyl_plus_duydxl(INDEX_IJK) + c26*duydyl(INDEX_IJK) + &
             c56*duzdxl_plus_duxdzl(INDEX_IJK) + c46*duzdyl_plus_duydzl(INDEX_IJK) + c36*duzdzl(INDEX_IJK)

    sigma_xz(INDEX_IJK) = c15*duxdxl(INDEX_IJK) + c56*duxdyl_plus_duydxl(INDEX_IJK) + c25*duydyl(INDEX_IJK) + &
             c55*duzdxl_plus_duxdzl(INDEX_IJK) + c45*duzdyl_plus_duydzl(INDEX_IJK) + c35*duzdzl(INDEX_IJK)

    sigma_yz(INDEX_IJK) = c14*duxdxl(INDEX_IJK) + c46*duxdyl_plus_duydxl(INDEX_IJK) + c24*duydyl(INDEX_IJK) + &
             c45*duzdxl_plus_duxdzl(INDEX_IJK) + c44*duzdyl_plus_duydzl(INDEX_IJK) + c34*duzdzl(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! attenuation contribution to stress
  ! subtract memory variables if attenuation
  if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
    call compute_element_stress_attenuation_contrib(R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                                    R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
  endif

  ! define symmetric components of sigma (to be general in case of gravity)
  DO_LOOP_IJK
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! compute non-symmetric terms for gravity
  if (GRAVITY_VAL) then
    call compute_element_gravity(ispec,NSPEC,NGLOB,ibool,jacobianl,wgll_cube, &
                                 gravity_pre_store,gravity_H, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product of stress tensor with test vector, non-symmetric form
  call compute_element_dot_product_stress(deriv(:,:,:,:,ispec),jacobianl, &
                                                sigma_xx,sigma_yy,sigma_zz, &
                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)

  end subroutine compute_element_aniso



!
!--------------------------------------------------------------------------------------------

  subroutine compute_element_aniso_ic(ispec, &
                                      gravity_pre_store,gravity_H, &
                                      deriv, &
                                      wgll_cube, &
                                      c11store,c12store,c13store,c33store,c44store, &
                                      ibool, &
                                      R_xx,R_yy,R_xy,R_xz,R_yz, &
                                      epsilon_trace_over_3, &
                                      tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                      dummyx_loc,dummyy_loc,dummyz_loc, &
                                      epsilondev_loc,rho_s_H)

! anisotropic element in inner core with hexagonal symmetry (and vertical symmetry axis)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,N_SLS

  use constants_solver, only: &
    NSPEC => NSPEC_INNER_CORE, &
    NGLOB => NGLOB_INNER_CORE, &
    NSPECMAX_ANISO => NSPECMAX_ANISO_IC, &
    NSPEC_ATTENUATION => NSPEC_INNER_CORE_ATTENUATION, &
    NSPEC_STRAIN_ONLY => NSPEC_INNER_CORE_STRAIN_ONLY, &
    ATTENUATION_VAL, &
    PARTIAL_PHYS_DISPERSION_ONLY_VAL,GRAVITY_VAL

  use specfem_par, only: COMPUTE_AND_STORE_STRAIN

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! element id
  integer,intent(in) :: ispec

  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: deriv

  ! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO),intent(in) :: &
        c11store,c12store,c13store,c33store,c44store

  !real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: muvstore

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_ATTENUATION),intent(in) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(in) :: gravity_pre_store
  real(kind=CUSTOM_REAL),dimension(6,NGLOB),intent(in) :: gravity_H

  ! element info
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(inout) :: epsilondev_loc

  ! local parameters
  real(kind=CUSTOM_REAL) :: c11,c12,c13,c22,c23,c33,c44,c55,c66

  ! local element arrays
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl, duydyl, duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#else
  integer :: i,j,k
#endif

  !  anisotropic elements

  ! precomputes factors
  call compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                          deriv(:,:,:,:,ispec),jacobianl, &
                                          duxdxl,duydyl,duzdzl, &
                                          duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                                          duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)

  ! compute deviatoric strain
  if (COMPUTE_AND_STORE_STRAIN) then
    call compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
                                           duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
                                           ispec,NSPEC_STRAIN_ONLY, &
                                           epsilon_trace_over_3,epsilondev_loc)
  endif

  !
  ! compute anisotropic elements
  !
  DO_LOOP_IJK
    ! elastic tensor for hexagonal symmetry in reduced notation:
    !
    !      c11 c12 c13  0   0        0
    !      c12 c11 c13  0   0        0
    !      c13 c13 c33  0   0        0
    !       0   0   0  c44  0        0
    !       0   0   0   0  c44       0
    !       0   0   0   0   0  (c11-c12)/2
    !
    !       in terms of the A, C, L, N and F of Love (1927):
    !
    !       c11 = A
    !       c12 = A-2N
    !       c13 = F
    !       c33 = C
    !       c44 = L
    !
    ! transversely isotropic
    !       C11 = A = rho * vph**2
    !       C33 = C = rho * vpv**2
    !       C44 = L = rho * vsv**2
    !       C13 = F = eta * (A - 2*L)
    !       C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
    !       C22 = C11 = A
    !       C23 = C13 = F
    !       C55 = C44 = L
    !       C66 = N = rho * vsh**2 = (C11-C12)/2
    !
    ! isotropic
    ! Lame parameters: mu = rho * vs**2
    !                  lambda = rho * (vp**2 - 2 vs**2) = rho * vp**2 - 2 mu
    !
    !            then: C11 = C22 = C33 = lambda + 2mu
    !                  C12 = C13 = C23 = lambda
    !                  C44 = C55 = C66 = mu
    c11 = c11store(INDEX_IJK,ispec)
    c12 = c12store(INDEX_IJK,ispec)
    c13 = c13store(INDEX_IJK,ispec)
    c33 = c33store(INDEX_IJK,ispec)
    c44 = c44store(INDEX_IJK,ispec)
    ! tiso symmetry
    c22 = c11
    c23 = c13
    c55 = c44
    c66 = 0.5_CUSTOM_REAL*(c11-c12)

    sigma_xx(INDEX_IJK) = c11*duxdxl(INDEX_IJK) + c12*duydyl(INDEX_IJK) + c13*duzdzl(INDEX_IJK)
    sigma_yy(INDEX_IJK) = c12*duxdxl(INDEX_IJK) + c22*duydyl(INDEX_IJK) + c23*duzdzl(INDEX_IJK)
    sigma_zz(INDEX_IJK) = c13*duxdxl(INDEX_IJK) + c23*duydyl(INDEX_IJK) + c33*duzdzl(INDEX_IJK)

    sigma_xy(INDEX_IJK) = c66*duxdyl_plus_duydxl(INDEX_IJK)
    sigma_xz(INDEX_IJK) = c55*duzdxl_plus_duxdzl(INDEX_IJK)
    sigma_yz(INDEX_IJK) = c44*duzdyl_plus_duydzl(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! attenuation contribution to stress
  ! subtract memory variables if attenuation
  if (ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL) then
    call compute_element_stress_attenuation_contrib(R_xx(1,1,1,1,ispec),R_yy(1,1,1,1,ispec),R_xy(1,1,1,1,ispec), &
                                                    R_xz(1,1,1,1,ispec),R_yz(1,1,1,1,ispec), &
                                                    sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)
  endif

  ! define symmetric components of sigma (to be general in case of gravity)
  DO_LOOP_IJK
    sigma_yx(INDEX_IJK) = sigma_xy(INDEX_IJK)
    sigma_zx(INDEX_IJK) = sigma_xz(INDEX_IJK)
    sigma_zy(INDEX_IJK) = sigma_yz(INDEX_IJK)
  ENDDO_LOOP_IJK

  ! compute non-symmetric terms for gravity
  if (GRAVITY_VAL) then
    call compute_element_gravity(ispec,NSPEC,NGLOB,ibool,jacobianl,wgll_cube, &
                                 gravity_pre_store,gravity_H, &
                                 dummyx_loc,dummyy_loc,dummyz_loc, &
                                 sigma_xx,sigma_yy,sigma_zz, &
                                 sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                 rho_s_H)
  endif

  ! dot product of stress tensor with test vector, non-symmetric form
  call compute_element_dot_product_stress(deriv(:,:,:,:,ispec),jacobianl, &
                                                sigma_xx,sigma_yy,sigma_zz, &
                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)

  end subroutine compute_element_aniso_ic



!--------------------------------------------------------------------------------------------
!
! helper functions
!
!--------------------------------------------------------------------------------------------

!
! please leave this routine in this file, to help compilers inlining this function...
!

  subroutine compute_element_stress_attenuation_contrib(R_xx_loc,R_yy_loc,R_xy_loc,R_xz_loc,R_yz_loc, &
                                                        sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES INLINE :: compute_element_stress_attenuation_contrib
#else
! cray
!DIR$ INLINEALWAYS compute_element_stress_attenuation_contrib
#endif

! updates stress with attenuation coontribution

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,N_SLS

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS),intent(in) :: R_xx_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS),intent(in) :: R_yy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS),intent(in) :: R_xy_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS),intent(in) :: R_xz_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,N_SLS),intent(in) :: R_yz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xy,sigma_xz,sigma_yz

  ! local parameters
  real(kind=CUSTOM_REAL) :: R_xx_val,R_yy_val

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
  integer :: i_SLS
#endif

#ifdef FORCE_VECTORIZATION

  ! here we assume that N_SLS == 3 in order to be able to unroll and suppress the loop
  ! in order to vectorize the outer loop
  DO_LOOP_IJK
    R_xx_val = R_xx_loc(INDEX_IJK,1)
    R_yy_val = R_yy_loc(INDEX_IJK,1)
    sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) - R_xx_val
    sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) - R_yy_val
    sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + R_xx_val + R_yy_val
    sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - R_xy_loc(INDEX_IJK,1)
    sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - R_xz_loc(INDEX_IJK,1)
    sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - R_yz_loc(INDEX_IJK,1)

    R_xx_val = R_xx_loc(INDEX_IJK,2)
    R_yy_val = R_yy_loc(INDEX_IJK,2)
    sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) - R_xx_val
    sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) - R_yy_val
    sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + R_xx_val + R_yy_val
    sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - R_xy_loc(INDEX_IJK,2)
    sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - R_xz_loc(INDEX_IJK,2)
    sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - R_yz_loc(INDEX_IJK,2)

    R_xx_val = R_xx_loc(INDEX_IJK,3)
    R_yy_val = R_yy_loc(INDEX_IJK,3)
    sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) - R_xx_val
    sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) - R_yy_val
    sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + R_xx_val + R_yy_val
    sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - R_xy_loc(INDEX_IJK,3)
    sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - R_xz_loc(INDEX_IJK,3)
    sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - R_yz_loc(INDEX_IJK,3)
  ENDDO_LOOP_IJK

#else

  ! loops over standard linear solids
  do i_SLS = 1,N_SLS
    DO_LOOP_IJK
      R_xx_val = R_xx_loc(INDEX_IJK,i_SLS)
      R_yy_val = R_yy_loc(INDEX_IJK,i_SLS)
      sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) - R_xx_val
      sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) - R_yy_val
      sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + R_xx_val + R_yy_val
      sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - R_xy_loc(INDEX_IJK,i_SLS)
      sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - R_xz_loc(INDEX_IJK,i_SLS)
      sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - R_yz_loc(INDEX_IJK,i_SLS)
    ENDDO_LOOP_IJK
  enddo

#endif

  end subroutine compute_element_stress_attenuation_contrib

!
!--------------------------------------------------------------------------------------------
!

! please leave this routine in this file, to help compilers inlining this function...

  subroutine compute_element_precompute_factors(tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                                                deriv_loc,jacobianl, &
                                                duxdxl,duydyl,duzdzl, &
                                                duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
                                                duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES INLINE :: compute_element_precompute_factors
#else
! cray
!DIR$ INLINEALWAYS compute_element_precompute_factors
#endif

! precomputes factors

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: tempz1,tempz2,tempz3

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ),intent(in) :: deriv_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: jacobianl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: duxdxl,duydyl,duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl, &
    duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  ! local parameters
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobian
  real(kind=CUSTOM_REAL) :: duxdyl,duxdzl,duydxl,duydzl,duzdxl,duzdyl
  real(kind=CUSTOM_REAL) :: x1,x2,x3,y1,y2,y3,z1,z2,z3

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! precomputes factors
  DO_LOOP_IJK
    ! get derivatives of ux, uy and uz with respect to x, y and z
    xixl = deriv_loc(1,INDEX_IJK)
    xiyl = deriv_loc(2,INDEX_IJK)
    xizl = deriv_loc(3,INDEX_IJK)
    etaxl = deriv_loc(4,INDEX_IJK)
    etayl = deriv_loc(5,INDEX_IJK)
    etazl = deriv_loc(6,INDEX_IJK)
    gammaxl = deriv_loc(7,INDEX_IJK)
    gammayl = deriv_loc(8,INDEX_IJK)
    gammazl = deriv_loc(9,INDEX_IJK)

    ! compute the Jacobian
    jacobian = (xixl*(etayl*gammazl-etazl*gammayl) &
              - xiyl*(etaxl*gammazl-etazl*gammaxl) &
              + xizl*(etaxl*gammayl-etayl*gammaxl))
    if (jacobian <= 0.0_CUSTOM_REAL) stop 'Error invalid jacobian in compute_element_precompute_factors()'

    jacobianl(INDEX_IJK) = 1.0_CUSTOM_REAL / jacobian

    x1 = tempx1(INDEX_IJK)
    x2 = tempx2(INDEX_IJK)
    x3 = tempx3(INDEX_IJK)

    y1 = tempy1(INDEX_IJK)
    y2 = tempy2(INDEX_IJK)
    y3 = tempy3(INDEX_IJK)

    z1 = tempz1(INDEX_IJK)
    z2 = tempz2(INDEX_IJK)
    z3 = tempz3(INDEX_IJK)

    duxdxl(INDEX_IJK) = xixl*x1 + etaxl*x2 + gammaxl*x3
    duxdyl            = xiyl*x1 + etayl*x2 + gammayl*x3
    duxdzl            = xizl*x1 + etazl*x2 + gammazl*x3

    duydxl            = xixl*y1 + etaxl*y2 + gammaxl*y3
    duydyl(INDEX_IJK) = xiyl*y1 + etayl*y2 + gammayl*y3
    duydzl            = xizl*y1 + etazl*y2 + gammazl*y3

    duzdxl            = xixl*z1 + etaxl*z2 + gammaxl*z3
    duzdyl            = xiyl*z1 + etayl*z2 + gammayl*z3
    duzdzl(INDEX_IJK) = xizl*z1 + etazl*z2 + gammazl*z3

    ! precompute some sums to save CPU time
    duxdxl_plus_duydyl(INDEX_IJK) = duxdxl(INDEX_IJK) + duydyl(INDEX_IJK)
    duxdxl_plus_duzdzl(INDEX_IJK) = duxdxl(INDEX_IJK) + duzdzl(INDEX_IJK)
    duydyl_plus_duzdzl(INDEX_IJK) = duydyl(INDEX_IJK) + duzdzl(INDEX_IJK)
    duxdyl_plus_duydxl(INDEX_IJK) = duxdyl            + duydxl
    duzdxl_plus_duxdzl(INDEX_IJK) = duzdxl            + duxdzl
    duzdyl_plus_duydzl(INDEX_IJK) = duzdyl            + duydzl
  ENDDO_LOOP_IJK

  end subroutine compute_element_precompute_factors

!
!--------------------------------------------------------------------------------------------
!

! please leave this routine in this file, to help compilers inlining this function...

  subroutine compute_element_deviatoric_strain(duxdxl,duydyl,duzdzl, &
                                               duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl, &
                                               ispec,NSPEC_STRAIN_ONLY, &
                                               epsilon_trace_over_3,epsilondev_loc)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES INLINE :: compute_element_deviatoric_strain
#else
! cray
!DIR$ INLINEALWAYS compute_element_deviatoric_strain
#endif

! computes deviatoric strain

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,ONE_THIRD

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: duxdxl,duydyl,duzdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  integer, intent(in) :: ispec,NSPEC_STRAIN_ONLY
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_STRAIN_ONLY),intent(inout) :: epsilon_trace_over_3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5),intent(inout) :: epsilondev_loc

  ! local parameters
  real(kind=CUSTOM_REAL) :: templ

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

!ZN beware, here the expression differs from the strain used in memory variable equation (6) in D. Komatitsch and J. Tromp 1999,
!ZN here Brian Savage uses the engineering strain which are epsilon = 1/2*(grad U + (grad U)^T),
!ZN where U is the displacement vector and grad the gradient operator, i.e. there is a 1/2 factor difference between the two.
!ZN Both expressions are fine, but we need to keep in mind that if we put the 1/2 factor here there we need to remove it
!ZN from the expression in which we use the strain later in the code.

  ! compute deviatoric strain
  if (NSPEC_STRAIN_ONLY == 1) then

    DO_LOOP_IJK
      templ = ONE_THIRD * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
      epsilondev_loc(INDEX_IJK,1) = duxdxl(INDEX_IJK) - templ
      epsilondev_loc(INDEX_IJK,2) = duydyl(INDEX_IJK) - templ
      epsilondev_loc(INDEX_IJK,3) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl(INDEX_IJK)
      epsilondev_loc(INDEX_IJK,4) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl(INDEX_IJK)
      epsilondev_loc(INDEX_IJK,5) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl(INDEX_IJK)
    ENDDO_LOOP_IJK

    if (ispec == 1) then
      DO_LOOP_IJK
        templ = ONE_THIRD * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
        epsilon_trace_over_3(INDEX_IJK,1) = templ
      ENDDO_LOOP_IJK
    endif

  else

    DO_LOOP_IJK
      templ = ONE_THIRD * (duxdxl(INDEX_IJK) + duydyl(INDEX_IJK) + duzdzl(INDEX_IJK))
      epsilon_trace_over_3(INDEX_IJK,ispec) = templ
      epsilondev_loc(INDEX_IJK,1) = duxdxl(INDEX_IJK) - templ
      epsilondev_loc(INDEX_IJK,2) = duydyl(INDEX_IJK) - templ
      epsilondev_loc(INDEX_IJK,3) = 0.5_CUSTOM_REAL * duxdyl_plus_duydxl(INDEX_IJK)
      epsilondev_loc(INDEX_IJK,4) = 0.5_CUSTOM_REAL * duzdxl_plus_duxdzl(INDEX_IJK)
      epsilondev_loc(INDEX_IJK,5) = 0.5_CUSTOM_REAL * duzdyl_plus_duydzl(INDEX_IJK)
    ENDDO_LOOP_IJK

  endif

  end subroutine compute_element_deviatoric_strain

!
!--------------------------------------------------------------------------------------------
!

! please leave this routine in this file, to help compilers inlining this function...

  subroutine compute_element_dot_product_stress(deriv_loc,jacobianl, &
                                                sigma_xx,sigma_yy,sigma_zz, &
                                                sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                                tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ATTRIBUTES INLINE :: compute_element_dot_product_stress
#else
! cray
!DIR$ INLINEALWAYS compute_element_dot_product_stress
#endif

! computes dot product of stress tensor with test vector, non-symmetric form

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  real(kind=CUSTOM_REAL), dimension(9,NGLLX,NGLLY,NGLLZ),intent(in) :: deriv_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: jacobianl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: tempy1,tempy2,tempy3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(out) :: tempz1,tempz2,tempz3

  ! local parameters
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: sxx,syy,szz,sxy,sxz,syz,syx,szx,szy

  real(kind=CUSTOM_REAL) :: fac

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! dot product of stress tensor with test vector, non-symmetric form
  DO_LOOP_IJK
    ! reloads derivatives of ux, uy and uz with respect to x, y and z
    xixl = deriv_loc(1,INDEX_IJK)
    xiyl = deriv_loc(2,INDEX_IJK)
    xizl = deriv_loc(3,INDEX_IJK)
    etaxl = deriv_loc(4,INDEX_IJK)
    etayl = deriv_loc(5,INDEX_IJK)
    etazl = deriv_loc(6,INDEX_IJK)
    gammaxl = deriv_loc(7,INDEX_IJK)
    gammayl = deriv_loc(8,INDEX_IJK)
    gammazl = deriv_loc(9,INDEX_IJK)

    ! common factor
    fac = jacobianl(INDEX_IJK)

    sxx = sigma_xx(INDEX_IJK)
    syy = sigma_yy(INDEX_IJK)
    szz = sigma_zz(INDEX_IJK)
    sxy = sigma_xy(INDEX_IJK)
    sxz = sigma_xz(INDEX_IJK)
    syz = sigma_yz(INDEX_IJK)
    syx = sigma_yx(INDEX_IJK)
    szx = sigma_zx(INDEX_IJK)
    szy = sigma_zy(INDEX_IJK)

    ! form dot product with test vector, non-symmetric form

    ! this goes to accel_x
    tempx1(INDEX_IJK) = fac * (sxx*xixl + syx*xiyl + szx*xizl)
    ! this goes to accel_y
    tempy1(INDEX_IJK) = fac * (sxy*xixl + syy*xiyl + szy*xizl)
    ! this goes to accel_z
    tempz1(INDEX_IJK) = fac * (sxz*xixl + syz*xiyl + szz*xizl)

    ! this goes to accel_x
    tempx2(INDEX_IJK) = fac * (sxx*etaxl + syx*etayl + szx*etazl)
    ! this goes to accel_y
    tempy2(INDEX_IJK) = fac * (sxy*etaxl + syy*etayl + szy*etazl)
    ! this goes to accel_z
    tempz2(INDEX_IJK) = fac * (sxz*etaxl + syz*etayl + szz*etazl)

    ! this goes to accel_x
    tempx3(INDEX_IJK) = fac * (sxx*gammaxl + syx*gammayl + szx*gammazl)
    ! this goes to accel_y
    tempy3(INDEX_IJK) = fac * (sxy*gammaxl + syy*gammayl + szy*gammazl)
    ! this goes to accel_z
    tempz3(INDEX_IJK) = fac * (sxz*gammaxl + syz*gammayl + szz*gammazl)
  ENDDO_LOOP_IJK

  end subroutine compute_element_dot_product_stress

!
!--------------------------------------------------------------------------------------------
!

! please leave this routine in this file, to help compilers inlining this function...

  subroutine compute_element_gravity(ispec,NSPEC,NGLOB,ibool,jacobianl,wgll_cube, &
                                     gravity_pre_store,gravity_H, &
                                     dummyx_loc,dummyy_loc,dummyz_loc, &
                                     sigma_xx,sigma_yy,sigma_zz, &
                                     sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
                                     rho_s_H)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES INLINE :: compute_element_gravity
#else
! cray
!DIR$ INLINEALWAYS compute_element_gravity
#endif

! computes non-symmetric stress terms for gravity

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

#ifdef FORCE_VECTORIZATION
  use constants, only: NGLLCUBE
#endif

  implicit none

  integer,intent(in) :: ispec
  integer,intent(in) :: NSPEC,NGLOB

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
!  real(kind=CUSTOM_REAL), dimension(3,NGLOB),intent(in) :: rstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: jacobianl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube

  ! gravity
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB),intent(in) :: gravity_pre_store
  real(kind=CUSTOM_REAL),dimension(6,NGLOB),intent(in) :: gravity_H

!  double precision, dimension(NRAD_GRAVITY),intent(in) :: minus_gravity_table,density_table,minus_deriv_gravity_table

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xx,sigma_yy,sigma_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(inout) :: rho_s_H

  ! local parameters
  ! for gravity
!  double precision :: dphi,dtheta
!  double precision :: radius,rho,minus_g,minus_dg
!  double precision :: minus_g_over_radius,minus_dg_plus_g_over_radius
!  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
!  double precision :: cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  real(kind=CUSTOM_REAL) :: factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  real(kind=CUSTOM_REAL) :: Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl

!  integer :: int_radius
  integer :: iglob

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! minimum radius in inner core (to avoid zero radius)
  !double precision, parameter :: MINIMUM_RADIUS_INNER_CORE = 100.d0 / R_PLANET

  ! computes non-symmetric terms for gravity
  DO_LOOP_IJK
    ! use mesh coordinates to get theta and phi
    ! x y and z contain r theta and phi
    iglob = ibool(INDEX_IJK,ispec)

    ! Cartesian components of the gravitational acceleration
    gxl = gravity_pre_store(1,iglob) ! minus_g*sin_theta*cos_phi * rho
    gyl = gravity_pre_store(2,iglob) ! minus_g*sin_theta*sin_phi * rho
    gzl = gravity_pre_store(3,iglob) ! minus_g*cos_theta * rho

    ! Cartesian components of gradient of gravitational acceleration
    ! get displacement and multiply by density to compute G tensor
    sx_l = dummyx_loc(INDEX_IJK)
    sy_l = dummyy_loc(INDEX_IJK)
    sz_l = dummyz_loc(INDEX_IJK)

    ! compute G tensor from s . g and add to sigma (not symmetric)
    sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) + sy_l * gyl + sz_l * gzl
    sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) + sx_l * gxl + sz_l * gzl
    sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + sx_l * gxl + sy_l * gyl

    sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - sx_l * gyl
    sigma_yx(INDEX_IJK) = sigma_yx(INDEX_IJK) - sy_l * gxl

    sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - sx_l * gzl
    sigma_zx(INDEX_IJK) = sigma_zx(INDEX_IJK) - sz_l * gxl

    sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - sy_l * gzl
    sigma_zy(INDEX_IJK) = sigma_zy(INDEX_IJK) - sz_l * gyl

    Hxxl = gravity_H(1,iglob) ! minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq * rho
    Hyyl = gravity_H(2,iglob) ! minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq * rho
    Hzzl = gravity_H(3,iglob) ! cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq * rho
    Hxyl = gravity_H(4,iglob) ! cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq * rho
    Hxzl = gravity_H(5,iglob) ! cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta * rho
    Hyzl = gravity_H(6,iglob) ! cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta * rho

    ! precompute vector
    factor = jacobianl(INDEX_IJK) * wgll_cube(INDEX_IJK)

    rho_s_H(1,INDEX_IJK) = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
    rho_s_H(2,INDEX_IJK) = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
    rho_s_H(3,INDEX_IJK) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)
  ENDDO_LOOP_IJK

  end subroutine compute_element_gravity

! left for reference: original routine...
!
! subroutine compute_element_gravity(ispec,NSPEC,NGLOB,ibool,rstore,jacobianl,wgll_cube, &
!                                     minus_gravity_table,minus_deriv_gravity_table,density_table, &
!                                     dummyx_loc,dummyy_loc,dummyz_loc, &
!                                     sigma_xx,sigma_yy,sigma_zz, &
!                                     sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy, &
!                                     rho_s_H)
!
!! computes non-symmetric stress terms for gravity
!
!  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,NRAD_GRAVITY,R_PLANET,R_PLANET_KM
!
!#ifdef FORCE_VECTORIZATION
!  use constants, only: NGLLCUBE
!#endif
!
!  implicit none
!
!  integer,intent(in) :: ispec
!  integer,intent(in) :: NSPEC,NGLOB
!
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
!  real(kind=CUSTOM_REAL), dimension(3,NGLOB),intent(in) :: rstore
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: jacobianl
!
!  double precision, dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: wgll_cube
!
!  ! gravity
!  double precision, dimension(NRAD_GRAVITY),intent(in) :: minus_gravity_table,density_table,minus_deriv_gravity_table
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: dummyx_loc,dummyy_loc,dummyz_loc
!
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xx,sigma_yy,sigma_zz
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ),intent(inout) :: sigma_xy,sigma_xz,sigma_yz,sigma_yx,sigma_zx,sigma_zy
!
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ),intent(out) :: rho_s_H
!
!  ! local parameters
!  ! for gravity
!  double precision :: dphi,dtheta
!  double precision :: radius,rho,minus_g,minus_dg
!  double precision :: minus_g_over_radius,minus_dg_plus_g_over_radius
!  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
!  double precision :: cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
!  double precision :: factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
!  double precision :: Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
!
!  integer :: int_radius
!  integer :: iglob
!
!#ifdef FORCE_VECTORIZATION
!  integer :: ijk
!#else
!  integer :: i,j,k
!#endif
!
!  ! minimum radius in inner core (to avoid zero radius)
!  double precision, parameter :: MINIMUM_RADIUS_INNER_CORE = 100.d0 / R_PLANET
!
!  ! computes non-symmetric terms for gravity
!  DO_LOOP_IJK
!
!    ! use mesh coordinates to get theta and phi
!    ! x y and z contain r theta and phi
!    iglob = ibool(INDEX_IJK,ispec)
!
!    radius = dble(rstore(1,iglob))
!    dtheta = dble(rstore(2,iglob))
!    dphi = dble(rstore(3,iglob))
!
!    ! make sure radius is never zero even for points at center of cube
!    ! because we later divide by radius
!    if (radius < MINIMUM_RADIUS_INNER_CORE) radius = MINIMUM_RADIUS_INNER_CORE
!
!    cos_theta = dcos(dtheta)
!    sin_theta = dsin(dtheta)
!    cos_phi = dcos(dphi)
!    sin_phi = dsin(dphi)
!
!    cos_theta_sq = cos_theta*cos_theta
!    sin_theta_sq = sin_theta*sin_theta
!    cos_phi_sq = cos_phi*cos_phi
!    sin_phi_sq = sin_phi*sin_phi
!
!    ! get g, rho and dg/dr=dg
!    ! spherical components of the gravitational acceleration
!
!    ! for efficiency replace with lookup table every 100 m in radial direction
!    !int_radius = nint(10.d0 * radius * R_PLANET_KM )
!
!    ! make sure we never use zero for point exactly at the center of the Earth
!    int_radius = max(1,nint(10.d0 * radius * R_PLANET_KM))
!
!    minus_g = minus_gravity_table(int_radius)
!    minus_dg = minus_deriv_gravity_table(int_radius)
!    rho = density_table(int_radius)
!
!    ! Cartesian components of the gravitational acceleration
!    gxl = minus_g*sin_theta*cos_phi
!    gyl = minus_g*sin_theta*sin_phi
!    gzl = minus_g*cos_theta
!
!    ! Cartesian components of gradient of gravitational acceleration
!    ! obtained from spherical components
!    minus_g_over_radius = minus_g / radius
!    minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius
!
!    Hxxl = minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
!    Hyyl = minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
!    Hzzl = cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
!    Hxyl = cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
!    Hxzl = cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
!    Hyzl = cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta
!
!    ! get displacement and multiply by density to compute G tensor
!    sx_l = rho * dble(dummyx_loc(INDEX_IJK)) ! dble(displ_crust_mantle(1,iglob))
!    sy_l = rho * dble(dummyy_loc(INDEX_IJK)) ! dble(displ_crust_mantle(2,iglob))
!    sz_l = rho * dble(dummyz_loc(INDEX_IJK)) ! dble(displ_crust_mantle(3,iglob))
!
!    ! compute G tensor from s . g and add to sigma (not symmetric)
!    sigma_xx(INDEX_IJK) = sigma_xx(INDEX_IJK) + real(sy_l*gyl + sz_l*gzl, kind=CUSTOM_REAL)
!    sigma_yy(INDEX_IJK) = sigma_yy(INDEX_IJK) + real(sx_l*gxl + sz_l*gzl, kind=CUSTOM_REAL)
!    sigma_zz(INDEX_IJK) = sigma_zz(INDEX_IJK) + real(sx_l*gxl + sy_l*gyl, kind=CUSTOM_REAL)
!
!    sigma_xy(INDEX_IJK) = sigma_xy(INDEX_IJK) - real(sx_l * gyl, kind=CUSTOM_REAL)
!    sigma_yx(INDEX_IJK) = sigma_yx(INDEX_IJK) - real(sy_l * gxl, kind=CUSTOM_REAL)
!
!    sigma_xz(INDEX_IJK) = sigma_xz(INDEX_IJK) - real(sx_l * gzl, kind=CUSTOM_REAL)
!    sigma_zx(INDEX_IJK) = sigma_zx(INDEX_IJK) - real(sz_l * gxl, kind=CUSTOM_REAL)
!
!    sigma_yz(INDEX_IJK) = sigma_yz(INDEX_IJK) - real(sy_l * gzl, kind=CUSTOM_REAL)
!    sigma_zy(INDEX_IJK) = sigma_zy(INDEX_IJK) - real(sz_l * gyl, kind=CUSTOM_REAL)
!
!    ! precompute vector
!    factor = dble(jacobianl(INDEX_IJK)) * wgll_cube(INDEX_IJK)
!
!    rho_s_H(1,INDEX_IJK) = real(factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl), kind=CUSTOM_REAL)
!    rho_s_H(2,INDEX_IJK) = real(factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl), kind=CUSTOM_REAL)
!    rho_s_H(3,INDEX_IJK) = real(factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl), kind=CUSTOM_REAL)
!
!  ENDDO_LOOP_IJK
!
!  end subroutine compute_element_gravity

