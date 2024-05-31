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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



  subroutine compute_forces_outer_core_Dev(timeval,deltat,two_omega_earth, &
                                           NSPEC_ROT,NGLOB, &
                                           A_array_rotation,B_array_rotation, &
                                           A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                           displfluid,accelfluid, &
                                           div_displfluid,iphase,sum_terms, &
                                           pgrav_outer_core)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver, only: &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLCUBE,NDIM, &
    NSPEC_OUTER_CORE_3DMOVIE,NSPEC_OUTER_CORE, &
    ROTATION_VAL,GRAVITY_VAL, &
    NSTAGE,ALPHA_LDDRK,BETA_LDDRK, &
    FULL_GRAVITY_VAL,DISCARD_GCONTRIB, &
    m1,m2

  use specfem_par, only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgll_cube, &
    gravity_pre_store => gravity_pre_store_outer_core, &
    MOVIE_VOLUME, &
    USE_LDDRK,istage

  use specfem_par_outercore, only: &
    deriv => deriv_mapping_outer_core, &
    ibool => ibool_outer_core, &
    phase_ispec_inner => phase_ispec_inner_outer_core, &
    nspec_outer => nspec_outer_outer_core, &
    nspec_inner => nspec_inner_outer_core

  use specfem_par, only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

#ifdef FORCE_VECTORIZATION
  use specfem_par_outercore, only: &
    ibool_inv_tbl => ibool_inv_tbl_outer_core, &
    ibool_inv_st => ibool_inv_st_outer_core, &
    num_globs => num_globs_outer_core, &
    phase_iglob => phase_iglob_outer_core
#endif

  ! full gravity
  use specfem_par_full_gravity, only: &
    gravity_rho_g_over_kappa => gravity_rho_g_over_kappa_outer_core

  implicit none

  integer,intent(in) :: NSPEC_ROT,NGLOB

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL),intent(in) :: timeval,deltat,two_omega_earth

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ROT),intent(inout) :: &
    A_array_rotation,B_array_rotation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_ROT),intent(inout) :: &
    A_array_rotation_lddrk,B_array_rotation_lddrk

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: displfluid
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(inout) :: accelfluid

  ! divergence of displacement
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE),intent(out) :: div_displfluid

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),intent(out) :: sum_terms

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: pgrav_outer_core

  ! local parameters
  ! for gravity
  real(kind=CUSTOM_REAL) :: gravity_term,vec_x,vec_y,vec_z

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t
  real(kind=CUSTOM_REAL) :: ux_rotation,uy_rotation
  real(kind=CUSTOM_REAL) :: source_euler_A,source_euler_B
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: A_rotation,B_rotation

  integer :: ispec,iglob
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobianl

  ! Deville
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: chi_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: temp1,temp2,temp3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtemp1,newtemp2,newtemp3

  integer :: num_elements,ispec_p
  integer :: i,j,k
#ifdef FORCE_VECTORIZATION
  integer :: ijk
  integer :: ijk_spec, ip, iglob_p
#endif
  real(kind=CUSTOM_REAL), dimension(NSTAGE) :: MYALPHA_LDDRK,MYBETA_LDDRK

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

  if (MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE > 1 .and. (iphase == 1)) then
    div_displfluid(:,:,:,:) = 0._CUSTOM_REAL
  endif

  if (iphase == 1) then
    ! outer elements
    num_elements = nspec_outer
  else
    num_elements = nspec_inner
  endif

  ! note: this is an OpenMP work-around for gfortran
  ! gfortran has problems with parameters inside a DEFAULT(NONE) OpenMP block
  ! for a fix see https://gcc.gnu.org/ml/Fortran/2014-10/msg00064.html
  !
  ! we use local variables to cirumvent it, another possibility would be to use DEFAULT(SHARED)
  !
  ! FIRSTPRIVATE declaration is chosen based on suggestion for performance optimization:
  ! see https://docs.oracle.com/cd/E19059-01/stud.10/819-0501/7_tuning.html

  if (USE_LDDRK) then
    MYALPHA_LDDRK(:) = ALPHA_LDDRK(:)
    MYBETA_LDDRK(:) = BETA_LDDRK(:)
  endif

! openmp solver
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED( deriv, &
!$OMP num_elements, phase_ispec_inner, iphase, ibool, &
!$OMP displfluid, accelfluid, &
!$OMP gravity_pre_store, gravity_rho_g_over_kappa, pgrav_outer_core, &
!$OMP deltat, two_omega_earth, timeval, &
!$OMP A_array_rotation, B_array_rotation, &
!$OMP A_array_rotation_lddrk, B_array_rotation_lddrk, &
!$OMP sum_terms, &
#ifdef FORCE_VECTORIZATION
!$OMP ibool_inv_tbl, ibool_inv_st, num_globs, phase_iglob, &
#endif
!$OMP istage, div_displfluid ) &
!$OMP PRIVATE( &
!$OMP ispec_p, ispec, i, j, k, iglob, chi_elem, &
!$OMP xixl, xiyl, xizl, etaxl, etayl, etazl, gammaxl, gammayl, gammazl, jacobianl, &
!$OMP dpotentialdxl, dpotentialdyl, dpotentialdzl, temp1, temp2, temp3, &
!$OMP two_omega_deltat, cos_two_omega_t, sin_two_omega_t, &
!$OMP source_euler_A, source_euler_B, A_rotation, B_rotation, ux_rotation, uy_rotation, &
!$OMP gravity_term,vec_x,vec_y,vec_z, &
#ifdef FORCE_VECTORIZATION
!$OMP ijk, &
!$OMP ijk_spec, ip, iglob_p, &
#endif
!$OMP newtemp1, newtemp2, newtemp3 ) &
!$OMP FIRSTPRIVATE( hprime_xx, hprime_xxT, hprimewgll_xxT, hprimewgll_xx, &
!$OMP wgllwgll_yz_3D, wgllwgll_xz_3D, wgllwgll_xy_3D, wgll_cube, &
!$OMP MYALPHA_LDDRK,MYBETA_LDDRK, &
!$OMP USE_LDDRK,ROTATION_VAL,GRAVITY_VAL,FULL_GRAVITY_VAL,MOVIE_VOLUME, &
!$OMP NSPEC_OUTER_CORE_3DMOVIE )

!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner(ispec_p,iphase)

    ! only compute element which belong to current phase (inner or outer elements)

    ! note: this loop will not fully vectorize because it contains a dependency (through indirect addressing with array ibool())
    !       thus, instead of DO_LOOP_IJK we use do k=..;do j=..;do i=..,
    !       which helps the compiler to unroll the innermost loop
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
          ! gets "displacement"
          chi_elem(i,j,k) = displfluid(iglob)
        enddo
      enddo
    enddo

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for temp1
    ! computes 2. matrix multiplication for temp2
    ! computes 3. matrix multiplication for temp3
    call mxm5_single(hprime_xx,m1,chi_elem,temp1,m2)
    call mxm5_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    call mxm5_single(chi_elem,m2,hprime_xxT,temp3,m1)

    ! note: this compute_forces_outer_core_Dev() routine is called for USE_DEVILLE_PRODUCTS_VAL == .true.
    !       which is only the case for NGLLX == NGLLY == NGLLZ == 5
    !
    ! for more general cases one could do the following:
    !select case (NGLLX)
    !case (5)
    !  call mxm5_single(hprime_xx,m1,chi_elem,temp1,m2)
    !  call mxm5_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    !  call mxm5_single(chi_elem,m2,hprime_xxT,temp3,m1)
    !case (6)
    !  call mxm6_single(hprime_xx,m1,chi_elem,temp1,m2)
    !  call mxm6_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    !  call mxm6_single(chi_elem,m2,hprime_xxT,temp3,m1)
    !case (7)
    !  call mxm7_single(hprime_xx,m1,chi_elem,temp1,m2)
    !  call mxm7_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    !  call mxm7_single(chi_elem,m2,hprime_xxT,temp3,m1)
    !case (8)
    !  call mxm8_single(hprime_xx,m1,chi_elem,temp1,m2)
    !  call mxm8_3dmat_single(chi_elem,m1,hprime_xxT,m1,temp2,NGLLX)
    !  call mxm8_single(chi_elem,m2,hprime_xxT,temp3,m1)
    !end select

    DO_LOOP_IJK
      ! get derivatives of potential with respect to x, y and z
      xixl = deriv(1,INDEX_IJK,ispec)
      xiyl = deriv(2,INDEX_IJK,ispec)
      xizl = deriv(3,INDEX_IJK,ispec)
      etaxl = deriv(4,INDEX_IJK,ispec)
      etayl = deriv(5,INDEX_IJK,ispec)
      etazl = deriv(6,INDEX_IJK,ispec)
      gammaxl = deriv(7,INDEX_IJK,ispec)
      gammayl = deriv(8,INDEX_IJK,ispec)
      gammazl = deriv(9,INDEX_IJK,ispec)

      ! compute the Jacobian
      jacobianl(INDEX_IJK) = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                                             - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                                             + xizl*(etaxl*gammayl-etayl*gammaxl))

      dpotentialdxl(INDEX_IJK) = xixl*temp1(INDEX_IJK) + etaxl*temp2(INDEX_IJK) + gammaxl*temp3(INDEX_IJK)
      dpotentialdyl(INDEX_IJK) = xiyl*temp1(INDEX_IJK) + etayl*temp2(INDEX_IJK) + gammayl*temp3(INDEX_IJK)
      dpotentialdzl(INDEX_IJK) = xizl*temp1(INDEX_IJK) + etazl*temp2(INDEX_IJK) + gammazl*temp3(INDEX_IJK)
    ENDDO_LOOP_IJK

    ! compute contribution of rotation and add to gradient of potential
    ! this term has no Z component
    if (ROTATION_VAL) then
      ! store the source for the Euler scheme for A_rotation and B_rotation
      two_omega_deltat = deltat * two_omega_earth

      cos_two_omega_t = cos(two_omega_earth*timeval)
      sin_two_omega_t = sin(two_omega_earth*timeval)

      ! update rotation term with Euler scheme
      !
      ! note: rotation with euler scheme is not exactly revertible;
      !       backpropagation/reconstruction of wavefield will lead to small differences
      if (USE_LDDRK) then
        ! use the source saved above
        DO_LOOP_IJK
          A_rotation(INDEX_IJK) = A_array_rotation(INDEX_IJK,ispec)
          B_rotation(INDEX_IJK) = B_array_rotation(INDEX_IJK,ispec)

          ! time step deltat of Euler scheme is included in the source
          source_euler_A = two_omega_deltat &
                * (cos_two_omega_t * dpotentialdyl(INDEX_IJK) + sin_two_omega_t * dpotentialdxl(INDEX_IJK))
          source_euler_B = two_omega_deltat &
                * (sin_two_omega_t * dpotentialdyl(INDEX_IJK) - cos_two_omega_t * dpotentialdxl(INDEX_IJK))

          A_array_rotation_lddrk(INDEX_IJK,ispec) = MYALPHA_LDDRK(istage) * A_array_rotation_lddrk(INDEX_IJK,ispec) &
                                                    + source_euler_A
          A_array_rotation(INDEX_IJK,ispec) = A_array_rotation(INDEX_IJK,ispec) &
                                              + MYBETA_LDDRK(istage) * A_array_rotation_lddrk(INDEX_IJK,ispec)

          B_array_rotation_lddrk(INDEX_IJK,ispec) = MYALPHA_LDDRK(istage) * B_array_rotation_lddrk(INDEX_IJK,ispec) &
                                                    + source_euler_B
          B_array_rotation(INDEX_IJK,ispec) = B_array_rotation(INDEX_IJK,ispec) &
                                              + MYBETA_LDDRK(istage) * B_array_rotation_lddrk(INDEX_IJK,ispec)
        ENDDO_LOOP_IJK
      else
        ! Euler scheme
        DO_LOOP_IJK
          A_rotation(INDEX_IJK) = A_array_rotation(INDEX_IJK,ispec)
          B_rotation(INDEX_IJK) = B_array_rotation(INDEX_IJK,ispec)

          ! time step deltat of Euler scheme is included in the source
          source_euler_A = two_omega_deltat &
                * (cos_two_omega_t * dpotentialdyl(INDEX_IJK) + sin_two_omega_t * dpotentialdxl(INDEX_IJK))
          source_euler_B = two_omega_deltat &
                * (sin_two_omega_t * dpotentialdyl(INDEX_IJK) - cos_two_omega_t * dpotentialdxl(INDEX_IJK))

          A_array_rotation(INDEX_IJK,ispec) = A_array_rotation(INDEX_IJK,ispec) + source_euler_A
          B_array_rotation(INDEX_IJK,ispec) = B_array_rotation(INDEX_IJK,ispec) + source_euler_B
        ENDDO_LOOP_IJK
      endif

      ! adds rotation
      DO_LOOP_IJK
        ux_rotation =   A_rotation(INDEX_IJK)*cos_two_omega_t + B_rotation(INDEX_IJK)*sin_two_omega_t
        uy_rotation = - A_rotation(INDEX_IJK)*sin_two_omega_t + B_rotation(INDEX_IJK)*cos_two_omega_t

        dpotentialdxl(INDEX_IJK) = dpotentialdxl(INDEX_IJK) + ux_rotation
        dpotentialdyl(INDEX_IJK) = dpotentialdyl(INDEX_IJK) + uy_rotation
      ENDDO_LOOP_IJK
    endif  ! end of section with rotation

    ! add (chi/rho)grad(rho) term in no gravity case
    if (.not. GRAVITY_VAL) then
      ! With regards to the non-gravitating case: we cannot set N^2 = 0 *and* let g = 0.
      ! We can *either* assume N^2 = 0 but keep gravity g, *or* we can assume that gravity
      ! is negligible to begin with, as in our GJI 2002a, in which case N does not arise.
      ! We get:
      !
      ! \ddot\chi = \rho^{-1}\kappa\bdel\cdot(\bdel\chi+\chi\bdel\ln\rho)
      !
      ! Then the displacement is
      !
      ! \bu = \bdel\chi+\chi\bdel\ln\rho = \rho^{-1}\bdel(\rho\chi)
      !
      ! and the pressure is
      !
      ! p = -\rho\ddot{\chi}
      !
      ! Thus in our 2002b GJI paper eqn (21) is wrong, and equation (41)
      ! in our AGU monograph is incorrect; these equations should be replaced by
      !
      ! \ddot\chi = \rho^{-1}\kappa\bdel\cdot(\bdel\chi+\chi\bdel\ln\rho)
      !
      ! Note that the fluid potential we use in GJI 2002a differs from the one used here:
      !
      ! \chi_GJI2002a = \rho\partial\t\chi
      !
      ! such that
      !
      ! \bv = \partial_t\bu=\rho^{-1}\bdel\chi_GJI2002a  (GJI 2002a eqn 20)
      !
      ! p = - \partial_t\chi_GJI2002a (GJI 2002a eqn 19)
      DO_LOOP_IJK
        iglob = ibool(INDEX_IJK,ispec)

        ! gradient of d ln(rho)/dr in Cartesian coordinates
        vec_x = gravity_pre_store(1,iglob)
        vec_y = gravity_pre_store(2,iglob)
        vec_z = gravity_pre_store(3,iglob)

        ! grad(rho)/rho in Cartesian components
        dpotentialdxl(INDEX_IJK) = dpotentialdxl(INDEX_IJK) + chi_elem(INDEX_IJK) * vec_x
        dpotentialdyl(INDEX_IJK) = dpotentialdyl(INDEX_IJK) + chi_elem(INDEX_IJK) * vec_y
        dpotentialdzl(INDEX_IJK) = dpotentialdzl(INDEX_IJK) + chi_elem(INDEX_IJK) * vec_z
      ENDDO_LOOP_IJK
    endif

    DO_LOOP_IJK
      ! reloads derivatives of ux, uy and uz with respect to x, y and z
      xixl = deriv(1,INDEX_IJK,ispec)
      xiyl = deriv(2,INDEX_IJK,ispec)
      xizl = deriv(3,INDEX_IJK,ispec)
      etaxl = deriv(4,INDEX_IJK,ispec)
      etayl = deriv(5,INDEX_IJK,ispec)
      etazl = deriv(6,INDEX_IJK,ispec)
      gammaxl = deriv(7,INDEX_IJK,ispec)
      gammayl = deriv(8,INDEX_IJK,ispec)
      gammazl = deriv(9,INDEX_IJK,ispec)

      temp1(INDEX_IJK) = jacobianl(INDEX_IJK) * &
                         (xixl*dpotentialdxl(INDEX_IJK) + xiyl*dpotentialdyl(INDEX_IJK) + xizl*dpotentialdzl(INDEX_IJK))
      temp2(INDEX_IJK) = jacobianl(INDEX_IJK) * &
                         (etaxl*dpotentialdxl(INDEX_IJK) + etayl*dpotentialdyl(INDEX_IJK) + etazl*dpotentialdzl(INDEX_IJK))
      temp3(INDEX_IJK) = jacobianl(INDEX_IJK) * &
                         (gammaxl*dpotentialdxl(INDEX_IJK) + gammayl*dpotentialdyl(INDEX_IJK) + gammazl*dpotentialdzl(INDEX_IJK))
    ENDDO_LOOP_IJK

    ! second double-loop over GLL to compute all the terms along the x,y,z directions and assemble the contributions

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtemp1
    ! computes 2. matrix multiplication for newtemp2
    ! computes 3. matrix multiplication for newtemp3
    call mxm5_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    call mxm5_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    call mxm5_single(temp3,m2,hprimewgll_xx,newtemp3,m1)

    ! note: this compute_forces_outer_core_Dev() routine is called for USE_DEVILLE_PRODUCTS_VAL == .true.
    !       which is only the case for NGLLX == NGLLY == NGLLZ == 5
    !
    ! for more general cases one could do the following:
    !select case (NGLLX)
    !case (5)
    !  call mxm5_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    !  call mxm5_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    !  call mxm5_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    !case (6)
    !  call mxm6_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    !  call mxm6_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    !  call mxm6_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    !case (7)
    !  call mxm7_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    !  call mxm7_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    !  call mxm7_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    !case (8)
    !  call mxm8_single(hprimewgll_xxT,m1,temp1,newtemp1,m2)
    !  call mxm8_3dmat_single(temp2,m1,hprimewgll_xx,m1,newtemp2,NGLLX)
    !  call mxm8_single(temp3,m2,hprimewgll_xx,newtemp3,m1)
    !end select

    ! sum contributions from each element to the global mesh and add gravity term
    DO_LOOP_IJK
      sum_terms(INDEX_IJK,ispec) = - ( wgllwgll_yz_3D(INDEX_IJK)*newtemp1(INDEX_IJK) &
                                     + wgllwgll_xz_3D(INDEX_IJK)*newtemp2(INDEX_IJK) &
                                     + wgllwgll_xy_3D(INDEX_IJK)*newtemp3(INDEX_IJK))
    ENDDO_LOOP_IJK

    ! adds gravity term
    if (GRAVITY_VAL) then
      DO_LOOP_IJK
        ! uses double precision calculations. needs testing if okay with single precision.
        iglob = ibool(INDEX_IJK,ispec)

        ! Cartesian components of the gravitational acceleration
        ! gravitational acceleration (integrated and multiply by rho / Kappa)
        vec_x = gravity_pre_store(1,iglob)
        vec_y = gravity_pre_store(2,iglob)
        vec_z = gravity_pre_store(3,iglob)

        ! compute divergence of displacement
        ! distinguish between single and double precision for reals
        gravity_term = jacobianl(INDEX_IJK) * wgll_cube(INDEX_IJK) &
                            * (dpotentialdxl(INDEX_IJK) * vec_x &
                             + dpotentialdyl(INDEX_IJK) * vec_y &
                             + dpotentialdzl(INDEX_IJK) * vec_z)

        ! full gravity contribution
        if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
          gravity_term = gravity_term &
              - gravity_rho_g_over_kappa(iglob) * jacobianl(INDEX_IJK) * wgll_cube(INDEX_IJK) * pgrav_outer_core(iglob)
        endif

        ! divergence of displacement field with gravity on
        ! note: these calculations are only considered for SIMULATION_TYPE == 1 .and. SAVE_FORWARD
        !          and one has set MOVIE_VOLUME_TYPE == 4 when MOVIE_VOLUME is .true.;
        !         in case of SIMULATION_TYPE == 3, it gets overwritten by compute_kernels_outer_core()
        if (MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE > 1) then
          div_displfluid(INDEX_IJK,ispec) = dpotentialdxl(INDEX_IJK) * vec_x &
                                          + dpotentialdyl(INDEX_IJK) * vec_y &
                                          + dpotentialdzl(INDEX_IJK) * vec_z
        endif

        ! adds gravity term
        sum_terms(INDEX_IJK,ispec) = sum_terms(INDEX_IJK,ispec) + gravity_term
      ENDDO_LOOP_IJK
    endif

  enddo   ! spectral element loop
!$OMP ENDDO

  ! updates acceleration
#ifdef FORCE_VECTORIZATION
  ! updates acceleration
!$OMP DO
  do iglob_p = 1,num_globs(iphase)
    iglob = phase_iglob(iglob_p,iphase)
    do ip = ibool_inv_st(iglob_p,iphase),ibool_inv_st(iglob_p+1,iphase)-1
      ijk_spec = ibool_inv_tbl(ip,iphase)
      accelfluid(iglob) = accelfluid(iglob) + sum_terms(ijk_spec,1,1,1)
    enddo
  enddo
!$OMP ENDDO
#else
  ! updates for non-vectorization case
!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements
    ispec = phase_ispec_inner(ispec_p,iphase)
    DO_LOOP_IJK
      iglob = ibool(INDEX_IJK,ispec)
!$OMP ATOMIC
      accelfluid(iglob) = accelfluid(iglob) + sum_terms(INDEX_IJK,ispec)
    ENDDO_LOOP_IJK
  enddo
!$OMP ENDDO
#endif
!$OMP END PARALLEL

  contains

!--------------------------------------------------------------------------------------------
!
! matrix-matrix multiplications
!
! subroutines adapted from Deville, Fischer and Mund, High-order methods
! for incompressible fluid flow, Cambridge University Press (2002),
! pages 386 and 389 and Figure 8.3.1
!
!--------------------------------------------------------------------------------------------
!
! note: the matrix-matrix multiplications are used for very small matrices ( 5 x 5 x 5 elements);
!       thus, calling external optimized libraries for these multiplications are in general slower
!
! please leave the routines here to help compilers inlining the code

  subroutine mxm5_single(A,n1,B,C,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_single
#else
! cray
! note: with Cray Fortran versions >= 14 on Frontier, inlining this routine together with optimization -O3 leads to problems.
!       for now, will avoid inlining by this directive INLINENEVER to allow for default compilation,
!       otherwise the compilation flag -hipa0 would need to be added to suppress all inlining as well.
!!DIR$ INLINEALWAYS mxm5_single
!DIR$ INLINENEVER mxm5_single
#endif

! 2-dimensional arrays (25,5)/(5,25)

  use constants_solver, only: CUSTOM_REAL

#ifdef USE_XSMM
  !use my_libxsmm, only: libxsmm_smm_5_25_5,libxsmm_smm_25_5_5
  use my_libxsmm, only: xmm1,xmm2,libxsmm_mmcall_abc => libxsmm_smmcall_abc
#endif

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C

  ! local parameters
  integer :: i,j

! note: on CPUs like Haswell or Sandy Bridge, the following will slow down computations
!       however, on Intel Phi (KNC) it is still helpful
#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  ! matrix-matrix multiplication C = alpha A * B + beta C
  if (n1 == 5) then
    ! with A(n1,n2) 5x5-matrix, B(n2,n3) 5x25-matrix and C(n1,n3) 5x25-matrix
    !call libxsmm_smm_5_25_5(a=A, b=B, c=C, pa=A, pb=B, pc=C)
    ! dispatch
    call libxsmm_mmcall_abc(xmm1, A, B, C)
  else if (n1 == 25) then
    ! with A(n1,n2) 25x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3) 25x5-matrix
    !call libxsmm_smm_25_5_5(a=A, b=B, c=C, pa=A, pb=B, pc=C)
    ! dispatch
    call libxsmm_mmcall_abc(xmm2, A, B, C)
  else
    stop 'Invalid n1 value for LibXSMM in mxm5_single routine'
  endif
  return
#endif

  ! matrix-matrix multiplication
  do j = 1,n3
    do i = 1,n1
      C(i,j) =  A(i,1) * B(1,j) &
              + A(i,2) * B(2,j) &
              + A(i,3) * B(3,j) &
              + A(i,4) * B(4,j) &
              + A(i,5) * B(5,j)
    enddo
  enddo

  end subroutine mxm5_single

  !-------------

! unused so far..
!
!  subroutine mxm6_single(A,n1,B,C,n3)
!
!! two-dimensional arrays (36,6)/(6,36)
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,6),intent(in) :: A
!  real(kind=CUSTOM_REAL),dimension(6,n3),intent(in) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j) &
!              + A(i,6) * B(6,j)
!    enddo
!  enddo
!
!  end subroutine mxm6_single

  !-------------

! unused so far..
!
!  subroutine mxm7_single(A,n1,B,C,n3)
!
!! two-dimensional arrays (49,7)/(7,49)
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,7),intent(in) :: A
!  real(kind=CUSTOM_REAL),dimension(7,n3),intent(in) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j) &
!              + A(i,6) * B(6,j) &
!              + A(i,7) * B(7,j)
!    enddo
!  enddo
!
!  end subroutine mxm7_single

  !-------------

! unused so far..
!
!  subroutine mxm8_single(A,n1,B,C,n3)
!
!! two-dimensional arrays (64,8)/(8,64)
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n3
!  real(kind=CUSTOM_REAL),dimension(n1,8),intent(in) :: A
!  real(kind=CUSTOM_REAL),dimension(8,n3),intent(in) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n3),intent(out) :: C
!
!  ! local parameters
!  integer :: i,j
!
!  ! matrix-matrix multiplication
!  do j = 1,n3
!    do i = 1,n1
!      C(i,j) =  A(i,1) * B(1,j) &
!              + A(i,2) * B(2,j) &
!              + A(i,3) * B(3,j) &
!              + A(i,4) * B(4,j) &
!              + A(i,5) * B(5,j) &
!              + A(i,6) * B(6,j) &
!              + A(i,7) * B(7,j) &
!              + A(i,8) * B(8,j)
!    enddo
!  enddo
!
!  end subroutine mxm8_single

!--------------------------------------------------------------------------------------------

  subroutine mxm5_3dmat_single(A,n1,B,n2,C,n3)

! we can force inlining (Intel compiler)
#if defined __INTEL_COMPILER
!DIR$ ATTRIBUTES FORCEINLINE :: mxm5_3dmat_single
#else
! cray
! note: with Cray Fortran versions >= 14 on Frontier, inlining this routine together with optimization -O3 leads to problems.
!       for now, will avoid inlining by this directive INLINENEVER to allow for default compilation,
!       otherwise the compilation flag -hipa0 would need to be added to suppress all inlining as well.
!!DIR$ INLINEALWAYS mxm5_3dmat_single
!DIR$ INLINENEVER mxm5_3dmat_single
#endif

! 3-dimensional arrays (5,5,5) for A and C

  use constants_solver, only: CUSTOM_REAL

! note: on CPUs like Haswell or Sandy Bridge, the following will slow down computations
!       however, on Intel Phi (KNC) it is still helpful (speedup +3%)
#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  !use my_libxsmm, only: libxsmm_smm_5_5_5
  use my_libxsmm, only: xmm3,libxsmm_mmcall_abc => libxsmm_smmcall_abc
#endif

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3),intent(in) :: A
  real(kind=CUSTOM_REAL),dimension(5,n2),intent(in) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C

  ! local parameters
  integer :: i,j,k

#if defined(XSMM_FORCE_EVEN_IF_SLOWER) || ( defined(XSMM) && defined(__MIC__) )
  ! matrix-matrix multiplication C = alpha A * B + beta C
  ! with A(n1,n2,n4) 5x5x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3,n4) 5x5x5-matrix
  ! static version
  !do k = 1,5
  !  call libxsmm_call(xmm3, A(:,:,k), B, C(:,:,k))
  !enddo
  ! dispatch
  do k = 1,5
    call libxsmm_mmcall_abc(xmm3, A(1,1,k), B, C(1,1,k))
  enddo
  ! unrolled
  !call libxsmm_smm_5_5_5(a=A(1,1,1), b=B, c=C(1,1,1),pa=A(1,1,1+1), pb=B, pc=C(1,1,1+1))
  !call libxsmm_smm_5_5_5(a=A(1,1,2), b=B, c=C(1,1,2),pa=A(1,1,2+1), pb=B, pc=C(1,1,2+1))
  !call libxsmm_smm_5_5_5(a=A(1,1,3), b=B, c=C(1,1,3),pa=A(1,1,3+1), pb=B, pc=C(1,1,3+1))
  !call libxsmm_smm_5_5_5(a=A(1,1,4), b=B, c=C(1,1,4),pa=A(1,1,4+1), pb=B, pc=C(1,1,4+1))
  !call libxsmm_smm_5_5_5(a=A(1,1,5), b=B, c=C(1,1,5),pa=A(1,1,1), pb=B, pc=C(1,1,1))
  return
#endif

  ! matrix-matrix multiplication
  do j = 1,n2
    do i = 1,n1
      ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
      do k = 1,n3
        C(i,j,k) =  A(i,1,k) * B(1,j) &
                  + A(i,2,k) * B(2,j) &
                  + A(i,3,k) * B(3,j) &
                  + A(i,4,k) * B(4,j) &
                  + A(i,5,k) * B(5,j)
      enddo
    enddo
  enddo

  end subroutine mxm5_3dmat_single

  !-------------

! unused so far..
!
!  subroutine mxm6_3dmat_single(A,n1,B,n2,C,n3)
!
!! three-dimensional arrays (6,6,6) for A and C
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,6,n3),intent(in) :: A
!  real(kind=CUSTOM_REAL),dimension(6,n2),intent(in) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C
!
!  ! local parameters
!  integer :: i,j,k
!
!  ! matrix-matrix multiplication
!  do k = 1,n3
!    do j = 1,n2
!      do i = 1,n1
!        C(i,j,k) =  A(i,1,k) * B(1,j) &
!                  + A(i,2,k) * B(2,j) &
!                  + A(i,3,k) * B(3,j) &
!                  + A(i,4,k) * B(4,j) &
!                  + A(i,5,k) * B(5,j) &
!                  + A(i,6,k) * B(6,j)
!      enddo
!    enddo
!  enddo
!
!  end subroutine mxm6_3dmat_single

  !-------------

! unused so far..
!
!  subroutine mxm7_3dmat_single(A,n1,B,n2,C,n3)
!
!! three-dimensional arrays (7,7,7) for A and C
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,7,n3),intent(in) :: A
!  real(kind=CUSTOM_REAL),dimension(7,n2),intent(in) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C
!
!  ! local parameters
!  integer :: i,j,k
!
!  ! matrix-matrix multiplication
!  do k = 1,n3
!    do j = 1,n2
!      do i = 1,n1
!        C(i,j,k) =  A(i,1,k) * B(1,j) &
!                  + A(i,2,k) * B(2,j) &
!                  + A(i,3,k) * B(3,j) &
!                  + A(i,4,k) * B(4,j) &
!                  + A(i,5,k) * B(5,j) &
!                  + A(i,6,k) * B(6,j) &
!                  + A(i,7,k) * B(7,j)
!      enddo
!    enddo
!  enddo
!
!  end subroutine mxm7_3dmat_single

  !-------------

! unused so far..
!
!  subroutine mxm8_3dmat_single(A,n1,B,n2,C,n3)
!
!! three-dimensional arrays (8,8,8) for A and C
!
!  use constants, only: CUSTOM_REAL
!
!  implicit none
!
!  integer,intent(in) :: n1,n2,n3
!  real(kind=CUSTOM_REAL),dimension(n1,8,n3),intent(in) :: A
!  real(kind=CUSTOM_REAL),dimension(8,n2),intent(in) :: B
!  real(kind=CUSTOM_REAL),dimension(n1,n2,n3),intent(out) :: C
!
!  ! local parameters
!  integer :: i,j,k
!
!  ! matrix-matrix multiplication
!  do k = 1,n3
!    do j = 1,n2
!      do i = 1,n1
!        C(i,j,k) =  A(i,1,k) * B(1,j) &
!                  + A(i,2,k) * B(2,j) &
!                  + A(i,3,k) * B(3,j) &
!                  + A(i,4,k) * B(4,j) &
!                  + A(i,5,k) * B(5,j) &
!                  + A(i,6,k) * B(6,j) &
!                  + A(i,7,k) * B(7,j) &
!                  + A(i,8,k) * B(8,j)
!      enddo
!    enddo
!  enddo
!
!  end subroutine mxm8_3dmat_single

  end subroutine compute_forces_outer_core_Dev

