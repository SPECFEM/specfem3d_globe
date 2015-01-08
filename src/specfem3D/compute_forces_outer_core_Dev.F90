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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"



  subroutine compute_forces_outer_core_Dev(timeval,deltat,two_omega_earth, &
                                           NSPEC,NGLOB, &
                                           A_array_rotation,B_array_rotation, &
                                           A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                           displfluid,accelfluid, &
                                           div_displfluid,phase_is_inner)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver

  use specfem_par,only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT, &
    wgll_cube, &
    minus_rho_g_over_kappa_fluid,d_ln_density_dr_table, &
    MOVIE_VOLUME, &
    USE_LDDRK,istage

  use specfem_par_outercore,only: &
    xstore => xstore_outer_core,ystore => ystore_outer_core,zstore => zstore_outer_core, &
    xix => xix_outer_core,xiy => xiy_outer_core,xiz => xiz_outer_core, &
    etax => etax_outer_core,etay => etay_outer_core,etaz => etaz_outer_core, &
    gammax => gammax_outer_core,gammay => gammay_outer_core,gammaz => gammaz_outer_core, &
    ibool => ibool_outer_core, &
    phase_ispec_inner => phase_ispec_inner_outer_core, &
    nspec_outer => nspec_outer_outer_core, &
    nspec_inner => nspec_inner_outer_core

  use specfem_par,only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D

  implicit none

  integer :: NSPEC,NGLOB

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) timeval,deltat,two_omega_earth

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    A_array_rotation,B_array_rotation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    A_array_rotation_lddrk,B_array_rotation_lddrk

  ! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB) :: displfluid,accelfluid

  ! divergence of displacement
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_3DMOVIE) :: div_displfluid

  ! inner/outer element run flag
  logical :: phase_is_inner

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sum_terms

  ! for gravity
  integer :: int_radius
  double precision :: radius,theta,phi,gxl,gyl,gzl
  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: gravity_term

  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: temp_gxl,temp_gyl,temp_gzl

  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rotation,B_rotation, &
       ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer :: ispec,iglob
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl

  ! Deville
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3


  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    displ_times_grad_x_ln_rho,displ_times_grad_y_ln_rho,displ_times_grad_z_ln_rho

  integer :: num_elements,ispec_p
  integer :: iphase

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif
  real(kind=CUSTOM_REAL), dimension(NSTAGE) :: MYALPHA_LDDRK,MYBETA_LDDRK

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

  if (MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE /= 1 .and. (.not. phase_is_inner)) then
    div_displfluid(:,:,:,:) = 0._CUSTOM_REAL
  endif

! computed_elements = 0
  if (.not. phase_is_inner) then
    iphase = 1
    num_elements = nspec_outer
  else
    iphase = 2
    num_elements = nspec_inner
  endif

  ! note: this is an OpenMP work-around for gfortran
  ! gfortran has problems with parameters inside a DEFAULT(NONE) OpenMP block
  ! for a fix see https://gcc.gnu.org/ml/fortran/2014-10/msg00064.html
  !
  ! we use local variables to cirumvent it, another possibility would be to use DEFAULT(SHARED)
  !
  ! FIRSTPRIVATE declaration is chosen based on suggestion for performance optimization:
  ! see https://docs.oracle.com/cd/E19059-01/stud.10/819-0501/7_tuning.html

  if (USE_LDDRK) then
    MYALPHA_LDDRK(:) = ALPHA_LDDRK(:)
    MYBETA_LDDRK(:) = BETA_LDDRK(:)
  endif

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED( &
!$OMP num_elements, phase_ispec_inner, iphase, ibool, displfluid, xstore, ystore, zstore, &
!$OMP d_ln_density_dr_table, hprime_xx, hprime_xxT, xix, xiy, xiz,  etax, etay, etaz, &
!$OMP gammax, gammay, gammaz, deltat, two_omega_earth, timeval, A_array_rotation, B_array_rotation, &
!$OMP minus_rho_g_over_kappa_fluid, wgll_cube, MOVIE_VOLUME, hprimewgll_xxT, hprimewgll_xx, &
!$OMP wgllwgll_yz_3D, wgllwgll_xz_3D, wgllwgll_xy_3D, accelfluid, USE_LDDRK, A_array_rotation_lddrk, &
!$OMP istage, B_array_rotation_lddrk, div_displfluid ) &
!$OMP PRIVATE( &
!$OMP ispec_p, ispec, iglob, dummyx_loc, radius, theta, phi, &
!$OMP cos_theta, sin_theta, cos_phi, sin_phi, int_radius, &
!$OMP displ_times_grad_x_ln_rho, displ_times_grad_y_ln_rho, displ_times_grad_z_ln_rho, &
!$OMP temp_gxl, temp_gyl, temp_gzl, tempx2, &
!$OMP xixl, xiyl, xizl, etaxl, etayl, etazl, gammaxl, gammayl, gammazl, jacobianl, &
!$OMP dpotentialdxl, tempx1, tempx3, dpotentialdyl, dpotentialdzl, two_omega_deltat, cos_two_omega_t, &
!$OMP sin_two_omega_t, source_euler_A, source_euler_B, A_rotation, B_rotation, ux_rotation, uy_rotation, &
!$OMP dpotentialdx_with_rot, dpotentialdy_with_rot, gxl, gyl, gzl, gravity_term, &
!$OMP sum_terms, newtempx1, newtempx3, newtempx2 ) &
!$OMP FIRSTPRIVATE( MYALPHA_LDDRK,MYBETA_LDDRK )

!$OMP DO SCHEDULE(GUIDED)
  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner(ispec_p,iphase)

    ! only compute element which belong to current phase (inner or outer elements)

    DO_LOOP_IJK

      iglob = ibool(INDEX_IJK,ispec)

      ! stores "displacement"
      dummyx_loc(INDEX_IJK) = displfluid(iglob)

      ! pre-computes factors
      ! use mesh coordinates to get theta and phi
      ! x y z contain r theta phi
      radius = dble(xstore(iglob))
      theta = dble(ystore(iglob))
      phi = dble(zstore(iglob))

      cos_theta = dcos(theta)
      sin_theta = dsin(theta)
      cos_phi = dcos(phi)
      sin_phi = dsin(phi)

      int_radius = nint(radius * R_EARTH_KM * 10.d0)

      if (.not. GRAVITY_VAL) then
        ! grad(rho)/rho in Cartesian components
        displ_times_grad_x_ln_rho(INDEX_IJK) = dummyx_loc(INDEX_IJK) &
              * sngl(sin_theta * cos_phi * d_ln_density_dr_table(int_radius))
        displ_times_grad_y_ln_rho(INDEX_IJK) = dummyx_loc(INDEX_IJK) &
              * sngl(sin_theta * sin_phi * d_ln_density_dr_table(int_radius))
        displ_times_grad_z_ln_rho(INDEX_IJK) = dummyx_loc(INDEX_IJK) &
              * sngl(cos_theta * d_ln_density_dr_table(int_radius))
      else
        ! Cartesian components of the gravitational acceleration
        ! integrate and multiply by rho / Kappa
        temp_gxl(INDEX_IJK) = sin_theta*cos_phi
        temp_gyl(INDEX_IJK) = sin_theta*sin_phi
        temp_gzl(INDEX_IJK) = cos_theta
      endif

    ENDDO_LOOP_IJK

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for tempx1,..
    call mxm5_single(hprime_xx,m1,dummyx_loc,tempx1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm5_3dmat_single(dummyx_loc,m1,hprime_xxT,m1,tempx2,NGLLX)
    ! computes 3. matrix multiplication for tempx1,..
    call mxm5_single(dummyx_loc,m2,hprime_xxT,tempx3,m1)

    DO_LOOP_IJK

      ! get derivatives of velocity potential with respect to x, y and z
      xixl = xix(INDEX_IJK,ispec)
      xiyl = xiy(INDEX_IJK,ispec)
      xizl = xiz(INDEX_IJK,ispec)
      etaxl = etax(INDEX_IJK,ispec)
      etayl = etay(INDEX_IJK,ispec)
      etazl = etaz(INDEX_IJK,ispec)
      gammaxl = gammax(INDEX_IJK,ispec)
      gammayl = gammay(INDEX_IJK,ispec)
      gammazl = gammaz(INDEX_IJK,ispec)

      ! compute the Jacobian
      jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))

      dpotentialdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
      dpotentialdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
      dpotentialdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

      ! compute contribution of rotation and add to gradient of potential
      ! this term has no Z component
      if (ROTATION_VAL) then

        ! store the source for the Euler scheme for A_rotation and B_rotation
        two_omega_deltat = deltat * two_omega_earth

        cos_two_omega_t = cos(two_omega_earth*timeval)
        sin_two_omega_t = sin(two_omega_earth*timeval)

        ! time step deltat of Euler scheme is included in the source
        source_euler_A(INDEX_IJK) = two_omega_deltat &
              * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
        source_euler_B(INDEX_IJK) = two_omega_deltat &
              * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

        A_rotation = A_array_rotation(INDEX_IJK,ispec)
        B_rotation = B_array_rotation(INDEX_IJK,ispec)

        ux_rotation =   A_rotation*cos_two_omega_t + B_rotation*sin_two_omega_t
        uy_rotation = - A_rotation*sin_two_omega_t + B_rotation*cos_two_omega_t

        dpotentialdx_with_rot = dpotentialdxl + ux_rotation
        dpotentialdy_with_rot = dpotentialdyl + uy_rotation

      else

        dpotentialdx_with_rot = dpotentialdxl
        dpotentialdy_with_rot = dpotentialdyl

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

        ! use mesh coordinates to get theta and phi
        ! x y z contain r theta phi
        dpotentialdx_with_rot = dpotentialdx_with_rot + displ_times_grad_x_ln_rho(INDEX_IJK)
        dpotentialdy_with_rot = dpotentialdy_with_rot + displ_times_grad_y_ln_rho(INDEX_IJK)
        dpotentialdzl = dpotentialdzl + displ_times_grad_z_ln_rho(INDEX_IJK)

     else
        ! if gravity is turned on

        ! compute divergence of displacement
        gxl = temp_gxl(INDEX_IJK)
        gyl = temp_gyl(INDEX_IJK)
        gzl = temp_gzl(INDEX_IJK)

        ! distinguish between single and double precision for reals
        gravity_term(INDEX_IJK) = &
                real(minus_rho_g_over_kappa_fluid(int_radius) &
                   * dble(jacobianl) * wgll_cube(INDEX_IJK) &
                   * (dble(dpotentialdx_with_rot) * gxl &
                    + dble(dpotentialdy_with_rot) * gyl &
                    + dble(dpotentialdzl)         * gzl), kind=CUSTOM_REAL)

        ! divergence of displacement field with gravity on
        ! note: these calculations are only considered for SIMULATION_TYPE == 1 .and. SAVE_FORWARD
        !          and one has set MOVIE_VOLUME_TYPE == 4 when MOVIE_VOLUME is .true.;
        !         in case of SIMULATION_TYPE == 3, it gets overwritten by compute_kernels_outer_core()
        if (MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE /= 1) then
          div_displfluid(INDEX_IJK,ispec) = &
                    real(minus_rho_g_over_kappa_fluid(int_radius) &
                       * (dpotentialdx_with_rot * gxl &
                        + dpotentialdy_with_rot * gyl &
                        + dpotentialdzl         * gzl), kind=CUSTOM_REAL)
        endif

      endif

      tempx1(INDEX_IJK) = jacobianl*(xixl*dpotentialdx_with_rot &
                               + xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
      tempx2(INDEX_IJK) = jacobianl*(etaxl*dpotentialdx_with_rot &
                               + etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
      tempx3(INDEX_IJK) = jacobianl*(gammaxl*dpotentialdx_with_rot &
                               + gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

    ENDDO_LOOP_IJK

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

    ! computes 1. matrix multiplication for newtempx1,..
    call mxm5_single(hprimewgll_xxT,m1,tempx1,newtempx1,m2)
    ! computes 2. matrix multiplication for tempx2,..
    call mxm5_3dmat_single(tempx2,m1,hprimewgll_xx,m1,newtempx2,NGLLX)
    ! computes 3. matrix multiplication for newtempx3,..
    call mxm5_single(tempx3,m2,hprimewgll_xx,newtempx3,m1)

    ! sum contributions from each element to the global mesh and add gravity term

    DO_LOOP_IJK

          sum_terms(INDEX_IJK) = - ( wgllwgll_yz_3D(INDEX_IJK)*newtempx1(INDEX_IJK) &
                                   + wgllwgll_xz_3D(INDEX_IJK)*newtempx2(INDEX_IJK) &
                                   + wgllwgll_xy_3D(INDEX_IJK)*newtempx3(INDEX_IJK))

    ENDDO_LOOP_IJK

    ! adds gravity
    if (GRAVITY_VAL) then

      DO_LOOP_IJK

        sum_terms(INDEX_IJK) = sum_terms(INDEX_IJK) + gravity_term(INDEX_IJK)

      ENDDO_LOOP_IJK

    endif

    ! updates acceleration

#ifdef FORCE_VECTORIZATION
#ifndef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP CRITICAL
#endif
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
! IBM, Portland PGI, and Intel and Cray syntax (Intel and Cray are the same)
!IBM* ASSERT (NODEPS)
!pgi$ ivdep
!DIR$ IVDEP
    do ijk = 1,NGLLCUBE
#else
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
#endif
          iglob = ibool(INDEX_IJK,ispec)
#ifdef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP ATOMIC
#endif
          accelfluid(iglob) = accelfluid(iglob) + sum_terms(INDEX_IJK)

#ifdef FORCE_VECTORIZATION
    enddo
#ifndef USE_OPENMP_ATOMIC_INSTEAD_OF_CRITICAL
!$OMP END CRITICAL
#endif
#else
        enddo
      enddo
    enddo
#endif

    ! update rotation term with Euler scheme

    if (ROTATION_VAL) then

      if (USE_LDDRK) then

        ! use the source saved above
        DO_LOOP_IJK

          A_array_rotation_lddrk(INDEX_IJK,ispec) = MYALPHA_LDDRK(istage) * A_array_rotation_lddrk(INDEX_IJK,ispec) &
                                                    + source_euler_A(INDEX_IJK)
          A_array_rotation(INDEX_IJK,ispec) = A_array_rotation(INDEX_IJK,ispec) &
                                              + MYBETA_LDDRK(istage) * A_array_rotation_lddrk(INDEX_IJK,ispec)

          B_array_rotation_lddrk(INDEX_IJK,ispec) = MYALPHA_LDDRK(istage) * B_array_rotation_lddrk(INDEX_IJK,ispec) &
                                                    + source_euler_B(INDEX_IJK)
          B_array_rotation(INDEX_IJK,ispec) = B_array_rotation(INDEX_IJK,ispec) &
                                              + MYBETA_LDDRK(istage) * B_array_rotation_lddrk(INDEX_IJK,ispec)

        ENDDO_LOOP_IJK

      else

        DO_LOOP_IJK

          A_array_rotation(INDEX_IJK,ispec) = A_array_rotation(INDEX_IJK,ispec) + source_euler_A(INDEX_IJK)
          B_array_rotation(INDEX_IJK,ispec) = B_array_rotation(INDEX_IJK,ispec) + source_euler_B(INDEX_IJK)

        ENDDO_LOOP_IJK

      endif
    endif

  enddo   ! spectral element loop
!$OMP enddo
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

! 2-dimensional arrays (25,5)/(5,25)

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n3
  real(kind=CUSTOM_REAL),dimension(n1,5) :: A
  real(kind=CUSTOM_REAL),dimension(5,n3) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n3) :: C

  ! local parameters
  integer :: i,j

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


!--------------------------------------------------------------------------------------------

  subroutine mxm5_3dmat_single(A,n1,B,n2,C,n3)

! 3-dimensional arrays (5,5,5) for A and C

  use constants_solver,only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: n1,n2,n3
  real(kind=CUSTOM_REAL),dimension(n1,5,n3) :: A
  real(kind=CUSTOM_REAL),dimension(5,n2) :: B
  real(kind=CUSTOM_REAL),dimension(n1,n2,n3) :: C

  ! local parameters
  integer :: i,j,k

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

  end subroutine compute_forces_outer_core_Dev

