!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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



  subroutine compute_forces_outer_core_Dev(time,deltat,two_omega_earth, &
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
  real(kind=CUSTOM_REAL) time,deltat,two_omega_earth

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
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: newtempx1,newtempx2,newtempx3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: sum_terms

  ! for gravity
  integer :: int_radius
  double precision :: radius,theta,phi,gxl,gyl,gzl
  double precision :: cos_theta,sin_theta,cos_phi,sin_phi
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: gravity_term
  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) :: two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rotation,B_rotation, &
       ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer :: ispec,iglob
  integer :: i,j,k
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: dpotentialdxl,dpotentialdyl,dpotentialdzl

  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(newtempx1,E1_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: E1_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)

  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: temp_gxl,temp_gyl,temp_gzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    displ_times_grad_x_ln_rho,displ_times_grad_y_ln_rho,displ_times_grad_z_ln_rho

! integer :: computed_elements
  integer :: num_elements,ispec_p
  integer :: iphase

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

  if( MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE /= 1 .and. (.not. phase_is_inner) ) then
    div_displfluid(:,:,:,:) = 0._CUSTOM_REAL
  endif

! computed_elements = 0
  if( .not. phase_is_inner ) then
    iphase = 1
    num_elements = nspec_outer
  else
    iphase = 2
    num_elements = nspec_inner
  endif

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

      if( .not. GRAVITY_VAL ) then
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
    do j=1,m2
      do i=1,m1
        C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                                hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                                hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                                hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                                hprime_xx(i,5)*B1_m1_m2_5points(5,j)
      enddo
    enddo

    do k = 1,NGLLX
      do j=1,m1
        do i=1,m1
          tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                          dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                          dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                          dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                          dummyx_loc(i,5,k)*hprime_xxT(5,j)
        enddo
      enddo
    enddo

    do j=1,m1
      do i=1,m2
        C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                    A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                    A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                    A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                    A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
      enddo
    enddo

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

      ! compute the jacobian
      jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))

      dpotentialdxl = xixl*tempx1(INDEX_IJK) + etaxl*tempx2(INDEX_IJK) + gammaxl*tempx3(INDEX_IJK)
      dpotentialdyl = xiyl*tempx1(INDEX_IJK) + etayl*tempx2(INDEX_IJK) + gammayl*tempx3(INDEX_IJK)
      dpotentialdzl = xizl*tempx1(INDEX_IJK) + etazl*tempx2(INDEX_IJK) + gammazl*tempx3(INDEX_IJK)

      ! compute contribution of rotation and add to gradient of potential
      ! this term has no Z component
      if(ROTATION_VAL) then

        ! store the source for the Euler scheme for A_rotation and B_rotation
        two_omega_deltat = deltat * two_omega_earth

        cos_two_omega_t = cos(two_omega_earth*time)
        sin_two_omega_t = sin(two_omega_earth*time)

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
      if(.not. GRAVITY_VAL) then

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

        ! compute divergence of displacment
        gxl = temp_gxl(INDEX_IJK)
        gyl = temp_gyl(INDEX_IJK)
        gzl = temp_gzl(INDEX_IJK)

        ! distinguish between single and double precision for reals
        if(CUSTOM_REAL == SIZE_REAL) then
          gravity_term(INDEX_IJK) = &
                  sngl( minus_rho_g_over_kappa_fluid(int_radius) &
                  * dble(jacobianl) * wgll_cube(INDEX_IJK) &
                  * (dble(dpotentialdx_with_rot) * gxl  &
                     + dble(dpotentialdy_with_rot) * gyl &
                     + dble(dpotentialdzl) * gzl) )
        else
          gravity_term(INDEX_IJK) = minus_rho_g_over_kappa_fluid(int_radius) * &
                    jacobianl * wgll_cube(INDEX_IJK) &
                    * (dpotentialdx_with_rot * gxl  &
                      + dpotentialdy_with_rot * gyl &
                      + dpotentialdzl * gzl)
        endif

        ! divergence of displacement field with gravity on
        ! note: these calculations are only considered for SIMULATION_TYPE == 1 .and. SAVE_FORWARD
        !          and one has set MOVIE_VOLUME_TYPE == 4 when MOVIE_VOLUME is .true.;
        !         in case of SIMULATION_TYPE == 3, it gets overwritten by compute_kernels_outer_core()
        if( MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE /= 1 ) then
          div_displfluid(INDEX_IJK,ispec) =  &
                    minus_rho_g_over_kappa_fluid(int_radius) &
                    * (dpotentialdx_with_rot * gxl &
                     + dpotentialdy_with_rot * gyl &
                     + dpotentialdzl * gzl)
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
    do j=1,m2
      do i=1,m1
        E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
                                hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
                                hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
                                hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
                                hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)
      enddo
    enddo

    do k = 1,NGLLX
      do j=1,m1
        do i=1,m1
          newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
                             tempx2(i,2,k)*hprimewgll_xx(2,j) + &
                             tempx2(i,3,k)*hprimewgll_xx(3,j) + &
                             tempx2(i,4,k)*hprimewgll_xx(4,j) + &
                             tempx2(i,5,k)*hprimewgll_xx(5,j)
        enddo
      enddo
    enddo

    do j=1,m1
      do i=1,m2
        E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                    C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                    C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                    C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                    C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
      enddo
    enddo

    ! sum contributions from each element to the global mesh and add gravity term

    DO_LOOP_IJK

          sum_terms(INDEX_IJK) = - ( wgllwgll_yz_3D(INDEX_IJK)*newtempx1(INDEX_IJK) &
                                   + wgllwgll_xz_3D(INDEX_IJK)*newtempx2(INDEX_IJK) &
                                   + wgllwgll_xy_3D(INDEX_IJK)*newtempx3(INDEX_IJK))

    ENDDO_LOOP_IJK

    ! adds gravity
    if(GRAVITY_VAL) then

      DO_LOOP_IJK

        sum_terms(INDEX_IJK) = sum_terms(INDEX_IJK) + gravity_term(INDEX_IJK)

      ENDDO_LOOP_IJK

    endif

    ! updates acceleration

#ifdef FORCE_VECTORIZATION
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
! IBM, Portland PGI, and Intel and Cray syntax (Intel and Cray are the same)
!IBM* ASSERT (NODEPS)
!pgi$ ivdep
!DIR$ IVDEP
    do ijk = 1,NGLLCUBE
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
#endif
          iglob = ibool(INDEX_IJK,ispec)
          accelfluid(iglob) = accelfluid(iglob) + sum_terms(INDEX_IJK)

#ifdef FORCE_VECTORIZATION
    enddo
#else
        enddo
      enddo
    enddo
#endif

    ! update rotation term with Euler scheme

    if(ROTATION_VAL) then

      if(USE_LDDRK) then

        ! use the source saved above
        DO_LOOP_IJK

          A_array_rotation_lddrk(INDEX_IJK,ispec) = ALPHA_LDDRK(istage) * A_array_rotation_lddrk(INDEX_IJK,ispec) &
                                                    + source_euler_A(INDEX_IJK)
          A_array_rotation(INDEX_IJK,ispec) = A_array_rotation(INDEX_IJK,ispec) &
                                              + BETA_LDDRK(istage) * A_array_rotation_lddrk(INDEX_IJK,ispec)

          B_array_rotation_lddrk(INDEX_IJK,ispec) = ALPHA_LDDRK(istage) * B_array_rotation_lddrk(INDEX_IJK,ispec) &
                                                    + source_euler_B(INDEX_IJK)
          B_array_rotation(INDEX_IJK,ispec) = B_array_rotation(INDEX_IJK,ispec) &
                                              + BETA_LDDRK(istage) * B_array_rotation_lddrk(INDEX_IJK,ispec)

        ENDDO_LOOP_IJK

      else

        DO_LOOP_IJK

          A_array_rotation(INDEX_IJK,ispec) = A_array_rotation(INDEX_IJK,ispec) + source_euler_A(INDEX_IJK)
          B_array_rotation(INDEX_IJK,ispec) = B_array_rotation(INDEX_IJK,ispec) + source_euler_B(INDEX_IJK)

        ENDDO_LOOP_IJK

      endif
    endif

  enddo   ! spectral element loop

  end subroutine compute_forces_outer_core_Dev

