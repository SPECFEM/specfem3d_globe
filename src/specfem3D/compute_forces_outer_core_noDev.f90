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

  subroutine compute_forces_outer_core_noDev(timeval,deltat,two_omega_earth, &
                                             NSPEC_ROT,NGLOB, &
                                             A_array_rotation,B_array_rotation, &
                                             A_array_rotation_lddrk,B_array_rotation_lddrk, &
                                             displfluid,accelfluid, &
                                             div_displfluid,iphase, &
                                             pgrav_outer_core)

  use constants_solver

  use specfem_par, only: &
    hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
    gravity_pre_store => gravity_pre_store_outer_core, &
    MOVIE_VOLUME, &
    USE_LDDRK,istage

  use specfem_par_outercore, only: &
    xix => xix_outer_core,xiy => xiy_outer_core,xiz => xiz_outer_core, &
    etax => etax_outer_core,etay => etay_outer_core,etaz => etaz_outer_core, &
    gammax => gammax_outer_core,gammay => gammay_outer_core,gammaz => gammaz_outer_core, &
    ibool => ibool_outer_core, &
    phase_ispec_inner => phase_ispec_inner_outer_core, &
    nspec_outer => nspec_outer_outer_core, &
    nspec_inner => nspec_inner_outer_core

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

  ! inner/outer element run flag
  integer,intent(in) :: iphase

  ! full gravity
  real(kind=CUSTOM_REAL), dimension(NGLOB),intent(in) :: pgrav_outer_core

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3
  ! for gravity
  real(kind=CUSTOM_REAL) :: vec_x,vec_y,vec_z
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: gravity_term
  ! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rotation,B_rotation, &
       ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer ispec,iglob
  integer i,j,k,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l,sum_terms

!  integer :: computed_elements
  integer :: num_elements,ispec_p

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************

  if (MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE > 1 .and. (iphase == 1)) then
    div_displfluid(:,:,:,:) = 0._CUSTOM_REAL
  endif

!  computed_elements = 0
  if (iphase == 1) then
    ! outer elements
    num_elements = nspec_outer
  else
    ! inner elements
    num_elements = nspec_inner
  endif

  do ispec_p = 1,num_elements

    ispec = phase_ispec_inner(ispec_p,iphase)

    ! only compute element which belong to current phase (inner or outer elements)

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l = 1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + displfluid(ibool(l,j,k,ispec)) * hprime_xx(i,l)
            tempx2l = tempx2l + displfluid(ibool(i,l,k,ispec)) * hprime_yy(j,l)
            tempx3l = tempx3l + displfluid(ibool(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          ! get derivatives of potential with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! compute the Jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          dpotentialdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          dpotentialdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          dpotentialdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          ! compute contribution of rotation and add to gradient of potential
          ! this term has no Z component
          if (ROTATION_VAL) then

            ! store the source for the Euler scheme for A_rotation and B_rotation
            two_omega_deltat = deltat * two_omega_earth

            cos_two_omega_t = cos(two_omega_earth*timeval)
            sin_two_omega_t = sin(two_omega_earth*timeval)

            ! time step deltat of Euler scheme is included in the source
            source_euler_A(i,j,k) = two_omega_deltat &
                  * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
            source_euler_B(i,j,k) = two_omega_deltat &
                  * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

            A_rotation = A_array_rotation(i,j,k,ispec)
            B_rotation = B_array_rotation(i,j,k,ispec)

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
            iglob = ibool(i,j,k,ispec)

            vec_x = gravity_pre_store(1,iglob)
            vec_y = gravity_pre_store(2,iglob)
            vec_z = gravity_pre_store(3,iglob)

            ! adding (chi/rho)grad(rho)
            dpotentialdx_with_rot = dpotentialdx_with_rot + displfluid(iglob) * vec_x
            dpotentialdy_with_rot = dpotentialdy_with_rot + displfluid(iglob) * vec_y
            dpotentialdzl = dpotentialdzl + displfluid(iglob) * vec_z


         else
            ! if gravity is turned on

            ! compute divergence of displacement
            ! precompute and store gravity term

            ! use mesh coordinates to get theta and phi
            ! x y z contain r theta phi
            iglob = ibool(i,j,k,ispec)

            ! get g, rho and dg/dr=dg
            ! Cartesian components of the gravitational acceleration
            ! integrate and multiply by rho / Kappa
            vec_x = gravity_pre_store(1,iglob)
            vec_y = gravity_pre_store(2,iglob)
            vec_z = gravity_pre_store(3,iglob)

            gravity_term(i,j,k) = jacobianl * wgll_cube(i,j,k) &
                                  * (dpotentialdx_with_rot * vec_x &
                                   + dpotentialdy_with_rot * vec_y &
                                   + dpotentialdzl         * vec_z)

            ! full gravity contribution
            if (FULL_GRAVITY_VAL .and. .not. DISCARD_GCONTRIB) then
              gravity_term(i,j,k) = gravity_term(i,j,k) &
                    - gravity_rho_g_over_kappa(iglob) * jacobianl * wgll_cube(i,j,k) * pgrav_outer_core(iglob)
            endif

            ! divergence of displacement field with gravity on
            ! note: these calculations are only considered for SIMULATION_TYPE == 1 .and. SAVE_FORWARD
            !          and one has set MOVIE_VOLUME_TYPE == 4 when MOVIE_VOLUME is .true.;
            !         in case of SIMULATION_TYPE == 3, it gets overwritten by compute_kernels_outer_core()
            if (MOVIE_VOLUME .and. NSPEC_OUTER_CORE_3DMOVIE > 1) then
              div_displfluid(i,j,k,ispec) = dpotentialdx_with_rot * vec_x &
                                          + dpotentialdy_with_rot * vec_y &
                                          + dpotentialdzl * vec_z
            endif

          endif

          tempx1(i,j,k) = jacobianl*(xixl*dpotentialdx_with_rot + xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
          tempx2(i,j,k) = jacobianl*(etaxl*dpotentialdx_with_rot + etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
          tempx3(i,j,k) = jacobianl*(gammaxl*dpotentialdx_with_rot + gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

        enddo
      enddo
    enddo

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l = 1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + tempx1(l,j,k) * hprimewgll_xx(l,i)
            tempx2l = tempx2l + tempx2(i,l,k) * hprimewgll_yy(l,j)
            tempx3l = tempx3l + tempx3(i,j,l) * hprimewgll_zz(l,k)
          enddo

          ! sum contributions from each element to the global mesh and add gravity term
          sum_terms = - (wgllwgll_yz(j,k)*tempx1l + wgllwgll_xz(i,k)*tempx2l + wgllwgll_xy(i,j)*tempx3l)
          if (GRAVITY_VAL) sum_terms = sum_terms + gravity_term(i,j,k)

          accelfluid(ibool(i,j,k,ispec)) = accelfluid(ibool(i,j,k,ispec)) + sum_terms

        enddo
      enddo
    enddo

    ! update rotation term with Euler scheme
    if (ROTATION_VAL) then
      if (USE_LDDRK) then
        ! use the source saved above
        A_array_rotation_lddrk(:,:,:,ispec) = ALPHA_LDDRK(istage) * A_array_rotation_lddrk(:,:,:,ispec) + source_euler_A(:,:,:)
        A_array_rotation(:,:,:,ispec) = A_array_rotation(:,:,:,ispec) + BETA_LDDRK(istage) * A_array_rotation_lddrk(:,:,:,ispec)

        B_array_rotation_lddrk(:,:,:,ispec) = ALPHA_LDDRK(istage) * B_array_rotation_lddrk(:,:,:,ispec) + source_euler_B(:,:,:)
        B_array_rotation(:,:,:,ispec) = B_array_rotation(:,:,:,ispec) + BETA_LDDRK(istage) * B_array_rotation_lddrk(:,:,:,ispec)
      else
        ! use the source saved above
        A_array_rotation(:,:,:,ispec) = A_array_rotation(:,:,:,ispec) + source_euler_A(:,:,:)
        B_array_rotation(:,:,:,ispec) = B_array_rotation(:,:,:,ispec) + source_euler_B(:,:,:)
      endif
    endif

  enddo   ! spectral element loop

  end subroutine compute_forces_outer_core_noDev

