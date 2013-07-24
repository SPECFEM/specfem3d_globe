!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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

  subroutine compute_forces_outer_core(time,deltat,two_omega_earth, &
                          A_array_rotation,B_array_rotation, &
                          d_ln_density_dr_table, &
                          minus_rho_g_over_kappa_fluid,displfluid,accelfluid, &
                          div_displfluid, &
                          xstore,ystore,zstore, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
          is_on_a_slice_edge_outer_core, &
          myrank,iproc_xi,iproc_eta,ichunk,addressing, &
          iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
          npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
          iboolfaces_outer_core,iboolcorner_outer_core, &
          iprocfrom_faces,iprocto_faces, &
          iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
          buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
          buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar,iphase,icall, &
                          hprime_xx,hprime_yy,hprime_zz, &
                          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
                          ibool,MOVIE_VOLUME,&
                          istage,A_array_rotation_lddrk,B_array_rotation_lddrk)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displfluid,accelfluid

! divergence of displacement
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: div_displfluid

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec_outer_core) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_outer_core) :: xix,xiy,xiz, &
                      etax,etay,etaz,gammax,gammay,gammaz

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  logical :: MOVIE_VOLUME

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3

! for gravity
  integer int_radius
  double precision radius,theta,phi,gxl,gyl,gzl
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision, dimension(NRAD_GRAVITY) :: minus_rho_g_over_kappa_fluid
  double precision, dimension(NRAD_GRAVITY) :: d_ln_density_dr_table
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: gravity_term
  real(kind=CUSTOM_REAL), dimension(nglob_outer_core) :: xstore,ystore,zstore

! for the Euler scheme for rotation
  real(kind=CUSTOM_REAL) time,deltat,two_omega_earth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION) :: &
    A_array_rotation,B_array_rotation

  real(kind=CUSTOM_REAL) two_omega_deltat,cos_two_omega_t,sin_two_omega_t,A_rotation,B_rotation, &
       ux_rotation,uy_rotation,dpotentialdx_with_rot,dpotentialdy_with_rot
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: source_euler_A,source_euler_B

  integer ispec,iglob
  integer i,j,k,l

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) dpotentialdxl,dpotentialdyl,dpotentialdzl
  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l,sum_terms

  double precision grad_x_ln_rho,grad_y_ln_rho,grad_z_ln_rho

! this for non blocking MPI
  integer :: ichunk,iproc_xi,iproc_eta,myrank

  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_OC) :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_OC) :: iboolleft_eta_outer_core,iboolright_eta_outer_core

  integer npoin2D_faces_outer_core(NUMFACES_SHARED)
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_outer_core,npoin2D_eta_outer_core

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: iprocfrom_faces,iprocto_faces

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

! indirect addressing for each message for faces and corners of the chunks
! a given slice can belong to at most one corner and at most two faces
  integer, dimension(NGLOB2DMAX_XY_OC_VAL,NUMFACES_SHARED) :: iboolfaces_outer_core

! buffers for send and receive between faces of the slices and the chunks
! we use the same buffers to assemble scalars and vectors because vectors are
! always three times bigger and therefore scalars can use the first part
! of the vector buffer in memory even if it has an additional index here
  integer :: npoin2D_max_all_CM_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED) :: buffer_send_faces,buffer_received_faces

  integer, dimension(NGLOB1D_RADIAL_OC,NUMCORNERS_SHARED) :: iboolcorner_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL_OC) :: buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar

  logical, dimension(NSPEC_OUTER_CORE) :: is_on_a_slice_edge_outer_core

  integer :: iphase,icall

  integer :: computed_elements

! for LDDRK
  integer :: istage
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION) :: &
            A_array_rotation_lddrk,B_array_rotation_lddrk

! ****************************************************
!   big loop over all spectral elements in the fluid
! ****************************************************
  if(istage == 1) then
    if (NSPEC_OUTER_CORE_ADJOINT /= 1 .and. SIMULATION_TYPE == 1 .and. icall == 1) div_displfluid(:,:,:,:) = 0._CUSTOM_REAL
  endif

  computed_elements = 0

  do ispec = 1,NSPEC_OUTER_CORE

! hide communications by computing the edges first
    if((icall == 2 .and. is_on_a_slice_edge_outer_core(ispec)) .or. &
       (icall == 1 .and. .not. is_on_a_slice_edge_outer_core(ispec))) cycle

! process the non-blocking communications every ELEMENTS_NONBLOCKING elements
    computed_elements = computed_elements + 1
    if (icall == 2 .and. mod(computed_elements,ELEMENTS_NONBLOCKING_OC) == 0 .and. iphase <= 7) &
      call assemble_MPI_scalar(myrank,accelfluid,NGLOB_OUTER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
            npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
            iboolfaces_outer_core,iboolcorner_outer_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES_VAL,NCORNERSCHUNKS_VAL, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL_OC, &
            NGLOB2DMAX_XMIN_XMAX_OC,NGLOB2DMAX_YMIN_YMAX_OC, &
            NGLOB2DMAX_XY_OC_VAL,NCHUNKS_VAL,iphase)

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + displfluid(ibool(l,j,k,ispec)) * hprime_xx(i,l)
            tempx2l = tempx2l + displfluid(ibool(i,l,k,ispec)) * hprime_yy(j,l)
            tempx3l = tempx3l + displfluid(ibool(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          ! get derivatives of velocity potential with respect to x, y and z
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          dpotentialdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          dpotentialdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          dpotentialdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          ! compute contribution of rotation and add to gradient of potential
          ! this term has no Z component
          if(ROTATION_VAL) then

            ! store the source for the Euler scheme for A_rotation and B_rotation
            two_omega_deltat = deltat * two_omega_earth

            cos_two_omega_t = cos(two_omega_earth*time)
            sin_two_omega_t = sin(two_omega_earth*time)

            ! time step deltat of Euler scheme is included in the source
            source_euler_A(i,j,k) = two_omega_deltat * (cos_two_omega_t * dpotentialdyl + sin_two_omega_t * dpotentialdxl)
            source_euler_B(i,j,k) = two_omega_deltat * (sin_two_omega_t * dpotentialdyl - cos_two_omega_t * dpotentialdxl)

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
            iglob = ibool(i,j,k,ispec)

            radius = dble(xstore(iglob))
            theta = dble(ystore(iglob))
            phi = dble(zstore(iglob))

            cos_theta = dcos(theta)
            sin_theta = dsin(theta)
            cos_phi = dcos(phi)
            sin_phi = dsin(phi)

            int_radius = nint(radius * R_EARTH_KM * 10.d0)

            ! grad(rho)/rho in Cartesian components
            grad_x_ln_rho = sin_theta * cos_phi * d_ln_density_dr_table(int_radius)
            grad_y_ln_rho = sin_theta * sin_phi * d_ln_density_dr_table(int_radius)
            grad_z_ln_rho = cos_theta * d_ln_density_dr_table(int_radius)

            ! adding (chi/rho)grad(rho)
            dpotentialdx_with_rot = dpotentialdx_with_rot + displfluid(iglob) * grad_x_ln_rho
            dpotentialdy_with_rot = dpotentialdy_with_rot + displfluid(iglob) * grad_y_ln_rho
            dpotentialdzl = dpotentialdzl + displfluid(iglob) * grad_z_ln_rho


         else  ! if gravity is turned on

            ! compute divergence of displacment
            ! precompute and store gravity term

            ! use mesh coordinates to get theta and phi
            ! x y z contain r theta phi
            iglob = ibool(i,j,k,ispec)

            radius = dble(xstore(iglob))
            theta = dble(ystore(iglob))
            phi = dble(zstore(iglob))

            cos_theta = dcos(theta)
            sin_theta = dsin(theta)
            cos_phi = dcos(phi)
            sin_phi = dsin(phi)

            ! get g, rho and dg/dr=dg
            ! spherical components of the gravitational acceleration
            ! for efficiency replace with lookup table every 100 m in radial direction
            int_radius = nint(radius * R_EARTH_KM * 10.d0)

            ! Cartesian components of the gravitational acceleration
            ! integrate and multiply by rho / Kappa
            gxl = sin_theta*cos_phi
            gyl = sin_theta*sin_phi
            gzl = cos_theta

            ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then
              gravity_term(i,j,k) = &
                sngl(minus_rho_g_over_kappa_fluid(int_radius) * &
                dble(jacobianl) * wgll_cube(i,j,k) * &
               (dble(dpotentialdx_with_rot) * gxl + &
                dble(dpotentialdy_with_rot) * gyl + dble(dpotentialdzl) * gzl))
            else
              gravity_term(i,j,k) = minus_rho_g_over_kappa_fluid(int_radius) * &
                 jacobianl * wgll_cube(i,j,k) * (dpotentialdx_with_rot * gxl + &
                 dpotentialdy_with_rot * gyl + dpotentialdzl * gzl)
            endif

            if(istage == 1) then
              ! divergence of displacement field with gravity on
              ! note: these calculations are only considered for SIMULATION_TYPE == 1 .and. SAVE_FORWARD
              !          and one has set MOVIE_VOLUME_TYPE == 4 when MOVIE_VOLUME is .true.;
              !         in case of SIMULATION_TYPE == 3, it gets overwritten by compute_kernels_outer_core()
              if (NSPEC_OUTER_CORE_ADJOINT /= 1 .and. SIMULATION_TYPE == 1 .and. MOVIE_VOLUME) then
                div_displfluid(i,j,k,ispec) =  &
                   minus_rho_g_over_kappa_fluid(int_radius) * (dpotentialdx_with_rot * gxl + &
                   dpotentialdy_with_rot * gyl + dpotentialdzl * gzl)
              endif
            endif

          endif

          tempx1(i,j,k) = jacobianl*(xixl*dpotentialdx_with_rot + xiyl*dpotentialdy_with_rot + xizl*dpotentialdzl)
          tempx2(i,j,k) = jacobianl*(etaxl*dpotentialdx_with_rot + etayl*dpotentialdy_with_rot + etazl*dpotentialdzl)
          tempx3(i,j,k) = jacobianl*(gammaxl*dpotentialdx_with_rot + gammayl*dpotentialdy_with_rot + gammazl*dpotentialdzl)

        enddo
      enddo
    enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            !!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo
            tempx1l = tempx1l + tempx1(l,j,k) * hprimewgll_xx(l,i)
            tempx2l = tempx2l + tempx2(i,l,k) * hprimewgll_yy(l,j)
            tempx3l = tempx3l + tempx3(i,j,l) * hprimewgll_zz(l,k)
          enddo

          ! sum contributions from each element to the global mesh and add gravity term
          sum_terms = - (wgllwgll_yz(j,k)*tempx1l + wgllwgll_xz(i,k)*tempx2l + wgllwgll_xy(i,j)*tempx3l)
          if(GRAVITY_VAL) sum_terms = sum_terms + gravity_term(i,j,k)

          accelfluid(ibool(i,j,k,ispec)) = accelfluid(ibool(i,j,k,ispec)) + sum_terms

        enddo
      enddo
    enddo

    ! update rotation term with Euler scheme
    if(ROTATION_VAL) then
      if(USE_LDDRK)then
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

  end subroutine compute_forces_outer_core

