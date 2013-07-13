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

  subroutine compute_kernels_crust_mantle(ibool_crust_mantle, &
                          rho_kl_crust_mantle,beta_kl_crust_mantle, &
                          alpha_kl_crust_mantle,cijkl_kl_crust_mantle, &
                          accel_crust_mantle,b_displ_crust_mantle, &
                          deltat,displ_crust_mantle,hprime_xx,hprime_xxT,&
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,ANISOTROPIC_KL,&
                          RECOMPUTE_STRAIN_DO_NOT_STORE,& 
                          epsilondev_crust_mantle,eps_trace_over_3_crust_mantle) 

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    rho_kl_crust_mantle, beta_kl_crust_mantle, alpha_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL) :: &
    cijkl_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     accel_crust_mantle,displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle

  real(kind=CUSTOM_REAL) deltat
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  logical :: ANISOTROPIC_KL

  logical :: RECOMPUTE_STRAIN_DO_NOT_STORE 
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: epsilondev_crust_mantle 
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: eps_trace_over_3_crust_mantle 

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(21) :: prod
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_matrix,b_epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_matrix,&
                                                          b_eps_trace_over_3_loc_matrix
  integer :: i,j,k,ispec,iglob

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE

    if(COMPUTE_AND_STORE_STRAIN .and. RECOMPUTE_STRAIN_DO_NOT_STORE)then  
        eps_trace_over_3_loc_matrix(:,:,:) = eps_trace_over_3_crust_mantle(:,:,:,ispec)
        epsilondev_loc_matrix(1,:,:,:) = epsilondev_crust_mantle(1,:,:,:,ispec)
        epsilondev_loc_matrix(2,:,:,:) = epsilondev_crust_mantle(2,:,:,:,ispec)
        epsilondev_loc_matrix(3,:,:,:) = epsilondev_crust_mantle(3,:,:,:,ispec)
        epsilondev_loc_matrix(4,:,:,:) = epsilondev_crust_mantle(4,:,:,:,ispec)
        epsilondev_loc_matrix(5,:,:,:) = epsilondev_crust_mantle(5,:,:,:,ispec)
    else
      call compute_element_strain_undo_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE,&
                                             displ_CRUST_MANTLE,ibool_crust_mantle,hprime_xx,hprime_xxT,&
                                             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                                             epsilondev_loc_matrix,eps_trace_over_3_loc_matrix)
    endif

    call compute_element_strain_undo_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE,&
                                             b_displ_CRUST_MANTLE,ibool_crust_mantle,hprime_xx,hprime_xxT,&
                                             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                                             b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)


    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)

          ! density kernel: see e.g. Tromp et al.(2005), equation (14)
          !                         b_displ_crust_mantle is the backward/reconstructed wavefield, that is s(x,t) in eq. (14),
          !                         accel_crust_mantle is the adjoint wavefield, that corresponds to s_dagger(x,T-t)
          !
          !                         note with respect to eq. (14) the second time derivative is applied to the
          !                         adjoint wavefield here rather than the backward/reconstructed wavefield.
          !                         this is a valid operation and the resultant kernel identical to the eq. (14).
          !
          !                         reason for this is that the adjoint wavefield is in general smoother
          !                         since the adjoint sources normally are obtained for filtered traces.
          !                         numerically, the time derivative by a finite-difference scheme should
          !                         behave better for smoother wavefields, thus containing less numerical artefacts.
          rho_kl_crust_mantle(i,j,k,ispec) =  rho_kl_crust_mantle(i,j,k,ispec) &
             + deltat * (accel_crust_mantle(1,iglob) * b_displ_crust_mantle(1,iglob) &
             + accel_crust_mantle(2,iglob) * b_displ_crust_mantle(2,iglob) &
             + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob) )

          epsilondev_loc(:) = epsilondev_loc_matrix(:,i,j,k)
          b_epsilondev_loc(:) = b_epsilondev_loc_matrix(:,i,j,k)

          ! For anisotropic kernels
          if (ANISOTROPIC_KL) then

            call compute_strain_product(prod,eps_trace_over_3_loc_matrix(i,j,k),epsilondev_loc, &
                                        b_eps_trace_over_3_loc_matrix(i,j,k),b_epsilondev_loc)
            cijkl_kl_crust_mantle(:,i,j,k,ispec) = cijkl_kl_crust_mantle(:,i,j,k,ispec) + deltat * prod(:)

          else

            ! kernel for shear modulus, see e.g. Tromp et al. (2005), equation (17)
            ! note: multiplication with 2*mu(x) will be done after the time loop
            beta_kl_crust_mantle(i,j,k,ispec) =  beta_kl_crust_mantle(i,j,k,ispec) &
               + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
               + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
               + 2.d0 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) + &
                epsilondev_loc(5)*b_epsilondev_loc(5)) )


            ! kernel for bulk modulus, see e.g. Tromp et al. (2005), equation (18)
            ! note: multiplication with kappa(x) will be done after the time loop
            alpha_kl_crust_mantle(i,j,k,ispec) = alpha_kl_crust_mantle(i,j,k,ispec) &
               + deltat * (9.d0 * eps_trace_over_3_loc_matrix(i,j,k) &
                                * b_eps_trace_over_3_loc_matrix(i,j,k))

          endif

        enddo
      enddo
    enddo
  enddo

  end subroutine compute_kernels_crust_mantle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_outer_core(ibool_outer_core, &
                        xix_outer_core,xiy_outer_core,xiz_outer_core, &
                        etax_outer_core,etay_outer_core,etaz_outer_core, &
                        gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                        hprime_xx,hprime_yy,hprime_zz, &
                        displ_outer_core,accel_outer_core, &
                        b_displ_outer_core,b_accel_outer_core, &
                        vector_accel_outer_core,vector_displ_outer_core, &
                        b_vector_displ_outer_core, &
                        div_displ_outer_core, &
                        rhostore_outer_core,kappavstore_outer_core, &
                        rho_kl_outer_core,alpha_kl_outer_core, &
                        deltat)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        xix_outer_core,xiy_outer_core,xiz_outer_core,&
        etax_outer_core,etay_outer_core,etaz_outer_core, &
        gammax_outer_core,gammay_outer_core,gammaz_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    displ_outer_core,accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_displ_outer_core,b_accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: vector_accel_outer_core,&
             vector_displ_outer_core, b_vector_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        rhostore_outer_core,kappavstore_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: &
    rho_kl_outer_core,alpha_kl_outer_core

  real(kind=CUSTOM_REAL) deltat

  ! local parameters
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,kappal
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: b_div_displ_outer_core

  integer :: i,j,k,l,ispec,iglob

  ! outer core -- compute the actual displacement and acceleration
  do ispec = 1, NSPEC_OUTER_CORE

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          xixl = xix_outer_core(i,j,k,ispec)
          xiyl = xiy_outer_core(i,j,k,ispec)
          xizl = xiz_outer_core(i,j,k,ispec)
          etaxl = etax_outer_core(i,j,k,ispec)
          etayl = etay_outer_core(i,j,k,ispec)
          etazl = etaz_outer_core(i,j,k,ispec)
          gammaxl = gammax_outer_core(i,j,k,ispec)
          gammayl = gammay_outer_core(i,j,k,ispec)
          gammazl = gammaz_outer_core(i,j,k,ispec)

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            tempx1l = tempx1l + b_displ_outer_core(ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
          enddo

          do l=1,NGLLY
            tempx2l = tempx2l + b_displ_outer_core(ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
          enddo

          do l=1,NGLLZ
            tempx3l = tempx3l + b_displ_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          b_vector_displ_outer_core(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          b_vector_displ_outer_core(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          b_vector_displ_outer_core(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            tempx1l = tempx1l + accel_outer_core(ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
          enddo

          do l=1,NGLLY
            tempx2l = tempx2l + accel_outer_core(ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
          enddo

          do l=1,NGLLZ
            tempx3l = tempx3l + accel_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          vector_accel_outer_core(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          vector_accel_outer_core(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          vector_accel_outer_core(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            tempx1l = tempx1l + displ_outer_core(ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
          enddo

          do l=1,NGLLY
            tempx2l = tempx2l + displ_outer_core(ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
          enddo

          do l=1,NGLLZ
            tempx3l = tempx3l + displ_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          vector_displ_outer_core(1,i,j,k) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          vector_displ_outer_core(2,i,j,k) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          vector_displ_outer_core(3,i,j,k) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          rho_kl_outer_core(i,j,k,ispec) = rho_kl_outer_core(i,j,k,ispec) &
             + deltat * dot_product(vector_accel_outer_core(:,i,j,k), b_vector_displ_outer_core(:,i,j,k))

          kappal = rhostore_outer_core(i,j,k,ispec)/kappavstore_outer_core(i,j,k,ispec)

          iglob = ibool_outer_core(i,j,k,ispec)

          div_displ_outer_core(i,j,k,ispec) =  kappal * accel_outer_core(iglob)
          b_div_displ_outer_core =  kappal * b_accel_outer_core(iglob)

          alpha_kl_outer_core(i,j,k,ispec) = alpha_kl_outer_core(i,j,k,ispec) &
             + deltat * div_displ_outer_core(i,j,k,ispec) * b_div_displ_outer_core

        enddo
      enddo
    enddo

  enddo

  end subroutine compute_kernels_outer_core

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_inner_core(ibool_inner_core, &
                          rho_kl_inner_core,beta_kl_inner_core, &
                          alpha_kl_inner_core, &
                          accel_inner_core,b_displ_inner_core, &
                          deltat,displ_inner_core,hprime_xx,hprime_xxT,&
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                          RECOMPUTE_STRAIN_DO_NOT_STORE,& 
                          epsilondev_inner_core,eps_trace_over_3_inner_core) 

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    rho_kl_inner_core, beta_kl_inner_core, alpha_kl_inner_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: &
     accel_inner_core,displ_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core

  real(kind=CUSTOM_REAL) deltat
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  logical :: RECOMPUTE_STRAIN_DO_NOT_STORE 
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev_inner_core 
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: eps_trace_over_3_inner_core 

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_matrix,b_epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_matrix,&
                                                          b_eps_trace_over_3_loc_matrix

  integer :: ispec,iglob

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! inner_core
  do ispec = 1, NSPEC_INNER_CORE

    if(COMPUTE_AND_STORE_STRAIN .and. RECOMPUTE_STRAIN_DO_NOT_STORE)then  
        eps_trace_over_3_loc_matrix(:,:,:) = eps_trace_over_3_inner_core(:,:,:,ispec)
        epsilondev_loc_matrix(1,:,:,:) = epsilondev_inner_core(1,:,:,:,ispec)
        epsilondev_loc_matrix(2,:,:,:) = epsilondev_inner_core(2,:,:,:,ispec)
        epsilondev_loc_matrix(3,:,:,:) = epsilondev_inner_core(3,:,:,:,ispec)
        epsilondev_loc_matrix(4,:,:,:) = epsilondev_inner_core(4,:,:,:,ispec)
        epsilondev_loc_matrix(5,:,:,:) = epsilondev_inner_core(5,:,:,:,ispec)
    else
      call compute_element_strain_undo_att_Dev(ispec,NGLOB_inner_core,NSPEC_inner_core,&
                                             displ_inner_core,ibool_inner_core,hprime_xx,hprime_xxT,&
                                             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                                             epsilondev_loc_matrix,eps_trace_over_3_loc_matrix)
    endif

    call compute_element_strain_undo_att_Dev(ispec,NGLOB_inner_core,NSPEC_inner_core,&
                                             b_displ_inner_core,ibool_inner_core,hprime_xx,hprime_xxT,&
                                             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                                             b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)

#ifdef FORCE_VECTORIZATION
        do ijk = 1,NGLLCUBE
          iglob = ibool_inner_core(ijk,1,1,ispec)

          rho_kl_inner_core(ijk,1,1,ispec) =  rho_kl_inner_core(ijk,1,1,ispec) &
             + deltat * (accel_inner_core(1,iglob) * b_displ_inner_core(1,iglob) &
             + accel_inner_core(2,iglob) * b_displ_inner_core(2,iglob) &
             + accel_inner_core(3,iglob) * b_displ_inner_core(3,iglob))

          beta_kl_inner_core(ijk,1,1,ispec) =  beta_kl_inner_core(ijk,1,1,ispec) &
             + deltat * (epsilondev_loc_matrix(1,ijk,1,1)*b_epsilondev_loc_matrix(1,ijk,1,1) &
                + epsilondev_loc_matrix(2,ijk,1,1)*b_epsilondev_loc_matrix(2,ijk,1,1) &
                + (epsilondev_loc_matrix(1,ijk,1,1)+epsilondev_loc_matrix(2,ijk,1,1)) * &
                  (b_epsilondev_loc_matrix(1,ijk,1,1)+b_epsilondev_loc_matrix(2,ijk,1,1)) &
                + 2.d0 * (epsilondev_loc_matrix(3,ijk,1,1)*b_epsilondev_loc_matrix(3,ijk,1,1) &
                + epsilondev_loc_matrix(4,ijk,1,1)*b_epsilondev_loc_matrix(4,ijk,1,1) &
                + epsilondev_loc_matrix(5,ijk,1,1)*b_epsilondev_loc_matrix(5,ijk,1,1)))

          alpha_kl_inner_core(ijk,1,1,ispec) = alpha_kl_inner_core(ijk,1,1,ispec) &
                + deltat * (9.d0 * eps_trace_over_3_loc_matrix(ijk,1,1) * b_eps_trace_over_3_loc_matrix(ijk,1,1))
        enddo
#else
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_inner_core(i,j,k,ispec)

          rho_kl_inner_core(i,j,k,ispec) =  rho_kl_inner_core(i,j,k,ispec) &
             + deltat * (accel_inner_core(1,iglob) * b_displ_inner_core(1,iglob) &
             + accel_inner_core(2,iglob) * b_displ_inner_core(2,iglob) &
             + accel_inner_core(3,iglob) * b_displ_inner_core(3,iglob))

          beta_kl_inner_core(i,j,k,ispec) =  beta_kl_inner_core(i,j,k,ispec) &
             + deltat * (epsilondev_loc_matrix(1,i,j,k)*b_epsilondev_loc_matrix(1,i,j,k) &
                + epsilondev_loc_matrix(2,i,j,k)*b_epsilondev_loc_matrix(2,i,j,k) &
                + (epsilondev_loc_matrix(1,i,j,k)+epsilondev_loc_matrix(2,i,j,k)) * &
                  (b_epsilondev_loc_matrix(1,i,j,k)+b_epsilondev_loc_matrix(2,i,j,k)) &
                + 2.d0 * (epsilondev_loc_matrix(3,i,j,k)*b_epsilondev_loc_matrix(3,i,j,k) &
                + epsilondev_loc_matrix(4,i,j,k)*b_epsilondev_loc_matrix(4,i,j,k) &
                + epsilondev_loc_matrix(5,i,j,k)*b_epsilondev_loc_matrix(5,i,j,k)))

          alpha_kl_inner_core(i,j,k,ispec) = alpha_kl_inner_core(i,j,k,ispec) &
                + deltat * (9.d0 * eps_trace_over_3_loc_matrix(i,j,k) * b_eps_trace_over_3_loc_matrix(i,j,k))
        enddo
      enddo
    enddo
#endif

  enddo

  end subroutine compute_kernels_inner_core


!
!-------------------------------------------------------------------------------------------------
!
! Subroutines to compute the kernels for the 21 elastic coefficients

  subroutine compute_strain_product(prod,eps_trace_over_3,epsdev,&
                                    b_eps_trace_over_3,b_epsdev)

  ! Purpose: compute the 21 strain products at a grid point
  ! (ispec,i,j,k fixed) and at a time t to compute then the kernels cij_kl (Voigt notation)
  ! (eq. 15 of Tromp et al., 2005)
  ! prod(1)=eps11*eps11 -> c11, prod(2)=eps11eps22 -> c12, prod(3)=eps11eps33 -> c13, ...
  ! prod(7)=eps22*eps22 -> c22, prod(8)=eps22eps33 -> c23, prod(9)=eps22eps23 -> c24, ...
  ! prod(19)=eps13*eps13 -> c55, prod(20)=eps13eps12 -> c56, prod(21)=eps12eps12 -> c66
  ! This then gives how the 21 kernels are organized

  implicit none
  include  "constants.h"

  real(kind=CUSTOM_REAL),dimension(21) :: prod
  real(kind=CUSTOM_REAL) :: eps_trace_over_3,b_eps_trace_over_3
  real(kind=CUSTOM_REAL),dimension(5) :: epsdev,b_epsdev

  real(kind=CUSTOM_REAL), dimension(6) :: eps,b_eps
  integer :: p,i,j

  ! Building of the local matrix of the strain tensor
  ! for the adjoint field and the regular backward field
  eps(1:2)=epsdev(1:2)+eps_trace_over_3           ! eps11 and eps22
  eps(3)=-(eps(1)+eps(2))+3*eps_trace_over_3      ! eps33
  eps(4)=epsdev(5)                                ! eps23
  eps(5)=epsdev(4)                                ! eps13
  eps(6)=epsdev(3)                                ! eps12

  b_eps(1:2)=b_epsdev(1:2)+b_eps_trace_over_3
  b_eps(3)=-(b_eps(1)+b_eps(2))+3*b_eps_trace_over_3
  b_eps(4)=b_epsdev(5)
  b_eps(5)=b_epsdev(4)
  b_eps(6)=b_epsdev(3)

  ! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
  p=1
  do i=1,6
    do j=i,6
      prod(p)=eps(i)*b_eps(j)
      if(j>i) then
        prod(p)=prod(p)+eps(j)*b_eps(i)
        if(j>3 .and. i<4) prod(p) = prod(p) * 2.0_CUSTOM_REAL
      endif
      if(i>3) prod(p) = prod(p) * 4.0_CUSTOM_REAL
      p=p+1
    enddo
  enddo

  end subroutine compute_strain_product

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_hessian(ibool_crust_mantle, &
                                    hess_kl_crust_mantle, &
                                    accel_crust_mantle,b_accel_crust_mantle, &
                                    deltat)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    hess_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
      b_accel_crust_mantle

  real(kind=CUSTOM_REAL) deltat

  ! local parameters
  integer :: ispec,iglob

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

! approximates Hessian term with adjoint acceleration and backward/reconstructed acceleration

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE

#ifdef FORCE_VECTORIZATION
        do ijk = 1,NGLLCUBE
          iglob = ibool_crust_mantle(ijk,1,1,ispec)
          hess_kl_crust_mantle(ijk,1,1,ispec) =  hess_kl_crust_mantle(ijk,1,1,ispec) &
             + deltat * (accel_crust_mantle(1,iglob) * b_accel_crust_mantle(1,iglob) &
             + accel_crust_mantle(2,iglob) * b_accel_crust_mantle(2,iglob) &
             + accel_crust_mantle(3,iglob) * b_accel_crust_mantle(3,iglob))
        enddo
#else
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          hess_kl_crust_mantle(i,j,k,ispec) =  hess_kl_crust_mantle(i,j,k,ispec) &
             + deltat * (accel_crust_mantle(1,iglob) * b_accel_crust_mantle(1,iglob) &
             + accel_crust_mantle(2,iglob) * b_accel_crust_mantle(2,iglob) &
             + accel_crust_mantle(3,iglob) * b_accel_crust_mantle(3,iglob))
        enddo
      enddo
    enddo
#endif

  enddo

  end subroutine compute_kernels_hessian
