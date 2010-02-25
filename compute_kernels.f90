!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology / Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology / Princeton University and University of Pau / CNRS / INRIA
!                            March 2010
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
                          epsilondev_crust_mantle,b_epsilondev_crust_mantle, &
                          eps_trace_over_3_crust_mantle,b_eps_trace_over_3_crust_mantle, &
                          deltat)
  
  implicit none
  
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    rho_kl_crust_mantle, beta_kl_crust_mantle, alpha_kl_crust_mantle
    
  real(kind=CUSTOM_REAL), dimension(21,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    cijkl_kl_crust_mantle

  
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT) :: &
    epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    b_epsilondev_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY) :: &
    eps_trace_over_3_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    b_eps_trace_over_3_crust_mantle

  real(kind=CUSTOM_REAL) deltat

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(21) :: prod !, cijkl_kl_local  
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc
  integer :: i,j,k,ispec,iglob

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE
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

          epsilondev_loc(:) = epsilondev_crust_mantle(:,i,j,k,ispec)
          b_epsilondev_loc(:) = b_epsilondev_crust_mantle(:,i,j,k,ispec)

          ! For anisotropic kernels
          if (ANISOTROPIC_KL) then

            call compute_strain_product(prod,eps_trace_over_3_crust_mantle(i,j,k,ispec),epsilondev_loc, &
                                        b_eps_trace_over_3_crust_mantle(i,j,k,ispec),b_epsilondev_loc)
            cijkl_kl_crust_mantle(:,i,j,k,ispec) = cijkl_kl_crust_mantle(:,i,j,k,ispec) + deltat * prod(:)

          else
            
            ! kernel for shear modulus, see e.g. Tromp et al. (2005), equation (17)
            ! note: multiplication with 2*mu(x) will be done after the time loop
            beta_kl_crust_mantle(i,j,k,ispec) =  beta_kl_crust_mantle(i,j,k,ispec) &
               + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
               + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
               + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) + &
                epsilondev_loc(5)*b_epsilondev_loc(5)) )


            ! kernel for bulk modulus, see e.g. Tromp et al. (2005), equation (18)
            ! note: multiplication with kappa(x) will be done after the time loop
            alpha_kl_crust_mantle(i,j,k,ispec) = alpha_kl_crust_mantle(i,j,k,ispec) &
               + deltat * (9 * eps_trace_over_3_crust_mantle(i,j,k,ispec) &
                             * b_eps_trace_over_3_crust_mantle(i,j,k,ispec))

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
                        div_displ_outer_core,b_div_displ_outer_core, &
                        rhostore_outer_core,kappavstore_outer_core, &
                        rho_kl_outer_core,alpha_kl_outer_core, &
                        deviatoric_outercore,nspec_beta_kl_outer_core,beta_kl_outer_core, &
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

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_OUTER_CORE) :: vector_accel_outer_core,&
             vector_displ_outer_core, b_vector_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: div_displ_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: b_div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        rhostore_outer_core,kappavstore_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: &
    rho_kl_outer_core,alpha_kl_outer_core

  integer nspec_beta_kl_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_beta_kl_outer_core) :: &
    beta_kl_outer_core
  logical deviatoric_outercore
  
  real(kind=CUSTOM_REAL) deltat
  
  ! local parameters
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,kappal
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) :: tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc

  integer :: i,j,k,l,ispec,iglob

  ! outer_core -- compute the actual displacement and acceleration (NDIM,NGLOBMAX_OUTER_CORE)
  do ispec = 1, NSPEC_OUTER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_outer_core(i,j,k,ispec)

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
            tempx3l = tempx3l +  b_displ_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
          enddo

          b_vector_displ_outer_core(1,iglob) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          b_vector_displ_outer_core(2,iglob) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          b_vector_displ_outer_core(3,iglob) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l


          !deviatoric kernel check
          if( deviatoric_outercore ) then

            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempx3l = 0._CUSTOM_REAL

            tempy1l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL

            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL

            ! assumes NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              tempx1l = tempx1l + b_vector_displ_outer_core(1,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              tempy1l = tempy1l + b_vector_displ_outer_core(2,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              tempz1l = tempz1l + b_vector_displ_outer_core(3,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)

              tempx2l = tempx2l + b_vector_displ_outer_core(1,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              tempy2l = tempy2l + b_vector_displ_outer_core(2,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              tempz2l = tempz2l + b_vector_displ_outer_core(3,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)

              tempx3l = tempx3l +  b_vector_displ_outer_core(1,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              tempy3l = tempy3l +  b_vector_displ_outer_core(2,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              tempz3l = tempz3l +  b_vector_displ_outer_core(3,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
            enddo


            !deviatoric strain
            b_epsilondev_loc(1) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l  &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )
                              
            b_epsilondev_loc(2) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l  &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

            b_epsilondev_loc(3) = 0.5*( xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l  &
                                      + xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l ) &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

            b_epsilondev_loc(4) = 0.5*( xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l  &
                                      + xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l ) &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

            b_epsilondev_loc(5) = 0.5*( xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l  &
                                      + xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l ) &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )
                          
          endif !deviatoric kernel check


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

          vector_accel_outer_core(1,iglob) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          vector_accel_outer_core(2,iglob) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          vector_accel_outer_core(3,iglob) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

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

          vector_displ_outer_core(1,iglob) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          vector_displ_outer_core(2,iglob) = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          vector_displ_outer_core(3,iglob) = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l


          !deviatoric kernel check
          if( deviatoric_outercore ) then

            tempx1l = 0._CUSTOM_REAL
            tempx2l = 0._CUSTOM_REAL
            tempx3l = 0._CUSTOM_REAL

            tempy1l = 0._CUSTOM_REAL
            tempy2l = 0._CUSTOM_REAL
            tempy3l = 0._CUSTOM_REAL

            tempz1l = 0._CUSTOM_REAL
            tempz2l = 0._CUSTOM_REAL
            tempz3l = 0._CUSTOM_REAL

            ! assumes NGLLX = NGLLY = NGLLZ
            do l=1,NGLLX
              tempx1l = tempx1l + vector_displ_outer_core(1,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              tempy1l = tempy1l + vector_displ_outer_core(2,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              tempz1l = tempz1l + vector_displ_outer_core(3,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)

              tempx2l = tempx2l + vector_displ_outer_core(1,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              tempy2l = tempy2l + vector_displ_outer_core(2,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              tempz2l = tempz2l + vector_displ_outer_core(3,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)

              tempx3l = tempx3l + vector_displ_outer_core(1,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              tempy3l = tempy3l + vector_displ_outer_core(2,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              tempz3l = tempz3l + vector_displ_outer_core(3,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
            enddo


            !deviatoric strain
            epsilondev_loc(1) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l  &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )
                              
            epsilondev_loc(2) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l  &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

            epsilondev_loc(3) = 0.5*( xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l  &
                                      + xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l ) &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

            epsilondev_loc(4) = 0.5*( xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l  &
                                      + xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l ) &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

            epsilondev_loc(5) = 0.5*( xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l  &
                                      + xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l ) &
                - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                              + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                              + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

          beta_kl_outer_core(i,j,k,ispec) =  beta_kl_outer_core(i,j,k,ispec) &
               + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
               + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
               + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) + &
                epsilondev_loc(5)*b_epsilondev_loc(5)) )
                          
          endif !deviatoric kernel check



          rho_kl_outer_core(i,j,k,ispec) = rho_kl_outer_core(i,j,k,ispec) &
             + deltat * dot_product(vector_accel_outer_core(:,iglob), b_vector_displ_outer_core(:,iglob))

          kappal = rhostore_outer_core(i,j,k,ispec)/kappavstore_outer_core(i,j,k,ispec)
          
          div_displ_outer_core(i,j,k,ispec) =  kappal * accel_outer_core(iglob)
          b_div_displ_outer_core(i,j,k,ispec) =  kappal * b_accel_outer_core(iglob)

          alpha_kl_outer_core(i,j,k,ispec) = alpha_kl_outer_core(i,j,k,ispec) &
             + deltat * div_displ_outer_core(i,j,k,ispec) * b_div_displ_outer_core(i,j,k,ispec)


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
                          epsilondev_inner_core,b_epsilondev_inner_core, &
                          eps_trace_over_3_inner_core,b_eps_trace_over_3_inner_core, &
                          deltat)
  

  implicit none
  
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    rho_kl_inner_core, beta_kl_inner_core, alpha_kl_inner_core
    
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: &
     accel_inner_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE_ADJOINT) :: &
    b_displ_inner_core

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STR_OR_ATT) :: &
    epsilondev_inner_core
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    b_epsilondev_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_STRAIN_ONLY) :: &
    eps_trace_over_3_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ADJOINT) :: &
    b_eps_trace_over_3_inner_core

  real(kind=CUSTOM_REAL) deltat

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
  
  integer :: i,j,k,ispec,iglob
  

  ! inner_core
  do ispec = 1, NSPEC_INNER_CORE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool_inner_core(i,j,k,ispec)

          rho_kl_inner_core(i,j,k,ispec) =  rho_kl_inner_core(i,j,k,ispec) &
             + deltat * (accel_inner_core(1,iglob) * b_displ_inner_core(1,iglob) &
             + accel_inner_core(2,iglob) * b_displ_inner_core(2,iglob) &
             + accel_inner_core(3,iglob) * b_displ_inner_core(3,iglob) )

          epsilondev_loc(:) = epsilondev_inner_core(:,i,j,k,ispec)
          b_epsilondev_loc(:) = b_epsilondev_inner_core(:,i,j,k,ispec)
          beta_kl_inner_core(i,j,k,ispec) =  beta_kl_inner_core(i,j,k,ispec) &
             + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
                + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
                + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) &
                + epsilondev_loc(5)*b_epsilondev_loc(5)) )

          alpha_kl_inner_core(i,j,k,ispec) = alpha_kl_inner_core(i,j,k,ispec) &
             + deltat * (9 * eps_trace_over_3_inner_core(i,j,k,ispec) * b_eps_trace_over_3_inner_core(i,j,k,ispec))
        enddo
      enddo
    enddo
  enddo

  end subroutine compute_kernels_inner_core
  

!
!-------------------------------------------------------------------------------------------------
!
! Subroutines to compute the kernels for the 21 elastic coefficients
! Last modified 19/04/2007

!-------------------------------------------------------------------
  subroutine compute_strain_product(prod,eps_trace_over_3,epsdev,&
                          b_eps_trace_over_3,b_epsdev)

  ! Purpose : compute the 21 strain products at a grid point
  ! (ispec,i,j,k fixed) and at a time t to compute then the kernels cij_kl (Voigt notation)
  ! (eq. 15 of Tromp et al., 2005)
  ! prod(1)=eps11*eps11 -> c11, prod(2)=eps11eps22 -> c12, prod(3)=eps11eps33 -> c13, ...
  ! prod(7)=eps22*eps22 -> c22, prod(8)=eps22eps33 -> c23, prod(9)=eps22eps23 -> c24, ...
  ! prod(19)=eps13*eps13 -> c55, prod(20)=eps13eps12 -> c56, prod(21)=eps12eps12 -> c66
  ! This then gives how the 21 kernels are organized
  ! For crust_mantle

  ! Modif 09/11/2005

  implicit none
  include  "constants.h"

  real(kind=CUSTOM_REAL),dimension(21) :: prod
  real(kind=CUSTOM_REAL) :: eps_trace_over_3,b_eps_trace_over_3
  real(kind=CUSTOM_REAL),dimension(5) :: epsdev,b_epsdev

  real(kind=CUSTOM_REAL), dimension(6) :: eps,b_eps
  integer :: p,i,j

  ! Building of the local matrix of the strain tensor
  ! for the adjoint field and the regular backward field
  eps(1:2)=epsdev(1:2)+eps_trace_over_3           !eps11 et eps22
  eps(3)=-(eps(1)+eps(2))+3*eps_trace_over_3     !eps33
  eps(4)=epsdev(5)                                !eps23
  eps(5)=epsdev(4)                                !eps13
  eps(6)=epsdev(3)                                !eps12

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
            if(j>3 .and. i<4) prod(p)=prod(p)*2
       endif
       if(i>3) prod(p)=prod(p)*4
       p=p+1
       enddo
  enddo

  end subroutine compute_strain_product

!
!-------------------------------------------------------------------------------------------------
!

  subroutine rotate_kernels_dble(cij_kl,cij_kll,theta_in,phi_in)

! Purpose : compute the kernels in r,theta,phi (cij_kll)
! from the kernels in x,y,z (cij_kl) (x,y,z <-> r,theta,phi)
! At r,theta,phi fixed
! theta and phi are in radians

! Coeff from Min's routine rotate_anisotropic_tensor
! with the help of Collect[Expand[cij],{dij}] in Mathematica

! Definition of the output array cij_kll :
! cij_kll(1) = C11 ; cij_kll(2) = C12 ; cij_kll(3) = C13
! cij_kll(4) = C14 ; cij_kll(5) = C15 ; cij_kll(6) = C16
! cij_kll(7) = C22 ; cij_kll(8) = C23 ; cij_kll(9) = C24
! cij_kll(10) = C25 ; cij_kll(11) = C26 ; cij_kll(12) = C33
! cij_kll(13) = C34 ; cij_kll(14) = C35 ; cij_kll(15) = C36
! cij_kll(16) = C44 ; cij_kll(17) = C45 ; cij_kll(18) = C46
! cij_kll(19) = C55 ; cij_kll(20) = C56 ; cij_kll(21) = C66
! where the Cij (Voigt's notation) are defined as function of
! the components of the elastic tensor in spherical coordinates
! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

  implicit none
  include  "constants.h"

  real(kind=CUSTOM_REAL) :: theta_in,phi_in
  real(kind=CUSTOM_REAL),dimension(21) :: cij_kll,cij_kl

  double precision :: theta,phi
  double precision :: costheta,sintheta,cosphi,sinphi
  double precision :: costhetasq,sinthetasq,cosphisq,sinphisq
  double precision :: costwotheta,sintwotheta,costwophi,sintwophi
  double precision :: cosfourtheta,sinfourtheta,cosfourphi,sinfourphi
  double precision :: costhetafour,sinthetafour,cosphifour,sinphifour
  double precision :: sintwophisq,sintwothetasq
  double precision :: costhreetheta,sinthreetheta,costhreephi,sinthreephi


   if (CUSTOM_REAL == SIZE_REAL) then
      theta = dble(theta_in)
      phi = dble(phi_in)
    else
      theta = theta_in
      phi = phi_in
    endif

  costheta = dcos(theta)
  sintheta = dsin(theta)
  cosphi = dcos(phi)
  sinphi = dsin(phi)

  costhetasq = costheta * costheta
  sinthetasq = sintheta * sintheta
  cosphisq = cosphi * cosphi
  sinphisq = sinphi * sinphi

  costhetafour = costhetasq * costhetasq
  sinthetafour = sinthetasq * sinthetasq
  cosphifour = cosphisq * cosphisq
  sinphifour = sinphisq * sinphisq

  costwotheta = dcos(2.d0*theta)
  sintwotheta = dsin(2.d0*theta)
  costwophi = dcos(2.d0*phi)
  sintwophi = dsin(2.d0*phi)

  costhreetheta=dcos(3.d0*theta)
  sinthreetheta=dsin(3.d0*theta)
  costhreephi=dcos(3.d0*phi)
  sinthreephi=dsin(3.d0*phi)

  cosfourtheta = dcos(4.d0*theta)
  sinfourtheta = dsin(4.d0*theta)
  cosfourphi = dcos(4.d0*phi)
  sinfourphi = dsin(4.d0*phi)
  sintwothetasq = sintwotheta * sintwotheta
  sintwophisq = sintwophi * sintwophi


 cij_kll(1) = 1.d0/16.d0* (cij_kl(16) - cij_kl(16)* costwophi + &
     16.d0* cosphi*cosphisq* costhetafour* (cij_kl(1)* cosphi + cij_kl(6)* sinphi) + &
     2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq - &
     2.d0* (cij_kl(16)* cosfourtheta* sinphisq + &
     2.d0* costhetafour* (-4* cij_kl(7)* sinphifour - &
     (cij_kl(2) + cij_kl(21))* sintwophisq) + &
     8.d0* cij_kl(5)* cosphi*cosphisq* costheta*costhetasq* sintheta - &
     8.d0* cij_kl(8)* costhetasq* sinphisq* sinthetasq - &
     8.d0* cij_kl(12)* sinthetafour + &
     8.d0* cosphisq* costhetasq* sintheta* ((cij_kl(4) + &
     cij_kl(20))* costheta* sinphi - &
     (cij_kl(3) + cij_kl(19))*sintheta) + &
     8.d0* cosphi* costheta* (-cij_kl(11)* costheta*costhetasq* &
     sinphi*sinphisq + (cij_kl(10) + cij_kl(18))* costhetasq* sinphisq* sintheta + &
     cij_kl(14)* sintheta*sinthetasq) + 2.d0* sinphi* (cij_kl(13) + &
     cij_kl(9)* sinphisq)* sintwotheta + &
     sinphi* (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta))

 cij_kll(2) = 1.d0/4.d0* (costhetasq* (cij_kl(1) + 3.d0* cij_kl(2) + cij_kl(7) - &
      cij_kl(21) + (-cij_kl(1) + cij_kl(2) - cij_kl(7) + &
      cij_kl(21))* cosfourphi + (-cij_kl(6) + cij_kl(11))* sinfourphi) + &
      4.d0* (cij_kl(8)* cosphisq - cij_kl(15)* cosphi* sinphi + &
      cij_kl(3)* sinphisq)* sinthetasq - &
      2.d0* (cij_kl(10)* cosphisq*cosphi + &
      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
      cij_kl(4)* sinphisq*sinphi)* sintwotheta)

 cij_kll(3) = 1.d0/8.d0* (sintwophi* (3.d0* cij_kl(15) - cij_kl(17) + &
     4.d0* (cij_kl(2) + cij_kl(21))* costhetasq* sintwophi* sinthetasq) + &
     4.d0* cij_kl(12)* sintwothetasq + 4.d0* cij_kl(1)* cosphifour* sintwothetasq + &
     2.d0* cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
     cij_kl(5)* sinfourtheta) + 2.d0* cosphisq* (3.d0* cij_kl(3) -  cij_kl(19) + &
     (cij_kl(3) + cij_kl(19))* cosfourtheta + &
     (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
     2.d0* sinphi* (sinphi* (3.d0* cij_kl(8) - &
     cij_kl(16) + (cij_kl(8) + cij_kl(16))* cosfourtheta + &
     2.d0* cij_kl(7)* sinphisq* sintwothetasq)+ &
     (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta)+ &
     2.d0* cosphi* ((cij_kl(15) + cij_kl(17))* cosfourtheta* sinphi + &
     8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
     (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)*sinfourtheta))

 cij_kll(4) = 1.d0/8.d0* (cosphi* costheta *(5.d0* cij_kl(4) - &
     cij_kl(9) + 4.d0* cij_kl(13) - &
     3.d0* cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
     4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
     1.d0/2.d0* (cij_kl(4) - cij_kl(9) + &
     cij_kl(20))* costhreephi * (costheta + 3.d0* costhreetheta) - &
     costheta* (-cij_kl(5) + 5.d0* cij_kl(10) + &
     4.d0* cij_kl(14) - 3.d0* cij_kl(18) + &
     (3.d0* cij_kl(5) + cij_kl(10) - &
     4.d0* cij_kl(14) + cij_kl(18))* costwotheta)* sinphi - &
     1.d0/2.d0* (cij_kl(5) - cij_kl(10) - cij_kl(18))* (costheta + &
     3.d0* costhreetheta)* sinthreephi + &
     4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* costhetasq* sintheta - &
     4.d0* (cij_kl(1) + cij_kl(3) - cij_kl(7) - cij_kl(8) + cij_kl(16) - cij_kl(19) + &
     (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + &
     cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi* sintheta - &
     4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
     cij_kl(21))* costhetasq* sinfourphi* sintheta + &
     costwophi* ((cij_kl(6) + cij_kl(11) + 6.d0* cij_kl(15) - &
     2.d0* cij_kl(17))* sintheta + &
     (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))

 cij_kll(5) = 1.d0/4.d0* (2.d0* (cij_kl(4) + &
     cij_kl(20))* cosphisq* (costwotheta + cosfourtheta)* sinphi + &
     2.d0* cij_kl(9)* (costwotheta + cosfourtheta)* sinphi*sinphisq + &
     16.d0* cij_kl(1)* cosphifour* costheta*costhetasq* sintheta + &
     4.d0* costheta*costhetasq* (-2.d0* cij_kl(8)* sinphisq + &
     4.d0* cij_kl(7)* sinphifour + &
     (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta + &
     4.d0* cij_kl(13)* (1.d0 + 2.d0* costwotheta)* sinphi* sinthetasq + &
     8.d0* costheta* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta*sinthetasq + &
     2.d0* cosphi*cosphisq* (cij_kl(5)* (costwotheta + cosfourtheta) + &
     8.d0* cij_kl(6)* costheta*costhetasq* sinphi* sintheta) + &
     2.d0* cosphi* (cosfourtheta* (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
     costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
     8.d0* cij_kl(11)* costheta*costhetasq* sinphi*sinphisq* sintheta) - &
     (cij_kl(3) + cij_kl(16) + cij_kl(19) + &
     (cij_kl(3) - cij_kl(16) + cij_kl(19))* costwophi + &
     (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)

 cij_kll(6) = 1.d0/2.d0* costheta*costhetasq* ((cij_kl(6) + cij_kl(11))* costwophi + &
      (cij_kl(6) - cij_kl(11))* cosfourphi + 2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi) + &
      1.d0/4.d0* costhetasq* (-(cij_kl(4) + 3* cij_kl(9) + cij_kl(20))* cosphi - &
      3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
      3.d0* (cij_kl(5) - cij_kl(10) - cij_kl(18))* sinthreephi)* sintheta + &
      costheta* ((cij_kl(15) + cij_kl(17))* costwophi + &
      (-cij_kl(3) + cij_kl(8) + cij_kl(16) - cij_kl(19))* sintwophi)* sinthetasq + &
      (-cij_kl(13)* cosphi + cij_kl(14)* sinphi)* sintheta*sinthetasq

 cij_kll(7) = cij_kl(7)* cosphifour - cij_kl(11)* cosphi*cosphisq* sinphi + &
      (cij_kl(2) + cij_kl(21))* cosphisq* sinphisq - &
      cij_kl(6)* cosphi* sinphi*sinphisq + &
      cij_kl(1)* sinphifour

 cij_kll(8) = 1.d0/2.d0* (2.d0* costhetasq* sinphi* (-cij_kl(15)* cosphi + &
      cij_kl(3)* sinphi) + 2.d0* cij_kl(2)* cosphifour* sinthetasq + &
      (2.d0* cij_kl(2)* sinphifour + &
      (cij_kl(1) + cij_kl(7) - cij_kl(21))* sintwophisq)* sinthetasq + &
      cij_kl(4)* sinphi*sinphisq* sintwotheta + &
      cosphi*cosphisq* (2.d0* (-cij_kl(6) + cij_kl(11))* sinphi* sinthetasq + &
      cij_kl(10)* sintwotheta) + cosphi* sinphisq* (2.d0* (cij_kl(6) - &
      cij_kl(11))* sinphi* sinthetasq + &
      (cij_kl(5) - cij_kl(18))* sintwotheta) + &
      cosphisq* (2.d0* cij_kl(8)* costhetasq + &
      (cij_kl(9) - cij_kl(20))* sinphi* sintwotheta))

 cij_kll(9) = cij_kl(11)* cosphifour* sintheta - sinphi*sinphisq* (cij_kl(5)* costheta + &
      cij_kl(6)* sinphi* sintheta) +  cosphisq* sinphi* (-(cij_kl(10) + &
      cij_kl(18))* costheta + &
      3.d0* (cij_kl(6) - cij_kl(11))* sinphi* sintheta) + &
      cosphi* sinphisq* ((cij_kl(4) + cij_kl(20))* costheta + &
      2.d0* (-2.d0* cij_kl(1) + cij_kl(2) + cij_kl(21))* sinphi* sintheta) + &
      cosphi*cosphisq* (cij_kl(9)* costheta - 2.d0* (cij_kl(2) - 2.d0* cij_kl(7) + &
      cij_kl(21))* sinphi* sintheta)

 cij_kll(10) = 1.d0/4.d0* (4.d0* costwotheta* (cij_kl(10)* cosphi*cosphisq + &
      (cij_kl(9) - cij_kl(20))* cosphisq* sinphi + &
      (cij_kl(5) - cij_kl(18))* cosphi* sinphisq + &
      cij_kl(4)* sinphi*sinphisq) + (cij_kl(1) + 3.d0* cij_kl(2) - &
      2.d0* cij_kl(3) + cij_kl(7) - &
      2.d0* cij_kl(8) - cij_kl(21) + 2.d0* (cij_kl(3) - cij_kl(8))* costwophi + &
      (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
      2.d0* cij_kl(15)* sintwophi + &
      (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)

 cij_kll(11) = 1.d0/4.d0* (2.d0* costheta* ((cij_kl(6) + cij_kl(11))* costwophi + &
      (-cij_kl(6) + cij_kl(11))* cosfourphi + &
      2.d0* (-cij_kl(1) + cij_kl(7))* sintwophi + &
      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(21))* sinfourphi) + &
      (-(cij_kl(4) + 3.d0* cij_kl(9) + cij_kl(20))* cosphi + &
      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi + &
      (3.d0* cij_kl(5) + cij_kl(10) + cij_kl(18))* sinphi + &
      (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sintheta)

 cij_kll(12) = 1.d0/16.d0* (cij_kl(16) - 2.d0* cij_kl(16)* cosfourtheta* sinphisq + &
      costwophi* (-cij_kl(16) + 8.d0* costheta* sinthetasq* ((cij_kl(3) - &
      cij_kl(8) + cij_kl(19))* costheta + &
      (cij_kl(5) - cij_kl(10) - cij_kl(18))* cosphi* sintheta)) + &
      2.d0* (cij_kl(15) + cij_kl(17))* sintwophi* sintwothetasq + &
      2.d0* (8.d0* cij_kl(12)* costhetafour + &
      8.d0* cij_kl(14)* cosphi* costheta*costhetasq* sintheta + &
      4.d0* cosphi* costheta* (cij_kl(5) + cij_kl(10) + cij_kl(18) + &
      (cij_kl(4) + cij_kl(20))* sintwophi)* &
      sintheta*sinthetasq + 8.d0* cij_kl(1)* cosphifour* sinthetafour + &
      8.d0* cij_kl(6)* cosphi*cosphisq* sinphi* sinthetafour + &
      8.d0* cij_kl(11)* cosphi* sinphi*sinphisq* sinthetafour + &
      8.d0* cij_kl(7)* sinphifour* sinthetafour + &
      2.d0* cij_kl(2)* sintwophisq* sinthetafour + &
      2.d0* cij_kl(21)* sintwophisq* sinthetafour + &
      2.d0* cij_kl(13)* sinphi* sintwotheta + &
      2.d0* cij_kl(9)* sinphi*sinphisq* sintwotheta + &
      cij_kl(3)* sintwothetasq + cij_kl(8)* sintwothetasq + &
      cij_kl(19)* sintwothetasq + cij_kl(13)* sinphi* sinfourtheta - &
      cij_kl(9)* sinphi*sinphisq* sinfourtheta))

 cij_kll(13) = 1.d0/8.d0* (cosphi* costheta* (cij_kl(4) + 3.d0* cij_kl(9) + &
      4.d0* cij_kl(13) + cij_kl(20) - (cij_kl(4) + 3.d0* cij_kl(9) - &
      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + 4.d0* (-cij_kl(1) - &
      cij_kl(3) + cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19) + &
      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
      cij_kl(19))* costwotheta)* sintwophi* sintheta + &
      4.d0* (cij_kl(6) - cij_kl(11))* cosfourphi* sinthetasq*sintheta - &
      4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* sinfourphi* sinthetasq*sintheta + &
      costheta* ((-3.d0* cij_kl(5) - cij_kl(10) - 4.d0* cij_kl(14) - &
      cij_kl(18) + (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + &
      cij_kl(18))* costwotheta)* sinphi + 6.d0* ((cij_kl(4) - cij_kl(9) + &
      cij_kl(20))* costhreephi + (-cij_kl(5) + cij_kl(10) + &
      cij_kl(18))* sinthreephi)* sinthetasq) + costwophi* ((3* cij_kl(6) + &
      3.d0* cij_kl(11) + 2.d0* (cij_kl(15) + cij_kl(17)))* sintheta - &
      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
      cij_kl(17)))* sinthreetheta))

 cij_kll(14) = 1.d0/4.d0* (2.d0* cij_kl(13)* (costwotheta + cosfourtheta)* sinphi + &
      8.d0* costheta*costhetasq* (-2.d0* cij_kl(12) + cij_kl(8)* sinphisq)* sintheta + &
      4.d0* (cij_kl(4) + cij_kl(20))* cosphisq* (1.d0 + &
      2.d0* costwotheta)* sinphi* sinthetasq + &
      4.d0* cij_kl(9)* (1.d0 + 2.d0* costwotheta)* sinphi*sinphisq* sinthetasq + &
      16.d0* cij_kl(1)* cosphifour* costheta* sintheta*sinthetasq + &
      4.d0* costheta* (-2.d0* cij_kl(8)* sinphisq + 4.d0* cij_kl(7)* sinphifour + &
      (cij_kl(2) + cij_kl(21))* sintwophisq)* sintheta*sinthetasq + &
      4.d0* cosphi*cosphisq* sinthetasq* (cij_kl(5) + 2.d0* cij_kl(5)* costwotheta + &
      4.d0* cij_kl(6)* costheta* sinphi* sintheta) + &
      2.d0* cosphi* (cosfourtheta* (cij_kl(14) - (cij_kl(10) + cij_kl(18))* sinphisq) + &
      costwotheta* (cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq) + &
      8.d0* cij_kl(11)* costheta* sinphi*sinphisq* sintheta*sinthetasq) + &
      (cij_kl(3) + cij_kl(16) + cij_kl(19) + (cij_kl(3) - cij_kl(16) + &
      cij_kl(19))* costwophi + (cij_kl(15) + cij_kl(17))* sintwophi)* sinfourtheta)

 cij_kll(15) = costwophi* costheta* (-cij_kl(17) + (cij_kl(15) + cij_kl(17))* costhetasq) + &
       1.d0/16.d0* (-((11.d0* cij_kl(4) + cij_kl(9) + 4.d0* cij_kl(13) - &
       5.d0* cij_kl(20))* cosphi + (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (cij_kl(5) + 11.d0* cij_kl(10) + 4.d0* cij_kl(14) - &
       5.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
       cij_kl(18))* sinthreephi)* sintheta + &
       8.d0* costheta* ((-cij_kl(1) - cij_kl(3) + cij_kl(7) + cij_kl(8) - cij_kl(16) +&
       cij_kl(19) + (cij_kl(1) - cij_kl(3) - &
       cij_kl(7) + cij_kl(8) + cij_kl(16) - cij_kl(19))* costwotheta)* sintwophi +&
       ((cij_kl(6) + cij_kl(11))* costwophi + &
       (cij_kl(6) - cij_kl(11))* cosfourphi + (-cij_kl(1) + cij_kl(2) - cij_kl(7) +&
       cij_kl(21))* sinfourphi)* sinthetasq) +&
       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)

 cij_kll(16) = 1.d0/4.d0*(cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
       cij_kl(19) + cij_kl(21) + 2.d0*(cij_kl(16) - cij_kl(19))*costwophi* costhetasq + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(16) + &
       cij_kl(19) - cij_kl(21))*costwotheta - 2.d0* cij_kl(17)* costhetasq* sintwophi + &
       2.d0* ((-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sinthetasq + ((cij_kl(5) - cij_kl(10) +&
       cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) + cij_kl(18))* costhreephi +&
       (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - &
       (cij_kl(4) - cij_kl(9) + cij_kl(20))* sinthreephi)* sintwotheta)

 cij_kll(17) = 1.d0/8.d0* (4.d0* costwophi* costheta* (cij_kl(6) + cij_kl(11) - &
       2.d0* cij_kl(15) - (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + &
       cij_kl(17)))* costwotheta) - (2.d0* cosphi* (-3.d0* cij_kl(4) +&
       cij_kl(9) + 2.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) - cij_kl(9) + &
       cij_kl(20))* costwophi) - (cij_kl(5) - 5.d0* cij_kl(10) + &
       4.d0* cij_kl(14) + 3.d0* cij_kl(18))* sinphi + (-cij_kl(5) + cij_kl(10) + &
       cij_kl(18))* sinthreephi)* sintheta + &
       8.d0* costheta* ((-cij_kl(1) + cij_kl(3) + cij_kl(7) - cij_kl(8) + &
       (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
       cij_kl(19))* costwotheta)* sintwophi + ((cij_kl(6) - cij_kl(11))* cosfourphi + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* sinfourphi)* sinthetasq) +&
       ((cij_kl(4) + 3.d0* cij_kl(9) - 4.d0* cij_kl(13) + cij_kl(20))* cosphi + &
       3.d0* (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi - &
       (3.d0* cij_kl(5) + cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))* sinphi + &
       3.d0* (-cij_kl(5) + cij_kl(10) + cij_kl(18))* sinthreephi)* sinthreetheta)

 cij_kll(18) = 1.d0/2.d0* ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi* costwotheta - &
       (cij_kl(5) - cij_kl(10) - cij_kl(18))* costhreephi* costwotheta - &
       2.d0* (cij_kl(4) - cij_kl(9) + &
       (cij_kl(4) - cij_kl(9) + cij_kl(20))* costwophi)* costwotheta* sinphi + &
       (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + cij_kl(21) + &
       (-cij_kl(16) + cij_kl(19))* costwophi + &
       (-cij_kl(1) + cij_kl(2) - cij_kl(7) + cij_kl(21))* cosfourphi + &
       cij_kl(17)* sintwophi + &
       (-cij_kl(6) + cij_kl(11))* sinfourphi)* sintwotheta)

 cij_kll(19) = 1.d0/4.d0* (cij_kl(16) - cij_kl(16)* costwophi + &
      (-cij_kl(15) + cij_kl(17))* sintwophi + &
      4.d0* cij_kl(12)* sintwothetasq + &
      2.d0* (2.d0* cij_kl(1)* cosphifour* sintwothetasq + &
      cosphi*cosphisq* (8.d0* cij_kl(6)* costhetasq* sinphi* sinthetasq + &
      cij_kl(5)* sinfourtheta) + cosphisq* (-cij_kl(3) + cij_kl(19) + (cij_kl(3) +&
      cij_kl(19))* cosfourtheta + (cij_kl(4) + cij_kl(20))* sinphi* sinfourtheta) + &
      sinphi* (cosfourtheta* ((cij_kl(15) + cij_kl(17))* cosphi + &
      cij_kl(16)* sinphi) + (cij_kl(2) + cij_kl(7) - 2.d0* cij_kl(8) + cij_kl(21) + &
      (cij_kl(2) - cij_kl(7) + cij_kl(21))* costwophi)* sinphi* sintwothetasq + &
      (-cij_kl(13) + cij_kl(9)* sinphisq)* sinfourtheta) + &
      cosphi* (8.d0* cij_kl(11)* costhetasq* sinphi*sinphisq* sinthetasq + &
      (-cij_kl(14) + (cij_kl(10) + cij_kl(18))* sinphisq)* sinfourtheta)))

 cij_kll(20) = 1.d0/8.d0* (2.d0* cosphi* costheta* (-3.d0* cij_kl(4) - cij_kl(9) + &
      4.d0* cij_kl(13) + cij_kl(20) + (cij_kl(4) + 3.d0* cij_kl(9) - &
      4.d0* cij_kl(13) + cij_kl(20))* costwotheta) + &
      (cij_kl(4) - cij_kl(9) + cij_kl(20))* costhreephi* (costheta + &
      3.d0* costhreetheta) - &
      2.d0* costheta* (-cij_kl(5) - 3.d0* cij_kl(10) + 4.d0* cij_kl(14) + &
      cij_kl(18) + (3.d0* cij_kl(5) + &
      cij_kl(10) - 4.d0* cij_kl(14) + cij_kl(18))*costwotheta)* sinphi - &
      (cij_kl(5) - cij_kl(10) - cij_kl(18))* &
      (costheta + 3.d0* costhreetheta)* sinthreephi + 8.d0* (cij_kl(6) - &
      cij_kl(11))* cosfourphi* costhetasq* sintheta - 8.d0* (cij_kl(1) - &
      cij_kl(3) - cij_kl(7) + cij_kl(8) + &
      (cij_kl(1) - cij_kl(3) - cij_kl(7) + cij_kl(8) + cij_kl(16) - &
      cij_kl(19))* costwotheta)* sintwophi* sintheta - &
      8.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* costhetasq* sinfourphi* sintheta + &
      2.d0* costwophi* ((cij_kl(6) + cij_kl(11) - 2.d0* cij_kl(15) + &
      2.d0* cij_kl(17))* sintheta + &
      (cij_kl(6) + cij_kl(11) - 2.d0* (cij_kl(15) + cij_kl(17)))* sinthreetheta))

 cij_kll(21) = 1.d0/4.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) + cij_kl(16) + &
      cij_kl(19) + cij_kl(21) - 2.d0* (cij_kl(1) - cij_kl(2) + cij_kl(7) - &
      cij_kl(21))* cosfourphi* costhetasq + &
      (cij_kl(1) - cij_kl(2) + cij_kl(7) - cij_kl(16) - cij_kl(19) + &
      cij_kl(21))* costwotheta + &
      2.d0* (-cij_kl(6) + cij_kl(11))* costhetasq* sinfourphi - &
      2.d0* ((-cij_kl(16) + cij_kl(19))* costwophi + cij_kl(17)* sintwophi)* sinthetasq - &
      ((cij_kl(5) - cij_kl(10) + cij_kl(18))* cosphi + (-cij_kl(5) + cij_kl(10) +&
      cij_kl(18))* costhreephi + &
      (-cij_kl(4) + cij_kl(9) + cij_kl(20))* sinphi - (cij_kl(4) - cij_kl(9) + &
      cij_kl(20))* sinthreephi)* sintwotheta)

  end subroutine rotate_kernels_dble

!-----------------------------------------------------------------------------



