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
                          deltat,hprime_xx,hprime_xxT,&
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,ANISOTROPIC_KL,&
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
     accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_displ_crust_mantle

  real(kind=CUSTOM_REAL) deltat
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  logical :: ANISOTROPIC_KL

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: epsilondev_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: eps_trace_over_3_crust_mantle

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: b_epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: b_eps_trace_over_3_loc_matrix
  integer :: ispec,iglob
  real(kind=CUSTOM_REAL) :: eps1, eps2, eps3, eps4, eps5, eps6, b_eps1, b_eps2, b_eps3, b_eps4, b_eps5, b_eps6
  real(kind=CUSTOM_REAL) :: eps_trace_over_3,b_eps_trace_over_3
  real(kind=CUSTOM_REAL) :: prod1, prod2, prod3, prod4, prod5, prod6, prod7, prod8, prod9, &
                            prod10, prod11, prod12, prod13, prod14, prod15, prod16, prod17, prod18, prod19, prod20, prod21

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE

! in principle there should also probably be a _noDev() call here as well
! and a "if(USE_DEVILLE_PRODUCTS_VAL) then" test, but for now it is not implemented
! by lack of time (and nobody uses NGLL /= 5 anyway, thus in practice USE_DEVILLE_PRODUCTS_VAL is always true...)
    call compute_element_strain_undo_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE,&
                                             b_displ_CRUST_MANTLE,ibool_crust_mantle,hprime_xx,hprime_xxT,&
                                             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                                             b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)

  ! For anisotropic kernels
  if (ANISOTROPIC_KL) then

#ifdef FORCE_VECTORIZATION
        do ijk = 1, NGLLCUBE
          iglob = ibool_crust_mantle(ijk,1,1,ispec)

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
          rho_kl_crust_mantle(ijk,1,1,ispec) =  rho_kl_crust_mantle(ijk,1,1,ispec) &
             + deltat * (accel_crust_mantle(1,iglob) * b_displ_crust_mantle(1,iglob) &
             + accel_crust_mantle(2,iglob) * b_displ_crust_mantle(2,iglob) &
             + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob))

! compute the 21 strain products

!! DK DK July 2013: inlined the call to this subroutine below to speed up the code, thus the subroutine is now unused
!       call compute_strain_product(prod,eps_trace_over_3_crust_mantle(i,j,k,ispec),epsilondev_crust_mantle(:,i,j,k,ispec), &
!                                   b_eps_trace_over_3_loc_matrix(i,j,k),b_epsilondev_loc_matrix(:,i,j,k))

  ! Purpose: compute the 21 strain products at a grid point
  ! (ispec,i,j,k fixed) and at a time t to compute then the kernels cij_kl (Voigt notation)
  ! (eq. 15 of Tromp et al., 2005)
  ! prod(1)=eps11*eps11 -> c11, prod(2)=eps11eps22 -> c12, prod(3)=eps11eps33 -> c13, ...
  ! prod(7)=eps22*eps22 -> c22, prod(8)=eps22eps33 -> c23, prod(9)=eps22eps23 -> c24, ...
  ! prod(19)=eps13*eps13 -> c55, prod(20)=eps13eps12 -> c56, prod(21)=eps12eps12 -> c66
  ! This then gives how the 21 kernels are organized

  ! Building of the local matrix of the strain tensor
  ! for the adjoint field and the regular backward field
  eps_trace_over_3 = eps_trace_over_3_crust_mantle(ijk,1,1,ispec)
  b_eps_trace_over_3 = b_eps_trace_over_3_loc_matrix(ijk,1,1)

  eps1=epsilondev_crust_mantle(1,ijk,1,1,ispec)+eps_trace_over_3  ! eps11
  eps2=epsilondev_crust_mantle(2,ijk,1,1,ispec)+eps_trace_over_3  ! eps22
  eps3=-(eps1+eps2)+3._CUSTOM_REAL*eps_trace_over_3               ! eps33
  eps4=epsilondev_crust_mantle(5,ijk,1,1,ispec)                   ! eps23
  eps5=epsilondev_crust_mantle(4,ijk,1,1,ispec)                   ! eps13
  eps6=epsilondev_crust_mantle(3,ijk,1,1,ispec)                   ! eps12

  b_eps1=b_epsilondev_loc_matrix(1,ijk,1,1)+b_eps_trace_over_3
  b_eps2=b_epsilondev_loc_matrix(2,ijk,1,1)+b_eps_trace_over_3
  b_eps3=-(b_eps1+b_eps2)+3._CUSTOM_REAL*b_eps_trace_over_3
  b_eps4=b_epsilondev_loc_matrix(5,ijk,1,1)
  b_eps5=b_epsilondev_loc_matrix(4,ijk,1,1)
  b_eps6=b_epsilondev_loc_matrix(3,ijk,1,1)

  ! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
!! DK DK July 2013: manually unrolled the calculations to speed up the code
  prod1 = eps1*b_eps1
  prod2 = eps1*b_eps2 + eps2*b_eps1
  prod3 = eps1*b_eps3 + eps3*b_eps1
  prod4 = (eps1*b_eps4 + eps4*b_eps1)*2._CUSTOM_REAL
  prod5 = (eps1*b_eps5 + eps5*b_eps1)*2._CUSTOM_REAL
  prod6 = (eps1*b_eps6 + eps6*b_eps1)*2._CUSTOM_REAL
  prod7 = eps2*b_eps2
  prod8 = eps2*b_eps3 + eps3*b_eps2
  prod9 = (eps2*b_eps4 + eps4*b_eps2)*2._CUSTOM_REAL
  prod10 = (eps2*b_eps5 + eps5*b_eps2)*2._CUSTOM_REAL
  prod11 = (eps2*b_eps6 + eps6*b_eps2)*2._CUSTOM_REAL
  prod12 = eps3*b_eps3
  prod13 = (eps3*b_eps4 + eps4*b_eps3)*2._CUSTOM_REAL
  prod14 = (eps3*b_eps5 + eps5*b_eps3)*2._CUSTOM_REAL
  prod15 = (eps3*b_eps6 + eps6*b_eps3)*2._CUSTOM_REAL
  prod16 = eps4*b_eps4*4._CUSTOM_REAL
  prod17 = (eps4*b_eps5 + eps5*b_eps4)*4._CUSTOM_REAL
  prod18 = (eps4*b_eps6 + eps6*b_eps4)*4._CUSTOM_REAL
  prod19 = eps5*b_eps5*4._CUSTOM_REAL
  prod20 = (eps5*b_eps6 + eps6*b_eps5)*4._CUSTOM_REAL
  prod21 = eps6*b_eps6*4._CUSTOM_REAL

          ! do not use a ":" array syntax for the first index below otherwise
          ! most compilers will not be able to vectorize the outer loop and the code will be slower
          cijkl_kl_crust_mantle( 1,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 1,ijk,1,1,ispec) + deltat * prod1
          cijkl_kl_crust_mantle( 2,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 2,ijk,1,1,ispec) + deltat * prod2
          cijkl_kl_crust_mantle( 3,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 3,ijk,1,1,ispec) + deltat * prod3
          cijkl_kl_crust_mantle( 4,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 4,ijk,1,1,ispec) + deltat * prod4
          cijkl_kl_crust_mantle( 5,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 5,ijk,1,1,ispec) + deltat * prod5
          cijkl_kl_crust_mantle( 6,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 6,ijk,1,1,ispec) + deltat * prod6
          cijkl_kl_crust_mantle( 7,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 7,ijk,1,1,ispec) + deltat * prod7
          cijkl_kl_crust_mantle( 8,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 8,ijk,1,1,ispec) + deltat * prod8
          cijkl_kl_crust_mantle( 9,ijk,1,1,ispec) = cijkl_kl_crust_mantle( 9,ijk,1,1,ispec) + deltat * prod9

          cijkl_kl_crust_mantle(10,ijk,1,1,ispec) = cijkl_kl_crust_mantle(10,ijk,1,1,ispec) + deltat * prod10
          cijkl_kl_crust_mantle(11,ijk,1,1,ispec) = cijkl_kl_crust_mantle(11,ijk,1,1,ispec) + deltat * prod11
          cijkl_kl_crust_mantle(12,ijk,1,1,ispec) = cijkl_kl_crust_mantle(12,ijk,1,1,ispec) + deltat * prod12
          cijkl_kl_crust_mantle(13,ijk,1,1,ispec) = cijkl_kl_crust_mantle(13,ijk,1,1,ispec) + deltat * prod13
          cijkl_kl_crust_mantle(14,ijk,1,1,ispec) = cijkl_kl_crust_mantle(14,ijk,1,1,ispec) + deltat * prod14
          cijkl_kl_crust_mantle(15,ijk,1,1,ispec) = cijkl_kl_crust_mantle(15,ijk,1,1,ispec) + deltat * prod15
          cijkl_kl_crust_mantle(16,ijk,1,1,ispec) = cijkl_kl_crust_mantle(16,ijk,1,1,ispec) + deltat * prod16
          cijkl_kl_crust_mantle(17,ijk,1,1,ispec) = cijkl_kl_crust_mantle(17,ijk,1,1,ispec) + deltat * prod17
          cijkl_kl_crust_mantle(18,ijk,1,1,ispec) = cijkl_kl_crust_mantle(18,ijk,1,1,ispec) + deltat * prod18
          cijkl_kl_crust_mantle(19,ijk,1,1,ispec) = cijkl_kl_crust_mantle(19,ijk,1,1,ispec) + deltat * prod19

          cijkl_kl_crust_mantle(20,ijk,1,1,ispec) = cijkl_kl_crust_mantle(20,ijk,1,1,ispec) + deltat * prod20
          cijkl_kl_crust_mantle(21,ijk,1,1,ispec) = cijkl_kl_crust_mantle(21,ijk,1,1,ispec) + deltat * prod21
        enddo
#else
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
             + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob))

! compute the 21 strain products

!! DK DK July 2013: inlined the call to this subroutine below to speed up the code, thus the subroutine is now unused
!       call compute_strain_product(prod,eps_trace_over_3_crust_mantle(i,j,k,ispec),epsilondev_crust_mantle(:,i,j,k,ispec), &
!                                   b_eps_trace_over_3_loc_matrix(i,j,k),b_epsilondev_loc_matrix(:,i,j,k))

  ! Purpose: compute the 21 strain products at a grid point
  ! (ispec,i,j,k fixed) and at a time t to compute then the kernels cij_kl (Voigt notation)
  ! (eq. 15 of Tromp et al., 2005)
  ! prod(1)=eps11*eps11 -> c11, prod(2)=eps11eps22 -> c12, prod(3)=eps11eps33 -> c13, ...
  ! prod(7)=eps22*eps22 -> c22, prod(8)=eps22eps33 -> c23, prod(9)=eps22eps23 -> c24, ...
  ! prod(19)=eps13*eps13 -> c55, prod(20)=eps13eps12 -> c56, prod(21)=eps12eps12 -> c66
  ! This then gives how the 21 kernels are organized

  ! Building of the local matrix of the strain tensor
  ! for the adjoint field and the regular backward field
  eps_trace_over_3 = eps_trace_over_3_crust_mantle(i,j,k,ispec)
  b_eps_trace_over_3 = b_eps_trace_over_3_loc_matrix(i,j,k)

  eps1=epsilondev_crust_mantle(1,i,j,k,ispec)+eps_trace_over_3  ! eps11
  eps2=epsilondev_crust_mantle(2,i,j,k,ispec)+eps_trace_over_3  ! eps22
  eps3=-(eps1+eps2)+3._CUSTOM_REAL*eps_trace_over_3             ! eps33
  eps4=epsilondev_crust_mantle(5,i,j,k,ispec)                   ! eps23
  eps5=epsilondev_crust_mantle(4,i,j,k,ispec)                   ! eps13
  eps6=epsilondev_crust_mantle(3,i,j,k,ispec)                   ! eps12

  b_eps1=b_epsilondev_loc_matrix(1,i,j,k)+b_eps_trace_over_3
  b_eps2=b_epsilondev_loc_matrix(2,i,j,k)+b_eps_trace_over_3
  b_eps3=-(b_eps1+b_eps2)+3._CUSTOM_REAL*b_eps_trace_over_3
  b_eps4=b_epsilondev_loc_matrix(5,i,j,k)
  b_eps5=b_epsilondev_loc_matrix(4,i,j,k)
  b_eps6=b_epsilondev_loc_matrix(3,i,j,k)

  ! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
!! DK DK July 2013: manually unrolled the calculations to speed up the code
  prod1 = eps1*b_eps1
  prod2 = eps1*b_eps2 + eps2*b_eps1
  prod3 = eps1*b_eps3 + eps3*b_eps1
  prod4 = (eps1*b_eps4 + eps4*b_eps1)*2._CUSTOM_REAL
  prod5 = (eps1*b_eps5 + eps5*b_eps1)*2._CUSTOM_REAL
  prod6 = (eps1*b_eps6 + eps6*b_eps1)*2._CUSTOM_REAL
  prod7 = eps2*b_eps2
  prod8 = eps2*b_eps3 + eps3*b_eps2
  prod9 = (eps2*b_eps4 + eps4*b_eps2)*2._CUSTOM_REAL
  prod10 = (eps2*b_eps5 + eps5*b_eps2)*2._CUSTOM_REAL
  prod11 = (eps2*b_eps6 + eps6*b_eps2)*2._CUSTOM_REAL
  prod12 = eps3*b_eps3
  prod13 = (eps3*b_eps4 + eps4*b_eps3)*2._CUSTOM_REAL
  prod14 = (eps3*b_eps5 + eps5*b_eps3)*2._CUSTOM_REAL
  prod15 = (eps3*b_eps6 + eps6*b_eps3)*2._CUSTOM_REAL
  prod16 = eps4*b_eps4*4._CUSTOM_REAL
  prod17 = (eps4*b_eps5 + eps5*b_eps4)*4._CUSTOM_REAL
  prod18 = (eps4*b_eps6 + eps6*b_eps4)*4._CUSTOM_REAL
  prod19 = eps5*b_eps5*4._CUSTOM_REAL
  prod20 = (eps5*b_eps6 + eps6*b_eps5)*4._CUSTOM_REAL
  prod21 = eps6*b_eps6*4._CUSTOM_REAL

          ! do not use a ":" array syntax for the first index below otherwise
          ! most compilers will not be able to vectorize the outer loop and the code will be slower
          cijkl_kl_crust_mantle( 1,i,j,k,ispec) = cijkl_kl_crust_mantle( 1,i,j,k,ispec) + deltat * prod1
          cijkl_kl_crust_mantle( 2,i,j,k,ispec) = cijkl_kl_crust_mantle( 2,i,j,k,ispec) + deltat * prod2
          cijkl_kl_crust_mantle( 3,i,j,k,ispec) = cijkl_kl_crust_mantle( 3,i,j,k,ispec) + deltat * prod3
          cijkl_kl_crust_mantle( 4,i,j,k,ispec) = cijkl_kl_crust_mantle( 4,i,j,k,ispec) + deltat * prod4
          cijkl_kl_crust_mantle( 5,i,j,k,ispec) = cijkl_kl_crust_mantle( 5,i,j,k,ispec) + deltat * prod5
          cijkl_kl_crust_mantle( 6,i,j,k,ispec) = cijkl_kl_crust_mantle( 6,i,j,k,ispec) + deltat * prod6
          cijkl_kl_crust_mantle( 7,i,j,k,ispec) = cijkl_kl_crust_mantle( 7,i,j,k,ispec) + deltat * prod7
          cijkl_kl_crust_mantle( 8,i,j,k,ispec) = cijkl_kl_crust_mantle( 8,i,j,k,ispec) + deltat * prod8
          cijkl_kl_crust_mantle( 9,i,j,k,ispec) = cijkl_kl_crust_mantle( 9,i,j,k,ispec) + deltat * prod9

          cijkl_kl_crust_mantle(10,i,j,k,ispec) = cijkl_kl_crust_mantle(10,i,j,k,ispec) + deltat * prod10
          cijkl_kl_crust_mantle(11,i,j,k,ispec) = cijkl_kl_crust_mantle(11,i,j,k,ispec) + deltat * prod11
          cijkl_kl_crust_mantle(12,i,j,k,ispec) = cijkl_kl_crust_mantle(12,i,j,k,ispec) + deltat * prod12
          cijkl_kl_crust_mantle(13,i,j,k,ispec) = cijkl_kl_crust_mantle(13,i,j,k,ispec) + deltat * prod13
          cijkl_kl_crust_mantle(14,i,j,k,ispec) = cijkl_kl_crust_mantle(14,i,j,k,ispec) + deltat * prod14
          cijkl_kl_crust_mantle(15,i,j,k,ispec) = cijkl_kl_crust_mantle(15,i,j,k,ispec) + deltat * prod15
          cijkl_kl_crust_mantle(16,i,j,k,ispec) = cijkl_kl_crust_mantle(16,i,j,k,ispec) + deltat * prod16
          cijkl_kl_crust_mantle(17,i,j,k,ispec) = cijkl_kl_crust_mantle(17,i,j,k,ispec) + deltat * prod17
          cijkl_kl_crust_mantle(18,i,j,k,ispec) = cijkl_kl_crust_mantle(18,i,j,k,ispec) + deltat * prod18
          cijkl_kl_crust_mantle(19,i,j,k,ispec) = cijkl_kl_crust_mantle(19,i,j,k,ispec) + deltat * prod19

          cijkl_kl_crust_mantle(20,i,j,k,ispec) = cijkl_kl_crust_mantle(20,i,j,k,ispec) + deltat * prod20
          cijkl_kl_crust_mantle(21,i,j,k,ispec) = cijkl_kl_crust_mantle(21,i,j,k,ispec) + deltat * prod21
        enddo
      enddo
    enddo
#endif

  else

#ifdef FORCE_VECTORIZATION
        do ijk = 1, NGLLCUBE
          iglob = ibool_crust_mantle(ijk,1,1,ispec)

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
          rho_kl_crust_mantle(ijk,1,1,ispec) =  rho_kl_crust_mantle(ijk,1,1,ispec) &
             + deltat * (accel_crust_mantle(1,iglob) * b_displ_crust_mantle(1,iglob) &
             + accel_crust_mantle(2,iglob) * b_displ_crust_mantle(2,iglob) &
             + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob))

          ! kernel for shear modulus, see e.g. Tromp et al. (2005), equation (17)
          ! note: multiplication with 2*mu(x) will be done after the time loop
          beta_kl_crust_mantle(ijk,1,1,ispec) =  beta_kl_crust_mantle(ijk,1,1,ispec) &
             + deltat * (epsilondev_crust_mantle(1,ijk,1,1,ispec)*b_epsilondev_loc_matrix(1,ijk,1,1) &
             + epsilondev_crust_mantle(2,ijk,1,1,ispec)*b_epsilondev_loc_matrix(2,ijk,1,1) &
             + (epsilondev_crust_mantle(1,ijk,1,1,ispec)+epsilondev_crust_mantle(2,ijk,1,1,ispec)) &
             * (b_epsilondev_loc_matrix(1,ijk,1,1)+b_epsilondev_loc_matrix(2,ijk,1,1)) &
             + 2.d0 * (epsilondev_crust_mantle(3,ijk,1,1,ispec)*b_epsilondev_loc_matrix(3,ijk,1,1) &
             + epsilondev_crust_mantle(4,ijk,1,1,ispec)*b_epsilondev_loc_matrix(4,ijk,1,1) + &
              epsilondev_crust_mantle(5,ijk,1,1,ispec)*b_epsilondev_loc_matrix(5,ijk,1,1)))

          ! kernel for bulk modulus, see e.g. Tromp et al. (2005), equation (18)
          ! note: multiplication with kappa(x) will be done after the time loop
          alpha_kl_crust_mantle(ijk,1,1,ispec) = alpha_kl_crust_mantle(ijk,1,1,ispec) &
             + deltat * (9.d0 * eps_trace_over_3_crust_mantle(ijk,1,1,ispec) &
                              * b_eps_trace_over_3_loc_matrix(ijk,1,1))
        enddo
#else
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
             + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob))

          ! kernel for shear modulus, see e.g. Tromp et al. (2005), equation (17)
          ! note: multiplication with 2*mu(x) will be done after the time loop
          beta_kl_crust_mantle(i,j,k,ispec) =  beta_kl_crust_mantle(i,j,k,ispec) &
             + deltat * (epsilondev_crust_mantle(1,i,j,k,ispec)*b_epsilondev_loc_matrix(1,i,j,k) &
             + epsilondev_crust_mantle(2,i,j,k,ispec)*b_epsilondev_loc_matrix(2,i,j,k) &
             + (epsilondev_crust_mantle(1,i,j,k,ispec)+epsilondev_crust_mantle(2,i,j,k,ispec)) &
             * (b_epsilondev_loc_matrix(1,i,j,k)+b_epsilondev_loc_matrix(2,i,j,k)) &
             + 2.d0 * (epsilondev_crust_mantle(3,i,j,k,ispec)*b_epsilondev_loc_matrix(3,i,j,k) &
             + epsilondev_crust_mantle(4,i,j,k,ispec)*b_epsilondev_loc_matrix(4,i,j,k) + &
              epsilondev_crust_mantle(5,i,j,k,ispec)*b_epsilondev_loc_matrix(5,i,j,k)))

          ! kernel for bulk modulus, see e.g. Tromp et al. (2005), equation (18)
          ! note: multiplication with kappa(x) will be done after the time loop
          alpha_kl_crust_mantle(i,j,k,ispec) = alpha_kl_crust_mantle(i,j,k,ispec) &
             + deltat * (9.d0 * eps_trace_over_3_crust_mantle(i,j,k,ispec) &
                              * b_eps_trace_over_3_loc_matrix(i,j,k))
        enddo
      enddo
    enddo
#endif

  endif ! of if ANISOTROPIC_KL

  enddo

  end subroutine compute_kernels_crust_mantle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_outer_core(ibool_outer_core, &
                        xix_outer_core,xiy_outer_core,xiz_outer_core, &
                        etax_outer_core,etay_outer_core,etaz_outer_core, &
                        gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                        hprime_xx,hprime_xxT, &
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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT

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
  real(kind=CUSTOM_REAL) :: b_div_displ_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: tempx1,tempx2,tempx3

  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: A1_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: C1_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)

  integer :: i,j,k,ispec,iglob

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

! in principle there should also probably be a _noDev() call here as well
! and a "if(USE_DEVILLE_PRODUCTS_VAL) then" test, but for now it is not implemented
! by lack of time (and nobody uses NGLL /= 5 anyway, thus in practice USE_DEVILLE_PRODUCTS_VAL is always true...)

  ! outer core, compute the actual displacement and acceleration
  do ispec = 1, NSPEC_OUTER_CORE

!----------------------------------------------------------------------------------------

    ! store the field locally
#ifdef FORCE_VECTORIZATION
        do ijk=1,NGLLCUBE
          dummyx_loc(ijk,1,1) = b_displ_outer_core(ibool_outer_core(ijk,1,1,ispec))
        enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          dummyx_loc(i,j,k) = b_displ_outer_core(ibool_outer_core(i,j,k,ispec))
        enddo
      enddo
    enddo
#endif

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

    ! get derivatives of velocity potential with respect to x, y and z
#ifdef FORCE_VECTORIZATION
        do ijk=1,NGLLCUBE
          xixl = xix_outer_core(ijk,1,1,ispec)
          xiyl = xiy_outer_core(ijk,1,1,ispec)
          xizl = xiz_outer_core(ijk,1,1,ispec)
          etaxl = etax_outer_core(ijk,1,1,ispec)
          etayl = etay_outer_core(ijk,1,1,ispec)
          etazl = etaz_outer_core(ijk,1,1,ispec)
          gammaxl = gammax_outer_core(ijk,1,1,ispec)
          gammayl = gammay_outer_core(ijk,1,1,ispec)
          gammazl = gammaz_outer_core(ijk,1,1,ispec)

          b_vector_displ_outer_core(1,ijk,1,1) = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
          b_vector_displ_outer_core(2,ijk,1,1) = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
          b_vector_displ_outer_core(3,ijk,1,1) = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)
        enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          xixl = xix_outer_core(i,j,k,ispec)
          xiyl = xiy_outer_core(i,j,k,ispec)
          xizl = xiz_outer_core(i,j,k,ispec)
          etaxl = etax_outer_core(i,j,k,ispec)
          etayl = etay_outer_core(i,j,k,ispec)
          etazl = etaz_outer_core(i,j,k,ispec)
          gammaxl = gammax_outer_core(i,j,k,ispec)
          gammayl = gammay_outer_core(i,j,k,ispec)
          gammazl = gammaz_outer_core(i,j,k,ispec)

          b_vector_displ_outer_core(1,i,j,k) = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
          b_vector_displ_outer_core(2,i,j,k) = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
          b_vector_displ_outer_core(3,i,j,k) = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)
        enddo
      enddo
    enddo
#endif

!----------------------------------------------------------------------------------------

    ! store the field locally
#ifdef FORCE_VECTORIZATION
        do ijk=1,NGLLCUBE
          dummyx_loc(ijk,1,1) = accel_outer_core(ibool_outer_core(ijk,1,1,ispec))
        enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          dummyx_loc(i,j,k) = accel_outer_core(ibool_outer_core(i,j,k,ispec))
        enddo
      enddo
    enddo
#endif

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

    ! get derivatives of velocity potential with respect to x, y and z
#ifdef FORCE_VECTORIZATION
        do ijk=1,NGLLCUBE
          xixl = xix_outer_core(ijk,1,1,ispec)
          xiyl = xiy_outer_core(ijk,1,1,ispec)
          xizl = xiz_outer_core(ijk,1,1,ispec)
          etaxl = etax_outer_core(ijk,1,1,ispec)
          etayl = etay_outer_core(ijk,1,1,ispec)
          etazl = etaz_outer_core(ijk,1,1,ispec)
          gammaxl = gammax_outer_core(ijk,1,1,ispec)
          gammayl = gammay_outer_core(ijk,1,1,ispec)
          gammazl = gammaz_outer_core(ijk,1,1,ispec)

          vector_accel_outer_core(1,ijk,1,1) = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
          vector_accel_outer_core(2,ijk,1,1) = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
          vector_accel_outer_core(3,ijk,1,1) = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)
        enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          xixl = xix_outer_core(i,j,k,ispec)
          xiyl = xiy_outer_core(i,j,k,ispec)
          xizl = xiz_outer_core(i,j,k,ispec)
          etaxl = etax_outer_core(i,j,k,ispec)
          etayl = etay_outer_core(i,j,k,ispec)
          etazl = etaz_outer_core(i,j,k,ispec)
          gammaxl = gammax_outer_core(i,j,k,ispec)
          gammayl = gammay_outer_core(i,j,k,ispec)
          gammazl = gammaz_outer_core(i,j,k,ispec)

          vector_accel_outer_core(1,i,j,k) = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
          vector_accel_outer_core(2,i,j,k) = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
          vector_accel_outer_core(3,i,j,k) = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)
        enddo
      enddo
    enddo
#endif

!----------------------------------------------------------------------------------------

    ! store the field locally
#ifdef FORCE_VECTORIZATION
        do ijk=1,NGLLCUBE
          dummyx_loc(ijk,1,1) = displ_outer_core(ibool_outer_core(ijk,1,1,ispec))
        enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          dummyx_loc(i,j,k) = displ_outer_core(ibool_outer_core(i,j,k,ispec))
        enddo
      enddo
    enddo
#endif

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

    ! get derivatives of velocity potential with respect to x, y and z
#ifdef FORCE_VECTORIZATION
        do ijk=1,NGLLCUBE
          xixl = xix_outer_core(ijk,1,1,ispec)
          xiyl = xiy_outer_core(ijk,1,1,ispec)
          xizl = xiz_outer_core(ijk,1,1,ispec)
          etaxl = etax_outer_core(ijk,1,1,ispec)
          etayl = etay_outer_core(ijk,1,1,ispec)
          etazl = etaz_outer_core(ijk,1,1,ispec)
          gammaxl = gammax_outer_core(ijk,1,1,ispec)
          gammayl = gammay_outer_core(ijk,1,1,ispec)
          gammazl = gammaz_outer_core(ijk,1,1,ispec)

          vector_displ_outer_core(1,ijk,1,1) = xixl*tempx1(ijk,1,1) + etaxl*tempx2(ijk,1,1) + gammaxl*tempx3(ijk,1,1)
          vector_displ_outer_core(2,ijk,1,1) = xiyl*tempx1(ijk,1,1) + etayl*tempx2(ijk,1,1) + gammayl*tempx3(ijk,1,1)
          vector_displ_outer_core(3,ijk,1,1) = xizl*tempx1(ijk,1,1) + etazl*tempx2(ijk,1,1) + gammazl*tempx3(ijk,1,1)
        enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          xixl = xix_outer_core(i,j,k,ispec)
          xiyl = xiy_outer_core(i,j,k,ispec)
          xizl = xiz_outer_core(i,j,k,ispec)
          etaxl = etax_outer_core(i,j,k,ispec)
          etayl = etay_outer_core(i,j,k,ispec)
          etazl = etaz_outer_core(i,j,k,ispec)
          gammaxl = gammax_outer_core(i,j,k,ispec)
          gammayl = gammay_outer_core(i,j,k,ispec)
          gammazl = gammaz_outer_core(i,j,k,ispec)

          vector_displ_outer_core(1,i,j,k) = xixl*tempx1(i,j,k) + etaxl*tempx2(i,j,k) + gammaxl*tempx3(i,j,k)
          vector_displ_outer_core(2,i,j,k) = xiyl*tempx1(i,j,k) + etayl*tempx2(i,j,k) + gammayl*tempx3(i,j,k)
          vector_displ_outer_core(3,i,j,k) = xizl*tempx1(i,j,k) + etazl*tempx2(i,j,k) + gammazl*tempx3(i,j,k)
        enddo
      enddo
    enddo
#endif

!----------------------------------------------------------------------------------------

#ifdef FORCE_VECTORIZATION
        do ijk = 1, NGLLCUBE
          rho_kl_outer_core(ijk,1,1,ispec) = rho_kl_outer_core(ijk,1,1,ispec) &
             + deltat * (vector_accel_outer_core(1,ijk,1,1) * b_vector_displ_outer_core(1,ijk,1,1) + &
                         vector_accel_outer_core(2,ijk,1,1) * b_vector_displ_outer_core(2,ijk,1,1) + &
                         vector_accel_outer_core(3,ijk,1,1) * b_vector_displ_outer_core(3,ijk,1,1))

          kappal = rhostore_outer_core(ijk,1,1,ispec) / kappavstore_outer_core(ijk,1,1,ispec)

          iglob = ibool_outer_core(ijk,1,1,ispec)

          div_displ_outer_core(ijk,1,1,ispec) =  kappal * accel_outer_core(iglob)
          b_div_displ_outer_core =  kappal * b_accel_outer_core(iglob)

          alpha_kl_outer_core(ijk,1,1,ispec) = alpha_kl_outer_core(ijk,1,1,ispec) &
             + deltat * div_displ_outer_core(ijk,1,1,ispec) * b_div_displ_outer_core
        enddo
#else
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          rho_kl_outer_core(i,j,k,ispec) = rho_kl_outer_core(i,j,k,ispec) &
!            + deltat * dot_product(vector_accel_outer_core(:,i,j,k), b_vector_displ_outer_core(:,i,j,k))
!! DK DK July 2013: replaces dot_product() with an unrolled expression, otherwise most compilers
!! DK DK July 2013: will try to vectorize this rather than the outer loop, resulting in a much slower code
             + deltat * (vector_accel_outer_core(1,i,j,k) * b_vector_displ_outer_core(1,i,j,k) + &
                         vector_accel_outer_core(2,i,j,k) * b_vector_displ_outer_core(2,i,j,k) + &
                         vector_accel_outer_core(3,i,j,k) * b_vector_displ_outer_core(3,i,j,k))

          kappal = rhostore_outer_core(i,j,k,ispec) / kappavstore_outer_core(i,j,k,ispec)

          iglob = ibool_outer_core(i,j,k,ispec)

          div_displ_outer_core(i,j,k,ispec) =  kappal * accel_outer_core(iglob)
          b_div_displ_outer_core =  kappal * b_accel_outer_core(iglob)

          alpha_kl_outer_core(i,j,k,ispec) = alpha_kl_outer_core(i,j,k,ispec) &
             + deltat * div_displ_outer_core(i,j,k,ispec) * b_div_displ_outer_core
        enddo
      enddo
    enddo
#endif

  enddo

  end subroutine compute_kernels_outer_core

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_inner_core(ibool_inner_core, &
                          rho_kl_inner_core,beta_kl_inner_core, &
                          alpha_kl_inner_core, &
                          accel_inner_core,b_displ_inner_core, &
                          deltat,hprime_xx,hprime_xxT,&
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,&
                          epsilondev_inner_core,eps_trace_over_3_inner_core)

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

  real(kind=CUSTOM_REAL) deltat
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprime_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: eps_trace_over_3_inner_core

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: b_epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: b_eps_trace_over_3_loc_matrix

  integer :: ispec,iglob

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! inner_core
  do ispec = 1, NSPEC_INNER_CORE

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
             + deltat * (epsilondev_inner_core(1,ijk,1,1,ispec)*b_epsilondev_loc_matrix(1,ijk,1,1) &
                + epsilondev_inner_core(2,ijk,1,1,ispec)*b_epsilondev_loc_matrix(2,ijk,1,1) &
                + (epsilondev_inner_core(1,ijk,1,1,ispec)+epsilondev_inner_core(2,ijk,1,1,ispec)) * &
                  (b_epsilondev_loc_matrix(1,ijk,1,1)+b_epsilondev_loc_matrix(2,ijk,1,1)) &
                + 2.d0 * (epsilondev_inner_core(3,ijk,1,1,ispec)*b_epsilondev_loc_matrix(3,ijk,1,1) &
                + epsilondev_inner_core(4,ijk,1,1,ispec)*b_epsilondev_loc_matrix(4,ijk,1,1) &
                + epsilondev_inner_core(5,ijk,1,1,ispec)*b_epsilondev_loc_matrix(5,ijk,1,1)))

          alpha_kl_inner_core(ijk,1,1,ispec) = alpha_kl_inner_core(ijk,1,1,ispec) &
                + deltat * (9.d0 * eps_trace_over_3_inner_core(ijk,1,1,ispec) * b_eps_trace_over_3_loc_matrix(ijk,1,1))
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
             + deltat * (epsilondev_inner_core(1,i,j,k,ispec)*b_epsilondev_loc_matrix(1,i,j,k) &
                + epsilondev_inner_core(2,i,j,k,ispec)*b_epsilondev_loc_matrix(2,i,j,k) &
                + (epsilondev_inner_core(1,i,j,k,ispec)+epsilondev_inner_core(2,i,j,k,ispec)) * &
                  (b_epsilondev_loc_matrix(1,i,j,k)+b_epsilondev_loc_matrix(2,i,j,k)) &
                + 2.d0 * (epsilondev_inner_core(3,i,j,k,ispec)*b_epsilondev_loc_matrix(3,i,j,k) &
                + epsilondev_inner_core(4,i,j,k,ispec)*b_epsilondev_loc_matrix(4,i,j,k) &
                + epsilondev_inner_core(5,i,j,k,ispec)*b_epsilondev_loc_matrix(5,i,j,k)))

          alpha_kl_inner_core(i,j,k,ispec) = alpha_kl_inner_core(i,j,k,ispec) &
                + deltat * (9.d0 * eps_trace_over_3_inner_core(i,j,k,ispec) * b_eps_trace_over_3_loc_matrix(i,j,k))
        enddo
      enddo
    enddo
#endif

  enddo

  end subroutine compute_kernels_inner_core

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

!!
!!-------------------------------------------------------------------------------------------------
!!
!! compute the 21 strain products

!! DK DK July 2013: inlined the call to this subroutine above to speed up the code, thus the subroutine is now unused

!  subroutine compute_strain_product(prod,eps_trace_over_3,epsdev,&
!                                    b_eps_trace_over_3,b_epsdev)

!  ! Purpose: compute the 21 strain products at a grid point
!  ! (ispec,i,j,k fixed) and at a time t to compute then the kernels cij_kl (Voigt notation)
!  ! (eq. 15 of Tromp et al., 2005)
!  ! prod(1)=eps11*eps11 -> c11, prod(2)=eps11eps22 -> c12, prod(3)=eps11eps33 -> c13, ...
!  ! prod(7)=eps22*eps22 -> c22, prod(8)=eps22eps33 -> c23, prod(9)=eps22eps23 -> c24, ...
!  ! prod(19)=eps13*eps13 -> c55, prod(20)=eps13eps12 -> c56, prod(21)=eps12eps12 -> c66
!  ! This then gives how the 21 kernels are organized

!  implicit none
!  include  "constants.h"

!  real(kind=CUSTOM_REAL),dimension(21) :: prod
!  real(kind=CUSTOM_REAL) :: eps_trace_over_3,b_eps_trace_over_3
!  real(kind=CUSTOM_REAL),dimension(5) :: epsdev,b_epsdev

!! integer :: p,i,j
!  real(kind=CUSTOM_REAL) :: eps1, eps2, eps3, eps4, eps5, eps6, b_eps1, b_eps2, b_eps3, b_eps4, b_eps5, b_eps6

!  ! Building of the local matrix of the strain tensor
!  ! for the adjoint field and the regular backward field
!  eps1=epsdev(1)+eps_trace_over_3                    ! eps11
!  eps2=epsdev(2)+eps_trace_over_3                    ! eps22
!  eps3=-(eps1+eps2)+3._CUSTOM_REAL*eps_trace_over_3  ! eps33
!  eps4=epsdev(5)                                     ! eps23
!  eps5=epsdev(4)                                     ! eps13
!  eps6=epsdev(3)                                     ! eps12

!  b_eps1=b_epsdev(1)+b_eps_trace_over_3
!  b_eps2=b_epsdev(2)+b_eps_trace_over_3
!  b_eps3=-(b_eps1+b_eps2)+3._CUSTOM_REAL*b_eps_trace_over_3
!  b_eps4=b_epsdev(5)
!  b_eps5=b_epsdev(4)
!  b_eps6=b_epsdev(3)

!  ! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
!! p=1
!! do i=1,6
!!   do j=i,6
!!     prod(p)=eps(i)*b_eps(j)
!!     if(j>i) then
!!       prod(p)=prod(p)+eps(j)*b_eps(i)
!!       if(j>3 .and. i<4) prod(p) = prod(p) * 2.0_CUSTOM_REAL
!!     endif
!!     if(i>3) prod(p) = prod(p) * 4.0_CUSTOM_REAL
!!     p=p+1
!!   enddo
!! enddo
!!! DK DK July 2013: manually unrolled the calculations to speed up the code
!  prod(1) = eps1*b_eps1
!  prod(2) = eps1*b_eps2 + eps2*b_eps1
!  prod(3) = eps1*b_eps3 + eps3*b_eps1
!  prod(4) = (eps1*b_eps4 + eps4*b_eps1)*2._CUSTOM_REAL
!  prod(5) = (eps1*b_eps5 + eps5*b_eps1)*2._CUSTOM_REAL
!  prod(6) = (eps1*b_eps6 + eps6*b_eps1)*2._CUSTOM_REAL
!  prod(7) = eps2*b_eps2
!  prod(8) = eps2*b_eps3 + eps3*b_eps2
!  prod(9) = (eps2*b_eps4 + eps4*b_eps2)*2._CUSTOM_REAL
!  prod(10) = (eps2*b_eps5 + eps5*b_eps2)*2._CUSTOM_REAL
!  prod(11) = (eps2*b_eps6 + eps6*b_eps2)*2._CUSTOM_REAL
!  prod(12) = eps3*b_eps3
!  prod(13) = (eps3*b_eps4 + eps4*b_eps3)*2._CUSTOM_REAL
!  prod(14) = (eps3*b_eps5 + eps5*b_eps3)*2._CUSTOM_REAL
!  prod(15) = (eps3*b_eps6 + eps6*b_eps3)*2._CUSTOM_REAL
!  prod(16) = eps4*b_eps4*4._CUSTOM_REAL
!  prod(17) = (eps4*b_eps5 + eps5*b_eps4)*4._CUSTOM_REAL
!  prod(18) = (eps4*b_eps6 + eps6*b_eps4)*4._CUSTOM_REAL
!  prod(19) = eps5*b_eps5*4._CUSTOM_REAL
!  prod(20) = (eps5*b_eps6 + eps6*b_eps5)*4._CUSTOM_REAL
!  prod(21) = eps6*b_eps6*4._CUSTOM_REAL

!  end subroutine compute_strain_product

