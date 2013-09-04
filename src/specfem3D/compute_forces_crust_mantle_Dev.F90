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

  subroutine compute_forces_crust_mantle_Dev( NSPEC,NGLOB,NSPEC_ATT, &
                                              deltat, &
                                              displ_crust_mantle, &
                                              veloc_crust_mantle, &
                                              accel_crust_mantle, &
                                              phase_is_inner, &
                                              R_xx,R_yy,R_xy,R_xz,R_yz, &
                                              R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                              epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                              epsilondev_xz,epsilondev_yz, &
                                              epsilon_trace_over_3, &
                                              alphaval,betaval,gammaval, &
                                              factor_common,vnspec,is_backward_field)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  use constants_solver

  use specfem_par,only: &
    hprime_xx,hprime_xxT,hprimewgll_xx,hprimewgll_xxT,wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
    minus_gravity_table,density_table,minus_deriv_gravity_table, &
    COMPUTE_AND_STORE_STRAIN,USE_LDDRK

  use specfem_par_crustmantle,only: &
    xstore => xstore_crust_mantle,ystore => ystore_crust_mantle,zstore => zstore_crust_mantle, &
    xix => xix_crust_mantle,xiy => xiy_crust_mantle,xiz => xiz_crust_mantle, &
    etax => etax_crust_mantle,etay => etay_crust_mantle,etaz => etaz_crust_mantle, &
    gammax => gammax_crust_mantle,gammay => gammay_crust_mantle,gammaz => gammaz_crust_mantle, &
    kappavstore => kappavstore_crust_mantle,kappahstore => kappahstore_crust_mantle, &
    muvstore => muvstore_crust_mantle,muhstore => muhstore_crust_mantle, &
    eta_anisostore => eta_anisostore_crust_mantle, &
    c11store => c11store_crust_mantle,c12store => c12store_crust_mantle,c13store => c13store_crust_mantle, &
    c14store => c14store_crust_mantle,c15store => c15store_crust_mantle,c16store => c16store_crust_mantle, &
    c22store => c22store_crust_mantle,c23store => c23store_crust_mantle,c24store => c24store_crust_mantle, &
    c25store => c25store_crust_mantle,c26store => c26store_crust_mantle,c33store => c33store_crust_mantle, &
    c34store => c34store_crust_mantle,c35store => c35store_crust_mantle,c36store => c36store_crust_mantle, &
    c44store => c44store_crust_mantle,c45store => c45store_crust_mantle,c46store => c46store_crust_mantle, &
    c55store => c55store_crust_mantle,c56store => c56store_crust_mantle,c66store => c66store_crust_mantle, &
    ibool => ibool_crust_mantle, &
    ispec_is_tiso => ispec_is_tiso_crust_mantle, &
    one_minus_sum_beta => one_minus_sum_beta_crust_mantle, &
    phase_ispec_inner => phase_ispec_inner_crust_mantle, &
    nspec_outer => nspec_outer_crust_mantle, &
    nspec_inner => nspec_inner_crust_mantle

#ifdef FORCE_VECTORIZATION
  use specfem_par,only: wgllwgll_xy_3D,wgllwgll_xz_3D,wgllwgll_yz_3D
#endif

!daniel: att - debug
!  use specfem_par,only: it,NSTEP

  implicit none

  integer :: NSPEC,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL) deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: veloc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_crust_mantle

  ! variable sized array variables
  integer :: vnspec

  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATT) :: R_xx,R_yy,R_xy,R_xz,R_yz
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATT) :: &
    R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: epsilon_trace_over_3

  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  ! inner/outer element run flag
  logical :: phase_is_inner

  logical :: is_backward_field

  ! local parameters

  ! Deville
  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points,E2_m1_m2_5points,E3_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  equivalence(newtempx1,E1_m1_m2_5points)
  equivalence(newtempy1,E2_m1_m2_5points)
  equivalence(newtempz1,E3_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)
  equivalence(newtempy3,E2_mxm_m2_m1_5points)
  equivalence(newtempz3,E3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  ! for gravity
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H

  integer :: ispec,i,j,k,iglob

  integer :: num_elements,ispec_p
  integer :: iphase

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#endif

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

!  computed_elements = 0
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

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1

#ifdef FORCE_VECTORIZATION
    do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      dummyx_loc(ijk,1,1) = displ_crust_mantle(1,iglob)
      dummyy_loc(ijk,1,1) = displ_crust_mantle(2,iglob)
      dummyz_loc(ijk,1,1) = displ_crust_mantle(3,iglob)
    enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          dummyx_loc(i,j,k) = displ_crust_mantle(1,iglob)
          dummyy_loc(i,j,k) = displ_crust_mantle(2,iglob)
          dummyz_loc(i,j,k) = displ_crust_mantle(3,iglob)
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

          C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
                                  hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
                                  hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
                                  hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
                                  hprime_xx(i,5)*B2_m1_m2_5points(5,j)

          C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
                                  hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
                                  hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
                                  hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
                                  hprime_xx(i,5)*B3_m1_m2_5points(5,j)
       enddo
    enddo

    do j=1,m1
       do i=1,m1
          ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
          do k = 1,NGLLX
             tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                             dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                             dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                             dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                             dummyx_loc(i,5,k)*hprime_xxT(5,j)

             tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
                             dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
                             dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
                             dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
                             dummyy_loc(i,5,k)*hprime_xxT(5,j)

             tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
                             dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
                             dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
                             dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
                             dummyz_loc(i,5,k)*hprime_xxT(5,j)
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

          C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                      A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                      A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                      A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                      A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

          C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                      A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                      A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                      A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                      A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
       enddo
    enddo

    !
    ! compute either isotropic, transverse isotropic or anisotropic elements
    !
    if(ANISOTROPIC_3D_MANTLE_VAL) then
       ! anisotropic element
       call compute_element_aniso(ispec, &
            minus_gravity_table,density_table,minus_deriv_gravity_table, &
            xstore,ystore,zstore, &
            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
            wgll_cube, &
            c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
            c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
            c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
            ibool, &
            R_xx,R_yy,R_xy,R_xz,R_yz, &
            epsilon_trace_over_3, &
            one_minus_sum_beta,vnspec, &
            tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
            dummyx_loc,dummyy_loc,dummyz_loc, &
            epsilondev_loc, &
            rho_s_H,is_backward_field)
    else
       if(.not. ispec_is_tiso(ispec)) then
          ! isotropic element
          call compute_element_iso(ispec, &
               minus_gravity_table,density_table,minus_deriv_gravity_table, &
               xstore,ystore,zstore, &
               xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
               wgll_cube, &
               kappavstore,muvstore, &
               ibool, &
               R_xx,R_yy,R_xy,R_xz,R_yz, &
               epsilon_trace_over_3, &
               one_minus_sum_beta,vnspec, &
               tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
               dummyx_loc,dummyy_loc,dummyz_loc, &
               epsilondev_loc, &
               rho_s_H,is_backward_field)
       else
          ! transverse isotropic element
          call compute_element_tiso(ispec, &
               minus_gravity_table,density_table,minus_deriv_gravity_table, &
               xstore,ystore,zstore, &
               xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
               wgll_cube, &
               kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
               ibool, &
               R_xx,R_yy,R_xy,R_xz,R_yz, &
               epsilon_trace_over_3, &
               one_minus_sum_beta,vnspec, &
               tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
               dummyx_loc,dummyy_loc,dummyz_loc, &
               epsilondev_loc, &
               rho_s_H,is_backward_field)
       endif ! .not. ispec_is_tiso
    endif

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

        E2_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C2_m1_m2_5points(1,j) + &
                                hprimewgll_xxT(i,2)*C2_m1_m2_5points(2,j) + &
                                hprimewgll_xxT(i,3)*C2_m1_m2_5points(3,j) + &
                                hprimewgll_xxT(i,4)*C2_m1_m2_5points(4,j) + &
                                hprimewgll_xxT(i,5)*C2_m1_m2_5points(5,j)

        E3_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C3_m1_m2_5points(1,j) + &
                                hprimewgll_xxT(i,2)*C3_m1_m2_5points(2,j) + &
                                hprimewgll_xxT(i,3)*C3_m1_m2_5points(3,j) + &
                                hprimewgll_xxT(i,4)*C3_m1_m2_5points(4,j) + &
                                hprimewgll_xxT(i,5)*C3_m1_m2_5points(5,j)
      enddo
    enddo

    do i=1,m1
      do j=1,m1
        ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
        do k = 1,NGLLX
          newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
                             tempx2(i,2,k)*hprimewgll_xx(2,j) + &
                             tempx2(i,3,k)*hprimewgll_xx(3,j) + &
                             tempx2(i,4,k)*hprimewgll_xx(4,j) + &
                             tempx2(i,5,k)*hprimewgll_xx(5,j)

          newtempy2(i,j,k) = tempy2(i,1,k)*hprimewgll_xx(1,j) + &
                             tempy2(i,2,k)*hprimewgll_xx(2,j) + &
                             tempy2(i,3,k)*hprimewgll_xx(3,j) + &
                             tempy2(i,4,k)*hprimewgll_xx(4,j) + &
                             tempy2(i,5,k)*hprimewgll_xx(5,j)

          newtempz2(i,j,k) = tempz2(i,1,k)*hprimewgll_xx(1,j) + &
                             tempz2(i,2,k)*hprimewgll_xx(2,j) + &
                             tempz2(i,3,k)*hprimewgll_xx(3,j) + &
                             tempz2(i,4,k)*hprimewgll_xx(4,j) + &
                             tempz2(i,5,k)*hprimewgll_xx(5,j)
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

        E2_mxm_m2_m1_5points(i,j) = C2_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                    C2_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                    C2_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                    C2_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                    C2_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

        E3_mxm_m2_m1_5points(i,j) = C3_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                    C3_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                    C3_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                    C3_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                    C3_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
      enddo
    enddo

    ! sum contributions
#ifdef FORCE_VECTORIZATION
    do ijk=1,NGLLCUBE
      fac1 = wgllwgll_yz_3D(ijk,1,1)
      fac2 = wgllwgll_xz_3D(ijk,1,1)
      fac3 = wgllwgll_xy_3D(ijk,1,1)
      sum_terms(1,ijk,1,1) = - (fac1*newtempx1(ijk,1,1) + fac2*newtempx2(ijk,1,1) + fac3*newtempx3(ijk,1,1))
      sum_terms(2,ijk,1,1) = - (fac1*newtempy1(ijk,1,1) + fac2*newtempy2(ijk,1,1) + fac3*newtempy3(ijk,1,1))
      sum_terms(3,ijk,1,1) = - (fac1*newtempz1(ijk,1,1) + fac2*newtempz2(ijk,1,1) + fac3*newtempz3(ijk,1,1))
    enddo
    ! adds gravity terms
    if(GRAVITY_VAL) then
      do ijk = 1,NDIM*NGLLCUBE
        sum_terms(ijk,1,1,1) = sum_terms(ijk,1,1,1) + rho_s_H(ijk,1,1,1)
      enddo
    endif
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        fac1 = wgllwgll_yz(j,k)
        do i=1,NGLLX
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)
          ! sums contributions
          sum_terms(1,i,j,k) = - (fac1*newtempx1(i,j,k) + fac2*newtempx2(i,j,k) + fac3*newtempx3(i,j,k))
          sum_terms(2,i,j,k) = - (fac1*newtempy1(i,j,k) + fac2*newtempy2(i,j,k) + fac3*newtempy3(i,j,k))
          sum_terms(3,i,j,k) = - (fac1*newtempz1(i,j,k) + fac2*newtempz2(i,j,k) + fac3*newtempz3(i,j,k))
          ! adds gravity terms
          if(GRAVITY_VAL) sum_terms(:,i,j,k) = sum_terms(:,i,j,k) + rho_s_H(:,i,j,k)
        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ
#endif

    ! sum contributions from each element to the global mesh and add gravity terms
#ifdef FORCE_VECTORIZATION
! we can force vectorization using a compiler directive here because we know that there is no dependency
! inside a given spectral element, since all the global points of a local elements are different by definition
! (only common points between different elements can be the same)
! IBM, Portland PGI, and Intel and Cray syntax (Intel and Cray are the same)
!IBM* ASSERT (NODEPS)
!pgi$ ivdep
!DIR$ IVDEP
    do ijk = 1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      ! do NOT use array syntax ":" for the three statements below otherwise most compilers
      ! will not be able to vectorize the outer loop
      accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,ijk,1,1)
      accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,ijk,1,1)
      accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,ijk,1,1)
    enddo
#else
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,i,j,k)
          accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,i,j,k)
          accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,i,j,k)
        enddo
      enddo
    enddo
#endif
    ! update memory variables based upon the Runge-Kutta scheme
    ! convention for attenuation
    ! term in xx = 1
    ! term in yy = 2
    ! term in xy = 3
    ! term in xz = 4
    ! term in yz = 5
    ! term in zz not computed since zero trace
    ! This is because we only implement Q_\mu attenuation and not Q_\kappa.
    ! Note that this does *NOT* imply that there is no attenuation for P waves
    ! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
    ! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
    ! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
    ! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu

    if( ATTENUATION_VAL .and. .not. PARTIAL_PHYS_DISPERSION_ONLY_VAL ) then
      ! updates R_memory
      if( USE_LDDRK ) then
        call compute_element_att_memory_cm_lddrk(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                                 R_xx_lddrk,R_yy_lddrk,R_xy_lddrk,R_xz_lddrk,R_yz_lddrk, &
                                                 ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec,factor_common, &
                                                 c44store,muvstore, &
                                                 epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                                 epsilondev_xz,epsilondev_yz, &
                                                 epsilondev_loc, &
                                                 deltat)
      else
        call compute_element_att_memory_cm(ispec,R_xx,R_yy,R_xy,R_xz,R_yz, &
                                           ATT1_VAL,ATT2_VAL,ATT3_VAL,vnspec,factor_common, &
                                           alphaval,betaval,gammaval, &
                                           c44store,muvstore, &
                                           epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                           epsilondev_xz,epsilondev_yz, &
                                           epsilondev_loc,is_backward_field)
      endif
    endif

    ! save deviatoric strain for Runge-Kutta scheme
    if(COMPUTE_AND_STORE_STRAIN) then
      epsilondev_xx(:,:,:,ispec) = epsilondev_loc(1,:,:,:)
      epsilondev_yy(:,:,:,ispec) = epsilondev_loc(2,:,:,:)
      epsilondev_xy(:,:,:,ispec) = epsilondev_loc(3,:,:,:)
      epsilondev_xz(:,:,:,ispec) = epsilondev_loc(4,:,:,:)
      epsilondev_yz(:,:,:,ispec) = epsilondev_loc(5,:,:,:)
    endif

  enddo ! of spectral element loop NSPEC_CRUST_MANTLE

  end subroutine compute_forces_crust_mantle_Dev

