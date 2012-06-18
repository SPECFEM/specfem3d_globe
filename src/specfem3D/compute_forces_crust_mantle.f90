!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine compute_forces_crust_mantle(NSPEC,NGLOB,NSPEC_ATT, &
                                        deltat, &
                                        displ_crust_mantle, &
                                        veloc_crust_mantle, &
                                        accel_crust_mantle, &
                                        phase_is_inner, &
                                        R_xx,R_yy,R_xy,R_xz,R_yz, &
                                        epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                                        epsilondev_xz,epsilondev_yz, &
                                        epsilon_trace_over_3, &
                                        alphaval,betaval,gammaval, &
                                        factor_common,vx,vy,vz,vnspec)

  use constants_solver

  use specfem_par,only: &
    hprime_xx,hprime_yy,hprime_zz,hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
    wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
    minus_gravity_table,density_table,minus_deriv_gravity_table

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

  implicit none

  integer :: NSPEC,NGLOB,NSPEC_ATT

  ! time step
  real(kind=CUSTOM_REAL) deltat

  ! displacement, velocity and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: veloc_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB) :: accel_crust_mantle

  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
!  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory
  real(kind=CUSTOM_REAL), dimension(N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_ATT) :: &
    R_xx,R_yy,R_xy,R_xz,R_yz


!  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: &
    epsilondev_xx,epsilondev_yy,epsilondev_xy,epsilondev_xz,epsilondev_yz

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,NSPEC) :: epsilon_trace_over_3

  ! variable sized array variables for one_minus_sum_beta and factor_common
  integer vx, vy, vz, vnspec
  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  ! inner/outer element run flag
  logical :: phase_is_inner

  ! local parameters

  ! for attenuation
  real(kind=CUSTOM_REAL) one_minus_sum_beta_use,minus_sum_beta
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: factor_common_c44_muv
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

  integer ispec,iglob,ispec_strain
  integer i,j,k,l
  integer i_SLS

  ! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  real(kind=CUSTOM_REAL) two_rhovsvsq,two_rhovshsq ! two_rhovpvsq,two_rhovphsq
  real(kind=CUSTOM_REAL) four_rhovsvsq,four_rhovshsq ! four_rhovpvsq,four_rhovphsq,

  real(kind=CUSTOM_REAL) twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) two_eta_aniso,four_eta_aniso,six_eta_aniso

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal,kappavl,kappahl,muvl,muhl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) tempx1l_att,tempx2l_att,tempx3l_att
  real(kind=CUSTOM_REAL) tempy1l_att,tempy2l_att,tempy3l_att
  real(kind=CUSTOM_REAL) tempz1l_att,tempz2l_att,tempz3l_att

  real(kind=CUSTOM_REAL) duxdxl_att,duxdyl_att,duxdzl_att,duydxl_att
  real(kind=CUSTOM_REAL) duydyl_att,duydzl_att,duzdxl_att,duzdyl_att,duzdzl_att;
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl_att,duzdxl_plus_duxdzl_att,duzdyl_plus_duydzl_att;

  ! for gravity
  integer int_radius
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H

!  integer :: computed_elements
  integer :: num_elements,ispec_p
  integer :: iphase

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

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool(l,j,k,ispec)
            tempx1l = tempx1l + displ_crust_mantle(1,iglob)*hp1
            tempy1l = tempy1l + displ_crust_mantle(2,iglob)*hp1
            tempz1l = tempz1l + displ_crust_mantle(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            hp2 = hprime_yy(j,l)
            iglob = ibool(i,l,k,ispec)
            tempx2l = tempx2l + displ_crust_mantle(1,iglob)*hp2
            tempy2l = tempy2l + displ_crust_mantle(2,iglob)*hp2
            tempz2l = tempz2l + displ_crust_mantle(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool(i,j,l,ispec)
            tempx3l = tempx3l + displ_crust_mantle(1,iglob)*hp3
            tempy3l = tempy3l + displ_crust_mantle(2,iglob)*hp3
            tempz3l = tempz3l + displ_crust_mantle(3,iglob)*hp3
          enddo

          if( ATTENUATION_VAL .and. COMPUTE_AND_STORE_STRAIN ) then
             ! temporary variables used for fixing attenuation in a consistent way

             tempx1l_att = tempx1l
             tempx2l_att = tempx2l
             tempx3l_att = tempx3l

             tempy1l_att = tempy1l
             tempy2l_att = tempy2l
             tempy3l_att = tempy3l
             
             tempz1l_att = tempz1l
             tempz2l_att = tempz2l
             tempz3l_att = tempz3l

             if(ATTENUATION_NEW_VAL) then
                ! takes new routines
                ! use first order Taylor expansion of displacement for local storage of stresses 
                ! at this current time step, to fix attenuation in a consistent way
                do l=1,NGLLX
                   hp1 = hprime_xx(i,l)
                   iglob = ibool(l,j,k,ispec)
                   tempx1l_att = tempx1l_att + deltat*veloc_crust_mantle(1,iglob)*hp1
                   tempy1l_att = tempy1l_att + deltat*veloc_crust_mantle(2,iglob)*hp1
                   tempz1l_att = tempz1l_att + deltat*veloc_crust_mantle(3,iglob)*hp1

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
                   hp2 = hprime_yy(j,l)
                   iglob = ibool(i,l,k,ispec)
                   tempx2l_att = tempx2l_att + deltat*veloc_crust_mantle(1,iglob)*hp2
                   tempy2l_att = tempy2l_att + deltat*veloc_crust_mantle(2,iglob)*hp2
                   tempz2l_att = tempz2l_att + deltat*veloc_crust_mantle(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
                   hp3 = hprime_zz(k,l)
                   iglob = ibool(i,j,l,ispec)
                   tempx3l_att = tempx3l_att + deltat*veloc_crust_mantle(1,iglob)*hp3
                   tempy3l_att = tempy3l_att + deltat*veloc_crust_mantle(2,iglob)*hp3
                   tempz3l_att = tempz3l_att + deltat*veloc_crust_mantle(3,iglob)*hp3
                enddo
             endif
          endif

!         get derivatives of ux, uy and uz with respect to x, y and z

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

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

          if( ATTENUATION_VAL .and. COMPUTE_AND_STORE_STRAIN ) then
             ! temporary variables used for fixing attenuation in a consistent way
             duxdxl_att = xixl*tempx1l_att + etaxl*tempx2l_att + gammaxl*tempx3l_att
             duxdyl_att = xiyl*tempx1l_att + etayl*tempx2l_att + gammayl*tempx3l_att
             duxdzl_att = xizl*tempx1l_att + etazl*tempx2l_att + gammazl*tempx3l_att

             duydxl_att = xixl*tempy1l_att + etaxl*tempy2l_att + gammaxl*tempy3l_att
             duydyl_att = xiyl*tempy1l_att + etayl*tempy2l_att + gammayl*tempy3l_att
             duydzl_att = xizl*tempy1l_att + etazl*tempy2l_att + gammazl*tempy3l_att

             duzdxl_att = xixl*tempz1l_att + etaxl*tempz2l_att + gammaxl*tempz3l_att
             duzdyl_att = xiyl*tempz1l_att + etayl*tempz2l_att + gammayl*tempz3l_att
             duzdzl_att = xizl*tempz1l_att + etazl*tempz2l_att + gammazl*tempz3l_att

             ! precompute some sums to save CPU time
             duxdyl_plus_duydxl_att = duxdyl_att + duydxl_att
             duzdxl_plus_duxdzl_att = duzdxl_att + duxdzl_att
             duzdyl_plus_duydzl_att = duzdyl_att + duydzl_att

             ! compute deviatoric strain
             if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
                ispec_strain = 1
             else
                ispec_strain = ispec
             endif
             epsilon_trace_over_3(i,j,k,ispec_strain) = ONE_THIRD * (duxdxl_att + duydyl_att + duzdzl_att)
             epsilondev_loc(1,i,j,k) = duxdxl_att - epsilon_trace_over_3(i,j,k,ispec_strain)
             epsilondev_loc(2,i,j,k) = duydyl_att - epsilon_trace_over_3(i,j,k,ispec_strain)
             epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl_att
             epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl_att
             epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl_att
          else
             ! compute deviatoric strain
             if (COMPUTE_AND_STORE_STRAIN) then
                if(NSPEC_CRUST_MANTLE_STRAIN_ONLY == 1) then
                   ispec_strain = 1
                else
                   ispec_strain = ispec
                endif
                epsilon_trace_over_3(i,j,k,ispec_strain) = ONE_THIRD * (duxdxl + duydyl + duzdzl)
                epsilondev_loc(1,i,j,k) = duxdxl - epsilon_trace_over_3(i,j,k,ispec_strain)
                epsilondev_loc(2,i,j,k) = duydyl - epsilon_trace_over_3(i,j,k,ispec_strain)
                epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
                epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
                epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
             endif
          endif

          ! precompute terms for attenuation if needed
          if(ATTENUATION_VAL) then
            one_minus_sum_beta_use = one_minus_sum_beta(i,j,k,ispec)
            minus_sum_beta =  one_minus_sum_beta_use - 1.0
          endif

        !
        ! compute either isotropic or anisotropic elements
        !

          if(ANISOTROPIC_3D_MANTLE_VAL) then

            c11 = c11store(i,j,k,ispec)
            c12 = c12store(i,j,k,ispec)
            c13 = c13store(i,j,k,ispec)
            c14 = c14store(i,j,k,ispec)
            c15 = c15store(i,j,k,ispec)
            c16 = c16store(i,j,k,ispec)
            c22 = c22store(i,j,k,ispec)
            c23 = c23store(i,j,k,ispec)
            c24 = c24store(i,j,k,ispec)
            c25 = c25store(i,j,k,ispec)
            c26 = c26store(i,j,k,ispec)
            c33 = c33store(i,j,k,ispec)
            c34 = c34store(i,j,k,ispec)
            c35 = c35store(i,j,k,ispec)
            c36 = c36store(i,j,k,ispec)
            c44 = c44store(i,j,k,ispec)
            c45 = c45store(i,j,k,ispec)
            c46 = c46store(i,j,k,ispec)
            c55 = c55store(i,j,k,ispec)
            c56 = c56store(i,j,k,ispec)
            c66 = c66store(i,j,k,ispec)

            if(ATTENUATION_VAL) then
              mul = c44
              c11 = c11 + FOUR_THIRDS * minus_sum_beta * mul
              c12 = c12 - TWO_THIRDS * minus_sum_beta * mul
              c13 = c13 - TWO_THIRDS * minus_sum_beta * mul
              c22 = c22 + FOUR_THIRDS * minus_sum_beta * mul
              c23 = c23 - TWO_THIRDS * minus_sum_beta * mul
              c33 = c33 + FOUR_THIRDS * minus_sum_beta * mul
              c44 = c44 + minus_sum_beta * mul
              c55 = c55 + minus_sum_beta * mul
              c66 = c66 + minus_sum_beta * mul
            endif

            sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                       c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

            sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                       c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

            sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                       c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

            sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                       c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

            sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                       c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

            sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                       c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

          else

        ! do not use transverse isotropy except if element is between d220 and Moho
!            if(.not. (TRANSVERSE_ISOTROPY_VAL .and. (idoubling(ispec)==IFLAG_220_80 .or. idoubling(ispec)==IFLAG_80_MOHO))) then
            if( .not. ispec_is_tiso(ispec) ) then
        ! layer with no transverse isotropy, use kappav and muv
              kappal = kappavstore(i,j,k,ispec)
              mul = muvstore(i,j,k,ispec)

        ! use unrelaxed parameters if attenuation
              if(ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use

              lambdalplus2mul = kappal + FOUR_THIRDS * mul
              lambdal = lambdalplus2mul - 2.*mul

        ! compute stress sigma

              sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
              sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
              sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

              sigma_xy = mul*duxdyl_plus_duydxl
              sigma_xz = mul*duzdxl_plus_duxdzl
              sigma_yz = mul*duzdyl_plus_duydzl

            else

        ! use Kappa and mu from transversely isotropic model
              kappavl = kappavstore(i,j,k,ispec)
              muvl = muvstore(i,j,k,ispec)

              kappahl = kappahstore(i,j,k,ispec)
              muhl = muhstore(i,j,k,ispec)

        ! use unrelaxed parameters if attenuation
        ! eta does not need to be shifted since it is a ratio
              if(ATTENUATION_VAL) then
                muvl = muvl * one_minus_sum_beta_use
                muhl = muhl * one_minus_sum_beta_use
              endif

              rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
              rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A

              rhovsvsq = muvl  !!! that is L
              rhovshsq = muhl  !!! that is N

              eta_aniso = eta_anisostore(i,j,k,ispec)  !!! that is  F / (A - 2 L)

        ! use mesh coordinates to get theta and phi
        ! ystore and zstore contain theta and phi

              iglob = ibool(i,j,k,ispec)
              theta = ystore(iglob)
              phi = zstore(iglob)

              costheta = cos(theta)
              sintheta = sin(theta)
              cosphi = cos(phi)
              sinphi = sin(phi)

              costhetasq = costheta * costheta
              sinthetasq = sintheta * sintheta
              cosphisq = cosphi * cosphi
              sinphisq = sinphi * sinphi

              costhetafour = costhetasq * costhetasq
              sinthetafour = sinthetasq * sinthetasq
              cosphifour = cosphisq * cosphisq
              sinphifour = sinphisq * sinphisq

              costwotheta = cos(2.*theta)
              sintwotheta = sin(2.*theta)
              costwophi = cos(2.*phi)
              sintwophi = sin(2.*phi)

              cosfourtheta = cos(4.*theta)
              cosfourphi = cos(4.*phi)

              costwothetasq = costwotheta * costwotheta

              costwophisq = costwophi * costwophi
              sintwophisq = sintwophi * sintwophi

              etaminone = eta_aniso - 1.
              twoetaminone = 2. * eta_aniso - 1.

        ! precompute some products to reduce the CPU time

              two_eta_aniso = 2.*eta_aniso
              four_eta_aniso = 4.*eta_aniso
              six_eta_aniso = 6.*eta_aniso

              !two_rhovpvsq = 2.*rhovpvsq
              !two_rhovphsq = 2.*rhovphsq
              two_rhovsvsq = 2.*rhovsvsq
              two_rhovshsq = 2.*rhovshsq

              !four_rhovpvsq = 4.*rhovpvsq
              !four_rhovphsq = 4.*rhovphsq
              four_rhovsvsq = 4.*rhovsvsq
              four_rhovshsq = 4.*rhovshsq

        ! the 21 anisotropic coefficients computed using Mathematica

             c11 = rhovphsq*sinphifour + 2.*cosphisq*sinphisq* &
               (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  sinthetasq) + cosphifour* &
               (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  costhetasq*sinthetasq + rhovpvsq*sinthetafour)

             c12 = ((rhovphsq - two_rhovshsq)*(3. + cosfourphi)*costhetasq)/4. - &
              four_rhovshsq*cosphisq*costhetasq*sinphisq + &
              (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. + &
              eta_aniso*(rhovphsq - two_rhovsvsq)*(cosphifour + &
                 2.*cosphisq*costhetasq*sinphisq + sinphifour)*sinthetasq + &
              rhovpvsq*cosphisq*sinphisq*sinthetafour - &
              rhovsvsq*sintwophisq*sinthetafour

             c13 = (cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - &
                   12.*eta_aniso*rhovsvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - &
                      four_eta_aniso*rhovsvsq)*cosfourtheta))/8. + &
              sinphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
                 (rhovphsq - two_rhovshsq)*sinthetasq)

             c14 = costheta*sinphi*((cosphisq* &
                   (-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
                     (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                        four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
                (etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))*sinphisq)* sintheta

             c15 = cosphi*costheta*((cosphisq* (-rhovphsq + rhovpvsq + &
                     (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                      costwotheta))/2. + etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sintheta

             c16 = (cosphi*sinphi*(cosphisq* (-rhovphsq + rhovpvsq + &
                     (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                        four_eta_aniso*rhovsvsq)*costwotheta) + &
                  2.*etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sinthetasq)/2.

             c22 = rhovphsq*cosphifour + 2.*cosphisq*sinphisq* &
               (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  sinthetasq) + sinphifour* &
               (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                  costhetasq*sinthetasq + rhovpvsq*sinthetafour)

             c23 = ((rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - 12.*eta_aniso*rhovsvsq + &
                   (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    cosfourtheta)*sinphisq)/8. + &
              cosphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
                 (rhovphsq - two_rhovshsq)*sinthetasq)

             c24 = costheta*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
                ((-rhovphsq + rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + &
                        four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

             c25 = cosphi*costheta*((etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))* &
                 cosphisq + ((-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
                     (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                        four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

             c26 = (cosphi*sinphi*(2.*etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
                  (-rhovphsq + rhovpvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                        four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*sinthetasq)/2.

             c33 = rhovpvsq*costhetafour + 2.*(eta_aniso*(rhovphsq - two_rhovsvsq) + two_rhovsvsq)* &
               costhetasq*sinthetasq + rhovphsq*sinthetafour

             c34 = -((rhovphsq - rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq &
                       - four_eta_aniso*rhovsvsq)*costwotheta)*sinphi*sintwotheta)/4.

             c35 = -(cosphi*(rhovphsq - rhovpvsq + &
                   (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    costwotheta)*sintwotheta)/4.

             c36 = -((rhovphsq - rhovpvsq - four_rhovshsq + four_rhovsvsq + &
                   (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                    costwotheta)*sintwophi*sinthetasq)/4.

             c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
              sinphisq*(rhovsvsq*costwothetasq + &
                 (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

             c45 = ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                  four_eta_aniso*rhovsvsq + (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + &
                     4.*etaminone*rhovsvsq)*costwotheta)*sintwophi*sinthetasq)/4.

             c46 = -(cosphi*costheta*((rhovshsq - rhovsvsq)*cosphisq - &
                  ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                       four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
                          four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)* sintheta)

             c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
              cosphisq*(rhovsvsq*costwothetasq + &
                 (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

             c56 = costheta*sinphi*((cosphisq* &
                   (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                     four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
                        four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
                (-rhovshsq + rhovsvsq)*sinphisq)*sintheta

             c66 = rhovshsq*costwophisq*costhetasq - &
              2.*(rhovphsq - two_rhovshsq)*cosphisq*costhetasq*sinphisq + &
              (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. - &
              (rhovsvsq*(-6. - 2.*cosfourphi + cos(4.*phi - 2.*theta) - 2.*costwotheta + &
                   cos(2.*(2.*phi + theta)))*sinthetasq)/8. + &
              rhovpvsq*cosphisq*sinphisq*sinthetafour - &
              (eta_aniso*(rhovphsq - two_rhovsvsq)*sintwophisq*sinthetafour)/2.

        ! general expression of stress tensor for full Cijkl with 21 coefficients

             sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                       c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

             sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                       c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

             sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                       c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

             sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                       c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

             sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                       c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

             sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                       c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

            endif

          endif   ! end of test whether isotropic or anisotropic element

        ! subtract memory variables if attenuation
          if(ATTENUATION_VAL .and. ( USE_ATTENUATION_MIMIC .eqv. .false. ) ) then
            do i_SLS = 1,N_SLS
              R_xx_val = R_xx(i_SLS,i,j,k,ispec)
              R_yy_val = R_yy(i_SLS,i,j,k,ispec)
              sigma_xx = sigma_xx - R_xx_val
              sigma_yy = sigma_yy - R_yy_val
              sigma_zz = sigma_zz + R_xx_val + R_yy_val
              sigma_xy = sigma_xy - R_xy(i_SLS,i,j,k,ispec)
              sigma_xz = sigma_xz - R_xz(i_SLS,i,j,k,ispec)
              sigma_yz = sigma_yz - R_yz(i_SLS,i,j,k,ispec)
            enddo
          endif

        ! define symmetric components of sigma for gravity
          sigma_yx = sigma_xy
          sigma_zx = sigma_xz
          sigma_zy = sigma_yz

        ! compute non-symmetric terms for gravity
          if(GRAVITY_VAL) then

        ! use mesh coordinates to get theta and phi
        ! x y and z contain r theta and phi

            iglob = ibool(i,j,k,ispec)
            radius = dble(xstore(iglob))
            theta = ystore(iglob)
            phi = zstore(iglob)

            cos_theta = dcos(dble(theta))
            sin_theta = dsin(dble(theta))
            cos_phi = dcos(dble(phi))
            sin_phi = dsin(dble(phi))

        ! get g, rho and dg/dr=dg
        ! spherical components of the gravitational acceleration
        ! for efficiency replace with lookup table every 100 m in radial direction
            int_radius = nint(radius * R_EARTH_KM * 10.d0)
            minus_g = minus_gravity_table(int_radius)
            minus_dg = minus_deriv_gravity_table(int_radius)
            rho = density_table(int_radius)

        ! Cartesian components of the gravitational acceleration
            gxl = minus_g*sin_theta*cos_phi
            gyl = minus_g*sin_theta*sin_phi
            gzl = minus_g*cos_theta

        ! Cartesian components of gradient of gravitational acceleration
        ! obtained from spherical components

            minus_g_over_radius = minus_g / radius
            minus_dg_plus_g_over_radius = minus_dg - minus_g_over_radius

            cos_theta_sq = cos_theta**2
            sin_theta_sq = sin_theta**2
            cos_phi_sq = cos_phi**2
            sin_phi_sq = sin_phi**2

            Hxxl = minus_g_over_radius*(cos_phi_sq*cos_theta_sq + sin_phi_sq) + cos_phi_sq*minus_dg*sin_theta_sq
            Hyyl = minus_g_over_radius*(cos_phi_sq + cos_theta_sq*sin_phi_sq) + minus_dg*sin_phi_sq*sin_theta_sq
            Hzzl = cos_theta_sq*minus_dg + minus_g_over_radius*sin_theta_sq
            Hxyl = cos_phi*minus_dg_plus_g_over_radius*sin_phi*sin_theta_sq
            Hxzl = cos_phi*cos_theta*minus_dg_plus_g_over_radius*sin_theta
            Hyzl = cos_theta*minus_dg_plus_g_over_radius*sin_phi*sin_theta

            iglob = ibool(i,j,k,ispec)

        ! distinguish between single and double precision for reals
            if(CUSTOM_REAL == SIZE_REAL) then

        ! get displacement and multiply by density to compute G tensor
              sx_l = rho * dble(displ_crust_mantle(1,iglob))
              sy_l = rho * dble(displ_crust_mantle(2,iglob))
              sz_l = rho * dble(displ_crust_mantle(3,iglob))

        ! compute G tensor from s . g and add to sigma (not symmetric)
              sigma_xx = sigma_xx + sngl(sy_l*gyl + sz_l*gzl)
              sigma_yy = sigma_yy + sngl(sx_l*gxl + sz_l*gzl)
              sigma_zz = sigma_zz + sngl(sx_l*gxl + sy_l*gyl)

              sigma_xy = sigma_xy - sngl(sx_l * gyl)
              sigma_yx = sigma_yx - sngl(sy_l * gxl)

              sigma_xz = sigma_xz - sngl(sx_l * gzl)
              sigma_zx = sigma_zx - sngl(sz_l * gxl)

              sigma_yz = sigma_yz - sngl(sy_l * gzl)
              sigma_zy = sigma_zy - sngl(sz_l * gyl)

        ! precompute vector
              factor = dble(jacobianl) * wgll_cube(i,j,k)
              rho_s_H(1,i,j,k) = sngl(factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl))
              rho_s_H(2,i,j,k) = sngl(factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl))
              rho_s_H(3,i,j,k) = sngl(factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl))

            else

        ! get displacement and multiply by density to compute G tensor
              sx_l = rho * displ_crust_mantle(1,iglob)
              sy_l = rho * displ_crust_mantle(2,iglob)
              sz_l = rho * displ_crust_mantle(3,iglob)

        ! compute G tensor from s . g and add to sigma (not symmetric)
              sigma_xx = sigma_xx + sy_l*gyl + sz_l*gzl
              sigma_yy = sigma_yy + sx_l*gxl + sz_l*gzl
              sigma_zz = sigma_zz + sx_l*gxl + sy_l*gyl

              sigma_xy = sigma_xy - sx_l * gyl
              sigma_yx = sigma_yx - sy_l * gxl

              sigma_xz = sigma_xz - sx_l * gzl
              sigma_zx = sigma_zx - sz_l * gxl

              sigma_yz = sigma_yz - sy_l * gzl
              sigma_zy = sigma_zy - sz_l * gyl

        ! precompute vector
              factor = jacobianl * wgll_cube(i,j,k)
              rho_s_H(1,i,j,k) = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
              rho_s_H(2,i,j,k) = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
              rho_s_H(3,i,j,k) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)

            endif

          endif  ! end of section with gravity terms

        ! form dot product with test vector, non-symmetric form
          tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl)
          tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl)
          tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

          tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl)
          tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl)
          tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

          tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl)
          tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl)
          tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempy1l = 0._CUSTOM_REAL
          tempz1l = 0._CUSTOM_REAL

          tempx2l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL

          tempx3l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1l = tempx1l + tempx1(l,j,k)*fac1
            tempy1l = tempy1l + tempy1(l,j,k)*fac1
            tempz1l = tempz1l + tempz1(l,j,k)*fac1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            fac2 = hprimewgll_yy(l,j)
            tempx2l = tempx2l + tempx2(i,l,k)*fac2
            tempy2l = tempy2l + tempy2(i,l,k)*fac2
            tempz2l = tempz2l + tempz2(i,l,k)*fac2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            fac3 = hprimewgll_zz(l,k)
            tempx3l = tempx3l + tempx3(i,j,l)*fac3
            tempy3l = tempy3l + tempy3(i,j,l)*fac3
            tempz3l = tempz3l + tempz3(i,j,l)*fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          sum_terms(1,i,j,k) = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
          sum_terms(2,i,j,k) = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
          sum_terms(3,i,j,k) = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

          if(GRAVITY_VAL) sum_terms(:,i,j,k) = sum_terms(:,i,j,k) + rho_s_H(:,i,j,k)

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

! sum contributions from each element to the global mesh and add gravity terms
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool(i,j,k,ispec)
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) + sum_terms(:,i,j,k)
        enddo
      enddo
    enddo

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

    if(ATTENUATION_VAL .and. ( USE_ATTENUATION_MIMIC .eqv. .false. )) then

! use Runge-Kutta scheme to march in time
      do i_SLS = 1,N_SLS

! get coefficients for that standard linear solid
! IMPROVE we use mu_v here even if there is some anisotropy
! IMPROVE we should probably use an average value instead

        ! reformatted R_memory to handle large factor_common and reduced [alpha,beta,gamma]val
        factor_common_c44_muv = factor_common(i_SLS,:,:,:,ispec)
        if(ANISOTROPIC_3D_MANTLE_VAL) then
          factor_common_c44_muv = factor_common_c44_muv * c44store(:,:,:,ispec)
        else
          factor_common_c44_muv = factor_common_c44_muv * muvstore(:,:,:,ispec)
        endif

!        do i_memory = 1,5
!          R_memory(i_memory,i_SLS,:,:,:,ispec) = alphaval(i_SLS) * &
!                    R_memory(i_memory,i_SLS,:,:,:,ispec) + &
!                    factor_common_c44_muv * &
!                    (betaval(i_SLS) * epsilondev(i_memory,:,:,:,ispec) + &
!                    gammaval(i_SLS) * epsilondev_loc(i_memory,:,:,:))
!        enddo

        R_xx(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xx(i_SLS,:,:,:,ispec) + factor_common_c44_muv * &
                    (betaval(i_SLS) * epsilondev_xx(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(1,:,:,:))

        R_yy(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_yy(i_SLS,:,:,:,ispec) + factor_common_c44_muv * &
                    (betaval(i_SLS) * epsilondev_yy(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(2,:,:,:))

        R_xy(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xy(i_SLS,:,:,:,ispec) + factor_common_c44_muv * &
                    (betaval(i_SLS) * epsilondev_xy(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(3,:,:,:))

        R_xz(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_xz(i_SLS,:,:,:,ispec) + factor_common_c44_muv * &
                    (betaval(i_SLS) * epsilondev_xz(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(4,:,:,:))

        R_yz(i_SLS,:,:,:,ispec) = alphaval(i_SLS) * R_yz(i_SLS,:,:,:,ispec) + factor_common_c44_muv * &
                    (betaval(i_SLS) * epsilondev_yz(:,:,:,ispec) + gammaval(i_SLS) * epsilondev_loc(5,:,:,:))

      enddo

    endif

! save deviatoric strain for Runge-Kutta scheme
    if(COMPUTE_AND_STORE_STRAIN) then
      !epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            epsilondev_xx(i,j,k,ispec) = epsilondev_loc(1,i,j,k)
            epsilondev_yy(i,j,k,ispec) = epsilondev_loc(2,i,j,k)
            epsilondev_xy(i,j,k,ispec) = epsilondev_loc(3,i,j,k)
            epsilondev_xz(i,j,k,ispec) = epsilondev_loc(4,i,j,k)
            epsilondev_yz(i,j,k,ispec) = epsilondev_loc(5,i,j,k)
          enddo
        enddo
      enddo
    endif

  enddo   ! spectral element loop NSPEC_CRUST_MANTLE

  end subroutine compute_forces_crust_mantle

