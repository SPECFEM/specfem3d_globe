!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine compute_forces_crust_mantle(ell_d80,minus_gravity_table,density_table,minus_deriv_gravity_table, &
          iter,nspec,displ,accel,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgll_cube,wgllwgll_yz_no_i,wgllwgll_xz_no_j,wgllwgll_xy_no_k, &
          kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
          c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
          c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
          c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
          ibool,idoubling,R_memory,epsilondev, &
          one_minus_sum_beta,alphaval,betaval,gammaval,factor_common,update_dof,nglob,index_i,index_k,index_dim)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  integer iter,nspec,nglob

! for ellipticity for d80 attenuation
  real(kind=CUSTOM_REAL) ell_d80,p20,cost

! array with the local to global mapping per slice
  integer, dimension(NSPECMAX_CRUST_MANTLE) :: idoubling

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOBMAX_CRUST_MANTLE) :: displ,accel

! global points in the matching regions
  logical, dimension(NGLOBMAX_CRUST_MANTLE) :: update_dof

! memory variables for attenuation
! memory variables R_ij are stored at the local rather than global level
! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL) one_minus_sum_beta_use,minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION) :: one_minus_sum_beta
  double precision dist
  integer iregion_selected

! for attenuation
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE_ATTENUAT,N_SLS) :: R_memory
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE_ATTENUAT) :: epsilondev
  real(kind=CUSTOM_REAL), dimension(NUM_REGIONS_ATTENUATION,N_SLS) :: alphaval,betaval,gammaval,factor_common
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL) epsilon_trace_over_3

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE) :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: &
    temp1,temp2,temp3,temp1inline,temp2inline,temp3inline
  real(kind=CUSTOM_REAL) temp11,temp12,temp13,temp21,temp22,temp23,temp31,temp32,temp33

! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOBMAX_CRUST_MANTLE) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        kappavstore,muvstore

! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        kappahstore,muhstore,eta_anisostore

! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  integer ispec,iglob
  integer i,j,k,d,ijd,ijk
  integer, dimension(NGLLSQUARE_NDIM) :: index_i,index_k,index_dim

  logical update_elem_iterations

! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  real(kind=CUSTOM_REAL) two_rhovpvsq,two_rhovphsq,two_rhovsvsq,two_rhovshsq
  real(kind=CUSTOM_REAL) four_rhovpvsq,four_rhovphsq,four_rhovsvsq,four_rhovshsq

  real(kind=CUSTOM_REAL) twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) two_eta_aniso,four_eta_aniso,six_eta_aniso

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal,kappavl,kappahl,muvl,muhl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: disploc

! for gravity
  integer int_radius
  real(kind=CUSTOM_REAL) sigma_yx,sigma_zx,sigma_zy
  double precision radius,rho,minus_g,minus_dg
  double precision minus_g_over_radius,minus_dg_plus_g_over_radius
  double precision cos_theta,sin_theta,cos_phi,sin_phi
  double precision cos_theta_sq,sin_theta_sq,cos_phi_sq,sin_phi_sq
  double precision factor,sx_l,sy_l,sz_l,gxl,gyl,gzl
  double precision Hxxl,Hyyl,Hzzl,Hxyl,Hxzl,Hyzl
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: wgllwgll_yz_no_i,wgllwgll_xz_no_j,wgllwgll_xy_no_k

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

! set acceleration to zero where needed
  do i=1,nglob
    if(iter == 1 .or. update_dof(i)) then
      accel(1,i) = 0._CUSTOM_REAL
      accel(2,i) = 0._CUSTOM_REAL
      accel(3,i) = 0._CUSTOM_REAL
    endif
  enddo

!CDIR NOVECTOR
  do ispec = 1,nspec

! define flag to selectively update elements at second iteration
    if(idoubling(ispec) == IFLAG_BOTTOM_MANTLE &
        .or. idoubling(ispec) == IFLAG_BOTTOM_MANTLE_LEV2) then
      update_elem_iterations = .true.
    else
      update_elem_iterations = .false.
    endif

! only matching layers if not first iteration
    if(iter == 1 .or. update_elem_iterations) then

! copy global displacement to local
    do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      disploc(1,ijk,1,1) = displ(1,iglob)
      disploc(2,ijk,1,1) = displ(2,iglob)
      disploc(3,ijk,1,1) = displ(3,iglob)
    enddo

! inlined first matrix product

!CDIR NODEP(temp1,temp2)
  do ijd = 1,NGLLSQUARE_NDIM

    i = index_i(ijd)
    j = i
    k = index_k(ijd)
    d = index_dim(ijd)

!---
  temp1(d,1,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,1) + disploc(d,2,j,k)*hprime_xx(2,1) + &
    disploc(d,3,j,k)*hprime_xx(3,1) + disploc(d,4,j,k)*hprime_xx(4,1) + &
    disploc(d,5,j,k)*hprime_xx(5,1) + disploc(d,6,j,k)*hprime_xx(6,1) + &
    disploc(d,7,j,k)*hprime_xx(7,1) + disploc(d,8,j,k)*hprime_xx(8,1) + &
    disploc(d,9,j,k)*hprime_xx(9,1)

!---
  temp2(d,i,1,k) = &
    disploc(d,i,1,k)*hprime_yy(1,1) + disploc(d,i,2,k)*hprime_yy(2,1) + &
    disploc(d,i,3,k)*hprime_yy(3,1) + disploc(d,i,4,k)*hprime_yy(4,1) + &
    disploc(d,i,5,k)*hprime_yy(5,1) + disploc(d,i,6,k)*hprime_yy(6,1) + &
    disploc(d,i,7,k)*hprime_yy(7,1) + disploc(d,i,8,k)*hprime_yy(8,1) + &
    disploc(d,i,9,k)*hprime_yy(9,1)

!---
  temp3(ijd,1,1,1) = &
    disploc(ijd,1,1,1)*hprime_zz(1,1) + disploc(ijd,1,1,2)*hprime_zz(2,1) + &
    disploc(ijd,1,1,3)*hprime_zz(3,1) + disploc(ijd,1,1,4)*hprime_zz(4,1) + &
    disploc(ijd,1,1,5)*hprime_zz(5,1) + disploc(ijd,1,1,6)*hprime_zz(6,1) + &
    disploc(ijd,1,1,7)*hprime_zz(7,1) + disploc(ijd,1,1,8)*hprime_zz(8,1) + &
    disploc(ijd,1,1,9)*hprime_zz(9,1)

!---
  temp1(d,2,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,2) + disploc(d,2,j,k)*hprime_xx(2,2) + &
    disploc(d,3,j,k)*hprime_xx(3,2) + disploc(d,4,j,k)*hprime_xx(4,2) + &
    disploc(d,5,j,k)*hprime_xx(5,2) + disploc(d,6,j,k)*hprime_xx(6,2) + &
    disploc(d,7,j,k)*hprime_xx(7,2) + disploc(d,8,j,k)*hprime_xx(8,2) + &
    disploc(d,9,j,k)*hprime_xx(9,2)

!---
  temp2(d,i,2,k) = &
    disploc(d,i,1,k)*hprime_yy(1,2) + disploc(d,i,2,k)*hprime_yy(2,2) + &
    disploc(d,i,3,k)*hprime_yy(3,2) + disploc(d,i,4,k)*hprime_yy(4,2) + &
    disploc(d,i,5,k)*hprime_yy(5,2) + disploc(d,i,6,k)*hprime_yy(6,2) + &
    disploc(d,i,7,k)*hprime_yy(7,2) + disploc(d,i,8,k)*hprime_yy(8,2) + &
    disploc(d,i,9,k)*hprime_yy(9,2)

!---
  temp3(ijd,1,1,2) = &
    disploc(ijd,1,1,1)*hprime_zz(1,2) + disploc(ijd,1,1,2)*hprime_zz(2,2) + &
    disploc(ijd,1,1,3)*hprime_zz(3,2) + disploc(ijd,1,1,4)*hprime_zz(4,2) + &
    disploc(ijd,1,1,5)*hprime_zz(5,2) + disploc(ijd,1,1,6)*hprime_zz(6,2) + &
    disploc(ijd,1,1,7)*hprime_zz(7,2) + disploc(ijd,1,1,8)*hprime_zz(8,2) + &
    disploc(ijd,1,1,9)*hprime_zz(9,2)

!---
  temp1(d,3,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,3) + disploc(d,2,j,k)*hprime_xx(2,3) + &
    disploc(d,3,j,k)*hprime_xx(3,3) + disploc(d,4,j,k)*hprime_xx(4,3) + &
    disploc(d,5,j,k)*hprime_xx(5,3) + disploc(d,6,j,k)*hprime_xx(6,3) + &
    disploc(d,7,j,k)*hprime_xx(7,3) + disploc(d,8,j,k)*hprime_xx(8,3) + &
    disploc(d,9,j,k)*hprime_xx(9,3)

!---
  temp2(d,i,3,k) = &
    disploc(d,i,1,k)*hprime_yy(1,3) + disploc(d,i,2,k)*hprime_yy(2,3) + &
    disploc(d,i,3,k)*hprime_yy(3,3) + disploc(d,i,4,k)*hprime_yy(4,3) + &
    disploc(d,i,5,k)*hprime_yy(5,3) + disploc(d,i,6,k)*hprime_yy(6,3) + &
    disploc(d,i,7,k)*hprime_yy(7,3) + disploc(d,i,8,k)*hprime_yy(8,3) + &
    disploc(d,i,9,k)*hprime_yy(9,3)

!---
  temp3(ijd,1,1,3) = &
    disploc(ijd,1,1,1)*hprime_zz(1,3) + disploc(ijd,1,1,2)*hprime_zz(2,3) + &
    disploc(ijd,1,1,3)*hprime_zz(3,3) + disploc(ijd,1,1,4)*hprime_zz(4,3) + &
    disploc(ijd,1,1,5)*hprime_zz(5,3) + disploc(ijd,1,1,6)*hprime_zz(6,3) + &
    disploc(ijd,1,1,7)*hprime_zz(7,3) + disploc(ijd,1,1,8)*hprime_zz(8,3) + &
    disploc(ijd,1,1,9)*hprime_zz(9,3)

!---
  temp1(d,4,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,4) + disploc(d,2,j,k)*hprime_xx(2,4) + &
    disploc(d,3,j,k)*hprime_xx(3,4) + disploc(d,4,j,k)*hprime_xx(4,4) + &
    disploc(d,5,j,k)*hprime_xx(5,4) + disploc(d,6,j,k)*hprime_xx(6,4) + &
    disploc(d,7,j,k)*hprime_xx(7,4) + disploc(d,8,j,k)*hprime_xx(8,4) + &
    disploc(d,9,j,k)*hprime_xx(9,4)

!---
  temp2(d,i,4,k) = &
    disploc(d,i,1,k)*hprime_yy(1,4) + disploc(d,i,2,k)*hprime_yy(2,4) + &
    disploc(d,i,3,k)*hprime_yy(3,4) + disploc(d,i,4,k)*hprime_yy(4,4) + &
    disploc(d,i,5,k)*hprime_yy(5,4) + disploc(d,i,6,k)*hprime_yy(6,4) + &
    disploc(d,i,7,k)*hprime_yy(7,4) + disploc(d,i,8,k)*hprime_yy(8,4) + &
    disploc(d,i,9,k)*hprime_yy(9,4)

!---
  temp3(ijd,1,1,4) = &
    disploc(ijd,1,1,1)*hprime_zz(1,4) + disploc(ijd,1,1,2)*hprime_zz(2,4) + &
    disploc(ijd,1,1,3)*hprime_zz(3,4) + disploc(ijd,1,1,4)*hprime_zz(4,4) + &
    disploc(ijd,1,1,5)*hprime_zz(5,4) + disploc(ijd,1,1,6)*hprime_zz(6,4) + &
    disploc(ijd,1,1,7)*hprime_zz(7,4) + disploc(ijd,1,1,8)*hprime_zz(8,4) + &
    disploc(ijd,1,1,9)*hprime_zz(9,4)

!---
  temp1(d,5,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,5) + disploc(d,2,j,k)*hprime_xx(2,5) + &
    disploc(d,3,j,k)*hprime_xx(3,5) + disploc(d,4,j,k)*hprime_xx(4,5) + &
    disploc(d,5,j,k)*hprime_xx(5,5) + disploc(d,6,j,k)*hprime_xx(6,5) + &
    disploc(d,7,j,k)*hprime_xx(7,5) + disploc(d,8,j,k)*hprime_xx(8,5) + &
    disploc(d,9,j,k)*hprime_xx(9,5)

!---
  temp2(d,i,5,k) = &
    disploc(d,i,1,k)*hprime_yy(1,5) + disploc(d,i,2,k)*hprime_yy(2,5) + &
    disploc(d,i,3,k)*hprime_yy(3,5) + disploc(d,i,4,k)*hprime_yy(4,5) + &
    disploc(d,i,5,k)*hprime_yy(5,5) + disploc(d,i,6,k)*hprime_yy(6,5) + &
    disploc(d,i,7,k)*hprime_yy(7,5) + disploc(d,i,8,k)*hprime_yy(8,5) + &
    disploc(d,i,9,k)*hprime_yy(9,5)

!---
  temp3(ijd,1,1,5) = &
    disploc(ijd,1,1,1)*hprime_zz(1,5) + disploc(ijd,1,1,2)*hprime_zz(2,5) + &
    disploc(ijd,1,1,3)*hprime_zz(3,5) + disploc(ijd,1,1,4)*hprime_zz(4,5) + &
    disploc(ijd,1,1,5)*hprime_zz(5,5) + disploc(ijd,1,1,6)*hprime_zz(6,5) + &
    disploc(ijd,1,1,7)*hprime_zz(7,5) + disploc(ijd,1,1,8)*hprime_zz(8,5) + &
    disploc(ijd,1,1,9)*hprime_zz(9,5)

!---
  temp1(d,6,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,6) + disploc(d,2,j,k)*hprime_xx(2,6) + &
    disploc(d,3,j,k)*hprime_xx(3,6) + disploc(d,4,j,k)*hprime_xx(4,6) + &
    disploc(d,5,j,k)*hprime_xx(5,6) + disploc(d,6,j,k)*hprime_xx(6,6) + &
    disploc(d,7,j,k)*hprime_xx(7,6) + disploc(d,8,j,k)*hprime_xx(8,6) + &
    disploc(d,9,j,k)*hprime_xx(9,6)

!---
  temp2(d,i,6,k) = &
    disploc(d,i,1,k)*hprime_yy(1,6) + disploc(d,i,2,k)*hprime_yy(2,6) + &
    disploc(d,i,3,k)*hprime_yy(3,6) + disploc(d,i,4,k)*hprime_yy(4,6) + &
    disploc(d,i,5,k)*hprime_yy(5,6) + disploc(d,i,6,k)*hprime_yy(6,6) + &
    disploc(d,i,7,k)*hprime_yy(7,6) + disploc(d,i,8,k)*hprime_yy(8,6) + &
    disploc(d,i,9,k)*hprime_yy(9,6)

!---
  temp3(ijd,1,1,6) = &
    disploc(ijd,1,1,1)*hprime_zz(1,6) + disploc(ijd,1,1,2)*hprime_zz(2,6) + &
    disploc(ijd,1,1,3)*hprime_zz(3,6) + disploc(ijd,1,1,4)*hprime_zz(4,6) + &
    disploc(ijd,1,1,5)*hprime_zz(5,6) + disploc(ijd,1,1,6)*hprime_zz(6,6) + &
    disploc(ijd,1,1,7)*hprime_zz(7,6) + disploc(ijd,1,1,8)*hprime_zz(8,6) + &
    disploc(ijd,1,1,9)*hprime_zz(9,6)

!---
  temp1(d,7,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,7) + disploc(d,2,j,k)*hprime_xx(2,7) + &
    disploc(d,3,j,k)*hprime_xx(3,7) + disploc(d,4,j,k)*hprime_xx(4,7) + &
    disploc(d,5,j,k)*hprime_xx(5,7) + disploc(d,6,j,k)*hprime_xx(6,7) + &
    disploc(d,7,j,k)*hprime_xx(7,7) + disploc(d,8,j,k)*hprime_xx(8,7) + &
    disploc(d,9,j,k)*hprime_xx(9,7)

!---
  temp2(d,i,7,k) = &
    disploc(d,i,1,k)*hprime_yy(1,7) + disploc(d,i,2,k)*hprime_yy(2,7) + &
    disploc(d,i,3,k)*hprime_yy(3,7) + disploc(d,i,4,k)*hprime_yy(4,7) + &
    disploc(d,i,5,k)*hprime_yy(5,7) + disploc(d,i,6,k)*hprime_yy(6,7) + &
    disploc(d,i,7,k)*hprime_yy(7,7) + disploc(d,i,8,k)*hprime_yy(8,7) + &
    disploc(d,i,9,k)*hprime_yy(9,7)

!---
  temp3(ijd,1,1,7) = &
    disploc(ijd,1,1,1)*hprime_zz(1,7) + disploc(ijd,1,1,2)*hprime_zz(2,7) + &
    disploc(ijd,1,1,3)*hprime_zz(3,7) + disploc(ijd,1,1,4)*hprime_zz(4,7) + &
    disploc(ijd,1,1,5)*hprime_zz(5,7) + disploc(ijd,1,1,6)*hprime_zz(6,7) + &
    disploc(ijd,1,1,7)*hprime_zz(7,7) + disploc(ijd,1,1,8)*hprime_zz(8,7) + &
    disploc(ijd,1,1,9)*hprime_zz(9,7)

!---
  temp1(d,8,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,8) + disploc(d,2,j,k)*hprime_xx(2,8) + &
    disploc(d,3,j,k)*hprime_xx(3,8) + disploc(d,4,j,k)*hprime_xx(4,8) + &
    disploc(d,5,j,k)*hprime_xx(5,8) + disploc(d,6,j,k)*hprime_xx(6,8) + &
    disploc(d,7,j,k)*hprime_xx(7,8) + disploc(d,8,j,k)*hprime_xx(8,8) + &
    disploc(d,9,j,k)*hprime_xx(9,8)

!---
  temp2(d,i,8,k) = &
    disploc(d,i,1,k)*hprime_yy(1,8) + disploc(d,i,2,k)*hprime_yy(2,8) + &
    disploc(d,i,3,k)*hprime_yy(3,8) + disploc(d,i,4,k)*hprime_yy(4,8) + &
    disploc(d,i,5,k)*hprime_yy(5,8) + disploc(d,i,6,k)*hprime_yy(6,8) + &
    disploc(d,i,7,k)*hprime_yy(7,8) + disploc(d,i,8,k)*hprime_yy(8,8) + &
    disploc(d,i,9,k)*hprime_yy(9,8)

!---
  temp3(ijd,1,1,8) = &
    disploc(ijd,1,1,1)*hprime_zz(1,8) + disploc(ijd,1,1,2)*hprime_zz(2,8) + &
    disploc(ijd,1,1,3)*hprime_zz(3,8) + disploc(ijd,1,1,4)*hprime_zz(4,8) + &
    disploc(ijd,1,1,5)*hprime_zz(5,8) + disploc(ijd,1,1,6)*hprime_zz(6,8) + &
    disploc(ijd,1,1,7)*hprime_zz(7,8) + disploc(ijd,1,1,8)*hprime_zz(8,8) + &
    disploc(ijd,1,1,9)*hprime_zz(9,8)

!---
  temp1(d,9,j,k) = &
    disploc(d,1,j,k)*hprime_xx(1,9) + disploc(d,2,j,k)*hprime_xx(2,9) + &
    disploc(d,3,j,k)*hprime_xx(3,9) + disploc(d,4,j,k)*hprime_xx(4,9) + &
    disploc(d,5,j,k)*hprime_xx(5,9) + disploc(d,6,j,k)*hprime_xx(6,9) + &
    disploc(d,7,j,k)*hprime_xx(7,9) + disploc(d,8,j,k)*hprime_xx(8,9) + &
    disploc(d,9,j,k)*hprime_xx(9,9)

!---
  temp2(d,i,9,k) = &
    disploc(d,i,1,k)*hprime_yy(1,9) + disploc(d,i,2,k)*hprime_yy(2,9) + &
    disploc(d,i,3,k)*hprime_yy(3,9) + disploc(d,i,4,k)*hprime_yy(4,9) + &
    disploc(d,i,5,k)*hprime_yy(5,9) + disploc(d,i,6,k)*hprime_yy(6,9) + &
    disploc(d,i,7,k)*hprime_yy(7,9) + disploc(d,i,8,k)*hprime_yy(8,9) + &
    disploc(d,i,9,k)*hprime_yy(9,9)

!---
  temp3(ijd,1,1,9) = &
    disploc(ijd,1,1,1)*hprime_zz(1,9) + disploc(ijd,1,1,2)*hprime_zz(2,9) + &
    disploc(ijd,1,1,3)*hprime_zz(3,9) + disploc(ijd,1,1,4)*hprime_zz(4,9) + &
    disploc(ijd,1,1,5)*hprime_zz(5,9) + disploc(ijd,1,1,6)*hprime_zz(6,9) + &
    disploc(ijd,1,1,7)*hprime_zz(7,9) + disploc(ijd,1,1,8)*hprime_zz(8,9) + &
    disploc(ijd,1,1,9)*hprime_zz(9,9)

  enddo

! fused the three loops for inlined version
    do ijk=1,NGLLCUBE

!         get derivatives of ux, uy and uz with respect to x, y and z

          xixl = xix(ijk,1,1,ispec)
          xiyl = xiy(ijk,1,1,ispec)
          xizl = xiz(ijk,1,1,ispec)
          etaxl = etax(ijk,1,1,ispec)
          etayl = etay(ijk,1,1,ispec)
          etazl = etaz(ijk,1,1,ispec)
          gammaxl = gammax(ijk,1,1,ispec)
          gammayl = gammay(ijk,1,1,ispec)
          gammazl = gammaz(ijk,1,1,ispec)
          jacobianl = jacobian(ijk,1,1,ispec)

          temp11 = temp1(1,ijk,1,1)
          temp12 = temp1(2,ijk,1,1)
          temp13 = temp1(3,ijk,1,1)

          temp21 = temp2(1,ijk,1,1)
          temp22 = temp2(2,ijk,1,1)
          temp23 = temp2(3,ijk,1,1)

          temp31 = temp3(1,ijk,1,1)
          temp32 = temp3(2,ijk,1,1)
          temp33 = temp3(3,ijk,1,1)

          duxdxl = xixl*temp11 + etaxl*temp21 + gammaxl*temp31
          duxdyl = xiyl*temp11 + etayl*temp21 + gammayl*temp31
          duxdzl = xizl*temp11 + etazl*temp21 + gammazl*temp31

          duydxl = xixl*temp12 + etaxl*temp22 + gammaxl*temp32
          duydyl = xiyl*temp12 + etayl*temp22 + gammayl*temp32
          duydzl = xizl*temp12 + etazl*temp22 + gammazl*temp32

          duzdxl = xixl*temp13 + etaxl*temp23 + gammaxl*temp33
          duzdyl = xiyl*temp13 + etayl*temp23 + gammayl*temp33
          duzdzl = xizl*temp13 + etazl*temp23 + gammazl*temp33

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! precompute terms for attenuation if needed
  if(ATTENUATION_VAL) then

! compute deviatoric strain
    epsilon_trace_over_3 = ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilondev_loc(1,ijk,1,1) = duxdxl - epsilon_trace_over_3
    epsilondev_loc(2,ijk,1,1) = duydyl - epsilon_trace_over_3
    epsilondev_loc(3,ijk,1,1) = 0.5 * duxdyl_plus_duydxl
    epsilondev_loc(4,ijk,1,1) = 0.5 * duzdxl_plus_duxdzl
    epsilondev_loc(5,ijk,1,1) = 0.5 * duzdyl_plus_duydzl

! distinguish regions in the mantle, including case of the d80

  if(idoubling(ispec) == IFLAG_DOUBLING_670 .or. &
     idoubling(ispec) == IFLAG_MANTLE_NORMAL .or. &
     idoubling(ispec) == IFLAG_BOTTOM_MANTLE_LEV2 .or. &
     idoubling(ispec) == IFLAG_BOTTOM_MANTLE) then

    iregion_selected = IREGION_ATTENUATION_CMB_670

  else if(idoubling(ispec) == IFLAG_670_220) then

    iregion_selected = IREGION_ATTENUATION_670_220

  else

! particular case of d80 which is not honored by the mesh
! xstore contains the radius
    iglob = ibool(ijk,1,1,ispec)
    dist = xstore(iglob)

! map ellipticity back for d80 detection
! ystore contains theta
    if(ELLIPTICITY_VAL) then
      theta = ystore(iglob)
      cost = cos(theta)
      p20 = 0.5*(3.*cost*cost-1.)
      dist = dist*(1.+(2./3.)*ell_d80*p20)
    endif

    if(dist > R80/R_EARTH) then
      iregion_selected = IREGION_ATTENUATION_80_SURFACE
    else
      iregion_selected = IREGION_ATTENUATION_220_80
    endif

   endif

    one_minus_sum_beta_use = one_minus_sum_beta(iregion_selected)
    minus_sum_beta =  one_minus_sum_beta_use - 1.

  endif

!
! compute either isotropic or anisotropic elements
!

  if(ANISOTROPIC_MANTLE_VAL) then

    c11 = c11store(ijk,1,1,ispec)
    c12 = c12store(ijk,1,1,ispec)
    c13 = c13store(ijk,1,1,ispec)
    c14 = c14store(ijk,1,1,ispec)
    c15 = c15store(ijk,1,1,ispec)
    c16 = c16store(ijk,1,1,ispec)
    c22 = c22store(ijk,1,1,ispec)
    c23 = c23store(ijk,1,1,ispec)
    c24 = c24store(ijk,1,1,ispec)
    c25 = c25store(ijk,1,1,ispec)
    c26 = c26store(ijk,1,1,ispec)
    c33 = c33store(ijk,1,1,ispec)
    c34 = c34store(ijk,1,1,ispec)
    c35 = c35store(ijk,1,1,ispec)
    c36 = c36store(ijk,1,1,ispec)
    c44 = c44store(ijk,1,1,ispec)
    c45 = c45store(ijk,1,1,ispec)
    c46 = c46store(ijk,1,1,ispec)
    c55 = c55store(ijk,1,1,ispec)
    c56 = c56store(ijk,1,1,ispec)
    c66 = c66store(ijk,1,1,ispec)

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
  if(.not. (TRANSVERSE_ISOTROPY_VAL .and. idoubling(ispec) == IFLAG_220_MOHO)) then

! layer with no transverse isotropy, use kappav and muv
          kappal = kappavstore(ijk,1,1,ispec)
          mul = muvstore(ijk,1,1,ispec)

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
      kappavl = kappavstore(ijk,1,1,ispec)
      muvl = muvstore(ijk,1,1,ispec)

      kappahl = kappahstore(ijk,1,1,ispec)
      muhl = muhstore(ijk,1,1,ispec)

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

  eta_aniso = eta_anisostore(ijk,1,1,ispec)  !!! that is  F / (A - 2 L)

! use mesh coordinates to get theta and phi
! ystore and zstore contain theta and phi

  iglob = ibool(ijk,1,1,ispec)
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

      two_rhovpvsq = 2.*rhovpvsq
      two_rhovphsq = 2.*rhovphsq
      two_rhovsvsq = 2.*rhovsvsq
      two_rhovshsq = 2.*rhovshsq

      four_rhovpvsq = 4.*rhovpvsq
      four_rhovphsq = 4.*rhovphsq
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
  if(ATTENUATION_VAL) then

    R_xx_val = R_memory(1,ijk,1,1,ispec,1)
    R_yy_val = R_memory(2,ijk,1,1,ispec,1)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_memory(3,ijk,1,1,ispec,1)
    sigma_xz = sigma_xz - R_memory(4,ijk,1,1,ispec,1)
    sigma_yz = sigma_yz - R_memory(5,ijk,1,1,ispec,1)

    R_xx_val = R_memory(1,ijk,1,1,ispec,2)
    R_yy_val = R_memory(2,ijk,1,1,ispec,2)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_memory(3,ijk,1,1,ispec,2)
    sigma_xz = sigma_xz - R_memory(4,ijk,1,1,ispec,2)
    sigma_yz = sigma_yz - R_memory(5,ijk,1,1,ispec,2)

    R_xx_val = R_memory(1,ijk,1,1,ispec,3)
    R_yy_val = R_memory(2,ijk,1,1,ispec,3)
    sigma_xx = sigma_xx - R_xx_val
    sigma_yy = sigma_yy - R_yy_val
    sigma_zz = sigma_zz + R_xx_val + R_yy_val
    sigma_xy = sigma_xy - R_memory(3,ijk,1,1,ispec,3)
    sigma_xz = sigma_xz - R_memory(4,ijk,1,1,ispec,3)
    sigma_yz = sigma_yz - R_memory(5,ijk,1,1,ispec,3)

  endif

! define symmetric components of sigma for gravity
  sigma_yx = sigma_xy
  sigma_zx = sigma_xz
  sigma_zy = sigma_yz

! compute non-symmetric terms for gravity
  if(GRAVITY_VAL) then

! use mesh coordinates to get theta and phi
! x y and z contain r theta and phi

  iglob = ibool(ijk,1,1,ispec)
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

! distinguish whether single or double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then

! get displacement and multiply by density to compute G tensor
    sx_l = rho * dble(disploc(1,ijk,1,1))
    sy_l = rho * dble(disploc(2,ijk,1,1))
    sz_l = rho * dble(disploc(3,ijk,1,1))

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
    factor = dble(jacobianl) * wgll_cube(ijk,1,1)
    rho_s_H(1,ijk,1,1) = sngl(factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl))
    rho_s_H(2,ijk,1,1) = sngl(factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl))
    rho_s_H(3,ijk,1,1) = sngl(factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl))

  else

! get displacement and multiply by density to compute G tensor
    sx_l = rho * disploc(1,ijk,1,1)
    sy_l = rho * disploc(2,ijk,1,1)
    sz_l = rho * disploc(3,ijk,1,1)

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
    factor = jacobianl * wgll_cube(ijk,1,1)
    rho_s_H(1,ijk,1,1) = factor * (sx_l * Hxxl + sy_l * Hxyl + sz_l * Hxzl)
    rho_s_H(2,ijk,1,1) = factor * (sx_l * Hxyl + sy_l * Hyyl + sz_l * Hyzl)
    rho_s_H(3,ijk,1,1) = factor * (sx_l * Hxzl + sy_l * Hyzl + sz_l * Hzzl)

  endif

  endif  ! end of section with gravity terms

! form dot product with test vector, non-symmetric form
      temp1(1,ijk,1,1) = jacobianl * (sigma_xx*xixl + sigma_yx*xiyl + sigma_zx*xizl)
      temp1(2,ijk,1,1) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_zy*xizl)
      temp1(3,ijk,1,1) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

      temp2(1,ijk,1,1) = jacobianl * (sigma_xx*etaxl + sigma_yx*etayl + sigma_zx*etazl)
      temp2(2,ijk,1,1) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_zy*etazl)
      temp2(3,ijk,1,1) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

      temp3(1,ijk,1,1) = jacobianl * (sigma_xx*gammaxl + sigma_yx*gammayl + sigma_zx*gammazl)
      temp3(2,ijk,1,1) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_zy*gammazl)
      temp3(3,ijk,1,1) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)

      enddo

! inlined second matrix product

!CDIR NODEP(temp1inline,temp2inline)
  do ijd = 1,NGLLSQUARE_NDIM

    i = index_i(ijd)
    j = i
    k = index_k(ijd)
    d = index_dim(ijd)

!---
  temp1inline(d,1,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(1,1) + temp1(d,2,j,k)*hprimewgll_xx(1,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(1,3) + temp1(d,4,j,k)*hprimewgll_xx(1,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(1,5) + temp1(d,6,j,k)*hprimewgll_xx(1,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(1,7) + temp1(d,8,j,k)*hprimewgll_xx(1,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(1,9)

!---
  temp2inline(d,i,1,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(1,1) + temp2(d,i,2,k)*hprimewgll_yy(1,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(1,3) + temp2(d,i,4,k)*hprimewgll_yy(1,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(1,5) + temp2(d,i,6,k)*hprimewgll_yy(1,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(1,7) + temp2(d,i,8,k)*hprimewgll_yy(1,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(1,9)

!---
  temp3inline(ijd,1,1,1) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(1,1) + temp3(ijd,1,1,2)*hprimewgll_zz(1,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(1,3) + temp3(ijd,1,1,4)*hprimewgll_zz(1,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(1,5) + temp3(ijd,1,1,6)*hprimewgll_zz(1,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(1,7) + temp3(ijd,1,1,8)*hprimewgll_zz(1,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(1,9)

!---
  temp1inline(d,2,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(2,1) + temp1(d,2,j,k)*hprimewgll_xx(2,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(2,3) + temp1(d,4,j,k)*hprimewgll_xx(2,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(2,5) + temp1(d,6,j,k)*hprimewgll_xx(2,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(2,7) + temp1(d,8,j,k)*hprimewgll_xx(2,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(2,9)

!---
  temp2inline(d,i,2,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(2,1) + temp2(d,i,2,k)*hprimewgll_yy(2,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(2,3) + temp2(d,i,4,k)*hprimewgll_yy(2,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(2,5) + temp2(d,i,6,k)*hprimewgll_yy(2,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(2,7) + temp2(d,i,8,k)*hprimewgll_yy(2,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(2,9)

!---
  temp3inline(ijd,1,1,2) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(2,1) + temp3(ijd,1,1,2)*hprimewgll_zz(2,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(2,3) + temp3(ijd,1,1,4)*hprimewgll_zz(2,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(2,5) + temp3(ijd,1,1,6)*hprimewgll_zz(2,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(2,7) + temp3(ijd,1,1,8)*hprimewgll_zz(2,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(2,9)

!---
  temp1inline(d,3,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(3,1) + temp1(d,2,j,k)*hprimewgll_xx(3,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(3,3) + temp1(d,4,j,k)*hprimewgll_xx(3,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(3,5) + temp1(d,6,j,k)*hprimewgll_xx(3,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(3,7) + temp1(d,8,j,k)*hprimewgll_xx(3,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(3,9)

!---
  temp2inline(d,i,3,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(3,1) + temp2(d,i,2,k)*hprimewgll_yy(3,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(3,3) + temp2(d,i,4,k)*hprimewgll_yy(3,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(3,5) + temp2(d,i,6,k)*hprimewgll_yy(3,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(3,7) + temp2(d,i,8,k)*hprimewgll_yy(3,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(3,9)

!---
  temp3inline(ijd,1,1,3) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(3,1) + temp3(ijd,1,1,2)*hprimewgll_zz(3,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(3,3) + temp3(ijd,1,1,4)*hprimewgll_zz(3,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(3,5) + temp3(ijd,1,1,6)*hprimewgll_zz(3,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(3,7) + temp3(ijd,1,1,8)*hprimewgll_zz(3,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(3,9)

!---
  temp1inline(d,4,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(4,1) + temp1(d,2,j,k)*hprimewgll_xx(4,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(4,3) + temp1(d,4,j,k)*hprimewgll_xx(4,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(4,5) + temp1(d,6,j,k)*hprimewgll_xx(4,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(4,7) + temp1(d,8,j,k)*hprimewgll_xx(4,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(4,9)

!---
  temp2inline(d,i,4,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(4,1) + temp2(d,i,2,k)*hprimewgll_yy(4,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(4,3) + temp2(d,i,4,k)*hprimewgll_yy(4,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(4,5) + temp2(d,i,6,k)*hprimewgll_yy(4,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(4,7) + temp2(d,i,8,k)*hprimewgll_yy(4,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(4,9)

!---
  temp3inline(ijd,1,1,4) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(4,1) + temp3(ijd,1,1,2)*hprimewgll_zz(4,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(4,3) + temp3(ijd,1,1,4)*hprimewgll_zz(4,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(4,5) + temp3(ijd,1,1,6)*hprimewgll_zz(4,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(4,7) + temp3(ijd,1,1,8)*hprimewgll_zz(4,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(4,9)

!---
  temp1inline(d,5,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(5,1) + temp1(d,2,j,k)*hprimewgll_xx(5,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(5,3) + temp1(d,4,j,k)*hprimewgll_xx(5,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(5,5) + temp1(d,6,j,k)*hprimewgll_xx(5,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(5,7) + temp1(d,8,j,k)*hprimewgll_xx(5,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(5,9)

!---
  temp2inline(d,i,5,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(5,1) + temp2(d,i,2,k)*hprimewgll_yy(5,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(5,3) + temp2(d,i,4,k)*hprimewgll_yy(5,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(5,5) + temp2(d,i,6,k)*hprimewgll_yy(5,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(5,7) + temp2(d,i,8,k)*hprimewgll_yy(5,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(5,9)

!---
  temp3inline(ijd,1,1,5) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(5,1) + temp3(ijd,1,1,2)*hprimewgll_zz(5,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(5,3) + temp3(ijd,1,1,4)*hprimewgll_zz(5,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(5,5) + temp3(ijd,1,1,6)*hprimewgll_zz(5,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(5,7) + temp3(ijd,1,1,8)*hprimewgll_zz(5,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(5,9)

!---
  temp1inline(d,6,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(6,1) + temp1(d,2,j,k)*hprimewgll_xx(6,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(6,3) + temp1(d,4,j,k)*hprimewgll_xx(6,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(6,5) + temp1(d,6,j,k)*hprimewgll_xx(6,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(6,7) + temp1(d,8,j,k)*hprimewgll_xx(6,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(6,9)

!---
  temp2inline(d,i,6,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(6,1) + temp2(d,i,2,k)*hprimewgll_yy(6,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(6,3) + temp2(d,i,4,k)*hprimewgll_yy(6,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(6,5) + temp2(d,i,6,k)*hprimewgll_yy(6,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(6,7) + temp2(d,i,8,k)*hprimewgll_yy(6,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(6,9)

!---
  temp3inline(ijd,1,1,6) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(6,1) + temp3(ijd,1,1,2)*hprimewgll_zz(6,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(6,3) + temp3(ijd,1,1,4)*hprimewgll_zz(6,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(6,5) + temp3(ijd,1,1,6)*hprimewgll_zz(6,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(6,7) + temp3(ijd,1,1,8)*hprimewgll_zz(6,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(6,9)

!---
  temp1inline(d,7,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(7,1) + temp1(d,2,j,k)*hprimewgll_xx(7,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(7,3) + temp1(d,4,j,k)*hprimewgll_xx(7,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(7,5) + temp1(d,6,j,k)*hprimewgll_xx(7,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(7,7) + temp1(d,8,j,k)*hprimewgll_xx(7,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(7,9)

!---
  temp2inline(d,i,7,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(7,1) + temp2(d,i,2,k)*hprimewgll_yy(7,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(7,3) + temp2(d,i,4,k)*hprimewgll_yy(7,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(7,5) + temp2(d,i,6,k)*hprimewgll_yy(7,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(7,7) + temp2(d,i,8,k)*hprimewgll_yy(7,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(7,9)

!---
  temp3inline(ijd,1,1,7) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(7,1) + temp3(ijd,1,1,2)*hprimewgll_zz(7,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(7,3) + temp3(ijd,1,1,4)*hprimewgll_zz(7,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(7,5) + temp3(ijd,1,1,6)*hprimewgll_zz(7,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(7,7) + temp3(ijd,1,1,8)*hprimewgll_zz(7,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(7,9)

!---
  temp1inline(d,8,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(8,1) + temp1(d,2,j,k)*hprimewgll_xx(8,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(8,3) + temp1(d,4,j,k)*hprimewgll_xx(8,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(8,5) + temp1(d,6,j,k)*hprimewgll_xx(8,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(8,7) + temp1(d,8,j,k)*hprimewgll_xx(8,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(8,9)

!---
  temp2inline(d,i,8,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(8,1) + temp2(d,i,2,k)*hprimewgll_yy(8,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(8,3) + temp2(d,i,4,k)*hprimewgll_yy(8,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(8,5) + temp2(d,i,6,k)*hprimewgll_yy(8,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(8,7) + temp2(d,i,8,k)*hprimewgll_yy(8,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(8,9)

!---
  temp3inline(ijd,1,1,8) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(8,1) + temp3(ijd,1,1,2)*hprimewgll_zz(8,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(8,3) + temp3(ijd,1,1,4)*hprimewgll_zz(8,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(8,5) + temp3(ijd,1,1,6)*hprimewgll_zz(8,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(8,7) + temp3(ijd,1,1,8)*hprimewgll_zz(8,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(8,9)

!---
  temp1inline(d,9,j,k) = &
    temp1(d,1,j,k)*hprimewgll_xx(9,1) + temp1(d,2,j,k)*hprimewgll_xx(9,2) + &
    temp1(d,3,j,k)*hprimewgll_xx(9,3) + temp1(d,4,j,k)*hprimewgll_xx(9,4) + &
    temp1(d,5,j,k)*hprimewgll_xx(9,5) + temp1(d,6,j,k)*hprimewgll_xx(9,6) + &
    temp1(d,7,j,k)*hprimewgll_xx(9,7) + temp1(d,8,j,k)*hprimewgll_xx(9,8) + &
    temp1(d,9,j,k)*hprimewgll_xx(9,9)

!---
  temp2inline(d,i,9,k) = &
    temp2(d,i,1,k)*hprimewgll_yy(9,1) + temp2(d,i,2,k)*hprimewgll_yy(9,2) + &
    temp2(d,i,3,k)*hprimewgll_yy(9,3) + temp2(d,i,4,k)*hprimewgll_yy(9,4) + &
    temp2(d,i,5,k)*hprimewgll_yy(9,5) + temp2(d,i,6,k)*hprimewgll_yy(9,6) + &
    temp2(d,i,7,k)*hprimewgll_yy(9,7) + temp2(d,i,8,k)*hprimewgll_yy(9,8) + &
    temp2(d,i,9,k)*hprimewgll_yy(9,9)

!---
  temp3inline(ijd,1,1,9) = &
   temp3(ijd,1,1,1)*hprimewgll_zz(9,1) + temp3(ijd,1,1,2)*hprimewgll_zz(9,2) + &
   temp3(ijd,1,1,3)*hprimewgll_zz(9,3) + temp3(ijd,1,1,4)*hprimewgll_zz(9,4) + &
   temp3(ijd,1,1,5)*hprimewgll_zz(9,5) + temp3(ijd,1,1,6)*hprimewgll_zz(9,6) + &
   temp3(ijd,1,1,7)*hprimewgll_zz(9,7) + temp3(ijd,1,1,8)*hprimewgll_zz(9,8) + &
   temp3(ijd,1,1,9)*hprimewgll_zz(9,9)

  enddo

  if(GRAVITY_VAL) then
    do ijk=1,NGLLCUBE_NDIM
      sum_terms(ijk,1,1,1) = &
           - (wgllwgll_yz_no_i(ijk,1,1,1)*temp1inline(ijk,1,1,1) + &
              wgllwgll_xz_no_j(ijk,1,1,1)*temp2inline(ijk,1,1,1) + &
              wgllwgll_xy_no_k(ijk,1,1,1)*temp3inline(ijk,1,1,1)) + &
              rho_s_H(ijk,1,1,1)
    enddo
  else
    do ijk=1,NGLLCUBE_NDIM
      sum_terms(ijk,1,1,1) = &
           - (wgllwgll_yz_no_i(ijk,1,1,1)*temp1inline(ijk,1,1,1) + &
              wgllwgll_xz_no_j(ijk,1,1,1)*temp2inline(ijk,1,1,1) + &
              wgllwgll_xy_no_k(ijk,1,1,1)*temp3inline(ijk,1,1,1))
    enddo
  endif

! sum contributions from each element to the global mesh and add gravity terms
!CDIR NODEP(accel)
    do ijk=1,NGLLCUBE
      iglob = ibool(ijk,1,1,ispec)
      if(iter == 1 .or. update_dof(iglob)) then
        accel(1,iglob) = accel(1,iglob) + sum_terms(1,ijk,1,1)
        accel(2,iglob) = accel(2,iglob) + sum_terms(2,ijk,1,1)
        accel(3,iglob) = accel(3,iglob) + sum_terms(3,ijk,1,1)
      endif
    enddo

! update memory variables based upon the Runge-Kutta scheme
! update at first iteration everywhere except in the matching layer at the CMB
! convention for attenuation
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace

  if(ATTENUATION_VAL) then

    if(iter == 2 .or. (iter == 1 .and. .not. update_elem_iterations)) then

! use Runge-Kutta scheme to march in time

! get coefficients for that standard linear solid
! IMPROVE we use mu_v here even if there is some anisotropy
! IMPROVE we should probably use an average value instead
  if(ANISOTROPIC_MANTLE_VAL) then
    do ijk = 1,NGLLCUBE

!---
    R_memory(1,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(1,ijk,1,1,ispec,1) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(1,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(2,ijk,1,1,ispec,1) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(2,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(3,ijk,1,1,ispec,1) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(3,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(4,ijk,1,1,ispec,1) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(4,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(5,ijk,1,1,ispec,1) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(5,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(5,ijk,1,1))

!---
    R_memory(1,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(1,ijk,1,1,ispec,2) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(1,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(2,ijk,1,1,ispec,2) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(2,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(3,ijk,1,1,ispec,2) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(3,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(4,ijk,1,1,ispec,2) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(4,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(5,ijk,1,1,ispec,2) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(5,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(5,ijk,1,1))

!---
    R_memory(1,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(1,ijk,1,1,ispec,3) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(1,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(2,ijk,1,1,ispec,3) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(2,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(3,ijk,1,1,ispec,3) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(3,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(4,ijk,1,1,ispec,3) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(4,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(5,ijk,1,1,ispec,3) + c44store(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(5,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(5,ijk,1,1))

    enddo
  else

    do ijk = 1,NGLLCUBE

!---
    R_memory(1,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(1,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(1,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(2,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(2,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(3,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(3,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(4,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(4,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,1) = alphaval(iregion_selected,1) * &
      R_memory(5,ijk,1,1,ispec,1) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,1) * &
      (betaval(iregion_selected,1) * epsilondev(5,ijk,1,1,ispec) + &
       gammaval(iregion_selected,1) * epsilondev_loc(5,ijk,1,1))

!---
    R_memory(1,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(1,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(1,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(2,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(2,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(3,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(3,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(4,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(4,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,2) = alphaval(iregion_selected,2) * &
      R_memory(5,ijk,1,1,ispec,2) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,2) * &
      (betaval(iregion_selected,2) * epsilondev(5,ijk,1,1,ispec) + &
       gammaval(iregion_selected,2) * epsilondev_loc(5,ijk,1,1))

!---
    R_memory(1,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(1,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(1,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(1,ijk,1,1))

    R_memory(2,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(2,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(2,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(2,ijk,1,1))

    R_memory(3,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(3,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(3,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(3,ijk,1,1))

    R_memory(4,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(4,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(4,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(4,ijk,1,1))

    R_memory(5,ijk,1,1,ispec,3) = alphaval(iregion_selected,3) * &
      R_memory(5,ijk,1,1,ispec,3) + muvstore(ijk,1,1,ispec) * &
      factor_common(iregion_selected,3) * &
      (betaval(iregion_selected,3) * epsilondev(5,ijk,1,1,ispec) + &
       gammaval(iregion_selected,3) * epsilondev_loc(5,ijk,1,1))

    enddo
  endif

  endif

! save deviatoric strain for Runge-Kutta scheme
    do ijk = 1,5*NGLLCUBE
      epsilondev(ijk,1,1,1,ispec) = epsilondev_loc(ijk,1,1,1)
    enddo

  endif

  endif   ! end test only matching layers if not first iteration

  enddo   ! spectral element loop

  end subroutine compute_forces_crust_mantle

