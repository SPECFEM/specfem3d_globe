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

  subroutine get_model(myrank,iregion_code,nspec,iproc_xi,iproc_eta, &
    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,rhostore, &
    nspec_ani, &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    xelm,yelm,zelm,shape3D,ispec, &
    iboun,iMPIcut_xi,iMPIcut_eta,rmin,rmax,ichunk,idoubling, &
    beta_montagner,pro_montagner,npar1_montagner, &
    dvs_a,dvs_b,dvp_a,dvp_b,spknt,qq0,qq,abbreviation,code,thlr,velocp,velocs,dens, &
    NPROC_XI,NPROC_ETA, &
    TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,THREE_D, &
    CRUSTAL,ONE_CRUST)

  implicit none

  include "constants.h"

  integer ispec,nspec,ichunk,idoubling,iregion_code,myrank
  integer NPROC_XI,NPROC_ETA
  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,THREE_D,CRUSTAL,ONE_CRUST

  logical iboun(6,nspec)
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

! mantle coefficients for 3D models
  double precision dvs_a(0:NK,0:NS,0:NS),dvs_b(0:NK,0:NS,0:NS)
  double precision dvp_a(0:NK,0:NS,0:NS),dvp_b(0:NK,0:NS,0:NS)
  double precision spknt(NK+1),qq0(NK+1,NK+1),qq(3,NK+1,NK+1)

! for Montagner model
  integer npar1_montagner
  double precision beta_montagner(14,34,37,73),pro_montagner(47)

! for crustal model
  double precision thlr(NKEYS_CRUST,NLAYERS_CRUST),velocp(NKEYS_CRUST,NLAYERS_CRUST)
  double precision velocs(NKEYS_CRUST,NLAYERS_CRUST),dens(NKEYS_CRUST,NLAYERS_CRUST)
  character(len=2) code(NKEYS_CRUST),abbreviation(NCAP_CRUST/2,NCAP_CRUST)

  double precision rmin,rmax

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) kappahstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muhstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec)

  integer nspec_ani

! the 21 coefficients for an anisotropic medium in reduced notation
  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store, &
    c33store,c34store,c35store,c36store, &
    c44store,c45store,c46store,c55store,c56store,c66store

  integer iproc_xi,iproc_eta

  double precision xmesh,ymesh,zmesh

  integer i,j,k,ia
  double precision rho,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision dvp,dvs,drho
  double precision xstore(NGLLX,NGLLY,NGLLZ)
  double precision ystore(NGLLX,NGLLY,NGLLZ)
  double precision zstore(NGLLX,NGLLY,NGLLZ)
  double precision r,r_prem,r_moho,r_dummy,theta,phi
  double precision lat,lon
  double precision vpc,vsc,rhoc,moho

  logical found_crust

  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
       xmesh = ZERO
       ymesh = ZERO
       zmesh = ZERO
       do ia=1,NGNOD
         xmesh = xmesh + shape3D(ia,i,j,k)*xelm(ia)
         ymesh = ymesh + shape3D(ia,i,j,k)*yelm(ia)
         zmesh = zmesh + shape3D(ia,i,j,k)*zelm(ia)
       enddo
       r = dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh)
       xstore(i,j,k) = xmesh
       ystore(i,j,k) = ymesh
       zstore(i,j,k) = zmesh

!      make sure we are within the right shell in PREM to honor discontinuities
!      use small geometrical tolerance
       r_prem = r
       if(r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
       if(r >= rmax*0.999999d0) r_prem = rmax*0.999999d0

!      get the anisotropic PREM parameters
       if(TRANSVERSE_ISOTROPY) then
         call prem_aniso(myrank,r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                                 Qkappa,Qmu,idoubling,CRUSTAL,ONE_CRUST)
       else
         call prem_iso(myrank,r_prem,rho,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,ONE_CRUST,.true.)
         vpv = vp
         vph = vp
         vsv = vs
         vsh = vs
         eta_aniso = 1.d0
       endif

!      get the 3-D model parameters
       if(THREE_D) then
         if(r_prem > RCMB/R_EARTH .and. r_prem < RMOHO/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           dvs = ZERO
           dvp = ZERO
           drho = ZERO
           call three_d_model(r,theta,phi,dvs,dvp,drho, &
                            dvs_a,dvs_b,dvp_a,dvp_b,spknt,qq0,qq)
           vpv=vpv*(1.0d0+dvp)
           vph=vph*(1.0d0+dvp)
           vsv=vsv*(1.0d0+dvs)
           vsh=vsh*(1.0d0+dvs)
           rho=rho*(1.0d0+drho)

! extend 3-D mantle model above the Moho to the surface before adding the crust
         else if(r_prem >= RMOHO/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           dvs = ZERO
           dvp = ZERO
           drho = ZERO
           r_moho = RMOHO/R_EARTH
           call three_d_model(r_moho,theta,phi,dvs,dvp,drho, &
                            dvs_a,dvs_b,dvp_a,dvp_b,spknt,qq0,qq)
           vpv=vpv*(1.0d0+dvp)
           vph=vph*(1.0d0+dvp)
           vsv=vsv*(1.0d0+dvs)
           vsh=vsh*(1.0d0+dvs)
           rho=rho*(1.0d0+drho)

         endif
       endif

       if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) &
           call anisotropic_inner_core_model(r_prem,c11,c33,c12,c13,c44)

       if(ANISOTROPIC_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

! Montagner's model between the Moho and 670 km
         if(r_prem < RMOHO/R_EARTH .and. r_prem > R670/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           call read_montagner_model(r_prem,theta,phi,beta_montagner,pro_montagner,npar1_montagner, &
              rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,myrank)
! extend 3-D mantle model above the Moho to the surface before adding the crust
         elseif(r_prem >= RMOHO/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           r_moho = RMOHO/R_EARTH
           call read_montagner_model(r_moho,theta,phi,beta_montagner,pro_montagner,npar1_montagner, &
              rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,myrank)
! fill the rest of the mantle with the isotropic model
         else
           c11 = rho*vpv*vpv
           c12 = rho*(vpv*vpv-2.*vsv*vsv)
           c13 = c12
           c14 = 0.
           c15 = 0.
           c16 = 0.
           c22 = c11
           c23 = c12
           c24 = 0.
           c25 = 0.
           c26 = 0.
           c33 = c11
           c34 = 0.
           c35 = 0.
           c36 = 0.
           c44 = rho*vsv*vsv
           c45 = 0.
           c46 = 0.
           c55 = c44
           c56 = 0.
           c66 = c44
         endif
       endif

!      get the 3-D crustal model
       if(CRUSTAL) then
         if(r > R_DEEPEST_CRUST) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)

           lat=(PI/2.0d0-theta)*180.0d0/PI
           lon=phi*180.0d0/PI
           if(lon>180.0d0) lon=lon-360.0d0
           call get_crust(lat,lon,r,vpc,vsc,rhoc,moho,found_crust, &
             abbreviation,code,thlr,velocp,velocs,dens)
           if(found_crust) then
             vpv=vpc
             vph=vpc
             vsv=vsc
             vsh=vsc
             rho=rhoc
             if(ANISOTROPIC_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
               c11 = rho*vpv*vpv
               c12 = rho*(vpv*vpv-2.*vsv*vsv)
               c13 = c12
               c14 = 0.
               c15 = 0.
               c16 = 0.
               c22 = c11
               c23 = c12
               c24 = 0.
               c25 = 0.
               c26 = 0.
               c33 = c11
               c34 = 0.
               c35 = 0.
               c36 = 0.
               c44 = rho*vsv*vsv
               c45 = 0.
               c46 = 0.
               c55 = c44
               c56 = 0.
               c66 = c44
             endif
           endif
         endif
       endif

! define elastic parameters in the model

! distinguish whether single or double precision for reals
       if(CUSTOM_REAL == SIZE_REAL) then
         rhostore(i,j,k,ispec) = sngl(rho)
         kappavstore(i,j,k,ispec) = sngl(rho*(vpv*vpv - 4.d0*vsv*vsv/3.d0))
         kappahstore(i,j,k,ispec) = sngl(rho*(vph*vph - 4.d0*vsh*vsh/3.d0))
         muvstore(i,j,k,ispec) = sngl(rho*vsv*vsv)
         muhstore(i,j,k,ispec) = sngl(rho*vsh*vsh)
         eta_anisostore(i,j,k,ispec) = sngl(eta_aniso)

         if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
           c11store(i,j,k,ispec) = sngl(c11)
           c33store(i,j,k,ispec) = sngl(c33)
           c12store(i,j,k,ispec) = sngl(c12)
           c13store(i,j,k,ispec) = sngl(c13)
           c44store(i,j,k,ispec) = sngl(c44)
         endif

         if(ANISOTROPIC_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
           c11store(i,j,k,ispec) = sngl(c11)
           c12store(i,j,k,ispec) = sngl(c12)
           c13store(i,j,k,ispec) = sngl(c13)
           c14store(i,j,k,ispec) = sngl(c14)
           c15store(i,j,k,ispec) = sngl(c15)
           c16store(i,j,k,ispec) = sngl(c16)
           c22store(i,j,k,ispec) = sngl(c22)
           c23store(i,j,k,ispec) = sngl(c23)
           c24store(i,j,k,ispec) = sngl(c24)
           c25store(i,j,k,ispec) = sngl(c25)
           c26store(i,j,k,ispec) = sngl(c26)
           c33store(i,j,k,ispec) = sngl(c33)
           c34store(i,j,k,ispec) = sngl(c34)
           c35store(i,j,k,ispec) = sngl(c35)
           c36store(i,j,k,ispec) = sngl(c36)
           c44store(i,j,k,ispec) = sngl(c44)
           c45store(i,j,k,ispec) = sngl(c45)
           c46store(i,j,k,ispec) = sngl(c46)
           c55store(i,j,k,ispec) = sngl(c55)
           c56store(i,j,k,ispec) = sngl(c56)
           c66store(i,j,k,ispec) = sngl(c66)
         endif

       else
         rhostore(i,j,k,ispec) = rho
         kappavstore(i,j,k,ispec) = rho*(vpv*vpv - 4.d0*vsv*vsv/3.d0)
         kappahstore(i,j,k,ispec) = rho*(vph*vph - 4.d0*vsh*vsh/3.d0)
         muvstore(i,j,k,ispec) = rho*vsv*vsv
         muhstore(i,j,k,ispec) = rho*vsh*vsh
         eta_anisostore(i,j,k,ispec) = eta_aniso

         if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
           c11store(i,j,k,ispec) = c11
           c33store(i,j,k,ispec) = c33
           c12store(i,j,k,ispec) = c12
           c13store(i,j,k,ispec) = c13
           c44store(i,j,k,ispec) = c44
         endif

         if(ANISOTROPIC_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
           c11store(i,j,k,ispec) = c11
           c12store(i,j,k,ispec) = c12
           c13store(i,j,k,ispec) = c13
           c14store(i,j,k,ispec) = c14
           c15store(i,j,k,ispec) = c15
           c16store(i,j,k,ispec) = c16
           c22store(i,j,k,ispec) = c22
           c23store(i,j,k,ispec) = c23
           c24store(i,j,k,ispec) = c24
           c25store(i,j,k,ispec) = c25
           c26store(i,j,k,ispec) = c26
           c33store(i,j,k,ispec) = c33
           c34store(i,j,k,ispec) = c34
           c35store(i,j,k,ispec) = c35
           c36store(i,j,k,ispec) = c36
           c44store(i,j,k,ispec) = c44
           c45store(i,j,k,ispec) = c45
           c46store(i,j,k,ispec) = c46
           c55store(i,j,k,ispec) = c55
           c56store(i,j,k,ispec) = c56
           c66store(i,j,k,ispec) = c66
         endif

       endif

     enddo
   enddo
 enddo

 call get_flags_boundaries(myrank,iregion_code,nspec,iproc_xi,iproc_eta,ispec,xstore,ystore,zstore, &
        iboun,iMPIcut_xi,iMPIcut_eta,ichunk,idoubling,NPROC_XI,NPROC_ETA)

 end subroutine get_model

