!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
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
    rho_vp,rho_vs,nspec_stacey, &
    NPROC_XI,NPROC_ETA, &
    TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE, &
    CRUSTAL,ONE_CRUST, &
    crustal_model,mantle_model,aniso_mantle_model, &
    aniso_inner_core_model,rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD,&
    attenuation_model,ATTENUATION,ATTENUATION_3D,tau_s,tau_e_store,Qmu_store,T_c_source,vx,vy,vz,vnspec, &
    NCHUNKS,INFLATE_CENTRAL_CUBE,ABSORBING_CONDITIONS,IASPEI, &
    R_CENTRAL_CUBE,RCMB,RICB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN)

  implicit none

  include "constants.h"

  external mantle_model,crustal_model,aniso_mantle_model, &
       aniso_inner_core_model,attenuation_model

  integer ispec,nspec,ichunk,idoubling,iregion_code,myrank,nspec_stacey
  integer NPROC_XI,NPROC_ETA,NCHUNKS

  logical ATTENUATION,ATTENUATION_3D,ABSORBING_CONDITIONS,INFLATE_CENTRAL_CUBE,IASPEI
  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST

  logical iboun(6,nspec)
  logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision rmin,rmax,R_CENTRAL_CUBE,RCMB,RICB,R670,RMOHO, &
    RTOPDDOUBLEPRIME,R600,R220,R771,R400,R80,RMIDDLE_CRUST,ROCEAN

  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) kappahstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) muhstore(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec)

  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey),rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey)

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

  ! Attenuation values
  integer vx, vy, vz, vnspec
  double precision, dimension(N_SLS)                     :: tau_s, tau_e
  double precision, dimension(vx, vy, vz, vnspec)        :: Qmu_store
  double precision, dimension(N_SLS, vx, vy, vz, vnspec) :: tau_e_store
  double precision  T_c_source

  logical found_crust

  double precision ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

! rotation matrix from Euler angles
  double precision rotation_matrix(3,3)

  tau_s(:)   = 0.0d0
  tau_e(:)   = 0.0d0
  T_c_source = 0.0d0

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
           Qkappa,Qmu,idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
           R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
       else
         call prem_iso(myrank,r_prem,rho,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
           ONE_CRUST,.true.,IASPEI,RICB,RCMB,RTOPDDOUBLEPRIME, &
           R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
         vpv = vp
         vph = vp
         vsv = vs
         vsh = vs
         eta_aniso = 1.d0
       endif

!      get the 3-D model parameters
       if(ISOTROPIC_3D_MANTLE) then
         if(r_prem > RCMB/R_EARTH .and. r_prem < RMOHO/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           dvs = ZERO
           dvp = ZERO
           drho = ZERO
           call mantle_model(r,theta,phi,dvs,dvp,drho)
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
           call mantle_model(r_moho,theta,phi,dvs,dvp,drho)
           vpv=vpv*(1.0d0+dvp)
           vph=vph*(1.0d0+dvp)
           vsv=vsv*(1.0d0+dvs)
           vsh=vsh*(1.0d0+dvs)
           rho=rho*(1.0d0+drho)

         endif
       endif

       if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) &
           call aniso_inner_core_model(r_prem,c11,c33,c12,c13,c44,IASPEI)

       if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

! Montagner's model between the Moho and 670 km
         if(r_prem < RMOHO/R_EARTH .and. r_prem > R670/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           call aniso_mantle_model(myrank,r_prem,theta,phi, &
              rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
! extend 3-D mantle model above the Moho to the surface before adding the crust
         elseif(r_prem >= RMOHO/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           r_moho = RMOHO/R_EARTH
           call aniso_mantle_model(myrank,r_moho,theta,phi, &
              rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)
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

       if(ATTENUATION .and. ATTENUATION_3D) then
         call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
         call reduce(theta,phi)
         lat=(PI/2.0d0-theta)*180.0d0/PI
         lon=phi*180.0d0/PI
         if(lon > 180.0d0) lon = lon - 360.0d0
         call attenuation_model(myrank, lat, lon, r, Qmu, tau_s, tau_e, T_c_source,RICB,RCMB, &
      RTOPDDOUBLEPRIME,R600,R670,R220,R771,R400,R80)
       endif

!      get the 3-D crustal model
       if(CRUSTAL) then
         if(r > R_DEEPEST_CRUST) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)

           lat=(PI/2.0d0-theta)*180.0d0/PI
           lon=phi*180.0d0/PI
           if(lon>180.0d0) lon=lon-360.0d0
           call crustal_model(lat,lon,r,vpc,vsc,rhoc,moho,found_crust)
           if (found_crust) then
             vpv=vpc
             vph=vpc
             vsv=vsc
             vsh=vsc
             rho=rhoc
             if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
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

! distinguish between single and double precision for reals
       if(CUSTOM_REAL == SIZE_REAL) then
         rhostore(i,j,k,ispec) = sngl(rho)
         kappavstore(i,j,k,ispec) = sngl(rho*(vpv*vpv - 4.d0*vsv*vsv/3.d0))
         kappahstore(i,j,k,ispec) = sngl(rho*(vph*vph - 4.d0*vsh*vsh/3.d0))
         muvstore(i,j,k,ispec) = sngl(rho*vsv*vsv)
         muhstore(i,j,k,ispec) = sngl(rho*vsh*vsh)
         eta_anisostore(i,j,k,ispec) = sngl(eta_aniso)

         if(ABSORBING_CONDITIONS) then
           if(iregion_code == IREGION_OUTER_CORE) then
! we need just vp in the outer core for STacey conditions
             rho_vp(i,j,k,ispec) = sngl(vph)
             rho_vs(i,j,k,ispec) = sngl(0.d0)
           else
             rho_vp(i,j,k,ispec) = sngl(rho*vph)
             rho_vs(i,j,k,ispec) = sngl(rho*vsh)
           endif
         endif

         if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
           c11store(i,j,k,ispec) = sngl(c11)
           c33store(i,j,k,ispec) = sngl(c33)
           c12store(i,j,k,ispec) = sngl(c12)
           c13store(i,j,k,ispec) = sngl(c13)
           c44store(i,j,k,ispec) = sngl(c44)
         endif

         if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
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

         if(ABSORBING_CONDITIONS) then
           if(iregion_code == IREGION_OUTER_CORE) then
! we need just vp in the outer core for STacey conditions
             rho_vp(i,j,k,ispec) = vph
             rho_vs(i,j,k,ispec) = 0.d0
           else
             rho_vp(i,j,k,ispec) = rho*vph
             rho_vs(i,j,k,ispec) = rho*vsh
           endif
         endif

         if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
           c11store(i,j,k,ispec) = c11
           c33store(i,j,k,ispec) = c33
           c12store(i,j,k,ispec) = c12
           c13store(i,j,k,ispec) = c13
           c44store(i,j,k,ispec) = c44
         endif

         if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
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

       if(ATTENUATION .and. ATTENUATION_3D) then
          tau_e_store(:,i,j,k,ispec) = tau_e(:)
          Qmu_store(i,j,k,ispec)     = Qmu
       endif

     enddo
   enddo
 enddo

 call get_flags_boundaries(myrank,iregion_code,nspec,iproc_xi,iproc_eta,ispec,xstore,ystore,zstore, &
        iboun,iMPIcut_xi,iMPIcut_eta,ichunk,idoubling,NPROC_XI,NPROC_ETA, &
        rotation_matrix,ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD, &
             NCHUNKS,INFLATE_CENTRAL_CUBE,R_CENTRAL_CUBE,RCMB,RICB)

 end subroutine get_model

