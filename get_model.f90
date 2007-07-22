!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine get_model(myrank,iregion_code,nspec, &
    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore,rhostore, &
    nspec_ani, &
    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
    xelm,yelm,zelm,shape3D,ispec, &
    rmin,rmax,idoubling, &
    rho_vp,rho_vs,nspec_stacey, &
    TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE, &
    CRUSTAL,ONE_CRUST, &
    crustal_model,mantle_model,aniso_mantle_model, &
    aniso_inner_core_model,&
    attenuation_model,ATTENUATION,ATTENUATION_3D,tau_s,tau_e_store,Qmu_store,T_c_source,vx,vy,vz,vnspec, &
    ABSORBING_CONDITIONS,REFERENCE_1D_MODEL,THREE_D_MODEL, &
    RCMB,RICB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN,&
    AMM_V, AM_V, M1066a_V, Mak135_V,Mref_V,D3MM_V,CM_V, AM_S, AS_V)

  implicit none

  include "constants.h"

  external mantle_model,crustal_model,aniso_mantle_model, &
       aniso_inner_core_model,attenuation_model

! aniso_mantle_model_variables
  type aniso_mantle_model_variables
    sequence
    double precision beta(14,34,37,73)
    double precision pro(47)
    integer npar1
  end type aniso_mantle_model_variables

  type (aniso_mantle_model_variables) AMM_V
! aniso_mantle_model_variables

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

  type (attenuation_model_variables) AM_V
! attenuation_model_variables

! model_1066a_variables
  type model_1066a_variables
    sequence
      double precision, dimension(NR_1066A) :: radius_1066a
      double precision, dimension(NR_1066A) :: density_1066a
      double precision, dimension(NR_1066A) :: vp_1066a
      double precision, dimension(NR_1066A) :: vs_1066a
      double precision, dimension(NR_1066A) :: Qkappa_1066a
      double precision, dimension(NR_1066A) :: Qmu_1066a
  end type model_1066a_variables

  type (model_1066a_variables) M1066a_V
! model_1066a_variables

! model_ak135_variables
  type model_ak135_variables
    sequence
    double precision, dimension(NR_AK135) :: radius_ak135
    double precision, dimension(NR_AK135) :: density_ak135
    double precision, dimension(NR_AK135) :: vp_ak135
    double precision, dimension(NR_AK135) :: vs_ak135
    double precision, dimension(NR_AK135) :: Qkappa_ak135
    double precision, dimension(NR_AK135) :: Qmu_ak135
  end type model_ak135_variables

 type (model_ak135_variables) Mak135_V
! model_ak135_variables

! model_ref_variables
  type model_ref_variables
    sequence
      double precision, dimension(NR_REF) :: radius_ref
      double precision, dimension(NR_REF) :: density_ref
      double precision, dimension(NR_REF) :: vpv_ref
      double precision, dimension(NR_REF) :: vph_ref
      double precision, dimension(NR_REF) :: vsv_ref
      double precision, dimension(NR_REF) :: vsh_ref
      double precision, dimension(NR_REF) :: eta_ref
      double precision, dimension(NR_REF) :: Qkappa_ref
      double precision, dimension(NR_REF) :: Qmu_ref
  end type model_ref_variables

  type (model_ref_variables) Mref_V
! model_ref_variables

! three_d_mantle_model_variables
  type three_d_mantle_model_variables
    sequence
    double precision dvs_a(0:NK,0:NS,0:NS)
    double precision dvs_b(0:NK,0:NS,0:NS)
    double precision dvp_a(0:NK,0:NS,0:NS)
    double precision dvp_b(0:NK,0:NS,0:NS)
    double precision spknt(NK+1)
    double precision qq0(NK+1,NK+1)
    double precision qq(3,NK+1,NK+1)
  end type three_d_mantle_model_variables

  type (three_d_mantle_model_variables) D3MM_V
! three_d_mantle_model_variables

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables

! attenuation_model_storage
  type attenuation_model_storage
    sequence
    integer Q_resolution
    integer Q_max
    double precision, dimension(:,:), pointer :: tau_e_storage
    double precision, dimension(:), pointer :: Qmu_storage
  end type attenuation_model_storage

  type (attenuation_model_storage) AM_S
! attenuation_model_storage

! attenuation_simplex_variables
  type attenuation_simplex_variables
    sequence
    integer nf          ! nf    = Number of Frequencies
    integer nsls        ! nsls  = Number of Standard Linear Solids
    double precision Q  ! Q     = Desired Value of Attenuation or Q
    double precision iQ ! iQ    = 1/Q
    double precision, dimension(:), pointer ::  f
    ! f = Frequencies at which to evaluate the solution
    double precision, dimension(:), pointer :: tau_s
    ! tau_s = Tau_sigma defined by the frequency range and
    !             number of standard linear solids
  end type attenuation_simplex_variables

  type(attenuation_simplex_variables) AS_V
! attenuation_simplex_variables


  integer ispec,nspec,idoubling,iregion_code,myrank,nspec_stacey
  integer REFERENCE_1D_MODEL,THREE_D_MODEL

  logical ATTENUATION,ATTENUATION_3D,ABSORBING_CONDITIONS
  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ISOTROPIC_3D_MANTLE,CRUSTAL,ONE_CRUST

!   logical iboun(6,nspec)
!   logical iMPIcut_xi(2,nspec),iMPIcut_eta(2,nspec)

  double precision shape3D(NGNOD,NGLLX,NGLLY,NGLLZ)

  double precision xelm(NGNOD)
  double precision yelm(NGNOD)
  double precision zelm(NGNOD)

  double precision rmin,rmax,RCMB,RICB,R670,RMOHO, &
    RTOPDDOUBLEPRIME,R600,R220,R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN

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

  double precision xmesh,ymesh,zmesh

  integer i,j,k,ia
  double precision rho,drhodr,vp,vs,Qkappa,Qmu
  double precision vpv,vph,vsv,vsh,eta_aniso
  double precision dvp,dvs,drho
  double precision xstore(NGLLX,NGLLY,NGLLZ)
  double precision ystore(NGLLX,NGLLY,NGLLZ)
  double precision zstore(NGLLX,NGLLY,NGLLZ)
  double precision r,r_prem,r_moho,r_dummy,theta,phi
  double precision lat,lon
  double precision vpc,vsc,rhoc,moho

! attenuation values
  integer vx, vy, vz, vnspec
  double precision, dimension(N_SLS)                     :: tau_s, tau_e
  double precision, dimension(vx, vy, vz, vnspec)        :: Qmu_store
  double precision, dimension(N_SLS, vx, vy, vz, vnspec) :: tau_e_store
  double precision  T_c_source

  logical found_crust
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

         if(REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
           call model_iasp91(myrank,r_prem,rho,vp,vs,Qkappa,Qmu,idoubling,ONE_CRUST, &
             .true.,RICB,RCMB,RTOPDDOUBLEPRIME,R771,R670,R400,R220,R120,RMOHO,RMIDDLE_CRUST)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM) then
           call prem_iso(myrank,r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
             ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
             R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
           call model_1066a(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,M1066a_V)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135) then
           call model_ak135(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code,Mak135_V)

         else if(REFERENCE_1D_MODEL == REFERENCE_MODEL_REF) then
           call model_ref(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,Mref_V)

         else
           stop 'unknown 1D reference Earth model in get_model'
         endif

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
           call mantle_model(r,theta,phi,dvs,dvp,drho,D3MM_V)
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
           call mantle_model(r_moho,theta,phi,dvs,dvp,drho,D3MM_V)
           vpv=vpv*(1.0d0+dvp)
           vph=vph*(1.0d0+dvp)
           vsv=vsv*(1.0d0+dvs)
           vsh=vsh*(1.0d0+dvs)
           rho=rho*(1.0d0+drho)

         endif
       endif

       if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) &
           call aniso_inner_core_model(r_prem,c11,c33,c12,c13,c44,REFERENCE_1D_MODEL)

       if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

! anisotropic model between the Moho and 670 km (change to CMB if desired)
         if(r_prem < RMOHO/R_EARTH .and. r_prem > R670/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           call aniso_mantle_model(r_prem,theta,phi,rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,AMM_V)
! extend 3-D mantle model above the Moho to the surface before adding the crust
         elseif(r_prem >= RMOHO/R_EARTH) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)
           r_moho = RMOHO/R_EARTH
           call aniso_mantle_model(r_moho,theta,phi,rho,c11,c12,c13,c14,c15,c16, &
              c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,AMM_V)
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
         call attenuation_model(r_prem, lat, lon, Qmu)
         call attenuation_conversion(Qmu, T_c_source, tau_s, tau_e, AM_V, AM_S, AS_V)
       endif

!      get the 3-D crustal model
       if(CRUSTAL) then
         if(r > R_DEEPEST_CRUST) then
           call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
           call reduce(theta,phi)

           lat=(PI/2.0d0-theta)*180.0d0/PI
           lon=phi*180.0d0/PI
           if(lon>180.0d0) lon=lon-360.0d0
           call crustal_model(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,CM_V)
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

! we need just vp in the outer core for Stacey conditions
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
! we need just vp in the outer core for Stacey conditions
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

 end subroutine get_model

