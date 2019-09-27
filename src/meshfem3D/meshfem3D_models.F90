!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine meshfem3D_models_broadcast()

! preparing model parameter coefficients on all processes

  use shared_parameters, only: MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,LOCAL_PATH,R80,R220,R670,RCMB,RICB

  use meshfem3D_models_par

  implicit none

  ! local parameters
  integer :: ier

  ! sets up spline coefficients for ellipticity
  if (ELLIPTICITY) call make_ellipticity(nspl,rspl,espl,espl2,ONE_CRUST)

  ! read topography and bathymetry file
  if (TOPOGRAPHY) then
    ! arrays for elevations
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if (ier /= 0) stop 'Error allocating ibathy_topo array'

    ! initializes
    ibathy_topo(:,:) = 0

    ! sets up topo/bathy
    call model_topo_bathy_broadcast(ibathy_topo,LOCAL_PATH)
  endif

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! reads 1D reference models
  call meshfem3D_reference_model_broadcast()

  ! reads in 3D mantle models
  call meshfem3D_mantle_broadcast()

  ! reads in crustal model
  call meshfem3D_crust_broadcast()

  ! attenuation
  if (ATTENUATION) then
    call model_attenuation_broadcast(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD)

    ! 3D attenuation
    if (ATTENUATION_3D) then
      ! Colleen's model defined originally between 24.4km and 650km
      call model_atten3D_QRFSI12_broadcast()
    else
      ! sets up attenuation coefficients according to the chosen, "pure" 1D model
      ! (including their 1D-crustal profiles)
      call model_attenuation_setup(REFERENCE_1D_MODEL,RICB,RCMB,R670,R220,R80,CRUSTAL)
    endif

  endif

  end subroutine meshfem3D_models_broadcast


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_reference_model_broadcast()

! preparing model parameter coefficients on all processes for reference models

  use meshfem3D_models_par

  implicit none

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! re-defines/initializes models 1066a and ak135 and ref
  ! ( with possible external crustal model: if CRUSTAL is set to true
  !    it strips the 1-D crustal profile and replaces it with mantle properties)
  select case (REFERENCE_1D_MODEL)
    case (REFERENCE_MODEL_1066A)
      call model_1066a_broadcast(CRUSTAL)

    case (REFERENCE_MODEL_AK135F_NO_MUD)
      call model_ak135_broadcast(CRUSTAL)

    case (REFERENCE_MODEL_1DREF)
      call model_1dref_broadcast(CRUSTAL)

    case (REFERENCE_MODEL_SEA1D)
      call model_sea1d_broadcast(CRUSTAL)
  end select

  end subroutine meshfem3D_reference_model_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_mantle_broadcast()

! preparing model parameter coefficients on all processes for mantle models

  use constants, only: myrank

  use meshfem3D_models_par

  implicit none

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! GLL model uses S29EA as reference 3D model
  if (THREE_D_MODEL == THREE_D_MODEL_GLL) then
    ! sets to initial reference model from which iterations started
    THREE_D_MODEL = GLL_REFERENCE_MODEL
  endif

  ! 3D mantle models
  if (MODEL_3D_MANTLE_PERTUBATIONS) then

    select case (THREE_D_MODEL)

      case (THREE_D_MODEL_S20RTS)
        call model_s20rts_broadcast()

      case (THREE_D_MODEL_S40RTS)
        call model_s40rts_broadcast()

      case(THREE_D_MODEL_MANTLE_SH)
        call model_mantle_sh_broadcast()

      case (THREE_D_MODEL_SEA99_JP3D)
        ! the variables read are declared and stored in structure model_sea99_s_par and model_jp3d_par
        call model_sea99_s_broadcast()
        call model_jp3d_broadcast()

      case (THREE_D_MODEL_SEA99)
        ! the variables read are declared and stored in structure model_sea99_s_par
        call model_sea99_s_broadcast()

      case (THREE_D_MODEL_JP3D)
        ! the variables read are declared and stored in structure model_jp3d_par
        call model_jp3d_broadcast()

      case (THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
            THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
        ! the variables read are declared and stored in structure model_s362ani_par
        call model_s362ani_broadcast(THREE_D_MODEL)

      case (THREE_D_MODEL_PPM)
        ! Point Profile Models
        ! the variables read are declared and stored in structure model_ppm_par
        call model_ppm_broadcast()

      case (THREE_D_MODEL_GAPP2)
        ! GAP model
        call model_gapp2_broadcast()

      case (THREE_D_MODEL_SGLOBE,THREE_D_MODEL_SGLOBE_ISO)
        ! SGLOBE-rani model
        call model_sglobe_broadcast()

      case default
        call exit_MPI(myrank,'3D model not defined')

    end select

  endif

  ! arbitrary mantle models
  ! adds additional perturbations on top of a reference 3D model
  if (HETEROGEN_3D_MANTLE) &
    call model_heterogen_mntl_broadcast()

  ! anisotropic mantle
  if (ANISOTROPIC_3D_MANTLE) &
    call model_aniso_mantle_broadcast()

  ! Enclose this in an ifdef so we don't link to netcdf
  ! if we don't need it.
#ifdef USE_CEM
  if (CEM_REQUEST .or. CEM_ACCEPT) &
    call model_cem_broadcast()
#endif

  ! GLL model
  if (MODEL_GLL) &
    call model_gll_broadcast()

  end subroutine meshfem3D_mantle_broadcast


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_crust_broadcast()

! preparing model parameter coefficients on all processes for crustal models

  use meshfem3D_models_par

  implicit none

  ! checks if anything to do
  if (.not. CRUSTAL) return

!---
!
! ADD YOUR MODEL HERE
!
!---

  select case (REFERENCE_CRUSTAL_MODEL)

    case (ICRUST_CRUST1)
      ! crust 1.0
      call model_crust_1_0_broadcast()

    case (ICRUST_CRUST2)
      ! default
      ! crust 2.0
      call model_crust_2_0_broadcast()

    case (ICRUST_CRUSTMAPS)
      ! general crustmaps
      call model_crustmaps_broadcast()

    case (ICRUST_EPCRUST)
      ! EPcrust (regional crustal model for Europe)
      call model_epcrust_broadcast()
      ! by default crust 1.0 (global coverage)
      call model_crust_1_0_broadcast()

    case (ICRUST_CRUST_SH)
      ! SH crustmaps
      call model_crust_sh_broadcast()

    case (ICRUST_EUCRUST)
      ! EUcrust07 Vp crustal structure (regional crustal model)
      call model_eucrust_broadcast()
      ! by default (vs,rho,eta,moho) from crust 1.0 (global coverage)
      call model_crust_1_0_broadcast()

    case default
      stop 'crustal model type not defined'

  end select

  end subroutine meshfem3D_crust_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_get1D_val(iregion_code,idoubling, &
                              r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                              Qkappa,Qmu,RICB,RCMB, &
                              RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                              RMOHO,RMIDDLE_CRUST,ROCEAN)
! reference model values
!
! for a given location radius (r_prem, which is the point's radius with tolerance factor),
! this calculates density and velocities
!
! note: if CRUSTAL is set, it strips the 1-D crustal profile and mantle gets expanded
!          up to the surface.
!          only exception is JP1D...
!
! routine returns: rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu

  use meshfem3D_models_par

  implicit none

  integer :: iregion_code,idoubling
  double precision :: r_prem,rho
  double precision :: vpv,vph,vsv,vsh,eta_aniso
  double precision :: Qkappa,Qmu
  double precision :: RICB,RCMB,RTOPDDOUBLEPRIME,R80,R120,R220,R400, &
    R600,R670,R771,RMOHO,RMIDDLE_CRUST,ROCEAN

  ! local parameters
  double precision :: drhodr,vp,vs

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! gets 1-D reference model parameters
  select case (REFERENCE_1D_MODEL)

    case (REFERENCE_MODEL_PREM)
      ! PREM (by Dziewonski & Anderson) - used also as background for 3D models
      if (TRANSVERSE_ISOTROPY) then
        ! default PREM:
        !   gets anisotropic PREM parameters, with radial anisotropic extension (from moho to surface for crustal model)
        call model_prem_aniso(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                  Qkappa,Qmu,idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
                  R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

        !daniel todo:
        ! specific 3D models with PREM references which would become too fast at shorter periods ( < 40s Love waves)
        !
        !select case (THREE_D_MODEL)
        !
        ! eventually sgloberani, check...
        !case (THREE_D_MODEL_SGLOBE,THREE_D_MODEL_SGLOBE_ISO)
        !  ! gets anisotropic PREM parameters, with isotropic extension (from moho to surface for crustal model)
        !  call model_prem_aniso_extended_isotropic(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu, &
        !            idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
        !            R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
        !
        ! eventually also Ritsema models, check...
        !case (THREE_D_MODEL_S20RTS,THREE_D_MODEL_S40RTS)
        !  ! gets anisotropic PREM parameters, with isotropic extension (from moho to surface for crustal model)
        !  call model_prem_aniso_extended_isotropic(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu, &
        !            idoubling,CRUSTAL,ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
        !            R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
        !
        !case default
        !  continue
        !end select

      else
        ! isotropic PREM model
        call model_prem_iso(r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
                  ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                  R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
      endif

    case (REFERENCE_MODEL_1DREF)
      ! 1D-REF also known as STW105 (by Kustowski et al.) - used also as background for 3D models
      call model_1dref(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL)
      if (.not. TRANSVERSE_ISOTROPY) then
        if (.not. MODEL_3D_MANTLE_PERTUBATIONS) then
          ! this case here is only executed for 1D_ref_iso
          ! calculates isotropic values
          vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                    + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
          vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                    + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
        endif
      endif

    case (REFERENCE_MODEL_1066A)
      ! 1066A (by Gilbert & Dziewonski) - pure isotropic model, used in 1D model mode only
      call model_1066a(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)

    case (REFERENCE_MODEL_AK135F_NO_MUD)
      ! AK135 (by Kennett et al. ) - pure isotropic model, used in 1D model mode only
      call model_ak135(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)

    case (REFERENCE_MODEL_IASP91)
      ! IASP91 (by Kennett & Engdahl) - pure isotropic model, used in 1D model mode only
      call model_iasp91(r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
                    ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                    R771,R670,R400,R220,R120,RMOHO,RMIDDLE_CRUST)

    case (REFERENCE_MODEL_JP1D)
      !JP1D (by Zhao et al.) - pure isotropic model, used also as background for 3D models
      call model_jp1d(r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
                      .true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                      R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST)

    case (REFERENCE_MODEL_SEA1D)
      ! SEA1D (by Lebedev & Nolet) - pure isotropic model, used also as background for 3D models
      call model_sea1d(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)

    case default
      stop 'unknown 1D reference Earth model in meshfem3D_models_get1D_val()'

  end select

  ! needs to set vpv,vph,vsv,vsh and eta_aniso for isotropic models
  if (.not. TRANSVERSE_ISOTROPY) then
     ! in the case of s362iso we want to save the anisotropic constants for the Voigt average
     if (.not. (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF .and. MODEL_3D_MANTLE_PERTUBATIONS)) then
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0
     endif
  endif ! TRANSVERSE_ISOTROPY

  end subroutine meshfem3D_models_get1D_val


!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho, &
                              vpv,vph,vsv,vsh,eta_aniso, &
                              RCMB,RMOHO, &
                              xmesh,ymesh,zmesh,r, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              ispec,i,j,k)

  use meshfem3D_models_par

  implicit none

  integer, intent (in) :: iregion_code
  double precision, intent(in) :: r_prem
  double precision, intent(out) :: rho
  double precision, intent(out) :: vpv,vph,vsv,vsh,eta_aniso

  double precision :: RCMB,RMOHO
  double precision :: xmesh,ymesh,zmesh,r

  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! heterogen model and CEM needs these (CEM to determine iglob)
  integer, intent (in) :: ispec, i, j, k

  ! local parameters
  double precision :: r_used,r_dummy,theta,phi
  double precision :: dvp,dvs,drho,vp,vs
  double precision :: dvpv,dvph,dvsv,dvsh,deta
  double precision :: lat,lon
  double precision :: A,C,L,N,F

  real(kind=4) :: xcolat,xlon,xrad
  real(kind=4) :: xdvpv,xdvph,xdvsv,xdvsh

  logical :: found_crust,suppress_mantle_extension,is_inside_region

  ! initializes perturbation values
  dvs = ZERO
  dvp = ZERO
  drho = ZERO

  dvpv = ZERO
  dvph = ZERO
  dvsv = ZERO
  dvsh = ZERO
  deta = ZERO

  xdvpv = 0.d0
  xdvph = 0.d0
  xdvsv = 0.d0
  xdvsh = 0.d0

  r_used = ZERO
  suppress_mantle_extension = .false.

  ! gets point's theta/phi
  call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
  call reduce(theta,phi)

!---
!
! ADD YOUR MODEL HERE
!
!---

  if (iregion_code == IREGION_CRUST_MANTLE) then
    ! crust/mantle
    ! sets flag when mantle should not be extended to surface
    if (r_prem >= RMOHO/R_EARTH .and. .not. CRUSTAL) then
      suppress_mantle_extension = .true.
    endif

    ! gets parameters for isotropic 3D mantle model
    !
    ! note: there can be transverse isotropy in the mantle, but only Lame parameters
    !           like kappav,kappah,muv,muh and eta_aniso are used for these simulations
    !
    ! note: in general, models here make use of perturbation values with respect to their
    !          corresponding 1-D reference models
    if (MODEL_3D_MANTLE_PERTUBATIONS .and. r_prem > RCMB/R_EARTH .and. .not. suppress_mantle_extension) then

      ! extend 3-D mantle model above the Moho to the surface before adding the crust
      if (r_prem > RCMB/R_EARTH .and. r_prem < RMOHO/R_EARTH) then
        ! GLL point is in mantle region, takes exact location
        r_used = r
      else ! else if (r_prem >= RMOHO/R_EARTH) then
        if (CRUSTAL) then
          ! GLL point is above moho
          ! takes radius slightly below moho radius, this will then "extend the mantle up to the surface";
          ! crustal values will be superimposed later on
          !
          ! note: this assumes that all the following mantle models are defined below RMOHO.
          !       this is in general true, almost all mantle models fix the moho depth at 24.4 km (set as RMOHO reference)
          !       and invert their velocity models for structures below.
          !
          !       however, in case a mantle models extends to shallower depths, those velocity regions will be ignored.
          !
          ! to do in future: we might want to let the mantle routines decide where to this upper cut-off radius value
          r_used = 0.999999d0*RMOHO/R_EARTH
        endif
      endif

      ! gets model parameters
      select case (THREE_D_MODEL)

        case (THREE_D_MODEL_S20RTS)
          ! s20rts
          call mantle_s20rts(r_used,theta,phi,dvs,dvp,drho)
          vpv = vpv*(1.0d0+dvp)
          vph = vph*(1.0d0+dvp)
          vsv = vsv*(1.0d0+dvs)
          vsh = vsh*(1.0d0+dvs)
          rho = rho*(1.0d0+drho)

        case (THREE_D_MODEL_S40RTS)
          ! s40rts
          call mantle_s40rts(r_used,theta,phi,dvs,dvp,drho)
          vpv = vpv*(1.0d0+dvp)
          vph = vph*(1.0d0+dvp)
          vsv = vsv*(1.0d0+dvs)
          vsh = vsh*(1.0d0+dvs)
          rho = rho*(1.0d0+drho)

        case(THREE_D_MODEL_MANTLE_SH)
          ! full_sh model
          lat = (PI/2.0d0-theta)*180.0d0/PI
          lon = phi*180.0d0/PI
          if (lon > 180.0d0) lon = lon - 360.0d0

          call mantle_sh(lat,lon,r_used,dvpv,dvph,dvsv,dvsh,deta,drho)
          vpv = vpv*(1.0d0+dvpv)
          vph = vph*(1.0d0+dvph)
          vsv = vsv*(1.0d0+dvsv)
          vsh = vsh*(1.0d0+dvsh)
          eta_aniso = eta_aniso*(1.0d0+deta)
          rho = rho*(1.0d0+drho)

        case (THREE_D_MODEL_SEA99_JP3D)
          ! sea99 + jp3d1994
          call model_sea99_s(r_used,theta,phi,dvs)
          vsv = vsv*(1.0d0+dvs)
          vsh = vsh*(1.0d0+dvs)
          ! use Lebedev model sea99 as background and add vp & vs perturbation from Zhao 1994 model jp3d
          call model_jp3d_iso_zhao(r_used,theta,phi,vp,vs,dvp,dvs,rho,found_crust,is_inside_region)
          if (is_inside_region) then
            vpv = vpv*(1.0d0+dvp)
            vph = vph*(1.0d0+dvp)
            vsv = vsv*(1.0d0+dvs)
            vsh = vsh*(1.0d0+dvs)
          endif

        case (THREE_D_MODEL_SEA99)
          ! sea99 Vs-only
          call model_sea99_s(r_used,theta,phi,dvs)
          vsv = vsv*(1.0d0+dvs)
          vsh = vsh*(1.0d0+dvs)

        case (THREE_D_MODEL_JP3D)
          ! jp3d1994
          call model_jp3d_iso_zhao(r_used,theta,phi,vp,vs,dvp,dvs,rho,found_crust,is_inside_region)
          if (is_inside_region) then
            vpv = vpv*(1.0d0+dvp)
            vph = vph*(1.0d0+dvp)
            vsv = vsv*(1.0d0+dvs)
            vsh = vsh*(1.0d0+dvs)
          endif

        case (THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
              THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
          ! 3D Harvard models s362ani, s362wmani, s362ani_prem and s2.9ea
          xcolat = sngl(theta*180.0d0/PI)
          xlon = sngl(phi*180.0d0/PI)
          xrad = sngl(r_used*R_EARTH_KM)
          call model_s362ani_subshsv(xcolat,xlon,xrad,xdvsh,xdvsv,xdvph,xdvpv)

          ! to use speed values from the 1D reference model but with 3D mesh variations
          if (USE_1D_REFERENCE) then
            ! sets all 3D variations in the mantle to zero
            xdvpv = 0.d0
            xdvph = 0.d0
            xdvsv = 0.d0
            xdvsh = 0.d0
          endif

          if (TRANSVERSE_ISOTROPY) then
            ! tiso perturbation
            vpv = vpv*(1.0d0+dble(xdvpv))
            vph = vph*(1.0d0+dble(xdvph))
            vsv = vsv*(1.0d0+dble(xdvsv))
            vsh = vsh*(1.0d0+dble(xdvsh))
          else
            ! isotropic model
            vpv = vpv+xdvpv
            vph = vph+xdvph
            vsv = vsv+xdvsv
            vsh = vsh+xdvsh
            ! isotropic average (considers anisotropic parameterization eta,vsv,vsh,vpv,vph)
            vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                      + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
            vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                      + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
            vpv = vp
            vph = vp
            vsv = vs
            vsh = vs
            eta_aniso = 1.0d0
          endif

        case (THREE_D_MODEL_PPM)
          ! point profile model
          call model_PPM(r_used,theta,phi,dvs,dvp,drho)
          vpv = vpv*(1.0d0+dvp)
          vph = vph*(1.0d0+dvp)
          vsv = vsv*(1.0d0+dvs)
          vsh = vsh*(1.0d0+dvs)
          rho = rho*(1.0d0+drho)

        case (THREE_D_MODEL_GAPP2)
          ! 3D GAP model (Obayashi)
          call mantle_gapmodel(r_used,theta,phi,dvs,dvp,drho)
          vpv = vpv*(1.0d0+dvp)
          vph = vph*(1.0d0+dvp)
          vsv = vsv*(1.0d0+dvs)
          vsh = vsh*(1.0d0+dvs)
          rho = rho*(1.0d0+drho)

        case (THREE_D_MODEL_SGLOBE,THREE_D_MODEL_SGLOBE_ISO)
          ! 3D SGLOBE-rani model (Chang)

          ! normally mantle perturbations are taken from 24.4km (R_MOHO) up.
          ! we need to add the if statement for sgloberani_iso or sgloberani_aniso to take from 50km up:
          if (r_prem > RCMB/R_EARTH .and. r_prem < 6321000.d0/R_EARTH) then
            r_used = r
          else   ! if (r_prem >= 6321000.d0/R_EARTH) then
            ! this will then "extend the mantle up to the surface" from 50km depth
            r_used = 6321000.d0/R_EARTH
          endif

          call mantle_sglobe(r_used,theta,phi,dvsv,dvsh,dvp,drho)

          if (TRANSVERSE_ISOTROPY) then
            ! tiso perturbation
            vpv = vpv*(1.0d0+dvp)
            vph = vph*(1.0d0+dvp)
            vsv = vsv*(1.0d0+dvsv)
            vsh = vsh*(1.0d0+dvsh)
          else
            ! isotropic model
            vpv = vpv*(1.0d0+dvp)
            vph = vph*(1.0d0+dvp)
            vsv = vsv*(1.0d0+dvsv)
            vsh = vsh*(1.0d0+dvsh)
            ! Voigt average
            vp = sqrt( (2.d0*vpv**2 + vph**2)/3.d0 )
            vs = sqrt( (2.d0*vsv**2 + vsh**2)/3.d0 )
            vph = vp
            vpv = vp
            vsh = vs
            vsv = vs
            eta_aniso = 1.d0
          endif
          rho = rho*(1.0d0+drho)

        case (THREE_D_MODEL_ANISO_MANTLE)
          ! Montagner anisotropic model
          call model_aniso_mantle(r_used,theta,phi,vpv,vph,vsv,vsh,eta_aniso,rho, &
                                  c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                  c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

        case default
          print *,'Error: do not recognize value for THREE_D_MODEL ',THREE_D_MODEL
          stop 'unknown 3D Earth model in meshfem3D_models_get3Dmntl_val(), please check... '

      end select ! THREE_D_MODEL

    endif ! MODEL_3D_MANTLE_PERTUBATIONS

    ! heterogen model
    ! adds additional mantle perturbations
    ! (based on density variations) on top of reference 3D model
    if (HETEROGEN_3D_MANTLE .and. .not. suppress_mantle_extension) then
      ! gets spherical coordinates of actual point location
      call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_used,theta,phi)
      call reduce(theta,phi)
      ! adds hetergeneous perturbations (isotropic)
      call model_heterogen_mantle(ispec,i,j,k,r_used,theta,phi,dvs,dvp,drho)
      vpv = vpv*(1.0d0+dvp)
      vph = vpv*(1.0d0+dvp)
      vsv = vsv*(1.0d0+dvs)
      vsh = vsh*(1.0d0+dvs)
      rho = rho*(1.0d0+drho)
    endif ! HETEROGEN_3D_MANTLE

  else if (iregion_code == IREGION_INNER_CORE) then
    ! inner core 3D models
    select case (THREE_D_MODEL_IC)
    case (THREE_D_MODEL_INNER_CORE_ISHII)
      ! Ishii et al. (2002) model
      call model_aniso_inner_core(r_prem,c11,c33,c12,c13,c44,REFERENCE_1D_MODEL, &
                                  vpv,vph,vsv,vsh,rho,eta_aniso)
    case default
      ! nothing to do yet, takes inner core velocities from 1D reference model...
      continue

    end select

  else if (iregion_code == IREGION_OUTER_CORE) then
    ! no 3D outer core models implemented yet... nothing to do
    continue

  else
    ! case should not occur
    print *,'Error invalid iregion_code',iregion_code
    stop 'Invalid iregion_code case in meshfem3D_models_get3Dmntl_val() routine'
  endif

  ! collaborative earth model for whole globe
#ifdef USE_CEM
  if (CEM_ACCEPT) then
    ! over-imposes velocity values for all regions (crust/mantle,outer core,inner core)
    call request_cem (vsh, vsv, vph, vpv, rho, iregion_code, ispec, i, j, k)
  endif
#endif

  ! parameterization
  !
  ! converts isotropic/tiso parameters to fully anisotropic parameters
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

    ! anisotropic model between the Moho and 670 km (change to CMB if desired)
    if (THREE_D_MODEL == THREE_D_MODEL_ANISO_MANTLE) then
      ! special case for model_aniso_mantle()
      ! c11,.. already set in model_aniso_mantle() routine
      ! nothing to do left
      continue
    else
      ! convert isotropic/tiso parameters to full cij
      if (TRANSVERSE_ISOTROPY) then
        ! assume that transverse isotropy is given for a radial symmetry axis
        ! we will need to rotate from this local (radial) axis to the SPECFEM global reference
        !
        ! Love parameterization
        A = rho * vph**2
        C = rho * vpv**2
        L = rho * vsv**2
        N = rho * vsh**2
        F = eta_aniso * (A - 2.d0 * L)

        ! local (radial) coordinate system to global SPECFEM reference
        call rotate_tensor_Love_to_global(theta,phi, &
                                          A,C,N,L,F, &
                                          c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                          c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

      else
        ! calculates isotropic values from given (transversely isotropic) reference values
        ! (are non-dimensionalized)
        !
        ! note: in case vpv == vph and vsv == vsh and eta == 1,
        !       this reduces to vp == vpv and vs == vsv
        vp = sqrt( ((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                  + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0 )
        vs = sqrt( ((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                          + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0 )

        ! fills the rest of the mantle with the isotropic model
        !
        ! note: no rotation needed as for isotropic case, there is no pre-defined symmetry axis
        ! (cij becomes rotation invariant)
        c11 = rho*vp*vp
        c12 = rho*(vp*vp - 2.d0*vs*vs)
        c13 = c12
        c14 = 0.d0
        c15 = 0.d0
        c16 = 0.d0
        c22 = c11
        c23 = c12
        c24 = 0.d0
        c25 = 0.d0
        c26 = 0.d0
        c33 = c11
        c34 = 0.d0
        c35 = 0.d0
        c36 = 0.d0
        c44 = rho*vs*vs
        c45 = 0.d0
        c46 = 0.d0
        c55 = c44
        c56 = 0.d0
        c66 = c44
      endif
    endif
  endif ! ANISOTROPIC_3D_MANTLE

  if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
    ! calculates isotropic values from given (transversely isotropic) reference values
    ! (are non-dimensionalized)
    !
    ! note: in case vpv == vph and vsv == vsh and eta == 1,
    !       this reduces to vp == vpv and vs == vsv
    vp = sqrt( ((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
              + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0 )
    vs = sqrt( ((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                      + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0 )

    ! elastic tensor for hexagonal symmetry in reduced notation:
    !
    !      c11 c12 c13  0   0        0
    !      c12 c11 c13  0   0        0
    !      c13 c13 c33  0   0        0
    !       0   0   0  c44  0        0
    !       0   0   0   0  c44       0
    !       0   0   0   0   0  c66=(c11-c12)/2
    !
    !       in terms of the A, C, L, N and F of Love (1927):
    !
    !       c11 = A
    !       c33 = C
    !       c12 = A-2N
    !       c13 = F
    !       c44 = L
    !       c66 = N
    !
    !       isotropic equivalent:
    !
    !       c11 = lambda+2mu
    !       c33 = lambda+2mu
    !       c12 = lambda
    !       c13 = lambda
    !       c44 = mu
    !       c66 = mu

    ! fills the rest of the mantle with the isotropic model
    !
    ! note: no rotation needed as for isotropic case, there is no pre-defined symmetry axis
    c11 = rho*vp*vp
    c12 = rho*(vp*vp - 2.d0*vs*vs)
    c13 = c12
    c22 = c11
    c23 = c12
    c33 = c11
    c44 = rho*vs*vs
    c55 = c44
    c66 = c44
  endif

  end subroutine meshfem3D_models_get3Dmntl_val

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_get3Dcrust_val(iregion_code,xmesh,ymesh,zmesh,r, &
                                             vpv,vph,vsv,vsh,rho,eta_aniso, &
                                             c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                             c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                             elem_in_crust,moho)

! returns velocities and density for points in 3D crustal region

  use meshfem3D_models_par

  implicit none

  integer,intent(in) :: iregion_code
  ! note: r is the exact radius (and not r_prem with tolerance)
  double precision,intent(in) :: xmesh,ymesh,zmesh,r
  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,eta_aniso

  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  logical,intent(in) :: elem_in_crust
  double precision,intent(out) :: moho

  ! local parameters
  double precision :: r_dummy,theta,phi
  double precision :: lat,lon
  double precision :: vpvc,vphc,vsvc,vshc,etac
  double precision :: vpc,vsc,rhoc !vpc_eu

  double precision :: dvs,dvp
  logical :: found_crust,is_inside_region

  ! checks if anything to do, that is, there is nothing to do
  ! for point radius smaller than deepest possible crust radius (~80 km depth)
  if (r < R_DEEPEST_CRUST ) return

  ! gets point's position theta/phi, lat/lon
  call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
  call reduce(theta,phi)

  ! lat/lon in degrees (range lat/lon = [-90,90] / [-180,180]
  lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  lon = phi * RADIANS_TO_DEGREES
  if (lon > 180.0d0 ) lon = lon - 360.0d0

!---
!
! ADD YOUR MODEL HERE
!
!---
  found_crust = .false.

  ! crustal model can vary for different 3-D models
  select case (THREE_D_MODEL)

    case (THREE_D_MODEL_SEA99_JP3D,THREE_D_MODEL_JP3D)
      ! tries to use Zhao's model of the crust
      call model_jp3d_iso_zhao(r,theta,phi,vpc,vsc,dvp,dvs,rhoc,found_crust,is_inside_region)
      if (.not. is_inside_region) then
        ! uses default crust outside of model region
        call meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,found_crust,elem_in_crust)
      endif

    case default
      ! default crust
      call meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,found_crust,elem_in_crust)

  end select

  ! sets crustal values
  if (found_crust) then
    vpv = vpvc
    vph = vphc
    vsv = vsvc
    vsh = vshc
    rho = rhoc
    eta_aniso = etac

    ! sets anisotropy in crustal region as well
    if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
      ! equivalent with an isotropic elastic tensor (given vpv and vsv as isotropic wave speeds)
      ! note: todo - this could be written as a transversely isotropic tensor (given vphc,vpvc,vshc,vsvc and etac from above)
      c11 = rho * vpv*vpv
      c12 = rho * (vpv*vpv - 2.d0*vsv*vsv)
      c13 = c12
      c14 = 0.d0
      c15 = 0.d0
      c16 = 0.d0
      c22 = c11
      c23 = c12
      c24 = 0.d0
      c25 = 0.d0
      c26 = 0.d0
      c33 = c11
      c34 = 0.d0
      c35 = 0.d0
      c36 = 0.d0
      c44 = rho * vsv*vsv
      c45 = 0.d0
      c46 = 0.d0
      c55 = c44
      c56 = 0.d0
      c66 = c44
    endif
  endif

  end subroutine meshfem3D_models_get3Dcrust_val

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,found_crust,elem_in_crust)

! returns velocity/density for default crust

  use meshfem3D_models_par

  implicit none

  double precision,intent(in) :: lat,lon,r
  double precision,intent(out) :: vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision,intent(out) :: moho
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust

  ! local parameters
  ! for isotropic crust
  double precision :: vpc,vsc
  double precision :: vpc_area,vsc_area,rhoc_area,moho_area
  logical :: found_crust_area,point_in_area

  ! initializes
  vpvc = 0.d0
  vphc = 0.d0
  vsvc = 0.d0
  vshc = 0.d0

  vpc = 0.d0
  vsc = 0.d0
  rhoc = 0.d0

  ! isotropic by default
  etac = 1.d0

  ! moho depth
  moho = 0.d0

  ! flag to indicate if position inside crust
  found_crust = .false.

!---
!
! ADD YOUR MODEL HERE
!
!---
  ! lat/lon range: [-90,90] / [-180,180]

  select case (REFERENCE_CRUSTAL_MODEL)

    case (ICRUST_CRUST1)
      ! crust 1.0
      call model_crust_1_0(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_CRUST2)
      ! default
      ! crust 2.0
      call model_crust_2_0(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_CRUSTMAPS)
      ! general crustmaps
      call model_crustmaps(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_EPCRUST)
      ! if defined within lat/lon-range, takes vp/vs/rho/moho from eucrust07
      call model_epcrust(lat,lon,r,vpc_area,vsc_area,rhoc_area,moho_area,found_crust_area,elem_in_crust,point_in_area)
      if (point_in_area) then
        vpvc = vpc_area
        vphc = vpc_area
        vsvc = vsc_area
        vshc = vsc_area
        rhoc = rhoc_area
        moho = moho_area
        found_crust = found_crust_area
      else
        ! by default takes Crust1.0 values
        call model_crust_1_0(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)
        vpvc = vpc
        vphc = vpc
        vsvc = vsc
        vshc = vsc
      endif

    case (ICRUST_CRUST_SH)
      ! SH crust: provides TI crust
      call crust_sh(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,found_crust,elem_in_crust)

    case (ICRUST_EUCRUST)
      ! by default takes Crust1.0 values for vs/vp/rho/moho
      call model_crust_1_0(lat,lon,r,vpc,vsc,rhoc,moho,found_crust,elem_in_crust)
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc
      ! if defined within lat/lon-range, takes vp/moho from eucrust07
      call model_eucrust(lat,lon,r,vpc_area,moho_area,found_crust_area,point_in_area)
      if (point_in_area) then
        vpvc = vpc_area
        vphc = vpc_area
        moho = moho_area
        found_crust = found_crust_area
      endif

    case default
      stop 'crustal model type not defined'

  end select


  end subroutine meshfem3D_model_crust

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_getatten_val(idoubling,xmesh,ymesh,zmesh,r_prem, &
                                           tau_e,tau_s,T_c_source, &
                                           moho,Qmu,Qkappa,elem_in_crust)

! sets attenuation values tau_e and Qmu for a given point
!
! note:  only Qmu attenuation considered, Qkappa attenuation not used so far in solver...

  use constants, only: N_SLS

  use meshfem3D_models_par

  implicit none

  integer :: idoubling

  double precision :: xmesh,ymesh,zmesh

  double precision :: r_prem
  double precision :: moho

  ! attenuation values
  double precision :: Qkappa,Qmu
  double precision, dimension(N_SLS) :: tau_s, tau_e
  double precision  :: T_c_source

  logical,intent(in) :: elem_in_crust

  ! local parameters
  double precision :: r_dummy,theta,phi,theta_degrees,phi_degrees
  double precision :: r_used
  double precision, parameter :: rmoho_prem = R_EARTH_KM - 24.4d0

  ! initializes
  tau_e(:)   = 0.0d0

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! Get the value of Qmu (Attenuation) dependent on
  ! the radius (r_prem) and idoubling flag
  if (ATTENUATION_3D) then
    ! used for models: s362ani_3DQ, s362iso_3DQ, 3D_attenuation

    ! gets spherical coordinates
    call xyz_2_rthetaphi_dble(xmesh,ymesh,zmesh,r_dummy,theta,phi)
    call reduce(theta,phi)
    theta_degrees = theta / DEGREES_TO_RADIANS
    phi_degrees = phi / DEGREES_TO_RADIANS

    ! in case models incorporate a 3D crust, attenuation values for mantle
    ! get expanded up to surface, and for the crustal points Qmu for PREM crust is imposed
    r_used = r_prem*R_EARTH_KM
    if (CRUSTAL) then
      if (r_prem > (ONE-moho) .or. elem_in_crust) then
        ! points in actual crust: puts point radius into prem crust
        r_used = rmoho_prem*1.0001
      else if (r_prem*R_EARTH_KM >= rmoho_prem) then
        ! points below actual crust (e.g. oceanic crust case), but above prem moho:
        ! puts point slightly below prem moho to expand mantle values at that depth
        r_used = rmoho_prem*0.99999
      endif
    endif ! CRUSTAL

    ! gets value according to radius/theta/phi location and idoubling flag
    call model_atten3D_QRFSI12(r_used,theta_degrees,phi_degrees,Qmu,idoubling)

  else

    select case (REFERENCE_1D_MODEL)

      ! case (REFERENCE_MODEL_PREM)
      ! this case is probably not needed since Qmu is 600. between R80 and surface
      !   call model_attenuation_1D_PREM(r_prem, Qmu)

      case (REFERENCE_MODEL_1DREF)
        ! 1D Ref changes Qmu at moho depth of 24.4km
        ! we take the crustal value and assign it to points only inside actual crust,
        ! otherwise the mantle values is taken
        ! makes sense especially for points below thin oceanic and thick continental crust
        if (CRUSTAL) then
          ! takes crustal Q value only if point is in actual crust
          if (r_prem > (ONE-moho) .or. elem_in_crust) then
            ! reference from 1D-REF aka STW105
            Qmu=300.0d0
            Qkappa=57822.5d0 !  not used so far...
          endif
        endif ! CRUSTAL

      case (REFERENCE_MODEL_SEA1D)
        ! SEA1D changes Qmu at 25km (moho) depth. we take the crustal value
        ! for points only inside actual crust
        if (CRUSTAL) then
          ! takes crustal Q value only if point is in actual crust
          if (r_prem > (ONE-moho) .or. elem_in_crust) then
            ! reference from Sea1D
            Qmu = 300.0d0
            Qkappa = 57822.5d0  ! not used so far...
          endif
        endif

    end select

  endif

  ! Get tau_e from tau_s and Qmu
  call model_attenuation_getstored_tau(Qmu, T_c_source, tau_s, tau_e)

  end subroutine meshfem3D_models_getatten_val


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_impose_val(iregion_code,xmesh,ymesh,zmesh,ispec,i,j,k, &
                                         vpv,vph,vsv,vsh,rho,eta_aniso, &
                                         c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                         c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! overwrites values with updated model values (from iteration step) here, given at all GLL points

  use meshfem3D_models_par

  implicit none

  integer,intent(in) :: iregion_code,ispec,i,j,k
  double precision,intent(in) :: xmesh,ymesh,zmesh

  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,eta_aniso
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66


  ! checks if anything to do
  if (.not. MODEL_GLL) return

  ! only valid for crust/mantle region at the moment...
  if (iregion_code /= IREGION_CRUST_MANTLE) return

  ! over-impose values from model GLL values
  call model_gll_impose_val(xmesh,ymesh,zmesh,ispec,i,j,k, &
                            vpv,vph,vsv,vsh,rho,eta_aniso, &
                            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine meshfem3D_models_impose_val


