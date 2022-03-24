!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

  use shared_parameters, only: LOCAL_PATH,SAVE_MESH_FILES
  use meshfem3D_models_par

  implicit none

  ! local parameters
  integer :: ier

  ! sets up spline coefficients for ellipticity
  if (ELLIPTICITY) call make_ellipticity(nspl,rspl,ellipicity_spline,ellipicity_spline2)

  ! read topography and bathymetry file
  if (TOPOGRAPHY) then
    ! arrays for elevations
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if (ier /= 0) stop 'Error allocating ibathy_topo array'

    ! initializes
    ibathy_topo(:,:) = 0

    ! sets up topo/bathy
    call model_topo_bathy_broadcast(ibathy_topo,LOCAL_PATH)

    ! saves VTK output
    if (SAVE_MESH_FILES) call meshfem3D_plot_VTK_topo_bathy()
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
    call model_attenuation_broadcast()

    ! 3D attenuation
    if (ATTENUATION_GLL) then
      call model_attenuation_gll_broadcast()
    else if (ATTENUATION_3D) then
      ! Colleen's model defined originally between 24.4km and 650km
      call model_atten3D_QRFSI12_broadcast()
    else
      ! sets up attenuation coefficients according to the chosen, "pure" 1D model
      ! (including their 1D-crustal profiles)
      call model_attenuation_setup(REFERENCE_1D_MODEL,CRUSTAL)
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

    case (REFERENCE_MODEL_CASE65TAY)
      call model_case65TAY_broadcast(CRUSTAL)

    case (REFERENCE_MODEL_MARS_1D)
      call model_mars_1D_broadcast(CRUSTAL)

    case (REFERENCE_MODEL_CCREM)
      call model_ccrem_broadcast(CRUSTAL)

  end select

  end subroutine meshfem3D_reference_model_broadcast

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_mantle_broadcast()

! preparing model parameter coefficients on all processes for mantle models

  use constants, only: myrank,IMAIN
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

      case (THREE_D_MODEL_ANISO_MANTLE)
        ! Montagner anisotropic model
        call model_aniso_mantle_broadcast()

      case (THREE_D_MODEL_BKMNS_GLAD)
        ! GLAD model expansion on block-mantle-spherical-harmonics
        if (myrank == 0) write(IMAIN,*) 'setting up both s362ani and bkmns models...'
        ! needs s362ani topo for 410/660
        call model_s362ani_broadcast(THREE_D_MODEL_S362ANI)
        ! block-mantle-spherical-harmonics model
        call model_bkmns_mantle_broadcast()

      case (THREE_D_MODEL_SPIRAL)
        ! SPiRal model
        call model_mantle_spiral_broadcast()

      case (THREE_D_MODEL_HETEROGEN_PREM)
        !chris modif checker 02/20/21
        call model_heterogen_mntl_broadcast()

      case default
        call exit_MPI(myrank,'3D model not defined')

    end select

  endif

  ! adds additional perturbations on top of a reference 3D model
  ! using "arbitrary" mantle models (for example density models derived from geodynamics modeling)
  if (HETEROGEN_3D_MANTLE) &
    call model_heterogen_mntl_broadcast()

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

  use constants, only: myrank,IMAIN
  use shared_parameters, only: SAVE_MESH_FILES
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

    case (ICRUST_SGLOBECRUST)
      ! modified crust2.0 for SGLOBE-rani
      call model_sglobecrust_broadcast()

    case (ICRUST_BKMNS_GLAD)
      ! GLAD model expansion on block-mantle-spherical-harmonics
      if (myrank == 0) write(IMAIN,*) 'setting up both crust 2.0 and bkmns crust models ...'
      ! needs crust 2.0 for moho depths
      call model_crust_2_0_broadcast()
      ! crustal structure in blocks between topography and 80km depth
      call model_bkmns_crust_broadcast()

    case (ICRUST_SPIRAL)
      ! anisotropic crust from SPiRaL
      call model_crust_spiral_broadcast()

    case default
      stop 'crustal model type not defined'

  end select

  ! plot moho as VTK file output
  if (SAVE_MESH_FILES) call meshfem3D_plot_VTK_crust_moho()

  end subroutine meshfem3D_crust_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_get1D_val(iregion_code,idoubling, &
                                        r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                                        Qkappa,Qmu,RICB,RCMB, &
                                        RTOPDDOUBLEPRIME,R80,R120,R220,R400,R670,R771, &
                                        RMOHO,RMIDDLE_CRUST)
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

  integer,intent(in) :: iregion_code,idoubling
  double precision,intent(in) :: r_prem
  double precision,intent(inout) :: vpv,vph,vsv,vsh,eta_aniso,rho
  double precision,intent(inout) :: Qkappa,Qmu
  double precision,intent(in) :: RICB,RCMB,RTOPDDOUBLEPRIME,R80,R120,R220,R400, &
    R670,R771,RMOHO,RMIDDLE_CRUST

  ! local parameters
  double precision :: drhodr,vp,vs

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! gets 1-D reference model parameters
  select case (REFERENCE_1D_MODEL)

    case (REFERENCE_MODEL_PREM,REFERENCE_MODEL_PREM2)
      ! PREM (by Dziewonski & Anderson) - used also as background for 3D models
      if (TRANSVERSE_ISOTROPY) then
        ! default PREM:
        !   gets anisotropic PREM parameters, with radial anisotropic extension (from moho to surface for crustal model)
        call model_prem_aniso(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,idoubling,CRUSTAL)

        ! calculates isotropic values
        vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                  + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
        if (iregion_code == IREGION_OUTER_CORE) then
            ! fluid with zero shear speed
            vs = 0.d0
        else
            vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                  + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
        endif

        !daniel todo:
        ! specific 3D models with PREM references which would become too fast at shorter periods ( < 40s Love waves)
        !
        !select case (THREE_D_MODEL)
        !
        ! eventually sgloberani, check...
        !case (THREE_D_MODEL_SGLOBE,THREE_D_MODEL_SGLOBE_ISO)
        !  ! gets anisotropic PREM parameters, with isotropic extension (from moho to surface for crustal model)
        !  call model_prem_aniso_extended_isotropic(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu, &
        !                                           idoubling,CRUSTAL)
        !
        ! eventually also Ritsema models, check...
        !case (THREE_D_MODEL_S20RTS,THREE_D_MODEL_S40RTS)
        !  ! gets anisotropic PREM parameters, with isotropic extension (from moho to surface for crustal model)
        !  call model_prem_aniso_extended_isotropic(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu, &
        !                                           idoubling,CRUSTAL)
        !
        !case default
        !  continue
        !end select

      else
        ! isotropic PREM model
        call model_prem_iso(r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,.true.)
        vpv = vp
        vph = vp
        vsv = vs
        vsh = vs
        eta_aniso = 1.d0
      endif

    case (REFERENCE_MODEL_1DREF)
      ! 1D-REF also known as STW105 (by Kustowski et al.) - used also as background for 3D models
      call model_1dref(r_prem,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,iregion_code,CRUSTAL)

      ! calculates isotropic values
      vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
      if (iregion_code == IREGION_OUTER_CORE) then
          ! fluid with zero shear speed
          vs = 0.d0
      else
          vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)
      endif
      if (.not. TRANSVERSE_ISOTROPY) then
        if (.not. (MODEL_3D_MANTLE_PERTUBATIONS .and. iregion_code == IREGION_CRUST_MANTLE)) then
          ! this case here is only executed for 1D_ref_iso
          vpv = vp
          vph = vp
          vsv = vs
          vsh = vs
          eta_aniso = 1.d0
        endif
      endif

    case (REFERENCE_MODEL_1066A)
      ! 1066A (by Gilbert & Dziewonski) - pure isotropic model, used in 1D model mode only
      call model_1066a(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_AK135F_NO_MUD)
      ! AK135 (by Kennett et al. ) - pure isotropic model, used in 1D model mode only
      call model_ak135(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_IASP91)
      ! IASP91 (by Kennett & Engdahl) - pure isotropic model, used in 1D model mode only
      call model_iasp91(r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
                        ONE_CRUST,.true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                        R771,R670,R400,R220,R120,RMOHO,RMIDDLE_CRUST)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_JP1D)
      !JP1D (by Zhao et al.) - pure isotropic model, used also as background for 3D models
      call model_jp1d(r_prem,rho,vp,vs,Qkappa,Qmu,idoubling, &
                      .true.,RICB,RCMB,RTOPDDOUBLEPRIME, &
                      R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_SEA1D)
      ! SEA1D (by Lebedev & Nolet) - pure isotropic model, used also as background for 3D models
      call model_sea1d(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_CCREM)
      ! CCREM (by Ma & Tkalcic) - pure isotropic model
      call model_ccrem(r_prem,rho,vp,vs,Qkappa,Qmu,idoubling,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    ! Mars 1D models
    case (REFERENCE_MODEL_SOHL)
      ! Mars
      call model_Sohl(r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,.true.)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_CASE65TAY)
      ! Mars
      call model_case65TAY(r_prem,rho,vp,vs,Qkappa,Qmu,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case (REFERENCE_MODEL_MARS_1D)
      ! Mars
      call model_mars_1D(r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,.true.,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    ! Moon 1D models
    case (REFERENCE_MODEL_VPREMOON)
      ! Moon
      call model_vpremoon(r_prem,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,.true.,iregion_code)
      vpv = vp
      vph = vp
      vsv = vs
      vsh = vs
      eta_aniso = 1.d0

    case default
      stop 'unknown 1D reference Earth model in meshfem3D_models_get1D_val()'

  end select

  ! needs to set vpv,vph,vsv,vsh and eta_aniso for isotropic models
  if (.not. TRANSVERSE_ISOTROPY) then
     ! in the case of s362iso we want to save the anisotropic constants for the Voigt average
     if (.not. (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF &
                .and. MODEL_3D_MANTLE_PERTUBATIONS &
                .and. iregion_code == IREGION_CRUST_MANTLE)) then
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
                                            r,theta,phi, &
                                            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                            ispec,i,j,k)

  use meshfem3D_models_par

  implicit none

  integer, intent(in) :: iregion_code
  double precision, intent(in) :: r_prem
  double precision, intent(inout) :: rho
  double precision, intent(inout) :: vpv,vph,vsv,vsh,eta_aniso

  double precision,intent(in) :: RCMB,RMOHO
  double precision,intent(in) :: r,theta,phi

  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! heterogen model and CEM needs these (CEM to determine iglob)
  integer, intent(in) :: ispec, i, j, k

  ! local parameters
  double precision :: r_used
  double precision :: dvp,dvs,drho,vp,vs,moho,sediment
  double precision :: dvpv,dvph,dvsv,dvsh,deta
  double precision :: lat,lon
  double precision :: A,C,L,N,F

  real(kind=4) :: xcolat,xlon,xrad
  real(kind=4) :: xdvpv,xdvph,xdvsv,xdvsh

  logical :: found_crust,suppress_mantle_extension,is_inside_region,convert_tiso_to_cij

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

!---
!
! ADD YOUR MODEL HERE
!
!---

  if (iregion_code == IREGION_CRUST_MANTLE) then
    ! crust/mantle
    ! sets flag when mantle should not be extended to surface
    if (r_prem >= RMOHO/R_PLANET .and. .not. CRUSTAL) then
      suppress_mantle_extension = .true.
    endif

    ! gets parameters for isotropic 3D mantle model
    !
    ! note: there can be transverse isotropy in the mantle, but only Lame parameters
    !           like kappav,kappah,muv,muh and eta_aniso are used for these simulations
    !
    ! note: in general, models here make use of perturbation values with respect to their
    !          corresponding 1-D reference models
    if (MODEL_3D_MANTLE_PERTUBATIONS .and. r_prem > RCMB/R_PLANET .and. .not. suppress_mantle_extension) then

      ! extend 3-D mantle model above the Moho to the surface before adding the crust
      if (r_prem > RCMB/R_PLANET .and. r_prem < RMOHO/R_PLANET) then
        ! GLL point is in mantle region, takes exact location
        r_used = r
      else ! else if (r_prem >= RMOHO/R_PLANET) then
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
          r_used = 0.999999d0*RMOHO/R_PLANET
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
          call model_jp3d_iso_zhao(r_used,theta,phi,vp,vs,dvp,dvs,rho,moho,sediment,found_crust,is_inside_region)
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
          call model_jp3d_iso_zhao(r_used,theta,phi,vp,vs,dvp,dvs,rho,moho,sediment,found_crust,is_inside_region)
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
          xrad = sngl(r_used*R_PLANET_KM)
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
          if (r_prem > RCMB/R_PLANET .and. r_prem < 6321000.d0/R_PLANET) then
            r_used = r
          else   ! if (r_prem >= 6321000.d0/R_PLANET) then
            ! this will then "extend the mantle up to the surface" from 50km depth
            r_used = 6321000.d0/R_PLANET
          endif

          call mantle_sglobe(r_used,theta,phi,vsv,vsh,dvsv,dvsh,dvp,drho)

          ! updates velocities from reference model
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

        case (THREE_D_MODEL_BKMNS_GLAD)
          ! GLAD model expansion on block-mantle-spherical-harmonics
          ! model takes actual position (includes topography between 80km depth and surface topo)
          r_used = r
          call model_bkmns_mantle(r_used,theta,phi,vpv,vph,vsv,vsh,eta_aniso,rho)

        case (THREE_D_MODEL_SPIRAL)
          ! SPiRaL v1.4 anisotropic model
          lat = (PI/2.0d0-theta)*180.0d0/PI
          lon = phi*180.0d0/PI
          if (lon > 180.0d0) lon = lon - 360.0d0
          call model_mantle_spiral(r_used,lat,lon,vpv,vph,vsv,vsh,eta_aniso,rho, &
                                   c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                   c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

        case (THREE_D_MODEL_HETEROGEN_PREM)
          ! chris modif checkers 02/20/21
          call model_heterogen_mantle(ispec,i,j,k,r_prem,theta,phi,dvs,dvp,drho)
          ! vpv = vpv*(1.0d0+dvp) ! correct format with initial heterogenous model
          vpv = dvp
          vph = dvp
          vsv = dvs
          vsh = dvs
          rho = drho

        case default
          print *,'Error: do not recognize value for THREE_D_MODEL ',THREE_D_MODEL
          stop 'unknown 3D Earth model in meshfem3D_models_get3Dmntl_val(), please check... '

      end select ! THREE_D_MODEL

    endif ! MODEL_3D_MANTLE_PERTUBATIONS

    ! heterogen model
    ! adds additional mantle perturbations
    ! (based on density variations) on top of reference 3D model
    if (HETEROGEN_3D_MANTLE &
        .and. .not. suppress_mantle_extension &
        .and. THREE_D_MODEL /= THREE_D_MODEL_HETEROGEN_PREM) then
      ! gets spherical coordinates of actual point location
      r_used = r
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
      call model_aniso_inner_core(r_prem,c11,c12,c13,c33,c44,REFERENCE_1D_MODEL, &
                                  vpv,vph,vsv,vsh,rho,eta_aniso)
      ! converts to isotropic values
      if (.not. TRANSVERSE_ISOTROPY) then
        ! converts further to iso
        vp = sqrt( ((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                  + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0 )
        vs = sqrt( ((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                  + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0 )
        vpv = vp
        vph = vp
        vsv = vs
        vsh = vs
        eta_aniso = 1.d0
      endif
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
    ! checks models if we need to convert tiso values to full cij
    convert_tiso_to_cij = .true.

    ! special case for fully anisotropic models
    ! anisotropic models defined between the Moho and 670 km (change to CMB if desired)
    if (THREE_D_MODEL == THREE_D_MODEL_ANISO_MANTLE) then
      ! special case for model_aniso_mantle()
      if (suppress_mantle_extension) then
        ! point is in crust, and CRUSTAL == .false. (no 3D crustal model); model_aniso_mantle() hasn't been called.
        ! we still need to convert the 1D reference crustal values (vpv,vph,..) to full cij coefficients
        convert_tiso_to_cij = .true.
      else
        ! c11,.. already set in model_aniso_mantle() routine
        ! nothing left to do
        convert_tiso_to_cij = .false.
      endif
    endif

    ! SPiRal model
    if (THREE_D_MODEL == THREE_D_MODEL_SPIRAL) then
      ! special case for appended _1Dcrust cases
      if (suppress_mantle_extension) then
        ! point is in crust, and CRUSTAL == .false. (no 3D crustal model); model_mantle_spiral() hasn't been called.
        ! we still need to convert the 1D reference crustal values (vpv,vph,..) to full cij coefficients
        convert_tiso_to_cij = .true.
      else
        ! sets c11,.. in model_mantle_spiral() routine
        ! nothing left to do
        convert_tiso_to_cij = .false.
      endif
    endif

    if (convert_tiso_to_cij) then
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
        ! we need to convert tiso/iso values to full cij coefficients
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
    if (THREE_D_MODEL_IC == THREE_D_MODEL_INNER_CORE_ISHII) then
      ! c11,c33,c12,c13,c44 provided by Ishii model
      c11 = c11
      c12 = c12
      c13 = c13
      c33 = c33
      c44 = c44
      ! setting additional ones
      c23 = c13
      c55 = c44
      c66 = 0.5d0 * (c11 - c12)
    else
      if (TRANSVERSE_ISOTROPY) then
        ! transversly isotropic
        ! C11 = A = rho * vph**2
        ! C33 = C = rho * vpv**2
        ! C44 = L = rho * vsv**2
        ! C13 = F = eta * (A - 2*L)
        ! C12 = C11 - 2 C66 = A - 2*N = rho * (vph**2 - 2 * vsh**2)
        ! C22 = C11 = A
        ! C23 = C13 = F
        ! C55 = C44 = L
        ! C66 = N = rho * vsh**2 = (C11-C12)/2
        c11 = rho*vph*vph
        c12 = rho*(vph*vph - 2.d0*vsh*vsh)
        c33 = rho*vpv*vpv
        c44 = rho*vsv*vsv
        c13 = eta_aniso * (c11 - 2.d0 * c44)
        c22 = c11
        c23 = c13
        c55 = c44
        c66 = rho*vsh*vsh
      else
        ! isotropic
        ! calculates isotropic values from given (transversely isotropic) reference values
        ! (are non-dimensionalized)
        !
        ! note: in case vpv == vph and vsv == vsh and eta == 1,
        !       this reduces to vp == vpv and vs == vsv
        vp = sqrt( ((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
                  + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0 )
        vs = sqrt( ((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                  + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0 )
        c11 = rho*vp*vp
        c12 = rho*(vp*vp - 2.d0*vs*vs)
        c13 = c12
        c22 = c11
        c23 = c13
        c33 = c11
        c44 = rho*vs*vs
        c55 = c44
        c66 = c44
      endif
    endif
  endif ! ANISOTROPIC_INNER_CORE

  end subroutine meshfem3D_models_get3Dmntl_val

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_get3Dcrust_val(iregion_code,r,theta,phi, &
                                             vpv,vph,vsv,vsh,rho,eta_aniso, &
                                             c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                             c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                             elem_in_crust,moho,sediment)

! returns velocities and density for points in 3D crustal region

  use meshfem3D_models_par

  implicit none

  integer,intent(in) :: iregion_code
  ! note: r is the exact radius (and not r_prem with tolerance)
  !       theta in [0,PI], phi in [0,2PI]
  double precision,intent(in) :: r,theta,phi
  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,eta_aniso

  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  logical,intent(in) :: elem_in_crust
  double precision,intent(inout) :: moho,sediment

  ! local parameters
  double precision :: lat,lon
  double precision :: vpvc,vphc,vsvc,vshc,etac
  double precision :: vpc,vsc,rhoc !vpc_eu
  double precision :: c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                      c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c
  double precision :: A,C,L,N,F
  double precision :: dvs,dvp
  logical :: found_crust,is_inside_region,moho_only

  ! checks if anything to do, that is, there is nothing to do
  ! for point radius smaller than deepest possible crust radius (~80 km depth)
  if (r < R_DEEPEST_CRUST) return

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
  moho_only = .false.

  ! crustal model can vary for different 3-D models
  select case (THREE_D_MODEL)

    case (THREE_D_MODEL_SEA99_JP3D,THREE_D_MODEL_JP3D)
      ! partial/regional crustal model only
      ! tries to use Zhao's model of the crust
      call model_jp3d_iso_zhao(r,theta,phi,vpc,vsc,dvp,dvs,rhoc,moho,sediment,found_crust,is_inside_region)

      if (.not. is_inside_region) then
        ! uses default crust outside of model region
        call meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only, &
                                   c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                                   c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c)
      endif

    case default
      ! default crust
      call meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only, &
                                 c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                                 c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c)

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

      if (THREE_D_MODEL == THREE_D_MODEL_SPIRAL .and. REFERENCE_CRUSTAL_MODEL == ICRUST_SPIRAL) then
        ! SPiRal has fully anisotropic crustal parameters
        ! c11,.. already set in model_crust_spiral() routine
        c11 = c11c
        c12 = c12c
        c13 = c13c
        c14 = c14c
        c15 = c15c
        c16 = c16c
        c22 = c22c
        c23 = c23c
        c24 = c24c
        c25 = c25c
        c26 = c26c
        c33 = c33c
        c34 = c34c
        c35 = c35c
        c36 = c36c
        c44 = c44c
        c45 = c45c
        c46 = c46c
        c55 = c55c
        c56 = c56c
        c66 = c66c
      else
        if (TRANSVERSE_ISOTROPY) then
          ! converts from a tiso crust parameterization to a fully anisotropic one
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
          ! equivalent with an isotropic elastic tensor (given vpv and vsv as isotropic wave speeds)
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
    endif
  endif

  end subroutine meshfem3D_models_get3Dcrust_val

!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_model_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc, &
                                   moho,sediment,found_crust,elem_in_crust,moho_only, &
                                   c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                                   c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c)

! returns velocity/density for default crust

  use meshfem3D_models_par

  implicit none

  ! lat/lon  - in degrees (range lat/lon = [-90,90] / [-180,180]
  ! radius r - normalized by globe radius [0,1.x]
  double precision,intent(in) :: lat,lon,r

  double precision,intent(out) :: vpvc,vphc,vsvc,vshc,etac,rhoc
  double precision,intent(out) :: moho,sediment
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust,moho_only
  double precision,intent(out) :: c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                                  c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c
  ! local parameters
  ! for isotropic crust
  double precision :: vpc,vsc
  double precision :: vpc_area,vsc_area,rhoc_area,moho_area,sediment_area
  logical :: found_crust_area,point_in_area

  ! initializes
  if (moho_only) then
    ! moho depth
    moho = 0.d0

    ! sediment depth
    sediment = 0.d0
  else
    vpvc = 0.d0
    vphc = 0.d0
    vsvc = 0.d0
    vshc = 0.d0

    vpc = 0.d0
    vsc = 0.d0
    rhoc = 0.d0

    ! isotropic by default
    etac = 1.d0

    ! anisotropy
    c11c = 0.d0; c12c = 0.d0; c13c = 0.d0
    c14c = 0.d0; c15c = 0.d0; c16c = 0.d0
    c22c = 0.d0; c23c = 0.d0; c24c = 0.d0
    c25c = 0.d0; c26c = 0.d0; c33c = 0.d0
    c34c = 0.d0; c35c = 0.d0; c36c = 0.d0
    c44c = 0.d0; c45c = 0.d0; c46c = 0.d0
    c55c = 0.d0; c56c = 0.d0; c66c = 0.d0

    ! moho depth
    moho = 0.d0

    ! sediment depth
    sediment = 0.d0
  endif

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
      call model_crust_1_0(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
      if (moho_only) return
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_CRUST2)
      ! default
      ! crust 2.0
      call model_crust_2_0(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
      if (moho_only) return
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_CRUSTMAPS)
      ! general crustmaps
      call model_crustmaps(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
      if (moho_only) return
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_EPCRUST)
      ! if defined within lat/lon-range, takes vp/vs/rho/moho from eucrust07
      call model_epcrust(lat,lon,r,vpc_area,vsc_area,rhoc_area,moho_area,sediment_area, &
                         found_crust_area,elem_in_crust,point_in_area)
      if (point_in_area) then
        vpvc = vpc_area
        vphc = vpc_area
        vsvc = vsc_area
        vshc = vsc_area
        rhoc = rhoc_area
        moho = moho_area
        sediment = sediment_area
        found_crust = found_crust_area
      else
        ! by default takes Crust1.0 values
        call model_crust_1_0(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
        vpvc = vpc
        vphc = vpc
        vsvc = vsc
        vshc = vsc
      endif

    case (ICRUST_CRUST_SH)
      ! SH crust: provides TI crust
      call crust_sh(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)

    case (ICRUST_EUCRUST)
      ! by default takes Crust1.0 values for vs/vp/rho/moho
      call model_crust_1_0(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc
      ! if defined within lat/lon-range, takes vp/moho from eucrust07
      call model_eucrust(lat,lon,r,vpc_area,moho_area,sediment_area,found_crust_area,point_in_area)
      if (point_in_area) then
        vpvc = vpc_area
        vphc = vpc_area
        moho = moho_area
        sediment = sediment_area
        found_crust = found_crust_area
      endif

    case (ICRUST_SGLOBECRUST)
      ! modified crust 2.0 for SGLOBE-rani
      call model_sglobecrust(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
      if (moho_only) return
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc

    case (ICRUST_BKMNS_GLAD)
      ! GLAD model expansion on block-mantle-spherical-harmonics
      ! gets moho/sediment depth from Crust2.0 model
      call model_crust_2_0(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)
      if (moho_only) return
      vpvc = vpc
      vphc = vpc
      vsvc = vsc
      vshc = vsc
      ! overimposes model values taken at actual position (with topography) between 80km depth and surface topo
      call model_bkmns_crust(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc)

    case (ICRUST_SPIRAL)
      ! anisotropic crust from SPiRaL
      ! unless we use SPiRaL anisotropic mantle, provdies TI crust only
      call model_crust_spiral(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment, &
                              c11c,c12c,c13c,c14c,c15c,c16c,c22c,c23c,c24c,c25c,c26c, &
                              c33c,c34c,c35c,c36c,c44c,c45c,c46c,c55c,c56c,c66c, &
                              found_crust,elem_in_crust,moho_only)

    case default
      stop 'crustal model type not defined'

  end select


  end subroutine meshfem3D_model_crust

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_models_getatten_val(idoubling,r_prem,theta,phi, &
                                           ispec, i, j, k, &
                                           tau_e,tau_s, &
                                           moho,Qmu,Qkappa,elem_in_crust)

! sets attenuation values tau_e and Qmu for a given point
!
! note:  only Qmu attenuation considered, Qkappa attenuation not used so far in solver...

  use constants, only: N_SLS,SPONGE_MIN_Q,SPONGE_WIDTH_IN_DEGREES

  use meshfem3D_models_par

  use model_prem_par, only: PREM_RMOHO

  use shared_parameters, only: ABSORB_USING_GLOBAL_SPONGE, &
          SPONGE_LATITUDE_IN_DEGREES,SPONGE_LONGITUDE_IN_DEGREES,SPONGE_RADIUS_IN_DEGREES

  implicit none

  integer,intent(in) :: idoubling

  double precision,intent(in) :: r_prem
  double precision,intent(in) :: theta,phi

  integer,intent(in) :: ispec,i,j,k

  double precision,intent(in) :: moho

  ! attenuation values
  double precision,intent(inout) :: Qkappa,Qmu
  double precision, dimension(N_SLS),intent(inout) :: tau_s, tau_e

  logical,intent(in) :: elem_in_crust

  ! local parameters
  double precision :: theta_degrees,phi_degrees
  double precision :: r_used

  ! geographical values
  double precision :: dist, theta_c, phi_c, dist_c, edge, sponge

  ! initializes
  tau_e(:)   = 0.0d0

!---
!
! ADD YOUR MODEL HERE
!
!---

  ! Get the value of Qmu (Attenuation) dependent on
  ! the radius (r_prem) and idoubling flag
  if (ATTENUATION_GLL) then
    ! GLL models with attenuation
    call model_attenuation_gll(ispec, i, j, k, Qmu)

  else if (ATTENUATION_3D) then
    ! used for models: s362ani_3DQ, s362iso_3DQ, 3D_attenuation, SPiRal

    ! gets spherical coordinates in degrees
    theta_degrees = theta / DEGREES_TO_RADIANS
    phi_degrees = phi / DEGREES_TO_RADIANS

    ! in case models incorporate a 3D crust, attenuation values for mantle
    ! get expanded up to surface, and for the crustal points Qmu for PREM crust is imposed
    r_used = r_prem*R_PLANET_KM
    if (CRUSTAL) then
      if (r_prem > (ONE-moho) .or. elem_in_crust) then
        ! points in actual crust: puts point radius into prem crust
        r_used = PREM_RMOHO/1000.d0 * 1.0001
      else if (r_prem*R_PLANET_KM >= PREM_RMOHO/1000.d0) then
        ! points below actual crust (e.g. oceanic crust case), but above prem moho:
        ! puts point slightly below prem moho to expand mantle values at that depth
        r_used = PREM_RMOHO/1000.d0 * 0.99999
      endif
    endif ! CRUSTAL

    ! gets value according to radius/theta/phi location and idoubling flag
    call model_atten3D_QRFSI12(r_used,theta_degrees,phi_degrees,Qmu,idoubling)

  else

    select case (REFERENCE_1D_MODEL)

      ! case (REFERENCE_MODEL_PREM,REFERENCE_MODEL_PREM2)
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
            Qmu = 300.0d0
            Qkappa = 57822.5d0 !  not used so far...
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

  ! sponge layer
  if (ABSORB_USING_GLOBAL_SPONGE) then
    ! get distance to chunk center
    call lat_2_geocentric_colat_dble(SPONGE_LATITUDE_IN_DEGREES, theta_c)
    phi_c = SPONGE_LONGITUDE_IN_DEGREES * DEGREES_TO_RADIANS
    call reduce(theta_c, phi_c)

    dist = acos(cos(theta)*cos(theta_c) + sin(theta)*sin(theta_c)*cos(phi-phi_c))
    dist_c = SPONGE_RADIUS_IN_DEGREES * DEGREES_TO_RADIANS
    edge = SPONGE_WIDTH_IN_DEGREES * DEGREES_TO_RADIANS

    if (dist > dist_c .and. Qmu > SPONGE_MIN_Q) then
      if (dist - dist_c < edge) then
        sponge = (1 + cos((dist - dist_c) / edge * PI)) / 2
      else
        sponge = 0.d0
      endif
      Qmu = SPONGE_MIN_Q + (Qmu - SPONGE_MIN_Q) * sponge
    endif
  endif

  ! Get tau_e from tau_s and Qmu
  call model_attenuation_getstored_tau(Qmu, tau_s, tau_e)

  end subroutine meshfem3D_models_getatten_val


!
!-------------------------------------------------------------------------------------------------
!


  subroutine meshfem3D_models_impose_val(iregion_code,r,theta,phi,ispec,i,j,k, &
                                         vpv,vph,vsv,vsh,rho,eta_aniso, &
                                         c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                         c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

! overwrites values with updated model values (from iteration step) here, given at all GLL points

  use meshfem3D_models_par

  implicit none

  integer,intent(in) :: iregion_code,ispec,i,j,k
  double precision,intent(in) :: r,theta,phi

  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,eta_aniso
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66


  ! checks if anything to do
  if (.not. MODEL_GLL) return

  ! over-impose values from model GLL values
  call model_gll_impose_val(iregion_code,r,theta,phi,ispec,i,j,k, &
                            vpv,vph,vsv,vsh,rho,eta_aniso, &
                            c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                            c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine meshfem3D_models_impose_val

!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_plot_VTK_crust_moho()

! creates VTK output file for moho depths spanning full globe

  use constants, only: PI,IREGION_CRUST_MANTLE,MAX_STRING_LEN,IMAIN,myrank
  use shared_parameters, only: LOCAL_PATH,R_PLANET_KM

  implicit none

  ! local parameters
  double precision :: lat,lon,r,phi,theta
  double precision :: xmesh,ymesh,zmesh

  double precision :: vpv,vph,vsv,vsh,rho,eta_aniso
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                      c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision :: moho,sediment
  integer :: iregion_code
  logical :: elem_in_crust

  ! 2D grid
  integer, parameter :: NLAT = 180,NLON = 360
  double precision :: dlat,dlon

  integer :: i,j,ier,iglob,ispec
  integer :: icorner1,icorner2,icorner3,icorner4
  integer :: nspec,nglob
  double precision, dimension(:), allocatable :: moho_depth,sediment_depth
  double precision, dimension(:), allocatable :: tmp_x,tmp_y,tmp_z
  integer, dimension(:,:), allocatable :: ibool2D

  character(len=MAX_STRING_LEN) :: filename

  ! only main process writes file
  if (myrank /= 0) return

  write(IMAIN,*)
  write(IMAIN,*) '  VTK moho output: resolution = ',1,'deg'
  write(IMAIN,*) '                   NLAT = ',NLAT
  write(IMAIN,*) '                   NLON = ',NLON

  elem_in_crust = .true.
  iregion_code = IREGION_CRUST_MANTLE

  dlat = 180.d0/NLAT
  dlon = 360.d0/NLON

  nglob = NLON * NLAT
  allocate(moho_depth(nglob), &
           sediment_depth(nglob), &
           tmp_x(nglob),tmp_y(nglob),tmp_z(nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating moho_depth arrays'
  moho_depth(:) = 0.d0
  sediment_depth(:) = 0.d0
  tmp_x(:) = 0.d0
  tmp_y(:) = 0.d0
  tmp_z(:) = 0.d0

  ! loop in 1-degree steps over the globe
  do j = 1,NLAT
    do i = 1,NLON
      ! lat/lon in degrees (range lat/lon = [-90,90] / [-180,180]
      lat = 90.d0 - j*dlat + 0.5d0
      lon = -180.d0 + i*dlon - 0.5d0

      ! converts to colatitude theta/phi in radians
      theta = (90.d0 - lat) * PI/180.d0   ! colatitude between [0,pi]
      phi = lon * PI/180.d0               ! longitude between [-pi,pi]
      r = 1.0d0                           ! radius at surface (normalized)

      ! theta to [0,PI] and phi to [0,2PI]
      call reduce(theta,phi)

      ! gets moho
      call meshfem3D_models_get3Dcrust_val(iregion_code,r,theta,phi, &
                                           vpv,vph,vsv,vsh,rho,eta_aniso, &
                                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                                           elem_in_crust,moho,sediment)
      ! stores
      iglob = i + (j-1) * NLON
      moho_depth(iglob) = moho * R_PLANET_KM  ! dimensionalize moho depth to km
      sediment_depth(iglob) = sediment * R_PLANET_KM

      ! gets point's position x/y/z
      call rthetaphi_2_xyz_dble(xmesh,ymesh,zmesh,r,theta,phi)
      tmp_x(iglob) = xmesh
      tmp_y(iglob) = ymesh
      tmp_z(iglob) = zmesh

      ! debug
      !print *,'debug: lat/lon',lat,lon,theta,phi,'xyz',xmesh,ymesh,zmesh
    enddo
  enddo

  ! creates an ibool2D array for 2D quad mesh
  nspec = (NLAT-1) * (NLON-1)
  allocate(ibool2D(4,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating moho_depth ibool2D array'
  ibool2D(:,:) = 0

  ispec = 0
  do j = 1,NLAT-1
    do i = 1,NLON-1
      !     1       2        NLON
      !     * ----  * -//-- *
      !     | ispec |       |
      !     * ----  * -//-- *
      !     3       4

      ispec = ispec + 1
      icorner1 = i   + (j  -1) * NLON
      icorner2 = i+1 + (j  -1) * NLON
      icorner3 = i   + (j+1-1) * NLON
      icorner4 = i+1 + (j+1-1) * NLON

      ibool2D(1,ispec) = icorner1
      ibool2D(2,ispec) = icorner2
      ibool2D(3,ispec) = icorner3
      ibool2D(4,ispec) = icorner4
    enddo
  enddo
  if (ispec /= nspec) stop 'Error invalid ispec value in plot_crustal_model_moho()'

  filename = trim(LOCAL_PATH)//'/mesh_moho_depth'

  call write_VTK_2Ddata_dp(nspec,nglob,tmp_x,tmp_y,tmp_z,ibool2D,moho_depth,filename)

  write(IMAIN,*)
  write(IMAIN,*) '  moho depths written to file: ',trim(filename)//'.vtk'
  write(IMAIN,*) '  min/max = ',sngl(minval(moho_depth(:))),'/',sngl(maxval(moho_depth(:))),'(km)'
  write(IMAIN,*)

  filename = trim(LOCAL_PATH)//'/mesh_sediment_depth'

  call write_VTK_2Ddata_dp(nspec,nglob,tmp_x,tmp_y,tmp_z,ibool2D,sediment_depth,filename)

  write(IMAIN,*) '  sediment depths written to file: ',trim(filename)//'.vtk'
  write(IMAIN,*) '  min/max = ',sngl(minval(sediment_depth(:))),'/',sngl(maxval(sediment_depth(:))),'(km)'
  write(IMAIN,*)

  ! frees memory
  deallocate(moho_depth,sediment_depth,tmp_x,tmp_y,tmp_z,ibool2D)

  end subroutine meshfem3D_plot_VTK_crust_moho


!
!-------------------------------------------------------------------------------------------------
!

  subroutine meshfem3D_plot_VTK_topo_bathy()

! creates VTK output file for moho depths spanning full globe

  use constants, only: PI,MAX_STRING_LEN,IMAIN,myrank
  use shared_parameters, only: LOCAL_PATH,RESOLUTION_TOPO_FILE,R_PLANET_KM
  use meshfem3D_models_par, only: ibathy_topo

  implicit none

  ! local parameters
  double precision :: lat,lon,r,phi,theta,elevation
  double precision :: xmesh,ymesh,zmesh

  ! 2D grid
  integer, parameter :: NLAT_g = 180,NLON_g = 360  ! full sphere
  integer :: NLAT,NLON
  double precision :: dlat,dlon,samples_per_degree,dist

  integer :: i,j,ier,iglob,ispec
  integer :: icorner1,icorner2,icorner3,icorner4
  integer :: nspec,nglob
  double precision, dimension(:), allocatable :: topo,tmp_x,tmp_y,tmp_z
  integer, dimension(:,:), allocatable :: ibool2D

  character(len=MAX_STRING_LEN) :: filename

  ! only main process writes file
  if (myrank /= 0) return

  ! number of samples per degree
  samples_per_degree =  60.d0 / RESOLUTION_TOPO_FILE
  ! distance between sampling points
  dist = 2.d0 * PI * R_PLANET_KM / 360.d0 / samples_per_degree

  ! total nlat/nlon
  NLAT = int(NLAT_g * samples_per_degree)
  NLON = int(NLON_g * samples_per_degree)

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) '  VTK topo output: topo resolution in minutes = ',sngl(RESOLUTION_TOPO_FILE)
  write(IMAIN,*) '                   samples per degree         = ',sngl(samples_per_degree)
  write(IMAIN,*) '                   resolution distance        = ',sngl(dist),'(km)'
  write(IMAIN,*) '                   full globe NLAT = ',NLAT
  write(IMAIN,*) '                              NLON = ',NLON
  write(IMAIN,*) '                              total number of points NLAT x NLON = ',NLAT*NLON
  call flush_IMAIN()

  ! limits size of output file
  if (samples_per_degree > 2) then
    samples_per_degree = 2.d0
    NLAT = int(NLAT_g * samples_per_degree)
    NLON = int(NLON_g * samples_per_degree)
    ! info
    write(IMAIN,*) '                   limiting output to samples per degree         = ',int(samples_per_degree)
    call flush_IMAIN()
  endif

  ! integer 32-bit limitation (integer overflow)
  if (NLAT > int(2147483646.d0 / dble(NLON))) then
    print *,'Error: integer 32-bit overflow: nlat x nlon = ',NLAT,'x',NLON,' exceeds limit ',2147483646
    stop 'Error vtu (nlat x nlon) might exceed integer limit'
  endif

  nglob = NLON * NLAT
  allocate(topo(nglob), &
           tmp_x(nglob),tmp_y(nglob),tmp_z(nglob),stat=ier)
  if (ier /= 0) stop 'Error allocating topo_bathy arrays'
  topo(:) = 0.d0
  tmp_x(:) = 0.d0
  tmp_y(:) = 0.d0
  tmp_z(:) = 0.d0

  dlat = 180.d0/NLAT
  dlon = 360.d0/NLON

  ! loop in 1-degree steps over the globe
  do j = 1,NLAT
    do i = 1,NLON
      ! lat/lon in degrees (range lat/lon = [-90,90] / [-180,180]
      lat = 90.d0 - j*dlat + 0.5d0*dlat
      lon = -180.d0 + i*dlon - 0.5d0*dlon

      ! converts to colatitude theta/phi in radians
      theta = (90.d0 - lat) * PI/180.d0   ! colatitude between [0,pi]
      phi = lon * PI/180.d0               ! longitude between [-pi,pi]
      r = 1.0d0                           ! radius at surface (normalized)

      ! gets point's position
      call rthetaphi_2_xyz_dble(xmesh,ymesh,zmesh,r,theta,phi)

      ! debug
      !print *,'debug: lat/lon',lat,lon,theta,phi,'xyz',xmesh,ymesh,zmesh

      ! compute elevation at current point
      call get_topo_bathy(lat,lon,elevation,ibathy_topo)

      ! stores
      iglob = i + (j-1) * NLON
      topo(iglob) = elevation / 1000.d0  ! convert to km, original elevation given in m
      tmp_x(iglob) = xmesh
      tmp_y(iglob) = ymesh
      tmp_z(iglob) = zmesh
    enddo
  enddo

  ! creates an ibool2D array for 2D quad mesh
  nspec = (NLAT-1) * (NLON-1)
  allocate(ibool2D(4,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating topo_bathy ibool2D array'
  ibool2D(:,:) = 0

  ispec = 0
  do j = 1,NLAT-1
    do i = 1,NLON-1
      !     1       2        NLON
      !     * ----  * -//-- *
      !     | ispec |       |
      !     * ----  * -//-- *
      !     3       4

      ispec = ispec + 1
      icorner1 = i   + (j  -1) * NLON
      icorner2 = i+1 + (j  -1) * NLON
      icorner3 = i   + (j+1-1) * NLON
      icorner4 = i+1 + (j+1-1) * NLON

      ibool2D(1,ispec) = icorner1
      ibool2D(2,ispec) = icorner2
      ibool2D(3,ispec) = icorner3
      ibool2D(4,ispec) = icorner4
    enddo
  enddo
  if (ispec /= nspec) stop 'Error invalid ispec value in plot_VTK_topo_bathy()'

  filename = trim(LOCAL_PATH)//'/mesh_topo_bathy'

  call write_VTU_2Ddata_dp(nspec,nglob,tmp_x,tmp_y,tmp_z,ibool2D,topo,filename)

  write(IMAIN,*)
  write(IMAIN,*) '  elevations written to file: ',trim(filename)//'.vtk'
  write(IMAIN,*) '  min/max = ',sngl(minval(topo(:))),'/',sngl(maxval(topo(:))),'(km)'
  write(IMAIN,*)


  ! frees memory
  deallocate(topo,tmp_x,tmp_y,tmp_z,ibool2D)

  end subroutine meshfem3D_plot_VTK_topo_bathy


