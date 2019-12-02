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

  subroutine get_model_parameters()

  use constants
  use shared_parameters

  implicit none

  ! turns on/off corresponding 1-D/3-D model flags
  call get_model_parameters_flags()

  ! sets constants for planet
  call get_model_planet_constants()

  ! sets radius for each discontinuity and ocean density values
  call get_model_parameters_radii()

  end subroutine get_model_parameters


!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_model_parameters_flags()

  use constants

  use shared_parameters, only: MODEL, &
    REFERENCE_1D_MODEL,REFERENCE_CRUSTAL_MODEL, &
    THREE_D_MODEL,THREE_D_MODEL_IC, &
    MODEL_GLL,MODEL_GLL_TYPE, &
    ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
    ATTENUATION_GLL, CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE, &
    HONOR_1D_SPHERICAL_MOHO, MODEL_3D_MANTLE_PERTUBATIONS, &
    ONE_CRUST, TRANSVERSE_ISOTROPY, OCEANS,TOPOGRAPHY, &
    CEM_REQUEST,CEM_ACCEPT,GPU_MODE

  implicit none

  ! local parameters
  character(len=64) :: ending
  character(len=MAX_STRING_LEN) :: MODEL_ROOT,MODEL_L
  integer :: impose_crust
  integer :: irange,i

  ! defaults:
  !
  ! HONOR_1D_SPHERICAL_MOHO: honor PREM Moho or not: doing so drastically reduces
  ! the stability condition and therefore the time step, resulting in expensive
  ! calculations. If not, honor a fictitious Moho at the depth of 40 km
  ! in order to have even radial sampling from the d220 to the Earth surface.
  !
  ! ONE_CRUST: in order to increase stability and therefore to allow cheaper
  ! simulations (larger time step), 1D models can be run with just one average crustal
  ! layer instead of two.
  !
  ! CASE_3D : this flag allows the stretching of the elements in the crustal
  ! layers in the case of 3D models. The purpose of this stretching is to squeeze more
  ! GLL points per km in the upper part of the crust than in the lower part.
  !
  ! CRUSTAL : flag set to .true. if a 3D crustal model (e.g. Crust-2.0) will be used or
  !  to .false. for a 1D crustal model.

  ! converts all string characters to lowercase (to make user input case-insensitive)
  MODEL_L = MODEL
  irange = iachar('a') - iachar('A')
  do i = 1,len_trim(MODEL_L)
    if (lge(MODEL_L(i:i),'A') .and. lle(MODEL_L(i:i),'Z')) then
      MODEL_L(i:i) = achar(iachar(MODEL_L(i:i)) + irange)
    endif
  enddo
  MODEL_ROOT = MODEL_L ! sets root name of model to original one

  ! inner core anisotropy
  ! initializes inner core parameters
  ANISOTROPIC_INNER_CORE = .false.
  THREE_D_MODEL_IC = 0
  ! extract ending of model name
  ending = ''
  if (len_trim(MODEL_ROOT) > 4 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-3:len_trim(MODEL_ROOT))
  ! determines if the anisotropic inner core option should be turned on
  if (trim(ending) == '_aic') then
    ANISOTROPIC_INNER_CORE = .true.
    THREE_D_MODEL_IC = THREE_D_MODEL_INNER_CORE_ISHII ! since we only have a single inner core aniso model, assumes we take it
    ! in case it has an ending for the inner core, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1:len_trim(MODEL_ROOT)-4)
  endif

  ! crust/mantle anisotropy
  ! uses no full anisotropy calculations by default
  ANISOTROPIC_3D_MANTLE = .false.
  ! determines if the anisotropic mantle option should be turned on
  ending = ''
  if (len_trim(MODEL_ROOT) > 4 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-3:len_trim(MODEL_ROOT))
  if (trim(ending) == '_acm') then
    ! turn handling on of crust/mantle elements for fully anisotropic calculations
    ANISOTROPIC_3D_MANTLE = .true.
    ! in case it has an ending "_aniso", remove it from the name
    MODEL_ROOT = MODEL_ROOT(1:len_trim(MODEL_ROOT)-4)
  endif

  ! crustal option by name ending
  impose_crust = 0
  ending = ''
  ! 1D crust options
  ! checks with '_onecrust' option
  if (len_trim(MODEL_ROOT) > 9 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-8:len_trim(MODEL_ROOT))
  if (trim(ending) == '_onecrust') then
    impose_crust = -1 ! negative flag for 1Dcrust option with 1-layer
    ! in case it has an ending for the 1D crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-9)
  endif
  ! checks with '_1Dcrust' option
  if (len_trim(MODEL_ROOT) > 8 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-7:len_trim(MODEL_ROOT))
  if (trim(ending) == '_1dcrust') then
    impose_crust = -2 ! negative flag for 1Dcrust option with 2-layers
    ! in case it has an ending for the 1D crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-8)
  endif

  ! 3D crust options
  ! checks with '_crustmaps' option
  if (len_trim(MODEL_ROOT) > 10 ) &
    ending = MODEL_ROOT(len_trim(MODEL_ROOT)-9:len_trim(MODEL_ROOT))
  if (trim(ending) == '_crustmaps') then
    impose_crust = ICRUST_CRUSTMAPS
    ! in case it has an ending for the crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-10)
  endif
  ! checks with '_crust1.0' option
  if (len_trim(MODEL_ROOT) > 9 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-8:len_trim(MODEL_ROOT))
  if (trim(ending) == '_crust1.0') then
    impose_crust = ICRUST_CRUST1
    ! in case it has an ending for the crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-9)
  endif
  ! checks with '_crust2.0' option
  if (len_trim(MODEL_ROOT) > 9 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-8:len_trim(MODEL_ROOT))
  if (trim(ending) == '_crust2.0') then
    impose_crust = ICRUST_CRUST2
    ! in case it has an ending for the crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-9)
  endif
  ! checks with '_epcrust' option
  if (len_trim(MODEL_ROOT) > 8 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-7:len_trim(MODEL_ROOT))
  if (trim(ending) == '_epcrust') then
    impose_crust = ICRUST_EPCRUST
    ! in case it has an ending for the crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-8)
  endif
  ! checks with '_eucrust' option
  if (len_trim(MODEL_ROOT) > 8 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-7:len_trim(MODEL_ROOT))
  if (trim(ending) == '_eucrust') then
    impose_crust = ICRUST_EUCRUST
    ! in case it has an ending for the crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-8)
  endif
  ! checks with '_crustSH' option
  if (len_trim(MODEL_ROOT) > 8 ) ending = MODEL_ROOT(len_trim(MODEL_ROOT)-7:len_trim(MODEL_ROOT))
  if (trim(ending) == '_crustsh') then
    impose_crust = ICRUST_CRUST_SH
    ! in case it has an ending for the crust, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL_ROOT)-8)
  endif


!---
!
! ADD YOUR MODEL HERE
!
!---

  ! default values

  ! uses PREM (isotropic) as the 1D reference model by default
  REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM

  ! default crustal model
  ! (used Crust2.0 as default when CRUSTAL flag is set for simulation)
  REFERENCE_CRUSTAL_MODEL = ICRUST_CRUST2

  ! uses 1D attenuation model by default
  ATTENUATION_3D = .false.

  ! no crustal mesh stretching and 3D crust models by default
  CASE_3D = .false.
  CRUSTAL = .false.

  ! by default, crust will be split into upper and lower crust using 2 element layers
  ! (if set to true, this uses 1 element layer only for the crust and assign a single crustal material)
  ONE_CRUST = .false.

  ! uses no 3D heterogeneity mantle by default
  HETEROGEN_3D_MANTLE = .false.
  MODEL_3D_MANTLE_PERTUBATIONS = .false.
  HONOR_1D_SPHERICAL_MOHO = .false.
  MODEL_GLL = .false.
  MODEL_GLL_TYPE = 0

  ! no CEM by default
  CEM_REQUEST = .false.
  CEM_ACCEPT  = .false.

  ! no 3D model by default
  THREE_D_MODEL = 0
  TRANSVERSE_ISOTROPY = .false.

  ! model specifics
  select case (trim(MODEL_ROOT))

  ! 1-D models
  case ('1d_isotropic_prem')
    HONOR_1D_SPHERICAL_MOHO = .true.

  case ('1d_transversely_isotropic_prem')
    HONOR_1D_SPHERICAL_MOHO = .true.
    TRANSVERSE_ISOTROPY = .true.

  case ('1d_iasp91','1d_1066a','1d_ak135f_no_mud','1d_jp3d','1d_sea99')
    HONOR_1D_SPHERICAL_MOHO = .true.
    if (MODEL_ROOT == '1d_iasp91') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if (MODEL_ROOT == '1d_1066a') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if (MODEL_ROOT == '1d_ak135f_no_mud') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135F_NO_MUD
    else if (MODEL_ROOT == '1d_jp3d') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_JP1D
    else if (MODEL_ROOT == '1d_sea99') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    else
      stop 'reference 1D Earth model unknown'
    endif

  case ('1d_ref')
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    TRANSVERSE_ISOTROPY = .true.

  case ('1d_ref_iso')
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF

  case ('1d_sohl')
    ! Mars
    TRANSVERSE_ISOTROPY = .false. ! enforces isotropic model
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SOHL

  case('1d_sohl_3d_crust')
    ! Mars
    TRANSVERSE_ISOTROPY = .false. ! enforces isotropic model
    CASE_3D = .true.
    CRUSTAL = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SOHL
    REFERENCE_CRUSTAL_MODEL = ICRUST_CRUSTMAPS

  ! 3-D models
  case ('transversely_isotropic_prem_plus_3d_crust_2.0')
    CASE_3D = .true.
    CRUSTAL = .true.
    ONE_CRUST = .true.
    TRANSVERSE_ISOTROPY = .true.

  case ('transversely_isotropic_prem_plus_3d_crust_1.0')
    CASE_3D = .true.
    CRUSTAL = .true.
    ONE_CRUST = .true.
    TRANSVERSE_ISOTROPY = .true.
    impose_crust = ICRUST_CRUST1

  case ('transversely_isotropic_prem_plus_3d')
    CASE_3D = .true.
    CRUSTAL = .true.
    ONE_CRUST = .true.
    TRANSVERSE_ISOTROPY = .true.

  case ('s20rts')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_S20RTS
    TRANSVERSE_ISOTROPY = .true.

  case ('s40rts')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_S40RTS
    TRANSVERSE_ISOTROPY = .true.

  case ('full_sh')
    ! uses PREM by default, uncomment if necessary
    !!REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    ! uses SH crustal model by default
    REFERENCE_CRUSTAL_MODEL = ICRUST_CRUST_SH
    !ATTENUATION_3D = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_MANTLE_SH
    TRANSVERSE_ISOTROPY = .true.

  case ('sea99_jp3d1994')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    THREE_D_MODEL = THREE_D_MODEL_SEA99_JP3D

  case ('sea99')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    THREE_D_MODEL = THREE_D_MODEL_SEA99

  case ('jp3d1994')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_JP1D
    THREE_D_MODEL = THREE_D_MODEL_JP3D

  case ('s362ani')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.

  case ('s362iso')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  case ('s362wmani')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362WMANI
    TRANSVERSE_ISOTROPY = .true.

  case ('s362ani_prem')
    CASE_3D = .true.
    CRUSTAL = .true.
    TRANSVERSE_ISOTROPY = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_S362ANI_PREM

  case ('s362ani_3dq')
    ATTENUATION_3D = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.

  case ('s362ani_azi')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.    ! to include original tiso perturbations from reference model
    ANISOTROPIC_3D_MANTLE = .true.  ! to impose mantle elements as fully anisotropic in solver

  case ('s362iso_3dq')
    ATTENUATION_3D = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  case ('s29ea')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S29EA
    TRANSVERSE_ISOTROPY = .true.

  case ('sgloberani_aniso')
    ! anisotropic perturbations to PREM
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_SGLOBE
    TRANSVERSE_ISOTROPY = .true.

  case ('sgloberani_iso')
    ! isotropic perturbations to PREM
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_SGLOBE_ISO
    TRANSVERSE_ISOTROPY = .true. ! still based on transversely isotropic PREM

  case ('3d_attenuation')
    ! uses default model (PREM isotropic) as reference velocity model
    ! adds Daltons' 3D attenuation model
    ATTENUATION_3D = .true.
    CASE_3D = .true.
    ONE_CRUST = .true.

  case ('3d_anisotropic')
    ! CRUSTAL = .true. ! with 3D crust: depends on 3D mantle reference model
    CASE_3D = .true. ! crustal moho stretching
    ONE_CRUST = .true. ! 1 element layer in top crust region
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    TRANSVERSE_ISOTROPY = .true. ! to use transverse isotropic PREM 1D ref model
    ! REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135F_NO_MUD
    ! TRANSVERSE_ISOTROPY = .false. ! for AK135 ref model
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    THREE_D_MODEL = THREE_D_MODEL_ANISO_MANTLE
    ANISOTROPIC_3D_MANTLE = .true. ! treats mantle elements as fully anisotropic

  case ('heterogen')
    ONE_CRUST = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    HETEROGEN_3D_MANTLE = .true.  ! adds additional (dvp,dvs,drho) perturbations on top of reference 3D model
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.

#ifdef USE_CEM
  case ('cem_request')
    CEM_REQUEST         = .true.
    TRANSVERSE_ISOTROPY = .true.

  case ('cem_accept')
    CEM_ACCEPT          = .true.
    TRANSVERSE_ISOTROPY = .true.

  case ('cem_gll')
    MODEL_GLL = .true.
    THREE_D_MODEL = THREE_D_MODEL_GLL
    TRANSVERSE_ISOTROPY = .true.
#else
  case ('cem_request','cem_accept','cem_gll')
    print *,'Error model ',trim(MODEL),': package compiled without CEM model support. Please re-configure with --with-cem support.'
    stop 'Invalid CEM model requested, compiled without CEM support'
#endif

  case ('ppm')
    ! superimposed based on isotropic-prem
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_PPM
    TRANSVERSE_ISOTROPY = .true. ! to use transverse-isotropic prem

  case ('gll','gll_tiso')
    ! default GLL model:
    ! will be given on local basis, at all GLL points,
    ! as from meshfem3d output from routine save_arrays_solver()
    !
    ! transverse_isotropic should match the one from the reference model;
    ! reference model set in constants.h: GLL_REFERENCE_MODEL and GLL_REFERENCE_1D_MODEL
    ONE_CRUST = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    TRANSVERSE_ISOTROPY = .true.  ! same as reference model
    MODEL_GLL = .true.
    MODEL_GLL_TYPE = 2 ! (2 == tiso) input model files are tiso
    REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
    THREE_D_MODEL = THREE_D_MODEL_GLL
    ! note: after call to this routine we will reset
    !       THREE_D_MODEL = THREE_D_MODEL_S29EA
    !       to initialize 3D mesh structure based on the initial 3D model (like 420/660 topography,..)

  case ('gll_iso')
    ! isotropic GLL model
    ONE_CRUST = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    TRANSVERSE_ISOTROPY = .true.  ! same as reference model
    !TRANSVERSE_ISOTROPY = .false.  ! forces earth model to be considered isotropic, not tiso
    MODEL_GLL = .true.
    MODEL_GLL_TYPE = 1 ! (1 == iso) input model files are iso
    REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
    THREE_D_MODEL = THREE_D_MODEL_GLL

  case ('gll_azi')
    ! azimuthal anisotropy GLL model
    ONE_CRUST = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.    ! to read in 3D reference model setup (like s362ani)
    TRANSVERSE_ISOTROPY = .true.    ! to include original tiso perturbations from reference model
    ANISOTROPIC_3D_MANTLE = .true.  ! to impose mantle elements as fully anisotropic in solver
    MODEL_GLL = .true.
    MODEL_GLL_TYPE = 3 ! azimuthal type
    REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
    THREE_D_MODEL = THREE_D_MODEL_GLL

  case ('gll_qmu')
    ! default GLL model with Qmu values:
    ONE_CRUST = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    TRANSVERSE_ISOTROPY = .true.  ! same as reference model
    MODEL_GLL = .true.
    MODEL_GLL_TYPE = 2 ! (2 == tiso) input model files are tiso
    REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
    THREE_D_MODEL = THREE_D_MODEL_GLL
    ATTENUATION_GLL = .true.

  case ('gll_mars')
    ! Mars
    ! default GLL model for mars
    CASE_3D = .true.
    CRUSTAL = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SOHL
    REFERENCE_CRUSTAL_MODEL = ICRUST_CRUSTMAPS
    TRANSVERSE_ISOTROPY = .false. ! enforces isotropic model
    MODEL_3D_MANTLE_PERTUBATIONS = .false. ! not based on a 3D mantle model, but 1D model Sohl
    THREE_D_MODEL = 0
    MODEL_GLL = .true.
    MODEL_GLL_TYPE = 1 ! (1 == iso) input model files are iso (vp,vs,rho)

  case ('gapp2')
    CASE_3D = .true.
    CRUSTAL = .true.
    MODEL_3D_MANTLE_PERTUBATIONS = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    THREE_D_MODEL = THREE_D_MODEL_GAPP2
    TRANSVERSE_ISOTROPY = .true.

  case ('ishii')
    ! Ishii et al. (2002) inner core model
    THREE_D_MODEL_IC = THREE_D_MODEL_INNER_CORE_ISHII
    ! takes tiso PREM mantle model as reference
    TRANSVERSE_ISOTROPY = .true.

  case default
    print *
    print *,'Error model: ',trim(MODEL)
    stop 'model not implemented yet, edit get_model_parameters.f90, or ensure you have run ./configure correctly, and recompile'
  end select

  ! additional option for 3D mantle models:
  ! this takes crust from reference 1D model rather than a 3D crust;
  if (impose_crust < 0) then
    ! imposes 1Dcrust options
    ! no 3D crust
    CRUSTAL = .false.
    ! no crustal moho stretching
    CASE_3D = .false.
    ! mesh honors the 1D moho depth
    HONOR_1D_SPHERICAL_MOHO = .true.
    select case (impose_crust)
    case (-1)
      ! onecrust option, single layer
      ONE_CRUST = .true.
    case (-2)
      ! 1Dcrust option, default 2-element layers in top crust region (rather than just one)
      ONE_CRUST = .false.
    case default
      stop 'Invalid impose_crust value for 1D crust'
    end select
  endif

  ! 3D crustal maps
  if (impose_crust > 0) then
    ! sets 3D crustal model
    REFERENCE_CRUSTAL_MODEL = impose_crust
    ! 3D crustal structure flags
    HONOR_1D_SPHERICAL_MOHO = .false.
    CRUSTAL = .true.
    CASE_3D = .true.
    ONE_CRUST = .true.
  endif

  ! suppress the crustal layers (only mantle values, extended to surface)
  if (SUPPRESS_CRUSTAL_MESH) then
    CRUSTAL = .false.
    ONE_CRUST = .false.
    OCEANS = .false.
    TOPOGRAPHY = .false.
  endif

  ! checks flag consistency for crust
  if (HONOR_1D_SPHERICAL_MOHO .and. CRUSTAL ) &
    stop 'honor 1D spherical moho excludes having 3D crustal structure'

  if (HONOR_1D_SPHERICAL_MOHO .and. CASE_3D ) &
    stop 'honor 1D spherical moho excludes having 3D crustal mesh stretching'

  ! checks that IASP91, AK135, 1066A, JP1D or SEA1D is isotropic
  if ((REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91 .or. &
       REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD .or. &
       REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A .or. &
       REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D .or. &
       REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) .and. TRANSVERSE_ISOTROPY) &
        stop 'models IASP91, AK135, 1066A, JP1D and SEA1D are currently isotropic'

  ! Mars
  ! Mars 1D_Sohl is isotropic
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SOHL .and. TRANSVERSE_ISOTROPY) &
      stop 'model 1D_Sohl is currently isotropic'
  ! Mars has no ocean
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SOHL .and. OCEANS) &
    stop 'model 1D_Sohl cannot use an ocean approximation'
  ! Mars not implemented yet on GPU (missing correct gravity)
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SOHL .and. GPU_MODE) &
    stop 'model 1D_Sohl cannot use GPU_MODE'

  end subroutine get_model_parameters_flags

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_model_planet_constants()

  use constants
  use shared_parameters

  implicit none

! note: Please make sure to broadcast the values below in broadcast_computed_parameters.f90
!
!       here, only the master process is setting the new defaults.
!       thus, they need to be broadcast to all other processes.
!       this will be done in routine broadcast_computed_parameters().

  ! sets Planet constants
  select case(REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_SOHL)
    ! Mars
    ! sets planet
    PLANET_TYPE = IPLANET_MARS
    ! radius
    R_EARTH = MARS_R
    R_EARTH_KM = MARS_R_KM
    ! average density
    RHOAV = MARS_RHOAV
    ! gravity
    STANDARD_GRAVITY = MARS_STANDARD_GRAVITY
    ! flattening/eccentricity
    ONE_MINUS_F_SQUARED = MARS_ONE_MINUS_F_SQUARED
    ! topo
    PATHNAME_TOPO_FILE = MARS_PATHNAME_TOPO_FILE
    RESOLUTION_TOPO_FILE = MARS_RESOLUTION_TOPO_FILE
    NX_BATHY = MARS_NX_BATHY
    NY_BATHY = MARS_NY_BATHY
    TOPO_MINIMUM = MARS_TOPO_MINIMUM
    TOPO_MAXIMUM = MARS_TOPO_MAXIMUM
    ! crust
    R_DEEPEST_CRUST = MARS_R_DEEPEST_CRUST
    ! rotation
    HOURS_PER_DAY = MARS_HOURS_PER_DAY
    SECONDS_PER_HOUR = MARS_SECONDS_PER_HOUR
    ! mesh
    MAX_RATIO_CRUST_STRETCHING = MARS_MAX_RATIO_CRUST_STRETCHING
    RMOHO_STRETCH_ADJUSTMENT = MARS_RMOHO_STRETCH_ADJUSTMENT
    R80_STRETCH_ADJUSTMENT = MARS_R80_STRETCH_ADJUSTMENT
    REGIONAL_MOHO_MESH = MARS_REGIONAL_MOHO_MESH
    HONOR_DEEP_MOHO = MARS_HONOR_DEEP_MOHO

  case default
    ! Earth
    ! default: sets Earth as default for R_EARTH, RHOAV, ..
    PLANET_TYPE = IPLANET_EARTH
    ! radius
    R_EARTH = EARTH_R
    R_EARTH_KM = EARTH_R_KM
    ! average density
    RHOAV = EARTH_RHOAV
    ! gravity
    STANDARD_GRAVITY = EARTH_STANDARD_GRAVITY
    ! flattening/eccentricity
    ONE_MINUS_F_SQUARED = EARTH_ONE_MINUS_F_SQUARED
    ! topo
    PATHNAME_TOPO_FILE = EARTH_PATHNAME_TOPO_FILE
    RESOLUTION_TOPO_FILE = EARTH_RESOLUTION_TOPO_FILE
    NX_BATHY = EARTH_NX_BATHY
    NY_BATHY = EARTH_NY_BATHY
    TOPO_MINIMUM = EARTH_TOPO_MINIMUM
    TOPO_MAXIMUM = EARTH_TOPO_MAXIMUM
    ! crust
    R_DEEPEST_CRUST = EARTH_R_DEEPEST_CRUST
    ! rotation
    HOURS_PER_DAY = EARTH_HOURS_PER_DAY
    SECONDS_PER_HOUR = EARTH_SECONDS_PER_HOUR
    ! mesh
    MAX_RATIO_CRUST_STRETCHING = EARTH_MAX_RATIO_CRUST_STRETCHING
    RMOHO_STRETCH_ADJUSTMENT = EARTH_RMOHO_STRETCH_ADJUSTMENT
    R80_STRETCH_ADJUSTMENT = EARTH_R80_STRETCH_ADJUSTMENT
    REGIONAL_MOHO_MESH = EARTH_REGIONAL_MOHO_MESH
    HONOR_DEEP_MOHO = EARTH_HONOR_DEEP_MOHO
  end select

  end subroutine get_model_planet_constants

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_model_parameters_radii()

  use constants

  use shared_parameters, only: &
    ROCEAN,RMIDDLE_CRUST, &
    RMOHO,R80,R120,R220,R400,R600,R670,R771, &
    RTOPDDOUBLEPRIME,RCMB,RICB, &
    RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER, &
    R80_STRETCH_ADJUSTMENT,RMOHO_STRETCH_ADJUSTMENT, &
    RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
    R_EARTH

  use shared_parameters, only: &
    HONOR_1D_SPHERICAL_MOHO,CASE_3D,CRUSTAL,REFERENCE_1D_MODEL

  implicit none


! sets radii in PREM or IASP91 and normalized density at fluid-solid interface on fluid size for coupling
!
! ROCEAN: radius of the ocean (m)
! RMIDDLE_CRUST: radius of the middle crust (m)
! RMOHO: radius of the Moho (m)
! R80: radius of 80 km discontinuity (m)
! R120: radius of 120 km discontinuity (m) in IASP91
! R220: radius of 220 km discontinuity (m)
! R400: radius of 400 km discontinuity (m)
! R600: radius of 600 km 2nd order discontinuity (m)
! R670: radius of 670 km discontinuity (m)
! R771: radius of 771 km 2nd order discontinuity (m)
! RTOPDDOUBLEPRIME: radius of top of D" 2nd order discontinuity (m)
! RCMB: radius of CMB (m)
! RICB: radius of ICB (m)


!---
!
! ADD YOUR MODEL HERE
!
!---

  ! default: PREM
  ROCEAN = 6368000.d0         ! at 3km depth
  RMIDDLE_CRUST = 6356000.d0  ! at 15km depth
  RMOHO = 6346600.d0          ! at 24.4km depth
  R80  = 6291000.d0
  R120 = -1.d0                ! by default there is no d120 discontinuity, except in IASP91, therefore set to fictitious value
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0           ! at 670km depth
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

  ! density ocean
  RHO_OCEANS = 1020.0 / RHOAV   ! value common to all models
  ! densities fluid outer core (PREM)
  RHO_TOP_OC = 9903.4384 / RHOAV
  RHO_BOTTOM_OC = 12166.5885 / RHOAV

  ! differing 1-D model radii
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
    ! IASP91
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO = 6336000.d0          ! at 35km depth
    R80  = 6291000.d0
    R120 = 6251000.d0
    R220 = 6161000.d0
    R400 = 5961000.d0
    ! there is no d600 discontinuity in IASP91 therefore this value is useless
    ! but it needs to be there for compatibility with other subroutines
    R600 = R_EARTH - 600000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0

    RHO_TOP_OC = 9900.2379 / RHOAV
    RHO_BOTTOM_OC = 12168.6383 / RHOAV

  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD) then

!! DK DK values below entirely checked and fixed by Dimitri Komatitsch in December 2012.

    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO  = 6336000.d0         ! at 35km depth
    R80    = 6293500.d0
    R220   = 6161000.d0
    R400   = 5961000.d0
    R670   = 5711000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB   = 3479500.d0
    RICB   = 1217500.d0

    ! values for AK135F that are not discontinuities
    R600 = 5771000.d0
    R771 = 5600000.d0

    RHO_TOP_OC = 9914.5000 / RHOAV
    RHO_BOTTOM_OC = 12139.1000 / RHOAV

  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A) then
    ! values below corrected by Ying Zhou
    ! 1066A
    RMOHO = 6360000.d0 ! at 11km depth
    R400 = 5950000.d0
    R600 = 5781000.d0
    R670 = 5700000.d0
    RCMB = 3484300.d0
    RICB = 1229480.d0

    ! values for 1066A that are not discontinuities
    RTOPDDOUBLEPRIME = 3631000.d0
    R220 = 6161000.d0
    R771 = 5611000.d0
    ! RMIDDLE_CRUST used only for high resolution FFSW1C model, with 3 elements crust simulations
    ! mid_crust = 10 km
    RMIDDLE_CRUST = 6361000.d0
    R80 = 6291000.d0

    ! model 1066A has no oceans, therefore we use the radius of the Earth instead
    ROCEAN = R_EARTH

    RHO_TOP_OC = 9917.4500 / RHOAV
    RHO_BOTTOM_OC = 12160.6500 / RHOAV

  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_1DREF) then
    ! REF
    ROCEAN = 6368000.d0         ! at 3km depth
    RMIDDLE_CRUST = 6356000.d0  ! at 15km
    RMOHO = 6346600.d0          ! at 24.4km depth
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5961000.d0           ! 410km discontinuity
    R600 = 5771000.d0
    R670 = 5721000.d0           ! 650km discontinuity
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3479958.d0
    RICB = 1221491.d0

    RHO_TOP_OC = 9903.48 / RHOAV
    RHO_BOTTOM_OC = 12166.35 / RHOAV

  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D) then
    ! values below corrected by Min Chen
    ! jp1d
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6359000.d0
    RMOHO = 6345000.d0 ! at 26km depth
    R80 = 6291000.d0
    R220 = 6161000.d0
    R400 = 5949000.d0
    R600 = 5781000.d0
    R670 = 5711000.d0
    R771 = 5611000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB = 3482000.d0
    RICB = 1217000.d0

    RHO_TOP_OC = 9900.2379 / RHOAV
    RHO_BOTTOM_OC = 12168.6383 / RHOAV

  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) then
    ! SEA1D without the 2 km of mud layer or the 3km water layer
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6361000.d0
    RMOHO  = 6346000.d0 ! at 25km depth
    R80    = 6291000.d0
    R220   = 6161000.d0
    R400   = 5961000.d0

    R670   = 5711000.d0
    RTOPDDOUBLEPRIME = 3631000.d0
    RCMB   = 3485700.d0
    RICB   = 1217100.d0

    ! values for SEA1D that are not discontinuities
    R600 = 5771000.d0
    R771 = 5611000.d0

  else if ( REFERENCE_1D_MODEL == REFERENCE_MODEL_SOHL) then
    ! Mars
    ! Sohl & Spohn, 1997: Model A, Table 5, pg. 1623
    ROCEAN = 3390000.d0         ! (Physical surface)
    RMIDDLE_CRUST = 3340000.d0  ! 50 km
    RMOHO = 3280000.d0          ! (Crust-mantle boundary)  110 km average crustal thickness
    R80  = 3055500.d0           ! (Rheological lithosphere) 334.5 km deep, too small will cause negative Jacobian err
    R120 = -1.d0                ! by default there is no d120 discontinuity, except in IASP91, therefore set to fictitious value
    R220 = 2908000.d0           ! (Thermal lithosphere) d = 482 km
    R400 = 2655000.d0           ! d = 735 km
    R600 = 2455000.d0           ! d = 935 km
    R670 = 2360000.d0           ! (alpha-olivine-beta-spinel transition) d = 1030 km
    R771 = 2033000.d0           ! (beta-spinel-gamma-spinel transition) d = 1357 km, below which the second doubling implemented
    RTOPDDOUBLEPRIME = 1503000.d0  ! (lower thermal boundary layer) d = 1887 km, too thin for mesh ?
    RCMB = 1468000.d0           ! (Core-mantle boundary) d = 1922 km
    ! note: Sohl & Spohn assume an entirely molten iron alloy core.
    !       since SPECFEM assumes a structure of solid crust/mantle - fluid outer core - solid inner core,
    !       we set a very small inner core radius here.
    RICB = 515000.d0            ! d = 2875 km good for both stability and efficiency

    ! densities
    RHO_OCEANS = 1020.0 / MARS_RHOAV   ! will not be used
    ! densities fluid outer core (from modSOHL)
    RHO_TOP_OC = 6936.40 / MARS_RHOAV
    RHO_BOTTOM_OC = 7268.20 / MARS_RHOAV
  endif

  ! honor the PREM Moho or define a fictitious Moho in order to have even radial sampling
  ! from the d220 to the Earth surface
  if (HONOR_1D_SPHERICAL_MOHO) then
    ! 1D models: all honor their spherical moho
    RMOHO_FICTITIOUS_IN_MESHER = RMOHO
    R80_FICTITIOUS_IN_MESHER = R80
  else
    ! 3D models do not honor PREM moho but a fictitious moho at 40km depth:
    ! either to make simulation cheaper or to have a 3D crustal structure
    RMOHO_FICTITIOUS_IN_MESHER = (R80 + R_EARTH) / 2.0d0
    R80_FICTITIOUS_IN_MESHER = R80
    if (CRUSTAL .and. CASE_3D) then
      !> Hejun
      ! mesh will honor 3D crustal moho topography
      ! moves MOHO up 5km to honor moho topography deeper than 35 km
      ! moves R80 down to 120km depth in order to have less squeezing for elements below moho
      RMOHO_FICTITIOUS_IN_MESHER = RMOHO_FICTITIOUS_IN_MESHER + RMOHO_STRETCH_ADJUSTMENT
      R80_FICTITIOUS_IN_MESHER = R80_FICTITIOUS_IN_MESHER + R80_STRETCH_ADJUSTMENT
    endif
  endif

  ! Mars
  ! Ebru - a quick fix for MARS but should be done and checked again during 3D implementation
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_SOHL) then
    RMOHO_FICTITIOUS_IN_MESHER = RMOHO
    R80_FICTITIOUS_IN_MESHER = R80
  endif

  end subroutine get_model_parameters_radii

