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

  subroutine get_model_parameters(MODEL,REFERENCE_1D_MODEL,THREE_D_MODEL, &
                        ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
                        CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO, &
                        ISOTROPIC_3D_MANTLE,ONE_CRUST,TRANSVERSE_ISOTROPY, &
                        OCEANS,TOPOGRAPHY, &
                        ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R120,R220,R400,R600,R670,R771, &
                        RTOPDDOUBLEPRIME,RCMB,RICB,RMOHO_FICTITIOUS_IN_MESHER, &
                        R80_FICTITIOUS_IN_MESHER,RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                        CEM_REQUEST, CEM_ACCEPT)

  use constants

  implicit none

  character(len=MAX_STRING_LEN) MODEL

  integer REFERENCE_1D_MODEL,THREE_D_MODEL

  logical ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
    CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO,&
    ISOTROPIC_3D_MANTLE,ONE_CRUST,TRANSVERSE_ISOTROPY,CEM_REQUEST,CEM_ACCEPT

  logical OCEANS,TOPOGRAPHY

  double precision ROCEAN,RMIDDLE_CRUST, &
    RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
    RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER

  double precision RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS

  ! turns on/off corresponding 1-D/3-D model flags
  call get_model_parameters_flags(MODEL,REFERENCE_1D_MODEL,THREE_D_MODEL, &
                        ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
                        CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO, &
                        ISOTROPIC_3D_MANTLE,ONE_CRUST,TRANSVERSE_ISOTROPY, &
                        OCEANS,TOPOGRAPHY,CEM_REQUEST,CEM_ACCEPT)

  ! sets radius for each discontinuity and ocean density values
  call get_model_parameters_radii(REFERENCE_1D_MODEL,ROCEAN,RMIDDLE_CRUST, &
                                  RMOHO,R80,R120,R220,R400,R600,R670,R771, &
                                  RTOPDDOUBLEPRIME,RCMB,RICB, &
                                  RMOHO_FICTITIOUS_IN_MESHER, &
                                  R80_FICTITIOUS_IN_MESHER, &
                                  RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                                  HONOR_1D_SPHERICAL_MOHO,CASE_3D,CRUSTAL)


  end subroutine get_model_parameters


!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_model_parameters_flags(MODEL,REFERENCE_1D_MODEL,THREE_D_MODEL, &
                        ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
                        CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO, &
                        ISOTROPIC_3D_MANTLE,ONE_CRUST,TRANSVERSE_ISOTROPY, &
                        OCEANS,TOPOGRAPHY,CEM_REQUEST,CEM_ACCEPT)

  use constants

  implicit none

  character(len=MAX_STRING_LEN) MODEL

  integer REFERENCE_1D_MODEL,THREE_D_MODEL

  logical ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,ATTENUATION_3D, &
         CASE_3D,CRUSTAL,HETEROGEN_3D_MANTLE,HONOR_1D_SPHERICAL_MOHO,&
         ISOTROPIC_3D_MANTLE,ONE_CRUST,TRANSVERSE_ISOTROPY,&
         CEM_REQUEST,CEM_ACCEPT
  logical OCEANS,TOPOGRAPHY

  ! local parameters
  character(len=4) ending
  character(len=8) ending_1Dcrust

  character(len=MAX_STRING_LEN) MODEL_ROOT
  logical :: impose_1Dcrust

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

  ! extract ending of model name
  ending = ' '
  if (len_trim(MODEL) > 4 ) ending = MODEL(len_trim(MODEL)-3:len_trim(MODEL))

  ! determines if the anisotropic inner core option should be turned on
  if (ending == '_AIC') then
    ANISOTROPIC_INNER_CORE = .true.
    ! in case it has an ending for the inner core, remove it from the name
    MODEL_ROOT = MODEL(1: len_trim(MODEL)-4)
  else
    ANISOTROPIC_INNER_CORE = .false.
    ! sets root name of model to original one
    MODEL_ROOT = MODEL
  endif

  ! checks with '_1Dcrust' option
  impose_1Dcrust = .false.
  ending_1Dcrust = ' '
  if (len_trim(MODEL_ROOT) > 8 ) &
    ending_1Dcrust = MODEL_ROOT(len_trim(MODEL_ROOT)-7:len_trim(MODEL_ROOT))
  if (ending_1Dcrust == '_1Dcrust') then
    impose_1Dcrust = .true.
    ! in case it has an ending for the inner core, remove it from the name
    MODEL_ROOT = MODEL_ROOT(1: len_trim(MODEL)-8)
  endif


!---
!
! ADD YOUR MODEL HERE
!
!---

  ! default values

  ! uses PREM as the 1D reference model by default
  REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM

  ! uses no anisotropic 3D model by default
  ANISOTROPIC_3D_MANTLE = .false.

  ! uses 1D attenuation model by default
  ATTENUATION_3D = .false.

  ! no crustal mesh stretching and 3D crust models by default
  CASE_3D = .false.
  CRUSTAL = .false.
  ONE_CRUST = .false.

  ! uses no 3D heterogeneity mantle by default
  HETEROGEN_3D_MANTLE = .false.
  ISOTROPIC_3D_MANTLE = .false.
  HONOR_1D_SPHERICAL_MOHO = .false.

  ! no CEM by default
  CEM_REQUEST = .false.
  CEM_ACCEPT  = .false.

  ! no 3D model by default
  THREE_D_MODEL = 0
  TRANSVERSE_ISOTROPY = .false.

  ! model specifics

  ! 1-D models
  if (MODEL_ROOT == '1D_isotropic_prem') then
    HONOR_1D_SPHERICAL_MOHO = .true.

  else if (MODEL_ROOT == '1D_transversely_isotropic_prem') then
    HONOR_1D_SPHERICAL_MOHO = .true.
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == '1D_iasp91' .or. MODEL_ROOT == '1D_1066a' .or. &
          MODEL_ROOT == '1D_ak135f_no_mud' .or. MODEL_ROOT == '1D_jp3d' .or. &
          MODEL_ROOT == '1D_sea99') then
    HONOR_1D_SPHERICAL_MOHO = .true.
    if (MODEL_ROOT == '1D_iasp91') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if (MODEL_ROOT == '1D_1066a') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if (MODEL_ROOT == '1D_ak135f_no_mud') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135F_NO_MUD
    else if (MODEL_ROOT == '1D_jp3d') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_JP1D
    else if (MODEL_ROOT == '1D_sea99') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    else
      stop 'reference 1D Earth model unknown'
    endif

  else if (MODEL_ROOT == '1D_ref') then
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == '1D_ref_iso') then
    HONOR_1D_SPHERICAL_MOHO = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF

  else if (MODEL_ROOT == '1D_isotropic_prem_onecrust') then
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if (MODEL_ROOT == '1D_transversely_isotropic_prem_onecrust') then
    TRANSVERSE_ISOTROPY = .true.
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.

  else if (MODEL_ROOT == '1D_iasp91_onecrust' .or. MODEL_ROOT == '1D_1066a_onecrust' &
        .or. MODEL_ROOT == '1D_ak135f_no_mud_onecrust') then
    HONOR_1D_SPHERICAL_MOHO = .true.
    ONE_CRUST = .true.
    if (MODEL_ROOT == '1D_iasp91_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_IASP91
    else if (MODEL_ROOT == '1D_1066a_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_1066A
    else if (MODEL_ROOT == '1D_ak135f_no_mud_onecrust') then
      REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135F_NO_MUD
    else
      stop 'reference 1D Earth model unknown'
    endif

  ! 3-D models
  else if (MODEL_ROOT == 'transversely_isotropic_prem_plus_3D_crust_2.0') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ONE_CRUST = .true.
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 's20rts') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_S20RTS
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 's40rts') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_S40RTS
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 'sea99_jp3d1994') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    THREE_D_MODEL = THREE_D_MODEL_SEA99_JP3D

  else if (MODEL_ROOT == 'sea99') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_SEA1D
    THREE_D_MODEL = THREE_D_MODEL_SEA99

  else if (MODEL_ROOT == 'jp3d1994') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_JP1D
    THREE_D_MODEL = THREE_D_MODEL_JP3D

  else if (MODEL_ROOT == 's362ani') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 's362iso') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  else if (MODEL_ROOT == 's362wmani') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362WMANI
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 's362ani_prem') then
    CASE_3D = .true.
    CRUSTAL = .true.
    TRANSVERSE_ISOTROPY = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_S362ANI_PREM

  else if (MODEL_ROOT == 's362ani_3DQ') then
    ATTENUATION_3D = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.

 else if (MODEL_ROOT == 's362iso_3DQ') then
    ATTENUATION_3D = .true.
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI

  else if (MODEL_ROOT == 's29ea') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S29EA
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == '3D_attenuation') then
    ATTENUATION_3D = .true.
    CASE_3D = .true.
    ONE_CRUST = .true.

  else if (MODEL_ROOT == '3D_anisotropic') then
    ANISOTROPIC_3D_MANTLE = .true.
    CASE_3D = .true. ! crustal moho stretching
    ONE_CRUST = .true. ! 1 element layer in top crust region
    TRANSVERSE_ISOTROPY = .true. ! to use transverse isotropic PREM 1D ref model
    ! CRUSTAL = .true. ! with 3D crust: depends on 3D mantle reference model
    ! THREE_D_MODEL = 0 ! for default crustal model
    ! REFERENCE_1D_MODEL = REFERENCE_MODEL_AK135F_NO_MUD
    ! TRANSVERSE_ISOTROPY = .false. ! for AK135 ref model

  else if (MODEL_ROOT == 'heterogen') then
    CASE_3D = .true.
    CRUSTAL = .true.
    HETEROGEN_3D_MANTLE = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_1DREF
    THREE_D_MODEL = THREE_D_MODEL_S362ANI
    TRANSVERSE_ISOTROPY = .true.

#ifdef CEM
  else if (MODEL_ROOT == 'CEM_REQUEST') then
    CEM_REQUEST         = .true.
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 'CEM_ACCEPT') then
    CEM_ACCEPT          = .true.
    TRANSVERSE_ISOTROPY = .true.

  else if (MODEL_ROOT == 'CEM_GLL') then
    THREE_D_MODEL = THREE_D_MODEL_GLL
    TRANSVERSE_ISOTROPY = .true.
#endif

  else if (MODEL_ROOT == 'PPM') then
    ! superimposed based on isotropic-prem
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    THREE_D_MODEL = THREE_D_MODEL_PPM
    TRANSVERSE_ISOTROPY = .true. ! to use transverse-isotropic prem

  else if (MODEL_ROOT == 'GLL' .or. MODEL_ROOT == 'gll') then
    ! model will be given on local basis, at all GLL points,
    ! as from meshfem3d output from routine save_arrays_solver()
    ! based on model S29EA
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = GLL_REFERENCE_1D_MODEL
    THREE_D_MODEL = THREE_D_MODEL_GLL
    TRANSVERSE_ISOTROPY = .true.
    ! note: after call to this routine we will set
    ! mgll_v%model_gll flag and reset
    ! THREE_D_MODEL = THREE_D_MODEL_S29EA
    ! (not done here because we will use mgll_v%model_gll flag to identify this
    !  model, based upon the S29EA model, but putting mgll_v as parameter to this
    !  routine involves too many changes. )

  else if (MODEL == 'gapp2') then
    CASE_3D = .true.
    CRUSTAL = .true.
    ISOTROPIC_3D_MANTLE = .true.
    ONE_CRUST = .true.
    REFERENCE_1D_MODEL = REFERENCE_MODEL_PREM
    THREE_D_MODEL = THREE_D_MODEL_GAPP2
    TRANSVERSE_ISOTROPY = .true.

  else
    print*
    print*,'Error model: ',trim(MODEL)
    stop 'model not implemented yet, edit get_model_parameters.f90, or ensure you have run ./configure correctly, and recompile'
  endif

  ! suppress the crustal layers
  if (SUPPRESS_CRUSTAL_MESH) then
    CRUSTAL = .false.
    OCEANS = .false.
    ONE_CRUST = .false.
    TOPOGRAPHY = .false.
  endif

  ! additional option for 3D mantle models:
  ! this takes crust from reference 1D model rather than a 3D crust;
  if (impose_1Dcrust) then
    ! no 3D crust
    CRUSTAL = .false.
    ! no crustal moho stretching
    CASE_3D = .false.
    ! mesh honors the 1D moho depth
    HONOR_1D_SPHERICAL_MOHO = .true.
    ! 2 element layers in top crust region rather than just one
    ONE_CRUST = .false.
  endif

  ! checks flag consistency for crust
  if (HONOR_1D_SPHERICAL_MOHO .and. CRUSTAL ) &
    stop 'honor 1D spherical moho excludes having 3D crustal structure'

  ! checks that IASP91, AK135, 1066A, JP1D or SEA1D is isotropic
  if ((REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91 .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_AK135F_NO_MUD .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_1066A .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D .or. &
      REFERENCE_1D_MODEL == REFERENCE_MODEL_SEA1D) .and. TRANSVERSE_ISOTROPY) &
        stop 'models IASP91, AK135, 1066A, JP1D and SEA1D are currently isotropic'


  end subroutine get_model_parameters_flags

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_model_parameters_radii(REFERENCE_1D_MODEL,ROCEAN,RMIDDLE_CRUST, &
                                  RMOHO,R80,R120,R220,R400,R600,R670,R771, &
                                  RTOPDDOUBLEPRIME,RCMB,RICB, &
                                  RMOHO_FICTITIOUS_IN_MESHER, &
                                  R80_FICTITIOUS_IN_MESHER, &
                                  RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS, &
                                  HONOR_1D_SPHERICAL_MOHO,CASE_3D,CRUSTAL)


  use constants

  implicit none

! parameters read from parameter file
  integer REFERENCE_1D_MODEL

  double precision ROCEAN,RMIDDLE_CRUST, &
          RMOHO,R80,R120,R220,R400,R600,R670,R771,RTOPDDOUBLEPRIME,RCMB,RICB, &
          RMOHO_FICTITIOUS_IN_MESHER,R80_FICTITIOUS_IN_MESHER

  double precision RHO_TOP_OC,RHO_BOTTOM_OC,RHO_OCEANS

  logical HONOR_1D_SPHERICAL_MOHO,CASE_3D,CRUSTAL

! radii in PREM or IASP91
! and normalized density at fluid-solid interface on fluid size for coupling
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
  ROCEAN = 6368000.d0
  RMIDDLE_CRUST = 6356000.d0
  RMOHO = 6346600.d0
  R80  = 6291000.d0
  R120 = -1.d0   ! by default there is no d120 discontinuity, except in IASP91, therefore set to fictitious value
  R220 = 6151000.d0
  R400 = 5971000.d0
  R600 = 5771000.d0
  R670 = 5701000.d0
  R771 = 5600000.d0
  RTOPDDOUBLEPRIME = 3630000.d0
  RCMB = 3480000.d0
  RICB = 1221000.d0

  ! density ocean
  RHO_OCEANS = 1020.0 / RHOAV   ! value common to all models
  ! densities fluid outer core
  RHO_TOP_OC = 9903.4384 / RHOAV
  RHO_BOTTOM_OC = 12166.5885 / RHOAV

  ! differing 1-D model radii
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_IASP91) then
    ! IASP91
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6351000.d0
    RMOHO = 6336000.d0
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
    RMOHO  = 6336000.d0
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
    ! values below corrected by Ying Zhou <yingz@gps.caltech.edu>
    ! 1066A
    RMOHO = 6360000.d0
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
    ROCEAN = 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5961000.d0 ! 410km discontinuity
    R600 = 5771000.d0
    R670 = 5721000.d0 ! 650km discontinuity
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3479958.d0
    RICB = 1221491.d0

    RHO_TOP_OC = 9903.48 / RHOAV
    RHO_BOTTOM_OC = 12166.35 / RHOAV

  else if (REFERENCE_1D_MODEL == REFERENCE_MODEL_JP1D) then
    ! values below corrected by Min Chen <mchen@gps.caltech.edu>
    ! jp1d
    ROCEAN = 6371000.d0
    RMIDDLE_CRUST = 6359000.d0
    RMOHO = 6345000.d0
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
    RMOHO  = 6346000.d0
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

    RHO_TOP_OC = 9903.4384 / RHOAV
    RHO_BOTTOM_OC = 12166.5885 / RHOAV

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

  end subroutine get_model_parameters_radii

