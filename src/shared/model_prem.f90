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

!--------------------------------------------------------------------------------------------------
!
! PREM [Dziewonski and Anderson, 1981].
!
! A. M. Dziewonski and D. L. Anderson.
! Preliminary reference Earth model.
! Phys. Earth Planet. Inter., 25:297-356, 1981.
!
! Isotropic (iso) and transversely isotropic (aniso) version of the
! spherically symmetric Preliminary Reference Earth Model
!
!
! PREM2 [Song and Helmberger, 1995]
!
! X. Song and D. V. Helmberger.
! A P wave velocity model of Earth's core.
! J. Geophys. Res., 100, B7, pp. 9817-9830, 1995
!
! major feature of PREM2 are:
! (1) large velocity jump at the inner core boundary
! (2) near zero gradient for bottom 150km of the fluid core
! (3) nearly constant t*=0.35s for inner core attenuation
! PREM2 has an average velocity lower than PREM in the lowermost mantle.
!
!--------------------------------------------------------------------------------------------------

  module model_prem_par

  implicit none

  ! Preliminary reference Earth model (PREM)
  double precision, parameter :: PREM_RSURFACE = 6371000.d0
  ! uncomment line below for PREM without oceans
  !double precision, parameter :: PREM_RSURFACE = 6368000.d0

  ! PREM radii (in m)
  double precision, parameter :: PREM_ROCEAN = 6368000.d0           ! at   3 km depth
  double precision, parameter :: PREM_RMIDDLE_CRUST = 6356000.d0    !     15 km depth
  double precision, parameter :: PREM_RMOHO = PREM_RSURFACE - 24400.d0  ! 6346600.d0  at 24.4km depth
  double precision, parameter :: PREM_R80  = 6291000.d0             !     80 km
  double precision, parameter :: PREM_R220 = 6151000.d0
  double precision, parameter :: PREM_R400 = 5971000.d0
  double precision, parameter :: PREM_R600 = 5771000.d0
  double precision, parameter :: PREM_R670 = 5701000.d0             ! at 670 km depth
  double precision, parameter :: PREM_R771 = 5600000.d0
  double precision, parameter :: PREM_RTOPDDOUBLEPRIME = 3630000.d0
  double precision, parameter :: PREM_RCMB = 3480000.d0             !   2891 km depth
  ! note: SPECFEM versions up to 8.0 (Aug, 2021) used an inner core radius of 1221 km;
  !       based on Table 1 in Dziewonski&Anderson's PREM paper, the inner core radius is at 1221.5 km
  !double precision, parameter :: PREM_RICB = 1221000.d0            ! old versions
  double precision, parameter :: PREM_RICB = 1221500.d0

  ! PREM2 additional radii for modifications
  double precision, parameter :: PREM2_RDDOUBLEPRIME_UPPER = 3840000.d0 ! upper D'' region at 3840km radius
  double precision, parameter :: PREM2_ROC_LOWER = 1621500.d0     ! lower outer core at 1621.5km radius
  double precision, parameter :: PREM2_RIC_UPPER = 1010000.d0     ! upper inner core at 1010 km radius


  ! densities
  double precision, parameter :: PREM_RHO_OCEANS = 1020.0         ! for ocean load
  double precision, parameter :: PREM_RHO_TOP_OC = 9903.4384      ! for coupling with CMB
  double precision, parameter :: PREM_RHO_BOTTOM_OC = 12166.5885  ! for coupling with ICB

  end module model_prem_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_prem_iso(x,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,check_doubling_flag)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,ONE_CRUST,REFERENCE_1D_MODEL

  use model_prem_par

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  double precision,intent(in) :: x
  double precision,intent(inout) :: rho,drhodr,vp,vs,Qkappa,Qmu

  logical,intent(in) :: CRUSTAL,check_doubling_flag
  integer,intent(in) :: idoubling

  ! local parameters
  double precision :: r,scaleval

! compute real physical radius in meters
  r = x * R_PLANET

! check flags to make sure we correctly honor the discontinuities
! we use strict inequalities since r has been slightly changed in mesher

! note: using stop statements, not exit_mpi() calls to avoid the need for MPI libraries when linking xcreate_header_file

 if (check_doubling_flag) then
    !
    !--- inner core
    !
    if (r >= 0.d0 .and. r < PREM_RICB) then
      if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
         idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
         idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
         idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
         idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
           stop 'wrong doubling flag for inner core point in model_prem_iso()'
    !
    !--- outer core
    !
    else if (r > PREM_RICB .and. r < PREM_RCMB) then
      if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
        stop 'wrong doubling flag for outer core point in model_prem_iso()'
    !
    !--- D" at the base of the mantle
    !
    else if (r > PREM_RCMB .and. r < PREM_RTOPDDOUBLEPRIME) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for D" point in model_prem_iso()'
    !
    !--- mantle: from top of D" to d670
    !
    else if (r > PREM_RTOPDDOUBLEPRIME .and. r < PREM_R670) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for top D" to d670 point in model_prem_iso()'
    !
    !--- mantle: from d670 to d220
    !
    else if (r > PREM_R670 .and. r < PREM_R220) then
      if (idoubling /= IFLAG_670_220) &
        stop 'wrong doubling flag for d670 to d220 point in model_prem_iso()'
    !
    !--- mantle and crust: from d220 to MOHO and then to surface
    !
    else if (r > PREM_R220) then
      if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
        stop 'wrong doubling flag for d220 to Moho to surface point in model_prem_iso()'
    endif

  endif

!
!--- inner core
!
  if (r >= 0.d0 .and. r <= PREM_RICB) then
    drhodr = -2.0d0*8.8381d0*x
    rho = 13.0885d0 - 8.8381d0*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2
      if (r <= PREM2_RIC_UPPER) then
        ! lower inner core, same value as PREM
        vp = 11.2622d0 - 6.3640d0*x*x
      else
        ! upper inner core modification
        vp = 11.3041d0 - 1.2730d0*x
      endif
    else
      ! PREM
      vp = 11.2622d0 - 6.3640d0*x*x
    endif
    vs = 3.6678d0 - 4.4475d0*x*x
    Qmu = 84.6d0
    Qkappa = 1327.7d0
!
!--- outer core
!
  else if (r > PREM_RICB .and. r <= PREM_RCMB) then
    drhodr = -1.2638d0 - 2.0d0*3.6426d0*x - 3.0d0*5.5281d0*x*x
    rho = 12.5815d0 - 1.2638d0*x - 3.6426d0*x*x - 5.5281d0*x*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2
      if (r <= PREM2_ROC_LOWER) then
        ! lower outer core modification
        vp = 4.0354d0 + 82.0080d0*x - 347.7690d0*x*x + 468.7860d0*x*x*x
      else
        ! upper outer core, same value as PREM
        vp = 11.0487d0 - 4.0362d0*x + 4.8023d0*x*x - 13.5732d0*x*x*x
      endif
    else
      ! PREM
      vp = 11.0487d0 - 4.0362d0*x + 4.8023d0*x*x - 13.5732d0*x*x*x
    endif
    vs = 0.0d0
    Qmu = 0.0d0
    Qkappa = 57827.0d0
!
!--- D" at the base of the mantle
!
  else if (r > PREM_RCMB .and. r <= PREM_RTOPDDOUBLEPRIME) then
    drhodr = -6.4761d0 + 2.0d0*5.5283d0*x - 3.0d0*3.0807d0*x*x
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2 modification between RCMB and D'' upper region
      vp = 14.2743d0 - 1.3998d0*x
    else
      ! PREM
      vp = 15.3891d0 - 5.3181d0*x + 5.5242d0*x*x - 2.5514d0*x*x*x
    endif
    vs = 6.9254d0 + 1.4672d0*x - 2.0834d0*x*x + 0.9783d0*x*x*x
    Qmu = 312.0d0
    Qkappa = 57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if (r > PREM_RTOPDDOUBLEPRIME .and. r <= PREM_R771) then
    drhodr = -6.4761d0 + 2.0d0*5.5283d0*x - 3.0d0*3.0807d0*x*x
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2
      if (r <= PREM2_RDDOUBLEPRIME_UPPER) then
        ! PREM2 modification between RCMB and D'' upper region
        vp = 14.2743d0 - 1.3998d0*x
      else
        ! above D'' extended region, same value as PREM
        vp = 24.9520d0 - 40.4673d0*x + 51.4832d0*x*x - 26.6419d0*x*x*x
      endif
    else
      ! PREM
      vp = 24.9520d0 - 40.4673d0*x + 51.4832d0*x*x - 26.6419d0*x*x*x
    endif
    vs = 11.1671d0 - 13.7818d0*x + 17.4575d0*x*x - 9.2777d0*x*x*x
    Qmu = 312.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R771 .and. r <= PREM_R670) then
    drhodr = -6.4761d0 + 2.0d0*5.5283d0*x - 3.0d0*3.0807d0*x*x
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
    vp = 29.2766d0 - 23.6027d0*x + 5.5242d0*x*x - 2.5514d0*x*x*x
    vs = 22.3459d0 - 17.2473d0*x - 2.0834d0*x*x + 0.9783d0*x*x*x
    Qmu = 312.0d0
    Qkappa = 57827.0d0
!
!--- mantle: above d670
!
  else if (r > PREM_R670 .and. r <= PREM_R600) then
    drhodr = -1.4836d0
    rho = 5.3197d0 - 1.4836d0*x
    vp = 19.0957d0 - 9.8672d0*x
    vs = 9.9839d0 - 4.9324d0*x
    Qmu = 143.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R600 .and. r <= PREM_R400) then
    drhodr = -8.0298d0
    rho = 11.2494d0 - 8.0298d0*x
    vp = 39.7027d0 - 32.6166d0*x
    vs = 22.3512d0 - 18.5856d0*x
    Qmu = 143.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R400 .and. r <= PREM_R220) then
    drhodr = -3.8045d0
    rho = 7.1089d0 - 3.8045d0*x
    vp = 20.3926d0 - 12.2569d0*x
    vs = 8.9496d0 - 4.4597d0*x
    Qmu = 143.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R220 .and. r <= PREM_R80) then
    drhodr = 0.6924d0
    rho = 2.6910d0 + 0.6924d0*x
    vp = 4.1875d0 + 3.9382d0*x
    vs = 2.1519d0 + 2.3481d0*x
    Qmu = 80.0d0
    Qkappa = 57827.0d0
  else
    if (CRUSTAL .and. .not. SUPPRESS_CRUSTAL_MESH) then
! fill with PREM mantle and later add CRUST2.0
      if (r > PREM_R80) then
        ! density/velocity from mantle just below moho
        drhodr = 0.6924d0
        rho = 2.6910d0 + 0.6924d0*x
        vp = 4.1875d0 + 3.9382d0*x
        vs = 2.1519d0 + 2.3481d0*x
        ! shear attenuation for R80 to surface
        Qmu = 600.0d0
        Qkappa = 57827.0d0
      endif
    else
! use PREM crust
      if (r > PREM_R80 .and. r <= PREM_RMOHO) then
        drhodr = 0.6924d0
        rho = 2.6910d0 + 0.6924d0*x
        vp = 4.1875d0 + 3.9382d0*x
        vs = 2.1519d0 + 2.3481d0*x
        Qmu = 600.0d0
        Qkappa = 57827.0d0

      else if (SUPPRESS_CRUSTAL_MESH) then
!! DK DK extend the Moho up to the surface instead of the crust
        drhodr = 0.6924d0
        rho = 2.6910d0 + 0.6924d0*(PREM_RMOHO / R_PLANET)
        vp = 4.1875d0 + 3.9382d0*(PREM_RMOHO / R_PLANET)
        vs = 2.1519d0 + 2.3481d0*(PREM_RMOHO / R_PLANET)
        Qmu = 600.0d0
        Qkappa = 57827.0d0

      else if (r > PREM_RMOHO .and. r <= PREM_RMIDDLE_CRUST) then
        drhodr = 0.0d0
        rho = 2.9d0
        vp = 6.8d0
        vs = 3.9d0
        Qmu = 600.0d0
        Qkappa = 57827.0d0

! same properties everywhere in PREM crust if we decide to define only one layer in the crust
        if (ONE_CRUST) then
          drhodr = 0.0d0
          rho = 2.6d0
          vp = 5.8d0
          vs = 3.2d0
          Qmu = 600.0d0
          Qkappa = 57827.0d0
        endif

      else if (r > PREM_RMIDDLE_CRUST .and. r <= PREM_ROCEAN) then
        drhodr = 0.0d0
        rho = 2.6d0
        vp = 5.8d0
        vs = 3.2d0
        Qmu = 600.0d0
        Qkappa = 57827.0d0
! for density profile for gravity, we do not check that r <= R_PLANET
      else if (r > PREM_ROCEAN) then
        drhodr = 0.0d0
        rho = 2.6d0
        vp = 5.8d0
        vs = 3.2d0
        Qmu = 600.0d0
        Qkappa = 57827.0d0

      endif
    endif
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  drhodr = drhodr * 1000.0d0/RHOAV
  rho = rho * 1000.0d0/RHOAV
  vp = vp * 1000.0d0/(R_PLANET*scaleval)
  vs = vs * 1000.0d0/(R_PLANET*scaleval)

  end subroutine model_prem_iso

!
!=====================================================================
!

  subroutine model_prem_aniso(x,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,idoubling,CRUSTAL)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,ONE_CRUST,REFERENCE_1D_MODEL

  use model_prem_par

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  logical,intent(in) :: CRUSTAL

  integer,intent(in) :: idoubling

  double precision,intent(in) :: x
  double precision,intent(inout) :: rho,Qkappa,Qmu,vpv,vph,vsv,vsh,eta_aniso

  ! local parameters
  double precision :: r
  double precision :: scaleval

! compute real physical radius in meters
  r = x * R_PLANET

! check flags to make sure we correctly honor the discontinuities
! we use strict inequalities since r has been slightly changed in mesher

!
!--- inner core
!
  if (r >= 0.d0 .and. r < PREM_RICB) then
    if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
       idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
         stop 'wrong doubling flag for inner core point in model_prem_aniso()'
!
!--- outer core
!
  else if (r > PREM_RICB .and. r < PREM_RCMB) then
    if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
      stop 'wrong doubling flag for outer core point in model_prem_aniso()'
!
!--- D" at the base of the mantle
!
  else if (r > PREM_RCMB .and. r < PREM_RTOPDDOUBLEPRIME) then
    if (idoubling /= IFLAG_MANTLE_NORMAL) then
      print *,'Error dprime point:',r, PREM_RCMB,PREM_RTOPDDOUBLEPRIME,idoubling,IFLAG_MANTLE_NORMAL
      stop 'wrong doubling flag for D" point in model_prem_aniso()'
    endif
!
!--- mantle: from top of D" to d670
!
  else if (r > PREM_RTOPDDOUBLEPRIME .and. r < PREM_R670) then
    if (idoubling /= IFLAG_MANTLE_NORMAL) &
      stop 'wrong doubling flag for top D" to d670 point in model_prem_aniso()'

!
!--- mantle: from d670 to d220
!
  else if (r > PREM_R670 .and. r < PREM_R220) then
    if (idoubling /= IFLAG_670_220) &
      stop 'wrong doubling flag for d670 to d220 point in model_prem_aniso()'

!
!--- mantle and crust: from d220 to MOHO and then to surface
!
  else if (r > PREM_R220) then
    if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
      stop 'wrong doubling flag for d220 to Moho to surface point in model_prem_aniso()'

  endif

! no anisotropy by default
  eta_aniso = 1.d0

!
!--- inner core
!
  if (r >= 0.d0 .and. r <= PREM_RICB) then
    rho = 13.0885d0 - 8.8381d0*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2
      if (r <= PREM2_RIC_UPPER) then
        ! lower inner core, same value as PREM
        vpv = 11.2622d0 - 6.3640d0*x*x
      else
        ! upper inner core modification
        vpv = 11.3041d0 - 1.2730d0*x
      endif
    else
      ! PREM
      vpv = 11.2622d0 - 6.3640d0*x*x
    endif
    vsv = 3.6678d0 - 4.4475d0*x*x
    vph = vpv
    vsh = vsv
    Qmu = 84.6d0
    Qkappa = 1327.7d0
!
!--- outer core
!
  else if (r > PREM_RICB .and. r <= PREM_RCMB) then
    rho = 12.5815d0 - 1.2638d0*x - 3.6426d0*x*x - 5.5281d0*x*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2
      if (r <= PREM2_ROC_LOWER) then
        ! lower outer core modification
        vpv = 4.0354d0 + 82.0080d0*x - 347.7690d0*x*x + 468.7860d0*x*x*x
      else
        ! upper outer core, same value as PREM
        vpv = 11.0487d0 - 4.0362d0*x + 4.8023d0*x*x - 13.5732d0*x*x*x
      endif
    else
      ! PREM
      vpv = 11.0487d0 - 4.0362d0*x + 4.8023d0*x*x - 13.5732d0*x*x*x
    endif
    vsv = 0.0d0
    vph = vpv
    vsh = vsv
    Qmu = 0.0d0
    Qkappa = 57827.0d0
!
!--- D" at the base of the mantle
!
  else if (r > PREM_RCMB .and. r <= PREM_RTOPDDOUBLEPRIME) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2 modification between RCMB and D'' upper region
      vpv = 14.2743d0 - 1.3998d0*x
    else
      ! PREM
      vpv = 15.3891d0 - 5.3181d0*x + 5.5242d0*x*x - 2.5514d0*x*x*x
    endif
    vsv = 6.9254d0 + 1.4672d0*x - 2.0834d0*x*x + 0.9783d0*x*x*x
    vph = vpv
    vsh = vsv
    Qmu = 312.0d0
    Qkappa = 57827.0d0
!
!--- mantle: from top of D" to d670
!
  else if (r > PREM_RTOPDDOUBLEPRIME .and. r <= PREM_R771) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
    if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
      ! PREM2
      if (r <= PREM2_RDDOUBLEPRIME_UPPER) then
        ! PREM2 modification between RCMB and D'' upper region
        vpv = 14.2743d0 - 1.3998d0*x
      else
        ! above D'' extended region, same value as PREM
        vpv = 24.9520d0 - 40.4673d0*x + 51.4832d0*x*x - 26.6419d0*x*x*x
      endif
    else
      ! PREM
      vpv = 24.9520d0 - 40.4673d0*x + 51.4832d0*x*x - 26.6419d0*x*x*x
    endif
    vsv = 11.1671d0 - 13.7818d0*x + 17.4575d0*x*x - 9.2777d0*x*x*x
    vph = vpv
    vsh = vsv
    Qmu = 312.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R771 .and. r <= PREM_R670) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
    vpv = 29.2766d0 - 23.6027d0*x + 5.5242d0*x*x - 2.5514d0*x*x*x
    vsv = 22.3459d0 - 17.2473d0*x - 2.0834d0*x*x + 0.9783d0*x*x*x
    vph = vpv
    vsh = vsv
    Qmu = 312.0d0
    Qkappa = 57827.0d0
!
!--- mantle: above d670
!
  else if (r > PREM_R670 .and. r <= PREM_R600) then
    rho = 5.3197d0 - 1.4836d0*x
    vpv = 19.0957d0 - 9.8672d0*x
    vsv = 9.9839d0 - 4.9324d0*x
    vph = vpv
    vsh = vsv
    Qmu = 143.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R600 .and. r <= PREM_R400) then
    rho = 11.2494d0 - 8.0298d0*x
    vpv = 39.7027d0 - 32.6166d0*x
    vsv = 22.3512d0 - 18.5856d0*x
    vph = vpv
    vsh = vsv
    Qmu = 143.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R400 .and. r <= PREM_R220) then
    rho = 7.1089d0 - 3.8045d0*x
    vpv = 20.3926d0 - 12.2569d0*x
    vsv = 8.9496d0 - 4.4597d0*x
    vph = vpv
    vsh = vsv
    Qmu = 143.0d0
    Qkappa = 57827.0d0
  else if (r > PREM_R220 .and. r <= PREM_R80) then

! anisotropy in PREM only above 220 km

    rho = 2.6910d0 + 0.6924d0*x
    vpv = 0.8317d0 + 7.2180d0*x
    vph = 3.5908d0 + 4.6172d0*x
    vsv = 5.8582d0 - 1.4678d0*x
    vsh = -1.0839d0 + 5.7176d0*x
    eta_aniso = 3.3687d0 - 2.4778d0*x
    Qmu = 80.0d0
    Qkappa = 57827.0d0

  else
    if (CRUSTAL) then
! fill with PREM mantle and later add CRUST2.0
      if (r > PREM_R80) then
        rho = 2.6910d0 + 0.6924d0*x
        vpv = 0.8317d0 + 7.2180d0*x
        vph = 3.5908d0 + 4.6172d0*x
        vsv = 5.8582d0 - 1.4678d0*x
        vsh = -1.0839d0 + 5.7176d0*x
        eta_aniso = 3.3687d0 - 2.4778d0*x
        Qmu = 600.0d0
        Qkappa = 57827.0d0
      endif
    else
! use PREM crust
      if (r > PREM_R80 .and. r <= PREM_RMOHO) then

! anisotropy in PREM only above 220 km

        rho = 2.6910d0 + 0.6924d0*x
        vpv = 0.8317d0 + 7.2180d0*x
        vph = 3.5908d0 + 4.6172d0*x
        vsv = 5.8582d0 - 1.4678d0*x
        vsh = -1.0839d0 + 5.7176d0*x
        eta_aniso = 3.3687d0 - 2.4778d0*x
        Qmu = 600.0d0
        Qkappa = 57827.0d0

! no anisotropy in the crust in PREM

      else if (r > PREM_RMOHO .and. r <= PREM_RMIDDLE_CRUST) then
        rho = 2.9d0
        vpv = 6.8d0
        vsv = 3.9d0
        vph = vpv
        vsh = vsv
        Qmu = 600.0d0
        Qkappa = 57827.0d0

! same properties everywhere in PREM crust (only one layer in the crust)
        if (ONE_CRUST) then
          rho = 2.6d0
          vpv = 5.8d0
          vsv = 3.2d0
          vph = vpv
          vsh = vsv
          Qmu = 600.0d0
          Qkappa = 57827.0d0
        endif

      else if (r > PREM_RMIDDLE_CRUST .and. r <= PREM_ROCEAN) then
        rho = 2.6d0
        vpv = 5.8d0
        vsv = 3.2d0
        vph = vpv
        vsh = vsv
        Qmu = 600.0d0
        Qkappa = 57827.0d0
      else if (r > PREM_ROCEAN) then
        rho = 2.6d0
        vpv = 5.8d0
        vsv = 3.2d0
        vph = vpv
        vsh = vsv
        Qmu = 600.0d0
        Qkappa = 57827.0d0
      endif
    endif
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
! do not scale anisotropy parameter eta_aniso, which is dimensionless
  scaleval = dsqrt(PI*GRAV*RHOAV)
  rho = rho * 1000.0d0/RHOAV
  vpv = vpv * 1000.0d0/(R_PLANET*scaleval)
  vsv = vsv * 1000.0d0/(R_PLANET*scaleval)
  vph = vph * 1000.0d0/(R_PLANET*scaleval)
  vsh = vsh * 1000.0d0/(R_PLANET*scaleval)

  end subroutine model_prem_aniso

!
!=====================================================================
!

  subroutine model_prem_aniso_extended_isotropic(x,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu, &
                                                 idoubling,CRUSTAL)

! note: for 3D crustal models, we extend the mantle reference up to the surface and then superimpose the crustal values later.
!       however, PREM mantle is anisotropic (eta < 1 and vsh > vsv) and the extension continues with strong TISO, thus
!       extending it with anisotropic values up to the surface will create very fast, strong mantle elements
!       below thin oceanic crust. in particular, oceanic paths will see Love waves becoming too fast for mantle models
!       like S20RTS using PREM as reference model.
!
!       here, we will extend the mantle with isotropic values (averaged), and no radial anisotropy above moho depths.

  use constants
  use shared_parameters, only: R_PLANET,RHOAV

  use model_prem_par

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  double precision,intent(in) :: x
  double precision,intent(out) :: rho,Qkappa,Qmu,vpv,vph,vsv,vsh,eta_aniso

  logical,intent(in) :: CRUSTAL
  integer,intent(in) :: idoubling

  ! local parameters
  double precision :: r,vp,vs
  double precision :: scaleval

  ! gets default values
  call model_prem_aniso(x,rho,vpv,vph,vsv,vsh,eta_aniso,Qkappa,Qmu,idoubling,CRUSTAL)

  if (CRUSTAL) then
    ! adds crustal model like CRUST2.0 later on top

    ! compute real physical radius in meters
    r = x * R_PLANET

    ! fills with PREM mantle only till moho depth
    if (r > PREM_R80 .and. r <= PREM_RMOHO) then
      rho = 2.6910d0 + 0.6924d0*x

      vpv = 0.8317d0 + 7.2180d0*x
      vph = 3.5908d0 + 4.6172d0*x

      vsv = 5.8582d0 - 1.4678d0*x
      vsh = -1.0839d0 + 5.7176d0*x

      eta_aniso = 3.3687d0 - 2.4778d0*x

      ! non-dimensionalize
      ! time scaling (s^{-1}) is done with scaleval
      ! do not scale anisotropy parameter eta_aniso, which is dimensionless
      scaleval = dsqrt(PI*GRAV*RHOAV)
      rho = rho * 1000.0d0/RHOAV
      vpv = vpv * 1000.0d0/(R_PLANET*scaleval)
      vsv = vsv * 1000.0d0/(R_PLANET*scaleval)
      vph = vph * 1000.0d0/(R_PLANET*scaleval)
      vsh = vsh * 1000.0d0/(R_PLANET*scaleval)

    ! extension to surface
    else if (r > PREM_RMOHO) then
      ! extends the mantle up to the surface
      ! uses isotropic values (averaged), no anisotropy above moho depths

      ! continues mantle trend
      rho = 2.6910d0 + 0.6924d0*x

      ! P-velocity
      vpv = 0.8317d0 + 7.2180d0*x
      vph = 3.5908d0 + 4.6172d0*x

      ! Voigt average:  v_avg**2 = ( 2 vpv**2 + vph**2 ) / 3
      vp = sqrt( (2.d0 * vpv**2 + vph**2) / 3.d0 )
      vpv = vp
      vph = vp

      ! S-velocity
      vsv = 5.8582d0 - 1.4678d0*x
      vsh = -1.0839d0 + 5.7176d0*x

      ! Voigt average:  v_avg**2 = ( 2 vsv**2 + vsh**2 ) / 3
      vs = sqrt( (2.d0 * vsv**2 + vsh**2) / 3.d0 )
      vsv = vs
      vsh = vs

      ! no radial anisotropy by default
      eta_aniso = 1.d0

      ! non-dimensionalize
      ! time scaling (s^{-1}) is done with scaleval
      ! do not scale anisotropy parameter eta_aniso, which is dimensionless
      scaleval = dsqrt(PI*GRAV*RHOAV)
      rho = rho * 1000.0d0/RHOAV
      vpv = vpv * 1000.0d0/(R_PLANET*scaleval)
      vsv = vsv * 1000.0d0/(R_PLANET*scaleval)
      vph = vph * 1000.0d0/(R_PLANET*scaleval)
      vsh = vsh * 1000.0d0/(R_PLANET*scaleval)
    endif
  endif ! CRUSTAL

  end subroutine model_prem_aniso_extended_isotropic

!
!=====================================================================
!

  subroutine prem_display_outer_core(x,rho,vp,vs,Qkappa,Qmu,idoubling)

! routine used for AVS or DX display of stability condition
! and number of points per wavelength only in the fluid outer core

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,REFERENCE_1D_MODEL

  use model_prem_par, only: PREM2_ROC_LOWER

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  integer,intent(in) :: idoubling
  double precision,intent(in) :: x
  double precision,intent(inout) :: rho,vp,vs,Qkappa,Qmu

  ! local parameters
  double precision :: scaleval,r

  if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
    stop 'wrong doubling flag for outer core point in prem_display_outer_core()'

!
!--- outer core
!
  rho    = 12.5815d0 - 1.2638d0 * x - 3.6426d0 * x*x - 5.5281d0 * x*x*x
  if (REFERENCE_1D_MODEL == REFERENCE_MODEL_PREM2) then
    ! PREM2
    r = x * R_PLANET
    if (r <= PREM2_ROC_LOWER) then
      ! lower outer core modification
      vp = 4.0354d0 + 82.0080d0*x - 347.7690d0*x*x + 468.7860d0*x*x*x
    else
      ! upper outer core, same value as PREM
      vp = 11.0487d0 - 4.0362d0*x + 4.8023d0*x*x - 13.5732d0*x*x*x
    endif
  else
    ! PREM
    vp     = 11.0487d0 - 4.0362d0 * x + 4.8023d0 * x*x - 13.5732d0 * x*x*x
  endif
  vs     = 0.0d0
  Qmu    = 0.0d0
  Qkappa = 57827.0d0

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  rho = rho*1000.0d0/RHOAV
  vp = vp*1000.0d0/(R_PLANET*scaleval)
  vs = vs*1000.0d0/(R_PLANET*scaleval)

  end subroutine prem_display_outer_core

!
!=====================================================================
!

  subroutine prem_density(x,rho)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,ONE_CRUST

  use model_prem_par

  implicit none

  double precision,intent(in) :: x
  double precision,intent(out) :: rho

  ! local parameters
  double precision :: r

  ! compute real physical radius in meters
  r = x * R_PLANET

  ! calculates density according to radius
  if (r <= PREM_RICB) then
    rho = 13.0885d0 - 8.8381d0*x*x
  else if (r > PREM_RICB .and. r <= PREM_RCMB) then
    rho = 12.5815d0 - 1.2638d0*x - 3.6426d0*x*x - 5.5281d0*x*x*x
  else if (r > PREM_RCMB .and. r <= PREM_RTOPDDOUBLEPRIME) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
  else if (r > PREM_RTOPDDOUBLEPRIME .and. r <= PREM_R771) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
  else if (r > PREM_R771 .and. r <= PREM_R670) then
    rho = 7.9565d0 - 6.4761d0*x + 5.5283d0*x*x - 3.0807d0*x*x*x
  else if (r > PREM_R670 .and. r <= PREM_R600) then
    rho = 5.3197d0 - 1.4836d0*x
  else if (r > PREM_R600 .and. r <= PREM_R400) then
    rho = 11.2494d0 - 8.0298d0*x
  else if (r > PREM_R400 .and. r <= PREM_R220) then
    rho = 7.1089d0 - 3.8045d0*x
  else if (r > PREM_R220 .and. r <= PREM_R80) then
    rho = 2.6910d0 + 0.6924d0*x
  else
    if (r > PREM_R80 .and. r <= PREM_RMOHO) then
      rho = 2.6910d0 + 0.6924d0*x
    else if (r > PREM_RMOHO .and. r <= PREM_RMIDDLE_CRUST) then
      if (ONE_CRUST) then
        rho = 2.6d0  ! takes upper crust value
      else
        rho = 2.9d0
      endif
    else if (r > PREM_RMIDDLE_CRUST .and. r <= PREM_ROCEAN) then
      rho = 2.6d0
    else if (r > PREM_ROCEAN) then
      rho = 2.6d0 ! extends upper crust
    endif
  endif

  ! non-dimensionalizes
  rho = rho * 1000.0d0 / RHOAV

  end subroutine prem_density

