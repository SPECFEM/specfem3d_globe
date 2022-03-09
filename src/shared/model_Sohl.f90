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
! Isotropic version of the spherically symmetric Reference Mars Model
! model A in the reference
!
! Sohl, F., and T. Spohn (1997), The interior structure of Mars: Implications
! from SNC meteorites, J. Geophys. Res., 102(E1), 1613-1635,
! doi:10.1029/96JE03419.
!
!
!--------------------------------------------------------------------------------------------------


  module model_sohl_par

  implicit none

  ! radius at physical surface
  double precision, parameter :: SOHL_RSURFACE = 3390000.d0
  ! no ocean, same as surface radius (Physical surface)
  double precision, parameter :: SOHL_ROCEAN = SOHL_RSURFACE
  ! some mid-crust layer for meshing
  double precision, parameter :: SOHL_RMIDDLE_CRUST = 3340000.d0  ! depth = 50 km

  ! Sohl & Spohn, 1997, propose two models:
  !   * model A: mostly constrained by geophysical moment of inertia
  !   * model B: mostly constrained by global Fe/Si composition by geochemical requirement
  ! we will use model A as default
  !
  ! Model A
  ! Sohl & Spohn 1997, Table 5, pg. 1623
  ! Crust-mantle boundary
  double precision, parameter :: SOHL_A_RMOHO = 3280000.d0      ! depth = 110 km average crustal thickness
  ! Rheological lithosphere
  double precision, parameter :: SOHL_A_R80  = 3055500.d0       ! depth = 334.5 km, if smaller will cause negative Jacobian err
  ! Thermal lithosphere
  double precision, parameter :: SOHL_A_R220 = 2908000.d0       ! depth = 482 km
  double precision, parameter :: SOHL_A_R400 = 2655000.d0       ! depth = 735 km (added for meshing)
  double precision, parameter :: SOHL_A_R600 = 2455000.d0       ! depth = 935 km (added for meshing)
  ! alpha-olivine-beta-spinel transition
  double precision, parameter :: SOHL_A_R670 = 2360000.d0       ! depth = 1030 km
  ! beta-spinel-gamma-spinel transition
  double precision, parameter :: SOHL_A_R771 = 2033000.d0       ! depth = 1357 km, below which the second doubling implemented
  ! lower thermal boundary layer
  double precision, parameter :: SOHL_A_RTOPDDOUBLEPRIME = 1503000.d0 ! depth = 1887 km, too thin for mesh ?
  ! Core-mantle boundary
  double precision, parameter :: SOHL_A_RCMB = 1468000.d0       ! depth = 1922 km

  ! Model B
  ! Sohl & Spohn 1997, Table 6, pg. 1623
  ! Crust-mantle boundary
  double precision, parameter :: SOHL_B_RMOHO = 3138000.d0      ! depth =  252 km
  ! Rheological lithosphere
  double precision, parameter :: SOHL_B_R80  = 3016000.d0       ! depth =  374 km
  ! Thermal lithosphere
  double precision, parameter :: SOHL_B_R220 = 2783000.d0       ! depth =  607 km
  double precision, parameter :: SOHL_B_R400 = 2655000.d0       ! depth =  735 km (added for meshing)
  double precision, parameter :: SOHL_B_R600 = 2455000.d0       ! depth =  935 km (added for meshing)
  ! alpha-olivine-beta-spinel transition
  double precision, parameter :: SOHL_B_R670 = 2332000.d0       ! depth = 1058 km
  ! beta-spinel-gamma-spinel transition
  double precision, parameter :: SOHL_B_R771 = 1974000.d0       ! depth = 1416 km
  ! lower thermal boundary layer
  double precision, parameter :: SOHL_B_RTOPDDOUBLEPRIME = 1700000.d0 ! depth = 1690 km
  ! Core-mantle boundary
  double precision, parameter :: SOHL_B_RCMB = 1667000.d0       ! depth = 1723 km

  ! note: Sohl & Spohn assume an entirely molten iron alloy core.
  !       since SPECFEM assumes a structure of solid crust/mantle - fluid outer core - solid inner core,
  !       we set a very small inner core radius here.
  double precision, parameter :: SOHL_RICB = 515000.d0          ! depth = 2875 km good for both stability and efficiency

  ! densities
  double precision, parameter :: SOHL_RHO_OCEANS = 1020.0         ! will not be used
  double precision, parameter :: SOHL_RHO_TOP_OC = 6936.40        ! densities fluid outer core (from modSOHL)
  double precision, parameter :: SOHL_RHO_BOTTOM_OC = 7268.20

  end module model_sohl_par

!
!-------------------------------------------------------------------------------------------------
!


  subroutine model_Sohl(x,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL,check_doubling_flag)

  use constants, only: PI,GRAV, &
    IFLAG_CRUST,IFLAG_220_80,IFLAG_670_220,IFLAG_80_MOHO,IFLAG_MANTLE_NORMAL, &
    IFLAG_TOP_CENTRAL_CUBE,IFLAG_INNER_CORE_NORMAL,IFLAG_OUTER_CORE_NORMAL, &
    IFLAG_BOTTOM_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE,IFLAG_MIDDLE_CENTRAL_CUBE, &
    SUPPRESS_CRUSTAL_MESH, &
    REFERENCE_MODEL_SOHL,REFERENCE_MODEL_SOHL_B,REFERENCE_MODEL_CASE65TAY

  use shared_parameters, only: R_PLANET,RHOAV,REFERENCE_1D_MODEL,ONE_CRUST

  use model_sohl_par

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  logical,intent(in) :: CRUSTAL,check_doubling_flag

  integer,intent(in) :: idoubling

  double precision,intent(in) :: x
  double precision,intent(out) :: rho,drhodr,vp,vs,Qkappa,Qmu

  ! local parameters
  double precision :: r,scaleval, x_moho
  double precision :: SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                      SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB

  ! compute real physical radius in meters
  r = x * R_PLANET

  ! check flags to make sure we correctly honor the discontinuities
  ! we use strict inequalities since r has been slightly changed in mesher

  ! gets corresponding Sohl&Spoon model radii
  call get_model_Sohl_radii(SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                            SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB)

  ! note: using stop statements, not exit_mpi() calls to avoid the need for MPI libraries when linking xcreate_header_file

  if (check_doubling_flag) then
    !
    !--- inner core
    !
    if (r >= 0.d0 .and. r < SOHL_RICB) then
      ! debug
      !if (myrank == 0) print *,'debug: model Sohl:',idoubling,IFLAG_INNER_CORE_NORMAL,IFLAG_MIDDLE_CENTRAL_CUBE, &
      !                            IFLAG_BOTTOM_CENTRAL_CUBE,IFLAG_TOP_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE

      ! checks with inner core flags
      if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
          idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
           stop 'wrong doubling flag for inner core point in model_Sohl()'
    !
    !--- outer core
    !
    else if (r > SOHL_RICB .and. r < SOHL_RCMB) then
      if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
        stop 'wrong doubling flag for outer core point in model_Sohl()'
    !
    !--- D" at the base of the mantle
    !
    else if (r > SOHL_RCMB .and. r < SOHL_RTOPDDOUBLEPRIME) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for D" point in model_Sohl()'
    !
    !--- mantle: from top of D" to d670
    !
    else if (r > SOHL_RTOPDDOUBLEPRIME .and. r < SOHL_R670) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for top D" - > d670 point in model_Sohl()'
    !
    !--- mantle: from d670 to d220
    !
    else if (r > SOHL_R670 .and. r < SOHL_R220) then
      if (idoubling /= IFLAG_670_220) &
        stop 'wrong doubling flag for d670 - > d220 point in model_Sohl()'
    !
    !--- mantle and crust: from d220 to MOHO and then to surface
    !
    else if (r > SOHL_R220) then
      if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
        stop 'wrong doubling flag for d220 - > Moho - > surface point in model_Sohl()'
    endif

  endif

! determines material properties

! note: this routine will be called for gravity/no-gravity to set density gradients and gravitational effects, also for
!       moon models other than Sohl, like 1D_case65TAY. in that case, we'll take the Sohl A model as reference.
!
  select case(REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_SOHL,REFERENCE_MODEL_CASE65TAY)
    ! Model A
    !
    !--- inner core
    !
    if (r >= 0.d0 .and. r <= SOHL_RICB) then
      drhodr=-13.4654d-3-2.0d0*1.6906d0*x-3.0d0*0.4225d0*x*x
      !rho=7.2942d0-13.4654d-3*x-1.6906d0*x*x-0.4225d0*x*x*x
      !vp=6.5395d0-22.5299d-3*x-2.3767d0*x*x-7.2716d-1*x*x*x
      !vs=2.4d0
      ! constant values for inner core
      ! todo: shouldn't we have therefore drhodr = 0 ? probably doesn't matter yet...
      rho=8.752d0
      vp=7.3d0
      vs=3d0
      Qmu=84.6d0
      Qkappa=1327.7d0
    !
    !--- outer core
    !
    else if (r > SOHL_RICB .and. r <= SOHL_RCMB) then
      drhodr=-13.4654d-3-2.0d0*1.6906d0*x-3.0d0*0.4225d0*x*x
      rho=7.2942d0-13.4654d-3*x-1.6906d0*x*x-0.4225d0*x*x*x
      vp=6.5395d0-22.5299d-3*x-2.3767d0*x*x-7.2716d-1*x*x*x
      vs=0.0d0
      Qmu=0.0d0
      Qkappa=57827.0d0
    !
    !--- D" at the base of the mantle
    !
    else if (r > SOHL_RCMB .and. r <= SOHL_RTOPDDOUBLEPRIME) then
      drhodr=-0.5184d0-2.0d0*3.2190d-2*x-3.0d0*1.3738d-1*x*x
      rho=4.5706d0-0.5184d0*x-3.2190d-2*x*x-1.3738d-1*x*x*x
      vp=11.4166d0-9.0421d-1*x-2.6380d0*x*x+9.4287d-1*x*x*x
      vs=6.7644d0-2.3139d0*x+1.7809d0*x*x-1.5312d0*x*x*x
      Qmu=312.0d0
      Qkappa=57827.0d0
    !
    !--- mantle: from top of D" to d670
    !
    else if (r > SOHL_RTOPDDOUBLEPRIME .and. r <= SOHL_R771) then
      drhodr=-0.5184d0-2.0d0*3.2190d-2*x-3.0d0*1.3738d-1*x*x
      rho=4.5706d0-0.5184d0*x-3.2190d-2*x*x-1.3738d-1*x*x*x
      vp=11.4166d0-9.0421d-1*x-2.6380d0*x*x+9.4287d-1*x*x*x
      vs=6.7644d0-2.3139d0*x+1.7809d0*x*x-1.5312d0*x*x*x
      Qmu=312.0d0
      Qkappa=57827.0d0

    else if (r > SOHL_R771 .and. r <= SOHL_R670) then
      drhodr=-5.8162d-1+2.0d0*1.7083d-1*x-3.0d0*3.0264d-1*x*x
      rho=4.4067d0-5.8162d-1*x+1.7083d-1*x*x-3.0264d-1*x*x*x
      vp=11.8365d0-4.1713d0*x+3.1366d0*x*x-2.5691d0*x*x*x
      vs=6.5172d0-1.8055d0*x+8.0800d-1*x*x-9.5676d-1*x*x*x
      Qmu=312.0d0
      Qkappa=57827.0d0
    !
    !--- mantle: above d670
    !
    else if (r > SOHL_R670 .and. r <= SOHL_R600) then
      drhodr=-0.8166d0+2.0d0*0.4844d0*x-3.0d0*0.4175d0*x*x
      rho=4.1601d0-0.8166d0*x+0.4844d0*x*x-0.4175d0*x*x*x
      vp=11.9300d0-4.8826d0*x+3.5016d0*x*x-2.5538d0*x*x*x
      vs=6.5847d0-2.5543d0*x+1.6024d0*x*x-1.2675d0*x*x*x
      Qmu=143.0d0
      Qkappa=57827.0d0

    else if (r > SOHL_R600 .and. r <= SOHL_R400) then
      drhodr=-0.8166d0+2.0d0*0.4844d0*x-3.0d0*0.4175d0*x*x
      rho=4.1601d0-0.8166d0*x+0.4844d0*x*x-0.4175d0*x*x*x
      vp=11.9300d0-4.8826d0*x+3.5016d0*x*x-2.5538d0*x*x*x
      vs=6.5847d0-2.5543d0*x+1.6024d0*x*x-1.2675d0*x*x*x
      Qmu=143.0d0
      Qkappa=57827.0d0

    else if (r > SOHL_R400 .and. r <= SOHL_R220) then
      drhodr=-0.8166d0+2.0d0*0.4844d0*x-3.0d0*0.4175d0*x*x
      rho=4.1601d0-0.8166d0*x+0.4844d0*x*x-0.4175d0*x*x*x
      vp=11.9300d0-4.8826d0*x+3.5016d0*x*x-2.5538d0*x*x*x
      vs=6.5847d0-2.5543d0*x+1.6024d0*x*x-1.2675d0*x*x*x
      Qmu=143.0d0
      Qkappa=57827.0d0

    else if (r > SOHL_R220 .and. r <= SOHL_R80) then
      drhodr=-3.4678d01+2.0d0*4.0167d01*x-3.0d0*1.5572d01*x*x
      rho=1.3572d01-3.4678d01*x+4.0167d01*x*x-1.5572d01*x*x*x
      vp=1.5559d01-1.7460d01*x+1.8004d01*x*x-8.1157d0*x*x*x
      vs=2.2395d01-5.7011d01*x+6.3447d01*x*x-2.4406d01*x*x*x
      Qmu=80.0d0
      Qkappa=57827.0d0
    else
      if (CRUSTAL .and. .not. SUPPRESS_CRUSTAL_MESH) then
        ! fill with mantle and later add CRUST2.0
        if (r > SOHL_R80) then
          ! density/velocity from mantle just below moho
          drhodr=-3.4678d01+2.0d0*4.0167d01*x-3.0d0*1.5572d01*x*x
          rho=1.3572d01-3.4678d01*x+4.0167d01*x*x-1.5572d01*x*x*x
          vp=1.5559d01-1.7460d01*x+1.8004d01*x*x-8.1157d0*x*x*x
          vs=2.2395d01-5.7011d01*x+6.3447d01*x*x-2.4406d01*x*x*x
          ! shear attenuation for R80 to surface
          Qmu=600.0d0
          Qkappa=57827.0d0
        endif
      else
        if (r > SOHL_R80 .and. r <= SOHL_RMOHO) then
          drhodr=1.7831d0+2.0d0*1.2168d0*x-3.0d0*1.7209d0*x*x
          rho=2.2095d0+1.7831d0*x+1.2168d0*x*x-1.7209d0*x*x*x
          vp=1.4372d01-1.3395d01*x+1.3353d01*x*x-6.3398d0*x*x*x
          vs=8.6113d0-1.1932d01*x+1.4301d01*x*x-6.5441d0*x*x*x
          Qmu=600.0d0
          Qkappa=57827.0d0

        else if (SUPPRESS_CRUSTAL_MESH) then
          ! DK DK extend the Moho up to the surface instead of the crust
          x_moho = SOHL_RMOHO / R_PLANET
          drhodr=1.7831d0+2.0d0*1.2168d0*x_moho-3.0d0*1.7209d0*x_moho*x_moho
          rho=2.2095d0+1.7831d0*x_moho+1.2168d0*x_moho*x_moho-1.7209d0*x_moho*x_moho*x_moho
          vp=1.4372d01-1.3395d01*x_moho+1.3353d01*x_moho*x_moho-6.3398d0*x_moho*x_moho*x_moho
          vs=8.6113d0-1.1932d01*x_moho+1.4301d01*x_moho*x_moho-6.5441d0*x_moho*x_moho*x_moho
          Qmu=600.0d0
          Qkappa=57827.0d0

        else if (r > SOHL_RMOHO .and. r <= SOHL_RMIDDLE_CRUST) then
          drhodr=1.2062d02-2.0d0*1.2140d02*x+3.0d0*4.0519d01*x*x
          rho=-3.6937d01+1.2062d02*x-1.2140d02*x*x+4.0519d01*x*x*x
          vp=5.0537d01-1.2368d02*x+1.2381d02*x*x-4.3012d01*x*x*x
          vs=-4.3819d01+1.4853d02*x-1.5121d02*x*x+5.0526d01*x*x*x
          Qmu=600.0d0
          Qkappa=57827.0d0

          ! same properties everywhere in Mars crust if we decide to
          ! define only one layer in the crust. Use the value at
          ! 60 km depth (r = 3330 km) in modSOHL from IPGP
          if (ONE_CRUST) then
            drhodr=0.0d0
            rho=2.815d0
            vp=7.747d0
            vs=4.069d0
            Qmu=600.0d0
            Qkappa=57827.0d0
          endif

        else if (r > SOHL_RMIDDLE_CRUST .and. r <= SOHL_ROCEAN) then
          drhodr=1.2062d02-2.0d0*1.2140d02*x+3.0d0*4.0519d01*x*x
          rho=-3.6937d01+1.2062d02*x-1.2140d02*x*x+4.0519d01*x*x*x
          vp=5.0537d01-1.2368d02*x+1.2381d02*x*x-4.3012d01*x*x*x
          vs=-4.3819d01+1.4853d02*x-1.5121d02*x*x+5.0526d01*x*x*x
          Qmu=600.0d0
          Qkappa=57827.0d0

        ! for density profile for gravity, we do not check that r <= R_MARS
        ! Obviously will not be used for Mars but could keep it for regolith
        else if (r > SOHL_ROCEAN) then
          drhodr=1.2062d02-2.0d0*1.2140d02*x+3.0d0*4.0519d01*x*x
          rho=-3.6937d01+1.2062d02*x-1.2140d02*x*x+4.0519d01*x*x*x
          vp=5.0537d01-1.2368d02*x+1.2381d02*x*x-4.3012d01*x*x*x
          vs=-4.3819d01+1.4853d02*x-1.5121d02*x*x+5.0526d01*x*x*x
          Qmu=600.0d0
          Qkappa=57827.0d0

        endif

      endif

    endif

  case (REFERENCE_MODEL_SOHL_B)
    ! Model B
    stop 'Error model_Sohl() routine: material properties for model B not implemented yet'

  case default
    stop 'Invalid reference model Sohl in model_Sohl()'
  end select

! note: we will set RHOAV and R_PLANET to MARS_RHOAV and mars radius in get_model_parameters
!       to have these non-dimensionalization factors account for a different planet. however, it might also work
!       with earth factor as long as we are consistent in dimensionalize/non-dimensionalize

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  drhodr = drhodr*1000.0d0/RHOAV
  rho = rho*1000.0d0/RHOAV
  vp = vp*1000.0d0/(R_PLANET*scaleval)
  vs = vs*1000.0d0/(R_PLANET*scaleval)

  end subroutine model_Sohl

!
!=====================================================================
!

! not used yet, mostly in case for AVS_DX mesh outputs

  subroutine Sohl_display_outer_core(x,rho,vp,vs,Qkappa,Qmu,idoubling)

! routine used for AVS or DX display of stability condition
! and number of points per wavelength only in the fluid outer core

  use constants, only: IFLAG_OUTER_CORE_NORMAL,PI,GRAV
  use shared_parameters, only: R_PLANET,RHOAV

  implicit none

! given a normalized radius x, gives the non-dimensionalized density rho,
! speeds vp and vs, and the quality factors Qkappa and Qmu

  integer, intent(in) :: idoubling
  double precision,intent(in) :: x
  double precision,intent(out) :: rho,vp,vs,Qkappa,Qmu

  ! local parameters
  double precision :: scaleval

  if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
    stop 'wrong doubling flag for outer core point in Sohl_display_outer_core()'

!
!--- outer core
!
  rho = 7.2942d0-13.4654d-3*x-1.6906d0*x*x-0.4225d0*x*x*x
  vp = 6.5395d0-22.5299d-3*x-2.3767d0*x*x-7.2716d-1*x*x*x
  vs = 0.0d0
  Qmu = 0.0d0
  Qkappa = 57827.0d0

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  rho = rho*1000.0d0/RHOAV
  vp = vp*1000.0d0/(R_PLANET*scaleval)
  vs = vs*1000.0d0/(R_PLANET*scaleval)

  end subroutine Sohl_display_outer_core

!
!=====================================================================
!

  subroutine Sohl_density(x,rho)

  use constants, only: R_MARS,REFERENCE_MODEL_SOHL,REFERENCE_MODEL_SOHL_B,REFERENCE_MODEL_CASE65TAY

  use shared_parameters, only: RHOAV,REFERENCE_1D_MODEL,ONE_CRUST

  use model_sohl_par

  implicit none

  double precision, intent(in) :: x
  double precision, intent(out) :: rho

  ! local parameters
  double precision :: r
  double precision :: SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                      SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB

  ! note: this routine will also be called in case ellipticity is set and the reference model is REFERENCE_MODEL_CASE65TAY.
  !       to compute ellipicity, we take the density profile of the Sohl model A reference model in that case.

  ! compute real physical radius in meters
  r = x * R_MARS

  ! gets corresponding Sohl&Spoon model radii
  call get_model_Sohl_radii(SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                            SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB)

  ! calculates density according to radius
  select case(REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_SOHL,REFERENCE_MODEL_CASE65TAY)
    ! Model A
    if (r <= SOHL_RICB) then
      rho = 7.2942d0-13.4654d-3*x-1.6906d0*x*x-0.4225d0*x*x*x
    else if (r > SOHL_RICB .and. r <= SOHL_RCMB) then
      rho = 7.2942d0-13.4654d-2*x-1.6906d0*x*x-0.4225d0*x*x*x
    else if (r > SOHL_RCMB .and. r <= SOHL_RTOPDDOUBLEPRIME) then
      rho = 4.5706d0-0.5184d0*x-3.2190d-2*x*x-1.3738d-1*x*x*x
    else if (r > SOHL_RTOPDDOUBLEPRIME .and. r <= SOHL_R771) then
      rho = 4.5706d0-0.5184d0*x-3.2190d-2*x*x-1.3738d-1*x*x*x
    else if (r > SOHL_R771 .and. r <= SOHL_R670) then
      rho = 4.4067d0-5.8162d-1*x+1.7083d-1*x*x-3.0264d-1*x*x*x
    else if (r > SOHL_R670 .and. r <= SOHL_R600) then
      rho = 4.1601d0-0.8166d0*x+0.4844d0*x*x-0.4175d0*x*x*x
    else if (r > SOHL_R600 .and. r <= SOHL_R400) then
      rho = 4.1601d0-0.8166d0*x+0.4844d0*x*x-0.4175d0*x*x*x
    else if (r > SOHL_R400 .and. r <= SOHL_R220) then
      rho = 4.1601d0-0.8166d0*x+0.4844d0*x*x-0.4175d0*x*x*x
    else if (r > SOHL_R220 .and. r <= SOHL_R80) then
      rho = 1.3572d01-3.4678d01*x+4.0167d01*x*x-1.5572d01*x*x*x
    else
      if (r > SOHL_R80 .and. r <= SOHL_RMOHO) then
        rho = 2.2095d0+1.7831d0*x+1.2168d0*x*x-1.7209d0*x*x*x
      else if (r > SOHL_RMOHO .and. r <= SOHL_RMIDDLE_CRUST) then
        if (ONE_CRUST) then
          rho = 2.815d0  ! constant upper crust value
        else
          rho = -3.6937d01+1.2062d02*x-1.2140d02*x*x+4.0519d01*x*x*x
        endif
      else if (r > SOHL_RMIDDLE_CRUST .and. r <= SOHL_ROCEAN) then
        rho = -3.6937d01+1.2062d02*x-1.2140d02*x*x+4.0519d01*x*x*x
      else if (r > SOHL_ROCEAN) then
        rho = -3.6937d01+1.2062d02*x-1.2140d02*x*x+4.0519d01*x*x*x ! no ocean, same as upper crust
      endif
    endif

  case (REFERENCE_MODEL_SOHL_B)
    ! Model B
    stop 'Error Sohl_density() routine: material properties for model B not implemented yet'

  case default
    stop 'Invalid 1D reference model in routine Sohl_density()'
  end select

  ! non-dimensionalize
  rho = rho * 1000.0d0 / RHOAV

  end subroutine Sohl_density

!
!=====================================================================
!

  subroutine get_model_Sohl_radii(SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                                  SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB)

! determines radii for Model A or Model B

  use constants, only: REFERENCE_MODEL_SOHL,REFERENCE_MODEL_SOHL_B,REFERENCE_MODEL_CASE65TAY

  use shared_parameters, only: REFERENCE_1D_MODEL

  use model_sohl_par

  implicit none

  double precision,intent(out) :: SOHL_RMOHO,SOHL_R80,SOHL_R220,SOHL_R400,SOHL_R600,SOHL_R670, &
                                  SOHL_R771,SOHL_RTOPDDOUBLEPRIME,SOHL_RCMB

  ! note: this routine will also be called in case ellipticity is set and the reference model is REFERENCE_MODEL_CASE65TAY.
  !       to compute ellipicity, we take the density profile of the Sohl model A reference model in that case.

  select case(REFERENCE_1D_MODEL)
  case (REFERENCE_MODEL_SOHL,REFERENCE_MODEL_CASE65TAY)
    ! Model A
    SOHL_RMOHO = SOHL_A_RMOHO
    SOHL_R80 = SOHL_A_R80
    SOHL_R220 = SOHL_A_R220
    SOHL_R400 = SOHL_A_R400
    SOHL_R600 = SOHL_A_R600
    SOHL_R670 = SOHL_A_R670
    SOHL_R771 = SOHL_A_R771
    SOHL_RTOPDDOUBLEPRIME = SOHL_A_RTOPDDOUBLEPRIME
    SOHL_RCMB = SOHL_A_RCMB

  case (REFERENCE_MODEL_SOHL_B)
    ! Model B
    SOHL_RMOHO = SOHL_B_RMOHO
    SOHL_R80 = SOHL_B_R80
    SOHL_R220 = SOHL_B_R220
    SOHL_R400 = SOHL_B_R400
    SOHL_R600 = SOHL_B_R600
    SOHL_R670 = SOHL_B_R670
    SOHL_R771 = SOHL_B_R771
    SOHL_RTOPDDOUBLEPRIME = SOHL_B_RTOPDDOUBLEPRIME
    SOHL_RCMB = SOHL_B_RCMB

  case default
    stop 'Invalid 1D reference model in routine get_model_Sohl_radii()'
  end select

  end subroutine get_model_Sohl_radii



!=====================================================================
!
! Mars 1D model - generic model from file in DATA/mars/mars_1D.dat
!
!=====================================================================

! mars 1D model defined by table in file mars_1D.dat
!
! similar implementation and file format as VPREMOON model (in DATA/moon/vpremoon.dat)

  module model_mars_1D_par

  implicit none

  ! radius at physical surface
  double precision, parameter :: MARS_1D_RSURFACE = 1737100.d0

  ! radii
  ! no ocean, same as surface radius (Physical surface)
  double precision, parameter :: MARS_1D_ROCEAN = MARS_1D_RSURFACE ! no ocean
  double precision, parameter :: MARS_1D_RMIDDLE_CRUST = 1725100.d0
  ! Crust-mantle boundary: at 28km depth corresponds to seismic model, not geodesic model
  double precision, parameter :: MARS_1D_RMOHO = 1709100.d0
  ! Rheological lithosphere
  double precision, parameter :: MARS_1D_R80  =  1607100.d0
  ! Thermal lithosphere
  double precision, parameter :: MARS_1D_R220 = 1517100.d0
  double precision, parameter :: MARS_1D_R400 = 1337100.d0
  double precision, parameter :: MARS_1D_R600 = 1137100.d0
  ! alpha-olivine-beta-spinel transition
  double precision, parameter :: MARS_1D_R670 = 1067100.d0
  ! beta-spinel-gamma-spinel transition
  double precision, parameter :: MARS_1D_R771 =  966100.d0
  ! lower thermal boundary layer
  double precision, parameter :: MARS_1D_RTOPDDOUBLEPRIME = 480000.d0
  ! Core-mantle boundary
  double precision, parameter :: MARS_1D_RCMB = 380000.d0
  ! inner core radius
  double precision, parameter :: MARS_1D_RICB = 240000.d0

  ! densities [g/cm^3]
  double precision, parameter :: MARS_1D_RHO_OCEANS    = 1020.d0  ! will not be used
  double precision, parameter :: MARS_1D_RHO_TOP_OC    = 5171.d0  ! outer core density fluid outer core
  double precision, parameter :: MARS_1D_RHO_BOTTOM_OC = 5171.d0  ! outer core density (constant value)

  ! number of layers
  integer, parameter :: NR_mars_1D_layers = 78  ! check with mars_1D.dat in DATA/mars/ folder

  ! model arrays
  double precision, dimension(:),allocatable :: mars_1D_radius,mars_1D_density, &
    mars_1D_vp,mars_1D_vs,mars_1D_Qkappa,mars_1D_Qmu

  ! storage of initial Qmu profile (used for attenuation setup)
  double precision, dimension(:),allocatable :: MARS_1D_Qmu_original

  end module model_mars_1D_par

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_mars_1D(x,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,CRUSTAL, &
                            check_doubling_flag,iregion_code)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,ONE_CRUST

  use model_mars_1D_par

  implicit none

! input:
! radius r: meters

! output:
! density rho: kg/m^3
! compressional wave speed vp: km/s
! shear wave speed vs: km/s

! note: although we provide Qkappa, it is not further used in the code.
!       SPECFEM3D_GLOBE only implements shear attenuation at the moment.

  double precision,intent(in) :: x
  double precision,intent(inout) :: rho,drhodr
  double precision,intent(inout) :: vp,vs,Qmu,Qkappa
  integer,intent(in) :: idoubling
  integer,intent(in) :: iregion_code
  logical,intent(in) :: CRUSTAL,check_doubling_flag

  ! local parameters
  double precision :: r,frac,scaleval,dr,drho
  double precision :: vp_tmp,vs_tmp,rho_tmp,Qmu_tmp,Qkappa_tmp
  integer :: i

  ! first allocates and defines model arrays (will be done only once)
  if (.not. allocated(MARS_1D_radius)) then
    ! all processes will define same parameters
    call define_model_mars_1D()

    ! ONE_CRUST uses a single layer value for the crust
    if (ONE_CRUST) then
      ! upper crust is in (3:4)
      rho_tmp = mars_1D_density(NR_mars_1D_layers-3)
      vp_tmp = mars_1D_vp(NR_mars_1D_layers-3)
      vs_tmp = mars_1D_vs(NR_mars_1D_layers-3)
      Qkappa_tmp = mars_1D_Qkappa(NR_mars_1D_layers-3)
      Qmu_tmp = mars_1D_Qmu(NR_mars_1D_layers-3)
      ! assign all crust (1:6) values to upper crust
      mars_1D_density(NR_mars_1D_layers-5:NR_mars_1D_layers) = rho_tmp
      mars_1D_vp(NR_mars_1D_layers-5:NR_mars_1D_layers) = vp_tmp
      mars_1D_vs(NR_mars_1D_layers-5:NR_mars_1D_layers) = vs_tmp
      mars_1D_Qkappa(NR_mars_1D_layers-5:NR_mars_1D_layers) = Qkappa_tmp
      mars_1D_Qmu(NR_mars_1D_layers-5:NR_mars_1D_layers) = Qmu_tmp
    endif

    ! in case an external crustal model will be superimposed on top, we extend mantle values to the surface
    !
    ! strip the crust and replace it by mantle if we use an external crustal model
    ! note: assumes that the crust is given by 6 layers
    !       see in mars_1D.dat: regolith layer (1:2), upper crust (3:4), lower crust (5:6), the mantle
    if (SUPPRESS_CRUSTAL_MESH .or. CRUSTAL) then
      do i = NR_mars_1D_layers-5,NR_mars_1D_layers
        mars_1D_density(i) = mars_1D_density(NR_mars_1D_layers-6)
        mars_1D_vp(i) = mars_1D_vp(NR_mars_1D_layers-6)
        mars_1D_vs(i) = mars_1D_vs(NR_mars_1D_layers-6)
        mars_1D_Qkappa(i) = mars_1D_Qkappa(NR_mars_1D_layers-6)
        mars_1D_Qmu(i) = mars_1D_Qmu(NR_mars_1D_layers-6)
      enddo
    endif
  endif

  ! compute real physical radius in meters
  r = x * R_PLANET

  ! check flags to make sure we correctly honor the discontinuities
  ! we use strict inequalities since r has been slightly changed in mesher

  ! note: using stop statements, not exit_mpi() calls to avoid the need for MPI libraries when linking xcreate_header_file

  if (check_doubling_flag) then
    !
    !--- inner core
    !
    if (r >= 0.d0 .and. r < MARS_1D_RICB) then
      ! debug
      !if (myrank == 0) print *,'debug: model Sohl:',idoubling,IFLAG_INNER_CORE_NORMAL,IFLAG_MIDDLE_CENTRAL_CUBE, &
      !                            IFLAG_BOTTOM_CENTRAL_CUBE,IFLAG_TOP_CENTRAL_CUBE,IFLAG_IN_FICTITIOUS_CUBE

      ! checks with inner core flags
      if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
          idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
          idoubling /= IFLAG_IN_FICTITIOUS_CUBE) &
           stop 'wrong doubling flag for inner core point in model_mars_1D()'
    !
    !--- outer core
    !
    else if (r > MARS_1D_RICB .and. r < MARS_1D_RCMB) then
      if (idoubling /= IFLAG_OUTER_CORE_NORMAL) &
        stop 'wrong doubling flag for outer core point in model_mars_1D()'
    !
    !--- D" at the base of the mantle
    !
    else if (r > MARS_1D_RCMB .and. r < MARS_1D_RTOPDDOUBLEPRIME) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for D" point in model_mars_1D()'
    !
    !--- mantle: from top of D" to d670
    !
    else if (r > MARS_1D_RTOPDDOUBLEPRIME .and. r < MARS_1D_R670) then
      if (idoubling /= IFLAG_MANTLE_NORMAL) &
        stop 'wrong doubling flag for top D" - > d670 point in model_mars_1D()'
    !
    !--- mantle: from d670 to d220
    !
    else if (r > MARS_1D_R670 .and. r < MARS_1D_R220) then
      if (idoubling /= IFLAG_670_220) &
        stop 'wrong doubling flag for d670 - > d220 point in model_mars_1D()'
    !
    !--- mantle and crust: from d220 to MOHO and then to surface
    !
    else if (r > MARS_1D_R220) then
      if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) &
        stop 'wrong doubling flag for d220 - > Moho - > surface point in model_mars_1D()'
    endif
  endif

  ! searches radius index
  i = 1
  do while(r >= MARS_1D_radius(i) .and. i /= NR_mars_1D_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_mars_1D_layers) i = NR_mars_1D_layers

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (iregion_code /= 0) then
    ! check
    if (iregion_code /= IREGION_INNER_CORE .and. &
        iregion_code /= IREGION_OUTER_CORE .and. &
        iregion_code /= IREGION_CRUST_MANTLE) stop 'Invalid iregion_code in model_mars_1D() routine'

    ! inner core range [1:2]
    if (iregion_code == IREGION_INNER_CORE .and. i > 2) i = 2

    ! outer core range [3:4]
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_OUTER_CORE .and. i <= 3) i = 3 + 1
    if (iregion_code == IREGION_OUTER_CORE .and. i > 4) i = 4

    ! crust/mantle range [5:78] NR_mars_1D_layers == 78
    ! due to the interpolation below, we add +1 to the bottom
    if (iregion_code == IREGION_CRUST_MANTLE .and. i <= 5) i = 5 + 1
    if (iregion_code == IREGION_CRUST_MANTLE .and. i > NR_mars_1D_layers) i = NR_mars_1D_layers
  endif

  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    drhodr = 0.d0
    rho = mars_1D_density(i)
    vp = mars_1D_vp(i)
    vs = mars_1D_vs(i)
    Qmu = mars_1D_Qmu(i)
    Qkappa = mars_1D_Qkappa(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = mars_1D_radius(i) - mars_1D_radius(i-1)
    if (abs(dr) > 1.d-9) then
      drho = (mars_1D_density(i)-mars_1D_density(i-1))
      drhodr = drho/dr
      frac = (r - mars_1D_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      drhodr = 0.d0
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = mars_1D_density(i-1) + frac * (mars_1D_density(i)-mars_1D_density(i-1))
    vp = mars_1D_vp(i-1) + frac * (mars_1D_vp(i)-mars_1D_vp(i-1))
    vs = mars_1D_vs(i-1) + frac * (mars_1D_vs(i)-mars_1D_vs(i-1))
    Qmu = mars_1D_Qmu(i-1) + frac * (mars_1D_Qmu(i)-mars_1D_Qmu(i-1))
    Qkappa = mars_1D_Qkappa(i-1) + frac * (mars_1D_Qkappa(i)-mars_1D_Qkappa(i-1))
  endif

! make sure Vs is zero in the outer core even if roundoff errors on depth
! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if (iregion_code == IREGION_OUTER_CORE) then
    vs = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval = dsqrt(PI*GRAV*RHOAV)
  rho = rho*1000.0d0/RHOAV
  vp = vp*1000.0d0/(R_PLANET*scaleval)
  vs = vs*1000.0d0/(R_PLANET*scaleval)

  end subroutine model_mars_1D

!
!-------------------------------------------------------------------------------------------------
!

  subroutine define_model_mars_1D()

  use constants, only: myrank,SUPPRESS_CRUSTAL_MESH,IIN
  use shared_parameters, only: RICB,RCMB,ROCEAN
  use model_mars_1D_par

  implicit none

  ! local parameters
  integer :: nlines,ir,ier
  double precision :: r,vp,vs,rho
  double precision :: Qp,Qs,Qkappa,Qmu
  double precision :: Qkappa_inv,Qp_inv,dr
  character(len=256) :: datafile,line

  double precision, parameter :: TOL = 1.d-9
  ! debugging
  logical, parameter :: DEBUG = .false.

  ! file name
  datafile = 'DATA/mars/mars_1D.dat'

  ! reads in data file
  open(IIN,file=trim(datafile),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(datafile)
    stop 'Error opening mars_1D.dat file'
  endif

  ! checks number of line
  nlines = 0
  do while(ier == 0)
    read(IIN,"(a256)",iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    line = adjustl(line)        ! suppress leading white spaces, if any
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    nlines = nlines + 1
  enddo
  rewind(IIN)

  ! safety check
  if (nlines /= NR_mars_1D_layers) then
    print *,'Error: number of data entries ',nlines,' should match number of layers ',NR_mars_1D_layers,' for mars_1D model'
    print *,'Please check file: ',trim(datafile)
    stop 'Invalid number of data line entries in mars_1D.dat'
  endif

  ! allocates model arrays
  allocate(mars_1D_radius(NR_mars_1D_layers), &
           mars_1D_density(NR_mars_1D_layers), &
           mars_1D_vp(NR_mars_1D_layers), &
           mars_1D_vs(NR_mars_1D_layers), &
           mars_1D_Qkappa(NR_mars_1D_layers), &
           mars_1D_Qmu(NR_mars_1D_layers), &
           mars_1D_Qmu_original(NR_mars_1D_layers), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating mars_1D arrays'
  mars_1D_radius(:) = 0.d0; mars_1D_density(:) = 0.d0
  mars_1D_vp(:) = 0.d0; mars_1D_vs(:) = 0.d0
  mars_1D_Qkappa(:) = 0.d0; mars_1D_Qmu(:) = 0.d0; mars_1D_Qmu_original(:) = 0.d0

  ! read in layering
  ! note: table 6 starts from top surface to inner core
  !       for the model arrays here, we choose an opposite order, with index 1 being the center of Moon
  !       and then increasing radius.
  nlines = 0
  do while(ier == 0)
    ! read line
    read(IIN,"(a256)",iostat=ier) line
    if (ier /= 0) exit

    ! skip empty/comment lines
    line = adjustl(line)        ! suppress leading white spaces, if any
    if (len_trim(line) == 0) cycle
    if (line(1:1) == '#' .or. line(1:1) == '!') cycle

    ! read values from line
    ! format: # radius(km) #vp(km/s) #vs(km/s) #density(kg/cm^3) #Qp  #Qs
    read(line,*,iostat=ier) r,vp,vs,rho,Qp,Qs

    if (ier /= 0) then
      print *,'Error: reading mars_1D.dat data, on line ',nlines
      print *,'       data line: ',trim(line)
      print *,'Please check data file mars_1D.dat'
      stop 'Error reading data in mars_1D data'
    endif

    if (r < 0.d0) stop 'Error found negative r value'
    if (vs < 0.d0) stop 'Error found negative vs value'
    if (vp <= 0.d0) stop 'Error found negative or zero vp value'
    if (Qp <= 0.d0) stop 'Error found negative or zero Qp value'
    if (Qs <= 0.d0) stop 'Error found negative or zero Qs value'

    ! converts Qp,Qs to Qkappa,Qmu
    ! see for example Dahlen and Tromp 1998, formulas 9.59 and 9.60
    Qmu = Qs

    if (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2) <= 0.d0) then
      print *,'Error: invalid vp,vs = ',vp,vs,' leads to factor ',(1.d0 - 4.d0/3.d0*(vs**2)/(vp**2))
      print *,'       data line: ',trim(line)
      stop 'Error invalid vp,vs values'
    endif

    ! Qkappa converted from (Qp,Qs)
    ! note: in case Qp == Qs then it follows that Qkappa = Qs
    Qkappa_inv = (1.d0/Qp - 4.d0/3.d0*(vs**2)/(vp**2)/Qs) / (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2))
    if (Qkappa_inv <= 0.d0) stop 'Error got negative or zero 1/Qkappa value'

    Qkappa = 1.d0/Qkappa_inv

    ! check back with Qp
    if (DEBUG) then
      Qp_inv = (1.d0 - 4.d0/3.d0*(vs**2)/(vp**2)) * 1.d0/Qkappa + 4.d0/3.d0*(vs**2)/(vp**2) * 1.d0/Qmu
      ! debug output
      if (myrank == 0) &
        print *,'debug: mars_1D values ',r,vp,vs,rho,Qp,Qs
        print *,'debug: Qkappa,Qmu = ',Qkappa,Qmu,' check Qp',Qp,1.d0/Qp_inv
      ! check
      if (abs(Qp - 1.d0/Qp_inv) > 1.d-9) stop 'Error invalid Qp check'
    endif

    ! increases layer number count
    nlines = nlines + 1

    ! define all the values in the model
    ! change ordering from core to top
    ir = NR_mars_1D_layers - (nlines-1)
    if (ir < 1) stop 'Error reading in mars_1D data lines'

    ! we need radius in [m], density in [kg/cm^3], vp,vs in [km/s] and Qkappa and Qmu
    mars_1D_radius(ir) = r * 1000.d0
    mars_1D_density(ir) = rho
    mars_1D_vp(ir) = vp
    mars_1D_vs(ir) = vs
    mars_1D_Qkappa(ir) = Qkappa
    mars_1D_Qmu(ir) = Qmu
  enddo
  close(IIN)

  ! debug
  if (DEBUG) then
    if (myrank == 0) then
      print *,'debug: model mars_1D'
      print *,'format: #line #r #vp #vs #rho #Qkappa #Qmu'
      do ir = 1,NR_mars_1D_layers
        print *,ir,sngl(mars_1D_radius(ir)),sngl(mars_1D_vp(ir)),sngl(mars_1D_vs(ir)),sngl(mars_1D_density(ir)), &
                   sngl(mars_1D_Qkappa(ir)),sngl(mars_1D_Qmu(ir))
      enddo
      print *
    endif
  endif

  ! checks radius is increasing
  do ir = 2,NR_mars_1D_layers
    dr = mars_1D_radius(ir) - mars_1D_radius(ir-1)
    ! slightly negative steps should not occur
    if (dr < - TOL) then
      print *,'Error layering inversed: radius should decrease in file mars_1D.dat and start from top down to center.'
      print *,'  wrong layering: ',ir,mars_1D_radius(ir),'and ',ir-1,mars_1D_radius(ir-1)
      stop 'Error layering in mars_1D.dat has wrong ordering'
    endif
  enddo

  ! checks data table
  ! ICB at radius(33)
  if (abs(RICB - mars_1D_radius(2)) > TOL) &
    stop 'Error: mars_1D radius RICB and model radius index not matching'
  if (abs(RCMB - mars_1D_radius(4)) > TOL) &
    stop 'Error: mars_1D radius RCMB and model radius index not matching'
  if (abs(ROCEAN - mars_1D_radius(NR_mars_1D_layers)) > TOL) &
    stop 'Error: mars_1D radius ROCEAN and model radius index not matching'

  ! stores Qmu original values (without CRUSTAL modifications) for further use in attenuation routine
  mars_1D_Qmu_original(:) = mars_1D_Qmu(:)

  end subroutine define_model_mars_1D

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_mars_1D_density(x,rho)

  use constants
  use shared_parameters, only: R_PLANET,RHOAV,ONE_CRUST

  use model_mars_1D_par

  implicit none

  double precision,intent(in) :: x
  double precision,intent(out) :: rho

  ! local parameters
  integer :: i
  double precision :: r,dr,frac
  double precision :: vp_tmp,vs_tmp,rho_tmp,Qmu_tmp,Qkappa_tmp

  ! first allocates and defines model arrays (will be done only once)
  if (.not. allocated(mars_1D_radius)) then
    ! all processes will define same parameters
    call define_model_mars_1D()

    ! ONE_CRUST uses a single layer value for the crust
    if (ONE_CRUST) then
      ! upper crust is in (3:4)
      rho_tmp = mars_1D_density(NR_mars_1D_layers-3)
      vp_tmp = mars_1D_vp(NR_mars_1D_layers-3)
      vs_tmp = mars_1D_vs(NR_mars_1D_layers-3)
      Qkappa_tmp = mars_1D_Qkappa(NR_mars_1D_layers-3)
      Qmu_tmp = mars_1D_Qmu(NR_mars_1D_layers-3)
      ! assign all crust (1:6) values to upper crust
      mars_1D_density(NR_mars_1D_layers-5:NR_mars_1D_layers) = rho_tmp
      mars_1D_vp(NR_mars_1D_layers-5:NR_mars_1D_layers) = vp_tmp
      mars_1D_vs(NR_mars_1D_layers-5:NR_mars_1D_layers) = vs_tmp
      mars_1D_Qkappa(NR_mars_1D_layers-5:NR_mars_1D_layers) = Qkappa_tmp
      mars_1D_Qmu(NR_mars_1D_layers-5:NR_mars_1D_layers) = Qmu_tmp
    endif
  endif

  ! initializes
  rho = 0.d0

  ! compute real physical radius in meters
  r = x * R_PLANET

  i = 1
  do while(r >= mars_1D_radius(i) .and. i /= NR_mars_1D_layers)
    i = i + 1
  enddo

  ! just double check
  if (i > NR_mars_1D_layers) i = NR_mars_1D_layers

  ! note: this routine helps to sample the density profile of the 1D model
  !       we might assign the "upper" value for points on interfaces like ICB, CMB etc.
  !       this effect should not matter much though since the sampling rate is quite high along the profile.

  ! calculates density according to radius
  ! interpolates model values
  if (i == 1) then
    ! no interpolation with previous layer possible
    rho = mars_1D_density(i)
  else
    ! interpolate from radius(i-1) to r using the values at i-1 and i
    dr = mars_1D_radius(i) - mars_1D_radius(i-1)
    ! fraction for linear interpolation
    if (abs(dr) > 1.d-9) then
      frac = (r - mars_1D_radius(i-1))/dr
    else
      ! both current and lower layer have same radius
      frac = 1.d0 ! will take values from current layer i
    endif
    rho = mars_1D_density(i-1) + frac * (mars_1D_density(i)-mars_1D_density(i-1))
  endif

  ! non-dimensionalize
  rho = rho*1000.0d0/RHOAV

  end subroutine model_mars_1D_density
