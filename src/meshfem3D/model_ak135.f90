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

!--------------------------------------------------------------------------------------------------
!
! Modified AK135 model:
!
! Spherically symmetric isotropic AK135 model (Kennett et al., 1995).
! modified to use the density and Q attenuation models of Montagner and Kennett (1995).
! That modified model is traditionally called AK135-F,
! see http://rses.anu.edu.au/seismology/ak135/ak135f.html for more details.
! As we do not want to use the 300 m-thick mud layer from that model nor the ocean layer,
! above the d120 discontinuity we switch back to the classical AK135 model of Kennett et al. (1995),
! i.e., we use AK135-F below and AK135 above.

! B. L. N. Kennett, E. R. Engdahl and R. Buland,
! Constraints on seismic velocities in the Earth from traveltimes,
! Geophysical Journal International, volume 122, issue 1, pages 108-124 (1995),
! DOI: 10.1111/j.1365-246X.1995.tb03540.x
!--------------------------------------------------------------------------------------------------

! J. P. Montagner and B. L. N. Kennett,
! How to reconcile body-wave and normal-mode reference Earth models?,
! Geophysical Journal International, volume 122, issue 1, pages 229-248 (1995)

!! DK DK values below entirely checked and fixed by Dimitri Komatitsch in December 2012.

!--------------------------------------------------------------------------------------------------

  module model_ak135_par

  ! number of layers in DATA/ak135/ak135.dat
  integer, parameter :: NR_AK135F_NO_MUD = 136

  ! model_ak135_variables
  double precision, dimension(:), allocatable :: &
    Mak135_V_radius_ak135,Mak135_V_density_ak135, &
    Mak135_V_vp_ak135,Mak135_V_vs_ak135, &
    Mak135_V_Qkappa_ak135,Mak135_V_Qmu_ak135

  end module model_ak135_par

!
!--------------------------------------------------------------------------------------------------
!


  subroutine model_ak135_broadcast(CRUSTAL)

! standard routine to setup model

  use constants, only: myrank
  use model_ak135_par

  implicit none

  logical :: CRUSTAL

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(Mak135_V_radius_ak135(NR_AK135F_NO_MUD), &
          Mak135_V_density_ak135(NR_AK135F_NO_MUD), &
          Mak135_V_vp_ak135(NR_AK135F_NO_MUD), &
          Mak135_V_vs_ak135(NR_AK135F_NO_MUD), &
          Mak135_V_Qkappa_ak135(NR_AK135F_NO_MUD), &
          Mak135_V_Qmu_ak135(NR_AK135F_NO_MUD), &
          stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating Mak135_V arrays')

  ! all processes will define same parameters
  call define_model_ak135(CRUSTAL)

  end subroutine model_ak135_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_ak135(x,rho,vp,vs,Qkappa,Qmu,iregion_code)

  use constants
  use model_ak135_par

  implicit none

! input:
! radius r: meters

! output:
! density rho: kg/m^3
! compressional wave speed vp: km/s
! shear wave speed vs: km/s

  double precision :: x,rho,vp,vs,Qmu,Qkappa
  integer :: iregion_code

  ! local parameters
  double precision :: r,frac,scaleval
  integer :: i

  ! compute real physical radius in meters
  r = x * R_EARTH

  i = 1
  do while(r >= Mak135_V_radius_ak135(i) .and. i /= NR_AK135F_NO_MUD)
    i = i + 1
  enddo

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (iregion_code == IREGION_INNER_CORE .and. i > 24) i = 24

  if (iregion_code == IREGION_OUTER_CORE .and. i < 26) i = 26
  if (iregion_code == IREGION_OUTER_CORE .and. i > 69) i = 69

  if (iregion_code == IREGION_CRUST_MANTLE .and. i < 71) i = 71

  if (i == 1) then
    rho = Mak135_V_density_ak135(i)
    vp = Mak135_V_vp_ak135(i)
    vs = Mak135_V_vs_ak135(i)
    Qmu = Mak135_V_Qmu_ak135(i)
    Qkappa = Mak135_V_Qkappa_ak135(i)
  else

! interpolate from radius_ak135(i-1) to r using the values at i-1 and i
    frac = (r-Mak135_V_radius_ak135(i-1))/(Mak135_V_radius_ak135(i)-Mak135_V_radius_ak135(i-1))

    rho = Mak135_V_density_ak135(i-1) + frac * (Mak135_V_density_ak135(i)-Mak135_V_density_ak135(i-1))
    vp = Mak135_V_vp_ak135(i-1) + frac * (Mak135_V_vp_ak135(i)-Mak135_V_vp_ak135(i-1))
    vs = Mak135_V_vs_ak135(i-1) + frac * (Mak135_V_vs_ak135(i)-Mak135_V_vs_ak135(i-1))
    Qmu = Mak135_V_Qmu_ak135(i-1) + frac * (Mak135_V_Qmu_ak135(i)-Mak135_V_Qmu_ak135(i-1))
    Qkappa = Mak135_V_Qkappa_ak135(i-1) + frac * (Mak135_V_Qkappa_ak135(i)-Mak135_V_Qkappa_ak135(i-1))

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
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)

  end subroutine model_ak135

!-------------------

  subroutine define_model_ak135(USE_EXTERNAL_CRUSTAL_MODEL)

  use constants
  use model_ak135_par

  implicit none

  logical :: USE_EXTERNAL_CRUSTAL_MODEL

! define all the values in the model

  Mak135_V_radius_ak135(  1) =  0.d0
  Mak135_V_radius_ak135(  2) =   50710.d0
  Mak135_V_radius_ak135(  3) =   101430.d0
  Mak135_V_radius_ak135(  4) =   152140.d0
  Mak135_V_radius_ak135(  5) =   202850.d0
  Mak135_V_radius_ak135(  6) =   253560.d0
  Mak135_V_radius_ak135(  7) =   304280.d0
  Mak135_V_radius_ak135(  8) =   354990.d0
  Mak135_V_radius_ak135(  9) =   405700.d0
  Mak135_V_radius_ak135( 10) =   456410.d0
  Mak135_V_radius_ak135( 11) =   507130.d0
  Mak135_V_radius_ak135( 12) =   557840.d0
  Mak135_V_radius_ak135( 13) =   659260.d0
  Mak135_V_radius_ak135( 14) =   710000.d0
  Mak135_V_radius_ak135( 15) =   760690.d0
  Mak135_V_radius_ak135( 16) =   811400.d0
  Mak135_V_radius_ak135( 17) =   862110.d0
  Mak135_V_radius_ak135( 18) =   912830.d0
  Mak135_V_radius_ak135( 19) =   963540.d0
  Mak135_V_radius_ak135( 20) =   1014250.d0
  Mak135_V_radius_ak135( 21) =   1064960.d0
  Mak135_V_radius_ak135( 22) =   1115680.d0
  Mak135_V_radius_ak135( 23) =   1166390.d0
  Mak135_V_radius_ak135( 24) =   1217500.d0
  Mak135_V_radius_ak135( 25) =   1217500.d0
  Mak135_V_radius_ak135( 26) =   1267430.d0
  Mak135_V_radius_ak135( 27) =   1317760.d0
  Mak135_V_radius_ak135( 28) =   1368090.d0
  Mak135_V_radius_ak135( 29) =   1418420.d0
  Mak135_V_radius_ak135( 30) =   1468760.d0
  Mak135_V_radius_ak135( 31) =   1519090.d0
  Mak135_V_radius_ak135( 32) =   1569420.d0
  Mak135_V_radius_ak135( 33) =   1670080.d0
  Mak135_V_radius_ak135( 34) =   1720410.d0
  Mak135_V_radius_ak135( 35) =   1770740.d0
  Mak135_V_radius_ak135( 36) =   1821070.d0
  Mak135_V_radius_ak135( 37) =   1871400.d0
  Mak135_V_radius_ak135( 38) =   1921740.d0
  Mak135_V_radius_ak135( 39) =   1972070.d0
  Mak135_V_radius_ak135( 40) =   2022400.d0
  Mak135_V_radius_ak135( 41) =   2072730.d0
  Mak135_V_radius_ak135( 42) =   2123060.d0
  Mak135_V_radius_ak135( 43) =   2173390.d0
  Mak135_V_radius_ak135( 44) =   2223720.d0
  Mak135_V_radius_ak135( 45) =   2274050.d0
  Mak135_V_radius_ak135( 46) =   2324380.d0
  Mak135_V_radius_ak135( 47) =   2374720.d0
  Mak135_V_radius_ak135( 48) =   2425050.d0
  Mak135_V_radius_ak135( 49) =   2475380.d0
  Mak135_V_radius_ak135( 50) =   2525710.d0
  Mak135_V_radius_ak135( 51) =   2576040.d0
  Mak135_V_radius_ak135( 52) =   2626370.d0
  Mak135_V_radius_ak135( 53) =   2676700.d0
  Mak135_V_radius_ak135( 54) =   2727030.d0
  Mak135_V_radius_ak135( 55) =   2777360.d0
  Mak135_V_radius_ak135( 56) =   2827700.d0
  Mak135_V_radius_ak135( 57) =   2878030.d0
  Mak135_V_radius_ak135( 58) =   2928360.d0
  Mak135_V_radius_ak135( 59) =   2978690.d0
  Mak135_V_radius_ak135( 60) =   3029020.d0
  Mak135_V_radius_ak135( 61) =   3079350.d0
  Mak135_V_radius_ak135( 62) =   3129680.d0
  Mak135_V_radius_ak135( 63) =   3180010.d0
  Mak135_V_radius_ak135( 64) =   3230340.d0
  Mak135_V_radius_ak135( 65) =   3280680.d0
  Mak135_V_radius_ak135( 66) =   3331010.d0
  Mak135_V_radius_ak135( 67) =   3381340.d0
  Mak135_V_radius_ak135( 68) =   3431670.d0
  Mak135_V_radius_ak135( 69) =   3479500.d0
  Mak135_V_radius_ak135( 70) =   3479500.d0
  Mak135_V_radius_ak135( 71) =   3531670.d0
  Mak135_V_radius_ak135( 72) =   3581330.d0
  Mak135_V_radius_ak135( 73) =   3631000.d0
  Mak135_V_radius_ak135( 74) =   3631000.d0
  Mak135_V_radius_ak135( 75) =   3681000.d0
  Mak135_V_radius_ak135( 76) =   3731000.d0
  Mak135_V_radius_ak135( 77) =   3779500.d0
  Mak135_V_radius_ak135( 78) =   3829000.d0
  Mak135_V_radius_ak135( 79) =   3878500.d0
  Mak135_V_radius_ak135( 80) =   3928000.d0
  Mak135_V_radius_ak135( 81) =   3977500.d0
  Mak135_V_radius_ak135( 82) =   4027000.d0
  Mak135_V_radius_ak135( 83) =   4076500.d0
  Mak135_V_radius_ak135( 84) =   4126000.d0
  Mak135_V_radius_ak135( 85) =   4175500.d0
  Mak135_V_radius_ak135( 86) =   4225000.d0
  Mak135_V_radius_ak135( 87) =   4274500.d0
  Mak135_V_radius_ak135( 88) =   4324000.d0
  Mak135_V_radius_ak135( 89) =   4373500.d0
  Mak135_V_radius_ak135( 90) =   4423000.d0
  Mak135_V_radius_ak135( 91) =   4472500.d0
  Mak135_V_radius_ak135( 92) =   4522000.d0
  Mak135_V_radius_ak135( 93) =   4571500.d0
  Mak135_V_radius_ak135( 94) =   4621000.d0
  Mak135_V_radius_ak135( 95) =   4670500.d0
  Mak135_V_radius_ak135( 96) =   4720000.d0
  Mak135_V_radius_ak135( 97) =   4769500.d0
  Mak135_V_radius_ak135( 98) =   4819000.d0
  Mak135_V_radius_ak135( 99) =   4868500.d0
  Mak135_V_radius_ak135(100) =   4918000.d0
  Mak135_V_radius_ak135(101) =   4967500.d0
  Mak135_V_radius_ak135(102) =   5017000.d0
  Mak135_V_radius_ak135(103) =   5066500.d0
  Mak135_V_radius_ak135(104) =   5116000.d0
  Mak135_V_radius_ak135(105) =   5165500.d0
  Mak135_V_radius_ak135(106) =   5215000.d0
  Mak135_V_radius_ak135(107) =   5264500.d0
  Mak135_V_radius_ak135(108) =   5314000.d0
  Mak135_V_radius_ak135(109) =   5363500.d0
  Mak135_V_radius_ak135(110) =   5413000.d0
  Mak135_V_radius_ak135(111) =   5462500.d0
  Mak135_V_radius_ak135(112) =   5512000.d0
  Mak135_V_radius_ak135(113) =   5561500.d0
  Mak135_V_radius_ak135(114) =   5611000.d0
  Mak135_V_radius_ak135(115) =   5661000.d0
  Mak135_V_radius_ak135(116) =   5711000.d0
  Mak135_V_radius_ak135(117) =   5711000.d0
  Mak135_V_radius_ak135(118) =   5761000.d0
  Mak135_V_radius_ak135(119) =   5811000.d0
  Mak135_V_radius_ak135(120) =   5861000.d0
  Mak135_V_radius_ak135(121) =   5911000.d0
  Mak135_V_radius_ak135(122) =   5961000.d0
  Mak135_V_radius_ak135(123) =   5961000.d0
  Mak135_V_radius_ak135(124) =   6011000.d0
  Mak135_V_radius_ak135(125) =   6061000.d0
  Mak135_V_radius_ak135(126) =   6111000.d0
  Mak135_V_radius_ak135(127) =   6161000.d0
  Mak135_V_radius_ak135(128) =   6161000.d0
  Mak135_V_radius_ak135(129) =   6206000.d0
  Mak135_V_radius_ak135(130) =   6251000.d0
  Mak135_V_radius_ak135(131) =   6293500.d0
  Mak135_V_radius_ak135(132) =   6336000.d0
  Mak135_V_radius_ak135(133) =   6336000.d0
  Mak135_V_radius_ak135(134) =   6351000.d0
  Mak135_V_radius_ak135(135) =   6351000.d0
  Mak135_V_radius_ak135(136) =   6371000.d0

  Mak135_V_density_ak135(  1) =   13.0122d0
  Mak135_V_density_ak135(  2) =   13.0117d0
  Mak135_V_density_ak135(  3) =   13.01d0
  Mak135_V_density_ak135(  4) =   13.0074d0
  Mak135_V_density_ak135(  5) =   13.0036d0
  Mak135_V_density_ak135(  6) =   12.9988d0
  Mak135_V_density_ak135(  7) =   12.9929d0
  Mak135_V_density_ak135(  8) =   12.9859d0
  Mak135_V_density_ak135(  9) =   12.9779d0
  Mak135_V_density_ak135( 10) =   12.9688d0
  Mak135_V_density_ak135( 11) =   12.9586d0
  Mak135_V_density_ak135( 12) =   12.9474d0
  Mak135_V_density_ak135( 13) =   12.9217d0
  Mak135_V_density_ak135( 14) =   12.907d0
  Mak135_V_density_ak135( 15) =   12.8917d0
  Mak135_V_density_ak135( 16) =   12.8751d0
  Mak135_V_density_ak135( 17) =   12.8574d0
  Mak135_V_density_ak135( 18) =   12.8387d0
  Mak135_V_density_ak135( 19) =   12.8188d0
  Mak135_V_density_ak135( 20) =   12.798d0
  Mak135_V_density_ak135( 21) =   12.776d0
  Mak135_V_density_ak135( 22) =   12.753d0
  Mak135_V_density_ak135( 23) =   12.7289d0
  Mak135_V_density_ak135( 24) =   12.7037d0
  Mak135_V_density_ak135( 25) =   12.1391d0
  Mak135_V_density_ak135( 26) =   12.1133d0
  Mak135_V_density_ak135( 27) =   12.0867d0
  Mak135_V_density_ak135( 28) =   12.0593d0
  Mak135_V_density_ak135( 29) =   12.0311d0
  Mak135_V_density_ak135( 30) =   12.0001d0
  Mak135_V_density_ak135( 31) =   11.9722d0
  Mak135_V_density_ak135( 32) =   11.9414d0
  Mak135_V_density_ak135( 33) =   11.8772d0
  Mak135_V_density_ak135( 34) =   11.8437d0
  Mak135_V_density_ak135( 35) =   11.8092d0
  Mak135_V_density_ak135( 36) =   11.7737d0
  Mak135_V_density_ak135( 37) =   11.7373d0
  Mak135_V_density_ak135( 38) =   11.6998d0
  Mak135_V_density_ak135( 39) =   11.6612d0
  Mak135_V_density_ak135( 40) =   11.6216d0
  Mak135_V_density_ak135( 41) =   11.5809d0
  Mak135_V_density_ak135( 42) =   11.5391d0
  Mak135_V_density_ak135( 43) =   11.4962d0
  Mak135_V_density_ak135( 44) =   11.4521d0
  Mak135_V_density_ak135( 45) =   11.4069d0
  Mak135_V_density_ak135( 46) =   11.3604d0
  Mak135_V_density_ak135( 47) =   11.3127d0
  Mak135_V_density_ak135( 48) =   11.2639d0
  Mak135_V_density_ak135( 49) =   11.2137d0
  Mak135_V_density_ak135( 50) =   11.1623d0
  Mak135_V_density_ak135( 51) =   11.1095d0
  Mak135_V_density_ak135( 52) =   11.0555d0
  Mak135_V_density_ak135( 53) =   11.0001d0
  Mak135_V_density_ak135( 54) =   10.9434d0
  Mak135_V_density_ak135( 55) =   10.8852d0
  Mak135_V_density_ak135( 56) =   10.8257d0
  Mak135_V_density_ak135( 57) =   10.7647d0
  Mak135_V_density_ak135( 58) =   10.7023d0
  Mak135_V_density_ak135( 59) =   10.6385d0
  Mak135_V_density_ak135( 60) =   10.5731d0
  Mak135_V_density_ak135( 61) =   10.5062d0
  Mak135_V_density_ak135( 62) =   10.4378d0
  Mak135_V_density_ak135( 63) =   10.3679d0
  Mak135_V_density_ak135( 64) =   10.2964d0
  Mak135_V_density_ak135( 65) =   10.2233d0
  Mak135_V_density_ak135( 66) =   10.1485d0
  Mak135_V_density_ak135( 67) =   10.0722d0
  Mak135_V_density_ak135( 68) =   9.9942d0
  Mak135_V_density_ak135( 69) =   9.9145d0
  Mak135_V_density_ak135( 70) =   5.7721d0
  Mak135_V_density_ak135( 71) =   5.7458d0
  Mak135_V_density_ak135( 72) =   5.7196d0
  Mak135_V_density_ak135( 73) =   5.6934d0
  Mak135_V_density_ak135( 74) =   5.4387d0
  Mak135_V_density_ak135( 75) =   5.4176d0
  Mak135_V_density_ak135( 76) =   5.3962d0
  Mak135_V_density_ak135( 77) =   5.3748d0
  Mak135_V_density_ak135( 78) =   5.3531d0
  Mak135_V_density_ak135( 79) =   5.3313d0
  Mak135_V_density_ak135( 80) =   5.3092d0
  Mak135_V_density_ak135( 81) =   5.287d0
  Mak135_V_density_ak135( 82) =   5.2646d0
  Mak135_V_density_ak135( 83) =   5.242d0
  Mak135_V_density_ak135( 84) =   5.2192d0
  Mak135_V_density_ak135( 85) =   5.1963d0
  Mak135_V_density_ak135( 86) =   5.1732d0
  Mak135_V_density_ak135( 87) =   5.1499d0
  Mak135_V_density_ak135( 88) =   5.1264d0
  Mak135_V_density_ak135( 89) =   5.1027d0
  Mak135_V_density_ak135( 90) =   5.0789d0
  Mak135_V_density_ak135( 91) =   5.0548d0
  Mak135_V_density_ak135( 92) =   5.0306d0
  Mak135_V_density_ak135( 93) =   5.0062d0
  Mak135_V_density_ak135( 94) =   4.9817d0
  Mak135_V_density_ak135( 95) =   4.957d0
  Mak135_V_density_ak135( 96) =   4.9321d0
  Mak135_V_density_ak135( 97) =   4.9069d0
  Mak135_V_density_ak135( 98) =   4.8817d0
  Mak135_V_density_ak135( 99) =   4.8562d0
  Mak135_V_density_ak135(100) =   4.8307d0
  Mak135_V_density_ak135(101) =   4.805d0
  Mak135_V_density_ak135(102) =   4.779d0
  Mak135_V_density_ak135(103) =   4.7528d0
  Mak135_V_density_ak135(104) =   4.7266d0
  Mak135_V_density_ak135(105) =   4.7001d0
  Mak135_V_density_ak135(106) =   4.6735d0
  Mak135_V_density_ak135(107) =   4.6467d0
  Mak135_V_density_ak135(108) =   4.6198d0
  Mak135_V_density_ak135(109) =   4.5926d0
  Mak135_V_density_ak135(110) =   4.5654d0
  Mak135_V_density_ak135(111) =   4.5162d0
  Mak135_V_density_ak135(112) =   4.465d0
  Mak135_V_density_ak135(113) =   4.4118d0
  Mak135_V_density_ak135(114) =   4.3565d0
  Mak135_V_density_ak135(115) =   4.2986d0
  Mak135_V_density_ak135(116) =   4.2387d0
  Mak135_V_density_ak135(117) =   3.9201d0
  Mak135_V_density_ak135(118) =   3.9206d0
  Mak135_V_density_ak135(119) =   3.9218d0
  Mak135_V_density_ak135(120) =   3.9233d0
  Mak135_V_density_ak135(121) =   3.9273d0
  Mak135_V_density_ak135(122) =   3.9317d0
  Mak135_V_density_ak135(123) =   3.5068d0
  Mak135_V_density_ak135(124) =   3.4577d0
  Mak135_V_density_ak135(125) =   3.411d0
  Mak135_V_density_ak135(126) =   3.3663d0
  Mak135_V_density_ak135(127) =   3.3243d0
  Mak135_V_density_ak135(128) =   3.3243d0
  Mak135_V_density_ak135(129) =   3.3711d0
  Mak135_V_density_ak135(130) =   3.4268d0
  Mak135_V_density_ak135(131) =   3.345d0
  Mak135_V_density_ak135(132) =   3.32d0
  Mak135_V_density_ak135(133) =   2.92d0
  Mak135_V_density_ak135(134) =   2.92d0
  Mak135_V_density_ak135(135) =   2.72d0
  Mak135_V_density_ak135(136) =   2.72d0

  Mak135_V_vp_ak135(  1) =   11.2622d0
  Mak135_V_vp_ak135(  2) =   11.2618d0
  Mak135_V_vp_ak135(  3) =   11.2606d0
  Mak135_V_vp_ak135(  4) =   11.2586d0
  Mak135_V_vp_ak135(  5) =   11.2557d0
  Mak135_V_vp_ak135(  6) =   11.2521d0
  Mak135_V_vp_ak135(  7) =   11.2477d0
  Mak135_V_vp_ak135(  8) =   11.2424d0
  Mak135_V_vp_ak135(  9) =   11.2364d0
  Mak135_V_vp_ak135( 10) =   11.2295d0
  Mak135_V_vp_ak135( 11) =   11.2219d0
  Mak135_V_vp_ak135( 12) =   11.2134d0
  Mak135_V_vp_ak135( 13) =   11.1941d0
  Mak135_V_vp_ak135( 14) =   11.183d0
  Mak135_V_vp_ak135( 15) =   11.1715d0
  Mak135_V_vp_ak135( 16) =   11.159d0
  Mak135_V_vp_ak135( 17) =   11.1457d0
  Mak135_V_vp_ak135( 18) =   11.1316d0
  Mak135_V_vp_ak135( 19) =   11.1166d0
  Mak135_V_vp_ak135( 20) =   11.0983d0
  Mak135_V_vp_ak135( 21) =   11.085d0
  Mak135_V_vp_ak135( 22) =   11.0718d0
  Mak135_V_vp_ak135( 23) =   11.0585d0
  Mak135_V_vp_ak135( 24) =   11.0427d0
  Mak135_V_vp_ak135( 25) =   10.289d0
  Mak135_V_vp_ak135( 26) =   10.2854d0
  Mak135_V_vp_ak135( 27) =   10.2745d0
  Mak135_V_vp_ak135( 28) =   10.2565d0
  Mak135_V_vp_ak135( 29) =   10.2329d0
  Mak135_V_vp_ak135( 30) =   10.2049d0
  Mak135_V_vp_ak135( 31) =   10.1739d0
  Mak135_V_vp_ak135( 32) =   10.1415d0
  Mak135_V_vp_ak135( 33) =   10.0768d0
  Mak135_V_vp_ak135( 34) =   10.0439d0
  Mak135_V_vp_ak135( 35) =   10.0103d0
  Mak135_V_vp_ak135( 36) =   9.9761d0
  Mak135_V_vp_ak135( 37) =   9.941d0
  Mak135_V_vp_ak135( 38) =   9.9051d0
  Mak135_V_vp_ak135( 39) =   9.8682d0
  Mak135_V_vp_ak135( 40) =   9.8304d0
  Mak135_V_vp_ak135( 41) =   9.7914d0
  Mak135_V_vp_ak135( 42) =   9.7513d0
  Mak135_V_vp_ak135( 43) =   9.71d0
  Mak135_V_vp_ak135( 44) =   9.6673d0
  Mak135_V_vp_ak135( 45) =   9.6232d0
  Mak135_V_vp_ak135( 46) =   9.5777d0
  Mak135_V_vp_ak135( 47) =   9.5306d0
  Mak135_V_vp_ak135( 48) =   9.4814d0
  Mak135_V_vp_ak135( 49) =   9.4297d0
  Mak135_V_vp_ak135( 50) =   9.376d0
  Mak135_V_vp_ak135( 51) =   9.3205d0
  Mak135_V_vp_ak135( 52) =   9.2634d0
  Mak135_V_vp_ak135( 53) =   9.2042d0
  Mak135_V_vp_ak135( 54) =   9.1426d0
  Mak135_V_vp_ak135( 55) =   9.0792d0
  Mak135_V_vp_ak135( 56) =   9.0138d0
  Mak135_V_vp_ak135( 57) =   8.9461d0
  Mak135_V_vp_ak135( 58) =   8.8761d0
  Mak135_V_vp_ak135( 59) =   8.8036d0
  Mak135_V_vp_ak135( 60) =   8.7283d0
  Mak135_V_vp_ak135( 61) =   8.6496d0
  Mak135_V_vp_ak135( 62) =   8.5692d0
  Mak135_V_vp_ak135( 63) =   8.4861d0
  Mak135_V_vp_ak135( 64) =   8.4001d0
  Mak135_V_vp_ak135( 65) =   8.3122d0
  Mak135_V_vp_ak135( 66) =   8.2213d0
  Mak135_V_vp_ak135( 67) =   8.1283d0
  Mak135_V_vp_ak135( 68) =   8.0382d0
  Mak135_V_vp_ak135( 69) =   8.d0
  Mak135_V_vp_ak135( 70) =   13.6601d0
  Mak135_V_vp_ak135( 71) =   13.657d0
  Mak135_V_vp_ak135( 72) =   13.6533d0
  Mak135_V_vp_ak135( 73) =   13.6498d0
  Mak135_V_vp_ak135( 74) =   13.6498d0
  Mak135_V_vp_ak135( 75) =   13.5899d0
  Mak135_V_vp_ak135( 76) =   13.5311d0
  Mak135_V_vp_ak135( 77) =   13.4741d0
  Mak135_V_vp_ak135( 78) =   13.4156d0
  Mak135_V_vp_ak135( 79) =   13.3584d0
  Mak135_V_vp_ak135( 80) =   13.3017d0
  Mak135_V_vp_ak135( 81) =   13.2465d0
  Mak135_V_vp_ak135( 82) =   13.1895d0
  Mak135_V_vp_ak135( 83) =   13.1337d0
  Mak135_V_vp_ak135( 84) =   13.0786d0
  Mak135_V_vp_ak135( 85) =   13.0226d0
  Mak135_V_vp_ak135( 86) =   12.9663d0
  Mak135_V_vp_ak135( 87) =   12.9093d0
  Mak135_V_vp_ak135( 88) =   12.8524d0
  Mak135_V_vp_ak135( 89) =   12.7956d0
  Mak135_V_vp_ak135( 90) =   12.7384d0
  Mak135_V_vp_ak135( 91) =   12.6807d0
  Mak135_V_vp_ak135( 92) =   12.6226d0
  Mak135_V_vp_ak135( 93) =   12.5638d0
  Mak135_V_vp_ak135( 94) =   12.503d0
  Mak135_V_vp_ak135( 95) =   12.4427d0
  Mak135_V_vp_ak135( 96) =   12.3813d0
  Mak135_V_vp_ak135( 97) =   12.3181d0
  Mak135_V_vp_ak135( 98) =   12.2558d0
  Mak135_V_vp_ak135( 99) =   12.1912d0
  Mak135_V_vp_ak135(100) =   12.1247d0
  Mak135_V_vp_ak135(101) =   12.0571d0
  Mak135_V_vp_ak135(102) =   11.9891d0
  Mak135_V_vp_ak135(103) =   11.9208d0
  Mak135_V_vp_ak135(104) =   11.8491d0
  Mak135_V_vp_ak135(105) =   11.7768d0
  Mak135_V_vp_ak135(106) =   11.702d0
  Mak135_V_vp_ak135(107) =   11.6265d0
  Mak135_V_vp_ak135(108) =   11.5493d0
  Mak135_V_vp_ak135(109) =   11.4704d0
  Mak135_V_vp_ak135(110) =   11.3897d0
  Mak135_V_vp_ak135(111) =   11.3068d0
  Mak135_V_vp_ak135(112) =   11.2228d0
  Mak135_V_vp_ak135(113) =   11.1355d0
  Mak135_V_vp_ak135(114) =   11.0553d0
  Mak135_V_vp_ak135(115) =   10.9222d0
  Mak135_V_vp_ak135(116) =   10.7909d0
  Mak135_V_vp_ak135(117) =   10.2d0
  Mak135_V_vp_ak135(118) =   10.032d0
  Mak135_V_vp_ak135(119) =   9.864d0
  Mak135_V_vp_ak135(120) =   9.6962d0
  Mak135_V_vp_ak135(121) =   9.528d0
  Mak135_V_vp_ak135(122) =   9.3601d0
  Mak135_V_vp_ak135(123) =   9.0302d0
  Mak135_V_vp_ak135(124) =   8.8476d0
  Mak135_V_vp_ak135(125) =   8.665d0
  Mak135_V_vp_ak135(126) =   8.4822d0
  Mak135_V_vp_ak135(127) =   8.3007d0
  Mak135_V_vp_ak135(128) =   8.3007d0
  Mak135_V_vp_ak135(129) =   8.175d0
  Mak135_V_vp_ak135(130) =   8.0505d0
  Mak135_V_vp_ak135(131) =   8.045d0
  Mak135_V_vp_ak135(132) =   8.04d0
  Mak135_V_vp_ak135(133) =   6.5d0
  Mak135_V_vp_ak135(134) =   6.5d0
  Mak135_V_vp_ak135(135) =   5.8d0
  Mak135_V_vp_ak135(136) =   5.8d0

  Mak135_V_vs_ak135(  1) =   3.6678d0
  Mak135_V_vs_ak135(  2) =   3.6675d0
  Mak135_V_vs_ak135(  3) =   3.6667d0
  Mak135_V_vs_ak135(  4) =   3.6653d0
  Mak135_V_vs_ak135(  5) =   3.6633d0
  Mak135_V_vs_ak135(  6) =   3.6608d0
  Mak135_V_vs_ak135(  7) =   3.6577d0
  Mak135_V_vs_ak135(  8) =   3.654d0
  Mak135_V_vs_ak135(  9) =   3.6498d0
  Mak135_V_vs_ak135( 10) =   3.645d0
  Mak135_V_vs_ak135( 11) =   3.6396d0
  Mak135_V_vs_ak135( 12) =   3.6337d0
  Mak135_V_vs_ak135( 13) =   3.6202d0
  Mak135_V_vs_ak135( 14) =   3.613d0
  Mak135_V_vs_ak135( 15) =   3.6044d0
  Mak135_V_vs_ak135( 16) =   3.5957d0
  Mak135_V_vs_ak135( 17) =   3.5864d0
  Mak135_V_vs_ak135( 18) =   3.5765d0
  Mak135_V_vs_ak135( 19) =   3.5661d0
  Mak135_V_vs_ak135( 20) =   3.5551d0
  Mak135_V_vs_ak135( 21) =   3.5435d0
  Mak135_V_vs_ak135( 22) =   3.5314d0
  Mak135_V_vs_ak135( 23) =   3.5187d0
  Mak135_V_vs_ak135( 24) =   3.5043d0
  Mak135_V_vs_ak135( 25) =  0.d0
  Mak135_V_vs_ak135( 26) =  0.d0
  Mak135_V_vs_ak135( 27) =  0.d0
  Mak135_V_vs_ak135( 28) =  0.d0
  Mak135_V_vs_ak135( 29) =  0.d0
  Mak135_V_vs_ak135( 30) =  0.d0
  Mak135_V_vs_ak135( 31) =  0.d0
  Mak135_V_vs_ak135( 32) =  0.d0
  Mak135_V_vs_ak135( 33) =  0.d0
  Mak135_V_vs_ak135( 34) =  0.d0
  Mak135_V_vs_ak135( 35) =  0.d0
  Mak135_V_vs_ak135( 36) =  0.d0
  Mak135_V_vs_ak135( 37) =  0.d0
  Mak135_V_vs_ak135( 38) =  0.d0
  Mak135_V_vs_ak135( 39) =  0.d0
  Mak135_V_vs_ak135( 40) =  0.d0
  Mak135_V_vs_ak135( 41) =  0.d0
  Mak135_V_vs_ak135( 42) =  0.d0
  Mak135_V_vs_ak135( 43) =  0.d0
  Mak135_V_vs_ak135( 44) =  0.d0
  Mak135_V_vs_ak135( 45) =  0.d0
  Mak135_V_vs_ak135( 46) =  0.d0
  Mak135_V_vs_ak135( 47) =  0.d0
  Mak135_V_vs_ak135( 48) =  0.d0
  Mak135_V_vs_ak135( 49) =  0.d0
  Mak135_V_vs_ak135( 50) =  0.d0
  Mak135_V_vs_ak135( 51) =  0.d0
  Mak135_V_vs_ak135( 52) =  0.d0
  Mak135_V_vs_ak135( 53) =  0.d0
  Mak135_V_vs_ak135( 54) =  0.d0
  Mak135_V_vs_ak135( 55) =  0.d0
  Mak135_V_vs_ak135( 56) =  0.d0
  Mak135_V_vs_ak135( 57) =  0.d0
  Mak135_V_vs_ak135( 58) =  0.d0
  Mak135_V_vs_ak135( 59) =  0.d0
  Mak135_V_vs_ak135( 60) =  0.d0
  Mak135_V_vs_ak135( 61) =  0.d0
  Mak135_V_vs_ak135( 62) =  0.d0
  Mak135_V_vs_ak135( 63) =  0.d0
  Mak135_V_vs_ak135( 64) =  0.d0
  Mak135_V_vs_ak135( 65) =  0.d0
  Mak135_V_vs_ak135( 66) =  0.d0
  Mak135_V_vs_ak135( 67) =  0.d0
  Mak135_V_vs_ak135( 68) =  0.d0
  Mak135_V_vs_ak135( 69) =  0.d0
  Mak135_V_vs_ak135( 70) =   7.2817d0
  Mak135_V_vs_ak135( 71) =   7.27d0
  Mak135_V_vs_ak135( 72) =   7.2593d0
  Mak135_V_vs_ak135( 73) =   7.2485d0
  Mak135_V_vs_ak135( 74) =   7.2485d0
  Mak135_V_vs_ak135( 75) =   7.2253d0
  Mak135_V_vs_ak135( 76) =   7.2031d0
  Mak135_V_vs_ak135( 77) =   7.1804d0
  Mak135_V_vs_ak135( 78) =   7.1584d0
  Mak135_V_vs_ak135( 79) =   7.1368d0
  Mak135_V_vs_ak135( 80) =   7.1144d0
  Mak135_V_vs_ak135( 81) =   7.0932d0
  Mak135_V_vs_ak135( 82) =   7.0722d0
  Mak135_V_vs_ak135( 83) =   7.0504d0
  Mak135_V_vs_ak135( 84) =   7.0286d0
  Mak135_V_vs_ak135( 85) =   7.0069d0
  Mak135_V_vs_ak135( 86) =   6.9852d0
  Mak135_V_vs_ak135( 87) =   6.9625d0
  Mak135_V_vs_ak135( 88) =   6.9416d0
  Mak135_V_vs_ak135( 89) =   6.9194d0
  Mak135_V_vs_ak135( 90) =   6.8972d0
  Mak135_V_vs_ak135( 91) =   6.8743d0
  Mak135_V_vs_ak135( 92) =   6.8517d0
  Mak135_V_vs_ak135( 93) =   6.8289d0
  Mak135_V_vs_ak135( 94) =   6.8056d0
  Mak135_V_vs_ak135( 95) =   6.782d0
  Mak135_V_vs_ak135( 96) =   6.7579d0
  Mak135_V_vs_ak135( 97) =   6.7323d0
  Mak135_V_vs_ak135( 98) =   6.707d0
  Mak135_V_vs_ak135( 99) =   6.6813d0
  Mak135_V_vs_ak135(100) =   6.6554d0
  Mak135_V_vs_ak135(101) =   6.6285d0
  Mak135_V_vs_ak135(102) =   6.6009d0
  Mak135_V_vs_ak135(103) =   6.5728d0
  Mak135_V_vs_ak135(104) =   6.5431d0
  Mak135_V_vs_ak135(105) =   6.5131d0
  Mak135_V_vs_ak135(106) =   6.4822d0
  Mak135_V_vs_ak135(107) =   6.4514d0
  Mak135_V_vs_ak135(108) =   6.4182d0
  Mak135_V_vs_ak135(109) =   6.386d0
  Mak135_V_vs_ak135(110) =   6.3519d0
  Mak135_V_vs_ak135(111) =   6.3164d0
  Mak135_V_vs_ak135(112) =   6.2799d0
  Mak135_V_vs_ak135(113) =   6.2424d0
  Mak135_V_vs_ak135(114) =   6.21d0
  Mak135_V_vs_ak135(115) =   6.0898d0
  Mak135_V_vs_ak135(116) =   5.9607d0
  Mak135_V_vs_ak135(117) =   5.6104d0
  Mak135_V_vs_ak135(118) =   5.5047d0
  Mak135_V_vs_ak135(119) =   5.3989d0
  Mak135_V_vs_ak135(120) =   5.2922d0
  Mak135_V_vs_ak135(121) =   5.1864d0
  Mak135_V_vs_ak135(122) =   5.0806d0
  Mak135_V_vs_ak135(123) =   4.8702d0
  Mak135_V_vs_ak135(124) =   4.7832d0
  Mak135_V_vs_ak135(125) =   4.6964d0
  Mak135_V_vs_ak135(126) =   4.6094d0
  Mak135_V_vs_ak135(127) =   4.5184d0
  Mak135_V_vs_ak135(128) =   4.5184d0
  Mak135_V_vs_ak135(129) =   4.509d0
  Mak135_V_vs_ak135(130) =   4.5d0
  Mak135_V_vs_ak135(131) =   4.49d0
  Mak135_V_vs_ak135(132) =   4.48d0
  Mak135_V_vs_ak135(133) =   3.85d0
  Mak135_V_vs_ak135(134) =   3.85d0
  Mak135_V_vs_ak135(135) =   3.46d0
  Mak135_V_vs_ak135(136) =   3.46d0

  Mak135_V_Qkappa_ak135(  1) =   601.27d0
  Mak135_V_Qkappa_ak135(  2) =   601.32d0
  Mak135_V_Qkappa_ak135(  3) =   601.46d0
  Mak135_V_Qkappa_ak135(  4) =   601.7d0
  Mak135_V_Qkappa_ak135(  5) =   602.05d0
  Mak135_V_Qkappa_ak135(  6) =   602.49d0
  Mak135_V_Qkappa_ak135(  7) =   603.04d0
  Mak135_V_Qkappa_ak135(  8) =   603.69d0
  Mak135_V_Qkappa_ak135(  9) =   604.44d0
  Mak135_V_Qkappa_ak135( 10) =   605.28d0
  Mak135_V_Qkappa_ak135( 11) =   606.26d0
  Mak135_V_Qkappa_ak135( 12) =   607.31d0
  Mak135_V_Qkappa_ak135( 13) =   609.74d0
  Mak135_V_Qkappa_ak135( 14) =   611.18d0
  Mak135_V_Qkappa_ak135( 15) =   612.62d0
  Mak135_V_Qkappa_ak135( 16) =   614.21d0
  Mak135_V_Qkappa_ak135( 17) =   615.93d0
  Mak135_V_Qkappa_ak135( 18) =   617.78d0
  Mak135_V_Qkappa_ak135( 19) =   619.71d0
  Mak135_V_Qkappa_ak135( 20) =   621.5d0
  Mak135_V_Qkappa_ak135( 21) =   624.08d0
  Mak135_V_Qkappa_ak135( 22) =   626.87d0
  Mak135_V_Qkappa_ak135( 23) =   629.89d0
  Mak135_V_Qkappa_ak135( 24) =   633.26d0
  Mak135_V_Qkappa_ak135( 25) =   57822.d0
  Mak135_V_Qkappa_ak135( 26) =   57822.d0
  Mak135_V_Qkappa_ak135( 27) =   57822.d0
  Mak135_V_Qkappa_ak135( 28) =   57822.d0
  Mak135_V_Qkappa_ak135( 29) =   57822.d0
  Mak135_V_Qkappa_ak135( 30) =   57822.d0
  Mak135_V_Qkappa_ak135( 31) =   57822.d0
  Mak135_V_Qkappa_ak135( 32) =   57822.d0
  Mak135_V_Qkappa_ak135( 33) =   57822.d0
  Mak135_V_Qkappa_ak135( 34) =   57822.d0
  Mak135_V_Qkappa_ak135( 35) =   57822.d0
  Mak135_V_Qkappa_ak135( 36) =   57822.d0
  Mak135_V_Qkappa_ak135( 37) =   57822.d0
  Mak135_V_Qkappa_ak135( 38) =   57822.d0
  Mak135_V_Qkappa_ak135( 39) =   57822.d0
  Mak135_V_Qkappa_ak135( 40) =   57822.d0
  Mak135_V_Qkappa_ak135( 41) =   57822.d0
  Mak135_V_Qkappa_ak135( 42) =   57822.d0
  Mak135_V_Qkappa_ak135( 43) =   57822.d0
  Mak135_V_Qkappa_ak135( 44) =   57822.d0
  Mak135_V_Qkappa_ak135( 45) =   57822.d0
  Mak135_V_Qkappa_ak135( 46) =   57822.d0
  Mak135_V_Qkappa_ak135( 47) =   57822.d0
  Mak135_V_Qkappa_ak135( 48) =   57822.d0
  Mak135_V_Qkappa_ak135( 49) =   57822.d0
  Mak135_V_Qkappa_ak135( 50) =   57822.d0
  Mak135_V_Qkappa_ak135( 51) =   57822.d0
  Mak135_V_Qkappa_ak135( 52) =   57822.d0
  Mak135_V_Qkappa_ak135( 53) =   57822.d0
  Mak135_V_Qkappa_ak135( 54) =   57822.d0
  Mak135_V_Qkappa_ak135( 55) =   57822.d0
  Mak135_V_Qkappa_ak135( 56) =   57822.d0
  Mak135_V_Qkappa_ak135( 57) =   57822.d0
  Mak135_V_Qkappa_ak135( 58) =   57822.d0
  Mak135_V_Qkappa_ak135( 59) =   57822.d0
  Mak135_V_Qkappa_ak135( 60) =   57822.d0
  Mak135_V_Qkappa_ak135( 61) =   57822.d0
  Mak135_V_Qkappa_ak135( 62) =   57822.d0
  Mak135_V_Qkappa_ak135( 63) =   57822.d0
  Mak135_V_Qkappa_ak135( 64) =   57822.d0
  Mak135_V_Qkappa_ak135( 65) =   57822.d0
  Mak135_V_Qkappa_ak135( 66) =   57822.d0
  Mak135_V_Qkappa_ak135( 67) =   57822.d0
  Mak135_V_Qkappa_ak135( 68) =   57822.d0
  Mak135_V_Qkappa_ak135( 69) =   57822.d0
  Mak135_V_Qkappa_ak135( 70) =   723.12d0
  Mak135_V_Qkappa_ak135( 71) =   725.11d0
  Mak135_V_Qkappa_ak135( 72) =   726.87d0
  Mak135_V_Qkappa_ak135( 73) =   722.73d0
  Mak135_V_Qkappa_ak135( 74) =   933.21d0
  Mak135_V_Qkappa_ak135( 75) =   940.88d0
  Mak135_V_Qkappa_ak135( 76) =   952.d0
  Mak135_V_Qkappa_ak135( 77) =   960.36d0
  Mak135_V_Qkappa_ak135( 78) =   968.46d0
  Mak135_V_Qkappa_ak135( 79) =   976.81d0
  Mak135_V_Qkappa_ak135( 80) =   985.63d0
  Mak135_V_Qkappa_ak135( 81) =   990.77d0
  Mak135_V_Qkappa_ak135( 82) =   999.44d0
  Mak135_V_Qkappa_ak135( 83) =   1008.79d0
  Mak135_V_Qkappa_ak135( 84) =   1018.38d0
  Mak135_V_Qkappa_ak135( 85) =   1032.14d0
  Mak135_V_Qkappa_ak135( 86) =   1042.07d0
  Mak135_V_Qkappa_ak135( 87) =   1048.09d0
  Mak135_V_Qkappa_ak135( 88) =   1058.03d0
  Mak135_V_Qkappa_ak135( 89) =   1064.23d0
  Mak135_V_Qkappa_ak135( 90) =   1070.38d0
  Mak135_V_Qkappa_ak135( 91) =   1085.97d0
  Mak135_V_Qkappa_ak135( 92) =   1097.16d0
  Mak135_V_Qkappa_ak135( 93) =   1108.58d0
  Mak135_V_Qkappa_ak135( 94) =   1120.09d0
  Mak135_V_Qkappa_ak135( 95) =   1127.02d0
  Mak135_V_Qkappa_ak135( 96) =   1134.01d0
  Mak135_V_Qkappa_ak135( 97) =   1141.32d0
  Mak135_V_Qkappa_ak135( 98) =   1148.76d0
  Mak135_V_Qkappa_ak135( 99) =   1156.04d0
  Mak135_V_Qkappa_ak135(100) =   1163.16d0
  Mak135_V_Qkappa_ak135(101) =   1170.53d0
  Mak135_V_Qkappa_ak135(102) =   1178.19d0
  Mak135_V_Qkappa_ak135(103) =   1186.06d0
  Mak135_V_Qkappa_ak135(104) =   1193.99d0
  Mak135_V_Qkappa_ak135(105) =   1202.04d0
  Mak135_V_Qkappa_ak135(106) =   1210.02d0
  Mak135_V_Qkappa_ak135(107) =   1217.91d0
  Mak135_V_Qkappa_ak135(108) =   1226.52d0
  Mak135_V_Qkappa_ak135(109) =   1234.54d0
  Mak135_V_Qkappa_ak135(110) =   1243.02d0
  Mak135_V_Qkappa_ak135(111) =   1251.69d0
  Mak135_V_Qkappa_ak135(112) =   1260.68d0
  Mak135_V_Qkappa_ak135(113) =   1269.44d0
  Mak135_V_Qkappa_ak135(114) =   1277.93d0
  Mak135_V_Qkappa_ak135(115) =   1311.17d0
  Mak135_V_Qkappa_ak135(116) =   1350.54d0
  Mak135_V_Qkappa_ak135(117) =   428.69d0
  Mak135_V_Qkappa_ak135(118) =   425.51d0
  Mak135_V_Qkappa_ak135(119) =   422.55d0
  Mak135_V_Qkappa_ak135(120) =   419.94d0
  Mak135_V_Qkappa_ak135(121) =   417.32d0
  Mak135_V_Qkappa_ak135(122) =   413.66d0
  Mak135_V_Qkappa_ak135(123) =   377.93d0
  Mak135_V_Qkappa_ak135(124) =   366.34d0
  Mak135_V_Qkappa_ak135(125) =   355.85d0
  Mak135_V_Qkappa_ak135(126) =   346.37d0
  Mak135_V_Qkappa_ak135(127) =   338.47d0
  Mak135_V_Qkappa_ak135(128) =   200.97d0
  Mak135_V_Qkappa_ak135(129) =   188.72d0
  Mak135_V_Qkappa_ak135(130) =   182.57d0
  Mak135_V_Qkappa_ak135(131) =   182.03d0
  Mak135_V_Qkappa_ak135(132) =   182.03d0
  Mak135_V_Qkappa_ak135(133) =   972.77d0
  Mak135_V_Qkappa_ak135(134) =   972.77d0
  Mak135_V_Qkappa_ak135(135) =   1368.02d0
  Mak135_V_Qkappa_ak135(136) =   1368.02d0

  Mak135_V_Qmu_ak135(  1) =   85.03d0
  Mak135_V_Qmu_ak135(  2) =   85.03d0
  Mak135_V_Qmu_ak135(  3) =   85.03d0
  Mak135_V_Qmu_ak135(  4) =   85.03d0
  Mak135_V_Qmu_ak135(  5) =   85.03d0
  Mak135_V_Qmu_ak135(  6) =   85.03d0
  Mak135_V_Qmu_ak135(  7) =   85.03d0
  Mak135_V_Qmu_ak135(  8) =   85.03d0
  Mak135_V_Qmu_ak135(  9) =   85.03d0
  Mak135_V_Qmu_ak135( 10) =   85.03d0
  Mak135_V_Qmu_ak135( 11) =   85.03d0
  Mak135_V_Qmu_ak135( 12) =   85.03d0
  Mak135_V_Qmu_ak135( 13) =   85.03d0
  Mak135_V_Qmu_ak135( 14) =   85.03d0
  Mak135_V_Qmu_ak135( 15) =   85.03d0
  Mak135_V_Qmu_ak135( 16) =   85.03d0
  Mak135_V_Qmu_ak135( 17) =   85.03d0
  Mak135_V_Qmu_ak135( 18) =   85.03d0
  Mak135_V_Qmu_ak135( 19) =   85.03d0
  Mak135_V_Qmu_ak135( 20) =   85.03d0
  Mak135_V_Qmu_ak135( 21) =   85.03d0
  Mak135_V_Qmu_ak135( 22) =   85.03d0
  Mak135_V_Qmu_ak135( 23) =   85.03d0
  Mak135_V_Qmu_ak135( 24) =   85.03d0
  Mak135_V_Qmu_ak135( 25) =  0.d0
  Mak135_V_Qmu_ak135( 26) =  0.d0
  Mak135_V_Qmu_ak135( 27) =  0.d0
  Mak135_V_Qmu_ak135( 28) =  0.d0
  Mak135_V_Qmu_ak135( 29) =  0.d0
  Mak135_V_Qmu_ak135( 30) =  0.d0
  Mak135_V_Qmu_ak135( 31) =  0.d0
  Mak135_V_Qmu_ak135( 32) =  0.d0
  Mak135_V_Qmu_ak135( 33) =  0.d0
  Mak135_V_Qmu_ak135( 34) =  0.d0
  Mak135_V_Qmu_ak135( 35) =  0.d0
  Mak135_V_Qmu_ak135( 36) =  0.d0
  Mak135_V_Qmu_ak135( 37) =  0.d0
  Mak135_V_Qmu_ak135( 38) =  0.d0
  Mak135_V_Qmu_ak135( 39) =  0.d0
  Mak135_V_Qmu_ak135( 40) =  0.d0
  Mak135_V_Qmu_ak135( 41) =  0.d0
  Mak135_V_Qmu_ak135( 42) =  0.d0
  Mak135_V_Qmu_ak135( 43) =  0.d0
  Mak135_V_Qmu_ak135( 44) =  0.d0
  Mak135_V_Qmu_ak135( 45) =  0.d0
  Mak135_V_Qmu_ak135( 46) =  0.d0
  Mak135_V_Qmu_ak135( 47) =  0.d0
  Mak135_V_Qmu_ak135( 48) =  0.d0
  Mak135_V_Qmu_ak135( 49) =  0.d0
  Mak135_V_Qmu_ak135( 50) =  0.d0
  Mak135_V_Qmu_ak135( 51) =  0.d0
  Mak135_V_Qmu_ak135( 52) =  0.d0
  Mak135_V_Qmu_ak135( 53) =  0.d0
  Mak135_V_Qmu_ak135( 54) =  0.d0
  Mak135_V_Qmu_ak135( 55) =  0.d0
  Mak135_V_Qmu_ak135( 56) =  0.d0
  Mak135_V_Qmu_ak135( 57) =  0.d0
  Mak135_V_Qmu_ak135( 58) =  0.d0
  Mak135_V_Qmu_ak135( 59) =  0.d0
  Mak135_V_Qmu_ak135( 60) =  0.d0
  Mak135_V_Qmu_ak135( 61) =  0.d0
  Mak135_V_Qmu_ak135( 62) =  0.d0
  Mak135_V_Qmu_ak135( 63) =  0.d0
  Mak135_V_Qmu_ak135( 64) =  0.d0
  Mak135_V_Qmu_ak135( 65) =  0.d0
  Mak135_V_Qmu_ak135( 66) =  0.d0
  Mak135_V_Qmu_ak135( 67) =  0.d0
  Mak135_V_Qmu_ak135( 68) =  0.d0
  Mak135_V_Qmu_ak135( 69) =  0.d0
  Mak135_V_Qmu_ak135( 70) =   273.97d0
  Mak135_V_Qmu_ak135( 71) =   273.97d0
  Mak135_V_Qmu_ak135( 72) =   273.97d0
  Mak135_V_Qmu_ak135( 73) =   271.74d0
  Mak135_V_Qmu_ak135( 74) =   350.88d0
  Mak135_V_Qmu_ak135( 75) =   354.61d0
  Mak135_V_Qmu_ak135( 76) =   359.71d0
  Mak135_V_Qmu_ak135( 77) =   363.64d0
  Mak135_V_Qmu_ak135( 78) =   367.65d0
  Mak135_V_Qmu_ak135( 79) =   371.75d0
  Mak135_V_Qmu_ak135( 80) =   375.94d0
  Mak135_V_Qmu_ak135( 81) =   378.79d0
  Mak135_V_Qmu_ak135( 82) =   383.14d0
  Mak135_V_Qmu_ak135( 83) =   387.6d0
  Mak135_V_Qmu_ak135( 84) =   392.16d0
  Mak135_V_Qmu_ak135( 85) =   398.41d0
  Mak135_V_Qmu_ak135( 86) =   403.23d0
  Mak135_V_Qmu_ak135( 87) =   406.5d0
  Mak135_V_Qmu_ak135( 88) =   411.52d0
  Mak135_V_Qmu_ak135( 89) =   414.94d0
  Mak135_V_Qmu_ak135( 90) =   418.41d0
  Mak135_V_Qmu_ak135( 91) =   425.53d0
  Mak135_V_Qmu_ak135( 92) =   431.03d0
  Mak135_V_Qmu_ak135( 93) =   436.68d0
  Mak135_V_Qmu_ak135( 94) =   442.48d0
  Mak135_V_Qmu_ak135( 95) =   446.43d0
  Mak135_V_Qmu_ak135( 96) =   450.45d0
  Mak135_V_Qmu_ak135( 97) =   454.55d0
  Mak135_V_Qmu_ak135( 98) =   458.72d0
  Mak135_V_Qmu_ak135( 99) =   462.96d0
  Mak135_V_Qmu_ak135(100) =   467.29d0
  Mak135_V_Qmu_ak135(101) =   471.7d0
  Mak135_V_Qmu_ak135(102) =   476.19d0
  Mak135_V_Qmu_ak135(103) =   480.77d0
  Mak135_V_Qmu_ak135(104) =   485.44d0
  Mak135_V_Qmu_ak135(105) =   490.2d0
  Mak135_V_Qmu_ak135(106) =   495.05d0
  Mak135_V_Qmu_ak135(107) =   500.d0
  Mak135_V_Qmu_ak135(108) =   505.05d0
  Mak135_V_Qmu_ak135(109) =   510.2d0
  Mak135_V_Qmu_ak135(110) =   515.46d0
  Mak135_V_Qmu_ak135(111) =   520.83d0
  Mak135_V_Qmu_ak135(112) =   526.32d0
  Mak135_V_Qmu_ak135(113) =   531.91d0
  Mak135_V_Qmu_ak135(114) =   537.63d0
  Mak135_V_Qmu_ak135(115) =   543.48d0
  Mak135_V_Qmu_ak135(116) =   549.45d0
  Mak135_V_Qmu_ak135(117) =   172.93d0
  Mak135_V_Qmu_ak135(118) =   170.82d0
  Mak135_V_Qmu_ak135(119) =   168.78d0
  Mak135_V_Qmu_ak135(120) =   166.8d0
  Mak135_V_Qmu_ak135(121) =   164.87d0
  Mak135_V_Qmu_ak135(122) =   162.5d0
  Mak135_V_Qmu_ak135(123) =   146.57d0
  Mak135_V_Qmu_ak135(124) =   142.76d0
  Mak135_V_Qmu_ak135(125) =   139.38d0
  Mak135_V_Qmu_ak135(126) =   136.38d0
  Mak135_V_Qmu_ak135(127) =   133.72d0
  Mak135_V_Qmu_ak135(128) =   79.4d0
  Mak135_V_Qmu_ak135(129) =   76.55d0
  Mak135_V_Qmu_ak135(130) =   76.06d0
  Mak135_V_Qmu_ak135(131) =   75.6d0
  Mak135_V_Qmu_ak135(132) =   75.6d0
  Mak135_V_Qmu_ak135(133) =   403.93d0
  Mak135_V_Qmu_ak135(134) =   403.93d0
  Mak135_V_Qmu_ak135(135) =   599.99d0
  Mak135_V_Qmu_ak135(136) =   599.99d0

  ! strip the crust and replace it with mantle
  if (SUPPRESS_CRUSTAL_MESH .or. USE_EXTERNAL_CRUSTAL_MODEL) then
    Mak135_V_vp_ak135(133:136) = Mak135_V_vp_ak135(132)
    Mak135_V_vs_ak135(133:136) = Mak135_V_vs_ak135(132)
    Mak135_V_density_ak135(133:136) = Mak135_V_density_ak135(132)
    Mak135_V_Qkappa_ak135(133:136) = Mak135_V_Qkappa_ak135(132)
    Mak135_V_Qmu_ak135(133:136) = Mak135_V_Qmu_ak135(132)
  endif

  end subroutine define_model_ak135

