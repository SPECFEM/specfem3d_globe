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

!--------------------------------------------------------------------------------------------------
! 1066A
!
! Spherically symmetric earth model 1066A [Gilbert and Dziewonski, 1975].
!
! When ATTENUATION is on, it uses an unpublished 1D attenuation model from Scripps.
!--------------------------------------------------------------------------------------------------

  module model_1066a_par

  ! number of layers in DATA/1066a/1066a.dat
  integer, parameter :: NR_1066A = 160

  ! model_1066a_variables
  double precision, dimension(:),allocatable :: &
    M1066a_V_radius_1066a,M1066a_V_density_1066a, &
    M1066a_V_vp_1066a,M1066a_V_vs_1066a, &
    M1066a_V_Qkappa_1066a,M1066a_V_Qmu_1066a

  end module model_1066a_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_1066a_broadcast(CRUSTAL)

! standard routine to setup model

  use constants, only: myrank
  use model_1066a_par

  implicit none

  logical :: CRUSTAL

  ! local parameters
  integer :: ier

  ! allocates model arrays
  allocate(M1066a_V_radius_1066a(NR_1066A), &
          M1066a_V_density_1066a(NR_1066A), &
          M1066a_V_vp_1066a(NR_1066A), &
          M1066a_V_vs_1066a(NR_1066A), &
          M1066a_V_Qkappa_1066a(NR_1066A), &
          M1066a_V_Qmu_1066a(NR_1066A), &
          stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating M1066a_V arrays')

  ! all processes will define same parameters
  call define_model_1066a(CRUSTAL)

  end subroutine model_1066a_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_1066a(x,rho,vp,vs,Qkappa,Qmu,iregion_code)

  use constants
  use model_1066a_par

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
  do while(r >= M1066a_V_radius_1066a(i) .and. i /= NR_1066A)
    i = i + 1
  enddo

! make sure we stay in the right region and never take a point above
! and a point below the ICB or the CMB and interpolate between them,
! which would lead to a wrong value (keeping in mind that we interpolate
! between points i-1 and i below)
  if (iregion_code == IREGION_INNER_CORE .and. i > 33) i = 33

  if (iregion_code == IREGION_OUTER_CORE .and. i < 35) i = 35
  if (iregion_code == IREGION_OUTER_CORE .and. i > 66) i = 66

  if (iregion_code == IREGION_CRUST_MANTLE .and. i < 68) i = 68

  if (i == 1) then
    rho = M1066a_V_density_1066a(i)
    vp = M1066a_V_vp_1066a(i)
    vs = M1066a_V_vs_1066a(i)
    Qmu = M1066a_V_Qmu_1066a(i)
    Qkappa = M1066a_V_Qkappa_1066a(i)
  else

! interpolate from radius_1066a(i-1) to r using the values at i-1 and i
    frac = (r-M1066a_V_radius_1066a(i-1))/(M1066a_V_radius_1066a(i)-M1066a_V_radius_1066a(i-1))

    rho = M1066a_V_density_1066a(i-1) + frac * (M1066a_V_density_1066a(i)-M1066a_V_density_1066a(i-1))
    vp = M1066a_V_vp_1066a(i-1) + frac * (M1066a_V_vp_1066a(i)-M1066a_V_vp_1066a(i-1))
    vs = M1066a_V_vs_1066a(i-1) + frac * (M1066a_V_vs_1066a(i)-M1066a_V_vs_1066a(i-1))
    Qmu = M1066a_V_Qmu_1066a(i-1) + frac * (M1066a_V_Qmu_1066a(i)-M1066a_V_Qmu_1066a(i-1))
    Qkappa = M1066a_V_Qkappa_1066a(i-1) + frac * (M1066a_V_Qkappa_1066a(i)-M1066a_V_Qkappa_1066a(i-1))

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

  end subroutine model_1066a

!-------------------

  subroutine define_model_1066a(USE_EXTERNAL_CRUSTAL_MODEL)

  use constants
  use model_1066a_par

  implicit none

  logical :: USE_EXTERNAL_CRUSTAL_MODEL

  ! local parameters
  integer :: i

! define all the values in the model

  M1066a_V_radius_1066a(  1) =  0.d0
  M1066a_V_radius_1066a(  2) =   38400.d0
  M1066a_V_radius_1066a(  3) =   76810.d0
  M1066a_V_radius_1066a(  4) =   115210.d0
  M1066a_V_radius_1066a(  5) =   153610.d0
  M1066a_V_radius_1066a(  6) =   192020.d0
  M1066a_V_radius_1066a(  7) =   230420.d0
  M1066a_V_radius_1066a(  8) =   268820.d0
  M1066a_V_radius_1066a(  9) =   307220.d0
  M1066a_V_radius_1066a( 10) =   345630.d0
  M1066a_V_radius_1066a( 11) =   384030.d0
  M1066a_V_radius_1066a( 12) =   422430.d0
  M1066a_V_radius_1066a( 13) =   460840.d0
  M1066a_V_radius_1066a( 14) =   499240.d0
  M1066a_V_radius_1066a( 15) =   537640.d0
  M1066a_V_radius_1066a( 16) =   576050.d0
  M1066a_V_radius_1066a( 17) =   614450.d0
  M1066a_V_radius_1066a( 18) =   652850.d0
  M1066a_V_radius_1066a( 19) =   691260.d0
  M1066a_V_radius_1066a( 20) =   729660.d0
  M1066a_V_radius_1066a( 21) =   768060.d0
  M1066a_V_radius_1066a( 22) =   806460.d0
  M1066a_V_radius_1066a( 23) =   844870.d0
  M1066a_V_radius_1066a( 24) =   883270.d0
  M1066a_V_radius_1066a( 25) =   921670.d0
  M1066a_V_radius_1066a( 26) =   960080.d0
  M1066a_V_radius_1066a( 27) =   998480.d0
  M1066a_V_radius_1066a( 28) =   1036880.d0
  M1066a_V_radius_1066a( 29) =   1075290.d0
  M1066a_V_radius_1066a( 30) =   1113690.d0
  M1066a_V_radius_1066a( 31) =   1152090.d0
  M1066a_V_radius_1066a( 32) =   1190500.d0
  M1066a_V_radius_1066a( 33) =   1229480.d0
  M1066a_V_radius_1066a( 34) =   1229480.d0
  M1066a_V_radius_1066a( 35) =   1299360.d0
  M1066a_V_radius_1066a( 36) =   1369820.d0
  M1066a_V_radius_1066a( 37) =   1440280.d0
  M1066a_V_radius_1066a( 38) =   1510740.d0
  M1066a_V_radius_1066a( 39) =   1581190.d0
  M1066a_V_radius_1066a( 40) =   1651650.d0
  M1066a_V_radius_1066a( 41) =   1722110.d0
  M1066a_V_radius_1066a( 42) =   1792570.d0
  M1066a_V_radius_1066a( 43) =   1863030.d0
  M1066a_V_radius_1066a( 44) =   1933490.d0
  M1066a_V_radius_1066a( 45) =   2003950.d0
  M1066a_V_radius_1066a( 46) =   2074410.d0
  M1066a_V_radius_1066a( 47) =   2144870.d0
  M1066a_V_radius_1066a( 48) =   2215330.d0
  M1066a_V_radius_1066a( 49) =   2285790.d0
  M1066a_V_radius_1066a( 50) =   2356240.d0
  M1066a_V_radius_1066a( 51) =   2426700.d0
  M1066a_V_radius_1066a( 52) =   2497160.d0
  M1066a_V_radius_1066a( 53) =   2567620.d0
  M1066a_V_radius_1066a( 54) =   2638080.d0
  M1066a_V_radius_1066a( 55) =   2708540.d0
  M1066a_V_radius_1066a( 56) =   2779000.d0
  M1066a_V_radius_1066a( 57) =   2849460.d0
  M1066a_V_radius_1066a( 58) =   2919920.d0
  M1066a_V_radius_1066a( 59) =   2990380.d0
  M1066a_V_radius_1066a( 60) =   3060840.d0
  M1066a_V_radius_1066a( 61) =   3131300.d0
  M1066a_V_radius_1066a( 62) =   3201750.d0
  M1066a_V_radius_1066a( 63) =   3272210.d0
  M1066a_V_radius_1066a( 64) =   3342670.d0
  M1066a_V_radius_1066a( 65) =   3413130.d0
  M1066a_V_radius_1066a( 66) =   3484300.d0
  M1066a_V_radius_1066a( 67) =   3484300.d0
  M1066a_V_radius_1066a( 68) =   3518220.d0
  M1066a_V_radius_1066a( 69) =   3552850.d0
  M1066a_V_radius_1066a( 70) =   3587490.d0
  M1066a_V_radius_1066a( 71) =   3622120.d0
  M1066a_V_radius_1066a( 72) =   3656750.d0
  M1066a_V_radius_1066a( 73) =   3691380.d0
  M1066a_V_radius_1066a( 74) =   3726010.d0
  M1066a_V_radius_1066a( 75) =   3760640.d0
  M1066a_V_radius_1066a( 76) =   3795270.d0
  M1066a_V_radius_1066a( 77) =   3829910.d0
  M1066a_V_radius_1066a( 78) =   3864540.d0
  M1066a_V_radius_1066a( 79) =   3899170.d0
  M1066a_V_radius_1066a( 80) =   3933800.d0
  M1066a_V_radius_1066a( 81) =   3968430.d0
  M1066a_V_radius_1066a( 82) =   4003060.d0
  M1066a_V_radius_1066a( 83) =   4037690.d0
  M1066a_V_radius_1066a( 84) =   4072330.d0
  M1066a_V_radius_1066a( 85) =   4106960.d0
  M1066a_V_radius_1066a( 86) =   4141590.d0
  M1066a_V_radius_1066a( 87) =   4176220.d0
  M1066a_V_radius_1066a( 88) =   4210850.d0
  M1066a_V_radius_1066a( 89) =   4245480.d0
  M1066a_V_radius_1066a( 90) =   4280110.d0
  M1066a_V_radius_1066a( 91) =   4314740.d0
  M1066a_V_radius_1066a( 92) =   4349380.d0
  M1066a_V_radius_1066a( 93) =   4384010.d0
  M1066a_V_radius_1066a( 94) =   4418640.d0
  M1066a_V_radius_1066a( 95) =   4453270.d0
  M1066a_V_radius_1066a( 96) =   4487900.d0
  M1066a_V_radius_1066a( 97) =   4522530.d0
  M1066a_V_radius_1066a( 98) =   4557160.d0
  M1066a_V_radius_1066a( 99) =   4591800.d0
  M1066a_V_radius_1066a(100) =   4626430.d0
  M1066a_V_radius_1066a(101) =   4661060.d0
  M1066a_V_radius_1066a(102) =   4695690.d0
  M1066a_V_radius_1066a(103) =   4730320.d0
  M1066a_V_radius_1066a(104) =   4764950.d0
  M1066a_V_radius_1066a(105) =   4799580.d0
  M1066a_V_radius_1066a(106) =   4834220.d0
  M1066a_V_radius_1066a(107) =   4868850.d0
  M1066a_V_radius_1066a(108) =   4903480.d0
  M1066a_V_radius_1066a(109) =   4938110.d0
  M1066a_V_radius_1066a(110) =   4972740.d0
  M1066a_V_radius_1066a(111) =   5007370.d0
  M1066a_V_radius_1066a(112) =   5042000.d0
  M1066a_V_radius_1066a(113) =   5076640.d0
  M1066a_V_radius_1066a(114) =   5111270.d0
  M1066a_V_radius_1066a(115) =   5145900.d0
  M1066a_V_radius_1066a(116) =   5180530.d0
  M1066a_V_radius_1066a(117) =   5215160.d0
  M1066a_V_radius_1066a(118) =   5249790.d0
  M1066a_V_radius_1066a(119) =   5284420.d0
  M1066a_V_radius_1066a(120) =   5319060.d0
  M1066a_V_radius_1066a(121) =   5353690.d0
  M1066a_V_radius_1066a(122) =   5388320.d0
  M1066a_V_radius_1066a(123) =   5422950.d0
  M1066a_V_radius_1066a(124) =   5457580.d0
  M1066a_V_radius_1066a(125) =   5492210.d0
  M1066a_V_radius_1066a(126) =   5526840.d0
  M1066a_V_radius_1066a(127) =   5561470.d0
  M1066a_V_radius_1066a(128) =   5596110.d0
  M1066a_V_radius_1066a(129) =   5630740.d0
  M1066a_V_radius_1066a(130) =   5665370.d0
  M1066a_V_radius_1066a(131) =   5700000.d0
  M1066a_V_radius_1066a(132) =   5700000.d0
  M1066a_V_radius_1066a(133) =   5731250.d0
  M1066a_V_radius_1066a(134) =   5762500.d0
  M1066a_V_radius_1066a(135) =   5793750.d0
  M1066a_V_radius_1066a(136) =   5825000.d0
  M1066a_V_radius_1066a(137) =   5856250.d0
  M1066a_V_radius_1066a(138) =   5887500.d0
  M1066a_V_radius_1066a(139) =   5918750.d0
  M1066a_V_radius_1066a(140) =   5950000.d0
  M1066a_V_radius_1066a(141) =   5950000.d0
  M1066a_V_radius_1066a(142) =   5975630.d0
  M1066a_V_radius_1066a(143) =   6001250.d0
  M1066a_V_radius_1066a(144) =   6026880.d0
  M1066a_V_radius_1066a(145) =   6052500.d0
  M1066a_V_radius_1066a(146) =   6078130.d0
  M1066a_V_radius_1066a(147) =   6103750.d0
  M1066a_V_radius_1066a(148) =   6129380.d0
  M1066a_V_radius_1066a(149) =   6155000.d0
  M1066a_V_radius_1066a(150) =   6180630.d0
  M1066a_V_radius_1066a(151) =   6206250.d0
  M1066a_V_radius_1066a(152) =   6231880.d0
  M1066a_V_radius_1066a(153) =   6257500.d0
  M1066a_V_radius_1066a(154) =   6283130.d0
  M1066a_V_radius_1066a(155) =   6308750.d0
  M1066a_V_radius_1066a(156) =   6334380.d0
  M1066a_V_radius_1066a(157) =   6360000.d0
  M1066a_V_radius_1066a(158) =   6360000.d0
  M1066a_V_radius_1066a(159) =   6365500.d0
  M1066a_V_radius_1066a(160) =   6371000.d0

  M1066a_V_density_1066a(  1) =   13.42903d0
  M1066a_V_density_1066a(  2) =   13.42563d0
  M1066a_V_density_1066a(  3) =   13.41913d0
  M1066a_V_density_1066a(  4) =   13.41353d0
  M1066a_V_density_1066a(  5) =   13.40723d0
  M1066a_V_density_1066a(  6) =   13.40032d0
  M1066a_V_density_1066a(  7) =   13.39292d0
  M1066a_V_density_1066a(  8) =   13.38471d0
  M1066a_V_density_1066a(  9) =   13.3754d0
  M1066a_V_density_1066a( 10) =   13.3649d0
  M1066a_V_density_1066a( 11) =   13.35279d0
  M1066a_V_density_1066a( 12) =   13.33898d0
  M1066a_V_density_1066a( 13) =   13.32387d0
  M1066a_V_density_1066a( 14) =   13.30785d0
  M1066a_V_density_1066a( 15) =   13.29144d0
  M1066a_V_density_1066a( 16) =   13.27503d0
  M1066a_V_density_1066a( 17) =   13.25891d0
  M1066a_V_density_1066a( 18) =   13.2431d0
  M1066a_V_density_1066a( 19) =   13.22758d0
  M1066a_V_density_1066a( 20) =   13.21236d0
  M1066a_V_density_1066a( 21) =   13.19725d0
  M1066a_V_density_1066a( 22) =   13.18233d0
  M1066a_V_density_1066a( 23) =   13.16751d0
  M1066a_V_density_1066a( 24) =   13.15278d0
  M1066a_V_density_1066a( 25) =   13.13826d0
  M1066a_V_density_1066a( 26) =   13.12394d0
  M1066a_V_density_1066a( 27) =   13.10952d0
  M1066a_V_density_1066a( 28) =   13.09539d0
  M1066a_V_density_1066a( 29) =   13.08116d0
  M1066a_V_density_1066a( 30) =   13.06704d0
  M1066a_V_density_1066a( 31) =   13.05251d0
  M1066a_V_density_1066a( 32) =   13.03858d0
  M1066a_V_density_1066a( 33) =   13.02875d0
  M1066a_V_density_1066a( 34) =   12.16065d0
  M1066a_V_density_1066a( 35) =   12.11699d0
  M1066a_V_density_1066a( 36) =   12.07483d0
  M1066a_V_density_1066a( 37) =   12.03307d0
  M1066a_V_density_1066a( 38) =   11.9916d0
  M1066a_V_density_1066a( 39) =   11.95073d0
  M1066a_V_density_1066a( 40) =   11.91046d0
  M1066a_V_density_1066a( 41) =   11.86938d0
  M1066a_V_density_1066a( 42) =   11.82481d0
  M1066a_V_density_1066a( 43) =   11.77532d0
  M1066a_V_density_1066a( 44) =   11.72204d0
  M1066a_V_density_1066a( 45) =   11.66655d0
  M1066a_V_density_1066a( 46) =   11.60856d0
  M1066a_V_density_1066a( 47) =   11.54696d0
  M1066a_V_density_1066a( 48) =   11.48096d0
  M1066a_V_density_1066a( 49) =   11.41166d0
  M1066a_V_density_1066a( 50) =   11.34116d0
  M1066a_V_density_1066a( 51) =   11.27055d0
  M1066a_V_density_1066a( 52) =   11.19824d0
  M1066a_V_density_1066a( 53) =   11.12142d0
  M1066a_V_density_1066a( 54) =   11.03841d0
  M1066a_V_density_1066a( 55) =   10.95119d0
  M1066a_V_density_1066a( 56) =   10.86316d0
  M1066a_V_density_1066a( 57) =   10.77703d0
  M1066a_V_density_1066a( 58) =   10.6925d0
  M1066a_V_density_1066a( 59) =   10.60767d0
  M1066a_V_density_1066a( 60) =   10.52073d0
  M1066a_V_density_1066a( 61) =   10.4312d0
  M1066a_V_density_1066a( 62) =   10.33775d0
  M1066a_V_density_1066a( 63) =   10.23961d0
  M1066a_V_density_1066a( 64) =   10.13786d0
  M1066a_V_density_1066a( 65) =   10.0323d0
  M1066a_V_density_1066a( 66) =   9.91745d0
  M1066a_V_density_1066a( 67) =   5.53205d0
  M1066a_V_density_1066a( 68) =   5.52147d0
  M1066a_V_density_1066a( 69) =   5.50959d0
  M1066a_V_density_1066a( 70) =   5.49821d0
  M1066a_V_density_1066a( 71) =   5.48673d0
  M1066a_V_density_1066a( 72) =   5.47495d0
  M1066a_V_density_1066a( 73) =   5.46297d0
  M1066a_V_density_1066a( 74) =   5.45049d0
  M1066a_V_density_1066a( 75) =   5.43741d0
  M1066a_V_density_1066a( 76) =   5.42382d0
  M1066a_V_density_1066a( 77) =   5.40934d0
  M1066a_V_density_1066a( 78) =   5.39375d0
  M1066a_V_density_1066a( 79) =   5.37717d0
  M1066a_V_density_1066a( 80) =   5.35958d0
  M1066a_V_density_1066a( 81) =   5.34079d0
  M1066a_V_density_1066a( 82) =   5.321d0
  M1066a_V_density_1066a( 83) =   5.30031d0
  M1066a_V_density_1066a( 84) =   5.27902d0
  M1066a_V_density_1066a( 85) =   5.25733d0
  M1066a_V_density_1066a( 86) =   5.23554d0
  M1066a_V_density_1066a( 87) =   5.21375d0
  M1066a_V_density_1066a( 88) =   5.19196d0
  M1066a_V_density_1066a( 89) =   5.17056d0
  M1066a_V_density_1066a( 90) =   5.14937d0
  M1066a_V_density_1066a( 91) =   5.12827d0
  M1066a_V_density_1066a( 92) =   5.10758d0
  M1066a_V_density_1066a( 93) =   5.08728d0
  M1066a_V_density_1066a( 94) =   5.06738d0
  M1066a_V_density_1066a( 95) =   5.04769d0
  M1066a_V_density_1066a( 96) =   5.02809d0
  M1066a_V_density_1066a( 97) =   5.00869d0
  M1066a_V_density_1066a( 98) =   4.98929d0
  M1066a_V_density_1066a( 99) =   4.96968d0
  M1066a_V_density_1066a(100) =   4.95008d0
  M1066a_V_density_1066a(101) =   4.93048d0
  M1066a_V_density_1066a(102) =   4.91128d0
  M1066a_V_density_1066a(103) =   4.89257d0
  M1066a_V_density_1066a(104) =   4.87447d0
  M1066a_V_density_1066a(105) =   4.85716d0
  M1066a_V_density_1066a(106) =   4.84095d0
  M1066a_V_density_1066a(107) =   4.82554d0
  M1066a_V_density_1066a(108) =   4.81084d0
  M1066a_V_density_1066a(109) =   4.79683d0
  M1066a_V_density_1066a(110) =   4.78312d0
  M1066a_V_density_1066a(111) =   4.76951d0
  M1066a_V_density_1066a(112) =   4.7553d0
  M1066a_V_density_1066a(113) =   4.74008d0
  M1066a_V_density_1066a(114) =   4.72317d0
  M1066a_V_density_1066a(115) =   4.70426d0
  M1066a_V_density_1066a(116) =   4.68264d0
  M1066a_V_density_1066a(117) =   4.65863d0
  M1066a_V_density_1066a(118) =   4.63351d0
  M1066a_V_density_1066a(119) =   4.60859d0
  M1066a_V_density_1066a(120) =   4.58538d0
  M1066a_V_density_1066a(121) =   4.56536d0
  M1066a_V_density_1066a(122) =   4.55044d0
  M1066a_V_density_1066a(123) =   4.54072d0
  M1066a_V_density_1066a(124) =   4.5348d0
  M1066a_V_density_1066a(125) =   4.53478d0
  M1066a_V_density_1066a(126) =   4.53275d0
  M1066a_V_density_1066a(127) =   4.50893d0
  M1066a_V_density_1066a(128) =   4.46541d0
  M1066a_V_density_1066a(129) =   4.40098d0
  M1066a_V_density_1066a(130) =   4.31686d0
  M1066a_V_density_1066a(131) =   4.20553d0
  M1066a_V_density_1066a(132) =   4.20553d0
  M1066a_V_density_1066a(133) =   4.10272d0
  M1066a_V_density_1066a(134) =   4.0225d0
  M1066a_V_density_1066a(135) =   3.95789d0
  M1066a_V_density_1066a(136) =   3.89997d0
  M1066a_V_density_1066a(137) =   3.84675d0
  M1066a_V_density_1066a(138) =   3.80144d0
  M1066a_V_density_1066a(139) =   3.76072d0
  M1066a_V_density_1066a(140) =   3.7084d0
  M1066a_V_density_1066a(141) =   3.7084d0
  M1066a_V_density_1066a(142) =   3.6537d0
  M1066a_V_density_1066a(143) =   3.5964d0
  M1066a_V_density_1066a(144) =   3.54731d0
  M1066a_V_density_1066a(145) =   3.50511d0
  M1066a_V_density_1066a(146) =   3.46861d0
  M1066a_V_density_1066a(147) =   3.43851d0
  M1066a_V_density_1066a(148) =   3.41471d0
  M1066a_V_density_1066a(149) =   3.39751d0
  M1066a_V_density_1066a(150) =   3.3882d0
  M1066a_V_density_1066a(151) =   3.382d0
  M1066a_V_density_1066a(152) =   3.3745d0
  M1066a_V_density_1066a(153) =   3.3671d0
  M1066a_V_density_1066a(154) =   3.3598d0
  M1066a_V_density_1066a(155) =   3.35259d0
  M1066a_V_density_1066a(156) =   3.34549d0
  M1066a_V_density_1066a(157) =   3.33828d0
  M1066a_V_density_1066a(158) =   2.17798d0
  M1066a_V_density_1066a(159) =   2.17766d0
  M1066a_V_density_1066a(160) =   2.17734d0

  M1066a_V_vp_1066a(  1) =   11.3383d0
  M1066a_V_vp_1066a(  2) =   11.3374d0
  M1066a_V_vp_1066a(  3) =   11.3347d0
  M1066a_V_vp_1066a(  4) =   11.3301d0
  M1066a_V_vp_1066a(  5) =   11.3237d0
  M1066a_V_vp_1066a(  6) =   11.3155d0
  M1066a_V_vp_1066a(  7) =   11.3056d0
  M1066a_V_vp_1066a(  8) =   11.294d0
  M1066a_V_vp_1066a(  9) =   11.281d0
  M1066a_V_vp_1066a( 10) =   11.2666d0
  M1066a_V_vp_1066a( 11) =   11.2512d0
  M1066a_V_vp_1066a( 12) =   11.2349d0
  M1066a_V_vp_1066a( 13) =   11.2181d0
  M1066a_V_vp_1066a( 14) =   11.201d0
  M1066a_V_vp_1066a( 15) =   11.184d0
  M1066a_V_vp_1066a( 16) =   11.1672d0
  M1066a_V_vp_1066a( 17) =   11.1508d0
  M1066a_V_vp_1066a( 18) =   11.1351d0
  M1066a_V_vp_1066a( 19) =   11.1201d0
  M1066a_V_vp_1066a( 20) =   11.1059d0
  M1066a_V_vp_1066a( 21) =   11.0924d0
  M1066a_V_vp_1066a( 22) =   11.0798d0
  M1066a_V_vp_1066a( 23) =   11.0678d0
  M1066a_V_vp_1066a( 24) =   11.0564d0
  M1066a_V_vp_1066a( 25) =   11.0455d0
  M1066a_V_vp_1066a( 26) =   11.035d0
  M1066a_V_vp_1066a( 27) =   11.0248d0
  M1066a_V_vp_1066a( 28) =   11.0149d0
  M1066a_V_vp_1066a( 29) =   11.0051d0
  M1066a_V_vp_1066a( 30) =   10.9953d0
  M1066a_V_vp_1066a( 31) =   10.9857d0
  M1066a_V_vp_1066a( 32) =   10.9756d0
  M1066a_V_vp_1066a( 33) =   10.9687d0
  M1066a_V_vp_1066a( 34) =   10.414d0
  M1066a_V_vp_1066a( 35) =   10.3518d0
  M1066a_V_vp_1066a( 36) =   10.2922d0
  M1066a_V_vp_1066a( 37) =   10.2351d0
  M1066a_V_vp_1066a( 38) =   10.1808d0
  M1066a_V_vp_1066a( 39) =   10.1297d0
  M1066a_V_vp_1066a( 40) =   10.0788d0
  M1066a_V_vp_1066a( 41) =   10.0284d0
  M1066a_V_vp_1066a( 42) =   9.9788d0
  M1066a_V_vp_1066a( 43) =   9.9307d0
  M1066a_V_vp_1066a( 44) =   9.8836d0
  M1066a_V_vp_1066a( 45) =   9.8353d0
  M1066a_V_vp_1066a( 46) =   9.7825d0
  M1066a_V_vp_1066a( 47) =   9.7211d0
  M1066a_V_vp_1066a( 48) =   9.6521d0
  M1066a_V_vp_1066a( 49) =   9.5806d0
  M1066a_V_vp_1066a( 50) =   9.5115d0
  M1066a_V_vp_1066a( 51) =   9.4465d0
  M1066a_V_vp_1066a( 52) =   9.3828d0
  M1066a_V_vp_1066a( 53) =   9.3166d0
  M1066a_V_vp_1066a( 54) =   9.2442d0
  M1066a_V_vp_1066a( 55) =   9.1658d0
  M1066a_V_vp_1066a( 56) =   9.0833d0
  M1066a_V_vp_1066a( 57) =   8.9987d0
  M1066a_V_vp_1066a( 58) =   8.9116d0
  M1066a_V_vp_1066a( 59) =   8.8201d0
  M1066a_V_vp_1066a( 60) =   8.7223d0
  M1066a_V_vp_1066a( 61) =   8.6171d0
  M1066a_V_vp_1066a( 62) =   8.503d0
  M1066a_V_vp_1066a( 63) =   8.3807d0
  M1066a_V_vp_1066a( 64) =   8.2556d0
  M1066a_V_vp_1066a( 65) =   8.1318d0
  M1066a_V_vp_1066a( 66) =   8.0112d0
  M1066a_V_vp_1066a( 67) =   13.7172d0
  M1066a_V_vp_1066a( 68) =   13.7134d0
  M1066a_V_vp_1066a( 69) =   13.7089d0
  M1066a_V_vp_1066a( 70) =   13.6806d0
  M1066a_V_vp_1066a( 71) =   13.6517d0
  M1066a_V_vp_1066a( 72) =   13.6251d0
  M1066a_V_vp_1066a( 73) =   13.5916d0
  M1066a_V_vp_1066a( 74) =   13.5564d0
  M1066a_V_vp_1066a( 75) =   13.5165d0
  M1066a_V_vp_1066a( 76) =   13.4725d0
  M1066a_V_vp_1066a( 77) =   13.4248d0
  M1066a_V_vp_1066a( 78) =   13.3742d0
  M1066a_V_vp_1066a( 79) =   13.3216d0
  M1066a_V_vp_1066a( 80) =   13.2679d0
  M1066a_V_vp_1066a( 81) =   13.2142d0
  M1066a_V_vp_1066a( 82) =   13.1619d0
  M1066a_V_vp_1066a( 83) =   13.1114d0
  M1066a_V_vp_1066a( 84) =   13.0631d0
  M1066a_V_vp_1066a( 85) =   13.0174d0
  M1066a_V_vp_1066a( 86) =   12.9745d0
  M1066a_V_vp_1066a( 87) =   12.9346d0
  M1066a_V_vp_1066a( 88) =   12.8977d0
  M1066a_V_vp_1066a( 89) =   12.8635d0
  M1066a_V_vp_1066a( 90) =   12.8318d0
  M1066a_V_vp_1066a( 91) =   12.8022d0
  M1066a_V_vp_1066a( 92) =   12.7739d0
  M1066a_V_vp_1066a( 93) =   12.7463d0
  M1066a_V_vp_1066a( 94) =   12.7186d0
  M1066a_V_vp_1066a( 95) =   12.6903d0
  M1066a_V_vp_1066a( 96) =   12.661d0
  M1066a_V_vp_1066a( 97) =   12.6302d0
  M1066a_V_vp_1066a( 98) =   12.5978d0
  M1066a_V_vp_1066a( 99) =   12.5637d0
  M1066a_V_vp_1066a(100) =   12.5276d0
  M1066a_V_vp_1066a(101) =   12.4893d0
  M1066a_V_vp_1066a(102) =   12.4485d0
  M1066a_V_vp_1066a(103) =   12.4052d0
  M1066a_V_vp_1066a(104) =   12.3592d0
  M1066a_V_vp_1066a(105) =   12.3105d0
  M1066a_V_vp_1066a(106) =   12.2596d0
  M1066a_V_vp_1066a(107) =   12.2072d0
  M1066a_V_vp_1066a(108) =   12.1538d0
  M1066a_V_vp_1066a(109) =   12.0998d0
  M1066a_V_vp_1066a(110) =   12.0458d0
  M1066a_V_vp_1066a(111) =   11.992d0
  M1066a_V_vp_1066a(112) =   11.9373d0
  M1066a_V_vp_1066a(113) =   11.8804d0
  M1066a_V_vp_1066a(114) =   11.82d0
  M1066a_V_vp_1066a(115) =   11.7554d0
  M1066a_V_vp_1066a(116) =   11.6844d0
  M1066a_V_vp_1066a(117) =   11.6079d0
  M1066a_V_vp_1066a(118) =   11.5308d0
  M1066a_V_vp_1066a(119) =   11.4579d0
  M1066a_V_vp_1066a(120) =   11.3935d0
  M1066a_V_vp_1066a(121) =   11.3418d0
  M1066a_V_vp_1066a(122) =   11.3085d0
  M1066a_V_vp_1066a(123) =   11.2938d0
  M1066a_V_vp_1066a(124) =   11.2915d0
  M1066a_V_vp_1066a(125) =   11.3049d0
  M1066a_V_vp_1066a(126) =   11.3123d0
  M1066a_V_vp_1066a(127) =   11.2643d0
  M1066a_V_vp_1066a(128) =   11.1635d0
  M1066a_V_vp_1066a(129) =   11.0063d0
  M1066a_V_vp_1066a(130) =   10.7959d0
  M1066a_V_vp_1066a(131) =   10.5143d0
  M1066a_V_vp_1066a(132) =   10.5143d0
  M1066a_V_vp_1066a(133) =   10.2513d0
  M1066a_V_vp_1066a(134) =   10.0402d0
  M1066a_V_vp_1066a(135) =   9.8648d0
  M1066a_V_vp_1066a(136) =   9.7086d0
  M1066a_V_vp_1066a(137) =   9.5681d0
  M1066a_V_vp_1066a(138) =   9.4512d0
  M1066a_V_vp_1066a(139) =   9.351d0
  M1066a_V_vp_1066a(140) =   9.2283d0
  M1066a_V_vp_1066a(141) =   9.2283d0
  M1066a_V_vp_1066a(142) =   9.1087d0
  M1066a_V_vp_1066a(143) =   8.9823d0
  M1066a_V_vp_1066a(144) =   8.8592d0
  M1066a_V_vp_1066a(145) =   8.7386d0
  M1066a_V_vp_1066a(146) =   8.6193d0
  M1066a_V_vp_1066a(147) =   8.5018d0
  M1066a_V_vp_1066a(148) =   8.3871d0
  M1066a_V_vp_1066a(149) =   8.2736d0
  M1066a_V_vp_1066a(150) =   8.1585d0
  M1066a_V_vp_1066a(151) =   8.054d0
  M1066a_V_vp_1066a(152) =   7.9652d0
  M1066a_V_vp_1066a(153) =   7.8734d0
  M1066a_V_vp_1066a(154) =   7.7972d0
  M1066a_V_vp_1066a(155) =   7.7391d0
  M1066a_V_vp_1066a(156) =   7.7134d0
  M1066a_V_vp_1066a(157) =   7.7046d0
  M1066a_V_vp_1066a(158) =   4.7022d0
  M1066a_V_vp_1066a(159) =   4.7001d0
  M1066a_V_vp_1066a(160) =   4.6979d0

  M1066a_V_vs_1066a(  1) =   3.6298d0
  M1066a_V_vs_1066a(  2) =   3.6297d0
  M1066a_V_vs_1066a(  3) =   3.6294d0
  M1066a_V_vs_1066a(  4) =   3.6288d0
  M1066a_V_vs_1066a(  5) =   3.6281d0
  M1066a_V_vs_1066a(  6) =   3.6271d0
  M1066a_V_vs_1066a(  7) =   3.6259d0
  M1066a_V_vs_1066a(  8) =   3.6244d0
  M1066a_V_vs_1066a(  9) =   3.6228d0
  M1066a_V_vs_1066a( 10) =   3.6209d0
  M1066a_V_vs_1066a( 11) =   3.6187d0
  M1066a_V_vs_1066a( 12) =   3.6163d0
  M1066a_V_vs_1066a( 13) =   3.6137d0
  M1066a_V_vs_1066a( 14) =   3.6108d0
  M1066a_V_vs_1066a( 15) =   3.6076d0
  M1066a_V_vs_1066a( 16) =   3.6042d0
  M1066a_V_vs_1066a( 17) =   3.6004d0
  M1066a_V_vs_1066a( 18) =   3.5965d0
  M1066a_V_vs_1066a( 19) =   3.5922d0
  M1066a_V_vs_1066a( 20) =   3.5876d0
  M1066a_V_vs_1066a( 21) =   3.5828d0
  M1066a_V_vs_1066a( 22) =   3.5777d0
  M1066a_V_vs_1066a( 23) =   3.5724d0
  M1066a_V_vs_1066a( 24) =   3.5668d0
  M1066a_V_vs_1066a( 25) =   3.561d0
  M1066a_V_vs_1066a( 26) =   3.5551d0
  M1066a_V_vs_1066a( 27) =   3.549d0
  M1066a_V_vs_1066a( 28) =   3.5428d0
  M1066a_V_vs_1066a( 29) =   3.5365d0
  M1066a_V_vs_1066a( 30) =   3.5301d0
  M1066a_V_vs_1066a( 31) =   3.5238d0
  M1066a_V_vs_1066a( 32) =   3.5172d0
  M1066a_V_vs_1066a( 33) =   3.5118d0
  M1066a_V_vs_1066a( 34) =  0.d0
  M1066a_V_vs_1066a( 35) =  0.d0
  M1066a_V_vs_1066a( 36) =  0.d0
  M1066a_V_vs_1066a( 37) =  0.d0
  M1066a_V_vs_1066a( 38) =  0.d0
  M1066a_V_vs_1066a( 39) =  0.d0
  M1066a_V_vs_1066a( 40) =  0.d0
  M1066a_V_vs_1066a( 41) =  0.d0
  M1066a_V_vs_1066a( 42) =  0.d0
  M1066a_V_vs_1066a( 43) =  0.d0
  M1066a_V_vs_1066a( 44) =  0.d0
  M1066a_V_vs_1066a( 45) =  0.d0
  M1066a_V_vs_1066a( 46) =  0.d0
  M1066a_V_vs_1066a( 47) =  0.d0
  M1066a_V_vs_1066a( 48) =  0.d0
  M1066a_V_vs_1066a( 49) =  0.d0
  M1066a_V_vs_1066a( 50) =  0.d0
  M1066a_V_vs_1066a( 51) =  0.d0
  M1066a_V_vs_1066a( 52) =  0.d0
  M1066a_V_vs_1066a( 53) =  0.d0
  M1066a_V_vs_1066a( 54) =  0.d0
  M1066a_V_vs_1066a( 55) =  0.d0
  M1066a_V_vs_1066a( 56) =  0.d0
  M1066a_V_vs_1066a( 57) =  0.d0
  M1066a_V_vs_1066a( 58) =  0.d0
  M1066a_V_vs_1066a( 59) =  0.d0
  M1066a_V_vs_1066a( 60) =  0.d0
  M1066a_V_vs_1066a( 61) =  0.d0
  M1066a_V_vs_1066a( 62) =  0.d0
  M1066a_V_vs_1066a( 63) =  0.d0
  M1066a_V_vs_1066a( 64) =  0.d0
  M1066a_V_vs_1066a( 65) =  0.d0
  M1066a_V_vs_1066a( 66) =  0.d0
  M1066a_V_vs_1066a( 67) =   7.2498d0
  M1066a_V_vs_1066a( 68) =   7.2376d0
  M1066a_V_vs_1066a( 69) =   7.2239d0
  M1066a_V_vs_1066a( 70) =   7.21d0
  M1066a_V_vs_1066a( 71) =   7.1964d0
  M1066a_V_vs_1066a( 72) =   7.183d0
  M1066a_V_vs_1066a( 73) =   7.1699d0
  M1066a_V_vs_1066a( 74) =   7.1571d0
  M1066a_V_vs_1066a( 75) =   7.1445d0
  M1066a_V_vs_1066a( 76) =   7.132d0
  M1066a_V_vs_1066a( 77) =   7.1196d0
  M1066a_V_vs_1066a( 78) =   7.1074d0
  M1066a_V_vs_1066a( 79) =   7.0953d0
  M1066a_V_vs_1066a( 80) =   7.0832d0
  M1066a_V_vs_1066a( 81) =   7.0712d0
  M1066a_V_vs_1066a( 82) =   7.0592d0
  M1066a_V_vs_1066a( 83) =   7.0471d0
  M1066a_V_vs_1066a( 84) =   7.0347d0
  M1066a_V_vs_1066a( 85) =   7.0219d0
  M1066a_V_vs_1066a( 86) =   7.0086d0
  M1066a_V_vs_1066a( 87) =   6.9947d0
  M1066a_V_vs_1066a( 88) =   6.9803d0
  M1066a_V_vs_1066a( 89) =   6.9651d0
  M1066a_V_vs_1066a( 90) =   6.9493d0
  M1066a_V_vs_1066a( 91) =   6.9329d0
  M1066a_V_vs_1066a( 92) =   6.9162d0
  M1066a_V_vs_1066a( 93) =   6.8991d0
  M1066a_V_vs_1066a( 94) =   6.882d0
  M1066a_V_vs_1066a( 95) =   6.8652d0
  M1066a_V_vs_1066a( 96) =   6.849d0
  M1066a_V_vs_1066a( 97) =   6.8334d0
  M1066a_V_vs_1066a( 98) =   6.8182d0
  M1066a_V_vs_1066a( 99) =   6.8036d0
  M1066a_V_vs_1066a(100) =   6.7891d0
  M1066a_V_vs_1066a(101) =   6.7744d0
  M1066a_V_vs_1066a(102) =   6.7589d0
  M1066a_V_vs_1066a(103) =   6.7427d0
  M1066a_V_vs_1066a(104) =   6.7255d0
  M1066a_V_vs_1066a(105) =   6.7073d0
  M1066a_V_vs_1066a(106) =   6.6881d0
  M1066a_V_vs_1066a(107) =   6.6684d0
  M1066a_V_vs_1066a(108) =   6.6485d0
  M1066a_V_vs_1066a(109) =   6.6288d0
  M1066a_V_vs_1066a(110) =   6.6095d0
  M1066a_V_vs_1066a(111) =   6.5911d0
  M1066a_V_vs_1066a(112) =   6.5731d0
  M1066a_V_vs_1066a(113) =   6.5548d0
  M1066a_V_vs_1066a(114) =   6.5351d0
  M1066a_V_vs_1066a(115) =   6.5133d0
  M1066a_V_vs_1066a(116) =   6.4881d0
  M1066a_V_vs_1066a(117) =   6.4594d0
  M1066a_V_vs_1066a(118) =   6.4286d0
  M1066a_V_vs_1066a(119) =   6.3976d0
  M1066a_V_vs_1066a(120) =   6.3684d0
  M1066a_V_vs_1066a(121) =   6.3428d0
  M1066a_V_vs_1066a(122) =   6.3235d0
  M1066a_V_vs_1066a(123) =   6.3114d0
  M1066a_V_vs_1066a(124) =   6.3041d0
  M1066a_V_vs_1066a(125) =   6.3052d0
  M1066a_V_vs_1066a(126) =   6.3021d0
  M1066a_V_vs_1066a(127) =   6.2643d0
  M1066a_V_vs_1066a(128) =   6.1947d0
  M1066a_V_vs_1066a(129) =   6.0912d0
  M1066a_V_vs_1066a(130) =   5.9555d0
  M1066a_V_vs_1066a(131) =   5.7755d0
  M1066a_V_vs_1066a(132) =   5.7755d0
  M1066a_V_vs_1066a(133) =   5.6083d0
  M1066a_V_vs_1066a(134) =   5.4752d0
  M1066a_V_vs_1066a(135) =   5.3653d0
  M1066a_V_vs_1066a(136) =   5.2665d0
  M1066a_V_vs_1066a(137) =   5.1762d0
  M1066a_V_vs_1066a(138) =   5.0996d0
  M1066a_V_vs_1066a(139) =   5.0322d0
  M1066a_V_vs_1066a(140) =   4.9488d0
  M1066a_V_vs_1066a(141) =   4.9488d0
  M1066a_V_vs_1066a(142) =   4.8667d0
  M1066a_V_vs_1066a(143) =   4.7806d0
  M1066a_V_vs_1066a(144) =   4.6995d0
  M1066a_V_vs_1066a(145) =   4.6211d0
  M1066a_V_vs_1066a(146) =   4.5479d0
  M1066a_V_vs_1066a(147) =   4.4882d0
  M1066a_V_vs_1066a(148) =   4.4421d0
  M1066a_V_vs_1066a(149) =   4.4084d0
  M1066a_V_vs_1066a(150) =   4.3874d0
  M1066a_V_vs_1066a(151) =   4.3795d0
  M1066a_V_vs_1066a(152) =   4.3904d0
  M1066a_V_vs_1066a(153) =   4.4331d0
  M1066a_V_vs_1066a(154) =   4.483d0
  M1066a_V_vs_1066a(155) =   4.5389d0
  M1066a_V_vs_1066a(156) =   4.604d0
  M1066a_V_vs_1066a(157) =   4.6487d0
  M1066a_V_vs_1066a(158) =   2.5806d0
  M1066a_V_vs_1066a(159) =   2.5814d0
  M1066a_V_vs_1066a(160) =   2.5822d0

  if (SUPPRESS_CRUSTAL_MESH) then
    M1066a_V_vp_1066a(158:160) = M1066a_V_vp_1066a(157)
    M1066a_V_vs_1066a(158:160) = M1066a_V_vs_1066a(157)
    M1066a_V_density_1066a(158:160) = M1066a_V_density_1066a(157)
  endif

  M1066a_V_Qkappa_1066a(  1) =   156900.d0
  M1066a_V_Qkappa_1066a(  2) =   156900.d0
  M1066a_V_Qkappa_1066a(  3) =   156900.d0
  M1066a_V_Qkappa_1066a(  4) =   156900.d0
  M1066a_V_Qkappa_1066a(  5) =   156900.d0
  M1066a_V_Qkappa_1066a(  6) =   156900.d0
  M1066a_V_Qkappa_1066a(  7) =   156900.d0
  M1066a_V_Qkappa_1066a(  8) =   156900.d0
  M1066a_V_Qkappa_1066a(  9) =   156900.d0
  M1066a_V_Qkappa_1066a( 10) =   156900.d0
  M1066a_V_Qkappa_1066a( 11) =   156900.d0
  M1066a_V_Qkappa_1066a( 12) =   156900.d0
  M1066a_V_Qkappa_1066a( 13) =   156900.d0
  M1066a_V_Qkappa_1066a( 14) =   156900.d0
  M1066a_V_Qkappa_1066a( 15) =   156900.d0
  M1066a_V_Qkappa_1066a( 16) =   156900.d0
  M1066a_V_Qkappa_1066a( 17) =   156900.d0
  M1066a_V_Qkappa_1066a( 18) =   156900.d0
  M1066a_V_Qkappa_1066a( 19) =   156900.d0
  M1066a_V_Qkappa_1066a( 20) =   156900.d0
  M1066a_V_Qkappa_1066a( 21) =   156900.d0
  M1066a_V_Qkappa_1066a( 22) =   156900.d0
  M1066a_V_Qkappa_1066a( 23) =   156900.d0
  M1066a_V_Qkappa_1066a( 24) =   156900.d0
  M1066a_V_Qkappa_1066a( 25) =   156900.d0
  M1066a_V_Qkappa_1066a( 26) =   156900.d0
  M1066a_V_Qkappa_1066a( 27) =   156900.d0
  M1066a_V_Qkappa_1066a( 28) =   156900.d0
  M1066a_V_Qkappa_1066a( 29) =   156900.d0
  M1066a_V_Qkappa_1066a( 30) =   156900.d0
  M1066a_V_Qkappa_1066a( 31) =   156900.d0
  M1066a_V_Qkappa_1066a( 32) =   156900.d0
  M1066a_V_Qkappa_1066a( 33) =   156900.d0
  M1066a_V_Qkappa_1066a( 34) =  0.d0
  M1066a_V_Qkappa_1066a( 35) =  0.d0
  M1066a_V_Qkappa_1066a( 36) =  0.d0
  M1066a_V_Qkappa_1066a( 37) =  0.d0
  M1066a_V_Qkappa_1066a( 38) =  0.d0
  M1066a_V_Qkappa_1066a( 39) =  0.d0
  M1066a_V_Qkappa_1066a( 40) =  0.d0
  M1066a_V_Qkappa_1066a( 41) =  0.d0
  M1066a_V_Qkappa_1066a( 42) =  0.d0
  M1066a_V_Qkappa_1066a( 43) =  0.d0
  M1066a_V_Qkappa_1066a( 44) =  0.d0
  M1066a_V_Qkappa_1066a( 45) =  0.d0
  M1066a_V_Qkappa_1066a( 46) =  0.d0
  M1066a_V_Qkappa_1066a( 47) =  0.d0
  M1066a_V_Qkappa_1066a( 48) =  0.d0
  M1066a_V_Qkappa_1066a( 49) =  0.d0
  M1066a_V_Qkappa_1066a( 50) =  0.d0
  M1066a_V_Qkappa_1066a( 51) =  0.d0
  M1066a_V_Qkappa_1066a( 52) =  0.d0
  M1066a_V_Qkappa_1066a( 53) =  0.d0
  M1066a_V_Qkappa_1066a( 54) =  0.d0
  M1066a_V_Qkappa_1066a( 55) =  0.d0
  M1066a_V_Qkappa_1066a( 56) =  0.d0
  M1066a_V_Qkappa_1066a( 57) =  0.d0
  M1066a_V_Qkappa_1066a( 58) =  0.d0
  M1066a_V_Qkappa_1066a( 59) =  0.d0
  M1066a_V_Qkappa_1066a( 60) =  0.d0
  M1066a_V_Qkappa_1066a( 61) =  0.d0
  M1066a_V_Qkappa_1066a( 62) =  0.d0
  M1066a_V_Qkappa_1066a( 63) =  0.d0
  M1066a_V_Qkappa_1066a( 64) =  0.d0
  M1066a_V_Qkappa_1066a( 65) =  0.d0
  M1066a_V_Qkappa_1066a( 66) =  0.d0
  M1066a_V_Qkappa_1066a( 67) =   16600.d0
  M1066a_V_Qkappa_1066a( 68) =   16600.d0
  M1066a_V_Qkappa_1066a( 69) =   16600.d0
  M1066a_V_Qkappa_1066a( 70) =   16600.d0
  M1066a_V_Qkappa_1066a( 71) =   16600.d0
  M1066a_V_Qkappa_1066a( 72) =   16600.d0
  M1066a_V_Qkappa_1066a( 73) =   16600.d0
  M1066a_V_Qkappa_1066a( 74) =   16600.d0
  M1066a_V_Qkappa_1066a( 75) =   16600.d0
  M1066a_V_Qkappa_1066a( 76) =   16600.d0
  M1066a_V_Qkappa_1066a( 77) =   16600.d0
  M1066a_V_Qkappa_1066a( 78) =   16600.d0
  M1066a_V_Qkappa_1066a( 79) =   16600.d0
  M1066a_V_Qkappa_1066a( 80) =   16600.d0
  M1066a_V_Qkappa_1066a( 81) =   16600.d0
  M1066a_V_Qkappa_1066a( 82) =   16600.d0
  M1066a_V_Qkappa_1066a( 83) =   16600.d0
  M1066a_V_Qkappa_1066a( 84) =   16600.d0
  M1066a_V_Qkappa_1066a( 85) =   16600.d0
  M1066a_V_Qkappa_1066a( 86) =   16600.d0
  M1066a_V_Qkappa_1066a( 87) =   16600.d0
  M1066a_V_Qkappa_1066a( 88) =   16600.d0
  M1066a_V_Qkappa_1066a( 89) =   16600.d0
  M1066a_V_Qkappa_1066a( 90) =   16600.d0
  M1066a_V_Qkappa_1066a( 91) =   16600.d0
  M1066a_V_Qkappa_1066a( 92) =   16600.d0
  M1066a_V_Qkappa_1066a( 93) =   16600.d0
  M1066a_V_Qkappa_1066a( 94) =   16600.d0
  M1066a_V_Qkappa_1066a( 95) =   16600.d0
  M1066a_V_Qkappa_1066a( 96) =   16600.d0
  M1066a_V_Qkappa_1066a( 97) =   16600.d0
  M1066a_V_Qkappa_1066a( 98) =   16600.d0
  M1066a_V_Qkappa_1066a( 99) =   16600.d0
  M1066a_V_Qkappa_1066a(100) =   16600.d0
  M1066a_V_Qkappa_1066a(101) =   16600.d0
  M1066a_V_Qkappa_1066a(102) =   16600.d0
  M1066a_V_Qkappa_1066a(103) =   16600.d0
  M1066a_V_Qkappa_1066a(104) =   16600.d0
  M1066a_V_Qkappa_1066a(105) =   16600.d0
  M1066a_V_Qkappa_1066a(106) =   16600.d0
  M1066a_V_Qkappa_1066a(107) =   16600.d0
  M1066a_V_Qkappa_1066a(108) =   16600.d0
  M1066a_V_Qkappa_1066a(109) =   16600.d0
  M1066a_V_Qkappa_1066a(110) =   16600.d0
  M1066a_V_Qkappa_1066a(111) =   16600.d0
  M1066a_V_Qkappa_1066a(112) =   16600.d0
  M1066a_V_Qkappa_1066a(113) =   16600.d0
  M1066a_V_Qkappa_1066a(114) =   16600.d0
  M1066a_V_Qkappa_1066a(115) =   16600.d0
  M1066a_V_Qkappa_1066a(116) =   16600.d0
  M1066a_V_Qkappa_1066a(117) =   16600.d0
  M1066a_V_Qkappa_1066a(118) =   16600.d0
  M1066a_V_Qkappa_1066a(119) =   16600.d0
  M1066a_V_Qkappa_1066a(120) =   16600.d0
  M1066a_V_Qkappa_1066a(121) =   16600.d0
  M1066a_V_Qkappa_1066a(122) =   16600.d0
  M1066a_V_Qkappa_1066a(123) =   16600.d0
  M1066a_V_Qkappa_1066a(124) =   16600.d0
  M1066a_V_Qkappa_1066a(125) =   16600.d0
  M1066a_V_Qkappa_1066a(126) =   16600.d0
  M1066a_V_Qkappa_1066a(127) =   16600.d0
  M1066a_V_Qkappa_1066a(128) =   16600.d0
  M1066a_V_Qkappa_1066a(129) =   16600.d0
  M1066a_V_Qkappa_1066a(130) =   16600.d0
  M1066a_V_Qkappa_1066a(131) =   16600.d0
  M1066a_V_Qkappa_1066a(132) =   13840.d0
  M1066a_V_Qkappa_1066a(133) =   13840.d0
  M1066a_V_Qkappa_1066a(134) =   13840.d0
  M1066a_V_Qkappa_1066a(135) =   13840.d0
  M1066a_V_Qkappa_1066a(136) =   13840.d0
  M1066a_V_Qkappa_1066a(137) =   13840.d0
  M1066a_V_Qkappa_1066a(138) =   13840.d0
  M1066a_V_Qkappa_1066a(139) =   13840.d0
  M1066a_V_Qkappa_1066a(140) =   13840.d0
  M1066a_V_Qkappa_1066a(141) =   5893.d0
  M1066a_V_Qkappa_1066a(142) =   5893.d0
  M1066a_V_Qkappa_1066a(143) =   5893.d0
  M1066a_V_Qkappa_1066a(144) =   5893.d0
  M1066a_V_Qkappa_1066a(145) =   5893.d0
  M1066a_V_Qkappa_1066a(146) =   5893.d0
  M1066a_V_Qkappa_1066a(147) =   5893.d0
  M1066a_V_Qkappa_1066a(148) =   5893.d0
  M1066a_V_Qkappa_1066a(149) =   5893.d0
  M1066a_V_Qkappa_1066a(150) =   5893.d0
  M1066a_V_Qkappa_1066a(151) =   5893.d0
  M1066a_V_Qkappa_1066a(152) =   5893.d0
  M1066a_V_Qkappa_1066a(153) =   5893.d0
  M1066a_V_Qkappa_1066a(154) =   5893.d0
  M1066a_V_Qkappa_1066a(155) =   5893.d0
  M1066a_V_Qkappa_1066a(156) =   5893.d0
  M1066a_V_Qkappa_1066a(157) =   5893.d0
  M1066a_V_Qkappa_1066a(158) =   5893.d0
  M1066a_V_Qkappa_1066a(159) =   5893.d0
  M1066a_V_Qkappa_1066a(160) =   5893.d0

  M1066a_V_Qmu_1066a(  1) =   3138.d0
  M1066a_V_Qmu_1066a(  2) =   3138.d0
  M1066a_V_Qmu_1066a(  3) =   3138.d0
  M1066a_V_Qmu_1066a(  4) =   3138.d0
  M1066a_V_Qmu_1066a(  5) =   3138.d0
  M1066a_V_Qmu_1066a(  6) =   3138.d0
  M1066a_V_Qmu_1066a(  7) =   3138.d0
  M1066a_V_Qmu_1066a(  8) =   3138.d0
  M1066a_V_Qmu_1066a(  9) =   3138.d0
  M1066a_V_Qmu_1066a( 10) =   3138.d0
  M1066a_V_Qmu_1066a( 11) =   3138.d0
  M1066a_V_Qmu_1066a( 12) =   3138.d0
  M1066a_V_Qmu_1066a( 13) =   3138.d0
  M1066a_V_Qmu_1066a( 14) =   3138.d0
  M1066a_V_Qmu_1066a( 15) =   3138.d0
  M1066a_V_Qmu_1066a( 16) =   3138.d0
  M1066a_V_Qmu_1066a( 17) =   3138.d0
  M1066a_V_Qmu_1066a( 18) =   3138.d0
  M1066a_V_Qmu_1066a( 19) =   3138.d0
  M1066a_V_Qmu_1066a( 20) =   3138.d0
  M1066a_V_Qmu_1066a( 21) =   3138.d0
  M1066a_V_Qmu_1066a( 22) =   3138.d0
  M1066a_V_Qmu_1066a( 23) =   3138.d0
  M1066a_V_Qmu_1066a( 24) =   3138.d0
  M1066a_V_Qmu_1066a( 25) =   3138.d0
  M1066a_V_Qmu_1066a( 26) =   3138.d0
  M1066a_V_Qmu_1066a( 27) =   3138.d0
  M1066a_V_Qmu_1066a( 28) =   3138.d0
  M1066a_V_Qmu_1066a( 29) =   3138.d0
  M1066a_V_Qmu_1066a( 30) =   3138.d0
  M1066a_V_Qmu_1066a( 31) =   3138.d0
  M1066a_V_Qmu_1066a( 32) =   3138.d0
  M1066a_V_Qmu_1066a( 33) =   3138.d0
  M1066a_V_Qmu_1066a( 34) =  0.d0
  M1066a_V_Qmu_1066a( 35) =  0.d0
  M1066a_V_Qmu_1066a( 36) =  0.d0
  M1066a_V_Qmu_1066a( 37) =  0.d0
  M1066a_V_Qmu_1066a( 38) =  0.d0
  M1066a_V_Qmu_1066a( 39) =  0.d0
  M1066a_V_Qmu_1066a( 40) =  0.d0
  M1066a_V_Qmu_1066a( 41) =  0.d0
  M1066a_V_Qmu_1066a( 42) =  0.d0
  M1066a_V_Qmu_1066a( 43) =  0.d0
  M1066a_V_Qmu_1066a( 44) =  0.d0
  M1066a_V_Qmu_1066a( 45) =  0.d0
  M1066a_V_Qmu_1066a( 46) =  0.d0
  M1066a_V_Qmu_1066a( 47) =  0.d0
  M1066a_V_Qmu_1066a( 48) =  0.d0
  M1066a_V_Qmu_1066a( 49) =  0.d0
  M1066a_V_Qmu_1066a( 50) =  0.d0
  M1066a_V_Qmu_1066a( 51) =  0.d0
  M1066a_V_Qmu_1066a( 52) =  0.d0
  M1066a_V_Qmu_1066a( 53) =  0.d0
  M1066a_V_Qmu_1066a( 54) =  0.d0
  M1066a_V_Qmu_1066a( 55) =  0.d0
  M1066a_V_Qmu_1066a( 56) =  0.d0
  M1066a_V_Qmu_1066a( 57) =  0.d0
  M1066a_V_Qmu_1066a( 58) =  0.d0
  M1066a_V_Qmu_1066a( 59) =  0.d0
  M1066a_V_Qmu_1066a( 60) =  0.d0
  M1066a_V_Qmu_1066a( 61) =  0.d0
  M1066a_V_Qmu_1066a( 62) =  0.d0
  M1066a_V_Qmu_1066a( 63) =  0.d0
  M1066a_V_Qmu_1066a( 64) =  0.d0
  M1066a_V_Qmu_1066a( 65) =  0.d0
  M1066a_V_Qmu_1066a( 66) =  0.d0
  M1066a_V_Qmu_1066a( 67) =   332.d0
  M1066a_V_Qmu_1066a( 68) =   332.d0
  M1066a_V_Qmu_1066a( 69) =   332.d0
  M1066a_V_Qmu_1066a( 70) =   332.d0
  M1066a_V_Qmu_1066a( 71) =   332.d0
  M1066a_V_Qmu_1066a( 72) =   332.d0
  M1066a_V_Qmu_1066a( 73) =   332.d0
  M1066a_V_Qmu_1066a( 74) =   332.d0
  M1066a_V_Qmu_1066a( 75) =   332.d0
  M1066a_V_Qmu_1066a( 76) =   332.d0
  M1066a_V_Qmu_1066a( 77) =   332.d0
  M1066a_V_Qmu_1066a( 78) =   332.d0
  M1066a_V_Qmu_1066a( 79) =   332.d0
  M1066a_V_Qmu_1066a( 80) =   332.d0
  M1066a_V_Qmu_1066a( 81) =   332.d0
  M1066a_V_Qmu_1066a( 82) =   332.d0
  M1066a_V_Qmu_1066a( 83) =   332.d0
  M1066a_V_Qmu_1066a( 84) =   332.d0
  M1066a_V_Qmu_1066a( 85) =   332.d0
  M1066a_V_Qmu_1066a( 86) =   332.d0
  M1066a_V_Qmu_1066a( 87) =   332.d0
  M1066a_V_Qmu_1066a( 88) =   332.d0
  M1066a_V_Qmu_1066a( 89) =   332.d0
  M1066a_V_Qmu_1066a( 90) =   332.d0
  M1066a_V_Qmu_1066a( 91) =   332.d0
  M1066a_V_Qmu_1066a( 92) =   332.d0
  M1066a_V_Qmu_1066a( 93) =   332.d0
  M1066a_V_Qmu_1066a( 94) =   332.d0
  M1066a_V_Qmu_1066a( 95) =   332.d0
  M1066a_V_Qmu_1066a( 96) =   332.d0
  M1066a_V_Qmu_1066a( 97) =   332.d0
  M1066a_V_Qmu_1066a( 98) =   332.d0
  M1066a_V_Qmu_1066a( 99) =   332.d0
  M1066a_V_Qmu_1066a(100) =   332.d0
  M1066a_V_Qmu_1066a(101) =   332.d0
  M1066a_V_Qmu_1066a(102) =   332.d0
  M1066a_V_Qmu_1066a(103) =   332.d0
  M1066a_V_Qmu_1066a(104) =   332.d0
  M1066a_V_Qmu_1066a(105) =   332.d0
  M1066a_V_Qmu_1066a(106) =   332.d0
  M1066a_V_Qmu_1066a(107) =   332.d0
  M1066a_V_Qmu_1066a(108) =   332.d0
  M1066a_V_Qmu_1066a(109) =   332.d0
  M1066a_V_Qmu_1066a(110) =   332.d0
  M1066a_V_Qmu_1066a(111) =   332.d0
  M1066a_V_Qmu_1066a(112) =   332.d0
  M1066a_V_Qmu_1066a(113) =   332.d0
  M1066a_V_Qmu_1066a(114) =   332.d0
  M1066a_V_Qmu_1066a(115) =   332.d0
  M1066a_V_Qmu_1066a(116) =   332.d0
  M1066a_V_Qmu_1066a(117) =   332.d0
  M1066a_V_Qmu_1066a(118) =   332.d0
  M1066a_V_Qmu_1066a(119) =   332.d0
  M1066a_V_Qmu_1066a(120) =   332.d0
  M1066a_V_Qmu_1066a(121) =   332.d0
  M1066a_V_Qmu_1066a(122) =   332.d0
  M1066a_V_Qmu_1066a(123) =   332.d0
  M1066a_V_Qmu_1066a(124) =   332.d0
  M1066a_V_Qmu_1066a(125) =   332.d0
  M1066a_V_Qmu_1066a(126) =   332.d0
  M1066a_V_Qmu_1066a(127) =   332.d0
  M1066a_V_Qmu_1066a(128) =   332.d0
  M1066a_V_Qmu_1066a(129) =   332.d0
  M1066a_V_Qmu_1066a(130) =   332.d0
  M1066a_V_Qmu_1066a(131) =   332.d0
  M1066a_V_Qmu_1066a(132) =   276.8d0
  M1066a_V_Qmu_1066a(133) =   276.8d0
  M1066a_V_Qmu_1066a(134) =   276.8d0
  M1066a_V_Qmu_1066a(135) =   276.8d0
  M1066a_V_Qmu_1066a(136) =   276.8d0
  M1066a_V_Qmu_1066a(137) =   276.8d0
  M1066a_V_Qmu_1066a(138) =   276.8d0
  M1066a_V_Qmu_1066a(139) =   276.8d0
  M1066a_V_Qmu_1066a(140) =   276.8d0
  M1066a_V_Qmu_1066a(141) =   117.9d0
  M1066a_V_Qmu_1066a(142) =   117.9d0
  M1066a_V_Qmu_1066a(143) =   117.9d0
  M1066a_V_Qmu_1066a(144) =   117.9d0
  M1066a_V_Qmu_1066a(145) =   117.9d0
  M1066a_V_Qmu_1066a(146) =   117.9d0
  M1066a_V_Qmu_1066a(147) =   117.9d0
  M1066a_V_Qmu_1066a(148) =   117.9d0
  M1066a_V_Qmu_1066a(149) =   117.9d0
  M1066a_V_Qmu_1066a(150) =   117.9d0
  M1066a_V_Qmu_1066a(151) =   117.9d0
  M1066a_V_Qmu_1066a(152) =   117.9d0
  M1066a_V_Qmu_1066a(153) =   117.9d0
  M1066a_V_Qmu_1066a(154) =   117.9d0
  M1066a_V_Qmu_1066a(155) =   117.9d0
  M1066a_V_Qmu_1066a(156) =   117.9d0
  M1066a_V_Qmu_1066a(157) =   117.9d0
  M1066a_V_Qmu_1066a(158) =   117.9d0
  M1066a_V_Qmu_1066a(159) =   117.9d0
  M1066a_V_Qmu_1066a(160) =   117.9d0

! strip the crust and replace it by mantle if we use an external crustal model
  if (SUPPRESS_CRUSTAL_MESH .or. USE_EXTERNAL_CRUSTAL_MODEL) then
    do i=NR_1066A-3,NR_1066A
      M1066a_V_density_1066a(i) = M1066a_V_density_1066a(NR_1066A-4)
      M1066a_V_vp_1066a(i) = M1066a_V_vp_1066a(NR_1066A-4)
      M1066a_V_vs_1066a(i) = M1066a_V_vs_1066a(NR_1066A-4)
      M1066a_V_Qkappa_1066a(i) = M1066a_V_Qkappa_1066a(NR_1066A-4)
      M1066a_V_Qmu_1066a(i) = M1066a_V_Qmu_1066a(NR_1066A-4)
    enddo
  endif

  end subroutine define_model_1066a

