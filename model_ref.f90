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


  subroutine model_ref(x,rho,vpv,vph,vsv,vsh,eta,Qkappa,Qmu,iregion_code,Mref_V)

  implicit none

  include "constants.h"

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

! input:
! dimensionless radius x

! output: non-dimensionalized
! mass density rho
! compressional wave speed vpv
! compressional wave speed vph
! shear wave speed vsv
! shear wave speed vsh
! dimensionless parameter eta
! shear quality factor Qmu
! bulk quality factor Qkappa

  integer iregion_code

  double precision x,rho,vpv,vph,vsv,vsh,eta,Qmu,Qkappa

  integer i

  double precision r,frac,scaleval

! compute real physical radius in meters
  r = x * R_EARTH

  i = 1
  do while(r >= Mref_V%radius_ref(i) .and. i /= NR_REF)
    i = i + 1
  enddo

! make sure we stay in the right region
  if(iregion_code == IREGION_INNER_CORE .and. i > 180) i = 180

  if(iregion_code == IREGION_OUTER_CORE .and. i < 181) i = 181
  if(iregion_code == IREGION_OUTER_CORE .and. i > 358) i = 358

  if(iregion_code == IREGION_CRUST_MANTLE .and. i < 359) i = 359

  if(i == 1) then
    rho = Mref_V%density_ref(i)
    vpv = Mref_V%vpv_ref(i)
    vph = Mref_V%vph_ref(i)
    vsv = Mref_V%vsv_ref(i)
    vsh = Mref_V%vsh_ref(i)
    eta = Mref_V%eta_ref(i)
    Qkappa = Mref_V%Qkappa_ref(i)
    Qmu = Mref_V%Qmu_ref(i)
  else

! interpolate from radius_ref(i-1) to r using the values at i-1 and i
    frac = (r-Mref_V%radius_ref(i-1))/(Mref_V%radius_ref(i)-Mref_V%radius_ref(i-1))

    rho = Mref_V%density_ref(i-1) + frac * (Mref_V%density_ref(i)-Mref_V%density_ref(i-1))
    vpv = Mref_V%vpv_ref(i-1) + frac * (Mref_V%vpv_ref(i)-Mref_V%vpv_ref(i-1))
    vph = Mref_V%vph_ref(i-1) + frac * (Mref_V%vph_ref(i)-Mref_V%vph_ref(i-1))
    vsv = Mref_V%vsv_ref(i-1) + frac * (Mref_V%vsv_ref(i)-Mref_V%vsv_ref(i-1))
    vsh = Mref_V%vsh_ref(i-1) + frac * (Mref_V%vsh_ref(i)-Mref_V%vsh_ref(i-1))
    eta = Mref_V%eta_ref(i-1) + frac * (Mref_V%eta_ref(i)-Mref_V%eta_ref(i-1))
    Qkappa = Mref_V%Qkappa_ref(i-1) + frac * (Mref_V%Qkappa_ref(i)-Mref_V%Qkappa_ref(i-1))
    Qmu = Mref_V%Qmu_ref(i-1) + frac * (Mref_V%Qmu_ref(i)-Mref_V%Qmu_ref(i-1))

  endif

! make sure Vs is zero in the outer core even if roundoff errors on depth
! also set fictitious attenuation to a very high value (attenuation is not used in the fluid)
  if(iregion_code == IREGION_OUTER_CORE) then
    vsv = 0.d0
    vsh = 0.d0
    Qkappa = 3000.d0
    Qmu = 3000.d0
  endif

! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho/RHOAV
  vpv=vpv/(R_EARTH*scaleval)
  vph=vph/(R_EARTH*scaleval)
  vsv=vsv/(R_EARTH*scaleval)
  vsh=vsh/(R_EARTH*scaleval)

  end subroutine model_ref

!-------------------

  subroutine define_model_ref(USE_EXTERNAL_CRUSTAL_MODEL,Mref_V)

  implicit none
  include "constants.h"

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

  logical USE_EXTERNAL_CRUSTAL_MODEL

  integer i

! define the 1D REF model of Kustowski et al. (2007)

  Mref_V%radius_ref( 1 ) = 0.000000000000000E+000
  Mref_V%density_ref( 1 ) = 13088.4800000000
  Mref_V%vpv_ref( 1 ) = 11262.2000000000
  Mref_V%vsv_ref( 1 ) = 3667.80000000000
  Mref_V%Qkappa_ref( 1 ) = 1327.60000000000
  Mref_V%Qmu_ref( 1 ) = 104.000000000000
  Mref_V%vph_ref( 1 ) = 11262.2000000000
  Mref_V%vsh_ref( 1 ) = 3667.80000000000
  Mref_V%eta_ref( 1 ) = 1.00000000000000

  Mref_V%radius_ref( 2 ) = 6824.00000000000
  Mref_V%density_ref( 2 ) = 13088.4700000000
  Mref_V%vpv_ref( 2 ) = 11262.2000000000
  Mref_V%vsv_ref( 2 ) = 3667.79000000000
  Mref_V%Qkappa_ref( 2 ) = 1327.60000000000
  Mref_V%Qmu_ref( 2 ) = 104.000000000000
  Mref_V%vph_ref( 2 ) = 11262.2000000000
  Mref_V%vsh_ref( 2 ) = 3667.79000000000
  Mref_V%eta_ref( 2 ) = 1.00000000000000

  Mref_V%radius_ref( 3 ) = 13648.0000000000
  Mref_V%density_ref( 3 ) = 13088.4400000000
  Mref_V%vpv_ref( 3 ) = 11262.1800000000
  Mref_V%vsv_ref( 3 ) = 3667.78000000000
  Mref_V%Qkappa_ref( 3 ) = 1327.60000000000
  Mref_V%Qmu_ref( 3 ) = 104.000000000000
  Mref_V%vph_ref( 3 ) = 11262.1800000000
  Mref_V%vsh_ref( 3 ) = 3667.78000000000
  Mref_V%eta_ref( 3 ) = 1.00000000000000

  Mref_V%radius_ref( 4 ) = 20472.0000000000
  Mref_V%density_ref( 4 ) = 13088.3900000000
  Mref_V%vpv_ref( 4 ) = 11262.1400000000
  Mref_V%vsv_ref( 4 ) = 3667.75000000000
  Mref_V%Qkappa_ref( 4 ) = 1327.60000000000
  Mref_V%Qmu_ref( 4 ) = 104.000000000000
  Mref_V%vph_ref( 4 ) = 11262.1400000000
  Mref_V%vsh_ref( 4 ) = 3667.75000000000
  Mref_V%eta_ref( 4 ) = 1.00000000000000

  Mref_V%radius_ref( 5 ) = 27296.0000000000
  Mref_V%density_ref( 5 ) = 13088.3200000000
  Mref_V%vpv_ref( 5 ) = 11262.0900000000
  Mref_V%vsv_ref( 5 ) = 3667.72000000000
  Mref_V%Qkappa_ref( 5 ) = 1327.60000000000
  Mref_V%Qmu_ref( 5 ) = 104.000000000000
  Mref_V%vph_ref( 5 ) = 11262.0900000000
  Mref_V%vsh_ref( 5 ) = 3667.72000000000
  Mref_V%eta_ref( 5 ) = 1.00000000000000

  Mref_V%radius_ref( 6 ) = 34120.0000000000
  Mref_V%density_ref( 6 ) = 13088.2200000000
  Mref_V%vpv_ref( 6 ) = 11262.0200000000
  Mref_V%vsv_ref( 6 ) = 3667.67000000000
  Mref_V%Qkappa_ref( 6 ) = 1327.60000000000
  Mref_V%Qmu_ref( 6 ) = 104.000000000000
  Mref_V%vph_ref( 6 ) = 11262.0200000000
  Mref_V%vsh_ref( 6 ) = 3667.67000000000
  Mref_V%eta_ref( 6 ) = 1.00000000000000

  Mref_V%radius_ref( 7 ) = 40944.0000000000
  Mref_V%density_ref( 7 ) = 13088.1100000000
  Mref_V%vpv_ref( 7 ) = 11261.9400000000
  Mref_V%vsv_ref( 7 ) = 3667.62000000000
  Mref_V%Qkappa_ref( 7 ) = 1327.60000000000
  Mref_V%Qmu_ref( 7 ) = 104.000000000000
  Mref_V%vph_ref( 7 ) = 11261.9400000000
  Mref_V%vsh_ref( 7 ) = 3667.62000000000
  Mref_V%eta_ref( 7 ) = 1.00000000000000

  Mref_V%radius_ref( 8 ) = 47768.0000000000
  Mref_V%density_ref( 8 ) = 13087.9800000000
  Mref_V%vpv_ref( 8 ) = 11261.8500000000
  Mref_V%vsv_ref( 8 ) = 3667.55000000000
  Mref_V%Qkappa_ref( 8 ) = 1327.60000000000
  Mref_V%Qmu_ref( 8 ) = 104.000000000000
  Mref_V%vph_ref( 8 ) = 11261.8500000000
  Mref_V%vsh_ref( 8 ) = 3667.55000000000
  Mref_V%eta_ref( 8 ) = 1.00000000000000

  Mref_V%radius_ref( 9 ) = 54592.0000000000
  Mref_V%density_ref( 9 ) = 13087.8300000000
  Mref_V%vpv_ref( 9 ) = 11261.7400000000
  Mref_V%vsv_ref( 9 ) = 3667.47000000000
  Mref_V%Qkappa_ref( 9 ) = 1327.60000000000
  Mref_V%Qmu_ref( 9 ) = 104.000000000000
  Mref_V%vph_ref( 9 ) = 11261.7400000000
  Mref_V%vsh_ref( 9 ) = 3667.47000000000
  Mref_V%eta_ref( 9 ) = 1.00000000000000

  Mref_V%radius_ref( 10 ) = 61416.0000000000
  Mref_V%density_ref( 10 ) = 13087.6600000000
  Mref_V%vpv_ref( 10 ) = 11261.6100000000
  Mref_V%vsv_ref( 10 ) = 3667.39000000000
  Mref_V%Qkappa_ref( 10 ) = 1327.60000000000
  Mref_V%Qmu_ref( 10 ) = 104.000000000000
  Mref_V%vph_ref( 10 ) = 11261.6100000000
  Mref_V%vsh_ref( 10 ) = 3667.39000000000
  Mref_V%eta_ref( 10 ) = 1.00000000000000

  Mref_V%radius_ref( 11 ) = 68240.0000000000
  Mref_V%density_ref( 11 ) = 13087.4600000000
  Mref_V%vpv_ref( 11 ) = 11261.4700000000
  Mref_V%vsv_ref( 11 ) = 3667.29000000000
  Mref_V%Qkappa_ref( 11 ) = 1327.60000000000
  Mref_V%Qmu_ref( 11 ) = 104.000000000000
  Mref_V%vph_ref( 11 ) = 11261.4700000000
  Mref_V%vsh_ref( 11 ) = 3667.29000000000
  Mref_V%eta_ref( 11 ) = 1.00000000000000

  Mref_V%radius_ref( 12 ) = 75064.0000000000
  Mref_V%density_ref( 12 ) = 13087.2500000000
  Mref_V%vpv_ref( 12 ) = 11261.3200000000
  Mref_V%vsv_ref( 12 ) = 3667.18000000000
  Mref_V%Qkappa_ref( 12 ) = 1327.60000000000
  Mref_V%Qmu_ref( 12 ) = 104.000000000000
  Mref_V%vph_ref( 12 ) = 11261.3200000000
  Mref_V%vsh_ref( 12 ) = 3667.18000000000
  Mref_V%eta_ref( 12 ) = 1.00000000000000

  Mref_V%radius_ref( 13 ) = 81888.0000000000
  Mref_V%density_ref( 13 ) = 13087.0200000000
  Mref_V%vpv_ref( 13 ) = 11261.1500000000
  Mref_V%vsv_ref( 13 ) = 3667.06000000000
  Mref_V%Qkappa_ref( 13 ) = 1327.60000000000
  Mref_V%Qmu_ref( 13 ) = 104.000000000000
  Mref_V%vph_ref( 13 ) = 11261.1500000000
  Mref_V%vsh_ref( 13 ) = 3667.06000000000
  Mref_V%eta_ref( 13 ) = 1.00000000000000

  Mref_V%radius_ref( 14 ) = 88712.0000000000
  Mref_V%density_ref( 14 ) = 13086.7600000000
  Mref_V%vpv_ref( 14 ) = 11260.9700000000
  Mref_V%vsv_ref( 14 ) = 3666.94000000000
  Mref_V%Qkappa_ref( 14 ) = 1327.60000000000
  Mref_V%Qmu_ref( 14 ) = 104.000000000000
  Mref_V%vph_ref( 14 ) = 11260.9700000000
  Mref_V%vsh_ref( 14 ) = 3666.94000000000
  Mref_V%eta_ref( 14 ) = 1.00000000000000

  Mref_V%radius_ref( 15 ) = 95536.0000000000
  Mref_V%density_ref( 15 ) = 13086.4900000000
  Mref_V%vpv_ref( 15 ) = 11260.7700000000
  Mref_V%vsv_ref( 15 ) = 3666.80000000000
  Mref_V%Qkappa_ref( 15 ) = 1327.60000000000
  Mref_V%Qmu_ref( 15 ) = 104.000000000000
  Mref_V%vph_ref( 15 ) = 11260.7700000000
  Mref_V%vsh_ref( 15 ) = 3666.80000000000
  Mref_V%eta_ref( 15 ) = 1.00000000000000

  Mref_V%radius_ref( 16 ) = 102360.000000000
  Mref_V%density_ref( 16 ) = 13086.2000000000
  Mref_V%vpv_ref( 16 ) = 11260.5600000000
  Mref_V%vsv_ref( 16 ) = 3666.65000000000
  Mref_V%Qkappa_ref( 16 ) = 1327.60000000000
  Mref_V%Qmu_ref( 16 ) = 104.000000000000
  Mref_V%vph_ref( 16 ) = 11260.5600000000
  Mref_V%vsh_ref( 16 ) = 3666.65000000000
  Mref_V%eta_ref( 16 ) = 1.00000000000000

  Mref_V%radius_ref( 17 ) = 109184.000000000
  Mref_V%density_ref( 17 ) = 13085.8800000000
  Mref_V%vpv_ref( 17 ) = 11260.3400000000
  Mref_V%vsv_ref( 17 ) = 3666.49000000000
  Mref_V%Qkappa_ref( 17 ) = 1327.60000000000
  Mref_V%Qmu_ref( 17 ) = 104.000000000000
  Mref_V%vph_ref( 17 ) = 11260.3400000000
  Mref_V%vsh_ref( 17 ) = 3666.49000000000
  Mref_V%eta_ref( 17 ) = 1.00000000000000

  Mref_V%radius_ref( 18 ) = 116007.000000000
  Mref_V%density_ref( 18 ) = 13085.5500000000
  Mref_V%vpv_ref( 18 ) = 11260.0900000000
  Mref_V%vsv_ref( 18 ) = 3666.32000000000
  Mref_V%Qkappa_ref( 18 ) = 1327.60000000000
  Mref_V%Qmu_ref( 18 ) = 104.000000000000
  Mref_V%vph_ref( 18 ) = 11260.0900000000
  Mref_V%vsh_ref( 18 ) = 3666.32000000000
  Mref_V%eta_ref( 18 ) = 1.00000000000000

  Mref_V%radius_ref( 19 ) = 122831.000000000
  Mref_V%density_ref( 19 ) = 13085.1900000000
  Mref_V%vpv_ref( 19 ) = 11259.8400000000
  Mref_V%vsv_ref( 19 ) = 3666.15000000000
  Mref_V%Qkappa_ref( 19 ) = 1327.60000000000
  Mref_V%Qmu_ref( 19 ) = 104.000000000000
  Mref_V%vph_ref( 19 ) = 11259.8400000000
  Mref_V%vsh_ref( 19 ) = 3666.15000000000
  Mref_V%eta_ref( 19 ) = 1.00000000000000

  Mref_V%radius_ref( 20 ) = 129655.000000000
  Mref_V%density_ref( 20 ) = 13084.8200000000
  Mref_V%vpv_ref( 20 ) = 11259.5700000000
  Mref_V%vsv_ref( 20 ) = 3665.96000000000
  Mref_V%Qkappa_ref( 20 ) = 1327.60000000000
  Mref_V%Qmu_ref( 20 ) = 104.000000000000
  Mref_V%vph_ref( 20 ) = 11259.5700000000
  Mref_V%vsh_ref( 20 ) = 3665.96000000000
  Mref_V%eta_ref( 20 ) = 1.00000000000000

  Mref_V%radius_ref( 21 ) = 136479.000000000
  Mref_V%density_ref( 21 ) = 13084.4200000000
  Mref_V%vpv_ref( 21 ) = 11259.2800000000
  Mref_V%vsv_ref( 21 ) = 3665.76000000000
  Mref_V%Qkappa_ref( 21 ) = 1327.60000000000
  Mref_V%Qmu_ref( 21 ) = 104.000000000000
  Mref_V%vph_ref( 21 ) = 11259.2800000000
  Mref_V%vsh_ref( 21 ) = 3665.76000000000
  Mref_V%eta_ref( 21 ) = 1.00000000000000

  Mref_V%radius_ref( 22 ) = 143303.000000000
  Mref_V%density_ref( 22 ) = 13084.0100000000
  Mref_V%vpv_ref( 22 ) = 11258.9900000000
  Mref_V%vsv_ref( 22 ) = 3665.55000000000
  Mref_V%Qkappa_ref( 22 ) = 1327.60000000000
  Mref_V%Qmu_ref( 22 ) = 104.000000000000
  Mref_V%vph_ref( 22 ) = 11258.9900000000
  Mref_V%vsh_ref( 22 ) = 3665.55000000000
  Mref_V%eta_ref( 22 ) = 1.00000000000000

  Mref_V%radius_ref( 23 ) = 150127.000000000
  Mref_V%density_ref( 23 ) = 13083.5700000000
  Mref_V%vpv_ref( 23 ) = 11258.6700000000
  Mref_V%vsv_ref( 23 ) = 3665.33000000000
  Mref_V%Qkappa_ref( 23 ) = 1327.60000000000
  Mref_V%Qmu_ref( 23 ) = 104.000000000000
  Mref_V%vph_ref( 23 ) = 11258.6700000000
  Mref_V%vsh_ref( 23 ) = 3665.33000000000
  Mref_V%eta_ref( 23 ) = 1.00000000000000

  Mref_V%radius_ref( 24 ) = 156951.000000000
  Mref_V%density_ref( 24 ) = 13083.1100000000
  Mref_V%vpv_ref( 24 ) = 11258.3400000000
  Mref_V%vsv_ref( 24 ) = 3665.10000000000
  Mref_V%Qkappa_ref( 24 ) = 1327.60000000000
  Mref_V%Qmu_ref( 24 ) = 104.000000000000
  Mref_V%vph_ref( 24 ) = 11258.3400000000
  Mref_V%vsh_ref( 24 ) = 3665.10000000000
  Mref_V%eta_ref( 24 ) = 1.00000000000000

  Mref_V%radius_ref( 25 ) = 163775.000000000
  Mref_V%density_ref( 25 ) = 13082.6400000000
  Mref_V%vpv_ref( 25 ) = 11258.0000000000
  Mref_V%vsv_ref( 25 ) = 3664.86000000000
  Mref_V%Qkappa_ref( 25 ) = 1327.60000000000
  Mref_V%Qmu_ref( 25 ) = 104.000000000000
  Mref_V%vph_ref( 25 ) = 11258.0000000000
  Mref_V%vsh_ref( 25 ) = 3664.86000000000
  Mref_V%eta_ref( 25 ) = 1.00000000000000

  Mref_V%radius_ref( 26 ) = 170599.000000000
  Mref_V%density_ref( 26 ) = 13082.1400000000
  Mref_V%vpv_ref( 26 ) = 11257.6400000000
  Mref_V%vsv_ref( 26 ) = 3664.61000000000
  Mref_V%Qkappa_ref( 26 ) = 1327.60000000000
  Mref_V%Qmu_ref( 26 ) = 104.000000000000
  Mref_V%vph_ref( 26 ) = 11257.6400000000
  Mref_V%vsh_ref( 26 ) = 3664.61000000000
  Mref_V%eta_ref( 26 ) = 1.00000000000000

  Mref_V%radius_ref( 27 ) = 177423.000000000
  Mref_V%density_ref( 27 ) = 13081.6200000000
  Mref_V%vpv_ref( 27 ) = 11257.2700000000
  Mref_V%vsv_ref( 27 ) = 3664.35000000000
  Mref_V%Qkappa_ref( 27 ) = 1327.60000000000
  Mref_V%Qmu_ref( 27 ) = 104.000000000000
  Mref_V%vph_ref( 27 ) = 11257.2700000000
  Mref_V%vsh_ref( 27 ) = 3664.35000000000
  Mref_V%eta_ref( 27 ) = 1.00000000000000

  Mref_V%radius_ref( 28 ) = 184247.000000000
  Mref_V%density_ref( 28 ) = 13081.0900000000
  Mref_V%vpv_ref( 28 ) = 11256.8800000000
  Mref_V%vsv_ref( 28 ) = 3664.08000000000
  Mref_V%Qkappa_ref( 28 ) = 1327.60000000000
  Mref_V%Qmu_ref( 28 ) = 104.000000000000
  Mref_V%vph_ref( 28 ) = 11256.8800000000
  Mref_V%vsh_ref( 28 ) = 3664.08000000000
  Mref_V%eta_ref( 28 ) = 1.00000000000000

  Mref_V%radius_ref( 29 ) = 191071.000000000
  Mref_V%density_ref( 29 ) = 13080.5300000000
  Mref_V%vpv_ref( 29 ) = 11256.4800000000
  Mref_V%vsv_ref( 29 ) = 3663.80000000000
  Mref_V%Qkappa_ref( 29 ) = 1327.60000000000
  Mref_V%Qmu_ref( 29 ) = 104.000000000000
  Mref_V%vph_ref( 29 ) = 11256.4800000000
  Mref_V%vsh_ref( 29 ) = 3663.80000000000
  Mref_V%eta_ref( 29 ) = 1.00000000000000

  Mref_V%radius_ref( 30 ) = 197895.000000000
  Mref_V%density_ref( 30 ) = 13079.9500000000
  Mref_V%vpv_ref( 30 ) = 11256.0600000000
  Mref_V%vsv_ref( 30 ) = 3663.51000000000
  Mref_V%Qkappa_ref( 30 ) = 1327.60000000000
  Mref_V%Qmu_ref( 30 ) = 104.000000000000
  Mref_V%vph_ref( 30 ) = 11256.0600000000
  Mref_V%vsh_ref( 30 ) = 3663.51000000000
  Mref_V%eta_ref( 30 ) = 1.00000000000000

  Mref_V%radius_ref( 31 ) = 204719.000000000
  Mref_V%density_ref( 31 ) = 13079.3500000000
  Mref_V%vpv_ref( 31 ) = 11255.6300000000
  Mref_V%vsv_ref( 31 ) = 3663.21000000000
  Mref_V%Qkappa_ref( 31 ) = 1327.60000000000
  Mref_V%Qmu_ref( 31 ) = 104.000000000000
  Mref_V%vph_ref( 31 ) = 11255.6300000000
  Mref_V%vsh_ref( 31 ) = 3663.21000000000
  Mref_V%eta_ref( 31 ) = 1.00000000000000

  Mref_V%radius_ref( 32 ) = 211543.000000000
  Mref_V%density_ref( 32 ) = 13078.7300000000
  Mref_V%vpv_ref( 32 ) = 11255.1900000000
  Mref_V%vsv_ref( 32 ) = 3662.90000000000
  Mref_V%Qkappa_ref( 32 ) = 1327.60000000000
  Mref_V%Qmu_ref( 32 ) = 104.000000000000
  Mref_V%vph_ref( 32 ) = 11255.1900000000
  Mref_V%vsh_ref( 32 ) = 3662.90000000000
  Mref_V%eta_ref( 32 ) = 1.00000000000000

  Mref_V%radius_ref( 33 ) = 218367.000000000
  Mref_V%density_ref( 33 ) = 13078.0900000000
  Mref_V%vpv_ref( 33 ) = 11254.7300000000
  Mref_V%vsv_ref( 33 ) = 3662.57000000000
  Mref_V%Qkappa_ref( 33 ) = 1327.60000000000
  Mref_V%Qmu_ref( 33 ) = 104.000000000000
  Mref_V%vph_ref( 33 ) = 11254.7300000000
  Mref_V%vsh_ref( 33 ) = 3662.57000000000
  Mref_V%eta_ref( 33 ) = 1.00000000000000

  Mref_V%radius_ref( 34 ) = 225191.000000000
  Mref_V%density_ref( 34 ) = 13077.4400000000
  Mref_V%vpv_ref( 34 ) = 11254.2500000000
  Mref_V%vsv_ref( 34 ) = 3662.24000000000
  Mref_V%Qkappa_ref( 34 ) = 1327.60000000000
  Mref_V%Qmu_ref( 34 ) = 104.000000000000
  Mref_V%vph_ref( 34 ) = 11254.2500000000
  Mref_V%vsh_ref( 34 ) = 3662.24000000000
  Mref_V%eta_ref( 34 ) = 1.00000000000000

  Mref_V%radius_ref( 35 ) = 232015.000000000
  Mref_V%density_ref( 35 ) = 13076.7600000000
  Mref_V%vpv_ref( 35 ) = 11253.7600000000
  Mref_V%vsv_ref( 35 ) = 3661.90000000000
  Mref_V%Qkappa_ref( 35 ) = 1327.60000000000
  Mref_V%Qmu_ref( 35 ) = 104.000000000000
  Mref_V%vph_ref( 35 ) = 11253.7600000000
  Mref_V%vsh_ref( 35 ) = 3661.90000000000
  Mref_V%eta_ref( 35 ) = 1.00000000000000

  Mref_V%radius_ref( 36 ) = 238839.000000000
  Mref_V%density_ref( 36 ) = 13076.0600000000
  Mref_V%vpv_ref( 36 ) = 11253.2600000000
  Mref_V%vsv_ref( 36 ) = 3661.55000000000
  Mref_V%Qkappa_ref( 36 ) = 1327.60000000000
  Mref_V%Qmu_ref( 36 ) = 104.000000000000
  Mref_V%vph_ref( 36 ) = 11253.2600000000
  Mref_V%vsh_ref( 36 ) = 3661.55000000000
  Mref_V%eta_ref( 36 ) = 1.00000000000000

  Mref_V%radius_ref( 37 ) = 245663.000000000
  Mref_V%density_ref( 37 ) = 13075.3400000000
  Mref_V%vpv_ref( 37 ) = 11252.7400000000
  Mref_V%vsv_ref( 37 ) = 3661.19000000000
  Mref_V%Qkappa_ref( 37 ) = 1327.60000000000
  Mref_V%Qmu_ref( 37 ) = 104.000000000000
  Mref_V%vph_ref( 37 ) = 11252.7400000000
  Mref_V%vsh_ref( 37 ) = 3661.19000000000
  Mref_V%eta_ref( 37 ) = 1.00000000000000

  Mref_V%radius_ref( 38 ) = 252487.000000000
  Mref_V%density_ref( 38 ) = 13074.6000000000
  Mref_V%vpv_ref( 38 ) = 11252.2100000000
  Mref_V%vsv_ref( 38 ) = 3660.81000000000
  Mref_V%Qkappa_ref( 38 ) = 1327.60000000000
  Mref_V%Qmu_ref( 38 ) = 104.000000000000
  Mref_V%vph_ref( 38 ) = 11252.2100000000
  Mref_V%vsh_ref( 38 ) = 3660.81000000000
  Mref_V%eta_ref( 38 ) = 1.00000000000000

  Mref_V%radius_ref( 39 ) = 259311.000000000
  Mref_V%density_ref( 39 ) = 13073.8400000000
  Mref_V%vpv_ref( 39 ) = 11251.6600000000
  Mref_V%vsv_ref( 39 ) = 3660.43000000000
  Mref_V%Qkappa_ref( 39 ) = 1327.60000000000
  Mref_V%Qmu_ref( 39 ) = 104.000000000000
  Mref_V%vph_ref( 39 ) = 11251.6600000000
  Mref_V%vsh_ref( 39 ) = 3660.43000000000
  Mref_V%eta_ref( 39 ) = 1.00000000000000

  Mref_V%radius_ref( 40 ) = 266135.000000000
  Mref_V%density_ref( 40 ) = 13073.0600000000
  Mref_V%vpv_ref( 40 ) = 11251.1000000000
  Mref_V%vsv_ref( 40 ) = 3660.04000000000
  Mref_V%Qkappa_ref( 40 ) = 1327.60000000000
  Mref_V%Qmu_ref( 40 ) = 104.000000000000
  Mref_V%vph_ref( 40 ) = 11251.1000000000
  Mref_V%vsh_ref( 40 ) = 3660.04000000000
  Mref_V%eta_ref( 40 ) = 1.00000000000000

  Mref_V%radius_ref( 41 ) = 272959.000000000
  Mref_V%density_ref( 41 ) = 13072.2500000000
  Mref_V%vpv_ref( 41 ) = 11250.5200000000
  Mref_V%vsv_ref( 41 ) = 3659.64000000000
  Mref_V%Qkappa_ref( 41 ) = 1327.60000000000
  Mref_V%Qmu_ref( 41 ) = 104.000000000000
  Mref_V%vph_ref( 41 ) = 11250.5200000000
  Mref_V%vsh_ref( 41 ) = 3659.64000000000
  Mref_V%eta_ref( 41 ) = 1.00000000000000

  Mref_V%radius_ref( 42 ) = 279783.000000000
  Mref_V%density_ref( 42 ) = 13071.4300000000
  Mref_V%vpv_ref( 42 ) = 11249.9300000000
  Mref_V%vsv_ref( 42 ) = 3659.22000000000
  Mref_V%Qkappa_ref( 42 ) = 1327.60000000000
  Mref_V%Qmu_ref( 42 ) = 104.000000000000
  Mref_V%vph_ref( 42 ) = 11249.9300000000
  Mref_V%vsh_ref( 42 ) = 3659.22000000000
  Mref_V%eta_ref( 42 ) = 1.00000000000000

  Mref_V%radius_ref( 43 ) = 286607.000000000
  Mref_V%density_ref( 43 ) = 13070.5900000000
  Mref_V%vpv_ref( 43 ) = 11249.3300000000
  Mref_V%vsv_ref( 43 ) = 3658.80000000000
  Mref_V%Qkappa_ref( 43 ) = 1327.60000000000
  Mref_V%Qmu_ref( 43 ) = 104.000000000000
  Mref_V%vph_ref( 43 ) = 11249.3300000000
  Mref_V%vsh_ref( 43 ) = 3658.80000000000
  Mref_V%eta_ref( 43 ) = 1.00000000000000

  Mref_V%radius_ref( 44 ) = 293431.000000000
  Mref_V%density_ref( 44 ) = 13069.7300000000
  Mref_V%vpv_ref( 44 ) = 11248.7100000000
  Mref_V%vsv_ref( 44 ) = 3658.36000000000
  Mref_V%Qkappa_ref( 44 ) = 1327.60000000000
  Mref_V%Qmu_ref( 44 ) = 104.000000000000
  Mref_V%vph_ref( 44 ) = 11248.7100000000
  Mref_V%vsh_ref( 44 ) = 3658.36000000000
  Mref_V%eta_ref( 44 ) = 1.00000000000000

  Mref_V%radius_ref( 45 ) = 300255.000000000
  Mref_V%density_ref( 45 ) = 13068.8500000000
  Mref_V%vpv_ref( 45 ) = 11248.0700000000
  Mref_V%vsv_ref( 45 ) = 3657.92000000000
  Mref_V%Qkappa_ref( 45 ) = 1327.60000000000
  Mref_V%Qmu_ref( 45 ) = 104.000000000000
  Mref_V%vph_ref( 45 ) = 11248.0700000000
  Mref_V%vsh_ref( 45 ) = 3657.92000000000
  Mref_V%eta_ref( 45 ) = 1.00000000000000

  Mref_V%radius_ref( 46 ) = 307079.000000000
  Mref_V%density_ref( 46 ) = 13067.9500000000
  Mref_V%vpv_ref( 46 ) = 11247.4200000000
  Mref_V%vsv_ref( 46 ) = 3657.47000000000
  Mref_V%Qkappa_ref( 46 ) = 1327.60000000000
  Mref_V%Qmu_ref( 46 ) = 104.000000000000
  Mref_V%vph_ref( 46 ) = 11247.4200000000
  Mref_V%vsh_ref( 46 ) = 3657.47000000000
  Mref_V%eta_ref( 46 ) = 1.00000000000000

  Mref_V%radius_ref( 47 ) = 313903.000000000
  Mref_V%density_ref( 47 ) = 13067.0200000000
  Mref_V%vpv_ref( 47 ) = 11246.7600000000
  Mref_V%vsv_ref( 47 ) = 3657.00000000000
  Mref_V%Qkappa_ref( 47 ) = 1327.60000000000
  Mref_V%Qmu_ref( 47 ) = 104.000000000000
  Mref_V%vph_ref( 47 ) = 11246.7600000000
  Mref_V%vsh_ref( 47 ) = 3657.00000000000
  Mref_V%eta_ref( 47 ) = 1.00000000000000

  Mref_V%radius_ref( 48 ) = 320727.000000000
  Mref_V%density_ref( 48 ) = 13066.0800000000
  Mref_V%vpv_ref( 48 ) = 11246.0800000000
  Mref_V%vsv_ref( 48 ) = 3656.53000000000
  Mref_V%Qkappa_ref( 48 ) = 1327.60000000000
  Mref_V%Qmu_ref( 48 ) = 104.000000000000
  Mref_V%vph_ref( 48 ) = 11246.0800000000
  Mref_V%vsh_ref( 48 ) = 3656.53000000000
  Mref_V%eta_ref( 48 ) = 1.00000000000000

  Mref_V%radius_ref( 49 ) = 327551.000000000
  Mref_V%density_ref( 49 ) = 13065.1200000000
  Mref_V%vpv_ref( 49 ) = 11245.3800000000
  Mref_V%vsv_ref( 49 ) = 3656.04000000000
  Mref_V%Qkappa_ref( 49 ) = 1327.60000000000
  Mref_V%Qmu_ref( 49 ) = 104.000000000000
  Mref_V%vph_ref( 49 ) = 11245.3800000000
  Mref_V%vsh_ref( 49 ) = 3656.04000000000
  Mref_V%eta_ref( 49 ) = 1.00000000000000

  Mref_V%radius_ref( 50 ) = 334375.000000000
  Mref_V%density_ref( 50 ) = 13064.1300000000
  Mref_V%vpv_ref( 50 ) = 11244.6700000000
  Mref_V%vsv_ref( 50 ) = 3655.55000000000
  Mref_V%Qkappa_ref( 50 ) = 1327.60000000000
  Mref_V%Qmu_ref( 50 ) = 104.000000000000
  Mref_V%vph_ref( 50 ) = 11244.6700000000
  Mref_V%vsh_ref( 50 ) = 3655.55000000000
  Mref_V%eta_ref( 50 ) = 1.00000000000000

  Mref_V%radius_ref( 51 ) = 341199.000000000
  Mref_V%density_ref( 51 ) = 13063.1300000000
  Mref_V%vpv_ref( 51 ) = 11243.9500000000
  Mref_V%vsv_ref( 51 ) = 3655.04000000000
  Mref_V%Qkappa_ref( 51 ) = 1327.60000000000
  Mref_V%Qmu_ref( 51 ) = 104.000000000000
  Mref_V%vph_ref( 51 ) = 11243.9500000000
  Mref_V%vsh_ref( 51 ) = 3655.04000000000
  Mref_V%eta_ref( 51 ) = 1.00000000000000

  Mref_V%radius_ref( 52 ) = 348022.000000000
  Mref_V%density_ref( 52 ) = 13062.1000000000
  Mref_V%vpv_ref( 52 ) = 11243.2100000000
  Mref_V%vsv_ref( 52 ) = 3654.53000000000
  Mref_V%Qkappa_ref( 52 ) = 1327.60000000000
  Mref_V%Qmu_ref( 52 ) = 104.000000000000
  Mref_V%vph_ref( 52 ) = 11243.2100000000
  Mref_V%vsh_ref( 52 ) = 3654.53000000000
  Mref_V%eta_ref( 52 ) = 1.00000000000000

  Mref_V%radius_ref( 53 ) = 354846.000000000
  Mref_V%density_ref( 53 ) = 13061.0600000000
  Mref_V%vpv_ref( 53 ) = 11242.4600000000
  Mref_V%vsv_ref( 53 ) = 3654.00000000000
  Mref_V%Qkappa_ref( 53 ) = 1327.60000000000
  Mref_V%Qmu_ref( 53 ) = 104.000000000000
  Mref_V%vph_ref( 53 ) = 11242.4600000000
  Mref_V%vsh_ref( 53 ) = 3654.00000000000
  Mref_V%eta_ref( 53 ) = 1.00000000000000

  Mref_V%radius_ref( 54 ) = 361670.000000000
  Mref_V%density_ref( 54 ) = 13060.0000000000
  Mref_V%vpv_ref( 54 ) = 11241.7000000000
  Mref_V%vsv_ref( 54 ) = 3653.47000000000
  Mref_V%Qkappa_ref( 54 ) = 1327.60000000000
  Mref_V%Qmu_ref( 54 ) = 104.000000000000
  Mref_V%vph_ref( 54 ) = 11241.7000000000
  Mref_V%vsh_ref( 54 ) = 3653.47000000000
  Mref_V%eta_ref( 54 ) = 1.00000000000000

  Mref_V%radius_ref( 55 ) = 368494.000000000
  Mref_V%density_ref( 55 ) = 13058.9100000000
  Mref_V%vpv_ref( 55 ) = 11240.9100000000
  Mref_V%vsv_ref( 55 ) = 3652.92000000000
  Mref_V%Qkappa_ref( 55 ) = 1327.60000000000
  Mref_V%Qmu_ref( 55 ) = 104.000000000000
  Mref_V%vph_ref( 55 ) = 11240.9100000000
  Mref_V%vsh_ref( 55 ) = 3652.92000000000
  Mref_V%eta_ref( 55 ) = 1.00000000000000

  Mref_V%radius_ref( 56 ) = 375318.000000000
  Mref_V%density_ref( 56 ) = 13057.8100000000
  Mref_V%vpv_ref( 56 ) = 11240.1200000000
  Mref_V%vsv_ref( 56 ) = 3652.36000000000
  Mref_V%Qkappa_ref( 56 ) = 1327.60000000000
  Mref_V%Qmu_ref( 56 ) = 104.000000000000
  Mref_V%vph_ref( 56 ) = 11240.1200000000
  Mref_V%vsh_ref( 56 ) = 3652.36000000000
  Mref_V%eta_ref( 56 ) = 1.00000000000000

  Mref_V%radius_ref( 57 ) = 382142.000000000
  Mref_V%density_ref( 57 ) = 13056.6800000000
  Mref_V%vpv_ref( 57 ) = 11239.3100000000
  Mref_V%vsv_ref( 57 ) = 3651.80000000000
  Mref_V%Qkappa_ref( 57 ) = 1327.60000000000
  Mref_V%Qmu_ref( 57 ) = 104.000000000000
  Mref_V%vph_ref( 57 ) = 11239.3100000000
  Mref_V%vsh_ref( 57 ) = 3651.80000000000
  Mref_V%eta_ref( 57 ) = 1.00000000000000

  Mref_V%radius_ref( 58 ) = 388966.000000000
  Mref_V%density_ref( 58 ) = 13055.5300000000
  Mref_V%vpv_ref( 58 ) = 11238.4800000000
  Mref_V%vsv_ref( 58 ) = 3651.22000000000
  Mref_V%Qkappa_ref( 58 ) = 1327.60000000000
  Mref_V%Qmu_ref( 58 ) = 104.000000000000
  Mref_V%vph_ref( 58 ) = 11238.4800000000
  Mref_V%vsh_ref( 58 ) = 3651.22000000000
  Mref_V%eta_ref( 58 ) = 1.00000000000000

  Mref_V%radius_ref( 59 ) = 395790.000000000
  Mref_V%density_ref( 59 ) = 13054.3700000000
  Mref_V%vpv_ref( 59 ) = 11237.6400000000
  Mref_V%vsv_ref( 59 ) = 3650.63000000000
  Mref_V%Qkappa_ref( 59 ) = 1327.60000000000
  Mref_V%Qmu_ref( 59 ) = 104.000000000000
  Mref_V%vph_ref( 59 ) = 11237.6400000000
  Mref_V%vsh_ref( 59 ) = 3650.63000000000
  Mref_V%eta_ref( 59 ) = 1.00000000000000

  Mref_V%radius_ref( 60 ) = 402614.000000000
  Mref_V%density_ref( 60 ) = 13053.1800000000
  Mref_V%vpv_ref( 60 ) = 11236.7900000000
  Mref_V%vsv_ref( 60 ) = 3650.04000000000
  Mref_V%Qkappa_ref( 60 ) = 1327.60000000000
  Mref_V%Qmu_ref( 60 ) = 104.000000000000
  Mref_V%vph_ref( 60 ) = 11236.7900000000
  Mref_V%vsh_ref( 60 ) = 3650.04000000000
  Mref_V%eta_ref( 60 ) = 1.00000000000000

  Mref_V%radius_ref( 61 ) = 409438.000000000
  Mref_V%density_ref( 61 ) = 13051.9800000000
  Mref_V%vpv_ref( 61 ) = 11235.9200000000
  Mref_V%vsv_ref( 61 ) = 3649.43000000000
  Mref_V%Qkappa_ref( 61 ) = 1327.60000000000
  Mref_V%Qmu_ref( 61 ) = 104.000000000000
  Mref_V%vph_ref( 61 ) = 11235.9200000000
  Mref_V%vsh_ref( 61 ) = 3649.43000000000
  Mref_V%eta_ref( 61 ) = 1.00000000000000

  Mref_V%radius_ref( 62 ) = 416262.000000000
  Mref_V%density_ref( 62 ) = 13050.7500000000
  Mref_V%vpv_ref( 62 ) = 11235.0400000000
  Mref_V%vsv_ref( 62 ) = 3648.81000000000
  Mref_V%Qkappa_ref( 62 ) = 1327.60000000000
  Mref_V%Qmu_ref( 62 ) = 104.000000000000
  Mref_V%vph_ref( 62 ) = 11235.0400000000
  Mref_V%vsh_ref( 62 ) = 3648.81000000000
  Mref_V%eta_ref( 62 ) = 1.00000000000000

  Mref_V%radius_ref( 63 ) = 423086.000000000
  Mref_V%density_ref( 63 ) = 13049.5000000000
  Mref_V%vpv_ref( 63 ) = 11234.1400000000
  Mref_V%vsv_ref( 63 ) = 3648.19000000000
  Mref_V%Qkappa_ref( 63 ) = 1327.60000000000
  Mref_V%Qmu_ref( 63 ) = 104.000000000000
  Mref_V%vph_ref( 63 ) = 11234.1400000000
  Mref_V%vsh_ref( 63 ) = 3648.19000000000
  Mref_V%eta_ref( 63 ) = 1.00000000000000

  Mref_V%radius_ref( 64 ) = 429910.000000000
  Mref_V%density_ref( 64 ) = 13048.2300000000
  Mref_V%vpv_ref( 64 ) = 11233.2300000000
  Mref_V%vsv_ref( 64 ) = 3647.55000000000
  Mref_V%Qkappa_ref( 64 ) = 1327.60000000000
  Mref_V%Qmu_ref( 64 ) = 104.000000000000
  Mref_V%vph_ref( 64 ) = 11233.2300000000
  Mref_V%vsh_ref( 64 ) = 3647.55000000000
  Mref_V%eta_ref( 64 ) = 1.00000000000000

  Mref_V%radius_ref( 65 ) = 436734.000000000
  Mref_V%density_ref( 65 ) = 13046.9500000000
  Mref_V%vpv_ref( 65 ) = 11232.3000000000
  Mref_V%vsv_ref( 65 ) = 3646.90000000000
  Mref_V%Qkappa_ref( 65 ) = 1327.60000000000
  Mref_V%Qmu_ref( 65 ) = 104.000000000000
  Mref_V%vph_ref( 65 ) = 11232.3000000000
  Mref_V%vsh_ref( 65 ) = 3646.90000000000
  Mref_V%eta_ref( 65 ) = 1.00000000000000

  Mref_V%radius_ref( 66 ) = 443558.000000000
  Mref_V%density_ref( 66 ) = 13045.6400000000
  Mref_V%vpv_ref( 66 ) = 11231.3600000000
  Mref_V%vsv_ref( 66 ) = 3646.24000000000
  Mref_V%Qkappa_ref( 66 ) = 1327.60000000000
  Mref_V%Qmu_ref( 66 ) = 104.000000000000
  Mref_V%vph_ref( 66 ) = 11231.3600000000
  Mref_V%vsh_ref( 66 ) = 3646.24000000000
  Mref_V%eta_ref( 66 ) = 1.00000000000000

  Mref_V%radius_ref( 67 ) = 450382.000000000
  Mref_V%density_ref( 67 ) = 13044.3100000000
  Mref_V%vpv_ref( 67 ) = 11230.4000000000
  Mref_V%vsv_ref( 67 ) = 3645.57000000000
  Mref_V%Qkappa_ref( 67 ) = 1327.60000000000
  Mref_V%Qmu_ref( 67 ) = 104.000000000000
  Mref_V%vph_ref( 67 ) = 11230.4000000000
  Mref_V%vsh_ref( 67 ) = 3645.57000000000
  Mref_V%eta_ref( 67 ) = 1.00000000000000

  Mref_V%radius_ref( 68 ) = 457206.000000000
  Mref_V%density_ref( 68 ) = 13042.9600000000
  Mref_V%vpv_ref( 68 ) = 11229.4300000000
  Mref_V%vsv_ref( 68 ) = 3644.89000000000
  Mref_V%Qkappa_ref( 68 ) = 1327.60000000000
  Mref_V%Qmu_ref( 68 ) = 104.000000000000
  Mref_V%vph_ref( 68 ) = 11229.4300000000
  Mref_V%vsh_ref( 68 ) = 3644.89000000000
  Mref_V%eta_ref( 68 ) = 1.00000000000000

  Mref_V%radius_ref( 69 ) = 464030.000000000
  Mref_V%density_ref( 69 ) = 13041.5900000000
  Mref_V%vpv_ref( 69 ) = 11228.4400000000
  Mref_V%vsv_ref( 69 ) = 3644.21000000000
  Mref_V%Qkappa_ref( 69 ) = 1327.60000000000
  Mref_V%Qmu_ref( 69 ) = 104.000000000000
  Mref_V%vph_ref( 69 ) = 11228.4400000000
  Mref_V%vsh_ref( 69 ) = 3644.21000000000
  Mref_V%eta_ref( 69 ) = 1.00000000000000

  Mref_V%radius_ref( 70 ) = 470854.000000000
  Mref_V%density_ref( 70 ) = 13040.2000000000
  Mref_V%vpv_ref( 70 ) = 11227.4400000000
  Mref_V%vsv_ref( 70 ) = 3643.51000000000
  Mref_V%Qkappa_ref( 70 ) = 1327.60000000000
  Mref_V%Qmu_ref( 70 ) = 104.000000000000
  Mref_V%vph_ref( 70 ) = 11227.4400000000
  Mref_V%vsh_ref( 70 ) = 3643.51000000000
  Mref_V%eta_ref( 70 ) = 1.00000000000000

  Mref_V%radius_ref( 71 ) = 477678.000000000
  Mref_V%density_ref( 71 ) = 13038.7900000000
  Mref_V%vpv_ref( 71 ) = 11226.4300000000
  Mref_V%vsv_ref( 71 ) = 3642.80000000000
  Mref_V%Qkappa_ref( 71 ) = 1327.60000000000
  Mref_V%Qmu_ref( 71 ) = 104.000000000000
  Mref_V%vph_ref( 71 ) = 11226.4300000000
  Mref_V%vsh_ref( 71 ) = 3642.80000000000
  Mref_V%eta_ref( 71 ) = 1.00000000000000

  Mref_V%radius_ref( 72 ) = 484502.000000000
  Mref_V%density_ref( 72 ) = 13037.3600000000
  Mref_V%vpv_ref( 72 ) = 11225.4000000000
  Mref_V%vsv_ref( 72 ) = 3642.08000000000
  Mref_V%Qkappa_ref( 72 ) = 1327.60000000000
  Mref_V%Qmu_ref( 72 ) = 104.000000000000
  Mref_V%vph_ref( 72 ) = 11225.4000000000
  Mref_V%vsh_ref( 72 ) = 3642.08000000000
  Mref_V%eta_ref( 72 ) = 1.00000000000000

  Mref_V%radius_ref( 73 ) = 491326.000000000
  Mref_V%density_ref( 73 ) = 13035.9100000000
  Mref_V%vpv_ref( 73 ) = 11224.3600000000
  Mref_V%vsv_ref( 73 ) = 3641.35000000000
  Mref_V%Qkappa_ref( 73 ) = 1327.60000000000
  Mref_V%Qmu_ref( 73 ) = 104.000000000000
  Mref_V%vph_ref( 73 ) = 11224.3600000000
  Mref_V%vsh_ref( 73 ) = 3641.35000000000
  Mref_V%eta_ref( 73 ) = 1.00000000000000

  Mref_V%radius_ref( 74 ) = 498150.000000000
  Mref_V%density_ref( 74 ) = 13034.4400000000
  Mref_V%vpv_ref( 74 ) = 11223.3000000000
  Mref_V%vsv_ref( 74 ) = 3640.61000000000
  Mref_V%Qkappa_ref( 74 ) = 1327.60000000000
  Mref_V%Qmu_ref( 74 ) = 104.000000000000
  Mref_V%vph_ref( 74 ) = 11223.3000000000
  Mref_V%vsh_ref( 74 ) = 3640.61000000000
  Mref_V%eta_ref( 74 ) = 1.00000000000000

  Mref_V%radius_ref( 75 ) = 504974.000000000
  Mref_V%density_ref( 75 ) = 13032.9500000000
  Mref_V%vpv_ref( 75 ) = 11222.2200000000
  Mref_V%vsv_ref( 75 ) = 3639.86000000000
  Mref_V%Qkappa_ref( 75 ) = 1327.60000000000
  Mref_V%Qmu_ref( 75 ) = 104.000000000000
  Mref_V%vph_ref( 75 ) = 11222.2200000000
  Mref_V%vsh_ref( 75 ) = 3639.86000000000
  Mref_V%eta_ref( 75 ) = 1.00000000000000

  Mref_V%radius_ref( 76 ) = 511798.000000000
  Mref_V%density_ref( 76 ) = 13031.4400000000
  Mref_V%vpv_ref( 76 ) = 11221.1400000000
  Mref_V%vsv_ref( 76 ) = 3639.10000000000
  Mref_V%Qkappa_ref( 76 ) = 1327.60000000000
  Mref_V%Qmu_ref( 76 ) = 104.000000000000
  Mref_V%vph_ref( 76 ) = 11221.1400000000
  Mref_V%vsh_ref( 76 ) = 3639.10000000000
  Mref_V%eta_ref( 76 ) = 1.00000000000000

  Mref_V%radius_ref( 77 ) = 518622.000000000
  Mref_V%density_ref( 77 ) = 13029.9100000000
  Mref_V%vpv_ref( 77 ) = 11220.0300000000
  Mref_V%vsv_ref( 77 ) = 3638.33000000000
  Mref_V%Qkappa_ref( 77 ) = 1327.60000000000
  Mref_V%Qmu_ref( 77 ) = 104.000000000000
  Mref_V%vph_ref( 77 ) = 11220.0300000000
  Mref_V%vsh_ref( 77 ) = 3638.33000000000
  Mref_V%eta_ref( 77 ) = 1.00000000000000

  Mref_V%radius_ref( 78 ) = 525446.000000000
  Mref_V%density_ref( 78 ) = 13028.3600000000
  Mref_V%vpv_ref( 78 ) = 11218.9200000000
  Mref_V%vsv_ref( 78 ) = 3637.55000000000
  Mref_V%Qkappa_ref( 78 ) = 1327.60000000000
  Mref_V%Qmu_ref( 78 ) = 104.000000000000
  Mref_V%vph_ref( 78 ) = 11218.9200000000
  Mref_V%vsh_ref( 78 ) = 3637.55000000000
  Mref_V%eta_ref( 78 ) = 1.00000000000000

  Mref_V%radius_ref( 79 ) = 532270.000000000
  Mref_V%density_ref( 79 ) = 13026.7900000000
  Mref_V%vpv_ref( 79 ) = 11217.7800000000
  Mref_V%vsv_ref( 79 ) = 3636.76000000000
  Mref_V%Qkappa_ref( 79 ) = 1327.60000000000
  Mref_V%Qmu_ref( 79 ) = 104.000000000000
  Mref_V%vph_ref( 79 ) = 11217.7800000000
  Mref_V%vsh_ref( 79 ) = 3636.76000000000
  Mref_V%eta_ref( 79 ) = 1.00000000000000

  Mref_V%radius_ref( 80 ) = 539094.000000000
  Mref_V%density_ref( 80 ) = 13025.2000000000
  Mref_V%vpv_ref( 80 ) = 11216.6400000000
  Mref_V%vsv_ref( 80 ) = 3635.96000000000
  Mref_V%Qkappa_ref( 80 ) = 1327.60000000000
  Mref_V%Qmu_ref( 80 ) = 104.000000000000
  Mref_V%vph_ref( 80 ) = 11216.6400000000
  Mref_V%vsh_ref( 80 ) = 3635.96000000000
  Mref_V%eta_ref( 80 ) = 1.00000000000000

  Mref_V%radius_ref( 81 ) = 545918.000000000
  Mref_V%density_ref( 81 ) = 13023.5800000000
  Mref_V%vpv_ref( 81 ) = 11215.4800000000
  Mref_V%vsv_ref( 81 ) = 3635.14000000000
  Mref_V%Qkappa_ref( 81 ) = 1327.60000000000
  Mref_V%Qmu_ref( 81 ) = 104.000000000000
  Mref_V%vph_ref( 81 ) = 11215.4800000000
  Mref_V%vsh_ref( 81 ) = 3635.14000000000
  Mref_V%eta_ref( 81 ) = 1.00000000000000

  Mref_V%radius_ref( 82 ) = 552742.000000000
  Mref_V%density_ref( 82 ) = 13021.9500000000
  Mref_V%vpv_ref( 82 ) = 11214.3000000000
  Mref_V%vsv_ref( 82 ) = 3634.32000000000
  Mref_V%Qkappa_ref( 82 ) = 1327.60000000000
  Mref_V%Qmu_ref( 82 ) = 104.000000000000
  Mref_V%vph_ref( 82 ) = 11214.3000000000
  Mref_V%vsh_ref( 82 ) = 3634.32000000000
  Mref_V%eta_ref( 82 ) = 1.00000000000000

  Mref_V%radius_ref( 83 ) = 559566.000000000
  Mref_V%density_ref( 83 ) = 13020.3000000000
  Mref_V%vpv_ref( 83 ) = 11213.1100000000
  Mref_V%vsv_ref( 83 ) = 3633.49000000000
  Mref_V%Qkappa_ref( 83 ) = 1327.60000000000
  Mref_V%Qmu_ref( 83 ) = 104.000000000000
  Mref_V%vph_ref( 83 ) = 11213.1100000000
  Mref_V%vsh_ref( 83 ) = 3633.49000000000
  Mref_V%eta_ref( 83 ) = 1.00000000000000

  Mref_V%radius_ref( 84 ) = 566390.000000000
  Mref_V%density_ref( 84 ) = 13018.6300000000
  Mref_V%vpv_ref( 84 ) = 11211.9100000000
  Mref_V%vsv_ref( 84 ) = 3632.65000000000
  Mref_V%Qkappa_ref( 84 ) = 1327.60000000000
  Mref_V%Qmu_ref( 84 ) = 104.000000000000
  Mref_V%vph_ref( 84 ) = 11211.9100000000
  Mref_V%vsh_ref( 84 ) = 3632.65000000000
  Mref_V%eta_ref( 84 ) = 1.00000000000000

  Mref_V%radius_ref( 85 ) = 573214.000000000
  Mref_V%density_ref( 85 ) = 13016.9300000000
  Mref_V%vpv_ref( 85 ) = 11210.6900000000
  Mref_V%vsv_ref( 85 ) = 3631.80000000000
  Mref_V%Qkappa_ref( 85 ) = 1327.60000000000
  Mref_V%Qmu_ref( 85 ) = 104.000000000000
  Mref_V%vph_ref( 85 ) = 11210.6900000000
  Mref_V%vsh_ref( 85 ) = 3631.80000000000
  Mref_V%eta_ref( 85 ) = 1.00000000000000

  Mref_V%radius_ref( 86 ) = 580037.000000000
  Mref_V%density_ref( 86 ) = 13015.2200000000
  Mref_V%vpv_ref( 86 ) = 11209.4500000000
  Mref_V%vsv_ref( 86 ) = 3630.93000000000
  Mref_V%Qkappa_ref( 86 ) = 1327.60000000000
  Mref_V%Qmu_ref( 86 ) = 104.000000000000
  Mref_V%vph_ref( 86 ) = 11209.4500000000
  Mref_V%vsh_ref( 86 ) = 3630.93000000000
  Mref_V%eta_ref( 86 ) = 1.00000000000000

  Mref_V%radius_ref( 87 ) = 586861.000000000
  Mref_V%density_ref( 87 ) = 13013.4900000000
  Mref_V%vpv_ref( 87 ) = 11208.2100000000
  Mref_V%vsv_ref( 87 ) = 3630.06000000000
  Mref_V%Qkappa_ref( 87 ) = 1327.60000000000
  Mref_V%Qmu_ref( 87 ) = 104.000000000000
  Mref_V%vph_ref( 87 ) = 11208.2100000000
  Mref_V%vsh_ref( 87 ) = 3630.06000000000
  Mref_V%eta_ref( 87 ) = 1.00000000000000

  Mref_V%radius_ref( 88 ) = 593685.000000000
  Mref_V%density_ref( 88 ) = 13011.7300000000
  Mref_V%vpv_ref( 88 ) = 11206.9400000000
  Mref_V%vsv_ref( 88 ) = 3629.18000000000
  Mref_V%Qkappa_ref( 88 ) = 1327.60000000000
  Mref_V%Qmu_ref( 88 ) = 104.000000000000
  Mref_V%vph_ref( 88 ) = 11206.9400000000
  Mref_V%vsh_ref( 88 ) = 3629.18000000000
  Mref_V%eta_ref( 88 ) = 1.00000000000000

  Mref_V%radius_ref( 89 ) = 600509.000000000
  Mref_V%density_ref( 89 ) = 13009.9600000000
  Mref_V%vpv_ref( 89 ) = 11205.6700000000
  Mref_V%vsv_ref( 89 ) = 3628.29000000000
  Mref_V%Qkappa_ref( 89 ) = 1327.60000000000
  Mref_V%Qmu_ref( 89 ) = 104.000000000000
  Mref_V%vph_ref( 89 ) = 11205.6700000000
  Mref_V%vsh_ref( 89 ) = 3628.29000000000
  Mref_V%eta_ref( 89 ) = 1.00000000000000

  Mref_V%radius_ref( 90 ) = 607333.000000000
  Mref_V%density_ref( 90 ) = 13008.1600000000
  Mref_V%vpv_ref( 90 ) = 11204.3700000000
  Mref_V%vsv_ref( 90 ) = 3627.38000000000
  Mref_V%Qkappa_ref( 90 ) = 1327.60000000000
  Mref_V%Qmu_ref( 90 ) = 104.000000000000
  Mref_V%vph_ref( 90 ) = 11204.3700000000
  Mref_V%vsh_ref( 90 ) = 3627.38000000000
  Mref_V%eta_ref( 90 ) = 1.00000000000000

  Mref_V%radius_ref( 91 ) = 614157.000000000
  Mref_V%density_ref( 91 ) = 13006.3500000000
  Mref_V%vpv_ref( 91 ) = 11203.0700000000
  Mref_V%vsv_ref( 91 ) = 3626.47000000000
  Mref_V%Qkappa_ref( 91 ) = 1327.60000000000
  Mref_V%Qmu_ref( 91 ) = 104.000000000000
  Mref_V%vph_ref( 91 ) = 11203.0700000000
  Mref_V%vsh_ref( 91 ) = 3626.47000000000
  Mref_V%eta_ref( 91 ) = 1.00000000000000

  Mref_V%radius_ref( 92 ) = 620981.000000000
  Mref_V%density_ref( 92 ) = 13004.5100000000
  Mref_V%vpv_ref( 92 ) = 11201.7400000000
  Mref_V%vsv_ref( 92 ) = 3625.55000000000
  Mref_V%Qkappa_ref( 92 ) = 1327.60000000000
  Mref_V%Qmu_ref( 92 ) = 104.000000000000
  Mref_V%vph_ref( 92 ) = 11201.7400000000
  Mref_V%vsh_ref( 92 ) = 3625.55000000000
  Mref_V%eta_ref( 92 ) = 1.00000000000000

  Mref_V%radius_ref( 93 ) = 627805.000000000
  Mref_V%density_ref( 93 ) = 13002.6600000000
  Mref_V%vpv_ref( 93 ) = 11200.4100000000
  Mref_V%vsv_ref( 93 ) = 3624.61000000000
  Mref_V%Qkappa_ref( 93 ) = 1327.60000000000
  Mref_V%Qmu_ref( 93 ) = 104.000000000000
  Mref_V%vph_ref( 93 ) = 11200.4100000000
  Mref_V%vsh_ref( 93 ) = 3624.61000000000
  Mref_V%eta_ref( 93 ) = 1.00000000000000

  Mref_V%radius_ref( 94 ) = 634629.000000000
  Mref_V%density_ref( 94 ) = 13000.7800000000
  Mref_V%vpv_ref( 94 ) = 11199.0600000000
  Mref_V%vsv_ref( 94 ) = 3623.67000000000
  Mref_V%Qkappa_ref( 94 ) = 1327.60000000000
  Mref_V%Qmu_ref( 94 ) = 104.000000000000
  Mref_V%vph_ref( 94 ) = 11199.0600000000
  Mref_V%vsh_ref( 94 ) = 3623.67000000000
  Mref_V%eta_ref( 94 ) = 1.00000000000000

  Mref_V%radius_ref( 95 ) = 641453.000000000
  Mref_V%density_ref( 95 ) = 12998.8800000000
  Mref_V%vpv_ref( 95 ) = 11197.6900000000
  Mref_V%vsv_ref( 95 ) = 3622.71000000000
  Mref_V%Qkappa_ref( 95 ) = 1327.60000000000
  Mref_V%Qmu_ref( 95 ) = 104.000000000000
  Mref_V%vph_ref( 95 ) = 11197.6900000000
  Mref_V%vsh_ref( 95 ) = 3622.71000000000
  Mref_V%eta_ref( 95 ) = 1.00000000000000

  Mref_V%radius_ref( 96 ) = 648277.000000000
  Mref_V%density_ref( 96 ) = 12996.9700000000
  Mref_V%vpv_ref( 96 ) = 11196.3100000000
  Mref_V%vsv_ref( 96 ) = 3621.75000000000
  Mref_V%Qkappa_ref( 96 ) = 1327.60000000000
  Mref_V%Qmu_ref( 96 ) = 104.000000000000
  Mref_V%vph_ref( 96 ) = 11196.3100000000
  Mref_V%vsh_ref( 96 ) = 3621.75000000000
  Mref_V%eta_ref( 96 ) = 1.00000000000000

  Mref_V%radius_ref( 97 ) = 655101.000000000
  Mref_V%density_ref( 97 ) = 12995.0300000000
  Mref_V%vpv_ref( 97 ) = 11194.9200000000
  Mref_V%vsv_ref( 97 ) = 3620.78000000000
  Mref_V%Qkappa_ref( 97 ) = 1327.60000000000
  Mref_V%Qmu_ref( 97 ) = 104.000000000000
  Mref_V%vph_ref( 97 ) = 11194.9200000000
  Mref_V%vsh_ref( 97 ) = 3620.78000000000
  Mref_V%eta_ref( 97 ) = 1.00000000000000

  Mref_V%radius_ref( 98 ) = 661925.000000000
  Mref_V%density_ref( 98 ) = 12993.0700000000
  Mref_V%vpv_ref( 98 ) = 11193.5100000000
  Mref_V%vsv_ref( 98 ) = 3619.79000000000
  Mref_V%Qkappa_ref( 98 ) = 1327.60000000000
  Mref_V%Qmu_ref( 98 ) = 104.000000000000
  Mref_V%vph_ref( 98 ) = 11193.5100000000
  Mref_V%vsh_ref( 98 ) = 3619.79000000000
  Mref_V%eta_ref( 98 ) = 1.00000000000000

  Mref_V%radius_ref( 99 ) = 668749.000000000
  Mref_V%density_ref( 99 ) = 12991.1000000000
  Mref_V%vpv_ref( 99 ) = 11192.0900000000
  Mref_V%vsv_ref( 99 ) = 3618.80000000000
  Mref_V%Qkappa_ref( 99 ) = 1327.60000000000
  Mref_V%Qmu_ref( 99 ) = 104.000000000000
  Mref_V%vph_ref( 99 ) = 11192.0900000000
  Mref_V%vsh_ref( 99 ) = 3618.80000000000
  Mref_V%eta_ref( 99 ) = 1.00000000000000

  Mref_V%radius_ref( 100 ) = 675573.000000000
  Mref_V%density_ref( 100 ) = 12989.1000000000
  Mref_V%vpv_ref( 100 ) = 11190.6500000000
  Mref_V%vsv_ref( 100 ) = 3617.79000000000
  Mref_V%Qkappa_ref( 100 ) = 1327.60000000000
  Mref_V%Qmu_ref( 100 ) = 104.000000000000
  Mref_V%vph_ref( 100 ) = 11190.6500000000
  Mref_V%vsh_ref( 100 ) = 3617.79000000000
  Mref_V%eta_ref( 100 ) = 1.00000000000000

  Mref_V%radius_ref( 101 ) = 682397.000000000
  Mref_V%density_ref( 101 ) = 12987.0800000000
  Mref_V%vpv_ref( 101 ) = 11189.1900000000
  Mref_V%vsv_ref( 101 ) = 3616.78000000000
  Mref_V%Qkappa_ref( 101 ) = 1327.60000000000
  Mref_V%Qmu_ref( 101 ) = 104.000000000000
  Mref_V%vph_ref( 101 ) = 11189.1900000000
  Mref_V%vsh_ref( 101 ) = 3616.78000000000
  Mref_V%eta_ref( 101 ) = 1.00000000000000

  Mref_V%radius_ref( 102 ) = 689221.000000000
  Mref_V%density_ref( 102 ) = 12985.0400000000
  Mref_V%vpv_ref( 102 ) = 11187.7300000000
  Mref_V%vsv_ref( 102 ) = 3615.75000000000
  Mref_V%Qkappa_ref( 102 ) = 1327.60000000000
  Mref_V%Qmu_ref( 102 ) = 104.000000000000
  Mref_V%vph_ref( 102 ) = 11187.7300000000
  Mref_V%vsh_ref( 102 ) = 3615.75000000000
  Mref_V%eta_ref( 102 ) = 1.00000000000000

  Mref_V%radius_ref( 103 ) = 696045.000000000
  Mref_V%density_ref( 103 ) = 12982.9900000000
  Mref_V%vpv_ref( 103 ) = 11186.2400000000
  Mref_V%vsv_ref( 103 ) = 3614.71000000000
  Mref_V%Qkappa_ref( 103 ) = 1327.60000000000
  Mref_V%Qmu_ref( 103 ) = 104.000000000000
  Mref_V%vph_ref( 103 ) = 11186.2400000000
  Mref_V%vsh_ref( 103 ) = 3614.71000000000
  Mref_V%eta_ref( 103 ) = 1.00000000000000

  Mref_V%radius_ref( 104 ) = 702869.000000000
  Mref_V%density_ref( 104 ) = 12980.9100000000
  Mref_V%vpv_ref( 104 ) = 11184.7500000000
  Mref_V%vsv_ref( 104 ) = 3613.67000000000
  Mref_V%Qkappa_ref( 104 ) = 1327.60000000000
  Mref_V%Qmu_ref( 104 ) = 104.000000000000
  Mref_V%vph_ref( 104 ) = 11184.7500000000
  Mref_V%vsh_ref( 104 ) = 3613.67000000000
  Mref_V%eta_ref( 104 ) = 1.00000000000000

  Mref_V%radius_ref( 105 ) = 709693.000000000
  Mref_V%density_ref( 105 ) = 12978.8100000000
  Mref_V%vpv_ref( 105 ) = 11183.2400000000
  Mref_V%vsv_ref( 105 ) = 3612.61000000000
  Mref_V%Qkappa_ref( 105 ) = 1327.60000000000
  Mref_V%Qmu_ref( 105 ) = 104.000000000000
  Mref_V%vph_ref( 105 ) = 11183.2400000000
  Mref_V%vsh_ref( 105 ) = 3612.61000000000
  Mref_V%eta_ref( 105 ) = 1.00000000000000

  Mref_V%radius_ref( 106 ) = 716517.000000000
  Mref_V%density_ref( 106 ) = 12976.6900000000
  Mref_V%vpv_ref( 106 ) = 11181.7100000000
  Mref_V%vsv_ref( 106 ) = 3611.55000000000
  Mref_V%Qkappa_ref( 106 ) = 1327.60000000000
  Mref_V%Qmu_ref( 106 ) = 104.000000000000
  Mref_V%vph_ref( 106 ) = 11181.7100000000
  Mref_V%vsh_ref( 106 ) = 3611.55000000000
  Mref_V%eta_ref( 106 ) = 1.00000000000000

  Mref_V%radius_ref( 107 ) = 723341.000000000
  Mref_V%density_ref( 107 ) = 12974.5500000000
  Mref_V%vpv_ref( 107 ) = 11180.1700000000
  Mref_V%vsv_ref( 107 ) = 3610.47000000000
  Mref_V%Qkappa_ref( 107 ) = 1327.60000000000
  Mref_V%Qmu_ref( 107 ) = 104.000000000000
  Mref_V%vph_ref( 107 ) = 11180.1700000000
  Mref_V%vsh_ref( 107 ) = 3610.47000000000
  Mref_V%eta_ref( 107 ) = 1.00000000000000

  Mref_V%radius_ref( 108 ) = 730165.000000000
  Mref_V%density_ref( 108 ) = 12972.3900000000
  Mref_V%vpv_ref( 108 ) = 11178.6100000000
  Mref_V%vsv_ref( 108 ) = 3609.38000000000
  Mref_V%Qkappa_ref( 108 ) = 1327.60000000000
  Mref_V%Qmu_ref( 108 ) = 104.000000000000
  Mref_V%vph_ref( 108 ) = 11178.6100000000
  Mref_V%vsh_ref( 108 ) = 3609.38000000000
  Mref_V%eta_ref( 108 ) = 1.00000000000000

  Mref_V%radius_ref( 109 ) = 736989.000000000
  Mref_V%density_ref( 109 ) = 12970.2100000000
  Mref_V%vpv_ref( 109 ) = 11177.0400000000
  Mref_V%vsv_ref( 109 ) = 3608.28000000000
  Mref_V%Qkappa_ref( 109 ) = 1327.60000000000
  Mref_V%Qmu_ref( 109 ) = 104.000000000000
  Mref_V%vph_ref( 109 ) = 11177.0400000000
  Mref_V%vsh_ref( 109 ) = 3608.28000000000
  Mref_V%eta_ref( 109 ) = 1.00000000000000

  Mref_V%radius_ref( 110 ) = 743813.000000000
  Mref_V%density_ref( 110 ) = 12968.0100000000
  Mref_V%vpv_ref( 110 ) = 11175.4600000000
  Mref_V%vsv_ref( 110 ) = 3607.18000000000
  Mref_V%Qkappa_ref( 110 ) = 1327.60000000000
  Mref_V%Qmu_ref( 110 ) = 104.000000000000
  Mref_V%vph_ref( 110 ) = 11175.4600000000
  Mref_V%vsh_ref( 110 ) = 3607.18000000000
  Mref_V%eta_ref( 110 ) = 1.00000000000000

  Mref_V%radius_ref( 111 ) = 750637.000000000
  Mref_V%density_ref( 111 ) = 12965.7900000000
  Mref_V%vpv_ref( 111 ) = 11173.8600000000
  Mref_V%vsv_ref( 111 ) = 3606.06000000000
  Mref_V%Qkappa_ref( 111 ) = 1327.60000000000
  Mref_V%Qmu_ref( 111 ) = 104.000000000000
  Mref_V%vph_ref( 111 ) = 11173.8600000000
  Mref_V%vsh_ref( 111 ) = 3606.06000000000
  Mref_V%eta_ref( 111 ) = 1.00000000000000

  Mref_V%radius_ref( 112 ) = 757461.000000000
  Mref_V%density_ref( 112 ) = 12963.5500000000
  Mref_V%vpv_ref( 112 ) = 11172.2500000000
  Mref_V%vsv_ref( 112 ) = 3604.93000000000
  Mref_V%Qkappa_ref( 112 ) = 1327.60000000000
  Mref_V%Qmu_ref( 112 ) = 104.000000000000
  Mref_V%vph_ref( 112 ) = 11172.2500000000
  Mref_V%vsh_ref( 112 ) = 3604.93000000000
  Mref_V%eta_ref( 112 ) = 1.00000000000000

  Mref_V%radius_ref( 113 ) = 764285.000000000
  Mref_V%density_ref( 113 ) = 12961.2900000000
  Mref_V%vpv_ref( 113 ) = 11170.6200000000
  Mref_V%vsv_ref( 113 ) = 3603.79000000000
  Mref_V%Qkappa_ref( 113 ) = 1327.60000000000
  Mref_V%Qmu_ref( 113 ) = 104.000000000000
  Mref_V%vph_ref( 113 ) = 11170.6200000000
  Mref_V%vsh_ref( 113 ) = 3603.79000000000
  Mref_V%eta_ref( 113 ) = 1.00000000000000

  Mref_V%radius_ref( 114 ) = 771109.000000000
  Mref_V%density_ref( 114 ) = 12959.0100000000
  Mref_V%vpv_ref( 114 ) = 11168.9800000000
  Mref_V%vsv_ref( 114 ) = 3602.65000000000
  Mref_V%Qkappa_ref( 114 ) = 1327.60000000000
  Mref_V%Qmu_ref( 114 ) = 104.000000000000
  Mref_V%vph_ref( 114 ) = 11168.9800000000
  Mref_V%vsh_ref( 114 ) = 3602.65000000000
  Mref_V%eta_ref( 114 ) = 1.00000000000000

  Mref_V%radius_ref( 115 ) = 777933.000000000
  Mref_V%density_ref( 115 ) = 12956.7000000000
  Mref_V%vpv_ref( 115 ) = 11167.3200000000
  Mref_V%vsv_ref( 115 ) = 3601.49000000000
  Mref_V%Qkappa_ref( 115 ) = 1327.60000000000
  Mref_V%Qmu_ref( 115 ) = 104.000000000000
  Mref_V%vph_ref( 115 ) = 11167.3200000000
  Mref_V%vsh_ref( 115 ) = 3601.49000000000
  Mref_V%eta_ref( 115 ) = 1.00000000000000

  Mref_V%radius_ref( 116 ) = 784757.000000000
  Mref_V%density_ref( 116 ) = 12954.3800000000
  Mref_V%vpv_ref( 116 ) = 11165.6500000000
  Mref_V%vsv_ref( 116 ) = 3600.32000000000
  Mref_V%Qkappa_ref( 116 ) = 1327.60000000000
  Mref_V%Qmu_ref( 116 ) = 104.000000000000
  Mref_V%vph_ref( 116 ) = 11165.6500000000
  Mref_V%vsh_ref( 116 ) = 3600.32000000000
  Mref_V%eta_ref( 116 ) = 1.00000000000000

  Mref_V%radius_ref( 117 ) = 791581.000000000
  Mref_V%density_ref( 117 ) = 12952.0400000000
  Mref_V%vpv_ref( 117 ) = 11163.9600000000
  Mref_V%vsv_ref( 117 ) = 3599.14000000000
  Mref_V%Qkappa_ref( 117 ) = 1327.60000000000
  Mref_V%Qmu_ref( 117 ) = 104.000000000000
  Mref_V%vph_ref( 117 ) = 11163.9600000000
  Mref_V%vsh_ref( 117 ) = 3599.14000000000
  Mref_V%eta_ref( 117 ) = 1.00000000000000

  Mref_V%radius_ref( 118 ) = 798405.000000000
  Mref_V%density_ref( 118 ) = 12949.6800000000
  Mref_V%vpv_ref( 118 ) = 11162.2600000000
  Mref_V%vsv_ref( 118 ) = 3597.95000000000
  Mref_V%Qkappa_ref( 118 ) = 1327.60000000000
  Mref_V%Qmu_ref( 118 ) = 104.000000000000
  Mref_V%vph_ref( 118 ) = 11162.2600000000
  Mref_V%vsh_ref( 118 ) = 3597.95000000000
  Mref_V%eta_ref( 118 ) = 1.00000000000000

  Mref_V%radius_ref( 119 ) = 805229.000000000
  Mref_V%density_ref( 119 ) = 12947.2900000000
  Mref_V%vpv_ref( 119 ) = 11160.5400000000
  Mref_V%vsv_ref( 119 ) = 3596.75000000000
  Mref_V%Qkappa_ref( 119 ) = 1327.60000000000
  Mref_V%Qmu_ref( 119 ) = 104.000000000000
  Mref_V%vph_ref( 119 ) = 11160.5400000000
  Mref_V%vsh_ref( 119 ) = 3596.75000000000
  Mref_V%eta_ref( 119 ) = 1.00000000000000

  Mref_V%radius_ref( 120 ) = 812052.000000000
  Mref_V%density_ref( 120 ) = 12944.8900000000
  Mref_V%vpv_ref( 120 ) = 11158.8100000000
  Mref_V%vsv_ref( 120 ) = 3595.54000000000
  Mref_V%Qkappa_ref( 120 ) = 1327.60000000000
  Mref_V%Qmu_ref( 120 ) = 104.000000000000
  Mref_V%vph_ref( 120 ) = 11158.8100000000
  Mref_V%vsh_ref( 120 ) = 3595.54000000000
  Mref_V%eta_ref( 120 ) = 1.00000000000000

  Mref_V%radius_ref( 121 ) = 818876.000000000
  Mref_V%density_ref( 121 ) = 12942.4700000000
  Mref_V%vpv_ref( 121 ) = 11157.0700000000
  Mref_V%vsv_ref( 121 ) = 3594.32000000000
  Mref_V%Qkappa_ref( 121 ) = 1327.60000000000
  Mref_V%Qmu_ref( 121 ) = 104.000000000000
  Mref_V%vph_ref( 121 ) = 11157.0700000000
  Mref_V%vsh_ref( 121 ) = 3594.32000000000
  Mref_V%eta_ref( 121 ) = 1.00000000000000

  Mref_V%radius_ref( 122 ) = 825700.000000000
  Mref_V%density_ref( 122 ) = 12940.0200000000
  Mref_V%vpv_ref( 122 ) = 11155.3100000000
  Mref_V%vsv_ref( 122 ) = 3593.10000000000
  Mref_V%Qkappa_ref( 122 ) = 1327.60000000000
  Mref_V%Qmu_ref( 122 ) = 104.000000000000
  Mref_V%vph_ref( 122 ) = 11155.3100000000
  Mref_V%vsh_ref( 122 ) = 3593.10000000000
  Mref_V%eta_ref( 122 ) = 1.00000000000000

  Mref_V%radius_ref( 123 ) = 832524.000000000
  Mref_V%density_ref( 123 ) = 12937.5600000000
  Mref_V%vpv_ref( 123 ) = 11153.5400000000
  Mref_V%vsv_ref( 123 ) = 3591.86000000000
  Mref_V%Qkappa_ref( 123 ) = 1327.60000000000
  Mref_V%Qmu_ref( 123 ) = 104.000000000000
  Mref_V%vph_ref( 123 ) = 11153.5400000000
  Mref_V%vsh_ref( 123 ) = 3591.86000000000
  Mref_V%eta_ref( 123 ) = 1.00000000000000

  Mref_V%radius_ref( 124 ) = 839348.000000000
  Mref_V%density_ref( 124 ) = 12935.0800000000
  Mref_V%vpv_ref( 124 ) = 11151.7500000000
  Mref_V%vsv_ref( 124 ) = 3590.61000000000
  Mref_V%Qkappa_ref( 124 ) = 1327.60000000000
  Mref_V%Qmu_ref( 124 ) = 104.000000000000
  Mref_V%vph_ref( 124 ) = 11151.7500000000
  Mref_V%vsh_ref( 124 ) = 3590.61000000000
  Mref_V%eta_ref( 124 ) = 1.00000000000000

  Mref_V%radius_ref( 125 ) = 846172.000000000
  Mref_V%density_ref( 125 ) = 12932.5700000000
  Mref_V%vpv_ref( 125 ) = 11149.9400000000
  Mref_V%vsv_ref( 125 ) = 3589.34000000000
  Mref_V%Qkappa_ref( 125 ) = 1327.60000000000
  Mref_V%Qmu_ref( 125 ) = 104.000000000000
  Mref_V%vph_ref( 125 ) = 11149.9400000000
  Mref_V%vsh_ref( 125 ) = 3589.34000000000
  Mref_V%eta_ref( 125 ) = 1.00000000000000

  Mref_V%radius_ref( 126 ) = 852996.000000000
  Mref_V%density_ref( 126 ) = 12930.0500000000
  Mref_V%vpv_ref( 126 ) = 11148.1300000000
  Mref_V%vsv_ref( 126 ) = 3588.07000000000
  Mref_V%Qkappa_ref( 126 ) = 1327.60000000000
  Mref_V%Qmu_ref( 126 ) = 104.000000000000
  Mref_V%vph_ref( 126 ) = 11148.1300000000
  Mref_V%vsh_ref( 126 ) = 3588.07000000000
  Mref_V%eta_ref( 126 ) = 1.00000000000000

  Mref_V%radius_ref( 127 ) = 859820.000000000
  Mref_V%density_ref( 127 ) = 12927.5000000000
  Mref_V%vpv_ref( 127 ) = 11146.2900000000
  Mref_V%vsv_ref( 127 ) = 3586.79000000000
  Mref_V%Qkappa_ref( 127 ) = 1327.60000000000
  Mref_V%Qmu_ref( 127 ) = 104.000000000000
  Mref_V%vph_ref( 127 ) = 11146.2900000000
  Mref_V%vsh_ref( 127 ) = 3586.79000000000
  Mref_V%eta_ref( 127 ) = 1.00000000000000

  Mref_V%radius_ref( 128 ) = 866644.000000000
  Mref_V%density_ref( 128 ) = 12924.9400000000
  Mref_V%vpv_ref( 128 ) = 11144.4500000000
  Mref_V%vsv_ref( 128 ) = 3585.50000000000
  Mref_V%Qkappa_ref( 128 ) = 1327.60000000000
  Mref_V%Qmu_ref( 128 ) = 104.000000000000
  Mref_V%vph_ref( 128 ) = 11144.4500000000
  Mref_V%vsh_ref( 128 ) = 3585.50000000000
  Mref_V%eta_ref( 128 ) = 1.00000000000000

  Mref_V%radius_ref( 129 ) = 873468.000000000
  Mref_V%density_ref( 129 ) = 12922.3500000000
  Mref_V%vpv_ref( 129 ) = 11142.5800000000
  Mref_V%vsv_ref( 129 ) = 3584.20000000000
  Mref_V%Qkappa_ref( 129 ) = 1327.60000000000
  Mref_V%Qmu_ref( 129 ) = 104.000000000000
  Mref_V%vph_ref( 129 ) = 11142.5800000000
  Mref_V%vsh_ref( 129 ) = 3584.20000000000
  Mref_V%eta_ref( 129 ) = 1.00000000000000

  Mref_V%radius_ref( 130 ) = 880292.000000000
  Mref_V%density_ref( 130 ) = 12919.7500000000
  Mref_V%vpv_ref( 130 ) = 11140.7100000000
  Mref_V%vsv_ref( 130 ) = 3582.89000000000
  Mref_V%Qkappa_ref( 130 ) = 1327.60000000000
  Mref_V%Qmu_ref( 130 ) = 104.000000000000
  Mref_V%vph_ref( 130 ) = 11140.7100000000
  Mref_V%vsh_ref( 130 ) = 3582.89000000000
  Mref_V%eta_ref( 130 ) = 1.00000000000000

  Mref_V%radius_ref( 131 ) = 887116.000000000
  Mref_V%density_ref( 131 ) = 12917.1200000000
  Mref_V%vpv_ref( 131 ) = 11138.8200000000
  Mref_V%vsv_ref( 131 ) = 3581.57000000000
  Mref_V%Qkappa_ref( 131 ) = 1327.60000000000
  Mref_V%Qmu_ref( 131 ) = 104.000000000000
  Mref_V%vph_ref( 131 ) = 11138.8200000000
  Mref_V%vsh_ref( 131 ) = 3581.57000000000
  Mref_V%eta_ref( 131 ) = 1.00000000000000

  Mref_V%radius_ref( 132 ) = 893940.000000000
  Mref_V%density_ref( 132 ) = 12914.4700000000
  Mref_V%vpv_ref( 132 ) = 11136.9100000000
  Mref_V%vsv_ref( 132 ) = 3580.24000000000
  Mref_V%Qkappa_ref( 132 ) = 1327.60000000000
  Mref_V%Qmu_ref( 132 ) = 104.000000000000
  Mref_V%vph_ref( 132 ) = 11136.9100000000
  Mref_V%vsh_ref( 132 ) = 3580.24000000000
  Mref_V%eta_ref( 132 ) = 1.00000000000000

  Mref_V%radius_ref( 133 ) = 900764.000000000
  Mref_V%density_ref( 133 ) = 12911.8100000000
  Mref_V%vpv_ref( 133 ) = 11134.9900000000
  Mref_V%vsv_ref( 133 ) = 3578.90000000000
  Mref_V%Qkappa_ref( 133 ) = 1327.60000000000
  Mref_V%Qmu_ref( 133 ) = 104.000000000000
  Mref_V%vph_ref( 133 ) = 11134.9900000000
  Mref_V%vsh_ref( 133 ) = 3578.90000000000
  Mref_V%eta_ref( 133 ) = 1.00000000000000

  Mref_V%radius_ref( 134 ) = 907588.000000000
  Mref_V%density_ref( 134 ) = 12909.1200000000
  Mref_V%vpv_ref( 134 ) = 11133.0600000000
  Mref_V%vsv_ref( 134 ) = 3577.54000000000
  Mref_V%Qkappa_ref( 134 ) = 1327.60000000000
  Mref_V%Qmu_ref( 134 ) = 104.000000000000
  Mref_V%vph_ref( 134 ) = 11133.0600000000
  Mref_V%vsh_ref( 134 ) = 3577.54000000000
  Mref_V%eta_ref( 134 ) = 1.00000000000000

  Mref_V%radius_ref( 135 ) = 914412.000000000
  Mref_V%density_ref( 135 ) = 12906.4100000000
  Mref_V%vpv_ref( 135 ) = 11131.1100000000
  Mref_V%vsv_ref( 135 ) = 3576.18000000000
  Mref_V%Qkappa_ref( 135 ) = 1327.60000000000
  Mref_V%Qmu_ref( 135 ) = 104.000000000000
  Mref_V%vph_ref( 135 ) = 11131.1100000000
  Mref_V%vsh_ref( 135 ) = 3576.18000000000
  Mref_V%eta_ref( 135 ) = 1.00000000000000

  Mref_V%radius_ref( 136 ) = 921236.000000000
  Mref_V%density_ref( 136 ) = 12903.6800000000
  Mref_V%vpv_ref( 136 ) = 11129.1400000000
  Mref_V%vsv_ref( 136 ) = 3574.81000000000
  Mref_V%Qkappa_ref( 136 ) = 1327.60000000000
  Mref_V%Qmu_ref( 136 ) = 104.000000000000
  Mref_V%vph_ref( 136 ) = 11129.1400000000
  Mref_V%vsh_ref( 136 ) = 3574.81000000000
  Mref_V%eta_ref( 136 ) = 1.00000000000000

  Mref_V%radius_ref( 137 ) = 928060.000000000
  Mref_V%density_ref( 137 ) = 12900.9400000000
  Mref_V%vpv_ref( 137 ) = 11127.1600000000
  Mref_V%vsv_ref( 137 ) = 3573.43000000000
  Mref_V%Qkappa_ref( 137 ) = 1327.60000000000
  Mref_V%Qmu_ref( 137 ) = 104.000000000000
  Mref_V%vph_ref( 137 ) = 11127.1600000000
  Mref_V%vsh_ref( 137 ) = 3573.43000000000
  Mref_V%eta_ref( 137 ) = 1.00000000000000

  Mref_V%radius_ref( 138 ) = 934884.000000000
  Mref_V%density_ref( 138 ) = 12898.1700000000
  Mref_V%vpv_ref( 138 ) = 11125.1700000000
  Mref_V%vsv_ref( 138 ) = 3572.03000000000
  Mref_V%Qkappa_ref( 138 ) = 1327.60000000000
  Mref_V%Qmu_ref( 138 ) = 104.000000000000
  Mref_V%vph_ref( 138 ) = 11125.1700000000
  Mref_V%vsh_ref( 138 ) = 3572.03000000000
  Mref_V%eta_ref( 138 ) = 1.00000000000000

  Mref_V%radius_ref( 139 ) = 941708.000000000
  Mref_V%density_ref( 139 ) = 12895.3800000000
  Mref_V%vpv_ref( 139 ) = 11123.1600000000
  Mref_V%vsv_ref( 139 ) = 3570.63000000000
  Mref_V%Qkappa_ref( 139 ) = 1327.60000000000
  Mref_V%Qmu_ref( 139 ) = 104.000000000000
  Mref_V%vph_ref( 139 ) = 11123.1600000000
  Mref_V%vsh_ref( 139 ) = 3570.63000000000
  Mref_V%eta_ref( 139 ) = 1.00000000000000

  Mref_V%radius_ref( 140 ) = 948532.000000000
  Mref_V%density_ref( 140 ) = 12892.5700000000
  Mref_V%vpv_ref( 140 ) = 11121.1400000000
  Mref_V%vsv_ref( 140 ) = 3569.22000000000
  Mref_V%Qkappa_ref( 140 ) = 1327.60000000000
  Mref_V%Qmu_ref( 140 ) = 104.000000000000
  Mref_V%vph_ref( 140 ) = 11121.1400000000
  Mref_V%vsh_ref( 140 ) = 3569.22000000000
  Mref_V%eta_ref( 140 ) = 1.00000000000000

  Mref_V%radius_ref( 141 ) = 955356.000000000
  Mref_V%density_ref( 141 ) = 12889.7400000000
  Mref_V%vpv_ref( 141 ) = 11119.1000000000
  Mref_V%vsv_ref( 141 ) = 3567.79000000000
  Mref_V%Qkappa_ref( 141 ) = 1327.60000000000
  Mref_V%Qmu_ref( 141 ) = 104.000000000000
  Mref_V%vph_ref( 141 ) = 11119.1000000000
  Mref_V%vsh_ref( 141 ) = 3567.79000000000
  Mref_V%eta_ref( 141 ) = 1.00000000000000

  Mref_V%radius_ref( 142 ) = 962180.000000000
  Mref_V%density_ref( 142 ) = 12886.8900000000
  Mref_V%vpv_ref( 142 ) = 11117.0500000000
  Mref_V%vsv_ref( 142 ) = 3566.36000000000
  Mref_V%Qkappa_ref( 142 ) = 1327.60000000000
  Mref_V%Qmu_ref( 142 ) = 104.000000000000
  Mref_V%vph_ref( 142 ) = 11117.0500000000
  Mref_V%vsh_ref( 142 ) = 3566.36000000000
  Mref_V%eta_ref( 142 ) = 1.00000000000000

  Mref_V%radius_ref( 143 ) = 969004.000000000
  Mref_V%density_ref( 143 ) = 12884.0200000000
  Mref_V%vpv_ref( 143 ) = 11114.9900000000
  Mref_V%vsv_ref( 143 ) = 3564.91000000000
  Mref_V%Qkappa_ref( 143 ) = 1327.60000000000
  Mref_V%Qmu_ref( 143 ) = 104.000000000000
  Mref_V%vph_ref( 143 ) = 11114.9900000000
  Mref_V%vsh_ref( 143 ) = 3564.91000000000
  Mref_V%eta_ref( 143 ) = 1.00000000000000

  Mref_V%radius_ref( 144 ) = 975828.000000000
  Mref_V%density_ref( 144 ) = 12881.1300000000
  Mref_V%vpv_ref( 144 ) = 11112.9000000000
  Mref_V%vsv_ref( 144 ) = 3563.46000000000
  Mref_V%Qkappa_ref( 144 ) = 1327.60000000000
  Mref_V%Qmu_ref( 144 ) = 104.000000000000
  Mref_V%vph_ref( 144 ) = 11112.9000000000
  Mref_V%vsh_ref( 144 ) = 3563.46000000000
  Mref_V%eta_ref( 144 ) = 1.00000000000000

  Mref_V%radius_ref( 145 ) = 982652.000000000
  Mref_V%density_ref( 145 ) = 12878.2200000000
  Mref_V%vpv_ref( 145 ) = 11110.8100000000
  Mref_V%vsv_ref( 145 ) = 3562.00000000000
  Mref_V%Qkappa_ref( 145 ) = 1327.60000000000
  Mref_V%Qmu_ref( 145 ) = 104.000000000000
  Mref_V%vph_ref( 145 ) = 11110.8100000000
  Mref_V%vsh_ref( 145 ) = 3562.00000000000
  Mref_V%eta_ref( 145 ) = 1.00000000000000

  Mref_V%radius_ref( 146 ) = 989476.000000000
  Mref_V%density_ref( 146 ) = 12875.2900000000
  Mref_V%vpv_ref( 146 ) = 11108.7000000000
  Mref_V%vsv_ref( 146 ) = 3560.52000000000
  Mref_V%Qkappa_ref( 146 ) = 1327.60000000000
  Mref_V%Qmu_ref( 146 ) = 104.000000000000
  Mref_V%vph_ref( 146 ) = 11108.7000000000
  Mref_V%vsh_ref( 146 ) = 3560.52000000000
  Mref_V%eta_ref( 146 ) = 1.00000000000000

  Mref_V%radius_ref( 147 ) = 996300.000000000
  Mref_V%density_ref( 147 ) = 12872.3400000000
  Mref_V%vpv_ref( 147 ) = 11106.5700000000
  Mref_V%vsv_ref( 147 ) = 3559.04000000000
  Mref_V%Qkappa_ref( 147 ) = 1327.60000000000
  Mref_V%Qmu_ref( 147 ) = 104.000000000000
  Mref_V%vph_ref( 147 ) = 11106.5700000000
  Mref_V%vsh_ref( 147 ) = 3559.04000000000
  Mref_V%eta_ref( 147 ) = 1.00000000000000

  Mref_V%radius_ref( 148 ) = 1003124.00000000
  Mref_V%density_ref( 148 ) = 12869.3700000000
  Mref_V%vpv_ref( 148 ) = 11104.4400000000
  Mref_V%vsv_ref( 148 ) = 3557.54000000000
  Mref_V%Qkappa_ref( 148 ) = 1327.60000000000
  Mref_V%Qmu_ref( 148 ) = 104.000000000000
  Mref_V%vph_ref( 148 ) = 11104.4400000000
  Mref_V%vsh_ref( 148 ) = 3557.54000000000
  Mref_V%eta_ref( 148 ) = 1.00000000000000

  Mref_V%radius_ref( 149 ) = 1009948.00000000
  Mref_V%density_ref( 149 ) = 12866.3800000000
  Mref_V%vpv_ref( 149 ) = 11102.2800000000
  Mref_V%vsv_ref( 149 ) = 3556.04000000000
  Mref_V%Qkappa_ref( 149 ) = 1327.60000000000
  Mref_V%Qmu_ref( 149 ) = 104.000000000000
  Mref_V%vph_ref( 149 ) = 11102.2800000000
  Mref_V%vsh_ref( 149 ) = 3556.04000000000
  Mref_V%eta_ref( 149 ) = 1.00000000000000

  Mref_V%radius_ref( 150 ) = 1016772.00000000
  Mref_V%density_ref( 150 ) = 12863.3700000000
  Mref_V%vpv_ref( 150 ) = 11100.1100000000
  Mref_V%vsv_ref( 150 ) = 3554.52000000000
  Mref_V%Qkappa_ref( 150 ) = 1327.60000000000
  Mref_V%Qmu_ref( 150 ) = 104.000000000000
  Mref_V%vph_ref( 150 ) = 11100.1100000000
  Mref_V%vsh_ref( 150 ) = 3554.52000000000
  Mref_V%eta_ref( 150 ) = 1.00000000000000

  Mref_V%radius_ref( 151 ) = 1023596.00000000
  Mref_V%density_ref( 151 ) = 12860.3400000000
  Mref_V%vpv_ref( 151 ) = 11097.9300000000
  Mref_V%vsv_ref( 151 ) = 3553.00000000000
  Mref_V%Qkappa_ref( 151 ) = 1327.60000000000
  Mref_V%Qmu_ref( 151 ) = 104.000000000000
  Mref_V%vph_ref( 151 ) = 11097.9300000000
  Mref_V%vsh_ref( 151 ) = 3553.00000000000
  Mref_V%eta_ref( 151 ) = 1.00000000000000

  Mref_V%radius_ref( 152 ) = 1030420.00000000
  Mref_V%density_ref( 152 ) = 12857.2900000000
  Mref_V%vpv_ref( 152 ) = 11095.7300000000
  Mref_V%vsv_ref( 152 ) = 3551.46000000000
  Mref_V%Qkappa_ref( 152 ) = 1327.60000000000
  Mref_V%Qmu_ref( 152 ) = 104.000000000000
  Mref_V%vph_ref( 152 ) = 11095.7300000000
  Mref_V%vsh_ref( 152 ) = 3551.46000000000
  Mref_V%eta_ref( 152 ) = 1.00000000000000

  Mref_V%radius_ref( 153 ) = 1037244.00000000
  Mref_V%density_ref( 153 ) = 12854.2100000000
  Mref_V%vpv_ref( 153 ) = 11093.5200000000
  Mref_V%vsv_ref( 153 ) = 3549.91000000000
  Mref_V%Qkappa_ref( 153 ) = 1327.60000000000
  Mref_V%Qmu_ref( 153 ) = 104.000000000000
  Mref_V%vph_ref( 153 ) = 11093.5200000000
  Mref_V%vsh_ref( 153 ) = 3549.91000000000
  Mref_V%eta_ref( 153 ) = 1.00000000000000

  Mref_V%radius_ref( 154 ) = 1044067.00000000
  Mref_V%density_ref( 154 ) = 12851.1200000000
  Mref_V%vpv_ref( 154 ) = 11091.2900000000
  Mref_V%vsv_ref( 154 ) = 3548.36000000000
  Mref_V%Qkappa_ref( 154 ) = 1327.60000000000
  Mref_V%Qmu_ref( 154 ) = 104.000000000000
  Mref_V%vph_ref( 154 ) = 11091.2900000000
  Mref_V%vsh_ref( 154 ) = 3548.36000000000
  Mref_V%eta_ref( 154 ) = 1.00000000000000

  Mref_V%radius_ref( 155 ) = 1050891.00000000
  Mref_V%density_ref( 155 ) = 12848.0100000000
  Mref_V%vpv_ref( 155 ) = 11089.0500000000
  Mref_V%vsv_ref( 155 ) = 3546.79000000000
  Mref_V%Qkappa_ref( 155 ) = 1327.60000000000
  Mref_V%Qmu_ref( 155 ) = 104.000000000000
  Mref_V%vph_ref( 155 ) = 11089.0500000000
  Mref_V%vsh_ref( 155 ) = 3546.79000000000
  Mref_V%eta_ref( 155 ) = 1.00000000000000

  Mref_V%radius_ref( 156 ) = 1057715.00000000
  Mref_V%density_ref( 156 ) = 12844.8800000000
  Mref_V%vpv_ref( 156 ) = 11086.8000000000
  Mref_V%vsv_ref( 156 ) = 3545.21000000000
  Mref_V%Qkappa_ref( 156 ) = 1327.60000000000
  Mref_V%Qmu_ref( 156 ) = 104.000000000000
  Mref_V%vph_ref( 156 ) = 11086.8000000000
  Mref_V%vsh_ref( 156 ) = 3545.21000000000
  Mref_V%eta_ref( 156 ) = 1.00000000000000

  Mref_V%radius_ref( 157 ) = 1064539.00000000
  Mref_V%density_ref( 157 ) = 12841.7200000000
  Mref_V%vpv_ref( 157 ) = 11084.5300000000
  Mref_V%vsv_ref( 157 ) = 3543.63000000000
  Mref_V%Qkappa_ref( 157 ) = 1327.60000000000
  Mref_V%Qmu_ref( 157 ) = 104.000000000000
  Mref_V%vph_ref( 157 ) = 11084.5300000000
  Mref_V%vsh_ref( 157 ) = 3543.63000000000
  Mref_V%eta_ref( 157 ) = 1.00000000000000

  Mref_V%radius_ref( 158 ) = 1071363.00000000
  Mref_V%density_ref( 158 ) = 12838.5500000000
  Mref_V%vpv_ref( 158 ) = 11082.2400000000
  Mref_V%vsv_ref( 158 ) = 3542.03000000000
  Mref_V%Qkappa_ref( 158 ) = 1327.60000000000
  Mref_V%Qmu_ref( 158 ) = 104.000000000000
  Mref_V%vph_ref( 158 ) = 11082.2400000000
  Mref_V%vsh_ref( 158 ) = 3542.03000000000
  Mref_V%eta_ref( 158 ) = 1.00000000000000

  Mref_V%radius_ref( 159 ) = 1078187.00000000
  Mref_V%density_ref( 159 ) = 12835.3500000000
  Mref_V%vpv_ref( 159 ) = 11079.9400000000
  Mref_V%vsv_ref( 159 ) = 3540.42000000000
  Mref_V%Qkappa_ref( 159 ) = 1327.60000000000
  Mref_V%Qmu_ref( 159 ) = 104.000000000000
  Mref_V%vph_ref( 159 ) = 11079.9400000000
  Mref_V%vsh_ref( 159 ) = 3540.42000000000
  Mref_V%eta_ref( 159 ) = 1.00000000000000

  Mref_V%radius_ref( 160 ) = 1085011.00000000
  Mref_V%density_ref( 160 ) = 12832.1400000000
  Mref_V%vpv_ref( 160 ) = 11077.6300000000
  Mref_V%vsv_ref( 160 ) = 3538.81000000000
  Mref_V%Qkappa_ref( 160 ) = 1327.60000000000
  Mref_V%Qmu_ref( 160 ) = 104.000000000000
  Mref_V%vph_ref( 160 ) = 11077.6300000000
  Mref_V%vsh_ref( 160 ) = 3538.81000000000
  Mref_V%eta_ref( 160 ) = 1.00000000000000

  Mref_V%radius_ref( 161 ) = 1091835.00000000
  Mref_V%density_ref( 161 ) = 12828.9100000000
  Mref_V%vpv_ref( 161 ) = 11075.3000000000
  Mref_V%vsv_ref( 161 ) = 3537.18000000000
  Mref_V%Qkappa_ref( 161 ) = 1327.60000000000
  Mref_V%Qmu_ref( 161 ) = 104.000000000000
  Mref_V%vph_ref( 161 ) = 11075.3000000000
  Mref_V%vsh_ref( 161 ) = 3537.18000000000
  Mref_V%eta_ref( 161 ) = 1.00000000000000

  Mref_V%radius_ref( 162 ) = 1098659.00000000
  Mref_V%density_ref( 162 ) = 12825.6500000000
  Mref_V%vpv_ref( 162 ) = 11072.9500000000
  Mref_V%vsv_ref( 162 ) = 3535.54000000000
  Mref_V%Qkappa_ref( 162 ) = 1327.60000000000
  Mref_V%Qmu_ref( 162 ) = 104.000000000000
  Mref_V%vph_ref( 162 ) = 11072.9500000000
  Mref_V%vsh_ref( 162 ) = 3535.54000000000
  Mref_V%eta_ref( 162 ) = 1.00000000000000

  Mref_V%radius_ref( 163 ) = 1105483.00000000
  Mref_V%density_ref( 163 ) = 12822.3800000000
  Mref_V%vpv_ref( 163 ) = 11070.5900000000
  Mref_V%vsv_ref( 163 ) = 3533.89000000000
  Mref_V%Qkappa_ref( 163 ) = 1327.60000000000
  Mref_V%Qmu_ref( 163 ) = 104.000000000000
  Mref_V%vph_ref( 163 ) = 11070.5900000000
  Mref_V%vsh_ref( 163 ) = 3533.89000000000
  Mref_V%eta_ref( 163 ) = 1.00000000000000

  Mref_V%radius_ref( 164 ) = 1112307.00000000
  Mref_V%density_ref( 164 ) = 12819.0800000000
  Mref_V%vpv_ref( 164 ) = 11068.2200000000
  Mref_V%vsv_ref( 164 ) = 3532.23000000000
  Mref_V%Qkappa_ref( 164 ) = 1327.60000000000
  Mref_V%Qmu_ref( 164 ) = 104.000000000000
  Mref_V%vph_ref( 164 ) = 11068.2200000000
  Mref_V%vsh_ref( 164 ) = 3532.23000000000
  Mref_V%eta_ref( 164 ) = 1.00000000000000

  Mref_V%radius_ref( 165 ) = 1119131.00000000
  Mref_V%density_ref( 165 ) = 12815.7600000000
  Mref_V%vpv_ref( 165 ) = 11065.8300000000
  Mref_V%vsv_ref( 165 ) = 3530.57000000000
  Mref_V%Qkappa_ref( 165 ) = 1327.60000000000
  Mref_V%Qmu_ref( 165 ) = 104.000000000000
  Mref_V%vph_ref( 165 ) = 11065.8300000000
  Mref_V%vsh_ref( 165 ) = 3530.57000000000
  Mref_V%eta_ref( 165 ) = 1.00000000000000

  Mref_V%radius_ref( 166 ) = 1125955.00000000
  Mref_V%density_ref( 166 ) = 12812.4300000000
  Mref_V%vpv_ref( 166 ) = 11063.4300000000
  Mref_V%vsv_ref( 166 ) = 3528.89000000000
  Mref_V%Qkappa_ref( 166 ) = 1327.60000000000
  Mref_V%Qmu_ref( 166 ) = 104.000000000000
  Mref_V%vph_ref( 166 ) = 11063.4300000000
  Mref_V%vsh_ref( 166 ) = 3528.89000000000
  Mref_V%eta_ref( 166 ) = 1.00000000000000

  Mref_V%radius_ref( 167 ) = 1132779.00000000
  Mref_V%density_ref( 167 ) = 12809.0700000000
  Mref_V%vpv_ref( 167 ) = 11061.0200000000
  Mref_V%vsv_ref( 167 ) = 3527.20000000000
  Mref_V%Qkappa_ref( 167 ) = 1327.60000000000
  Mref_V%Qmu_ref( 167 ) = 104.000000000000
  Mref_V%vph_ref( 167 ) = 11061.0200000000
  Mref_V%vsh_ref( 167 ) = 3527.20000000000
  Mref_V%eta_ref( 167 ) = 1.00000000000000

  Mref_V%radius_ref( 168 ) = 1139603.00000000
  Mref_V%density_ref( 168 ) = 12805.7000000000
  Mref_V%vpv_ref( 168 ) = 11058.5800000000
  Mref_V%vsv_ref( 168 ) = 3525.50000000000
  Mref_V%Qkappa_ref( 168 ) = 1327.60000000000
  Mref_V%Qmu_ref( 168 ) = 104.000000000000
  Mref_V%vph_ref( 168 ) = 11058.5800000000
  Mref_V%vsh_ref( 168 ) = 3525.50000000000
  Mref_V%eta_ref( 168 ) = 1.00000000000000

  Mref_V%radius_ref( 169 ) = 1146427.00000000
  Mref_V%density_ref( 169 ) = 12802.3000000000
  Mref_V%vpv_ref( 169 ) = 11056.1400000000
  Mref_V%vsv_ref( 169 ) = 3523.79000000000
  Mref_V%Qkappa_ref( 169 ) = 1327.60000000000
  Mref_V%Qmu_ref( 169 ) = 104.000000000000
  Mref_V%vph_ref( 169 ) = 11056.1400000000
  Mref_V%vsh_ref( 169 ) = 3523.79000000000
  Mref_V%eta_ref( 169 ) = 1.00000000000000

  Mref_V%radius_ref( 170 ) = 1153251.00000000
  Mref_V%density_ref( 170 ) = 12798.8800000000
  Mref_V%vpv_ref( 170 ) = 11053.6800000000
  Mref_V%vsv_ref( 170 ) = 3522.07000000000
  Mref_V%Qkappa_ref( 170 ) = 1327.60000000000
  Mref_V%Qmu_ref( 170 ) = 104.000000000000
  Mref_V%vph_ref( 170 ) = 11053.6800000000
  Mref_V%vsh_ref( 170 ) = 3522.07000000000
  Mref_V%eta_ref( 170 ) = 1.00000000000000

  Mref_V%radius_ref( 171 ) = 1160075.00000000
  Mref_V%density_ref( 171 ) = 12795.4400000000
  Mref_V%vpv_ref( 171 ) = 11051.2000000000
  Mref_V%vsv_ref( 171 ) = 3520.34000000000
  Mref_V%Qkappa_ref( 171 ) = 1327.60000000000
  Mref_V%Qmu_ref( 171 ) = 104.000000000000
  Mref_V%vph_ref( 171 ) = 11051.2000000000
  Mref_V%vsh_ref( 171 ) = 3520.34000000000
  Mref_V%eta_ref( 171 ) = 1.00000000000000

  Mref_V%radius_ref( 172 ) = 1166899.00000000
  Mref_V%density_ref( 172 ) = 12791.9900000000
  Mref_V%vpv_ref( 172 ) = 11048.7100000000
  Mref_V%vsv_ref( 172 ) = 3518.60000000000
  Mref_V%Qkappa_ref( 172 ) = 1327.60000000000
  Mref_V%Qmu_ref( 172 ) = 104.000000000000
  Mref_V%vph_ref( 172 ) = 11048.7100000000
  Mref_V%vsh_ref( 172 ) = 3518.60000000000
  Mref_V%eta_ref( 172 ) = 1.00000000000000

  Mref_V%radius_ref( 173 ) = 1173723.00000000
  Mref_V%density_ref( 173 ) = 12788.5100000000
  Mref_V%vpv_ref( 173 ) = 11046.2100000000
  Mref_V%vsv_ref( 173 ) = 3516.85000000000
  Mref_V%Qkappa_ref( 173 ) = 1327.60000000000
  Mref_V%Qmu_ref( 173 ) = 104.000000000000
  Mref_V%vph_ref( 173 ) = 11046.2100000000
  Mref_V%vsh_ref( 173 ) = 3516.85000000000
  Mref_V%eta_ref( 173 ) = 1.00000000000000

  Mref_V%radius_ref( 174 ) = 1180547.00000000
  Mref_V%density_ref( 174 ) = 12785.0100000000
  Mref_V%vpv_ref( 174 ) = 11043.6900000000
  Mref_V%vsv_ref( 174 ) = 3515.09000000000
  Mref_V%Qkappa_ref( 174 ) = 1327.60000000000
  Mref_V%Qmu_ref( 174 ) = 104.000000000000
  Mref_V%vph_ref( 174 ) = 11043.6900000000
  Mref_V%vsh_ref( 174 ) = 3515.09000000000
  Mref_V%eta_ref( 174 ) = 1.00000000000000

  Mref_V%radius_ref( 175 ) = 1187371.00000000
  Mref_V%density_ref( 175 ) = 12781.4900000000
  Mref_V%vpv_ref( 175 ) = 11041.1600000000
  Mref_V%vsv_ref( 175 ) = 3513.32000000000
  Mref_V%Qkappa_ref( 175 ) = 1327.60000000000
  Mref_V%Qmu_ref( 175 ) = 104.000000000000
  Mref_V%vph_ref( 175 ) = 11041.1600000000
  Mref_V%vsh_ref( 175 ) = 3513.32000000000
  Mref_V%eta_ref( 175 ) = 1.00000000000000

  Mref_V%radius_ref( 176 ) = 1194195.00000000
  Mref_V%density_ref( 176 ) = 12777.9500000000
  Mref_V%vpv_ref( 176 ) = 11038.6100000000
  Mref_V%vsv_ref( 176 ) = 3511.54000000000
  Mref_V%Qkappa_ref( 176 ) = 1327.60000000000
  Mref_V%Qmu_ref( 176 ) = 104.000000000000
  Mref_V%vph_ref( 176 ) = 11038.6100000000
  Mref_V%vsh_ref( 176 ) = 3511.54000000000
  Mref_V%eta_ref( 176 ) = 1.00000000000000

  Mref_V%radius_ref( 177 ) = 1201019.00000000
  Mref_V%density_ref( 177 ) = 12774.4000000000
  Mref_V%vpv_ref( 177 ) = 11036.0500000000
  Mref_V%vsv_ref( 177 ) = 3509.75000000000
  Mref_V%Qkappa_ref( 177 ) = 1327.60000000000
  Mref_V%Qmu_ref( 177 ) = 104.000000000000
  Mref_V%vph_ref( 177 ) = 11036.0500000000
  Mref_V%vsh_ref( 177 ) = 3509.75000000000
  Mref_V%eta_ref( 177 ) = 1.00000000000000

  Mref_V%radius_ref( 178 ) = 1207843.00000000
  Mref_V%density_ref( 178 ) = 12770.8200000000
  Mref_V%vpv_ref( 178 ) = 11033.4700000000
  Mref_V%vsv_ref( 178 ) = 3507.95000000000
  Mref_V%Qkappa_ref( 178 ) = 1327.60000000000
  Mref_V%Qmu_ref( 178 ) = 104.000000000000
  Mref_V%vph_ref( 178 ) = 11033.4700000000
  Mref_V%vsh_ref( 178 ) = 3507.95000000000
  Mref_V%eta_ref( 178 ) = 1.00000000000000

  Mref_V%radius_ref( 179 ) = 1214667.00000000
  Mref_V%density_ref( 179 ) = 12767.2200000000
  Mref_V%vpv_ref( 179 ) = 11030.8800000000
  Mref_V%vsv_ref( 179 ) = 3506.13000000000
  Mref_V%Qkappa_ref( 179 ) = 1327.60000000000
  Mref_V%Qmu_ref( 179 ) = 104.000000000000
  Mref_V%vph_ref( 179 ) = 11030.8800000000
  Mref_V%vsh_ref( 179 ) = 3506.13000000000
  Mref_V%eta_ref( 179 ) = 1.00000000000000

  Mref_V%radius_ref( 180 ) = 1221491.00000000
  Mref_V%density_ref( 180 ) = 12763.6000000000
  Mref_V%vpv_ref( 180 ) = 11028.2700000000
  Mref_V%vsv_ref( 180 ) = 3504.31000000000
  Mref_V%Qkappa_ref( 180 ) = 1327.60000000000
  Mref_V%Qmu_ref( 180 ) = 104.000000000000
  Mref_V%vph_ref( 180 ) = 11028.2700000000
  Mref_V%vsh_ref( 180 ) = 3504.31000000000
  Mref_V%eta_ref( 180 ) = 1.00000000000000

  Mref_V%radius_ref( 181 ) = 1221491.00000000
  Mref_V%density_ref( 181 ) = 12166.3500000000
  Mref_V%vpv_ref( 181 ) = 10355.6900000000
  Mref_V%vsv_ref( 181 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 181 ) = 57822.5000000000
  Mref_V%Qmu_ref( 181 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 181 ) = 10355.6900000000
  Mref_V%vsh_ref( 181 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 181 ) = 1.00000000000000

  Mref_V%radius_ref( 182 ) = 1234250.00000000
  Mref_V%density_ref( 182 ) = 12159.7700000000
  Mref_V%vpv_ref( 182 ) = 10348.2800000000
  Mref_V%vsv_ref( 182 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 182 ) = 57822.5000000000
  Mref_V%Qmu_ref( 182 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 182 ) = 10348.2800000000
  Mref_V%vsh_ref( 182 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 182 ) = 1.00000000000000

  Mref_V%radius_ref( 183 ) = 1247010.00000000
  Mref_V%density_ref( 183 ) = 12153.1400000000
  Mref_V%vpv_ref( 183 ) = 10340.8500000000
  Mref_V%vsv_ref( 183 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 183 ) = 57822.5000000000
  Mref_V%Qmu_ref( 183 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 183 ) = 10340.8500000000
  Mref_V%vsh_ref( 183 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 183 ) = 1.00000000000000

  Mref_V%radius_ref( 184 ) = 1259770.00000000
  Mref_V%density_ref( 184 ) = 12146.4500000000
  Mref_V%vpv_ref( 184 ) = 10333.3900000000
  Mref_V%vsv_ref( 184 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 184 ) = 57822.5000000000
  Mref_V%Qmu_ref( 184 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 184 ) = 10333.3900000000
  Mref_V%vsh_ref( 184 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 184 ) = 1.00000000000000

  Mref_V%radius_ref( 185 ) = 1272530.00000000
  Mref_V%density_ref( 185 ) = 12139.7100000000
  Mref_V%vpv_ref( 185 ) = 10325.9100000000
  Mref_V%vsv_ref( 185 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 185 ) = 57822.5000000000
  Mref_V%Qmu_ref( 185 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 185 ) = 10325.9100000000
  Mref_V%vsh_ref( 185 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 185 ) = 1.00000000000000

  Mref_V%radius_ref( 186 ) = 1285289.00000000
  Mref_V%density_ref( 186 ) = 12132.9100000000
  Mref_V%vpv_ref( 186 ) = 10318.4000000000
  Mref_V%vsv_ref( 186 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 186 ) = 57822.5000000000
  Mref_V%Qmu_ref( 186 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 186 ) = 10318.4000000000
  Mref_V%vsh_ref( 186 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 186 ) = 1.00000000000000

  Mref_V%radius_ref( 187 ) = 1298049.00000000
  Mref_V%density_ref( 187 ) = 12126.0500000000
  Mref_V%vpv_ref( 187 ) = 10310.8700000000
  Mref_V%vsv_ref( 187 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 187 ) = 57822.5000000000
  Mref_V%Qmu_ref( 187 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 187 ) = 10310.8700000000
  Mref_V%vsh_ref( 187 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 187 ) = 1.00000000000000

  Mref_V%radius_ref( 188 ) = 1310809.00000000
  Mref_V%density_ref( 188 ) = 12119.1400000000
  Mref_V%vpv_ref( 188 ) = 10303.3000000000
  Mref_V%vsv_ref( 188 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 188 ) = 57822.5000000000
  Mref_V%Qmu_ref( 188 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 188 ) = 10303.3000000000
  Mref_V%vsh_ref( 188 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 188 ) = 1.00000000000000

  Mref_V%radius_ref( 189 ) = 1323568.00000000
  Mref_V%density_ref( 189 ) = 12112.1800000000
  Mref_V%vpv_ref( 189 ) = 10295.7100000000
  Mref_V%vsv_ref( 189 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 189 ) = 57822.5000000000
  Mref_V%Qmu_ref( 189 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 189 ) = 10295.7100000000
  Mref_V%vsh_ref( 189 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 189 ) = 1.00000000000000

  Mref_V%radius_ref( 190 ) = 1336328.00000000
  Mref_V%density_ref( 190 ) = 12105.1500000000
  Mref_V%vpv_ref( 190 ) = 10288.0900000000
  Mref_V%vsv_ref( 190 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 190 ) = 57822.5000000000
  Mref_V%Qmu_ref( 190 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 190 ) = 10288.0900000000
  Mref_V%vsh_ref( 190 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 190 ) = 1.00000000000000

  Mref_V%radius_ref( 191 ) = 1349088.00000000
  Mref_V%density_ref( 191 ) = 12098.0700000000
  Mref_V%vpv_ref( 191 ) = 10280.4400000000
  Mref_V%vsv_ref( 191 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 191 ) = 57822.5000000000
  Mref_V%Qmu_ref( 191 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 191 ) = 10280.4400000000
  Mref_V%vsh_ref( 191 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 191 ) = 1.00000000000000

  Mref_V%radius_ref( 192 ) = 1361847.00000000
  Mref_V%density_ref( 192 ) = 12090.9300000000
  Mref_V%vpv_ref( 192 ) = 10272.7600000000
  Mref_V%vsv_ref( 192 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 192 ) = 57822.5000000000
  Mref_V%Qmu_ref( 192 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 192 ) = 10272.7600000000
  Mref_V%vsh_ref( 192 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 192 ) = 1.00000000000000

  Mref_V%radius_ref( 193 ) = 1374607.00000000
  Mref_V%density_ref( 193 ) = 12083.7300000000
  Mref_V%vpv_ref( 193 ) = 10265.0400000000
  Mref_V%vsv_ref( 193 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 193 ) = 57822.5000000000
  Mref_V%Qmu_ref( 193 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 193 ) = 10265.0400000000
  Mref_V%vsh_ref( 193 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 193 ) = 1.00000000000000

  Mref_V%radius_ref( 194 ) = 1387367.00000000
  Mref_V%density_ref( 194 ) = 12076.4800000000
  Mref_V%vpv_ref( 194 ) = 10257.3000000000
  Mref_V%vsv_ref( 194 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 194 ) = 57822.5000000000
  Mref_V%Qmu_ref( 194 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 194 ) = 10257.3000000000
  Mref_V%vsh_ref( 194 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 194 ) = 1.00000000000000

  Mref_V%radius_ref( 195 ) = 1400127.00000000
  Mref_V%density_ref( 195 ) = 12069.1700000000
  Mref_V%vpv_ref( 195 ) = 10249.5200000000
  Mref_V%vsv_ref( 195 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 195 ) = 57822.5000000000
  Mref_V%Qmu_ref( 195 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 195 ) = 10249.5200000000
  Mref_V%vsh_ref( 195 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 195 ) = 1.00000000000000

  Mref_V%radius_ref( 196 ) = 1412886.00000000
  Mref_V%density_ref( 196 ) = 12061.8000000000
  Mref_V%vpv_ref( 196 ) = 10241.7100000000
  Mref_V%vsv_ref( 196 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 196 ) = 57822.5000000000
  Mref_V%Qmu_ref( 196 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 196 ) = 10241.7100000000
  Mref_V%vsh_ref( 196 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 196 ) = 1.00000000000000

  Mref_V%radius_ref( 197 ) = 1425646.00000000
  Mref_V%density_ref( 197 ) = 12054.3700000000
  Mref_V%vpv_ref( 197 ) = 10233.8600000000
  Mref_V%vsv_ref( 197 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 197 ) = 57822.5000000000
  Mref_V%Qmu_ref( 197 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 197 ) = 10233.8600000000
  Mref_V%vsh_ref( 197 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 197 ) = 1.00000000000000

  Mref_V%radius_ref( 198 ) = 1438406.00000000
  Mref_V%density_ref( 198 ) = 12046.8800000000
  Mref_V%vpv_ref( 198 ) = 10225.9800000000
  Mref_V%vsv_ref( 198 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 198 ) = 57822.5000000000
  Mref_V%Qmu_ref( 198 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 198 ) = 10225.9800000000
  Mref_V%vsh_ref( 198 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 198 ) = 1.00000000000000

  Mref_V%radius_ref( 199 ) = 1451165.00000000
  Mref_V%density_ref( 199 ) = 12039.3300000000
  Mref_V%vpv_ref( 199 ) = 10218.0600000000
  Mref_V%vsv_ref( 199 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 199 ) = 57822.5000000000
  Mref_V%Qmu_ref( 199 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 199 ) = 10218.0600000000
  Mref_V%vsh_ref( 199 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 199 ) = 1.00000000000000

  Mref_V%radius_ref( 200 ) = 1463925.00000000
  Mref_V%density_ref( 200 ) = 12031.7200000000
  Mref_V%vpv_ref( 200 ) = 10210.1100000000
  Mref_V%vsv_ref( 200 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 200 ) = 57822.5000000000
  Mref_V%Qmu_ref( 200 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 200 ) = 10210.1100000000
  Mref_V%vsh_ref( 200 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 200 ) = 1.00000000000000

  Mref_V%radius_ref( 201 ) = 1476685.00000000
  Mref_V%density_ref( 201 ) = 12024.0500000000
  Mref_V%vpv_ref( 201 ) = 10202.1200000000
  Mref_V%vsv_ref( 201 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 201 ) = 57822.5000000000
  Mref_V%Qmu_ref( 201 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 201 ) = 10202.1200000000
  Mref_V%vsh_ref( 201 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 201 ) = 1.00000000000000

  Mref_V%radius_ref( 202 ) = 1489444.00000000
  Mref_V%density_ref( 202 ) = 12016.3300000000
  Mref_V%vpv_ref( 202 ) = 10194.1000000000
  Mref_V%vsv_ref( 202 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 202 ) = 57822.5000000000
  Mref_V%Qmu_ref( 202 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 202 ) = 10194.1000000000
  Mref_V%vsh_ref( 202 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 202 ) = 1.00000000000000

  Mref_V%radius_ref( 203 ) = 1502204.00000000
  Mref_V%density_ref( 203 ) = 12008.5400000000
  Mref_V%vpv_ref( 203 ) = 10186.0400000000
  Mref_V%vsv_ref( 203 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 203 ) = 57822.5000000000
  Mref_V%Qmu_ref( 203 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 203 ) = 10186.0400000000
  Mref_V%vsh_ref( 203 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 203 ) = 1.00000000000000

  Mref_V%radius_ref( 204 ) = 1514964.00000000
  Mref_V%density_ref( 204 ) = 12000.6900000000
  Mref_V%vpv_ref( 204 ) = 10177.9400000000
  Mref_V%vsv_ref( 204 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 204 ) = 57822.5000000000
  Mref_V%Qmu_ref( 204 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 204 ) = 10177.9400000000
  Mref_V%vsh_ref( 204 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 204 ) = 1.00000000000000

  Mref_V%radius_ref( 205 ) = 1527724.00000000
  Mref_V%density_ref( 205 ) = 11992.7800000000
  Mref_V%vpv_ref( 205 ) = 10169.7900000000
  Mref_V%vsv_ref( 205 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 205 ) = 57822.5000000000
  Mref_V%Qmu_ref( 205 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 205 ) = 10169.7900000000
  Mref_V%vsh_ref( 205 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 205 ) = 1.00000000000000

  Mref_V%radius_ref( 206 ) = 1540483.00000000
  Mref_V%density_ref( 206 ) = 11984.8100000000
  Mref_V%vpv_ref( 206 ) = 10161.6100000000
  Mref_V%vsv_ref( 206 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 206 ) = 57822.5000000000
  Mref_V%Qmu_ref( 206 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 206 ) = 10161.6100000000
  Mref_V%vsh_ref( 206 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 206 ) = 1.00000000000000

  Mref_V%radius_ref( 207 ) = 1553243.00000000
  Mref_V%density_ref( 207 ) = 11976.7800000000
  Mref_V%vpv_ref( 207 ) = 10153.3900000000
  Mref_V%vsv_ref( 207 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 207 ) = 57822.5000000000
  Mref_V%Qmu_ref( 207 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 207 ) = 10153.3900000000
  Mref_V%vsh_ref( 207 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 207 ) = 1.00000000000000

  Mref_V%radius_ref( 208 ) = 1566003.00000000
  Mref_V%density_ref( 208 ) = 11968.6800000000
  Mref_V%vpv_ref( 208 ) = 10145.1300000000
  Mref_V%vsv_ref( 208 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 208 ) = 57822.5000000000
  Mref_V%Qmu_ref( 208 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 208 ) = 10145.1300000000
  Mref_V%vsh_ref( 208 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 208 ) = 1.00000000000000

  Mref_V%radius_ref( 209 ) = 1578762.00000000
  Mref_V%density_ref( 209 ) = 11960.5300000000
  Mref_V%vpv_ref( 209 ) = 10136.8300000000
  Mref_V%vsv_ref( 209 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 209 ) = 57822.5000000000
  Mref_V%Qmu_ref( 209 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 209 ) = 10136.8300000000
  Mref_V%vsh_ref( 209 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 209 ) = 1.00000000000000

  Mref_V%radius_ref( 210 ) = 1591522.00000000
  Mref_V%density_ref( 210 ) = 11952.3100000000
  Mref_V%vpv_ref( 210 ) = 10128.4800000000
  Mref_V%vsv_ref( 210 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 210 ) = 57822.5000000000
  Mref_V%Qmu_ref( 210 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 210 ) = 10128.4800000000
  Mref_V%vsh_ref( 210 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 210 ) = 1.00000000000000

  Mref_V%radius_ref( 211 ) = 1604282.00000000
  Mref_V%density_ref( 211 ) = 11944.0300000000
  Mref_V%vpv_ref( 211 ) = 10120.0900000000
  Mref_V%vsv_ref( 211 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 211 ) = 57822.5000000000
  Mref_V%Qmu_ref( 211 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 211 ) = 10120.0900000000
  Mref_V%vsh_ref( 211 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 211 ) = 1.00000000000000

  Mref_V%radius_ref( 212 ) = 1617041.00000000
  Mref_V%density_ref( 212 ) = 11935.6900000000
  Mref_V%vpv_ref( 212 ) = 10111.6600000000
  Mref_V%vsv_ref( 212 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 212 ) = 57822.5000000000
  Mref_V%Qmu_ref( 212 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 212 ) = 10111.6600000000
  Mref_V%vsh_ref( 212 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 212 ) = 1.00000000000000

  Mref_V%radius_ref( 213 ) = 1629801.00000000
  Mref_V%density_ref( 213 ) = 11927.2800000000
  Mref_V%vpv_ref( 213 ) = 10103.1800000000
  Mref_V%vsv_ref( 213 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 213 ) = 57822.5000000000
  Mref_V%Qmu_ref( 213 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 213 ) = 10103.1800000000
  Mref_V%vsh_ref( 213 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 213 ) = 1.00000000000000

  Mref_V%radius_ref( 214 ) = 1642561.00000000
  Mref_V%density_ref( 214 ) = 11918.8100000000
  Mref_V%vpv_ref( 214 ) = 10094.6600000000
  Mref_V%vsv_ref( 214 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 214 ) = 57822.5000000000
  Mref_V%Qmu_ref( 214 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 214 ) = 10094.6600000000
  Mref_V%vsh_ref( 214 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 214 ) = 1.00000000000000

  Mref_V%radius_ref( 215 ) = 1655321.00000000
  Mref_V%density_ref( 215 ) = 11910.2800000000
  Mref_V%vpv_ref( 215 ) = 10086.0900000000
  Mref_V%vsv_ref( 215 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 215 ) = 57822.5000000000
  Mref_V%Qmu_ref( 215 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 215 ) = 10086.0900000000
  Mref_V%vsh_ref( 215 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 215 ) = 1.00000000000000

  Mref_V%radius_ref( 216 ) = 1668080.00000000
  Mref_V%density_ref( 216 ) = 11901.6800000000
  Mref_V%vpv_ref( 216 ) = 10077.4800000000
  Mref_V%vsv_ref( 216 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 216 ) = 57822.5000000000
  Mref_V%Qmu_ref( 216 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 216 ) = 10077.4800000000
  Mref_V%vsh_ref( 216 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 216 ) = 1.00000000000000

  Mref_V%radius_ref( 217 ) = 1680840.00000000
  Mref_V%density_ref( 217 ) = 11893.0200000000
  Mref_V%vpv_ref( 217 ) = 10068.8200000000
  Mref_V%vsv_ref( 217 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 217 ) = 57822.5000000000
  Mref_V%Qmu_ref( 217 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 217 ) = 10068.8200000000
  Mref_V%vsh_ref( 217 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 217 ) = 1.00000000000000

  Mref_V%radius_ref( 218 ) = 1693600.00000000
  Mref_V%density_ref( 218 ) = 11884.3000000000
  Mref_V%vpv_ref( 218 ) = 10060.1100000000
  Mref_V%vsv_ref( 218 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 218 ) = 57822.5000000000
  Mref_V%Qmu_ref( 218 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 218 ) = 10060.1100000000
  Mref_V%vsh_ref( 218 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 218 ) = 1.00000000000000

  Mref_V%radius_ref( 219 ) = 1706359.00000000
  Mref_V%density_ref( 219 ) = 11875.5100000000
  Mref_V%vpv_ref( 219 ) = 10051.3500000000
  Mref_V%vsv_ref( 219 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 219 ) = 57822.5000000000
  Mref_V%Qmu_ref( 219 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 219 ) = 10051.3500000000
  Mref_V%vsh_ref( 219 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 219 ) = 1.00000000000000

  Mref_V%radius_ref( 220 ) = 1719119.00000000
  Mref_V%density_ref( 220 ) = 11866.6600000000
  Mref_V%vpv_ref( 220 ) = 10042.5400000000
  Mref_V%vsv_ref( 220 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 220 ) = 57822.5000000000
  Mref_V%Qmu_ref( 220 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 220 ) = 10042.5400000000
  Mref_V%vsh_ref( 220 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 220 ) = 1.00000000000000

  Mref_V%radius_ref( 221 ) = 1731879.00000000
  Mref_V%density_ref( 221 ) = 11857.7400000000
  Mref_V%vpv_ref( 221 ) = 10033.6900000000
  Mref_V%vsv_ref( 221 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 221 ) = 57822.5000000000
  Mref_V%Qmu_ref( 221 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 221 ) = 10033.6900000000
  Mref_V%vsh_ref( 221 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 221 ) = 1.00000000000000

  Mref_V%radius_ref( 222 ) = 1744638.00000000
  Mref_V%density_ref( 222 ) = 11848.7500000000
  Mref_V%vpv_ref( 222 ) = 10024.7800000000
  Mref_V%vsv_ref( 222 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 222 ) = 57822.5000000000
  Mref_V%Qmu_ref( 222 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 222 ) = 10024.7800000000
  Mref_V%vsh_ref( 222 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 222 ) = 1.00000000000000

  Mref_V%radius_ref( 223 ) = 1757398.00000000
  Mref_V%density_ref( 223 ) = 11839.7000000000
  Mref_V%vpv_ref( 223 ) = 10015.8200000000
  Mref_V%vsv_ref( 223 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 223 ) = 57822.5000000000
  Mref_V%Qmu_ref( 223 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 223 ) = 10015.8200000000
  Mref_V%vsh_ref( 223 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 223 ) = 1.00000000000000

  Mref_V%radius_ref( 224 ) = 1770158.00000000
  Mref_V%density_ref( 224 ) = 11830.5800000000
  Mref_V%vpv_ref( 224 ) = 10006.8200000000
  Mref_V%vsv_ref( 224 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 224 ) = 57822.5000000000
  Mref_V%Qmu_ref( 224 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 224 ) = 10006.8200000000
  Mref_V%vsh_ref( 224 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 224 ) = 1.00000000000000

  Mref_V%radius_ref( 225 ) = 1782918.00000000
  Mref_V%density_ref( 225 ) = 11821.4000000000
  Mref_V%vpv_ref( 225 ) = 9997.75000000000
  Mref_V%vsv_ref( 225 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 225 ) = 57822.5000000000
  Mref_V%Qmu_ref( 225 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 225 ) = 9997.75000000000
  Mref_V%vsh_ref( 225 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 225 ) = 1.00000000000000

  Mref_V%radius_ref( 226 ) = 1795677.00000000
  Mref_V%density_ref( 226 ) = 11812.1500000000
  Mref_V%vpv_ref( 226 ) = 9988.64000000000
  Mref_V%vsv_ref( 226 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 226 ) = 57822.5000000000
  Mref_V%Qmu_ref( 226 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 226 ) = 9988.64000000000
  Mref_V%vsh_ref( 226 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 226 ) = 1.00000000000000

  Mref_V%radius_ref( 227 ) = 1808437.00000000
  Mref_V%density_ref( 227 ) = 11802.8400000000
  Mref_V%vpv_ref( 227 ) = 9979.47000000000
  Mref_V%vsv_ref( 227 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 227 ) = 57822.5000000000
  Mref_V%Qmu_ref( 227 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 227 ) = 9979.47000000000
  Mref_V%vsh_ref( 227 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 227 ) = 1.00000000000000

  Mref_V%radius_ref( 228 ) = 1821197.00000000
  Mref_V%density_ref( 228 ) = 11793.4500000000
  Mref_V%vpv_ref( 228 ) = 9970.25000000000
  Mref_V%vsv_ref( 228 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 228 ) = 57822.5000000000
  Mref_V%Qmu_ref( 228 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 228 ) = 9970.25000000000
  Mref_V%vsh_ref( 228 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 228 ) = 1.00000000000000

  Mref_V%radius_ref( 229 ) = 1833956.00000000
  Mref_V%density_ref( 229 ) = 11784.0100000000
  Mref_V%vpv_ref( 229 ) = 9960.97000000000
  Mref_V%vsv_ref( 229 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 229 ) = 57822.5000000000
  Mref_V%Qmu_ref( 229 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 229 ) = 9960.97000000000
  Mref_V%vsh_ref( 229 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 229 ) = 1.00000000000000

  Mref_V%radius_ref( 230 ) = 1846716.00000000
  Mref_V%density_ref( 230 ) = 11774.4900000000
  Mref_V%vpv_ref( 230 ) = 9951.64000000000
  Mref_V%vsv_ref( 230 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 230 ) = 57822.5000000000
  Mref_V%Qmu_ref( 230 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 230 ) = 9951.64000000000
  Mref_V%vsh_ref( 230 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 230 ) = 1.00000000000000

  Mref_V%radius_ref( 231 ) = 1859476.00000000
  Mref_V%density_ref( 231 ) = 11764.9000000000
  Mref_V%vpv_ref( 231 ) = 9942.25000000000
  Mref_V%vsv_ref( 231 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 231 ) = 57822.5000000000
  Mref_V%Qmu_ref( 231 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 231 ) = 9942.25000000000
  Mref_V%vsh_ref( 231 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 231 ) = 1.00000000000000

  Mref_V%radius_ref( 232 ) = 1872235.00000000
  Mref_V%density_ref( 232 ) = 11755.2500000000
  Mref_V%vpv_ref( 232 ) = 9932.81000000000
  Mref_V%vsv_ref( 232 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 232 ) = 57822.5000000000
  Mref_V%Qmu_ref( 232 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 232 ) = 9932.81000000000
  Mref_V%vsh_ref( 232 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 232 ) = 1.00000000000000

  Mref_V%radius_ref( 233 ) = 1884995.00000000
  Mref_V%density_ref( 233 ) = 11745.5300000000
  Mref_V%vpv_ref( 233 ) = 9923.31000000000
  Mref_V%vsv_ref( 233 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 233 ) = 57822.5000000000
  Mref_V%Qmu_ref( 233 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 233 ) = 9923.31000000000
  Mref_V%vsh_ref( 233 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 233 ) = 1.00000000000000

  Mref_V%radius_ref( 234 ) = 1897755.00000000
  Mref_V%density_ref( 234 ) = 11735.7400000000
  Mref_V%vpv_ref( 234 ) = 9913.75000000000
  Mref_V%vsv_ref( 234 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 234 ) = 57822.5000000000
  Mref_V%Qmu_ref( 234 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 234 ) = 9913.75000000000
  Mref_V%vsh_ref( 234 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 234 ) = 1.00000000000000

  Mref_V%radius_ref( 235 ) = 1910515.00000000
  Mref_V%density_ref( 235 ) = 11725.8800000000
  Mref_V%vpv_ref( 235 ) = 9904.13000000000
  Mref_V%vsv_ref( 235 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 235 ) = 57822.5000000000
  Mref_V%Qmu_ref( 235 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 235 ) = 9904.13000000000
  Mref_V%vsh_ref( 235 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 235 ) = 1.00000000000000

  Mref_V%radius_ref( 236 ) = 1923274.00000000
  Mref_V%density_ref( 236 ) = 11715.9500000000
  Mref_V%vpv_ref( 236 ) = 9894.45000000000
  Mref_V%vsv_ref( 236 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 236 ) = 57822.5000000000
  Mref_V%Qmu_ref( 236 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 236 ) = 9894.45000000000
  Mref_V%vsh_ref( 236 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 236 ) = 1.00000000000000

  Mref_V%radius_ref( 237 ) = 1936034.00000000
  Mref_V%density_ref( 237 ) = 11705.9500000000
  Mref_V%vpv_ref( 237 ) = 9884.71000000000
  Mref_V%vsv_ref( 237 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 237 ) = 57822.5000000000
  Mref_V%Qmu_ref( 237 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 237 ) = 9884.71000000000
  Mref_V%vsh_ref( 237 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 237 ) = 1.00000000000000

  Mref_V%radius_ref( 238 ) = 1948794.00000000
  Mref_V%density_ref( 238 ) = 11695.8900000000
  Mref_V%vpv_ref( 238 ) = 9874.91000000000
  Mref_V%vsv_ref( 238 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 238 ) = 57822.5000000000
  Mref_V%Qmu_ref( 238 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 238 ) = 9874.91000000000
  Mref_V%vsh_ref( 238 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 238 ) = 1.00000000000000

  Mref_V%radius_ref( 239 ) = 1961553.00000000
  Mref_V%density_ref( 239 ) = 11685.7500000000
  Mref_V%vpv_ref( 239 ) = 9865.05000000000
  Mref_V%vsv_ref( 239 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 239 ) = 57822.5000000000
  Mref_V%Qmu_ref( 239 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 239 ) = 9865.05000000000
  Mref_V%vsh_ref( 239 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 239 ) = 1.00000000000000

  Mref_V%radius_ref( 240 ) = 1974313.00000000
  Mref_V%density_ref( 240 ) = 11675.5400000000
  Mref_V%vpv_ref( 240 ) = 9855.13000000000
  Mref_V%vsv_ref( 240 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 240 ) = 57822.5000000000
  Mref_V%Qmu_ref( 240 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 240 ) = 9855.13000000000
  Mref_V%vsh_ref( 240 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 240 ) = 1.00000000000000

  Mref_V%radius_ref( 241 ) = 1987073.00000000
  Mref_V%density_ref( 241 ) = 11665.2600000000
  Mref_V%vpv_ref( 241 ) = 9845.14000000000
  Mref_V%vsv_ref( 241 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 241 ) = 57822.5000000000
  Mref_V%Qmu_ref( 241 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 241 ) = 9845.14000000000
  Mref_V%vsh_ref( 241 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 241 ) = 1.00000000000000

  Mref_V%radius_ref( 242 ) = 1999832.00000000
  Mref_V%density_ref( 242 ) = 11654.9200000000
  Mref_V%vpv_ref( 242 ) = 9835.09000000000
  Mref_V%vsv_ref( 242 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 242 ) = 57822.5000000000
  Mref_V%Qmu_ref( 242 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 242 ) = 9835.09000000000
  Mref_V%vsh_ref( 242 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 242 ) = 1.00000000000000

  Mref_V%radius_ref( 243 ) = 2012592.00000000
  Mref_V%density_ref( 243 ) = 11644.5000000000
  Mref_V%vpv_ref( 243 ) = 9824.98000000000
  Mref_V%vsv_ref( 243 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 243 ) = 57822.5000000000
  Mref_V%Qmu_ref( 243 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 243 ) = 9824.98000000000
  Mref_V%vsh_ref( 243 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 243 ) = 1.00000000000000

  Mref_V%radius_ref( 244 ) = 2025352.00000000
  Mref_V%density_ref( 244 ) = 11634.0100000000
  Mref_V%vpv_ref( 244 ) = 9814.80000000000
  Mref_V%vsv_ref( 244 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 244 ) = 57822.5000000000
  Mref_V%Qmu_ref( 244 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 244 ) = 9814.80000000000
  Mref_V%vsh_ref( 244 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 244 ) = 1.00000000000000

  Mref_V%radius_ref( 245 ) = 2038112.00000000
  Mref_V%density_ref( 245 ) = 11623.4400000000
  Mref_V%vpv_ref( 245 ) = 9804.56000000000
  Mref_V%vsv_ref( 245 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 245 ) = 57822.5000000000
  Mref_V%Qmu_ref( 245 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 245 ) = 9804.56000000000
  Mref_V%vsh_ref( 245 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 245 ) = 1.00000000000000

  Mref_V%radius_ref( 246 ) = 2050871.00000000
  Mref_V%density_ref( 246 ) = 11612.8100000000
  Mref_V%vpv_ref( 246 ) = 9794.25000000000
  Mref_V%vsv_ref( 246 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 246 ) = 57822.5000000000
  Mref_V%Qmu_ref( 246 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 246 ) = 9794.25000000000
  Mref_V%vsh_ref( 246 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 246 ) = 1.00000000000000

  Mref_V%radius_ref( 247 ) = 2063631.00000000
  Mref_V%density_ref( 247 ) = 11602.1000000000
  Mref_V%vpv_ref( 247 ) = 9783.87000000000
  Mref_V%vsv_ref( 247 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 247 ) = 57822.5000000000
  Mref_V%Qmu_ref( 247 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 247 ) = 9783.87000000000
  Mref_V%vsh_ref( 247 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 247 ) = 1.00000000000000

  Mref_V%radius_ref( 248 ) = 2076391.00000000
  Mref_V%density_ref( 248 ) = 11591.3200000000
  Mref_V%vpv_ref( 248 ) = 9773.43000000000
  Mref_V%vsv_ref( 248 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 248 ) = 57822.5000000000
  Mref_V%Qmu_ref( 248 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 248 ) = 9773.43000000000
  Mref_V%vsh_ref( 248 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 248 ) = 1.00000000000000

  Mref_V%radius_ref( 249 ) = 2089150.00000000
  Mref_V%density_ref( 249 ) = 11580.4700000000
  Mref_V%vpv_ref( 249 ) = 9762.92000000000
  Mref_V%vsv_ref( 249 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 249 ) = 57822.5000000000
  Mref_V%Qmu_ref( 249 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 249 ) = 9762.92000000000
  Mref_V%vsh_ref( 249 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 249 ) = 1.00000000000000

  Mref_V%radius_ref( 250 ) = 2101910.00000000
  Mref_V%density_ref( 250 ) = 11569.5500000000
  Mref_V%vpv_ref( 250 ) = 9752.34000000000
  Mref_V%vsv_ref( 250 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 250 ) = 57822.5000000000
  Mref_V%Qmu_ref( 250 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 250 ) = 9752.34000000000
  Mref_V%vsh_ref( 250 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 250 ) = 1.00000000000000

  Mref_V%radius_ref( 251 ) = 2114670.00000000
  Mref_V%density_ref( 251 ) = 11558.5500000000
  Mref_V%vpv_ref( 251 ) = 9741.69000000000
  Mref_V%vsv_ref( 251 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 251 ) = 57822.5000000000
  Mref_V%Qmu_ref( 251 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 251 ) = 9741.69000000000
  Mref_V%vsh_ref( 251 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 251 ) = 1.00000000000000

  Mref_V%radius_ref( 252 ) = 2127429.00000000
  Mref_V%density_ref( 252 ) = 11547.4800000000
  Mref_V%vpv_ref( 252 ) = 9730.97000000000
  Mref_V%vsv_ref( 252 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 252 ) = 57822.5000000000
  Mref_V%Qmu_ref( 252 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 252 ) = 9730.97000000000
  Mref_V%vsh_ref( 252 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 252 ) = 1.00000000000000

  Mref_V%radius_ref( 253 ) = 2140189.00000000
  Mref_V%density_ref( 253 ) = 11536.3400000000
  Mref_V%vpv_ref( 253 ) = 9720.18000000000
  Mref_V%vsv_ref( 253 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 253 ) = 57822.5000000000
  Mref_V%Qmu_ref( 253 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 253 ) = 9720.18000000000
  Mref_V%vsh_ref( 253 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 253 ) = 1.00000000000000

  Mref_V%radius_ref( 254 ) = 2152949.00000000
  Mref_V%density_ref( 254 ) = 11525.1200000000
  Mref_V%vpv_ref( 254 ) = 9709.32000000000
  Mref_V%vsv_ref( 254 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 254 ) = 57822.5000000000
  Mref_V%Qmu_ref( 254 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 254 ) = 9709.32000000000
  Mref_V%vsh_ref( 254 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 254 ) = 1.00000000000000

  Mref_V%radius_ref( 255 ) = 2165709.00000000
  Mref_V%density_ref( 255 ) = 11513.8300000000
  Mref_V%vpv_ref( 255 ) = 9698.39000000000
  Mref_V%vsv_ref( 255 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 255 ) = 57822.5000000000
  Mref_V%Qmu_ref( 255 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 255 ) = 9698.39000000000
  Mref_V%vsh_ref( 255 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 255 ) = 1.00000000000000

  Mref_V%radius_ref( 256 ) = 2178468.00000000
  Mref_V%density_ref( 256 ) = 11502.4600000000
  Mref_V%vpv_ref( 256 ) = 9687.38000000000
  Mref_V%vsv_ref( 256 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 256 ) = 57822.5000000000
  Mref_V%Qmu_ref( 256 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 256 ) = 9687.38000000000
  Mref_V%vsh_ref( 256 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 256 ) = 1.00000000000000

  Mref_V%radius_ref( 257 ) = 2191228.00000000
  Mref_V%density_ref( 257 ) = 11491.0200000000
  Mref_V%vpv_ref( 257 ) = 9676.31000000000
  Mref_V%vsv_ref( 257 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 257 ) = 57822.5000000000
  Mref_V%Qmu_ref( 257 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 257 ) = 9676.31000000000
  Mref_V%vsh_ref( 257 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 257 ) = 1.00000000000000

  Mref_V%radius_ref( 258 ) = 2203988.00000000
  Mref_V%density_ref( 258 ) = 11479.5100000000
  Mref_V%vpv_ref( 258 ) = 9665.15000000000
  Mref_V%vsv_ref( 258 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 258 ) = 57822.5000000000
  Mref_V%Qmu_ref( 258 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 258 ) = 9665.15000000000
  Mref_V%vsh_ref( 258 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 258 ) = 1.00000000000000

  Mref_V%radius_ref( 259 ) = 2216747.00000000
  Mref_V%density_ref( 259 ) = 11467.9100000000
  Mref_V%vpv_ref( 259 ) = 9653.93000000000
  Mref_V%vsv_ref( 259 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 259 ) = 57822.5000000000
  Mref_V%Qmu_ref( 259 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 259 ) = 9653.93000000000
  Mref_V%vsh_ref( 259 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 259 ) = 1.00000000000000

  Mref_V%radius_ref( 260 ) = 2229507.00000000
  Mref_V%density_ref( 260 ) = 11456.2500000000
  Mref_V%vpv_ref( 260 ) = 9642.63000000000
  Mref_V%vsv_ref( 260 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 260 ) = 57822.5000000000
  Mref_V%Qmu_ref( 260 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 260 ) = 9642.63000000000
  Mref_V%vsh_ref( 260 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 260 ) = 1.00000000000000

  Mref_V%radius_ref( 261 ) = 2242267.00000000
  Mref_V%density_ref( 261 ) = 11444.5000000000
  Mref_V%vpv_ref( 261 ) = 9631.25000000000
  Mref_V%vsv_ref( 261 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 261 ) = 57822.5000000000
  Mref_V%Qmu_ref( 261 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 261 ) = 9631.25000000000
  Mref_V%vsh_ref( 261 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 261 ) = 1.00000000000000

  Mref_V%radius_ref( 262 ) = 2255026.00000000
  Mref_V%density_ref( 262 ) = 11432.6900000000
  Mref_V%vpv_ref( 262 ) = 9619.80000000000
  Mref_V%vsv_ref( 262 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 262 ) = 57822.5000000000
  Mref_V%Qmu_ref( 262 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 262 ) = 9619.80000000000
  Mref_V%vsh_ref( 262 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 262 ) = 1.00000000000000

  Mref_V%radius_ref( 263 ) = 2267786.00000000
  Mref_V%density_ref( 263 ) = 11420.7900000000
  Mref_V%vpv_ref( 263 ) = 9608.27000000000
  Mref_V%vsv_ref( 263 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 263 ) = 57822.5000000000
  Mref_V%Qmu_ref( 263 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 263 ) = 9608.27000000000
  Mref_V%vsh_ref( 263 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 263 ) = 1.00000000000000

  Mref_V%radius_ref( 264 ) = 2280546.00000000
  Mref_V%density_ref( 264 ) = 11408.8200000000
  Mref_V%vpv_ref( 264 ) = 9596.66000000000
  Mref_V%vsv_ref( 264 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 264 ) = 57822.5000000000
  Mref_V%Qmu_ref( 264 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 264 ) = 9596.66000000000
  Mref_V%vsh_ref( 264 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 264 ) = 1.00000000000000

  Mref_V%radius_ref( 265 ) = 2293306.00000000
  Mref_V%density_ref( 265 ) = 11396.7700000000
  Mref_V%vpv_ref( 265 ) = 9584.97000000000
  Mref_V%vsv_ref( 265 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 265 ) = 57822.5000000000
  Mref_V%Qmu_ref( 265 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 265 ) = 9584.97000000000
  Mref_V%vsh_ref( 265 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 265 ) = 1.00000000000000

  Mref_V%radius_ref( 266 ) = 2306065.00000000
  Mref_V%density_ref( 266 ) = 11384.6400000000
  Mref_V%vpv_ref( 266 ) = 9573.20000000000
  Mref_V%vsv_ref( 266 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 266 ) = 57822.5000000000
  Mref_V%Qmu_ref( 266 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 266 ) = 9573.20000000000
  Mref_V%vsh_ref( 266 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 266 ) = 1.00000000000000

  Mref_V%radius_ref( 267 ) = 2318825.00000000
  Mref_V%density_ref( 267 ) = 11372.4400000000
  Mref_V%vpv_ref( 267 ) = 9561.36000000000
  Mref_V%vsv_ref( 267 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 267 ) = 57822.5000000000
  Mref_V%Qmu_ref( 267 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 267 ) = 9561.36000000000
  Mref_V%vsh_ref( 267 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 267 ) = 1.00000000000000

  Mref_V%radius_ref( 268 ) = 2331585.00000000
  Mref_V%density_ref( 268 ) = 11360.1600000000
  Mref_V%vpv_ref( 268 ) = 9549.43000000000
  Mref_V%vsv_ref( 268 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 268 ) = 57822.5000000000
  Mref_V%Qmu_ref( 268 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 268 ) = 9549.43000000000
  Mref_V%vsh_ref( 268 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 268 ) = 1.00000000000000

  Mref_V%radius_ref( 269 ) = 2344344.00000000
  Mref_V%density_ref( 269 ) = 11347.8000000000
  Mref_V%vpv_ref( 269 ) = 9537.43000000000
  Mref_V%vsv_ref( 269 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 269 ) = 57822.5000000000
  Mref_V%Qmu_ref( 269 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 269 ) = 9537.43000000000
  Mref_V%vsh_ref( 269 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 269 ) = 1.00000000000000

  Mref_V%radius_ref( 270 ) = 2357104.00000000
  Mref_V%density_ref( 270 ) = 11335.3700000000
  Mref_V%vpv_ref( 270 ) = 9525.34000000000
  Mref_V%vsv_ref( 270 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 270 ) = 57822.5000000000
  Mref_V%Qmu_ref( 270 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 270 ) = 9525.34000000000
  Mref_V%vsh_ref( 270 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 270 ) = 1.00000000000000

  Mref_V%radius_ref( 271 ) = 2369864.00000000
  Mref_V%density_ref( 271 ) = 11322.8500000000
  Mref_V%vpv_ref( 271 ) = 9513.17000000000
  Mref_V%vsv_ref( 271 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 271 ) = 57822.5000000000
  Mref_V%Qmu_ref( 271 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 271 ) = 9513.17000000000
  Mref_V%vsh_ref( 271 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 271 ) = 1.00000000000000

  Mref_V%radius_ref( 272 ) = 2382623.00000000
  Mref_V%density_ref( 272 ) = 11310.2600000000
  Mref_V%vpv_ref( 272 ) = 9500.91000000000
  Mref_V%vsv_ref( 272 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 272 ) = 57822.5000000000
  Mref_V%Qmu_ref( 272 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 272 ) = 9500.91000000000
  Mref_V%vsh_ref( 272 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 272 ) = 1.00000000000000

  Mref_V%radius_ref( 273 ) = 2395383.00000000
  Mref_V%density_ref( 273 ) = 11297.5800000000
  Mref_V%vpv_ref( 273 ) = 9488.57000000000
  Mref_V%vsv_ref( 273 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 273 ) = 57822.5000000000
  Mref_V%Qmu_ref( 273 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 273 ) = 9488.57000000000
  Mref_V%vsh_ref( 273 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 273 ) = 1.00000000000000

  Mref_V%radius_ref( 274 ) = 2408143.00000000
  Mref_V%density_ref( 274 ) = 11284.8300000000
  Mref_V%vpv_ref( 274 ) = 9476.15000000000
  Mref_V%vsv_ref( 274 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 274 ) = 57822.5000000000
  Mref_V%Qmu_ref( 274 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 274 ) = 9476.15000000000
  Mref_V%vsh_ref( 274 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 274 ) = 1.00000000000000

  Mref_V%radius_ref( 275 ) = 2420903.00000000
  Mref_V%density_ref( 275 ) = 11272.0000000000
  Mref_V%vpv_ref( 275 ) = 9463.64000000000
  Mref_V%vsv_ref( 275 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 275 ) = 57822.5000000000
  Mref_V%Qmu_ref( 275 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 275 ) = 9463.64000000000
  Mref_V%vsh_ref( 275 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 275 ) = 1.00000000000000

  Mref_V%radius_ref( 276 ) = 2433662.00000000
  Mref_V%density_ref( 276 ) = 11259.0900000000
  Mref_V%vpv_ref( 276 ) = 9451.05000000000
  Mref_V%vsv_ref( 276 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 276 ) = 57822.5000000000
  Mref_V%Qmu_ref( 276 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 276 ) = 9451.05000000000
  Mref_V%vsh_ref( 276 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 276 ) = 1.00000000000000

  Mref_V%radius_ref( 277 ) = 2446422.00000000
  Mref_V%density_ref( 277 ) = 11246.1000000000
  Mref_V%vpv_ref( 277 ) = 9438.37000000000
  Mref_V%vsv_ref( 277 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 277 ) = 57822.5000000000
  Mref_V%Qmu_ref( 277 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 277 ) = 9438.37000000000
  Mref_V%vsh_ref( 277 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 277 ) = 1.00000000000000

  Mref_V%radius_ref( 278 ) = 2459182.00000000
  Mref_V%density_ref( 278 ) = 11233.0300000000
  Mref_V%vpv_ref( 278 ) = 9425.61000000000
  Mref_V%vsv_ref( 278 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 278 ) = 57822.5000000000
  Mref_V%Qmu_ref( 278 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 278 ) = 9425.61000000000
  Mref_V%vsh_ref( 278 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 278 ) = 1.00000000000000

  Mref_V%radius_ref( 279 ) = 2471941.00000000
  Mref_V%density_ref( 279 ) = 11219.8700000000
  Mref_V%vpv_ref( 279 ) = 9412.75000000000
  Mref_V%vsv_ref( 279 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 279 ) = 57822.5000000000
  Mref_V%Qmu_ref( 279 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 279 ) = 9412.75000000000
  Mref_V%vsh_ref( 279 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 279 ) = 1.00000000000000

  Mref_V%radius_ref( 280 ) = 2484701.00000000
  Mref_V%density_ref( 280 ) = 11206.6400000000
  Mref_V%vpv_ref( 280 ) = 9399.81000000000
  Mref_V%vsv_ref( 280 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 280 ) = 57822.5000000000
  Mref_V%Qmu_ref( 280 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 280 ) = 9399.81000000000
  Mref_V%vsh_ref( 280 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 280 ) = 1.00000000000000

  Mref_V%radius_ref( 281 ) = 2497461.00000000
  Mref_V%density_ref( 281 ) = 11193.3300000000
  Mref_V%vpv_ref( 281 ) = 9386.78000000000
  Mref_V%vsv_ref( 281 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 281 ) = 57822.5000000000
  Mref_V%Qmu_ref( 281 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 281 ) = 9386.78000000000
  Mref_V%vsh_ref( 281 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 281 ) = 1.00000000000000

  Mref_V%radius_ref( 282 ) = 2510220.00000000
  Mref_V%density_ref( 282 ) = 11179.9300000000
  Mref_V%vpv_ref( 282 ) = 9373.66000000000
  Mref_V%vsv_ref( 282 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 282 ) = 57822.5000000000
  Mref_V%Qmu_ref( 282 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 282 ) = 9373.66000000000
  Mref_V%vsh_ref( 282 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 282 ) = 1.00000000000000

  Mref_V%radius_ref( 283 ) = 2522980.00000000
  Mref_V%density_ref( 283 ) = 11166.4500000000
  Mref_V%vpv_ref( 283 ) = 9360.45000000000
  Mref_V%vsv_ref( 283 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 283 ) = 57822.5000000000
  Mref_V%Qmu_ref( 283 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 283 ) = 9360.45000000000
  Mref_V%vsh_ref( 283 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 283 ) = 1.00000000000000

  Mref_V%radius_ref( 284 ) = 2535740.00000000
  Mref_V%density_ref( 284 ) = 11152.8900000000
  Mref_V%vpv_ref( 284 ) = 9347.15000000000
  Mref_V%vsv_ref( 284 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 284 ) = 57822.5000000000
  Mref_V%Qmu_ref( 284 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 284 ) = 9347.15000000000
  Mref_V%vsh_ref( 284 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 284 ) = 1.00000000000000

  Mref_V%radius_ref( 285 ) = 2548500.00000000
  Mref_V%density_ref( 285 ) = 11139.2500000000
  Mref_V%vpv_ref( 285 ) = 9333.76000000000
  Mref_V%vsv_ref( 285 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 285 ) = 57822.5000000000
  Mref_V%Qmu_ref( 285 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 285 ) = 9333.76000000000
  Mref_V%vsh_ref( 285 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 285 ) = 1.00000000000000

  Mref_V%radius_ref( 286 ) = 2561259.00000000
  Mref_V%density_ref( 286 ) = 11125.5300000000
  Mref_V%vpv_ref( 286 ) = 9320.27000000000
  Mref_V%vsv_ref( 286 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 286 ) = 57822.5000000000
  Mref_V%Qmu_ref( 286 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 286 ) = 9320.27000000000
  Mref_V%vsh_ref( 286 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 286 ) = 1.00000000000000

  Mref_V%radius_ref( 287 ) = 2574019.00000000
  Mref_V%density_ref( 287 ) = 11111.7200000000
  Mref_V%vpv_ref( 287 ) = 9306.70000000000
  Mref_V%vsv_ref( 287 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 287 ) = 57822.5000000000
  Mref_V%Qmu_ref( 287 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 287 ) = 9306.70000000000
  Mref_V%vsh_ref( 287 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 287 ) = 1.00000000000000

  Mref_V%radius_ref( 288 ) = 2586779.00000000
  Mref_V%density_ref( 288 ) = 11097.8300000000
  Mref_V%vpv_ref( 288 ) = 9293.03000000000
  Mref_V%vsv_ref( 288 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 288 ) = 57822.5000000000
  Mref_V%Qmu_ref( 288 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 288 ) = 9293.03000000000
  Mref_V%vsh_ref( 288 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 288 ) = 1.00000000000000

  Mref_V%radius_ref( 289 ) = 2599538.00000000
  Mref_V%density_ref( 289 ) = 11083.8600000000
  Mref_V%vpv_ref( 289 ) = 9279.26000000000
  Mref_V%vsv_ref( 289 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 289 ) = 57822.5000000000
  Mref_V%Qmu_ref( 289 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 289 ) = 9279.26000000000
  Mref_V%vsh_ref( 289 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 289 ) = 1.00000000000000

  Mref_V%radius_ref( 290 ) = 2612298.00000000
  Mref_V%density_ref( 290 ) = 11069.8000000000
  Mref_V%vpv_ref( 290 ) = 9265.40000000000
  Mref_V%vsv_ref( 290 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 290 ) = 57822.5000000000
  Mref_V%Qmu_ref( 290 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 290 ) = 9265.40000000000
  Mref_V%vsh_ref( 290 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 290 ) = 1.00000000000000

  Mref_V%radius_ref( 291 ) = 2625058.00000000
  Mref_V%density_ref( 291 ) = 11055.6600000000
  Mref_V%vpv_ref( 291 ) = 9251.45000000000
  Mref_V%vsv_ref( 291 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 291 ) = 57822.5000000000
  Mref_V%Qmu_ref( 291 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 291 ) = 9251.45000000000
  Mref_V%vsh_ref( 291 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 291 ) = 1.00000000000000

  Mref_V%radius_ref( 292 ) = 2637818.00000000
  Mref_V%density_ref( 292 ) = 11041.4400000000
  Mref_V%vpv_ref( 292 ) = 9237.40000000000
  Mref_V%vsv_ref( 292 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 292 ) = 57822.5000000000
  Mref_V%Qmu_ref( 292 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 292 ) = 9237.40000000000
  Mref_V%vsh_ref( 292 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 292 ) = 1.00000000000000

  Mref_V%radius_ref( 293 ) = 2650577.00000000
  Mref_V%density_ref( 293 ) = 11027.1300000000
  Mref_V%vpv_ref( 293 ) = 9223.25000000000
  Mref_V%vsv_ref( 293 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 293 ) = 57822.5000000000
  Mref_V%Qmu_ref( 293 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 293 ) = 9223.25000000000
  Mref_V%vsh_ref( 293 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 293 ) = 1.00000000000000

  Mref_V%radius_ref( 294 ) = 2663337.00000000
  Mref_V%density_ref( 294 ) = 11012.7400000000
  Mref_V%vpv_ref( 294 ) = 9209.00000000000
  Mref_V%vsv_ref( 294 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 294 ) = 57822.5000000000
  Mref_V%Qmu_ref( 294 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 294 ) = 9209.00000000000
  Mref_V%vsh_ref( 294 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 294 ) = 1.00000000000000

  Mref_V%radius_ref( 295 ) = 2676097.00000000
  Mref_V%density_ref( 295 ) = 10998.2600000000
  Mref_V%vpv_ref( 295 ) = 9194.66000000000
  Mref_V%vsv_ref( 295 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 295 ) = 57822.5000000000
  Mref_V%Qmu_ref( 295 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 295 ) = 9194.66000000000
  Mref_V%vsh_ref( 295 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 295 ) = 1.00000000000000

  Mref_V%radius_ref( 296 ) = 2688856.00000000
  Mref_V%density_ref( 296 ) = 10983.7000000000
  Mref_V%vpv_ref( 296 ) = 9180.22000000000
  Mref_V%vsv_ref( 296 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 296 ) = 57822.5000000000
  Mref_V%Qmu_ref( 296 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 296 ) = 9180.22000000000
  Mref_V%vsh_ref( 296 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 296 ) = 1.00000000000000

  Mref_V%radius_ref( 297 ) = 2701616.00000000
  Mref_V%density_ref( 297 ) = 10969.0500000000
  Mref_V%vpv_ref( 297 ) = 9165.68000000000
  Mref_V%vsv_ref( 297 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 297 ) = 57822.5000000000
  Mref_V%Qmu_ref( 297 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 297 ) = 9165.68000000000
  Mref_V%vsh_ref( 297 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 297 ) = 1.00000000000000

  Mref_V%radius_ref( 298 ) = 2714376.00000000
  Mref_V%density_ref( 298 ) = 10954.3200000000
  Mref_V%vpv_ref( 298 ) = 9151.03000000000
  Mref_V%vsv_ref( 298 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 298 ) = 57822.5000000000
  Mref_V%Qmu_ref( 298 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 298 ) = 9151.03000000000
  Mref_V%vsh_ref( 298 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 298 ) = 1.00000000000000

  Mref_V%radius_ref( 299 ) = 2727135.00000000
  Mref_V%density_ref( 299 ) = 10939.5000000000
  Mref_V%vpv_ref( 299 ) = 9136.29000000000
  Mref_V%vsv_ref( 299 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 299 ) = 57822.5000000000
  Mref_V%Qmu_ref( 299 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 299 ) = 9136.29000000000
  Mref_V%vsh_ref( 299 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 299 ) = 1.00000000000000

  Mref_V%radius_ref( 300 ) = 2739895.00000000
  Mref_V%density_ref( 300 ) = 10924.5900000000
  Mref_V%vpv_ref( 300 ) = 9121.45000000000
  Mref_V%vsv_ref( 300 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 300 ) = 57822.5000000000
  Mref_V%Qmu_ref( 300 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 300 ) = 9121.45000000000
  Mref_V%vsh_ref( 300 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 300 ) = 1.00000000000000

  Mref_V%radius_ref( 301 ) = 2752655.00000000
  Mref_V%density_ref( 301 ) = 10909.6000000000
  Mref_V%vpv_ref( 301 ) = 9106.50000000000
  Mref_V%vsv_ref( 301 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 301 ) = 57822.5000000000
  Mref_V%Qmu_ref( 301 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 301 ) = 9106.50000000000
  Mref_V%vsh_ref( 301 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 301 ) = 1.00000000000000

  Mref_V%radius_ref( 302 ) = 2765415.00000000
  Mref_V%density_ref( 302 ) = 10894.5200000000
  Mref_V%vpv_ref( 302 ) = 9091.46000000000
  Mref_V%vsv_ref( 302 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 302 ) = 57822.5000000000
  Mref_V%Qmu_ref( 302 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 302 ) = 9091.46000000000
  Mref_V%vsh_ref( 302 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 302 ) = 1.00000000000000

  Mref_V%radius_ref( 303 ) = 2778174.00000000
  Mref_V%density_ref( 303 ) = 10879.3500000000
  Mref_V%vpv_ref( 303 ) = 9076.30000000000
  Mref_V%vsv_ref( 303 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 303 ) = 57822.5000000000
  Mref_V%Qmu_ref( 303 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 303 ) = 9076.30000000000
  Mref_V%vsh_ref( 303 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 303 ) = 1.00000000000000

  Mref_V%radius_ref( 304 ) = 2790934.00000000
  Mref_V%density_ref( 304 ) = 10864.1000000000
  Mref_V%vpv_ref( 304 ) = 9061.05000000000
  Mref_V%vsv_ref( 304 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 304 ) = 57822.5000000000
  Mref_V%Qmu_ref( 304 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 304 ) = 9061.05000000000
  Mref_V%vsh_ref( 304 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 304 ) = 1.00000000000000

  Mref_V%radius_ref( 305 ) = 2803694.00000000
  Mref_V%density_ref( 305 ) = 10848.7600000000
  Mref_V%vpv_ref( 305 ) = 9045.69000000000
  Mref_V%vsv_ref( 305 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 305 ) = 57822.5000000000
  Mref_V%Qmu_ref( 305 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 305 ) = 9045.69000000000
  Mref_V%vsh_ref( 305 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 305 ) = 1.00000000000000

  Mref_V%radius_ref( 306 ) = 2816453.00000000
  Mref_V%density_ref( 306 ) = 10833.3300000000
  Mref_V%vpv_ref( 306 ) = 9030.23000000000
  Mref_V%vsv_ref( 306 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 306 ) = 57822.5000000000
  Mref_V%Qmu_ref( 306 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 306 ) = 9030.23000000000
  Mref_V%vsh_ref( 306 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 306 ) = 1.00000000000000

  Mref_V%radius_ref( 307 ) = 2829213.00000000
  Mref_V%density_ref( 307 ) = 10817.8100000000
  Mref_V%vpv_ref( 307 ) = 9014.65000000000
  Mref_V%vsv_ref( 307 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 307 ) = 57822.5000000000
  Mref_V%Qmu_ref( 307 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 307 ) = 9014.65000000000
  Mref_V%vsh_ref( 307 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 307 ) = 1.00000000000000

  Mref_V%radius_ref( 308 ) = 2841973.00000000
  Mref_V%density_ref( 308 ) = 10802.2100000000
  Mref_V%vpv_ref( 308 ) = 8998.98000000000
  Mref_V%vsv_ref( 308 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 308 ) = 57822.5000000000
  Mref_V%Qmu_ref( 308 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 308 ) = 8998.98000000000
  Mref_V%vsh_ref( 308 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 308 ) = 1.00000000000000

  Mref_V%radius_ref( 309 ) = 2854732.00000000
  Mref_V%density_ref( 309 ) = 10786.5100000000
  Mref_V%vpv_ref( 309 ) = 8983.19000000000
  Mref_V%vsv_ref( 309 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 309 ) = 57822.5000000000
  Mref_V%Qmu_ref( 309 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 309 ) = 8983.19000000000
  Mref_V%vsh_ref( 309 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 309 ) = 1.00000000000000

  Mref_V%radius_ref( 310 ) = 2867492.00000000
  Mref_V%density_ref( 310 ) = 10770.7300000000
  Mref_V%vpv_ref( 310 ) = 8967.30000000000
  Mref_V%vsv_ref( 310 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 310 ) = 57822.5000000000
  Mref_V%Qmu_ref( 310 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 310 ) = 8967.30000000000
  Mref_V%vsh_ref( 310 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 310 ) = 1.00000000000000

  Mref_V%radius_ref( 311 ) = 2880252.00000000
  Mref_V%density_ref( 311 ) = 10754.8600000000
  Mref_V%vpv_ref( 311 ) = 8951.30000000000
  Mref_V%vsv_ref( 311 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 311 ) = 57822.5000000000
  Mref_V%Qmu_ref( 311 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 311 ) = 8951.30000000000
  Mref_V%vsh_ref( 311 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 311 ) = 1.00000000000000

  Mref_V%radius_ref( 312 ) = 2893012.00000000
  Mref_V%density_ref( 312 ) = 10738.9000000000
  Mref_V%vpv_ref( 312 ) = 8935.19000000000
  Mref_V%vsv_ref( 312 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 312 ) = 57822.5000000000
  Mref_V%Qmu_ref( 312 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 312 ) = 8935.19000000000
  Mref_V%vsh_ref( 312 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 312 ) = 1.00000000000000

  Mref_V%radius_ref( 313 ) = 2905771.00000000
  Mref_V%density_ref( 313 ) = 10722.8500000000
  Mref_V%vpv_ref( 313 ) = 8918.97000000000
  Mref_V%vsv_ref( 313 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 313 ) = 57822.5000000000
  Mref_V%Qmu_ref( 313 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 313 ) = 8918.97000000000
  Mref_V%vsh_ref( 313 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 313 ) = 1.00000000000000

  Mref_V%radius_ref( 314 ) = 2918531.00000000
  Mref_V%density_ref( 314 ) = 10706.7100000000
  Mref_V%vpv_ref( 314 ) = 8902.64000000000
  Mref_V%vsv_ref( 314 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 314 ) = 57822.5000000000
  Mref_V%Qmu_ref( 314 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 314 ) = 8902.64000000000
  Mref_V%vsh_ref( 314 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 314 ) = 1.00000000000000

  Mref_V%radius_ref( 315 ) = 2931291.00000000
  Mref_V%density_ref( 315 ) = 10690.4800000000
  Mref_V%vpv_ref( 315 ) = 8886.20000000000
  Mref_V%vsv_ref( 315 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 315 ) = 57822.5000000000
  Mref_V%Qmu_ref( 315 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 315 ) = 8886.20000000000
  Mref_V%vsh_ref( 315 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 315 ) = 1.00000000000000

  Mref_V%radius_ref( 316 ) = 2944050.00000000
  Mref_V%density_ref( 316 ) = 10674.1600000000
  Mref_V%vpv_ref( 316 ) = 8869.64000000000
  Mref_V%vsv_ref( 316 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 316 ) = 57822.5000000000
  Mref_V%Qmu_ref( 316 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 316 ) = 8869.64000000000
  Mref_V%vsh_ref( 316 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 316 ) = 1.00000000000000

  Mref_V%radius_ref( 317 ) = 2956810.00000000
  Mref_V%density_ref( 317 ) = 10657.7500000000
  Mref_V%vpv_ref( 317 ) = 8852.98000000000
  Mref_V%vsv_ref( 317 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 317 ) = 57822.5000000000
  Mref_V%Qmu_ref( 317 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 317 ) = 8852.98000000000
  Mref_V%vsh_ref( 317 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 317 ) = 1.00000000000000

  Mref_V%radius_ref( 318 ) = 2969570.00000000
  Mref_V%density_ref( 318 ) = 10641.2400000000
  Mref_V%vpv_ref( 318 ) = 8836.20000000000
  Mref_V%vsv_ref( 318 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 318 ) = 57822.5000000000
  Mref_V%Qmu_ref( 318 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 318 ) = 8836.20000000000
  Mref_V%vsh_ref( 318 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 318 ) = 1.00000000000000

  Mref_V%radius_ref( 319 ) = 2982329.00000000
  Mref_V%density_ref( 319 ) = 10624.6500000000
  Mref_V%vpv_ref( 319 ) = 8819.31000000000
  Mref_V%vsv_ref( 319 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 319 ) = 57822.5000000000
  Mref_V%Qmu_ref( 319 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 319 ) = 8819.31000000000
  Mref_V%vsh_ref( 319 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 319 ) = 1.00000000000000

  Mref_V%radius_ref( 320 ) = 2995089.00000000
  Mref_V%density_ref( 320 ) = 10607.9600000000
  Mref_V%vpv_ref( 320 ) = 8802.30000000000
  Mref_V%vsv_ref( 320 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 320 ) = 57822.5000000000
  Mref_V%Qmu_ref( 320 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 320 ) = 8802.30000000000
  Mref_V%vsh_ref( 320 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 320 ) = 1.00000000000000

  Mref_V%radius_ref( 321 ) = 3007849.00000000
  Mref_V%density_ref( 321 ) = 10591.1900000000
  Mref_V%vpv_ref( 321 ) = 8785.18000000000
  Mref_V%vsv_ref( 321 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 321 ) = 57822.5000000000
  Mref_V%Qmu_ref( 321 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 321 ) = 8785.18000000000
  Mref_V%vsh_ref( 321 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 321 ) = 1.00000000000000

  Mref_V%radius_ref( 322 ) = 3020609.00000000
  Mref_V%density_ref( 322 ) = 10574.3200000000
  Mref_V%vpv_ref( 322 ) = 8767.94000000000
  Mref_V%vsv_ref( 322 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 322 ) = 57822.5000000000
  Mref_V%Qmu_ref( 322 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 322 ) = 8767.94000000000
  Mref_V%vsh_ref( 322 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 322 ) = 1.00000000000000

  Mref_V%radius_ref( 323 ) = 3033368.00000000
  Mref_V%density_ref( 323 ) = 10557.3600000000
  Mref_V%vpv_ref( 323 ) = 8750.59000000000
  Mref_V%vsv_ref( 323 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 323 ) = 57822.5000000000
  Mref_V%Qmu_ref( 323 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 323 ) = 8750.59000000000
  Mref_V%vsh_ref( 323 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 323 ) = 1.00000000000000

  Mref_V%radius_ref( 324 ) = 3046128.00000000
  Mref_V%density_ref( 324 ) = 10540.3000000000
  Mref_V%vpv_ref( 324 ) = 8733.12000000000
  Mref_V%vsv_ref( 324 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 324 ) = 57822.5000000000
  Mref_V%Qmu_ref( 324 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 324 ) = 8733.12000000000
  Mref_V%vsh_ref( 324 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 324 ) = 1.00000000000000

  Mref_V%radius_ref( 325 ) = 3058888.00000000
  Mref_V%density_ref( 325 ) = 10523.1600000000
  Mref_V%vpv_ref( 325 ) = 8715.53000000000
  Mref_V%vsv_ref( 325 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 325 ) = 57822.5000000000
  Mref_V%Qmu_ref( 325 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 325 ) = 8715.53000000000
  Mref_V%vsh_ref( 325 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 325 ) = 1.00000000000000

  Mref_V%radius_ref( 326 ) = 3071647.00000000
  Mref_V%density_ref( 326 ) = 10505.9200000000
  Mref_V%vpv_ref( 326 ) = 8697.82000000000
  Mref_V%vsv_ref( 326 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 326 ) = 57822.5000000000
  Mref_V%Qmu_ref( 326 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 326 ) = 8697.82000000000
  Mref_V%vsh_ref( 326 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 326 ) = 1.00000000000000

  Mref_V%radius_ref( 327 ) = 3084407.00000000
  Mref_V%density_ref( 327 ) = 10488.5800000000
  Mref_V%vpv_ref( 327 ) = 8680.00000000000
  Mref_V%vsv_ref( 327 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 327 ) = 57822.5000000000
  Mref_V%Qmu_ref( 327 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 327 ) = 8680.00000000000
  Mref_V%vsh_ref( 327 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 327 ) = 1.00000000000000

  Mref_V%radius_ref( 328 ) = 3097167.00000000
  Mref_V%density_ref( 328 ) = 10471.1500000000
  Mref_V%vpv_ref( 328 ) = 8662.05000000000
  Mref_V%vsv_ref( 328 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 328 ) = 57822.5000000000
  Mref_V%Qmu_ref( 328 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 328 ) = 8662.05000000000
  Mref_V%vsh_ref( 328 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 328 ) = 1.00000000000000

  Mref_V%radius_ref( 329 ) = 3109926.00000000
  Mref_V%density_ref( 329 ) = 10453.6300000000
  Mref_V%vpv_ref( 329 ) = 8643.99000000000
  Mref_V%vsv_ref( 329 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 329 ) = 57822.5000000000
  Mref_V%Qmu_ref( 329 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 329 ) = 8643.99000000000
  Mref_V%vsh_ref( 329 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 329 ) = 1.00000000000000

  Mref_V%radius_ref( 330 ) = 3122686.00000000
  Mref_V%density_ref( 330 ) = 10436.0200000000
  Mref_V%vpv_ref( 330 ) = 8625.80000000000
  Mref_V%vsv_ref( 330 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 330 ) = 57822.5000000000
  Mref_V%Qmu_ref( 330 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 330 ) = 8625.80000000000
  Mref_V%vsh_ref( 330 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 330 ) = 1.00000000000000

  Mref_V%radius_ref( 331 ) = 3135446.00000000
  Mref_V%density_ref( 331 ) = 10418.3100000000
  Mref_V%vpv_ref( 331 ) = 8607.49000000000
  Mref_V%vsv_ref( 331 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 331 ) = 57822.5000000000
  Mref_V%Qmu_ref( 331 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 331 ) = 8607.49000000000
  Mref_V%vsh_ref( 331 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 331 ) = 1.00000000000000

  Mref_V%radius_ref( 332 ) = 3148206.00000000
  Mref_V%density_ref( 332 ) = 10400.5100000000
  Mref_V%vpv_ref( 332 ) = 8589.06000000000
  Mref_V%vsv_ref( 332 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 332 ) = 57822.5000000000
  Mref_V%Qmu_ref( 332 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 332 ) = 8589.06000000000
  Mref_V%vsh_ref( 332 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 332 ) = 1.00000000000000

  Mref_V%radius_ref( 333 ) = 3160965.00000000
  Mref_V%density_ref( 333 ) = 10382.6100000000
  Mref_V%vpv_ref( 333 ) = 8570.51000000000
  Mref_V%vsv_ref( 333 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 333 ) = 57822.5000000000
  Mref_V%Qmu_ref( 333 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 333 ) = 8570.51000000000
  Mref_V%vsh_ref( 333 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 333 ) = 1.00000000000000

  Mref_V%radius_ref( 334 ) = 3173725.00000000
  Mref_V%density_ref( 334 ) = 10364.6100000000
  Mref_V%vpv_ref( 334 ) = 8551.83000000000
  Mref_V%vsv_ref( 334 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 334 ) = 57822.5000000000
  Mref_V%Qmu_ref( 334 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 334 ) = 8551.83000000000
  Mref_V%vsh_ref( 334 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 334 ) = 1.00000000000000

  Mref_V%radius_ref( 335 ) = 3186485.00000000
  Mref_V%density_ref( 335 ) = 10346.5200000000
  Mref_V%vpv_ref( 335 ) = 8533.03000000000
  Mref_V%vsv_ref( 335 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 335 ) = 57822.5000000000
  Mref_V%Qmu_ref( 335 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 335 ) = 8533.03000000000
  Mref_V%vsh_ref( 335 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 335 ) = 1.00000000000000

  Mref_V%radius_ref( 336 ) = 3199244.00000000
  Mref_V%density_ref( 336 ) = 10328.3400000000
  Mref_V%vpv_ref( 336 ) = 8514.10000000000
  Mref_V%vsv_ref( 336 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 336 ) = 57822.5000000000
  Mref_V%Qmu_ref( 336 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 336 ) = 8514.10000000000
  Mref_V%vsh_ref( 336 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 336 ) = 1.00000000000000

  Mref_V%radius_ref( 337 ) = 3212004.00000000
  Mref_V%density_ref( 337 ) = 10310.0500000000
  Mref_V%vpv_ref( 337 ) = 8495.05000000000
  Mref_V%vsv_ref( 337 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 337 ) = 57822.5000000000
  Mref_V%Qmu_ref( 337 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 337 ) = 8495.05000000000
  Mref_V%vsh_ref( 337 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 337 ) = 1.00000000000000

  Mref_V%radius_ref( 338 ) = 3224764.00000000
  Mref_V%density_ref( 338 ) = 10291.6800000000
  Mref_V%vpv_ref( 338 ) = 8475.87000000000
  Mref_V%vsv_ref( 338 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 338 ) = 57822.5000000000
  Mref_V%Qmu_ref( 338 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 338 ) = 8475.87000000000
  Mref_V%vsh_ref( 338 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 338 ) = 1.00000000000000

  Mref_V%radius_ref( 339 ) = 3237523.00000000
  Mref_V%density_ref( 339 ) = 10273.2000000000
  Mref_V%vpv_ref( 339 ) = 8456.57000000000
  Mref_V%vsv_ref( 339 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 339 ) = 57822.5000000000
  Mref_V%Qmu_ref( 339 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 339 ) = 8456.57000000000
  Mref_V%vsh_ref( 339 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 339 ) = 1.00000000000000

  Mref_V%radius_ref( 340 ) = 3250283.00000000
  Mref_V%density_ref( 340 ) = 10254.6300000000
  Mref_V%vpv_ref( 340 ) = 8437.14000000000
  Mref_V%vsv_ref( 340 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 340 ) = 57822.5000000000
  Mref_V%Qmu_ref( 340 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 340 ) = 8437.14000000000
  Mref_V%vsh_ref( 340 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 340 ) = 1.00000000000000

  Mref_V%radius_ref( 341 ) = 3263043.00000000
  Mref_V%density_ref( 341 ) = 10235.9600000000
  Mref_V%vpv_ref( 341 ) = 8417.58000000000
  Mref_V%vsv_ref( 341 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 341 ) = 57822.5000000000
  Mref_V%Qmu_ref( 341 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 341 ) = 8417.58000000000
  Mref_V%vsh_ref( 341 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 341 ) = 1.00000000000000

  Mref_V%radius_ref( 342 ) = 3275803.00000000
  Mref_V%density_ref( 342 ) = 10217.2000000000
  Mref_V%vpv_ref( 342 ) = 8397.89000000000
  Mref_V%vsv_ref( 342 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 342 ) = 57822.5000000000
  Mref_V%Qmu_ref( 342 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 342 ) = 8397.89000000000
  Mref_V%vsh_ref( 342 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 342 ) = 1.00000000000000

  Mref_V%radius_ref( 343 ) = 3288562.00000000
  Mref_V%density_ref( 343 ) = 10198.3300000000
  Mref_V%vpv_ref( 343 ) = 8378.07000000000
  Mref_V%vsv_ref( 343 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 343 ) = 57822.5000000000
  Mref_V%Qmu_ref( 343 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 343 ) = 8378.07000000000
  Mref_V%vsh_ref( 343 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 343 ) = 1.00000000000000

  Mref_V%radius_ref( 344 ) = 3301322.00000000
  Mref_V%density_ref( 344 ) = 10179.3700000000
  Mref_V%vpv_ref( 344 ) = 8358.12000000000
  Mref_V%vsv_ref( 344 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 344 ) = 57822.5000000000
  Mref_V%Qmu_ref( 344 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 344 ) = 8358.12000000000
  Mref_V%vsh_ref( 344 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 344 ) = 1.00000000000000

  Mref_V%radius_ref( 345 ) = 3314082.00000000
  Mref_V%density_ref( 345 ) = 10160.3100000000
  Mref_V%vpv_ref( 345 ) = 8338.04000000000
  Mref_V%vsv_ref( 345 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 345 ) = 57822.5000000000
  Mref_V%Qmu_ref( 345 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 345 ) = 8338.04000000000
  Mref_V%vsh_ref( 345 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 345 ) = 1.00000000000000

  Mref_V%radius_ref( 346 ) = 3326841.00000000
  Mref_V%density_ref( 346 ) = 10141.1500000000
  Mref_V%vpv_ref( 346 ) = 8317.83000000000
  Mref_V%vsv_ref( 346 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 346 ) = 57822.5000000000
  Mref_V%Qmu_ref( 346 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 346 ) = 8317.83000000000
  Mref_V%vsh_ref( 346 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 346 ) = 1.00000000000000

  Mref_V%radius_ref( 347 ) = 3339601.00000000
  Mref_V%density_ref( 347 ) = 10121.9000000000
  Mref_V%vpv_ref( 347 ) = 8297.49000000000
  Mref_V%vsv_ref( 347 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 347 ) = 57822.5000000000
  Mref_V%Qmu_ref( 347 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 347 ) = 8297.49000000000
  Mref_V%vsh_ref( 347 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 347 ) = 1.00000000000000

  Mref_V%radius_ref( 348 ) = 3352361.00000000
  Mref_V%density_ref( 348 ) = 10102.5400000000
  Mref_V%vpv_ref( 348 ) = 8277.01000000000
  Mref_V%vsv_ref( 348 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 348 ) = 57822.5000000000
  Mref_V%Qmu_ref( 348 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 348 ) = 8277.01000000000
  Mref_V%vsh_ref( 348 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 348 ) = 1.00000000000000

  Mref_V%radius_ref( 349 ) = 3365120.00000000
  Mref_V%density_ref( 349 ) = 10083.0900000000
  Mref_V%vpv_ref( 349 ) = 8256.41000000000
  Mref_V%vsv_ref( 349 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 349 ) = 57822.5000000000
  Mref_V%Qmu_ref( 349 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 349 ) = 8256.41000000000
  Mref_V%vsh_ref( 349 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 349 ) = 1.00000000000000

  Mref_V%radius_ref( 350 ) = 3377880.00000000
  Mref_V%density_ref( 350 ) = 10063.5300000000
  Mref_V%vpv_ref( 350 ) = 8235.66000000000
  Mref_V%vsv_ref( 350 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 350 ) = 57822.5000000000
  Mref_V%Qmu_ref( 350 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 350 ) = 8235.66000000000
  Mref_V%vsh_ref( 350 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 350 ) = 1.00000000000000

  Mref_V%radius_ref( 351 ) = 3390640.00000000
  Mref_V%density_ref( 351 ) = 10043.8800000000
  Mref_V%vpv_ref( 351 ) = 8214.79000000000
  Mref_V%vsv_ref( 351 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 351 ) = 57822.5000000000
  Mref_V%Qmu_ref( 351 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 351 ) = 8214.79000000000
  Mref_V%vsh_ref( 351 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 351 ) = 1.00000000000000

  Mref_V%radius_ref( 352 ) = 3403400.00000000
  Mref_V%density_ref( 352 ) = 10024.1200000000
  Mref_V%vpv_ref( 352 ) = 8193.77000000000
  Mref_V%vsv_ref( 352 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 352 ) = 57822.5000000000
  Mref_V%Qmu_ref( 352 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 352 ) = 8193.77000000000
  Mref_V%vsh_ref( 352 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 352 ) = 1.00000000000000

  Mref_V%radius_ref( 353 ) = 3416159.00000000
  Mref_V%density_ref( 353 ) = 10004.2700000000
  Mref_V%vpv_ref( 353 ) = 8172.62000000000
  Mref_V%vsv_ref( 353 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 353 ) = 57822.5000000000
  Mref_V%Qmu_ref( 353 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 353 ) = 8172.62000000000
  Mref_V%vsh_ref( 353 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 353 ) = 1.00000000000000

  Mref_V%radius_ref( 354 ) = 3428919.00000000
  Mref_V%density_ref( 354 ) = 9984.32000000000
  Mref_V%vpv_ref( 354 ) = 8151.34000000000
  Mref_V%vsv_ref( 354 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 354 ) = 57822.5000000000
  Mref_V%Qmu_ref( 354 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 354 ) = 8151.34000000000
  Mref_V%vsh_ref( 354 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 354 ) = 1.00000000000000

  Mref_V%radius_ref( 355 ) = 3441679.00000000
  Mref_V%density_ref( 355 ) = 9964.26000000000
  Mref_V%vpv_ref( 355 ) = 8129.92000000000
  Mref_V%vsv_ref( 355 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 355 ) = 57822.5000000000
  Mref_V%Qmu_ref( 355 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 355 ) = 8129.92000000000
  Mref_V%vsh_ref( 355 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 355 ) = 1.00000000000000

  Mref_V%radius_ref( 356 ) = 3454438.00000000
  Mref_V%density_ref( 356 ) = 9944.10000000000
  Mref_V%vpv_ref( 356 ) = 8108.36000000000
  Mref_V%vsv_ref( 356 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 356 ) = 57822.5000000000
  Mref_V%Qmu_ref( 356 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 356 ) = 8108.36000000000
  Mref_V%vsh_ref( 356 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 356 ) = 1.00000000000000

  Mref_V%radius_ref( 357 ) = 3467198.00000000
  Mref_V%density_ref( 357 ) = 9923.84000000000
  Mref_V%vpv_ref( 357 ) = 8086.66000000000
  Mref_V%vsv_ref( 357 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 357 ) = 57822.5000000000
  Mref_V%Qmu_ref( 357 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 357 ) = 8086.66000000000
  Mref_V%vsh_ref( 357 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 357 ) = 1.00000000000000

  Mref_V%radius_ref( 358 ) = 3479958.00000000
  Mref_V%density_ref( 358 ) = 9903.48000000000
  Mref_V%vpv_ref( 358 ) = 8064.82000000000
  Mref_V%vsv_ref( 358 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 358 ) = 57822.5000000000
  Mref_V%Qmu_ref( 358 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 358 ) = 8064.82000000000
  Mref_V%vsh_ref( 358 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 358 ) = 1.00000000000000

  Mref_V%radius_ref( 359 ) = 3479958.00000000
  Mref_V%density_ref( 359 ) = 5566.45000000000
  Mref_V%vpv_ref( 359 ) = 13716.6000000000
  Mref_V%vsv_ref( 359 ) = 7264.66000000000
  Mref_V%Qkappa_ref( 359 ) = 57822.5000000000
  Mref_V%Qmu_ref( 359 ) = 355.000000000000
  Mref_V%vph_ref( 359 ) = 13716.6000000000
  Mref_V%vsh_ref( 359 ) = 7264.66000000000
  Mref_V%eta_ref( 359 ) = 1.00000000000000

  Mref_V%radius_ref( 360 ) = 3489335.00000000
  Mref_V%density_ref( 360 ) = 5561.75000000000
  Mref_V%vpv_ref( 360 ) = 13714.2900000000
  Mref_V%vsv_ref( 360 ) = 7264.75000000000
  Mref_V%Qkappa_ref( 360 ) = 57822.5000000000
  Mref_V%Qmu_ref( 360 ) = 355.000000000000
  Mref_V%vph_ref( 360 ) = 13714.2900000000
  Mref_V%vsh_ref( 360 ) = 7264.75000000000
  Mref_V%eta_ref( 360 ) = 1.00000000000000

  Mref_V%radius_ref( 361 ) = 3498713.00000000
  Mref_V%density_ref( 361 ) = 5557.05000000000
  Mref_V%vpv_ref( 361 ) = 13712.0000000000
  Mref_V%vsv_ref( 361 ) = 7264.85000000000
  Mref_V%Qkappa_ref( 361 ) = 57822.5000000000
  Mref_V%Qmu_ref( 361 ) = 355.000000000000
  Mref_V%vph_ref( 361 ) = 13712.0000000000
  Mref_V%vsh_ref( 361 ) = 7264.85000000000
  Mref_V%eta_ref( 361 ) = 1.00000000000000

  Mref_V%radius_ref( 362 ) = 3508091.00000000
  Mref_V%density_ref( 362 ) = 5552.36000000000
  Mref_V%vpv_ref( 362 ) = 13709.7000000000
  Mref_V%vsv_ref( 362 ) = 7264.94000000000
  Mref_V%Qkappa_ref( 362 ) = 57822.5000000000
  Mref_V%Qmu_ref( 362 ) = 355.000000000000
  Mref_V%vph_ref( 362 ) = 13709.7000000000
  Mref_V%vsh_ref( 362 ) = 7264.94000000000
  Mref_V%eta_ref( 362 ) = 1.00000000000000

  Mref_V%radius_ref( 363 ) = 3517468.00000000
  Mref_V%density_ref( 363 ) = 5547.66000000000
  Mref_V%vpv_ref( 363 ) = 13707.4200000000
  Mref_V%vsv_ref( 363 ) = 7265.03000000000
  Mref_V%Qkappa_ref( 363 ) = 57822.5000000000
  Mref_V%Qmu_ref( 363 ) = 355.000000000000
  Mref_V%vph_ref( 363 ) = 13707.4200000000
  Mref_V%vsh_ref( 363 ) = 7265.03000000000
  Mref_V%eta_ref( 363 ) = 1.00000000000000

  Mref_V%radius_ref( 364 ) = 3526846.00000000
  Mref_V%density_ref( 364 ) = 5542.97000000000
  Mref_V%vpv_ref( 364 ) = 13705.1400000000
  Mref_V%vsv_ref( 364 ) = 7265.12000000000
  Mref_V%Qkappa_ref( 364 ) = 57822.5000000000
  Mref_V%Qmu_ref( 364 ) = 355.000000000000
  Mref_V%vph_ref( 364 ) = 13705.1400000000
  Mref_V%vsh_ref( 364 ) = 7265.12000000000
  Mref_V%eta_ref( 364 ) = 1.00000000000000

  Mref_V%radius_ref( 365 ) = 3536224.00000000
  Mref_V%density_ref( 365 ) = 5538.28000000000
  Mref_V%vpv_ref( 365 ) = 13702.8600000000
  Mref_V%vsv_ref( 365 ) = 7265.21000000000
  Mref_V%Qkappa_ref( 365 ) = 57822.5000000000
  Mref_V%Qmu_ref( 365 ) = 355.000000000000
  Mref_V%vph_ref( 365 ) = 13702.8600000000
  Mref_V%vsh_ref( 365 ) = 7265.21000000000
  Mref_V%eta_ref( 365 ) = 1.00000000000000

  Mref_V%radius_ref( 366 ) = 3545601.00000000
  Mref_V%density_ref( 366 ) = 5533.59000000000
  Mref_V%vpv_ref( 366 ) = 13700.5900000000
  Mref_V%vsv_ref( 366 ) = 7265.29000000000
  Mref_V%Qkappa_ref( 366 ) = 57822.5000000000
  Mref_V%Qmu_ref( 366 ) = 355.000000000000
  Mref_V%vph_ref( 366 ) = 13700.5900000000
  Mref_V%vsh_ref( 366 ) = 7265.29000000000
  Mref_V%eta_ref( 366 ) = 1.00000000000000

  Mref_V%radius_ref( 367 ) = 3554979.00000000
  Mref_V%density_ref( 367 ) = 5528.90000000000
  Mref_V%vpv_ref( 367 ) = 13698.3300000000
  Mref_V%vsv_ref( 367 ) = 7265.38000000000
  Mref_V%Qkappa_ref( 367 ) = 57822.5000000000
  Mref_V%Qmu_ref( 367 ) = 355.000000000000
  Mref_V%vph_ref( 367 ) = 13698.3300000000
  Mref_V%vsh_ref( 367 ) = 7265.38000000000
  Mref_V%eta_ref( 367 ) = 1.00000000000000

  Mref_V%radius_ref( 368 ) = 3564357.00000000
  Mref_V%density_ref( 368 ) = 5524.21000000000
  Mref_V%vpv_ref( 368 ) = 13696.0700000000
  Mref_V%vsv_ref( 368 ) = 7265.46000000000
  Mref_V%Qkappa_ref( 368 ) = 57822.5000000000
  Mref_V%Qmu_ref( 368 ) = 355.000000000000
  Mref_V%vph_ref( 368 ) = 13696.0700000000
  Mref_V%vsh_ref( 368 ) = 7265.46000000000
  Mref_V%eta_ref( 368 ) = 1.00000000000000

  Mref_V%radius_ref( 369 ) = 3573734.00000000
  Mref_V%density_ref( 369 ) = 5519.53000000000
  Mref_V%vpv_ref( 369 ) = 13693.8200000000
  Mref_V%vsv_ref( 369 ) = 7265.54000000000
  Mref_V%Qkappa_ref( 369 ) = 57822.5000000000
  Mref_V%Qmu_ref( 369 ) = 355.000000000000
  Mref_V%vph_ref( 369 ) = 13693.8200000000
  Mref_V%vsh_ref( 369 ) = 7265.54000000000
  Mref_V%eta_ref( 369 ) = 1.00000000000000

  Mref_V%radius_ref( 370 ) = 3583112.00000000
  Mref_V%density_ref( 370 ) = 5514.85000000000
  Mref_V%vpv_ref( 370 ) = 13691.5700000000
  Mref_V%vsv_ref( 370 ) = 7265.62000000000
  Mref_V%Qkappa_ref( 370 ) = 57822.5000000000
  Mref_V%Qmu_ref( 370 ) = 355.000000000000
  Mref_V%vph_ref( 370 ) = 13691.5700000000
  Mref_V%vsh_ref( 370 ) = 7265.62000000000
  Mref_V%eta_ref( 370 ) = 1.00000000000000

  Mref_V%radius_ref( 371 ) = 3592489.00000000
  Mref_V%density_ref( 371 ) = 5510.16000000000
  Mref_V%vpv_ref( 371 ) = 13689.3300000000
  Mref_V%vsv_ref( 371 ) = 7265.69000000000
  Mref_V%Qkappa_ref( 371 ) = 57822.5000000000
  Mref_V%Qmu_ref( 371 ) = 355.000000000000
  Mref_V%vph_ref( 371 ) = 13689.3300000000
  Mref_V%vsh_ref( 371 ) = 7265.69000000000
  Mref_V%eta_ref( 371 ) = 1.00000000000000

  Mref_V%radius_ref( 372 ) = 3601867.00000000
  Mref_V%density_ref( 372 ) = 5505.48000000000
  Mref_V%vpv_ref( 372 ) = 13687.0900000000
  Mref_V%vsv_ref( 372 ) = 7265.76000000000
  Mref_V%Qkappa_ref( 372 ) = 57822.5000000000
  Mref_V%Qmu_ref( 372 ) = 355.000000000000
  Mref_V%vph_ref( 372 ) = 13687.0900000000
  Mref_V%vsh_ref( 372 ) = 7265.76000000000
  Mref_V%eta_ref( 372 ) = 1.00000000000000

  Mref_V%radius_ref( 373 ) = 3611245.00000000
  Mref_V%density_ref( 373 ) = 5500.81000000000
  Mref_V%vpv_ref( 373 ) = 13684.8600000000
  Mref_V%vsv_ref( 373 ) = 7265.84000000000
  Mref_V%Qkappa_ref( 373 ) = 57822.5000000000
  Mref_V%Qmu_ref( 373 ) = 355.000000000000
  Mref_V%vph_ref( 373 ) = 13684.8600000000
  Mref_V%vsh_ref( 373 ) = 7265.84000000000
  Mref_V%eta_ref( 373 ) = 1.00000000000000

  Mref_V%radius_ref( 374 ) = 3620622.00000000
  Mref_V%density_ref( 374 ) = 5496.13000000000
  Mref_V%vpv_ref( 374 ) = 13682.6300000000
  Mref_V%vsv_ref( 374 ) = 7265.91000000000
  Mref_V%Qkappa_ref( 374 ) = 57822.5000000000
  Mref_V%Qmu_ref( 374 ) = 355.000000000000
  Mref_V%vph_ref( 374 ) = 13682.6300000000
  Mref_V%vsh_ref( 374 ) = 7265.91000000000
  Mref_V%eta_ref( 374 ) = 1.00000000000000

  Mref_V%radius_ref( 375 ) = 3630000.00000000
  Mref_V%density_ref( 375 ) = 5491.45000000000
  Mref_V%vpv_ref( 375 ) = 13680.4100000000
  Mref_V%vsv_ref( 375 ) = 7265.97000000000
  Mref_V%Qkappa_ref( 375 ) = 57822.5000000000
  Mref_V%Qmu_ref( 375 ) = 355.000000000000
  Mref_V%vph_ref( 375 ) = 13680.4100000000
  Mref_V%vsh_ref( 375 ) = 7265.97000000000
  Mref_V%eta_ref( 375 ) = 1.00000000000000

  Mref_V%radius_ref( 376 ) = 3630000.00000000
  Mref_V%density_ref( 376 ) = 5491.45000000000
  Mref_V%vpv_ref( 376 ) = 13680.4100000000
  Mref_V%vsv_ref( 376 ) = 7265.97000000000
  Mref_V%Qkappa_ref( 376 ) = 57822.5000000000
  Mref_V%Qmu_ref( 376 ) = 355.000000000000
  Mref_V%vph_ref( 376 ) = 13680.4100000000
  Mref_V%vsh_ref( 376 ) = 7265.97000000000
  Mref_V%eta_ref( 376 ) = 1.00000000000000

  Mref_V%radius_ref( 377 ) = 3639471.00000000
  Mref_V%density_ref( 377 ) = 5486.73000000000
  Mref_V%vpv_ref( 377 ) = 13668.9000000000
  Mref_V%vsv_ref( 377 ) = 7261.63000000000
  Mref_V%Qkappa_ref( 377 ) = 57822.5000000000
  Mref_V%Qmu_ref( 377 ) = 355.000000000000
  Mref_V%vph_ref( 377 ) = 13668.9000000000
  Mref_V%vsh_ref( 377 ) = 7261.63000000000
  Mref_V%eta_ref( 377 ) = 1.00000000000000

  Mref_V%radius_ref( 378 ) = 3648942.00000000
  Mref_V%density_ref( 378 ) = 5482.01000000000
  Mref_V%vpv_ref( 378 ) = 13657.4300000000
  Mref_V%vsv_ref( 378 ) = 7257.29000000000
  Mref_V%Qkappa_ref( 378 ) = 57822.5000000000
  Mref_V%Qmu_ref( 378 ) = 355.000000000000
  Mref_V%vph_ref( 378 ) = 13657.4300000000
  Mref_V%vsh_ref( 378 ) = 7257.29000000000
  Mref_V%eta_ref( 378 ) = 1.00000000000000

  Mref_V%radius_ref( 379 ) = 3658413.00000000
  Mref_V%density_ref( 379 ) = 5477.29000000000
  Mref_V%vpv_ref( 379 ) = 13645.9700000000
  Mref_V%vsv_ref( 379 ) = 7252.97000000000
  Mref_V%Qkappa_ref( 379 ) = 57822.5000000000
  Mref_V%Qmu_ref( 379 ) = 355.000000000000
  Mref_V%vph_ref( 379 ) = 13645.9700000000
  Mref_V%vsh_ref( 379 ) = 7252.97000000000
  Mref_V%eta_ref( 379 ) = 1.00000000000000

  Mref_V%radius_ref( 380 ) = 3667885.00000000
  Mref_V%density_ref( 380 ) = 5472.57000000000
  Mref_V%vpv_ref( 380 ) = 13634.5400000000
  Mref_V%vsv_ref( 380 ) = 7248.64000000000
  Mref_V%Qkappa_ref( 380 ) = 57822.5000000000
  Mref_V%Qmu_ref( 380 ) = 355.000000000000
  Mref_V%vph_ref( 380 ) = 13634.5400000000
  Mref_V%vsh_ref( 380 ) = 7248.64000000000
  Mref_V%eta_ref( 380 ) = 1.00000000000000

  Mref_V%radius_ref( 381 ) = 3677356.00000000
  Mref_V%density_ref( 381 ) = 5467.85000000000
  Mref_V%vpv_ref( 381 ) = 13623.1400000000
  Mref_V%vsv_ref( 381 ) = 7244.33000000000
  Mref_V%Qkappa_ref( 381 ) = 57822.5000000000
  Mref_V%Qmu_ref( 381 ) = 355.000000000000
  Mref_V%vph_ref( 381 ) = 13623.1400000000
  Mref_V%vsh_ref( 381 ) = 7244.33000000000
  Mref_V%eta_ref( 381 ) = 1.00000000000000

  Mref_V%radius_ref( 382 ) = 3686827.00000000
  Mref_V%density_ref( 382 ) = 5463.13000000000
  Mref_V%vpv_ref( 382 ) = 13611.7600000000
  Mref_V%vsv_ref( 382 ) = 7240.01000000000
  Mref_V%Qkappa_ref( 382 ) = 57822.5000000000
  Mref_V%Qmu_ref( 382 ) = 355.000000000000
  Mref_V%vph_ref( 382 ) = 13611.7600000000
  Mref_V%vsh_ref( 382 ) = 7240.01000000000
  Mref_V%eta_ref( 382 ) = 1.00000000000000

  Mref_V%radius_ref( 383 ) = 3696298.00000000
  Mref_V%density_ref( 383 ) = 5458.42000000000
  Mref_V%vpv_ref( 383 ) = 13600.4000000000
  Mref_V%vsv_ref( 383 ) = 7235.71000000000
  Mref_V%Qkappa_ref( 383 ) = 57822.5000000000
  Mref_V%Qmu_ref( 383 ) = 355.000000000000
  Mref_V%vph_ref( 383 ) = 13600.4000000000
  Mref_V%vsh_ref( 383 ) = 7235.71000000000
  Mref_V%eta_ref( 383 ) = 1.00000000000000

  Mref_V%radius_ref( 384 ) = 3705769.00000000
  Mref_V%density_ref( 384 ) = 5453.70000000000
  Mref_V%vpv_ref( 384 ) = 13589.0700000000
  Mref_V%vsv_ref( 384 ) = 7231.41000000000
  Mref_V%Qkappa_ref( 384 ) = 57822.5000000000
  Mref_V%Qmu_ref( 384 ) = 355.000000000000
  Mref_V%vph_ref( 384 ) = 13589.0700000000
  Mref_V%vsh_ref( 384 ) = 7231.41000000000
  Mref_V%eta_ref( 384 ) = 1.00000000000000

  Mref_V%radius_ref( 385 ) = 3715240.00000000
  Mref_V%density_ref( 385 ) = 5448.99000000000
  Mref_V%vpv_ref( 385 ) = 13577.7600000000
  Mref_V%vsv_ref( 385 ) = 7227.12000000000
  Mref_V%Qkappa_ref( 385 ) = 57822.5000000000
  Mref_V%Qmu_ref( 385 ) = 355.000000000000
  Mref_V%vph_ref( 385 ) = 13577.7600000000
  Mref_V%vsh_ref( 385 ) = 7227.12000000000
  Mref_V%eta_ref( 385 ) = 1.00000000000000

  Mref_V%radius_ref( 386 ) = 3724712.00000000
  Mref_V%density_ref( 386 ) = 5444.27000000000
  Mref_V%vpv_ref( 386 ) = 13566.4700000000
  Mref_V%vsv_ref( 386 ) = 7222.83000000000
  Mref_V%Qkappa_ref( 386 ) = 57822.5000000000
  Mref_V%Qmu_ref( 386 ) = 355.000000000000
  Mref_V%vph_ref( 386 ) = 13566.4700000000
  Mref_V%vsh_ref( 386 ) = 7222.83000000000
  Mref_V%eta_ref( 386 ) = 1.00000000000000

  Mref_V%radius_ref( 387 ) = 3734183.00000000
  Mref_V%density_ref( 387 ) = 5439.56000000000
  Mref_V%vpv_ref( 387 ) = 13555.2000000000
  Mref_V%vsv_ref( 387 ) = 7218.55000000000
  Mref_V%Qkappa_ref( 387 ) = 57822.5000000000
  Mref_V%Qmu_ref( 387 ) = 355.000000000000
  Mref_V%vph_ref( 387 ) = 13555.2000000000
  Mref_V%vsh_ref( 387 ) = 7218.55000000000
  Mref_V%eta_ref( 387 ) = 1.00000000000000

  Mref_V%radius_ref( 388 ) = 3743654.00000000
  Mref_V%density_ref( 388 ) = 5434.85000000000
  Mref_V%vpv_ref( 388 ) = 13543.9500000000
  Mref_V%vsv_ref( 388 ) = 7214.27000000000
  Mref_V%Qkappa_ref( 388 ) = 57822.5000000000
  Mref_V%Qmu_ref( 388 ) = 355.000000000000
  Mref_V%vph_ref( 388 ) = 13543.9500000000
  Mref_V%vsh_ref( 388 ) = 7214.27000000000
  Mref_V%eta_ref( 388 ) = 1.00000000000000

  Mref_V%radius_ref( 389 ) = 3753125.00000000
  Mref_V%density_ref( 389 ) = 5430.13000000000
  Mref_V%vpv_ref( 389 ) = 13532.7200000000
  Mref_V%vsv_ref( 389 ) = 7210.00000000000
  Mref_V%Qkappa_ref( 389 ) = 57822.5000000000
  Mref_V%Qmu_ref( 389 ) = 355.000000000000
  Mref_V%vph_ref( 389 ) = 13532.7200000000
  Mref_V%vsh_ref( 389 ) = 7210.00000000000
  Mref_V%eta_ref( 389 ) = 1.00000000000000

  Mref_V%radius_ref( 390 ) = 3762596.00000000
  Mref_V%density_ref( 390 ) = 5425.42000000000
  Mref_V%vpv_ref( 390 ) = 13521.5100000000
  Mref_V%vsv_ref( 390 ) = 7205.73000000000
  Mref_V%Qkappa_ref( 390 ) = 57822.5000000000
  Mref_V%Qmu_ref( 390 ) = 355.000000000000
  Mref_V%vph_ref( 390 ) = 13521.5100000000
  Mref_V%vsh_ref( 390 ) = 7205.73000000000
  Mref_V%eta_ref( 390 ) = 1.00000000000000

  Mref_V%radius_ref( 391 ) = 3772067.00000000
  Mref_V%density_ref( 391 ) = 5420.71000000000
  Mref_V%vpv_ref( 391 ) = 13510.3200000000
  Mref_V%vsv_ref( 391 ) = 7201.47000000000
  Mref_V%Qkappa_ref( 391 ) = 57822.5000000000
  Mref_V%Qmu_ref( 391 ) = 355.000000000000
  Mref_V%vph_ref( 391 ) = 13510.3200000000
  Mref_V%vsh_ref( 391 ) = 7201.47000000000
  Mref_V%eta_ref( 391 ) = 1.00000000000000

  Mref_V%radius_ref( 392 ) = 3781538.00000000
  Mref_V%density_ref( 392 ) = 5416.00000000000
  Mref_V%vpv_ref( 392 ) = 13499.1400000000
  Mref_V%vsv_ref( 392 ) = 7197.21000000000
  Mref_V%Qkappa_ref( 392 ) = 57822.5000000000
  Mref_V%Qmu_ref( 392 ) = 355.000000000000
  Mref_V%vph_ref( 392 ) = 13499.1400000000
  Mref_V%vsh_ref( 392 ) = 7197.21000000000
  Mref_V%eta_ref( 392 ) = 1.00000000000000

  Mref_V%radius_ref( 393 ) = 3791010.00000000
  Mref_V%density_ref( 393 ) = 5411.29000000000
  Mref_V%vpv_ref( 393 ) = 13487.9900000000
  Mref_V%vsv_ref( 393 ) = 7192.95000000000
  Mref_V%Qkappa_ref( 393 ) = 57822.5000000000
  Mref_V%Qmu_ref( 393 ) = 355.000000000000
  Mref_V%vph_ref( 393 ) = 13487.9900000000
  Mref_V%vsh_ref( 393 ) = 7192.95000000000
  Mref_V%eta_ref( 393 ) = 1.00000000000000

  Mref_V%radius_ref( 394 ) = 3800481.00000000
  Mref_V%density_ref( 394 ) = 5406.57000000000
  Mref_V%vpv_ref( 394 ) = 13476.8500000000
  Mref_V%vsv_ref( 394 ) = 7188.70000000000
  Mref_V%Qkappa_ref( 394 ) = 57822.5000000000
  Mref_V%Qmu_ref( 394 ) = 355.000000000000
  Mref_V%vph_ref( 394 ) = 13476.8500000000
  Mref_V%vsh_ref( 394 ) = 7188.70000000000
  Mref_V%eta_ref( 394 ) = 1.00000000000000

  Mref_V%radius_ref( 395 ) = 3809952.00000000
  Mref_V%density_ref( 395 ) = 5401.86000000000
  Mref_V%vpv_ref( 395 ) = 13465.7300000000
  Mref_V%vsv_ref( 395 ) = 7184.45000000000
  Mref_V%Qkappa_ref( 395 ) = 57822.5000000000
  Mref_V%Qmu_ref( 395 ) = 355.000000000000
  Mref_V%vph_ref( 395 ) = 13465.7300000000
  Mref_V%vsh_ref( 395 ) = 7184.45000000000
  Mref_V%eta_ref( 395 ) = 1.00000000000000

  Mref_V%radius_ref( 396 ) = 3819423.00000000
  Mref_V%density_ref( 396 ) = 5397.15000000000
  Mref_V%vpv_ref( 396 ) = 13454.6300000000
  Mref_V%vsv_ref( 396 ) = 7180.21000000000
  Mref_V%Qkappa_ref( 396 ) = 57822.5000000000
  Mref_V%Qmu_ref( 396 ) = 355.000000000000
  Mref_V%vph_ref( 396 ) = 13454.6300000000
  Mref_V%vsh_ref( 396 ) = 7180.21000000000
  Mref_V%eta_ref( 396 ) = 1.00000000000000

  Mref_V%radius_ref( 397 ) = 3828894.00000000
  Mref_V%density_ref( 397 ) = 5392.44000000000
  Mref_V%vpv_ref( 397 ) = 13443.5400000000
  Mref_V%vsv_ref( 397 ) = 7175.97000000000
  Mref_V%Qkappa_ref( 397 ) = 57822.5000000000
  Mref_V%Qmu_ref( 397 ) = 355.000000000000
  Mref_V%vph_ref( 397 ) = 13443.5400000000
  Mref_V%vsh_ref( 397 ) = 7175.97000000000
  Mref_V%eta_ref( 397 ) = 1.00000000000000

  Mref_V%radius_ref( 398 ) = 3838365.00000000
  Mref_V%density_ref( 398 ) = 5387.73000000000
  Mref_V%vpv_ref( 398 ) = 13432.4600000000
  Mref_V%vsv_ref( 398 ) = 7171.73000000000
  Mref_V%Qkappa_ref( 398 ) = 57822.5000000000
  Mref_V%Qmu_ref( 398 ) = 355.000000000000
  Mref_V%vph_ref( 398 ) = 13432.4600000000
  Mref_V%vsh_ref( 398 ) = 7171.73000000000
  Mref_V%eta_ref( 398 ) = 1.00000000000000

  Mref_V%radius_ref( 399 ) = 3847837.00000000
  Mref_V%density_ref( 399 ) = 5383.02000000000
  Mref_V%vpv_ref( 399 ) = 13421.4100000000
  Mref_V%vsv_ref( 399 ) = 7167.50000000000
  Mref_V%Qkappa_ref( 399 ) = 57822.5000000000
  Mref_V%Qmu_ref( 399 ) = 355.000000000000
  Mref_V%vph_ref( 399 ) = 13421.4100000000
  Mref_V%vsh_ref( 399 ) = 7167.50000000000
  Mref_V%eta_ref( 399 ) = 1.00000000000000

  Mref_V%radius_ref( 400 ) = 3857308.00000000
  Mref_V%density_ref( 400 ) = 5378.30000000000
  Mref_V%vpv_ref( 400 ) = 13410.3600000000
  Mref_V%vsv_ref( 400 ) = 7163.27000000000
  Mref_V%Qkappa_ref( 400 ) = 57822.5000000000
  Mref_V%Qmu_ref( 400 ) = 355.000000000000
  Mref_V%vph_ref( 400 ) = 13410.3600000000
  Mref_V%vsh_ref( 400 ) = 7163.27000000000
  Mref_V%eta_ref( 400 ) = 1.00000000000000

  Mref_V%radius_ref( 401 ) = 3866779.00000000
  Mref_V%density_ref( 401 ) = 5373.59000000000
  Mref_V%vpv_ref( 401 ) = 13399.3300000000
  Mref_V%vsv_ref( 401 ) = 7159.04000000000
  Mref_V%Qkappa_ref( 401 ) = 57822.5000000000
  Mref_V%Qmu_ref( 401 ) = 355.000000000000
  Mref_V%vph_ref( 401 ) = 13399.3300000000
  Mref_V%vsh_ref( 401 ) = 7159.04000000000
  Mref_V%eta_ref( 401 ) = 1.00000000000000

  Mref_V%radius_ref( 402 ) = 3876250.00000000
  Mref_V%density_ref( 402 ) = 5368.88000000000
  Mref_V%vpv_ref( 402 ) = 13388.3100000000
  Mref_V%vsv_ref( 402 ) = 7154.81000000000
  Mref_V%Qkappa_ref( 402 ) = 57822.5000000000
  Mref_V%Qmu_ref( 402 ) = 355.000000000000
  Mref_V%vph_ref( 402 ) = 13388.3100000000
  Mref_V%vsh_ref( 402 ) = 7154.81000000000
  Mref_V%eta_ref( 402 ) = 1.00000000000000

  Mref_V%radius_ref( 403 ) = 3885721.00000000
  Mref_V%density_ref( 403 ) = 5364.17000000000
  Mref_V%vpv_ref( 403 ) = 13377.3100000000
  Mref_V%vsv_ref( 403 ) = 7150.59000000000
  Mref_V%Qkappa_ref( 403 ) = 57822.5000000000
  Mref_V%Qmu_ref( 403 ) = 355.000000000000
  Mref_V%vph_ref( 403 ) = 13377.3100000000
  Mref_V%vsh_ref( 403 ) = 7150.59000000000
  Mref_V%eta_ref( 403 ) = 1.00000000000000

  Mref_V%radius_ref( 404 ) = 3895192.00000000
  Mref_V%density_ref( 404 ) = 5359.45000000000
  Mref_V%vpv_ref( 404 ) = 13366.3100000000
  Mref_V%vsv_ref( 404 ) = 7146.37000000000
  Mref_V%Qkappa_ref( 404 ) = 57822.5000000000
  Mref_V%Qmu_ref( 404 ) = 355.000000000000
  Mref_V%vph_ref( 404 ) = 13366.3100000000
  Mref_V%vsh_ref( 404 ) = 7146.37000000000
  Mref_V%eta_ref( 404 ) = 1.00000000000000

  Mref_V%radius_ref( 405 ) = 3904663.00000000
  Mref_V%density_ref( 405 ) = 5354.74000000000
  Mref_V%vpv_ref( 405 ) = 13355.3300000000
  Mref_V%vsv_ref( 405 ) = 7142.15000000000
  Mref_V%Qkappa_ref( 405 ) = 57822.5000000000
  Mref_V%Qmu_ref( 405 ) = 355.000000000000
  Mref_V%vph_ref( 405 ) = 13355.3300000000
  Mref_V%vsh_ref( 405 ) = 7142.15000000000
  Mref_V%eta_ref( 405 ) = 1.00000000000000

  Mref_V%radius_ref( 406 ) = 3914135.00000000
  Mref_V%density_ref( 406 ) = 5350.02000000000
  Mref_V%vpv_ref( 406 ) = 13344.3600000000
  Mref_V%vsv_ref( 406 ) = 7137.93000000000
  Mref_V%Qkappa_ref( 406 ) = 57822.5000000000
  Mref_V%Qmu_ref( 406 ) = 355.000000000000
  Mref_V%vph_ref( 406 ) = 13344.3600000000
  Mref_V%vsh_ref( 406 ) = 7137.93000000000
  Mref_V%eta_ref( 406 ) = 1.00000000000000

  Mref_V%radius_ref( 407 ) = 3923606.00000000
  Mref_V%density_ref( 407 ) = 5345.31000000000
  Mref_V%vpv_ref( 407 ) = 13333.4000000000
  Mref_V%vsv_ref( 407 ) = 7133.71000000000
  Mref_V%Qkappa_ref( 407 ) = 57822.5000000000
  Mref_V%Qmu_ref( 407 ) = 355.000000000000
  Mref_V%vph_ref( 407 ) = 13333.4000000000
  Mref_V%vsh_ref( 407 ) = 7133.71000000000
  Mref_V%eta_ref( 407 ) = 1.00000000000000

  Mref_V%radius_ref( 408 ) = 3933077.00000000
  Mref_V%density_ref( 408 ) = 5340.59000000000
  Mref_V%vpv_ref( 408 ) = 13322.4500000000
  Mref_V%vsv_ref( 408 ) = 7129.50000000000
  Mref_V%Qkappa_ref( 408 ) = 57822.5000000000
  Mref_V%Qmu_ref( 408 ) = 355.000000000000
  Mref_V%vph_ref( 408 ) = 13322.4500000000
  Mref_V%vsh_ref( 408 ) = 7129.50000000000
  Mref_V%eta_ref( 408 ) = 1.00000000000000

  Mref_V%radius_ref( 409 ) = 3942548.00000000
  Mref_V%density_ref( 409 ) = 5335.87000000000
  Mref_V%vpv_ref( 409 ) = 13311.5100000000
  Mref_V%vsv_ref( 409 ) = 7125.29000000000
  Mref_V%Qkappa_ref( 409 ) = 57822.5000000000
  Mref_V%Qmu_ref( 409 ) = 355.000000000000
  Mref_V%vph_ref( 409 ) = 13311.5100000000
  Mref_V%vsh_ref( 409 ) = 7125.29000000000
  Mref_V%eta_ref( 409 ) = 1.00000000000000

  Mref_V%radius_ref( 410 ) = 3952019.00000000
  Mref_V%density_ref( 410 ) = 5331.16000000000
  Mref_V%vpv_ref( 410 ) = 13300.5800000000
  Mref_V%vsv_ref( 410 ) = 7121.07000000000
  Mref_V%Qkappa_ref( 410 ) = 57822.5000000000
  Mref_V%Qmu_ref( 410 ) = 355.000000000000
  Mref_V%vph_ref( 410 ) = 13300.5800000000
  Mref_V%vsh_ref( 410 ) = 7121.07000000000
  Mref_V%eta_ref( 410 ) = 1.00000000000000

  Mref_V%radius_ref( 411 ) = 3961490.00000000
  Mref_V%density_ref( 411 ) = 5326.44000000000
  Mref_V%vpv_ref( 411 ) = 13289.6600000000
  Mref_V%vsv_ref( 411 ) = 7116.86000000000
  Mref_V%Qkappa_ref( 411 ) = 57822.5000000000
  Mref_V%Qmu_ref( 411 ) = 355.000000000000
  Mref_V%vph_ref( 411 ) = 13289.6600000000
  Mref_V%vsh_ref( 411 ) = 7116.86000000000
  Mref_V%eta_ref( 411 ) = 1.00000000000000

  Mref_V%radius_ref( 412 ) = 3970962.00000000
  Mref_V%density_ref( 412 ) = 5321.72000000000
  Mref_V%vpv_ref( 412 ) = 13278.7400000000
  Mref_V%vsv_ref( 412 ) = 7112.65000000000
  Mref_V%Qkappa_ref( 412 ) = 57822.5000000000
  Mref_V%Qmu_ref( 412 ) = 355.000000000000
  Mref_V%vph_ref( 412 ) = 13278.7400000000
  Mref_V%vsh_ref( 412 ) = 7112.65000000000
  Mref_V%eta_ref( 412 ) = 1.00000000000000

  Mref_V%radius_ref( 413 ) = 3980433.00000000
  Mref_V%density_ref( 413 ) = 5317.00000000000
  Mref_V%vpv_ref( 413 ) = 13267.8400000000
  Mref_V%vsv_ref( 413 ) = 7108.44000000000
  Mref_V%Qkappa_ref( 413 ) = 57822.5000000000
  Mref_V%Qmu_ref( 413 ) = 355.000000000000
  Mref_V%vph_ref( 413 ) = 13267.8400000000
  Mref_V%vsh_ref( 413 ) = 7108.44000000000
  Mref_V%eta_ref( 413 ) = 1.00000000000000

  Mref_V%radius_ref( 414 ) = 3989904.00000000
  Mref_V%density_ref( 414 ) = 5312.28000000000
  Mref_V%vpv_ref( 414 ) = 13256.9300000000
  Mref_V%vsv_ref( 414 ) = 7104.23000000000
  Mref_V%Qkappa_ref( 414 ) = 57822.5000000000
  Mref_V%Qmu_ref( 414 ) = 355.000000000000
  Mref_V%vph_ref( 414 ) = 13256.9300000000
  Mref_V%vsh_ref( 414 ) = 7104.23000000000
  Mref_V%eta_ref( 414 ) = 1.00000000000000

  Mref_V%radius_ref( 415 ) = 3999375.00000000
  Mref_V%density_ref( 415 ) = 5307.55000000000
  Mref_V%vpv_ref( 415 ) = 13246.0400000000
  Mref_V%vsv_ref( 415 ) = 7100.02000000000
  Mref_V%Qkappa_ref( 415 ) = 57822.5000000000
  Mref_V%Qmu_ref( 415 ) = 355.000000000000
  Mref_V%vph_ref( 415 ) = 13246.0400000000
  Mref_V%vsh_ref( 415 ) = 7100.02000000000
  Mref_V%eta_ref( 415 ) = 1.00000000000000

  Mref_V%radius_ref( 416 ) = 4008846.00000000
  Mref_V%density_ref( 416 ) = 5302.83000000000
  Mref_V%vpv_ref( 416 ) = 13235.1500000000
  Mref_V%vsv_ref( 416 ) = 7095.81000000000
  Mref_V%Qkappa_ref( 416 ) = 57822.5000000000
  Mref_V%Qmu_ref( 416 ) = 355.000000000000
  Mref_V%vph_ref( 416 ) = 13235.1500000000
  Mref_V%vsh_ref( 416 ) = 7095.81000000000
  Mref_V%eta_ref( 416 ) = 1.00000000000000

  Mref_V%radius_ref( 417 ) = 4018317.00000000
  Mref_V%density_ref( 417 ) = 5298.10000000000
  Mref_V%vpv_ref( 417 ) = 13224.2700000000
  Mref_V%vsv_ref( 417 ) = 7091.60000000000
  Mref_V%Qkappa_ref( 417 ) = 57822.5000000000
  Mref_V%Qmu_ref( 417 ) = 355.000000000000
  Mref_V%vph_ref( 417 ) = 13224.2700000000
  Mref_V%vsh_ref( 417 ) = 7091.60000000000
  Mref_V%eta_ref( 417 ) = 1.00000000000000

  Mref_V%radius_ref( 418 ) = 4027788.00000000
  Mref_V%density_ref( 418 ) = 5293.38000000000
  Mref_V%vpv_ref( 418 ) = 13213.3900000000
  Mref_V%vsv_ref( 418 ) = 7087.39000000000
  Mref_V%Qkappa_ref( 418 ) = 57822.5000000000
  Mref_V%Qmu_ref( 418 ) = 355.000000000000
  Mref_V%vph_ref( 418 ) = 13213.3900000000
  Mref_V%vsh_ref( 418 ) = 7087.39000000000
  Mref_V%eta_ref( 418 ) = 1.00000000000000

  Mref_V%radius_ref( 419 ) = 4037260.00000000
  Mref_V%density_ref( 419 ) = 5288.65000000000
  Mref_V%vpv_ref( 419 ) = 13202.5100000000
  Mref_V%vsv_ref( 419 ) = 7083.18000000000
  Mref_V%Qkappa_ref( 419 ) = 57822.5000000000
  Mref_V%Qmu_ref( 419 ) = 355.000000000000
  Mref_V%vph_ref( 419 ) = 13202.5100000000
  Mref_V%vsh_ref( 419 ) = 7083.18000000000
  Mref_V%eta_ref( 419 ) = 1.00000000000000

  Mref_V%radius_ref( 420 ) = 4046731.00000000
  Mref_V%density_ref( 420 ) = 5283.92000000000
  Mref_V%vpv_ref( 420 ) = 13191.6400000000
  Mref_V%vsv_ref( 420 ) = 7078.96000000000
  Mref_V%Qkappa_ref( 420 ) = 57822.5000000000
  Mref_V%Qmu_ref( 420 ) = 355.000000000000
  Mref_V%vph_ref( 420 ) = 13191.6400000000
  Mref_V%vsh_ref( 420 ) = 7078.96000000000
  Mref_V%eta_ref( 420 ) = 1.00000000000000

  Mref_V%radius_ref( 421 ) = 4056202.00000000
  Mref_V%density_ref( 421 ) = 5279.19000000000
  Mref_V%vpv_ref( 421 ) = 13180.7800000000
  Mref_V%vsv_ref( 421 ) = 7074.75000000000
  Mref_V%Qkappa_ref( 421 ) = 57822.5000000000
  Mref_V%Qmu_ref( 421 ) = 355.000000000000
  Mref_V%vph_ref( 421 ) = 13180.7800000000
  Mref_V%vsh_ref( 421 ) = 7074.75000000000
  Mref_V%eta_ref( 421 ) = 1.00000000000000

  Mref_V%radius_ref( 422 ) = 4065673.00000000
  Mref_V%density_ref( 422 ) = 5274.45000000000
  Mref_V%vpv_ref( 422 ) = 13169.9100000000
  Mref_V%vsv_ref( 422 ) = 7070.54000000000
  Mref_V%Qkappa_ref( 422 ) = 57822.5000000000
  Mref_V%Qmu_ref( 422 ) = 355.000000000000
  Mref_V%vph_ref( 422 ) = 13169.9100000000
  Mref_V%vsh_ref( 422 ) = 7070.54000000000
  Mref_V%eta_ref( 422 ) = 1.00000000000000

  Mref_V%radius_ref( 423 ) = 4075144.00000000
  Mref_V%density_ref( 423 ) = 5269.72000000000
  Mref_V%vpv_ref( 423 ) = 13159.0500000000
  Mref_V%vsv_ref( 423 ) = 7066.32000000000
  Mref_V%Qkappa_ref( 423 ) = 57822.5000000000
  Mref_V%Qmu_ref( 423 ) = 355.000000000000
  Mref_V%vph_ref( 423 ) = 13159.0500000000
  Mref_V%vsh_ref( 423 ) = 7066.32000000000
  Mref_V%eta_ref( 423 ) = 1.00000000000000

  Mref_V%radius_ref( 424 ) = 4084615.00000000
  Mref_V%density_ref( 424 ) = 5264.98000000000
  Mref_V%vpv_ref( 424 ) = 13148.1900000000
  Mref_V%vsv_ref( 424 ) = 7062.10000000000
  Mref_V%Qkappa_ref( 424 ) = 57822.5000000000
  Mref_V%Qmu_ref( 424 ) = 355.000000000000
  Mref_V%vph_ref( 424 ) = 13148.1900000000
  Mref_V%vsh_ref( 424 ) = 7062.10000000000
  Mref_V%eta_ref( 424 ) = 1.00000000000000

  Mref_V%radius_ref( 425 ) = 4094087.00000000
  Mref_V%density_ref( 425 ) = 5260.25000000000
  Mref_V%vpv_ref( 425 ) = 13137.3300000000
  Mref_V%vsv_ref( 425 ) = 7057.88000000000
  Mref_V%Qkappa_ref( 425 ) = 57822.5000000000
  Mref_V%Qmu_ref( 425 ) = 355.000000000000
  Mref_V%vph_ref( 425 ) = 13137.3300000000
  Mref_V%vsh_ref( 425 ) = 7057.88000000000
  Mref_V%eta_ref( 425 ) = 1.00000000000000

  Mref_V%radius_ref( 426 ) = 4103558.00000000
  Mref_V%density_ref( 426 ) = 5255.51000000000
  Mref_V%vpv_ref( 426 ) = 13126.4700000000
  Mref_V%vsv_ref( 426 ) = 7053.66000000000
  Mref_V%Qkappa_ref( 426 ) = 57822.5000000000
  Mref_V%Qmu_ref( 426 ) = 355.000000000000
  Mref_V%vph_ref( 426 ) = 13126.4700000000
  Mref_V%vsh_ref( 426 ) = 7053.66000000000
  Mref_V%eta_ref( 426 ) = 1.00000000000000

  Mref_V%radius_ref( 427 ) = 4113029.00000000
  Mref_V%density_ref( 427 ) = 5250.77000000000
  Mref_V%vpv_ref( 427 ) = 13115.6100000000
  Mref_V%vsv_ref( 427 ) = 7049.44000000000
  Mref_V%Qkappa_ref( 427 ) = 57822.5000000000
  Mref_V%Qmu_ref( 427 ) = 355.000000000000
  Mref_V%vph_ref( 427 ) = 13115.6100000000
  Mref_V%vsh_ref( 427 ) = 7049.44000000000
  Mref_V%eta_ref( 427 ) = 1.00000000000000

  Mref_V%radius_ref( 428 ) = 4122500.00000000
  Mref_V%density_ref( 428 ) = 5246.02000000000
  Mref_V%vpv_ref( 428 ) = 13104.7500000000
  Mref_V%vsv_ref( 428 ) = 7045.22000000000
  Mref_V%Qkappa_ref( 428 ) = 57822.5000000000
  Mref_V%Qmu_ref( 428 ) = 355.000000000000
  Mref_V%vph_ref( 428 ) = 13104.7500000000
  Mref_V%vsh_ref( 428 ) = 7045.22000000000
  Mref_V%eta_ref( 428 ) = 1.00000000000000

  Mref_V%radius_ref( 429 ) = 4131971.00000000
  Mref_V%density_ref( 429 ) = 5241.28000000000
  Mref_V%vpv_ref( 429 ) = 13093.8900000000
  Mref_V%vsv_ref( 429 ) = 7040.99000000000
  Mref_V%Qkappa_ref( 429 ) = 57822.5000000000
  Mref_V%Qmu_ref( 429 ) = 355.000000000000
  Mref_V%vph_ref( 429 ) = 13093.8900000000
  Mref_V%vsh_ref( 429 ) = 7040.99000000000
  Mref_V%eta_ref( 429 ) = 1.00000000000000

  Mref_V%radius_ref( 430 ) = 4141442.00000000
  Mref_V%density_ref( 430 ) = 5236.53000000000
  Mref_V%vpv_ref( 430 ) = 13083.0200000000
  Mref_V%vsv_ref( 430 ) = 7036.76000000000
  Mref_V%Qkappa_ref( 430 ) = 57822.5000000000
  Mref_V%Qmu_ref( 430 ) = 355.000000000000
  Mref_V%vph_ref( 430 ) = 13083.0200000000
  Mref_V%vsh_ref( 430 ) = 7036.76000000000
  Mref_V%eta_ref( 430 ) = 1.00000000000000

  Mref_V%radius_ref( 431 ) = 4150913.00000000
  Mref_V%density_ref( 431 ) = 5231.78000000000
  Mref_V%vpv_ref( 431 ) = 13072.1600000000
  Mref_V%vsv_ref( 431 ) = 7032.52000000000
  Mref_V%Qkappa_ref( 431 ) = 57822.5000000000
  Mref_V%Qmu_ref( 431 ) = 355.000000000000
  Mref_V%vph_ref( 431 ) = 13072.1600000000
  Mref_V%vsh_ref( 431 ) = 7032.52000000000
  Mref_V%eta_ref( 431 ) = 1.00000000000000

  Mref_V%radius_ref( 432 ) = 4160385.00000000
  Mref_V%density_ref( 432 ) = 5227.03000000000
  Mref_V%vpv_ref( 432 ) = 13061.2900000000
  Mref_V%vsv_ref( 432 ) = 7028.29000000000
  Mref_V%Qkappa_ref( 432 ) = 57822.5000000000
  Mref_V%Qmu_ref( 432 ) = 355.000000000000
  Mref_V%vph_ref( 432 ) = 13061.2900000000
  Mref_V%vsh_ref( 432 ) = 7028.29000000000
  Mref_V%eta_ref( 432 ) = 1.00000000000000

  Mref_V%radius_ref( 433 ) = 4169856.00000000
  Mref_V%density_ref( 433 ) = 5222.27000000000
  Mref_V%vpv_ref( 433 ) = 13050.4200000000
  Mref_V%vsv_ref( 433 ) = 7024.05000000000
  Mref_V%Qkappa_ref( 433 ) = 57822.5000000000
  Mref_V%Qmu_ref( 433 ) = 355.000000000000
  Mref_V%vph_ref( 433 ) = 13050.4200000000
  Mref_V%vsh_ref( 433 ) = 7024.05000000000
  Mref_V%eta_ref( 433 ) = 1.00000000000000

  Mref_V%radius_ref( 434 ) = 4179327.00000000
  Mref_V%density_ref( 434 ) = 5217.52000000000
  Mref_V%vpv_ref( 434 ) = 13039.5500000000
  Mref_V%vsv_ref( 434 ) = 7019.81000000000
  Mref_V%Qkappa_ref( 434 ) = 57822.5000000000
  Mref_V%Qmu_ref( 434 ) = 355.000000000000
  Mref_V%vph_ref( 434 ) = 13039.5500000000
  Mref_V%vsh_ref( 434 ) = 7019.81000000000
  Mref_V%eta_ref( 434 ) = 1.00000000000000

  Mref_V%radius_ref( 435 ) = 4188798.00000000
  Mref_V%density_ref( 435 ) = 5212.76000000000
  Mref_V%vpv_ref( 435 ) = 13028.6700000000
  Mref_V%vsv_ref( 435 ) = 7015.56000000000
  Mref_V%Qkappa_ref( 435 ) = 57822.5000000000
  Mref_V%Qmu_ref( 435 ) = 355.000000000000
  Mref_V%vph_ref( 435 ) = 13028.6700000000
  Mref_V%vsh_ref( 435 ) = 7015.56000000000
  Mref_V%eta_ref( 435 ) = 1.00000000000000

  Mref_V%radius_ref( 436 ) = 4198269.00000000
  Mref_V%density_ref( 436 ) = 5208.00000000000
  Mref_V%vpv_ref( 436 ) = 13017.7800000000
  Mref_V%vsv_ref( 436 ) = 7011.31000000000
  Mref_V%Qkappa_ref( 436 ) = 57822.5000000000
  Mref_V%Qmu_ref( 436 ) = 355.000000000000
  Mref_V%vph_ref( 436 ) = 13017.7800000000
  Mref_V%vsh_ref( 436 ) = 7011.31000000000
  Mref_V%eta_ref( 436 ) = 1.00000000000000

  Mref_V%radius_ref( 437 ) = 4207740.00000000
  Mref_V%density_ref( 437 ) = 5203.23000000000
  Mref_V%vpv_ref( 437 ) = 13006.9000000000
  Mref_V%vsv_ref( 437 ) = 7007.06000000000
  Mref_V%Qkappa_ref( 437 ) = 57822.5000000000
  Mref_V%Qmu_ref( 437 ) = 355.000000000000
  Mref_V%vph_ref( 437 ) = 13006.9000000000
  Mref_V%vsh_ref( 437 ) = 7007.06000000000
  Mref_V%eta_ref( 437 ) = 1.00000000000000

  Mref_V%radius_ref( 438 ) = 4217212.00000000
  Mref_V%density_ref( 438 ) = 5198.47000000000
  Mref_V%vpv_ref( 438 ) = 12996.0000000000
  Mref_V%vsv_ref( 438 ) = 7002.80000000000
  Mref_V%Qkappa_ref( 438 ) = 57822.5000000000
  Mref_V%Qmu_ref( 438 ) = 355.000000000000
  Mref_V%vph_ref( 438 ) = 12996.0000000000
  Mref_V%vsh_ref( 438 ) = 7002.80000000000
  Mref_V%eta_ref( 438 ) = 1.00000000000000

  Mref_V%radius_ref( 439 ) = 4226683.00000000
  Mref_V%density_ref( 439 ) = 5193.70000000000
  Mref_V%vpv_ref( 439 ) = 12985.1000000000
  Mref_V%vsv_ref( 439 ) = 6998.54000000000
  Mref_V%Qkappa_ref( 439 ) = 57822.5000000000
  Mref_V%Qmu_ref( 439 ) = 355.000000000000
  Mref_V%vph_ref( 439 ) = 12985.1000000000
  Mref_V%vsh_ref( 439 ) = 6998.54000000000
  Mref_V%eta_ref( 439 ) = 1.00000000000000

  Mref_V%radius_ref( 440 ) = 4236154.00000000
  Mref_V%density_ref( 440 ) = 5188.93000000000
  Mref_V%vpv_ref( 440 ) = 12974.1900000000
  Mref_V%vsv_ref( 440 ) = 6994.27000000000
  Mref_V%Qkappa_ref( 440 ) = 57822.5000000000
  Mref_V%Qmu_ref( 440 ) = 355.000000000000
  Mref_V%vph_ref( 440 ) = 12974.1900000000
  Mref_V%vsh_ref( 440 ) = 6994.27000000000
  Mref_V%eta_ref( 440 ) = 1.00000000000000

  Mref_V%radius_ref( 441 ) = 4245625.00000000
  Mref_V%density_ref( 441 ) = 5184.15000000000
  Mref_V%vpv_ref( 441 ) = 12963.2800000000
  Mref_V%vsv_ref( 441 ) = 6990.00000000000
  Mref_V%Qkappa_ref( 441 ) = 57822.5000000000
  Mref_V%Qmu_ref( 441 ) = 355.000000000000
  Mref_V%vph_ref( 441 ) = 12963.2800000000
  Mref_V%vsh_ref( 441 ) = 6990.00000000000
  Mref_V%eta_ref( 441 ) = 1.00000000000000

  Mref_V%radius_ref( 442 ) = 4255096.00000000
  Mref_V%density_ref( 442 ) = 5179.38000000000
  Mref_V%vpv_ref( 442 ) = 12952.3600000000
  Mref_V%vsv_ref( 442 ) = 6985.72000000000
  Mref_V%Qkappa_ref( 442 ) = 57822.5000000000
  Mref_V%Qmu_ref( 442 ) = 355.000000000000
  Mref_V%vph_ref( 442 ) = 12952.3600000000
  Mref_V%vsh_ref( 442 ) = 6985.72000000000
  Mref_V%eta_ref( 442 ) = 1.00000000000000

  Mref_V%radius_ref( 443 ) = 4264567.00000000
  Mref_V%density_ref( 443 ) = 5174.60000000000
  Mref_V%vpv_ref( 443 ) = 12941.4200000000
  Mref_V%vsv_ref( 443 ) = 6981.44000000000
  Mref_V%Qkappa_ref( 443 ) = 57822.5000000000
  Mref_V%Qmu_ref( 443 ) = 355.000000000000
  Mref_V%vph_ref( 443 ) = 12941.4200000000
  Mref_V%vsh_ref( 443 ) = 6981.44000000000
  Mref_V%eta_ref( 443 ) = 1.00000000000000

  Mref_V%radius_ref( 444 ) = 4274038.00000000
  Mref_V%density_ref( 444 ) = 5169.82000000000
  Mref_V%vpv_ref( 444 ) = 12930.4800000000
  Mref_V%vsv_ref( 444 ) = 6977.15000000000
  Mref_V%Qkappa_ref( 444 ) = 57822.5000000000
  Mref_V%Qmu_ref( 444 ) = 355.000000000000
  Mref_V%vph_ref( 444 ) = 12930.4800000000
  Mref_V%vsh_ref( 444 ) = 6977.15000000000
  Mref_V%eta_ref( 444 ) = 1.00000000000000

  Mref_V%radius_ref( 445 ) = 4283510.00000000
  Mref_V%density_ref( 445 ) = 5165.03000000000
  Mref_V%vpv_ref( 445 ) = 12919.5400000000
  Mref_V%vsv_ref( 445 ) = 6972.86000000000
  Mref_V%Qkappa_ref( 445 ) = 57822.5000000000
  Mref_V%Qmu_ref( 445 ) = 355.000000000000
  Mref_V%vph_ref( 445 ) = 12919.5400000000
  Mref_V%vsh_ref( 445 ) = 6972.86000000000
  Mref_V%eta_ref( 445 ) = 1.00000000000000

  Mref_V%radius_ref( 446 ) = 4292981.00000000
  Mref_V%density_ref( 446 ) = 5160.24000000000
  Mref_V%vpv_ref( 446 ) = 12908.5800000000
  Mref_V%vsv_ref( 446 ) = 6968.57000000000
  Mref_V%Qkappa_ref( 446 ) = 57822.5000000000
  Mref_V%Qmu_ref( 446 ) = 355.000000000000
  Mref_V%vph_ref( 446 ) = 12908.5800000000
  Mref_V%vsh_ref( 446 ) = 6968.57000000000
  Mref_V%eta_ref( 446 ) = 1.00000000000000

  Mref_V%radius_ref( 447 ) = 4302452.00000000
  Mref_V%density_ref( 447 ) = 5155.45000000000
  Mref_V%vpv_ref( 447 ) = 12897.6100000000
  Mref_V%vsv_ref( 447 ) = 6964.26000000000
  Mref_V%Qkappa_ref( 447 ) = 57822.5000000000
  Mref_V%Qmu_ref( 447 ) = 355.000000000000
  Mref_V%vph_ref( 447 ) = 12897.6100000000
  Mref_V%vsh_ref( 447 ) = 6964.26000000000
  Mref_V%eta_ref( 447 ) = 1.00000000000000

  Mref_V%radius_ref( 448 ) = 4311923.00000000
  Mref_V%density_ref( 448 ) = 5150.65000000000
  Mref_V%vpv_ref( 448 ) = 12886.6300000000
  Mref_V%vsv_ref( 448 ) = 6959.95000000000
  Mref_V%Qkappa_ref( 448 ) = 57822.5000000000
  Mref_V%Qmu_ref( 448 ) = 355.000000000000
  Mref_V%vph_ref( 448 ) = 12886.6300000000
  Mref_V%vsh_ref( 448 ) = 6959.95000000000
  Mref_V%eta_ref( 448 ) = 1.00000000000000

  Mref_V%radius_ref( 449 ) = 4321394.00000000
  Mref_V%density_ref( 449 ) = 5145.86000000000
  Mref_V%vpv_ref( 449 ) = 12875.6300000000
  Mref_V%vsv_ref( 449 ) = 6955.64000000000
  Mref_V%Qkappa_ref( 449 ) = 57822.5000000000
  Mref_V%Qmu_ref( 449 ) = 355.000000000000
  Mref_V%vph_ref( 449 ) = 12875.6300000000
  Mref_V%vsh_ref( 449 ) = 6955.64000000000
  Mref_V%eta_ref( 449 ) = 1.00000000000000

  Mref_V%radius_ref( 450 ) = 4330865.00000000
  Mref_V%density_ref( 450 ) = 5141.06000000000
  Mref_V%vpv_ref( 450 ) = 12864.6300000000
  Mref_V%vsv_ref( 450 ) = 6951.32000000000
  Mref_V%Qkappa_ref( 450 ) = 57822.5000000000
  Mref_V%Qmu_ref( 450 ) = 355.000000000000
  Mref_V%vph_ref( 450 ) = 12864.6300000000
  Mref_V%vsh_ref( 450 ) = 6951.32000000000
  Mref_V%eta_ref( 450 ) = 1.00000000000000

  Mref_V%radius_ref( 451 ) = 4340337.00000000
  Mref_V%density_ref( 451 ) = 5136.25000000000
  Mref_V%vpv_ref( 451 ) = 12853.6100000000
  Mref_V%vsv_ref( 451 ) = 6946.99000000000
  Mref_V%Qkappa_ref( 451 ) = 57822.5000000000
  Mref_V%Qmu_ref( 451 ) = 355.000000000000
  Mref_V%vph_ref( 451 ) = 12853.6100000000
  Mref_V%vsh_ref( 451 ) = 6946.99000000000
  Mref_V%eta_ref( 451 ) = 1.00000000000000

  Mref_V%radius_ref( 452 ) = 4349808.00000000
  Mref_V%density_ref( 452 ) = 5131.44000000000
  Mref_V%vpv_ref( 452 ) = 12842.5800000000
  Mref_V%vsv_ref( 452 ) = 6942.66000000000
  Mref_V%Qkappa_ref( 452 ) = 57822.5000000000
  Mref_V%Qmu_ref( 452 ) = 355.000000000000
  Mref_V%vph_ref( 452 ) = 12842.5800000000
  Mref_V%vsh_ref( 452 ) = 6942.66000000000
  Mref_V%eta_ref( 452 ) = 1.00000000000000

  Mref_V%radius_ref( 453 ) = 4359279.00000000
  Mref_V%density_ref( 453 ) = 5126.63000000000
  Mref_V%vpv_ref( 453 ) = 12831.5400000000
  Mref_V%vsv_ref( 453 ) = 6938.31000000000
  Mref_V%Qkappa_ref( 453 ) = 57822.5000000000
  Mref_V%Qmu_ref( 453 ) = 355.000000000000
  Mref_V%vph_ref( 453 ) = 12831.5400000000
  Mref_V%vsh_ref( 453 ) = 6938.31000000000
  Mref_V%eta_ref( 453 ) = 1.00000000000000

  Mref_V%radius_ref( 454 ) = 4368750.00000000
  Mref_V%density_ref( 454 ) = 5121.82000000000
  Mref_V%vpv_ref( 454 ) = 12820.4800000000
  Mref_V%vsv_ref( 454 ) = 6933.97000000000
  Mref_V%Qkappa_ref( 454 ) = 57822.5000000000
  Mref_V%Qmu_ref( 454 ) = 355.000000000000
  Mref_V%vph_ref( 454 ) = 12820.4800000000
  Mref_V%vsh_ref( 454 ) = 6933.97000000000
  Mref_V%eta_ref( 454 ) = 1.00000000000000

  Mref_V%radius_ref( 455 ) = 4378221.00000000
  Mref_V%density_ref( 455 ) = 5117.00000000000
  Mref_V%vpv_ref( 455 ) = 12809.4100000000
  Mref_V%vsv_ref( 455 ) = 6929.61000000000
  Mref_V%Qkappa_ref( 455 ) = 57822.5000000000
  Mref_V%Qmu_ref( 455 ) = 355.000000000000
  Mref_V%vph_ref( 455 ) = 12809.4100000000
  Mref_V%vsh_ref( 455 ) = 6929.61000000000
  Mref_V%eta_ref( 455 ) = 1.00000000000000

  Mref_V%radius_ref( 456 ) = 4387692.00000000
  Mref_V%density_ref( 456 ) = 5112.18000000000
  Mref_V%vpv_ref( 456 ) = 12798.3200000000
  Mref_V%vsv_ref( 456 ) = 6925.25000000000
  Mref_V%Qkappa_ref( 456 ) = 57822.5000000000
  Mref_V%Qmu_ref( 456 ) = 355.000000000000
  Mref_V%vph_ref( 456 ) = 12798.3200000000
  Mref_V%vsh_ref( 456 ) = 6925.25000000000
  Mref_V%eta_ref( 456 ) = 1.00000000000000

  Mref_V%radius_ref( 457 ) = 4397163.00000000
  Mref_V%density_ref( 457 ) = 5107.35000000000
  Mref_V%vpv_ref( 457 ) = 12787.2200000000
  Mref_V%vsv_ref( 457 ) = 6920.88000000000
  Mref_V%Qkappa_ref( 457 ) = 57822.5000000000
  Mref_V%Qmu_ref( 457 ) = 355.000000000000
  Mref_V%vph_ref( 457 ) = 12787.2200000000
  Mref_V%vsh_ref( 457 ) = 6920.88000000000
  Mref_V%eta_ref( 457 ) = 1.00000000000000

  Mref_V%radius_ref( 458 ) = 4406635.00000000
  Mref_V%density_ref( 458 ) = 5102.52000000000
  Mref_V%vpv_ref( 458 ) = 12776.1000000000
  Mref_V%vsv_ref( 458 ) = 6916.50000000000
  Mref_V%Qkappa_ref( 458 ) = 57822.5000000000
  Mref_V%Qmu_ref( 458 ) = 355.000000000000
  Mref_V%vph_ref( 458 ) = 12776.1000000000
  Mref_V%vsh_ref( 458 ) = 6916.50000000000
  Mref_V%eta_ref( 458 ) = 1.00000000000000

  Mref_V%radius_ref( 459 ) = 4416106.00000000
  Mref_V%density_ref( 459 ) = 5097.69000000000
  Mref_V%vpv_ref( 459 ) = 12764.9600000000
  Mref_V%vsv_ref( 459 ) = 6912.11000000000
  Mref_V%Qkappa_ref( 459 ) = 57822.5000000000
  Mref_V%Qmu_ref( 459 ) = 355.000000000000
  Mref_V%vph_ref( 459 ) = 12764.9600000000
  Mref_V%vsh_ref( 459 ) = 6912.11000000000
  Mref_V%eta_ref( 459 ) = 1.00000000000000

  Mref_V%radius_ref( 460 ) = 4425577.00000000
  Mref_V%density_ref( 460 ) = 5092.85000000000
  Mref_V%vpv_ref( 460 ) = 12753.8100000000
  Mref_V%vsv_ref( 460 ) = 6907.72000000000
  Mref_V%Qkappa_ref( 460 ) = 57822.5000000000
  Mref_V%Qmu_ref( 460 ) = 355.000000000000
  Mref_V%vph_ref( 460 ) = 12753.8100000000
  Mref_V%vsh_ref( 460 ) = 6907.72000000000
  Mref_V%eta_ref( 460 ) = 1.00000000000000

  Mref_V%radius_ref( 461 ) = 4435048.00000000
  Mref_V%density_ref( 461 ) = 5088.01000000000
  Mref_V%vpv_ref( 461 ) = 12742.6300000000
  Mref_V%vsv_ref( 461 ) = 6903.32000000000
  Mref_V%Qkappa_ref( 461 ) = 57822.5000000000
  Mref_V%Qmu_ref( 461 ) = 355.000000000000
  Mref_V%vph_ref( 461 ) = 12742.6300000000
  Mref_V%vsh_ref( 461 ) = 6903.32000000000
  Mref_V%eta_ref( 461 ) = 1.00000000000000

  Mref_V%radius_ref( 462 ) = 4444519.00000000
  Mref_V%density_ref( 462 ) = 5083.16000000000
  Mref_V%vpv_ref( 462 ) = 12731.4400000000
  Mref_V%vsv_ref( 462 ) = 6898.91000000000
  Mref_V%Qkappa_ref( 462 ) = 57822.5000000000
  Mref_V%Qmu_ref( 462 ) = 355.000000000000
  Mref_V%vph_ref( 462 ) = 12731.4400000000
  Mref_V%vsh_ref( 462 ) = 6898.91000000000
  Mref_V%eta_ref( 462 ) = 1.00000000000000

  Mref_V%radius_ref( 463 ) = 4453990.00000000
  Mref_V%density_ref( 463 ) = 5078.31000000000
  Mref_V%vpv_ref( 463 ) = 12720.2400000000
  Mref_V%vsv_ref( 463 ) = 6894.49000000000
  Mref_V%Qkappa_ref( 463 ) = 57822.5000000000
  Mref_V%Qmu_ref( 463 ) = 355.000000000000
  Mref_V%vph_ref( 463 ) = 12720.2400000000
  Mref_V%vsh_ref( 463 ) = 6894.49000000000
  Mref_V%eta_ref( 463 ) = 1.00000000000000

  Mref_V%radius_ref( 464 ) = 4463462.00000000
  Mref_V%density_ref( 464 ) = 5073.46000000000
  Mref_V%vpv_ref( 464 ) = 12709.0100000000
  Mref_V%vsv_ref( 464 ) = 6890.06000000000
  Mref_V%Qkappa_ref( 464 ) = 57822.5000000000
  Mref_V%Qmu_ref( 464 ) = 355.000000000000
  Mref_V%vph_ref( 464 ) = 12709.0100000000
  Mref_V%vsh_ref( 464 ) = 6890.06000000000
  Mref_V%eta_ref( 464 ) = 1.00000000000000

  Mref_V%radius_ref( 465 ) = 4472933.00000000
  Mref_V%density_ref( 465 ) = 5068.60000000000
  Mref_V%vpv_ref( 465 ) = 12697.7600000000
  Mref_V%vsv_ref( 465 ) = 6885.62000000000
  Mref_V%Qkappa_ref( 465 ) = 57822.5000000000
  Mref_V%Qmu_ref( 465 ) = 355.000000000000
  Mref_V%vph_ref( 465 ) = 12697.7600000000
  Mref_V%vsh_ref( 465 ) = 6885.62000000000
  Mref_V%eta_ref( 465 ) = 1.00000000000000

  Mref_V%radius_ref( 466 ) = 4482404.00000000
  Mref_V%density_ref( 466 ) = 5063.74000000000
  Mref_V%vpv_ref( 466 ) = 12686.4900000000
  Mref_V%vsv_ref( 466 ) = 6881.17000000000
  Mref_V%Qkappa_ref( 466 ) = 57822.5000000000
  Mref_V%Qmu_ref( 466 ) = 355.000000000000
  Mref_V%vph_ref( 466 ) = 12686.4900000000
  Mref_V%vsh_ref( 466 ) = 6881.17000000000
  Mref_V%eta_ref( 466 ) = 1.00000000000000

  Mref_V%radius_ref( 467 ) = 4491875.00000000
  Mref_V%density_ref( 467 ) = 5058.87000000000
  Mref_V%vpv_ref( 467 ) = 12675.2000000000
  Mref_V%vsv_ref( 467 ) = 6876.72000000000
  Mref_V%Qkappa_ref( 467 ) = 57822.5000000000
  Mref_V%Qmu_ref( 467 ) = 355.000000000000
  Mref_V%vph_ref( 467 ) = 12675.2000000000
  Mref_V%vsh_ref( 467 ) = 6876.72000000000
  Mref_V%eta_ref( 467 ) = 1.00000000000000

  Mref_V%radius_ref( 468 ) = 4501346.00000000
  Mref_V%density_ref( 468 ) = 5054.00000000000
  Mref_V%vpv_ref( 468 ) = 12663.8900000000
  Mref_V%vsv_ref( 468 ) = 6872.25000000000
  Mref_V%Qkappa_ref( 468 ) = 57822.5000000000
  Mref_V%Qmu_ref( 468 ) = 355.000000000000
  Mref_V%vph_ref( 468 ) = 12663.8900000000
  Mref_V%vsh_ref( 468 ) = 6872.25000000000
  Mref_V%eta_ref( 468 ) = 1.00000000000000

  Mref_V%radius_ref( 469 ) = 4510817.00000000
  Mref_V%density_ref( 469 ) = 5049.13000000000
  Mref_V%vpv_ref( 469 ) = 12652.5600000000
  Mref_V%vsv_ref( 469 ) = 6867.78000000000
  Mref_V%Qkappa_ref( 469 ) = 57822.5000000000
  Mref_V%Qmu_ref( 469 ) = 355.000000000000
  Mref_V%vph_ref( 469 ) = 12652.5600000000
  Mref_V%vsh_ref( 469 ) = 6867.78000000000
  Mref_V%eta_ref( 469 ) = 1.00000000000000

  Mref_V%radius_ref( 470 ) = 4520288.00000000
  Mref_V%density_ref( 470 ) = 5044.25000000000
  Mref_V%vpv_ref( 470 ) = 12641.2000000000
  Mref_V%vsv_ref( 470 ) = 6863.29000000000
  Mref_V%Qkappa_ref( 470 ) = 57822.5000000000
  Mref_V%Qmu_ref( 470 ) = 355.000000000000
  Mref_V%vph_ref( 470 ) = 12641.2000000000
  Mref_V%vsh_ref( 470 ) = 6863.29000000000
  Mref_V%eta_ref( 470 ) = 1.00000000000000

  Mref_V%radius_ref( 471 ) = 4529760.00000000
  Mref_V%density_ref( 471 ) = 5039.36000000000
  Mref_V%vpv_ref( 471 ) = 12629.8200000000
  Mref_V%vsv_ref( 471 ) = 6858.80000000000
  Mref_V%Qkappa_ref( 471 ) = 57822.5000000000
  Mref_V%Qmu_ref( 471 ) = 355.000000000000
  Mref_V%vph_ref( 471 ) = 12629.8200000000
  Mref_V%vsh_ref( 471 ) = 6858.80000000000
  Mref_V%eta_ref( 471 ) = 1.00000000000000

  Mref_V%radius_ref( 472 ) = 4539231.00000000
  Mref_V%density_ref( 472 ) = 5034.47000000000
  Mref_V%vpv_ref( 472 ) = 12618.4200000000
  Mref_V%vsv_ref( 472 ) = 6854.29000000000
  Mref_V%Qkappa_ref( 472 ) = 57822.5000000000
  Mref_V%Qmu_ref( 472 ) = 355.000000000000
  Mref_V%vph_ref( 472 ) = 12618.4200000000
  Mref_V%vsh_ref( 472 ) = 6854.29000000000
  Mref_V%eta_ref( 472 ) = 1.00000000000000

  Mref_V%radius_ref( 473 ) = 4548702.00000000
  Mref_V%density_ref( 473 ) = 5029.58000000000
  Mref_V%vpv_ref( 473 ) = 12606.9900000000
  Mref_V%vsv_ref( 473 ) = 6849.78000000000
  Mref_V%Qkappa_ref( 473 ) = 57822.5000000000
  Mref_V%Qmu_ref( 473 ) = 355.000000000000
  Mref_V%vph_ref( 473 ) = 12606.9900000000
  Mref_V%vsh_ref( 473 ) = 6849.78000000000
  Mref_V%eta_ref( 473 ) = 1.00000000000000

  Mref_V%radius_ref( 474 ) = 4558173.00000000
  Mref_V%density_ref( 474 ) = 5024.68000000000
  Mref_V%vpv_ref( 474 ) = 12595.5400000000
  Mref_V%vsv_ref( 474 ) = 6845.25000000000
  Mref_V%Qkappa_ref( 474 ) = 57822.5000000000
  Mref_V%Qmu_ref( 474 ) = 355.000000000000
  Mref_V%vph_ref( 474 ) = 12595.5400000000
  Mref_V%vsh_ref( 474 ) = 6845.25000000000
  Mref_V%eta_ref( 474 ) = 1.00000000000000

  Mref_V%radius_ref( 475 ) = 4567644.00000000
  Mref_V%density_ref( 475 ) = 5019.78000000000
  Mref_V%vpv_ref( 475 ) = 12584.0600000000
  Mref_V%vsv_ref( 475 ) = 6840.71000000000
  Mref_V%Qkappa_ref( 475 ) = 57822.5000000000
  Mref_V%Qmu_ref( 475 ) = 355.000000000000
  Mref_V%vph_ref( 475 ) = 12584.0600000000
  Mref_V%vsh_ref( 475 ) = 6840.71000000000
  Mref_V%eta_ref( 475 ) = 1.00000000000000

  Mref_V%radius_ref( 476 ) = 4577115.00000000
  Mref_V%density_ref( 476 ) = 5014.87000000000
  Mref_V%vpv_ref( 476 ) = 12572.5600000000
  Mref_V%vsv_ref( 476 ) = 6836.16000000000
  Mref_V%Qkappa_ref( 476 ) = 57822.5000000000
  Mref_V%Qmu_ref( 476 ) = 355.000000000000
  Mref_V%vph_ref( 476 ) = 12572.5600000000
  Mref_V%vsh_ref( 476 ) = 6836.16000000000
  Mref_V%eta_ref( 476 ) = 1.00000000000000

  Mref_V%radius_ref( 477 ) = 4586587.00000000
  Mref_V%density_ref( 477 ) = 5009.96000000000
  Mref_V%vpv_ref( 477 ) = 12561.0300000000
  Mref_V%vsv_ref( 477 ) = 6831.60000000000
  Mref_V%Qkappa_ref( 477 ) = 57822.5000000000
  Mref_V%Qmu_ref( 477 ) = 355.000000000000
  Mref_V%vph_ref( 477 ) = 12561.0300000000
  Mref_V%vsh_ref( 477 ) = 6831.60000000000
  Mref_V%eta_ref( 477 ) = 1.00000000000000

  Mref_V%radius_ref( 478 ) = 4596058.00000000
  Mref_V%density_ref( 478 ) = 5005.04000000000
  Mref_V%vpv_ref( 478 ) = 12549.4800000000
  Mref_V%vsv_ref( 478 ) = 6827.03000000000
  Mref_V%Qkappa_ref( 478 ) = 57822.5000000000
  Mref_V%Qmu_ref( 478 ) = 355.000000000000
  Mref_V%vph_ref( 478 ) = 12549.4800000000
  Mref_V%vsh_ref( 478 ) = 6827.03000000000
  Mref_V%eta_ref( 478 ) = 1.00000000000000

  Mref_V%radius_ref( 479 ) = 4605529.00000000
  Mref_V%density_ref( 479 ) = 5000.12000000000
  Mref_V%vpv_ref( 479 ) = 12537.8900000000
  Mref_V%vsv_ref( 479 ) = 6822.45000000000
  Mref_V%Qkappa_ref( 479 ) = 57822.5000000000
  Mref_V%Qmu_ref( 479 ) = 355.000000000000
  Mref_V%vph_ref( 479 ) = 12537.8900000000
  Mref_V%vsh_ref( 479 ) = 6822.45000000000
  Mref_V%eta_ref( 479 ) = 1.00000000000000

  Mref_V%radius_ref( 480 ) = 4615000.00000000
  Mref_V%density_ref( 480 ) = 4995.19000000000
  Mref_V%vpv_ref( 480 ) = 12526.2800000000
  Mref_V%vsv_ref( 480 ) = 6817.85000000000
  Mref_V%Qkappa_ref( 480 ) = 57822.5000000000
  Mref_V%Qmu_ref( 480 ) = 355.000000000000
  Mref_V%vph_ref( 480 ) = 12526.2800000000
  Mref_V%vsh_ref( 480 ) = 6817.85000000000
  Mref_V%eta_ref( 480 ) = 1.00000000000000

  Mref_V%radius_ref( 481 ) = 4624471.00000000
  Mref_V%density_ref( 481 ) = 4990.26000000000
  Mref_V%vpv_ref( 481 ) = 12514.6400000000
  Mref_V%vsv_ref( 481 ) = 6813.25000000000
  Mref_V%Qkappa_ref( 481 ) = 57822.5000000000
  Mref_V%Qmu_ref( 481 ) = 355.000000000000
  Mref_V%vph_ref( 481 ) = 12514.6400000000
  Mref_V%vsh_ref( 481 ) = 6813.25000000000
  Mref_V%eta_ref( 481 ) = 1.00000000000000

  Mref_V%radius_ref( 482 ) = 4633942.00000000
  Mref_V%density_ref( 482 ) = 4985.32000000000
  Mref_V%vpv_ref( 482 ) = 12502.9800000000
  Mref_V%vsv_ref( 482 ) = 6808.63000000000
  Mref_V%Qkappa_ref( 482 ) = 57822.5000000000
  Mref_V%Qmu_ref( 482 ) = 355.000000000000
  Mref_V%vph_ref( 482 ) = 12502.9800000000
  Mref_V%vsh_ref( 482 ) = 6808.63000000000
  Mref_V%eta_ref( 482 ) = 1.00000000000000

  Mref_V%radius_ref( 483 ) = 4643413.00000000
  Mref_V%density_ref( 483 ) = 4980.38000000000
  Mref_V%vpv_ref( 483 ) = 12491.2800000000
  Mref_V%vsv_ref( 483 ) = 6804.00000000000
  Mref_V%Qkappa_ref( 483 ) = 57822.5000000000
  Mref_V%Qmu_ref( 483 ) = 355.000000000000
  Mref_V%vph_ref( 483 ) = 12491.2800000000
  Mref_V%vsh_ref( 483 ) = 6804.00000000000
  Mref_V%eta_ref( 483 ) = 1.00000000000000

  Mref_V%radius_ref( 484 ) = 4652885.00000000
  Mref_V%density_ref( 484 ) = 4975.43000000000
  Mref_V%vpv_ref( 484 ) = 12479.5500000000
  Mref_V%vsv_ref( 484 ) = 6799.35000000000
  Mref_V%Qkappa_ref( 484 ) = 57822.5000000000
  Mref_V%Qmu_ref( 484 ) = 355.000000000000
  Mref_V%vph_ref( 484 ) = 12479.5500000000
  Mref_V%vsh_ref( 484 ) = 6799.35000000000
  Mref_V%eta_ref( 484 ) = 1.00000000000000

  Mref_V%radius_ref( 485 ) = 4662356.00000000
  Mref_V%density_ref( 485 ) = 4970.47000000000
  Mref_V%vpv_ref( 485 ) = 12467.7900000000
  Mref_V%vsv_ref( 485 ) = 6794.70000000000
  Mref_V%Qkappa_ref( 485 ) = 57822.5000000000
  Mref_V%Qmu_ref( 485 ) = 355.000000000000
  Mref_V%vph_ref( 485 ) = 12467.7900000000
  Mref_V%vsh_ref( 485 ) = 6794.70000000000
  Mref_V%eta_ref( 485 ) = 1.00000000000000

  Mref_V%radius_ref( 486 ) = 4671827.00000000
  Mref_V%density_ref( 486 ) = 4965.51000000000
  Mref_V%vpv_ref( 486 ) = 12456.0100000000
  Mref_V%vsv_ref( 486 ) = 6790.03000000000
  Mref_V%Qkappa_ref( 486 ) = 57822.5000000000
  Mref_V%Qmu_ref( 486 ) = 355.000000000000
  Mref_V%vph_ref( 486 ) = 12456.0100000000
  Mref_V%vsh_ref( 486 ) = 6790.03000000000
  Mref_V%eta_ref( 486 ) = 1.00000000000000

  Mref_V%radius_ref( 487 ) = 4681298.00000000
  Mref_V%density_ref( 487 ) = 4960.55000000000
  Mref_V%vpv_ref( 487 ) = 12444.1900000000
  Mref_V%vsv_ref( 487 ) = 6785.34000000000
  Mref_V%Qkappa_ref( 487 ) = 57822.5000000000
  Mref_V%Qmu_ref( 487 ) = 355.000000000000
  Mref_V%vph_ref( 487 ) = 12444.1900000000
  Mref_V%vsh_ref( 487 ) = 6785.34000000000
  Mref_V%eta_ref( 487 ) = 1.00000000000000

  Mref_V%radius_ref( 488 ) = 4690769.00000000
  Mref_V%density_ref( 488 ) = 4955.58000000000
  Mref_V%vpv_ref( 488 ) = 12432.3300000000
  Mref_V%vsv_ref( 488 ) = 6780.65000000000
  Mref_V%Qkappa_ref( 488 ) = 57822.5000000000
  Mref_V%Qmu_ref( 488 ) = 355.000000000000
  Mref_V%vph_ref( 488 ) = 12432.3300000000
  Mref_V%vsh_ref( 488 ) = 6780.65000000000
  Mref_V%eta_ref( 488 ) = 1.00000000000000

  Mref_V%radius_ref( 489 ) = 4700240.00000000
  Mref_V%density_ref( 489 ) = 4950.60000000000
  Mref_V%vpv_ref( 489 ) = 12420.4500000000
  Mref_V%vsv_ref( 489 ) = 6775.94000000000
  Mref_V%Qkappa_ref( 489 ) = 57822.5000000000
  Mref_V%Qmu_ref( 489 ) = 355.000000000000
  Mref_V%vph_ref( 489 ) = 12420.4500000000
  Mref_V%vsh_ref( 489 ) = 6775.94000000000
  Mref_V%eta_ref( 489 ) = 1.00000000000000

  Mref_V%radius_ref( 490 ) = 4709712.00000000
  Mref_V%density_ref( 490 ) = 4945.62000000000
  Mref_V%vpv_ref( 490 ) = 12408.5300000000
  Mref_V%vsv_ref( 490 ) = 6771.22000000000
  Mref_V%Qkappa_ref( 490 ) = 57822.5000000000
  Mref_V%Qmu_ref( 490 ) = 355.000000000000
  Mref_V%vph_ref( 490 ) = 12408.5300000000
  Mref_V%vsh_ref( 490 ) = 6771.22000000000
  Mref_V%eta_ref( 490 ) = 1.00000000000000

  Mref_V%radius_ref( 491 ) = 4719183.00000000
  Mref_V%density_ref( 491 ) = 4940.63000000000
  Mref_V%vpv_ref( 491 ) = 12396.5800000000
  Mref_V%vsv_ref( 491 ) = 6766.48000000000
  Mref_V%Qkappa_ref( 491 ) = 57822.5000000000
  Mref_V%Qmu_ref( 491 ) = 355.000000000000
  Mref_V%vph_ref( 491 ) = 12396.5800000000
  Mref_V%vsh_ref( 491 ) = 6766.48000000000
  Mref_V%eta_ref( 491 ) = 1.00000000000000

  Mref_V%radius_ref( 492 ) = 4728654.00000000
  Mref_V%density_ref( 492 ) = 4935.64000000000
  Mref_V%vpv_ref( 492 ) = 12384.6000000000
  Mref_V%vsv_ref( 492 ) = 6761.73000000000
  Mref_V%Qkappa_ref( 492 ) = 57822.5000000000
  Mref_V%Qmu_ref( 492 ) = 355.000000000000
  Mref_V%vph_ref( 492 ) = 12384.6000000000
  Mref_V%vsh_ref( 492 ) = 6761.73000000000
  Mref_V%eta_ref( 492 ) = 1.00000000000000

  Mref_V%radius_ref( 493 ) = 4738125.00000000
  Mref_V%density_ref( 493 ) = 4930.64000000000
  Mref_V%vpv_ref( 493 ) = 12372.5800000000
  Mref_V%vsv_ref( 493 ) = 6756.97000000000
  Mref_V%Qkappa_ref( 493 ) = 57822.5000000000
  Mref_V%Qmu_ref( 493 ) = 355.000000000000
  Mref_V%vph_ref( 493 ) = 12372.5800000000
  Mref_V%vsh_ref( 493 ) = 6756.97000000000
  Mref_V%eta_ref( 493 ) = 1.00000000000000

  Mref_V%radius_ref( 494 ) = 4747596.00000000
  Mref_V%density_ref( 494 ) = 4925.63000000000
  Mref_V%vpv_ref( 494 ) = 12360.5200000000
  Mref_V%vsv_ref( 494 ) = 6752.19000000000
  Mref_V%Qkappa_ref( 494 ) = 57822.5000000000
  Mref_V%Qmu_ref( 494 ) = 355.000000000000
  Mref_V%vph_ref( 494 ) = 12360.5200000000
  Mref_V%vsh_ref( 494 ) = 6752.19000000000
  Mref_V%eta_ref( 494 ) = 1.00000000000000

  Mref_V%radius_ref( 495 ) = 4757067.00000000
  Mref_V%density_ref( 495 ) = 4920.62000000000
  Mref_V%vpv_ref( 495 ) = 12348.4300000000
  Mref_V%vsv_ref( 495 ) = 6747.40000000000
  Mref_V%Qkappa_ref( 495 ) = 57822.5000000000
  Mref_V%Qmu_ref( 495 ) = 355.000000000000
  Mref_V%vph_ref( 495 ) = 12348.4300000000
  Mref_V%vsh_ref( 495 ) = 6747.40000000000
  Mref_V%eta_ref( 495 ) = 1.00000000000000

  Mref_V%radius_ref( 496 ) = 4766538.00000000
  Mref_V%density_ref( 496 ) = 4915.60000000000
  Mref_V%vpv_ref( 496 ) = 12336.3000000000
  Mref_V%vsv_ref( 496 ) = 6742.59000000000
  Mref_V%Qkappa_ref( 496 ) = 57822.5000000000
  Mref_V%Qmu_ref( 496 ) = 355.000000000000
  Mref_V%vph_ref( 496 ) = 12336.3000000000
  Mref_V%vsh_ref( 496 ) = 6742.59000000000
  Mref_V%eta_ref( 496 ) = 1.00000000000000

  Mref_V%radius_ref( 497 ) = 4776010.00000000
  Mref_V%density_ref( 497 ) = 4910.58000000000
  Mref_V%vpv_ref( 497 ) = 12324.1400000000
  Mref_V%vsv_ref( 497 ) = 6737.76000000000
  Mref_V%Qkappa_ref( 497 ) = 57822.5000000000
  Mref_V%Qmu_ref( 497 ) = 355.000000000000
  Mref_V%vph_ref( 497 ) = 12324.1400000000
  Mref_V%vsh_ref( 497 ) = 6737.76000000000
  Mref_V%eta_ref( 497 ) = 1.00000000000000

  Mref_V%radius_ref( 498 ) = 4785481.00000000
  Mref_V%density_ref( 498 ) = 4905.55000000000
  Mref_V%vpv_ref( 498 ) = 12311.9400000000
  Mref_V%vsv_ref( 498 ) = 6732.93000000000
  Mref_V%Qkappa_ref( 498 ) = 57822.5000000000
  Mref_V%Qmu_ref( 498 ) = 355.000000000000
  Mref_V%vph_ref( 498 ) = 12311.9400000000
  Mref_V%vsh_ref( 498 ) = 6732.93000000000
  Mref_V%eta_ref( 498 ) = 1.00000000000000

  Mref_V%radius_ref( 499 ) = 4794952.00000000
  Mref_V%density_ref( 499 ) = 4900.51000000000
  Mref_V%vpv_ref( 499 ) = 12299.7000000000
  Mref_V%vsv_ref( 499 ) = 6728.07000000000
  Mref_V%Qkappa_ref( 499 ) = 57822.5000000000
  Mref_V%Qmu_ref( 499 ) = 355.000000000000
  Mref_V%vph_ref( 499 ) = 12299.7000000000
  Mref_V%vsh_ref( 499 ) = 6728.07000000000
  Mref_V%eta_ref( 499 ) = 1.00000000000000

  Mref_V%radius_ref( 500 ) = 4804423.00000000
  Mref_V%density_ref( 500 ) = 4895.47000000000
  Mref_V%vpv_ref( 500 ) = 12287.4200000000
  Mref_V%vsv_ref( 500 ) = 6723.21000000000
  Mref_V%Qkappa_ref( 500 ) = 57822.5000000000
  Mref_V%Qmu_ref( 500 ) = 355.000000000000
  Mref_V%vph_ref( 500 ) = 12287.4200000000
  Mref_V%vsh_ref( 500 ) = 6723.21000000000
  Mref_V%eta_ref( 500 ) = 1.00000000000000

  Mref_V%radius_ref( 501 ) = 4813894.00000000
  Mref_V%density_ref( 501 ) = 4890.42000000000
  Mref_V%vpv_ref( 501 ) = 12275.1100000000
  Mref_V%vsv_ref( 501 ) = 6718.32000000000
  Mref_V%Qkappa_ref( 501 ) = 57822.5000000000
  Mref_V%Qmu_ref( 501 ) = 355.000000000000
  Mref_V%vph_ref( 501 ) = 12275.1100000000
  Mref_V%vsh_ref( 501 ) = 6718.32000000000
  Mref_V%eta_ref( 501 ) = 1.00000000000000

  Mref_V%radius_ref( 502 ) = 4823365.00000000
  Mref_V%density_ref( 502 ) = 4885.37000000000
  Mref_V%vpv_ref( 502 ) = 12262.7500000000
  Mref_V%vsv_ref( 502 ) = 6713.42000000000
  Mref_V%Qkappa_ref( 502 ) = 57822.5000000000
  Mref_V%Qmu_ref( 502 ) = 355.000000000000
  Mref_V%vph_ref( 502 ) = 12262.7500000000
  Mref_V%vsh_ref( 502 ) = 6713.42000000000
  Mref_V%eta_ref( 502 ) = 1.00000000000000

  Mref_V%radius_ref( 503 ) = 4832837.00000000
  Mref_V%density_ref( 503 ) = 4880.31000000000
  Mref_V%vpv_ref( 503 ) = 12250.3500000000
  Mref_V%vsv_ref( 503 ) = 6708.51000000000
  Mref_V%Qkappa_ref( 503 ) = 57822.5000000000
  Mref_V%Qmu_ref( 503 ) = 355.000000000000
  Mref_V%vph_ref( 503 ) = 12250.3500000000
  Mref_V%vsh_ref( 503 ) = 6708.51000000000
  Mref_V%eta_ref( 503 ) = 1.00000000000000

  Mref_V%radius_ref( 504 ) = 4842308.00000000
  Mref_V%density_ref( 504 ) = 4875.24000000000
  Mref_V%vpv_ref( 504 ) = 12237.9200000000
  Mref_V%vsv_ref( 504 ) = 6703.57000000000
  Mref_V%Qkappa_ref( 504 ) = 57822.5000000000
  Mref_V%Qmu_ref( 504 ) = 355.000000000000
  Mref_V%vph_ref( 504 ) = 12237.9200000000
  Mref_V%vsh_ref( 504 ) = 6703.57000000000
  Mref_V%eta_ref( 504 ) = 1.00000000000000

  Mref_V%radius_ref( 505 ) = 4851779.00000000
  Mref_V%density_ref( 505 ) = 4870.16000000000
  Mref_V%vpv_ref( 505 ) = 12225.4400000000
  Mref_V%vsv_ref( 505 ) = 6698.62000000000
  Mref_V%Qkappa_ref( 505 ) = 57822.5000000000
  Mref_V%Qmu_ref( 505 ) = 355.000000000000
  Mref_V%vph_ref( 505 ) = 12225.4400000000
  Mref_V%vsh_ref( 505 ) = 6698.62000000000
  Mref_V%eta_ref( 505 ) = 1.00000000000000

  Mref_V%radius_ref( 506 ) = 4861250.00000000
  Mref_V%density_ref( 506 ) = 4865.08000000000
  Mref_V%vpv_ref( 506 ) = 12212.9200000000
  Mref_V%vsv_ref( 506 ) = 6693.66000000000
  Mref_V%Qkappa_ref( 506 ) = 57822.5000000000
  Mref_V%Qmu_ref( 506 ) = 355.000000000000
  Mref_V%vph_ref( 506 ) = 12212.9200000000
  Mref_V%vsh_ref( 506 ) = 6693.66000000000
  Mref_V%eta_ref( 506 ) = 1.00000000000000

  Mref_V%radius_ref( 507 ) = 4870721.00000000
  Mref_V%density_ref( 507 ) = 4859.99000000000
  Mref_V%vpv_ref( 507 ) = 12200.3600000000
  Mref_V%vsv_ref( 507 ) = 6688.68000000000
  Mref_V%Qkappa_ref( 507 ) = 57822.5000000000
  Mref_V%Qmu_ref( 507 ) = 355.000000000000
  Mref_V%vph_ref( 507 ) = 12200.3600000000
  Mref_V%vsh_ref( 507 ) = 6688.68000000000
  Mref_V%eta_ref( 507 ) = 1.00000000000000

  Mref_V%radius_ref( 508 ) = 4880192.00000000
  Mref_V%density_ref( 508 ) = 4854.90000000000
  Mref_V%vpv_ref( 508 ) = 12187.7600000000
  Mref_V%vsv_ref( 508 ) = 6683.68000000000
  Mref_V%Qkappa_ref( 508 ) = 57822.5000000000
  Mref_V%Qmu_ref( 508 ) = 355.000000000000
  Mref_V%vph_ref( 508 ) = 12187.7600000000
  Mref_V%vsh_ref( 508 ) = 6683.68000000000
  Mref_V%eta_ref( 508 ) = 1.00000000000000

  Mref_V%radius_ref( 509 ) = 4889663.00000000
  Mref_V%density_ref( 509 ) = 4849.80000000000
  Mref_V%vpv_ref( 509 ) = 12175.1100000000
  Mref_V%vsv_ref( 509 ) = 6678.66000000000
  Mref_V%Qkappa_ref( 509 ) = 57822.5000000000
  Mref_V%Qmu_ref( 509 ) = 355.000000000000
  Mref_V%vph_ref( 509 ) = 12175.1100000000
  Mref_V%vsh_ref( 509 ) = 6678.66000000000
  Mref_V%eta_ref( 509 ) = 1.00000000000000

  Mref_V%radius_ref( 510 ) = 4899135.00000000
  Mref_V%density_ref( 510 ) = 4844.69000000000
  Mref_V%vpv_ref( 510 ) = 12162.4300000000
  Mref_V%vsv_ref( 510 ) = 6673.63000000000
  Mref_V%Qkappa_ref( 510 ) = 57822.5000000000
  Mref_V%Qmu_ref( 510 ) = 355.000000000000
  Mref_V%vph_ref( 510 ) = 12162.4300000000
  Mref_V%vsh_ref( 510 ) = 6673.63000000000
  Mref_V%eta_ref( 510 ) = 1.00000000000000

  Mref_V%radius_ref( 511 ) = 4908606.00000000
  Mref_V%density_ref( 511 ) = 4839.57000000000
  Mref_V%vpv_ref( 511 ) = 12149.6900000000
  Mref_V%vsv_ref( 511 ) = 6668.58000000000
  Mref_V%Qkappa_ref( 511 ) = 57822.5000000000
  Mref_V%Qmu_ref( 511 ) = 355.000000000000
  Mref_V%vph_ref( 511 ) = 12149.6900000000
  Mref_V%vsh_ref( 511 ) = 6668.58000000000
  Mref_V%eta_ref( 511 ) = 1.00000000000000

  Mref_V%radius_ref( 512 ) = 4918077.00000000
  Mref_V%density_ref( 512 ) = 4834.45000000000
  Mref_V%vpv_ref( 512 ) = 12136.9100000000
  Mref_V%vsv_ref( 512 ) = 6663.51000000000
  Mref_V%Qkappa_ref( 512 ) = 57822.5000000000
  Mref_V%Qmu_ref( 512 ) = 355.000000000000
  Mref_V%vph_ref( 512 ) = 12136.9100000000
  Mref_V%vsh_ref( 512 ) = 6663.51000000000
  Mref_V%eta_ref( 512 ) = 1.00000000000000

  Mref_V%radius_ref( 513 ) = 4927548.00000000
  Mref_V%density_ref( 513 ) = 4829.32000000000
  Mref_V%vpv_ref( 513 ) = 12124.0900000000
  Mref_V%vsv_ref( 513 ) = 6658.43000000000
  Mref_V%Qkappa_ref( 513 ) = 57822.5000000000
  Mref_V%Qmu_ref( 513 ) = 355.000000000000
  Mref_V%vph_ref( 513 ) = 12124.0900000000
  Mref_V%vsh_ref( 513 ) = 6658.43000000000
  Mref_V%eta_ref( 513 ) = 1.00000000000000

  Mref_V%radius_ref( 514 ) = 4937019.00000000
  Mref_V%density_ref( 514 ) = 4824.18000000000
  Mref_V%vpv_ref( 514 ) = 12111.2200000000
  Mref_V%vsv_ref( 514 ) = 6653.32000000000
  Mref_V%Qkappa_ref( 514 ) = 57822.5000000000
  Mref_V%Qmu_ref( 514 ) = 355.000000000000
  Mref_V%vph_ref( 514 ) = 12111.2200000000
  Mref_V%vsh_ref( 514 ) = 6653.32000000000
  Mref_V%eta_ref( 514 ) = 1.00000000000000

  Mref_V%radius_ref( 515 ) = 4946490.00000000
  Mref_V%density_ref( 515 ) = 4819.04000000000
  Mref_V%vpv_ref( 515 ) = 12098.3100000000
  Mref_V%vsv_ref( 515 ) = 6648.20000000000
  Mref_V%Qkappa_ref( 515 ) = 57822.5000000000
  Mref_V%Qmu_ref( 515 ) = 355.000000000000
  Mref_V%vph_ref( 515 ) = 12098.3100000000
  Mref_V%vsh_ref( 515 ) = 6648.20000000000
  Mref_V%eta_ref( 515 ) = 1.00000000000000

  Mref_V%radius_ref( 516 ) = 4955962.00000000
  Mref_V%density_ref( 516 ) = 4813.88000000000
  Mref_V%vpv_ref( 516 ) = 12085.3400000000
  Mref_V%vsv_ref( 516 ) = 6643.06000000000
  Mref_V%Qkappa_ref( 516 ) = 57822.5000000000
  Mref_V%Qmu_ref( 516 ) = 355.000000000000
  Mref_V%vph_ref( 516 ) = 12085.3400000000
  Mref_V%vsh_ref( 516 ) = 6643.06000000000
  Mref_V%eta_ref( 516 ) = 1.00000000000000

  Mref_V%radius_ref( 517 ) = 4965433.00000000
  Mref_V%density_ref( 517 ) = 4808.73000000000
  Mref_V%vpv_ref( 517 ) = 12072.3400000000
  Mref_V%vsv_ref( 517 ) = 6637.90000000000
  Mref_V%Qkappa_ref( 517 ) = 57822.5000000000
  Mref_V%Qmu_ref( 517 ) = 355.000000000000
  Mref_V%vph_ref( 517 ) = 12072.3400000000
  Mref_V%vsh_ref( 517 ) = 6637.90000000000
  Mref_V%eta_ref( 517 ) = 1.00000000000000

  Mref_V%radius_ref( 518 ) = 4974904.00000000
  Mref_V%density_ref( 518 ) = 4803.56000000000
  Mref_V%vpv_ref( 518 ) = 12059.2800000000
  Mref_V%vsv_ref( 518 ) = 6632.73000000000
  Mref_V%Qkappa_ref( 518 ) = 57822.5000000000
  Mref_V%Qmu_ref( 518 ) = 355.000000000000
  Mref_V%vph_ref( 518 ) = 12059.2800000000
  Mref_V%vsh_ref( 518 ) = 6632.73000000000
  Mref_V%eta_ref( 518 ) = 1.00000000000000

  Mref_V%radius_ref( 519 ) = 4984375.00000000
  Mref_V%density_ref( 519 ) = 4798.39000000000
  Mref_V%vpv_ref( 519 ) = 12046.1700000000
  Mref_V%vsv_ref( 519 ) = 6627.53000000000
  Mref_V%Qkappa_ref( 519 ) = 57822.5000000000
  Mref_V%Qmu_ref( 519 ) = 355.000000000000
  Mref_V%vph_ref( 519 ) = 12046.1700000000
  Mref_V%vsh_ref( 519 ) = 6627.53000000000
  Mref_V%eta_ref( 519 ) = 1.00000000000000

  Mref_V%radius_ref( 520 ) = 4993846.00000000
  Mref_V%density_ref( 520 ) = 4793.20000000000
  Mref_V%vpv_ref( 520 ) = 12033.0200000000
  Mref_V%vsv_ref( 520 ) = 6622.31000000000
  Mref_V%Qkappa_ref( 520 ) = 57822.5000000000
  Mref_V%Qmu_ref( 520 ) = 355.000000000000
  Mref_V%vph_ref( 520 ) = 12033.0200000000
  Mref_V%vsh_ref( 520 ) = 6622.31000000000
  Mref_V%eta_ref( 520 ) = 1.00000000000000

  Mref_V%radius_ref( 521 ) = 5003317.00000000
  Mref_V%density_ref( 521 ) = 4788.02000000000
  Mref_V%vpv_ref( 521 ) = 12019.8200000000
  Mref_V%vsv_ref( 521 ) = 6617.08000000000
  Mref_V%Qkappa_ref( 521 ) = 57822.5000000000
  Mref_V%Qmu_ref( 521 ) = 355.000000000000
  Mref_V%vph_ref( 521 ) = 12019.8200000000
  Mref_V%vsh_ref( 521 ) = 6617.08000000000
  Mref_V%eta_ref( 521 ) = 1.00000000000000

  Mref_V%radius_ref( 522 ) = 5012788.00000000
  Mref_V%density_ref( 522 ) = 4782.82000000000
  Mref_V%vpv_ref( 522 ) = 12006.5600000000
  Mref_V%vsv_ref( 522 ) = 6611.82000000000
  Mref_V%Qkappa_ref( 522 ) = 57822.5000000000
  Mref_V%Qmu_ref( 522 ) = 355.000000000000
  Mref_V%vph_ref( 522 ) = 12006.5600000000
  Mref_V%vsh_ref( 522 ) = 6611.82000000000
  Mref_V%eta_ref( 522 ) = 1.00000000000000

  Mref_V%radius_ref( 523 ) = 5022260.00000000
  Mref_V%density_ref( 523 ) = 4777.61000000000
  Mref_V%vpv_ref( 523 ) = 11993.2600000000
  Mref_V%vsv_ref( 523 ) = 6606.55000000000
  Mref_V%Qkappa_ref( 523 ) = 57822.5000000000
  Mref_V%Qmu_ref( 523 ) = 355.000000000000
  Mref_V%vph_ref( 523 ) = 11993.2600000000
  Mref_V%vsh_ref( 523 ) = 6606.55000000000
  Mref_V%eta_ref( 523 ) = 1.00000000000000

  Mref_V%radius_ref( 524 ) = 5031731.00000000
  Mref_V%density_ref( 524 ) = 4772.40000000000
  Mref_V%vpv_ref( 524 ) = 11979.9000000000
  Mref_V%vsv_ref( 524 ) = 6601.26000000000
  Mref_V%Qkappa_ref( 524 ) = 57822.5000000000
  Mref_V%Qmu_ref( 524 ) = 355.000000000000
  Mref_V%vph_ref( 524 ) = 11979.9000000000
  Mref_V%vsh_ref( 524 ) = 6601.26000000000
  Mref_V%eta_ref( 524 ) = 1.00000000000000

  Mref_V%radius_ref( 525 ) = 5041202.00000000
  Mref_V%density_ref( 525 ) = 4767.18000000000
  Mref_V%vpv_ref( 525 ) = 11966.5000000000
  Mref_V%vsv_ref( 525 ) = 6595.94000000000
  Mref_V%Qkappa_ref( 525 ) = 57822.5000000000
  Mref_V%Qmu_ref( 525 ) = 355.000000000000
  Mref_V%vph_ref( 525 ) = 11966.5000000000
  Mref_V%vsh_ref( 525 ) = 6595.94000000000
  Mref_V%eta_ref( 525 ) = 1.00000000000000

  Mref_V%radius_ref( 526 ) = 5050673.00000000
  Mref_V%density_ref( 526 ) = 4761.95000000000
  Mref_V%vpv_ref( 526 ) = 11953.0400000000
  Mref_V%vsv_ref( 526 ) = 6590.61000000000
  Mref_V%Qkappa_ref( 526 ) = 57822.5000000000
  Mref_V%Qmu_ref( 526 ) = 355.000000000000
  Mref_V%vph_ref( 526 ) = 11953.0400000000
  Mref_V%vsh_ref( 526 ) = 6590.61000000000
  Mref_V%eta_ref( 526 ) = 1.00000000000000

  Mref_V%radius_ref( 527 ) = 5060144.00000000
  Mref_V%density_ref( 527 ) = 4756.72000000000
  Mref_V%vpv_ref( 527 ) = 11939.5300000000
  Mref_V%vsv_ref( 527 ) = 6584.91000000000
  Mref_V%Qkappa_ref( 527 ) = 57822.5000000000
  Mref_V%Qmu_ref( 527 ) = 355.000000000000
  Mref_V%vph_ref( 527 ) = 11939.5300000000
  Mref_V%vsh_ref( 527 ) = 6584.91000000000
  Mref_V%eta_ref( 527 ) = 1.00000000000000

  Mref_V%radius_ref( 528 ) = 5069615.00000000
  Mref_V%density_ref( 528 ) = 4751.47000000000
  Mref_V%vpv_ref( 528 ) = 11925.9700000000
  Mref_V%vsv_ref( 528 ) = 6579.51000000000
  Mref_V%Qkappa_ref( 528 ) = 57822.5000000000
  Mref_V%Qmu_ref( 528 ) = 355.000000000000
  Mref_V%vph_ref( 528 ) = 11925.9700000000
  Mref_V%vsh_ref( 528 ) = 6579.51000000000
  Mref_V%eta_ref( 528 ) = 1.00000000000000

  Mref_V%radius_ref( 529 ) = 5079087.00000000
  Mref_V%density_ref( 529 ) = 4746.22000000000
  Mref_V%vpv_ref( 529 ) = 11912.3500000000
  Mref_V%vsv_ref( 529 ) = 6574.11000000000
  Mref_V%Qkappa_ref( 529 ) = 57822.5000000000
  Mref_V%Qmu_ref( 529 ) = 355.000000000000
  Mref_V%vph_ref( 529 ) = 11912.3500000000
  Mref_V%vsh_ref( 529 ) = 6574.11000000000
  Mref_V%eta_ref( 529 ) = 1.00000000000000

  Mref_V%radius_ref( 530 ) = 5088558.00000000
  Mref_V%density_ref( 530 ) = 4740.95000000000
  Mref_V%vpv_ref( 530 ) = 11898.6900000000
  Mref_V%vsv_ref( 530 ) = 6568.67000000000
  Mref_V%Qkappa_ref( 530 ) = 57822.5000000000
  Mref_V%Qmu_ref( 530 ) = 355.000000000000
  Mref_V%vph_ref( 530 ) = 11898.6900000000
  Mref_V%vsh_ref( 530 ) = 6568.67000000000
  Mref_V%eta_ref( 530 ) = 1.00000000000000

  Mref_V%radius_ref( 531 ) = 5098029.00000000
  Mref_V%density_ref( 531 ) = 4735.68000000000
  Mref_V%vpv_ref( 531 ) = 11884.9600000000
  Mref_V%vsv_ref( 531 ) = 6563.22000000000
  Mref_V%Qkappa_ref( 531 ) = 57822.5000000000
  Mref_V%Qmu_ref( 531 ) = 355.000000000000
  Mref_V%vph_ref( 531 ) = 11884.9600000000
  Mref_V%vsh_ref( 531 ) = 6563.22000000000
  Mref_V%eta_ref( 531 ) = 1.00000000000000

  Mref_V%radius_ref( 532 ) = 5107500.00000000
  Mref_V%density_ref( 532 ) = 4730.40000000000
  Mref_V%vpv_ref( 532 ) = 11871.1900000000
  Mref_V%vsv_ref( 532 ) = 6557.74000000000
  Mref_V%Qkappa_ref( 532 ) = 57822.5000000000
  Mref_V%Qmu_ref( 532 ) = 355.000000000000
  Mref_V%vph_ref( 532 ) = 11871.1900000000
  Mref_V%vsh_ref( 532 ) = 6557.74000000000
  Mref_V%eta_ref( 532 ) = 1.00000000000000

  Mref_V%radius_ref( 533 ) = 5116971.00000000
  Mref_V%density_ref( 533 ) = 4725.10000000000
  Mref_V%vpv_ref( 533 ) = 11857.3700000000
  Mref_V%vsv_ref( 533 ) = 6552.24000000000
  Mref_V%Qkappa_ref( 533 ) = 57822.5000000000
  Mref_V%Qmu_ref( 533 ) = 355.000000000000
  Mref_V%vph_ref( 533 ) = 11857.3700000000
  Mref_V%vsh_ref( 533 ) = 6552.24000000000
  Mref_V%eta_ref( 533 ) = 1.00000000000000

  Mref_V%radius_ref( 534 ) = 5126442.00000000
  Mref_V%density_ref( 534 ) = 4719.80000000000
  Mref_V%vpv_ref( 534 ) = 11843.4800000000
  Mref_V%vsv_ref( 534 ) = 6546.73000000000
  Mref_V%Qkappa_ref( 534 ) = 57822.5000000000
  Mref_V%Qmu_ref( 534 ) = 355.000000000000
  Mref_V%vph_ref( 534 ) = 11843.4800000000
  Mref_V%vsh_ref( 534 ) = 6546.73000000000
  Mref_V%eta_ref( 534 ) = 1.00000000000000

  Mref_V%radius_ref( 535 ) = 5135913.00000000
  Mref_V%density_ref( 535 ) = 4714.48000000000
  Mref_V%vpv_ref( 535 ) = 11829.5500000000
  Mref_V%vsv_ref( 535 ) = 6541.19000000000
  Mref_V%Qkappa_ref( 535 ) = 57822.5000000000
  Mref_V%Qmu_ref( 535 ) = 355.000000000000
  Mref_V%vph_ref( 535 ) = 11829.5500000000
  Mref_V%vsh_ref( 535 ) = 6541.19000000000
  Mref_V%eta_ref( 535 ) = 1.00000000000000

  Mref_V%radius_ref( 536 ) = 5145385.00000000
  Mref_V%density_ref( 536 ) = 4709.15000000000
  Mref_V%vpv_ref( 536 ) = 11815.5700000000
  Mref_V%vsv_ref( 536 ) = 6535.63000000000
  Mref_V%Qkappa_ref( 536 ) = 57822.5000000000
  Mref_V%Qmu_ref( 536 ) = 355.000000000000
  Mref_V%vph_ref( 536 ) = 11815.5700000000
  Mref_V%vsh_ref( 536 ) = 6535.63000000000
  Mref_V%eta_ref( 536 ) = 1.00000000000000

  Mref_V%radius_ref( 537 ) = 5154856.00000000
  Mref_V%density_ref( 537 ) = 4703.81000000000
  Mref_V%vpv_ref( 537 ) = 11801.5300000000
  Mref_V%vsv_ref( 537 ) = 6530.05000000000
  Mref_V%Qkappa_ref( 537 ) = 57822.5000000000
  Mref_V%Qmu_ref( 537 ) = 355.000000000000
  Mref_V%vph_ref( 537 ) = 11801.5300000000
  Mref_V%vsh_ref( 537 ) = 6530.05000000000
  Mref_V%eta_ref( 537 ) = 1.00000000000000

  Mref_V%radius_ref( 538 ) = 5164327.00000000
  Mref_V%density_ref( 538 ) = 4698.44000000000
  Mref_V%vpv_ref( 538 ) = 11787.4400000000
  Mref_V%vsv_ref( 538 ) = 6524.44000000000
  Mref_V%Qkappa_ref( 538 ) = 57822.5000000000
  Mref_V%Qmu_ref( 538 ) = 355.000000000000
  Mref_V%vph_ref( 538 ) = 11787.4400000000
  Mref_V%vsh_ref( 538 ) = 6524.44000000000
  Mref_V%eta_ref( 538 ) = 1.00000000000000

  Mref_V%radius_ref( 539 ) = 5173798.00000000
  Mref_V%density_ref( 539 ) = 4693.08000000000
  Mref_V%vpv_ref( 539 ) = 11773.3000000000
  Mref_V%vsv_ref( 539 ) = 6518.82000000000
  Mref_V%Qkappa_ref( 539 ) = 57822.5000000000
  Mref_V%Qmu_ref( 539 ) = 355.000000000000
  Mref_V%vph_ref( 539 ) = 11773.3000000000
  Mref_V%vsh_ref( 539 ) = 6518.82000000000
  Mref_V%eta_ref( 539 ) = 1.00000000000000

  Mref_V%radius_ref( 540 ) = 5183269.00000000
  Mref_V%density_ref( 540 ) = 4687.69000000000
  Mref_V%vpv_ref( 540 ) = 11759.1000000000
  Mref_V%vsv_ref( 540 ) = 6513.17000000000
  Mref_V%Qkappa_ref( 540 ) = 57822.5000000000
  Mref_V%Qmu_ref( 540 ) = 355.000000000000
  Mref_V%vph_ref( 540 ) = 11759.1000000000
  Mref_V%vsh_ref( 540 ) = 6513.17000000000
  Mref_V%eta_ref( 540 ) = 1.00000000000000

  Mref_V%radius_ref( 541 ) = 5192740.00000000
  Mref_V%density_ref( 541 ) = 4682.29000000000
  Mref_V%vpv_ref( 541 ) = 11744.8500000000
  Mref_V%vsv_ref( 541 ) = 6507.50000000000
  Mref_V%Qkappa_ref( 541 ) = 57822.5000000000
  Mref_V%Qmu_ref( 541 ) = 355.000000000000
  Mref_V%vph_ref( 541 ) = 11744.8500000000
  Mref_V%vsh_ref( 541 ) = 6507.50000000000
  Mref_V%eta_ref( 541 ) = 1.00000000000000

  Mref_V%radius_ref( 542 ) = 5202212.00000000
  Mref_V%density_ref( 542 ) = 4676.87000000000
  Mref_V%vpv_ref( 542 ) = 11730.5500000000
  Mref_V%vsv_ref( 542 ) = 6501.80000000000
  Mref_V%Qkappa_ref( 542 ) = 57822.5000000000
  Mref_V%Qmu_ref( 542 ) = 355.000000000000
  Mref_V%vph_ref( 542 ) = 11730.5500000000
  Mref_V%vsh_ref( 542 ) = 6501.80000000000
  Mref_V%eta_ref( 542 ) = 1.00000000000000

  Mref_V%radius_ref( 543 ) = 5211683.00000000
  Mref_V%density_ref( 543 ) = 4671.44000000000
  Mref_V%vpv_ref( 543 ) = 11716.1800000000
  Mref_V%vsv_ref( 543 ) = 6496.09000000000
  Mref_V%Qkappa_ref( 543 ) = 57822.5000000000
  Mref_V%Qmu_ref( 543 ) = 355.000000000000
  Mref_V%vph_ref( 543 ) = 11716.1800000000
  Mref_V%vsh_ref( 543 ) = 6496.09000000000
  Mref_V%eta_ref( 543 ) = 1.00000000000000

  Mref_V%radius_ref( 544 ) = 5221154.00000000
  Mref_V%density_ref( 544 ) = 4665.99000000000
  Mref_V%vpv_ref( 544 ) = 11701.7800000000
  Mref_V%vsv_ref( 544 ) = 6490.35000000000
  Mref_V%Qkappa_ref( 544 ) = 57822.5000000000
  Mref_V%Qmu_ref( 544 ) = 355.000000000000
  Mref_V%vph_ref( 544 ) = 11701.7800000000
  Mref_V%vsh_ref( 544 ) = 6490.35000000000
  Mref_V%eta_ref( 544 ) = 1.00000000000000

  Mref_V%radius_ref( 545 ) = 5230625.00000000
  Mref_V%density_ref( 545 ) = 4660.52000000000
  Mref_V%vpv_ref( 545 ) = 11687.3100000000
  Mref_V%vsv_ref( 545 ) = 6484.59000000000
  Mref_V%Qkappa_ref( 545 ) = 57822.5000000000
  Mref_V%Qmu_ref( 545 ) = 355.000000000000
  Mref_V%vph_ref( 545 ) = 11687.3100000000
  Mref_V%vsh_ref( 545 ) = 6484.59000000000
  Mref_V%eta_ref( 545 ) = 1.00000000000000

  Mref_V%radius_ref( 546 ) = 5240096.00000000
  Mref_V%density_ref( 546 ) = 4655.03000000000
  Mref_V%vpv_ref( 546 ) = 11672.8000000000
  Mref_V%vsv_ref( 546 ) = 6478.80000000000
  Mref_V%Qkappa_ref( 546 ) = 57822.5000000000
  Mref_V%Qmu_ref( 546 ) = 355.000000000000
  Mref_V%vph_ref( 546 ) = 11672.8000000000
  Mref_V%vsh_ref( 546 ) = 6478.80000000000
  Mref_V%eta_ref( 546 ) = 1.00000000000000

  Mref_V%radius_ref( 547 ) = 5249567.00000000
  Mref_V%density_ref( 547 ) = 4649.52000000000
  Mref_V%vpv_ref( 547 ) = 11658.2300000000
  Mref_V%vsv_ref( 547 ) = 6472.99000000000
  Mref_V%Qkappa_ref( 547 ) = 57822.5000000000
  Mref_V%Qmu_ref( 547 ) = 355.000000000000
  Mref_V%vph_ref( 547 ) = 11658.2300000000
  Mref_V%vsh_ref( 547 ) = 6472.99000000000
  Mref_V%eta_ref( 547 ) = 1.00000000000000

  Mref_V%radius_ref( 548 ) = 5259038.00000000
  Mref_V%density_ref( 548 ) = 4644.00000000000
  Mref_V%vpv_ref( 548 ) = 11643.6000000000
  Mref_V%vsv_ref( 548 ) = 6467.16000000000
  Mref_V%Qkappa_ref( 548 ) = 57822.5000000000
  Mref_V%Qmu_ref( 548 ) = 355.000000000000
  Mref_V%vph_ref( 548 ) = 11643.6000000000
  Mref_V%vsh_ref( 548 ) = 6467.16000000000
  Mref_V%eta_ref( 548 ) = 1.00000000000000

  Mref_V%radius_ref( 549 ) = 5268510.00000000
  Mref_V%density_ref( 549 ) = 4638.46000000000
  Mref_V%vpv_ref( 549 ) = 11628.9200000000
  Mref_V%vsv_ref( 549 ) = 6461.30000000000
  Mref_V%Qkappa_ref( 549 ) = 57822.5000000000
  Mref_V%Qmu_ref( 549 ) = 355.000000000000
  Mref_V%vph_ref( 549 ) = 11628.9200000000
  Mref_V%vsh_ref( 549 ) = 6461.30000000000
  Mref_V%eta_ref( 549 ) = 1.00000000000000

  Mref_V%radius_ref( 550 ) = 5277981.00000000
  Mref_V%density_ref( 550 ) = 4632.89000000000
  Mref_V%vpv_ref( 550 ) = 11614.1900000000
  Mref_V%vsv_ref( 550 ) = 6455.42000000000
  Mref_V%Qkappa_ref( 550 ) = 57822.5000000000
  Mref_V%Qmu_ref( 550 ) = 355.000000000000
  Mref_V%vph_ref( 550 ) = 11614.1900000000
  Mref_V%vsh_ref( 550 ) = 6455.42000000000
  Mref_V%eta_ref( 550 ) = 1.00000000000000

  Mref_V%radius_ref( 551 ) = 5287452.00000000
  Mref_V%density_ref( 551 ) = 4627.31000000000
  Mref_V%vpv_ref( 551 ) = 11599.4000000000
  Mref_V%vsv_ref( 551 ) = 6449.51000000000
  Mref_V%Qkappa_ref( 551 ) = 57822.5000000000
  Mref_V%Qmu_ref( 551 ) = 355.000000000000
  Mref_V%vph_ref( 551 ) = 11599.4000000000
  Mref_V%vsh_ref( 551 ) = 6449.51000000000
  Mref_V%eta_ref( 551 ) = 1.00000000000000

  Mref_V%radius_ref( 552 ) = 5296923.00000000
  Mref_V%density_ref( 552 ) = 4621.70000000000
  Mref_V%vpv_ref( 552 ) = 11584.5700000000
  Mref_V%vsv_ref( 552 ) = 6443.58000000000
  Mref_V%Qkappa_ref( 552 ) = 57822.5000000000
  Mref_V%Qmu_ref( 552 ) = 355.000000000000
  Mref_V%vph_ref( 552 ) = 11584.5700000000
  Mref_V%vsh_ref( 552 ) = 6443.58000000000
  Mref_V%eta_ref( 552 ) = 1.00000000000000

  Mref_V%radius_ref( 553 ) = 5306394.00000000
  Mref_V%density_ref( 553 ) = 4616.08000000000
  Mref_V%vpv_ref( 553 ) = 11569.6800000000
  Mref_V%vsv_ref( 553 ) = 6437.63000000000
  Mref_V%Qkappa_ref( 553 ) = 57822.5000000000
  Mref_V%Qmu_ref( 553 ) = 355.000000000000
  Mref_V%vph_ref( 553 ) = 11569.6800000000
  Mref_V%vsh_ref( 553 ) = 6437.63000000000
  Mref_V%eta_ref( 553 ) = 1.00000000000000

  Mref_V%radius_ref( 554 ) = 5315865.00000000
  Mref_V%density_ref( 554 ) = 4610.44000000000
  Mref_V%vpv_ref( 554 ) = 11554.7200000000
  Mref_V%vsv_ref( 554 ) = 6431.65000000000
  Mref_V%Qkappa_ref( 554 ) = 57822.5000000000
  Mref_V%Qmu_ref( 554 ) = 355.000000000000
  Mref_V%vph_ref( 554 ) = 11554.7200000000
  Mref_V%vsh_ref( 554 ) = 6431.65000000000
  Mref_V%eta_ref( 554 ) = 1.00000000000000

  Mref_V%radius_ref( 555 ) = 5325337.00000000
  Mref_V%density_ref( 555 ) = 4604.76000000000
  Mref_V%vpv_ref( 555 ) = 11539.7200000000
  Mref_V%vsv_ref( 555 ) = 6425.65000000000
  Mref_V%Qkappa_ref( 555 ) = 57822.5000000000
  Mref_V%Qmu_ref( 555 ) = 355.000000000000
  Mref_V%vph_ref( 555 ) = 11539.7200000000
  Mref_V%vsh_ref( 555 ) = 6425.65000000000
  Mref_V%eta_ref( 555 ) = 1.00000000000000

  Mref_V%radius_ref( 556 ) = 5334808.00000000
  Mref_V%density_ref( 556 ) = 4599.08000000000
  Mref_V%vpv_ref( 556 ) = 11524.6700000000
  Mref_V%vsv_ref( 556 ) = 6419.61000000000
  Mref_V%Qkappa_ref( 556 ) = 57822.5000000000
  Mref_V%Qmu_ref( 556 ) = 355.000000000000
  Mref_V%vph_ref( 556 ) = 11524.6700000000
  Mref_V%vsh_ref( 556 ) = 6419.61000000000
  Mref_V%eta_ref( 556 ) = 1.00000000000000

  Mref_V%radius_ref( 557 ) = 5344279.00000000
  Mref_V%density_ref( 557 ) = 4593.36000000000
  Mref_V%vpv_ref( 557 ) = 11509.5600000000
  Mref_V%vsv_ref( 557 ) = 6413.56000000000
  Mref_V%Qkappa_ref( 557 ) = 57822.5000000000
  Mref_V%Qmu_ref( 557 ) = 355.000000000000
  Mref_V%vph_ref( 557 ) = 11509.5600000000
  Mref_V%vsh_ref( 557 ) = 6413.56000000000
  Mref_V%eta_ref( 557 ) = 1.00000000000000

  Mref_V%radius_ref( 558 ) = 5353750.00000000
  Mref_V%density_ref( 558 ) = 4587.63000000000
  Mref_V%vpv_ref( 558 ) = 11494.3900000000
  Mref_V%vsv_ref( 558 ) = 6407.48000000000
  Mref_V%Qkappa_ref( 558 ) = 57822.5000000000
  Mref_V%Qmu_ref( 558 ) = 355.000000000000
  Mref_V%vph_ref( 558 ) = 11494.3900000000
  Mref_V%vsh_ref( 558 ) = 6407.48000000000
  Mref_V%eta_ref( 558 ) = 1.00000000000000

  Mref_V%radius_ref( 559 ) = 5363221.00000000
  Mref_V%density_ref( 559 ) = 4581.86000000000
  Mref_V%vpv_ref( 559 ) = 11479.1700000000
  Mref_V%vsv_ref( 559 ) = 6401.37000000000
  Mref_V%Qkappa_ref( 559 ) = 57822.5000000000
  Mref_V%Qmu_ref( 559 ) = 355.000000000000
  Mref_V%vph_ref( 559 ) = 11479.1700000000
  Mref_V%vsh_ref( 559 ) = 6401.37000000000
  Mref_V%eta_ref( 559 ) = 1.00000000000000

  Mref_V%radius_ref( 560 ) = 5372692.00000000
  Mref_V%density_ref( 560 ) = 4576.07000000000
  Mref_V%vpv_ref( 560 ) = 11463.8900000000
  Mref_V%vsv_ref( 560 ) = 6395.25000000000
  Mref_V%Qkappa_ref( 560 ) = 57822.5000000000
  Mref_V%Qmu_ref( 560 ) = 355.000000000000
  Mref_V%vph_ref( 560 ) = 11463.8900000000
  Mref_V%vsh_ref( 560 ) = 6395.25000000000
  Mref_V%eta_ref( 560 ) = 1.00000000000000

  Mref_V%radius_ref( 561 ) = 5382163.00000000
  Mref_V%density_ref( 561 ) = 4570.26000000000
  Mref_V%vpv_ref( 561 ) = 11448.5500000000
  Mref_V%vsv_ref( 561 ) = 6389.09000000000
  Mref_V%Qkappa_ref( 561 ) = 57822.5000000000
  Mref_V%Qmu_ref( 561 ) = 355.000000000000
  Mref_V%vph_ref( 561 ) = 11448.5500000000
  Mref_V%vsh_ref( 561 ) = 6389.09000000000
  Mref_V%eta_ref( 561 ) = 1.00000000000000

  Mref_V%radius_ref( 562 ) = 5391635.00000000
  Mref_V%density_ref( 562 ) = 4564.43000000000
  Mref_V%vpv_ref( 562 ) = 11433.1700000000
  Mref_V%vsv_ref( 562 ) = 6382.91000000000
  Mref_V%Qkappa_ref( 562 ) = 57822.5000000000
  Mref_V%Qmu_ref( 562 ) = 355.000000000000
  Mref_V%vph_ref( 562 ) = 11433.1700000000
  Mref_V%vsh_ref( 562 ) = 6382.91000000000
  Mref_V%eta_ref( 562 ) = 1.00000000000000

  Mref_V%radius_ref( 563 ) = 5401106.00000000
  Mref_V%density_ref( 563 ) = 4558.56000000000
  Mref_V%vpv_ref( 563 ) = 11417.7300000000
  Mref_V%vsv_ref( 563 ) = 6376.70000000000
  Mref_V%Qkappa_ref( 563 ) = 57822.5000000000
  Mref_V%Qmu_ref( 563 ) = 355.000000000000
  Mref_V%vph_ref( 563 ) = 11417.7300000000
  Mref_V%vsh_ref( 563 ) = 6376.70000000000
  Mref_V%eta_ref( 563 ) = 1.00000000000000

  Mref_V%radius_ref( 564 ) = 5410577.00000000
  Mref_V%density_ref( 564 ) = 4552.67000000000
  Mref_V%vpv_ref( 564 ) = 11402.2300000000
  Mref_V%vsv_ref( 564 ) = 6370.46000000000
  Mref_V%Qkappa_ref( 564 ) = 57822.5000000000
  Mref_V%Qmu_ref( 564 ) = 355.000000000000
  Mref_V%vph_ref( 564 ) = 11402.2300000000
  Mref_V%vsh_ref( 564 ) = 6370.46000000000
  Mref_V%eta_ref( 564 ) = 1.00000000000000

  Mref_V%radius_ref( 565 ) = 5420048.00000000
  Mref_V%density_ref( 565 ) = 4546.76000000000
  Mref_V%vpv_ref( 565 ) = 11386.6800000000
  Mref_V%vsv_ref( 565 ) = 6364.20000000000
  Mref_V%Qkappa_ref( 565 ) = 57822.5000000000
  Mref_V%Qmu_ref( 565 ) = 355.000000000000
  Mref_V%vph_ref( 565 ) = 11386.6800000000
  Mref_V%vsh_ref( 565 ) = 6364.20000000000
  Mref_V%eta_ref( 565 ) = 1.00000000000000

  Mref_V%radius_ref( 566 ) = 5429519.00000000
  Mref_V%density_ref( 566 ) = 4540.82000000000
  Mref_V%vpv_ref( 566 ) = 11371.0700000000
  Mref_V%vsv_ref( 566 ) = 6357.91000000000
  Mref_V%Qkappa_ref( 566 ) = 57822.5000000000
  Mref_V%Qmu_ref( 566 ) = 355.000000000000
  Mref_V%vph_ref( 566 ) = 11371.0700000000
  Mref_V%vsh_ref( 566 ) = 6357.91000000000
  Mref_V%eta_ref( 566 ) = 1.00000000000000

  Mref_V%radius_ref( 567 ) = 5438990.00000000
  Mref_V%density_ref( 567 ) = 4534.84000000000
  Mref_V%vpv_ref( 567 ) = 11355.4100000000
  Mref_V%vsv_ref( 567 ) = 6351.59000000000
  Mref_V%Qkappa_ref( 567 ) = 57822.5000000000
  Mref_V%Qmu_ref( 567 ) = 355.000000000000
  Mref_V%vph_ref( 567 ) = 11355.4100000000
  Mref_V%vsh_ref( 567 ) = 6351.59000000000
  Mref_V%eta_ref( 567 ) = 1.00000000000000

  Mref_V%radius_ref( 568 ) = 5448462.00000000
  Mref_V%density_ref( 568 ) = 4528.85000000000
  Mref_V%vpv_ref( 568 ) = 11339.6900000000
  Mref_V%vsv_ref( 568 ) = 6345.25000000000
  Mref_V%Qkappa_ref( 568 ) = 57822.5000000000
  Mref_V%Qmu_ref( 568 ) = 355.000000000000
  Mref_V%vph_ref( 568 ) = 11339.6900000000
  Mref_V%vsh_ref( 568 ) = 6345.25000000000
  Mref_V%eta_ref( 568 ) = 1.00000000000000

  Mref_V%radius_ref( 569 ) = 5457933.00000000
  Mref_V%density_ref( 569 ) = 4522.81000000000
  Mref_V%vpv_ref( 569 ) = 11323.9100000000
  Mref_V%vsv_ref( 569 ) = 6338.88000000000
  Mref_V%Qkappa_ref( 569 ) = 57822.5000000000
  Mref_V%Qmu_ref( 569 ) = 355.000000000000
  Mref_V%vph_ref( 569 ) = 11323.9100000000
  Mref_V%vsh_ref( 569 ) = 6338.88000000000
  Mref_V%eta_ref( 569 ) = 1.00000000000000

  Mref_V%radius_ref( 570 ) = 5467404.00000000
  Mref_V%density_ref( 570 ) = 4516.76000000000
  Mref_V%vpv_ref( 570 ) = 11308.0900000000
  Mref_V%vsv_ref( 570 ) = 6332.49000000000
  Mref_V%Qkappa_ref( 570 ) = 57822.5000000000
  Mref_V%Qmu_ref( 570 ) = 355.000000000000
  Mref_V%vph_ref( 570 ) = 11308.0900000000
  Mref_V%vsh_ref( 570 ) = 6332.49000000000
  Mref_V%eta_ref( 570 ) = 1.00000000000000

  Mref_V%radius_ref( 571 ) = 5476875.00000000
  Mref_V%density_ref( 571 ) = 4510.67000000000
  Mref_V%vpv_ref( 571 ) = 11292.2000000000
  Mref_V%vsv_ref( 571 ) = 6326.05000000000
  Mref_V%Qkappa_ref( 571 ) = 57822.5000000000
  Mref_V%Qmu_ref( 571 ) = 355.000000000000
  Mref_V%vph_ref( 571 ) = 11292.2000000000
  Mref_V%vsh_ref( 571 ) = 6326.05000000000
  Mref_V%eta_ref( 571 ) = 1.00000000000000

  Mref_V%radius_ref( 572 ) = 5486346.00000000
  Mref_V%density_ref( 572 ) = 4504.56000000000
  Mref_V%vpv_ref( 572 ) = 11276.2500000000
  Mref_V%vsv_ref( 572 ) = 6319.60000000000
  Mref_V%Qkappa_ref( 572 ) = 57822.5000000000
  Mref_V%Qmu_ref( 572 ) = 355.000000000000
  Mref_V%vph_ref( 572 ) = 11276.2500000000
  Mref_V%vsh_ref( 572 ) = 6319.60000000000
  Mref_V%eta_ref( 572 ) = 1.00000000000000

  Mref_V%radius_ref( 573 ) = 5495817.00000000
  Mref_V%density_ref( 573 ) = 4498.41000000000
  Mref_V%vpv_ref( 573 ) = 11260.2500000000
  Mref_V%vsv_ref( 573 ) = 6313.13000000000
  Mref_V%Qkappa_ref( 573 ) = 57822.5000000000
  Mref_V%Qmu_ref( 573 ) = 355.000000000000
  Mref_V%vph_ref( 573 ) = 11260.2500000000
  Mref_V%vsh_ref( 573 ) = 6313.13000000000
  Mref_V%eta_ref( 573 ) = 1.00000000000000

  Mref_V%radius_ref( 574 ) = 5505288.00000000
  Mref_V%density_ref( 574 ) = 4492.23000000000
  Mref_V%vpv_ref( 574 ) = 11244.1900000000
  Mref_V%vsv_ref( 574 ) = 6306.62000000000
  Mref_V%Qkappa_ref( 574 ) = 57822.5000000000
  Mref_V%Qmu_ref( 574 ) = 355.000000000000
  Mref_V%vph_ref( 574 ) = 11244.1900000000
  Mref_V%vsh_ref( 574 ) = 6306.62000000000
  Mref_V%eta_ref( 574 ) = 1.00000000000000

  Mref_V%radius_ref( 575 ) = 5514760.00000000
  Mref_V%density_ref( 575 ) = 4486.02000000000
  Mref_V%vpv_ref( 575 ) = 11228.0800000000
  Mref_V%vsv_ref( 575 ) = 6300.08000000000
  Mref_V%Qkappa_ref( 575 ) = 57822.5000000000
  Mref_V%Qmu_ref( 575 ) = 355.000000000000
  Mref_V%vph_ref( 575 ) = 11228.0800000000
  Mref_V%vsh_ref( 575 ) = 6300.08000000000
  Mref_V%eta_ref( 575 ) = 1.00000000000000

  Mref_V%radius_ref( 576 ) = 5524231.00000000
  Mref_V%density_ref( 576 ) = 4479.78000000000
  Mref_V%vpv_ref( 576 ) = 11211.9000000000
  Mref_V%vsv_ref( 576 ) = 6293.52000000000
  Mref_V%Qkappa_ref( 576 ) = 57822.5000000000
  Mref_V%Qmu_ref( 576 ) = 355.000000000000
  Mref_V%vph_ref( 576 ) = 11211.9000000000
  Mref_V%vsh_ref( 576 ) = 6293.52000000000
  Mref_V%eta_ref( 576 ) = 1.00000000000000

  Mref_V%radius_ref( 577 ) = 5533702.00000000
  Mref_V%density_ref( 577 ) = 4473.51000000000
  Mref_V%vpv_ref( 577 ) = 11195.6700000000
  Mref_V%vsv_ref( 577 ) = 6286.92000000000
  Mref_V%Qkappa_ref( 577 ) = 57822.5000000000
  Mref_V%Qmu_ref( 577 ) = 355.000000000000
  Mref_V%vph_ref( 577 ) = 11195.6700000000
  Mref_V%vsh_ref( 577 ) = 6286.92000000000
  Mref_V%eta_ref( 577 ) = 1.00000000000000

  Mref_V%radius_ref( 578 ) = 5543173.00000000
  Mref_V%density_ref( 578 ) = 4467.20000000000
  Mref_V%vpv_ref( 578 ) = 11179.3800000000
  Mref_V%vsv_ref( 578 ) = 6280.29000000000
  Mref_V%Qkappa_ref( 578 ) = 57822.5000000000
  Mref_V%Qmu_ref( 578 ) = 355.000000000000
  Mref_V%vph_ref( 578 ) = 11179.3800000000
  Mref_V%vsh_ref( 578 ) = 6280.29000000000
  Mref_V%eta_ref( 578 ) = 1.00000000000000

  Mref_V%radius_ref( 579 ) = 5552644.00000000
  Mref_V%density_ref( 579 ) = 4460.87000000000
  Mref_V%vpv_ref( 579 ) = 11163.0400000000
  Mref_V%vsv_ref( 579 ) = 6273.64000000000
  Mref_V%Qkappa_ref( 579 ) = 57822.5000000000
  Mref_V%Qmu_ref( 579 ) = 355.000000000000
  Mref_V%vph_ref( 579 ) = 11163.0400000000
  Mref_V%vsh_ref( 579 ) = 6273.64000000000
  Mref_V%eta_ref( 579 ) = 1.00000000000000

  Mref_V%radius_ref( 580 ) = 5562115.00000000
  Mref_V%density_ref( 580 ) = 4454.49000000000
  Mref_V%vpv_ref( 580 ) = 11146.6300000000
  Mref_V%vsv_ref( 580 ) = 6266.96000000000
  Mref_V%Qkappa_ref( 580 ) = 57822.5000000000
  Mref_V%Qmu_ref( 580 ) = 355.000000000000
  Mref_V%vph_ref( 580 ) = 11146.6300000000
  Mref_V%vsh_ref( 580 ) = 6266.96000000000
  Mref_V%eta_ref( 580 ) = 1.00000000000000

  Mref_V%radius_ref( 581 ) = 5571587.00000000
  Mref_V%density_ref( 581 ) = 4448.08000000000
  Mref_V%vpv_ref( 581 ) = 11130.1800000000
  Mref_V%vsv_ref( 581 ) = 6260.25000000000
  Mref_V%Qkappa_ref( 581 ) = 57822.5000000000
  Mref_V%Qmu_ref( 581 ) = 355.000000000000
  Mref_V%vph_ref( 581 ) = 11130.1800000000
  Mref_V%vsh_ref( 581 ) = 6260.25000000000
  Mref_V%eta_ref( 581 ) = 1.00000000000000

  Mref_V%radius_ref( 582 ) = 5581058.00000000
  Mref_V%density_ref( 582 ) = 4441.63000000000
  Mref_V%vpv_ref( 582 ) = 11113.6700000000
  Mref_V%vsv_ref( 582 ) = 6253.51000000000
  Mref_V%Qkappa_ref( 582 ) = 57822.5000000000
  Mref_V%Qmu_ref( 582 ) = 355.000000000000
  Mref_V%vph_ref( 582 ) = 11113.6700000000
  Mref_V%vsh_ref( 582 ) = 6253.51000000000
  Mref_V%eta_ref( 582 ) = 1.00000000000000

  Mref_V%radius_ref( 583 ) = 5590529.00000000
  Mref_V%density_ref( 583 ) = 4435.13000000000
  Mref_V%vpv_ref( 583 ) = 11097.1100000000
  Mref_V%vsv_ref( 583 ) = 6246.75000000000
  Mref_V%Qkappa_ref( 583 ) = 57822.5000000000
  Mref_V%Qmu_ref( 583 ) = 355.000000000000
  Mref_V%vph_ref( 583 ) = 11097.1100000000
  Mref_V%vsh_ref( 583 ) = 6246.75000000000
  Mref_V%eta_ref( 583 ) = 1.00000000000000

  Mref_V%radius_ref( 584 ) = 5600000.00000000
  Mref_V%density_ref( 584 ) = 4428.60000000000
  Mref_V%vpv_ref( 584 ) = 11080.5100000000
  Mref_V%vsv_ref( 584 ) = 6239.95000000000
  Mref_V%Qkappa_ref( 584 ) = 57822.5000000000
  Mref_V%Qmu_ref( 584 ) = 355.000000000000
  Mref_V%vph_ref( 584 ) = 11080.5100000000
  Mref_V%vsh_ref( 584 ) = 6239.95000000000
  Mref_V%eta_ref( 584 ) = 1.00000000000000

  Mref_V%radius_ref( 585 ) = 5600000.00000000
  Mref_V%density_ref( 585 ) = 4428.59000000000
  Mref_V%vpv_ref( 585 ) = 11080.5100000000
  Mref_V%vsv_ref( 585 ) = 6239.95000000000
  Mref_V%Qkappa_ref( 585 ) = 57822.5000000000
  Mref_V%Qmu_ref( 585 ) = 355.000000000000
  Mref_V%vph_ref( 585 ) = 11080.5100000000
  Mref_V%vsh_ref( 585 ) = 6239.95000000000
  Mref_V%eta_ref( 585 ) = 1.00000000000000

  Mref_V%radius_ref( 586 ) = 5607562.00000000
  Mref_V%density_ref( 586 ) = 4423.32000000000
  Mref_V%vpv_ref( 586 ) = 11063.0100000000
  Mref_V%vsv_ref( 586 ) = 6219.68000000000
  Mref_V%Qkappa_ref( 586 ) = 57822.5000000000
  Mref_V%Qmu_ref( 586 ) = 355.000000000000
  Mref_V%vph_ref( 586 ) = 11063.0100000000
  Mref_V%vsh_ref( 586 ) = 6219.68000000000
  Mref_V%eta_ref( 586 ) = 1.00000000000000

  Mref_V%radius_ref( 587 ) = 5615125.00000000
  Mref_V%density_ref( 587 ) = 4418.01000000000
  Mref_V%vpv_ref( 587 ) = 11045.2200000000
  Mref_V%vsv_ref( 587 ) = 6200.29000000000
  Mref_V%Qkappa_ref( 587 ) = 57822.5000000000
  Mref_V%Qmu_ref( 587 ) = 355.000000000000
  Mref_V%vph_ref( 587 ) = 11045.2200000000
  Mref_V%vsh_ref( 587 ) = 6200.29000000000
  Mref_V%eta_ref( 587 ) = 1.00000000000000

  Mref_V%radius_ref( 588 ) = 5622688.00000000
  Mref_V%density_ref( 588 ) = 4412.67000000000
  Mref_V%vpv_ref( 588 ) = 11026.8200000000
  Mref_V%vsv_ref( 588 ) = 6181.16000000000
  Mref_V%Qkappa_ref( 588 ) = 57822.5000000000
  Mref_V%Qmu_ref( 588 ) = 355.000000000000
  Mref_V%vph_ref( 588 ) = 11026.8200000000
  Mref_V%vsh_ref( 588 ) = 6181.16000000000
  Mref_V%eta_ref( 588 ) = 1.00000000000000

  Mref_V%radius_ref( 589 ) = 5630250.00000000
  Mref_V%density_ref( 589 ) = 4407.30000000000
  Mref_V%vpv_ref( 589 ) = 11008.4700000000
  Mref_V%vsv_ref( 589 ) = 6162.04000000000
  Mref_V%Qkappa_ref( 589 ) = 57822.5000000000
  Mref_V%Qmu_ref( 589 ) = 355.000000000000
  Mref_V%vph_ref( 589 ) = 11008.4700000000
  Mref_V%vsh_ref( 589 ) = 6162.04000000000
  Mref_V%eta_ref( 589 ) = 1.00000000000000

  Mref_V%radius_ref( 590 ) = 5637812.00000000
  Mref_V%density_ref( 590 ) = 4401.90000000000
  Mref_V%vpv_ref( 590 ) = 10989.0400000000
  Mref_V%vsv_ref( 590 ) = 6143.01000000000
  Mref_V%Qkappa_ref( 590 ) = 57822.5000000000
  Mref_V%Qmu_ref( 590 ) = 355.000000000000
  Mref_V%vph_ref( 590 ) = 10989.0400000000
  Mref_V%vsh_ref( 590 ) = 6143.01000000000
  Mref_V%eta_ref( 590 ) = 1.00000000000000

  Mref_V%radius_ref( 591 ) = 5645375.00000000
  Mref_V%density_ref( 591 ) = 4396.45000000000
  Mref_V%vpv_ref( 591 ) = 10969.6300000000
  Mref_V%vsv_ref( 591 ) = 6123.98000000000
  Mref_V%Qkappa_ref( 591 ) = 57822.5000000000
  Mref_V%Qmu_ref( 591 ) = 355.000000000000
  Mref_V%vph_ref( 591 ) = 10969.6300000000
  Mref_V%vsh_ref( 591 ) = 6123.98000000000
  Mref_V%eta_ref( 591 ) = 1.00000000000000

  Mref_V%radius_ref( 592 ) = 5652938.00000000
  Mref_V%density_ref( 592 ) = 4390.96000000000
  Mref_V%vpv_ref( 592 ) = 10948.7600000000
  Mref_V%vsv_ref( 592 ) = 6103.71000000000
  Mref_V%Qkappa_ref( 592 ) = 57822.5000000000
  Mref_V%Qmu_ref( 592 ) = 355.000000000000
  Mref_V%vph_ref( 592 ) = 10948.7600000000
  Mref_V%vsh_ref( 592 ) = 6103.71000000000
  Mref_V%eta_ref( 592 ) = 1.00000000000000

  Mref_V%radius_ref( 593 ) = 5660500.00000000
  Mref_V%density_ref( 593 ) = 4385.40000000000
  Mref_V%vpv_ref( 593 ) = 10928.0200000000
  Mref_V%vsv_ref( 593 ) = 6083.53000000000
  Mref_V%Qkappa_ref( 593 ) = 57822.5000000000
  Mref_V%Qmu_ref( 593 ) = 355.000000000000
  Mref_V%vph_ref( 593 ) = 10928.0200000000
  Mref_V%vsh_ref( 593 ) = 6083.53000000000
  Mref_V%eta_ref( 593 ) = 1.00000000000000

  Mref_V%radius_ref( 594 ) = 5668062.00000000
  Mref_V%density_ref( 594 ) = 4379.81000000000
  Mref_V%vpv_ref( 594 ) = 10907.4200000000
  Mref_V%vsv_ref( 594 ) = 6063.45000000000
  Mref_V%Qkappa_ref( 594 ) = 57822.5000000000
  Mref_V%Qmu_ref( 594 ) = 355.000000000000
  Mref_V%vph_ref( 594 ) = 10907.4200000000
  Mref_V%vsh_ref( 594 ) = 6063.45000000000
  Mref_V%eta_ref( 594 ) = 1.00000000000000

  Mref_V%radius_ref( 595 ) = 5675625.00000000
  Mref_V%density_ref( 595 ) = 4374.16000000000
  Mref_V%vpv_ref( 595 ) = 10886.9400000000
  Mref_V%vsv_ref( 595 ) = 6043.44000000000
  Mref_V%Qkappa_ref( 595 ) = 57822.5000000000
  Mref_V%Qmu_ref( 595 ) = 355.000000000000
  Mref_V%vph_ref( 595 ) = 10886.9400000000
  Mref_V%vsh_ref( 595 ) = 6043.44000000000
  Mref_V%eta_ref( 595 ) = 1.00000000000000

  Mref_V%radius_ref( 596 ) = 5683188.00000000
  Mref_V%density_ref( 596 ) = 4368.47000000000
  Mref_V%vpv_ref( 596 ) = 10866.6000000000
  Mref_V%vsv_ref( 596 ) = 6023.52000000000
  Mref_V%Qkappa_ref( 596 ) = 57822.5000000000
  Mref_V%Qmu_ref( 596 ) = 355.000000000000
  Mref_V%vph_ref( 596 ) = 10866.6000000000
  Mref_V%vsh_ref( 596 ) = 6023.52000000000
  Mref_V%eta_ref( 596 ) = 1.00000000000000

  Mref_V%radius_ref( 597 ) = 5690750.00000000
  Mref_V%density_ref( 597 ) = 4362.72000000000
  Mref_V%vpv_ref( 597 ) = 10846.4100000000
  Mref_V%vsv_ref( 597 ) = 6003.73000000000
  Mref_V%Qkappa_ref( 597 ) = 57822.5000000000
  Mref_V%Qmu_ref( 597 ) = 355.000000000000
  Mref_V%vph_ref( 597 ) = 10846.4100000000
  Mref_V%vsh_ref( 597 ) = 6003.73000000000
  Mref_V%eta_ref( 597 ) = 1.00000000000000

  Mref_V%radius_ref( 598 ) = 5698312.00000000
  Mref_V%density_ref( 598 ) = 4356.92000000000
  Mref_V%vpv_ref( 598 ) = 10826.3500000000
  Mref_V%vsv_ref( 598 ) = 5984.03000000000
  Mref_V%Qkappa_ref( 598 ) = 57822.5000000000
  Mref_V%Qmu_ref( 598 ) = 355.000000000000
  Mref_V%vph_ref( 598 ) = 10826.3500000000
  Mref_V%vsh_ref( 598 ) = 5984.03000000000
  Mref_V%eta_ref( 598 ) = 1.00000000000000

  Mref_V%radius_ref( 599 ) = 5705875.00000000
  Mref_V%density_ref( 599 ) = 4351.08000000000
  Mref_V%vpv_ref( 599 ) = 10806.4200000000
  Mref_V%vsv_ref( 599 ) = 5964.38000000000
  Mref_V%Qkappa_ref( 599 ) = 57822.5000000000
  Mref_V%Qmu_ref( 599 ) = 355.000000000000
  Mref_V%vph_ref( 599 ) = 10806.4200000000
  Mref_V%vsh_ref( 599 ) = 5964.38000000000
  Mref_V%eta_ref( 599 ) = 1.00000000000000

  Mref_V%radius_ref( 600 ) = 5713438.00000000
  Mref_V%density_ref( 600 ) = 4345.18000000000
  Mref_V%vpv_ref( 600 ) = 10786.6100000000
  Mref_V%vsv_ref( 600 ) = 5944.81000000000
  Mref_V%Qkappa_ref( 600 ) = 57822.5000000000
  Mref_V%Qmu_ref( 600 ) = 355.000000000000
  Mref_V%vph_ref( 600 ) = 10786.6100000000
  Mref_V%vsh_ref( 600 ) = 5944.81000000000
  Mref_V%eta_ref( 600 ) = 1.00000000000000

  Mref_V%radius_ref( 601 ) = 5721000.00000000
  Mref_V%density_ref( 601 ) = 4339.24000000000
  Mref_V%vpv_ref( 601 ) = 10766.9000000000
  Mref_V%vsv_ref( 601 ) = 5925.27000000000
  Mref_V%Qkappa_ref( 601 ) = 57822.5000000000
  Mref_V%Qmu_ref( 601 ) = 355.000000000000
  Mref_V%vph_ref( 601 ) = 10766.9000000000
  Mref_V%vsh_ref( 601 ) = 5925.27000000000
  Mref_V%eta_ref( 601 ) = 1.00000000000000

  Mref_V%radius_ref( 602 ) = 5721000.00000000
  Mref_V%density_ref( 602 ) = 4047.01000000000
  Mref_V%vpv_ref( 602 ) = 10278.8800000000
  Mref_V%vsv_ref( 602 ) = 5550.32000000000
  Mref_V%Qkappa_ref( 602 ) = 943.000000000000
  Mref_V%Qmu_ref( 602 ) = 165.000000000000
  Mref_V%vph_ref( 602 ) = 10278.8800000000
  Mref_V%vsh_ref( 602 ) = 5550.32000000000
  Mref_V%eta_ref( 602 ) = 1.00000000000000

  Mref_V%radius_ref( 603 ) = 5724572.00000000
  Mref_V%density_ref( 603 ) = 4042.50000000000
  Mref_V%vpv_ref( 603 ) = 10261.8700000000
  Mref_V%vsv_ref( 603 ) = 5541.20000000000
  Mref_V%Qkappa_ref( 603 ) = 943.000000000000
  Mref_V%Qmu_ref( 603 ) = 165.000000000000
  Mref_V%vph_ref( 603 ) = 10261.8700000000
  Mref_V%vsh_ref( 603 ) = 5541.20000000000
  Mref_V%eta_ref( 603 ) = 1.00000000000000

  Mref_V%radius_ref( 604 ) = 5728143.00000000
  Mref_V%density_ref( 604 ) = 4037.98000000000
  Mref_V%vpv_ref( 604 ) = 10244.8400000000
  Mref_V%vsv_ref( 604 ) = 5532.08000000000
  Mref_V%Qkappa_ref( 604 ) = 943.000000000000
  Mref_V%Qmu_ref( 604 ) = 165.000000000000
  Mref_V%vph_ref( 604 ) = 10244.8400000000
  Mref_V%vsh_ref( 604 ) = 5532.08000000000
  Mref_V%eta_ref( 604 ) = 1.00000000000000

  Mref_V%radius_ref( 605 ) = 5731714.00000000
  Mref_V%density_ref( 605 ) = 4033.48000000000
  Mref_V%vpv_ref( 605 ) = 10227.8200000000
  Mref_V%vsv_ref( 605 ) = 5522.96000000000
  Mref_V%Qkappa_ref( 605 ) = 943.000000000000
  Mref_V%Qmu_ref( 605 ) = 165.000000000000
  Mref_V%vph_ref( 605 ) = 10227.8200000000
  Mref_V%vsh_ref( 605 ) = 5522.96000000000
  Mref_V%eta_ref( 605 ) = 1.00000000000000

  Mref_V%radius_ref( 606 ) = 5735286.00000000
  Mref_V%density_ref( 606 ) = 4028.95000000000
  Mref_V%vpv_ref( 606 ) = 10210.8000000000
  Mref_V%vsv_ref( 606 ) = 5513.83000000000
  Mref_V%Qkappa_ref( 606 ) = 943.000000000000
  Mref_V%Qmu_ref( 606 ) = 165.000000000000
  Mref_V%vph_ref( 606 ) = 10210.8000000000
  Mref_V%vsh_ref( 606 ) = 5513.83000000000
  Mref_V%eta_ref( 606 ) = 1.00000000000000

  Mref_V%radius_ref( 607 ) = 5738857.00000000
  Mref_V%density_ref( 607 ) = 4024.45000000000
  Mref_V%vpv_ref( 607 ) = 10193.7800000000
  Mref_V%vsv_ref( 607 ) = 5504.71000000000
  Mref_V%Qkappa_ref( 607 ) = 943.000000000000
  Mref_V%Qmu_ref( 607 ) = 165.000000000000
  Mref_V%vph_ref( 607 ) = 10193.7800000000
  Mref_V%vsh_ref( 607 ) = 5504.71000000000
  Mref_V%eta_ref( 607 ) = 1.00000000000000

  Mref_V%radius_ref( 608 ) = 5742428.00000000
  Mref_V%density_ref( 608 ) = 4019.93000000000
  Mref_V%vpv_ref( 608 ) = 10176.7700000000
  Mref_V%vsv_ref( 608 ) = 5495.59000000000
  Mref_V%Qkappa_ref( 608 ) = 943.000000000000
  Mref_V%Qmu_ref( 608 ) = 165.000000000000
  Mref_V%vph_ref( 608 ) = 10176.7700000000
  Mref_V%vsh_ref( 608 ) = 5495.59000000000
  Mref_V%eta_ref( 608 ) = 1.00000000000000

  Mref_V%radius_ref( 609 ) = 5746000.00000000
  Mref_V%density_ref( 609 ) = 4015.42000000000
  Mref_V%vpv_ref( 609 ) = 10159.7400000000
  Mref_V%vsv_ref( 609 ) = 5486.47000000000
  Mref_V%Qkappa_ref( 609 ) = 943.000000000000
  Mref_V%Qmu_ref( 609 ) = 165.000000000000
  Mref_V%vph_ref( 609 ) = 10159.7400000000
  Mref_V%vsh_ref( 609 ) = 5486.47000000000
  Mref_V%eta_ref( 609 ) = 1.00000000000000

  Mref_V%radius_ref( 610 ) = 5749572.00000000
  Mref_V%density_ref( 610 ) = 4010.90000000000
  Mref_V%vpv_ref( 610 ) = 10142.7200000000
  Mref_V%vsv_ref( 610 ) = 5477.35000000000
  Mref_V%Qkappa_ref( 610 ) = 943.000000000000
  Mref_V%Qmu_ref( 610 ) = 165.000000000000
  Mref_V%vph_ref( 610 ) = 10142.7200000000
  Mref_V%vsh_ref( 610 ) = 5477.35000000000
  Mref_V%eta_ref( 610 ) = 1.00000000000000

  Mref_V%radius_ref( 611 ) = 5753143.00000000
  Mref_V%density_ref( 611 ) = 4006.38000000000
  Mref_V%vpv_ref( 611 ) = 10125.7100000000
  Mref_V%vsv_ref( 611 ) = 5468.22000000000
  Mref_V%Qkappa_ref( 611 ) = 943.000000000000
  Mref_V%Qmu_ref( 611 ) = 165.000000000000
  Mref_V%vph_ref( 611 ) = 10125.7100000000
  Mref_V%vsh_ref( 611 ) = 5468.22000000000
  Mref_V%eta_ref( 611 ) = 1.00000000000000

  Mref_V%radius_ref( 612 ) = 5756714.00000000
  Mref_V%density_ref( 612 ) = 4001.87000000000
  Mref_V%vpv_ref( 612 ) = 10108.7000000000
  Mref_V%vsv_ref( 612 ) = 5459.10000000000
  Mref_V%Qkappa_ref( 612 ) = 943.000000000000
  Mref_V%Qmu_ref( 612 ) = 165.000000000000
  Mref_V%vph_ref( 612 ) = 10108.7000000000
  Mref_V%vsh_ref( 612 ) = 5459.10000000000
  Mref_V%eta_ref( 612 ) = 1.00000000000000

  Mref_V%radius_ref( 613 ) = 5760286.00000000
  Mref_V%density_ref( 613 ) = 3997.35000000000
  Mref_V%vpv_ref( 613 ) = 10091.6800000000
  Mref_V%vsv_ref( 613 ) = 5449.97000000000
  Mref_V%Qkappa_ref( 613 ) = 943.000000000000
  Mref_V%Qmu_ref( 613 ) = 165.000000000000
  Mref_V%vph_ref( 613 ) = 10091.6800000000
  Mref_V%vsh_ref( 613 ) = 5449.97000000000
  Mref_V%eta_ref( 613 ) = 1.00000000000000

  Mref_V%radius_ref( 614 ) = 5763857.00000000
  Mref_V%density_ref( 614 ) = 3992.84000000000
  Mref_V%vpv_ref( 614 ) = 10074.6800000000
  Mref_V%vsv_ref( 614 ) = 5440.84000000000
  Mref_V%Qkappa_ref( 614 ) = 943.000000000000
  Mref_V%Qmu_ref( 614 ) = 165.000000000000
  Mref_V%vph_ref( 614 ) = 10074.6800000000
  Mref_V%vsh_ref( 614 ) = 5440.84000000000
  Mref_V%eta_ref( 614 ) = 1.00000000000000

  Mref_V%radius_ref( 615 ) = 5767428.00000000
  Mref_V%density_ref( 615 ) = 3988.32000000000
  Mref_V%vpv_ref( 615 ) = 10057.6800000000
  Mref_V%vsv_ref( 615 ) = 5431.71000000000
  Mref_V%Qkappa_ref( 615 ) = 943.000000000000
  Mref_V%Qmu_ref( 615 ) = 165.000000000000
  Mref_V%vph_ref( 615 ) = 10057.6800000000
  Mref_V%vsh_ref( 615 ) = 5431.71000000000
  Mref_V%eta_ref( 615 ) = 1.00000000000000

  Mref_V%radius_ref( 616 ) = 5771000.00000000
  Mref_V%density_ref( 616 ) = 3983.80000000000
  Mref_V%vpv_ref( 616 ) = 10040.6400000000
  Mref_V%vsv_ref( 616 ) = 5422.57000000000
  Mref_V%Qkappa_ref( 616 ) = 943.000000000000
  Mref_V%Qmu_ref( 616 ) = 165.000000000000
  Mref_V%vph_ref( 616 ) = 10040.6400000000
  Mref_V%vsh_ref( 616 ) = 5422.57000000000
  Mref_V%eta_ref( 616 ) = 1.00000000000000

  Mref_V%radius_ref( 617 ) = 5771000.00000000
  Mref_V%density_ref( 617 ) = 3983.80000000000
  Mref_V%vpv_ref( 617 ) = 10040.6700000000
  Mref_V%vsv_ref( 617 ) = 5422.59000000000
  Mref_V%Qkappa_ref( 617 ) = 943.000000000000
  Mref_V%Qmu_ref( 617 ) = 165.000000000000
  Mref_V%vph_ref( 617 ) = 10040.6700000000
  Mref_V%vsh_ref( 617 ) = 5422.59000000000
  Mref_V%eta_ref( 617 ) = 1.00000000000000

  Mref_V%radius_ref( 618 ) = 5777334.00000000
  Mref_V%density_ref( 618 ) = 3975.79000000000
  Mref_V%vpv_ref( 618 ) = 10010.5200000000
  Mref_V%vsv_ref( 618 ) = 5406.39000000000
  Mref_V%Qkappa_ref( 618 ) = 943.000000000000
  Mref_V%Qmu_ref( 618 ) = 165.000000000000
  Mref_V%vph_ref( 618 ) = 10010.5200000000
  Mref_V%vsh_ref( 618 ) = 5406.39000000000
  Mref_V%eta_ref( 618 ) = 1.00000000000000

  Mref_V%radius_ref( 619 ) = 5783666.00000000
  Mref_V%density_ref( 619 ) = 3967.77000000000
  Mref_V%vpv_ref( 619 ) = 9980.51000000000
  Mref_V%vsv_ref( 619 ) = 5390.30000000000
  Mref_V%Qkappa_ref( 619 ) = 943.000000000000
  Mref_V%Qmu_ref( 619 ) = 165.000000000000
  Mref_V%vph_ref( 619 ) = 9980.51000000000
  Mref_V%vsh_ref( 619 ) = 5390.30000000000
  Mref_V%eta_ref( 619 ) = 1.00000000000000

  Mref_V%radius_ref( 620 ) = 5790000.00000000
  Mref_V%density_ref( 620 ) = 3959.76000000000
  Mref_V%vpv_ref( 620 ) = 9950.64000000000
  Mref_V%vsv_ref( 620 ) = 5374.34000000000
  Mref_V%Qkappa_ref( 620 ) = 943.000000000000
  Mref_V%Qmu_ref( 620 ) = 165.000000000000
  Mref_V%vph_ref( 620 ) = 9950.64000000000
  Mref_V%vsh_ref( 620 ) = 5374.34000000000
  Mref_V%eta_ref( 620 ) = 1.00000000000000

  Mref_V%radius_ref( 621 ) = 5796334.00000000
  Mref_V%density_ref( 621 ) = 3951.75000000000
  Mref_V%vpv_ref( 621 ) = 9920.91000000000
  Mref_V%vsv_ref( 621 ) = 5358.52000000000
  Mref_V%Qkappa_ref( 621 ) = 943.000000000000
  Mref_V%Qmu_ref( 621 ) = 165.000000000000
  Mref_V%vph_ref( 621 ) = 9920.91000000000
  Mref_V%vsh_ref( 621 ) = 5358.52000000000
  Mref_V%eta_ref( 621 ) = 1.00000000000000

  Mref_V%radius_ref( 622 ) = 5802666.00000000
  Mref_V%density_ref( 622 ) = 3943.73000000000
  Mref_V%vpv_ref( 622 ) = 9891.35000000000
  Mref_V%vsv_ref( 622 ) = 5342.83000000000
  Mref_V%Qkappa_ref( 622 ) = 943.000000000000
  Mref_V%Qmu_ref( 622 ) = 165.000000000000
  Mref_V%vph_ref( 622 ) = 9891.35000000000
  Mref_V%vsh_ref( 622 ) = 5342.83000000000
  Mref_V%eta_ref( 622 ) = 1.00000000000000

  Mref_V%radius_ref( 623 ) = 5809000.00000000
  Mref_V%density_ref( 623 ) = 3935.71000000000
  Mref_V%vpv_ref( 623 ) = 9861.96000000000
  Mref_V%vsv_ref( 623 ) = 5327.31000000000
  Mref_V%Qkappa_ref( 623 ) = 943.000000000000
  Mref_V%Qmu_ref( 623 ) = 165.000000000000
  Mref_V%vph_ref( 623 ) = 9861.96000000000
  Mref_V%vsh_ref( 623 ) = 5327.31000000000
  Mref_V%eta_ref( 623 ) = 1.00000000000000

  Mref_V%radius_ref( 624 ) = 5815334.00000000
  Mref_V%density_ref( 624 ) = 3927.69000000000
  Mref_V%vpv_ref( 624 ) = 9832.79000000000
  Mref_V%vsv_ref( 624 ) = 5311.92000000000
  Mref_V%Qkappa_ref( 624 ) = 943.000000000000
  Mref_V%Qmu_ref( 624 ) = 165.000000000000
  Mref_V%vph_ref( 624 ) = 9832.79000000000
  Mref_V%vsh_ref( 624 ) = 5311.92000000000
  Mref_V%eta_ref( 624 ) = 1.00000000000000

  Mref_V%radius_ref( 625 ) = 5821666.00000000
  Mref_V%density_ref( 625 ) = 3919.67000000000
  Mref_V%vpv_ref( 625 ) = 9803.79000000000
  Mref_V%vsv_ref( 625 ) = 5296.73000000000
  Mref_V%Qkappa_ref( 625 ) = 943.000000000000
  Mref_V%Qmu_ref( 625 ) = 165.000000000000
  Mref_V%vph_ref( 625 ) = 9803.79000000000
  Mref_V%vsh_ref( 625 ) = 5296.73000000000
  Mref_V%eta_ref( 625 ) = 1.00000000000000

  Mref_V%radius_ref( 626 ) = 5828000.00000000
  Mref_V%density_ref( 626 ) = 3911.65000000000
  Mref_V%vpv_ref( 626 ) = 9774.98000000000
  Mref_V%vsv_ref( 626 ) = 5281.71000000000
  Mref_V%Qkappa_ref( 626 ) = 943.000000000000
  Mref_V%Qmu_ref( 626 ) = 165.000000000000
  Mref_V%vph_ref( 626 ) = 9774.98000000000
  Mref_V%vsh_ref( 626 ) = 5281.71000000000
  Mref_V%eta_ref( 626 ) = 1.00000000000000

  Mref_V%radius_ref( 627 ) = 5834334.00000000
  Mref_V%density_ref( 627 ) = 3903.61000000000
  Mref_V%vpv_ref( 627 ) = 9746.41000000000
  Mref_V%vsv_ref( 627 ) = 5266.86000000000
  Mref_V%Qkappa_ref( 627 ) = 943.000000000000
  Mref_V%Qmu_ref( 627 ) = 165.000000000000
  Mref_V%vph_ref( 627 ) = 9746.41000000000
  Mref_V%vsh_ref( 627 ) = 5266.86000000000
  Mref_V%eta_ref( 627 ) = 1.00000000000000

  Mref_V%radius_ref( 628 ) = 5840666.00000000
  Mref_V%density_ref( 628 ) = 3895.58000000000
  Mref_V%vpv_ref( 628 ) = 9718.08000000000
  Mref_V%vsv_ref( 628 ) = 5252.21000000000
  Mref_V%Qkappa_ref( 628 ) = 943.000000000000
  Mref_V%Qmu_ref( 628 ) = 165.000000000000
  Mref_V%vph_ref( 628 ) = 9718.08000000000
  Mref_V%vsh_ref( 628 ) = 5252.21000000000
  Mref_V%eta_ref( 628 ) = 1.00000000000000

  Mref_V%radius_ref( 629 ) = 5847000.00000000
  Mref_V%density_ref( 629 ) = 3887.56000000000
  Mref_V%vpv_ref( 629 ) = 9689.96000000000
  Mref_V%vsv_ref( 629 ) = 5237.78000000000
  Mref_V%Qkappa_ref( 629 ) = 943.000000000000
  Mref_V%Qmu_ref( 629 ) = 165.000000000000
  Mref_V%vph_ref( 629 ) = 9689.96000000000
  Mref_V%vsh_ref( 629 ) = 5237.78000000000
  Mref_V%eta_ref( 629 ) = 1.00000000000000

  Mref_V%radius_ref( 630 ) = 5853334.00000000
  Mref_V%density_ref( 630 ) = 3879.53000000000
  Mref_V%vpv_ref( 630 ) = 9662.10000000000
  Mref_V%vsv_ref( 630 ) = 5223.55000000000
  Mref_V%Qkappa_ref( 630 ) = 943.000000000000
  Mref_V%Qmu_ref( 630 ) = 165.000000000000
  Mref_V%vph_ref( 630 ) = 9662.10000000000
  Mref_V%vsh_ref( 630 ) = 5223.55000000000
  Mref_V%eta_ref( 630 ) = 1.00000000000000

  Mref_V%radius_ref( 631 ) = 5859666.00000000
  Mref_V%density_ref( 631 ) = 3871.50000000000
  Mref_V%vpv_ref( 631 ) = 9634.47000000000
  Mref_V%vsv_ref( 631 ) = 5209.54000000000
  Mref_V%Qkappa_ref( 631 ) = 943.000000000000
  Mref_V%Qmu_ref( 631 ) = 165.000000000000
  Mref_V%vph_ref( 631 ) = 9634.47000000000
  Mref_V%vsh_ref( 631 ) = 5209.54000000000
  Mref_V%eta_ref( 631 ) = 1.00000000000000

  Mref_V%radius_ref( 632 ) = 5866000.00000000
  Mref_V%density_ref( 632 ) = 3863.46000000000
  Mref_V%vpv_ref( 632 ) = 9607.11000000000
  Mref_V%vsv_ref( 632 ) = 5195.72000000000
  Mref_V%Qkappa_ref( 632 ) = 943.000000000000
  Mref_V%Qmu_ref( 632 ) = 165.000000000000
  Mref_V%vph_ref( 632 ) = 9607.11000000000
  Mref_V%vsh_ref( 632 ) = 5195.72000000000
  Mref_V%eta_ref( 632 ) = 1.00000000000000

  Mref_V%radius_ref( 633 ) = 5872334.00000000
  Mref_V%density_ref( 633 ) = 3855.42000000000
  Mref_V%vpv_ref( 633 ) = 9579.97000000000
  Mref_V%vsv_ref( 633 ) = 5182.10000000000
  Mref_V%Qkappa_ref( 633 ) = 943.000000000000
  Mref_V%Qmu_ref( 633 ) = 165.000000000000
  Mref_V%vph_ref( 633 ) = 9579.97000000000
  Mref_V%vsh_ref( 633 ) = 5182.10000000000
  Mref_V%eta_ref( 633 ) = 1.00000000000000

  Mref_V%radius_ref( 634 ) = 5878666.00000000
  Mref_V%density_ref( 634 ) = 3847.38000000000
  Mref_V%vpv_ref( 634 ) = 9553.08000000000
  Mref_V%vsv_ref( 634 ) = 5168.69000000000
  Mref_V%Qkappa_ref( 634 ) = 943.000000000000
  Mref_V%Qmu_ref( 634 ) = 165.000000000000
  Mref_V%vph_ref( 634 ) = 9553.08000000000
  Mref_V%vsh_ref( 634 ) = 5168.69000000000
  Mref_V%eta_ref( 634 ) = 1.00000000000000

  Mref_V%radius_ref( 635 ) = 5885000.00000000
  Mref_V%density_ref( 635 ) = 3839.33000000000
  Mref_V%vpv_ref( 635 ) = 9526.38000000000
  Mref_V%vsv_ref( 635 ) = 5155.42000000000
  Mref_V%Qkappa_ref( 635 ) = 943.000000000000
  Mref_V%Qmu_ref( 635 ) = 165.000000000000
  Mref_V%vph_ref( 635 ) = 9526.38000000000
  Mref_V%vsh_ref( 635 ) = 5155.42000000000
  Mref_V%eta_ref( 635 ) = 1.00000000000000

  Mref_V%radius_ref( 636 ) = 5891334.00000000
  Mref_V%density_ref( 636 ) = 3831.27000000000
  Mref_V%vpv_ref( 636 ) = 9499.78000000000
  Mref_V%vsv_ref( 636 ) = 5142.22000000000
  Mref_V%Qkappa_ref( 636 ) = 943.000000000000
  Mref_V%Qmu_ref( 636 ) = 165.000000000000
  Mref_V%vph_ref( 636 ) = 9499.78000000000
  Mref_V%vsh_ref( 636 ) = 5142.22000000000
  Mref_V%eta_ref( 636 ) = 1.00000000000000

  Mref_V%radius_ref( 637 ) = 5897666.00000000
  Mref_V%density_ref( 637 ) = 3823.22000000000
  Mref_V%vpv_ref( 637 ) = 9473.25000000000
  Mref_V%vsv_ref( 637 ) = 5129.05000000000
  Mref_V%Qkappa_ref( 637 ) = 943.000000000000
  Mref_V%Qmu_ref( 637 ) = 165.000000000000
  Mref_V%vph_ref( 637 ) = 9473.25000000000
  Mref_V%vsh_ref( 637 ) = 5129.05000000000
  Mref_V%eta_ref( 637 ) = 1.00000000000000

  Mref_V%radius_ref( 638 ) = 5904000.00000000
  Mref_V%density_ref( 638 ) = 3815.16000000000
  Mref_V%vpv_ref( 638 ) = 9446.74000000000
  Mref_V%vsv_ref( 638 ) = 5115.84000000000
  Mref_V%Qkappa_ref( 638 ) = 943.000000000000
  Mref_V%Qmu_ref( 638 ) = 165.000000000000
  Mref_V%vph_ref( 638 ) = 9446.74000000000
  Mref_V%vsh_ref( 638 ) = 5115.84000000000
  Mref_V%eta_ref( 638 ) = 1.00000000000000

  Mref_V%radius_ref( 639 ) = 5910334.00000000
  Mref_V%density_ref( 639 ) = 3807.09000000000
  Mref_V%vpv_ref( 639 ) = 9420.19000000000
  Mref_V%vsv_ref( 639 ) = 5102.55000000000
  Mref_V%Qkappa_ref( 639 ) = 943.000000000000
  Mref_V%Qmu_ref( 639 ) = 165.000000000000
  Mref_V%vph_ref( 639 ) = 9420.19000000000
  Mref_V%vsh_ref( 639 ) = 5102.55000000000
  Mref_V%eta_ref( 639 ) = 1.00000000000000

  Mref_V%radius_ref( 640 ) = 5916666.00000000
  Mref_V%density_ref( 640 ) = 3799.01000000000
  Mref_V%vpv_ref( 640 ) = 9393.55000000000
  Mref_V%vsv_ref( 640 ) = 5089.14000000000
  Mref_V%Qkappa_ref( 640 ) = 943.000000000000
  Mref_V%Qmu_ref( 640 ) = 165.000000000000
  Mref_V%vph_ref( 640 ) = 9393.55000000000
  Mref_V%vsh_ref( 640 ) = 5089.14000000000
  Mref_V%eta_ref( 640 ) = 1.00000000000000

  Mref_V%radius_ref( 641 ) = 5923000.00000000
  Mref_V%density_ref( 641 ) = 3790.94000000000
  Mref_V%vpv_ref( 641 ) = 9366.75000000000
  Mref_V%vsv_ref( 641 ) = 5075.50000000000
  Mref_V%Qkappa_ref( 641 ) = 943.000000000000
  Mref_V%Qmu_ref( 641 ) = 165.000000000000
  Mref_V%vph_ref( 641 ) = 9366.75000000000
  Mref_V%vsh_ref( 641 ) = 5075.50000000000
  Mref_V%eta_ref( 641 ) = 1.00000000000000

  Mref_V%radius_ref( 642 ) = 5929334.00000000
  Mref_V%density_ref( 642 ) = 3782.85000000000
  Mref_V%vpv_ref( 642 ) = 9339.76000000000
  Mref_V%vsv_ref( 642 ) = 5061.63000000000
  Mref_V%Qkappa_ref( 642 ) = 943.000000000000
  Mref_V%Qmu_ref( 642 ) = 165.000000000000
  Mref_V%vph_ref( 642 ) = 9339.76000000000
  Mref_V%vsh_ref( 642 ) = 5061.63000000000
  Mref_V%eta_ref( 642 ) = 1.00000000000000

  Mref_V%radius_ref( 643 ) = 5935666.00000000
  Mref_V%density_ref( 643 ) = 3774.78000000000
  Mref_V%vpv_ref( 643 ) = 9312.50000000000
  Mref_V%vsv_ref( 643 ) = 5047.46000000000
  Mref_V%Qkappa_ref( 643 ) = 943.000000000000
  Mref_V%Qmu_ref( 643 ) = 165.000000000000
  Mref_V%vph_ref( 643 ) = 9312.50000000000
  Mref_V%vsh_ref( 643 ) = 5047.46000000000
  Mref_V%eta_ref( 643 ) = 1.00000000000000

  Mref_V%radius_ref( 644 ) = 5942000.00000000
  Mref_V%density_ref( 644 ) = 3766.66000000000
  Mref_V%vpv_ref( 644 ) = 9284.96000000000
  Mref_V%vsv_ref( 644 ) = 5032.93000000000
  Mref_V%Qkappa_ref( 644 ) = 943.000000000000
  Mref_V%Qmu_ref( 644 ) = 165.000000000000
  Mref_V%vph_ref( 644 ) = 9284.96000000000
  Mref_V%vsh_ref( 644 ) = 5032.93000000000
  Mref_V%eta_ref( 644 ) = 1.00000000000000

  Mref_V%radius_ref( 645 ) = 5948334.00000000
  Mref_V%density_ref( 645 ) = 3758.56000000000
  Mref_V%vpv_ref( 645 ) = 9257.04000000000
  Mref_V%vsv_ref( 645 ) = 5018.03000000000
  Mref_V%Qkappa_ref( 645 ) = 943.000000000000
  Mref_V%Qmu_ref( 645 ) = 165.000000000000
  Mref_V%vph_ref( 645 ) = 9257.04000000000
  Mref_V%vsh_ref( 645 ) = 5018.03000000000
  Mref_V%eta_ref( 645 ) = 1.00000000000000

  Mref_V%radius_ref( 646 ) = 5954666.00000000
  Mref_V%density_ref( 646 ) = 3750.45000000000
  Mref_V%vpv_ref( 646 ) = 9228.73000000000
  Mref_V%vsv_ref( 646 ) = 5002.66000000000
  Mref_V%Qkappa_ref( 646 ) = 943.000000000000
  Mref_V%Qmu_ref( 646 ) = 165.000000000000
  Mref_V%vph_ref( 646 ) = 9228.73000000000
  Mref_V%vsh_ref( 646 ) = 5002.66000000000
  Mref_V%eta_ref( 646 ) = 1.00000000000000

  Mref_V%radius_ref( 647 ) = 5961000.00000000
  Mref_V%density_ref( 647 ) = 3742.34000000000
  Mref_V%vpv_ref( 647 ) = 9199.94000000000
  Mref_V%vsv_ref( 647 ) = 4986.77000000000
  Mref_V%Qkappa_ref( 647 ) = 943.000000000000
  Mref_V%Qmu_ref( 647 ) = 165.000000000000
  Mref_V%vph_ref( 647 ) = 9199.94000000000
  Mref_V%vsh_ref( 647 ) = 4986.77000000000
  Mref_V%eta_ref( 647 ) = 1.00000000000000

  Mref_V%radius_ref( 648 ) = 5961000.00000000
  Mref_V%density_ref( 648 ) = 3554.91000000000
  Mref_V%vpv_ref( 648 ) = 8940.94000000000
  Mref_V%vsv_ref( 648 ) = 4802.15000000000
  Mref_V%Qkappa_ref( 648 ) = 943.000000000000
  Mref_V%Qmu_ref( 648 ) = 165.000000000000
  Mref_V%vph_ref( 648 ) = 8940.94000000000
  Mref_V%vsh_ref( 648 ) = 4803.78000000000
  Mref_V%eta_ref( 648 ) = 1.00000000000000

  Mref_V%radius_ref( 649 ) = 5967334.00000000
  Mref_V%density_ref( 649 ) = 3551.00000000000
  Mref_V%vpv_ref( 649 ) = 8930.61000000000
  Mref_V%vsv_ref( 649 ) = 4798.23000000000
  Mref_V%Qkappa_ref( 649 ) = 943.000000000000
  Mref_V%Qmu_ref( 649 ) = 165.000000000000
  Mref_V%vph_ref( 649 ) = 8930.61000000000
  Mref_V%vsh_ref( 649 ) = 4800.54000000000
  Mref_V%eta_ref( 649 ) = 1.00000000000000

  Mref_V%radius_ref( 650 ) = 5973666.00000000
  Mref_V%density_ref( 650 ) = 3547.07000000000
  Mref_V%vpv_ref( 650 ) = 8920.22000000000
  Mref_V%vsv_ref( 650 ) = 4794.28000000000
  Mref_V%Qkappa_ref( 650 ) = 943.000000000000
  Mref_V%Qmu_ref( 650 ) = 165.000000000000
  Mref_V%vph_ref( 650 ) = 8920.22000000000
  Mref_V%vsh_ref( 650 ) = 4797.28000000000
  Mref_V%eta_ref( 650 ) = 1.00000000000000

  Mref_V%radius_ref( 651 ) = 5980000.00000000
  Mref_V%density_ref( 651 ) = 3543.16000000000
  Mref_V%vpv_ref( 651 ) = 8909.68000000000
  Mref_V%vsv_ref( 651 ) = 4790.38000000000
  Mref_V%Qkappa_ref( 651 ) = 943.000000000000
  Mref_V%Qmu_ref( 651 ) = 165.000000000000
  Mref_V%vph_ref( 651 ) = 8909.68000000000
  Mref_V%vsh_ref( 651 ) = 4793.96000000000
  Mref_V%eta_ref( 651 ) = 1.00000000000000

  Mref_V%radius_ref( 652 ) = 5986334.00000000
  Mref_V%density_ref( 652 ) = 3539.23000000000
  Mref_V%vpv_ref( 652 ) = 8898.47000000000
  Mref_V%vsv_ref( 652 ) = 4785.95000000000
  Mref_V%Qkappa_ref( 652 ) = 943.000000000000
  Mref_V%Qmu_ref( 652 ) = 165.000000000000
  Mref_V%vph_ref( 652 ) = 8898.47000000000
  Mref_V%vsh_ref( 652 ) = 4790.18000000000
  Mref_V%eta_ref( 652 ) = 0.999990000000000

  Mref_V%radius_ref( 653 ) = 5992666.00000000
  Mref_V%density_ref( 653 ) = 3535.32000000000
  Mref_V%vpv_ref( 653 ) = 8886.28000000000
  Mref_V%vsv_ref( 653 ) = 4780.83000000000
  Mref_V%Qkappa_ref( 653 ) = 943.000000000000
  Mref_V%Qmu_ref( 653 ) = 165.000000000000
  Mref_V%vph_ref( 653 ) = 8886.28000000000
  Mref_V%vsh_ref( 653 ) = 4785.78000000000
  Mref_V%eta_ref( 653 ) = 0.999970000000000

  Mref_V%radius_ref( 654 ) = 5999000.00000000
  Mref_V%density_ref( 654 ) = 3531.39000000000
  Mref_V%vpv_ref( 654 ) = 8873.03000000000
  Mref_V%vsv_ref( 654 ) = 4775.01000000000
  Mref_V%Qkappa_ref( 654 ) = 943.000000000000
  Mref_V%Qmu_ref( 654 ) = 165.000000000000
  Mref_V%vph_ref( 654 ) = 8873.03000000000
  Mref_V%vsh_ref( 654 ) = 4780.71000000000
  Mref_V%eta_ref( 654 ) = 0.999950000000000

  Mref_V%radius_ref( 655 ) = 6005334.00000000
  Mref_V%density_ref( 655 ) = 3527.46000000000
  Mref_V%vpv_ref( 655 ) = 8858.58000000000
  Mref_V%vsv_ref( 655 ) = 4768.45000000000
  Mref_V%Qkappa_ref( 655 ) = 943.000000000000
  Mref_V%Qmu_ref( 655 ) = 165.000000000000
  Mref_V%vph_ref( 655 ) = 8858.58000000000
  Mref_V%vsh_ref( 655 ) = 4775.00000000000
  Mref_V%eta_ref( 655 ) = 0.999940000000000

  Mref_V%radius_ref( 656 ) = 6011666.00000000
  Mref_V%density_ref( 656 ) = 3523.57000000000
  Mref_V%vpv_ref( 656 ) = 8842.82000000000
  Mref_V%vsv_ref( 656 ) = 4761.12000000000
  Mref_V%Qkappa_ref( 656 ) = 943.000000000000
  Mref_V%Qmu_ref( 656 ) = 165.000000000000
  Mref_V%vph_ref( 656 ) = 8842.82000000000
  Mref_V%vsh_ref( 656 ) = 4768.58000000000
  Mref_V%eta_ref( 656 ) = 0.999900000000000

  Mref_V%radius_ref( 657 ) = 6018000.00000000
  Mref_V%density_ref( 657 ) = 3519.67000000000
  Mref_V%vpv_ref( 657 ) = 8825.64000000000
  Mref_V%vsv_ref( 657 ) = 4752.97000000000
  Mref_V%Qkappa_ref( 657 ) = 943.000000000000
  Mref_V%Qmu_ref( 657 ) = 165.000000000000
  Mref_V%vph_ref( 657 ) = 8825.64000000000
  Mref_V%vsh_ref( 657 ) = 4761.41000000000
  Mref_V%eta_ref( 657 ) = 0.999860000000000

  Mref_V%radius_ref( 658 ) = 6024334.00000000
  Mref_V%density_ref( 658 ) = 3515.77000000000
  Mref_V%vpv_ref( 658 ) = 8806.94000000000
  Mref_V%vsv_ref( 658 ) = 4744.01000000000
  Mref_V%Qkappa_ref( 658 ) = 943.000000000000
  Mref_V%Qmu_ref( 658 ) = 165.000000000000
  Mref_V%vph_ref( 658 ) = 8806.94000000000
  Mref_V%vsh_ref( 658 ) = 4753.51000000000
  Mref_V%eta_ref( 658 ) = 0.999800000000000

  Mref_V%radius_ref( 659 ) = 6030666.00000000
  Mref_V%density_ref( 659 ) = 3511.91000000000
  Mref_V%vpv_ref( 659 ) = 8786.67000000000
  Mref_V%vsv_ref( 659 ) = 4734.25000000000
  Mref_V%Qkappa_ref( 659 ) = 943.000000000000
  Mref_V%Qmu_ref( 659 ) = 165.000000000000
  Mref_V%vph_ref( 659 ) = 8786.67000000000
  Mref_V%vsh_ref( 659 ) = 4744.86000000000
  Mref_V%eta_ref( 659 ) = 0.999740000000000

  Mref_V%radius_ref( 660 ) = 6037000.00000000
  Mref_V%density_ref( 660 ) = 3508.06000000000
  Mref_V%vpv_ref( 660 ) = 8764.85000000000
  Mref_V%vsv_ref( 660 ) = 4723.77000000000
  Mref_V%Qkappa_ref( 660 ) = 943.000000000000
  Mref_V%Qmu_ref( 660 ) = 165.000000000000
  Mref_V%vph_ref( 660 ) = 8764.85000000000
  Mref_V%vsh_ref( 660 ) = 4735.64000000000
  Mref_V%eta_ref( 660 ) = 0.999660000000000

  Mref_V%radius_ref( 661 ) = 6043334.00000000
  Mref_V%density_ref( 661 ) = 3504.21000000000
  Mref_V%vpv_ref( 661 ) = 8741.49000000000
  Mref_V%vsv_ref( 661 ) = 4712.70000000000
  Mref_V%Qkappa_ref( 661 ) = 943.000000000000
  Mref_V%Qmu_ref( 661 ) = 165.000000000000
  Mref_V%vph_ref( 661 ) = 8741.49000000000
  Mref_V%vsh_ref( 661 ) = 4725.88000000000
  Mref_V%eta_ref( 661 ) = 0.999570000000000

  Mref_V%radius_ref( 662 ) = 6049666.00000000
  Mref_V%density_ref( 662 ) = 3500.39000000000
  Mref_V%vpv_ref( 662 ) = 8716.63000000000
  Mref_V%vsv_ref( 662 ) = 4701.12000000000
  Mref_V%Qkappa_ref( 662 ) = 943.000000000000
  Mref_V%Qmu_ref( 662 ) = 165.000000000000
  Mref_V%vph_ref( 662 ) = 8716.63000000000
  Mref_V%vsh_ref( 662 ) = 4715.76000000000
  Mref_V%eta_ref( 662 ) = 0.999470000000000

  Mref_V%radius_ref( 663 ) = 6056000.00000000
  Mref_V%density_ref( 663 ) = 3496.58000000000
  Mref_V%vpv_ref( 663 ) = 8690.30000000000
  Mref_V%vsv_ref( 663 ) = 4689.11000000000
  Mref_V%Qkappa_ref( 663 ) = 943.000000000000
  Mref_V%Qmu_ref( 663 ) = 165.000000000000
  Mref_V%vph_ref( 663 ) = 8690.30000000000
  Mref_V%vsh_ref( 663 ) = 4705.34000000000
  Mref_V%eta_ref( 663 ) = 0.999340000000000

  Mref_V%radius_ref( 664 ) = 6062334.00000000
  Mref_V%density_ref( 664 ) = 3492.80000000000
  Mref_V%vpv_ref( 664 ) = 8662.50000000000
  Mref_V%vsv_ref( 664 ) = 4676.77000000000
  Mref_V%Qkappa_ref( 664 ) = 943.000000000000
  Mref_V%Qmu_ref( 664 ) = 165.000000000000
  Mref_V%vph_ref( 664 ) = 8662.50000000000
  Mref_V%vsh_ref( 664 ) = 4694.74000000000
  Mref_V%eta_ref( 664 ) = 0.999200000000000

  Mref_V%radius_ref( 665 ) = 6068666.00000000
  Mref_V%density_ref( 665 ) = 3489.05000000000
  Mref_V%vpv_ref( 665 ) = 8633.28000000000
  Mref_V%vsv_ref( 665 ) = 4664.20000000000
  Mref_V%Qkappa_ref( 665 ) = 943.000000000000
  Mref_V%Qmu_ref( 665 ) = 165.000000000000
  Mref_V%vph_ref( 665 ) = 8633.28000000000
  Mref_V%vsh_ref( 665 ) = 4684.08000000000
  Mref_V%eta_ref( 665 ) = 0.999040000000000

  Mref_V%radius_ref( 666 ) = 6075000.00000000
  Mref_V%density_ref( 666 ) = 3485.32000000000
  Mref_V%vpv_ref( 666 ) = 8602.66000000000
  Mref_V%vsv_ref( 666 ) = 4651.49000000000
  Mref_V%Qkappa_ref( 666 ) = 943.000000000000
  Mref_V%Qmu_ref( 666 ) = 165.000000000000
  Mref_V%vph_ref( 666 ) = 8602.66000000000
  Mref_V%vsh_ref( 666 ) = 4673.46000000000
  Mref_V%eta_ref( 666 ) = 0.998860000000000

  Mref_V%radius_ref( 667 ) = 6081334.00000000
  Mref_V%density_ref( 667 ) = 3481.61000000000
  Mref_V%vpv_ref( 667 ) = 8570.81000000000
  Mref_V%vsv_ref( 667 ) = 4638.69000000000
  Mref_V%Qkappa_ref( 667 ) = 943.000000000000
  Mref_V%Qmu_ref( 667 ) = 165.000000000000
  Mref_V%vph_ref( 667 ) = 8570.81000000000
  Mref_V%vsh_ref( 667 ) = 4662.94000000000
  Mref_V%eta_ref( 667 ) = 0.998640000000000

  Mref_V%radius_ref( 668 ) = 6087666.00000000
  Mref_V%density_ref( 668 ) = 3477.88000000000
  Mref_V%vpv_ref( 668 ) = 8538.06000000000
  Mref_V%vsv_ref( 668 ) = 4625.88000000000
  Mref_V%Qkappa_ref( 668 ) = 943.000000000000
  Mref_V%Qmu_ref( 668 ) = 165.000000000000
  Mref_V%vph_ref( 668 ) = 8538.06000000000
  Mref_V%vsh_ref( 668 ) = 4652.61000000000
  Mref_V%eta_ref( 668 ) = 0.998320000000000

  Mref_V%radius_ref( 669 ) = 6094000.00000000
  Mref_V%density_ref( 669 ) = 3474.16000000000
  Mref_V%vpv_ref( 669 ) = 8504.66000000000
  Mref_V%vsv_ref( 669 ) = 4613.07000000000
  Mref_V%Qkappa_ref( 669 ) = 943.000000000000
  Mref_V%Qmu_ref( 669 ) = 165.000000000000
  Mref_V%vph_ref( 669 ) = 8504.66000000000
  Mref_V%vsh_ref( 669 ) = 4642.55000000000
  Mref_V%eta_ref( 669 ) = 0.997900000000000

  Mref_V%radius_ref( 670 ) = 6100334.00000000
  Mref_V%density_ref( 670 ) = 3470.41000000000
  Mref_V%vpv_ref( 670 ) = 8470.92000000000
  Mref_V%vsv_ref( 670 ) = 4600.31000000000
  Mref_V%Qkappa_ref( 670 ) = 943.000000000000
  Mref_V%Qmu_ref( 670 ) = 165.000000000000
  Mref_V%vph_ref( 670 ) = 8470.92000000000
  Mref_V%vsh_ref( 670 ) = 4632.81000000000
  Mref_V%eta_ref( 670 ) = 0.997320000000000

  Mref_V%radius_ref( 671 ) = 6106666.00000000
  Mref_V%density_ref( 671 ) = 3466.59000000000
  Mref_V%vpv_ref( 671 ) = 8437.13000000000
  Mref_V%vsv_ref( 671 ) = 4587.67000000000
  Mref_V%Qkappa_ref( 671 ) = 943.000000000000
  Mref_V%Qmu_ref( 671 ) = 165.000000000000
  Mref_V%vph_ref( 671 ) = 8437.13000000000
  Mref_V%vsh_ref( 671 ) = 4623.51000000000
  Mref_V%eta_ref( 671 ) = 0.996540000000000

  Mref_V%radius_ref( 672 ) = 6113000.00000000
  Mref_V%density_ref( 672 ) = 3462.74000000000
  Mref_V%vpv_ref( 672 ) = 8403.52000000000
  Mref_V%vsv_ref( 672 ) = 4575.18000000000
  Mref_V%Qkappa_ref( 672 ) = 943.000000000000
  Mref_V%Qmu_ref( 672 ) = 165.000000000000
  Mref_V%vph_ref( 672 ) = 8403.52000000000
  Mref_V%vsh_ref( 672 ) = 4614.68000000000
  Mref_V%eta_ref( 672 ) = 0.995530000000000

  Mref_V%radius_ref( 673 ) = 6119334.00000000
  Mref_V%density_ref( 673 ) = 3458.78000000000
  Mref_V%vpv_ref( 673 ) = 8370.42000000000
  Mref_V%vsv_ref( 673 ) = 4562.88000000000
  Mref_V%Qkappa_ref( 673 ) = 943.000000000000
  Mref_V%Qmu_ref( 673 ) = 165.000000000000
  Mref_V%vph_ref( 673 ) = 8370.42000000000
  Mref_V%vsh_ref( 673 ) = 4606.39000000000
  Mref_V%eta_ref( 673 ) = 0.994260000000000

  Mref_V%radius_ref( 674 ) = 6125666.00000000
  Mref_V%density_ref( 674 ) = 3454.75000000000
  Mref_V%vpv_ref( 674 ) = 8338.11000000000
  Mref_V%vsv_ref( 674 ) = 4550.85000000000
  Mref_V%Qkappa_ref( 674 ) = 943.000000000000
  Mref_V%Qmu_ref( 674 ) = 165.000000000000
  Mref_V%vph_ref( 674 ) = 8338.11000000000
  Mref_V%vsh_ref( 674 ) = 4598.73000000000
  Mref_V%eta_ref( 674 ) = 0.992680000000000

  Mref_V%radius_ref( 675 ) = 6132000.00000000
  Mref_V%density_ref( 675 ) = 3450.61000000000
  Mref_V%vpv_ref( 675 ) = 8306.25000000000
  Mref_V%vsv_ref( 675 ) = 4539.08000000000
  Mref_V%Qkappa_ref( 675 ) = 943.000000000000
  Mref_V%Qmu_ref( 675 ) = 165.000000000000
  Mref_V%vph_ref( 675 ) = 8307.42000000000
  Mref_V%vsh_ref( 675 ) = 4591.76000000000
  Mref_V%eta_ref( 675 ) = 0.990750000000000

  Mref_V%radius_ref( 676 ) = 6138334.00000000
  Mref_V%density_ref( 676 ) = 3446.33000000000
  Mref_V%vpv_ref( 676 ) = 8275.42000000000
  Mref_V%vsv_ref( 676 ) = 4527.67000000000
  Mref_V%Qkappa_ref( 676 ) = 943.000000000000
  Mref_V%Qmu_ref( 676 ) = 165.000000000000
  Mref_V%vph_ref( 676 ) = 8278.36000000000
  Mref_V%vsh_ref( 676 ) = 4585.56000000000
  Mref_V%eta_ref( 676 ) = 0.988430000000000

  Mref_V%radius_ref( 677 ) = 6144666.00000000
  Mref_V%density_ref( 677 ) = 3441.91000000000
  Mref_V%vpv_ref( 677 ) = 8241.77000000000
  Mref_V%vsv_ref( 677 ) = 4516.65000000000
  Mref_V%Qkappa_ref( 677 ) = 943.000000000000
  Mref_V%Qmu_ref( 677 ) = 165.000000000000
  Mref_V%vph_ref( 677 ) = 8255.33000000000
  Mref_V%vsh_ref( 677 ) = 4580.21000000000
  Mref_V%eta_ref( 677 ) = 0.985710000000000

  Mref_V%radius_ref( 678 ) = 6151000.00000000
  Mref_V%density_ref( 678 ) = 3437.35000000000
  Mref_V%vpv_ref( 678 ) = 8207.37000000000
  Mref_V%vsv_ref( 678 ) = 4506.09000000000
  Mref_V%Qkappa_ref( 678 ) = 943.000000000000
  Mref_V%Qmu_ref( 678 ) = 165.000000000000
  Mref_V%vph_ref( 678 ) = 8236.90000000000
  Mref_V%vsh_ref( 678 ) = 4575.75000000000
  Mref_V%eta_ref( 678 ) = 0.982550000000000

  Mref_V%radius_ref( 679 ) = 6151000.00000000
  Mref_V%density_ref( 679 ) = 3437.34000000000
  Mref_V%vpv_ref( 679 ) = 8207.01000000000
  Mref_V%vsv_ref( 679 ) = 4506.00000000000
  Mref_V%Qkappa_ref( 679 ) = 943.000000000000
  Mref_V%Qmu_ref( 679 ) = 70.0000000000000
  Mref_V%vph_ref( 679 ) = 8236.81000000000
  Mref_V%vsh_ref( 679 ) = 4575.74000000000
  Mref_V%eta_ref( 679 ) = 0.982500000000000

  Mref_V%radius_ref( 680 ) = 6157087.00000000
  Mref_V%density_ref( 680 ) = 3432.81000000000
  Mref_V%vpv_ref( 680 ) = 8174.32000000000
  Mref_V%vsv_ref( 680 ) = 4496.29000000000
  Mref_V%Qkappa_ref( 680 ) = 943.000000000000
  Mref_V%Qmu_ref( 680 ) = 70.0000000000000
  Mref_V%vph_ref( 680 ) = 8222.27000000000
  Mref_V%vsh_ref( 680 ) = 4572.27000000000
  Mref_V%eta_ref( 680 ) = 0.979070000000000

  Mref_V%radius_ref( 681 ) = 6163174.00000000
  Mref_V%density_ref( 681 ) = 3428.15000000000
  Mref_V%vpv_ref( 681 ) = 8141.99000000000
  Mref_V%vsv_ref( 681 ) = 4487.00000000000
  Mref_V%Qkappa_ref( 681 ) = 943.000000000000
  Mref_V%Qmu_ref( 681 ) = 70.0000000000000
  Mref_V%vph_ref( 681 ) = 8210.47000000000
  Mref_V%vsh_ref( 681 ) = 4569.53000000000
  Mref_V%eta_ref( 681 ) = 0.975310000000000

  Mref_V%radius_ref( 682 ) = 6169261.00000000
  Mref_V%density_ref( 682 ) = 3423.37000000000
  Mref_V%vpv_ref( 682 ) = 8110.40000000000
  Mref_V%vsv_ref( 682 ) = 4478.17000000000
  Mref_V%Qkappa_ref( 682 ) = 943.000000000000
  Mref_V%Qmu_ref( 682 ) = 70.0000000000000
  Mref_V%vph_ref( 682 ) = 8201.14000000000
  Mref_V%vsh_ref( 682 ) = 4567.46000000000
  Mref_V%eta_ref( 682 ) = 0.971280000000000

  Mref_V%radius_ref( 683 ) = 6175348.00000000
  Mref_V%density_ref( 683 ) = 3418.47000000000
  Mref_V%vpv_ref( 683 ) = 8079.71000000000
  Mref_V%vsv_ref( 683 ) = 4469.83000000000
  Mref_V%Qkappa_ref( 683 ) = 943.000000000000
  Mref_V%Qmu_ref( 683 ) = 70.0000000000000
  Mref_V%vph_ref( 683 ) = 8193.99000000000
  Mref_V%vsh_ref( 683 ) = 4566.02000000000
  Mref_V%eta_ref( 683 ) = 0.967040000000000

  Mref_V%radius_ref( 684 ) = 6181435.00000000
  Mref_V%density_ref( 684 ) = 3413.47000000000
  Mref_V%vpv_ref( 684 ) = 8050.15000000000
  Mref_V%vsv_ref( 684 ) = 4462.00000000000
  Mref_V%Qkappa_ref( 684 ) = 943.000000000000
  Mref_V%Qmu_ref( 684 ) = 70.0000000000000
  Mref_V%vph_ref( 684 ) = 8188.67000000000
  Mref_V%vsh_ref( 684 ) = 4565.10000000000
  Mref_V%eta_ref( 684 ) = 0.962680000000000

  Mref_V%radius_ref( 685 ) = 6187522.00000000
  Mref_V%density_ref( 685 ) = 3408.35000000000
  Mref_V%vpv_ref( 685 ) = 8021.89000000000
  Mref_V%vsv_ref( 685 ) = 4454.69000000000
  Mref_V%Qkappa_ref( 685 ) = 943.000000000000
  Mref_V%Qmu_ref( 685 ) = 70.0000000000000
  Mref_V%vph_ref( 685 ) = 8184.92000000000
  Mref_V%vsh_ref( 685 ) = 4564.66000000000
  Mref_V%eta_ref( 685 ) = 0.958230000000000

  Mref_V%radius_ref( 686 ) = 6193609.00000000
  Mref_V%density_ref( 686 ) = 3403.15000000000
  Mref_V%vpv_ref( 686 ) = 7995.08000000000
  Mref_V%vsv_ref( 686 ) = 4447.94000000000
  Mref_V%Qkappa_ref( 686 ) = 943.000000000000
  Mref_V%Qmu_ref( 686 ) = 70.0000000000000
  Mref_V%vph_ref( 686 ) = 8182.39000000000
  Mref_V%vsh_ref( 686 ) = 4564.65000000000
  Mref_V%eta_ref( 686 ) = 0.953780000000000

  Mref_V%radius_ref( 687 ) = 6199696.00000000
  Mref_V%density_ref( 687 ) = 3397.87000000000
  Mref_V%vpv_ref( 687 ) = 7969.97000000000
  Mref_V%vsv_ref( 687 ) = 4441.76000000000
  Mref_V%Qkappa_ref( 687 ) = 943.000000000000
  Mref_V%Qmu_ref( 687 ) = 70.0000000000000
  Mref_V%vph_ref( 687 ) = 8180.79000000000
  Mref_V%vsh_ref( 687 ) = 4564.99000000000
  Mref_V%eta_ref( 687 ) = 0.949380000000000

  Mref_V%radius_ref( 688 ) = 6205783.00000000
  Mref_V%density_ref( 688 ) = 3392.50000000000
  Mref_V%vpv_ref( 688 ) = 7946.70000000000
  Mref_V%vsv_ref( 688 ) = 4436.18000000000
  Mref_V%Qkappa_ref( 688 ) = 943.000000000000
  Mref_V%Qmu_ref( 688 ) = 70.0000000000000
  Mref_V%vph_ref( 688 ) = 8179.83000000000
  Mref_V%vsh_ref( 688 ) = 4565.62000000000
  Mref_V%eta_ref( 688 ) = 0.945090000000000

  Mref_V%radius_ref( 689 ) = 6211870.00000000
  Mref_V%density_ref( 689 ) = 3387.07000000000
  Mref_V%vpv_ref( 689 ) = 7925.45000000000
  Mref_V%vsv_ref( 689 ) = 4431.20000000000
  Mref_V%Qkappa_ref( 689 ) = 943.000000000000
  Mref_V%Qmu_ref( 689 ) = 70.0000000000000
  Mref_V%vph_ref( 689 ) = 8179.17000000000
  Mref_V%vsh_ref( 689 ) = 4566.47000000000
  Mref_V%eta_ref( 689 ) = 0.940980000000000

  Mref_V%radius_ref( 690 ) = 6217957.00000000
  Mref_V%density_ref( 690 ) = 3381.58000000000
  Mref_V%vpv_ref( 690 ) = 7906.44000000000
  Mref_V%vsv_ref( 690 ) = 4426.83000000000
  Mref_V%Qkappa_ref( 690 ) = 943.000000000000
  Mref_V%Qmu_ref( 690 ) = 70.0000000000000
  Mref_V%vph_ref( 690 ) = 8178.54000000000
  Mref_V%vsh_ref( 690 ) = 4567.46000000000
  Mref_V%eta_ref( 690 ) = 0.937120000000000

  Mref_V%radius_ref( 691 ) = 6224043.00000000
  Mref_V%density_ref( 691 ) = 3376.03000000000
  Mref_V%vpv_ref( 691 ) = 7889.80000000000
  Mref_V%vsv_ref( 691 ) = 4423.12000000000
  Mref_V%Qkappa_ref( 691 ) = 943.000000000000
  Mref_V%Qmu_ref( 691 ) = 70.0000000000000
  Mref_V%vph_ref( 691 ) = 8177.64000000000
  Mref_V%vsh_ref( 691 ) = 4568.58000000000
  Mref_V%eta_ref( 691 ) = 0.933560000000000

  Mref_V%radius_ref( 692 ) = 6230130.00000000
  Mref_V%density_ref( 692 ) = 3370.45000000000
  Mref_V%vpv_ref( 692 ) = 7875.56000000000
  Mref_V%vsv_ref( 692 ) = 4420.09000000000
  Mref_V%Qkappa_ref( 692 ) = 943.000000000000
  Mref_V%Qmu_ref( 692 ) = 70.0000000000000
  Mref_V%vph_ref( 692 ) = 8176.30000000000
  Mref_V%vsh_ref( 692 ) = 4569.70000000000
  Mref_V%eta_ref( 692 ) = 0.930340000000000

  Mref_V%radius_ref( 693 ) = 6236217.00000000
  Mref_V%density_ref( 693 ) = 3364.87000000000
  Mref_V%vpv_ref( 693 ) = 7863.64000000000
  Mref_V%vsv_ref( 693 ) = 4417.81000000000
  Mref_V%Qkappa_ref( 693 ) = 943.000000000000
  Mref_V%Qmu_ref( 693 ) = 70.0000000000000
  Mref_V%vph_ref( 693 ) = 8174.55000000000
  Mref_V%vsh_ref( 693 ) = 4570.85000000000
  Mref_V%eta_ref( 693 ) = 0.927430000000000

  Mref_V%radius_ref( 694 ) = 6242304.00000000
  Mref_V%density_ref( 694 ) = 3359.31000000000
  Mref_V%vpv_ref( 694 ) = 7853.87000000000
  Mref_V%vsv_ref( 694 ) = 4416.30000000000
  Mref_V%Qkappa_ref( 694 ) = 943.000000000000
  Mref_V%Qmu_ref( 694 ) = 70.0000000000000
  Mref_V%vph_ref( 694 ) = 8172.42000000000
  Mref_V%vsh_ref( 694 ) = 4571.91000000000
  Mref_V%eta_ref( 694 ) = 0.924830000000000

  Mref_V%radius_ref( 695 ) = 6248391.00000000
  Mref_V%density_ref( 695 ) = 3353.79000000000
  Mref_V%vpv_ref( 695 ) = 7846.17000000000
  Mref_V%vsv_ref( 695 ) = 4415.67000000000
  Mref_V%Qkappa_ref( 695 ) = 943.000000000000
  Mref_V%Qmu_ref( 695 ) = 70.0000000000000
  Mref_V%vph_ref( 695 ) = 8169.91000000000
  Mref_V%vsh_ref( 695 ) = 4572.83000000000
  Mref_V%eta_ref( 695 ) = 0.922510000000000

  Mref_V%radius_ref( 696 ) = 6254478.00000000
  Mref_V%density_ref( 696 ) = 3348.37000000000
  Mref_V%vpv_ref( 696 ) = 7840.38000000000
  Mref_V%vsv_ref( 696 ) = 4415.93000000000
  Mref_V%Qkappa_ref( 696 ) = 943.000000000000
  Mref_V%Qmu_ref( 696 ) = 70.0000000000000
  Mref_V%vph_ref( 696 ) = 8167.05000000000
  Mref_V%vsh_ref( 696 ) = 4573.60000000000
  Mref_V%eta_ref( 696 ) = 0.920460000000000

  Mref_V%radius_ref( 697 ) = 6260565.00000000
  Mref_V%density_ref( 697 ) = 3343.03000000000
  Mref_V%vpv_ref( 697 ) = 7836.39000000000
  Mref_V%vsv_ref( 697 ) = 4417.15000000000
  Mref_V%Qkappa_ref( 697 ) = 943.000000000000
  Mref_V%Qmu_ref( 697 ) = 70.0000000000000
  Mref_V%vph_ref( 697 ) = 8163.88000000000
  Mref_V%vsh_ref( 697 ) = 4574.16000000000
  Mref_V%eta_ref( 697 ) = 0.918670000000000

  Mref_V%radius_ref( 698 ) = 6266652.00000000
  Mref_V%density_ref( 698 ) = 3337.85000000000
  Mref_V%vpv_ref( 698 ) = 7834.11000000000
  Mref_V%vsv_ref( 698 ) = 4419.42000000000
  Mref_V%Qkappa_ref( 698 ) = 943.000000000000
  Mref_V%Qmu_ref( 698 ) = 70.0000000000000
  Mref_V%vph_ref( 698 ) = 8160.37000000000
  Mref_V%vsh_ref( 698 ) = 4574.44000000000
  Mref_V%eta_ref( 698 ) = 0.917110000000000

  Mref_V%radius_ref( 699 ) = 6272739.00000000
  Mref_V%density_ref( 699 ) = 3332.83000000000
  Mref_V%vpv_ref( 699 ) = 7833.38000000000
  Mref_V%vsv_ref( 699 ) = 4422.78000000000
  Mref_V%Qkappa_ref( 699 ) = 943.000000000000
  Mref_V%Qmu_ref( 699 ) = 70.0000000000000
  Mref_V%vph_ref( 699 ) = 8156.58000000000
  Mref_V%vsh_ref( 699 ) = 4574.42000000000
  Mref_V%eta_ref( 699 ) = 0.915770000000000

  Mref_V%radius_ref( 700 ) = 6278826.00000000
  Mref_V%density_ref( 700 ) = 3328.01000000000
  Mref_V%vpv_ref( 700 ) = 7834.11000000000
  Mref_V%vsv_ref( 700 ) = 4427.25000000000
  Mref_V%Qkappa_ref( 700 ) = 943.000000000000
  Mref_V%Qmu_ref( 700 ) = 70.0000000000000
  Mref_V%vph_ref( 700 ) = 8152.57000000000
  Mref_V%vsh_ref( 700 ) = 4574.04000000000
  Mref_V%eta_ref( 700 ) = 0.914650000000000

  Mref_V%radius_ref( 701 ) = 6284913.00000000
  Mref_V%density_ref( 701 ) = 3323.39000000000
  Mref_V%vpv_ref( 701 ) = 7836.11000000000
  Mref_V%vsv_ref( 701 ) = 4432.88000000000
  Mref_V%Qkappa_ref( 701 ) = 943.000000000000
  Mref_V%Qmu_ref( 701 ) = 70.0000000000000
  Mref_V%vph_ref( 701 ) = 8148.41000000000
  Mref_V%vsh_ref( 701 ) = 4573.36000000000
  Mref_V%eta_ref( 701 ) = 0.913710000000000

  Mref_V%radius_ref( 702 ) = 6291000.00000000
  Mref_V%density_ref( 702 ) = 3319.04000000000
  Mref_V%vpv_ref( 702 ) = 7839.12000000000
  Mref_V%vsv_ref( 702 ) = 4439.57000000000
  Mref_V%Qkappa_ref( 702 ) = 943.000000000000
  Mref_V%Qmu_ref( 702 ) = 70.0000000000000
  Mref_V%vph_ref( 702 ) = 8144.20000000000
  Mref_V%vsh_ref( 702 ) = 4572.41000000000
  Mref_V%eta_ref( 702 ) = 0.912960000000000

  Mref_V%radius_ref( 703 ) = 6291000.00000000
  Mref_V%density_ref( 703 ) = 3319.03000000000
  Mref_V%vpv_ref( 703 ) = 7839.37000000000
  Mref_V%vsv_ref( 703 ) = 4439.74000000000
  Mref_V%Qkappa_ref( 703 ) = 943.000000000000
  Mref_V%Qmu_ref( 703 ) = 191.000000000000
  Mref_V%vph_ref( 703 ) = 8144.32000000000
  Mref_V%vsh_ref( 703 ) = 4572.46000000000
  Mref_V%eta_ref( 703 ) = 0.912940000000000

  Mref_V%radius_ref( 704 ) = 6294971.00000000
  Mref_V%density_ref( 704 ) = 3316.33000000000
  Mref_V%vpv_ref( 704 ) = 7841.82000000000
  Mref_V%vsv_ref( 704 ) = 4444.71000000000
  Mref_V%Qkappa_ref( 704 ) = 943.000000000000
  Mref_V%Qmu_ref( 704 ) = 191.000000000000
  Mref_V%vph_ref( 704 ) = 8141.60000000000
  Mref_V%vsh_ref( 704 ) = 4571.71000000000
  Mref_V%eta_ref( 704 ) = 0.912540000000000

  Mref_V%radius_ref( 705 ) = 6298943.00000000
  Mref_V%density_ref( 705 ) = 3313.75000000000
  Mref_V%vpv_ref( 705 ) = 7844.77000000000
  Mref_V%vsv_ref( 705 ) = 4450.28000000000
  Mref_V%Qkappa_ref( 705 ) = 943.000000000000
  Mref_V%Qmu_ref( 705 ) = 191.000000000000
  Mref_V%vph_ref( 705 ) = 8139.01000000000
  Mref_V%vsh_ref( 705 ) = 4570.93000000000
  Mref_V%eta_ref( 705 ) = 0.912210000000000

  Mref_V%radius_ref( 706 ) = 6302914.00000000
  Mref_V%density_ref( 706 ) = 3311.30000000000
  Mref_V%vpv_ref( 706 ) = 7848.07000000000
  Mref_V%vsv_ref( 706 ) = 4456.35000000000
  Mref_V%Qkappa_ref( 706 ) = 943.000000000000
  Mref_V%Qmu_ref( 706 ) = 191.000000000000
  Mref_V%vph_ref( 706 ) = 8136.50000000000
  Mref_V%vsh_ref( 706 ) = 4570.06000000000
  Mref_V%eta_ref( 706 ) = 0.911930000000000

  Mref_V%radius_ref( 707 ) = 6306886.00000000
  Mref_V%density_ref( 707 ) = 3308.97000000000
  Mref_V%vpv_ref( 707 ) = 7851.72000000000
  Mref_V%vsv_ref( 707 ) = 4462.89000000000
  Mref_V%Qkappa_ref( 707 ) = 943.000000000000
  Mref_V%Qmu_ref( 707 ) = 191.000000000000
  Mref_V%vph_ref( 707 ) = 8134.11000000000
  Mref_V%vsh_ref( 707 ) = 4569.16000000000
  Mref_V%eta_ref( 707 ) = 0.911710000000000

  Mref_V%radius_ref( 708 ) = 6310857.00000000
  Mref_V%density_ref( 708 ) = 3306.73000000000
  Mref_V%vpv_ref( 708 ) = 7855.75000000000
  Mref_V%vsv_ref( 708 ) = 4469.94000000000
  Mref_V%Qkappa_ref( 708 ) = 943.000000000000
  Mref_V%Qmu_ref( 708 ) = 191.000000000000
  Mref_V%vph_ref( 708 ) = 8131.82000000000
  Mref_V%vsh_ref( 708 ) = 4568.21000000000
  Mref_V%eta_ref( 708 ) = 0.911550000000000

  Mref_V%radius_ref( 709 ) = 6314829.00000000
  Mref_V%density_ref( 709 ) = 3304.58000000000
  Mref_V%vpv_ref( 709 ) = 7860.14000000000
  Mref_V%vsv_ref( 709 ) = 4477.40000000000
  Mref_V%Qkappa_ref( 709 ) = 943.000000000000
  Mref_V%Qmu_ref( 709 ) = 191.000000000000
  Mref_V%vph_ref( 709 ) = 8129.66000000000
  Mref_V%vsh_ref( 709 ) = 4567.22000000000
  Mref_V%eta_ref( 709 ) = 0.911420000000000

  Mref_V%radius_ref( 710 ) = 6318800.00000000
  Mref_V%density_ref( 710 ) = 3302.53000000000
  Mref_V%vpv_ref( 710 ) = 7864.89000000000
  Mref_V%vsv_ref( 710 ) = 4485.33000000000
  Mref_V%Qkappa_ref( 710 ) = 943.000000000000
  Mref_V%Qmu_ref( 710 ) = 191.000000000000
  Mref_V%vph_ref( 710 ) = 8127.60000000000
  Mref_V%vsh_ref( 710 ) = 4566.21000000000
  Mref_V%eta_ref( 710 ) = 0.911340000000000

  Mref_V%radius_ref( 711 ) = 6322771.00000000
  Mref_V%density_ref( 711 ) = 3300.55000000000
  Mref_V%vpv_ref( 711 ) = 7870.01000000000
  Mref_V%vsv_ref( 711 ) = 4493.69000000000
  Mref_V%Qkappa_ref( 711 ) = 943.000000000000
  Mref_V%Qmu_ref( 711 ) = 191.000000000000
  Mref_V%vph_ref( 711 ) = 8125.65000000000
  Mref_V%vsh_ref( 711 ) = 4565.16000000000
  Mref_V%eta_ref( 711 ) = 0.911300000000000

  Mref_V%radius_ref( 712 ) = 6326743.00000000
  Mref_V%density_ref( 712 ) = 3298.63000000000
  Mref_V%vpv_ref( 712 ) = 7875.49000000000
  Mref_V%vsv_ref( 712 ) = 4502.48000000000
  Mref_V%Qkappa_ref( 712 ) = 943.000000000000
  Mref_V%Qmu_ref( 712 ) = 191.000000000000
  Mref_V%vph_ref( 712 ) = 8123.87000000000
  Mref_V%vsh_ref( 712 ) = 4564.11000000000
  Mref_V%eta_ref( 712 ) = 0.911290000000000

  Mref_V%radius_ref( 713 ) = 6330714.00000000
  Mref_V%density_ref( 713 ) = 3296.79000000000
  Mref_V%vpv_ref( 713 ) = 7881.33000000000
  Mref_V%vsv_ref( 713 ) = 4511.66000000000
  Mref_V%Qkappa_ref( 713 ) = 943.000000000000
  Mref_V%Qmu_ref( 713 ) = 191.000000000000
  Mref_V%vph_ref( 713 ) = 8122.23000000000
  Mref_V%vsh_ref( 713 ) = 4563.05000000000
  Mref_V%eta_ref( 713 ) = 0.911300000000000

  Mref_V%radius_ref( 714 ) = 6334686.00000000
  Mref_V%density_ref( 714 ) = 3295.00000000000
  Mref_V%vpv_ref( 714 ) = 7887.54000000000
  Mref_V%vsv_ref( 714 ) = 4521.24000000000
  Mref_V%Qkappa_ref( 714 ) = 943.000000000000
  Mref_V%Qmu_ref( 714 ) = 191.000000000000
  Mref_V%vph_ref( 714 ) = 8120.74000000000
  Mref_V%vsh_ref( 714 ) = 4562.00000000000
  Mref_V%eta_ref( 714 ) = 0.911350000000000

  Mref_V%radius_ref( 715 ) = 6338657.00000000
  Mref_V%density_ref( 715 ) = 3293.25000000000
  Mref_V%vpv_ref( 715 ) = 7894.13000000000
  Mref_V%vsv_ref( 715 ) = 4531.23000000000
  Mref_V%Qkappa_ref( 715 ) = 943.000000000000
  Mref_V%Qmu_ref( 715 ) = 191.000000000000
  Mref_V%vph_ref( 715 ) = 8119.38000000000
  Mref_V%vsh_ref( 715 ) = 4560.94000000000
  Mref_V%eta_ref( 715 ) = 0.911400000000000

  Mref_V%radius_ref( 716 ) = 6342629.00000000
  Mref_V%density_ref( 716 ) = 3291.54000000000
  Mref_V%vpv_ref( 716 ) = 7901.10000000000
  Mref_V%vsv_ref( 716 ) = 4541.57000000000
  Mref_V%Qkappa_ref( 716 ) = 943.000000000000
  Mref_V%Qmu_ref( 716 ) = 191.000000000000
  Mref_V%vph_ref( 716 ) = 8118.22000000000
  Mref_V%vsh_ref( 716 ) = 4559.94000000000
  Mref_V%eta_ref( 716 ) = 0.911470000000000

  Mref_V%radius_ref( 717 ) = 6346600.00000000
  Mref_V%density_ref( 717 ) = 3289.84000000000
  Mref_V%vpv_ref( 717 ) = 7908.24000000000
  Mref_V%vsv_ref( 717 ) = 4552.08000000000
  Mref_V%Qkappa_ref( 717 ) = 943.000000000000
  Mref_V%Qmu_ref( 717 ) = 191.000000000000
  Mref_V%vph_ref( 717 ) = 8117.13000000000
  Mref_V%vsh_ref( 717 ) = 4558.94000000000
  Mref_V%eta_ref( 717 ) = 0.911550000000000

  Mref_V%radius_ref( 718 ) = 6346600.00000000
  Mref_V%density_ref( 718 ) = 2900.00000000000
  Mref_V%vpv_ref( 718 ) = 6800.00000000000
  Mref_V%vsv_ref( 718 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 718 ) = 57822.5000000000
  Mref_V%Qmu_ref( 718 ) = 300.000000000000
  Mref_V%vph_ref( 718 ) = 6800.00000000000
  Mref_V%vsh_ref( 718 ) = 3900.00000000000
  Mref_V%eta_ref( 718 ) = 1.00000000000000

  Mref_V%radius_ref( 719 ) = 6347540.00000000
  Mref_V%density_ref( 719 ) = 2900.00000000000
  Mref_V%vpv_ref( 719 ) = 6800.00000000000
  Mref_V%vsv_ref( 719 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 719 ) = 57822.5000000000
  Mref_V%Qmu_ref( 719 ) = 300.000000000000
  Mref_V%vph_ref( 719 ) = 6800.00000000000
  Mref_V%vsh_ref( 719 ) = 3900.00000000000
  Mref_V%eta_ref( 719 ) = 1.00000000000000

  Mref_V%radius_ref( 720 ) = 6348480.00000000
  Mref_V%density_ref( 720 ) = 2900.00000000000
  Mref_V%vpv_ref( 720 ) = 6800.00000000000
  Mref_V%vsv_ref( 720 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 720 ) = 57822.5000000000
  Mref_V%Qmu_ref( 720 ) = 300.000000000000
  Mref_V%vph_ref( 720 ) = 6800.00000000000
  Mref_V%vsh_ref( 720 ) = 3900.00000000000
  Mref_V%eta_ref( 720 ) = 1.00000000000000

  Mref_V%radius_ref( 721 ) = 6349420.00000000
  Mref_V%density_ref( 721 ) = 2900.00000000000
  Mref_V%vpv_ref( 721 ) = 6800.00000000000
  Mref_V%vsv_ref( 721 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 721 ) = 57822.5000000000
  Mref_V%Qmu_ref( 721 ) = 300.000000000000
  Mref_V%vph_ref( 721 ) = 6800.00000000000
  Mref_V%vsh_ref( 721 ) = 3900.00000000000
  Mref_V%eta_ref( 721 ) = 1.00000000000000

  Mref_V%radius_ref( 722 ) = 6350360.00000000
  Mref_V%density_ref( 722 ) = 2900.00000000000
  Mref_V%vpv_ref( 722 ) = 6800.00000000000
  Mref_V%vsv_ref( 722 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 722 ) = 57822.5000000000
  Mref_V%Qmu_ref( 722 ) = 300.000000000000
  Mref_V%vph_ref( 722 ) = 6800.00000000000
  Mref_V%vsh_ref( 722 ) = 3900.00000000000
  Mref_V%eta_ref( 722 ) = 1.00000000000000

  Mref_V%radius_ref( 723 ) = 6351300.00000000
  Mref_V%density_ref( 723 ) = 2900.00000000000
  Mref_V%vpv_ref( 723 ) = 6800.00000000000
  Mref_V%vsv_ref( 723 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 723 ) = 57822.5000000000
  Mref_V%Qmu_ref( 723 ) = 300.000000000000
  Mref_V%vph_ref( 723 ) = 6800.00000000000
  Mref_V%vsh_ref( 723 ) = 3900.00000000000
  Mref_V%eta_ref( 723 ) = 1.00000000000000

  Mref_V%radius_ref( 724 ) = 6352240.00000000
  Mref_V%density_ref( 724 ) = 2900.00000000000
  Mref_V%vpv_ref( 724 ) = 6800.00000000000
  Mref_V%vsv_ref( 724 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 724 ) = 57822.5000000000
  Mref_V%Qmu_ref( 724 ) = 300.000000000000
  Mref_V%vph_ref( 724 ) = 6800.00000000000
  Mref_V%vsh_ref( 724 ) = 3900.00000000000
  Mref_V%eta_ref( 724 ) = 1.00000000000000

  Mref_V%radius_ref( 725 ) = 6353180.00000000
  Mref_V%density_ref( 725 ) = 2900.00000000000
  Mref_V%vpv_ref( 725 ) = 6800.00000000000
  Mref_V%vsv_ref( 725 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 725 ) = 57822.5000000000
  Mref_V%Qmu_ref( 725 ) = 300.000000000000
  Mref_V%vph_ref( 725 ) = 6800.00000000000
  Mref_V%vsh_ref( 725 ) = 3900.00000000000
  Mref_V%eta_ref( 725 ) = 1.00000000000000

  Mref_V%radius_ref( 726 ) = 6354120.00000000
  Mref_V%density_ref( 726 ) = 2900.00000000000
  Mref_V%vpv_ref( 726 ) = 6800.00000000000
  Mref_V%vsv_ref( 726 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 726 ) = 57822.5000000000
  Mref_V%Qmu_ref( 726 ) = 300.000000000000
  Mref_V%vph_ref( 726 ) = 6800.00000000000
  Mref_V%vsh_ref( 726 ) = 3900.00000000000
  Mref_V%eta_ref( 726 ) = 1.00000000000000

  Mref_V%radius_ref( 727 ) = 6355060.00000000
  Mref_V%density_ref( 727 ) = 2900.00000000000
  Mref_V%vpv_ref( 727 ) = 6800.00000000000
  Mref_V%vsv_ref( 727 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 727 ) = 57822.5000000000
  Mref_V%Qmu_ref( 727 ) = 300.000000000000
  Mref_V%vph_ref( 727 ) = 6800.00000000000
  Mref_V%vsh_ref( 727 ) = 3900.00000000000
  Mref_V%eta_ref( 727 ) = 1.00000000000000

  Mref_V%radius_ref( 728 ) = 6356000.00000000
  Mref_V%density_ref( 728 ) = 2900.00000000000
  Mref_V%vpv_ref( 728 ) = 6800.00000000000
  Mref_V%vsv_ref( 728 ) = 3900.00000000000
  Mref_V%Qkappa_ref( 728 ) = 57822.5000000000
  Mref_V%Qmu_ref( 728 ) = 300.000000000000
  Mref_V%vph_ref( 728 ) = 6800.00000000000
  Mref_V%vsh_ref( 728 ) = 3900.00000000000
  Mref_V%eta_ref( 728 ) = 1.00000000000000

  Mref_V%radius_ref( 729 ) = 6356000.00000000
  Mref_V%density_ref( 729 ) = 2600.00000000000
  Mref_V%vpv_ref( 729 ) = 5800.00000000000
  Mref_V%vsv_ref( 729 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 729 ) = 57822.5000000000
  Mref_V%Qmu_ref( 729 ) = 300.000000000000
  Mref_V%vph_ref( 729 ) = 5800.00000000000
  Mref_V%vsh_ref( 729 ) = 3200.00000000000
  Mref_V%eta_ref( 729 ) = 1.00000000000000

  Mref_V%radius_ref( 730 ) = 6357200.00000000
  Mref_V%density_ref( 730 ) = 2600.00000000000
  Mref_V%vpv_ref( 730 ) = 5800.00000000000
  Mref_V%vsv_ref( 730 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 730 ) = 57822.5000000000
  Mref_V%Qmu_ref( 730 ) = 300.000000000000
  Mref_V%vph_ref( 730 ) = 5800.00000000000
  Mref_V%vsh_ref( 730 ) = 3200.00000000000
  Mref_V%eta_ref( 730 ) = 1.00000000000000

  Mref_V%radius_ref( 731 ) = 6358400.00000000
  Mref_V%density_ref( 731 ) = 2600.00000000000
  Mref_V%vpv_ref( 731 ) = 5800.00000000000
  Mref_V%vsv_ref( 731 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 731 ) = 57822.5000000000
  Mref_V%Qmu_ref( 731 ) = 300.000000000000
  Mref_V%vph_ref( 731 ) = 5800.00000000000
  Mref_V%vsh_ref( 731 ) = 3200.00000000000
  Mref_V%eta_ref( 731 ) = 1.00000000000000

  Mref_V%radius_ref( 732 ) = 6359600.00000000
  Mref_V%density_ref( 732 ) = 2600.00000000000
  Mref_V%vpv_ref( 732 ) = 5800.00000000000
  Mref_V%vsv_ref( 732 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 732 ) = 57822.5000000000
  Mref_V%Qmu_ref( 732 ) = 300.000000000000
  Mref_V%vph_ref( 732 ) = 5800.00000000000
  Mref_V%vsh_ref( 732 ) = 3200.00000000000
  Mref_V%eta_ref( 732 ) = 1.00000000000000

  Mref_V%radius_ref( 733 ) = 6360800.00000000
  Mref_V%density_ref( 733 ) = 2600.00000000000
  Mref_V%vpv_ref( 733 ) = 5800.00000000000
  Mref_V%vsv_ref( 733 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 733 ) = 57822.5000000000
  Mref_V%Qmu_ref( 733 ) = 300.000000000000
  Mref_V%vph_ref( 733 ) = 5800.00000000000
  Mref_V%vsh_ref( 733 ) = 3200.00000000000
  Mref_V%eta_ref( 733 ) = 1.00000000000000

  Mref_V%radius_ref( 734 ) = 6362000.00000000
  Mref_V%density_ref( 734 ) = 2600.00000000000
  Mref_V%vpv_ref( 734 ) = 5800.00000000000
  Mref_V%vsv_ref( 734 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 734 ) = 57822.5000000000
  Mref_V%Qmu_ref( 734 ) = 300.000000000000
  Mref_V%vph_ref( 734 ) = 5800.00000000000
  Mref_V%vsh_ref( 734 ) = 3200.00000000000
  Mref_V%eta_ref( 734 ) = 1.00000000000000

  Mref_V%radius_ref( 735 ) = 6363200.00000000
  Mref_V%density_ref( 735 ) = 2600.00000000000
  Mref_V%vpv_ref( 735 ) = 5800.00000000000
  Mref_V%vsv_ref( 735 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 735 ) = 57822.5000000000
  Mref_V%Qmu_ref( 735 ) = 300.000000000000
  Mref_V%vph_ref( 735 ) = 5800.00000000000
  Mref_V%vsh_ref( 735 ) = 3200.00000000000
  Mref_V%eta_ref( 735 ) = 1.00000000000000

  Mref_V%radius_ref( 736 ) = 6364400.00000000
  Mref_V%density_ref( 736 ) = 2600.00000000000
  Mref_V%vpv_ref( 736 ) = 5800.00000000000
  Mref_V%vsv_ref( 736 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 736 ) = 57822.5000000000
  Mref_V%Qmu_ref( 736 ) = 300.000000000000
  Mref_V%vph_ref( 736 ) = 5800.00000000000
  Mref_V%vsh_ref( 736 ) = 3200.00000000000
  Mref_V%eta_ref( 736 ) = 1.00000000000000

  Mref_V%radius_ref( 737 ) = 6365600.00000000
  Mref_V%density_ref( 737 ) = 2600.00000000000
  Mref_V%vpv_ref( 737 ) = 5800.00000000000
  Mref_V%vsv_ref( 737 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 737 ) = 57822.5000000000
  Mref_V%Qmu_ref( 737 ) = 300.000000000000
  Mref_V%vph_ref( 737 ) = 5800.00000000000
  Mref_V%vsh_ref( 737 ) = 3200.00000000000
  Mref_V%eta_ref( 737 ) = 1.00000000000000

  Mref_V%radius_ref( 738 ) = 6366800.00000000
  Mref_V%density_ref( 738 ) = 2600.00000000000
  Mref_V%vpv_ref( 738 ) = 5800.00000000000
  Mref_V%vsv_ref( 738 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 738 ) = 57822.5000000000
  Mref_V%Qmu_ref( 738 ) = 300.000000000000
  Mref_V%vph_ref( 738 ) = 5800.00000000000
  Mref_V%vsh_ref( 738 ) = 3200.00000000000
  Mref_V%eta_ref( 738 ) = 1.00000000000000

  Mref_V%radius_ref( 739 ) = 6368000.00000000
  Mref_V%density_ref( 739 ) = 2600.00000000000
  Mref_V%vpv_ref( 739 ) = 5800.00000000000
  Mref_V%vsv_ref( 739 ) = 3200.00000000000
  Mref_V%Qkappa_ref( 739 ) = 57822.5000000000
  Mref_V%Qmu_ref( 739 ) = 300.000000000000
  Mref_V%vph_ref( 739 ) = 5800.00000000000
  Mref_V%vsh_ref( 739 ) = 3200.00000000000
  Mref_V%eta_ref( 739 ) = 1.00000000000000

  Mref_V%radius_ref( 740 ) = 6368000.00000000
  Mref_V%density_ref( 740 ) = 1020.00000000000
  Mref_V%vpv_ref( 740 ) = 1450.00000000000
  Mref_V%vsv_ref( 740 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 740 ) = 57822.5000000000
  Mref_V%Qmu_ref( 740 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 740 ) = 1450.00000000000
  Mref_V%vsh_ref( 740 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 740 ) = 1.00000000000000

  Mref_V%radius_ref( 741 ) = 6368300.00000000
  Mref_V%density_ref( 741 ) = 1020.00000000000
  Mref_V%vpv_ref( 741 ) = 1450.00000000000
  Mref_V%vsv_ref( 741 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 741 ) = 57822.5000000000
  Mref_V%Qmu_ref( 741 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 741 ) = 1450.00000000000
  Mref_V%vsh_ref( 741 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 741 ) = 1.00000000000000

  Mref_V%radius_ref( 742 ) = 6368600.00000000
  Mref_V%density_ref( 742 ) = 1020.00000000000
  Mref_V%vpv_ref( 742 ) = 1450.00000000000
  Mref_V%vsv_ref( 742 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 742 ) = 57822.5000000000
  Mref_V%Qmu_ref( 742 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 742 ) = 1450.00000000000
  Mref_V%vsh_ref( 742 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 742 ) = 1.00000000000000

  Mref_V%radius_ref( 743 ) = 6368900.00000000
  Mref_V%density_ref( 743 ) = 1020.00000000000
  Mref_V%vpv_ref( 743 ) = 1450.00000000000
  Mref_V%vsv_ref( 743 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 743 ) = 57822.5000000000
  Mref_V%Qmu_ref( 743 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 743 ) = 1450.00000000000
  Mref_V%vsh_ref( 743 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 743 ) = 1.00000000000000

  Mref_V%radius_ref( 744 ) = 6369200.00000000
  Mref_V%density_ref( 744 ) = 1020.00000000000
  Mref_V%vpv_ref( 744 ) = 1450.00000000000
  Mref_V%vsv_ref( 744 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 744 ) = 57822.5000000000
  Mref_V%Qmu_ref( 744 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 744 ) = 1450.00000000000
  Mref_V%vsh_ref( 744 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 744 ) = 1.00000000000000

  Mref_V%radius_ref( 745 ) = 6369500.00000000
  Mref_V%density_ref( 745 ) = 1020.00000000000
  Mref_V%vpv_ref( 745 ) = 1450.00000000000
  Mref_V%vsv_ref( 745 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 745 ) = 57822.5000000000
  Mref_V%Qmu_ref( 745 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 745 ) = 1450.00000000000
  Mref_V%vsh_ref( 745 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 745 ) = 1.00000000000000

  Mref_V%radius_ref( 746 ) = 6369800.00000000
  Mref_V%density_ref( 746 ) = 1020.00000000000
  Mref_V%vpv_ref( 746 ) = 1450.00000000000
  Mref_V%vsv_ref( 746 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 746 ) = 57822.5000000000
  Mref_V%Qmu_ref( 746 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 746 ) = 1450.00000000000
  Mref_V%vsh_ref( 746 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 746 ) = 1.00000000000000

  Mref_V%radius_ref( 747 ) = 6370100.00000000
  Mref_V%density_ref( 747 ) = 1020.00000000000
  Mref_V%vpv_ref( 747 ) = 1450.00000000000
  Mref_V%vsv_ref( 747 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 747 ) = 57822.5000000000
  Mref_V%Qmu_ref( 747 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 747 ) = 1450.00000000000
  Mref_V%vsh_ref( 747 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 747 ) = 1.00000000000000

  Mref_V%radius_ref( 748 ) = 6370400.00000000
  Mref_V%density_ref( 748 ) = 1020.00000000000
  Mref_V%vpv_ref( 748 ) = 1450.00000000000
  Mref_V%vsv_ref( 748 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 748 ) = 57822.5000000000
  Mref_V%Qmu_ref( 748 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 748 ) = 1450.00000000000
  Mref_V%vsh_ref( 748 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 748 ) = 1.00000000000000

  Mref_V%radius_ref( 749 ) = 6370700.00000000
  Mref_V%density_ref( 749 ) = 1020.00000000000
  Mref_V%vpv_ref( 749 ) = 1450.00000000000
  Mref_V%vsv_ref( 749 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 749 ) = 57822.5000000000
  Mref_V%Qmu_ref( 749 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 749 ) = 1450.00000000000
  Mref_V%vsh_ref( 749 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 749 ) = 1.00000000000000

  Mref_V%radius_ref( 750 ) = 6371000.00000000
  Mref_V%density_ref( 750 ) = 1020.00000000000
  Mref_V%vpv_ref( 750 ) = 1450.00000000000
  Mref_V%vsv_ref( 750 ) = 0.000000000000000E+000
  Mref_V%Qkappa_ref( 750 ) = 57822.5000000000
  Mref_V%Qmu_ref( 750 ) = 0.000000000000000E+000
  Mref_V%vph_ref( 750 ) = 1450.00000000000
  Mref_V%vsh_ref( 750 ) = 0.000000000000000E+000
  Mref_V%eta_ref( 750 ) = 1.00000000000000

! strip the crust and replace it by mantle if we use an external crustal model
  if(USE_EXTERNAL_CRUSTAL_MODEL) then
    do i=NR_REF-10,NR_REF
      Mref_V%density_ref(i) = Mref_V%density_ref(NR_REF-11)
      Mref_V%vpv_ref(i) = Mref_V%vpv_ref(NR_REF-11)
      Mref_V%vph_ref(i) = Mref_V%vph_ref(NR_REF-11)
      Mref_V%vsv_ref(i) = Mref_V%vsv_ref(NR_REF-11)
      Mref_V%vsh_ref(i) = Mref_V%vsh_ref(NR_REF-11)
      Mref_V%eta_ref(i) = Mref_V%eta_ref(NR_REF-11)
      Mref_V%Qkappa_ref(i) = Mref_V%Qkappa_ref(NR_REF-11)
      Mref_V%Qmu_ref(i) = Mref_V%Qmu_ref(NR_REF-11)
    enddo
  endif

  end subroutine define_model_ref

