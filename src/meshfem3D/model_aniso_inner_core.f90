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
! based on scaling factors by Ishii et al. (2002)
!
! one should add an MPI_BCAST in meshfem3D_models.f90 if one
! adds a 3D model or a read_aniso_inner_core_model subroutine
!--------------------------------------------------------------------------------------------------

  subroutine model_aniso_inner_core(x,c11,c33,c12,c13,c44,REFERENCE_1D_MODEL, &
                                    vpv,vph,vsv,vsh,rho,eta_aniso)

  use constants

  implicit none

! given a normalized radius x, gives non-dimensionalized c11,c33,c12,c13,c44

  integer REFERENCE_1D_MODEL

  double precision x,c11,c33,c12,c13,c44
  double precision rho,vpv,vph,vsv,vsh,eta_aniso

  ! local parameters
  double precision vp,vs
  double precision vpc,vsc,rhoc
  double precision vp0,vs0,rho0,A0
  double precision c66
  double precision scale_fac

  ! calculates isotropic values from given (transversely isotropic) reference values
  ! (are non-dimensionalized)
  vp = sqrt(((8.d0+4.d0*eta_aniso)*vph*vph + 3.d0*vpv*vpv &
            + (8.d0 - 8.d0*eta_aniso)*vsv*vsv)/15.d0)
  vs = sqrt(((1.d0-2.d0*eta_aniso)*vph*vph + vpv*vpv &
                    + 5.d0*vsh*vsh + (6.d0+4.d0*eta_aniso)*vsv*vsv)/15.d0)

  ! scale to dimensions (e.g. used in prem model)
  scale_fac = R_EARTH*dsqrt(PI*GRAV*RHOAV)/1000.d0
  vp = vp * scale_fac
  vs = vs * scale_fac
  rho = rho * RHOAV/1000.d0

  select case (REFERENCE_1D_MODEL)

    case (REFERENCE_MODEL_IASP91)
      vpc=11.24094d0-4.09689d0*x*x
      vsc=3.56454d0-3.45241d0*x*x
      rhoc=13.0885d0-8.8381d0*x*x
      ! checks with given values
      if (abs(vpc-vp) > TINYVAL .or. abs(vsc-vs) > TINYVAL .or. abs(rhoc-rho) > TINYVAL) then
        stop 'Error isotropic IASP91 values in model_aniso_inner_core() '
      endif

      ! values at center
      vp0=11.24094d0
      vs0=3.56454d0
      rho0=13.0885d0

    case (REFERENCE_MODEL_PREM)
      vpc=11.2622d0-6.3640d0*x*x
      vsc=3.6678d0-4.4475d0*x*x
      rhoc=13.0885d0-8.8381d0*x*x
      ! checks
      if (abs(vpc-vp) > TINYVAL .or. abs(vsc-vs) > TINYVAL .or. abs(rhoc-rho) > TINYVAL) then
        stop 'Error isotropic PREM values in model_aniso_inner_core() '
      endif

      ! values at center
      vp0=11.2622d0
      vs0=3.6678d0
      rho0=13.0885d0

    case (REFERENCE_MODEL_1DREF)
      ! values at center
      vp0 = 11262.20 / 1000.0d0
      vs0 = 3667.800 / 1000.0d0
      rho0 = 13088.480 / 1000.0d0

    case (REFERENCE_MODEL_1066A)
      ! values at center
      vp0 = 11.33830
      vs0 = 3.62980
      rho0 = 13.429030

    case (REFERENCE_MODEL_AK135F_NO_MUD)
      ! values at center
      vp0 = 11.26220
      vs0 = 3.667800
      rho0 = 13.01220

    case (REFERENCE_MODEL_JP1D)
      ! values at center
      vp0 = 11.24094
      vs0 = 3.56454
      rho0 = 13.0885d0

    case (REFERENCE_MODEL_SEA1D)
      ! values at center
      vp0 = 11.240940
      vs0 = 3.564540
      rho0 = 13.012190

    case default
      stop 'unknown 1D reference Earth model in anisotropic inner core'

  end select

! non-dimensionalization of elastic parameters (GPa--[g/cm^3][(km/s)^2])
  scale_fac = RHOAV*R_EARTH*R_EARTH*PI*GRAV*RHOAV
  scale_fac = 1.d9 / scale_fac

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

! Ishii et al. (2002):
!
! alpha = 3.490 % = (C-A)/A0    = (c33-c11)/A0
!  beta = 0.988 % = (L-N)/A0    = (c44-c66)/A0
! gamma = 0.881 % = (A-2N-F)/A0    = (c12-c13)/A0
! where A0 is A at the Earth's center
!
! assume c11 = lambda+2mu
!        c66 = (c11-c12)/2 = mu
!
! then   c33 = c11 + alpha*A0
!        c44 = c66 + beta*A0
!        c13 = c12 - gamma*A0
!        and c12 = c11 - 2*c66
!
! Note: The value vs0, while set above, is unnecessary for the
! calculation below.
!
! Steinle-Neumann (2002):
!
!  r    T    rho    c11   c12  c13  c33  c44 KS   mu
! (km) (K) (Mg/m3) (GPa)
! 0    5735 13.09   1693 1253 1364 1813 154 1457 184
! 200  5729 13.08   1689 1251 1362 1809 154 1455 184
! 400  5711 13.05   1676 1243 1353 1795 151 1444 181
! 600  5682 13.01   1661 1232 1341 1779 150 1432 180
! 800  5642 12.95   1638 1214 1321 1755 148 1411 178
! 1000 5590 12.87   1606 1190 1295 1720 146 1383 175
! 1200 5527 12.77   1559 1155 1257 1670 141 1343 169
!
  c11 = rho*vp*vp*scale_fac
  c66 = rho*vs*vs*scale_fac
  c12 = c11 - 2.0d0*c66

  A0 = rho0*vp0*vp0*scale_fac

  c33 = c11 + 0.0349d0*A0
  c44 = c66 + 0.00988d0*A0
  c13 = c12 - 0.00881d0*A0

  end subroutine model_aniso_inner_core

