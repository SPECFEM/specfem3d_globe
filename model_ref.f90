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

! 1D REF model of Kustowski et al. (2007)
  open(unit=10,file='DATA/s362ani/REF',status='old',action='read')
! skip the first three lines
  read(10,*)
  read(10,*)
  read(10,*)
! read the card deck model
  do i=1,NR_REF
    read(10,*) Mref_V%radius_ref(i),Mref_V%density_ref(i),Mref_V%vpv_ref(i), &
               Mref_V%vsv_ref(i),Mref_V%Qkappa_ref(i),Mref_V%Qmu_ref(i), &
               Mref_V%vph_ref(i),Mref_V%vsh_ref(i),Mref_V%eta_ref(i)
  enddo
  close(10)

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

