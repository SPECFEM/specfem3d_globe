!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine mesh_radial(myrank,rn,rmin,rmax,NER,NER_TOP_CENTRAL_CUBE_ICB, &
       NER_CMB_TOPDDOUBLEPRIME,NER_TOPDDOUBLEPRIME_771, &
       NER_771_670,NER_670_600,NER_600_400,NER_400_220,NER_220_MOHO,NER_CRUST, &
       NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC, &
       CRUSTAL,TOPOGRAPHY,ONE_CRUST)

! create the radial mesh, honoring the major discontinuities in PREM
! we also honor the top of D" and the other second-order discontinuities

  implicit none

  include "constants.h"

! in the case of a very fine mesh, move the bottom of crustal elements
! below the PREM Moho, otherwise the elements become too distorted
  double precision, parameter :: RMOHO_FICTITIOUS_2ELEMS = 6330000.d0
  double precision, parameter :: RMOHO_FICTITIOUS_4ELEMS = 6330000.d0

  integer myrank

  integer NER
  integer NER_TOP_CENTRAL_CUBE_ICB,NER_CMB_TOPDDOUBLEPRIME, &
       NER_TOPDDOUBLEPRIME_771,NER_771_670,NER_670_600,NER_600_400, &
       NER_400_220,NER_220_MOHO,NER_CRUST,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB

  logical CRUSTAL,TOPOGRAPHY,ONE_CRUST

  double precision RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC
  double precision rn(0:2*NER),rmin(0:2*NER),rmax(0:2*NER)

  integer npr,ir
  integer NER_CRUST_LOWER,NER_CRUST_UPPER

  double precision R_0,r_moho,r_220,r_400,r_600,r_670,r_771
  double precision r_cmb,r_bottomdblfluid,r_topdblfluid,r_icb
  double precision r_top_ddoubleprime,r_surface,r_middle_crust

! normalized radii
  R_0 = R_UNIT_SPHERE
  r_middle_crust = RMIDDLE_CRUST / R_EARTH
  r_220 = R220 / R_EARTH
  r_400 = R400 / R_EARTH
  r_600 = R600 / R_EARTH
  r_670 = R670 / R_EARTH
  r_771 = R771 / R_EARTH
  r_top_ddoubleprime = RTOPDDOUBLEPRIME / R_EARTH
  r_cmb = RCMB / R_EARTH
  r_topdblfluid = (RICB + RATIO_TOP_DBL_OC*(RCMB-RICB)) / R_EARTH
  r_bottomdblfluid = (RICB + RATIO_BOTTOM_DBL_OC*(RCMB-RICB)) / R_EARTH
  r_icb = RICB / R_EARTH

! in the case of a very fine mesh, move the bottom of crustal elements
! below the PREM Moho, otherwise the elements become too distorted
!! DK DK this is not very clean, should write something more general one day
!! DK DK do not change anything if model is spherically symmetric (PREM)
  r_moho = RMOHO / R_EARTH
  if(CRUSTAL .or. TOPOGRAPHY) then
    if(NER_CRUST == 1) then
      r_moho = RMOHO / R_EARTH
    else if(NER_CRUST == 2) then
      r_moho = RMOHO_FICTITIOUS_2ELEMS / R_EARTH
    else if(NER_CRUST == 3 .or. NER_CRUST == 4) then
      r_moho = RMOHO_FICTITIOUS_4ELEMS / R_EARTH
    else
      stop 'incorrect definition of r_moho in mesh_radial'
    endif
  endif

! define surface for mapping to [0,1]
  r_surface = R_0

  npr = -1

!
!--- bottom of the mesh (R_CENTRAL_CUBE) to ICB
!
  do ir=0,2*NER_TOP_CENTRAL_CUBE_ICB-1
    npr=npr+1
    rn(npr)=((r_icb-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_TOP_CENTRAL_CUBE_ICB)
    rmin(npr)=R_CENTRAL_CUBE
    rmax(npr)=r_icb
  enddo

!
!--- ICB to bottom doubling fluid
!
  do ir=0,2*NER_ICB_BOTTOMDBL-1
    npr=npr+1
    rn(npr)=(r_icb-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_bottomdblfluid-r_icb)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_ICB_BOTTOMDBL)
    rmin(npr)=r_icb
    rmax(npr)=r_cmb
  enddo

!
!--- bottom doubling fluid to top doubling fluid
!
  do ir=0,2*NER_BOTTOMDBL_TOPDBL-1
    npr=npr+1
    rn(npr)=(r_bottomdblfluid-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_topdblfluid-r_bottomdblfluid)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_BOTTOMDBL_TOPDBL)
    rmin(npr)=r_icb
    rmax(npr)=r_cmb
  enddo

!
!--- top doubling fluid to CMB
!
  do ir=0,2*NER_TOPDBL_CMB-1
    npr=npr+1
    rn(npr)=(r_topdblfluid-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_cmb-r_topdblfluid)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_TOPDBL_CMB)
    rmin(npr)=r_icb
    rmax(npr)=r_cmb
  enddo

!
!--- CMB to top of D"
!
  do ir=0,2*NER_CMB_TOPDDOUBLEPRIME-1
    npr=npr+1
    rn(npr)=(r_cmb-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_top_ddoubleprime-r_cmb)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_CMB_TOPDDOUBLEPRIME)
    rmin(npr)=r_cmb
    rmax(npr)=r_top_ddoubleprime
  enddo
!
!--- top of D" to d771
!
  do ir=0,2*NER_TOPDDOUBLEPRIME_771-1
    npr=npr+1
    rn(npr)=(r_top_ddoubleprime-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_771-r_top_ddoubleprime)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_TOPDDOUBLEPRIME_771)
    rmin(npr)=r_top_ddoubleprime
    rmax(npr)=r_771
  enddo
!
!--- d771 to d670
!
  do ir=0,2*NER_771_670-1
    npr=npr+1
    rn(npr)=(r_771-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_670-r_771)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_771_670)
    rmin(npr)=r_771
    rmax(npr)=r_670
  enddo
!
!--- d670 to d600
!
  do ir=0,2*NER_670_600-1
    npr=npr+1
    rn(npr)=(r_670-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_600-r_670)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_670_600)
    rmin(npr)=r_670
    rmax(npr)=r_600
  enddo
!
!--- d600 to d400
!
  do ir=0,2*NER_600_400-1
    npr=npr+1
    rn(npr)=(r_600-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_400-r_600)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_600_400)
    rmin(npr)=r_600
    rmax(npr)=r_400
  enddo
!
!--- d400 to d220
!
  do ir=0,2*NER_400_220-1
    npr=npr+1
    rn(npr)=(r_400-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_220-r_400)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_400_220)
    rmin(npr)=r_400
    rmax(npr)=r_220
  enddo

!! DK DK avoid problem with bathymetry trenches on the ES
!! DK DK do not honor the fictitious Moho if high-res mesh with topography
!! DK DK also use regular mesh if regional code
  if(REGIONAL_CODE .or. (NER_CRUST > 1 .and. (CRUSTAL .or. TOPOGRAPHY))) then

!! DK DK uniform radial mesh from d220 to surface if high-res 3D model

! also create last point exactly at the surface
! other regions above stop one point below
  do ir=0,2*(NER_220_MOHO+NER_CRUST) - 0
    npr=npr+1
    rn(npr)=(r_220-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
              +((R_0-r_220)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*(NER_220_MOHO+NER_CRUST))
    rmin(npr)=r_220
    rmax(npr)=R_0
  enddo

  else

!
!--- d220 to Moho
!
  do ir=0,2*NER_220_MOHO-1
    npr=npr+1
    rn(npr)=(r_220-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
            +((r_moho-r_220)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_220_MOHO)
    rmin(npr)=r_220
    rmax(npr)=r_moho
  enddo
!
!--- Moho to surface
!
! also create last point exactly at the surface
! other regions above stop one point below

!! DK DK this to compute PREM at high-resolution
!! DK DK this is not very clean, should write something more general one day
  if(.not. ONE_CRUST .and. .not. CRUSTAL .and. .not. TOPOGRAPHY .and. NER_CRUST > 1) then

!! DK DK this is not flexible, should write something more general one day
!! DK DK NER_CRUST = 3 on the ES in Japan, impose 2 in upper and 1 in lower
    NER_CRUST_UPPER = 2
    NER_CRUST_LOWER = NER_CRUST - NER_CRUST_UPPER

! lower crust
    do ir=0,2*NER_CRUST_LOWER - 1
      npr=npr+1
      rn(npr)=(r_moho-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
              +((r_middle_crust-r_moho)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_CRUST_LOWER)
      rmin(npr)=r_moho
      rmax(npr)=R_0
    enddo

! upper crust
    do ir=0,2*NER_CRUST_UPPER - 0
      npr=npr+1
      rn(npr)=(r_middle_crust-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
              +((R_0-r_middle_crust)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_CRUST_UPPER)
      rmin(npr)=r_moho
      rmax(npr)=R_0
    enddo

  else
! regular mesh here
    do ir=0,2*NER_CRUST - 0
      npr=npr+1
      rn(npr)=(r_moho-R_CENTRAL_CUBE)/(r_surface-R_CENTRAL_CUBE) &
              +((R_0-r_moho)/(r_surface-R_CENTRAL_CUBE))*dble(ir)/dble(2*NER_CRUST)
      rmin(npr)=r_moho
      rmax(npr)=R_0
    enddo
  endif

  endif

! check that the mesh that has been generated is correct
! radial distribution is mapped to [0,1]
  if(npr /= 2*NER .or. rn(0) < ZERO .or. rn(2*NER) > ONE) &
      call exit_MPI(myrank,'incorrect intervals for global Earth')

! check that radial spacing makes sense
  do ir=0,2*NER-1
    if(rn(ir+1) < rn(ir)) call exit_MPI(myrank,'incorrect radial spacing for global Earth')
  enddo

  end subroutine mesh_radial

