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

  subroutine make_gravity(nspl,rspl,gspl,gspl2,ONE_CRUST)

! creates a spline for the gravity profile in PREM
! radius and density are non-dimensional

  use constants, only: NR
  use shared_parameters, only: PLANET_TYPE,IPLANET_MARS

  implicit none

  integer,intent(out) :: nspl

  double precision,intent(inout) :: rspl(NR),gspl(NR),gspl2(NR)

  logical,intent(in) :: ONE_CRUST

  ! local parameters
  integer :: i
  double precision :: ROCEAN,RMIDDLE_CRUST,RMOHO,R80,R220,R400,R600,R670, &
                      R771,RTOPDDOUBLEPRIME,RCMB,RICB,RSURFACE
  double precision :: r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision :: r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision :: r(NR),rho(NR),g(NR),i_rho
  double precision :: s1(NR),s2(NR),s3(NR)
  double precision :: yp1,ypn

  ! Earth
  ! radius of the Earth for gravity calculation
  double precision, parameter :: R_EARTH_GRAVITY = 6371000.d0
  ! radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_GRAVITY = 6368000.d0
  ! Mars
  ! radius of Mars for gravity calculation
  double precision, parameter :: R_MARS_GRAVITY = 3390000.d0
  ! radius of the ocean floor for gravity calculation
  double precision, parameter :: ROCEAN_MARS_GRAVITY = 3390000.d0

  ! sets radii
  select case (PLANET_TYPE)
  case (IPLANET_MARS)
    ! Mars
    ! Sohl & Spohn (Mars Model A)
    RSURFACE = R_MARS_GRAVITY
    ROCEAN = ROCEAN_MARS_GRAVITY
    RMIDDLE_CRUST = 3340000.d0
    RMOHO = 3280000.d0 ! moho depth at 110 km
    R80  = 3055000.d0
    R220 = 2908000.d0
    R400 = 2655000.d0
    R600 = 2455000.d0
    R670 = 2360000.d0
    R771 = 2033000.d0
    RTOPDDOUBLEPRIME = 1503000.d0
    RCMB = 1468000.d0
    RICB = 515000.d0

  case default
    ! Earth default
    ! PREM
    RSURFACE = R_EARTH_GRAVITY
    ROCEAN = ROCEAN_GRAVITY ! PREM defines this as 6368000.d0
    RMIDDLE_CRUST = 6356000.d0
    RMOHO = 6346600.d0 ! PREM moho depth at 24.4 km
    R80  = 6291000.d0
    R220 = 6151000.d0
    R400 = 5971000.d0
    R600 = 5771000.d0
    R670 = 5701000.d0
    R771 = 5600000.d0
    RTOPDDOUBLEPRIME = 3630000.d0
    RCMB = 3480000.d0
    RICB = 1221000.d0
  end select

  ! non-dimensionalize
  r_icb = RICB/RSURFACE
  r_cmb = RCMB/RSURFACE
  r_topddoubleprime = RTOPDDOUBLEPRIME/RSURFACE
  r_771 = R771/RSURFACE
  r_670 = R670/RSURFACE
  r_600 = R600/RSURFACE
  r_400 = R400/RSURFACE
  r_220 = R220/RSURFACE
  r_80 = R80/RSURFACE
  r_moho = RMOHO/RSURFACE
  r_middle_crust = RMIDDLE_CRUST/RSURFACE
  r_ocean = ROCEAN/RSURFACE
  r_0 = 1.d0

! note: for Mars
!       discretize the radius point to fit the spline, inherited from PREM
! To do: (may not be necessary)
!        Find number of points at each layer to have better integrations

  do i=1,163
    r(i) = r_icb*dble(i-1)/dble(162)
  enddo
  do i=164,323
    r(i) = r_icb+(r_cmb-r_icb)*dble(i-164)/dble(159)
  enddo
  do i=324,336
    r(i) = r_cmb+(r_topddoubleprime-r_cmb)*dble(i-324)/dble(12)
  enddo
  do i=337,517
    r(i) = r_topddoubleprime+(r_771-r_topddoubleprime)*dble(i-337)/dble(180)
  enddo
  do i=518,530
    r(i) = r_771+(r_670-r_771)*dble(i-518)/dble(12)
  enddo
  do i=531,540
    r(i) = r_670+(r_600-r_670)*dble(i-531)/dble(9)
  enddo
  do i=541,565
    r(i) = r_600+(r_400-r_600)*dble(i-541)/dble(24)
  enddo
  do i=566,590
    r(i) = r_400+(r_220-r_400)*dble(i-566)/dble(24)
  enddo
  do i=591,609
    r(i) = r_220+(r_80-r_220)*dble(i-591)/dble(18)
  enddo
  do i=610,619
    r(i) = r_80+(r_moho-r_80)*dble(i-610)/dble(9)
  enddo
  do i=620,626
    r(i) = r_moho+(r_middle_crust-r_moho)*dble(i-620)/dble(6)
  enddo
  do i=627,633
    r(i) = r_middle_crust+(r_ocean-r_middle_crust)*dble(i-627)/dble(6)
  enddo
  do i=634,NR
    r(i) = r_ocean+(r_0-r_ocean)*dble(i-634)/dble(6)
  enddo

  ! density profile
  if (PLANET_TYPE == IPLANET_MARS) then
    ! Mars
    ! Sohn & Spohn Model A
    ! No Ocean
    do i = 627,NR
      r(i) = r_middle_crust+(r_0-r_middle_crust)*dble(i-627)/dble(12)
    enddo
    ! use Sohl & Spohn model A to get the density profile for gravity
    do i = 1,NR
      call sohl_density(r(i),rho(i),ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
                        R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)
    enddo

  else
    ! Earth
    ! use PREM to get the density profile for ellipticity (fine for other 1D reference models)
    do i = 1,NR
      call prem_density(r(i),rho(i),ONE_CRUST)
    enddo
  endif

  g(1) = 0.0d0
  do i = 2,NR
    call intgrl(i_rho,r,1,i,rho,s1,s2,s3)
    g(i) = 4.0d0*i_rho/(r(i)*r(i))
  enddo

!
! get ready to spline g
!
  nspl = 1
  rspl(1) = r(1)
  gspl(1) = g(1)
  do i=2,NR
    if (r(i) /= r(i-1)) then
      nspl = nspl+1
      rspl(nspl) = r(i)
      gspl(nspl) = g(i)
    endif
  enddo
  yp1 = (4.0d0/3.0d0)*rho(1)
  ypn = 4.0d0*rho(NR)-2.0d0*g(NR)/r(NR)

  call spline_construction(rspl,gspl,nspl,yp1,ypn,gspl2)

  end subroutine make_gravity

