!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine make_gravity(nspl,rspl,gspl,gspl2,ONE_CRUST)

! creates a spline for the gravity profile in PREM
! radius and density are non-dimensional

  implicit none

  include "constants.h"

  integer nspl
  logical ONE_CRUST
  double precision rspl(NR),gspl(NR),gspl2(NR)

  integer i
  double precision r_icb,r_cmb,r_topddoubleprime,r_771,r_670,r_600
  double precision r_400,r_220,r_80,r_moho,r_middle_crust,r_ocean,r_0
  double precision r(NR),rho(NR),g(NR),i_rho
  double precision s1(NR),s2(NR),s3(NR)
  double precision yp1,ypn

  r_icb = RICB/R_EARTH
  r_cmb = RCMB/R_EARTH
  r_topddoubleprime = RTOPDDOUBLEPRIME/R_EARTH
  r_771 = R771/R_EARTH
  r_670 = R670/R_EARTH
  r_600 = R600/R_EARTH
  r_400 = R400/R_EARTH
  r_220 = R220/R_EARTH
  r_80 = R80/R_EARTH
  r_moho = RMOHO/R_EARTH
  r_middle_crust = RMIDDLE_CRUST/R_EARTH
  r_ocean = ROCEAN/R_EARTH
  r_0 = 1.d0

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

  do i=1,NR
    call prem_density(r(i),rho(i),ONE_CRUST)
  enddo

  g(1)=0.0d0
  do i=2,NR
    call intgrl(i_rho,r,1,i,rho,s1,s2,s3)
    g(i)=4.0d0*i_rho/(r(i)*r(i))
  enddo

!
! get ready to spline g
!
  nspl=1
  rspl(1)=r(1)
  gspl(1)=g(1)
  do i=2,NR
    if(r(i)/=r(i-1)) then
      nspl=nspl+1
      rspl(nspl)=r(i)
      gspl(nspl)=g(i)
    endif
  enddo
  yp1=(4.0d0/3.0d0)*rho(1)
  ypn=4.0d0*rho(NR)-2.0d0*g(NR)/r(NR)
  call spline(rspl,gspl,nspl,yp1,ypn,gspl2)

  end subroutine make_gravity

