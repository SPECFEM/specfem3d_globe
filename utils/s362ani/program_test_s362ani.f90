!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2013
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


  program xtest_s362ani

  implicit none

  include "constants.h"

! ocean-continent function maximum spherical harmonic degree
  integer,parameter :: NL_OCEAN_CONTINENT = 12

! spherical harmonic coefficients of the ocean-continent function (km)
  double precision  A_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT),B_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT)

  common /ocf/ A_lm,B_lm

! minimum and maximum allowable topography of the Moho (km) for stretching purposes
  double precision,parameter :: EMAX = 1.0d0/R_EARTH_KM
  double precision,parameter :: RMOHO = 6346600.d0

! model_ref_variables
  type model_ref_variables
    sequence
      double precision,dimension(NR_REF) :: radius_ref
      double precision,dimension(NR_REF) :: density_ref
      double precision,dimension(NR_REF) :: vpv_ref
      double precision,dimension(NR_REF) :: vph_ref
      double precision,dimension(NR_REF) :: vsv_ref
      double precision,dimension(NR_REF) :: vsh_ref
      double precision,dimension(NR_REF) :: eta_ref
      double precision,dimension(NR_REF) :: Qkappa_ref
      double precision,dimension(NR_REF) :: Qmu_ref
  end type model_ref_variables

! type (model_ref_variables) Mref_V
! model_ref_variables

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision,dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision,dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision,dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision,dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
    character(len=2) abbreviation(NCAP_CRUST/2,NCAP_CRUST)
    character(len=2) code(NKEYS_CRUST)
  end type crustal_model_variables

  type (crustal_model_variables) CM_V
! crustal_model_variables


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

  integer ilat,ilon
! integer i,iregion_code
! integer THREE_D_MODEL

! real(kind=4) xlat,xcolat,xlon,xdep,xrad
! real(kind=4) dvsh,dvsv,dvph,dvpv
! real(kind=4) topo410,topo650
! double precision scaleval
! double precision xr,vpv,vph,vsv,vsh,eta,Qmu,Qkappa
  double precision rho
  double precision lat,lon,r,vp,vs,moho,r_moho,topomoho,elevation
  logical found_crust

  integer l,m
  double precision theta,phi,sint,cost,x(2*NL_OCEAN_CONTINENT+1),dx(2*NL_OCEAN_CONTINENT+1)
  double precision ocf

!  call define_model_ref(.false.,Mref_V)
!
!  scaleval=dsqrt(PI*GRAV*RHOAV)
!  do i=1,NR_REF
!    x = dble(i-1)/dble(NR_REF-1)
!    call model_ref(xr,rho,vpv,vph,vsv,vsh,eta,Qkappa,Qmu,iregion_code,Mref_V)
!    print *,i,sngl(x*R_EARTH),sngl(rho*RHOAV),sngl(vpv*(R_EARTH*scaleval)),sngl(vph*(R_EARTH*scaleval)), &
!            sngl(vsv*(R_EARTH*scaleval)),sngl(vsh*(R_EARTH*scaleval)),sngl(eta),sngl(Qmu)
!
! find and check discontinuities
!    if (Mref_V%radius_ref(i) == Mref_V%radius_ref(i+1)) then
!      print *,i,Mref_V%radius_ref(i)
!    endif
!
!  enddo
!
!  print *,"xlat: "
!  read(5,*) xlat
!  xcolat=90.0-xlat
!  print *,"xlon: "
!  read(5,*) xlon
!  print *,"xdep: "
!  read(5,*) xdep
!  xrad=6371.0-xdep
!  print *,"THREE_D_MODEL: "
!  read(5,*) THREE_D_MODEL
!!
!  call read_model_s362ani(THREE_D_MODEL,THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
!               THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA, &
!               numker,numhpa,ihpa,lmxhpa,itypehpa,ihpakern,numcoe,ivarkern,itpspl, &
!               xlaspl,xlospl,radspl,coe,hsplfl,dskker,kerstr,varstr,refmdl)
!
!  call subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv, &
!               numker,numhpa,numcof,ihpa,lmax,nylm, &
!               lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
!               nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
!               coe,vercof,vercofd,ylmcof,wk1,wk2,wk3,kerstr,varstr)
!  write(*,"('    dvsh      dvsv      dvph      dvpv    ')")
!  write(*,"(6f10.5)") 100.0*dvsh,100.0*dvsv,100.0*dvph,100.0*dvpv
!
!  call subtopo(xcolat,xlon,topo410,topo650, &
!               numker,numhpa,numcof,ihpa,lmax,nylm, &
!               lmxhpa,itypehpa,ihpakern,numcoe,ivarkern, &
!               nconpt,iver,iconpt,conpt,xlaspl,xlospl,radspl, &
!               coe,ylmcof,wk1,wk2,wk3,varstr)
!  write(*,"('   topo410    topo650 ')")
!  write(*,"(2f11.5)") topo410,topo650
!

! ocean-continent function (km)
  call read_smooth_moho

! test moho stretching
  call read_crustal_model(CM_V)

  r = 0.9999d0
  do ilat = -89,89,5
    do ilon = -179,179,5

      lat = dble(ilat)
      lon = dble(ilon)

      theta = (90.0d0-lat)*PI/180.0d0
      phi = lon*PI/180.0d0

      ocf=0.0d0
      do l=0,NL_OCEAN_CONTINENT
        sint = dsin(theta)
        cost = dcos(theta)
        call lgndr(l,cost,sint,x,dx)
        m = 0
        ocf = ocf +  A_lm(l,m)*x(m+1)
        do m=1,l
          ocf = ocf + (A_lm(l,m)*dcos(dble(m)*phi)+B_lm(l,m)*dsin(dble(m)*phi))*x(m+1)
        enddo
      enddo
      ocf = -ocf

! get Moho depth
      call crustal_model(lat,lon,r,vp,vs,rho,moho,found_crust,CM_V)
      if (.not. found_crust) stop 'failed to find Moho in add_topography_moho'
! Moho radius
      r_moho = 1.0d0 - moho

! Moho topography: negative for a depression,positive for an elevation
      topomoho = r_moho - RMOHO/R_EARTH
      print *,sngl(lat),sngl(lon),sngl(moho*R_EARTH_KM),sngl(topomoho*R_EARTH_KM),sngl(ocf)

! smoothed version of topographic variations
     elevation = 0.5d0*EMAX
     if (topomoho < EMAX) elevation = 0.5d0 * topomoho

!      print *,sngl(lat),sngl(lon),sngl(moho*R_EARTH_KM),sngl(r_moho*R_EARTH_KM-RMOHO/1000.0d0),sngl(elevation*R_EARTH_KM)
    enddo
  enddo

  end program xtest_s362ani

  subroutine read_smooth_moho

  implicit none

! ocean-continent function maximum spherical harmonic degree
  integer,parameter :: NL_OCEAN_CONTINENT = 12

! spherical harmonic coefficients of the ocean-continent function (km)
  double precision  A_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT),B_lm(0:NL_OCEAN_CONTINENT,0:NL_OCEAN_CONTINENT)

  common /ocf/ A_lm,B_lm

  A_lm(0,0) = -3.8201999E-04
  B_lm(0,0) = 0.
  A_lm(1,0) = 13.88800
  B_lm(1,0) = 0.
  A_lm(1,1) = -15.24000
  B_lm(1,1) = -9.187200
  A_lm(2,0) = 11.21500
  B_lm(2,0) = 0.
  A_lm(2,1) = -6.754500
  B_lm(2,1) = -8.516700
  A_lm(2,2) = -8.327800
  B_lm(2,2) = -5.029200
  A_lm(3,0) = -3.614500
  B_lm(3,0) = 0.
  A_lm(3,1) = 5.394800
  B_lm(3,1) = -0.9220800
  A_lm(3,2) = -10.05100
  B_lm(3,2) = 13.98100
  A_lm(3,3) = -2.711200
  B_lm(3,3) = -13.57100
  A_lm(4,0) = 7.523300
  B_lm(4,0) = 0.
  A_lm(4,1) = 5.156100
  B_lm(4,1) = 2.184400
  A_lm(4,2) = -10.67300
  B_lm(4,2) = 2.640600
  A_lm(4,3) = -7.786300
  B_lm(4,3) = 0.3674500
  A_lm(4,4) = -3.076400
  B_lm(4,4) = 16.83000
  A_lm(5,0) = -9.681000
  B_lm(5,0) = 0.
  A_lm(5,1) = 0.5026800
  B_lm(5,1) = 2.111300
  A_lm(5,2) = -2.931000
  B_lm(5,2) = -4.329000
  A_lm(5,3) = -1.766800
  B_lm(5,3) = -3.621200
  A_lm(5,4) = 16.08200
  B_lm(5,4) = -4.493900
  A_lm(5,5) = -0.3705800
  B_lm(5,5) = -5.574500
  A_lm(6,0) = 4.407900
  B_lm(6,0) = 0.
  A_lm(6,1) = 0.3799000
  B_lm(6,1) = 1.589400
  A_lm(6,2) = -1.886400
  B_lm(6,2) = -0.5686300
  A_lm(6,3) = -0.9816800
  B_lm(6,3) = -5.827800
  A_lm(6,4) = 3.620600
  B_lm(6,4) = -2.713100
  A_lm(6,5) = 1.445600
  B_lm(6,5) = 3.964100
  A_lm(6,6) = 1.167400
  B_lm(6,6) = 2.134100
  A_lm(7,0) = -4.086100
  B_lm(7,0) = 0.
  A_lm(7,1) = 0.5462000
  B_lm(7,1) = -4.488100
  A_lm(7,2) = 3.116400
  B_lm(7,2) = 1.793600
  A_lm(7,3) = 2.594600
  B_lm(7,3) = -2.129100
  A_lm(7,4) = -5.445000
  B_lm(7,4) = 0.5381500
  A_lm(7,5) = -2.178100
  B_lm(7,5) = 1.766700
  A_lm(7,6) = -1.040000
  B_lm(7,6) = -5.541000
  A_lm(7,7) = 1.536500
  B_lm(7,7) = 3.700600
  A_lm(8,0) = -2.562200
  B_lm(8,0) = 0.
  A_lm(8,1) = 0.3736200
  B_lm(8,1) = 1.488000
  A_lm(8,2) = 1.347500
  B_lm(8,2) = 0.5288200
  A_lm(8,3) = -0.8493700
  B_lm(8,3) = -1.626500
  A_lm(8,4) = 0.2423400
  B_lm(8,4) = 4.202800
  A_lm(8,5) = 2.052200
  B_lm(8,5) = 0.6880400
  A_lm(8,6) = 2.838500
  B_lm(8,6) = 2.835700
  A_lm(8,7) = -4.981400
  B_lm(8,7) = -1.883100
  A_lm(8,8) = -1.102800
  B_lm(8,8) = -1.951700
  A_lm(9,0) = -1.202100
  B_lm(9,0) = 0.
  A_lm(9,1) = 1.020300
  B_lm(9,1) = 1.371000
  A_lm(9,2) = -0.3430100
  B_lm(9,2) = 0.8782800
  A_lm(9,3) = -0.4462500
  B_lm(9,3) = -0.3046100
  A_lm(9,4) = 0.7750700
  B_lm(9,4) = 2.351600
  A_lm(9,5) = -2.092600
  B_lm(9,5) = -2.377100
  A_lm(9,6) = 0.3126900
  B_lm(9,6) = 4.996000
  A_lm(9,7) = -2.284000
  B_lm(9,7) = 1.183700
  A_lm(9,8) = 1.445900
  B_lm(9,8) = 1.080000
  A_lm(9,9) = 1.146700
  B_lm(9,9) = 1.457800
  A_lm(10,0) = -2.516900
  B_lm(10,0) = 0.
  A_lm(10,1) = -0.9739500
  B_lm(10,1) = -0.7195500
  A_lm(10,2) = -2.846000
  B_lm(10,2) = -1.464700
  A_lm(10,3) = 2.720100
  B_lm(10,3) = 0.8241400
  A_lm(10,4) = -1.247800
  B_lm(10,4) = 1.220300
  A_lm(10,5) = -1.638500
  B_lm(10,5) = -1.099500
  A_lm(10,6) = 3.043000
  B_lm(10,6) = -1.976400
  A_lm(10,7) = -1.007300
  B_lm(10,7) = -1.604900
  A_lm(10,8) = 0.6620500
  B_lm(10,8) = -1.135000
  A_lm(10,9) = -3.576800
  B_lm(10,9) = 0.5554900
  A_lm(10,10) = 2.418700
  B_lm(10,10) = -1.482200
  A_lm(11,0) = 0.7158800
  B_lm(11,0) = 0.
  A_lm(11,1) = -3.694800
  B_lm(11,1) = 0.8491400
  A_lm(11,2) = 9.3208998E-02
  B_lm(11,2) = -1.276000
  A_lm(11,3) = 1.575600
  B_lm(11,3) = 0.1972100
  A_lm(11,4) = 0.8989600
  B_lm(11,4) = -1.063000
  A_lm(11,5) = -0.6301000
  B_lm(11,5) = -1.329400
  A_lm(11,6) = 1.389000
  B_lm(11,6) = 1.184100
  A_lm(11,7) = 0.5640700
  B_lm(11,7) = 2.286200
  A_lm(11,8) = 1.530300
  B_lm(11,8) = 0.7677500
  A_lm(11,9) = 0.8495500
  B_lm(11,9) = 0.7247500
  A_lm(11,10) = 2.106800
  B_lm(11,10) = 0.6588000
  A_lm(11,11) = 0.6067800
  B_lm(11,11) = 0.1366800
  A_lm(12,0) = -2.598700
  B_lm(12,0) = 0.
  A_lm(12,1) = -1.150500
  B_lm(12,1) = -0.8425700
  A_lm(12,2) = -0.1593300
  B_lm(12,2) = -1.241400
  A_lm(12,3) = 1.508600
  B_lm(12,3) = 0.3385500
  A_lm(12,4) = -1.941200
  B_lm(12,4) = 1.120000
  A_lm(12,5) = -0.4630500
  B_lm(12,5) = -6.4753003E-02
  A_lm(12,6) = 0.8967000
  B_lm(12,6) = 4.7417998E-02
  A_lm(12,7) = 4.5407999E-02
  B_lm(12,7) = 0.8876400
  A_lm(12,8) = -2.444400
  B_lm(12,8) = 1.172500
  A_lm(12,9) = -2.593400
  B_lm(12,9) = 0.1703700
  A_lm(12,10) = 0.5662700
  B_lm(12,10) = 0.7050800
  A_lm(12,11) = -0.1930000
  B_lm(12,11) = -2.008100
  A_lm(12,12) = -3.187900
  B_lm(12,12) = -1.672000

  end subroutine read_smooth_moho

