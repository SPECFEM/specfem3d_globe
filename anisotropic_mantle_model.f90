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

!=====================================================================
!
!       Jean-Paul Montagner, January 2002
!       modified by Min Chen, Caltech, February 2002
!
!  input is (r, theta, phi), output is the matrix cij(6x6)
!  0 <= r <= 1, 0 <= theta <= pi, 0 <= phi <= 2 pi
!
!  returns non-dimensionalized cij
!
!  creates parameters p(i=1,14,r,theta,phi)
!  from model glob-prem3sm01 globpreman3sm01 (Montagner, 2002)
!
!======================================================================

module aniso_mantle_model_constants

  implicit none

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: DEGREES_TO_RADIANS = PI / 180.d0
  double precision, parameter :: R_EARTH = 6371000.d0
  double precision, parameter :: R_EARTH_KM = R_EARTH / 1000.d0
  double precision, parameter :: ZERO = 0.d0,ONE = 1.d0
  double precision, parameter :: R_UNIT_SPHERE = ONE
  double precision, parameter :: GRAV = 6.6723d-11

! average density in the full Earth to normalize equation
  double precision, parameter :: RHOAV = 5514.3d0

end module aniso_mantle_model_constants

!------------

module aniso_mantle_model_variables

  use aniso_mantle_model_constants

  implicit none

  double precision beta(14,34,37,73),pro(47)
  integer npar1

end module aniso_mantle_model_variables

!-------------------

  subroutine aniso_mantle_model(myrank,r,theta,phi,rho, &
    c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  use aniso_mantle_model_variables

  implicit none

  integer myrank

  double precision r,theta,phi
  double precision rho
  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                   c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision d11,d12,d13,d16,d22,d23,d26,d33,d36,d44,d45,d55,d66
  double precision colat,lon

! uncomment this line to suppress anisotropic mantle models
!  call exit_MPI(myrank,'anisotropic mantle models not implemented yet')

  lon = phi / DEGREES_TO_RADIANS
  colat = theta / DEGREES_TO_RADIANS

  call build_cij(pro,npar1,rho,beta,r,colat,lon,d11,d12,d13,d16,&
        d22,d23,d26,d33,d36,d44,d45,d55,d66,myrank)

  call rotate_aniso_tensor(theta,phi,d11,d12,d13,d16,d22,d23,d26,&
                           d33,d36,d44,d45,d55,d66,&
                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  end subroutine aniso_mantle_model

!--------------------------------------------------------------------

  subroutine build_cij(pro,npar1,rho,beta,r,theta,phi,d11,d12,d13,d16,&
        d22,d23,d26,d33,d36,d44,d45,d55,d66,myrank)

  use aniso_mantle_model_constants

  implicit none

  integer npar1,ndepth,idep,ipar,itheta,ilon,icz0,nx0,ny0,nz0,&
          ict0,ict1,icp0,icp1,icz1,myrank

  double precision d11,d12,d13,d16,d22,d23,d26,d33,d36,d44,d45,d55,d66
  double precision r,theta,phi,rho,depth,tei,tet,ph,fi,x0,y0,pxy0
  double precision d1,d2,d3,d4,sd,thickness,dprof1,dprof2,eps,pc1,pc2,pc3,pc4,&
                   dpr1,dpr2,param,scale_GPa,scaleval
  double precision A,C,F,AL,AN,BC,BS,GC,GS,HC,HS,EC,ES
  double precision beta(14,34,37,73),pro(47)
  double precision anispara(14,2,4),elpar(14)

  ndepth = npar1
  pxy0 = 5.
  x0 = 0.
  y0 = 0.
  nx0 = 37
  ny0 = 73
  nz0 = 34

! avoid edge effects
  if(theta==0.0d0) theta=0.000001d0
  if(theta==180.d0) theta=0.999999d0*theta
  if(phi==0.0d0) phi=0.000001d0
  if(phi==360.d0) phi=0.999999d0*phi

! dimensionalize
  depth = R_EARTH_KM*(R_UNIT_SPHERE - r)
  if(depth <= pro(nz0) .or. depth >= pro(1)) call exit_MPI(myrank,'r out of range in build_cij')
  itheta = int(theta + pxy0)/pxy0
  ilon = int(phi + pxy0)/pxy0
  tet = theta
  ph = phi

  icz0 = 0
  do idep = 1,ndepth
    if(pro(idep) > depth) icz0 = icz0 + 1
  enddo

!
! Interpolation for depth between dep1(iz0) and dep2(iz1)
!
!    1 (ict0,icp0)      2 (ict0,icp1)
!    3 (ict1,icp0)      4 (ict1,icp1)
!

  ict0 = itheta
  ict1 = ict0 + 1
  icp0 = ilon
  icp1 = icp0 + 1
  icz1 = icz0 + 1

! check that parameters make sense
  if(ict0 < 1 .or. ict0 > nx0) call exit_MPI(myrank,'ict0 out of range')
  if(ict1 < 1 .or. ict1 > nx0) call exit_MPI(myrank,'ict1 out of range')
  if(icp0 < 1 .or. icp0 > ny0) call exit_MPI(myrank,'icp0 out of range')
  if(icp1 < 1 .or. icp1 > ny0) call exit_MPI(myrank,'icp1 out of range')
  if(icz0 < 1 .or. icz0 > nz0) call exit_MPI(myrank,'icz0 out of range')
  if(icz1 < 1 .or. icz1 > nz0) call exit_MPI(myrank,'icz1 out of range')

  do ipar = 1,14
    anispara(ipar,1,1) = beta(ipar,icz0,ict0,icp0)
    anispara(ipar,2,1) = beta(ipar,icz1,ict0,icp0)
    anispara(ipar,1,2) = beta(ipar,icz0,ict0,icp1)
    anispara(ipar,2,2) = beta(ipar,icz1,ict0,icp1)
    anispara(ipar,1,3) = beta(ipar,icz0,ict1,icp0)
    anispara(ipar,2,3) = beta(ipar,icz1,ict1,icp0)
    anispara(ipar,1,4) = beta(ipar,icz0,ict1,icp1)
    anispara(ipar,2,4) = beta(ipar,icz1,ict1,icp1)
  enddo

!
! calculation of distances between the selected point and grid points
!
  tei = pxy0*ict0 + x0 - pxy0
  fi = pxy0*icp0 + y0 - pxy0

!***  d1=de(tet,ph,tei,fi)

  d1 = dsqrt(((tei - tet)**2) + ((fi - ph)**2)*(dsin((tet + tei)*DEGREES_TO_RADIANS/2.)**2))

!***  d2=de(tet,ph,tei+pxy0,fi)

  d2 = dsqrt(((tei - tet + pxy0)**2) + ((fi - ph)**2)*(dsin((tet + tei + pxy0)*DEGREES_TO_RADIANS/2.)**2))

!***  d3=de(tet,ph,tei,fi+pxy0)

  d3 = dsqrt(((tei - tet)**2) + ((fi - ph + pxy0)**2)*(dsin((tet + tei)*DEGREES_TO_RADIANS/2.)**2))

!***  d4=de(tet,ph,tei+pxy0,fi+pxy0)

  d4 = dsqrt(((tei - tet + pxy0)**2) + ((fi - ph + pxy0)**2)*(dsin((tet + tei + pxy0)*DEGREES_TO_RADIANS/2.)**2))

  sd = d2*d3*d4 + d1*d2*d4 + d1*d3*d4 + d1*d2*d3
  thickness = pro(icz0) - pro(icz1)
  dprof1 = pro(icz0) - depth
  dprof2 = depth - pro(icz1)
  eps = 0.01

  do ipar = 1,14
     if(thickness < eps)then
      pc1 = anispara(ipar,1,1)
      pc2 = anispara(ipar,1,2)
      pc3 = anispara(ipar,1,3)
      pc4 = anispara(ipar,1,4)
     else
      dpr1 = dprof1/thickness
      dpr2 = dprof2/thickness
      pc1 = anispara(ipar,1,1)*dpr2+anispara(ipar,2,1)*dpr1
      pc2 = anispara(ipar,1,2)*dpr2+anispara(ipar,2,2)*dpr1
      pc3 = anispara(ipar,1,3)*dpr2+anispara(ipar,2,3)*dpr1
      pc4 = anispara(ipar,1,4)*dpr2+anispara(ipar,2,4)*dpr1
     endif
     param = pc1*d2*d3*d4 + pc2*d1*d3*d4 + pc3*d1*d2*d4 + pc4*d1*d2*d3
     param = param/sd
     elpar(ipar) = param
  enddo

  d11 = ZERO
  d12 = ZERO
  d13 = ZERO
  d22 = ZERO
  d33 = ZERO
  d23 = ZERO
  d16 = ZERO
  d26 = ZERO
  d36 = ZERO
  d44 = ZERO
  d45 = ZERO
  d55 = ZERO
  d66 = ZERO

!
!   create dij
!
  rho = elpar(1)
  A = elpar(2)
  C = elpar(3)
  F = elpar(4)
  AL = elpar(5)
  AN = elpar(6)
  BC = elpar(7)
  BS = elpar(8)
  GC = elpar(9)
  GS = elpar(10)
  HC = elpar(11)
  HS = elpar(12)
  EC = elpar(13)
  ES = elpar(14)
  d11 = A + BC + EC
  d12 = A - 2.*AN - EC
  d13 = F + HC
  d22 = A - BC + EC
  d33 = C
  d23 = F - HC
  d16 = BS/2. + ES
  d26 = BS/2. - ES
  d36 = HS
  d44 = AL - GC
  d45 = GS
  d55 = AL + GC
  d66 = AN - EC

! convert these elastic coefficients
! into conventional local Cartesian coordinate system
  d45 = - d45
  d16 = - d16
  d26 = - d26
  d36 = - d36

! non-dimensionalize the elastic coefficients using
! the scale of GPa--[g/cm^3][(km/s)^2]
  scaleval = dsqrt(PI*GRAV*RHOAV)
  scale_GPa =(RHOAV/1000.d0)*((R_EARTH*scaleval/1000.d0)**2)
  d11 = d11/scale_GPa
  d12 = d12/scale_GPa
  d13 = d13/scale_GPa
  d22 = d22/scale_GPa
  d33 = d33/scale_GPa
  d23 = d23/scale_GPa
  d16 = d16/scale_GPa
  d26 = d26/scale_GPa
  d36 = d36/scale_GPa
  d44 = d44/scale_GPa
  d45 = d45/scale_GPa
  d55 = d55/scale_GPa
  d66 = d66/scale_GPa

! non-dimensionalize
  rho = rho*1000.d0/RHOAV

  end subroutine build_cij

!--------------------------------------------------------------

  subroutine read_aniso_mantle_model

  use aniso_mantle_model_variables
  implicit none

  integer nx,ny,np1,np2,ipar,ipa1,ipa,ilat,ilon,il,idep,nfin,nfi0,nf,nri
  double precision xinf,yinf,pxy,ppp,angle,A,A2L,AL,af
  double precision ra(47),pari(14,37)
  double precision bet2(14,34,37,73)
  double precision alph(73,37),ph(73,37)

  np1 = 1
  np2 = 34
  npar1 = (np2 - np1 + 1)

!
! glob-prem3sm01: model with rho,A,L,xi-1,1-phi,eta
!
  open(19,file='DATA/Montagner_model/glob-prem3sm01',status='old')

!
! read the models
!
! reference model: PREM or ACY400
!
  call lecmod(nri,pari,ra)
!
! read tomographic model (equivalent T.I. model)
!
  ipa = 0
  nfi0 = 6
  nfin = 14
  do nf = 1,nfi0
    ipa = ipa + 1
    do idep = 1,npar1
      il = idep + np1 - 1
      read(19,2011,end = 88) xinf,yinf,nx,ny,pxy

      ppp = 1.
      read(19,2013,end = 88) pro(idep),ppp

      if(nf == 1) pari(nf,il) = ppp
      if(nf == 2) pari(nf,il) = ppp
      if(nf == 3) pari(nf,il) = ppp
      if(nf == 4) ppp = pari(nf,il)
      if(nf == 5) ppp = pari(nf,il)
      do ilat = 1,nx
        read(19,2015,end = 88) (beta(ipa,idep,ilat,ilon),ilon = 1,ny)
!
! calculation of A,C,F,L,N
!
! bet2(1,...)=rho, bet2(2,...)=A,bet2(3,...)=L,bet2(4,...)=xi
! bet2(5,...)=phi=C/A, bet2(6,...)=eta=F/A-2L
! bet2(7,...)=Bc, bet2(8,...)=Bs,bet2(9,...)=Gc,bet2(10,...)=Gs
! bet2(11,...)=Hc, bet2(12,...)=Hs,bet2(13,...)=Ec,bet2(14,...)=Es
!
        do ilon = 1,ny
          if(nf <= 3 .or. nf >= 6)then
            bet2(ipa,idep,ilat,ilon) = beta(ipa,idep,ilat,ilon)*0.01*ppp + ppp
          else
            if(nf == 4)bet2(ipa,idep,ilat,ilon) = beta(ipa,idep,ilat,ilon)*0.01 + 1.
            if(nf == 5)bet2(ipa,idep,ilat,ilon) = - beta(ipa,idep,ilat,ilon)*0.01 + 1.
          endif
        enddo

       enddo
     enddo
   enddo
88 close(19)

!
! read anisotropic azimuthal parameters
!

!
! beta(ipa,idep,ilat,ilon) are sorted in (amplitude, phase)
! normalized, in percents: 100 G/L
!
  open(unit=15,file='DATA/Montagner_model/globpreman3sm01',status='old')

  do nf = 7,nfin,2
    ipa = nf
    ipa1 = ipa + 1
    do idep = 1,npar1
      il = idep + np1 - 1
      read(15,2011,end = 888) xinf,yinf,nx,ny,pxy
      read(15,2013,end = 888) pro(idep),ppp
      if(nf == 7) ppp = pari(2,il)
      if(nf == 9) ppp = pari(3,il)
      af = pari(6,il)*(pari(2,il) - 2.*pari(3,il))
      if(nf == 11) ppp = af
      if(nf == 13) ppp = (pari(4,il) + 1.)*pari(3,il)

      do ilat = 1,nx
        read(15,2015,end = 888) (alph(ilon,ilat),ilon = 1,ny)
      enddo

      do ilat=1,nx
        read(15,2015,end = 888) (ph(ilon,ilat),ilon = 1,ny)
      enddo

      do ilat = 1,nx
        do ilon = 1,ny
          angle = 2.*DEGREES_TO_RADIANS*ph(ilon,ilat)
          beta(ipa,idep,ilat,ilon) = alph(ilon,ilat)*ppp*0.01d0
          beta(ipa1,idep,ilat,ilon) = ph(ilon,ilat)
          bet2(ipa,idep,ilat,ilon) = alph(ilon,ilat)*dcos(angle)*ppp*0.01d0
          bet2(ipa1,idep,ilat,ilon) = alph(ilon,ilat)*dsin(angle)*ppp*0.01d0
        enddo
      enddo

    enddo
  enddo

888 close(15)

  do idep = 1,npar1
    do ilat = 1,nx
      do ilon = 1,ny

! rho
        beta(1,idep,ilat,ilon) = bet2(1,idep,ilat,ilon)

! A
        beta(2,idep,ilat,ilon) = bet2(2,idep,ilat,ilon)
        A=bet2(2,idep,ilat,ilon)

!  C
        beta(3,idep,ilat,ilon) = bet2(5,idep,ilat,ilon)*A

!  F
        A2L = A - 2.*bet2(3,idep,ilat,ilon)
        beta(4,idep,ilat,ilon) = bet2(6,idep,ilat,ilon)*A2L

!  L
        beta(5,idep,ilat,ilon) = bet2(3,idep,ilat,ilon)
        AL = bet2(3,idep,ilat,ilon)

!  N
        beta(6,idep,ilat,ilon) = bet2(4,idep,ilat,ilon)*AL

!  azimuthal terms
        do ipar = 7,14
          beta(ipar,idep,ilat,ilon) = bet2(ipar,idep,ilat,ilon)
        enddo

      enddo
    enddo
  enddo

 2011 format(2f4.0,2i3,f4.0)
 2013 format(f5.0,f8.4)
 2015 format(17f7.2)

 end subroutine read_aniso_mantle_model


!--------------------------------------------------------------------

  subroutine lecmod(nri,pari,ra)

  implicit none

! read the reference Earth model: rho, Vph, Vsv, XI, PHI, ETA
! array par(i,nlayer)
! output: array pari(ipar, nlayer): rho, A, L, xi-1, phi-1, eta-1

  integer i,j,k,ip,ifanis,idum1,idum2,idum3,nlayer,nout,neff,&
          nband,nri,minlay,moho,kiti
  double precision pari(14,47),qkappa(47),qshear(47),par(6,47)
  double precision epa(14,47),ra(47),dcori(47),ri(47)
  double precision corpar(21,47)
  double precision aa,an,al,af,ac,vpv,vph,vsv,vsh,rho,red,a2l
  character(len=80) null

     ifanis = 1
     nri = 47

     open(unit=13,file='DATA/Montagner_model/Adrem119',status='old')
     read(13,*,end = 77) nlayer,minlay,moho,nout,neff,nband,kiti,null

     if(kiti == 0) read(13,1000,end = 77) idum1
     read(13,1000,end = 77) idum2
     read(13,1000,end = 77) idum3

     do i = 1,nlayer
       read(13,1203,end = 77) ra(i),(par(k,i),k = 1,6),qshear(i),qkappa(i)
     enddo

     do i = 1,nlayer
       ri(i) = 0.001*ra(i)
     enddo

     do i = 1,nlayer
       rho = par(1,i)
       pari(1,i) = rho
!    A : pari(2,i)
       pari(2,i) = rho*(par(2,i)**2)
       aa = pari(2,i)
!    L : pari(3,i)
       pari(3,i) = rho*(par(3,i)**2)
       al = pari(3,i)
!    Xi : pari(4,i)= (N-L)/L
       an = al*par(4,i)
       pari(4,i) = 0.
       pari(4,i) = par(4,i) - 1.
!    Phi : pari(5,i)=(a-c)/a
       pari(5,i) = - par(5,i) + 1.
       ac = par(5,i)*aa
!    f : pari(4,i)
       af = par(6,i)*(aa - 2.*al)
       pari(6,i) = par(6,i)
       do ip = 7,14
         pari(ip,i) = 0.
       enddo
       vsv = 0.
       vsh = 0.
       if(al < 0.0001 .or. an < 0.0001) goto 12
       vsv = dsqrt(al/rho)
       vsh = dsqrt(an/rho)
 12    vpv = dsqrt(ac/rho)
       vph = dsqrt(aa/rho)
     enddo

  red = 1.
  do i = 1,nlayer
    read(13,1205,end = 77) (epa(j,i),j = 1,6),dcori(i)
    epa(7,i) = epa(2,i)
    epa(8,i) = epa(2,i)
    epa(9,i) = epa(3,i)
    epa(10,i) = epa(3,i)

    a2l = pari(2,i) - 2.*pari(3,i)
    epa(11,i) = epa(6,i)*a2l
    epa(12,i) = epa(6,i)*a2l
    epa(13,i) = epa(3,i)
    epa(14,i) = epa(3,i)

    do j = 1,14
      epa(j,i) = red*epa(j,i)
    enddo

    read(13,1208,end = 77) (corpar(j,i),j = 1,21)

  enddo

77 close(13)

 1000 format(20a4)
 1203 format(4x,f11.1,8d12.5)
 1205 format(15X,6e12.5,f11.1)
 1208 format(21f5.3)

  end subroutine lecmod

!--------------------------------------------------------------------

  subroutine rotate_aniso_tensor(theta,phi,d11,d12,d13,d16,d22,d23,d26,&
                           d33,d36,d44,d45,d55,d66,&
                           c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,&
                           c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  implicit none


  double precision theta,phi
  double precision c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                   c33,c34,c35,c36,c44,c45,c46,c55,c56,c66
  double precision d11,d12,d13,d16,d22,d23,d26,d33,d36,d44,d45,d55,d66
  double precision costheta,sintheta,cosphi,sinphi
  double precision costhetasq,sinthetasq,cosphisq,sinphisq
  double precision costwotheta,sintwotheta,costwophi,sintwophi
  double precision cosfourtheta,costhetafour,sinthetafour,cosphifour,sinphifour
  double precision sintwophisq,sintwothetasq

  costheta = dcos(theta)
  sintheta = dsin(theta)
  cosphi = dcos(phi)
  sinphi = dsin(phi)

  costhetasq = costheta * costheta
  sinthetasq = sintheta * sintheta
  cosphisq = cosphi * cosphi
  sinphisq = sinphi * sinphi

  costhetafour = costhetasq * costhetasq
  sinthetafour = sinthetasq * sinthetasq
  cosphifour = cosphisq * cosphisq
  sinphifour = sinphisq * sinphisq

  costwotheta = dcos(2.d0*theta)
  sintwotheta = dsin(2.d0*theta)
  costwophi = dcos(2.d0*phi)
  sintwophi = dsin(2.d0*phi)

  cosfourtheta = dcos(4.d0*theta)
  sintwothetasq = sintwotheta * sintwotheta
  sintwophisq = sintwophi * sintwophi

! recompute 21 anisotropic coefficients for Montagner's model using Mathematica

  c11 = cosphisq*(d12 + d23 + 2.*(d44 + d66) + &
        (d12 - d23 - 2.*d44 + 2.*d66)*costwotheta)*sinphisq - &
        2.*d26*sintwophi*costheta*sinphisq + d22*sinphifour - &
        2.*cosphisq*costheta*sintwophi* &
        (d16*costhetasq + (d36 + 2.*d45)*sinthetasq) + &
        cosphifour*(d11*costhetafour + &
        (d13 + 2.*d55)*sintwothetasq/2. + d33*sinthetafour)

  c12 = costwophi*costheta*(d16 - 2.*d26 + d36 + 2.*d45 + &
        (d16 - d36 - 2.*d45)*costwotheta)*sintwophi/2. + &
        (cosphifour + sinphifour)*(d12*costhetasq + d23*sinthetasq) + &
        sintwophisq*(d22 + d11*costhetafour - 4.*d44*sinthetasq + &
        d33*sinthetafour + costhetasq* &
        (-4.*d66 + 2.*(d13 + 2.*d55)*sinthetasq))/4.

  c13 = (cosphisq*(d11 + 6.*d13 + d33 - 4.*d55 - &
        (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta))/8. + &
        sinphisq*(d23*costhetasq + d12*sinthetasq) - &
        costheta*sintwophi* &
        (d36*costhetasq + (d16 - 2.*d45)*sinthetasq)

  c14 = -(sintheta*((cosphisq*costheta* &
        (d11 - d33 + 4.*d44 - 4.*d66 + &
        (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sinphi)/2. + &
        cosphi*(-d16 + d26 + d36 + (-d16 + d36 + 2.*d45)*costwotheta)* &
        sinphisq + (d12 - d23)*costheta*sinphi*sinphisq + &
        cosphi*cosphisq*((d16 - 2.*d45)*costhetasq + d36*sinthetasq)))

  c15 = (-(cosphi*cosphisq*costheta*(d11 - d33 + &
        (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta))/2. + &
        (cosphisq*(3.*d16 - d36 - 2.*d45 + &
        3.*(d16 - d36 - 2.*d45)*costwotheta)*sinphi)/2. - &
        (d12 - d23 - 2.*d44 + 2.*d66)*cosphi*costheta*sinphisq + &
        d26*sinphi*sinphisq)*sintheta

  c16 = (-3.*cosphisq*costheta*(d16 - 2.*d26 + d36 + 2.*d45 + &
        (d16 - d36 - 2.*d45)*costwotheta)*sinphisq)/2. - &
        d26*costheta*sinphifour + &
        sintwophi*sinphisq*(-d22 + (d12 + 2.*d66)*costhetasq + &
        (d23 + 2.*d44)*sinthetasq)/2. + &
        cosphifour*costheta*(d16*costhetasq + (d36 + 2.*d45)*sinthetasq) + &
        cosphisq*sintwophi*(d11*costhetafour - (d23 + 2.*d44)*sinthetasq + &
        d33*sinthetafour - costhetasq* &
        (d12 + 2.*d66 - 2.*(d13 + 2.*d55)*sinthetasq))/2.

  c22 = d22*cosphifour + 2.*d26*cosphisq*costheta*sintwophi + &
        (d12 + d23 + 2.*(d44 + d66) + &
        (d12 - d23 - 2.*d44 + 2.*d66)*costwotheta)*sintwophisq/4. + &
        costheta*(d16 + d36 + 2.*d45 + &
        (d16 - d36 - 2.*d45)*costwotheta)*sintwophi*sinphisq + &
        sinphifour*(d11*costhetafour + &
        (d13 + 2.*d55)*sintwothetasq/2. + d33*sinthetafour)

  c23 = costhetasq*(d23*cosphisq + &
        costheta*(d13*costheta*sinphisq + d36*sintwophi)) + &
        (d12*cosphisq - 2.*d45*costheta*sintwophi + &
        costheta*((d11 + d33 - 4.*d55)*costheta*sinphisq + &
        d16*sintwophi))*sinthetasq + d13*sinphisq*sinthetafour

  c24 = -((d26*cosphi*cosphisq + &
        (d12 - d23 - 2.*d44 + 2.*d66)*cosphisq*costheta* &
        sinphi + (cosphi*(3.*d16 - d36 - 2.*d45 + &
        3.*(d16 - d36 - 2.*d45)*costwotheta)*sinphisq)/2. + &
        (costheta*(d11 - d33 + (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)* &
        sinphi*sinphisq)/2.)*sintheta)

  c25 = sintheta*(-((d12 - d23)*cosphi*cosphisq*costheta) + &
        cosphisq*(-d16 + d26 + d36 + (-d16 + d36 + 2.*d45)*costwotheta)* &
        sinphi - (cosphi*costheta*(d11 - d33 + 4.*d44 - 4.*d66 + &
        (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sinphisq)/2. + &
        sinphi*sinphisq*((d16 - 2.*d45)*costhetasq + d36*sinthetasq))

  c26 = d26*cosphifour*costheta + (3.*cosphisq*costheta* &
        (d16 - 2.*d26 + d36 + 2.*d45 + (d16 - d36 - 2.*d45)*costwotheta)* &
        sinphisq)/2. + cosphisq*sintwophi* &
        (-d22 + (d12 + 2.*d66)*costhetasq + (d23 + 2.*d44)*sinthetasq)/2. - &
        costheta*sinphifour*(d16*costhetasq + (d36 + 2.*d45)*sinthetasq) + &
        sintwophi*sinphisq*(d11*costhetafour - (d23 + 2.*d44)*sinthetasq + &
        d33*sinthetafour - costhetasq* &
        (d12 + 2.*d66 - 2.*(d13 + 2.*d55)*sinthetasq))/2.

  c33 = d33*costhetafour + (d13 + 2.*d55)*sintwothetasq/2. + d11*sinthetafour

  c34 = -(sintheta*(-(costheta*(-d11 + d33 + &
        (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sinphi)/2. + &
        cosphi*((d36 + 2.*d45)*costhetasq + d16*sinthetasq)))

  c35 = sintheta*((cosphi*costheta* &
        (-d11 + d33 + (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta))/2. + &
        sinphi*((d36 + 2.*d45)*costhetasq + d16*sinthetasq))

  c36 = -((-d11 + 4.*d12 - 6.*d13 + 4.*d23 - d33 + 4.*d55 + &
        4.*(-d12 + d23)*costwotheta + &
        (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta)*sintwophi)/16. + &
        costwophi*costheta*(d36*costhetasq + (d16 - 2.*d45)*sinthetasq)

  c44 = costheta*costhetasq*(d55*costheta*sinphisq + d45*sintwophi) + &
        costheta*((d11 - 2.*d13 + d33 - 2.*d55)*costheta*sinphisq + &
        (d16 - d36 - d45)*sintwophi)*sinthetasq + &
        d55*sinphisq*sinthetafour + &
        cosphisq*(d44*costhetasq + d66*sinthetasq)

  c45 = (cosphisq*costheta*(d16 - d36 + &
        (-d16 + d36 + 2.*d45)*costwotheta))/2. + &
        ((d11 - 2.*d13 + d33 - 4.*(d44 - d55 + d66) + &
        4.*(-d44 + d66)*costwotheta - &
        (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta)*sintwophi)/16. - &
        costheta*sinphisq*(d45*costwotheta + (d16 - d36)*sinthetasq)

  c46 = sintheta*((d44 - d66)*cosphi*cosphisq*costheta + &
        cosphisq*(-d16 + d26 + d45 + (-d16 + d36 + 2.*d45)*costwotheta)* &
        sinphi - (cosphi*costheta* &
        (d11 -2.*d12 + 2.*d23 - d33 + 2.*d44 - 2.*d66 + &
        (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sinphisq)/2. +&
        sinphi*sinphisq*((d16 - d36 - d45)*costhetasq + d45*sinthetasq))

  c55 = (cosphisq*(d11 - 2.*d13 + d33 + 4.*d55 - &
        (d11 - 2.*d13 + d33 - 4.*d55)*cosfourtheta))/8. + &
        d44*costhetasq*sinphisq + &
        (d66*sinphisq + (d36 + d45)*costheta*sintwophi)*sinthetasq - &
        costheta*sintwophi*(d45*costhetasq + d16*sinthetasq)

  c56 = sintheta*(-(cosphisq*costheta* &
        (d11 - 2.*d12 + 2.*d23 - d33 + 2.*d44 - 2.*d66 + &
        (d11 - 2.*d13 + d33 - 4.*d55)*costwotheta)*sinphi)/2. - &
        cosphi*(-d16 + d26 + d45 + (-d16 + d36 + 2.*d45)*costwotheta)* &
        sinphisq + (d44 - d66)*costheta*sinphi*sinphisq + &
        cosphi*cosphisq*((-d16 + d36 + d45)*costhetasq - d45*sinthetasq))

  c66 = costheta*costwophi*(d16 - 2.*d26 + d36 + 2.*d45 + &
        (d16 - d36 - 2.*d45)*costwotheta)*sintwophi/2. + &
        (cosphifour + sinphifour)*(d66*costhetasq + d44*sinthetasq) + &
        sintwophisq*(d22 + d11*costhetafour - &
        2.*(d23 + d44)*sinthetasq + d33*sinthetafour - &
        2.*costhetasq*(d12 + d66 - (d13 + 2.*d55)*sinthetasq))/4.

  end subroutine rotate_aniso_tensor

