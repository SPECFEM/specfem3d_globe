!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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
! JP3D
!
! 3D japan Vp velocity model
!
! based on:
!
!          Program ----- veljp3d.f -----
!
!       This program is used to calculate 3-D P-wave velocity
!    distribution beneath the Japan Islands which is obtained
!    by a simultaneous inversion of arrival time data from local,
!    regional and teleseismic events.  For details, see "Deep
!    structure of the Japan subduction zone as derived from local,
!    regional, and teleseismic events" by Zhao, Hasegawa & Kanamori,
!    JGR, 99, 22313-22329, 1994.
!
!       The meaningful range of this model is as follows:
!        latitude : 32 - 45 N
!        longitude: 130-145 E
!        depth    : 0  - 500 km
!
!                            Dapeng Zhao
!                            Dept. of Earth & Planet. Sci
!                            Washington University
!                            St. Louis, MO 63130
!                            U.S.A.
!                            dapeng@izu.wustl.edu
!
!
!         Last Time Modified by Min Chen, Caltech, 03/14/2008
!
!--------------------------------------------------------------------------------------------------

  module model_jp3d_par

  ! Japan 3D model (Zhao, 1994) constants
  integer, parameter :: MPA=42,MRA=48,MHA=21,MPB=42,MRB=48,MHB=18
  integer, parameter :: MKA=2101,MKB=2101

  ! JP3DM_V arrays: vmod3d
  double precision,dimension(:),allocatable :: JP3DM_PNA,JP3DM_RNA,JP3DM_HNA, &
    JP3DM_PNB,JP3DM_RNB,JP3DM_HNB
  double precision,dimension(:,:,:),allocatable :: JP3DM_VELAP,JP3DM_VELBP
  ! discon
  double precision,dimension(:),allocatable :: JP3DM_PN,JP3DM_RRN
  double precision,dimension(:,:),allocatable :: JP3DM_DEPA,JP3DM_DEPB,JP3DM_DEPC
  ! locate
  double precision :: JP3DM_PLA
  double precision :: JP3DM_RLA
  double precision :: JP3DM_HLA
  double precision :: JP3DM_PLB
  double precision :: JP3DM_RLB
  double precision :: JP3DM_HLB
  ! weight
  double precision,dimension(:),allocatable :: JP3DM_WV
  ! prhfd
  double precision :: JP3DM_P
  double precision :: JP3DM_R
  double precision :: JP3DM_H
  double precision :: JP3DM_PF
  double precision :: JP3DM_RF
  double precision :: JP3DM_HF
  double precision :: JP3DM_PF1
  double precision :: JP3DM_RF1
  double precision :: JP3DM_HF1
  double precision :: JP3DM_PD
  double precision :: JP3DM_RD
  double precision :: JP3DM_HD
  ! jpmodv
  double precision,dimension(:),allocatable :: JP3DM_VP,JP3DM_VS,JP3DM_RA,JP3DM_DEPJ
  ! locate integers
  integer,dimension(:),allocatable :: JP3DM_IPLOCA,JP3DM_IRLOCA,JP3DM_IHLOCA, &
    JP3DM_IPLOCB,JP3DM_IRLOCB,JP3DM_IHLOCB

  ! vmod3D integers
  integer :: JP3DM_NPA
  integer :: JP3DM_NRA
  integer :: JP3DM_NHA
  integer :: JP3DM_NPB
  integer :: JP3DM_NRB
  integer :: JP3DM_NHB
  ! weight integers
  integer :: JP3DM_IP
  integer :: JP3DM_JP
  integer :: JP3DM_KP
  integer :: JP3DM_IP1
  integer :: JP3DM_JP1
  integer :: JP3DM_KP1

  end module model_jp3d_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_jp3d_broadcast(myrank)

! standard routine to setup model

  use model_jp3d_par

  implicit none

  include "constants.h"
  ! standard include of the MPI library
  include 'mpif.h'

  integer :: myrank
  integer :: ier

  ! allocates arrays
! model_jp3d_variables
  allocate(JP3DM_PNA(MPA),JP3DM_RNA(MRA),JP3DM_HNA(MHA), &
          JP3DM_PNB(MPB),JP3DM_RNB(MRB),JP3DM_HNB(MHB), &
          JP3DM_VELAP(MPA,MRA,MHA), &
          JP3DM_VELBP(MPB,MRB,MHB), &
          JP3DM_PN(51),JP3DM_RRN(63), &
          JP3DM_DEPA(51,63),JP3DM_DEPB(51,63),JP3DM_DEPC(51,63), &
          JP3DM_WV(8),JP3DM_VP(29),JP3DM_VS(29), &
          JP3DM_RA(29),JP3DM_DEPJ(29), &
          JP3DM_IPLOCA(MKA),JP3DM_IRLOCA(MKA), &
          JP3DM_IHLOCA(MKA),JP3DM_IPLOCB(MKB), &
          JP3DM_IRLOCB(MKB),JP3DM_IHLOCB(MKB), &
          stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating JP3D arrays')

  ! master reads in values
  if(myrank == 0) call read_jp3d_iso_zhao_model()

  ! JP3DM_V
  call MPI_BCAST(JP3DM_NPA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_NRA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_NHA,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_NPB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_NRB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_NHB,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PNA,MPA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RNA,MRA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HNA,MHA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PNB,MPB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RNB,MRB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HNB,MHB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_VELAP,MPA*MRA*MHA,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_VELBP,MPB*MRB*MHB,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PN,51,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RRN,63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_DEPA,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_DEPB,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_DEPC,51*63,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IPLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IRLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IHLOCA,MKA,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IPLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IRLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IHLOCB,MKB,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HLA,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HLB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_JP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_KP,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_IP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_JP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_KP1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_WV,8,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_P,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_R,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_H,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HF,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HF1,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_PD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_HD,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_VP,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_VS,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_RA,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(JP3DM_DEPJ,29,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)


  end subroutine model_jp3d_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_jp3d_iso_zhao_model()

  use model_jp3d_par

  implicit none

  include "constants.h"

  OPEN(2,FILE="DATA/Zhao_JP_model/m3d1341")
  OPEN(3,FILE="DATA/Zhao_JP_model/datadis")

  CALL INPUTJP()
  CALL INPUT1()
  CALL INPUT2()

  end subroutine read_jp3d_iso_zhao_model

!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_jp3d_iso_zhao(radius,theta,phi,vp,vs,dvp,dvs,rho,found_crust)

  use model_jp3d_par

  implicit none

  include "constants.h"

  logical found_crust
  double precision :: radius,theta,phi,vp,vs,dvs,dvp,rho
  double precision :: PE,RE,HE,H1,H2,H3,scaleval
  integer :: LAY


  found_crust = .false.

  PE = theta
  RE = phi
  HE = (ONE - radius)*R_EARTH_KM
!  calculate depths of the Conrad, the Moho and
!  the plate boundary beneath the location (PHI,RAM)
  CALL HLAY(PE,RE,H1,1)
  CALL HLAY(PE,RE,H2,2)
  CALL HLAY(PE,RE,H3,3)
!   when LAY = 1, the focus is in the upper crust;
!   when LAY = 2, the focus is in the lower crust;
!   when LAY = 3, the focus is in the mantle wedge;
!   when LAY = 4, the focus is beneath the plate boundary.
  IF(HE.LE.H1)                   THEN
     LAY = 1
     found_crust = .true.
  ELSE IF(HE.GT.H1.AND.HE.LE.H2) THEN
     LAY = 2
     found_crust = .true.
  ELSE IF(HE.GT.H2.AND.HE.LE.H3) THEN
     LAY = 3
  ELSE
     LAY = 4
  END IF

  CALL VEL1D(HE,vp,LAY,1)
  CALL VEL1D(HE,vs,LAY,2)
  CALL VEL3(PE,RE,HE,dvp,LAY)

  dvp = 0.01d0*dvp
  dvs = 1.5d0*dvp
  vp = vp*(1.0d0+dvp)
  vs = vs*(1.0d0+dvs)

! determine rho
  if(LAY .eq. 1) then
     rho=2.6
  endif
  if(LAY .eq. 2) then
     rho=2.9
  endif
  if(LAY .GT. 2) then
     rho=3.3+(vs-4.4)*0.66667
  endif
! non-dimensionalize
! time scaling (s^{-1}) is done with scaleval
  scaleval=dsqrt(PI*GRAV*RHOAV)
  rho=rho*1000.0d0/RHOAV
  vp=vp*1000.0d0/(R_EARTH*scaleval)
  vs=vs*1000.0d0/(R_EARTH*scaleval)

  END subroutine model_jp3d_iso_zhao

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE INPUT1()

  use model_jp3d_par

  implicit none

  include "constants.h"

100     FORMAT(3I3)

  READ(2,100)  JP3DM_NPA,JP3DM_NRA,JP3DM_NHA

  CALL PUT1(JP3DM_NPA,JP3DM_NRA,JP3DM_NHA,JP3DM_PNA,JP3DM_RNA,JP3DM_HNA,JP3DM_VELAP)

  READ(2,100)  JP3DM_NPB,JP3DM_NRB,JP3DM_NHB

  CALL PUT1(JP3DM_NPB,JP3DM_NRB,JP3DM_NHB,JP3DM_PNB,JP3DM_RNB,JP3DM_HNB,JP3DM_VELBP)
  CALL BLDMAP()

  RETURN
  END SUBROUTINE INPUT1

!
!-------------------------------------------------------------------------------------------------
!
  SUBROUTINE PUT1(NPX,NRX,NHX,PNX,RNX,HNX,VELXP)

  implicit none

  integer :: NPX,NRX,NHX,K,I,J
  double precision ::  VELXP(NPX,NRX,NHX), &
                PNX(NPX),RNX(NRX),HNX(NHX)

  READ(2,110) (PNX(I),I=1,NPX)
  READ(2,110) (RNX(I),I=1,NRX)
  READ(2,120) (HNX(I),I=1,NHX)
      DO K = 1,NHX
         DO I = 1,NPX
            READ(2,140) (VELXP(I,J,K),J=1,NRX)
110         FORMAT(6(9F7.2/))
120         FORMAT(3(8F7.2/))
140         FORMAT(4(14F5.2/))
         enddo
      enddo
  END SUBROUTINE PUT1

!
!---------------------------------------------------------------------------------------------
!
  SUBROUTINE INPUT2()

  use model_jp3d_par

  implicit none

  include "constants.h"

  integer :: NP,NNR,I,J

  READ(3,100)  NP,NNR
  READ(3,110) (JP3DM_PN(I),I=1,NP)
  READ(3,120) (JP3DM_RRN(I),I=1,NNR)
  DO 1  I = NP,1,-1
      READ(3,130) (JP3DM_DEPA(I,J),J=1,NNR)
1     CONTINUE
  DO 2  I = NP,1,-1
      READ(3,130) (JP3DM_DEPB(I,J),J=1,NNR)
2     CONTINUE
  DO 3  I = NP,1,-1
      READ(3,130) (JP3DM_DEPC(I,J),J=1,NNR)
3     CONTINUE

100   FORMAT(2I6)
110   FORMAT(5(10F7.2/),F7.2)
120   FORMAT(6(10F7.2/),3F7.2)
130   FORMAT(6(10F7.1/),3F7.1)

  RETURN
  END

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE BLDMAP()

  use model_jp3d_par

  implicit none

  include "constants.h"

  CALL LOCX(JP3DM_PNA,JP3DM_RNA,JP3DM_HNA,JP3DM_NPA,JP3DM_NRA,JP3DM_NHA,MKA, &
           JP3DM_PLA,JP3DM_RLA,JP3DM_HLA,JP3DM_IPLOCA,JP3DM_IRLOCA,JP3DM_IHLOCA)
  CALL LOCX(JP3DM_PNB,JP3DM_RNB,JP3DM_HNB,JP3DM_NPB,JP3DM_NRB,JP3DM_NHB,MKB, &
           JP3DM_PLB,JP3DM_RLB,JP3DM_HLB,JP3DM_IPLOCB,JP3DM_IRLOCB,JP3DM_IHLOCB)

  RETURN
  END

!
!-------------------------------------------------------------------------------------------------
!
      SUBROUTINE LOCX(PNX,RNX,HNX,NPX,NRX,NHX,MKX, &
                 PLX,RLX,HLX,IPLOCX,IRLOCX,IHLOCX)
     integer ::  NPX,NRX,NHX,MKX,IPLOCX(MKX),IRLOCX(MKX),IHLOCX(MKX)
     integer ::  IPMAX,IP,IP1,IRMAX,IR,IR1,IH1,IH,IHMAX,I
     double precision :: PNX(NPX),RNX(NRX),HNX(NHX)
     double precision :: PLX,RLX,HLX,PNOW,RNOW,HNOW
      PLX      = 1.0-PNX(1)*100.0
      IPMAX    = IDNINT(PNX(NPX)*100.0+PLX)
      IP       = 1
      DO 10 I  = 1,IPMAX
      IP1      = IP+1
      PNOW     = (FLOAT(I)-PLX)/100.0
      IF(PNOW.GE.PNX(IP1))   IP = IP1
      IPLOCX(I)= IP
10    CONTINUE
      RLX      = 1.0-RNX(1)*100.0
      IRMAX    = IDNINT(RNX(NRX)*100.0+RLX)
      IR       = 1
      DO 20 I  = 1,IRMAX
      IR1      = IR+1
      RNOW     = (FLOAT(I)-RLX)/100.0
      IF(RNOW.GE.RNX(IR1))   IR = IR1
      IRLOCX(I)= IR
20    CONTINUE
      HLX      = 1.0-HNX(1)
      IHMAX    = IDNINT(HNX(NHX)+HLX)
      IH       = 1
      DO 30 I  = 1,IHMAX
      IH1      = IH+1
      HNOW     = FLOAT(I)-HLX
      IF(HNOW.GE.HNX(IH1))   IH = IH1
      IHLOCX(I)= IH
30    CONTINUE
      RETURN
      END

!
!-------------------------------------------------------------------------------------------
!

  SUBROUTINE VEL3(PE,RE,HE,V,LAY)

  use model_jp3d_par

  implicit none

  include "constants.h"

  double precision :: PE,RE,HE,V

  integer :: LAY

  JP3DM_P     = 90.0-PE/DEGREES_TO_RADIANS
  JP3DM_R     = RE/DEGREES_TO_RADIANS
  JP3DM_H     = HE
  IF(LAY.LE.3)       THEN
     CALL PRHF(JP3DM_IPLOCA,JP3DM_IRLOCA,JP3DM_IHLOCA,JP3DM_PLA,JP3DM_RLA,JP3DM_HLA, &
          JP3DM_PNA,JP3DM_RNA,JP3DM_HNA,MPA,MRA,MHA,MKA)
  ELSE IF(LAY.EQ.4)  THEN
     CALL PRHF(JP3DM_IPLOCB,JP3DM_IRLOCB,JP3DM_IHLOCB,JP3DM_PLB,JP3DM_RLB,JP3DM_HLB, &
          JP3DM_PNB,JP3DM_RNB,JP3DM_HNB,MPB,MRB,MHB,MKB)
  ELSE
  END IF

  JP3DM_WV(1) = JP3DM_PF1*JP3DM_RF1*JP3DM_HF1
  JP3DM_WV(2) = JP3DM_PF*JP3DM_RF1*JP3DM_HF1
  JP3DM_WV(3) = JP3DM_PF1*JP3DM_RF*JP3DM_HF1
  JP3DM_WV(4) = JP3DM_PF*JP3DM_RF*JP3DM_HF1
  JP3DM_WV(5) = JP3DM_PF1*JP3DM_RF1*JP3DM_HF
  JP3DM_WV(6) = JP3DM_PF*JP3DM_RF1*JP3DM_HF
  JP3DM_WV(7) = JP3DM_PF1*JP3DM_RF*JP3DM_HF
  JP3DM_WV(8) = JP3DM_PF*JP3DM_RF*JP3DM_HF

  !   calculate velocity
  IF(LAY.LE.3)      THEN
     CALL VABPS(MPA,MRA,MHA,JP3DM_VELAP,V)
  ELSE IF(LAY.EQ.4) THEN
     CALL VABPS(MPB,MRB,MHB,JP3DM_VELBP,V)
  ELSE
  END IF

  RETURN
  END SUBROUTINE VEL3

!
!---------------------------------------------------------------------------------------
!

  SUBROUTINE VABPS(MP,MR,MH,V,VEL)

  use model_jp3d_par

  implicit none

  include "constants.h"

  double precision :: VEL
  integer :: MP,MR,MH
  double precision :: V(MP,MR,MH)

  VEL = JP3DM_WV(1)*V(JP3DM_IP,JP3DM_JP,JP3DM_KP)  + JP3DM_WV(2)*V(JP3DM_IP1,JP3DM_JP,JP3DM_KP) &
      + JP3DM_WV(3)*V(JP3DM_IP,JP3DM_JP1,JP3DM_KP) + JP3DM_WV(4)*V(JP3DM_IP1,JP3DM_JP1,JP3DM_KP) &
      + JP3DM_WV(5)*V(JP3DM_IP,JP3DM_JP,JP3DM_KP1) + JP3DM_WV(6)*V(JP3DM_IP1,JP3DM_JP,JP3DM_KP1) &
      + JP3DM_WV(7)*V(JP3DM_IP,JP3DM_JP1,JP3DM_KP1)+ JP3DM_WV(8)*V(JP3DM_IP1,JP3DM_JP1,JP3DM_KP1)

  RETURN
  END

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE INTMAP(R,IRLOC,NNR,RL,IR)

  implicit none
  integer :: NNR,IRLOC(NNR),IS,IR
  double precision :: R,RL

  IS      = IDNINT(R+RL)
  IR      = IRLOC(IS)

  RETURN
  END

!
!------------------------------------------------------------------------------------------------
!

  SUBROUTINE PRHF(IPLOCX,IRLOCX,IHLOCX,PLX,RLX,HLX, &
                      PNX,RNX,HNX,MPX,MRX,MHX,MKX)

  use model_jp3d_par

  implicit none

  include "constants.h"

  integer :: MPX,MRX,MHX,MKX
  integer ::  IPLOCX(MKX),IRLOCX(MKX),IHLOCX(MKX)
  double precision :: PNX(MPX),RNX(MRX),HNX(MHX)
  double precision :: PLX,RLX,HLX

  CALL LIMIT(PNX(1),PNX(MPX),JP3DM_P)
  CALL LIMIT(RNX(1),RNX(MRX),JP3DM_R)
  CALL LIMIT(HNX(1),HNX(MHX),JP3DM_H)
  CALL INTMAP(JP3DM_P*100.0,IPLOCX,MKX,PLX,JP3DM_IP)
  CALL INTMAP(JP3DM_R*100.0,IRLOCX,MKX,RLX,JP3DM_JP)
  CALL INTMAP(JP3DM_H,IHLOCX,MKX,HLX,JP3DM_KP)
  JP3DM_IP1   = JP3DM_IP+1
  JP3DM_JP1   = JP3DM_JP+1
  JP3DM_KP1   = JP3DM_KP+1
  JP3DM_PD    = PNX(JP3DM_IP1)-PNX(JP3DM_IP)
  JP3DM_RD    = RNX(JP3DM_JP1)-RNX(JP3DM_JP)
  JP3DM_HD    = HNX(JP3DM_KP1)-HNX(JP3DM_KP)
  JP3DM_PF    = (JP3DM_P-PNX(JP3DM_IP))/JP3DM_PD
  JP3DM_RF    = (JP3DM_R-RNX(JP3DM_JP))/JP3DM_RD
  JP3DM_HF    = (JP3DM_H-HNX(JP3DM_KP))/JP3DM_HD
  JP3DM_PF1   = 1.0-JP3DM_PF
  JP3DM_RF1   = 1.0-JP3DM_RF
  JP3DM_HF1   = 1.0-JP3DM_HF
  RETURN
  END

!
!----------------------------------------------------------------------------------------------
!

  SUBROUTINE HLAY(PE,RE,HE,IJK)

  use model_jp3d_par

  implicit none

  include "constants.h"

  double precision :: PE,RE,HE,WV1,WV2,WV3,WV4,P,R,PF,RF,PF1,RF1
  integer :: IJK,J,J1,I,I1

  P = 90.0-PE/DEGREES_TO_RADIANS
  R = RE/DEGREES_TO_RADIANS

  CALL LIMIT(JP3DM_PN(1),JP3DM_PN(51),P)
  CALL LIMIT(JP3DM_RRN(1),JP3DM_RRN(63),R)

        DO 1 I = 1,50
           I1     = I+1
           IF(P.GE.JP3DM_PN(I).AND.P.LT.JP3DM_PN(I1)) GO TO 11
1          CONTINUE
11         CONTINUE
           DO 2 J = 1,62
              J1     = J+1
              IF(R.GE.JP3DM_RRN(J).AND.R.LT.JP3DM_RRN(J1)) GO TO 22
2             CONTINUE
22            CONTINUE
              PF    = (P-JP3DM_PN(I))/(JP3DM_PN(I1)-JP3DM_PN(I))
              RF    = (R-JP3DM_RRN(J))/(JP3DM_RRN(J1)-JP3DM_RRN(J))
              PF1   = 1.0-PF
              RF1   = 1.0-RF
              WV1   = PF1*RF1
              WV2   = PF*RF1
              WV3   = PF1*RF
              WV4   = PF*RF
              IF(IJK.EQ.1)       THEN
                 HE  = WV1*JP3DM_DEPA(I,J)  + WV2*JP3DM_DEPA(I1,J) &
                      + WV3*JP3DM_DEPA(I,J1) + WV4*JP3DM_DEPA(I1,J1)
              ELSE IF(IJK.EQ.2)  THEN
                 HE  = WV1*JP3DM_DEPB(I,J)  + WV2*JP3DM_DEPB(I1,J) &
                      + WV3*JP3DM_DEPB(I,J1) + WV4*JP3DM_DEPB(I1,J1)
              ELSE IF(IJK.EQ.3)  THEN
                 HE  = WV1*JP3DM_DEPC(I,J)  + WV2*JP3DM_DEPC(I1,J) &
                      + WV3*JP3DM_DEPC(I,J1) + WV4*JP3DM_DEPC(I1,J1)
              ELSE
              END IF
  RETURN
  END SUBROUTINE HLAY

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE LIMIT(C1,C2,C)

  implicit none
  double precision :: A1,A2,C1,C2,C

  A1    = dmin1(C1,C2)
  A2    = dmax1(C1,C2)
  IF(C.LT.A1)   C = A1
  IF(C.GT.A2)   C = A2

  END SUBROUTINE LIMIT

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE VEL1D(HE,V,LAY,IPS)

  use model_jp3d_par

  implicit none

  include "constants.h"

  integer :: IPS,LAY
  double precision :: HE,V,VM,HM

  IF(LAY.EQ.1)      THEN
    V    = 6.0
    IF(IPS.EQ.2)    V = 3.5
  ELSE IF(LAY.EQ.2) THEN
    V    = 6.7
    IF(IPS.EQ.2)    V = 3.8
  ELSE IF(LAY.GE.3) THEN
    HM   = 40.0
    IF(HE.LT.HM)    THEN
      CALL JPMODEL(IPS,HM,VM)
      V  = VM-(HM-HE)*0.003
    ELSE
      CALL JPMODEL(IPS,HE,V)
    END IF
  ELSE
  END IF

  RETURN
  END

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE INPUTJP()

  use model_jp3d_par

  implicit none

  include "constants.h"

  double precision :: VP1(29),VS1(29),RA1(29)
  integer :: L

  DATA VP1/7.75, 7.94, 8.13, 8.33, 8.54, 8.75, 8.97, &
           9.50, 9.91,10.26,10.55,10.99,11.29,11.50, &
          11.67,11.85,12.03,12.20,12.37,12.54,12.71, &
          12.87,13.02,13.16,13.32,13.46,13.60,13.64,13.64/
  DATA VS1/4.353,4.444,4.539,4.638,4.741,4.850,4.962, &
           5.227,5.463,5.670,5.850,6.125,6.295,6.395, &
           6.483,6.564,6.637,6.706,6.770,6.833,6.893, &
           6.953,7.012,7.074,7.137,7.199,7.258,7.314,7.304/
  DATA RA1/1.00,0.99,0.98,0.97,0.96,0.95,0.94,0.93, &
           0.92,0.91,0.90,0.88,0.86,0.84,0.82,0.80, &
           0.78,0.76,0.74,0.72,0.70,0.68,0.66,0.64, &
           0.62,0.60,0.58,0.56,0.55/

  DO 1 L  = 1,29
    JP3DM_VP(L)   = VP1(L)
    JP3DM_VS(L)   = VS1(L)
    JP3DM_RA(L)   = RA1(L)
    JP3DM_DEPJ(L) = 40.0+6325.59*(1.0-RA1(L))
1     CONTINUE

  RETURN
  END

!
!-------------------------------------------------------------------------------------------------
!

  SUBROUTINE JPMODEL(IPS,H,V)

  use model_jp3d_par

  implicit none

  include "constants.h"

  integer :: IPS,K,K1
  double precision :: H1,H2,H12,H,V

  DO 2 K = 1,28
      K1     = K+1
      H1     = JP3DM_DEPJ(K)
      H2     = JP3DM_DEPJ(K1)
      IF(H.GE.H1.AND.H.LT.H2) GO TO 3
2     CONTINUE
3     CONTINUE

  H12    = (H-H1)/(H2-H1)
  IF(IPS.EQ.1)  THEN
     V   = (JP3DM_VP(K1)-JP3DM_VP(K))*H12+JP3DM_VP(K)
  ELSE
     V   = (JP3DM_VS(K1)-JP3DM_VS(K))*H12+JP3DM_VS(K)
  END IF

  RETURN
  END

