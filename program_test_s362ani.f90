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

  program xtest_s362ani

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

  integer i,iregion_code
  integer THREE_D_MODEL

  real xlat,xcolat,xlon,xdep,xrad
  real vshout,vsvout,vphout,vpvout,etaout,rhoout
  real topo410out,topo650out
  integer ifknowmodel
  double precision scaleval
  double precision x,rho,vpv,vph,vsv,vsh,eta,Qmu,Qkappa

!  call define_model_ref(.FALSE.,Mref_V)
!
!  scaleval=dsqrt(PI*GRAV*RHOAV)
!  do i=1,NR_REF
!    x = dble(i-1)/dble(NR_REF-1)
!    call model_ref(x,rho,vpv,vph,vsv,vsh,eta,Qkappa,Qmu,iregion_code,Mref_V)
!    print *,i,sngl(x*R_EARTH),sngl(rho*RHOAV),sngl(vpv*(R_EARTH*scaleval)),sngl(vph*(R_EARTH*scaleval)), &
!            sngl(vsv*(R_EARTH*scaleval)),sngl(vsh*(R_EARTH*scaleval)),sngl(eta),sngl(Qmu)
!
! find and check discontinuities
!    if(Mref_V%radius_ref(i) == Mref_V%radius_ref(i+1)) then
!      print *,i,Mref_V%radius_ref(i)
!    endif
!
!  enddo

  print *,"xlat: "
  read(5,*) xlat
  print *,"xlon: "
  read(5,*) xlon
  print *,"xdep: "
  read(5,*) xdep
  print *,"THREE_D_MODEL: "
  read(5,*) THREE_D_MODEL

  xcolat=90.0-xlat
  xrad=6371.0-xdep
  ifknowmodel=0
!  THREE_D_MODEL=THREE_D_MODEL_S362ANI

  call subshsv(xcolat,xlon,xrad,vshout,vsvout,vphout,vpvout,etaout,rhoout,ifknowmodel, &
               THREE_D_MODEL,THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
               THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
  write(6,"('    vsh       vsv       vph       vpv       eta       rho    ')") 
  write(6,"(6f10.5)") vshout,vsvout,vphout,vpvout,etaout,rhoout

  call subtopo(xcolat,xlon,topo410out,topo650out,ifknowmodel, &
               THREE_D_MODEL,THREE_D_MODEL_S362ANI,THREE_D_MODEL_S362WMANI, &
               THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
  write(6,"('   topo410    topo650 ')") 
  write(6,"(2f11.5)") topo410out,topo650out

  end program xtest_s362ani
