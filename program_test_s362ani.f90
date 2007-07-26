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

! minimum and maximum allowable topography of the Moho (km) for stretching purposes
  double precision, parameter :: TOPO_MAX = 5.0d0
  double precision, parameter :: TOPO_MIN = -5.0d0
  double precision, parameter :: TOPO_A = 0.5d0*(TOPO_MAX - TOPO_MIN)/R_EARTH_KM
  double precision, parameter :: TOPO_B = 0.5d0*(TOPO_MAX + TOPO_MIN)/R_EARTH_KM
  double precision, parameter :: RMOHO = 6346600.d0

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

! crustal_model_variables
  type crustal_model_variables
    sequence
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: thlr
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocp
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: velocs
    double precision, dimension(NKEYS_CRUST,NLAYERS_CRUST) :: dens
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

  integer i,iregion_code,ilat,ilon
  integer THREE_D_MODEL

  real(kind=4) xlat,xcolat,xlon,xdep,xrad
  real(kind=4) dvsh,dvsv,dvph,dvpv
  real(kind=4) topo410,topo650
  double precision scaleval
  double precision x,rho,vpv,vph,vsv,vsh,eta,Qmu,Qkappa
  double precision lat,lon,r,vp,vs,moho,r_moho,topomoho
  logical found_crust

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
!               THREE_D_MODEL_S362ANI_PREM,THREE_D_MODEL_S29EA)
!
!  call subshsv(xcolat,xlon,xrad,dvsh,dvsv,dvph,dvpv)
!  write(6,"('    dvsh      dvsv      dvph      dvpv    ')") 
!  write(6,"(6f10.5)") 100.0*dvsh,100.0*dvsv,100.0*dvph,100.0*dvpv
!
!  call subtopo(xcolat,xlon,topo410,topo650)
!  write(6,"('   topo410    topo650 ')") 
!  write(6,"(2f11.5)") topo410,topo650
!

! test moho stretching
  call read_crustal_model(CM_V)

  r = 0.9999d0
  do ilat = -89,89,5
    do ilon = -179,179,

      lat = dble(ilat)
      lon = dble(ilon)

! get Moho depth
      call crustal_model(lat,lon,r,vp,vs,rho,moho,found_crust,CM_V)
     if (.not. found_crust) stop 'failed to find Moho in add_topography_moho'
! Moho radius
      r_moho = 1.0d0 - moho

! Moho topography: negative for a depression, positive for an elevation
      topomoho = r_moho - RMOHO/R_EARTH

      topomoho = TOPO_B+TOPO_A*2.0d0*atan(10.0d0*(topomoho-TOPO_B)/TOPO_A)/PI

      print *, sngl(lat),sngl(lon),sngl(moho*R_EARTH_KM),sngl(r_moho*R_EARTH_KM-RMOHO/1000.0d0),sngl(topomoho*R_EARTH_KM)
    enddo
  enddo

  end program xtest_s362ani
