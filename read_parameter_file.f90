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

  subroutine read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
        NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
        NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB,NER_DOUBLING_OUTER_CORE, &
        NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP, &
        DT,TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,OCEANS,ELLIPTICITY, &
        GRAVITY,ONE_CRUST,ATTENUATION, &
        ROTATION,THREE_D,TOPOGRAPHY,LOCAL_PATH,NSOURCES,&
        MOVIE_SURFACE,MOVIE_VOLUME,NMOVIE,HDUR_MIN_MOVIES, &
        NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC)

  implicit none

  include "constants.h"

  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_DOUBLING_OUTER_CORE, &
             NER_CRUST,NER_220_MOHO,NER_400_220, &
             NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
             NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
             NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP

  double precision DT,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC,HDUR_MIN_MOVIES

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY, &
             GRAVITY,ONE_CRUST,ROTATION, &
             THREE_D,TOPOGRAPHY,ATTENUATION,OCEANS, &
             MOVIE_SURFACE,MOVIE_VOLUME
  integer NSOURCES,NMOVIE,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB

  character(len=150) LOCAL_PATH

  integer i

! first 27 characters of each line in the file are a comment
  character(len=27) junk

  open(unit=IIN,file='DATA/Par_file',status='old')

! ignore header
  do i=1,11
    read(IIN,*)
  enddo

  read(IIN,1) junk,NEX_XI
  read(IIN,1) junk,NEX_ETA
  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NPROC_XI
  read(IIN,1) junk,NPROC_ETA
  read(IIN,*)
  read(IIN,*)
  read(IIN,2) junk,DT
  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NSTEP
  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,OCEANS
  read(IIN,3) junk,ELLIPTICITY
  read(IIN,3) junk,TOPOGRAPHY
  read(IIN,3) junk,THREE_D
  read(IIN,3) junk,GRAVITY
  read(IIN,3) junk,ROTATION
  read(IIN,3) junk,TRANSVERSE_ISOTROPY
  read(IIN,3) junk,ANISOTROPIC_MANTLE
  read(IIN,3) junk,ANISOTROPIC_INNER_CORE
  read(IIN,3) junk,CRUSTAL
  read(IIN,3) junk,ONE_CRUST
  read(IIN,3) junk,ATTENUATION
  read(IIN,1) junk,MIN_ATTENUATION_PERIOD
  read(IIN,1) junk,MAX_ATTENUATION_PERIOD
  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NER_CRUST
  read(IIN,1) junk,NER_220_MOHO
  read(IIN,1) junk,NER_400_220
  read(IIN,1) junk,NER_600_400
  read(IIN,1) junk,NER_670_600
  read(IIN,1) junk,NER_771_670
  read(IIN,1) junk,NER_TOPDDOUBLEPRIME_771
  read(IIN,1) junk,NER_CMB_TOPDDOUBLEPRIME
  read(IIN,2) junk,RATIO_TOP_DBL_OC
  read(IIN,2) junk,RATIO_BOTTOM_DBL_OC
  read(IIN,1) junk,NER_TOPDBL_CMB
  read(IIN,1) junk,NER_ICB_BOTTOMDBL
  read(IIN,1) junk,NER_TOP_CENTRAL_CUBE_ICB

! scale radial mesh parameters according to definitions used in mesher
! in order to implement mesh doubling
  NER_220_MOHO = NER_220_MOHO * 2
  NER_400_220 = NER_400_220 * 2
  NER_600_400 = NER_600_400 * 2
  NER_670_600 = NER_670_600 * 2

  NER_771_670 = NER_771_670 * 4
  NER_TOPDDOUBLEPRIME_771 = NER_TOPDDOUBLEPRIME_771 * 4
  NER_CMB_TOPDDOUBLEPRIME = NER_CMB_TOPDDOUBLEPRIME * 4
  NER_TOPDBL_CMB = NER_TOPDBL_CMB * 4
  NER_ICB_BOTTOMDBL = NER_ICB_BOTTOMDBL * 4
  NER_TOP_CENTRAL_CUBE_ICB = NER_TOP_CENTRAL_CUBE_ICB * 4

  NER_ICB_CMB = NER_ICB_BOTTOMDBL + NER_BOTTOMDBL_TOPDBL + NER_TOPDBL_CMB
  NER_DOUBLING_OUTER_CORE = NER_TOP_CENTRAL_CUBE_ICB + NER_ICB_BOTTOMDBL + NER_BOTTOMDBL_TOPDBL

  read(IIN,*)
  read(IIN,*)
  read(IIN,4) junk,LOCAL_PATH

! ignore name of machine file (used by scripts but not by mesher nor solver)
  read(IIN,*)
  read(IIN,*)
  read(IIN,*)

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NSEIS

  read(IIN,*)
  read(IIN,*)
  read(IIN,1) junk,NSOURCES

  read(IIN,*)
  read(IIN,*)
  read(IIN,3) junk,MOVIE_SURFACE
  read(IIN,3) junk,MOVIE_VOLUME
  read(IIN,1) junk,NMOVIE
  read(IIN,2) junk,HDUR_MIN_MOVIES

! close parameter file
  close(IIN)

! formats

 1 format(a,i8)
 2 format(a,f12.5)
 3 format(a,l8)
 4 format(a,a)

  end subroutine read_parameter_file

