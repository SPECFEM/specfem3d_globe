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

! create file OUTPUT_FILES/values_from_mesher.h based upon DATA/Par_file
! in order to compile the solver with the right array sizes

  program create_header_file

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
             NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
             NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
             NEX_ETA,NEX_XI,NER_DOUBLING_OUTER_CORE, &
             NPROC_ETA,NPROC_XI,NSEIS,NSTEP

  double precision DT

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,ELLIPTICITY, &
             GRAVITY,ONE_CRUST,ROTATION, &
             THREE_D,TOPOGRAPHY,ATTENUATION,OCEANS, &
! BS
!            MOVIE_SURFACE,MOVIE_VOLUME
             MOVIE_SURFACE,MOVIE_VOLUME,ATTENUATION_3D
! BS END
  integer NSOURCES,NMOVIE,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB
  double precision RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC,HDUR_MIN_MOVIES

  character(len=150) LOCAL_PATH

! parameters deduced from parameters read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

! this for all the regions
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
               NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
               NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
               NSPEC2D_BOTTOM,NSPEC2D_TOP, &
               NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
               NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
               nglob_AB,nglob_AC,nglob_BC

! topology of the elements
  integer iaddx(NGNOD)
  integer iaddy(NGNOD)
  integer iaddz(NGNOD)

! auxiliary variables to generate the mesh
  integer ix,iy,ir,ir1,ir2,dir
  integer ix1,ix2,dix,iy1,iy2,diy
  integer iax,iay,iar
  integer ichunk,isubregion,nsubregions,doubling_index

  integer ispec_aniso
  integer npx,npy

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'creating file OUTPUT_FILES/values_from_mesher.h to compile solver with correct values'

! read the parameter file
  call read_parameter_file(MIN_ATTENUATION_PERIOD,MAX_ATTENUATION_PERIOD,NER_CRUST,NER_220_MOHO,NER_400_220, &
        NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
        NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB,NER_DOUBLING_OUTER_CORE, &
        NEX_ETA,NEX_XI,NPROC_ETA,NPROC_XI,NSEIS,NSTEP, &
        DT,TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE,CRUSTAL,OCEANS,ELLIPTICITY, &
        GRAVITY,ONE_CRUST,ATTENUATION, &
        ROTATION,THREE_D,TOPOGRAPHY,LOCAL_PATH,NSOURCES, &
        MOVIE_SURFACE,MOVIE_VOLUME,NMOVIE,HDUR_MIN_MOVIES, &
! BS
!       NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC)
        NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB,RATIO_BOTTOM_DBL_OC,RATIO_TOP_DBL_OC,ATTENUATION_3D)
! BS END

! compute other parameters based upon values read
  call compute_parameters(NER_CRUST,NER_220_MOHO,NER_400_220, &
      NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
      NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
      NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB, &
      NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA, &
      NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NSPEC_AB,NSPEC_AC,NSPEC_BC, &
      NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
      NGLOB_AB,NGLOB_AC,NGLOB_BC,NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB)

! count number of elements in the largest chunk
  ichunk = CHUNK_BC

! number of elements in each slice in a given chunk
  npx = 2*NEX_PER_PROC_XI
  npy = 2*NEX_PER_PROC_ETA

! generate the elements in all the regions of the mesh
  ispec_aniso = 0

! define number of subregions in current region
  nsubregions = 32

  do isubregion = 1,nsubregions

! define shape of elements in current subregion
    call define_subregions_crust_mantle(isubregion,ichunk,iaddx,iaddy,iaddz, &
      ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
      doubling_index,npx,npy,NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400, &
      NER_400_220,NER_220_MOHO,NER_CRUST)

  do ir = ir1,ir2,dir
    do iy = iy1,iy2,diy
      do ix = ix1,ix2,dix

! count number of anisotropic elements
        if(doubling_index == IFLAG_220_MOHO) ispec_aniso = ispec_aniso + 1

      enddo
    enddo
  enddo

! end of loop on all the subregions of the current region the mesh
  enddo

  print *
  print *,'found ',ispec_aniso,' tranversely isotropic elements in the mantle'

! create include file for the solver
  call save_header_file(NSPEC_AB,NSPEC_AC,NSPEC_BC, &
        nglob_AB,nglob_AC,nglob_BC,NEX_XI,NEX_ETA,ispec_aniso,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE, &
! BS
!       ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION)
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION,ATTENUATION_3D)
! BS END

  print *
  print *,'edit file OUTPUT_FILES/values_from_mesher.h to see some statistics about the mesh'
  print *
  print *,'on NEC SX-5 and Earth Simulator, make sure "loopcnt=" parameter'
! use fused loops on the ES
  print *,'in Makefile is greater than max vector length = ',nglob_BC(IREGION_CRUST_MANTLE)*NDIM

  print *
  print *,'done'
  print *

  end program create_header_file

