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

  subroutine compute_parameters(NER_CRUST,NER_220_MOHO,NER_400_220, &
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

  implicit none

  include "constants.h"

! parameters read from parameter file
  integer NER_CRUST,NER_220_MOHO,NER_400_220, &
          NER_600_400,NER_670_600,NER_771_670,NER_TOPDDOUBLEPRIME_771, &
          NER_CMB_TOPDDOUBLEPRIME,NER_ICB_CMB,NER_TOP_CENTRAL_CUBE_ICB, &
          NER_ICB_BOTTOMDBL,NER_TOPDBL_CMB
  integer NEX_XI,NEX_ETA,NPROC_XI,NPROC_ETA

! parameters to be computed based upon parameters above read from file
  integer NPROC,NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
      NER,NER_CMB_670,NER_670_400,NER_CENTRAL_CUBE_CMB

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
      NSPEC2D_A_XI,NSPEC2D_B_XI,NSPEC2D_C_XI, &
      NSPEC2D_A_ETA,NSPEC2D_B_ETA,NSPEC2D_C_ETA, &
      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
      NSPEC1D_RADIAL,NPOIN1D_RADIAL, &
      NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX, &
      NGLOB_AB,NGLOB_AC,NGLOB_BC

  integer NEX_DOUBLING_MOHO_XI,NEX_DOUBLING_MOHO_ETA
  integer NEX_DOUBLING_MOHO_PER_PROC_XI,NEX_DOUBLING_MOHO_PER_PROC_ETA
  integer NSPEC2D_DOUBLING_A_XI,NSPEC2D_DOUBLING_A_ETA
  integer NSPEC2D_DOUBLING_B_XI,NSPEC2D_DOUBLING_B_ETA
  integer NSPEC2D_DOUBLING_C_XI,NSPEC2D_DOUBLING_C_ETA
  integer NSPEC2D_DOUBLING_A_XI_OC,NSPEC2D_DOUBLING_A_ETA_OC
  integer NSPEC2D_DOUBLING_B_XI_OC,NSPEC2D_DOUBLING_B_ETA_OC
  integer NSPEC2D_DOUBLING_C_XI_OC,NSPEC2D_DOUBLING_C_ETA_OC
  integer NSPEC_DOUBLING_AB,NSPEC_DOUBLING_AC,NSPEC_DOUBLING_BC
  integer NSPEC_DOUBLING_AB_OC,NSPEC_DOUBLING_AC_OC,NSPEC_DOUBLING_BC_OC
  integer NSPEC_DOUBLING_BC_ANISO
  integer NUM_DOUBLING_BRICKS,NUM_DOUBLING_BRICKS_ANISO,NUM_DOUBLING_BRICKS_OC
  integer NUM2D_DOUBLING_BRICKS_XI,NUM2D_DOUBLING_BRICKS_ETA
  integer NUM2D_DOUBLING_BRICKS_XI_OC,NUM2D_DOUBLING_BRICKS_ETA_OC
  integer nglob_no_doubling_volume,nglob_no_doubling_surface
  integer nblocks_xi,nblocks_eta
  integer nglob_surface_typeA,nglob_surface_typeB,nglob_surface_typeC

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_NO_DOUBLING,NSPEC2D_NO_DOUBLING_XI,NSPEC2D_NO_DOUBLING_ETA

! compute parameters for the Earth

! number of elements in radial direction
  NER_670_400 = NER_670_600 + NER_600_400
  NER_CMB_670 = NER_CMB_TOPDDOUBLEPRIME + NER_TOPDDOUBLEPRIME_771 + NER_771_670
  NER_CENTRAL_CUBE_CMB = NER_TOP_CENTRAL_CUBE_ICB + NER_ICB_CMB

! total number of spectral elements along radius
  NER = NER_CENTRAL_CUBE_CMB + NER_CMB_670 + NER_670_400 + NER_400_220 &
                    + NER_220_MOHO + NER_CRUST

! number of elements horizontally in each slice (i.e. per processor)
! these two values MUST be equal in all cases
  NEX_PER_PROC_XI = NEX_XI/NPROC_XI
  NEX_PER_PROC_ETA = NEX_ETA/NPROC_ETA

! total number of processors in each of the six chunks
  NPROC = NPROC_XI * NPROC_ETA

! total number of processors in the full Earth composed of the six chunks
  NPROCTOT = NCHUNKS * NPROC

! number of spectral elements at the bottom of the doubling below the moho
  NEX_DOUBLING_MOHO_XI=NEX_XI/2
  NEX_DOUBLING_MOHO_ETA=NEX_ETA/2
  NEX_DOUBLING_MOHO_PER_PROC_XI=NEX_PER_PROC_XI/2
  NEX_DOUBLING_MOHO_PER_PROC_ETA=NEX_PER_PROC_ETA/2

! exact number of spectral elements for a chunk without doubling layers

! in the crust and mantle
  NSPEC_NO_DOUBLING(IREGION_CRUST_MANTLE) = NEX_XI*NEX_ETA*NER_CRUST + ( &
    +NEX_DOUBLING_MOHO_XI*NEX_DOUBLING_MOHO_ETA*(NER_220_MOHO/2-3) &
    +NEX_DOUBLING_MOHO_XI*NEX_DOUBLING_MOHO_ETA*(NER_400_220+NER_670_400)/2 &
    +(NEX_XI/4)*(NEX_ETA/4)*(NER_CMB_670/4-3))

! in the outer core
  NSPEC_NO_DOUBLING(IREGION_OUTER_CORE) = &
    (NEX_XI/4)*(NEX_ETA/4)*(NER_TOPDBL_CMB/4) + (NEX_XI/8)*(NEX_ETA/8)*(NER_ICB_BOTTOMDBL/4)

! in the inner core
  NSPEC_NO_DOUBLING(IREGION_INNER_CORE) = (NEX_XI/(4*2))*(NEX_ETA/(4*2))*NER_TOP_CENTRAL_CUBE_ICB/4

! exact number of spectral elements in the doubling regions

! number of elementary bricks in the two regions with doubling
  NUM_DOUBLING_BRICKS = ((NEX_XI/4)*(NEX_ETA/4) &
        +NEX_DOUBLING_MOHO_XI*NEX_DOUBLING_MOHO_ETA)/4

! for type AB, each doubling brick contains 40 elements on 3 levels
  NSPEC_DOUBLING_AB=40*NUM_DOUBLING_BRICKS

! for type AC, each doubling brick contains 44 elements on 3 levels
  NSPEC_DOUBLING_AC=44*NUM_DOUBLING_BRICKS

! for type BC, each doubling brick contains 52 elements on 3 levels
  NSPEC_DOUBLING_BC=52*NUM_DOUBLING_BRICKS

! special case of elements with anisotropic properties, max found in BC chunks
  NUM_DOUBLING_BRICKS_ANISO = NEX_DOUBLING_MOHO_XI*NEX_DOUBLING_MOHO_ETA/4
  NSPEC_DOUBLING_BC_ANISO=52*NUM_DOUBLING_BRICKS_ANISO

! doubling in the outer core

! number of elementary bricks in the two regions with doubling
  NUM_DOUBLING_BRICKS_OC = ((NEX_XI/8)*(NEX_ETA/8))/4

! for type AB, each doubling brick contains 40 elements on 3 levels
  NSPEC_DOUBLING_AB_OC=40*NUM_DOUBLING_BRICKS_OC

! for type AC, each doubling brick contains 44 elements on 3 levels
  NSPEC_DOUBLING_AC_OC=44*NUM_DOUBLING_BRICKS_OC

! for type BC, each doubling brick contains 52 elements on 3 levels
  NSPEC_DOUBLING_BC_OC=52*NUM_DOUBLING_BRICKS_OC


! %%%%%%%%%%%%%% surface elements %%%%%%%%%%%%%%%%%%%

! exact number of surface elements for a chunk without doubling layers

! in the crust and mantle
  NSPEC2D_NO_DOUBLING_XI(IREGION_CRUST_MANTLE) = NEX_PER_PROC_XI*NER_CRUST + &
      NEX_DOUBLING_MOHO_PER_PROC_XI*(NER_220_MOHO/2-3) &
     +NEX_DOUBLING_MOHO_PER_PROC_XI*(NER_400_220+NER_670_400)/2 &
     +(NEX_PER_PROC_XI/4)*(NER_CMB_670/4-3)
  NSPEC2D_NO_DOUBLING_ETA(IREGION_CRUST_MANTLE) = NEX_PER_PROC_ETA*NER_CRUST + &
      NEX_DOUBLING_MOHO_PER_PROC_ETA*(NER_220_MOHO/2-3) &
     +NEX_DOUBLING_MOHO_PER_PROC_ETA*(NER_400_220+NER_670_400)/2 &
     +(NEX_PER_PROC_ETA/4)*(NER_CMB_670/4-3)

! in the outer core
  NSPEC2D_NO_DOUBLING_XI(IREGION_OUTER_CORE) = &
    (NEX_PER_PROC_XI/4)*(NER_TOPDBL_CMB/4) + (NEX_PER_PROC_XI/8)*(NER_ICB_BOTTOMDBL/4)
  NSPEC2D_NO_DOUBLING_ETA(IREGION_OUTER_CORE) = &
    (NEX_PER_PROC_ETA/4)*(NER_TOPDBL_CMB/4) + (NEX_PER_PROC_ETA/8)*(NER_ICB_BOTTOMDBL/4)

! in the inner core
  NSPEC2D_NO_DOUBLING_XI(IREGION_INNER_CORE) = (NEX_PER_PROC_XI/(4*2))*NER_TOP_CENTRAL_CUBE_ICB/4
  NSPEC2D_NO_DOUBLING_ETA(IREGION_INNER_CORE) = (NEX_PER_PROC_ETA/(4*2))*NER_TOP_CENTRAL_CUBE_ICB/4
  if(INCLUDE_CENTRAL_CUBE) then
    NSPEC2D_NO_DOUBLING_XI(IREGION_INNER_CORE) = NSPEC2D_NO_DOUBLING_XI(IREGION_INNER_CORE) + &
      (NEX_PER_PROC_XI / 8) * (NEX_XI / 8)
    NSPEC2D_NO_DOUBLING_ETA(IREGION_INNER_CORE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_INNER_CORE) + &
      (NEX_PER_PROC_ETA / 8) * (NEX_XI / 8)
  endif

! exact number of surface elements in the doubling regions

! number of elementary bricks in the two regions with doubling
  NUM2D_DOUBLING_BRICKS_XI = ((NEX_PER_PROC_XI/4) &
        +NEX_DOUBLING_MOHO_PER_PROC_XI)/2

  NUM2D_DOUBLING_BRICKS_ETA = ((NEX_PER_PROC_ETA/4) &
        +NEX_DOUBLING_MOHO_PER_PROC_ETA)/2

! for type A, each doubling brick contains 10 elements on 3 levels
  NSPEC2D_DOUBLING_A_XI=10*NUM2D_DOUBLING_BRICKS_XI
  NSPEC2D_DOUBLING_A_ETA=10*NUM2D_DOUBLING_BRICKS_ETA

! for type B, each doubling brick contains 12 elements on 3 levels
  NSPEC2D_DOUBLING_B_XI=12*NUM2D_DOUBLING_BRICKS_XI
  NSPEC2D_DOUBLING_B_ETA=12*NUM2D_DOUBLING_BRICKS_ETA

! for type C, each doubling brick contains 14 elements on 3 levels
  NSPEC2D_DOUBLING_C_XI=14*NUM2D_DOUBLING_BRICKS_XI
  NSPEC2D_DOUBLING_C_ETA=14*NUM2D_DOUBLING_BRICKS_ETA

! doubling in outer core

! number of elementary bricks in the two regions with doubling
  NUM2D_DOUBLING_BRICKS_XI_OC = ((NEX_PER_PROC_XI/8))/2

  NUM2D_DOUBLING_BRICKS_ETA_OC = ((NEX_PER_PROC_ETA/8))/2

! for type A, each doubling brick contains 10 elements on 3 levels
  NSPEC2D_DOUBLING_A_XI_OC=10*NUM2D_DOUBLING_BRICKS_XI_OC
  NSPEC2D_DOUBLING_A_ETA_OC=10*NUM2D_DOUBLING_BRICKS_ETA_OC

! for type B, each doubling brick contains 12 elements on 3 levels
  NSPEC2D_DOUBLING_B_XI_OC=12*NUM2D_DOUBLING_BRICKS_XI_OC
  NSPEC2D_DOUBLING_B_ETA_OC=12*NUM2D_DOUBLING_BRICKS_ETA_OC

! for type C, each doubling brick contains 14 elements on 3 levels
  NSPEC2D_DOUBLING_C_XI_OC=14*NUM2D_DOUBLING_BRICKS_XI_OC
  NSPEC2D_DOUBLING_C_ETA_OC=14*NUM2D_DOUBLING_BRICKS_ETA_OC

! exact number of surface elements for faces A, B and C along XI and ETA
  NSPEC2D_A_XI(:) = NSPEC2D_NO_DOUBLING_XI(:)
  NSPEC2D_B_XI(:) = NSPEC2D_A_XI(:)
  NSPEC2D_C_XI(:) = NSPEC2D_A_XI(:)
  NSPEC2D_A_ETA(:) = NSPEC2D_NO_DOUBLING_ETA(:)
  NSPEC2D_B_ETA(:) = NSPEC2D_A_ETA(:)
  NSPEC2D_C_ETA(:) = NSPEC2D_A_ETA(:)

! exact number of spectral elements in crust and mantle
  NSPEC_AB(IREGION_CRUST_MANTLE) = (NSPEC_NO_DOUBLING(IREGION_CRUST_MANTLE) + NSPEC_DOUBLING_AB) / NPROC
  NSPEC_AC(IREGION_CRUST_MANTLE) = (NSPEC_NO_DOUBLING(IREGION_CRUST_MANTLE) + NSPEC_DOUBLING_AC) / NPROC
  NSPEC_BC(IREGION_CRUST_MANTLE) = (NSPEC_NO_DOUBLING(IREGION_CRUST_MANTLE) + NSPEC_DOUBLING_BC) / NPROC

! exact number of spectral elements in outer core
  NSPEC_AB(IREGION_OUTER_CORE) = (NSPEC_NO_DOUBLING(IREGION_OUTER_CORE) + NSPEC_DOUBLING_AB_OC) / NPROC
  NSPEC_AC(IREGION_OUTER_CORE) = (NSPEC_NO_DOUBLING(IREGION_OUTER_CORE) + NSPEC_DOUBLING_AC_OC) / NPROC
  NSPEC_BC(IREGION_OUTER_CORE) = (NSPEC_NO_DOUBLING(IREGION_OUTER_CORE) + NSPEC_DOUBLING_BC_OC) / NPROC

! exact number of spectral elements in inner core
  NSPEC_AB(IREGION_INNER_CORE) = NSPEC_NO_DOUBLING(IREGION_INNER_CORE) / NPROC
  if(INCLUDE_CENTRAL_CUBE) NSPEC_AB(IREGION_INNER_CORE) = NSPEC_AB(IREGION_INNER_CORE) + &
         (NEX_PER_PROC_XI / 8) * (NEX_PER_PROC_ETA / 8) * (NEX_XI / 8)
  NSPEC_AC(IREGION_INNER_CORE) = NSPEC_AB(IREGION_INNER_CORE)
  NSPEC_BC(IREGION_INNER_CORE) = NSPEC_AB(IREGION_INNER_CORE)

! exact number of surface elements for faces A, B and C along XI and ETA for doubling region
  NSPEC2D_A_XI(IREGION_CRUST_MANTLE) = NSPEC2D_NO_DOUBLING_XI(IREGION_CRUST_MANTLE) + NSPEC2D_DOUBLING_A_XI
  NSPEC2D_B_XI(IREGION_CRUST_MANTLE) = NSPEC2D_NO_DOUBLING_XI(IREGION_CRUST_MANTLE) + NSPEC2D_DOUBLING_B_XI
  NSPEC2D_C_XI(IREGION_CRUST_MANTLE) = NSPEC2D_NO_DOUBLING_XI(IREGION_CRUST_MANTLE) + NSPEC2D_DOUBLING_C_XI
  NSPEC2D_A_ETA(IREGION_CRUST_MANTLE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_CRUST_MANTLE) + NSPEC2D_DOUBLING_A_ETA
  NSPEC2D_B_ETA(IREGION_CRUST_MANTLE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_CRUST_MANTLE) + NSPEC2D_DOUBLING_B_ETA
  NSPEC2D_C_ETA(IREGION_CRUST_MANTLE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_CRUST_MANTLE) + NSPEC2D_DOUBLING_C_ETA

! added doubling region in the outer core
  NSPEC2D_A_XI(IREGION_OUTER_CORE) = NSPEC2D_NO_DOUBLING_XI(IREGION_OUTER_CORE) + NSPEC2D_DOUBLING_A_XI_OC
  NSPEC2D_B_XI(IREGION_OUTER_CORE) = NSPEC2D_NO_DOUBLING_XI(IREGION_OUTER_CORE) + NSPEC2D_DOUBLING_B_XI_OC
  NSPEC2D_C_XI(IREGION_OUTER_CORE) = NSPEC2D_NO_DOUBLING_XI(IREGION_OUTER_CORE) + NSPEC2D_DOUBLING_C_XI_OC
  NSPEC2D_A_ETA(IREGION_OUTER_CORE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_OUTER_CORE) + NSPEC2D_DOUBLING_A_ETA_OC
  NSPEC2D_B_ETA(IREGION_OUTER_CORE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_OUTER_CORE) + NSPEC2D_DOUBLING_B_ETA_OC
  NSPEC2D_C_ETA(IREGION_OUTER_CORE) = NSPEC2D_NO_DOUBLING_ETA(IREGION_OUTER_CORE) + NSPEC2D_DOUBLING_C_ETA_OC

! exact number of surface elements on the bottom and top boundaries
! and theoretical number of spectral elements in radial direction

! in the crust and mantle
  NSPEC2D_TOP(IREGION_CRUST_MANTLE) = NEX_XI*NEX_ETA/NPROC
  NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE) = (NEX_XI/4)*(NEX_ETA/4)/NPROC
  NSPEC1D_RADIAL(IREGION_CRUST_MANTLE) = NER_CRUST + (NER_220_MOHO+NER_400_220+NER_670_400)/2 + NER_CMB_670/4

! in the outer core with mesh doubling
  NSPEC2D_TOP(IREGION_OUTER_CORE) = (NEX_XI/4)*(NEX_ETA/4)/NPROC
  NSPEC2D_BOTTOM(IREGION_OUTER_CORE) = (NEX_XI/4/2)*(NEX_ETA/4/2)/NPROC
  NSPEC1D_RADIAL(IREGION_OUTER_CORE) = NER_ICB_CMB/4 - 3

! in the top of the inner core
  NSPEC2D_TOP(IREGION_INNER_CORE) = (NEX_XI/(4*2))*(NEX_ETA/(4*2))/NPROC
  NSPEC2D_BOTTOM(IREGION_INNER_CORE) = NSPEC2D_TOP(IREGION_INNER_CORE)
  NSPEC1D_RADIAL(IREGION_INNER_CORE) = NER_TOP_CENTRAL_CUBE_ICB/4

! maximum number of surface elements on vertical boundaries of the slices
  NSPEC2DMAX_XMIN_XMAX(:) = NSPEC2D_C_ETA(:)
  NSPEC2DMAX_YMIN_YMAX(:) = NSPEC2D_C_XI(:)

! theoretical number of Gauss-Lobatto points in radial direction
  NPOIN1D_RADIAL(:)=NSPEC1D_RADIAL(:)*(NGLLZ-1)+1

! 2-D addressing and buffers for summation between slices
! we add one to number of points because of the flag after the last point
  NPOIN2DMAX_XMIN_XMAX(:) = NSPEC2DMAX_XMIN_XMAX(:)*NGLLY*NGLLZ + 1
  NPOIN2DMAX_YMIN_YMAX(:) = NSPEC2DMAX_YMIN_YMAX(:)*NGLLX*NGLLZ + 1

! exact number of global points in each region

! initialize array
  NGLOB_AB(:) = 0

! in the inner core
  if(INCLUDE_CENTRAL_CUBE) then
    NGLOB_AB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/8) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/8) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB/4 + NEX_XI / 8)*(NGLLZ-1)+1)
  else
    NGLOB_AB(IREGION_INNER_CORE) = ((NEX_PER_PROC_XI/8) &
      *(NGLLX-1)+1)*((NEX_PER_PROC_ETA/8) &
      *(NGLLY-1)+1)*((NER_TOP_CENTRAL_CUBE_ICB/4)*(NGLLZ-1)+1)
  endif

! blocks in the inner core have the same number of points
! but not blocks in the mantle or in the outer core because of mesh doubling
  NGLOB_AC(:) = NGLOB_AB(:)
  NGLOB_BC(:) = NGLOB_AB(:)

! ---
! case of the mantle and crust that have doubling regions
! values computed using Mathematica for the basic 3D doubling brick
! ---

! compute number of points in blocks with no doubling
! exclude the three surfaces in contact with the doubling regions
  nglob_no_doubling_volume = &
     (4*(NGLLX-1)+1)*(4*(NGLLX-1)+1)*((NER_220_MOHO/2-3 + (NER_400_220+NER_670_400)/2)*(NGLLX-1)-1) &
    + (2*(NGLLX-1)+1)*(2*(NGLLX-1)+1)*((NER_CMB_670/4-3)*(NGLLX-1)+0) &
    + (8*(NGLLX-1)+1)*(8*(NGLLX-1)+1)*(NER_CRUST*(NGLLX-1)+0)

! number of basic blocks in each slice
  nblocks_xi = NEX_PER_PROC_XI / 8
  nblocks_eta = NEX_PER_PROC_ETA / 8

! in the crust and mantle
  NGLOB_AB(IREGION_CRUST_MANTLE) = &
    nblocks_xi*nblocks_eta*(200*NGLLX**3 - 484*NGLLX**2 + 392*NGLLX - 106 + nglob_no_doubling_volume)
  NGLOB_AC(IREGION_CRUST_MANTLE) = &
    nblocks_xi*nblocks_eta*(220*NGLLX**3 - 538*NGLLX**2 + 440*NGLLX - 120 + nglob_no_doubling_volume)
  NGLOB_BC(IREGION_CRUST_MANTLE) = &
    nblocks_xi*nblocks_eta*(260*NGLLX**3 - 652*NGLLX**2 + 548*NGLLX - 154 + nglob_no_doubling_volume)

! same thing for 2D surfaces for the three types of faces
  nglob_no_doubling_surface = (8*(NGLLX-1)+1)*(NER_CRUST*(NGLLX-1)+0) &
    + (4*(NGLLX-1)+1)*((NER_220_MOHO/2-3 + (NER_400_220+NER_670_400)/2)*(NGLLX-1)-1) &
    + (2*(NGLLX-1)+1)*((NER_CMB_670/4-3)*(NGLLX-1)+0)

  nglob_surface_typeA = 30*NGLLX**2 - 45 * NGLLX + 17
  nglob_surface_typeB = 36*NGLLX**2 - 57 * NGLLX + 23
  nglob_surface_typeC = 42*NGLLX**2 - 69 * NGLLX + 29

! final number of points in volume obtained by removing planes counted twice
  NGLOB_AB(IREGION_CRUST_MANTLE) = nglob_AB(IREGION_CRUST_MANTLE) &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeA + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeB + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL(IREGION_CRUST_MANTLE)

  NGLOB_AC(IREGION_CRUST_MANTLE) = nglob_AC(IREGION_CRUST_MANTLE) &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeA + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeC + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL(IREGION_CRUST_MANTLE)

  NGLOB_BC(IREGION_CRUST_MANTLE) = nglob_BC(IREGION_CRUST_MANTLE) &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeB + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeC + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL(IREGION_CRUST_MANTLE)

! ---
! case of the outer core that also has a mesh doubling region
! values computed using Mathematica for the basic 3D doubling brick
! ---

! compute number of points in blocks with no doubling
! exclude the three surfaces in contact with the doubling regions
  nglob_no_doubling_volume = (4*(NGLLX-1)+1)*(4*(NGLLX-1)+1)*((NER_TOPDBL_CMB/4)*(NGLLX-1)+0) + &
    (2*(NGLLX-1)+1)*(2*(NGLLX-1)+1)*((NER_ICB_BOTTOMDBL/4)*(NGLLX-1)+0)

! number of basic blocks in each slice
  nblocks_xi = NEX_PER_PROC_XI / 16
  nblocks_eta = NEX_PER_PROC_ETA / 16

  NGLOB_AB(IREGION_OUTER_CORE) = nblocks_xi*nblocks_eta*(40*NGLLX**3 - 88*NGLLX**2 + 65*NGLLX - 16 + nglob_no_doubling_volume)
  NGLOB_AC(IREGION_OUTER_CORE) = nblocks_xi*nblocks_eta*(44*NGLLX**3 - 98*NGLLX**2 + 73*NGLLX - 18 + nglob_no_doubling_volume)
  NGLOB_BC(IREGION_OUTER_CORE) = nblocks_xi*nblocks_eta*(52*NGLLX**3 - 120*NGLLX**2 + 93*NGLLX - 24 + nglob_no_doubling_volume)

! same thing for 2D surfaces for the three types of faces
  nglob_no_doubling_surface = (4*(NGLLX-1)+1)*((NER_TOPDBL_CMB/4)*(NGLLX-1)+0) + &
    (2*(NGLLX-1)+1)*((NER_ICB_BOTTOMDBL/4)*(NGLLX-1)+0)

  nglob_surface_typeA = 10*(NGLLX-2)**2 + 26*(NGLLX-2) + 17
  nglob_surface_typeB = 12*(NGLLX-2)**2 + 30*(NGLLX-2) + 19
  nglob_surface_typeC = 14*(NGLLX-2)**2 + 34*(NGLLX-2) + 21

! final number of points in volume obtained by removing planes counted twice
  NGLOB_AB(IREGION_OUTER_CORE) = nglob_AB(IREGION_OUTER_CORE) &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeA + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeB + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL(IREGION_OUTER_CORE)

  NGLOB_AC(IREGION_OUTER_CORE) = nglob_AC(IREGION_OUTER_CORE) &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeA + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeC + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL(IREGION_OUTER_CORE)

  NGLOB_BC(IREGION_OUTER_CORE) = nglob_BC(IREGION_OUTER_CORE) &
     - (nblocks_xi-1)*nblocks_eta*(nglob_surface_typeB + nglob_no_doubling_surface) &
     - (nblocks_eta-1)*nblocks_xi*(nglob_surface_typeC + nglob_no_doubling_surface) &
     + (nblocks_eta-1)*(nblocks_xi-1)*NPOIN1D_RADIAL(IREGION_OUTER_CORE)

  end subroutine compute_parameters

