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

  subroutine define_subregions_crust_mantle(isubregion,ichunk,iaddx,iaddy,iaddz, &
        ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
        doubling_index,npx,npy,NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400, &
        NER_400_220,NER_220_MOHO,NER_CRUST)

! define shape of elements in current subregion of the mesh

! do not move the order of the subregions here since
! we have to make sure anisotropic elements are created first
! (region number 1 at the end of the routine)
! in order to be able to store the anisotropic material properties
! for a contiguous subset of elements to save memory in the solver

  implicit none

  include "constants.h"

  integer ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
  integer iax,iay,iar
  integer isubregion,ichunk,doubling_index
  integer npx,npy

  integer NER_CENTRAL_CUBE_CMB,NER_CMB_670,NER_670_400,NER_400_220,NER_220_MOHO,NER_CRUST

! topology of the elements
  integer iaddx(NGNOD)
  integer iaddy(NGNOD)
  integer iaddz(NGNOD)

! **************

! this last region is the crust (elements above the Moho)

  if(isubregion == 32) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-2
    diy=2

    ix1=0
    ix2=npx-2
    dix=2

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)
    ir2=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO+NER_CRUST)-2
    dir=2

    iax=1
    iay=1
    iar=1

    doubling_index=IFLAG_CRUST

! all other regions below are the mantle (elements below the Moho)

  else if(isubregion == 31) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*NER_CENTRAL_CUBE_CMB  + 16
    ir2=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-32

    dir=8

    iax=4
    iay=4
    iar=4

    doubling_index=IFLAG_MANTLE_NORMAL

  else if(isubregion == 30) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*NER_CENTRAL_CUBE_CMB
    ir2=2*NER_CENTRAL_CUBE_CMB
    dir=8

    iax=4
    iay=4
    iar=4

    doubling_index=IFLAG_BOTTOM_MANTLE

  else if(isubregion == 29) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*NER_CENTRAL_CUBE_CMB + 8
    ir2=2*NER_CENTRAL_CUBE_CMB + 8
    dir=8

    iax=4
    iay=4
    iar=4

    doubling_index=IFLAG_MANTLE_NORMAL

  else if(isubregion == 28) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 1 of the mesh doubling below 670

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
    dir=8

    iax=4
    iay=4
    iar=4
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-4
    dix=4

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
    dir=8

    iax=2
    iay=4
    iar=4
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-4
    dix=4

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
    dir=8

    iax=2
    iay=2
    iar=4
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 27) then

    call unusual_hex_nodes1(iaddx,iaddy,iaddz)

! generating stage 2 of the mesh doubling below 670

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-16
    dix=16

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 26) then

    call unusual_hex_nodes1p(iaddx,iaddy,iaddz)

! generating stage 3 of the mesh doubling below 670

    iy1=0
    iy2=npy-8
    diy=8

    ix1=8
    ix2=npx-8
    dix=16

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 25) then

    call unusual_hex_nodes2(iaddx,iaddy,iaddz)

! generating stage 4 of the mesh doubling below 670

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-16
    dix=16

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-24
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 24) then

    call unusual_hex_nodes2p(iaddx,iaddy,iaddz)

! generating stage 5 of the mesh doubling below 670

    iy1=0
    iy2=npy-8
    diy=8

    ix1=12
    ix2=npx-4
    dix=16

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-12
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-20
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-20
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 23) then

    call unusual_hex_nodes3(iaddx,iaddy,iaddz)

! generating stage 6 of the mesh doubling below 670

    iy1=0
    iy2=npy-8
    diy=8

    ix1=4
    ix2=npx-12
    dix=16

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-12
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-20
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-20
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 22) then

    call unusual_hex_nodes3(iaddx,iaddy,iaddz)

! generating stage 7 of the mesh doubling below 670

    iy1=0
    iy2=npy-8
    diy=8

    ix1=8
    ix2=npx-8
    dix=16

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-12
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-20
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-20
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 21) then

    call unusual_hex_nodes4(iaddx,iaddy,iaddz)

! generating stage 8 of the mesh doubling below 670

    iy1=8
    iy2=npy-8
    diy=16

    ix1=0
    ix2=npx-4
    dix=4

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 20) then

    call unusual_hex_nodes4p(iaddx,iaddy,iaddz)

! generating stage 9 of the mesh doubling below 670

    iy1=0
    iy2=npy-16
    diy=16

    ix1=0
    ix2=npx-4
    dix=4

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 19) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 10 of the mesh doubling below 670

    iy1=8
    iy2=npy-8
    diy=16

    ix1=0
    ix2=npx-4
    dix=4

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-4
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-4
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-12
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 18) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 11 of the mesh doubling below 670

    iy1=4
    iy2=npy-12
    diy=16

    ix1=0
    ix2=npx-4
    dix=4

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-4
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-4
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-12
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 17) then

    call unusual_hex_nodes6(iaddx,iaddy,iaddz)

! generating stage 12 of the mesh doubling below 670

    iy1=12
    iy2=npy-4
    diy=16

    ix1=0
    ix2=npx-4
    dix=4

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-4
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-4
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-12
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 16) then

    call unusual_hex_nodes6p(iaddx,iaddy,iaddz)

! generating stage 13 of the mesh doubling below 670

    iy1=0
    iy2=npy-16
    diy=16

    ix1=0
    ix2=npx-4
    dix=4

    dir=8

    iax=2
    iay=2
    iar=2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-8
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)-16
    ir2=ir1
  endif

    doubling_index=IFLAG_DOUBLING_670

  else if(isubregion == 15) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-4
    dix=4

! honor discontinuity at 220 km for accuracy and for anisotropy
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670)
    ir2=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220)-4
    dir=4

    iax=2
    iay=2
    iar=2

    doubling_index=IFLAG_670_220

  else if(isubregion == 14) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-4
    dix=4

! honor discontinuity at 220 km for accuracy and for anisotropy
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220)
    ir2=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-16
    dir=4

    iax=2
    iay=2
    iar=2

    doubling_index=IFLAG_220_MOHO


  else if(isubregion == 13) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 1 of the mesh doubling below the Moho

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-4
    dix=4

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
    dir=4

    iax=2
    iay=2
    iar=2
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-2
    dix=2

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
    dir=2

    iax=1
    iay=2
    iar=2
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    iy1=0
    iy2=npy-2
    diy=2

    ix1=0
    ix2=npx-2
    dix=2

    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
    dir=2

    iax=1
    iay=1
    iar=2
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 12) then

    call unusual_hex_nodes1(iaddx,iaddy,iaddz)

! generating stage 2 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-8
    dix=8

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 11) then

    call unusual_hex_nodes1p(iaddx,iaddy,iaddz)

! generating stage 3 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-4
    diy=4

    ix1=4
    ix2=npx-4
    dix=8

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 10) then

    call unusual_hex_nodes2(iaddx,iaddy,iaddz)

! generating stage 4 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-4
    diy=4

    ix1=0
    ix2=npx-8
    dix=8

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-12
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 9) then

    call unusual_hex_nodes2p(iaddx,iaddy,iaddz)

! generating stage 5 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-4
    diy=4

    ix1=6
    ix2=npx-2
    dix=8

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-6
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-10
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-10
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 8) then

    call unusual_hex_nodes3(iaddx,iaddy,iaddz)

! generating stage 6 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-4
    diy=4

    ix1=2
    ix2=npx-6
    dix=8

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-6
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-10
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-10
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 7) then

    call unusual_hex_nodes3(iaddx,iaddy,iaddz)

! generating stage 7 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-4
    diy=4

    ix1=4
    ix2=npx-4
    dix=8

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-6
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-10
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-10
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 6) then

    call unusual_hex_nodes4(iaddx,iaddy,iaddz)

! generating stage 8 of the mesh doubling below the Moho

    iy1=4
    iy2=npy-4
    diy=8

    ix1=0
    ix2=npx-2
    dix=2

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 5) then

    call unusual_hex_nodes4p(iaddx,iaddy,iaddz)

! generating stage 9 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-2
    dix=2

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 4) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 10 of the mesh doubling below the Moho

    iy1=4
    iy2=npy-4
    diy=8

    ix1=0
    ix2=npx-2
    dix=2

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-6
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 3) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 11 of the mesh doubling below the Moho

    iy1=2
    iy2=npy-6
    diy=8

    ix1=0
    ix2=npx-2
    dix=2

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-6
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 2) then

    call unusual_hex_nodes6(iaddx,iaddy,iaddz)

! generating stage 12 of the mesh doubling below the Moho

    iy1=6
    iy2=npy-2
    diy=8

    ix1=0
    ix2=npx-2
    dix=2

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-6
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else if(isubregion == 1) then

    call unusual_hex_nodes6p(iaddx,iaddy,iaddz)

! generating stage 13 of the mesh doubling below the Moho

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-2
    dix=2

    dir=4

    iax=1
    iay=1
    iar=1

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-4
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*(NER_CENTRAL_CUBE_CMB+NER_CMB_670+NER_670_400+NER_400_220+NER_220_MOHO)-8
    ir2=ir1
  endif

    doubling_index=IFLAG_220_MOHO

  else

    stop 'incorrect subregion code'

  endif

  end subroutine define_subregions_crust_mantle

