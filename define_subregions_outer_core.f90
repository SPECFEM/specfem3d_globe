!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine define_subregions_outer_core(isubregion,ichunk,iaddx,iaddy,iaddz, &
        ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
        doubling_index,npx,npy,NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_DOUBLING_OUTER_CORE,NER_ICB_BOTTOMDBL)

! define shape of elements in current subregion of the mesh
! the position of the top of the doubling region in the outer core
! is controlled by NER_DOUBLING_OUTER_CORE

  implicit none

  include "constants.h"

  integer ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
  integer iax,iay,iar
  integer isubregion,doubling_index
  integer npx,npy

  integer NER_TOP_CENTRAL_CUBE_ICB,NER_CENTRAL_CUBE_CMB,NER_DOUBLING_OUTER_CORE,NER_ICB_BOTTOMDBL

! topology of the elements
  integer iaddx(NGNOD)
  integer iaddy(NGNOD)
  integer iaddz(NGNOD)

  integer ichunk

! **************

  if(isubregion == 1) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*NER_DOUBLING_OUTER_CORE
    ir2=2*NER_CENTRAL_CUBE_CMB-24
    dir=8

    iax=4
    iay=4
    iar=4

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 2) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

    ir1=2*NER_TOP_CENTRAL_CUBE_ICB + 8
    ir2=ir1
    dir=8

    iax=4*2
    iay=4*2
    iar=4

    doubling_index = IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 3) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

    ir1=2*NER_TOP_CENTRAL_CUBE_ICB
    ir2=ir1
    dir=8

    iax=4*2
    iay=4*2
    iar=4

    doubling_index = IFLAG_BOTTOM_OUTER_CORE

  else if(isubregion == 4) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*NER_CENTRAL_CUBE_CMB-16
    ir2=2*NER_CENTRAL_CUBE_CMB-16
    dir=8

    iax=4
    iay=4
    iar=4

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 5) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8
    diy=8

    ix1=0
    ix2=npx-8
    dix=8

    ir1=2*NER_CENTRAL_CUBE_CMB-8
    ir2=2*NER_CENTRAL_CUBE_CMB-8
    dir=8

    iax=4
    iay=4
    iar=4

    doubling_index=IFLAG_TOP_OUTER_CORE

  else if(isubregion == 6) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

! subtract the first two layers, which are coupled to the inner core
    ir1=2*NER_TOP_CENTRAL_CUBE_ICB + 16
    ir2=2*NER_TOP_CENTRAL_CUBE_ICB + 16 + 8*(NER_ICB_BOTTOMDBL/4 -2 -1)
    dir=8

    iax=4*2
    iay=4*2
    iar=4

    doubling_index=IFLAG_OUTER_CORE_NORMAL

! added doubling at the bottom of the outer core

  else if(isubregion == 7) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 1 of the mesh doubling in the outer core

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
    dir=8*2

    iax=4*2
    iay=4*2
    iar=4*2
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
    dir=8*2

    iax=2*2
    iay=4*2
    iar=4*2
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    iy1=0
    iy2=npy-4*2
    diy=4*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
    dir=8*2

    iax=2*2
    iay=2*2
    iar=4*2
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 8) then

    call unusual_hex_nodes1(iaddx,iaddy,iaddz)

! generating stage 2 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-16*2
    dix=16*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 9) then

    call unusual_hex_nodes1p(iaddx,iaddy,iaddz)

! generating stage 3 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=8*2
    ix2=npx-8*2
    dix=16*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 10) then

    call unusual_hex_nodes2(iaddx,iaddy,iaddz)

! generating stage 4 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-16*2
    dix=16*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-24*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 11) then

    call unusual_hex_nodes2p(iaddx,iaddy,iaddz)

! generating stage 5 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=12*2
    ix2=npx-4*2
    dix=16*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-12*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-20*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-20*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 12) then

    call unusual_hex_nodes3(iaddx,iaddy,iaddz)

! generating stage 6 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=4*2
    ix2=npx-12*2
    dix=16*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-12*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-20*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-20*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 13) then

    call unusual_hex_nodes3(iaddx,iaddy,iaddz)

! generating stage 7 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=8*2
    ix2=npx-8*2
    dix=16*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-12*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-20*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-20*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 14) then

    call unusual_hex_nodes4(iaddx,iaddy,iaddz)

! generating stage 8 of the mesh doubling in the outer core

    iy1=8*2
    iy2=npy-8*2
    diy=16*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 15) then

    call unusual_hex_nodes4p(iaddx,iaddy,iaddz)

! generating stage 9 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-16*2
    diy=16*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 16) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 10 of the mesh doubling in the outer core

    iy1=8*2
    iy2=npy-8*2
    diy=16*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-4*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-4*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-12*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 17) then

    call usual_hex_nodes(iaddx,iaddy,iaddz)

! generating stage 11 of the mesh doubling in the outer core

    iy1=4*2
    iy2=npy-12*2
    diy=16*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-4*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-4*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-12*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 18) then

    call unusual_hex_nodes6(iaddx,iaddy,iaddz)

! generating stage 12 of the mesh doubling in the outer core

    iy1=12*2
    iy2=npy-4*2
    diy=16*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-4*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-4*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-12*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else if(isubregion == 19) then

    call unusual_hex_nodes6p(iaddx,iaddy,iaddz)

! generating stage 13 of the mesh doubling in the outer core

    iy1=0
    iy2=npy-16*2
    diy=16*2

    ix1=0
    ix2=npx-4*2
    dix=4*2

    dir=8*2

    iax=2*2
    iay=2*2
    iar=2*2

  if(ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
  else if(ichunk == CHUNK_AC .or. ichunk == CHUNK_AC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-8*2
    ir2=ir1
  else if(ichunk == CHUNK_BC .or. ichunk == CHUNK_BC_ANTIPODE) then
    ir1=2*NER_DOUBLING_OUTER_CORE-16*2
    ir2=ir1
  endif

    doubling_index=IFLAG_OUTER_CORE_NORMAL

  else

    stop 'incorrect subregion code'

  endif

  end subroutine define_subregions_outer_core

