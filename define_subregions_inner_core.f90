!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
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

  subroutine define_subregions_inner_core(isubregion,iaddx,iaddy,iaddz, &
        ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir,iax,iay,iar, &
        doubling_index,npx,npy,NER_TOP_CENTRAL_CUBE_ICB)

! define shape of elements in current subregion of the mesh

  implicit none

  include "constants.h"

  integer ix1,ix2,dix,iy1,iy2,diy,ir1,ir2,dir
  integer iax,iay,iar
  integer isubregion,doubling_index
  integer npx,npy,NER_TOP_CENTRAL_CUBE_ICB

! topology of the elements
  integer iaddx(NGNOD)
  integer iaddy(NGNOD)
  integer iaddz(NGNOD)

! **************

! elements are always regular in inner core
  call usual_hex_nodes(iaddx,iaddy,iaddz)

! top of the inner core (one layer of elements)
  if(isubregion == 1) then

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

    ir1=2*NER_TOP_CENTRAL_CUBE_ICB - 8
    ir2=ir1
    dir=8

    iax=4*2
    iay=4*2
    iar=4

    doubling_index = IFLAG_TOP_INNER_CORE

! second coupling layer (one layer of elements)
  else if(isubregion == 2) then

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

    ir1=2*NER_TOP_CENTRAL_CUBE_ICB - 16
    ir2=ir1
    dir=8

    iax=4*2
    iay=4*2
    iar=4

    doubling_index = IFLAG_TOP_INNER_CORE_LEV2

! layer in contact with central cube
  else if(isubregion == 3) then

    iy1=0
    iy2=npy-8*2
    diy=8*2

    ix1=0
    ix2=npx-8*2
    dix=8*2

    ir1=0
    ir2=2*NER_TOP_CENTRAL_CUBE_ICB - 24
    dir=8

    iax=4*2
    iay=4*2
    iar=4

    doubling_index = IFLAG_INNER_CORE_NORMAL

  else

    stop 'incorrect subregion code'

  endif

  end subroutine define_subregions_inner_core

