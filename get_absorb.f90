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

  subroutine get_absorb(myrank,prname,iboun,nspec, &
        nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! Stacey, define flags for absorbing boundaries

  implicit none

  include "constants.h"

  integer nspec,myrank

  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM

  integer nimin(2,NSPEC2DMAX_YMIN_YMAX),nimax(2,NSPEC2DMAX_YMIN_YMAX)
  integer njmin(2,NSPEC2DMAX_XMIN_XMAX),njmax(2,NSPEC2DMAX_XMIN_XMAX)
  integer nkmin_xi(2,NSPEC2DMAX_XMIN_XMAX),nkmin_eta(2,NSPEC2DMAX_YMIN_YMAX)

  logical iboun(6,nspec)

! global element numbering
  integer ispecg

! counters to keep track of the number of elements on each of the
! five absorbing boundaries
  integer ispecb1,ispecb2,ispecb3,ispecb4,ispecb5

! processor identification
  character(len=150) prname

  ispecb1=0
  ispecb2=0
  ispecb3=0
  ispecb4=0
  ispecb5=0

  do ispecg=1,nspec

! determine if the element falls on an absorbing boundary

  if(iboun(1,ispecg)) then

!   on boundary 1: xmin
    ispecb1=ispecb1+1

! this is useful even if it is constant because it can be zero inside the slices
    njmin(1,ispecb1)=1
    njmax(1,ispecb1)=NGLLY

!   check for ovelap with other boundaries
    nkmin_xi(1,ispecb1)=1
    if(iboun(5,ispecg)) nkmin_xi(1,ispecb1)=2
  endif

  if(iboun(2,ispecg)) then

!   on boundary 2: xmax
    ispecb2=ispecb2+1

! this is useful even if it is constant because it can be zero inside the slices
    njmin(2,ispecb2)=1
    njmax(2,ispecb2)=NGLLY

!   check for ovelap with other boundaries
    nkmin_xi(2,ispecb2)=1
    if(iboun(5,ispecg)) nkmin_xi(2,ispecb2)=2
  endif

  if(iboun(3,ispecg)) then

!   on boundary 3: ymin
    ispecb3=ispecb3+1

!   check for ovelap with other boundaries
    nimin(1,ispecb3)=1
    if(iboun(1,ispecg)) nimin(1,ispecb3)=2
    nimax(1,ispecb3)=NGLLX
    if(iboun(2,ispecg)) nimax(1,ispecb3)=NGLLX-1
    nkmin_eta(1,ispecb3)=1
    if(iboun(5,ispecg)) nkmin_eta(1,ispecb3)=2
  endif

  if(iboun(4,ispecg)) then

!   on boundary 4: ymax
    ispecb4=ispecb4+1

!   check for ovelap with other boundaries
    nimin(2,ispecb4)=1
    if(iboun(1,ispecg)) nimin(2,ispecb4)=2
    nimax(2,ispecb4)=NGLLX
    if(iboun(2,ispecg)) nimax(2,ispecb4)=NGLLX-1
    nkmin_eta(2,ispecb4)=1
    if(iboun(5,ispecg)) nkmin_eta(2,ispecb4)=2
  endif

! on boundary 5: bottom
  if(iboun(5,ispecg)) ispecb5=ispecb5+1

  enddo

! check theoretical value of elements at the bottom
  if(ispecb5 /= NSPEC2D_BOTTOM) &
    call exit_MPI(myrank,'ispecb5 should equal NSPEC2D_BOTTOM in absorbing boundary detection')

  end subroutine get_absorb

