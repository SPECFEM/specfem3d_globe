!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine reduce(theta,phi)

! bring theta between 0 and PI, and phi between 0 and 2*PI

  implicit none

  include "constants.h"

  double precision theta,phi

  integer i
  double precision th,ph

  th=theta
  ph=phi
  i=abs(int(ph/TWO_PI))
  if(ph<ZERO) then
    ph=ph+(i+1)*TWO_PI
  else
    if(ph>TWO_PI) ph=ph-i*TWO_PI
  endif
  phi=ph
  if(th<ZERO .or. th>PI) then
    i=int(th/PI)
    if(th>ZERO) then
      if(mod(i,2) /= 0) then
        th=(i+1)*PI-th
        if(ph<PI) then
          ph=ph+PI
        else
          ph=ph-PI
        endif
      else
        th=th-i*PI
      endif
    else
      if(mod(i,2) == 0) then
        th=-th+i*PI
        if(ph<PI) then
          ph=ph+PI
        else
          ph=ph-PI
        endif
      else
        th=th-i*PI
      endif
    endif
    theta=th
    phi=ph
  endif

  if(theta<ZERO .or. theta>PI) stop 'theta out of range in reduce'

  if(phi<ZERO .or. phi>TWO_PI) stop 'phi out of range in reduce'

  end subroutine reduce

