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

  subroutine get_cmt(yr,jda,ho,mi,sec,t_cmt,hdur,elat,elon,depth,moment_tensor,DT,isource)

  implicit none

  include "constants.h"

  integer yr,jda,ho,mi,isource
  double precision sec,t_cmt,hdur,elat,elon,depth
  double precision moment_tensor(6)
  double precision DT

  integer mo,da,julian_day,iread
  double precision scaleM
  character(len=5) datasource
  character(len=150) string

!
!---- read hypocenter info
!
  open(unit=1,file='DATA/CMTSOLUTION',status='old')

! read source number isource
  do iread=1,isource

! read header with event information
  read(1,"(a4,i5,i3,i3,i3,i3,f6.2)") datasource,yr,mo,da,ho,mi,sec
  jda=julian_day(yr,mo,da)

! ignore line with event name
  read(1,"(a)") string

! read time shift
  read(1,"(a)") string
  read(string(12:len_trim(string)),*) t_cmt

! read half duration
  read(1,"(a)") string
  read(string(15:len_trim(string)),*) hdur

! read latitude
  read(1,"(a)") string
  read(string(10:len_trim(string)),*) elat

! read longitude
  read(1,"(a)") string
  read(string(11:len_trim(string)),*) elon

! read depth
  read(1,"(a)") string
  read(string(7:len_trim(string)),*) depth

! read Mrr
  read(1,"(a)") string
  read(string(5:len_trim(string)),*) moment_tensor(1)

! read Mtt
  read(1,"(a)") string
  read(string(5:len_trim(string)),*) moment_tensor(2)

! read Mpp
  read(1,"(a)") string
  read(string(5:len_trim(string)),*) moment_tensor(3)

! read Mrt
  read(1,"(a)") string
  read(string(5:len_trim(string)),*) moment_tensor(4)

! read Mrp
  read(1,"(a)") string
  read(string(5:len_trim(string)),*) moment_tensor(5)

! read Mtp
  read(1,"(a)") string
  read(string(5:len_trim(string)),*) moment_tensor(6)

  enddo

  close(1)

! null half-duration indicates a Heaviside
! replace with very short error function
  if(hdur < 5. * DT) hdur = 5. * DT

!
! scale the moment-tensor (dimensions dyn-cm)
!
  scaleM = 1.d7 * RHOAV * (R_EARTH**5) * PI*GRAV*RHOAV
  moment_tensor(:) = moment_tensor(:) / scaleM

  end subroutine get_cmt

! ------------------------------------------------------------------

  integer function julian_day(yr,mo,da)

  implicit none

  integer yr,mo,da

  integer mon(12)
  integer lpyr
  data mon /0,31,59,90,120,151,181,212,243,273,304,334/

  julian_day = da + mon(mo)
  if(mo>2) julian_day = julian_day + lpyr(yr)

  end function julian_day

! ------------------------------------------------------------------

  integer function lpyr(yr)

  implicit none

  integer yr
!
!---- returns 1 if leap year
!
  lpyr=0
  if(mod(yr,400) == 0) then
    lpyr=1
  else if(mod(yr,4) == 0) then
    lpyr=1
    if(mod(yr,100) == 0) lpyr=0
  endif

  end function lpyr

