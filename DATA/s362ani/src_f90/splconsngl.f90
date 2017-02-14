
  subroutine splcon(xlat,xlon,nver,verlat,verlon,verrad,ncon,icon,con)

  implicit none

  integer icon(1)

  real(kind=4) verlat(1)
  real(kind=4) verlon(1)
  real(kind=4) verrad(1)
  real(kind=4) con(1)

  double precision dd
  double precision rn
  double precision dr
  double precision xrad
  double precision ver8
  double precision xla8

  integer :: ncon,iver,nver

  real(kind=4) :: xlat,xlon

  xrad=3.14159265358979/180.d0

  ncon=0

  do iver=1,nver
  if (xlat > verlat(iver)-2.*verrad(iver)) then
    if (xlat < verlat(iver)+2.*verrad(iver)) then
      ver8=xrad*(verlat(iver))
      xla8=xrad*(xlat)
      dd=sin(ver8)*sin(xla8)
      dd=dd+cos(ver8)*cos(xla8)* cos(xrad*(xlon-verlon(iver)))
      dd=acos(dd)/xrad
      if (dd > (verrad(iver))*2.d0) then
      else
        ncon=ncon+1
        icon(ncon)=iver
        rn=dd/(verrad(iver))
        dr=rn-1.d0
        if (rn <= 1.d0) then
          con(ncon)=(0.75d0*rn-1.5d0)*(rn**2)+1.d0
        else if (rn > 1.d0) then
          con(ncon)=((-0.25d0*dr+0.75d0)*dr-0.75d0)*dr+0.25d0
        else
          con(ncon)=0.
        endif
      endif
    endif
  endif
  enddo

  end subroutine splcon

