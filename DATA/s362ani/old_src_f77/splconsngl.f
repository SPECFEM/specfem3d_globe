      subroutine splcon(xlat,xlon,nver,verlat,verlon,verrad,
     #                  ncon,icon,con)
      dimension verlat(1)
      dimension verlon(1)
      dimension verrad(1)
      dimension icon(1)
      dimension con(1)
c
      real*8 dd
      real*8 rn
      real*8 dr
      real*8 xrad
      real*8 ver8
      real*8 xla8
c
      xrad=3.14159265358979/180.d0
      ncon=0
      do iver=1,nver
        if(xlat.gt.verlat(iver)-2.*verrad(iver)) then
          if(xlat.lt.verlat(iver)+2.*verrad(iver)) then
            ver8=xrad*(verlat(iver))
            xla8=xrad*(xlat)
            dd=sin(ver8)*sin(xla8)
            dd=dd+cos(ver8)*cos(xla8)*
     #         cos(xrad*(xlon-verlon(iver)))
            dd=acos(dd)/xrad
            if(dd.gt.(verrad(iver))*2.d0) then
            else
              ncon=ncon+1
              icon(ncon)=iver
              rn=dd/(verrad(iver))
              dr=rn-1.d0
              if(rn.le.1.d0) then
                con(ncon)=(0.75d0*rn-1.5d0)*(rn**2)+1.d0
              else if(rn.gt.1.d0) then
                con(ncon)=((-0.25d0*dr+0.75d0)*dr-0.75d0)*dr+0.25d0
              else
                con(ncon)=0.
              endif
            endif
          endif
        endif
      enddo
      return
      end
