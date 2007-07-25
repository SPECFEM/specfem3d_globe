      subroutine vbspl(x,np,xarr,splcon,splcond)
c
c---- this subroutine returns the spline contributions at a particular value of x
c
      dimension xarr(np)
      dimension splcon(np)
      dimension splcond(np)
c
c---- iflag=1 ==>> second derivative is 0 at end points
c---- iflag=0 ==>> first derivative is 0 at end points
c
      iflag=1
c
c---- first, find out within which interval x falls
c
      interval=0
      ik=1
      do while(interval.eq.0.and.ik.lt.np)
        ik=ik+1
        if(x.ge.xarr(ik-1).and.x.le.xarr(ik)) interval=ik-1
      enddo
      if(x.gt.xarr(np)) then
        interval=np
      endif
c
      if(interval.eq.0) then
cc        write(6,"('low value:',2f10.3)") x,xarr(1)
      else if(interval.gt.0.and.interval.lt.np) then
cc        write(6,"('bracket:',i5,3f10.3)") interval,xarr(interval),x,xarr(interval+1)
      else
cc        write(6,"('high value:',2f10.3)") xarr(np),x
      endif
c
      do ib=1,np
        val=0.
        vald=0.
        if(ib.eq.1) then
c
          r1=(x-xarr(1))/(xarr(2)-xarr(1))
          r2=(xarr(3)-x)/(xarr(3)-xarr(1))
          r4=(xarr(2)-x)/(xarr(2)-xarr(1))
          r5=(x-xarr(1))/(xarr(2)-xarr(1))
          r6=(xarr(3)-x)/(xarr(3)-xarr(1))
         r10=(xarr(2)-x)/(xarr(2)-xarr(1))
         r11=(x-xarr(1))  /(xarr(2)-xarr(1))
         r12=(xarr(3)-x)/(xarr(3)-xarr(2))
         r13=(xarr(2)-x)/(xarr(2)-xarr(1))
c
          r1d=1./(xarr(2)-xarr(1))
          r2d=-1./(xarr(3)-xarr(1))
          r4d=-1./(xarr(2)-xarr(1))
          r5d=1./(xarr(2)-xarr(1))
          r6d=-1./(xarr(3)-xarr(1))
         r10d=-1./(xarr(2)-xarr(1))
         r11d=1./(xarr(2)-xarr(1))
         r12d=-1./(xarr(3)-xarr(2))
         r13d=-1./(xarr(2)-xarr(1))
c

          if(interval.eq.ib.or.interval.eq.0) then
               if(iflag.eq.0) then 
                 val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11 +r13**3
                 vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
                 vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
                 vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
                 vald=vald+3.*r13d*r13**2
               else if(iflag.eq.1) then
                 val=0.6667*(r1*r4*r10 + r2*r5*r10 + r2*r6*r11 
     #			+ 1.5*r13**3)
                 vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
                 vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
                 vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
                 vald=vald+4.5*r13d*r13**2
                 vald=0.6667*vald
               endif
          else if(interval.eq.ib+1) then
               if(iflag.eq.0) then
                 val=r2*r6*r12
                 vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
               else if(iflag.eq.1) then
                 val=0.6667*r2*r6*r12
                 vald=0.6667*(r2d*r6*r12+r2*r6d*r12+r2*r6*r12d)
               endif
          else
            val=0.
          endif
c
        else if(ib.eq.2) then
c
          rr1=(x-xarr(1))/(xarr(2)-xarr(1))
          rr2=(xarr(3)-x)/(xarr(3)-xarr(1))
          rr4=(xarr(2)-x)/(xarr(2)-xarr(1))
          rr5=(x-xarr(1))/(xarr(2)-xarr(1))
          rr6=(xarr(3)-x)/(xarr(3)-xarr(1))
         rr10=(xarr(2)-x)/(xarr(2)-xarr(1))
         rr11=(x-xarr(1))  /(xarr(2)-xarr(1))
         rr12=(xarr(3)-x)/(xarr(3)-xarr(2))
c
          rr1d=1./(xarr(2)-xarr(1))
          rr2d=-1./(xarr(3)-xarr(1))
          rr4d=-1./(xarr(2)-xarr(1))
          rr5d=1./(xarr(2)-xarr(1))
          rr6d=-1./(xarr(3)-xarr(1))
         rr10d=-1./(xarr(2)-xarr(1))
         rr11d=1./(xarr(2)-xarr(1))
         rr12d=-1./(xarr(3)-xarr(2))
c
          r1=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
          r3=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
          r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
          r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
          r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
         r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
         r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
         r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))
c
          r1d=1./(xarr(ib+1)-xarr(ib-1))
          r2d=-1./(xarr(ib+2)-xarr(ib-1))
          r3d=1./(xarr(ib)-xarr(ib-1))
          r4d=-1./(xarr(ib+1)-xarr(ib-1))
          r5d=1./(xarr(ib+1)-xarr(ib-1))
          r6d=-1./(xarr(ib+2)-xarr(ib))
          r8d=-1./  (xarr(ib)-xarr(ib-1))
          r9d=1./(xarr(ib)-xarr(ib-1))
         r10d=-1./(xarr(ib+1)-xarr(ib))
         r11d=1./(xarr(ib+1)-xarr(ib))
         r12d=-1./(xarr(ib+2)-xarr(ib+1))
c
          if(interval.eq.ib-1.or.interval.eq.0) then
               val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
               vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
               vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
               vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
               if(iflag.eq.1) then
                 val=val+0.3333*(rr1*rr4*rr10 + rr2*rr5*rr10 +
     #			 rr2*rr6*rr11)
                 vald=vald+0.3333*(rr1d*rr4*rr10+rr1*rr4d*rr10+
     #			rr1*rr4*rr10d)
                 vald=vald+0.3333*(rr2d*rr5*rr10+rr2*rr5d*rr10+
     #			rr2*rr5*rr10d)
                 vald=vald+0.3333*(rr2d*rr6*rr11+rr2*rr6d*rr11+
     #			rr2*rr6*rr11d)
               endif
          else if(interval.eq.ib) then
               val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
               vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
               vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
               vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
               if(iflag.eq.1) then
                 val=val+0.3333*rr2*rr6*rr12
                 vald=vald+0.3333*(rr2d*rr6*rr12+rr2*rr6d*rr12+
     #			rr2*rr6*rr12d)
               endif
          else if(interval.eq.ib+1) then
               val=r2*r6*r12
               vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
          else
               val=0.
          endif
        else if(ib.eq.np-1) then
c
          rr1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          rr2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          rr3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          rr4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          rr5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          rr7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
          rr8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
          rr9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
c
          rr1d=1./(xarr(np)-xarr(np-2))
          rr2d=-1./(xarr(np)-xarr(np-1))
          rr3d=1./(xarr(np)-xarr(np-2))
          rr4d=-1./(xarr(np)-xarr(np-1))
          rr5d=1./(xarr(np)-xarr(np-1))
          rr7d=1./(xarr(np-1)-xarr(np-2))
          rr8d=-1./  (xarr(np)-xarr(np-1))
          rr9d=1./(xarr(np)-xarr(np-1))
c
          r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
          r2=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
          r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r6=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
          r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
          r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
          r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
         r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
         r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
c
          r1d=1./(xarr(ib+1)-xarr(ib-2))
          r2d=-1./(xarr(ib+1)-xarr(ib-1))
          r3d=1./(xarr(ib)-xarr(ib-2))
          r4d=-1./(xarr(ib+1)-xarr(ib-1))
          r5d=1./(xarr(ib+1)-xarr(ib-1))
          r6d=-1./(xarr(ib+1)-xarr(ib))
          r7d=1./(xarr(ib-1)-xarr(ib-2))
          r8d=-1./(xarr(ib)-xarr(ib-1))
          r9d=1./(xarr(ib)-xarr(ib-1))
         r10d=-1./(xarr(ib+1)-xarr(ib))
         r11d=1./(xarr(ib+1)-xarr(ib))
c
          if(interval.eq.ib-2) then
               val=r1*r3*r7
               vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
          else if(interval.eq.ib-1) then
               val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
               vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
               vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
               vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
               if(iflag.eq.1) then
                 val=val+0.3333*rr1*rr3*rr7
                 vald=vald+0.3333*(rr1d*rr3*rr7+rr1*rr3d*rr7+
     #			rr1*rr3*rr7d)
               endif
          else if(interval.eq.ib.or.interval.eq.np) then
               val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
               vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
               vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
               vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
               if(iflag.eq.1) then
                 val=val+0.3333*(rr1*rr3*rr8 + rr1*rr4*rr9 +
     #			 rr2*rr5*rr9)
                 vald=vald+0.3333*(rr1d*rr3*rr8+rr1*rr3d*rr8+
     #			rr1*rr3*rr8d)
                 vald=vald+0.3333*(rr1d*rr4*rr9+rr1*rr4d*rr9+
     #			rr1*rr4*rr9d)
                 vald=vald+0.3333*(rr2d*rr5*rr9+rr2*rr5d*rr9+
     #			rr2*rr5*rr9d)
               endif
          else
            val=0.
          endif
        else if(ib.eq.np) then
c
          r1=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          r2=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          r3=(x-xarr(np-2))/(xarr(np)-xarr(np-2))
          r4=(xarr(np)-x)/(xarr(np)-xarr(np-1))
          r5=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          r7=(x-xarr(np-2))/(xarr(np-1)-xarr(np-2))
          r8=(xarr(np)-x)/  (xarr(np)-xarr(np-1))
          r9=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
          r13=(x-xarr(np-1))/(xarr(np)-xarr(np-1))
c
          r1d=1./(xarr(np)-xarr(np-2))
          r2d=-1./(xarr(np)-xarr(np-1))
          r3d=1./(xarr(np)-xarr(np-2))
          r4d=-1./(xarr(np)-xarr(np-1))
          r5d=1./(xarr(np)-xarr(np-1))
          r7d=1./(xarr(np-1)-xarr(np-2))
          r8d=-1./  (xarr(np)-xarr(np-1))
          r9d=1./(xarr(np)-xarr(np-1))
          r13d=1./(xarr(np)-xarr(np-1))
c
          if(interval.eq.np-2) then
               if(iflag.eq.0) then
                 val=r1*r3*r7
                 vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
               else if(iflag.eq.1) then
                 val=0.6667*r1*r3*r7
                 vald=0.6667*(r1d*r3*r7+r1*r3d*r7+r1*r3*r7d)
               endif
          else if(interval.eq.np-1.or.interval.eq.np) then
               if(iflag.eq.0) then
                 val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9 + r13**3
                 vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
                 vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
                 vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
                 vald=vald+3.*r13d*r13**2
               else if(iflag.eq.1) then
                 val=0.6667*(r1*r3*r8 + r1*r4*r9 + r2*r5*r9 +
     #			 1.5*r13**3)
                 vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
                 vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
                 vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
                 vald=vald+4.5*r13d*r13**2
                 vald=0.6667*vald
               endif
          else
            val=0.
          endif
        else
c
          r1=(x-xarr(ib-2))/(xarr(ib+1)-xarr(ib-2))
          r2=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib-1))
          r3=(x-xarr(ib-2))/(xarr(ib)-xarr(ib-2))
          r4=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib-1))
          r5=(x-xarr(ib-1))/(xarr(ib+1)-xarr(ib-1))
          r6=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib))
          r7=(x-xarr(ib-2))/(xarr(ib-1)-xarr(ib-2))
          r8=(xarr(ib)-x)/  (xarr(ib)-xarr(ib-1))
          r9=(x-xarr(ib-1))/(xarr(ib)-xarr(ib-1))
         r10=(xarr(ib+1)-x)/(xarr(ib+1)-xarr(ib))
         r11=(x-xarr(ib))  /(xarr(ib+1)-xarr(ib))
         r12=(xarr(ib+2)-x)/(xarr(ib+2)-xarr(ib+1))
c
          r1d=1./(xarr(ib+1)-xarr(ib-2))
          r2d=-1./(xarr(ib+2)-xarr(ib-1))
          r3d=1./(xarr(ib)-xarr(ib-2))
          r4d=-1./(xarr(ib+1)-xarr(ib-1))
          r5d=1./(xarr(ib+1)-xarr(ib-1))
          r6d=-1./(xarr(ib+2)-xarr(ib))
          r7d=1./(xarr(ib-1)-xarr(ib-2))
          r8d=-1./  (xarr(ib)-xarr(ib-1))
          r9d=1./(xarr(ib)-xarr(ib-1))
         r10d=-1./(xarr(ib+1)-xarr(ib))
         r11d=1./(xarr(ib+1)-xarr(ib))
         r12d=-1./(xarr(ib+2)-xarr(ib+1))
c
          if(interval.eq.ib-2) then
               val=r1*r3*r7
               vald=r1d*r3*r7+r1*r3d*r7+r1*r3*r7d
          else if(interval.eq.ib-1) then
               val=r1*r3*r8 + r1*r4*r9 + r2*r5*r9
               vald=r1d*r3*r8+r1*r3d*r8+r1*r3*r8d
               vald=vald+r1d*r4*r9+r1*r4d*r9+r1*r4*r9d
               vald=vald+r2d*r5*r9+r2*r5d*r9+r2*r5*r9d
          else if(interval.eq.ib) then
               val=r1*r4*r10 + r2*r5*r10 + r2*r6*r11
               vald=r1d*r4*r10+r1*r4d*r10+r1*r4*r10d
               vald=vald+r2d*r5*r10+r2*r5d*r10+r2*r5*r10d
               vald=vald+r2d*r6*r11+r2*r6d*r11+r2*r6*r11d
          else if(interval.eq.ib+1) then
               val=r2*r6*r12
               vald=r2d*r6*r12+r2*r6d*r12+r2*r6*r12d
          else
            val=0.
          endif
        endif
        splcon(ib)=val
        splcond(ib)=vald
      enddo
      return
      end
