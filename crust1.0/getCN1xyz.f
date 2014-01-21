      parameter(np=9,nlo=360,nla=180)      

      character*2 cmap*7,dmap*10
      dimension vp(np,nla,nlo),vs(np,nla,nlo),rho(np,nla,nlo)
      dimension bnd(np,nla,nlo)
      dimension thk(np,nla,nlo),thc(nla,nlo),ths(nla,nlo)
        
      open(51,file='crust1.vp')
      open(52,file='crust1.vs')
      open(53,file='crust1.rho')
      open(54,file='crust1.bnds')

      print*,' .... reading all maps ... ' 

      do j=1,nla
         do i=1,nlo
            read(51,*)(vp(k,j,i),k=1,np)
            read(52,*)(vs(k,j,i),k=1,np)
            read(53,*)(rho(k,j,i),k=1,np)
            read(54,*)(bnd(k,j,i),k=1,np)
         enddo
      enddo
      close(51)
      close(52)
      close(53)
      close(54)

c write out maps
      do k=1,9
         write(cmap,'(a,i1)')'xyz-vp',k
         open(66,file=cmap)
         do j=1,nla
         flat=90.-j+0.5
         do i=1,nlo
            flon=-180.+i-0.5
            write(66,104)flon,flat,vp(k,j,i)
         enddo
         enddo
         close(66)
         write(cmap,'(a,i1)')'xyz-vs',k
         open(66,file=cmap)
         do j=1,nla
         flat=90.-j+0.5
         do i=1,nlo
            flon=-180.+i-0.5
            write(66,104)flon,flat,vs(k,j,i)
         enddo
         enddo
         close(66)
         write(cmap,'(a,i1)')'xyz-ro',k
         open(66,file=cmap)
         do j=1,nla
         flat=90.-j+0.5
         do i=1,nlo
            flon=-180.+i-0.5
            write(66,104)flon,flat,rho(k,j,i)
         enddo
         enddo
         close(66)
         write(cmap,'(a,i1)')'xyz-bd',k
         open(66,file=cmap)
         do j=1,nla
         flat=90.-j+0.5
         do i=1,nlo
            flon=-180.+i-0.5
            write(66,105)flon,flat,bnd(k,j,i)
         enddo
         enddo
         close(66)
      enddo
      do k=1,8
         write(cmap,'(a,i1)')'xyz-th',k
         open(66,file=cmap)
         do j=1,nla
         flat=90.-j+0.5
         do i=1,nlo
            flon=-180.+i-0.5
            thk(k,j,i)=-(bnd(k+1,j,i)-bnd(k,j,i))
            write(66,105)flon,flat,thk(k,j,i)
         enddo
         enddo
         close(66)
      enddo
      write(dmap,'(a)')'sedthk.xyz'
      open(64,file=dmap)
      write(dmap,'(a)')'crsthk.xyz'
      open(65,file=dmap)
      do j=1,nla
         flat=90.-j+0.5
         do i=1,nlo
            flon=-180.+i-0.5
            ths(j,i)=thk(3,j,i)+thk(4,j,i)+thk(5,j,i)
            thc(j,i)=thk(2,j,i)+ths(j,i)+
     +      thk(6,j,i)+thk(7,j,i)+thk(8,j,i)
            write(64,105)flon,flat,ths(j,i)
            write(65,105)flon,flat,thc(j,i)
         enddo
      enddo
      close(64)
      close(65)
 104  format(2f6.1,1x,f5.2)
 105  format(2f6.1,1x,f7.2)
      end
