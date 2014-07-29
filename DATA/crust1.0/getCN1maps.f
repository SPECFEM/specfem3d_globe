      parameter(np=9,nlo=360,nla=180)      

      character*2 cmap*7,dmap*6
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
         write(cmap,'(a,i1)')'map-vp',k
         open(66,file=cmap)
         write(66,104)((vp(k,j,i),i=1,nlo),j=1,nla)         
         close(66)
         write(cmap,'(a,i1)')'map-vs',k
         open(66,file=cmap)
         write(66,104)((vs(k,j,i),i=1,nlo),j=1,nla)         
         close(66)
         write(cmap,'(a,i1)')'map-ro',k
         open(66,file=cmap)
         write(66,104)((rho(k,j,i),i=1,nlo),j=1,nla)         
         close(66)
         write(cmap,'(a,i1)')'map-bd',k
         open(66,file=cmap)
         write(66,105)((bnd(k,j,i),i=1,nlo),j=1,nla)         
         close(66)
      enddo
      do k=1,8
         write(cmap,'(a,i1)')'map-th',k
         open(66,file=cmap)
         do j=1,nla
            do i=1,nlo
               thk(k,j,i)=-(bnd(k+1,j,i)-bnd(k,j,i))
            enddo
         enddo
         write(66,105)((thk(k,j,i),i=1,nlo),j=1,nla)         
         close(66)
      enddo
      write(dmap,'(a)')'sedthk'
      open(64,file=dmap)
      write(dmap,'(a)')'crsthk'
      open(65,file=dmap)
      do j=1,nla
         do i=1,nlo
            ths(j,i)=thk(3,j,i)+thk(4,j,i)+thk(5,j,i)
            thc(j,i)=thk(2,j,i)+ths(j,i)+
     +      thk(6,j,i)+thk(7,j,i)+thk(8,j,i)
         enddo
      enddo
      write(64,105)((ths(j,i),i=1,nlo),j=1,nla)         
      write(65,105)((thc(j,i),i=1,nlo),j=1,nla)         
      close(64)
      close(65)
 104  format(72f5.2)
 105  format(72f7.2)
      end
