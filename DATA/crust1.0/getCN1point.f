      parameter(np=9,nlo=360,nla=180)      

      dimension vp(np,nla,nlo),vs(np,nla,nlo),rho(np,nla,nlo)
      dimension bnd(np,nla,nlo)
        
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
 1    continue
      print*,'enter center lat, long of desired tile (q to quit)' 
      read(*,*,err=99)flat,flon
c make sure longitudes go from -180 to 180
      if(flon.gt.180.)flon=flon-360.
      if(flon.lt.-180.)flon=flon+360.

      ilat=90.-flat+1
      ilon=180.+flon+1
      print 102,'ilat,ilon,crustal type: ',ilat,ilon
 102  format(a,2i4)
      print*,'topography: ',bnd(1,ilat,ilon)
      print*,' layers: vp,vs,rho,bottom'
 103  format(4f7.2)
      do i=1,np-1
         print 103,vp(i,ilat,ilon),vs(i,ilat,ilon),rho(i,ilat,ilon),
     +          bnd(i+1,ilat,ilon) 
      enddo
      print 104,' pn,sn,rho-mantle: ',vp(9,ilat,ilon),
     +         vs(9,ilat,ilon),rho(9,ilat,ilon)
 104  format(a,3f7.2)
      goto 1
 99   continue
      end
