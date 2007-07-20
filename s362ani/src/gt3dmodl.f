      subroutine gt3dmodl(lu,targetfile,
     #    maxhpa,maxker,maxcoe,
     #    numhpa,numker,numcoe,lmxhpa,
     #    ihpakern,itypehpa,coe,
     #    itpspl,xlatspl,xlonspl,radispl,
     #    numvar,ivarkern,varstr,
     #    refmdl,kerstr,hsplfl,dskker,ierror)
c
      character*80 targetfile
      integer numhpa
      integer numker
      dimension numcoe(maxhpa)
      dimension lmxhpa(maxhpa)
      dimension ihpakern(maxker)
      dimension itypehpa(maxhpa)
      dimension coe(maxcoe,maxker)
      dimension itpspl(maxcoe,maxhpa)
      dimension xlatspl(maxcoe,maxhpa)
      dimension xlonspl(maxcoe,maxhpa)
      dimension radispl(maxcoe,maxhpa)
      character*80 refmdl
      character*80 kerstr
      character*80 hsplfl(maxhpa)
      character*40 dskker(maxker)
c
      character*40 string
      character*40 varstr(maxker)
      dimension ivarkern(maxker)
      integer numvar
c
      include '3dmodl.h'
c
      ierror=0
      call rd3dmodl(lu,targetfile,ierror)
c
      if(nhorpar.le.maxhpa) then
        numhpa=nhorpar
      else
        ierror=ierror+1
      endif
c
      if(nmodkern.le.maxker) then
        numker=nmodkern
      else
        ierror=ierror+1
      endif
c
      do i=1,nmodkern
        ihpakern(i)=ihorpar(i)
        dskker(i)=desckern(i)
        do j=1,ncoefhor(ihpakern(i))
          coe(j,i)=coef(j,i)
c          if(j.eq.1) then
c            write(6,"(e12.4)") coe(j,i)
c          endif
        enddo
      enddo
c
      do i=1,nhorpar
        numcoe(i)=ncoefhor(i)
        lmxhpa(i)=lmaxhor(i)
        itypehpa(i)=ityphpar(i)
        if(itypehpa(i).eq.2) then
          do j=1,ncoefhor(i)
            itpspl(j,i)=ixlspl(j,i)
            xlatspl(j,i)=xlaspl(j,i)
            xlonspl(j,i)=xlospl(j,i)
            radispl(j,i)=xraspl(j,i)
          enddo
        endif
        hsplfl(i)=hsplfile(i)
      enddo
c
      numvar=0
      do i=1,nmodkern
        string=dskker(i)
        lstr=lnblnk(string)
        j=1
        do while(string(j:j).ne.','.and.j.lt.lstr)
          j=j+1
        enddo
        ivarkern(i)=0
        do k=1,numvar
          if(string(1:j).eq.varstr(k)(1:j)) then
            ivarkern(i)=k
          endif
        enddo
        if(ivarkern(i).eq.0) then
          numvar=numvar+1
          varstr(numvar)=string(1:j)
          ivarkern(i)=numvar
        endif
      enddo
c
      refmdl=refmodel
      kerstr=kernstri
c
      return
      end

