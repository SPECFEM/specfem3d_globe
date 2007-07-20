      subroutine rd3dmodl(lu,filename,ierror)
c
      common /plevel/iprtlv
c
      character*(*) filename
c
      character*128 string
      character*128 substr
c
      include '3dmodl.h'
c
      open(lu,file=filename,iostat=ios)
      if(ios.ne.0) then
        stop 'error opening 3-d model'
      endif
      do while (ios.eq.0)
        read(lu,"(a)",iostat=ios) string
        lstr=lnblnk(string)
        if(ios.eq.0) then
          if(string(1:16).eq.'REFERENCE MODEL:') then
            substr=string(17:lstr)
            ifst=1
            ilst=lnblnk(substr)
            do while (substr(ifst:ifst).eq.' '.and.ifst.lt.ilst)
              ifst=ifst+1
            enddo
            if(ilst-ifst.le.0) then
              stop 'error reading model 1'
            else
              refmodel=substr(ifst:ilst)
            endif
          else if(string(1:11).eq.'KERNEL SET:') then
            substr=string(12:lnblnk(string))
            ifst=1
            ilst=lnblnk(substr)
            do while (substr(ifst:ifst).eq.' '.and.ifst.lt.ilst)
              ifst=ifst+1
            enddo
            if(ilst-ifst.le.0) then
              stop 'error reading model 2'
            else
              kernstri=substr(ifst:ilst)
            endif
          else if(string(1:25).eq.'RADIAL STRUCTURE KERNELS:') then
            substr=string(26:lnblnk(string))
            read(substr,*,iostat=ierror) nmodkern
            if(ierror.ne.0) then
              stop 'error reading model 3'
            endif
          else if(string(1:4).eq.'DESC'.and.string(9:9).eq.':') then
            read(string(5:8),"(i4)") idummy
            substr=string(10:lnblnk(string))
            ifst=1
            ilst=lnblnk(substr)
            do while (substr(ifst:ifst).eq.' '.and.ifst.lt.ilst)
              ifst=ifst+1
            enddo
            if(ilst-ifst.le.0) then
              stop 'error reading model 4'
            else
              desckern(idummy)=substr(ifst:ilst)
            endif
          else if(string(1:29).eq.'HORIZONTAL PARAMETERIZATIONS:') then
            substr=string(30:lnblnk(string))
            read(substr,*,iostat=ierror) nhorpar
            if(ierror.ne.0) then
              stop 'error reading model 5'
            endif
          else if(string(1:4).eq.'HPAR'.and.string(9:9).eq.':') then
            read(string(5:8),"(i4)") idummy
            ifst=10
            ilst=lnblnk(string)
            do while (string(ifst:ifst).eq.' '.and.ifst.lt.ilst)
              ifst=ifst+1
            enddo
            if(ilst-ifst.le.0) then
              stop 'error reading model 6'
            else if(string(ifst:ifst+19).eq.'SPHERICAL HARMONICS,') then
              substr=string(20+ifst:lnblnk(string))
              read(substr,*) lmax
              ityphpar(idummy)=1
              lmaxhor(idummy)=lmax
              ncoefhor(idummy)=(lmax+1)**2
            else if(string(ifst:ifst+17).eq.'SPHERICAL SPLINES,') then
              ifst1=ifst+18
              ifst=lnblnk(string)
              ilst=lnblnk(string)
              do while(string(ifst:ifst).ne.',') 
                ifst=ifst-1
              enddo
              read(string(ifst+1:ilst),*) ncoef
              substr=string(ifst1:ifst-1)
              do while (string(ifst1:ifst1).eq.' '.and.ifst1.lt.ifst)
                ifst1=ifst1+1
              enddo
              hsplfile(idummy)=string(ifst1:ifst-1)
              ityphpar(idummy)=2
              lmaxhor(idummy)=0
              ncoefhor(idummy)=ncoef
              do i=1,ncoef
                read(lu,*) ixlspl(i,idummy),xlaspl(i,idummy),
     #           xlospl(i,idummy),xraspl(i,idummy)
              enddo
            endif
          else if(string(1:4).eq.'STRU'.and.string(9:9).eq.':') then
            read(string(5:8),"(i4)") idummy
            substr=string(10:lnblnk(string))
            read(substr,*) ihor
            ihorpar(idummy)=ihor
            ncoef=ncoefhor(ihor)
            read(lu,"(6e12.4)") (coef(i,idummy),i=1,ncoef)
          endif
        endif
      enddo
      close(lu)
      return
      end
        
