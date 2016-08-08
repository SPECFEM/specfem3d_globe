
  subroutine rd3dmodl(lu,filename,ierror)

  implicit none

  include '3dmodl.h'

  character(len=*) filename

  character(len=128) string
  character(len=128) substr

  integer :: lu,ierror

  integer :: ncoef,i,ihor,iprtlv,ifst,ilst,ifst1,ios,lstr,nmodkern,idummy,nhorpar,lmax

  open(lu,file=filename,iostat=ios)
  if (ios /= 0) then
  stop 'error opening 3-d model'
  endif
  do while (ios == 0)
  read(lu,"(a)",iostat=ios) string
  lstr=len_trim(string)
  if (ios == 0) then
    if (string(1:16) == 'REFERENCE MODEL:') then
      substr=string(17:lstr)
      ifst=1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst=ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'error reading model 1'
      else
        refmodel=substr(ifst:ilst)
      endif
    else if (string(1:11) == 'KERNEL SET:') then
      substr=string(12:len_trim(string))
      ifst=1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst=ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'error reading model 2'
      else
        kernstri=substr(ifst:ilst)
      endif
    else if (string(1:25) == 'RADIAL STRUCTURE KERNELS:') then
      substr=string(26:len_trim(string))
      read(substr,*,iostat=ierror) nmodkern
      if (ierror /= 0) then
        stop 'error reading model 3'
      endif
    else if (string(1:4) == 'DESC' .and. string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      substr=string(10:len_trim(string))
      ifst=1
      ilst=len_trim(substr)
      do while (substr(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst=ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'error reading model 4'
      else
        desckern(idummy)=substr(ifst:ilst)
      endif
    else if (string(1:29) == 'HORIZONTAL PARAMETERIZATIONS:') then
      substr=string(30:len_trim(string))
      read(substr,*,iostat=ierror) nhorpar
      if (ierror /= 0) then
        stop 'error reading model 5'
      endif
    else if (string(1:4) == 'HPAR' .and. string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      ifst=10
      ilst=len_trim(string)
      do while (string(ifst:ifst) == ' ' .and. ifst < ilst)
        ifst=ifst+1
      enddo
      if (ilst-ifst <= 0) then
        stop 'error reading model 6'
      else if (string(ifst:ifst+19) == 'SPHERICAL HARMONICS,') then
        substr=string(20+ifst:len_trim(string))
        read(substr,*) lmax
        ityphpar(idummy)=1
        lmaxhor(idummy)=lmax
        ncoefhor(idummy)=(lmax+1)**2
      else if (string(ifst:ifst+17) == 'SPHERICAL SPLINES,') then
        ifst1=ifst+18
        ifst=len_trim(string)
        ilst=len_trim(string)
        do while(string(ifst:ifst) /= ',')
          ifst=ifst-1
        enddo
        read(string(ifst+1:ilst),*) ncoef
        substr=string(ifst1:ifst-1)
        do while (string(ifst1:ifst1) == ' ' .and. ifst1 < ifst)
          ifst1=ifst1+1
        enddo
        hsplfile(idummy)=string(ifst1:ifst-1)
        ityphpar(idummy)=2
        lmaxhor(idummy)=0
        ncoefhor(idummy)=ncoef
        do i=1,ncoef
          read(lu,*) ixlspl(i,idummy),xlaspl(i,idummy), &
             xlospl(i,idummy),xraspl(i,idummy)
        enddo
      endif
    else if (string(1:4) == 'STRU' .and. string(9:9) == ':') then
      read(string(5:8),"(i4)") idummy
      substr=string(10:len_trim(string))
      read(substr,*) ihor
      ihorpar(idummy)=ihor
      ncoef=ncoefhor(ihor)
      read(lu,"(6e12.4)") (coef(i,idummy),i=1,ncoef)
    endif
  endif
  enddo
  close(lu)

  end subroutine rd3dmodl

