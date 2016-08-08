
  subroutine gt3dmodl(lu,targetfile, &
      maxhpa,maxker,maxcoe, &
      numhpa,numker,numcoe,lmxhpa, &
      ihpakern,itypehpa,coe, &
      itpspl,xlatspl,xlonspl,radispl, &
      numvar,ivarkern,varstr, &
      refmdl,kerstr,hsplfl,dskker,ierror)

  implicit none

  include '3dmodl.h'

  character(len=80) targetfile

  integer numhpa,numker,maxhpa,maxker,maxcoe

  integer numcoe(maxhpa)
  integer lmxhpa(maxhpa)
  integer ihpakern(maxker)
  integer itypehpa(maxhpa)
  integer itpspl(maxcoe,maxhpa)
  integer ivarkern(maxker)

  real(kind=4) coe(maxcoe,maxker)
  real(kind=4) xlatspl(maxcoe,maxhpa)
  real(kind=4) xlonspl(maxcoe,maxhpa)
  real(kind=4) radispl(maxcoe,maxhpa)

  character(len=80) refmdl
  character(len=80) kerstr
  character(len=80) hsplfl(maxhpa)
  character(len=40) dskker(maxker)
  character(len=40) string
  character(len=40) varstr(maxker)

  integer numvar,ierror,lu,nhorpar,nmodkern,i,j,lstr,k

  ierror=0
  call rd3dmodl(lu,targetfile,ierror)

  if (nhorpar <= maxhpa) then
  numhpa=nhorpar
  else
  ierror=ierror+1
  endif

  if (nmodkern <= maxker) then
  numker=nmodkern
  else
  ierror=ierror+1
  endif

  do i=1,nmodkern
  ihpakern(i)=ihorpar(i)
  dskker(i)=desckern(i)
  do j=1,ncoefhor(ihpakern(i))
    coe(j,i)=coef(j,i)
!          if (j==1) then
!            write(6,"(e12.4)") coe(j,i)
!          endif
  enddo
  enddo

  do i=1,nhorpar
  numcoe(i)=ncoefhor(i)
  lmxhpa(i)=lmaxhor(i)
  itypehpa(i)=ityphpar(i)
  if (itypehpa(i) == 2) then
    do j=1,ncoefhor(i)
      itpspl(j,i)=ixlspl(j,i)
      xlatspl(j,i)=xlaspl(j,i)
      xlonspl(j,i)=xlospl(j,i)
      radispl(j,i)=xraspl(j,i)
    enddo
  endif
  hsplfl(i)=hsplfile(i)
  enddo

  numvar=0
  do i=1,nmodkern
  string=dskker(i)
  lstr=len_trim(string)
  j=1
  do while(string(j:j) /= ',' .and. j < lstr)
    j=j+1
  enddo
  ivarkern(i)=0
  do k=1,numvar
    if (string(1:j) == varstr(k)(1:j)) then
      ivarkern(i)=k
    endif
  enddo
  if (ivarkern(i) == 0) then
    numvar=numvar+1
    varstr(numvar)=string(1:j)
    ivarkern(i)=numvar
  endif
  enddo

  refmdl=refmodel
  kerstr=kernstri

  end subroutine gt3dmodl

