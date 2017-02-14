
  subroutine putcurrmodels(model,imod)

  implicit none

  integer :: imod

  character(len=*) model
  character(len=80) string(10)

  integer :: i,ios

  open(13,file='current_3d_model')
  ios=0
  do i=1,10
  string(i)=''
  if (ios == 0) read(13,"(a)",iostat=ios) string(i)
  enddo
  string(imod)=model
  rewind(13)
  do i=1,10
  write(13,"(a)") string(i)
  enddo
  close(13)

  end subroutine putcurrmodels

