
  subroutine getcurrmodels(model,imod)

  implicit none

  character(len=*) model
  integer :: imod
  logical exists

  integer :: ios,i

  inquire(file='current_3d_model',exist=exists)
  model=''
  if (exists) then
  open(13,file='current_3d_model')
  ios=0
  do i=1,imod
  if (ios == 0) read(13,"(a)",iostat=ios) model
  enddo
  close(13)
  if (ios == 0) then
  inquire(file=model,exist=exists)
  if (exists) then
  else
    model=''
    endif
  else
  model=''
  endif
  endif

  end subroutine getcurrmodels
