      subroutine putcurrmodels(model,imod)
      character*(*) model
      character*80 string(10)
      logical exists
      open(13,file='current_3d_model')
      ios=0
      do i=1,10
	string(i)=''
	if(ios.eq.0) read(13,"(a)",iostat=ios) string(i)
      enddo
      string(imod)=model
      rewind(13)
      do i=1,10
	write(13,"(a)") string(i)
      enddo
      close(13)
      return
      end
