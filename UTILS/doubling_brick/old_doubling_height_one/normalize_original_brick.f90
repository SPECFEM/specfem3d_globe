
 program normalize_brick

! normalize the tripling brick given by Emanuele Casarotti

 implicit none

 integer i

! there are NPOINTS points in the original tripling brick
 integer, parameter :: NPOINTS = 32

 double precision minix,maxix,miniy,maxiy,miniz,maxiz

 double precision, dimension(NPOINTS) :: x,y,z

! read original brick
 do i = 1,NPOINTS
   read(*,*) x(i),y(i),z(i)
 enddo

 minix = minval(x)
 maxix = maxval(x)
 print *,' for x, min max delta = ',minix,maxix,maxix-minix

 miniy = minval(y)
 maxiy = maxval(y)
 print *,' for y, min max delta = ',miniy,maxiy,maxiy-miniy

 miniz = minval(z)
 maxiz = maxval(z)
 print *,' for z, min max delta = ',miniz,maxiz,maxiz-miniz

! write normalized brick
 do i = 1,NPOINTS
   write(*,*) sngl((x(i)-minix)/(maxix-minix)),sngl((y(i)-miniy)/(maxiy-miniy)),sngl((z(i)-miniz)/(maxiz-miniz))
 enddo

 end program normalize_brick

