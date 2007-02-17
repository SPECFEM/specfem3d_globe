
 program create_new_brick

! create new doubling brick based upon modification of Emanuele Casarotti's tripling brick

! to execute it:  xcreate_new_brick < original_points_emanuele_DX.txt

 implicit none

! there are NPOINTS points in the original tripling brick
 integer, parameter :: NPOINTS = 32

 integer i

 double precision minix,maxix,miniy,maxiy,miniz,maxiz,xnew,ynew,znew

 double precision, dimension(NPOINTS) :: x,y,z

! read original brick
 do i = 1,NPOINTS
   read(*,*) x(i),y(i),z(i)
 enddo

 minix = minval(x)
 maxix = maxval(x)
 print *,' for x, min max delta before = ',minix,maxix,maxix-minix

 miniy = minval(y)
 maxiy = maxval(y)
 print *,' for y, min max delta before = ',miniy,maxiy,maxiy-miniy

 miniz = minval(z)
 maxiz = maxval(z)
 print *,' for z, min max delta before = ',miniz,maxiz,maxiz-miniz

! create new brick by truncating X and Y in Xmin and Ymin
 do i = 1,NPOINTS
   xnew = (x(i)-minix)/(maxix-minix)
   ynew = (y(i)-miniy)/(maxiy-miniy)
   znew = (z(i)-miniz)/(maxiz-miniz)
   if(xnew < 0.001d0) xnew = 0.3333333333333d0
   if(ynew < 0.001d0) ynew = 0.3333333333333d0

! change Z from (O, 0.50, 0.75, 1) to (0, 0.33, 0.66, 1) to have more even elements
   if(abs(znew - 0.50d0) < 0.001d0) znew = 1.d0 / 3.d0
   if(abs(znew - 0.75d0) < 0.001d0) znew = 2.d0 / 3.d0

   x(i) = xnew
   y(i) = ynew
   z(i) = znew
 enddo

 minix = minval(x)
 maxix = maxval(x)
 print *,' for x, min max delta after = ',minix,maxix,maxix-minix

 miniy = minval(y)
 maxiy = maxval(y)
 print *,' for y, min max delta after = ',miniy,maxiy,maxiy-miniy

 miniz = minval(z)
 maxiz = maxval(z)
 print *,' for z, min max delta after = ',miniz,maxiz,maxiz-miniz

! write new brick
 do i = 1,NPOINTS
   write(*,*) sngl((x(i)-minix)/(maxix-minix)),sngl((y(i)-miniy)/(maxiy-miniy)),sngl((z(i)-miniz)/(maxiz-miniz))
 enddo

 end program create_new_brick

