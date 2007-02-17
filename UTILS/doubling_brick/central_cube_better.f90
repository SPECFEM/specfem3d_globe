
  program cubed_sphere
!
!=======================================================================
!
!     "cubed sphere"
!
! ======================================================================
!

  implicit none

! nombre d'elements spectraux dans une direction d'une face, et en radial
  integer, parameter :: nx=20
  integer, parameter :: ny=nx
  integer, parameter :: nz=nx

  integer, parameter :: npoinsurf=5*(nx+1)*(ny+1)
  integer, parameter :: nspecsurf=5*nx*ny

! choisir le rayon de la sphere
  double precision, parameter :: radius_sphere_surf = 1.d0

  integer, parameter :: ngnodside=4

  double precision, parameter :: pi=3.14159265d0, pioverfour=pi/4.d0

  double precision, parameter :: zero=0.d0, one=1.d0

! coef to distort central cube
  double precision, parameter :: COEF1=1.d0
  double precision, parameter :: COEF2=1.35d0/3.d0

  double precision inv_sqrt3

!---- surface

  double precision x1surf(0:nx,0:nx)
  double precision y1surf(0:nx,0:nx)
  double precision z1surf(0:nx,0:nx)

  double precision x2surf(0:nx,0:nx)
  double precision y2surf(0:nx,0:nx)
  double precision z2surf(0:nx,0:nx)

  double precision x3surf(0:nx,0:nx)
  double precision y3surf(0:nx,0:nx)
  double precision z3surf(0:nx,0:nx)

  double precision x4surf(0:nx,0:nx)
  double precision y4surf(0:nx,0:nx)
  double precision z4surf(0:nx,0:nx)

  double precision x5surf(0:nx,0:nx)
  double precision y5surf(0:nx,0:nx)
  double precision z5surf(0:nx,0:nx)

  double precision x6surf(0:nx,0:nx)
  double precision y6surf(0:nx,0:nx)
  double precision z6surf(0:nx,0:nx)

  double precision x1central(0:nx,0:nx)
  double precision y1central(0:nx,0:nx)
  double precision z1central(0:nx,0:nx)

  double precision x2central(0:nx,0:nx)
  double precision y2central(0:nx,0:nx)
  double precision z2central(0:nx,0:nx)

  double precision x3central(0:nx,0:nx)
  double precision y3central(0:nx,0:nx)
  double precision z3central(0:nx,0:nx)

  double precision x4central(0:nx,0:nx)
  double precision y4central(0:nx,0:nx)
  double precision z4central(0:nx,0:nx)

  double precision x5central(0:nx,0:nx)
  double precision y5central(0:nx,0:nx)
  double precision z5central(0:nx,0:nx)

  double precision x6central(0:nx,0:nx)
  double precision y6central(0:nx,0:nx)
  double precision z6central(0:nx,0:nx)

  integer ix,iy
  integer nbdonnees,ispec,numat

  double precision xicoord,etacoord,xival,etaval
  double precision Xbig,Ybig
  double precision delta,deltainvsqrt
  double precision rad_delta_inv_surf

  inv_sqrt3=1.d0/dsqrt(3.d0)

! generation sphere cubique

  do ix=0,nx
  do iy=0,ny

  xicoord  = dble(ix)/dble(nx)
  etacoord = dble(iy)/dble(ny)

  xival = pioverfour * xicoord + ( - pioverfour ) * (one - xicoord)
  etaval = pioverfour * etacoord + ( - pioverfour ) * (one - etacoord)

  Xbig = tan(xival)
  Ybig = tan(etaval)

  delta = one + Xbig*Xbig + Ybig*Ybig
  deltainvsqrt = one / dsqrt(delta)
  rad_delta_inv_surf = radius_sphere_surf * deltainvsqrt

! coordonnees cartesiennes correspondantes

!---- surface

! face 1

  x1surf(ix,iy) = rad_delta_inv_surf
  y1surf(ix,iy) = Xbig * rad_delta_inv_surf
  z1surf(ix,iy) = Ybig * rad_delta_inv_surf

! face 2

  x2surf(ix,iy) = - Xbig * rad_delta_inv_surf
  y2surf(ix,iy) = rad_delta_inv_surf
  z2surf(ix,iy) = Ybig * rad_delta_inv_surf

! face 3

  x3surf(ix,iy) = - rad_delta_inv_surf
  y3surf(ix,iy) = - Xbig * rad_delta_inv_surf
  z3surf(ix,iy) = Ybig * rad_delta_inv_surf

! face 4

  x4surf(ix,iy) = Xbig * rad_delta_inv_surf
  y4surf(ix,iy) = - rad_delta_inv_surf
  z4surf(ix,iy) = Ybig * rad_delta_inv_surf

! face 5

  x5surf(ix,iy) = - Ybig * rad_delta_inv_surf
  y5surf(ix,iy) = Xbig * rad_delta_inv_surf
  z5surf(ix,iy) = rad_delta_inv_surf

! face 6

  x6surf(ix,iy) = - Ybig * rad_delta_inv_surf
  y6surf(ix,iy) = - Xbig * rad_delta_inv_surf
  z6surf(ix,iy) = - rad_delta_inv_surf

!!! DK DK central cube

! face 1

  x1central(ix,iy) = inv_sqrt3
  y1central(ix,iy) = Xbig * inv_sqrt3
  z1central(ix,iy) = Ybig * inv_sqrt3

! face 2

  x2central(ix,iy) = - Xbig * inv_sqrt3
  y2central(ix,iy) = inv_sqrt3
  z2central(ix,iy) = Ybig * inv_sqrt3

! face 3

  x3central(ix,iy) = - inv_sqrt3
  y3central(ix,iy) = - Xbig * inv_sqrt3
  z3central(ix,iy) = Ybig * inv_sqrt3

! face 4

  x4central(ix,iy) = Xbig * inv_sqrt3
  y4central(ix,iy) = - inv_sqrt3
  z4central(ix,iy) = Ybig * inv_sqrt3

! face 5

  x5central(ix,iy) = - Ybig * inv_sqrt3
  y5central(ix,iy) = Xbig * inv_sqrt3
  z5central(ix,iy) = inv_sqrt3

! face 6

  x6central(ix,iy) = - Ybig * inv_sqrt3
  y6central(ix,iy) = - Xbig * inv_sqrt3
  z6central(ix,iy) = - inv_sqrt3

  enddo
  enddo

! representation sous AVS du volume

  open(unit=20,file='cubed_sphere_surface.inp',status='unknown')

! nb de noeuds, de cellules, de donnees par cellule
  nbdonnees = 1
  write(20,*) npoinsurf,nspecsurf,' 1 0 0'

! numero et coordonnees des noeuds
  do ix=0,nx
  do iy=0,ny

    write(20,1200) '1',ix,iy, &
      x1central(ix,iy)+COEF1*(x1surf(ix,iy)-x1central(ix,iy)), &
      y1central(ix,iy)+COEF1*(y1surf(ix,iy)-y1central(ix,iy)), &
      z1central(ix,iy)+COEF1*(z1surf(ix,iy)-z1central(ix,iy))
    write(20,1200) '2',ix,iy, &
      x2central(ix,iy)+COEF2*(x2surf(ix,iy)-x2central(ix,iy)), &
      y2central(ix,iy)+COEF2*(y2surf(ix,iy)-y2central(ix,iy)), &
      z2central(ix,iy)+COEF2*(z2surf(ix,iy)-z2central(ix,iy))
    write(20,1200) '3',ix,iy, &
      x3central(ix,iy)+COEF2*(x3surf(ix,iy)-x3central(ix,iy)), &
      y3central(ix,iy)+COEF2*(y3surf(ix,iy)-y3central(ix,iy)), &
      z3central(ix,iy)+COEF2*(z3surf(ix,iy)-z3central(ix,iy))
    write(20,1200) '4',ix,iy, &
      x4central(ix,iy)+COEF2*(x4surf(ix,iy)-x4central(ix,iy)), &
      y4central(ix,iy)+COEF2*(y4surf(ix,iy)-y4central(ix,iy)), &
      z4central(ix,iy)+COEF2*(z4surf(ix,iy)-z4central(ix,iy))
    write(20,1200) '5',ix,iy, &
      x5central(ix,iy)+COEF2*(x5surf(ix,iy)-x5central(ix,iy)), &
      y5central(ix,iy)+COEF2*(y5surf(ix,iy)-y5central(ix,iy)), &
      z5central(ix,iy)+COEF2*(z5surf(ix,iy)-z5central(ix,iy))

  enddo
  enddo

! numero et coordonnees des cellules, definition materiau
  ispec = 0
  numat = 1
  do ix=0,nx-1
  do iy=0,ny-1
      ispec = ispec + 1
      write(20,1310) ispec,numat,'1',ix,iy,'1',ix+1,iy,'1',ix+1,iy+1,'1',ix,iy+1
      write(20,1310) ispec,numat,'2',ix,iy,'2',ix+1,iy,'2',ix+1,iy+1,'2',ix,iy+1
      write(20,1310) ispec,numat,'3',ix,iy,'3',ix+1,iy,'3',ix+1,iy+1,'3',ix,iy+1
      write(20,1310) ispec,numat,'4',ix,iy,'4',ix+1,iy,'4',ix+1,iy+1,'4',ix,iy+1
      write(20,1310) ispec,numat,'5',ix,iy,'5',ix+1,iy,'5',ix+1,iy+1,'5',ix,iy+1
  enddo
  enddo

! structure des donnees aux noeuds
  write(20,*) '1 1'
  write(20,*) 'Zcoord, meters'

! donnees aux noeuds
  do ix=0,nx
  do iy=0,ny
      write(20,1210) '1',ix,iy,z1surf(ix,iy)
      write(20,1210) '2',ix,iy,z2surf(ix,iy)
      write(20,1210) '3',ix,iy,z3surf(ix,iy)
      write(20,1210) '4',ix,iy,z4surf(ix,iy)
      write(20,1210) '5',ix,iy,z5surf(ix,iy)
  enddo
  enddo

  close(20)

 1200  format(a1,i3.3,i3.3,1x,e12.5,1x,e12.5,1x,e12.5)
 1210  format(a1,i3.3,i3.3,1x,e12.5)
 1310  format(i6,1x,i1,' quad ',a1,i3.3,i3.3,1x,a1,i3.3,i3.3,1x, &
              a1,i3.3,i3.3,1x,a1,i3.3,i3.3)

  end program cubed_sphere

