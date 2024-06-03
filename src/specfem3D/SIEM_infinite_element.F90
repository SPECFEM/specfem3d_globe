!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================


! this module contains infinite-element routines

module siem_infinite_element

  use siem_gll_library, only: kdble

  implicit none

  private

  public :: shape_function_infiniteGLHEX8ZW_GLLR

contains

! TODO: compute 1D lagrange shape function iusing GEN rotuine since we take
! equidistant interpolation points along infinite direction. But now I have
! changed to GLL points so not necessary!
! this subroutine computes GLL (along finite directions) and Radau (along
! infinite direction) quadrature points and weights for 3D

  subroutine shape_function_infiniteGLHEX8ZW_GLLR(ndim,ngllx,nglly,ngllz,ngll,nip, &
                                                  iface,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)

  use siem_gll_library, only: lagrange1dGLLAS,lagrange1dGENAS,zwgljd

  implicit none

  integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,nip,iface
  !of decay
  !integer,parameter :: ngllinf=ngll-nglly*ngllz
  real(kind=kdble),dimension(nip,8),intent(out) :: shape_infinite
  real(kind=kdble),dimension(ndim,nip,8),intent(out) :: dshape_infinite
  real(kind=kdble),dimension(nip,ngll),intent(out) :: lagrange_gl
  real(kind=kdble),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
  real(kind=kdble),intent(out) :: GLw(nip)

  real(kind=kdble),parameter :: jacobi_alpha=0.0_kdble,jacobi_beta=0.0_kdble, &
                                one = 1.0_kdble!,five = 5.0_kdble

  integer :: i,ii,j,k,n,i1,j1,k1,nipx(ndim)
  real(kind=kdble) :: ddir,xi(ndim) !,eta,zeta
  real(kind=kdble),dimension(ngllx) :: gllpx,gllwx,igllpx,igllwx ! GLL points and weights
  real(kind=kdble),dimension(nglly) :: gllpy,gllwy,igllpy,igllwy ! GLL points and weights
  real(kind=kdble),dimension(ngllz) :: gllpz,gllwz,igllpz,igllwz ! GLL points and weights
  real(kind=kdble),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
  real(kind=kdble),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
  integer :: iINF

  iINF = 0
  ddir = one

  if (iface == 1) then
    iINF = 2; ddir = -one
  else if (iface == 2) then
    iINF = 1; ddir = one
  else if (iface == 3) then
    iINF = 2; ddir = one
  else if (iface == 4) then
    iINF = 1; ddir = -one
  else if (iface == 5) then
    iINF = 3; ddir = -one
  else if (iface == 6) then
    iINF = 3; ddir = one
  endif

  nipx(1) = ngllx
  nipx(2) = nglly
  nipx(3) = ngllz

  ! compute everything in indexed order
  ! get GLL points
  ! for alpha=beta=0, jacobi polynomial is legendre polynomial
  ! for ngllx=nglly=ngllz, need to call only once
  ! get GLL points and weights
  call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
  call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
  call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)

  ! integration points are the GLL points
  igllpx = gllpx; igllpy = gllpy; igllpz = gllpz
  igllwx = gllwx; igllwy = gllwy; igllwz = gllwz;

  ! overwrite GLL points/weights with radau counterpart along infinite direction
  if (iINF == 1) call radau_quadrature(ngllx,igllpx,igllwx)
  if (iINF == 2) call radau_quadrature(nglly,igllpy,igllwy)
  if (iINF == 3) call radau_quadrature(ngllz,igllpz,igllwz)

  ii = 0
  do k1 = 1,nipx(3)
    do j1 = 1,nipx(2)
      do i1 = 1,nipx(1)
        ii = ii+1

        ! integration points
        xi(1) = igllpx(i1) !xi,   gll_points(1,ii)
        xi(2) = igllpy(j1) !eta,  gll_points(2,ii)
        xi(3) = igllpz(k1) !zeta, gll_points(3,ii)

        ! integration weights
        GLw(ii) = igllwx(i1)*igllwy(j1)*igllwz(k1)

        call lagrange1dGLLAS(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
        call lagrange1dGLLAS(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
        call lagrange1dGLLAS(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))

        ! interpolation functions
        n = 0
        do k = 1,ngllz
          do j = 1,nglly
            do i = 1,ngllx
              n = n+1
              lagrange_gl(ii,n) = lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
              dlagrange_gl(1,ii,n) = lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
              dlagrange_gl(2,ii,n) = lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
              dlagrange_gl(3,ii,n) = lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
            enddo
          enddo
        enddo

        ! shape functions for HEX8
        ! compute 1d lagrange polynomials
        call lagrange1dGENAS(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
        call lagrange1dGENAS(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
        call lagrange1dGENAS(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))

        ! consider 3 nodes but compute only at 2 nodes
        call lagrange1d_infiniteZWAS(3,xi(iINF),lagrangeINF_x(iINF,:),lagrangeINF_dx(iINF,:))

        n = 0
        do k = 1,2
          do j = 1,2
            do i = 1,2
              n = n+1
              shape_infinite(ii,n) = lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(1,ii,n) = lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(2,ii,n) = lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(3,ii,n) = lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
            enddo
          enddo
        enddo

      enddo
    enddo
  enddo

  end subroutine shape_function_infiniteGLHEX8ZW_GLLR

!
!===========================================
!
! this subroutine computes Gauss quadrature points and weights for 3D
! TODO: compute 1D lagrange shape function iusing GEN rotuine since we take
! equidistant interpolation points along infinite direction.  But now I have
! changed to GLL points so not necessary!

! not used so far...
!
!  subroutine shape_function_infiniteGLHEX8ZW_GQ(ndim,ngllx,nglly,ngllz,ngll,nipx, &
!                                                nip,iface,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)
!
!  use siem_gll_library, only: lagrange1dGLL,lagrange1dGEN,zwgjd,zwgljd
!  implicit none
!  integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,nipx,nip,iface
!  real(kind=kdble),dimension(nip,8),intent(out) :: shape_infinite
!  real(kind=kdble),dimension(ndim,nip,8),intent(out) :: dshape_infinite
!  real(kind=kdble),dimension(nip,ngll),intent(out) :: lagrange_gl
!  real(kind=kdble),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
!  real(kind=kdble),intent(out) :: GLw(nip)
!
!  real(kind=kdble),parameter :: jacobi_alpha=0.0_kdble,jacobi_beta=0.0_kdble,one=1.0_kdble
!
!  integer :: i,ii,j,k,n,i1,j1,k1
!  real(kind=kdble) :: ddir,xi(ndim)
!  real(kind=kdble),dimension(ngllx) :: gllpx,gllwx ! GLL points and weights
!  real(kind=kdble),dimension(nglly) :: gllpy,gllwy ! GLL points and weights
!  real(kind=kdble),dimension(ngllz) :: gllpz,gllwz ! GLL points and weights
!  real(kind=kdble),dimension(nipx) :: igllpx,igllwx ! GLL points and weights
!  real(kind=kdble),dimension(nipx) :: igllpy,igllwy ! GLL points and weights
!  real(kind=kdble),dimension(nipx) :: igllpz,igllwz ! GLL points and weights
!  real(kind=kdble),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
!  real(kind=kdble),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
!  integer :: iINF
!
!  iINF = 0
!  ddir = one
!  if (iface == 1) then
!    iINF = 2; ddir = -one
!  else if (iface == 2) then
!    iINF = 1; ddir = one
!  else if (iface == 3) then
!    iINF = 2; ddir = one
!  else if (iface == 4) then
!    iINF = 1; ddir = -one
!  else if (iface == 5) then
!    iINF = 3; ddir = -one
!  else if (iface == 6) then
!    iINF = 3; ddir = one
!  endif
!
!  ! interpolation points
!  ! compute everything in indexed order
!  ! get GLL points
!  ! for alpha=beta=0, jacobi polynomial is legendre polynomial
!  ! for ngllx=nglly=ngllz, need to call only once
!  call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
!  call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
!  call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)
!
!  ! integration points
!  ! gauss-jacobi or gauss-legendre points and weights
!  call zwgjd(igllpx,igllwx,nipx,jacobi_alpha,jacobi_beta)
!
!  if (nip /= 8.and.nip /= 27) then
!    print *,'ERROR: illegal number of Gauss points:',nip,'!'
!    stop
!  endif
!
!  igllpy = igllpx; igllpz = igllpx;
!  igllwy = igllwx; igllwz = igllwx;
!
!  ii = 0
!  do k1 = 1,nipx
!    do j1 = 1,nipx
!      do i1 = 1,nipx
!        ii = ii+1
!
!        ! integration points
!        xi(1) = igllpx(i1) !xi,   gll_points(1,ii)
!        xi(2) = igllpy(j1) !eta,  gll_points(2,ii)
!        xi(3) = igllpz(k1) !zeta, gll_points(3,ii)
!
!        ! integration weights
!        GLw(ii) = igllwx(i1)*igllwy(j1)*igllwz(k1)
!
!        call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
!        call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
!        call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))
!
!        ! interpolation functions
!        n = 0
!        do k = 1,ngllz
!          do j = 1,nglly
!            do i = 1,ngllx
!              n = n+1
!              lagrange_gl(ii,n) = lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
!              dlagrange_gl(1,ii,n) = lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
!              dlagrange_gl(2,ii,n) = lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
!              dlagrange_gl(3,ii,n) = lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
!            enddo
!          enddo
!        enddo
!
!        ! shape functions for HEX8
!        ! compute 1d lagrange polynomials
!        call lagrange1dGEN(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
!        call lagrange1dGEN(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
!        call lagrange1dGEN(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
!
!        ! consider 3 nodes but compute only at 2 nodes
!        call lagrange1d_infiniteZW(3,xi(iINF),lagrangeINF_x(iINF,:),lagrangeINF_dx(iINF,:))
!
!        n = 0
!        do k = 1,2
!          do j = 1,2
!            do i = 1,2
!              n = n+1
!              shape_infinite(ii,n) = lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
!              dshape_infinite(1,ii,n) = lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
!              dshape_infinite(2,ii,n) = lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
!              dshape_infinite(3,ii,n) = lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
!            enddo
!          enddo
!        enddo
!
!      enddo
!    enddo
!  enddo
!
!  end subroutine shape_function_infiniteGLHEX8ZW_GQ

!
!===========================================
!
! this subroutine extracts the nodes for HEX8 of the finite region of an infinite element

! not used yet...

!  subroutine get_gnodinfHEX8(ndim,ngllx,nglly,ngllz,nginf,iface,gnodinf)
!
!  implicit none
!  integer,intent(in) :: ndim,ngllx,nglly,ngllz,nginf,iface
!  integer,intent(out) :: gnodinf(nginf)
!  integer :: i,j,k,inum
!  integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF
!  real(kind=kdble) :: ddir
!  real(kind=kdble),parameter :: one=1.0_kdble
!
!  if (iface < 1.or.iface > 6) then
!    print *,'ERROR: illegal outer face ID:',iface
!    stop
!  endif
!
!  ! initialize ngllINF indices
!  ngllxINF0 = 1
!  ngllxINF(1) = ngllx
!  ngllxINF(2) = nglly
!  ngllxINF(3) = ngllz
!
!  iINF = 0
!  ddir = one
!
!  if (iface == 1) then
!    iINF = 2; ddir = -one
!  else if (iface == 2) then
!    iINF = 1; ddir = one
!  else if (iface == 3) then
!    iINF = 2; ddir = one
!  else if (iface == 4) then
!    iINF = 1; ddir = -one
!  else if (iface == 5) then
!    iINF = 3; ddir = -one
!  else if (iface == 6) then
!    iINF = 3; ddir = one
!  endif
!
!  if (ddir < 0) then
!    ngllxINF0(iINF) = 2
!  else
!    ngllxINF(iINF) = ngllxINF(iINF)-1
!  endif
!
!  ! extract only the corner nodes
!  inc = ngllxINF - ngllxINF0
!  inum = 0
!  do k = ngllxINF0(3),ngllxINF(3),inc(3)
!    do j = ngllxINF0(2),ngllxINF(2),inc(2)
!      do i = ngllxINF0(1),ngllxINF(1),inc(1)
!        inum = inum+1
!        gnodinf(inum) = nglly*ngllx*(k-1) + ngllx*(j-1)+i
!      enddo
!    enddo
!  enddo
!
!  end subroutine get_gnodinfHEX8

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

! not used yet...
!
!  subroutine lagrange1d_infiniteMO(nenod,nd,xi,ddir,phi,dphi_dxi)
!
!  implicit none
!  integer,intent(in) :: nenod ! number of nodes in an 1d element
!  real(kind=kdble),intent(in) :: nd,xi,ddir ! xi: point where to calculate lagrange function and
!  !its derivative
!  real(kind=kdble),dimension(nenod-1),intent(out) :: phi,dphi_dxi
!  real(kind=kdble),dimension(nenod) :: xii
!  real(kind=kdble) :: fac !dx
!  real(kind=kdble),parameter :: one=1.0_kdble,two=2.0_kdble
!
!  if (nenod /= 3) then
!    print *,'ERROR: infinite element is currently implemented only for 3 nodes!'
!    stop
!  endif
!
!  fac = one/(one-xi)
!
!  phi(1) = -two*xi*fac
!  phi(2) = one - phi(1)
!
!  dphi_dxi(1) = -two*fac*fac
!  dphi_dxi(2) = two*fac*fac
!
!  end subroutine lagrange1d_infiniteMO

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
! Assumed Shape array: pass pointer, subarray or allocatable array

  subroutine lagrange1d_infiniteZWAS(nenod,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  real(kind=kdble),intent(in) :: xi ! xi: point where to calculate lagrange
  !function and its derivative
  real(kind=kdble),dimension(:),intent(out) :: phi,dphi_dxi !,dimension(nenod-1)
  real(kind=kdble) :: fac
  real(kind=kdble),parameter :: one=1.0_kdble

  if (nenod /= 3) then
    stop 'ERROR: infinite element is currently implemented only for 3 nodes!'
  endif

  fac = one/(one-xi)

  phi(1) = -xi*fac
  phi(2) = fac !one-phi(1) !+xi*fac

  dphi_dxi(1) = -fac*fac
  dphi_dxi(2) = fac*fac

  end subroutine lagrange1d_infiniteZWAS

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

! not used ...

!  subroutine lagrange1d_infiniteZW(nenod,xi,phi,dphi_dxi)
!
!  implicit none
!  integer,intent(in) :: nenod ! number of nodes in an 1d element
!  real(kind=kdble),intent(in) :: xi ! xi: point where to calculate lagrange
!  !function and its derivative
!
!  real(kind=kdble),dimension(nenod-1),intent(out) :: phi,dphi_dxi
!  real(kind=kdble) :: fac !dx
!  real(kind=kdble),parameter :: one=1.0_kdble
!
!  if (nenod /= 3) then
!    stop 'ERROR: infinite element is currently implemented only for 3 nodes!'
!  endif
!
!  fac = one/(one-xi)
!
!  phi(1) = -xi*fac
!  phi(2) = fac !one-phi(1) !+xi*fac
!
!  dphi_dxi(1) = -fac*fac
!  dphi_dxi(2) = fac*fac
!
!  end subroutine lagrange1d_infiniteZW

!
!===========================================
!

! Revision:
!   HNG, Apr 19,2012
! RADAU_COMPUTE computes a Radau quadrature rule.
! the Radau rule is distinguished by the fact that the left endpoint
! (-1) is always an abscissa.
!
! the integral:
! integral ( -1 <= x <= 1 ) f(x) dx
!
! the quadrature rule:
! sum ( 1 <= i <= n ) w(i) * f ( x(i) )
!
! the quadrature rule will integrate exactly all polynomials up to
! X**(2*N-2).
!
! Licensing:!
! this code is distributed under the GNU LGPL license.
!
! Modified:
! 06 February 2007
!
! Author:
!    Original MATLAB code by Greg von Winckel.
!    This MATLAB version by John Burkardt.
!
! References:
!    Milton Abramowitz, Irene Stegun,
!    Handbook of Mathematical functions,
!    National Bureau of Standards, 1964,
!    ISBN: 0-486-61272-4,
!    LC: QA47.A34.
!
!    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
!    Spectral Methods in Fluid Dynamics,
!    Springer, 1993,
!    ISNB13: 978-3540522058,
!    LC: QA377.S676.
!
!    Francis Hildebrand,
!    Section 8.11,
!    Introduction to Numerical Analysis,
!    Dover, 1987,
!    ISBN13: 978-0486653631,
!    LC: QA300.H5.
!
!    Arthur Stroud, Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966,
!    LC: QA299.4G3S7.
!
!    Daniel Zwillinger, editor,
!    CRC Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996,
!    ISBN: 0-8493-2479-3,
!    LC: QA47.M315.
!
! Input:
! N: the order or the number of integration points (>0, integer)
! Output:
! X(N): the abscissas
! W(N): the weights

  subroutine radau_quadrature(n,x,w)

  implicit none
  integer,intent(in) :: n
  integer :: i,j
  real(kind=kdble),intent(out) :: x(n),w(n)
  real(kind=kdble) :: rj,rn,fac,p(n,n+1),xold(n)

  real(kind=kdble),parameter :: one=1.0_kdble,pi=3.141592653589793_kdble, &
                                two = 2.0_kdble,tol = 1.0e-12_kdble,zero = 0.0_kdble,zerotol = 1.0e-12_kdble

  if (n < 1) then
    stop 'ERROR: number of quadrature points must be > 1!'
  endif

  x = zero; w = zero
  rn = real(n,kdble)

  ! initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
  fac = two*pi/(two*rn-one)

  ! initialize the Legendre Vandermonde matrix.
  p = zero
  p(2:n,1) = one;
  do i = 1,n
    x(i) = -cos(fac*real(i-1,kdble))
    p(1,i) = (-one)**(i-1)
  enddo
  p(1,n+1) = (-one)**(n)

  ! compute P using the recursion relation.
  ! compute its first and second derivatives and
  ! update X using the Newton-Raphson method.
  xold = two
  do i = 1,100
    if (maxval(abs(x-xold)) <= zerotol) exit
    if (i >= 100) then
      stop 'ERROR: Legendre Vandermonde matrix does not converge!'
    endif

    xold = x;
    p(2:n,2) = x(2:n);
    do j = 2,n
      rj = real(j,kdble)
      p(2:n,j+1) = ((two*rj-one)*x(2:n)*p(2:n,j)+(-rj+one)*p(2:n,j-1))/rj
    enddo
    x(2:n) = xold(2:n)-((one-xold(2:n))/rn)*(p(2:n,n)+p(2:n,n+1))/(p(2:n,n)-p(2:n,n+1))
  enddo

  ! compute the weights.
  w = zero
  w(1) = two/(rn*rn)
  w(2:n) = (one-x(2:n))/(rn*p(2:n,n)*rn*p(2:n,n))

  end subroutine radau_quadrature

end module siem_infinite_element

