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

! TODO: full gravity is not working yet, needs to fully implement solver...
#ifdef USE_PETSC_NOT_WORKING_YET

! this module contains infinite-element routines
! REVISION
!   HNG, Apr 11,2012; HNG, Jul 12,2011; HNG, Apr 09,2010

module infinite_element

  integer,parameter :: kdble=selected_real_kind(15) ! double precision

contains

  !! this subroutine add a layer of infinite mesh outside the model given the
  !! reference surface and infinite element information.
  !! REVISION
  !!   HNG, Jul 12,2011; ; HNG, Apr 09,2010
  !subroutine add_infmesh(ismpi,myid,nproc,errcode,errtag)
  !use global
  !use math_constants, only: zero
  !use math_library, only: distance,i_uniinv
  !implicit none
  !logical,intent(in) :: ismpi
  !integer,intent(in) :: myid,nproc
  !integer,intent(out) :: errcode
  !character(len=250),intent(out) :: errtag
  !integer :: bctype,i,ios
  !integer :: ielmt,iface,nelpart,i_elpart
  !integer,dimension(6,4) :: node_face ! local node numbers in each face
  !
  !real(kind=kdble) :: gaminf,x0(ndim),val
  !real(kind=kdble),parameter :: one=1.0_kdble
  !integer :: g_numOLD(8,nelmt),mat_idOLD(nelmt),nelmtOLD,nnodeOLD
  !real(kind=kdble) :: r1,g_coordOLD(ndim,nnode)
  !real(kind=kdble),allocatable :: xs(:,:),mirxs(:,:)
  !
  !integer :: n1,n2,nelmtINF,nnode_inf,nsnode,nsnode_all
  !integer,allocatable :: nodelist(:),inode_order(:),g_numinf(:),iface_elmt(:)
  !logical,allocatable :: isnode(:)
  !
  !character(len=20) :: format_str,ptail
  !character(len=250) :: fname
  !character(len=150) :: data_path,strline
  !
  !integer :: ipart ! partition ID
  !
  !errtag="ERROR: unknown!"
  !errcode=-1
  !! set data path
  !if (ismpi) then
  !  data_path=trim(part_path)
  !else
  !  data_path=trim(inp_path)
  !endif
  !
  !ipart=myid-1 ! partition ID starts from 0
  !if (ismpi) then
  !  write(format_str,*)ceiling(log10(real(nproc)+1))
  !  format_str='(a,i'//trim(adjustl(format_str))//'.'//trim(adjustl(format_str))//')'
  !  write(ptail, fmt=format_str)'_proc',ipart
  !else
  !  ptail=""
  !endif
  !
  !! local node numbering in each face CUBIT/EXODUS convention
  !node_face(1,:)=(/1,2,6,5/) ! counterclockwise w.r.t. outer normal
  !node_face(2,:)=(/2,3,7,6/) ! counterclockwise w.r.t. outer normal
  !node_face(3,:)=(/4,3,7,8/) ! clockwise w.r.t. outer normal
  !node_face(4,:)=(/1,4,8,5/) ! clockwise w.r.t. outer normal
  !node_face(5,:)=(/1,2,3,4/) ! clockwise w.r.t. outer normal
  !node_face(6,:)=(/5,6,7,8/) ! counterclockwise w.r.t. outer normal
  !
  !fname=trim(data_path)//trim(infrfile)//trim(ptail)
  !!print *,fname
  !open(unit=11,file=trim(fname),status='old',action='read',iostat = ios)
  !if (ios /= 0) then
  !  write(errtag,*)'ERROR: file "'//trim(fname)//'" cannot be opened!'
  !  return
  !endif
  !
  !read(11,*,iostat=ios)nelpart
  !if (ios /= 0) then
  !  write(errtag,*)'ERROR: bad infrfile!'
  !  return
  !endif
  !nelmtINF=nelpart
  !allocate(iface_elmt(nelmtINF))
  !nsnode_all=4*nelmtINF
  !allocate(nodelist(nsnode_all),inode_order(nsnode_all))
  !n1=1; n2=4
  !do i_elpart=1,nelpart
  !  !print *,n1,n2
  !  read(11,*)ielmt,iface ! This will read a line and proceed to next line
  !  iface_elmt(i_elpart)=iface
  !  nodelist(n1:n2)=g_num(node_face(iface,:),ielmt)
  !  n1=n2+1; n2=n1+3
  !enddo
  !close(11)
  !!stop
  !!print *,n2,nsnode_all;! stop
  !call i_uniinv(nodelist,inode_order)
  !
  !!print *,inode_order;
  !!print *,minval(nodelist),maxval(nodelist);
  !!print *,minval(abs(inode_order)),maxval(abs(inode_order)); stop
  !
  !nsnode=maxval(inode_order)
  !allocate(isnode(nsnode),xs(ndim,nsnode),mirxs(ndim,nsnode),g_numinf(nsnode))
  !
  !isnode=.false.
  !! assign xs
  !xs(:,inode_order(1))=g_coord(:,nodelist(1))
  !isnode(inode_order(1))=.true.
  !do i=2,nsnode_all
  !  if (.not.isnode(inode_order(i))) then
  !     xs(:,inode_order(i))=g_coord(:,nodelist(i))
  !     isnode(inode_order(i))=.true.
  !  endif
  !enddo
  !deallocate(isnode)
  !!print *,'Hi:',inod,nsnode_all
  !!stop
  !! pole specific to the spherical body which has the center at (0,0,0)
  !x0=zero
  !
  !! compute mirror nodes
  !do i=1,nsnode
  !  r1=distance(x0,xs(:,i),ndim)
  !  if (rinf <= r1) then
  !    write(errtag,*)'ERROR: reference infinite radius is smaller than the model!'
  !    return
  !  endif
  !  gaminf=r1/(rinf-r1)
  !  !print *,one+one/gaminf
  !  ! division formula
  !  mirxs(:,i)=((gaminf+one)*xs(:,i)-x0)/gaminf
  !  g_numinf(i)=nnode+i
  !enddo
  !!stop
  !deallocate(xs)
  !g_numOLD=g_num;
  !g_coordOLD=g_coord;
  !mat_idOLD=mat_id
  !deallocate(g_num,g_coord,mat_id)
  !
  !nelmtOLD=nelmt; nnodeOLD=nnode
  !
  !nelmt=nelmtOLD+nelmtINF
  !nnode=nnodeOLD+nsnode
  !
  !! reallocate global node - and element-arrays
  !allocate(g_num(8,nelmt),g_coord(ndim,nnode),mat_id(nelmt))
  !
  !! update connectivity, coordinates, and material IDs
  !g_num(:,1:nelmtOLD)=g_numOLD
  !g_coord(:,1:nnodeOLD)=g_coordOLD
  !mat_id(1:nelmtOLD)=mat_idOLD
  !! add to material list
  !mat_id(nelmtOLD+1:nelmt)=mat_idINF
  !!print *,mat_id,'infmat:',mat_idINF; stop
  !! add to global node list
  !g_coord(:,nnodeOLD+1:nnode)=mirxs
  !deallocate(mirxs)
  !
  !! add to global element list
  !! to have face 5 of new element always clockwise w.r.t. outer normal
  !! face 3 or 4 or 5 of reference element has to be reordered
  !n1=1; n2=4
  !do i=1,nelmtINF
  !  if (iface_elmt(i)==3.or.iface_elmt(i)==4.or.iface_elmt(i)==5) then
  !    g_num(1:4,nelmtOLD+i)=nodelist(n2:n1:-1)
  !    g_num(5:8,nelmtOLD+i)=g_numinf(inode_order(n2:n1:-1))
  !  else
  !    g_num(1:4,nelmtOLD+i)=nodelist(n1:n2)
  !    g_num(5:8,nelmtOLD+i)=g_numinf(inode_order(n1:n2))
  !  endif
  !  n1=n2+1; n2=n1+3
  !enddo
  !deallocate(nodelist,inode_order,g_numinf,iface_elmt)
  !
  !ielmtINF1=nelmtOLD+1; ielmtINF2=nelmt
  !ifaceINF=6
  !!sync all
  !!call stop_all()
  !
  !! compute nodal to global
  !errcode=0
  !return
  !
  !end subroutine add_infmesh

!
!===========================================
!

! TODO: compute 1D lagrange shape function iusing GEN rotuine since we take
! equidistant interpolation points along infinite direction. But now I have
! changed to GLL points so not necessary!
! this subroutine computes GLL (along finite directions) and Radau (along
! infinite direction) quadrature points and weights for 3D

  subroutine shape_function_infiniteGLHEX8ZW_GLLR(ndim,ngllx,nglly,ngllz,ngll,nip, &
                                                  iface,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)

  use gll_library1, only: lagrange1dGLLAS,lagrange1dGENAS,zwgljd
  implicit none
  integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,nip,iface
  !of decay
  !integer,parameter :: ngllinf=ngll-nglly*ngllz
  real(kind=kdble),intent(in) :: gam,a!,nd !nd: order
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
  !real(kind=kdble),dimension(nglly) :: lagrange_y,lagrange_dy
  !real(kind=kdble),dimension(ngllz) :: lagrange_z,lagrange_dz
  integer :: iINF !,ngllxINF(ndim)
  !print *,iface; stop
  !ngllxINF(1)=ngllx
  !ngllxINF(2)=nglly
  !ngllxINF(3)=ngllz
  !print *,nip; stop
  !print *,size(shape_infinite),size(dshape_infinite)
  !stop
  ddir = one
  if (iface == 1) then
    iINF = 2; ddir=-one
  else if (iface == 2) then
    iINF = 1; ddir = one
  else if (iface == 3) then
    iINF = 2; ddir = one
  else if (iface == 4) then
    iINF = 1; ddir=-one
  else if (iface == 5) then
    iINF = 3; ddir=-one
  else if (iface == 6) then
    iINF = 3; ddir = one
  endif
  !print *,ddir; stop
  !ngllxINF(iINF)=ngllxINF(iINF)-1

  nipx(1)=ngllx
  nipx(2)=nglly
  nipx(3)=ngllz
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

  !print *,gllpz
  !print *,gllwz
  !stop
  ! gauss-jacobi or gauss-legendre points and weights
  !call zwgjd(gllpx,gllwx,nipx,jacobi_alpha,jacobi_beta)
  !print *,nipx
  !print *,gllpx
  !print *,gllwx
  ! for an infinite element we use Gauss-Legendre quadrature
  !if (nip==8) then
  !  gllpx(1)=-one/sqrt(3.0_kdble); gllpx(2)=-gllpx(1)
  !  gllwx(1)=one; gllwx(2)=one;
  !else if (nip==27) then
  !  gllpx(1)=-sqrt(3.0_kdble/five); gllpx(2)=0.0_kdble; gllpx(3)=-gllpx(1)
  !  gllwx(1)=five/9.0_kdble; gllwx(2)=8.0_kdble/9.0_kdble; gllwx(3)=gllwx(1);
  !else
  !if (nip /= 8.and.nip /= 27) then
  !  print *,'ERROR: illegal number of Gauss points:',nip,'!'
  !  stop
  !endif
  !print *,gllpx
  !print *,gllwx
  !stop
  !gllpy=gllpx; gllpz=gllpx;
  !gllwy=gllwx; gllwz=gllwx;

  !print *,gllpz
  !print *,gllwz
  !stop
  ii = 0
  do k1 = 1,nipx(3)
    do j1 = 1,nipx(2)
      do i1 = 1,nipx(1)
        ii = ii+1
        !do ii=1,ngll ! ngllx*nglly*ngllz

        ! integration points
        xi(1)=igllpx(i1) !xi,   gll_points(1,ii)
        xi(2)=igllpy(j1) !eta,  gll_points(2,ii)
        xi(3)=igllpz(k1) !zeta, gll_points(3,ii)

        !xi(iINF)=ddir*xi(iINF)

        ! integration weights
        GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)

        call lagrange1dGLLAS(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
        call lagrange1dGLLAS(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
        call lagrange1dGLLAS(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))

        ! interpolation functions
        n = 0
        do k = 1,ngllz
          do j = 1,nglly
            do i = 1,ngllx
              n = n+1
              lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
              dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
              dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
              dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
            enddo
          enddo
        enddo

        ! shape functions for HEX8
        ! compute 1d lagrange polynomials
        call lagrange1dGENAS(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
        call lagrange1dGENAS(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
        call lagrange1dGENAS(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
        ! consider 3 nodes but compute only at 2 nodes
        !call lagrange_infinite(3,nd,xi(iINF),ddir,gam,a,lagrangeINF_x(iINF,:),lagrangeINF_dx(iINF,:))
        call lagrange1d_infiniteZWAS(3,xi(iINF),lagrangeINF_x(iINF,:),lagrangeINF_dx(iINF,:))
        !call lagrange_infinite(ngllx,nd,1.0_kdble,gam,a,lagrangeINF_x,lagrangeINF_dx)
        !print *,lagrangeINF_x,lagrangeINF_dx; stop
        !call lagrange1d_infinite(ngllx,xi,lagrangeINF_x,lagrangeINF_dx)
        !call lagrange1d_infiniteZW(ngllx,xi,lagrangeINF_x,lagrangeINF_dx)
        n = 0
        do k = 1,2
          do j = 1,2
            do i = 1,2
              n = n+1
              shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
            enddo
          enddo
        enddo

        !enddo
      enddo
    enddo
  enddo
  !stop 'Hi!'
  return
  end subroutine shape_function_infiniteGLHEX8ZW_GLLR

!
!===========================================
!
! this subroutine computes Gauss quadrature points and weights for 3D
! TODO: compute 1D lagrange shape function iusing GEN rotuine since we take
! equidistant interpolation points along infinite direction.  But now I have
! changed to GLL points so not necessary!

  subroutine shape_function_infiniteGLHEX8ZW_GQ(ndim,ngllx,nglly,ngllz,ngll,nipx, &
                                                nip,iface,nd,gam,a,shape_infinite,dshape_infinite,lagrange_gl,dlagrange_gl,GLw)

  use gll_library1, only: lagrange1dGLL,lagrange1dGEN,zwgjd,zwgljd
  implicit none
  integer,intent(in) :: ndim,ngllx,nglly,ngllz,ngll,nipx,nip,iface
  real(kind=kdble),intent(in) :: gam,a,nd !nd: order
  real(kind=kdble),dimension(nip,8),intent(out) :: shape_infinite
  real(kind=kdble),dimension(ndim,nip,8),intent(out) :: dshape_infinite
  real(kind=kdble),dimension(nip,ngll),intent(out) :: lagrange_gl
  real(kind=kdble),dimension(ndim,nip,ngll),intent(out) :: dlagrange_gl
  real(kind=kdble),intent(out) :: GLw(nip)
  real(kind=kdble),parameter :: jacobi_alpha=0.0_kdble,jacobi_beta=0.0_kdble,one=1.0_kdble
  integer :: i,ii,j,k,n,i1,j1,k1
  real(kind=kdble) :: ddir,xi(ndim),tmp !,eta,zeta
  real(kind=kdble),dimension(ngllx) :: gllpx,gllwx ! GLL points and weights
  real(kind=kdble),dimension(nglly) :: gllpy,gllwy ! GLL points and weights
  real(kind=kdble),dimension(ngllz) :: gllpz,gllwz ! GLL points and weights
  real(kind=kdble),dimension(nipx) :: igllpx,igllwx ! GLL points and weights
  real(kind=kdble),dimension(nipx) :: igllpy,igllwy ! GLL points and weights
  real(kind=kdble),dimension(nipx) :: igllpz,igllwz ! GLL points and weights
  real(kind=kdble),dimension(ndim,ngllx) :: lagrange_x,lagrange_dx
  real(kind=kdble),dimension(ndim,2) :: lagrangeINF_x,lagrangeINF_dx
  !real(kind=kdble),dimension(nglly) :: lagrange_y,lagrange_dy
  !real(kind=kdble),dimension(ngllz) :: lagrange_z,lagrange_dz
  integer :: iINF !,ngllxINF(ndim)
  !print *,iface; stop
  !ngllxINF(1)=ngllx
  !ngllxINF(2)=nglly
  !ngllxINF(3)=ngllz

  ddir = one
  if (iface == 1) then
    iINF = 2; ddir=-one
  else if (iface == 2) then
    iINF = 1; ddir = one
  else if (iface == 3) then
    iINF = 2; ddir = one
  else if (iface == 4) then
    iINF = 1; ddir=-one
  else if (iface == 5) then
    iINF = 3; ddir=-one
  else if (iface == 6) then
    iINF = 3; ddir = one
  endif
  !print *,iface,iINF,ddir
  !ngllxINF(iINF)=ngllxINF(iINF)-1

  ! interpolation points
  ! compute everything in indexed order
  ! get GLL points
  ! for alpha=beta=0, jacobi polynomial is legendre polynomial
  ! for ngllx=nglly=ngllz, need to call only once
  call zwgljd(gllpx,gllwx,ngllx,jacobi_alpha,jacobi_beta)
  call zwgljd(gllpy,gllwy,nglly,jacobi_alpha,jacobi_beta)
  call zwgljd(gllpz,gllwz,ngllz,jacobi_alpha,jacobi_beta)


  ! integration points
  ! gauss-jacobi or gauss-legendre points and weights
  call zwgjd(igllpx,igllwx,nipx,jacobi_alpha,jacobi_beta)
  !print *,nipx
  !print *,gllpx
  !print *,gllwx; stop
  ! for an infinite element we use Gauss-Legendre quadrature
  !if (nip==8) then
  !  gllpx(1)=-one/sqrt(3.0_kdble); gllpx(2)=-gllpx(1)
  !  gllwx(1)=one; gllwx(2)=one;
  !else if (nip==27) then
  !  gllpx(1)=-sqrt(3.0_kdble/five); gllpx(2)=0.0_kdble; gllpx(3)=-gllpx(1)
  !  gllwx(1)=five/9.0_kdble; gllwx(2)=8.0_kdble/9.0_kdble; gllwx(3)=gllwx(1);
  !else
  if (nip /= 8.and.nip /= 27) then
    print *,'ERROR: illegal number of Gauss points:',nip,'!'
    stop
  endif
  !print *,igllpx
  !print *,igllwx
  !stop
  igllpy = igllpx; igllpz = igllpx;
  igllwy = igllwx; igllwz = igllwx;

  !print *,iface,iINF,ddir; stop
  !print *,gllpz
  !print *,gllwz
  !stop
  ii = 0
  do k1 = 1,nipx
    do j1 = 1,nipx
      do i1 = 1,nipx
        ii = ii+1
        !do ii=1,ngll ! ngllx*nglly*ngllz

        ! integration points
        xi(1)=igllpx(i1) !xi,   gll_points(1,ii)
        xi(2)=igllpy(j1) !eta,  gll_points(2,ii)
        xi(3)=igllpz(k1) !zeta, gll_points(3,ii)

        !xi(iINF)=ddir*xi(iINF)

        ! integration weights
        GLw(ii)=igllwx(i1)*igllwy(j1)*igllwz(k1)

        call lagrange1dGLL(ngllx,gllpx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
        call lagrange1dGLL(nglly,gllpy,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
        call lagrange1dGLL(ngllz,gllpz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))

        !call lagrange1dGEN(ngllx,xi(1),lagrange_x(1,:),lagrange_dx(1,:))
        !call lagrange1dGEN(nglly,xi(2),lagrange_x(2,:),lagrange_dx(2,:))
        !call lagrange1dGEN(ngllz,xi(3),lagrange_x(3,:),lagrange_dx(3,:))

        ! interpolation functions
        n = 0
        do k = 1,ngllz
          do j = 1,nglly
            do i = 1,ngllx
              n = n+1
              lagrange_gl(ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
              dlagrange_gl(1,ii,n)=lagrange_dx(1,i)*lagrange_x(2,j)*lagrange_x(3,k)
              dlagrange_gl(2,ii,n)=lagrange_x(1,i)*lagrange_dx(2,j)*lagrange_x(3,k)
              dlagrange_gl(3,ii,n)=lagrange_x(1,i)*lagrange_x(2,j)*lagrange_dx(3,k)
            enddo
          enddo
        enddo

        ! shape functions for HEX8
        ! compute 1d lagrange polynomials
        call lagrange1dGEN(2,xi(1),lagrangeINF_x(1,:),lagrangeINF_dx(1,:))
        call lagrange1dGEN(2,xi(2),lagrangeINF_x(2,:),lagrangeINF_dx(2,:))
        call lagrange1dGEN(2,xi(3),lagrangeINF_x(3,:),lagrangeINF_dx(3,:))
        ! consider 3 nodes but compute only at 2 nodes
        call lagrange1d_infiniteZW(3,xi(iINF),lagrangeINF_x(iINF,:),lagrangeINF_dx(iINF,:))
        !call lagrange_infinite(ngllx,nd,1.0_kdble,gam,a,lagrangeINF_x,lagrangeINF_dx)
        !print *,lagrangeINF_x,lagrangeINF_dx; stop
        !call lagrange1d_infinite(ngllx,xi,lagrangeINF_x,lagrangeINF_dx)
        !call lagrange1d_infiniteZW(ngllx,xi,lagrangeINF_x,lagrangeINF_dx)
        n = 0
        do k = 1,2
          do j = 1,2
            do i = 1,2
              n = n+1
              shape_infinite(ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(1,ii,n)=lagrangeINF_dx(1,i)*lagrangeINF_x(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(2,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_dx(2,j)*lagrangeINF_x(3,k)
              dshape_infinite(3,ii,n)=lagrangeINF_x(1,i)*lagrangeINF_x(2,j)*lagrangeINF_dx(3,k)
            enddo
          enddo
        enddo

        !enddo
      enddo
    enddo
  enddo

  return
  end subroutine shape_function_infiniteGLHEX8ZW_GQ

!
!===========================================
!
! this subroutine extracts the nodes for HEX8 of the finite region of an infinite element

  subroutine get_gnodinfHEX8(ndim,ngllx,nglly,ngllz,nginf,iface,gnodinf)

  implicit none
  integer,intent(in) :: ndim,ngllx,nglly,ngllz,nginf,iface
  integer,intent(out) :: gnodinf(nginf)
  integer :: i,j,k,inum
  integer :: inc(ndim),ngllxINF0(ndim),ngllxINF(ndim),iINF
  real(kind=kdble) :: ddir
  real(kind=kdble),parameter :: one=1.0_kdble

  if (iface < 1.or.iface > 6) then
    write(*,*) 'ERROR: illegal outer face ID:',iface
    stop
  endif

  ! initialize ngllINF indices
  ngllxINF0 = 1
  ngllxINF(1)=ngllx
  ngllxINF(2)=nglly
  ngllxINF(3)=ngllz

  if (iface == 1) then
    iINF = 2; ddir=-one
  else if (iface == 2) then
    iINF = 1; ddir = one
  else if (iface == 3) then
    iINF = 2; ddir = one
  else if (iface == 4) then
    iINF = 1; ddir=-one
  else if (iface == 5) then
    iINF = 3; ddir=-one
  else if (iface == 6) then
    iINF = 3; ddir = one
  endif

  if (ddir < 0) then
    ngllxINF0(iINF)=2
  else
    ngllxINF(iINF)=ngllxINF(iINF)-1
  endif

  ! extract only the corner nodes
  inc = ngllxINF-ngllxINF0
  inum = 0
  do k = ngllxINF0(3),ngllxINF(3),inc(3)
    do j = ngllxINF0(2),ngllxINF(2),inc(2)
      do i = ngllxINF0(1),ngllxINF(1),inc(1)
        inum = inum+1
        gnodinf(inum)=nglly*ngllx*(k-1)+ngllx*(j-1)+i
      enddo
    enddo
  enddo
  !print *,inum
  end subroutine get_gnodinfHEX8

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

  subroutine lagrange1d_infiniteMO(nenod,nd,xi,ddir,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  !integer :: i,j,k
  real(kind=kdble),intent(in) :: nd,xi,ddir ! xi: point where to calculate lagrange function and
  !its derivative
  !real(kind=kdble),intent(in) :: gam,a
  real(kind=kdble),dimension(nenod-1),intent(out) :: phi,dphi_dxi
  real(kind=kdble),dimension(nenod) :: xii
  real(kind=kdble) :: fac !dx
  real(kind=kdble),parameter :: one=1.0_kdble,two=2.0_kdble

  if (nenod /= 3) then
    write(*,*) 'ERROR: infinite element is currently implemented only for 3 nodes!'
    stop
  endif
  !! compute natural coordnates
  !dx=2.0_kdble/real((nenod-1),kdble)! length = 2.0 as xi is taken -1 to +1
  !do i=1,nenod
  !  ! coordinates when origin is in the left
  !  xii(i)=real((i-1),kdble)*dx
  !enddo

  !! origin is tranformed to mid point
  !xii=xii-1.0_kdble

  fac=one/(one-xi)

  phi(1)=-two*xi*fac
  phi(2)=one-phi(1)

  dphi_dxi(1)=-two*fac*fac
  dphi_dxi(2)=two*fac*fac

  return
  end subroutine lagrange1d_infiniteMO

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.
! Assumed Shape array: pass pointer, subarray or allocatable array

  subroutine lagrange1d_infiniteZWAS(nenod,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  !integer :: i,j,k
  real(kind=kdble),intent(in) :: xi ! xi: point where to calculate lagrange
  !function and its derivative
  real(kind=kdble),dimension(:),intent(out) :: phi,dphi_dxi !,dimension(nenod-1)
  real(kind=kdble) :: fac
  real(kind=kdble),parameter :: one=1.0_kdble

  if (nenod /= 3) then
    write(*,*) 'ERROR: infinite element is currently implemented only for 3 nodes!'
    stop
  endif

  fac=one/(one-xi)

  phi(1)=-xi*fac
  phi(2)=fac !one-phi(1) !+xi*fac

  dphi_dxi(1)=-fac*fac
  dphi_dxi(2)=fac*fac

  return
  end subroutine lagrange1d_infiniteZWAS

!
!===========================================
!

! this subroutine computes the 1d lagrange interpolation functions and their
! derivatives at a given point xi.

  subroutine lagrange1d_infiniteZW(nenod,xi,phi,dphi_dxi)

  implicit none
  integer,intent(in) :: nenod ! number of nodes in an 1d element
  !integer :: i,j,k
  real(kind=kdble),intent(in) :: xi ! xi: point where to calculate lagrange
  !function and its derivative

  real(kind=kdble),dimension(nenod-1),intent(out) :: phi,dphi_dxi
  real(kind=kdble) :: fac !dx
  real(kind=kdble),parameter :: one=1.0_kdble

  if (nenod /= 3) then
    write(*,*) 'ERROR: infinite element is currently implemented only for 3 nodes!'
    stop
  endif
  !! compute natural coordnates
  !dx=2.0_kdble/real((nenod-1),kdble)! length = 2.0 as xi is taken -1 to +1
  !do i=1,nenod
  !  ! coordinates when origin is in the left
  !  xii(i)=real((i-1),kdble)*dx
  !enddo

  !! origin is tranformed to mid point
  !xii=xii-1.0_kdble

  fac=one/(one-xi)

  phi(1)=-xi*fac
  phi(2)=fac !one-phi(1) !+xi*fac

  dphi_dxi(1)=-fac*fac
  dphi_dxi(2)=fac*fac

  return
  end subroutine lagrange1d_infiniteZW

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
  real(kind=kdble) :: rj,rn,twopi,fac,p(n,n+1),xold(n)
  real(kind=kdble),parameter :: one=1.0_kdble,pi=3.141592653589793_kdble, &
  two = 2.0_kdble,tol = 1.0e-12_kdble,zero = 0.0_kdble,zerotol = 1.0e-12_kdble

  if (n < 1) then
    write(*,*) 'ERROR: number of quadrature points must be > 1!'
    stop
  endif

  x = zero; w = zero
  rn=real(n,kdble)

  ! initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
  fac=two*pi/(two*rn-one)

  ! initialize the Legendre Vandermonde matrix.
  p = zero
  p(2:n,1) = one;
  do i = 1,n
    x(i)=-cos(fac*real(i-1,kdble))
    p(1,i) = (-one)**(i-1)
  enddo
  p(1,n+1)=(-one)**(n)

  ! compute P using the recursion relation.
  ! compute its first and second derivatives and
  ! update X using the Newton-Raphson method.
  xold = two
  do i = 1,100
    if (maxval(abs(x-xold)) <= zerotol)exit
    if (i >= 100) then
      write(*,*) 'ERROR: Legendre Vandermonde matrix does not converge!'
      stop
    endif
    xold = x;
    p(2:n,2) = x(2:n);
    do j = 2,n
      rj=real(j,kdble)
      p(2:n,j+1) = ((two*rj-one)*x(2:n)*p(2:n,j)+(-rj+one)*p(2:n,j-1))/rj
    enddo
    x(2:n) = xold(2:n)-((one-xold(2:n))/rn)*(p(2:n,n)+p(2:n,n+1))/(p(2:n,n)-p(2:n,n+1))
  enddo

  ! compute the weights.
  w = zero
  w(1) = two/(rn*rn)
  w(2:n)=(one-x(2:n))/(rn*p(2:n,n)*rn*p(2:n,n))
  return
  end subroutine radau_quadrature

end module infinite_element

#endif

