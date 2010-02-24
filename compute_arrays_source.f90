!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine compute_arrays_source(ispec_selected_source, &
             xi_source,eta_source,gamma_source,sourcearray, &
             Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
             xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
             xigll,yigll,zigll,nspec)

  implicit none

  include "constants.h"

  integer ispec_selected_source,nspec

  double precision xi_source,eta_source,gamma_source
  double precision Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,l,m

! calculate G_ij for general source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX

        xixd    = dble(xix(k,l,m,ispec_selected_source))
        xiyd    = dble(xiy(k,l,m,ispec_selected_source))
        xizd    = dble(xiz(k,l,m,ispec_selected_source))
        etaxd   = dble(etax(k,l,m,ispec_selected_source))
        etayd   = dble(etay(k,l,m,ispec_selected_source))
        etazd   = dble(etaz(k,l,m,ispec_selected_source))
        gammaxd = dble(gammax(k,l,m,ispec_selected_source))
        gammayd = dble(gammay(k,l,m,ispec_selected_source))
        gammazd = dble(gammaz(k,l,m,ispec_selected_source))

        G11(k,l,m) = Mxx*xixd+Mxy*xiyd+Mxz*xizd
        G12(k,l,m) = Mxx*etaxd+Mxy*etayd+Mxz*etazd
        G13(k,l,m) = Mxx*gammaxd+Mxy*gammayd+Mxz*gammazd
        G21(k,l,m) = Mxy*xixd+Myy*xiyd+Myz*xizd
        G22(k,l,m) = Mxy*etaxd+Myy*etayd+Myz*etazd
        G23(k,l,m) = Mxy*gammaxd+Myy*gammayd+Myz*gammazd
        G31(k,l,m) = Mxz*xixd+Myz*xiyd+Mzz*xizd
        G32(k,l,m) = Mxz*etaxd+Myz*etayd+Mzz*etazd
        G33(k,l,m) = Mxz*gammaxd+Myz*gammayd+Mzz*gammazd

      enddo
    enddo
  enddo

! compute Lagrange polynomials at the source location
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

! calculate source array
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX
        call multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)
      enddo
    enddo
  enddo

! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    sourcearray(:,:,:,:) = sngl(sourcearrayd(:,:,:,:))
  else
    sourcearray(:,:,:,:) = sourcearrayd(:,:,:,:)
  endif

  end subroutine compute_arrays_source

!================================================================

! we put these multiplications in a separate routine because otherwise
! some compilers try to unroll the six loops above and take forever to compile
  subroutine multiply_arrays_source(sourcearrayd,G11,G12,G13,G21,G22,G23, &
                  G31,G32,G33,hxis,hpxis,hetas,hpetas,hgammas,hpgammas,k,l,m)

  implicit none

  include "constants.h"

! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,l,m

  integer ir,it,iv

  sourcearrayd(:,k,l,m) = ZERO

  do iv=1,NGLLZ
    do it=1,NGLLY
      do ir=1,NGLLX

        sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                           *(G11(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                           +G12(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                           +G13(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

        sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                           *(G21(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                           +G22(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                           +G23(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

        sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + hxis(ir)*hetas(it)*hgammas(iv) &
                           *(G31(ir,it,iv)*hpxis(k)*hetas(l)*hgammas(m) &
                           +G32(ir,it,iv)*hxis(k)*hpetas(l)*hgammas(m) &
                           +G33(ir,it,iv)*hxis(k)*hetas(l)*hpgammas(m))

      enddo
    enddo
  enddo

  end subroutine multiply_arrays_source

!================================================================

subroutine compute_arrays_source_adjoint(myrank, adj_source_file, &
      xi_receiver,eta_receiver,gamma_receiver, nu,adj_sourcearray, &
      xigll,yigll,zigll,NSTEP_BLOCK,iadjsrc,it_sub_adj,NSTEP_SUB_ADJ, &
      NTSTEP_BETWEEN_READ_ADJSRC)

  implicit none

  include 'constants.h'

! input -- notice here NSTEP_BLOCK is different from the NSTEP in the main program
! instead NSTEP_BLOCK = iadjsrc_len(it_sub_adj), the length of this specific block

  integer myrank, NSTEP_BLOCK

  double precision xi_receiver, eta_receiver, gamma_receiver

  character(len=*) adj_source_file

  ! Vala added
  integer it_sub_adj,NSTEP_SUB_ADJ,NTSTEP_BETWEEN_READ_ADJSRC
  integer, dimension(NSTEP_SUB_ADJ,2) :: iadjsrc

  ! output
  real(kind=CUSTOM_REAL) :: adj_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NTSTEP_BETWEEN_READ_ADJSRC)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  double precision, dimension(NDIM,NDIM) :: nu

  double precision,parameter :: scale_displ_inv = 1.d0/R_EARTH

  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
  real(kind=CUSTOM_REAL) :: adj_src(NDIM,NSTEP_BLOCK),adj_src_u(NDIM,NSTEP_BLOCK)

  integer icomp, itime, i, j, k, ios
  integer it_start,it_end,index_i
  real(kind=CUSTOM_REAL) :: junk
  character(len=3),dimension(NDIM) :: comp = (/ "LHN", "LHE", "LHZ" /)
  character(len=150) :: filename
  
  ! (sub)trace start and end
  ! reading starts in chunks of NSTEP_BLOCK from the end of the trace,
  ! i.e. as an example: total length NSTEP = 3000, chunk length NSTEP_BLOCK= 1000  
  !                                then it will read in first it_start=2001 to it_end=3000, 
  !                                second time, it will be it_start=1001 to it_end=2000 and so on...
  it_start = iadjsrc(it_sub_adj,1)
  it_end = iadjsrc(it_sub_adj,1)+NSTEP_BLOCK-1

  
  ! unfortunately, things become more tricky because of the Newark time scheme at
  ! the very beginning of the time loop.
  !
  ! see the comment on where we add the adjoint source (compute_add_sources_adjoint()).
  !
  ! we will thus shift this indices by minus 1, to read in the adjoint source trace between 0 to 2999.
  ! since 0 index is out of bounds, we put that adjoint source displacement artifically to zero
  !
  ! that is e.g., it_start is now 2000 and it_end = 2999, then 1000 to 1999, then 0 to 999.
  it_start = it_start - 1
  it_end = it_end - 1
  
  adj_src = 0._CUSTOM_REAL
  do icomp = 1, NDIM

    ! opens adjoint component file
    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)
    if (ios /= 0) cycle ! cycles to next file
     
    ! jumps over unused trace length
    do itime =1,it_start-1
      read(IIN,*,iostat=ios) junk,junk
      if( ios /= 0) &
        call exit_MPI(myrank,&
          'file '//trim(filename)//' has wrong length, please check with your simulation duration')
    enddo
     
    ! reads in (sub)trace
    do itime = it_start,it_end

      ! index will run from 1 to NSTEP_BLOCK
      index_i = itime - it_start + 1
    
      ! skips read and sets source artifically to zero if out of bounds, see comments above
      if( it_start == 0 .and. itime == 0 ) then 
        adj_src(icomp,1) = 0._CUSTOM_REAL
        cycle
      endif
      
      ! reads in adjoint source trace
      !read(IIN,*,iostat=ios) junk, adj_src(icomp,itime-it_start+1)
      read(IIN,*,iostat=ios) junk, adj_src(icomp,index_i)
      
      if( ios /= 0) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration')
    enddo
    
    close(IIN)

  enddo

  ! non-dimensionalize 
  adj_src = adj_src*scale_displ_inv

  ! rotates to cartesian
  do itime = 1, NSTEP_BLOCK
    adj_src_u(:,itime) = nu(1,:) * adj_src(1,itime) &
                       + nu(2,:) * adj_src(2,itime) &
                       + nu(3,:) * adj_src(3,itime)
  enddo

  ! receiver interpolators  
  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_receiver,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)

  ! adds interpolated source contribution to all GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        do itime = 1, NSTEP_BLOCK
          adj_sourcearray(:,i,j,k,itime) = hxir(i) * hetar(j) * hgammar(k) * adj_src_u(:,itime)
        enddo
      enddo
    enddo
  enddo


end subroutine compute_arrays_source_adjoint

! =======================================================================

! compute the integrated derivatives of source parameters (M_jk and X_s)

subroutine compute_adj_source_frechet(displ_s,Mxx,Myy,Mzz,Mxy,Mxz,Myz,eps_s,eps_m_s, &
           hxir,hetar,hgammar,hpxir,hpetar,hpgammar, hprime_xx,hprime_yy,hprime_zz, &
           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  implicit none

  include 'constants.h'

  ! input
  real(kind=CUSTOM_REAL) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ)
  double precision :: Mxx, Myy, Mzz, Mxy, Mxz, Myz
  ! output
  real(kind=CUSTOM_REAL) :: eps_s(NDIM,NDIM), eps_m_s(NDIM)

  ! auxilliary
  double precision :: hxir(NGLLX), hetar(NGLLY), hgammar(NGLLZ), &
             hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

! local variables
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l, tempy1l,tempy2l,tempy3l, &
             tempz1l,tempz2l,tempz3l, hp1, hp2, hp3, &
             xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, &
             duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl, &
             xix_s,xiy_s,xiz_s,etax_s,etay_s,etaz_s,gammax_s,gammay_s,gammaz_s, &
             hlagrange_xi, hlagrange_eta, hlagrange_gamma, hlagrange

  real(kind=CUSTOM_REAL) :: eps(NDIM,NDIM), eps_array(NDIM,NDIM,NGLLX,NGLLY,NGLLZ), &
             eps_m_array(NGLLX,NGLLY,NGLLZ)

  integer i,j,k,l


! first compute the strain at all the GLL points of the source element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX

        tempx1l = 0._CUSTOM_REAL
        tempx2l = 0._CUSTOM_REAL
        tempx3l = 0._CUSTOM_REAL

        tempy1l = 0._CUSTOM_REAL
        tempy2l = 0._CUSTOM_REAL
        tempy3l = 0._CUSTOM_REAL

        tempz1l = 0._CUSTOM_REAL
        tempz2l = 0._CUSTOM_REAL
        tempz3l = 0._CUSTOM_REAL

        do l=1,NGLLX
          hp1 = hprime_xx(i,l)
          tempx1l = tempx1l + displ_s(1,l,j,k)*hp1
          tempy1l = tempy1l + displ_s(2,l,j,k)*hp1
          tempz1l = tempz1l + displ_s(3,l,j,k)*hp1

          hp2 = hprime_yy(j,l)
          tempx2l = tempx2l + displ_s(1,i,l,k)*hp2
          tempy2l = tempy2l + displ_s(2,i,l,k)*hp2
          tempz2l = tempz2l + displ_s(3,i,l,k)*hp2

          hp3 = hprime_zz(k,l)
          tempx3l = tempx3l + displ_s(1,i,j,l)*hp3
          tempy3l = tempy3l + displ_s(2,i,j,l)*hp3
          tempz3l = tempz3l + displ_s(3,i,j,l)*hp3
        enddo

! dudx
        xixl = xix(i,j,k)
        xiyl = xiy(i,j,k)
        xizl = xiz(i,j,k)
        etaxl = etax(i,j,k)
        etayl = etay(i,j,k)
        etazl = etaz(i,j,k)
        gammaxl = gammax(i,j,k)
        gammayl = gammay(i,j,k)
        gammazl = gammaz(i,j,k)

        duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
        duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
        duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

        duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
        duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
        duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

        duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
        duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
        duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! strain eps_jk
        eps(1,1) = duxdxl
        eps(1,2) = (duxdyl + duydxl) / 2
        eps(1,3) = (duxdzl + duzdxl) / 2
        eps(2,2) = duydyl
        eps(2,3) = (duydzl + duzdyl) / 2
        eps(3,3) = duzdzl
        eps(2,1) = eps(1,2)
        eps(3,1) = eps(1,3)
        eps(3,2) = eps(2,3)

        eps_array(:,:,i,j,k) = eps(:,:)

! Mjk eps_jk
        eps_m_array(i,j,k) = Mxx * eps(1,1) + Myy * eps(2,2) + Mzz * eps(3,3) + &
                   2 * (Mxy * eps(1,2) + Mxz * eps(1,3) + Myz * eps(2,3))

      enddo
    enddo
  enddo

  ! interpolate the strain eps_s(:,:) from eps_array(:,:,i,j,k)
  eps_s = 0.
  xix_s = 0.;  xiy_s = 0.;  xiz_s = 0.
  etax_s = 0.; etay_s = 0.; etaz_s = 0.
  gammax_s = 0.; gammay_s = 0.; gammaz_s = 0.

  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        hlagrange = hxir(i)*hetar(j)*hgammar(k)

        eps_s(1,1) = eps_s(1,1) + eps_array(1,1,i,j,k)*hlagrange
        eps_s(1,2) = eps_s(1,2) + eps_array(1,2,i,j,k)*hlagrange
        eps_s(1,3) = eps_s(1,3) + eps_array(1,3,i,j,k)*hlagrange
        eps_s(2,2) = eps_s(2,2) + eps_array(2,2,i,j,k)*hlagrange
        eps_s(2,3) = eps_s(2,3) + eps_array(2,3,i,j,k)*hlagrange
        eps_s(3,3) = eps_s(3,3) + eps_array(3,3,i,j,k)*hlagrange

        xix_s = xix_s + xix(i,j,k)*hlagrange
        xiy_s = xiy_s + xiy(i,j,k)*hlagrange
        xiz_s = xiz_s + xiz(i,j,k)*hlagrange
        etax_s = etax_s + etax(i,j,k)*hlagrange
        etay_s = etay_s + etay(i,j,k)*hlagrange
        etaz_s = etaz_s + etaz(i,j,k)*hlagrange
        gammax_s = gammax_s + gammax(i,j,k)*hlagrange
        gammay_s = gammay_s + gammay(i,j,k)*hlagrange
        gammaz_s = gammaz_s + gammaz(i,j,k)*hlagrange

      enddo
    enddo
  enddo

! for completion purpose, not used in specfem3D.f90
  eps_s(2,1) = eps_s(1,2)
  eps_s(3,1) = eps_s(1,3)
  eps_s(3,2) = eps_s(2,3)

! compute the gradient of M_jk * eps_jk, and then interpolate it

  eps_m_s = 0.
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        hlagrange_xi = hpxir(i)*hetar(j)*hgammar(k)
        hlagrange_eta = hxir(i)*hpetar(j)*hgammar(k)
        hlagrange_gamma = hxir(i)*hetar(j)*hpgammar(k)

        eps_m_s(1) = eps_m_s(1) +  eps_m_array(i,j,k) * (hlagrange_xi * xix_s &
                   + hlagrange_eta * etax_s + hlagrange_gamma * gammax_s)
        eps_m_s(2) = eps_m_s(2) +  eps_m_array(i,j,k) * (hlagrange_xi * xiy_s &
                   + hlagrange_eta * etay_s + hlagrange_gamma * gammay_s)
        eps_m_s(3) = eps_m_s(3) +  eps_m_array(i,j,k) * (hlagrange_xi * xiz_s &
                   + hlagrange_eta * etaz_s + hlagrange_gamma * gammaz_s)

      enddo
    enddo
  enddo

end subroutine compute_adj_source_frechet

!================================================================
!
! deprecated...
!
!subroutine compute_arrays_adjoint_source(myrank, adj_source_file, &
!      xi_receiver,eta_receiver,gamma_receiver, nu,adj_sourcearray, &
!      xigll,yigll,zigll,NSTEP)
!
!  implicit none
!
!  include 'constants.h'
!
!! input
!  integer myrank, NSTEP
!
!  double precision xi_receiver, eta_receiver, gamma_receiver
!
!  character(len=*) adj_source_file
!
!! output
!  real(kind=CUSTOM_REAL) :: adj_sourcearray(NSTEP,NDIM,NGLLX,NGLLY,NGLLZ)
!
!! Gauss-Lobatto-Legendre points of integration and weights
!  double precision, dimension(NGLLX) :: xigll
!  double precision, dimension(NGLLY) :: yigll
!  double precision, dimension(NGLLZ) :: zigll
!
!  double precision, dimension(NDIM,NDIM) :: nu
!
!  double precision scale_displ
!
!  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
!        hgammar(NGLLZ), hpgammar(NGLLZ)
!  real(kind=CUSTOM_REAL) :: adj_src(NSTEP,NDIM),adj_src_u(NSTEP,NDIM)
!
!  integer icomp, itime, i, j, k, ios
!  double precision :: junk
!  character(len=3) :: comp(NDIM)
!  character(len=150) :: filename
!
!  scale_displ = R_EARTH
!
!  call lagrange_any(xi_receiver,NGLLX,xigll,hxir,hpxir)
!  call lagrange_any(eta_receiver,NGLLY,yigll,hetar,hpetar)
!  call lagrange_any(gamma_receiver,NGLLZ,zigll,hgammar,hpgammar)
!
!  adj_sourcearray(:,:,:,:,:) = 0.
!
!  comp = (/"LHN", "LHE", "LHZ"/)
!
!  do icomp = 1, NDIM
!
!    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
!    open(unit = IIN, file = trim(filename), iostat = ios)
!    if (ios /= 0) call exit_MPI(myrank, ' file '//trim(filename)//' does not exist')
!    do itime = 1, NSTEP
!      read(IIN,*) junk, adj_src(itime,icomp)
!    enddo
!    close(IIN)
!
!  enddo
!
!  adj_src = adj_src/scale_displ
!
!  do itime = 1, NSTEP
!    adj_src_u(itime,:) = nu(1,:) * adj_src(itime,1) + nu(2,:) * adj_src(itime,2) + nu(3,:) * adj_src(itime,3)
!  enddo
!
!  do k = 1, NGLLZ
!    do j = 1, NGLLY
!      do i = 1, NGLLX
!        adj_sourcearray(:,:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src_u(:,:)
!      enddo
!    enddo
!  enddo
!
!
!end subroutine compute_arrays_adjoint_source
!
