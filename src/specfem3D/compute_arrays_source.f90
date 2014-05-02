!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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

  subroutine compute_arrays_source(sourcearray, &
                                   xi_source,eta_source,gamma_source, &
                                   Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                                   xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                   xigll,yigll,zigll)

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  double precision :: xi_source,eta_source,gamma_source
  double precision :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: xix,xiy,xiz,etax,etay,etaz, &
        gammax,gammay,gammaz

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  ! local parameters
  double precision :: xixd,xiyd,xizd,etaxd,etayd,etazd,gammaxd,gammayd,gammazd
  ! source arrays
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: G11,G12,G13,G21,G22,G23,G31,G32,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer :: k,l,m

  ! calculate G_ij for general source location
  ! the source does not necessarily correspond to a Gauss-Lobatto point
  do m=1,NGLLZ
    do l=1,NGLLY
      do k=1,NGLLX

        if( CUSTOM_REAL == SIZE_REAL ) then
          xixd    = dble(xix(k,l,m))
          xiyd    = dble(xiy(k,l,m))
          xizd    = dble(xiz(k,l,m))
          etaxd   = dble(etax(k,l,m))
          etayd   = dble(etay(k,l,m))
          etazd   = dble(etaz(k,l,m))
          gammaxd = dble(gammax(k,l,m))
          gammayd = dble(gammay(k,l,m))
          gammazd = dble(gammaz(k,l,m))
        else
          xixd    = xix(k,l,m)
          xiyd    = xiy(k,l,m)
          xizd    = xiz(k,l,m)
          etaxd   = etax(k,l,m)
          etayd   = etay(k,l,m)
          etazd   = etaz(k,l,m)
          gammaxd = gammax(k,l,m)
          gammayd = gammay(k,l,m)
          gammazd = gammaz(k,l,m)
        endif

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

  subroutine compute_arrays_source_adjoint(myrank, adj_source_file, &
                                           xi_receiver,eta_receiver,gamma_receiver, nu,adj_sourcearray, &
                                           xigll,yigll,zigll,NSTEP_BLOCK,iadjsrc,it_sub_adj,NSTEP_SUB_ADJ, &
                                           NTSTEP_BETWEEN_READ_ADJSRC,DT)

  use constants,only: CUSTOM_REAL,SIZE_REAL,NDIM,NGLLX,NGLLY,NGLLZ,IIN_ADJ,R_EARTH
  use write_seismograms_mod, only: band_instrument_code

  implicit none

! input -- notice here NSTEP_BLOCK is different from the NSTEP in the main program
! instead NSTEP_BLOCK = iadjsrc_len(it_sub_adj), the length of this specific block

  integer myrank, NSTEP_BLOCK

  double precision xi_receiver, eta_receiver, gamma_receiver
  double precision DT

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
  double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd

  real(kind=CUSTOM_REAL) :: adj_src(NDIM,NSTEP_BLOCK)
  double precision, dimension(NDIM,NSTEP_BLOCK) :: adj_src_u

  integer icomp, itime, ios
  integer it_start,it_end,index_i
  real(kind=CUSTOM_REAL) :: junk
  character(len=3),dimension(NDIM) :: comp
  character(len=150) :: filename
  character(len=2) :: bic

! by Ebru
  call band_instrument_code(DT,bic)
  comp(1) = bic(1:2)//'N'
  comp(2) = bic(1:2)//'E'
  comp(3) = bic(1:2)//'Z'
!

  ! (sub)trace start and end
  ! reading starts in chunks of NSTEP_BLOCK from the end of the trace,
  ! i.e. as an example: total length NSTEP = 3000, chunk length NSTEP_BLOCK= 1000
  !                                then it will read in first it_start=2001 to it_end=3000,
  !                                second time, it will be it_start=1001 to it_end=2000 and so on...
  it_start = iadjsrc(it_sub_adj,1)
  it_end = iadjsrc(it_sub_adj,1)+NSTEP_BLOCK-1


  ! unfortunately, things become more tricky because of the Newmark time scheme at
  ! the very beginning of the time loop. however, when we read in the backward/reconstructed
  ! wavefields at the end of the first time loop, we can use the adjoint source index from 3000 down to 1.
  !
  ! see the comment on where we add the adjoint source (compute_add_sources_adjoint()).
  !
  ! otherwise,
  ! we would have to shift this indices by minus 1, to read in the adjoint source trace between 0 to 2999.
  ! since 0 index is out of bounds, we would have to put that adjoint source displacement artificially to zero
  !
  ! here now, it_start is now 2001 and it_end = 3000, then 1001 to 2000, then 1 to 1000.
  it_start = it_start
  it_end = it_end

  adj_src = 0._CUSTOM_REAL
  do icomp = 1, NDIM

    ! opens adjoint component file
    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
    open(unit=IIN_ADJ,file=trim(filename),status='old',action='read',iostat=ios)

    ! note: adjoint source files must be available for all three components E/N/Z, even
    !          if a component is just zeroed out
    if (ios /= 0) then
      ! adjoint source file not found
      ! stops simulation
      call exit_MPI(myrank,&
          'file '//trim(filename)//' not found, please check with your STATIONS_ADJOINT file')
    endif
    !if (ios /= 0) cycle ! cycles to next file - this is too error prone and users might easily end up with wrong results

    ! jumps over unused trace length
    do itime =1,it_start-1
      read(IIN_ADJ,*,iostat=ios) junk,junk
      if( ios /= 0) &
        call exit_MPI(myrank,&
          'file '//trim(filename)//' has wrong length, please check with your simulation duration')
    enddo

    ! reads in (sub)trace
    do itime = it_start,it_end

      ! index will run from 1 to NSTEP_BLOCK
      index_i = itime - it_start + 1

      ! would skip read and set source artificially to zero if out of bounds, see comments above
      if( it_start == 0 .and. itime == 0 ) then
        adj_src(icomp,1) = 0._CUSTOM_REAL
        cycle
      endif

      ! reads in adjoint source trace
      !read(IIN_ADJ,*,iostat=ios) junk, adj_src(icomp,itime-it_start+1)
      read(IIN_ADJ,*,iostat=ios) junk, adj_src(icomp,index_i)

      if( ios /= 0) &
        call exit_MPI(myrank, &
          'file '//trim(filename)//' has wrong length, please check with your simulation duration')
    enddo

    close(IIN_ADJ)

  enddo

  ! non-dimensionalize
  adj_src = adj_src*scale_displ_inv

  ! rotates to Cartesian
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
  do itime = 1, NSTEP_BLOCK

    ! multiply with interpolators
    call multiply_arrays_adjoint(sourcearrayd,hxir,hetar,hgammar,adj_src_u(:,itime))

    ! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      adj_sourcearray(:,:,:,:,itime) = sngl(sourcearrayd(:,:,:,:))
    else
      adj_sourcearray(:,:,:,:,itime) = sourcearrayd(:,:,:,:)
    endif

  enddo

  end subroutine compute_arrays_source_adjoint
