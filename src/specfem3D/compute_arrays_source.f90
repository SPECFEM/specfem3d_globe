!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
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
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLY) :: hetas,hpetas
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  double precision :: hlagrange
  double precision :: dsrc_dx, dsrc_dy, dsrc_dz
  double precision :: dxis_dx, detas_dx, dgammas_dx
  double precision :: dxis_dy, detas_dy, dgammas_dy
  double precision :: dxis_dz, detas_dz, dgammas_dz

  integer :: k,l,m

! compute Lagrange polynomials at the source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(eta_source,NGLLY,yigll,hetas,hpetas)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

  dxis_dx = ZERO
  dxis_dy = ZERO
  dxis_dz = ZERO
  detas_dx = ZERO
  detas_dy = ZERO
  detas_dz = ZERO
  dgammas_dx = ZERO
  dgammas_dy = ZERO
  dgammas_dz = ZERO

  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           xixd    = dble(xix(k,l,m))
           xiyd    = dble(xiy(k,l,m))
           xizd    = dble(xiz(k,l,m))
           etaxd   = dble(etax(k,l,m))
           etayd   = dble(etay(k,l,m))
           etazd   = dble(etaz(k,l,m))
           gammaxd = dble(gammax(k,l,m))
           gammayd = dble(gammay(k,l,m))
           gammazd = dble(gammaz(k,l,m))

           hlagrange = hxis(k) * hetas(l) * hgammas(m)

           dxis_dx = dxis_dx + hlagrange * xixd
           dxis_dy = dxis_dy + hlagrange * xiyd
           dxis_dz = dxis_dz + hlagrange * xizd

           detas_dx = detas_dx + hlagrange * etaxd
           detas_dy = detas_dy + hlagrange * etayd
           detas_dz = detas_dz + hlagrange * etazd

           dgammas_dx = dgammas_dx + hlagrange * gammaxd
           dgammas_dy = dgammas_dy + hlagrange * gammayd
           dgammas_dz = dgammas_dz + hlagrange * gammazd

       enddo
     enddo
  enddo

! calculate source array
  sourcearrayd(:,:,:,:) = ZERO
  do m = 1,NGLLZ
     do l = 1,NGLLY
        do k = 1,NGLLX

           dsrc_dx = (hpxis(k)*dxis_dx)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dx)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dx)
           dsrc_dy = (hpxis(k)*dxis_dy)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dy)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dy)
           dsrc_dz = (hpxis(k)*dxis_dz)*hetas(l)*hgammas(m) + hxis(k)*(hpetas(l)*detas_dz)*hgammas(m) + &
                                                                hxis(k)*hetas(l)*(hpgammas(m)*dgammas_dz)

           sourcearrayd(1,k,l,m) = sourcearrayd(1,k,l,m) + (Mxx*dsrc_dx + Mxy*dsrc_dy + Mxz*dsrc_dz)
           sourcearrayd(2,k,l,m) = sourcearrayd(2,k,l,m) + (Mxy*dsrc_dx + Myy*dsrc_dy + Myz*dsrc_dz)
           sourcearrayd(3,k,l,m) = sourcearrayd(3,k,l,m) + (Mxz*dsrc_dx + Myz*dsrc_dy + Mzz*dsrc_dz)

       enddo
     enddo
  enddo

  ! distinguish between single and double precision for reals
  sourcearray(:,:,:,:) = real(sourcearrayd(:,:,:,:), kind=CUSTOM_REAL)

  end subroutine compute_arrays_source

!================================================================

  subroutine compute_arrays_source_adjoint(myrank, adj_source_file, &
                                           nu,source_adjoint, &
                                           NSTEP_BLOCK,iadjsrc,it_sub_adj,NSTEP_SUB_ADJ, &
                                           NTSTEP_BETWEEN_READ_ADJSRC,DT)

  use constants, only: CUSTOM_REAL,SIZE_REAL,NDIM,NGLLX,NGLLY,NGLLZ,IIN_ADJ,R_EARTH,MAX_STRING_LEN

  use specfem_par, only: scale_displ_inv, NUMBER_OF_SIMULTANEOUS_RUNS, READ_ADJSRC_ASDF, mygroup

  use iso_c_binding, only: C_NULL_CHAR

  implicit none

! input -- notice here NSTEP_BLOCK is different from the NSTEP in the main program
! instead NSTEP_BLOCK = iadjsrc_len(it_sub_adj), the length of this specific block

  integer,intent(in) :: myrank
  character(len=*),intent(in) :: adj_source_file

  double precision, dimension(NDIM,NDIM),intent(in) :: nu


  ! output
  integer,intent(in) :: NTSTEP_BETWEEN_READ_ADJSRC
  real(kind=CUSTOM_REAL),intent(out) :: source_adjoint(NDIM,NTSTEP_BETWEEN_READ_ADJSRC)

  integer,intent(in) :: NSTEP_SUB_ADJ
  integer, dimension(NSTEP_SUB_ADJ,2),intent(in) :: iadjsrc

  integer,intent(in) :: NSTEP_BLOCK
  integer,intent(in) :: it_sub_adj

  double precision,intent(in) :: DT

  ! local parameters
  double precision, dimension(NDIM,NSTEP_BLOCK) :: adj_src_u

  real(kind=CUSTOM_REAL), dimension(NDIM,NSTEP_BLOCK) :: adj_src
  real(kind=CUSTOM_REAL), dimension(NSTEP_BLOCK) :: adj_source_asdf
  real(kind=CUSTOM_REAL) :: junk

  integer :: icomp, itime, ios
  integer :: index_start,index_end,index_i
  character(len=3),dimension(NDIM) :: comp
  character(len=MAX_STRING_LEN) :: filename, path_to_add
  character(len=80) :: adj_source_name
  character(len=2) :: bic

  call band_instrument_code(DT,bic)
  comp(1) = bic(1:2)//'N'
  comp(2) = bic(1:2)//'E'
  comp(3) = bic(1:2)//'Z'

  ! (sub)trace start and end
  ! reading starts in chunks of NSTEP_BLOCK from the end of the trace,
  ! i.e. as an example: total length NSTEP = 3000, chunk length NSTEP_BLOCK= 1000
  !                                then it will read in first index_start=2001 to index_end=3000,
  !                                second time, it will be index_start=1001 to index_end=2000 and so on...
  index_start = iadjsrc(it_sub_adj,1)
  index_end = iadjsrc(it_sub_adj,1)+NSTEP_BLOCK-1

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
  ! here now, index_start is now 2001 and index_end = 3000, then 1001 to 2000, then 1 to 1000.
  index_start = index_start
  index_end = index_end

  itime = 0
  adj_src(:,:) = 0._CUSTOM_REAL

  if (READ_ADJSRC_ASDF) then

    do icomp = 1, NDIM ! 3 components

      ! print *, "READING ADJOINT SOURCES USING ASDF"

      adj_source_name = trim(adj_source_file) // '_' // comp(icomp)

      ! would skip read and set source artificially to zero if out of bounds, see comments above
      if (index_start == 0 .and. itime == 0) then
        adj_src(icomp,1) = 0._CUSTOM_REAL
        cycle
      endif

      call read_adjoint_sources_ASDF(adj_source_name, adj_source_asdf, index_start, index_end)

      adj_src(icomp,:) = real(adj_source_asdf(1:NSTEP_BLOCK))

    enddo

  else

    do icomp = 1, NDIM

      ! opens adjoint component file
      filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'

      if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
        write(path_to_add,"('run',i4.4,'/')") mygroup + 1
        filename = path_to_add(1:len_trim(path_to_add))//filename(1:len_trim(filename))
      endif

      open(unit=IIN_ADJ,file=trim(filename),status='old',action='read',iostat=ios)

      ! note: adjoint source files must be available for all three components E/N/Z, even
      !          if a component is just zeroed out
      if (ios /= 0) then
        ! adjoint source file not found
        ! stops simulation
        call exit_MPI(myrank, &
            'file '//trim(filename)//' not found, please check with your STATIONS_ADJOINT file')
      endif
      !if (ios /= 0) cycle ! cycles to next file - this is too error prone and users might easily end up with wrong results

      ! jumps over unused trace length
      do itime  = 1,index_start-1
        read(IIN_ADJ,*,iostat=ios) junk,junk
        if (ios /= 0) &
          call exit_MPI(myrank, &
            'file '//trim(filename)//' has wrong length, please check with your simulation duration')
      enddo

      ! reads in (sub)trace
      do itime = index_start,index_end

        ! index will run from 1 to NSTEP_BLOCK
        index_i = itime - index_start + 1

        ! would skip read and set source artificially to zero if out of bounds, see comments above
        if (index_start == 0 .and. itime == 0) then
          adj_src(icomp,1) = 0._CUSTOM_REAL
          cycle
        endif

        ! reads in adjoint source trace
        !read(IIN_ADJ,*,iostat=ios) junk, adj_src(icomp,itime-index_start+1)
        read(IIN_ADJ,*,iostat=ios) junk, adj_src(icomp,index_i)

        if (ios /= 0) then
          print *,'Error reading adjoint source: ',trim(filename)
          print *,'rank ',myrank,' - time step: ',itime,' index_start: ',index_start,' index_end: ',index_end
          print *,'  ',trim(filename)//'has wrong length, please check with your simulation duration'
          call exit_MPI(myrank,'file '//trim(filename)//' has wrong length, please check with your simulation duration')
        endif
      enddo

      close(IIN_ADJ)

    enddo
  endif

  ! non-dimensionalize
  adj_src(:,:) = adj_src(:,:) * scale_displ_inv

  ! rotates to Cartesian
  do itime = 1, NSTEP_BLOCK
    adj_src_u(:,itime) = nu(1,:) * adj_src(1,itime) &
                       + nu(2,:) * adj_src(2,itime) &
                       + nu(3,:) * adj_src(3,itime)
  enddo

  do icomp = 1, NDIM
    source_adjoint(icomp,:) = adj_src_u(icomp,:)
  enddo

  contains

    subroutine multiply_arrays_adjoint(sourcearrayd,hxir,hetar,hgammar,adj_src_ud)

    use constants

    implicit none

    double precision, dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearrayd
    double precision, dimension(NGLLX) :: hxir
    double precision, dimension(NGLLY) :: hetar
    double precision, dimension(NGLLZ) :: hgammar
    double precision, dimension(NDIM) :: adj_src_ud

    integer :: i,j,k

    ! adds interpolated source contribution to all GLL points within this element
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          sourcearrayd(:,i,j,k) = hxir(i) * hetar(j) * hgammar(k) * adj_src_ud(:)
        enddo
      enddo
    enddo

    end subroutine multiply_arrays_adjoint

  end subroutine compute_arrays_source_adjoint
