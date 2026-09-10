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

!----
!---- A runnable tour of the gf3d Fortran API, on the global example.
!----
!---- The Python counterpart of this program is
!---- ../python/gf3d_api_demo.py; it prints the same numbers, because it
!---- drives the same library. This one is what a Fortran caller sees:
!----
!----     use gf3d
!----
!---- and nothing else. Everything below -- the types, the routines, the
!---- error codes, the names and units of the partial derivatives -- comes
!---- from that one module. If it did not, this program would not compile,
!---- which is rather the point of writing it.
!----
!---- Build (see the Makefile):
!----
!----     gfortran gf3d_api_demo.f90 -I<repo>/include \
!----              -L<repo>/lib -lgf3d -Wl,-rpath,<repo>/lib
!----
!---- No HDF5 flags: lib/libgf3d.so records HDF5 as a dependency of its own
!---- and carries the library path with it, so a caller does not have to
!---- know where this tree's HDF5 came from. Linking the static
!---- lib/libgf3d.a instead works too, but then the full HDF5 link line is
!---- the caller's problem.
!----
!---- Five sections, each printing what it measured:
!----
!----   1. open the database and say what is in it
!----   2. read a CMTSOLUTION and locate it in the mesh
!----   3. plan the output axis and the source time function conversion
!----   4. extract seismograms and their ten partial derivatives
!----   5. use a centroid partial as a derivative, and show it converges
!----
!---- Section 5 is the one worth reading. Sections 1 to 4 can be had from
!---- bin/xgf3d and a pile of SAC files; what the library is for is section
!---- 5, where the partials are Frechet derivatives to be used rather than
!---- traces to be written.
!----

  program gf3d_api_demo

  use gf3d

  implicit none

  ! paths, relative to this directory
  character(len=*), parameter :: DB_PATH = '../../GFDB'
  character(len=*), parameter :: CMT_PATH = '../../validation_data/CMTSOLUTION'

  ! the relocation used in section 5, in degrees of latitude, and how many
  ! times it is halved
  double precision, parameter :: STEP0 = 0.05d0
  integer, parameter :: NHALVE = 2

  type(t_gfdb) :: db
  type(t_gf_source) :: src,moved
  type(t_gf_location) :: loc
  type(t_gf_taxis) :: tax
  type(t_gf_stf) :: stf

  double precision, dimension(:,:,:), allocatable :: synt,truth
  double precision, dimension(:,:,:,:), allocatable :: dp
  double precision, dimension(:), allocatable :: t

  ! scratch outputs for the itypsokern = 0 calls of section 5: they are
  ! allocated zero-sized and never read, but get_seismograms allocates
  ! through them, so they have to be variables rather than expressions
  double precision, dimension(:,:,:,:), allocatable :: dp_unused
  double precision, dimension(:), allocatable :: t_unused

  double precision :: t_open,t_locate,t_locate2,t_extract,t_total
  double precision, dimension(0:NHALVE) :: t_reloc,step,err
  double precision :: peak,worst,ratio
  integer :: ierr,i,ista,icomp,ip,k,nsta,nt

  character(len=34) :: label

  ! for tic()/toc(), by host association. Both are integer(8) so that
  ! system_clock reports nanoseconds rather than milliseconds -- some of
  ! the calls below take well under one of the latter.
  integer(kind=8) :: clock_start,clock_rate

  logical :: ok

  ok = .true.
  t_total = 0.d0

  write(*,*)
  write(*,*) '========================================================'
  write(*,*) ' gf3d Fortran API demonstration'
  write(*,*) '========================================================'
  write(*,'(a,a)') '  library version : ',trim(GF3D_VERSION)

!
!--- 1. the database ---------------------------------------------------
!

  call banner('1. the database')

  ! check_completion = .false. skips verifying that every element/station
  ! file the index promises is actually present. That check is one stat()
  ! per file -- 328 of them here, about half a second -- and is what
  ! `xgf3d --info` does; an extraction does not need it.
  call tic()
  call gf_open(DB_PATH,db,ierr,check_completion = .false.)
  call toc(t_open)

  if (ierr /= GF_OK) then
    write(*,*)
    write(*,'(a,a)') '  could not open ',DB_PATH
    write(*,'(a,a)') '  ',trim(gf_errmsg)
    write(*,*)
    write(*,*) '  The example database is built by the workflow, and is'
    write(*,*) '  gitignored. From EXAMPLES/green_function_database/global:'
    write(*,*) '      snakemake -j1'
    stop 1
  endif

  write(*,'(a,f9.1,a)') '  gf_open took ',t_open*1.d3,' ms (metadata only;' // &
                        ' elements are read on demand)'
  write(*,*)

  ! straight out of the handle: this is what a caller reads
  write(*,'(a,a)')                '  path             : ',trim(db%path)
  write(*,'(a,i0)')               '  elements         : ',db%nelem
  write(*,'(a,i0,a,i0,a)')        '  stored samples   : ',db%nt_subsampled, &
                                  ' of ',db%nstep,' solver steps'
  write(*,'(a,f6.3,a,i0,a,f7.3,a)') '  sample spacing   : ',db%dt,' s x ', &
                                  db%subsample_step,' = ', &
                                  db%dt*dble(db%subsample_step),' s'
  write(*,'(a,l1,a,l1,a,l1)')     '  topography       : ',db%topography, &
                                  '   ellipticity ',db%ellipticity, &
                                  '   rotation ',db%rotation
  write(*,'(a,i0)')               '  stations         : ',db%nstations

  write(*,*)
  write(*,'(a)') '    station         latitude    longitude   burial [m]'
  do i = 1,db%nstations
    write(*,'(a,a12,2f13.4,f12.1)') '    ',trim(db%stations(i)%id), &
      db%stations(i)%latitude,db%stations(i)%longitude,db%stations(i)%depth
  enddo

  write(*,*)
  write(*,*) '  (gf_print_info(db,6,.false.) prints all of the metadata,'
  write(*,*) '   which is what `xgf3d --info` shows)'

  nsta = db%nstations

!
!--- 2. the source ------------------------------------------------------
!

  call banner('2. the source, and where it sits')

  ! The file route. A Fortran caller may use the solver's own CMTSOLUTION
  ! reader directly -- it is re-exported here -- where the C and Python
  ! interfaces deliberately cannot, because get_cmt() stops on malformed
  ! input and a stop inside a shared object would kill the caller. Reading
  ! the file is safe from a program that owns its own process.
  call gf_read_cmt_source(CMT_PATH,db%dt,src,ierr)
  if (ierr /= GF_OK) then
    write(*,'(a,a)') '  could not read the CMTSOLUTION: ',trim(gf_errmsg)
    stop 1
  endif

  call gf_print_source(src,6)

  call tic()
  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  call toc(t_locate)

  if (ierr /= GF_OK) then
    write(*,'(a,a)') '  could not locate the source: ',trim(gf_errmsg)
    stop 1
  endif

  ! again, to show what the first one paid for
  call tic()
  call gf_locate_source(db,src%latitude,src%longitude,src%depth,loc,ierr)
  call toc(t_locate2)

  write(*,*)
  write(*,'(a,i0,a,a)')   '  element        : ',loc%ielem,'  ',trim(loc%morton_hex)
  write(*,'(a,3f14.9)')   '  xi, eta, gamma : ',loc%xi,loc%eta,loc%gamma
  write(*,'(a,es12.4,a)') '  |mapped - requested| : ',loc%distance_km,' km'
  write(*,'(a,es12.4)')   '  27-anchor residual   : ',loc%anchor_err
  write(*,'(a,f9.1,a)')   '  gf_locate_source took ',t_locate*1.d3, &
                          ' ms the first time,'
  write(*,'(a,f9.1,a)')   '                        ',t_locate2*1.d3, &
                          ' ms the second'
  write(*,*) '  (the first call builds the element search tree and loads the'
  write(*,*) '   topography grid; every later one reuses both)'

  if (max(abs(loc%xi),abs(loc%eta),abs(loc%gamma)) > 1.d0) then
    write(*,*)
    write(*,*) '  Note: one coordinate is outside [-1,1], so the source sits'
    write(*,*) '  just outside the element that holds it. That is not a bug:'
    write(*,*) '  it is where specfem itself placed this source, and'
    write(*,*) '  reproducing that run is what the extraction is judged on.'
  endif

!
!--- 3. the plan --------------------------------------------------------
!

  call banner('3. the output axis and the source time function')

  call gf_seis_plan(db,src,-1.d0,tax,stf,ierr)
  if (ierr /= GF_OK) then
    write(*,'(a,a)') '  could not plan: ',trim(gf_errmsg)
    stop 1
  endif

  ! -1 above asks for specfem's own rule for a forward run, 1.5*hdur before
  ! the centroid time for a moment tensor. Pass a positive number to choose.
  call gf_print_stf(stf,tax,6)

  nt = tax%nt

!
!--- 4. extraction ------------------------------------------------------
!

  call banner('4. seismograms and partial derivatives')

  ! get_seismograms() is GF3DF's name and argument shape: it locates, plans
  ! and extracts in one call, and allocates its own outputs. itypsokern = 2
  ! asks for all ten partials; 1 gives the six moment-tensor ones, 0 none.
  call tic()
  call get_seismograms(db,src,synt,dp,2,t,ierr)
  call toc(t_extract)

  if (ierr /= GF_OK) then
    write(*,'(a,a)') '  extraction failed: ',trim(gf_errmsg)
    stop 1
  endif

  write(*,'(a,f9.1,a)') '  get_seismograms took ',t_extract*1.d3,' ms'
  write(*,*)
  write(*,'(a,i0,a,i0,a,i0,a)') '  synt(', size(synt,1),',',size(synt,2),',', &
                                size(synt,3),')   stations, components N/E/Z, samples'
  write(*,'(a,i0,a,i0,a,i0,a,i0,a)') '  dp  (', size(dp,1),',',size(dp,2),',', &
                                size(dp,3),',',size(dp,4),')  parameters first'
  write(*,'(a,i0,a,f12.4,a,f12.4,a)') '  t   (',size(t),')   from ',t(1), &
                                ' s to ',t(size(t)),' s'
  write(*,*) '  (t = 0 is the centroid time: the CMTSOLUTION time shift is'
  write(*,*) '   header metadata and is not in the trace)'

  write(*,*)
  write(*,*) '  peak displacement per station and component [m]'
  write(*,'(a)') '      station               N            E            Z'
  do ista = 1,nsta
    write(*,'(a,a12,3es13.4)') '    ',adjustr(db%stations(ista)%id(1:12)), &
      (maxval(abs(synt(ista,icomp,:))),icomp = 1,3)
  enddo

  write(*,*)
  write(*,*) '  the ten partials, their units, and their peaks'
  do ip = 1,GF_NDP_LOC
    write(*,'(a,i2,a,a3,a,a9,a,es12.4)') '    ',ip,'   ',GF_DP_NAME(ip), &
      '   ',GF_DP_UNIT(ip),'   ',maxval(abs(dp(ip,:,:,:)))
  enddo

  ! -- the identity the moment-tensor partials satisfy --------------------
  !
  ! The first six are per dyne-cm and the seismogram is linear in the moment
  ! tensor, so contracting them with the CMTSOLUTION's own components has to
  ! give the seismogram back. src%moment_tensor is non-dimensional --
  ! get_cmt divided by scaleM -- so scale_moment puts it back into dyne-cm.
  worst = 0.d0
  peak = maxval(abs(synt))
  do ista = 1,nsta
    do icomp = 1,3
      do i = 1,nt
        worst = max(worst,abs(sum(src%moment_tensor(1:6)*src%scale_moment &
                                  *dp(1:6,ista,icomp,i)) - synt(ista,icomp,i)))
      enddo
    enddo
  enddo

  write(*,*)
  write(*,*) '  sum over v of M_v * dp(v) against the seismogram:'
  write(*,'(a,es12.4,a)') '    max difference ',worst/peak,' of the trace peak'
  write(*,*) '  So a new mechanism at the same hypocentre is a contraction of'
  write(*,*) '  arrays already in memory, not another extraction.'

  if (worst/peak > 1.d-12) ok = .false.

!
!--- 5. the centroid partials are derivatives ---------------------------
!

  call banner('5. the centroid partials are derivatives')

  write(*,*) '  Relocating the source changes the seismogram non-linearly, so'
  write(*,*) '  here the partial is a genuine derivative and one Taylor term'
  write(*,*) '  is only an approximation:'
  write(*,*)
  write(*,*) '      u(lat + h)  ~=  u(lat) + h * du/dlat'
  write(*,*)
  write(*,*) '  Its error must fall as h^2, so halving h divides the error by'
  write(*,*) '  four. A partial with the wrong sign, units or scaling cannot'
  write(*,*) '  do that, however plausible its traces look.'
  write(*,*)

  do k = 0,NHALVE
    step(k) = STEP0 / dble(2**k)

    moved = src
    moved%latitude = src%latitude + step(k)

    call tic()
    call get_seismograms(db,moved,truth,dp_unused,0,t_unused,ierr)
    call toc(t_reloc(k))

    if (ierr /= GF_OK) then
      write(*,'(a,f9.4,a,a)') '    step ',step(k),' could not be extracted: ', &
        trim(gf_error_string(ierr))
      write(*,*) '    (the relocated source left the database)'
      err(k) = -1.d0
      cycle
    endif

    ! predicted from the *unperturbed* run: synt + h * dp(lat)
    err(k) = 0.d0
    do ista = 1,nsta
      do icomp = 1,3
        do i = 1,nt
          err(k) = max(err(k),abs(synt(ista,icomp,i) &
                   + step(k)*dp(GF_DP_LAT,ista,icomp,i) - truth(ista,icomp,i)))
        enddo
      enddo
    enddo
    err(k) = err(k)/peak

    deallocate(truth)
  enddo

  write(*,'(a)') '      step [deg]    error / peak    ratio'
  do k = 0,NHALVE
    if (err(k) < 0.d0) cycle
    if (k > 0 .and. err(max(k-1,0)) > 0.d0) then
      ratio = err(k-1)/err(k)
      write(*,'(a,f11.5,es16.4,f10.2)') '    ',step(k),err(k),ratio
    else
      write(*,'(a,f11.5,es16.4)') '    ',step(k),err(k)
    endif
  enddo

  if (err(0) > 0.d0 .and. err(NHALVE) > 0.d0) then
    ratio = (err(0)/err(NHALVE))**(1.d0/dble(NHALVE))
    write(*,*)
    if (ratio > 3.d0 .and. ratio < 5.d0) then
      write(*,'(a,f5.2,a)') '  mean ratio ',ratio, &
        ': second order, as a first derivative must be'
    else
      write(*,'(a,f5.2,a)') '  mean ratio ',ratio, &
        ': NOT second order -- something is wrong'
      ok = .false.
    endif
  endif

!
!--- timings ------------------------------------------------------------
!

  call banner('timings')

  t_total = t_open + t_locate + t_locate2 + t_extract + sum(t_reloc)

  call timing_line('gf_open',t_open)
  call timing_line('gf_locate_source, first',t_locate)
  call timing_line('gf_locate_source, again',t_locate2)
  call timing_line('get_seismograms, 10 partials',t_extract)
  do k = 0,NHALVE
    write(label,'(a,f7.4,a)') 'get_seismograms, at +',step(k),' deg'
    call timing_line(label,t_reloc(k))
  enddo
  write(*,'(a)') '                                    ----------'
  call timing_line('in the library',t_total)
  write(*,*)
  write(*,'(a,i0,a,i0,a)') '  Every extraction above produced ',nsta, &
    ' stations x 3 components x ',nt,' samples.'
  write(*,*) '  The database is opened once and stays open; the reciprocal'
  write(*,*) '  simulations that built it took about forty minutes.'

!
!--- done ---------------------------------------------------------------
!

  deallocate(synt,dp,t)

  ! the pair a long-lived caller needs: gf_close alone leaves the search
  ! tree allocated, because gf_locate cannot be unwound from gf_database
  call gf_release(db)

  write(*,*)
  if (ok) then
    write(*,*) 'gf3d_api_demo: done'
  else
    write(*,*) 'gf3d_api_demo: something did not agree, see above'
    stop 1
  endif
  write(*,*)

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine timing_line(label,seconds)

! one row of the timing table, in a fixed column

  implicit none

  character(len=*), intent(in) :: label
  double precision, intent(in) :: seconds

  ! local parameters
  character(len=34) :: padded

  padded = label

  write(*,'(a,a,f10.1,a)') '  ',padded,seconds*1.d3,' ms'

  end subroutine timing_line

!
!-------------------------------------------------------------------------------------------------
!

  subroutine banner(title)

  implicit none

  character(len=*), intent(in) :: title

  write(*,*)
  write(*,'(a,a)') ' ',title
  write(*,'(a,a)') ' ',repeat('-',len_trim(title))

  end subroutine banner

!
!-------------------------------------------------------------------------------------------------
!

  subroutine tic()

  implicit none

  call system_clock(clock_start,clock_rate)

  end subroutine tic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine toc(seconds)

  implicit none

  double precision, intent(out) :: seconds

  ! local parameters
  integer(kind=8) :: clock_end

  call system_clock(clock_end)

  seconds = dble(clock_end - clock_start)/dble(clock_rate)

  end subroutine toc

  end program gf3d_api_demo
