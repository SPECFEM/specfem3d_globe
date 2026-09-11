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
!---- Seismograms for several events from one open database.
!----
!---- The Fortran counterpart of ../python/gf3d_events_demo.py.
!----
!---- The global example's database was not built for a single earthquake.
!---- Its db_base/DATA/GF_LOCATIONS names three hypocentres -- the 2010
!---- Chile (Maule) event, the 2019 Peru event, and a point near the centre
!---- of the chunk -- and the reciprocal runs stored the strain around all
!---- of them. That is what a Green function database is for: the expensive
!---- simulations are done per *station*, once, and any source inside the
!---- covered volume is then a look-up.
!----
!---- So this program opens the database once and extracts for every
!---- hypocentre the database itself declares, timing each, and prints the
!---- peak amplitudes. The moment tensor is held fixed at the validation
!---- event's, so the three records differ only by where the source sits.
!---- A real catalogue would vary the mechanism too, and that costs nothing
!---- extra: the moment-tensor partials are exact, so a new mechanism at the
!---- same hypocentre is a contraction of arrays already in memory rather
!---- than another extraction. gf3d_api_demo shows that identity.
!----

  program gf3d_events_demo

  use gf3d

  implicit none

  character(len=*), parameter :: DB_PATH = '../../GFDB'
  character(len=*), parameter :: CMT_PATH = '../../validation_data/CMTSOLUTION'
  character(len=*), parameter :: LOC_PATH = '../../db_base/DATA/GF_LOCATIONS'

  integer, parameter :: MAX_EVENTS = 64

  type(t_gfdb) :: db
  type(t_gf_source) :: src,event
  type(t_gf_location) :: loc

  double precision, dimension(:,:,:), allocatable :: synt
  double precision, dimension(:,:,:,:), allocatable :: dp_unused
  double precision, dimension(:), allocatable :: t

  ! the hypocentres, as GF_LOCATIONS gives them
  character(len=64), dimension(MAX_EVENTS) :: ev_name
  double precision, dimension(MAX_EVENTS) :: ev_lat,ev_lon,ev_dep
  integer :: nev

  ! what each one measured
  double precision, dimension(MAX_EVENTS) :: t_locate,t_extract,peak
  double precision, dimension(:,:), allocatable :: peak_sta
  integer, dimension(MAX_EVENTS) :: ielem
  logical, dimension(MAX_EVENTS) :: done

  character(len=32) :: hdr

  double precision :: t_open,sum_extract
  integer :: ierr,k,ista,nsta,nt,nok

  integer(kind=8) :: clock_start,clock_rate

  write(*,*)
  write(*,*) '========================================================'
  write(*,*) ' gf3d: several events, one open database'
  write(*,*) '========================================================'
  write(*,'(a,a)') '  library version : ',trim(GF3D_VERSION)

!
!--- the hypocentres the database declares ------------------------------
!

  call read_locations(LOC_PATH,ev_name,ev_lat,ev_lon,ev_dep,nev,ierr)
  if (ierr /= 0) then
    write(*,'(a,a)') '  could not read ',LOC_PATH
    stop 1
  endif
  if (nev < 1) then
    write(*,'(a,a)') '  no hypocentres declared in ',LOC_PATH
    stop 1
  endif

  write(*,'(a,i0,a,a)') '  events          : ',nev,' from ',LOC_PATH

!
!--- open, once ---------------------------------------------------------
!

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

  nsta = db%nstations
  allocate(peak_sta(nev,nsta))
  peak_sta(:,:) = 0.d0

  write(*,'(a,i0,a,10(a,a))') '  stations        : ',nsta,':  ', &
    (trim(db%stations(ista)%id),'  ',ista = 1,min(nsta,10))
  write(*,'(a,f9.1,a)') '  gf_open took ',t_open*1.d3,' ms'

!
!--- the mechanism, from the validation event ---------------------------
!

  call gf_read_cmt_source(CMT_PATH,db%dt,src,ierr)
  if (ierr /= GF_OK) then
    write(*,'(a,a)') '  could not read the CMTSOLUTION: ',trim(gf_errmsg)
    stop 1
  endif

  write(*,'(a,a,a)') '  mechanism held fixed at ',trim(src%event_name), &
                     '''s, so only the hypocentre differs'

!
!--- one extraction per hypocentre --------------------------------------
!

  ! a left-justified header: Fortran right-justifies a character value in a
  ! wider aW field, which would put this over the numbers instead
  hdr = 'event'

  write(*,*)
  write(*,'(a,a32,3a11,a7,2a11,a13)') '  ',hdr,'lat','lon','depth', &
    'elem','locate','extract','peak [m]'
  write(*,'(a)') '  ' // repeat('-',107)

  nok = 0
  nt = 0

  do k = 1,nev
    done(k) = .false.

    ! same source, moved. t_gf_source is a plain derived type, so this is
    ! the whole of "make me another event at a different place".
    event = src
    event%latitude = ev_lat(k)
    event%longitude = ev_lon(k)
    event%depth = ev_dep(k)

    call tic()
    call gf_locate_source(db,ev_lat(k),ev_lon(k),ev_dep(k),loc,ierr)
    call toc(t_locate(k))

    if (ierr /= GF_OK) then
      write(*,'(a,a32,3f11.3,a,a)') '  ',ev_name(k),ev_lat(k),ev_lon(k),ev_dep(k), &
        '      -- ',trim(gf_error_string(ierr))
      cycle
    endif

    ! itypsokern = 0: seismograms only, no partials
    call tic()
    call get_seismograms(db,event,synt,dp_unused,0,t,ierr)
    call toc(t_extract(k))

    if (ierr /= GF_OK) then
      write(*,'(a,a32,3f11.3,a,a)') '  ',ev_name(k),ev_lat(k),ev_lon(k),ev_dep(k), &
        '      -- ',trim(gf_error_string(ierr))
      cycle
    endif

    ielem(k) = loc%ielem
    peak(k) = maxval(abs(synt))
    do ista = 1,nsta
      peak_sta(k,ista) = maxval(abs(synt(ista,:,:)))
    enddo
    nt = size(t)
    done(k) = .true.
    nok = nok + 1

    write(*,'(a,a32,3f11.3,i7,f8.1,a,f8.1,a,es13.4)') '  ',ev_name(k), &
      ev_lat(k),ev_lon(k),ev_dep(k),ielem(k), &
      t_locate(k)*1.d3,' ms',t_extract(k)*1.d3,' ms',peak(k)

    deallocate(synt,t)
    if (allocated(dp_unused)) deallocate(dp_unused)
  enddo

  if (nok == 0) then
    write(*,*)
    write(*,*) '  none of the declared hypocentres is inside this database'
    call gf_release(db)
    stop 1
  endif

!
!--- peaks per station --------------------------------------------------
!

  write(*,*)
  write(*,*) '  peak amplitude per station [m]'
  write(*,'(a,a32,10a14)') '    ',hdr, &
    (trim(db%stations(ista)%id),ista = 1,nsta)
  do k = 1,nev
    if (.not. done(k)) cycle
    write(*,'(a,a32,10es14.4)') '    ',ev_name(k),(peak_sta(k,ista),ista = 1,nsta)
  enddo

!
!--- what it cost -------------------------------------------------------
!

  sum_extract = 0.d0
  do k = 1,nev
    if (done(k)) sum_extract = sum_extract + t_extract(k)
  enddo

  write(*,*)
  write(*,'(a,i0,a,f6.1,a)') '  ',nok,' events: ',t_open*1.d3,' ms to open,'
  write(*,'(a,f6.1,a)')      '  ',t_locate(1)*1.d3, &
    ' ms for the first locate (it builds the element search'
  write(*,*) '  tree and loads the topography grid), then'
  write(*,'(a,f6.1,a,i0,a,i0,a)') '  ',sum_extract*1.d3/dble(nok), &
    ' ms per event on average, each ',nsta,' stations x 3 components x ',nt, &
    ' samples.'
  write(*,*)
  write(*,*) '  The reciprocal simulations that built this database took about'
  write(*,*) '  forty minutes. Every source inside it is now a look-up.'

  deallocate(peak_sta)

  call gf_release(db)

  write(*,*)
  write(*,*) 'gf3d_events_demo: done'
  write(*,*)

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_locations(filename,name,lat,lon,dep,n,ierr)

! the hypocentres a GF_LOCATIONS file declares, with their labels
!
! Three columns -- latitude longitude depth_km -- with `#` comments. The
! comment immediately above a row is taken as that row's name, which is how
! the shipped file labels its events.

  implicit none

  character(len=*), intent(in) :: filename
  character(len=64), dimension(:), intent(out) :: name
  double precision, dimension(:), intent(out) :: lat,lon,dep
  integer, intent(out) :: n,ierr

  ! local parameters
  integer, parameter :: IIN_LOC = 71
  character(len=256) :: line
  character(len=64) :: label
  integer :: ios
  double precision :: a,b,c

  n = 0
  ierr = 0
  label = ''

  open(unit=IIN_LOC,file=trim(filename),status='old',action='read',iostat=ios)
  if (ios /= 0) then
    ierr = 1
    return
  endif

  do
    read(IIN_LOC,'(a)',iostat=ios) line
    if (ios /= 0) exit

    line = adjustl(line)
    if (len_trim(line) == 0) cycle

    if (line(1:1) == '#') then
      ! remember it: if a data row follows, this is its name
      label = adjustl(line(2:))
      cycle
    endif

    read(line,*,iostat=ios) a,b,c
    if (ios /= 0) cycle

    if (n >= size(lat)) exit

    n = n + 1
    lat(n) = a
    lon(n) = b
    dep(n) = c
    if (len_trim(label) > 0) then
      name(n) = label
    else
      write(name(n),'(a,i0)') 'event ',n
    endif
    label = ''
  enddo

  close(IIN_LOC)

  end subroutine read_locations

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

  end program gf3d_events_demo
