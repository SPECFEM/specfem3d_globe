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
!---- xgf3d — the Green function extraction tool.
!----
!---- One executable with several modes rather than a family of small
!---- binaries, so that every mode survives as a debugging entry point:
!----
!----   xgf3d --version
!----   xgf3d --info <GFDB> [--topo] [--no-check]
!----
!---- Later stages add --locate, --seis, --dump and the SAC-writing
!---- default mode.
!----

  program xgf3d

  use gf_par, only: t_gfdb,gf_errmsg,gf_error_string,GF_OK,GF3D_VERSION
  use gf_database, only: gf_open,gf_close,gf_print_info

  use constants, only: MAX_STRING_LEN,NGLLX,NGNOD

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: arg,mode,dbpath
  type(t_gfdb) :: db
  integer :: nargs,iarg,ierr
  logical :: with_topo,do_check

  ! standard error, used for anything that is not the requested output
  integer, parameter :: ISTDERR = 0

  nargs = command_argument_count()

  if (nargs < 1) then
    call print_usage(ISTDERR)
    stop 1
  endif

  call get_command_argument(1,mode)

  select case (trim(mode))

  case ('--version','-v')
    write(*,'(a)') 'xgf3d version '//GF3D_VERSION//' (specfem3d_globe)'
#ifdef USE_HDF5
    write(*,'(a)') '  HDF5 support : yes'
#else
    write(*,'(a)') '  HDF5 support : no  (re-run configure with --with-hdf5)'
#endif
    write(*,'(a,i0,a,i0)') '  built with     NGLLX = ',NGLLX,'  NGNOD = ',NGNOD
    stop

  case ('--help','-h')
    call print_usage(6)
    stop

  case ('--info')

    if (nargs < 2) then
      write(ISTDERR,'(a)') 'Error: --info needs a database directory'
      call print_usage(ISTDERR)
      stop 1
    endif

    call get_command_argument(2,dbpath)

    with_topo = .false.
    do_check = .true.
    do iarg = 3,nargs
      call get_command_argument(iarg,arg)
      select case (trim(arg))
      case ('--topo')
        with_topo = .true.
      case ('--no-check')
        do_check = .false.
      case default
        write(ISTDERR,'(a)') 'Error: unknown option for --info: '//trim(arg)
        stop 1
      end select
    enddo

    call gf_open(dbpath,db,ierr,check_completion=do_check)
    if (ierr /= GF_OK) then
      write(ISTDERR,'(a)') 'Error opening the Green function database'
      write(ISTDERR,'(a)') '  '//trim(gf_error_string(ierr))//': '//trim(gf_errmsg)
      stop 1
    endif

    call gf_print_info(db,6,with_topo)

    call gf_close(db)

  case default
    write(ISTDERR,'(a)') 'Error: unknown mode: '//trim(mode)
    call print_usage(ISTDERR)
    stop 1

  end select

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine print_usage(iunit)

  implicit none

  integer, intent(in) :: iunit

  write(iunit,'(a)') 'usage: xgf3d <mode> [arguments]'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  --version                       print the version and build settings'
  write(iunit,'(a)') '  --help                          print this message'
  write(iunit,'(a)') '  --info <GFDB> [options]         describe a Green function database'
  write(iunit,'(a)') ''
  write(iunit,'(a)') '  options for --info:'
  write(iunit,'(a)') '    --topo       also load the topography grid and probe it at each station'
  write(iunit,'(a)') '    --no-check   skip the per-element completion scan'
  write(iunit,'(a)') ''

  end subroutine print_usage

  end program xgf3d
