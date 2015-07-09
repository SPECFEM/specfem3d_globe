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


  subroutine create_addressing(myrank,NCHUNKS,NPROC,NPROC_ETA,NPROC_XI,NPROCTOT, &
                        addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
                        OUTPUT_FILES)

  use constants

  implicit none

  integer :: myrank,NCHUNKS,NPROC,NPROC_ETA,NPROC_XI,NPROCTOT

  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1) :: addressing
  integer, dimension(0:NPROCTOT-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice

  character(len=MAX_STRING_LEN) :: OUTPUT_FILES

  ! local parameters
  integer :: ichunk,iproc_eta,iproc_xi,iprocnum,ier

  ! initializes
  addressing(:,:,:) = 0
  ichunk_slice(:) = 0
  iproc_xi_slice(:) = 0
  iproc_eta_slice(:) = 0

  ! loop on all the chunks to create global slice addressing for solver
  if (myrank == 0) then
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/addressing.txt',status='unknown',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening addressing.txt')
    write(IMAIN,*) 'creating global slice addressing'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  do ichunk = 1,NCHUNKS
    do iproc_eta = 0,NPROC_ETA-1
      do iproc_xi = 0,NPROC_XI-1
        iprocnum = (ichunk-1)*NPROC + iproc_eta * NPROC_XI + iproc_xi
        addressing(ichunk,iproc_xi,iproc_eta) = iprocnum
        ichunk_slice(iprocnum) = ichunk
        iproc_xi_slice(iprocnum) = iproc_xi
        iproc_eta_slice(iprocnum) = iproc_eta
        if (myrank == 0) write(IOUT,*) iprocnum,ichunk,iproc_xi,iproc_eta
      enddo
    enddo
  enddo

  if (myrank == 0) close(IOUT)

  ! output a topology map of slices - fix 20x by nproc
  if (myrank == 0) then
    if (NCHUNKS == 6 .and. NPROCTOT < 1000) then
      write(IMAIN,*) 'Spatial distribution of the slices'
      do iproc_xi = NPROC_XI-1, 0, -1
        write(IMAIN,'(20x)',advance='no')
        do iproc_eta = NPROC_ETA -1, 0, -1
          ichunk = CHUNK_AB
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      do iproc_xi = NPROC_XI-1, 0, -1
        write(IMAIN,'(1x)',advance='no')
        do iproc_eta = NPROC_ETA -1, 0, -1
          ichunk = CHUNK_BC
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(3x)',advance='no')
        do iproc_eta = NPROC_ETA -1, 0, -1
          ichunk = CHUNK_AC
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(3x)',advance='no')
        do iproc_eta = NPROC_ETA -1, 0, -1
          ichunk = CHUNK_BC_ANTIPODE
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      do iproc_xi = NPROC_XI-1, 0, -1
        write(IMAIN,'(20x)',advance='no')
        do iproc_eta = NPROC_ETA -1, 0, -1
          ichunk = CHUNK_AB_ANTIPODE
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      do iproc_xi = NPROC_XI-1, 0, -1
        write(IMAIN,'(20x)',advance='no')
        do iproc_eta = NPROC_ETA -1, 0, -1
          ichunk = CHUNK_AC_ANTIPODE
          write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
        enddo
        write(IMAIN,'(1x)',advance='yes')
      enddo
      write(IMAIN, *) ' '
      call flush_IMAIN()
    endif
  endif

  end subroutine create_addressing
