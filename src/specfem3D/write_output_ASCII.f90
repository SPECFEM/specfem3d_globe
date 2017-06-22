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

  subroutine write_output_ASCII(seismogram_tmp,iorientation,sisname,sisname_big_file)

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

  use constants, only: CUSTOM_REAL,SIZE_REAL,IOUT,MAX_STRING_LEN

  use specfem_par, only: &
    DT,t0,NSTEP, &
    seismo_offset,seismo_current, &
    NTSTEP_BETWEEN_OUTPUT_SEISMOS,OUTPUT_FILES,SIMULATION_TYPE, &
    SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE, &
    myrank

  implicit none

  real(kind=CUSTOM_REAL), dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismogram_tmp

  integer :: iorientation
  character(len=MAX_STRING_LEN) :: sisname,sisname_big_file

  ! local parameters
  integer :: it
  integer :: ier,isample
  double precision :: value
  double precision :: timeval
  character(len=MAX_STRING_LEN) :: sisname_2

  ! add .ascii extension to seismogram file name for ASCII seismograms
  write(sisname_2,"('/',a,'.ascii')") trim(sisname)

  ! save all seismograms in one large combined file instead of one file per seismogram
  ! to avoid overloading shared non-local file systems

  ! create one large file instead of one small file per station to avoid file system overload
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE) then
    if (USE_BINARY_FOR_LARGE_FILE) then
      write(IOUT) sisname_big_file
    else
      write(IOUT,*) sisname_big_file(1:len_trim(sisname_big_file))
    endif
  else
    if (seismo_offset == 0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2), &
            status='unknown',action='write',iostat=ier)
    else
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2), &
            status='old',position='append',action='write',iostat=ier)
    endif
    if (ier /= 0) call exit_mpi(myrank,'Error opening file:'//trim(OUTPUT_FILES)//trim(sisname_2))
  endif

  ! subtract half duration of the source to make sure travel time is correct
  do isample = 1,seismo_current

    ! seismogram value
    value = dble(seismogram_tmp(iorientation,isample))

    ! current time increment
    it = seismo_offset + isample

    ! current time
    if (SIMULATION_TYPE == 3) then
      timeval = dble(NSTEP-it)*DT - t0
    else
      timeval = dble(it-1)*DT - t0
    endif

    ! writes out to file
    if (SAVE_ALL_SEISMOS_IN_ONE_FILE .and. USE_BINARY_FOR_LARGE_FILE) then
      ! distinguish between single and double precision for reals
      write(IOUT) real(timeval, kind=CUSTOM_REAL), real(value, kind=CUSTOM_REAL)
    else
      ! distinguish between single and double precision for reals
      write(IOUT,*) real(timeval, kind=CUSTOM_REAL), ' ', real(value, kind=CUSTOM_REAL)
    endif
  enddo

  if (.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

  end subroutine write_output_ASCII
