!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2011
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

  subroutine write_output_ASCII(seismogram_tmp, &
              DT,hdur,OUTPUT_FILES, &
              NTSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current, &
              SAVE_ALL_SEISMOS_IN_ONE_FILE,USE_BINARY_FOR_LARGE_FILE,myrank, &
              iorientation,sisname,sisname_big_file)

! save seismograms in text format with no subsampling.
! Because we do not subsample the output, this can result in large files
! if the simulation uses many time steps. However, subsampling the output
! here would result in a loss of accuracy when one later convolves
! the results with the source time function

  implicit none

  include "constants.h"

  integer :: seismo_offset, seismo_current, NTSTEP_BETWEEN_OUTPUT_SEISMOS

  real(kind=CUSTOM_REAL), dimension(5,NTSTEP_BETWEEN_OUTPUT_SEISMOS) :: seismogram_tmp

  integer myrank
  double precision hdur,DT

  integer iorientation

  character(len=256) sisname,sisname_big_file
  character(len=150) OUTPUT_FILES

  ! save all seismograms in one large combined file instead of one file per seismogram
  ! to avoid overloading shared non-local file systems such as GPFS for instance
  logical SAVE_ALL_SEISMOS_IN_ONE_FILE
  logical USE_BINARY_FOR_LARGE_FILE

  ! local parameters
  integer ier,isample
  character(len=256) sisname_2
  double precision value


  ! add .ascii extension to seismogram file name for ASCII seismograms
  write(sisname_2,"('/',a,'.ascii')") trim(sisname)

  ! create one large file instead of one small file per station to avoid file system overload
  if(SAVE_ALL_SEISMOS_IN_ONE_FILE) then
    if(USE_BINARY_FOR_LARGE_FILE) then
      write(IOUT) sisname_big_file
    else
      write(IOUT,*) sisname_big_file(1:len_trim(sisname_big_file))
    endif
  else
    if (seismo_offset==0) then
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2), &
            status='unknown',action='write',iostat=ier)
    else
      open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(sisname_2), &
            status='old',position='append',action='write',iostat=ier)
    endif
    if( ier /= 0 ) call exit_mpi(myrank,'error opening file:'//trim(OUTPUT_FILES)//trim(sisname_2))
  endif

  ! subtract half duration of the source to make sure travel time is correct
  do isample = 1,seismo_current
    value = dble(seismogram_tmp(iorientation,isample))

    if(SAVE_ALL_SEISMOS_IN_ONE_FILE .and. USE_BINARY_FOR_LARGE_FILE) then
      ! distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        write(IOUT) sngl(dble(seismo_offset+isample-1)*DT - hdur),sngl(value)
      else
        write(IOUT) dble(seismo_offset+isample-1)*DT - hdur,value
      endif
    else
      ! distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        write(IOUT,*) sngl(dble(seismo_offset+isample-1)*DT - hdur),' ',sngl(value)
      else
        write(IOUT,*) dble(seismo_offset+isample-1)*DT - hdur,' ',value
      endif
    endif
  enddo

  if(.not. SAVE_ALL_SEISMOS_IN_ONE_FILE) close(IOUT)

  end subroutine write_output_ASCII
