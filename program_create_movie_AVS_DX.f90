!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
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

!
!---  create a movie of radial component of surface displacement
!---  in AVS or OpenDX format
!

  program xcreate_movie_AVS_DX

  implicit none

  integer it1,it2
  integer iformat

! parameters read from parameter file
  integer NEX_XI,NEX_ETA
  integer NSTEP,NTSTEP_BETWEEN_FRAMES,NCHUNKS
  integer NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  logical MOVIE_SURFACE

! ************** PROGRAM STARTS HERE **************

  call read_AVS_DX_parameters(NEX_XI,NEX_ETA, &
           NSTEP,NTSTEP_BETWEEN_FRAMES, &
           NCHUNKS,MOVIE_SURFACE, &
           NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  if(.not. MOVIE_SURFACE) stop 'movie frames were not saved by the solver'

  print *,'1 = create files in OpenDX format'
  print *,'2 = create files in AVS UCD format with individual files'
  print *,'3 = create files in AVS UCD format with one time-dependent file'
  print *,'4 = create files in GMT xyz Ascii long/lat/Uz format'
  print *,'any other value = exit'
  print *
  print *,'enter value:'
  read(5,*) iformat
  if(iformat<1 .or. iformat>4) stop 'exiting...'

  print *,'movie frames have been saved every ',NTSTEP_BETWEEN_FRAMES,' time steps'
  print *

  print *,'enter first time step of movie (e.g. 1)'
  read(5,*) it1

  print *,'enter last time step of movie (e.g. ',NSTEP,')'
  read(5,*) it2

! run the main program
  call create_movie_AVS_DX(iformat,it1,it2, &
           NEX_XI,NEX_ETA, &
           NSTEP,NTSTEP_BETWEEN_FRAMES, &
           NCHUNKS, &
           NPROCTOT,NEX_PER_PROC_XI,NEX_PER_PROC_ETA)

  end program xcreate_movie_AVS_DX
