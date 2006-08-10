!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 5
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology July 2004
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
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
