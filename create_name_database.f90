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

  subroutine create_name_database(prname,iproc,iregion_code,LOCAL_PATH)

! create the name of the database for the mesher and the solver

  implicit none

  integer iproc,iregion_code

! name of the database file
  character(len=150) prname,procname,LOCAL_PATH,clean_LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,10) iproc,iregion_code
 10 format('/proc',i4.4,'_reg',i1,'_')

! suppress white spaces if any
  clean_LOCAL_PATH = adjustl(LOCAL_PATH)

! create full name with path
  prname = clean_LOCAL_PATH(1:len_trim(clean_LOCAL_PATH)) // procname

  end subroutine create_name_database

