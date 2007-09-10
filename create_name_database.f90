!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, October 2007
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
  character(len=150) prname,procname,LOCAL_PATH

! create the name for the database of the current slide and region
  write(procname,"('/proc',i6.6,'_reg',i1,'_')") iproc,iregion_code

! create full name with path
  prname = trim(LOCAL_PATH) // procname

  end subroutine create_name_database

