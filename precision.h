!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology August 2003
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
! solver in single or double precision depending on the machine
!
!  ALSO CHANGE FILE  constants.h ACCORDINGLY
!
! uncomment this to run in single precision
  integer, parameter :: CUSTOM_MPI_TYPE = MPI_REAL
! uncomment this to run in double precision
! integer, parameter :: CUSTOM_MPI_TYPE = MPI_DOUBLE_PRECISION

