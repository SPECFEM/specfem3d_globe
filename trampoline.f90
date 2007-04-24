!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program trampoline

! Perform Fortran mojo, and then run the Python script.

! With ifort v9 in particular, this function (i.e., MAIN__) will call
! the undocumented function __intel_new_proc_init or
! __intel_new_proc_init_P.  Without this, SPECFEM runs several
! times slower (!).

  call FC_PY_MAIN()

  end program trampoline
