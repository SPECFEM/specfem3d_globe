!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  double precision function comp_source_spectrum(om,hdur)

  implicit none

  double precision om,hdur

! source decay rate (also change in source time function if needed)
  double precision, parameter :: decay_rate = 2.628d0

  comp_source_spectrum = dexp(-0.25d0*(om*hdur/decay_rate)**2)

  end function comp_source_spectrum
