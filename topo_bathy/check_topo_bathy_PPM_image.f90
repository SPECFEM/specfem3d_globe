!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 2
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!          (c) California Institute of Technology July 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  program check_topo_bathy_PPM_image
!
!--- read topography and bathymetry ASCII file and create PPM image to check it
!
!--- compile and link with file topo_bathy.f90 from main SEM code
!
  implicit none

  include "../../constants.h"

! use integer array to store values
  integer ibathy_topo(NX_BATHY,NY_BATHY)

  integer grey(NX_BATHY,NY_BATHY)
  integer ivalue,icurrent_rec,ix,iz

! read the topography file
  print *,'reading topo file'
  print *,'file used has a resolution of ',RESOLUTION_TOPO_FILE
  call read_topo_bathy_file(ibathy_topo)

  print *,'min and max of topography = ',minval(ibathy_topo),maxval(ibathy_topo)

! creating image with grey levels
  grey = 255 * (ibathy_topo - minval(ibathy_topo)) / (maxval(ibathy_topo) - minval(ibathy_topo))
  where(grey < 1) grey = 1
  where(grey > 255) grey = 255

! store image in PPM format with grey levels

! creating the header
  open(unit=27,file='____tutu1',status='unknown')
  write(27,100)
  write(27,101)
  write(27,102) NX_BATHY,NY_BATHY
  write(27,103)
  close(27)

 100 format('P5')
 101 format('# creator DK')
 102 format(i6,' ',i6)
 103 format('255')

! create the PPM image
  print *,'creating PPM image'
  open(unit=27,file='____tutu2',status='unknown',access='direct',recl=1)

  icurrent_rec = 1

  do iz = 1,NY_BATHY
  do ix = 1,NX_BATHY

! write grey
      write(27,rec=icurrent_rec) achar(grey(ix,iz))
      icurrent_rec = icurrent_rec + 1

  enddo
  enddo

  close(27)

  call system('cat ____tutu1 ____tutu2 > image_topo_bathy.ppm ; rm -f ____tutu1 ____tutu2')

  end program check_topo_bathy_PPM_image

