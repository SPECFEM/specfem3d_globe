!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 4
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

! save header file OUTPUT_FILES/values_from_mesher.h

  subroutine save_header_file(NSPEC_AB,NSPEC_AC,NSPEC_BC, &
        nglob_AB,nglob_AC,nglob_BC,NEX_XI,NEX_ETA, &
        nspec_aniso_mantle,NPROC,NPROCTOT, &
        TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE, &
        ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION)

  implicit none

  include "constants.h"

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC_AB,NSPEC_AC,NSPEC_BC, &
               nglob_AB,nglob_AC,nglob_BC

  integer NEX_XI,NEX_ETA,NPROC,NPROCTOT
  integer nspec_aniso_mantle

  logical TRANSVERSE_ISOTROPY,ANISOTROPIC_MANTLE,ANISOTROPIC_INNER_CORE, &
          ELLIPTICITY,GRAVITY,ROTATION,ATTENUATION

  integer subtract_central_cube_elems
  double precision subtract_central_cube_points

! copy number of elements and points in an include file for the solver
  open(unit=IOUT,file='OUTPUT_FILES/values_from_mesher.h',status='unknown')
  write(IOUT,*)

  write(IOUT,*) '!'
  write(IOUT,*) '! this is the parameter file for static compilation of the solver'
  write(IOUT,*) '!'
  write(IOUT,*) '! mesh statistics:'
  write(IOUT,*) '! ---------------'
  write(IOUT,*) '!'

! the central cube is counted 6 times, therefore remove 5 times
  if(INCLUDE_CENTRAL_CUBE) then
    write(IOUT,*) '! these statistics include the central cube'
    subtract_central_cube_elems = 5 * (NEX_XI/8)**3
    subtract_central_cube_points = 5.d0 * (dble(NEX_XI/8)*dble(NGLLX-1)+1.d0)**3
  else
    write(IOUT,*) '! these statistics do not include the central cube'
    subtract_central_cube_elems = 0
    subtract_central_cube_points = 0.d0
  endif

  write(IOUT,*) '!'
  write(IOUT,*) '! number of processors = ',NPROCTOT
  write(IOUT,*) '!'
  write(IOUT,*) '! number of ES nodes = ',real(NPROCTOT)/8.
  write(IOUT,*) '! percentage of total 640 ES nodes = ',100.*(real(NPROCTOT)/8.)/640.,' %'
  write(IOUT,*) '! total memory available on these ES nodes (Gb) = ',16.*real(NPROCTOT)/8.

  write(IOUT,*) '!'
  write(IOUT,*) '! max points in largest region = ',nglob_BC(IREGION_CRUST_MANTLE)
! use fused loops on the ES
  write(IOUT,*) '! max vector length = ',nglob_BC(IREGION_CRUST_MANTLE)*NDIM
  write(IOUT,*) '! min vector length = ',NGLLSQUARE
  write(IOUT,*) '! min critical vector length = ',NGLLSQUARE_NDIM
  write(IOUT,*) '!'
  write(IOUT,*) '! on ES and SX-5, make sure "loopcnt=" parameter'
! use fused loops on the ES
  write(IOUT,*) '! in Makefile is greater than ',nglob_BC(IREGION_CRUST_MANTLE)*NDIM
  write(IOUT,*) '!'

  write(IOUT,*) '! total elements per AB slice = ',sum(NSPEC_AB)
  write(IOUT,*) '! total points per AB slice = ',sum(nglob_AB)
  write(IOUT,*) '!'
  write(IOUT,*) '! total elements per AC slice = ',sum(NSPEC_AC)
  write(IOUT,*) '! total points per AC slice = ',sum(nglob_AC)
  write(IOUT,*) '!'
  write(IOUT,*) '! total elements per BC slice = ',sum(NSPEC_BC)
  write(IOUT,*) '! total points per BC slice = ',sum(nglob_BC)
  write(IOUT,*) '!'
  write(IOUT,*) '! load balancing AB/BC for points = ',100.*real(sum(nglob_AB))/real(sum(nglob_BC)),' %'
  write(IOUT,*) '! load balancing AB/BC for elements = ',100.*real(sum(NSPEC_AB))/real(sum(NSPEC_BC)),' %'
  write(IOUT,*) '!'
  write(IOUT,*) '! load balancing AC/BC for points = ',100.*real(sum(nglob_AC))/real(sum(nglob_BC)),' %'
  write(IOUT,*) '! load balancing AC/BC for elements = ',100.*real(sum(NSPEC_AC))/real(sum(NSPEC_BC)),' %'
  write(IOUT,*) '!'

  write(IOUT,*) '! total for full 6-chunk mesh:'
  write(IOUT,*) '! ---------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! exact total number of spectral elements in entire mesh = '
  write(IOUT,*) '! ',2*NPROC*(sum(NSPEC_AB) + sum(NSPEC_AC) + sum(NSPEC_BC)) &
    - subtract_central_cube_elems
  write(IOUT,*) '! approximate total number of points in entire mesh = '
  write(IOUT,*) '! ',2.d0*dble(NPROC)*(dble(sum(nglob_AB)) + dble(sum(nglob_AC)) + dble(sum(nglob_BC))) &
    - subtract_central_cube_points
! there are 3 DOFs in solid regions, but only 1 in fluid outer core
  write(IOUT,*) '! approximate total number of degrees of freedom in entire mesh = '
  write(IOUT,*) '! ',2.d0*dble(NPROC)*(3.d0*(dble(sum(nglob_AB)) + &
    dble(sum(nglob_AC)) + dble(sum(nglob_BC))) &
    - 2.d0*dble(nglob_AB(IREGION_OUTER_CORE) + nglob_AC(IREGION_OUTER_CORE) + nglob_BC(IREGION_OUTER_CORE))) &
    - 3.d0*subtract_central_cube_points
  write(IOUT,*) '!'

  write(IOUT,*) '! resolution of the mesh at the surface:'
  write(IOUT,*) '! -------------------------------------'
  write(IOUT,*) '!'
  write(IOUT,*) '! spectral elements along a great circle = ',4*NEX_XI
  write(IOUT,*) '! GLL points along a great circle = ',4*NEX_XI*(NGLLX-1)
  write(IOUT,*) '! average distance between points in degrees = ',360./real(4*NEX_XI*(NGLLX-1))
  write(IOUT,*) '! average distance between points in km = ',real(TWO_PI*R_EARTH/1000.d0)/real(4*NEX_XI*(NGLLX-1))
  write(IOUT,*) '!'
  write(IOUT,*)

  if(NCHUNKS == 1) write(IOUT,*) '! values for AC and BC below undefined for one chunk'
  write(IOUT,*) 'integer, parameter :: NEX_XI_VAL = ',NEX_XI
  write(IOUT,*) 'integer, parameter :: NEX_ETA_VAL = ',NEX_ETA
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_AB = ',NSPEC_AB(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_AC = ',NSPEC_AC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_CRUST_MANTLE_BC = ',NSPEC_BC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_AB = ',NSPEC_AB(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_AC = ',NSPEC_AC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_OUTER_CORE_BC = ',NSPEC_BC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE = ',NSPEC_AB(IREGION_INNER_CORE)
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_AB = ',nglob_AB(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_AC = ',nglob_AC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_CRUST_MANTLE_BC = ',nglob_BC(IREGION_CRUST_MANTLE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_AB = ',nglob_AB(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_AC = ',nglob_AC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_OUTER_CORE_BC = ',nglob_BC(IREGION_OUTER_CORE)
  write(IOUT,*) 'integer, parameter :: NGLOB_INNER_CORE = ',nglob_AB(IREGION_INNER_CORE)
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPECMAX_CRUST_MANTLE = NSPEC_CRUST_MANTLE_BC'
  write(IOUT,*) 'integer, parameter :: NGLOBMAX_CRUST_MANTLE = NGLOB_CRUST_MANTLE_BC'
  write(IOUT,*)
  write(IOUT,*) 'integer, parameter :: NSPECMAX_OUTER_CORE = NSPEC_OUTER_CORE_BC'
  write(IOUT,*) 'integer, parameter :: NGLOBMAX_OUTER_CORE = NGLOB_OUTER_CORE_BC'
  write(IOUT,*)

  if(ANISOTROPIC_INNER_CORE) then
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_IC = ',NSPEC_AB(IREGION_INNER_CORE)
  else
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_IC = 1'
  endif

  if(ANISOTROPIC_MANTLE) then
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ISO_MANTLE = ',1
    write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',1
    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_MANTLE = NSPECMAX_CRUST_MANTLE'
  else

    write(IOUT,*) 'integer, parameter :: NSPECMAX_ISO_MANTLE = NSPECMAX_CRUST_MANTLE'
    if(TRANSVERSE_ISOTROPY) then
      write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',nspec_aniso_mantle
    else
      write(IOUT,*) 'integer, parameter :: NSPECMAX_TISO_MANTLE = ',1
    endif

    write(IOUT,*) 'integer, parameter :: NSPECMAX_ANISO_MANTLE = 1'
  endif

  write(IOUT,*)

! if attenuation is off, set dummy size of arrays to one
  if(ATTENUATION) then
    write(IOUT,*) 'integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = NSPECMAX_CRUST_MANTLE'
    write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = NSPEC_INNER_CORE'
  else
    write(IOUT,*) 'integer, parameter :: NSPECMAX_CRUST_MANTLE_ATTENUAT = 1'
    write(IOUT,*) 'integer, parameter :: NSPEC_INNER_CORE_ATTENUATION = 1'
  endif

! this to allow for code elimination by compiler in solver for performance
  write(IOUT,*)

  if(TRANSVERSE_ISOTROPY) then
    write(IOUT,*) 'logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: TRANSVERSE_ISOTROPY_VAL = .false.'
  endif
  write(IOUT,*)

  if(ANISOTROPIC_MANTLE) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_MANTLE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_MANTLE_VAL = .false.'
  endif
  write(IOUT,*)

  if(ANISOTROPIC_INNER_CORE) then
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ANISOTROPIC_INNER_CORE_VAL = .false.'
  endif
  write(IOUT,*)

  if(ATTENUATION) then
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ATTENUATION_VAL = .false.'
  endif
  write(IOUT,*)

  if(ELLIPTICITY) then
    write(IOUT,*) 'logical, parameter :: ELLIPTICITY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: ELLIPTICITY_VAL = .false.'
  endif
  write(IOUT,*)

  if(GRAVITY) then
    write(IOUT,*) 'logical, parameter :: GRAVITY_VAL = .true.'
  else
    write(IOUT,*) 'logical, parameter :: GRAVITY_VAL = .false.'
  endif
  write(IOUT,*)

  if(ROTATION) then
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .true.'
    write(IOUT,*) 'integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = NSPECMAX_OUTER_CORE'
  else
    write(IOUT,*) 'logical, parameter :: ROTATION_VAL = .false.'
    write(IOUT,*) 'integer, parameter :: NSPECMAX_OUTER_CORE_ROTATION = 1'
  endif
  write(IOUT,*)

  close(IOUT)

  end subroutine save_header_file

