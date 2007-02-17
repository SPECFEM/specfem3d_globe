
! extract only the points that are needed for the new doubling brick
! from the modified version of Emanuele's tripling brick

  program convert_OpenDX_SPECFEM_format

  implicit none

! number of points in the superbrick
  integer, parameter :: NPOINTS = 67

! number of spectral elements in the superbrick
  integer, parameter :: NSPEC = 32

! bricks are made of eight points
  integer, parameter :: NGNOD = 8

  integer :: ispec,ipoin,ignod

! coordinates of all the mesh points
  double precision, dimension(NPOINTS) :: x,y,z

! point numbers of all the nodes of each spectral element
  integer, dimension(NGNOD) :: iglob

! read the basic brick
  open(unit=27,file='dimitri_superbrick_without_duplicates.dx',status='old',action='read')

! skip header
  read(27,*)

  print *
  print *,' subroutine define_superbrick(x_superbrick,y_superbrick,z_superbrick,ibool_superbrick)'
  print *
  print *,' implicit none'
  print *
  print *,' include "constants_modified.h"'
  print *
  print *,' integer, dimension(NGNOD_DOUBLING_SUPERBRICK,NSPEC_DOUBLING_SUPERBRICK) :: ibool_superbrick'
  print *,' double precision, dimension(NGLOB_DOUBLING_SUPERBRICK) :: x_superbrick,y_superbrick,z_superbrick'
  print *

! read the points and output their coordinates in Fortran90 format
! output all the points, without duplicates
  do ipoin = 1,NPOINTS

    read(27,*) x(ipoin),y(ipoin),z(ipoin)

! write in a nice format with full precision
    if(abs(x(ipoin) - 0.d0) < 0.0001d0) then
      print *,' x_superbrick(',ipoin,') = 0.d0'
    else if(abs(x(ipoin) - 0.5d0) < 0.0001d0) then
      print *,' x_superbrick(',ipoin,') = 1.d0 / 2.d0'
    else if(abs(x(ipoin) - 1.d0) < 0.0001d0) then
      print *,' x_superbrick(',ipoin,') = 1.d0'
    else if(abs(x(ipoin) - 1.5d0) < 0.0001d0) then
      print *,' x_superbrick(',ipoin,') = 3.d0 / 2.d0'
    else if(abs(x(ipoin) - 2.d0) < 0.0001d0) then
      print *,' x_superbrick(',ipoin,') = 2.d0'
    else
      stop 'unknown coordinate found for x'
    endif

    if(abs(y(ipoin) - 0.d0) < 0.0001d0) then
      print *,' y_superbrick(',ipoin,') = 0.d0'
    else if(abs(y(ipoin) - 0.5d0) < 0.0001d0) then
      print *,' y_superbrick(',ipoin,') = 1.d0 / 2.d0'
    else if(abs(y(ipoin) - 1.d0) < 0.0001d0) then
      print *,' y_superbrick(',ipoin,') = 1.d0'
    else if(abs(y(ipoin) - 1.5d0) < 0.0001d0) then
      print *,' y_superbrick(',ipoin,') = 3.d0 / 2.d0'
    else if(abs(y(ipoin) - 2.d0) < 0.0001d0) then
      print *,' y_superbrick(',ipoin,') = 2.d0'
    else
      stop 'unknown coordinate found for y'
    endif

    if(abs(z(ipoin) - 0.d0) < 0.0001d0) then
      print *,' z_superbrick(',ipoin,') = 0.d0'
    else if(abs(z(ipoin) - 0.5d0) < 0.0001d0) then
      print *,' z_superbrick(',ipoin,') = 1.d0 / 2.d0'
    else if(abs(z(ipoin) - 1.d0) < 0.0001d0) then
      print *,' z_superbrick(',ipoin,') = 1.d0'
    else if(abs(z(ipoin) - 1.5d0) < 0.0001d0) then
      print *,' z_superbrick(',ipoin,') = 3.d0 / 2.d0'
    else if(abs(z(ipoin) - 2.d0) < 0.0001d0) then
      print *,' z_superbrick(',ipoin,') = 2.d0'
    else
      stop 'unknown coordinate found for z'
    endif

    print *

  enddo

! skip header
  read(27,*)

! read the points of each of the NSPEC elements
  do ispec = 1, NSPEC

! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in SPECFEM3D, therefore swap the points
    read(27,"(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
            iglob(4),iglob(1),iglob(8),iglob(5),iglob(3),iglob(2),iglob(7),iglob(6)

    do ignod = 1,NGNOD
! in the case of OpenDX, node numbers start at zero
      print *,' ibool_superbrick(',ignod,',',ispec,') = ',iglob(ignod) + 1
    enddo
    print *

  enddo

  close(27)

  print *,' end subroutine define_superbrick'
  print *

  end program convert_OpenDX_SPECFEM_format

