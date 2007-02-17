
! extract only the points that are needed for the new doubling brick
! from the modified version of Emanuele's tripling brick

  program create_superbrick_no_duplicates

  implicit none

! number of points in the superbrick
  integer, parameter :: OLD_NPOINTS = 108, NEW_NPOINTS = 67

! number of spectral elements in the superbrick
  integer, parameter :: NSPEC = 32

! bricks are made of eight points
  integer, parameter :: NGNOD = 8

  integer :: ispec,ipoin,ignod,ipoin_found,ipoin_previous
  double precision :: distance
  logical :: already_found_before

! coordinates of all the mesh points
  double precision, dimension(OLD_NPOINTS) :: x,y,z
  double precision, dimension(NEW_NPOINTS) :: xnew,ynew,znew

! point numbers of all the nodes of each spectral element
  integer, dimension(NGNOD) :: point_number

! read the basic brick
  open(unit=27,file='dimitri_superbrick_with_duplicates.dx',status='old',action='read')

! skip header
  read(27,*)

! create an OpenDX file with the new superbrick without duplicates

! write OpenDX header with element data
  open(unit=11,file='dimitri_superbrick_without_duplicates.dx',status='unknown')
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',NEW_NPOINTS,' data follows'

! read the points, with duplicates
! output all the points, without duplicates
  ipoin_found = 0
  do ipoin = 1,OLD_NPOINTS

    read(27,*) x(ipoin),y(ipoin),z(ipoin)

! check if this point has previously been found, print only if not
    already_found_before = .false.
    if(ipoin > 1) then
      do ipoin_previous = 1,ipoin - 1
        distance = dsqrt((x(ipoin) - x(ipoin_previous))**2 + (y(ipoin) - y(ipoin_previous))**2 + &
                         (z(ipoin) - z(ipoin_previous))**2)
! it is the same point as before if distance with a previous point is zero (down to roundoff error)
        if(distance < 0.00001d0) then
          already_found_before = .true.
! we can exit the loop if we have already found this point in the list
          exit
        endif
      enddo
    endif
    if(.not. already_found_before) then
      ipoin_found = ipoin_found + 1
      xnew(ipoin_found) = x(ipoin)
      ynew(ipoin_found) = y(ipoin)
      znew(ipoin_found) = z(ipoin)
      write(11,*) sngl(xnew(ipoin_found)),sngl(ynew(ipoin_found)),sngl(znew(ipoin_found))
    endif

  enddo

  if(ipoin_found /= NEW_NPOINTS) then
    print *,'ipoin_found,OLD_NPOINTS,NEW_NPOINTS = ',ipoin_found,OLD_NPOINTS,NEW_NPOINTS
    stop 'incorrect number of unique points found'
  endif

! skip header
  read(27,*)

! write element header
  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',NSPEC,' data follows'

! read the NGNOD points of each of the NSPEC elements
  do ispec = 1, NSPEC

    read(27,*) (point_number(ignod), ignod = 1, NGNOD)

! change old point numbers to new
    do ignod = 1, NGNOD

! in the case of OpenDX, node numbers start at zero
      ipoin = point_number(ignod) + 1

      already_found_before = .false.
      do ipoin_found = 1,NEW_NPOINTS
        distance = dsqrt((x(ipoin) - xnew(ipoin_found))**2 + (y(ipoin) - ynew(ipoin_found))**2 + &
                         (z(ipoin) - znew(ipoin_found))**2)
! it is the same point as before if distance with a previous point is zero (down to roundoff error)
        if(distance < 0.00001d0) then
          already_found_before = .true.
! we can exit the loop if we have already found this point in the list
          exit
        endif
      enddo

! problem if this point is never found
      if(.not. already_found_before) stop 'problem, point not found'

! in the case of OpenDX, node numbers start at zero
      point_number(ignod) = ipoin_found - 1
      if(point_number(ignod) < 0 .or. point_number(ignod) > NEW_NPOINTS - 1) stop 'incorrect point number found'

    enddo

! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
    write(11,"(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") (point_number(ignod), ignod = 1, NGNOD)

  enddo

  close(27)

! output OpenDX header for data
! label for hexahedra in OpenDX is "cubes"
  write(11,*) 'attribute "element type" string "cubes"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',NSPEC,' data follows'

! write element data (use a different color for each element)
  do ispec=1,NSPEC
    write(11,*) ispec
  enddo

! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

  close(11)

  end program create_superbrick_no_duplicates

