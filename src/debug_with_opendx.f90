
  subroutine debug_with_opendx(is_on_a_slice_edge_crust_mantle,xstore_crust_mantle, &
                              ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle)

  implicit none

  include "constants.h"

!! DK DK for the merged version
! include values created by the mesher
  include "values_from_mesher.h"

  logical, dimension(NSPEC_CRUST_MANTLE) :: is_on_a_slice_edge_crust_mantle

! mesh parameters
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  integer :: ntotspecAVS_DX,ioffset,ntotpoinAVS_DX,i,j,k,ispec,iglob
  real :: xval,yval,zval

!! DK DK create an OpenDX file showing all the elements computed in the first step
!! DK DK of non blocking overlapping

 ntotspecAVS_DX = count(is_on_a_slice_edge_crust_mantle(:))
 ntotpoinAVS_DX = 8 * ntotspecAVS_DX

    open(unit=11,file='DX_fullmesh.dx',status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ntotpoinAVS_DX,' data follows'
  do ispec=1,nspec_crust_mantle
if(is_on_a_slice_edge_crust_mantle(ispec)) then

i = 1
j = 1
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = 1
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = NGLLX
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = 1
j = NGLLX
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = 1
j = 1
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = 1
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = NGLLX
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = 1
j = NGLLX
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

endif
  enddo

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',ntotspecAVS_DX,' data follows'

! read local elements in this slice and output global AVS or DX elements
ioffset = 0
  do ispec=1,nspec_crust_mantle
if(is_on_a_slice_edge_crust_mantle(ispec)) then
! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
write(11,"(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
            ioffset+3,ioffset+0,ioffset+7,ioffset+4,ioffset+2,ioffset+1,ioffset+6,ioffset+5
ioffset = ioffset + 8
endif
  enddo


! ************* generate element data values ******************

! output AVS or DX header for data
    write(11,*) 'attribute "element type" string "cubes"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec_crust_mantle
if(is_on_a_slice_edge_crust_mantle(ispec)) then
        write(11,*) '100'
endif
  enddo

! define OpenDX field
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'

  close(11)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! DK DK create an OpenDX file showing all the elements that seems to have a problem
!! DK DK in the second step of non blocking overlapping

!! DK DK suppressed for now
 goto 777

 ntotspecAVS_DX = 2
 ntotpoinAVS_DX = 8 * ntotspecAVS_DX

    open(unit=11,file='DX_fullmesh2.dx',status='unknown')
    write(11,*) 'object 1 class array type float rank 1 shape 3 items ',ntotpoinAVS_DX,' data follows'
  do ispec=1,nspec_crust_mantle
if(ispec == 201 .or. ispec == 203) then

i = 1
j = 1
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = 1
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = NGLLX
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = 1
j = NGLLX
k = 1
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = 1
j = 1
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = 1
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = NGLLX
j = NGLLX
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

i = 1
j = NGLLX
k = NGLLZ
iglob = ibool_crust_mantle(i,j,k,ispec)
call rthetaphi_2_xyz(xval,yval,zval,xstore_crust_mantle(iglob), &
    ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))
        write(11,"(f10.7,1x,f10.7,1x,f10.7)") xval,yval,zval

endif
  enddo

! ************* generate elements ******************

  write(11,*) 'object 2 class array type int rank 1 shape 8 items ',ntotspecAVS_DX,' data follows'

! read local elements in this slice and output global AVS or DX elements
ioffset = 0
  do ispec=1,nspec_crust_mantle
if(ispec == 201 .or. ispec == 203) then
! point order in OpenDX is 4,1,8,5,3,2,7,6, *not* 1,2,3,4,5,6,7,8 as in AVS
! in the case of OpenDX, node numbers start at zero
write(11,"(i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6,1x,i6)") &
            ioffset+3,ioffset+0,ioffset+7,ioffset+4,ioffset+2,ioffset+1,ioffset+6,ioffset+5
ioffset = ioffset + 8
endif
  enddo


! ************* generate element data values ******************

! output AVS or DX header for data
    write(11,*) 'attribute "element type" string "cubes"'
    write(11,*) 'attribute "ref" string "positions"'
    write(11,*) 'object 3 class array type float rank 0 items ',ntotspecAVS_DX,' data follows'

! read local elements in this slice and output global AVS or DX elements
  do ispec=1,nspec_crust_mantle
if(ispec == 201 .or. ispec == 203) then
        write(11,*) '200'
endif
  enddo

! define OpenDX field
    write(11,*) 'attribute "dep" string "connections"'
    write(11,*) 'object "irregular positions irregular connections" class field'
    write(11,*) 'component "positions" value 1'
    write(11,*) 'component "connections" value 2'
    write(11,*) 'component "data" value 3'
    write(11,*) 'end'

  close(11)

  777 continue

  end subroutine debug_with_opendx

