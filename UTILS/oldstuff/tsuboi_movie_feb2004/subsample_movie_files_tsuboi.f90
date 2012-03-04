
!! DK DK subsample high-resolution OpenDX movies to low-resolution movies for ES in Japan

  program subsample_movie_files_tsuboi

  implicit none

!! DK DK this read manually from one of the DX files (total number of points and elements of high-res movie files)
  integer, parameter :: npoin = 2457602
  integer, parameter :: nelem = 2457600

!! DK DK subsampling from high-res to low-res by taking the corners of the cell
  integer, parameter :: NGLL = 5
  integer, parameter :: NSUBSAMP_FACTOR = (NGLL-1)*(NGLL-1)

  double precision, dimension(npoin) :: x,y,z,dataval
  integer, dimension(nelem) :: ibool1,ibool2,ibool3,ibool4

!! DK DK new numbering
  integer, dimension(npoin) :: inewpointnumber
  logical, dimension(npoin) :: mask_numbered

  integer ipoin,ielem,ielem_base_highres,npoin_lowres,nelem_lowres,ipoin_current

! skip header (split on two lines because of line length) in high-res OpenDX file
  read(*,*)
  read(*,*)

  do ipoin=1,npoin
    read(*,*) x(ipoin),y(ipoin),z(ipoin)
  enddo

! skip header
  read(*,*)

  do ielem=1,nelem
    read(*,*) ibool1(ielem),ibool2(ielem),ibool3(ielem),ibool4(ielem)
  enddo

! skip header
  read(*,*)
  read(*,*)
  read(*,*)

  do ipoin=1,npoin
    read(*,*) dataval(ipoin)
  enddo

!! DK DK write smaller file

! compute number of elements in low-res file
  nelem_lowres = nelem / NSUBSAMP_FACTOR

! mask data points that are used
  inewpointnumber(:) = 0
  mask_numbered(:) = .false.
  ipoin_current = 0
  ielem_base_highres = 0

  do ielem=1,nelem_lowres

!! DK DK point numbers start at 0 in OpenDX, therefore add 1 to index value
    if(.not. mask_numbered(ibool1(ielem_base_highres+1) + 1)) then
      inewpointnumber(ibool1(ielem_base_highres+1) + 1) = ipoin_current
      ipoin_current = ipoin_current + 1
      mask_numbered(ibool1(ielem_base_highres+1) + 1) = .true.
    endif

    if(.not. mask_numbered(ibool2(ielem_base_highres+13) + 1)) then
      inewpointnumber(ibool2(ielem_base_highres+13) + 1) = ipoin_current
      ipoin_current = ipoin_current + 1
      mask_numbered(ibool2(ielem_base_highres+13) + 1) = .true.
    endif

    if(.not. mask_numbered(ibool3(ielem_base_highres+4) + 1)) then
      inewpointnumber(ibool3(ielem_base_highres+4) + 1) = ipoin_current
      ipoin_current = ipoin_current + 1
      mask_numbered(ibool3(ielem_base_highres+4) + 1) = .true.
    endif

    if(.not. mask_numbered(ibool4(ielem_base_highres+16) + 1)) then
      inewpointnumber(ibool4(ielem_base_highres+16) + 1) = ipoin_current
      ipoin_current = ipoin_current + 1
      mask_numbered(ibool4(ielem_base_highres+16) + 1) = .true.
    endif

    ielem_base_highres = ielem_base_highres + NSUBSAMP_FACTOR
  enddo

! compute number of points in low-res file
  npoin_lowres = ipoin_current

! points
  write(*,*) 'object 1 class array type float rank 1 shape 3 items ',npoin_lowres,' data follows'

  do ipoin=1,npoin
    if(mask_numbered(ipoin)) write(*,"(f8.5,1x,f8.5,1x,f8.5)") x(ipoin),y(ipoin),z(ipoin)
  enddo

! elements
  write(*,*) 'object 2 class array type int rank 1 shape 4 items ',nelem_lowres,' data follows'

! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
  ielem_base_highres = 0
  do ielem=1,nelem_lowres
!! DK DK extract the four corners of the NGLL x NGLL high-res OpenDX element
    write(*,210) inewpointnumber(ibool1(ielem_base_highres+1) + 1),inewpointnumber(ibool2(ielem_base_highres+13) + 1), &
                 inewpointnumber(ibool3(ielem_base_highres+4) + 1),inewpointnumber(ibool4(ielem_base_highres+16) + 1)
    ielem_base_highres = ielem_base_highres + NSUBSAMP_FACTOR
  enddo

! data values
  write(*,*) 'attribute "element type" string "quads"'
  write(*,*) 'attribute "ref" string "positions"'
  write(*,*) 'object 3 class array type float rank 0 items ',npoin_lowres,' data follows'
  do ipoin=1,npoin
    if(mask_numbered(ipoin)) write(*,501) dataval(ipoin)
  enddo

  write(*,*) 'attribute "dep" string "positions"'
  write(*,*) 'object "irregular positions irregular connections" class field'
  write(*,*) 'component "positions" value 1'
  write(*,*) 'component "connections" value 2'
  write(*,*) 'component "data" value 3'
  write(*,*) 'end'

 210 format(i10,1x,i10,1x,i10,1x,i10)
 501 format(f7.2)

 end program subsample_movie_files_tsuboi

