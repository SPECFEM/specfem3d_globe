!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
!          --------------------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

!--------------------------------------------------------------------------------------------------
! ETOPO
!
! Global Gridded Elevation Data
!
! by default (constants.h), it uses a smoothed ETOPO 4 dataset
!--------------------------------------------------------------------------------------------------

  subroutine model_topo_bathy_broadcast(myrank,ibathy_topo,LOCAL_PATH)

! standard routine to setup model

  use constants

  implicit none

  integer :: myrank

  ! bathymetry and topography: use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  if (myrank == 0) then
    ! user output
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating topography'
    call flush_IMAIN()

    ! read/save topo file on master
    call read_topo_bathy_file(ibathy_topo)

    if (.not. GRAVITY_INTEGRALS) call save_topo_bathy_database(ibathy_topo,LOCAL_PATH)
  endif

  ! broadcast the information read on the master to the nodes
  call bcast_all_i(ibathy_topo,NX_BATHY*NY_BATHY)

  end subroutine model_topo_bathy_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_topo_bathy_file(ibathy_topo)
!
!---- read topography and bathymetry file once and for all
!
  use constants

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! local parameters
  integer(kind=8) :: filesize
  integer(kind=2) :: ival
  integer :: indx,itopo_x,itopo_y
  logical :: byteswap
  integer(kind=2) :: HEADER_IS_BYTE_SWAPPED
  data HEADER_IS_BYTE_SWAPPED/z'3412'/

  ! reads in topography values from file
  filesize = NX_BATHY * NY_BATHY * 2 + 2
  call open_file_abs_r(10, trim(PATHNAME_TOPO_FILE), len_trim(PATHNAME_TOPO_FILE), filesize)

  ! checks byte ordering
  indx = 1
  call read_abs(10, ival, 2, indx)
  byteswap = (ival == HEADER_IS_BYTE_SWAPPED)

  ! reads in topography array
  do itopo_y = 1,NY_BATHY
    do itopo_x = 1,NX_BATHY
      indx = indx + 1
      call read_abs(10, ival, 2, indx)
      if (byteswap) then
        ival = ishftc(ival, 8, 16)
      endif

      ! checks values
      if (ival < TOPO_MINIMUM .or. ival > TOPO_MAXIMUM) then
        print *,'Error read topo_bathy: ival = ',ival,'ix,iy = ',itopo_x,itopo_y
        print *,'topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
        call exit_mpi(0,'Error reading topo_bathy file')
      endif

      ! stores in array
      ibathy_topo(itopo_x,itopo_y) = ival
    enddo
  enddo
  call close_file_abs(10)

  ! user output
  write(IMAIN,*) "  topography/bathymetry: min/max = ",minval(ibathy_topo),maxval(ibathy_topo)

  ! plots image
  call plot_topo_bathy_pnm(ibathy_topo)

  end subroutine read_topo_bathy_file

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_topo_bathy_database(ibathy_topo,LOCAL_PATH)

  use constants

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  character(len=MAX_STRING_LEN) :: prname
  integer :: ier

  ! create the name for the database of the current slide and region
  ! only master needs to save this
  call create_name_database(prname,0,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! saves topography and bathymetry file for solver
  open(unit=IOUT,file=prname(1:len_trim(prname))//'topo.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)

  if (ier /= 0) then
    ! inform about missing database topo file
    print *,'TOPOGRAPHY problem:'
    print *,'Error opening file: ',prname(1:len_trim(prname))//'topo.bin'
    print *,'please check if path exists and rerun mesher'
    call exit_mpi(0,'Error opening file for database topo')
  endif

  write(IOUT) ibathy_topo
  close(IOUT)

  end subroutine save_topo_bathy_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_topo_bathy_database(ibathy_topo,LOCAL_PATH)

  use constants

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  character(len=MAX_STRING_LEN) :: prname
  integer :: ier

  ! create the name for the database of the current slide and region
  ! only master needs to save this
  call create_name_database(prname,0,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! reads topography and bathymetry file from saved database file
  open(unit=IIN,file=prname(1:len_trim(prname))//'topo.bin', &
        status='unknown',form='unformatted',action='read',iostat=ier)

  if (ier /= 0) then
    ! inform user
    print *,'TOPOGRAPHY problem:'
    print *,'Error opening file: ',prname(1:len_trim(prname))//'topo.bin'
    !print *,'please check if file exists and rerun solver'
    !call exit_mpi(0,'Error opening file for database topo')

    ! read by original file
    print *,'trying original topography file...'
    call read_topo_bathy_file(ibathy_topo)

    ! saves database topo file for next time
    if (.not. GRAVITY_INTEGRALS) call save_topo_bathy_database(ibathy_topo,LOCAL_PATH)

  else
    ! database topo file exists
    read(IIN) ibathy_topo
    close(IIN)

    ! user output
    write(IMAIN,*) "  topography/bathymetry: min/max = ",minval(ibathy_topo),maxval(ibathy_topo)
    call flush_IMAIN()

  endif

  ! plots image
  call plot_topo_bathy_pnm(ibathy_topo)

  end subroutine read_topo_bathy_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_bathy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

  use constants

  implicit none

  ! location latitude/longitude (in degree)
  double precision,intent(in):: xlat,xlon

  ! returns elevation (in meters)
  double precision,intent(out):: value

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY),intent(in) :: ibathy_topo

  ! local parameters
  integer :: iadd1,iel1
  double precision :: samples_per_degree_topo
  double precision :: xlo
  double precision :: lon_corner,lat_corner,ratio_lon,ratio_lat

  ! initializes elevation
  value = ZERO

  ! longitude within range [0,360] degrees
  xlo = xlon
  if (xlo < 0.d0) xlo = xlo + 360.d0
  if (xlo > 360.d0) xlo = xlo - 360.d0

  ! compute number of samples per degree
  samples_per_degree_topo = dble(RESOLUTION_TOPO_FILE) / 60.d0

  ! compute offset in data file and avoid edge effects
  iadd1 = 1 + int((90.d0-xlat)/samples_per_degree_topo)
  if (iadd1 < 1) iadd1 = 1
  if (iadd1 > NY_BATHY) iadd1 = NY_BATHY

  iel1 = int(xlo/samples_per_degree_topo)
  if (iel1 <= 0 .or. iel1 > NX_BATHY) iel1 = NX_BATHY

! Use bilinear interpolation rather nearest point interpolation

  ! convert integer value to double precision
  !  value = dble(ibathy_topo(iel1,iadd1))

  lon_corner=iel1*samples_per_degree_topo
  lat_corner=90.d0-iadd1*samples_per_degree_topo

  ratio_lon = (xlo-lon_corner)/samples_per_degree_topo
  ratio_lat = (xlat-lat_corner)/samples_per_degree_topo

  if (ratio_lon<0.0) ratio_lon=0.0
  if (ratio_lon>1.0) ratio_lon=1.0
  if (ratio_lat<0.0) ratio_lat=0.0
  if (ratio_lat>1.0) ratio_lat=1.0

  ! convert integer value to double precision
  if (iadd1 <= NY_BATHY-1 .and. iel1 <= NX_BATHY-1) then
    ! interpolates for points within boundaries
    value = dble(ibathy_topo(iel1,iadd1))*(1-ratio_lon)*(1.-ratio_lat) &
            + dble(ibathy_topo(iel1+1,iadd1))*ratio_lon*(1.-ratio_lat) &
            + dble(ibathy_topo(iel1+1,iadd1+1))*ratio_lon*ratio_lat &
            + dble(ibathy_topo(iel1,iadd1+1))*(1.-ratio_lon)*ratio_lat

  else if (iadd1 <= NY_BATHY-1 .and. iel1 == NX_BATHY) then
    ! interpolates for points on longitude border
    value = dble(ibathy_topo(iel1,iadd1))*(1-ratio_lon)*(1.-ratio_lat) &
            + dble(ibathy_topo(1,iadd1))*ratio_lon*(1.-ratio_lat) &
            + dble(ibathy_topo(1,iadd1+1))*ratio_lon*ratio_lat &
            + dble(ibathy_topo(iel1,iadd1+1))*(1.-ratio_lon)*ratio_lat

  else
    ! for points on latitude boundaries
    value = dble(ibathy_topo(iel1,iadd1))
  endif

  end subroutine get_topo_bathy

!
!-------------------------------------------------------------------------------------------------
!

  subroutine plot_topo_bathy_pnm(ibathy_topo)

! stores topo_bathy image in PNM format with grey levels

  use constants,only: NX_BATHY,NY_BATHY,IOUT,IMAIN
  use shared_input_parameters, only: OUTPUT_FILES

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY),intent(in) :: ibathy_topo

  ! local parameters
  integer :: ix,iy,ival,ier
  integer :: minvalue,maxvalue

  !----------------------------------------------------------------------

  ! for debugging: plots pnm-image showing used topography
  !                file can become fairly big for large topo-files, e.g. ETOPO1 creates a ~2.7 GB pnm-image
  logical,parameter :: DO_IMAGE_PLOT = .false.

  !----------------------------------------------------------------------

  ! checks if anything to do
  if (.not. DO_IMAGE_PLOT) return

  ! gets min and max
  minvalue = minval(ibathy_topo)
  maxvalue = maxval(ibathy_topo)

  ! creates the PNM image
  write(IMAIN,*) '  plotting PNM image ',trim(OUTPUT_FILES)//'/'//'image_topo_bathy.pnm'
  write(IMAIN,*)

  ! creating the header
  open(unit=IOUT,file=trim(OUTPUT_FILES)//'/'//'image_topo_bathy.pnm',status='unknown',iostat=ier)
  if (ier /= 0) stop 'Error opening file image_topo_bathy.pnm'

  write(IOUT,'(a)') 'P3'
  write(IOUT,'(i6,1x,i6)') NX_BATHY,NY_BATHY
  write(IOUT,'(i3)') 255

  ! creates image with grey levels
  do iy = 1,NY_BATHY
    do ix = 1,NX_BATHY
      if (minvalue == maxvalue) then
        ival = 128
      else
        ival = 255 * (ibathy_topo(ix,iy) - minvalue) / (maxvalue - minvalue)
      endif

      if(ival < 1) ival = 1
      if(ival > 255) ival = 255

      ! write data value (red = green = blue to produce grey levels)
      write(IOUT,'(i3)') ival
      write(IOUT,'(i3)') ival
      write(IOUT,'(i3)') ival
    enddo
  enddo

  close(IOUT)

  end subroutine plot_topo_bathy_pnm

