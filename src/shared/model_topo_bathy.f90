!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine model_topo_bathy_broadcast(ibathy_topo,LOCAL_PATH)

! standard routine to setup model

  use constants, only: myrank,MAX_STRING_LEN,IMAIN,GRAVITY_INTEGRALS, &
                       PLOT_PNM_IMAGE_TOPO_BATHY

  use shared_parameters, only: NX_BATHY,NY_BATHY,PATHNAME_TOPO_FILE,RESOLUTION_TOPO_FILE

  implicit none

  ! bathymetry and topography: use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo
  character(len=MAX_STRING_LEN) :: LOCAL_PATH

  ! local parameters
  character(len=4) :: ending

  ! timer MPI
  double precision :: time1,tCPU
  double precision, external :: wtime

  ending = ''
  if (len_trim(PATHNAME_TOPO_FILE) > 4) ending = PATHNAME_TOPO_FILE(len_trim(PATHNAME_TOPO_FILE)-3:len_trim(PATHNAME_TOPO_FILE))

  if (myrank == 0) then
    ! user output
    write(IMAIN,*)
    write(IMAIN,*) 'incorporating topography'
    write(IMAIN,*) '  topo file            : ',trim(PATHNAME_TOPO_FILE)
    write(IMAIN,*) '  resolution in minutes: ',sngl(RESOLUTION_TOPO_FILE)
    write(IMAIN,*)
    call flush_IMAIN()

    ! get MPI starting time
    time1 = wtime()

    ! read/save topo file on main proc
    if (ending == '.dat') then
      call read_topo_bathy_file_dat_text(ibathy_topo)
    else
      call read_topo_bathy_file(ibathy_topo)
    endif

    ! elapsed time
    tCPU = wtime() - time1

    ! user output
    write(IMAIN,*) '  topography/bathymetry: min/max = ',minval(ibathy_topo),maxval(ibathy_topo)
    write(IMAIN,*)
    write(IMAIN,*) '  Elapsed time for reading in seconds = ',tCPU
    write(IMAIN,*)
    call flush_IMAIN()

    if (.not. GRAVITY_INTEGRALS) call save_topo_bathy_database(ibathy_topo,LOCAL_PATH)

    ! plots image
    if (PLOT_PNM_IMAGE_TOPO_BATHY) call plot_topo_bathy_pnm(ibathy_topo)
  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_i(ibathy_topo,NX_BATHY*NY_BATHY)

  end subroutine model_topo_bathy_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_topo_bathy_file(ibathy_topo)

! reads topography and bathymetry file (given in binary format, name ending in *.bin)
!
! from Elliot, 2014:
!    This expects a file containing signed 16-bit integers, as generated
!    using the ascii2bin.py script in DATA/topo_bathy/. Byte swapping should
!    be handled automatically if necessary.

  use constants
  use shared_parameters, only: NX_BATHY,NY_BATHY,TOPO_MINIMUM,TOPO_MAXIMUM,PATHNAME_TOPO_FILE

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! local parameters
  integer(kind=8) :: filesize   ! 8-bytes / 64-bits
  integer(kind=2) :: ival       ! 2-bytes / 16-bits
  integer(kind=2),dimension(NX_BATHY) :: ival_array
  integer :: indx,itopo_x,itopo_y,itmp
  logical :: byteswap

  integer(kind=2) :: HEADER_IS_BYTE_SWAPPED
  data HEADER_IS_BYTE_SWAPPED/z'3412'/

  ! reads in topography values from file
  filesize = NX_BATHY * NY_BATHY * 2 + 2
  call open_file_abs_r(10, trim(PATHNAME_TOPO_FILE), len_trim(PATHNAME_TOPO_FILE), filesize)

  ! checks byte ordering
  !
  ! first entry in file is a byte marker
  ! (see file convert_etopo_files_from_specfem_ASCII_to_binary_.**.py in DATA/topo_bathy/ directory)
  indx = 1
  call read_abs(10, ival, 2, indx)
  byteswap = (ival == HEADER_IS_BYTE_SWAPPED)

  !debug
  !print *,'topo bathy: byteswap ',byteswap,ival

  ! reads in topography array
  if (byteswap) then
    ! swapping byte from little- to big-endian or vice verse
    do itopo_y = 1,NY_BATHY
      do itopo_x = 1,NX_BATHY
        indx = indx + 1
        call read_abs(10, ival, 2, indx)

        ! note: ibm's xlf compiler warns about ishftc() with integer(2) input. ival should have type integer.
        !       other compilers would use iishift for integer(2) types.
        !ival = ishftc(ival, 8, 16)
        ! work-around
        itmp = ival
        ival = int(ishftc(itmp, 8, 16),kind=2)

        ! stores in array
        ibathy_topo(itopo_x,itopo_y) = ival

        ! checks values
        if (ival < TOPO_MINIMUM .or. ival > TOPO_MAXIMUM) then
          print *,'Error read topo_bathy: ival = ',ival,'at ix/iy = ',itopo_x,itopo_y,'exceeds min/max topography bounds'
          print *,'topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
          call exit_mpi(0,'Error reading topo_bathy file value exceeds min/max bounds')
        endif
      enddo
    enddo
  else
    ! keeps same byte order
    do itopo_y = 1,NY_BATHY
      ! fast way: reads all values in one direction (shifted by 2 bytes due to first marker entry)
      call read_abs_shifted(10,ival_array,2 * NX_BATHY,itopo_y,2)

      ! stores one-by-one
      do itopo_x = 1,NX_BATHY
        ! slow way: reads 1 entry at a time
        !indx = indx + 1
        !call read_abs(10, ival, 2, indx)
        !
        ! fast way: gets topo value from array
        ival = ival_array(itopo_x)

        ! stores in array
        ibathy_topo(itopo_x,itopo_y) = ival

        ! checks values
        if (ival < TOPO_MINIMUM .or. ival > TOPO_MAXIMUM) then
          print *,'Error read topo_bathy: ival = ',ival,'at ix/iy = ',itopo_x,itopo_y,'exceeds min/max topography bounds'
          print *,'topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
          call exit_mpi(0,'Error reading topo_bathy file value exceeds min/max bounds')
        endif
      enddo
    enddo
  endif ! byteswap

  ! closes file
  call close_file_abs(10)

  ! debug
  !print *,'ibathy_topo min/max = ',minval(ibathy_topo),maxval(ibathy_topo)
  !print *,'ibathy_topo ',ibathy_topo(1:10,5)

  end subroutine read_topo_bathy_file

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_topo_bathy_file_dat_text(ibathy_topo)

! reads topography and bathymetry file (given in ASCII format, name ending in *.dat)
! (older formats used for version 6.0)

  use constants
  use shared_parameters, only: NX_BATHY,NY_BATHY,TOPO_MINIMUM,TOPO_MAXIMUM,PATHNAME_TOPO_FILE

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  ! local parameters
  real :: val
  integer :: ival
  integer :: itopo_x,itopo_y,ier

  ! reads in topography values from file
  open(unit=IIN,file=trim(PATHNAME_TOPO_FILE),status='old',action='read',iostat=ier)
  if ( ier /= 0 ) then
    print *,'Error opening:',trim(PATHNAME_TOPO_FILE)
    call exit_mpi(0,'Error opening topography data file')
  endif

  ! reads in topography array
  do itopo_y = 1,NY_BATHY
    do itopo_x = 1,NX_BATHY
      read(IIN,*,iostat=ier) val

      ! checks
      if ( ier /= 0 ) then
        print *,'error read topo_bathy: ix,iy = ',itopo_x,itopo_y,val
        print *,'topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
        call exit_mpi(0,'error reading topo_bathy file')
      endif

      ! converts to integer
      ival = nint(val)

      ! checks values
      if ( ival < TOPO_MINIMUM .or. ival > TOPO_MAXIMUM ) then
        print *,'Error read topo_bathy: ival = ',ival,val,'ix,iy = ',itopo_x,itopo_y
        print *,' topo_bathy dimension: nx,ny = ',NX_BATHY,NY_BATHY
        call exit_mpi(0,'Error reading topo_bathy file')
      endif

      ! stores in array
      ibathy_topo(itopo_x,itopo_y) = ival

    enddo
  enddo
  close(IIN)

  end subroutine read_topo_bathy_file_dat_text


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_topo_bathy_database(ibathy_topo,LOCAL_PATH)

  use constants
  use shared_parameters, only: NX_BATHY,NY_BATHY

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY),intent(in) :: ibathy_topo
  character(len=MAX_STRING_LEN),intent(in) :: LOCAL_PATH

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer :: ier

  filename = trim(LOCAL_PATH) // '/' // 'mesh_topo_bathy.bin'

  ! saves topography and bathymetry file for solver
  open(unit=IOUT,file=trim(filename), &
        status='unknown',form='unformatted',action='write',iostat=ier)

  if (ier /= 0) then
    ! inform about missing database topo file
    print *,'TOPOGRAPHY problem:'
    print *,'Error opening file: ',trim(filename)
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
  use shared_parameters, only: NX_BATHY,NY_BATHY

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY),intent(inout) :: ibathy_topo
  character(len=MAX_STRING_LEN),intent(in) :: LOCAL_PATH

  ! local parameters
  character(len=MAX_STRING_LEN) :: filename
  integer :: ier

  filename = trim(LOCAL_PATH) // '/' // 'mesh_topo_bathy.bin'

  ! reads topography and bathymetry file from saved database file
  open(unit=IIN,file=trim(filename), &
        status='unknown',form='unformatted',action='read',iostat=ier)

  if (ier /= 0) then
    ! inform user
    print *,'TOPOGRAPHY problem:'
    print *,'Error opening file: ',trim(filename)
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
  if (PLOT_PNM_IMAGE_TOPO_BATHY) call plot_topo_bathy_pnm(ibathy_topo)

  end subroutine read_topo_bathy_database

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_topo_bathy(xlat,xlon,value,ibathy_topo)

!
!---- get elevation or ocean depth in meters at a given latitude and longitude
!

  use constants
  use shared_parameters, only: NX_BATHY,NY_BATHY,RESOLUTION_TOPO_FILE

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

  lon_corner = iel1 * samples_per_degree_topo
  lat_corner = 90.d0 - iadd1 * samples_per_degree_topo

  ratio_lon = (xlo-lon_corner)/samples_per_degree_topo
  ratio_lat = (xlat-lat_corner)/samples_per_degree_topo

  if (ratio_lon < 0.d0) ratio_lon = 0.d0
  if (ratio_lon > 1.d0) ratio_lon = 1.d0
  if (ratio_lat < 0.d0) ratio_lat = 0.d0
  if (ratio_lat > 1.d0) ratio_lat = 1.d0

  ! convert integer value to double precision
  if (iadd1 <= NY_BATHY-1 .and. iel1 <= NX_BATHY-1) then
    ! interpolates for points within boundaries
    value = dble(ibathy_topo(iel1,iadd1))     * (1.d0-ratio_lon) * (1.d0-ratio_lat) &
          + dble(ibathy_topo(iel1+1,iadd1))   * ratio_lon * (1.d0-ratio_lat) &
          + dble(ibathy_topo(iel1+1,iadd1+1)) * ratio_lon * ratio_lat &
          + dble(ibathy_topo(iel1,iadd1+1))   * (1.d0-ratio_lon) * ratio_lat

  else if (iadd1 <= NY_BATHY-1 .and. iel1 == NX_BATHY) then
    ! interpolates for points on longitude border
    value = dble(ibathy_topo(iel1,iadd1))   * (1.d0-ratio_lon)*(1.d0-ratio_lat) &
          + dble(ibathy_topo(1,iadd1))      * ratio_lon*(1.d0-ratio_lat) &
          + dble(ibathy_topo(1,iadd1+1))    * ratio_lon*ratio_lat &
          + dble(ibathy_topo(iel1,iadd1+1)) * (1.d0-ratio_lon)*ratio_lat

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

  use constants, only: IOUT,IMAIN,PLOT_PNM_IMAGE_TOPO_BATHY
  use shared_parameters, only: OUTPUT_FILES,NX_BATHY,NY_BATHY

  implicit none

  ! use integer array to store values
  integer, dimension(NX_BATHY,NY_BATHY),intent(in) :: ibathy_topo

  ! local parameters
  integer :: ix,iy,ival,ier
  integer :: minvalue,maxvalue

  ! checks if anything to do
  if (.not. PLOT_PNM_IMAGE_TOPO_BATHY) return

  ! gets min and max
  minvalue = minval(ibathy_topo)
  maxvalue = maxval(ibathy_topo)

  ! creates the PNM image
  write(IMAIN,*) '  plotting PNM image ',trim(OUTPUT_FILES)//'/'//'image_topo_bathy.pnm'
  write(IMAIN,*)
  call flush_IMAIN()

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

      if (ival < 0) ival = 0
      if (ival > 255) ival = 255

      ! write data value (red = green = blue to produce grey levels)
      write(IOUT,'(i3,1x,i3,1x,i3)') ival,ival,ival
    enddo
  enddo

  close(IOUT)

  end subroutine plot_topo_bathy_pnm

