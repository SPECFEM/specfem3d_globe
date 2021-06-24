!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
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

! Block Mantle Spherical Harmonics (BKMNS) models
!
! based on model format by Caio H. Ciardelli
!
!-----------------------
! Model parametrization
!-----------------------
!
! This parametrization includes Vp (km/s), Vs (km/s), Rho (g/cm^3), Vpv (km/s), Vph (km/s), Vsv (km/s), Vsh (km/s)
! and Eta (km/s).
!
! The model is separated in four zones:
!
! Zone 1: Represents the crust using a block model. The horizontal resolution is 0.25 degrees both in latitude
!         and longitude. The vertical resolution is 1 km. This zone stretches from 6 km above sea level to 80 km depth.
!
! Zone 2: Represents the upper mantle using a combination of spherical harmonics up to degree 100 and 6 cubic
!         b-splines. This zone stretches from 80 km to 410 km depth.
!
! Zone 3: Represents the transition zone using a combination of spherical harmonics up to degree 90 and 5 cubic
!         b-splines. This zone stretches from 410 km to 650 km depth.
!
! Zone 4: Represents the lower mantle using a combination of spherical harmonics up to degree 80 and 14 cubic
!         b-splines. This zone stretches from 650 km to 2891 km depth.
!
! This model format combines a block model at the top, between +6km to 80km depth, and the spherical harmonics expansion for
! the mantle below 80km to CMB. It is based on a spherical Earth, but with topography/bathymetry and moho variations.
!
! The block model includes topography/bathymetry. For points above the topo/bathy surface, velocities are set to zero.
! It also includes Moho variations, based on Crust2.0 used for the original GLAD models, thus includes crustal
! and mantle velocities for the actual point positions, but without ellipticity. Note the maximum Moho depth is around 75km,
! therefore the block model expands down to 80km depth to capture its variations.
!
! The mantle expansion into spherical harmonics relies on a spherical mesh, without topography
! of the 410/660/CMB discontinuities.
!
! The model parameterization uses mainly transverse isotropy, i.e., (vpv,vph,vsv,vsh,eta,rho) values, for both block and
! spherical harmonics parts. (and provides additionally derived vp and vs velocities based on a Voigt average).
!


  module model_bkmns_par

  implicit none

  !! constants
  !integer,parameter :: MAX_STRING_LEN = 200         ! Buffer size for file names

  integer,parameter :: NPT            = 100000      ! Number of interpolation points in the radial basis

  ! total number of blocks: N = NP_b * NT_b * NR_b
  integer,parameter :: NR_b           = 87          ! Number of points of the block model in the radial direction
  integer,parameter :: NT_b           = 721         ! Number of points of the block model in the polar direction
  integer,parameter :: NP_b           = 1441        ! Number of points of the block model in the azimuthal direction

  double precision,parameter :: ANGLE_TOLERANCE = 1d-10  ! Water level to convert from Cartesian to spherical coordinates
  double precision,parameter :: WATER_LEVEL     = 1d-15  ! Water level to prevent divisions by zero

  ! Pi and factors to convert radians to degrees and vice-versa
  !double precision,parameter :: PI              = 3.14159265358979323846d0
  !double precision,parameter :: R2D             = 180.d0 / 3.14159265358979323846d0
  !double precision,parameter :: D2R             = 3.14159265358979323846d0 / 180.d0
  !double precision,parameter :: EARTH_R         = 6371.d0     ! Earth radius for models

  double precision,parameter :: NORM_CMB_R      = 0.54622508  ! Normalized core-mantle boundary radius
  double precision,parameter :: NORM_MOHO_R     = 0.98744310  ! Normalized Moho radius
  double precision,parameter :: NORM_MAX_H_R    = 1.00094177  ! Normalized maximum distance from the Earth's centre
                                                              ! to the highest point on the surface of the mesh

  ! K1, K2 and K3 are constantes used to avoid underflow on the Legendre associated funtions
  double precision,parameter :: LEGENDRE_K1              = 1d140
  double precision,parameter :: LEGENDRE_K2              = 1d280
  double precision,parameter :: LEGENDRE_K3              = 1d-140

  ! number of parameters
  integer :: BKMNS_NPAR = 6  ! tiso: vpv,vph,vsv,vsh,eta,rho

  ! block model arrays (crust/upper mantle)
  real,dimension(:,:,:),allocatable, target :: BKMNS_crust_rho,BKMNS_crust_eta, &
    BKMNS_crust_vsh,BKMNS_crust_vsv, &
    BKMNS_crust_vph,BKMNS_crust_vpv

  ! values at surface
  real,dimension(:,:,:), allocatable :: BKMNS_crust_surface_value

  ! mantle zones 2,3,4
  ! radius range
  double precision :: r2min,r2max,r3min,r3max,r4min,r4max

  ! note: using double precision for mantle sh arrays since these are smaller arrays

  ! radial splines positions
  integer :: num_spline_positions_zone2,num_spline_positions_zone3,num_spline_positions_zone4
  double precision,dimension(:), allocatable, target :: T_zone2,T_zone3,T_zone4

  ! spherical harmonics degree
  integer :: degree_N_zone2,degree_N_zone3,degree_N_zone4
  ! spline degrees
  integer :: num_spline_NS_zone2,num_spline_NS_zone3,num_spline_NS_zone4
  double precision, dimension(:,:), allocatable, target :: Rbasis_zone2,Rbasis_zone3,Rbasis_zone4

  ! factors for Legendre polynomials
  integer :: degree_NMAX
  double precision, dimension(:), allocatable :: Legendre_degree_factor

  ! coefficients
  integer :: num_sh_coeffs_zone2,num_sh_coeffs_zone3,num_sh_coeffs_zone4
  double precision,dimension(:,:,:),allocatable, target :: A_zone2, B_zone2, A_zone3, B_zone3, A_zone4, B_zone4
  double precision, dimension(:), allocatable, target :: nFactors_zone2,nFactors_zone3,nFactors_zone4

  end module model_bkmns_par


!
!-----------------------------------------------------------------------------------------
!
! model setup
!
!-----------------------------------------------------------------------------------------
!

  subroutine model_bkmns_crust_broadcast()

! standard routine to setup model

  use constants, only: IMAIN,myrank
  use model_bkmns_par

  implicit none

  ! local parameters
  integer :: ier

  ! number of blocks N = NP_b * NT_b * NR_b
  !                    = 1441 * 721  * 87
  !                    = 90,389,607
  !
  ! for real == 4 bytes values:
  ! array size in MB: N * real / 1024. / 1024. = 344 MB
  !
  ! total size for vsh,..,rho around 6 * 344 MB ~ 2 GB
  ! which should fit memory per core around 4-8 GB

  ! allocate memory
  allocate (BKMNS_crust_vsh(NP_b,NT_b,NR_b), &
            BKMNS_crust_vsv(NP_b,NT_b,NR_b), &
            BKMNS_crust_vph(NP_b,NT_b,NR_b), &
            BKMNS_crust_vpv(NP_b,NT_b,NR_b), &
            BKMNS_crust_eta(NP_b,NT_b,NR_b), &
            BKMNS_crust_rho(NP_b,NT_b,NR_b), &
            BKMNS_crust_surface_value(NP_b,NT_b,BKMNS_NPAR), stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating crust arrays')
  BKMNS_crust_rho(:,:,:) = 0.0; BKMNS_crust_eta(:,:,:) = 0.0
  BKMNS_crust_vph(:,:,:) = 0.0; BKMNS_crust_vpv(:,:,:) = 0.0
  BKMNS_crust_vsh(:,:,:) = 0.0; BKMNS_crust_vsv(:,:,:) = 0.0
  BKMNS_crust_surface_value(:,:,:) = 0.0

  ! the variables read are declared and stored in crust/ subfolder
  if (myrank == 0) call read_crust_bkmns_model()

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_r(BKMNS_crust_vsh,NP_b * NT_b * NR_b)
  call bcast_all_r(BKMNS_crust_vsv,NP_b * NT_b * NR_b)
  call bcast_all_r(BKMNS_crust_vph,NP_b * NT_b * NR_b)
  call bcast_all_r(BKMNS_crust_vpv,NP_b * NT_b * NR_b)
  call bcast_all_r(BKMNS_crust_eta,NP_b * NT_b * NR_b)
  call bcast_all_r(BKMNS_crust_rho,NP_b * NT_b * NR_b)
  call bcast_all_r(BKMNS_crust_surface_value,NP_b * NT_b * BKMNS_NPAR)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done bkmns crust'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine model_bkmns_crust_broadcast

!
!-----------------------------------------------------------------------------------------
!

  subroutine read_crust_bkmns_model()

  use constants, only: myrank,IIN,IMAIN,MAX_STRING_LEN,HUGEVAL,EARTH_R_KM
  use model_bkmns_par

  implicit none

  ! local parameters
  integer :: ipar,ier
  character(len=MAX_STRING_LEN) :: filename
  character(len=MAX_STRING_LEN) :: rootdir
  real,dimension(:,:),allocatable :: values_at_surface

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: bkmns '
  write(IMAIN,*) '  number of blocks: radial     NR_b = ',NR_b
  write(IMAIN,*) '                    latitudes  NT_b = ',NT_b
  write(IMAIN,*) '                    longitudes NP_b = ',NP_b
  write(IMAIN,*)
  call flush_IMAIN()

  ! user output
  write(IMAIN,*) '  reading crustal model from DATA/gladm15/crust/ folder:'
  ! min/max depths of the block model
  write(IMAIN,*) '    zone 1: depth min/max = ',(1.d0 - NORM_MAX_H_R)*EARTH_R_KM,'/',(1.d0-NORM_MOHO_R)*EARTH_R_KM,'(km)'
  write(IMAIN,*)
  call flush_IMAIN()

  ! temporary array
  allocate(values_at_surface(NP_b,NT_b),stat=ier)
  if (ier /= 0) stop 'Error allocating values_at_surface array'
  values_at_surface(:,:) = 0.0

  ! root directory
  rootdir = 'DATA/gladm15/crust/'

  do ipar = 1,BKMNS_NPAR
    ! reads in parameter in block format
    select case(ipar)
    case (1)
      ! VPV
      filename = trim(rootdir)//'vpv.bin'
      call read_crust_bkmns_model_file(filename,BKMNS_crust_vpv,values_at_surface)
    case (2)
      ! VPH
      filename = trim(rootdir)//'vph.bin'
      call read_crust_bkmns_model_file(filename,BKMNS_crust_vph,values_at_surface)
    case (3)
      ! VSV
      filename = trim(rootdir)//'vsv.bin'
      call read_crust_bkmns_model_file(filename,BKMNS_crust_vsv,values_at_surface)
    case (4)
      ! VSH
      filename = trim(rootdir)//'vsh.bin'
      call read_crust_bkmns_model_file(filename,BKMNS_crust_vsh,values_at_surface)
    case (5)
      ! ETA
      filename = trim(rootdir)//'eta.bin'
      call read_crust_bkmns_model_file(filename,BKMNS_crust_eta,values_at_surface)
    case (6)
      ! RHO
      filename = trim(rootdir)//'rho.bin'
      call read_crust_bkmns_model_file(filename,BKMNS_crust_rho,values_at_surface)
    case default
      call exit_MPI(myrank,'Invalid parameter for crust file in bkmns')
    end select

    ! stores surface values
    BKMNS_crust_surface_value(:,:,ipar) = values_at_surface(:,:)
  enddo

  ! free arrays
  deallocate(values_at_surface)

  end subroutine read_crust_bkmns_model

!
!-----------------------------------------------------------------------------------------
!

  subroutine read_crust_bkmns_model_file(filename,array,values_at_surface)

! reads crustal file (given in binary format, name ending in *.bin)
!
! This expects a file containing 32-bit floats values.
! Byte swapping should be handled automatically if necessary.

  use constants, only: IMAIN,MAX_STRING_LEN,myrank
  use model_bkmns_par

  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: filename
  real, dimension(NP_b,NT_b,NR_b),intent(inout) :: array
  real, dimension(NP_b,NT_b), intent(out) :: values_at_surface

  ! local parameters
  integer :: indx,i,j,k
  real :: val,val_swap
  real,dimension(NR_b) :: val_array
  real,dimension(1) :: val_single
  integer(kind=8) :: filesize
  logical :: file_exists

  ! byte-order endianess
  logical :: byteswap
  character(len=1),parameter :: HEADER_IS_BYTE_SWAPPED = '>'
  character(len=1) :: c   ! 1-byte marker

  real,external :: float_swap

  ! initializes
  values_at_surface(:,:) = 0.0

  ! user output
  write(IMAIN,*) '  crustal file: ',trim(filename)
  call flush_IMAIN()

  ! opens model file
  inquire(file=trim(filename),exist=file_exists)
  if (.not. file_exists) then
    print *,'Error file: ',trim(filename),'does not exist'
    call exit_MPI(myrank,'Error crustal file for bkmns model does not exist')
  endif

  ! reads in binary values from file
  filesize = NP_b * NT_b * NR_b * 4 + 1 ! total filesize = array size (float) + 1 char
  call open_file_abs_r(11, trim(filename), len_trim(filename), filesize)

  ! checks byte ordering
  !
  ! first entry in file is a byte marker (char)
  indx = 1
  call read_abs(11, c, 1, indx)
  byteswap = (c == HEADER_IS_BYTE_SWAPPED)

  write(IMAIN,*) '    byte swap: ',byteswap
  call flush_IMAIN()

  !debug
  !print *,'debug: bkmns crust file ',trim(filename),' - byteswap ',byteswap,'***',c,'***',HEADER_IS_BYTE_SWAPPED,'***'

  ! reads in topography array
  indx = 0
  if (byteswap) then
    ! swapping byte from little- to big-endian or vice verse
    do i = 1,NP_b
      do j = 1,NT_b
        do k = 1,NR_b
          ! reads single real value
          indx = indx + 1
          call read_abs_shifted(11, val_single, 4, indx, 1)  ! 4-byte read into val_single, 1-byte shifted
          val = val_single(1)
          ! stores swapped real value
          val_swap = float_swap(val)
          array(i,j,k) = val_swap
          ! stores last value greater than zero as the value at surface
          ! (for points above topo surface, values will be 0.0)
          if (val_swap > WATER_LEVEL) values_at_surface(i,j) = val_swap
        enddo
      enddo
    enddo
  else
    ! keeps same byte order
    do i = 1,NP_b
      do j = 1,NT_b
        ! slow way: reads single entries
        ! stores one-by-one
        !do k = 1,NR_b
        !  indx = indx + 1
        !  call read_abs_shifted(11, val_single, 4, indx, 1)  ! 4-byte read into val, 1-byte shifted
        !  array(i,j,k) = val
        !  ! stores last value greater than zero as the value at surface
        !  ! (for points above topo surface, values will be 0.0)
        !  if (val > WATER_LEVEL) values_at_surface(i,j) = val
        !enddo

        ! fast way: reads all values in one direction (shifted by 1 byte due to first marker entry)
        ! reads line with length NR_b
        indx = indx + 1
        call read_abs_shifted(11, val_array, 4*NR_b, indx, 1)  ! 4-byte read into val, 1-byte shifted
        ! stores one-by-one
        do k = 1,NR_b
          val = val_array(k)
          array(i,j,k) = val
          ! stores last value greater than zero as the value at surface
          ! (for points above topo surface, values will be 0.0)
          if (val > WATER_LEVEL) values_at_surface(i,j) = val
        enddo
      enddo
    enddo
  endif ! byteswap

  ! closes file
  call close_file_abs(11)

  ! user output
  write(IMAIN,*) '    min/max block values   = ',minval(array),maxval(array)
  write(IMAIN,*) '    min/max surface values = ',minval(values_at_surface),maxval(values_at_surface)
  call flush_IMAIN()

  !debug
  !print *,'debug: bkmns crust array: ',trim(filename),' min/max = ',minval(array),maxval(array)

  end subroutine read_crust_bkmns_model_file

!
!-----------------------------------------------------------------------------------------
!

  real function float_swap(r_in)

! swaps endianness of a real
! based on: http://www.cgd.ucar.edu/cas/software/endian.html

  implicit none

  real,intent(in) :: r_in

  ! local parameters
  integer :: i_element,i_element_br

  ! transfer 32 bits of generic 32 bit integer space
  i_element = transfer(r_in, 0)

  ! reverses order of 4 bytes
  call mvbits( i_element, 24, 8, i_element_br, 0  )
  call mvbits( i_element, 16, 8, i_element_br, 8  )
  call mvbits( i_element,  8, 8, i_element_br, 16 )
  call mvbits( i_element,  0, 8, i_element_br, 24 )

  ! transfers reversed order bytes to 32 bit REAL space
  float_swap = transfer( i_element_br, 0.0 )

  end function float_swap


!
!-----------------------------------------------------------------------------------------
!

  subroutine model_bkmns_mantle_broadcast()

  ! standard routine to setup model

  use constants, only: IMAIN,myrank
  use model_bkmns_par

  implicit none

  integer :: ier

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: bkmns mantle'
    call flush_IMAIN()
  endif

  ! model files only read by main process
  if (myrank == 0) call read_mantle_bkmns_model()

  ! broadcasts array sizes
  call bcast_all_singlei(num_spline_positions_zone2)
  call bcast_all_singlei(num_spline_positions_zone3)
  call bcast_all_singlei(num_spline_positions_zone4)

  call bcast_all_singlei(num_sh_coeffs_zone2)
  call bcast_all_singlei(num_sh_coeffs_zone3)
  call bcast_all_singlei(num_sh_coeffs_zone4)

  call bcast_all_singlei(degree_N_zone2)
  call bcast_all_singlei(degree_N_zone3)
  call bcast_all_singlei(degree_N_zone4)
  call bcast_all_singlei(degree_NMAX)

  call bcast_all_singlei(num_spline_NS_zone2)
  call bcast_all_singlei(num_spline_NS_zone3)
  call bcast_all_singlei(num_spline_NS_zone4)

  ! allocate mantle arrays for all other processes
  if (myrank /= 0) then
    allocate(T_zone2(0:num_spline_positions_zone2-1), &
             T_zone3(0:num_spline_positions_zone3-1), &
             T_zone4(0:num_spline_positions_zone4-1),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating bkmns mantle spline arrays')
    T_zone2(:) = 0.d0; T_zone3(:) = 0.d0; T_zone4(:) = 0.d0

    allocate(A_zone2(num_spline_NS_zone2,num_sh_coeffs_zone2,BKMNS_NPAR), &
             B_zone2(num_spline_NS_zone2,num_sh_coeffs_zone2,BKMNS_NPAR), &
             A_zone3(num_spline_NS_zone3,num_sh_coeffs_zone3,BKMNS_NPAR), &
             B_zone3(num_spline_NS_zone3,num_sh_coeffs_zone3,BKMNS_NPAR), &
             A_zone4(num_spline_NS_zone4,num_sh_coeffs_zone4,BKMNS_NPAR), &
             B_zone4(num_spline_NS_zone4,num_sh_coeffs_zone4,BKMNS_NPAR),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating bkmns mantle sh arrays')
    A_zone2(:,:,:) = 0.d0; B_zone2(:,:,:) = 0.d0
    A_zone3(:,:,:) = 0.d0; B_zone3(:,:,:) = 0.d0
    A_zone4(:,:,:) = 0.d0; B_zone4(:,:,:) = 0.d0

    allocate(Rbasis_zone2(NPT,num_spline_NS_zone2), &
             Rbasis_zone3(NPT,num_spline_NS_zone3), &
             Rbasis_zone4(NPT,num_spline_NS_zone4),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating Rbasis array')
    Rbasis_zone2(:,:) = 0.d0; Rbasis_zone3(:,:) = 0.d0; Rbasis_zone4(:,:) = 0.d0

    !long double *nF = malloc (sizeof (long double[num_sh_coeffs]));
    allocate(nFactors_zone2(0:num_sh_coeffs_zone2-1), &
             nFactors_zone3(0:num_sh_coeffs_zone3-1), &
             nFactors_zone4(0:num_sh_coeffs_zone4-1),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating nFactors array')
    nFactors_zone2(:) = 0.d0; nFactors_zone3(:) = 0.d0; nFactors_zone4(:) = 0.d0

    allocate(Legendre_degree_factor(0:degree_NMAX-1),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating Legendre factor array')
    Legendre_degree_factor(:) = 0.d0
  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(T_zone2,num_spline_positions_zone2)
  call bcast_all_dp(T_zone3,num_spline_positions_zone3)
  call bcast_all_dp(T_zone4,num_spline_positions_zone4)

  call bcast_all_dp(A_zone2,num_spline_NS_zone2 * num_sh_coeffs_zone2 * BKMNS_NPAR)
  call bcast_all_dp(B_zone2,num_spline_NS_zone2 * num_sh_coeffs_zone2 * BKMNS_NPAR)

  call bcast_all_dp(A_zone3,num_spline_NS_zone3 * num_sh_coeffs_zone3 * BKMNS_NPAR)
  call bcast_all_dp(B_zone3,num_spline_NS_zone3 * num_sh_coeffs_zone3 * BKMNS_NPAR)

  call bcast_all_dp(A_zone4,num_spline_NS_zone4 * num_sh_coeffs_zone4 * BKMNS_NPAR)
  call bcast_all_dp(B_zone4,num_spline_NS_zone4 * num_sh_coeffs_zone4 * BKMNS_NPAR)

  call bcast_all_dp(Rbasis_zone2,NPT * num_spline_NS_zone2)
  call bcast_all_dp(Rbasis_zone3,NPT * num_spline_NS_zone3)
  call bcast_all_dp(Rbasis_zone4,NPT * num_spline_NS_zone4)

  call bcast_all_dp(nFactors_zone2,num_sh_coeffs_zone2)
  call bcast_all_dp(nFactors_zone3,num_sh_coeffs_zone3)
  call bcast_all_dp(nFactors_zone4,num_sh_coeffs_zone4)

  call bcast_all_dp(Legendre_degree_factor,degree_NMAX)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine model_bkmns_mantle_broadcast

!
!-----------------------------------------------------------------------------------------
!

  subroutine read_mantle_bkmns_model()

! reads mantle files with spline positions and spherical harmonic expansions

  use constants, only: MAX_STRING_LEN,IIN,myrank
  use model_bkmns_par

  implicit none

  ! local parameters
  integer :: i,i1,i2,idx,ipar,s,n,m,ier,izone
  integer :: nnodes,degree,num_spline_NS,num_spline_positions
  integer :: degree_N,num_sh_coeffs
  double precision :: rmin,rmax
  double precision :: sgn,double_factorial_r

  double precision, dimension(:), pointer :: T_zone
  double precision, dimension(:,:), pointer :: Rbasis
  double precision, dimension(:), pointer :: nFactors

  character(len=MAX_STRING_LEN) :: filename,rootdir
  character(len=MAX_STRING_LEN) :: line,substring
  character(len=1) :: str_zone

  ! root directory
  rootdir = 'DATA/gladm15/mantle/'

  ! zone range within [2,4]
  do izone = 2,4
    ! zone string
    write(str_zone,'(i1)') izone

    ! reads spline file with node positions
    filename = trim(rootdir) // 'nodes_Z' // str_zone // '.dat'

    ! opens model file
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      call exit_MPI(myrank,'Error opening file BKMNS Bspline model file')
    endif
    ! Basis splines description (from top to bottom):
    ! Number of nodes, degree, minimum radius, maximum radius, nodes positions
    read(IIN,*) line ! comment line
    read(IIN,*) nnodes
    read(IIN,*) degree
    read(IIN,*) rmin
    read(IIN,*) rmax

    ! total number of spline positions
    num_spline_positions = nnodes + 2 * degree

    ! allocates spline position array
    select case(izone)
    case (2)
      ! total number of positions
      num_spline_positions_zone2 = num_spline_positions

      allocate(T_zone2(0:num_spline_positions_zone2-1),stat=ier)
      if (ier /= 0) call exit_MPI(0,'Error allocating T_zone2 array')
      T_zone2(:) = 0.0

      ! sets current zone pointer
      !
      ! note: here we assign
      !          T_zone => T_zone2
      !       rather than
      !          T_zone => T_zone2(:)
      !       as we would loose the range from (0:..) bounds in the latter and T_zone(0) would become an error.
      !       another way would be to use
      !          T_zone(0:) => T_zone2(:)
      !       to start the bounds indexing at 0 for the pointer.
      !       this might however not be recognized & implemented by all Fortran compilers.
      !
      T_zone => T_zone2

    case (3)
      ! total number of positions
      num_spline_positions_zone3 = num_spline_positions

      allocate(T_zone3(0:num_spline_positions_zone3-1),stat=ier)
      if (ier /= 0) call exit_MPI(0,'Error allocating T3 array')
      T_zone3(:) = 0.0

      ! sets current zone pointer
      T_zone => T_zone3

    case (4)
      ! total number of positions
      num_spline_positions_zone4 = num_spline_positions

      allocate(T_zone4(0:num_spline_positions_zone4-1),stat=ier)
      if (ier /= 0) call exit_MPI(0,'Error allocating T4 array')
      T_zone4(:) = 0.0

      ! sets current zone pointers
      T_zone => T_zone4

    case default
      call exit_MPI(myrank,'Error invalid zone in read_mantle_bkmns_model() routine')
    end select

    ! set initial positions
    do i = 0,degree
      T_zone(i) = rmin
      T_zone(num_spline_positions - 1 - i) = rmax
    enddo

    ! reads in additional splines node position
    do i = 1,nnodes-2
      idx = degree + i
      read(IIN,*) T_zone(idx)
    enddo

    close(IIN)

    !debug
    !print *,'debug: BKMNS zone ',izone,'spline node positions degree/nnodes = ',degree,nnodes,' T_zone: ',T_zone

    ! coefficient file
    do ipar = 1,BKMNS_NPAR
      ! parameter file name
      select case(ipar)
      case (1)
        ! vpv
        filename = trim(rootdir) // 'mns_Z' // str_zone // '_vpv.dat'
      case (2)
        ! vph
        filename = trim(rootdir) // 'mns_Z' // str_zone // '_vph.dat'
      case (3)
        ! vsv
        filename = trim(rootdir) // 'mns_Z' // str_zone // '_vsv.dat'
      case (4)
        ! vsh
        filename = trim(rootdir) // 'mns_Z' // str_zone // '_vsh.dat'
      case (5)
        ! eta
        filename = trim(rootdir) // 'mns_Z' // str_zone // '_eta.dat'
      case (6)
        ! rho
        filename = trim(rootdir) // 'mns_Z' // str_zone // '_rho.dat'
      case default
        call exit_MPI(myrank,'Invalid parameter for mantle file in bkmns')
      end select

      !debug
      !print *,'debug: BKMNS ipar ',ipar

      ! opens model coefficient file
      open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening file: ',trim(filename)
        call exit_MPI(0,'Error opening file BKMNS coefficient model file')
      endif

      ! file format:
      ! #radial and spherical harmonics coefficients  (A^m_n)_s  (B^m_n)_s  (m <= n <= N)  (s < NS)  N = 100  NS = 6
      ! (A^0_0)_0       =  3.169672E+01
      ! (B^0_0)_0       =  0.000000E+00
      ! (A^0_1)_0       =  4.603058E-02
      ! (B^0_1)_0       =  0.000000E+00
      ! (A^1_1)_0       =  1.566959E-02
      ! (B^1_1)_0       =  2.517482E-02
      ! ..

      read(IIN,'(a160)') line ! comment line
      ! remove leading and trailing whitespace
      line = trim(adjustl(line))

      ! debug
      !print *,'debug: bkmns line: ****',trim(line),'****'

      ! format:
      !#radial and spherical harmonics coefficients  (A^m_n)_s  (B^m_n)_s  (m <= n <= N)  (s < NS)  N = 100  NS = 6
      i1 = index(line,'N =')
      if (i1 == 0) then
        print *,'Error reading file: ',trim(filename), ' - section "N =" index',i1
        print *,'  line: ',line
        call exit_MPI(0,'Error: could not find degree N = info in file')
      endif
      i2 = index(line,'NS =')
      if (i2 == 0) then
        print *,'Error reading file: ',trim(filename), ' - section "NS =" index',i2
        print *,'  line: ',line
        call exit_MPI(0,'Error: could not find degree NS = info in file')
      endif

      ! reads N
      substring = ''
      substring = line(i1+3:i2-1)
      read(substring, *) degree_N                 ! spherical harmonic degree

      !debug
      !print *,'debug: BKMNS N substring ****',trim(substring),'**** degree N = ',degree_N

      ! reads NS
      substring = ''
      substring = line(i2+4:len_trim(line))
      read(substring, *) num_spline_NS            ! NS number of radial splines

      !debug
      !print *,'debug: BKMNS NS substring ****',trim(substring),'**** degree NS = ',num_spline_NS

      num_sh_coeffs = degree_N * (degree_N+1) / 2 + degree_N + 1

      ! allocate spherical harmonics coefficient arrays
      if (ipar == 1) then
        ! we asssume that for a given zone, all the parameter expansions (vpv,vph,..) will have the same degrees
        select case(izone)
        case (2)
          ! maximum number of coefficients (per degree NS)
          num_sh_coeffs_zone2 = num_sh_coeffs
          degree_N_zone2 = degree_N
          num_spline_NS_zone2 = num_spline_NS

          allocate(A_zone2(num_spline_NS_zone2,num_sh_coeffs_zone2,BKMNS_NPAR), &
                   B_zone2(num_spline_NS_zone2,num_sh_coeffs_zone2,BKMNS_NPAR),stat=ier)
          if (ier /= 0) call exit_MPI(0,'Error allocating BKMNS spherical harmonics zone 2 arrays')
          A_zone2(:,:,:) = 0.0
          B_zone2(:,:,:) = 0.0
        case (3)
          ! maximum number of coefficients (per degree NS)
          num_sh_coeffs_zone3 = num_sh_coeffs
          degree_N_zone3 = degree_N
          num_spline_NS_zone3 = num_spline_NS

          allocate(A_zone3(num_spline_NS_zone3,num_sh_coeffs_zone3,BKMNS_NPAR), &
                   B_zone3(num_spline_NS_zone3,num_sh_coeffs_zone3,BKMNS_NPAR),stat=ier)
          if (ier /= 0) call exit_MPI(0,'Error allocating BKMNS spherical harmonics zone 3 arrays')
          A_zone3(:,:,:) = 0.0
          B_zone3(:,:,:) = 0.0
        case (4)
          ! maximum number of coefficients (per degree NS)
          num_sh_coeffs_zone4 = num_sh_coeffs
          degree_N_zone4 = degree_N
          num_spline_NS_zone4 = num_spline_NS

          allocate(A_zone4(num_spline_NS_zone4,num_sh_coeffs_zone4,BKMNS_NPAR), &
                   B_zone4(num_spline_NS_zone4,num_sh_coeffs_zone4,BKMNS_NPAR),stat=ier)
          if (ier /= 0) call exit_MPI(0,'Error allocating BKMNS spherical harmonics zone 4 arrays')
          A_zone4(:,:,:) = 0.0
          B_zone4(:,:,:) = 0.0
        end select
      else
        ! checks degree for ipar /= 1
        select case(izone)
        case(2)
          if (degree_N /= degree_N_zone2 .or. num_spline_NS /= num_spline_NS_zone2) then
            print *,'Error: BKNS mantle model parameter ',ipar,'has different expansion degrees N/NS = ', &
                     degree_N,'/',num_spline_NS
            print *,'  zone 2 should have: ',degree_N_zone2,'/',num_spline_NS_zone2
            call exit_MPI(0,'Error reading BKMNS mantle model parameter zone 2')
          endif
        case(3)
          if (degree_N /= degree_N_zone3 .or. num_spline_NS /= num_spline_NS_zone3) then
            print *,'Error: BKNS mantle model parameter ',ipar,'has different expansion degrees N/NS = ', &
                     degree_N,'/',num_spline_NS
            print *,'  zone 3 should have: ',degree_N_zone3,'/',num_spline_NS_zone3
            call exit_MPI(0,'Error reading BKMNS mantle model parameter zone 3')
          endif
        case(4)
          if (degree_N /= degree_N_zone4 .or. num_spline_NS /= num_spline_NS_zone4) then
            print *,'Error: BKNS mantle model parameter ',ipar,'has different expansion degrees N/NS = ', &
                     degree_N,'/',num_spline_NS
            print *,'  zone 4 should have: ',degree_N_zone4,'/',num_spline_NS_zone4
            call exit_MPI(0,'Error reading BKMNS mantle model parameter zone 4')
          endif
        end select
      endif

      ! reads spherical harmonic coefficients
      do s = 0,num_spline_NS-1   ! number of radial splines
        do n = 0,degree_N         ! spherical harmonic degree n
          do m = 0,n              ! order m
            !mn2Index(m,n)
            idx = n * (n + 1) / 2 + m

            ! format: (A^m_n)_s
            ! (A^5_6)_0       = -6.547580E-03
            ! ..

            ! A coefficient
            read(IIN,'(a80)') line

            ! debug
            !print *,'debug: bkmns read coeff line ****',trim(line),'****'

            i = index(line,'=')
            if (i == 0) then
              print *,'Error reading file: ',trim(filename),' - coefficient A m/n/s = ',m,n,s
              call exit_MPI(myrank,'Error reading A coefficient')
            endif

            ! debug
            !print *,'debug: A coeff ',idx,'line ****',line(i:len_trim(line)),'****'

            select case(izone)
            case (2)
              read(line(i+1:len_trim(line)),*) A_zone2(s+1,idx+1,ipar)  ! radial
            case (3)
              read(line(i+1:len_trim(line)),*) A_zone3(s+1,idx+1,ipar)  ! radial
            case (4)
              read(line(i+1:len_trim(line)),*) A_zone4(s+1,idx+1,ipar)  ! radial
            end select

            ! B coefficient
            read(IIN,'(a80)') line

            i = index(line,'=')
            if (i == 0) then
              print *,'Error reading file: ',trim(filename),' - coefficient B m/n/s = ',m,n,s
              call exit_MPI(myrank,'Error reading B coefficient')
            endif

            select case(izone)
            case (2)
              read(line(i+1:len_trim(line)),*) B_zone2(s+1,idx+1,ipar)  ! spherical
            case (3)
              read(line(i+1:len_trim(line)),*) B_zone3(s+1,idx+1,ipar)  ! spherical
            case (4)
              read(line(i+1:len_trim(line)),*) B_zone4(s+1,idx+1,ipar)  ! spherical
            end select

          enddo
        enddo
      enddo
      close(IIN)
    enddo ! ipar

    ! allocates radial spline values and spherical harmonics factors
    select case(izone)
    case (2)
      !double (*Rb)[ns] = malloc (sizeof (double[NPT][ns]));
      allocate(Rbasis_zone2(NPT,num_spline_NS),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating Rbasis array')
      Rbasis_zone2(:,:) = 0.d0

      !long double *nF = malloc (sizeof (long double[num_sh_coeffs]));
      allocate(nFactors_zone2(0:num_sh_coeffs-1),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating nFactors array')
      nFactors_zone2(:) = 0.d0

      ! sets current zone pointers
      Rbasis => Rbasis_zone2
      nFactors => nFactors_zone2

    case (3)
      !double (*Rb)[ns] = malloc (sizeof (double[NPT][ns]));
      allocate(Rbasis_zone3(NPT,num_spline_NS),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating Rbasis array')
      Rbasis_zone3(:,:) = 0.d0

      !long double *nF = malloc (sizeof (long double[num_sh_coeffs]));
      allocate(nFactors_zone3(0:num_sh_coeffs-1),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating nFactors array')
      nFactors_zone3(:) = 0.d0

      ! sets current zone pointers
      Rbasis => Rbasis_zone3
      nFactors => nFactors_zone3

    case (4)
      !double (*Rb)[ns] = malloc (sizeof (double[NPT][ns]));
      allocate(Rbasis_zone4(NPT,num_spline_NS),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating Rbasis array')
      Rbasis_zone4(:,:) = 0.d0

      !long double *nF = malloc (sizeof (long double[num_sh_coeffs]));
      allocate(nFactors_zone4(0:num_sh_coeffs-1),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating nFactors array')
      nFactors_zone4(:) = 0.d0

      ! sets current zone pointers
      Rbasis => Rbasis_zone4
      nFactors => nFactors_zone4

    end select

    !radialBasis (ns, dg, nnt, r1, r2, T, Rb);
    call get_radial_Basis(num_spline_NS, degree, num_spline_positions, rmin, rmax, T_zone, Rbasis)

    !nmlFactors (N, nlg, nF);
    call get_nml_Factors(degree_N, num_sh_coeffs, nFactors)

  enddo ! izone

  ! pre-computes Legendre polynomial factors
  degree_NMAX = max(degree_N_zone2,degree_N_zone3,degree_N_zone4)

  ! safety check
  if (degree_NMAX > 150) then
    print *,'Error: maximum degree NMAX ',degree_NMAX,' too big, Legendre factors will exceed double precision range'
    call exit_MPI(myrank,'Error maximum degree N exceeds double precision range')
  endif

  ! allocates array (main process only)
  allocate(Legendre_degree_factor(0:degree_NMAX-1),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating Legendre factor array')
  Legendre_degree_factor(:) = 0.d0

  do n = 0,degree_NMAX
    ! from P(i_n + 1,N+2) = sgn (n) * doubleFactorial (2 * n - 1) * pow (sint, n)
    ! pre-computes factors sgn (n) * doubleFactorial (2 * n - 1)
    ! sgn (n)
    if (mod(n,2) == 1) then
      sgn = -1.d0
    else
      sgn = 1.d0
    endif
    ! doubleFactorial(2 * n - 1)
    double_factorial_r = 1.d0
    do i = 2 * n - 1,0,-2
      if (i == 0 .or. i == 1) then
        exit
      else
        double_factorial_r = double_factorial_r * i
      endif
    enddo
    ! stores factor
    Legendre_degree_factor(n) = sgn * double_factorial_r
  enddo
  ! debug
  !print *,'debug: Legendre degree ',degree_NMAX,' factors: ',Legendre_degree_factor(:)

  end subroutine read_mantle_bkmns_model

!
!-----------------------------------------------------------------------------------------
!

  subroutine get_radial_Basis(num_spline_NS, degree, num_spline_positions, rmin, rmax, T_zone, Rbasis)

  use model_bkmns_par, only: NPT

  implicit none

  integer, intent(in) :: num_spline_NS, degree, num_spline_positions
  double precision, intent(in) :: rmin,rmax
  double precision, dimension(0:num_spline_positions-1),intent(in) :: T_zone

  double precision, dimension(NPT,num_spline_NS),intent(inout) :: Rbasis

  ! local parameters
  double precision :: dr,rval
  double precision, external :: get_bSplines

  integer :: i,j

  ! radius increment
  dr = (rmax - rmin) / (NPT - 1)

  ! loops over spline degrees
  do i = 0,num_spline_NS-1
    rval = rmin
    do j = 1,NPT
      ! radial basis using cubic B-splines
      Rbasis(j,i+1) = get_bSplines(rval, num_spline_positions, degree, T_zone, i)
      rval = rval + dr
    enddo
  enddo
  Rbasis(NPT,num_spline_NS) = 1.d0

  !debug
  !do j = 1,NPT
  !  if (Rbasis(j,1) > 0.d0) print *,'debug: radial splines ',Rbasis(j,1),j,num_spline_NS,dr,NPT,rmin,rmax
  !enddo
  !stop 'debug'

  end subroutine get_radial_Basis

!
!-----------------------------------------------------------------------------------------
!

  recursive function get_bSplines(x, num_spline_positions, degree, T_zone, i) result(bSpline)

  implicit none

  double precision, intent(in) :: x
  integer, intent(in) :: num_spline_positions, degree, i
  double precision, dimension(0:num_spline_positions-1), intent(in) :: T_zone

  double precision :: bSpline

  ! local parameters
  double precision :: d1,d2,k1,k2
  double precision :: val,val1,val2

  ! initializes
  bSpline = 0.d0

  ! determines spline value
  if (degree == 0) then
    ! degree-0 spline
    if (T_zone(i) <= x .and. x < T_zone(i+1)) then
      val = 1.d0
    else
      val = 0.d0
    endif
    bSpline = val
  else
    ! spline degrees > 0
    d1 = T_zone(i + degree) - T_zone(i)
    d2 = T_zone(i + degree + 1) - T_zone(i + 1)

    if (d1 > 0.d0) then
      k1 = (x - T_zone(i)) / d1
    else
      k1 = 0.d0
    endif
    if (d2 > 0.d0) then
      k2 = (T_zone(i + degree + 1) - x) / d2
    else
      k2 = 0.d0
    endif
    val1 = get_bSplines(x, num_spline_positions, degree - 1, T_zone, i)
    val2 = get_bSplines(x, num_spline_positions, degree - 1, T_zone, i + 1)

    val = k1 * val1 + k2 * val2

    bSpline = val
  endif

  ! debug
  !if (bSpline > 0.d0) &
  !  print *,'debug: get bsplines ',x,i,degree,num_spline_positions,T_zone(i),T_zone(i+1),' bspline = ',bSpline

  end function get_bSplines

!
!-----------------------------------------------------------------------------------------
!

  subroutine get_nml_Factors(degree_N, num_sh_coeffs, nFactors)

  use constants, only: PI
  use model_bkmns_par, only: LEGENDRE_K1,LEGENDRE_K2

  implicit none

  integer,intent(in) :: degree_N,num_sh_coeffs
  double precision, dimension(0:num_sh_coeffs-1),intent(inout) :: nFactors

  ! local parameters
  integer :: n,m,idx,i,i0
  double precision :: val,sgn,factor
  real(kind=16) :: factorial_r

  do n = 0,degree_N
    val = sqrt ((2.d0 * n + 1) / (4.d0 * PI))
    do m = 0,n
      !mn2Index(m,n)
      idx = n * (n + 1) / 2 + m

      ! sign
      !sgn(m)
      if (mod(m,2) == 1) then
        sgn = -1.d0
      else
        sgn = 1.d0
      endif

      ! factor (m , n)
      if (m == 0) then
        factor = LEGENDRE_K1
      else
        ! factorial (n+m)!/(n-m)!
        ! note: the factorial number can become very large for n > 30 and overflow double precision;
        !       we use here real(kind=16), the maximum size for a floating point number by Fortran standards
        i0 = n - m + 1
        factorial_r = i0
        do i = i0 + 1, n + m
          factorial_r = factorial_r * i
        enddo
        factor = 2.d0 * real(LEGENDRE_K2 / factorial_r,kind=16)
        factor = sqrt(factor)
      endif

      ! stores factor
      !nF(i) = sgn(m) * val * factor (m , n)
      nFactors(idx) = sgn * val * factor
    enddo
  enddo

  end subroutine get_nml_Factors

!
!-----------------------------------------------------------------------------------------
!
! model value evaluations
!
!-----------------------------------------------------------------------------------------
!

  subroutine model_bkmns_crust(lat,lon,radius,vpvc,vphc,vsvc,vshc,etac,rhoc)

! gets crustal value for location lat/lon/r

  use constants, only: PI
  use model_bkmns_par

  implicit none

  ! lat/lon  - in degrees (range lat/lon = [-90,90] / [-180,180]
  ! radius   - normalized by globe radius [0,1.x]
  double precision,intent(in) :: lat,lon,radius

  double precision,intent(inout) :: vpvc,vphc,vsvc,vshc,etac,rhoc

  ! local parameters
  double precision :: phi,theta
  double precision :: Bk_vpv,Bk_vph,Bk_vsv,Bk_vsh,Bk_eta,Bk_rho

  ! converts to colatitude theta/phi in radians
  theta = (90.d0 - lat) * PI/180.d0   ! colatitude between [0,pi]
  phi = lon * PI/180.d0               ! longitude between [-pi,pi]

  ! gets crustal values
  if (radius >= NORM_MOHO_R .and. radius <= NORM_MAX_H_R) then
    ! position within block model range [-80km,topo]
    call bkmns_block2crust(phi,theta,radius,Bk_vpv,Bk_vph,Bk_vsv,Bk_vsh,Bk_eta,Bk_rho)

    ! sets crustal value from block values if non-zero
    if (Bk_vpv > WATER_LEVEL) then
      vpvc = Bk_vpv
      vphc = Bk_vph
      vsvc = Bk_vsv
      vshc = Bk_vsh
      etac = Bk_eta
      rhoc = Bk_rho
    endif
  endif

  end subroutine model_bkmns_crust

!
!-----------------------------------------------------------------------------------------
!

  subroutine bkmns_block2crust(phi,theta,radius,vpvc,vphc,vsvc,vshc,etac,rhoc)

  use constants, only: PI,myrank,EARTH_R_KM,GRAV
  use shared_parameters, only: R_PLANET,RHOAV

  use model_bkmns_par

  implicit none

  ! longitude between [-pi,pi], colatitude between [0,pi], r normalized [0,1.x]
  double precision, intent(in) :: phi,theta,radius
  double precision, intent(inout) :: vpvc,vphc,vsvc,vshc,etac,rhoc

  ! local parameters
  integer :: i0,j0,k0,i1,j1,k1
  integer :: ipar
  real,dimension(:,:,:), pointer :: BKMNS_crust_par

  double precision :: p,t,r
  double precision :: dp,dt,dr
  double precision :: p0,t0,r0,p1,t1,r1
  double precision :: scaleval_rho,scaleval_vel

  ! interpolation
  integer :: n
  double precision :: f000,f100,f010,f110,f001,f101,f011,f111
  double precision :: k000,k100,k010,k110,k001,k101,k011,k111
  double precision :: fac000,fac100,fac010,fac110,fac001,fac101,fac011,fac111
  double precision :: sum,mean,val

  ! initializes
  vpvc = 0.d0 !km/s
  vphc = 0.d0 !km/s
  vsvc = 0.d0 !km/s
  vshc = 0.d0 !km/s
  etac = 0.d0
  rhoc = 0.d0 !g/cm^3

  ! input position
  p = phi
  t = theta
  r = radius

  ! limits input position
  ! phi in [-PI,PI]
  if (p < -PI) p = -PI
  if (p > PI) p = PI
  ! theta in [0,PI]
  if (t < 0.d0) t = 0.d0
  if (t > PI) t = PI
  ! radius
  if (r < NORM_MOHO_R)  r = NORM_MOHO_R
  if (r > NORM_MAX_H_R) r = NORM_MAX_H_R

  ! re-positions ranges
  p = p + PI                ! range [0,2PI]
  t = t                     ! range [0,PI]
  r = r - NORM_MOHO_R       ! range [0,1.x] with 0 being at the base depth of -80km

  ! increments
  dp = 2.d0 * PI / (NP_b - 1)
  dt = PI / (NT_b - 1)
  dr = (NORM_MAX_H_R - NORM_MOHO_R) / (NR_b - 1)

  ! block index
  ! note: index will start at 0
  i0 = int(p / dp)
  j0 = int(t / dt)
  k0 = int(r / dr)

  ! block position
  p0 = i0 * dp
  t0 = j0 * dt
  r0 = k0 * dr

  ! next block position
  p1 = p0 + dp
  t1 = t0 + dt
  r1 = r0 + dr

  ! puts index to range [1,N] for Fortran array
  i0 = i0 + 1
  j0 = j0 + 1
  k0 = k0 + 1

  ! checks
  if (i0 == 0 .or. j0 == 0 .or. k0 == 0) then
    print *,'Error: invalid block index ',i0,j0,k0
    call exit_MPI(myrank,'Error block2crust invalid block position')
  endif
  if (i0 > NP_b .or. j0 > NT_b .or. k0 > NR_b) then
    print *,'Error: invalid block index ',i0,j0,k0
    call exit_MPI(myrank,'Error block2crust invalid block position')
  endif

  ! next block index (for interpolation)
  ! note: index must be in range [1,N]
  if (i0 < NP_b) then
    i1 = i0 + 1
  else
    i1 = i0
  endif
  if (j0 < NT_b) then
    j1 = j0 + 1
  else
    j1 = j0
  endif
  if (k0 < NR_b) then
    k1 = k0 + 1
  else
    k1 = k0
  endif

  ! linear interpolation factors
  fac000 = (p1 - p) * (t1 - t) * (r1 - r)
  fac100 = (p - p0) * (t1 - t) * (r1 - r)
  fac010 = (p1 - p) * (t - t0) * (r1 - r)
  fac110 = (p - p0) * (t - t0) * (r1 - r)
  fac001 = (p1 - p) * (t1 - t) * (r - r0)
  fac101 = (p - p0) * (t1 - t) * (r - r0)
  fac011 = (p1 - p) * (t - t0) * (r - r0)
  fac111 = (p - p0) * (t - t0) * (r - r0)

  ! non-dimensionalize (from km/s)
  scaleval_rho = 1000.0d0 / RHOAV
  scaleval_vel = 1000.0d0 / (R_PLANET * sqrt(PI*GRAV*RHOAV))

  ! loops over parameters
  do ipar = 1,6
    ! gets block values (neighboring blocks used for interpolation)
    select case(ipar)
    case(1)
      ! vpv
      BKMNS_crust_par => BKMNS_crust_vpv
    case(2)
      ! vph
      BKMNS_crust_par => BKMNS_crust_vph
    case(3)
      ! vsv
      BKMNS_crust_par => BKMNS_crust_vsv
    case(4)
      ! vsh
      BKMNS_crust_par => BKMNS_crust_vsh
    case(5)
      ! eta
      BKMNS_crust_par => BKMNS_crust_eta
    case(6)
      ! rho
      BKMNS_crust_par => BKMNS_crust_rho
    end select

    f000 = BKMNS_crust_par(i0,j0,k0); f100 = BKMNS_crust_par(i1,j0,k0)
    f010 = BKMNS_crust_par(i0,j1,k0); f110 = BKMNS_crust_par(i1,j1,k0)
    f001 = BKMNS_crust_par(i0,j0,k1); f101 = BKMNS_crust_par(i1,j0,k1)
    f011 = BKMNS_crust_par(i0,j1,k1); f111 = BKMNS_crust_par(i1,j1,k1)

    ! interpolation
    sum = 0.d0
    n = 0
    if (f000 >= WATER_LEVEL) then
      sum = sum + f000; n = n + 1
    endif
    if (f100 >= WATER_LEVEL) then
      sum = sum + f100; n = n + 1
    endif
    if (f010 >= WATER_LEVEL) then
      sum = sum + f010; n = n + 1
    endif
    if (f110 >= WATER_LEVEL) then
      sum = sum + f110; n = n + 1
    endif
    if (f001 >= WATER_LEVEL) then
      sum = sum + f001; n = n + 1
    endif
    if (f101 >= WATER_LEVEL) then
      sum = sum + f101; n = n + 1
    endif
    if (f011 >= WATER_LEVEL) then
      sum = sum + f011; n = n + 1
    endif
    if (f111 >= WATER_LEVEL) then
      sum = sum + f111; n = n + 1
    endif

    ! checks if values found
    if (n == 0) then
      ! all block values are zero - position might be above surface

      !debug
      !print *,'Error: rank ',myrank,' block2crust position has no valid values'
      !print *,'  position lon/lat/r = ',phi*180.0/PI,90.0-theta*180.0/PI,radius*EARTH_R_KM,' (km)',' phi/theta/radius = ',p,t,r
      !print *,'  parameter ',ipar,' f values: ',f000,f100,f010,f110,f001,f101,f011,f111
      !print *,'  block position : lon/lat,r = ',(p0-PI)*180.0/PI,90.0-t0*180.0/PI,r0*EARTH_R_KM,' index i/j/k',i0,j0,k0
      !print *,'  taking surface value: ',BKMNS_crust_surface_value(i0,j0,ipar)
      !do k1 = k0,1,-1
      !  f000 = BKMNS_crust_par(i0,j0,k1)
      !  print *,'    depth profile index ',k1,' radius ',k1 * dr * EARTH_R_KM,' par = ',f000
      !enddo
      !call exit_MPI(myrank,'Error block2crust invalid position value')

      ! sets single block value to surface value
      f000 = BKMNS_crust_surface_value(i0,j0,ipar)
      sum = sum + f000; n = n + 1
    endif

    ! takes average for missing blocks
    if (n < 8) then
      mean = sum / n
      if (f000 < WATER_LEVEL) f000 = mean
      if (f100 < WATER_LEVEL) f100 = mean
      if (f010 < WATER_LEVEL) f010 = mean
      if (f110 < WATER_LEVEL) f110 = mean
      if (f001 < WATER_LEVEL) f001 = mean
      if (f101 < WATER_LEVEL) f101 = mean
      if (f011 < WATER_LEVEL) f011 = mean
      if (f111 < WATER_LEVEL) f111 = mean
    endif

    ! linear interpolation with position
    k000 = f000 * fac000
    k100 = f100 * fac100
    k010 = f010 * fac010
    k110 = f110 * fac110
    k001 = f001 * fac001
    k101 = f101 * fac101
    k011 = f011 * fac011
    k111 = f111 * fac111
    val = (k000 + k100 + k010 + k110 + k001 + k101 + k011 + k111) / (dp * dt * dr)

    ! sets corresponding parameter value
    select case(ipar)
    case(1)
      vpvc = val * scaleval_vel
    case(2)
      vphc = val * scaleval_vel
    case(3)
      vsvc = val * scaleval_vel
    case(4)
      vshc = val * scaleval_vel
    case(5)
      etac = val  ! eta is non-dimensional already
    case(6)
      rhoc = val * scaleval_rho
    end select
  enddo

  end subroutine bkmns_block2crust

!
!-----------------------------------------------------------------------------------------
!

  subroutine model_bkmns_mantle(lat,lon,radius,vpv,vph,vsv,vsh,eta,rho)

  use constants, only: PI
  use model_bkmns_par

  implicit none

  ! lat/lon  - in degrees (range lat/lon = [-90,90] / [-180,180]
  ! radius   - normalized by globe radius [0,1.x]
  double precision,intent(in) :: lat,lon,radius

  ! absolute values, not perturbations
  double precision,intent(inout) :: vpv,vph,vsv,vsh,eta,rho

  ! local parameters
  double precision :: r,theta,phi
  double precision :: M_vpv,M_vph,M_vsv,M_vsh,M_eta,M_rho

  ! initializes
  M_vpv = 0.d0
  M_vph = 0.d0
  M_vsv = 0.d0
  M_vsh = 0.d0
  M_eta = 0.d0
  M_rho = 0.d0

  ! point location
  ! converts to colatitude theta/phi in radians
  theta = (90.d0 - lat) * PI/180.d0   ! colatitude between [0,pi]
  phi = lon * PI/180.d0               ! longitude between [-pi,pi]
  r = radius                          ! radius (normalized)

  ! gets crustal values
  if (r >= NORM_MOHO_R .and. r <= NORM_MAX_H_R) then
    ! position within block model range [-80km,topo]
    call bkmns_block2crust(phi,theta,r,M_vpv,M_vph,M_vsv,M_vsh,M_eta,M_rho)
  else
    ! position in mantle zones
    call bkmns_get_mantle_value(phi,theta,r,M_vpv,M_vph,M_vsv,M_vsh,M_eta,M_rho)
  endif

  ! returns mantle values if non-zero
  if (M_vpv > WATER_LEVEL) then
    vpv = M_vpv
    vph = M_vph
    vsv = M_vsv
    vsh = M_vsh
    eta = M_eta
    rho = M_rho
  endif

  end subroutine model_bkmns_mantle

!
!-----------------------------------------------------------------------------------------
!

  subroutine bkmns_get_mantle_value(phi,theta,radius,vpv,vph,vsv,vsh,eta,rho)

  use constants, only: PI,EARTH_R_KM,GRAV
  use shared_parameters, only: R_PLANET,RHOAV
  use model_bkmns_par

  implicit none

  ! longitude between [-pi,pi], colatitude between [0,pi], r normalized [0,1.x]
  double precision,intent(in) :: phi,theta,radius

  ! absolute values, not perturbations
  double precision,intent(inout) :: vpv,vph,vsv,vsh,eta,rho

  ! local parameters
  integer :: ri,n,m,s,idx
  integer :: izone,ipar
  integer :: num_sh_coeffs,degree_N,num_spline_NS

  double precision :: rmin,rmax
  double precision :: Cmn,Smn
  double precision, dimension(0:degree_NMAX) :: C_coeffs,S_coeffs
  double precision, dimension(0:degree_NMAX+1) :: P,Pnormalized

  double precision, dimension(:,:), pointer :: Rbasis
  double precision, dimension(:), pointer :: nFactors
  double precision, dimension(:,:,:), pointer  :: A_zone,B_zone

  ! mantle parameters (1==vpv, 2==vph, ..)
  double precision, dimension(BKMNS_NPAR) :: M_par
  double precision :: scaleval_vel,scaleval_rho

  ! initializes mantle parameter values
  M_par(:) = 0.d0

  ! loops over mantle zones
  do izone = 2,4
    ! gets min/max for shell radius
    select case(izone)
    case (2)
      rmin = T_zone2(0)
      rmax = T_zone2(num_spline_positions_zone2-1)
    case (3)
      rmin = T_zone3(0)
      rmax = T_zone3(num_spline_positions_zone3-1)
    case (4)
      rmin = T_zone4(0)
      rmax = T_zone4(num_spline_positions_zone4-1)
    end select

    ! checks if in zone shell min/max radius
    if (radius >= rmin .and. radius <= rmax) then
      ! selects corresponding zone arrays
      select case(izone)
      case (2)
        A_zone => A_zone2
        B_zone => B_zone2
        Rbasis => Rbasis_zone2
        nFactors => nFactors_zone2
        degree_N = degree_N_zone2
        num_spline_NS = num_spline_NS_zone2
        num_sh_coeffs = num_sh_coeffs_zone2
      case (3)
        A_zone => A_zone3
        B_zone => B_zone3
        Rbasis => Rbasis_zone3
        nFactors => nFactors_zone3
        degree_N = degree_N_zone3
        num_spline_NS = num_spline_NS_zone3
        num_sh_coeffs = num_sh_coeffs_zone3
      case (4)
        A_zone => A_zone4
        B_zone => B_zone4
        Rbasis => Rbasis_zone4
        nFactors => nFactors_zone4
        degree_N = degree_N_zone4
        num_spline_NS = num_spline_NS_zone4
        num_sh_coeffs = num_sh_coeffs_zone4
      end select

      ! radial index
      ri = int( (radius - rmin)/(rmax - rmin) * (NPT - 1) + 0.5 ) + 1   ! + 1 due to array starting at 1

      ! initialize Legendre polynomials and coefficients
      P(:) = 0.d0
      Pnormalized(:) = 0.d0
      C_coeffs(:) = 0.d0
      S_coeffs(:) = 0.d0

      ! azimutal Basis
      call get_azimuthal_Basis(phi,degree_N,C_coeffs,S_coeffs)

      do n = 0,degree_N
        ! polar Basis
        call get_polar_Basis(theta,n,num_sh_coeffs,nFactors,P,Pnormalized)

        do m = 0,n
          Cmn = Pnormalized(m + 1) * C_coeffs(m)
          Smn = Pnormalized(m + 1) * S_coeffs(m)

          !mn2Index(m,n)
          idx = n * (n + 1) / 2 + m + 1   ! + 1 due to arrays starting at 1

          !debug
          !print *,'debug: position lon/lat/r ',phi*180.0/PI,90.0-theta*180.0/PI,radius*EARTH_R_KM, &
          !        ' index ri ',ri,' Rbasis ',Rbasis(ri,1),' Cmn/Smn ',Cmn,Smn, &
          !        ' idx ',idx,' A/B zone ',A_zone(1,idx,1),B_zone(1,idx,1)

          ! spline contributions
          do s = 1,num_spline_NS              ! starting at [1,NS] due to arrays starting at 1
            if (Rbasis(ri,s) > 0.d0) then
              do ipar = 1,BKMNS_NPAR
                M_par(ipar) = M_par(ipar) + Rbasis(ri,s) * (Cmn * A_zone(s,idx,ipar) + Smn * B_zone(s,idx,ipar))
              enddo
            endif
          enddo
        enddo
      enddo

      ! not used: calculates perturbations
      !if (dvv) M[i] = log (M[i] / meanModel (r, nl, Mm));

      !debug
      !print *,'debug: position lon/lat/r ',phi*180.0/PI,90.0-theta*180.0/PI,radius*EARTH_R_KM, &
      !        ' vpv ',M_par(1),' vph ',M_par(2),' vsv ',M_par(3),' vsh ',M_par(4),' eta ',M_par(5),' rho ',M_par(6)

      ! non-dimensionalize (from km/s)
      scaleval_rho = 1000.0d0 / RHOAV
      scaleval_vel = 1000.0d0 / (R_PLANET * sqrt(PI*GRAV*RHOAV))

      ! sets corresponding parameter value
      vpv = M_par(1) * scaleval_vel
      vph = M_par(2) * scaleval_vel
      vsv = M_par(3) * scaleval_vel
      vsh = M_par(4) * scaleval_vel
      eta = M_par(5)                 ! eta is non-dimensional already
      rho = M_par(6) * scaleval_rho

      ! found mantle values, all done
      return

    endif ! r

  enddo ! izone

  end subroutine bkmns_get_mantle_value


!
!-----------------------------------------------------------------------------------------
!

  subroutine get_azimuthal_Basis(phi,degree_N,C_coeffs,S_coeffs)

  use model_bkmns_par, only: degree_NMAX

  implicit none

  double precision, intent(in) :: phi
  integer, intent(in) :: degree_N
  double precision,dimension(0:degree_NMAX),intent(inout) :: C_coeffs,S_coeffs

  ! local parameters
  integer :: m

  do m = 0,degree_N
    C_coeffs(m) = cos(m * phi)
    S_coeffs(m) = sin(m * phi)
  enddo

  end subroutine get_azimuthal_Basis

!
!-----------------------------------------------------------------------------------------
!

  subroutine get_polar_Basis(theta,n,num_sh_coeffs,nFactors,P,Pnormalized)

  use model_bkmns_par, only: LEGENDRE_K3,Legendre_degree_factor,degree_NMAX

  implicit none

  double precision, intent(in) :: theta
  integer, intent(in) :: n,num_sh_coeffs
  double precision,dimension(0:num_sh_coeffs-1),intent(in) :: nFactors

  double precision,dimension(0:degree_NMAX+1),intent(inout) :: P,Pnormalized

  ! local parameters
  integer :: m,idx
  double precision :: z,sint,cott
  double precision :: P0,factor

  ! calculates normalized spherical harmonics based on associated Legendre Polynomials
  ! Dahlen & Tromp, 1998

  ! spherical harmonic P
  ! Pmn (cos (t), n, N, P);
  z = cos(theta)

  if (abs(z) - 1.d0 < 1.d-9) then
    sint = 0.d0
    cott = 0.d0
  else
    sint = sqrt (1.d0 - z*z)
    cott = z / sint
  endif

  ! associated Legrendre polynomials Pmn
  if (n == 0) then
    P(1) = 1.d0
  else if (n == 1) then
    P(0) = P(1)
    P(1) = z
    P(2) = -sint
  else
    P0 = P(0)
    P(0) = P(1)
    P(1) = ((2 * n - 1) * z * P(1) - (n - 1) * P0) / n
    !P(i_n + 1,N+2) = sgn (n) * doubleFactorial (2 * n - 1) * pow (sint, n)
    factor = Legendre_degree_factor(n)  ! sgn (n) * doubleFactorial (2 * n - 1)
    P(n + 1) = factor * sint**n

    P(n) = P(n) * z * (2 * (n - 1) + 1)

    do m = n - 1, 1+1, -1
      P(m) = (-2.d0 * m * cott * P(m + 1) - P(m + 2)) / ((n + m) * (n - m + 1))
    enddo
  endif

  ! normalized spherical harmonic Pnormalized
  ! normalize (n, nlg, N, nF, P, nP);
  do m = 0,n
    !mn2Index(m,n)
    idx = n * (n + 1) / 2 + m
    Pnormalized(m + 1) = LEGENDRE_K3 * nFactors(idx) * P(m + 1)
  enddo

  end subroutine get_polar_Basis

