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
!
! Spherical harmonic model for Mars
!
! spherical harmonics format given by Ina's evolution model output
! uses the case65TAY 1D reference model defined in model_case65TAY.f90
!
! based on:
!   A.-C. Plesa, E. Bozdag, A. Rivoldini, M. Knapmeyer, S. M. McLennan, S. Padovan, N. Tosi, D. Breuer,
!   D. Peter, S. Staehler, M. A. Wieczorek, M. van Driel, A. Khan, T. Spohn, 2021.
!   Seismic Velocity Variations in a 3D Martian Mantle: Implications for the InSight Measurements,
!   JGR Planets, 126, 6, e2020JE006755.
!   https://doi.org/10.1029/2020JE006755
!
!
!--------------------------------------------------------------------------------------------------


  module model_sh_mars_par

  ! spherical harmonics mars model for 3D crust/mantle velocity variations
  ! model extend from the CMB at 1,550 km depth to surface, where a radius of 3,389.5 km is originally used for Mars

  ! Mars radius (in km)
  double precision, parameter :: MARS_R_KM               = 3389.5d0
  ! Mars CMB radius (in km)
  double precision, parameter :: MARS_CMB_R_KM           = 1839.5d0   ! CMB at 1,550 km depth

  ! model root directory
  character(len=160),parameter :: SH_model_rootdir = 'DATA/mars/SH_model/'

  ! model file specification (specifies where to find density, vp and vs spherical-harmonic files)
  character(len=160),parameter :: model_specs = 'SH_model_files.dat'

  ! number of model variables: isotropic rho,vp,vs
  integer,parameter :: NPAR = 3

  ! number of shell layers
  integer :: N_shells

  ! spherical harmonic degree
  integer :: degree_N

  ! number of spherical harmonic degree coefficients
  integer :: num_sh_coeffs

  ! Water level to prevent divisions by zero
  double precision,parameter :: WATER_LEVEL     = 1d-15

  ! spherical harmonics coefficients (for each shell)
  double precision,dimension(:,:,:),allocatable, target :: A_shell, B_shell

  ! shell radii
  double precision, dimension(:), allocatable, target :: R_shell

  ! Legendre polynomials normalization factors
  double precision, dimension(:), allocatable, target :: nFactors

  ! factors for Legendre polynomials
  double precision, dimension(:), allocatable :: Legendre_degree_factor

  end module model_sh_mars_par


!
!-----------------------------------------------------------------------------------------
!

  subroutine model_SH_mars_broadcast()

  ! standard routine to setup model

  use constants, only: IMAIN,myrank
  use model_sh_mars_par

  implicit none

  integer :: ier

  ! checks if already done
  if (allocated(A_shell)) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: SH mars model (crust/mantle)'
    call flush_IMAIN()
  endif

  ! model files only read by main process
  if (myrank == 0) call read_SH_mars_model()

  ! broadcasts array sizes
  call bcast_all_singlei(N_shells)
  call bcast_all_singlei(degree_N)
  call bcast_all_singlei(num_sh_coeffs)

  ! allocate SH mantle arrays for all other processes
  if (myrank /= 0) then
    allocate(A_shell(NPAR,num_sh_coeffs,N_shells), &
             B_shell(NPAR,num_sh_coeffs,N_shells),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating SH mars model arrays')
    A_shell(:,:,:) = 0.d0; B_shell(:,:,:) = 0.d0

    allocate(R_shell(N_shells),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating shell array')
    R_shell(:) = 0.d0

    !long double *nF = malloc (sizeof (long double[num_sh_coeffs]));
    allocate(nFactors(0:num_sh_coeffs-1),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating nFactors array')
    nFactors(:) = 0.d0

    allocate(Legendre_degree_factor(0:degree_N),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating Legendre factor array')
    Legendre_degree_factor(:) = 0.d0
  endif

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(A_shell,N_shells * num_sh_coeffs * NPAR)
  call bcast_all_dp(B_shell,N_shells * num_sh_coeffs * NPAR)
  call bcast_all_dp(R_shell,N_shells)
  call bcast_all_dp(nFactors,num_sh_coeffs)
  call bcast_all_dp(Legendre_degree_factor,degree_N+1)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine model_SH_mars_broadcast

!
!-----------------------------------------------------------------------------------------
!

  subroutine read_SH_mars_model()

! reads SH files with spherical harmonic expansions

  use constants, only: MAX_STRING_LEN,TINYVAL,IIN,IMAIN,myrank
  use model_sh_mars_par

  implicit none

  ! local parameters
  integer :: i,ipar,ier,ilayer,n,m,idx
  integer :: degree_N_read,N_shells_read,n_read,m_read
  double precision :: A_coeff,B_coeff
  double precision :: shell_radius,r_norm

  character(len=MAX_STRING_LEN) :: filename,filename_rho,filename_vp,filename_vs
  character(len=MAX_STRING_LEN) :: line,substring

  ! reads in model file specifications
  filename = trim(SH_model_rootdir) // trim(model_specs)
  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call exit_MPI(0,'Error opening file SH mars model file specifications')
  endif

  ! reads filenames
  filename_rho = ''
  filename_vp = ''
  filename_vs = ''

  do while (ier == 0)
    read(IIN,'(a512)',iostat=ier) line ! comment line
    if (ier == 0) then
      ! remove leading and trailing whitespace
      line = trim(adjustl(line))
      if (len_trim(line) == 0) cycle
      if (line(1:1) == '#' .or. line(1:1) == '!') cycle

      ! no comment, then exit loop
      ! format: #parameter : #filename
      ! for example:
      !density: case65/SH_Case65_TAY13_Density_2102
      !vp     : case65/SH_Case65_TAY13_Vp_2102
      !vs     : case65/SH_Case65_TAY13_Vs_2102
      i = index(line,":")
      if (i > 1) then
        substring = line(1:i-1)
        filename = line(i+1:len_trim(line))
      else
        print *,'Error reading SH model specification: ',trim(model_specs)
        print *,'  line: ',trim(line)
        call exit_MPI(0,'Error reading SH model specification file')
      endif

      ! remove leading and trailing whitespace from filename
      filename = trim(adjustl(filename))

      ! parameter type
      ! converts all string characters to lowercase (to make user input case-insensitive)
      call convert_to_lowercase(substring,substring)
      !debug
      !print *,'debug: SH mars model: parameter = ',trim(substring),' file = ',trim(filename)

      ! sets corresponding filename
      select case(trim(substring))
      case('density')
        filename_rho = trim(SH_model_rootdir) // trim(filename)
      case('vp')
        filename_vp = trim(SH_model_rootdir) // trim(filename)
      case('vs')
        filename_vs = trim(SH_model_rootdir) // trim(filename)
      case default
        ! not recognized
        print *,'Error reading SH model specification: ',trim(model_specs)
        print *,'  line: ',trim(line)
        print *,'  parameter: ',trim(substring)
        print *,'  filename : ',trim(filename)
        call exit_MPI(0,'Error reading SH model specification parameter file')
      end select
    endif
  enddo

  ! checks if all parameter files read in
  if (len_trim(filename_rho) == 0 .or. len_trim(filename_vp) == 0 .or. len_trim(filename_vs) == 0) then
    print *,'Error reading SH model files from ',trim(model_specs)
    call exit_MPI(0,'Error reading SH mars model parameter file names')
  endif

  ! user output
  write(IMAIN,*) '  reading spherical coefficients files: '
  call flush_IMAIN()

  ! reads in coefficient file
  do ipar = 1,NPAR
    ! parameter file name
    select case(ipar)
    case (1)
      ! Density
      filename = trim(filename_rho)
    case (2)
      ! Vp
      filename = trim(filename_vp)
    case (3)
      ! Vs
      filename = trim(filename_vs)
    case default
      call exit_MPI(myrank,'Invalid parameter for SH mars file in read_SH_mars_model() routine')
    end select

    ! user output
    write(IMAIN,*) '    ',trim(filename)
    call flush_IMAIN()

    ! opens model coefficient file
    open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(filename)
      call exit_MPI(0,'Error opening file SH mars coefficient file')
    endif

    ! file format:
    !# Simulation ID: Case65_TAY13
    !# Crustal thickness model: 3200_1_DWTh2Ref1
    !# Grid file: Mars_ico66_res15_24p6_d1550.grid
    !# Simulation iteration: 2102
    !# Simulation time: 4.503460760464611 Gyr
    !# Number of shells of SH expansion: 68
    !# Maximum degree of SH expansion: 20
    !68 20
    !1.839500e+06
    !0 0 4.041686e+03 0.000000e+00
    !1 0 -1.450610e+00 0.000000e+00
    !1 1 -6.312329e-01 9.913989e-02
    !..

    ! skips header lines
    do while (ier == 0)
      read(IIN,'(a160)',iostat=ier) line
      if (ier == 0) then
          ! remove leading and trailing whitespace
          line = trim(adjustl(line))
          if (len_trim(line) == 0) cycle
          if (line(1:1) == '#' .or. line(1:1) == '!') cycle
          ! no comment, then exit loop
          exit
      endif
    enddo

    ! reads expansion numbers
    ! format: #number_of_spherical_shells  #SH_degree
    read(line,*) N_shells_read,degree_N_read

    ! checks
    if (N_shells_read <= 1) then
      print *,'Error: number of shells ',N_shells_read,' must be at least 2.'
      call exit_MPI(0,'Error invalid number of shells')
    endif
    if (degree_N_read < 0) then
      print *,'Error: invalid spherical harmonic degree ',degree_N_read
      call exit_MPI(0,'Error invalid spherical harmonic degree')
    endif

    ! allocates arrays
    ! note: we assume that all parameter files (rho,vp,vs) will have the same spherical harmonic degree
    !       and number of shells
    if (ipar == 1) then
      ! sets model parameters
      degree_N = degree_N_read
      N_shells = N_shells_read

      ! number of SH coefficients
      num_sh_coeffs = degree_N * (degree_N+1) / 2 + degree_N + 1

      ! allocates spherical harmonics coefficient arrays
      allocate(A_shell(NPAR,num_sh_coeffs,N_shells), &
               B_shell(NPAR,num_sh_coeffs,N_shells),stat=ier)
      if (ier /= 0) call exit_MPI(0,'Error allocating spherical harmonics arrays')
      A_shell(:,:,:) = 0.0
      B_shell(:,:,:) = 0.0

      ! non-dimensionalizes shell radii
      allocate(R_shell(N_shells),stat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error allocating shell array')
      R_shell(:) = 0.d0
    else
      ! checks if degree match between different parameter files
      if (degree_N_read /= degree_N) then
        print *,'Error: invalid SH degree N ',degree_N_read,', should match ',degree_N
        print *,'  parameter file ',trim(filename)
        print *,'  All SH mars model parameter files must have same SH degree, please check with density file.'
        call exit_MPI(0,'Invalid SH degree in parameter file')
      endif
      if (N_shells_read /= N_shells) then
        print *,'Error: invalid number of shells ',N_shells_read,', should match ',N_shells
        print *,'  parameter file ',trim(filename)
        print *,'  All SH mars model parameter files must have same number of shells, please check with density file.'
        call exit_MPI(0,'Invalid number of shells in parameter file')
      endif
    endif

    ! loops of spherical shells
    do ilayer = 1,N_shells
      ! reads radius
      read(IIN,*) shell_radius ! shell radius (in m), e.g. 1.839500e+06

      ! stores non-dimensionalized radius
      ! note: this will scale the radii from the normalizes CMB radius to 0 [CMB_R/MARS_R, 1]
      !       therefore, even if the GLL mesh uses a different mars radius, like 3390km, the normalized
      !       GLL position will go up to 1, as if the shells here are stretched between the CMB and the surface radius.
      r_norm = shell_radius / (MARS_R_KM * 1000.d0)

      ! stores radii
      if (ipar == 1) then
        R_shell(ilayer) = r_norm
      else
        ! checks if all parameter files have the same shell layering, otherwise something must be off
        if (abs(shell_radius - R_shell(ilayer) * MARS_R_KM * 1000.d0) > TINYVAL) then
          print *,'Error: invalid shell radius ',shell_radius,', should match ',R_shell(ilayer) * MARS_R_KM * 1000.d0
          print *,'  parameter file ',trim(filename)
          print *,'  All SH mars model parameter files must have same shell layer radii, please check with density file.'
          call exit_MPI(0,'Invalid parameter file shell radius')
        endif
      endif

      ! checks that shell layering must be given in increasing radius order
      if (ilayer > 1) then
        if (r_norm < R_shell(ilayer-1)) then
          print *,'Error: invalid shell radius ',shell_radius,', should be bigger than ', &
                   R_shell(ilayer-1) * MARS_R_KM * 1000.d0
          print *,'  parameter file ',trim(filename)
          call exit_MPI(0,'Invalid parameter file shell radius not in increasing order')
        endif
      endif

      ! reads in coefficients
      ! format:
      !0 0 4.041686e+03 0.000000e+00
      !1 0 -1.450610e+00 0.000000e+00
      !1 1 -6.312329e-01 9.913989e-02
      !2 0 4.808198e-01 0.000000e+00
      !..

      ! reads spherical harmonic coefficients
      do n = 0,degree_N         ! spherical harmonic degree n
        do m = 0,n              ! order m
          ! format: #n #m #A #B
          read(IIN,*) n_read,m_read,A_coeff,B_coeff

          ! checks if n,m match
          if (n_read /= n .or. m_read /= m) then
            print *,'Error: invalid (n,m) on coefficient line ',n_read,m_read,A_coeff,B_coeff
            print *,'  (n,m) should be ',n,m
            call exit_MPI(0,'Invalid coefficient line')
          endif

          ! stores coefficients
          !mn2Index(m,n)
          idx = n * (n + 1) / 2 + m

          A_shell(ipar,idx+1,ilayer) = A_coeff  ! radial
          B_shell(ipar,idx+1,ilayer) = B_coeff ! spherical
        enddo
      enddo
    enddo ! ilayer

    ! closes parameter file
    close(IIN)

  enddo ! ipar

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) '  mars SH model read:'
  write(IMAIN,*) '    number of spherical harmonic degrees: N        = ',degree_N
  write(IMAIN,*) '    number of shells                    : N_shells = ',N_shells
  write(IMAIN,*)
  write(IMAIN,*) '    number of spherical harmonic coefficients (per shell): coeffs  = ',num_sh_coeffs
  write(IMAIN,*)
  write(IMAIN,*) '    shell radius range (in km) = ',sngl(R_shell(1) * MARS_R_KM),'/', &
                                                     sngl(R_shell(N_shells) * MARS_R_KM)
  write(IMAIN,*) '    shell depth range  (in km) = ',sngl((1.d0-R_shell(1)) * MARS_R_KM),'/', &
                                                     sngl((1.d0-R_shell(N_shells)) * MARS_R_KM)
  write(IMAIN,*)
  call flush_IMAIN()

  ! user output
  write(IMAIN,*) '  setting up Legendre polynomial factors'
  write(IMAIN,*)
  call flush_IMAIN()

  ! allocates spherical harmonics factors
  ! Legendre polynomials normalization factors
  !long double *nF = malloc (sizeof (long double[num_sh_coeffs]));
  allocate(nFactors(0:num_sh_coeffs-1),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating nFactors array')
  nFactors(:) = 0.d0

  ! note: Legendre polynomials can differ by a normalization factor 1/(4 PI)
  !       while the bkmns model uses this normalization factor in routine get_nml_Factors(),
  !       here we use polynomials without this normalization.
  !
  !nmlFactors (N, nlg, nF);
  ! (similar to get_nml_Factors() routine, but without the sqrt(1/(4 PI) normalization, implemented in model_bkmns.f90)
  call get_nml_Factors_without_4PI_norm(degree_N, num_sh_coeffs, nFactors)

  ! Legendre degree factors
  allocate(Legendre_degree_factor(0:degree_N),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating Legendre factor array')
  Legendre_degree_factor(:) = 0.d0

  ! pre-computes Legendre polynomial factors
  ! (using same routine as in model_bkmns.f90)
  call get_Legendre_degree_Factors(degree_N,Legendre_degree_factor)

  end subroutine read_SH_mars_model


!
!-----------------------------------------------------------------------------------------
!

  subroutine model_SH_mars_crust(lat,lon,r,vpc,vsc,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)

  ! for crustal velocities

  use constants, only: ZERO,ONE,TWO_PI,PI_OVER_TWO,DEGREES_TO_RADIANS
  use shared_parameters, only: R_PLANET
  use model_sh_mars_par

  implicit none

  double precision,intent(in) :: lat,lon,r
  double precision,intent(out) :: vpc,vsc,rhoc,moho,sediment
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust,moho_only

  ! local parameters
  double precision :: vpv,vph,vsv,vsh,eta,rho
  double precision :: theta,phi

  ! initializes
  vpv = ZERO
  vph = ZERO
  vsv = ZERO
  vsh = ZERO
  rho = ZERO
  eta = ONE

  moho = ZERO
  sediment = ZERO
  found_crust = .true.

  ! determine moho based on 1D moho from case65TAY model
  ! RMOHO is set at 3280000.d0, that is at a depth of ~110 km.
  ! moho thickness
  moho = 110.d0
  ! non-dimensionalize thicknesses (given in km)
  moho = moho * 1000.d0 / R_PLANET

  ! checks if anything further to do
  if (moho_only) return

  ! lat/lon range: [-90,90] / [-180,180]
  ! lat/lon in degrees (range lat/lon = [-90,90] / [-180,180]
  !lat = (PI_OVER_TWO - theta) * RADIANS_TO_DEGREES
  !lon = phi * RADIANS_TO_DEGREES
  !if (lon > 180.0d0 ) lon = lon - 360.0d0

  ! theta/phi  - colatitude/longitude in rad (range theta/phi = [0,pi] / [0,2pi]
  theta = PI_OVER_TWO - lat * DEGREES_TO_RADIANS
  phi = lon * DEGREES_TO_RADIANS
  if (phi < 0.0d0 ) phi = phi + TWO_PI

  ! determine crustal velocities from SH model (defined for both crust & mantle, between CMB and surface)
  call model_SH_mars_crustmantle(r,theta,phi,vpv,vph,vsv,vsh,eta,rho)

  if (r > ONE - moho .or. elem_in_crust) then
    ! returns isotropic values
    rhoc = rho
    vpc = vpv   ! vpv == vph as crustmantle model is isotropic
    vsc = vsv
  else
    ! note: if x is exactly the moho depth this will return false
    found_crust = .false.
  endif

  end subroutine model_SH_mars_crust
!
!-----------------------------------------------------------------------------------------
!

  subroutine model_SH_mars_crustmantle(radius_in,theta_in,phi_in,vpv_out,vph_out,vsv_out,vsh_out,eta_out,rho_out)

  use constants, only: PI,PI_OVER_TWO,TWO_PI,EARTH_R_KM,GRAV
  use shared_parameters, only: R_PLANET,RHOAV

  use model_sh_mars_par

  implicit none

  ! radius     - normalized by globe radius [0,1.x]
  ! theta/phi  - colatitude/longitude in rad (range theta/phi = [0,pi] / [0,2pi]
  double precision,intent(in) :: radius_in,theta_in,phi_in

  ! absolute values, not perturbations
  double precision,intent(inout) :: vpv_out,vph_out,vsv_out,vsh_out,eta_out,rho_out

  ! local parameters
  double precision :: vp,vs,rho
  double precision :: radius,theta,phi

  integer :: n,m,idx
  integer :: ipar,ilayer
  integer :: ishell,ishell_lower,ishell_upper

  double precision :: x
  double precision :: c_interpolate
  double precision :: Cmn,Smn
  double precision :: A_coeff,B_coeff,A_coeff_lower,B_coeff_lower,A_coeff_upper,B_coeff_upper
  double precision, dimension(0:degree_N) :: C_coeffs,S_coeffs
  double precision, dimension(0:degree_N+1) :: P,Pnormalized

  ! mantle parameters (1==rho, 2==vp, 3==vs)
  double precision, dimension(NPAR) :: M_par
  double precision :: scaleval_vel,scaleval_rho

  ! takes position
  radius = radius_in

  ! note: theta,phi values here are outputs from reduce(theta,phi), which limits theta to [0,PI] and phi to [0,2PI].
  !        however, below we expect theta to go from south to north pole, and phi shifted by (PI/2 + PI)
  !        to have its meridian at PI/2 (instead of 0) and on the opposite hemisphere.
  !
  ! SH mars files use input definition:
  !  - colatitude: south pole at 0, north pole at PI
  !  - longitude : shifted by +PI/2
  theta = PI - theta_in              ! switches theta to have 0 at south pole
  phi = phi_in + PI_OVER_TWO + PI    ! shifts meridian & hemisphere)

  ! initializes parameter values
  vp = 0.d0
  vs = 0.d0
  rho = 0.d0
  M_par(:) = 0.d0

  ! cos(theta) needed for polar basis
  x = cos(theta)

  ! determines radial shell index and interpolation coefficient
  ! note: for position r between two shells, we interpolate between the two shell values
  ! (r2Index (NS, r1, r2, r, R[el][i][j][k], &c); from sh2gll.c)
  c_interpolate = 0.d0
  ishell = 0

  if (radius < R_shell(1)) then
    ishell = 1
    c_interpolate = 0.d0
  else if (radius >= R_shell(N_shells)) then
    ishell = N_shells
    c_interpolate = 1.d0
  else
    do ilayer = 1,N_shells-1
      if (radius >= R_shell(ilayer) .and. radius < R_shell(ilayer+1)) then
        ! interpolation factor between shells
        c_interpolate = (radius - R_shell(ilayer)) / (R_shell(ilayer+1) - R_shell(ilayer))
        ishell = ilayer
        ! done searching shell layer
        exit
      endif
    enddo
  endif

  ! checks ishell
  if (ishell <= 0 .or. ishell > N_shells) call exit_MPI(0,'Error invalid ishell index')

  ! sets lower/upper shell index for interpolation
  ishell_lower = ishell
  ishell_upper = ishell + 1
  if (ishell_upper > N_shells) ishell_upper = N_shells

  !debug
  !if (ishell < N_shells) print *,'debug: ishell ',ishell,N_shells,'radius',radius,R_shell(ishell_lower),R_shell(ishell_upper)

  ! determines model values based on spherical harmonics expansion
  ! (same as in model_bkmns.f90 mantle routine)

  ! initialize Legendre polynomials and coefficients
  P(:) = 0.d0
  Pnormalized(:) = 0.d0
  C_coeffs(:) = 0.d0
  S_coeffs(:) = 0.d0

  ! azimutal Basis
  call get_azimuthal_Basis(phi,degree_N,degree_N,C_coeffs,S_coeffs)

  do n = 0,degree_N
    ! polar Basis
    call get_polar_Basis(x,n,num_sh_coeffs,nFactors,degree_N,P,Pnormalized,Legendre_degree_factor)

    do m = 0,n
      Cmn = Pnormalized(m + 1) * C_coeffs(m)
      Smn = Pnormalized(m + 1) * S_coeffs(m)

      !mn2Index(m,n)
      idx = n * (n + 1) / 2 + m + 1   ! + 1 due to arrays starting at 1

      ! loops over parameters (rho,vp,vs)
      do ipar = 1,3
        ! lower shell coefficients
        A_coeff_lower = A_shell(ipar,idx,ishell_lower)
        B_coeff_lower = B_shell(ipar,idx,ishell_lower)
        ! upper shell coefficients
        A_coeff_upper = A_shell(ipar,idx,ishell_upper)
        B_coeff_upper = B_shell(ipar,idx,ishell_upper)
        ! interpolates between shells
        A_coeff = (1.d0 - c_interpolate) * A_coeff_lower + c_interpolate * A_coeff_upper
        B_coeff = (1.d0 - c_interpolate) * B_coeff_lower + c_interpolate * B_coeff_upper
        ! sets model parameter value
        M_par(ipar) = M_par(ipar) + (Cmn * A_coeff + Smn * B_coeff)
      enddo
    enddo
    !debug
    !if (n == 0) print *,'debug: rho/vp/vs=',M_par(1),M_par(2),M_par(3),' Cmn/Smn=',C_coeffs(0),S_coeffs(0), &
    !                    ' A_coeff/B_coeff=',A_coeff,B_coeff, &
    !                    ' Pnormalized ',Pnormalized(0),Pnormalized(1)
  enddo

  ! non-dimensionalize
  scaleval_rho = 1.0d0 / RHOAV                             ! from kg/m3
  scaleval_vel = 1.0d0 / (R_PLANET * sqrt(PI*GRAV*RHOAV))  ! from m/s

  ! sets corresponding parameter value
  rho = M_par(1) * scaleval_rho
  vp = M_par(2) * scaleval_vel
  vs = M_par(3) * scaleval_vel

  ! returns crust/mantle values if non-zero
  if (vp > WATER_LEVEL) then
    ! converts isotropic values to transverse isotropic
    vpv_out = vp
    vph_out = vp
    vsv_out = vs
    vsh_out = vs
    eta_out = 1.d0
    rho_out = rho
  endif

  end subroutine model_SH_mars_crustmantle

