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


! von Karman perturbation distribution
!
! apply to velocity models for investigating scattering effects

  module model_scattering_par

  use constants, only: CUSTOM_REAL
  use shared_parameters, only: SCATTERING_STRENGTH,SCATTERING_CORRELATION
  implicit none

  ! apply perturbation to Earth regions
  logical, parameter :: USE_CRUST_MANTLE_SCATTERING = .true.
  logical, parameter :: USE_INNER_CORE_SCATTERING = .false.
  logical, parameter :: USE_OUTER_CORE_SCATTERING = .false.

  ! perturbation array (regular x/y/z grid)
  real(kind=CUSTOM_REAL), dimension(:,:,:),allocatable :: perturbation_grid

  ! SPECFEM3D_GLOBE: normalized radius in range [-1,1] (for spherical mesh)
  double precision, parameter :: grid_length = 2.d0        ! model dimension [-1,1]

  ! for mesh interpolations
  ! grid origin location (-1), grid spacing
  double precision :: grid_origin,grid_dx
  ! fft size
  integer :: grid_N

  ! note: not working yet, mesh point positions are determined on the fly and not yet known at the onset
  !       of velocity determination get_model() routine - todo for future use: shrink fft grid size per process
  ! indexing of partial grid to cover single process slice
  !integer :: part_grid_imin,part_grid_imax,part_grid_jmin,part_grid_jmax,part_grid_kmin,part_grid_kmax
  !integer :: part_grid_size_i,part_grid_size_j,part_grid_size_k

  ! special functions
  ! log2 interface
  interface log2
    module procedure rlog2,dlog2
  end interface

  contains

    !---------------------------------------------------
    ! logarithm base 2 (single precision)
    real function rlog2(x)

    implicit none
    real, intent(in) :: x

    rlog2 = log(x) / log(2.0)

    end function

    !---------------------------------------------------
    ! logarithm base 2 (double precision)
    double precision function dlog2(x)

    implicit none
    double precision, intent(in) :: x

    dlog2 = log(x) / log(2.d0)

    end function

  end module model_scattering_par

!
!--------------------------------------------------------------------------------------------------
!

  subroutine model_scattering_broadcast()

  ! standard routine to setup model

  use constants
  use model_scattering_par

  implicit none

  ! user info
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'adding scattering perturbations:'
    call flush_IMAIN()
  endif

  ! generate perturbation distribution
  call generate_perturbations()

  end subroutine model_scattering_broadcast

!
!--------------------------------------------------------------------------------------------------
!


  real function psd_vonKarman_3D(a_in,kx,ky,kz)

  ! von Karman distribution
  !
  ! 3D power spectral density: Imperatori & Mai (2012)
  !
  ! P(k) = sigma**2 ( 2 sqrt(pi * a) )**E gamma(H + E/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + E/2) )
  !      = sigma**2 ( 2 sqrt(pi * a) )**3 gamma(H + 3/2) / ( gamma(H) * (1 + k**2 a**2 )**(H + 3/2) )
  !
  ! E = euclidian dimension -> E = 3 for 3D medium
  !
  !from scipy.special import gamma

  use constants, only: CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: a_in,kx,ky,kz

  ! local parameters
  real(kind=CUSTOM_REAL) :: a,H,sigma
  real(kind=CUSTOM_REAL) :: g1,g2,amp,k,psd

  real(kind=CUSTOM_REAL), parameter :: PI = real(acos(-1.d0),kind=CUSTOM_REAL)
  real(kind=CUSTOM_REAL), parameter :: CONST_3HALF = real(3.d0 / 2.d0,kind=CUSTOM_REAL)

  a = a_in    ! correlation length
  H = 0.3     ! Hurst exponent: Imperatori & Mai use H = 0.1 and 0.3
  sigma = 0.1 ! standard deviation: 10%

  ! gamma function - a standard intrinsic function for Fortran 2008 and later
  !
  ! coefficients
  ! see: https://en.wikipedia.org/wiki/Particular_values_of_the_gamma_function
  ! with gamma(1/2) = sqrt(pi)
  !      gamma(1)   = 1
  !      gamma(3/2) = sqrt(pi)/2
  !
  g1 = gamma(H + CONST_3HALF)                            ! gamma function: gamma(H+3/2)
  g2 = gamma(H)                                          ! gamma function: gamma(H)
  amp = sigma * sigma * ( 2.0 * sqrt(PI * a) )**3        ! coefficient   : sigma**2 * ( 2 * sqrt(pi * a) )**3

  ! wavenumber k = sqrt(kx**2 + ky**2 + kz**2)
  k = sqrt(kx*kx + ky*ky + kz*kz)

  ! 3D von Karman power spectral density
  ! P(k) = sigma**2 ( 2 sqrt(pi * a) )**3 *  gamma(H + 3/2) / (gamma(H) * (1 + k**2 a**2 )**(H + 3/2)
  !      =             amp                *        g1       / (   g2    * (1 + k**2 a**2 )**(H + 3/2)
  psd = amp * g1 / ( g2 * (1.0 + k**2 * a**2)**(H + CONST_3HALF) )

  ! return value
  psd_vonKarman_3D = psd
  return

  end function

!
!--------------------------------------------------------------------------------------------------
!

  subroutine generate_perturbations()

  use constants
  use shared_parameters, only: R_PLANET_KM,T_min_period,estimated_min_wavelength,SAVE_MESH_FILES
  use model_scattering_par

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: lambda_min,lambda_min_norm
  real(kind=CUSTOM_REAL) :: dx,length,min_length,a_corr
  real(kind=CUSTOM_REAL) :: dk,k_min,k_max,k_lambda,kx,ky,kz
  real(kind=CUSTOM_REAL) :: psd,rand_phase
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: freqs

  integer :: N,npower_of_2,index_k_lambda
  integer :: i,j,k,ik,ier

  integer, parameter :: CUSTOM_CMPLX = 8
  complex(kind=CUSTOM_CMPLX),dimension(:,:,:),allocatable :: kxyz_dist
  complex(kind=CUSTOM_CMPLX) :: k_random

  ! helper array for FFT
  real(kind=CUSTOM_REAL) :: mpow(30)
  real(kind=CUSTOM_REAL) :: val_avg,val_max

  ! random seed
  integer :: num_seed
  integer,dimension(:),allocatable :: myseed

  ! debug output
  character(len=MAX_STRING_LEN) :: name

  ! external function
  real,external :: psd_vonKarman_3D

  ! timing
  double precision, external :: wtime
  double precision :: time_start,tCPU

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  perturbation distribution       : ','von Karman'
    write(IMAIN,*) '  perturbation correlation factor : ',SCATTERING_CORRELATION
    write(IMAIN,*) '  perturbation maximum amplitude  : ',SCATTERING_STRENGTH
    write(IMAIN,*)
    write(IMAIN,*) '  planet radius                   : ',R_PLANET_KM,'(km)'
    write(IMAIN,*) '  estimated minimum period        : ',sngl(T_min_period),'(s)'
    write(IMAIN,*) '  estimated minimum wavelength    : ',sngl(estimated_min_wavelength),'(km)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

  ! estimated minimum wavelength resolved by mesh
  lambda_min = real(estimated_min_wavelength,kind=CUSTOM_REAL)

  ! normalized wavelength
  lambda_min_norm = real(lambda_min / R_PLANET_KM,kind=CUSTOM_REAL)

  ! quarter of wavelength for grid estimate (at least 5 grid points per wavelength)
  min_length = lambda_min_norm / 4.0_CUSTOM_REAL

  ! total length
  length = grid_length

  ! next power of 2 for FFT
  npower_of_2 = int(log2(length / min_length)) + 1
  N = 2**npower_of_2

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  estimated number of grid points : ',N
    call flush_IMAIN()
  endif

  ! limits minimum/maximum number of points for FFT
  !if (N < 1024) N = 1024  ! number of minimum points along x-direction (power of 2 for fft) 2**(10)
  !if (N < 512) N = 512    ! number of minimum points along x-direction (power of 2 for fft) 2**(9)
  if (N < 256) N = 256    ! number of minimum points along x-direction (power of 2 for fft) 2**(8)

  !if (N > 8192) N = 8192  ! number of maximum points along x-direction (power of 2 for fft) 2**(13)
  !if (N > 4096) N = 4096  ! number of maximum points along x-direction (power of 2 for fft) 2**(12)
  if (N > 2048) N = 2048  ! number of minimum points along x-direction (power of 2 for fft) 2**(11)

  ! power of 2
  npower_of_2 = ceiling(log2(real(N)))

  ! grid space increment
  dx = length / (N-1)   ! dx == dy == dz

  ! for mesh interpolations (uses double precision values)
  grid_origin  = - grid_length / 2.d0         ! origin at -1.0 for x,y,z positions
  grid_dx      = grid_length / (N-1)
  grid_N       = N

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  FFT using number of grid points : N           = ',N
    write(IMAIN,*) '  FFT using power of 2            : npower_of_2 = ',npower_of_2
    write(IMAIN,*) '  grid spacing                    : dx          = ',dx * R_PLANET_KM,'(km)'
    write(IMAIN,*) '  grid memory size                : ',dble(N) * dble(N) * dble(N) * dble(CUSTOM_REAL) / 1024.d0 / 1024.d0,'MB'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! wavenumbers
  ! min/max: k_max = 2 pi / (2 dx)
  k_min = real(2.d0 * PI / (N * dx),kind=CUSTOM_REAL)
  k_max = real(2.d0 * PI / (2.0 * dx),kind=CUSTOM_REAL)

  ! wavenumber increment
  dk = real(2.d0 * PI / (N * dx),kind=CUSTOM_REAL)

  ! maximum index for k_max
  !index_k_max = N / 2 # k_max / dk = (2 pi / (2 dx) ) / (2 pi / (N dx) ) = (N * dx) / (2 * dx) = N / 2

  ! index for k_lambda_min
  k_lambda = real(2.d0 * PI / (lambda_min_norm),kind=CUSTOM_REAL)
  index_k_lambda = int(k_lambda / dk)

  ! correlation length
  ! such that k * a ~ 1 (for correlation factor == 1)
  a_corr = real(lambda_min_norm / (2.d0 * PI) * SCATTERING_CORRELATION,kind=CUSTOM_REAL)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  wavenumbers: k min/max          = ',k_min,'/',k_max
    write(IMAIN,*) '               dk increment       = ',dk
    write(IMAIN,*) '               k lambda_min       = ',k_lambda
    write(IMAIN,*) '               index k lambda_min = ',index_k_lambda
    write(IMAIN,*)
    write(IMAIN,*) '  correlation length              : ',a_corr,' - in km: ',a_corr * R_PLANET_KM
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! setup wavenumber arrays
  if (myrank == 0) then
    ! only main process creates pertubation grid, then distributes it to all others
    ! to have the same perturbation grid for all processes.

    ! 3D arrays
    !allocate(freqs(N), wavenumbers_kx(N,N,N), wavenumbers_ky(N,N,N), wavenumbers_kz(N,N,N),stat=ier)
    allocate(freqs(N),stat=ier)
    if (ier /= 0) stop 'Error allocating wavenumbers arrays'

    ! fft indexing
    freqs(:) = 0.0
    do i = 1,N
      ! numpy-like
      !if (i <= N/2+1) then
      !  ! indices: 0,1,..,N/2
      !  ik = i - 1
      !else
      !  ! indices: -N/2-1,-N/2-2,..,-1
      !  ik = - (N - (i - 1))
      !endif

      ! here positive values only; will take conjugates later for second half of array (when calling to apply symmetries)
      if (i <= N/2+1) then
        ! indices: 0,1,..,N/2
        ik = i - 1
      else
        ! indices: N/2-1,N/2-2,..,1
        ik = N - (i - 1)
      endif
      freqs(i) = real(ik,kind=CUSTOM_REAL)     ! for d = 1.0/N: freqs = ik / (d*N) = ik / ((1.0 / N) * N) = ik / ( 1.0 ) = ik
    enddo
    ! debug
    !print *,'debug: freqs ',freqs(:); print *
    ! same as:
    !freqs = (/0.0,(real(i,kind=CUSTOM_REAL),i=1,N/2)/)
    !freqs = (/freqs(1:N/2+1),(real(i,kind=CUSTOM_REAL),i=N/2-1,1,-1)/)
    !print *,'debug: new freqs ',freqs(:); print *
    ! debug
    !print *,'maximum wavenumber = ',maxval(freqs(:)) * dk ! maximum wavenumber

    ! sets up helper array for FFTs
    do i = 1,npower_of_2
      mpow(i) = 2**(npower_of_2-i)
    enddo

    ! initializes random number generator
    ! seed with fixed value to make successive runs repeatable
    call random_seed(size=num_seed)
    allocate(myseed(num_seed))
    myseed(1:num_seed) = 12345
    call random_seed(put=myseed)

    ! debug
    !call random_number(rand_phase)
    !print *,'debug: rank ',myrank,'random number 1: ',rand_phase

    ! perturbation grid array
    allocate(perturbation_grid(N,N,N),stat=ier)
    if (ier /= 0) stop 'Error allocating perturbation_grid array'
    perturbation_grid(:,:,:) = 0.0_CUSTOM_REAL

    ! wavenumber distribution
    allocate(kxyz_dist(N,N,N),stat=ier)
    if (ier /= 0) stop 'Error allocating kxyz_dist array'
    kxyz_dist(:,:,:) = cmplx(0.0,0.0)

    ! applies amplitudes, which follow defined power spectral density (psd), to random phases
    do k = 1,N
      do j = 1,N
        do i = 1,N
          ! wavenumbers
          kx = freqs(i) * dk  ! (2.0 * np.pi * freqs[icol]) / L
          ky = freqs(j) * dk  ! (2.0 * np.pi * freqs[irow]) / L
          kz = freqs(k) * dk  ! (2.0 * np.pi * freqs[iz]) / L

          ! amplitudes w/ von Karman distribution
          psd = psd_vonKarman_3D(a_corr,kx,ky,kz)

          ! random phase
          call random_number(rand_phase)
          ! range [0,2pi]
          rand_phase = real(rand_phase * 2.d0 * PI,kind=CUSTOM_REAL)
          k_random = cmplx( cos(rand_phase), sin(rand_phase) )

          ! stores wavenumber distribution
          kxyz_dist(i,j,k) = k_random * sqrt(psd)
        enddo
      enddo
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  starting 3D FFTs'
      call flush_IMAIN()
    endif

    ! define symmetry conditions for 3D FFT
    call fft_apply_3D_symmetry(kxyz_dist,N)

    ! FFT arrays
    ! example: 1D fft
    !allocate(k_line(N),x_FFT(N),stat=ier)
    !if (ier /= 0) stop 'Error allocating x_FFT array'
    !k_line(:) = cmplx(0.0,0.0)
    !x_FFT(:) = 0.0_CUSTOM_REAL
    ! inverse Fourier transform
    ! w/ 1D FFTs
    !do k = 1,N
    !  do j = 1,N
    !    ! takes 1D line
    !    k_line(:) = kxyz_dist(:,j,k)
    !    ! pad negative k
    !    do ii = 2, N/2
    !      ! fills from N,N-1,..,N/2+2
    !      k_line(N+2-ii) = conjg(k_line(ii))
    !    enddo
    !    ! 1D inverse FFT
    !    call FFTinv(npower_of_2, k_line(:), 1.0_CUSTOM_REAL, dk, x_FFT(:), mpow) ! inverse FFT, outputs real array x_FFT
    !
    !    !call rspec(k_line(:),N/2)                                    ! restructuring
    !    !call FFT(npower_of_2, k_line(:), -1.0_CUSTOM_REAL, dk, mpow) ! inverse FFT, outputs complex array k_line
    !    !x_FFT(1:N) = real(k_line(1:N))                               ! takes the real part
    !
    !    ! stores perturbations array
    !    perturbation_grid(:,j,k) = x_FFT(:)
    !  enddo
    !enddo

    ! 3D FFT
    call FFT_3D(N, npower_of_2, kxyz_dist, -1.0_CUSTOM_REAL, dk, mpow) ! inverse 3D FFT

    ! stores real part
    do k = 1,N
      do j = 1,N
        do i = 1,N
          ! stores perturbations array
          perturbation_grid(i,j,k) = real(kxyz_dist(i,j,k),kind=CUSTOM_REAL)
        enddo
      enddo
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  applying scaling: maximum amplitude = ',SCATTERING_STRENGTH
      call flush_IMAIN()
    endif

    ! debug
    !print *,'before scattering perturbation: min/max = ',minval(perturbation_grid),'/',maxval(perturbation_grid)
    !print *,'                                average = ',sum(perturbation_grid)/(N*N*N)

    ! applies scaling
    ! makes sure it has a zero average
    val_avg = sum(perturbation_grid) / (N*N*N)
    perturbation_grid(:,:,:) = perturbation_grid(:,:,:) - val_avg

    ! normalizes to range [-1,1]
    val_max = maxval(abs(perturbation_grid))
    if (val_max > 0.0_CUSTOM_REAL) then
      perturbation_grid(:,:,:) = perturbation_grid(:,:,:) / val_max
    endif

    ! scales with maximum strength
    perturbation_grid(:,:,:) = real(perturbation_grid(:,:,:) * SCATTERING_STRENGTH,kind=CUSTOM_REAL)

    ! debug
    !print *,'scattering perturbation: min/max = ',minval(perturbation_grid),'/',maxval(perturbation_grid)
    !print *,'                         average = ',sum(perturbation_grid)/(N*N*N)
    !print *,'debug: loop done'

    ! free memory
    deallocate(freqs)
    deallocate(kxyz_dist)
  endif

  call synchronize_all()

  ! allocates grid on secondary processes
  if (myrank /= 0) then
    ! perturbation grid array
    allocate(perturbation_grid(N,N,N),stat=ier)
    if (ier /= 0) stop 'Error allocating perturbation_grid array'
    perturbation_grid(:,:,:) = 0.0_CUSTOM_REAL
  endif

  ! main process distributes grid to all other arrays
  call bcast_all_cr(perturbation_grid,N*N*N)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  perturbations: min/max = ',minval(perturbation_grid),'/',maxval(perturbation_grid)
    write(IMAIN,*) '                 average = ',sum(perturbation_grid)/(N*N*N)
    write(IMAIN,*)
    ! timing
    tCPU = wtime() - time_start
    write(IMAIN,*) '  Elapsed time for perturbation setup in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! visualization output
  if (SAVE_MESH_FILES .and. myrank == 0) then
    ! output VTU file for visual inspection
    name = 'perturbation_grid'
    call plot_grid_data(perturbation_grid,N,length,dx,name)
  endif

  end subroutine generate_perturbations

!
!-------------------------------------------------------------------------------------------------
!

  subroutine plot_grid_data(array,N,length,dx,name)

  use constants, only: myrank,IMAIN,CUSTOM_REAL,MAX_STRING_LEN

  implicit none

  integer, intent(in) :: N
  real(kind=CUSTOM_REAL),dimension(N,N,N),intent(in) :: array

  real(kind=CUSTOM_REAL),intent(in) :: length,dx
  character(len=MAX_STRING_LEN) :: name

  ! local parameters
  integer :: ne,np,ixyz,n1,n2,n3,n4,n5,n6,n7,n8
  integer :: i,j,k,ier

  ! global point data
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: total_dat
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: total_dat_xyz
  integer,dimension(:,:),allocatable :: total_dat_con
  character(len=MAX_STRING_LEN) :: mesh_file,var_name

  ! regular grid
  np = N * N * N             ! total number of points
  ne = (N-1) * (N-1) * (N-1) ! total number of elements

  ! creates array to hold point data
  allocate(total_dat(np),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat array'
  total_dat(:) = 0.0_CUSTOM_REAL
  allocate(total_dat_xyz(3,np),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_xyz array'
  total_dat_xyz(:,:) = 0.0_CUSTOM_REAL
  allocate(total_dat_con(8,ne),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_con array'
  total_dat_con(:,:) = 0

  ! regular grid points values
  np = 0
  do k = 1,N
    do j = 1,N
      do i = 1,N
        ! number counter
        np = np + 1

        ! checks point index
        ixyz = i + (j-1) * (N) + (k-1) * (N * N)
        if (np /= ixyz) stop 'Invalid grid point'

        ! array value
        total_dat(np) = array(i,j,k)

        ! grid point position [-L/2,L/2]
        total_dat_xyz(1,np) = -length/2.0 + (i-1) * dx
        total_dat_xyz(2,np) = -length/2.0 + (j-1) * dx
        total_dat_xyz(3,np) = -length/2.0 + (k-1) * dx
      enddo
    enddo
  enddo
  ! element connections
  ne = 0
  do k = 1,N-1
    do j = 1,N-1
      do i = 1,N-1
        ! element counter
        ne = ne + 1

        ! point index
        !ixyz = i + (j-1) * (N) + (k-1) * (N * N)

        ! element corners
        n1 = i     + (j-1) * (N) + (k-1) * (N * N)  ! i  ,j  ,k
        n2 = (i+1) + (j-1) * (N) + (k-1) * (N * N)  ! i+1,j  ,k
        n3 = (i+1) + (j  ) * (N) + (k-1) * (N * N)  ! i+1,j+1,k
        n4 = i     + (j  ) * (N) + (k-1) * (N * N)  ! i  ,j+1,k
        n5 = i     + (j-1) * (N) + (k  ) * (N * N)  ! i  ,j  ,k+1
        n6 = (i+1) + (j-1) * (N) + (k  ) * (N * N)  ! i+1,j  ,k+1
        n7 = (i+1) + (j  ) * (N) + (k  ) * (N * N)  ! i+1,j+1,k+1
        n8 = i     + (j  ) * (N) + (k  ) * (N * N)  ! i  ,j+1,k+1

        ! note: indices for VTK start at 0
        total_dat_con(1,ne) = n1 - 1
        total_dat_con(2,ne) = n2 - 1
        total_dat_con(3,ne) = n3 - 1
        total_dat_con(4,ne) = n4 - 1
        total_dat_con(5,ne) = n5 - 1
        total_dat_con(6,ne) = n6 - 1
        total_dat_con(7,ne) = n7 - 1
        total_dat_con(8,ne) = n8 - 1
      enddo
    enddo
  enddo

  ! VTU binary format
  mesh_file = 'OUTPUT_FILES/' // trim(name) // '.vtu'
  var_name = 'val'
  call write_VTU_movie_data_binary(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) '  VTU file written: ',trim(mesh_file)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! free memory
  deallocate(total_dat,total_dat_xyz,total_dat_con)

  end subroutine plot_grid_data


!
!-------------------------------------------------------------------------------------------------
!

  subroutine model_scattering_add_perturbations(iregion_code,xmesh,ymesh,zmesh, &
                                                vpv,vph,vsv,vsh,rho,eta_aniso, &
                                                c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                                c33,c34,c35,c36,c44,c45,c46,c55,c56,c66)

  ! overwrites values with updated model values (from iteration step) here, given at all GLL points

  use constants, only: myrank,IREGION_CRUST_MANTLE,IREGION_INNER_CORE,IREGION_OUTER_CORE
  use meshfem_models_par, only: ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE
  use model_scattering_par

  implicit none

  integer,intent(in) :: iregion_code
  double precision,intent(in) :: xmesh,ymesh,zmesh

  double precision,intent(inout) :: vpv,vph,vsv,vsh,rho,eta_aniso
  double precision,intent(inout) :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                                    c33,c34,c35,c36,c44,c45,c46,c55,c56,c66

  ! local parameters
  double precision :: pert_val,scaling
  integer :: ix,iy,iz,N
  ! positioning
  double precision :: spac_x,spac_y,spac_z
  double precision :: gamma_interp_x,gamma_interp_y,gamma_interp_z
  double precision :: val1,val2,val3,val4,val5,val6,val7,val8

  ! checks if anything to do
  if ((.not. USE_CRUST_MANTLE_SCATTERING) .and. (.not. USE_OUTER_CORE_SCATTERING) .and. (.not. USE_INNER_CORE_SCATTERING)) return

  ! determine spacing and cell for linear interpolation
  spac_x = (xmesh - grid_origin) / grid_dx
  spac_y = (ymesh - grid_origin) / grid_dx
  spac_z = (zmesh - grid_origin) / grid_dx

  ix = int(spac_x)
  iy = int(spac_y)
  iz = int(spac_z)

  gamma_interp_x = spac_x - dble(ix)
  gamma_interp_y = spac_y - dble(iy)
  gamma_interp_z = spac_z - dble(iz)

  ! shifting indices to start at 1
  ix = ix + 1
  iy = iy + 1
  iz = iz + 1

  ! number of grid points (in x,y,z direction)
  N = grid_N

  ! suppress edge effects for points outside of the model SPOSTARE DOPO
  if (ix < 1) then
     ix = 1
     gamma_interp_x = 0.d0
  endif
  if (ix > N-1) then
     ix = N-1
     gamma_interp_x = 1.d0
  endif

  if (iy < 1) then
     iy = 1
     gamma_interp_y = 0.d0
  endif
  if (iy > N-1) then
     iy = N-1
     gamma_interp_y = 1.d0
  endif

  if (iz < 1) then
     iz = 1
     gamma_interp_z = 0.d0
  endif
  if (iz > N-1) then
     iz = N-1
     gamma_interp_z = 1.d0
  endif

  ! checks
  if (ix < 0 .or. iy < 0 .or. iz < 0) then
    print *,'Error: scattering position has invalid index: '
    print *,'  rank        : ',myrank
    print *,'  corner index: ',ix,iy,iz
    print *,'  location    : ',sngl(xmesh),sngl(ymesh),sngl(zmesh)
    print *,'  origin      : ',sngl(grid_origin)
    call exit_MPI(myrank,'Error corner index in model_scattering routine')
  endif

  ! perturbation values
  val1 = perturbation_grid(ix  ,iy  ,iz  )
  val2 = perturbation_grid(ix+1,iy  ,iz  )
  val3 = perturbation_grid(ix+1,iy+1,iz  )
  val4 = perturbation_grid(ix  ,iy+1,iz  )
  val5 = perturbation_grid(ix  ,iy  ,iz+1)
  val6 = perturbation_grid(ix+1,iy  ,iz+1)
  val7 = perturbation_grid(ix+1,iy+1,iz+1)
  val8 = perturbation_grid(ix  ,iy+1,iz+1)

  ! interpolation rule
  ! use trilinear interpolation in cell to define perturbation value
  pert_val =  &
       val1 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val2 * gamma_interp_x        * (1.d0-gamma_interp_y) * (1.d0-gamma_interp_z) + &
       val3 * gamma_interp_x        * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val4 * (1.d0-gamma_interp_x) * gamma_interp_y        * (1.d0-gamma_interp_z) + &
       val5 * (1.d0-gamma_interp_x) * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val6 * gamma_interp_x        * (1.d0-gamma_interp_y) * gamma_interp_z + &
       val7 * gamma_interp_x        * gamma_interp_y        * gamma_interp_z + &
       val8 * (1.d0-gamma_interp_x) * gamma_interp_y        * gamma_interp_z

  !debug
  !if (myrank == 0) &
  !  print *,'debug: loc ',sngl(xmesh),sngl(ymesh),sngl(zmesh),'ixyz',ix,iy,iz, &
  !          'perturbation ',pert_val,'corners',val1,val2,val3,val4,val5,val6,val7,val8

  ! adds perturbation
  scaling = 1.d0 + pert_val

  ! crust/mantle
  if (USE_CRUST_MANTLE_SCATTERING .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! define new elastic parameters in the model
    rho = scaling * rho
    vpv = scaling * vpv
    vph = scaling * vph
    vsv = scaling * vsv
    vsh = scaling * vsh
    eta_aniso = scaling * eta_aniso

    ! anisotropy
    if (ANISOTROPIC_3D_MANTLE) then
      c11 = scaling * c11
      c12 = scaling * c12
      c13 = scaling * c13
      c14 = scaling * c14
      c15 = scaling * c15
      c16 = scaling * c16
      c22 = scaling * c22
      c23 = scaling * c23
      c24 = scaling * c24
      c25 = scaling * c25
      c26 = scaling * c26
      c33 = scaling * c33
      c34 = scaling * c34
      c35 = scaling * c35
      c36 = scaling * c36
      c44 = scaling * c44
      c45 = scaling * c45
      c46 = scaling * c46
      c55 = scaling * c55
      c56 = scaling * c56
      c66 = scaling * c66
    endif
  endif ! crust/mantle

  ! outer core
  if (USE_OUTER_CORE_SCATTERING .and. iregion_code == IREGION_OUTER_CORE) then
    ! define new elastic parameters in the model
    rho = scaling * rho
    vpv = scaling * vpv
    vph = scaling * vph
    vsv = scaling * vsv
    vsh = scaling * vsh
    eta_aniso = scaling * eta_aniso
  endif ! outer core

  ! inner core
  if (USE_INNER_CORE_SCATTERING .and. iregion_code == IREGION_INNER_CORE) then
    ! define new elastic parameters in the model
    rho = scaling * rho
    vpv = scaling * vpv
    vph = scaling * vph
    vsv = scaling * vsv
    vsh = scaling * vsh
    eta_aniso = scaling * eta_aniso

    ! anisotropy
    if (ANISOTROPIC_INNER_CORE) then
      c11 = scaling * c11
      c12 = scaling * c12
      c13 = scaling * c13
      c33 = scaling * c33
      c44 = scaling * c44
    endif
  endif ! inner core

  end subroutine model_scattering_add_perturbations
