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

! XSMOOTH_SEM
!
! USAGE
!   mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR
!
!
! COMMAND LINE ARGUMENTS
!   SIGMA_H                - horizontal smoothing radius
!   SIGMA_V                - vertical smoothing radius
!   KERNEL_NAME            - kernel name, e.g. alpha_kernel
!   INPUT_DIR              - directory from which arrays are read
!   OUTPUT_DIR             - directory to which smoothed array are written
!
! DESCRIPTION
!   Reads kernels from INPUT_DIR, smooths by convolution with a Gaussian, and
!   write the resulting smoothed kernels to OUTPUT_DIR.
!
!   Files written to OUTPUT_DIR have the suffix 'smooth' appended,
!   e.g. proc***alpha_kernel.bin becomes proc***alpha_kernel_smooth.bin
!
!   This program's primary use case is to smooth kernels. It can be used though on
!   any "reg1" array, i.e. any array of dimension
!   (NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE).
!
!   This is a parrallel program -- it must be invoked with mpirun or other
!   appropriate utility.  Operations are performed in embarassingly-parallel
!   fashion.
!
!   Because this routine uses constants.h and values_from_mesher_globe.h,
!   you need to compile it for your specific case.
!
!   We switch between vectorized and non-vectorized version by using pre-processor
!   flag FORCE_VECTORIZATION and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK
!   defined in config.fh

#include "config.fh"

program smooth_sem_globe

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,IIN,IOUT, &
    GAUSSALPHA,GAUSSBETA,PI,TWO_PI,R_EARTH_KM,MAX_STRING_LEN,DEGREES_TO_RADIANS,NGLLCUBE,myrank

  use postprocess_par, only: &
    NCHUNKS_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,NPROCTOT_VAL,NEX_XI_VAL,NEX_ETA_VAL, &
    ANGULAR_WIDTH_XI_IN_DEGREES_VAL,ANGULAR_WIDTH_ETA_IN_DEGREES_VAL, &
    NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,MAX_KERNEL_NAMES,LOCAL_PATH

  use kdtree_search

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use adios_helpers_mod
  use manager_adios
#endif

  implicit none

  !-------------------------------------------------------------
  ! Parameters

  ! copy from static compilation (depends on Par_file values)
  integer, parameter :: NPROC_XI  = NPROC_XI_VAL
  integer, parameter :: NPROC_ETA = NPROC_ETA_VAL
  integer, parameter :: NCHUNKS   = NCHUNKS_VAL

  !takes region 1 kernels
  integer, parameter :: NSPEC_AB = NSPEC_CRUST_MANTLE
  integer, parameter :: NGLOB_AB = NGLOB_CRUST_MANTLE

  ! only include the neighboring 3 x 3 slices
  integer, parameter :: NSLICES = 3
  integer ,parameter :: NSLICES2 = NSLICES * NSLICES
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer, parameter :: NARGS = 6
  character(len=*), parameter :: reg_name = 'reg1/'
#else
  integer, parameter :: NARGS = 5
  character(len=*), parameter :: reg_name = '_reg1_'
#endif


  !-------------------------------------------------------------

  integer :: islice(NSLICES2), islice0(NSLICES2)

  integer, dimension(:,:,:,:), allocatable :: ibool
  integer, dimension(:), allocatable :: idoubling
  logical, dimension(:), allocatable :: ispec_is_tiso

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore, ystore, zstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: kernel, kernel_smooth,tk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: bk, jacobian
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xx0, yy0, zz0, xx, yy, zz

  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx0, cy0, cz0, cx, cy, cz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: exp_val
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: factor

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: norm, norm_h, norm_v, element_size_m
  real(kind=CUSTOM_REAL) :: x0, y0, z0
  real(kind=CUSTOM_REAL) :: center_x0, center_y0, center_z0
  real(kind=CUSTOM_REAL) :: center_x, center_y, center_z

  real(kind=CUSTOM_REAL),dimension(MAX_KERNEL_NAMES) :: min_old,max_old
  real(kind=CUSTOM_REAL) :: max_new, min_new
  real(kind=CUSTOM_REAL) :: max_old_all, max_new_all, min_old_all, min_new_all

  integer :: nspec, nglob
  integer :: ier,ichunk, ixi, ieta, iglob,nums
  integer :: ispec,iproc,ispec2,inum
  integer :: i,j,k
  integer :: sizeprocs

  character(len=MAX_STRING_LEN),dimension(NARGS) :: arg
  character(len=MAX_STRING_LEN),dimension(MAX_KERNEL_NAMES) :: kernel_names
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: kernel_name, topo_dir
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  character(len=MAX_STRING_LEN) :: input_file, solver_file
  character(len=MAX_STRING_LEN) :: varname
#else
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  character(len=MAX_STRING_LEN) :: local_data_file
#endif
  character(len=MAX_STRING_LEN) :: output_file
  character(len=MAX_STRING_LEN) :: filename

  integer :: nker,iker

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll

  ! array with all the weights in the cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: element_size
  real(kind=CUSTOM_REAL) :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! search elements
  integer :: ielem,num_elem_local

  ! tree nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: r_search,r_search_dist_v,r_search_dist_h

  ! debugging
  integer, dimension(:),allocatable :: ispec_flag
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tmp_bk

  ! timing
  double precision :: time_start
  double precision :: tCPU
  double precision, external :: wtime

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer(kind=8) :: group_size_inc,local_dim
  logical :: is_kernel
#endif

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#endif

  !------------------------------------------------------------------------------------------
  ! USER PARAMETERS

  ! number of steps to reach 100 percent, i.e. 10 outputs info for every 10 percent
  integer,parameter :: NSTEP_PERCENT_INFO = 10

  ! brute-force search for neighboring elements, if set to .false. a kd-tree search will be used
  logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

  ! kd-tree search uses either a spherical or ellipsoid search
  logical,parameter :: DO_SEARCH_ELLIP = .true.

  ! debugging
  logical, parameter :: DEBUG = .false.
  ! boundary point: 910, interior point: 2382, surface point: 5128
  integer ,parameter :: tmp_ispec_dbg = 5128

  !------------------------------------------------------------------------------------------

  ! ============ program starts here =====================

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! check command line arguments
  if (command_argument_count() /= NARGS) then
    if (myrank == 0) then
#ifdef USE_ADIOS_INSTEAD_OF_MESH
      print *,'Usage: mpirun -np NPROC bin/xsmooth_sem_adios SIGMA_H SIGMA_V KERNEL_NAME INPUT_FILE SOLVER_FILE OUTPUT_FILE'
      print *,'   with'
      print *,'     SIGMA_H, SIGMA_V - horizontal and vertical smoothing lenghts'
      print *,'     KERNEL_NAME      - comma-separated kernel names (e.g., alpha_kernel,beta_kernel)'
      print *,'     INPUT_FILE       - ADIOS file with kernel values (e.g., kernels.bp)'
      print *,'     SOLVER_FILE      - ADIOS file with mesh arrays (e.g., DATABASES_MPI/solver_data.bp)'
      print *,'     OUTPUT_FILE      - ADIOS file for smoothed output'
      print *
#else
      print *, 'Usage: mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR'
      print *,'   with'
      print *,'     SIGMA_H, SIGMA_V - horizontal and vertical smoothing lenghts'
      print *,'     KERNEL_NAME      - comma-separated kernel names (e.g., alpha_kernel,beta_kernel)'
      print *,'     INPUT_DIR        - directory with kernel files (e.g., proc***_alpha_kernel.bin)'
      print *,'     OUTPUT_DIR       - directory for smoothed output files'
      print *
#endif
      stop ' Please check command line arguments'
    endif
  endif
  call synchronize_all()

  ! check number of MPI processes
  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print *
      print *,'Expected number of MPI processes: ', NPROCTOT_VAL
      print *,'Actual number of MPI processes: ', sizeprocs
      print *
    endif
    call synchronize_all()
    stop 'Error wrong number of MPI processes'
  endif
  call synchronize_all()

  if (myrank == 0) then
    print *, 'Running XSMOOTH_SEM'
    print *
  endif
  call synchronize_all()

  ! parse command line arguments
  do i = 1, NARGS
    call get_command_argument(i,arg(i))
    if (i <= 5 .and. trim(arg(i)) == '') then
      stop ' Please check command line arguments'
    endif
  enddo
  read(arg(1),*) sigma_h
  read(arg(2),*) sigma_v
  kernel_names_comma_delimited = arg(3)

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS arguments
  input_file = arg(4)
  solver_file = arg(5)
  output_file = arg(6)
#else
  input_dir = arg(4)
  output_dir = arg(5)
#endif
  call synchronize_all()

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  if (nker > MAX_KERNEL_NAMES) stop 'number of kernel_names exceeds MAX_KERNEL_NAMES'

  if (myrank == 0) then
    ! The machinery for reading multiple names from the command line is in
    ! place,
    ! but the smoothing routines themselves have not yet been modified to work
    !  on multiple arrays.
    if (myrank == 0) print *, 'Smoothing list: ', trim(kernel_names_comma_delimited),' - total: ',nker
    if (myrank == 0) print *
  endif
  call synchronize_all()

  call read_parameter_file()
  topo_dir = trim(LOCAL_PATH)//'/'

  ! checks if basin code or global code: global code uses nchunks /= 0
  if (NCHUNKS == 0) stop 'Error nchunks'
  if (sizeprocs /= NPROCTOT_VAL) call exit_mpi(myrank,'Error total number of slices')

  ! estimates mesh element size
  ! note: this estimation is for global meshes valid only
  ! see values_from_mesher.h:
  !   average size of a spectral element in km = ...
  !   e.g. nproc 12x12, nex 192: element_size = 52.122262
  if (NCHUNKS_VAL == 6) then
    element_size = real(TWO_PI / dble(4) * R_EARTH_KM / dble(NEX_XI_VAL),kind=CUSTOM_REAL)
  else
    ANGULAR_WIDTH_XI_RAD = ANGULAR_WIDTH_XI_IN_DEGREES_VAL * DEGREES_TO_RADIANS
    ANGULAR_WIDTH_ETA_RAD = ANGULAR_WIDTH_ETA_IN_DEGREES_VAL * DEGREES_TO_RADIANS
    element_size = max( ANGULAR_WIDTH_XI_RAD/NEX_XI_VAL,ANGULAR_WIDTH_ETA_RAD/NEX_ETA_VAL ) * real(R_EARTH_KM,kind=CUSTOM_REAL)
  endif

  ! user output
  if (myrank == 0) then
    print *,"defaults:"
    print *,"  NPROC_XI , NPROC_ETA        : ",NPROC_XI,NPROC_ETA
    print *,"  NCHUNKS                     : ",NCHUNKS
    print *,"  element size on surface (km): ",element_size
    print *
    print *,"  smoothing sigma_h , sigma_v (km)                : ",sigma_h,sigma_v
    ! scalelength: approximately S ~ sigma * sqrt(8.0) for a Gaussian smoothing
    print *,"  smoothing scalelengths horizontal, vertical (km): ",sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *
    print *,"  data name      : ",trim(kernel_names_comma_delimited)
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS arguments
    print *,"  input file     : ",trim(input_file)
    print *,"  solver file    : ",trim(solver_file)
    print *,"  output file    : ",trim(output_file)
#else
    print *,"  input dir      : ",trim(input_dir)
    print *,"  output dir     : ",trim(output_dir)
#endif
    print *
    print *,"number of elements per slice: ",NSPEC_AB
    print *
  endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  call synchronize_all()
  ! initializes ADIOS
  if (myrank == 0) then
    print *, 'initializing ADIOS...'
    print *, ' '
  endif
  call initialize_adios()
#endif

  ! synchronizes
  call synchronize_all()

  ! initializes lengths (non-dimensionalizes)
  element_size_m = element_size / real(R_EARTH_KM,kind=CUSTOM_REAL) ! e.g. 9 km on the surface, 36 km at CMB

  sigma_h = sigma_h / real(R_EARTH_KM,kind=CUSTOM_REAL) ! scale
  sigma_v = sigma_v / real(R_EARTH_KM,kind=CUSTOM_REAL) ! scale

  sigma_h2 = 2.0 * sigma_h ** 2  ! factor two for Gaussian distribution with standard variance sigma
  sigma_v2 = 2.0 * sigma_v ** 2

  ! checks
  if (sigma_h2 < 1.e-18) stop 'Error sigma_h2 zero, must be non-zero'
  if (sigma_v2 < 1.e-18) stop 'Error sigma_v2 zero, must be non-zero'

  ! search radius
  ! note: crust/mantle region has at least 2 doubling layers, we enlarge the search radius according to element size at CMB;
  !       we always smooth among several elements, thus jumps at discontinuities will get averaged
  sigma_h3 = 3.0  * sigma_h + 4.0 * element_size_m
  sigma_v3 = 3.0  * sigma_v + 4.0 * element_size_m

  if (.not. DO_BRUTE_FORCE_SEARCH) then
    ! kd-tree search uses a slightly larger constant search radius
    r_search = max( 1.1d0 * sigma_h3, 1.1d0 * sigma_v3 )
    r_search_dist_v = sigma_v3
    r_search_dist_h = sigma_h3
  else
    r_search = 0.d0
  endif

  ! theoretic normal value
  ! (see integral over -inf to +inf of exp[- x*x/(2*sigma) ] = sigma * sqrt(2*pi) )

! note: smoothing is using a Gaussian (ellipsoid for sigma_h /= sigma_v),
!          but in spherical coordinates, we use horizontal distance as epicentral distance
!          and vertical distance as radial distance?

! not squared since epicentral distance is taken? values from bk seem to be closer to squared ones...
  norm_h = sqrt(2.0*PI) * sigma_h
  norm_h = norm_h * norm_h ! squared since 2 horizontal directions
  norm_v = sqrt(2.0*PI) * sigma_v
  norm   = norm_h * norm_v
  !norm = (sqrt(2.0*PI) * sigma) ** 3 ! for sigma_h = sigma_v = sigma

  ! GLL points weights
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        wgll_cube(i,j,k) = real(wxgll(i)*wygll(j)*wzgll(k),kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  ! ---- figure out the neighboring 8 or 7 slices: (ichunk,ixi,ieta) index start at 0------
  ichunk = myrank / (NPROC_XI * NPROC_ETA)
  ieta = (myrank - ichunk * NPROC_XI * NPROC_ETA) / NPROC_XI
  ixi = myrank - ichunk * NPROC_XI * NPROC_ETA - ieta * NPROC_XI

  ! get the neighboring slices:
  call get_all_eight_slices(ichunk,ixi,ieta, &
                            islice0(1),islice0(2),islice0(3),islice0(4),islice0(5),islice0(6),islice0(7),islice0(8), &
                            NPROC_XI,NPROC_ETA)

  ! remove the repeated slices (only 8 for corner slices in global case)
  islice(1) = myrank; j = 1
  do i = 1, 8
    if (.not. any(islice(1:i) == islice0(i)) .and. islice0(i) < sizeprocs) then
      j = j + 1
      islice(j) = islice0(i)
    endif
  enddo
  nums = j

  if (myrank == 0) then
    print *,'slices:',nums
    print *,'  rank:',myrank,'  smoothing slices'
    print *,'  ',islice(1:nums)
    print *
  endif

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           idoubling(NSPEC_AB), &
           ispec_is_tiso(NSPEC_AB), &
           xstore(NGLOB_AB), &
           ystore(NGLOB_AB), &
           zstore(NGLOB_AB), &
           xix(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiy(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammax(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammay(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammaz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating mesh arrays'

  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           kernel_smooth(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           tk(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           jacobian(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xx0(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           yy0(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           zz0(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xx(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           yy(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           zz(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           cx0(NSPEC_AB), &
           cy0(NSPEC_AB), &
           cz0(NSPEC_AB), &
           cx(NSPEC_AB), &
           cy(NSPEC_AB), &
           cz(NSPEC_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'

  ! initializes
  kernel(:,:,:,:,:) = 0.0_CUSTOM_REAL
  kernel_smooth(:,:,:,:,:) = 0.0_CUSTOM_REAL
  tk(:,:,:,:,:) = 0.0_CUSTOM_REAL
  bk(:,:,:,:) = 0.0_CUSTOM_REAL

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! user output
  if (myrank == 0) print *, 'reading in ADIOS solver file: ',trim(solver_file)

  ! opens file for reading
  call init_adios_group(myadios_group, "MeshReader")
  call open_file_adios_read_and_init_method(myadios_file, myadios_group, trim(solver_file))

  call read_adios_scalar(myadios_file, myadios_group, myrank, trim(reg_name)//"nspec", nspec)
  call read_adios_scalar(myadios_file, myadios_group, myrank, trim(reg_name)//"nglob", nglob)

  if (nspec /= NSPEC_AB) call exit_mpi(myrank,'Error invalid nspec value in solver_data.bp')
  if (nglob /= NGLOB_AB) call exit_mpi(myrank,'Error invalid nglob value in solver_data.bp')

  ! reads mesh arrays
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "ibool", ibool(:, :, :, :))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xixstore", xix(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xiystore", xiy(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xizstore", xiz(:, :, :, :))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammaxstore", gammax(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammaystore", gammay(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammazstore", gammaz(:, :, :, :))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etaxstore", etax(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etaystore", etay(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etazstore", etaz(:, :, :, :))

  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "x_global", xstore(:))
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "y_global", ystore(:))
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "z_global", zstore(:))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "idoubling", idoubling(:))
#else
  ! read in the topology files of the current and neighboring slices
  ! point locations
  write(filename,'(a,i6.6,a)') trim(topo_dir)//'/proc',myrank,trim(reg_name)//'solver_data.bin'
  open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    call exit_mpi(myrank,'Error opening solver_data.bin file')
  endif

  ! checks nspec
  read(IIN) nspec
  if (nspec /= NSPEC_AB) call exit_mpi(myrank,'Error invalid nspec value in solver_data.bin')

  ! checks nglob
  read(IIN) nglob
  if (nglob /= NGLOB_AB) call exit_mpi(myrank,'Error invalid nglob value in solver_data.bin')

  ! node locations
  read(IIN) xstore
  read(IIN) ystore
  read(IIN) zstore

  read(IIN) ibool
  read(IIN) idoubling
  read(IIN) ispec_is_tiso

  read(IIN) xix
  read(IIN) xiy
  read(IIN) xiz
  read(IIN) etax
  read(IIN) etay
  read(IIN) etaz
  read(IIN) gammax
  read(IIN) gammay
  read(IIN) gammaz
  close(IIN)
#endif

  ! synchronizes
  call synchronize_all()

  ! get the location of the center of the elements
  do ispec = 1, NSPEC_AB

    DO_LOOP_IJK

      iglob = ibool(INDEX_IJK,ispec)
      xx0(INDEX_IJK,ispec) = xstore(iglob)
      yy0(INDEX_IJK,ispec) = ystore(iglob)
      zz0(INDEX_IJK,ispec) = zstore(iglob)

      ! build jacobian
      ! get derivatives of ux, uy and uz with respect to x, y and z
      xixl = xix(INDEX_IJK,ispec)
      xiyl = xiy(INDEX_IJK,ispec)
      xizl = xiz(INDEX_IJK,ispec)
      etaxl = etax(INDEX_IJK,ispec)
      etayl = etay(INDEX_IJK,ispec)
      etazl = etaz(INDEX_IJK,ispec)
      gammaxl = gammax(INDEX_IJK,ispec)
      gammayl = gammay(INDEX_IJK,ispec)
      gammazl = gammaz(INDEX_IJK,ispec)
      ! compute the jacobian
      jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))

      jacobian(INDEX_IJK,ispec) = jacobianl

    ENDDO_LOOP_IJK

    cx0(ispec) = (xx0(1,1,1,ispec) + xx0(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cy0(ispec) = (yy0(1,1,1,ispec) + yy0(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cz0(ispec) = (zz0(1,1,1,ispec) + zz0(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
  enddo

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS single file opening
  ! user output
  if (myrank == 0) print *, 'reading in ADIOS input file : ',trim(input_file)
  ! open kernel file before looping over slices
  ! note: the number of neighbor slices might be different for different procs,
  !       thus no synchronization should occur within the loop over neighbor slices.
  !       adios however requires that all processes are synchronized for opening/closing files, so we do it here.
  call init_adios_group(myadios_val_group, "ValReader")
  call open_file_adios_read(myadios_val_file, myadios_val_group, trim(input_file))
#endif

  ! element search
  if (.not. DO_BRUTE_FORCE_SEARCH) then
    ! search by kd-tree
    ! user output
    if (myrank == 0) then
      print *
      print *,'using kd-tree search:'
      if (DO_SEARCH_ELLIP) then
        print *,'  search radius horizontal: ',r_search_dist_h * R_EARTH_KM,'km'
        print *,'  search radius vertical  : ',r_search_dist_v * R_EARTH_KM,'km'
      else
        print *,'  search sphere radius: ',r_search * R_EARTH_KM,'km'
      endif
      print *
    endif

    ! set number of tree nodes
    kdtree_num_nodes = NSPEC_AB

    ! checks
    if (kdtree_num_nodes <= 0) stop 'Error number of kd-tree nodes is zero, please check NSPEC_AB'

    ! allocates tree arrays
    allocate(kdtree_nodes_location(NDIM,kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
    allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
    if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

    ! tree verbosity
    if (myrank == 0) call kdtree_set_verbose()

  else
    ! brute-force search
    ! user output
    if (myrank == 0) then
      print *
      print *,'using brute-force search:'
      print *,'  search radius horizontal: ',sigma_h3 * R_EARTH_KM,'km'
      print *,'  search radius vertical  : ',sigma_v3 * R_EARTH_KM,'km'
      print *
    endif
  endif

  ! synchronizes
  call synchronize_all()

  if (myrank == 0) print *, 'start looping over elements and points for smoothing ...'

  ! synchronizes
  call synchronize_all()

  ! loop over all slices
  do inum = 1, nums

    ! timing
    time_start = wtime()

    iproc = islice(inum)

    ! user output
    if (myrank == 0) then
      print *
      print *,'slice ',inum,'out of ',nums
      print *,'  reading slice proc:',iproc
      print *
    endif

    ! debugging
    if (DEBUG .and. myrank == 0) then
      allocate(tmp_bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
      if (ier /= 0) stop 'Error allocating array tmp_bk'
      tmp_bk(:,:,:,:) = 0.0_CUSTOM_REAL
    endif

    ! sets up mesh locations
    if (iproc == myrank) then
      ! element centers
      cx(:) = cx0(:)
      cy(:) = cy0(:)
      cz(:) = cz0(:)
      ! locations
      xx(:,:,:,:) = xx0(:,:,:,:)
      yy(:,:,:,:) = yy0(:,:,:,:)
      zz(:,:,:,:) = zz0(:,:,:,:)
    else
      ! read in neighbor slice
#ifdef USE_ADIOS_INSTEAD_OF_MESH
      ! ADIOS
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "ibool", ibool(:, :, :, :))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "xixstore", xix(:, :, :, :))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "xiystore", xiy(:, :, :, :))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "xizstore", xiz(:, :, :, :))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "etaxstore", etax(:, :, :, :))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "etaystore", etay(:, :, :, :))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "etazstore", etaz(:, :, :, :))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "gammaxstore", gammax(:, :, :, :))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "gammaystore", gammay(:, :, :, :))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "gammazstore", gammaz(:, :, :, :))

      call read_adios_array(myadios_file, myadios_group, iproc, nglob, trim(reg_name) // "x_global", xstore(:))
      call read_adios_array(myadios_file, myadios_group, iproc, nglob, trim(reg_name) // "y_global", ystore(:))
      call read_adios_array(myadios_file, myadios_group, iproc, nglob, trim(reg_name) // "z_global", zstore(:))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "idoubling", idoubling(:))
#else
      ! neighbor database file
      write(filename,'(a,i6.6,a)') trim(topo_dir)//'/proc',iproc,trim(reg_name)//'solver_data.bin'

      ! read in the topology, kernel files, calculate center of elements
      ! point locations
      ! given in Cartesian coordinates
      open(IIN,file=trim(filename),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error could not open database file: ',trim(filename)
        call exit_mpi(myrank,'Error opening slices in solver_data.bin file')
      endif

      read(IIN) nspec
      read(IIN) nglob
      read(IIN) xstore
      read(IIN) ystore
      read(IIN) zstore

      read(IIN) ibool
      read(IIN) idoubling
      read(IIN) ispec_is_tiso

      read(IIN) xix
      read(IIN) xiy
      read(IIN) xiz
      read(IIN) etax
      read(IIN) etay
      read(IIN) etaz
      read(IIN) gammax
      read(IIN) gammay
      read(IIN) gammaz
      close(IIN)
#endif

      ! get the location of the center of the elements
      do ispec = 1, NSPEC_AB

        DO_LOOP_IJK

          iglob = ibool(INDEX_IJK,ispec)
          xx(INDEX_IJK,ispec) = xstore(iglob)
          yy(INDEX_IJK,ispec) = ystore(iglob)
          zz(INDEX_IJK,ispec) = zstore(iglob)

          ! build jacobian
          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xix(INDEX_IJK,ispec)
          xiyl = xiy(INDEX_IJK,ispec)
          xizl = xiz(INDEX_IJK,ispec)
          etaxl = etax(INDEX_IJK,ispec)
          etayl = etay(INDEX_IJK,ispec)
          etazl = etaz(INDEX_IJK,ispec)
          gammaxl = gammax(INDEX_IJK,ispec)
          gammayl = gammay(INDEX_IJK,ispec)
          gammazl = gammaz(INDEX_IJK,ispec)
          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))
          jacobian(INDEX_IJK,ispec) = jacobianl

        ENDDO_LOOP_IJK

        ! calculate element center location
        cx(ispec) = (xx(1,1,1,ispec) + xx(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
        cy(ispec) = (yy(1,1,1,ispec) + yy(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
        cz(ispec) = (zz(1,1,1,ispec) + zz(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
      enddo
    endif

    ! loops over input kernels
    do iker = 1,nker
      kernel_name = kernel_names(iker)
      ! user output
      if (myrank == 0) then
        print *,'  kernel ',iker,'out of ',nker
        print *,'  reading data file: proc = ',iproc,' name = ',trim(kernel_name)
        print *
      endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
      ! ADIOS array name
      ! determines if parameter name is for a kernel
      if (len_trim(kernel_name) > 3) then
        if (kernel_name(len_trim(kernel_name)-2:len_trim(kernel_name)) == '_kl') then
          is_kernel = .true.
        endif
      endif
      if (is_kernel) then
        ! NOTE: reg1 => crust_mantle, others are not implemented
        varname = trim(kernel_name) // "_crust_mantle"
      else
        varname = trim(reg_name) // trim(kernel_name)
      endif
      ! user output
      if (myrank == 0) then
        print *, '  data: ADIOS ',trim(kernel_name), " is_kernel = ", is_kernel
        print *, '  data: ADIOS using array name = ',trim(varname)
        print *
      endif
      ! reads kernel values
      call read_adios_array(myadios_val_file, myadios_val_group, iproc, nspec, trim(varname), kernel(:, :, :, :, iker))
#else
      ! data file
      write(local_data_file,'(a,i6.6,a)') &
        trim(input_dir)//'/proc',iproc,trim(reg_name)//trim(kernel_name)//'.bin'

      open(IIN,file=trim(local_data_file),status='old',action='read',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print *,'Error opening data file: ',trim(local_data_file)
        call exit_mpi(myrank,'Error opening data file')
      endif

      read(IIN) kernel(:,:,:,:,iker)
      close(IIN)
#endif

      ! statistics
      ! get the global maximum value of the original kernel file
      if (iproc == myrank) then
        min_old(iker) = minval(kernel(:,:,:,:,iker))
        max_old(iker) = maxval(kernel(:,:,:,:,iker))
      endif
    enddo
    if (myrank == 0) print *

    ! search setup
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      ! search by kd-tree

      ! prepares search arrays, each element takes its midpoint for tree search
      kdtree_nodes_index(:) = 0
      kdtree_nodes_location(:,:) = 0.0

      ! fills kd-tree arrays
      do ispec = 1,NSPEC_AB
        ! adds node index ( index points to same ispec for all internal GLL points)
        kdtree_nodes_index(ispec) = ispec

        ! adds node location
        kdtree_nodes_location(1,ispec) = cx(ispec)
        kdtree_nodes_location(2,ispec) = cy(ispec)
        kdtree_nodes_location(3,ispec) = cz(ispec)
      enddo

      ! creates kd-tree for searching
      call kdtree_setup()
    endif

    ! loops over elements to be smoothed in the current slice
    do ispec = 1, NSPEC_AB

      ! user info about progress
      if (myrank == 0) then
        tCPU = wtime() - time_start
        if (mod(ispec-1,max(int(1.0*NSPEC_AB/NSTEP_PERCENT_INFO),1)) == 0 &
            .and. ispec < (NSPEC_AB - 0.5*NSPEC_AB/NSTEP_PERCENT_INFO)) then
          print *,'    ',int((ispec-1) / max(int(1.0*NSPEC_AB/NSTEP_PERCENT_INFO),1)) * (100.0 / NSTEP_PERCENT_INFO), &
                 ' % elements done - Elapsed time in seconds = ',tCPU
        endif
        if (ispec == NSPEC_AB) then
          print *,'    ',100.0,' % elements done - Elapsed time in seconds = ',tCPU
        endif
      endif

      ! initializes
      num_elem_local = 0

      ! element center position
      center_x0 = cx0(ispec)
      center_y0 = cy0(ispec)
      center_z0 = cz0(ispec)

      ! sets number of elements to loop over
      if (.not. DO_BRUTE_FORCE_SEARCH) then
        xyz_target(1) = center_x0
        xyz_target(2) = center_y0
        xyz_target(3) = center_z0

        ! counts nearest neighbors
        if (DO_SEARCH_ELLIP) then
          ! (within ellipsoid)
          call kdtree_count_nearest_n_neighbors_ellip(xyz_target,r_search_dist_v,r_search_dist_h,num_elem_local)
        else
          ! (within search sphere)
          call kdtree_count_nearest_n_neighbors(xyz_target,r_search,num_elem_local)
        endif

        ! checks that at least a single element was choosen
        if (iproc == myrank) then
          if (num_elem_local < 1) stop 'Error no local search element found'
        endif

        ! sets n-search number of nodes
        kdtree_search_num_nodes = num_elem_local

        ! allocates search array
        if (kdtree_search_num_nodes > 0) then
          allocate(kdtree_search_index(kdtree_search_num_nodes),stat=ier)
          if (ier /= 0) stop 'Error allocating array kdtree_search_index'

          ! finds closest point in target chunk
          if (DO_SEARCH_ELLIP) then
            ! (within ellipsoid)
            call kdtree_get_nearest_n_neighbors_ellip(xyz_target,r_search_dist_v,r_search_dist_h,num_elem_local)
          else
            ! (within search sphere)
            call kdtree_get_nearest_n_neighbors(xyz_target,r_search,num_elem_local)
          endif
        endif

        ! debug output
        if (DEBUG .and. myrank == 0) then
          ! user info
          if (ispec < 10) then
            print *,'  total number of search elements: ',num_elem_local,'ispec',ispec
          endif
          ! file output
          if (ispec == tmp_ispec_dbg) then
            ! outputs search elements with integer flags
            allocate(ispec_flag(NSPEC_AB))
            ispec_flag(:) = 0
            ! debug element (tmp_ispec_dbg)
            do ielem = 1,num_elem_local
              i = kdtree_search_index(ielem)
              ispec_flag(i) = ielem
            enddo
            ! writes out vtk file for debugging
            write(filename,'(a,i4.4,a,i6.6)') trim(LOCAL_PATH)//'/search_elem',tmp_ispec_dbg,'_proc',iproc
            call write_VTK_data_elem_i(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore, &
                                       ibool,ispec_flag,filename)
            print *,'file written: ',trim(filename)//'.vtk'
            deallocate(ispec_flag)
          endif
        endif
      else
        ! brute-force search always loops over whole mesh slice
        num_elem_local = NSPEC_AB
      endif

      ! --- only double loop over the elements in the search radius ---
      do ielem = 1, num_elem_local
        ! search element
        if (.not. DO_BRUTE_FORCE_SEARCH) then
          ! kd-tree search elements
          ispec2 = kdtree_search_index(ielem)
        else
          ! brute-force search
          ispec2 = ielem
        endif

        ! search element center position
        center_x = cx(ispec2)
        center_y = cy(ispec2)
        center_z = cz(ispec2)

        ! takes only elements in search radius
        ! calculates horizontal and vertical distance between two element centers
        call get_distance_vec(dist_h,dist_v,center_x0,center_y0,center_z0,center_x,center_y,center_z)

        ! checks distance between centers of elements
        ! note: distances and sigmah, sigmav are normalized by R_EARTH_KM
        if ( dist_h > sigma_h3 .or. dist_v > sigma_v3 ) cycle

        ! integration factors:
        ! uses volume assigned to GLL points
        factor(:,:,:) = jacobian(:,:,:,ispec2) * wgll_cube(:,:,:)
        ! no volume
        !factor(:,:,:) = 1.0_CUSTOM_REAL

        ! loops over all GLL points of reference element in current slice (ispec) and search elements (ispec2)
        DO_LOOP_IJK

          ! reference location
          ! current point (i,j,k,ispec) location, Cartesian coordinates
          x0 = xx0(INDEX_IJK,ispec)
          y0 = yy0(INDEX_IJK,ispec)
          z0 = zz0(INDEX_IJK,ispec)

          ! calculate weights based on Gaussian smoothing
          call smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val, &
                                     xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

          ! debug
          !if (DEBUG .and. myrank == 0) then
          !  if (ispec == tmp_ispec_dbg) then
#ifdef FORCE_VECTORIZATION
          !    if (ijk == (2 + 3*NGLLX + 4*NGLLY*NGLLX)) then
#else
          !    if (i == 2 .and. j == 3 .and. k == 4) then
#endif
          !      tmp_bk(:,:,:,ispec2) = exp_val(:,:,:)
          !      print *,'debug',myrank,'ispec',ispec,'ispec2',ispec2,'dist:',dist_h,dist_v
          !      print *,'debug exp',minval(exp_val),maxval(exp_val)
          !      print *,'debug kernel',minval(kernel(:,:,:,ispec2,1)),maxval(kernel(:,:,:,ispec2,1))
          !    endif
          !  endif
          !endif

          ! adds GLL integration weights
          exp_val(:,:,:) = exp_val(:,:,:) * factor(:,:,:)

          ! normalization, integrated values of Gaussian smoothing function
          bk(INDEX_IJK,ispec) = bk(INDEX_IJK,ispec) + sum(exp_val(:,:,:))

          ! adds contribution of element ispec2 to smoothed kernel values
          do iker = 1,nker
            tk(INDEX_IJK,ispec,iker) = tk(INDEX_IJK,ispec,iker) + sum(exp_val(:,:,:) * kernel(:,:,:,ispec2,iker))
          enddo

          ! checks number
          !if (isNaN(tk(INDEX_IJK,ispec,1))) then
          !  print *,'Error tk NaN: ',tk(INDEX_IJK,ispec,1)
          !  print *,'rank:',myrank
          !  print *,'INDEX_IJK,ispec:',INDEX_IJK,ispec
          !  print *,'tk: ',tk(INDEX_IJK,ispec,1),'bk:',bk(INDEX_IJK,ispec)
          !  print *,'sum exp_val: ',sum(exp_val(:,:,:)),'sum factor:',sum(factor(:,:,:))
          !  print *,'sum kernel:',sum(kernel(:,:,:,ispec2))
          !  call exit_mpi(myrank, 'Error NaN')
          !endif

        ENDDO_LOOP_IJK

      enddo ! (ielem)

      ! frees search results
      if (.not. DO_BRUTE_FORCE_SEARCH) then
        if (kdtree_search_num_nodes > 0) then
          deallocate(kdtree_search_index)
          kdtree_search_num_nodes = 0
        endif
      endif

      ! debug output
      if (DEBUG .and. myrank == 0) then
        ! debug element (tmp_ispec_dbg)
        if (ispec == tmp_ispec_dbg) then
          ! outputs Gaussian weighting function
          ! writes out vtk file for debugging
          write(filename,'(a,i4.4,a,i6.6)') trim(LOCAL_PATH)//'/search_elem',tmp_ispec_dbg,'_Gaussian_proc',iproc
          call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,tmp_bk,filename)
          print *,'file written: ',trim(filename)//'.vtk'
        endif
      endif

    enddo   ! (ispec)

    ! debug output
    if (DEBUG .and. myrank == 0) then
      deallocate(tmp_bk)
    endif

    ! deletes search tree nodes
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      call kdtree_delete()
    endif
    if (myrank == 0) print *

  enddo     ! islice

  if (myrank == 0) print *

  ! frees memory
  if (.not. DO_BRUTE_FORCE_SEARCH) then
    deallocate(kdtree_nodes_location)
    deallocate(kdtree_nodes_index)
  endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! closes kernel file
  call close_file_adios_read(myadios_val_file)
#endif

  ! synchronizes
  call synchronize_all()

  ! normalizes/scaling factor
  if (myrank == 0) then
    print *, 'Scaling values:'
    print *, '  bk min/max = ',minval(bk),maxval(bk)
    print *, '  theoretical norm = ',norm
    do iker = 1,nker
      print *, '  ',trim(kernel_names(iker)),' tk min/max = ',minval(tk(:,:,:,:,iker)),maxval(tk(:,:,:,:,iker))
    enddo
    print *
  endif

  ! compute the smoothed kernel values
  do ispec = 1, NSPEC_AB

    ! avoids division by zero
    DO_LOOP_IJK
      if (bk(INDEX_IJK,ispec) == 0.0_CUSTOM_REAL) bk(INDEX_IJK,ispec) = 1.0_CUSTOM_REAL
    ENDDO_LOOP_IJK

    ! loops over kernels
    do iker = 1,nker

      DO_LOOP_IJK
        ! checks the normalization criterion
        ! e.g. sigma_h 160km, sigma_v 40km:
        !     norm (not squared sigma_h ) ~ 0.001
        !     norm ( squared sigma_h) ~ 6.23 * e-5
        !if (abs(bk(INDEX_IJK,ispec) - norm) > 1.e-4) then
        !  print *, 'Problem norm here --- ', myrank, ispec, i, j, k, bk(INDEX_IJK,ispec), norm
        !  !call exit_mpi(myrank, 'Error computing Gaussian function on the grid')
        !endif

        !debug
        !if (tk(INDEX_IJK,ispec) == 0.0_CUSTOM_REAL .and. myrank == 0) then
        !  print *,myrank,'zero tk: ',INDEX_IJK,ispec,'tk:',tk(INDEX_IJK,ispec),'bk:',bk(INDEX_IJK,ispec)
        !endif

        ! normalizes smoothed kernel values by integral value of Gaussian weighting
        kernel_smooth(INDEX_IJK,ispec,iker) = tk(INDEX_IJK,ispec,iker) / bk(INDEX_IJK,ispec)

        ! checks number (isNaN check)
        if (kernel_smooth(INDEX_IJK,ispec,iker) /= kernel_smooth(INDEX_IJK,ispec,iker)) then
          print *,'Error kernel_smooth value not a number: ',kernel_smooth(INDEX_IJK,ispec,iker),trim(kernel_names(iker))
          print *,'rank:',myrank
          print *,'INDEX_IJK,ispec,iker:',INDEX_IJK,ispec,iker
          print *,'tk: ',tk(INDEX_IJK,ispec,iker),'bk:',bk(INDEX_IJK,ispec),'norm:',norm
          call exit_mpi(myrank, 'Error kernel value is NaN')
        endif

      ENDDO_LOOP_IJK

    enddo
  enddo

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! user output
  if (myrank == 0) print *, 'writing to ADIOS output file: ',trim(output_file)
  ! determines group size
  call init_adios_group(myadios_val_group, "ValWriter")
  group_size_inc = 0
  call define_adios_scalar(myadios_val_group, group_size_inc, '', trim(reg_name)//"nspec", nspec)
  call define_adios_scalar(myadios_val_group, group_size_inc, '', trim(reg_name)//"nglob", nglob)
  do iker = 1,nker
    ! name
    kernel_name = kernel_names(iker)
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, trim(reg_name), &
                                     trim(kernel_name), kernel_smooth(:, :, :, :, iker))
  enddo
  ! opens output files
  call open_file_adios_write(myadios_val_file, myadios_val_group, trim(output_file), "ValWriter")
  call set_adios_group_size(myadios_val_file,group_size_inc)
  ! writes to file
  call write_adios_scalar(myadios_val_file, myadios_val_group, trim(reg_name)//"nspec", nspec)
  call write_adios_scalar(myadios_val_file, myadios_val_group, trim(reg_name)//"nglob", nglob)
#endif

  ! output
  do iker = 1,nker
    ! name
    kernel_name = kernel_names(iker)
    if (myrank == 0) then
      print *,'smoothed: ',trim(kernel_name)
    endif

    ! file output
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    ! smoothed kernel values
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(reg_name) // trim(kernel_name), kernel_smooth(:, :, :, :, iker))
#else
    ! smoothed kernel file name
    write(output_file,'(a,i6.6,a)') trim(output_dir)//'/proc', myrank, trim(reg_name)//trim(kernel_name)//'_smooth.bin'

    open(IOUT,file=trim(output_file),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'Error opening smoothed kernel file')

    ! Note: output the following instead of kernel_smooth(:,:,:,1:NSPEC_AB) to create files of the same sizes
    write(IOUT) kernel_smooth(:,:,:,:,iker)
    close(IOUT)
    if (myrank == 0) print *,'  written: ',trim(output_file)
#endif

    ! statistics
    min_new = minval(kernel_smooth(:,:,:,:,iker))
    max_new = maxval(kernel_smooth(:,:,:,:,iker))

    ! the minimum/maximum value for the smoothed kernel
    call min_all_cr(min_old(iker), min_old_all)
    call min_all_cr(min_new, min_new_all)
    call max_all_cr(max_old(iker), max_old_all)
    call max_all_cr(max_new, max_new_all)

    if (myrank == 0) then
      print *
      print *, '  Minimum data value before smoothing = ', min_old_all
      print *, '  Minimum data value after smoothing  = ', min_new_all
      print *
      print *, '  Maximum data value before smoothing = ', max_old_all
      print *, '  Maximum data value after smoothing  = ', max_new_all
      print *
    endif
    ! synchronizes
    call synchronize_all()
  enddo

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! actual writing
  ! (note: done at the very end, otherwise it will reset path and we would need to re-initiate a group)
  call write_adios_perform(myadios_val_file)
  ! closes adios files
  call close_file_adios(myadios_val_file)
  call close_file_adios_read_and_finalize_method(myadios_file)
  ! finalizes adios
  call finalize_adios()
#endif

  ! user output
  if (myrank == 0) print *, 'all done'

  ! stop all the processes, and exit
  call finalize_mpi()

end program smooth_sem_globe
