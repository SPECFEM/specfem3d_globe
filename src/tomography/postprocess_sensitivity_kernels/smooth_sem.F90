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

! XSMOOTH_SEM
!
! USAGE
!   mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR USE_GPU
!
!
! COMMAND LINE ARGUMENTS
!   SIGMA_H                - horizontal smoothing radius
!   SIGMA_V                - vertical smoothing radius
!   KERNEL_NAME            - kernel name, e.g. alpha_kernel
!   INPUT_DIR              - directory from which arrays are read
!   OUTPUT_DIR             - directory to which smoothed array are written
!   USE_GPU                - (optional) use GPUs for computation
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

  use constants, only: myrank

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLCUBE,NDIM,IIN,IOUT, &
    GAUSSALPHA,GAUSSBETA,PI,TWO_PI,MAX_STRING_LEN,DEGREES_TO_RADIANS, &
    USE_QUADRATURE_RULE_FOR_SMOOTHING,USE_VECTOR_DISTANCE_FOR_SMOOTHING

  use shared_parameters, only: R_PLANET_KM

  ! GPU
  use shared_parameters, only: GPU_RUNTIME,GPU_PLATFORM,GPU_DEVICE

  use postprocess_par, only: &
    NCHUNKS_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,NPROCTOT_VAL,NEX_XI_VAL,NEX_ETA_VAL, &
    ANGULAR_WIDTH_XI_IN_DEGREES_VAL,ANGULAR_WIDTH_ETA_IN_DEGREES_VAL, &
    NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,MAX_KERNEL_NAMES,LOCAL_PATH

  use kdtree_search

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use adios_helpers_mod
  use manager_adios
#endif

  use iso_c_binding, only: C_NULL_CHAR

  implicit none

  !-------------------------------------------------------------
  ! Parameters
  integer,parameter :: NGLL3 = NGLLX * NGLLY * NGLLZ

  ! only include the neighboring 3 x 3 slices
  integer, parameter :: NSLICES = 3
  integer ,parameter :: NSLICES2 = NSLICES * NSLICES

  ! arguments
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer, parameter :: NARGS = 7
  character(len=*), parameter :: reg_name = 'reg1/'
#else
  integer, parameter :: NARGS = 6
  character(len=*), parameter :: reg_name = '_reg1_'
#endif
  !-------------------------------------------------------------

  integer :: islice(NSLICES2), islice0(NSLICES2)

  integer, dimension(:,:,:,:), allocatable :: ibool
  integer, dimension(:), allocatable :: idoubling
  logical, dimension(:), allocatable :: ispec_is_tiso

  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore, ystore, zstore
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: kernel, kernel_smooth,tk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: bk
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: xx0, yy0, zz0, xx, yy, zz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: integ_factor
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: cx0, cy0, cz0, cx, cy, cz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: jacobian

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) :: sigma_h, sigma_h2, sigma_h3, sigma_v, sigma_v2, sigma_v3
  real(kind=CUSTOM_REAL) :: norm, norm_h, norm_v, element_size_m
  real(kind=CUSTOM_REAL) :: center_x0, center_y0, center_z0

  real(kind=CUSTOM_REAL),dimension(:),allocatable :: min_old,max_old
  real(kind=CUSTOM_REAL) :: max_new, min_new
  real(kind=CUSTOM_REAL) :: max_old_all, max_new_all, min_old_all, min_new_all

  ! local copies of mesh parameters
  integer :: NPROC_XI
  integer :: NPROC_ETA
  integer :: NCHUNKS
  integer :: NSPEC_AB
  integer :: NGLOB_AB

  integer :: nspec, nglob
  integer :: ier,ichunk, ixi, ieta, iglob, num_slices
  integer :: ispec,iproc,inum
  integer :: i,j,k
  integer :: sizeprocs

  ! GPU
  integer :: ngpu_devices
  ! smooth structure pointer
  integer(kind=8) :: Container

  ! arguments
  character(len=MAX_STRING_LEN),dimension(NARGS) :: arg
  character(len=MAX_STRING_LEN),dimension(:),allocatable :: kernel_names
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: kernel_name, topo_dir
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  character(len=MAX_STRING_LEN) :: kernel_output_name
  character(len=MAX_STRING_LEN),dimension(:),allocatable :: kernel_output_names
  character(len=MAX_STRING_LEN) :: input_file, solver_file
  character(len=MAX_STRING_LEN) :: varname
#else
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  character(len=MAX_STRING_LEN) :: local_data_file
#endif
  character(len=MAX_STRING_LEN) :: output_file
  character(len=MAX_STRING_LEN) :: filename
  logical :: USE_GPU

  integer :: nker,iker

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll, wxgll
  double precision, dimension(NGLLY) :: yigll, wygll
  double precision, dimension(NGLLZ) :: zigll, wzgll

  ! array with all the weights in the cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube

  real(kind=CUSTOM_REAL) :: element_size
  real(kind=CUSTOM_REAL) :: ANGULAR_WIDTH_XI_RAD,ANGULAR_WIDTH_ETA_RAD

  ! search elements
  integer :: ielem,num_elem_local,num_elem_local_max
  integer, dimension(:), allocatable :: nsearch_elements

  ! tree nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: r_search,r_search_dist_v,r_search_dist_h
  logical :: use_kdtree_search

  ! debugging
  integer, dimension(:),allocatable :: ispec_flag
  !real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tmp_bk
  double precision :: memory_size

  ! timing
  integer :: ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start_all
  double precision :: tCPU
  double precision, external :: wtime

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer(kind=8) :: group_size_inc,local_dim
  logical :: is_kernel = .false.
#endif

#ifdef FORCE_VECTORIZATION
! in this vectorized version we have to assume that N_SLS == 3 in order to be able to unroll and thus suppress
! an inner loop that would otherwise prevent vectorization; this is safe in practice in all cases because N_SLS == 3
! in all known applications, and in the main program we check that N_SLS == 3 if FORCE_VECTORIZATION is used and we stop
  integer :: ijk
#endif

  !------------------------------------------------------------------------------------------
  ! USER PARAMETERS

  ! number of steps to reach 100 percent, i.e. 20 outputs info for every 5 percent
  integer,parameter :: NSTEP_PERCENT_INFO = 20     ! or 10,..

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
  if (command_argument_count() < NARGS-1) then
    if (myrank == 0) then
#ifdef USE_ADIOS_INSTEAD_OF_MESH
      print *,'Usage: mpirun -np NPROC bin/xsmooth_sem_adios SIGMA_H SIGMA_V KERNEL_NAME INPUT_FILE SOLVER_FILE OUTPUT_FILE'
      print *,'   with'
      print *,'     SIGMA_H, SIGMA_V - horizontal and vertical smoothing lenghts'
      print *,'     KERNEL_NAME      - comma-separated kernel names (e.g., alpha_kernel,beta_kernel)'
      print *,'     INPUT_FILE       - ADIOS file with kernel values (e.g., kernels.bp)'
      print *,'     SOLVER_FILE      - ADIOS file with mesh arrays (e.g., DATABASES_MPI/solver_data.bp)'
      print *,'     OUTPUT_FILE      - ADIOS file for smoothed output'
#else
      print *, 'Usage: mpirun -np NPROC bin/xsmooth_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR'
      print *,'   with'
      print *,'     SIGMA_H, SIGMA_V - horizontal and vertical smoothing lenghts'
      print *,'     KERNEL_NAME      - comma-separated kernel names (e.g., alpha_kernel,beta_kernel)'
      print *,'     INPUT_DIR        - directory with kernel files (e.g., proc***_alpha_kernel.bin)'
      print *,'     OUTPUT_DIR       - directory for smoothed output files'
#endif
      print *,'     GPU_MODE         - (optional) set to .true. to use GPU, otherwise set to .false. for CPU run (default off)'
      print *
      stop ' Please check command line arguments'
    endif
  endif
  call synchronize_all()

  if (myrank == 0) then
    print *, 'Running XSMOOTH_SEM'
    print *
  endif
  call synchronize_all()

  ! timing
  time_start_all = wtime()

  ! allocates array
  allocate(kernel_names(MAX_KERNEL_NAMES), &
           min_old(MAX_KERNEL_NAMES), &
           max_old(MAX_KERNEL_NAMES),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel_names array'
  kernel_names(:) = ''
  min_old(:) = 0._CUSTOM_REAL
  max_old(:) = 0._CUSTOM_REAL
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  allocate(kernel_output_names(MAX_KERNEL_NAMES), stat=ier)
  if (ier /= 0) stop 'Error allocating kernel_output_names array'
  kernel_output_names(:) = ''
#endif

  ! parse command line arguments
  do i = 1, NARGS
    if (command_argument_count() >= i) then
      call get_command_argument(i,arg(i), status=ier)
    else
      ! optional argument
      arg(i) = ''
    endif
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

  if (command_argument_count() == NARGS) then
    read(arg(NARGS),*) USE_GPU
  else
    USE_GPU = .false.
  endif

  ! sets flag to use kdtree
  use_kdtree_search = (.not. DO_BRUTE_FORCE_SEARCH) .and. (.not. USE_GPU)

  call synchronize_all()

  call parse_kernel_names(kernel_names_comma_delimited,kernel_names,nker)
  if (nker > MAX_KERNEL_NAMES) stop 'number of kernel_names exceeds MAX_KERNEL_NAMES'

  if (myrank == 0) then
    print *, 'Smoothing list: ', trim(kernel_names_comma_delimited),' - total: ',nker
    print *
  endif
  call synchronize_all()

  ! reads in Par_file and sets compute parameters
  call read_compute_parameters()

  ! reads mesh parameters
  if (myrank == 0) then
    ! reads mesh_parameters.bin file from LOCAL_PATH
    call read_mesh_parameters()
  endif
  ! broadcast parameters to all processes
  call bcast_mesh_parameters()

  topo_dir = trim(LOCAL_PATH)//'/'

  ! user output
  if (myrank == 0) then
    print *,'mesh parameters (from local_path parameter):'
    print *,'  LOCAL_PATH         = ',trim(LOCAL_PATH)
    print *,'  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
    print *,'  NPROCTOT           = ',NPROCTOT_VAL
    print *
  endif

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

  ! copy from static compilation (depends on Par_file values)
  NPROC_XI  = NPROC_XI_VAL
  NPROC_ETA = NPROC_ETA_VAL
  NCHUNKS   = NCHUNKS_VAL

  !takes region 1 kernels
  NSPEC_AB = NSPEC_CRUST_MANTLE
  NGLOB_AB = NGLOB_CRUST_MANTLE

  ! checks if basin code or global code: global code uses nchunks /= 0
  if (NCHUNKS == 0) stop 'Error NCHUNKS'

  ! estimates mesh element size
  ! note: this estimation is for global meshes valid only
  ! see values_from_mesher.h:
  !   average size of a spectral element in km = ...
  !   e.g. nproc 12x12, nex 192: element_size = 52.122262
  if (NCHUNKS_VAL == 6) then
    element_size = real(TWO_PI / dble(4) * R_PLANET_KM / dble(NEX_XI_VAL),kind=CUSTOM_REAL)
  else
    ANGULAR_WIDTH_XI_RAD = real(ANGULAR_WIDTH_XI_IN_DEGREES_VAL * DEGREES_TO_RADIANS,kind=CUSTOM_REAL)
    ANGULAR_WIDTH_ETA_RAD = real(ANGULAR_WIDTH_ETA_IN_DEGREES_VAL * DEGREES_TO_RADIANS,kind=CUSTOM_REAL)
    element_size = max( ANGULAR_WIDTH_XI_RAD/NEX_XI_VAL,ANGULAR_WIDTH_ETA_RAD/NEX_ETA_VAL ) * real(R_PLANET_KM,kind=CUSTOM_REAL)
  endif

  ! user output
  if (myrank == 0) then
    print *,'defaults:'
    print *,'  NPROC_XI / NPROC_ETA         : ',NPROC_XI,NPROC_ETA
    print *,'  NEX_XI   / NEX_ETA           : ',NEX_XI_VAL,NEX_ETA_VAL
    print *,'  NGLLX    / NGLLY     / NGLLZ : ',NGLLX,NGLLY,NGLLZ
    print *,'  NCHUNKS                      : ',NCHUNKS
    print *,'  element size on surface (km) : ',element_size
    print *
    print *,'  smoothing sigma_h , sigma_v (km)                : ',sigma_h,sigma_v
    ! scalelength: approximately S ~ sigma * sqrt(8.0) ~ sigma * 2.8 for a Gaussian smoothing
    print *,'  smoothing scalelengths horizontal, vertical (km): ',sigma_h*sqrt(8.0),sigma_v*sqrt(8.0)
    print *
    print *,'  data name      : ',trim(kernel_names_comma_delimited)
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS arguments
    print *,'  input file     : ',trim(input_file)
    print *,'  solver file    : ',trim(solver_file)
    print *,'  output file    : ',trim(output_file)
#else
    print *,'  input dir      : ',trim(input_dir)
    print *,'  output dir     : ',trim(output_dir)
#endif
    print *
    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      print *,'  using quadrature rule for smoothing integration'
    else
      print *,'  using distance weighted smoothing'
    endif
    if (USE_VECTOR_DISTANCE_FOR_SMOOTHING) then
      print *,'  using vector approximation for distance calculation'
    else
      print *,'  using epicentral distance for distance calculation'
    endif
    if (USE_GPU) then
      print *,'  using GPU'
    endif
    print *
    print *,'number of elements per slice: ',NSPEC_AB
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

  ! sets GPU devices
  if (USE_GPU) then
    call initialize_gpu_device(GPU_RUNTIME,trim(GPU_PLATFORM)//C_NULL_CHAR,trim(GPU_DEVICE)//C_NULL_CHAR,myrank,ngpu_devices)
    ! checks if found
    if (ngpu_devices == 0) stop 'no GPU devices found, please switch to CPU smoothing'
  endif

  ! synchronizes
  call synchronize_all()

  ! initializes lengths (non-dimensionalizes)
  element_size_m = element_size / real(R_PLANET_KM,kind=CUSTOM_REAL) ! e.g. 9 km on the surface, 36 km at CMB

  sigma_h = sigma_h / real(R_PLANET_KM,kind=CUSTOM_REAL) ! scale
  sigma_v = sigma_v / real(R_PLANET_KM,kind=CUSTOM_REAL) ! scale

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

  if (use_kdtree_search) then
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
  norm_h = real(sqrt(2.d0*PI) * sigma_h,kind=CUSTOM_REAL)
  norm_h = norm_h * norm_h ! squared since 2 horizontal directions
  norm_v = real(sqrt(2.d0*PI) * sigma_v,kind=CUSTOM_REAL)
  norm   = norm_h * norm_v
  !norm = (sqrt(2.0*PI) * sigma) ** 3 ! for sigma_h = sigma_v = sigma

  ! GLL points weights
  if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
    call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
    call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
    call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          wgll_cube(i,j,k) = real(wxgll(i)*wygll(j)*wzgll(k),kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  endif

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
  num_slices = j

  ! user output
  if (myrank == 0) then
    print *,'slices:',num_slices
    print *,'  rank:',myrank,'  smoothing slices'
    print *,'  ',islice(1:num_slices)
    print *
  endif

  ! debug
  !do iproc = 0,NPROCTOT_VAL-1
  !  if (iproc == myrank) then
  !    print *,'debug: rank:',myrank,' has smoothing slices',islice(1:num_slices)
  !    print *
  !  endif
  !enddo

  ! note: for the globe version, all slices have the same number of elements NSPEC_AB
  !       we can thus allocate all mesh arrays here already, valid for all slices to loop over

  ! mesh
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           idoubling(NSPEC_AB), &
           ispec_is_tiso(NSPEC_AB), &
           xstore(NGLOB_AB), &
           ystore(NGLOB_AB), &
           zstore(NGLOB_AB), &
           xixstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xiystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           xizstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etaystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           etazstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammaxstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammaystore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           gammazstore(NGLLX,NGLLY,NGLLZ,NSPEC_AB), stat=ier)
  if (ier /= 0) stop 'Error allocating mesh arrays'

  ibool(:,:,:,:) = 0; idoubling(:) = 0; ispec_is_tiso(:) = .false.
  xstore(:) = 0.0_CUSTOM_REAL; ystore(:) = 0.0_CUSTOM_REAL; zstore(:) = 0.0_CUSTOM_REAL
  xixstore(:,:,:,:) = 0.0_CUSTOM_REAL; xiystore(:,:,:,:) = 0.0_CUSTOM_REAL; xizstore(:,:,:,:) = 0.0_CUSTOM_REAL
  etaxstore(:,:,:,:) = 0.0_CUSTOM_REAL; etaystore(:,:,:,:) = 0.0_CUSTOM_REAL; etazstore(:,:,:,:) = 0.0_CUSTOM_REAL
  gammaxstore(:,:,:,:) = 0.0_CUSTOM_REAL; gammaystore(:,:,:,:) = 0.0_CUSTOM_REAL; gammazstore(:,:,:,:) = 0.0_CUSTOM_REAL

  ! weights and positions
  allocate(bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
           integ_factor(NGLLX,NGLLY,NGLLZ,NSPEC_AB), &
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
  if (ier /= 0) stop 'Error allocating position arrays'

  ! initializes
  bk(:,:,:,:) = 0.0_CUSTOM_REAL
  integ_factor(:,:,:,:) = 0.0_CUSTOM_REAL
  xx0(:,:,:,:) = 0.0_CUSTOM_REAL; yy0(:,:,:,:) = 0.0_CUSTOM_REAL; zz0(:,:,:,:) = 0.0_CUSTOM_REAL
  xx(:,:,:,:) = 0.0_CUSTOM_REAL; yy(:,:,:,:) = 0.0_CUSTOM_REAL; zz(:,:,:,:) = 0.0_CUSTOM_REAL
  cx0(:) = 0.0_CUSTOM_REAL; cy0(:) = 0.0_CUSTOM_REAL; cz0(:) = 0.0_CUSTOM_REAL
  cx(:) = 0.0_CUSTOM_REAL; cy(:) = 0.0_CUSTOM_REAL; cz(:) = 0.0_CUSTOM_REAL

  ! Gaussian kernels
  allocate(kernel(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           kernel_smooth(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), &
           tk(NGLLX,NGLLY,NGLLZ,NSPEC_AB,nker), stat=ier)
  if (ier /= 0) stop 'Error allocating kernel arrays'
  kernel(:,:,:,:,:) = 0.0_CUSTOM_REAL
  kernel_smooth(:,:,:,:,:) = 0.0_CUSTOM_REAL
  tk(:,:,:,:,:) = 0.0_CUSTOM_REAL

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! user output
  if (myrank == 0) print *, 'reading in ADIOS solver file: ',trim(solver_file)

  ! opens file for reading
  call init_adios_group(myadios_group, "MeshReader")
  call open_file_adios_read_and_init_method(myadios_file, myadios_group, trim(solver_file))

  call read_adios_scalar(myadios_file, myadios_group, myrank, trim(reg_name)//"nspec", nspec)
  call read_adios_scalar(myadios_file, myadios_group, myrank, trim(reg_name)//"nglob", nglob)

  if (nspec /= NSPEC_AB) call exit_mpi(myrank,'Error invalid nspec value in adios solver_data file')
  if (nglob /= NGLOB_AB) call exit_mpi(myrank,'Error invalid nglob value in adios solver_data file')

  ! reads mesh arrays
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "ibool", ibool(:,:,:,:))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xixstore", xixstore(:,:,:,:))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xiystore", xiystore(:,:,:,:))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xizstore", xizstore(:,:,:,:))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etaxstore", etaxstore(:,:,:,:))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etaystore", etaystore(:,:,:,:))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etazstore", etazstore(:,:,:,:))

  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammaxstore", gammaxstore(:,:,:,:))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammaystore", gammaystore(:,:,:,:))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammazstore", gammazstore(:,:,:,:))

  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "x_global", xstore(:))
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "y_global", ystore(:))
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "z_global", zstore(:))
  ! not needed so far
  !call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "idoubling", idoubling(:))
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

  read(IIN) xixstore
  read(IIN) xiystore
  read(IIN) xizstore
  read(IIN) etaxstore
  read(IIN) etaystore
  read(IIN) etazstore
  read(IIN) gammaxstore
  read(IIN) gammaystore
  read(IIN) gammazstore
  close(IIN)
#endif

  ! synchronizes
  call synchronize_all()

  ! get the location of the center of the elements
  do ispec = 1, NSPEC_AB

    ! gets point locations for this element
    DO_LOOP_IJK
      iglob = ibool(INDEX_IJK,ispec)
      xx0(INDEX_IJK,ispec) = xstore(iglob)
      yy0(INDEX_IJK,ispec) = ystore(iglob)
      zz0(INDEX_IJK,ispec) = zstore(iglob)
    ENDDO_LOOP_IJK

    cx0(ispec) = (xx0(1,1,1,ispec) + xx0(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cy0(ispec) = (yy0(1,1,1,ispec) + yy0(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
    cz0(ispec) = (zz0(1,1,1,ispec) + zz0(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL

    ! build jacobian
    DO_LOOP_IJK
      ! get derivatives of ux, uy and uz with respect to x, y and z
      xixl = xixstore(INDEX_IJK,ispec)
      xiyl = xiystore(INDEX_IJK,ispec)
      xizl = xizstore(INDEX_IJK,ispec)
      etaxl = etaxstore(INDEX_IJK,ispec)
      etayl = etaystore(INDEX_IJK,ispec)
      etazl = etazstore(INDEX_IJK,ispec)
      gammaxl = gammaxstore(INDEX_IJK,ispec)
      gammayl = gammaystore(INDEX_IJK,ispec)
      gammazl = gammazstore(INDEX_IJK,ispec)
      ! compute the jacobian
      jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                    - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                    + xizl*(etaxl*gammayl-etayl*gammaxl))
      jacobian(INDEX_IJK) = jacobianl
    ENDDO_LOOP_IJK

    ! integration factors
    if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
      ! adds GLL integration weights
      ! uses volume assigned to GLL points
      integ_factor(:,:,:,ispec) = jacobian(:,:,:) * wgll_cube(:,:,:)
    else
      ! no volume, only distance weights
      integ_factor(:,:,:,ispec) = 1.0_CUSTOM_REAL
    endif
  enddo

  ! user info
  if (myrank == 0) then
    ! mesh stats
    print *,'slice 0:'
    print *,'  mesh dimensions x min/max = ',minval(xstore),'/',maxval(xstore)
    print *,'  mesh dimensions y min/max = ',minval(ystore),'/',maxval(ystore)
    print *,'  mesh dimensions z min/max = ',minval(zstore),'/',maxval(zstore)
    print *
    print *,'  integration factor min/max =',minval(integ_factor),'/',maxval(integ_factor)
    print *
  endif

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
  allocate(nsearch_elements(NSPEC_AB),stat=ier)
  if (ier /= 0) stop 'Error allocating nsearch_elements array'
  nsearch_elements(:) = 0

  if (use_kdtree_search) then
    ! search by kd-tree
    ! user output
    if (myrank == 0) then
      print *
      print *,'using kd-tree search:'
      if (DO_SEARCH_ELLIP) then
        print *,'  search radius horizontal: ',r_search_dist_h * R_PLANET_KM,'km'
        print *,'  search radius vertical  : ',r_search_dist_v * R_PLANET_KM,'km'
      else
        print *,'  search sphere radius: ',r_search * R_PLANET_KM,'km'
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
      print *,'  search radius horizontal: ',sigma_h3 * R_PLANET_KM,'km'
      print *,'  search radius vertical  : ',sigma_v3 * R_PLANET_KM,'km'
      print *
    endif
  endif

  ! GPU setup
  if (USE_GPU) then
    ! user output
    if (myrank == 0) then
      print *
      print *,'preparing GPU arrays:'
      ! memory estimate
      ! data_smooth
      memory_size = NGLL3 * NSPEC_AB * nker * dble(CUSTOM_REAL)
      ! x_me/y_me/z_me + x_other/y_other/z_other +
      memory_size = memory_size + 6.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      ! normalisation + integ_factor
      memory_size = memory_size + 2.d0 * NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      ! kernel
      memory_size = memory_size + NGLL3 * NSPEC_AB * dble(CUSTOM_REAL)
      print *,'  minimum memory requested     : ',sngl(memory_size / 1024.d0 / 1024.d0),'MB per process'
      print *
    endif
    ! copies arrays onto GPU
    call prepare_smooth_gpu(Container,xx0,yy0,zz0,sigma_h2,sigma_v2,sigma_h3,sigma_v3,NSPEC_AB,nker)
  endif

  ! synchronizes
  call synchronize_all()

  if (myrank == 0) print *, 'start looping over elements and points for smoothing ...'

  ! synchronizes
  call synchronize_all()

  ! loop over all slices
  do inum = 1, num_slices

    iproc = islice(inum)

    ! user output
    if (myrank == 0) then
      print *
      print *,'slice ',inum,'out of ',num_slices
      print *,'  reading slice proc:',iproc
      print *
    endif

    ! debugging
    !if (DEBUG) then
    !  if (myrank == 0) then
    !    allocate(tmp_bk(NGLLX,NGLLY,NGLLZ,NSPEC_AB),stat=ier)
    !    if (ier /= 0) stop 'Error allocating array tmp_bk'
    !    tmp_bk(:,:,:,:) = 0.0_CUSTOM_REAL
    !  endif
    !endif

    ! sets up mesh locations
    if (iproc == myrank) then
      ! positions already read in, just copy
      ! element centers
      cx(:) = cx0(:)
      cy(:) = cy0(:)
      cz(:) = cz0(:)
      ! locations
      xx(:,:,:,:) = xx0(:,:,:,:)
      yy(:,:,:,:) = yy0(:,:,:,:)
      zz(:,:,:,:) = zz0(:,:,:,:)
      ! integration factor
      ! not needed as the first slice to loop over is set be always the myrank slice, and integ_factor has already been setup
      !integ_factor(:,:,:,:) = integ_factor0(:,:,:,:)
    else
      ! read in neighbor slice
#ifdef USE_ADIOS_INSTEAD_OF_MESH
      ! ADIOS
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "ibool", ibool(:,:,:,:))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "xixstore", xixstore(:,:,:,:))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "xiystore", xiystore(:,:,:,:))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "xizstore", xizstore(:,:,:,:))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "etaxstore", etaxstore(:,:,:,:))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "etaystore", etaystore(:,:,:,:))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "etazstore", etazstore(:,:,:,:))

      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "gammaxstore", gammaxstore(:,:,:,:))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "gammaystore", gammaystore(:,:,:,:))
      call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "gammazstore", gammazstore(:,:,:,:))

      call read_adios_array(myadios_file, myadios_group, iproc, nglob, trim(reg_name) // "x_global", xstore(:))
      call read_adios_array(myadios_file, myadios_group, iproc, nglob, trim(reg_name) // "y_global", ystore(:))
      call read_adios_array(myadios_file, myadios_group, iproc, nglob, trim(reg_name) // "z_global", zstore(:))
      ! not needed so far
      !call read_adios_array(myadios_file, myadios_group, iproc, nspec, trim(reg_name) // "idoubling", idoubling(:))
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

      read(IIN) xixstore
      read(IIN) xiystore
      read(IIN) xizstore
      read(IIN) etaxstore
      read(IIN) etaystore
      read(IIN) etazstore
      read(IIN) gammaxstore
      read(IIN) gammaystore
      read(IIN) gammazstore
      close(IIN)
#endif

      ! get the location of the center of the elements
      do ispec = 1, NSPEC_AB

        ! gets point locations for this element
        DO_LOOP_IJK
          iglob = ibool(INDEX_IJK,ispec)
          xx(INDEX_IJK,ispec) = xstore(iglob)
          yy(INDEX_IJK,ispec) = ystore(iglob)
          zz(INDEX_IJK,ispec) = zstore(iglob)
        ENDDO_LOOP_IJK

        ! calculate element center location
        cx(ispec) = (xx(1,1,1,ispec) + xx(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
        cy(ispec) = (yy(1,1,1,ispec) + yy(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL
        cz(ispec) = (zz(1,1,1,ispec) + zz(NGLLX,NGLLY,NGLLZ,ispec)) * 0.5_CUSTOM_REAL

        ! build jacobian
        DO_LOOP_IJK
          ! get derivatives of ux, uy and uz with respect to x, y and z
          xixl = xixstore(INDEX_IJK,ispec)
          xiyl = xiystore(INDEX_IJK,ispec)
          xizl = xizstore(INDEX_IJK,ispec)
          etaxl = etaxstore(INDEX_IJK,ispec)
          etayl = etaystore(INDEX_IJK,ispec)
          etazl = etazstore(INDEX_IJK,ispec)
          gammaxl = gammaxstore(INDEX_IJK,ispec)
          gammayl = gammaystore(INDEX_IJK,ispec)
          gammazl = gammazstore(INDEX_IJK,ispec)
          ! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))
          jacobian(INDEX_IJK) = jacobianl
        ENDDO_LOOP_IJK

        ! integration factors
        if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
          ! adds GLL integration weights
          ! uses volume assigned to GLL points
          integ_factor(:,:,:,ispec) = jacobian(:,:,:) * wgll_cube(:,:,:)
        else
          ! no volume, only distance weights
          integ_factor(:,:,:,ispec) = 1.0_CUSTOM_REAL
        endif
      enddo

      ! user info
      if (myrank == 0) then
        ! mesh stats
        print *,'slice ',iproc,':'
        print *,'  mesh dimensions x min/max = ',minval(xstore),'/',maxval(xstore)
        print *,'  mesh dimensions y min/max = ',minval(ystore),'/',maxval(ystore)
        print *,'  mesh dimensions z min/max = ',minval(zstore),'/',maxval(zstore)
        print *
        print *,'  integration factor min/max =',minval(integ_factor),'/',maxval(integ_factor)
        print *
      endif
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
      kernel_output_names(iker) = varname
      ! reads kernel values
      call read_adios_array(myadios_val_file, myadios_val_group, iproc, nspec, trim(varname), kernel(:,:,:,:,iker))
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
    if (use_kdtree_search) then
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

    ! pre-determines number of search elements
    if (use_kdtree_search) then
      ! gets search range
      num_elem_local_max = 0
      do ispec = 1, NSPEC_AB
        ! element center position
        center_x0 = cx0(ispec)
        center_y0 = cy0(ispec)
        center_z0 = cz0(ispec)

        ! initializes
        num_elem_local = 0
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

        ! debug output
        if (DEBUG) then
          if (myrank == 0) then
            ! user info
            if (ispec < 10) then
              print *,'  total number of search elements: ',num_elem_local,'ispec',ispec
            endif
          endif
        endif

        ! sets n-search number of nodes
        nsearch_elements(ispec) = num_elem_local

        ! maximum
        if (num_elem_local_max < num_elem_local) num_elem_local_max = num_elem_local
      enddo

      ! user output
      if (myrank == 0) then
        print *,'  maximum number of local search elements: ',num_elem_local_max
        print *
      endif

      ! allocates search array
      allocate(kdtree_search_index(num_elem_local_max),stat=ier)
      if (ier /= 0) stop 'Error allocating array kdtree_search_index'
      kdtree_search_index(:) = 0

      ! debug output
      if (DEBUG) then
        if (myrank == 0) then
          ! file output
          ! outputs search elements with integer flags
          allocate(ispec_flag(NSPEC_AB))
          ispec_flag(:) = 0

          ! debug element (tmp_ispec_dbg)
          num_elem_local = nsearch_elements(tmp_ispec_dbg)
          xyz_target(1) = cx0(ispec)
          xyz_target(2) = cy0(ispec)
          xyz_target(3) = cz0(ispec)

          ! finds closest point in target chunk
          kdtree_search_index(:) = 0
          if (DO_SEARCH_ELLIP) then
            ! (within ellipsoid)
            call kdtree_get_nearest_n_neighbors_ellip(xyz_target,r_search_dist_v,r_search_dist_h,num_elem_local)
          else
            ! (within search sphere)
            call kdtree_get_nearest_n_neighbors(xyz_target,r_search,num_elem_local)
          endif

          ! sets flags
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
      ! brute-force search
      ! always loop over whole mesh slice
      nsearch_elements(:) = NSPEC_AB
    endif

    ! smooth
    if (USE_GPU) then
      ! on GPU
      call compute_smooth_gpu(Container,kernel,integ_factor,xx,yy,zz,NSPEC_AB,USE_VECTOR_DISTANCE_FOR_SMOOTHING)
    else
      ! on CPU
      call compute_smooth(nker,bk,tk, &
                          kernel,integ_factor, &
                          cx0,cy0,cz0,cx,cy,cz, &
                          xx0,yy0,zz0,xx,yy,zz, &
                          sigma_v2,sigma_h2, &
                          sigma_v3,sigma_h3, &
                          r_search,r_search_dist_v,r_search_dist_h, &
                          nsearch_elements,use_kdtree_search, &
                          NSTEP_PERCENT_INFO)
    endif

    ! debug output
    !if (DEBUG) then
    !  if (myrank == 0) deallocate(tmp_bk)
    !endif

    ! deletes search tree nodes
    if (use_kdtree_search) then
      ! frees search index array
      deallocate(kdtree_search_index)
      call kdtree_delete()
    endif

    if (myrank == 0) print *

  enddo     ! islice

  ! user output
  if (myrank == 0) then
    print *
    print *,'normalizing smoothed kernel values...'
    print *
  endif

  ! frees memory
  if (use_kdtree_search) then
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

  ! outputs infos for normalizes/scaling factor
  if (myrank == 0) then
    if (.not. USE_GPU) then
      ! arrays bk,tk only available for CPU-only calculations
      print *, 'Scaling values:'
      print *, '  bk min/max = ',minval(bk),maxval(bk)
      print *, '  theoretical norm = ',norm
      do iker = 1,nker
        print *, '  ',trim(kernel_names(iker)),' tk min/max = ',minval(tk(:,:,:,:,iker)),maxval(tk(:,:,:,:,iker))
      enddo
      print *
    endif
  endif

  ! compute the smoothed kernel values
  if (USE_GPU) then
    ! on GPU
    call get_smooth_gpu(Container,kernel_smooth)
  else
    ! on CPU
    do ispec = 1, NSPEC_AB

      ! avoids division by zero
      DO_LOOP_IJK
        if (abs(bk(INDEX_IJK,ispec)) < 1.d-24) bk(INDEX_IJK,ispec) = 1.0_CUSTOM_REAL
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
  endif

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
    kernel_output_name = kernel_output_names(iker)
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call define_adios_global_array1D(myadios_val_group, group_size_inc, local_dim, '', &
                                     trim(kernel_output_name), kernel_smooth(:, :, :, :, iker))
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
    kernel_output_name = kernel_output_names(iker)
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(myadios_val_file, myadios_val_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(kernel_output_name), kernel_smooth(:,:,:,:,iker))
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

  ! frees arrays
  deallocate(ibool,idoubling,ispec_is_tiso)
  deallocate(xstore,ystore,zstore)
  deallocate(xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore)
  deallocate(bk,integ_factor)
  deallocate(xx0,yy0,zz0,xx,yy,zz,cx0,cy0,cz0,cx,cy,cz)
  deallocate(kernel,kernel_smooth,tk)

  ! timing
  if (myrank == 0) then
    tCPU = wtime() - time_start_all
    ! format time
    int_tCPU = int(tCPU)
    ihours = int_tCPU / 3600
    iminutes = (int_tCPU - 3600*ihours) / 60
    iseconds = int_tCPU - 3600*ihours - 60*iminutes
    write(*,*)
    write(*,*) 'Elapsed time in seconds  = ',tCPU
    write(*,"(' Elapsed time in hh:mm:ss = ',i6,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
    write(*,*)
  endif

  ! user output
  if (myrank == 0) print *, 'all done'

  ! stop all the processes, and exit
  call finalize_mpi()

end program smooth_sem_globe

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_smooth(nker,bk,tk, &
                            kernel,integ_factor, &
                            cx0,cy0,cz0,cx,cy,cz, &
                            xx0,yy0,zz0,xx,yy,zz, &
                            sigma_v2,sigma_h2, &
                            sigma_v3,sigma_h3, &
                            r_search,r_search_dist_v,r_search_dist_h, &
                            nsearch_elements,use_kdtree_search, &
                            NSTEP_PERCENT_INFO)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGLLCUBE,PI,myrank, &
    USE_QUADRATURE_RULE_FOR_SMOOTHING,USE_VECTOR_DISTANCE_FOR_SMOOTHING

  use postprocess_par, only: NSPEC_CRUST_MANTLE

  use kdtree_search

  implicit none

  integer, intent(in) :: nker
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),intent(inout) :: bk        ! normalization factors
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,nker),intent(inout) ::  tk  ! kernel values

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,nker), intent(in) :: kernel
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), intent(in) :: integ_factor

  ! center positions
  real(kind=CUSTOM_REAL), dimension(NSPEC_CRUST_MANTLE), intent(in) :: cx0,cy0,cz0,cx,cy,cz
  ! GLL point positions
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), intent(in) :: xx0,yy0,zz0,xx,yy,zz

  ! smoothing radii
  real(kind=CUSTOM_REAL), intent(in) :: sigma_v2,sigma_h2
  real(kind=CUSTOM_REAL), intent(in) :: sigma_h3,sigma_v3
  ! search radii
  double precision,intent(in) :: r_search,r_search_dist_v,r_search_dist_h

  ! number of search elements
  integer,dimension(NSPEC_CRUST_MANTLE), intent(in) :: nsearch_elements
  logical, intent(in) :: use_kdtree_search

  ! user output infos
  integer, intent(in) :: NSTEP_PERCENT_INFO

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: exp_val

  real(kind=CUSTOM_REAL) :: val,val_Gaussian
  real(kind=CUSTOM_REAL) :: center_x0, center_y0, center_z0
  real(kind=CUSTOM_REAL) :: center_x, center_y, center_z
  real(kind=CUSTOM_REAL) :: contrib_weights,contrib_kernel

  ! for distance calc
  real(kind=CUSTOM_REAL) :: sigma_h2_inv,sigma_v2_inv
  real(kind=CUSTOM_REAL) :: sigma_h3_sq,sigma_v3_sq
  real(kind=CUSTOM_REAL) :: dist_h,dist_v
  real(kind=CUSTOM_REAL) :: x0, y0, z0
  real(kind=CUSTOM_REAL) :: x1, y1, z1
  real(kind=CUSTOM_REAL) :: r0,r1
  real(kind=CUSTOM_REAL) :: r0_squared,r1_squared
  real(kind=CUSTOM_REAL) :: theta,ratio,alpha
  real(kind=CUSTOM_REAL) :: vx,vy,vz

  double precision,dimension(3) :: xyz_target

  integer :: ispec,ispec2,iker
  ! search elements
  integer :: ielem,num_elem_local

  integer :: NELEM_OUTPUT_INFO,NELEM_MAX_OUTPUT_INFO
  integer :: NSPEC_AB

  ! timing
  double precision :: time_start
  double precision :: tCPU
  double precision, external :: wtime

! switches do-loops between: do k=1,NGLLZ; do j=1,NGLLY; do i=1,NGLLX <-> do ijk=1,NGLLCUBE
#ifdef FORCE_VECTORIZATION
  integer :: ijk
  integer :: ijk2
#  define INDEX_IJK2  ijk2,1,1
#  define DO_LOOP_IJK2  do ijk2 = 1,NGLLCUBE
#  define ENDDO_LOOP_IJK2  enddo                  ! NGLLCUBE
#else
  integer :: i,j,k
  integer :: ii,jj,kk
#  define INDEX_IJK2  ii,jj,kk
#  define DO_LOOP_IJK2  do kk = 1,NGLLZ; do jj = 1,NGLLY; do ii = 1,NGLLX
#  define ENDDO_LOOP_IJK2  enddo; enddo; enddo    ! NGLLZ,NGLLY,NGLLX
#endif

  !---------------------
  ! Parameters

  ! kd-tree search uses either a spherical or ellipsoid search
  logical,parameter :: DO_SEARCH_ELLIP = .true.

  real(kind=CUSTOM_REAL),parameter :: PI2 = real(PI * PI,kind=CUSTOM_REAL) ! squared

  ! debugging
  !logical, parameter :: DEBUG = .false.
  ! boundary point: 910, interior point: 2382, surface point: 5128
  !integer ,parameter :: tmp_ispec_dbg = 5128
  !---------------------

  !takes region 1 kernels
  NSPEC_AB = NSPEC_CRUST_MANTLE

  ! for time output
  NELEM_OUTPUT_INFO = max(int(1.0*NSPEC_AB/NSTEP_PERCENT_INFO),1)
  NELEM_MAX_OUTPUT_INFO = NSPEC_AB - int(0.5*NSPEC_AB/NSTEP_PERCENT_INFO)

  ! helper variables
  sigma_h2_inv = 1.0_CUSTOM_REAL / sigma_h2
  sigma_v2_inv = 1.0_CUSTOM_REAL / sigma_v2

  sigma_h3_sq = sigma_h3 * sigma_h3  ! squared
  sigma_v3_sq = sigma_v3 * sigma_v3

  ! timing
  time_start = wtime()

  ! loops over elements to be smoothed in the current slice
  do ispec = 1, NSPEC_AB

    ! user info about progress
    if (myrank == 0) then
      if (ispec == NSPEC_AB) then
        tCPU = wtime() - time_start
        print *,'    ',100.0,' % elements done - Elapsed time in seconds = ',tCPU
      else if (mod(ispec-1,NELEM_OUTPUT_INFO) == 0) then
        if (ispec < NELEM_MAX_OUTPUT_INFO) then
          tCPU = wtime() - time_start
          print *,'    ',int((ispec-1) / NELEM_OUTPUT_INFO) * (100.0 / NSTEP_PERCENT_INFO), &
                  ' % elements done - Elapsed time in seconds = ',tCPU
        endif
      endif
    endif

    ! search range
    num_elem_local = nsearch_elements(ispec)

    ! element center position
    center_x0 = cx0(ispec)
    center_y0 = cy0(ispec)
    center_z0 = cz0(ispec)

    ! sets number of elements to loop over
    if (use_kdtree_search) then
      ! sets search array
      if (num_elem_local > 0) then
        ! sets n-search number of nodes
        kdtree_search_num_nodes = num_elem_local
        ! finds closest point in target chunk
        xyz_target(1) = center_x0
        xyz_target(2) = center_y0
        xyz_target(3) = center_z0
        if (DO_SEARCH_ELLIP) then
          ! (within ellipsoid)
          call kdtree_get_nearest_n_neighbors_ellip(xyz_target,r_search_dist_v,r_search_dist_h,num_elem_local)
        else
          ! (within search sphere)
          call kdtree_get_nearest_n_neighbors(xyz_target,r_search,num_elem_local)
        endif
      endif
    endif

    ! --- only double loop over the elements in the search radius ---
    do ielem = 1, num_elem_local
      ! search element
      if (use_kdtree_search) then
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

      ! calculates horizontal and vertical distance between two element centers
      call get_distance_vec_squared(dist_h,dist_v,center_x0,center_y0,center_z0,center_x,center_y,center_z)

      ! takes only elements in search radius
      ! note: distances and sigmah, sigmav are normalized by R_PLANET_KM
      if (dist_h > sigma_h3_sq .or. dist_v > sigma_v3_sq) cycle

      ! loops over all GLL points of reference element in current slice (ispec) and search elements (ispec2)
      DO_LOOP_IJK

        ! reference location
        ! current point (i,j,k,ispec) location, Cartesian coordinates
        x0 = xx0(INDEX_IJK,ispec)
        y0 = yy0(INDEX_IJK,ispec)
        z0 = zz0(INDEX_IJK,ispec)

        ! calculate weights based on Gaussian smoothing

        ! calculates horizontal and vertical distance between two element centers
        !call smoothing_weights_vec(x0,y0,z0,sigma_h2,sigma_v2,exp_val, &
        !                           xx(:,:,:,ispec2),yy(:,:,:,ispec2),zz(:,:,:,ispec2))

        ! same as calling routine smoothing_weights_vec(),
        ! but inlined here to improve speed - a little..

        ! length of first position vector (squared)
        r0_squared = x0*x0 + y0*y0 + z0*z0

        if (USE_VECTOR_DISTANCE_FOR_SMOOTHING) then
          ! vector approximation (fast computation): neglects curvature

          ! length of first position vector
          r0 = sqrt( r0_squared )

          DO_LOOP_IJK2
            ! point in second slice
            x1 = xx(INDEX_IJK2,ispec2)
            y1 = yy(INDEX_IJK2,ispec2)
            z1 = zz(INDEX_IJK2,ispec2)

            ! calculates horizontal and vertical distance between two element centers
            !call get_distance_vec_squared(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

            ! same as calling routine get_distance_vec_squared(),
            ! but to help compiler optimize the loop calculation...

            ! without explicit function calls to help compiler optimize loops
            ! length of position vector
            r1_squared = x1*x1 + y1*y1 + z1*z1
            r1 = sqrt( r1_squared )

            ! vertical distance (squared)
            dist_v = (r1 - r0)*(r1 - r0)

            ! horizontal distance
            ! length of vector from point 0 to point 1
            ! assuming small earth curvature  (since only for neighboring elements)
            ! scales r0 to have same length as r1
            alpha = r1 / r0
            vx = alpha * x0
            vy = alpha * y0
            vz = alpha * z0

            ! vector in horizontal between new r0 and r1
            vx = x1 - vx
            vy = y1 - vy
            vz = z1 - vz

            ! distance is vector length
            !dist_h = sqrt( vx*vx + vy*vy + vz*vz )
            ! or squared:
            dist_h = vx*vx + vy*vy + vz*vz

            ! Gaussian function
            val = - dist_h * sigma_h2_inv - dist_v * sigma_v2_inv

            ! limits to single precision
            if (val < - 86.0_CUSTOM_REAL) then
              ! smaller than numerical precision: exp(-86) < 1.e-37
              val_Gaussian = 0.0_CUSTOM_REAL
            else
              val_Gaussian = exp(val)
            endif

            exp_val(INDEX_IJK2) = val_Gaussian

          ENDDO_LOOP_IJK2

        else
          ! w/ exact epicentral distance calculation
          DO_LOOP_IJK2
            ! point in second slice
            x1 = xx(INDEX_IJK2,ispec2)
            y1 = yy(INDEX_IJK2,ispec2)
            z1 = zz(INDEX_IJK2,ispec2)

            ! calculates horizontal and vertical distance between two element centers
            !call get_distance_vec_squared(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

            ! same as calling routine get_distance_vec_squared(),
            ! but to help compiler optimize the loop calculation...

            ! without explicit function calls to help compiler optimize loops
            ! length of position vector
            r1_squared = x1*x1 + y1*y1 + z1*z1

            ! vertical distance (squared)
            ! dist_v = (r1 - r0)*(r1 - r0)
            !        = r1**2 + r0**2 - 2 * r0 * r1
            !        = r1**2 + r0**2 - 2 * sqrt( r0**2 ) * sqrt( r1**2 )
            !        = r1**2 + r0**2 - 2 * sqrt( r0**2 * r1**2 )
            !        = r1**2 + r0**2 - 2 * alpha
            !          with alpha = sqrt( r0**2 * r1**2 ) = r0 * r1
            ! this avoids using sqrt() function too often which is costly
            alpha = sqrt( r0_squared * r1_squared )
            dist_v = r1_squared + r0_squared - 2.0_CUSTOM_REAL * alpha

            ! only for flat earth with z in depth:
            !  dist_v = sqrt( (cz(ispec2)-cz0(ispec))** 2)

            ! epicentral distance
            ! (accounting for spherical curvature)
            ! calculates distance of circular segment
            ! angle between r0 and r1 in radian
            ! given by dot-product of two vectors
            !ratio = (x0*x1 + y0*y1 + z0*z1) / (r0 * r1)
            !      = (x0*x1 + y0*y1 + z0*z1) / sqrt(r0**2 * r1**2) -> alpha = r0 * r1
            if (alpha > 0.0_CUSTOM_REAL) then
              ratio = (x0*x1 + y0*y1 + z0*z1) / alpha
            else
              ratio = 1.0_CUSTOM_REAL
            endif

            ! checks boundaries of ratio (due to numerical inaccuracies)
            if (ratio >= 1.0_CUSTOM_REAL) then
              ! ratio = 1.0_CUSTOM_REAL
              ! -> acos( 1 ) = 0
              ! -> dist_h = 0
              dist_h = 0.0_CUSTOM_REAL
            else if (ratio <= -1.0_CUSTOM_REAL) then
              ! ratio = -1.0_CUSTOM_REAL
              ! -> acos ( -1 ) = PI
              ! -> dist_h = (r1 * pi) * (r1 * pi ) = r1**2 * PI**2
              dist_h = r1_squared * PI2
            else
              theta = acos( ratio )
              ! segment length at heigth of r1 (squared)
              ! dist_h = (r1 * theta)*(r1 * theta) = r1**2 * theta**2
              dist_h = r1_squared * (theta*theta)
            endif

            ! Gaussian function
            val = - dist_h * sigma_h2_inv - dist_v * sigma_v2_inv

            ! limits to single precision
            if (val < - 86.0_CUSTOM_REAL) then
              ! smaller than numerical precision: exp(-86) < 1.e-37
              val_Gaussian = 0.0_CUSTOM_REAL
            else
              val_Gaussian = exp(val)
            endif

            exp_val(INDEX_IJK2) = val_Gaussian

          ENDDO_LOOP_IJK2
        endif

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
        if (USE_QUADRATURE_RULE_FOR_SMOOTHING) then
          !exp_val(:,:,:) = exp_val(:,:,:) * integ_factor(:,:,:,ispec2)
          ! explicit loop
          DO_LOOP_IJK2
            exp_val(INDEX_IJK2) = exp_val(INDEX_IJK2) * integ_factor(INDEX_IJK2,ispec2)
          ENDDO_LOOP_IJK2
        endif

        ! intrinsic sum
        !contrib_weights = sum(exp_val(:,:,:))
        ! explicit loop
        contrib_weights = 0.0_CUSTOM_REAL
        DO_LOOP_IJK2
          contrib_weights = contrib_weights + exp_val(INDEX_IJK2)
        ENDDO_LOOP_IJK2

        ! normalization, integrated values of Gaussian smoothing function
        bk(INDEX_IJK,ispec) = bk(INDEX_IJK,ispec) + contrib_weights

        ! adds contribution of element ispec2 to smoothed kernel values
        do iker = 1,nker
          ! kernel contributions
          !contrib_kernel = sum(exp_val(:,:,:) * kernel(:,:,:,ispec2,iker))
          ! explicit loop
          contrib_kernel = 0.0_CUSTOM_REAL
          DO_LOOP_IJK2
            contrib_kernel = contrib_kernel + exp_val(INDEX_IJK2) * kernel(INDEX_IJK2,ispec2,iker)
          ENDDO_LOOP_IJK2

          ! new kernel value
          tk(INDEX_IJK,ispec,iker) = tk(INDEX_IJK,ispec,iker) + contrib_kernel
        enddo

        ! checks number
        !if (isNaN(tk(INDEX_IJK,ispec,1))) then
        !  print *,'Error tk NaN: ',tk(INDEX_IJK,ispec,1)
        !  print *,'rank:',myrank
        !  print *,'INDEX_IJK,ispec:',INDEX_IJK,ispec
        !  print *,'tk: ',tk(INDEX_IJK,ispec,1),'bk:',bk(INDEX_IJK,ispec)
        !  print *,'sum exp_val: ',sum(exp_val(:,:,:)),'sum factor:',sum(integ_factor(:,:,:,ispec2))
        !  print *,'sum kernel:',sum(kernel(:,:,:,ispec2))
        !  call exit_mpi(myrank, 'Error NaN')
        !endif

      ENDDO_LOOP_IJK

    enddo ! (ielem)

    ! debug output
    !if (DEBUG) then
    !  if (myrank == 0) then
    !    ! debug element (tmp_ispec_dbg)
    !    if (ispec == tmp_ispec_dbg) then
    !      ! outputs Gaussian weighting function
    !      ! writes out vtk file for debugging
    !      write(filename,'(a,i4.4,a,i6.6)') trim(LOCAL_PATH)//'/search_elem',tmp_ispec_dbg,'_Gaussian_proc',iproc
    !      call write_VTK_data_gll_cr(NSPEC_AB,NGLOB_AB,xstore,ystore,zstore,ibool,tmp_bk,filename)
    !      print *,'file written: ',trim(filename)//'.vtk'
    !    endif
    !  endif
    !endif

  enddo   ! (ispec)

  end subroutine compute_smooth

!
! -----------------------------------------------------------------------------
!

! subroutine put in this file to help compiler inlining code...
! not used, but here for comparison

  subroutine get_distance_vec_squared(dist_h,dist_v,x0,y0,z0,x1,y1,z1)

! returns distances in radial and horizontal direction (squared)

  use constants, only: CUSTOM_REAL,PI,USE_VECTOR_DISTANCE_FOR_SMOOTHING

  implicit none

  real(kind=CUSTOM_REAL),intent(out) :: dist_h,dist_v
  real(kind=CUSTOM_REAL),intent(in) :: x0,y0,z0,x1,y1,z1

  ! local parameters
  real(kind=CUSTOM_REAL) :: r0_squared,r1_squared
  real(kind=CUSTOM_REAL) :: theta,ratio,alpha
  real(kind=CUSTOM_REAL) :: vx,vy,vz,r0,r1

  real(kind=CUSTOM_REAL),parameter :: PI2 = real(PI * PI,kind=CUSTOM_REAL) ! squared

  ! note: instead of distance we use distance squared to avoid too many sqrt() operations

  ! spherical Earth: distance in radial direction
  ! radius (squared)
  r0_squared = x0*x0 + y0*y0 + z0*z0    ! length of first position vector (squared)
  r1_squared = x1*x1 + y1*y1 + z1*z1

  if (USE_VECTOR_DISTANCE_FOR_SMOOTHING) then
    ! vector approximation (fast computation): neglects curvature

    ! only for flat earth with z in depth: dist_v = sqrt( (z1-z0)** 2)
    !dist_v = sqrt( (center_z-center_z0)** 2)
    ! or squared:
    !dist_v = (center_z - center_z0)** 2

    ! spherical Earth: distance in radial direction
    ! radius
    r0 = sqrt( r0_squared )    ! length of first position vector
    r1 = sqrt( r1_squared )

    ! vertical distance (squared)
    dist_v = (r1 - r0)*(r1 - r0)

    ! horizontal distance
    ! length of vector from point 0 to point 1
    ! assuming small earth curvature  (since only for neighboring elements)
    ! scales r0 to have same length as r1
    alpha = r1 / r0
    vx = alpha * x0
    vy = alpha * y0
    vz = alpha * z0

    ! vector in horizontal between new r0 and r1
    vx = x1 - vx
    vy = y1 - vy
    vz = z1 - vz

    ! distance is vector length
    !dist_h = sqrt( vx*vx + vy*vy + vz*vz )
    ! or squared:
    dist_h = vx*vx + vy*vy + vz*vz
  else
    ! epicentral distance
    ! (accounting for spherical curvature)

    ! vertical distance (squared)
    ! dist_v = (r1 - r0)*(r1 - r0)
    !        = r1**2 + r0**2 - 2 * alpha
    !          with alpha = sqrt( r0**2 * r1**2 ) = r0 * r1
    ! this avoids using sqrt() function too often which is costly
    alpha = sqrt( r0_squared * r1_squared )
    dist_v = r1_squared + r0_squared - 2.0_CUSTOM_REAL * alpha

    ! epicentral dsitance
    ! calculates distance of circular segment
    ! angle between r0 and r1 in radian
    ! given by dot-product of two vectors
    !ratio = (center_x0*center_x + center_y0*center_y + center_z0*center_z)/(r0 * r1)
    if (alpha > 0.0_CUSTOM_REAL) then
      ratio = (x0*x1 + y0*y1 + z0*z1) / alpha
    else
      ratio = 1.0_CUSTOM_REAL
    endif

    ! checks boundaries of ratio (due to numerical inaccuracies)
    if (ratio >= 1.0_CUSTOM_REAL) then
      ! ratio = 1.0_CUSTOM_REAL
      ! -> acos(1) = 0
      ! -> dist_h = 0
      dist_h = 0.0_CUSTOM_REAL
    else if (ratio <= -1.0_CUSTOM_REAL) then
      ! ratio = -1.0_CUSTOM_REAL
      ! -> acos(-1) = PI
      ! -> dist_h = r1**2 * PI**2
      dist_h = r1_squared * PI2
    else
      theta = acos( ratio )
      ! segment length at heigth of r1 (squared)
      dist_h = r1_squared * (theta*theta)
    endif
  endif

  end subroutine get_distance_vec_squared
