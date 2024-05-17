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

!------------------------------------------------------------------------------
! creates a horizontal/vertical cross-section at a given depth/great-circle path
!
! inputs:
!   - model parameter                       e.g. vsv, vsh   (-> proc***_reg1_vsv.bin)
!
!   - source model directory                 e.g. MODEL_M10/ (directory where proc***_reg1_vsv.bin lies)
!
! outputs:
!   - cross-section ascii-file into output_dir/
!
!
! needs:
!   - topo files from databases mesh            e.g. DATABASES_MPI/proc000001_reg1_solver_data.bin
!
!
! usage: xcreate_cross_section param section-param mesh-dir/ model-dir/ output-dir/ topoography-flag ellipticity-flag
!
!   for example
!
!   > mpirun -np 24 ./bin/xcreate_cross_section vpv -H 100.0,1.0 DATABASES_MPI/ DATABASES_MPI/ OUTPUT_FILES/ 1 1
!
!   creates a horizontal cross-section for mesh file proc**_vpv.bin at 100 km depth,
!   with a 1 x 1 degree regular latitude-longitude grid, taking into account topography and ellipticity
!
!------------------------------------------------------------------------------


  program cross_section

  use constants, only: myrank

  use constants, only: SIZE_INTEGER, &
    TWO_PI,R_UNIT_SPHERE, &
    NGNOD,CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
    GAUSSALPHA,GAUSSBETA, &
    IIN,IOUT,MAX_STRING_LEN, &
    RADIANS_TO_DEGREES,SMALLVAL,HUGEVAL

  use shared_parameters, only: R_PLANET_KM,LOCAL_PATH

  use postprocess_par, only: MAX_KERNEL_NAMES, &
    NCHUNKS_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,NPROCTOT_VAL,NEX_XI_VAL, &
    NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE

  use kdtree_search, only: kdtree_setup,kdtree_set_verbose,kdtree_delete,kdtree_find_nearest_neighbor, &
    kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

  implicit none

  !-------------------------------------------------------------------
  ! USER PARAMETERS

  ! region code, by default uses crust/mantle region (reg1)
  character(len=16) :: reg = '_reg1_'

  !-------------------------------------------------------------------

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLY) :: wygll
  double precision, dimension(NGLLZ) :: wzgll

  ! topology of the control points of the surface element
  integer :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! typical element size in old mesh used as search radius
  double precision :: typical_size

  ! mesh
  ! combined slices mesh
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: xstore,ystore,zstore
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: model
  integer,dimension(:,:,:,:),allocatable :: ibool

  ! target mesh points (on cross-section surface)
  integer :: nglob_target
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x2, y2, z2
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: model2
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: model_distance2
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: model_diff,model_pert

  ! input arguments
  character(len=MAX_STRING_LEN) :: arg
  character(len=MAX_STRING_LEN) :: dir_topo
  character(len=MAX_STRING_LEN) :: input_model_dir,output_dir

  integer :: i,j,k,iglob,ispec,ier,iflag
  integer :: nspec, nglob
  integer(kind=8) :: arraysize

  ! model
  character(len=16) :: fname
  character(len=MAX_STRING_LEN) :: m_file,filename
  character(len=MAX_STRING_LEN) :: solver_file

  ! MPI processes
  integer :: sizeprocs

  ! nodes search
  integer :: inodes
  real(kind=CUSTOM_REAL) :: model_maxdiff,min_val,max_val
  real(kind=CUSTOM_REAL) :: val

  ! horizontal latitude-longitude grid
  integer :: nlat,nlon
  ! vertical latitude-longitude grid
  integer :: ndepth,nsec

  double precision :: lat0,dlat
  double precision :: lon0,dlon
  double precision :: depth0,ddepth,depth_min,depth_max
  double precision :: dincr
  double precision :: r_top,r_bottom

  double precision :: v_lat1,v_lon1,v_lat2,v_lon2
  double precision :: theta1,phi1,theta2,phi2
  double precision :: epidist

  ! statistics
  double precision :: m_avg_total,point_avg_total

  ! cross-section infos
  integer :: section_type
  integer :: nsection_params
  character(len=MAX_STRING_LEN),dimension(:),allocatable :: param_args

  ! depth with respect to with topography
  logical :: TOPOGRAPHY

  ! depth considering ellipticity
  logical :: ELLIPTICITY

  ! timer MPI
  double precision :: tCPU,time_start
  double precision, external :: wtime

  ! minimum radius for CMB (from AK135, see get_model_parameters.F90)
  double precision, parameter :: RCMB_LIMIT_KM   = 3479.5d0

  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks program arguments
  if (myrank == 0) then
    do i = 1,8
      call get_command_argument(i,arg)
      ! usage info
      if (len_trim(arg) == 0) then
        if (myrank == 0) then
          print *
          print *,' Usage: xcreate_cross_section param section-param mesh-dir/ model-dir/ output-dir/ ' // &
                  'topography-flag ellipticity-flag'
          print *
          print *,' with'
          print *,'   param         - model parameter name (e.g. vpv)'
          print *
          print *,'   section-param - cross-section parameters, use format:'
          print *,'        -H depth,delta-incr'
          print *,'        for horizontal cross-section at depth (in km) '
          print *,'        with increment delta-incr for regular latitude/longitudes (in degrees)'
          print *
          print *,'        -V lat1,lon1,lat2,lon2,delta-incr,depth-incr,(depth_min),(depth_max)'
          print *,'        for vertical cross-section along great-circle path '
          print *,'        at points (lat1,lon2) and (lat2,lon2) '
          print *,'        with increment delta-incr for regular great-circle section'
          print *,'        and depth-incr for depth increments (in km)'
          print *,'        optionally, depth_min and depth_max can be set for minimum/maximum depths (given in km)'
          print *
          print *,'   mesh-dir/     - mesh directory with topology files (e.g. proc***_solver_data.bin)'
          print *,'   model-dir/    - directoy which holds model files (e.g. proc***_vpv.bin)'
          print *,'   output-dir/   - output directory with topology files (e.g. proc***_solver_data.bin)'
          print *
          print *,'   topography-flag  - depth will be taken with respect to surface topography (0 == off / 1 == on);'
          print *,'                      if no topography is used, then depth is with respect to sea-level'
          print *,'   ellipticity-flag - depth will consider Earth ellipticity (0 == off / 1 == on);'
          print *,'                      if no ellipticity is used, Earth shape is assumed to be perfectly spherical'
          print *
        endif
        stop ' Reenter command line options'
      endif
    enddo
  endif
  call synchronize_all()

  ! allocates arrays
  allocate(param_args(MAX_KERNEL_NAMES),stat=ier)
  if (ier /= 0) stop 'Error allocating param_args array'
  param_args(:) = ''

  ! initializes cross section type
  section_type = 0
  depth0 = 0.d0
  ddepth = 0.d0
  dincr = 0.d0
  v_lat1 = 0.d0
  v_lon1 = 0.d0
  v_lat2 = 0.d0
  v_lon2 = 0.d0
  depth_min = 0.d0
  depth_max = R_PLANET_KM - RCMB_LIMIT_KM

  ! reads input arguments
  do i = 1, 8
    call get_command_argument(i,arg)
    ! assignes values
    select case (i)
    case (1)
      fname = trim(arg)
    case (2)
      ! section parameter format:
      !   -H depth,delta-incr                 for horizontal x-section
      !   -V lat1,lon1,lat2,lon2,delta-incr   for vertical x-section
      if (trim(arg) == "-H") then
        section_type = 0
      else if (trim(arg) == "-V") then
        section_type = 1
      else
        stop 'Error section-param not recognized, please use -H or -V ..'
      endif
    case (3)
      ! section parameters, comma-separated
      call parse_kernel_names(arg, param_args, nsection_params)
      if (section_type == 0) then
        ! horizontal cross-section
        if (nsection_params /= 2) &
          stop 'Error horizontal section format: -H depth,delta-incr'
        read(param_args(1),*) depth0
        read(param_args(2),*) dincr
        ! check
        if (depth0 > R_PLANET_KM) stop 'Error cross-section depth (given in km) must be smaller than Earth radius'
        if (dincr > 360.0) stop 'Error lat/lon increment (given in degrees) must be smaller than 360'
      else
        ! vertical cross-section
        if (nsection_params /= 6 .and. nsection_params /= 8) &
          stop 'Error vertical section format: -V lat1,lon1,lat2,lon2,delta-incr,depth-incr,(depth_min),(depth_max)'
        read(param_args(1),*) v_lat1
        read(param_args(2),*) v_lon1
        read(param_args(3),*) v_lat2
        read(param_args(4),*) v_lon2
        read(param_args(5),*) dincr
        read(param_args(6),*) ddepth
        if (nsection_params == 8) then
          read(param_args(7),*) depth_min
          read(param_args(8),*) depth_max
        endif
      endif
    case (4)
      dir_topo = trim(arg)
    case (5)
      input_model_dir = trim(arg)
    case (6)
      output_dir = trim(arg)
    case (7)
      read(arg,*) iflag
      if (iflag == 0) then
        TOPOGRAPHY = .false.
      else if (iflag == 1) then
        TOPOGRAPHY = .true.
      else
        stop 'Error topography-flag value is invalid, please use 0 for turning off, or 1 for on'
      endif
    case (8)
      read(arg,*) iflag
      if (iflag == 0) then
        ELLIPTICITY = .false.
      else if (iflag == 1) then
        ELLIPTICITY = .true.
      else
        stop 'Error ellipticity-flag value is invalid, please use 0 for turning off, or 1 for on'
      endif
    end select
  enddo

  ! kdtree search:
  ! searches closest element using internal GLL points

  ! console output
  if (myrank == 0) then
    print *
    if (section_type == 0) then
      print *,'horizontal cross-section:'
    else if (section_type == 1) then
      print *,'vertical cross-section:'
    else
      stop 'Error cross-section type invalid'
    endif
    print *
    print *,'model parameter: ',trim(fname)
    print *
    print *,'  mesh directory: ',trim(dir_topo)
    print *,' model directory: ',trim(input_model_dir)
    print *,'output directory: ',trim(output_dir)
    print *
  endif

  ! reads mesh parameters
  if (myrank == 0) then
    ! reads mesh_parameters.bin file from input1dir/
    LOCAL_PATH = dir_topo
    call read_mesh_parameters()
  endif
  ! broadcast parameters to all processes
  call bcast_mesh_parameters()

  ! user output
  if (myrank == 0) then
    print *,'mesh parameters (from mesh directory):'
    print *,'  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
    print *,'  NPROCTOT           = ',NPROCTOT_VAL
    print *
  endif

  ! checks number of processes
  ! note: must run with same number of process as new mesh was created
  if (sizeprocs /= NPROCTOT_VAL) then
    ! usage info
    if (myrank == 0) then
      print *, "this program must be executed in parallel with NPROCTOT_VAL = ",NPROCTOT_VAL,"processes"
      print *, "Invalid number of processes used: ", sizeprocs, " procs"
      print *
      print *, "Please run: mpirun -np ",NPROCTOT_VAL," ./bin/xcreate_cross_section .."
    endif
    call abort_mpi()
  endif

  ! console output
  if (myrank == 0) then
    print *,'mesh:  '
    print *,'  processors = ',NPROCTOT_VAL
    print *,'  nproc_eta / nproc_xi = ',NPROC_ETA_VAL,NPROC_XI_VAL
    print *,'  nex        = ',NEX_XI_VAL
    print *,'  nspec      = ',NSPEC_CRUST_MANTLE
    print *,'  nglob      = ',NGLOB_CRUST_MANTLE
    print *
    print *,'array size:'
    ! array size in bytes (note: the multiplication is split into two line to avoid integer arithmetic overflow)
    arraysize = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE
    arraysize = arraysize * NPROC_ETA_VAL * NPROC_XI_VAL * SIZE_INTEGER
    print *,'  ibool   = ',dble(arraysize)/1024./1024.,'MB'
    arraysize = NGLOB_CRUST_MANTLE
    arraysize = arraysize * NPROC_ETA_VAL * NPROC_XI_VAL * CUSTOM_REAL
    print *,'  x,y,z   = ',dble(arraysize)/1024./1024.,'MB'
    print *
    arraysize = NGLLX * NGLLY * NGLLZ * NSPEC_CRUST_MANTLE
    arraysize = arraysize * NPROC_ETA_VAL * NPROC_XI_VAL * CUSTOM_REAL
    print *,'  model   = ',dble(arraysize)/1024./1024.,'MB'
    print *
    print *,'total MPI processes: ',sizeprocs
    print *
    print *,'location search by : kd-tree search'
    print *,'  uses internal GLL points'
    print *
    if (TOPOGRAPHY) then
      print *,'uses topography'
      print *,'  depth will be determined from (mesh) elevation'
    else
      print *,'no topography'
      print *,'  depth given with respect to sea-level'
    endif
    print *
    if (ELLIPTICITY) then
      print *,'uses ellipticity'
      print *,'  radius and latitude of points will be determined for elliptical Earth'
    else
      print *,'no ellipticity'
      print *,'  radius and latitude of points with respect to a spherical Earth'
    endif
    print *
  endif
  call synchronize_all()

  ! checks temporary file creation, to see if we could write out new model
  if (myrank == 0) then
    write(m_file,'(a,i6.6,a)') trim(output_dir)// '/cross_section_'//trim(fname)//'.tmp'
    open(IOUT,file=trim(m_file),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(m_file)
      stop 'Error opening new output model file, please check if output directory exists...'
    endif
    close(IOUT,status='delete')
  endif
  call synchronize_all()

  ! timing
  time_start = wtime()

  ! GLL points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! define topology of the control element
  call hex_nodes_anchor_ijk(anchor_iax,anchor_iay,anchor_iaz)

  ! source mesh
  ! compute typical size of elements at the surface
  typical_size = TWO_PI * R_UNIT_SPHERE / ( 4.0 * NEX_XI_VAL )
  ! use 2 times the distance as a criterion for point detection
  typical_size = 2.0 * typical_size

  ! user output
  if (myrank == 0) then
    print *
    print *,'creating target cross-section points... '
    print *
  endif

  ! cross-section points
  if (section_type == 0) then
    ! horizontal cross-section
    ! user output
    if (myrank == 0) then
      print *,'  horizontal cross-section'
      print *
    endif

    ! uses a regular lat/lon grid
    !
    ! latitudes range [-90,90]
    lat0 = -90.d0
    dlat = dincr
    nlat = int(180.d0 / dlat) + 1
    ! longitudes range [0,360[
    lon0 = 0.d0
    dlon = dincr
    nlon = int(360.d0 / dlon)
    ! depths (in km)
    depth0 = depth0
    ddepth = 10.d0
    ndepth = 1

    ! total number of points
    ! (north/south pole will count as single points)
    nglob_target = (nlat - 2) * nlon * ndepth + 2

    ! user output
    if (myrank == 0) then
      print *,'  number of lat/lon points      = ',nlat,'/',nlon
      print *,'  total number of target points = ',nglob_target
      print *
    endif

    ! allocates arrays for target points
    allocate(x2(nglob_target), &
             y2(nglob_target), &
             z2(nglob_target), &
             model_distance2(nglob_target), &
             model2(nglob_target),stat=ier)
    if (ier /= 0) stop 'Error allocating target model point arrays'
    x2(:) = 0.0_CUSTOM_REAL; y2(:) = 0.0_CUSTOM_REAL; z2(:) = 0.0_CUSTOM_REAL
    model_distance2(:) = real(HUGEVAL,kind=CUSTOM_REAL)
    model2(:) = 0.0_CUSTOM_REAL

    ! creates cross-section points
    call set_horiz_cross_section_points(myrank,nglob_target,x2,y2,z2, &
                                        lat0,dlat,nlat,lon0,dlon,nlon,depth0,ddepth,ndepth, &
                                        TOPOGRAPHY,ELLIPTICITY)

  else
    ! vertical cross-section
    ! user output
    if (myrank == 0) then
      print *,'  vertical cross-section'
      print *
    endif

    ! checks depth input
    if (depth_min > depth_max) then
      ! switch depths
      depth0 = depth_min
      depth_min = depth_max
      depth_max = depth0
    endif

    ! number of depths between min/max radius (by default, surface down to CMB)
    r_top = R_PLANET_KM - depth_min
    r_bottom = R_PLANET_KM - depth_max
    ndepth = int((r_top - r_bottom)/ddepth) + 1

    ! starting depth
    depth0 = depth_min

    ! determines number of sections along great-circle arc
    ! gets geocentric locations
    call latlon_2_geocentric_thetaphi_dble(v_lat1,v_lon1,theta1,phi1,ELLIPTICITY)
    call latlon_2_geocentric_thetaphi_dble(v_lat2,v_lon2,theta2,phi2,ELLIPTICITY)

    ! computes epicentral distance (in radians)
    call get_greatcircle_distance(theta1,phi1,theta2,phi2,epidist)

    ! distance in degrees
    epidist = epidist * RADIANS_TO_DEGREES

    if (epidist <= SMALLVAL) then
      stop 'Error great-circle points too close to each other'
    endif

    ! determines number of steps along great-circle arc
    if (NCHUNKS_VAL == 6) then
      ! full circle
      nsec = int(360.d0 / dincr)
    else
      nsec = int(epidist / dincr) + 1
    endif

    ! user output
    if (myrank == 0) then
      print *,'  number of nsec/ndepth         = ',nlat,'/',nlon
      print *,'  total number of target points = ',nglob_target
      print *
    endif

    ! total number of points
    nglob_target = ndepth * nsec

    ! allocates arrays for target points
    allocate(x2(nglob_target), &
             y2(nglob_target), &
             z2(nglob_target), &
             model_distance2(nglob_target), &
             model2(nglob_target),stat=ier)
    if (ier /= 0) stop 'Error allocating target model point arrays'
    x2(:) = 0.0_CUSTOM_REAL; y2(:) = 0.0_CUSTOM_REAL; z2(:) = 0.0_CUSTOM_REAL
    model_distance2(:) = real(HUGEVAL,kind=CUSTOM_REAL)
    model2(:) = 0.0_CUSTOM_REAL

    ! creates cross-section points
    call set_vertical_cross_section_points(myrank,nglob_target,x2,y2,z2, &
                                           v_lat1,v_lon1,v_lat2,v_lon2,dincr,depth0,ddepth,ndepth,nsec,epidist, &
                                           TOPOGRAPHY,ELLIPTICITY)

  endif

  ! initializes interpolated model values
  model2(:) = 0.0_CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    print *
    print *, 'loading source mesh ... '
  endif

  ! mesh arrays for a single slice
  allocate(xstore(NGLOB_CRUST_MANTLE), &
           ystore(NGLOB_CRUST_MANTLE), &
           zstore(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating locations'
  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating ibool'

  ! reads in model and locations of old, source mesh
  ! combines all slices for this whole chunk
  xstore(:) = 0.0_CUSTOM_REAL
  ystore(:) = 0.0_CUSTOM_REAL
  zstore(:) = 0.0_CUSTOM_REAL
  ibool(:,:,:,:) = 0

  ! reads in mesh for this slice
  ! source mesh locations
  write(solver_file,'(a,i6.6,a)') trim(dir_topo)//'proc',myrank,trim(reg)//'solver_data.bin'
  open(IIN,file=solver_file,status='old',form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening solver_data file: ',trim(solver_file)
    stop 'Error opening old solver_data.bin file'
  endif
  read(IIN) nspec
  read(IIN) nglob

  ! checks dimensions
  if (nspec /= NSPEC_CRUST_MANTLE .or. nglob /= NGLOB_CRUST_MANTLE) then
    print *,'Error dimension of old, source mesh: solver_data nspec/nglob = ',nspec,nglob
    stop 'Error new mesh dimensions'
  endif

  read(IIN) xstore(:)
  read(IIN) ystore(:)
  read(IIN) zstore(:)
  read(IIN) ibool(:,:,:,:)
  close(IIN)

  ! user output
  if (myrank == 0) then
    print *
    print *,'  source mesh chunk read successfully'
    print *
  endif
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *,'loading source model ... '
    print *
    print *,'  parameter: ',trim(fname)
    print *
  endif

  ! model files
  allocate(model(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating initial model'

  ! reads in old model files
  model(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reads in model slice
  ! debug user output
  !if (myrank == 0) print *, '  for parameter: ',trim(fname)

  ! opens model file
  write(m_file,'(a,i6.6,a)') trim(input_model_dir)//'proc',myrank,trim(reg)//trim(fname)//'.bin'
  open(IIN,file=trim(m_file),status='old',form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening model file: ',trim(m_file)
    stop 'Error opening old model file'
  endif
  read(IIN) model(:,:,:,:)
  close(IIN)

  ! user output
  if (myrank == 0) then
    print *
    print *,'  source model slice read successfully'
    print *
  endif
  call synchronize_all()

  ! statistics
  model_maxdiff = 0.0_CUSTOM_REAL

  ! builds search tree
  ! counts total number of points in this slice
  ! takes all internal GLL points ( 2 to NGLLX-1 ) to build kd-tree points
  kdtree_num_nodes = NSPEC_CRUST_MANTLE * (NGLLX-2) * (NGLLY-2) * (NGLLZ-2)

  ! debug
  !print *,'kdtree nodes: ',kdtree_num_nodes

  ! allocates tree arrays
  allocate(kdtree_nodes_location(3,kdtree_num_nodes),stat=ier)
  if (ier /= 0) stop 'Error allocating kdtree_nodes_location arrays'
  allocate(kdtree_nodes_index(kdtree_num_nodes),stat=ier)
  if (ier /= 0) stop 'Error allocating kdtree_nodes_index arrays'

  ! tree verbosity
  if (myrank == 0) call kdtree_set_verbose()

  ! prepares search arrays, each element takes its internal GLL points for tree search
  kdtree_nodes_index(:) = 0
  kdtree_nodes_location(:,:) = 0.0

  ! fills all local nodes into tree array
  inodes = 0
  do ispec = 1,NSPEC_CRUST_MANTLE
    ! sets up tree nodes
    ! all internal GLL points
    do k = 2,NGLLZ-1
      do j = 2,NGLLY-1
        do i = 2,NGLLX-1
          iglob = ibool(i,j,k,ispec)

          ! counts nodes
          inodes = inodes + 1
          if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

          ! adds node index ( index points to same ispec for all internal GLL points)
          kdtree_nodes_index(inodes) = ispec

          ! adds node location
          kdtree_nodes_location(1,inodes) = xstore(iglob)
          kdtree_nodes_location(2,inodes) = ystore(iglob)
          kdtree_nodes_location(3,inodes) = zstore(iglob)
        enddo
      enddo
    enddo

  enddo
  if (inodes /= kdtree_num_nodes ) stop 'Error index inodes does not match nnodes_local'
  call synchronize_all()

  ! creates kd-tree for searching
  ! serial way
  !do i = 0,NPROCTOT_VAL-1
  !  if (myrank == i) then
  !    print *,'kd-tree setup for process: ',myrank
  !    call kdtree_setup()
  !  endif
  !  call synchronize_all()
  !enddo
  ! parallel way
  call kdtree_setup()

  ! synchronizes all mpi-processes
  call synchronize_all()

  ! gets model values
  ! kdtree search
  call get_model_values_cross_section(nglob_target,x2,y2,z2,model2,model_distance2, &
                                      NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,ibool,xstore,ystore,zstore,model, &
                                      anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll, &
                                      typical_size,myrank,model_maxdiff)

  ! frees tree memory
  ! deletes tree arrays
  deallocate(kdtree_nodes_location)
  deallocate(kdtree_nodes_index)

  ! deletes search tree nodes
  call kdtree_delete()

  ! frees memory
  deallocate(xstore,ystore,zstore)
  deallocate(ibool)
  deallocate(model)

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    print *,'elapsed time in seconds = ',tCPU
    print *
  endif

  ! collects best points on main proc
  call collect_closest_point_values(myrank,NPROCTOT_VAL,nglob_target,model2,model_distance2)

  ! statistics
  call max_all_cr(model_maxdiff,val)
  if (myrank == 0) then
    ! model min/max values
    min_val = minval(model2(:),mask = model_distance2(:) < 2 * typical_size)
    max_val = maxval(model2(:),mask = model_distance2(:) < 2 * typical_size)
    ! user output
    print *,'search statistics:','  parameter ',trim(fname)
    print *,'  min/max values = ',min_val,max_val
    print *

    ! distance
    ! note: the cross-section is over the whole lat/lon range [-90,90/0,360]
    !       even if the mesh is only regional; thus having cross-section points with huge distance can happen.
    !       we thus apply the mask to get the distance for those with valid model points
    if (NCHUNKS_VAL /= 6) then
      max_val = maxval(model_distance2(:),mask = model_distance2(:) < 2 * typical_size)
    else
      max_val = maxval(model_distance2(:))
    endif
    ! user output
    print *,'  maximum distance to target point = ',max_val * R_PLANET_KM,'(km)'
    print *,'  maximum model value difference between closest GLL point = ',val
    print *
  endif
  call synchronize_all()

  ! main process only
  if (myrank == 0) then
    ! allocates arrays for statistics
    allocate(model_diff(nglob_target), &
             model_pert(nglob_target),stat=ier)
    if (ier /= 0) stop 'Error allocating statistics arrays'

    ! gets statistics values
    call get_cross_section_avg(nglob_target,x2,y2,z2,model2,model_distance2, &
                               typical_size,dlat,dlon,nlat,lat0,dincr,ddepth,section_type,ELLIPTICITY, &
                               model_diff,model_pert,m_avg_total,point_avg_total)

    ! user output
    print *,'writing out cross-section:'
    print *, '  for parameter: ',trim(fname)

    ! cross-section file name
    m_file = '/cross_section_'//trim(fname)//'.dat'
    write(filename,'(a,i6.6,a)') trim(output_dir) // trim(m_file)

    ! writes out cross section
    call write_cross_section(nglob_target,x2,y2,z2,model2,model_distance2, &
                             model_diff,model_pert,m_avg_total,point_avg_total, &
                             typical_size,depth0,section_type,filename,fname)

    ! frees temporary arrays
    deallocate(model_diff,model_pert)
  endif

  ! frees memory
  deallocate(x2,y2,z2,model2,model_distance2)

  ! synchronizes MPI processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    ! elapsed time since beginning of the simulation
    tCPU = wtime() - time_start
    print *,'elapsed time in seconds = ',tCPU
    print *
    print *, 'done successfully'
    print *
  endif

  ! exiting MPI processes
  call finalize_mpi()

  end program cross_section

!
!------------------------------------------------------------------------------
!


  subroutine set_horiz_cross_section_points(myrank,nglob_target,x2,y2,z2, &
                                            lat0,dlat,nlat,lon0,dlon,nlon,depth0,ddepth,ndepth, &
                                            TOPOGRAPHY,ELLIPTICITY)

! creates point locations of horizontal cross-section points

  use constants, only: CUSTOM_REAL,NR_DENSITY,R_UNIT_SPHERE
  use shared_parameters, only: R_PLANET,NX_BATHY,NY_BATHY

  ! for testing
  use constants, only: PI_OVER_TWO

  implicit none

  integer,intent(in) :: myrank,nglob_target

  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(out) :: x2, y2, z2

  ! latitude-longitude grid
  integer,intent(in) :: nlat,nlon,ndepth
  double precision,intent(in) :: lat0,dlat
  double precision,intent(in) :: lon0,dlon
  double precision,intent(in) :: depth0,ddepth

  logical,intent(in) :: TOPOGRAPHY,ELLIPTICITY

  ! local parameters
  ! topography/bathymetry & oceans
  ! use integer array to store values
  integer, dimension(:,:),allocatable :: ibathy_topo

  ! for ellipticity
  integer :: nspl
  double precision,dimension(NR_DENSITY) :: rspl,ellipicity_spline,ellipicity_spline2

  double precision :: r0
  double precision :: elevation

  integer :: iglob,ier
  integer :: ilat,ilon,idep

  double precision :: theta,phi
  double precision :: lat,lon,depth
  double precision :: x_target,y_target,z_target

  ! debugging
  logical, parameter :: DO_TEST = .false.

  ! user output
  if (myrank == 0) then
    print *,'  depth: ',depth0,'km'
    print *,'  lat/lon increments: ',sngl(dlat),'/',sngl(dlon),'degrees'
    print *
    print *,'  total number of points = ',nglob_target
    print *,'  total size for point arrays = ',5 * dble(nglob_target) * dble(CUSTOM_REAL) / 1024. / 1024.,'MB'
    print *
    if (TOPOGRAPHY) print *,'  considering topography'
    if (ELLIPTICITY) print *,'  considering ellipticity'
    print *
  endif

  ! make ellipticity
  if (ELLIPTICITY) then
    ! splines used for locating exact positions
    call make_ellipticity(nspl,rspl,ellipicity_spline,ellipicity_spline2)
  endif

  ! read topography and bathymetry file
  if (TOPOGRAPHY) then
    ! allocates topography array
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating ibathy_topo array')

    ! initializes
    ibathy_topo(:,:) = 0

    ! main reads file
    if (myrank == 0) then
      ! user output
      print *,'topography:'

      ! reads topo file
      call read_topo_bathy_file(ibathy_topo)

      print *,"  topography/bathymetry: min/max = ",minval(ibathy_topo),maxval(ibathy_topo)
    endif

    ! broadcast the information read on the main node to all the nodes
    call bcast_all_i(ibathy_topo,NX_BATHY*NY_BATHY)
  endif
  call synchronize_all()

  ! get values in elements for this new slice
  iglob = 0
  do idep = 1,ndepth
    ! depth (in km)
    depth = depth0 + dble(idep-1) * ddepth
    ! converts to meters
    depth = depth * 1000.d0

    ! loops over lat/lon points
    do ilat = 1,nlat
      do ilon = 1,nlon
        ! north/south pole: uses only first longitude to have single point
        ! checks if at south pole
        if (ilat == 1 .and. ilon > 1) cycle
        ! checks if at north pole
        if (ilat == nlat .and. ilon > 1) cycle

        ! lat/lon/depth
        lat = lat0 + dble(ilat-1) * dlat
        lon = lon0 + dble(ilon-1) * dlon

        !if (myrank == 0) print *,'depth/lat/lon = ',depth,lat,lon

        ! get colatitude/longitude
        ! limits longitude to [0.0,360.0]
        if (lon < 0.d0 ) lon = lon + 360.d0
        if (lon > 360.d0 ) lon = lon - 360.d0

        ! converts geographic latitude lat (degrees) to geocentric colatitude theta (radians)
        ! (reduce theta/phi to range [0,PI] and [0,2PI])
        call latlon_2_geocentric_thetaphi_dble(lat,lon,theta,phi,ELLIPTICITY)

        ! normalized receiver radius
        r0 = R_UNIT_SPHERE

        ! finds elevation
        if (TOPOGRAPHY) then
          ! gets elevation in meters
          call get_topo_bathy(lat,lon,elevation,ibathy_topo)
          ! adds to spherical radius
          r0 = r0 + elevation/R_PLANET
        endif

        ! ellipticity
        if (ELLIPTICITY) then
          ! adds ellipticity factor to radius
          call add_ellipticity_rtheta(r0,theta,nspl,rspl,ellipicity_spline,ellipicity_spline2)
        endif

        ! subtracts desired depth (in meters)
        r0 = r0 - depth/R_PLANET

        !if (myrank == 0) print *,'  elevation = ',elevation,'ellip = ',ELLIPTICITY,'radius = ',r0 * R_PLANET_KM

        ! compute the Cartesian position of the receiver
        x_target = r0 * sin(theta) * cos(phi)
        y_target = r0 * sin(theta) * sin(phi)
        z_target = r0 * cos(theta)

        ! adds point
        iglob = iglob + 1
        if (iglob > nglob_target) stop 'Error iglob exceeds total size'

        x2(iglob) = real(x_target,kind=CUSTOM_REAL)
        y2(iglob) = real(y_target,kind=CUSTOM_REAL)
        z2(iglob) = real(z_target,kind=CUSTOM_REAL)
      enddo
    enddo
  enddo
  ! checks point count
  if (iglob /= nglob_target) stop 'Error iglob count invalid'

  ! double check conversion from geographic to geocentric positions
  if (DO_TEST .and. (myrank == 0)) then
    ! test position
    lat = 68.d0        ! over Alaska
    lon = 200.d0
    depth = 0.d0

    ! initializes
    theta = 0.d0
    phi = 0.d0
    r0 = 0.d0
    elevation = 0.d0  ! for topo

    ! converts geographic latitude lat (degrees) to geocentric colatitude theta (radians)
    ! (reduce theta/phi to range [0,PI] and [0,2PI])
    call latlon_2_geocentric_thetaphi_dble(lat,lon,theta,phi,ELLIPTICITY)

    ! normalized receiver radius
    r0 = R_UNIT_SPHERE

    print *
    print *,'test: lat/lon/r   = ',lat,lon,r0
    print *,'test: theta/phi/r = ',theta,phi,r0

    ! finds elevation
    if (TOPOGRAPHY) then
      ! gets elevation in meters
      call get_topo_bathy(lat,lon,elevation,ibathy_topo)
      ! adds to spherical radius
      r0 = r0 + elevation/R_PLANET
    endif

    print *,'test: after topo : theta/phi/r = ',theta,phi,r0

    ! ellipticity
    if (ELLIPTICITY) then
      ! adds ellipticity factor to radius
      call add_ellipticity_rtheta(r0,theta,nspl,rspl,ellipicity_spline,ellipicity_spline2)
    endif

    print *,'test: after ellip: theta/phi/r = ',theta,phi,r0

    ! subtracts desired depth (in meters)
    r0 = r0 - depth/R_PLANET

    print *,'test: elevation = ',elevation,'ellip = ',ELLIPTICITY,'radius = ',r0 * R_PLANET / 1000.0,'(km)'

    ! compute the Cartesian position of the receiver
    x_target = r0 * sin(theta) * cos(phi)
    y_target = r0 * sin(theta) * sin(phi)
    z_target = r0 * cos(theta)

    print *,'test: x/y/z = ',x_target,y_target,z_target

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    ! for ellipticity, this is correcting latitude
    call xyz_2_rlatlon_dble(x_target,y_target,z_target,r0,lat,lon,ELLIPTICITY)

    print *,'test: back to geographic lat/lon/r = ',lat,lon,r0
    print *,'test: done'
    print *
  endif ! test

  ! frees topography array
  if (TOPOGRAPHY) deallocate(ibathy_topo)

  end subroutine set_horiz_cross_section_points



!
!------------------------------------------------------------------------------
!


  subroutine set_vertical_cross_section_points(myrank,nglob_target,x2,y2,z2, &
                                               lat1,lon1,lat2,lon2,dincr,depth0,ddepth,ndepth,nsec,epidist, &
                                               TOPOGRAPHY,ELLIPTICITY)

! creates point locations of horizontal cross-section points

  use constants, only: R_UNIT_SPHERE,CUSTOM_REAL,NR_DENSITY, &
    DEGREES_TO_RADIANS,RADIANS_TO_DEGREES, &
    PI_OVER_TWO,TWO_PI,SMALLVAL

  use shared_parameters, only: R_PLANET,NX_BATHY,NY_BATHY

  implicit none

  integer,intent(in) :: myrank,nglob_target

  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(out) :: x2, y2, z2

  ! latitude-longitude grid
  integer,intent(in) :: nsec,ndepth
  double precision,intent(in) :: lat1,lon1
  double precision,intent(in) :: lat2,lon2
  double precision,intent(in) :: dincr
  double precision,intent(in) :: depth0,ddepth
  double precision,intent(in) :: epidist

  logical,intent(in) :: TOPOGRAPHY,ELLIPTICITY

  ! local parameters
  ! topography/bathymetry & oceans
  ! use integer array to store values
  integer, dimension(:,:),allocatable :: ibathy_topo

  ! for ellipticity
  integer :: nspl
  double precision,dimension(NR_DENSITY) :: rspl,ellipicity_spline,ellipicity_spline2

  double precision :: r0
  double precision :: elevation

  integer :: iglob,ier
  integer :: isec,idep

  double precision :: theta,phi
  double precision :: theta1,phi1,theta2,phi2
  double precision :: thn,phn,thr_new,phr_new,ths_new,phs_new,thp,php

  double precision :: lat,lon,depth
  double precision :: x_target,y_target,z_target

  logical :: use_positive_lon

  ! user output
  if (myrank == 0) then
    print *,'  starting depth         : ',sngl(depth0),'km'
    print *,'  reaching maximum depth : ',sngl(depth0 + (ndepth-1) * ddepth),'km'
    print *,'  depth increment        : ',ddepth,'km'
    print *
    print *,'  start location lat/lon = ',sngl(lat1),sngl(lon1)
    print *,'    end location lat/lon = ',sngl(lat2),sngl(lon2)
    print *
    print *,'  epicentral distance    : ',sngl(epidist),'degrees'
    print *
    print *,'  great-circle increments: ',sngl(dincr),'degrees'
    print *,'  total number of great-circle increments = ',nsec
    print *,'  total path length of great-circle = ',sngl((nsec-1)*dincr),'degrees'
    print *
    print *,'  total number of points = ',nglob_target
    print *,'  total size for point arrays = ',5 * dble(nglob_target) * dble(CUSTOM_REAL) / 1024. / 1024.,'MB'
    print *
    if (TOPOGRAPHY) print *,'  considering topography'
    if (ELLIPTICITY) print *,'  considering ellipticity'
    print *
  endif

  ! make ellipticity
  if (ELLIPTICITY) then
    ! splines used for locating exact positions
    call make_ellipticity(nspl,rspl,ellipicity_spline,ellipicity_spline2)
  endif

  ! read topography and bathymetry file
  if (TOPOGRAPHY) then
    ! allocates topography array
    allocate(ibathy_topo(NX_BATHY,NY_BATHY),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating ibathy_topo array')

    ! initializes
    ibathy_topo(:,:) = 0

    ! main reads file
    if (myrank == 0) then
      ! user output
      print *,'topography:'

      ! reads topo file
      call read_topo_bathy_file(ibathy_topo)

      print *,"  topography/bathymetry: min/max = ",minval(ibathy_topo),maxval(ibathy_topo)
    endif

    ! broadcast the information read on the main node to all the nodes
    call bcast_all_i(ibathy_topo,NX_BATHY*NY_BATHY)
  endif
  call synchronize_all()

  ! determines if longitudes in range [-180,180] or [0,360]
  if (lon1 < 0.d0 .or. lon2 < 0.d0) then
    use_positive_lon = .false.
  else
    use_positive_lon = .true.
  endif

  ! sets up great-circle points
  ! converts to radians and uses co-latitude
  call latlon_2_geocentric_thetaphi_dble(lat1,lon1,theta1,phi1,ELLIPTICITY)
  call latlon_2_geocentric_thetaphi_dble(lat2,lon2,theta2,phi2,ELLIPTICITY)

  ! reduce lat/lon to corresponding range
  call set_thetaphi_range(theta1,phi1,use_positive_lon)
  call set_thetaphi_range(theta2,phi2,use_positive_lon)

  ! debug
  !if (myrank == 0) print *,'colat/lon: ',theta1*RADIANS_TO_DEGREES,phi1*RADIANS_TO_DEGREES, &
  !                          '/',theta2*RADIANS_TO_DEGREES,phi2*RADIANS_TO_DEGREES

  ! great circle
  ! we will rotate the source/receiver plane into the equatorial one
  ! to easily increment points along the zero-meridian, and rotate them back
  ! figure out the normal of source-receiver plane, and rotate them to the equator
  call tp2norm(theta1,phi1,theta2,phi2,thn,phn)
  call set_thetaphi_range(thn,phn,use_positive_lon)

  ! debug
  !if (myrank == 0) print *,'normal: ',thn*RADIANS_TO_DEGREES,phn*RADIANS_TO_DEGREES

  ! new source position
  call norm_rot(thn,phn,theta1,phi1,ths_new,phs_new)
  call set_thetaphi_range(ths_new,phs_new,use_positive_lon)

  ! new receiver position
  call norm_rot(thn,phn,theta2,phi2,thr_new,phr_new)
  call set_thetaphi_range(thr_new,phr_new,use_positive_lon)

  ! debug
  !if (myrank == 0) print *,'rotated colat,lon = ',ths_new*RADIANS_TO_DEGREES,phs_new*RADIANS_TO_DEGREES, &
  !                          '/',thr_new*RADIANS_TO_DEGREES,phr_new*RADIANS_TO_DEGREES

  ! checks colatitude
  if ((ths_new - PI_OVER_TWO) > SMALLVAL) then
    print *,'New lat of first point is not on equator but, colat/lon = ',ths_new*RADIANS_TO_DEGREES,phs_new*RADIANS_TO_DEGREES
    stop 'Error rotating great-circle points'
  endif
  ! check colatitude
  if ((thr_new - PI_OVER_TWO) > SMALLVAL) then
    print *,'New lat of second point is not on equator but, colat/lon = ',thr_new*RADIANS_TO_DEGREES,phr_new*RADIANS_TO_DEGREES
    stop 'Error rotating great-circle points'
  endif

  ! ordering
  if (use_positive_lon) then
    if (phs_new > phr_new) then
      print *,'New rotating points switched longitude orders:',phs_new*RADIANS_TO_DEGREES,phr_new*RADIANS_TO_DEGREES
      stop 'Error rotating points switched longitude orders in range [0,2PI]'
    endif
  else
    if (phs_new > phr_new + TWO_PI) then
      print *,'New rotating points switched longitude orders:',phs_new*RADIANS_TO_DEGREES,phr_new*RADIANS_TO_DEGREES
      stop 'Error rotating points switched longitude orders in range [-PI,PI]'
    endif
  endif

  ! get values in elements for this new slice
  iglob = 0
  do idep = 1,ndepth
    ! depth (in km)
    depth = depth0 + dble(idep - 1) * ddepth
    ! converts to meters
    depth = depth * 1000.d0

    ! loops over great-circle sections
    do isec = 1,nsec
      ! new point location on equator
      thp = PI_OVER_TWO
      php = phs_new + dble(isec - 1) * dincr * DEGREES_TO_RADIANS

      ! actual point location
      call norm_rot_back(thn,phn,thp,php,theta,phi)
      call set_thetaphi_range(theta,phi,use_positive_lon)

      ! lat/lon in degrees
      lat = 90.d0 - theta * RADIANS_TO_DEGREES
      lon = phi * RADIANS_TO_DEGREES

      ! debug
      !if (myrank == 0) print *,'depth/lat/lon = ',depth,lat,lon

      ! normalized receiver radius
      r0 = R_UNIT_SPHERE

      ! finds elevation
      if (TOPOGRAPHY) then
        ! gets elevation in meters
        call get_topo_bathy(lat,lon,elevation,ibathy_topo)
        ! adds to spherical radius
        r0 = r0 + elevation/R_PLANET
      endif

      ! ellipticity
      if (ELLIPTICITY) then
        ! adds ellipticity factor to radius
        call add_ellipticity_rtheta(r0,theta,nspl,rspl,ellipicity_spline,ellipicity_spline2)
      endif

      ! subtracts desired depth (in meters)
      r0 = r0 - depth/R_PLANET

      !if (myrank == 0) print *,'  elevation = ',elevation,'ellip = ',ELLIPTICITY,'radius = ',r0 * R_PLANET_KM

      ! compute the Cartesian position of the receiver
      x_target = r0 * sin(theta) * cos(phi)
      y_target = r0 * sin(theta) * sin(phi)
      z_target = r0 * cos(theta)

      ! adds point
      iglob = iglob + 1
      if (iglob > nglob_target) stop 'Error iglob exceeds total size'

      x2(iglob) = real(x_target,kind=CUSTOM_REAL)
      y2(iglob) = real(y_target,kind=CUSTOM_REAL)
      z2(iglob) = real(z_target,kind=CUSTOM_REAL)
    enddo
  enddo
  ! checks point count
  if (iglob /= nglob_target) stop 'Error iglob count invalid'

  ! frees topography array
  if (TOPOGRAPHY) deallocate(ibathy_topo)

  end subroutine set_vertical_cross_section_points

!
!------------------------------------------------------------------------------
!


  subroutine collect_closest_point_values(myrank,nproc,nglob_target,model2,model_distance2)

  use constants, only: CUSTOM_REAL

  implicit none

  integer,intent(in) :: myrank,nproc
  integer,intent(in) :: nglob_target
  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(inout) :: model2,model_distance2

  ! local parameters
  integer :: iglob,iproc,ier
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: buffer
  real(kind=CUSTOM_REAL) :: dist,val
  ! MPI communication tag
  integer, parameter :: itag = 11

  ! temporary array
  allocate(buffer(2,nglob_target),stat=ier)
  if (ier /= 0) stop 'Error allocating buffer'
  buffer(:,:) = 0.0_CUSTOM_REAL

  ! main gets point distances
  if (myrank == 0) then
    ! main process collects info
    do iproc = 1,nproc - 1
      ! gets buffer arrays from secondary procs
      call recv_cr(buffer, 2 * nglob_target, iproc, itag)

      ! checks if closer point found
      do iglob = 1,nglob_target
        ! extracts info
        dist = buffer(1,iglob)
        val = buffer(2,iglob)

        ! debug
        !print *,'collect from ',iproc,'iglob = ',iglob,'value = ',val,'dist = ',dist

        ! fills in best points
        if (dist < model_distance2(iglob)) then
          model_distance2(iglob) = dist
          model2(iglob) = val
        endif
      enddo
    enddo
  else
    ! fills buffer, format: (distance,model-value)
    do iglob = 1,nglob_target
      buffer(1,iglob) = model_distance2(iglob)
      buffer(2,iglob) = model2(iglob)
    enddo

    ! secondary process sends its (distance,value) buffer to main proc
    call send_cr(buffer, 2 * nglob_target, 0, itag)
  endif

  ! synchronizes all mpi-processes
  call synchronize_all()

  ! free temporary array
  deallocate(buffer)

  end subroutine collect_closest_point_values

!
!------------------------------------------------------------------------------
!


  subroutine get_model_values_cross_section(nglob_target,x2,y2,z2,model2,model_distance2, &
                                            nspec,nglob,ibool,xstore,ystore,zstore,model, &
                                            anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll, &
                                            typical_size,myrank,model_maxdiff)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,HUGEVAL,SMALLVAL
  use shared_parameters, only: R_PLANET_KM

  use kdtree_search, only: kdtree_find_nearest_neighbor,kdtree_nodes_location, &
    kdtree_count_nearest_n_neighbors,kdtree_get_nearest_n_neighbors, &
    kdtree_search_index,kdtree_search_num_nodes,kdtree_nodes_index

  implicit none

  ! target mesh
  integer,intent(in) :: nglob_target
  real(kind=CUSTOM_REAL),dimension(nglob_target),intent(in) :: x2,y2,z2
  real(kind=CUSTOM_REAL),dimension(nglob_target),intent(inout) :: model2
  real(kind=CUSTOM_REAL),dimension(nglob_target),intent(inout) :: model_distance2

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(nglob),intent(in) :: xstore,ystore,zstore
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: model

  ! indices of the control points of the surface element
  integer,intent(in) :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL),intent(inout) :: model_maxdiff

  ! local parameters
  integer :: iglob
  ! interpolated point location
  integer :: ispec_selected
  double precision :: xi,eta,gamma
  ! point location
  real(kind=CUSTOM_REAL) :: x_target,y_target,z_target

  ! nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: dist_min,dist
  double precision :: element_size
  integer :: iglob_min
  integer :: i_selected,j_selected,k_selected
  logical :: is_inside_element

  double precision :: r_search
  integer :: nsearch_points,ielem,i,ier

  integer,dimension(:),allocatable :: search_elements
  integer :: num_search_elem
  logical :: is_listed

  ! locations
  real(kind=CUSTOM_REAL) :: x_found,y_found,z_found
  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax

  integer :: ipoin,nslice_points
  integer,dimension(:),allocatable :: slice_points

  ! debug warning about large model value differences
  logical,parameter :: DO_WARNING = .false.

  ! user output
  if (myrank == 0) print *,'looping over target points ...'
  call synchronize_all()

  ! initializes closest distance
  model_distance2(:) = real(HUGEVAL,kind=CUSTOM_REAL)

  ! normalized search radius (around target)
  r_search = 1.d0 * typical_size

  ! gets cross-section points in this slice
  ! slice dimensions
  xmin = minval(xstore)
  xmax = maxval(xstore)
  ymin = minval(ystore)
  ymax = maxval(ystore)
  zmin = minval(zstore)
  zmax = maxval(zstore)

  ! temporary array
  allocate(slice_points(nglob_target),stat=ier)
  if (ier /= 0) stop 'Error allocating slice_points array'
  slice_points(:) = 0

  ! counts points in this slice
  nslice_points = 0
  do iglob = 1, nglob_target
    ! target point location
    x_target = x2(iglob)
    y_target = y2(iglob)
    z_target = z2(iglob)

    ! check with slice dimensions
    if (x_target < xmin - r_search .or. x_target > xmax + r_search) cycle
    if (y_target < ymin - r_search .or. y_target > ymax + r_search) cycle
    if (z_target < zmin - r_search .or. z_target > zmax + r_search) cycle

    ! adds point index to cross-section points in this slice
    nslice_points = nslice_points + 1
    slice_points(nslice_points) = iglob
  enddo

  ! user output
  if (myrank == 0) then
    print *,'  using search radius: ',sngl(r_search * R_PLANET_KM),'(km)'
    print *
    print *,'  points contained in main slice: ',nslice_points,' out of ',nglob_target
    print *
  endif

  ! loops over all cross-section points in this slice
  do ipoin = 1, nslice_points
    ! gets global target point index
    iglob = slice_points(ipoin)

    ! user output
    if (myrank == 0) then
      if (ipoin == 1 .or. mod(ipoin,int(0.1*nslice_points)) == 0 .or. ipoin == nslice_points) then
        print *,'  point',ipoin,' out of ',nslice_points
      endif
    endif

    ! target point location
    x_target = x2(iglob)
    y_target = y2(iglob)
    z_target = z2(iglob)

    ! kdtree search for each single GLL point
    xyz_target(1) = x_target
    xyz_target(2) = y_target
    xyz_target(3) = z_target

    ! gets number of tree points within search radius
    ! (within search sphere)
    call kdtree_count_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

    ! debug
    !if (myrank == 0) print *,'search radius: ',r_search*R_PLANET_KM,'has number of neighbors = ',nsearch_points

    ! checks if points near target locations have been found, if not then it is outside of mesh
    if (nsearch_points == 0) cycle

    ! first, best guess element considered
    !
    ! finds closest point in source chunk
    call kdtree_find_nearest_neighbor(xyz_target,iglob_min,dist_min)

    ! selected closest element associated with tree node
    ispec_selected = iglob_min

    ! checks if point was found properly
    call check_point_result()

    ! debug
    !if (myrank == 0 .and. iglob < 100) &
    !  print *,'dist_min kdtree: ',dist_min * R_PLANET_KM,'(km)',typical_size * R_PLANET_KM

    ! restores original target point location for locating/interpolating
    ! (re-assignment might help for compiler optimizations with better locality)
    !x_target = x2(iglob)
    !y_target = y2(iglob)
    !z_target = z2(iglob)

    ! gets interpolated position within selected element
    call locate_single_point(x_target,y_target,z_target, &
                             xi,eta,gamma, &
                             ispec_selected,nspec,nglob, &
                             ibool,xstore,ystore,zstore, &
                             anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size, &
                             i_selected,j_selected,k_selected,dist_min, &
                             is_inside_element,element_size)

    ! checks point location
    call check_location()

    ! checks if point outside of mesh slice
    if (.not. is_inside_element .and. dist_min > 2.d0 * r_search) cycle

    ! checks if point is inside element
    if (is_inside_element .and. dist_min < SMALLVAL) then
      ! sets new interpolated model value
      call set_interpolated_value(dist_min,nglob_target,model2,model_distance2, &
                                  xi,eta,gamma,ispec_selected,i_selected,j_selected,k_selected,iglob, &
                                  nspec,model,xigll,yigll,zigll, &
                                  model_maxdiff,DO_WARNING,xyz_target,x_found,y_found,z_found,myrank)

      ! debug
      !if (myrank == 0) print *,'search first guess = ',dist_min*R_PLANET_KM,'point',ipoin

      ! no more searching needed, switch to next point
      cycle
    endif

    ! search continues with more element, looks within all closest elements
    !
    ! allocates search index
    allocate(kdtree_search_index(nsearch_points),stat=ier)
    if (ier /= 0) stop 'Error allocating array kdtree_search_index'
    kdtree_search_index(:) = 0

    ! gets closest n points in target chunk
    ! (within search sphere)
    kdtree_search_num_nodes = nsearch_points
    call kdtree_get_nearest_n_neighbors(xyz_target,r_search,nsearch_points)

    ! debug
    !if (myrank == 0) print *,'search nodes: ',kdtree_search_index(:)

    ! allocates search elements
    allocate(search_elements(nsearch_points),stat=ier)
    if (ier /= 0) stop 'Error allocating array search_elements array'
    search_elements(:) = 0

    ! gets all elements holding target point
    num_search_elem = 0
    do i = 1,nsearch_points
      ! gets search point/element index
      iglob_min = kdtree_search_index(i)
      ispec_selected = kdtree_nodes_index(iglob_min)

      ! debug
      !if (myrank == 0) print *,'search iglob,ispec = ',iglob_min,ispec_selected

      ! checks global point and element index
      if (iglob_min < 1 .or. iglob_min > nglob) stop 'Error global index is invalid'
      if (ispec_selected < 1 .or. ispec_selected > nspec) stop 'Error element index is invalid'

      ! checks if already listed
      is_listed = .false.
      do ielem = 1,num_search_elem
        if (search_elements(ielem) == ispec_selected) then
          is_listed = .true.
          exit
        endif
      enddo

      ! adds as new search element
      if (.not. is_listed) then
        num_search_elem = num_search_elem + 1
        search_elements(num_search_elem) = ispec_selected
      endif
    enddo

    ! checks all elements close to this target point
    nsearch_points = 0
    do ielem = 1,num_search_elem
      ! gets search element
      ispec_selected = search_elements(ielem)

      ! gets interpolated position within selected element
      call locate_single_point(x_target,y_target,z_target, &
                               xi,eta,gamma, &
                               ispec_selected,nspec,nglob, &
                               ibool,xstore,ystore,zstore, &
                               anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size, &
                               i_selected,j_selected,k_selected,dist_min, &
                               is_inside_element,element_size)

      ! checks if points is in element
      if (is_inside_element) then
        ! checks point location
        call check_location()

        ! debug
        !if (myrank == 0) print *,'search element: ',ispec_selected,'dist:',dist_min*R_PLANET_KM,'inside:',is_inside_element

        ! sets new interpolated model value
        call set_interpolated_value(dist_min,nglob_target,model2,model_distance2, &
                                    xi,eta,gamma,ispec_selected,i_selected,j_selected,k_selected,iglob, &
                                    nspec,model,xigll,yigll,zigll, &
                                    model_maxdiff,DO_WARNING,xyz_target,x_found,y_found,z_found,myrank)

        ! counter adds element to new search list
        nsearch_points = nsearch_points + 1
      endif

    enddo

    ! debug
    !if (myrank == 0) print *,'search elements = ',num_search_elem,' contained = ',nsearch_points,'radius:',r_search*R_PLANET_KM

    ! frees search arrays
    deallocate(kdtree_search_index)
    deallocate(search_elements)

  enddo

  ! free temporary array
  deallocate(slice_points)

  ! done looping over target locations
  if (myrank == 0) print *
  call synchronize_all()

  contains

    !-------------------------------------------

    subroutine check_point_result()

    implicit none

    ! checks valid iglob
    if (iglob_min < 1 .or. iglob_min > nspec) then
      print *,'Error iglob_min :',iglob_min
      print *,'nspec / nproc :',nspec
      stop 'Error invalid iglob_min index'
    endif

    ! checks valid ispec
    if (ispec_selected < 1 .or. ispec_selected > nspec) then
      print *,'Error rank:',myrank,'invalid selected ispec ',ispec_selected
      print *,'iglob_min:',iglob_min,'nspec:',nspec
      print *,'target location:',xyz_target(:)
      print *,'dist_min: ',dist_min * R_PLANET_KM,'(km)'
      stop 'Error specifying closest ispec element'
    endif

    ! checks minimum distance within a typical element size
    if (DO_WARNING) then
      if (dist_min > 2 * typical_size) then
        print *,'Warning: rank ',myrank,' - large dist_min: ',dist_min * R_PLANET_KM,'(km)', &
               'element size:',typical_size * R_PLANET_KM
        print *,'element selected ispec:',ispec_selected,'iglob_min:',iglob_min
        print *,'typical element size:',typical_size * 0.5 * R_PLANET_KM
        print *,'target location:',xyz_target(:)
        print *,'target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_PLANET_KM,'(km)'
        print *,'found location :',kdtree_nodes_location(:,iglob_min)
        print *,'found radius   :',sqrt(kdtree_nodes_location(1,iglob_min)**2 &
                                     + kdtree_nodes_location(2,iglob_min)**2 &
                                     + kdtree_nodes_location(3,iglob_min)**2 ) * R_PLANET_KM,'(km)'
        !debug
        !stop 'Error dist_min too large in check_point_result() routine'
      endif
    endif

    end subroutine check_point_result

    !-------------------------------------------

    subroutine check_location()

    implicit none

    ! checks closest GLL point
    iglob_min = ibool(i_selected,j_selected,k_selected,ispec_selected)
    x_found = xstore(iglob_min)
    y_found = ystore(iglob_min)
    z_found = zstore(iglob_min)

    ! checks distance
    if (DO_WARNING) then
      ! distance to closest GLL (still normalized)
      dist = sqrt((x_found-x_target)**2 + (y_found-y_target)**2 + (z_found-z_target)**2)
      if (dist > 2 * typical_size) then
        print *,'Warning: rank ',myrank,' - large GLL distance: ',dist * R_PLANET_KM,'(km)', &
               'element size:',element_size*R_PLANET_KM,'typical size:',typical_size * R_PLANET_KM
        print *,'target location:',xyz_target(:)
        print *,'target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_PLANET_KM,'(km)'
        print *,'gll location   :',x_found,y_found,z_found
        print *,'gll radius     :',sqrt(x_found**2 + y_found**2 + z_found**2) * R_PLANET_KM,'(km)'
        print *,'minimum distance gll:',dist_min,'(km)'
        ! debug
        !stop 'Error GLL model value invalid'
      endif
      ! debug
      !if (myrank == 0 .and. iglob_min < 100) &
      !  print *,'dist_min GLL point: ',dist_min * R_PLANET_KM,'(km)',typical_size * R_PLANET_KM
    endif

    end subroutine check_location

  end subroutine get_model_values_cross_section


!
!------------------------------------------------------------------------------
!

  subroutine set_interpolated_value(dist_min,nglob_target,model2,model_distance2, &
                                    xi,eta,gamma,ispec_selected,i_selected,j_selected,k_selected,iglob, &
                                    nspec,model,xigll,yigll,zigll, &
                                    model_maxdiff,DO_WARNING,xyz_target,x_found,y_found,z_found,myrank)

  use constants, only: NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use shared_parameters, only: R_PLANET_KM

  implicit none

  double precision,intent(in) :: dist_min
  ! new mesh
  integer,intent(in) :: nglob_target
  real(kind=CUSTOM_REAL),dimension(nglob_target),intent(inout) :: model2
  real(kind=CUSTOM_REAL),dimension(nglob_target),intent(inout) :: model_distance2

  ! point location
  double precision,intent(in) :: xi,eta,gamma
  integer,intent(in) :: ispec_selected
  integer,intent(in) :: i_selected,j_selected,k_selected
  integer,intent(in) :: iglob

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: model

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  real(kind=CUSTOM_REAL),intent(inout) :: model_maxdiff
  logical, intent(in) :: DO_WARNING

  ! nodes search
  double precision,dimension(3),intent(in) :: xyz_target
  real(kind=CUSTOM_REAL),intent(in) :: x_found,y_found,z_found
  integer, intent(in) :: myrank

  ! local parameters
  real(kind=CUSTOM_REAL) :: val,val_initial

  ! checks if new point distance is better than stored one
  if (dist_min < model_distance2(iglob)) then

    ! sets new minimum distance
    model_distance2(iglob) = real(dist_min,kind=CUSTOM_REAL)

    ! interpolate model values
    call interpolate_element_value(xi,eta,gamma,ispec_selected, &
                                   nspec,model(:,:,:,ispec_selected), &
                                   val,xigll,yigll,zigll)

    ! sets new model value
    model2(iglob) = val

    ! checks difference of interpolated value with model value of closest GLL point
    val_initial = model(i_selected,j_selected,k_selected,ispec_selected)
    if (abs(val - val_initial ) > abs(model_maxdiff)) model_maxdiff = val - val_initial

    ! checks model difference
    if (DO_WARNING) then
      ! note: warns for top elements, probably due to crustal structure
      if (abs(val - val_initial ) > abs( 0.2 * val_initial )) then
        print *,'Warning: model ',' value:',val,'is very different from initial value ',val_initial
        print *,'  rank ',myrank,' - dist_min: ',dist_min * R_PLANET_KM,'(km)'
        print *,'  selected ispec:',ispec_selected
        print *,'  interpolation i,j,k :',i_selected,j_selected,k_selected
        print *,'  interpolation       :',xi,eta,gamma
        print *,'  target location:',xyz_target(:)
        print *,'  target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_PLANET_KM,'(km)'
        print *,'  GLL location   :',x_found,y_found,z_found
        print *,'  GLL radius     :',sqrt(x_found**2 + y_found**2 + z_found**2) * R_PLANET_KM,'(km)'
        print *,'  distance gll:',dist_min * R_PLANET_KM,'(km)'
        !stop 'Error model value invalid'
      endif
    endif

    ! debug
    !if (myrank == 0 .and. iglob < 100) &
    !  print *,'new model ',': value ',val,'initial ',val_initial,'diff ',(val - val_initial)/val_initial*100.0,'(%)'

  endif

  end subroutine set_interpolated_value

!
!------------------------------------------------------------------------------
!


  subroutine locate_single_point(x_target,y_target,z_target, &
                                 xi_target,eta_target,gamma_target, &
                                 ispec_selected,nspec,nglob, &
                                 ibool,xstore,ystore,zstore, &
                                 anchor_iax,anchor_iay,anchor_iaz, &
                                 xigll,yigll,zigll,typical_size, &
                                 i_selected,j_selected,k_selected,dist_min, &
                                 is_inside_element,element_size)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,HUGEVAL,NUM_ITER
  use shared_parameters, only: R_PLANET_KM

  implicit none

  ! point location
  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target,z_target

  integer,intent(in) :: ispec_selected

  ! interpolated point location
  double precision,intent(out) :: xi_target,eta_target,gamma_target

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec,nglob

  ! arrays containing coordinates of the points
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(nglob),intent(in) :: xstore,ystore,zstore

  ! indices of the control points of the surface element
  integer,intent(in) :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size

  ! best guess for i,j,k
  integer,intent(out) :: i_selected,j_selected,k_selected

  ! minimum distance to target location (normalized)
  double precision,intent(out) :: dist_min

  logical,intent(out) :: is_inside_element
  double precision,intent(out) :: element_size

  ! local parameters
  ! use integer array to store values
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: i,j,k,iglob !,ispec,rank
  double precision :: dist
  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma
  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)
  integer :: iter_loop
  integer :: ia
  !integer :: ier

  double precision :: x,y,z
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz
  double precision :: final_distance

  integer :: iglob1,iglob2

  ! debugging
  logical,parameter :: DO_WARNING = .false.

  ! xi/eta/gamma coordinate limit for point within reference element
  double precision,parameter :: COORD_REFERENCE_LIMIT = 1.001d0

  ! checks
  if (ispec_selected < 1 ) stop 'Error locating: specifying closest ispec element'

  ! set distance to huge initial value
  dist_min = HUGEVAL
  ix_initial_guess = 0
  iy_initial_guess = 0
  iz_initial_guess = 0

  ! finds closest interior GLL point
  ! interior points should be good enough for iterative point search below
  do k = 2,NGLLZ-1
    do j = 2,NGLLY-1
      do i = 2,NGLLX-1
        iglob = ibool(i,j,k,ispec_selected)

        dist = (x_target - xstore(iglob))*(x_target - xstore(iglob)) &
              +(y_target - ystore(iglob))*(y_target - ystore(iglob)) &
              +(z_target - zstore(iglob))*(z_target - zstore(iglob))

        ! keep this point if it is closer to the receiver
        if (dist < dist_min) then
          dist_min = dist
          ix_initial_guess = i
          iy_initial_guess = j
          iz_initial_guess = k
        endif
      enddo
    enddo
  enddo

  ! checks
  if (ix_initial_guess == 0 .or. iy_initial_guess == 0 .or. iz_initial_guess == 0 ) &
    stop 'Error locating: no initial guess'

  ! initial minimum distance (still normalized)
  dist_min = sqrt(dist_min)

  ! debug
  !print *,'dist min = ',sngl(dist_min * R_PLANET_KM),'(km)'

  ! find the best (xi,eta)
  ! use initial guess in xi, eta and gamma from closest point found
  xi = xigll(ix_initial_guess)
  eta = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  ! iterate to solve the non linear system to find the exact position within the element

  ! initializes
  final_distance = HUGEVAL

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD
    iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec_selected)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))
  enddo

  dxi = 0.d0
  deta = 0.d0
  dgamma = 0.d0

  ! loop leads to invalid jacobian... probably some GLL points are too far outside of the selected element
  do iter_loop = 1,2*NUM_ITER

    ! recompute jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

    ! compute distance to target location
    dx = - (x - x_target)
    dy = - (y - y_target)
    dz = - (z - z_target)

    ! compute increments
    dxi  = xix*dx + xiy*dy + xiz*dz
    deta = etax*dx + etay*dy + etaz*dz
    dgamma = gammax*dx + gammay*dy + gammaz*dz

    ! decreases step length if step is large
    if ((dxi*dxi + deta*deta + dgamma*dgamma) > 1.0d0) then
      dxi = dxi * 0.33333333333d0
      deta = deta * 0.33333333333d0
      dgamma = dgamma * 0.33333333333d0
    endif
    ! alternative: impose limit on increments (seems to result in slightly less accurate locations)
    !if (abs(dxi) > 0.1d0 ) dxi = sign(1.0d0,dxi)*0.1d0
    !if (abs(deta) > 0.1d0 ) deta = sign(1.0d0,deta)*0.1d0
    !if (abs(dgamma) > 0.1d0 ) dgamma = sign(1.0d0,dgamma)*0.1d0

    ! update values
    xi = xi + dxi
    eta = eta + deta
    gamma = gamma + dgamma

    ! impose that we stay in that element
    ! (useful if user gives a receiver outside the mesh for instance)
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    ! can be useful for convergence of iterative scheme with distorted elements
    ! point leaves element, stay to closest guess
    if (xi > 1.20d0 .or. xi < -1.20d0 .or. eta > 1.20d0 .or. eta < -1.20d0 .or. gamma > 1.20d0 .or. gamma < -1.20d0) then
      ! uses previous guess
      xi = xi - dxi
      eta = eta - deta
      gamma = gamma - dgamma
      ! exits loop
      exit
    endif
  enddo

  ! compute final coordinates of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  ! found interpolated position
  ! compute final distance between asked and found (still normalized)
  final_distance = sqrt((x_target-x)*(x_target-x) + (y_target-y)*(y_target-y) + (z_target-z)*(z_target-z))

  ! debug
  !if (final_distance > 5.0 ) &
  !  print *,'final distance = ',sngl(final_distance * R_PLANET_KM),'(km)',dist_min * R_PLANET_KM,'xi/eta/gamma = ',xi,eta,gamma

  ! checks if location improved
  if (final_distance < dist_min ) then
    ! closest distance improved
    dist_min = final_distance
  else
    ! uses initial guess
    xi = xigll(ix_initial_guess)
    eta = yigll(iy_initial_guess)
    gamma = zigll(iz_initial_guess)
    final_distance = dist_min

    ! recomputes location x,y,z
    if (DO_WARNING) then
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
    endif
  endif

  ! checks valid distance
  if (final_distance == HUGEVAL) stop 'Error locating location'

  ! return xi,eta and gamma of point found
  xi_target = xi
  eta_target = eta
  gamma_target = gamma

  i_selected = ix_initial_guess
  j_selected = iy_initial_guess
  k_selected = iz_initial_guess

  ! element size (across diagonal)
  iglob1 = ibool(1,1,1,ispec_selected)
  iglob2 = ibool(NGLLX,NGLLY,NGLLZ,ispec_selected)
  element_size = sqrt((xstore(iglob1) - xstore(iglob2))**2 &
                    + (ystore(iglob1) - ystore(iglob2))**2 &
                    + (zstore(iglob1) - zstore(iglob2))**2)

  ! flag to indicate if point found lies inside or outside of element
  is_inside_element = .false.

  ! checks if interpolation weights small enough for point within element
  if (abs(xi_target) < COORD_REFERENCE_LIMIT &
      .and. abs(eta_target) < COORD_REFERENCE_LIMIT &
      .and. abs(gamma_target) < COORD_REFERENCE_LIMIT) then
    ! checks if distance close enough
    if (dist_min < 0.1 * element_size) then
      is_inside_element = .true.
      ! debug
      !print *,'xi/eta/gamma:',xi_target,eta_target,gamma_target,'dist',dist_min * R_PLANET_KM,'elemsize',element_size*R_PLANET_KM
    endif
  endif

  ! add warning if estimate is poor
  ! (usually means receiver outside the mesh given by the user)
  if (DO_WARNING) then
    if (dist_min > element_size) then
      print *, '*****************************************************************'
      print *, '***** WARNING: location estimate is poor                    *****'
      print *, '*****************************************************************'
      print *, 'closest estimate found: ',sngl(dist_min * R_PLANET_KM),'km away, not within:',typical_size * R_PLANET_KM
      print *, ' in element ',ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess
      print *, ' at xi,eta,gamma coordinates = ',xi,eta,gamma
      print *, ' at radius ',sqrt(x**2 + y**2 + z**2) * R_PLANET_KM,'(km)'
      print *, ' initial distance :',dist_min * R_PLANET_KM,'(km)'
    else
      print *,'point inside: ',is_inside_element, &
              'distances = ',dist_min * R_PLANET_KM,element_size * R_PLANET_KM,typical_size*R_PLANET_KM
      print *, ' at xi,eta,gamma coordinates = ',xi,eta,gamma
    endif

  endif

  end subroutine locate_single_point



!
!------------------------------------------------------------------------------
!

  subroutine write_cross_section(nglob_target,x2,y2,z2,model2,model_distance2, &
                                 model_diff,model_pert,m_avg_total,point_avg_total, &
                                 typical_size,depth0,section_type,filename,fname)

  use constants, only: CUSTOM_REAL,IOUT,MAX_STRING_LEN
  use shared_parameters, only: R_PLANET_KM

  implicit none

  integer,intent(in) :: nglob_target
  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(inout) :: x2, y2, z2
  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(in) :: model2,model_distance2
  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(in) :: model_diff,model_pert

  double precision,intent(in) :: m_avg_total,point_avg_total

  double precision,intent(in) :: typical_size
  double precision,intent(in) :: depth0
  integer,intent(in) :: section_type

  character(len=MAX_STRING_LEN),intent(in) :: filename
  character(len=16),intent(in) :: fname

  ! local parameters
  double precision :: r,lat,lon
  double precision :: val,pert,diff
  integer :: iglob,ier
  real(kind=CUSTOM_REAL) :: dist_min,distance_limit

  ! flag to chose to plot either all targets points (even if location is outside of mesh), or only close ones
  logical, parameter :: PLOT_ALL_POINTS = .false.

  ! minimum distance allowed
  distance_limit = real(1.d0 * typical_size,kind=CUSTOM_REAL)

  ! opens file
  open(IOUT,file=trim(filename),status='unknown', action='write',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(filename)
    stop 'Error opening output model file'
  endif

  ! note: to avoid line breaks after 80 characters, we specify a format in the write-statement:
  !         write(IOUT,'(a)') ..
  !       instead of
  !         write(IOUT,*) ..
  !       this is because compilers like ifort will add line breaks in the free formatted (list-directed) version,
  !       and then the files becomes unreadable, e.g., for a plotting script.

  ! header info
  write(IOUT,'(a)') '# cross-section informations:'
  if (section_type == 0) then
    write(IOUT,'(a)') '#   horizontal cross-section'
    write(IOUT,'(a)') '#'
    write(IOUT,'(a,f14.6)') '#   depth = ',depth0
    write(IOUT,'(a)') '#'
  else
    write(IOUT,'(a)') '#   vertical cross-section'
    write(IOUT,'(a)') '#'
  endif
  write(IOUT,'(a)') '#   parameter                           = ' // trim(fname)
  write(IOUT,'(a,f14.6)') '#   cross-section average value (m_avg) = ',m_avg_total
  write(IOUT,'(a,f14.6)') '#   point average value                 = ',point_avg_total
  write(IOUT,'(a)') '#'
  ! format: #lon #lat #parameter-value #actual-radius #minimum-distance(km)
  write(IOUT,'(a)') '#lon(degrees)   #lat(degrees)    #radius(km)   #parameter-value  #perturbation ln(m/m_avg)  ' // &
                    '#difference (m - m_avg)  #closest-point-distance(km)'

  ! writes out data points
  do iglob = 1,nglob_target
    ! gets r/lat/lon from converted arrays
    lat = x2(iglob)
    lon = y2(iglob)
    r = z2(iglob)

    ! distance to closest mesh point
    dist_min = model_distance2(iglob)

    ! plots points, or only points within close distance to mesh points will be considered
    if (PLOT_ALL_POINTS .or. dist_min < distance_limit) then

      ! model value
      val = model2(iglob)

      ! relative perturbations
      pert = model_pert(iglob)

      ! difference
      diff = model_diff(iglob)

      ! outputs cross-section points to file
      write(IOUT,'(7f14.6)') sngl(lon),sngl(lat),sngl(r*R_PLANET_KM),sngl(val), &
                             sngl(pert),sngl(diff),sngl(dist_min*R_PLANET_KM)
    endif
    !debug
    !print *,'lat/lon/r/param = ',sngl(lat),sngl(lon),sngl(r*R_PLANET_KM),model2(iglob),model_distance2(iglob)*R_PLANET_KM
  enddo

  ! closes file
  close(IOUT)

  ! user output
  print *
  print *, '  written cross-section file: ',trim(filename)
  print *

  end subroutine write_cross_section


!
!------------------------------------------------------------------------------
!

  subroutine get_cross_section_avg(nglob_target,x2,y2,z2,model2,model_distance2, &
                                   typical_size,dlat,dlon,nlat,lat0,dincr,ddepth,section_type,ELLIPTICITY, &
                                   model_diff,model_pert,m_avg_total,point_avg_total)

  use constants, only: CUSTOM_REAL,RADIANS_TO_DEGREES, &
    TINYVAL,HUGEVAL,PI,PI_OVER_TWO
  use shared_parameters, only: R_PLANET_KM

  implicit none

  integer,intent(in) :: nglob_target
  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(inout) :: x2, y2, z2
  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(in) :: model2,model_distance2

  double precision,intent(in) :: typical_size
  double precision,intent(in) :: dlat,dlon,lat0
  double precision,intent(in) :: dincr,ddepth
  integer,intent(in) :: nlat,section_type
  logical,intent(in) :: ELLIPTICITY

  real(kind=CUSTOM_REAL), dimension(nglob_target),intent(out) :: model_diff,model_pert
  double precision,intent(out) :: m_avg_total,point_avg_total

  ! local parameters
  double precision :: x,y,z
  double precision :: r,lat,lon

  integer :: iglob
  real(kind=CUSTOM_REAL) :: dist_min,distance_limit

  ! surface integral
  double precision :: area,total_spherical_area
  double precision :: model_integral
  double precision :: m_val,pert,diff
  double precision :: pert_min,pert_max
  integer :: ipoints_integral

  ! note: only main rank is computing this on collected array values
  !
  ! user output
  print *,'cross-section statistics:'

  ! average values
  m_avg_total = 0.d0
  point_avg_total = 0.d0

  ! surface integral
  model_integral = 0.d0
  total_spherical_area = 0.d0
  ipoints_integral = 0

  ! minimum distance allowed
  distance_limit = real(1.d0 * typical_size,kind=CUSTOM_REAL)

  ! loops over all cross-section points
  do iglob = 1,nglob_target
    ! converts x/y/z to radius/lat/lon
    x = x2(iglob)
    y = y2(iglob)
    z = z2(iglob)

    ! converts geocentric x/y/z to geographic lat/lon
    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

    ! stores lat/lon/r in x2,y2,z2 arrays
    ! such that we don't need to recompute these when writing out
    x2(iglob) = real(lat,kind=CUSTOM_REAL)
    y2(iglob) = real(lon,kind=CUSTOM_REAL)
    z2(iglob) = real(r,kind=CUSTOM_REAL)

    ! distance to closest mesh point
    dist_min = model_distance2(iglob)

    ! only points within close distance to mesh points will be considered for statistics
    if (dist_min < distance_limit) then
      ! associated area around point
      if (section_type == 0) then
        ! horizontal cross-section
        call get_regular_lat_lon_cell_area(r,lat,lon,dlat,dlon,nlat,lat0,area)
      else
        ! vertical cross-section
        call get_vertical_lat_lon_section_area(r,lat,lon,dincr,ddepth,area)
      endif

      ! model value
      m_val = model2(iglob)

      ! adds to model integral value
      model_integral = model_integral + area * m_val

      ! adds to total area
      total_spherical_area = total_spherical_area + area

      ! adds to point average
      point_avg_total = point_avg_total + m_val

      ! point counter
      ipoints_integral = ipoints_integral + 1
    endif
  enddo

  ! checks if some points got counted
  if (ipoints_integral == 0) then
    print *
    print *,'  Warning: no cross-section points within mesh volume'
    call exit_mpi(0,'No cross-section points found')
  endif
  if (total_spherical_area <= 0.d0) then
    print *
    print *,'  Warning: integrated surface area is invalid'
    print *
    call exit_mpi(0,'No valid cross-section area found')
  endif

  ! mean or average value over cross-section surface
  if (total_spherical_area /= 0.d0) then
    m_avg_total = model_integral / total_spherical_area
  endif
  if (ipoints_integral /= 0) then
    point_avg_total = point_avg_total / ipoints_integral
  endif

  ! initializes
  pert_min = HUGEVAL
  pert_max = - HUGEVAL
  model_pert(:) = 0._CUSTOM_REAL
  model_diff(:) = 0._CUSTOM_REAL

  ! determines difference to average and relative perturbation for all data points
  do iglob = 1,nglob_target
    ! distance to closest mesh point
    dist_min = model_distance2(iglob)

    ! model value
    m_val = model2(iglob)

    ! difference
    diff = m_val - m_avg_total
    model_diff(iglob) = real(diff,kind=CUSTOM_REAL)

    ! relative perturbations
    ! logarithmic perturbation: log( m_new) - log( m_avg) = log( m_new / m_avg )
    if (m_avg_total == 0.d0) then
      ! assumes that model values are already perturbations
      pert = m_val
    else
      if (m_avg_total > 0.d0 .and. m_val > 0.d0) then
        ! material values both positive
        pert = log( m_val / m_avg_total )
      else
        ! contains negative values
        pert = (m_val - m_avg_total) / abs(m_avg_total)
      endif
    endif
    model_pert(iglob) = real(pert,kind=CUSTOM_REAL)

    ! only points within close distance to mesh points will be considered for statistics
    if (dist_min < distance_limit) then
      ! statistics
      if (pert < pert_min) pert_min = pert
      if (pert > pert_max) pert_max = pert
    endif
  enddo

  ! user output
  print *,'  distance limit of close points       = ',sngl(distance_limit * R_PLANET_KM),'(km)'
  print *,'  number of close cross-section points = ',ipoints_integral
  print *
  print *,'  full Earth surface area    = ',sngl(4.0 * PI * R_PLANET_KM**2),'(km^2)'
  print *,'  cross-section surface area = ',sngl(total_spherical_area * R_PLANET_KM**2),'(km^2)'
  print *,'                             = ',sngl(total_spherical_area / (4.0 * PI) * 100.d0) ,'%'
  print *
  print *,'  cross-section average value = ',sngl(m_avg_total)
  print *,'          point average value = ',sngl(point_avg_total)
  print *
  print *,'  perturbation value min/max = ',sngl(pert_min),sngl(pert_max)
  print *

  end subroutine get_cross_section_avg

!
!------------------------------------------------------------------------------
!

  subroutine get_regular_lat_lon_cell_area(r,lat,lon,dlat,dlon,nlat,lat0,area)

! determines cell area around lat/lon location for a regular latitude-longitude grid point
! with dlat/dlon increments
!
! grid will have a single point at north/south pole

  use constants, only: DEGREES_TO_RADIANS

  implicit none

  ! point location (lat/lon in degrees)
  double precision,intent(in) :: r,lat,lon
  ! grid increments
  double precision,intent(in) :: dlat,dlon
  ! number of grid points along latitude
  integer,intent(in) :: nlat
  ! start latitude at -90. degrees
  double precision,intent(in) :: lat0

  double precision,intent(out) :: area

  ! local parameters
  double precision :: phi1,phi2
  double precision :: theta1,theta2
  double precision :: theta,phi
  double precision :: dtheta_half,dphi_half

  ! initializes
  area = 0.d0

  ! converts latitude (to radian)
  theta = lat * DEGREES_TO_RADIANS
  dtheta_half = 0.5d0 * dlat * DEGREES_TO_RADIANS

  ! converts longitude
  phi = lon * DEGREES_TO_RADIANS
  dphi_half = 0.5d0 * dlon * DEGREES_TO_RADIANS

  ! determines associated cell to lat/lon point
  if (lat > (lat0 + 0.5*dlat) .and. lat < (lat0 + (nlat - 1 - 0.5)*dlat)) then
    ! regular cell
    ! longitude cell borders
    phi1 = phi - dphi_half
    phi2 = phi + dphi_half

    ! latitude cell borders
    theta1 = theta - dtheta_half
    theta2 = theta + dtheta_half

    !debug
    !print *,'  phi1/phi2 = ',phi1/DEGREES_TO_RADIANS,phi2/DEGREES_TO_RADIANS, &
    !        'theta1/theta2 = ',theta1/DEGREES_TO_RADIANS,theta2/DEGREES_TO_RADIANS

    ! get cell area
    call get_area_spherical_cell(r,phi1,phi2,theta1,theta2,area)

  else if (lat < (lat0 + 0.5*dlat)) then
    ! south pole cap at -PI/2
    theta1 = theta + dtheta_half

    !debug
    !print *,'south pole :',theta1/DEGREES_TO_RADIANS

    ! get spherical cap area
    call get_area_spherical_cap(r,theta1,area)

  else if (lat > (lat0 + (nlat - 1 - 0.5)*dlat)) then
    ! north pole cap at PI/2
    theta1 = theta - dtheta_half

    !debug
    !print *,'north pole :',theta1/DEGREES_TO_RADIANS

    ! get spherical cap area
    call get_area_spherical_cap(r,theta1,area)
  else
    stop 'Error in selecting latitude region'
  endif

  end subroutine get_regular_lat_lon_cell_area

!
!------------------------------------------------------------------------------
!

  subroutine get_area_spherical_cell(r,phi1,phi2,theta1,theta2,area)

! returns the surface area of a grid cell at radius r included
! between phi1/phi2 longitudes [0,2*PI] and theta1/theta2 latitudes [-PI/2,PI/2] (given in radians)
!
! see: https://badc.nerc.ac.uk/help/coordinates/cell-surf-area.html

  implicit none

  double precision,intent(in) :: r
  double precision,intent(in) :: phi1,phi2
  double precision,intent(in) :: theta1,theta2
  double precision,intent(out) :: area

  ! local parameters
  double precision :: p1,p2,t1,t2

  ! longitudes
  p1 = phi1
  p2 = phi2
  ! checks order
  if (phi1 > phi2) then
    p1 = phi2
    p2 = phi1
  endif

  ! latitudes
  t1 = theta1
  t2 = theta2
  ! checks order
  if (theta1 > theta2) then
    t1 = theta2
    t2 = theta1
  endif

  ! cell area
  area = r**2 * (p2 - p1) * (sin(t2) - sin(t1))

  end subroutine get_area_spherical_cell

!
!------------------------------------------------------------------------------
!

  subroutine get_area_spherical_cap(r,theta,area)

! returns the surface area of a spherical cap between closest pole and given latitude,
! e.g. north pole at PI/2 and lat between [0,PI/2],
! of a sphere with radius r
!
! see: http://mathforum.org/library/drmath/view/63767.html

  use constants, only: PI

  implicit none

  double precision,intent(in) :: r
  double precision,intent(in) :: theta
  double precision,intent(out) :: area

  ! local parameters
  double precision :: t

  ! latitudes
  t = theta
  ! checks if in lower hemisphere
  if (theta < 0) then
    ! swaps hemisphere
    t = - theta
  endif

  ! spherical cap area
  area = r**2 * 2.d0 * PI * (1.d0 - sin(t))

  end subroutine get_area_spherical_cap

!
!------------------------------------------------------------------------------
!

  subroutine get_vertical_lat_lon_section_area(r,lat,lon,dincr,ddepth,area)

! determines area around lat/lon location for a regular, vertical latitude-longitude grid
! with dincr/ddepth increments

  use constants, only: R_UNIT_SPHERE,DEGREES_TO_RADIANS
  use shared_parameters, only: R_PLANET_KM

  implicit none

  ! point location (lat/lon in degrees)
  double precision,intent(in) :: r,lat,lon
  ! grid increments
  double precision,intent(in) :: dincr,ddepth

  double precision,intent(out) :: area

  ! local parameters
  double precision :: r_top,r_bot
  double precision :: theta,phi
  double precision :: dincr_rad,ddepth_half

  ! minimum radius for CMB (from AK135, see get_model_parameters.F90)
  double precision, parameter :: RCMB_LIMIT_KM   = 3479.5d0

  logical,parameter :: DEBUG = .false.

  ! initializes
  area = 0.d0

  ! converts latitude (to radian)
  theta = lat * DEGREES_TO_RADIANS
  ! converts longitude
  phi = lon * DEGREES_TO_RADIANS

  ! increment along great-circle path (in radians
  dincr_rad = dincr * DEGREES_TO_RADIANS

  ! normalized half vertical increment
  ddepth_half = 0.5d0 * ddepth / R_PLANET_KM

  ! radius bottom/top
  if (r + ddepth_half >= R_UNIT_SPHERE) then
    ! point at surface
    r_top = r
  else
    r_top = r + ddepth_half
  endif
  if (r - ddepth_half <= RCMB_LIMIT_KM) then
    ! point at crust/mantle bottom
    r_bot = r
  else
    r_bot = r - ddepth_half
  endif

  ! area of circle section around grid point (section between r_top and r_bot)
  area = 0.5 * (r_top**2 - r_bot**2) * dincr_rad

  !debug
  if (DEBUG) then
    print *,'  r/theta/phi = ',r*R_PLANET_KM,theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS,'area = ',area * R_PLANET_KM**2
  endif

  end subroutine get_vertical_lat_lon_section_area

!
!------------------------------------------------------------------------------
!

  subroutine set_thetaphi_range(theta,phi,use_positive_lon)

! returns co-latitude and longitude in range [0,PI]
! and [0,2PI] (use_positive_lon == true) / [-PI,PI] (use_positive_lon == false)

  use constants, only: PI,TWO_PI

  implicit none

  double precision,intent(inout) :: theta,phi
  logical, intent(in) :: use_positive_lon

  ! reduce to range [0,PI] and [0,2PI]
  call reduce(theta,phi)

  ! for range [-PI,PI]
  if (.not. use_positive_lon) then
    if (phi > PI) then
      phi = phi - TWO_PI
    endif
  endif

  end subroutine set_thetaphi_range

!
!------------------------------------------------------------------------------
!

  subroutine tp2norm(ths,phs,thr,phr,thn,phn)

! find out the normal of the plane constructed
!  by source (ths,phs), receiver (thr,phr), and the origin using
!  cross-product.

  implicit none

  double precision,intent(in) :: ths,phs,thr,phr
  double precision,intent(out) :: thn,phn

  ! local parameters
  double precision :: xs,ys,zs,xr,yr,zr
  double precision :: nx,ny,nz

  xs = 0.d0
  ys = 0.d0
  zs = 0.d0
  xr = 0.d0
  yr = 0.d0
  zr = 0.d0

  call tp2xyz(ths,phs,xs,ys,zs)
  call tp2xyz(thr,phr,xr,yr,zr)

  ! normal
  nx = ys*zr - zs*yr
  ny = -(xs*zr - zs*xr)
  nz = xs*yr - ys*xr

  call xyz2tp(nx,ny,nz,thn,phn)

  end subroutine tp2norm

!
!------------------------------------------------------------------------------
!

  subroutine tp2xyz(th,ph,x,y,z)

! convert (th,ph) to (x,y,z) on unit sphere

  implicit none

  double precision,intent(in) :: th,ph
  double precision,intent(out) :: x,y,z

  x = sin(th) * cos(ph)
  y = sin(th) * sin(ph)
  z = cos(th)

  end subroutine tp2xyz

!
!------------------------------------------------------------------------------
!

  subroutine xyz2tp(x,y,z,th,ph)

! convert x,y,z to a point on unit sphere

  implicit none

  double precision,intent(in) :: x,y,z
  double precision,intent(out) :: th,ph

  ph = atan2(y,x)
  th = atan2(sqrt(x*x+y*y),z)

  end subroutine xyz2tp

!
!------------------------------------------------------------------------------
!

  subroutine norm_rot(thn,phn,th,ph,th_new,ph_new)

! coordinates change from (th,ph) to (th_new,ph_new)
! according to a rotation that converts (thn,phn) to
! z axis

  implicit none

  double precision,intent(in) :: thn,phn,th,ph
  double precision,intent(out) :: th_new,ph_new

  ! local parameters
  double precision,dimension(3,3) :: rot
  double precision :: x,y,z
  double precision :: x_new,y_new,z_new

  ! sets up rotation matrix
  rot(1,1) = cos(thn)*cos(phn)
  rot(1,2) = cos(thn)*sin(phn)
  rot(1,3) = -sin(thn)
  rot(2,1) = -sin(phn)
  rot(2,2) = cos(phn)
  rot(2,3) = 0.d0
  rot(3,1) = sin(thn)*cos(phn)
  rot(3,2) = sin(thn)*sin(phn)
  rot(3,3) = cos(thn)

  call tp2xyz(th,ph,x,y,z)

  x_new = rot(1,1) * x + rot(1,2) * y + rot(1,3) * z
  y_new = rot(2,1) * x + rot(2,2) * y + rot(2,3) * z
  z_new = rot(3,1) * x + rot(3,2) * y + rot(3,3) * z

  call xyz2tp(x_new,y_new,z_new,th_new,ph_new)

  end subroutine norm_rot

!
!------------------------------------------------------------------------------
!

  subroutine norm_rot_back(thn,phn,th,ph,th_old,ph_old)

! coordinates change from (th,ph) to (th_old,ph_old)
! according to a rotation that converts (thn,phn) to z axis

  implicit none

  double precision,intent(in) :: thn,phn,th,ph
  double precision,intent(out) :: th_old,ph_old

  ! local parameters
  double precision,dimension(3,3) :: rot
  double precision :: x,y,z
  double precision :: x_old,y_old,z_old

  !rotation matrix
  rot(1,1) = cos(thn)*cos(phn)
  rot(1,2) = cos(thn)*sin(phn)
  rot(1,3) = -sin(thn)
  rot(2,1) = -sin(phn)
  rot(2,2) = cos(phn)
  rot(2,3) = 0.d0
  rot(3,1) = sin(thn)*cos(phn)
  rot(3,2) = sin(thn)*sin(phn)
  rot(3,3) = cos(thn)

  call tp2xyz(th,ph,x,y,z)

  x_old = rot(1,1) * x + rot(2,1) * y + rot(3,1) * z
  y_old = rot(1,2) * x + rot(2,2) * y + rot(3,2) * z
  z_old = rot(1,3) * x + rot(2,3) * y + rot(3,3) * z

  call xyz2tp(x_old,y_old,z_old,th_old,ph_old)

  end subroutine norm_rot_back

