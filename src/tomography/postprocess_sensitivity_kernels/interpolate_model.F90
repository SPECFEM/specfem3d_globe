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
! interpolates from one mesh resolution to another
!
! inputs:
!   - model parameter                       e.g. vsv, vsh   (-> proc***_reg1_vsv.bin)
!
!   - source model directory                 e.g. MODEL_M10/ (directory where proc***_reg1_vsv.bin lies)
!
! outputs:
!   - velocity model into output_model_dir/
!
!
! needs:
!   - topo files from first mesh            e.g. topo_1/proc000001_reg1_solver_data.bin
!
!   - topo files from target mesh           e.g. topo_2/proc000001_reg1_solver_data.bin
!
!
! usage: xinterpolate_model old-topo-dir/ old-model-dir/ new-topo-dir/ model-output-dir/ (midpoint-search) (nchunks_old) (iso/tiso)
!
! note on mid-point search option:
!  - mid-point-search == 1: looking for mid-points only is a good approach when changing number of processes (NPROC) only,
!                           while keeping the mesh resolution (NEX) the same
!                           (by default set to .true.)
!  - mid-point-search == 0: searching for each single GLL point is a good approach when changing resolution (NEX) of meshes;
!                           in general, interpolation suffers and might lead to differences at internal interfaces (e.g. 410)
!
!------------------------------------------------------------------------------


  program interpolate_model

  use constants, only: myrank

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
    GAUSSALPHA,GAUSSBETA, &
    SIZE_INTEGER,IIN,IOUT,MAX_STRING_LEN, &
    TWO_PI,R_UNIT_SPHERE,HUGEVAL_SNGL, &
    NGNOD,MIDX,MIDY,MIDZ, &
    IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80,IFLAG_670_220,IFLAG_MANTLE_NORMAL

  use shared_parameters, only: R_PLANET_KM,LOCAL_PATH

  use postprocess_par, only: &
    NCHUNKS_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,NPROCTOT_VAL,NEX_XI_VAL,NEX_ETA_VAL, &
    NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE

  use kdtree_search, only: kdtree_setup,kdtree_set_verbose,kdtree_delete,kdtree_find_nearest_neighbor, &
    kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use manager_adios
  use adios_helpers_mod
#endif

  implicit none

  !-------------------------------------------------------------------
  ! USER PARAMETERS

  ! model type
  ! isotropic model parameters (vp,vs,rho) or
  ! transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)
  ! transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)
  ! azimuthally anisotropic model parameters (vpv,vph,vsv,vsh,eta,rho,mu0,Gc_prime,Gs_prime)
  ! defaults: TI models
  logical :: USE_TRANSVERSE_ISOTROPY = .true.
  logical :: USE_AZIMUTHAL_ANISOTROPY = .false.

  ! shear attenuation
  logical :: USE_ATTENUATION_Q = .false.

  ! brute-force search for closest element
  ! by default set to .false., thus a tree search for initial guess element is used
  logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

  ! kd-tree setup
  ! uses all internal GLL points for search tree
  ! or only element mid-points (only one option must be true)
  logical,parameter :: TREE_INTERNAL_GLL_POINTS = .true.
  logical,parameter :: TREE_MID_POINTS = .false.

  ! special case for elements
  ! around 410-km discontinuity where internal topography distorts meshes
  logical,parameter :: DO_SEPARATION_410_650 = .true.
  ! around surface (due to moho-stretching)
  logical,parameter :: DO_SEPARATION_TOPO = .true.

  ! use closest point value in case of large differences
  logical,parameter :: USE_FALLBACK = .true.

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

  ! old, initial (source) mesh
  ! check constants defined for old, first mesh from file: 'values_from_mesher.h'
  ! NPROC_XI and NPROC_ETA
  integer :: nproc_eta_old,nproc_xi_old,nproctot_old
  integer :: nex_xi_old
  integer :: nchunks_old
  integer :: nspec_old,nglob_old
  ! combined slices, old mesh
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: x1,y1,z1
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:,:),allocatable :: model1
  integer,dimension(:,:,:,:,:),allocatable :: ibool1
  integer,dimension(:,:),allocatable :: idoubling1
  integer, dimension(:,:,:),allocatable :: addressing1

  ! new mesh
  integer :: nproc_eta_new,nproc_xi_new,nproctot_new
  integer :: nex_xi_new
  ! single slice, target mesh
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: x2, y2, z2
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: model2
  integer, dimension(:,:,:,:),allocatable :: ibool2
  integer, dimension(:),allocatable :: idoubling2
  integer, dimension(:,:,:),allocatable :: addressing2

  ! input arguments
  character(len=MAX_STRING_LEN) :: arg
  character(len=MAX_STRING_LEN) :: dir_topo1
  character(len=MAX_STRING_LEN) :: dir_topo2
  character(len=MAX_STRING_LEN) :: input_model_dir,output_model_dir

  ! kdtree search
  integer :: want_midpoint
  logical :: USE_MIDPOINT_SEARCH

  integer :: i,j,k,iglob,ispec,ier,iker
  integer :: nspec,nglob
  integer :: nspec_new,nglob_new
  integer :: rank,nproc_chunk1,nchunks_new
  integer :: ichunk,ichunk_selected
  integer :: iproc_eta,iproc_xi,iprocnum
  integer :: iproc_eta_selected,iproc_xi_selected
  integer :: ilayer
  integer :: total_nspec

  ! model
  integer :: nparams
  character(len=16),dimension(10) :: fname
  character(len=MAX_STRING_LEN) :: m_file
  character(len=MAX_STRING_LEN) :: solver_file

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer(kind=8) :: local_dim
  integer(kind=8) :: global_dim
  character(len=MAX_STRING_LEN) :: group_name
  integer(kind=8) :: group_size_inc
#endif

  ! MPI processes
  integer :: sizeprocs

  ! nodes search
  integer :: inodes
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: model_maxdiff
  real(kind=CUSTOM_REAL) :: val,val_all
  double precision :: sizeval
  logical :: use_single_process_per_chunk
  integer :: istart_xi,iend_xi,istart_eta,iend_eta
  integer :: istart_chunk,iend_chunk

  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: point_distance
  integer,dimension(4) :: loc_max
  logical :: is_updated,is_updated_all
  logical :: target_is_local_mesh

  ! timing
  double precision, external :: wtime
  double precision :: tstart,tCPU

  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! checks program arguments
  if (myrank == 0) then
    do i = 1,4
      call get_command_argument(i,arg)
      ! usage info
      if (len_trim(arg) == 0) then
        print *,' '
        print *,' Usage: xinterpolate_model old-topo-dir/ old-model-dir/ new-topo-dir/ model-output-dir/ (midpoint-search)' &
                // ' (nchunks) (model_p)'
        print *,' '
        print *,' with'
#ifdef USE_ADIOS_INSTEAD_OF_MESH
        print *,'   old-topo-dir/     - old mesh directory with topology files (e.g. solver_data.bp)'
        print *,'   old-model-dir/    - directoy which holds old model files (e.g. model_gll.bp)'
        print *,'   new-topo-dir/     - new mesh directory with topology files (e.g. solver_data.bp)'
#else
        print *,'   old-topo-dir/     - old mesh directory with topology files (e.g. proc***_solver_data.bin)'
        print *,'   old-model-dir/    - directoy which holds old model files (e.g. proc***_vpv.bin)'
        print *,'   new-topo-dir/     - new mesh directory with topology files (e.g. proc***_solver_data.bin)'
#endif
        print *,'   model-output-dir/ - output directory for interpolated model on new mesh'
        print *,'   (optional) midpoint-search = 0  - uses every single GLL point for search of closest element'
        print *,'                              = 1  - uses midpoints for search of closest element (default)'
        print *,'   (optional) nchunks - NCHUNKS (1,2,6) of old-model'
        print *,'   (optional) model_p - "iso" takes isotropic (vp,vs,rho) instead of transversely isotropic (vpv,vph,..) model,'
        print *,'                        "tiso_att" adds Q model (qmu) parameter, "iso_att" takes isotropic + Q model'
        print *,'                        "azi" takes azimuthally anisotropic (vpv, vph, ..., mu0, Gc_prime, Gs_prime) model,'
        print *,'                        "azi_att" adds Q model (qmu) parameter to azimuthally anisotropic model'
        print *,' '
        stop ' Reenter command line options'
      endif
    enddo
  endif
  call synchronize_all()

  ! initializes
  nchunks_old = 0
  nchunks_new = 0

  ! reads input arguments
  want_midpoint = 1
  do i = 1, 7
    call get_command_argument(i,arg)
    ! assignes values
    select case (i)
    case (1)
      dir_topo1 = trim(arg)
    case (2)
      input_model_dir = trim(arg)
    case (3)
      dir_topo2 = trim(arg)
    case (4)
      output_model_dir = trim(arg)
    end select
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      if (i == 5 .and. len_trim(arg) > 0) then
        read(arg(1:len_trim(arg)),*,iostat=ier) want_midpoint
        if (ier /= 0) then
          if (myrank == 0) print *,'Error reading in midpoint-search value, please check your arguments...'
          stop ' Reenter command line options'
        endif
      endif
    endif
    if (i == 6 .and. len_trim(arg) > 0) then
      read(arg(1:len_trim(arg)),*,iostat=ier) nchunks_old
    endif
    if (i == 7 .and. len_trim(arg) > 0) then
      if (trim(arg) == 'iso') then
        USE_TRANSVERSE_ISOTROPY = .false.
        USE_AZIMUTHAL_ANISOTROPY = .false.
        USE_ATTENUATION_Q = .false.
      else if (trim(arg) == 'iso_att') then
        USE_TRANSVERSE_ISOTROPY = .false.
        USE_AZIMUTHAL_ANISOTROPY = .false.
        USE_ATTENUATION_Q = .true.
      else if (trim(arg) == 'tiso_att') then
        USE_TRANSVERSE_ISOTROPY = .true.
        USE_AZIMUTHAL_ANISOTROPY = .false.
        USE_ATTENUATION_Q = .true.
      else if (trim(arg) == 'azi') then
        USE_TRANSVERSE_ISOTROPY = .false.
        USE_AZIMUTHAL_ANISOTROPY = .true.
        USE_ATTENUATION_Q = .false.
      else if (trim(arg) == 'azi_att') then
        USE_TRANSVERSE_ISOTROPY = .false.
        USE_AZIMUTHAL_ANISOTROPY = .true.
        USE_ATTENUATION_Q = .true.
      endif
    endif
  enddo

  ! kdtree search:
  ! searches closest element using mid-points only, rather than for every single GLL point
  ! note: looking for mid-points only is a good approach when changing number of processes (NPROC)
  !       while keeping the mesh resolution (NEX) the same;
  !       searching for each single GLL point is a good approach when changing resolution (NEX) of meshes;
  !       in general, interpolation suffers and might lead to differences at internal interfaces (e.g. 410);
  if (want_midpoint == 1) then
    USE_MIDPOINT_SEARCH = .true.
  else
    USE_MIDPOINT_SEARCH = .false.
  endif

  ! console output
  if (myrank == 0) then
    print *
    print *,'model interpolation:'
    print *,'  old topo  : ',trim(dir_topo1)
    print *,'  old model : ',trim(input_model_dir)
    print *,'  new topo  : ',trim(dir_topo2)
    print *,'  output dir: ',trim(output_model_dir)
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      print *,'  using midpoint search: ',want_midpoint
    else
      print *,'  using brute force search'
    endif
    if (nchunks_old /= 0) then
      print *,'  old topo NCHUNKS: ',nchunks_old
    endif
    print *
  endif

  ! reads mesh parameters
  if (myrank == 0) then
    ! reads mesh_parameters.bin file from input1dir/
    LOCAL_PATH = dir_topo2
    call read_mesh_parameters()
  endif
  ! broadcast parameters to all processes
  call bcast_mesh_parameters()

  ! user output
  if (myrank == 0) then
    print *,'mesh parameters (from new topo directory):'
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
      print *, "Please run: mpirun -np ",NPROCTOT_VAL," ./bin/xinterpolate_model .."
    endif
    call abort_mpi()
  endif

  ! safety check, assumes to have NEX_XI == NEX_ETA for now
  if (NEX_XI_VAL /= NEX_ETA_VAL) then
    print *,'Error: NEX_XI must be equal to NEX_ETA'
    stop 'Invalid NEX_XI/NEX_ETA values'
  endif

  ! initializes chunks
  ! new (target) mesh
  nchunks_new = NCHUNKS_VAL  ! from compilation
  nproctot_new = NPROCTOT_VAL
  nproc_xi_new = NPROC_XI_VAL
  nproc_eta_new = NPROC_ETA_VAL
  nex_xi_new = NEX_XI_VAL
  nspec_new = NSPEC_CRUST_MANTLE
  nglob_new = NGLOB_CRUST_MANTLE

  ! old (source) mesh
  if (nchunks_old == 0) then
    nchunks_old = nchunks_new  ! by default assumes same nchunks (e.g., global to global interpolation)
  endif
  nproctot_old = 0
  nproc_eta_old = 0
  nproc_xi_old = 0
  nex_xi_old = 0
  nspec_old = 0
  nglob_old = 0

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  call synchronize_all()
  ! initializes ADIOS
  if (myrank == 0) then
    print *, 'initializing ADIOS...'
    print *, ' '
  endif
  call initialize_adios()

  ! initializes i/o group
  call init_adios_group(myadios_group,"Interpolator")
#endif

  !  main process gets old, source mesh dimension
  if (myrank == 0) then

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    solver_file = get_adios_filename(trim(dir_topo1)//'/solver_data')

    ! user output
    print *, 'reading in ADIOS solver file: ',trim(solver_file)

    ! opens adios file
    call open_file_adios_read_only_rank(myadios_file,myadios_group,myrank,solver_file)

    ! reads in scalars
    print *,'ADIOS file only rank ',myrank,'reading scalars'
    call read_adios_scalar(myadios_file,myadios_group,myrank,"reg1/nspec",nspec_old)
    call read_adios_scalar(myadios_file,myadios_group,myrank,"reg1/nglob",nglob_old)
    print *,'  nspec/nglob = ',nspec_old,'/',nglob_old

    ! determines total number of processes
    call read_adios_scalar(myadios_file,myadios_group,myrank,"reg1/ibool/local_dim",local_dim)
    call read_adios_scalar(myadios_file,myadios_group,myrank,"reg1/ibool/global_dim",global_dim)
    print *,'  local_dim/global_dim = ',local_dim,'/',global_dim

    rank = 0
    ! note: the mod-function is only standard for integer, but an extension for integer(8).
    !       thus, instead of testing mod(global_dim,local_dim) == 0,
    !       we will just do what mod() would do and compute mod(a,b) = a - int(a/b) * b explicitly
    if ((global_dim - int(global_dim/local_dim) * local_dim) == 0) then
      ! sizeprocs
      rank = int(global_dim/local_dim)
    else
      print *,'Error invalid local_dim/global_dim ratio in file: ',trim(solver_file)
      stop 'Error adios array has invalid local_dim/global_dim ratio'
    endif
    ! closes file
    call close_file_adios_read_and_finalize_method_only_rank(myadios_file,myrank)

#else
    ! default binary
    ! user output
    print *, 'reading in binary solver files: ',trim(dir_topo1)//'/proc***_reg1_solver_data.bin'

    ! gets nspec/nglob
    ! opens file
    write(solver_file,'(a,i6.6,a)') trim(dir_topo1)//'/proc',myrank,'_reg1_'//'solver_data.bin'
    open(IIN,file=trim(solver_file),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(solver_file)
      stop 'Error opening old solver_data.bin file, please check arguments...'
    endif
    read(IIN) nspec_old
    read(IIN) nglob_old
    close(IIN)

    ! gets number of processes from old mesh (1 file per process)
    rank = 0
    do while (ier == 0)
      write(solver_file,'(a,i6.6,a)') trim(dir_topo1)//'/proc',rank,'_reg1_'//'solver_data.bin'
      open(IIN,file=trim(solver_file),status='old',action='read',form='unformatted',iostat=ier)
      if (ier == 0) then
        rank = rank + 1
        close(IIN)
      endif
    enddo
#endif

    ! assumes same number of chunks and nproc_eta = nproc_xi
    nproc_eta_old = int (sqrt ( dble(rank) / nchunks_old ))
    nproc_xi_old = int (sqrt ( dble(rank) / nchunks_old ))

    ! sets old nproc_xi (assumes equal nproc_xi/nproc_eta)
    !nproc_xi_old = nproc_eta_old

    ! total procs
    nproctot_old = nproc_xi_old * nproc_eta_old * nchunks_old
    ! try with global chunks if numbers don't match
    if (nproctot_old /= rank .and. nchunks_old /= 6) then
      ! check with 6 chunk mesh
      if (int (sqrt ( dble(rank) / 6 )) * int (sqrt ( dble(rank) / 6 )) * 6 == rank) then
        ! setup as global mesh
        nchunks_old = 6
        nproc_eta_old = int (sqrt ( dble(rank) / nchunks_old ))
        nproc_xi_old = int (sqrt ( dble(rank) / nchunks_old ))
        ! total procs
        nproctot_old = nproc_xi_old * nproc_eta_old * nchunks_old
      endif
    endif

    ! user output
    print *,'  estimated  number of total processes = ',rank
    print *,'  calculated number of total processes = ',nproc_xi_old * nproc_eta_old * nchunks_old
    print *

    ! NEX (must be multiple of 16 and 8 * multiple of nproc)
    nex_xi_old = 8 * nproc_xi_old
    ! finds maximum nex value to have nspec divisable by nex
    do while (mod(nspec_old,2*nex_xi_old) == 0)
      nex_xi_old = nex_xi_old * 2
    enddo
    ! unfortunatly, that doesn't work for all mesh setups...
    !
    ! the number of elements in the mesh is quite a complicated sum due to doubling layers in the mesh.
    ! so, there is no easy way to find out NEX based on nproc and nspec values.
    ! let's omit user output to avoid confusing in case the estimate is wrong
    !print *,'  estimated  NEX = ',nex_xi_old
    !print *
    ! this NEX value is only needed if we want estimate the typical size of an element in the old mesh.
    ! we can though base the estimate also on the new mesh where we have the NEX value through compilation values.

    ! checks
    if (rank == 0) then
      print *,'Error invalid number of processes found for old setup'
      stop 'Error invalid number of processes for old setup'
    endif
    if (nproc_xi_old * nproc_eta_old * nchunks_old /= rank) then
      print *,'Error: invalid old mesh setting, number of processes ',nproc_xi_old * nproc_eta_old * nchunks_old,' found ',rank
      print *,'Consider setting NCHUNKS value as input argument for old mesh. exiting...'
      stop 'Error invalid number of processes for old setup found'
    endif

  endif ! main

  ! main broadcasts to all other processes (assumes all slices have equal nspec/nglob values)
  call bcast_all_singlei(nspec_old)
  call bcast_all_singlei(nglob_old)
  call bcast_all_singlei(nproctot_old)
  call bcast_all_singlei(nproc_eta_old)
  call bcast_all_singlei(nproc_xi_old)
  call bcast_all_singlei(nchunks_old)

  ! check
  if (nproc_xi_old /= nproc_eta_old) then
    if (myrank == 0) then
      print *
      print *,'Warning: nproc_eta = ',nproc_eta_old,' not equal to nproc_xi = ',nproc_xi_old
      print *,'         routine works only by assuming nproc_xi == nproc_eta...'
      print *,'Please consider changing your old mesh setup!'
      print *
    endif
    stop 'Invalid NPROC_XI/NPROC_ETA for old mesh'
  endif

  ! total number of processes per chunk (source mesh)
  nproc_chunk1 = nproc_eta_old * nproc_xi_old

  ! note: in case nproc_xi and nproc_eta are the same for old and new mesh, then the slice will have the same mesh geometry.
  !       thus, we can only read in a single process mesh. otherwise, we will read in a full chunk and then search for
  !       the closest mesh point in the whole chunk mesh. however, this will require much more memory to store the whole chunk.
  if (nproc_xi_old == nproc_xi_new .and. nproc_eta_old == nproc_eta_new .and. nchunks_old == nchunks_new) then
    use_single_process_per_chunk = .true.
    nproc_chunk1 = 1
  else
    use_single_process_per_chunk = .false.
  endif
  ! out-comment to enforce full chunk search
  !use_single_process_per_chunk = .false.
  !nproc_chunk1 = nproc_eta_old * nproc_xi_old

  ! defines model parameters
  ! note: the fname setup is done here now to avoid corruption by the above adios routines. the adios calls above
  !       for some reason corrupt the fname list if this setup is done before...
  if (USE_TRANSVERSE_ISOTROPY) then
    ! transversly isotropic (TI) model
    nparams = 6
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    fname(1:6) = (/character(len=16) :: "reg1/vpv","reg1/vph","reg1/vsv","reg1/vsh","reg1/eta","reg1/rho"/)
#else
    fname(1:6) = (/character(len=16) :: "vpv","vph","vsv","vsh","eta","rho"/)
#endif

  else if (USE_AZIMUTHAL_ANISOTROPY) then
     nparams = 9

#ifdef USE_ADIOS_INSTEAD_OF_MESH
     fname(1:9) = (/character(len=16) :: "reg1/vpv","reg1/vph","reg1/vsv","reg1/vsh","reg1/eta","reg1/rho", &
                                         "reg1/mu0","reg1/Gc_prime","reg1/Gs_prime" /)
#else
     fname(1:9) = (/character(len=16) :: "vpv","vph","vsv","vsh","eta","rho", &
                                         "mu0","Gc_prime","Gs_prime" /)
#endif

  else
    ! isotropic model
    nparams = 3
    ! note: adds space at endings to equal number of characters
    !       to avoid compiler error: "Different shape for array assignment.."
    fname(1:3) = (/character(len=16) :: "vp ","vs ","rho"/)
  endif

  ! adds shear attenuation
  if (USE_ATTENUATION_Q) then
    nparams = nparams + 1
    fname(nparams) = "qmu"
  endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! adios only for model_gll.bp implemented which uses vpv,vph,..,rho
  if (nparams /= 6 .and. nparams /= 9) &
    stop 'ADIOS version only works for purely transversely isotropic and azimuthally anisotropic model file model_gll so far...'
#endif

  ! console output
  if (myrank == 0) then
    print *,'source mesh:  '
    print *,'  total number of processors = ',nproctot_old
    print *,'  nchunks    = ',nchunks_old
    print *,'  nproc_eta / nproc_xi = ',nproc_eta_old,nproc_xi_old
    print *,'  nspec      = ',nspec_old
    print *,'  nglob      = ',nglob_old
    print *
    print *,'target mesh:  '
    print *,'  total number of processors = ',nproctot_new
    print *,'  nchunks    = ',nchunks_new
    print *,'  nproc_eta / nproc_xi = ',nproc_eta_new,nproc_xi_new
    print *,'  nex        = ',nex_xi_new
    print *,'  nspec      = ',nspec_new
    print *,'  nglob      = ',nglob_new
    print *
    if (USE_TRANSVERSE_ISOTROPY) then
      print *,'model parameters:',nparams,' - transversely isotropic model'
    else if (USE_AZIMUTHAL_ANISOTROPY) then
      print *,'model parameters:',nparams,' - azimuthally anisotropic model'
    else
      print *,'model parameters:',nparams,' - isotropic model'
    endif
    if (USE_ATTENUATION_Q) then
      print *,'  includes qmu model parameter'
    endif
    print *,'  ( ',(trim(fname(iker))//" ",iker = 1,nparams),')'
    print *
    print *,'input model  directory: ',trim(input_model_dir)
    print *,'output model directory: ',trim(output_model_dir)
    print *
    print *,'array size:'
    print *,'  ibool1   = ',NGLLX*NGLLY*NGLLZ*nspec_old*nproc_chunk1*dble(SIZE_INTEGER)/1024./1024.,'MB'
    print *,'  x1,y1,z1 = ',nglob_old*nproc_chunk1*dble(CUSTOM_REAL)/1024./1024.,'MB'
    sizeval = NGLLX*NGLLY*NGLLZ*nspec_old*nproc_chunk1*nparams*dble(CUSTOM_REAL)
    print *,'  model1   = ',sngl(sizeval/1024./1024.),'MB'
    print *
    sizeval = NGLLX*NGLLY*NGLLZ*nspec_new*nparams*dble(CUSTOM_REAL)
    print *,'  model2   = ',sngl(sizeval/1024./1024.),'MB'
    print *
    print *,'total MPI processes: ',sizeprocs
    print *
    if (DO_BRUTE_FORCE_SEARCH) then
      print *,'location search by : brute-force approach'
    else
      print *,'location search by : kd-tree search'
      if (USE_MIDPOINT_SEARCH) then
        print *,'  uses midpoints of elements only'
      else
        print *,'  uses internal GLL points'
      endif
      if (DO_SEPARATION_410_650) then
        print *,'  uses element separation for 410-km/650-km discontinuity'
      endif
      if (DO_SEPARATION_TOPO) then
        print *,'  uses element separation for surface (moho) discontinuity'
      endif
      if (USE_FALLBACK) then
        print *,'  uses fall-back to model value of closest point in case of large differences'
      endif
    endif
    print *
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    print *,'file format: ADIOS files'
    print *
#endif
  endif
  call synchronize_all()

  ! checks
  if (sizeprocs /= nproctot_new) stop 'Error target mesh processors not equal to current total MPI processes'

  ! checks interpolation from regional to global model not supported yet
  if (nchunks_new == 6 .and. nchunks_old /= 6) then
    print *,'Error: interpolation from a regional mesh (old) to a global mesh (new) is not supported yet'
    stop 'Error interpolation from regional to global mesh'
  endif

  ! checks temporary file creation, to see if we could write out new model
  if (myrank == 0) then
    write(m_file,'(a,i6.6,a)') trim(output_model_dir)// '/tmp_proc',myrank,'.tmp'
    open(IOUT,file=trim(m_file),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(m_file)
      stop 'Error opening new output model file, please check if output directory exists...'
    endif
    close(IOUT,status='delete')
  endif
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *,'allocating search arrays:'
    print *,'  collecting number of processes per chunk = ',nproc_chunk1
    print *,'  use search by single process per chunk   = ',use_single_process_per_chunk
    print *
  endif

  ! collected mesh arrays for a single chunk
  allocate( x1(nglob_old,0:nproc_chunk1-1), &
            y1(nglob_old,0:nproc_chunk1-1), &
            z1(nglob_old,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating locations'
  x1(:,:) = 0.0_CUSTOM_REAL
  y1(:,:) = 0.0_CUSTOM_REAL
  z1(:,:) = 0.0_CUSTOM_REAL

  allocate( ibool1(NGLLX,NGLLY,NGLLZ,nspec_old,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating ibool1'
  ibool1(:,:,:,:,:) = 0

  allocate( idoubling1(nspec_old,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating idoubling1'
  idoubling1(:,:) = 0

  allocate( addressing1(nchunks_old,0:nproc_xi_old-1,0:nproc_eta_old-1),stat=ier )
  if (ier /= 0) stop 'Error allocating addressing1'
  addressing1(:,:,:) = 0

  ! model files
  allocate(addressing2(nchunks_new,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1), stat=ier)
  if (ier /= 0) stop 'Error allocating addressing2'
  addressing2(:,:,:) = 0

  allocate( model1(NGLLX,NGLLY,NGLLZ,nspec_old,nparams,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating initial model1'
  model1(:,:,:,:,:,:) = 0.0_CUSTOM_REAL

  allocate( model2(NGLLX,NGLLY,NGLLZ,nspec_new,nparams),stat=ier )
  if (ier /= 0) stop 'Error allocating target model2'
  model2(:,:,:,:,:) = 0.0_CUSTOM_REAL

  allocate( point_distance(NGLLX,NGLLY,NGLLZ,nspec_new),stat=ier )
  if (ier /= 0) stop 'Error allocating target model2'
  point_distance(:,:,:,:) = HUGEVAL_SNGL

  ! statistics
  allocate( model_maxdiff(nparams),stat=ier)
  if (ier /= 0) stop 'Error allocating model_maxdiff'
  model_maxdiff(:) = 0.0_CUSTOM_REAL

  ! new model topo
  allocate(x2(nglob_new), &
           y2(nglob_new), &
           z2(nglob_new), &
           ibool2(NGLLX,NGLLY,NGLLZ,nspec_new), &
           idoubling2(nspec_new), stat=ier)
  if (ier /= 0) stop 'Error allocating mesh arrays'
  x2(:) = 0.0_CUSTOM_REAL
  y2(:) = 0.0_CUSTOM_REAL
  z2(:) = 0.0_CUSTOM_REAL
  ibool2(:,:,:,:) = 0
  idoubling2(:) = 0

  ! synchronizes
  call synchronize_all()

  ! GLL points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! define indices of the control element
  call hex_nodes_anchor_ijk(anchor_iax,anchor_iay,anchor_iaz)

  ! mesh size
  ! compute typical size of elements at the surface
  ! note: we take NEX from new (target) mesh, as it is difficult to accurately determine NEX for the old mesh.
  !       since we only need this typical size as a search radius, having a rough estimate should be still fine.
  typical_size = TWO_PI * R_UNIT_SPHERE / ( 4.0 * nex_xi_new )

  ! user output
  if (myrank == 0) then
    print *,'  typical element size (new mesh at surface) = ',sngl(typical_size * R_PLANET_KM),'(km)'
    print *
  endif

  ! use 2 times the distance as a criterion for point detection
  typical_size = 2.0 * typical_size

  ! creates addressing
  ! old mesh addressing
  do ichunk = 1,nchunks_old
    do iproc_eta = 0,nproc_eta_old-1
      do iproc_xi = 0,nproc_xi_old-1
        iprocnum = (ichunk-1)*nproc_xi_old*nproc_eta_old + iproc_eta * nproc_xi_old + iproc_xi
        addressing1(ichunk,iproc_xi,iproc_eta) = iprocnum
      enddo
    enddo
  enddo
  ! new target mesh addressing
  do ichunk = 1,nchunks_new
    do iproc_eta = 0,NPROC_ETA_VAL-1
      do iproc_xi = 0,NPROC_XI_VAL-1
        iprocnum = (ichunk-1)*NPROC_XI_VAL*NPROC_ETA_VAL + iproc_eta * NPROC_XI_VAL + iproc_xi
        addressing2(ichunk,iproc_xi,iproc_eta) = iprocnum
      enddo
    enddo
  enddo

  ! gets associated chunk for this rank process
  ichunk_selected = 0
  iproc_eta_selected = 0
  iproc_xi_selected = 0
  do ichunk = 1,nchunks_new
    do iproc_eta = 0, NPROC_ETA_VAL - 1
      do iproc_xi = 0, NPROC_XI_VAL - 1
        ! gets slice number
        rank = addressing2(ichunk,iproc_xi,iproc_eta)

        ! only associated MPI process continues
        if (myrank == rank) then
          ichunk_selected = ichunk
          iproc_eta_selected = iproc_eta
          iproc_xi_selected = iproc_xi
          exit
        endif
      enddo
    enddo
  enddo
  if (ichunk_selected < 1 .or. ichunk_selected > nchunks_new ) stop 'Error selecting ichunk'
  ! debug
  !print *, 'selected chunk: ',ichunk_selected,' - eta/xi : ',iproc_eta_selected,iproc_xi_selected

  ! reads in the topology files of the target slices
  ! gets slice number
  rank = addressing2(ichunk_selected,iproc_xi_selected,iproc_eta_selected)

  ! only associated MPI process continues
  if (myrank /= rank) stop 'Error selected addressing rank'

  ! user output
  if (myrank == 0) then
    print *
    print *,'reading new mesh slice ... '
  endif

  ! checks new mesh locations
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! single adios file for all process slices
  ! re-initiate new group
  call init_adios_group(myadios_group,"InterpolatorNew")

  ! opens adios file
  solver_file = get_adios_filename(trim(dir_topo2)//'/solver_data')
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,solver_file)

  ! reads in scalars for rank
  call read_adios_scalar(myadios_file,myadios_group,rank,"reg1/nspec",nspec)
  call read_adios_scalar(myadios_file,myadios_group,rank,"reg1/nglob",nglob)

  ! checks dimensions
  if (nspec /= NSPEC_CRUST_MANTLE .or. nglob /= NGLOB_CRUST_MANTLE) then
    print *,'Error dimension of new mesh: solver_data nspec/nglob = ',nspec,nglob
    stop 'Error new mesh dimensions'
  endif
  call synchronize_all()

  ! reads in arrays
  call read_adios_array(myadios_file,myadios_group,rank,nglob,"reg1/x_global",x2(:))
  call read_adios_array(myadios_file,myadios_group,rank,nglob,"reg1/y_global",y2(:))
  call read_adios_array(myadios_file,myadios_group,rank,nglob,"reg1/z_global",z2(:))
  call read_adios_array(myadios_file,myadios_group,rank,nspec,"reg1/ibool",ibool2(:,:,:,:))
  call read_adios_array(myadios_file,myadios_group,rank,nspec,"reg1/idoubling",idoubling2(:))

  ! closes file
  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"InterpolatorNew")

#else
  ! opens binary file
  write(solver_file,'(a,i6.6,a)') trim(dir_topo2)//'/proc',rank,'_reg1_'//'solver_data.bin'
  open(IIN,file=solver_file,status='old',form='unformatted',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening file: ',trim(solver_file)
    stop 'Error opening new solver_data.bin file'
  endif
  read(IIN) nspec
  read(IIN) nglob

  ! checks dimensions
  if (nspec /= NSPEC_CRUST_MANTLE .or. nglob /= NGLOB_CRUST_MANTLE) then
    print *,'Error dimension of new mesh: solver_data nspec/nglob = ',nspec,nglob
    stop 'Error new mesh dimensions'
  endif

  ! locations
  read(IIN) x2(:)
  read(IIN) y2(:)
  read(IIN) z2(:)
  read(IIN) ibool2(:,:,:,:)
  read(IIN) idoubling2(:)
  close(IIN)
#endif
  if (myrank == 0) print *,'  locations done'

  ! user output
  if (myrank == 0) then
    print *
    print *, 'loading source mesh ... '
  endif

  ! reads in model and locations of old, source mesh
  ! combines all slices for this whole chunk

  ! sets loop bounds
  if (use_single_process_per_chunk) then
    istart_eta = iproc_eta_selected
    iend_eta = iproc_eta_selected

    istart_xi = iproc_xi_selected
    iend_xi = iproc_xi_selected

    istart_chunk = ichunk_selected
    iend_chunk = ichunk_selected
  else
    istart_eta = 0
    iend_eta = nproc_eta_old - 1

    istart_xi = 0
    iend_xi = nproc_xi_old - 1

    if (nchunks_old == nchunks_new) then
      istart_chunk = ichunk_selected
      iend_chunk = ichunk_selected
    else
      ! loops over all chunks
      istart_chunk = 1
      iend_chunk = nchunks_old
    endif
  endif

  ! timing
  tstart = wtime()

  ! loops over chunks
  do ichunk = istart_chunk,iend_chunk
    ! user output
    if (myrank == 0) then
      print *,'chunk number: ',ichunk,' out of ',iend_chunk - istart_chunk + 1
    endif
    call synchronize_all()

    ! initializes chunk mesh
    x1(:,:) = 0.0_CUSTOM_REAL
    y1(:,:) = 0.0_CUSTOM_REAL
    z1(:,:) = 0.0_CUSTOM_REAL
    ibool1(:,:,:,:,:) = 0
    idoubling1(:,:) = 0

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! single adios file for all old process slices
    ! re-initiate group
    call init_adios_group(myadios_group,"InterpolatorOld")
    ! opens adios file
    solver_file = get_adios_filename(trim(dir_topo1)//'/solver_data')
    call open_file_adios_read_and_init_method(myadios_file,myadios_group,solver_file)
#endif

    iprocnum = 0
    do iproc_eta = istart_eta, iend_eta
      do iproc_xi = istart_xi, iend_xi
        ! gets slice number
        rank = addressing1(ichunk,iproc_xi,iproc_eta)

        !debug
        !print *,'debug: chunk',ichunk,'myrank',myrank,'has process',rank,'iprocnum',iprocnum,'iproc xi/eta',iproc_xi,iproc_eta

        ! counter
        iprocnum = iprocnum + 1

        ! user output
        if (myrank == 0) then
          print *,'  slice number: ',iprocnum,' out of ',nproc_chunk1
        endif

        ! warning
        if (use_single_process_per_chunk) then
          if (rank /= myrank) then
            print *,'Warning: old mesh rank ',rank,' not identical to target rank ',myrank
          endif
        endif

        ! reads in old arrays
#ifdef USE_ADIOS_INSTEAD_OF_MESH
        ! ADIOS
        ! reads in scalars for rank
        call read_adios_scalar(myadios_file,myadios_group,rank,"reg1/nspec",nspec)
        call read_adios_scalar(myadios_file,myadios_group,rank,"reg1/nglob",nglob)

        ! checks dimensions
        if (nspec /= nspec_old .or. nglob /= nglob_old) then
          print *,'Error dimension of old, source mesh: solver_data nspec/nglob = ',nspec,nglob
          stop 'Error new mesh dimensions'
        endif

        ! reads in arrays
        call read_adios_array(myadios_file,myadios_group,rank,nglob,"reg1/x_global",x1(:,iprocnum-1))
        call read_adios_array(myadios_file,myadios_group,rank,nglob,"reg1/y_global",y1(:,iprocnum-1))
        call read_adios_array(myadios_file,myadios_group,rank,nglob,"reg1/z_global",z1(:,iprocnum-1))

        call read_adios_array(myadios_file,myadios_group,rank,nspec,"reg1/ibool",ibool1(:,:,:,:,iprocnum-1))
        call read_adios_array(myadios_file,myadios_group,rank,nspec,"reg1/idoubling",idoubling1(:,iprocnum-1))

#else
        ! old, source mesh locations
        write(solver_file,'(a,i6.6,a)') trim(dir_topo1)//'/proc',rank,'_reg1_'//'solver_data.bin'
        open(IIN,file=solver_file,status='old',form='unformatted',action='read',iostat=ier)
        if (ier /= 0) then
          print *,'Error opening file: ',trim(solver_file)
          stop 'Error opening old solver_data.bin file'
        endif
        read(IIN) nspec
        read(IIN) nglob

        ! checks dimensions
        if (nspec /= nspec_old .or. nglob /= nglob_old) then
          print *,'Error dimension of old, source mesh: solver_data nspec/nglob = ',nspec,nglob
          stop 'Error new mesh dimensions'
        endif

        read(IIN) x1(:,iprocnum-1)
        read(IIN) y1(:,iprocnum-1)
        read(IIN) z1(:,iprocnum-1)
        read(IIN) ibool1(:,:,:,:,iprocnum-1)
        read(IIN) idoubling1(:,iprocnum-1)
        close(IIN)
#endif
      enddo
    enddo

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! closes file
    call close_file_adios_read_and_finalize_method(myadios_file)
    call delete_adios_group(myadios_group,"InterpolatorOld")
#endif

    ! user output
    if (myrank == 0) then
      print *
      print *,'  source mesh chunk read successfully'
      print *
    endif
    call synchronize_all()

    ! user output
    if (myrank == 0) then
      print *, '  loading source model ... '
      !do iker = 1,nparams
      !  print *, '  for parameter: ',trim(fname(iker))
      !enddo
      !print *
    endif

    ! reads in old model files
    model1(:,:,:,:,:,:) = 0.0_CUSTOM_REAL

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! single adios file for all old model arrays
    ! opens adios file with old model values
    solver_file = get_adios_filename(trim(input_model_dir)//'/model_gll')
    call open_file_adios_read_and_init_method(myadios_val_file,myadios_val_group,solver_file)
#endif

    iprocnum = 0
    do iproc_eta = istart_eta, iend_eta
      do iproc_xi = istart_xi, iend_xi
        ! gets slice number
        rank = addressing1(ichunk,iproc_xi,iproc_eta)

        ! counter
        iprocnum = iprocnum + 1

        ! user output
        if (myrank == 0) then
          print *,'  slice number: ',iprocnum,' out of ',nproc_chunk1
        endif
        ! warning
        if (use_single_process_per_chunk) then
          if (rank /= myrank) then
            print *,'Warning: old mesh rank ',rank,' not identical to target rank ',myrank
          endif
        endif

        ! reads in model slices
        do iker = 1,nparams
          ! debug user output
          if (myrank == 0) print *, '  for parameter: ',trim(fname(iker))

#ifdef USE_ADIOS_INSTEAD_OF_MESH
          ! reads in array
          call read_adios_array(myadios_val_file,myadios_val_group,rank,nspec_old,fname(iker),model1(:,:,:,:,iker,iprocnum-1))
#else
          ! opens model file
          write(m_file,'(a,i6.6,a)') trim(input_model_dir)//'/proc',rank,'_reg1_'//trim(fname(iker))//'.bin'
          open(IIN,file=trim(m_file),status='old',form='unformatted',action='read',iostat=ier)
          if (ier /= 0) then
            print *,'Error opening file: ',trim(m_file)
            stop 'Error opening old model file'
          endif
          read(IIN) model1(:,:,:,:,iker,iprocnum-1)
          close(IIN)
#endif
        enddo ! nparams
      enddo
    enddo

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! closes file
    call close_file_adios_read_and_finalize_method(myadios_val_file)
#endif

    ! user output
    if (myrank == 0) then
      print *
      print *,'  source model chunk read successfully'
      print *
    endif
    call synchronize_all()

    ! determines if target mesh is a cut-off mesh with USE_LOCAL_MESH flag turned on
    ! note: For local meshes, all mesh elements belong to idoubling == IFLAG_CRUST == 1.
    !       Thus, interpolating from a "normal chunk mesh" with idoubling layers min/max
    !       from 1 to 5 (IFLAG_CRUST to IFLAG_MANTLE_NORMAL), and staying within the layer for the element search
    !       would only interpolate crustal values from the source mesh to deeper elements in the target mesh.
    !       We therefore need a flag to decide if we consider elements from all layers in the source mesh for these target meshes.
    !
    !       For cut-off meshes with only REGIONAL_MESH_CUTOFF set to .true., the idoubling flags min/max can go from 1 to 4
    !       (IFLAG_CRUST to IFLAG_670_220) depending on the REGIONAL_MESH_CUTOFF_DEPTH value.
    !       In these cut-off meshes, staying within layers for element search is still okay to properly account for discontinuities.
    if (minval(idoubling2) == 1 .and. maxval(idoubling2) == 1) then
      target_is_local_mesh = .true.
    else
      target_is_local_mesh = .false.
    endif

    ! checks that layers match
    if (myrank == 0) then
      print *,'mesh doubling layers:'
      print *,'  source mesh: min/max = ',minval(idoubling1),maxval(idoubling1)
      print *,'  target mesh: min/max = ',minval(idoubling2),maxval(idoubling2)
      print *
      if (target_is_local_mesh) then
        print *,'  target mesh: uses cut-off local mesh'
        print *
      endif
    endif

    ! source mesh must have at least all layers that the target mesh has (in case of cut-off meshes REGIONAL_MESH_CUTOFF)
    if (minval(idoubling1) > minval(idoubling2) .or. maxval(idoubling1) < maxval(idoubling2)) then
      print *,'Error idoubling range:'
      print *,'  idoubling 1:',minval(idoubling1),maxval(idoubling1)
      print *,'  idoubling 2:',minval(idoubling2),maxval(idoubling2)
      stop 'Error invalid idoubling range'
    endif

    ! user output
    if (myrank == 0) then
      print *,'Earth layers: ',minval(idoubling2),' to ',maxval(idoubling2)
      print *
    endif

    ! loops over layers (crust/moho-80/80-220/220-660/660-CMB)
    total_nspec = 0
    do ilayer = minval(idoubling2),maxval(idoubling2)

      ! user output
      if (myrank == 0) then
        print *,'layer: ',ilayer,' out of ',maxval(idoubling2)
        select case (ilayer)
        case (IFLAG_CRUST)
          print *,'layer: crust'
        case (IFLAG_80_MOHO)
          print *,'layer: 80 - MOHO'
        case (IFLAG_220_80)
          print *,'layer: 220 - 80'
        case (IFLAG_670_220)
          print *,'layer: 670 - 220'
        case (IFLAG_MANTLE_NORMAL)
          print *,'layer: mantle normal'
        end select
      endif

      ! statistics
      model_maxdiff(:) = 0.0_CUSTOM_REAL
      is_updated = .false.

      ! builds search tree
      if (.not. DO_BRUTE_FORCE_SEARCH) then
        ! counts total number of points in this layer in source mesh
        iprocnum = 0
        inodes = 0
        do iproc_eta = istart_eta, iend_eta
          do iproc_xi = istart_xi, iend_xi
            ! counter
            iprocnum = iprocnum + 1
            ! search elements
            do ispec = 1, nspec_old
              ! search elements
              if (.not. target_is_local_mesh) then
                ! skips elements outside of this layer
                if (idoubling1(ispec,iprocnum-1) /= ilayer ) cycle
              endif
              if (TREE_INTERNAL_GLL_POINTS) then
                ! all internal GLL points ( 2 to NGLLX-1 )
                inodes = inodes + (NGLLX-2)*(NGLLY-2)*(NGLLZ-2)
              endif
              if (TREE_MID_POINTS) then
                ! only element mid-points
                inodes = inodes + 1
              endif
            enddo
          enddo
        enddo
        ! checks
        if (inodes < 1 ) stop 'Error no search tree nodes in source mesh for this layer'
        ! checks maximum number of nodes
        k = 0
        if (TREE_INTERNAL_GLL_POINTS) k = nspec_old * nproc_chunk1 * (NGLLX-2)*(NGLLY-2)*(NGLLZ-2)
        if (TREE_MID_POINTS) k = nspec_old * nproc_chunk1 * 1
        if (inodes > k ) stop 'Error invalid number of search tree nodes in this layer'

        ! set number of tree nodes
        kdtree_num_nodes = inodes

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
        iprocnum = 0
        inodes = 0
        do iproc_eta = istart_eta, iend_eta
          do iproc_xi = istart_xi, iend_xi

            ! counter
            iprocnum = iprocnum + 1

            ! adds tree nodes
            do ispec = 1,nspec_old

              ! search elements
              if (.not. target_is_local_mesh) then
                ! skips elements outside of this layer
                if (idoubling1(ispec,iprocnum-1) /= ilayer ) cycle
              endif

              ! sets up tree nodes
              ! all internal GLL points
              if (TREE_INTERNAL_GLL_POINTS) then
                do k = 2,NGLLZ-1
                  do j = 2,NGLLY-1
                    do i = 2,NGLLX-1
                      iglob = ibool1(i,j,k,ispec,iprocnum-1)

                      ! counts nodes
                      inodes = inodes + 1
                      if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

                      ! adds node index ( index points to same ispec for all internal GLL points)
                      kdtree_nodes_index(inodes) = ispec + (iprocnum - 1) * nspec_old

                      ! adds node location
                      kdtree_nodes_location(1,inodes) = x1(iglob,iprocnum-1)
                      kdtree_nodes_location(2,inodes) = y1(iglob,iprocnum-1)
                      kdtree_nodes_location(3,inodes) = z1(iglob,iprocnum-1)
                    enddo
                  enddo
                enddo
              endif

              ! only element midpoints
              if (TREE_MID_POINTS) then
                iglob = ibool1(MIDX,MIDY,MIDZ,ispec,iprocnum-1)

                ! counts nodes
                inodes = inodes + 1
                if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

                ! adds node index ( index points to same ispec for all internal GLL points)
                kdtree_nodes_index(inodes) = ispec + (iprocnum - 1) * nspec_old

                ! adds node location
                kdtree_nodes_location(1,inodes) = x1(iglob,iprocnum-1)
                kdtree_nodes_location(2,inodes) = y1(iglob,iprocnum-1)
                kdtree_nodes_location(3,inodes) = z1(iglob,iprocnum-1)
              endif

            enddo
          enddo
        enddo
        if (inodes /= kdtree_num_nodes) stop 'Error index inodes does not match nnodes_local'
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

      endif

      ! user output
      if (myrank == 0) print *,'looping over elements:'

      ! loop over all elements (mainly those in this layer)
      do ispec = 1, nspec_new
        ! user output
        if (myrank == 0) then
          if (ispec == 1 .or. mod(ispec,int(0.1*nspec_new)) == 0 .or. ispec == nspec_new) then
            print *,'  ispec',ispec,' out of ',nspec_new
          endif
        endif

        ! skip elements out of this layer
        if (idoubling2(ispec) /= ilayer ) cycle

        ! increases element counter
        total_nspec = total_nspec + 1

        ! gets model values
        if (DO_BRUTE_FORCE_SEARCH) then
          ! brute-force search over all GLL points
          call get_model_values_bruteforce(ispec,nspec_new,nglob_new,ibool2,x2,y2,z2,nparams,model2, &
                                           nspec_old,nglob_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                           anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size,myrank, &
                                           model_maxdiff,point_distance,is_updated)
        else
          ! kdtree search
          call get_model_values_kdtree(ispec,nspec_new,nglob_new,ibool2,x2,y2,z2,nparams,model2, &
                                       nspec_old,nglob_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                       anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size,myrank, &
                                       model_maxdiff,point_distance,is_updated, &
                                       USE_MIDPOINT_SEARCH,DO_SEPARATION_410_650,DO_SEPARATION_TOPO,USE_FALLBACK)
        endif
      enddo ! ispec

      ! user output
      if (myrank == 0) print *

      ! frees tree memory
      if (.not. DO_BRUTE_FORCE_SEARCH) then
        ! deletes tree arrays
        deallocate(kdtree_nodes_location)
        deallocate(kdtree_nodes_index)
        ! deletes search tree nodes
        call kdtree_delete()
      endif

      ! statistics
      call any_all_l(is_updated,is_updated_all)
      if (is_updated_all) then
        ! user output
        if (myrank == 0) print *,'  statistics:'
        do iker = 1,nparams
          call max_all_cr(model_maxdiff(iker),val_all)
          if (myrank == 0) print *,'  parameter ',trim(fname(iker)),': maximum difference = ',val_all
        enddo
      else
        if (myrank == 0) print *,'  statistics: no updates'
      endif

      ! user output
      if (myrank == 0) print *
    enddo ! ilayer

    ! checks
    if (total_nspec /= nspec_new) then
      print *,'Error: invalid total number of elements',total_nspec,' should be ',nspec_new
      stop 'Error invalid total number of elements after loop'
    endif
    call synchronize_all()

    ! user output
    if (myrank == 0) then
      print *
      print *,'  search in chunk ',ichunk,' done'
      print *
    endif
    call synchronize_all()

  enddo  ! ichunks

  ! distance info
  val = maxval(point_distance(:,:,:,:))
  call max_all_all_cr(val,val_all)
  ! convert to km
  val_all = sqrt(val_all) * real(R_PLANET_KM,kind=CUSTOM_REAL)
  ! user output
  if (myrank == 0) then
    print *
    print *,'point distances have a maximum distance = ',val_all,'(km)'
    print *
  endif
  ! highlights maximum distances
  if (val_all >= 5.0) then
    if (myrank == 0) then
      print *,'***********************'
      print *,'WARNING: some point locations found have poor match'
      print *,'***********************'
      print *,'Due to topography of internal interfaces and doubling layers this might occur (often for coarse meshes).'
      print *,'Please double-check interpolation results.'
      print *
      print *,'More details:'
    endif
    ! outputs location with maximum distance
    do i = 0,nproctot_new-1
      if (i == myrank) then
        loc_max = maxloc(point_distance(:,:,:,:))
        print *,'poor location: rank ',myrank,'has maximum distance ',sqrt(val)*R_PLANET_KM,'(km)'
        print *,'     location: at index ',loc_max
        print *,'     location: doubling layer ',idoubling2(loc_max(4)),'(1=crust,2=80-MOHO,3=220-80,4=670-220,5=mantle-normal)'
        iglob = ibool2(loc_max(1),loc_max(2),loc_max(3),loc_max(4))
        val = sqrt(x2(iglob)*x2(iglob) + y2(iglob)*y2(iglob) + z2(iglob)*z2(iglob)) * real(R_PLANET_KM,kind=CUSTOM_REAL)
        print *,'     location: ',x2(iglob),y2(iglob),z2(iglob),'at radius ',val,'(km) depth',sngl(R_PLANET_KM - val),'(km)'
        print *
      endif
      call synchronize_all()
    enddo
  endif
  call synchronize_all()
  if (myrank == 0) print *

  ! frees memory
  deallocate(x1,y1,z1)
  deallocate(ibool1)
  deallocate(idoubling1)
  deallocate(model1)
  deallocate(model_maxdiff)
  deallocate(point_distance)
  deallocate(ibool2)
  deallocate(x2,y2,z2)

  ! user output
  if (myrank == 0) then
    print *, 'writing out new model files'
  endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! sets up adios group
  group_name = "MODELS_GROUP"
  call init_adios_group(myadios_val_group,group_name)

  ! defines group size
  group_size_inc = 0
  call define_adios_scalar(myadios_val_group, group_size_inc, '', "NSPEC", nspec_new)
  call define_adios_scalar(myadios_val_group, group_size_inc, '', "reg1/nspec", nspec_new)

  local_dim = size(model2(:,:,:,:,iker))
  do iker = 1,nparams
    call define_adios_global_array1D(myadios_val_group, group_size_inc,local_dim,'',trim(fname(iker)),model2(:,:,:,:,iker))
  enddo

  ! opens new adios model file
  solver_file = get_adios_filename(trim(output_model_dir) //'/model_gll_interpolated')
  call open_file_adios_write(myadios_val_file,myadios_val_group,solver_file,group_name)

  call set_adios_group_size(myadios_val_file,group_size_inc)

  ! writes nspec (for checking and backward compatibility)
  call write_adios_scalar(myadios_val_file,myadios_val_group,"NSPEC",nspec_new)
  call write_adios_scalar(myadios_val_file,myadios_val_group,"reg1/nspec",nspec_new)
#endif

  ! writes out new model
  do iker = 1,nparams
    ! user output
    if (myrank == 0) then
      print *, '  for parameter: ',trim(fname(iker))
      print *, '    slice rank 0 has min/max = ',minval(model2(:,:,:,:,iker)),'/',maxval(model2(:,:,:,:,iker))
    endif
    if (minval(model2(:,:,:,:,iker)) <= 0.001) then
      print *
      print *,'***********************'
      print *,'WARNING: minimum model value ',minval(model2(:,:,:,:,iker)),'is (almost) zero for slice ',myrank,' - check output'
      print *,'***********************'
      print *
    endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    ! writes previously defined ADIOS variables
    call write_adios_array_gll(myadios_val_file,myadios_val_group,myrank,sizeprocs_adios,nspec_new,trim(fname(iker)), &
                               model2(:,:,:,:,iker))
#else
    ! default binary
    write(m_file,'(a,i6.6,a)') trim(output_model_dir) // '/proc',myrank,'_reg1_'//trim(fname(iker))//'.bin'
    open(IOUT,file=trim(m_file),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(m_file)
      stop 'Error opening output model file'
    endif
    write(IOUT) model2(:,:,:,:,iker)
    close(IOUT)
#endif

  enddo ! nparams

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! Reset the path to its original value to avoid bugs.
  call write_adios_perform(myadios_val_file)
  ! closing performs actual write
  call close_file_adios(myadios_val_file)
#endif

  ! frees memory
  deallocate(model2)

  ! synchronizes MPI processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *
    print *, 'check new model files in directory: ',trim(output_model_dir)
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    print *, 'see file: ',trim(solver_file)
#else
    print *, 'see files: ',trim(output_model_dir) // '/proc***_reg1_**.bin'
#endif
    ! timing
    tCPU = wtime() - tstart
    print *
    print *,'Elapsed time for model interpolation in seconds = ',sngl(tCPU)
    print *
    print *, 'done successfully'
    print *
  endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! finalizes adios
  call finalize_adios()
#endif

  ! exiting MPI processes
  call finalize_mpi()

  end program interpolate_model

!
!------------------------------------------------------------------------------
!

  subroutine get_model_values_bruteforce(ispec,nspec_new,nglob_new,ibool2,x2,y2,z2,nparams,model2, &
                                         nspec_old,nglob_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                         anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size,myrank, &
                                         model_maxdiff,point_distance,is_updated)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD

  implicit none

  integer,intent(in) :: ispec
  integer,intent(in) :: nparams

  ! new, target mesh
  integer,intent(in) :: nspec_new,nglob_new
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec_new),intent(in) :: ibool2
  real(kind=CUSTOM_REAL),dimension(nglob_new),intent(in) :: x2,y2,z2
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_new,nparams),intent(inout) :: model2

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec_old,nglob_old,nproc_chunk1
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec_old,0:nproc_chunk1-1),intent(in) :: ibool1
  real(kind=CUSTOM_REAL),dimension(nglob_old,0:nproc_chunk1-1),intent(in) :: x1,y1,z1
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_old,nparams,0:nproc_chunk1-1),intent(in) :: model1

  ! indices of the control points of the surface element
  integer,intent(in) :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL),dimension(nparams),intent(inout) :: model_maxdiff
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_new),intent(inout) :: point_distance
  logical, intent(inout) :: is_updated

  ! local parameters
  integer :: i,j,k,iglob,iker
  ! interpolated point location
  integer :: ispec_selected,rank_selected
  integer :: i_selected,j_selected,k_selected
  double precision :: xi,eta,gamma
  ! point location
  real(kind=CUSTOM_REAL) :: x_target,y_target,z_target
  real(kind=CUSTOM_REAL) :: val,val_initial
  real(kind=CUSTOM_REAL) :: dist_sq

  ! checks given ispec
  if (ispec < 1 .or. ispec > nspec_new) then
    print *,'Error: rank ',myrank,' has invalid ispec'
    stop 'Error invalid ispec in get_model_values_bruteforce() routine'
  endif

  ! brute-force search over all points
  ! loops over all GLL points
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        ! target point location
        iglob = ibool2(i,j,k,ispec)
        x_target = x2(iglob)
        y_target = y2(iglob)
        z_target = z2(iglob)

        ! gets interpolated position
        call locate(x_target,y_target,z_target, &
                    xi,eta,gamma, &
                    ispec_selected,rank_selected, &
                    nspec_old,nglob_old,nproc_chunk1, &
                    ibool1,x1,y1,z1, &
                    anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size, &
                    i_selected,j_selected,k_selected, &
                    dist_sq)

        ! interpolate model values
        if (dist_sq < point_distance(i,j,k,ispec)) then
          ! better position found
          point_distance(i,j,k,ispec) = dist_sq
          is_updated = .true.

          ! new model value
          do iker = 1,nparams
            call interpolate_element_value(xi,eta,gamma,ispec_selected, &
                                           nspec_old,model1(:,:,:,ispec_selected,iker,rank_selected), &
                                           val,xigll,yigll,zigll)

            ! sets new model value
            model2(i,j,k,ispec,iker) = val
          enddo

          ! checks model value differences
          do iker = 1,nparams
            val_initial = model1(i_selected,j_selected,k_selected,ispec_selected,iker,rank_selected)
            val = model2(i,j,k,ispec,iker)

            ! statistics
            if (abs(val - val_initial ) > abs(model_maxdiff(iker))) model_maxdiff(iker) = val - val_initial
          enddo
        endif

      enddo
    enddo
  enddo

  end subroutine get_model_values_bruteforce


!
!------------------------------------------------------------------------------
!

  subroutine get_model_values_kdtree(ispec,nspec_new,nglob_new,ibool2,x2,y2,z2,nparams,model2, &
                                     nspec_old,nglob_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                     anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size,myrank, &
                                     model_maxdiff,point_distance,is_updated, &
                                     USE_MIDPOINT_SEARCH,DO_SEPARATION_410_650,DO_SEPARATION_TOPO,USE_FALLBACK)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,MIDX,MIDY,MIDZ

  use shared_parameters, only: R220,R400,R600,R670,R_PLANET,R_PLANET_KM

  use kdtree_search, only: kdtree_find_nearest_neighbor,kdtree_nodes_location

  implicit none

  integer,intent(in) :: ispec
  integer,intent(in) :: nparams

  ! new, target mesh
  integer,intent(in) :: nspec_new,nglob_new
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec_new),intent(in) :: ibool2
  real(kind=CUSTOM_REAL),dimension(nglob_new),intent(in) :: x2,y2,z2
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_new,nparams),intent(inout) :: model2

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec_old,nglob_old,nproc_chunk1
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec_old,0:nproc_chunk1-1),intent(in) :: ibool1
  real(kind=CUSTOM_REAL),dimension(nglob_old,0:nproc_chunk1-1),intent(in) :: x1,y1,z1
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_old,nparams,0:nproc_chunk1-1),intent(in) :: model1

  ! indices of the control points of the surface element
  integer,intent(in) :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL),dimension(nparams),intent(inout) :: model_maxdiff
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_new),intent(inout) :: point_distance
  logical, intent(inout) :: is_updated

  logical,intent(in) :: USE_MIDPOINT_SEARCH
  logical,intent(in) :: DO_SEPARATION_410_650,DO_SEPARATION_TOPO
  logical,intent(in) :: USE_FALLBACK

  ! local parameters
  integer :: i,j,k,iglob,iker,iglob1,iglob2
  ! interpolated point location
  integer :: ispec_selected,rank_selected
  double precision :: xi,eta,gamma
  ! point location
  real(kind=CUSTOM_REAL) :: x_target,y_target,z_target
  ! nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: dist_min
  double precision :: elem_height,r1,r2,r_sq
  integer :: iglob_min
  integer :: i_selected,j_selected,k_selected
  ! locations
  real(kind=CUSTOM_REAL) :: x_found,y_found,z_found
  real(kind=CUSTOM_REAL) :: val,val_initial
  real(kind=CUSTOM_REAL) :: dist_sq
  logical :: is_critical,search_internal
  integer :: ii,jj,kk


  !daniel todo: check if these margins would need to be adapted for mars
  !
  ! radius in m
  double precision :: RTOP
  ! note: for s362ani, 410-discontinuity is defined at R400, and 650-discontinuity at R670
  ! margins
  ! surface topography and moho stretching: uses a 120-km margin
  ! (moho max ~ 80 km below Himalayan)
  double precision,parameter :: RTOP_MARGIN = 120000.d0
  ! 410-km: uses a 50-km margin
  ! (s362ani: 410 topography perturbations have min/max ~ -13/+13 km)
  ! (element heights are around 60km below and 64km above the 410-km discontinuity)
  double precision,parameter :: R400_MARGIN = 50000.d0
  ! 650-km: uses a 50-km margin
  ! (s362ani: 650 topography perturbations have min/max ~ -14/+19 km)
  double precision,parameter :: R670_MARGIN = 50000.d0

  ! debug warning about large model value differences
  logical,parameter :: DO_WARNING = .false.
  ! debugging
  logical,parameter :: DEBUG = .false.

  ! checks given ispec
  if (ispec < 1 .or. ispec > nspec_new) then
    print *,'Error: rank ',myrank,' has invalid ispec'
    stop 'Error invalid ispec in get_model_values_kdtree() routine'
  endif

  ! initializes for 410/660/topo special case
  is_critical = .false.

  ! searches for element using mid-point
  if (USE_MIDPOINT_SEARCH .or. DO_SEPARATION_410_650 .or. DO_SEPARATION_TOPO) then
    ! element mid-point
    iglob = ibool2(MIDX,MIDY,MIDZ,ispec)

    xyz_target(1) = x2(iglob)
    xyz_target(2) = y2(iglob)
    xyz_target(3) = z2(iglob)

    ! finds closest point in source chunk
    if (USE_MIDPOINT_SEARCH) then
      call kdtree_find_nearest_neighbor(xyz_target,iglob_min,dist_min)

      ! selected source rank
      rank_selected = int((iglob_min-1)/nspec_old)
      ! selected closest element
      ispec_selected = iglob_min - rank_selected * nspec_old

      ! checks if mid-point was found properly
      call check_point_result()

      ! debug
      !if (myrank == 0 .and. iglob < 100) &
      !  print *,'dist_min kdtree midpoint: ',dist_min * R_PLANET_KM,'(km)',typical_size * R_PLANET_KM

    else
      ! uses all GLL points to find the best element
      !
      ! checks if element close to surface/410/660 interfaces
      !
      ! note: we want to make sure that for points on a discontinuity/interface,
      !       the interpolation is done from elements on the correct side.
      !
      !       we thus flag elements here as critical in case they are close to such an interface.
      !       elements with this flag will then avoid locating the closest GLL point for points being on a surface,
      !       but rather take only internal element GLL points for the search of the best element.
      !
      !       for midpoint search, we located the best element using element midpoints only.
      !       in this case, we are probably fine by finding the closest element on the correct side of an interface.
      !
      ! make sure elements close to an interface are flagged as "critical" for special care
      if (DO_SEPARATION_410_650 .or. DO_SEPARATION_TOPO) then
        ! point radius (squared)
        r_sq = xyz_target(1)*xyz_target(1) + xyz_target(2)*xyz_target(2) + xyz_target(3)*xyz_target(3)

        ! special case for 410-km/650-km discontinuity
        if (DO_SEPARATION_410_650) then
          ! surface
          if (r_sq >= (R220/R_PLANET)**2) then
            ! elements close to surface
            is_critical = .true.
          endif
          ! 410-km discontinuity
          if (r_sq >= (R600/R_PLANET)**2 .and. r_sq <= (R220/R_PLANET)**2) then
            ! elements within 220 - 600 km depth
            is_critical = .true.
          endif
          ! 650-km discontinuity
          if (r_sq >= ((R670 - R670_MARGIN)/R_PLANET)**2 .and. r_sq <= ((R670 + R670_MARGIN)/R_PLANET)**2) then
            ! elements within around 650 km depth
            is_critical = .true.
          endif
        endif

        ! special case for surface (moho) discontinuity
        if (DO_SEPARATION_TOPO) then
          ! surface
          if (r_sq >= (R220/R_PLANET)**2) then
            ! elements close to surface
            is_critical = .true.
          endif
        endif

        ! debugging
        if (DEBUG .and. is_critical) then
          ! element height: size along a vertical edge
          ! top point
          iglob1 = ibool2(1,1,NGLLZ,ispec)
          r1 = sqrt(x2(iglob1)*x2(iglob1) + y2(iglob1)*y2(iglob1) + z2(iglob1)*z2(iglob1))
          ! bottom point
          iglob2 = ibool2(1,1,1,ispec)
          r2 = sqrt(x2(iglob2)*x2(iglob2) + y2(iglob2)*y2(iglob2) + z2(iglob2)*z2(iglob2))
          ! element height
          elem_height = r1 - r2
          ! debug
          if (myrank == 0) print *,'element height: ',elem_height * R_PLANET_KM,'(km)','radius: ',sqrt(r_sq)*R_PLANET_KM
        endif
      endif
    endif
  endif

  ! loops over all element GLL points
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        ! target point location
        iglob = ibool2(i,j,k,ispec)

        ! target point location
        x_target = x2(iglob)
        y_target = y2(iglob)
        z_target = z2(iglob)

        ! kd-search for this single GLL point
        if (.not. USE_MIDPOINT_SEARCH) then
          ! finds closest element in original mesh
          ! search using only internal GLL points (avoid GLL points on the surface of the element)
          search_internal = .false.

          ! avoids getting values from "wrong" side on 410-km discontinuity,etc.
          if (DO_SEPARATION_410_650 .or. DO_SEPARATION_TOPO) then
            ! only points on surface of element
            if (i == 1 .or. i == NGLLX .or. j == 1 .or. j == NGLLY .or. k == 1 .or. k == NGLLZ) then
              ! further limits the number of internal searches
              ! checks each critical element if the points are close to an interface/discontinuity
              if (is_critical) then
                ! GLL point radius (squared)
                r_sq = x_target*x_target + y_target*y_target + z_target*z_target

                if (DO_SEPARATION_410_650) then
                  ! takes corresponding internal GLL point for element search
                  ! 410-km discontinuity
                  if (r_sq >= ((R400 - R400_MARGIN)/R_PLANET)**2 .and. &
                      r_sq <= ((R400 + R400_MARGIN)/R_PLANET)**2) search_internal = .true.
                  ! 650-km discontinuity
                  if (r_sq >= ((R670 - R670_MARGIN)/R_PLANET)**2 .and. &
                      r_sq <= ((R670 + R670_MARGIN)/R_PLANET)**2) search_internal = .true.
                endif

                if (DO_SEPARATION_TOPO) then
                  RTOP = R_PLANET
                  ! takes corresponding internal GLL point for element search
                  ! surface elements
                  if (r_sq >= ((RTOP - RTOP_MARGIN)/R_PLANET)**2) search_internal = .true.
                endif
              endif
            endif
          endif

          ! avoid using nodes at upper/lower/..outer surfaces
          if (search_internal) then
            ! new search point indices
            ii = i
            jj = j
            kk = k

            ! takes corresponding internal one for element search
            if (i == 1) then
              ii = 2
            else if (i == NGLLX) then
              ii = NGLLX - 1
            endif
            if (j == 1) then
              jj = 2
            else if (j == NGLLY) then
              jj = NGLLY - 1
            endif
            if (k == 1) then
              kk = 2
            else if (k == NGLLZ) then
              kk = NGLLZ - 1
            endif

            ! new target point location
            iglob = ibool2(ii,jj,kk,ispec)
            x_target = x2(iglob)
            y_target = y2(iglob)
            z_target = z2(iglob)
          endif

          ! kdtree search for each single GLL point
          xyz_target(1) = x_target
          xyz_target(2) = y_target
          xyz_target(3) = z_target

          ! finds closest point in source chunk
          call kdtree_find_nearest_neighbor(xyz_target,iglob_min,dist_min)

          ! selected source rank
          rank_selected = int((iglob_min-1)/nspec_old)
          ! selected closest element
          ispec_selected = iglob_min - rank_selected * nspec_old

          ! checks if point was found properly
          call check_point_result()

          ! debug
          !if (myrank == 0 .and. iglob < 100) &
          !  print *,'dist_min kdtree: ',dist_min * R_PLANET_KM,'(km)',typical_size * R_PLANET_KM

          ! restores original target point location for locating/interpolating
          iglob = ibool2(i,j,k,ispec)
          x_target = x2(iglob)
          y_target = y2(iglob)
          z_target = z2(iglob)
        endif

        ! gets interpolated position within selected element
        call locate_single(x_target,y_target,z_target, &
                           xi,eta,gamma, &
                           ispec_selected, &
                           nspec_old,nglob_old, &
                           ibool1(:,:,:,:,rank_selected),x1(:,rank_selected),y1(:,rank_selected),z1(:,rank_selected), &
                           anchor_iax,anchor_iay,anchor_iaz,xigll,yigll,zigll,typical_size, &
                           i_selected,j_selected,k_selected, &
                           dist_sq)

        ! interpolate model values
        if (dist_sq < point_distance(i,j,k,ispec)) then
          ! better position found
          point_distance(i,j,k,ispec) = dist_sq
          is_updated = .true.

          ! checks closest GLL point
          iglob = ibool1(i_selected,j_selected,k_selected,ispec_selected,rank_selected)
          x_found = x1(iglob,rank_selected)
          y_found = y1(iglob,rank_selected)
          z_found = z1(iglob,rank_selected)

          ! checks distance
          if (DO_WARNING) then
            dist_min = sqrt((x_found-x_target)**2 + (y_found-y_target)**2 + (z_found-z_target)**2)
            if (dist_min > 2.0 * typical_size) then
              print *,'Warning: rank ',myrank,' - large dist_min: ',dist_min * R_PLANET_KM,'(km)', &
                     'element size:',typical_size * R_PLANET_KM
              print *,'target location:',xyz_target(:)
              print *,'target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_PLANET_KM,'(km)'
              print *,'gll location   :',x_found,y_found,z_found
              print *,'gll radius     :',sqrt(x_found**2 + y_found**2 + z_found**2) * R_PLANET_KM,'(km)'
              print *,'distance gll:',dist_min * R_PLANET_KM,'(km)'
              ! debug
              !stop 'Error GLL model value invalid'
            endif
            ! debug
            !if (myrank == 0 .and. iglob < 100) &
            !  print *,'dist_min GLL point: ',dist_min * R_PLANET_KM,'(km)',typical_size * R_PLANET_KM
          endif

          ! new model values
          do iker = 1,nparams
            if (USE_FALLBACK) then
              call interpolate_limited(xi,eta,gamma,ispec_selected, &
                                       nspec_old,model1(:,:,:,ispec_selected,iker,rank_selected), &
                                       val,xigll,yigll,zigll, &
                                       i_selected,j_selected,k_selected)
            else
              call interpolate_element_value(xi,eta,gamma,ispec_selected, &
                                             nspec_old,model1(:,:,:,ispec_selected,iker,rank_selected), &
                                             val,xigll,yigll,zigll)
            endif

            ! sets new model value
            model2(i,j,k,ispec,iker) = val
          enddo

          ! checks model value differences
          do iker = 1,nparams
            val_initial = model1(i_selected,j_selected,k_selected,ispec_selected,iker,rank_selected)
            val = model2(i,j,k,ispec,iker)

            ! statistics
            if (abs(val - val_initial ) > abs(model_maxdiff(iker))) model_maxdiff(iker) = val - val_initial

            ! checks model difference
            if (DO_WARNING) then
              ! note: warns for top elements, probably due to crustal structure
              if (abs(val - val_initial ) > abs( 0.2 * val_initial )) then
                print *,'Warning: model ',iker,' value:',val,'is very different from initial value ',val_initial
                print *,'  rank ',myrank,' - dist_min: ',dist_min * R_PLANET_KM,'(km)'
                print *,'  element',ispec,'selected ispec:',ispec_selected,'in rank:',rank_selected,'iglob_min:',iglob_min
                print *,'  typical element size:',typical_size * 0.5 * R_PLANET_KM
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
            !  print *,'new model ',iker,': value ',val,'initial ',val_initial,'diff ',(val - val_initial)/val_initial*100.0,'(%)'
          enddo
        endif

      enddo
    enddo
  enddo

  contains

    subroutine check_point_result()

    implicit none

    ! checks valid iglob
    if (iglob_min < 1 .or. iglob_min > nspec_old * nproc_chunk1) then
      print *,'Error iglob_min :',iglob_min
      print *,'nspec / nproc :',nspec_old,nproc_chunk1
      stop 'Error invalid iglob_min index'
    endif

    ! checks valid rank
    if (rank_selected < 0 .or. rank_selected >= nproc_chunk1) then
      print *,'Error rank:',myrank,'invalid selected rank ',rank_selected,'for element',ispec
      print *,'target location:',xyz_target(:)
      stop 'Error specifying closest rank for element'
    endif

    ! checks valid ispec
    if (ispec_selected < 1 .or. ispec_selected > nspec_old) then
      print *,'Error rank:',myrank,'invalid selected ispec ',ispec_selected,'for element',ispec
      print *,'rank_selected:',rank_selected,'iglob_min:',iglob_min,'nspec_old:',nspec_old
      print *,'target location:',xyz_target(:)
      print *,'dist_min: ',dist_min * R_PLANET_KM,'(km)'
      stop 'Error specifying closest ispec element'
    endif

    ! checks minimum distance within a typical element size
    if (DO_WARNING) then
      if (dist_min > 2.0 * typical_size) then
        print *,'Warning: rank ',myrank,' - large dist_min: ',dist_min * R_PLANET_KM,'(km)', &
               'element size:',typical_size * R_PLANET_KM
        print *,'element',ispec,'selected ispec:',ispec_selected,'in rank:',rank_selected,'iglob_min:',iglob_min
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

  end subroutine get_model_values_kdtree

!
!------------------------------------------------------------------------------
!

  subroutine locate(x_target,y_target,z_target, &
                    xi_target,eta_target,gamma_target, &
                    ispec_selected,rank_selected, &
                    nspec,nglob,nproc, &
                    ibool,xstore,ystore,zstore, &
                    anchor_iax,anchor_iay,anchor_iaz, &
                    xigll,yigll,zigll,typical_size, &
                    i_selected,j_selected,k_selected, &
                    dist_sq)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,MIDX,MIDY,MIDZ,HUGEVAL

  implicit none

  ! point location
  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target,z_target

  ! interpolated point location
  integer,intent(out) :: ispec_selected,rank_selected
  double precision,intent(out) :: xi_target,eta_target,gamma_target

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec,nglob,nproc

  ! arrays containing coordinates of the points
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec,0:nproc-1),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(nglob,0:nproc-1),intent(in) :: xstore,ystore,zstore

  ! indices of the control points of the surface element
  integer,intent(in) :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size

  integer, intent(out) :: i_selected,j_selected,k_selected

  real(kind=CUSTOM_REAL),intent(inout) :: dist_sq

  ! local parameters
  integer :: i,j,k,ispec,iglob,rank
  double precision :: dist,distmin,typical_size_sq

  ! set distance to huge initial value
  distmin = HUGEVAL
  typical_size_sq = typical_size * typical_size

  ispec_selected = 0
  rank_selected = -1

  ! finds closest point
  do rank = 0,nproc-1
    do ispec = 1,nspec
      ! distance to cell center
      ! midpoint
      iglob = ibool(MIDX,MIDY,MIDZ,ispec,rank)
      ! distance to midpoint
      dist =   (x_target - dble(xstore(iglob,rank)))**2 &
              +(y_target - dble(ystore(iglob,rank)))**2 &
              +(z_target - dble(zstore(iglob,rank)))**2

      if (dist > typical_size_sq) cycle

      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do k = 2,NGLLZ-1
        do j = 2,NGLLY-1
          do i = 2,NGLLX-1
            iglob = ibool(i,j,k,ispec,rank)
            dist = (x_target - dble(xstore(iglob,rank)))**2 &
                  +(y_target - dble(ystore(iglob,rank)))**2 &
                  +(z_target - dble(zstore(iglob,rank)))**2

            ! keep this point if it is closer to the receiver
            if (dist < distmin) then
              distmin = dist
              rank_selected = rank
              ispec_selected = ispec
            endif
          enddo
        enddo
      enddo
    enddo ! ispec
  enddo ! rank

  ! checks if closest point found
  if (ispec_selected < 1 ) stop 'Error locating closest ispec element'
  if (rank_selected < 0 ) stop 'Error locating closest rank for element'

  ! find the best (xi,eta,gamma)
  call locate_single(x_target,y_target,z_target, &
                     xi_target,eta_target,gamma_target, &
                     ispec_selected, &
                     nspec,nglob, &
                     ibool(:,:,:,:,rank_selected),xstore(:,rank_selected),ystore(:,rank_selected),zstore(:,rank_selected), &
                     anchor_iax,anchor_iay,anchor_iaz, &
                     xigll,yigll,zigll,typical_size, &
                     i_selected,j_selected,k_selected, &
                     dist_sq)

  end subroutine locate


!
!------------------------------------------------------------------------------
!

  subroutine locate_single(x_target,y_target,z_target, &
                           xi_target,eta_target,gamma_target, &
                           ispec_selected, &
                           nspec,nglob,ibool,xstore,ystore,zstore, &
                           anchor_iax,anchor_iay,anchor_iaz, &
                           xigll,yigll,zigll,typical_size, &
                           i_selected,j_selected,k_selected, &
                           dist_sq)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,HUGEVAL,NUM_ITER, &
    MIDX,MIDY,MIDZ,myrank
  use shared_parameters, only: R_PLANET_KM

  implicit none

  ! point location
  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target,z_target

  integer,intent(inout) :: ispec_selected

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

  real(kind=CUSTOM_REAL),intent(inout) :: dist_sq

  ! local parameters
  ! use integer array to store values
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: i,j,k,iglob
  double precision :: dist
  double precision :: xi,eta,gamma
  double precision :: x,y,z
  double precision :: final_distance
  double precision :: distmin
  double precision :: typical_size_sq

  ! neighbor search
  integer :: ispec,ispec_ref
  integer :: iglob_corner(8)
  integer :: iglob_edge(2,4),iedge
  integer :: ix_n,iy_n,iz_n
  double precision :: xi_n,eta_n,gamma_n
  double precision :: dist_n,distmin_n
  integer :: i_neighbor

  !------------------------------------------------------

  ! exact position search
  logical,parameter :: DO_REFINE_LOCATION = .true.

  ! neighbor search location improvement
  ! slows downs search quite a bit, for example from 7.9 s -> 86.9 s
  ! thus, only use if really needed
  logical,parameter :: DO_ADJACENT_SEARCH = .false.

  ! debug warning
  logical,parameter :: DO_WARNING = .false.

  !------------------------------------------------------

  ! set distance to huge initial value
  distmin = HUGEVAL
  ix_initial_guess = 0
  iy_initial_guess = 0
  iz_initial_guess = 0

  typical_size_sq = typical_size * typical_size

  ! finds closest interior GLL point
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool(i,j,k,ispec_selected)

        dist = (x_target - xstore(iglob))**2 &
              +(y_target - ystore(iglob))**2 &
              +(z_target - zstore(iglob))**2

        ! keep this point if it is closer to the receiver
        if (dist < distmin) then
          distmin = dist
          ix_initial_guess = i
          iy_initial_guess = j
          iz_initial_guess = k
        endif
      enddo
    enddo
  enddo

  ! checks
  if (ix_initial_guess == 0 .or. iy_initial_guess == 0 .or. iz_initial_guess == 0 ) &
    stop 'Error locate_single: no initial guess'

  ! debug
  !print *,'distmin = ',sngl(sqrt(distmin) * R_PLANET_KM),'(km)'

  ! find the best (xi,eta)
  ! use initial guess in xi, eta and gamma from closest point found
  xi = xigll(ix_initial_guess)
  eta = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  ! iterate to solve the non linear system to find the exact position within the element
  if (DO_REFINE_LOCATION) then

    ! only if element is close enough
    if (distmin < 6.d0 * typical_size_sq) then

      ! finds refined position within selected element and rank slice
      call find_local_coordinates_refined(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                          ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                          nspec,nglob,ibool,xstore,ystore,zstore, &
                                          anchor_iax,anchor_iay,anchor_iaz)

      ! found interpolated position
      ! compute final distance between asked and found
      final_distance = (x_target-x)**2 + (y_target-y)**2 + (z_target-z)**2

      ! debug
      !if (final_distance > 1.e-2 ) &
      !  print *,'debug: final distance = ',sngl(final_distance),'distmin',distmin,'xi/eta/gamma',xi,eta,gamma

      ! checks if location improved
      if (final_distance <= distmin) then
        ! updates position distance
        distmin = final_distance
      else
        ! uses initial guess
        xi = xigll(ix_initial_guess)
        eta = yigll(iy_initial_guess)
        gamma = zigll(iz_initial_guess)
      endif

      ! loops over neighbors and tries to find better location
      if (DO_ADJACENT_SEARCH) then
        ! checks if position lies on an element boundary
        if (abs(xi) > 1.099d0 .or. abs(eta) > 1.099d0 .or. abs(gamma) > 1.099d0) then
          ! searches for better position in neighboring elements
          ispec_ref = ispec_selected

          ! takes corner indexes of elements left/right/front/back edges of reference one
          ! we won't take top/bottom edges, as they might be below/above a discontinuity and provide wrong interpolation values.
          ! common edges
          iglob_edge(1,1) = ibool(1,1,1,ispec_ref)
          iglob_edge(2,1) = ibool(1,1,NGLLZ,ispec_ref)

          iglob_edge(1,2) = ibool(NGLLX,1,1,ispec_ref)
          iglob_edge(2,2) = ibool(NGLLX,1,NGLLZ,ispec_ref)

          iglob_edge(1,3) = ibool(1,NGLLY,1,ispec_ref)
          iglob_edge(2,3) = ibool(1,NGLLY,NGLLZ,ispec_ref)

          iglob_edge(1,4) = ibool(NGLLX,NGLLY,1,ispec_ref)
          iglob_edge(2,4) = ibool(NGLLX,NGLLY,NGLLZ,ispec_ref)

          do ispec = 1,nspec
            ! skip reference element
            if (ispec == ispec_ref) cycle

            ! exclude elements that are too far from target
            if (((x_target - xstore(ibool(MIDX,MIDY,MIDZ,ispec)))**2 &
                +(y_target - ystore(ibool(MIDX,MIDY,MIDZ,ispec)))**2 &
                +(z_target - zstore(ibool(MIDX,MIDY,MIDZ,ispec)))**2) > 6.d0 * typical_size_sq) cycle

            ! the eight corners of the current element
            iglob_corner(1) = ibool(1,1,1,ispec)
            iglob_corner(2) = ibool(NGLLX,1,1,ispec)
            iglob_corner(3) = ibool(NGLLX,NGLLY,1,ispec)
            iglob_corner(4) = ibool(1,NGLLY,1,ispec)
            iglob_corner(5) = ibool(1,1,NGLLZ,ispec)
            iglob_corner(6) = ibool(NGLLX,1,NGLLZ,ispec)
            iglob_corner(7) = ibool(NGLLX,NGLLY,NGLLZ,ispec)
            iglob_corner(8) = ibool(1,NGLLY,NGLLZ,ispec)
            i_neighbor = 0

            do iedge = 1,4
              ! checks if corner also has reference element
              if (any(iglob_corner(:) == iglob_edge(1,iedge)) .and. &
                  any(iglob_corner(:) == iglob_edge(2,iedge))) then
                i_neighbor = iedge
                exit
              endif
            enddo

            ! search position in neighbor
            if (i_neighbor > 0) then
              ! starts refinement
              distmin_n = HUGEVAL
              ix_n = MIDX
              iy_n = MIDY
              iz_n = MIDZ
              ! finds closest interior GLL point
              do k = 1,NGLLZ
                do j = 1,NGLLY
                  do i = 1,NGLLX
                    ! keep this point if it is closer to the target
                    dist_n =  (x_target - xstore(ibool(i,j,k,ispec)))**2 &
                             +(y_target - ystore(ibool(i,j,k,ispec)))**2 &
                             +(z_target - zstore(ibool(i,j,k,ispec)))**2
                    if ( dist_n < distmin_n) then
                      distmin_n = dist_n
                      ix_n = i
                      iy_n = j
                      iz_n = k
                    endif
                  enddo
                enddo
              enddo
              ! sets initial xi,eta,gamma
              xi_n = xigll(ix_n)
              eta_n = yigll(ix_n)
              gamma_n = zigll(ix_n)

              ! finds refined position within selected element and rank slice
              call find_local_coordinates_refined(x_target,y_target,z_target,xi_n,eta_n,gamma_n,x,y,z, &
                                                  ispec,ix_n,iy_n,iz_n, &
                                                  nspec,nglob,ibool,xstore,ystore,zstore, &
                                                  anchor_iax,anchor_iay,anchor_iaz)

              ! compute final distance between asked and found
              final_distance = (x_target-x)**2 + (y_target-y)**2 + (z_target-z)**2

              ! checks if position lies inside element (which usually means that located position is accurate)
              if (abs(xi) < 1.099d0 .and. abs(eta) < 1.099d0 .and. abs(gamma) < 1.099d0) then
                ! checks if location improved
                if (final_distance < distmin) then
                  ! debug
                  !print *,'debug: neighbor',i_neighbor,' final distance = ',sngl(final_distance),'distmin',sngl(distmin), &
                  !        'xi/eta/gamma',sngl(xi_n),sngl(eta_n),sngl(gamma_n)

                  ! updates position distance
                  distmin = final_distance
                  xi = xi_n
                  eta = eta_n
                  gamma = gamma_n
                  ispec_selected = ispec
                  ix_initial_guess = ix_n
                  iy_initial_guess = iy_n
                  iz_initial_guess = iz_n

                  ! okay position, no need to look further
                  exit
                endif
              endif
            endif ! is_neighbor
          enddo ! ispec
        endif

      endif ! DO_ADJACENT_SEARCH

      ! checks valid distance
      if (distmin == HUGEVAL .or. distmin /= distmin) then
        print *,'Error: locating location ',x_target,y_target,z_target,'found',x,y,z,'dist',distmin,final_distance
        stop 'Error locating location'
      endif

    endif
  endif ! DO_REFINE_LOCATION

  ! return xi,eta,gamma of point found
  xi_target = xi
  eta_target = eta
  gamma_target = gamma

  i_selected = ix_initial_guess
  j_selected = iy_initial_guess
  k_selected = iz_initial_guess

  dist_sq = real(distmin,kind=CUSTOM_REAL)

  ! add warning if estimate is poor
  ! (usually means receiver outside the mesh given by the user)
  if (DO_WARNING) then
    if (dist_sq > typical_size_sq) then
      print *, '*****************************************************************'
      print *, '***** WARNING: location estimate is poor                    *****'
      print *, '*****************************************************************'
      print *, 'closest estimate found: ',sngl(sqrt(dist_sq)*R_PLANET_KM),'km away',' - not within:',typical_size * R_PLANET_KM
      print *, ' myrank ',myrank,' in element ',ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess
      print *, ' at xi,eta,gamma coordinates = ',xi,eta,gamma
      print *, ' at radius ',sqrt(x**2 + y**2 + z**2) * R_PLANET_KM,'(km)'
      print *, ' initial distance :',sqrt(distmin)*R_PLANET_KM,'(km)'
    endif
  endif

  end subroutine locate_single


!
!-------------------------------------------------------------------------------------------------
!

  subroutine find_local_coordinates_refined(x_target,y_target,z_target,xi,eta,gamma,x,y,z, &
                                            ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess, &
                                            nspec,nglob,ibool,xstore,ystore,zstore, &
                                            anchor_iax,anchor_iay,anchor_iaz)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,HUGEVAL,NUM_ITER
  use shared_parameters, only: R_PLANET_KM

  implicit none

  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target,z_target

  double precision,intent(inout) :: xi,eta,gamma
  double precision,intent(out) :: x,y,z

  integer,intent(in) :: ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess

  ! arrays containing coordinates of the points
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(nglob),intent(in) :: xstore,ystore,zstore

  ! indices of the control points of the surface element
  integer,intent(in) :: anchor_iax(NGNOD),anchor_iay(NGNOD),anchor_iaz(NGNOD)

  ! local parameters
  integer :: ia,iter_loop
  integer :: iglob

  ! coordinates of the control points of the surface element
  double precision :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  double precision :: dx,dy,dz,dx_min,dy_min,dz_min,d_min_sq
  double precision :: dxi,deta,dgamma
  double precision :: xix,xiy,xiz
  double precision :: etax,etay,etaz
  double precision :: gammax,gammay,gammaz
  ! debug
  double precision :: dist

  ! see original refinement version in routine find_local_coordinates() in file locate_point.f90

  ! define coordinates of the control points of the element
  do ia = 1,NGNOD
    iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec_selected)
    xelm(ia) = dble(xstore(iglob))
    yelm(ia) = dble(ystore(iglob))
    zelm(ia) = dble(zstore(iglob))
  enddo

  ! starts with initial guess in xi,eta,gamma
  d_min_sq = HUGEVAL
  dx_min = HUGEVAL
  dy_min = HUGEVAL
  dz_min = HUGEVAL

  ! iterate to solve the non linear system
  do iter_loop = 1,NUM_ITER

    ! recompute Jacobian for the new point
    call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                            xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

    ! debug
    if (.false.) then
      print *,'debug: jacobian error in locate_single(): '
      print *,'debug: jacobian error i,j,k,ispec :',ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected
      print *,'debug: jacobian error iter_loop   :',iter_loop
      dist = (x_target-x)**2 + (y_target-y)**2 + (z_target-z)**2
      print *,'debug: jacobian error dist        :',dist,sqrt(dist)*R_PLANET_KM,'(km)'
      ! uses initial guess again
      !xi = xigll(ix_initial_guess)
      !eta = yigll(iy_initial_guess)
      !gamma = zigll(iz_initial_guess)
      ! uses previous guess
      !xi = xi - dxi
      !eta = eta - deta
      !gamma = gamma - dgamma
      ! exits loop
      !exit
      !stop 'Error recomputing jacobian'
    endif

    ! compute distance to target location
    dx = - (x - x_target)
    dy = - (y - y_target)
    dz = - (z - z_target)

    !debug
    !print *,'  iter ',iter_loop,'dx',sngl(dx),sngl(dx_min),'dy',sngl(dy),sngl(dy_min),'dz',sngl(dz),sngl(dz_min),d_min_sq

    ! compute increments
    if ((dx**2 + dy**2 + dz**2) < d_min_sq) then
      d_min_sq = dx**2 + dy**2 + dz**2
      dx_min = dx
      dy_min = dy
      dz_min = dz

      dxi = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz
    else
      ! new position is worse than old one, no change necessary
      dxi = 0.d0
      deta = 0.d0
      dgamma = 0.d0
    endif

    ! decreases step length if step is large
    if ((dxi*dxi + deta*deta + dgamma*dgamma) > 1.0d0) then
      dxi = dxi * 0.33333333333d0
      deta = deta * 0.33333333333d0
      dgamma = dgamma * 0.33333333333d0
    endif
    ! alternative: impose limit on increments (seems to result in slightly less accurate locations)
    !if (abs(dxi) > 0.3d0 ) dxi = sign(1.0d0,dxi)*0.3d0
    !if (abs(deta) > 0.3d0 ) deta = sign(1.0d0,deta)*0.3d0
    !if (abs(dgamma) > 0.3d0 ) dgamma = sign(1.0d0,dgamma)*0.3d0

    !debug
    !print *,'  dxi/..',(dxi**2 + deta**2 + dgamma**2),dxi,deta,dgamma

    ! update values
    xi = xi + dxi
    eta = eta + deta
    gamma = gamma + dgamma

    ! impose that we stay in that element
    ! (useful if user gives a receiver outside the mesh for instance)
    ! we can go slightly outside the [1,1] segment since with finite elements
    ! the polynomial solution is defined everywhere
    ! can be useful for convergence of iterative scheme with distorted elements
    if (xi > 1.10d0) xi = 1.10d0
    if (xi < -1.10d0) xi = -1.10d0
    if (eta > 1.10d0) eta = 1.10d0
    if (eta < -1.10d0) eta = -1.10d0
    if (gamma > 1.10d0) gamma = 1.10d0
    if (gamma < -1.10d0) gamma = -1.10d0

  ! end of non linear iterations
  enddo

  ! compute final coordinates x,y,z (and xix,xiy,..) of point found
  call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

  end subroutine find_local_coordinates_refined

