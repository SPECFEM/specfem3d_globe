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
! usage: ./xinterpolate_model vsv MODEL_M10/
!
! note on mid-point search option:
!  - mid-point-search == 1: looking for mid-points only is a good approach when changing number of processes (NPROC) only,
!                           while keeping the mesh resolution (NEX) the same
!                           (by default set to .true.)
!  - mid-point-search == 0: searching for each single gll point is a good approach when changing resolution (NEX) of meshes;
!                           in general, interpolation suffers and might lead to differences at internal interfaces (e.g. 410)
!
!------------------------------------------------------------------------------


  program interpolate_model

  use constants,only: SIZE_INTEGER, &
    TWO_PI,R_UNIT_SPHERE, &
    NGNOD,MIDX,MIDY,MIDZ,R_EARTH_KM, &
    IFLAG_CRUST,IFLAG_80_MOHO,IFLAG_220_80,IFLAG_670_220,IFLAG_MANTLE_NORMAL

  use postprocess_par,only: &
    CUSTOM_REAL,NGLLX,NGLLY,NGLLZ, &
    GAUSSALPHA,GAUSSBETA,R_EARTH_KM, &
    IIN,IOUT,MAX_STRING_LEN, &
    NCHUNKS_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,NPROCTOT_VAL,NEX_XI_VAL, &
    NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE

  use kdtree_search, only: kdtree_setup,kdtree_set_verbose,kdtree_delete,kdtree_find_nearest_neighbor, &
    kdtree_num_nodes,kdtree_nodes_location,kdtree_nodes_index

  implicit none

  !-------------------------------------------------------------------
  ! USER PARAMETERS

  ! model type
  ! isotropic model parameters (vp,vs,rho) or
  ! transversely isotropic model parameters (vpv,vph,vsv,vsh,eta,rho)
  ! defaults: TI models
  logical,parameter :: USE_TRANSVERSE_ISOTROPY = .true.

  ! shear attenuation
  logical,parameter :: USE_ATTENUATION_Q = .false.

  ! brute-force search for closest element
  ! by default set to .false., thus a tree search for initial guess element is used
  logical,parameter :: DO_BRUTE_FORCE_SEARCH = .false.

  ! kd-tree setup
  ! uses all internal gll points for search tree
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
  integer :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! typical element size in old mesh used as search radius
  double precision :: typical_size

  ! old, initial (source) mesh
  ! check constants defined for old, first mesh from file: 'values_from_mesher.h'
  ! NPROC_XI and NPROC_ETA
  integer :: nproc_eta_old,nproc_xi_old
  ! NSPEC_CRUST_MANTLE
  integer :: nspec_max_old
  ! NGLOB_CRUST_MANTLE
  integer :: nglob_max_old
  ! combined slices, old mesh
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: x1,y1,z1
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:,:),allocatable :: model1
  integer,dimension(:,:,:,:,:),allocatable :: ibool1
  integer,dimension(:,:),allocatable :: idoubling1
  integer, dimension(:,:,:),allocatable :: addressing1

  ! new mesh
  ! single slice, target mesh
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: x2, y2, z2
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:),allocatable :: model2
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool2
  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling2
  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing2

  ! input arguments
  character(len=MAX_STRING_LEN) :: arg
  character(len=MAX_STRING_LEN) :: dir_topo1
  character(len=MAX_STRING_LEN) :: dir_topo2
  character(len=MAX_STRING_LEN) :: input_model_dir,output_model_dir

  ! kdtree search
  integer :: want_midpoint
  logical :: USE_MIDPOINT_SEARCH

  integer :: i,j,k,iglob,ispec,ier,iker
  integer :: nspec, nglob, rank
  integer :: nproc_chunk1
  integer :: ichunk,ichunk_selected
  integer :: iproc_eta,iproc_xi,iprocnum
  integer :: iproc_eta_selected,iproc_xi_selected
  integer :: ilayer
  integer :: total_nspec

  ! model
  integer :: nparams
  character(len=16) :: fname(7)
#ifdef ADIOS_INPUT
  character(len=16) :: model_name(7)
#endif
  character(len=MAX_STRING_LEN) :: m_file
  character(len=MAX_STRING_LEN) :: solver_file

  ! mpi parameters
  integer :: sizeprocs,myrank

  ! nodes search
  integer :: inodes
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: model_maxdiff
  real(kind=CUSTOM_REAL) :: val

  ! starts mpi
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

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


  ! checks program arguments
  if (myrank == 0) then
    do i = 1,4
      call get_command_argument(i,arg)
      ! usage info
      if (len_trim(arg) == 0) then
        if (myrank == 0) then
          print *,' '
          print *,' Usage: xinterpolate_model old-topo-dir/ old-model-dir/ new-topo-dir/ model-output-dir/ (midpoint-search)'
          print *,' '
          print *,' with'
          print *,'   old-topo-dir/     - old mesh directory with topology files (e.g. proc***_solver_data.bin)'
          print *,'   old-model-dir/    - directoy which holds old model files (e.g. proc***_vpv.bin)'
          print *,'   new-topo-dir/     - new mesh directory with topology files (e.g. proc***_solver_data.bin)'
          print *,'   model-output-dir/ - output directory for interpolated model on new mesh'
          print *,'   (optional) midpoint-search = 0  - uses every single gll point for search of closest element'
          print *,'                              = 1  - uses midpoints for search of closest element (default)'
          print *,' '
        endif
        call synchronize_all()
        stop ' Reenter command line options'
      endif
    enddo
  endif
  call synchronize_all()

  ! reads input arguments
  want_midpoint = 1
  do i = 1, 5
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
          if (myrank == 0) print*,'Error reading in midpoint-search value, please check your arguments...'
          stop ' Reenter command line options'
        endif
      endif
    endif
  enddo

  ! kdtree search:
  ! searches closest element using mid-points only, rather than for every single gll point
  ! note: looking for mid-points only is a good approach when changing number of processes (NPROC)
  !       while keeping the mesh resolution (NEX) the same;
  !       searching for each single gll point is a good approach when changing resolution (NEX) of meshes;
  !       in general, interpolation suffers and might lead to differences at internal interfaces (e.g. 410);
  if (want_midpoint == 1) then
    USE_MIDPOINT_SEARCH = .true.
  else
    USE_MIDPOINT_SEARCH = .false.
  endif

  ! console output
  if (myrank == 0) then
    print*
    print*,'model interpolation:'
    print*
  endif

#ifdef ADIOS_INPUT
  stop 'safety stop: ADIOS support not fully implemented yet...'

  ! initializes ADIOS
  if (myrank == 0) then
    print *, 'initializing ADIOS...'
    print *, ' '
  endif
  call adios_setup()
#endif

  ! defines model parameters
  if (USE_TRANSVERSE_ISOTROPY) then
    ! transversly isotropic (TI) model
    nparams = 6
#ifdef ADIOS_INPUT
    model_name(1:6) = (/character(len=16) :: "reg1/vpv","reg1/vph","reg1/vsv","reg1/vsh","reg1/eta","reg1/rho"/)
#endif
    fname(1:6) = (/character(len=16) :: "vpv","vph","vsv","vsh","eta","rho"/)
  else
    ! isotropic model
    nparams = 3
    ! note: adds space at endings to equal number of characters
    !       to avoid compiler error: "Different shape for array assignment.."
#ifdef ADIOS_INPUT
    model_name(1:3) = (/character(len=16) :: "reg1/vp ","reg1/vs ","reg1/rho"/)
#endif
    fname(1:3) = (/character(len=16) :: "vp ","vs ","rho"/)
  endif
  ! adds shear attenuation
  if (USE_ATTENUATION_Q) then
    nparams = nparams + 1
#ifdef ADIOS_INPUT
    model_name(nparams) = "reg1/qmu"
#endif
    fname(nparams) = "qmu"
  endif

  ! master process gets old, source mesh dimension
  if (myrank == 0) then
    ! gets nspec/nglob
    write(solver_file,'(a,i6.6,a)') trim(dir_topo1)//'proc',myrank,'_reg1_'//'solver_data.bin'
    open(IIN,file=trim(solver_file),status='old',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening file: ',trim(solver_file)
      stop 'Error opening old solver_data.bin file, please check arguments...'
    endif
    read(IIN) nspec_max_old
    read(IIN) nglob_max_old
    close(IIN)

    ! gets number of processes from old mesh (1 file per process)
    rank = 0
    do while (ier == 0)
      write(solver_file,'(a,i6.6,a)') trim(dir_topo1)//'proc',rank,'_reg1_'//'solver_data.bin'
      open(IIN,file=trim(solver_file),status='old',form='unformatted',iostat=ier)
      if (ier == 0) then
        rank = rank + 1
        close(IIN)
      endif
    enddo
    ! assumes same number of chunks and nproc_eta = nproc_xi
    nproc_eta_old = int (sqrt ( dble(rank) / NCHUNKS_VAL ))
    nproc_xi_old = int (sqrt ( dble(rank) / NCHUNKS_VAL ))
  endif
  ! broadcasts to all other processes (assumes all slices have equal nspec/nglob values)
  call bcast_all_singlei(nspec_max_old)
  call bcast_all_singlei(nglob_max_old)
  call bcast_all_singlei(nproc_eta_old)
  call bcast_all_singlei(nproc_xi_old)

  ! sets old nproc_xi (assumes equal nproc_xi/nproc_eta)
  nproc_xi_old = nproc_eta_old

  ! console output
  if (myrank == 0) then
    print*,'source mesh:  '
    print*,'  processors = ',nproc_eta_old * nproc_xi_old * NCHUNKS_VAL
    print*,'  nproc_eta / nproc_xi = ',nproc_eta_old,nproc_xi_old
    print*,'  nspec      = ',nspec_max_old
    print*,'  nglob      = ',nglob_max_old
    print*
    print*,'target mesh:  '
    print*,'  processors = ',NPROCTOT_VAL
    print*,'  nproc_eta / nproc_xi = ',NPROC_ETA_VAL,NPROC_XI_VAL
    print*,'  nex        = ',NEX_XI_VAL
    print*,'  nspec      = ',NSPEC_CRUST_MANTLE
    print*,'  nglob      = ',NGLOB_CRUST_MANTLE
    print*
    if (USE_TRANSVERSE_ISOTROPY) then
      print*,'model parameters:',nparams,' - transversely isotropic model'
    else
      print*,'model parameters:',nparams,' - isotropic model'
    endif
    if (USE_ATTENUATION_Q) then
      print*,'  includes qmu model parameter'
    endif
    print*,'  ( ',(trim(fname(i))//" ",i=1,nparams),')'
    print*
    print*,'input model  directory: ',trim(input_model_dir)
    print*,'output model directory: ',trim(output_model_dir)
    print*
    print*,'array size:'
    print*,'  ibool1   = ',NGLLX*NGLLY*NGLLZ*nspec_max_old*nproc_eta_old*nproc_xi_old*dble(SIZE_INTEGER)/1024./1024.,'MB'
    print*,'  x1,y1,z1 = ',nglob_max_old*nproc_eta_old*nproc_xi_old*dble(CUSTOM_REAL)/1024./1024.,'MB'
    print*
    print*,'  model1   = ',NGLLX*NGLLY*NGLLZ*nspec_max_old*nparams*nproc_eta_old*nproc_xi_old*dble(CUSTOM_REAL)/1024./1024.,'MB'
    print*,'  model2   = ',NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE*nparams*dble(CUSTOM_REAL)/1024./1024.,'MB'
    print*
    print*,'total mpi processes: ',sizeprocs
    print*
    if (DO_BRUTE_FORCE_SEARCH) then
      print*,'location search by : brute-force approach'
    else
      print*,'location search by : kd-tree search'
      if (USE_MIDPOINT_SEARCH) then
        print*,'location search by : uses midpoints of elements only'
      else
        print*,'  uses internal gll points'
      endif
      if (DO_SEPARATION_410_650) then
        print*,'  uses element separation for 410-km/650-km discontinuity'
      endif
      if (DO_SEPARATION_TOPO) then
        print*,'  uses element separation for surface (moho) discontinuity'
      endif
      if (USE_FALLBACK) then
        print*,'  uses fall-back to model value of closest point in case of large differences'
      endif
    endif
    print*
  endif
  call synchronize_all()

  ! checks
  if (sizeprocs /= NPROCTOT_VAL) stop 'Error target mesh processors not equal to current total mpi processes'

  ! checks temporary file creation, to see if we could write out new model
  if (myrank == 0) then
    write(m_file,'(a,i6.6,a)') trim(output_model_dir)// '/proc',myrank,'_reg1_'//trim(fname(1))//'.tmp'
    open(IOUT,file=trim(m_file),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening file: ',trim(m_file)
      stop 'Error opening new output model file, please check if output directory exists...'
    endif
    close(IOUT,status='delete')
  endif
  call synchronize_all()

  ! total number of processes per chunk (source mesh)
  nproc_chunk1 = nproc_eta_old * nproc_xi_old

  ! collected mesh arrays for a single chunk
  allocate( x1(nglob_max_old,0:nproc_chunk1-1), &
            y1(nglob_max_old,0:nproc_chunk1-1), &
            z1(nglob_max_old,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating locations'
  allocate( ibool1(NGLLX,NGLLY,NGLLZ,nspec_max_old,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating ibool1'
  allocate( idoubling1(nspec_max_old,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating idoubling1'
  allocate( addressing1(NCHUNKS_VAL,0:nproc_xi_old-1,0:nproc_eta_old-1),stat=ier )
  if (ier /= 0) stop 'Error allocating addressing1'

  ! model files
  allocate( model1(NGLLX,NGLLY,NGLLZ,nspec_max_old,nparams,0:nproc_chunk1-1),stat=ier )
  if (ier /= 0) stop 'Error allocating initial model1'
  allocate( model2(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE,nparams),stat=ier )
  if (ier /= 0) stop 'Error allocating target model2'
  ! statistics
  allocate( model_maxdiff(nparams),stat=ier)
  if (ier /= 0) stop 'Error allocating model_maxdiff'

  ! synchronizes
  call synchronize_all()

  ! GLL points
  call zwgljd(xigll,wxgll,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ,GAUSSALPHA,GAUSSBETA)

  ! define topology of the control element
  call hex_nodes(iaddx,iaddy,iaddr)

  ! source mesh
  ! compute typical size of elements at the surface
  typical_size = TWO_PI * R_UNIT_SPHERE / ( 4.0 * NEX_XI_VAL )
  ! use 2 times the distance as a criterion for point detection
  typical_size = 2.0 * typical_size

  ! creates addressing
  addressing1(:,:,:) = 0
  addressing2(:,:,:) = 0
  do ichunk = 1,NCHUNKS_VAL
    ! old mesh addressing
    do iproc_eta = 0,nproc_eta_old-1
      do iproc_xi = 0,nproc_xi_old-1
        iprocnum = (ichunk-1)*nproc_xi_old*nproc_eta_old + iproc_eta * nproc_xi_old + iproc_xi
        addressing1(ichunk,iproc_xi,iproc_eta) = iprocnum
      enddo
    enddo

    ! new target mesh addressing
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
  do ichunk = 1,NCHUNKS_VAL
    do iproc_eta = 0, NPROC_ETA_VAL - 1
      do iproc_xi = 0, NPROC_XI_VAL - 1
        ! gets slice number
        rank = addressing2(ichunk,iproc_xi,iproc_eta)

        ! only associated mpi process continues
        if (myrank == rank) then
          ichunk_selected = ichunk
          iproc_eta_selected = iproc_eta
          iproc_xi_selected = iproc_xi
          exit
        endif
      enddo
    enddo
  enddo
  if (ichunk_selected < 1 .or. ichunk_selected > NCHUNKS_VAL ) stop 'Error selecting ichunk'
  ! debug
  !print*, 'selected chunk: ',ichunk_selected,' - eta/xi : ',iproc_eta_selected,iproc_xi_selected

  ! user output
  if (myrank == 0) then
    print*
    print*, 'loading source mesh ... '
  endif

  ! reads in model and locations of old, source mesh
  ! combines all slices for this whole chunk
  x1(:,:) = 0.0_CUSTOM_REAL
  y1(:,:) = 0.0_CUSTOM_REAL
  z1(:,:) = 0.0_CUSTOM_REAL
  ibool1(:,:,:,:,:) = 0

  iprocnum = 0
  do iproc_eta = 0, nproc_eta_old - 1
    do iproc_xi = 0, nproc_xi_old - 1
      ! gets slice number
      rank = addressing1(ichunk_selected,iproc_xi,iproc_eta)

      ! counter
      iprocnum = iprocnum + 1

      ! user output
      if (myrank == 0) then
        print*,'  slice number: ',iprocnum,' out of ',nproc_chunk1
      endif

      ! old, source mesh locations
      write(solver_file,'(a,i6.6,a)') trim(dir_topo1)//'proc',rank,'_reg1_'//'solver_data.bin'
      open(IIN,file=solver_file,status='old',form='unformatted',iostat=ier)
      if (ier /= 0) then
        print*,'Error opening file: ',trim(solver_file)
        stop 'Error opening old solver_data.bin file'
      endif
      read(IIN) nspec
      read(IIN) nglob

      ! checks dimensions
      if (nspec /= nspec_max_old .or. nglob /= nglob_max_old) then
        print*,'Error dimension of old, source mesh: solver_data nspec/nglob = ',nspec,nglob
        stop 'Error new mesh dimensions'
      endif

      read(IIN) x1(:,iprocnum-1)
      read(IIN) y1(:,iprocnum-1)
      read(IIN) z1(:,iprocnum-1)
      read(IIN) ibool1(:,:,:,:,iprocnum-1)
      read(IIN) idoubling1(:,iprocnum-1)
      close(IIN)
    enddo
  enddo
  ! user output
  if (myrank == 0) then
    print*
    print*,'  source mesh chunk read successfully'
    print*
  endif
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print*, 'loading source model ... '
    do iker = 1,nparams
      print*, '  for parameter: ',trim(fname(iker))
    enddo
    print*
  endif

  ! reads in old model files
  model1(:,:,:,:,:,:) = 0.0_CUSTOM_REAL
  iprocnum = 0
  do iproc_eta = 0, nproc_eta_old - 1
    do iproc_xi = 0, nproc_xi_old - 1
      ! gets slice number
      rank = addressing1(ichunk_selected,iproc_xi,iproc_eta)

      ! counter
      iprocnum = iprocnum + 1

      ! user output
      if (myrank == 0) then
        print*,'  slice number: ',iprocnum,' out of ',nproc_chunk1
      endif

      ! reads in model slices
      do iker = 1,nparams
        ! debug user output
        !if (myrank == 0) print *, '  for parameter: ',trim(fname(iker))
        ! opens model file
        write(m_file,'(a,i6.6,a)') trim(input_model_dir)//'proc',rank,'_reg1_'//trim(fname(iker))//'.bin'
        open(IIN,file=trim(m_file),status='old',form='unformatted',iostat=ier)
        if (ier /= 0) then
          print*,'Error opening file: ',trim(m_file)
          stop 'Error opening old model file'
        endif
        read(IIN) model1(:,:,:,:,iker,iprocnum-1)
        close(IIN)
      enddo
    enddo
  enddo
  ! user output
  if (myrank == 0) then
    print*
    print*,'  source model chunk read successfully'
    print*
  endif
  call synchronize_all()

  ! reads in the topology files of the target slices
  ! gets slice number
  rank = addressing2(ichunk_selected,iproc_xi_selected,iproc_eta_selected)

  ! only associated mpi process continues
  if (myrank /= rank) stop 'Error selected addressing rank'

  ! user output
  if (myrank == 0) then
    print*
    print*,'reading new mesh slice ... '
    print*
  endif

  ! checks new mesh locations
  write(solver_file,'(a,i6.6,a)') trim(dir_topo2)//'proc',rank,'_reg1_'//'solver_data.bin'
  open(IIN,file=solver_file,status='old',form='unformatted',iostat=ier)
  if (ier /= 0) then
    print*,'Error opening file: ',trim(solver_file)
    stop 'Error opening new solver_data.bin file'
  endif
  read(IIN) nspec
  read(IIN) nglob

  ! checks dimensions
  if (nspec /= NSPEC_CRUST_MANTLE .or. nglob /= NGLOB_CRUST_MANTLE) then
    print*,'Error dimension of new mesh: solver_data nspec/nglob = ',nspec,nglob
    stop 'Error new mesh dimensions'
  endif

  ! locations
  read(IIN) x2(:)
  read(IIN) y2(:)
  read(IIN) z2(:)
  read(IIN) ibool2(:,:,:,:)
  read(IIN) idoubling2(:)
  close(IIN)

  ! checks that layers match
  if (minval(idoubling1) /= minval(idoubling2) .or. maxval(idoubling1) /= maxval(idoubling2)) then
    print*,'Error idoubling range:'
    print*,'idoubling 1:',minval(idoubling1),maxval(idoubling1)
    print*,'idoubling 2:',minval(idoubling2),maxval(idoubling2)
    stop 'Error invalid idoubling range'
  endif
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print*,'Earth layers: ',minval(idoubling2),' to ',maxval(idoubling2)
    print*
  endif

  ! get values in elements for this new slice
  model2(:,:,:,:,:) = 0.0_CUSTOM_REAL

  ! loops over layers (crust/moho-80/80-220/220-660/660-CMB)
  total_nspec = 0
  do ilayer = minval(idoubling2),maxval(idoubling2)

    ! user output
    if (myrank == 0) then
      print*,'layer: ',ilayer,' out of ',maxval(idoubling2)
      select case (ilayer)
      case (IFLAG_CRUST)
        print*,'layer: crust'
      case (IFLAG_80_MOHO)
        print*,'layer: 80 - MOHO'
      case (IFLAG_220_80)
        print*,'layer: 220 - 80'
      case (IFLAG_670_220)
        print*,'layer: 670 - 220'
      case (IFLAG_MANTLE_NORMAL)
        print*,'layer: mantle normal'
      end select
    endif

    ! statistics
    model_maxdiff(:) = 0.0_CUSTOM_REAL

    ! builds search tree
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      ! counts total number of points in this layer in source mesh
      iprocnum = 0
      inodes = 0
      do iproc_eta = 0, nproc_eta_old - 1
        do iproc_xi = 0, nproc_xi_old - 1
          ! counter
          iprocnum = iprocnum + 1
          ! all elements
          do ispec = 1, nspec_max_old
            if (idoubling1(ispec,iprocnum-1) == ilayer) then
              if (TREE_INTERNAL_GLL_POINTS) then
                ! all internal gll points ( 2 to NGLLX-1 )
                inodes = inodes + (NGLLX-2)*(NGLLY-2)*(NGLLZ-2)
              endif
              if (TREE_MID_POINTS) then
                ! only element mid-points
                inodes = inodes + 1
              endif
            endif
          enddo
        enddo
      enddo
      ! checks
      if (inodes < 1 ) stop 'Error no search tree nodes in source mesh for this layer'
      ! checks maximum number of nodes
      k = 0
      if (TREE_INTERNAL_GLL_POINTS) k = nspec_max_old * nproc_chunk1 * (NGLLX-2)*(NGLLY-2)*(NGLLZ-2)
      if (TREE_MID_POINTS) k = nspec_max_old * nproc_chunk1 * 1
      if (inodes > k ) stop 'Error invalid number of search tree nodes in this layer'

      ! set number of tree nodes
      kdtree_num_nodes = inodes

      ! debug
      !print*,'kdtree nodes: ',kdtree_num_nodes

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
      do iproc_eta = 0, nproc_eta_old - 1
        do iproc_xi = 0, nproc_xi_old - 1

          ! counter
          iprocnum = iprocnum + 1

          ! adds tree nodes
          do ispec = 1,nspec_max_old

            ! skips elements outside of this layer
            if (idoubling1(ispec,iprocnum-1) /= ilayer ) cycle

            ! sets up tree nodes
            ! all internal gll points
            if (TREE_INTERNAL_GLL_POINTS) then
              do k = 2,NGLLZ-1
                do j = 2,NGLLY-1
                  do i = 2,NGLLX-1
                    iglob = ibool1(i,j,k,ispec,iprocnum-1)

                    ! counts nodes
                    inodes = inodes + 1
                    if (inodes > kdtree_num_nodes ) stop 'Error index inodes bigger than kdtree_num_nodes'

                    ! adds node index ( index points to same ispec for all internal gll points)
                    kdtree_nodes_index(inodes) = ispec + (iprocnum - 1) * nspec_max_old

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

              ! adds node index ( index points to same ispec for all internal gll points)
              kdtree_nodes_index(inodes) = ispec + (iprocnum - 1) * nspec_max_old

              ! adds node location
              kdtree_nodes_location(1,inodes) = x1(iglob,iprocnum-1)
              kdtree_nodes_location(2,inodes) = y1(iglob,iprocnum-1)
              kdtree_nodes_location(3,inodes) = z1(iglob,iprocnum-1)
            endif

          enddo
        enddo
      enddo
      if (inodes /= kdtree_num_nodes ) stop 'Error index inodes does not match nnodes_local'
      call synchronize_all()

      ! creates kd-tree for searching
      ! serial way
      !do i = 0,NPROCTOT_VAL-1
      !  if (myrank == i) then
      !    print*,'kd-tree setup for process: ',myrank
      !    call kdtree_setup()
      !  endif
      !  call synchronize_all()
      !enddo
      ! parallel way
      call kdtree_setup()

    endif
    call synchronize_all()

    ! user output
    if (myrank == 0) print*,'looping over elements:'
    call synchronize_all()

    ! loop over all elements (mainly those in this layer)
    do ispec = 1, nspec
      ! user output
      if (myrank == 0) then
        if (ispec == 1 .or. mod(ispec,int(0.1*nspec)) == 0 .or. ispec == nspec) then
          print*,'  ispec',ispec,' out of ',nspec
        endif
      endif

      ! skip elements out of this layer
      if (idoubling2(ispec) /= ilayer ) cycle

      ! increases element counter
      total_nspec = total_nspec + 1

      ! gets model values
      if (DO_BRUTE_FORCE_SEARCH) then
        ! brute-force search over all gll points
        call get_model_values_bruteforce(ispec,nspec,nglob,ibool2,x2,y2,z2,nparams,model2, &
                                         nspec_max_old,nglob_max_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                         iaddx,iaddy,iaddr,xigll,yigll,zigll,typical_size,myrank,model_maxdiff)
      else
        ! kdtree search
        call get_model_values_kdtree(ispec,nspec,nglob,ibool2,x2,y2,z2,nparams,model2, &
                                     nspec_max_old,nglob_max_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                     iaddx,iaddy,iaddr,xigll,yigll,zigll,typical_size,myrank,model_maxdiff, &
                                     USE_MIDPOINT_SEARCH,DO_SEPARATION_410_650,DO_SEPARATION_TOPO,USE_FALLBACK)
      endif
    enddo ! ispec
    if (myrank == 0) print*
    call synchronize_all()

    ! frees tree memory
    if (.not. DO_BRUTE_FORCE_SEARCH) then
      ! deletes tree arrays
      deallocate(kdtree_nodes_location)
      deallocate(kdtree_nodes_index)
      ! deletes search tree nodes
      call kdtree_delete()
    endif

    ! statistics
    ! user output
    if (myrank == 0) print*,'statistics:'
    do iker = 1,nparams
      call max_all_cr(model_maxdiff(iker),val)
      if (myrank == 0) print*,'  parameter ',trim(fname(iker)),': maximum difference = ',val
    enddo

    ! user output
    if (myrank == 0) print*
    call synchronize_all()
  enddo ! ilayer

  ! frees memory
  deallocate(x1,y1,z1)
  deallocate(ibool1)
  deallocate(idoubling1)
  deallocate(model1)
  deallocate(model_maxdiff)

  ! checks
  if (total_nspec /= nspec) stop 'Error invalid total number of elements after loop'

  ! user output
  if (myrank == 0) then
    print *, 'writing out new model files'
  endif

  ! writes out new model
  do iker = 1,nparams
    ! user output
    if (myrank == 0) then
      print *, '  for parameter: ',trim(fname(iker))
    endif

    write(m_file,'(a,i6.6,a)') trim(output_model_dir) // '/proc',rank,'_reg1_'//trim(fname(iker))//'.bin'
    open(IOUT,file=trim(m_file),status='unknown',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print*,'Error opening file: ',trim(m_file)
      stop 'Error opening output model file'
    endif
    write(IOUT) model2(:,:,:,:,iker)
    close(IOUT)
  enddo

  ! frees memory
  deallocate(model2)

  ! synchronizes MPI processes
  call synchronize_all()

  ! user output
  if (myrank == 0) then
    print *
    print *, 'check new model files in directory: ',trim(output_model_dir)
    print *, 'done successfully'
    print *
  endif

#ifdef ADIOS_INPUT
  call adios_finalize (myrank, ier)
#endif

  ! exiting MPI processes
  call finalize_mpi()

  end program interpolate_model

!
!------------------------------------------------------------------------------
!

  subroutine get_model_values_bruteforce(ispec,nspec,nglob,ibool2,x2,y2,z2,nparams,model2, &
                                         nspec_max_old,nglob_max_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                         iaddx,iaddy,iaddr,xigll,yigll,zigll,typical_size,myrank,model_maxdiff)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD

  implicit none

  integer,intent(in) :: ispec
  integer,intent(in) :: nparams

  ! new, target mesh
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool2
  real(kind=CUSTOM_REAL),dimension(nglob),intent(in) :: x2,y2,z2
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec,nparams),intent(inout) :: model2

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec_max_old,nglob_max_old,nproc_chunk1
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec_max_old,0:nproc_chunk1-1),intent(in) :: ibool1
  real(kind=CUSTOM_REAL),dimension(nglob_max_old,0:nproc_chunk1-1),intent(in) :: x1,y1,z1
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_max_old,nparams,0:nproc_chunk1-1),intent(in) :: model1

  ! topology of the control points of the surface element
  integer,intent(in) :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL),dimension(nparams),intent(inout) :: model_maxdiff

  ! local parameters
  integer :: i,j,k,iglob,iker
  ! interpolated point location
  integer :: ispec_selected,rank_selected
  integer :: i_selected,j_selected,k_selected
  double precision :: xi,eta,gamma
  ! point location
  real(kind=CUSTOM_REAL) :: x_target,y_target,z_target
  real(kind=CUSTOM_REAL) :: val,val_initial

  ! checks given ispec
  if (ispec < 1 .or. ispec > nspec) then
    print*,'Error: rank ',myrank,' has invalid ispec'
    stop 'Error invalid ispec in get_model_values_bruteforce() routine'
  endif

  ! brute-force search over all points
  ! loops over all gll points
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
                    xi,eta,gamma,&
                    ispec_selected,rank_selected, &
                    nspec_max_old,nglob_max_old,nproc_chunk1, &
                    ibool1,x1,y1,z1, &
                    iaddx,iaddy,iaddr,xigll,yigll,zigll,typical_size, &
                    i_selected,j_selected,k_selected)

        ! interpolate model values
        do iker = 1,nparams
          call interpolate(xi,eta,gamma,ispec_selected, &
                           nspec_max_old,model1(:,:,:,:,iker,rank_selected), &
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

      enddo
    enddo
  enddo

  end subroutine get_model_values_bruteforce


!
!------------------------------------------------------------------------------
!

  subroutine get_model_values_kdtree(ispec,nspec,nglob,ibool2,x2,y2,z2,nparams,model2, &
                                     nspec_max_old,nglob_max_old,nproc_chunk1,ibool1,x1,y1,z1,model1, &
                                     iaddx,iaddy,iaddr,xigll,yigll,zigll,typical_size,myrank,model_maxdiff, &
                                     USE_MIDPOINT_SEARCH,DO_SEPARATION_410_650,DO_SEPARATION_TOPO,USE_FALLBACK)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,MIDX,MIDY,MIDZ,R_EARTH_KM,R_EARTH

  use kdtree_search, only: kdtree_find_nearest_neighbor,kdtree_nodes_location

  implicit none

  integer,intent(in) :: ispec
  integer,intent(in) :: nparams

  ! new, target mesh
  integer,intent(in) :: nspec,nglob
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: ibool2
  real(kind=CUSTOM_REAL),dimension(nglob),intent(in) :: x2,y2,z2
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec,nparams),intent(inout) :: model2

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec_max_old,nglob_max_old,nproc_chunk1
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec_max_old,0:nproc_chunk1-1),intent(in) :: ibool1
  real(kind=CUSTOM_REAL),dimension(nglob_max_old,0:nproc_chunk1-1),intent(in) :: x1,y1,z1
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_max_old,nparams,0:nproc_chunk1-1),intent(in) :: model1

  ! topology of the control points of the surface element
  integer,intent(in) :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size
  integer, intent(in) :: myrank

  real(kind=CUSTOM_REAL),dimension(nparams),intent(inout) :: model_maxdiff
  logical,intent(in) :: USE_MIDPOINT_SEARCH
  logical,intent(in) :: DO_SEPARATION_410_650,DO_SEPARATION_TOPO
  logical,intent(in) :: USE_FALLBACK

  ! local parameters
  integer :: i,j,k,iglob,iker
  ! interpolated point location
  integer :: ispec_selected,rank_selected
  double precision :: xi,eta,gamma
  ! point location
  real(kind=CUSTOM_REAL) :: x_target,y_target,z_target
  ! nodes search
  double precision,dimension(3) :: xyz_target
  double precision :: dist_min
  double precision :: mid_radius,elem_height,r
  integer :: iglob_min
  integer :: i_selected,j_selected,k_selected
  ! locations
  real(kind=CUSTOM_REAL) :: x_found,y_found,z_found
  real(kind=CUSTOM_REAL) :: val,val_initial
  logical :: is_critical,search_internal
  integer :: ii,jj,kk
  ! Earth radius
  double precision,parameter :: RTOP = R_EARTH      ! 6371000.d0
  double precision,parameter :: R220 = 6151000.d0
  double precision,parameter :: R410 = 5961000.d0
  double precision,parameter :: R600 = 5771000.d0
  double precision,parameter :: R650 = 5721000.d0
  ! margins
  ! surface topography and moho stretching: uses a 120-km margin
  ! (moho max ~ 80 km below Himalayan)
  double precision,parameter :: RTOP_MARGIN = 120000.d0
  ! 410-km: uses a 50-km margin
  ! (s362ani: 410 topography perturbations have min/max ~ -13/+13 km)
  ! (element heights are around 60km below and 64km above the 410-km discontinuity)
  double precision,parameter :: R410_MARGIN = 50000.d0
  ! 650-km: uses a 50-km margin
  ! (s362ani: 650 topography perturbations have min/max ~ -14/+19 km)
  double precision,parameter :: R650_MARGIN = 50000.d0

  ! debug warning about large model value differences
  logical,parameter :: DO_WARNING = .false.

  ! checks given ispec
  if (ispec < 1 .or. ispec > nspec) then
    print*,'Error: rank ',myrank,' has invalid ispec'
    stop 'Error invalid ispec in get_model_values_kdtree() routine'
  endif

  ! initializes for 410 special case
  mid_radius = 0.d0
  elem_height = 0.d0
  is_critical = .false.

  ! searches for element using mid-point
  if (USE_MIDPOINT_SEARCH .or. DO_SEPARATION_410_650 .or. DO_SEPARATION_TOPO) then
    iglob = ibool2(MIDX,MIDY,MIDZ,ispec)

    xyz_target(1) = x2(iglob)
    xyz_target(2) = y2(iglob)
    xyz_target(3) = z2(iglob)

    ! finds closest point in source chunk
    call kdtree_find_nearest_neighbor(xyz_target,iglob_min,dist_min)

    ! selected source rank
    rank_selected = int((iglob_min-1)/nspec_max_old)
    ! selected closest element
    ispec_selected = iglob_min - rank_selected * nspec_max_old

    ! checks if mid-point was found properly
    call check_point_result()

    ! debug
    !if (myrank == 0 .and. iglob < 100) &
    !  print*,'dist_min kdtree midpoint: ',dist_min * R_EARTH_KM,'(km)',typical_size * R_EARTH_KM

    ! special case for 410-km/650-km discontinuity
    if (DO_SEPARATION_410_650) then
      ! point radius
      r = sqrt(xyz_target(1)*xyz_target(1) + xyz_target(2)*xyz_target(2) + xyz_target(3)*xyz_target(3))

      ! surface
      if (r >= R220/R_EARTH) then
        ! elements close to surface
        is_critical = .true.
      endif

      ! 410-km discontinuity
      if (r >= R600/R_EARTH .and. r <= R220/R_EARTH) then
        ! elements within 220 - 600 km depth
        is_critical = .true.
      endif

      ! 650-km discontinuity
      if (r >= (R650 - R650_MARGIN)/R_EARTH .and. r <= (R650 + R650_MARGIN)/R_EARTH) then
        ! elements within around 650 km depth
        is_critical = .true.
      endif
    endif

    ! special case for surface (moho) discontinuity
    if (DO_SEPARATION_TOPO) then
      ! point radius
      r = sqrt(xyz_target(1)*xyz_target(1) + xyz_target(2)*xyz_target(2) + xyz_target(3)*xyz_target(3))

      ! surface
      if (r >= R220/R_EARTH) then
        ! elements close to surface
        is_critical = .true.
      endif
    endif

    if (is_critical) then
      ! stores mid-point radius
      mid_radius = r

      ! element height: size along a vertical edge
      ! top point
      iglob = ibool2(1,1,NGLLZ,ispec)
      r = sqrt(x2(iglob)*x2(iglob) + y2(iglob)*y2(iglob) + z2(iglob)*z2(iglob))
      ! bottom point
      iglob = ibool2(1,1,1,ispec)
      elem_height = r - sqrt(x2(iglob)*x2(iglob) + y2(iglob)*y2(iglob) + z2(iglob)*z2(iglob))
      ! debug
      !if (myrank == 0) print*,'element height: ',elem_height * R_EARTH_KM,'(km)','radius: ',mid_radius*R_EARTH_KM
    endif

  endif

  ! loops over all element gll points
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
          search_internal = .false.

          ! avoids getting values from "wrong" side on 410-km discontinuity,etc.
          if (DO_SEPARATION_410_650) then
            if (is_critical) then
              ! gll point radius
              r = sqrt(x_target*x_target + y_target*y_target + z_target*z_target)

              ! takes corresponding internal gll point for element search
              ! 410-km discontinuity
              if (r >= (R410 - R410_MARGIN)/R_EARTH .and. r <= (R410 + R410_MARGIN)/R_EARTH) search_internal = .true.
              ! 650-km discontinuity
              if (r >= (R650 - R650_MARGIN)/R_EARTH .and. r <= (R650 + R650_MARGIN)/R_EARTH) search_internal = .true.
            endif
          endif

          if (DO_SEPARATION_TOPO) then
            if (is_critical) then
              ! gll point radius
              r = sqrt(x_target*x_target + y_target*y_target + z_target*z_target)

              ! takes corresponding internal gll point for element search
              ! surface elements
              if (r >= (RTOP - RTOP_MARGIN)/R_EARTH) search_internal = .true.
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

            ! target point location
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
          rank_selected = int((iglob_min-1)/nspec_max_old)
          ! selected closest element
          ispec_selected = iglob_min - rank_selected * nspec_max_old

          ! checks if point was found properly
          call check_point_result()

          ! debug
          !if (myrank == 0 .and. iglob < 100) &
          !  print*,'dist_min kdtree: ',dist_min * R_EARTH_KM,'(km)',typical_size * R_EARTH_KM

          ! restores original target point location for locating/interpolating
          iglob = ibool2(i,j,k,ispec)
          x_target = x2(iglob)
          y_target = y2(iglob)
          z_target = z2(iglob)
        endif

        ! gets interpolated position within selected element
        call locate_single(x_target,y_target,z_target, &
                           xi,eta,gamma,&
                           ispec_selected,rank_selected, &
                           nspec_max_old,nglob_max_old,nproc_chunk1, &
                           ibool1,x1,y1,z1, &
                           iaddx,iaddy,iaddr,xigll,yigll,zigll,typical_size, &
                           i_selected,j_selected,k_selected)

        ! checks closest gll point
        iglob = ibool1(i_selected,j_selected,k_selected,ispec_selected,rank_selected)
        x_found = x1(iglob,rank_selected)
        y_found = y1(iglob,rank_selected)
        z_found = z1(iglob,rank_selected)

        ! checks distance
        dist_min = sqrt((x_found-x_target)**2 + (y_found-y_target)**2 + (z_found-z_target)**2)
        if (dist_min > 2 * typical_size) then
          print*,'Warning: rank ',myrank,' - large dist_min: ',dist_min * R_EARTH_KM,'(km)', &
                 'element size:',typical_size * R_EARTH_KM
          print*,'target location:',xyz_target(:)
          print*,'target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_EARTH_KM,'(km)'
          print*,'gll location   :',x_found,y_found,z_found
          print*,'gll radius     :',sqrt(x_found**2 + y_found**2 + z_found**2) * R_EARTH_KM,'(km)'
          print*,'distance gll:',dist_min * R_EARTH_KM,'(km)'
          ! debug
          !stop 'Error gll model value invalid'
        endif
        ! debug
        !if (myrank == 0 .and. iglob < 100) &
        !  print*,'dist_min gll point: ',dist_min * R_EARTH_KM,'(km)',typical_size * R_EARTH_KM

        ! interpolate model values
        do iker = 1,nparams
          if (USE_FALLBACK) then
            call interpolate_limited(xi,eta,gamma,ispec_selected, &
                                     nspec_max_old,model1(:,:,:,:,iker,rank_selected), &
                                     val,xigll,yigll,zigll, &
                                     i_selected,j_selected,k_selected)
          else
            call interpolate(xi,eta,gamma,ispec_selected, &
                             nspec_max_old,model1(:,:,:,:,iker,rank_selected), &
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
              print*,'Warning: model ',iker,' value:',val,'is very different from initial value ',val_initial
              print*,'  rank ',myrank,' - dist_min: ',dist_min * R_EARTH_KM,'(km)'
              print*,'  element',ispec,'selected ispec:',ispec_selected,'in rank:',rank_selected,'iglob_min:',iglob_min
              print*,'  typical element size:',typical_size * 0.5 * R_EARTH_KM
              print*,'  interpolation i,j,k :',i_selected,j_selected,k_selected
              print*,'  interpolation       :',xi,eta,gamma
              print*,'  target location:',xyz_target(:)
              print*,'  target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_EARTH_KM,'(km)'
              print*,'  gll location   :',x_found,y_found,z_found
              print*,'  gll radius     :',sqrt(x_found**2 + y_found**2 + z_found**2) * R_EARTH_KM,'(km)'
              print*,'  distance gll:',dist_min * R_EARTH_KM,'(km)'
              !stop 'Error model value invalid'
            endif
          endif

          ! debug
          !if (myrank == 0 .and. iglob < 100) &
          !  print*,'new model ',iker,': value ',val,'initial ',val_initial,'diff ',(val - val_initial)/val_initial*100.0,'(%)'

        enddo

      enddo
    enddo
  enddo

  contains

    subroutine check_point_result()

    implicit none

    ! checks valid iglob
    if (iglob_min < 1 .or. iglob_min > nspec_max_old * nproc_chunk1) then
      print*,'Error iglob_min :',iglob_min
      print*,'nspec / nproc :',nspec_max_old,nproc_chunk1
      stop 'Error invalid iglob_min index'
    endif

    ! checks valid rank
    if (rank_selected < 0 .or. rank_selected >= nproc_chunk1) then
      print*,'Error rank:',myrank,'invalid selected rank ',rank_selected,'for element',ispec
      print*,'target location:',xyz_target(:)
      stop 'Error specifying closest rank for element'
    endif

    ! checks valid ispec
    if (ispec_selected < 1 .or. ispec_selected > nspec_max_old) then
      print*,'Error rank:',myrank,'invalid selected ispec ',ispec_selected,'for element',ispec
      print*,'rank_selected:',rank_selected,'iglob_min:',iglob_min,'nspec_max_old:',nspec_max_old
      print*,'target location:',xyz_target(:)
      print*,'dist_min: ',dist_min * R_EARTH_KM,'(km)'
      stop 'Error specifying closest ispec element'
    endif

    ! checks minimum distance within a typical element size
    if (dist_min > 2 * typical_size) then
      print*,'Warning: rank ',myrank,' - large dist_min: ',dist_min * R_EARTH_KM,'(km)', &
             'element size:',typical_size * R_EARTH_KM
      print*,'element',ispec,'selected ispec:',ispec_selected,'in rank:',rank_selected,'iglob_min:',iglob_min
      print*,'typical element size:',typical_size * 0.5 * R_EARTH_KM
      print*,'target location:',xyz_target(:)
      print*,'target radius  :',sqrt(xyz_target(1)**2 + xyz_target(2)**2 + xyz_target(3)**2) * R_EARTH_KM,'(km)'
      print*,'found location :',kdtree_nodes_location(:,iglob_min)
      print*,'found radius   :',sqrt(kdtree_nodes_location(1,iglob_min)**2 &
                                   + kdtree_nodes_location(2,iglob_min)**2 &
                                   + kdtree_nodes_location(3,iglob_min)**2 ) * R_EARTH_KM,'(km)'
      !debug
      !stop 'Error dist_min too large in check_point_result() routine'
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
                    iaddx,iaddy,iaddr, &
                    xigll,yigll,zigll,typical_size, &
                    i_selected,j_selected,k_selected)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,MIDX,MIDY,MIDZ,HUGEVAL,NUM_ITER,R_EARTH

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

  ! topology of the control points of the surface element
  integer,intent(in) :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size

  integer, intent(out) :: i_selected,j_selected,k_selected

  ! local parameters
  integer :: i,j,k,ispec,iglob,rank
  double precision :: dist,dist_typical_sq
  double precision :: distmin

  ! set distance to huge initial value
  distmin = HUGEVAL
  dist_typical_sq = typical_size * typical_size

  ispec_selected = 0
  rank_selected = -1

  ! finds closest point
  do rank=0,nproc-1
    do ispec=1,nspec
      ! distance to cell center
      ! midpoint
      iglob = ibool(MIDX,MIDY,MIDZ,ispec,rank)
      ! distance to midpoint
      dist =   (x_target - dble(xstore(iglob,rank)))**2 &
              +(y_target - dble(ystore(iglob,rank)))**2 &
              +(z_target - dble(zstore(iglob,rank)))**2

      if (dist > dist_typical_sq) cycle

      ! loop only on points inside the element
      ! exclude edges to ensure this point is not shared with other elements
      do k=2,NGLLZ-1
        do j=2,NGLLY-1
          do i=2,NGLLX-1
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
                     ispec_selected,rank_selected, &
                     nspec,nglob,nproc, &
                     ibool,xstore,ystore,zstore, &
                     iaddx,iaddy,iaddr, &
                     xigll,yigll,zigll,typical_size, &
                     i_selected,j_selected,k_selected)

  end subroutine locate


!
!------------------------------------------------------------------------------
!

  subroutine locate_single(x_target,y_target,z_target, &
                           xi_target,eta_target,gamma_target, &
                           ispec_selected,rank_selected, &
                           nspec,nglob,nproc, &
                           ibool,xstore,ystore,zstore, &
                           iaddx,iaddy,iaddr, &
                           xigll,yigll,zigll,typical_size, &
                           i_selected,j_selected,k_selected)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NGNOD,HUGEVAL,NUM_ITER,R_EARTH_KM

  implicit none

  ! point location
  real(kind=CUSTOM_REAL),intent(in) :: x_target,y_target,z_target

  integer,intent(in) :: ispec_selected,rank_selected

  ! interpolated point location
  double precision,intent(out) :: xi_target,eta_target,gamma_target

  ! for old, first mesh we interpolate on
  integer,intent(in) :: nspec,nglob,nproc

  ! arrays containing coordinates of the points
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec,0:nproc-1),intent(in) :: ibool
  real(kind=CUSTOM_REAL),dimension(nglob,0:nproc-1),intent(in) :: xstore,ystore,zstore

  ! topology of the control points of the surface element
  integer,intent(in) :: iaddx(NGNOD),iaddy(NGNOD),iaddr(NGNOD)

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! typical element size in old mesh used as search radius
  double precision,intent(in) :: typical_size

  ! best guess for i,j,k
  integer,intent(out) :: i_selected,j_selected,k_selected

  ! local parameters
  ! use integer array to store values
  integer :: ix_initial_guess,iy_initial_guess,iz_initial_guess
  integer :: i,j,k,iglob !,ispec,rank
  double precision :: dist
  double precision :: xi,eta,gamma,dx,dy,dz,dxi,deta,dgamma
  integer :: iax,iay,iaz
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
  double precision :: distmin

  !------------------------------------------------------

  ! exact position search
  logical,parameter :: DO_REFINE_LOCATION = .true.
  logical,parameter :: DO_WARNING = .true.

  !------------------------------------------------------

  ! checks
  if (ispec_selected < 1 ) stop 'Error locate_single: specifying closest ispec element'
  if (rank_selected < 0 ) stop 'Error locate_single: specifying closest rank for element'

  ! set distance to huge initial value
  distmin = HUGEVAL
  ix_initial_guess = 0
  iy_initial_guess = 0
  iz_initial_guess = 0

  ! finds closest interior gll point
  do k=1,NGLLZ
    do j=1,NGLLY
      do i=1,NGLLX
        iglob = ibool(i,j,k,ispec_selected,rank_selected)

        dist = (x_target - xstore(iglob,rank_selected))**2 &
              +(y_target - ystore(iglob,rank_selected))**2 &
              +(z_target - zstore(iglob,rank_selected))**2

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

  ! initial minimum distance in km
  distmin = sqrt(distmin) * R_EARTH_KM

  ! debug
  !print*,'distmin = ',sngl(distmin),'(km)'

  ! find the best (xi,eta)
  ! use initial guess in xi, eta and gamma from closest point found
  xi = xigll(ix_initial_guess)
  eta = yigll(iy_initial_guess)
  gamma = zigll(iz_initial_guess)

  ! iterate to solve the non linear system to find the exact position within the element
  if (DO_REFINE_LOCATION) then

    ! define coordinates of the control points of the element
    do ia=1,NGNOD
      if (iaddx(ia) == 0) then
        iax = 1
      else if (iaddx(ia) == 1) then
        iax = (NGLLX+1)/2
      else if (iaddx(ia) == 2) then
        iax = NGLLX
      else
        stop 'incorrect value of iaddx'
      endif
      if (iaddy(ia) == 0) then
        iay = 1
      else if (iaddy(ia) == 1) then
        iay = (NGLLY+1)/2
      else if (iaddy(ia) == 2) then
        iay = NGLLY
      else
        stop 'incorrect value of iaddy'
      endif
      if (iaddr(ia) == 0) then
        iaz = 1
      else if (iaddr(ia) == 1) then
        iaz = (NGLLZ+1)/2
      else if (iaddr(ia) == 2) then
        iaz = NGLLZ
      else
        stop 'incorrect value of iaddr'
      endif

      iglob = ibool(iax,iay,iaz,ispec_selected,rank_selected)

      xelm(ia) = dble(xstore(iglob,rank_selected))
      yelm(ia) = dble(ystore(iglob,rank_selected))
      zelm(ia) = dble(zstore(iglob,rank_selected))
    enddo

    dxi = 0.d0
    deta = 0.d0
    dgamma = 0.d0

    ! loop leads to invalid jacobian... probably some gll points are too far outside of the selected element
    do iter_loop = 1,2*NUM_ITER

      ! recompute jacobian for the new point
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)

! debug
!      if (ier /= 0) then
!        ! debug
!        if (.true.) then
!          print*,'jacobian error in locate_single(): '
!          print*,'jacobian error i,j,k,ispec :',ix_initial_guess,iy_initial_guess,iz_initial_guess,ispec_selected
!          print*,'jacobian error iter_loop   :',iter_loop
!          dist = sqrt((x_target-x)**2+(y_target-y)**2+(z_target-z)**2) * R_EARTH_KM
!          print*,'jacobian error dist        :',dist,'(km)',distmin
!        endif
!
!        ! uses initial guess again
!        xi = xigll(ix_initial_guess)
!        eta = yigll(iy_initial_guess)
!        gamma = zigll(iz_initial_guess)
!
!        ! uses previous guess
!        !xi = xi - dxi
!        !eta = eta - deta
!        !gamma = gamma - dgamma
!
!        ! exits loop
!        exit
!
!        !stop 'Error recomputing jacobian'
!      endif
!

      ! compute distance to target location
      dx = - (x - x_target)
      dy = - (y - y_target)
      dz = - (z - z_target)

      ! compute increments
      dxi  = xix*dx + xiy*dy + xiz*dz
      deta = etax*dx + etay*dy + etaz*dz
      dgamma = gammax*dx + gammay*dy + gammaz*dz

      ! impose limit on increments
      if (abs(dxi) > 0.1d0 ) dxi = sign(1.0d0,dxi)*0.1d0
      if (abs(deta) > 0.1d0 ) deta = sign(1.0d0,deta)*0.1d0
      if (abs(dgamma) > 0.1d0 ) dgamma = sign(1.0d0,dgamma)*0.1d0

      ! update values
      xi = xi + dxi
      eta = eta + deta
      gamma = gamma + dgamma

      ! impose that we stay in that element
      ! (useful if user gives a receiver outside the mesh for instance)
      ! we can go slightly outside the [1,1] segment since with finite elements
      ! the polynomial solution is defined everywhere
      ! can be useful for convergence of iterative scheme with distorted elements
      !if (xi > 1.10d0) xi = 1.10d0
      !if (xi < -1.10d0) xi = -1.10d0
      !if (eta > 1.10d0) eta = 1.10d0
      !if (eta < -1.10d0) eta = -1.10d0
      !if (gamma > 1.10d0) gamma = 1.10d0
      !if (gamma < -1.10d0) gamma = -1.10d0

      ! point leaves element, stay to closest guess
      if (xi > 1.10d0 .or. xi < -1.10d0 .or. eta > 1.10d0 .or. eta < -1.10d0 .or. gamma > 1.10d0 .or. gamma < -1.10d0) then
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
    ! compute final distance between asked and found (converted to km)
    final_distance = sqrt((x_target-x)**2+(y_target-y)**2+(z_target-z)**2) * R_EARTH_KM

    ! debug
    !if (final_distance > 5.0 ) &
    !  print*,'final distance = ',sngl(final_distance),'(km)',distmin,xi,eta,gamma

    ! checks if location improved
    if (distmin <= final_distance) then
      ! uses initial guess
      xi = xigll(ix_initial_guess)
      eta = yigll(iy_initial_guess)
      gamma = zigll(iz_initial_guess)
      final_distance = distmin
      call recompute_jacobian(xelm,yelm,zelm,xi,eta,gamma,x,y,z, &
                              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz)
    endif

    ! add warning if estimate is poor
    ! (usually means receiver outside the mesh given by the user)
    if (DO_WARNING) then
      if (final_distance > typical_size * R_EARTH_KM) then
        print*, '*****************************************************************'
        print*, '***** WARNING: location estimate is poor                    *****'
        print*, '*****************************************************************'
        print*, 'closest estimate found: ',sngl(final_distance),'km away',' - not within:',typical_size * R_EARTH_KM
        print*, ' in rank ',rank_selected,' in element ',ispec_selected,ix_initial_guess,iy_initial_guess,iz_initial_guess
        print*, ' at xi,eta,gamma coordinates = ',xi,eta,gamma
        print*, ' at radius ',sqrt(x**2 + y**2 + z**2) * R_EARTH_KM,'(km)'
        print*, ' initial distance :',distmin,'(km)'
      endif
    endif

    ! checks valid distance
    if (final_distance == HUGEVAL) stop 'Error locating location'
  endif

  ! return xi,eta and gamma of point found
  xi_target = xi
  eta_target = eta
  gamma_target = gamma

  i_selected = ix_initial_guess
  j_selected = iy_initial_guess
  k_selected = iz_initial_guess

  end subroutine locate_single


!
!------------------------------------------------------------------------------
!

  subroutine interpolate(xi,eta,gamma,ielem, &
                         nspec,model, &
                         val,xigll,yigll,zigll)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  double precision,intent(in):: xi,eta,gamma
  integer,intent(in):: ielem

  integer,intent(in):: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: model

  real(kind=CUSTOM_REAL),intent(out) :: val

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! local parameters
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
                      hgammar(NGLLZ), hpgammar(NGLLZ)
  integer:: i,j,k

  ! interpolation weights
  call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates value
  val = 0.0
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
          val = val +  hxir(i) * hetar(j) * hgammar(k) * model(i,j,k,ielem)
      enddo
    enddo
  enddo

  end subroutine interpolate

!
!------------------------------------------------------------------------------
!

  subroutine interpolate_limited(xi,eta,gamma,ielem, &
                                 nspec,model, &
                                 val,xigll,yigll,zigll, &
                                 i_selected,j_selected,k_selected)


  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  double precision,intent(in):: xi,eta,gamma
  integer,intent(in):: ielem

  integer,intent(in):: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(in) :: model

  real(kind=CUSTOM_REAL),intent(out) :: val

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  integer,intent(in):: i_selected,j_selected,k_selected

  ! local parameters
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
                      hgammar(NGLLZ), hpgammar(NGLLZ)
  integer:: i,j,k
  real(kind=CUSTOM_REAL) :: val_initial,val_avg,pert,pert_limit

  ! percentage
  real(kind=CUSTOM_REAL), parameter :: PERCENTAGE_LIMIT = 0.01

  ! interpolation weights
  call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates value
  val = 0.0_CUSTOM_REAL
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        val = val + hxir(i) * hetar(j) * hgammar(k) * model(i,j,k,ielem)
      enddo
    enddo
  enddo

  ! note: interpolation of values close to the surface or 3D moho encounters problems;
  !       this is a fall-back to the closest point value
  !
  ! uses average/closest point value if too far off

  ! closest point value
  val_initial = model(i_selected,j_selected,k_selected,ielem)

  ! average value
  val_avg = sum(model(:,:,:,ielem)) / NGLLX / NGLLY / NGLLZ

  ! initial difference
  pert = abs(val_initial - val_avg)

  ! upper/lower perturbation bound
  pert_limit = PERCENTAGE_LIMIT * abs(val_avg)
  if (pert > pert_limit) pert_limit = pert

  ! within a certain percentage range
  if (abs(val - val_avg ) > pert_limit) then
    val = val_initial
  endif

  end subroutine interpolate_limited

