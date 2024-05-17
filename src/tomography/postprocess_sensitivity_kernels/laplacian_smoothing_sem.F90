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

#include "config.fh"

program smooth_laplacian_sem

  use constants, only: myrank

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM,IIN,IOUT, &
       GAUSSALPHA,GAUSSBETA,MAX_STRING_LEN,GRAV,PI,TINYVAL

  use shared_parameters, only: R_PLANET_KM,CRUSTAL

  use postprocess_par, only: &
       NCHUNKS_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,NPROCTOT_VAL, &
       NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,MAX_KERNEL_NAMES,LOCAL_PATH

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  use adios_helpers_mod
  use manager_adios
#endif

  implicit none

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer, parameter :: NARGS = 7
  character(len=*), parameter :: reg_name = 'reg1/'
#else
  integer, parameter :: NARGS = 6
  character(len=*), parameter :: reg_name = '_reg1_'
#endif

  integer :: nspec, nglob, nker, niter_cg_max
  integer :: iker, i, j, k, idof, iel, i1, i2, ier, sizeprocs

  double precision    :: Lh, Lv, conv_crit, Lh2, Lv2
  double precision    :: x, y, z, rel_to_prem
  double precision    :: r, theta, phi
  double precision    :: rho,drhodr,vp,vs,Qkappa,Qmu

  real(kind=CUSTOM_REAL), dimension(:),       allocatable :: m, s
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: mo, so

  ! Mesh parameters
  real(kind=CUSTOM_REAL), dimension(:),       allocatable :: xglob, yglob, zglob, valglob
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dxsi_dx, deta_dx, dgam_dx
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dxsi_dy, deta_dy, dgam_dy
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: dxsi_dz, deta_dz, dgam_dz
  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: jacobian

  integer, dimension(:,:,:,:), allocatable :: ibool
  integer, dimension(:),       allocatable :: idoubling, ispec_is_tiso, my_neighbors, &
                                              nibool_interfaces
  integer, dimension(:,:),     allocatable :: ibool_interfaces

  integer :: max_nibool_interfaces, num_interfaces

  double precision, dimension(NGLLX)             :: wx, xgll
  double precision, dimension(NGLLY)             :: wy, ygll
  double precision, dimension(NGLLZ)             :: wz, zgll

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX)       :: dlagx,   dlagy,   dlagz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX)       :: dlagxwx, dlagywy, dlagzwz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX)       :: wxy,     wxz,     wyz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX,NGLLX) :: wxyz


  character(len=MAX_STRING_LEN),dimension(NARGS) :: arg
  character(len=MAX_STRING_LEN),dimension(:),allocatable :: kernel_names
  character(len=MAX_STRING_LEN) :: kernel_names_comma_delimited
  character(len=MAX_STRING_LEN) :: kernel_name, topo_dir

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  character(len=MAX_STRING_LEN) :: input_file, solver_file, solver_file_mpi, tmp_dir
  character(len=MAX_STRING_LEN) :: varname, tmp_str
  integer, dimension(:),       allocatable :: ibool_interfaces_tmp
#else
  character(len=MAX_STRING_LEN) :: input_dir, output_dir
  character(len=MAX_STRING_LEN) :: local_data_file
  character(len=MAX_STRING_LEN) :: filename
#endif

  character(len=MAX_STRING_LEN) :: output_file

  real(kind=CUSTOM_REAL) :: dxsi_dxl, deta_dxl, dgam_dxl
  real(kind=CUSTOM_REAL) :: dxsi_dyl, deta_dyl, dgam_dyl
  real(kind=CUSTOM_REAL) :: dxsi_dzl, deta_dzl, dgam_dzl, jacobianl

  real(kind=CUSTOM_REAL) :: norm_kerl, norm_ker

  ! timing
  integer :: ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start_all
  double precision :: tCPU
  double precision, external :: wtime

  ! Hessian
  logical :: is_hess
  logical :: is_kernel

  ! ADIOS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  integer(kind=8) :: group_size_inc,local_dim
#endif

  ! number of steps to reach 100 percent, i.e. 10 outputs info for every 10 percent
  integer,parameter :: NSTEP_PERCENT_INFO = 10

  double precision, external :: lagrange_deriv_GLL

  ! local copies of mesh parameters
  integer :: NPROC_XI
  integer :: NPROC_ETA
  integer :: NCHUNKS
  integer :: NSPEC_AB
  integer :: NGLOB_AB

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  ! Parameters for convergence
  conv_crit    = 1e-5    ! convergence criterion
  niter_cg_max = 1000    ! max number of iteration

 ! check command line arguments
  if (command_argument_count() /= NARGS .and. command_argument_count() /= NARGS-1 .and. command_argument_count() /= NARGS-2) then
     if (myrank == 0) then
#ifdef USE_ADIOS_INSTEAD_OF_MESH
        print *,'Usage: mpirun -np NPROC bin/xsmooth_laplacian_sem_adios SIGMA_H SIGMA_V KERNEL_NAME', &
               ' INPUT_FILE SOLVER_DIR OUTPUT_FILE'
        print *,'   with'
        print *,'     SIGMA_H, SIGMA_V - Horizontal and vertical smoothing lenghts'
        print *,'     KERNEL_NAME       - comma-separated kernel names (e.g., alpha_kernel,beta_kernel)'
        print *,'     INPUT_FILE        - ADIOS file with kernel values (e.g., kernels.bp)'
        print *,'     SOLVER_DIR        - directory w/ ADIOS file with mesh arrays (e.g., DATABASES_MPI/) containing', &
               ' solver_data.bp solver_data_mpi.bp'
        print *,'     OUTPUT_FILE       - ADIOS file for smoothed output'
        print *,'     REL_TO_PREM       - (optional) increase smoothing radius based on PREM model P-wave velocity'
        print *
#else
        print *,'Usage: mpirun -np NPROC bin/xsmooth_laplacian_sem SIGMA_H SIGMA_V KERNEL_NAME INPUT_DIR OUPUT_DIR'
        print *,'   with'
        print *,'     SIGMA_H, SIGMA_V - Horizontal and vertical smoothing lenghts'
        print *,'     KERNEL_NAME       - comma-separated kernel names (e.g., alpha_kernel,beta_kernel)'
        print *,'     INPUT_DIR         - directory with kernel files (e.g., proc***_alpha_kernel.bin)'
        print *,'     OUTPUT_DIR        - directory for smoothed output files'
        print *,'     REL_TO_PREM       - (optional) increase smoothing radius based on PREM model P-wave velocity'
        print *
#endif
        stop ' Please check command line arguments'
     endif
  endif
  call synchronize_all()

  if (myrank == 0) then
     print *, 'Running smooth_laplacian_SEM'
     print *
  endif
  call synchronize_all()

  ! timing
  time_start_all = wtime()

  ! allocates arrays
  allocate(kernel_names(MAX_KERNEL_NAMES),stat=ier)
  if (ier /= 0) stop 'Error allocating kernel_names array'
  kernel_names(:) = ''

  ! parse command line arguments
  do i = 1, NARGS
     call get_command_argument(i,arg(i))
     if (i <= 5 .and. trim(arg(i)) == '') then
        stop ' Please check command line arguments'
     endif
  enddo
  read(arg(1),*) Lh
  read(arg(2),*) Lv
  kernel_names_comma_delimited = arg(3)

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS arguments
  input_file = arg(4)
  tmp_dir     = arg(5)
  output_file = arg(6)
  solver_file = get_adios_filename(trim(tmp_dir)//'solver_data')
  solver_file_mpi = get_adios_filename(trim(tmp_dir)//'solver_data_mpi')
#else
  input_dir = arg(4)
  output_dir = arg(5)
#endif

  if (trim(arg(NARGS)) == '') then
      rel_to_prem = 0.0
  else
      read(arg(NARGS),*) rel_to_prem
      if (myrank == 0) print *, 'Increase smoothing length based on PREM model     :', rel_to_prem
  endif

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
    write(*,*) 'mesh parameters (from local_path parameter):'
    write(*,*) '  LOCAL_PATH         = ',trim(LOCAL_PATH)
    write(*,*) '  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
    write(*,*) '  NPROCTOT           = ',NPROCTOT_VAL
    write(*,*)
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

  ! user output
  if (myrank == 0) then
     print *,"defaults:"
     print *,"  NPROC_XI , NPROC_ETA        : ",NPROC_XI,NPROC_ETA
     print *,"  NCHUNKS                     : ",NCHUNKS
     print *
     print *,"  smoothing sigma_xy , sigma_z (km)                : ",Lh,Lv
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

  Lh = Lh / real(R_PLANET_KM,kind=CUSTOM_REAL) ! scale
  Lv = Lv / real(R_PLANET_KM,kind=CUSTOM_REAL) ! scale

  if (myrank == 0 .and. abs(Lh - Lv) > TINYVAL) then
    print *, 'WARNING: different smoothing length in XY and Z direction'
  endif
!   taper_vertical = taper_vertical / real(R_PLANET_KM,kind=CUSTOM_REAL) ! scale

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

  ! 1. Read inputs and prepare MPI / gpu
  if (.not. allocated(ibool))         allocate(ibool(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(valglob))       allocate(valglob(nglob_ab))
  if (.not. allocated(dxsi_dx))       allocate(dxsi_dx(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(dxsi_dy))       allocate(dxsi_dy(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(dxsi_dz))       allocate(dxsi_dz(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(deta_dx))       allocate(deta_dx(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(deta_dy))       allocate(deta_dy(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(deta_dz))       allocate(deta_dz(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(dgam_dx))       allocate(dgam_dx(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(dgam_dy))       allocate(dgam_dy(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(dgam_dz))       allocate(dgam_dz(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(jacobian))      allocate(jacobian(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(xglob))         allocate(xglob(nglob_ab))
  if (.not. allocated(yglob))         allocate(yglob(nglob_ab))
  if (.not. allocated(zglob))         allocate(zglob(nglob_ab))
  if (.not. allocated(idoubling))     allocate(idoubling(nspec_ab))
  if (.not. allocated(ispec_is_tiso)) allocate(ispec_is_tiso(nspec_ab))

  if (.not. allocated(m))  allocate(m(nglob_ab))
  if (.not. allocated(s))  allocate(s(nglob_ab))
  if (.not. allocated(mo)) allocate(mo(ngllx, nglly, ngllz, nspec_ab))
  if (.not. allocated(so)) allocate(so(ngllx, nglly, ngllz, nspec_ab))

  m(:) = 0
  s(:) = 0

  ! GLL points weights
  call zwgljd(xgll,wx,NGLLX,GAUSSALPHA,GAUSSBETA)
  call zwgljd(ygll,wy,NGLLY,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zgll,wz,NGLLZ,GAUSSALPHA,GAUSSBETA)
  do k = 1,NGLLZ
     do j = 1,NGLLY
        do i = 1,NGLLX
           wxyz(i,j,k) = real(wx(i)*wy(j)*wz(k),kind=CUSTOM_REAL)
        enddo
     enddo
  enddo
  ! Get derivative matrices
  do i1 = 1,NGLLX
     do i2 = 1,NGLLX
        dlagx(i2,i1)   = real(lagrange_deriv_GLL(i1-1,i2-1,xgll,NGLLX), kind=CUSTOM_REAL)
        dlagy(i2,i1)   = real(lagrange_deriv_GLL(i1-1,i2-1,xgll,NGLLX), kind=CUSTOM_REAL)
        dlagz(i2,i1)   = real(lagrange_deriv_GLL(i1-1,i2-1,xgll,NGLLX), kind=CUSTOM_REAL)
        dlagxwx(i2,i1) = real(lagrange_deriv_GLL(i1-1,i2-1,xgll,NGLLX)*wx(i2), kind=CUSTOM_REAL)
        dlagywy(i2,i1) = real(lagrange_deriv_GLL(i1-1,i2-1,xgll,NGLLX)*wy(i2), kind=CUSTOM_REAL)
        dlagzwz(i2,i1) = real(lagrange_deriv_GLL(i1-1,i2-1,xgll,NGLLX)*wz(i2), kind=CUSTOM_REAL)
        wxy(i1, i2)    = real(wx(i1)*wy(i2), kind=CUSTOM_REAL)
        wyz(i1, i2)    = real(wy(i1)*wz(i2), kind=CUSTOM_REAL)
        wxz(i1, i2)    = real(wx(i1)*wz(i2), kind=CUSTOM_REAL)
     enddo
  enddo


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
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "ibool", ibool(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xixstore", dxsi_dx(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xiystore", dxsi_dy(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "xizstore", dxsi_dz(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammaxstore", dgam_dx(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammaystore", dgam_dy(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "gammazstore", dgam_dz(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etaxstore", deta_dx(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etaystore", deta_dy(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nspec, trim(reg_name) // "etazstore", deta_dz(:, :, :, :))
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "x_global", xglob)
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "y_global", yglob)
  call read_adios_array(myadios_file, myadios_group, myrank, nglob, trim(reg_name) // "z_global", zglob)

  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"MeshReader")

  ! opens file for mesh
  if (myrank == 0) print *, 'reading in ADIOS solver file: ',trim(solver_file_mpi)

  call init_adios_group(myadios_group, "MeshMPIReader")
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,solver_file_mpi)

  ! MPI interfaces
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(reg_name) // "num_interfaces",num_interfaces)
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(reg_name) // "max_nibool_interfaces",max_nibool_interfaces)

  allocate(my_neighbors(num_interfaces), &
           nibool_interfaces(num_interfaces), stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_crust_mantle etc.')
  if (num_interfaces > 0) then
     allocate(ibool_interfaces(max_nibool_interfaces,num_interfaces), stat=ier)
     if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
     allocate(ibool_interfaces_tmp(max_nibool_interfaces * num_interfaces), stat=ier)
     if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle_tmp')
     call read_adios_array(myadios_file, myadios_group, myrank, num_interfaces, trim(reg_name) // "my_neighbors", &
             my_neighbors(:))
     call read_adios_array(myadios_file, myadios_group, myrank, num_interfaces, trim(reg_name) // "nibool_interfaces", &
             nibool_interfaces(:))
     call read_adios_array(myadios_file, myadios_group, myrank, num_interfaces*max_nibool_interfaces, &
        trim(reg_name) // "ibool_interfaces", ibool_interfaces_tmp(:))
     k = 0
     do j = 1, num_interfaces
        do i = 1, max_nibool_interfaces
            k = k + 1
            ibool_interfaces(i, j) = ibool_interfaces_tmp(k)
        enddo
     enddo
     deallocate(ibool_interfaces_tmp)
  else
     max_nibool_interfaces = 0
     allocate(ibool_interfaces(0,0),stat=ier)
     if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_crust_mantle')
  endif

  call close_file_adios_read_and_finalize_method(myadios_file)
  call delete_adios_group(myadios_group,"MeshMPIReader")

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
  read(IIN) xglob
  read(IIN) yglob
  read(IIN) zglob

  read(IIN) ibool
  read(IIN) idoubling
  read(IIN) ispec_is_tiso

  read(IIN) dxsi_dx
  read(IIN) dxsi_dy
  read(IIN) dxsi_dz
  read(IIN) deta_dx
  read(IIN) deta_dy
  read(IIN) deta_dz
  read(IIN) dgam_dx
  read(IIN) dgam_dy
  read(IIN) dgam_dz

  close(IIN)

  !! Read MPI
  write(filename,'(a,i6.6,a)') trim(topo_dir)//'/proc',myrank,trim(reg_name)//'solver_data_mpi.bin'
  open(unit=IIN,file=filename, &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')
  ! MPI interfaces
  read(IIN) num_interfaces
  allocate(my_neighbors(num_interfaces), &
           nibool_interfaces(num_interfaces), stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array my_neighbors_crust_mantle etc.')
  if (num_interfaces > 0) then
     read(IIN) max_nibool_interfaces
     allocate(ibool_interfaces(max_nibool_interfaces,num_interfaces), &
          stat=ier)
     if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
     read(IIN) my_neighbors
     read(IIN) nibool_interfaces
     read(IIN) ibool_interfaces
  else
     ! dummy array
     max_nibool_interfaces = 0
     allocate(ibool_interfaces(0,0),stat=ier)
     if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_crust_mantle')
  endif
  close(IIN)

  deallocate(idoubling, ispec_is_tiso)
#endif
  call synchronize_all()

  ! 2. Modify jacobian and inverse jacobian for non-dimensionalize smoothing
  !    note: correlation mengths Lx, Ly and Lz could be made depth or space dependent
  do iel = 1, nspec_ab
    do k = 1, ngllz
      do j = 1, nglly
        do i = 1, ngllx
          ! Get inverse jacobian
          dxsi_dxl = dxsi_dx(i,j,k,iel)
          deta_dxl = deta_dx(i,j,k,iel)
          dgam_dxl = dgam_dx(i,j,k,iel)
          dxsi_dyl = dxsi_dy(i,j,k,iel)
          deta_dyl = deta_dy(i,j,k,iel)
          dgam_dyl = dgam_dy(i,j,k,iel)
          dxsi_dzl = dxsi_dz(i,j,k,iel)
          deta_dzl = deta_dz(i,j,k,iel)
          dgam_dzl = dgam_dz(i,j,k,iel)
          ! Compute jacobian
          jacobianl = dxsi_dxl * (deta_dyl * dgam_dzl - deta_dzl * dgam_dyl) &
                    - dxsi_dyl * (deta_dxl * dgam_dzl - deta_dzl * dgam_dxl) &
                    + dxsi_dzl * (deta_dxl * dgam_dyl - deta_dyl * dgam_dxl)
          jacobianl = 1. / jacobianl
          ! convert Lv, Lh to Lx, Ly, Lz
          idof = ibool(i, j, k, iel)
          x = xglob(idof)
          y = yglob(idof)
          z = zglob(idof)

          ! converts x/y/z to geocentric coordinates r/theta/phi
          call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)

          ! theta/phi not used any further...
          ! note: converts the geocentric colatitude to a geographic colatitude by ellipticity correction
          !if (ELLIPTICITY) call geocentric_2_geographic_dble(theta,theta)

          ! determine smoothing radius
          Lh2 = Lh
          Lv2 = Lv
          !   if (taper_vertical > 0) then
          !       if ((1 - r) < taper_vertical) then
          !          Lv2 = Lv + (Lh - Lv) * (1 - cos(PI * (1 - r) / taper_vertical)) / 2
          !       else
          !          Lv2 = Lh
          !       endif
          !   else
          !       Lv2 = Lv
          !   endif
          ! increase radius based on PREM model velocity
          if (rel_to_prem > TINYVAL) then
            call model_prem_iso(r,rho,drhodr,vp,vs,Qkappa,Qmu,0,CRUSTAL,.false.)
            if (vp > 1.0) then
              Lv2 = Lv2 * vp ** rel_to_prem
              Lh2 = Lh2 * vp ** rel_to_prem
            endif
          endif
          !   ! convert Lv, Lh to Lx, Ly, Lz
          !   e2 = 1 - (Lv2 / Lh2) ** 2
          !   Lz = Lh2 * (1 - e2 * cos(theta) ** 2)
          !   Lx = Lh2 * (1 - e2 * (sin(theta) * cos(phi)) ** 2 )
          !   Ly = Lh2 * (1 - e2 * (sin(theta) * sin(phi)) ** 2 )

          ! Apply scaling
          dxsi_dx(i,j,k,iel)  = real(dxsi_dxl * Lh2,kind=CUSTOM_REAL)
          deta_dx(i,j,k,iel)  = real(deta_dxl * Lh2,kind=CUSTOM_REAL)
          dgam_dx(i,j,k,iel)  = real(dgam_dxl * Lv2,kind=CUSTOM_REAL)
          dxsi_dy(i,j,k,iel)  = real(dxsi_dyl * Lh2,kind=CUSTOM_REAL)
          deta_dy(i,j,k,iel)  = real(deta_dyl * Lh2,kind=CUSTOM_REAL)
          dgam_dy(i,j,k,iel)  = real(dgam_dyl * Lv2,kind=CUSTOM_REAL)
          dxsi_dz(i,j,k,iel)  = real(dxsi_dzl * Lh2,kind=CUSTOM_REAL)
          deta_dz(i,j,k,iel)  = real(deta_dzl * Lh2,kind=CUSTOM_REAL)
          dgam_dz(i,j,k,iel)  = real(dgam_dzl * Lv2,kind=CUSTOM_REAL)
          jacobian(i,j,k,iel) = real(jacobianl / (Lh*Lh*Lv),kind=CUSTOM_REAL)
        enddo
      enddo
    enddo
  enddo
  call synchronize_all()


  !! PREPARE ADIOS OUTPUTS
#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! user output
  if (myrank == 0) then
    print *
    print *, 'writing to ADIOS output file: ',trim(output_file)
    print *
  endif
  ! determines group size
  call init_adios_group(myadios_group, "ValWriter")
  group_size_inc = 0
  call define_adios_scalar(myadios_group, group_size_inc, '',trim(reg_name)//"nspec", nspec)
  call define_adios_scalar(myadios_group, group_size_inc, '',trim(reg_name)//"nglob", nglob)
  do iker = 1, nker
    kernel_name = kernel_names(iker)
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    !!! warning different conventions
    write(tmp_str,'(a,a)')trim(kernel_name),'_crust_mantle'
    !write(tmp_str,'(a,a)')trim(reg_name),trim(kernel_name)
    call define_adios_global_array1D(myadios_group, group_size_inc,local_dim, '', &
                                     trim(tmp_str), so(:, :, :, :))
  enddo
  ! opens output files
  call open_file_adios_write(myadios_file, myadios_group,trim(output_file), "ValWriter")
  call set_adios_group_size(myadios_file,group_size_inc)
  ! writes to file
  call write_adios_scalar(myadios_file, myadios_group, trim(reg_name)//"nspec", nspec)
  call write_adios_scalar(myadios_file, myadios_group, trim(reg_name)//"nglob", nglob)
#endif


  ! 3. Smooth kernels / models
  !    solves inverse filtering problem for each kernel / model
  !         must solve A A s = M m
  !         where A = (M + K)
  !    done with two conjugate gradients
  is_hess = .false.
  is_kernel = .false.

  do iker = 1, nker
    !! Read input kernels
    kernel_name = kernel_names(iker)

    ! determines if parameter name is for a kernel
    is_kernel = .false.
    if (len_trim(kernel_name) > 3) then
      if (kernel_name(len_trim(kernel_name)-2:len_trim(kernel_name)) == '_kl') then
        is_kernel = .true.
      endif
    endif
    is_hess = .false.
    if (len_trim(kernel_name) > 5) then
      if (kernel_name(1:5) == 'hess_') then
        is_hess = .true.
      endif
    endif

#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS single file opening
    ! user output
    if (myrank == 0) then
      print *, 'reading in ADIOS input file : ',trim(input_file)
    endif
    call init_adios_group(myadios_val_group, "ValReader")
    call open_file_adios_read(myadios_val_file, myadios_val_group, trim(input_file))
    ! ADIOS array name
    if (is_kernel) then
      ! NOTE: reg1 => crust_mantle, others are not implemented
      varname = trim(kernel_name) // "_crust_mantle"
    else
      varname = trim(reg_name) // trim(kernel_name)
    endif
    ! user output
    if (myrank == 0) then
      print *, '  data: ADIOS ',trim(kernel_name), " is_kernel = ", is_kernel, " is_hess = ", is_hess
      print *, '  data: ADIOS using array name = ',trim(varname)
      print *
    endif
    ! reads kernel values
    call read_adios_array(myadios_val_file, myadios_val_group, myrank, nspec, trim(varname), mo(:, :, :, :))
    call synchronize_all()
    call close_file_adios_read(myadios_val_file)
    call delete_adios_group(myadios_val_group, "ValReader")
#else
    ! user output
    if (myrank == 0) then
      print *,'  kernel ',iker,'out of ',nker
      print *,'  reading data file: proc = ',myrank,' name = ',trim(kernel_name)
      print *
    endif
    ! data file
    write(local_data_file,'(a,i6.6,a)') &
          trim(input_dir)//'/proc',myrank,trim(reg_name)//trim(kernel_name)//'.bin'

    open(IIN,file=trim(local_data_file),status='old',action='read',form='unformatted',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening data file: ',trim(local_data_file)
      call exit_mpi(myrank,'Error opening data file')
    endif
    read(IIN) mo(:,:,:,:)
    close(IIN)
#endif
    ! Get normalization factor
    if (is_hess) then
      mo = sqrt(abs(mo))
    endif
    norm_ker = 0.
    norm_kerl = sum(mo**2)/nglob
    call sum_all_cr(norm_kerl, norm_ker)
    norm_ker = sqrt(norm_ker)
    call bcast_all_singlecr(norm_ker)
    mo = mo / norm_ker
    call synchronize_all()

    !! First pass
    call apply_mass_matrix_gll(mo, m)
    call solve_laplace_linear_system_cg(m, s) ! Solve A y = m

    !! Second pass
    call apply_mass_matrix_glob(s, m)
    call solve_laplace_linear_system_cg(m, s) ! Solve A s = y

    !! Save kernels
    call model_glob_to_gll(s, so)
    so = so * norm_ker
    if (is_hess) then
      so = so * so
    endif

    ! smoothed kernel file name
#ifdef USE_ADIOS_INSTEAD_OF_MESH
    ! ADIOS
    write(tmp_str,'(a,a)')trim(kernel_name),'_crust_mantle'
    !write(tmp_str,'(a,a)')trim(reg_name),trim(kernel_name)
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    call write_adios_global_1d_array(myadios_file, myadios_group, myrank, sizeprocs_adios, local_dim, &
                                     trim(tmp_str), so(:, :, :, :))
#else
    write(output_file,'(a,i6.6,a)') trim(output_dir)//'/proc', myrank, trim(reg_name)//trim(kernel_name)//'_smooth.bin'
    open(IOUT,file=trim(output_file),status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0) call exit_mpi(myrank,'Error opening smoothed kernel file')
    ! Note: output the following instead of kernel_smooth(:,:,:,1:NSPEC_AB) to create files of the same sizes
    write(IOUT) so(:,:,:,:)
    close(IOUT)
#endif
    ! user output
    if (myrank == 0) then
      print *
      print *,'written: ',trim(output_file)
      print *
    endif
    ! synchronizes
    call synchronize_all()
  enddo

#ifdef USE_ADIOS_INSTEAD_OF_MESH
  ! ADIOS
  ! actual writing
  ! (note: done at the very end, otherwise it will reset path and we would need
  ! to re-initiate a group)
  call write_adios_perform(myadios_file)
  ! closes adios files
  call close_file_adios(myadios_file)
  call finalize_adios()
#endif

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


contains

!-------------------------------------------------------------------------------------------------

  subroutine solve_laplace_linear_system_cg(m, s)

    implicit none

    integer :: iter_cg
    real(kind=CUSTOM_REAL)    :: res_ini, res_norm, res_norm_new, res_norm_inf
    real(kind=CUSTOM_REAL)    :: pAp, alpha, beta

    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in)    :: m
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: s
    real(kind=CUSTOM_REAL), dimension(:), allocatable                :: r, p, Ap


    ! Initializations and allocations
    if (.not. allocated(r))  allocate(r(nglob))
    if (.not. allocated(p))  allocate(p(nglob))
    if (.not. allocated(Ap)) allocate(Ap(nglob))
    r(:)  = 0.
    p(:)  = 0.
    Ap(:) = 0.
    !! niter_cg_max = nglob !300 !nglob   ! could be changed, here it's the theoretical maximum value
                                          ! warning should be bcast from main...

    ! Compute As and form residuals, use them as fir direction
    ! As = 0 because we assume s=0
    s(:) = 0.
    r(:) = m(:)
    p(:) = r(:)

    ! Compute initial squared norm of residuals
    call compute_scalar_product(r, r, res_norm)
    res_ini = res_norm
    if (myrank == 0) then
        write(*,*) 'Initial residual: ',res_ini
    endif
    call synchronize_all()

    ! Conjugate gradient iterations
    do iter_cg = 1, niter_cg_max

       ! Compute Ap product and scalar product
       call compute_Ax_product(p, Ap)
       call compute_scalar_product(p, Ap, pAp)

       ! Compute alpha, then update smooth model and residual
       alpha = res_norm / pAp
       s     = s + real(alpha, kind=4) * p
       r     = r - real(alpha, kind=4) * Ap

       ! Compute new norm of residual and check convergence
       call compute_scalar_product(r, r, res_norm_new)
       call compute_infinite_norm(r, res_norm_inf)
       !if ( res_norm_inf < conv_crit) then
       if ( res_norm_new < res_ini*conv_crit) then
          exit
       endif

       ! Compute beta and update conjugate direction
       beta     = res_norm_new / res_norm
       res_norm = res_norm_new
       p        = r + real(beta, kind=4) * p

       ! Print infos
       if (myrank == 0) then
           write(*,*) 'Iterations ',iter_cg,', max residual ',res_norm_inf,' l2 norm ',res_norm_new
       endif

    enddo

    if (myrank == 0) write(*,*) 'Iterations ',iter_cg,', max residual ',res_norm_inf


  end subroutine solve_laplace_linear_system_cg

!-------------------------------------------------------------------------------------------------

  subroutine compute_Ax_product(s, As)

    implicit none

    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: sl, stif, mass
    real(kind=CUSTOM_REAL), dimension(NGLLX, NGLLY, NGLLZ) :: grad_xsil, grad_etal, grad_gaml

    integer :: i, j, k, l, iel, idof

    real(kind=CUSTOM_REAL)    :: lapla_x, lapla_y, lapla_z, jacobianl
    real(kind=CUSTOM_REAL)    :: dxsi_dxl, deta_dxl, dgam_dxl
    real(kind=CUSTOM_REAL)    :: dxsi_dyl, deta_dyl, dgam_dyl
    real(kind=CUSTOM_REAL)    :: dxsi_dzl, deta_dzl, dgam_dzl
    real(kind=CUSTOM_REAL)    :: dsl_dxsi, dsl_deta, dsl_dgam
    real(kind=CUSTOM_REAL)    :: dsl_dxl, dsl_dyl, dsl_dzl

    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in)    :: s
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: As

    As(:) = 0

    ! Loop over elements
    do iel = 1, nspec

       ! Get element values
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx
                idof      = ibool(i, j, k, iel)
                sl(i,j,k) = s(idof)
             enddo
          enddo
       enddo

       ! First loop
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx

                ! Compute derivatives at elemental level
                dsl_dxsi = 0.
                dsl_deta = 0.
                dsl_dgam = 0.
                do l = 1,ngllx
                   dsl_dxsi = dsl_dxsi + sl(l, j, k) * dlagx(i, l)
                   dsl_deta = dsl_deta + sl(i, l, k) * dlagy(j, l)
                   dsl_dgam = dsl_dgam + sl(i, j, l) * dlagz(k, l)
                enddo

                ! Get derivatives informations
                dxsi_dxl  = dxsi_dx(i, j, k, iel)
                dxsi_dyl  = dxsi_dy(i, j, k, iel)
                dxsi_dzl  = dxsi_dz(i, j, k, iel)
                deta_dxl  = deta_dx(i, j, k, iel)
                deta_dyl  = deta_dy(i, j, k, iel)
                deta_dzl  = deta_dz(i, j, k, iel)
                dgam_dxl  = dgam_dx(i, j, k, iel)
                dgam_dyl  = dgam_dy(i, j, k, iel)
                dgam_dzl  = dgam_dz(i, j, k, iel)
                jacobianl = jacobian(i, j, k, iel)

                ! Get physical derivatives
                dsl_dxl = dsl_dxsi * dxsi_dxl + dsl_deta * deta_dxl + dsl_dgam * dgam_dxl
                dsl_dyl = dsl_dxsi * dxsi_dyl + dsl_deta * deta_dyl + dsl_dgam * dgam_dyl
                dsl_dzl = dsl_dxsi * dxsi_dzl + dsl_deta * deta_dzl + dsl_dgam * dgam_dzl

                ! Start product
                grad_xsil(i, j, k) = jacobianl * (dxsi_dxl * dsl_dxl + dxsi_dyl * dsl_dyl + dxsi_dzl * dsl_dzl)
                grad_etal(i, j, k) = jacobianl * (deta_dxl * dsl_dxl + deta_dyl * dsl_dyl + deta_dzl * dsl_dzl)
                grad_gaml(i, j, k) = jacobianl * (dgam_dxl * dsl_dxl + dgam_dyl * dsl_dyl + dgam_dzl * dsl_dzl)

                ! Compute mass term
                mass(i, j, k) = wxyz(i, j, k) * jacobianl *  sl(i, j, k)

             enddo
          enddo
       enddo

       ! Second loop
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx

                ! Compute derivatives at elemental level
                lapla_x = 0.
                lapla_y = 0.
                lapla_z = 0.
                do l = 1,ngllx
                   lapla_x = lapla_x + grad_xsil(l, j, k) * dlagxwx(l, i)
                   lapla_y = lapla_y + grad_etal(i, l, k) * dlagywy(l, j)
                   lapla_z = lapla_z + grad_gaml(i, j, l) * dlagzwz(l, k)
                enddo

                ! Stiffness
                stif (i,j,k) = wyz(j,k)*lapla_x + wxz(i,k)*lapla_y + wxy(i,j)*lapla_z

             enddo
          enddo
       enddo

      ! Go back to dof
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx
                idof     = ibool(i, j, k, iel)
                As(idof) = As(idof) + stif (i, j, k) + mass(i, j, k)
             enddo
          enddo
       enddo

    enddo


    call assemble_MPI_scalar(sizeprocs, nglob, As, &
                             num_interfaces, max_nibool_interfaces, &
                             nibool_interfaces, ibool_interfaces, my_neighbors)


  end subroutine compute_Ax_product

!-------------------------------------------------------------------------------------------------

  ! subroutine model_gll_to_glob(mgll, mglob)

  !   implicit none

  !   integer :: i, j, k, iel, idof

  !   real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(in)    :: mgll
  !   real(kind=CUSTOM_REAL), dimension(:),       allocatable, intent(inout) :: mglob

  !   mglob(:) = 0.
  !   do iel = 1, nspec
  !      do k = 1, ngllz
  !         do j = 1, nglly
  !            do i = 1, ngllx
  !               idof        = ibool(i, j, k, iel)
  !               mglob(idof) = mglob(idof) + mgll(i, j, k, iel) / valglob(idof)
  !            enddo
  !         enddo
  !      enddo
  !   enddo

  !   call assemble_MPI_scalar(sizeprocs, nglob_ab, mglob, &
  !                            num_interfaces, max_nibool_interfaces, &
  !                            nibool_interfaces, ibool_interfaces, my_neighbors)

  ! end subroutine model_gll_to_glob

  subroutine model_glob_to_gll(mglob, mgll)

    implicit none

    integer :: i, j, k, iel, idof

    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(inout) :: mgll
    real(kind=CUSTOM_REAL), dimension(:),       allocatable, intent(inout)    :: mglob

    do iel = 1, nspec
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx
                idof               = ibool(i, j, k, iel)
                mgll(i, j, k, iel) =  mglob(idof)
             enddo
          enddo
       enddo
    enddo

  end subroutine model_glob_to_gll

!-------------------------------------------------------------------------------------------------

  subroutine apply_mass_matrix_glob(m, s)

    implicit none

    integer :: i, j, k, iel, idof

    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in)    :: m
    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(inout) :: s

    s(:) = 0.
    do iel = 1, nspec
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx
                idof    = ibool(i, j, k, iel)
                s(idof) = s(idof) + m(idof) * jacobian(i, j, k, iel) * wxyz(i,j,k)
             enddo
          enddo
       enddo
    enddo

    call assemble_MPI_scalar(sizeprocs, nglob, s, &
                             num_interfaces, max_nibool_interfaces, &
                             nibool_interfaces, ibool_interfaces, my_neighbors)

  end subroutine apply_mass_matrix_glob

!-------------------------------------------------------------------------------------------------

  subroutine apply_mass_matrix_gll(m, s)

    implicit none

    integer :: i, j, k, iel, idof

    real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable, intent(in)    :: m
    real(kind=CUSTOM_REAL), dimension(:),       allocatable, intent(inout) :: s

    s(:) = 0.
    do iel = 1, nspec
       do k = 1, ngllz
          do j = 1, nglly
             do i = 1, ngllx
                idof    = ibool(i, j, k, iel)
                s(idof) = s(idof) + m(i, j, k, iel) * jacobian(i, j, k, iel) * wxyz(i,j,k)
                !s(idof) = m(i, j, k, iel) * jacobian(i, j, k, iel) * wxyz(i,j,k)
             enddo
          enddo
       enddo
    enddo

    call assemble_MPI_scalar(sizeprocs, nglob, s, &
                             num_interfaces, max_nibool_interfaces, &
                             nibool_interfaces, ibool_interfaces, my_neighbors)

  end subroutine apply_mass_matrix_gll

!-------------------------------------------------------------------------------------------------

  subroutine compute_scalar_product(a, b, val)

    implicit none

    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in) :: a, b
    real(kind=CUSTOM_REAL), intent(out) :: val
    real(kind=CUSTOM_REAL) :: valloc
    integer :: i
    valloc = 0
    val    = 0
    if (size(a) /= size(b)) then
       stop
    endif
    do i = 1,size(a)
       valloc = valloc + a(i) * b(i)
    enddo
    call sum_all_cr(valloc, val)
    call bcast_all_singlecr(val)
    !print *,'scal prod',myrank,valloc,val

  end subroutine compute_scalar_product

!-------------------------------------------------------------------------------------------------

  subroutine compute_infinite_norm(a, val)

    implicit none

    real(kind=CUSTOM_REAL), dimension(:), allocatable, intent(in) :: a
    integer :: i
    real(kind=CUSTOM_REAL), intent(inout) :: val
    real(kind=CUSTOM_REAL) :: valloc

    valloc = 0
    val    = 0
    do i = 1,size(a)
       valloc = max(valloc,abs(a(i)))
    enddo
    call max_all_cr(valloc, val)
    call bcast_all_singlecr(val)
    !print *,'inf norm ',myrank,valloc,val

  end subroutine compute_infinite_norm

end program smooth_laplacian_sem
