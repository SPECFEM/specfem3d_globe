!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine read_mesh_databases()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  include 'mpif.h'
  
  ! local parameters
  integer :: ier

  ! get MPI starting time
  time_start = MPI_WTIME()

  ! allocates temporary arrays
  allocate( is_on_a_slice_edge_crust_mantle(NSPEC_CRUST_MANTLE), &
           is_on_a_slice_edge_inner_core(NSPEC_INNER_CORE), &
           is_on_a_slice_edge_outer_core(NSPEC_OUTER_CORE), &
           stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating temporary is_on_a_slice_edge arrays')

  ! start reading the databases
  ! read arrays created by the mesher

  ! reads "solver_data_1.bin" & "solver_data_2.bin" files for crust and mantle
  call read_mesh_databases_CM()
  
  ! reads "solver_data_1.bin" & "solver_data_2.bin" files for outer core
  call read_mesh_databases_OC()
  
  ! reads "solver_data_1.bin" & "solver_data_2.bin" files for inner core
  call read_mesh_databases_IC()

  ! reads "boundary.bin" files to couple mantle with outer core and inner core boundaries
  call read_mesh_databases_coupling()
  
  ! reads "addressing.txt" 2-D addressing for summation between slices with MPI
  call read_mesh_databases_addressing()

  ! reads "iboolleft_..txt", "iboolright_..txt" (and "list_messages_..txt", "buffer_...txt") files and sets up MPI buffers
  call read_mesh_databases_MPIbuffers()

  ! sets up MPI interfaces
  call read_mesh_databases_MPIinter()

  ! sets up inner/outer element arrays
  call read_mesh_databases_InnerOuter()

  ! absorbing boundaries
  if(ABSORBING_CONDITIONS) then
    ! reads "stacey.bin" files and sets up arrays for Stacey conditions
    call read_mesh_databases_stacey()
  endif

  ! user output
  if( myrank == 0 ) then
    ! elapsed time since beginning of mesh generation
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for reading mesh in seconds = ',sngl(tCPU)
    write(IMAIN,*)
  endif
  
  ! frees temporary allocated arrays
  deallocate(is_on_a_slice_edge_crust_mantle, &
            is_on_a_slice_edge_outer_core, &
            is_on_a_slice_edge_inner_core)

  end subroutine read_mesh_databases

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_CM()

! mesh for CRUST MANTLE region

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
  integer :: nspec_iso,nspec_tiso,nspec_ani
  logical :: READ_KAPPA_MU,READ_TISO

  ! dummy array that does not need to be actually read
  integer, dimension(NSPEC_CRUST_MANTLE) :: dummy_i

  ! crust and mantle
  if(ANISOTROPIC_3D_MANTLE_VAL) then
    READ_KAPPA_MU = .false.
    READ_TISO = .false.
    nspec_iso = 1
    nspec_tiso = 1
    nspec_ani = NSPEC_CRUST_MANTLE
  else
    nspec_iso = NSPEC_CRUST_MANTLE
    if(TRANSVERSE_ISOTROPY_VAL) then
      nspec_tiso = NSPECMAX_TISO_MANTLE
    else
      nspec_tiso = 1
    endif
    nspec_ani = 1
    READ_KAPPA_MU = .true.
    READ_TISO = .true.
  endif
  
  call read_arrays_solver(IREGION_CRUST_MANTLE,myrank, &
            rho_vp_crust_mantle,rho_vs_crust_mantle, &
            xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
            kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
            nspec_iso,nspec_tiso,nspec_ani, &
            c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
            c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
            c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
            c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
            c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
            c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
            c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
            ibool_crust_mantle,dummy_i, &
            ispec_is_tiso_crust_mantle, &
            is_on_a_slice_edge_crust_mantle,rmass_crust_mantle,rmass_ocean_load, &
            NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  ! check that the number of points in this slice is correct
  if(minval(ibool_crust_mantle(:,:,:,:)) /= 1 .or. &
    maxval(ibool_crust_mantle(:,:,:,:)) /= NGLOB_CRUST_MANTLE) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in crust and mantle')


  end subroutine read_mesh_databases_CM

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_OC()

! mesh for OUTER CORE region

  use specfem_par
  use specfem_par_outercore
  implicit none

  ! local parameters
  integer :: nspec_iso,nspec_tiso,nspec_ani
  logical :: READ_KAPPA_MU,READ_TISO
  
  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array

  logical, dimension(:),allocatable:: dummy_ispec_is_tiso
  integer, dimension(:),allocatable :: dummy_idoubling_outer_core
  
  ! outer core (no anisotropy nor S velocity)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the outer core
  READ_KAPPA_MU = .false.
  READ_TISO = .false.
  nspec_iso = NSPEC_OUTER_CORE
  nspec_tiso = 1
  nspec_ani = 1
  
  ! dummy allocation
  allocate(dummy_ispec_is_tiso(NSPEC_OUTER_CORE))
  allocate(dummy_idoubling_outer_core(NSPEC_OUTER_CORE))

  call read_arrays_solver(IREGION_OUTER_CORE,myrank, &
            vp_outer_core,dummy_array, &
            xstore_outer_core,ystore_outer_core,zstore_outer_core, &
            xix_outer_core,xiy_outer_core,xiz_outer_core, &
            etax_outer_core,etay_outer_core,etaz_outer_core, &
            gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
            rhostore_outer_core,kappavstore_outer_core,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            nspec_iso,nspec_tiso,nspec_ani, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            ibool_outer_core,dummy_idoubling_outer_core,dummy_ispec_is_tiso, &
            is_on_a_slice_edge_outer_core,rmass_outer_core,rmass_ocean_load, &
            NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  deallocate(dummy_idoubling_outer_core)
  deallocate(dummy_ispec_is_tiso)

  ! check that the number of points in this slice is correct
  if(minval(ibool_outer_core(:,:,:,:)) /= 1 .or. &
     maxval(ibool_outer_core(:,:,:,:)) /= NGLOB_OUTER_CORE) &
    call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in outer core')

  end subroutine read_mesh_databases_OC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_IC()

! mesh for INNER CORE region

  use specfem_par
  use specfem_par_innercore
  implicit none

  ! local parameters
  integer :: nspec_iso,nspec_tiso,nspec_ani
  logical :: READ_KAPPA_MU,READ_TISO

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  logical, dimension(:),allocatable:: dummy_ispec_is_tiso
  
  ! inner core (no anisotropy)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the inner core
  READ_KAPPA_MU = .true.
  READ_TISO = .false.
  nspec_iso = NSPEC_INNER_CORE
  nspec_tiso = 1
  if(ANISOTROPIC_INNER_CORE_VAL) then
    nspec_ani = NSPEC_INNER_CORE
  else
    nspec_ani = 1
  endif
  
  ! dummy allocation
  allocate(dummy_ispec_is_tiso(NSPEC_INNER_CORE))
  
  call read_arrays_solver(IREGION_INNER_CORE,myrank, &
            dummy_array,dummy_array, &
            xstore_inner_core,ystore_inner_core,zstore_inner_core, &
            xix_inner_core,xiy_inner_core,xiz_inner_core, &
            etax_inner_core,etay_inner_core,etaz_inner_core, &
            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
            dummy_array,dummy_array,dummy_array, &
            nspec_iso,nspec_tiso,nspec_ani, &
            c11store_inner_core,c12store_inner_core,c13store_inner_core, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,c33store_inner_core, &
            dummy_array,dummy_array,dummy_array, &
            c44store_inner_core,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            ibool_inner_core,idoubling_inner_core,dummy_ispec_is_tiso, &
            is_on_a_slice_edge_inner_core,rmass_inner_core,rmass_ocean_load, &
            NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  deallocate(dummy_ispec_is_tiso)

  ! check that the number of points in this slice is correct
  if(minval(ibool_inner_core(:,:,:,:)) /= 1 .or. maxval(ibool_inner_core(:,:,:,:)) /= NGLOB_INNER_CORE) &
    call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in inner core')

  end subroutine read_mesh_databases_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_coupling()

! to couple mantle with outer core

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  include 'mpif.h'
  
  ! local parameters
  integer njunk1,njunk2,njunk3

  ! crust and mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! Stacey put back
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read')
  read(27) nspec2D_xmin_crust_mantle
  read(27) nspec2D_xmax_crust_mantle
  read(27) nspec2D_ymin_crust_mantle
  read(27) nspec2D_ymax_crust_mantle
  read(27) njunk1
  read(27) njunk2

! boundary parameters
  read(27) ibelm_xmin_crust_mantle
  read(27) ibelm_xmax_crust_mantle
  read(27) ibelm_ymin_crust_mantle
  read(27) ibelm_ymax_crust_mantle
  read(27) ibelm_bottom_crust_mantle
  read(27) ibelm_top_crust_mantle

  read(27) normal_xmin_crust_mantle
  read(27) normal_xmax_crust_mantle
  read(27) normal_ymin_crust_mantle
  read(27) normal_ymax_crust_mantle
  read(27) normal_bottom_crust_mantle
  read(27) normal_top_crust_mantle

  read(27) jacobian2D_xmin_crust_mantle
  read(27) jacobian2D_xmax_crust_mantle
  read(27) jacobian2D_ymin_crust_mantle
  read(27) jacobian2D_ymax_crust_mantle
  read(27) jacobian2D_bottom_crust_mantle
  read(27) jacobian2D_top_crust_mantle
  close(27)


  ! read parameters to couple fluid and solid regions
  !
  ! outer core

  ! create name of database
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  ! boundary parameters

  ! Stacey put back
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read')
  read(27) nspec2D_xmin_outer_core
  read(27) nspec2D_xmax_outer_core
  read(27) nspec2D_ymin_outer_core
  read(27) nspec2D_ymax_outer_core
  read(27) njunk1
  read(27) njunk2

  nspec2D_zmin_outer_core = NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

  read(27) ibelm_xmin_outer_core
  read(27) ibelm_xmax_outer_core
  read(27) ibelm_ymin_outer_core
  read(27) ibelm_ymax_outer_core
  read(27) ibelm_bottom_outer_core
  read(27) ibelm_top_outer_core

  read(27) normal_xmin_outer_core
  read(27) normal_xmax_outer_core
  read(27) normal_ymin_outer_core
  read(27) normal_ymax_outer_core
  read(27) normal_bottom_outer_core
  read(27) normal_top_outer_core

  read(27) jacobian2D_xmin_outer_core
  read(27) jacobian2D_xmax_outer_core
  read(27) jacobian2D_ymin_outer_core
  read(27) jacobian2D_ymax_outer_core
  read(27) jacobian2D_bottom_outer_core
  read(27) jacobian2D_top_outer_core
  close(27)


  !
  ! inner core
  !

  ! create name of database
  call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_PATH)

  ! read info for vertical edges for central cube matching in inner core
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read')
  read(27) nspec2D_xmin_inner_core
  read(27) nspec2D_xmax_inner_core
  read(27) nspec2D_ymin_inner_core
  read(27) nspec2D_ymax_inner_core
  read(27) njunk1
  read(27) njunk2

  ! boundary parameters
  read(27) ibelm_xmin_inner_core
  read(27) ibelm_xmax_inner_core
  read(27) ibelm_ymin_inner_core
  read(27) ibelm_ymax_inner_core
  read(27) ibelm_bottom_inner_core
  read(27) ibelm_top_inner_core
  close(27)


  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then

    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

    open(unit=27,file=prname(1:len_trim(prname))//'boundary_disc.bin', &
          status='old',form='unformatted',action='read')
    read(27) njunk1,njunk2,njunk3
    if (njunk1 /= NSPEC2D_MOHO .and. njunk2 /= NSPEC2D_400 .and. njunk3 /= NSPEC2D_670) &
               call exit_mpi(myrank, 'Error reading ibelm_disc.bin file')
    read(27) ibelm_moho_top
    read(27) ibelm_moho_bot
    read(27) ibelm_400_top
    read(27) ibelm_400_bot
    read(27) ibelm_670_top
    read(27) ibelm_670_bot
    read(27) normal_moho
    read(27) normal_400
    read(27) normal_670
    close(27)

    k_top = 1
    k_bot = NGLLZ

    ! initialization
    moho_kl = 0.; d400_kl = 0.; d670_kl = 0.; cmb_kl = 0.; icb_kl = 0.
  endif

  end subroutine read_mesh_databases_coupling

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_mesh_databases_addressing()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  include 'mpif.h'

  ! local parameters
  integer :: ier,iproc,iproc_read
  integer :: NUM_FACES,NPROC_ONE_DIRECTION

  ! open file with global slice number addressing
  if(myrank == 0) then
    open(unit=IIN,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening addressing.txt')
    do iproc = 0,NPROCTOT_VAL-1
      read(IIN,*) iproc_read,ichunk,iproc_xi,iproc_eta
      if(iproc_read /= iproc) call exit_MPI(myrank,'incorrect slice number read')
      addressing(ichunk,iproc_xi,iproc_eta) = iproc
      ichunk_slice(iproc) = ichunk
      iproc_xi_slice(iproc) = iproc_xi
      iproc_eta_slice(iproc) = iproc_eta
    enddo
    close(IIN)
  endif

  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(addressing,NCHUNKS_VAL*NPROC_XI_VAL*NPROC_ETA_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ichunk_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iproc_xi_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iproc_eta_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! output a topology map of slices - fix 20x by nproc
  if (myrank == 0 .and. NCHUNKS_VAL == 6) then
    write(IMAIN,*) 'Spatial distribution of the slices'
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(20x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AB
        write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(1x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_BC
        write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(3x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AC
        write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(3x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_BC_ANTIPODE
        write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(20x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AB_ANTIPODE
        write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(20x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AC_ANTIPODE
        write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
  endif

  ! determine chunk number and local slice coordinates using addressing
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  ! define maximum size for message buffers
  ! use number of elements found in the mantle since it is the largest region
  NGLOB2DMAX_XY = max(NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE))

  ! number of corners and faces shared between chunks and number of message types
  if(NCHUNKS_VAL == 1 .or. NCHUNKS_VAL == 2) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 1
  else if(NCHUNKS_VAL == 3) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 3
  else if(NCHUNKS_VAL == 6) then
    NCORNERSCHUNKS = 8
    NUM_FACES = 4
    NUM_MSG_TYPES = 3
  else
    call exit_MPI(myrank,'number of chunks must be either 1, 2, 3 or 6')
  endif
  ! if more than one chunk then same number of processors in each direction
  NPROC_ONE_DIRECTION = NPROC_XI_VAL
  ! total number of messages corresponding to these common faces
  NUMMSGS_FACES = NPROC_ONE_DIRECTION*NUM_FACES*NUM_MSG_TYPES

  ! debug checks with compiled value
  !if( NUMMSGS_FACES /= NUMMSGS_FACES_VAL ) then
  !  print*,'check: NUMMSGS_FACES',NUMMSGS_FACES,NUMMSGS_FACES_VAL
  !  call exit_mpi(myrank,'error NUMMSGS_FACES_VAL, please recompile solver')
  !endif

  end subroutine read_mesh_databases_addressing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPIbuffers()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore  
  implicit none

  ! local parameters
  integer :: ier
  character(len=150) :: filename

  ! read 2-D addressing for summation between slices with MPI

  ! mantle and crust
  if(myrank == 0) write(IMAIN,*) 'crust/mantle region:'
  
  call read_arrays_buffers_solver(IREGION_CRUST_MANTLE,myrank,iboolleft_xi_crust_mantle, &
     iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
     npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
     iboolfaces_crust_mantle,npoin2D_faces_crust_mantle, &
     iboolcorner_crust_mantle, &
     NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE), &
     NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_XY,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,LOCAL_PATH,NCHUNKS_VAL)

  ! outer core
  if(myrank == 0) write(IMAIN,*) 'outer core region:'

  call read_arrays_buffers_solver(IREGION_OUTER_CORE,myrank, &
     iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
     npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
     iboolfaces_outer_core,npoin2D_faces_outer_core, &
     iboolcorner_outer_core, &
     NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE), &
     NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE),NGLOB2DMAX_XY,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,LOCAL_PATH,NCHUNKS_VAL)

  ! inner core
  if(myrank == 0) write(IMAIN,*) 'inner core region:'
  
  call read_arrays_buffers_solver(IREGION_INNER_CORE,myrank, &
     iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
     npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
     iboolfaces_inner_core,npoin2D_faces_inner_core, &
     iboolcorner_inner_core, &
     NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE), &
     NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE),NGLOB2DMAX_XY,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,LOCAL_PATH,NCHUNKS_VAL)

  ! added this to reduce the size of the buffers
  ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_crust_mantle(:) + npoin2D_xi_inner_core(:)), &
                        maxval(npoin2D_eta_crust_mantle(:) + npoin2D_eta_inner_core(:)))

  allocate(buffer_send_faces(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED), &
          buffer_received_faces(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating mpi buffer')

  allocate(b_buffer_send_faces(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED), &
          b_buffer_received_faces(NDIM,npoin2D_max_all_CM_IC,NUMFACES_SHARED),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating mpi b_buffer')

  ! central cube buffers
  if(INCLUDE_CENTRAL_CUBE) then

    if(myrank == 0) write(IMAIN,*) 'including central cube'

    ! compute number of messages to expect in cube as well as their size
    call comp_central_cube_buffer_size(iproc_xi,iproc_eta,ichunk, &
                NPROC_XI_VAL,NPROC_ETA_VAL,NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                nb_msgs_theor_in_cube,npoin2D_cube_from_slices)

    ! this value is used for dynamic memory allocation, therefore make sure it is never zero
    if(nb_msgs_theor_in_cube > 0) then
      non_zero_nb_msgs_theor_in_cube = nb_msgs_theor_in_cube
    else
      non_zero_nb_msgs_theor_in_cube = 1
    endif

    ! allocate buffers for cube and slices
    allocate(sender_from_slices_to_cube(non_zero_nb_msgs_theor_in_cube), &
            buffer_all_cube_from_slices(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM), &
            b_buffer_all_cube_from_slices(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM), &
            buffer_slices(npoin2D_cube_from_slices,NDIM), &
            b_buffer_slices(npoin2D_cube_from_slices,NDIM), &
            buffer_slices2(npoin2D_cube_from_slices,NDIM), &
            ibool_central_cube(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating cube buffers')

    ! handles the communications with the central cube if it was included in the mesh
    ! create buffers to assemble with the central cube
    call create_central_cube_buffers(myrank,iproc_xi,iproc_eta,ichunk, &
               NPROC_XI_VAL,NPROC_ETA_VAL,NCHUNKS_VAL,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
               NSPEC2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NSPEC2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
               NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
               addressing,ibool_inner_core,idoubling_inner_core, &
               xstore_inner_core,ystore_inner_core,zstore_inner_core, &
               nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
               nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
               ibelm_xmin_inner_core,ibelm_xmax_inner_core, &
               ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
               nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices, &
               receiver_cube_from_slices,sender_from_slices_to_cube,ibool_central_cube, &
               buffer_slices,buffer_slices2,buffer_all_cube_from_slices)

    if(myrank == 0) write(IMAIN,*) ''

  else

    ! allocate fictitious buffers for cube and slices with a dummy size
    ! just to be able to use them as arguments in subroutine calls
    allocate(sender_from_slices_to_cube(1), &
            buffer_all_cube_from_slices(1,1,1), &
            b_buffer_all_cube_from_slices(1,1,1), &
            buffer_slices(1,1), &
            b_buffer_slices(1,1), &
            buffer_slices2(1,1), &
            ibool_central_cube(1,1),stat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error allocating dummy buffers')

  endif

  ! note: fix_... routines below update is_on_a_slice_edge_.. arrays:
  !          assign flags for each element which is on a rim of the slice
  !          thus, they include elements on top and bottom not shared with other MPI partitions  
  !
  !          we will re-set these flags when setting up inner/outer elements, but will
  !          use these arrays for now as initial guess for the search for elements which share a global point 
  !          between different MPI processes  
  call fix_non_blocking_slices(is_on_a_slice_edge_crust_mantle,iboolright_xi_crust_mantle, &
         iboolleft_xi_crust_mantle,iboolright_eta_crust_mantle,iboolleft_eta_crust_mantle, &
         npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle,ibool_crust_mantle, &
         mask_ibool,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,NGLOB2DMAX_XMIN_XMAX_CM,NGLOB2DMAX_YMIN_YMAX_CM)

  call fix_non_blocking_slices(is_on_a_slice_edge_outer_core,iboolright_xi_outer_core, &
         iboolleft_xi_outer_core,iboolright_eta_outer_core,iboolleft_eta_outer_core, &
         npoin2D_xi_outer_core,npoin2D_eta_outer_core,ibool_outer_core, &
         mask_ibool,NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,NGLOB2DMAX_XMIN_XMAX_OC,NGLOB2DMAX_YMIN_YMAX_OC)

  call fix_non_blocking_slices(is_on_a_slice_edge_inner_core,iboolright_xi_inner_core, &
         iboolleft_xi_inner_core,iboolright_eta_inner_core,iboolleft_eta_inner_core, &
         npoin2D_xi_inner_core,npoin2D_eta_inner_core,ibool_inner_core, &
         mask_ibool,NSPEC_INNER_CORE,NGLOB_INNER_CORE,NGLOB2DMAX_XMIN_XMAX_IC,NGLOB2DMAX_YMIN_YMAX_IC)

  if(INCLUDE_CENTRAL_CUBE) then
    ! updates flags for elements on slice boundaries
    call fix_non_blocking_central_cube(is_on_a_slice_edge_inner_core, &
         ibool_inner_core,NSPEC_INNER_CORE,NGLOB_INNER_CORE,nb_msgs_theor_in_cube,ibelm_bottom_inner_core, &
         idoubling_inner_core,npoin2D_cube_from_slices,ibool_central_cube, &
         NSPEC2D_BOTTOM(IREGION_INNER_CORE),ichunk)
  endif
  
  ! debug: saves element flags
  ! crust mantle
  !write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_crust_mantle_proc',myrank  
  !call write_VTK_data_elem_l(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
  !                          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                          ibool_crust_mantle, &
  !                          is_on_a_slice_edge_crust_mantle,filename)  
  ! outer core
  !write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_outer_core_proc',myrank  
  !call write_VTK_data_elem_l(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
  !                          xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                          ibool_outer_core, &
  !                          is_on_a_slice_edge_outer_core,filename)
!daniel  
  ! inner core
  write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_inner_core_proc',myrank  
  call write_VTK_data_elem_l(NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                            xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                            ibool_inner_core, &
                            is_on_a_slice_edge_inner_core,filename)
  
  end subroutine read_mesh_databases_MPIbuffers
  
  
!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPIinter()

! sets up interfaces between MPI processes

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore  
  implicit none

  include 'mpif.h'
  
  ! local parameters  
  integer :: ier,ndim_assemble
  character(len=150) :: filename

  ! temporary buffers for send and receive between faces of the slices and the chunks
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC) ::  &
    buffer_send_faces_scalar,buffer_received_faces_scalar
  
  ! assigns initial maximum arrays
  ! for global slices, maximum number of neighbor is around 17 ( 8 horizontal, max of 8 on bottom )
  integer, parameter :: MAX_NEIGHBOURS = 8 + NCORNERSCHUNKS_VAL
  integer, dimension(MAX_NEIGHBOURS) :: my_neighbours,nibool_neighbours
  integer, dimension(:,:),allocatable :: ibool_neighbours
  integer :: max_nibool  
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: test_flag
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: test_flag_cc
  integer :: i,j,k,ispec,iglob
  
  ! estimates initial maximum ibool array
  max_nibool = npoin2D_max_all_CM_IC * NUMFACES_SHARED &
               + non_zero_nb_msgs_theor_in_cube*npoin2D_cube_from_slices

  allocate(ibool_neighbours(max_nibool,MAX_NEIGHBOURS), stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating ibool_neighbours')


! sets up MPI interfaces  
! crust mantle region  
  if( myrank == 0 ) write(IMAIN,*) 'crust mantle mpi:'
  allocate(test_flag(NGLOB_CRUST_MANTLE), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating test_flag')
  
  ! sets flag to rank id (+1 to avoid problems with zero rank)
  test_flag(:) = myrank + 1.0

  ! assembles values
  call assemble_MPI_scalar_block(myrank,test_flag, &
            NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_XY,NCHUNKS_VAL)
    
  ! removes own myrank id (+1)
  test_flag(:) = test_flag(:) - ( myrank + 1.0)

  ! debug: saves array
  !write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_test_flag_crust_mantle_proc',myrank
  !call write_VTK_glob_points(NGLOB_CRUST_MANTLE, &
  !                      xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                      test_flag,filename)
  
  ! determines neighbor rank for shared faces 
  call rmd_get_MPI_interfaces(myrank,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                            test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                            num_interfaces_crust_mantle,max_nibool_interfaces_crust_mantle, &
                            max_nibool,MAX_NEIGHBOURS, &
                            ibool_crust_mantle,&
                            is_on_a_slice_edge_crust_mantle, &
                            IREGION_CRUST_MANTLE,.false.)  

  deallocate(test_flag)
    
  ! stores MPI interfaces informations
  allocate(my_neighbours_crust_mantle(num_interfaces_crust_mantle), &
          nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &          
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array my_neighbours_crust_mantle etc.')

  ! copies interfaces arrays
  if( num_interfaces_crust_mantle > 0 ) then
    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_crust_mantle,num_interfaces_crust_mantle), &
           stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array ibool_interfaces_crust_mantle')

    ! ranks of neighbour processes
    my_neighbours_crust_mantle(:) = my_neighbours(1:num_interfaces_crust_mantle)
    ! number of global ibool entries on each interface
    nibool_interfaces_crust_mantle(:) = nibool_neighbours(1:num_interfaces_crust_mantle)
    ! global iglob point ids on each interface
    ibool_interfaces_crust_mantle(:,:) = ibool_neighbours(1:max_nibool_interfaces_crust_mantle,1:num_interfaces_crust_mantle)    
  else
    ! dummy allocation (fortran90 should allow allocate statement with zero array size)
    max_nibool_interfaces_crust_mantle = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
  endif
  
  ! debug: saves 1. MPI interface
  !if( num_interfaces_crust_mantle >= 1 ) then
  !  write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_1_points_crust_mantle_proc',myrank
  !  call write_VTK_data_points(NGLOB_CRUST_MANTLE, &
  !                      xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                      ibool_interfaces_crust_mantle(1:nibool_interfaces_crust_mantle(1),1), &
  !                      nibool_interfaces_crust_mantle(1),filename)
  !endif


! outer core region  
  if( myrank == 0 ) write(IMAIN,*) 'outer core mpi:'

  allocate(test_flag(NGLOB_OUTER_CORE), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating test_flag outer core')
  
  ! sets flag to rank id (+1 to avoid problems with zero rank)
  test_flag(:) = myrank + 1.0

  ! assembles values
  call assemble_MPI_scalar_block(myrank,test_flag, &
            NGLOB_OUTER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
            npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
            iboolfaces_outer_core,iboolcorner_outer_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE),NGLOB2DMAX_XY,NCHUNKS_VAL)

    
  ! removes own myrank id (+1)
  test_flag(:) = test_flag(:) - ( myrank + 1.0)

  ! debug: saves array
  !write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_test_flag_outer_core_proc',myrank
  !call write_VTK_glob_points(NGLOB_OUTER_CORE, &
  !                      xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                      test_flag,filename)

  
  ! determines neighbor rank for shared faces 
  call rmd_get_MPI_interfaces(myrank,NGLOB_OUTER_CORE,NSPEC_OUTER_CORE, &
                            test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                            num_interfaces_outer_core,max_nibool_interfaces_outer_core, &
                            max_nibool,MAX_NEIGHBOURS, &
                            ibool_outer_core,&
                            is_on_a_slice_edge_outer_core, &
                            IREGION_OUTER_CORE,.false.)  

  deallocate(test_flag)
    
  ! stores MPI interfaces informations
  allocate(my_neighbours_outer_core(num_interfaces_outer_core), &
          nibool_interfaces_outer_core(num_interfaces_outer_core), &          
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array my_neighbours_outer_core etc.')

  ! copies interfaces arrays
  if( num_interfaces_outer_core > 0 ) then
    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_outer_core,num_interfaces_outer_core), &
           stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array ibool_interfaces_outer_core')

    ! ranks of neighbour processes
    my_neighbours_outer_core(:) = my_neighbours(1:num_interfaces_outer_core)
    ! number of global ibool entries on each interface
    nibool_interfaces_outer_core(:) = nibool_neighbours(1:num_interfaces_outer_core)
    ! global iglob point ids on each interface
    ibool_interfaces_outer_core(:,:) = ibool_neighbours(1:max_nibool_interfaces_outer_core,1:num_interfaces_outer_core)    
  else
    ! dummy allocation (fortran90 should allow allocate statement with zero array size)
    max_nibool_interfaces_outer_core = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ier)
  endif

  ! debug: saves 1. MPI interface
  !if( num_interfaces_outer_core >= 1 ) then
  !  write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_1_points_outer_core_proc',myrank
  !  call write_VTK_data_points(NGLOB_OUTER_CORE, &
  !                      xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                      ibool_interfaces_outer_core(1:nibool_interfaces_outer_core(1),1), &
  !                      nibool_interfaces_outer_core(1),filename)
  !endif

! inner core
  if( myrank == 0 ) write(IMAIN,*) 'inner core mpi:'

  allocate(test_flag(NGLOB_INNER_CORE), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating test_flag inner core')
  
  ! sets flag to rank id (+1 to avoid problems with zero rank)
  test_flag(:) = 0.0
  do ispec=1,NSPEC_INNER_CORE
    ! suppress fictitious elements in central cube
    if(idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle    
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_inner_core(i,j,k,ispec)          
          test_flag(iglob) = myrank + 1.0
        enddo
      enddo
    enddo
  enddo

  ! assembles values
  call assemble_MPI_scalar_block(myrank,test_flag, &
            NGLOB_INNER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE),NGLOB2DMAX_XY,NCHUNKS_VAL)

  ! removes own myrank id (+1)
  test_flag(:) = test_flag(:) - ( myrank + 1.0)  
  where( test_flag(:) < 0.0 ) test_flag(:) = 0.0

  ! debug: saves array
  write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_test_flag_inner_core_A_proc',myrank
  call write_VTK_glob_points(NGLOB_INNER_CORE, &
                        xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                        test_flag,filename)

!  ! gets new interfaces for inner_core without central cube yet
!  ! determines neighbor rank for shared faces 
!  call rmd_get_MPI_interfaces(myrank,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
!                            test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
!                            num_interfaces_inner_core,max_nibool_interfaces_inner_core, &
!                            max_nibool,MAX_NEIGHBOURS, &
!                            ibool_inner_core,&
!                            is_on_a_slice_edge_inner_core, &
!                            IREGION_INNER_CORE,.false.,idoubling_inner_core)  

  ! including central cube
  if(INCLUDE_CENTRAL_CUBE) then
    if( myrank == 0 ) write(IMAIN,*) 'inner core with central cube mpi:'

    allocate(test_flag_cc(NGLOB_INNER_CORE), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating test_flag_cc inner core')
  
    ! re-sets flag to rank id (+1 to avoid problems with zero rank)
    test_flag_cc(:) = 0.0
    do ispec=1,NSPEC_INNER_CORE
      ! suppress fictitious elements in central cube
      if(idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle    
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)          
            test_flag_cc(iglob) = myrank + 1.0
          enddo
        enddo
      enddo
    enddo  
  
    ! test_flag is a scalar, not a vector
    ndim_assemble = 1
    ! use central cube buffers to assemble the inner core mass matrix with the central cube
    call assemble_MPI_central_cube_block(ichunk,nb_msgs_theor_in_cube, sender_from_slices_to_cube, &
                 npoin2D_cube_from_slices, buffer_all_cube_from_slices, &
                 buffer_slices, buffer_slices2, ibool_central_cube, &
                 receiver_cube_from_slices, ibool_inner_core, &
                 idoubling_inner_core, NSPEC_INNER_CORE, &
                 ibelm_bottom_inner_core, NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                 NGLOB_INNER_CORE, &
                 test_flag_cc,ndim_assemble)
                 

    ! removes own myrank id (+1)
    test_flag_cc(:) = test_flag_cc(:) - ( myrank + 1.0)  
    where( test_flag_cc(:) < 0.0 ) test_flag_cc(:) = 0.0
                 
    write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_test_flag_inner_core_B_proc',myrank
    call write_VTK_glob_points(NGLOB_INNER_CORE, &
                        xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                        test_flag_cc,filename)

!    ! adds additional inner core points
!    call rmd_get_MPI_interfaces(myrank,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
!                            test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
!                            num_interfaces_inner_core,max_nibool_interfaces_inner_core, &
!                            max_nibool,MAX_NEIGHBOURS, &
!                            ibool_inner_core,&
!                            is_on_a_slice_edge_inner_core, &
!                            IREGION_INNER_CORE,.true.,idoubling_inner_core)
     
    ! adds both together 
    test_flag(:) = test_flag(:) + test_flag_cc(:)

    deallocate(test_flag_cc)    

    ! debug: saves array
    write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_test_flag_inner_core_C_proc',myrank
    call write_VTK_glob_points(NGLOB_INNER_CORE, &
                              xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                              test_flag,filename)
  
  endif

  ! gets new interfaces for inner_core without central cube yet
  ! determines neighbor rank for shared faces 
  call rmd_get_MPI_interfaces(myrank,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                            test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                            num_interfaces_inner_core,max_nibool_interfaces_inner_core, &
                            max_nibool,MAX_NEIGHBOURS, &
                            ibool_inner_core,&
                            is_on_a_slice_edge_inner_core, &
                            IREGION_INNER_CORE,.false.,idoubling_inner_core)  
      
  deallocate(test_flag)
  
  ! stores MPI interfaces informations
  allocate(my_neighbours_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &          
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array my_neighbours_inner_core etc.')

  ! copies interfaces arrays
  if( num_interfaces_inner_core > 0 ) then
    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_inner_core,num_interfaces_inner_core), &
           stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array ibool_interfaces_inner_core')

    ! ranks of neighbour processes
    my_neighbours_inner_core(:) = my_neighbours(1:num_interfaces_inner_core)
    ! number of global ibool entries on each interface
    nibool_interfaces_inner_core(:) = nibool_neighbours(1:num_interfaces_inner_core)
    ! global iglob point ids on each interface
    ibool_interfaces_inner_core(:,:) = ibool_neighbours(1:max_nibool_interfaces_inner_core,1:num_interfaces_inner_core)    
  else
    ! dummy allocation (fortran90 should allow allocate statement with zero array size)
    max_nibool_interfaces_inner_core = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ier)
  endif

  ! debug: saves 1. MPI interface
  if( myrank == 0 .and. num_interfaces_inner_core >= 1 ) then
    write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_1_points_inner_core_proc',myrank  
    call write_VTK_data_points(NGLOB_INNER_CORE, &
                        xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                        ibool_interfaces_inner_core(1:nibool_interfaces_inner_core(1),1), &
                        nibool_interfaces_inner_core(1),filename)
    !print*,'saved: ',trim(filename)//'.vtk'
  endif  

  ! synchronizes MPI processes
  call sync_all()

  ! frees temporary array
  deallocate(ibool_neighbours)  


  ! allocates MPI buffers
  ! crust mantle
  allocate(buffer_send_vector_crust_mantle(NDIM,max_nibool_interfaces_crust_mantle,num_interfaces_crust_mantle), &
          buffer_recv_vector_crust_mantle(NDIM,max_nibool_interfaces_crust_mantle,num_interfaces_crust_mantle), &
          request_send_vector_crust_mantle(num_interfaces_crust_mantle), &
          request_recv_vector_crust_mantle(num_interfaces_crust_mantle), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array buffer_send_vector_crust_mantle etc.')
  if( SIMULATION_TYPE == 3 ) then
    allocate(b_buffer_send_vector_crust_mantle(NDIM,max_nibool_interfaces_crust_mantle,num_interfaces_crust_mantle), &
            b_buffer_recv_vector_crust_mantle(NDIM,max_nibool_interfaces_crust_mantle,num_interfaces_crust_mantle), &
            b_request_send_vector_crust_mantle(num_interfaces_crust_mantle), &
            b_request_recv_vector_crust_mantle(num_interfaces_crust_mantle), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array b_buffer_send_vector_crust_mantle etc.')
  endif

  ! outer core
  allocate(buffer_send_scalar_outer_core(max_nibool_interfaces_outer_core,num_interfaces_outer_core), &
          buffer_recv_scalar_outer_core(max_nibool_interfaces_outer_core,num_interfaces_outer_core), &
          request_send_scalar_outer_core(num_interfaces_outer_core), &
          request_recv_scalar_outer_core(num_interfaces_outer_core), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array buffer_send_vector_outer_core etc.')
  if( SIMULATION_TYPE == 3 ) then
    allocate(b_buffer_send_scalar_outer_core(max_nibool_interfaces_outer_core,num_interfaces_outer_core), &
            b_buffer_recv_scalar_outer_core(max_nibool_interfaces_outer_core,num_interfaces_outer_core), &
            b_request_send_scalar_outer_core(num_interfaces_outer_core), &
            b_request_recv_scalar_outer_core(num_interfaces_outer_core), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array b_buffer_send_vector_outer_core etc.')
  endif
  
  ! inner core
  allocate(buffer_send_vector_inner_core(NDIM,max_nibool_interfaces_inner_core,num_interfaces_inner_core), &
          buffer_recv_vector_inner_core(NDIM,max_nibool_interfaces_inner_core,num_interfaces_inner_core), &
          request_send_vector_inner_core(num_interfaces_inner_core), &
          request_recv_vector_inner_core(num_interfaces_inner_core), &
          stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array buffer_send_vector_inner_core etc.')
  if( SIMULATION_TYPE == 3 ) then
    allocate(b_buffer_send_vector_inner_core(NDIM,max_nibool_interfaces_inner_core,num_interfaces_inner_core), &
            b_buffer_recv_vector_inner_core(NDIM,max_nibool_interfaces_inner_core,num_interfaces_inner_core), &
            b_request_send_vector_inner_core(num_interfaces_inner_core), &
            b_request_recv_vector_inner_core(num_interfaces_inner_core), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array b_buffer_send_vector_inner_core etc.')
  endif
  
  end subroutine read_mesh_databases_MPIinter


!
!-------------------------------------------------------------------------------------------------
!

  subroutine rmd_get_MPI_interfaces(myrank,NGLOB,NSPEC, &
                                    test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                                    num_interfaces,max_nibool_interfaces, &
                                    max_nibool,MAX_NEIGHBOURS, &
                                    ibool,&
                                    is_on_a_slice_edge, &
                                    IREGION,add_central_cube,idoubling)

  use constants
  
  implicit none

  include 'mpif.h'

  integer,intent(in) :: myrank,NGLOB,NSPEC

  real(kind=CUSTOM_REAL),dimension(NGLOB),intent(in) :: test_flag

  integer,intent(in) :: max_nibool
  integer,intent(in) :: MAX_NEIGHBOURS
  integer, dimension(MAX_NEIGHBOURS),intent(inout) :: my_neighbours,nibool_neighbours
  integer, dimension(max_nibool,MAX_NEIGHBOURS),intent(inout) :: ibool_neighbours

  integer,intent(out) :: num_interfaces,max_nibool_interfaces

  integer,dimension(NGLLX,NGLLY,NGLLZ,NSPEC),intent(in) :: ibool
  
  logical,dimension(NSPEC),intent(inout) :: is_on_a_slice_edge
  
  integer,intent(in) :: IREGION  
  logical,intent(in) :: add_central_cube
  integer,dimension(NSPEC),optional:: idoubling
  
  ! local parameters  
  integer :: ispec,iglob,j,k
  integer :: iface,iedge,icorner
  integer :: ii,iinterface,icurrent,rank
  integer :: npoin
  logical :: is_done,ispec_is_outer
  integer,dimension(NGLOB) :: work_test_flag
  logical,dimension(NSPEC) :: work_ispec_is_outer
  
  ! initializes 
  if( add_central_cube) then
    ! adds points to existing inner_core interfaces
    iinterface = num_interfaces
    work_ispec_is_outer(:) = is_on_a_slice_edge(:)    
  else
    ! creates new interfaces
    iinterface = 0
    num_interfaces = 0
    max_nibool_interfaces = 0
    my_neighbours(:) = -1
    nibool_neighbours(:) = 0
    ibool_neighbours(:,:) = 0
    work_ispec_is_outer(:) = .false.
  endif
  
  ! makes working copy (converted to nearest integers)
  work_test_flag(:) = nint( test_flag(:) )

  ! loops over all elements
  do ispec = 1,NSPEC

    ! exclude elements in inner part of slice
    !if( .not. is_on_a_slice_edge(ispec) ) cycle

    ! exclude elements in fictitious core
    if( IREGION == IREGION_INNER_CORE) then
      if( idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE ) cycle
    endif

    ! sets flag if element has global points shared with other processes    
    ispec_is_outer = .false.    

    ! 1. finds neighbours which share a whole face with this process
    ! (faces are shared only with 1 other neighbour process)
        
    ! loops over all faces of element
    do iface = 1, 6
    
      ! chooses a point inside face
      select case( iface )
      case( 1 )
        ! face I == 1    
        iglob = ibool(1,2,2,ispec) 
      case( 2 )
        ! face I == NGLLX
        iglob = ibool(NGLLX,2,2,ispec) 
      case( 3 )
        ! face J == 1
        iglob = ibool(2,1,2,ispec) 
      case( 4 )
        ! face J == NGLLY
        iglob = ibool(2,NGLLY,2,ispec) 
      case( 5 )
        ! face K == 1
        iglob = ibool(2,2,1,ispec)
      case( 6 )
        ! face K == NGLLZ
        iglob = ibool(2,2,NGLLZ,ispec)
      end select      

      ! checks assembled flag on global point
      if( work_test_flag(iglob) > 0 ) then
        ispec_is_outer = .true.
        
        ! rank of neighbor process
        rank = work_test_flag(iglob) - 1
        
        ! checks ranks range
        if( rank < 0 .or. rank >= NPROCTOT_VAL ) then
          print*,'error face rank: ',myrank,'ispec=',ispec
          print*,'  neighbor rank = ',rank,'exceeds total nproc:',NPROCTOT_VAL
          print*,'  face ',iface
          call exit_mpi(myrank,'error face neighbor mpi rank')
        endif
        
        ! checks if already stored
        icurrent = 0
        is_done = .false.
        do ii = 1,iinterface
          if( rank == my_neighbours(ii) ) then 
            icurrent = ii
            is_done = .true.
            exit
          endif
        enddo
        
        ! updates interfaces array
        if( .not. is_done ) then
          iinterface = iinterface + 1
          if( iinterface > MAX_NEIGHBOURS ) &
            call exit_mpi(myrank,'interface face exceeds MAX_NEIGHBOURS range')
          ! adds as neighbor new interface
          my_neighbours(iinterface) = rank
          icurrent = iinterface
        endif
        if( icurrent == 0 ) &
          call exit_mpi(myrank,'could not find current interface for this neighbor, please check my_neighbours')
        
        ! adds interface points and removes neighbor flag from face
        ! assumes NGLLX == NGLLY == NGLLZ
        do k=1,NGLLX
          do j=1,NGLLX
            select case( iface )
            case( 1 )
              ! face I == 1    
              iglob = ibool(1,j,k,ispec) 
            case( 2 )
              ! face I == NGLLX
              iglob = ibool(NGLLX,j,k,ispec) 
            case( 3 )
              ! face J == 1
              iglob = ibool(j,1,k,ispec) 
            case( 4 )
              ! face J == NGLLY
              iglob = ibool(j,NGLLY,k,ispec) 
            case( 5 )
              ! face K == 1
              iglob = ibool(j,k,1,ispec)
            case( 6 )
              ! face K == NGLLZ
              iglob = ibool(j,k,NGLLZ,ispec)
            end select      
            
            ! checks that we take each global point (on edges and corners) only once
            if( work_test_flag(iglob) <= 0 ) cycle ! continues to next point
            
            ! increases number of total points on this interface
            nibool_neighbours(icurrent) = nibool_neighbours(icurrent) + 1
            if( nibool_neighbours(icurrent) > max_nibool) &
              call exit_mpi(myrank,'interface face exceeds max_nibool range')
              
            ! stores interface iglob index
            ibool_neighbours( nibool_neighbours(icurrent),icurrent ) = iglob
            
            ! re-sets flag
            work_test_flag(iglob) = work_test_flag(iglob) - ( rank + 1 )
            ! debug
            if( work_test_flag(iglob) < 0 ) then
              print*,'error face flag:',myrank,'ispec=',ispec,'rank=',rank
              print*,'  flag=',work_test_flag(iglob),'iface jk=',iface,j,k
              call exit_mpi(myrank,'error face flag')
            endif
          enddo
        enddo
      endif    
    enddo ! iface

    ! 2. finds neighbours which share a single edge with this process
    ! note: by now, faces have subtracted their neighbours, edges can hold only one more process info

    ! loops over all edges of element
    do iedge = 1, 12
    
      ! chooses a point inside edge but not corner
      select case( iedge )
      case( 1 )
        ! face I == 1, J == 1
        iglob = ibool(1,1,2,ispec) 
      case( 2 )
        ! face I == 1, J == NGLLY
        iglob = ibool(1,NGLLY,2,ispec)
      case( 3 )
        ! face I == 1, K == 1
        iglob = ibool(1,2,1,ispec)
      case( 4 )
        ! face I == 1, K == NGLLZ
        iglob = ibool(1,2,NGLLZ,ispec)
      case( 5 )
        ! face I == NGLLX, J == 1
        iglob = ibool(NGLLX,1,2,ispec)
      case( 6 )
        ! face I == NGLLX, J == NGLLY
        iglob = ibool(NGLLX,NGLLY,2,ispec)
      case( 7 )
        ! face I == NGLLX, K == 1
        iglob = ibool(NGLLX,2,1,ispec)
      case( 8 )
        ! face I == NGLLX, K == NGLLZ
        iglob = ibool(NGLLX,2,NGLLZ,ispec)
      case( 9 )
        ! face J == 1, K == 1
        iglob = ibool(2,1,1,ispec)
      case( 10 )
        ! face J == 1, K == NGLLZ
        iglob = ibool(2,1,NGLLZ,ispec)
      case( 11 )
        ! face J == NGLLY, K == 1
        iglob = ibool(2,NGLLY,1,ispec)
      case( 12 )
        ! face J == NGLLY, K == NGLLZ
        iglob = ibool(2,NGLLY,NGLLZ,ispec)
      end select

      ! checks assembled flag on global point
      if( work_test_flag(iglob) > 0 ) then
        ispec_is_outer = .true.
        
        ! rank of neighbor process
        rank = work_test_flag(iglob) - 1
        
        ! checks ranks range
        if( rank < 0 .or. rank >= NPROCTOT_VAL ) then
          print*,'error egde rank: ',myrank
          print*,'  neighbor rank = ',rank,'exceeds total nproc:',NPROCTOT_VAL
          print*,'  edge ',iedge
          call exit_mpi(myrank,'error edge neighbor mpi rank')
        endif
        
        ! checks if already stored
        icurrent = 0
        is_done = .false.
        do ii = 1,iinterface
          if( rank == my_neighbours(ii) ) then 
            icurrent = ii
            is_done = .true.
            exit
          endif
        enddo
        
        ! updates interfaces array
        if( .not. is_done ) then
          iinterface = iinterface + 1
          if( iinterface > MAX_NEIGHBOURS ) &
            call exit_mpi(myrank,'interface edge exceeds MAX_NEIGHBOURS range')
          ! adds as neighbor new interface
          my_neighbours(iinterface) = rank
          icurrent = iinterface
        endif
        if( icurrent == 0 ) &
          call exit_mpi(myrank,'could not find current interface for this neighbor, please check my_neighbours')
        
        ! adds interface points and removes neighbor flag from edge
        ! assumes NGLLX == NGLLY == NGLLZ
        do k = 1,NGLLX
          select case( iedge )
          case( 1 )
            ! face I == 1, J == 1
            iglob = ibool(1,1,k,ispec) 
          case( 2 )
            ! face I == 1, J == NGLLY
            iglob = ibool(1,NGLLY,k,ispec)
          case( 3 )
            ! face I == 1, K == 1
            iglob = ibool(1,k,1,ispec)
          case( 4 )
            ! face I == 1, K == NGLLZ
            iglob = ibool(1,k,NGLLZ,ispec)
          case( 5 )
            ! face I == NGLLX, J == 1
            iglob = ibool(NGLLX,1,k,ispec)
          case( 6 )
            ! face I == NGLLX, J == NGLLY
            iglob = ibool(NGLLX,NGLLY,k,ispec)
          case( 7 )
            ! face I == NGLLX, K == 1
            iglob = ibool(NGLLX,k,1,ispec)
          case( 8 )
            ! face I == NGLLX, K == NGLLZ
            iglob = ibool(NGLLX,k,NGLLZ,ispec)
          case( 9 )
            ! face J == 1, K == 1
            iglob = ibool(k,1,1,ispec)
          case( 10 )
            ! face J == 1, K == NGLLZ
            iglob = ibool(k,1,NGLLZ,ispec)
          case( 11 )
            ! face J == NGLLY, K == 1
            iglob = ibool(k,NGLLY,1,ispec)
          case( 12 )
            ! face J == NGLLY, K == NGLLZ
            iglob = ibool(k,NGLLY,NGLLZ,ispec)
          end select

          ! checks that we take each global point (on edges and corners) only once
          if( work_test_flag(iglob) <= 0 ) cycle ! continues to next point

          ! increases number of total points on this interface
          nibool_neighbours(icurrent) = nibool_neighbours(icurrent) + 1
          if( nibool_neighbours(icurrent) > max_nibool) &
            call exit_mpi(myrank,'interface edge exceeds max_nibool range')
            
          ! stores interface iglob index          
          ibool_neighbours( nibool_neighbours(icurrent),icurrent ) = iglob
            
          ! re-sets flag
          work_test_flag(iglob) = work_test_flag(iglob) - ( rank + 1 )

          ! debug
          if( work_test_flag(iglob) < 0 ) call exit_mpi(myrank,'error edge flag')
          
        enddo        
      endif    
    enddo ! iedge


    ! 3. finds neighbours which share a single corner with this process
    ! note: faces and edges have subtracted their neighbors, only one more process left possible

    ! loops over all corners of element
    do icorner = 1, 8
    
      ! chooses a corner point
      select case( icorner )
      case( 1 )
        ! face I == 1
        iglob = ibool(1,1,1,ispec)
      case( 2 )
        ! face I == 1
        iglob = ibool(1,NGLLY,1,ispec)
      case( 3 )
        ! face I == 1
        iglob = ibool(1,1,NGLLZ,ispec)
      case( 4 )
        ! face I == 1
        iglob = ibool(1,NGLLY,NGLLZ,ispec)
      case( 5 )
        ! face I == NGLLX
        iglob = ibool(NGLLX,1,1,ispec)
      case( 6 )
        ! face I == NGLLX
        iglob = ibool(NGLLX,NGLLY,1,ispec)
      case( 7 )
        ! face I == NGLLX
        iglob = ibool(NGLLX,1,NGLLZ,ispec)
      case( 8 )
        ! face I == NGLLX
        iglob = ibool(NGLLX,NGLLY,NGLLZ,ispec)
      end select

      ! makes sure that all elements on mpi interfaces are included
      ! uses original test_flag array, since the working copy reduces values
      ! note: there can be elements which have an edge or corner shared with
      !          other mpi partitions, but have the work_test_flag value already set to zero
      !          since the iglob point was found before.
      !          also, this check here would suffice to determine the outer flag, but we also include the 
      !          check everywhere we encounter it too
      if( test_flag(iglob) > 0.5 ) then
        ispec_is_outer = .true.
      endif

      ! checks assembled flag on global point
      if( work_test_flag(iglob) > 0 ) then
        ispec_is_outer = .true.
        
        ! rank of neighbor process
        rank = work_test_flag(iglob) - 1
        
        ! checks ranks range
        if( rank < 0 .or. rank >= NPROCTOT_VAL ) then
          print*,'error corner: ',myrank
          print*,'  neighbor rank = ',rank,'exceeds total nproc:',NPROCTOT_VAL
          print*,'  corner ',icorner
          call exit_mpi(myrank,'error corner neighbor mpi rank')
        endif
        
        ! checks if already stored
        icurrent = 0
        is_done = .false.
        do ii = 1,iinterface
          if( rank == my_neighbours(ii) ) then 
            icurrent = ii
            is_done = .true.
            exit
          endif
        enddo
        
        ! updates interfaces array
        if( .not. is_done ) then
          iinterface = iinterface + 1
          if( iinterface > MAX_NEIGHBOURS ) &
            call exit_mpi(myrank,'interface corner exceed MAX_NEIGHBOURS range')
          ! adds as neighbor new interface
          my_neighbours(iinterface) = rank
          icurrent = iinterface
        endif
        if( icurrent == 0 ) &
          call exit_mpi(myrank,'could not find current interface for this neighbor, please check my_neighbours')
        
        ! adds this corner as interface point and removes neighbor flag from face
        ! increases number of total points on this interface
        nibool_neighbours(icurrent) = nibool_neighbours(icurrent) + 1
        if( nibool_neighbours(icurrent) > max_nibool) &
          call exit_mpi(myrank,'interface corner exceeds max_nibool range')
          
        ! stores interface iglob index        
        ibool_neighbours( nibool_neighbours(icurrent),icurrent ) = iglob
            
        ! re-sets flag
        work_test_flag(iglob) = work_test_flag(iglob) - ( rank + 1 )

        ! debug
        if( work_test_flag(iglob) < 0 ) call exit_mpi(myrank,'error corner flag')
        
      endif
      
    enddo ! icorner
        
    ! stores flags for outer elements when recognized as such 
    ! (inner/outer elements separated for non-blocking mpi communications)
    if( ispec_is_outer ) then
      work_ispec_is_outer(ispec) = .true.
    endif
    
  enddo 

  ! number of outer elements (on MPI interfaces)
  npoin = count( work_ispec_is_outer )
    
  ! debug: user output
  if( myrank == 0 ) then
    write(IMAIN,*) '  interfaces : ',iinterface
    write(IMAIN,*) '  my_neighbours: ',my_neighbours(1:iinterface)  
    write(IMAIN,*) '  nibool_neighbours: ',nibool_neighbours(1:iinterface)
    write(IMAIN,*) '  test flag min/max: ',minval(work_test_flag),maxval(work_test_flag)
    write(IMAIN,*) '  outer elements: ',npoin
    write(IMAIN,*)
  endif
  call sync_all()

  ! checks if all points were recognized
  if( maxval(work_test_flag) > 0 ) then
    print*,'error mpi interface rank: ',myrank
    print*,'  work_test_flag min/max :',minval(work_test_flag),maxval(work_test_flag)
    call exit_mpi(myrank,'error: mpi points remain unrecognized, please check mesh interfaces')
  endif

  ! checks if all points were taken only once
  if( minval(work_test_flag) < 0 ) then
    print*,'error mpi interface rank: ',myrank
    print*,'  work_test_flag min/max :',minval(work_test_flag),maxval(work_test_flag)
    call exit_mpi(myrank,'error: mpi points counted more than once, please check mesh interfaces')
  endif

  ! sets interfaces infos
  num_interfaces = iinterface
  max_nibool_interfaces = maxval( nibool_neighbours(1:num_interfaces) )

  ! optional: ibool usually is already sorted,
  !                this makes sure ibool_neighbours arrays are still sorted 
  !               (iglob indices in increasing order; we will access acceleration fields accel(:,iglob), 
  !                thus it helps if iglob strides are short and accesses are close-by)
  do iinterface = 1,num_interfaces
    npoin = nibool_neighbours(iinterface)    
    call heap_sort( npoin, ibool_neighbours(1:npoin,iinterface) )
    
    ! debug: checks if unique set of iglob values
    do j=1,npoin-1
      if( ibool_neighbours(j,iinterface) == ibool_neighbours(j+1,iinterface) ) then
        print*,'error mpi interface rank:',myrank
        print*,'  interface: ',my_neighbours(iinterface),'point: ',j,'of',npoin
        call exit_mpi(myrank,'error: mpi points not unique on interface')
      endif
    enddo
  enddo


  ! re-sets flags for outer elements
  is_on_a_slice_edge(:) = work_ispec_is_outer(:)

  end subroutine rmd_get_MPI_interfaces

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_InnerOuter()

! sets up inner/outer elements for non-blocking MPI communication

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore  
  implicit none
  
  ! local parameters
  real :: percentage_edge
  integer :: ier,ispec,iinner,iouter
  character(len=150) :: filename

  ! stores inner / outer elements
  !
  ! note: arrays is_on_a_slice_edge_.. have flags set for elements which need to
  !         communicate with other MPI processes
  
  ! crust_mantle
  nspec_outer_crust_mantle = count( is_on_a_slice_edge_crust_mantle )
  nspec_inner_crust_mantle = NSPEC_CRUST_MANTLE - nspec_outer_crust_mantle
  
  num_phase_ispec_crust_mantle = max(nspec_inner_crust_mantle,nspec_outer_crust_mantle)
  
  allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array phase_ispec_inner_crust_mantle')
  
  phase_ispec_inner_crust_mantle(:,:) = 0
  iinner = 0
  iouter = 0
  do ispec=1,NSPEC_CRUST_MANTLE
    if( is_on_a_slice_edge_crust_mantle(ispec) ) then
      ! outer element
      iouter = iouter + 1
      phase_ispec_inner_crust_mantle(iouter,1) = ispec
    else
      ! inner element
      iinner = iinner + 1
      phase_ispec_inner_crust_mantle(iinner,2) = ispec
    endif
  enddo

  ! outer_core
  nspec_outer_outer_core = count( is_on_a_slice_edge_outer_core )
  nspec_inner_outer_core = NSPEC_OUTER_CORE - nspec_outer_outer_core

  num_phase_ispec_outer_core = max(nspec_inner_outer_core,nspec_outer_outer_core)
  
  allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array phase_ispec_inner_outer_core')
  
  phase_ispec_inner_outer_core(:,:) = 0
  iinner = 0
  iouter = 0
  do ispec=1,NSPEC_OUTER_CORE
    if( is_on_a_slice_edge_outer_core(ispec) ) then
      ! outer element
      iouter = iouter + 1
      phase_ispec_inner_outer_core(iouter,1) = ispec
    else
      ! inner element
      iinner = iinner + 1
      phase_ispec_inner_outer_core(iinner,2) = ispec
    endif
  enddo

  ! inner_core
  nspec_outer_inner_core = count( is_on_a_slice_edge_inner_core )
  nspec_inner_inner_core = NSPEC_INNER_CORE - nspec_outer_inner_core

  num_phase_ispec_inner_core = max(nspec_inner_inner_core,nspec_outer_inner_core)
  
  allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating array phase_ispec_inner_inner_core')
  
  phase_ispec_inner_inner_core(:,:) = 0
  iinner = 0
  iouter = 0
  do ispec=1,NSPEC_INNER_CORE
    if( is_on_a_slice_edge_inner_core(ispec) ) then
      ! outer element
      iouter = iouter + 1
      phase_ispec_inner_inner_core(iouter,1) = ispec
    else
      ! inner element
      iinner = iinner + 1
      phase_ispec_inner_inner_core(iinner,2) = ispec
    endif
  enddo
  
  ! user output
  if(myrank == 0) then

    write(IMAIN,*) 'for overlapping of communications with calculations:'
    write(IMAIN,*)

    percentage_edge = 100.0 - 100. * nspec_outer_crust_mantle / real(NSPEC_CRUST_MANTLE)
    write(IMAIN,*) 'percentage of edge elements in crust/mantle ',percentage_edge,'%'
    write(IMAIN,*) 'percentage of volume elements in crust/mantle ',100. - percentage_edge,'%'
    write(IMAIN,*)

    percentage_edge = 100.0 - 100.* nspec_outer_outer_core / real(NSPEC_OUTER_CORE)
    write(IMAIN,*) 'percentage of edge elements in outer core ',percentage_edge,'%'
    write(IMAIN,*) 'percentage of volume elements in outer core ',100. - percentage_edge,'%'
    write(IMAIN,*)

    percentage_edge = 100.0 - 100. * nspec_outer_inner_core / real(NSPEC_INNER_CORE)
    write(IMAIN,*) 'percentage of edge elements in inner core ',percentage_edge,'%'
    write(IMAIN,*) 'percentage of volume elements in inner core ',100. - percentage_edge,'%'
    write(IMAIN,*)

  endif

  ! debug: saves element flags
  ! crust mantle
  !write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_crust_mantle_proc',myrank  
  !call write_VTK_data_elem_l(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
  !                          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
  !                          ibool_crust_mantle, &
  !                          is_on_a_slice_edge_crust_mantle,filename)  
  ! outer core
  !write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_outer_core_proc',myrank  
  !call write_VTK_data_elem_l(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
  !                          xstore_outer_core,ystore_outer_core,zstore_outer_core, &
  !                          ibool_outer_core, &
  !                          is_on_a_slice_edge_outer_core,filename)
  ! inner core
  write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_innerouter_inner_core_proc',myrank  
  call write_VTK_data_elem_l(NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                            xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                            ibool_inner_core, &
                            is_on_a_slice_edge_inner_core,filename)

  end subroutine read_mesh_databases_InnerOuter


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_stacey()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer(kind=8) :: filesize
  integer :: ier
  integer :: nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm
  integer :: nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc,nabs_zmin_oc

  ! sets up absorbing boundary buffer arrays

  ! crust_mantle
  if (nspec2D_xmin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmin_cm = nspec2D_xmin_crust_mantle
  else
    nabs_xmin_cm = 1
  endif
  !daniel: not sure why dimensions are (..,..,..,..,8)  ?
  !allocate(absorb_xmin_crust_mantle5(NDIM,NGLLY,NGLLZ,nabs_xmin_cm,8),stat=ier)
  allocate(absorb_xmin_crust_mantle(NDIM,NGLLY,NGLLZ,nabs_xmin_cm),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb xmin')

  if (nspec2D_xmax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmax_cm = nspec2D_xmax_crust_mantle
  else
    nabs_xmax_cm = 1
  endif
  !daniel: not sure why dimensions are (..,..,..,..,8)
  !allocate(absorb_xmax_crust_mantle5(NDIM,NGLLY,NGLLZ,nabs_xmax_cm,8),stat=ier)
  allocate(absorb_xmax_crust_mantle(NDIM,NGLLY,NGLLZ,nabs_xmax_cm),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb xmax')

  if (nspec2D_ymin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymin_cm = nspec2D_ymin_crust_mantle
  else
    nabs_ymin_cm = 1
  endif
  !daniel: not sure why dimensions are (..,..,..,..,8)
  !allocate(absorb_ymin_crust_mantle5(NDIM,NGLLX,NGLLZ,nabs_ymin_cm,8),stat=ier)
  allocate(absorb_ymin_crust_mantle(NDIM,NGLLX,NGLLZ,nabs_ymin_cm),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb ymin')

  if (nspec2D_ymax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymax_cm = nspec2D_ymax_crust_mantle
  else
    nabs_ymax_cm = 1
  endif
  !daniel: not sure why dimensions are (..,..,..,..,8)
  !allocate(absorb_ymax_crust_mantle5(NDIM,NGLLX,NGLLZ,nabs_ymax_cm,8),stat=ier)
  allocate(absorb_ymax_crust_mantle(NDIM,NGLLX,NGLLZ,nabs_ymax_cm),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb ymax')

  ! outer_core
  if (nspec2D_xmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmin_oc = nspec2D_xmin_outer_core
  else
    nabs_xmin_oc = 1
  endif
  allocate(absorb_xmin_outer_core(NGLLY,NGLLZ,nabs_xmin_oc),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb xmin')

  if (nspec2D_xmax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_xmax_oc = nspec2D_xmax_outer_core
  else
    nabs_xmax_oc = 1
  endif
  allocate(absorb_xmax_outer_core(NGLLY,NGLLZ,nabs_xmax_oc),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb xmax')

  if (nspec2D_ymin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymin_oc = nspec2D_ymin_outer_core
  else
    nabs_ymin_oc = 1
  endif
  allocate(absorb_ymin_outer_core(NGLLX,NGLLZ,nabs_ymin_oc),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb ymin')

  if (nspec2D_ymax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_ymax_oc = nspec2D_ymax_outer_core
  else
    nabs_ymax_oc = 1
  endif
  allocate(absorb_ymax_outer_core(NGLLX,NGLLZ,nabs_ymax_oc),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb ymax')

  if (nspec2D_zmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    nabs_zmin_oc = nspec2D_zmin_outer_core
  else
    nabs_zmin_oc = 1
  endif
  allocate(absorb_zmin_outer_core(NGLLX,NGLLY,nabs_zmin_oc),stat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error allocating absorb zmin')


  ! crust and mantle
  
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! read arrays for Stacey conditions
  open(unit=27,file=prname(1:len_trim(prname))//'stacey.bin', &
        status='old',form='unformatted',action='read')
  read(27) nimin_crust_mantle
  read(27) nimax_crust_mantle
  read(27) njmin_crust_mantle
  read(27) njmax_crust_mantle
  read(27) nkmin_xi_crust_mantle
  read(27) nkmin_eta_crust_mantle
  close(27)

  if (nspec2D_xmin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmin_crust_mantle)

    ! total file size
    filesize = reclen_xmin_crust_mantle
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=51,file=trim(prname)//'absorb_xmin.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmin_crust_mantle+2*4)
!    else
!      open(unit=51,file=trim(prname)//'absorb_xmin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_xmin_crust_mantle+2*4)

      call open_file_abs_r(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_xmax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmax_crust_mantle)

    ! total file size
    filesize = reclen_xmax_crust_mantle
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=52,file=trim(prname)//'absorb_xmax.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmax_crust_mantle+2*4)
!    else
!      open(unit=52,file=trim(prname)//'absorb_xmax.bin', &
!            status='unknown',form='unformatted',access='direct', &
!            recl=reclen_xmax_crust_mantle+2*4)

      call open_file_abs_r(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_ymin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymin_crust_mantle)

    ! total file size
    filesize = reclen_ymin_crust_mantle
    filesize = filesize*NSTEP


    if (SIMULATION_TYPE == 3) then
!      open(unit=53,file=trim(prname)//'absorb_ymin.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymin_crust_mantle+2*4)
!    else
!      open(unit=53,file=trim(prname)//'absorb_ymin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymin_crust_mantle+2*4)

      call open_file_abs_r(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_ymax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymax_crust_mantle)

    ! total file size
    filesize = reclen_ymax_crust_mantle
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=54,file=trim(prname)//'absorb_ymax.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymax_crust_mantle+2*4)
!    else
!      open(unit=54,file=trim(prname)//'absorb_ymax.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymax_crust_mantle+2*4)

      call open_file_abs_r(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif


  ! outer core
  
  ! create name of database
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  ! read arrays for Stacey conditions
  open(unit=27,file=prname(1:len_trim(prname))//'stacey.bin', &
        status='old',form='unformatted',action='read')
  read(27) nimin_outer_core
  read(27) nimax_outer_core
  read(27) njmin_outer_core
  read(27) njmax_outer_core
  read(27) nkmin_xi_outer_core
  read(27) nkmin_eta_outer_core
  close(27)

  if (nspec2D_xmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmin_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmin_outer_core)

    ! total file size
    filesize = reclen_xmin_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=61,file=trim(prname)//'absorb_xmin.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmin_outer_core+2*4)
!    else
!      open(unit=61,file=trim(prname)//'absorb_xmin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_xmin_outer_core+2*4)

      call open_file_abs_r(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_xmax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmax_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmax_outer_core)

    ! total file size
    filesize = reclen_xmax_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=62,file=trim(prname)//'absorb_xmax.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmax_outer_core+2*4)
!    else
!      open(unit=62,file=trim(prname)//'absorb_xmax.bin', &
!            status='unknown',form='unformatted',access='direct', &
!            recl=reclen_xmax_outer_core+2*4)

      call open_file_abs_r(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
   endif

  endif

  if (nspec2D_ymin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymin_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymin_outer_core)

    ! total file size
    filesize = reclen_ymin_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=63,file=trim(prname)//'absorb_ymin.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymin_outer_core+2*4)
!    else
!      open(unit=63,file=trim(prname)//'absorb_ymin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymin_outer_core+2*4)

      call open_file_abs_r(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)

    endif
  endif

  if (nspec2D_ymax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymax_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymax_outer_core)

    ! total file size
    filesize = reclen_ymax_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=64,file=trim(prname)//'absorb_ymax.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymax_outer_core+2*4)
!    else
!      open(unit=64,file=trim(prname)//'absorb_ymax.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymax_outer_core+2*4)

      call open_file_abs_r(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)

    endif
  endif

  if (nspec2D_zmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)))then

    ! size of single record
    reclen_zmin = CUSTOM_REAL * (NGLLX * NGLLY * nspec2D_zmin_outer_core)

    ! total file size
    filesize = reclen_zmin
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=65,file=trim(prname)//'absorb_zmin.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_zmin+2*4)
!    else
!      open(unit=65,file=trim(prname)//'absorb_zmin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_zmin+2*4)

      call open_file_abs_r(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    endif
  endif

  end subroutine read_mesh_databases_stacey


!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!

!daniel: TODO - place this auxiliary function...

  subroutine heap_sort( N, array )

! heap sort algorithm
! sorts integer array (in increasing order, like 1 - 5 - 6 - 9 - 12 - 13 - 14 -...)
  
  implicit none
  integer,intent(in) :: N
  integer,dimension(N),intent(inout) :: array
  
  ! local parameters
  integer :: tmp
  integer :: i
  
  ! checks if anything to do
  if( N < 2 ) return
  
  ! builds heap
  do i = N/2, 1, -1  
    call heap_sort_siftdown(N,array,i,N)
  enddo
    
  ! sorts array
  do i = N, 2, -1
    ! swaps last and first entry in this section
    tmp = array(1)
    array(1) = array(i)
    array(i) = tmp
    call heap_sort_siftdown(N,array,1,i-1)
  enddo
  
  end subroutine heap_sort

!
!----
!  

  subroutine heap_sort_siftdown(N,array,start,bottom)

  implicit none
  
  integer,intent(in):: N
  integer,dimension(N),intent(inout) :: array
  integer :: start,bottom
  
  ! local parameters
  integer :: i,j
  integer :: tmp
  
  i = start
  tmp = array(i)  
  j = 2*i
  do while( j <= bottom )
    ! chooses larger value first in this section
    if( j < bottom ) then
      if( array(j) <= array(j+1) ) j = j + 1
    endif
    
    ! checks if section already smaller than inital value
    if( array(j) < tmp ) exit
        
    array(i) = array(j)
    i = j
    j = 2*i
  enddo
  
  array(i) = tmp      
  return
  
  end subroutine heap_sort_siftdown

!
!-------------------------------------------------------------------------------------------------
!

!daniel: TODO - place this auxiliary function...

! external mesh routine for saving vtk files for points locations

  subroutine write_VTK_data_points(nglob, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  points_globalindices,num_points_globalindices, &
                                  prname_file)

  implicit none

  include "constants.h"

  integer :: nglob

! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! gll data values array
  integer :: num_points_globalindices
  integer, dimension(num_points_globalindices) :: points_globalindices

! file name
  character(len=150) prname_file

  integer :: i,iglob

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', num_points_globalindices, ' float'
  do i=1,num_points_globalindices
    iglob = points_globalindices(i)
    if( iglob <= 0 .or. iglob > nglob ) then
      print*,'error: '//prname_file(1:len_trim(prname_file))//'.vtk'
      print*,'error global index: ',iglob,i
      close(IOVTK)
      stop 'error vtk points file'
    endif

    write(IOVTK,'(3e18.6)') xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOVTK,*) ""
  
  close(IOVTK)

  end subroutine write_VTK_data_points

!
!-------------------------------------------------------------------------------------------------
!

! external mesh routine for saving vtk files for points locations

  subroutine write_VTK_glob_points(nglob, &
                                  xstore_dummy,ystore_dummy,zstore_dummy, &
                                  glob_values, &
                                  prname_file)

  implicit none

  include "constants.h"

  integer :: nglob

  ! global coordinates
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

  ! gll data values array
  real(kind=CUSTOM_REAL), dimension(nglob) :: glob_values

  ! file name
  character(len=150) prname_file

  integer :: iglob

  ! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do iglob=1,nglob
    write(IOVTK,*) xstore_dummy(iglob),ystore_dummy(iglob),zstore_dummy(iglob)
  enddo
  write(IOVTK,*) ""

  ! writes out gll-data (velocity) for each element point
  write(IOVTK,'(a,i12)') "POINT_DATA ",nglob
  write(IOVTK,'(a)') "SCALARS glob_data float"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do iglob=1,nglob
    write(IOVTK,*) glob_values(iglob)
  enddo
  write(IOVTK,*) ""
  
  close(IOVTK)

  end subroutine write_VTK_glob_points



!
!-------------------------------------------------------------------------------------------------
!

! routine for saving vtk file holding logical flag on each spectral element

  subroutine write_VTK_data_elem_l(nspec,nglob, &
                        xstore_dummy,ystore_dummy,zstore_dummy,ibool, &
                        elem_flag,prname_file)


  implicit none

  include "constants.h"

  integer :: nspec,nglob

! global coordinates
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore_dummy,ystore_dummy,zstore_dummy

! element flag array
  logical, dimension(nspec) :: elem_flag
  integer :: ispec,i

! file name
  character(len=150) prname_file

! write source and receiver VTK files for Paraview
  !debug
  !write(IMAIN,*) '  vtk file: '
  !write(IMAIN,*) '    ',prname_file(1:len_trim(prname_file))//'.vtk'

  open(IOVTK,file=prname_file(1:len_trim(prname_file))//'.vtk',status='unknown')
  write(IOVTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOVTK,'(a)') 'material model VTK file'
  write(IOVTK,'(a)') 'ASCII'
  write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOVTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i=1,nglob
    write(IOVTK,'(3e18.6)') xstore_dummy(i),ystore_dummy(i),zstore_dummy(i)
  enddo
  write(IOVTK,*) ""

  ! note: indices for vtk start at 0
  write(IOVTK,'(a,i12,i12)') "CELLS ",nspec,nspec*9
  do ispec=1,nspec
    write(IOVTK,'(9i12)') 8,ibool(1,1,1,ispec)-1,ibool(NGLLX,1,1,ispec)-1,ibool(NGLLX,NGLLY,1,ispec)-1,ibool(1,NGLLY,1,ispec)-1,&
          ibool(1,1,NGLLZ,ispec)-1,ibool(NGLLX,1,NGLLZ,ispec)-1,ibool(NGLLX,NGLLY,NGLLZ,ispec)-1,ibool(1,NGLLY,NGLLZ,ispec)-1
  enddo
  write(IOVTK,*) ""

  ! type: hexahedrons
  write(IOVTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOVTK,*) (12,ispec=1,nspec)
  write(IOVTK,*) ""

  write(IOVTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOVTK,'(a)') "SCALARS elem_flag integer"
  write(IOVTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    if( elem_flag(ispec) .eqv. .true. ) then
      write(IOVTK,*) 1
    else
      write(IOVTK,*) 0
    endif
  enddo
  write(IOVTK,*) ""
  close(IOVTK)


  end subroutine write_VTK_data_elem_l

