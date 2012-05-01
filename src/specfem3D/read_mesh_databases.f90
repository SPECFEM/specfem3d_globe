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

  ! get MPI starting time
  time_start = MPI_WTIME()

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

  ! sets up MPI interfaces, inner/outer elements and mesh coloring
  call read_mesh_databases_MPI()

  ! absorbing boundaries
  if(ABSORBING_CONDITIONS) then
    ! reads "stacey.bin" files and sets up arrays for Stacey conditions
    call read_mesh_databases_stacey()
  endif

  ! user output
  call sync_all()
  if( myrank == 0 ) then
    ! elapsed time since beginning of mesh generation
    tCPU = MPI_WTIME() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for reading mesh in seconds = ',sngl(tCPU)
    write(IMAIN,*)
  endif

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
  integer, dimension(:),allocatable :: dummy_i
  integer :: ier

  ! allocates dummy array
  allocate(dummy_i(NSPEC_CRUST_MANTLE),stat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error allocating dummy_i array in read_mesh_databases_CM')

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

  ! sets number of top elements for surface movies & noise tomography
  NSPEC_TOP = NSPEC2D_TOP(IREGION_CRUST_MANTLE)

  ! reads databases file
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
            rmass_crust_mantle,rmass_ocean_load, &
            NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  ! check that the number of points in this slice is correct
  if(minval(ibool_crust_mantle(:,:,:,:)) /= 1 .or. &
    maxval(ibool_crust_mantle(:,:,:,:)) /= NGLOB_CRUST_MANTLE) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in crust and mantle')

  deallocate(dummy_i)

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
            rmass_outer_core,rmass_ocean_load, &
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
            rmass_inner_core,rmass_ocean_load, &
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
  integer :: njunk1,njunk2,njunk3
  integer :: ier

  ! crust and mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! Stacey put back
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening crust_mantle boundary.bin file')

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
        status='old',form='unformatted',action='read',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening outer_core boundary.bin file')

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
        status='old',form='unformatted',action='read',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening inner_core boundary.bin file')

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
          status='old',form='unformatted',action='read',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening boundary_disc.bin file')

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
  if (myrank == 0 .and. NCHUNKS_VAL == 6 .and. NPROCTOT_VAL < 1000 ) then
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
  ! (needed for stacey conditions)
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  end subroutine read_mesh_databases_addressing


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  ! local parameters
  real :: percentage_edge
  integer :: ier

  ! read MPI interfaces from file

  ! crust mantle
  call read_mesh_databases_MPI_CM()

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
  call read_mesh_databases_MPI_OC()

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
  call read_mesh_databases_MPI_IC()

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


  ! user output
  if(myrank == 0) then
    write(IMAIN,*) 'for overlapping of communications with calculations:'
    write(IMAIN,*)

    percentage_edge = 100. * nspec_outer_crust_mantle / real(NSPEC_CRUST_MANTLE)
    write(IMAIN,*) 'percentage of edge elements in crust/mantle ',percentage_edge,'%'
    write(IMAIN,*) 'percentage of volume elements in crust/mantle ',100. - percentage_edge,'%'
    write(IMAIN,*)

    percentage_edge = 100.* nspec_outer_outer_core / real(NSPEC_OUTER_CORE)
    write(IMAIN,*) 'percentage of edge elements in outer core ',percentage_edge,'%'
    write(IMAIN,*) 'percentage of volume elements in outer core ',100. - percentage_edge,'%'
    write(IMAIN,*)

    percentage_edge = 100. * nspec_outer_inner_core / real(NSPEC_INNER_CORE)
    write(IMAIN,*) 'percentage of edge elements in inner core ',percentage_edge,'%'
    write(IMAIN,*) 'percentage of volume elements in inner core ',100. - percentage_edge,'%'
    write(IMAIN,*)
  endif
  ! synchronizes MPI processes
  call sync_all()

  end subroutine read_mesh_databases_MPI

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI_CM()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
  integer :: ier

  ! crust mantle region

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_crust_mantle
  allocate(my_neighbours_crust_mantle(num_interfaces_crust_mantle), &
          nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
          stat=ier)
  if( ier /= 0 ) &
    call exit_mpi(myrank,'error allocating array my_neighbours_crust_mantle etc.')

  if( num_interfaces_crust_mantle > 0 ) then
    read(IIN) max_nibool_interfaces_crust_mantle
    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_crust_mantle,num_interfaces_crust_mantle), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array ibool_interfaces_crust_mantle')

    read(IIN) my_neighbours_crust_mantle
    read(IIN) nibool_interfaces_crust_mantle
    read(IIN) ibool_interfaces_crust_mantle
  else
    ! dummy array
    max_nibool_interfaces_crust_mantle = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array dummy ibool_interfaces_crust_mantle')
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_crust_mantle,nspec_outer_crust_mantle
  read(IIN) num_phase_ispec_crust_mantle
  if( num_phase_ispec_crust_mantle < 0 ) &
    call exit_mpi(myrank,'error num_phase_ispec_crust_mantle is < zero')

  allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),&
          stat=ier)
  if( ier /= 0 ) &
    call exit_mpi(myrank,'error allocating array phase_ispec_inner_crust_mantle')

  if(num_phase_ispec_crust_mantle > 0 ) read(IIN) phase_ispec_inner_crust_mantle

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    ! colors
    read(IIN) num_colors_outer_crust_mantle,num_colors_inner_crust_mantle

    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
            stat=ier)
    if( ier /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_crust_mantle array')

    read(IIN) num_elem_colors_crust_mantle
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
            stat=ier)
    if( ier /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_crust_mantle array')
  endif

  close(IIN)

  end subroutine read_mesh_databases_MPI_CM


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI_OC()

  use specfem_par
  use specfem_par_outercore
  implicit none

  ! local parameters
  integer :: ier

  ! crust mantle region

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_outer_core
  allocate(my_neighbours_outer_core(num_interfaces_outer_core), &
          nibool_interfaces_outer_core(num_interfaces_outer_core), &
          stat=ier)
  if( ier /= 0 ) &
    call exit_mpi(myrank,'error allocating array my_neighbours_outer_core etc.')

  if( num_interfaces_outer_core > 0 ) then
    read(IIN) max_nibool_interfaces_outer_core
    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_outer_core,num_interfaces_outer_core), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array ibool_interfaces_outer_core')

    read(IIN) my_neighbours_outer_core
    read(IIN) nibool_interfaces_outer_core
    read(IIN) ibool_interfaces_outer_core
  else
    ! dummy array
    max_nibool_interfaces_outer_core = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array dummy ibool_interfaces_outer_core')
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_outer_core,nspec_outer_outer_core
  read(IIN) num_phase_ispec_outer_core
  if( num_phase_ispec_outer_core < 0 ) &
    call exit_mpi(myrank,'error num_phase_ispec_outer_core is < zero')

  allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),&
          stat=ier)
  if( ier /= 0 ) &
    call exit_mpi(myrank,'error allocating array phase_ispec_inner_outer_core')

  if(num_phase_ispec_outer_core > 0 ) read(IIN) phase_ispec_inner_outer_core

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    ! colors
    read(IIN) num_colors_outer_outer_core,num_colors_inner_outer_core

    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
            stat=ier)
    if( ier /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_outer_core array')

    read(IIN) num_elem_colors_outer_core
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
            stat=ier)
    if( ier /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_outer_core array')
  endif

  close(IIN)

  end subroutine read_mesh_databases_MPI_OC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI_IC()

  use specfem_par
  use specfem_par_innercore
  implicit none

  ! local parameters
  integer :: ier

  ! crust mantle region

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_inner_core
  allocate(my_neighbours_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &
          stat=ier)
  if( ier /= 0 ) &
    call exit_mpi(myrank,'error allocating array my_neighbours_inner_core etc.')

  if( num_interfaces_inner_core > 0 ) then
    read(IIN) max_nibool_interfaces_inner_core
    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_inner_core,num_interfaces_inner_core), &
            stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array ibool_interfaces_inner_core')

    read(IIN) my_neighbours_inner_core
    read(IIN) nibool_interfaces_inner_core
    read(IIN) ibool_interfaces_inner_core
  else
    ! dummy array
    max_nibool_interfaces_inner_core = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error allocating array dummy ibool_interfaces_inner_core')
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_inner_core,nspec_outer_inner_core
  read(IIN) num_phase_ispec_inner_core
  if( num_phase_ispec_inner_core < 0 ) &
    call exit_mpi(myrank,'error num_phase_ispec_inner_core is < zero')

  allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),&
          stat=ier)
  if( ier /= 0 ) &
    call exit_mpi(myrank,'error allocating array phase_ispec_inner_inner_core')

  if(num_phase_ispec_inner_core > 0 ) read(IIN) phase_ispec_inner_inner_core

  ! mesh coloring for GPUs
  if( USE_MESH_COLORING_GPU ) then
    ! colors
    read(IIN) num_colors_outer_inner_core,num_colors_inner_inner_core

    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
            stat=ier)
    if( ier /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_inner_core array')

    read(IIN) num_elem_colors_inner_core
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
            stat=ier)
    if( ier /= 0 ) &
      call exit_mpi(myrank,'error allocating num_elem_colors_inner_core array')
  endif

  close(IIN)

  end subroutine read_mesh_databases_MPI_IC


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
      call open_file_abs_r(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    endif
  endif

  end subroutine read_mesh_databases_stacey

