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

  subroutine read_mesh_databases()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  ! timing
  double precision :: tCPU
  double precision, external :: wtime

#ifdef USE_SERIAL_CASCADE_FOR_IOs
  logical,parameter :: SYNC_READING = .false.
#else
  logical,parameter :: SYNC_READING = .true.
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'mesh databases:'
    call flush_IMAIN()
  endif

  ! get MPI starting time
  time_start = wtime()

#ifdef USE_SERIAL_CASCADE_FOR_IOs
  ! serial i/o
  you_can_start_doing_IOs = .false.
  if (myrank > 0) call recv_singlel(you_can_start_doing_IOs,myrank-1,itag)
#endif

  ! start reading the databases
  ! read arrays created by the mesher

  ! reads "solver_data.bin" files for crust and mantle
  if (SYNC_READING ) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in crust/mantle databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_CM()

  ! reads "solver_data.bin" files for outer core
  if (SYNC_READING ) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in outer core databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_OC()

  ! reads "solver_data.bin" files for inner core
  if (SYNC_READING ) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in inner core databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_IC()

  ! reads "boundary.bin" files to couple mantle with outer core and inner core boundaries
  if (SYNC_READING ) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in coupling surface databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_coupling()

  ! reads "addressing.txt" 2-D addressing (needed for Stacey boundaries and
  ! regular grid kernels)
  if (SYNC_READING ) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in addressing...'
    call flush_IMAIN()
  endif
  allocate(addressing(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1))
  call read_mesh_databases_addressing()
  if (.not.SAVE_REGULAR_KL) then
    deallocate(addressing)
  endif

  ! sets up MPI interfaces, inner/outer elements and mesh coloring
  if (SYNC_READING ) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in MPI databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_MPI()

  ! absorbing boundaries
  if (ABSORBING_CONDITIONS) then
    ! reads "stacey.bin" files
    if (SYNC_READING ) call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  reading in Stacey databases...'
      call flush_IMAIN()
    endif
    call read_mesh_databases_stacey()
  endif

  ! kernels on regular grids
  if (SAVE_REGULAR_KL) then
    if (SYNC_READING ) call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  reading in regular kernel databases...'
      call flush_IMAIN()
    endif
    call read_mesh_databases_regular_kl()
    deallocate(addressing)
  endif

#ifdef USE_SERIAL_CASCADE_FOR_IOs
  ! serial i/o
  you_can_start_doing_IOs = .true.
  if (myrank < NPROC_XI_VAL*NPROC_ETA_VAL-1) call send_singlel(you_can_start_doing_IOs,myrank+1,itag)
#endif

  ! user output
  call synchronize_all()
  if (myrank == 0) then
    ! elapsed time since beginning of mesh generation
    tCPU = wtime() - time_start
    write(IMAIN,*)
    write(IMAIN,*) 'Elapsed time for reading mesh in seconds = ',sngl(tCPU)
    write(IMAIN,*)
    call flush_IMAIN()
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
  integer, dimension(:),allocatable :: dummy_idoubling
  integer :: ier

  ! crust and mantle

  if (ANISOTROPIC_3D_MANTLE_VAL) then
    READ_KAPPA_MU = .false.
    READ_TISO = .false.
    nspec_iso = NSPECMAX_ISO_MANTLE ! 1
    nspec_tiso = NSPECMAX_TISO_MANTLE ! 1
    nspec_ani = NSPEC_CRUST_MANTLE
  else
    READ_KAPPA_MU = .true.
    nspec_iso = NSPEC_CRUST_MANTLE
    if (TRANSVERSE_ISOTROPY_VAL) then
      nspec_tiso = NSPECMAX_TISO_MANTLE
    else
      nspec_tiso = 1
    endif
    nspec_ani = NSPECMAX_ANISO_MANTLE ! 1
    READ_TISO = .true.
  endif

  ! sets number of top elements for surface movies & noise tomography
  NSPEC_TOP = NSPEC2D_TOP(IREGION_CRUST_MANTLE)

  ! allocates dummy array
  allocate(dummy_idoubling(NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy idoubling in crust_mantle')

  ! allocates mass matrices in this slice (will be fully assembled in the solver)
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete

  ! allocates mass matrices
  allocate(rmassx_crust_mantle(NGLOB_XY_CM), &
           rmassy_crust_mantle(NGLOB_XY_CM),stat=ier)
  if (ier /= 0) stop 'Error allocating rmassx, rmassy in crust_mantle'

  ! b_rmassx and b_rmassy will be different to rmassx and rmassy
  ! needs new arrays
  allocate(b_rmassx_crust_mantle(NGLOB_XY_CM), &
           b_rmassy_crust_mantle(NGLOB_XY_CM),stat=ier)
  if (ier /= 0) stop 'Error allocating b_rmassx, b_rmassy in crust_mantle'

  ! reads databases file
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_CRUST_MANTLE,myrank, &
                                    NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,NGLOB_XY_CM, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    rho_vp_crust_mantle,rho_vs_crust_mantle, &
                                    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                                    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                    rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                    kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                    c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                                    c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
                                    c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
                                    c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
                                    c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
                                    c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                    c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                    ibool_crust_mantle,dummy_idoubling,ispec_is_tiso_crust_mantle, &
                                    rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                    NGLOB_CRUST_MANTLE_OCEANS,rmass_ocean_load, &
                                    READ_KAPPA_MU,READ_TISO, &
                                    b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    else
      call read_arrays_solver(IREGION_CRUST_MANTLE,myrank, &
                              NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,NGLOB_XY_CM, &
                              nspec_iso,nspec_tiso,nspec_ani, &
                              rho_vp_crust_mantle,rho_vs_crust_mantle, &
                              xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                              xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                              etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                              gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                              rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                              kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                              c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                              c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
                              c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
                              c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
                              c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
                              c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                              c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                              ibool_crust_mantle,dummy_idoubling,ispec_is_tiso_crust_mantle, &
                              rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                              NGLOB_CRUST_MANTLE_OCEANS,rmass_ocean_load, &
                              READ_KAPPA_MU,READ_TISO, &
                              b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    endif
  endif
  call bcast_mesh_databases_CM()

  ! check that the number of points in this slice is correct
  if (minval(ibool_crust_mantle(:,:,:,:)) /= 1) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmin is not equal to 1 in crust and mantle')
  if (maxval(ibool_crust_mantle(:,:,:,:)) /= NGLOB_CRUST_MANTLE) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in crust and mantle')

  deallocate(dummy_idoubling)

  ! mass matrix corrections
  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION)) then
    ! mass matrices differ for rmassx,rmassy
    ! continue
  else
    ! uses single mass matrix without correction
    ! frees pointer memory
    deallocate(rmassx_crust_mantle,rmassy_crust_mantle)
    ! re-associates with corresponding rmassz
    rmassx_crust_mantle => rmassz_crust_mantle(:)
    rmassy_crust_mantle => rmassz_crust_mantle(:)
  endif

  ! kernel simulations
  if (SIMULATION_TYPE == 3) then
    ! associates mass matrix used for backward/reconstructed wavefields
    b_rmassz_crust_mantle => rmassz_crust_mantle
    ! checks if we can take rmassx and rmassy (only differs for rotation correction)
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
      ! mass matrices differ for b_rmassx,b_rmassy
      ! continue
    else
      ! frees pointer memory
      deallocate(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
      ! re-associates with corresponding rmassx,rmassy
      b_rmassx_crust_mantle => rmassx_crust_mantle(:)
      b_rmassy_crust_mantle => rmassy_crust_mantle(:)
    endif
  else
    ! b_rmassx,b_rmassy not used anymore
    deallocate(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    nullify(b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    nullify(b_rmassz_crust_mantle)
  endif

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
  integer :: nspec_iso,nspec_tiso,nspec_ani,NGLOB_XY_dummy
  logical :: READ_KAPPA_MU,READ_TISO
  integer :: ier

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dummy_rmass

  logical, dimension(:), allocatable :: dummy_ispec_is_tiso
  integer, dimension(:), allocatable :: dummy_idoubling_outer_core

  ! outer core (no anisotropy nor S velocity)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the outer core
  READ_KAPPA_MU = .false.
  READ_TISO = .false.
  nspec_iso = NSPEC_OUTER_CORE
  nspec_tiso = 1
  nspec_ani = 1

  ! dummy allocation
  NGLOB_XY_dummy = 1

  allocate(dummy_rmass(NGLOB_XY_dummy), &
          dummy_ispec_is_tiso(NSPEC_OUTER_CORE), &
          dummy_idoubling_outer_core(NSPEC_OUTER_CORE), &
          stat=ier)
  if (ier /= 0) stop 'Error allocating dummy rmass and dummy ispec/idoubling in outer core'

  ! reads in mesh arrays
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_OUTER_CORE,myrank, &
                                    NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,NGLOB_XY_dummy, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    vp_outer_core,dummy_array, &
                                    xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                                    xix_outer_core,xiy_outer_core,xiz_outer_core, &
                                    etax_outer_core,etay_outer_core,etaz_outer_core, &
                                    gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                                    rhostore_outer_core,kappavstore_outer_core,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    ibool_outer_core,dummy_idoubling_outer_core,dummy_ispec_is_tiso, &
                                    dummy_rmass,dummy_rmass,rmass_outer_core, &
                                    1,dummy_array, &
                                    READ_KAPPA_MU,READ_TISO, &
                                    dummy_rmass,dummy_rmass)
    else
      call read_arrays_solver(IREGION_OUTER_CORE,myrank, &
                              NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,NGLOB_XY_dummy, &
                              nspec_iso,nspec_tiso,nspec_ani, &
                              vp_outer_core,dummy_array, &
                              xstore_outer_core,ystore_outer_core,zstore_outer_core, &
                              xix_outer_core,xiy_outer_core,xiz_outer_core, &
                              etax_outer_core,etay_outer_core,etaz_outer_core, &
                              gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
                              rhostore_outer_core,kappavstore_outer_core,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              ibool_outer_core,dummy_idoubling_outer_core,dummy_ispec_is_tiso, &
                              dummy_rmass,dummy_rmass,rmass_outer_core, &
                              1, dummy_array, &
                              READ_KAPPA_MU,READ_TISO, &
                              dummy_rmass,dummy_rmass)
    endif
  endif
  call bcast_mesh_databases_OC()

  deallocate(dummy_idoubling_outer_core,dummy_ispec_is_tiso,dummy_rmass)

  ! check that the number of points in this slice is correct
  ! check that the number of points in this slice is correct
  if (minval(ibool_outer_core(:,:,:,:)) /= 1) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmin is not equal to 1 in outer core')
  if (maxval(ibool_outer_core(:,:,:,:)) /= NGLOB_OUTER_CORE) then
    call exit_MPI(myrank, 'incorrect global numbering: iboolmax does not equal nglob in outer core')
  endif

  ! kernel simulations
  if (SIMULATION_TYPE == 3) then
    ! associates mass matrix used for backward/reconstructed wavefields
    b_rmass_outer_core => rmass_outer_core
  else
    nullify(b_rmass_outer_core)
  endif

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
  integer :: ier

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  logical, dimension(:),allocatable:: dummy_ispec_is_tiso

  ! inner core (no anisotropy)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the inner core
  READ_KAPPA_MU = .true. ! (muvstore needed for attenuation)
  READ_TISO = .false.
  nspec_iso = NSPEC_INNER_CORE
  nspec_tiso = 1
  if (ANISOTROPIC_INNER_CORE_VAL) then
    nspec_ani = NSPEC_INNER_CORE
  else
    nspec_ani = 1
  endif

  allocate(dummy_ispec_is_tiso(NSPEC_INNER_CORE), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating dummy ispec in inner core'

  ! allocates mass matrices in this slice (will be fully assembled in the solver)
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  allocate(rmassx_inner_core(NGLOB_XY_IC), &
           rmassy_inner_core(NGLOB_XY_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating rmassx, rmassy in inner_core'

  ! b_rmassx and b_rmassy maybe different to rmassx,rmassy
  allocate(b_rmassx_inner_core(NGLOB_XY_IC), &
           b_rmassy_inner_core(NGLOB_XY_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating b_rmassx, b_rmassy in inner_core'

  ! reads in arrays
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_INNER_CORE,myrank, &
                                    NSPEC_INNER_CORE,NGLOB_INNER_CORE,NGLOB_XY_IC, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    dummy_array,dummy_array, &
                                    xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                                    xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                    etax_inner_core,etay_inner_core,etaz_inner_core, &
                                    gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                    rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                    dummy_array,dummy_array,dummy_array, &
                                    c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,c33store_inner_core, &
                                    dummy_array,dummy_array,dummy_array, &
                                    c44store_inner_core,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    ibool_inner_core,idoubling_inner_core,dummy_ispec_is_tiso, &
                                    rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                                    1,dummy_array, &
                                    READ_KAPPA_MU,READ_TISO, &
                                    b_rmassx_inner_core,b_rmassy_inner_core)
    else
      call read_arrays_solver(IREGION_INNER_CORE,myrank, &
                              NSPEC_INNER_CORE,NGLOB_INNER_CORE,NGLOB_XY_IC, &
                              nspec_iso,nspec_tiso,nspec_ani, &
                              dummy_array,dummy_array, &
                              xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                              xix_inner_core,xiy_inner_core,xiz_inner_core, &
                              etax_inner_core,etay_inner_core,etaz_inner_core, &
                              gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                              rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                              dummy_array,dummy_array,dummy_array, &
                              c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              dummy_array,dummy_array,c33store_inner_core, &
                              dummy_array,dummy_array,dummy_array, &
                              c44store_inner_core,dummy_array,dummy_array, &
                              dummy_array,dummy_array,dummy_array, &
                              ibool_inner_core,idoubling_inner_core,dummy_ispec_is_tiso, &
                              rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                              1,dummy_array, &
                              READ_KAPPA_MU,READ_TISO, &
                              b_rmassx_inner_core,b_rmassy_inner_core)
    endif
  endif
  call bcast_mesh_databases_IC()

  deallocate(dummy_ispec_is_tiso)

  ! check that the number of points in this slice is correct
  if (minval(ibool_inner_core(:,:,:,:)) /= 1 .or. maxval(ibool_inner_core(:,:,:,:)) /= NGLOB_INNER_CORE) &
    call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in inner core')

  ! mass matrix corrections
  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    ! uses corrected mass matrices
    ! continue
  else
    ! uses single mass matrix without correction
    ! frees pointer memory
    deallocate(rmassx_inner_core,rmassy_inner_core)
    ! re-associates with corresponding rmassz
    rmassx_inner_core => rmassz_inner_core(:)
    rmassy_inner_core => rmassz_inner_core(:)
  endif

  ! kernel simulations
  if (SIMULATION_TYPE == 3) then
    ! associates mass matrix used for backward/reconstructed wavefields
    b_rmassz_inner_core => rmassz_inner_core
    ! checks if we can take rmassx and rmassy (only differs for rotation correction)
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
      ! uses corrected mass matrices
      ! continue
    else
      ! frees pointer memory
      deallocate(b_rmassx_inner_core,b_rmassy_inner_core)
      ! re-associates with corresponding rmassx,rmassy
      b_rmassx_inner_core => rmassx_inner_core
      b_rmassy_inner_core => rmassy_inner_core
    endif
  else
    ! b_rmassx,b_rmassy not used anymore
    deallocate(b_rmassx_inner_core,b_rmassy_inner_core)
    nullify(b_rmassx_inner_core,b_rmassy_inner_core)
    ! use dummy pointers used for passing as function arguments
    ! associates mass matrix used for backward/reconstructed wavefields
    !b_rmassz_inner_core => rmassz_inner_core
    !b_rmassx_inner_core => rmassz_inner_core
    !b_rmassy_inner_core => rmassz_inner_core
    nullify(b_rmassz_inner_core)
  endif

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

  ! local parameters
  integer :: njunk1,njunk2,njunk3
  integer :: ier

  ! reads in arrays
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_mesh_databases_coupling_adios()
    else
      ! crust and mantle
      ! create name of database
      call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

      ! Stacey put back
      open(unit=IIN,file=prname(1:len_trim(prname))//'boundary.bin', &
            status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening crust_mantle boundary.bin file')

      read(IIN) nspec2D_xmin_crust_mantle
      read(IIN) nspec2D_xmax_crust_mantle
      read(IIN) nspec2D_ymin_crust_mantle
      read(IIN) nspec2D_ymax_crust_mantle
      read(IIN) njunk1
      read(IIN) njunk2

    ! boundary parameters
      read(IIN) ibelm_xmin_crust_mantle
      read(IIN) ibelm_xmax_crust_mantle
      read(IIN) ibelm_ymin_crust_mantle
      read(IIN) ibelm_ymax_crust_mantle
      read(IIN) ibelm_bottom_crust_mantle
      read(IIN) ibelm_top_crust_mantle

      read(IIN) normal_xmin_crust_mantle
      read(IIN) normal_xmax_crust_mantle
      read(IIN) normal_ymin_crust_mantle
      read(IIN) normal_ymax_crust_mantle
      read(IIN) normal_bottom_crust_mantle
      read(IIN) normal_top_crust_mantle

      read(IIN) jacobian2D_xmin_crust_mantle
      read(IIN) jacobian2D_xmax_crust_mantle
      read(IIN) jacobian2D_ymin_crust_mantle
      read(IIN) jacobian2D_ymax_crust_mantle
      read(IIN) jacobian2D_bottom_crust_mantle
      read(IIN) jacobian2D_top_crust_mantle
      close(IIN)

      ! read parameters to couple fluid and solid regions
      !
      ! outer core

      ! create name of database
      call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

      ! boundary parameters

      ! Stacey put back
      open(unit=IIN,file=prname(1:len_trim(prname))//'boundary.bin', &
            status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening outer_core boundary.bin file')

      read(IIN) nspec2D_xmin_outer_core
      read(IIN) nspec2D_xmax_outer_core
      read(IIN) nspec2D_ymin_outer_core
      read(IIN) nspec2D_ymax_outer_core
      read(IIN) njunk1
      read(IIN) njunk2

      nspec2D_zmin_outer_core = NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

      read(IIN) ibelm_xmin_outer_core
      read(IIN) ibelm_xmax_outer_core
      read(IIN) ibelm_ymin_outer_core
      read(IIN) ibelm_ymax_outer_core
      read(IIN) ibelm_bottom_outer_core
      read(IIN) ibelm_top_outer_core

      read(IIN) normal_xmin_outer_core
      read(IIN) normal_xmax_outer_core
      read(IIN) normal_ymin_outer_core
      read(IIN) normal_ymax_outer_core
      read(IIN) normal_bottom_outer_core
      read(IIN) normal_top_outer_core

      read(IIN) jacobian2D_xmin_outer_core
      read(IIN) jacobian2D_xmax_outer_core
      read(IIN) jacobian2D_ymin_outer_core
      read(IIN) jacobian2D_ymax_outer_core
      read(IIN) jacobian2D_bottom_outer_core
      read(IIN) jacobian2D_top_outer_core
      close(IIN)

      !
      ! inner core
      !

      ! create name of database
      call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_PATH)

      ! read info for vertical edges for central cube matching in inner core
      open(unit=IIN,file=prname(1:len_trim(prname))//'boundary.bin', &
            status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening inner_core boundary.bin file')

      read(IIN) nspec2D_xmin_inner_core
      read(IIN) nspec2D_xmax_inner_core
      read(IIN) nspec2D_ymin_inner_core
      read(IIN) nspec2D_ymax_inner_core
      read(IIN) njunk1
      read(IIN) njunk2

      ! boundary parameters
      read(IIN) ibelm_xmin_inner_core
      read(IIN) ibelm_xmax_inner_core
      read(IIN) ibelm_ymin_inner_core
      read(IIN) ibelm_ymax_inner_core
      read(IIN) ibelm_bottom_inner_core
      read(IIN) ibelm_top_inner_core
      close(IIN)

      ! -- Boundary Mesh for crust and mantle ---
      if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
        call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

        open(unit=IIN,file=prname(1:len_trim(prname))//'boundary_disc.bin', &
              status='old',form='unformatted',action='read',iostat=ier)
        if (ier /= 0 ) call exit_mpi(myrank,'Error opening boundary_disc.bin file')

        read(IIN) njunk1,njunk2,njunk3
        if (njunk1 /= NSPEC2D_MOHO .and. njunk2 /= NSPEC2D_400 .and. njunk3 /= NSPEC2D_670) &
                   call exit_mpi(myrank, 'Error reading ibelm_disc.bin file')
        read(IIN) ibelm_moho_top
        read(IIN) ibelm_moho_bot
        read(IIN) ibelm_400_top
        read(IIN) ibelm_400_bot
        read(IIN) ibelm_670_top
        read(IIN) ibelm_670_bot
        read(IIN) normal_moho
        read(IIN) normal_400
        read(IIN) normal_670
        close(IIN)
      endif

    endif ! ADIOS
  endif
  call bcast_mesh_databases_coupling()

  ! checks dimensions
  ! crust mantle
  if (nspec2d_xmin_crust_mantle < 0 .or. nspec2d_xmin_crust_mantle > NSPEC2DMAX_XMIN_XMAX_CM .or. &
      nspec2d_xmax_crust_mantle < 0 .or. nspec2d_xmax_crust_mantle > NSPEC2DMAX_XMIN_XMAX_CM .or. &
      nspec2d_ymin_crust_mantle < 0 .or. nspec2d_ymin_crust_mantle > NSPEC2DMAX_YMIN_YMAX_CM .or. &
      nspec2d_ymax_crust_mantle < 0 .or. nspec2d_ymax_crust_mantle > NSPEC2DMAX_YMIN_YMAX_CM ) &
      call exit_mpi(myrank,'Error reading crust/mantle boundary')
  ! outer core
  if (nspec2D_xmin_outer_core < 0 .or. nspec2d_xmin_outer_core > NSPEC2DMAX_XMIN_XMAX_OC .or. &
      nspec2D_xmax_outer_core < 0 .or. nspec2d_xmax_outer_core > NSPEC2DMAX_XMIN_XMAX_OC .or. &
      nspec2D_ymin_outer_core < 0 .or. nspec2d_ymin_outer_core > NSPEC2DMAX_YMIN_YMAX_OC .or. &
      nspec2D_ymax_outer_core < 0 .or. nspec2d_ymax_outer_core > NSPEC2DMAX_YMIN_YMAX_OC ) &
    call exit_mpi(myrank, 'Error reading outer core boundary')
  ! inner core
  if (nspec2D_xmin_inner_core < 0 .or. nspec2d_xmin_inner_core > NSPEC2DMAX_XMIN_XMAX_IC .or. &
      nspec2D_xmax_inner_core < 0 .or. nspec2d_xmax_inner_core > NSPEC2DMAX_XMIN_XMAX_IC .or. &
      nspec2D_ymin_inner_core < 0 .or. nspec2d_ymin_inner_core > NSPEC2DMAX_YMIN_YMAX_IC .or. &
      nspec2D_ymax_inner_core < 0 .or. nspec2d_ymax_inner_core > NSPEC2DMAX_YMIN_YMAX_IC ) &
    call exit_mpi(myrank, 'Error reading inner core boundary')

  ! initializes
  nspec2D_zmin_outer_core = NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

  ! Boundary Mesh for crust and mantle
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
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

  use constants

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer, dimension(0:NPROCTOT_VAL-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer :: ier,iproc,iproc_read,iproc_xi,iproc_eta

  if (I_should_read_the_database) then
    ! open file with global slice number addressing
    if (myrank == 0) then
      open(unit=IIN,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read',iostat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error opening addressing.txt')

      do iproc = 0,NPROCTOT_VAL-1
        read(IIN,*) iproc_read,ichunk,iproc_xi,iproc_eta

        if (iproc_read /= iproc) call exit_MPI(myrank,'incorrect slice number read')

        addressing(ichunk,iproc_xi,iproc_eta) = iproc
        ichunk_slice(iproc) = ichunk
        iproc_xi_slice(iproc) = iproc_xi
        iproc_eta_slice(iproc) = iproc_eta
      enddo
      close(IIN)
    endif

    ! broadcast the information read on the master to the nodes
    call bcast_all_i(addressing,NCHUNKS_VAL*NPROC_XI_VAL*NPROC_ETA_VAL)
    call bcast_all_i(ichunk_slice,NPROCTOT_VAL)
    call bcast_all_i(iproc_xi_slice,NPROCTOT_VAL)
    call bcast_all_i(iproc_eta_slice,NPROCTOT_VAL)

    ! output a topology map of slices - fix 20x by nproc
    if (myrank == 0) then
      if (NCHUNKS_VAL == 6 .and. NPROCTOT_VAL < 1000) then
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
    endif

    ! determine chunk number and local slice coordinates using addressing
    ! (needed for Stacey conditions)
    ichunk = ichunk_slice(myrank)
  endif ! I_should_read_the_database

  call bcast_all_i_for_database(ichunk, 1)
  call bcast_all_i_for_database(addressing(1,0,0), size(addressing))

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
  if (I_should_read_the_database) then
    if (ADIOS_FOR_MPI_ARRAYS) then
      call read_mesh_databases_MPI_CM_adios()
    else
      call read_mesh_databases_MPI_CM()
    endif
    !call synchronize_all()
    !print *, "read_mesh_databases_MPI_CM done"
    !call flush()
  endif
  call bcast_mesh_databases_MPI_CM()

  allocate(buffer_send_vector_crust_mantle(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
           buffer_recv_vector_crust_mantle(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
           request_send_vector_cm(num_interfaces_crust_mantle), &
           request_recv_vector_cm(num_interfaces_crust_mantle), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_vector_crust_mantle etc.')

  if (SIMULATION_TYPE == 3) then
    allocate(b_buffer_send_vector_cm(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
             b_buffer_recv_vector_cm(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
             b_request_send_vector_cm(num_interfaces_crust_mantle), &
             b_request_recv_vector_cm(num_interfaces_crust_mantle), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_vector_cm etc.')
  endif

  ! outer core
  if (I_should_read_the_database) then
    if (ADIOS_FOR_MPI_ARRAYS) then
      call read_mesh_databases_MPI_OC_adios()
    else
      call read_mesh_databases_MPI_OC()
    endif
  endif
  call bcast_mesh_databases_MPI_OC()

  allocate(buffer_send_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
           buffer_recv_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
           request_send_scalar_oc(num_interfaces_outer_core), &
           request_recv_scalar_oc(num_interfaces_outer_core), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_vector_outer_core etc.')

  if (SIMULATION_TYPE == 3) then
    allocate(b_buffer_send_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
             b_buffer_recv_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
             b_request_send_scalar_oc(num_interfaces_outer_core), &
             b_request_recv_scalar_oc(num_interfaces_outer_core), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_vector_outer_core etc.')
  endif

  ! inner core
  if (I_should_read_the_database) then
    if (ADIOS_FOR_MPI_ARRAYS) then
      call read_mesh_databases_MPI_IC_adios()
    else
      call read_mesh_databases_MPI_IC()
    endif
  endif
  call bcast_mesh_databases_MPI_IC()

  allocate(buffer_send_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core), &
           buffer_recv_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core), &
           request_send_vector_ic(num_interfaces_inner_core), &
           request_recv_vector_ic(num_interfaces_inner_core), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_vector_inner_core etc.')

  if (SIMULATION_TYPE == 3) then
    allocate(b_buffer_send_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core), &
             b_buffer_recv_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core), &
             b_request_send_vector_ic(num_interfaces_inner_core), &
             b_request_recv_vector_ic(num_interfaces_inner_core), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_vector_inner_core etc.')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  for overlapping of communications with calculations:'
    write(IMAIN,*)

    percentage_edge = 100. * nspec_outer_crust_mantle / real(NSPEC_CRUST_MANTLE)
    write(IMAIN,*) '  percentage of edge elements in crust/mantle ',percentage_edge,'%'
    write(IMAIN,*) '  percentage of volume elements in crust/mantle ',100. - percentage_edge,'%'
    write(IMAIN,*)

    percentage_edge = 100.* nspec_outer_outer_core / real(NSPEC_OUTER_CORE)
    write(IMAIN,*) '  percentage of edge elements in outer core ',percentage_edge,'%'
    write(IMAIN,*) '  percentage of volume elements in outer core ',100. - percentage_edge,'%'
    write(IMAIN,*)

    percentage_edge = 100. * nspec_outer_inner_core / real(NSPEC_INNER_CORE)
    write(IMAIN,*) '  percentage of edge elements in inner core ',percentage_edge,'%'
    write(IMAIN,*) '  percentage of volume elements in inner core ',100. - percentage_edge,'%'
    write(IMAIN,*)
    call flush_IMAIN()
  endif
  ! synchronizes MPI processes
  call synchronize_all()

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
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_crust_mantle
  allocate(my_neighbours_crust_mantle(num_interfaces_crust_mantle), &
          nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
          stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array my_neighbours_crust_mantle etc.')

  if (num_interfaces_crust_mantle > 0) then
    read(IIN) max_nibool_interfaces_cm
    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')

    read(IIN) my_neighbours_crust_mantle
    read(IIN) nibool_interfaces_crust_mantle
    read(IIN) ibool_interfaces_crust_mantle
  else
    ! dummy array
    max_nibool_interfaces_cm = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_crust_mantle')
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_crust_mantle,nspec_outer_crust_mantle
  read(IIN) num_phase_ispec_crust_mantle
  if (num_phase_ispec_crust_mantle < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_crust_mantle is < zero')

  allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),&
          stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')

  if (num_phase_ispec_crust_mantle > 0 ) read(IIN) phase_ispec_inner_crust_mantle

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_crust_mantle,num_colors_inner_crust_mantle

    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')

    read(IIN) num_elem_colors_crust_mantle
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
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
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_outer_core
  allocate(my_neighbours_outer_core(num_interfaces_outer_core), &
          nibool_interfaces_outer_core(num_interfaces_outer_core), &
          stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array my_neighbours_outer_core etc.')

  if (num_interfaces_outer_core > 0) then
    read(IIN) max_nibool_interfaces_oc
    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_outer_core')

    read(IIN) my_neighbours_outer_core
    read(IIN) nibool_interfaces_outer_core
    read(IIN) ibool_interfaces_outer_core
  else
    ! dummy array
    max_nibool_interfaces_oc = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_outer_core')
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_outer_core,nspec_outer_outer_core
  read(IIN) num_phase_ispec_outer_core
  if (num_phase_ispec_outer_core < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_outer_core is < zero')

  allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),&
          stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')

  if (num_phase_ispec_outer_core > 0 ) read(IIN) phase_ispec_inner_outer_core

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_outer_core,num_colors_inner_outer_core

    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')

    read(IIN) num_elem_colors_outer_core
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
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
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_inner_core
  allocate(my_neighbours_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &
          stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array my_neighbours_inner_core etc.')

  if (num_interfaces_inner_core > 0) then
    read(IIN) max_nibool_interfaces_ic
    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic,num_interfaces_inner_core), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_inner_core')

    read(IIN) my_neighbours_inner_core
    read(IIN) nibool_interfaces_inner_core
    read(IIN) ibool_interfaces_inner_core
  else
    ! dummy array
    max_nibool_interfaces_ic = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_inner_core')
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_inner_core,nspec_outer_inner_core
  read(IIN) num_phase_ispec_inner_core
  if (num_phase_ispec_inner_core < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_inner_core is < zero')

  allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),&
          stat=ier)
  if (ier /= 0 ) &
    call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')

  if (num_phase_ispec_inner_core > 0 ) read(IIN) phase_ispec_inner_inner_core

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_inner_core,num_colors_inner_inner_core

    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')

    read(IIN) num_elem_colors_inner_core
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
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
  integer :: ier

  ! reads in arrays
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_mesh_databases_stacey_adios()
    else
      ! crust and mantle

      ! create name of database
      call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

      ! read arrays for Stacey conditions
      open(unit=IIN,file=prname(1:len_trim(prname))//'stacey.bin', &
            status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey.bin file for crust mantle')

      read(IIN) nimin_crust_mantle
      read(IIN) nimax_crust_mantle
      read(IIN) njmin_crust_mantle
      read(IIN) njmax_crust_mantle
      read(IIN) nkmin_xi_crust_mantle
      read(IIN) nkmin_eta_crust_mantle
      close(IIN)

      ! outer core

      ! create name of database
      call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

      ! read arrays for Stacey conditions
      open(unit=IIN,file=prname(1:len_trim(prname))//'stacey.bin', &
            status='old',form='unformatted',action='read',iostat=ier)
      if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey.bin file for outer core')

      read(IIN) nimin_outer_core
      read(IIN) nimax_outer_core
      read(IIN) njmin_outer_core
      read(IIN) njmax_outer_core
      read(IIN) nkmin_xi_outer_core
      read(IIN) nkmin_eta_outer_core
      close(IIN)
    endif ! ADIOS
  endif ! I_should_read_the_database
  call bcast_mesh_databases_stacey()

  end subroutine read_mesh_databases_stacey

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_regular_kl()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer, dimension(:), allocatable :: slice_number
  integer :: i,isp,ier
  ! grid parameters
  type kl_reg_grid_variables
    sequence
    real dlat
    real dlon
    integer nlayer
    real rlayer(NM_KL_REG_LAYER)
    integer ndoubling(NM_KL_REG_LAYER)
    integer nlat(NM_KL_REG_LAYER)
    integer nlon(NM_KL_REG_LAYER)
    integer npts_total
    integer npts_before_layer(NM_KL_REG_LAYER+1)
  end type kl_reg_grid_variables
  type (kl_reg_grid_variables) KL_REG_GRID

  call read_kl_regular_grid(KL_REG_GRID)

  if (myrank == 0) then
    ! master process
    allocate(slice_number(KL_REG_GRID%npts_total),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating slice_number array')

    ! print *, 'slice npts =', KL_REG_GRID%npts_total
    call find_regular_grid_slice_number(slice_number, KL_REG_GRID)

    do i = NPROCTOT_VAL-1,0,-1
      npoints_slice = 0
      do isp = 1,KL_REG_GRID%npts_total
        if (slice_number(isp) == i) then
          npoints_slice = npoints_slice + 1
          if (npoints_slice > NM_KL_REG_PTS) stop 'Exceeding NM_KL_REG_PTS limit'
          points_slice(npoints_slice) = isp
        endif
      enddo

      if (i /= 0) then
        call send_singlei(npoints_slice,i,i)
        if (npoints_slice > 0) then
          call send_i(points_slice,npoints_slice,i,2*i)
        endif
      endif
    enddo

    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/kl_grid_slice.txt',status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file kl_grid_slice.txt for writing')
    write(IOUT,*) slice_number
    close(IOUT)

    deallocate(slice_number)
  else
    ! slave processes
    call recv_singlei(npoints_slice,0,myrank)
    if (npoints_slice > 0) then
      call recv_i(points_slice,npoints_slice,0,2*myrank)
    endif
  endif

  ! this is the core part that takes up most of the computation time,
  ! and presumably the more processors involved the faster.
  if (npoints_slice > 0) then
    call locate_regular_points(npoints_slice, points_slice, KL_REG_GRID, &
                               NSPEC_CRUST_MANTLE, &
                               xstore_crust_mantle, ystore_crust_mantle, zstore_crust_mantle, &
                               ibool_crust_mantle, &
                               xigll, yigll, zigll, &
                               ispec_reg, hxir_reg, hetar_reg, hgammar_reg)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) ' '
    write(IMAIN,*) 'Finished locating kernel output regular grid'
    write(IMAIN,*) ' '
    call flush_IMAIN()
  endif

  end subroutine read_mesh_databases_regular_kl

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_CM()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  !call bcast_all_i_for_database(NSPEC_CRUST_MANTLE, 1)
  !call bcast_all_i_for_database(NGLOB_CRUST_MANTLE, 1)
  !call bcast_all_i_for_database(NGLOB_XY_CM, 1)
  call bcast_all_cr_for_database(rho_vp_crust_mantle(1,1,1,1), size(rho_vp_crust_mantle))
  call bcast_all_cr_for_database(rho_vs_crust_mantle(1,1,1,1), size(rho_vs_crust_mantle))

  call bcast_all_cr_for_database(xstore_crust_mantle(1), size(xstore_crust_mantle))
  call bcast_all_cr_for_database(ystore_crust_mantle(1), size(ystore_crust_mantle))
  call bcast_all_cr_for_database(zstore_crust_mantle(1), size(zstore_crust_mantle))
  call bcast_all_cr_for_database(xix_crust_mantle(1,1,1,1), size(xix_crust_mantle))
  call bcast_all_cr_for_database(xiy_crust_mantle(1,1,1,1), size(xiy_crust_mantle))
  call bcast_all_cr_for_database(xiz_crust_mantle(1,1,1,1), size(xiz_crust_mantle))
  call bcast_all_cr_for_database(etax_crust_mantle(1,1,1,1), size(etax_crust_mantle))
  call bcast_all_cr_for_database(etay_crust_mantle(1,1,1,1), size(etay_crust_mantle))
  call bcast_all_cr_for_database(etaz_crust_mantle(1,1,1,1), size(etaz_crust_mantle))
  call bcast_all_cr_for_database(gammax_crust_mantle(1,1,1,1), size(gammax_crust_mantle))
  call bcast_all_cr_for_database(gammay_crust_mantle(1,1,1,1), size(gammay_crust_mantle))
  call bcast_all_cr_for_database(gammaz_crust_mantle(1,1,1,1), size(gammaz_crust_mantle))
  call bcast_all_cr_for_database(rhostore_crust_mantle(1,1,1,1), size(rhostore_crust_mantle))
  call bcast_all_cr_for_database(kappavstore_crust_mantle(1,1,1,1), size(kappavstore_crust_mantle))
  call bcast_all_cr_for_database(muvstore_crust_mantle(1,1,1,1), size(muvstore_crust_mantle))
  call bcast_all_cr_for_database(kappahstore_crust_mantle(1,1,1,1), size(kappahstore_crust_mantle))
  call bcast_all_cr_for_database(muhstore_crust_mantle(1,1,1,1), size(muhstore_crust_mantle))
  call bcast_all_cr_for_database(eta_anisostore_crust_mantle(1,1,1,1), size(eta_anisostore_crust_mantle))
  call bcast_all_cr_for_database(c11store_crust_mantle(1,1,1,1), size(c11store_crust_mantle))
  call bcast_all_cr_for_database(c12store_crust_mantle(1,1,1,1), size(c12store_crust_mantle))
  call bcast_all_cr_for_database(c13store_crust_mantle(1,1,1,1), size(c13store_crust_mantle))
  call bcast_all_cr_for_database(c14store_crust_mantle(1,1,1,1), size(c14store_crust_mantle))
  call bcast_all_cr_for_database(c15store_crust_mantle(1,1,1,1), size(c15store_crust_mantle))
  call bcast_all_cr_for_database(c16store_crust_mantle(1,1,1,1), size(c16store_crust_mantle))
  call bcast_all_cr_for_database(c22store_crust_mantle(1,1,1,1), size(c22store_crust_mantle))
  call bcast_all_cr_for_database(c23store_crust_mantle(1,1,1,1), size(c23store_crust_mantle))
  call bcast_all_cr_for_database(c24store_crust_mantle(1,1,1,1), size(c24store_crust_mantle))
  call bcast_all_cr_for_database(c25store_crust_mantle(1,1,1,1), size(c25store_crust_mantle))
  call bcast_all_cr_for_database(c26store_crust_mantle(1,1,1,1), size(c26store_crust_mantle))
  call bcast_all_cr_for_database(c33store_crust_mantle(1,1,1,1), size(c33store_crust_mantle))
  call bcast_all_cr_for_database(c34store_crust_mantle(1,1,1,1), size(c34store_crust_mantle))
  call bcast_all_cr_for_database(c35store_crust_mantle(1,1,1,1), size(c35store_crust_mantle))
  call bcast_all_cr_for_database(c36store_crust_mantle(1,1,1,1), size(c36store_crust_mantle))
  call bcast_all_cr_for_database(c44store_crust_mantle(1,1,1,1), size(c44store_crust_mantle))
  call bcast_all_cr_for_database(c45store_crust_mantle(1,1,1,1), size(c45store_crust_mantle))
  call bcast_all_cr_for_database(c46store_crust_mantle(1,1,1,1), size(c46store_crust_mantle))
  call bcast_all_cr_for_database(c55store_crust_mantle(1,1,1,1), size(c55store_crust_mantle))
  call bcast_all_cr_for_database(c56store_crust_mantle(1,1,1,1), size(c56store_crust_mantle))
  call bcast_all_cr_for_database(c66store_crust_mantle(1,1,1,1), size(c66store_crust_mantle))
  call bcast_all_i_for_database(ibool_crust_mantle(1,1,1,1), size(ibool_crust_mantle))
  call bcast_all_l_for_database(ispec_is_tiso_crust_mantle(1), size(ispec_is_tiso_crust_mantle))
  call bcast_all_cr_for_database(rmassx_crust_mantle(1), size(rmassx_crust_mantle))
  call bcast_all_cr_for_database(rmassy_crust_mantle(1), size(rmassy_crust_mantle))
  call bcast_all_cr_for_database(rmassz_crust_mantle(1), size(rmassz_crust_mantle))
  !call bcast_all_i_for_database(NGLOB_CRUST_MANTLE_OCEANS, 1)
  call bcast_all_cr_for_database(rmass_ocean_load(1), size(rmass_ocean_load))
  call bcast_all_cr_for_database(b_rmassx_crust_mantle(1), size(b_rmassx_crust_mantle))
  call bcast_all_cr_for_database(b_rmassy_crust_mantle(1), size(b_rmassy_crust_mantle))

  end subroutine bcast_mesh_databases_CM

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_OC()

  use specfem_par
  use specfem_par_outercore
  implicit none

  !call bcast_all_i_for_database(NSPEC_OUTER_CORE, 1)
  !call bcast_all_i_for_database(NGLOB_OUTER_CORE, 1)
  call bcast_all_cr_for_database(vp_outer_core(1,1,1,1), size(vp_outer_core))
  call bcast_all_cr_for_database(xstore_outer_core(1), size(xstore_outer_core))
  call bcast_all_cr_for_database(ystore_outer_core(1), size(ystore_outer_core))
  call bcast_all_cr_for_database(zstore_outer_core(1), size(zstore_outer_core))
  call bcast_all_cr_for_database(xix_outer_core(1,1,1,1), size(xix_outer_core))
  call bcast_all_cr_for_database(xiy_outer_core(1,1,1,1), size(xiy_outer_core))
  call bcast_all_cr_for_database(xiz_outer_core(1,1,1,1), size(xiz_outer_core))
  call bcast_all_cr_for_database(etax_outer_core(1,1,1,1), size(etax_outer_core))
  call bcast_all_cr_for_database(etay_outer_core(1,1,1,1), size(etay_outer_core))
  call bcast_all_cr_for_database(etaz_outer_core(1,1,1,1), size(etaz_outer_core))
  call bcast_all_cr_for_database(gammax_outer_core(1,1,1,1), size(gammax_outer_core))
  call bcast_all_cr_for_database(gammay_outer_core(1,1,1,1), size(gammay_outer_core))
  call bcast_all_cr_for_database(gammaz_outer_core(1,1,1,1), size(gammaz_outer_core))
  call bcast_all_cr_for_database(rhostore_outer_core(1,1,1,1), size(rhostore_outer_core))
  call bcast_all_cr_for_database(kappavstore_outer_core(1,1,1,1), size(kappavstore_outer_core))
  call bcast_all_i_for_database(ibool_outer_core(1,1,1,1), size(ibool_outer_core))
  call bcast_all_cr_for_database(rmass_outer_core(1), size(rmass_outer_core))

  end subroutine bcast_mesh_databases_OC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_IC()

  use specfem_par
  use specfem_par_innercore
  implicit none

  !call bcast_all_i_for_database(NSPEC_INNER_CORE, 1)
  !call bcast_all_i_for_database(NGLOB_INNER_CORE, 1)
  !call bcast_all_i_for_database(NGLOB_XY_IC, 1)
  call bcast_all_cr_for_database(xstore_inner_core(1), size(xstore_inner_core))
  call bcast_all_cr_for_database(ystore_inner_core(1), size(ystore_inner_core))
  call bcast_all_cr_for_database(zstore_inner_core(1), size(zstore_inner_core))
  call bcast_all_cr_for_database(xix_inner_core(1,1,1,1), size(xix_inner_core))
  call bcast_all_cr_for_database(xiy_inner_core(1,1,1,1), size(xiy_inner_core))
  call bcast_all_cr_for_database(xiz_inner_core(1,1,1,1), size(xiz_inner_core))
  call bcast_all_cr_for_database(etax_inner_core(1,1,1,1), size(etax_inner_core))
  call bcast_all_cr_for_database(etay_inner_core(1,1,1,1), size(etay_inner_core))
  call bcast_all_cr_for_database(etaz_inner_core(1,1,1,1), size(etaz_inner_core))
  call bcast_all_cr_for_database(gammax_inner_core(1,1,1,1), size(gammax_inner_core))
  call bcast_all_cr_for_database(gammay_inner_core(1,1,1,1), size(gammay_inner_core))
  call bcast_all_cr_for_database(gammaz_inner_core(1,1,1,1), size(gammaz_inner_core))
  call bcast_all_cr_for_database(rhostore_inner_core(1,1,1,1), size(rhostore_inner_core))
  call bcast_all_cr_for_database(kappavstore_inner_core(1,1,1,1), size(kappavstore_inner_core))
  call bcast_all_cr_for_database(muvstore_inner_core(1,1,1,1), size(muvstore_inner_core))
  call bcast_all_cr_for_database(c11store_inner_core(1,1,1,1), size(c11store_inner_core))
  call bcast_all_cr_for_database(c12store_inner_core(1,1,1,1), size(c12store_inner_core))
  call bcast_all_cr_for_database(c13store_inner_core(1,1,1,1), size(c13store_inner_core))
  call bcast_all_cr_for_database(c33store_inner_core(1,1,1,1), size(c33store_inner_core))
  call bcast_all_cr_for_database(c44store_inner_core(1,1,1,1), size(c44store_inner_core))
  call bcast_all_i_for_database(ibool_inner_core(1,1,1,1), size(ibool_inner_core))
  call bcast_all_i_for_database(idoubling_inner_core(1), size(idoubling_inner_core))
  call bcast_all_cr_for_database(rmassx_inner_core(1), size(rmassx_inner_core))
  call bcast_all_cr_for_database(rmassy_inner_core(1), size(rmassy_inner_core))
  call bcast_all_cr_for_database(rmassz_inner_core(1), size(rmassz_inner_core))
  call bcast_all_cr_for_database(b_rmassx_inner_core(1), size(b_rmassx_inner_core))
  call bcast_all_cr_for_database(b_rmassy_inner_core(1), size(b_rmassy_inner_core))

  end subroutine bcast_mesh_databases_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_coupling()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  call bcast_all_i_for_database(nspec2D_xmin_crust_mantle, 1)
  call bcast_all_i_for_database(nspec2D_xmax_crust_mantle, 1)
  call bcast_all_i_for_database(nspec2D_ymin_crust_mantle, 1)
  call bcast_all_i_for_database(nspec2D_ymax_crust_mantle, 1)

  call bcast_all_i_for_database(ibelm_xmin_crust_mantle(1), size(ibelm_xmin_crust_mantle))
  call bcast_all_i_for_database(ibelm_xmax_crust_mantle(1), size(ibelm_xmax_crust_mantle))
  call bcast_all_i_for_database(ibelm_ymin_crust_mantle(1), size(ibelm_ymin_crust_mantle))
  call bcast_all_i_for_database(ibelm_ymax_crust_mantle(1), size(ibelm_ymax_crust_mantle))
  call bcast_all_i_for_database(ibelm_bottom_crust_mantle(1), size(ibelm_bottom_crust_mantle))
  call bcast_all_i_for_database(ibelm_top_crust_mantle(1), size(ibelm_top_crust_mantle))

  call bcast_all_cr_for_database(normal_xmin_crust_mantle(1,1,1,1), size(normal_xmin_crust_mantle))
  call bcast_all_cr_for_database(normal_xmax_crust_mantle(1,1,1,1), size(normal_xmax_crust_mantle))
  call bcast_all_cr_for_database(normal_ymin_crust_mantle(1,1,1,1), size(normal_ymin_crust_mantle))
  call bcast_all_cr_for_database(normal_ymax_crust_mantle(1,1,1,1), size(normal_ymax_crust_mantle))
  call bcast_all_cr_for_database(normal_bottom_crust_mantle(1,1,1,1), size(normal_bottom_crust_mantle))
  call bcast_all_cr_for_database(normal_top_crust_mantle(1,1,1,1), size(normal_top_crust_mantle))

  call bcast_all_cr_for_database(jacobian2D_xmin_crust_mantle(1,1,1), size(jacobian2D_xmin_crust_mantle))
  call bcast_all_cr_for_database(jacobian2D_xmax_crust_mantle(1,1,1), size(jacobian2D_xmax_crust_mantle))
  call bcast_all_cr_for_database(jacobian2D_ymin_crust_mantle(1,1,1), size(jacobian2D_ymin_crust_mantle))
  call bcast_all_cr_for_database(jacobian2D_ymax_crust_mantle(1,1,1), size(jacobian2D_ymax_crust_mantle))
  call bcast_all_cr_for_database(jacobian2D_bottom_crust_mantle(1,1,1), size(jacobian2D_bottom_crust_mantle))
  call bcast_all_cr_for_database(jacobian2D_top_crust_mantle(1,1,1), size(jacobian2D_top_crust_mantle))

  call bcast_all_i_for_database(nspec2D_xmin_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_xmax_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_ymin_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_ymax_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_zmin_outer_core, 1)

  call bcast_all_i_for_database(ibelm_xmin_outer_core(1), size(ibelm_xmin_outer_core))
  call bcast_all_i_for_database(ibelm_xmax_outer_core(1), size(ibelm_xmax_outer_core))
  call bcast_all_i_for_database(ibelm_ymin_outer_core(1), size(ibelm_ymin_outer_core))
  call bcast_all_i_for_database(ibelm_ymax_outer_core(1), size(ibelm_ymax_outer_core))
  call bcast_all_i_for_database(ibelm_bottom_outer_core(1), size(ibelm_bottom_outer_core))
  call bcast_all_i_for_database(ibelm_top_outer_core(1), size(ibelm_top_outer_core))

  call bcast_all_cr_for_database(normal_xmin_outer_core(1,1,1,1), size(normal_xmin_outer_core))
  call bcast_all_cr_for_database(normal_xmax_outer_core(1,1,1,1), size(normal_xmax_outer_core))
  call bcast_all_cr_for_database(normal_ymin_outer_core(1,1,1,1), size(normal_ymin_outer_core))
  call bcast_all_cr_for_database(normal_ymax_outer_core(1,1,1,1), size(normal_ymax_outer_core))
  call bcast_all_cr_for_database(normal_bottom_outer_core(1,1,1,1), size(normal_bottom_outer_core))
  call bcast_all_cr_for_database(normal_top_outer_core(1,1,1,1), size(normal_top_outer_core))

  call bcast_all_cr_for_database(jacobian2D_xmin_outer_core(1,1,1), size(jacobian2D_xmin_outer_core))
  call bcast_all_cr_for_database(jacobian2D_xmax_outer_core(1,1,1), size(jacobian2D_xmax_outer_core))
  call bcast_all_cr_for_database(jacobian2D_ymin_outer_core(1,1,1), size(jacobian2D_ymin_outer_core))
  call bcast_all_cr_for_database(jacobian2D_ymax_outer_core(1,1,1), size(jacobian2D_ymax_outer_core))
  call bcast_all_cr_for_database(jacobian2D_bottom_outer_core(1,1,1), size(jacobian2D_bottom_outer_core))
  call bcast_all_cr_for_database(jacobian2D_top_outer_core(1,1,1), size(jacobian2D_top_outer_core))

  call bcast_all_i_for_database(nspec2D_xmin_inner_core, 1)
  call bcast_all_i_for_database(nspec2D_xmax_inner_core, 1)
  call bcast_all_i_for_database(nspec2D_ymin_inner_core, 1)
  call bcast_all_i_for_database(nspec2D_ymax_inner_core, 1)

    ! boundary parameters
  call bcast_all_i_for_database(ibelm_xmin_inner_core(1), size(ibelm_xmin_inner_core))
  call bcast_all_i_for_database(ibelm_xmax_inner_core(1), size(ibelm_xmax_inner_core))
  call bcast_all_i_for_database(ibelm_ymin_inner_core(1), size(ibelm_ymin_inner_core))
  call bcast_all_i_for_database(ibelm_ymax_inner_core(1), size(ibelm_ymax_inner_core))
  call bcast_all_i_for_database(ibelm_bottom_inner_core(1), size(ibelm_bottom_inner_core))
  call bcast_all_i_for_database(ibelm_top_inner_core(1), size(ibelm_top_inner_core))

  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
    call bcast_all_i_for_database(ibelm_moho_top(1), size(ibelm_moho_top))
    call bcast_all_i_for_database(ibelm_moho_bot(1), size(ibelm_moho_bot))
    call bcast_all_i_for_database(ibelm_400_top(1), size(ibelm_400_top))
    call bcast_all_i_for_database(ibelm_400_bot(1), size(ibelm_400_bot))
    call bcast_all_i_for_database(ibelm_670_top(1), size(ibelm_670_top))
    call bcast_all_i_for_database(ibelm_670_bot(1), size(ibelm_670_bot))
    call bcast_all_cr_for_database(normal_moho(1,1,1,1), size(normal_moho))
    call bcast_all_cr_for_database(normal_400(1,1,1,1), size(normal_400))
    call bcast_all_cr_for_database(normal_670(1,1,1,1), size(normal_670))
  endif

  end subroutine bcast_mesh_databases_coupling

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_CM()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_crust_mantle, 1)

  ! could also test for not allocated, only reader processes have
  ! allocated these arrays.
  if (.not. I_should_read_the_database) then
    allocate(my_neighbours_crust_mantle(num_interfaces_crust_mantle), &
            nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating array my_neighbours_crust_mantle etc.')
  endif
  if (num_interfaces_crust_mantle > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_cm, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
              stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
    endif

    call bcast_all_i_for_database(my_neighbours_crust_mantle(1), size(my_neighbours_crust_mantle))
    call bcast_all_i_for_database(nibool_interfaces_crust_mantle(1), size(nibool_interfaces_crust_mantle))
    call bcast_all_i_for_database(ibool_interfaces_crust_mantle(1,1), size(ibool_interfaces_crust_mantle))
  else
    ! dummy array
    max_nibool_interfaces_cm = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_crust_mantle')
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_crust_mantle, 1)
  call bcast_all_i_for_database(nspec_outer_crust_mantle, 1)
  call bcast_all_i_for_database(num_phase_ispec_crust_mantle, 1)
  if (num_phase_ispec_crust_mantle < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_crust_mantle is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),&
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')
  endif

  if (num_phase_ispec_crust_mantle > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_crust_mantle(1,1), size(phase_ispec_inner_crust_mantle))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_crust_mantle, 1)
    call bcast_all_i_for_database(num_colors_inner_crust_mantle, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
              stat=ier)
      if (ier /= 0 ) &
        call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
    endif

    call bcast_all_i_for_database(num_elem_colors_crust_mantle(1), size(num_elem_colors_crust_mantle))
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
              stat=ier)
      if (ier /= 0 ) &
        call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_CM

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_OC()

  use specfem_par
  use specfem_par_outercore
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_outer_core, 1)
  if (.not. I_should_read_the_database) then
    allocate(my_neighbours_outer_core(num_interfaces_outer_core), &
            nibool_interfaces_outer_core(num_interfaces_outer_core), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating array my_neighbours_outer_core etc.')
  endif

  if (num_interfaces_outer_core > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_oc, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
              stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_outer_core')
    endif
  endif

  if (num_interfaces_outer_core > 0) then
    call bcast_all_i_for_database(my_neighbours_outer_core(1), size(my_neighbours_outer_core))
    call bcast_all_i_for_database(nibool_interfaces_outer_core(1), size(nibool_interfaces_outer_core))
    call bcast_all_i_for_database(ibool_interfaces_outer_core(1,1), size(ibool_interfaces_outer_core))
  else
    ! dummy array
    max_nibool_interfaces_oc = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_outer_core(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_outer_core')
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_outer_core, 1)
  call bcast_all_i_for_database(nspec_outer_outer_core, 1)
  call bcast_all_i_for_database(num_phase_ispec_outer_core, 1)
  if (num_phase_ispec_outer_core < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_outer_core is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),&
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')
  endif

  if (num_phase_ispec_outer_core > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_outer_core(1,1), size(phase_ispec_inner_outer_core))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_outer_core, 1)
    call bcast_all_i_for_database(num_colors_inner_outer_core, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
              stat=ier)
      if (ier /= 0 ) &
        call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
    endif

    call bcast_all_i_for_database(num_elem_colors_outer_core(1), size(num_elem_colors_outer_core))
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
              stat=ier)
      if (ier /= 0 ) &
        call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_OC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_IC()

  use specfem_par
  use specfem_par_innercore
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_inner_core, 1)
  if (.not. I_should_read_the_database) then
    allocate(my_neighbours_inner_core(num_interfaces_inner_core), &
            nibool_interfaces_inner_core(num_interfaces_inner_core), &
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating array my_neighbours_inner_core etc.')
  endif

  if (num_interfaces_inner_core > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_ic, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic,num_interfaces_inner_core), &
              stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_inner_core')
    endif

    call bcast_all_i_for_database(my_neighbours_inner_core(1), size(my_neighbours_inner_core))
    call bcast_all_i_for_database(nibool_interfaces_inner_core(1), size(nibool_interfaces_inner_core))
    call bcast_all_i_for_database(ibool_interfaces_inner_core(1,1), size(ibool_interfaces_inner_core))
  else
    ! dummy array
    max_nibool_interfaces_ic = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_inner_core(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_inner_core')
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_inner_core, 1)
  call bcast_all_i_for_database(nspec_outer_inner_core, 1)
  call bcast_all_i_for_database(num_phase_ispec_inner_core, 1)
  if (num_phase_ispec_inner_core < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_inner_core is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),&
            stat=ier)
    if (ier /= 0 ) &
      call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')
  endif

  if (num_phase_ispec_inner_core > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_inner_core(1,1), size(phase_ispec_inner_inner_core))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_inner_core, 1)
    call bcast_all_i_for_database(num_colors_inner_inner_core, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
              stat=ier)
      if (ier /= 0 ) &
        call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
    endif

    call bcast_all_i_for_database(num_elem_colors_inner_core(1), size(num_elem_colors_inner_core))
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
              stat=ier)
      if (ier /= 0 ) &
        call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_stacey()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! crust and mantle
  call bcast_all_i_for_database(nimin_crust_mantle(1,1), size(nimin_crust_mantle))
  call bcast_all_i_for_database(nimax_crust_mantle(1,1), size(nimax_crust_mantle))
  call bcast_all_i_for_database(njmin_crust_mantle(1,1), size(njmin_crust_mantle))
  call bcast_all_i_for_database(njmax_crust_mantle(1,1), size(njmax_crust_mantle))
  call bcast_all_i_for_database(nkmin_xi_crust_mantle(1,1), size(nkmin_xi_crust_mantle))
  call bcast_all_i_for_database(nkmin_eta_crust_mantle(1,1), size(nkmin_eta_crust_mantle))

  ! outer core
  call bcast_all_i_for_database(nimin_outer_core(1,1), size(nimin_outer_core))
  call bcast_all_i_for_database(nimax_outer_core(1,1), size(nimax_outer_core))
  call bcast_all_i_for_database(njmin_outer_core(1,1), size(njmin_outer_core))
  call bcast_all_i_for_database(njmax_outer_core(1,1), size(njmax_outer_core))
  call bcast_all_i_for_database(nkmin_xi_outer_core(1,1), size(nkmin_xi_outer_core))
  call bcast_all_i_for_database(nkmin_eta_outer_core(1,1), size(nkmin_eta_outer_core))

  end subroutine bcast_mesh_databases_stacey
