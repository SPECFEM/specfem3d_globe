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
  if (SYNC_READING) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in crust/mantle databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_CM()

  ! reads "solver_data.bin" files for outer core
  if (SYNC_READING) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in outer core databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_OC()

  ! reads "solver_data.bin" files for inner core
  if (SYNC_READING) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in inner core databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_IC()

  ! full gravity arrays
  if (FULL_GRAVITY_VAL) then
    ! transition-to-infinite region
    if (ADD_TRINF) then
      if (SYNC_READING) call synchronize_all()
      if (myrank == 0) then
        write(IMAIN,*) '  reading in transition-to-infinite databases...'
        call flush_IMAIN()
      endif
      call read_mesh_databases_TRINF()
    endif

    ! infinite region
    if (SYNC_READING) call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  reading in infinite databases...'
      call flush_IMAIN()
    endif
    call read_mesh_databases_INF()
  endif

  ! reads "boundary.bin" files to couple mantle with outer core and inner core boundaries
  if (SYNC_READING) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in coupling surface databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_coupling()

  ! sets up MPI interfaces, inner/outer elements and mesh coloring
  if (SYNC_READING) call synchronize_all()
  if (myrank == 0) then
    write(IMAIN,*) '  reading in MPI databases...'
    call flush_IMAIN()
  endif
  call read_mesh_databases_MPI()

  ! absorbing boundaries
  if (ABSORBING_CONDITIONS) then
    ! reads "stacey.bin" files
    if (SYNC_READING) call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  reading in Stacey databases...'
      call flush_IMAIN()
    endif
    call read_mesh_databases_stacey()
  endif

  ! kernels on regular grids
  if (SAVE_REGULAR_KL .and. SIMULATION_TYPE == 3) then
    if (SYNC_READING) call synchronize_all()
    if (myrank == 0) then
      write(IMAIN,*) '  reading in regular kernel databases...'
      call flush_IMAIN()
    endif
    call read_mesh_databases_regular_kl()
  endif

#ifdef USE_SERIAL_CASCADE_FOR_IOs
  ! serial i/o
  you_can_start_doing_IOs = .true.
  if (myrank < NPROC_XI_VAL*NPROC_ETA_VAL-1) call send_singlel(you_can_start_doing_IOs,myrank+1,itag)
#endif

  ! frees arrays
  if (allocated(addressing)) deallocate(addressing)

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
  ! dummy array that does not need to be actually read
  integer, dimension(:),allocatable :: dummy_idoubling
  integer :: ier

  ! crust and mantle
  nspec_iso = NSPEC_CRUST_MANTLE
  if (ANISOTROPIC_3D_MANTLE_VAL) then
    nspec_tiso = NSPECMAX_TISO_MANTLE ! 1
    nspec_ani = NSPEC_CRUST_MANTLE
  else
    if (TRANSVERSE_ISOTROPY_VAL) then
      nspec_tiso = NSPECMAX_TISO_MANTLE
    else
      nspec_tiso = 0
    endif
    nspec_ani = NSPECMAX_ANISO_MANTLE ! 1
  endif

  ! sets number of top elements for surface movies & noise tomography
  NSPEC_TOP = NSPEC2D_TOP(IREGION_CRUST_MANTLE)

  ! allocates dummy array
  allocate(dummy_idoubling(NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy idoubling in crust_mantle')
  dummy_idoubling(:) = 0

  ! allocates mass matrices in this slice (will be fully assembled in the solver)
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete

  allocate(rmassz_crust_mantle(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating rmassz in crust_mantle'
  rmassz_crust_mantle(:) = 0.0_CUSTOM_REAL

  ! allocates mass matrices
  allocate(rmassx_crust_mantle(NGLOB_XY_CM), &
           rmassy_crust_mantle(NGLOB_XY_CM),stat=ier)
  if (ier /= 0) stop 'Error allocating rmassx, rmassy in crust_mantle'
  rmassx_crust_mantle(:) = 0.0_CUSTOM_REAL
  rmassy_crust_mantle(:) = 0.0_CUSTOM_REAL

  ! b_rmassx and b_rmassy will be different to rmassx and rmassy
  ! needs new arrays
  allocate(b_rmassx_crust_mantle(NGLOB_XY_CM), &
           b_rmassy_crust_mantle(NGLOB_XY_CM),stat=ier)
  if (ier /= 0) stop 'Error allocating b_rmassx, b_rmassy in crust_mantle'
  b_rmassx_crust_mantle(:) = 0.0_CUSTOM_REAL
  b_rmassy_crust_mantle(:) = 0.0_CUSTOM_REAL

  ! x/y/z locations
  allocate(xstore_crust_mantle(NGLOB_CRUST_MANTLE), &
           ystore_crust_mantle(NGLOB_CRUST_MANTLE), &
           zstore_crust_mantle(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating x/y/zstore in crust_mantle'
  xstore_crust_mantle(:) = 0.0_CUSTOM_REAL
  ystore_crust_mantle(:) = 0.0_CUSTOM_REAL
  zstore_crust_mantle(:) = 0.0_CUSTOM_REAL

  allocate(rmass_ocean_load(NGLOB_CRUST_MANTLE_OCEANS),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rmass_ocean_load,..'
  rmass_ocean_load(:) = 0.0_CUSTOM_REAL

  allocate(ibool_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           xix_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           xiy_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           xiz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           etax_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           etay_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           etaz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           gammax_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           gammay_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           gammaz_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibool_crust_mantle,..'
  ibool_crust_mantle(:,:,:,:) = 0
  xix_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  xiy_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  xiz_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  etax_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  etay_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  etaz_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  gammax_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  gammay_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  gammaz_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(rhostore_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           kappavstore_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           muvstore_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rhostore_crust_mantle,..'
  rhostore_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  kappavstore_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  muvstore_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(kappahstore_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
           muhstore_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
           eta_anisostore_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays kappahstore_crust_mantle,..'
  kappahstore_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  muhstore_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  eta_anisostore_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  if (.not. ANISOTROPIC_3D_MANTLE_VAL) then
    ! allocates c11stores,.. for tiso elements
    allocate(c11store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c12store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c13store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c14store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c15store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c16store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c22store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c23store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c24store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c25store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c26store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c33store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c34store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c35store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c36store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c44store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c45store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c46store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c55store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c56store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE), &
             c66store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays c11store_crust_mantle,..'
  else
    ! allocates c11stores,.. for aniso elements
    allocate(c11store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c12store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c13store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c14store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c15store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c16store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c22store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c23store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c24store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c25store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c26store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c33store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c34store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c35store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c36store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c44store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c45store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c46store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c55store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c56store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE), &
             c66store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays c11store_crust_mantle,..'
  endif
  c11store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c12store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c13store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c14store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c15store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c16store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c22store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c23store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c24store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c25store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c26store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c33store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c34store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c35store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c36store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c44store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c45store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c46store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c55store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL; c56store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  c66store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  ! for azimuthal
  allocate(mu0store_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating mu0 array'
  mu0store_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  allocate(ispec_is_tiso_crust_mantle(NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating array ispec_is_tiso_crust_mantle'
  ispec_is_tiso_crust_mantle(:) = .false.

  allocate(rho_vp_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY), &
           rho_vs_crust_mantle(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rho_vp_crust_mantle,..'
  rho_vp_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL
  rho_vs_crust_mantle(:,:,:,:) = 0.0_CUSTOM_REAL

  ! reads databases file
  if (I_should_read_the_database .and. NSPEC_CRUST_MANTLE > 0) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_CRUST_MANTLE, &
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
                                    mu0store_crust_mantle, &
                                    ibool_crust_mantle,dummy_idoubling,ispec_is_tiso_crust_mantle, &
                                    rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                                    NGLOB_CRUST_MANTLE_OCEANS,rmass_ocean_load, &
                                    b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    else
      call read_arrays_solver(IREGION_CRUST_MANTLE, &
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
                              mu0store_crust_mantle, &
                              ibool_crust_mantle,dummy_idoubling,ispec_is_tiso_crust_mantle, &
                              rmassx_crust_mantle,rmassy_crust_mantle,rmassz_crust_mantle, &
                              NGLOB_CRUST_MANTLE_OCEANS,rmass_ocean_load, &
                              b_rmassx_crust_mantle,b_rmassy_crust_mantle)
    endif
  endif
  call bcast_mesh_databases_CM()

  ! check that the number of points in this slice is correct
  if (NSPEC_CRUST_MANTLE > 0) then
    if (minval(ibool_crust_mantle(:,:,:,:)) /= 1) then
      print *,'Error: rank ',myrank,' has invalid crust_mantle ibool indexing min/max = ', &
              minval(ibool_crust_mantle(:,:,:,:)),maxval(ibool_crust_mantle(:,:,:,:)),'max glob = ',NGLOB_CRUST_MANTLE
      call exit_MPI(myrank,'incorrect global numbering: iboolmin is not equal to 1 in crust and mantle')
    endif
    if (maxval(ibool_crust_mantle(:,:,:,:)) /= NGLOB_CRUST_MANTLE) then
      print *,'Error: rank ',myrank,' has invalid ibool indexing min/max = ', &
              minval(ibool_crust_mantle(:,:,:,:)),maxval(ibool_crust_mantle(:,:,:,:)), 'max glob =',NGLOB_CRUST_MANTLE
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in crust and mantle')
    endif
  endif

  deallocate(dummy_idoubling)

  ! mass matrix corrections
  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL)) then
    ! mass matrices differ for rmassx,rmassy in case Newmark time scheme is used
    ! note: uses different, allocated rmassx/y** arrays
    !       we will also use them for LDDRK schemes even though they don't need an additional mass term for absorbing conditions
    continue
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
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      ! mass matrices differ for b_rmassx,b_rmassy in case Newmark time scheme is used
      ! note: uses different, allocated rmassx/y** arrays
      !       we will also use them for LDDRK schemes even though they don't need an additional mass term for absorbing conditions
      continue
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
    ! re-associates with corresponding rmassx,rmassy
    ! (uses dummy pointers used for passing as function arguments.
    !  associates mass matrix used for backward/reconstructed wavefields)
    b_rmassx_crust_mantle => rmassx_crust_mantle(:)
    b_rmassy_crust_mantle => rmassy_crust_mantle(:)
    nullify(b_rmassz_crust_mantle)
  endif

  ! only needed for azimuthal kernels
  if (.not. (SIMULATION_TYPE == 3 .and. (ANISOTROPIC_KL .and. SAVE_AZIMUTHAL_ANISO_KL_ONLY))) then
    deallocate(mu0store_crust_mantle)
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
  integer :: ier

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: dummy_rmass

  logical, dimension(:), allocatable :: dummy_ispec_is_tiso
  integer, dimension(:), allocatable :: dummy_idoubling_outer_core

  ! outer core (no anisotropy nor S velocity)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the outer core
  nspec_iso = NSPEC_OUTER_CORE
  nspec_tiso = 0
  nspec_ani = 0

  ! dummy allocation
  NGLOB_XY_dummy = 0

  allocate(dummy_rmass(NGLOB_XY_dummy), &
           dummy_ispec_is_tiso(NSPEC_OUTER_CORE), &
           dummy_idoubling_outer_core(NSPEC_OUTER_CORE), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating dummy rmass and dummy ispec/idoubling in outer core'
  dummy_rmass(:) = 0.0; dummy_ispec_is_tiso(:) = .false.; dummy_idoubling_outer_core(:) = 0

  ! x/y/z locations
  allocate(xstore_outer_core(NGLOB_OUTER_CORE), &
           ystore_outer_core(NGLOB_OUTER_CORE), &
           zstore_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating x/y/zstore in outer core'
  xstore_outer_core(:) = 0.0; ystore_outer_core(:) = 0.0; zstore_outer_core(:) = 0.0

  allocate(ibool_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           xix_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           xiy_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           xiz_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           etax_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           etay_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           etaz_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           gammax_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           gammay_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           gammaz_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibool_outer_core,..'
  ibool_outer_core(:,:,:,:) = 0
  xix_outer_core(:,:,:,:) = 0.0; xiy_outer_core(:,:,:,:) = 0.0; xiz_outer_core(:,:,:,:) = 0.0
  etax_outer_core(:,:,:,:) = 0.0; etay_outer_core(:,:,:,:) = 0.0; etaz_outer_core(:,:,:,:) = 0.0
  gammax_outer_core(:,:,:,:) = 0.0; gammay_outer_core(:,:,:,:) = 0.0; gammaz_outer_core(:,:,:,:) = 0.0

  allocate(rhostore_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
           kappavstore_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rhostore_outer_core,..'
  rhostore_outer_core(:,:,:,:) = 0.0; kappavstore_outer_core(:,:,:,:) = 0.0

  allocate(rmass_outer_core(NGLOB_OUTER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating array rmass_outer_core'
  rmass_outer_core(:) = 0.0

  allocate(vp_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_STACEY),stat=ier)
  if (ier /= 0) stop 'Error allocating array vp_outer_core'
  vp_outer_core(:,:,:,:) = 0.0

  ! reads in mesh arrays
  if (I_should_read_the_database .and. NSPEC_OUTER_CORE > 0) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_OUTER_CORE, &
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
                                    dummy_array, &
                                    ibool_outer_core,dummy_idoubling_outer_core,dummy_ispec_is_tiso, &
                                    dummy_rmass,dummy_rmass,rmass_outer_core, &
                                    1,dummy_array, &
                                    dummy_rmass,dummy_rmass)
    else
      call read_arrays_solver(IREGION_OUTER_CORE, &
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
                              dummy_array, &
                              ibool_outer_core,dummy_idoubling_outer_core,dummy_ispec_is_tiso, &
                              dummy_rmass,dummy_rmass,rmass_outer_core, &
                              1, dummy_array, &
                              dummy_rmass,dummy_rmass)
    endif
  endif
  call bcast_mesh_databases_OC()

  deallocate(dummy_idoubling_outer_core,dummy_ispec_is_tiso,dummy_rmass)

  ! check that the number of points in this slice is correct
  if (NSPEC_OUTER_CORE > 0) then
    if (minval(ibool_outer_core(:,:,:,:)) /= 1) then
      print *,'Error: rank ',myrank,' has invalid outer_core ibool indexing min/max = ', &
              minval(ibool_outer_core(:,:,:,:)),maxval(ibool_outer_core(:,:,:,:)),'max glob = ',NGLOB_OUTER_CORE
      call exit_MPI(myrank,'incorrect global numbering: iboolmin is not equal to 1 in outer core')
    endif
    if (maxval(ibool_outer_core(:,:,:,:)) /= NGLOB_OUTER_CORE) then
      print *,'Error: rank ',myrank,' has invalid ibool indexing min/max = ', &
              minval(ibool_outer_core(:,:,:,:)),maxval(ibool_outer_core(:,:,:,:)),'max glob = ',NGLOB_OUTER_CORE
      call exit_MPI(myrank, 'incorrect global numbering: iboolmax does not equal nglob in outer core')
    endif
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
  integer :: ier

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  logical, dimension(:),allocatable:: dummy_ispec_is_tiso

  ! inner core (no anisotropy)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the inner core
  nspec_iso = NSPEC_INNER_CORE
  nspec_tiso = 0
  if (ANISOTROPIC_INNER_CORE_VAL) then
    nspec_ani = NSPEC_INNER_CORE
  else
    nspec_ani = 0
  endif

  allocate(dummy_ispec_is_tiso(NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating dummy ispec in inner core'
  dummy_ispec_is_tiso(:) = .false.

  ! allocates mass matrices in this slice (will be fully assembled in the solver)
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  allocate(rmassz_inner_core(NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating rmassz in inner_core'
  rmassz_inner_core(:) = 0.0

  allocate(rmassx_inner_core(NGLOB_XY_IC), &
           rmassy_inner_core(NGLOB_XY_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating rmassx, rmassy in inner_core'
  rmassx_inner_core(:) = 0.0; rmassy_inner_core(:) = 0.0

  ! b_rmassx and b_rmassy maybe different to rmassx,rmassy
  allocate(b_rmassx_inner_core(NGLOB_XY_IC), &
           b_rmassy_inner_core(NGLOB_XY_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating b_rmassx, b_rmassy in inner_core'
  b_rmassx_inner_core(:) = 0.0; b_rmassy_inner_core(:) = 0.0

  ! x/y/z locations
  allocate(xstore_inner_core(NGLOB_INNER_CORE), &
           ystore_inner_core(NGLOB_INNER_CORE), &
           zstore_inner_core(NGLOB_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating x/y/zstore in inner core'
  xstore_inner_core(:) = 0.0; ystore_inner_core(:) = 0.0; zstore_inner_core(:) = 0.0

  allocate(ibool_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           xix_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           xiy_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           xiz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           etax_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           etay_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           etaz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           gammax_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           gammay_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           gammaz_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibool_inner_core,..'
  ibool_inner_core(:,:,:,:) = 0
  xix_inner_core(:,:,:,:) = 0.0; xiy_inner_core(:,:,:,:) = 0.0; xiz_inner_core(:,:,:,:) = 0.0
  etax_inner_core(:,:,:,:) = 0.0; etay_inner_core(:,:,:,:) = 0.0; etaz_inner_core(:,:,:,:) = 0.0
  gammax_inner_core(:,:,:,:) = 0.0; gammay_inner_core(:,:,:,:) = 0.0; gammaz_inner_core(:,:,:,:) = 0.0

  allocate(rhostore_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           kappavstore_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
           muvstore_inner_core(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays rhostore_inner_core,..'
  rhostore_inner_core(:,:,:,:) = 0.0; kappavstore_inner_core(:,:,:,:) = 0.0; muvstore_inner_core(:,:,:,:) = 0.0

  allocate(c11store_inner_core(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC), &
           c33store_inner_core(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC), &
           c12store_inner_core(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC), &
           c13store_inner_core(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC), &
           c44store_inner_core(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays c11store_inner_core,..'
  c11store_inner_core(:,:,:,:) = 0.0; c33store_inner_core(:,:,:,:) = 0.0; c12store_inner_core(:,:,:,:) = 0.0
  c13store_inner_core(:,:,:,:) = 0.0; c44store_inner_core(:,:,:,:) = 0.0

  allocate(idoubling_inner_core(NSPEC_INNER_CORE),stat=ier)
  if (ier /= 0) stop 'Error allocating array idoubling_inner_core'
  idoubling_inner_core(:) = 0

  ! reads in arrays
  if (I_should_read_the_database .and. NSPEC_INNER_CORE > 0) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_INNER_CORE, &
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
                                    dummy_array, &
                                    ibool_inner_core,idoubling_inner_core,dummy_ispec_is_tiso, &
                                    rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                                    1,dummy_array, &
                                    b_rmassx_inner_core,b_rmassy_inner_core)
    else
      call read_arrays_solver(IREGION_INNER_CORE, &
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
                              dummy_array, &
                              ibool_inner_core,idoubling_inner_core,dummy_ispec_is_tiso, &
                              rmassx_inner_core,rmassy_inner_core,rmassz_inner_core, &
                              1,dummy_array, &
                              b_rmassx_inner_core,b_rmassy_inner_core)
    endif
  endif
  call bcast_mesh_databases_IC()

  deallocate(dummy_ispec_is_tiso)

  ! check that the number of points in this slice is correct
  if (NSPEC_INNER_CORE > 0) then
    if (minval(ibool_inner_core(:,:,:,:)) /= 1 .or. maxval(ibool_inner_core(:,:,:,:)) /= NGLOB_INNER_CORE) then
      print *,'Error: rank ',myrank,' has invalid inner_core ibool indexing min/max = ', &
              minval(ibool_inner_core(:,:,:,:)),maxval(ibool_inner_core(:,:,:,:)),'max glob = ',NGLOB_INNER_CORE
      call exit_MPI(myrank,'incorrect global numbering: iboolmin/max does not equal 1/nglob in inner core')
    endif
  endif

  ! mass matrix corrections
  if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
    ! uses corrected mass matrices in case Newmark time scheme is used
    ! note: uses different, allocated rmassx/y** arrays
    !       we will also use them for LDDRK schemes even though they don't need an additional mass term for absorbing conditions
    continue
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
    if (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) then
      ! uses corrected mass matrices in case Newmark time scheme is used
      ! note: uses different, allocated rmassx/y** arrays
      !       we will also use them for LDDRK schemes even though they don't need an additional mass term for absorbing conditions
      continue
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
    ! re-associates with corresponding rmassx,rmassy
    ! (uses dummy pointers used for passing as function arguments.
    !  associates mass matrix used for backward/reconstructed wavefields)
    b_rmassx_inner_core => rmassx_inner_core
    b_rmassy_inner_core => rmassy_inner_core
    nullify(b_rmassz_inner_core)
  endif

  end subroutine read_mesh_databases_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_TRINF()

! mesh for TRINFINITE region

  use specfem_par
  use specfem_par_trinfinite
  implicit none

  ! local parameters
  integer :: nspec_iso,nspec_tiso,nspec_ani,NGLOB_XY_dummy
  integer :: ier

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  real(kind=CUSTOM_REAL), dimension(1) :: dummy_rmass

  logical, dimension(:), allocatable :: dummy_ispec_is_tiso
  integer, dimension(:), allocatable :: dummy_idoubling

  ! transition-to-infinite region
  ! (no isotropy, no tiso, no anisotropy)
  ! (no rmass)
  nspec_iso = 0
  nspec_tiso = 0
  nspec_ani = 0

  ! dummy allocation
  allocate(dummy_ispec_is_tiso(NSPEC_TRINFINITE), &
           dummy_idoubling(NSPEC_TRINFINITE), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating dummy rmass and dummy ispec/idoubling in trinfinite region'
  dummy_ispec_is_tiso(:) = .false.; dummy_idoubling(:) = 0

  ! x/y/z locations
  allocate(xstore_trinfinite(NGLOB_TRINFINITE), &
           ystore_trinfinite(NGLOB_TRINFINITE), &
           zstore_trinfinite(NGLOB_TRINFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating x/y/zstore in trinfinite'
  xstore_trinfinite(:) = 0.0; ystore_trinfinite(:) = 0.0; zstore_trinfinite(:) = 0.0

  allocate(ibool_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           xix_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           xiy_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           xiz_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           etax_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           etay_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           etaz_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           gammax_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           gammay_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE), &
           gammaz_trinfinite(NGLLX,NGLLY,NGLLZ,NSPEC_TRINFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibool_trinfinite,..'
  ibool_trinfinite(:,:,:,:) = 0
  xix_trinfinite(:,:,:,:) = 0.0; xiy_trinfinite(:,:,:,:) = 0.0; xiz_trinfinite(:,:,:,:) = 0.0
  etax_trinfinite(:,:,:,:) = 0.0; etay_trinfinite(:,:,:,:) = 0.0; etaz_trinfinite(:,:,:,:) = 0.0
  gammax_trinfinite(:,:,:,:) = 0.0; gammay_trinfinite(:,:,:,:) = 0.0; gammaz_trinfinite(:,:,:,:) = 0.0

  ! dummy rmass
  NGLOB_XY_dummy = 0
  dummy_rmass(:) = 0.0

  ! reads in arrays
  if (I_should_read_the_database .and. NSPEC_TRINFINITE > 0) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_TRINFINITE, &
                                    NSPEC_TRINFINITE,NGLOB_TRINFINITE,NGLOB_XY_dummy, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    dummy_array,dummy_array, &
                                    xstore_trinfinite,ystore_trinfinite,zstore_trinfinite, &
                                    xix_trinfinite,xiy_trinfinite,xiz_trinfinite, &
                                    etax_trinfinite,etay_trinfinite,etaz_trinfinite, &
                                    gammax_trinfinite,gammay_trinfinite,gammaz_trinfinite, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array, &
                                    ibool_trinfinite,dummy_idoubling,dummy_ispec_is_tiso, &
                                    dummy_rmass,dummy_rmass,dummy_rmass, &
                                    1,dummy_array, &
                                    dummy_rmass,dummy_rmass)
    else
      call read_arrays_solver(IREGION_TRINFINITE, &
                                    NSPEC_TRINFINITE,NGLOB_TRINFINITE,NGLOB_XY_dummy, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    dummy_array,dummy_array, &
                                    xstore_trinfinite,ystore_trinfinite,zstore_trinfinite, &
                                    xix_trinfinite,xiy_trinfinite,xiz_trinfinite, &
                                    etax_trinfinite,etay_trinfinite,etaz_trinfinite, &
                                    gammax_trinfinite,gammay_trinfinite,gammaz_trinfinite, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array, &
                                    ibool_trinfinite,dummy_idoubling,dummy_ispec_is_tiso, &
                                    dummy_rmass,dummy_rmass,dummy_rmass, &
                                    1,dummy_array, &
                                    dummy_rmass,dummy_rmass)
    endif
  endif
  ! broadcast to other processes
  call bcast_mesh_databases_TRINF()

  ! free dummy arrays
  deallocate(dummy_ispec_is_tiso,dummy_idoubling)

  ! check
  if (NSPEC_TRINFINITE > 0) then
    ! check that the number of points in this slice is correct
    if (minval(ibool_trinfinite(:,:,:,:)) /= 1 .or. &
        maxval(ibool_trinfinite(:,:,:,:)) /= NGLOB_TRINFINITE) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in transition infinite region')
  endif

  end subroutine read_mesh_databases_TRINF

!
!-------------------------------------------------------------------------------------------------
!

subroutine read_mesh_databases_INF()

! mesh for INFINITE region

  use specfem_par
  use specfem_par_infinite
  implicit none

  ! local parameters
  integer :: nspec_iso,nspec_tiso,nspec_ani,NGLOB_XY_dummy
  integer :: ier

  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  real(kind=CUSTOM_REAL), dimension(1) :: dummy_rmass

  logical, dimension(:), allocatable :: dummy_ispec_is_tiso
  integer, dimension(:), allocatable :: dummy_idoubling

  ! infinite region
  ! (no isotropy, no tiso, no anisotropy)
  ! (no rmass)
  nspec_iso = 0
  nspec_tiso = 0
  nspec_ani = 0

  ! dummy allocation
  allocate(dummy_ispec_is_tiso(NSPEC_INFINITE), &
           dummy_idoubling(NSPEC_INFINITE), &
           stat=ier)
  if (ier /= 0) stop 'Error allocating dummy rmass and dummy ispec/idoubling in infinite region'
  dummy_ispec_is_tiso(:) = .false.; dummy_idoubling(:) = 0

  ! x/y/z locations
  allocate(xstore_infinite(NGLOB_INFINITE), &
           ystore_infinite(NGLOB_INFINITE), &
           zstore_infinite(NGLOB_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating x/y/zstore in infinite'
  xstore_infinite(:) = 0.0; ystore_infinite(:) = 0.0; zstore_infinite(:) = 0.0

  allocate(ibool_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           xix_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           xiy_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           xiz_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           etax_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           etay_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           etaz_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           gammax_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           gammay_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE), &
           gammaz_infinite(NGLLX,NGLLY,NGLLZ,NSPEC_INFINITE),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibool_infinite,..'
  ibool_infinite(:,:,:,:) = 0
  xix_infinite(:,:,:,:) = 0.0; xiy_infinite(:,:,:,:) = 0.0; xiz_infinite(:,:,:,:) = 0.0
  etax_infinite(:,:,:,:) = 0.0; etay_infinite(:,:,:,:) = 0.0; etaz_infinite(:,:,:,:) = 0.0
  gammax_infinite(:,:,:,:) = 0.0; gammay_infinite(:,:,:,:) = 0.0; gammaz_infinite(:,:,:,:) = 0.0

  ! dummy rmass
  NGLOB_XY_dummy = 0
  dummy_rmass(:) = 0.0

  ! reads in arrays
  if (I_should_read_the_database .and. NSPEC_INFINITE > 0) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_arrays_solver_adios(IREGION_INFINITE, &
                                    NSPEC_INFINITE,NGLOB_INFINITE,NGLOB_XY_dummy, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    dummy_array,dummy_array, &
                                    xstore_infinite,ystore_infinite,zstore_infinite, &
                                    xix_infinite,xiy_infinite,xiz_infinite, &
                                    etax_infinite,etay_infinite,etaz_infinite, &
                                    gammax_infinite,gammay_infinite,gammaz_infinite, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array, &
                                    ibool_infinite,dummy_idoubling,dummy_ispec_is_tiso, &
                                    dummy_rmass,dummy_rmass,dummy_rmass, &
                                    1,dummy_array, &
                                    dummy_rmass,dummy_rmass)
    else
      call read_arrays_solver(IREGION_INFINITE, &
                                    NSPEC_INFINITE,NGLOB_INFINITE,NGLOB_XY_dummy, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    dummy_array,dummy_array, &
                                    xstore_infinite,ystore_infinite,zstore_infinite, &
                                    xix_infinite,xiy_infinite,xiz_infinite, &
                                    etax_infinite,etay_infinite,etaz_infinite, &
                                    gammax_infinite,gammay_infinite,gammaz_infinite, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array,dummy_array,dummy_array, &
                                    dummy_array, &
                                    ibool_infinite,dummy_idoubling,dummy_ispec_is_tiso, &
                                    dummy_rmass,dummy_rmass,dummy_rmass, &
                                    1,dummy_array, &
                                    dummy_rmass,dummy_rmass)
    endif
  endif
  ! broadcast to other processes
  call bcast_mesh_databases_INF()

  ! free dummy arrays
  deallocate(dummy_ispec_is_tiso,dummy_idoubling)

  ! check
  if (NSPEC_INFINITE > 0) then
    ! check that the number of points in this slice is correct
    if (minval(ibool_infinite(:,:,:,:)) /= 1 .or. &
       maxval(ibool_infinite(:,:,:,:)) /= NGLOB_INFINITE) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in infinite region')
  endif

  end subroutine read_mesh_databases_INF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_coupling()

! to couple mantle with outer core

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use specfem_par_trinfinite
  use specfem_par_infinite

  implicit none

  ! local parameters
  integer :: njunk1,njunk2,njunk3
  integer :: ier

  ! allocate arrays
  allocate(ibelm_moho_top(NSPEC2D_MOHO),ibelm_moho_bot(NSPEC2D_MOHO), &
           ibelm_400_top(NSPEC2D_400),ibelm_400_bot(NSPEC2D_400), &
           ibelm_670_top(NSPEC2D_670),ibelm_670_bot(NSPEC2D_670), &
           normal_moho(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO), &
           normal_400(NDIM,NGLLX,NGLLY,NSPEC2D_400), &
           normal_670(NDIM,NGLLX,NGLLY,NSPEC2D_670),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibelm_moho_top,..'
  ibelm_moho_top(:) = 0; ibelm_moho_bot(:) = 0
  ibelm_400_top(:) = 0; ibelm_400_bot(:) = 0
  ibelm_670_top(:) = 0; ibelm_670_bot(:) = 0
  normal_moho(:,:,:,:) = 0.0; normal_400(:,:,:,:) = 0.0; normal_670(:,:,:,:) = 0.0

  allocate(ibelm_xmin_crust_mantle(NSPEC2DMAX_XMIN_XMAX_CM),ibelm_xmax_crust_mantle(NSPEC2DMAX_XMIN_XMAX_CM), &
           ibelm_ymin_crust_mantle(NSPEC2DMAX_YMIN_YMAX_CM),ibelm_ymax_crust_mantle(NSPEC2DMAX_YMIN_YMAX_CM), &
           ibelm_bottom_crust_mantle(NSPEC2D_BOTTOM_CM),ibelm_top_crust_mantle(NSPEC2D_TOP_CM), &
           normal_xmin_crust_mantle(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM), &
           normal_xmax_crust_mantle(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM), &
           normal_ymin_crust_mantle(NDIM,NGLLX,NGLLY,NSPEC2DMAX_YMIN_YMAX_CM), &
           normal_ymax_crust_mantle(NDIM,NGLLX,NGLLY,NSPEC2DMAX_YMIN_YMAX_CM), &
           normal_bottom_crust_mantle(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_CM), &
           normal_top_crust_mantle(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_CM), &
           jacobian2D_bottom_crust_mantle(NGLLX,NGLLY,NSPEC2D_BOTTOM_CM), &
           jacobian2D_top_crust_mantle(NGLLX,NGLLY,NSPEC2D_TOP_CM), &
           jacobian2D_xmin_crust_mantle(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM), &
           jacobian2D_xmax_crust_mantle(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM), &
           jacobian2D_ymin_crust_mantle(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_CM), &
           jacobian2D_ymax_crust_mantle(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_CM),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibelm_xmin_crust_mantle,..'
  ibelm_xmin_crust_mantle(:) = 0; ibelm_xmax_crust_mantle(:) = 0
  ibelm_ymin_crust_mantle(:) = 0; ibelm_ymax_crust_mantle(:) = 0
  ibelm_bottom_crust_mantle(:) = 0; ibelm_top_crust_mantle(:) = 0
  normal_xmin_crust_mantle(:,:,:,:) = 0.0; normal_xmax_crust_mantle(:,:,:,:) = 0.0
  normal_ymin_crust_mantle(:,:,:,:) = 0.0; normal_ymax_crust_mantle(:,:,:,:) = 0.0
  normal_bottom_crust_mantle(:,:,:,:) = 0.0; normal_top_crust_mantle(:,:,:,:) = 0.0
  jacobian2D_xmin_crust_mantle(:,:,:) = 0.0; jacobian2D_xmax_crust_mantle(:,:,:) = 0.0
  jacobian2D_ymin_crust_mantle(:,:,:) = 0.0; jacobian2D_ymax_crust_mantle(:,:,:) = 0.0
  jacobian2D_bottom_crust_mantle(:,:,:) = 0.0; jacobian2D_top_crust_mantle(:,:,:) = 0.0

  allocate(ibelm_xmin_outer_core(NSPEC2DMAX_XMIN_XMAX_OC),ibelm_xmax_outer_core(NSPEC2DMAX_XMIN_XMAX_OC), &
           ibelm_ymin_outer_core(NSPEC2DMAX_YMIN_YMAX_OC),ibelm_ymax_outer_core(NSPEC2DMAX_YMIN_YMAX_OC), &
           ibelm_bottom_outer_core(NSPEC2D_BOTTOM_OC),ibelm_top_outer_core(NSPEC2D_TOP_OC), &
           normal_xmin_outer_core(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC), &
           normal_xmax_outer_core(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC), &
           normal_ymin_outer_core(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC), &
           normal_ymax_outer_core(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC), &
           normal_bottom_outer_core(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC), &
           normal_top_outer_core(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC), &
           jacobian2D_bottom_outer_core(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC), &
           jacobian2D_top_outer_core(NGLLX,NGLLY,NSPEC2D_TOP_OC), &
           jacobian2D_xmin_outer_core(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC), &
           jacobian2D_xmax_outer_core(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC), &
           jacobian2D_ymin_outer_core(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC), &
           jacobian2D_ymax_outer_core(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibelm_xmin_outer_core,..'
  ibelm_xmin_outer_core(:) = 0; ibelm_xmax_outer_core(:) = 0
  ibelm_ymin_outer_core(:) = 0; ibelm_ymax_outer_core(:) = 0
  ibelm_bottom_outer_core(:) = 0; ibelm_top_outer_core(:) = 0
  normal_xmin_outer_core(:,:,:,:) = 0.0; normal_xmax_outer_core(:,:,:,:) = 0.0
  normal_ymin_outer_core(:,:,:,:) = 0.0; normal_ymax_outer_core(:,:,:,:) = 0.0
  normal_bottom_outer_core(:,:,:,:) = 0.0; normal_top_outer_core(:,:,:,:) = 0.0
  jacobian2D_xmin_outer_core(:,:,:) = 0.0; jacobian2D_xmax_outer_core(:,:,:) = 0.0
  jacobian2D_ymin_outer_core(:,:,:) = 0.0; jacobian2D_ymax_outer_core(:,:,:) = 0.0
  jacobian2D_bottom_outer_core(:,:,:) = 0.0; jacobian2D_top_outer_core(:,:,:) = 0.0

  allocate(ibelm_xmin_inner_core(NSPEC2DMAX_XMIN_XMAX_IC),ibelm_xmax_inner_core(NSPEC2DMAX_XMIN_XMAX_IC), &
           ibelm_ymin_inner_core(NSPEC2DMAX_YMIN_YMAX_IC),ibelm_ymax_inner_core(NSPEC2DMAX_YMIN_YMAX_IC), &
           ibelm_bottom_inner_core(NSPEC2D_BOTTOM_IC),ibelm_top_inner_core(NSPEC2D_TOP_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ibelm_xmin_inner_core,..'
  ibelm_xmin_inner_core(:) = 0; ibelm_xmax_inner_core(:) = 0
  ibelm_ymin_inner_core(:) = 0; ibelm_ymax_inner_core(:) = 0
  ibelm_bottom_inner_core(:) = 0; ibelm_top_inner_core(:) = 0

  nspec2D_xmin_crust_mantle = 0; nspec2D_xmax_crust_mantle = 0
  nspec2D_ymin_crust_mantle = 0; nspec2D_ymax_crust_mantle = 0

  nspec2D_xmin_outer_core = 0; nspec2D_xmax_outer_core = 0
  nspec2D_ymin_outer_core = 0; nspec2D_ymax_outer_core = 0

  nspec2D_xmin_inner_core = 0; nspec2D_xmax_inner_core = 0
  nspec2D_ymin_inner_core = 0; nspec2D_ymax_inner_core = 0

  ! infinite regions
  if (FULL_GRAVITY) then
    ! transition-to-infinite
    if (ADD_TRINF) then
      allocate(ibelm_xmin_trinfinite(NSPEC2DMAX_XMIN_XMAX_TRINF),ibelm_xmax_trinfinite(NSPEC2DMAX_XMIN_XMAX_TRINF), &
               ibelm_ymin_trinfinite(NSPEC2DMAX_YMIN_YMAX_TRINF),ibelm_ymax_trinfinite(NSPEC2DMAX_YMIN_YMAX_TRINF), &
               ibelm_bottom_trinfinite(NSPEC2D_BOTTOM_TRINF),ibelm_top_trinfinite(NSPEC2D_TOP_TRINF),stat=ier)
      if (ier /= 0) stop 'Error allocating arrays ibelm_xmin_trinfinite,..'
      ibelm_xmin_trinfinite(:) = 0; ibelm_xmax_trinfinite(:) = 0
      ibelm_ymin_trinfinite(:) = 0; ibelm_ymax_trinfinite(:) = 0
      ibelm_bottom_trinfinite(:) = 0; ibelm_top_trinfinite(:) = 0
    endif
    nspec2D_xmin_trinfinite = 0; nspec2D_xmax_trinfinite = 0
    nspec2D_ymin_trinfinite = 0; nspec2D_ymax_trinfinite = 0

    ! infinite
    allocate(ibelm_xmin_infinite(NSPEC2DMAX_XMIN_XMAX_INF),ibelm_xmax_infinite(NSPEC2DMAX_XMIN_XMAX_INF), &
             ibelm_ymin_infinite(NSPEC2DMAX_YMIN_YMAX_INF),ibelm_ymax_infinite(NSPEC2DMAX_YMIN_YMAX_INF), &
             ibelm_bottom_infinite(NSPEC2D_BOTTOM_INF),ibelm_top_infinite(NSPEC2D_TOP_INF),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays ibelm_xmin_infinite,..'
    ibelm_xmin_infinite(:) = 0; ibelm_xmax_infinite(:) = 0
    ibelm_ymin_infinite(:) = 0; ibelm_ymax_infinite(:) = 0
    ibelm_bottom_infinite(:) = 0; ibelm_top_infinite(:) = 0

    nspec2D_xmin_infinite = 0; nspec2D_xmax_infinite = 0
    nspec2D_ymin_infinite = 0; nspec2D_ymax_infinite = 0
  endif

  ! reads in arrays
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_mesh_databases_coupling_adios()
    else
      ! crust and mantle
      if (NSPEC_CRUST_MANTLE > 0) then
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

        !if (REGIONAL_MESH_CUTOFF) &
        !  nspec2D_zmin_crust_mantle = NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)

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
      endif

      ! read parameters to couple fluid and solid regions
      !
      ! outer core
      if (NSPEC_OUTER_CORE > 0) then
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
      endif

      !
      ! inner core
      !
      if (NSPEC_INNER_CORE > 0) then
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
      endif

      ! full gravity
      if (FULL_GRAVITY) then
        if (ADD_TRINF) then
          !
          ! transition-to-infinite
          !
          if (NSPEC_TRINFINITE > 0) then
            ! create name of database
            call create_name_database(prname,myrank,IREGION_TRINFINITE,LOCAL_PATH)

            ! read info
            open(unit=IIN,file=prname(1:len_trim(prname))//'boundary.bin', &
                  status='old',form='unformatted',action='read',iostat=ier)
            if (ier /= 0 ) call exit_mpi(myrank,'Error opening trinfinite boundary.bin file')

            read(IIN) nspec2D_xmin_trinfinite
            read(IIN) nspec2D_xmax_trinfinite
            read(IIN) nspec2D_ymin_trinfinite
            read(IIN) nspec2D_ymax_trinfinite
            read(IIN) njunk1
            read(IIN) njunk2

            ! boundary parameters
            read(IIN) ibelm_xmin_trinfinite
            read(IIN) ibelm_xmax_trinfinite
            read(IIN) ibelm_ymin_trinfinite
            read(IIN) ibelm_ymax_trinfinite
            read(IIN) ibelm_bottom_trinfinite
            read(IIN) ibelm_top_trinfinite
            close(IIN)
          endif
        endif

        !
        ! infinite
        !
        if (NSPEC_INFINITE > 0) then
          ! create name of database
          call create_name_database(prname,myrank,IREGION_INFINITE,LOCAL_PATH)

          ! read info
          open(unit=IIN,file=prname(1:len_trim(prname))//'boundary.bin', &
                status='old',form='unformatted',action='read',iostat=ier)
          if (ier /= 0 ) call exit_mpi(myrank,'Error opening infinite boundary.bin file')

          read(IIN) nspec2D_xmin_infinite
          read(IIN) nspec2D_xmax_infinite
          read(IIN) nspec2D_ymin_infinite
          read(IIN) nspec2D_ymax_infinite
          read(IIN) njunk1
          read(IIN) njunk2

          ! boundary parameters
          read(IIN) ibelm_xmin_infinite
          read(IIN) ibelm_xmax_infinite
          read(IIN) ibelm_ymin_infinite
          read(IIN) ibelm_ymax_infinite
          read(IIN) ibelm_bottom_infinite
          read(IIN) ibelm_top_infinite
          close(IIN)
        endif
      endif

      ! -- Boundary Mesh for crust and mantle ---
      if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
        if (NSPEC_CRUST_MANTLE > 0) then
          call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

          open(unit=IIN,file=prname(1:len_trim(prname))//'boundary_disc.bin', &
                status='old',form='unformatted',action='read',iostat=ier)
          if (ier /= 0 ) call exit_mpi(myrank,'Error opening boundary_disc.bin file')

          read(IIN) njunk1,njunk2,njunk3

          ! checks setup
          if (njunk1 /= NSPEC2D_MOHO .or. njunk2 /= NSPEC2D_400 .or. njunk3 /= NSPEC2D_670) then
            print *,'Error: invalid NSPEC2D values read in for solver: ',njunk1,njunk2,njunk3,'(boundary_disc.bin)'
            print *,'       should be MOHO/400/670 : ',NSPEC2D_MOHO,NSPEC2D_400,NSPEC2D_670,'(mesh_parameters.bin)'
            call exit_mpi(myrank, 'Error reading boundary_disc.bin file')
          endif

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

  if (FULL_GRAVITY) then
    if (ADD_TRINF) then
      ! transition-to-infinite
      if (nspec2D_xmin_trinfinite < 0 .or. nspec2d_xmin_trinfinite > NSPEC2DMAX_XMIN_XMAX_TRINF .or. &
          nspec2D_xmax_trinfinite < 0 .or. nspec2d_xmax_trinfinite > NSPEC2DMAX_XMIN_XMAX_TRINF .or. &
          nspec2D_ymin_trinfinite < 0 .or. nspec2d_ymin_trinfinite > NSPEC2DMAX_YMIN_YMAX_TRINF .or. &
          nspec2D_ymax_trinfinite < 0 .or. nspec2d_ymax_trinfinite > NSPEC2DMAX_YMIN_YMAX_TRINF ) &
        call exit_mpi(myrank, 'Error reading trinfinite boundary')
    endif
    ! infinite
    if (nspec2D_xmin_infinite < 0 .or. nspec2d_xmin_infinite > NSPEC2DMAX_XMIN_XMAX_INF .or. &
        nspec2D_xmax_infinite < 0 .or. nspec2d_xmax_infinite > NSPEC2DMAX_XMIN_XMAX_INF .or. &
        nspec2D_ymin_infinite < 0 .or. nspec2d_ymin_infinite > NSPEC2DMAX_YMIN_YMAX_INF .or. &
        nspec2D_ymax_infinite < 0 .or. nspec2d_ymax_infinite > NSPEC2DMAX_YMIN_YMAX_INF ) &
      call exit_mpi(myrank, 'Error reading infinite boundary')
  endif

  ! Boundary Mesh for crust and mantle
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
    allocate(moho_kl(NGLLX,NGLLY,NSPEC2D_MOHO), &
             moho_kl_top(NGLLX,NGLLY,NSPEC2D_MOHO), &
             moho_kl_bot(NGLLX,NGLLY,NSPEC2D_MOHO), &
             d400_kl(NGLLX,NGLLY,NSPEC2D_400), &
             d400_kl_top(NGLLX,NGLLY,NSPEC2D_400), &
             d400_kl_bot(NGLLX,NGLLY,NSPEC2D_400), &
             d670_kl(NGLLX,NGLLY,NSPEC2D_670), &
             d670_kl_top(NGLLX,NGLLY,NSPEC2D_670), &
             d670_kl_bot(NGLLX,NGLLY,NSPEC2D_670), &
             cmb_kl(NGLLX,NGLLY,NSPEC2D_CMB), &
             cmb_kl_top(NGLLX,NGLLY,NSPEC2D_CMB), &
             cmb_kl_bot(NGLLX,NGLLY,NSPEC2D_CMB), &
             icb_kl(NGLLX,NGLLY,NSPEC2D_ICB), &
             icb_kl_top(NGLLX,NGLLY,NSPEC2D_ICB), &
             icb_kl_bot(NGLLX,NGLLY,NSPEC2D_ICB),stat=ier)
    if (ier /= 0) stop 'Error allocating arrays moho_kl,.. '
    ! initialization
    moho_kl(:,:,:) = 0._CUSTOM_REAL; moho_kl_top(:,:,:) = 0._CUSTOM_REAL; moho_kl_bot(:,:,:) = 0._CUSTOM_REAL
    d400_kl(:,:,:) = 0._CUSTOM_REAL; d400_kl_top(:,:,:) = 0._CUSTOM_REAL; d400_kl_bot(:,:,:) = 0._CUSTOM_REAL
    d670_kl(:,:,:) = 0._CUSTOM_REAL; d670_kl_top(:,:,:) = 0._CUSTOM_REAL; d670_kl_bot(:,:,:) = 0._CUSTOM_REAL
    cmb_kl(:,:,:) = 0._CUSTOM_REAL; cmb_kl_top(:,:,:) = 0._CUSTOM_REAL; cmb_kl_bot(:,:,:) = 0._CUSTOM_REAL
    icb_kl(:,:,:) = 0._CUSTOM_REAL; icb_kl_top(:,:,:) = 0._CUSTOM_REAL; icb_kl_bot(:,:,:) = 0._CUSTOM_REAL
  endif

  end subroutine read_mesh_databases_coupling

!
!-------------------------------------------------------------------------------------------------
!
!
! not used anymore... left here for reference
!
!  subroutine read_mesh_databases_addressing()
!
!  use constants
!
!  use specfem_par
!  use specfem_par_crustmantle
!  use specfem_par_innercore
!  use specfem_par_outercore
!
!  implicit none
!
!  ! local parameters
!  integer, dimension(0:NPROCTOT_VAL-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
!  integer :: ier,iproc,iproc_read,iproc_xi,iproc_eta
!
!  if (I_should_read_the_database) then
!    ! open file with global slice number addressing
!    if (myrank == 0) then
!      open(unit=IIN,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read',iostat=ier)
!      if (ier /= 0 ) call exit_mpi(myrank,'Error opening addressing.txt')
!
!      do iproc = 0,NPROCTOT_VAL-1
!        read(IIN,*) iproc_read,ichunk,iproc_xi,iproc_eta
!
!        if (iproc_read /= iproc) call exit_MPI(myrank,'incorrect slice number read')
!
!        addressing(ichunk,iproc_xi,iproc_eta) = iproc
!        ichunk_slice(iproc) = ichunk
!        iproc_xi_slice(iproc) = iproc_xi
!        iproc_eta_slice(iproc) = iproc_eta
!      enddo
!      close(IIN)
!    endif
!
!    ! broadcast the information read on the main to the nodes
!    call bcast_all_i(addressing,NCHUNKS_VAL*NPROC_XI_VAL*NPROC_ETA_VAL)
!    call bcast_all_i(ichunk_slice,NPROCTOT_VAL)
!    call bcast_all_i(iproc_xi_slice,NPROCTOT_VAL)
!    call bcast_all_i(iproc_eta_slice,NPROCTOT_VAL)
!
!    ! output a topology map of slices - fix 20x by nproc
!    if (myrank == 0) then
!      if (NCHUNKS_VAL == 6 .and. NPROCTOT_VAL < 1000) then
!        write(IMAIN,*) 'Spatial distribution of the slices'
!        do iproc_xi = NPROC_XI_VAL-1, 0, -1
!          write(IMAIN,'(20x)',advance='no')
!          do iproc_eta = NPROC_ETA_VAL -1, 0, -1
!            ichunk = CHUNK_AB
!            write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
!          enddo
!          write(IMAIN,'(1x)',advance='yes')
!        enddo
!        write(IMAIN, *) ' '
!        do iproc_xi = NPROC_XI_VAL-1, 0, -1
!          write(IMAIN,'(1x)',advance='no')
!          do iproc_eta = NPROC_ETA_VAL -1, 0, -1
!            ichunk = CHUNK_BC
!            write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
!          enddo
!          write(IMAIN,'(3x)',advance='no')
!          do iproc_eta = NPROC_ETA_VAL -1, 0, -1
!            ichunk = CHUNK_AC
!            write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
!          enddo
!          write(IMAIN,'(3x)',advance='no')
!          do iproc_eta = NPROC_ETA_VAL -1, 0, -1
!            ichunk = CHUNK_BC_ANTIPODE
!            write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
!          enddo
!          write(IMAIN,'(1x)',advance='yes')
!        enddo
!        write(IMAIN, *) ' '
!        do iproc_xi = NPROC_XI_VAL-1, 0, -1
!          write(IMAIN,'(20x)',advance='no')
!          do iproc_eta = NPROC_ETA_VAL -1, 0, -1
!            ichunk = CHUNK_AB_ANTIPODE
!            write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
!          enddo
!          write(IMAIN,'(1x)',advance='yes')
!        enddo
!        write(IMAIN, *) ' '
!        do iproc_xi = NPROC_XI_VAL-1, 0, -1
!          write(IMAIN,'(20x)',advance='no')
!          do iproc_eta = NPROC_ETA_VAL -1, 0, -1
!            ichunk = CHUNK_AC_ANTIPODE
!            write(IMAIN,'(i5)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
!          enddo
!          write(IMAIN,'(1x)',advance='yes')
!        enddo
!        write(IMAIN, *) ' '
!      endif
!    endif
!
!    ! determine chunk number and local slice coordinates using addressing
!    ! (needed for Stacey conditions)
!    ichunk = ichunk_slice(myrank)
!  endif ! I_should_read_the_database
!
!  call bcast_all_i_for_database(ichunk, 1)
!  call bcast_all_i_for_database(addressing(1,0,0), size(addressing,kind=4))
!
!  end subroutine read_mesh_databases_addressing
!
!
!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_trinfinite
  use specfem_par_infinite
  use specfem_par_full_gravity
  implicit none

  ! local parameters
  real :: percentage_edge
  integer :: ier,i,num_poin,num_elements
  integer :: iglob_min,iglob_max
  integer :: ispec_min,ispec_max
  integer :: shape_2d(2), shape_3d(3)

  ! read MPI interfaces from file
  num_interfaces_crust_mantle = 0
  num_interfaces_outer_core = 0
  num_interfaces_inner_core = 0
  num_interfaces_trinfinite = 0
  num_interfaces_infinite = 0

  max_nibool_interfaces_cm = 0
  max_nibool_interfaces_oc = 0
  max_nibool_interfaces_ic = 0
  max_nibool_interfaces_trinfinite = 0
  max_nibool_interfaces_infinite = 0

  !
  ! crust mantle
  !
  if (I_should_read_the_database .and. NSPEC_CRUST_MANTLE > 0) then
    if (ADIOS_FOR_MPI_ARRAYS) then
      call read_mesh_databases_MPI_adios(IREGION_CRUST_MANTLE)
    else
      call read_mesh_databases_MPI_CM()
    endif
  endif
  call bcast_mesh_databases_MPI_CM()

  ! checks interface values read
  if (NSPEC_CRUST_MANTLE > 0) then
    if (minval(my_neighbors_crust_mantle) < 0 .or. maxval(my_neighbors_crust_mantle) >= NPROCTOT) then
      print *,'Error: invalid MPI neighbors min/max',minval(my_neighbors_crust_mantle),maxval(my_neighbors_crust_mantle),NPROCTOT
      call exit_mpi(myrank,'Error invalid MPI neighbors crust_mantle')
    endif
  endif

  do i = 1,num_interfaces_crust_mantle
    ! number of points on interface
    num_poin = nibool_interfaces_crust_mantle(i)
    if (num_poin <= 0 .or. num_poin > NGLOB_CRUST_MANTLE) then
      print *,'Error: invalid nibool_interfaces_crust_mantle ',num_poin,'interface',i,'nglob',NGLOB_CRUST_MANTLE
      call exit_mpi(myrank,'Error invalid nibool_interfaces_crust_mantle')
    endif
    ! iglob min/max
    iglob_min = minval(ibool_interfaces_crust_mantle(1:num_poin,i))
    iglob_max = maxval(ibool_interfaces_crust_mantle(1:num_poin,i))
    if (iglob_min <= 0 .or. iglob_max > NGLOB_CRUST_MANTLE) then
      print *,'Error: invalid ibool_interfaces_crust_mantle min/max ',iglob_min,iglob_max,'interface',i,'nglob',NGLOB_CRUST_MANTLE
      call exit_mpi(myrank,'Error invalid ibool_interfaces_crust_mantle')
    endif
  enddo
  ! checks min/max
  if (num_phase_ispec_crust_mantle > 0) then
    do i = 1,2
      if (i == 1) then
        num_elements = nspec_outer_crust_mantle
      else
        num_elements = nspec_inner_crust_mantle
      endif
      ispec_min = minval(phase_ispec_inner_crust_mantle(1:num_elements,i))
      ispec_max = maxval(phase_ispec_inner_crust_mantle(1:num_elements,i))
      if (ispec_min <= 0 .or. ispec_max > NSPEC_CRUST_MANTLE) then
        print *,'Error: invalid phase_ispec_inner_crust_mantle min/max ',ispec_min,ispec_max,'phase',i, &
                nspec_outer_crust_mantle,nspec_inner_crust_mantle,'nspec',NSPEC_CRUST_MANTLE
        call exit_mpi(myrank,'Error invalid phase_ispec_inner_crust_mantle')
      endif
    enddo
  endif

  ! MPI buffers
  if (USE_CUDA_AWARE_MPI) then
    ! allocates buffers on GPU
    shape_3d(1) = NDIM
    shape_3d(2) = max_nibool_interfaces_cm
    shape_3d(3) = num_interfaces_crust_mantle
    call allocate_gpu_buffer_3d(buffer_send_vector_crust_mantle,shape_3d)
    call allocate_gpu_buffer_3d(buffer_recv_vector_crust_mantle,shape_3d)
  else
    ! allocates buffers on CPU
    allocate(buffer_send_vector_crust_mantle(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
             buffer_recv_vector_crust_mantle(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_vector_crust_mantle etc.')
    buffer_send_vector_crust_mantle(:,:,:) = 0.0; buffer_recv_vector_crust_mantle(:,:,:) = 0.0
  endif
  ! request buffers
  allocate(request_send_vector_cm(num_interfaces_crust_mantle), &
           request_recv_vector_cm(num_interfaces_crust_mantle),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array request_send_vector_cm etc.')
  request_send_vector_cm(:) = 0; request_recv_vector_cm(:) = 0

  if (SIMULATION_TYPE == 3) then
    ! MPI buffers
    if (USE_CUDA_AWARE_MPI) then
      ! allocates buffers on GPU
      shape_3d(1) = NDIM
      shape_3d(2) = max_nibool_interfaces_cm
      shape_3d(3) = num_interfaces_crust_mantle
      call allocate_gpu_buffer_3d(b_buffer_send_vector_cm,shape_3d)
      call allocate_gpu_buffer_3d(b_buffer_recv_vector_cm,shape_3d)
    else
      ! allocates buffers on CPU
      allocate(b_buffer_send_vector_cm(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
               b_buffer_recv_vector_cm(NDIM,max_nibool_interfaces_cm,num_interfaces_crust_mantle),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_vector_cm etc.')
      b_buffer_send_vector_cm(:,:,:) = 0.0; b_buffer_recv_vector_cm(:,:,:) = 0.0
    endif
    ! request buffers
    allocate(b_request_send_vector_cm(num_interfaces_crust_mantle), &
             b_request_recv_vector_cm(num_interfaces_crust_mantle),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_request_send_vector_cm etc.')
    b_request_send_vector_cm(:) = 0; b_request_recv_vector_cm(:) = 0
  endif

  !
  ! outer core
  !
  if (I_should_read_the_database .and. NSPEC_OUTER_CORE > 0) then
    if (ADIOS_FOR_MPI_ARRAYS) then
      call read_mesh_databases_MPI_adios(IREGION_OUTER_CORE)
    else
      call read_mesh_databases_MPI_OC()
    endif
  endif
  call bcast_mesh_databases_MPI_OC()

  ! checks interface values read
  if (NSPEC_OUTER_CORE > 0) then
    if (minval(my_neighbors_outer_core) < 0 .or. maxval(my_neighbors_outer_core) >= NPROCTOT) then
      print *,'Error: invalid MPI neighbors min/max',minval(my_neighbors_outer_core),maxval(my_neighbors_outer_core),NPROCTOT
      call exit_mpi(myrank,'Error invalid MPI neighbors outer_core')
    endif
  endif

  do i = 1,num_interfaces_outer_core
    ! number of points on interface
    num_poin = nibool_interfaces_outer_core(i)
    if (num_poin <= 0 .or. num_poin > NGLOB_outer_core) then
      print *,'Error: invalid nibool_interfaces_outer_core ',num_poin,'interface',i,'nglob',NGLOB_outer_core
      call exit_mpi(myrank,'Error invalid nibool_interfaces_outer_core')
    endif
    ! iglob min/max
    iglob_min = minval(ibool_interfaces_outer_core(1:num_poin,i))
    iglob_max = maxval(ibool_interfaces_outer_core(1:num_poin,i))
    if (iglob_min <= 0 .or. iglob_max > NGLOB_outer_core) then
      print *,'Error: invalid ibool_interfaces_outer_core min/max ',iglob_min,iglob_max,'interface',i,'nglob',NGLOB_outer_core
      call exit_mpi(myrank,'Error invalid ibool_interfaces_outer_core')
    endif
  enddo
  ! checks min/max
  if (num_phase_ispec_outer_core > 0) then
    do i = 1,2
      if (i == 1) then
        num_elements = nspec_outer_outer_core
      else
        num_elements = nspec_inner_outer_core
      endif
      ispec_min = minval(phase_ispec_inner_outer_core(1:num_elements,i))
      ispec_max = maxval(phase_ispec_inner_outer_core(1:num_elements,i))
      if (ispec_min <= 0 .or. ispec_max > NSPEC_OUTER_CORE) then
        print *,'Error: invalid phase_ispec_inner_outer_core min/max ',ispec_min,ispec_max,'phase',i, &
                nspec_outer_outer_core,nspec_inner_outer_core,'nspec',NSPEC_OUTER_CORE
        call exit_mpi(myrank,'Error invalid phase_ispec_inner_outer_core')
      endif
    enddo
  endif

  ! MPI buffers
  if (USE_CUDA_AWARE_MPI) then
    ! allocates buffers on GPU
    shape_2d(1) = max_nibool_interfaces_oc
    shape_2d(2) = num_interfaces_outer_core
    call allocate_gpu_buffer_2d(buffer_send_scalar_outer_core,shape_2d)
    call allocate_gpu_buffer_2d(buffer_recv_scalar_outer_core,shape_2d)
  else
    ! allocates buffers on CPU
    allocate(buffer_send_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
             buffer_recv_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_vector_outer_core etc.')
    buffer_send_scalar_outer_core(:,:) = 0.0; buffer_recv_scalar_outer_core(:,:) = 0.0
  endif
  ! request buffers
  allocate(request_send_scalar_oc(num_interfaces_outer_core), &
           request_recv_scalar_oc(num_interfaces_outer_core),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array request_send_scalar_oc etc.')
  request_send_scalar_oc(:) = 0; request_recv_scalar_oc(:) = 0

  if (SIMULATION_TYPE == 3) then
    if (USE_CUDA_AWARE_MPI) then
      ! allocates buffers on GPU
      shape_2d(1) = max_nibool_interfaces_oc
      shape_2d(2) = num_interfaces_outer_core
      call allocate_gpu_buffer_2d(b_buffer_send_scalar_outer_core,shape_2d)
      call allocate_gpu_buffer_2d(b_buffer_recv_scalar_outer_core,shape_2d)
    else
      ! allocates buffers on CPU
      allocate(b_buffer_send_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
               b_buffer_recv_scalar_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_vector_outer_core etc.')
      b_buffer_send_scalar_outer_core(:,:) = 0.0; b_buffer_recv_scalar_outer_core(:,:) = 0.0
    endif
    ! request buffers
    allocate(b_request_send_scalar_oc(num_interfaces_outer_core), &
             b_request_recv_scalar_oc(num_interfaces_outer_core),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_request_send_scalar_oc etc.')
    b_request_send_scalar_oc(:) = 0; b_request_recv_scalar_oc(:) = 0
  endif

  !
  ! inner core
  !
  if (I_should_read_the_database .and. NSPEC_INNER_CORE > 0) then
    if (ADIOS_FOR_MPI_ARRAYS) then
      call read_mesh_databases_MPI_adios(IREGION_INNER_CORE)
    else
      call read_mesh_databases_MPI_IC()
    endif
  endif
  call bcast_mesh_databases_MPI_IC()

  ! checks interface values read
  if (NSPEC_INNER_CORE > 0) then
    if (minval(my_neighbors_inner_core) < 0 .or. maxval(my_neighbors_inner_core) >= NPROCTOT) then
      print *,'Error: invalid MPI neighbors min/max',minval(my_neighbors_inner_core),maxval(my_neighbors_inner_core),NPROCTOT
      call exit_mpi(myrank,'Error invalid MPI neighbors inner_core')
    endif
  endif

  do i = 1,num_interfaces_inner_core
    ! number of points on interface
    num_poin = nibool_interfaces_inner_core(i)
    if (num_poin <= 0 .or. num_poin > NGLOB_inner_core) then
      print *,'Error: invalid nibool_interfaces_inner_core ',num_poin,'interface',i,'nglob',NGLOB_inner_core
      call exit_mpi(myrank,'Error invalid nibool_interfaces_inner_core')
    endif
    ! iglob min/max
    iglob_min = minval(ibool_interfaces_inner_core(1:num_poin,i))
    iglob_max = maxval(ibool_interfaces_inner_core(1:num_poin,i))
    if (iglob_min <= 0 .or. iglob_max > NGLOB_inner_core) then
      print *,'Error: invalid ibool_interfaces_inner_core min/max ',iglob_min,iglob_max,'interface',i,'nglob',NGLOB_inner_core
      call exit_mpi(myrank,'Error invalid ibool_interfaces_inner_core')
    endif
  enddo
  ! checks min/max
  if (num_phase_ispec_inner_core > 0) then
    do i = 1,2
      if (i == 1) then
        num_elements = nspec_outer_inner_core
      else
        num_elements = nspec_inner_inner_core
      endif
      ispec_min = minval(phase_ispec_inner_inner_core(1:num_elements,i))
      ispec_max = maxval(phase_ispec_inner_inner_core(1:num_elements,i))
      if (ispec_min <= 0 .or. ispec_max > NSPEC_INNER_CORE) then
        print *,'Error: invalid phase_ispec_inner_inner_core min/max ',ispec_min,ispec_max,'phase',i, &
                nspec_outer_inner_core,nspec_inner_inner_core,'nspec',NSPEC_INNER_CORE
        call exit_mpi(myrank,'Error invalid phase_ispec_inner_inner_core')
      endif
    enddo
  endif

  ! MPI buffers
  if (USE_CUDA_AWARE_MPI) then
    ! allocates buffers on GPU
    ! allocates buffers on GPU
    shape_3d(1) = NDIM
    shape_3d(2) = max_nibool_interfaces_ic
    shape_3d(3) = num_interfaces_inner_core
    call allocate_gpu_buffer_3d(buffer_send_vector_inner_core,shape_3d)
    call allocate_gpu_buffer_3d(buffer_recv_vector_inner_core,shape_3d)
  else
    ! allocates buffers on CPU
    allocate(buffer_send_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core), &
             buffer_recv_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_vector_inner_core etc.')
    buffer_send_vector_inner_core(:,:,:) = 0.0; buffer_recv_vector_inner_core(:,:,:) = 0.0
  endif
  ! request buffers
  allocate(request_send_vector_ic(num_interfaces_inner_core), &
           request_recv_vector_ic(num_interfaces_inner_core),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array request_send_vector_ic etc.')
  request_send_vector_ic(:) = 0; request_recv_vector_ic(:) = 0

  if (SIMULATION_TYPE == 3) then
    if (USE_CUDA_AWARE_MPI) then
      ! allocates buffers on GPU
      shape_3d(1) = NDIM
      shape_3d(2) = max_nibool_interfaces_ic
      shape_3d(3) = num_interfaces_inner_core
      call allocate_gpu_buffer_3d(b_buffer_send_vector_inner_core,shape_3d)
      call allocate_gpu_buffer_3d(b_buffer_recv_vector_inner_core,shape_3d)
    else
      ! allocates buffers on CPU
      allocate(b_buffer_send_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core), &
               b_buffer_recv_vector_inner_core(NDIM,max_nibool_interfaces_ic,num_interfaces_inner_core),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_vector_inner_core etc.')
      b_buffer_send_vector_inner_core(:,:,:) = 0.0; b_buffer_recv_vector_inner_core(:,:,:) = 0.0
    endif
    ! request buffers
    allocate(b_request_send_vector_ic(num_interfaces_inner_core), &
             b_request_recv_vector_ic(num_interfaces_inner_core),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_request_send_vector_ic etc.')
    b_request_send_vector_ic(:) = 0; b_request_recv_vector_ic(:) = 0
  endif

  ! full gravity
  if (FULL_GRAVITY) then
    if (ADD_TRINF) then
      !
      ! transition-to-infinite
      !
      if (I_should_read_the_database .and. NSPEC_TRINFINITE > 0) then
        if (ADIOS_FOR_MPI_ARRAYS) then
          call read_mesh_databases_MPI_adios(IREGION_TRINFINITE)
        else
          call read_mesh_databases_MPI_TRINF()
        endif
      endif
      call bcast_mesh_databases_MPI_TRINF()

      ! checks interface values read
      if (NSPEC_TRINFINITE > 0) then
        if (minval(my_neighbors_trinfinite) < 0 .or. maxval(my_neighbors_trinfinite) >= NPROCTOT) then
          print *,'Error: invalid MPI neighbors min/max',minval(my_neighbors_trinfinite),maxval(my_neighbors_trinfinite),NPROCTOT
          call exit_mpi(myrank,'Error invalid MPI neighbors trinfinite')
        endif
      endif

      do i = 1,num_interfaces_trinfinite
        ! number of points on interface
        num_poin = nibool_interfaces_trinfinite(i)
        if (num_poin <= 0 .or. num_poin > NGLOB_TRINFINITE) then
          print *,'Error: invalid nibool_interfaces_trinfinite ',num_poin,'interface',i,'nglob',NGLOB_TRINFINITE
          call exit_mpi(myrank,'Error invalid nibool_interfaces_trinfinite')
        endif
        ! iglob min/max
        iglob_min = minval(ibool_interfaces_trinfinite(1:num_poin,i))
        iglob_max = maxval(ibool_interfaces_trinfinite(1:num_poin,i))
        if (iglob_min <= 0 .or. iglob_max > NGLOB_TRINFINITE) then
          print *,'Error: invalid ibool_interfaces_trinfinite min/max ',iglob_min,iglob_max,'interface',i,'nglob',NGLOB_TRINFINITE
          call exit_mpi(myrank,'Error invalid ibool_interfaces_trinfinite')
        endif
      enddo
      ! checks min/max
      if (num_phase_ispec_trinfinite > 0) then
        do i = 1,2
          if (i == 1) then
            num_elements = nspec_outer_trinfinite
          else
            num_elements = nspec_inner_trinfinite
          endif
          ispec_min = minval(phase_ispec_inner_trinfinite(1:num_elements,i))
          ispec_max = maxval(phase_ispec_inner_trinfinite(1:num_elements,i))
          if (ispec_min <= 0 .or. ispec_max > NSPEC_TRINFINITE) then
            print *,'Error: invalid phase_ispec_inner_trinfinite min/max ',ispec_min,ispec_max,'phase',i, &
                    nspec_outer_trinfinite,nspec_inner_trinfinite,'nspec',NSPEC_TRINFINITE
            call exit_mpi(myrank,'Error invalid phase_ispec_inner_trinfinite')
          endif
        enddo
      endif

      ! MPI buffers
      if (USE_CUDA_AWARE_MPI) then
        ! allocates buffers on GPU
        shape_2d(1) = max_nibool_interfaces_trinfinite
        shape_2d(2) = num_interfaces_trinfinite
        call allocate_gpu_buffer_2d(buffer_send_scalar_trinfinite,shape_2d)
        call allocate_gpu_buffer_2d(buffer_recv_scalar_trinfinite,shape_2d)
      else
        ! allocates buffers on CPU
        allocate(buffer_send_scalar_trinfinite(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite), &
                 buffer_recv_scalar_trinfinite(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite),stat=ier)
        if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_scalar_trinfinite etc.')
        buffer_send_scalar_trinfinite(:,:) = 0.0; buffer_recv_scalar_trinfinite(:,:) = 0.0
      endif
      ! request buffers
      allocate(request_send_scalar_trinfinite(num_interfaces_trinfinite), &
               request_recv_scalar_trinfinite(num_interfaces_trinfinite),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array request_send_scalar_trinf etc.')
      request_send_scalar_trinfinite(:) = 0; request_recv_scalar_trinfinite(:) = 0

      if (SIMULATION_TYPE == 3) then
        if (USE_CUDA_AWARE_MPI) then
          ! allocates buffers on GPU
          shape_2d(1) = max_nibool_interfaces_trinfinite
          shape_2d(2) = num_interfaces_trinfinite
          call allocate_gpu_buffer_2d(b_buffer_send_scalar_trinfinite,shape_2d)
          call allocate_gpu_buffer_2d(b_buffer_recv_scalar_trinfinite,shape_2d)
        else
          ! allocates buffers on CPU
          allocate(b_buffer_send_scalar_trinfinite(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite), &
                   b_buffer_recv_scalar_trinfinite(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite),stat=ier)
          if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_scalar_trinfinite etc.')
          b_buffer_send_scalar_trinfinite(:,:) = 0.0; b_buffer_recv_scalar_trinfinite(:,:) = 0.0
        endif
        ! request buffers
        allocate(b_request_send_scalar_trinfinite(num_interfaces_trinfinite), &
                 b_request_recv_scalar_trinfinite(num_interfaces_trinfinite),stat=ier)
        if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_request_send_scalar_trinfinite etc.')
        b_request_send_scalar_trinfinite(:) = 0; b_request_recv_scalar_trinfinite(:) = 0
      endif
    endif

    !
    ! infinite
    !
    if (I_should_read_the_database .and. NSPEC_INFINITE > 0) then
      if (ADIOS_FOR_MPI_ARRAYS) then
        call read_mesh_databases_MPI_adios(IREGION_INFINITE)
      else
        call read_mesh_databases_MPI_INF()
      endif
    endif
    call bcast_mesh_databases_MPI_INF()

    ! checks interface values read
    if (NSPEC_INFINITE > 0) then
      if (minval(my_neighbors_infinite) < 0 .or. maxval(my_neighbors_infinite) >= NPROCTOT) then
        print *,'Error: invalid MPI neighbors min/max',minval(my_neighbors_infinite),maxval(my_neighbors_infinite),NPROCTOT
        call exit_mpi(myrank,'Error invalid MPI neighbors infinite')
      endif
    endif

    do i = 1,num_interfaces_infinite
      ! number of points on interface
      num_poin = nibool_interfaces_infinite(i)
      if (num_poin <= 0 .or. num_poin > NGLOB_INFINITE) then
        print *,'Error: invalid nibool_interfaces_infinite ',num_poin,'interface',i,'nglob',NGLOB_INFINITE
        call exit_mpi(myrank,'Error invalid nibool_interfaces_infinite')
      endif
      ! iglob min/max
      iglob_min = minval(ibool_interfaces_infinite(1:num_poin,i))
      iglob_max = maxval(ibool_interfaces_infinite(1:num_poin,i))
      if (iglob_min <= 0 .or. iglob_max > NGLOB_INFINITE) then
        print *,'Error: invalid ibool_interfaces_infinite min/max ',iglob_min,iglob_max,'interface',i,'nglob',NGLOB_INFINITE
        call exit_mpi(myrank,'Error invalid ibool_interfaces_infinite')
      endif
    enddo
    ! checks min/max
    if (num_phase_ispec_infinite > 0) then
      do i = 1,2
        if (i == 1) then
          num_elements = nspec_outer_infinite
        else
          num_elements = nspec_inner_infinite
        endif
        ispec_min = minval(phase_ispec_inner_infinite(1:num_elements,i))
        ispec_max = maxval(phase_ispec_inner_infinite(1:num_elements,i))
        if (ispec_min <= 0 .or. ispec_max > NSPEC_INFINITE) then
          print *,'Error: invalid phase_ispec_inner_infinite min/max ',ispec_min,ispec_max,'phase',i, &
                  nspec_outer_infinite,nspec_inner_infinite,'nspec',NSPEC_INFINITE
          call exit_mpi(myrank,'Error invalid phase_ispec_inner_infinite')
        endif
      enddo
    endif

    ! MPI buffers
    if (USE_CUDA_AWARE_MPI) then
      ! allocates buffers on GPU
      shape_2d(1) = max_nibool_interfaces_infinite
      shape_2d(2) = num_interfaces_infinite
      call allocate_gpu_buffer_2d(buffer_send_scalar_infinite,shape_2d)
      call allocate_gpu_buffer_2d(buffer_recv_scalar_infinite,shape_2d)
    else
      ! allocates buffers on CPU
      allocate(buffer_send_scalar_infinite(max_nibool_interfaces_infinite,num_interfaces_infinite), &
               buffer_recv_scalar_infinite(max_nibool_interfaces_infinite,num_interfaces_infinite),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array buffer_send_scalar_infinite etc.')
      buffer_send_scalar_infinite(:,:) = 0.0; buffer_recv_scalar_infinite(:,:) = 0.0
    endif
    ! request buffers
    allocate(request_send_scalar_infinite(num_interfaces_infinite), &
             request_recv_scalar_infinite(num_interfaces_infinite),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array request_send_scalar_infinite etc.')
    request_send_scalar_infinite(:) = 0; request_recv_scalar_infinite(:) = 0

    if (SIMULATION_TYPE == 3) then
      if (USE_CUDA_AWARE_MPI) then
        ! allocates buffers on GPU
        shape_2d(1) = max_nibool_interfaces_infinite
        shape_2d(2) = num_interfaces_infinite
        call allocate_gpu_buffer_2d(b_buffer_send_scalar_infinite,shape_2d)
        call allocate_gpu_buffer_2d(b_buffer_recv_scalar_infinite,shape_2d)
      else
        ! allocates buffers on CPU
        allocate(b_buffer_send_scalar_infinite(max_nibool_interfaces_infinite,num_interfaces_infinite), &
                 b_buffer_recv_scalar_infinite(max_nibool_interfaces_infinite,num_interfaces_infinite),stat=ier)
        if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_buffer_send_scalar_infinite etc.')
        b_buffer_send_scalar_infinite(:,:) = 0.0; b_buffer_recv_scalar_infinite(:,:) = 0.0
      endif
      ! request buffers
      allocate(b_request_send_scalar_infinite(num_interfaces_infinite), &
               b_request_recv_scalar_infinite(num_interfaces_infinite),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array b_request_send_scalar_infinite etc.')
      b_request_send_scalar_infinite(:) = 0; b_request_recv_scalar_infinite(:) = 0
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  for overlapping of communications with calculations:'
    write(IMAIN,*)
    if (NSPEC_CRUST_MANTLE > 0) then
      percentage_edge = 100. * nspec_outer_crust_mantle / real(NSPEC_CRUST_MANTLE)
      write(IMAIN,*) '  percentage of edge elements in crust/mantle ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in crust/mantle ',100. - percentage_edge,'%'
      write(IMAIN,*)
    endif
    if (NSPEC_OUTER_CORE > 0) then
      percentage_edge = 100.* nspec_outer_outer_core / real(NSPEC_OUTER_CORE)
      write(IMAIN,*) '  percentage of edge elements in outer core ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in outer core ',100. - percentage_edge,'%'
      write(IMAIN,*)
    endif
    if (NSPEC_INNER_CORE > 0) then
      percentage_edge = 100. * nspec_outer_inner_core / real(NSPEC_INNER_CORE)
      write(IMAIN,*) '  percentage of edge elements in inner core ',percentage_edge,'%'
      write(IMAIN,*) '  percentage of volume elements in inner core ',100. - percentage_edge,'%'
      write(IMAIN,*)
    endif

    if (FULL_GRAVITY) then
      if (ADD_TRINF .and. NSPEC_TRINFINITE > 0) then
        percentage_edge = 100. * nspec_outer_trinfinite / real(NSPEC_TRINFINITE)
        write(IMAIN,*) '  percentage of edge elements in trinfinite layer ',percentage_edge,'%'
        write(IMAIN,*) '  percentage of volume elements in trinfinite layer ',100. - percentage_edge,'%'
        write(IMAIN,*)
      endif
      if (NSPEC_INFINITE > 0) then
        percentage_edge = 100. * nspec_outer_infinite / real(NSPEC_INFINITE)
        write(IMAIN,*) '  percentage of edge elements in infinite layer ',percentage_edge,'%'
        write(IMAIN,*) '  percentage of volume elements in infinite layer ',100. - percentage_edge,'%'
        write(IMAIN,*)
      endif
    endif

    call flush_IMAIN()
  endif
  ! synchronizes MPI processes
  call synchronize_all()

  contains

    subroutine allocate_gpu_buffer_2d(buffer, shape)

    use iso_c_binding, only: c_ptr,c_f_pointer

    implicit none
    real(kind=CUSTOM_REAL), dimension(:,:), pointer, intent(inout) :: buffer
    integer, intent(in) :: shape(2)
    ! local parameters
    type(c_ptr) :: c_ptr_buffer

    call allocate_gpu_buffer(c_ptr_buffer, shape(1) * shape(2))
    call c_f_pointer(c_ptr_buffer, buffer, shape)

    end subroutine allocate_gpu_buffer_2d

    !--------

    subroutine allocate_gpu_buffer_3d(buffer, shape)

    use iso_c_binding, only: c_ptr,c_f_pointer

    implicit none
    real(kind=CUSTOM_REAL), dimension(:,:,:), pointer, intent(inout) :: buffer
    integer, intent(in) :: shape(3)
    ! local parameters
    type(c_ptr) :: c_ptr_buffer

    call allocate_gpu_buffer(c_ptr_buffer, shape(1) * shape(2) * shape(3))
    call c_f_pointer(c_ptr_buffer, buffer, shape)

    end subroutine allocate_gpu_buffer_3d

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
  allocate(my_neighbors_crust_mantle(num_interfaces_crust_mantle), &
           nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_crust_mantle etc.')
  my_neighbors_crust_mantle(:) = 0; nibool_interfaces_crust_mantle(:) = 0

  if (num_interfaces_crust_mantle > 0) then
    read(IIN) max_nibool_interfaces_cm
    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
    ibool_interfaces_crust_mantle(:,:) = 0

    read(IIN) my_neighbors_crust_mantle
    read(IIN) nibool_interfaces_crust_mantle
    read(IIN) ibool_interfaces_crust_mantle
  else
    ! dummy array
    max_nibool_interfaces_cm = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_crust_mantle')
    ibool_interfaces_crust_mantle(:,:) = 0
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_crust_mantle,nspec_outer_crust_mantle
  read(IIN) num_phase_ispec_crust_mantle
  if (num_phase_ispec_crust_mantle < 0 ) call exit_mpi(myrank,'Error num_phase_ispec_crust_mantle is < zero')

  allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')
  phase_ispec_inner_crust_mantle(:,:) = 0

  if (num_phase_ispec_crust_mantle > 0 ) read(IIN) phase_ispec_inner_crust_mantle

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_crust_mantle,num_colors_inner_crust_mantle

    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
    num_elem_colors_crust_mantle(:) = 0

    read(IIN) num_elem_colors_crust_mantle
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
    num_elem_colors_crust_mantle(:) = 0
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
  allocate(my_neighbors_outer_core(num_interfaces_outer_core), &
           nibool_interfaces_outer_core(num_interfaces_outer_core), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_outer_core etc.')
  my_neighbors_outer_core(:) = 0; nibool_interfaces_outer_core(:) = 0

  if (num_interfaces_outer_core > 0) then
    read(IIN) max_nibool_interfaces_oc
    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_outer_core')
    ibool_interfaces_outer_core(:,:) = 0

    read(IIN) my_neighbors_outer_core
    read(IIN) nibool_interfaces_outer_core
    read(IIN) ibool_interfaces_outer_core
  else
    ! dummy array
    max_nibool_interfaces_oc = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_outer_core')
    ibool_interfaces_outer_core(:,:) = 0
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_outer_core,nspec_outer_outer_core
  read(IIN) num_phase_ispec_outer_core
  if (num_phase_ispec_outer_core < 0 ) call exit_mpi(myrank,'Error num_phase_ispec_outer_core is < zero')

  allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')
  phase_ispec_inner_outer_core(:,:) = 0

  if (num_phase_ispec_outer_core > 0 ) read(IIN) phase_ispec_inner_outer_core

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_outer_core,num_colors_inner_outer_core

    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
    num_elem_colors_outer_core(:) = 0

    read(IIN) num_elem_colors_outer_core
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
    num_elem_colors_outer_core(:) = 0
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
  allocate(my_neighbors_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &
          stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_inner_core etc.')
  my_neighbors_inner_core(:) = 0; nibool_interfaces_inner_core(:) = 0

  if (num_interfaces_inner_core > 0) then
    read(IIN) max_nibool_interfaces_ic
    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic,num_interfaces_inner_core), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_inner_core')
    ibool_interfaces_inner_core(:,:) = 0

    read(IIN) my_neighbors_inner_core
    read(IIN) nibool_interfaces_inner_core
    read(IIN) ibool_interfaces_inner_core
  else
    ! dummy array
    max_nibool_interfaces_ic = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_inner_core')
    ibool_interfaces_inner_core(:,:) = 0
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_inner_core,nspec_outer_inner_core
  read(IIN) num_phase_ispec_inner_core
  if (num_phase_ispec_inner_core < 0 ) call exit_mpi(myrank,'Error num_phase_ispec_inner_core is < zero')

  allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')
  phase_ispec_inner_inner_core(:,:) = 0

  if (num_phase_ispec_inner_core > 0 ) read(IIN) phase_ispec_inner_inner_core

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_inner_core,num_colors_inner_inner_core

    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
    num_elem_colors_inner_core(:) = 0

    read(IIN) num_elem_colors_inner_core
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
    num_elem_colors_inner_core(:) = 0
  endif

  close(IIN)

  end subroutine read_mesh_databases_MPI_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI_TRINF()

  use specfem_par
  use specfem_par_trinfinite
  use specfem_par_full_gravity
  implicit none

  ! local parameters
  integer :: ier

  ! transition-to-infinite region

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,IREGION_TRINFINITE,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_trinfinite
  allocate(my_neighbors_trinfinite(num_interfaces_trinfinite), &
           nibool_interfaces_trinfinite(num_interfaces_trinfinite), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_trinfinite etc.')
  my_neighbors_trinfinite(:) = 0; nibool_interfaces_trinfinite(:) = 0

  if (num_interfaces_trinfinite > 0) then
    read(IIN) max_nibool_interfaces_trinfinite
    allocate(ibool_interfaces_trinfinite(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_trinfinite')
    ibool_interfaces_trinfinite(:,:) = 0

    read(IIN) my_neighbors_trinfinite
    read(IIN) nibool_interfaces_trinfinite
    read(IIN) ibool_interfaces_trinfinite
  else
    ! dummy array
    max_nibool_interfaces_trinfinite = 0
    allocate(ibool_interfaces_trinfinite(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_trinfinite')
    ibool_interfaces_trinfinite(:,:) = 0
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_trinfinite,nspec_outer_trinfinite
  read(IIN) num_phase_ispec_trinfinite
  if (num_phase_ispec_trinfinite < 0 ) call exit_mpi(myrank,'Error num_phase_ispec_trinfinite is < zero')

  allocate(phase_ispec_inner_trinfinite(num_phase_ispec_trinfinite,2),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_trinfinite')
  phase_ispec_inner_trinfinite(:,:) = 0

  if (num_phase_ispec_trinfinite > 0 ) read(IIN) phase_ispec_inner_trinfinite

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_trinfinite,num_colors_inner_trinfinite

    allocate(num_elem_colors_trinfinite(num_colors_outer_trinfinite + num_colors_inner_trinfinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_trinfinite array')
    num_elem_colors_trinfinite(:) = 0

    read(IIN) num_elem_colors_trinfinite
  else
    ! allocates dummy arrays
    num_colors_outer_trinfinite = 0
    num_colors_inner_trinfinite = 0
    allocate(num_elem_colors_trinfinite(num_colors_outer_trinfinite + num_colors_inner_trinfinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_trinfinite array')
    num_elem_colors_trinfinite(:) = 0
  endif

  close(IIN)

  end subroutine read_mesh_databases_MPI_TRINF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_MPI_INF()

  use specfem_par
  use specfem_par_infinite
  use specfem_par_full_gravity
  implicit none

  ! local parameters
  integer :: ier

  ! infinite region

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,IREGION_INFINITE,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_mpi.bin', &
       status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data_mpi.bin')

  ! MPI interfaces
  read(IIN) num_interfaces_infinite
  allocate(my_neighbors_infinite(num_interfaces_infinite), &
           nibool_interfaces_infinite(num_interfaces_infinite), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_infinite etc.')
  my_neighbors_infinite(:) = 0; nibool_interfaces_infinite(:) = 0

  if (num_interfaces_infinite > 0) then
    read(IIN) max_nibool_interfaces_infinite
    allocate(ibool_interfaces_infinite(max_nibool_interfaces_infinite,num_interfaces_infinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_infinite')
    ibool_interfaces_infinite(:,:) = 0

    read(IIN) my_neighbors_infinite
    read(IIN) nibool_interfaces_infinite
    read(IIN) ibool_interfaces_infinite
  else
    ! dummy array
    max_nibool_interfaces_infinite = 0
    allocate(ibool_interfaces_infinite(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_infinite')
    ibool_interfaces_infinite(:,:) = 0
  endif

  ! inner / outer elements
  read(IIN) nspec_inner_infinite,nspec_outer_infinite
  read(IIN) num_phase_ispec_infinite
  if (num_phase_ispec_infinite < 0 ) call exit_mpi(myrank,'Error num_phase_ispec_infinite is < zero')

  allocate(phase_ispec_inner_infinite(num_phase_ispec_infinite,2),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_infinite')
  phase_ispec_inner_infinite(:,:) = 0

  if (num_phase_ispec_infinite > 0 ) read(IIN) phase_ispec_inner_infinite

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    read(IIN) num_colors_outer_infinite,num_colors_inner_infinite

    allocate(num_elem_colors_infinite(num_colors_outer_infinite + num_colors_inner_infinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_infinite array')
    num_elem_colors_infinite(:) = 0

    read(IIN) num_elem_colors_infinite
  else
    ! allocates dummy arrays
    num_colors_outer_infinite = 0
    num_colors_inner_infinite = 0
    allocate(num_elem_colors_infinite(num_colors_outer_infinite + num_colors_inner_infinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_infinite array')
    num_elem_colors_infinite(:) = 0
  endif

  close(IIN)

  end subroutine read_mesh_databases_MPI_INF

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

  ! initializes
  num_abs_boundary_faces_crust_mantle = 0
  num_abs_boundary_faces_outer_core = 0

  ! reads in arrays
  if (I_should_read_the_database) then
    if (ADIOS_FOR_ARRAYS_SOLVER) then
      call read_mesh_databases_stacey_adios()
    else
      ! crust and mantle
      if (NSPEC_CRUST_MANTLE > 0) then
        ! create name of database
        call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

        ! read arrays for Stacey conditions
        open(unit=IIN,file=prname(1:len_trim(prname))//'stacey.bin', &
              status='old',form='unformatted',action='read',iostat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey.bin file for crust mantle')

        read(IIN) num_abs_boundary_faces_crust_mantle

        ! allocates absorbing boundary arrays
        if (num_abs_boundary_faces_crust_mantle > 0) then
          allocate(abs_boundary_ispec_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
          allocate(abs_boundary_ijk_crust_mantle(3,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
          allocate(abs_boundary_jacobian2Dw_crust_mantle(NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
          allocate(abs_boundary_normal_crust_mantle(NDIM,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
          allocate(abs_boundary_npoin_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
          if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

          ! reads in arrays
          read(IIN) abs_boundary_ispec_crust_mantle
          read(IIN) abs_boundary_npoin_crust_mantle
          read(IIN) abs_boundary_ijk_crust_mantle
          read(IIN) abs_boundary_jacobian2Dw_crust_mantle
          read(IIN) abs_boundary_normal_crust_mantle
        else
          ! dummy arrays
          allocate(abs_boundary_ispec_crust_mantle(1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
          allocate(abs_boundary_ijk_crust_mantle(1,1,1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
          allocate(abs_boundary_jacobian2Dw_crust_mantle(1,1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
          allocate(abs_boundary_normal_crust_mantle(1,1,1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
          allocate(abs_boundary_npoin_crust_mantle(1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
          abs_boundary_ispec_crust_mantle(:) = 0; abs_boundary_npoin_crust_mantle(:) = 0
          abs_boundary_ijk_crust_mantle(:,:,:) = 0
          abs_boundary_jacobian2Dw_crust_mantle(:,:) = 0.0; abs_boundary_normal_crust_mantle(:,:,:) = 0.0
        endif

        close(IIN)
      endif

      ! outer core
      if (NSPEC_OUTER_CORE > 0) then
        ! create name of database
        call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

        ! read arrays for Stacey conditions
        open(unit=IIN,file=prname(1:len_trim(prname))//'stacey.bin', &
              status='old',form='unformatted',action='read',iostat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey.bin file for outer core')

        read(IIN) num_abs_boundary_faces_outer_core

        ! allocates absorbing boundary arrays
        if (num_abs_boundary_faces_outer_core > 0) then
          allocate(abs_boundary_ispec_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
          allocate(abs_boundary_ijk_outer_core(3,NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
          allocate(abs_boundary_jacobian2Dw_outer_core(NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
          allocate(abs_boundary_npoin_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
          if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

          ! reads in arrays
          read(IIN) abs_boundary_ispec_outer_core
          read(IIN) abs_boundary_npoin_outer_core
          read(IIN) abs_boundary_ijk_outer_core
          read(IIN) abs_boundary_jacobian2Dw_outer_core
        else
          ! dummy arrays
          allocate(abs_boundary_ispec_outer_core(1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
          allocate(abs_boundary_ijk_outer_core(1,1,1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
          allocate(abs_boundary_jacobian2Dw_outer_core(1,1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
          allocate(abs_boundary_npoin_outer_core(1),stat=ier)
          if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
          abs_boundary_ispec_outer_core(:) = 0; abs_boundary_npoin_outer_core(:) = 0
          abs_boundary_ijk_outer_core(:,:,:) = 0
          abs_boundary_jacobian2Dw_outer_core(:,:) = 0.0
        endif
        close(IIN)
      endif
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

  ! checks if anything to do; setup only needed for kernel simulations
  if (SIMULATION_TYPE /= 3) return

  ! checks setup
  if (ADIOS_FOR_KERNELS) &
    call exit_mpi(myrank,'saving regular kernels in ADIOS file format is not supported yet')
  ! assuming 6 chunks full global simulations right now
  if (NCHUNKS_VAL /= 6 .or. NPROC_XI_VAL /= NPROC_ETA_VAL) &
    call exit_MPI(myrank, 'Only deal with 6 chunks at this moment')

  ! allocates arrays
  allocate(points_slice_reg(NM_KL_REG_PTS), &
           ispec_reg(NM_KL_REG_PTS), &
           hxir_reg(NGLLX, NM_KL_REG_PTS), &
           hetar_reg(NGLLY, NM_KL_REG_PTS), &
           hgammar_reg(NGLLZ, NM_KL_REG_PTS),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays points_slice_reg,..'
  points_slice_reg(:) = 0; ispec_reg(:) = 0
  hxir_reg(:,:) = 0.0; hetar_reg(:,:) = 0.0; hgammar_reg(:,:) = 0.0

  allocate(kl_reg_grid%rlayer(NM_KL_REG_LAYER), &
           kl_reg_grid%ndoubling(NM_KL_REG_LAYER), &
           kl_reg_grid%nlat(NM_KL_REG_LAYER), &
           kl_reg_grid%nlon(NM_KL_REG_LAYER), &
           kl_reg_grid%npts_before_layer(NM_KL_REG_LAYER+1),stat=ier)
  if (ier /= 0) stop 'Error allocating kl_reg_grid arrays'
  kl_reg_grid%rlayer(:) = 0.0; kl_reg_grid%ndoubling(:) = 0
  kl_reg_grid%nlat(:) = 0; kl_reg_grid%nlon(:) = 0
  kl_reg_grid%npts_before_layer(:) = 0

  ! reads in mesh inputs from file (see: PATHNAME_KL_REG = 'DATA/kl_reg_grid.txt' in constants.h)
  call read_kl_regular_grid(kl_reg_grid)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  locating regular kernel grid points...'
    call flush_IMAIN()
  endif

  if (myrank == 0) then
    ! main process
    allocate(slice_number(kl_reg_grid%npts_total),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating slice_number array')
    slice_number(:) = 0

    ! print *, 'slice npts =', kl_reg_grid%npts_total
    call find_regular_grid_slice_number(slice_number, kl_reg_grid)

    do i = NPROCTOT_VAL-1,0,-1
      npoints_slice_reg = 0
      do isp = 1,kl_reg_grid%npts_total
        if (slice_number(isp) == i) then
          npoints_slice_reg = npoints_slice_reg + 1
          if (npoints_slice_reg > NM_KL_REG_PTS) stop 'Exceeding NM_KL_REG_PTS limit'
          points_slice_reg(npoints_slice_reg) = isp
        endif
      enddo

      if (i /= 0) then
        call send_singlei(npoints_slice_reg,i,i)
        if (npoints_slice_reg > 0) then
          call send_i(points_slice_reg,npoints_slice_reg,i,2*i)
        endif
      endif
    enddo

    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/kl_grid_slice.txt',status='unknown',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening file kl_grid_slice.txt for writing')
    write(IOUT,*) kl_reg_grid%npts_total
    write(IOUT,*) slice_number(:)
    close(IOUT)

    deallocate(slice_number)
  else
    ! secondary processes
    call recv_singlei(npoints_slice_reg,0,myrank)
    if (npoints_slice_reg > 0) then
      call recv_i(points_slice_reg,npoints_slice_reg,0,2*myrank)
    endif
  endif

  ! this is the core part that takes up most of the computation time,
  ! and presumably the more processors involved the faster.
  if (npoints_slice_reg > 0) then
    call locate_regular_points(npoints_slice_reg, points_slice_reg, kl_reg_grid, &
                               NSPEC_CRUST_MANTLE, &
                               xstore_crust_mantle, ystore_crust_mantle, zstore_crust_mantle, &
                               ibool_crust_mantle, &
                               xigll, yigll, zigll, &
                               ispec_reg, hxir_reg, hetar_reg, hgammar_reg)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Finished locating kernel output regular grid'
    write(IMAIN,*) ' '
    call flush_IMAIN()
  endif

  end subroutine read_mesh_databases_regular_kl
