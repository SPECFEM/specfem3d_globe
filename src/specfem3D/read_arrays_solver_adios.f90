!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

!===============================================================================
!> \brief Read adios arrays created by the mesher (file: regX_solver_data.bp)
subroutine read_arrays_solver_adios(iregion_code,myrank, &
              nspec,nglob,nglob_xy, &
              nspec_iso,nspec_tiso,nspec_ani, &
              rho_vp,rho_vs,xstore,ystore,zstore, &
              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
              rhostore, kappavstore,muvstore,kappahstore,muhstore,eta_anisostore, &
              c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
              c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
              c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
              ibool,idoubling,ispec_is_tiso, &
              rmassx,rmassy,rmassz,rmass_ocean_load, &
              READ_KAPPA_MU,READ_TISO, &
              b_rmassx,b_rmassy)

  use mpi
  use adios_read_mod

  use constants_solver
  use specfem_par,only: &
    ABSORBING_CONDITIONS,TRANSVERSE_ISOTROPY, &
    ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS,LOCAL_PATH,ABSORBING_CONDITIONS,&
    EXACT_MASS_MATRIX_FOR_ROTATION

  implicit none

  integer :: iregion_code,myrank
  integer :: nspec,nglob,nglob_xy
  integer :: nspec_iso,nspec_tiso,nspec_ani

  ! Stacey
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec):: rho_vp,rho_vs

  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! material properties
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_iso) :: &
    rhostore,kappavstore,muvstore

  ! additional arrays for anisotropy stored only where needed to save memory
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_tiso) :: &
    kappahstore,muhstore,eta_anisostore

  ! additional arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store,c33store,c34store, &
    c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! global addressing
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: idoubling
  logical, dimension(nspec) :: ispec_is_tiso

  ! mass matrices and additional ocean load mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob_xy) :: rmassx,rmassy
  real(kind=CUSTOM_REAL), dimension(nglob_xy) :: b_rmassx,b_rmassy

  real(kind=CUSTOM_REAL), dimension(nglob)    :: rmassz
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load

  ! flags to know if we should read Vs and anisotropy arrays
  logical :: READ_KAPPA_MU,READ_TISO

  character(len=150) :: file_name

  ! local parameters
  integer :: ierr, comm, lnspec, lnglob, local_dim
  ! processor identification
  character(len=150) :: prname
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_group, adios_handle, varid
  integer(kind=8)         :: adios_groupsize, adios_totalsize
  integer :: vars_count, attrs_count, current_step, last_step, vsteps
  character(len=128), dimension(:), allocatable :: adios_names
  integer(kind=8), dimension(1) :: start, count

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num, i
  integer(kind=8), pointer :: sel => null()

  sel_num = 0

  ! create a prefix for the file name such as LOCAL_PATH/regX_
  call create_name_database_adios(prname, iregion_code, LOCAL_PATH)

  ! Postpend the actual file name.
  file_name= trim(prname) // "solver_data.bp"
  call MPI_Comm_dup (MPI_COMM_WORLD, comm, ierr)

  ! Setup the ADIOS library to read the file
  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, &
      "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_open_file (adios_handle, file_name, 0, comm, ierr)
  call check_adios_err(myrank,adios_err)

  ! read coordinates of the mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, "nspec", 0, 1, &
     lnspec, adios_err)
  call adios_schedule_read(adios_handle, sel, "nglob", 0, 1, &
     lnglob, adios_err)
  !call adios_get_scalar(adios_handle, "nspec", lnspec, adios_err)
  !call adios_get_scalar(adios_handle, "nglob", lnglob, adios_err)
  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)


  ! mesh coordinates
  local_dim = nglob
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "xstore/array", 0, 1, &
      xstore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ystore/array", 0, 1, &
      ystore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "zstore/array", 0, 1, &
      zstore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "rmassz/array", 0, 1, &
      rmassz, adios_err)
  call check_adios_err(myrank,adios_err)

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_iso
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "rhostore/array", 0, 1, &
      rhostore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "kappavstore/array", 0, 1, &
      kappavstore, adios_err)
  call check_adios_err(myrank,adios_err)
  if(READ_KAPPA_MU) then
    call adios_schedule_read(adios_handle, sel, "muvstore/array", 0, 1, &
        muvstore, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  if(TRANSVERSE_ISOTROPY_VAL .and. READ_TISO) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_tiso
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, "kappahstore/array", 0, 1, &
        kappahstore, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "muhstore/array", 0, 1, &
        muhstore, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "eta_anisostore/array", 0, 1, &
        eta_anisostore, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  local_dim = nspec
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "idoubling/array", 0, 1, &
      idoubling, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "ispec_is_tiso/array", 0, 1, &
      ispec_is_tiso, adios_err)
  call check_adios_err(myrank,adios_err)

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, "ibool/array", 0, 1, &
      ibool, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "xixstore/array", 0, 1, &
      xix, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "xiystore/array", 0, 1, &
      xiy, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "xizstore/array", 0, 1, &
      xiz, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "etaxstore/array", 0, 1, &
      etax, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "etaystore/array", 0, 1, &
      etay, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "etazstore/array", 0, 1, &
      etaz, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "gammaxstore/array", 0, 1, &
      gammax, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "gammaystore/array", 0, 1, &
      gammay, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, "gammazstore/array", 0, 1, &
      gammaz, adios_err)
  call check_adios_err(myrank,adios_err)

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)


  if(ANISOTROPIC_INNER_CORE_VAL .and. iregion_code == IREGION_INNER_CORE) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, "c11store/array", 0, 1, &
        c11store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c12store/array", 0, 1, &
        c12store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c13store/array", 0, 1, &
        c13store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c33store/array", 0, 1, &
        c33store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c44store/array", 0, 1, &
        c44store, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  if(ANISOTROPIC_3D_MANTLE_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, "c11store/array", 0, 1, &
        c11store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c12store/array", 0, 1, &
        c12store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c13store/array", 0, 1, &
        c13store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c14store/array", 0, 1, &
        c14store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c15store/array", 0, 1, &
        c15store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c16store/array", 0, 1, &
        c16store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c22store/array", 0, 1, &
        c22store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c23store/array", 0, 1, &
        c23store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c24store/array", 0, 1, &
        c24store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c25store/array", 0, 1, &
        c25store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c26store/array", 0, 1, &
        c26store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c33store/array", 0, 1, &
        c33store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c34store/array", 0, 1, &
        c34store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c35store/array", 0, 1, &
        c35store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c36store/array", 0, 1, &
        c36store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c44store/array", 0, 1, &
        c44store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c45store/array", 0, 1, &
        c45store, adios_err)
    call adios_schedule_read(adios_handle, sel, "c46store/array", 0, 1, &
        c46store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c55store/array", 0, 1, &
        c55store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c56store/array", 0, 1, &
        c56store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "c66store/array", 0, 1, &
        c66store, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  ! Stacey
  if(ABSORBING_CONDITIONS) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec ! nspec_stacey in meshfem3D
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    if(iregion_code == IREGION_CRUST_MANTLE) then
      call adios_schedule_read(adios_handle, sel, "rho_vp/array", 0, 1, &
          rho_vp, adios_err)
      call check_adios_err(myrank,adios_err)
      call adios_schedule_read(adios_handle, sel, "rho_vs/array", 0, 1, &
          rho_vs, adios_err)
      call check_adios_err(myrank,adios_err)
    else if(iregion_code == IREGION_OUTER_CORE) then
      call adios_schedule_read(adios_handle, sel, "rho_vp/array", 0, 1, &
          rho_vp, adios_err)
      call check_adios_err(myrank,adios_err)
    endif

  endif

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to
  ! the mass matrix on Stacey edges for the crust_mantle and outer_core regions
  ! but not for the inner_core region thus the mass matrix must be replaced by
  ! three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix
  ! is needed for the sake of performance, only "rmassz" array will be filled
  ! and "rmassx" & "rmassy" will be obsolete
  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  if( (NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then

    local_dim = nglob_xy
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, "rmassx/array", 0, 1, &
        rmassx, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "rmassy/array", 0, 1, &
        rmassy, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  local_dim = nglob
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)

  call adios_schedule_read(adios_handle, sel, "rmassz/array", 0, 1, &
      rmassz, adios_err)
  call check_adios_err(myrank,adios_err)

  if( (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE))then
    local_dim = nglob_xy
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, "b_rmassx/array", 0, 1, &
        b_rmassx, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, "b_rmassy/array", 0, 1, &
        b_rmassy, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  !call adios_perform_reads(adios_handle, adios_err)
  !call check_adios_err(myrank,adios_err)

  ! read additional ocean load mass matrix
  if(OCEANS_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = NGLOB_CRUST_MANTLE_OCEANS ! nglob_oceans
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, "rmass_ocean_load/array", &
        0, 1, rmass_ocean_load, adios_err)
    call check_adios_err(myrank,adios_err)

    !call adios_perform_reads(adios_handle, adios_err)
    !call check_adios_err(myrank,adios_err)
  endif


  call adios_perform_reads(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)


  ! Clean everything and close the ADIOS file
  do i = 1, sel_num
    sel => selections(i)
    call adios_selection_delete(sel)
  enddo
  call adios_read_close(adios_handle, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_read_finalize_method(ADIOS_READ_METHOD_BP, adios_err)
  call check_adios_err(myrank,adios_err)

  call synchronize_all()
  ! checks dimensions
  if( lnspec /= nspec ) then
    print*,'error file dimension: nspec in file = ',lnspec, &
        ' but nspec desired:',nspec
    print*,'please check file ', file_name
    call exit_mpi(myrank,'error dimensions in solver_data.bp')
  endif
  if( lnglob /= nglob ) then
    print*,'error file dimension: nglob in file = ',lnglob, &
        ' but nglob desired:',nglob
    print*,'please check file ', file_name
    call exit_mpi(myrank,'error dimensions in solver_data.bp')
  endif

end subroutine read_arrays_solver_adios
