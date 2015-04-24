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
                                    rmassx,rmassy,rmassz, &
                                    nglob_oceans,rmass_ocean_load, &
                                    READ_KAPPA_MU,READ_TISO, &
                                    b_rmassx,b_rmassy)

  use adios_read_mod
  use adios_helpers_mod, only: check_adios_err

  use constants_solver
  use specfem_par,only: &
    ABSORBING_CONDITIONS, &
    LOCAL_PATH,ABSORBING_CONDITIONS,&
    EXACT_MASS_MATRIX_FOR_ROTATION, NUMBER_OF_SIMULTANEOUS_RUNS, mygroup

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

  integer :: nglob_oceans
  real(kind=CUSTOM_REAL), dimension(nglob_oceans) :: rmass_ocean_load

  ! flags to know if we should read Vs and anisotropy arrays
  logical :: READ_KAPPA_MU,READ_TISO

  character(len=MAX_STRING_LEN) :: file_name, path_to_add

  ! local parameters
  integer :: comm, lnspec, lnglob, local_dim
  ! ADIOS variables
  integer                 :: adios_err
  integer(kind=8)         :: adios_handle
  integer(kind=8), dimension(1) :: start, count

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num, i
  integer(kind=8), pointer :: sel => null()

  character(len=128)      :: region_name, region_name_scalar

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  sel_num = 0

  file_name= trim(LOCAL_PATH) // "/solver_data.bp"
  ! Append the actual file name.
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    file_name = path_to_add(1:len_trim(path_to_add))//file_name(1:len_trim(file_name))
  endif

  call world_duplicate(comm)

  ! Setup the ADIOS library to read the file
  call adios_read_init_method (ADIOS_READ_METHOD_BP, comm, "verbose=1", adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_read_open_file (adios_handle, file_name, 0, comm, adios_err)
  if (adios_err /= 0) then
    print*,'Error rank ',myrank,' opening adios file: ',trim(file_name)
    call check_adios_err(myrank,adios_err)
  endif

  ! read coordinates of the mesh
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_writeblock(sel, myrank)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "nspec", 0, 1, &
     lnspec, adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "nglob", 0, 1, &
     lnglob, adios_err)

  ! mesh coordinates
  local_dim = nglob
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "x_global/array", 0, 1, &
      xstore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "y_global/array", 0, 1, &
      ystore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "z_global/array", 0, 1, &
      zstore, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_iso
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "rhostore/array", 0, 1, &
      rhostore, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "kappavstore/array", 0, 1, &
      kappavstore, adios_err)
  call check_adios_err(myrank,adios_err)
  if (READ_KAPPA_MU) then
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "muvstore/array", 0, 1, &
        muvstore, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  if (TRANSVERSE_ISOTROPY_VAL .and. READ_TISO) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_tiso
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "kappahstore/array", 0, 1, &
        kappahstore, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "muhstore/array", 0, 1, &
        muhstore, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "eta_anisostore/array", 0, 1, &
        eta_anisostore, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  local_dim = nspec
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "idoubling/array", 0, 1, &
      idoubling, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "ispec_is_tiso/array", 0, 1, &
      ispec_is_tiso, adios_err)
  call check_adios_err(myrank,adios_err)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)

  call adios_schedule_read(adios_handle, sel, trim(region_name) // "ibool/array", 0, 1, &
      ibool, adios_err)
  call check_adios_err(myrank,adios_err)

  call adios_schedule_read(adios_handle, sel, trim(region_name) // "xixstore/array", 0, 1, &
      xix, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "xiystore/array", 0, 1, &
      xiy, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "xizstore/array", 0, 1, &
      xiz, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "etaxstore/array", 0, 1, &
      etax, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "etaystore/array", 0, 1, &
      etay, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "etazstore/array", 0, 1, &
      etaz, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "gammaxstore/array", 0, 1, &
      gammax, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "gammaystore/array", 0, 1, &
      gammay, adios_err)
  call check_adios_err(myrank,adios_err)
  call adios_schedule_read(adios_handle, sel, trim(region_name) // "gammazstore/array", 0, 1, &
      gammaz, adios_err)
  call check_adios_err(myrank,adios_err)

  if (ANISOTROPIC_INNER_CORE_VAL .and. iregion_code == IREGION_INNER_CORE) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c11store/array", 0, 1, &
        c11store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c12store/array", 0, 1, &
        c12store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c13store/array", 0, 1, &
        c13store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c33store/array", 0, 1, &
        c33store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c44store/array", 0, 1, &
        c44store, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  if (ANISOTROPIC_3D_MANTLE_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c11store/array", 0, 1, &
        c11store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c12store/array", 0, 1, &
        c12store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c13store/array", 0, 1, &
        c13store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c14store/array", 0, 1, &
        c14store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c15store/array", 0, 1, &
        c15store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c16store/array", 0, 1, &
        c16store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c22store/array", 0, 1, &
        c22store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c23store/array", 0, 1, &
        c23store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c24store/array", 0, 1, &
        c24store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c25store/array", 0, 1, &
        c25store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c26store/array", 0, 1, &
        c26store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c33store/array", 0, 1, &
        c33store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c34store/array", 0, 1, &
        c34store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c35store/array", 0, 1, &
        c35store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c36store/array", 0, 1, &
        c36store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c44store/array", 0, 1, &
        c44store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c45store/array", 0, 1, &
        c45store, adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c46store/array", 0, 1, &
        c46store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c55store/array", 0, 1, &
        c55store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c56store/array", 0, 1, &
        c56store, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "c66store/array", 0, 1, &
        c66store, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! Stacey
  if (ABSORBING_CONDITIONS) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec ! nspec_stacey in meshfem3D
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    if (iregion_code == IREGION_CRUST_MANTLE) then
      call adios_schedule_read(adios_handle, sel, trim(region_name) // "rho_vp/array", 0, 1, &
          rho_vp, adios_err)
      call check_adios_err(myrank,adios_err)
      call adios_schedule_read(adios_handle, sel, trim(region_name) // "rho_vs/array", 0, 1, &
          rho_vs, adios_err)
      call check_adios_err(myrank,adios_err)
    else if (iregion_code == IREGION_OUTER_CORE) then
      call adios_schedule_read(adios_handle, sel, trim(region_name) // "rho_vp/array", 0, 1, &
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
  if ((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then

    local_dim = nglob_xy
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, trim(region_name) // "rmassx/array", 0, 1, &
        rmassx, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "rmassy/array", 0, 1, &
        rmassy, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  local_dim = nglob
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call adios_selection_boundingbox (sel , 1, start, count)

  call adios_schedule_read(adios_handle, sel, trim(region_name) // "rmassz/array", 0, 1, &
                           rmassz, adios_err)
  call check_adios_err(myrank,adios_err)


  if ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION .and. iregion_code == IREGION_INNER_CORE)) then
    local_dim = nglob_xy
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, trim(region_name) // "b_rmassx/array", 0, 1, &
                             b_rmassx, adios_err)
    call check_adios_err(myrank,adios_err)
    call adios_schedule_read(adios_handle, sel, trim(region_name) // "b_rmassy/array", 0, 1, &
                             b_rmassy, adios_err)
    call check_adios_err(myrank,adios_err)
  endif

  ! read additional ocean load mass matrix
  if (OCEANS_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = NGLOB_CRUST_MANTLE_OCEANS ! nglob_oceans
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call adios_selection_boundingbox (sel , 1, start, count)

    call adios_schedule_read(adios_handle, sel, trim(region_name) // "rmass_ocean_load/array", &
        0, 1, rmass_ocean_load, adios_err)
    call check_adios_err(myrank,adios_err)
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

  call synchronize_all_comm(comm)

  ! checks dimensions
  if (lnspec /= nspec) then
    print*,'Error file dimension: nspec in file = ',lnspec, &
        ' but nspec desired:',nspec
    print*,'please check file ', file_name
    call exit_mpi(myrank,'Error dimensions in solver_data.bp')
  endif
  if (lnglob /= nglob) then
    print*,'Error file dimension: nglob in file = ',lnglob, &
        ' but nglob desired:',nglob
    print*,'please check file ', file_name
    call exit_mpi(myrank,'Error dimensions in solver_data.bp')
  endif

end subroutine read_arrays_solver_adios
