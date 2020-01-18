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

!===============================================================================
!> \brief Read adios arrays created by the mesher (file: regX_solver_data.bp)
subroutine read_arrays_solver_adios(iregion_code, &
                                    nspec,nglob,nglob_xy, &
                                    nspec_iso,nspec_tiso,nspec_ani, &
                                    rho_vp,rho_vs,xstore,ystore,zstore, &
                                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                    rhostore, kappavstore,muvstore,kappahstore,muhstore,eta_anisostore, &
                                    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                    mu0store, &
                                    ibool,idoubling,ispec_is_tiso, &
                                    rmassx,rmassy,rmassz, &
                                    nglob_oceans,rmass_ocean_load, &
                                    b_rmassx,b_rmassy)


  use constants_solver
  use specfem_par, only: ABSORBING_CONDITIONS,LOCAL_PATH,ABSORBING_CONDITIONS

  use adios_helpers_mod
  use manager_adios

  implicit none

  integer :: iregion_code
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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: mu0store

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

  character(len=MAX_STRING_LEN) :: file_name

  ! local parameters
  integer :: lnspec, lnglob, local_dim
  ! ADIOS variables
  integer(kind=8), dimension(1) :: start, count

  integer(kind=8), dimension(256),target :: selections
  integer :: sel_num, i
  integer(kind=8), pointer :: sel => null()

  character(len=128)      :: region_name, region_name_scalar

  ! user output
  if (myrank == 0) then
#if defined(USE_ADIOS)
    write(IMAIN,*) '  reading arrays solver in ADIOS 1 file format'
#elif defined(USE_ADIOS2)
    write(IMAIN,*) '  reading arrays solver in ADIOS 2 file format'
#endif
    call flush_IMAIN()
  endif

  write(region_name,"('reg',i1, '/')") iregion_code
  write(region_name_scalar,"('reg',i1)") iregion_code

  sel_num = 0

  file_name= trim(LOCAL_PATH) // "/solver_data.bp"

  ! Setup the ADIOS library to read the file
  call open_file_adios_read_and_init_method(myadios_file,myadios_group,file_name)

  ! debug
  !call show_adios_file_variables(myadios_file,myadios_group,file_name)

  ! read coordinates of the mesh
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nspec",lnspec)
  call read_adios_scalar(myadios_file,myadios_group,myrank,trim(region_name) // "nglob",lnglob)

  ! mesh coordinates
  local_dim = nglob
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "x_global/array", xstore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "y_global/array", ystore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "z_global/array", zstore)

  local_dim = nspec
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "idoubling/array", idoubling)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "ispec_is_tiso/array", ispec_is_tiso)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "ibool/array", ibool)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "ibool/array", ibool)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "xixstore/array", xix)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "xiystore/array", xiy)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "xizstore/array", xiz)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "etaxstore/array", etax)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "etaystore/array", etay)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "etazstore/array", etaz)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "gammaxstore/array", gammax)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "gammaystore/array", gammay)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "gammazstore/array", gammaz)

  local_dim = NGLLX * NGLLY * NGLLZ * nspec_iso
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "rhostore/array", rhostore)
  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                 trim(region_name) // "kappavstore/array", kappavstore)
  if (iregion_code /= IREGION_OUTER_CORE) then
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "muvstore/array", muvstore)
  endif

  ! solid regions
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle
    if (ANISOTROPIC_3D_MANTLE_VAL) then
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
      start(1) = local_dim*myrank; count(1) = local_dim
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c11store/array", c11store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c12store/array", c12store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c13store/array", c13store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c14store/array", c14store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c15store/array", c15store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c16store/array", c16store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c22store/array", c22store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c23store/array", c23store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c24store/array", c24store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c25store/array", c25store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c26store/array", c26store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c33store/array", c33store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c34store/array", c34store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c35store/array", c35store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c36store/array", c36store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c44store/array", c44store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c45store/array", c45store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c46store/array", c46store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c55store/array", c55store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c56store/array", c56store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c66store/array", c66store)
    else
      if (TRANSVERSE_ISOTROPY_VAL) then
        local_dim = NGLLX * NGLLY * NGLLZ * nspec_tiso
        start(1) = local_dim*myrank; count(1) = local_dim
        sel_num = sel_num+1
        sel => selections(sel_num)
        call set_selection_boundingbox(sel, start, count)

        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                       trim(region_name) // "kappahstore/array", kappahstore)
        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                       trim(region_name) // "muhstore/array", muhstore)
        call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                       trim(region_name) // "eta_anisostore/array", eta_anisostore)
      endif
    endif

    ! for azimuthal
    local_dim = NGLLX * NGLLY * NGLLZ * nspec
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "mu0store/array", mu0store)

  case (IREGION_INNER_CORE)
    ! inner core
    if (ANISOTROPIC_INNER_CORE_VAL) then
      local_dim = NGLLX * NGLLY * NGLLZ * nspec_ani
      start(1) = local_dim*myrank; count(1) = local_dim
      sel_num = sel_num+1
      sel => selections(sel_num)
      call set_selection_boundingbox(sel, start, count)

      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c11store/array", c11store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c12store/array", c12store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c13store/array", c13store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c33store/array", c33store)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "c44store/array", c44store)
    endif
  end select

  ! Stacey
  if (ABSORBING_CONDITIONS) then
    local_dim = NGLLX * NGLLY * NGLLZ * nspec ! nspec_stacey in meshfem3D
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    if (iregion_code == IREGION_CRUST_MANTLE) then
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "rho_vp/array", rho_vp)
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "rho_vs/array", rho_vs)
    else if (iregion_code == IREGION_OUTER_CORE) then
      call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                     trim(region_name) // "rho_vp/array", rho_vp)
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
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL .and. iregion_code == IREGION_INNER_CORE)) then

    local_dim = nglob_xy
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "rmassx/array", rmassx)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "rmassy/array", rmassy)
  endif

  local_dim = nglob
  start(1) = local_dim*myrank; count(1) = local_dim
  sel_num = sel_num+1
  sel => selections(sel_num)
  call set_selection_boundingbox(sel, start, count)

  call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, trim(region_name) // "rmassz/array", rmassz)


  if ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL .and. iregion_code == IREGION_INNER_CORE)) then
    local_dim = nglob_xy
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "b_rmassx/array", b_rmassx)
    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "b_rmassy/array", b_rmassy)
  endif

  ! read additional ocean load mass matrix
  if (OCEANS_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
    local_dim = NGLOB_CRUST_MANTLE_OCEANS ! nglob_oceans
    start(1) = local_dim*myrank; count(1) = local_dim
    sel_num = sel_num+1
    sel => selections(sel_num)
    call set_selection_boundingbox(sel, start, count)

    call read_adios_schedule_array(myadios_file, myadios_group, sel, start, count, &
                                   trim(region_name) // "rmass_ocean_load/array", rmass_ocean_load)
  endif

  ! perform actual reading
  call read_adios_perform(myadios_file)

  ! Clean everything and close the ADIOS file
  do i = 1, sel_num
    sel => selections(i)
    call delete_adios_selection(sel)
  enddo

  ! closes adios file
  call close_file_adios_read_and_finalize_method(myadios_file)

  ! checks dimensions
  if (lnspec /= nspec) then
    print *,'Error file dimension: nspec in file = ',lnspec, &
        ' but nspec desired:',nspec
    print *,'please check file ', file_name
    call exit_mpi(myrank,'Error dimensions in solver_data.bp')
  endif
  if (lnglob /= nglob) then
    print *,'Error file dimension: nglob in file = ',lnglob, &
        ' but nglob desired:',nglob
    print *,'please check file ', file_name
    call exit_mpi(myrank,'Error dimensions in solver_data.bp')
  endif

end subroutine read_arrays_solver_adios
