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

! read arrays created by the mesher

  subroutine read_arrays_solver(iregion_code, &
                                nspec,nglob,nglob_xy, &
                                nspec_iso,nspec_tiso,nspec_ani, &
                                rho_vp,rho_vs, &
                                xstore,ystore,zstore, &
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
  use specfem_par, only: &
    ABSORBING_CONDITIONS, &
    LOCAL_PATH,ABSORBING_CONDITIONS

  implicit none

  integer,intent(in) :: iregion_code
  integer,intent(in) :: nspec,nglob,nglob_xy
  integer,intent(in) :: nspec_iso,nspec_tiso,nspec_ani

  ! Stacey
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: rho_vp,rho_vs

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! material properties
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_iso),intent(inout) :: &
    rhostore,kappavstore,muvstore

  ! additional arrays for anisotropy stored only where needed to save memory
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_tiso),intent(inout) :: &
    kappahstore,muhstore,eta_anisostore

  ! additional arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani),intent(inout) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store,c33store,c34store, &
    c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: mu0store

  ! global addressing
  integer,dimension(NGLLX,NGLLY,NGLLZ,nspec),intent(inout) :: ibool
  integer, dimension(nspec),intent(inout) :: idoubling
  logical, dimension(nspec),intent(inout) :: ispec_is_tiso

  ! mass matrices and additional ocean load mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob_xy),intent(inout) :: rmassx,rmassy
  real(kind=CUSTOM_REAL), dimension(nglob_xy),intent(inout) :: b_rmassx,b_rmassy

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout)    :: rmassz

  integer,intent(in) :: nglob_oceans
  real(kind=CUSTOM_REAL), dimension(nglob_oceans),intent(inout) :: rmass_ocean_load

  ! local parameters
  integer :: ier,lnspec,lnglob
  ! processor identification
  character(len=MAX_STRING_LEN) :: prname

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data.bin', &
        status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data.bin')

  ! read coordinates of the mesh

  read(IIN) lnspec
  read(IIN) lnglob

  ! checks dimensions
  if (lnspec /= nspec) then
    close(IIN)
    print *,'Error file dimension: nspec in file = ',lnspec,' but nspec desired:',nspec
    print *,'please check file ',prname(1:len_trim(prname))//'solver_data.bin'
    call exit_mpi(myrank,'Error dimensions in solver_data.bin')
  endif
  if (lnglob /= nglob) then
    close(IIN)
    print *,'Error file dimension: nglob in file = ',lnglob,' but nglob desired:',nglob
    print *,'please check file ',prname(1:len_trim(prname))//'solver_data.bin'
    call exit_mpi(myrank,'Error dimensions in solver_data.bin')
  endif

  ! mesh coordinates
  read(IIN) xstore
  read(IIN) ystore
  read(IIN) zstore

  read(IIN) ibool
  read(IIN) idoubling
  read(IIN) ispec_is_tiso

  ! local GLL points
  read(IIN) xix
  read(IIN) xiy
  read(IIN) xiz
  read(IIN) etax
  read(IIN) etay
  read(IIN) etaz
  read(IIN) gammax
  read(IIN) gammay
  read(IIN) gammaz

  ! checks if anything else to do for infinite regions
  if (iregion_code == IREGION_TRINFINITE .or. iregion_code == IREGION_INFINITE) then
    close(IIN)
    return
  endif

  ! model arrays
  read(IIN) rhostore
  read(IIN) kappavstore
  if (iregion_code /= IREGION_OUTER_CORE) then
    read(IIN) muvstore
  endif

  ! solid regions
  select case(iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle
    if (ANISOTROPIC_3D_MANTLE_VAL) then
      read(IIN) c11store
      read(IIN) c12store
      read(IIN) c13store
      read(IIN) c14store
      read(IIN) c15store
      read(IIN) c16store
      read(IIN) c22store
      read(IIN) c23store
      read(IIN) c24store
      read(IIN) c25store
      read(IIN) c26store
      read(IIN) c33store
      read(IIN) c34store
      read(IIN) c35store
      read(IIN) c36store
      read(IIN) c44store
      read(IIN) c45store
      read(IIN) c46store
      read(IIN) c55store
      read(IIN) c56store
      read(IIN) c66store
    else
      ! for anisotropy, gravity and rotation
      if (TRANSVERSE_ISOTROPY_VAL) then
        read(IIN) kappahstore
        read(IIN) muhstore
        read(IIN) eta_anisostore
      endif
    endif
    ! for azi kernels
    read(IIN) mu0store

  case (IREGION_INNER_CORE)
    ! inner core
    if (ANISOTROPIC_INNER_CORE_VAL) then
      read(IIN) c11store
      read(IIN) c12store
      read(IIN) c13store
      read(IIN) c33store
      read(IIN) c44store
    endif
  end select

  ! Stacey
  if (ABSORBING_CONDITIONS) then
    if (iregion_code == IREGION_CRUST_MANTLE) then
      read(IIN) rho_vp
      read(IIN) rho_vs
    else if (iregion_code == IREGION_OUTER_CORE) then
      read(IIN) rho_vp
    endif
  endif

  ! mass matrices
  !
  ! in the case of Stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  if (((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_INNER_CORE)) then
    read(IIN) rmassx
    read(IIN) rmassy
  endif

  read(IIN) rmassz

  if (((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION_VAL) .and. iregion_code == IREGION_INNER_CORE)) then
    read(IIN) b_rmassx
    read(IIN) b_rmassy
  endif

  ! read additional ocean load mass matrix
  if (OCEANS_VAL .and. iregion_code == IREGION_CRUST_MANTLE) read(IIN) rmass_ocean_load

  close(IIN) ! solver_data.bin

  end subroutine read_arrays_solver

