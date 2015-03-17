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

! read arrays created by the mesher

  subroutine read_arrays_solver(iregion_code,myrank, &
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

  use constants_solver
  use specfem_par,only: &
    ABSORBING_CONDITIONS, &
    LOCAL_PATH,ABSORBING_CONDITIONS,&
    EXACT_MASS_MATRIX_FOR_ROTATION, mygroup, NUMBER_OF_SIMULTANEOUS_RUNS

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


  ! local parameters
  integer :: ier,lnspec,lnglob
  ! processor identification
  character(len=MAX_STRING_LEN) :: prname, path_to_add

  ! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)
  if (NUMBER_OF_SIMULTANEOUS_RUNS > 1 .and. mygroup >= 0) then
    write(path_to_add,"('run',i4.4,'/')") mygroup + 1
    prname = path_to_add(1:len_trim(path_to_add))//prname(1:len_trim(prname))
  endif

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data.bin', &
        status='old',action='read',form='unformatted',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening solver_data.bin')

  ! read coordinates of the mesh

  read(IIN) lnspec
  read(IIN) lnglob

  ! checks dimensions
  if (lnspec /= nspec) then
    close(IIN)
    print*,'Error file dimension: nspec in file = ',lnspec,' but nspec desired:',nspec
    print*,'please check file ',prname(1:len_trim(prname))//'solver_data.bin'
    call exit_mpi(myrank,'Error dimensions in solver_data.bin')
  endif
  if (lnglob /= nglob) then
    close(IIN)
    print*,'Error file dimension: nglob in file = ',lnglob,' but nglob desired:',nglob
    print*,'please check file ',prname(1:len_trim(prname))//'solver_data.bin'
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

  ! model arrays
  read(IIN) rhostore
  read(IIN) kappavstore

  if (READ_KAPPA_MU) read(IIN) muvstore

  ! for anisotropy, gravity and rotation
  if (TRANSVERSE_ISOTROPY_VAL .and. READ_TISO) then
    read(IIN) kappahstore
    read(IIN) muhstore
    read(IIN) eta_anisostore
  endif

  if (ANISOTROPIC_INNER_CORE_VAL .and. iregion_code == IREGION_INNER_CORE) then
    read(IIN) c11store
    read(IIN) c12store
    read(IIN) c13store
    read(IIN) c33store
    read(IIN) c44store
  endif

  if (ANISOTROPIC_3D_MANTLE_VAL .and. iregion_code == IREGION_CRUST_MANTLE) then
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
  endif

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
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) .and. iregion_code == IREGION_INNER_CORE)) then
    read(IIN) rmassx
    read(IIN) rmassy
  endif

  read(IIN) rmassz

  if (((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) .and. iregion_code == IREGION_CRUST_MANTLE) .or. &
      ((ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION) .and. iregion_code == IREGION_INNER_CORE)) then
    read(IIN) b_rmassx
    read(IIN) b_rmassy
  endif

  ! read additional ocean load mass matrix
  if (OCEANS_VAL .and. iregion_code == IREGION_CRUST_MANTLE) read(IIN) rmass_ocean_load

  close(IIN) ! solver_data.bin


  end subroutine read_arrays_solver

