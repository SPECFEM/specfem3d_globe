!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, November 2007
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
              rho_vp,rho_vs,xstore,ystore,zstore, &
              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
              rhostore, kappavstore,muvstore,kappahstore,muhstore,eta_anisostore, &
              nspec_iso,nspec_tiso,nspec_ani, &
              c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
              c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
              c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
              ibool,idoubling,rmass,rmass_ocean_load,nspec,nglob, &
              READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY, &
              ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS,LOCAL_PATH,ABSORBING_CONDITIONS)

  implicit none

  include "constants.h"

  include "OUTPUT_FILES/values_from_mesher.h"

  integer iregion_code,myrank

! flags to know if we should read Vs and anisotropy arrays
  logical READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS,ABSORBING_CONDITIONS

  character(len=150) LOCAL_PATH

  integer nspec,nglob

  integer nspec_iso,nspec_tiso,nspec_ani

  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

! material properties
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec_iso)
  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec_iso)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec_iso)

! additional arrays for anisotropy stored only where needed to save memory
  real(kind=CUSTOM_REAL) kappahstore(NGLLX,NGLLY,NGLLZ,nspec_tiso)
  real(kind=CUSTOM_REAL) muhstore(NGLLX,NGLLY,NGLLZ,nspec_tiso)
  real(kind=CUSTOM_REAL) eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec_tiso)

! additional arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store,c33store,c34store, &
    c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store

! Stacey
  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,nspec)

! mass matrix and additional ocean load mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass,rmass_ocean_load

! global addressing
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

! processor identification
  character(len=150) prname

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_1.bin',status='old',action='read',form='unformatted')

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

  if(READ_KAPPA_MU) read(IIN) muvstore

! for anisotropy, gravity and rotation

  if(TRANSVERSE_ISOTROPY .and. READ_TISO) then
    read(IIN) kappahstore
    read(IIN) muhstore
    read(IIN) eta_anisostore
  endif

  if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
    read(IIN) c11store
    read(IIN) c12store
    read(IIN) c13store
    read(IIN) c33store
    read(IIN) c44store
  endif

  if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
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
  if(ABSORBING_CONDITIONS) then

    if(iregion_code == IREGION_CRUST_MANTLE) then
      read(IIN) rho_vp
      read(IIN) rho_vs
    else if(iregion_code == IREGION_OUTER_CORE) then
      read(IIN) rho_vp
    endif

  endif

! mass matrix
  read(IIN) rmass

! read additional ocean load mass matrix
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) read(IIN) rmass_ocean_load
  
  close(IIN)

! read coordinates of the mesh

  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_2.bin',status='old',action='read',form='unformatted')
  read(IIN) xstore
  read(IIN) ystore
  read(IIN) zstore

  read(IIN) ibool

  read(IIN) idoubling

  close(IIN)

  end subroutine read_arrays_solver

