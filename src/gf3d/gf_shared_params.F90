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

!----
!---- The `shared_parameters` bootstrap.
!----
!---- src/gf3d/ reuses routines from src/shared/ rather than re-porting
!---- specfem's coordinate and mesh mathematics. Those routines read
!---- module state out of `shared_parameters`, which in the solver is
!---- filled by read_parameter_file() from DATA/Par_file. A standalone
!---- library has no Par_file, so nothing initialises them.
!----
!---- Most of the affected variables have sane declaration-time defaults.
!---- Three groups do not, and each of them fails *silently*:
!----
!----   * NUMBER_OF_SIMULTANEOUS_RUNS and NOISE_TOMOGRAPHY have no
!----     initialiser at all (shared_par.f90:161 and :88), so they hold
!----     whatever the loader left there. get_cmt.f90:79 branches on the
!----     first and get_cmt.f90:387 zeroes every half duration when the
!----     second is non-zero.
!----
!----   * NX_BATHY / NY_BATHY / RESOLUTION_TOPO_FILE index the topography
!----     grid inside get_topo_bathy(). Their Earth defaults do not match
!----     a database written for another planet, or with another topo
!----     file, and get_topo_bathy() does not range-check: it returns
!----     plausible-looking garbage.
!----
!----   * PATHNAME_TOPO_FILE is compared against PATHNAME_TOPO_FILE_BERKELEY
!----     inside get_topo_bathy() to pick between the plain bilinear
!----     interpolation and the Berkeley smoothed branch. Equal, and the
!----     library silently takes a branch meant for a different grid.
!----
!---- Everything here is set from the database, so the library stays
!---- planet-agnostic instead of being quietly Earth-only.
!----

  module gf_shared_params

  use gf_par, only: t_gfdb,gf_set_error,GF_OK,GF_ERR_MISMATCH

  implicit none

  private

  public :: gf_init_shared_params

  contains

!
!-------------------------------------------------------------------------------------------------
!

  subroutine gf_init_shared_params(db,ierr)

! fills the parts of `shared_parameters` that reused src/shared/ routines rely on
!
! Called by gf_open() once the database metadata has been read.

  use constants, only: PATHNAME_TOPO_FILE_BERKELEY,EARTH_PATHNAME_TOPO_FILE,myrank

  use shared_parameters, only: NUMBER_OF_SIMULTANEOUS_RUNS,NOISE_TOMOGRAPHY, &
                               NX_BATHY,NY_BATHY,RESOLUTION_TOPO_FILE,PATHNAME_TOPO_FILE, &
                               TOPO_MINIMUM,TOPO_MAXIMUM, &
                               OUTPUT_FILES,R_PLANET,R_PLANET_KM,R_EARTH,R_EARTH_KM,RHOAV

  implicit none

  type(t_gfdb), intent(in) :: db
  integer, intent(out) :: ierr

  ! this is a serial library: there is exactly one process, and it is rank 0
  myrank = 0

  ! get_cmt.f90:79 branches on this; it has no initialiser in shared_par.f90
  NUMBER_OF_SIMULTANEOUS_RUNS = 1

  ! get_cmt.f90:387 zeroes every source half duration when this is non-zero;
  ! it has no initialiser in shared_par.f90 either
  NOISE_TOMOGRAPHY = 0

  ! used by incidental writers in reused routines
  OUTPUT_FILES = '.'

  ! planet constants, from the database rather than the Earth defaults
  if (db%R_PLANET > 0.d0) then
    R_PLANET    = db%R_PLANET
    R_PLANET_KM = db%R_PLANET / 1000.d0
    ! several reused routines still read the older R_EARTH names
    R_EARTH     = db%R_PLANET
    R_EARTH_KM  = db%R_PLANET / 1000.d0
  endif

  ! note: RHOAV is not written to mesh_info.h5 by the current writer, so
  !       db%RHOAV is zero for the shipped example databases and the Earth
  !       default declared in shared_par.f90 stands. It first matters in
  !       Stage 4, where get_cmt() forms scaleM from it. gf_open() reads the
  !       attribute if a future writer adds it.
  if (db%RHOAV > 0.d0) RHOAV = db%RHOAV

  ! topography grid geometry: get_topo_bathy() indexes ibathy_topo with these
  if (db%topography) then
    if (db%NX_BATHY <= 0 .or. db%NY_BATHY <= 0 .or. db%RESOLUTION_TOPO_FILE <= 0.d0) then
      call gf_set_error(ierr,GF_ERR_MISMATCH, &
        'database has TOPOGRAPHY set but no usable NX_BATHY/NY_BATHY/RESOLUTION_TOPO_FILE')
      return
    endif
    NX_BATHY = db%NX_BATHY
    NY_BATHY = db%NY_BATHY
    RESOLUTION_TOPO_FILE = db%RESOLUTION_TOPO_FILE
  endif

  ! get_topo_bathy() takes the Berkeley smoothed branch when
  ! PATHNAME_TOPO_FILE == PATHNAME_TOPO_FILE_BERKELEY. The topography grid
  ! comes out of mesh_info.h5, never off disk, so any other value will do —
  ! but it must not be that one.
  PATHNAME_TOPO_FILE = EARTH_PATHNAME_TOPO_FILE
  if (PATHNAME_TOPO_FILE == PATHNAME_TOPO_FILE_BERKELEY) then
    call gf_set_error(ierr,GF_ERR_MISMATCH, &
      'PATHNAME_TOPO_FILE default matches the Berkeley smoothed topography path; ' // &
      'get_topo_bathy would take the wrong branch')
    return
  endif

  ! range used when validating a topography grid read off disk; we never read
  ! one, but keep the values consistent with the Earth grid in the database
  if (TOPO_MINIMUM >= TOPO_MAXIMUM) then
    call gf_set_error(ierr,GF_ERR_MISMATCH,'inconsistent TOPO_MINIMUM/TOPO_MAXIMUM defaults')
    return
  endif

  ierr = GF_OK

  end subroutine gf_init_shared_params

  end module gf_shared_params
