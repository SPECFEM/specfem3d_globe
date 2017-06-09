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


  subroutine prepare_movie()

  use specfem_par
  implicit none

  ! allocate files to save surface movies
  if (MOVIE_SURFACE) call prepare_movie_surface()

  ! output point and element information for 3D movies
  if (MOVIE_VOLUME) call prepare_movie_volume()

  end subroutine prepare_movie

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_movie_surface()

  use specfem_par
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: ier

  ! checks if anything to do
  if (.not. MOVIE_SURFACE) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing movie surface"
    call flush_IMAIN()
  endif

  ! only output corners
  ! note: for noise tomography, must NOT be coarse (have to be saved on all GLL points)
  if (MOVIE_COARSE) then
    ! checks setup
    if (NGLLX /= NGLLY) &
      call exit_MPI(myrank,'MOVIE_COARSE together with MOVIE_SURFACE requires NGLLX=NGLLY')
    ! number of points
    nmovie_points = 2 * 2 * NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    NIT = NGLLX - 1
  else
    ! number of points
    nmovie_points = NGLLX * NGLLY * NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    NIT = 1
  endif

  ! checks exact number of points nmovie_points
  call movie_surface_count_points()

  ! those arrays are not necessary for noise tomography, so only allocate them in MOVIE_SURFACE case
  ! writes out movie point locations to file
  call write_movie_surface_mesh()

  ! allocates movie surface arrays for wavefield values
  allocate(store_val_ux(nmovie_points), &
           store_val_uy(nmovie_points), &
           store_val_uz(nmovie_points),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface arrays')

  ! allocates arrays for gathering wavefield values
  if (myrank == 0) then
    ! only master needs full arrays
    allocate(store_val_ux_all(nmovie_points,0:NPROCTOT_VAL-1), &
             store_val_uy_all(nmovie_points,0:NPROCTOT_VAL-1), &
             store_val_uz_all(nmovie_points,0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface all arrays')
  else
    ! slave processes only need dummy arrays
    allocate(store_val_ux_all(1,1), &
             store_val_uy_all(1,1), &
             store_val_uz_all(1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface all arrays')
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  Movie surface:'
    write(IMAIN,*) '  Writing to moviedata*** files in output directory'
    if (MOVIE_VOLUME_TYPE == 5) then
      write(IMAIN,*) '  movie output: displacement'
    else
      write(IMAIN,*) '  movie output: velocity'
    endif
    write(IMAIN,*) '  time steps every: ',NTSTEP_BETWEEN_FRAMES
    call flush_IMAIN()
  endif

  call synchronize_all()

  end subroutine prepare_movie_surface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_movie_volume()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie
  implicit none

  ! local parameters
  integer :: ier

  ! checks if anything to do
  if (.not. MOVIE_VOLUME) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing movie volume"
    call flush_IMAIN()
  endif

  ! checks
  ! the following has to be true for the the array dimensions of eps to match with those of rstore etc..
  ! note that epsilondev and eps_trace_over_3 don't have the same dimensions.. could cause trouble
  if (NSPEC_CRUST_MANTLE_STR_OR_ATT /= NSPEC_CRUST_MANTLE) &
    stop 'NSPEC_CRUST_MANTLE_STRAINS_ATT /= NSPEC_CRUST_MANTLE'
  if (NSPEC_CRUST_MANTLE_STRAIN_ONLY /= NSPEC_CRUST_MANTLE) &
    stop 'NSPEC_CRUST_MANTLE_STRAIN_ONLY /= NSPEC_CRUST_MANTLE'
  ! checks movie type
  if (MOVIE_VOLUME_TYPE < 1 .or. MOVIE_VOLUME_TYPE > 9) &
    stop 'MOVIE_VOLUME_TYPE has to be in range from 1 to 9'

  allocate(mask_ibool(NGLOB_CRUST_MANTLE_3DMOVIE), &
           num_ibool_3dmovie(NGLOB_CRUST_MANTLE_3DMOVIE), &
           mask_3dmovie(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE), &
           muvstore_crust_mantle_3dmovie(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE),stat=ier)
  if (ier /= 0 ) stop 'Error allocating arrays muvstore_crust_mantle_3dmovie,..'

  ! counts total number of points for movie file output
  call movie_volume_count_points()

  allocate(nu_3dmovie(3,3,npoints_3dmovie),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating nu for 3D movie')

  call write_movie_volume_mesh(nu_3dmovie,num_ibool_3dmovie,mask_3dmovie,mask_ibool, &
                                          muvstore_crust_mantle_3dmovie,npoints_3dmovie)

  if (myrank == 0) then
    write(IMAIN,*) '  Movie volume:'
    write(IMAIN,*) '  Writing to movie3D*** files on local disk databases directory'
    select case (MOVIE_VOLUME_TYPE)
    case (1)
      write(IMAIN,*) '  movie output: strains'
    case (2)
      write(IMAIN,*) '  movie output: time integral of strains'
    case (3)
      write(IMAIN,*) '  movie output: potency or integral of strain'
    case (4)
      write(IMAIN,*) '  movie output: divergence and curl'
    case (5)
      write(IMAIN,*) '  movie output: displacement'
    case (6)
      write(IMAIN,*) '  movie output: velocity'
    case (7)
      write(IMAIN,*) '  movie output: norm of displacement'
    case (8)
      write(IMAIN,*) '  movie output: norm of velocity'
    case (9)
      write(IMAIN,*) '  movie output: norm of acceleration'
    case default
      call exit_MPI(myrank, 'MOVIE_VOLUME_TYPE has to be in range from 1 to 9')
    end select
    write(IMAIN,*)
    write(IMAIN,*) '  depth(T,B):',MOVIE_TOP,MOVIE_BOTTOM
    write(IMAIN,*) '  lon(W,E)  :',MOVIE_WEST,MOVIE_EAST
    write(IMAIN,*) '  lat(S,N)  :',MOVIE_SOUTH,MOVIE_NORTH
    write(IMAIN,*) '  Starting at time step:',MOVIE_START, 'ending at:',MOVIE_STOP,'every: ',NTSTEP_BETWEEN_FRAMES
    call flush_IMAIN()
  endif

  call synchronize_all()

  end subroutine prepare_movie_volume

