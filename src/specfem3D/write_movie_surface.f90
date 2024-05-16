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


  subroutine movie_surface_count_points()

  use specfem_par
  use specfem_par_crustmantle, only: NSPEC_TOP
  use specfem_par_movie, only: NIT,nmovie_points

  implicit none

  ! local parameters
  integer :: ispec2D,i,j,npoin

  ! gets number of points on surface mesh
  npoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        npoin = npoin + 1
      enddo
    enddo
  enddo

  ! checks
  if (npoin /= nmovie_points) then
    print *,'Error: movie points collected ',npoin,'not equal to calculated :',nmovie_points
    call exit_mpi(myrank,'Error confusing number of movie points')
  endif

  end subroutine movie_surface_count_points

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_surface_mesh()

! writes out movie point locations to file

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: store_val_x,store_val_y,store_val_z
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: store_val_x_all,store_val_y_all,store_val_z_all
  integer :: ipoin,ispec2D,ispec,i,j,k,ier,iglob,npoin
  character(len=MAX_STRING_LEN) :: outputname

  ! allocates movie surface arrays
  allocate(store_val_x(nmovie_points), &
           store_val_y(nmovie_points), &
           store_val_z(nmovie_points),stat=ier)
  if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface location arrays')

  ! allocates arrays for gathering movie point locations
  if (myrank == 0) then
    ! only main needs full arrays
    allocate(store_val_x_all(nmovie_points,0:NPROCTOT_VAL-1), &
             store_val_y_all(nmovie_points,0:NPROCTOT_VAL-1), &
             store_val_z_all(nmovie_points,0:NPROCTOT_VAL-1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface all arrays')
  else
    ! secondary processes only need dummy arrays
    allocate(store_val_x_all(1,1), &
             store_val_y_all(1,1), &
             store_val_z_all(1,1),stat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error allocating movie surface all arrays')
  endif

  ! gets coordinates of surface mesh
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
    ! the cubed sphere, the mapping changes for different chunks,
    ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
    ! for future consideration, like in create_movie_GMT_global.f90 ...
    k = NGLLZ
    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        ipoin = ipoin + 1
        ! stores values
        iglob = ibool_crust_mantle(i,j,k,ispec)
        store_val_x(ipoin) = rstore_crust_mantle(1,iglob) ! radius r (normalized)
        store_val_y(ipoin) = rstore_crust_mantle(2,iglob) ! colatitude theta (in radian)
        store_val_z(ipoin) = rstore_crust_mantle(3,iglob) ! longitude phi (in radian)
      enddo
    enddo
  enddo
  npoin = ipoin
  if (npoin /= nmovie_points ) call exit_mpi(myrank,'Error number of movie points not equal to nmovie_points')

  ! gather info on main proc
  call gather_all_cr(store_val_x,nmovie_points,store_val_x_all,nmovie_points,NPROCTOT_VAL)
  call gather_all_cr(store_val_y,nmovie_points,store_val_y_all,nmovie_points,NPROCTOT_VAL)
  call gather_all_cr(store_val_z,nmovie_points,store_val_z_all,nmovie_points,NPROCTOT_VAL)

  ! save movie data locations to disk in home directory
  if (myrank == 0) then

    ! outputs movie point locations to moviedata_xyz.bin file
    outputname = "/moviedata_xyz.bin"
    open(unit=IOUT,file=trim(OUTPUT_FILES)//trim(outputname), &
         status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening moviedata_xyz.bin file')

    ! point coordinates
    ! (given as r theta phi for geocentric coordinate system)
    write(IOUT) store_val_x_all
    write(IOUT) store_val_y_all
    write(IOUT) store_val_z_all
    close(IOUT)
  endif

  deallocate(store_val_x_all,store_val_y_all,store_val_z_all)
  deallocate(store_val_x,store_val_y,store_val_z)

  end subroutine write_movie_surface_mesh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine write_movie_surface()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_movie

  implicit none

  ! local parameters
  character(len=MAX_STRING_LEN) :: outputname
  integer :: ipoin,ispec2D,ispec,i,j,k,ier,iglob

  ! by default: save velocity here to avoid static offset on displacement for movies

  ! gets coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    ! in case of global, NCHUNKS_VAL == 6 simulations, be aware that for
    ! the cubed sphere, the mapping changes for different chunks,
    ! i.e. e.g. x(1,1) and x(5,5) flip left and right sides of the elements in geographical coordinates.
    ! for future consideration, like in create_movie_GMT_global.f90 ...
    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)

        ! wavefield values
        if (MOVIE_VOLUME_TYPE == 5) then
          ! stores displacement
          store_val_ux(ipoin) = displ_crust_mantle(1,iglob) * real(scale_displ,kind=CUSTOM_REAL)
          store_val_uy(ipoin) = displ_crust_mantle(2,iglob) * real(scale_displ,kind=CUSTOM_REAL)
          store_val_uz(ipoin) = displ_crust_mantle(3,iglob) * real(scale_displ,kind=CUSTOM_REAL)
        else
          ! stores velocity
          store_val_ux(ipoin) = veloc_crust_mantle(1,iglob) * real(scale_veloc,kind=CUSTOM_REAL)
          store_val_uy(ipoin) = veloc_crust_mantle(2,iglob) * real(scale_veloc,kind=CUSTOM_REAL)
          store_val_uz(ipoin) = veloc_crust_mantle(3,iglob) * real(scale_veloc,kind=CUSTOM_REAL)
        endif

      enddo
    enddo
  enddo

  ! gather info on main proc
  ! wavefield
  call gather_all_cr(store_val_ux,nmovie_points,store_val_ux_all,nmovie_points,NPROCTOT_VAL)
  call gather_all_cr(store_val_uy,nmovie_points,store_val_uy_all,nmovie_points,NPROCTOT_VAL)
  call gather_all_cr(store_val_uz,nmovie_points,store_val_uz_all,nmovie_points,NPROCTOT_VAL)

  ! save movie data to disk in home directory
  if (myrank == 0) then
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname, &
         status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening moviedata file')

    ! -> find movie point locations stored in file moviedata_xyz.bin

    ! wavefield values
    write(IOUT) store_val_ux_all
    write(IOUT) store_val_uy_all
    write(IOUT) store_val_uz_all

    close(IOUT)
  endif

  end subroutine write_movie_surface
