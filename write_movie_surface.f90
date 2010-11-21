!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            December 2010
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

  subroutine write_movie_surface(myrank,nmovie_points,scale_veloc,veloc_crust_mantle, &
                    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                    store_val_x,store_val_y,store_val_z, &
                    store_val_x_all,store_val_y_all,store_val_z_all, &
                    store_val_ux,store_val_uy,store_val_uz, &
                    store_val_ux_all,store_val_uy_all,store_val_uz_all, &
                    ibelm_top_crust_mantle,ibool_crust_mantle,nspec_top, &
                    NIT,it,OUTPUT_FILES)

  implicit none

  include 'mpif.h'
  include "precision.h"
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,nmovie_points
  double precision :: scale_veloc

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     veloc_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: &
      store_val_x,store_val_y,store_val_z, &
      store_val_ux,store_val_uy,store_val_uz

  real(kind=CUSTOM_REAL), dimension(nmovie_points,0:NPROCTOT_VAL-1) :: &
      store_val_x_all,store_val_y_all,store_val_z_all, &
      store_val_ux_all,store_val_uy_all,store_val_uz_all

  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer nspec_top,NIT,it
  character(len=150) OUTPUT_FILES

  ! local parameters
  character(len=150) :: outputname
  integer :: ipoin,ispec2D,ispec,i,j,k,ier,iglob

  ! save velocity here to avoid static offset on displacement for movies


  ! get coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
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
        store_val_x(ipoin) = xstore_crust_mantle(iglob)
        store_val_y(ipoin) = ystore_crust_mantle(iglob)
        store_val_z(ipoin) = zstore_crust_mantle(iglob)
        store_val_ux(ipoin) = veloc_crust_mantle(1,iglob)*scale_veloc
        store_val_uy(ipoin) = veloc_crust_mantle(2,iglob)*scale_veloc
        store_val_uz(ipoin) = veloc_crust_mantle(3,iglob)*scale_veloc
      enddo
    enddo

  enddo

  ! gather info on master proc
  ispec = nmovie_points
  call MPI_GATHER(store_val_x,ispec,CUSTOM_MPI_TYPE,store_val_x_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_y,ispec,CUSTOM_MPI_TYPE,store_val_y_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_z,ispec,CUSTOM_MPI_TYPE,store_val_z_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_ux,ispec,CUSTOM_MPI_TYPE,store_val_ux_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_uy,ispec,CUSTOM_MPI_TYPE,store_val_uy_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_uz,ispec,CUSTOM_MPI_TYPE,store_val_uz_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  ! save movie data to disk in home directory
  if(myrank == 0) then
    write(outputname,"('/moviedata',i6.6)") it
    open(unit=IOUT,file=trim(OUTPUT_FILES)//outputname,status='unknown',form='unformatted',action='write')
    write(IOUT) store_val_x_all
    write(IOUT) store_val_y_all
    write(IOUT) store_val_z_all
    write(IOUT) store_val_ux_all
    write(IOUT) store_val_uy_all
    write(IOUT) store_val_uz_all
    close(IOUT)
  endif

  end subroutine write_movie_surface
