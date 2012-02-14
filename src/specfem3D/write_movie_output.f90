!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            April 2011
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

  subroutine write_movie_output()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie
  implicit none

  ! local parameters
!daniel: debugging
!  character(len=256) :: filename
!  logical, parameter :: SNAPSHOT_INNER_CORE = .true.

  ! save movie on surface
  if( MOVIE_SURFACE ) then
    if( mod(it,NTSTEP_BETWEEN_FRAMES) == 0) then
      ! save velocity here to avoid static offset on displacement for movies
      call write_movie_surface(myrank,nmovie_points,scale_veloc,veloc_crust_mantle, &
                    scale_displ,displ_crust_mantle, &
                    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                    store_val_x,store_val_y,store_val_z, &
                    store_val_x_all,store_val_y_all,store_val_z_all, &
                    store_val_ux,store_val_uy,store_val_uz, &
                    store_val_ux_all,store_val_uy_all,store_val_uz_all, &
                    ibelm_top_crust_mantle,ibool_crust_mantle, &
                    NSPEC2D_TOP(IREGION_CRUST_MANTLE), &
                    NIT,it,OUTPUT_FILES,MOVIE_VOLUME_TYPE)
    endif
  endif


  ! save movie in full 3D mesh
  if(MOVIE_VOLUME ) then
    if( mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0  &
      .and. it >= MOVIE_START .and. it <= MOVIE_STOP) then

      if (MOVIE_VOLUME_TYPE == 1) then  ! output strains

        call  write_movie_volume_strains(myrank,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,eps_trace_over_3_crust_mantle, &
                    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                    muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)

      else if (MOVIE_VOLUME_TYPE == 2 .or. MOVIE_VOLUME_TYPE == 3) then
        ! output the Time Integral of Strain, or \mu*TIS
        call  write_movie_volume_strains(myrank,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE,MOVIE_COARSE, &
                    it,Ieps_trace_over_3_crust_mantle,Iepsilondev_crust_mantle, &
                    muvstore_crust_mantle_3dmovie, &
                    mask_3dmovie,nu_3dmovie)

      else if (MOVIE_VOLUME_TYPE == 4) then ! output divergence and curl in whole volume

        call write_movie_volume_divcurl(myrank,it,eps_trace_over_3_crust_mantle,&
                        div_displ_outer_core, &
                        accel_outer_core,kappavstore_outer_core,rhostore_outer_core,ibool_outer_core, &
                        eps_trace_over_3_inner_core, &
                        epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                        epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                        epsilondev_xx_inner_core,epsilondev_yy_inner_core,epsilondev_xy_inner_core, &
                        epsilondev_xz_inner_core,epsilondev_yz_inner_core, &
                        LOCAL_PATH, &
                        displ_crust_mantle,displ_inner_core,displ_outer_core, &
                        veloc_crust_mantle,veloc_inner_core,veloc_outer_core, &
                        accel_crust_mantle,accel_inner_core, &
                        ibool_crust_mantle,ibool_inner_core)

      else if (MOVIE_VOLUME_TYPE == 5) then !output displacement
        scalingval = scale_displ
        call write_movie_volume_vector(myrank,it,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE, &
                    MOVIE_COARSE,ibool_crust_mantle,displ_crust_mantle, &
                    scalingval,mask_3dmovie,nu_3dmovie)

      else if (MOVIE_VOLUME_TYPE == 6) then !output velocity
        scalingval = scale_veloc
        call write_movie_volume_vector(myrank,it,npoints_3dmovie, &
                    LOCAL_PATH,MOVIE_VOLUME_TYPE, &
                    MOVIE_COARSE,ibool_crust_mantle,veloc_crust_mantle, &
                    scalingval,mask_3dmovie,nu_3dmovie)

      else

        call exit_MPI(myrank, 'MOVIE_VOLUME_TYPE has to be 1,2,3,4,5 or 6')

      endif ! MOVIE_VOLUME_TYPE
    endif
  endif ! MOVIE_VOLUME

!daniel: debugging
!  if( SNAPSHOT_INNER_CORE .and. mod(it-MOVIE_START,NTSTEP_BETWEEN_FRAMES) == 0  &
!      .and. it >= MOVIE_START .and. it <= MOVIE_STOP) then
!    ! VTK file output
!    ! displacement values
!    !write(prname,'(a,i6.6,a)') trim(LOCAL_PATH)//'/'//'proc',myrank,'_'
!    !write(filename,'(a,a,i6.6)') prname(1:len_trim(prname)),'reg_3_displ_',it
!    !call write_VTK_data_cr(idoubling_inner_core,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
!    !                    xstore_inner_core,ystore_inner_core,zstore_inner_core,ibool_inner_core, &
!    !                    displ_inner_core,filename)
!
!    write(prname,'(a)') 'OUTPUT_FILES/snapshot_all_'
!    write(filename,'(a,a,i6.6)') prname(1:len_trim(prname)),'reg_3_displ_',it
!    call write_VTK_data_cr_all(myrank,idoubling_inner_core, &
!                        NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
!                        xstore_inner_core,ystore_inner_core,zstore_inner_core,ibool_inner_core, &
!                        displ_inner_core,filename)
!
!  endif


  end subroutine write_movie_output
