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

!-------------------------------------------------------------------------------------------------
!
! compute the strain in the whole crust/mantle and inner core domains
!
!-------------------------------------------------------------------------------------------------

  subroutine compute_strain_att()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none
  ! local parameters
  integer :: ispec

  ! computes strain based on forward wavefield displ
  if (.not. GPU_MODE) then

    ! checks
    if (USE_DEVILLE_PRODUCTS_VAL) then

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_Dev(ispec,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                                            displ_inner_core,veloc_inner_core,0._CUSTOM_REAL, &
                                            ibool_inner_core, &
                                            hprime_xx,hprime_xxT, &
                                            deriv_mapping_inner_core, &
                                            epsilondev_xx_inner_core, &
                                            epsilondev_yy_inner_core, &
                                            epsilondev_xy_inner_core, &
                                            epsilondev_xz_inner_core, &
                                            epsilondev_yz_inner_core, &
                                            NSPEC_INNER_CORE_STRAIN_ONLY,eps_trace_over_3_inner_core)
      enddo
      ! crust mantle
      do ispec = 1, NSPEC_crust_mantle
        call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                            displ_crust_mantle,veloc_crust_mantle,0._CUSTOM_REAL, &
                                            ibool_crust_mantle, &
                                            hprime_xx,hprime_xxT, &
                                            deriv_mapping_crust_mantle, &
                                            epsilondev_xx_crust_mantle, &
                                            epsilondev_yy_crust_mantle, &
                                            epsilondev_xy_crust_mantle, &
                                            epsilondev_xz_crust_mantle, &
                                            epsilondev_yz_crust_mantle, &
                                            NSPEC_CRUST_MANTLE_STRAIN_ONLY,eps_trace_over_3_crust_mantle)
      enddo

    else

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_noDev(ispec,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                                              displ_inner_core,veloc_inner_core,0._CUSTOM_REAL, &
                                              ibool_inner_core, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                              etax_inner_core,etay_inner_core,etaz_inner_core, &
                                              gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                              epsilondev_xx_inner_core, &
                                              epsilondev_yy_inner_core, &
                                              epsilondev_xy_inner_core, &
                                              epsilondev_xz_inner_core, &
                                              epsilondev_yz_inner_core, &
                                              NSPEC_INNER_CORE_STRAIN_ONLY,eps_trace_over_3_inner_core)
      enddo
      ! crust mantle
      do ispec = 1, NSPEC_CRUST_MANTLE
        call compute_element_strain_att_noDev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                              displ_crust_mantle,veloc_crust_mantle,0._CUSTOM_REAL, &
                                              ibool_crust_mantle, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                              etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                              gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                              epsilondev_xx_crust_mantle, &
                                              epsilondev_yy_crust_mantle, &
                                              epsilondev_xy_crust_mantle, &
                                              epsilondev_xz_crust_mantle, &
                                              epsilondev_yz_crust_mantle, &
                                              NSPEC_CRUST_MANTLE_STRAIN_ONLY,eps_trace_over_3_crust_mantle)
      enddo
    endif

  else

    ! calculates strains on GPU
    ! note: deltat is zero, thus strain is computed based on < displ(:,:) > rather than < displ(:,:) + deltat * veloc(:,:) >
    !       nevertheless, we implement < displ(:,:) + deltat * veloc(:,:) > in order to have a more general calculation
    !       as done in the CPU routine as well
    call compute_strain_gpu(Mesh_pointer,0._CUSTOM_REAL,1)

  endif ! GPU_MODE

  end subroutine compute_strain_att

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_strain_att_backward()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore

  implicit none

  ! local parameters
  integer :: ispec

  ! computes strain based on backward/reconstructed wavefield b_displ
  if (.not. GPU_MODE) then

    ! checks
    if (USE_DEVILLE_PRODUCTS_VAL) then

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_Dev(ispec,NGLOB_INNER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                                            b_displ_inner_core,b_veloc_inner_core,0._CUSTOM_REAL, &
                                            ibool_inner_core, &
                                            hprime_xx,hprime_xxT, &
                                            deriv_mapping_inner_core, &
                                            b_epsilondev_xx_inner_core, &
                                            b_epsilondev_yy_inner_core, &
                                            b_epsilondev_xy_inner_core, &
                                            b_epsilondev_xz_inner_core, &
                                            b_epsilondev_yz_inner_core, &
                                            NSPEC_INNER_CORE_STRAIN_ONLY,b_eps_trace_over_3_inner_core)
      enddo

      ! crust mantle
      do ispec = 1, NSPEC_CRUST_MANTLE
        call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE_ADJOINT,NSPEC_CRUST_MANTLE_ADJOINT, &
                                            b_displ_crust_mantle,b_veloc_crust_mantle,0._CUSTOM_REAL, &
                                            ibool_crust_mantle, &
                                            hprime_xx,hprime_xxT, &
                                            deriv_mapping_crust_mantle, &
                                            b_epsilondev_xx_crust_mantle, &
                                            b_epsilondev_yy_crust_mantle, &
                                            b_epsilondev_xy_crust_mantle, &
                                            b_epsilondev_xz_crust_mantle, &
                                            b_epsilondev_yz_crust_mantle, &
                                            NSPEC_CRUST_MANTLE_STRAIN_ONLY,b_eps_trace_over_3_crust_mantle)
      enddo

    else

      ! inner core
      do ispec = 1, NSPEC_INNER_CORE
        call compute_element_strain_att_noDev(ispec,NGLOB_INNER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                                              b_displ_inner_core,b_veloc_inner_core,0._CUSTOM_REAL, &
                                              ibool_inner_core, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                              etax_inner_core,etay_inner_core,etaz_inner_core, &
                                              gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                              b_epsilondev_xx_inner_core, &
                                              b_epsilondev_yy_inner_core, &
                                              b_epsilondev_xy_inner_core, &
                                              b_epsilondev_xz_inner_core, &
                                              b_epsilondev_yz_inner_core, &
                                              NSPEC_INNER_CORE_STRAIN_ONLY,b_eps_trace_over_3_inner_core)
      enddo
      ! crust mantle
      do ispec = 1, NSPEC_crust_mantle
        call compute_element_strain_att_noDev(ispec,NGLOB_CRUST_MANTLE_ADJOINT,NSPEC_CRUST_MANTLE_ADJOINT, &
                                              b_displ_crust_mantle,b_veloc_crust_mantle,0._CUSTOM_REAL, &
                                              ibool_crust_mantle, &
                                              hprime_xx,hprime_yy,hprime_zz, &
                                              xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                              etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                              gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                              b_epsilondev_xx_crust_mantle, &
                                              b_epsilondev_yy_crust_mantle, &
                                              b_epsilondev_xy_crust_mantle, &
                                              b_epsilondev_xz_crust_mantle, &
                                              b_epsilondev_yz_crust_mantle, &
                                              NSPEC_CRUST_MANTLE_STRAIN_ONLY,b_eps_trace_over_3_crust_mantle)
      enddo
    endif

  else

    ! calculates strains on GPU
    call compute_strain_gpu(Mesh_pointer,0._CUSTOM_REAL,3)

  endif ! GPU_MODE

  end subroutine compute_strain_att_backward

