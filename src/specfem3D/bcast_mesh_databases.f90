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


  subroutine bcast_mesh_databases_CM()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

! note: the size(..) function returns either integer(kind=4) or integer(kind=8)
!       depending on compiler flags (-mcmedium), thus adding a kind argument to have integer(kind=4) output

  !call bcast_all_i_for_database(NSPEC_CRUST_MANTLE, 1)
  !call bcast_all_i_for_database(NGLOB_CRUST_MANTLE, 1)
  !call bcast_all_i_for_database(NGLOB_XY_CM, 1)

  if (size(rho_vp_crust_mantle) > 0) then
    call bcast_all_cr_for_database(rho_vp_crust_mantle(1,1,1,1), size(rho_vp_crust_mantle,kind=4))
    call bcast_all_cr_for_database(rho_vs_crust_mantle(1,1,1,1), size(rho_vs_crust_mantle,kind=4))
  endif

  if (size(xstore_crust_mantle) > 0) then
    call bcast_all_cr_for_database(xstore_crust_mantle(1), size(xstore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(ystore_crust_mantle(1), size(ystore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(zstore_crust_mantle(1), size(zstore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(xix_crust_mantle(1,1,1,1), size(xix_crust_mantle,kind=4))
    call bcast_all_cr_for_database(xiy_crust_mantle(1,1,1,1), size(xiy_crust_mantle,kind=4))
    call bcast_all_cr_for_database(xiz_crust_mantle(1,1,1,1), size(xiz_crust_mantle,kind=4))
    call bcast_all_cr_for_database(etax_crust_mantle(1,1,1,1), size(etax_crust_mantle,kind=4))
    call bcast_all_cr_for_database(etay_crust_mantle(1,1,1,1), size(etay_crust_mantle,kind=4))
    call bcast_all_cr_for_database(etaz_crust_mantle(1,1,1,1), size(etaz_crust_mantle,kind=4))
    call bcast_all_cr_for_database(gammax_crust_mantle(1,1,1,1), size(gammax_crust_mantle,kind=4))
    call bcast_all_cr_for_database(gammay_crust_mantle(1,1,1,1), size(gammay_crust_mantle,kind=4))
    call bcast_all_cr_for_database(gammaz_crust_mantle(1,1,1,1), size(gammaz_crust_mantle,kind=4))
  endif

  if (size(rhostore_crust_mantle) > 0) then
    call bcast_all_cr_for_database(rhostore_crust_mantle(1,1,1,1), size(rhostore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(kappavstore_crust_mantle(1,1,1,1), size(kappavstore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(muvstore_crust_mantle(1,1,1,1), size(muvstore_crust_mantle,kind=4))
  endif

  if (size(kappahstore_crust_mantle) > 0) then
    call bcast_all_cr_for_database(kappahstore_crust_mantle(1,1,1,1), size(kappahstore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(muhstore_crust_mantle(1,1,1,1), size(muhstore_crust_mantle,kind=4))
    call bcast_all_cr_for_database(eta_anisostore_crust_mantle(1,1,1,1), size(eta_anisostore_crust_mantle,kind=4))
  endif

  if (size(c11store_crust_mantle) > 0) then
    call bcast_all_cr_for_database(c11store_crust_mantle(1,1,1,1), size(c11store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c12store_crust_mantle(1,1,1,1), size(c12store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c13store_crust_mantle(1,1,1,1), size(c13store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c14store_crust_mantle(1,1,1,1), size(c14store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c15store_crust_mantle(1,1,1,1), size(c15store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c16store_crust_mantle(1,1,1,1), size(c16store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c22store_crust_mantle(1,1,1,1), size(c22store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c23store_crust_mantle(1,1,1,1), size(c23store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c24store_crust_mantle(1,1,1,1), size(c24store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c25store_crust_mantle(1,1,1,1), size(c25store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c26store_crust_mantle(1,1,1,1), size(c26store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c33store_crust_mantle(1,1,1,1), size(c33store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c34store_crust_mantle(1,1,1,1), size(c34store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c35store_crust_mantle(1,1,1,1), size(c35store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c36store_crust_mantle(1,1,1,1), size(c36store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c44store_crust_mantle(1,1,1,1), size(c44store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c45store_crust_mantle(1,1,1,1), size(c45store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c46store_crust_mantle(1,1,1,1), size(c46store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c55store_crust_mantle(1,1,1,1), size(c55store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c56store_crust_mantle(1,1,1,1), size(c56store_crust_mantle,kind=4))
    call bcast_all_cr_for_database(c66store_crust_mantle(1,1,1,1), size(c66store_crust_mantle,kind=4))
  endif

  if (size(mu0store_crust_mantle) > 0) then
    call bcast_all_cr_for_database(mu0store_crust_mantle(1,1,1,1), size(mu0store_crust_mantle,kind=4))
  endif

  if (size(ibool_crust_mantle) > 0) then
    call bcast_all_i_for_database(ibool_crust_mantle(1,1,1,1), size(ibool_crust_mantle,kind=4))
    call bcast_all_l_for_database(ispec_is_tiso_crust_mantle(1), size(ispec_is_tiso_crust_mantle,kind=4))
  endif

  if (size(rmassx_crust_mantle) > 0) then
    call bcast_all_cr_for_database(rmassx_crust_mantle(1), size(rmassx_crust_mantle,kind=4))
    call bcast_all_cr_for_database(rmassy_crust_mantle(1), size(rmassy_crust_mantle,kind=4))
  endif

  if (size(rmassz_crust_mantle) > 0) then
    call bcast_all_cr_for_database(rmassz_crust_mantle(1), size(rmassz_crust_mantle,kind=4))
  endif

  !call bcast_all_i_for_database(NGLOB_CRUST_MANTLE_OCEANS, 1)
  if (size(rmass_ocean_load) > 0) then
    call bcast_all_cr_for_database(rmass_ocean_load(1), size(rmass_ocean_load,kind=4))
  endif

  if (size(b_rmassx_crust_mantle) > 0) then
    call bcast_all_cr_for_database(b_rmassx_crust_mantle(1), size(b_rmassx_crust_mantle,kind=4))
    call bcast_all_cr_for_database(b_rmassy_crust_mantle(1), size(b_rmassy_crust_mantle,kind=4))
  endif

  end subroutine bcast_mesh_databases_CM

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_OC()

  use specfem_par
  use specfem_par_outercore
  implicit none

  !call bcast_all_i_for_database(NSPEC_OUTER_CORE, 1)
  !call bcast_all_i_for_database(NGLOB_OUTER_CORE, 1)

  if (size(vp_outer_core) > 0) then
    call bcast_all_cr_for_database(vp_outer_core(1,1,1,1), size(vp_outer_core,kind=4))
  endif

  if (size(xstore_outer_core) > 0) then
    call bcast_all_cr_for_database(xstore_outer_core(1), size(xstore_outer_core,kind=4))
    call bcast_all_cr_for_database(ystore_outer_core(1), size(ystore_outer_core,kind=4))
    call bcast_all_cr_for_database(zstore_outer_core(1), size(zstore_outer_core,kind=4))
  endif

  if (size(xix_outer_core) > 0) then
    call bcast_all_cr_for_database(xix_outer_core(1,1,1,1), size(xix_outer_core,kind=4))
    call bcast_all_cr_for_database(xiy_outer_core(1,1,1,1), size(xiy_outer_core,kind=4))
    call bcast_all_cr_for_database(xiz_outer_core(1,1,1,1), size(xiz_outer_core,kind=4))
    call bcast_all_cr_for_database(etax_outer_core(1,1,1,1), size(etax_outer_core,kind=4))
    call bcast_all_cr_for_database(etay_outer_core(1,1,1,1), size(etay_outer_core,kind=4))
    call bcast_all_cr_for_database(etaz_outer_core(1,1,1,1), size(etaz_outer_core,kind=4))
    call bcast_all_cr_for_database(gammax_outer_core(1,1,1,1), size(gammax_outer_core,kind=4))
    call bcast_all_cr_for_database(gammay_outer_core(1,1,1,1), size(gammay_outer_core,kind=4))
    call bcast_all_cr_for_database(gammaz_outer_core(1,1,1,1), size(gammaz_outer_core,kind=4))
  endif

  if (size(rhostore_outer_core) > 0) then
    call bcast_all_cr_for_database(rhostore_outer_core(1,1,1,1), size(rhostore_outer_core,kind=4))
    call bcast_all_cr_for_database(kappavstore_outer_core(1,1,1,1), size(kappavstore_outer_core,kind=4))
  endif

  if (size(ibool_outer_core) > 0) then
    call bcast_all_i_for_database(ibool_outer_core(1,1,1,1), size(ibool_outer_core,kind=4))
  endif

  if (size(rmass_outer_core) > 0) then
    call bcast_all_cr_for_database(rmass_outer_core(1), size(rmass_outer_core,kind=4))
  endif

  end subroutine bcast_mesh_databases_OC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_IC()

  use specfem_par
  use specfem_par_innercore
  implicit none

  !call bcast_all_i_for_database(NSPEC_INNER_CORE, 1)
  !call bcast_all_i_for_database(NGLOB_INNER_CORE, 1)
  !call bcast_all_i_for_database(NGLOB_XY_IC, 1)

  if (size(xstore_inner_core) > 0) then
    call bcast_all_cr_for_database(xstore_inner_core(1), size(xstore_inner_core,kind=4))
    call bcast_all_cr_for_database(ystore_inner_core(1), size(ystore_inner_core,kind=4))
    call bcast_all_cr_for_database(zstore_inner_core(1), size(zstore_inner_core,kind=4))
  endif

  if (size(xix_inner_core) > 0) then
    call bcast_all_cr_for_database(xix_inner_core(1,1,1,1), size(xix_inner_core,kind=4))
    call bcast_all_cr_for_database(xiy_inner_core(1,1,1,1), size(xiy_inner_core,kind=4))
    call bcast_all_cr_for_database(xiz_inner_core(1,1,1,1), size(xiz_inner_core,kind=4))
    call bcast_all_cr_for_database(etax_inner_core(1,1,1,1), size(etax_inner_core,kind=4))
    call bcast_all_cr_for_database(etay_inner_core(1,1,1,1), size(etay_inner_core,kind=4))
    call bcast_all_cr_for_database(etaz_inner_core(1,1,1,1), size(etaz_inner_core,kind=4))
    call bcast_all_cr_for_database(gammax_inner_core(1,1,1,1), size(gammax_inner_core,kind=4))
    call bcast_all_cr_for_database(gammay_inner_core(1,1,1,1), size(gammay_inner_core,kind=4))
    call bcast_all_cr_for_database(gammaz_inner_core(1,1,1,1), size(gammaz_inner_core,kind=4))
  endif

  if (size(rhostore_inner_core) > 0) then
    call bcast_all_cr_for_database(rhostore_inner_core(1,1,1,1), size(rhostore_inner_core,kind=4))
    call bcast_all_cr_for_database(kappavstore_inner_core(1,1,1,1), size(kappavstore_inner_core,kind=4))
    call bcast_all_cr_for_database(muvstore_inner_core(1,1,1,1), size(muvstore_inner_core,kind=4))
  endif

  if (size(c11store_inner_core) > 0) then
    call bcast_all_cr_for_database(c11store_inner_core(1,1,1,1), size(c11store_inner_core,kind=4))
    call bcast_all_cr_for_database(c12store_inner_core(1,1,1,1), size(c12store_inner_core,kind=4))
    call bcast_all_cr_for_database(c13store_inner_core(1,1,1,1), size(c13store_inner_core,kind=4))
    call bcast_all_cr_for_database(c33store_inner_core(1,1,1,1), size(c33store_inner_core,kind=4))
    call bcast_all_cr_for_database(c44store_inner_core(1,1,1,1), size(c44store_inner_core,kind=4))
  endif

  if (size(ibool_inner_core) > 0) then
    call bcast_all_i_for_database(ibool_inner_core(1,1,1,1), size(ibool_inner_core,kind=4))
    call bcast_all_i_for_database(idoubling_inner_core(1), size(idoubling_inner_core,kind=4))
  endif

  if (size(rmassx_inner_core) > 0) then
    call bcast_all_cr_for_database(rmassx_inner_core(1), size(rmassx_inner_core,kind=4))
    call bcast_all_cr_for_database(rmassy_inner_core(1), size(rmassy_inner_core,kind=4))
  endif

  if (size(rmassz_inner_core) > 0) then
    call bcast_all_cr_for_database(rmassz_inner_core(1), size(rmassz_inner_core,kind=4))
  endif

  if (size(b_rmassx_inner_core) > 0) then
    call bcast_all_cr_for_database(b_rmassx_inner_core(1), size(b_rmassx_inner_core,kind=4))
    call bcast_all_cr_for_database(b_rmassy_inner_core(1), size(b_rmassy_inner_core,kind=4))
  endif

  end subroutine bcast_mesh_databases_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_TRINF()

  use specfem_par
  use specfem_par_trinfinite
  implicit none

  if (size(xstore_trinfinite) > 0) then
    call bcast_all_cr_for_database(xstore_trinfinite(1), size(xstore_trinfinite,kind=4))
    call bcast_all_cr_for_database(ystore_trinfinite(1), size(ystore_trinfinite,kind=4))
    call bcast_all_cr_for_database(zstore_trinfinite(1), size(zstore_trinfinite,kind=4))
  endif

  if (size(xix_trinfinite) > 0) then
    call bcast_all_cr_for_database(xix_trinfinite(1,1,1,1), size(xix_trinfinite,kind=4))
    call bcast_all_cr_for_database(xiy_trinfinite(1,1,1,1), size(xiy_trinfinite,kind=4))
    call bcast_all_cr_for_database(xiz_trinfinite(1,1,1,1), size(xiz_trinfinite,kind=4))
    call bcast_all_cr_for_database(etax_trinfinite(1,1,1,1), size(etax_trinfinite,kind=4))
    call bcast_all_cr_for_database(etay_trinfinite(1,1,1,1), size(etay_trinfinite,kind=4))
    call bcast_all_cr_for_database(etaz_trinfinite(1,1,1,1), size(etaz_trinfinite,kind=4))
    call bcast_all_cr_for_database(gammax_trinfinite(1,1,1,1), size(gammax_trinfinite,kind=4))
    call bcast_all_cr_for_database(gammay_trinfinite(1,1,1,1), size(gammay_trinfinite,kind=4))
    call bcast_all_cr_for_database(gammaz_trinfinite(1,1,1,1), size(gammaz_trinfinite,kind=4))
  endif

  if (size(ibool_trinfinite) > 0) then
    call bcast_all_i_for_database(ibool_trinfinite(1,1,1,1), size(ibool_trinfinite,kind=4))
  endif

  end subroutine bcast_mesh_databases_TRINF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_INF()

  use specfem_par
  use specfem_par_infinite
  implicit none

  if (size(xstore_infinite) > 0) then
    call bcast_all_cr_for_database(xstore_infinite(1), size(xstore_infinite,kind=4))
    call bcast_all_cr_for_database(ystore_infinite(1), size(ystore_infinite,kind=4))
    call bcast_all_cr_for_database(zstore_infinite(1), size(zstore_infinite,kind=4))
  endif

  if (size(xix_infinite) > 0) then
    call bcast_all_cr_for_database(xix_infinite(1,1,1,1), size(xix_infinite,kind=4))
    call bcast_all_cr_for_database(xiy_infinite(1,1,1,1), size(xiy_infinite,kind=4))
    call bcast_all_cr_for_database(xiz_infinite(1,1,1,1), size(xiz_infinite,kind=4))
    call bcast_all_cr_for_database(etax_infinite(1,1,1,1), size(etax_infinite,kind=4))
    call bcast_all_cr_for_database(etay_infinite(1,1,1,1), size(etay_infinite,kind=4))
    call bcast_all_cr_for_database(etaz_infinite(1,1,1,1), size(etaz_infinite,kind=4))
    call bcast_all_cr_for_database(gammax_infinite(1,1,1,1), size(gammax_infinite,kind=4))
    call bcast_all_cr_for_database(gammay_infinite(1,1,1,1), size(gammay_infinite,kind=4))
    call bcast_all_cr_for_database(gammaz_infinite(1,1,1,1), size(gammaz_infinite,kind=4))
  endif

  if (size(ibool_infinite) > 0) then
    call bcast_all_i_for_database(ibool_infinite(1,1,1,1), size(ibool_infinite,kind=4))
  endif

  end subroutine bcast_mesh_databases_INF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_coupling()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  use specfem_par_trinfinite
  use specfem_par_infinite

  implicit none

  ! crust/mantle
  call bcast_all_i_for_database(nspec2D_xmin_crust_mantle, 1)
  call bcast_all_i_for_database(nspec2D_xmax_crust_mantle, 1)
  call bcast_all_i_for_database(nspec2D_ymin_crust_mantle, 1)
  call bcast_all_i_for_database(nspec2D_ymax_crust_mantle, 1)
  !call bcast_all_i_for_database(nspec2D_zmin_crust_mantle, 1)

  if (size(ibelm_xmin_crust_mantle) > 0) then
    call bcast_all_i_for_database(ibelm_xmin_crust_mantle(1), size(ibelm_xmin_crust_mantle,kind=4))
    call bcast_all_i_for_database(ibelm_xmax_crust_mantle(1), size(ibelm_xmax_crust_mantle,kind=4))
    call bcast_all_i_for_database(ibelm_ymin_crust_mantle(1), size(ibelm_ymin_crust_mantle,kind=4))
    call bcast_all_i_for_database(ibelm_ymax_crust_mantle(1), size(ibelm_ymax_crust_mantle,kind=4))
    call bcast_all_i_for_database(ibelm_bottom_crust_mantle(1), size(ibelm_bottom_crust_mantle,kind=4))
    call bcast_all_i_for_database(ibelm_top_crust_mantle(1), size(ibelm_top_crust_mantle,kind=4))
  endif

  if (size(normal_xmin_crust_mantle) > 0) then
    call bcast_all_cr_for_database(normal_xmin_crust_mantle(1,1,1,1), size(normal_xmin_crust_mantle,kind=4))
    call bcast_all_cr_for_database(normal_xmax_crust_mantle(1,1,1,1), size(normal_xmax_crust_mantle,kind=4))
    call bcast_all_cr_for_database(normal_ymin_crust_mantle(1,1,1,1), size(normal_ymin_crust_mantle,kind=4))
    call bcast_all_cr_for_database(normal_ymax_crust_mantle(1,1,1,1), size(normal_ymax_crust_mantle,kind=4))
    call bcast_all_cr_for_database(normal_bottom_crust_mantle(1,1,1,1), size(normal_bottom_crust_mantle,kind=4))
    call bcast_all_cr_for_database(normal_top_crust_mantle(1,1,1,1), size(normal_top_crust_mantle,kind=4))

    call bcast_all_cr_for_database(jacobian2D_xmin_crust_mantle(1,1,1), size(jacobian2D_xmin_crust_mantle,kind=4))
    call bcast_all_cr_for_database(jacobian2D_xmax_crust_mantle(1,1,1), size(jacobian2D_xmax_crust_mantle,kind=4))
    call bcast_all_cr_for_database(jacobian2D_ymin_crust_mantle(1,1,1), size(jacobian2D_ymin_crust_mantle,kind=4))
    call bcast_all_cr_for_database(jacobian2D_ymax_crust_mantle(1,1,1), size(jacobian2D_ymax_crust_mantle,kind=4))
    call bcast_all_cr_for_database(jacobian2D_bottom_crust_mantle(1,1,1), size(jacobian2D_bottom_crust_mantle,kind=4))
    call bcast_all_cr_for_database(jacobian2D_top_crust_mantle(1,1,1), size(jacobian2D_top_crust_mantle,kind=4))
  endif

  ! outer core
  call bcast_all_i_for_database(nspec2D_xmin_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_xmax_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_ymin_outer_core, 1)
  call bcast_all_i_for_database(nspec2D_ymax_outer_core, 1)

  if (size(ibelm_xmin_outer_core) > 0) then
    call bcast_all_i_for_database(ibelm_xmin_outer_core(1), size(ibelm_xmin_outer_core,kind=4))
    call bcast_all_i_for_database(ibelm_xmax_outer_core(1), size(ibelm_xmax_outer_core,kind=4))
    call bcast_all_i_for_database(ibelm_ymin_outer_core(1), size(ibelm_ymin_outer_core,kind=4))
    call bcast_all_i_for_database(ibelm_ymax_outer_core(1), size(ibelm_ymax_outer_core,kind=4))
    call bcast_all_i_for_database(ibelm_bottom_outer_core(1), size(ibelm_bottom_outer_core,kind=4))
    call bcast_all_i_for_database(ibelm_top_outer_core(1), size(ibelm_top_outer_core,kind=4))
  endif

  if (size(normal_xmin_outer_core) > 0) then
    call bcast_all_cr_for_database(normal_xmin_outer_core(1,1,1,1), size(normal_xmin_outer_core,kind=4))
    call bcast_all_cr_for_database(normal_xmax_outer_core(1,1,1,1), size(normal_xmax_outer_core,kind=4))
    call bcast_all_cr_for_database(normal_ymin_outer_core(1,1,1,1), size(normal_ymin_outer_core,kind=4))
    call bcast_all_cr_for_database(normal_ymax_outer_core(1,1,1,1), size(normal_ymax_outer_core,kind=4))
    call bcast_all_cr_for_database(normal_bottom_outer_core(1,1,1,1), size(normal_bottom_outer_core,kind=4))
    call bcast_all_cr_for_database(normal_top_outer_core(1,1,1,1), size(normal_top_outer_core,kind=4))

    call bcast_all_cr_for_database(jacobian2D_xmin_outer_core(1,1,1), size(jacobian2D_xmin_outer_core,kind=4))
    call bcast_all_cr_for_database(jacobian2D_xmax_outer_core(1,1,1), size(jacobian2D_xmax_outer_core,kind=4))
    call bcast_all_cr_for_database(jacobian2D_ymin_outer_core(1,1,1), size(jacobian2D_ymin_outer_core,kind=4))
    call bcast_all_cr_for_database(jacobian2D_ymax_outer_core(1,1,1), size(jacobian2D_ymax_outer_core,kind=4))
    call bcast_all_cr_for_database(jacobian2D_bottom_outer_core(1,1,1), size(jacobian2D_bottom_outer_core,kind=4))
    call bcast_all_cr_for_database(jacobian2D_top_outer_core(1,1,1), size(jacobian2D_top_outer_core,kind=4))
  endif

  ! inner core
  call bcast_all_i_for_database(nspec2D_xmin_inner_core, 1)
  call bcast_all_i_for_database(nspec2D_xmax_inner_core, 1)
  call bcast_all_i_for_database(nspec2D_ymin_inner_core, 1)
  call bcast_all_i_for_database(nspec2D_ymax_inner_core, 1)

  ! boundary parameters
  if (size(ibelm_xmin_inner_core) > 0) then
    call bcast_all_i_for_database(ibelm_xmin_inner_core(1), size(ibelm_xmin_inner_core,kind=4))
    call bcast_all_i_for_database(ibelm_xmax_inner_core(1), size(ibelm_xmax_inner_core,kind=4))
    call bcast_all_i_for_database(ibelm_ymin_inner_core(1), size(ibelm_ymin_inner_core,kind=4))
    call bcast_all_i_for_database(ibelm_ymax_inner_core(1), size(ibelm_ymax_inner_core,kind=4))
    call bcast_all_i_for_database(ibelm_bottom_inner_core(1), size(ibelm_bottom_inner_core,kind=4))
    call bcast_all_i_for_database(ibelm_top_inner_core(1), size(ibelm_top_inner_core,kind=4))
  endif

  if (FULL_GRAVITY) then
    if (ADD_TRINF) then
      ! transition-to-infinite
      call bcast_all_i_for_database(nspec2D_xmin_trinfinite, 1)
      call bcast_all_i_for_database(nspec2D_xmax_trinfinite, 1)
      call bcast_all_i_for_database(nspec2D_ymin_trinfinite, 1)
      call bcast_all_i_for_database(nspec2D_ymax_trinfinite, 1)
      if (size(ibelm_xmin_trinfinite) > 0) then
        call bcast_all_i_for_database(ibelm_xmin_trinfinite(1), size(ibelm_xmin_trinfinite,kind=4))
        call bcast_all_i_for_database(ibelm_xmax_trinfinite(1), size(ibelm_xmax_trinfinite,kind=4))
        call bcast_all_i_for_database(ibelm_ymin_trinfinite(1), size(ibelm_ymin_trinfinite,kind=4))
        call bcast_all_i_for_database(ibelm_ymax_trinfinite(1), size(ibelm_ymax_trinfinite,kind=4))
        call bcast_all_i_for_database(ibelm_bottom_trinfinite(1), size(ibelm_bottom_trinfinite,kind=4))
        call bcast_all_i_for_database(ibelm_top_trinfinite(1), size(ibelm_top_trinfinite,kind=4))
      endif
    endif
    ! infinite
    call bcast_all_i_for_database(nspec2D_xmin_infinite, 1)
    call bcast_all_i_for_database(nspec2D_xmax_infinite, 1)
    call bcast_all_i_for_database(nspec2D_ymin_infinite, 1)
    call bcast_all_i_for_database(nspec2D_ymax_infinite, 1)
    if (size(ibelm_xmin_infinite) > 0) then
      call bcast_all_i_for_database(ibelm_xmin_infinite(1), size(ibelm_xmin_infinite,kind=4))
      call bcast_all_i_for_database(ibelm_xmax_infinite(1), size(ibelm_xmax_infinite,kind=4))
      call bcast_all_i_for_database(ibelm_ymin_infinite(1), size(ibelm_ymin_infinite,kind=4))
      call bcast_all_i_for_database(ibelm_ymax_infinite(1), size(ibelm_ymax_infinite,kind=4))
      call bcast_all_i_for_database(ibelm_bottom_infinite(1), size(ibelm_bottom_infinite,kind=4))
      call bcast_all_i_for_database(ibelm_top_infinite(1), size(ibelm_top_infinite,kind=4))
    endif
  endif

  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then
    if (size(ibelm_moho_top) > 0) then
      call bcast_all_i_for_database(ibelm_moho_top(1), size(ibelm_moho_top,kind=4))
      call bcast_all_i_for_database(ibelm_moho_bot(1), size(ibelm_moho_bot,kind=4))
      call bcast_all_cr_for_database(normal_moho(1,1,1,1), size(normal_moho,kind=4))
    endif
    if (size(ibelm_400_top) > 0) then
      call bcast_all_i_for_database(ibelm_400_top(1), size(ibelm_400_top,kind=4))
      call bcast_all_i_for_database(ibelm_400_bot(1), size(ibelm_400_bot,kind=4))
      call bcast_all_cr_for_database(normal_400(1,1,1,1), size(normal_400,kind=4))
    endif
    if (size(ibelm_670_top) > 0) then
      call bcast_all_i_for_database(ibelm_670_top(1), size(ibelm_670_top,kind=4))
      call bcast_all_i_for_database(ibelm_670_bot(1), size(ibelm_670_bot,kind=4))
      call bcast_all_cr_for_database(normal_670(1,1,1,1), size(normal_670,kind=4))
    endif
  endif

  end subroutine bcast_mesh_databases_coupling

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_CM()

  use specfem_par
  use specfem_par_crustmantle
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_crust_mantle, 1)

  ! could also test for not allocated, only reader processes have
  ! allocated these arrays.
  if (.not. I_should_read_the_database) then
    allocate(my_neighbors_crust_mantle(num_interfaces_crust_mantle), &
             nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_crust_mantle etc.')
    my_neighbors_crust_mantle(:) = 0; nibool_interfaces_crust_mantle(:) = 0
  endif
  if (num_interfaces_crust_mantle > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_cm, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
      ibool_interfaces_crust_mantle(:,:) = 0
    endif

    call bcast_all_i_for_database(my_neighbors_crust_mantle(1), size(my_neighbors_crust_mantle,kind=4))
    call bcast_all_i_for_database(nibool_interfaces_crust_mantle(1), size(nibool_interfaces_crust_mantle,kind=4))
    call bcast_all_i_for_database(ibool_interfaces_crust_mantle(1,1), size(ibool_interfaces_crust_mantle,kind=4))
  else
    ! dummy array
    max_nibool_interfaces_cm = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_crust_mantle')
      ibool_interfaces_crust_mantle(:,:) = 0
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_crust_mantle, 1)
  call bcast_all_i_for_database(nspec_outer_crust_mantle, 1)
  call bcast_all_i_for_database(num_phase_ispec_crust_mantle, 1)
  if (num_phase_ispec_crust_mantle < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_crust_mantle is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_crust_mantle(num_phase_ispec_crust_mantle,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_crust_mantle')
    phase_ispec_inner_crust_mantle(:,:) = 0
  endif

  if (num_phase_ispec_crust_mantle > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_crust_mantle(1,1), size(phase_ispec_inner_crust_mantle,kind=4))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_crust_mantle, 1)
    call bcast_all_i_for_database(num_colors_inner_crust_mantle, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
      num_elem_colors_crust_mantle(:) = 0
    endif

    call bcast_all_i_for_database(num_elem_colors_crust_mantle(1), size(num_elem_colors_crust_mantle,kind=4))
  else
    ! allocates dummy arrays
    num_colors_outer_crust_mantle = 0
    num_colors_inner_crust_mantle = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_crust_mantle(num_colors_outer_crust_mantle + num_colors_inner_crust_mantle), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_crust_mantle array')
      num_elem_colors_crust_mantle(:) = 0
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_CM

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_OC()

  use specfem_par
  use specfem_par_outercore
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_outer_core, 1)
  if (.not. I_should_read_the_database) then
    allocate(my_neighbors_outer_core(num_interfaces_outer_core), &
             nibool_interfaces_outer_core(num_interfaces_outer_core), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_outer_core etc.')
    my_neighbors_outer_core(:) = 0; nibool_interfaces_outer_core(:) = 0
  endif

  if (num_interfaces_outer_core > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_oc, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_outer_core')
      ibool_interfaces_outer_core(:,:) = 0
    endif
  endif

  if (num_interfaces_outer_core > 0) then
    call bcast_all_i_for_database(my_neighbors_outer_core(1), size(my_neighbors_outer_core,kind=4))
    call bcast_all_i_for_database(nibool_interfaces_outer_core(1), size(nibool_interfaces_outer_core,kind=4))
    call bcast_all_i_for_database(ibool_interfaces_outer_core(1,1), size(ibool_interfaces_outer_core,kind=4))
  else
    ! dummy array
    max_nibool_interfaces_oc = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_outer_core(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_outer_core')
      ibool_interfaces_outer_core(:,:) = 0
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_outer_core, 1)
  call bcast_all_i_for_database(nspec_outer_outer_core, 1)
  call bcast_all_i_for_database(num_phase_ispec_outer_core, 1)
  if (num_phase_ispec_outer_core < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_outer_core is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_outer_core(num_phase_ispec_outer_core,2),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_outer_core')
    phase_ispec_inner_outer_core(:,:) = 0
  endif

  if (num_phase_ispec_outer_core > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_outer_core(1,1), size(phase_ispec_inner_outer_core,kind=4))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_outer_core, 1)
    call bcast_all_i_for_database(num_colors_inner_outer_core, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
      num_elem_colors_outer_core(:) = 0
    endif

    call bcast_all_i_for_database(num_elem_colors_outer_core(1), size(num_elem_colors_outer_core,kind=4))
  else
    ! allocates dummy arrays
    num_colors_outer_outer_core = 0
    num_colors_inner_outer_core = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_outer_core(num_colors_outer_outer_core + num_colors_inner_outer_core), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_outer_core array')
      num_elem_colors_outer_core(:) = 0
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_OC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_IC()

  use specfem_par
  use specfem_par_innercore
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_inner_core, 1)
  if (.not. I_should_read_the_database) then
    allocate(my_neighbors_inner_core(num_interfaces_inner_core), &
             nibool_interfaces_inner_core(num_interfaces_inner_core), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_inner_core etc.')
    my_neighbors_inner_core(:) = 0; nibool_interfaces_inner_core(:) = 0
  endif

  if (num_interfaces_inner_core > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_ic, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic,num_interfaces_inner_core), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_inner_core')
      ibool_interfaces_inner_core(:,:) = 0
    endif

    call bcast_all_i_for_database(my_neighbors_inner_core(1), size(my_neighbors_inner_core,kind=4))
    call bcast_all_i_for_database(nibool_interfaces_inner_core(1), size(nibool_interfaces_inner_core,kind=4))
    call bcast_all_i_for_database(ibool_interfaces_inner_core(1,1), size(ibool_interfaces_inner_core,kind=4))
  else
    ! dummy array
    max_nibool_interfaces_ic = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_inner_core(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_inner_core')
      ibool_interfaces_inner_core(:,:) = 0
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_inner_core, 1)
  call bcast_all_i_for_database(nspec_outer_inner_core, 1)
  call bcast_all_i_for_database(num_phase_ispec_inner_core, 1)
  if (num_phase_ispec_inner_core < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_inner_core is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_inner_core(num_phase_ispec_inner_core,2), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_inner_core')
    phase_ispec_inner_inner_core(:,:) = 0
  endif

  if (num_phase_ispec_inner_core > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_inner_core(1,1), size(phase_ispec_inner_inner_core,kind=4))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_inner_core, 1)
    call bcast_all_i_for_database(num_colors_inner_inner_core, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
      num_elem_colors_inner_core(:) = 0
    endif

    call bcast_all_i_for_database(num_elem_colors_inner_core(1), size(num_elem_colors_inner_core,kind=4))
  else
    ! allocates dummy arrays
    num_colors_outer_inner_core = 0
    num_colors_inner_inner_core = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_inner_core(num_colors_outer_inner_core + num_colors_inner_inner_core), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_inner_core array')
      num_elem_colors_inner_core(:) = 0
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_IC

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_TRINF()

  use specfem_par
  use specfem_par_trinfinite
  use specfem_par_full_gravity
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_trinfinite, 1)
  if (.not. I_should_read_the_database) then
    allocate(my_neighbors_trinfinite(num_interfaces_trinfinite), &
             nibool_interfaces_trinfinite(num_interfaces_trinfinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_trinfinite etc.')
    my_neighbors_trinfinite(:) = 0; nibool_interfaces_trinfinite(:) = 0
  endif

  if (num_interfaces_trinfinite > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_trinfinite, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_trinfinite(max_nibool_interfaces_trinfinite,num_interfaces_trinfinite), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_trinfinite')
      ibool_interfaces_trinfinite(:,:) = 0
    endif

    call bcast_all_i_for_database(my_neighbors_trinfinite(1), size(my_neighbors_trinfinite,kind=4))
    call bcast_all_i_for_database(nibool_interfaces_trinfinite(1), size(nibool_interfaces_trinfinite,kind=4))
    call bcast_all_i_for_database(ibool_interfaces_trinfinite(1,1), size(ibool_interfaces_trinfinite,kind=4))
  else
    ! dummy array
    max_nibool_interfaces_trinfinite = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_trinfinite(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_trinfinite')
      ibool_interfaces_trinfinite(:,:) = 0
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_trinfinite, 1)
  call bcast_all_i_for_database(nspec_outer_trinfinite, 1)
  call bcast_all_i_for_database(num_phase_ispec_trinfinite, 1)
  if (num_phase_ispec_trinfinite < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_trinfinite is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_trinfinite(num_phase_ispec_trinfinite,2), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_trinfinite')
    phase_ispec_inner_trinfinite(:,:) = 0
  endif

  if (num_phase_ispec_trinfinite > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_trinfinite(1,1), size(phase_ispec_inner_trinfinite,kind=4))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_trinfinite, 1)
    call bcast_all_i_for_database(num_colors_inner_trinfinite, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_trinfinite(num_colors_outer_trinfinite + num_colors_inner_trinfinite), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_trinfinite array')
      num_elem_colors_trinfinite(:) = 0
    endif

    call bcast_all_i_for_database(num_elem_colors_trinfinite(1), size(num_elem_colors_trinfinite,kind=4))
  else
    ! allocates dummy arrays
    num_colors_outer_trinfinite = 0
    num_colors_inner_trinfinite = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_trinfinite(num_colors_outer_trinfinite + num_colors_inner_trinfinite), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_trinfinite array')
      num_elem_colors_trinfinite(:) = 0
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_TRINF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_MPI_INF()

  use specfem_par
  use specfem_par_infinite
  use specfem_par_full_gravity
  implicit none

  ! local parameters
  integer :: ier

  ! MPI interfaces
  call bcast_all_i_for_database(num_interfaces_infinite, 1)
  if (.not. I_should_read_the_database) then
    allocate(my_neighbors_infinite(num_interfaces_infinite), &
             nibool_interfaces_infinite(num_interfaces_infinite), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbors_infinite etc.')
    my_neighbors_infinite(:) = 0; nibool_interfaces_infinite(:) = 0
  endif

  if (num_interfaces_infinite > 0) then
    call bcast_all_i_for_database(max_nibool_interfaces_infinite, 1)
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_infinite(max_nibool_interfaces_infinite,num_interfaces_infinite), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_infinite')
      ibool_interfaces_infinite(:,:) = 0
    endif

    call bcast_all_i_for_database(my_neighbors_infinite(1), size(my_neighbors_infinite,kind=4))
    call bcast_all_i_for_database(nibool_interfaces_infinite(1), size(nibool_interfaces_infinite,kind=4))
    call bcast_all_i_for_database(ibool_interfaces_infinite(1,1), size(ibool_interfaces_infinite,kind=4))
  else
    ! dummy array
    max_nibool_interfaces_infinite = 0
    if (.not. I_should_read_the_database) then
      allocate(ibool_interfaces_infinite(0,0),stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array dummy ibool_interfaces_infinite')
      ibool_interfaces_infinite(:,:) = 0
    endif
  endif

  ! inner / outer elements
  call bcast_all_i_for_database(nspec_inner_infinite, 1)
  call bcast_all_i_for_database(nspec_outer_infinite, 1)
  call bcast_all_i_for_database(num_phase_ispec_infinite, 1)
  if (num_phase_ispec_infinite < 0 ) &
    call exit_mpi(myrank,'Error num_phase_ispec_infinite is < zero')

  if (.not. I_should_read_the_database) then
    allocate(phase_ispec_inner_infinite(num_phase_ispec_infinite,2), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array phase_ispec_inner_infinite')
    phase_ispec_inner_infinite(:,:) = 0
  endif

  if (num_phase_ispec_infinite > 0) then
    call bcast_all_i_for_database(phase_ispec_inner_infinite(1,1), size(phase_ispec_inner_infinite,kind=4))
  endif

  ! mesh coloring for GPUs
  if (USE_MESH_COLORING_GPU) then
    ! colors
    call bcast_all_i_for_database(num_colors_outer_infinite, 1)
    call bcast_all_i_for_database(num_colors_inner_infinite, 1)

    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_infinite(num_colors_outer_infinite + num_colors_inner_infinite), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_infinite array')
      num_elem_colors_infinite(:) = 0
    endif

    call bcast_all_i_for_database(num_elem_colors_infinite(1), size(num_elem_colors_infinite,kind=4))
  else
    ! allocates dummy arrays
    num_colors_outer_infinite = 0
    num_colors_inner_infinite = 0
    if (.not. I_should_read_the_database) then
      allocate(num_elem_colors_infinite(num_colors_outer_infinite + num_colors_inner_infinite), &
               stat=ier)
      if (ier /= 0 ) call exit_mpi(myrank,'Error allocating num_elem_colors_infinite array')
      num_elem_colors_infinite(:) = 0
    endif
  endif

  end subroutine bcast_mesh_databases_MPI_INF

!
!-------------------------------------------------------------------------------------------------
!

  subroutine bcast_mesh_databases_stacey()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ier

  ! crust and mantle
  if (NSPEC_CRUST_MANTLE > 0) then
    call bcast_all_i_for_database(num_abs_boundary_faces_crust_mantle,1)

    if (.not. I_should_read_the_database) then
      ! allocates absorbing boundary arrays
      if (num_abs_boundary_faces_crust_mantle > 0) then
        allocate(abs_boundary_ispec_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
        allocate(abs_boundary_ijk_crust_mantle(3,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
        allocate(abs_boundary_jacobian2Dw_crust_mantle(NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
        allocate(abs_boundary_normal_crust_mantle(NDIM,NGLLSQUARE,num_abs_boundary_faces_crust_mantle),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
        allocate(abs_boundary_npoin_crust_mantle(num_abs_boundary_faces_crust_mantle),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
        if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'
      else
        ! dummy arrays
        allocate(abs_boundary_ispec_crust_mantle(1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
        allocate(abs_boundary_ijk_crust_mantle(1,1,1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
        allocate(abs_boundary_jacobian2Dw_crust_mantle(1,1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
        allocate(abs_boundary_normal_crust_mantle(1,1,1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_normal')
        allocate(abs_boundary_npoin_crust_mantle(1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      endif
      abs_boundary_ispec_crust_mantle(:) = 0; abs_boundary_npoin_crust_mantle(:) = 0
      abs_boundary_ijk_crust_mantle(:,:,:) = 0
      abs_boundary_jacobian2Dw_crust_mantle(:,:) = 0.0; abs_boundary_normal_crust_mantle(:,:,:) = 0.0
    endif
    call bcast_all_i_for_database(abs_boundary_ispec_crust_mantle(1), size(abs_boundary_ispec_crust_mantle,kind=4))
    call bcast_all_i_for_database(abs_boundary_npoin_crust_mantle(1), size(abs_boundary_npoin_crust_mantle,kind=4))
    call bcast_all_i_for_database(abs_boundary_ijk_crust_mantle(1,1,1), size(abs_boundary_ijk_crust_mantle,kind=4))
    call bcast_all_cr_for_database(abs_boundary_jacobian2Dw_crust_mantle(1,1), size(abs_boundary_jacobian2Dw_crust_mantle,kind=4))
    call bcast_all_cr_for_database(abs_boundary_normal_crust_mantle(1,1,1), size(abs_boundary_normal_crust_mantle,kind=4))
  endif

  ! outer core
  if (NSPEC_OUTER_CORE > 0) then
    call bcast_all_i_for_database(num_abs_boundary_faces_outer_core,1)

    if (.not. I_should_read_the_database) then
      ! allocates absorbing boundary arrays
      if (num_abs_boundary_faces_outer_core > 0) then
        allocate(abs_boundary_ispec_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
        allocate(abs_boundary_ijk_outer_core(3,NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
        allocate(abs_boundary_jacobian2Dw_outer_core(NGLLSQUARE,num_abs_boundary_faces_outer_core),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
        allocate(abs_boundary_npoin_outer_core(num_abs_boundary_faces_outer_core),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
        if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'
      else
        ! dummy arrays
        allocate(abs_boundary_ispec_outer_core(1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ispec')
        allocate(abs_boundary_ijk_outer_core(1,1,1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_ijk')
        allocate(abs_boundary_jacobian2Dw_outer_core(1,1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_jacobian2Dw')
        allocate(abs_boundary_npoin_outer_core(1),stat=ier)
        if (ier /= 0) call exit_mpi(myrank,'Error allocating array abs_boundary_npoin')
      endif
      abs_boundary_ispec_outer_core(:) = 0; abs_boundary_npoin_outer_core(:) = 0
      abs_boundary_ijk_outer_core(:,:,:) = 0
      abs_boundary_jacobian2Dw_outer_core(:,:) = 0.0
    endif
    call bcast_all_i_for_database(abs_boundary_ispec_outer_core(1), size(abs_boundary_ispec_outer_core,kind=4))
    call bcast_all_i_for_database(abs_boundary_npoin_outer_core(1), size(abs_boundary_npoin_outer_core,kind=4))
    call bcast_all_i_for_database(abs_boundary_ijk_outer_core(1,1,1), size(abs_boundary_ijk_outer_core,kind=4))
    call bcast_all_cr_for_database(abs_boundary_jacobian2Dw_outer_core(1,1), size(abs_boundary_jacobian2Dw_outer_core,kind=4))
  endif

  end subroutine bcast_mesh_databases_stacey

