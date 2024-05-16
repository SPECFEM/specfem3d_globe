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

! compute the approximate amount of static memory needed to run the solver

  subroutine memory_eval(NEX_PER_PROC_XI,NEX_PER_PROC_ETA, &
                         NPROCTOT, &
                         NSPEC_REGIONS,NGLOB_REGIONS, &
                         NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
                         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
                         NSPEC_INNER_CORE_ATTENUATION, &
                         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
                         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
                         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
                         NSPEC_CRUST_MANTLE_ADJOINT, &
                         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
                         NSPEC_TRINFINITE_ADJOINT,NSPEC_INFINITE_ADJOINT, &
                         NGLOB_CRUST_MANTLE_ADJOINT, &
                         NGLOB_OUTER_CORE_ADJOINT,NGLOB_INNER_CORE_ADJOINT, &
                         NGLOB_TRINFINITE_ADJOINT,NGLOB_INFINITE_ADJOINT, &
                         NSPEC_OUTER_CORE_ROT_ADJOINT, &
                         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
                         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION, &
                         NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                         static_memory_size)

  use constants
  use shared_parameters, only: &
    ATT1,ATT2,ATT3, &
    APPROXIMATE_HESS_KL,ANISOTROPIC_KL,NOISE_TOMOGRAPHY, &
    EXACT_MASS_MATRIX_FOR_ROTATION, &
    OCEANS,ABSORBING_CONDITIONS,ATTENUATION,ANISOTROPIC_3D_MANTLE, &
    TRANSVERSE_ISOTROPY,ANISOTROPIC_INNER_CORE,ROTATION,TOPOGRAPHY,GRAVITY, &
    FULL_GRAVITY, &
    ONE_CRUST,NCHUNKS, &
    SIMULATION_TYPE,SAVE_FORWARD, &
    MOVIE_VOLUME,MOVIE_VOLUME_TYPE, &
    NX_BATHY,NY_BATHY

  use shared_parameters, only: ner_mesh_layers,ratio_sampling_array,doubling_index,this_region_has_a_doubling

  implicit none

  ! input
  integer, dimension(MAX_NUM_REGIONS), intent(in) :: NSPEC_REGIONS, NGLOB_REGIONS, NSPEC2D_BOTTOM, NSPEC2D_TOP
  integer, intent(in) :: NEX_PER_PROC_XI,NEX_PER_PROC_ETA
  integer, intent(in) :: NPROCTOT

  ! output
  double precision, intent(out) :: static_memory_size

  integer, intent(out) :: NSPECMAX_ANISO_IC,NSPECMAX_ISO_MANTLE,NSPECMAX_TISO_MANTLE, &
         NSPECMAX_ANISO_MANTLE,NSPEC_CRUST_MANTLE_ATTENUATION, &
         NSPEC_INNER_CORE_ATTENUATION, &
         NSPEC_CRUST_MANTLE_STR_OR_ATT,NSPEC_INNER_CORE_STR_OR_ATT, &
         NSPEC_CRUST_MANTLE_STR_AND_ATT,NSPEC_INNER_CORE_STR_AND_ATT, &
         NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_INNER_CORE_STRAIN_ONLY, &
         NSPEC_CRUST_MANTLE_ADJOINT, &
         NSPEC_OUTER_CORE_ADJOINT,NSPEC_INNER_CORE_ADJOINT, &
         NGLOB_CRUST_MANTLE_ADJOINT,NGLOB_OUTER_CORE_ADJOINT, &
         NGLOB_INNER_CORE_ADJOINT,NSPEC_OUTER_CORE_ROT_ADJOINT, &
         NSPEC_CRUST_MANTLE_STACEY,NSPEC_OUTER_CORE_STACEY, &
         NGLOB_CRUST_MANTLE_OCEANS,NSPEC_OUTER_CORE_ROTATION

  integer, intent(out) :: NSPEC_TRINFINITE_ADJOINT,NSPEC_INFINITE_ADJOINT, &
                          NGLOB_TRINFINITE_ADJOINT,NGLOB_INFINITE_ADJOINT

  ! local variables
  integer :: ilayer,NUMBER_OF_MESH_LAYERS,ner_without_doubling,ispec_aniso

  integer :: NSPEC_CRUST_MANTLE_ADJOINT_HESS,NSPEC_CRUST_MANTLE_ADJOINT_NOISE, &
             NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL, &
             NGLOB_XY_CM,NGLOB_XY_IC

  ! generate the elements in all the regions of the mesh
  ispec_aniso = 0

  if (ONE_CRUST) then
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS - 1
  else
    NUMBER_OF_MESH_LAYERS = MAX_NUMBER_OF_MESH_LAYERS
  endif

  ! count anisotropic elements
  do ilayer = 1, NUMBER_OF_MESH_LAYERS
    if (doubling_index(ilayer) == IFLAG_220_80 .or. doubling_index(ilayer) == IFLAG_80_MOHO) then
      ner_without_doubling = ner_mesh_layers(ilayer)
      if (this_region_has_a_doubling(ilayer)) then
        ner_without_doubling = ner_without_doubling - 2
        ispec_aniso = ispec_aniso + &
            (NSPEC_DOUBLING_SUPERBRICK*(NEX_PER_PROC_XI/ratio_sampling_array(ilayer)/2)* &
            (NEX_PER_PROC_ETA/ratio_sampling_array(ilayer)/2))
      endif
      ispec_aniso = ispec_aniso + &
        ((NEX_PER_PROC_XI/ratio_sampling_array(ilayer))*(NEX_PER_PROC_ETA/ratio_sampling_array(ilayer))*ner_without_doubling)
    endif
  enddo

  ! define static size of the arrays whose size depends on logical tests
  if (ANISOTROPIC_INNER_CORE) then
    NSPECMAX_ANISO_IC = NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPECMAX_ANISO_IC = 0
  endif

  NSPECMAX_ISO_MANTLE = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
  if (ANISOTROPIC_3D_MANTLE) then
    NSPECMAX_TISO_MANTLE = 0
    NSPECMAX_ANISO_MANTLE = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
  else
    if (TRANSVERSE_ISOTROPY) then
      ! note: the number of transverse isotropic elements is ispec_aniso
      !          however for transverse isotropic kernels, the arrays muhstore,kappahstore,eta_anisostore,
      !          will be needed for the crust_mantle region everywhere still...
      !          originally: NSPECMAX_TISO_MANTLE = ispec_aniso
      NSPECMAX_TISO_MANTLE = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    else
      NSPECMAX_TISO_MANTLE = 0
    endif
    NSPECMAX_ANISO_MANTLE = 0
  endif

  ! if attenuation is off, set dummy size of arrays to one
  if (ATTENUATION) then
    NSPEC_CRUST_MANTLE_ATTENUATION = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_ATTENUATION = NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_ATTENUATION = 0
    NSPEC_INNER_CORE_ATTENUATION = 0
  endif

  if (ATTENUATION .or. SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    NSPEC_CRUST_MANTLE_STR_OR_ATT = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STR_OR_ATT = NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STR_OR_ATT = 0
    NSPEC_INNER_CORE_STR_OR_ATT = 0
  endif

  if (ATTENUATION .and. &
    ( SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then
    NSPEC_CRUST_MANTLE_STR_AND_ATT = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STR_AND_ATT = NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STR_AND_ATT = 0
    NSPEC_INNER_CORE_STR_AND_ATT = 0
  endif

  if (SIMULATION_TYPE /= 1 .or. SAVE_FORWARD .or. (MOVIE_VOLUME .and. SIMULATION_TYPE /= 3)) then
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_INNER_CORE_STRAIN_ONLY = NSPEC_REGIONS(IREGION_INNER_CORE)
  else
    NSPEC_CRUST_MANTLE_STRAIN_ONLY = 0
    NSPEC_INNER_CORE_STRAIN_ONLY = 0
  endif

  ! adjoint sizes
  if ((SIMULATION_TYPE == 1 .and. SAVE_FORWARD) .or. SIMULATION_TYPE == 3) then
    NSPEC_CRUST_MANTLE_ADJOINT = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_OUTER_CORE_ADJOINT = NSPEC_REGIONS(IREGION_OUTER_CORE)
    NSPEC_INNER_CORE_ADJOINT = NSPEC_REGIONS(IREGION_INNER_CORE)
    NSPEC_TRINFINITE_ADJOINT = NSPEC_REGIONS(IREGION_TRINFINITE)
    NSPEC_INFINITE_ADJOINT = NSPEC_REGIONS(IREGION_INFINITE)

    NGLOB_CRUST_MANTLE_ADJOINT = NGLOB_REGIONS(IREGION_CRUST_MANTLE)
    NGLOB_OUTER_CORE_ADJOINT = NGLOB_REGIONS(IREGION_OUTER_CORE)
    NGLOB_INNER_CORE_ADJOINT = NGLOB_REGIONS(IREGION_INNER_CORE)
    NGLOB_TRINFINITE_ADJOINT = NGLOB_REGIONS(IREGION_TRINFINITE)
    NGLOB_INFINITE_ADJOINT = NGLOB_REGIONS(IREGION_INFINITE)

    if (ROTATION) then
      NSPEC_OUTER_CORE_ROT_ADJOINT = NSPEC_REGIONS(IREGION_OUTER_CORE)
    else
      NSPEC_OUTER_CORE_ROT_ADJOINT = 0
    endif
  else
    NSPEC_CRUST_MANTLE_ADJOINT = 0
    NSPEC_OUTER_CORE_ADJOINT = 0
    NSPEC_INNER_CORE_ADJOINT = 0
    NSPEC_TRINFINITE_ADJOINT = 0
    NSPEC_INFINITE_ADJOINT = 0

    NGLOB_CRUST_MANTLE_ADJOINT = 0
    NGLOB_OUTER_CORE_ADJOINT = 0
    NGLOB_INNER_CORE_ADJOINT = 0
    NGLOB_TRINFINITE_ADJOINT = 0
    NGLOB_INFINITE_ADJOINT = 0

    NSPEC_OUTER_CORE_ROT_ADJOINT = 0
  endif

  ! if absorbing conditions are off, set dummy size of arrays to one
  if (ABSORBING_CONDITIONS) then
    NSPEC_CRUST_MANTLE_STACEY = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NSPEC_OUTER_CORE_STACEY = NSPEC_REGIONS(IREGION_OUTER_CORE)
  else
    NSPEC_CRUST_MANTLE_STACEY = 0
    NSPEC_OUTER_CORE_STACEY = 0
  endif

  ! if oceans are off, set dummy size of arrays to one
  if (OCEANS) then
    NGLOB_CRUST_MANTLE_OCEANS = NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  else
    NGLOB_CRUST_MANTLE_OCEANS = 0
  endif

  if (ROTATION) then
    NSPEC_OUTER_CORE_ROTATION = NSPEC_REGIONS(IREGION_OUTER_CORE)
  else
    NSPEC_OUTER_CORE_ROTATION = 0
  endif

! add size of each set of static arrays multiplied by the number of such arrays

  static_memory_size = 0.d0

! crust/mantle

  ! ibool_crust_mantle
  static_memory_size = static_memory_size + &
    dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_CRUST_MANTLE)*dble(SIZE_INTEGER)

  ! xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle
  ! etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,
  ! gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  static_memory_size = static_memory_size + &
    9.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

  ! xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle,rmass_crust_mantle
  if (NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then
     ! three mass matrices for the crust and mantle region: rmassx, rmassy and rmassz
     static_memory_size = static_memory_size + &
          6.d0*NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)
  else
     ! one only keeps one mass matrix for the calculations: rmassz
     static_memory_size = static_memory_size + &
          4.d0*NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)
  endif

  ! rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle
  static_memory_size = static_memory_size + &
    3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ISO_MANTLE*dble(CUSTOM_REAL)

  ! kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle
  static_memory_size = static_memory_size + &
    3.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_TISO_MANTLE*dble(CUSTOM_REAL)

  ! c11,.. store for tiso elements
  static_memory_size = static_memory_size + &
    21.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_TISO_MANTLE*dble(CUSTOM_REAL)

  ! for aniso elements
  ! c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,
  ! c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle,
  ! c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle,
  ! c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle,
  ! c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle,
  ! c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle,
  ! c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle
  static_memory_size = static_memory_size + &
    21.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ANISO_MANTLE*dble(CUSTOM_REAL)

  ! mu0
  static_memory_size = static_memory_size + &
    1.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ANISO_MANTLE*dble(CUSTOM_REAL)

  ! ispec_is_tiso_crust_mantle
  static_memory_size = static_memory_size + NSPEC_REGIONS(IREGION_CRUST_MANTLE)*dble(SIZE_LOGICAL)

  ! displ_crust_mantle,veloc_crust_mantle,accel_crust_mantle
  static_memory_size = static_memory_size + &
    3.d0*dble(NDIM)*NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

  ! attenuation arrays
  ! one_minus_sum_beta_crust_mantle, factor_scale_crust_mantle
  static_memory_size = static_memory_size + &
      2.d0*dble(ATT1)*dble(ATT2)*dble(ATT3)*NSPEC_CRUST_MANTLE_ATTENUATION*dble(CUSTOM_REAL)

  ! factor_common_crust_mantle
  static_memory_size = static_memory_size + &
      dble(N_SLS)*dble(ATT1)*dble(ATT2)*dble(ATT3)*NSPEC_CRUST_MANTLE_ATTENUATION*dble(CUSTOM_REAL)

  ! R_memory_crust_mantle (R_xx, R_yy, ..)
  static_memory_size = static_memory_size + &
    5.d0*dble(N_SLS)*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ATTENUATION*dble(CUSTOM_REAL)

  ! add arrays used to save strain for attenuation or for adjoint runs
  ! epsilondev_crust_mantle
  static_memory_size = static_memory_size + &
    5.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_STR_OR_ATT*dble(CUSTOM_REAL)

  ! eps_trace_over_3_crust_mantle
  static_memory_size = static_memory_size + &
    dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_STRAIN_ONLY*dble(CUSTOM_REAL)

  ! normal_bottom_crust_mantle
  static_memory_size = static_memory_size + &
    dble(NDIM)*dble(NGLLX)*dble(NGLLY)*NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

  ! normal_top_crust_mantle
  static_memory_size = static_memory_size + &
    dble(NDIM)*dble(NGLLX)*dble(NGLLY)*NSPEC2D_TOP(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

  ! rho_vp_crust_mantle,rho_vs_crust_mantle
  static_memory_size = static_memory_size + &
    2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_STACEY*dble(CUSTOM_REAL)

! inner core

  ! ibool_inner_core
  static_memory_size = static_memory_size + &
    dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_INNER_CORE)*dble(SIZE_INTEGER)

  ! xix_inner_core,xiy_inner_core,xiz_inner_core,
  ! etax_inner_core,etay_inner_core,etaz_inner_core,
  ! gammax_inner_core,gammay_inner_core,gammaz_inner_core,
  ! rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core
  static_memory_size = static_memory_size + &
    12.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

  ! xstore_inner_core,ystore_inner_core,zstore_inner_core,rmassz_inner_core
  static_memory_size = static_memory_size + &
    4.d0*NGLOB_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

  ! c11store_inner_core,c33store_inner_core,c12store_inner_core,c13store_inner_core,c44store_inner_core
  static_memory_size = static_memory_size + &
    5.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPECMAX_ANISO_IC*dble(CUSTOM_REAL)

  ! idoubling_inner_core
  static_memory_size = static_memory_size + NSPEC_REGIONS(IREGION_INNER_CORE)*dble(SIZE_INTEGER)

  ! displ_inner_core,veloc_inner_core,accel_inner_core
  static_memory_size = static_memory_size + &
    3.d0*dble(NDIM)*NGLOB_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)

  ! attenuation arrays
  ! one_minus_sum_beta_inner_core, factor_scale_inner_core
  static_memory_size = static_memory_size + &
    2.d0*dble(ATT1)*dble(ATT2)*dble(ATT3)*NSPEC_INNER_CORE_ATTENUATION*dble(CUSTOM_REAL)

  ! factor_common_inner_core
  static_memory_size = static_memory_size + &
    dble(N_SLS)*dble(ATT1)*dble(ATT2)*dble(ATT3)*NSPEC_INNER_CORE_ATTENUATION*dble(CUSTOM_REAL)

  ! R_memory_inner_core
  static_memory_size = static_memory_size + &
    5.d0*dble(N_SLS)*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_ATTENUATION*dble(CUSTOM_REAL)

  ! add arrays used to save strain for attenuation or for adjoint runs
  ! epsilondev_inner_core
  static_memory_size = static_memory_size + &
    5.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_STR_OR_ATT*dble(CUSTOM_REAL)

  ! eps_trace_over_3_inner_core
  static_memory_size = static_memory_size + &
    dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_STRAIN_ONLY*dble(CUSTOM_REAL)

! outer core

  ! ibool_outer_core
  static_memory_size = static_memory_size + &
    dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_OUTER_CORE)*dble(SIZE_INTEGER)

  ! xix_outer_core,xiy_outer_core,xiz_outer_core,
  ! etax_outer_core,etay_outer_core,etaz_outer_core,
  ! gammax_outer_core,gammay_outer_core,gammaz_outer_core
  ! rhostore_outer_core,kappavstore_outer_core
  static_memory_size = static_memory_size + &
    11.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

  ! xstore_outer_core, ystore_outer_core, zstore_outer_core, rmass_outer_core,
  ! displ_outer_core, veloc_outer_core, accel_outer_core
  static_memory_size = static_memory_size + &
    7.d0*NGLOB_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

  ! vp_outer_core
  static_memory_size = static_memory_size + &
    dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_STACEY*dble(CUSTOM_REAL)

  ! ispec_is_tiso_outer_core
  static_memory_size = static_memory_size + NSPEC_REGIONS(IREGION_OUTER_CORE)*dble(SIZE_LOGICAL)

! additional arrays

  ! A_array_rotation,B_array_rotation
  static_memory_size = static_memory_size + &
    2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ROTATION*dble(CUSTOM_REAL)

  ! GRAVITY
  if (GRAVITY) then
    ! minus_gravity_table, &
    ! minus_deriv_gravity_table,density_table,d_ln_density_dr_table,minus_rho_g_over_kappa_fluid
    static_memory_size = static_memory_size + &
      5.d0*NRAD_GRAVITY*dble(SIZE_DOUBLE)

    ! gravity_pre_store_crust_mantle,gravity_H_crust_mantle
    static_memory_size = static_memory_size + &
      (3.d0 + 6.d0)*NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

    ! gravity_pre_store_inner_core,gravity_H_inner_core
    static_memory_size = static_memory_size + &
      (3.d0 + 6.d0)*NGLOB_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)
  endif

  ! gravity_pre_store_outer_core
  static_memory_size = static_memory_size + &
    3.d0*NGLOB_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)

  ! ELLIPTICITY
  ! rspl,ellipicity_spline,ellipicity_spline2
  static_memory_size = static_memory_size + &
    3.d0*NR_DENSITY*dble(SIZE_DOUBLE)

  ! OCEANS
  ! rmass_ocean_load
  static_memory_size = static_memory_size + &
    NGLOB_CRUST_MANTLE_OCEANS*dble(CUSTOM_REAL)

  ! not accounted for yet (npoin_oceans unknown yet): rmass_ocean_load_selected,normal_ocean_load,ibool_ocean_load

  ! ichunk_slice,iproc_xi_slice,iproc_eta_slice,addressing
  static_memory_size = static_memory_size + &
      4.d0*NPROCTOT*dble(SIZE_INTEGER)

  if (TOPOGRAPHY .or. OCEANS) then
    ! ibathy_topo
    static_memory_size = static_memory_size + &
      NX_BATHY*NY_BATHY*dble(SIZE_INTEGER)
  endif

  if (MOVIE_VOLUME) then
    ! iepsilondev_.. crust_mantle
    static_memory_size = static_memory_size + &
      5.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

    ! muvstore_crust_mantle_3dmovie
    static_memory_size = static_memory_size + &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)

    ! mask_ibool_3dmovie
    static_memory_size = static_memory_size + &
      NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(SIZE_LOGICAL)

  endif

  ! div_displ_outer_core
  if (MOVIE_VOLUME .and. MOVIE_VOLUME_TYPE == 4) then
    static_memory_size = static_memory_size + &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)
  endif

! add arrays used for adjoint runs only (LQY: not very accurate)

  ! b_R_memory_crust_mantle
  ! b_epsilondev_crust_mantle
  ! b_eps_trace_over_3_crust_mantle
  ! rho_kl_crust_mantle,beta_kl_crust_mantle, alpha_kl_crust_mantle
  static_memory_size = static_memory_size + (5.d0*dble(N_SLS) + 9.d0)* &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ADJOINT*dble(CUSTOM_REAL)

  ! rho_kl_outer_core,alpha_kl_outer_core
  static_memory_size = static_memory_size + &
    2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ADJOINT*dble(CUSTOM_REAL)

  ! b_R_memory_inner_core
  ! b_epsilondev_inner_core
  ! b_eps_trace_over_3_inner_core
  ! rho_kl_inner_core,beta_kl_inner_core, alpha_kl_inner_core
  static_memory_size = static_memory_size + (5.d0*dble(N_SLS) + 9.d0)* &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_INNER_CORE_ADJOINT*dble(CUSTOM_REAL)

  ! b_displ_crust_mantle,b_veloc_crust_mantle,b_accel_crust_mantle
  static_memory_size = static_memory_size + &
    3.d0*dble(NDIM)*NGLOB_CRUST_MANTLE_ADJOINT*dble(CUSTOM_REAL)

  ! b_displ_outer_core,b_veloc_outer_core,b_accel_outer_core
  static_memory_size = static_memory_size + &
    3.d0*NGLOB_OUTER_CORE_ADJOINT*dble(CUSTOM_REAL)

  ! vector_accel_outer_core,vector_displ_outer_core,b_vector_displ_outer_core
  static_memory_size = static_memory_size + &
    3.d0*dble(NDIM)*NGLOB_OUTER_CORE_ADJOINT*dble(CUSTOM_REAL)

  ! b_displ_inner_core,b_veloc_inner_core,b_accel_inner_core
  static_memory_size = static_memory_size + &
    3.d0*dble(NDIM)*NGLOB_INNER_CORE_ADJOINT*dble(CUSTOM_REAL)

  ! b_A_array_rotation,b_B_array_rotation
  static_memory_size = static_memory_size + &
    2.d0*dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_OUTER_CORE_ROT_ADJOINT*dble(CUSTOM_REAL)

  ! cijkl_kl_crust_mantle (full anisotropic kernels with 21 coefficients)
  if (ANISOTROPIC_KL) then
    NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL = NSPEC_CRUST_MANTLE_ADJOINT
  else
    NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL = 0
  endif
  static_memory_size = static_memory_size + 21.d0* &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ADJOINT_ANISO_KL*dble(CUSTOM_REAL)

  ! hess_kl_crust_mantle
  if (APPROXIMATE_HESS_KL) then
    NSPEC_CRUST_MANTLE_ADJOINT_HESS = NSPEC_CRUST_MANTLE_ADJOINT
  else
    NSPEC_CRUST_MANTLE_ADJOINT_HESS = 0
  endif
  static_memory_size = static_memory_size + &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ADJOINT_HESS*dble(CUSTOM_REAL)

  ! sigma_kl_crust_mantle
  if (NOISE_TOMOGRAPHY > 0) then
    NSPEC_CRUST_MANTLE_ADJOINT_NOISE = NSPEC_CRUST_MANTLE_ADJOINT
  else
    NSPEC_CRUST_MANTLE_ADJOINT_NOISE = 0
  endif
  static_memory_size = static_memory_size + &
      dble(NGLLX)*dble(NGLLY)*dble(NGLLZ)*NSPEC_CRUST_MANTLE_ADJOINT_NOISE*dble(CUSTOM_REAL)

  ! noise_surface_movie_buffer
  if (NOISE_TOMOGRAPHY > 0) then
    ! noise_surface_movie_buffer
    static_memory_size = static_memory_size + &
        dble(NDIM)*dble(NGLLX)*dble(NGLLY)*dble(NSPEC2D_TOP(IREGION_CRUST_MANTLE))*dble(CUSTOM_REAL)
    ! noise buffer size not known yet
  endif

  ! in the case of Stacey boundary conditions, add C*delta/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be fictitious / unused
  if (NCHUNKS /= 6 .and. ABSORBING_CONDITIONS) then
     NGLOB_XY_CM = NGLOB_REGIONS(IREGION_CRUST_MANTLE)
  else
     NGLOB_XY_CM = 0
  endif
  NGLOB_XY_IC = 0

  if (ROTATION .and. EXACT_MASS_MATRIX_FOR_ROTATION) then
    NGLOB_XY_CM = NGLOB_REGIONS(IREGION_CRUST_MANTLE)
    NGLOB_XY_IC = NGLOB_REGIONS(IREGION_INNER_CORE)
  endif

  ! rmassx_crust_mantle,rmassy_crust_mantle for EXACT_MASS_MATRIX_FOR_ROTATION and/or ABSORBING_CONDITIONS
  static_memory_size = static_memory_size + 2.d0*NGLOB_XY_CM*4.d0*dble(CUSTOM_REAL)

  if (SIMULATION_TYPE == 3) then
    ! b_rmassx_crust_mantle,b_rmassy_crust_mantle for EXACT_MASS_MATRIX_FOR_ROTATION and/or ABSORBING_CONDITIONS
    static_memory_size = static_memory_size + 2.d0*NGLOB_XY_CM*4.d0*dble(CUSTOM_REAL)
  endif

  ! rmassx_inner_core,rmassy_inner_core for EXACT_MASS_MATRIX_FOR_ROTATION and/or ABSORBING_CONDITIONS
  static_memory_size = static_memory_size + 2.d0*NGLOB_XY_IC*4.d0*dble(CUSTOM_REAL)

  if (SIMULATION_TYPE == 3) then
    ! b_rmassx_inner_core,b_rmassy_inner_core for EXACT_MASS_MATRIX_FOR_ROTATION and/or ABSORBING_CONDITIONS
    static_memory_size = static_memory_size + 2.d0*NGLOB_XY_IC*4.d0*dble(CUSTOM_REAL)
  endif

  ! full gravity
  if (FULL_GRAVITY) then
    !TODO: add all full gravity array sizes?
    ! pgrav_ic,pgrav_oc,pgrav_cm,pgrav_trinf,pgrav_inf
    static_memory_size = static_memory_size + NGLOB_REGIONS(IREGION_INNER_CORE)*dble(CUSTOM_REAL)
    static_memory_size = static_memory_size + NGLOB_REGIONS(IREGION_OUTER_CORE)*dble(CUSTOM_REAL)
    static_memory_size = static_memory_size + NGLOB_REGIONS(IREGION_CRUST_MANTLE)*dble(CUSTOM_REAL)
    static_memory_size = static_memory_size + NGLOB_REGIONS(IREGION_TRINFINITE)*dble(CUSTOM_REAL)
    static_memory_size = static_memory_size + NGLOB_REGIONS(IREGION_INFINITE)*dble(CUSTOM_REAL)

    if (SIMULATION_TYPE == 3) then
      ! b_pgrav_ic,b_pgrav_oc,b_pgrav_cm,b_pgrav_trinf,b_pgrav_inf
      static_memory_size = static_memory_size + NGLOB_INNER_CORE_ADJOINT*dble(CUSTOM_REAL)
      static_memory_size = static_memory_size + NGLOB_OUTER_CORE_ADJOINT*dble(CUSTOM_REAL)
      static_memory_size = static_memory_size + NGLOB_CRUST_MANTLE_ADJOINT*dble(CUSTOM_REAL)
      static_memory_size = static_memory_size + NGLOB_TRINFINITE_ADJOINT*dble(CUSTOM_REAL)
      static_memory_size = static_memory_size + NGLOB_INFINITE_ADJOINT*dble(CUSTOM_REAL)
    endif
  endif

  end subroutine memory_eval

