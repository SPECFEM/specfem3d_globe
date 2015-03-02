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

  subroutine compute_boundary_kernels()

! kernel calculations

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  implicit none

  ! local parameters
  ! dummy array that does not need to be actually read
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  logical,dimension(:),allocatable:: dummy_ispec_is_tiso
  integer:: ispec,i,j,k,l,iglob
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: gradx,grady,gradz

  logical,dimension(:),allocatable :: mask_ibool
  integer :: ier

  ! transfers wavefields onto CPU
  if (GPU_MODE) then

    ! crust/mantle
    call transfer_accel_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,accel_crust_mantle,Mesh_pointer)
    call transfer_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,displ_crust_mantle,Mesh_pointer)
    call transfer_b_displ_cm_from_device(NDIM*NGLOB_CRUST_MANTLE,b_displ_crust_mantle,Mesh_pointer)
    ! inner core
    call transfer_accel_ic_from_device(NDIM*NGLOB_INNER_CORE,accel_inner_core,Mesh_pointer)
    call transfer_displ_ic_from_device(NDIM*NGLOB_INNER_CORE,displ_inner_core,Mesh_pointer)
    call transfer_b_displ_ic_from_device(NDIM*NGLOB_INNER_CORE,b_displ_inner_core,Mesh_pointer)
    ! outer core
    call transfer_accel_oc_from_device(NGLOB_OUTER_CORE,accel_outer_core,Mesh_pointer)
    call transfer_displ_oc_from_device(NGLOB_OUTER_CORE,displ_outer_core,Mesh_pointer)
    call transfer_b_displ_oc_from_device(NGLOB_OUTER_CORE,b_displ_outer_core,Mesh_pointer)

    ! pre-calculates gradients in outer core
    ! note: for CPU, this is already done in compute_kernels_outer_core() routine

    allocate(mask_ibool(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating mask_ibool array in routine compute_boundary_kernels()')
    mask_ibool(:) = .false.

    ! pre-calculates gradients
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

            ! global index
            iglob = ibool_outer_core(i,j,k,ispec)

            ! only calculates gradients once for shared nodes
            if (.not. mask_ibool(iglob)) then

              ! masks this global point
              mask_ibool(iglob) = .true.

              xixl = xix_outer_core(i,j,k,ispec)
              xiyl = xiy_outer_core(i,j,k,ispec)
              xizl = xiz_outer_core(i,j,k,ispec)
              etaxl = etax_outer_core(i,j,k,ispec)
              etayl = etay_outer_core(i,j,k,ispec)
              etazl = etaz_outer_core(i,j,k,ispec)
              gammaxl = gammax_outer_core(i,j,k,ispec)
              gammayl = gammay_outer_core(i,j,k,ispec)
              gammazl = gammaz_outer_core(i,j,k,ispec)

              ! calculates gradient grad(b_displ)
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempx3l = 0._CUSTOM_REAL
              do l = 1,NGLLX
                tempx1l = tempx1l + b_displ_outer_core(ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              enddo
              do l = 1,NGLLY
                tempx2l = tempx2l + b_displ_outer_core(ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              enddo
              do l = 1,NGLLZ
                tempx3l = tempx3l + b_displ_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              enddo
              gradx = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
              grady = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
              gradz = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

              ! assigns and stores gradient on global nodes
              b_vector_displ_outer_core(1,iglob) = gradx
              b_vector_displ_outer_core(2,iglob) = grady
              b_vector_displ_outer_core(3,iglob) = gradz

              ! calculates gradient grad(accel)
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempx3l = 0._CUSTOM_REAL
              do l = 1,NGLLX
                tempx1l = tempx1l + accel_outer_core(ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              enddo
              do l = 1,NGLLY
                tempx2l = tempx2l + accel_outer_core(ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              enddo
              do l = 1,NGLLZ
                tempx3l = tempx3l + accel_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              enddo
              gradx = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
              grady = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
              gradz = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

              vector_accel_outer_core(1,iglob) = gradx
              vector_accel_outer_core(2,iglob) = grady
              vector_accel_outer_core(3,iglob) = gradz

              ! calculates gradient grad(displ)
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempx3l = 0._CUSTOM_REAL
              do l = 1,NGLLX
                tempx1l = tempx1l + displ_outer_core(ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
              enddo
              do l = 1,NGLLY
                tempx2l = tempx2l + displ_outer_core(ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
              enddo
              do l = 1,NGLLZ
                tempx3l = tempx3l + displ_outer_core(ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              enddo
              gradx = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
              grady = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
              gradz = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

              vector_displ_outer_core(1,iglob) = gradx
              vector_displ_outer_core(2,iglob) = grady
              vector_displ_outer_core(3,iglob) = gradz
            endif

          enddo
        enddo
      enddo
    enddo

    ! frees memory
    deallocate(mask_ibool)

  endif ! GPU_MODE

  ! updates kernels on CPU
  fluid_solid_boundary = .false.
  iregion_code = IREGION_CRUST_MANTLE

  ! Moho
  if (.not. SUPPRESS_CRUSTAL_MESH .and. HONOR_1D_SPHERICAL_MOHO) then
    call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                       nspec_crust_mantle,iregion_code, &
                                       ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                       xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                       etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                       gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                       hprime_xx,hprime_yy,hprime_zz, &
                                       rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                       kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                       c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                       c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                       c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                       c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                       c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                       c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                       k_top,ibelm_moho_top,normal_moho,moho_kl_top,fluid_solid_boundary,NSPEC2D_MOHO)

    call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                       nspec_crust_mantle,iregion_code, &
                                       ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                       xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                       etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                       gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                       hprime_xx,hprime_yy,hprime_zz, &
                                       rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                       kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                       c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                       c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                       c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                       c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                       c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                       c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                       k_bot,ibelm_moho_bot,normal_moho,moho_kl_bot,fluid_solid_boundary,NSPEC2D_MOHO)

    moho_kl = moho_kl + (moho_kl_top - moho_kl_bot) * deltat
  endif

  ! 400
  call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                     nspec_crust_mantle,iregion_code, &
                                     ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                     c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                     c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                     c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                     c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     k_top,ibelm_400_top,normal_400,d400_kl_top,fluid_solid_boundary,NSPEC2D_400)

  call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                     nspec_crust_mantle,iregion_code, &
                                     ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                     c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                     c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                     c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                     c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     k_bot,ibelm_400_bot,normal_400,d400_kl_bot,fluid_solid_boundary,NSPEC2D_400)

  d400_kl = d400_kl + (d400_kl_top - d400_kl_bot) * deltat

  ! 670
  call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                     nspec_crust_mantle,iregion_code, &
                                     ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                     c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                     c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                     c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                     c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     k_top,ibelm_670_top,normal_670,d670_kl_top,fluid_solid_boundary,NSPEC2D_670)

  call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                     nspec_crust_mantle,iregion_code, &
                                     ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                     c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                     c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                     c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                     c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     k_bot,ibelm_670_bot,normal_670,d670_kl_bot,fluid_solid_boundary,NSPEC2D_670)

  d670_kl = d670_kl + (d670_kl_top - d670_kl_bot) * deltat

  ! CMB
  fluid_solid_boundary = .true.
  iregion_code = IREGION_CRUST_MANTLE

  call compute_boundary_kernel_depth(displ_crust_mantle,accel_crust_mantle,b_displ_crust_mantle, &
                                     nspec_crust_mantle,iregion_code, &
                                     ystore_crust_mantle,zstore_crust_mantle,ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
                                     xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                     etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle,&
                                     gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                     hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_crust_mantle,kappavstore_crust_mantle, muvstore_crust_mantle, &
                                     kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
                                     c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle,c14store_crust_mantle, &
                                     c15store_crust_mantle,c16store_crust_mantle,c22store_crust_mantle, &
                                     c23store_crust_mantle,c24store_crust_mantle,c25store_crust_mantle,c26store_crust_mantle, &
                                     c33store_crust_mantle,c34store_crust_mantle,c35store_crust_mantle, &
                                     c36store_crust_mantle,c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
                                     c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
                                     k_top,ibelm_bottom_crust_mantle,normal_top_outer_core, &
                                     cmb_kl_top,fluid_solid_boundary,NSPEC2D_CMB)

  iregion_code = IREGION_OUTER_CORE

  ! dummy allocation
  allocate(dummy_ispec_is_tiso(NSPEC_OUTER_CORE))
  dummy_ispec_is_tiso(:) = .false.

  call compute_boundary_kernel_depth(vector_displ_outer_core,vector_accel_outer_core,b_vector_displ_outer_core, &
                                     nspec_outer_core, &
                                     iregion_code,ystore_outer_core,zstore_outer_core,ibool_outer_core,dummy_ispec_is_tiso, &
                                     xix_outer_core,xiy_outer_core,xiz_outer_core, &
                                     etax_outer_core,etay_outer_core,etaz_outer_core,&
                                     gammax_outer_core,gammay_outer_core,gammaz_outer_core,hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_outer_core,kappavstore_outer_core,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     k_bot,ibelm_top_outer_core,normal_top_outer_core, &
                                     cmb_kl_bot,fluid_solid_boundary,NSPEC2D_CMB)

  cmb_kl = cmb_kl + (cmb_kl_top - cmb_kl_bot) * deltat

  ! ICB
  fluid_solid_boundary = .true.
  call compute_boundary_kernel_depth(vector_displ_outer_core,vector_accel_outer_core,b_vector_displ_outer_core, &
                                     nspec_outer_core, &
                                     iregion_code,ystore_outer_core,zstore_outer_core,ibool_outer_core,dummy_ispec_is_tiso, &
                                     xix_outer_core,xiy_outer_core,xiz_outer_core, &
                                     etax_outer_core,etay_outer_core,etaz_outer_core,&
                                     gammax_outer_core,gammay_outer_core,gammaz_outer_core,hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_outer_core,kappavstore_outer_core,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     k_top,ibelm_bottom_outer_core,normal_bottom_outer_core, &
                                     icb_kl_top,fluid_solid_boundary,NSPEC2D_ICB)

  deallocate(dummy_ispec_is_tiso)

  ! ICB
  iregion_code = IREGION_INNER_CORE

  ! dummy allocation
  allocate(dummy_ispec_is_tiso(NSPEC_INNER_CORE))
  dummy_ispec_is_tiso(:) = .false.

  call compute_boundary_kernel_depth(displ_inner_core,accel_inner_core,b_displ_inner_core, &
                                     nspec_inner_core,iregion_code, &
                                     ystore_inner_core,zstore_inner_core,ibool_inner_core,dummy_ispec_is_tiso, &
                                     xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                     etax_inner_core,etay_inner_core,etaz_inner_core,&
                                     gammax_inner_core,gammay_inner_core,gammaz_inner_core,hprime_xx,hprime_yy,hprime_zz, &
                                     rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
                                     dummy_array,dummy_array,dummy_array, &
                                     c11store_inner_core,c12store_inner_core,c13store_inner_core,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array,dummy_array, &
                                     c33store_inner_core,dummy_array,dummy_array, &
                                     dummy_array,c44store_inner_core,dummy_array,dummy_array, &
                                     dummy_array,dummy_array,dummy_array, &
                                     k_bot,ibelm_top_inner_core,normal_bottom_outer_core, &
                                     icb_kl_bot,fluid_solid_boundary,NSPEC2D_ICB)
  deallocate(dummy_ispec_is_tiso)

  icb_kl = icb_kl + (icb_kl_top - icb_kl_bot) * deltat

  end subroutine compute_boundary_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_boundary_kernel_depth(displ,accel,b_displ, &
                                           nspec,iregion_code, &
                                           ystore,zstore,ibool,ispec_is_tiso, &
                                           xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                                           hprime_xx,hprime_yy,hprime_zz, &
                                           rhostore,kappavstore,muvstore,kappahstore,muhstore,eta_anisostore, &
                                           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                                           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                                           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                                           k_disc,ibelm_disc,normal_disc,b_kl,fluid_solid_boundary,NSPEC2D_DISC)

  use constants

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,*) :: displ,accel,b_displ
  integer nspec, iregion_code
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool
  logical, dimension(*) :: ispec_is_tiso
  real(kind=CUSTOM_REAL), dimension(*) :: ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,*) :: xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,*) :: rhostore, kappavstore,muvstore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,*) :: kappahstore,muhstore,eta_anisostore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,*) :: c11store,c12store,c13store,c14store,c15store,c16store, &
             c22store,c23store,c24store,c25store,c26store,c33store, c34store,c35store,c36store, &
             c44store,c45store,c46store,c55store,c56store,c66store

  integer NSPEC2D_DISC, k_disc
  integer :: ibelm_disc(NSPEC2D_DISC)
  real(kind=CUSTOM_REAL) :: normal_disc(NDIM,NGLLX,NGLLY,NSPEC2D_DISC)
  real(kind=CUSTOM_REAL) :: b_kl(NGLLX,NGLLY,NSPEC2D_DISC)
  logical :: fluid_solid_boundary

  ! --- local variables ---
  integer ispec2D,i,j,k,iglob,ispec
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: displl, accell, b_displl, Kdvect
  real(kind=CUSTOM_REAL), dimension(NDIM) :: normal, temp1, temp2, temp3
  real(kind=CUSTOM_REAL) :: xixl, xiyl, xizl, etaxl, etayl, etazl, gammaxl, gammayl, gammazl
  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ) :: dsdx, sigma, b_dsdx, b_sigma
  real(kind=CUSTOM_REAL) :: b_kl_2(NGLLX,NGLLY)
  real(kind=CUSTOM_REAL) :: dKdx(NDIM,NDIM)

  ! ------

  ! initialization
  b_kl = 0.

  do ispec2D = 1, NSPEC2D_DISC

    ! assign local matrices
    ispec = ibelm_disc(ispec2D)
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX
          iglob = ibool(i,j,k,ispec)
          displl(:,i,j,k) = displ(:,iglob)
          accell(:,i,j,k) = accel(:,iglob)
          b_displl(:,i,j,k) = b_displ(:,iglob)
        enddo
      enddo
    enddo

    ! strain and stress
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          normal(:) = normal_disc(:,i,j,ispec2D)
          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! ----- adjoint strain ------
          temp1(:) = matmul(displl(:,:,j,k), hprime_xx(i,:))
          temp2(:) = matmul(displl(:,i,:,k), hprime_yy(j,:))
          temp3(:) = matmul(displl(:,i,j,:), hprime_zz(k,:))

          dsdx(1,1,i,j,k) = xixl*temp1(1) + etaxl*temp2(1) + gammaxl*temp3(1)
          dsdx(1,2,i,j,k) = xiyl*temp1(1) + etayl*temp2(1) + gammayl*temp3(1)
          dsdx(1,3,i,j,k) = xizl*temp1(1) + etazl*temp2(1) + gammazl*temp3(1)

          dsdx(2,1,i,j,k) = xixl*temp1(2) + etaxl*temp2(2) + gammaxl*temp3(2)
          dsdx(2,2,i,j,k) = xiyl*temp1(2) + etayl*temp2(2) + gammayl*temp3(2)
          dsdx(2,3,i,j,k) = xizl*temp1(2) + etazl*temp2(2) + gammazl*temp3(2)

          dsdx(3,1,i,j,k) = xixl*temp1(3) + etaxl*temp2(3) + gammaxl*temp3(3)
          dsdx(3,2,i,j,k) = xiyl*temp1(3) + etayl*temp2(3) + gammayl*temp3(3)
          dsdx(3,3,i,j,k) = xizl*temp1(3) + etazl*temp2(3) + gammazl*temp3(3)

          ! ------ adjoint stress -------
          call compute_stress_from_strain(dsdx(:,:,i,j,k),sigma(:,:,i,j,k),i,j,k,ispec,iregion_code, &
                     kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                     c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                     c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                     c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                     ystore,zstore,ibool,ispec_is_tiso)

          ! ----- forward strain -------
          temp1(:) = matmul(b_displl(:,:,j,k), hprime_xx(i,:))
          temp2(:) = matmul(b_displl(:,i,:,k), hprime_yy(j,:))
          temp3(:) = matmul(b_displl(:,i,j,:), hprime_zz(k,:))

          b_dsdx(1,1,i,j,k) = xixl*temp1(1) + etaxl*temp2(1) + gammaxl*temp3(1)
          b_dsdx(1,2,i,j,k) = xiyl*temp1(1) + etayl*temp2(1) + gammayl*temp3(1)
          b_dsdx(1,3,i,j,k) = xizl*temp1(1) + etazl*temp2(1) + gammazl*temp3(1)

          b_dsdx(2,1,i,j,k) = xixl*temp1(2) + etaxl*temp2(2) + gammaxl*temp3(2)
          b_dsdx(2,2,i,j,k) = xiyl*temp1(2) + etayl*temp2(2) + gammayl*temp3(2)
          b_dsdx(2,3,i,j,k) = xizl*temp1(2) + etazl*temp2(2) + gammazl*temp3(2)

          b_dsdx(3,1,i,j,k) = xixl*temp1(3) + etaxl*temp2(3) + gammaxl*temp3(3)
          b_dsdx(3,2,i,j,k) = xiyl*temp1(3) + etayl*temp2(3) + gammayl*temp3(3)
          b_dsdx(3,3,i,j,k) = xizl*temp1(3) + etazl*temp2(3) + gammazl*temp3(3)

          ! ----- forward stress ---------
          call compute_stress_from_strain(b_dsdx(:,:,i,j,k),b_sigma(:,:,i,j,k),i,j,k,ispec,iregion_code, &
                     kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                     c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                     c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                     c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                     ystore,zstore,ibool,ispec_is_tiso)

          ! ---- precompute K_d for F-S boundaries ----
          if (fluid_solid_boundary) then
            Kdvect(:,i,j,k) = dot_product( normal(:), matmul(sigma(:,:,i,j,k),normal(:)) ) * b_displl(:,i,j,k) &
                       + dot_product( normal(:), matmul(b_sigma(:,:,i,j,k),normal(:)) ) * displl(:,i,j,k)
            ! important: take only the surface part of the Kdvect
            Kdvect(:,i,j,k) = Kdvect(:,i,j,k) - normal(:) * dot_product(Kdvect(:,i,j,k),normal(:))
          endif


          ! ----- kernel contributions from all boundaries (S-S and F-S)-----
          if (k == k_disc) then
            b_kl(i,j,ispec2D) = rhostore(i,j,k,ispec) * dot_product(b_displl(:,i,j,k),accell(:,i,j,k)) &
                       + dot_product(b_dsdx(1,:,i,j,k),sigma(1,:,i,j,k)) &
                       + dot_product(b_dsdx(2,:,i,j,k),sigma(2,:,i,j,k)) &
                       + dot_product(b_dsdx(3,:,i,j,k),sigma(3,:,i,j,k)) &
                       - dot_product( matmul(b_dsdx(:,:,i,j,k), normal(:)), matmul(sigma(:,:,i,j,k),normal(:)) ) &
                       - dot_product( matmul(dsdx(:,:,i,j,k), normal(:)), matmul(b_sigma(:,:,i,j,k),normal(:)) )
          endif

        enddo
      enddo
    enddo

    ! ---- compute surface gradient of K_h for the surface element ----
    if (fluid_solid_boundary) then

      k = k_disc
      do j = 1, NGLLY
        do i = 1, NGLLX

          normal(:) = normal_disc(:,i,j,ispec2D)

          xixl = xix(i,j,k,ispec)
          xiyl = xiy(i,j,k,ispec)
          xizl = xiz(i,j,k,ispec)
          etaxl = etax(i,j,k,ispec)
          etayl = etay(i,j,k,ispec)
          etazl = etaz(i,j,k,ispec)
          gammaxl = gammax(i,j,k,ispec)
          gammayl = gammay(i,j,k,ispec)
          gammazl = gammaz(i,j,k,ispec)

          ! ----- gradient of vector boundary kernel K_h------
          temp1(:) = matmul(Kdvect(:,:,j,k), hprime_xx(i,:))
          temp2(:) = matmul(Kdvect(:,i,:,k), hprime_yy(j,:))
          temp3(:) = matmul(Kdvect(:,i,j,:), hprime_zz(k,:))

          dKdx(1,1) = xixl*temp1(1) + etaxl*temp2(1) + gammaxl*temp3(1)
          dKdx(1,2) = xiyl*temp1(1) + etayl*temp2(1) + gammayl*temp3(1)
          dKdx(1,3) = xizl*temp1(1) + etazl*temp2(1) + gammazl*temp3(1)

          dKdx(2,1) = xixl*temp1(2) + etaxl*temp2(2) + gammaxl*temp3(2)
          dKdx(2,2) = xiyl*temp1(2) + etayl*temp2(2) + gammayl*temp3(2)
          dKdx(2,3) = xizl*temp1(2) + etazl*temp2(2) + gammazl*temp3(2)

          dKdx(3,1) = xixl*temp1(3) + etaxl*temp2(3) + gammaxl*temp3(3)
          dKdx(3,2) = xiyl*temp1(3) + etayl*temp2(3) + gammayl*temp3(3)
          dKdx(3,3) = xizl*temp1(3) + etazl*temp2(3) + gammazl*temp3(3)

          ! ----- extra boundary kernel contribution for F-S ------
          b_kl_2(i,j) = dKdx(1,1) + dKdx(2,2) + dKdx(3,3) + &
                     dot_product( normal(:),matmul(dKdx(:,:),normal(:)) )
          enddo
        enddo

        b_kl(:,:,ispec2D) = b_kl(:,:,ispec2D) - b_kl_2(:,:)
      endif

    enddo

  end subroutine compute_boundary_kernel_depth


! ==========================================================================================

  subroutine compute_stress_from_strain(dsdx,sigma,i,j,k,ispec,iregion_code, &
           kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
           c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
           c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
           c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
           ystore,zstore,ibool,ispec_is_tiso)

  use constants_solver

  implicit none

  real(kind=CUSTOM_REAL) :: dsdx(NDIM,NDIM), sigma(NDIM,NDIM)
  integer i, j, k, ispec, iregion_code
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,*) :: kappavstore,muvstore, &
        kappahstore,muhstore,eta_anisostore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,*) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store
  real(kind=CUSTOM_REAL), dimension(*) :: ystore,zstore
  integer, dimension(NGLLX,NGLLY,NGLLZ,*) :: ibool
  logical, dimension(*) :: ispec_is_tiso

! --- local variables ---
  real(kind=CUSTOM_REAL) :: duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) :: duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl
  real(kind=CUSTOM_REAL) :: duxdxl,duydyl,duzdzl,duxdxl_plus_duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56
  real(kind=CUSTOM_REAL) kappal,mul,kappavl,kappahl,muvl,muhl,lambdal,lambdalplus2mul
  real(kind=CUSTOM_REAL) rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq,&
             cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
             costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
             sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta
  real(kind=CUSTOM_REAL) two_rhovpvsq,two_rhovphsq,two_rhovsvsq,two_rhovshsq
  real(kind=CUSTOM_REAL) four_rhovpvsq,four_rhovphsq,four_rhovsvsq,four_rhovshsq
  real(kind=CUSTOM_REAL) twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) two_eta_aniso,four_eta_aniso,six_eta_aniso

  integer :: iglob

  ! --- precompute sum ---

  duxdxl_plus_duydyl = dsdx(1,1) + dsdx(2,2)
  duxdxl_plus_duzdzl = dsdx(1,1) + dsdx(3,3)
  duydyl_plus_duzdzl = dsdx(2,2) + dsdx(3,3)
  duxdyl_plus_duydxl = dsdx(1,2) + dsdx(2,1)
  duzdxl_plus_duxdzl = dsdx(3,1) + dsdx(1,3)
  duzdyl_plus_duydzl = dsdx(3,2) + dsdx(2,3)
  duxdxl = dsdx(1,1)
  duydyl = dsdx(2,2)
  duzdzl = dsdx(3,3)

  ! ----------------- mantle-----------------------

  if (iregion_code == IREGION_CRUST_MANTLE) then

    if (ANISOTROPIC_3D_MANTLE_VAL) then

      c11 = c11store(i,j,k,ispec)
      c12 = c12store(i,j,k,ispec)
      c13 = c13store(i,j,k,ispec)
      c14 = c14store(i,j,k,ispec)
      c15 = c15store(i,j,k,ispec)
      c16 = c16store(i,j,k,ispec)
      c22 = c22store(i,j,k,ispec)
      c23 = c23store(i,j,k,ispec)
      c24 = c24store(i,j,k,ispec)
      c25 = c25store(i,j,k,ispec)
      c26 = c26store(i,j,k,ispec)
      c33 = c33store(i,j,k,ispec)
      c34 = c34store(i,j,k,ispec)
      c35 = c35store(i,j,k,ispec)
      c36 = c36store(i,j,k,ispec)
      c44 = c44store(i,j,k,ispec)
      c45 = c45store(i,j,k,ispec)
      c46 = c46store(i,j,k,ispec)
      c55 = c55store(i,j,k,ispec)
      c56 = c56store(i,j,k,ispec)
      c66 = c66store(i,j,k,ispec)

     sigma(1,1) = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
               c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

     sigma(2,2) = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
               c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

     sigma(3,3) = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
               c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

     sigma(1,2) = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
               c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

     sigma(1,3) = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
               c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

     sigma(2,3) = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
               c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

   else if (.not. ispec_is_tiso(ispec)) then

     ! isotropic elements

     kappal = kappavstore(i,j,k,ispec)
     mul = muvstore(i,j,k,ispec)

     lambdalplus2mul = kappal + FOUR_THIRDS * mul
     lambdal = lambdalplus2mul - 2.*mul

     sigma(1,1) = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
     sigma(2,2) = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
     sigma(3,3) = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

     sigma(1,2) = mul*duxdyl_plus_duydxl
     sigma(1,3) = mul*duzdxl_plus_duxdzl
     sigma(2,3) = mul*duzdyl_plus_duydzl

   else

     ! transverse isotropic elements

     kappavl = kappavstore(i,j,k,ispec)
     muvl = muvstore(i,j,k,ispec)

     kappahl = kappahstore(i,j,k,ispec)
     muhl = muhstore(i,j,k,ispec)

     rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
     rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A

     rhovsvsq = muvl  !!! that is L
     rhovshsq = muhl  !!! that is N

     eta_aniso = eta_anisostore(i,j,k,ispec)  !!! that is  F / (A - 2 L)

     ! use mesh coordinates to get theta and phi
     ! ystore and zstore contain theta and phi

     iglob = ibool(i,j,k,ispec)
     theta = ystore(iglob)
     phi = zstore(iglob)

     costheta = cos(theta)
     sintheta = sin(theta)
     cosphi = cos(phi)
     sinphi = sin(phi)

     costhetasq = costheta * costheta
     sinthetasq = sintheta * sintheta
     cosphisq = cosphi * cosphi
     sinphisq = sinphi * sinphi

     costhetafour = costhetasq * costhetasq
     sinthetafour = sinthetasq * sinthetasq
     cosphifour = cosphisq * cosphisq
     sinphifour = sinphisq * sinphisq

     costwotheta = cos(2.*theta)
     sintwotheta = sin(2.*theta)
     costwophi = cos(2.*phi)
     sintwophi = sin(2.*phi)

     cosfourtheta = cos(4.*theta)
     cosfourphi = cos(4.*phi)

     costwothetasq = costwotheta * costwotheta

     costwophisq = costwophi * costwophi
     sintwophisq = sintwophi * sintwophi

     etaminone = eta_aniso - 1.
     twoetaminone = 2. * eta_aniso - 1.

     ! precompute some products to reduce the CPU time

     two_eta_aniso = 2.*eta_aniso
     four_eta_aniso = 4.*eta_aniso
     six_eta_aniso = 6.*eta_aniso

     two_rhovpvsq = 2.*rhovpvsq
     two_rhovphsq = 2.*rhovphsq
     two_rhovsvsq = 2.*rhovsvsq
     two_rhovshsq = 2.*rhovshsq

     four_rhovpvsq = 4.*rhovpvsq
     four_rhovphsq = 4.*rhovphsq
     four_rhovsvsq = 4.*rhovsvsq
     four_rhovshsq = 4.*rhovshsq

     ! the 21 anisotropic coefficients computed using Mathematica

     c11 = rhovphsq*sinphifour + 2.*cosphisq*sinphisq* &
                (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                sinthetasq) + cosphifour* &
                (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                costhetasq*sinthetasq + rhovpvsq*sinthetafour)

     c12 = ((rhovphsq - two_rhovshsq)*(3. + cosfourphi)*costhetasq)/4. - &
                four_rhovshsq*cosphisq*costhetasq*sinphisq + &
                (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. + &
                eta_aniso*(rhovphsq - two_rhovsvsq)*(cosphifour + &
                2.*cosphisq*costhetasq*sinphisq + sinphifour)*sinthetasq + &
                rhovpvsq*cosphisq*sinphisq*sinthetafour - &
                rhovsvsq*sintwophisq*sinthetafour

     c13 = (cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - &
                12.*eta_aniso*rhovsvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - &
                four_eta_aniso*rhovsvsq)*cosfourtheta))/8. + &
                sinphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
                (rhovphsq - two_rhovshsq)*sinthetasq)

     c14 = costheta*sinphi*((cosphisq* &
                (-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
                (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
                (etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))*sinphisq)* sintheta

     c15 = cosphi*costheta*((cosphisq* (-rhovphsq + rhovpvsq + &
                (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                costwotheta))/2. + etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sintheta

     c16 = (cosphi*sinphi*(cosphisq* (-rhovphsq + rhovpvsq + &
                (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                four_eta_aniso*rhovsvsq)*costwotheta) + &
                2.*etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sinthetasq)/2.

     c22 = rhovphsq*cosphifour + 2.*cosphisq*sinphisq* &
                (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                sinthetasq) + sinphifour* &
                (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
                costhetasq*sinthetasq + rhovpvsq*sinthetafour)

     c23 = ((rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - 12.*eta_aniso*rhovsvsq + &
                (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                cosfourtheta)*sinphisq)/8. + &
                cosphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
                (rhovphsq - two_rhovshsq)*sinthetasq)

     c24 = costheta*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
                ((-rhovphsq + rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + &
                four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

     c25 = cosphi*costheta*((etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))* &
                cosphisq + ((-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
                (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

     c26 = (cosphi*sinphi*(2.*etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
                (-rhovphsq + rhovpvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
                four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*sinthetasq)/2.

     c33 = rhovpvsq*costhetafour + 2.*(eta_aniso*(rhovphsq - two_rhovsvsq) + two_rhovsvsq)* &
                costhetasq*sinthetasq + rhovphsq*sinthetafour

     c34 = -((rhovphsq - rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq &
                - four_eta_aniso*rhovsvsq)*costwotheta)*sinphi*sintwotheta)/4.

     c35 = -(cosphi*(rhovphsq - rhovpvsq + &
                (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                costwotheta)*sintwotheta)/4.

     c36 = -((rhovphsq - rhovpvsq - four_rhovshsq + four_rhovsvsq + &
                (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
                costwotheta)*sintwophi*sinthetasq)/4.

     c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
                sinphisq*(rhovsvsq*costwothetasq + &
                (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

     c45 = ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                four_eta_aniso*rhovsvsq + (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + &
                4.*etaminone*rhovsvsq)*costwotheta)*sintwophi*sinthetasq)/4.

     c46 = -(cosphi*costheta*((rhovshsq - rhovsvsq)*cosphisq - &
                ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
                four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)* sintheta)

     c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
                cosphisq*(rhovsvsq*costwothetasq + &
                (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

     c56 = costheta*sinphi*((cosphisq* &
                (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
                four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
                four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
                (-rhovshsq + rhovsvsq)*sinphisq)*sintheta

     c66 = rhovshsq*costwophisq*costhetasq - &
                2.*(rhovphsq - two_rhovshsq)*cosphisq*costhetasq*sinphisq + &
                (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. - &
                (rhovsvsq*(-6. - 2.*cosfourphi + cos(4.*phi - 2.*theta) - 2.*costwotheta + &
                cos(2.*(2.*phi + theta)))*sinthetasq)/8. + &
                rhovpvsq*cosphisq*sinphisq*sinthetafour - &
                (eta_aniso*(rhovphsq - two_rhovsvsq)*sintwophisq*sinthetafour)/2.

     ! general expression of stress tensor for full Cijkl with 21 coefficients

     sigma(1,1) = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
                c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

     sigma(2,2) = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
                c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

     sigma(3,3) = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
                c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

     sigma(1,2) = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
                c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

     sigma(1,3) = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
                c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

     sigma(2,3) = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
                c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

   endif ! end of test whether isotropic or anisotropic element for the mantle

! ------------------- outer  core --------------------------

 else if (iregion_code == IREGION_OUTER_CORE) then

   kappal = kappavstore(i,j,k,ispec)
   duxdxl_plus_duydyl_plus_duzdzl = duxdxl+duydyl_plus_duzdzl

   sigma(1,1) = kappal * duxdxl_plus_duydyl_plus_duzdzl
   sigma(2,2) = sigma(1,1)
   sigma(3,3) = sigma(1,1)

   sigma(1,2) = 0
   sigma(1,3) = 0
   sigma(2,3) = 0

! ------------------ inner core -------------------------

 else if (iregion_code == IREGION_INNER_CORE) then

   if (ANISOTROPIC_INNER_CORE_VAL) then

! elastic tensor for hexagonal symmetry in reduced notation:
!
!      c11 c12 c13  0   0        0
!      c12 c11 c13  0   0        0
!      c13 c13 c33  0   0        0
!       0   0   0  c44  0        0
!       0   0   0   0  c44       0
!       0   0   0   0   0  (c11-c12)/2
!
!       in terms of the A, C, L, N and F of Love (1927):
!
!       c11 = A
!       c12 = A-2N
!       c13 = F
!       c33 = C
!       c44 = L
! notice this is already in global coordinates

     c11 = c11store(i,j,k,ispec)
     c12 = c12store(i,j,k,ispec)
     c13 = c13store(i,j,k,ispec)
     c33 = c33store(i,j,k,ispec)
     c44 = c44store(i,j,k,ispec)

     sigma(1,1) = c11*duxdxl + c12*duydyl + c13*duzdzl
     sigma(2,2) = c12*duxdxl + c11*duydyl + c13*duzdzl
     sigma(3,3) = c13*duxdxl + c13*duydyl + c33*duzdzl
     sigma(1,2) = 0.5*(c11-c12)*duxdyl_plus_duydxl
     sigma(1,3) = c44*duzdxl_plus_duxdzl
     sigma(2,3) = c44*duzdyl_plus_duydzl
   else

! inner core with no anisotropy, use kappav and muv for instance
! layer with no anisotropy, use kappav and muv for instance
     kappal = kappavstore(i,j,k,ispec)
     mul = muvstore(i,j,k,ispec)

     lambdalplus2mul = kappal + FOUR_THIRDS * mul
     lambdal = lambdalplus2mul - 2.*mul

! compute stress sigma

     sigma(1,1) = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
     sigma(2,2) = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
     sigma(3,3) = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

     sigma(1,2) = mul*duxdyl_plus_duydxl
     sigma(1,3) = mul*duzdxl_plus_duxdzl
     sigma(2,3) = mul*duzdyl_plus_duydzl

   endif

 endif

! define symmetric components of sigma for gravity
  sigma(2,1) = sigma(1,2)
  sigma(3,1) = sigma(1,3)
  sigma(3,2) = sigma(2,3)

end subroutine compute_stress_from_strain
