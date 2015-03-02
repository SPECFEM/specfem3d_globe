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

  subroutine compute_kernels()

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore
  use specfem_par_movie

  implicit none

  ! crust mantle
  call compute_kernels_crust_mantle()

  ! only information to compute the crust_mantle kernels was saved to disk
  if (EXACT_UNDOING_TO_DISK) return

  ! outer core
  call compute_kernels_outer_core(vector_displ_outer_core,vector_accel_outer_core,b_vector_displ_outer_core, &
              displ_outer_core,accel_outer_core,b_displ_outer_core,b_accel_outer_core, &
              rhostore_outer_core,kappavstore_outer_core, &
              rho_kl_outer_core,alpha_kl_outer_core,beta_kl_outer_core, &
              xix_outer_core, xiy_outer_core, xiz_outer_core, etax_outer_core, etay_outer_core, etaz_outer_core, &
              gammax_outer_core, gammay_outer_core, gammaz_outer_core, ibool_outer_core, &
              nspec_beta_kl_outer_core,deviatoric_outercore)

  ! inner core
  call compute_kernels_inner_core()

  ! NOISE TOMOGRAPHY --- source strength kernel
  if (NOISE_TOMOGRAPHY == 3 ) call compute_kernels_strength_noise()

  ! boundary kernels
  if (SAVE_BOUNDARY_MESH ) call compute_boundary_kernels()

  ! approximate hessian
  if (APPROXIMATE_HESS_KL ) call compute_kernels_hessian()

  end subroutine compute_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_crust_mantle()

  use constants_solver

  use specfem_par,only: deltat,GPU_MODE,Mesh_pointer,ANISOTROPIC_KL,UNDO_ATTENUATION, &
    hprime_xx,hprime_xxT,hprime_yy,hprime_zz

  use specfem_par_crustmantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(21) :: prod
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: b_epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: b_eps_trace_over_3_loc_matrix

  integer :: i,j,k,ispec,iglob

  if (.not. GPU_MODE) then

    ! on CPU
    ! crust_mantle
    do ispec = 1, NSPEC_CRUST_MANTLE

      ! simulations with UNDO_ATTENUATION save as much memory as possible;
      ! backward/reconstructed wavefield strain will be re-computed locally here
      if (UNDO_ATTENUATION) then
        if (USE_DEVILLE_PRODUCTS_VAL) then
          call compute_element_strain_undoatt_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE,&
                                                  b_displ_crust_mantle,ibool_crust_mantle, &
                                                  hprime_xx,hprime_xxT,&
                                                  xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                                  etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                                  gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle,&
                                                  b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)

        else
          call compute_element_strain_undoatt_noDev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                                    b_displ_crust_mantle, &
                                                    hprime_xx,hprime_yy,hprime_zz,ibool_crust_mantle, &
                                                    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                                    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                                    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                                    b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)
        endif
      else
        ! backward/reconstructed strain arrays
        b_eps_trace_over_3_loc_matrix(:,:,:) = b_eps_trace_over_3_crust_mantle(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,1) = b_epsilondev_xx_crust_mantle(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,2) = b_epsilondev_yy_crust_mantle(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,3) = b_epsilondev_xy_crust_mantle(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,4) = b_epsilondev_xz_crust_mantle(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,5) = b_epsilondev_yz_crust_mantle(:,:,:,ispec)
      endif

      ! For anisotropic kernels
      if (ANISOTROPIC_KL) then

        ! computes fully anisotropic kernel cijkl_kl
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool_crust_mantle(i,j,k,ispec)

              ! density kernel: see e.g. Tromp et al.(2005), equation (14)
              !                         b_displ_crust_mantle is the backward/reconstructed wavefield, that is s(x,t) in eq. (14),
              !                         accel_crust_mantle is the adjoint wavefield, that corresponds to s_dagger(x,T-t)
              !
              !                         note with respect to eq. (14) the second time derivative is applied to the
              !                         adjoint wavefield here rather than the backward/reconstructed wavefield.
              !                         this is a valid operation and the resultant kernel identical to the eq. (14).
              !
              !                         reason for this is that the adjoint wavefield is in general smoother
              !                         since the adjoint sources normally are obtained for filtered traces.
              !                         numerically, the time derivative by a finite-difference scheme should
              !                         behave better for smoother wavefields, thus containing less numerical artifacts.
              rho_kl_crust_mantle(i,j,k,ispec) =  rho_kl_crust_mantle(i,j,k,ispec) &
                 + deltat * (accel_crust_mantle(1,iglob) * b_displ_crust_mantle(1,iglob) &
                           + accel_crust_mantle(2,iglob) * b_displ_crust_mantle(2,iglob) &
                           + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob) )

              ! fully anisotropic kernel
              ! temporary arrays
              epsilondev_loc(1) = epsilondev_xx_crust_mantle(i,j,k,ispec)
              epsilondev_loc(2) = epsilondev_yy_crust_mantle(i,j,k,ispec)
              epsilondev_loc(3) = epsilondev_xy_crust_mantle(i,j,k,ispec)
              epsilondev_loc(4) = epsilondev_xz_crust_mantle(i,j,k,ispec)
              epsilondev_loc(5) = epsilondev_yz_crust_mantle(i,j,k,ispec)

              b_epsilondev_loc(1) = b_epsilondev_loc_matrix(i,j,k,1)
              b_epsilondev_loc(2) = b_epsilondev_loc_matrix(i,j,k,2)
              b_epsilondev_loc(3) = b_epsilondev_loc_matrix(i,j,k,3)
              b_epsilondev_loc(4) = b_epsilondev_loc_matrix(i,j,k,4)
              b_epsilondev_loc(5) = b_epsilondev_loc_matrix(i,j,k,5)

              call compute_strain_product(prod,eps_trace_over_3_crust_mantle(i,j,k,ispec),epsilondev_loc, &
                                          b_eps_trace_over_3_loc_matrix(i,j,k),b_epsilondev_loc)

              cijkl_kl_crust_mantle(:,i,j,k,ispec) = cijkl_kl_crust_mantle(:,i,j,k,ispec) + deltat * prod(:)

            enddo
          enddo
        enddo

      else

        ! isotropic kernels

        ! if EXACT_UNDOING_TO_DISK is set, currently the rho and beta kernels below will be computed but will be equal
        ! to zero everywhere because the backward field is then not computed, and only the alpha kernel
        ! will be computed correctly based on the values of the forward run saved to disk in the first run
        ! and read back at the beginning of this routine

        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              iglob = ibool_crust_mantle(i,j,k,ispec)

              ! density kernel: see e.g. Tromp et al.(2005), equation (14)
              !                         b_displ_crust_mantle is the backward/reconstructed wavefield, that is s(x,t) in eq. (14),
              !                         accel_crust_mantle is the adjoint wavefield, that corresponds to s_dagger(x,T-t)
              !
              !                         note with respect to eq. (14) the second time derivative is applied to the
              !                         adjoint wavefield here rather than the backward/reconstructed wavefield.
              !                         this is a valid operation and the resultant kernel identical to the eq. (14).
              !
              !                         reason for this is that the adjoint wavefield is in general smoother
              !                         since the adjoint sources normally are obtained for filtered traces.
              !                         numerically, the time derivative by a finite-difference scheme should
              !                         behave better for smoother wavefields, thus containing less numerical artifacts.
              rho_kl_crust_mantle(i,j,k,ispec) =  rho_kl_crust_mantle(i,j,k,ispec) &
                 + deltat * (accel_crust_mantle(1,iglob) * b_displ_crust_mantle(1,iglob) &
                           + accel_crust_mantle(2,iglob) * b_displ_crust_mantle(2,iglob) &
                           + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob) )

              ! isotropic kernels
              ! temporary arrays
              epsilondev_loc(1) = epsilondev_xx_crust_mantle(i,j,k,ispec)
              epsilondev_loc(2) = epsilondev_yy_crust_mantle(i,j,k,ispec)
              epsilondev_loc(3) = epsilondev_xy_crust_mantle(i,j,k,ispec)
              epsilondev_loc(4) = epsilondev_xz_crust_mantle(i,j,k,ispec)
              epsilondev_loc(5) = epsilondev_yz_crust_mantle(i,j,k,ispec)

              b_epsilondev_loc(1) = b_epsilondev_loc_matrix(i,j,k,1)
              b_epsilondev_loc(2) = b_epsilondev_loc_matrix(i,j,k,2)
              b_epsilondev_loc(3) = b_epsilondev_loc_matrix(i,j,k,3)
              b_epsilondev_loc(4) = b_epsilondev_loc_matrix(i,j,k,4)
              b_epsilondev_loc(5) = b_epsilondev_loc_matrix(i,j,k,5)

              ! kernel for shear modulus, see e.g. Tromp et al. (2005), equation (17)
              ! note: multiplication with 2*mu(x) will be done after the time loop
              beta_kl_crust_mantle(i,j,k,ispec) =  beta_kl_crust_mantle(i,j,k,ispec) &
                 + deltat * &
                 (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
                  + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
                  + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) + &
                  epsilondev_loc(5)*b_epsilondev_loc(5)) )

              ! kernel for bulk modulus, see e.g. Tromp et al. (2005), equation (18)
              ! note: multiplication with kappa(x) will be done after the time loop
              alpha_kl_crust_mantle(i,j,k,ispec) = alpha_kl_crust_mantle(i,j,k,ispec) &
                 + deltat * (9 * eps_trace_over_3_crust_mantle(i,j,k,ispec) &
                               * b_eps_trace_over_3_loc_matrix(i,j,k))

            enddo
          enddo
        enddo

      endif ! ANISOTROPIC_KL

    enddo

  else
    ! updates kernel contribution on GPU

    ! computes contribution to density and isotropic/anisotropic kernels
    call compute_kernels_cm_gpu(Mesh_pointer,deltat)

  endif

  end subroutine compute_kernels_crust_mantle

!
!-------------------------------------------------------------------------------------------------
!

!! DK DK put the list of parameters back here to avoid a warning / error from the gfortran compiler
!! DK DK about undefined behavior when aggressive loop vectorization is used by the compiler
  subroutine compute_kernels_outer_core(vector_displ_outer_core,vector_accel_outer_core,b_vector_displ_outer_core, &
              displ_outer_core,accel_outer_core,b_displ_outer_core,b_accel_outer_core, &
              rhostore_outer_core,kappavstore_outer_core, &
              rho_kl_outer_core,alpha_kl_outer_core,beta_kl_outer_core, &
              xix_outer_core, xiy_outer_core, xiz_outer_core, etax_outer_core, etay_outer_core, etaz_outer_core, &
              gammax_outer_core, gammay_outer_core, gammaz_outer_core, ibool_outer_core, &
              nspec_beta_kl_outer_core,deviatoric_outercore)

  use constants_solver
  use specfem_par,only: deltat,hprime_xx,hprime_yy,hprime_zz,myrank
  use specfem_par,only: GPU_MODE,Mesh_pointer

  implicit none

  integer :: nspec_beta_kl_outer_core

  logical :: deviatoric_outercore

  ! velocity potential
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: displ_outer_core,accel_outer_core

  ! ADJOINT
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: b_displ_outer_core,b_accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_OUTER_CORE_ADJOINT) :: &
    vector_accel_outer_core,vector_displ_outer_core,b_vector_displ_outer_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
    xix_outer_core,xiy_outer_core,xiz_outer_core,&
    etax_outer_core,etay_outer_core,etaz_outer_core, &
    gammax_outer_core,gammay_outer_core,gammaz_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: rhostore_outer_core,kappavstore_outer_core

  ! adjoint kernels
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ADJOINT) :: rho_kl_outer_core,alpha_kl_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_beta_kl_outer_core) :: beta_kl_outer_core

  ! local parameters
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,kappal
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) :: tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) :: tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc
  real(kind=CUSTOM_REAL) :: div_displ,b_div_displ
  real(kind=CUSTOM_REAL) :: gradx,grady,gradz
  integer :: i,j,k,l,ispec,iglob

  logical,dimension(:),allocatable :: mask_ibool
  integer :: ier

  ! outer_core -- compute the actual displacement and acceleration (NDIM,NGLOBMAX_OUTER_CORE)

  if (.not. GPU_MODE) then
    ! on CPU

    allocate(mask_ibool(NGLOB_OUTER_CORE),stat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error allocating mask_ibool array in routine compute_boundary_kernels()')
    mask_ibool(:) = .false.

    ! pre-calculates gradients in outer core on CPU
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

            ! global index
            iglob = ibool_outer_core(i,j,k,ispec)

            ! only calculate the gradients once for shared nodes
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

              ! assigns gradient field on global points
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

              ! assigns gradient field on global points
              vector_accel_outer_core(1,iglob) = gradx
              vector_accel_outer_core(2,iglob) = grady
              vector_accel_outer_core(3,iglob) = gradz

              ! calculates gradient grad(displ) (also needed for boundary kernels)
              if (SAVE_BOUNDARY_MESH .or. deviatoric_outercore) then
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

            endif ! mask_ibool

          enddo
        enddo
      enddo
    enddo

    ! frees memory
    deallocate(mask_ibool)

    ! acoustic kernels
    do ispec = 1, NSPEC_OUTER_CORE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX

            ! global index
            iglob = ibool_outer_core(i,j,k,ispec)

            gradx = vector_accel_outer_core(1,iglob) * b_vector_displ_outer_core(1,iglob)
            grady = vector_accel_outer_core(2,iglob) * b_vector_displ_outer_core(2,iglob)
            gradz = vector_accel_outer_core(3,iglob) * b_vector_displ_outer_core(3,iglob)

            ! density kernel
            ! note: we replace dot_product() with an unrolled expression, otherwise most compilers
            !       will try to vectorize this rather than the outer loop, resulting in a much slower code
            rho_kl_outer_core(i,j,k,ispec) = rho_kl_outer_core(i,j,k,ispec) + deltat * (gradx + grady + gradz)

            ! bulk modulus kernel
            kappal = rhostore_outer_core(i,j,k,ispec)/kappavstore_outer_core(i,j,k,ispec)
            div_displ =  kappal * accel_outer_core(iglob)
            b_div_displ =  kappal * b_accel_outer_core(iglob)

            alpha_kl_outer_core(i,j,k,ispec) = alpha_kl_outer_core(i,j,k,ispec) + deltat * div_displ * b_div_displ

          enddo
        enddo
      enddo
    enddo

    !deviatoric kernel check
    if (deviatoric_outercore) then

      do ispec = 1, NSPEC_OUTER_CORE
        do k = 1, NGLLZ
          do j = 1, NGLLY
            do i = 1, NGLLX
              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempx3l = 0._CUSTOM_REAL

              tempy1l = 0._CUSTOM_REAL
              tempy2l = 0._CUSTOM_REAL
              tempy3l = 0._CUSTOM_REAL

              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              tempz3l = 0._CUSTOM_REAL

              ! assumes NGLLX = NGLLY = NGLLZ
              do l = 1,NGLLX
                tempx1l = tempx1l + b_vector_displ_outer_core(1,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
                tempy1l = tempy1l + b_vector_displ_outer_core(2,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
                tempz1l = tempz1l + b_vector_displ_outer_core(3,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)

                tempx2l = tempx2l + b_vector_displ_outer_core(1,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
                tempy2l = tempy2l + b_vector_displ_outer_core(2,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
                tempz2l = tempz2l + b_vector_displ_outer_core(3,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)

                tempx3l = tempx3l +  b_vector_displ_outer_core(1,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
                tempy3l = tempy3l +  b_vector_displ_outer_core(2,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
                tempz3l = tempz3l +  b_vector_displ_outer_core(3,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              enddo

              !deviatoric strain
              b_epsilondev_loc(1) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              b_epsilondev_loc(2) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              b_epsilondev_loc(3) = 0.5*( xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l &
                                        + xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l ) &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              b_epsilondev_loc(4) = 0.5*( xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l &
                                        + xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l ) &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              b_epsilondev_loc(5) = 0.5*( xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l &
                                        + xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l ) &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              tempx1l = 0._CUSTOM_REAL
              tempx2l = 0._CUSTOM_REAL
              tempx3l = 0._CUSTOM_REAL

              tempy1l = 0._CUSTOM_REAL
              tempy2l = 0._CUSTOM_REAL
              tempy3l = 0._CUSTOM_REAL

              tempz1l = 0._CUSTOM_REAL
              tempz2l = 0._CUSTOM_REAL
              tempz3l = 0._CUSTOM_REAL

              ! assumes NGLLX = NGLLY = NGLLZ
              do l = 1,NGLLX
                tempx1l = tempx1l + vector_displ_outer_core(1,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
                tempy1l = tempy1l + vector_displ_outer_core(2,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)
                tempz1l = tempz1l + vector_displ_outer_core(3,ibool_outer_core(l,j,k,ispec)) * hprime_xx(i,l)

                tempx2l = tempx2l + vector_displ_outer_core(1,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
                tempy2l = tempy2l + vector_displ_outer_core(2,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)
                tempz2l = tempz2l + vector_displ_outer_core(3,ibool_outer_core(i,l,k,ispec)) * hprime_yy(j,l)

                tempx3l = tempx3l + vector_displ_outer_core(1,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
                tempy3l = tempy3l + vector_displ_outer_core(2,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
                tempz3l = tempz3l + vector_displ_outer_core(3,ibool_outer_core(i,j,l,ispec)) * hprime_zz(k,l)
              enddo

              !deviatoric strain
              epsilondev_loc(1) = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              epsilondev_loc(2) = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              epsilondev_loc(3) = 0.5*( xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l &
                                        + xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l ) &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              epsilondev_loc(4) = 0.5*( xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l &
                                        + xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l ) &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              epsilondev_loc(5) = 0.5*( xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l &
                                        + xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l ) &
                  - ONE_THIRD* (xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l &
                                + xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l &
                                + xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l )

              beta_kl_outer_core(i,j,k,ispec) =  beta_kl_outer_core(i,j,k,ispec) &
                 + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
                 + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
                 + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + epsilondev_loc(4)*b_epsilondev_loc(4) + &
                  epsilondev_loc(5)*b_epsilondev_loc(5)) )

            enddo
          enddo
        enddo
      enddo

    endif !deviatoric kernel check

  else
    ! updates kernel contribution on GPU
    if (deviatoric_outercore ) call exit_mpi(myrank,'deviatoric kernel on GPU not supported yet')

    ! computes contribution to density and bulk modulus kernel
    call compute_kernels_oc_gpu(Mesh_pointer,deltat)

  endif

  end subroutine compute_kernels_outer_core

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_inner_core()

  use constants_solver

  use specfem_par,only: deltat,GPU_MODE,Mesh_pointer,UNDO_ATTENUATION, &
    hprime_xx,hprime_xxT,hprime_yy,hprime_zz

  use specfem_par_innercore

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(5) :: b_epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5) :: epsilondev_loc

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: b_epsilondev_loc_matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: b_eps_trace_over_3_loc_matrix

  integer :: i,j,k,ispec,iglob

  if (.not. GPU_MODE) then
    ! on CPU
    ! inner_core
    do ispec = 1, NSPEC_INNER_CORE

      ! gets element strain
      if (UNDO_ATTENUATION) then
        if (USE_DEVILLE_PRODUCTS_VAL) then
          call compute_element_strain_undoatt_Dev(ispec,NGLOB_inner_core,NSPEC_inner_core, &
                                                  b_displ_inner_core,ibool_inner_core, &
                                                  hprime_xx,hprime_xxT, &
                                                  xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                                  etax_inner_core,etay_inner_core,etaz_inner_core, &
                                                  gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                                  b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)
        else
          call compute_element_strain_undoatt_noDev(ispec,NGLOB_inner_core,NSPEC_inner_core, &
                                                    b_displ_inner_core, &
                                                    hprime_xx,hprime_yy,hprime_zz,ibool_inner_core, &
                                                    xix_inner_core,xiy_inner_core,xiz_inner_core, &
                                                    etax_inner_core,etay_inner_core,etaz_inner_core, &
                                                    gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
                                                    b_epsilondev_loc_matrix,b_eps_trace_over_3_loc_matrix)
        endif
      else
        ! backward/reconstructed strain arrays
        b_eps_trace_over_3_loc_matrix(:,:,:) = b_eps_trace_over_3_inner_core(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,1) = b_epsilondev_xx_inner_core(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,2) = b_epsilondev_yy_inner_core(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,3) = b_epsilondev_xy_inner_core(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,4) = b_epsilondev_xz_inner_core(:,:,:,ispec)
        b_epsilondev_loc_matrix(:,:,:,5) = b_epsilondev_yz_inner_core(:,:,:,ispec)
      endif

      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)

            rho_kl_inner_core(i,j,k,ispec) =  rho_kl_inner_core(i,j,k,ispec) &
               + deltat * (accel_inner_core(1,iglob) * b_displ_inner_core(1,iglob) + &
                           accel_inner_core(2,iglob) * b_displ_inner_core(2,iglob) + &
                           accel_inner_core(3,iglob) * b_displ_inner_core(3,iglob) )

            epsilondev_loc(1) = epsilondev_xx_inner_core(i,j,k,ispec)
            epsilondev_loc(2) = epsilondev_yy_inner_core(i,j,k,ispec)
            epsilondev_loc(3) = epsilondev_xy_inner_core(i,j,k,ispec)
            epsilondev_loc(4) = epsilondev_xz_inner_core(i,j,k,ispec)
            epsilondev_loc(5) = epsilondev_yz_inner_core(i,j,k,ispec)

            b_epsilondev_loc(1) = b_epsilondev_loc_matrix(i,j,k,1)
            b_epsilondev_loc(2) = b_epsilondev_loc_matrix(i,j,k,2)
            b_epsilondev_loc(3) = b_epsilondev_loc_matrix(i,j,k,3)
            b_epsilondev_loc(4) = b_epsilondev_loc_matrix(i,j,k,4)
            b_epsilondev_loc(5) = b_epsilondev_loc_matrix(i,j,k,5)

            beta_kl_inner_core(i,j,k,ispec) =  beta_kl_inner_core(i,j,k,ispec) &
               + deltat * (epsilondev_loc(1)*b_epsilondev_loc(1) + epsilondev_loc(2)*b_epsilondev_loc(2) &
                        + (epsilondev_loc(1)+epsilondev_loc(2)) * (b_epsilondev_loc(1)+b_epsilondev_loc(2)) &
                        + 2 * (epsilondev_loc(3)*b_epsilondev_loc(3) + &
                               epsilondev_loc(4)*b_epsilondev_loc(4) + &
                               epsilondev_loc(5)*b_epsilondev_loc(5)) )

            alpha_kl_inner_core(i,j,k,ispec) = alpha_kl_inner_core(i,j,k,ispec) &
               + deltat * (9 * eps_trace_over_3_inner_core(i,j,k,ispec) * &
                               b_eps_trace_over_3_loc_matrix(i,j,k))
          enddo
        enddo
      enddo
    enddo

  else
    ! updates kernel contribution on GPU

    ! computes contribution to density and bulk and shear modulus kernel
    call compute_kernels_ic_gpu(Mesh_pointer,deltat)

  endif

  end subroutine compute_kernels_inner_core


!
!-------------------------------------------------------------------------------------------------
!

! Subroutines to compute the kernels for the 21 elastic coefficients

  subroutine compute_strain_product(prod,eps_trace_over_3,epsdev,&
                                    b_eps_trace_over_3,b_epsdev)

  ! Purpose: compute the 21 strain products at a grid point
  ! (ispec,i,j,k fixed) and at a time t to compute then the kernels cij_kl (Voigt notation)
  ! (eq. 15 of Tromp et al., 2005)
  ! prod(1)=eps11*eps11 -> c11, prod(2)=eps11eps22 -> c12, prod(3)=eps11eps33 -> c13, ...
  ! prod(7)=eps22*eps22 -> c22, prod(8)=eps22eps33 -> c23, prod(9)=eps22eps23 -> c24, ...
  ! prod(19)=eps13*eps13 -> c55, prod(20)=eps13eps12 -> c56, prod(21)=eps12eps12 -> c66
  ! This then gives how the 21 kernels are organized
  ! For crust_mantle

  use constants

  implicit none

  real(kind=CUSTOM_REAL),dimension(21) :: prod
  real(kind=CUSTOM_REAL) :: eps_trace_over_3,b_eps_trace_over_3
  real(kind=CUSTOM_REAL),dimension(5) :: epsdev,b_epsdev

  real(kind=CUSTOM_REAL), dimension(6) :: eps,b_eps
  integer :: p,i,j

  ! Building of the local matrix of the strain tensor
  ! for the adjoint field and the regular backward field
  eps(1:2) = epsdev(1:2)+eps_trace_over_3           !eps11 et eps22
  eps(3) = -(eps(1)+eps(2))+3*eps_trace_over_3     !eps33
  eps(4) = epsdev(5)                                !eps23
  eps(5) = epsdev(4)                                !eps13
  eps(6) = epsdev(3)                                !eps12

  b_eps(1:2) = b_epsdev(1:2)+b_eps_trace_over_3
  b_eps(3) = -(b_eps(1)+b_eps(2))+3*b_eps_trace_over_3
  b_eps(4) = b_epsdev(5)
  b_eps(5) = b_epsdev(4)
  b_eps(6) = b_epsdev(3)

  ! Computing the 21 strain products without assuming eps(i)*b_eps(j) = eps(j)*b_eps(i)
  p = 1
  do i = 1,6
    do j = i,6
      prod(p) = eps(i)*b_eps(j)
      if (j>i) then
        prod(p) = prod(p)+eps(j)*b_eps(i)
        if (j>3 .and. i<4) prod(p) = prod(p) * 2.0_CUSTOM_REAL
      endif
      if (i>3) prod(p) = prod(p) * 4.0_CUSTOM_REAL
      p = p+1
    enddo
  enddo

  end subroutine compute_strain_product

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_kernels_hessian()

  use constants_solver
  use specfem_par,only: deltat
  use specfem_par,only: GPU_MODE,Mesh_pointer
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,iglob

  if (.not. GPU_MODE) then
    ! on CPU
    ! crust_mantle
    do ispec = 1, NSPEC_CRUST_MANTLE
      do k = 1, NGLLZ
        do j = 1, NGLLY
          do i = 1, NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)

            ! approximates hessian
            ! term with adjoint acceleration and backward/reconstructed acceleration
            hess_kl_crust_mantle(i,j,k,ispec) =  hess_kl_crust_mantle(i,j,k,ispec) &
               + deltat * (accel_crust_mantle(1,iglob) * b_accel_crust_mantle(1,iglob) &
               + accel_crust_mantle(2,iglob) * b_accel_crust_mantle(2,iglob) &
               + accel_crust_mantle(3,iglob) * b_accel_crust_mantle(3,iglob) )

          enddo
        enddo
      enddo
    enddo

  else
    ! updates kernel contribution on GPU

    ! computes contribution to density and bulk modulus kernel
    call compute_kernels_hess_gpu(Mesh_pointer,deltat)

  endif

  end subroutine compute_kernels_hessian
