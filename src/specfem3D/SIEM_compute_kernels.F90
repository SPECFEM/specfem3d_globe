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


  subroutine SIEM_compute_gravity_kernels()

  use specfem_par, only: FULL_GRAVITY,SIMULATION_TYPE,myrank

  implicit none

  ! checks if anything to do
  if (.not. FULL_GRAVITY) return
  if (SIMULATION_TYPE /= 3) return

  !debug
  if (myrank == 0) print *,'Calculate 1st gravity kernel...'
  call synchronize_all()

  call calculate_first_gravity_kernel()

  !debug
  if (myrank == 0) print *,'Calculate 2nd gravity kernel...'
  call synchronize_all()

  call calculate_second_gravity_kernel()

  !debug
  if (myrank == 0) print *,'DONE!'
  call synchronize_all()

  end subroutine SIEM_compute_gravity_kernels

!
!-------------------------------------------------------------------------------
!

  subroutine calculate_first_gravity_kernel()

  ! ------------------- FIRST GRAVITY KERNEL -------------------

  use specfem_par
  use specfem_par_crustmantle, only: ibool_crust_mantle, &
                                     xix_crust_mantle, etax_crust_mantle, gammax_crust_mantle, &
                                     xiy_crust_mantle, etay_crust_mantle, gammay_crust_mantle, &
                                     xiz_crust_mantle, etaz_crust_mantle, gammaz_crust_mantle

  use specfem_par_full_gravity, only: gdof_cm1, inode_elmt_cm, inode_map_cm, nmir_cm, nnode_cm1, &
    gknl1, &
    gravload1, pgrav1, pgrav_cm1, pgrav_cm, neq1, ndscale1, dprecon1, &
    is_active_gll,igll_active_on, &
    CG_SCALING

  use siem_math_library_mpi, only: maxvec
  use siem_poisson, only: compute_grav_kl1_load
  use siem_solver_petsc, only: petsc_set_vector1, petsc_zero_initialguess1, petsc_solve1
  use siem_solver_mpi, only: cg_solver3, diagpcg_solver3, interpolate3to5

  implicit none

  ! Local variables
  integer :: ispec,i,j,k,l,ier
  integer :: icomponent
  real(kind=CUSTOM_REAL) :: tempx1l, tempx2l, tempx3l, tempy1l, tempy2l, tempy3l, tempz1l, tempz2l, tempz3l
  real(kind=CUSTOM_REAL) :: maxload
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: gknl1_cm ! (3,NGLOB_CRUST_MANTLE)

  ! Allocate 1st grav kernel array :
  allocate(gknl1(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating gknl1'
  gknl1(:,:,:,:) = 0.0_CUSTOM_REAL

  ! temporary
  allocate(gknl1_cm(3,NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating gknl1_cm'
  gknl1_cm(:,:) = 0.0_CUSTOM_REAL

  ! Solve a modified poisson equation for each component
  do icomponent = 1,3

    !debug
    if (myrank == 0) print *,'icomponent: ', icomponent

    ! Calculate the RHS (gravload1)
    call compute_grav_kl1_load(icomponent)

    !debug
    maxload = maxvec(abs(gravload1))
    if (myrank == 0) print *,'  -- Max load: ', maxload

    if (POISSON_SOLVER == ISOLVER_BUILTIN) then
      ! builtin solver
      pgrav1(:) = 0.0_CUSTOM_REAL
      if (CG_SCALING) then
        gravload1(:) = ndscale1(:) * gravload1(:)
        call cg_solver3(myrank,neq1,pgrav1,gravload1)
        pgrav1(:) = ndscale1(:) * pgrav1(:)
      else
        call diagpcg_solver3(myrank,neq1,pgrav1,gravload1,dprecon1)
      endif
    else
      ! Petsc solver
      call petsc_set_vector1(gravload1)

      ! Need to zero guess since was previously for poisson eqn
      ! Do we need this - is zero guess set to PETSC_TRUE?
      call petsc_zero_initialguess1()

      pgrav1(:) = 0.0_CUSTOM_REAL
      ! Here we use pgrav1 as the vector the solution is put into, just to
      ! save on allocating another array. Note that we could allocate a
      ! separate temporary array but pgrav1 isnt used after this (apart from
      ! if the forward wavefield is being saved, but it shouldnt be for SIMTYPE 3)
      call petsc_solve1(pgrav1(1:))
    endif

    ! Now interpolate this component to the 5GLL setup (Crust mantle only)
    pgrav_cm1(:) = zero
    pgrav_cm1(:) = pgrav1(gdof_cm1(:))

    pgrav_cm(:) = 0.0_CUSTOM_REAL !initialise before interpolating to be safe

    call interpolate3to5(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,nnode_cm1, &
                         inode_elmt_cm,nmir_cm,inode_map_cm,is_active_gll,igll_active_on,pgrav_cm1,pgrav_cm)

    gknl1_cm(icomponent,:) = pgrav_cm(:)
    call synchronize_all()
  enddo

  if (myrank == 0) print *,'  -- Calculating divergence of gknl1_cm...'

  ! Now calculate NEGATIVE divergence of gknl1_cm to get kernel:
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          tempx1l = 0._CUSTOM_REAL;   tempx2l = 0._CUSTOM_REAL;   tempx3l = 0._CUSTOM_REAL
          tempy1l = 0._CUSTOM_REAL;   tempy2l = 0._CUSTOM_REAL;   tempy3l = 0._CUSTOM_REAL
          tempz1l = 0._CUSTOM_REAL;   tempz2l = 0._CUSTOM_REAL;   tempz3l = 0._CUSTOM_REAL
          do l = 1,NGLLX
            tempx1l = tempx1l + gknl1_cm(1,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempy1l = tempy1l + gknl1_cm(1,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempz1l = tempz1l + gknl1_cm(1,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)

            tempx2l = tempx2l + gknl1_cm(2,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempy2l = tempy2l + gknl1_cm(2,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempz2l = tempz2l + gknl1_cm(2,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)

            tempx3l = tempx3l + gknl1_cm(3,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempy3l = tempy3l + gknl1_cm(3,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempz3l = tempz3l + gknl1_cm(3,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
          enddo !l

          ! Note this gives you k_rho, not -k_rho - the final output we need is ( - rho * k_rho)
          ! So in save_kernels we will multiply this by -rho
          gknl1(i,j,k,ispec) = gknl1(i,j,k,ispec) - (xix_crust_mantle(i,j,k,ispec)    * tempx1l + &
                                                     etax_crust_mantle(i,j,k,ispec)   * tempy1l + &
                                                     gammax_crust_mantle(i,j,k,ispec) * tempz1l + &
                                                     xiy_crust_mantle(i,j,k,ispec)    * tempx2l + &
                                                     etay_crust_mantle(i,j,k,ispec)   * tempy2l + &
                                                     gammay_crust_mantle(i,j,k,ispec) * tempz2l + &
                                                     xiz_crust_mantle(i,j,k,ispec)    * tempx3l + &
                                                     etaz_crust_mantle(i,j,k,ispec)   * tempy3l + &
                                                     gammaz_crust_mantle(i,j,k,ispec) * tempz3l)
        enddo !i
      enddo !j
    enddo ! k
  enddo !ispec

  ! free temporary array
  deallocate(gknl1_cm)

  end subroutine calculate_first_gravity_kernel

!
!-------------------------------------------------------------------------------
!

  subroutine calculate_second_gravity_kernel()

  ! ------------------- SECOND GRAVITY KERNEL -------------------

  use specfem_par
  use specfem_par_crustmantle, only: ibool_crust_mantle, &
                                     xix_crust_mantle, xiy_crust_mantle, xiz_crust_mantle, &
                                     etax_crust_mantle, etay_crust_mantle, etaz_crust_mantle, &
                                     gammax_crust_mantle, gammay_crust_mantle, gammaz_crust_mantle

  use specfem_par_full_gravity, only: gdof_cm1, inode_elmt_cm, inode_map_cm, nmir_cm, nnode_cm1, &
    gknl2, &
    gravload1, pgrav1, pgrav_cm1, pgrav_cm, neq1, ndscale1, dprecon1, &
    is_active_gll,igll_active_on, &
    CG_SCALING

  use siem_math_library_mpi, only: maxvec
  use siem_poisson, only: compute_grav_kl2_load
  use siem_solver_petsc, only: petsc_set_vector1, petsc_zero_initialguess1, petsc_solve1
  use siem_solver_mpi, only: cg_solver3, diagpcg_solver3, interpolate3to5

  implicit none

  ! Local variables
  integer :: ispec,i,j,k,l,m,ier
  integer :: icomponent, jcomponent
  real(kind=CUSTOM_REAL) :: tempdiv1(3),tempdiv2(3),tempdiv3(3)
  real(kind=CUSTOM_REAL) :: maxload
  real(kind=CUSTOM_REAL),dimension(:,:,:), allocatable :: gknl2_cm
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: div_gknl2_cm

  ! Allocate 2nd gravity kernel array
  allocate(gknl2(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
  if (ier /= 0) stop 'Error allocating gknl2'
  gknl2(:,:,:,:) = 0.0_CUSTOM_REAL

  ! temporary
  allocate(gknl2_cm(3,3,NGLOB_CRUST_MANTLE), &
           div_gknl2_cm(3,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating gknl2_cm,.. arrays'
  gknl2_cm(:,:,:) = 0.0_CUSTOM_REAL
  div_gknl2_cm(:,:,:,:,:) = 0.0_CUSTOM_REAL

  ! Solve the equation 9 times!
  do icomponent = 1,3
    do jcomponent = 1,3
      !debug
      if (myrank == 0) print *,'component ', icomponent, jcomponent

      ! Calculate the RHS (gravload1)
      call compute_grav_kl2_load(icomponent, jcomponent)

      !debug
      maxload = maxvec(abs(gravload1))
      if (myrank == 0) print *,'  -- Max load: ', maxload

      if (POISSON_SOLVER == ISOLVER_BUILTIN) then
        ! builtin solver
        pgrav1(:) = 0.0_CUSTOM_REAL
        if (CG_SCALING) then
          gravload1(:) = ndscale1(:) * gravload1(:)
          call cg_solver3(myrank,neq1,pgrav1,gravload1)
          pgrav1(:) = ndscale1(:) * pgrav1(:)
        else
          call diagpcg_solver3(myrank,neq1,pgrav1,gravload1,dprecon1)
        endif
      else
        ! PETSc solver
        call petsc_set_vector1(gravload1)

        call petsc_zero_initialguess1()

        pgrav1(:) = 0.0_CUSTOM_REAL
        call petsc_solve1(pgrav1(1:))
      endif

      ! Now interpolate this component to the 5GLL setup (Crust mantle only)
      pgrav_cm1(:) = zero
      pgrav_cm1(:) = pgrav1(gdof_cm1(:))

      pgrav_cm(:)  = zero !initialise before interpolating to be safe

      call interpolate3to5(NSPEC_CRUST_MANTLE, NGLOB_CRUST_MANTLE, nnode_cm1, &
                           inode_elmt_cm, nmir_cm, inode_map_cm, is_active_gll, igll_active_on, pgrav_cm1, pgrav_cm)

      ! Store this component of the gravity kernel
      gknl2_cm(icomponent,jcomponent,:) = pgrav_cm(:)

    enddo !jcomponent
  enddo !icomponent

  ! Calculate the divergence of gknl2_cm
  ! note that in the manuscript div_gknl2_cm (a vector) is denoted by a roman b
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! Loop for each component of the resulting vector:
          do m = 1,3
            tempdiv1 = 0._CUSTOM_REAL
            tempdiv2 = 0._CUSTOM_REAL
            tempdiv3 = 0._CUSTOM_REAL

            do l = 1,NGLLX
              ! tempdiv1 is a 3 component array where each component (1,2,3) represents either \beta_{1,j}, \beta_{2,j} or \beta_{3,j}
              tempdiv1 = tempdiv1 + gknl2_cm(:,m,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
              tempdiv2 = tempdiv2 + gknl2_cm(:,m,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
              tempdiv3 = tempdiv3 + gknl2_cm(:,m,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
            enddo !l

            div_gknl2_cm(m,i,j,k,ispec) = div_gknl2_cm(m,i,j,k,ispec) + &
                                           xix_crust_mantle(i,j,k,ispec)*tempdiv1(1)    + &
                                           xiy_crust_mantle(i,j,k,ispec)*tempdiv1(2)    + &
                                           xiz_crust_mantle(i,j,k,ispec)*tempdiv1(3)    + &
                                           etax_crust_mantle(i,j,k,ispec)*tempdiv2(1)   + &
                                           etay_crust_mantle(i,j,k,ispec)*tempdiv2(2)   + &
                                           etaz_crust_mantle(i,j,k,ispec)*tempdiv2(3)   + &
                                           gammax_crust_mantle(i,j,k,ispec)*tempdiv3(1) + &
                                           gammay_crust_mantle(i,j,k,ispec)*tempdiv3(2) + &
                                           gammaz_crust_mantle(i,j,k,ispec)*tempdiv3(3)
          enddo !m
        enddo !i
      enddo !j
    enddo ! k
  enddo !ispec


  ! Finally, calculate the divergence of div_gknl2_cm
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          tempdiv1 = 0._CUSTOM_REAL
          tempdiv2 = 0._CUSTOM_REAL
          tempdiv3 = 0._CUSTOM_REAL

          do l = 1,NGLLX
            tempdiv1 = tempdiv1 + div_gknl2_cm(:,l,j,k,ispec)*hprime_xx(i,l)
            tempdiv2 = tempdiv2 + div_gknl2_cm(:,i,l,k,ispec)*hprime_yy(j,l)
            tempdiv3 = tempdiv3 + div_gknl2_cm(:,i,j,l,ispec)*hprime_zz(k,l)
          enddo !l

          ! This is multiplied by -rho in save_kernels for the final output
          gknl2(i,j,k,ispec) = gknl2(i,j,k,ispec) + (   xix_crust_mantle(i,j,k,ispec)*tempdiv1(1) +&
                                                        xiy_crust_mantle(i,j,k,ispec)*tempdiv1(2) +&
                                                        xiz_crust_mantle(i,j,k,ispec)*tempdiv1(3) +&
                                                       etax_crust_mantle(i,j,k,ispec)*tempdiv2(1) +&
                                                       etay_crust_mantle(i,j,k,ispec)*tempdiv2(2) +&
                                                       etaz_crust_mantle(i,j,k,ispec)*tempdiv2(3) +&
                                                     gammax_crust_mantle(i,j,k,ispec)*tempdiv3(1) +&
                                                     gammay_crust_mantle(i,j,k,ispec)*tempdiv3(2) +&
                                                     gammaz_crust_mantle(i,j,k,ispec)*tempdiv3(3)  )
        enddo !i
      enddo !j
    enddo ! k
  enddo !ispec

  ! free temporary array
  deallocate(gknl2_cm)
  deallocate(div_gknl2_cm)

  end subroutine calculate_second_gravity_kernel

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_compute_crust_mantle_kernels()

! adds full gravity contribution to density kernel (rho_kl) in crust/mantle region
! and computes 1. and 2. gravity-density kernels (rho1siem_kl,rho2siem_kl) contributions

  use constants_solver

  use specfem_par, only: deltat,GPU_MODE, &
    hprime_xx,hprime_yy,hprime_zz

  use specfem_par_crustmantle
  use specfem_par_full_gravity

  implicit none

  integer :: i,j,k,l,m,ispec,iglob,icont
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL) :: tempx1l_phi, tempx2l_phi, tempx3l_phi
  real(kind=CUSTOM_REAL) :: b_tempx1l_phi, b_tempx2l_phi, b_tempx3l_phi

  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl, rhol
  real(kind=CUSTOM_REAL) :: grav_cont, grav0_contraction, sdotphi1, sdotphi2, divS, rho_kl_cm_val

  real(kind=CUSTOM_REAL) :: gradS(3,3), gradphi(3), b_gradphi(3), sdagtmp(3), stmp(3), hmatloc(6)

  logical, save :: do_warn = .true.

  ! safety check
  if (.not. FULL_GRAVITY_VAL) return

  ! not implemented yet on GPU - warning
  if (GPU_MODE) then
    ! kernel contribution not implemented yet - issue a warning
    if (myrank == 0 .and. do_warn) then
      print *,'Warning - Full gravity density kernels not implemented for Crust/Mantle yet on GPUs.'
      ! only print warning once, turn off flag for subsequent calls
      do_warn = .false.
    endif
    ! all done
    return
  endif

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE
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
          !                         behave better for smoother wavefields, thus containing less numerical artefacts.

          ! Original rho kernel (no gravity) contribution value
          rho_kl_cm_val = deltat * (accel_crust_mantle(1,iglob) * b_displ_crust_mantle(1,iglob) &
                                  + accel_crust_mantle(2,iglob) * b_displ_crust_mantle(2,iglob) &
                                  + accel_crust_mantle(3,iglob) * b_displ_crust_mantle(3,iglob) )

          ! already done in compute_kernels_crust_mantle() routine:
          !  rho_kl_crust_mantle(i,j,k,ispec) =  rho_kl_crust_mantle(i,j,k,ispec) + rho_kl_cm_val

          ! Additions due to full gravity

          ! For calculating grad of adjoint displacement (perturb in grav acc) for GLL point:
          tempx1l       = 0._CUSTOM_REAL
          tempy1l       = 0._CUSTOM_REAL
          tempz1l       = 0._CUSTOM_REAL
          tempx2l       = 0._CUSTOM_REAL
          tempy2l       = 0._CUSTOM_REAL
          tempz2l       = 0._CUSTOM_REAL
          tempx3l       = 0._CUSTOM_REAL
          tempy3l       = 0._CUSTOM_REAL
          tempz3l       = 0._CUSTOM_REAL
          ! For calculating \grad\phi (perturb in grav acc) for GLL point (backward and adjoint):
          tempx1l_phi   = 0._CUSTOM_REAL
          tempx2l_phi   = 0._CUSTOM_REAL
          tempx3l_phi   = 0._CUSTOM_REAL
          b_tempx1l_phi = 0._CUSTOM_REAL
          b_tempx2l_phi = 0._CUSTOM_REAL
          b_tempx3l_phi = 0._CUSTOM_REAL

          ! requires NGLLX=NGLLY=NGLLZ otherwise need to be separate for each row
          do l = 1,NGLLX
            ! Adjoint grad phi
            tempx1l_phi = tempx1l_phi + pgrav_cm(ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempx2l_phi = tempx2l_phi + pgrav_cm(ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempx3l_phi = tempx3l_phi + pgrav_cm(ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
            ! Reconstructed forward grad phi
            b_tempx1l_phi = b_tempx1l_phi + b_pgrav_cm(ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            b_tempx2l_phi = b_tempx2l_phi + b_pgrav_cm(ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            b_tempx3l_phi = b_tempx3l_phi + b_pgrav_cm(ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)

            ! Gradient of adjoint displacement field:
            ! s1 * l'a  ;  s1 * l'b   ; s1 * l'c
            tempx1l = tempx1l + displ_crust_mantle(1,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempy1l = tempy1l + displ_crust_mantle(1,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempz1l = tempz1l + displ_crust_mantle(1,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
            ! s2 * l'a  ;  s2 * l'b   ; s2 * l'c
            tempx2l = tempx2l + displ_crust_mantle(2,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempy2l = tempy2l + displ_crust_mantle(2,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempz2l = tempz2l + displ_crust_mantle(2,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
            ! s3 * l'a  ;  s3 * l'b   ; s3 * l'c
            tempx3l = tempx3l + displ_crust_mantle(3,ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempy3l = tempy3l + displ_crust_mantle(3,ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempz3l = tempz3l + displ_crust_mantle(3,ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
          enddo

          ! Get partial derivatives of mapping (del xi/del x) etc
          xixl    = xix_crust_mantle(i,j,k,ispec)
          xiyl    = xiy_crust_mantle(i,j,k,ispec)
          xizl    = xiz_crust_mantle(i,j,k,ispec)
          etaxl   = etax_crust_mantle(i,j,k,ispec)
          etayl   = etay_crust_mantle(i,j,k,ispec)
          etazl   = etaz_crust_mantle(i,j,k,ispec)
          gammaxl = gammax_crust_mantle(i,j,k,ispec)
          gammayl = gammay_crust_mantle(i,j,k,ispec)
          gammazl = gammaz_crust_mantle(i,j,k,ispec)

          ! Adjoint \grad\phi^(dagger)
          gradphi(1) = (xixl * tempx1l_phi) + (etaxl * tempx2l_phi) + (gammaxl * tempx3l_phi)
          gradphi(2) = (xiyl * tempx1l_phi) + (etayl * tempx2l_phi) + (gammayl * tempx3l_phi)
          gradphi(3) = (xizl * tempx1l_phi) + (etazl * tempx2l_phi) + (gammazl * tempx3l_phi)

          ! Reconstructed forward \grad\phi
          b_gradphi(1) = (xixl * b_tempx1l_phi) + (etaxl * b_tempx2l_phi) + (gammaxl * b_tempx3l_phi)
          b_gradphi(2) = (xiyl * b_tempx1l_phi) + (etayl * b_tempx2l_phi) + (gammayl * b_tempx3l_phi)
          b_gradphi(3) = (xizl * b_tempx1l_phi) + (etazl * b_tempx2l_phi) + (gammazl * b_tempx3l_phi)

          ! Gradient of adjoint displacement (\nabla\s^dagger)
          gradS(1,1) = (xixl * tempx1l) + (etaxl * tempy1l) + (gammaxl * tempz1l) ! dS_1 / dx
          gradS(1,2) = (xiyl * tempx1l) + (etayl * tempy1l) + (gammayl * tempz1l) ! dS_1 / dy
          gradS(1,3) = (xizl * tempx1l) + (etazl * tempy1l) + (gammazl * tempz1l) ! dS_1 / dz
          gradS(2,1) = (xixl * tempx2l) + (etaxl * tempy2l) + (gammaxl * tempz2l) ! dS_2 / dx
          gradS(2,2) = (xiyl * tempx2l) + (etayl * tempy2l) + (gammayl * tempz2l) ! dS_2 / dy
          gradS(2,3) = (xizl * tempx2l) + (etazl * tempy2l) + (gammazl * tempz2l) ! dS_2 / dz
          gradS(3,1) = (xixl * tempx3l) + (etaxl * tempy3l) + (gammaxl * tempz3l) ! dS_3 / dx
          gradS(3,2) = (xiyl * tempx3l) + (etayl * tempy3l) + (gammayl * tempz3l) ! dS_3 / dy
          gradS(3,3) = (xizl * tempx3l) + (etazl * tempy3l) + (gammazl * tempz3l) ! dS_3 / dz

          ! Divergence of adjoint (backward) displacement wavefield:
          divS = gradS(1,1) + gradS(2,2) + gradS(3,3)

          ! Calculate gravity terms of density kernel
          grav_cont               = 0.0_CUSTOM_REAL
          grav0_contraction       = 0.0_CUSTOM_REAL
          sdotphi1                = 0.0_CUSTOM_REAL
          sdotphi2                = 0.0_CUSTOM_REAL

          ! Loop to calculate contractions in NDIM = 3
          do icont = 1,3
            ! (1) s_dagger \cdot \nabla \phi
            sdotphi1 = sdotphi1 + (displ_crust_mantle(icont,iglob) * b_gradphi(icont))

            ! (2) s \cdot \nabla \phi_dagger
            sdotphi2 = sdotphi2 + (b_displ_crust_mantle(icont,iglob) * gradphi(icont))

            ! (4) background grav dot adjoint wavefields
            !     \nabla\Phi0 \cdot ( s \cdot \nabla s_dagger  - s \nabla \cdot s_dagger )
            !     icont loop is sum over m
            !     nabla_m Phi ( s_1 nabla_1 s_dagger_m +  s_2 nabla_2 s_dagger_m + s_3 nabla_3 s_dagger_m )
            ! g_cm is the background gravity acc. vector in crustmantle
            grav0_contraction = grav0_contraction &
                                + g_cm(icont,iglob) * (b_displ_crust_mantle(1,iglob)*gradS(icont,1) &
                                                     + b_displ_crust_mantle(2,iglob)*gradS(icont,2) &
                                                     + b_displ_crust_mantle(3,iglob)*gradS(icont,3))

            grav0_contraction = grav0_contraction - (g_cm(icont,iglob) * b_displ_crust_mantle(icont,iglob) * divS)
          enddo ! icont

          ! Local values for the GLL point
          sdagtmp(:) = displ_crust_mantle(:,iglob)
          stmp(:)    = b_displ_crust_mantle(:,iglob)

          hmatloc(:) = gradg_cm(:, iglob)

          ! (3) Gravity contraction s \cdot \nabla\nabla\Phi \cdot \s_dagger
          ! hmatloc is [h11, h22, h33, h12, h13, h23] since symmetric
          grav_cont = grav_cont +  (stmp(1) * hmatloc(1) * sdagtmp(1))
          grav_cont = grav_cont +  (stmp(1) * hmatloc(4) * sdagtmp(2))
          grav_cont = grav_cont +  (stmp(1) * hmatloc(5) * sdagtmp(3))

          grav_cont = grav_cont +  (stmp(2) * hmatloc(4) * sdagtmp(1))
          grav_cont = grav_cont +  (stmp(2) * hmatloc(2) * sdagtmp(2))
          grav_cont = grav_cont +  (stmp(2) * hmatloc(6) * sdagtmp(3))

          grav_cont = grav_cont +  (stmp(3) * hmatloc(5) * sdagtmp(1))
          grav_cont = grav_cont +  (stmp(3) * hmatloc(6) * sdagtmp(2))
          grav_cont = grav_cont +  (stmp(3) * hmatloc(3) * sdagtmp(3))

          !! adds full gravity contribution to density kernel
          ! Combine the 4 gravity terms of the density kernel
          rho_kl_crust_mantle(i,j,k,ispec) = rho_kl_crust_mantle(i,j,k,ispec) &
              + deltat * (sdotphi1 + sdotphi2 + grav_cont + grav0_contraction)

          ! Calculate delta \Phi kernels (require global convolution)
          ! See Eaton paper.
          rhol = rhostore_crust_mantle(i,j,k,ispec) !local density at GLL point

          ! First gravity-density kernel: rho [s \cdot \nabla\s^dag   -  s \nabla\cdot\s^dagger  ]
          ! in index notation:
          ! kernel_l  =   rho * [s_m \nabla_m(\s^dag)_l   -  s_l \nabla_m (\s^dagger)_m ]
          do l = 1,3
            do m = 1,3
              rho1siem_kl_crust_mantle(l,i,j,k,ispec) = rho1siem_kl_crust_mantle(l,i,j,k,ispec) &
                                        + deltat * rhol * (stmp(m)*gradS(l,m) - stmp(l)*gradS(m,m))
            enddo
          enddo


          ! Second gravity-density kernel: \rho \s * \s^dagger
          ! Symmetric portion is taken so 0.5 rho  (s*sdag + sdag*s), where * is tensor product
          ! This is a rank 2 (s, sdagger not contracted)
          do l = 1,3
            do m = 1,3
              rho2siem_kl_crust_mantle(l,m,i,j,k,ispec) = rho2siem_kl_crust_mantle(l,m,i,j,k,ispec) &
                            + 0.5_CUSTOM_REAL * deltat * rhol * (stmp(l)*sdagtmp(m) + stmp(m)*sdagtmp(l))
            enddo
          enddo

          ! Store terms for debugging
          debug_rho_kl_cm(1,i,j,k,ispec) = debug_rho_kl_cm(1,i,j,k,ispec) + rho_kl_cm_val
          debug_rho_kl_cm(2,i,j,k,ispec) = debug_rho_kl_cm(2,i,j,k,ispec) + deltat*sdotphi1
          debug_rho_kl_cm(3,i,j,k,ispec) = debug_rho_kl_cm(3,i,j,k,ispec) + deltat*sdotphi2
          debug_rho_kl_cm(4,i,j,k,ispec) = debug_rho_kl_cm(4,i,j,k,ispec) + deltat*grav_cont
          debug_rho_kl_cm(5,i,j,k,ispec) = debug_rho_kl_cm(5,i,j,k,ispec) + deltat*grav0_contraction

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ
  enddo ! ispec

  end subroutine SIEM_compute_crust_mantle_kernels

!
!-------------------------------------------------------------------------------
!

  subroutine SIEM_save_crust_mantle_kernels()

  use specfem_par
  use specfem_par_crustmantle, only: rhostore_crust_mantle
  use specfem_par_full_gravity

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: scale_kl, scale_kl_int_rho1, scale_kl_int_rho2
  real(kind=CUSTOM_REAL) :: rhol
  real(kind=CUSTOM_REAL), dimension(:,:,:,:),allocatable :: tmpbinoutput

  integer :: ispec,i,j,k
  integer :: ier
  integer :: icomponent, jcomponent

  character(len=80) :: gravfname

  logical, parameter :: SAVE_INTERMEDIATE_GRAV_KERNELS_MESH_BINARIES = .true.

  ! ensight vis
  !logical, parameter :: SAVE_KERNELS_ENSIGHT_GOLD = .true.
  !logical, parameter :: SAVE_INTERMEDIATE_GRAV_KERNELS_MESH_ENSIGHT = .true.
  ! Ensight parameters
  !character(len=80) :: file_head
  !integer :: npart
  !character(len=80) :: destag ! this must be 80 characters long
  !character(len=250)::  grav_file
  !character(len=80),parameter :: ensight_etype='hexa8'

  ! interfaces for visualization routine
  interface
    subroutine write_ensight_pernodeSCALAS(out_fname,destag,npart,n,var)
    implicit none
    character(len=250),intent(in) :: out_fname
    character(len=80),intent(in) :: destag
    integer,intent(in) :: npart,n
    real,dimension(:),intent(in) :: var
    end subroutine write_ensight_pernodeSCALAS
    !============================================
    subroutine write_ensight_perelementAS(out_fname,etype,destag,npart,var)
    implicit none
    character(len=250),intent(in) :: out_fname
    character(len=20),intent(in) :: etype
    character(len=80),intent(in) :: destag
    integer,intent(in) :: npart
    real,dimension(:),intent(in) :: var
    end subroutine write_ensight_perelementAS
    !============================================
  end interface


  ! safety check
  if (.not. FULL_GRAVITY_VAL) return

  ! warning
  if (ADIOS_FOR_KERNELS) then
    if (myrank == 0) then
      print *,'Full gravity kernels: option ADIOS_FOR_KERNELS not implemented yet, saving them as binary files'
    endif
  endif

  ! scaling factors
  ! kernel unit [ s / km^3 ]
  scale_kl = real(scale_t * scale_displ_inv * 1.d9,kind=CUSTOM_REAL)

  ! intermediate scaling factors for gravity kernels:
  scale_kl_int_rho1 = real(scale_kl/(GRAV*scale_displ),kind=CUSTOM_REAL) ! vector rho kernel
  scale_kl_int_rho2 = real(scale_kl/GRAV,kind=CUSTOM_REAL)               ! matrix rho kernel

  ! crust_mantle
  do ispec = 1, NSPEC_CRUST_MANTLE
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          rhol = rhostore_crust_mantle(i,j,k,ispec)

          ! Gravity kernels:
          gknl1(i,j,k,ispec) = - rhol * gknl1(i,j,k,ispec) * scale_kl
          gknl2(i,j,k,ispec) = - rhol * gknl2(i,j,k,ispec) * scale_kl

          debug_rho_kl_cm(:,i,j,k,ispec) =  - rhol * debug_rho_kl_cm(:,i,j,k,ispec) * scale_kl
        enddo
      enddo
    enddo
  enddo

  ! 1. and 2. gravity-kernels
  ! Scale
  rho1siem_kl_crust_mantle = rho1siem_kl_crust_mantle * scale_kl_int_rho1
  rho2siem_kl_crust_mantle = rho2siem_kl_crust_mantle * scale_kl_int_rho2

  ! store files
  open(unit=IOUT,file=trim(prname)//'grav1_kernel.bin',status='unknown',form='unformatted',action='write')
  write(IOUT) gknl1
  close(IOUT)
  open(unit=IOUT,file=trim(prname)//'grav2_kernel.bin',status='unknown',form='unformatted',action='write')
  write(IOUT) gknl2
  close(IOUT)

  ! intermediate kernels
  if (SAVE_INTERMEDIATE_GRAV_KERNELS_MESH_BINARIES) then
    ! unsmoothed kernel (ie before SIEM integration)

    allocate(tmpbinoutput(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT),stat=ier)
    if (ier /= 0) stop 'Error allocating tmpbinoutput array'
    tmpbinoutput(:,:,:,:) = 0.0_CUSTOM_REAL

    ! first intermediate kernel
    do icomponent = 1,3
      write(gravfname,'(a,i1,a)')trim(prname)//'unconv', icomponent, '_grav1_kernel.bin'

      !debug
      if (myrank == 0) print *,'SIEM save kernels: Filename: ', gravfname

      tmpbinoutput = rho1siem_kl_crust_mantle(icomponent,:,:,:,:)

      open(unit=IOUT,file=gravfname,status='unknown',form='unformatted',action='write')
      write(IOUT) tmpbinoutput
      close(IOUT)
    enddo

    ! second intermediate kernel
    do icomponent = 1,3
      do jcomponent = 1,3
        write(gravfname,'(a,i1,i1,a)')trim(prname)//'unconv', icomponent,jcomponent, '_grav2_kernel.bin'

        !debug
        if (myrank == 0) print *,'SIEM save kernels: Filename: ', gravfname

        tmpbinoutput = rho2siem_kl_crust_mantle(icomponent,jcomponent,:,:,:,:)

        open(unit=IOUT,file=gravfname,status='unknown',form='unformatted',action='write')
        write(IOUT) tmpbinoutput
        close(IOUT)
      enddo
    enddo

    ! debugging files:
    do icomponent = 1,5
      write(gravfname,'(a,i1,a)')trim(prname)//'rho_debug', icomponent, '.bin'

      !debug
      if (myrank == 0) print *,'SIEM save kernels: Filename: ', gravfname

      tmpbinoutput =  debug_rho_kl_cm(icomponent,:,:,:,:)

      open(unit=IOUT,file=gravfname,status='unknown',form='unformatted',action='write')
      write(IOUT) tmpbinoutput
      close(IOUT)
    enddo

    ! free temporary array
    deallocate(tmpbinoutput)
  endif

  ! ensight vis
!  ! Saving as Ensight files
!  if (SAVE_KERNELS_ENSIGHT_GOLD) then
!    destag = 'unstructured meshes'
!    npart  = 1
!    allocate(ensight_tmp(NGLOB_CRUST_MANTLE))
!
!    ! Density (rhonotprime) kernel: rhonotprime_kl_crust_mantle
!    if (myrank==0) write(*,*) ' Writing rho kernel to Ensight...'
!    ! First need to reformat to global:
!    ensight_tmp(:) = 0.0_CUSTOM_REAL
!    do ispec=1,NSPEC_CRUST_MANTLE
!      do k = 1,NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = rhonotprime_kl_crust_mantle(i,j,k,ispec)
!          enddo !i
!        enddo !j
!      enddo !k
!    enddo !i
!    ! Save to ensight
!    write(file_head,'(a,i1,a)')'reg',1,'_proc'
!    grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.rhokernel'
!    call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!
!    ! Alpha kernel: alpha_kl_crust_mantle
!    if (myrank==0) write(*,*) ' Writing alpha kernel to Ensight...'
!    ! First need to reformat to global:
!    ensight_tmp(:) = 0.0_CUSTOM_REAL
!    do ispec=1,NSPEC_CRUST_MANTLE
!      do k = 1,NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = alpha_kl_crust_mantle(i,j,k,ispec)
!          enddo !i
!        enddo !j
!      enddo !k
!    enddo !i
!    ! Save to ensight
!    write(file_head,'(a,i1,a)')'reg',1,'_proc'
!    grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.alphakernel'
!    call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!
!    ! Beta kernel: beta_kl_crust_mantle
!    if (myrank==0) write(*,*) ' Writing beta kernel to Ensight...'
!    ! First need to reformat to global:
!    ensight_tmp(:) = 0.0_CUSTOM_REAL
!    do ispec=1,NSPEC_CRUST_MANTLE
!      do k = 1,NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = beta_kl_crust_mantle(i,j,k,ispec)
!          enddo !i
!        enddo !j
!      enddo !k
!    enddo !i
!    ! Save to ensight
!    write(file_head,'(a,i1,a)')'reg',1,'_proc'
!    grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.betakernel'
!    call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!
!    ! Ensight gravity kernels:
!    ! First gravity kernel:
!    if (myrank==0) write(*,*) ' Writing 1st gravity kernel to Ensight...'
!    ! First need to reformat to global:
!    ensight_tmp(:) = 0.0_CUSTOM_REAL
!    do ispec=1,NSPEC_CRUST_MANTLE
!      do k = 1,NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = gknl1(i,j,k,ispec)
!          enddo !i
!        enddo !j
!      enddo !k
!    enddo !i
!    ! Save to ensight
!    write(file_head,'(a,i1,a)')'reg',1,'_proc'
!    grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav1kernel'
!    call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!
!    ! Second gravity kernel:
!    if (myrank==0) write(*,*) ' Writing 2nd gravity kernel to Ensight...'
!    ! First need to reformat to global:
!    ensight_tmp(:) = 0.0_CUSTOM_REAL
!    do ispec=1,NSPEC_CRUST_MANTLE
!      do k = 1,NGLLZ
!        do j = 1, NGLLY
!          do i = 1, NGLLX
!            ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = gknl2(i,j,k,ispec)
!          enddo !i
!        enddo !j
!      enddo !k
!    enddo !i
!    ! Save to ensight
!    write(file_head,'(a,i1,a)')'reg',1,'_proc'
!    grav_file=trim(LOCAL_PATH)//'/'//trim(file_head)//trim(ptail)//'.grav2kernel'
!    call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!
!    if (SAVE_INTERMEDIATE_GRAV_KERNELS_MESH_ENSIGHT) then
!      if (myrank==0) write(*,*) ' Writing 1st intermediate gravity kernel to Ensight...'
!      ! loop each component to write out
!      do icomponent=1,3
!        ensight_tmp(:) = 0.0_CUSTOM_REAL
!        do ispec=1,NSPEC_CRUST_MANTLE
!          do k = 1,NGLLZ
!            do j = 1, NGLLY
!              do i = 1, NGLLX
!                ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = rho1siem_kl_crust_mantle(icomponent,i,j,k,ispec)
!              enddo !i
!            enddo !j
!          enddo !k
!        enddo !i
!        write(file_head,'(a,i1,a)')'reg',1,'_proc'
!        write(grav_file, '(a,a,a,a,a,i1)')trim(LOCAL_PATH), '/', trim(file_head), trim(ptail), '.rho1int', icomponent
!        call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!      enddo
!
!      if (myrank==0) write(*,*) ' Writing 2nd intermediate gravity kernel to Ensight...'
!      ! loop each component to write out
!      do jcomponent=1,3
!        do icomponent=1,3
!          ensight_tmp(:) = 0.0_CUSTOM_REAL
!          do ispec=1,NSPEC_CRUST_MANTLE
!            do k = 1,NGLLZ
!              do j = 1, NGLLY
!                do i = 1, NGLLX
!                  ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = rho2siem_kl_crust_mantle(icomponent,jcomponent,i,j,k,ispec)
!                enddo !i
!              enddo !j
!            enddo !k
!          enddo !i
!          write(file_head,'(a,i1,a)')'reg',1,'_proc'
!          write(grav_file, '(a,a,a,a,a,i1,i1)') trim(LOCAL_PATH), '/', trim(file_head), trim(ptail), &
!                                                '.rho2int', icomponent, jcomponent
!          call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!        enddo ! icomponent
!      enddo ! jcomponent
!
!      ! Debug files from first gravity kernel:
!      do icomponent=1,5
!        ensight_tmp(:) = 0.0_CUSTOM_REAL
!        do ispec=1,NSPEC_CRUST_MANTLE
!          do k = 1,NGLLZ
!            do j = 1, NGLLY
!              do i = 1, NGLLX
!                ensight_tmp(ibool_crust_mantle(i,j,k,ispec)) = debug_rho_kl_cm(icomponent,i,j,k,ispec)
!              enddo !i
!            enddo !j
!          enddo !k
!        enddo !i
!        write(file_head,'(a,i1,a)')'reg',1,'_proc'
!        write(grav_file, '(a,a,a,a,a,i1,i1)') trim(LOCAL_PATH), '/', trim(file_head), trim(ptail), &
!                                              '.debug_rho', icomponent
!        call write_ensight_pernodeSCALAS(grav_file,destag,npart,nnode4_cm,real(ensight_tmp(inode4_cm)))
!      enddo
!    endif ! SAVE_INTERMEDIATE_GRAV_KERNELS_MESH_ENSIGHT
!  endif ! SAVE_KERNELS_ENSIGHT_GOLD

  end subroutine SIEM_save_crust_mantle_kernels
