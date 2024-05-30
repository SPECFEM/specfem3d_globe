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
    gravload1, pgrav1, pgrav_cm1, pgrav_cm, &
    is_active_gll,igll_active_on

  use siem_math_library_mpi, only: maxvec
  use siem_poisson, only: compute_grav_kl1_load
  use siem_solver_petsc, only: petsc_set_vector1, petsc_zero_initialguess1, petsc_solve1
  use siem_solver_mpi, only: interpolate3to5

  implicit none

  ! Local variables
  integer :: ispec,i,j,k,l,ier
  integer :: icomponent
  real(kind=CUSTOM_REAL) :: tempx1l, tempx2l, tempx3l, tempy1l, tempy2l, tempy3l, tempz1l, tempz2l, tempz3l
  real(kind=CUSTOM_REAL) :: maxload
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: gknl1_cm ! (3,NGLOB_CRUST_MANTLE)

  ! safety check
  if (POISSON_SOLVER == ISOLVER_BUILTIN) then
    !TODO: full gravity builtin solver for kernels
    print *,'ERROR: builtin solver not setup for gravity kernels'
    call exit_MPI(myrank,'Error builtin solver not setup for gravity kernels')
  endif

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

    if (POISSON_SOLVER == ISOLVER_PETSC) then
      ! Petsc solver
      call petsc_set_vector1(gravload1)

      ! Need to zero guess since was previously for poisson eqn
      ! Do we need this - is zero guess set to PETSC_TRUE?
      call petsc_zero_initialguess1()

      pgrav1(:) = 0.0_CUSTOM_REAL
      ! Here we use pgrav as the vector the solution is put into, just to
      ! save on allocating another array. Note that we could allocate a
      ! separate temporary array but pgrav isnt used after this (apart from
      ! if the forward wavefield is being saved, but it shouldnt be for SIMTYPE 3)
      call petsc_solve1(pgrav1(1:))
    endif

    ! Now interpolate this component to the 5GLL setup (Crust mantle only):
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
    gravload1, pgrav1, pgrav_cm1, pgrav_cm, &
    is_active_gll,igll_active_on

  use siem_math_library_mpi, only: maxvec
  use siem_poisson, only: compute_grav_kl2_load
  use siem_solver_petsc, only: petsc_set_vector1, petsc_zero_initialguess1, petsc_solve1
  use siem_solver_mpi, only: interpolate3to5

  implicit none

  ! Local variables
  integer :: ispec,i,j,k,l,m,ier
  integer :: icomponent, jcomponent
  real(kind=CUSTOM_REAL) :: tempdiv1(3),tempdiv2(3),tempdiv3(3)
  real(kind=CUSTOM_REAL) :: maxload
  real(kind=CUSTOM_REAL),dimension(:,:,:), allocatable :: gknl2_cm
  real(kind=CUSTOM_REAL),dimension(:,:,:,:,:), allocatable :: div_gknl2_cm

  ! safety check
  if (POISSON_SOLVER == ISOLVER_BUILTIN) then
    !TODO: full gravity builtin solver for kernels
    print *,'ERROR: builtin solver not setup for gravity kernels'
    call exit_MPI(myrank,'Error builtin solver not setup for gravity kernels')
  endif

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

      if (POISSON_SOLVER == ISOLVER_PETSC) then
        ! PETSc solver
        call petsc_set_vector1(gravload1)

        call petsc_zero_initialguess1()

        pgrav1(:) = 0.0_CUSTOM_REAL
        call petsc_solve1(pgrav1(1:))
      endif

      pgrav_cm(:)  = zero
      pgrav_cm1(:) = zero
      pgrav_cm1(:) = pgrav1(gdof_cm1(:))

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

