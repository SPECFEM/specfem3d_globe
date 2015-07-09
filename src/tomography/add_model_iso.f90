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

! add_model_globe_iso
!
! this program can be used to update ISOTROPIC model files with
! (smoothed & summed) event kernels.
! the kernels are given for isotropic parameters (alpha,beta,rho) or ( bulk_c,beta,rho).
!
! the algorithm uses a steepest descent method with a step length
! determined by the given maximum update percentage.
!
! input:
!    - step_fac : step length to update the models, f.e. 0.03 for plusminus 3%
!
! setup:
!
!- INPUT_MODEL/  contains:
!       proc000***_reg1_vs.bin &
!       proc000***_reg1_vp.bin &
!       proc000***_reg1_rho.bin
!
!- INPUT_GRADIENT/ contains:
!       proc000***_reg1_bulk_c_kernel_smooth.bin &
!       proc000***_reg1_bulk_beta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!     or
!       proc000***_reg1_alpha_kernel_smooth.bin &
!       proc000***_reg1_beta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!
!- topo/ contains:
!       proc000***_reg1_solver_data.bin
!
! new models are stored in
!- OUTPUT_MODEL/ as
!   proc000***_reg1_vp_new.bin and
!   proc000***_reg1_vs_new.bin and
!   proc000***_reg1_rho_new.bin and
!
! USAGE: e.g. ./add_model_globe_iso 0.3


program add_model

  use tomography_model_iso
  use tomography_kernels_iso

  implicit none
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: beta1,beta0,rho1,rho0,alpha1,alpha0
  real(kind=CUSTOM_REAL) :: dbetaiso,dbulk

  ! ============ program starts here =====================

  ! initializes arrays
  call initialize()

  ! reads in parameters needed
  call read_parameters_tomo()

  ! user output
  if (myrank == 0) then
    print*
    print*,'***********'
    print*,'program add_model_iso: '
    print*,'  NPROC_XI , NPROC_ETA: ',nproc_xi_val,nproc_eta_val
    print*,'  NCHUNKS: ',nchunks_val
    print*
    print*,'model update for vs & vp & rho'
    print*,'  step_fac = ',step_fac
    print*
    if (USE_ALPHA_BETA_RHO) then
      print*,'kernel parameterization: (alpha,beta,rho)'
    else
      print*,'kernel parameterization: (bulk,beta,rho)'
    endif
    print*
    if (USE_RHO_SCALING) then
      print*,'scaling rho perturbations'
      print*
    endif
    print*,'***********'
    print*
  endif

  ! reads in current isotropic model files: vp & vs & rho
  call read_model_iso()

  ! reads in smoothed kernels: bulk, beta, rho
  call read_kernels_iso()

  ! calculates gradient
  ! steepest descent method
  call get_gradient_steepest_iso()

  ! computes new model values for alpha, beta and rho
  ! and stores new model files
  ! allocate new model arrays
  allocate(model_vp_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vs_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vp_new = 0.0_CUSTOM_REAL
  model_vs_new = 0.0_CUSTOM_REAL
  model_rho_new = 0.0_CUSTOM_REAL

  ! model update:
  !   isotropic update everywhere
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initial model values
          beta0 = model_vs(i,j,k,ispec)
          rho0 = model_rho(i,j,k,ispec)
          alpha0 = model_vp(i,j,k,ispec)

          beta1 = 0._CUSTOM_REAL
          rho1 = 0._CUSTOM_REAL
          alpha1 = 0._CUSTOM_REAL

          ! isotropic model update

          ! shear values
          dbetaiso = model_dbeta(i,j,k,ispec)
          beta1 = beta0 * exp( dbetaiso )

          ! density
          rho1 = rho0 * exp( model_drho(i,j,k,ispec) )

          ! alpha values
          dbulk = model_dbulk(i,j,k,ispec)
          if (USE_ALPHA_BETA_RHO) then
            ! new vp values use alpha model update
            alpha1 = alpha0 * exp( dbulk )
          else
            ! new vp values use bulk model update:
            ! this is based on vp_new = sqrt( bulk_new**2 + 4/3 vs_new**2 )
            alpha1 = sqrt( alpha0**2 * exp(2.0*dbulk) + FOUR_THIRDS * beta0**2 * ( &
                              exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
          endif

          ! stores new model values
          model_vp_new(i,j,k,ispec) = alpha1
          model_vs_new(i,j,k,ispec) = beta1
          model_rho_new(i,j,k,ispec) = rho1

        enddo
      enddo
    enddo
  enddo
  call synchronize_all()

  ! stores new model in files
  call write_new_model_iso()

  ! stores relative model perturbations
  call write_new_model_perturbations_iso()

  ! computes volume element associated with points, calculates kernel integral for statistics
  call compute_kernel_integral_iso()

  ! stop all the MPI processes, and exit
  call finalize_mpi()

end program add_model


!
!-------------------------------------------------------------------------------------------------
!

subroutine initialize()

! initializes arrays

  use tomography_par

  implicit none

  ! initialize the MPI communicator and start the NPROCTOT MPI processes
  call init_mpi()
  call world_size(sizeprocs)
  call world_rank(myrank)

  if (sizeprocs /= NPROCTOT_VAL) then
    if (myrank == 0) then
      print*, 'Error number of processors supposed to run on : ',NPROCTOT_VAL
      print*, 'Error number of MPI processors actually run on: ',sizeprocs
      print*
      print*, 'please rerun with: mpirun -np ',NPROCTOT_VAL,' bin/xadd_model .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! sets tomography array dimensions
  NSPEC = NSPEC_CRUST_MANTLE
  NGLOB = NGLOB_CRUST_MANTLE

end subroutine initialize

