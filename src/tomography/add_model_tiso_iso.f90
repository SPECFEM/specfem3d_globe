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

! add_model_globe_tiso_iso
!
! this program can be used to update TRANSVERSE ISOTROPIC model files
! based on smoothed, summed, ISOTROPIC event kernels.
! the kernels are given for isotropic parameters (bulk_c,bulk_beta,rho) or (alpha,beta,rho).
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
!       proc000***_reg1_vsv.bin &
!       proc000***_reg1_vsh.bin &
!       proc000***_reg1_vpv.bin &
!       proc000***_reg1_vph.bin &
!       proc000***_reg1_eta.bin &
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
!       proc000***_reg1_solver_data_1.bin
!
! new models are stored in
!- OUTPUT_MODEL/ as
!   proc000***_reg1_vpv_new.bin &
!   proc000***_reg1_vph_new.bin &
!   proc000***_reg1_vsv_new.bin &
!   proc000***_reg1_vsh_new.bin &
!   proc000***_reg1_eta_new.bin &
!   proc000***_reg1_rho_new.bin
!
! USAGE: ./add_model_globe_tiso_iso 0.3


program add_model

  use tomography_model_tiso
  use tomography_kernels_iso

  implicit none
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: betav1,betah1,betav0,betah0,rho1,rho0, &
    eta1,eta0,alphav1,alphav0,alphah1,alphah0
  real(kind=CUSTOM_REAL) :: dbetaiso,dbulk

  ! ============ program starts here =====================

  ! initializes arrays
  call initialize()

  ! reads in parameters needed
  call read_parameters_tomo()

  ! user output
  if (myrank == 0) then
    print *
    print *,'***********'
    print *,'program add_model_tiso_iso: '
    print *,'  NPROC_XI , NPROC_ETA: ',nproc_xi_val,nproc_eta_val
    print *,'  NCHUNKS: ',nchunks_val
    print *
    print *,'model update for vsv,vsh,vpv,vph,eta,rho:'
    print *,'  step_fac = ',step_fac
    print *
    if (USE_ALPHA_BETA_RHO) then
      print *,'kernel parameterization: (alpha,beta,rho)'
    else
      print *,'kernel parameterization: (bulk,beta,rho)'
    endif
    print *
    if (USE_RHO_SCALING) then
      print *,'scaling rho perturbations'
      print *
    endif
    print *,'***********'
    print *
  endif

  ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
  call read_model_tiso()

  ! reads in smoothed kernels: bulk, beta, rho
  call read_kernels_iso()

  ! calculates gradient
  ! steepest descent method
  call get_gradient_steepest_iso()

  ! compute new model in terms of alpha, beta, eta and rho
  ! (see also Carl's Latex notes)
  ! allocate new model arrays
  allocate(model_vpv_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vph_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vsv_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_vsh_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_eta_new(NGLLX,NGLLY,NGLLZ,NSPEC), &
           model_rho_new(NGLLX,NGLLY,NGLLZ,NSPEC),stat=ier)
  if (ier /= 0) stop 'Error allocating model arrays'

  ! initializes arrays
  model_vpv_new = 0.0_CUSTOM_REAL
  model_vph_new = 0.0_CUSTOM_REAL
  model_vsv_new = 0.0_CUSTOM_REAL
  model_vsh_new = 0.0_CUSTOM_REAL
  model_eta_new = 0.0_CUSTOM_REAL
  model_rho_new = 0.0_CUSTOM_REAL

  ! model update:
  !   isotropic update everywhere
  do ispec = 1, NSPEC
    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          ! initial model values
          eta0 = model_eta(i,j,k,ispec)
          betav0 = model_vsv(i,j,k,ispec)
          betah0 = model_vsh(i,j,k,ispec)
          rho0 = model_rho(i,j,k,ispec)
          alphav0 = model_vpv(i,j,k,ispec)
          alphah0 = model_vph(i,j,k,ispec)

          eta1 = 0._CUSTOM_REAL
          betav1 = 0._CUSTOM_REAL
          betah1 = 0._CUSTOM_REAL
          rho1 = 0._CUSTOM_REAL
          alphav1 = 0._CUSTOM_REAL
          alphah1 = 0._CUSTOM_REAL

          ! isotropic model update

          ! no eta perturbation for isotropic updates, takes the same as initial values
          eta1 = eta0

          ! shear values
          dbetaiso = model_dbeta(i,j,k,ispec)
          betav1 = betav0 * exp( dbetaiso )
          betah1 = betah0 * exp( dbetaiso )
          ! note: betah is probably not really used in isotropic layers
          !         (see SPECFEM3D_GLOBE/get_model.f90)

          ! density
          rho1 = rho0 * exp( model_drho(i,j,k,ispec) )

          ! alpha values
          dbulk = model_dbulk(i,j,k,ispec)
          if (USE_ALPHA_BETA_RHO) then
            ! new vp values use alpha model update
            alphav1 = alphav0 * exp( dbulk )
            alphah1 = alphah0 * exp( dbulk )
          else
            ! new vp values use bulk model update:
            ! this is based on vp_new = sqrt( bulk_new**2 + 4/3 vs_new**2 )
            alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betav0**2 * ( &
                              exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
            alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betah0**2 * ( &
                              exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
          endif
          ! note: alphah probably not used in SPECFEM3D_GLOBE

          ! stores new model values
          model_vpv_new(i,j,k,ispec) = alphav1
          model_vph_new(i,j,k,ispec) = alphah1
          model_vsv_new(i,j,k,ispec) = betav1
          model_vsh_new(i,j,k,ispec) = betah1
          model_eta_new(i,j,k,ispec) = eta1
          model_rho_new(i,j,k,ispec) = rho1

        enddo
      enddo
    enddo
  enddo

  ! stores new model in files
  call write_new_model_tiso()

  ! stores relative model perturbations
  call write_new_model_perturbations_tiso()

  ! computes volume element associated with points, calculates kernel integral for statistics
  call compute_kernel_integral_tiso_iso()

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
      print *, 'Error number of processors supposed to run on : ',NPROCTOT_VAL
      print *, 'Error number of MPI processors actually run on: ',sizeprocs
      print *
      print *, 'please rerun with: mpirun -np ',NPROCTOT_VAL,' bin/xadd_model .. '
    endif
    call exit_MPI(myrank,'Error wrong number of MPI processes')
  endif

  ! sets tomography array dimensions
  NSPEC = NSPEC_CRUST_MANTLE
  NGLOB = NGLOB_CRUST_MANTLE

end subroutine initialize

