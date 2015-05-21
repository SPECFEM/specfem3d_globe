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

! add_model_globe_tiso_cg
!
! this program can be used to update TRANSVERSE ISOTROPIC model files
! based on smoothed event kernels.
! the kernels are given for tranverse isotropic parameters (bulk_c,bulk_betav,bulk_betah,eta).
!
! the algorithm uses a conjugate gradient method with a step length
! limited by the given maximum update percentage.
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
!       proc000***_reg1_bulk_betav_kernel_smooth.bin &
!       proc000***_reg1_bulk_betah_kernel_smooth.bin &
!       proc000***_reg1_eta_kernel_smooth.bin &
!       proc000***_reg1_rho_kernel_smooth.bin
!
!- /tigress-hsm/dpeter/SPECFEM3D_GLOBE/KERNELS/OUTPUT_SUM.old/ contains old gradients:
!       proc000***_reg1_bulk_c_kernel_smooth.bin &
!       proc000***_reg1_bulk_betav_kernel_smooth.bin &
!       proc000***_reg1_bulk_betah_kernel_smooth.bin &
!       proc000***_reg1_eta_kernel_smooth.bin &
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
! usage: ./add_model_globe_tiso_cg 0.3
!
!
! NOTE: this routine uses previous model update in OUTPUT_SUM.old/
!             for a conjugate gradient update
!


program add_model

  use tomography_model_tiso
  use tomography_kernels_tiso
  use tomography_kernels_tiso_cg

  implicit none
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: betav1,betah1,betav0,betah0,rho1,rho0, &
    betaiso1,betaiso0,eta1,eta0,alphav1,alphav0,alphah1,alphah0
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
    print *,'program add_model_tiso_cg: '
    print *,'  NPROC_XI , NPROC_ETA: ',nproc_xi_val,nproc_eta_val
    print *,'  NCHUNKS: ',nchunks_val
    print *
    print *,'model update for vsv,vsh,vpv,vph,eta,rho:'
    print *,'  step_fac = ',step_fac
    print *
    print *,'***********'
    print *
  endif

  ! reads in current transverse isotropic model files: vpv.. & vsv.. & eta & rho
  call read_model_tiso()

  ! reads in smoothed kernels: bulk, betav, betah, eta
  call read_kernels_tiso()

  ! reads in old (former inversion) smoothed kernels: bulk, betav, betah, eta
  call read_kernels_cg_tiso_old()

  ! calculates gradient
  ! conjugate gradient method
  call get_gradient_cg_tiso()

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
  !   transverse isotropic update only in layer Moho to 220 (where SPECFEM3D_GLOBE considers TISO)
  !   everywhere else uses an isotropic update
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

          ! do not use transverse isotropy except if element is between d220 and Moho
!          if (.not. ( idoubling(ispec)==IFLAG_220_80 .or. idoubling(ispec)==IFLAG_80_MOHO)) then
          if (.not. ispec_is_tiso(ispec)) then
            ! isotropic model update

            ! no eta perturbation, since eta = 1 in isotropic media
            eta1 = eta0

            ! shear values
            ! isotropic kernel K_beta = K_betav + K_betah
            ! with same scaling step_length the model update dbeta_iso = dbetav + dbetah
            ! note:
            !   this step length can be twice as big as that given by the input
            dbetaiso = model_dbetav(i,j,k,ispec) + model_dbetah(i,j,k,ispec)
            betav1 = betav0 * exp( dbetaiso )
            betah1 = betah0 * exp( dbetaiso )
            ! note: betah is probably not really used in isotropic layers
            !         (see SPECFEM3D_GLOBE/get_model.f90)

            ! density: uses scaling relation with isotropic shear perturbations
            !               dln rho = RHO_SCALING * dln betaiso
            rho1 = rho0 * exp( RHO_SCALING * dbetaiso )

            ! alpha values
            dbulk = model_dbulk(i,j,k,ispec)
            alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betav0**2 * ( &
                                exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
            alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) + FOUR_THIRDS * betah0**2 * ( &
                                exp(2.0*dbetaiso) - exp(2.0*dbulk) ) )
            ! note: alphah probably not used in SPECFEM3D_GLOBE

          else

            ! transverse isotropic model update

            ! eta value : limits updated values for eta range constraint
            eta1 = eta0 * exp( model_deta(i,j,k,ispec) )
            if (eta1 < LIMIT_ETA_MIN ) eta1 = LIMIT_ETA_MIN
            if (eta1 > LIMIT_ETA_MAX ) eta1 = LIMIT_ETA_MAX

            ! shear values
            betav1 = betav0 * exp( model_dbetav(i,j,k,ispec) )
            betah1 = betah0 * exp( model_dbetah(i,j,k,ispec) )

            ! density: uses scaling relation with Voigt average of shear perturbations
            betaiso0 = sqrt(  ( 2.0 * betav0**2 + betah0**2 ) / 3.0 )
            betaiso1 = sqrt(  ( 2.0 * betav1**2 + betah1**2 ) / 3.0 )
            dbetaiso = log( betaiso1 / betaiso0 )
            rho1 = rho0 * exp( RHO_SCALING * dbetaiso )

            ! alpha values
            dbulk = model_dbulk(i,j,k,ispec)
            alphav1 = sqrt( alphav0**2 * exp(2.0*dbulk) &
                            + FOUR_THIRDS * betav0**2 * ( &
                                exp(2.0*model_dbetav(i,j,k,ispec)) - exp(2.0*dbulk) ) )
            alphah1 = sqrt( alphah0**2 * exp(2.0*dbulk) &
                            + FOUR_THIRDS * betah0**2 * ( &
                                exp(2.0*model_dbetah(i,j,k,ispec)) - exp(2.0*dbulk) ) )

          endif


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
  call compute_kernel_integral_tiso()

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

  ! needs density scaling
  if (.not. USE_RHO_SCALING) stop 'Error this model update needs density scaling to be turned on...'

end subroutine initialize

