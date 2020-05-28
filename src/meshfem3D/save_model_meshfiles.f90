

  subroutine save_model_meshfiles()

! outputs model files in binary format

  use constants
  use shared_parameters, only: R_PLANET,RHOAV

  use meshfem3D_par, only: nspec,iregion_code

  use meshfem3D_models_par, only: &
    TRANSVERSE_ISOTROPY,ATTENUATION,ATTENUATION_3D,ATTENUATION_1D_WITH_3D_STORAGE, &
    HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE

  use regions_mesh_par2, only: &
    rhostore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
    Qmu_store,Gc_prime_store,Gs_prime_store,mu0store, &
    prname

  implicit none

  ! local parameters
  integer :: i,j,k,ispec,ier
  real(kind=CUSTOM_REAL) :: scaleval1,scaleval2,scaleval,scale_GPa
  real(kind=CUSTOM_REAL),dimension(:,:,:,:),allocatable :: temp_store

  ! scaling factors to re-dimensionalize units
  scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_PLANET/1000.0d0) )
  scaleval2 = sngl( RHOAV/1000.0d0 )

  ! uses temporary array
  allocate(temp_store(NGLLX,NGLLY,NGLLZ,nspec),stat=ier)
  if (ier /= 0) stop 'Error allocating temp_store array'
  temp_store(:,:,:,:) = 0._CUSTOM_REAL

  ! isotropic model
  ! vp
  open(unit=IOUT,file=prname(1:len_trim(prname))//'vp.bin', &
       status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening vp.bin file')

  temp_store(:,:,:,:) = sqrt((kappavstore(:,:,:,:) + 4.0_CUSTOM_REAL * muvstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                        * scaleval1
  write(IOUT) temp_store
  close(IOUT)
  ! vs
  open(unit=IOUT,file=prname(1:len_trim(prname))//'vs.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening vs.bin file')

  temp_store(:,:,:,:) = sqrt( muvstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
  write(IOUT) temp_store
  close(IOUT)
  ! rho
  open(unit=IOUT,file=prname(1:len_trim(prname))//'rho.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error opening rho.bin file')

  temp_store(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
  write(IOUT) temp_store
  close(IOUT)

  ! transverse isotropic model
  if (TRANSVERSE_ISOTROPY) then
    ! vpv
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vpv.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vpv.bin file')

    temp_store(:,:,:,:) = sqrt((kappavstore(:,:,:,:) + 4.0_CUSTOM_REAL * muvstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! vph
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vph.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vph.bin file')

    temp_store(:,:,:,:) = sqrt((kappahstore(:,:,:,:) + 4.0_CUSTOM_REAL * muhstore(:,:,:,:)/3.0_CUSTOM_REAL)/rhostore(:,:,:,:)) &
                          * scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! vsv
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vsv.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vsv.bin file')

    temp_store(:,:,:,:) = sqrt( muvstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! vsh
    open(unit=IOUT,file=prname(1:len_trim(prname))//'vsh.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening vsh.bin file')

    temp_store(:,:,:,:) = sqrt( muhstore(:,:,:,:)/rhostore(:,:,:,:) )*scaleval1
    write(IOUT) temp_store
    close(IOUT)

    ! rho
    open(unit=IOUT,file=prname(1:len_trim(prname))//'rho.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening rho.bin file')

    temp_store(:,:,:,:) = rhostore(:,:,:,:) * scaleval2
    write(IOUT) temp_store
    close(IOUT)

    ! eta
    open(unit=IOUT,file=prname(1:len_trim(prname))//'eta.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening eta.bin file')
    write(IOUT) eta_anisostore
    close(IOUT)
  endif ! TRANSVERSE_ISOTROPY

  ! anisotropic values
  if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! the scale of GPa--[g/cm^3][(km/s)^2]
    scaleval = dsqrt(PI*GRAV*RHOAV)
    scale_GPa = (RHOAV/1000.d0)*((R_PLANET*scaleval/1000.d0)**2)

    ! Gc_prime
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Gc_prime.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening Gc_prime.bin file')
    write(IOUT) Gc_prime_store
    close(IOUT)

    ! Gs_prime
    open(unit=IOUT,file=prname(1:len_trim(prname))//'Gs_prime.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening Gs_prime.bin file')
    write(IOUT) Gs_prime_store
    close(IOUT)

    ! shear modulus mu0
    open(unit=IOUT,file=prname(1:len_trim(prname))//'mu0.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening mu0.bin file')
    temp_store(:,:,:,:) = mu0store(:,:,:,:) * scale_GPa
    write(IOUT) temp_store
    close(IOUT)
  endif

  ! shear attenuation
  if (ATTENUATION) then
    ! saves Qmu_store to full CUSTOM_REAL array
    if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then
      ! attenuation arrays are fully 3D
      temp_store(:,:,:,:) = Qmu_store(:,:,:,:)
    else
      ! attenuation array dimensions: Q_mustore(1,1,1,nspec)
      do ispec = 1,nspec
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              temp_store(i,j,k,ispec) = Qmu_store(1,1,1,ispec)
            enddo
          enddo
        enddo
      enddo
    endif

    ! Qmu
    open(unit=IOUT,file=prname(1:len_trim(prname))//'qmu.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error opening qmu.bin file')
    write(IOUT) temp_store
    close(IOUT)
  endif ! ATTENUATION

  ! frees temporary memory
  deallocate(temp_store)

  ! dvpstore file output
  if (HETEROGEN_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
    call model_heterogen_mantle_output_dvp(prname)
  endif

  end subroutine save_model_meshfiles
