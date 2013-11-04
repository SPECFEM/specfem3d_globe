
  program attenuation_output

! BS This program is intended to only read in the values created in the files:
!    tau_s.bin
!    tau_e.bin
!    T_c_source.bin
!    Q.bin
!    It is dependent on the LOCAL_PATH, process, iregion, and the number of nodes
!    which reside on that process
!    Do not rely on this code for processing
!    However, the values here should match those which are output from attenuation_prem.c
!    This code was used for testing the implementation of 3D attenuation
!    Brian Savage, 22/03/04

    use constants_solver

    implicit none

    integer i,j,k,ispec
    integer myrank, vnspec, process, iregion
    character(len=150) prname, LOCAL_PATH
    double precision, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_AC)       :: one_minus_sum_beta, scale_factor
    double precision, dimension(NGLLX,NGLLY,NGLLZ,N_SLS,NSPEC_CRUST_MANTLE_AC) :: factor_common
    double precision, dimension(N_SLS)                                         :: tau_s
    double precision T_c_source

    LOCAL_PATH = "/scratch/DATABASES_MPI_BRIAN/att"
    process    = 42
    iregion    = IREGION_CRUST_MANTLE

    call create_name_database(prname, process, iregion, LOCAL_PATH)

    open(unit=27, file=prname(1:len_trim(prname))//'tau_s.bin',status='old',form='unformatted')
    read(27) tau_s
    close(27)

    open(unit=27, file=prname(1:len_trim(prname))//'T_c_source.bin',status='old',form='unformatted')
    read(27) T_c_source
    close(27);

    open(unit=27, file=prname(1:len_trim(prname))//'Q.bin',status='old',form='unformatted')
    read(27) scale_factor
    close(27)

    open(unit=27, file=prname(1:len_trim(prname))//'tau_e.bin',status='old',form='unformatted')
    read(27) factor_common
    close(27)

    write(*,*)' T_c_source = 1000.d0 / ', T_c_source
    write(*,*)' tau_sigma(1) = ', tau_s(1)
    write(*,*)' tau_sigma(2) = ', tau_s(2)
    write(*,*)' tau_sigma(3) = ', tau_s(3)

    do ispec = 1, NSPEC_CRUST_MANTLE_AC
       do k = 1, NGLLZ
          do j = 1, NGLLY
             do i = 1, NGLLX
                write(*,*)' tau_mu(1) = ', factor_common(i,j,k,1,ispec)
                write(*,*)' tau_mu(2) = ', factor_common(i,j,k,2,ispec)
                write(*,*)' tau_mu(3) = ', factor_common(i,j,k,3,ispec)
                write(*,*)' Qmu = ', scale_factor(i,j,k,ispec)
                write(*,*)
             enddo
          enddo
       enddo
    enddo
  end program attenuation_output


