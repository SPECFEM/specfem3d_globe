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

  subroutine save_regular_kernels_cm()

  use specfem_par
  use specfem_par_crustmantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: &
    cijkl_kl_crust_mantle_reg
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    rho_kl_crust_mantle_reg, beta_kl_crust_mantle_reg, alpha_kl_crust_mantle_reg
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: &
    mu_kl_crust_mantle, kappa_kl_crust_mantle, rhonotprime_kl_crust_mantle
  real(kind=CUSTOM_REAL),dimension(21) ::  cijkl_kl_local
  real(kind=CUSTOM_REAL) :: scale_kl,scale_kl_ani,scale_kl_rho
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal,rho_kl,alpha_kl,beta_kl
  real(kind=CUSTOM_REAL) :: alphah_kl,alphav_kl,betah_kl,betav_kl,rhonotprime_kl
  integer :: ispec,i,j,k,iglob
  double precision :: hlagrange
  integer :: ipoint

  ! transverse isotropic parameters
  real(kind=CUSTOM_REAL), dimension(21) :: an_kl
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
    betav_kl_crust_mantle,betah_kl_crust_mantle, &
    eta_kl_crust_mantle

  ! bulk parameterization
  real(kind=CUSTOM_REAL), dimension(:),allocatable :: &
    bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
    bulk_betav_kl_crust_mantle,bulk_betah_kl_crust_mantle
  real(kind=CUSTOM_REAL) :: A,C,F,L,N,eta
  real(kind=CUSTOM_REAL) :: muvl,kappavl,muhl,kappahl
  real(kind=CUSTOM_REAL) :: alphav_sq,alphah_sq,betav_sq,betah_sq,bulk_sq

  ! scaling factors
  scale_kl = scale_t/scale_displ * 1.d9
  ! For anisotropic kernels
  ! final unit : [s km^(-3) GPa^(-1)]
  scale_kl_ani = scale_t**3 / (RHOAV*R_EARTH**3) * 1.d18
  ! final unit : [s km^(-3) (kg/m^3)^(-1)]
  scale_kl_rho = scale_t / scale_displ / RHOAV * 1.d9

  ! allocates temporary arrays
  allocate(rho_kl_crust_mantle_reg(npoints_slice), &
           beta_kl_crust_mantle_reg(npoints_slice), &
           alpha_kl_crust_mantle_reg(npoints_slice))

  if (ANISOTROPIC_KL) then
    allocate(cijkl_kl_crust_mantle_reg(21, npoints_slice))
    if (SAVE_TRANSVERSE_KL_ONLY) then
      ! transverse isotropic kernel arrays for file output
      allocate(alphav_kl_crust_mantle(npoints_slice), &
               alphah_kl_crust_mantle(npoints_slice), &
               betav_kl_crust_mantle(npoints_slice), &
               betah_kl_crust_mantle(npoints_slice), &
               eta_kl_crust_mantle(npoints_slice))

      ! isotropic kernel arrays for file output
      allocate(bulk_c_kl_crust_mantle(npoints_slice), &
               bulk_betav_kl_crust_mantle(npoints_slice), &
               bulk_betah_kl_crust_mantle(npoints_slice), &
               bulk_beta_kl_crust_mantle(npoints_slice))
    endif
  else
    ! allocates temporary isotropic kernel arrays for file output
    allocate(bulk_c_kl_crust_mantle(npoints_slice), &
             bulk_beta_kl_crust_mantle(npoints_slice))
    allocate(mu_kl_crust_mantle(npoints_slice), &
             kappa_kl_crust_mantle(npoints_slice), &
             rhonotprime_kl_crust_mantle(npoints_slice))
  endif

  ! crust_mantle
  do ipoint = 1, npoints_slice
    ispec = ispec_reg(ipoint)
    rho_kl_crust_mantle_reg(ipoint) = 0.0
    alpha_kl_crust_mantle_reg(ipoint) = 0.0
    beta_kl_crust_mantle_reg(ipoint) = 0.0
    if (ANISOTROPIC_KL) then
      cijkl_kl_crust_mantle_reg(:,ipoint) = 0.0
      if (SAVE_TRANSVERSE_KL_ONLY) then
        alphav_kl_crust_mantle(ipoint) = 0.0
        alphah_kl_crust_mantle(ipoint) = 0.0
        betav_kl_crust_mantle(ipoint) = 0.0
        betah_kl_crust_mantle(ipoint) = 0.0
        eta_kl_crust_mantle(ipoint) = 0.0
        bulk_c_kl_crust_mantle(ipoint) = 0.0
        bulk_betav_kl_crust_mantle(ipoint) = 0.0
        bulk_betah_kl_crust_mantle(ipoint) = 0.0
        bulk_beta_kl_crust_mantle(ipoint) = 0.0
      endif
    else
      rhonotprime_kl_crust_mantle(ipoint) = 0.0
      kappa_kl_crust_mantle(ipoint) = 0.0
      mu_kl_crust_mantle(ipoint) = 0.0

      bulk_c_kl_crust_mantle(ipoint) = 0.0
      bulk_beta_kl_crust_mantle(ipoint) = 0.0
    endif

    do k = 1, NGLLZ
      do j = 1, NGLLY
        do i = 1, NGLLX

          hlagrange = hxir_reg(i,ipoint)*hetar_reg(j,ipoint)*hgammar_reg(k,ipoint)

          if (ANISOTROPIC_KL) then

            ! For anisotropic kernels
            iglob = ibool_crust_mantle(i,j,k,ispec)

            ! The Cartesian global cijkl_kl are rotated into the spherical local cijkl_kl
            ! ystore and zstore are thetaval and phival (line 2252) -- dangerous
            call rotate_kernels_dble(cijkl_kl_crust_mantle(:,i,j,k,ispec),cijkl_kl_local, &
                                     ystore_crust_mantle(iglob),zstore_crust_mantle(iglob))

            cijkl_kl_crust_mantle_reg(:,ipoint) = cijkl_kl_crust_mantle_reg(:,ipoint) &
                                                + cijkl_kl_local * scale_kl_ani * hlagrange
            rho_kl = rho_kl_crust_mantle(i,j,k,ispec) * scale_kl_rho * hlagrange

            ! transverse isotropic kernel calculations
            if (SAVE_TRANSVERSE_KL_ONLY) then
              ! note: transverse isotropic kernels are calculated for all elements
              !
              !          however, the factors A,C,L,N,F are based only on transverse elements
              !          in between Moho and 220 km layer, otherwise they are taken from isotropic values

              rhol = rhostore_crust_mantle(i,j,k,ispec)

              ! transverse isotropic parameters from compute_force_crust_mantle.f90
              ! C=rhovpvsq A=rhovphsq L=rhovsvsq N=rhovshsq eta=F/(A - 2 L)

              ! Get A,C,F,L,N,eta from kappa,mu
              ! element can have transverse isotropy if between d220 and Moho
              if (.not. ispec_is_tiso_crust_mantle(ispec)) then

                ! layer with no transverse isotropy
                ! A,C,L,N,F from isotropic model

                mul = muvstore_crust_mantle(i,j,k,ispec)
                kappal = kappavstore_crust_mantle(i,j,k,ispec)
                muvl = mul
                muhl = mul

                A = kappal + FOUR_THIRDS * mul
                C = A
                L = mul
                N = mul
                F = kappal - 2._CUSTOM_REAL/3._CUSTOM_REAL * mul
                eta = 1._CUSTOM_REAL

              else

                ! A,C,L,N,F from transverse isotropic model
                kappavl = kappavstore_crust_mantle(i,j,k,ispec)
                kappahl = kappahstore_crust_mantle(i,j,k,ispec)
                muvl = muvstore_crust_mantle(i,j,k,ispec)
                muhl = muhstore_crust_mantle(i,j,k,ispec)
                kappal = kappavl

                A = kappahl + FOUR_THIRDS * muhl
                C = kappavl + FOUR_THIRDS * muvl
                L = muvl
                N = muhl
                eta = eta_anisostore_crust_mantle(i,j,k,ispec)  ! that is  F / (A - 2 L)
                F = eta * ( A - 2._CUSTOM_REAL * L )

              endif

              ! note: cijkl_kl_local() is fully anisotropic C_ij kernel components (non-dimensionalized)
              !          for GLL point at (i,j,k,ispec)

              ! Purpose : compute the kernels for the An coeffs (an_kl)
              ! from the kernels for Cij (cijkl_kl_local)
              ! At r,theta,phi fixed
              ! kernel def : dx = kij * dcij + krho * drho
              !                = kAn * dAn  + krho * drho

              ! Definition of the input array cij_kl :
              ! cij_kl(1) = C11 ; cij_kl(2) = C12 ; cij_kl(3) = C13
              ! cij_kl(4) = C14 ; cij_kl(5) = C15 ; cij_kl(6) = C16
              ! cij_kl(7) = C22 ; cij_kl(8) = C23 ; cij_kl(9) = C24
              ! cij_kl(10) = C25 ; cij_kl(11) = C26 ; cij_kl(12) = C33
              ! cij_kl(13) = C34 ; cij_kl(14) = C35 ; cij_kl(15) = C36
              ! cij_kl(16) = C44 ; cij_kl(17) = C45 ; cij_kl(18) = C46
              ! cij_kl(19) = C55 ; cij_kl(20) = C56 ; cij_kl(21) = C66
              ! where the Cij (Voigt's notation) are defined as function of
              ! the components of the elastic tensor in spherical coordinates
              ! by eq. (A.1) of Chen & Tromp, GJI 168 (2007)

              ! From the relations giving Cij in function of An
              ! Checked with Min Chen's results (routine build_cij)

              an_kl(1) = cijkl_kl_local(1)+cijkl_kl_local(2)+cijkl_kl_local(7)  !A
              an_kl(2) = cijkl_kl_local(12)                                     !C
              an_kl(3) = -2*cijkl_kl_local(2)+cijkl_kl_local(21)                !N
              an_kl(4) = cijkl_kl_local(16)+cijkl_kl_local(19)                  !L
              an_kl(5) = cijkl_kl_local(3)+cijkl_kl_local(8)                    !F

              ! not used yet
              !an_kl(6)=2*cijkl_kl_local(5)+2*cijkl_kl_local(10)+2*cijkl_kl_local(14)          !Jc
              !an_kl(7)=2*cijkl_kl_local(4)+2*cijkl_kl_local(9)+2*cijkl_kl_local(13)           !Js
              !an_kl(8)=-2*cijkl_kl_local(14)                                  !Kc
              !an_kl(9)=-2*cijkl_kl_local(13)                                  !Ks
              !an_kl(10)=-2*cijkl_kl_local(10)+cijkl_kl_local(18)                      !Mc
              !an_kl(11)=2*cijkl_kl_local(4)-cijkl_kl_local(20)                        !Ms
              !an_kl(12)=cijkl_kl_local(1)-cijkl_kl_local(7)                           !Bc
              !an_kl(13)=-1./2.*(cijkl_kl_local(6)+cijkl_kl_local(11))                 !Bs
              !an_kl(14)=cijkl_kl_local(3)-cijkl_kl_local(8)                           !Hc
              !an_kl(15)=-cijkl_kl_local(15)                                   !Hs
              !an_kl(16)=-cijkl_kl_local(16)+cijkl_kl_local(19)                        !Gc
              !an_kl(17)=-cijkl_kl_local(17)                                   !Gs
              !an_kl(18)=cijkl_kl_local(5)-cijkl_kl_local(10)-cijkl_kl_local(18)               !Dc
              !an_kl(19)=cijkl_kl_local(4)-cijkl_kl_local(9)+cijkl_kl_local(20)                !Ds
              !an_kl(20)=cijkl_kl_local(1)-cijkl_kl_local(2)+cijkl_kl_local(7)-cijkl_kl_local(21)      !Ec
              !an_kl(21)=-cijkl_kl_local(6)+cijkl_kl_local(11)                         !Es

              ! K_rho (primary kernel, for a parameterization (A,C,L,N,F,rho) )
              rhonotprime_kl = rhol * rho_kl / scale_kl_rho
              !rhonotprime_kl_crust_mantle(ipoint) = rhonotprime_kl_crust_mantle(ipoint) + rhonotprime_kl

              ! note: transverse isotropic kernels are calculated for ALL elements,
              !          and not just transverse elements
              !
              ! note: the kernels are for relative perturbations (delta ln (m_i) = (m_i - m_0)/m_i )
              !
              ! Gets transverse isotropic kernels
              ! (see Appendix B of Sieminski et al., GJI 171, 2007)

              ! for parameterization: ( alpha_v, alpha_h, beta_v, beta_h, eta, rho )
              ! K_alpha_v
              alphav_kl = 2*C*an_kl(2) * hlagrange
              alphav_kl_crust_mantle(ipoint) = alphav_kl_crust_mantle(ipoint) &
                                             + alphav_kl
              ! K_alpha_h
              alphah_kl = (2*A*an_kl(1) + 2*A*eta*an_kl(5)) * hlagrange
              alphah_kl_crust_mantle(ipoint) = alphah_kl_crust_mantle(ipoint) &
                                             + alphah_kl
              ! K_beta_v
              betav_kl = (2*L*an_kl(4) - 4*L*eta*an_kl(5)) * hlagrange
              betav_kl_crust_mantle(ipoint) = betav_kl_crust_mantle(ipoint) &
                                            + betav_kl
              ! K_beta_h
              betah_kl = 2*N*an_kl(3) * hlagrange
              betah_kl_crust_mantle(ipoint) = betah_kl_crust_mantle(ipoint) &
                                            + betah_kl
              ! K_eta
              eta_kl_crust_mantle(ipoint) = eta_kl_crust_mantle(ipoint) + F*an_kl(5) * hlagrange
              ! K_rhoprime  (for a parameterization (alpha_v, ..., rho) )
              rho_kl_crust_mantle_reg(ipoint) = rho_kl_crust_mantle_reg(ipoint) &
                                              + (A*an_kl(1) + C*an_kl(2) &
                                               + N*an_kl(3) + L*an_kl(4) + F*an_kl(5)) * hlagrange &
                                              + rhonotprime_kl

              ! for parameterization: ( bulk, beta_v, beta_h, eta, rho )
              ! where kappa_v = kappa_h = kappa and bulk c = sqrt( kappa / rho )
              betav_sq = muvl / rhol
              betah_sq = muhl / rhol
              alphav_sq = ( kappal + FOUR_THIRDS * muvl ) / rhol
              alphah_sq = ( kappal + FOUR_THIRDS * muhl ) / rhol
              bulk_sq = kappal / rhol

              bulk_c_kl_crust_mantle(ipoint) = bulk_c_kl_crust_mantle(ipoint) + &
                bulk_sq / alphav_sq * alphav_kl + &
                bulk_sq / alphah_sq * alphah_kl

              bulk_betah_kl_crust_mantle(ipoint) = bulk_betah_kl_crust_mantle(ipoint) + &
                betah_kl + FOUR_THIRDS * betah_sq / alphah_sq * alphah_kl

              bulk_betav_kl_crust_mantle(ipoint) = bulk_betav_kl_crust_mantle(ipoint) + &
                betav_kl + FOUR_THIRDS * betav_sq / alphav_sq * alphav_kl
              ! the rest, K_eta and K_rho are the same as above
            else

              rho_kl_crust_mantle_reg(ipoint) = rho_kl_crust_mantle_reg(ipoint) + rho_kl

            endif ! SAVE_TRANSVERSE_KL_ONLY

          else

            ! isotropic kernels

            rhol = rhostore_crust_mantle(i,j,k,ispec)
            mul = muvstore_crust_mantle(i,j,k,ispec)
            kappal = kappavstore_crust_mantle(i,j,k,ispec)

            ! kernel values for rho, kappa, mu (primary kernel values)
            rho_kl = - rhol * rho_kl_crust_mantle(i,j,k,ispec)
            alpha_kl = - kappal * alpha_kl_crust_mantle(i,j,k,ispec) ! note: alpha_kl corresponds to K_kappa
            beta_kl =  - 2 * mul * beta_kl_crust_mantle(i,j,k,ispec) ! note: beta_kl corresponds to K_mu

            ! for a parameterization: (rho,mu,kappa) "primary" kernels
            rhonotprime_kl_crust_mantle(ipoint) = rhonotprime_kl_crust_mantle(ipoint) &
                                                + rho_kl * scale_kl * hlagrange
            mu_kl_crust_mantle(ipoint) = mu_kl_crust_mantle(ipoint) + beta_kl * scale_kl * hlagrange
            kappa_kl_crust_mantle(ipoint) = kappa_kl_crust_mantle(ipoint) + alpha_kl * scale_kl * hlagrange

            ! for a parameterization: (rho,alpha,beta)
            ! kernels rho^prime, beta, alpha
            rho_kl_crust_mantle_reg(ipoint) = rho_kl_crust_mantle_reg(ipoint) &
                                            + (rho_kl + alpha_kl + beta_kl) * scale_kl * hlagrange
            beta_kl_crust_mantle_reg(ipoint) = beta_kl_crust_mantle_reg(ipoint) + &
              2._CUSTOM_REAL * (beta_kl - FOUR_THIRDS * mul * alpha_kl / kappal) * scale_kl * hlagrange
            alpha_kl_crust_mantle_reg(ipoint) = alpha_kl_crust_mantle_reg(ipoint) + &
              2._CUSTOM_REAL * (1._CUSTOM_REAL +  FOUR_THIRDS * mul / kappal) * alpha_kl * scale_kl * hlagrange

            ! for a parameterization: (rho,bulk, beta)
            ! where bulk wave speed is c = sqrt( kappa / rho)
            ! note: rhoprime is the same as for (rho,alpha,beta) parameterization
            bulk_c_kl_crust_mantle(ipoint) = bulk_c_kl_crust_mantle(ipoint) &
                                           + 2._CUSTOM_REAL * alpha_kl * scale_kl * hlagrange
            bulk_beta_kl_crust_mantle(ipoint) = bulk_beta_kl_crust_mantle(ipoint) &
                                              + 2._CUSTOM_REAL * beta_kl * scale_kl * hlagrange

          endif

        enddo
      enddo
    enddo

    ! do some transforms that are independent of GLL points
    if (ANISOTROPIC_KL) then
      if (SAVE_TRANSVERSE_KL_ONLY) then
        ! write the kernel in physical units
        !rhonotprime_kl_crust_mantle(ipoint) = - rhonotprime_kl_crust_mantle(ipoint) * scale_kl

        alphav_kl_crust_mantle(ipoint) = - alphav_kl_crust_mantle(ipoint) * scale_kl
        alphah_kl_crust_mantle(ipoint) = - alphah_kl_crust_mantle(ipoint) * scale_kl
        betav_kl_crust_mantle(ipoint) = - betav_kl_crust_mantle(ipoint) * scale_kl
        betah_kl_crust_mantle(ipoint) = - betah_kl_crust_mantle(ipoint) * scale_kl
        eta_kl_crust_mantle(ipoint) = - eta_kl_crust_mantle(ipoint) * scale_kl
        rho_kl_crust_mantle_reg(ipoint) = - rho_kl_crust_mantle_reg(ipoint) * scale_kl

        ! to check: isotropic kernels from transverse isotropic ones
        alpha_kl_crust_mantle_reg(ipoint) = alphav_kl_crust_mantle(ipoint) &
                                          + alphah_kl_crust_mantle(ipoint)
        beta_kl_crust_mantle_reg(ipoint) = betav_kl_crust_mantle(ipoint) &
                                         + betah_kl_crust_mantle(ipoint)
        !rho_kl_crust_mantle_reg(ipoint) = rho_kl_crust_mantle_reg(ipoint) &
        !                                + rhonotprime_kl_crust_mantle(ipoint) &
        !                                + alpha_kl_crust_mantle_reg(ipoint) &
        !                                + beta_kl_crust_mantle_reg(ipoint)
        bulk_beta_kl_crust_mantle(ipoint) = bulk_betah_kl_crust_mantle(ipoint) &
                                          + bulk_betav_kl_crust_mantle(ipoint)
      endif
    endif
  enddo

  ! writes out kernels to file
  if (ADIOS_FOR_KERNELS) then
    ! check implementation
    call exit_mpi(myrank,'saving regular kernels in ADIOS file format is not supported yet')
  else
    ! sets up database name
    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

    ! For anisotropic kernels
    if (ANISOTROPIC_KL) then

      ! outputs transverse isotropic kernels only
      if (SAVE_TRANSVERSE_KL_ONLY) then
        ! transverse isotropic kernels
        ! (alpha_v, alpha_h, beta_v, beta_h, eta, rho ) parameterization
        open(unit=IOUT,file=trim(prname)//'alphav_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) alphav_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'alphah_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) alphah_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'betav_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) betav_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'betah_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) betah_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'eta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) eta_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) rho_kl_crust_mantle_reg
        close(IOUT)

        ! in case one is interested in primary kernel K_rho
        !open(unit=IOUT,file=trim(prname)//'rhonotprime_kernel.bin',status='unknown',form='unformatted',action='write')
        !write(IOUT) rhonotprime_kl_crust_mantle
        !close(IOUT)

        ! (bulk, beta_v, beta_h, eta, rho ) parameterization: K_eta and K_rho same as above
        open(unit=IOUT,file=trim(prname)//'bulk_c_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) bulk_c_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'bulk_betav_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) bulk_betav_kl_crust_mantle
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'bulk_betah_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) bulk_betah_kl_crust_mantle
        close(IOUT)

        ! to check: isotropic kernels
        open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) alpha_kl_crust_mantle_reg
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) beta_kl_crust_mantle_reg
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'bulk_beta_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) bulk_beta_kl_crust_mantle
        close(IOUT)

      else

        ! fully anisotropic kernels
        ! note: the C_ij and density kernels are not for relative perturbations (delta ln( m_i) = delta m_i / m_i),
        !          but absolute perturbations (delta m_i = m_i - m_0)
        open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) - rho_kl_crust_mantle_reg
        close(IOUT)
        open(unit=IOUT,file=trim(prname)//'cijkl_kernel.bin',status='unknown',form='unformatted',action='write')
        write(IOUT) - cijkl_kl_crust_mantle_reg
        close(IOUT)

      endif

    else
      ! primary kernels: (rho,kappa,mu) parameterization
      open(unit=IOUT,file=trim(prname)//'rhonotprime_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) rhonotprime_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'kappa_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) kappa_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'mu_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) mu_kl_crust_mantle
      close(IOUT)

      ! (rho, alpha, beta ) parameterization
      open(unit=IOUT,file=trim(prname)//'rho_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) rho_kl_crust_mantle_reg
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'alpha_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) alpha_kl_crust_mantle_reg
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'beta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) beta_kl_crust_mantle_reg
      close(IOUT)

      ! (rho, bulk, beta ) parameterization, K_rho same as above
      open(unit=IOUT,file=trim(prname)//'bulk_c_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_c_kl_crust_mantle
      close(IOUT)
      open(unit=IOUT,file=trim(prname)//'bulk_beta_kernel.bin',status='unknown',form='unformatted',action='write')
      write(IOUT) bulk_beta_kl_crust_mantle
      close(IOUT)

    endif

  endif ! ADIOS_FOR_KERNELS

  ! cleans up temporary kernel arrays
  if (ANISOTROPIC_KL) then
    deallocate(cijkl_kl_crust_mantle_reg)
    if (SAVE_TRANSVERSE_KL_ONLY) then
      deallocate(alphav_kl_crust_mantle,alphah_kl_crust_mantle, &
                 betav_kl_crust_mantle,betah_kl_crust_mantle, &
                 eta_kl_crust_mantle)
      deallocate(bulk_c_kl_crust_mantle,bulk_betah_kl_crust_mantle, &
                 bulk_betav_kl_crust_mantle,bulk_beta_kl_crust_mantle)
    endif
  else
    deallocate(bulk_c_kl_crust_mantle,bulk_beta_kl_crust_mantle, &
               mu_kl_crust_mantle,kappa_kl_crust_mantle, &
               rhonotprime_kl_crust_mantle)
  endif
  deallocate(rho_kl_crust_mantle_reg, &
             beta_kl_crust_mantle_reg, &
             alpha_kl_crust_mantle_reg)

  end subroutine save_regular_kernels_cm

