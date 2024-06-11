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


  subroutine SIEM_compute_seismos()

  use specfem_par, only: SIMULATION_TYPE,ROTATION_VAL

  use specfem_par_crustmantle, only: displ_crust_mantle,veloc_crust_mantle

  use specfem_par_full_gravity, only: pgrav_cm, &
    seismograms_phi,seismograms_pgrav,seismograms_Hgrav,seismograms_grav,seismograms_corio

  implicit none

  select case (SIMULATION_TYPE)
  case (1)
    ! forward run
    ! perturbed gravitational potential
    call SIEM_compute_seismograms_phi(pgrav_cm,seismograms_phi)

    ! acceleration correction due to perturbed gravitational potential (\grad \phi)
    ! as well as gravity gradient for gravity strain measurements (\grad \grad \phi)
    call SIEM_compute_seismograms_pgrav(pgrav_cm,seismograms_pgrav,seismograms_Hgrav)

    ! acceleration correction due to background gravity
    ! s.\grad g (For vertical component, g=||g||)
    ! grad(g.s) (For vertical component)
    call SIEM_compute_seismograms_grav(displ_crust_mantle,seismograms_grav)

    if (ROTATION_VAL) then
      ! acceleration correction due to Coriolis acceleration (2\Omega x v)
      call SIEM_compute_seismograms_corio(veloc_crust_mantle,seismograms_corio)
    endif

  case (2)
    ! adjoint run
    ! nothing to do
    continue

  case (3)
    ! kernel run
    call SIEM_compute_seismograms_phi(pgrav_cm,seismograms_phi)

  end select

  end subroutine SIEM_compute_seismos

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_compute_seismograms_phi(phi_crust_mantle,seismograms)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ
  use constants_solver, only: NGLOB_CRUST_MANTLE

  use specfem_par, only: &
    nrec_local,ispec_selected_rec,number_receiver_global, &
    hxir_store,hetar_store,hgammar_store, &
    seismo_current,nlength_seismogram

  use specfem_par_crustmantle, only: ibool_crust_mantle

  use specfem_par_full_gravity, only: scale_pgrav

  implicit none

  real(kind=CUSTOM_REAL),dimension(1,NGLOB_CRUST_MANTLE),intent(in) :: phi_crust_mantle

  real(kind=CUSTOM_REAL),dimension(1,nrec_local,nlength_seismogram),intent(inout) :: seismograms

  ! local parameters
  ! phir = perturbed gravitational potential at receiver
  real(kind=CUSTOM_REAL) :: phir
  real(kind=CUSTOM_REAL) :: hlagrange
  integer :: i,j,k,ispec,iglob,irec_local,irec

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)
    ispec = ispec_selected_rec(irec)

    ! perform the general interpolation using Lagrange polynomials
    phir = 0.0_CUSTOM_REAL

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          hlagrange = real(hxir_store(i,irec_local)*hetar_store(j,irec_local)*hgammar_store(k,irec_local),kind=CUSTOM_REAL)
          phir = phir + phi_crust_mantle(1,iglob)*hlagrange
        enddo
      enddo
    enddo

    ! store North, East and Vertical components
    ! distinguish between single and double precision for reals
    seismograms(:,irec_local,seismo_current) = real(scale_pgrav * phir,kind=CUSTOM_REAL)
  enddo

  end subroutine SIEM_compute_seismograms_phi

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_compute_seismograms_pgrav(var_crust_mantle, seismograms, seismograms_Hgrav)

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use constants_solver, only: NGLOB_CRUST_MANTLE

  use specfem_par, only: hprime_xx, hprime_yy, hprime_zz, SAVE_SEISMOGRAMS_STRAIN

  use specfem_par, only: &
    nrec_local,nu_rec,ispec_selected_rec,number_receiver_global, &
    hxir_store,hetar_store,hgammar_store, &
    seismo_current,nlength_seismogram

  use specfem_par_crustmantle, only: ibool_crust_mantle, &
    xix_crust_mantle, xiy_crust_mantle, xiz_crust_mantle, &
    etax_crust_mantle, etay_crust_mantle, etaz_crust_mantle, &
    gammax_crust_mantle, gammay_crust_mantle, gammaz_crust_mantle

  use specfem_par_full_gravity, only: scale_pgrav

  implicit none

  real(kind=CUSTOM_REAL),dimension(NGLOB_CRUST_MANTLE),intent(in) :: var_crust_mantle

  real(kind=CUSTOM_REAL),dimension(NDIM,nrec_local,nlength_seismogram),intent(inout) :: seismograms
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM,nrec_local,nlength_seismogram),intent(inout) :: seismograms_Hgrav

  ! local parameters
  integer :: i,j,k,l,m, icomponent,ispec,irec_local,irec
  real(kind=CUSTOM_REAL) :: gradphi(NDIM), gll_gradphi(NDIM), gll_gstore(NDIM,NGLLX,NGLLY,NGLLZ)
  real(kind=CUSTOM_REAL) :: Hgrav_gll(NDIM,NDIM), Hgrav_rec(NDIM,NDIM)
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,hlagrange, &
                            tempx1l_phi, tempx2l_phi, tempx3l_phi
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) :: eps_loc_rot,eps_loc_new

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)
    ispec = ispec_selected_rec(irec)

    ! Initialise grad phi at receiver location
    gradphi(:) = 0.0_CUSTOM_REAL

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          ! Compute grad phi at GLL point
          tempx1l_phi   = 0._CUSTOM_REAL
          tempx2l_phi   = 0._CUSTOM_REAL
          tempx3l_phi   = 0._CUSTOM_REAL
          do l = 1,NGLLX
            tempx1l_phi = tempx1l_phi + var_crust_mantle(ibool_crust_mantle(l,j,k,ispec))*hprime_xx(i,l)
            tempx2l_phi = tempx2l_phi + var_crust_mantle(ibool_crust_mantle(i,l,k,ispec))*hprime_yy(j,l)
            tempx3l_phi = tempx3l_phi + var_crust_mantle(ibool_crust_mantle(i,j,l,ispec))*hprime_zz(k,l)
          enddo

          xixl    = xix_crust_mantle(i,j,k,ispec)
          xiyl    = xiy_crust_mantle(i,j,k,ispec)
          xizl    = xiz_crust_mantle(i,j,k,ispec)
          etaxl   = etax_crust_mantle(i,j,k,ispec)
          etayl   = etay_crust_mantle(i,j,k,ispec)
          etazl   = etaz_crust_mantle(i,j,k,ispec)
          gammaxl = gammax_crust_mantle(i,j,k,ispec)
          gammayl = gammay_crust_mantle(i,j,k,ispec)
          gammazl = gammaz_crust_mantle(i,j,k,ispec)

          gll_gradphi(1) = (xixl * tempx1l_phi) + (etaxl * tempx2l_phi) + (gammaxl * tempx3l_phi)
          gll_gradphi(2) = (xiyl * tempx1l_phi) + (etayl * tempx2l_phi) + (gammayl * tempx3l_phi)
          gll_gradphi(3) = (xizl * tempx1l_phi) + (etazl * tempx2l_phi) + (gammazl * tempx3l_phi)

          ! receiver interpolation
          hlagrange = real(hxir_store(i,irec_local)*hetar_store(j,irec_local)*hgammar_store(k,irec_local),kind=CUSTOM_REAL)

          do icomponent = 1,3
            gradphi(icomponent) = gradphi(icomponent) + gll_gradphi(icomponent)*hlagrange
          enddo !icomponent

          if (SAVE_SEISMOGRAMS_STRAIN) then
            ! Store the values of - \grad\phi for computation of gravity strain
            gll_gstore(:,i,j,k) = - gll_gradphi(:)
          endif

        enddo
      enddo
    enddo

    ! Calculate the gravity strain
    if (SAVE_SEISMOGRAMS_STRAIN) then
      ! interpolate the strain at actual receiver location
      Hgrav_rec(:,:) = 0.0_CUSTOM_REAL

      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX

            do m = 1,3 !component of g
              tempx1l_phi   = 0._CUSTOM_REAL
              tempx2l_phi   = 0._CUSTOM_REAL
              tempx3l_phi   = 0._CUSTOM_REAL
              do l = 1,NGLLX
                tempx1l_phi = tempx1l_phi + gll_gstore(m,l,j,k)*hprime_xx(i,l)
                tempx2l_phi = tempx2l_phi + gll_gstore(m,i,l,k)*hprime_yy(j,l)
                tempx3l_phi = tempx3l_phi + gll_gstore(m,i,j,l)*hprime_zz(k,l)
              enddo !l

              xixl    = xix_crust_mantle(i,j,k,ispec)
              xiyl    = xiy_crust_mantle(i,j,k,ispec)
              xizl    = xiz_crust_mantle(i,j,k,ispec)
              etaxl   = etax_crust_mantle(i,j,k,ispec)
              etayl   = etay_crust_mantle(i,j,k,ispec)
              etazl   = etaz_crust_mantle(i,j,k,ispec)
              gammaxl = gammax_crust_mantle(i,j,k,ispec)
              gammayl = gammay_crust_mantle(i,j,k,ispec)
              gammazl = gammaz_crust_mantle(i,j,k,ispec)

              ! Temporary Hgrav at each GLL point in the element
              Hgrav_gll(m,1) = (xixl * tempx1l_phi) + (etaxl * tempx2l_phi) + (gammaxl * tempx3l_phi) ! dg_m / dx
              Hgrav_gll(m,2) = (xiyl * tempx1l_phi) + (etayl * tempx2l_phi) + (gammayl * tempx3l_phi) ! dg_m / dy
              Hgrav_gll(m,3) = (xizl * tempx1l_phi) + (etazl * tempx2l_phi) + (gammazl * tempx3l_phi) ! dg_m / dz
            enddo

            ! Interpolate to the receiver location within element:
            hlagrange = real(hxir_store(i,irec_local)*hetar_store(j,irec_local)*hgammar_store(k,irec_local),kind=CUSTOM_REAL)

            Hgrav_rec(:,:) = Hgrav_rec(:,:) + Hgrav_gll(:,:)*hlagrange
          enddo
        enddo
      enddo
    endif


    ! store North, East and Vertical components

    ! distinguish between single and double precision for reals
    seismograms(:,irec_local,seismo_current) = real(scale_pgrav * (matmul(nu_rec(:,:,irec),gradphi)),kind=CUSTOM_REAL)

    if (SAVE_SEISMOGRAMS_STRAIN) then
      ! rotates from global x-y-z to the local coordinates (n-e-z):  eps_new = P*eps*P'
      ! nu is the rotation matrix from ECEF to local N-E-UP as defined.
      ! thus, if the nu is the rotation matrix that transforms coordinates from the global system (x,y,z) to the local
      ! coordinate system (N,E,V), e.g., a tensor is transformed as
      ! T_L = \nu * T_g * \nu^T
      !
      ! global -> local (n-e-up)
      ! eps_xx -> eps_nn
      ! eps_xy -> eps_ne
      ! eps_xz -> eps_nz
      ! eps_yx -> eps_en
      ! eps_yy -> eps_ee
      ! eps_yz -> eps_ez
      ! eps_zx -> eps_zn
      ! eps_zy -> eps_ze
      ! eps_zz -> eps_zz (z in radial direction up)
      eps_loc_rot(:,:) = real(matmul(nu_rec(:,:,irec), Hgrav_rec(:,:)),kind=CUSTOM_REAL)
      eps_loc_new(:,:) = real(matmul(eps_loc_rot(:,:), transpose(nu_rec(:,:,irec))),kind=CUSTOM_REAL)

      seismograms_Hgrav(:,:,irec_local,seismo_current) = real(scale_pgrav * eps_loc_new(:,:),kind=CUSTOM_REAL)

      ! original: misses P' projection for rotation
      !do m = 1,3
      !  seismogramsHgrav(:,m,irec_local,seismo_current) = &
      !                                        real(scale_pgrav*(matmul(nu_rec(:,:,irec),Hgrav_rec(:,m))),kind=CUSTOM_REAL)
      !enddo !m
    endif

  enddo !irec_local

  end subroutine SIEM_compute_seismograms_pgrav

!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_compute_seismograms_grav(var_crust_mantle, seismograms)

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,NGLLCUBE,CUSTOM_REAL
  use constants_solver, only: NGLOB_CRUST_MANTLE

  use specfem_par, only: &
    nrec_local,nu_rec,ispec_selected_rec,number_receiver_global, &
    hxir_store,hetar_store,hgammar_store, &
    seismo_current,nlength_seismogram

  use specfem_par_crustmantle, only: ibool_crust_mantle

  use specfem_par_full_gravity, only: scale_accel, &
    gradg_rec,g_spec_rec,storederiv_rec

  implicit none

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_CRUST_MANTLE),intent(in) :: var_crust_mantle

  real(kind=CUSTOM_REAL),dimension(NDIM,nrec_local,nlength_seismogram),intent(inout) :: seismograms

  ! local parameters
  integer :: i,j,k,igll,iglob,irec_local,irec,icomp,ispec
  real(kind=CUSTOM_REAL) :: av,ah1,ah2
  real(kind=CUSTOM_REAL) :: deriv(NDIM,NGLLCUBE),disp(NDIM),dispr(NDIM)
  real(kind=CUSTOM_REAL) :: g_dot_s(NGLLCUBE),grad_g_dot_s(NDIM)
  real(kind=CUSTOM_REAL) :: hlagrange

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)
    ispec = ispec_selected_rec(irec)

    ! perform the general interpolation using Lagrange polynomials
    dispr(:) = 0.0_CUSTOM_REAL
    g_dot_s(:) = 0.0_CUSTOM_REAL

    igll = 0
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          igll = igll+1
          iglob = ibool_crust_mantle(i,j,k,ispec)
          hlagrange = real(hxir_store(i,irec_local)*hetar_store(j,irec_local)*hgammar_store(k,irec_local),kind=CUSTOM_REAL)
          disp(:) = var_crust_mantle(:,iglob)
          g_dot_s(igll) = dot_product(g_spec_rec(:,igll,irec_local),disp)
          do icomp = 1,NDIM
            dispr(icomp) = dispr(icomp) + disp(icomp)*hlagrange
          enddo
        enddo
      enddo
    enddo

    ! vertical component
    av = dot_product(dispr,gradg_rec(:,irec_local))

    ! horizonal components
    deriv(:,:) = storederiv_rec(:,:,irec_local)
    grad_g_dot_s = matmul(deriv,g_dot_s)
    ! Note: nu is stored for nrec NOT for nrec_local
    ah1 = - real(dot_product(nu_rec(1,:,irec),grad_g_dot_s),kind=CUSTOM_REAL)
    ah2 = - real(dot_product(nu_rec(2,:,irec),grad_g_dot_s),kind=CUSTOM_REAL)

    ! store North, East and Vertical components
    ! distinguish between single and double precision for reals
    seismograms(1,irec_local,seismo_current) = real(scale_accel * ah1,kind=CUSTOM_REAL)
    seismograms(2,irec_local,seismo_current) = real(scale_accel * ah2,kind=CUSTOM_REAL)
    seismograms(3,irec_local,seismo_current) = real(scale_accel * av,kind=CUSTOM_REAL)
  enddo

  end subroutine SIEM_compute_seismograms_grav


!
!-------------------------------------------------------------------------------------------------
!

  subroutine SIEM_compute_seismograms_corio(var_crust_mantle, seismograms)

  use constants, only: NDIM,NGLLX,NGLLY,NGLLZ,CUSTOM_REAL
  use constants_solver, only: NGLOB_CRUST_MANTLE

  use specfem_par, only: two_omega_earth

  use specfem_par, only: &
    nrec_local,nu_rec,ispec_selected_rec,number_receiver_global, &
    hxir_store,hetar_store,hgammar_store, &
    seismo_current,nlength_seismogram

  use specfem_par_crustmantle, only: ibool_crust_mantle

  use specfem_par_full_gravity, only: scale_accel

  implicit none

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_CRUST_MANTLE),intent(in) :: var_crust_mantle

  real(kind=CUSTOM_REAL),dimension(NDIM,nrec_local,nlength_seismogram),intent(inout) :: seismograms

  ! local parameters
  real(kind=CUSTOM_REAL) :: acorio(NDIM),varx(NDIM)
  real(kind=CUSTOM_REAL) :: hlagrange
  integer :: i,j,k,ispec,iglob,irec_local,irec,icomp

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)
    ispec = ispec_selected_rec(irec)

    ! perform the general interpolation using Lagrange polynomials
    varx(:) = 0.0_CUSTOM_REAL

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          hlagrange = real(hxir_store(i,irec_local)*hetar_store(j,irec_local)*hgammar_store(k,irec_local),kind=CUSTOM_REAL)
          do icomp = 1,NDIM
            varx(icomp) = varx(icomp) + var_crust_mantle(icomp,iglob)*hlagrange
          enddo
        enddo
      enddo
    enddo
    ! store North, East and Vertical components

    ! \Omega x v
    acorio(:) = 0.d0
    acorio(1) = - two_omega_earth*varx(2)
    acorio(2) = two_omega_earth*varx(1)
    !acorio(3) = 0 ! because the axis of rotation is Z.

    ! distinguish between single and double precision for reals
    seismograms(:,irec_local,seismo_current) = real(scale_accel * (matmul(nu_rec(:,:,irec),acorio)),kind=CUSTOM_REAL)
  enddo

  end subroutine SIEM_compute_seismograms_corio

