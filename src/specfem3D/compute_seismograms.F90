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


  subroutine compute_seismograms(nglob,displ,seismo_current,seismograms)

  use constants_solver

  use specfem_par, only: &
    nlength_seismogram, &
    nrec_local,nu_rec,ispec_selected_rec,number_receiver_global, &
    scale_displ,hxir_store,hetar_store,hgammar_store

  use specfem_par_crustmantle, only: ibool_crust_mantle

  implicit none

  integer,intent(in) :: nglob
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ

  integer,intent(in) :: seismo_current

  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,nlength_seismogram),intent(out) :: &
    seismograms

  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  integer :: ispec,iglob,irec_local,irec
  integer :: i,j,k

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ispec = ispec_selected_rec(irec)

    ! perform the general interpolation using Lagrange polynomials
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO

    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)

          hlagrange = hxir_store(i,irec_local) * hetar_store(j,irec_local) * hgammar_store(k,irec_local)

          uxd = uxd + dble(displ(1,iglob))*hlagrange
          uyd = uyd + dble(displ(2,iglob))*hlagrange
          uzd = uzd + dble(displ(3,iglob))*hlagrange

        enddo
      enddo
    enddo

    ! store North, East and Vertical components
    ! distinguish between single and double precision for reals
    seismograms(:,irec_local,seismo_current) = real(scale_displ*(nu_rec(:,1,irec)*uxd + &
                                                                 nu_rec(:,2,irec)*uyd + &
                                                                 nu_rec(:,3,irec)*uzd), &
                                                    kind=CUSTOM_REAL)

  enddo

  end subroutine compute_seismograms

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_seismograms_adjoint(displ_crust_mantle, &
                                         eps_trace_over_3_crust_mantle, &
                                         epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
                                         epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle, &
                                         moment_der,sloc_der,stshift_der,shdur_der, &
                                         seismograms)

  use constants_solver, only: &
    CUSTOM_REAL,ZERO,PI,NGLLX,NGLLY,NGLLZ, &
    NDIM,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
    NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_CRUST_MANTLE_STR_OR_ATT

  use specfem_par, only: &
    NSTEP,NTSTEP_BETWEEN_OUTPUT_SAMPLE, &
    nlength_seismogram,seismo_current, &
    UNDO_ATTENUATION, &
    nrec_local, &
    nu_source,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
    hxir_store,hpxir_store,hetar_store,hpetar_store,hgammar_store,hpgammar_store, &
    tshift_src,hdur_Gaussian, &
    DT,t0,deltat,it, &
    scale_displ,scale_t, &
    hprime_xx,hprime_yy,hprime_zz, &
    ispec_selected_source,number_receiver_global

  use specfem_par_crustmantle, only: ibool_crust_mantle, &
    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

  implicit none

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE),intent(in) :: &
    displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STRAIN_ONLY),intent(in) :: &
    eps_trace_over_3_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STR_OR_ATT),intent(in) :: &
    epsilondev_xx_crust_mantle,epsilondev_yy_crust_mantle,epsilondev_xy_crust_mantle, &
    epsilondev_xz_crust_mantle,epsilondev_yz_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local),intent(inout) :: moment_der
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local),intent(inout) :: sloc_der
  real(kind=CUSTOM_REAL), dimension(nrec_local),intent(inout) :: stshift_der, shdur_der

  real(kind=CUSTOM_REAL), dimension(NDIM*NDIM,nrec_local,nlength_seismogram),intent(out) :: &
    seismograms

  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  double precision :: eps_trace,dxx,dyy,dxy,dxz,dyz
  double precision :: eps_loc(NDIM,NDIM), eps_loc_new(NDIM,NDIM)
  double precision :: stf,timeval

  real(kind=CUSTOM_REAL) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ)
  real(kind=CUSTOM_REAL) :: eps_s(NDIM,NDIM), eps_m_s, &
        eps_m_l_s(NDIM), stf_deltat, Kp_deltat, Hp_deltat
  integer :: iglob,irec_local,irec,ispec

  double precision, external :: get_stf_viscoelastic

  ! element strain
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_matrix

  double precision :: hxir(NGLLX), hetar(NGLLY), hgammar(NGLLZ), &
                      hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)
  integer :: i,j,k

  do irec_local = 1,nrec_local

    ! gets global number of that receiver
    irec = number_receiver_global(irec_local)

    ! gets element id
    ispec = ispec_selected_source(irec)

    ! adjoint strain
    if (UNDO_ATTENUATION) then
      ! recomputes strain
      call compute_element_strain_undoatt_noDev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                                                displ_crust_mantle, &
                                                hprime_xx,hprime_yy,hprime_zz,ibool_crust_mantle, &
                                                xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                                                etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                                                gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                                                epsilondev_loc_matrix,eps_trace_over_3_loc)
    else
      ! element adjoint strain
      eps_trace_over_3_loc(:,:,:) = eps_trace_over_3_crust_mantle(:,:,:,ispec)
      epsilondev_loc_matrix(1,:,:,:) = epsilondev_xx_crust_mantle(:,:,:,ispec)
      epsilondev_loc_matrix(2,:,:,:) = epsilondev_yy_crust_mantle(:,:,:,ispec)
      epsilondev_loc_matrix(3,:,:,:) = epsilondev_xy_crust_mantle(:,:,:,ispec)
      epsilondev_loc_matrix(4,:,:,:) = epsilondev_xz_crust_mantle(:,:,:,ispec)
      epsilondev_loc_matrix(5,:,:,:) = epsilondev_yz_crust_mantle(:,:,:,ispec)
    endif

    ! initializes
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO
    eps_trace = ZERO
    dxx = ZERO
    dyy = ZERO
    dxy = ZERO
    dxz = ZERO
    dyz = ZERO

    ! perform the general interpolation using Lagrange polynomials
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX

          iglob = ibool_crust_mantle(i,j,k,ispec)

          hlagrange = hxir_store(i,irec_local) * hetar_store(j,irec_local) * hgammar_store(k,irec_local)

          uxd = uxd + dble(displ_crust_mantle(1,iglob))*hlagrange
          uyd = uyd + dble(displ_crust_mantle(2,iglob))*hlagrange
          uzd = uzd + dble(displ_crust_mantle(3,iglob))*hlagrange

          eps_trace = eps_trace + dble(eps_trace_over_3_loc(i,j,k))*hlagrange

          dxx = dxx + dble(epsilondev_loc_matrix(1,i,j,k))*hlagrange
          dyy = dyy + dble(epsilondev_loc_matrix(2,i,j,k))*hlagrange
          dxy = dxy + dble(epsilondev_loc_matrix(3,i,j,k))*hlagrange
          dxz = dxz + dble(epsilondev_loc_matrix(4,i,j,k))*hlagrange
          dyz = dyz + dble(epsilondev_loc_matrix(5,i,j,k))*hlagrange

          displ_s(:,i,j,k) = displ_crust_mantle(:,iglob)
        enddo
      enddo
    enddo

    eps_loc(1,1) = eps_trace + dxx
    eps_loc(2,2) = eps_trace + dyy
    eps_loc(3,3) = eps_trace - dxx - dyy
    eps_loc(1,2) = dxy
    eps_loc(1,3) = dxz
    eps_loc(2,3) = dyz
    eps_loc(2,1) = dxy
    eps_loc(3,1) = dxz
    eps_loc(3,2) = dyz

    ! un-rotated
    !eps_loc_new(:,:) = eps_loc(:,:)
    !
    ! rotate to the local Cartesian coordinates (n-e-z):  eps_new = P*eps*P'
    eps_loc_new(:,:) = matmul(matmul(nu_source(:,:,irec),eps_loc(:,:)), transpose(nu_source(:,:,irec)))

    ! distinguish between single and double precision for reals
    seismograms(1,irec_local,seismo_current) = real(eps_loc_new(1,1), kind=CUSTOM_REAL)
    seismograms(2,irec_local,seismo_current) = real(eps_loc_new(2,2), kind=CUSTOM_REAL)
    seismograms(3,irec_local,seismo_current) = real(eps_loc_new(3,3), kind=CUSTOM_REAL)
    seismograms(4,irec_local,seismo_current) = real(eps_loc_new(1,2), kind=CUSTOM_REAL)
    seismograms(5,irec_local,seismo_current) = real(eps_loc_new(1,3), kind=CUSTOM_REAL)
    seismograms(6,irec_local,seismo_current) = real(eps_loc_new(2,3), kind=CUSTOM_REAL)

    seismograms(7:9,irec_local,seismo_current) = real(scale_displ*(nu_source(:,1,irec)*uxd + &
                                                                      nu_source(:,2,irec)*uyd + &
                                                                      nu_source(:,3,irec)*uzd), &
                                                                      kind=CUSTOM_REAL)

    ! interpolators
    ! note: we explicitly copy the store arrays to local temporary arrays here
    !       array indexing (irec_local,:) would be non-contiguous and compilers would have to do this anyway.
    !       since we now use (:,irec_local), we could however skip this... still, we keep it as an explicit way.
    hxir(:) = hxir_store(:,irec_local)
    hetar(:) = hetar_store(:,irec_local)
    hgammar(:) = hgammar_store(:,irec_local)
    hpxir(:) = hpxir_store(:,irec_local)
    hpetar(:) = hpetar_store(:,irec_local)
    hpgammar(:) = hpgammar_store(:,irec_local)

    ! Frechet derivatives of the source
    call compute_adj_source_frechet(displ_s,Mxx(irec),Myy(irec),Mzz(irec), &
                            Mxy(irec),Mxz(irec),Myz(irec), &
                            eps_s,eps_m_s,eps_m_l_s, &
                            hxir,hetar,hgammar, &
                            hpxir,hpetar,hpgammar, &
                            hprime_xx,hprime_yy,hprime_zz, &
                            xix_crust_mantle(1,1,1,ispec),xiy_crust_mantle(1,1,1,ispec),xiz_crust_mantle(1,1,1,ispec), &
                            etax_crust_mantle(1,1,1,ispec),etay_crust_mantle(1,1,1,ispec),etaz_crust_mantle(1,1,1,ispec), &
                            gammax_crust_mantle(1,1,1,ispec),gammay_crust_mantle(1,1,1,ispec),gammaz_crust_mantle(1,1,1,ispec))

    ! reverse time
    timeval = dble(NSTEP-it) * DT - t0 - tshift_src(irec)

    ! gets source-time function value
    !#TODO: double-check if time at it or NSTEP-it+1
    !       since adjoint simulation uses reversed adjoint sources, source time seems to correspond to NSTEP-(it-1)==NSTEP-it+1
    !       therefore, when reading in an external STF, the index would be `NSTEP-it+1` rather than `it`
    stf = get_stf_viscoelastic(timeval,irec,NSTEP-it+1)

    stf_deltat = real(stf * deltat * NTSTEP_BETWEEN_OUTPUT_SAMPLE,kind=CUSTOM_REAL)

    ! moment derivatives
    moment_der(:,:,irec_local) = moment_der(:,:,irec_local) + eps_s(:,:) * stf_deltat

    ! location derivatives
    sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_l_s(:) * stf_deltat

    ! derivatives for time shift and hduration
    Kp_deltat = real(-1.0d0/sqrt(PI)/hdur_Gaussian(irec) &
                     * exp(-(timeval/hdur_Gaussian(irec))**2) * deltat * scale_t,kind=CUSTOM_REAL)
    Hp_deltat = real(timeval/hdur_Gaussian(irec) * Kp_deltat,kind=CUSTOM_REAL)

    stshift_der(irec_local) = stshift_der(irec_local) + eps_m_s * Kp_deltat
    shdur_der(irec_local) = shdur_der(irec_local) + eps_m_s * Hp_deltat

  enddo

  end subroutine compute_seismograms_adjoint


!
!-------------------------------------------------------------------------------------------------
!

! unused...
!
!  subroutine compute_seismograms_undoatt()
!
!! re-orders seismogram entries
!
!  use specfem_par, only: CUSTOM_REAL,NDIM, &
!    NT_DUMP_ATTENUATION,NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
!    nrec_local,myrank, &
!    seismo_current,seismograms
!
!  implicit none
!
!  ! local parameters
!  integer :: i,j,k,irec_local
!  real(kind=CUSTOM_REAL), dimension(3) :: seismograms_temp
!
!  ! checks if anything to do
!  if (nrec_local == 0 ) return
!
!  if (mod(NT_DUMP_ATTENUATION,2) == 0) then
!
!    do irec_local = 1,nrec_local
!      do i = 1,seismo_current/NT_DUMP_ATTENUATION
!        do j = 1,NT_DUMP_ATTENUATION/2
!          do k = 1,NDIM
!            seismograms_temp(k) = seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + j)
!
!            seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + j) = &
!                          seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + (NT_DUMP_ATTENUATION-j+1))
!
!            seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + (NT_DUMP_ATTENUATION-j+1)) = seismograms_temp(k)
!          enddo
!        enddo
!      enddo
!    enddo
!
!  else
!
!    do irec_local = 1,nrec_local
!      do i = 1,seismo_current/NT_DUMP_ATTENUATION
!        do j = 1,(NT_DUMP_ATTENUATION-1)/2
!          do k = 1,NDIM
!            seismograms_temp(k) = seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + j)
!            seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + j) = &
!                  seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + (NT_DUMP_ATTENUATION-j+1))
!            seismograms(k,irec_local,(i-1)*NT_DUMP_ATTENUATION + (NT_DUMP_ATTENUATION-j+1)) = seismograms_temp(k)
!          enddo
!        enddo
!      enddo
!    enddo
!
!  endif
!
!  end subroutine compute_seismograms_undoatt

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_seismograms_strain(nglob,displ)

  use constants, only: CUSTOM_REAL,NDIM,NGLLX,NGLLY,NGLLZ

  use specfem_par, only: &
    nrec_local,nu_rec,ispec_selected_rec,number_receiver_global, &
    hxir_store,hetar_store,hgammar_store, &
    hprime_xx,hprime_yy,hprime_zz, &
    seismograms_eps, &
    seismo_current

  use specfem_par_crustmantle, only: &
    ibool_crust_mantle, &
    xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
    etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
    gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

  implicit none

  ! input
  integer,intent(in) :: nglob
  real(kind=CUSTOM_REAL),dimension(NDIM,nglob),intent(in) :: displ

  ! local parameters
  double precision :: hxir(NGLLX), hetar(NGLLY), hgammar(NGLLZ)
  real(kind=CUSTOM_REAL) :: tempx1l,tempx2l,tempx3l,tempy1l,tempy2l,tempy3l,tempz1l,tempz2l,tempz3l
  real(kind=CUSTOM_REAL) :: hp1,hp2,hp3
  real(kind=CUSTOM_REAL) :: xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl
  real(kind=CUSTOM_REAL) :: duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl
  real(kind=CUSTOM_REAL) :: eps_xx,eps_yy,eps_zz,eps_xy,eps_xz,eps_yz
  real(kind=CUSTOM_REAL) :: hlagrange

  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) :: eps_loc,eps_loc_new
  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM) :: eps_loc_rot

  real(kind=CUSTOM_REAL),dimension(NDIM,NDIM,NGLLX,NGLLY,NGLLZ) :: eps_array
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: displ_elem

  integer :: irec_local,irec
  integer :: iglob,ispec,i,j,k,l

  do irec_local = 1,nrec_local

    ! gets global number of that receiver
    irec = number_receiver_global(irec_local)

    ! gets element id
    ispec = ispec_selected_rec(irec)

    ! fetches local element displacements
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          displ_elem(:,i,j,k) = displ(:,iglob)
        enddo
      enddo
    enddo

    ! first compute the strain at all the GLL points of the source element
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

          do l = 1,NGLLX
            hp1 = hprime_xx(i,l)
            tempx1l = tempx1l + displ_elem(1,l,j,k)*hp1
            tempy1l = tempy1l + displ_elem(2,l,j,k)*hp1
            tempz1l = tempz1l + displ_elem(3,l,j,k)*hp1

            hp2 = hprime_yy(j,l)
            tempx2l = tempx2l + displ_elem(1,i,l,k)*hp2
            tempy2l = tempy2l + displ_elem(2,i,l,k)*hp2
            tempz2l = tempz2l + displ_elem(3,i,l,k)*hp2

            hp3 = hprime_zz(k,l)
            tempx3l = tempx3l + displ_elem(1,i,j,l)*hp3
            tempy3l = tempy3l + displ_elem(2,i,j,l)*hp3
            tempz3l = tempz3l + displ_elem(3,i,j,l)*hp3
          enddo

          ! derivatives dudx,..
          xixl = xix_crust_mantle(i,j,k,ispec)
          xiyl = xiy_crust_mantle(i,j,k,ispec)
          xizl = xiz_crust_mantle(i,j,k,ispec)
          etaxl = etax_crust_mantle(i,j,k,ispec)
          etayl = etay_crust_mantle(i,j,k,ispec)
          etazl = etaz_crust_mantle(i,j,k,ispec)
          gammaxl = gammax_crust_mantle(i,j,k,ispec)
          gammayl = gammay_crust_mantle(i,j,k,ispec)
          gammazl = gammaz_crust_mantle(i,j,k,ispec)

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

          ! strain eps_jk
          ! symmetric strain definition: \epsilon = 1/2 ( grad(u) + grad(u)^T )
          eps_xx = duxdxl                                 ! dx/dx
          eps_xy = 0.5_CUSTOM_REAL * (duxdyl + duydxl)    ! dx/dy
          eps_xz = 0.5_CUSTOM_REAL * (duxdzl + duzdxl)    ! dx/dz
          eps_yy = duydyl                                 ! dy/dy
          eps_yz = 0.5_CUSTOM_REAL * (duydzl + duzdyl)    ! dy/dz
          eps_zz = duzdzl                                 ! dz/dz

          ! stores local element strain array
          eps_array(1,1,i,j,k) = eps_xx
          eps_array(2,1,i,j,k) = eps_xy      ! symmetry dx/dy = dy/dx
          eps_array(3,1,i,j,k) = eps_xz      ! symmetry dx/dz = dz/dx
          eps_array(1,2,i,j,k) = eps_xy
          eps_array(2,2,i,j,k) = eps_yy
          eps_array(3,2,i,j,k) = eps_yz
          eps_array(1,3,i,j,k) = eps_xz
          eps_array(2,3,i,j,k) = eps_yz
          eps_array(3,3,i,j,k) = eps_zz
        enddo
      enddo
    enddo

    hxir(:) = hxir_store(:,irec_local)
    hetar(:) = hetar_store(:,irec_local)
    hgammar(:) = hgammar_store(:,irec_local)

    ! interpolate the strain at actual receiver locations within the element from eps_array(:,:,i,j,k)
    eps_loc(:,:) = 0.0_CUSTOM_REAL
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          hlagrange = real(hxir(i)*hetar(j)*hgammar(k),kind=CUSTOM_REAL)
          eps_loc(1,1) = eps_loc(1,1) + eps_array(1,1,i,j,k)*hlagrange
          eps_loc(1,2) = eps_loc(1,2) + eps_array(1,2,i,j,k)*hlagrange
          eps_loc(1,3) = eps_loc(1,3) + eps_array(1,3,i,j,k)*hlagrange
          eps_loc(2,2) = eps_loc(2,2) + eps_array(2,2,i,j,k)*hlagrange
          eps_loc(2,3) = eps_loc(2,3) + eps_array(2,3,i,j,k)*hlagrange
          eps_loc(3,3) = eps_loc(3,3) + eps_array(3,3,i,j,k)*hlagrange
        enddo
      enddo
    enddo
    ! for completion purpose, symmetries
    eps_loc(2,1) = eps_loc(1,2)
    eps_loc(3,1) = eps_loc(1,3)
    eps_loc(3,2) = eps_loc(2,3)

    ! stores North, East and Vertical components
    !
    ! un-rotated
    !eps_loc_new(:,:) = eps_loc(:,:)
    !
    ! rotates from global x-y-z to the local coordinates (n-e-z):  eps_new = P*eps*P'
    ! nu is the rotation matrix from ECEF to local N-E-UP as defined.
    ! thus, if the nu is the rotation matrix that transforms coordinates from the global system (x,y,z) to the local
    ! coordinate system (N,E,V), e.g., a tensor is transformed as
    ! T_L = \nu * T_g * \nu^T
    !
    ! global -> local (n-e-up)
    ! eps_xx -> eps_nn
    ! eps_yy -> eps_ee
    ! eps_zz -> eps_zz (z in radial direction up)
    ! eps_xy -> eps_ne
    ! eps_xz -> eps_nz
    ! eps_yz -> eps_ez
    eps_loc_rot(:,:) = real(matmul(nu_rec(:,:,irec),eps_loc(:,:)),kind=CUSTOM_REAL)
    eps_loc_new(:,:) = real(matmul(eps_loc_rot(:,:), transpose(nu_rec(:,:,irec))),kind=CUSTOM_REAL)

    ! distinguish between single and double precision for reals
    !
    ! note: strain is dimensionless, no scaling to dimensionalize it needed (as for example for displacement seismograms)
    seismograms_eps(1,irec_local,seismo_current) = real(eps_loc_new(1,1), kind=CUSTOM_REAL)  ! \eps_nn
    seismograms_eps(2,irec_local,seismo_current) = real(eps_loc_new(2,2), kind=CUSTOM_REAL)  ! \eps_ee
    seismograms_eps(3,irec_local,seismo_current) = real(eps_loc_new(3,3), kind=CUSTOM_REAL)  ! \eps_zz
    seismograms_eps(4,irec_local,seismo_current) = real(eps_loc_new(1,2), kind=CUSTOM_REAL)  ! \eps_ne
    seismograms_eps(5,irec_local,seismo_current) = real(eps_loc_new(1,3), kind=CUSTOM_REAL)  ! \eps_nz
    seismograms_eps(6,irec_local,seismo_current) = real(eps_loc_new(2,3), kind=CUSTOM_REAL)  ! \eps_ez

  enddo ! irec_local

  end subroutine compute_seismograms_strain
