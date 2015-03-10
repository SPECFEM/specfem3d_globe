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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"

  subroutine compute_seismograms(nglob,displ,seismo_current,seismograms)

  use constants_solver

  use specfem_par,only: &
    NTSTEP_BETWEEN_OUTPUT_SEISMOS, &
    nrec_local,nu,ispec_selected_rec,number_receiver_global, &
    scale_displ,hlagrange_store

  use specfem_par_crustmantle,only: ibool_crust_mantle

  implicit none

  integer,intent(in) :: nglob
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob),intent(in) :: displ

  integer,intent(in) :: seismo_current

  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(out) :: &
    seismograms

  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  integer :: ispec,iglob,irec_local,irec
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  do irec_local = 1,nrec_local

    ! get global number of that receiver
    irec = number_receiver_global(irec_local)

    ispec = ispec_selected_rec(irec)

    ! perform the general interpolation using Lagrange polynomials
    uxd = ZERO
    uyd = ZERO
    uzd = ZERO

    DO_LOOP_IJK

      iglob = ibool_crust_mantle(INDEX_IJK,ispec)

      hlagrange = hlagrange_store(INDEX_IJK,irec_local)

      uxd = uxd + dble(displ(1,iglob))*hlagrange
      uyd = uyd + dble(displ(2,iglob))*hlagrange
      uzd = uzd + dble(displ(3,iglob))*hlagrange

    ENDDO_LOOP_IJK

    ! store North, East and Vertical components
    ! distinguish between single and double precision for reals
    seismograms(:,irec_local,seismo_current) = real(scale_displ*(nu(:,1,irec)*uxd + &
                                                                 nu(:,2,irec)*uyd + &
                                                                 nu(:,3,irec)*uzd), &
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
                                         nit_written, &
                                         moment_der,sloc_der,stshift_der,shdur_der, &
                                         seismograms)

  use constants_solver,only: &
    CUSTOM_REAL,SIZE_REAL,ZERO,ONE,PI,GRAV,RHOAV,NGLLX,NGLLY,NGLLZ, NGLLCUBE, &
    NDIM,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
    NSPEC_CRUST_MANTLE_STRAIN_ONLY,NSPEC_CRUST_MANTLE_STR_OR_ATT

  use specfem_par,only: &
    NSTEP,NTSTEP_BETWEEN_OUTPUT_SEISMOS,UNDO_ATTENUATION, &
    nrec_local, &
    nu_source,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
    hlagrange_store, &
    hxir_store,hpxir_store,hetar_store,hpetar_store,hgammar_store,hpgammar_store, &
    tshift_cmt,hdur_gaussian, &
    DT,t0,deltat,it, &
    scale_displ, &
    hprime_xx,hprime_yy,hprime_zz, &
    ispec_selected_source,number_receiver_global

  use specfem_par_crustmantle,only: ibool_crust_mantle, &
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

  integer,intent(in) :: nit_written

  real(kind=CUSTOM_REAL), dimension(NDIM,NDIM,nrec_local),intent(inout) :: moment_der
  real(kind=CUSTOM_REAL), dimension(NDIM,nrec_local),intent(inout) :: sloc_der
  real(kind=CUSTOM_REAL), dimension(nrec_local),intent(inout) :: stshift_der, shdur_der

  real(kind=CUSTOM_REAL), dimension(NDIM*NDIM,nrec_local,NTSTEP_BETWEEN_OUTPUT_SEISMOS),intent(out) :: &
    seismograms

  ! local parameters
  double precision :: uxd,uyd,uzd,hlagrange
  double precision :: eps_trace,dxx,dyy,dxy,dxz,dyz
  double precision :: eps_loc(NDIM,NDIM), eps_loc_new(NDIM,NDIM)
  double precision :: stf
  double precision :: scale_t

  real(kind=CUSTOM_REAL) :: displ_s(NDIM,NGLLX,NGLLY,NGLLZ)
  real(kind=CUSTOM_REAL) :: eps_s(NDIM,NDIM), eps_m_s, &
        eps_m_l_s(NDIM), stf_deltat, Kp_deltat, Hp_deltat
  integer :: iglob,irec_local,irec,ispec

  double precision, external :: comp_source_time_function

  ! element strain
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: eps_trace_over_3_loc_cm
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,5) :: epsilondev_loc_crust_mantle

  double precision :: hxir(NGLLX), hetar(NGLLY), hgammar(NGLLZ), &
                      hpxir(NGLLX),hpetar(NGLLY),hpgammar(NGLLZ)
#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  do irec_local = 1,nrec_local

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
                                                epsilondev_loc_crust_mantle,eps_trace_over_3_loc_cm)
    else
      ! element adjoint strain
      eps_trace_over_3_loc_cm(:,:,:) = eps_trace_over_3_crust_mantle(:,:,:,ispec)
      epsilondev_loc_crust_mantle(:,:,:,1) = epsilondev_xx_crust_mantle(:,:,:,ispec)
      epsilondev_loc_crust_mantle(:,:,:,2) = epsilondev_yy_crust_mantle(:,:,:,ispec)
      epsilondev_loc_crust_mantle(:,:,:,3) = epsilondev_xy_crust_mantle(:,:,:,ispec)
      epsilondev_loc_crust_mantle(:,:,:,4) = epsilondev_xz_crust_mantle(:,:,:,ispec)
      epsilondev_loc_crust_mantle(:,:,:,5) = epsilondev_yz_crust_mantle(:,:,:,ispec)
    endif

    ! perform the general interpolation using Lagrange polynomials
    DO_LOOP_IJK

      iglob = ibool_crust_mantle(INDEX_IJK,ispec)

      hlagrange = hlagrange_store(INDEX_IJK,irec_local)

      uxd = uxd + dble(displ_crust_mantle(1,iglob))*hlagrange
      uyd = uyd + dble(displ_crust_mantle(2,iglob))*hlagrange
      uzd = uzd + dble(displ_crust_mantle(3,iglob))*hlagrange

      eps_trace = eps_trace + dble(eps_trace_over_3_loc_cm(INDEX_IJK))*hlagrange

      dxx = dxx + dble(epsilondev_loc_crust_mantle(INDEX_IJK,1))*hlagrange
      dyy = dyy + dble(epsilondev_loc_crust_mantle(INDEX_IJK,2))*hlagrange
      dxy = dxy + dble(epsilondev_loc_crust_mantle(INDEX_IJK,3))*hlagrange
      dxz = dxz + dble(epsilondev_loc_crust_mantle(INDEX_IJK,4))*hlagrange
      dyz = dyz + dble(epsilondev_loc_crust_mantle(INDEX_IJK,5))*hlagrange

      displ_s(:,INDEX_IJK) = displ_crust_mantle(:,iglob)

    ENDDO_LOOP_IJK

    eps_loc(1,1) = eps_trace + dxx
    eps_loc(2,2) = eps_trace + dyy
    eps_loc(3,3) = eps_trace - dxx - dyy
    eps_loc(1,2) = dxy
    eps_loc(1,3) = dxz
    eps_loc(2,3) = dyz
    eps_loc(2,1) = dxy
    eps_loc(3,1) = dxz
    eps_loc(3,2) = dyz

    eps_loc_new(:,:) = eps_loc(:,:)
    ! rotate to the local Cartesian coordinates (n-e-z):  eps_new=P*eps*P'
    eps_loc_new(:,:) = matmul(matmul(nu_source(:,:,irec),eps_loc(:,:)), transpose(nu_source(:,:,irec)))

    ! distinguish between single and double precision for reals
    seismograms(1,irec_local,it-nit_written) = real(eps_loc_new(1,1), kind=CUSTOM_REAL)
    seismograms(2,irec_local,it-nit_written) = real(eps_loc_new(2,2), kind=CUSTOM_REAL)
    seismograms(3,irec_local,it-nit_written) = real(eps_loc_new(3,3), kind=CUSTOM_REAL)
    seismograms(4,irec_local,it-nit_written) = real(eps_loc_new(1,2), kind=CUSTOM_REAL)
    seismograms(5,irec_local,it-nit_written) = real(eps_loc_new(1,3), kind=CUSTOM_REAL)
    seismograms(6,irec_local,it-nit_written) = real(eps_loc_new(2,3), kind=CUSTOM_REAL)
    seismograms(7:9,irec_local,it-nit_written) = real(scale_displ*(nu_source(:,1,irec)*uxd + &
                                                                   nu_source(:,2,irec)*uyd + &
                                                                   nu_source(:,3,irec)*uzd), &
                                                      kind=CUSTOM_REAL)

    ! interpolators
    ! note: we explicitly copy the store arrays to local temporary arrays here
    !       the array indexing (irec_local,:) is non-contiguous and compilers would have to do this anyway
    hxir(:) = hxir_store(irec_local,:)
    hetar(:) = hetar_store(irec_local,:)
    hgammar(:) = hgammar_store(irec_local,:)
    hpxir(:) = hpxir_store(irec_local,:)
    hpetar(:) = hpetar_store(irec_local,:)
    hpgammar(:) = hpgammar_store(irec_local,:)

    ! Frechet derivatives of the source
    call compute_adj_source_frechet(displ_s,Mxx(irec),Myy(irec),Mzz(irec), &
                Mxy(irec),Mxz(irec),Myz(irec),eps_s,eps_m_s,eps_m_l_s, &
                hxir,hetar,hgammar, &
                hpxir,hpetar,hpgammar, &
                hprime_xx,hprime_yy,hprime_zz, &
                xix_crust_mantle(1,1,1,ispec),xiy_crust_mantle(1,1,1,ispec),xiz_crust_mantle(1,1,1,ispec), &
                etax_crust_mantle(1,1,1,ispec),etay_crust_mantle(1,1,1,ispec),etaz_crust_mantle(1,1,1,ispec), &
                gammax_crust_mantle(1,1,1,ispec),gammay_crust_mantle(1,1,1,ispec),gammaz_crust_mantle(1,1,1,ispec))

    stf = comp_source_time_function(dble(NSTEP-it)*DT-t0-tshift_cmt(irec),hdur_gaussian(irec))
    stf_deltat = stf * deltat

    moment_der(:,:,irec_local) = moment_der(:,:,irec_local) + eps_s(:,:) * stf_deltat
    sloc_der(:,irec_local) = sloc_der(:,irec_local) + eps_m_l_s(:) * stf_deltat

    scale_t = ONE/dsqrt(PI*GRAV*RHOAV)

    Kp_deltat= -1.0d0/sqrt(PI)/hdur_gaussian(irec)*exp(-((dble(NSTEP-it)*DT-t0-tshift_cmt(irec))/hdur_gaussian(irec))**2) &
                       * deltat * scale_t
    Hp_deltat= (dble(NSTEP-it)*DT-t0-tshift_cmt(irec))/hdur_gaussian(irec)*Kp_deltat

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
!  use specfem_par,only: CUSTOM_REAL,NDIM, &
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

