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

  subroutine get_model(iregion_code,ispec,nspec,idoubling, &
                      kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                      rhostore,dvpstore,nspec_ani, &
                      c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                      c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                      c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                      nspec_stacey,rho_vp,rho_vs, &
                      xstore,ystore,zstore, &
                      rmin,rmax, &
                      elem_in_crust,elem_in_mantle)

  use meshfem3D_par, only: &
    RCMB,RICB,R670,RMOHO,RTOPDDOUBLEPRIME,R600,R220, &
    R771,R400,R120,R80,RMIDDLE_CRUST,ROCEAN, &
    ABSORBING_CONDITIONS

  use meshfem3D_models_par

  use regions_mesh_par2, only: &
    Qmu_store,tau_e_store,tau_s,T_c_source

  implicit none

  integer :: iregion_code,ispec,nspec,idoubling

  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec) :: kappavstore,kappahstore, &
    muvstore,muhstore,eta_anisostore,rhostore,dvpstore

  integer :: nspec_ani
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store, &
    c33store,c34store,c35store,c36store, &
    c44store,c45store,c46store,c55store,c56store,c66store

  integer :: nspec_stacey
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ,nspec_stacey):: rho_vp,rho_vs

  double precision, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: xstore,ystore,zstore

  double precision :: rmin,rmax
  logical,intent(in) :: elem_in_crust,elem_in_mantle

  ! local parameters
  double precision :: xmesh,ymesh,zmesh
  ! the 21 coefficients for an anisotropic medium in reduced notation
  double precision :: c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33, &
                   c34,c35,c36,c44,c45,c46,c55,c56,c66

  double precision :: Qkappa,Qmu
  double precision, dimension(N_SLS) :: tau_e

  double precision :: rho,dvp
  double precision :: vpv,vph,vsv,vsh,eta_aniso
  double precision :: r,r_prem,moho
  integer :: i,j,k,i_sls

  ! it is *CRUCIAL* to leave this initialization here, this was the cause of the "s362ani + attenuation" bug in 2013 and 2014
  ! thus please never remove the line below
  moho = 0.d0

  ! loops over all GLL points for this spectral element
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX

        ! initializes values
        rho = 0.d0
        vpv = 0.d0
        vph = 0.d0
        vsv = 0.d0
        vsh = 0.d0

        eta_aniso = 1.d0 ! default for isotropic element

        c11 = 0.d0
        c12 = 0.d0
        c13 = 0.d0
        c14 = 0.d0
        c15 = 0.d0
        c16 = 0.d0
        c22 = 0.d0
        c23 = 0.d0
        c24 = 0.d0
        c25 = 0.d0
        c26 = 0.d0
        c33 = 0.d0
        c34 = 0.d0
        c35 = 0.d0
        c36 = 0.d0
        c44 = 0.d0
        c45 = 0.d0
        c46 = 0.d0
        c55 = 0.d0
        c56 = 0.d0
        c66 = 0.d0

        Qmu = 0.d0
        Qkappa = 0.d0 ! not used, not stored so far...
        tau_e(:) = 0.d0
        dvp = 0.d0

        ! sets xyz coordinates of GLL point
        xmesh = xstore(i,j,k,ispec)
        ymesh = ystore(i,j,k,ispec)
        zmesh = zstore(i,j,k,ispec)

        ! exact point location radius
        r = dsqrt(xmesh*xmesh + ymesh*ymesh + zmesh*zmesh)

        ! make sure we are within the right shell in PREM to honor discontinuities
        ! use small geometrical tolerance
        r_prem = r
        if (r <= rmin*1.000001d0) r_prem = rmin*1.000001d0
        if (r >= rmax*0.999999d0) r_prem = rmax*0.999999d0
        ! checks r_prem,rmin/rmax and assigned idoubling
        call get_model_check_idoubling(r_prem,xmesh,ymesh,zmesh,rmin,rmax,idoubling, &
                            RICB,RCMB,RTOPDDOUBLEPRIME, &
                            R220,R670)

        ! gets reference model values: rho,vpv,vph,vsv,vsh and eta_aniso
        call meshfem3D_models_get1D_val(iregion_code,idoubling, &
                              r_prem,rho,vpv,vph,vsv,vsh,eta_aniso, &
                              Qkappa,Qmu,RICB,RCMB, &
                              RTOPDDOUBLEPRIME,R80,R120,R220,R400,R600,R670,R771, &
                              RMOHO,RMIDDLE_CRUST,ROCEAN)

        ! gets the 3-D model parameters for the mantle
        call meshfem3D_models_get3Dmntl_val(iregion_code,r_prem,rho,dvp, &
                              vpv,vph,vsv,vsh,eta_aniso, &
                              RCMB,R670,RMOHO, &
                              xmesh,ymesh,zmesh,r, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26, &
                              c33,c34,c35,c36,c44,c45,c46,c55,c56,c66 &
#ifdef CEM
                              ,ispec,i,j,k &
#endif
                              )

        ! gets the 3-D crustal model
        ! M.A. don't overwrite crust if using CEM.
        if (CRUSTAL .and. .not. CEM_ACCEPT) then
          if (.not. elem_in_mantle) &
            call meshfem3D_models_get3Dcrust_val(iregion_code,xmesh,ymesh,zmesh,r, &
                              vpv,vph,vsv,vsh,rho,eta_aniso,dvp, &
                              c11,c12,c13,c14,c15,c16,c22,c23,c24,c25, &
                              c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66, &
                              elem_in_crust,moho)
        endif

        ! overwrites with tomographic model values (from iteration step) here, given at all GLL points
        call meshfem3D_models_impose_val(vpv,vph,vsv,vsh,rho,dvp,eta_aniso,iregion_code,ispec,i,j,k)

        ! checks vpv: if close to zero then there is probably an error
        if (vpv < TINYVAL) then
          print *,'Error vpv: ',vpv,' vph:',vph,' vsv: ',vsv,' vsh: ',vsh,' rho:',rho
          print *,'radius:',r*R_EARTH_KM
          call exit_mpi(myrank,'Error get_model values')
        endif

        !> Hejun, New attenuation assignment
        ! Define 3D and 1D attenuation after Moho stretch
        ! and before TOPOGRAPHY / ELLIPTICITY
        !
        !note:  only Qmu attenuation considered, Qkappa attenuation not used so far...
        if (ATTENUATION ) &
          call meshfem3D_models_getatten_val(idoubling,xmesh,ymesh,zmesh,r_prem, &
                                             tau_e,tau_s,T_c_source, &
                                             moho,Qmu,Qkappa,elem_in_crust)

! define elastic parameters in the model

        rhostore(i,j,k,ispec) = real(rho, kind=CUSTOM_REAL)
        kappavstore(i,j,k,ispec) = real(rho*(vpv*vpv - 4.d0*vsv*vsv/3.d0), kind=CUSTOM_REAL)
        kappahstore(i,j,k,ispec) = real(rho*(vph*vph - 4.d0*vsh*vsh/3.d0), kind=CUSTOM_REAL)
        muvstore(i,j,k,ispec) = real(rho*vsv*vsv, kind=CUSTOM_REAL)
        muhstore(i,j,k,ispec) = real(rho*vsh*vsh, kind=CUSTOM_REAL)
        eta_anisostore(i,j,k,ispec) = real(eta_aniso, kind=CUSTOM_REAL)

        if (HETEROGEN_3D_MANTLE) then
          dvpstore(i,j,k,ispec) = real(dvp, kind=CUSTOM_REAL)
        endif

        if (ABSORBING_CONDITIONS) then
          if (iregion_code == IREGION_OUTER_CORE) then
            ! we need just vp in the outer core for Stacey conditions
            rho_vp(i,j,k,ispec) = real(vph, kind=CUSTOM_REAL)
            rho_vs(i,j,k,ispec) = real(0.d0, kind=CUSTOM_REAL)
          else
            rho_vp(i,j,k,ispec) = real(rho*vph, kind=CUSTOM_REAL)
            rho_vs(i,j,k,ispec) = real(rho*vsh, kind=CUSTOM_REAL)
          endif
        endif

        if (ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
          c11store(i,j,k,ispec) = real(c11, kind=CUSTOM_REAL)
          c33store(i,j,k,ispec) = real(c33, kind=CUSTOM_REAL)
          c12store(i,j,k,ispec) = real(c12, kind=CUSTOM_REAL)
          c13store(i,j,k,ispec) = real(c13, kind=CUSTOM_REAL)
          c44store(i,j,k,ispec) = real(c44, kind=CUSTOM_REAL)
        endif

        if (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
          c11store(i,j,k,ispec) = real(c11, kind=CUSTOM_REAL)
          c12store(i,j,k,ispec) = real(c12, kind=CUSTOM_REAL)
          c13store(i,j,k,ispec) = real(c13, kind=CUSTOM_REAL)
          c14store(i,j,k,ispec) = real(c14, kind=CUSTOM_REAL)
          c15store(i,j,k,ispec) = real(c15, kind=CUSTOM_REAL)
          c16store(i,j,k,ispec) = real(c16, kind=CUSTOM_REAL)
          c22store(i,j,k,ispec) = real(c22, kind=CUSTOM_REAL)
          c23store(i,j,k,ispec) = real(c23, kind=CUSTOM_REAL)
          c24store(i,j,k,ispec) = real(c24, kind=CUSTOM_REAL)
          c25store(i,j,k,ispec) = real(c25, kind=CUSTOM_REAL)
          c26store(i,j,k,ispec) = real(c26, kind=CUSTOM_REAL)
          c33store(i,j,k,ispec) = real(c33, kind=CUSTOM_REAL)
          c34store(i,j,k,ispec) = real(c34, kind=CUSTOM_REAL)
          c35store(i,j,k,ispec) = real(c35, kind=CUSTOM_REAL)
          c36store(i,j,k,ispec) = real(c36, kind=CUSTOM_REAL)
          c44store(i,j,k,ispec) = real(c44, kind=CUSTOM_REAL)
          c45store(i,j,k,ispec) = real(c45, kind=CUSTOM_REAL)
          c46store(i,j,k,ispec) = real(c46, kind=CUSTOM_REAL)
          c55store(i,j,k,ispec) = real(c55, kind=CUSTOM_REAL)
          c56store(i,j,k,ispec) = real(c56, kind=CUSTOM_REAL)
          c66store(i,j,k,ispec) = real(c66, kind=CUSTOM_REAL)
        endif

        if (ATTENUATION) then
          if (ATTENUATION_3D .or. ATTENUATION_1D_WITH_3D_STORAGE) then

            ! distinguish between single and double precision for reals
            do i_sls = 1,N_SLS
              tau_e_store(i,j,k,i_sls,ispec) = real(tau_e(i_sls), kind=CUSTOM_REAL)
            enddo
            Qmu_store(i,j,k,ispec) = real(Qmu, kind=CUSTOM_REAL)

          else

            ! distinguish between single and double precision for reals
            ! store values from mid-point for whole element
            if (i == NGLLX/2 .and. j == NGLLY/2 .and. k == NGLLZ/2) then
              do i_sls = 1,N_SLS
                tau_e_store(1,1,1,i_sls,ispec) = real(tau_e(i_sls), kind=CUSTOM_REAL)
              enddo
              Qmu_store(1,1,1,ispec) = real(Qmu, kind=CUSTOM_REAL)
            endif

          endif
        endif

      enddo
    enddo
  enddo

  end subroutine get_model


!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_model_check_idoubling(r_prem,x,y,z,rmin,rmax,idoubling, &
                            RICB,RCMB,RTOPDDOUBLEPRIME, &
                            R220,R670)

  use meshfem3D_models_par

  implicit none

  integer :: idoubling

  double precision :: r_prem,rmin,rmax,x,y,z

  double precision :: RICB,RCMB,RTOPDDOUBLEPRIME,R670,R220
  double precision :: r_m,r,theta,phi

  ! compute real physical radius in meters
  r_m = r_prem * R_EARTH

  ! checks layers
  if (abs(rmax - rmin ) < TINYVAL) then
    ! there's probably an error
    print *,'Error layer radius min/max:',rmin,rmax
    print *,'  point radius: ',r_prem
    call exit_mpi(myrank,'Error  in get_model_check_idoubling() layer radius')
  endif


  ! check flags to make sure we correctly honor the discontinuities
  ! we use strict inequalities since r has been slightly changed in mesher

  !
  !--- inner core
  !
  if (r_m >= 0.d0 .and. r_m < RICB) then
    if (idoubling /= IFLAG_INNER_CORE_NORMAL .and. &
       idoubling /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_TOP_CENTRAL_CUBE .and. &
       idoubling /= IFLAG_IN_FICTITIOUS_CUBE) then
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      print *,'Error point r/lat/lon:',r_m,90.0 - theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS
      print *,'  idoubling/IFLAG: ',idoubling,IFLAG_INNER_CORE_NORMAL,'-to-',IFLAG_IN_FICTITIOUS_CUBE
      call exit_MPI(myrank,'Error  in get_model_check_idoubling() wrong doubling flag for inner core point')
    endif
  !
  !--- outer core
  !
  else if (r_m > RICB .and. r_m < RCMB) then
    if (idoubling /= IFLAG_OUTER_CORE_NORMAL) then
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      print *,'Error point r/lat/lon:',r_m,90.0 - theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS
      print *,'  idoubling/IFLAG: ',idoubling,IFLAG_OUTER_CORE_NORMAL
      call exit_MPI(myrank,'Error  in get_model_check_idoubling() wrong doubling flag for outer core point')
    endif
  !
  !--- D" at the base of the mantle
  !
  else if (r_m > RCMB .and. r_m < RTOPDDOUBLEPRIME) then
    if (idoubling /= IFLAG_MANTLE_NORMAL) then
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      print *,'Error point r/lat/lon:',r_m,90.0 - theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS
      print *,'  dprime radius/RCMB/RTOPDDOUBLEPRIME:',r_m, RCMB,RTOPDDOUBLEPRIME
      print *,'  idoubling/IFLAG: ',idoubling,IFLAG_MANTLE_NORMAL
      call exit_MPI(myrank,'Error  in get_model_check_idoubling() wrong doubling flag for D" point')
    endif
  !
  !--- mantle: from top of D" to d670
  !
  else if (r_m > RTOPDDOUBLEPRIME .and. r_m < R670) then
    if (idoubling /= IFLAG_MANTLE_NORMAL) then
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      print *,'Error point r/lat/lon:',r_m,90.0 - theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS
      print *,'  idoubling/IFLAG: ',idoubling,IFLAG_MANTLE_NORMAL
      call exit_MPI(myrank,'Error  in get_model_check_idoubling() wrong doubling flag for top D" to d670 point')
    endif

  !
  !--- mantle: from d670 to d220
  !
  else if (r_m > R670 .and. r_m < R220) then
    if (idoubling /= IFLAG_670_220) then
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      print *,'Error point r/lat/lon:',r_m,90.0 - theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS
      print *,'  idoubling/IFLAG: ',idoubling,IFLAG_670_220
      call exit_MPI(myrank,'Error  in get_model_check_idoubling() wrong doubling flag for d670 to d220 point')
    endif

  !
  !--- mantle and crust: from d220 to MOHO and then to surface
  !
  else if (r_m > R220) then
    if (idoubling /= IFLAG_220_80 .and. idoubling /= IFLAG_80_MOHO .and. idoubling /= IFLAG_CRUST) then
      call xyz_2_rthetaphi_dble(x,y,z,r,theta,phi)
      print *,'Error point r/lat/lon:',r_m,90.0 - theta/DEGREES_TO_RADIANS,phi/DEGREES_TO_RADIANS
      print *,'  idoubling/IFLAG: ',idoubling,IFLAG_220_80,IFLAG_80_MOHO,IFLAG_CRUST
      call exit_MPI(myrank,'Error  in get_model_check_idoubling() wrong doubling flag for d220 to Moho to surface point')
    endif

  endif

  end subroutine get_model_check_idoubling
