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

!--------------------------------------------------------------------------------------------------
! Generic transversely isotropic perturbation models based on spherical harmonics expansion:
! SH lmax = 20 for mantle and 21 splines vertically,
!    lmax = 40 for crust with one layer
!
! PREM as a background model so far, 3D ATTENUATION is not yet incorporated.
!
! Jeannot Trampert, 2015
!--------------------------------------------------------------------------------------------------


!-----------------------------------------------------------------------------------------
!
! Crustal model
!
!-----------------------------------------------------------------------------------------

  module model_full_sh_crust_par

  ! three_d_mantle_model_constants
  integer, parameter :: NS_40 = 40
  integer, parameter :: NSH_40 = (NS_40+1)**2

  ! model_crust_variables
  double precision,dimension(:),allocatable :: CRUST_SH_V_rho, &
    CRUST_SH_V_vsh,CRUST_SH_V_vsv,CRUST_SH_V_vph,CRUST_SH_V_vpv,CRUST_SH_V_eta

  double precision,dimension(:),allocatable :: CRUST_SH_V_moho

  end module model_full_sh_crust_par

!
!-----------------------------------------------------------------------------------------
!

  subroutine model_crust_sh_broadcast()

! standard routine to setup model

  use constants, only: IMAIN,myrank

  use model_full_sh_crust_par

  implicit none

  ! local parameters
  integer :: ier

  ! allocate memory
  allocate (CRUST_SH_V_vsh(NSH_40), &   ! Full TI for one layer crust
            CRUST_SH_V_vsv(NSH_40), &
            CRUST_SH_V_vph(NSH_40), &
            CRUST_SH_V_vpv(NSH_40), &
            CRUST_SH_V_eta(NSH_40), &
            CRUST_SH_V_rho(NSH_40), &
            CRUST_SH_V_moho(NSH_40), &
            stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating crust arrays')
  CRUST_SH_V_vph(:) = 0.d0; CRUST_SH_V_vpv(:) = 0.d0
  CRUST_SH_V_vsh(:) = 0.d0; CRUST_SH_V_vsv(:) = 0.d0
  CRUST_SH_V_rho(:) = 0.d0; CRUST_SH_V_eta(:) = 0.d0
  CRUST_SH_V_moho(:) = 0.d0

  ! the variables read are declared and stored in structure CRUST_SH_V
  if (myrank == 0) call read_crust_sh_model()

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(CRUST_SH_V_vsh,NSH_40)
  call bcast_all_dp(CRUST_SH_V_vsv,NSH_40)
  call bcast_all_dp(CRUST_SH_V_vph,NSH_40)
  call bcast_all_dp(CRUST_SH_V_vpv,NSH_40)
  call bcast_all_dp(CRUST_SH_V_eta,NSH_40)
  call bcast_all_dp(CRUST_SH_V_rho,NSH_40)
  call bcast_all_dp(CRUST_SH_V_moho,NSH_40)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done full_sh crust'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine model_crust_sh_broadcast

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_crust_sh_model()

  use constants, only: ZERO,IIN,IMAIN,MAX_STRING_LEN,HUGEVAL

  use model_full_sh_crust_par

  implicit none

  ! local parameters
  integer :: l,lmax,nm,ipar,ier
  character(len=32) :: modelname
  character(len=MAX_STRING_LEN) :: rootdir
  ! for moho statistics
  integer :: i,j
  real(kind=4) :: shcof(NSH_40)
  real(kind=4) :: xlat,xlon
  double precision :: moho,moho_min,moho_max

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating crustal model: crust_sh '
  write(IMAIN,*)
  call flush_IMAIN()

  ! initialize
  do l = 1,NSH_40
    CRUST_SH_V_vsh(l) = ZERO
    CRUST_SH_V_vsv(l) = ZERO
    CRUST_SH_V_vph(l) = ZERO
    CRUST_SH_V_vpv(l) = ZERO
    CRUST_SH_V_eta(l) = ZERO
    CRUST_SH_V_rho(l) = ZERO
    CRUST_SH_V_moho(l) = ZERO
  enddo

  ! user output
  write(IMAIN,*) '  reading crustal spherical harmonic model from full_sh:'
  call flush_IMAIN()

  ! root directory
  rootdir = 'DATA/full_sphericalharmonic_model/'

  do ipar = 1,7

    select case(ipar)
    case (1)
      ! VSV degree 40 max model 1 layer
      modelname = 'CVSV.dat'
    case (2)
      ! VSH
      modelname = 'CVSH.dat'
    case (3)
      ! VPV
      modelname = 'CVPV.dat'
    case (4)
      ! VPH
      modelname = 'CVPH.dat'
    case (5)
      ! ETA
      modelname = 'CETA.dat'
    case (6)
      ! RHO
      modelname = 'CRHO.dat'
    case (7)
      ! MOHO degree 40 model
      modelname = 'MOHO.dat'
    case default
      stop 'Invalid parameter for crust file in full_sh'
    end select

    ! opens model file
    open(unit=IIN,file=trim(rootdir)//trim(modelname),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(rootdir)//trim(modelname)
      call exit_MPI(0,'Error opening file full_sphericalharmonic model file')
    endif

    read(IIN,*) lmax

    ! check maximum degree
    if (lmax > NS_40) then
      print *,'Error full_sh model file: lmax too big ',lmax,' for file ',trim(rootdir)//trim(modelname)
      call exit_MPI(0,'lmax too big in full_sh model file')
    endif

    ! user output
    write(IMAIN,*) '  name ',trim(modelname),': lmax = ',lmax
    call flush_IMAIN()

    ! reads in spherical harmonics coefficients
    nm = (lmax+1)**2

    select case(ipar)
    case (1)
      read(IIN,*) (CRUST_SH_V_vsv(l),l=1,nm)
    case (2)
      read(IIN,*) (CRUST_SH_V_vsh(l),l=1,nm)
    case (3)
      read(IIN,*) (CRUST_SH_V_vpv(l),l=1,nm)
    case (4)
      read(IIN,*) (CRUST_SH_V_vph(l),l=1,nm)
    case (5)
      read(IIN,*) (CRUST_SH_V_eta(l),l=1,nm)
    case (6)
      read(IIN,*) (CRUST_SH_V_rho(l),l=1,nm)
    case (7)
      read(IIN,*) (CRUST_SH_V_moho(l),l=1,nm)
    end select

    close(IIN)

  enddo

  ! moho thickness statistics
  ! gets min/max value for a set of sample locations (regular lat/lon grid points)
  moho_min = HUGEVAL
  moho_max = -HUGEVAL
  do i = -90,90
    do j = -180,180
      ! lat/lon in degrees
      xlat = i * 1.0
      xlon = j * 1.0
      call ylm(xlat,xlon,NS_40,shcof)

      ! gets moho depth (in km)
      moho = ZERO
      do l = 1,NSH_40
        moho = moho + CRUST_SH_V_moho(l)*dble(shcof(l))
      enddo

      ! checks
      if (moho <= ZERO) then
        print *,'Error in crustal sh model: invalid moho depth ',moho,' at lat/lon = ',xlat,xlon
        stop 'Error invalid moho depth in full_sh model'
      endif

      ! statistics
      if (moho > moho_max) moho_max = moho
      if (moho < moho_min) moho_min = moho
    enddo
  enddo

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) '  Moho crustal thickness min/max = ',sngl(moho_min),sngl(moho_max),' km'
  call flush_IMAIN()

  end subroutine read_crust_sh_model

!
!-----------------------------------------------------------------------------------------
!

  subroutine crust_sh(lat,lon,r,vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment,found_crust,elem_in_crust,moho_only)

! gets crustal value for location lat/lon/r

  use constants, only: ZERO,PI,GRAV
  use shared_parameters, only: R_PLANET,RHOAV

  use model_full_sh_crust_par

  implicit none

  ! lat/lon  - in degrees (range lat/lon = [-90,90] / [-180,180]
  ! radius r - normalized by globe radius [0,1.x]
  double precision,intent(in) :: lat,lon,r

  double precision,intent(out) :: vpvc,vphc,vsvc,vshc,etac,rhoc,moho,sediment
  logical,intent(out) :: found_crust
  logical,intent(in) :: elem_in_crust,moho_only

  ! local parameters
  double precision :: depth,scaleval
  integer :: l
  real(kind=4) :: shcof(NSH_40)
  real(kind=4) :: xlat,xlon

  ! initializes
  vsvc = ZERO !km/s
  vshc = ZERO !km/s
  vpvc = ZERO !km/s
  vphc = ZERO !km/s
  rhoc = ZERO !g/cm^3
  moho = ZERO !km
  sediment = ZERO
  etac = ZERO
  found_crust = .true.

  ! gets coefficient for location
  xlat = sngl(lat)
  xlon = sngl(lon)
  call ylm(xlat,xlon,NS_40,shcof)

  ! calculates model values
  do l = 1,NSH_40
    vsvc = vsvc + CRUST_SH_V_vsv(l)*dble(shcof(l))
    vshc = vshc + CRUST_SH_V_vsh(l)*dble(shcof(l))
    vpvc = vpvc + CRUST_SH_V_vpv(l)*dble(shcof(l))
    vphc = vphc + CRUST_SH_V_vph(l)*dble(shcof(l))
    etac = etac + CRUST_SH_V_eta(l)*dble(shcof(l))
    rhoc = rhoc + CRUST_SH_V_rho(l)*dble(shcof(l))
    moho = moho + CRUST_SH_V_moho(l)*dble(shcof(l))
  enddo

  ! scales (non-dimensionalizes) values
  moho = moho * 1000.d0/R_PLANET

  ! checks if anything further to do
  if (moho_only) return

  ! no sediment thickness information
  sediment = 0.d0

  ! check values
  if (vsvc <= ZERO) found_crust = .false.
  if (vshc <= ZERO) found_crust = .false.
  if (vpvc <= ZERO) found_crust = .false.
  if (vphc <= ZERO) found_crust = .false.
  if (etac < 0.6d0) found_crust = .false.
  if (etac > 1.4d0) found_crust = .false.
  if (rhoc <= ZERO) found_crust = .false.
  if (moho <= ZERO) found_crust = .false.

  ! checks if crustal value found
  if (.not. found_crust) then
    print *,'Error: crustal value not found in full_sh model! position is at lat/lon/r = ',lat,lon,r
    stop 'crustal value not found in full_sh model, please check...'
  endif

  ! checks if position is in crust
  found_crust = .false.
  if (elem_in_crust) then
    ! forces points to have crustal values
    found_crust = .true.
  else
    ! checks if depth of position above/below moho
    depth = (1.d0 - r)  ! non-dimensional
    !debug
    !print *,'sh crust: depth = ',depth,' moho = ',moho
    if (depth <= moho) then
      found_crust = .true.
    endif
  endif

  if (found_crust) then
    scaleval = dsqrt(PI*GRAV*RHOAV)
    vsvc = vsvc * 1000.d0/(R_PLANET*scaleval)
    vshc = vshc * 1000.d0/(R_PLANET*scaleval)
    vpvc = vpvc * 1000.d0/(R_PLANET*scaleval)
    vphc = vphc * 1000.d0/(R_PLANET*scaleval)
    rhoc = rhoc * 1000.d0/RHOAV
  endif

  end subroutine crust_sh



!-----------------------------------------------------------------------------------------
!
! Mantle model
!
!-----------------------------------------------------------------------------------------


  module model_full_sh_mantle_par

  ! three_d_mantle_model_constants
  integer, parameter :: NK_20 = 20
  integer, parameter :: NS_20 = 20
  integer, parameter :: NSH_20 = (NS_20+1)**2

  ! model_mantle_variables
  double precision,dimension(:,:),allocatable :: &
    MANTLE_SH_V_dvsh,MANTLE_SH_V_dvsv,MANTLE_SH_V_dvph,MANTLE_SH_V_dvpv,MANTLE_SH_V_deta,MANTLE_SH_V_drho

  double precision,dimension(:),allocatable :: &
    MANTLE_SH_V_t410,MANTLE_SH_V_t660,MANTLE_SH_V_tcmb

  ! splines
  double precision,dimension(:),allocatable :: MANTLE_SH_V_spknt
  double precision,dimension(:,:),allocatable :: MANTLE_SH_V_qq0
  double precision,dimension(:,:,:),allocatable :: MANTLE_SH_V_qq

  end module model_full_sh_mantle_par


!
!-----------------------------------------------------------------------------------------
!

  subroutine model_mantle_sh_broadcast()

! standard routine to setup model

  use constants, only: IMAIN,myrank

  use model_full_sh_mantle_par

  implicit none

  integer :: ier

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) 'broadcast model: full_sh'
    call flush_IMAIN()
  endif

  ! model_mantle_sh_variables
  allocate(MANTLE_SH_V_dvsh(0:NK_20,1:NSH_20), &   ! Full TI
           MANTLE_SH_V_dvsv(0:NK_20,1:NSH_20), &
           MANTLE_SH_V_dvph(0:NK_20,1:NSH_20), &
           MANTLE_SH_V_dvpv(0:NK_20,1:NSH_20), &
           MANTLE_SH_V_deta(0:NK_20,1:NSH_20), &
           MANTLE_SH_V_drho(0:NK_20,1:NSH_20), &
           MANTLE_SH_V_t410(1:NSH_20), &
           MANTLE_SH_V_t660(1:NSH_20), &
           MANTLE_SH_V_tcmb(1:NSH_20), &
           MANTLE_SH_V_spknt(NK_20+1), &
           MANTLE_SH_V_qq0(NK_20+1,NK_20+1), &
           MANTLE_SH_V_qq(3,NK_20+1,NK_20+1), &
           stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'Error allocating MANTLE_SH_V arrays')

  ! the variables read are declared and stored in structure S20RTS_V
  if (myrank == 0) call read_model_mantle_sh()

  ! broadcast the information read on the main node to all the nodes
  call bcast_all_dp(MANTLE_SH_V_dvsh,(NK_20+1)*(NSH_20))
  call bcast_all_dp(MANTLE_SH_V_dvsv,(NK_20+1)*(NSH_20))
  call bcast_all_dp(MANTLE_SH_V_dvph,(NK_20+1)*(NSH_20))
  call bcast_all_dp(MANTLE_SH_V_dvpv,(NK_20+1)*(NSH_20))
  call bcast_all_dp(MANTLE_SH_V_deta,(NK_20+1)*(NSH_20))
  call bcast_all_dp(MANTLE_SH_V_drho,(NK_20+1)*(NSH_20))

  call bcast_all_dp(MANTLE_SH_V_t410,NSH_20)
  call bcast_all_dp(MANTLE_SH_V_t660,NSH_20)
  call bcast_all_dp(MANTLE_SH_V_tcmb,NSH_20)

  call bcast_all_dp(MANTLE_SH_V_spknt,NK_20+1)
  call bcast_all_dp(MANTLE_SH_V_qq0,(NK_20+1)*(NK_20+1))
  call bcast_all_dp(MANTLE_SH_V_qq,3*(NK_20+1)*(NK_20+1))

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  done full_sh mantle'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine model_mantle_sh_broadcast


!
!-----------------------------------------------------------------------------------------
!

  subroutine read_model_mantle_sh()

  use constants, only: ZERO,IIN,IMAIN,MAX_STRING_LEN,HUGEVAL

  use model_full_sh_mantle_par

  implicit none

  ! local parameters
  integer :: k,l,kmax,lmax,nm,ipar,ier
  character(len=32) :: modelname
  character(len=MAX_STRING_LEN) :: rootdir
  ! for statistics
  integer :: i,j
  real(kind=4) :: shcof(NSH_20)
  real(kind=4) :: xlat,xlon
  double precision :: cmb,cmb_min,cmb_max
  double precision :: t410,t410_min,t410_max
  double precision :: t660,t660_min,t660_max

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) 'incorporating mantle model: full_sh '
  write(IMAIN,*)
  call flush_IMAIN()

  ! initialize
  do k = 0,NK_20
    do l = 1,NSH_20
      MANTLE_SH_V_dvsh(k,l) = ZERO
      MANTLE_SH_V_dvsv(k,l) = ZERO
      MANTLE_SH_V_dvph(k,l) = ZERO
      MANTLE_SH_V_dvpv(k,l) = ZERO
      MANTLE_SH_V_deta(k,l) = ZERO
      MANTLE_SH_V_drho(k,l) = ZERO
    enddo
  enddo
  do l = 1,NSH_20
    MANTLE_SH_V_t410(l) = ZERO
    MANTLE_SH_V_t660(l) = ZERO
    MANTLE_SH_V_tcmb(l) = ZERO
  enddo

  ! user output
  write(IMAIN,*) '  reading mantle spherical harmonic model:'
  call flush_IMAIN()

  ! root directory
  rootdir = 'DATA/full_sphericalharmonic_model/'

  ! mantle tiso values
  do ipar = 1,6

    select case(ipar)
    case (1)
      ! VSV degree 20 dSV model sph file structure
      modelname = 'MVSV.dat'
    case (2)
      ! VSH degree 20 dSH model sph file structure
      modelname = 'MVSH.dat'
    case (3)
      ! VPV degree 20 dPV model sph file structure
      modelname = 'MVPV.dat'
    case (4)
      ! VPH degree 20 dPH model sph file structure
      modelname = 'MVPH.dat'
    case (5)
      ! ETA degree 20 dETA model sph file structure
      modelname = 'META.dat'
    case (6)
      ! RHO degree 20 dRHO model sph file structure
      modelname = 'MRHO.dat'
    case default
      stop 'Invalid parameter for mantle file in full_sh'
    end select


    ! opens model file
    open(unit=IIN,file=trim(rootdir)//trim(modelname),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(rootdir)//trim(modelname)
      call exit_MPI(0,'Error opening full_sphericalharmonic model file')
    endif

    read(IIN,*) lmax,kmax

    ! check maximum degree
    if (kmax-1 > NK_20 .or. lmax > NS_20) then
      print *,'Error full_sh model file: lmax or kmax too big ',lmax,kmax,' for file ',trim(rootdir)//trim(modelname)
      call exit_MPI(0,'lmax or kmax too big in full_sh model file')
    endif

    ! user output
    write(IMAIN,*) '  name ',trim(modelname),': lmax = ',lmax,' kmax = ',kmax
    call flush_IMAIN()

    ! reads in spherical harmonics coefficients
    nm = (lmax+1)**2

    select case(ipar)
    case (1)
      do k = 0,kmax-1
        read(IIN,*) (MANTLE_SH_V_dvsv(k,l),l=1,nm)
      enddo
    case (2)
      do k = 0,kmax-1
        read(IIN,*) (MANTLE_SH_V_dvsh(k,l),l=1,nm)
      enddo
    case (3)
      do k = 0,kmax-1
        read(IIN,*) (MANTLE_SH_V_dvpv(k,l),l=1,nm)
      enddo
    case (4)
      do k = 0,kmax-1
        read(IIN,*) (MANTLE_SH_V_dvph(k,l),l=1,nm)
      enddo
    case (5)
      do k = 0,kmax-1
        read(IIN,*) (MANTLE_SH_V_deta(k,l),l=1,nm)
      enddo
    case (6)
      do k = 0,kmax-1
        read(IIN,*) (MANTLE_SH_V_drho(k,l),l=1,nm)
      enddo
    end select

    close(IIN)

  enddo

  ! internal topography
  do ipar = 1,3

    select case(ipar)
    case (1)
      ! T410 degree 20 model
      modelname = 'T410.dat'
    case (2)
      ! T660 degree 20 model
      modelname = 'T660.dat'
    case (3)
      ! TCMB degree 20 model
      modelname = 'TCMB.dat'
    case default
      stop 'Invalid parameter for internal topography in mantle_sh'
    end select


    ! opens model file
    open(unit=IIN,file=trim(rootdir)//trim(modelname),status='old',action='read',iostat=ier)
    if (ier /= 0) then
      print *,'Error opening file: ',trim(rootdir)//trim(modelname)
      call exit_MPI(0,'Error opening full_sphericalharmonic model file')
    endif

    read(IIN,*) lmax

    ! check maximum degree
    if (lmax > NS_20) then
      print *,'Error full_sh model file: lmax too big ',lmax,' for file ',trim(rootdir)//trim(modelname)
      call exit_MPI(0,'lmax too big in full_sh model file')
    endif

    ! user output
    write(IMAIN,*) '  name ',trim(modelname),': lmax = ',lmax
    call flush_IMAIN()

    ! reads in spherical harmonics coefficients
    nm = (lmax+1)**2

    select case(ipar)
    case (1)
      read(IIN,*) (MANTLE_SH_V_t410(l),l=1,nm)
    case (2)
      read(IIN,*) (MANTLE_SH_V_t660(l),l=1,nm)
    case (3)
      read(IIN,*) (MANTLE_SH_V_tcmb(l),l=1,nm)
    end select

    close(IIN)

  enddo

  ! set up the splines used as radial basis functions by Ritsema
  call mantle_sh_splhsetup()

  ! topography statistics
  ! gets min/max value for a set of sample locations (regular lat/lon grid points)
  t410_min = HUGEVAL
  t410_max = -HUGEVAL
  t660_min = HUGEVAL
  t660_max = -HUGEVAL
  cmb_min = HUGEVAL
  cmb_max = -HUGEVAL
  do i = -90,90
    do j = -180,180
      ! lat/lon in degrees
      xlat = i * 1.0
      xlon = j * 1.0
      call ylm(xlat,xlon,NS_20,shcof)

      ! gets topos (in km)
      t410 = ZERO
      t660 = ZERO
      cmb = ZERO
      do l = 1,NSH_20
        t410 = t410 + MANTLE_SH_V_t410(l)*dble(shcof(l))
        t660 = t660 + MANTLE_SH_V_t660(l)*dble(shcof(l))
        cmb = cmb + MANTLE_SH_V_tcmb(l)*dble(shcof(l))
      enddo

      ! statistics
      if (t410 > t410_max) t410_max = t410
      if (t410 < t410_min) t410_min = t410
      if (t660 > t660_max) t660_max = t660
      if (t660 < t660_min) t660_min = t660
      if (cmb > cmb_max) cmb_max = cmb
      if (cmb < cmb_min) cmb_min = cmb
    enddo
  enddo

  ! user output
  write(IMAIN,*)
  write(IMAIN,*) '  410 topography min/max = ',sngl(t410_min),sngl(t410_max),' km'
  write(IMAIN,*) '  660 topography min/max = ',sngl(t660_min),sngl(t660_max),' km'
  write(IMAIN,*) '  CMB topography min/max = ',sngl(cmb_min),sngl(cmb_max),' km'
  call flush_IMAIN()

  end subroutine read_model_mantle_sh


!
!-----------------------------------------------------------------------------------------
!

  subroutine mantle_sh(lat,lon,radius,dvpv,dvph,dvsv,dvsh,deta,drho)

  use constants, only: ZERO
  use shared_parameters, only: R_PLANET

  use model_full_sh_mantle_par

  implicit none

  double precision :: lat,lon,radius,dvpv,dvph,dvsv,dvsh,deta,drho

! factor to convert perturbations in shear speed to perturbations in density
  double precision, parameter :: SCALE_RHO = 0.40d0

  double precision, parameter :: RMOHO_ = 6346600.d0
  double precision, parameter :: RCMB_ = 3480000.d0

  integer :: l,k
  double precision :: r_moho,r_cmb,xr
  double precision :: mantle_sh_rsple,radial_basis(0:NK_20),xmap(NSH_20)
  real(kind=4) :: shcof(NSH_20)
  real(kind=4) :: xlat,xlon

  dvsv = ZERO
  dvsh = ZERO
  dvpv = ZERO
  dvph = ZERO
  drho = ZERO
  deta = ZERO

  r_moho = RMOHO_ / R_PLANET
  r_cmb = RCMB_ / R_PLANET
  if (radius >= r_moho .or. radius <= r_cmb) return

  ! get spherical harmonics coefficients
  xlat = sngl(lat)
  xlon = sngl(lon)
  call ylm(xlat,xlon,NS_20,shcof)

  xr = -1.0d0+2.0d0*(radius-r_cmb)/(r_moho-r_cmb)
  do k = 0,NK_20
    radial_basis(k) = mantle_sh_rsple(1,NK_20+1, &
                                      MANTLE_SH_V_spknt(1),MANTLE_SH_V_qq0(1,NK_20+1-k),MANTLE_SH_V_qq(1,1,NK_20+1-k),xr)
  enddo

  ! Vsv perturbation
  do l = 1,NSH_20
    xmap(l) = 0.0d0
    do k = 0,NK_20
      xmap(l) = xmap(l) + MANTLE_SH_V_dvsv(k,l)*radial_basis(k)
    enddo
  enddo
  do l = 1,NSH_20
    dvsv = dvsv + xmap(l)*dble(shcof(l))
  enddo

  ! Vpv perturbation
  do l = 1,NSH_20
    xmap(l) = 0.0d0
    do k = 0,NK_20
      xmap(l) = xmap(l) + MANTLE_SH_V_dvpv(k,l)*radial_basis(k)
    enddo
  enddo
  do l = 1,NSH_20
    dvpv = dvpv + xmap(l)*dble(shcof(l))
  enddo

  ! isotropic perturbations
  dvsh = dvsv
  dvph = dvpv

  ! density perturbation scaled from Vsv
  drho = SCALE_RHO*dvsv

  end subroutine mantle_sh


!
!-----------------------------------------------------------------------------------------
!

  subroutine mantle_sh_splhsetup()

! sets up (spknt,qq0,qq)

  use model_full_sh_mantle_par

  implicit none

  integer :: i,j
  double precision :: qqwk(3,NK_20+1)

  MANTLE_SH_V_spknt(1) = -1.00000d0
  MANTLE_SH_V_spknt(2) = -0.78631d0
  MANTLE_SH_V_spknt(3) = -0.59207d0
  MANTLE_SH_V_spknt(4) = -0.41550d0
  MANTLE_SH_V_spknt(5) = -0.25499d0
  MANTLE_SH_V_spknt(6) = -0.10909d0
  MANTLE_SH_V_spknt(7) = 0.02353d0
  MANTLE_SH_V_spknt(8) = 0.14409d0
  MANTLE_SH_V_spknt(9) = 0.25367d0
  MANTLE_SH_V_spknt(10) = 0.35329d0
  MANTLE_SH_V_spknt(11) = 0.44384d0
  MANTLE_SH_V_spknt(12) = 0.52615d0
  MANTLE_SH_V_spknt(13) = 0.60097d0
  MANTLE_SH_V_spknt(14) = 0.66899d0
  MANTLE_SH_V_spknt(15) = 0.73081d0
  MANTLE_SH_V_spknt(16) = 0.78701d0
  MANTLE_SH_V_spknt(17) = 0.83810d0
  MANTLE_SH_V_spknt(18) = 0.88454d0
  MANTLE_SH_V_spknt(19) = 0.92675d0
  MANTLE_SH_V_spknt(20) = 0.96512d0
  MANTLE_SH_V_spknt(21) = 1.00000d0

  do i = 1,NK_20+1
    do j = 1,NK_20+1
      if (i == j) then
        MANTLE_SH_V_qq0(j,i) = 1.0d0
      else
        MANTLE_SH_V_qq0(j,i) = 0.0d0
      endif
    enddo
  enddo
  do i = 1,NK_20+1
    call mantle_sh_rspln(1,NK_20+1,MANTLE_SH_V_spknt(1),MANTLE_SH_V_qq0(1,i),MANTLE_SH_V_qq(1,1,i),qqwk(1,1))
  enddo

  end subroutine mantle_sh_splhsetup


!
!-----------------------------------------------------------------------------------------
!

! changed the obsolecent f77 features in the two routines below
! now still awful Fortran, but at least conforms to f90 standard

  double precision function mantle_sh_rsple(I1,I2,X,Y,Q,S)

  implicit none

! rsple returns the value of the function y(x) evaluated at point S
! using the cubic spline coefficients computed by rspln and saved in Q.
! If S is outside the interval (x(i1),x(i2)) rsple extrapolates
! using the first or last interpolation polynomial. The arrays must
! be dimensioned at least - x(i2), y(i2), and q(3,i2).

  integer :: i1,i2
  double precision :: X(21),Y(21),Q(3,21),S

  integer :: i,ii
  double precision :: h

  i = 1
  II = I2-1

  !   GUARANTEE I WITHIN BOUNDS.
  I = MAX0(I,I1)
  I = MIN0(I,II)

  !   SEE IF X IS INCREASING OR DECREASING.
  if (X(I2)-X(I1) < 0) goto 1
  if (X(I2)-X(I1) >= 0) goto 2

  !   X IS DECREASING.  CHANGE I AS NECESSARY.
1  if (S-X(I) <= 0) goto 3
  if (S-X(I) > 0) goto 4

4  I = I-1

  if (I-I1 < 0) goto 11
  if (I-I1 == 0) goto 6
  if (I-I1 > 0) goto 1

3  if (S-X(I+1) < 0) goto 5
  if (S-X(I+1) >= 0) goto 6

5  I = I+1

  if (I-II < 0) goto 3
  if (I-II == 0) goto 6
  if (I-II > 0) goto 7

  !   X IS INCREASING.  CHANGE I AS NECESSARY.
2  if (S-X(I+1) <= 0) goto 8
  if (S-X(I+1) > 0) goto 9

9  I = I+1

  if (I-II < 0) goto 2
  if (I-II == 0) goto 6
  if (I-II > 0) goto 7

8  if (S-X(I) < 0) goto 10
  if (S-X(I) >= 0) goto 6

10 I = I-1
  if (I-I1 < 0) goto 11
  if (I-I1 == 0) goto 6
  if (I-I1 > 0) goto 8

7  I = II
  goto 6

11 I = I1

  !   CALCULATE RSPLE USING SPLINE COEFFICIENTS IN Y AND Q.
6  H = S-X(I)
  MANTLE_SH_RSPLE = Y(I)+H*(Q(1,I)+H*(Q(2,I)+H*Q(3,I)))

  end function mantle_sh_rsple


!
!-----------------------------------------------------------------------------------------
!

  subroutine mantle_sh_rspln(I1,I2,X,Y,Q,F)

  implicit none

! The subroutine rspln computes cubic spline interpolation coefficients
! for y(x) between grid points i1 and i2 saving them in q.The
! interpolation is continuous with continuous first and second
! derivatives. It agrees exactly with y at grid points and with the
! three point first derivatives at both end points (i1 and i2).
! X must be monotonic but if two successive values of x are equal
! a discontinuity is assumed and separate interpolation is done on
! each strictly monotonic segment. The arrays must be dimensioned at
! least - x(i2), y(i2), q(3,i2), and f(3,i2).
! F is working storage for rspln.

  integer :: i1,i2
  double precision :: X(21),Y(21),Q(3,21),F(3,21)

  integer :: i,j,k,j1,j2
  double precision :: y0,a0,b0,b1,h,h2,ha,h2a,h3a,h2b
  double precision :: YY(3),small

  equivalence (YY(1),Y0)

  data SMALL/1.0d-08/,YY/0.0d0,0.0d0,0.0d0/

  J1 = I1+1
  Y0 = 0.0d0

  !   BAIL OUT IF THERE ARE LESS THAN TWO POINTS TOTAL
  if (I2-I1 < 0) return
  if (I2-I1 == 0) goto 17
  if (I2-I1 > 0) goto 8

8  A0 = X(J1-1)

  !   SEARCH FOR DISCONTINUITIES.
  do 3 I = J1,I2
    B0 = A0
    A0 = X(I)
    if (DABS((A0-B0)/DMAX1(A0,B0)) < SMALL) goto 4

3   continue

17 J1 = J1-1
  J2 = I2-2
  goto 5

4  J1 = J1-1
  J2 = I-3

  !   SEE IF THERE ARE ENOUGH POINTS TO INTERPOLATE (AT LEAST THREE).
5  if (J2+1-J1 < 0) goto 9
  if (J2+1-J1 == 0) goto 10
  if (J2+1-J1 > 0) goto 11

  !   ONLY TWO POINTS.  USE LINEAR INTERPOLATION.
10  J2 = J2+2
  Y0 = (Y(J2)-Y(J1))/(X(J2)-X(J1))
  do J = 1,3
    Q(J,J1) = YY(J)
    Q(J,J2) = YY(J)
  enddo
  goto 12

    !   MORE THAN TWO POINTS.  DO SPLINE INTERPOLATION.
11  A0 = 0.
  H = X(J1+1)-X(J1)
  H2 = X(J1+2)-X(J1)
  Y0 = H*H2*(H2-H)
  H = H*H
  H2 = H2*H2

  ! CALCULATE DERIVITIVE AT NEAR END.
  B0 = (Y(J1)*(H-H2)+Y(J1+1)*H2-Y(J1+2)*H)/Y0
  B1 = B0

  !   EXPLICITLY REDUCE BANDED MATRIX TO AN UPPER BANDED MATRIX.
  do I = J1,J2
    H = X(I+1)-X(I)
    Y0 = Y(I+1)-Y(I)
    H2 = H*H
    HA = H-A0
    H2A = H-2.0d0*A0
    H3A = 2.0d0*H-3.0d0*A0
    H2B = H2*B0
    Q(1,I) = H2/HA
    Q(2,I) = -HA/(H2A*H2)
    Q(3,I) = -H*H2A/H3A
    F(1,I) = (Y0-H*B0)/(H*HA)
    F(2,I) = (H2B-Y0*(2.0d0*H-A0))/(H*H2*H2A)
    F(3,I) = -(H2B-3.0d0*Y0*HA)/(H*H3A)
    A0 = Q(3,I)
    B0 = F(3,I)
  enddo

  !   TAKE CARE OF LAST TWO ROWS.
  I = J2+1
  H = X(I+1)-X(I)
  Y0 = Y(I+1)-Y(I)
  H2 = H*H
  HA = H-A0
  H2A = H*HA
  H2B = H2*B0-Y0*(2.0d0*H-A0)
  Q(1,I) = H2/HA
  F(1,I) = (Y0-H*B0)/H2A
  HA = X(J2)-X(I+1)
  Y0 = -H*HA*(HA+H)
  HA = HA*HA

  !   CALCULATE DERIVATIVE AT FAR END.
  Y0 = (Y(I+1)*(H2-HA)+Y(I)*HA-Y(J2)*H2)/Y0
  Q(3,I) = (Y0*H2A+H2B)/(H*H2*(H-2.0d0*A0))
  Q(2,I) = F(1,I)-Q(1,I)*Q(3,I)

  !   SOLVE UPPER BANDED MATRIX BY REVERSE ITERATION.
  do J = J1,J2
    K = I-1
    Q(1,I) = F(3,K)-Q(3,K)*Q(2,I)
    Q(3,K) = F(2,K)-Q(2,K)*Q(1,I)
    Q(2,K) = F(1,K)-Q(1,K)*Q(3,K)
    I = K
  enddo
  Q(1,I) = B1

  !   FILL IN THE LAST POINT WITH A LINEAR EXTRAPOLATION.
9  J2 = J2+2
  do J = 1,3
    Q(J,J2) = YY(J)
  enddo

  !   SEE IF THIS DISCONTINUITY IS THE LAST.
12 if (J2-I2 < 0) then
    goto 6
  else
    return
  endif

  !   NO.  GO BACK FOR MORE.
6  J1 = J2+2
  if (J1-I2 <= 0) goto 8
  if (J1-I2 > 0) goto 7

  !   THERE IS ONLY ONE POINT LEFT AFTER THE LATEST DISCONTINUITY.
7 do J = 1,3
    Q(J,I2) = YY(J)
  enddo

  end subroutine mantle_sh_rspln

!
!-------------------------------------------------------------------------------------------------
!

! added by JT 2015

  subroutine add_topography_sh_mantle(xelm,yelm,zelm)

  use constants
  use shared_parameters, only: R_PLANET,ELLIPTICITY

  use meshfem_par, only: R220,R400,R670,R771

  implicit none

  double precision, intent(inout) :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  ! local parameters
  integer :: ia

  real(kind=4) :: topo410out,topo650out
  double precision :: topo410,topo650

  double precision :: r,lat,lon
  double precision :: gamma
  double precision :: x,y,z

  !statistics
  logical,parameter :: DEBUG_STATISTICS = .false.
  real(kind=CUSTOM_REAL), parameter :: HUGEVAL_REAL = real(HUGEVAL,kind=CUSTOM_REAL)
  real(kind=CUSTOM_REAL),save :: min_410 = HUGEVAL_REAL, max_410 = - HUGEVAL_REAL
  real(kind=CUSTOM_REAL),save :: min_650 = HUGEVAL_REAL, max_650 = - HUGEVAL_REAL
  real(kind=CUSTOM_REAL) :: min_410_all,max_410_all
  real(kind=CUSTOM_REAL) :: min_650_all,max_650_all

! note: adding topography to 410 and 650 strongly affects PcP, PKiKP, etc. phases,
!       we leave it in and check whether the stretching makes simulation unstable
!
! topography perturbations
!   410-km: minimum / maximum = -13.48 km / + 13.24 km
!   650-km: minimum / maximum = -14.34 km / + 19.19 km

! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! note: the topography on 410 and 650 is given in geographic colat/lon,
    !       thus we need to convert geocentric colatitude to geographic colatitudes
    !
    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    ! note: at this point, the mesh is still spherical (no need to correct latitude for ellipticity)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

    ! stretching occurs between 220 and 770
    if (r > R220/R_PLANET .or. r < R771/R_PLANET) cycle

    topo410 = 0.d0
    topo650 = 0.d0

    ! compute topography on 410 and 660 at current point
    call subtopo_sh(lat,lon,topo410,topo650)

    ! debug
    !write(*,*) 'ola',lat,lon,topo410,topo650

    if (topo410 == 0.d0 .and. topo650 == 0.d0) return

    ! min/max statistics
    topo410out = sngl(topo410)
    topo650out = sngl(topo650)
    if (DEBUG_STATISTICS) then
      if (topo410out < min_410) min_410 = topo410out
      if (topo410out > max_410) max_410 = topo410out
      if (topo650out < min_650) min_650 = topo650out
      if (topo650out > max_650) max_650 = topo650out
      ! debug
      !print *,'topo410 / topo650: ',r,lat,lon,topo410out,topo650out
    endif

    ! non-dimensionalize the topography, which is in km
    ! positive for a depression, so change the sign for a perturbation in radius
    topo410 = -(topo410) / (R_PLANET/1000.d0)
    topo650 = -(topo650) / (R_PLANET/1000.d0)

    gamma = 0.d0
    if (r >= R400/R_PLANET .and. r <= R220/R_PLANET) then
      ! stretching between R220 and R400
      gamma = (R220/R_PLANET - r) / (R220/R_PLANET - R400/R_PLANET)
      xelm(ia) = x*(ONE + gamma * topo410 / r)
      yelm(ia) = y*(ONE + gamma * topo410 / r)
      zelm(ia) = z*(ONE + gamma * topo410 / r)
    else if (r >= R771/R_PLANET .and. r <= R670/R_PLANET) then
      ! stretching between R771 and R670
      gamma = (r - R771/R_PLANET) / (R670/R_PLANET - R771/R_PLANET)
      xelm(ia) = x*(ONE + gamma * topo650 / r)
      yelm(ia) = y*(ONE + gamma * topo650 / r)
      zelm(ia) = z*(ONE + gamma * topo650 / r)
    else if (r > R670/R_PLANET .and. r < R400/R_PLANET) then
      ! stretching between R670 and R400
      gamma = (R400/R_PLANET - r) / (R400/R_PLANET - R670/R_PLANET)
      xelm(ia) = x*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      yelm(ia) = y*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
      zelm(ia) = z*(ONE + (topo410 + gamma * (topo650 - topo410)) / r)
    endif
    if (gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for 410-650 topography')

  enddo

  ! debug
  if (DEBUG_STATISTICS) then
    ! collects min/max on main
    call min_all_cr(min_410,min_410_all)
    call max_all_cr(max_410,max_410_all)
    call min_all_cr(min_650,min_650_all)
    call max_all_cr(max_650,max_650_all)
    if (myrank == 0) then
      if (r <= R220/R_PLANET .and. r >= R771/R_PLANET) then
        print *,'add_topography_410_650: min/max_410 = ',min_410_all,max_410_all,'min/max_650 = ',min_650_all,max_650_all
      endif
    endif
    !if (r <= R220/R_PLANET .and. r >= R771/R_PLANET) then
    !  print *,myrank,'add_topography_410_650: min/max_410 = ',min_410,max_410,'min/max_650 = ',min_650,max_650
    !  print *,myrank,'add_topography_410_650: depth = ',(1.d0 - r)*(R_PLANET/1000.d0), &
    !         ' 410-km = ',topo410out,' 650-km = ',topo650out
    !endif
  endif

  end subroutine add_topography_sh_mantle

!
!-------------------------------------------------------------------------------------------------
!

  subroutine subtopo_sh(lat,lon,topo410,topo660)

  use model_full_sh_mantle_par, only: NSH_20,NS_20,MANTLE_SH_V_t410,MANTLE_SH_V_t660

  implicit none

  double precision,intent(in) :: lat,lon
  double precision,intent(out) :: topo410,topo660

  ! local parameter
  integer :: l
  real(kind=4) :: shcof(NSH_20)
  real(kind=4) :: xlat,xlon

  xlat = sngl(lat)
  xlon = sngl(lon)
  call ylm(xlat,xlon,NS_20,shcof)

  topo410 = 0.d0
  topo660 = 0.d0
  do l = 1,NSH_20
    topo410 = topo410 + MANTLE_SH_V_t410(l)*dble(shcof(l))
    topo660 = topo660 + MANTLE_SH_V_t660(l)*dble(shcof(l))
  enddo

  end subroutine subtopo_sh

!
!-------------------------------------------------------------------------------------------------
!

  subroutine add_topography_sh_cmb(xelm,yelm,zelm)

! this is only a placeholder function, which is not used yet...user must supply the subtopo_cmb() routine

  use constants
  use shared_parameters, only: R_PLANET,RCMB,RTOPDDOUBLEPRIME,ELLIPTICITY

  implicit none

  double precision,intent(inout) :: xelm(NGNOD),yelm(NGNOD),zelm(NGNOD)

  ! PREM reference values
  double precision :: RTOPDDOUBLEPRIME_ = 3630000.d0
  double precision :: RCMB_ = 3480000.d0

  ! local parameters
  integer :: ia
  double precision :: r_start,topocmb
  double precision :: r,lat,lon
  double precision :: x,y,z
  double precision :: gamma

  ! checks if we have right reference values
  if (RCMB /= RCMB_) stop 'Error add_topography_sh_cmb has wrong RCMB reference'
  if (RTOPDDOUBLEPRIME /= RTOPDDOUBLEPRIME_) stop 'Error add_topography_sh_cmb has wrong RTOPDDOUBLEPRIME reference'

  ! we loop on all the points of the element
  do ia = 1,NGNOD

    x = xelm(ia)
    y = yelm(ia)
    z = zelm(ia)

    ! converts geocentric coordinates x/y/z to geographic radius/latitude/longitude (in degrees)
    ! note: at this point, the mesh is still spherical (no need to correct latitude for ellipticity)
    call xyz_2_rlatlon_dble(x,y,z,r,lat,lon,ELLIPTICITY)

    ! compute topography on CMB; routine subtopo_cmb needs to be supplied by the user
    call subtopo_sh_cmb(lat,lon,topocmb)

    ! checks if anything to do
    if (topocmb == 0.0d0) cycle

    ! non-dimensionalize the topography, which is in km
    ! positive for a depression, so change the sign for a perturbation in radius
    topocmb = -topocmb / (R_PLANET/1000.d0)

    ! start stretching a distance RTOPDDOUBLEPRIME - RCMB below the CMB
    ! and finish at RTOPDDOUBLEPRIME of D_double_prime
    r_start = (RCMB - (RTOPDDOUBLEPRIME - RCMB))/R_PLANET
    gamma = 0.0d0
    if (r >= RCMB/R_PLANET .and. r <= RTOPDDOUBLEPRIME/R_PLANET) then
      ! stretching between RCMB and RTOPDDOUBLEPRIME
      gamma = (RTOPDDOUBLEPRIME/R_PLANET - r) / (RTOPDDOUBLEPRIME/R_PLANET - RCMB/R_PLANET)
    else if (r >= r_start .and. r <= RCMB/R_PLANET) then
      ! stretching between r_start and RCMB
      gamma = (r - r_start) / (RCMB/R_PLANET - r_start)
    endif
    if (gamma < -0.0001 .or. gamma > 1.0001) call exit_MPI(myrank,'incorrect value of gamma for CMB topography')

    xelm(ia) = x*(ONE + gamma * topocmb / r)
    yelm(ia) = y*(ONE + gamma * topocmb / r)
    zelm(ia) = z*(ONE + gamma * topocmb / r)

  enddo

  end subroutine add_topography_sh_cmb

!
!-------------------------------------------------------------------------------------------------
!


!! JT dec 2015
  subroutine subtopo_sh_cmb(lat,lon,topocmb)

! note: uses reference values
!  RTOPDDOUBLEPRIME = 3630000.d0
!  RCMB = 3480000.d0

  use model_full_sh_mantle_par, only: NSH_20,NS_20,MANTLE_SH_V_tcmb

  implicit none

  double precision,intent(in) :: lat,lon
  double precision,intent(out) :: topocmb

  ! local parameters
  integer :: l
  real(kind=4) :: shcof(NSH_20)
  real(kind=4) :: xlat,xlon

  xlat = sngl(lat)
  xlon = sngl(lon)
  call ylm(xlat,xlon,NS_20,shcof)

  topocmb = 0.0d0
  do l = 1,NSH_20
    topocmb = topocmb + MANTLE_SH_V_tcmb(l)*dble(shcof(l))
  enddo

  end subroutine subtopo_sh_cmb

