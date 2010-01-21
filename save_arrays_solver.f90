!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

  subroutine save_arrays_solver(rho_vp,rho_vs,nspec_stacey, &
            prname,iregion_code,xixstore,xiystore,xizstore, &
            etaxstore,etaystore,etazstore, &
            gammaxstore,gammaystore,gammazstore, &
            xstore,ystore,zstore,rhostore,dvpstore, &
            kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
            nspec_ani, &
            c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
            c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
            c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
            ibool,idoubling,rmass,rmass_ocean_load,npointot_oceans, &
            ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
            nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
            normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
            jacobian2D_xmin,jacobian2D_xmax, &
            jacobian2D_ymin,jacobian2D_ymax, &
            jacobian2D_bottom,jacobian2D_top, &
            nspec,nglob, &
            NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
            TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS, &
            tau_s,tau_e_store,Qmu_store,T_c_source, &
            ATTENUATION,vx,vy,vz,vnspec, &
            NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NEX_XI,ichunk,NCHUNKS,ABSORBING_CONDITIONS)

  implicit none

  include "constants.h"

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(:), pointer   :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qomsb, Qomsb2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

!> Hejun
!  type (attenuation_model_variables) AM_V
! attenuation_model_variables
!< Hejun

  logical ATTENUATION,ABSORBING_CONDITIONS
!  logical ATTENUATION_3D

  integer NEX_PER_PROC_XI,NEX_PER_PROC_ETA,NEX_XI,ichunk

  integer nspec,nglob,nspec_stacey,NCHUNKS
  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP
  integer npointot_oceans

  logical TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS

! arrays with jacobian matrix
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore

! arrays with mesh parameters
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

! for anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    rhostore,dvpstore,kappavstore,kappahstore,muvstore,muhstore,eta_anisostore

  integer nspec_ani

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store,c33store,c34store, &
    c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store

! Stacey
  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey)
  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey)

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

! doubling mesh flag
  integer idoubling(nspec)

! mass matrix
  real(kind=CUSTOM_REAL) rmass(nglob)

! additional ocean load mass matrix
  real(kind=CUSTOM_REAL) rmass_ocean_load(npointot_oceans)

! boundary parameters locator
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)
  integer ibelm_bottom(NSPEC2D_BOTTOM),ibelm_top(NSPEC2D_TOP)

! normals
  real(kind=CUSTOM_REAL) normal_xmin(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_xmax(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) normal_ymin(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_ymax(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) normal_bottom(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) normal_top(NDIM,NGLLX,NGLLY,NSPEC2D_TOP)

! jacobian on 2D edges
  real(kind=CUSTOM_REAL) jacobian2D_xmin(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_xmax(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymin(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_ymax(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX)
  real(kind=CUSTOM_REAL) jacobian2D_bottom(NGLLX,NGLLY,NSPEC2D_BOTTOM)
  real(kind=CUSTOM_REAL) jacobian2D_top(NGLLX,NGLLY,NSPEC2D_TOP)

! number of elements on the boundaries
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  integer i,j,k,ispec,iglob,nspec1, nglob1

! attenuation
  integer vx, vy, vz, vnspec
  double precision  T_c_source
  double precision, dimension(N_SLS)                     :: tau_s
  double precision, dimension(vx, vy, vz, vnspec)        :: Qmu_store
  double precision, dimension(N_SLS, vx, vy, vz, vnspec) :: tau_e_store

! processor identification
  character(len=150) prname

  integer iregion_code

! save nspec and nglob, to be used in combine_paraview_data
  open(unit=27,file=prname(1:len_trim(prname))//'array_dims.txt',status='unknown',action='write')
  if (NCHUNKS == 6 .and. ichunk /= CHUNK_AB .and. iregion_code == IREGION_INNER_CORE) then
     nspec1 = nspec - (NEX_PER_PROC_XI/16) * (NEX_PER_PROC_ETA/16) * (NEX_XI/16)
     nglob1 = nglob -   ((NEX_PER_PROC_XI/16)*(NGLLX-1)+1) * ((NEX_PER_PROC_ETA/16)*(NGLLY-1)+1) &
       * (NEX_XI/16)*(NGLLZ-1)
  else
     nspec1 = nspec
     nglob1 = nglob
  endif
  write(27,*) nspec1
  write(27,*) nglob1
  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'solver_data_1.bin',status='unknown',form='unformatted',action='write')

  write(27) xixstore
  write(27) xiystore
  write(27) xizstore
  write(27) etaxstore
  write(27) etaystore
  write(27) etazstore
  write(27) gammaxstore
  write(27) gammaystore
  write(27) gammazstore

  write(27) rhostore
  write(27) kappavstore

  if(HETEROGEN_3D_MANTLE) then
     open(unit=29,file=prname(1:len_trim(prname))//'dvp.bin',status='unknown',form='unformatted',action='write')
     write(29) dvpstore
     close(29)
  endif

! other terms needed in the solid regions only
  if(iregion_code /= IREGION_OUTER_CORE) then

    if(.not. (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) write(27) muvstore

!   save anisotropy in the mantle only
    if(TRANSVERSE_ISOTROPY) then
      if(iregion_code == IREGION_CRUST_MANTLE .and. .not. ANISOTROPIC_3D_MANTLE) then
        write(27) kappahstore
        write(27) muhstore
        write(27) eta_anisostore
      endif
    endif

!   save anisotropy in the inner core only
    if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then
      write(27) c11store
      write(27) c33store
      write(27) c12store
      write(27) c13store
      write(27) c44store
    endif

    if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then
        write(27) c11store
        write(27) c12store
        write(27) c13store
        write(27) c14store
        write(27) c15store
        write(27) c16store
        write(27) c22store
        write(27) c23store
        write(27) c24store
        write(27) c25store
        write(27) c26store
        write(27) c33store
        write(27) c34store
        write(27) c35store
        write(27) c36store
        write(27) c44store
        write(27) c45store
        write(27) c46store
        write(27) c55store
        write(27) c56store
        write(27) c66store
    endif

  endif

! Stacey
  if(ABSORBING_CONDITIONS) then

    if(iregion_code == IREGION_CRUST_MANTLE) then
      write(27) rho_vp
      write(27) rho_vs
    else if(iregion_code == IREGION_OUTER_CORE) then
      write(27) rho_vp
    endif

  endif

! mass matrix
  write(27) rmass

! additional ocean load mass matrix if oceans and if we are in the crust
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) write(27) rmass_ocean_load

  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'solver_data_2.bin',status='unknown',form='unformatted',action='write')
! mesh arrays used in the solver to locate source and receivers
! and for anisotropy and gravity, save in single precision
! use rmass for temporary storage to perform conversion, since already saved

!--- x coordinate
  rmass(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglob) = sngl(xstore(i,j,k,ispec))
          else
            rmass(iglob) = xstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmass

!--- y coordinate
  rmass(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglob) = sngl(ystore(i,j,k,ispec))
          else
            rmass(iglob) = ystore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmass

!--- z coordinate
  rmass(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmass(iglob) = sngl(zstore(i,j,k,ispec))
          else
            rmass(iglob) = zstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo

  write(27) rmass

  write(27) ibool

  write(27) idoubling

  close(27)

! boundary parameters
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin',status='unknown',form='unformatted',action='write')

  write(27) nspec2D_xmin
  write(27) nspec2D_xmax
  write(27) nspec2D_ymin
  write(27) nspec2D_ymax
  write(27) NSPEC2D_BOTTOM
  write(27) NSPEC2D_TOP

  write(27) ibelm_xmin
  write(27) ibelm_xmax
  write(27) ibelm_ymin
  write(27) ibelm_ymax
  write(27) ibelm_bottom
  write(27) ibelm_top

  write(27) normal_xmin
  write(27) normal_xmax
  write(27) normal_ymin
  write(27) normal_ymax
  write(27) normal_bottom
  write(27) normal_top

  write(27) jacobian2D_xmin
  write(27) jacobian2D_xmax
  write(27) jacobian2D_ymin
  write(27) jacobian2D_ymax
  write(27) jacobian2D_bottom
  write(27) jacobian2D_top

  close(27)

!> Hejun
! No matter 1D or 3D Attenuation, we save value for gll points
!  if(ATTENUATION .and. ATTENUATION_3D) then
!     open(unit=27, file=prname(1:len_trim(prname))//'attenuation3D.bin', status='unknown', form='unformatted',action='write')
!     write(27) tau_s
!     write(27) tau_e_store
!     write(27) Qmu_store
!     write(27) T_c_source
!     close(27)
!  else if(ATTENUATION) then
!     call attenuation_save_arrays(prname, iregion_code, AM_V)
!  endif

  if(ATTENUATION) then
     open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', status='unknown', form='unformatted',action='write')
     write(27) tau_s
     write(27) tau_e_store
     write(27) Qmu_store
     write(27) T_c_source
     close(27)
  endif
!< Hejun

  end subroutine save_arrays_solver

