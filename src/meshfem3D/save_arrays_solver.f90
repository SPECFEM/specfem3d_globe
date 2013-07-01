!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            April 2011
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
                    etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore, &
                    xstore,ystore,zstore,rhostore,dvpstore, &
                    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    nspec_ani,c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    ibool,idoubling,is_on_a_slice_edge,nglob_xy,nglob, &
                    rmassx,rmassy,rmassz,rmass_ocean_load,npointot_oceans, &
                    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
                    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
                    normal_xmin,normal_xmax,normal_ymin,normal_ymax,normal_bottom,normal_top, &
                    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
                    jacobian2D_bottom,jacobian2D_top,nspec, &
                    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP, &
                    TRANSVERSE_ISOTROPY,HETEROGEN_3D_MANTLE,ANISOTROPIC_3D_MANTLE, &
                    ANISOTROPIC_INNER_CORE,OCEANS, &
                    tau_s,tau_e_store,Qmu_store,T_c_source,ATTENUATION,vx,vy,vz,vnspec, &
                    NCHUNKS,ABSORBING_CONDITIONS,SAVE_MESH_FILES,ispec_is_tiso)

  implicit none

  include "constants.h"

  logical ATTENUATION

  character(len=150) prname
  integer iregion_code,NCHUNKS

  integer nspec,nglob_xy,nglob,nspec_stacey
  integer npointot_oceans

! Stacey
  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec_stacey)
  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,nspec_stacey)

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

  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  ! doubling mesh flag
  integer, dimension(nspec) :: idoubling

! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

! mass matrices
! use rmassz as temporary storage for the rest of this routine
! once it is written to disk and thus not needed any more
  real(kind=CUSTOM_REAL), dimension(nglob_xy) :: rmassx,rmassy
  real(kind=CUSTOM_REAL), dimension(nglob)    :: rmassz

! additional ocean load mass matrix
  real(kind=CUSTOM_REAL) rmass_ocean_load(npointot_oceans)

! boundary parameters locator
  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM,NSPEC2D_TOP

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

! attenuation
  integer vx, vy, vz, vnspec
  double precision  T_c_source
  double precision, dimension(N_SLS)                     :: tau_s
  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec)        :: Qmu_store
  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: tau_e_store

  logical ABSORBING_CONDITIONS,SAVE_MESH_FILES

  logical, dimension(nspec) :: ispec_is_tiso

  ! local parameters
  integer i,j,k,ispec,iglob,nspec1, nglob1
  real(kind=CUSTOM_REAL) scaleval1,scaleval2

! save nspec and nglob, to be used in combine_paraview_data
  open(unit=27,file=prname(1:len_trim(prname))//'array_dims.txt',status='unknown',action='write')

  nspec1 = nspec
  nglob1 = nglob

  ! might be wrong, check...
  !if (NCHUNKS == 6 .and. ichunk /= CHUNK_AB .and. iregion_code == IREGION_INNER_CORE) then
  !  ! only chunk_AB contains inner core?
  !  ratio_divide_central_cube = 16
  !  ! corrects nspec/nglob
  !  nspec1 = nspec1 - (NEX_PER_PROC_XI/ratio_divide_central_cube) &
  !            * (NEX_PER_PROC_ETA/ratio_divide_central_cube) * (NEX_XI/ratio_divide_central_cube)
  !  nglob1 = nglob1 -   ((NEX_PER_PROC_XI/ratio_divide_central_cube)*(NGLLX-1)+1) &
  !            * ((NEX_PER_PROC_ETA/ratio_divide_central_cube)*(NGLLY-1)+1) &
  !            * (NEX_XI/ratio_divide_central_cube)*(NGLLZ-1)
  !endif

  write(27,*) nspec1
  write(27,*) nglob1
  close(27)

  open(unit=27,file=prname(1:len_trim(prname))//'solver_data_1.bin',status='unknown',form='unformatted',action='write')

  ! local GLL points
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

    ! note: muvstore needed for Q_mu shear attenuation in inner core
    if(.not. (ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE)) then
      write(27) muvstore
    endif

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

  ! mass matrices
  !
  ! in the case of stacey boundary conditions, add C*deltat/2 contribution to the mass matrix
  ! on the Stacey edges for the crust_mantle and outer_core regions but not for the inner_core region
  ! thus the mass matrix must be replaced by three mass matrices including the "C" damping matrix
  !
  ! if absorbing_conditions are not set or if NCHUNKS=6, only one mass matrix is needed
  ! for the sake of performance, only "rmassz" array will be filled and "rmassx" & "rmassy" will be obsolete
  if(NCHUNKS /= 6 .and. ABSORBING_CONDITIONS .and. iregion_code == IREGION_CRUST_MANTLE) then
     write(27) rmassx
     write(27) rmassy
  endif

  write(27) rmassz

  ! additional ocean load mass matrix if oceans and if we are in the crust
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) write(27) rmass_ocean_load

  close(27) ! solver_data.bin

  open(unit=27,file=prname(1:len_trim(prname))//'solver_data_2.bin',status='unknown',form='unformatted',action='write')

! mesh arrays used in the solver to locate source and receivers
! and for anisotropy and gravity, save in single precision
! use rmassz for temporary storage to perform conversion, since already saved

!--- x coordinate
  rmassz(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmassz(iglob) = sngl(xstore(i,j,k,ispec))
          else
            rmassz(iglob) = xstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmassz

!--- y coordinate
  rmassz(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmassz(iglob) = sngl(ystore(i,j,k,ispec))
          else
            rmassz(iglob) = ystore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmassz

!--- z coordinate
  rmassz(:) = 0._CUSTOM_REAL
  do ispec = 1,nspec
    do k = 1,NGLLZ
      do j = 1,NGLLY
        do i = 1,NGLLX
          iglob = ibool(i,j,k,ispec)
! distinguish between single and double precision for reals
          if(CUSTOM_REAL == SIZE_REAL) then
            rmassz(iglob) = sngl(zstore(i,j,k,ispec))
          else
            rmassz(iglob) = zstore(i,j,k,ispec)
          endif
        enddo
      enddo
    enddo
  enddo
  write(27) rmassz

  write(27) ibool

  write(27) idoubling

  write(27) is_on_a_slice_edge

  write(27) ispec_is_tiso

  close(27)

  ! absorbing boundary parameters
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_mpi(myrank,'error opening boundary.bin file')

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
  if(ATTENUATION) then
    open(unit=27, file=prname(1:len_trim(prname))//'attenuation.bin', &
          status='unknown', form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening attenuation.bin file')

    write(27) tau_s
    write(27) tau_e_store
    write(27) Qmu_store
    write(27) T_c_source
    close(27)
  endif

  ! uncomment for vp & vs model storage
  if( SAVE_MESH_FILES ) then
    scaleval1 = sngl( sqrt(PI*GRAV*RHOAV)*(R_EARTH/1000.0d0) )
    scaleval2 = sngl( RHOAV/1000.0d0 )

    ! isotropic model
    ! vp
    open(unit=27,file=prname(1:len_trim(prname))//'vp.bin',status='unknown',form='unformatted',action='write')
    write(27) sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1
    close(27)
    ! vs
    open(unit=27,file=prname(1:len_trim(prname))//'vs.bin',status='unknown',form='unformatted',action='write')
    write(27) sqrt( muvstore/rhostore )*scaleval1
    close(27)
    ! rho
    open(unit=27,file=prname(1:len_trim(prname))//'rho.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening rho.bin file')

    write(27) rhostore*scaleval2
    close(27)

    ! transverse isotropic model
    if( TRANSVERSE_ISOTROPY ) then
      ! vpv
      open(unit=27,file=prname(1:len_trim(prname))//'vpv.bin',status='unknown',form='unformatted',action='write')
      write(27) sqrt( (kappavstore+4.*muvstore/3.)/rhostore )*scaleval1
      close(27)
      ! vph
      open(unit=27,file=prname(1:len_trim(prname))//'vph.bin',status='unknown',form='unformatted',action='write')
      write(27) sqrt( (kappahstore+4.*muhstore/3.)/rhostore )*scaleval1
      close(27)
      ! vsv
      open(unit=27,file=prname(1:len_trim(prname))//'vsv.bin',status='unknown',form='unformatted',action='write')
      write(27) sqrt( muvstore/rhostore )*scaleval1
      close(27)
      ! vsh
      open(unit=27,file=prname(1:len_trim(prname))//'vsh.bin',status='unknown',form='unformatted',action='write')
      write(27) sqrt( muhstore/rhostore )*scaleval1
      close(27)
      ! rho
      open(unit=27,file=prname(1:len_trim(prname))//'rho.bin',status='unknown',form='unformatted',action='write')
      write(27) rhostore*scaleval2
      close(27)
      ! eta
      open(unit=27,file=prname(1:len_trim(prname))//'eta.bin',status='unknown',form='unformatted',action='write')
      write(27) eta_anisostore
      close(27)
    endif

  endif

  end subroutine save_arrays_solver




