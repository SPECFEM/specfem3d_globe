!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

! read arrays created by the mesher

  subroutine read_arrays_solver(iregion_code,myrank, &
              rho_vp,rho_vs,xstore,ystore,zstore, &
              xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian, &
              rhostore, kappavstore,muvstore,kappahstore,muhstore,eta_anisostore, &
              nspec_iso,nspec_tiso,nspec_ani, &
              c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
              c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
              c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
              ibool,idoubling,rmass,rmass_ocean_load,nspec,nglob, &
              READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY, &
              ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS,LOCAL_PATH,NCHUNKS)

  implicit none

  include "constants.h"

  include "OUTPUT_FILES/values_from_mesher.h"

  integer iregion_code,myrank,NCHUNKS

! flags to know if we should read Vs and anisotropy arrays
  logical READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY,ANISOTROPIC_3D_MANTLE,ANISOTROPIC_INNER_CORE,OCEANS

  character(len=150) LOCAL_PATH

  integer nspec,nglob

  integer nspec_iso,nspec_tiso,nspec_ani

  real(kind=CUSTOM_REAL), dimension(nglob) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec) :: &
    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,jacobian

! material properties
  real(kind=CUSTOM_REAL) rhostore(NGLLX,NGLLY,NGLLZ,nspec_iso)
  real(kind=CUSTOM_REAL) kappavstore(NGLLX,NGLLY,NGLLZ,nspec_iso)
  real(kind=CUSTOM_REAL) muvstore(NGLLX,NGLLY,NGLLZ,nspec_iso)

! additional arrays for anisotropy stored only where needed to save memory
  real(kind=CUSTOM_REAL) kappahstore(NGLLX,NGLLY,NGLLZ,nspec_tiso)
  real(kind=CUSTOM_REAL) muhstore(NGLLX,NGLLY,NGLLZ,nspec_tiso)
  real(kind=CUSTOM_REAL) eta_anisostore(NGLLX,NGLLY,NGLLZ,nspec_tiso)

! additional arrays for full anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,nspec_ani) :: &
    c11store,c12store,c13store,c14store,c15store,c16store, &
    c22store,c23store,c24store,c25store,c26store,c33store,c34store, &
    c35store,c36store,c44store,c45store,c46store,c55store,c56store,c66store

! Stacey
  real(kind=CUSTOM_REAL) rho_vp(NGLLX,NGLLY,NGLLZ,nspec)
  real(kind=CUSTOM_REAL) rho_vs(NGLLX,NGLLY,NGLLZ,nspec)

! mass matrix and additional ocean load mass matrix
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass,rmass_ocean_load

! global addressing
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

! processor identification
  character(len=150) prname

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

! xix
  open(unit=IIN,file=prname(1:len_trim(prname))//'solver_data_1.bin',status='old',action='read',form='unformatted')
  read(IIN) xix
!   close(IIN)

! xiy
!   open(unit=IIN,file=prname(1:len_trim(prname))//'xiy.bin',status='old',action='read',form='unformatted')
  read(IIN) xiy
!   close(IIN)

! xiz
!   open(unit=IIN,file=prname(1:len_trim(prname))//'xiz.bin',status='old',action='read',form='unformatted')
  read(IIN) xiz
!   close(IIN)

! etax
!   open(unit=IIN,file=prname(1:len_trim(prname))//'etax.bin',status='old',action='read',form='unformatted')
  read(IIN) etax
!   close(IIN)

! etay
!   open(unit=IIN,file=prname(1:len_trim(prname))//'etay.bin',status='old',action='read',form='unformatted')
  read(IIN) etay
!   close(IIN)

! etaz
!   open(unit=IIN,file=prname(1:len_trim(prname))//'etaz.bin',status='old',action='read',form='unformatted')
  read(IIN) etaz
!   close(IIN)

! gammax
!   open(unit=IIN,file=prname(1:len_trim(prname))//'gammax.bin',status='old',action='read',form='unformatted')
  read(IIN) gammax
!   close(IIN)

! gammay
!   open(unit=IIN,file=prname(1:len_trim(prname))//'gammay.bin',status='old',action='read',form='unformatted')
  read(IIN) gammay
!   close(IIN)

! gammaz
!   open(unit=IIN,file=prname(1:len_trim(prname))//'gammaz.bin',status='old',action='read',form='unformatted')
  read(IIN) gammaz
!   close(IIN)

! jacobian
!   open(unit=IIN,file=prname(1:len_trim(prname))//'jacobian.bin',status='old',action='read',form='unformatted')
  read(IIN) jacobian
!   close(IIN)

! model arrays

!   rho
!     open(unit=IIN,file=prname(1:len_trim(prname))//'rho.bin',status='old',action='read',form='unformatted')
    read(IIN) rhostore
!     close(IIN)

!   kappav
!     open(unit=IIN,file=prname(1:len_trim(prname))//'kappav.bin',status='old',action='read',form='unformatted')
    read(IIN) kappavstore
!     close(IIN)

  if(READ_KAPPA_MU) then

!   muv
!     open(unit=IIN,file=prname(1:len_trim(prname))//'muv.bin',status='old',action='read',form='unformatted')
    read(IIN) muvstore
!     close(IIN)

  endif


! for anisotropy, gravity and rotation

  if(TRANSVERSE_ISOTROPY .and. READ_TISO) then

!   kappah
!     open(unit=IIN,file=prname(1:len_trim(prname))//'kappah.bin',status='old',action='read',form='unformatted')
    read(IIN) kappahstore
!     close(IIN)

!   muh
!     open(unit=IIN,file=prname(1:len_trim(prname))//'muh.bin',status='old',action='read',form='unformatted')
    read(IIN) muhstore
!     close(IIN)

!   eta_aniso
!     open(unit=IIN,file=prname(1:len_trim(prname))//'eta_aniso.bin',status='old',action='read',form='unformatted')
    read(IIN) eta_anisostore
!     close(IIN)

  endif

  if(ANISOTROPIC_INNER_CORE .and. iregion_code == IREGION_INNER_CORE) then

!   model arrays

!   c11
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c11_inner_core.bin',status='old',action='read',form='unformatted')
    read(IIN) c11store
!     close(IIN)

!   c12
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c12_inner_core.bin',status='old',action='read',form='unformatted')
    read(IIN) c12store
!     close(IIN)

!   c13
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c13_inner_core.bin',status='old',action='read',form='unformatted')
    read(IIN) c13store
!     close(IIN)

!   c33
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c33_inner_core.bin',status='old',action='read',form='unformatted')
    read(IIN) c33store
!     close(IIN)

!   c44
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c44_inner_core.bin',status='old',action='read',form='unformatted')
    read(IIN) c44store
!     close(IIN)

  endif

  if(ANISOTROPIC_3D_MANTLE .and. iregion_code == IREGION_CRUST_MANTLE) then

!   model arrays

!   c11
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c11_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c11store
!     close(IIN)

!   c12
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c12_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c12store
!     close(IIN)

!   c13
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c13_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c13store
!     close(IIN)

!   c14
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c14_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c14store
!     close(IIN)

!   c15
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c15_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c15store
!     close(IIN)

!   c16
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c16_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c16store
!     close(IIN)

!   c22
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c22_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c22store
!     close(IIN)

!   c23
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c23_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c23store
!     close(IIN)

!   c24
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c24_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c24store
!     close(IIN)

!   c25
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c25_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c25store
!     close(IIN)

!   c26
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c26_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c26store
!     close(IIN)

!   c33
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c33_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c33store
!     close(IIN)

!   c34
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c34_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c34store
!     close(IIN)

!   c35
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c35_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c35store
!     close(IIN)

!   c36
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c36_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c36store
!     close(IIN)

!   c44
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c44_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c44store
!     close(IIN)

!   c45
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c45_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c45store
!     close(IIN)

!   c46
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c46_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c46store
!     close(IIN)

!   c55
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c55_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c55store
!     close(IIN)

!   c56
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c56_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c56store
!     close(IIN)

!   c66
!     open(unit=IIN,file=prname(1:len_trim(prname))//'c66_mantle.bin',status='old',action='read',form='unformatted')
    read(IIN) c66store
!     close(IIN)
  endif

! Stacey
  if(NCHUNKS /= 6) then

    if(iregion_code == IREGION_CRUST_MANTLE) then

! rho_vp
!       open(unit=IIN,file=prname(1:len_trim(prname))//'rho_vp_mantle.bin',status='old',action='read',form='unformatted')
      read(IIN) rho_vp
!       close(IIN)

! rho_vs
!       open(unit=IIN,file=prname(1:len_trim(prname))//'rho_vs_mantle.bin',status='old',action='read',form='unformatted')
      read(IIN) rho_vs
!       close(IIN)

    elseif(iregion_code == IREGION_OUTER_CORE) then

! vp
!       open(unit=IIN,file=prname(1:len_trim(prname))//'vp_outer_core.bin',status='old',action='read',form='unformatted')
      read(IIN) rho_vp
!       close(IIN)

    endif

  endif

! ibool
!   open(unit=IIN,file=prname(1:len_trim(prname))//'ibool.bin',status='old',action='read',form='unformatted')
  read(IIN) ibool
!   close(IIN)

! idoubling
!   open(unit=IIN,file=prname(1:len_trim(prname))//'idoubling.bin',status='old',action='read',form='unformatted')
  read(IIN) idoubling
!   close(IIN)

! mass matrix
!   open(unit=IIN,file=prname(1:len_trim(prname))//'rmass.bin',status='old',action='read',form='unformatted')
  read(IIN) rmass
!   close(IIN)

! read additional ocean load mass matrix
  if(OCEANS .and. iregion_code == IREGION_CRUST_MANTLE) then
!     open(unit=IIN,file=prname(1:len_trim(prname))//'rmass_ocean_load.bin',status='old',action='read',form='unformatted')
    read(IIN) rmass_ocean_load
!     close(IIN)
  endif

! read coordinates of the mesh
!   open(unit=IIN,file=prname(1:len_trim(prname))//'x.bin',status='old',action='read',form='unformatted')
  read(IIN) xstore
!   close(IIN)

!   open(unit=IIN,file=prname(1:len_trim(prname))//'y.bin',status='old',action='read',form='unformatted')
  read(IIN) ystore
!   close(IIN)

!   open(unit=IIN,file=prname(1:len_trim(prname))//'z.bin',status='old',action='read',form='unformatted')
  read(IIN) zstore
  close(IIN)

  end subroutine read_arrays_solver

