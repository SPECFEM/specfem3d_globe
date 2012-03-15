!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine read_mesh_databases(myrank,rho_vp_crust_mantle,rho_vs_crust_mantle, &
            xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
            kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
            nspec_iso,nspec_tiso,nspec_ani, &
            c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
            c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
            c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
            c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
            c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
            c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
            c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
            ibool_crust_mantle,ispec_is_tiso_crust_mantle, &
         ! -- idoubling_crust_mantle
            is_on_a_slice_edge_crust_mantle,rmass_crust_mantle,rmass_ocean_load, &
            vp_outer_core,xstore_outer_core,ystore_outer_core,zstore_outer_core, &
            xix_outer_core,xiy_outer_core,xiz_outer_core, &
            etax_outer_core,etay_outer_core,etaz_outer_core, &
            gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
            rhostore_outer_core,kappavstore_outer_core, &
            ibool_outer_core,idoubling_outer_core,ispec_is_tiso_outer_core, &
            is_on_a_slice_edge_outer_core,rmass_outer_core, &
            xstore_inner_core,ystore_inner_core,zstore_inner_core, &
            xix_inner_core,xiy_inner_core,xiz_inner_core, &
            etax_inner_core,etay_inner_core,etaz_inner_core, &
            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
            c11store_inner_core,c12store_inner_core,c13store_inner_core, &
            c33store_inner_core,c44store_inner_core, &
            ibool_inner_core,idoubling_inner_core,ispec_is_tiso_inner_core, &
            is_on_a_slice_edge_inner_core,rmass_inner_core, &
            ABSORBING_CONDITIONS,LOCAL_PATH)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! Stacey
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY) :: &
    rho_vp_crust_mantle,rho_vs_crust_mantle

  ! mesh parameters
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
        etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
        gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle
  ! arrays for anisotropic elements stored only where needed to save space
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle

  ! arrays for full anisotropy only when needed
  integer nspec_iso,nspec_tiso,nspec_ani
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
        c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
        c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
        c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
        c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
        c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
        c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
        c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

!  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling_crust_mantle
  logical, dimension(NSPEC_CRUST_MANTLE) :: ispec_is_tiso_crust_mantle

  ! mass matrix
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmass_crust_mantle
  ! additional mass matrix for ocean load
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load

  ! stacy outer core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_STACEY) :: vp_outer_core
  ! mesh parameters
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
        xstore_outer_core,ystore_outer_core,zstore_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        xix_outer_core,xiy_outer_core,xiz_outer_core,&
        etax_outer_core,etay_outer_core,etaz_outer_core, &
        gammax_outer_core,gammay_outer_core,gammaz_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: &
        rhostore_outer_core,kappavstore_outer_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core
  integer, dimension(NSPEC_OUTER_CORE) :: idoubling_outer_core
  logical, dimension(NSPEC_OUTER_CORE) :: ispec_is_tiso_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: rmass_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
        xix_inner_core,xiy_inner_core,xiz_inner_core,&
        etax_inner_core,etay_inner_core,etaz_inner_core, &
        gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
        rhostore_inner_core, kappavstore_inner_core,muvstore_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: &
        xstore_inner_core,ystore_inner_core,zstore_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC) :: &
        c11store_inner_core,c33store_inner_core,c12store_inner_core, &
        c13store_inner_core,c44store_inner_core
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core
  logical, dimension(NSPEC_INNER_CORE) :: ispec_is_tiso_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: rmass_inner_core

  logical ABSORBING_CONDITIONS
  character(len=150) LOCAL_PATH

  !local parameters
  logical READ_KAPPA_MU,READ_TISO
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,1) :: dummy_array
  integer, dimension(NSPEC_CRUST_MANTLE) :: dummy_i

! this for non blocking MPI
  logical, dimension(NSPEC_CRUST_MANTLE) :: is_on_a_slice_edge_crust_mantle
  logical, dimension(NSPEC_OUTER_CORE) :: is_on_a_slice_edge_outer_core
  logical, dimension(NSPEC_INNER_CORE) :: is_on_a_slice_edge_inner_core

  ! start reading the databases
  ! read arrays created by the mesher

  ! crust and mantle
  if(ANISOTROPIC_3D_MANTLE_VAL) then
    READ_KAPPA_MU = .false.
    READ_TISO = .false.
    nspec_iso = 1
    nspec_tiso = 1
    nspec_ani = NSPEC_CRUST_MANTLE
  else
    nspec_iso = NSPEC_CRUST_MANTLE
    if(TRANSVERSE_ISOTROPY_VAL) then
      nspec_tiso = NSPECMAX_TISO_MANTLE
    else
      nspec_tiso = 1
    endif
    nspec_ani = 1
    READ_KAPPA_MU = .true.
    READ_TISO = .true.
  endif
  call read_arrays_solver(IREGION_CRUST_MANTLE,myrank, &
            rho_vp_crust_mantle,rho_vs_crust_mantle, &
            xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
            xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
            etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
            gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
            rhostore_crust_mantle,kappavstore_crust_mantle,muvstore_crust_mantle, &
            kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle, &
            nspec_iso,nspec_tiso,nspec_ani, &
            c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
            c14store_crust_mantle,c15store_crust_mantle,c16store_crust_mantle, &
            c22store_crust_mantle,c23store_crust_mantle,c24store_crust_mantle, &
            c25store_crust_mantle,c26store_crust_mantle,c33store_crust_mantle, &
            c34store_crust_mantle,c35store_crust_mantle,c36store_crust_mantle, &
            c44store_crust_mantle,c45store_crust_mantle,c46store_crust_mantle, &
            c55store_crust_mantle,c56store_crust_mantle,c66store_crust_mantle, &
            ibool_crust_mantle,dummy_i, &
          ! --idoubling_crust_mantle,
            ispec_is_tiso_crust_mantle, &
            is_on_a_slice_edge_crust_mantle,rmass_crust_mantle,rmass_ocean_load, &
            NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  ! outer core (no anisotropy nor S velocity)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the outer core
  READ_KAPPA_MU = .false.
  READ_TISO = .false.
  nspec_iso = NSPEC_OUTER_CORE
  nspec_tiso = 1
  nspec_ani = 1

  call read_arrays_solver(IREGION_OUTER_CORE,myrank, &
            vp_outer_core,dummy_array, &
            xstore_outer_core,ystore_outer_core,zstore_outer_core, &
            xix_outer_core,xiy_outer_core,xiz_outer_core, &
            etax_outer_core,etay_outer_core,etaz_outer_core, &
            gammax_outer_core,gammay_outer_core,gammaz_outer_core, &
            rhostore_outer_core,kappavstore_outer_core,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            nspec_iso,nspec_tiso,nspec_ani, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            ibool_outer_core,idoubling_outer_core,ispec_is_tiso_outer_core, &
            is_on_a_slice_edge_outer_core,rmass_outer_core,rmass_ocean_load, &
            NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  ! inner core (no anisotropy)
  ! rmass_ocean_load is not used in this routine because it is meaningless in the inner core
  READ_KAPPA_MU = .true.
  READ_TISO = .false.
  nspec_iso = NSPEC_INNER_CORE
  nspec_tiso = 1
  if(ANISOTROPIC_INNER_CORE_VAL) then
    nspec_ani = NSPEC_INNER_CORE
  else
    nspec_ani = 1
  endif

  call read_arrays_solver(IREGION_INNER_CORE,myrank, &
            dummy_array,dummy_array, &
            xstore_inner_core,ystore_inner_core,zstore_inner_core, &
            xix_inner_core,xiy_inner_core,xiz_inner_core, &
            etax_inner_core,etay_inner_core,etaz_inner_core, &
            gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
            rhostore_inner_core,kappavstore_inner_core,muvstore_inner_core, &
            dummy_array,dummy_array,dummy_array, &
            nspec_iso,nspec_tiso,nspec_ani, &
            c11store_inner_core,c12store_inner_core,c13store_inner_core, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            dummy_array,dummy_array,c33store_inner_core, &
            dummy_array,dummy_array,dummy_array, &
            c44store_inner_core,dummy_array,dummy_array, &
            dummy_array,dummy_array,dummy_array, &
            ibool_inner_core,idoubling_inner_core,ispec_is_tiso_inner_core, &
            is_on_a_slice_edge_inner_core,rmass_inner_core,rmass_ocean_load, &
            NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
            READ_KAPPA_MU,READ_TISO,TRANSVERSE_ISOTROPY_VAL,ANISOTROPIC_3D_MANTLE_VAL, &
            ANISOTROPIC_INNER_CORE_VAL,OCEANS_VAL,LOCAL_PATH,ABSORBING_CONDITIONS)

  ! check that the number of points in this slice is correct
  if(minval(ibool_crust_mantle(:,:,:,:)) /= 1 .or. &
    maxval(ibool_crust_mantle(:,:,:,:)) /= NGLOB_CRUST_MANTLE) &
      call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in crust and mantle')

  if(minval(ibool_outer_core(:,:,:,:)) /= 1 .or. &
     maxval(ibool_outer_core(:,:,:,:)) /= NGLOB_OUTER_CORE) &
    call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in outer core')

  if(minval(ibool_inner_core(:,:,:,:)) /= 1 .or. maxval(ibool_inner_core(:,:,:,:)) /= NGLOB_INNER_CORE) &
    call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in inner core')

  end subroutine read_mesh_databases

!
!-------------------------------------------------------------------------------------------------
!


  subroutine read_mesh_databases_addressing(myrank, &
                    iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle, &
                    iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
                    npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
                    iboolfaces_crust_mantle,npoin2D_faces_crust_mantle, &
                    iboolcorner_crust_mantle, &
                    iboolleft_xi_outer_core,iboolright_xi_outer_core, &
                    iboolleft_eta_outer_core,iboolright_eta_outer_core, &
                    npoin2D_xi_outer_core,npoin2D_eta_outer_core,&
                    iboolfaces_outer_core,npoin2D_faces_outer_core, &
                    iboolcorner_outer_core, &
                    iboolleft_xi_inner_core,iboolright_xi_inner_core, &
                    iboolleft_eta_inner_core,iboolright_eta_inner_core, &
                    npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
                    iboolfaces_inner_core,npoin2D_faces_inner_core, &
                    iboolcorner_inner_core, &
                    iprocfrom_faces,iprocto_faces,imsg_type, &
                    iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
                    LOCAL_PATH,OUTPUT_FILES, &
                    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB1D_RADIAL, &
                    NGLOB2DMAX_XY,NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
                    addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
                    ichunk,iproc_xi,iproc_eta)

  implicit none

  include 'mpif.h'
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! 2-D addressing and buffers for summation between slices
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_OC) :: iboolleft_xi_outer_core,iboolright_xi_outer_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_OC) :: iboolleft_eta_outer_core,iboolright_eta_outer_core
  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_outer_core,npoin2D_eta_outer_core
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_inner_core,npoin2D_eta_inner_core

  integer, dimension(NGLOB2DMAX_XY_VAL,NUMFACES_SHARED) :: iboolfaces_crust_mantle, &
      iboolfaces_outer_core,iboolfaces_inner_core

  integer npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer npoin2D_faces_outer_core(NUMFACES_SHARED)
  integer npoin2D_faces_inner_core(NUMFACES_SHARED)

  ! indirect addressing for each corner of the chunks
  integer, dimension(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_OC,NUMCORNERS_SHARED) :: iboolcorner_outer_core
  integer, dimension(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED) :: iboolcorner_inner_core

  ! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: iprocfrom_faces,iprocto_faces,imsg_type
  ! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  character(len=150) LOCAL_PATH,OUTPUT_FILES

  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer NGLOB2DMAX_XY

  integer NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS

  ! for addressing of the slices
  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing
  integer, dimension(0:NPROCTOT_VAL-1) :: ichunk_slice,iproc_xi_slice,iproc_eta_slice
  integer ichunk,iproc_xi,iproc_eta

  ! local parameters
  integer :: ier,iproc,iproc_read
  integer :: NUM_FACES,NPROC_ONE_DIRECTION

  ! open file with global slice number addressing
  if(myrank == 0) then
    open(unit=IIN,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening addressing.txt')
    do iproc = 0,NPROCTOT_VAL-1
      read(IIN,*) iproc_read,ichunk,iproc_xi,iproc_eta
      if(iproc_read /= iproc) call exit_MPI(myrank,'incorrect slice number read')
      addressing(ichunk,iproc_xi,iproc_eta) = iproc
      ichunk_slice(iproc) = ichunk
      iproc_xi_slice(iproc) = iproc_xi
      iproc_eta_slice(iproc) = iproc_eta
    enddo
    close(IIN)
  endif

  ! broadcast the information read on the master to the nodes
  call MPI_BCAST(addressing,NCHUNKS_VAL*NPROC_XI_VAL*NPROC_ETA_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(ichunk_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iproc_xi_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iproc_eta_slice,NPROCTOT_VAL,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  ! output a topology map of slices
!!!!!!!  if (myrank == 0 .and. NCHUNKS_VAL == 6) then
!!!!!!! commented out because crashes when run on a very large machine
!!!!!!! because the records become too long
  if (.false.) then
    write(IMAIN,*) 'Spatial distribution of the slices'
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(20x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AB
        write(IMAIN,'(i10)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(1x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_BC
        write(IMAIN,'(i10)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(3x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AC
        write(IMAIN,'(i10)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(3x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_BC_ANTIPODE
        write(IMAIN,'(i10)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(20x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AB_ANTIPODE
        write(IMAIN,'(i10)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
    do iproc_xi = NPROC_XI_VAL-1, 0, -1
      write(IMAIN,'(20x)',advance='no')
      do iproc_eta = NPROC_ETA_VAL -1, 0, -1
        ichunk = CHUNK_AC_ANTIPODE
        write(IMAIN,'(i10)',advance='no') addressing(ichunk,iproc_xi,iproc_eta)
      enddo
      write(IMAIN,'(1x)',advance='yes')
    enddo
    write(IMAIN, *) ' '
  endif

  ! determine chunk number and local slice coordinates using addressing
  ichunk = ichunk_slice(myrank)
  iproc_xi = iproc_xi_slice(myrank)
  iproc_eta = iproc_eta_slice(myrank)

  ! define maximum size for message buffers
  ! use number of elements found in the mantle since it is the largest region
  NGLOB2DMAX_XY = max(NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE))

  ! number of corners and faces shared between chunks and number of message types
  if(NCHUNKS_VAL == 1 .or. NCHUNKS_VAL == 2) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 1
  else if(NCHUNKS_VAL == 3) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 3
  else if(NCHUNKS_VAL == 6) then
    NCORNERSCHUNKS = 8
    NUM_FACES = 4
    NUM_MSG_TYPES = 3
  else
    call exit_MPI(myrank,'number of chunks must be either 1, 2, 3 or 6')
  endif
  ! if more than one chunk then same number of processors in each direction
  NPROC_ONE_DIRECTION = NPROC_XI_VAL
  ! total number of messages corresponding to these common faces
  NUMMSGS_FACES = NPROC_ONE_DIRECTION*NUM_FACES*NUM_MSG_TYPES

  ! debug checks with compiled value
  !if( NUMMSGS_FACES /= NUMMSGS_FACES_VAL ) then
  !  print*,'check: NUMMSGS_FACES',NUMMSGS_FACES,NUMMSGS_FACES_VAL
  !  stop 'error NUMMSGS_FACES_VAL, please recompile solver'
  !endif

  ! read 2-D addressing for summation between slices with MPI

  ! mantle and crust
  call read_arrays_buffers_solver(IREGION_CRUST_MANTLE,myrank,iboolleft_xi_crust_mantle, &
     iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
     npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
     iboolfaces_crust_mantle,npoin2D_faces_crust_mantle, &
     iboolcorner_crust_mantle, &
     NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE), &
     NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_XY,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,LOCAL_PATH,NCHUNKS_VAL)

  ! outer core
  call read_arrays_buffers_solver(IREGION_OUTER_CORE,myrank, &
     iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
     npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
     iboolfaces_outer_core,npoin2D_faces_outer_core, &
     iboolcorner_outer_core, &
     NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE), &
     NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE),NGLOB2DMAX_XY,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,LOCAL_PATH,NCHUNKS_VAL)

  ! inner core
  call read_arrays_buffers_solver(IREGION_INNER_CORE,myrank, &
     iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
     npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
     iboolfaces_inner_core,npoin2D_faces_inner_core, &
     iboolcorner_inner_core, &
     NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE), &
     NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE),NGLOB2DMAX_XY,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT_VAL,NPROC_XI_VAL,NPROC_ETA_VAL,LOCAL_PATH,NCHUNKS_VAL)


  end subroutine read_mesh_databases_addressing


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_mesh_databases_coupling(myrank, &
              nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
              nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
              ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle,ibelm_ymin_crust_mantle, &
              ibelm_ymax_crust_mantle,ibelm_bottom_crust_mantle,ibelm_top_crust_mantle, &
              normal_xmin_crust_mantle,normal_xmax_crust_mantle,normal_ymin_crust_mantle, &
              normal_ymax_crust_mantle,normal_bottom_crust_mantle,normal_top_crust_mantle, &
              jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle,jacobian2D_ymin_crust_mantle, &
              jacobian2D_ymax_crust_mantle,jacobian2D_bottom_crust_mantle,jacobian2D_top_crust_mantle, &
              nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
              nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
              ibelm_xmin_outer_core,ibelm_xmax_outer_core,ibelm_ymin_outer_core, &
              ibelm_ymax_outer_core,ibelm_bottom_outer_core,ibelm_top_outer_core, &
              normal_xmin_outer_core,normal_xmax_outer_core,normal_ymin_outer_core, &
              normal_ymax_outer_core,normal_bottom_outer_core,normal_top_outer_core, &
              jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core,jacobian2D_ymin_outer_core, &
              jacobian2D_ymax_outer_core,jacobian2D_bottom_outer_core,jacobian2D_top_outer_core, &
              nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
              nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
              ibelm_xmin_inner_core,ibelm_xmax_inner_core,ibelm_ymin_inner_core, &
              ibelm_ymax_inner_core,ibelm_bottom_inner_core,ibelm_top_inner_core, &
              ibelm_moho_top,ibelm_moho_bot,ibelm_400_top,ibelm_400_bot, &
              ibelm_670_top,ibelm_670_bot,normal_moho,normal_400,normal_670, &
              k_top,k_bot,moho_kl,d400_kl,d670_kl,cmb_kl,icb_kl, &
              LOCAL_PATH,SIMULATION_TYPE)

! to couple mantle with outer core
  implicit none

  include 'mpif.h'
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! for crust/oceans coupling
  integer nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_CM) :: ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_CM) :: ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle
  integer, dimension(NSPEC2D_BOTTOM_CM) :: ibelm_bottom_crust_mantle
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_CM) :: jacobian2D_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_CM) :: jacobian2D_top_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    normal_xmin_crust_mantle,normal_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    normal_ymin_crust_mantle,normal_ymax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_CM) :: normal_bottom_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_CM) :: normal_top_crust_mantle

  ! arrays to couple with the fluid regions by pointwise matching
  integer nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_OC) :: ibelm_xmin_outer_core,ibelm_xmax_outer_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_OC) :: ibelm_ymin_outer_core,ibelm_ymax_outer_core
  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core
  integer, dimension(NSPEC2D_TOP_OC) :: ibelm_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_OC) :: jacobian2D_top_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: &
    jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: &
    jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: &
    normal_xmin_outer_core,normal_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: &
    normal_ymin_outer_core,normal_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: normal_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_TOP_OC) :: normal_top_outer_core

  ! inner core
  integer nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
    nspec2D_ymin_inner_core,nspec2D_ymax_inner_core
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_IC) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_IC) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_IC) :: ibelm_bottom_inner_core
  integer, dimension(NSPEC2D_TOP_IC) :: ibelm_top_inner_core

  ! boundary
  integer, dimension(NSPEC2D_MOHO) :: ibelm_moho_top,ibelm_moho_bot
  integer, dimension(NSPEC2D_400) :: ibelm_400_top,ibelm_400_bot
  integer, dimension(NSPEC2D_670) :: ibelm_670_top,ibelm_670_bot
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_MOHO) :: normal_moho
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_400) :: normal_400
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2D_670) :: normal_670

  integer k_top,k_bot

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_MOHO) :: moho_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_400) :: d400_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_670) ::  d670_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_CMB) :: cmb_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_ICB) :: icb_kl

  character(len=150) LOCAL_PATH
  integer SIMULATION_TYPE

  ! local parameters
  integer njunk1,njunk2,njunk3
  character(len=150) prname


  ! crust and mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! Stacey put back
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read')
  read(27) nspec2D_xmin_crust_mantle
  read(27) nspec2D_xmax_crust_mantle
  read(27) nspec2D_ymin_crust_mantle
  read(27) nspec2D_ymax_crust_mantle
  read(27) njunk1
  read(27) njunk2

! boundary parameters
  read(27) ibelm_xmin_crust_mantle
  read(27) ibelm_xmax_crust_mantle
  read(27) ibelm_ymin_crust_mantle
  read(27) ibelm_ymax_crust_mantle
  read(27) ibelm_bottom_crust_mantle
  read(27) ibelm_top_crust_mantle

  read(27) normal_xmin_crust_mantle
  read(27) normal_xmax_crust_mantle
  read(27) normal_ymin_crust_mantle
  read(27) normal_ymax_crust_mantle
  read(27) normal_bottom_crust_mantle
  read(27) normal_top_crust_mantle

  read(27) jacobian2D_xmin_crust_mantle
  read(27) jacobian2D_xmax_crust_mantle
  read(27) jacobian2D_ymin_crust_mantle
  read(27) jacobian2D_ymax_crust_mantle
  read(27) jacobian2D_bottom_crust_mantle
  read(27) jacobian2D_top_crust_mantle
  close(27)


  ! read parameters to couple fluid and solid regions
  !
  ! outer core

  ! create name of database
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  ! boundary parameters

  ! Stacey put back
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read')
  read(27) nspec2D_xmin_outer_core
  read(27) nspec2D_xmax_outer_core
  read(27) nspec2D_ymin_outer_core
  read(27) nspec2D_ymax_outer_core
  read(27) njunk1
  read(27) njunk2

  read(27) ibelm_xmin_outer_core
  read(27) ibelm_xmax_outer_core
  read(27) ibelm_ymin_outer_core
  read(27) ibelm_ymax_outer_core
  read(27) ibelm_bottom_outer_core
  read(27) ibelm_top_outer_core

  read(27) normal_xmin_outer_core
  read(27) normal_xmax_outer_core
  read(27) normal_ymin_outer_core
  read(27) normal_ymax_outer_core
  read(27) normal_bottom_outer_core
  read(27) normal_top_outer_core

  read(27) jacobian2D_xmin_outer_core
  read(27) jacobian2D_xmax_outer_core
  read(27) jacobian2D_ymin_outer_core
  read(27) jacobian2D_ymax_outer_core
  read(27) jacobian2D_bottom_outer_core
  read(27) jacobian2D_top_outer_core
  close(27)


  !
  ! inner core
  !

  ! create name of database
  call create_name_database(prname,myrank,IREGION_INNER_CORE,LOCAL_PATH)

  ! read info for vertical edges for central cube matching in inner core
  open(unit=27,file=prname(1:len_trim(prname))//'boundary.bin', &
        status='old',form='unformatted',action='read')
  read(27) nspec2D_xmin_inner_core
  read(27) nspec2D_xmax_inner_core
  read(27) nspec2D_ymin_inner_core
  read(27) nspec2D_ymax_inner_core
  read(27) njunk1
  read(27) njunk2

  ! boundary parameters
  read(27) ibelm_xmin_inner_core
  read(27) ibelm_xmax_inner_core
  read(27) ibelm_ymin_inner_core
  read(27) ibelm_ymax_inner_core
  read(27) ibelm_bottom_inner_core
  read(27) ibelm_top_inner_core
  close(27)


  ! -- Boundary Mesh for crust and mantle ---
  if (SAVE_BOUNDARY_MESH .and. SIMULATION_TYPE == 3) then

    call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

    open(unit=27,file=prname(1:len_trim(prname))//'boundary_disc.bin', &
          status='old',form='unformatted',action='read')
    read(27) njunk1,njunk2,njunk3
    if (njunk1 /= NSPEC2D_MOHO .and. njunk2 /= NSPEC2D_400 .and. njunk3 /= NSPEC2D_670) &
               call exit_mpi(myrank, 'Error reading ibelm_disc.bin file')
    read(27) ibelm_moho_top
    read(27) ibelm_moho_bot
    read(27) ibelm_400_top
    read(27) ibelm_400_bot
    read(27) ibelm_670_top
    read(27) ibelm_670_bot
    read(27) normal_moho
    read(27) normal_400
    read(27) normal_670
    close(27)

    k_top = 1
    k_bot = NGLLZ

    ! initialization
    moho_kl = 0.; d400_kl = 0.; d670_kl = 0.; cmb_kl = 0.; icb_kl = 0.
  endif

  end subroutine read_mesh_databases_coupling


!
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------
!


  subroutine read_mesh_databases_stacey(myrank, &
                      nimin_crust_mantle,nimax_crust_mantle,njmin_crust_mantle, &
                      njmax_crust_mantle,nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
                      nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
                      nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
                      reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
                      reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
                      nimin_outer_core,nimax_outer_core,njmin_outer_core, &
                      njmax_outer_core,nkmin_xi_outer_core,nkmin_eta_outer_core, &
                      nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                      nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                      reclen_xmin_outer_core,reclen_xmax_outer_core, &
                      reclen_ymin_outer_core,reclen_ymax_outer_core, &
                      reclen_zmin,NSPEC2D_BOTTOM, &
                      SIMULATION_TYPE,SAVE_FORWARD,LOCAL_PATH,NSTEP)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_CM) ::  &
    nimin_crust_mantle,nimax_crust_mantle,nkmin_eta_crust_mantle
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_CM) ::  &
    njmin_crust_mantle,njmax_crust_mantle,nkmin_xi_crust_mantle
  integer nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle
  integer reclen_xmin_crust_mantle, reclen_xmax_crust_mantle, reclen_ymin_crust_mantle, &
    reclen_ymax_crust_mantle

  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_OC) :: nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_OC) :: njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core
  integer nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core
  integer reclen_xmin_outer_core, reclen_xmax_outer_core,reclen_ymin_outer_core, &
    reclen_ymax_outer_core

  integer reclen_zmin
  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_BOTTOM

  integer SIMULATION_TYPE
  logical SAVE_FORWARD
  character(len=150) LOCAL_PATH
  integer NSTEP

  ! local parameters
  integer(kind=8) :: filesize
  character(len=150) prname


  ! crust and mantle
  ! create name of database
  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  ! read arrays for Stacey conditions
  open(unit=27,file=prname(1:len_trim(prname))//'stacey.bin', &
        status='old',form='unformatted',action='read')
  read(27) nimin_crust_mantle
  read(27) nimax_crust_mantle
  read(27) njmin_crust_mantle
  read(27) njmax_crust_mantle
  read(27) nkmin_xi_crust_mantle
  read(27) nkmin_eta_crust_mantle
  close(27)

  if (nspec2D_xmin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmin_crust_mantle)

    ! total file size
    filesize = reclen_xmin_crust_mantle
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=51,file=trim(prname)//'absorb_xmin.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmin_crust_mantle+2*4)
!    else
!      open(unit=51,file=trim(prname)//'absorb_xmin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_xmin_crust_mantle+2*4)

      call open_file_abs_r(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(0,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_xmin.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_xmax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLY * NGLLZ * nspec2D_xmax_crust_mantle)

    ! total file size
    filesize = reclen_xmax_crust_mantle
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=52,file=trim(prname)//'absorb_xmax.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmax_crust_mantle+2*4)
!    else
!      open(unit=52,file=trim(prname)//'absorb_xmax.bin', &
!            status='unknown',form='unformatted',access='direct', &
!            recl=reclen_xmax_crust_mantle+2*4)

      call open_file_abs_r(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(1,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_ymin_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymin_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymin_crust_mantle)

    ! total file size
    filesize = reclen_ymin_crust_mantle
    filesize = filesize*NSTEP


    if (SIMULATION_TYPE == 3) then
!      open(unit=53,file=trim(prname)//'absorb_ymin.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymin_crust_mantle+2*4)
!    else
!      open(unit=53,file=trim(prname)//'absorb_ymin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymin_crust_mantle+2*4)

      call open_file_abs_r(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(2,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_ymax_crust_mantle > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymax_crust_mantle = CUSTOM_REAL * (NDIM * NGLLX * NGLLZ * nspec2D_ymax_crust_mantle)

    ! total file size
    filesize = reclen_ymax_crust_mantle
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=54,file=trim(prname)//'absorb_ymax.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymax_crust_mantle+2*4)
!    else
!      open(unit=54,file=trim(prname)//'absorb_ymax.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymax_crust_mantle+2*4)

      call open_file_abs_r(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(3,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif


  ! outer core
  ! create name of database
  call create_name_database(prname,myrank,IREGION_OUTER_CORE,LOCAL_PATH)

  ! read arrays for Stacey conditions
  open(unit=27,file=prname(1:len_trim(prname))//'stacey.bin', &
        status='old',form='unformatted',action='read')
  read(27) nimin_outer_core
  read(27) nimax_outer_core
  read(27) njmin_outer_core
  read(27) njmax_outer_core
  read(27) nkmin_xi_outer_core
  read(27) nkmin_eta_outer_core
  close(27)

  if (nspec2D_xmin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmin_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmin_outer_core)

    ! total file size
    filesize = reclen_xmin_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=61,file=trim(prname)//'absorb_xmin.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmin_outer_core+2*4)
!    else
!      open(unit=61,file=trim(prname)//'absorb_xmin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_xmin_outer_core+2*4)

      call open_file_abs_r(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(4,trim(prname)//'absorb_xmin.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    endif
  endif

  if (nspec2D_xmax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_xmax_outer_core = CUSTOM_REAL * (NGLLY * NGLLZ * nspec2D_xmax_outer_core)

    ! total file size
    filesize = reclen_xmax_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=62,file=trim(prname)//'absorb_xmax.bin', &
!            status='old',action='read',form='unformatted',access='direct', &
!            recl=reclen_xmax_outer_core+2*4)
!    else
!      open(unit=62,file=trim(prname)//'absorb_xmax.bin', &
!            status='unknown',form='unformatted',access='direct', &
!            recl=reclen_xmax_outer_core+2*4)

      call open_file_abs_r(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
    else
      call open_file_abs_w(5,trim(prname)//'absorb_xmax.bin',len_trim(trim(prname)//'absorb_xmax.bin'), &
                          filesize)
   endif

  endif

  if (nspec2D_ymin_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymin_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymin_outer_core)

    ! total file size
    filesize = reclen_ymin_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=63,file=trim(prname)//'absorb_ymin.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymin_outer_core+2*4)
!    else
!      open(unit=63,file=trim(prname)//'absorb_ymin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymin_outer_core+2*4)

      call open_file_abs_r(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)
    else
      call open_file_abs_w(6,trim(prname)//'absorb_ymin.bin',len_trim(trim(prname)//'absorb_ymin.bin'), &
                          filesize)

    endif
  endif

  if (nspec2D_ymax_outer_core > 0 .and. (SIMULATION_TYPE == 3 &
    .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))) then

    ! size of single record
    reclen_ymax_outer_core = CUSTOM_REAL * (NGLLX * NGLLZ * nspec2D_ymax_outer_core)

    ! total file size
    filesize = reclen_ymax_outer_core
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=64,file=trim(prname)//'absorb_ymax.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_ymax_outer_core+2*4)
!    else
!      open(unit=64,file=trim(prname)//'absorb_ymax.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_ymax_outer_core+2*4)

      call open_file_abs_r(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)
    else
      call open_file_abs_w(7,trim(prname)//'absorb_ymax.bin',len_trim(trim(prname)//'absorb_ymax.bin'), &
                          filesize)

    endif
  endif

  if (NSPEC2D_BOTTOM(IREGION_OUTER_CORE) > 0 .and. &
     (SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)))then

    ! size of single record
    reclen_zmin = CUSTOM_REAL * (NGLLX * NGLLY * NSPEC2D_BOTTOM(IREGION_OUTER_CORE))

    ! total file size
    filesize = reclen_zmin
    filesize = filesize*NSTEP

    if (SIMULATION_TYPE == 3) then
!      open(unit=65,file=trim(prname)//'absorb_zmin.bin', &
!            status='old',action='read',form='unformatted',access='direct',&
!            recl=reclen_zmin+2*4)
!    else
!      open(unit=65,file=trim(prname)//'absorb_zmin.bin', &
!            status='unknown',form='unformatted',access='direct',&
!            recl=reclen_zmin+2*4)

      call open_file_abs_r(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    else
      call open_file_abs_w(8,trim(prname)//'absorb_zmin.bin',len_trim(trim(prname)//'absorb_zmin.bin'), &
                          filesize)
    endif
  endif

  end subroutine read_mesh_databases_stacey
