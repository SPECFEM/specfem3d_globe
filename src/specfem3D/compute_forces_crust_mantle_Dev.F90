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

  subroutine compute_forces_crust_mantle_Dev(minus_gravity_table,density_table,minus_deriv_gravity_table, &
          displ_crust_mantle,accel_crust_mantle,xstore,ystore,zstore, &
          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
            is_on_a_slice_edge_crust_mantle,icall, &
            accel_inner_core,ibool_inner_core,idoubling_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector,iphase, &
            nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
            npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
            receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_INNER_CORE,INCLUDE_CENTRAL_CUBE,iphase_CC, &
          hprime_xx,hprime_xxT, &
          hprimewgll_xx,hprimewgll_xxT, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube, &
          kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
          c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
          c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
          c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
          ibool,ispec_is_tiso, &
!ZN          R_memory,epsilondev,epsilon_trace_over_3,one_minus_sum_beta, &
          R_memory,one_minus_sum_beta,deltat,veloc_crust_mantle, &
          alphaval,betaval,gammaval,factor_common,vx,vy,vz,vnspec)

! this routine is optimized for NGLLX = NGLLY = NGLLZ = 5 using the Deville et al. (2002) inlined matrix-matrix products

  implicit none

  include "constants.h"

  ! include values created by the mesher
  ! done for performance only using static allocation to allow for loop unrolling
  include "OUTPUT_FILES/values_from_mesher.h"

  ! displacement and acceleration
  real(kind=CUSTOM_REAL) deltat  
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle,accel_crust_mantle,veloc_crust_mantle
  ! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool

  ! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore,ystore,zstore

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz

  ! array with derivatives of Lagrange polynomials and precalculated products
  double precision, dimension(NGLLX,NGLLY,NGLLZ) :: wgll_cube
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xxT,hprimewgll_xxT
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  ! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        kappahstore,muhstore,eta_anisostore
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        kappavstore,muvstore

  ! arrays for full anisotropy only when needed
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
        c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
        c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
        c36store,c44store,c45store,c46store,c55store,c56store,c66store

  ! attenuation
  ! memory variables for attenuation
  ! memory variables R_ij are stored at the local rather than global level
  ! to allow for optimization of cache access by compiler
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory
!ZN  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev
!ZN  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilon_trace_over_3

  integer :: vx,vy,vz,vnspec

  ! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(N_SLS, vx, vy, vz, vnspec) :: factor_common
  real(kind=CUSTOM_REAL), dimension(vx, vy, vz, vnspec) :: one_minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval

  ! array with the local to global mapping per slice
  logical, dimension(NSPEC_CRUST_MANTLE) :: ispec_is_tiso

  ! gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table,density_table,minus_deriv_gravity_table

! local parameters
  ! Deville
  ! manually inline the calls to the Deville et al. (2002) routines
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: dummyx_loc,dummyy_loc,dummyz_loc, &
    newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3
  real(kind=CUSTOM_REAL), dimension(NGLLX,m2) :: B1_m1_m2_5points,B2_m1_m2_5points,B3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: C1_m1_m2_5points,C2_m1_m2_5points,C3_m1_m2_5points
  real(kind=CUSTOM_REAL), dimension(m1,m2) :: E1_m1_m2_5points,E2_m1_m2_5points,E3_m1_m2_5points

  equivalence(dummyx_loc,B1_m1_m2_5points)
  equivalence(dummyy_loc,B2_m1_m2_5points)
  equivalence(dummyz_loc,B3_m1_m2_5points)
  equivalence(tempx1,C1_m1_m2_5points)
  equivalence(tempy1,C2_m1_m2_5points)
  equivalence(tempz1,C3_m1_m2_5points)
  equivalence(newtempx1,E1_m1_m2_5points)
  equivalence(newtempy1,E2_m1_m2_5points)
  equivalence(newtempz1,E3_m1_m2_5points)

  real(kind=CUSTOM_REAL), dimension(m2,NGLLX) :: &
    A1_mxm_m2_m1_5points,A2_mxm_m2_m1_5points,A3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    C1_mxm_m2_m1_5points,C2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points
  real(kind=CUSTOM_REAL), dimension(m2,m1) :: &
    E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points

  equivalence(dummyx_loc,A1_mxm_m2_m1_5points)
  equivalence(dummyy_loc,A2_mxm_m2_m1_5points)
  equivalence(dummyz_loc,A3_mxm_m2_m1_5points)
  equivalence(tempx3,C1_mxm_m2_m1_5points)
  equivalence(tempy3,C2_mxm_m2_m1_5points)
  equivalence(tempz3,C3_mxm_m2_m1_5points)
  equivalence(newtempx3,E1_mxm_m2_m1_5points)
  equivalence(newtempy3,E2_mxm_m2_m1_5points)
  equivalence(newtempz3,E3_mxm_m2_m1_5points)

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc_nplus1 !ZN
  real(kind=CUSTOM_REAL) fac1,fac2,fac3

  ! for gravity
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: rho_s_H

  integer :: ispec
  integer :: i,j,k
  integer :: iglob1


! this for non blocking MPI
  integer :: iphase,icall

  logical, dimension(NSPEC_CRUST_MANTLE) :: is_on_a_slice_edge_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: accel_inner_core

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core

  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core

  integer :: ichunk,iproc_xi,iproc_eta,myrank

  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer npoin2D_faces_inner_core(NUMFACES_SHARED)

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
       npoin2D_xi_inner_core,npoin2D_eta_inner_core

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: iprocfrom_faces,iprocto_faces

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  integer, dimension(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED) :: iboolcorner_inner_core

  integer, dimension(NGLOB2DMAX_XY_CM_VAL,NUMFACES_SHARED) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_IC_VAL,NUMFACES_SHARED) :: iboolfaces_inner_core

  integer :: npoin2D_max_all_CM_IC
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin2D_max_all_CM_IC) :: buffer_send_faces,buffer_received_faces

! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC) :: &
     buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector

! for matching with central cube in inner core
  integer nb_msgs_theor_in_cube, npoin2D_cube_from_slices,iphase_CC
  integer, dimension(nb_msgs_theor_in_cube) :: sender_from_slices_to_cube
  double precision, dimension(npoin2D_cube_from_slices,NDIM) :: buffer_slices
  double precision, dimension(npoin2D_cube_from_slices,NDIM,nb_msgs_theor_in_cube) :: buffer_all_cube_from_slices
  integer, dimension(nb_msgs_theor_in_cube,npoin2D_cube_from_slices):: ibool_central_cube
  integer receiver_cube_from_slices
  logical :: INCLUDE_CENTRAL_CUBE

! local to global mapping
  integer NSPEC2D_BOTTOM_INNER_CORE,iend,ispec_glob
  integer, dimension(NSPEC2D_BOTTOM_INNER_CORE) :: ibelm_bottom_inner_core

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

!$OMP PARALLEL DEFAULT(NONE) SHARED(xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
!$OMP one_minus_sum_beta,epsilon_trace_over_3,c11store,c12store,c13store,c14store,c15store, &
!$OMP c16store,c22store,c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
!$OMP c36store,c44store,c45store,c46store,c55store,c56store,c66store,ispec_is_tiso, &
!$OMP kappavstore,muvstore,kappahstore,muhstore,eta_anisostore,ibool,ystore,zstore, &
!$OMP R_memory,xstore,minus_gravity_table,minus_deriv_gravity_table,density_table, &
!$OMP displ_crust_mantle,wgll_cube,accel_inner_core,hprime_xxt,hprime_xx,idoubling_inner_core, &
!$OMP addressing,iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,vx,vy,vz,vnspec, &
!$OMP iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle,npoin2D_faces_crust_mantle, &
!$OMP npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle,iboolfaces_crust_mantle,iboolcorner_crust_mantle,iboolleft_xi_inner_core, &
!$OMP iboolright_xi_inner_core,ibool_inner_core, &
!$OMP iboolleft_eta_inner_core,iboolright_eta_inner_core,npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
!$OMP iboolfaces_inner_core,accel_crust_mantle, &
!$OMP iboolcorner_inner_core,iprocfrom_faces,iprocto_faces,iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
!$OMP buffer_send_faces,buffer_received_faces, &
!$OMP buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector, &
!$OMP sender_from_slices_to_cube,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
!$OMP ibelm_bottom_inner_core,hprimewgll_xx,hprimewgll_xxt,wgllwgll_xy,wgllwgll_xz, &
!$OMP wgllwgll_yz,INCLUDE_CENTRAL_CUBE,alphaval,betaval,epsilondev,gammaval,factor_common, &
!$OMP myrank,ichunk,iphase,iphase_CC,icall,iproc_xi,iproc_eta,is_on_a_slice_edge_crust_mantle, &
!$OMP npoin2D_cube_from_slices,nb_msgs_theor_in_cube,receiver_cube_from_slices, &
!$OMP npoin2D_max_all_CM_IC ) &
!$OMP PRIVATE(k,j,ispec,fac1,fac2,fac3,sum_terms,iend,ispec_glob, &
!$OMP C1_mxm_m2_m1_5points,A1_mxm_m2_m1_5points,B2_m1_m2_5points,C3_m1_m2_5points, &
!$OMP B3_m1_m2_5points,C2_mxm_m2_m1_5points,E1_m1_m2_5points,E2_m1_m2_5points,A2_mxm_m2_m1_5points,C3_mxm_m2_m1_5points, &
!$OMP A3_mxm_m2_m1_5points,B1_m1_m2_5points, &
!$OMP C2_m1_m2_5points,C1_m1_m2_5points,E3_m1_m2_5points,E1_mxm_m2_m1_5points,E2_mxm_m2_m1_5points,E3_mxm_m2_m1_5points, &
!$OMP tempx1,tempx2,tempx3, &
!$OMP newtempx1,newtempx2,newtempx3,newtempy1,newtempy2,newtempy3,newtempz1,newtempz2,newtempz3, &
!$OMP dummyx_loc,dummyy_loc,dummyz_loc,rho_s_H,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
!$OMP iglob1,epsilondev_loc)

  do ispec_glob = 1,NSPEC_CRUST_MANTLE,ELEMENTS_NONBLOCKING_CM_IC

! process the non-blocking communications every ELEMENTS_NONBLOCKING elements
    if (icall == 2) then

      if(iphase <= 7) then
!$OMP BARRIER
!$OMP MASTER
            call assemble_MPI_vector(myrank,accel_crust_mantle,accel_inner_core, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector, &
            NUMMSGS_FACES_VAL,NCORNERSCHUNKS_VAL, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL_CM, &
            NGLOB1D_RADIAL_IC,NCHUNKS_VAL,iphase)
!$OMP END MASTER
!$OMP BARRIER
      endif

      if(INCLUDE_CENTRAL_CUBE) then
          if(iphase > 7 .and. iphase_CC <= 4) then
!$OMP BARRIER
!$OMP MASTER
            call assemble_MPI_central_cube(ichunk,nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
                   npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
                   receiver_cube_from_slices,ibool_inner_core,idoubling_inner_core, &
                   ibelm_bottom_inner_core,NSPEC2D_BOTTOM_IC,accel_inner_core,NDIM,iphase_CC)
!$OMP END MASTER
!$OMP BARRIER
          endif
      endif

    endif

    iend = min(ispec_glob+ELEMENTS_NONBLOCKING_CM_IC-1,NSPEC_CRUST_MANTLE)

!$OMP DO SCHEDULE(runtime)
    do ispec=ispec_glob,iend

! hide communications by computing the edges first
    if((icall == 2 .and. is_on_a_slice_edge_crust_mantle(ispec)) .or. &
       (icall == 1 .and. .not. is_on_a_slice_edge_crust_mantle(ispec))) cycle

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
            iglob1 = ibool(i,j,k,ispec)
            dummyx_loc(i,j,k) = displ_crust_mantle(1,iglob1)
            dummyy_loc(i,j,k) = displ_crust_mantle(2,iglob1)
            dummyz_loc(i,j,k) = displ_crust_mantle(3,iglob1)
        enddo
      enddo
    enddo

    do j=1,m2
      do i=1,m1
        C1_m1_m2_5points(i,j) = hprime_xx(i,1)*B1_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B1_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B1_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B1_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B1_m1_m2_5points(5,j)

        C2_m1_m2_5points(i,j) = hprime_xx(i,1)*B2_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B2_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B2_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B2_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B2_m1_m2_5points(5,j)

        C3_m1_m2_5points(i,j) = hprime_xx(i,1)*B3_m1_m2_5points(1,j) + &
                              hprime_xx(i,2)*B3_m1_m2_5points(2,j) + &
                              hprime_xx(i,3)*B3_m1_m2_5points(3,j) + &
                              hprime_xx(i,4)*B3_m1_m2_5points(4,j) + &
                              hprime_xx(i,5)*B3_m1_m2_5points(5,j)
      enddo
    enddo

    do j=1,m1
      do i=1,m1
        ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
        do k = 1,NGLLX
          tempx2(i,j,k) = dummyx_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyx_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyx_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyx_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyx_loc(i,5,k)*hprime_xxT(5,j)

          tempy2(i,j,k) = dummyy_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyy_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyy_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyy_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyy_loc(i,5,k)*hprime_xxT(5,j)

          tempz2(i,j,k) = dummyz_loc(i,1,k)*hprime_xxT(1,j) + &
                        dummyz_loc(i,2,k)*hprime_xxT(2,j) + &
                        dummyz_loc(i,3,k)*hprime_xxT(3,j) + &
                        dummyz_loc(i,4,k)*hprime_xxT(4,j) + &
                        dummyz_loc(i,5,k)*hprime_xxT(5,j)
        enddo
      enddo
    enddo

    do j=1,m1
      do i=1,m2
        C1_mxm_m2_m1_5points(i,j) = A1_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A1_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A1_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A1_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A1_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

        C2_mxm_m2_m1_5points(i,j) = A2_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A2_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A2_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A2_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A2_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)

        C3_mxm_m2_m1_5points(i,j) = A3_mxm_m2_m1_5points(i,1)*hprime_xxT(1,j) + &
                                  A3_mxm_m2_m1_5points(i,2)*hprime_xxT(2,j) + &
                                  A3_mxm_m2_m1_5points(i,3)*hprime_xxT(3,j) + &
                                  A3_mxm_m2_m1_5points(i,4)*hprime_xxT(4,j) + &
                                  A3_mxm_m2_m1_5points(i,5)*hprime_xxT(5,j)
      enddo
    enddo

    !
    ! compute either isotropic, transverse isotropic or anisotropic elements
    !
    if(ANISOTROPIC_3D_MANTLE_VAL) then
      ! anisotropic element
      call compute_element_aniso(ispec, &
                    minus_gravity_table,density_table,minus_deriv_gravity_table, &
                    xstore,ystore,zstore, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    wgll_cube, &
                    c11store,c12store,c13store,c14store,c15store,c16store,c22store, &
                    c23store,c24store,c25store,c26store,c33store,c34store,c35store, &
                    c36store,c44store,c45store,c46store,c55store,c56store,c66store, &
                    ibool, &
!ZN                    R_memory,epsilon_trace_over_3, &
                    R_memory, &
                    one_minus_sum_beta,vx,vy,vz,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc,epsilondev_loc,rho_s_H)
    else
      if( .not. ispec_is_tiso(ispec) ) then
        ! isotropic element
        call compute_element_iso(ispec, &
                    minus_gravity_table,density_table,minus_deriv_gravity_table, &
                    xstore,ystore,zstore, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    wgll_cube, &
                    kappavstore,muvstore, &
                    ibool, &
!ZN                    R_memory,epsilon_trace_over_3, &
                    R_memory, & !ZN
                    one_minus_sum_beta,vx,vy,vz,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc,epsilondev_loc,rho_s_H)
      else
        ! transverse isotropic element
        call compute_element_tiso(ispec, &
                    minus_gravity_table,density_table,minus_deriv_gravity_table, &
                    xstore,ystore,zstore, &
                    xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz, &
                    wgll_cube, &
                    kappavstore,kappahstore,muvstore,muhstore,eta_anisostore, &
                    ibool, &
!ZN                    R_memory,epsilon_trace_over_3, &
                    R_memory, & !ZN
                    one_minus_sum_beta,vx,vy,vz,vnspec, &
                    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3, &
                    dummyx_loc,dummyy_loc,dummyz_loc,epsilondev_loc,rho_s_H)
      endif ! .not. ispec_is_tiso
    endif

    ! subroutines adapted from Deville, Fischer and Mund, High-order methods
    ! for incompressible fluid flow, Cambridge University Press (2002),
    ! pages 386 and 389 and Figure 8.3.1
    do j=1,m2
      do i=1,m1
        E1_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C1_m1_m2_5points(1,j) + &
                              hprimewgll_xxT(i,2)*C1_m1_m2_5points(2,j) + &
                              hprimewgll_xxT(i,3)*C1_m1_m2_5points(3,j) + &
                              hprimewgll_xxT(i,4)*C1_m1_m2_5points(4,j) + &
                              hprimewgll_xxT(i,5)*C1_m1_m2_5points(5,j)

        E2_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C2_m1_m2_5points(1,j) + &
                              hprimewgll_xxT(i,2)*C2_m1_m2_5points(2,j) + &
                              hprimewgll_xxT(i,3)*C2_m1_m2_5points(3,j) + &
                              hprimewgll_xxT(i,4)*C2_m1_m2_5points(4,j) + &
                              hprimewgll_xxT(i,5)*C2_m1_m2_5points(5,j)

        E3_m1_m2_5points(i,j) = hprimewgll_xxT(i,1)*C3_m1_m2_5points(1,j) + &
                              hprimewgll_xxT(i,2)*C3_m1_m2_5points(2,j) + &
                              hprimewgll_xxT(i,3)*C3_m1_m2_5points(3,j) + &
                              hprimewgll_xxT(i,4)*C3_m1_m2_5points(4,j) + &
                              hprimewgll_xxT(i,5)*C3_m1_m2_5points(5,j)
      enddo
    enddo

    do i=1,m1
      do j=1,m1
        ! for efficiency it is better to leave this loop on k inside, it leads to slightly faster code
        do k = 1,NGLLX
          newtempx2(i,j,k) = tempx2(i,1,k)*hprimewgll_xx(1,j) + &
                           tempx2(i,2,k)*hprimewgll_xx(2,j) + &
                           tempx2(i,3,k)*hprimewgll_xx(3,j) + &
                           tempx2(i,4,k)*hprimewgll_xx(4,j) + &
                           tempx2(i,5,k)*hprimewgll_xx(5,j)

          newtempy2(i,j,k) = tempy2(i,1,k)*hprimewgll_xx(1,j) + &
                           tempy2(i,2,k)*hprimewgll_xx(2,j) + &
                           tempy2(i,3,k)*hprimewgll_xx(3,j) + &
                           tempy2(i,4,k)*hprimewgll_xx(4,j) + &
                           tempy2(i,5,k)*hprimewgll_xx(5,j)

          newtempz2(i,j,k) = tempz2(i,1,k)*hprimewgll_xx(1,j) + &
                           tempz2(i,2,k)*hprimewgll_xx(2,j) + &
                           tempz2(i,3,k)*hprimewgll_xx(3,j) + &
                           tempz2(i,4,k)*hprimewgll_xx(4,j) + &
                           tempz2(i,5,k)*hprimewgll_xx(5,j)
        enddo
      enddo
    enddo

    do j=1,m1
      do i=1,m2
        E1_mxm_m2_m1_5points(i,j) = C1_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                  C1_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                  C1_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                  C1_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                  C1_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

        E2_mxm_m2_m1_5points(i,j) = C2_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                  C2_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                  C2_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                  C2_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                  C2_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)

        E3_mxm_m2_m1_5points(i,j) = C3_mxm_m2_m1_5points(i,1)*hprimewgll_xx(1,j) + &
                                  C3_mxm_m2_m1_5points(i,2)*hprimewgll_xx(2,j) + &
                                  C3_mxm_m2_m1_5points(i,3)*hprimewgll_xx(3,j) + &
                                  C3_mxm_m2_m1_5points(i,4)*hprimewgll_xx(4,j) + &
                                  C3_mxm_m2_m1_5points(i,5)*hprimewgll_xx(5,j)
      enddo
    enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        fac1 = wgllwgll_yz(j,k)
        do i=1,NGLLX
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          ! sum contributions
          sum_terms(1,i,j,k) = - (fac1*newtempx1(i,j,k) + fac2*newtempx2(i,j,k) + fac3*newtempx3(i,j,k))
          sum_terms(2,i,j,k) = - (fac1*newtempy1(i,j,k) + fac2*newtempy2(i,j,k) + fac3*newtempy3(i,j,k))
          sum_terms(3,i,j,k) = - (fac1*newtempz1(i,j,k) + fac2*newtempz2(i,j,k) + fac3*newtempz3(i,j,k))

          if(GRAVITY_VAL) sum_terms(:,i,j,k) = sum_terms(:,i,j,k) + rho_s_H(:,i,j,k)

        enddo ! NGLLX
      enddo ! NGLLY
    enddo ! NGLLZ

    ! sum contributions from each element to the global mesh and add gravity terms
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob1 = ibool(i,j,k,ispec)
          accel_crust_mantle(:,iglob1) = accel_crust_mantle(:,iglob1) + sum_terms(:,i,j,k)
        enddo
      enddo
    enddo

    ! update memory variables based upon the Runge-Kutta scheme
    ! convention for attenuation
    ! term in xx = 1
    ! term in yy = 2
    ! term in xy = 3
    ! term in xz = 4
    ! term in yz = 5
    ! term in zz not computed since zero trace
    ! This is because we only implement Q_\mu attenuation and not Q_\kappa.
    ! Note that this does *NOT* imply that there is no attenuation for P waves
    ! because for Q_\kappa = infinity one gets (see for instance Dahlen and Tromp (1998)
    ! equation (9.59) page 350): Q_\alpha = Q_\mu * 3 * (V_p/V_s)^2 / 4
    ! therefore Q_\alpha is not zero; for instance for V_p / V_s = sqrt(3)
    ! we get Q_\alpha = (9 / 4) * Q_\mu = 2.25 * Q_\mu

    if(ATTENUATION_VAL .and. ( PARTIAL_PHYS_DISPERSION_ONLY .eqv. .false. ) ) then

      call compute_element_strain_att_Dev(ispec,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE,displ_crust_mantle,veloc_crust_mantle,&
                                          deltat,ibool,hprime_xx,hprime_xxT,&
                                          xix,xiy,xiz,etax,etay,etaz,gammax,gammay,gammaz,epsilondev_loc_nplus1)
      ! updates R_memory
      call compute_element_att_memory_cr(ispec,R_memory, &
                                      vx,vy,vz,vnspec,factor_common, &
                                      alphaval,betaval,gammaval, &
                                      c44store,muvstore, &
!ZN                                      epsilondev,epsilondev_loc)
                                      epsilondev_loc_nplus1,epsilondev_loc)

    endif

    ! save deviatoric strain for Runge-Kutta scheme
!ZN    if(COMPUTE_AND_STORE_STRAIN) then
!ZN      epsilondev(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)
!ZN    endif
! end ispec loop
   enddo
!$OMP enddo

! end ispec_globe strided loop
  enddo   ! spectral element loop NSPEC_CRUST_MANTLE
!$OMP END PARALLEL

  end subroutine compute_forces_crust_mantle_Dev

