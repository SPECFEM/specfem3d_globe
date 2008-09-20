!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            August 2008
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

! compute forces in the crust/mantle and in the inner core (i.e., all the solid regions)
  subroutine compute_forces_CM_IC(displ_crust_mantle,accel_crust_mantle, &
          xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
          xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
          gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
          kappavstore_crust_mantle,kappahstore_crust_mantle,muvstore_crust_mantle,muhstore_crust_mantle, &
          eta_anisostore_crust_mantle, &
          ibool_crust_mantle,idoubling_crust_mantle,R_memory_crust_mantle,epsilondev_crust_mantle,one_minus_sum_beta_crust_mantle, &
          factor_common_crust_mantle,vx_crust_mantle,vy_crust_mantle,vz_crust_mantle,vnspec_crust_mantle, &
          displ_inner_core,accel_inner_core, &
          xix_inner_core,xiy_inner_core,xiz_inner_core,etax_inner_core,etay_inner_core,etaz_inner_core, &
          gammax_inner_core,gammay_inner_core,gammaz_inner_core, &
          kappavstore_inner_core,muvstore_inner_core,ibool_inner_core,idoubling_inner_core, &
          R_memory_inner_core,epsilondev_inner_core,&
          one_minus_sum_beta_inner_core,factor_common_inner_core, &
          vx_inner_core,vy_inner_core,vz_inner_core,vnspec_inner_core, &
          hprime_xx,hprime_yy,hprime_zz, &
          hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
          wgllwgll_xy,wgllwgll_xz,wgllwgll_yz, &
          alphaval,betaval,gammaval, &
#ifdef USE_MPI
            is_on_a_slice_edge_crust_mantle,is_on_a_slice_edge_inner_core, &
            myrank,iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all, &
            buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector,iphase, &
               nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
               npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
               receiver_cube_from_slices,ibelm_bottom_inner_core,NSPEC2D_BOTTOM_INNER_CORE,INCLUDE_CENTRAL_CUBE,iphase_CC, &
#endif
          COMPUTE_AND_STORE_STRAIN,AM_V,icall)

  implicit none

  include "constants.h"

! include values created by the mesher
! done for performance only using static allocation to allow for loop unrolling
  include "values_from_mesher.h"

  integer :: icall

! attenuation_model_variables
  type attenuation_model_variables
    sequence
    double precision min_period, max_period
    double precision                          :: QT_c_source        ! Source Frequency
    double precision, dimension(N_SLS)        :: Qtau_s             ! tau_sigma
    double precision, dimension(:), pointer   :: QrDisc             ! Discontinutitues Defined
    double precision, dimension(:), pointer   :: Qr                 ! Radius
    integer, dimension(:), pointer            :: interval_Q                 ! Steps
    double precision, dimension(:), pointer   :: Qmu                ! Shear Attenuation
    double precision, dimension(:,:), pointer :: Qtau_e             ! tau_epsilon
    double precision, dimension(:), pointer   :: Qone_minus_sum_beta, Qone_minus_sum_beta2      ! one_minus_sum_beta
    double precision, dimension(:,:), pointer :: Qfc, Qfc2          ! factor_common
    double precision, dimension(:), pointer   :: Qsf, Qsf2          ! scale_factor
    integer, dimension(:), pointer            :: Qrmin              ! Max and Mins of idoubling
    integer, dimension(:), pointer            :: Qrmax              ! Max and Mins of idoubling
    integer                                   :: Qn                 ! Number of points
  end type attenuation_model_variables

  type (attenuation_model_variables) AM_V
! attenuation_model_variables

! for forward or backward simulations
  logical COMPUTE_AND_STORE_STRAIN

! array with the local to global mapping per slice
  integer, dimension(NSPEC_CRUST_MANTLE) :: idoubling_crust_mantle

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: displ_crust_mantle,accel_crust_mantle

! memory variables for attenuation
! memory variables R_ij are stored at the local rather than global level
! to allow for optimization of cache access by compiler
  integer i_sls,i_memory

! variable lengths for factor_common and one_minus_sum_beta
  integer :: vx_crust_mantle,vy_crust_mantle,vz_crust_mantle,vnspec_crust_mantle

  real(kind=CUSTOM_REAL) one_minus_sum_beta_use,minus_sum_beta
  real(kind=CUSTOM_REAL), dimension(vx_crust_mantle,vy_crust_mantle,vz_crust_mantle,vnspec_crust_mantle) :: &
       one_minus_sum_beta_crust_mantle

  integer iregion_selected

! for attenuation
  real(kind=CUSTOM_REAL) R_xx_val,R_yy_val
  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ATTENUAT) :: R_memory_crust_mantle
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: epsilondev_crust_mantle

! [alpha,beta,gamma]val reduced to N_SLS and factor_common to N_SLS*NUM_NODES
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval,betaval,gammaval
  real(kind=CUSTOM_REAL), dimension(N_SLS,vx_crust_mantle,vy_crust_mantle,vz_crust_mantle,vnspec_crust_mantle) :: &
       factor_common_crust_mantle

  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ) :: epsilondev_loc

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
        gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

! array with derivatives of Lagrange polynomials and precalculated products
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLY) :: hprime_yy,hprimewgll_yy
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ) :: &
    tempx1,tempx2,tempx3,tempy1,tempy2,tempy3,tempz1,tempz2,tempz3

! x y and z contain r theta and phi
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        kappavstore_crust_mantle,muvstore_crust_mantle

! store anisotropic properties only where needed to save memory
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        kappahstore_crust_mantle,muhstore_crust_mantle,eta_anisostore_crust_mantle

  integer ispec,iglob
  integer i,j,k,l

! the 21 coefficients for an anisotropic medium in reduced notation
  real(kind=CUSTOM_REAL) c11,c22,c33,c44,c55,c66,c12,c13,c23,c14,c24,c34,c15,c25,c35,c45,c16,c26,c36,c46,c56

  real(kind=CUSTOM_REAL) rhovphsq,sinphifour,cosphisq,sinphisq,costhetasq,rhovsvsq,sinthetasq, &
        cosphifour,costhetafour,rhovpvsq,sinthetafour,rhovshsq,cosfourphi, &
        costwotheta,cosfourtheta,sintwophisq,costheta,sinphi,sintheta,cosphi, &
        sintwotheta,costwophi,sintwophi,costwothetasq,costwophisq,phi,theta

  real(kind=CUSTOM_REAL) two_rhovpvsq,two_rhovphsq,two_rhovsvsq,two_rhovshsq
  real(kind=CUSTOM_REAL) four_rhovpvsq,four_rhovphsq,four_rhovsvsq,four_rhovshsq

  real(kind=CUSTOM_REAL) twoetaminone,etaminone,eta_aniso
  real(kind=CUSTOM_REAL) two_eta_aniso,four_eta_aniso,six_eta_aniso

  real(kind=CUSTOM_REAL) xixl,xiyl,xizl,etaxl,etayl,etazl,gammaxl,gammayl,gammazl,jacobianl
  real(kind=CUSTOM_REAL) duxdxl,duxdyl,duxdzl,duydxl,duydyl,duydzl,duzdxl,duzdyl,duzdzl

  real(kind=CUSTOM_REAL) duxdxl_plus_duydyl,duxdxl_plus_duzdzl,duydyl_plus_duzdzl
  real(kind=CUSTOM_REAL) duxdyl_plus_duydxl,duzdxl_plus_duxdzl,duzdyl_plus_duydzl

  real(kind=CUSTOM_REAL) sigma_xx,sigma_yy,sigma_zz,sigma_xy,sigma_xz,sigma_yz

  real(kind=CUSTOM_REAL) hp1,hp2,hp3
  real(kind=CUSTOM_REAL) fac1,fac2,fac3
  real(kind=CUSTOM_REAL) lambdal,mul,lambdalplus2mul
  real(kind=CUSTOM_REAL) kappal,kappavl,kappahl,muvl,muhl
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sum_terms

  real(kind=CUSTOM_REAL) tempx1l,tempx2l,tempx3l
  real(kind=CUSTOM_REAL) tempy1l,tempy2l,tempy3l
  real(kind=CUSTOM_REAL) tempz1l,tempz2l,tempz3l

  real(kind=CUSTOM_REAL) radius_cr

!=====================================================================
!=====================================================================
!=====================================================================

! same attenuation everywhere in the inner core therefore no need to use Brian's routines
  integer, parameter :: iregion_selected_inner_core = 1

! displacement and acceleration
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_INNER_CORE) :: displ_inner_core,accel_inner_core

! variable lengths for factor_common and one_minus_sum_beta
  integer :: vx_inner_core,vy_inner_core,vz_inner_core,vnspec_inner_core

  real(kind=CUSTOM_REAL), dimension(vx_inner_core,vy_inner_core,vz_inner_core,vnspec_inner_core) :: one_minus_sum_beta_inner_core

  real(kind=CUSTOM_REAL), dimension(N_SLS,vx_inner_core,vy_inner_core,vz_inner_core,vnspec_inner_core) :: factor_common_inner_core

  real(kind=CUSTOM_REAL), dimension(5,N_SLS,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE_ATTENUATION) :: R_memory_inner_core
  real(kind=CUSTOM_REAL), dimension(5,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: epsilondev_inner_core

! array with the local to global mapping per slice
  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core

! arrays with mesh parameters per slice
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
       xix_inner_core,xiy_inner_core,xiz_inner_core,etax_inner_core,etay_inner_core,etaz_inner_core,&
       gammax_inner_core,gammay_inner_core,gammaz_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: kappavstore_inner_core,muvstore_inner_core

  integer :: computed_elements

#ifdef USE_MPI
  integer :: iphase

  logical, dimension(NSPEC_CRUST_MANTLE) :: is_on_a_slice_edge_crust_mantle
  logical, dimension(NSPEC_INNER_CORE) :: is_on_a_slice_edge_inner_core

  integer :: ichunk,iproc_xi,iproc_eta,myrank

  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_CM) :: iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_CM) :: iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle

  integer, dimension(NGLOB2DMAX_XMIN_XMAX_IC) :: iboolleft_xi_inner_core,iboolright_xi_inner_core
  integer, dimension(NGLOB2DMAX_YMIN_YMAX_IC) :: iboolleft_eta_inner_core,iboolright_eta_inner_core

  integer npoin2D_faces_crust_mantle(NUMFACES_SHARED)
  integer npoin2D_faces_inner_core(NUMFACES_SHARED)

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle
  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_inner_core,npoin2D_eta_inner_core

! communication pattern for faces between chunks
  integer, dimension(NUMMSGS_FACES_VAL) :: iprocfrom_faces,iprocto_faces

! communication pattern for corners between chunks
  integer, dimension(NCORNERSCHUNKS_VAL) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  integer, dimension(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED) :: iboolcorner_crust_mantle
  integer, dimension(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED) :: iboolcorner_inner_core

  integer, dimension(NGLOB2DMAX_XY_VAL_CM,NUMFACES_SHARED) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_VAL_IC,NUMFACES_SHARED) :: iboolfaces_inner_core

  integer :: npoin2D_max_all

  real(kind=CUSTOM_REAL), dimension(NDIM,npoin2D_max_all) :: buffer_send_faces,buffer_received_faces

! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC) :: &
     buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector

! for matching with central cube in inner core
  integer nb_msgs_theor_in_cube, npoin2D_cube_from_slices,iphase_CC
  integer, dimension(nb_msgs_theor_in_cube) :: sender_from_slices_to_cube
  double precision, dimension(npoin2D_cube_from_slices,NDIM) :: buffer_slices
  double precision, dimension(npoin2D_cube_from_slices,NDIM,nb_msgs_theor_in_cube) :: buffer_all_cube_from_slices
  integer, dimension(nb_msgs_theor_in_cube,npoin2D_cube_from_slices):: ibool_central_cube
  integer receiver_cube_from_slices
  logical :: INCLUDE_CENTRAL_CUBE

! local to global mapping
  integer NSPEC2D_BOTTOM_INNER_CORE
  integer, dimension(NSPEC2D_BOTTOM_INNER_CORE) :: ibelm_bottom_inner_core

#endif

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

! set acceleration to zero
  if(icall == 1) accel_crust_mantle(:,:) = 0._CUSTOM_REAL

  computed_elements = 0

  do ispec = 1,NSPEC_CRUST_MANTLE

#ifdef USE_MPI
! hide communications by computing the edges first
    if((icall == 2 .and. is_on_a_slice_edge_crust_mantle(ispec)) .or. &
       (icall == 1 .and. .not. is_on_a_slice_edge_crust_mantle(ispec))) cycle

! process the communications every ELEMENTS_NONBLOCKING elements
    computed_elements = computed_elements + 1
    if (USE_NONBLOCKING_COMMS .and. icall == 2 .and. mod(computed_elements,ELEMENTS_NONBLOCKING_CM_IC) == 0) then

      if(iphase <= 7) call assemble_MPI_vector(myrank,accel_crust_mantle,accel_inner_core, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle(1),npoin2D_eta_crust_mantle(1), &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core(1),npoin2D_eta_inner_core(1), &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all, &
            buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector, &
            NUMMSGS_FACES_VAL,NCORNERSCHUNKS_VAL, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL_CM, &
            NGLOB1D_RADIAL_IC,NCHUNKS_VAL,iphase)

      if(INCLUDE_CENTRAL_CUBE) then
          if(iphase > 7 .and. iphase_CC <= 4) &
            call assemble_MPI_central_cube(ichunk,nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
              npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
              receiver_cube_from_slices,ibool_inner_core,idoubling_inner_core, &
              ibelm_bottom_inner_core,NSPEC2D_BOTTOM_INNER_CORE,accel_inner_core,NDIM,iphase_CC)
      endif

    endif
#endif

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool_crust_mantle(l,j,k,ispec)
            tempx1l = tempx1l + displ_crust_mantle(1,iglob)*hp1
            tempy1l = tempy1l + displ_crust_mantle(2,iglob)*hp1
            tempz1l = tempz1l + displ_crust_mantle(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            hp2 = hprime_yy(j,l)
            iglob = ibool_crust_mantle(i,l,k,ispec)
            tempx2l = tempx2l + displ_crust_mantle(1,iglob)*hp2
            tempy2l = tempy2l + displ_crust_mantle(2,iglob)*hp2
            tempz2l = tempz2l + displ_crust_mantle(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool_crust_mantle(i,j,l,ispec)
            tempx3l = tempx3l + displ_crust_mantle(1,iglob)*hp3
            tempy3l = tempy3l + displ_crust_mantle(2,iglob)*hp3
            tempz3l = tempz3l + displ_crust_mantle(3,iglob)*hp3
          enddo

!         get derivatives of ux, uy and uz with respect to x, y and z

          xixl = xix_crust_mantle(i,j,k,ispec)
          xiyl = xiy_crust_mantle(i,j,k,ispec)
          xizl = xiz_crust_mantle(i,j,k,ispec)
          etaxl = etax_crust_mantle(i,j,k,ispec)
          etayl = etay_crust_mantle(i,j,k,ispec)
          etazl = etaz_crust_mantle(i,j,k,ispec)
          gammaxl = gammax_crust_mantle(i,j,k,ispec)
          gammayl = gammay_crust_mantle(i,j,k,ispec)
          gammazl = gammaz_crust_mantle(i,j,k,ispec)

! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! compute deviatoric strain
 if (COMPUTE_AND_STORE_STRAIN) then
    epsilondev_loc(1,i,j,k) = duxdxl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilondev_loc(2,i,j,k) = duydyl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
    epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
    epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
  endif

! precompute terms for attenuation if needed
  if(ATTENUATION_VAL) then
    radius_cr = xstore_crust_mantle(ibool_crust_mantle(i,j,k,ispec))
    call get_attenuation_index(idoubling_crust_mantle(ispec), dble(radius_cr), iregion_selected, .FALSE., AM_V)
    one_minus_sum_beta_use = one_minus_sum_beta_crust_mantle(1,1,1,iregion_selected)
    minus_sum_beta =  one_minus_sum_beta_use - 1.0
  endif

!
! compute either isotropic or anisotropic elements
!

  if(ANISOTROPIC_3D_MANTLE_VAL) then

  else

! do not use transverse isotropy except if element is between d220 and Moho
  if(.not. (TRANSVERSE_ISOTROPY_VAL .and. (idoubling_crust_mantle(ispec)==IFLAG_220_80 .or. &
                                           idoubling_crust_mantle(ispec)==IFLAG_80_MOHO))) then

! layer with no transverse isotropy, use kappav and muv
          kappal = kappavstore_crust_mantle(i,j,k,ispec)
          mul = muvstore_crust_mantle(i,j,k,ispec)

! use unrelaxed parameters if attenuation
    if(ATTENUATION_VAL) mul = mul * one_minus_sum_beta_use

          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2.*mul

! compute stress sigma

          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

          sigma_xy = mul*duxdyl_plus_duydxl
          sigma_xz = mul*duzdxl_plus_duxdzl
          sigma_yz = mul*duzdyl_plus_duydzl

    else

! use Kappa and mu from transversely isotropic model
      kappavl = kappavstore_crust_mantle(i,j,k,ispec)
      muvl = muvstore_crust_mantle(i,j,k,ispec)

      kappahl = kappahstore_crust_mantle(i,j,k,ispec)
      muhl = muhstore_crust_mantle(i,j,k,ispec)

! use unrelaxed parameters if attenuation
! eta does not need to be shifted since it is a ratio
  if(ATTENUATION_VAL) then
    muvl = muvl * one_minus_sum_beta_use
    muhl = muhl * one_minus_sum_beta_use
  endif

  rhovpvsq = kappavl + FOUR_THIRDS * muvl  !!! that is C
  rhovphsq = kappahl + FOUR_THIRDS * muhl  !!! that is A

  rhovsvsq = muvl  !!! that is L
  rhovshsq = muhl  !!! that is N

  eta_aniso = eta_anisostore_crust_mantle(i,j,k,ispec)  !!! that is  F / (A - 2 L)

! use mesh coordinates to get theta and phi
! ystore_crust_mantle and zstore_crust_mantle contain theta and phi

  iglob = ibool_crust_mantle(i,j,k,ispec)
  theta = ystore_crust_mantle(iglob)
  phi = zstore_crust_mantle(iglob)

  costheta = cos(theta)
  sintheta = sin(theta)
  cosphi = cos(phi)
  sinphi = sin(phi)

  costhetasq = costheta * costheta
  sinthetasq = sintheta * sintheta
  cosphisq = cosphi * cosphi
  sinphisq = sinphi * sinphi

  costhetafour = costhetasq * costhetasq
  sinthetafour = sinthetasq * sinthetasq
  cosphifour = cosphisq * cosphisq
  sinphifour = sinphisq * sinphisq

  costwotheta = cos(2.*theta)
  sintwotheta = sin(2.*theta)
  costwophi = cos(2.*phi)
  sintwophi = sin(2.*phi)

  cosfourtheta = cos(4.*theta)
  cosfourphi = cos(4.*phi)

  costwothetasq = costwotheta * costwotheta

  costwophisq = costwophi * costwophi
  sintwophisq = sintwophi * sintwophi

  etaminone = eta_aniso - 1.
  twoetaminone = 2. * eta_aniso - 1.

! precompute some products to reduce the CPU time

      two_eta_aniso = 2.*eta_aniso
      four_eta_aniso = 4.*eta_aniso
      six_eta_aniso = 6.*eta_aniso

      two_rhovpvsq = 2.*rhovpvsq
      two_rhovphsq = 2.*rhovphsq
      two_rhovsvsq = 2.*rhovsvsq
      two_rhovshsq = 2.*rhovshsq

      four_rhovpvsq = 4.*rhovpvsq
      four_rhovphsq = 4.*rhovphsq
      four_rhovsvsq = 4.*rhovsvsq
      four_rhovshsq = 4.*rhovshsq

! the 21 anisotropic coefficients computed using Mathematica

 c11 = rhovphsq*sinphifour + 2.*cosphisq*sinphisq* &
   (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
      sinthetasq) + cosphifour* &
   (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
      costhetasq*sinthetasq + rhovpvsq*sinthetafour)

 c12 = ((rhovphsq - two_rhovshsq)*(3. + cosfourphi)*costhetasq)/4. - &
  four_rhovshsq*cosphisq*costhetasq*sinphisq + &
  (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. + &
  eta_aniso*(rhovphsq - two_rhovsvsq)*(cosphifour + &
     2.*cosphisq*costhetasq*sinphisq + sinphifour)*sinthetasq + &
  rhovpvsq*cosphisq*sinphisq*sinthetafour - &
  rhovsvsq*sintwophisq*sinthetafour

 c13 = (cosphisq*(rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - &
       12.*eta_aniso*rhovsvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - &
          four_eta_aniso*rhovsvsq)*cosfourtheta))/8. + &
  sinphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
     (rhovphsq - two_rhovshsq)*sinthetasq)

 c14 = costheta*sinphi*((cosphisq* &
       (-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
         (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
            four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
    (etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))*sinphisq)* sintheta

 c15 = cosphi*costheta*((cosphisq* (-rhovphsq + rhovpvsq + &
         (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
          costwotheta))/2. + etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sintheta

 c16 = (cosphi*sinphi*(cosphisq* (-rhovphsq + rhovpvsq + &
         (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
            four_eta_aniso*rhovsvsq)*costwotheta) + &
      2.*etaminone*(rhovphsq - two_rhovsvsq)*sinphisq)*sinthetasq)/2.

 c22 = rhovphsq*cosphifour + 2.*cosphisq*sinphisq* &
   (rhovphsq*costhetasq + (eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
      sinthetasq) + sinphifour* &
   (rhovphsq*costhetafour + 2.*(eta_aniso*rhovphsq + two_rhovsvsq - two_eta_aniso*rhovsvsq)* &
      costhetasq*sinthetasq + rhovpvsq*sinthetafour)

 c23 = ((rhovphsq + six_eta_aniso*rhovphsq + rhovpvsq - four_rhovsvsq - 12.*eta_aniso*rhovsvsq + &
       (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
        cosfourtheta)*sinphisq)/8. + &
  cosphisq*(eta_aniso*(rhovphsq - two_rhovsvsq)*costhetasq + &
     (rhovphsq - two_rhovshsq)*sinthetasq)

 c24 = costheta*sinphi*(etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
    ((-rhovphsq + rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + &
            four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

 c25 = cosphi*costheta*((etaminone*rhovphsq + 2.*(rhovshsq - eta_aniso*rhovsvsq))* &
     cosphisq + ((-rhovphsq + rhovpvsq + four_rhovshsq - four_rhovsvsq + &
         (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
            four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)*sintheta

 c26 = (cosphi*sinphi*(2.*etaminone*(rhovphsq - two_rhovsvsq)*cosphisq + &
      (-rhovphsq + rhovpvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + four_rhovsvsq - &
            four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)*sinthetasq)/2.

 c33 = rhovpvsq*costhetafour + 2.*(eta_aniso*(rhovphsq - two_rhovsvsq) + two_rhovsvsq)* &
   costhetasq*sinthetasq + rhovphsq*sinthetafour

 c34 = -((rhovphsq - rhovpvsq + (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq &
           - four_eta_aniso*rhovsvsq)*costwotheta)*sinphi*sintwotheta)/4.

 c35 = -(cosphi*(rhovphsq - rhovpvsq + &
       (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
        costwotheta)*sintwotheta)/4.

 c36 = -((rhovphsq - rhovpvsq - four_rhovshsq + four_rhovsvsq + &
       (twoetaminone*rhovphsq - rhovpvsq + four_rhovsvsq - four_eta_aniso*rhovsvsq)* &
        costwotheta)*sintwophi*sinthetasq)/4.

 c44 = cosphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
  sinphisq*(rhovsvsq*costwothetasq + &
     (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

 c45 = ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
      four_eta_aniso*rhovsvsq + (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + &
         4.*etaminone*rhovsvsq)*costwotheta)*sintwophi*sinthetasq)/4.

 c46 = -(cosphi*costheta*((rhovshsq - rhovsvsq)*cosphisq - &
      ((rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
           four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
              four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta)*sinphisq)/2.)* sintheta)

 c55 = sinphisq*(rhovsvsq*costhetasq + rhovshsq*sinthetasq) + &
  cosphisq*(rhovsvsq*costwothetasq + &
     (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq + four_eta_aniso*rhovsvsq)*costhetasq* sinthetasq)

 c56 = costheta*sinphi*((cosphisq* &
       (rhovphsq - two_eta_aniso*rhovphsq + rhovpvsq - two_rhovshsq - two_rhovsvsq + &
         four_eta_aniso*rhovsvsq + (-rhovphsq + two_eta_aniso*rhovphsq - rhovpvsq + &
            four_rhovsvsq - four_eta_aniso*rhovsvsq)*costwotheta))/2. + &
    (-rhovshsq + rhovsvsq)*sinphisq)*sintheta

 c66 = rhovshsq*costwophisq*costhetasq - &
  2.*(rhovphsq - two_rhovshsq)*cosphisq*costhetasq*sinphisq + &
  (rhovphsq*(11. + 4.*costwotheta + cosfourtheta)*sintwophisq)/32. - &
  (rhovsvsq*(-6. - 2.*cosfourphi + cos(4.*phi - 2.*theta) - 2.*costwotheta + &
       cos(2.*(2.*phi + theta)))*sinthetasq)/8. + &
  rhovpvsq*cosphisq*sinphisq*sinthetafour - &
  (eta_aniso*(rhovphsq - two_rhovsvsq)*sintwophisq*sinthetafour)/2.

! general expression of stress tensor for full Cijkl with 21 coefficients

     sigma_xx = c11*duxdxl + c16*duxdyl_plus_duydxl + c12*duydyl + &
               c15*duzdxl_plus_duxdzl + c14*duzdyl_plus_duydzl + c13*duzdzl

     sigma_yy = c12*duxdxl + c26*duxdyl_plus_duydxl + c22*duydyl + &
               c25*duzdxl_plus_duxdzl + c24*duzdyl_plus_duydzl + c23*duzdzl

     sigma_zz = c13*duxdxl + c36*duxdyl_plus_duydxl + c23*duydyl + &
               c35*duzdxl_plus_duxdzl + c34*duzdyl_plus_duydzl + c33*duzdzl

     sigma_xy = c16*duxdxl + c66*duxdyl_plus_duydxl + c26*duydyl + &
               c56*duzdxl_plus_duxdzl + c46*duzdyl_plus_duydzl + c36*duzdzl

     sigma_xz = c15*duxdxl + c56*duxdyl_plus_duydxl + c25*duydyl + &
               c55*duzdxl_plus_duxdzl + c45*duzdyl_plus_duydzl + c35*duzdzl

     sigma_yz = c14*duxdxl + c46*duxdyl_plus_duydxl + c24*duydyl + &
               c45*duzdxl_plus_duxdzl + c44*duzdyl_plus_duydzl + c34*duzdzl

  endif

  endif   ! end of test whether isotropic or anisotropic element

! subtract memory variables if attenuation
  if(ATTENUATION_VAL) then
    do i_sls = 1,N_SLS
      R_xx_val = R_memory_crust_mantle(1,i_sls,i,j,k,ispec)
      R_yy_val = R_memory_crust_mantle(2,i_sls,i,j,k,ispec)
      sigma_xx = sigma_xx - R_xx_val
      sigma_yy = sigma_yy - R_yy_val
      sigma_zz = sigma_zz + R_xx_val + R_yy_val
      sigma_xy = sigma_xy - R_memory_crust_mantle(3,i_sls,i,j,k,ispec)
      sigma_xz = sigma_xz - R_memory_crust_mantle(4,i_sls,i,j,k,ispec)
      sigma_yz = sigma_yz - R_memory_crust_mantle(5,i_sls,i,j,k,ispec)
    enddo
  endif

! form dot product with test vector, symmetric form
      tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl)
      tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl)
      tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

      tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl)
      tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl)
      tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

      tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl)
      tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl)
      tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)

          enddo
        enddo
      enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempy1l = 0._CUSTOM_REAL
          tempz1l = 0._CUSTOM_REAL

          tempx2l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL

          tempx3l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1l = tempx1l + tempx1(l,j,k)*fac1
            tempy1l = tempy1l + tempy1(l,j,k)*fac1
            tempz1l = tempz1l + tempz1(l,j,k)*fac1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            fac2 = hprimewgll_yy(l,j)
            tempx2l = tempx2l + tempx2(i,l,k)*fac2
            tempy2l = tempy2l + tempy2(i,l,k)*fac2
            tempz2l = tempz2l + tempz2(i,l,k)*fac2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            fac3 = hprimewgll_zz(l,k)
            tempx3l = tempx3l + tempx3(i,j,l)*fac3
            tempy3l = tempy3l + tempy3(i,j,l)*fac3
            tempz3l = tempz3l + tempz3(i,j,l)*fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          sum_terms(1,i,j,k) = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
          sum_terms(2,i,j,k) = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
          sum_terms(3,i,j,k) = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

        enddo
      enddo
    enddo

! sum contributions from each element to the global mesh and add gravity terms
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec)
          accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sum_terms(1,i,j,k)
          accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sum_terms(2,i,j,k)
          accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sum_terms(3,i,j,k)
        enddo
      enddo
    enddo

! update memory variables based upon a Runge-Kutta scheme.
! convention for attenuation:
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace
    if(ATTENUATION_VAL) then
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            do i_sls = 1,N_SLS
              do i_memory = 1,5
! get coefficients for that standard linear solid
! IMPROVE we use mu_v here even if there is some anisotropy
! IMPROVE we should probably use an average value instead
                R_memory_crust_mantle(i_memory,i_sls,i,j,k,ispec) = alphaval(i_sls) * &
                  R_memory_crust_mantle(i_memory,i_sls,i,j,k,ispec) + &
                  factor_common_crust_mantle(i_sls,1,1,1,iregion_selected) * muvstore_crust_mantle(i,j,k,ispec) * &
                  (betaval(i_sls) * epsilondev_crust_mantle(i_memory,i,j,k,ispec) + &
                  gammaval(i_sls) * epsilondev_loc(i_memory,i,j,k))
              enddo
            enddo
          enddo
        enddo
      enddo
    endif

! save deviatoric strain for Runge-Kutta scheme
  if(COMPUTE_AND_STORE_STRAIN) epsilondev_crust_mantle(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)

  enddo   ! spectral element loop

!=====================================================================
!=====================================================================
!=====================================================================

! ****************************************************
!   big loop over all spectral elements in the solid
! ****************************************************

! set acceleration to zero
  if(icall == 1) accel_inner_core(:,:) = 0._CUSTOM_REAL

  computed_elements = 0

  do ispec = 1,NSPEC_INNER_CORE

#ifdef USE_MPI
! hide communications by computing the edges first
    if((icall == 2 .and. is_on_a_slice_edge_inner_core(ispec)) .or. &
       (icall == 1 .and. .not. is_on_a_slice_edge_inner_core(ispec))) cycle
#endif

! exclude fictitious elements in central cube
    if(idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

#ifdef USE_MPI
! process the communications every ELEMENTS_NONBLOCKING elements
    computed_elements = computed_elements + 1
    if (USE_NONBLOCKING_COMMS .and. icall == 2 .and. mod(computed_elements,ELEMENTS_NONBLOCKING_CM_IC) == 0) then

         if(iphase <= 7) call assemble_MPI_vector(myrank,accel_crust_mantle,accel_inner_core, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle(1),npoin2D_eta_crust_mantle(1), &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core(1),npoin2D_eta_inner_core(1), &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces,buffer_received_faces,npoin2D_max_all, &
            buffer_send_chunkcorners_vector,buffer_recv_chunkcorners_vector, &
            NUMMSGS_FACES_VAL,NCORNERSCHUNKS_VAL, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL_CM, &
            NGLOB1D_RADIAL_IC,NCHUNKS_VAL,iphase)

      if(INCLUDE_CENTRAL_CUBE) then
          if(iphase > 7 .and. iphase_CC <= 4) &
            call assemble_MPI_central_cube(ichunk,nb_msgs_theor_in_cube,sender_from_slices_to_cube, &
              npoin2D_cube_from_slices,buffer_all_cube_from_slices,buffer_slices,ibool_central_cube, &
              receiver_cube_from_slices,ibool_inner_core,idoubling_inner_core, &
              ibelm_bottom_inner_core,NSPEC2D_BOTTOM_INNER_CORE,accel_inner_core,NDIM,iphase_CC)
      endif

    endif
#endif

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempx2l = 0._CUSTOM_REAL
          tempx3l = 0._CUSTOM_REAL

          tempy1l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL

          tempz1l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            hp1 = hprime_xx(i,l)
            iglob = ibool_inner_core(l,j,k,ispec)
            tempx1l = tempx1l + displ_inner_core(1,iglob)*hp1
            tempy1l = tempy1l + displ_inner_core(2,iglob)*hp1
            tempz1l = tempz1l + displ_inner_core(3,iglob)*hp1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            hp2 = hprime_yy(j,l)
            iglob = ibool_inner_core(i,l,k,ispec)
            tempx2l = tempx2l + displ_inner_core(1,iglob)*hp2
            tempy2l = tempy2l + displ_inner_core(2,iglob)*hp2
            tempz2l = tempz2l + displ_inner_core(3,iglob)*hp2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            hp3 = hprime_zz(k,l)
            iglob = ibool_inner_core(i,j,l,ispec)
            tempx3l = tempx3l + displ_inner_core(1,iglob)*hp3
            tempy3l = tempy3l + displ_inner_core(2,iglob)*hp3
            tempz3l = tempz3l + displ_inner_core(3,iglob)*hp3
          enddo

!         get derivatives of ux, uy and uz with respect to x, y and z

          xixl = xix_inner_core(i,j,k,ispec)
          xiyl = xiy_inner_core(i,j,k,ispec)
          xizl = xiz_inner_core(i,j,k,ispec)
          etaxl = etax_inner_core(i,j,k,ispec)
          etayl = etay_inner_core(i,j,k,ispec)
          etazl = etaz_inner_core(i,j,k,ispec)
          gammaxl = gammax_inner_core(i,j,k,ispec)
          gammayl = gammay_inner_core(i,j,k,ispec)
          gammazl = gammaz_inner_core(i,j,k,ispec)

! compute the jacobian
          jacobianl = 1._CUSTOM_REAL / (xixl*(etayl*gammazl-etazl*gammayl) &
                        - xiyl*(etaxl*gammazl-etazl*gammaxl) &
                        + xizl*(etaxl*gammayl-etayl*gammaxl))

          duxdxl = xixl*tempx1l + etaxl*tempx2l + gammaxl*tempx3l
          duxdyl = xiyl*tempx1l + etayl*tempx2l + gammayl*tempx3l
          duxdzl = xizl*tempx1l + etazl*tempx2l + gammazl*tempx3l

          duydxl = xixl*tempy1l + etaxl*tempy2l + gammaxl*tempy3l
          duydyl = xiyl*tempy1l + etayl*tempy2l + gammayl*tempy3l
          duydzl = xizl*tempy1l + etazl*tempy2l + gammazl*tempy3l

          duzdxl = xixl*tempz1l + etaxl*tempz2l + gammaxl*tempz3l
          duzdyl = xiyl*tempz1l + etayl*tempz2l + gammayl*tempz3l
          duzdzl = xizl*tempz1l + etazl*tempz2l + gammazl*tempz3l

! precompute some sums to save CPU time
          duxdxl_plus_duydyl = duxdxl + duydyl
          duxdxl_plus_duzdzl = duxdxl + duzdzl
          duydyl_plus_duzdzl = duydyl + duzdzl
          duxdyl_plus_duydxl = duxdyl + duydxl
          duzdxl_plus_duxdzl = duzdxl + duxdzl
          duzdyl_plus_duydzl = duzdyl + duydzl

! compute deviatoric strain
  if (COMPUTE_AND_STORE_STRAIN) then
    epsilondev_loc(1,i,j,k) = duxdxl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilondev_loc(2,i,j,k) = duydyl - ONE_THIRD * (duxdxl + duydyl + duzdzl)
    epsilondev_loc(3,i,j,k) = 0.5 * duxdyl_plus_duydxl
    epsilondev_loc(4,i,j,k) = 0.5 * duzdxl_plus_duxdzl
    epsilondev_loc(5,i,j,k) = 0.5 * duzdyl_plus_duydzl
  endif

  if(ATTENUATION_VAL) then
! same attenuation everywhere in the inner core therefore no need to use Brian's routines
!!!!!        radius_cr = xstore(ibool_inner_core(i,j,k,ispec))
!!!!!        call get_attenuation_index(idoubling_inner_core(ispec), dble(radius_cr), iregion_selected_inner_core, .TRUE., AM_V)
        minus_sum_beta =  one_minus_sum_beta_inner_core(1,1,1,iregion_selected_inner_core) - 1.0
  endif ! ATTENUATION_VAL

       if(ANISOTROPIC_INNER_CORE_VAL) then

       else

! inner core with no anisotropy, use kappav and muv for instance
! layer with no anisotropy, use kappav and muv for instance
          kappal = kappavstore_inner_core(i,j,k,ispec)
          mul = muvstore_inner_core(i,j,k,ispec)

! use unrelaxed parameters if attenuation
  if(ATTENUATION_VAL) then
      mul = mul * one_minus_sum_beta_inner_core(1,1,1,iregion_selected_inner_core)
  endif

          lambdalplus2mul = kappal + FOUR_THIRDS * mul
          lambdal = lambdalplus2mul - 2.*mul

! compute stress sigma

          sigma_xx = lambdalplus2mul*duxdxl + lambdal*duydyl_plus_duzdzl
          sigma_yy = lambdalplus2mul*duydyl + lambdal*duxdxl_plus_duzdzl
          sigma_zz = lambdalplus2mul*duzdzl + lambdal*duxdxl_plus_duydyl

          sigma_xy = mul*duxdyl_plus_duydxl
          sigma_xz = mul*duzdxl_plus_duxdzl
          sigma_yz = mul*duzdyl_plus_duydzl

        endif

! subtract memory variables if attenuation
  if(ATTENUATION_VAL) then
    do i_sls = 1,N_SLS
      R_xx_val = R_memory_inner_core(1,i_sls,i,j,k,ispec)
      R_yy_val = R_memory_inner_core(2,i_sls,i,j,k,ispec)
      sigma_xx = sigma_xx - R_xx_val
      sigma_yy = sigma_yy - R_yy_val
      sigma_zz = sigma_zz + R_xx_val + R_yy_val
      sigma_xy = sigma_xy - R_memory_inner_core(3,i_sls,i,j,k,ispec)
      sigma_xz = sigma_xz - R_memory_inner_core(4,i_sls,i,j,k,ispec)
      sigma_yz = sigma_yz - R_memory_inner_core(5,i_sls,i,j,k,ispec)
    enddo
  endif

! form dot product with test vector, symmetric form
      tempx1(i,j,k) = jacobianl * (sigma_xx*xixl + sigma_xy*xiyl + sigma_xz*xizl)
      tempy1(i,j,k) = jacobianl * (sigma_xy*xixl + sigma_yy*xiyl + sigma_yz*xizl)
      tempz1(i,j,k) = jacobianl * (sigma_xz*xixl + sigma_yz*xiyl + sigma_zz*xizl)

      tempx2(i,j,k) = jacobianl * (sigma_xx*etaxl + sigma_xy*etayl + sigma_xz*etazl)
      tempy2(i,j,k) = jacobianl * (sigma_xy*etaxl + sigma_yy*etayl + sigma_yz*etazl)
      tempz2(i,j,k) = jacobianl * (sigma_xz*etaxl + sigma_yz*etayl + sigma_zz*etazl)

      tempx3(i,j,k) = jacobianl * (sigma_xx*gammaxl + sigma_xy*gammayl + sigma_xz*gammazl)
      tempy3(i,j,k) = jacobianl * (sigma_xy*gammaxl + sigma_yy*gammayl + sigma_yz*gammazl)
      tempz3(i,j,k) = jacobianl * (sigma_xz*gammaxl + sigma_yz*gammayl + sigma_zz*gammazl)

          enddo
        enddo
      enddo

    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          tempx1l = 0._CUSTOM_REAL
          tempy1l = 0._CUSTOM_REAL
          tempz1l = 0._CUSTOM_REAL

          tempx2l = 0._CUSTOM_REAL
          tempy2l = 0._CUSTOM_REAL
          tempz2l = 0._CUSTOM_REAL

          tempx3l = 0._CUSTOM_REAL
          tempy3l = 0._CUSTOM_REAL
          tempz3l = 0._CUSTOM_REAL

          do l=1,NGLLX
            fac1 = hprimewgll_xx(l,i)
            tempx1l = tempx1l + tempx1(l,j,k)*fac1
            tempy1l = tempy1l + tempy1(l,j,k)*fac1
            tempz1l = tempz1l + tempz1(l,j,k)*fac1
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLY
            fac2 = hprimewgll_yy(l,j)
            tempx2l = tempx2l + tempx2(i,l,k)*fac2
            tempy2l = tempy2l + tempy2(i,l,k)*fac2
            tempz2l = tempz2l + tempz2(i,l,k)*fac2
!!! can merge these loops because NGLLX = NGLLY = NGLLZ          enddo

!!! can merge these loops because NGLLX = NGLLY = NGLLZ          do l=1,NGLLZ
            fac3 = hprimewgll_zz(l,k)
            tempx3l = tempx3l + tempx3(i,j,l)*fac3
            tempy3l = tempy3l + tempy3(i,j,l)*fac3
            tempz3l = tempz3l + tempz3(i,j,l)*fac3
          enddo

          fac1 = wgllwgll_yz(j,k)
          fac2 = wgllwgll_xz(i,k)
          fac3 = wgllwgll_xy(i,j)

          sum_terms(1,i,j,k) = - (fac1*tempx1l + fac2*tempx2l + fac3*tempx3l)
          sum_terms(2,i,j,k) = - (fac1*tempy1l + fac2*tempy2l + fac3*tempy3l)
          sum_terms(3,i,j,k) = - (fac1*tempz1l + fac2*tempz2l + fac3*tempz3l)

        enddo
      enddo
    enddo

! sum contributions from each element to the global mesh and add gravity terms
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool_inner_core(i,j,k,ispec)
          accel_inner_core(:,iglob) = accel_inner_core(:,iglob) + sum_terms(:,i,j,k)
        enddo
      enddo
    enddo

! update memory variables based upon a Runge-Kutta scheme.
! convention for attenuation:
! term in xx = 1
! term in yy = 2
! term in xy = 3
! term in xz = 4
! term in yz = 5
! term in zz not computed since zero trace
    if(ATTENUATION_VAL) then
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            do i_sls = 1,N_SLS
              do i_memory = 1,5
                R_memory_inner_core(i_memory,i_sls,i,j,k,ispec) = &
                  alphaval(i_sls) * &
                  R_memory_inner_core(i_memory,i_sls,i,j,k,ispec) + muvstore_inner_core(i,j,k,ispec) * &
                  factor_common_inner_core(i_sls,1,1,1,iregion_selected_inner_core) * &
                  (betaval(i_sls) * &
                epsilondev_inner_core(i_memory,i,j,k,ispec) + gammaval(i_sls) * epsilondev_loc(i_memory,i,j,k))
              enddo
            enddo
          enddo
        enddo
      enddo
    endif

! save deviatoric strain for Runge-Kutta scheme
    if(COMPUTE_AND_STORE_STRAIN) epsilondev_inner_core(:,:,:,:,ispec) = epsilondev_loc(:,:,:,:)

  enddo ! spectral element loop

  end subroutine compute_forces_CM_IC

