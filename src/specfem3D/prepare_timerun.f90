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

  subroutine prepare_timerun_rmass(myrank,rmass_ocean_load,rmassx_crust_mantle,rmassy_crust_mantle, &
                      rmassz_crust_mantle,rmass_outer_core,rmass_inner_core, &
                      iproc_xi,iproc_eta,ichunk,addressing, &
                      iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle, &
                      iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
                      npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
                      iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
                      iboolleft_xi_outer_core,iboolright_xi_outer_core, &
                      iboolleft_eta_outer_core,iboolright_eta_outer_core, &
                      npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
                      iboolfaces_outer_core,iboolcorner_outer_core, &
                      iboolleft_xi_inner_core,iboolright_xi_inner_core, &
                      iboolleft_eta_inner_core,iboolright_eta_inner_core, &
                      npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
                      iboolfaces_inner_core,iboolcorner_inner_core, &
                      iprocfrom_faces,iprocto_faces,imsg_type, &
                      iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
                      buffer_send_faces_scalar,buffer_received_faces_scalar, &
                      buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
                      NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS,ABSORBING_CONDITIONS, &
                      NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,npoin2D_max_all_CM_IC, &
                      SIMULATION_TYPE,EXACT_MASS_MATRIX_FOR_ROTATION,USE_LDDRK,&
                      b_rmassx_crust_mantle,b_rmassy_crust_mantle,&
                      rmassx_inner_core,rmassy_inner_core,&
                      b_rmassx_inner_core,b_rmassy_inner_core)

  use mpi
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,npoin2D_max_all_CM_IC
  integer SIMULATION_TYPE

  logical ABSORBING_CONDITIONS,EXACT_MASS_MATRIX_FOR_ROTATION,USE_LDDRK

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE_OCEANS) :: rmass_ocean_load
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_CM) :: rmassx_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_CM) :: rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_CM_BACKWARD) :: b_rmassx_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_CM_BACKWARD) :: b_rmassy_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: rmassz_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: rmass_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_IC) :: rmassx_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_IC) :: rmassy_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_IC_BACKWARD) :: b_rmassx_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_XY_IC_BACKWARD) :: b_rmassy_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: rmass_inner_core

  integer ichunk,iproc_xi,iproc_eta
  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing

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

  integer, dimension(NGLOB2DMAX_XY_CM_VAL,NUMFACES_SHARED) :: iboolfaces_crust_mantle
  integer, dimension(NGLOB2DMAX_XY_OC_VAL,NUMFACES_SHARED) :: iboolfaces_outer_core
  integer, dimension(NGLOB2DMAX_XY_IC_VAL,NUMFACES_SHARED) :: iboolfaces_inner_core

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

  ! buffers for send and receive between faces of the slices and the chunks
  real(kind=CUSTOM_REAL), dimension(NGLOB2DMAX_XY_CM_VAL) ::  &
    buffer_send_faces_scalar,buffer_received_faces_scalar

  ! buffers for send and receive between corners of the chunks
  real(kind=CUSTOM_REAL), dimension(NGLOB1D_RADIAL_CM) :: &
    buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar

  integer NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS

  integer, dimension(MAX_NUM_REGIONS) :: NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX

  ! local parameters
  integer :: ier

  ! synchronize all the processes before assembling the mass matrix
  ! to make sure all the nodes have finished to read their databases
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

  ! the mass matrix needs to be assembled with MPI here once and for all

  ! ocean load
  if (OCEANS_VAL) then
    call assemble_MPI_scalar_block(myrank,rmass_ocean_load,NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XY_CM_VAL,NCHUNKS_VAL)
  endif

  ! crust and mantle
  if((NCHUNKS_VAL /= 6 .and. ABSORBING_CONDITIONS .or. &
     (ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION)) .and. NGLOB_CRUST_MANTLE > 0 &
      .and. (.not. USE_LDDRK)) then
     call assemble_MPI_scalar_block(myrank,rmassx_crust_mantle,NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XY_CM_VAL,NCHUNKS_VAL)

     call assemble_MPI_scalar_block(myrank,rmassy_crust_mantle,NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XY_CM_VAL,NCHUNKS_VAL)
  endif

  call assemble_MPI_scalar_block(myrank,rmassz_crust_mantle,NGLOB_CRUST_MANTLE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XY_CM_VAL,NCHUNKS_VAL)

 if(ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION &
    .and. (.not. USE_LDDRK) .and. NGLOB_XY_CM_BACKWARD > 0)then
    if(SIMULATION_TYPE == 3)then
       call assemble_MPI_scalar_block(myrank,b_rmassx_crust_mantle,NGLOB_XY_CM_BACKWARD, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XY_CM_VAL,NCHUNKS_VAL)
      call assemble_MPI_scalar_block(myrank,b_rmassy_crust_mantle,NGLOB_XY_CM_BACKWARD, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
            npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
            iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
            NGLOB2DMAX_XY_CM_VAL,NCHUNKS_VAL)
    endif
 endif

  ! outer core
  call assemble_MPI_scalar_block(myrank,rmass_outer_core,NGLOB_OUTER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_outer_core,iboolright_xi_outer_core,iboolleft_eta_outer_core,iboolright_eta_outer_core, &
            npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
            iboolfaces_outer_core,iboolcorner_outer_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
            NGLOB2DMAX_XY_OC_VAL,NCHUNKS_VAL)

  if(ROTATION_VAL .and. EXACT_MASS_MATRIX_FOR_ROTATION &
     .and. (.not. USE_LDDRK) .and. NGLOB_XY_IC > 0)then
    ! inner core
    call assemble_MPI_scalar_block(myrank,rmassx_inner_core,NGLOB_XY_IC, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
            NGLOB2DMAX_XY_IC_VAL,NCHUNKS_VAL)

   call assemble_MPI_scalar_block(myrank,rmassy_inner_core,NGLOB_XY_IC, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
            NGLOB2DMAX_XY_IC_VAL,NCHUNKS_VAL)

    if(SIMULATION_TYPE == 3  .and. NGLOB_XY_IC_BACKWARD > 0)then
      call assemble_MPI_scalar_block(myrank,b_rmassx_inner_core,NGLOB_XY_IC_BACKWARD, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
            NGLOB2DMAX_XY_IC_VAL,NCHUNKS_VAL)

      call assemble_MPI_scalar_block(myrank,b_rmassy_inner_core,NGLOB_XY_IC_BACKWARD, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
            NGLOB2DMAX_XY_IC_VAL,NCHUNKS_VAL)

    endif

  endif

  ! inner core
  call assemble_MPI_scalar_block(myrank,rmass_inner_core,NGLOB_INNER_CORE, &
            iproc_xi,iproc_eta,ichunk,addressing, &
            iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
            npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
            iboolfaces_inner_core,iboolcorner_inner_core, &
            iprocfrom_faces,iprocto_faces,imsg_type, &
            iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
            buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
            buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
            NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
            NPROC_XI_VAL,NPROC_ETA_VAL,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
            NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
            NGLOB2DMAX_XY_IC_VAL,NCHUNKS_VAL)

  if(myrank == 0) write(IMAIN,*) 'end assembling MPI mass matrix'

  end subroutine prepare_timerun_rmass

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_centralcube(myrank,rmass_inner_core, &
                      iproc_xi,iproc_eta,ichunk, &
                      NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM, &
                      addressing,ibool_inner_core,idoubling_inner_core, &
                      xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                      nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
                      nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
                      ibelm_xmin_inner_core,ibelm_xmax_inner_core, &
                      ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
                      nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube, &
                      npoin2D_cube_from_slices,receiver_cube_from_slices, &
                      sender_from_slices_to_cube,ibool_central_cube, &
                      buffer_slices,buffer_slices2,buffer_all_cube_from_slices)

  use mpi

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: rmass_inner_core

  integer ichunk,iproc_xi,iproc_eta

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM
  integer, dimension(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1) :: addressing

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: ibool_inner_core
  integer, dimension(NSPEC_INNER_CORE) :: idoubling_inner_core
  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE) :: &
        xstore_inner_core,ystore_inner_core,zstore_inner_core

  integer nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
    nspec2D_ymin_inner_core,nspec2D_ymax_inner_core

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_IC) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_IC) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_IC) :: ibelm_bottom_inner_core

  integer nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube, &
    npoin2D_cube_from_slices,receiver_cube_from_slices

  integer, dimension(non_zero_nb_msgs_theor_in_cube) :: sender_from_slices_to_cube
  integer, dimension(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices) :: ibool_central_cube
  double precision, dimension(npoin2D_cube_from_slices,NDIM) :: buffer_slices,buffer_slices2
  double precision, dimension(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM) :: &
    buffer_all_cube_from_slices

  ! local parameters
  integer :: ndim_assemble

  ! create buffers to assemble with the central cube
  call create_central_cube_buffers(myrank,iproc_xi,iproc_eta,ichunk, &
               NPROC_XI_VAL,NPROC_ETA_VAL,NCHUNKS_VAL,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
               NSPEC2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NSPEC2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
               NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
               addressing,ibool_inner_core,idoubling_inner_core, &
               xstore_inner_core,ystore_inner_core,zstore_inner_core, &
               nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
               nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
               ibelm_xmin_inner_core,ibelm_xmax_inner_core, &
               ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
               nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices, &
               receiver_cube_from_slices,sender_from_slices_to_cube,ibool_central_cube, &
               buffer_slices,buffer_slices2,buffer_all_cube_from_slices)

  if(myrank == 0) write(IMAIN,*) 'done including central cube'

  ! the mass matrix to assemble is a scalar, not a vector
  ndim_assemble = 1

  ! use these buffers to assemble the inner core mass matrix with the central cube
  call assemble_MPI_central_cube_block(ichunk,nb_msgs_theor_in_cube, sender_from_slices_to_cube, &
               npoin2D_cube_from_slices, buffer_all_cube_from_slices, &
               buffer_slices, buffer_slices2, ibool_central_cube, &
               receiver_cube_from_slices, ibool_inner_core, &
               idoubling_inner_core, NSPEC_INNER_CORE, &
               ibelm_bottom_inner_core, NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
               NGLOB_INNER_CORE,rmass_inner_core,ndim_assemble)

  ! suppress fictitious mass matrix elements in central cube
  ! because the slices do not compute all their spectral elements in the cube
  where(rmass_inner_core(:) <= 0.) rmass_inner_core = 1._CUSTOM_REAL

  end subroutine prepare_timerun_centralcube

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_constants(myrank,NSTEP, &
                    DT,t0,scale_t,scale_t_inv,scale_displ,scale_veloc, &
                    deltat,deltatover2,deltatsqover2, &
                    b_deltat,b_deltatover2,b_deltatsqover2, &
                    two_omega_earth,A_array_rotation,B_array_rotation, &
                    b_two_omega_earth, SIMULATION_TYPE)

! precomputes constants for time integration

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,NSTEP

  double precision DT
  double precision t0


  double precision scale_t,scale_t_inv,scale_displ,scale_veloc

  real(kind=CUSTOM_REAL) deltat,deltatover2,deltatsqover2
  real(kind=CUSTOM_REAL) b_deltat,b_deltatover2,b_deltatsqover2

  real(kind=CUSTOM_REAL) two_omega_earth
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_ROTATION) :: &
    A_array_rotation,B_array_rotation

  real(kind=CUSTOM_REAL) b_two_omega_earth

  integer SIMULATION_TYPE

  ! local parameters


  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '           time step: ',sngl(DT),' s'
    write(IMAIN,*) 'number of time steps: ',NSTEP
    write(IMAIN,*) 'total simulated time: ',sngl(((NSTEP-1)*DT-t0)/60.d0),' minutes'
    write(IMAIN,*) 'start time:',sngl(-t0),' seconds'
    write(IMAIN,*)
  endif

  ! define constants for the time integration
  ! scaling to make displacement in meters and velocity in meters per second
  scale_t = ONE/dsqrt(PI*GRAV*RHOAV)
  scale_t_inv = dsqrt(PI*GRAV*RHOAV)

  scale_displ = R_EARTH

  scale_veloc = scale_displ * scale_t_inv

  ! distinguish between single and double precision for reals
  if(CUSTOM_REAL == SIZE_REAL) then
    deltat = sngl(DT*scale_t_inv)
  else
    deltat = DT*scale_t_inv
  endif
  deltatover2 = 0.5d0*deltat
  deltatsqover2 = 0.5d0*deltat*deltat

  if (SIMULATION_TYPE == 3) then
    if(CUSTOM_REAL == SIZE_REAL) then
      b_deltat = - sngl(DT*scale_t_inv)
    else
      b_deltat = - DT*scale_t_inv
    endif
    b_deltatover2 = 0.5d0*b_deltat
    b_deltatsqover2 = 0.5d0*b_deltat*b_deltat
  endif

  ! non-dimensionalized rotation rate of the Earth times two
  if(ROTATION_VAL) then
    ! distinguish between single and double precision for reals
    if (SIMULATION_TYPE == 1) then
      if(CUSTOM_REAL == SIZE_REAL) then
        two_omega_earth = sngl(2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv))
      else
        two_omega_earth = 2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv)
      endif
    else
      if(CUSTOM_REAL == SIZE_REAL) then
        two_omega_earth = - sngl(2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv))
      else
        two_omega_earth = - 2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv)
      endif
    endif

    A_array_rotation = 0._CUSTOM_REAL
    B_array_rotation = 0._CUSTOM_REAL

    if (SIMULATION_TYPE == 3) then
      if(CUSTOM_REAL == SIZE_REAL) then
        b_two_omega_earth = sngl(2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv))
      else
        b_two_omega_earth = 2.d0 * TWO_PI / (HOURS_PER_DAY * 3600.d0 * scale_t_inv)
      endif
    endif
  else
    two_omega_earth = 0._CUSTOM_REAL
    if (SIMULATION_TYPE == 3) b_two_omega_earth = 0._CUSTOM_REAL
  endif

  end subroutine prepare_timerun_constants

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_gravity(myrank, &
                    minus_g_cmb,minus_g_icb, &
                    minus_gravity_table,minus_deriv_gravity_table, &
                    density_table,d_ln_density_dr_table,minus_rho_g_over_kappa_fluid, &
                    ONE_CRUST,RICB,RCMB,RTOPDDOUBLEPRIME, &
                    R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

! precomputes gravity factors

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  real(kind=CUSTOM_REAL) minus_g_cmb,minus_g_icb

  ! lookup table every km for gravity
  double precision, dimension(NRAD_GRAVITY) :: minus_gravity_table, &
    minus_deriv_gravity_table,density_table, &
    d_ln_density_dr_table,minus_rho_g_over_kappa_fluid

  logical ONE_CRUST

  double precision RICB,RCMB,RTOPDDOUBLEPRIME, &
    R80,R220,R400,R600,R670,R771,RMOHO,RMIDDLE_CRUST,ROCEAN

  ! local parameters
  double precision :: rspl_gravity(NR),gspl(NR),gspl2(NR)
  double precision :: radius,radius_km,g,dg
  double precision :: g_cmb_dble,g_icb_dble
  double precision :: rho,drhodr,vp,vs,Qkappa,Qmu
  integer :: int_radius,idoubling,nspl_gravity

  ! store g, rho and dg/dr=dg using normalized radius in lookup table every 100 m
  ! get density and velocity from PREM model using dummy doubling flag
  ! this assumes that the gravity perturbations are small and smooth
  ! and that we can neglect the 3D model and use PREM every 100 m in all cases
  ! this is probably a rather reasonable assumption
  if(GRAVITY_VAL) then
    call make_gravity(nspl_gravity,rspl_gravity,gspl,gspl2,ONE_CRUST)
    do int_radius = 1,NRAD_GRAVITY
      radius = dble(int_radius) / (R_EARTH_KM * 10.d0)
      call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g)

      ! use PREM density profile to calculate gravity (fine for other 1D models)
      idoubling = 0
      call model_prem_iso(myrank,radius,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,.false., &
          ONE_CRUST,.false.,RICB,RCMB,RTOPDDOUBLEPRIME, &
          R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

      dg = 4.0d0*rho - 2.0d0*g/radius

      minus_gravity_table(int_radius) = - g
      minus_deriv_gravity_table(int_radius) = - dg
      density_table(int_radius) = rho
      minus_rho_g_over_kappa_fluid(int_radius) = - g / vp**2
    enddo

    ! make sure fluid array is only assigned in outer core between 1222 and 3478 km
    ! lookup table is defined every 100 m
    do int_radius = 1,NRAD_GRAVITY
      radius_km = dble(int_radius) / 10.d0
      if(radius_km > RCMB/1000.d0 - 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RCMB/1000.d0 - 3.d0)*10.d0))
      if(radius_km < RICB/1000.d0 + 3.d0) &
        minus_rho_g_over_kappa_fluid(int_radius) = minus_rho_g_over_kappa_fluid(nint((RICB/1000.d0 + 3.d0)*10.d0))
    enddo

    ! compute gravity value at CMB and ICB once and for all
    radius = RCMB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_cmb_dble)

    radius = RICB / R_EARTH
    call spline_evaluation(rspl_gravity,gspl,gspl2,nspl_gravity,radius,g_icb_dble)

    ! distinguish between single and double precision for reals
    if(CUSTOM_REAL == SIZE_REAL) then
      minus_g_cmb = sngl(- g_cmb_dble)
      minus_g_icb = sngl(- g_icb_dble)
    else
      minus_g_cmb = - g_cmb_dble
      minus_g_icb = - g_icb_dble
    endif

  else

    ! tabulate d ln(rho)/dr needed for the no gravity fluid potential
    do int_radius = 1,NRAD_GRAVITY
       radius = dble(int_radius) / (R_EARTH_KM * 10.d0)
       idoubling = 0
       call model_prem_iso(myrank,radius,rho,drhodr,vp,vs,Qkappa,Qmu,idoubling,.false., &
           ONE_CRUST,.false.,RICB,RCMB,RTOPDDOUBLEPRIME, &
           R600,R670,R220,R771,R400,R80,RMOHO,RMIDDLE_CRUST,ROCEAN)

       d_ln_density_dr_table(int_radius) = drhodr/rho

    enddo

  endif

  end subroutine prepare_timerun_gravity

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_attenuation(myrank, &
                factor_scale_crust_mantle,one_minus_sum_beta_crust_mantle,factor_common_crust_mantle, &
                factor_scale_inner_core,one_minus_sum_beta_inner_core,factor_common_inner_core, &
                c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
                c22store_crust_mantle,c23store_crust_mantle, &
                c33store_crust_mantle,c44store_crust_mantle, &
                c55store_crust_mantle,c66store_crust_mantle, &
                muvstore_crust_mantle,muhstore_crust_mantle,ispec_is_tiso_crust_mantle, &
                muvstore_inner_core, &
                SIMULATION_TYPE,MOVIE_VOLUME,muvstore_crust_mantle_3dmovie, &
                c11store_inner_core,c12store_inner_core,c13store_inner_core, &
                c33store_inner_core,c44store_inner_core, &
                alphaval,betaval,gammaval,b_alphaval,b_betaval,b_gammaval, &
                deltat,b_deltat,LOCAL_PATH,tau_sigma_dble)

  ! precomputes attenuation factors

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank

  ! memory variables and standard linear solids for attenuation
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL) :: &
         one_minus_sum_beta_crust_mantle, factor_scale_crust_mantle
  real(kind=CUSTOM_REAL), dimension(ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL) :: &
         one_minus_sum_beta_inner_core, factor_scale_inner_core
  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT4_VAL) :: factor_common_crust_mantle
  real(kind=CUSTOM_REAL), dimension(N_SLS,ATT1_VAL,ATT2_VAL,ATT3_VAL,ATT5_VAL) :: factor_common_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_MANTLE) :: &
      c11store_crust_mantle,c12store_crust_mantle,c13store_crust_mantle, &
      c22store_crust_mantle,c23store_crust_mantle, &
      c33store_crust_mantle,c44store_crust_mantle, &
      c55store_crust_mantle,c66store_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ISO_MANTLE) :: &
        muvstore_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_TISO_MANTLE) :: &
        muhstore_crust_mantle
  logical, dimension(NSPEC_CRUST_MANTLE) :: ispec_is_tiso_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE) :: &
        muvstore_inner_core

  integer SIMULATION_TYPE
  logical MOVIE_VOLUME

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_3DMOVIE) :: muvstore_crust_mantle_3dmovie

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPECMAX_ANISO_IC) :: &
        c11store_inner_core,c33store_inner_core,c12store_inner_core, &
        c13store_inner_core,c44store_inner_core

  real(kind=CUSTOM_REAL), dimension(N_SLS) :: alphaval, betaval, gammaval
  real(kind=CUSTOM_REAL), dimension(N_SLS) :: b_alphaval, b_betaval, b_gammaval

  real(kind=CUSTOM_REAL) deltat,b_deltat

  character(len=150) LOCAL_PATH

  ! local parameters
  double precision, dimension(N_SLS) :: alphaval_dble, betaval_dble, gammaval_dble
  double precision, dimension(N_SLS) :: tau_sigma_dble

  double precision :: scale_factor,scale_factor_minus_one
  real(kind=CUSTOM_REAL) :: mul
  integer :: ispec,i,j,k
  character(len=150) :: prname

  ! get and store PREM attenuation model

  ! CRUST_MANTLE ATTENUATION
  call create_name_database(prname, myrank, IREGION_CRUST_MANTLE, LOCAL_PATH)
  call get_attenuation_model_3D(myrank, prname, one_minus_sum_beta_crust_mantle, &
           factor_common_crust_mantle,factor_scale_crust_mantle,tau_sigma_dble,NSPEC_CRUST_MANTLE)

  ! INNER_CORE ATTENUATION
  call create_name_database(prname, myrank, IREGION_INNER_CORE, LOCAL_PATH)
  call get_attenuation_model_3D(myrank, prname, one_minus_sum_beta_inner_core, &
           factor_common_inner_core,factor_scale_inner_core,tau_sigma_dble,NSPEC_INNER_CORE)

  ! if attenuation is on, shift PREM to right frequency
  ! rescale mu in PREM to average frequency for attenuation
  ! the formulas to implement the scaling can be found for instance in
  ! Liu, H. P., Anderson, D. L. and Kanamori, H., Velocity dispersion due to
  ! anelasticity: implications for seismology and mantle composition,
  ! Geophys. J. R. Astron. Soc., vol. 47, pp. 41-58 (1976)
  ! and in Aki, K. and Richards, P. G., Quantitative seismology, theory and methods,
  ! W. H. Freeman, (1980), second edition, sections 5.5 and 5.5.2, eq. (5.81) p. 170.
  ! Beware that in the book of Aki and Richards eq. (5.81) is given for velocities
  ! while we need an equation for "mu" and thus we have an additional factor of 2
  ! in the scaling factor below and in equation (49) of Komatitsch and Tromp, Geophys. J. Int. (2002) 149, 390-412,
  ! because "mu" is related to the square of velocity.

  ! rescale in crust and mantle

  do ispec = 1,NSPEC_CRUST_MANTLE
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          if(ATTENUATION_3D_VAL) then
            scale_factor = factor_scale_crust_mantle(i,j,k,ispec)
          else
            scale_factor = factor_scale_crust_mantle(1,1,1,ispec)
          endif

          if(ANISOTROPIC_3D_MANTLE_VAL) then
            scale_factor_minus_one = scale_factor - 1.d0

            mul = c44store_crust_mantle(i,j,k,ispec)
            c11store_crust_mantle(i,j,k,ispec) = c11store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c12store_crust_mantle(i,j,k,ispec) = c12store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c13store_crust_mantle(i,j,k,ispec) = c13store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c22store_crust_mantle(i,j,k,ispec) = c22store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c23store_crust_mantle(i,j,k,ispec) = c23store_crust_mantle(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c33store_crust_mantle(i,j,k,ispec) = c33store_crust_mantle(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c44store_crust_mantle(i,j,k,ispec) = c44store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
            c55store_crust_mantle(i,j,k,ispec) = c55store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
            c66store_crust_mantle(i,j,k,ispec) = c66store_crust_mantle(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
          else
            if(MOVIE_VOLUME .and. SIMULATION_TYPE==3) then
              ! store the original value of \mu to comput \mu*\eps
              muvstore_crust_mantle_3dmovie(i,j,k,ispec)=muvstore_crust_mantle(i,j,k,ispec)
            endif

            muvstore_crust_mantle(i,j,k,ispec) = muvstore_crust_mantle(i,j,k,ispec) * scale_factor

            ! scales transverse isotropic values for mu_h
            if( ispec_is_tiso_crust_mantle(ispec) ) then
              muhstore_crust_mantle(i,j,k,ispec) = muhstore_crust_mantle(i,j,k,ispec) * scale_factor
            endif
          endif

        enddo
      enddo
    enddo
  enddo ! enddo CRUST MANTLE

  ! rescale in inner core

  do ispec = 1,NSPEC_INNER_CORE
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX

          if(ATTENUATION_3D_VAL) then
            scale_factor_minus_one = factor_scale_inner_core(i,j,k,ispec) - 1.0
          else
            scale_factor_minus_one = factor_scale_inner_core(1,1,1,ispec) - 1.0
          endif

          if(ANISOTROPIC_INNER_CORE_VAL) then
            mul = muvstore_inner_core(i,j,k,ispec)
            c11store_inner_core(i,j,k,ispec) = c11store_inner_core(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c12store_inner_core(i,j,k,ispec) = c12store_inner_core(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c13store_inner_core(i,j,k,ispec) = c13store_inner_core(i,j,k,ispec) &
                    - TWO_THIRDS * scale_factor_minus_one * mul
            c33store_inner_core(i,j,k,ispec) = c33store_inner_core(i,j,k,ispec) &
                    + FOUR_THIRDS * scale_factor_minus_one * mul
            c44store_inner_core(i,j,k,ispec) = c44store_inner_core(i,j,k,ispec) &
                    + scale_factor_minus_one * mul
          endif

          if(ATTENUATION_3D_VAL) then
            muvstore_inner_core(i,j,k,ispec) = muvstore_inner_core(i,j,k,ispec) * factor_scale_inner_core(i,j,k,ispec)
          else
            muvstore_inner_core(i,j,k,ispec) = muvstore_inner_core(i,j,k,ispec) * factor_scale_inner_core(1,1,1,ispec)
          endif

        enddo
      enddo
    enddo
  enddo ! enddo INNER CORE

  ! precompute Runge-Kutta coefficients
  call get_attenuation_memory_values(tau_sigma_dble, deltat, alphaval_dble, betaval_dble, gammaval_dble)
  if(CUSTOM_REAL == SIZE_REAL) then
    alphaval = sngl(alphaval_dble)
    betaval  = sngl(betaval_dble)
    gammaval = sngl(gammaval_dble)
  else
    alphaval = alphaval_dble
    betaval  = betaval_dble
    gammaval = gammaval_dble
  endif

  if (SIMULATION_TYPE == 3) then
   call get_attenuation_memory_values(tau_sigma_dble, b_deltat, alphaval_dble, betaval_dble, gammaval_dble)
   if(CUSTOM_REAL == SIZE_REAL) then
     b_alphaval = sngl(alphaval_dble)
     b_betaval  = sngl(betaval_dble)
     b_gammaval = sngl(gammaval_dble)
   else
     b_alphaval = alphaval_dble
     b_betaval  = betaval_dble
     b_gammaval = gammaval_dble
   endif
  endif

  end subroutine prepare_timerun_attenuation
