!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
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

  subroutine compute_add_sources(myrank,NSOURCES, &
                                accel_crust_mantle,sourcearrays, &
                                DT,t0,t_cmt,hdur_gaussian,ibool_crust_mantle, &
                                islice_selected_source,ispec_selected_source,it)

  implicit none
  
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,NSOURCES

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    accel_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES) :: sourcearrays

  double precision, dimension(NSOURCES) :: t_cmt,hdur_gaussian

  double precision :: DT,t0

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source  
  integer :: it
  
  ! local parameters
  double precision :: stf
  real(kind=CUSTOM_REAL) :: stf_used
  integer :: isource,i,j,k,iglob,ispec
  double precision, external :: comp_source_time_function

  do isource = 1,NSOURCES

    
    ! add only if this proc carries the source
    if(myrank == islice_selected_source(isource)) then

      stf = comp_source_time_function(dble(it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

      !     distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        stf_used = sngl(stf)
      else
        stf_used = stf
      endif

      !     add source array
      ispec = ispec_selected_source(isource)      
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            
            accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
              + sourcearrays(:,i,j,k,isource)*stf_used
              
          enddo
        enddo
      enddo

    endif

  enddo
  
  end subroutine compute_add_sources
  
!
!-------------------------------------------------------------------------------------------------
!
  
  subroutine compute_add_sources_adjoint(myrank,nrec, &
                                nadj_rec_local,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC, &
                                accel_crust_mantle,adj_sourcearrays, &
                                nu,xi_receiver,eta_receiver,gamma_receiver, &
                                xigll,yigll,zigll,ibool_crust_mantle, &
                                islice_selected_rec,ispec_selected_rec, &
                                NSTEP_SUB_ADJ,iadjsrc_len,iadjsrc,iadj_vec, &
                                it,it_begin,station_name,network_name)

  implicit none
  
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,nrec,nadj_rec_local,NSTEP,NTSTEP_BETWEEN_READ_ADJSRC

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    accel_crust_mantle

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,NGLLZ,nadj_rec_local,NTSTEP_BETWEEN_READ_ADJSRC) :: &
    adj_sourcearrays
  
  double precision, dimension(NDIM,NDIM,nrec) :: nu
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec

  integer NSTEP_SUB_ADJ
  integer, dimension(NSTEP_SUB_ADJ) :: iadjsrc_len
  integer, dimension(NSTEP_SUB_ADJ,2) :: iadjsrc ! to read input in chunks
  integer, dimension(NSTEP) :: iadj_vec
  
  integer :: it,it_begin,itime

  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name
      
  ! local parameters  
  real(kind=CUSTOM_REAL), dimension(:,:,:,:,:), allocatable :: adj_sourcearray
  integer :: irec,irec_local,i,j,k,iglob,it_sub_adj
  character(len=150) :: adj_source_file
  logical :: ibool_read_adj_arrays
  
  ! figure out if we need to read in a chunk of the adjoint source at this timestep
  it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )   !chunk_number
  ibool_read_adj_arrays = (((it == it_begin) .or. (mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0)) &
                                                  .and. (nadj_rec_local > 0))

  ! needs to read in a new chunk/block of the adjoint source
  if(ibool_read_adj_arrays) then

    ! temporary source array
    allocate(adj_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NTSTEP_BETWEEN_READ_ADJSRC))
    adj_sourcearray = 0._CUSTOM_REAL
    
    irec_local = 0
    do irec = 1, nrec
      ! check that the source slice number is okay
      if(islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROCTOT_VAL-1) then
        if(islice_selected_rec(irec) < 0) call exit_MPI(myrank,'islice < 0')
        if(islice_selected_rec(irec) > NPROCTOT_VAL-1) call exit_MPI(myrank,'islice > NPROCTOT_VAL-1')
        call exit_MPI(myrank,'now: something is wrong with the source slice number in adjoint simulation')
      endif
      ! compute source arrays
      if(myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        
        ! reads in **sta**.**net**.**LH**.adj files
        adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
        call compute_arrays_source_adjoint(myrank,adj_source_file, &
                  xi_receiver(irec),eta_receiver(irec),gamma_receiver(irec), &
                  nu(:,:,irec),adj_sourcearray, xigll,yigll,zigll,iadjsrc_len(it_sub_adj), &
                  iadjsrc,it_sub_adj,NSTEP_SUB_ADJ,NTSTEP_BETWEEN_READ_ADJSRC)

        ! stores source array
        ! note: the adj_sourcearrays has a time stepping from 1 to NTSTEP_BETWEEN_READ_ADJSRC
        !          this gets overwritten every time a new block/chunk is read in
        do itime = 1,NTSTEP_BETWEEN_READ_ADJSRC
          adj_sourcearrays(:,:,:,:,irec_local,itime) = adj_sourcearray(:,:,:,:,itime)
        enddo
        
      endif
    enddo
    if(irec_local /= nadj_rec_local) &
      call exit_MPI(myrank,'irec_local /= nadj_rec_local in adjoint simulation')
    
    deallocate(adj_sourcearray)  
    
  endif

  irec_local = 0
  do irec = 1,nrec

    ! adds source (only if this proc carries the source)
    if(myrank == islice_selected_rec(irec)) then
      irec_local = irec_local + 1

      ! adds source contributions
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec_selected_rec(irec))
            
            
            ! adds adjoint source acting at this time step (it):
            !
            ! note: we use index iadj_vec(it) which is the corresponding time step
            !          for the adjoint source acting at this time step (it)
            !
            ! see routine: setup_sources_receivers_adjindx() how this adjoint index array is set up
            !
            !           e.g. total length NSTEP = 3000, chunk length NTSTEP_BETWEEN_READ_ADJSRC= 1000 
            !           then for it=1,..1000, first block has iadjsrc(1,1) with start = 2001 and end = 3000;
            !           corresponding iadj_vec(it) goes from 
            !           iadj_vec(1) = 1000, iadj_vec(2) = 999 to iadj_vec(1000) = 1, 
            !           that is, originally the idea was
            !           adj_sourcearrays(.. iadj_vec(1) ) corresponds to adjoint source trace at time index 3000
            !           adj_sourcearrays(.. iadj_vec(2) ) corresponds to adjoint source trace at time index 2999
            !           ..
            !           adj_sourcearrays(.. iadj_vec(1000) ) corresponds to adjoint source trace at time index 2001
            !           then a new block will be read, etc, and it is going down till to adjoint source trace at time index 1
            !           
            ! now comes the tricky part:
            !           adjoint source traces are based on the seismograms from the forward run;
            !           such seismograms have a time step index 1 which corresponds to time -t0
            !           then time step index 2 which corresponds to -t0 + DT, and
            !           the last time step in the file at time step NSTEP corresponds to time -t0 + (NSTEP-1)*DT
            !           (see how we add the sources to the simulation in compute_add_sources() and 
            !             how we write/save the seismograms and wavefields at the end of the time loop).
            !
            !           then you use that seismogram and take e.g. the velocity of it for a travetime adjoint source
            !
            !           now we read it in again, and remember the last time step in 
            !           the file at NSTEP corresponds to -t0 + (NSTEP-1)*DT
            !           
            !           the same time step is saved for the forward wavefields to reconstruct them;
            !           however, the Newark time scheme acts at the very beginning of this time loop
            !           such that we have the backward/reconstructed wavefield updated by
            !           a single time step into the direction -DT and b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT - DT
            !           after the Newark (predictor) time step update.
            !
            !           for the kernel calculations, we want: 
            !             adjoint wavefield at time t, starting from 0 to T
            !             and forward wavefield at time T-t, starting from T down to 0
            !           let's say time 0 corresponds to -t0 = -t0 + (it - 1)*DT at it=1
            !             and time T corresponds to -t0 + (NSTEP-1)*DT  at it = NSTEP
            !
            !           as seen before, the time for the forward wavefield b_displ(it=1) would then
            !           correspond to time -t0 + (NSTEP-1)*DT - DT, which is T - DT.
            !           the corresponding time for the adjoint wavefield thus would be 0 + DT
            !           and the adjoint source index would be iadj_vec(it+1) 
            !           however, iadj_vec(it+1) which would go from 999 down to 0. 0 is out of bounds.
            !           we thus have to read in the adjoint source trace beginning from 2999 down to 0.
            !           index 0 is not defined in the adjoint source trace, let's just set it to zero.
            !
            !           for this reason, we use iadj_vec(it) from 1000 down to 1, but 
            !           pointing to the indices in the adjoint source trace from 2999 down to 0 (virtually zero). 
            !           i guess by now everyone is confused. see how we read in the adjoint source trace
            !           in routine compute_arrays_source_adjoint() and note that there is a shift of indices.
            !           then the following should be correct...
            accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &         
                          + adj_sourcearrays(:,i,j,k,irec_local,iadj_vec(it))

          enddo
        enddo
      enddo
    endif

  enddo
  
  
  end subroutine compute_add_sources_adjoint

!
!-------------------------------------------------------------------------------------------------
!
  
  
  subroutine compute_add_sources_backward(myrank,NSOURCES,NSTEP, &
                                b_accel_crust_mantle,sourcearrays, &
                                DT,t0,t_cmt,hdur_gaussian,ibool_crust_mantle, &
                                islice_selected_source,ispec_selected_source,it)

  implicit none
  
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,NSOURCES,NSTEP

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_accel_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES) :: sourcearrays

  double precision, dimension(NSOURCES) :: t_cmt,hdur_gaussian

  double precision :: DT,t0

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source  
  integer :: it
  
  ! local parameters
  double precision :: stf
  real(kind=CUSTOM_REAL) :: stf_used
  integer :: isource,i,j,k,iglob,ispec
  double precision, external :: comp_source_time_function

  do isource = 1,NSOURCES

    !   add the source (only if this proc carries the source)
    if(myrank == islice_selected_source(isource)) then

! note on backward/reconstructed wavefields: 
!       time for b_displ( it ) corresponds to (NSTEP - it - 1 )*DT - t0  ...
!       as we start with saved wavefields b_displ( 0 ) = displ( NSTEP ) which correspond
!       to a time (NSTEP - 1)*DT - t0 
!       (see sources for simulation_type 1 and seismograms)
!
!       now, at the beginning of the time loop, the numerical Newark time scheme updates
!       the wavefields, that is b_displ( it=1) corresponds now to time (NSTEP -1 - 1)*DT - t0.
!       this leads to the timing (NSTEP-it-1)*DT-t0-t_cmt for the source time function here
      stf = comp_source_time_function(dble(NSTEP-it-1)*DT-t0-t_cmt(isource),hdur_gaussian(isource))

      !     distinguish between single and double precision for reals
      if(CUSTOM_REAL == SIZE_REAL) then
        stf_used = sngl(stf)
      else
        stf_used = stf
      endif

      !     add source array
      ispec = ispec_selected_source(isource)            
      do k=1,NGLLZ
        do j=1,NGLLY
          do i=1,NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)
            
            b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) &
              + sourcearrays(:,i,j,k,isource)*stf_used

          enddo
        enddo
      enddo

    endif

  enddo
  
  end subroutine compute_add_sources_backward
  
  