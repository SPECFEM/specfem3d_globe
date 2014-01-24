!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

  subroutine compute_add_sources()

  use specfem_par
  use specfem_par_crustmantle,only: accel_crust_mantle,ibool_crust_mantle

  implicit none

  ! local parameters
  double precision :: stf
  real(kind=CUSTOM_REAL) :: stf_used
  integer :: isource,i,j,k,iglob,ispec
  double precision :: f0
  double precision, dimension(NSOURCES) :: stf_pre_compute

  double precision, external :: comp_source_time_function
  double precision, external :: comp_source_time_function_rickr

  if( .not. GPU_MODE ) then
    ! on CPU
    do isource = 1,NSOURCES

      ! add only if this proc carries the source
      if(myrank == islice_selected_source(isource)) then

        if(USE_FORCE_POINT_SOURCE) then

          ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
          iglob = ibool_crust_mantle(nint(xi_source(isource)), &
                         nint(eta_source(isource)), &
                         nint(gamma_source(isource)), &
                         ispec_selected_source(isource))

          f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format

          ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
          if(USE_LDDRK)then
            stf_used = FACTOR_FORCE_SOURCE * &
                     comp_source_time_function_rickr(dble(it-1)*DT + dble(C_LDDRK(istage))*DT-t0-tshift_cmt(isource),f0)
          else
            stf_used = FACTOR_FORCE_SOURCE * comp_source_time_function_rickr(dble(it-1)*DT-t0-tshift_cmt(isource),f0)
          endif

          ! we use a force in a single direction along one of the components:
          !  x/y/z or E/N/Z-direction would correspond to 1/2/3 = COMPONENT_FORCE_SOURCE
          ! e.g. nu_source(3,:) here would be a source normal to the surface (z-direction).
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob)  &
                           + sngl( nu_source(COMPONENT_FORCE_SOURCE,:,isource) ) * stf_used

        else
          if(USE_LDDRK)then
            stf = comp_source_time_function(dble(it-1)*DT + &
                                            dble(C_LDDRK(istage))*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))
          else
            stf = comp_source_time_function(dble(it-1)*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))
          endif

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

        endif ! USE_FORCE_POINT_SOURCE

      endif

    enddo

  else
    ! on GPU
    ! prepares buffer with source time function values, to be copied onto GPU
    if(USE_FORCE_POINT_SOURCE) then
      do isource = 1,NSOURCES
        if(USE_LDDRK)then
          stf_pre_compute(isource) = FACTOR_FORCE_SOURCE * &
                     comp_source_time_function_rickr(dble(it-1)*DT + dble(C_LDDRK(istage))*DT-t0-tshift_cmt(isource),f0)
        else
          stf_pre_compute(isource) = &
                      FACTOR_FORCE_SOURCE * comp_source_time_function_rickr(dble(it-1)*DT-t0-tshift_cmt(isource),f0)
        endif
      enddo
    else
      do isource = 1,NSOURCES
        if(USE_LDDRK)then
          stf_pre_compute(isource) = comp_source_time_function(dble(it-1)*DT + &
                                            dble(C_LDDRK(istage))*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))
        else
          stf_pre_compute(isource) = &
            comp_source_time_function(dble(it-1)*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))
        endif
      enddo
    endif
    ! adds sources: only implements SIMTYPE=1 and NOISE_TOM=0
    call compute_add_sources_cuda(Mesh_pointer,NSOURCES,stf_pre_compute)
  endif


  end subroutine compute_add_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_adjoint()

  use specfem_par
  use specfem_par_crustmantle,only: accel_crust_mantle,ibool_crust_mantle

  implicit none

  ! local parameters
  integer :: irec,irec_local,i,j,k,iglob
  integer :: ivec_index
  logical :: ibool_read_adj_arrays

  ! note: we check if nadj_rec_local > 0 before calling this routine, but better be safe...
  if( nadj_rec_local == 0 ) return

  ! figure out if we need to read in a chunk of the adjoint source at this timestep
  ibool_read_adj_arrays = ( (it == it_begin) .or. (mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0) )

  ! needs to read in a new chunk/block of the adjoint source
  if(ibool_read_adj_arrays) then
    call read_adjoint_sources()
  endif

  ! adds adjoint sources
  if( .not. GPU_MODE ) then
    ! on CPU
    irec_local = 0
    do irec = 1,nrec

      ! adds source (only if this proc carries the source)
      if(myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1

        ! adjoint source array index
        ivec_index = iadj_vec(it)

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
              !           however, the Newmark time scheme acts at the very beginning of this time loop
              !           such that we have the backward/reconstructed wavefield updated by
              !           a single time step into the direction -DT and b_displ(it=1) would  corresponds to -t0 + (NSTEP-1)*DT - DT
              !           after the Newmark (predictor) time step update.
              !           however, we will read the backward/reconstructed wavefield at the end of the first time loop,
              !           such that b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT (which is the one saved in the files).
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
              !           we thus would have to read in the adjoint source trace beginning from 2999 down to 0.
              !           index 0 is not defined in the adjoint source trace, and would be set to zero.
              !
              !           however, since this complicates things, we read the backward/reconstructed
              !           wavefield at the end of the first time loop, such that b_displ(it=1) corresponds to -t0 + (NSTEP-1)*DT.
              !           assuming that until that end the backward/reconstructed wavefield and adjoint fields
              !           have a zero contribution to adjoint kernels.
              accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
                            + adj_sourcearrays(:,i,j,k,irec_local,ivec_index)

            enddo
          enddo
        enddo
      endif

    enddo

  else

    ! on GPU
    ! note: adjoint sourcearrays can become very big when used with many receiver stations
    !       we overlap here the memory transfer to GPUs

    ! current time index
    ivec_index = iadj_vec(it)

    if( GPU_ASYNC_COPY ) then
      ! only synchronuously transfers array at beginning or whenever new arrays were read in
      if( ibool_read_adj_arrays ) then
        ! transfers adjoint arrays to GPU device memory
        ! note: function call passes pointer to array adj_sourcearrays at corresponding time slice
        call transfer_adj_to_device(Mesh_pointer,nrec,adj_sourcearrays(1,1,1,1,1,ivec_index), &
                                    islice_selected_rec)
      endif
    else
      ! synchronuously transfers adjoint arrays to GPU device memory before adding adjoint sources on GPU
      call transfer_adj_to_device(Mesh_pointer,nrec,adj_sourcearrays(1,1,1,1,1,ivec_index), &
                                  islice_selected_rec)
    endif

    ! adds adjoint source contributions
    call compute_add_sources_adjoint_cuda(Mesh_pointer,nrec)

    if( GPU_ASYNC_COPY ) then
      ! starts asynchronuously transfer of next adjoint arrays to GPU device memory
      ! (making sure the next adj_sourcearrays values were already read in)
      if( (.not. ibool_read_adj_arrays) .and. &
          (.not. mod(it,NTSTEP_BETWEEN_READ_ADJSRC) == 0) .and. &
          (.not. it == it_end) ) then
        ! next time index
        ivec_index = iadj_vec(it+1)

        ! checks next index
        if( ivec_index < 1 .or. ivec_index > NTSTEP_BETWEEN_READ_ADJSRC ) then
          print*,'error iadj_vec bounds: rank',myrank,' it = ',it,' index = ',ivec_index, &
                 'out of bounds ',1,'to',NTSTEP_BETWEEN_READ_ADJSRC
          call exit_MPI(myrank,'error iadj_vec index bounds')
        endif

        ! asynchronuously transfers next time slice
        call transfer_adj_to_device_async(Mesh_pointer,nrec,adj_sourcearrays(1,1,1,1,1,ivec_index), &
                                          islice_selected_rec)
      endif
    endif

  endif

  end subroutine compute_add_sources_adjoint

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_backward()

  use specfem_par
  use specfem_par_crustmantle,only: b_accel_crust_mantle,ibool_crust_mantle
  implicit none

  ! local parameters
  double precision :: stf
  real(kind=CUSTOM_REAL) :: stf_used
  integer :: isource,i,j,k,iglob,ispec
  double precision :: f0
  double precision, dimension(NSOURCES) :: stf_pre_compute

  double precision, external :: comp_source_time_function
  double precision, external :: comp_source_time_function_rickr

  integer :: it_tmp

  ! iteration step
  if( UNDO_ATTENUATION ) then
    it_tmp = iteration_on_subset * NT_DUMP_ATTENUATION - it_of_this_subset + 1
  else
    it_tmp = it
  endif

  if( .not. GPU_MODE ) then
    ! on CPU
    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if(myrank == islice_selected_source(isource)) then

        ! note on backward/reconstructed wavefields:
        !       time for b_displ( it ) corresponds to (NSTEP - (it-1) - 1 )*DT - t0  ...
        !       as we start with saved wavefields b_displ( 1 ) = displ( NSTEP ) which correspond
        !       to a time (NSTEP - 1)*DT - t0
        !       (see sources for simulation_type 1 and seismograms)
        !
        !       now, at the beginning of the time loop, the numerical Newmark time scheme updates
        !       the wavefields, that is b_displ( it=1) would correspond to time (NSTEP -1 - 1)*DT - t0.
        !       however, we read in the backward/reconstructed wavefields at the end of the Newmark time scheme
        !       in the first (it=1) time loop.
        !       this leads to the timing (NSTEP-(it-1)-1)*DT-t0-tshift_cmt for the source time function here

        if(USE_FORCE_POINT_SOURCE) then

           ! note: for use_force_point_source xi/eta/gamma are in the range [1,NGLL*]
           iglob = ibool_crust_mantle(nint(xi_source(isource)), &
                         nint(eta_source(isource)), &
                         nint(gamma_source(isource)), &
                         ispec_selected_source(isource))

           f0 = hdur(isource) !! using hdur as a FREQUENCY just to avoid changing CMTSOLUTION file format

           !if (it == 1 .and. myrank == 0) then
           !   write(IMAIN,*) 'using a source of dominant frequency ',f0
           !   write(IMAIN,*) 'lambda_S at dominant frequency = ',3000./sqrt(3.)/f0
           !   write(IMAIN,*) 'lambda_S at highest significant frequency = ',3000./sqrt(3.)/(2.5*f0)
           !endif

           ! This is the expression of a Ricker; should be changed according maybe to the Par_file.
           stf_used = FACTOR_FORCE_SOURCE * comp_source_time_function_rickr(dble(NSTEP-it_tmp)*DT-t0-tshift_cmt(isource),f0)

           ! e.g. we use nu_source(3,:) here if we want a source normal to the surface.
           ! note: time step is now at NSTEP-it
           b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob)  &
                              + sngl( nu_source(COMPONENT_FORCE_SOURCE,:,isource) ) * stf_used

        else

          ! see note above: time step corresponds now to NSTEP-it
          stf = comp_source_time_function(dble(NSTEP-it_tmp)*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))

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

        endif ! USE_FORCE_POINT_SOURCE

      endif

    enddo

  else
    ! on GPU
    ! prepares buffer with source time function values, to be copied onto GPU
    if(USE_FORCE_POINT_SOURCE) then
      do isource = 1,NSOURCES
        stf_pre_compute(isource) = &
          FACTOR_FORCE_SOURCE * comp_source_time_function_rickr(dble(NSTEP-it_tmp)*DT-t0-tshift_cmt(isource),f0)
      enddo
    else
      do isource = 1,NSOURCES
        stf_pre_compute(isource) = &
          comp_source_time_function(dble(NSTEP-it_tmp)*DT-t0-tshift_cmt(isource),hdur_gaussian(isource))
      enddo
    endif
    ! adds sources: only implements SIMTYPE=3 (and NOISE_TOM=0)
    call compute_add_sources_backward_cuda(Mesh_pointer,NSOURCES,stf_pre_compute)
  endif

  end subroutine compute_add_sources_backward
