!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  8 . 0
!          --------------------------------------------------
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

  subroutine compute_add_sources()

  use specfem_par
  use specfem_par_crustmantle, only: accel_crust_mantle,ibool_crust_mantle

  implicit none

  ! local parameters
  integer :: isource,i,j,k,iglob,ispec
  double precision :: timeval,time_t
  double precision :: stf
  real(kind=CUSTOM_REAL) :: stf_used
  ! for gpu
  double precision, dimension(NSOURCES) :: stf_pre_compute

  double precision, external :: get_stf_viscoelastic

  ! checks if anything to do for noise simulation
  if (NOISE_TOMOGRAPHY /= 0) return

  ! sets current initial time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
    time_t = dble(it-1-1)*DT + dble(C_LDDRK(istage))*DT - t0
  else
    time_t = dble(it-1)*DT - t0
  endif

  if (.not. GPU_MODE) then
    ! on CPU
! openmp solver
!$OMP PARALLEL if (NSOURCES > 100) &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(isource,timeval,iglob,stf_used,stf,ispec,i,j,k)
!$OMP DO
    do isource = 1,NSOURCES

      ! add only if this proc carries the source
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        ! sets current time for this source
        timeval = time_t - tshift_src(isource)

        ! determines source time function value
        stf = get_stf_viscoelastic(timeval,isource)

        ! distinguishes between single and double precision for reals
        stf_used = real(stf,kind=CUSTOM_REAL)

        ! adds source contribution
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool_crust_mantle(i,j,k,ispec)
!$OMP ATOMIC
              accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + sourcearrays(1,i,j,k,isource)*stf_used
!$OMP ATOMIC
              accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + sourcearrays(2,i,j,k,isource)*stf_used
!$OMP ATOMIC
              accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + sourcearrays(3,i,j,k,isource)*stf_used
            enddo
          enddo
        enddo

      endif

    enddo
!$OMP ENDDO
!$OMP END PARALLEL

  else
    ! on GPU
    ! prepares buffer with source time function values, to be copied onto GPU
    do isource = 1,NSOURCES
      ! sets current time for this source
      timeval = time_t - tshift_src(isource)

      ! determines source time function value
      stf = get_stf_viscoelastic(timeval,isource)

      ! stores current stf values
      stf_pre_compute(isource) = stf
    enddo

    ! adds sources: only implements SIMTYPE=1 and NOISE_TOM = 0
    call compute_add_sources_gpu(Mesh_pointer,NSOURCES,stf_pre_compute)
  endif

  end subroutine compute_add_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_adjoint()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ,NDIM

  use specfem_par, only: myrank,it,it_begin,it_end,NTSTEP_BETWEEN_READ_ADJSRC, &
    nadj_rec_local,hxir_adjstore,hetar_adjstore,hgammar_adjstore,number_adjsources_global, &
    islice_selected_rec,ispec_selected_rec,nrec,iadj_vec, &
    GPU_MODE,GPU_ASYNC_COPY,Mesh_pointer

  use specfem_par_crustmantle, only: accel_crust_mantle,ibool_crust_mantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL),dimension(NDIM) :: stf_array
  real(kind=CUSTOM_REAL) :: hlagrange
  real(kind=CUSTOM_REAL),dimension(NDIM,nadj_rec_local) :: stf_array_adjoint

  integer :: irec,irec_local,i,j,k,iglob,ispec
  integer :: ivec_index
  logical :: ibool_read_adj_arrays

  ! note: we check if nadj_rec_local > 0 before calling this routine, but better be safe...
  if (nadj_rec_local == 0) return

  ! figure out if we need to read in a chunk of the adjoint source at this timestep
  ibool_read_adj_arrays = ( (it == it_begin) .or. (mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0) )

  ! needs to read in a new chunk/block of the adjoint source
  if (ibool_read_adj_arrays) then
    call read_adjoint_sources()
  endif

  ! adds adjoint sources
  if (.not. GPU_MODE) then
    ! on CPU

! work-around a cray compiler issue with the loop below.
! the issue occurs when running debugging flag -g without optimization specifier (-O0). the compiler still tries to
! optimize the routines, but runs into the following internal error:
!
!*** Optimization assertion failure:
!   'Analyze_aliases'
!
!   Error detected     ::  File 'pdgcs/v_df.c', line 8282
!   Initiated from     ::  Line 1562 (v_main.c)
!   Optimizer built    ::  2017-12-05 (production)
!
!   File               ::  src/specfem3D/compute_add_sources.f90
!   function           ::  compute_add_sources_adjoint
!   at or near line    ::  229
!
!problematic part:
!    irec_local = 0
!    do irec = 1,nrec
!
!      ! adds source (only if this proc carries the source)
!      if (myrank == islice_selected_rec(irec)) then
!        irec_local = irec_local + 1
!        do ..
!         ..
!        enddo
!      endif
!    enddo <-- this is the line mentioned in the error
!
! the work-around gets first all local receivers and then loops only over
! these local ones (without the need of the if-case)
!
! note: cray compilation with -g but without -O0 still fails for routines in compute_stacey_*.f90
!       cray version 8.6.x needs -O0 when debugging flags are used. as a work around, the following
!       debugging flags work for crayftn: -g -G0 -O0 -Rb -eF -rm -eC -eD -ec -en -eI -ea

    ! receivers act as sources
    do irec_local = 1, nadj_rec_local
      ! adjoint source time function (trace)
      call get_stf_adjoint_source(it,irec_local,stf_array)

      ! receiver location
      irec = number_adjsources_global(irec_local)

      ! element index
      ispec = ispec_selected_rec(irec)

      ! adds source contributions
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool_crust_mantle(i,j,k,ispec)

            hlagrange = hxir_adjstore(i,irec_local) * hetar_adjstore(j,irec_local) * hgammar_adjstore(k,irec_local)

            ! adds adjoint source acting at this time step (it):
            !
            ! note: we use index iadj_vec(it) which is the corresponding time step
            !          for the adjoint source acting at this time step (it)
            !
            ! see routine: setup_sources_receivers_adjindx() how this adjoint index array is set up
            !
            !           e.g. total length NSTEP = 3000, chunk length NTSTEP_BETWEEN_READ_ADJSRC= 1000
            !           then for it = 1,..1000, first block has iadjsrc(1,1) with start = 2001 and end = 3000;
            !           corresponding iadj_vec(it) goes from
            !           iadj_vec(1) = 1000, iadj_vec(2) = 999 to iadj_vec(1000) = 1,
            !           that is, originally the idea was
            !           source_adjoint(.. iadj_vec(1) ) corresponds to adjoint source trace at time index 3000
            !           source_adjoint(.. iadj_vec(2) ) corresponds to adjoint source trace at time index 2999
            !           ..
            !           source_adjoint(.. iadj_vec(1000) ) corresponds to adjoint source trace at time index 2001
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
            !           then you use that seismogram and take e.g. the velocity of it for a traveltime adjoint source
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
            accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) + stf_array(:) * hlagrange

          enddo ! NGLLX
        enddo ! NGLLY
      enddo ! NGLLZ

    enddo ! irec_local

  else

    ! on GPU
    ! note: adjoint sourcearrays can become very big when used with many receiver stations
    !       we overlap here the memory transfer to GPUs

    ! receivers act as sources
    ! determines STF contribution for each local adjoint source (taking account of LDDRK scheme stages)
    do irec_local = 1, nadj_rec_local
      ! adjoint source time function (trace)
      call get_stf_adjoint_source(it,irec_local,stf_array)
      ! stores stf for all local adjoint sources
      stf_array_adjoint(:,irec_local) = stf_array(:)
    enddo

    if (GPU_ASYNC_COPY) then
      ! only synchronously transfers array at beginning or whenever new arrays were read in
      if (ibool_read_adj_arrays) then
        ! transfers adjoint arrays to GPU device memory
        ! note: function call passes pointer to array source_adjoint at corresponding time slice
        call transfer_adj_to_device(Mesh_pointer,nrec,stf_array_adjoint,islice_selected_rec)
      endif
    else
      ! synchronously transfers adjoint arrays to GPU device memory before adding adjoint sources on GPU
      call transfer_adj_to_device(Mesh_pointer,nrec,stf_array_adjoint,islice_selected_rec)
    endif

    ! adds adjoint source contributions
    call compute_add_sources_adjoint_gpu(Mesh_pointer)

    if (GPU_ASYNC_COPY) then
      ! starts asynchronously transfer of next adjoint arrays to GPU device memory
      ! (making sure the next source_adjoint values were already read in)
      if ((.not. ibool_read_adj_arrays) .and. &
          (.not. mod(it,NTSTEP_BETWEEN_READ_ADJSRC) == 0) .and. &
          (.not. it == it_end)) then

        ! checks next index
        ivec_index = iadj_vec(it+1)
        if (ivec_index < 1 .or. ivec_index > NTSTEP_BETWEEN_READ_ADJSRC) then
          print *,'Error iadj_vec bounds: rank',myrank,' it = ',it,' index = ',ivec_index, &
                  'out of bounds ',1,'to',NTSTEP_BETWEEN_READ_ADJSRC
          call exit_MPI(myrank,'Error iadj_vec index bounds')
        endif

        ! next time index
        do irec_local = 1, nadj_rec_local
          ! adjoint source time function (trace)
          call get_stf_adjoint_source(it+1,irec_local,stf_array)
          ! stores stf for all local adjoint sources
          stf_array_adjoint(:,irec_local) = stf_array(:)
        enddo

        ! asynchronously transfers next time slice
        call transfer_adj_to_device_async(Mesh_pointer,nrec,stf_array_adjoint,islice_selected_rec)
      endif
    endif

  endif

  end subroutine compute_add_sources_adjoint

!
!-------------------------------------------------------------------------------------------------
!

  subroutine compute_add_sources_backward()

  use specfem_par
  use specfem_par_crustmantle, only: b_accel_crust_mantle,ibool_crust_mantle
  implicit none

  ! local parameters
  integer :: isource,i,j,k,iglob,ispec
  integer :: it_tmp
  double precision :: timeval,time_t
  double precision :: stf
  real(kind=CUSTOM_REAL) :: stf_used
  ! for gpu
  double precision, dimension(NSOURCES) :: stf_pre_compute

  double precision, external :: get_stf_viscoelastic

  ! checks if anything to do for noise simulation
  if (NOISE_TOMOGRAPHY /= 0) return

  ! iteration step
  if (UNDO_ATTENUATION) then
    ! example: NSTEP is a multiple of NT_DUMP_ATTENUATION
    !         NT_DUMP_ATTENUATION = 301, NSTEP = 1204, NSUBSET_ITERATIONS = 4, iteration_on_subset = 1 -> 4,
    !              1. subset, it_temp goes from 301 down to 1
    !              2. subset, it_temp goes from 602 down to 302
    !              3. subset, it_temp goes from 903 down to 603
    !              4. subset, it_temp goes from 1204 down to 904
    !valid for multiples only:
    !it_tmp = iteration_on_subset * NT_DUMP_ATTENUATION - it_of_this_subset + 1
    !
    ! example: NSTEP is **NOT** a multiple of NT_DUMP_ATTENUATION
    !          NT_DUMP_ATTENUATION = 301, NSTEP = 900, NSUBSET_ITERATIONS = 3, iteration_on_subset = 1 -> 3
    !              1. subset, it_temp goes from (900 - 602) = 298 down to 1
    !              2. subset, it_temp goes from (900 - 301) = 599 down to 299
    !              3. subset, it_temp goes from (900 - 0)   = 900 down to 600
    !works always:
    if (NSTEP_STEADY_STATE > 0) then
      it_tmp = NSTEP_STEADY_STATE - (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION - it_of_this_subset + 1
    else
      it_tmp = NSTEP - (NSUBSET_ITERATIONS - iteration_on_subset)*NT_DUMP_ATTENUATION - it_of_this_subset + 1
    endif
  else
    it_tmp = it
  endif

  !debug
  !if (myrank == 0 ) print *,'compute_add_sources_backward: it_tmp = ',it_tmp,it

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
  !       this leads to the timing (NSTEP-(it-1)-1)*DT-t0-tshift_src for the source time function here
  !
  ! sets current initial time
  if (USE_LDDRK) then
    ! LDDRK
    ! note: the LDDRK scheme updates displacement after the stiffness computations and
    !       after adding boundary/coupling/source terms.
    !       thus, at each time loop step it, displ(:) is still at (n) and not (n+1) like for the Newmark scheme
    !       when entering this routine. we therefore at an additional -DT to have the corresponding timing for the source.
    if (UNDO_ATTENUATION) then
      ! stepping moves forward from snapshot position
      time_t = dble(NSTEP-it_tmp-1)*DT + dble(C_LDDRK(istage))*DT - t0
    else
      ! stepping backwards
      time_t = dble(NSTEP-it_tmp-1)*DT - dble(C_LDDRK(istage))*DT - t0
    endif
  else
    time_t = dble(NSTEP-it_tmp)*DT - t0
  endif

  if (.not. GPU_MODE) then
    ! on CPU
    do isource = 1,NSOURCES

      !   add the source (only if this proc carries the source)
      if (myrank == islice_selected_source(isource)) then

        ispec = ispec_selected_source(isource)

        ! sets current time for this source
        timeval = time_t - tshift_src(isource)

        ! determines source time function value
        stf = get_stf_viscoelastic(timeval,isource)

        ! distinguishes between single and double precision for reals
        stf_used = real(stf,kind=CUSTOM_REAL)

        ! adds source contribution
        do k = 1,NGLLZ
          do j = 1,NGLLY
            do i = 1,NGLLX
              iglob = ibool_crust_mantle(i,j,k,ispec)

              b_accel_crust_mantle(:,iglob) = b_accel_crust_mantle(:,iglob) + sourcearrays(:,i,j,k,isource) * stf_used

            enddo
          enddo
        enddo

      endif

    enddo

  else
    ! on GPU
    ! prepares buffer with source time function values, to be copied onto GPU
    do isource = 1,NSOURCES
      ! sets current time for this source
      timeval = time_t - tshift_src(isource)

      ! determines source time function value
      stf = get_stf_viscoelastic(timeval,isource)

      ! stores current stf values
      stf_pre_compute(isource) = stf
    enddo

    ! adds sources: only implements SIMTYPE=3 (and NOISE_TOM = 0)
    call compute_add_sources_backward_gpu(Mesh_pointer,NSOURCES,stf_pre_compute)
  endif

  end subroutine compute_add_sources_backward


!
!-------------------------------------------------------------------------------------------------
!

  double precision function get_stf_viscoelastic(time_source_dble,isource)

! returns source time function value for specified time

  use specfem_par, only: USE_FORCE_POINT_SOURCE,USE_MONOCHROMATIC_CMT_SOURCE,force_stf,hdur,hdur_Gaussian

  implicit none

  double precision,intent(in) :: time_source_dble
  integer,intent(in) :: isource

  ! local parameters
  double precision :: stf,f0

  double precision, external :: comp_source_time_function
  double precision, external :: comp_source_time_function_rickr
  double precision, external :: comp_source_time_function_gauss
  double precision, external :: comp_source_time_function_gauss_2
  double precision, external :: comp_source_time_function_mono

  ! note: calling comp_source_time_function() includes the handling for external source time functions

  ! determines source time function value
  if (USE_FORCE_POINT_SOURCE) then
    ! single point force
    select case(force_stf(isource))
    case (0)
      ! Gaussian source time function value
      stf = comp_source_time_function_gauss(time_source_dble,hdur_Gaussian(isource))
    case (1)
      ! Ricker source time function
      f0 = hdur(isource) ! using hdur as a FREQUENCY just to avoid changing FORCESOLUTION file format
      stf = comp_source_time_function_rickr(time_source_dble,f0)
    case (2)
      ! Heaviside (step) source time function
      stf = comp_source_time_function(time_source_dble,hdur_Gaussian(isource))
    case (3)
      ! Monochromatic source time function
      f0 = 1.d0 / hdur(isource) ! using hdur as a PERIOD just to avoid changing FORCESOLUTION file format
      stf = comp_source_time_function_mono(time_source_dble,f0)
    case (4)
      ! Gaussian source time function by Meschede et al. (2011)
      stf = comp_source_time_function_gauss_2(time_source_dble,hdur(isource))
    case default
      stop 'unsupported force_stf value!'
    end select
  else
    ! moment-tensor
    ! Heaviside source time function
    if (USE_MONOCHROMATIC_CMT_SOURCE) then
      f0 = 1.d0 / hdur(isource) ! using half duration as a FREQUENCY just to avoid changing CMTSOLUTION file format
      stf = comp_source_time_function_mono(time_source_dble,f0)
    else
      stf = comp_source_time_function(time_source_dble,hdur_Gaussian(isource))
    endif
  endif

  ! return value
  get_stf_viscoelastic = stf

  end function get_stf_viscoelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine get_stf_adjoint_source(it_in,irec_local,stf_array)

  use constants, only: CUSTOM_REAL,NDIM,C_LDDRK

  use specfem_par, only: NSTEP,NTSTEP_BETWEEN_READ_ADJSRC, &
    source_adjoint,iadj_vec, &
    USE_LDDRK,istage

  !debug
  !use specfem_par, only: myrank,NSTAGE_TIME_SCHEME

  implicit none

! note: defining it as a subroutine here instead of a function. Fortran90/95 could return arrays as function, but would need
!       an interface definition in the subroutine which calls it (or put the function into a module).
!       so far, specfem doesn't use this technique and there are no such interface definitions within subroutines.
!       so let's just use a subroutine call for this.

  integer,intent(in) :: it_in,irec_local
  real(kind=CUSTOM_REAL),dimension(NDIM),intent(out) :: stf_array

  ! local parameters
  integer :: ivec_index
  double precision,dimension(NDIM) :: stf
  ! cubic interpolation
  integer :: ivec_index0,ivec_index1,ivec_index2,ivec_index3,idx
  double precision,dimension(NDIM) :: p0,p1,p2,p3
  double precision,dimension(NDIM) :: a,b,c,d
  double precision :: t

  ! gets stf for adjoint source
  if (USE_LDDRK) then
    ! LDDRK
    ! needs interpolation for different stages. the adjoint trace has only points stored at it,
    ! and not for the additional NSTAGES of the LDDRK scheme. we will interpolate between it and it+1 points.
    !
    ! note: due to the update step to displ(n+1) at the end of the compute_forces** routine as compared to the Newmark scheme
    !       (which does it before the compute_forces**), we have to shift the adjoint source by -1.
    !       iadj_vec(it) goes from iadj_vec(1) = 1000, iadj_vec(2) = 999 to iadj_vec(1000) = 1
    if (istage == 1) then
      ! exact position at time it
      if (it_in == 1) then
        idx = it_in
      else
        idx = it_in - 1
      endif
      ivec_index = iadj_vec(idx)
      ! adjoint source time function for 3 components
      stf(:) = source_adjoint(:,irec_local,ivec_index)
    else
      ! position at time it + ct where the fraction ct = C_LDDRK(istage) between ]0,1[
      ! thus, we will need to interpolate between it and it + 1 points of the adjoint source
      !
      ! Catmull-Rom (cubic) interpolation:
      ! needs 4 points p0,p1,p2,p3 and interpolates point p(t) between point p1 and p2
      ! see: https://www.iquilezles.org/www/articles/minispline/minispline.htm
      t = C_LDDRK(istage)
      ! shift index by -1 due to LDDRK scheme
      if (it_in == 1) then
        idx = 1
      else
        idx = it_in - 1
      endif
      ivec_index = iadj_vec(idx)

      ! interpolation points
      ! note: at the boundaries it == 1, it < NSTEP-1, having "no" control points p0,p3 or end point p2 is probably still fine
      !       as we usually taper the adjoint source at the end, so the interpolation could shrink down to just taking
      !       the value at it. still, we take the limited values and mimick interpolation even in these cases.
      if (idx == 1) then
        ivec_index0 = iadj_vec(idx)
        ivec_index1 = iadj_vec(idx)
        ivec_index2 = iadj_vec(idx+1)
        ivec_index3 = iadj_vec(idx+2)
      else if (idx == NSTEP - 1) then
        ivec_index0 = iadj_vec(idx-1)
        ivec_index1 = iadj_vec(idx)
        ivec_index2 = iadj_vec(idx+1)
        ivec_index3 = iadj_vec(idx+1)
      else if (idx == NSTEP) then
        ivec_index0 = iadj_vec(idx-1)
        ivec_index1 = iadj_vec(idx)
        ivec_index2 = iadj_vec(idx)
        ivec_index3 = iadj_vec(idx)
      else
        ! it > 1 .and. it < NSTEP-1
        ivec_index0 = iadj_vec(idx-1)
        ivec_index1 = iadj_vec(idx)
        ivec_index2 = iadj_vec(idx+1)
        ivec_index3 = iadj_vec(idx+2)
      endif
      ! checks bounds
      if (ivec_index0 < 1) ivec_index0 = 1
      if (ivec_index1 < 1) ivec_index1 = 1
      if (ivec_index2 < 1) ivec_index2 = 1
      if (ivec_index3 < 1) ivec_index3 = 1
      if (ivec_index0 > NTSTEP_BETWEEN_READ_ADJSRC) ivec_index0 = NTSTEP_BETWEEN_READ_ADJSRC
      if (ivec_index1 > NTSTEP_BETWEEN_READ_ADJSRC) ivec_index1 = NTSTEP_BETWEEN_READ_ADJSRC
      if (ivec_index2 > NTSTEP_BETWEEN_READ_ADJSRC) ivec_index2 = NTSTEP_BETWEEN_READ_ADJSRC
      if (ivec_index3 > NTSTEP_BETWEEN_READ_ADJSRC) ivec_index3 = NTSTEP_BETWEEN_READ_ADJSRC
      ! interpolation points
      p0(:) = source_adjoint(:,irec_local,ivec_index0)
      p1(:) = source_adjoint(:,irec_local,ivec_index1)
      p2(:) = source_adjoint(:,irec_local,ivec_index2)
      p3(:) = source_adjoint(:,irec_local,ivec_index3)
      ! Catmull-Rom interpolation
      a(:) = 2.d0 * p1(:)
      b(:) = p2(:) - p0(:)
      c(:) = 2.d0 * p0(:) - 5.d0 * p1(:) + 4.d0 * p2(:) - p3(:)
      d(:) = -p0(:) + 3.d0 * p1(:) - 3.d0 * p2(:) + p3(:)
      ! cubic polynomial: a + b * t + c * t^2 + d * t^3
      stf(:) = 0.5d0 * ( a(:) + (b(:) * t) + (c(:) * t * t) + (d(:) * t * t * t) )
    endif
  else
    ! Newmark
    ! has 1 stage, at index iadj_vec
    ! adjoint source array index
    ivec_index = iadj_vec(it_in)
    ! adjoint source time function for 3 components
    stf(:) = source_adjoint(:,irec_local,ivec_index)
  endif

  ! return value
  stf_array(:) = real(stf(:),kind=CUSTOM_REAL)

  !debug
  !if (irec_local == 1) print *,myrank,(it_in-1)*NSTAGE_TIME_SCHEME + istage,stf_array(1),stf_array(2),stf_array(3),'#i #stf'

  end subroutine get_stf_adjoint_source
