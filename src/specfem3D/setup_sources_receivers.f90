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

  subroutine setup_sources_receivers(NSOURCES,myrank,ibool_crust_mantle, &
                      xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                      xigll,yigll,zigll,TOPOGRAPHY, &
                      sec,tshift_cmt,theta_source,phi_source, &
                      NSTEP,DT,hdur,hdur_gaussian,t0,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                      islice_selected_source,ispec_selected_source, &
                      xi_source,eta_source,gamma_source,nu_source, &
                      rspl,espl,espl2,nspl,ibathy_topo,NEX_XI,PRINT_SOURCE_TIME_FUNCTION, &
                      rec_filename,nrec,islice_selected_rec,ispec_selected_rec, &
                      xi_receiver,eta_receiver,gamma_receiver,station_name,network_name, &
                      stlat,stlon,stele,stbur,nu, &
                      nrec_local,nadj_rec_local,nrec_simulation, &
                      SIMULATION_TYPE,RECEIVERS_CAN_BE_BURIED,MOVIE_SURFACE,MOVIE_VOLUME, &
                      HDUR_MOVIE,OUTPUT_FILES,LOCAL_PATH)


  implicit none

  include 'mpif.h'
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer NSOURCES,myrank

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  logical TOPOGRAPHY

  double precision sec,DT,t0,min_tshift_cmt_original

  double precision, dimension(NSOURCES) :: tshift_cmt,hdur,hdur_gaussian
  double precision, dimension(NSOURCES) :: theta_source,phi_source
  double precision, dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source,nu_source

  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  integer NSTEP

  ! for ellipticity
  integer nspl
  double precision rspl(NR),espl(NR),espl2(NR)

  integer, dimension(NX_BATHY,NY_BATHY) :: ibathy_topo

  integer NEX_XI
  logical PRINT_SOURCE_TIME_FUNCTION

  character(len=150) rec_filename

  integer nrec
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec

  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  double precision, dimension(nrec) :: stlat,stlon,stele,stbur
  double precision, dimension(NDIM,NDIM,nrec) :: nu

  integer nrec_local,nadj_rec_local,nrec_simulation

  integer SIMULATION_TYPE

  logical RECEIVERS_CAN_BE_BURIED,MOVIE_SURFACE,MOVIE_VOLUME

  double precision HDUR_MOVIE

  character(len=150) OUTPUT_FILES
  character(len=150) LOCAL_PATH

  ! local parameters
  double precision :: junk
  integer :: yr,jda,ho,mi
  integer :: irec,isource,nrec_tot_found,ier
  integer :: icomp,itime,nadj_files_found,nadj_files_found_tot
  character(len=3),dimension(NDIM) :: comp
  character(len=256) :: filename,adj_source_file,system_command,filename_new
  character(len=2) :: bic
  ! makes smaller hdur for movies
  logical,parameter :: USE_SMALLER_HDUR_MOVIE = .false.

! sources
  ! BS BS moved open statement and writing of first lines into sr.vtk before the
  ! call to locate_sources, where further write statements to that file follow
  if(myrank == 0) then
  ! write source and receiver VTK files for Paraview
    filename = trim(OUTPUT_FILES)//'/sr_tmp.vtk'
    open(IOVTK,file=trim(filename),status='unknown')
    write(IOVTK,'(a)') '# vtk DataFile Version 2.0'
    write(IOVTK,'(a)') 'Source and Receiver VTK file'
    write(IOVTK,'(a)') 'ASCII'
    write(IOVTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
    !  LQY -- won't be able to know NSOURCES+nrec at this point...
    write(IOVTK, '(a,i6,a)') 'POINTS ', NSOURCES, ' float'
  endif

  ! locate sources in the mesh
  call locate_sources(NSOURCES,myrank,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,ibool_crust_mantle, &
            xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
            xigll,yigll,zigll,NPROCTOT_VAL,ELLIPTICITY_VAL,TOPOGRAPHY, &
            sec,tshift_cmt,min_tshift_cmt_original,yr,jda,ho,mi,theta_source,phi_source, &
            NSTEP,DT,hdur,Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
            islice_selected_source,ispec_selected_source, &
            xi_source,eta_source,gamma_source, nu_source, &
            rspl,espl,espl2,nspl,ibathy_topo,NEX_XI,PRINT_SOURCE_TIME_FUNCTION, &
            LOCAL_PATH,SIMULATION_TYPE)

  if(abs(minval(tshift_cmt)) > TINYVAL) call exit_MPI(myrank,'one tshift_cmt must be zero, others must be positive')

  ! filter source time function by Gaussian with hdur = HDUR_MOVIE when outputing movies or shakemaps
  if (MOVIE_SURFACE .or. MOVIE_VOLUME ) then
    ! smaller hdur_movie will do
    if( USE_SMALLER_HDUR_MOVIE ) then
      ! hdur_movie gets assigned an automatic value based on the simulation resolution
      ! this will make that a bit smaller to have a higher-frequency movie output
      HDUR_MOVIE = 0.5* HDUR_MOVIE
    endif

    ! new hdur for simulation
    hdur = sqrt(hdur**2 + HDUR_MOVIE**2)
    if(myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Each source is being convolved with HDUR_MOVIE = ',HDUR_MOVIE
      write(IMAIN,*)
    endif
  endif

  ! convert the half duration for triangle STF to the one for gaussian STF
  hdur_gaussian(:) = hdur(:)/SOURCE_DECAY_MIMIC_TRIANGLE

  ! define t0 as the earliest start time
  t0 = - 1.5d0*minval( tshift_cmt(:) - hdur(:) )

  ! point force sources will start depending on the frequency given by hdur
  if( USE_FORCE_POINT_SOURCE ) then
    ! note: point force sources will give the dominant frequency in hdur,
    !          thus the main period is 1/hdur.
    !          also, these sources use a Ricker source time function instead of a gaussian.
    !          for a Ricker source time function, a start time ~1.2 * main_period is a good choice
    t0 = - 1.2d0 * minval(tshift_cmt(:) - 1.0d0/hdur(:))
  endif

  ! checks if user set USER_T0 to fix simulation start time
  ! note: USER_T0 has to be positive
  if( USER_T0 > 0.d0 ) then
    ! user cares about origin time and time shifts of the CMTSOLUTION
    ! and wants to fix simulation start time to a constant start time
    ! time 0 on time axis will correspond to given origin time

    ! notifies user
    if( myrank == 0 ) then
      write(IMAIN,*) 'USER_T0: ',USER_T0
      write(IMAIN,*) 't0: ',t0,'min_tshift_cmt_original: ',min_tshift_cmt_original
      write(IMAIN,*)
    endif

    ! checks if automatically set t0 is too small
    ! note: min_tshift_cmt_original can be a positive or negative time shift (minimum from all tshift)
    if( t0 <= USER_T0 + min_tshift_cmt_original ) then
      ! by default, tshift_cmt(:) holds relative time shifts with a minimum time shift set to zero
      ! re-adds (minimum) original time shift such that sources will kick in
      ! according to their absolute time shift
      tshift_cmt(:) = tshift_cmt(:) + min_tshift_cmt_original

      ! sets new simulation start time such that
      ! simulation starts at t = - t0 = - USER_T0
      t0 = USER_T0

      ! notifies user
      if( myrank == 0 ) then
        write(IMAIN,*) '  set new simulation start time: ', - t0
        write(IMAIN,*)
      endif
    else
      ! start time needs to be at least t0 for numerical stability
      ! notifies user
      if( myrank == 0 ) then
        write(IMAIN,*) 'error: USER_T0 is too small'
        write(IMAIN,*) '       must make one of three adjustements:'
        write(IMAIN,*) '       - increase USER_T0 to be at least: ',t0-min_tshift_cmt_original
        write(IMAIN,*) '       - decrease time shift in CMTSOLUTION file'
        write(IMAIN,*) '       - decrease hdur in CMTSOLUTION file'
      endif
      call exit_mpi(myrank,'error USER_T0 is set but too small')
    endif
  else if( USER_T0 < 0.d0 ) then
    if( myrank == 0 ) then
      write(IMAIN,*) 'error: USER_T0 is negative, must be set zero or positive!'
    endif
    call exit_mpi(myrank,'error negative USER_T0 parameter in constants.h')
  endif

  !  receivers
  if(myrank == 0) then
    write(IMAIN,*)
    if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
      write(IMAIN,*) 'Total number of receivers = ', nrec
    else
      write(IMAIN,*) 'Total number of adjoint sources = ', nrec
    endif
    write(IMAIN,*)
  endif

  ! locate receivers in the crust in the mesh
  call locate_receivers(myrank,DT,NSTEP,NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,ibool_crust_mantle, &
                      xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                      xigll,yigll,zigll,trim(rec_filename), &
                      nrec,islice_selected_rec,ispec_selected_rec, &
                      xi_receiver,eta_receiver,gamma_receiver,station_name,network_name, &
                      stlat,stlon,stele,stbur,nu, &
                      yr,jda,ho,mi,sec,NPROCTOT_VAL,ELLIPTICITY_VAL,TOPOGRAPHY, &
                      theta_source(1),phi_source(1),rspl,espl,espl2,nspl, &
                      ibathy_topo,RECEIVERS_CAN_BE_BURIED,NCHUNKS_VAL)


  ! count number of receivers located in this slice
  nrec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    nrec_simulation = nrec
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec)) nrec_local = nrec_local + 1
    enddo
  else
    nrec_simulation = NSOURCES
    do isource = 1, NSOURCES
      if(myrank == islice_selected_source(isource)) nrec_local = nrec_local + 1
    enddo
  endif

  ! counts receivers for adjoint simulations
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    ! by Ebru
    call band_instrument_code(DT,bic)
    comp(1) = bic(1:2)//'N'
    comp(2) = bic(1:2)//'E'
    comp(3) = bic(1:2)//'Z'

    ! counter for adjoint receiver stations in local slice, used to allocate adjoint source arrays
    nadj_rec_local = 0
    ! temporary counter to check if any files are found at all
    nadj_files_found = 0
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec))then
        ! adjoint receiver station in this process slice
        if(islice_selected_rec(irec) < 0 .or. islice_selected_rec(irec) > NPROCTOT_VAL-1) &
          call exit_MPI(myrank,'something is wrong with the source slice number in adjoint simulation')

        ! updates counter
        nadj_rec_local = nadj_rec_local + 1

        ! checks **sta**.**net**.**MX**.adj files for correct number of time steps
        adj_source_file = trim(station_name(irec))//'.'//trim(network_name(irec))
        do icomp = 1,NDIM

          ! opens adjoint source file for this component
          filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
          open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ier)

          if( ier == 0 ) then
            ! checks length of file
            itime = 0
            do while(ier == 0)
              read(IIN,*,iostat=ier) junk,junk
              if( ier == 0 ) itime = itime + 1
            enddo
            if( itime /= NSTEP) &
              call exit_MPI(myrank,&
                'file '//trim(filename)//' has wrong length, please check with your simulation duration')

            ! updates counter for found files
            nadj_files_found = nadj_files_found + 1
          else
            ! adjoint source file not found
            ! stops simulation
            call exit_MPI(myrank,&
                'file '//trim(filename)//' not found, please check with your STATIONS_ADJOINT file')
          endif
          close(IIN)
        enddo
      endif
    enddo

    ! checks if any adjoint source files found at all
    call MPI_REDUCE(nadj_files_found,nadj_files_found_tot,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)
    if( myrank == 0 ) then
      write(IMAIN,*)
      write(IMAIN,*) '    ',nadj_files_found_tot,' adjoint component traces found in all slices'
      if(nadj_files_found_tot == 0) &
        call exit_MPI(myrank,'no adjoint traces found, please check adjoint sources in directory SEM/')
    endif
  endif

  ! check that the sum of the number of receivers in each slice is nrec
  call MPI_REDUCE(nrec_local,nrec_tot_found,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ier)
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'found a total of ',nrec_tot_found,' receivers in all slices'
    if(nrec_tot_found /= nrec_simulation) then
      call exit_MPI(myrank,'problem when dispatching the receivers')
    else
      write(IMAIN,*) 'this total is okay'
    endif
  endif

  ! user output
  if(myrank == 0) then

    ! finishes vtk file
    write(IOVTK,*) ""
    close(IOVTK)

    !  we should know NSOURCES+nrec at this point...
    ! creates source/receiver location file
    filename = trim(OUTPUT_FILES)//'/sr_tmp.vtk'
    filename_new = trim(OUTPUT_FILES)//'/sr.vtk'
    write(system_command, &
  "('sed -e ',a1,'s/POINTS.*/POINTS',i6,' float/',a1,' < ',a,' > ',a)")&
      "'",NSOURCES + nrec,"'",trim(filename),trim(filename_new)
    call system(system_command)

    ! only extract receiver locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/receiver.vtk'
    write(system_command, &
  "('awk ',a1,'{if(NR<5) print $0;if(NR==6)print ',a1,'POINTS',i6,' float',a1,';if(NR>5+',i6,')print $0}',a1,' < ',a,' > ',a)")&
      "'",'"',nrec,'"',NSOURCES,"'",trim(filename),trim(filename_new)
    call system(system_command)

    ! only extract source locations and remove temporary file
    filename_new = trim(OUTPUT_FILES)//'/source.vtk'
    write(system_command, &
  "('awk ',a1,'{if(NR< 6 + ',i6,') print $0}END{print}',a1,' < ',a,' > ',a,'; rm -f ',a)")&
      "'",NSOURCES,"'",trim(filename),trim(filename_new),trim(filename)
    call system(system_command)


    write(IMAIN,*)
    write(IMAIN,*) 'Total number of samples for seismograms = ',NSTEP
    write(IMAIN,*)


    if(NSOURCES > 1) write(IMAIN,*) 'Using ',NSOURCES,' point sources'
  endif

  end subroutine setup_sources_receivers

!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_sources_receivers_srcarr(NSOURCES,myrank, &
                      ispec_selected_source,islice_selected_source, &
                      xi_source,eta_source,gamma_source, &
                      Mxx,Myy,Mzz,Mxy,Mxz,Myz, &
                      xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
                      etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
                      gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
                      xigll,yigll,zigll,sourcearrays)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer NSOURCES,myrank

  integer, dimension(NSOURCES) :: islice_selected_source,ispec_selected_source
  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(NSOURCES) :: Mxx,Myy,Mzz,Mxy,Mxz,Myz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: &
        xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle,&
        etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
        gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle

  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSOURCES) :: sourcearrays


  ! local parameters
  integer :: isource
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ) :: sourcearray

  do isource = 1,NSOURCES

    !   check that the source slice number is okay
    if(islice_selected_source(isource) < 0 .or. islice_selected_source(isource) > NPROCTOT_VAL-1) &
      call exit_MPI(myrank,'something is wrong with the source slice number')

    !   compute source arrays in source slice
    if(myrank == islice_selected_source(isource)) then
      call compute_arrays_source(ispec_selected_source(isource), &
             xi_source(isource),eta_source(isource),gamma_source(isource),sourcearray, &
             Mxx(isource),Myy(isource),Mzz(isource),Mxy(isource),Mxz(isource),Myz(isource), &
             xix_crust_mantle,xiy_crust_mantle,xiz_crust_mantle, &
             etax_crust_mantle,etay_crust_mantle,etaz_crust_mantle, &
             gammax_crust_mantle,gammay_crust_mantle,gammaz_crust_mantle, &
             xigll,yigll,zigll,NSPEC_CRUST_MANTLE)

      sourcearrays(:,:,:,:,isource) = sourcearray(:,:,:,:)

    endif
  enddo

  end subroutine setup_sources_receivers_srcarr


!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_sources_receivers_adjindx(NSTEP,NSTEP_SUB_ADJ, &
                      NTSTEP_BETWEEN_READ_ADJSRC, &
                      iadjsrc,iadjsrc_len,iadj_vec)

  implicit none

  include "constants.h"

  integer NSTEP,NSTEP_SUB_ADJ,NTSTEP_BETWEEN_READ_ADJSRC

  integer, dimension(NSTEP_SUB_ADJ,2) :: iadjsrc ! to read input in chunks
  integer, dimension(NSTEP_SUB_ADJ) :: iadjsrc_len
  integer, dimension(NSTEP) :: iadj_vec


  ! local parameters
  integer :: iadj_block,it,it_sub_adj

  iadj_block = 1  !counts blocks

  iadjsrc(:,:) = 0
  iadjsrc_len(:) = 0

  ! setting up chunks of NTSTEP_BETWEEN_READ_ADJSRC to read adjoint source traces
  ! i.e. as an example: total length NSTEP = 3000, chunk length NTSTEP_BETWEEN_READ_ADJSRC= 1000
  !                                then it will set first block from 2001 to 3000,
  !                                second block from 1001 to 2000 and so on...
  !
  ! see routine: compute_arrays_source_adjoint()
  !                     how we read in the adjoint source trace in blocks/chunk sizes
  !
  ! see routine: compute_add_sources_adjoint()
  !                     how the adjoint source is added to the (adjoint) acceleration field
  do it=1,NSTEP

    ! block number
    ! e.g. increases from 1 (case it=1-1000), 2 (case it=1001-2000) to 3 (case it=2001-3000)
    it_sub_adj = ceiling( dble(it)/dble(NTSTEP_BETWEEN_READ_ADJSRC) )

    ! we are at the edge of a block
    if(mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC) == 0) then
     ! block start time ( e.g. 2001)
     iadjsrc(iadj_block,1) = NSTEP-it_sub_adj*NTSTEP_BETWEEN_READ_ADJSRC+1
     ! block end time (e.g. 3000)
     iadjsrc(iadj_block,2) = NSTEP-(it_sub_adj-1)*NTSTEP_BETWEEN_READ_ADJSRC

     ! final adj src array
     ! e.g. will be from 1000 to 1, but doesn't go below 1 in cases where NSTEP isn't
     ! a multiple of NTSTEP_BETWEEN_READ_ADJSRC
     if(iadjsrc(iadj_block,1) < 0) iadjsrc(iadj_block,1) = 1

     ! actual block length
     iadjsrc_len(iadj_block) = iadjsrc(iadj_block,2)-iadjsrc(iadj_block,1)+1

     ! increases block number
     iadj_block = iadj_block+1
    endif

    ! time stepping for adjoint sources:
    ! adjoint time step that corresponds to time step in simulation (it).
    ! note, that adjoint source has to be time-reversed with respect to the forward wavefield
    ! e.g.: first block 1 has iadjsrc_len = 1000 with start at 2001 and end at 3000
    !         so iadj_vec(1) = 1000 - 0, iadj_vec(2) = 1000 - 1, ..., to iadj_vec(1000) = 1000 - 999 = 1
    !         then for block 2, iadjsrc_len = 1000 with start at 1001 and end at 2000
    !         so iadj_vec(1001) = 1000 - 0, iad_vec(1002) = 1000 - 1, .. and so on again down to 1
    !         then block 3 and your guess is right now... iadj_vec(2001) to iadj_vec(3000) is 1000 down to 1. :)
    iadj_vec(it) = iadjsrc_len(it_sub_adj) - mod(it-1,NTSTEP_BETWEEN_READ_ADJSRC)
  enddo

  end subroutine setup_sources_receivers_adjindx

!
!-------------------------------------------------------------------------------------------------
!


  subroutine setup_sources_receivers_intp(NSOURCES,myrank, &
                      islice_selected_source, &
                      xi_source,eta_source,gamma_source, &
                      xigll,yigll,zigll, &
                      SIMULATION_TYPE,nrec,nrec_local, &
                      islice_selected_rec,number_receiver_global, &
                      xi_receiver,eta_receiver,gamma_receiver, &
                      hxir_store,hetar_store,hgammar_store, &
                      nadj_hprec_local,hpxir_store,hpetar_store,hpgammar_store)

  implicit none

  include "constants.h"

  integer NSOURCES,myrank

  integer, dimension(NSOURCES) :: islice_selected_source

  double precision, dimension(NSOURCES) :: xi_source,eta_source,gamma_source
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll


  integer SIMULATION_TYPE

  integer nrec,nrec_local
  integer, dimension(nrec) :: islice_selected_rec
  integer, dimension(nrec_local) :: number_receiver_global
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver

  double precision, dimension(nrec_local,NGLLX) :: hxir_store
  double precision, dimension(nrec_local,NGLLY) :: hetar_store
  double precision, dimension(nrec_local,NGLLZ) :: hgammar_store

  integer nadj_hprec_local
  double precision, dimension(nadj_hprec_local,NGLLX) :: hpxir_store
  double precision, dimension(nadj_hprec_local,NGLLY) :: hpetar_store
  double precision, dimension(nadj_hprec_local,NGLLZ) :: hpgammar_store


  ! local parameters
  integer :: isource,irec,irec_local
  double precision, dimension(NGLLX) :: hxir,hpxir
  double precision, dimension(NGLLY) :: hpetar,hetar
  double precision, dimension(NGLLZ) :: hgammar,hpgammar


  ! select local receivers

  ! define local to global receiver numbering mapping
  irec_local = 0
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do irec = 1,nrec
      if(myrank == islice_selected_rec(irec)) then
        irec_local = irec_local + 1
        number_receiver_global(irec_local) = irec
      endif
    enddo
  else
    do isource = 1,NSOURCES
      if(myrank == islice_selected_source(isource)) then
        irec_local = irec_local + 1
        number_receiver_global(irec_local) = isource
      endif
    enddo
  endif

  ! define and store Lagrange interpolators at all the receivers
  if (SIMULATION_TYPE == 1 .or. SIMULATION_TYPE == 3) then
    do irec_local = 1,nrec_local
      irec = number_receiver_global(irec_local)
      call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_receiver(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
      hxir_store(irec_local,:) = hxir(:)
      hetar_store(irec_local,:) = hetar(:)
      hgammar_store(irec_local,:) = hgammar(:)
    enddo
  else
    do irec_local = 1,nrec_local
      irec = number_receiver_global(irec_local)
      call lagrange_any(xi_source(irec),NGLLX,xigll,hxir,hpxir)
      call lagrange_any(eta_source(irec),NGLLY,yigll,hetar,hpetar)
      call lagrange_any(gamma_source(irec),NGLLZ,zigll,hgammar,hpgammar)
      hxir_store(irec_local,:) = hxir(:)
      hetar_store(irec_local,:) = hetar(:)
      hgammar_store(irec_local,:) = hgammar(:)
      hpxir_store(irec_local,:) = hpxir(:)
      hpetar_store(irec_local,:) = hpetar(:)
      hpgammar_store(irec_local,:) = hpgammar(:)
    enddo
  endif

  end subroutine setup_sources_receivers_intp

