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

! subroutine for NOISE TOMOGRAPHY
! chracterize noise statistics
! for a given point (xcoord,ycoord,zcoord), specify the noise direction "normal_x/y/z_noise"
!     and noise distribution "mask_noise"
! USERS need to modify this subroutine for their own noise characteristics

  subroutine noise_distribution_direction(xcoord_in,ycoord_in,zcoord_in, &
                  normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                  mask_noise_out)
  implicit none
  include "constants.h"
  ! input parameters
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  ! output parameters
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  ! local parameters
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord


  ! coordinates "x/y/zcoord_in" actually contain r theta phi, therefore convert back to x y z
  call rthetaphi_2_xyz(xcoord,ycoord,zcoord, xcoord_in,ycoord_in,zcoord_in)
  ! NOTE that all coordinates are non-dimensionalized in GLOBAL package!
  ! USERS are free to choose which set to use,
  ! either "r theta phi" (xcoord_in,ycoord_in,zcoord_in)
  ! or     "x y z"       (xcoord,ycoord,zcoord)

  !*****************************************************************************************************************
  !******************************** change your noise characteristics below ****************************************
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! noise direction
  ! here, the noise is assumed to be vertical
  normal_x_noise_out = xcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  normal_y_noise_out = ycoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  normal_z_noise_out = zcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  noise distribution
  ! here, the noise is assumed to be uniform
  mask_noise_out = 1.0
  !******************************** change your noise characteristics above ****************************************
  !*****************************************************************************************************************

  end subroutine noise_distribution_direction

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! read parameters

  subroutine read_parameters_noise(myrank,nrec,NSTEP,nmovie_points, &
                                   islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu, &
                                   noise_sourcearray,xigll,yigll,zigll,nspec_top, &
                                   NIT, ibool_crust_mantle, ibelm_top_crust_mantle, &
                                   xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                                   irec_master_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise)
  implicit none

  include 'mpif.h'
  include "precision.h"
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: myrank, nrec, NSTEP, nmovie_points, nspec_top, NIT
  integer, dimension(nrec) :: islice_selected_rec
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  double precision, dimension(nrec)  :: xi_receiver,eta_receiver,gamma_receiver
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  double precision, dimension(NDIM,NDIM,nrec) :: nu
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  ! output parameters
  integer :: irec_master_noise
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise,mask_noise
  ! local parameters
  integer :: ipoin, ispec2D, ispec, i, j, k, iglob, ios, ier
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  character(len=150) :: filename
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: &
      store_val_x,store_val_y,store_val_z,  store_val_ux,store_val_uy,store_val_uz
  real(kind=CUSTOM_REAL), dimension(nmovie_points,0:NPROCTOT_VAL-1) :: &
      store_val_x_all,store_val_y_all,store_val_z_all, store_val_ux_all,store_val_uy_all,store_val_uz_all


  ! read master receiver ID -- the ID in DATA/STATIONS
  filename = 'OUTPUT_FILES/NOISE_TOMOGRAPHY/'//'irec_master_noise'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0)  &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file contains the ID of the master receiver')
  read(IIN_NOISE,*,iostat=ios) irec_master_noise
  close(IIN_NOISE)

  if (myrank == 0) then
    open(unit=IOUT_NOISE,file='OUTPUT_FILES/irec_master_noise', &
          status='unknown',action='write',iostat=ios)
    if( ios /= 0 ) call exit_MPI(myrank,'error opening output file irec_master_noise')

    WRITE(IOUT_NOISE,*) 'The master receiver is: (RECEIVER ID)', irec_master_noise
    close(IOUT_NOISE)
  endif

  ! checks master irec
  if( irec_master_noise < 1 .or. irec_master_noise > nrec ) then
    call exit_MPI(myrank,'error noise tomography: irec_master_noise is not in range of given number of receivers')
  endif

  ! compute source arrays for "ensemble forward source", which is source of "ensemble forward wavefield"
  if(myrank == islice_selected_rec(irec_master_noise) .OR. myrank == 0) then ! myrank == 0 is used for output only
    call compute_arrays_source_noise(myrank, &
              xi_receiver(irec_master_noise),eta_receiver(irec_master_noise),gamma_receiver(irec_master_noise), &
              nu(:,:,irec_master_noise),noise_sourcearray, xigll,yigll,zigll,NSTEP)
  endif

  ! noise distribution and noise direction
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ! this subroutine must be modified by USERS
        call noise_distribution_direction(xstore_crust_mantle(iglob), &
                  ystore_crust_mantle(iglob),zstore_crust_mantle(iglob), &
                  normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                  mask_noise_out)
        normal_x_noise(ipoin) = normal_x_noise_out
        normal_y_noise(ipoin) = normal_y_noise_out
        normal_z_noise(ipoin) = normal_z_noise_out
        mask_noise(ipoin)     = mask_noise_out
      enddo
    enddo

  enddo

  !!!BEGIN!!! save mask_noise for check, a file called "mask_noise" is saved in "./OUTPUT_FIELS/"
  ipoin = 0
  do ispec2D = 1, NSPEC_TOP ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    k = NGLLZ
    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)
        store_val_x(ipoin) = xstore_crust_mantle(iglob)
        store_val_y(ipoin) = ystore_crust_mantle(iglob)
        store_val_z(ipoin) = zstore_crust_mantle(iglob)
        store_val_ux(ipoin) = mask_noise(ipoin)
        store_val_uy(ipoin) = mask_noise(ipoin)
        store_val_uz(ipoin) = mask_noise(ipoin)
      enddo
    enddo
  enddo

  ! gather info on master proc
  ispec = nmovie_points
  call MPI_GATHER(store_val_x,ispec,CUSTOM_MPI_TYPE,store_val_x_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_y,ispec,CUSTOM_MPI_TYPE,store_val_y_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_z,ispec,CUSTOM_MPI_TYPE,store_val_z_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_ux,ispec,CUSTOM_MPI_TYPE,store_val_ux_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_uy,ispec,CUSTOM_MPI_TYPE,store_val_uy_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)
  call MPI_GATHER(store_val_uz,ispec,CUSTOM_MPI_TYPE,store_val_uz_all,ispec,CUSTOM_MPI_TYPE,0,MPI_COMM_WORLD,ier)

  ! save maks_noise data to disk in home directory
  ! this file can be viewed the same way as surface movie data (xcreate_movie_AVS_DX)
  ! create_movie_AVS_DX.f90 needs to be modified in order to do that,
  ! i.e., instead of showing the normal component, change it to either x, y or z component, or the norm.
  if(myrank == 0) then
    open(unit=IOUT_NOISE,file='OUTPUT_FILES/mask_noise', &
              status='unknown',form='unformatted',action='write',iostat=ios)
    if( ios /= 0 ) call exit_MPI(myrank,'error opening output file mask_noise')

    write(IOUT_NOISE) store_val_x_all
    write(IOUT_NOISE) store_val_y_all
    write(IOUT_NOISE) store_val_z_all
    write(IOUT_NOISE) store_val_ux_all
    write(IOUT_NOISE) store_val_uy_all
    write(IOUT_NOISE) store_val_uz_all

    close(IOUT_NOISE)
  endif
  !!!END!!! save mask_noise for check, a file called "mask_noise" is saved in "./OUTPUT_FIELS/"

  end subroutine read_parameters_noise

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! check for consistency of the parameters
  subroutine check_parameters_noise(myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                    NUMBER_OF_RUNS, NUMBER_OF_THIS_RUN,ROTATE_SEISMOGRAMS_RT, &
                                    SAVE_ALL_SEISMOS_IN_ONE_FILE, USE_BINARY_FOR_LARGE_FILE, &
                                    MOVIE_COARSE,LOCAL_PATH,NSPEC_TOP,NSTEP)
  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN,NSPEC_TOP,NSTEP
  logical :: SAVE_FORWARD,ROTATE_SEISMOGRAMS_RT,SAVE_ALL_SEISMOS_IN_ONE_FILE, &
            USE_BINARY_FOR_LARGE_FILE,MOVIE_COARSE
  character(len=150) :: LOCAL_PATH
  ! output parameters
  ! local parameters
  integer :: reclen,ier
  integer(kind=8) :: filesize
  character(len=150) :: outputname


  if (myrank == 0) then
     open(unit=IOUT_NOISE,file='OUTPUT_FILES/NOISE_SIMULATION', &
          status='unknown',action='write',iostat=ier)
     if( ier /= 0 ) call exit_MPI(myrank,'error opening output file NOISE_SIMULATION')

     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     WRITE(IOUT_NOISE,*) 'WARNING!!!!!!!!!!!!'
     WRITE(IOUT_NOISE,*) 'You are running simulations using NOISE TOMOGRAPHY techniques.'
     WRITE(IOUT_NOISE,*) 'Please make sure you understand the procedures before you have a try.'
     WRITE(IOUT_NOISE,*) 'Displacements everywhere at the free surface are saved every timestep,'
     WRITE(IOUT_NOISE,*) 'so make sure that LOCAL_PATH in DATA/Par_file is not global.'
     WRITE(IOUT_NOISE,*) 'Otherwise the disk storage may be a serious issue, as is the speed of I/O.'
     WRITE(IOUT_NOISE,*) 'Also make sure that NO earthquakes are included,'
     WRITE(IOUT_NOISE,*) 'i.e., set moment tensor to be ZERO in CMTSOLUTION'
     WRITE(IOUT_NOISE,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(IOUT_NOISE,*) 'If you just want a regular EARTHQUAKE simulation,'
     WRITE(IOUT_NOISE,*) 'set NOISE_TOMOGRAPHY=0 in DATA/Par_file'
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     WRITE(IOUT_NOISE,*) '*******************************************************************************'
     close(IOUT_NOISE)
  endif

  if (NUMBER_OF_RUNS/=1 .OR. NUMBER_OF_THIS_RUN/=1) &
     call exit_mpi(myrank,'NUMBER_OF_RUNS and NUMBER_OF_THIS_RUN must be 1 for NOISE TOMOGRAPHY! check DATA/Par_file')
  if (ROTATE_SEISMOGRAMS_RT) &
     call exit_mpi(myrank,'Do NOT rotate seismograms in the code, change ROTATE_SEISMOGRAMS_RT in DATA/Par_file')
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE .OR. USE_BINARY_FOR_LARGE_FILE) &
     call exit_mpi(myrank,'Please set SAVE_ALL_SEISMOS_IN_ONE_FILE and USE_BINARY_FOR_LARGE_FILE to be .false.')
  if (MOVIE_COARSE) &
     call exit_mpi(myrank,'Please set MOVIE_COARSE in DATA/Par_file to be .false.')


  if (NOISE_TOMOGRAPHY==1) then
     if (SIMULATION_TYPE/=1) &
        call exit_mpi(myrank,'NOISE_TOMOGRAPHY=1 requires SIMULATION_TYPE=1! check DATA/Par_file')
  else if (NOISE_TOMOGRAPHY==2) then
     if (SIMULATION_TYPE/=1) &
        call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SIMULATION_TYPE=1! check DATA/Par_file')
     if (.not. SAVE_FORWARD) &
        call exit_mpi(myrank,'NOISE_TOMOGRAPHY=2 requires SAVE_FORWARD=.true.! check DATA/Par_file')
  else if (NOISE_TOMOGRAPHY==3) then
     if (SIMULATION_TYPE/=3) &
        call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SIMULATION_TYPE=3! check DATA/Par_file')
     if (SAVE_FORWARD) &
        call exit_mpi(myrank,'NOISE_TOMOGRAPHY=3 requires SAVE_FORWARD=.false.! check DATA/Par_file')
  endif

  if (NOISE_TOMOGRAPHY/=0) then
     ! save/read the surface movie using the same c routine as we do for absorbing boundaries (file ID is 9)
     ! size of single record
     reclen=CUSTOM_REAL*NDIM*NGLLX*NGLLY*NSPEC_TOP
     ! total file size
     filesize = reclen
     filesize = filesize*NSTEP

     write(outputname,"('/proc',i6.6,'_surface_movie')") myrank
     if (NOISE_TOMOGRAPHY==1) call open_file_abs_w(9,trim(LOCAL_PATH)//trim(outputname), &
                                      len_trim(trim(LOCAL_PATH)//trim(outputname)), &
                                      filesize)
     if (NOISE_TOMOGRAPHY==2) call open_file_abs_r(9,trim(LOCAL_PATH)//trim(outputname), &
                                      len_trim(trim(LOCAL_PATH)//trim(outputname)), &
                                      filesize)
     if (NOISE_TOMOGRAPHY==3) call open_file_abs_r(9,trim(LOCAL_PATH)//trim(outputname), &
                                      len_trim(trim(LOCAL_PATH)//trim(outputname)), &
                                      filesize)
  endif

  end subroutine check_parameters_noise

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! read and construct the "source" (source time function based upon noise spectrum)
! for "ensemble forward source"

  subroutine compute_arrays_source_noise(myrank, &
                                         xi_noise,eta_noise,gamma_noise,nu_single,noise_sourcearray, &
                                         xigll,yigll,zigll,NSTEP)
  implicit none
  include 'constants.h'
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: myrank, NSTEP
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  double precision, dimension(NDIM,NDIM) :: nu_single  ! rotation matrix at the master receiver
  ! output parameters
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  ! local parameters
  integer itime, i, j, k, ios
  real(kind=CUSTOM_REAL) :: junk
  real(kind=CUSTOM_REAL) :: noise_src(NSTEP),noise_src_u(NDIM,NSTEP)
  double precision, dimension(NDIM) :: nu_master       ! component direction chosen at the master receiver
  double precision :: xi_noise, eta_noise, gamma_noise ! master receiver location
  double precision,parameter :: scale_displ_inv = 1.d0/R_EARTH ! non-dimesional scaling
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
  character(len=150) :: filename


  noise_src(:) = 0._CUSTOM_REAL
  ! noise file (source time function)
  filename = 'OUTPUT_FILES/NOISE_TOMOGRAPHY/'//'S_squared'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0)  &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file is generated by Matlab scripts')

  do itime =1,NSTEP
    read(IIN_NOISE,*,iostat=ios) junk, noise_src(itime)
    if( ios /= 0)  call exit_MPI(myrank,&
        'file '//trim(filename)//' has wrong length, please check your simulation duration')
  enddo
  close(IIN_NOISE)

  ! master receiver component direction, \nu_master
  filename = 'OUTPUT_FILES/NOISE_TOMOGRAPHY/'//'nu_master'
  open(unit=IIN_NOISE,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0)  call exit_MPI(myrank,&
        'file '//trim(filename)//' does NOT exist! nu_master is the component direction (NEZ) for master receiver')

  do itime =1,3
    read(IIN_NOISE,*,iostat=ios) nu_master(itime)
    if( ios /= 0) call exit_MPI(myrank,&
        'file '//trim(filename)//' has wrong length, the vector should have three components (NEZ)')
  enddo
  close(IIN_NOISE)

  if (myrank == 0) then
     open(unit=IOUT_NOISE,file='OUTPUT_FILES/nu_master',status='unknown',action='write')
     WRITE(IOUT_NOISE,*) 'The direction (NEZ) of selected component of master receiver is', nu_master
     close(IOUT_NOISE)
  endif

  ! rotates to cartesian
  do itime = 1, NSTEP
    noise_src_u(:,itime) = nu_single(1,:) * noise_src(itime) * nu_master(1) &
                         + nu_single(2,:) * noise_src(itime) * nu_master(2) &
                         + nu_single(3,:) * noise_src(itime) * nu_master(3)
  enddo

  ! receiver interpolators
  call lagrange_any(xi_noise,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta_noise,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma_noise,NGLLZ,zigll,hgammar,hpgammar)

  ! adds interpolated source contribution to all GLL points within this element
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        do itime = 1, NSTEP
          noise_sourcearray(:,i,j,k,itime) = hxir(i) * hetar(j) * hgammar(k) * noise_src_u(:,itime)
        enddo
      enddo
    enddo
  enddo

  end subroutine compute_arrays_source_noise

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! step 1: calculate the "ensemble forward source"
! add noise spectrum to the location of master receiver
  subroutine add_source_master_rec_noise(myrank,nrec, &
                                NSTEP,accel_crust_mantle,noise_sourcearray, &
                                ibool_crust_mantle,islice_selected_rec,ispec_selected_rec, &
                                it,irec_master_noise)
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: myrank,nrec,NSTEP, irec_master_noise
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP) :: noise_sourcearray
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle  ! both input and output
  ! output parameters

  ! local parameters
  integer :: i,j,k,iglob,it


  ! adds noise source (only if this proc carries the noise)
  if(myrank == islice_selected_rec(irec_master_noise)) then
    ! adds noise source contributions
    do k=1,NGLLZ
      do j=1,NGLLY
        do i=1,NGLLX
          iglob = ibool_crust_mantle(i,j,k,ispec_selected_rec(irec_master_noise))
          accel_crust_mantle(:,iglob) = accel_crust_mantle(:,iglob) &
                        + noise_sourcearray(:,i,j,k,it)
        enddo
      enddo
    enddo
  endif

  end subroutine add_source_master_rec_noise

!
!-------------------------------------------------------------------------------------------------
!

! subroutine for NOISE TOMOGRAPHY
! step 1: calculate the "ensemble forward source"
! save surface movie (displacement) at every time steps, for step 2 & 3.
!
! there are two subroutines --- noise_save_surface_movie_original & noise_save_surface_movie
!    noise_save_surface_movie_original is implemented at first, which creates one file at each time step
!    noise_save_surface_movie is implemented later, which utilizes 'src/shared/write_c_binary.c' for faster I/O,
!                                                   which creates one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage

!!!!! improved version !!!!!
  subroutine noise_save_surface_movie(displ_crust_mantle, &
                    ibelm_top_crust_mantle,ibool_crust_mantle, &
                    nspec_top,noise_surface_movie,it)
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: nspec_top,it
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) ::  displ_crust_mantle
  ! output parameters

  ! local parameters
  integer :: ispec2D,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL),dimension(NDIM,NGLLX,NGLLY,nspec_top) :: noise_surface_movie

  ! get coordinates of surface mesh and surface displacement
  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)
    k = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool_crust_mantle(i,j,k,ispec)
        noise_surface_movie(:,i,j,ispec2D) = displ_crust_mantle(:,iglob)
      enddo
    enddo
  enddo

  ! save surface motion to disk
  call write_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*nspec_top,it)

  end subroutine noise_save_surface_movie

!!!!! original implementation, not used anymore (but kept here for references) !!!!!
!  subroutine noise_save_surface_movie_original(myrank,nmovie_points,displ_crust_mantle, &
!                    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
!                    store_val_x,store_val_y,store_val_z, &
!                    store_val_ux,store_val_uy,store_val_uz, &
!                    ibelm_top_crust_mantle,ibool_crust_mantle,nspec_top, &
!                    NIT,it,LOCAL_PATH)
!  implicit none
!  include 'mpif.h'
!  include "precision.h"
!  include "constants.h"
!  include "OUTPUT_FILES/values_from_mesher.h"
!  ! input parameters
!  integer :: myrank,nmovie_points,nspec_top,NIT,it
!  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) ::  displ_crust_mantle
!  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
!        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
!  character(len=150) :: LOCAL_PATH
!  ! output parameters
!  ! local parameters
!  integer :: ipoin,ispec2D,ispec,i,j,k,iglob
!  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: &
!      store_val_x,store_val_y,store_val_z, &
!      store_val_ux,store_val_uy,store_val_uz
!  character(len=150) :: outputname
!
!
!  ! get coordinates of surface mesh and surface displacement
!  ipoin = 0
!  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
!    ispec = ibelm_top_crust_mantle(ispec2D)
!
!    k = NGLLZ
!
!    ! loop on all the points inside the element
!    do j = 1,NGLLY,NIT
!      do i = 1,NGLLX,NIT
!        ipoin = ipoin + 1
!        iglob = ibool_crust_mantle(i,j,k,ispec)
!        store_val_x(ipoin) = xstore_crust_mantle(iglob)
!        store_val_y(ipoin) = ystore_crust_mantle(iglob)
!        store_val_z(ipoin) = zstore_crust_mantle(iglob)
!        store_val_ux(ipoin) = displ_crust_mantle(1,iglob)
!        store_val_uy(ipoin) = displ_crust_mantle(2,iglob)
!        store_val_uz(ipoin) = displ_crust_mantle(3,iglob)
!      enddo
!    enddo
!
!  enddo
!
!  ! save surface motion to disk
!  ! LOCAL storage is better than GLOBAL, because we have to save the 'movie' at every time step
!  ! also note that the surface movie does NOT have to be shared with other nodes/CPUs
!  ! change LOCAL_PATH specified in "DATA/Par_file"
!    write(outputname,"('/proc',i6.6,'_surface_movie',i6.6)") myrank, it
!    open(unit=IOUT_NOISE,file=trim(LOCAL_PATH)//outputname,status='unknown',form='unformatted',action='write')
!    write(IOUT_NOISE) store_val_ux
!    write(IOUT_NOISE) store_val_uy
!    write(IOUT_NOISE) store_val_uz
!    close(IOUT_NOISE)
!  end subroutine noise_save_surface_movie_original

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 2/3: calculate/reconstructe the "ensemble forward wavefield"
! read surface movie (displacement) at every time steps, injected as the source of "ensemble forward wavefield"
! in step 2, call noise_read_add_surface_movie(..., NSTEP-it+1 ,...)
! in step 3, call noise_read_add_surface_movie(..., it ,...)
!
! there are two subroutines --- noise_read_add_surface_movie_original & noise_read_add_surface_movie
!    noise_read_add_surface_movie_original is implemented at first, which creates one file at each time step
!    noise_read_add_surface_movie is implemented later, which utilizes 'src/shared/write_c_binary.c' for faster I/O,
!                                                   which creates one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage

!!!!! improved version !!!!!
  subroutine noise_read_add_surface_movie(nmovie_points,accel_crust_mantle, &
                  normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                  ibelm_top_crust_mantle,ibool_crust_mantle, &
                  nspec_top,noise_surface_movie, &
                  it,jacobian2D_top_crust_mantle,wgllwgll_xy)
  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: nspec_top,it,nmovie_points
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_CM) :: jacobian2D_top_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle ! both input and output
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise, mask_noise
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  ! output parameters
  ! local parameters
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,nspec_top) :: noise_surface_movie

  ! read surface movie
  call read_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*nspec_top,it)

  ! get coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

      k = NGLLZ

      ! loop on all the points inside the element
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1
          iglob = ibool_crust_mantle(i,j,k,ispec)

          eta = noise_surface_movie(1,i,j,ispec2D) * normal_x_noise(ipoin) + &
                noise_surface_movie(2,i,j,ispec2D) * normal_y_noise(ipoin) + &
                noise_surface_movie(3,i,j,ispec2D) * normal_z_noise(ipoin)


        accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + eta * mask_noise(ipoin) * normal_x_noise(ipoin) &
                                                      * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
        accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + eta * mask_noise(ipoin) * normal_y_noise(ipoin) &
                                                      * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
        accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + eta * mask_noise(ipoin) * normal_z_noise(ipoin) &
                                                      * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
      enddo
    enddo

  enddo

  end subroutine noise_read_add_surface_movie

!!!!! original implementation, not used anymore (but kept here for references) !!!!!
!  subroutine noise_read_add_surface_movie_original(myrank,nmovie_points,accel_crust_mantle, &
!                  normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
!                  store_val_ux,store_val_uy,store_val_uz, &
!                  ibelm_top_crust_mantle,ibool_crust_mantle,nspec_top, &
!                  NIT,it,LOCAL_PATH,jacobian2D_top_crust_mantle,wgllwgll_xy)
!  implicit none
!  include 'mpif.h'
!  include "precision.h"
!  include "constants.h"
!  include "OUTPUT_FILES/values_from_mesher.h"
!  ! input parameters
!  integer :: myrank,nmovie_points,nspec_top,NIT,it
!  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_CM) :: jacobian2D_top_crust_mantle
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: accel_crust_mantle ! both input and output
!  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise, mask_noise
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
!  character(len=150) :: LOCAL_PATH
!  ! output parameters
!  ! local parameters
!  integer :: ipoin,ispec2D,ispec,i,j,k,iglob,ios
!  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: store_val_ux,store_val_uy,store_val_uz
!  real(kind=CUSTOM_REAL) :: eta
!  character(len=150) :: outputname
!
!
!  ! read surface movie
!  write(outputname,"('/proc',i6.6,'_surface_movie',i6.6)") myrank, it
!  open(unit=IIN_NOISE,file=trim(LOCAL_PATH)//outputname,status='old',form='unformatted',action='read',iostat=ios)
!  if( ios /= 0)  call exit_MPI(myrank,'file '//trim(outputname)//' does NOT exist!')
!  read(IIN_NOISE) store_val_ux
!  read(IIN_NOISE) store_val_uy
!  read(IIN_NOISE) store_val_uz
!  close(IIN_NOISE)
!
!  ! get coordinates of surface mesh and surface displacement
!  ipoin = 0
!  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
!    ispec = ibelm_top_crust_mantle(ispec2D)
!
!    k = NGLLZ
!
!    ! loop on all the points inside the element
!    do j = 1,NGLLY,NIT
!      do i = 1,NGLLX,NIT
!        ipoin = ipoin + 1
!        iglob = ibool_crust_mantle(i,j,k,ispec)
!
!        eta = store_val_ux(ipoin) * normal_x_noise(ipoin) + &
!              store_val_uy(ipoin) * normal_y_noise(ipoin) + &
!              store_val_uz(ipoin) * normal_z_noise(ipoin)
!
!        accel_crust_mantle(1,iglob) = accel_crust_mantle(1,iglob) + eta * mask_noise(ipoin) * normal_x_noise(ipoin) &
!                                                      * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
!        accel_crust_mantle(2,iglob) = accel_crust_mantle(2,iglob) + eta * mask_noise(ipoin) * normal_y_noise(ipoin) &
!                                                      * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
!        accel_crust_mantle(3,iglob) = accel_crust_mantle(3,iglob) + eta * mask_noise(ipoin) * normal_z_noise(ipoin) &
!                                                      * wgllwgll_xy(i,j) * jacobian2D_top_crust_mantle(i,j,ispec2D)
!      enddo
!    enddo
!
!  enddo
!
!  end subroutine noise_read_add_surface_movie_original

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 3: constructing noise source strength kernel
!
! there are two subroutines --- compute_kernels_strength_noise_original & compute_kernels_strength_noise
!    compute_kernels_strength_noise_original is implemented at first, which creates one file at each time step
!    compute_kernels_strength_noise is implemented later, which utilizes 'src/shared/write_c_binary.c' for faster I/O,
!                                                         which creates only one file for the all time steps
!
! by this modification, the efficiency is greatly improved
! and now, it should be OK to run NOISE_TOMOGRAPHY on a cluster with global storage

!!!!! improved version !!!!!
  subroutine compute_kernels_strength_noise(nmovie_points,ibool_crust_mantle, &
                          Sigma_kl_crust_mantle,displ_crust_mantle,deltat,it, &
                          normal_x_noise,normal_y_noise,normal_z_noise, &
                          nspec_top,noise_surface_movie, &
                          ibelm_top_crust_mantle)
  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer :: it,nspec_top,nmovie_points
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL) :: deltat
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: displ_crust_mantle
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise
  ! output parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    Sigma_kl_crust_mantle
  ! local parameters
  integer :: i,j,k,ispec,iglob,ipoin,ispec2D
  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,nspec_top) :: noise_surface_movie

  ! read surface movie, needed for Sigma_kl_crust_mantle
  call read_abs(9,noise_surface_movie,CUSTOM_REAL*NDIM*NGLLX*NGLLY*nspec_top,it)

  ! noise source strength kernel
  ! to keep similar structure to other kernels, the source strength kernel is saved as a volumetric kernel
  ! but only updated at the surface, because the noise is generated there
  ipoin = 0
  do ispec2D = 1, nspec_top
    ispec = ibelm_top_crust_mantle(ispec2D)

    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY
      do i = 1,NGLLX
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)

        eta = noise_surface_movie(1,i,j,ispec2D) * normal_x_noise(ipoin) + &
              noise_surface_movie(2,i,j,ispec2D) * normal_y_noise(ipoin) + &
              noise_surface_movie(3,i,j,ispec2D) * normal_z_noise(ipoin)

        Sigma_kl_crust_mantle(i,j,k,ispec) =  Sigma_kl_crust_mantle(i,j,k,ispec) &
           + deltat * eta * ( normal_x_noise(ipoin) * displ_crust_mantle(1,iglob) &
                            + normal_y_noise(ipoin) * displ_crust_mantle(2,iglob) &
                            + normal_z_noise(ipoin) * displ_crust_mantle(3,iglob) )
      enddo
    enddo

  enddo

  end subroutine compute_kernels_strength_noise

!!!!! original implementation, not used anymore (but kept here for references) !!!!!
!  subroutine compute_kernels_strength_noise_original(myrank,ibool_crust_mantle, &
!                          Sigma_kl_crust_mantle,displ_crust_mantle,deltat,it, &
!                          nmovie_points,normal_x_noise,normal_y_noise,normal_z_noise, &
!                          nspec_top,ibelm_top_crust_mantle,LOCAL_PATH)
!  implicit none
!  include "constants.h"
!  include "OUTPUT_FILES/values_from_mesher.h"
!  ! input parameters
!  integer :: myrank,nmovie_points,it,nspec_top
!  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
!  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
!  real(kind=CUSTOM_REAL) :: deltat
!  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: displ_crust_mantle
!  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise
!  character(len=150) :: LOCAL_PATH
!  ! output parameters
!  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
!    Sigma_kl_crust_mantle
!  ! local parameters
!  integer :: i,j,k,ispec,iglob,ipoin,ispec2D,ios
!  real(kind=CUSTOM_REAL) :: eta
!  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: store_val_ux,store_val_uy,store_val_uz
!  character(len=150) :: outputname
!
!
!  ! read surface movie, needed for Sigma_kl_crust_mantle
!  write(outputname,"('/proc',i6.6,'_surface_movie',i6.6)") myrank, it
!  open(unit=IIN_NOISE,file=trim(LOCAL_PATH)//outputname,status='old',form='unformatted',action='read',iostat=ios)
!  if( ios /= 0)  call exit_MPI(myrank,'file '//trim(outputname)//' does NOT exist!')
!
!  read(IIN_NOISE) store_val_ux
!  read(IIN_NOISE) store_val_uy
!  read(IIN_NOISE) store_val_uz
!  close(IIN_NOISE)
!
!  ! noise source strength kernel
!  ! to keep similar structure to other kernels, the source strength kernel is saved as a volumetric kernel
!  ! but only updated at the surface, because the noise is generated there
!  ipoin = 0
!  do ispec2D = 1, nspec_top
!    ispec = ibelm_top_crust_mantle(ispec2D)
!
!    k = NGLLZ
!
!    ! loop on all the points inside the element
!    do j = 1,NGLLY
!      do i = 1,NGLLX
!        ipoin = ipoin + 1
!        iglob = ibool_crust_mantle(i,j,k,ispec)
!
!        eta = store_val_ux(ipoin) * normal_x_noise(ipoin) + &
!              store_val_uy(ipoin) * normal_y_noise(ipoin) + &
!              store_val_uz(ipoin) * normal_z_noise(ipoin)
!
!        Sigma_kl_crust_mantle(i,j,k,ispec) =  Sigma_kl_crust_mantle(i,j,k,ispec) &
!           + deltat * eta * ( normal_x_noise(ipoin) * displ_crust_mantle(1,iglob) &
!                            + normal_y_noise(ipoin) * displ_crust_mantle(2,iglob) &
!                            + normal_z_noise(ipoin) * displ_crust_mantle(3,iglob) )
!      enddo
!    enddo
!
!  enddo
!
!  end subroutine compute_kernels_strength_noise_original

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 3: save noise source strength kernel
  subroutine save_kernels_strength_noise(myrank,LOCAL_PATH,Sigma_kl_crust_mantle)
  implicit none
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  ! input parameters
  integer myrank
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: Sigma_kl_crust_mantle
  character(len=150) :: LOCAL_PATH
  ! output parameters
  ! local parameters
  character(len=150) :: prname
  integer :: ier

  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

  open(unit=IOUT_NOISE,file=trim(prname)//'Sigma_kernel.bin', &
        status='unknown',form='unformatted',action='write',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error opening file Sigma_kernel.bin')

  write(IOUT_NOISE) Sigma_kl_crust_mantle     ! need to put dimensions back (not done yet)
  close(IOUT_NOISE)

  end subroutine save_kernels_strength_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================
