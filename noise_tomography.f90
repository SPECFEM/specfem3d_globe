!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
!                            March 2010
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

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! chracterize noise statistics
! for a given point (xcoord,ycoord,zcoord), specify the noise direction and noise distribution mask_noise
! USERS need to modify this subroutine for their own noise characteristics
  subroutine noise_distribution_direction(xcoord_in,ycoord_in,zcoord_in, &
                  normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                  mask_noise_out,NSTEP)
  implicit none
  include "constants.h"
  real(kind=CUSTOM_REAL) :: xcoord,ycoord,zcoord,normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  real(kind=CUSTOM_REAL) :: xcoord_in,ycoord_in,zcoord_in
  character(len=150) :: filename
  integer :: ios, itime, NSTEP
  real(kind=CUSTOM_REAL) :: noise_src(NSTEP)
  real(kind=CUSTOM_REAL) :: junk

  ! coordinates "x/y/zcoord" actually contain r theta phi, therefore convert back to x y z
  call rthetaphi_2_xyz(xcoord,ycoord,zcoord, xcoord_in,ycoord_in,zcoord_in)
  ! NOTE that all coordinates are non-dimensionalized in GLOBAL package!
  ! USERS are free to choose which set to use,
  ! either "r theta phi" (xcoord_in,ycoord_in,zcoord_in)
  ! or     "x y z"       (xcoord,ycoord,zcoord)

  !******************************** add your noise characteristics below ****************************************
  ! here, the noise is assumed to be vertical
  normal_x_noise_out = xcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  normal_y_noise_out = ycoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  normal_z_noise_out = zcoord / sqrt(xcoord**2 + ycoord**2 + zcoord**2)
  ! here, the noise is assumed to be uniform
  mask_noise_out = 1.0
  !******************************** add your noise characteristics above ****************************************

  ! here, the noise spectrum is based upon Peterson's model (Peterson, 1993)
  ! output the source time function based upon the noise's power spectrum
  ! NOT implemented yet. file "S_squared" is calculated before the simulation
  
  !filename = 'NOISE_TOMOGRAPHY/'//'S_squared'
  !if (myrank == 0) then
  !   open(unit=IOUT,file=trim(filename),status='unknown',action='write',iostat=ios)
  !   if( ios /= 0)  call exit_MPI(myrank, 'file '//trim(filename)//' cannot be created!')
  !   
  !   do itime =1,NSTEP
  !     WRITE(IOUT,*) itime, noise_src(itime)
  !   enddo
  !   close(IOUT)
  !endif

  end subroutine noise_distribution_direction

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! read parameters
  subroutine read_parameters_noise(myrank,nrec,NSTEP,nmovie_points, &
                                   islice_selected_rec,xi_receiver,eta_receiver,gamma_receiver,nu, &
                                   noise_sourcearray,xigll,yigll,zigll,nspec_top, & 
                                   NIT, ibool_crust_mantle, ibelm_top_crust_mantle, &
                                   xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                                   irec_master_noise,normal_x_noise,normal_y_noise,normal_z_noise,mask_noise)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise,mask_noise
  real(kind=CUSTOM_REAL) :: normal_x_noise_out,normal_y_noise_out,normal_z_noise_out,mask_noise_out
  integer, dimension(nrec) :: islice_selected_rec
  double precision, dimension(nrec) :: xi_receiver,eta_receiver,gamma_receiver
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob,ios,irec_master_noise,myrank,nmovie_points,nrec,NSTEP
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer nspec_top,NIT
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle
  double precision, dimension(NDIM,NDIM,nrec) :: nu
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll
  character(len=150) :: filename
  
  ! read master receiver ID -- the ID in DATA/STATIONS
  filename = 'NOISE_TOMOGRAPHY/'//'irec_master_noise'
  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0)  &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file contains the ID of the master receiver')
  read(IIN,*,iostat=ios) irec_master_noise
  close(IIN)
  if (myrank == 0) then
     WRITE(IMAIN,*) 'The master receiver is: (RECEIVER ID)', irec_master_noise
  endif

  ! compute source arrays for "ensemble forward source", which is source of "ensemble forward wavefield"
  if(myrank == islice_selected_rec(irec_master_noise) .OR. myrank==0) then ! myrank==0 is used for output only
    call compute_arrays_source_noise(myrank, &
              xi_receiver(irec_master_noise),eta_receiver(irec_master_noise),gamma_receiver(irec_master_noise), &
              nu(:,:,irec_master_noise),noise_sourcearray, xigll,yigll,zigll,NSTEP)
  endif

  ! noise distribution and noise direction
  ipoin = 0
  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ! this subroutine must be provided by USERS
        call noise_distribution_direction(xstore_crust_mantle(iglob), &
                  ystore_crust_mantle(iglob),zstore_crust_mantle(iglob), &
                  normal_x_noise_out,normal_y_noise_out,normal_z_noise_out, &
                  mask_noise_out,NSTEP)
        normal_x_noise(ipoin) = normal_x_noise_out
        normal_y_noise(ipoin) = normal_y_noise_out
        normal_z_noise(ipoin) = normal_z_noise_out
        mask_noise(ipoin)     = mask_noise_out
      enddo
    enddo

  enddo

  end subroutine read_parameters_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! check for consistency of parameters
  subroutine check_parameters_noise(myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                              NUMBER_OF_RUNS, NUMBER_OF_THIS_RUN,ROTATE_SEISMOGRAMS_RT, &
                                              SAVE_ALL_SEISMOS_IN_ONE_FILE, USE_BINARY_FOR_LARGE_FILE)

  implicit none

  include 'mpif.h'
  include "precision.h"
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"
  integer myrank,NOISE_TOMOGRAPHY,SIMULATION_TYPE,NUMBER_OF_RUNS,NUMBER_OF_THIS_RUN
  logical SAVE_FORWARD,ROTATE_SEISMOGRAMS_RT,SAVE_ALL_SEISMOS_IN_ONE_FILE, USE_BINARY_FOR_LARGE_FILE


  WRITE(IMAIN,*) '*******************************************************************************'
  WRITE(IMAIN,*) '*******************************************************************************'
  WRITE(IMAIN,*) 'IMPORTANT WARNING'
  WRITE(IMAIN,*) 'You are running simulations using NOISE TOMOGRAPHY techniques.'
  WRITE(IMAIN,*) 'Please make sure you understand the procedures before you have a try.'
  WRITE(IMAIN,*) 'Displacements everywhere at the free surface are saved every timestep,'
  WRITE(IMAIN,*) 'so make sure that LOCAL_PATH in DATA/Par_file is not global.'
  WRITE(IMAIN,*) 'Otherwise the disk storage may be a serious issue, as is the I/O.'
  WRITE(IMAIN,*) 'Also make sure that NO earthquakes are included,'
  WRITE(IMAIN,*) 'i.e., set moment tensor to be ZERO in CMTSOLUTION'
  WRITE(IMAIN,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  WRITE(IMAIN,*) 'If you just want a regular EARTHQUAKE simulation,'
  WRITE(IMAIN,*) 'set NOISE_TOMOGRAPHY=0 in DATA/Par_file'
  WRITE(IMAIN,*) '*******************************************************************************'
  WRITE(IMAIN,*) '*******************************************************************************'

  if (NUMBER_OF_RUNS/=1 .OR. NUMBER_OF_THIS_RUN/=1) &
     call exit_mpi(myrank,'NUMBER_OF_RUNS and NUMBER_OF_THIS_RUN must be 1 for NOISE TOMOGRAPHY! check DATA/Par_file')
  if (ROTATE_SEISMOGRAMS_RT) &
     call exit_mpi(myrank,'Do NOT rotate seismograms in the code, change ROTATE_SEISMOGRAMS_RT in DATA/Par_file')
  if (SAVE_ALL_SEISMOS_IN_ONE_FILE .OR. USE_BINARY_FOR_LARGE_FILE) &
     call exit_mpi(myrank,'Please set SAVE_ALL_SEISMOS_IN_ONE_FILE and USE_BINARY_FOR_LARGE_FILE to be .false.')

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
  end subroutine check_parameters_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! read and construct the "source" (source time function based upon noise spectrum) for "ensemble forward source"
  subroutine compute_arrays_source_noise(myrank, &
      xi_noise,eta_noise,gamma_noise, nu,noise_sourcearray, &
      xigll,yigll,zigll,NSTEP)

  implicit none

  include 'constants.h'
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank, NSTEP
  double precision xi_noise, eta_noise, gamma_noise

  ! output
  real(kind=CUSTOM_REAL) :: noise_sourcearray(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP)
  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLY) :: yigll
  double precision, dimension(NGLLZ) :: zigll

  double precision, dimension(NDIM,NDIM) :: nu
  double precision, dimension(NDIM) :: nu_master

  double precision,parameter :: scale_displ_inv = 1.d0/R_EARTH

  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
        hgammar(NGLLZ), hpgammar(NGLLZ)
  real(kind=CUSTOM_REAL) :: noise_src(NSTEP),noise_src_u(NDIM,NSTEP)

  integer itime, i, j, k, ios
  real(kind=CUSTOM_REAL) :: junk
  character(len=150) :: filename

  noise_src = 0._CUSTOM_REAL

  ! noise file
  filename = 'NOISE_TOMOGRAPHY/'//'S_squared'
  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0)  &
    call exit_MPI(myrank, 'file '//trim(filename)//' does NOT exist! This file is based upon noise spectrum')
  
  do itime =1,NSTEP
    read(IIN,*,iostat=ios) junk, noise_src(itime)
    if( ios /= 0)  call exit_MPI(myrank,&
        'file '//trim(filename)//' has wrong length, please check your simulation duration')
  enddo
  close(IIN)

  ! master receiver component direction, \nu_master
  filename = 'NOISE_TOMOGRAPHY/'//'nu_master'
  open(unit=IIN,file=trim(filename),status='old',action='read',iostat=ios)
  if( ios /= 0)  call exit_MPI(myrank,&
        'file '//trim(filename)//' does NOT exist! nu_master is the component direction for master receiver')
  
  do itime =1,3
    read(IIN,*,iostat=ios) nu_master(itime)
    if( ios /= 0) call exit_MPI(myrank,&
        'file '//trim(filename)//' has wrong length, the vector should have three components (ENZ)')
  enddo
  close(IIN)

  if (myrank == 0) &
     WRITE(IMAIN,*) 'The direction (ENZ) of selected component of master receiver is', nu_master
  
  ! rotates to cartesian
  do itime = 1, NSTEP
    noise_src_u(:,itime) = nu(1,:) * noise_src(itime) * nu_master(1) &
                         + nu(2,:) * noise_src(itime) * nu_master(2) &
                         + nu(3,:) * noise_src(itime) * nu_master(3)
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

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

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

  integer myrank,nrec,NSTEP

  real(kind=CUSTOM_REAL),dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    accel_crust_mantle

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(nrec) :: islice_selected_rec,ispec_selected_rec

  integer :: it

  ! local parameters
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NGLLZ,NSTEP) :: noise_sourcearray
  integer :: irec_master_noise,i,j,k,iglob

    ! adds noise source (only if this proc carries the noise)
    if(myrank == islice_selected_rec(irec_master_noise)) then
      ! adds nosie source contributions
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

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 1: calculate the "ensemble forward source"
! save surface movie (displacement) at every time steps, for step 2.
  subroutine noise_save_surface_movie(myrank,nmovie_points,displ_crust_mantle, &
                    xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle, &
                    store_val_x,store_val_y,store_val_z, &
                    store_val_ux,store_val_uy,store_val_uz, &
                    ibelm_top_crust_mantle,ibool_crust_mantle,nspec_top, &
                    NIT,it,LOCAL_PATH)

  implicit none

  include 'mpif.h'
  include "precision.h"
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,nmovie_points

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     displ_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: &
        xstore_crust_mantle,ystore_crust_mantle,zstore_crust_mantle

  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: &
      store_val_x,store_val_y,store_val_z, &
      store_val_ux,store_val_uy,store_val_uz

  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  integer nspec_top,NIT,it
  character(len=150) LOCAL_PATH

  ! local parameters
  character(len=150) :: outputname
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob

  ! get coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
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
        store_val_ux(ipoin) = displ_crust_mantle(1,iglob)
        store_val_uy(ipoin) = displ_crust_mantle(2,iglob)
        store_val_uz(ipoin) = displ_crust_mantle(3,iglob)
      enddo
    enddo

  enddo

  ! save surface motion to disk 
  ! LOCAL storage is better than GLOBAL, because we have to save the 'movie' at every time step
  ! also note that the surface movie needs not to be shared with other nodes/procs
  ! change LOCAL_PATH specified in "DATA/Par_file"
    write(outputname,"('/proc',i6.6,'_surface_movie',i6.6)") myrank, it
    open(unit=IOUT,file=trim(LOCAL_PATH)//outputname,status='unknown',form='unformatted',action='write')
    write(IOUT) store_val_ux
    write(IOUT) store_val_uy
    write(IOUT) store_val_uz
    close(IOUT)

  end subroutine noise_save_surface_movie 

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 2/3: calculate/reconstructe the "ensemble forward wavefield"
! read surface movie (displacement) at every time steps, injected as the source of "ensemble forward wavefield"
! in step 2, call noise_read_add_surface_movie(..., NSTEP-it+1 ,...)
! in step 3, call noise_read_add_surface_movie(..., it ,...)
  subroutine noise_read_add_surface_movie(myrank,nmovie_points,accel_crust_mantle, &
                  normal_x_noise,normal_y_noise,normal_z_noise,mask_noise, &
                  store_val_ux,store_val_uy,store_val_uz, &
                  ibelm_top_crust_mantle,ibool_crust_mantle,nspec_top, &
                  NIT,it,LOCAL_PATH,jacobian2D_top_crust_mantle,wgllwgll_xy)  

  implicit none

  include 'mpif.h'
  include "precision.h"
  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,nmovie_points

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
     accel_crust_mantle

  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: &
      store_val_ux,store_val_uy,store_val_uz

  real(kind=CUSTOM_REAL) :: eta
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise, mask_noise

  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_TOP_CM) :: jacobian2D_top_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy

  integer nspec_top,NIT,it

  character(len=150) :: outputname, LOCAL_PATH
  integer :: ipoin,ispec2D,ispec,i,j,k,iglob,ios


  ! read surface movie
  write(outputname,"('/proc',i6.6,'_surface_movie',i6.6)") myrank, it
  open(unit=IOUT,file=trim(LOCAL_PATH)//outputname,status='old',form='unformatted',action='read',iostat=ios)
  if( ios /= 0)  call exit_MPI(myrank,'file '//trim(outputname)//' does NOT exist!')
  read(IOUT) store_val_ux
  read(IOUT) store_val_uy
  read(IOUT) store_val_uz
  close(IOUT)

  ! get coordinates of surface mesh and surface displacement
  ipoin = 0
  do ispec2D = 1, nspec_top ! NSPEC2D_TOP(IREGION_CRUST_MANTLE)
    ispec = ibelm_top_crust_mantle(ispec2D)

    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY,NIT
      do i = 1,NGLLX,NIT
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)

        eta = store_val_ux(ipoin) * normal_x_noise(ipoin) + &
              store_val_uy(ipoin) * normal_y_noise(ipoin) + &
              store_val_uz(ipoin) * normal_z_noise(ipoin) 

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

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 3: constructing noise source strength kernel
  subroutine compute_kernels_strength_noise(myrank,ibool_crust_mantle, &
                          Sigma_kl_crust_mantle,displ_crust_mantle,deltat,it, &
                          nmovie_points,normal_x_noise,normal_y_noise,normal_z_noise, &
                          nspec_top,ibelm_top_crust_mantle,LOCAL_PATH)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle
  integer, dimension(NSPEC2D_TOP_CM) :: ibelm_top_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    Sigma_kl_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    displ_crust_mantle

  real(kind=CUSTOM_REAL) deltat,eta

  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: normal_x_noise,normal_y_noise,normal_z_noise
  real(kind=CUSTOM_REAL), dimension(nmovie_points) :: store_val_ux,store_val_uy,store_val_uz

  integer :: i,j,k,ispec,iglob, nmovie_points,it,nspec_top,ipoin,ispec2D,ios,myrank
  character(len=150) outputname,LOCAL_PATH

  ! read surface movie, needed for Sigma_kl_crust_mantle
  write(outputname,"('/proc',i6.6,'_surface_movie',i6.6)") myrank, it
  open(unit=IOUT,file=trim(LOCAL_PATH)//outputname,status='old',form='unformatted',action='read',iostat=ios)
  if( ios /= 0)  call exit_MPI(myrank,'file '//trim(outputname)//' does NOT exist!')

  read(IOUT) store_val_ux
  read(IOUT) store_val_uy
  read(IOUT) store_val_uz
  close(IOUT)

  ! noise source strength kernel (only updated at the surface, because the noise is generated there)
  ipoin = 0
  do ispec2D = 1, nspec_top 
    ispec = ibelm_top_crust_mantle(ispec2D)

    k = NGLLZ

    ! loop on all the points inside the element
    do j = 1,NGLLY
      do i = 1,NGLLX
        ipoin = ipoin + 1
        iglob = ibool_crust_mantle(i,j,k,ispec)

        eta = store_val_ux(ipoin) * normal_x_noise(ipoin) + &
              store_val_uy(ipoin) * normal_y_noise(ipoin) + &
              store_val_uz(ipoin) * normal_z_noise(ipoin) 

        Sigma_kl_crust_mantle(i,j,k,ispec) =  Sigma_kl_crust_mantle(i,j,k,ispec) &
           + deltat * eta * ( normal_x_noise(ipoin) * displ_crust_mantle(1,iglob) &
                            + normal_y_noise(ipoin) * displ_crust_mantle(2,iglob) &
                            + normal_z_noise(ipoin) * displ_crust_mantle(3,iglob) )
      enddo
    enddo

  enddo
  end subroutine compute_kernels_strength_noise

! =============================================================================================================
! =============================================================================================================
! =============================================================================================================

! subroutine for NOISE TOMOGRAPHY
! step 3: save noise source strength kernel
  subroutine save_kernels_strength_noise(myrank,LOCAL_PATH, &
                                        Sigma_kl_crust_mantle,scale_t,scale_displ)
  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  character(len=150) prname,LOCAL_PATH
  double precision :: scale_t,scale_displ,scale_Sigma_kl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_ADJOINT) :: &
    Sigma_kl_crust_mantle
  integer myrank

  call create_name_database(prname,myrank,IREGION_CRUST_MANTLE,LOCAL_PATH)

    open(unit=27,file=trim(prname)//'Sigma_kernel.bin',status='unknown',form='unformatted',action='write')
    write(27) Sigma_kl_crust_mantle
    close(27)
  end subroutine save_kernels_strength_noise
