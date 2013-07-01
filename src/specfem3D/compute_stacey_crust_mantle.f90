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

  subroutine compute_stacey_crust_mantle()

  use constants_solver

  use specfem_par,only: &
    ichunk,SIMULATION_TYPE,SAVE_FORWARD,NSTEP,it, &
    wgllwgll_xz,wgllwgll_yz

  use specfem_par,only: GPU_MODE,Mesh_pointer

  use specfem_par_crustmantle, only: &
    veloc_crust_mantle,accel_crust_mantle,b_accel_crust_mantle, &
    ibool_crust_mantle, &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
    normal_xmin_crust_mantle,normal_xmax_crust_mantle, &
    normal_ymin_crust_mantle,normal_ymax_crust_mantle, &
    rho_vp_crust_mantle,rho_vs_crust_mantle, &
    ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle, &
    ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle, &
    nimin_crust_mantle,nimax_crust_mantle, &
    njmin_crust_mantle,njmax_crust_mantle, &
    nkmin_xi_crust_mantle,nkmin_eta_crust_mantle, &
    nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle, &
    reclen_xmin_crust_mantle,reclen_xmax_crust_mantle, &
    reclen_ymin_crust_mantle,reclen_ymax_crust_mantle, &
    absorb_xmin_crust_mantle,absorb_xmax_crust_mantle, &
    absorb_ymin_crust_mantle,absorb_ymax_crust_mantle

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: weight
  real(kind=CUSTOM_REAL) :: vn,vx,vy,vz,nx,ny,nz,tx,ty,tz
  integer :: i,j,k,ispec,iglob,ispec2D
  !integer :: reclen1,reclen2

  ! note: we use c functions for I/O as they still have a better performance than
  !           fortran, unformatted file I/O. however, using -assume byterecl together with fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storages (every 8 step) and/or asynchronious
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)


  ! crust & mantle

  !   xmin
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then

    ! reads absorbing boundary values
    if (SIMULATION_TYPE == 3 .and. nspec2D_xmin_crust_mantle > 0)  then
      ! note: backward/reconstructed wavefields are read in after the Newmark time scheme in the first time loop
      !          this leads to a corresponding boundary condition at time index NSTEP - (it-1) = NSTEP - it + 1
      call read_abs(0,absorb_xmin_crust_mantle,reclen_xmin_crust_mantle,NSTEP-it+1)
    endif

    if ( .NOT. GPU_MODE) then
      ! on CPU
      do ispec2D=1,nspec2D_xmin_crust_mantle

        ispec=ibelm_xmin_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if(nkmin_xi_crust_mantle(1,ispec2D) == 0 .or. njmin_crust_mantle(1,ispec2D) == 0) cycle

        i=1
        do k=nkmin_xi_crust_mantle(1,ispec2D),NGLLZ
          do j=njmin_crust_mantle(1,ispec2D),njmax_crust_mantle(1,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            vx=veloc_crust_mantle(1,iglob)
            vy=veloc_crust_mantle(2,iglob)
            vz=veloc_crust_mantle(3,iglob)

            nx=normal_xmin_crust_mantle(1,j,k,ispec2D)
            ny=normal_xmin_crust_mantle(2,j,k,ispec2D)
            nz=normal_xmin_crust_mantle(3,j,k,ispec2D)

            vn=vx*nx+vy*ny+vz*nz

            tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
            ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
            tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

            weight=jacobian2D_xmin_crust_mantle(j,k,ispec2D)*wgllwgll_yz(j,k)

            accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
            accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
            accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

            if (SIMULATION_TYPE == 3) then
              b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_xmin_crust_mantle(:,j,k,ispec2D)
            else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              absorb_xmin_crust_mantle(1,j,k,ispec2D) = tx*weight
              absorb_xmin_crust_mantle(2,j,k,ispec2D) = ty*weight
              absorb_xmin_crust_mantle(3,j,k,ispec2D) = tz*weight
            endif
          enddo
        enddo
      enddo

    else
      ! on GPU
      if( nspec2D_xmin_crust_mantle > 0 ) call compute_stacey_elastic_cuda(Mesh_pointer, &
                                                                absorb_xmin_crust_mantle, &
                                                                0) ! <= xmin
    endif

    ! writes absorbing boundary values
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_xmin_crust_mantle > 0 ) then
      call write_abs(0,absorb_xmin_crust_mantle, reclen_xmin_crust_mantle,it)
    endif

  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC

  !   xmax
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then

    ! reads absorbing boundary values
    if (SIMULATION_TYPE == 3 .and. nspec2D_xmax_crust_mantle > 0)  then
      call read_abs(1,absorb_xmax_crust_mantle,reclen_xmax_crust_mantle,NSTEP-it+1)
    endif

    if(.NOT. GPU_MODE ) then
      ! on CPU
      do ispec2D=1,nspec2D_xmax_crust_mantle

        ispec=ibelm_xmax_crust_mantle(ispec2D)

        ! exclude elements that are not on absorbing edges
        if(nkmin_xi_crust_mantle(2,ispec2D) == 0 .or. njmin_crust_mantle(2,ispec2D) == 0) cycle

        i=NGLLX
        do k=nkmin_xi_crust_mantle(2,ispec2D),NGLLZ
          do j=njmin_crust_mantle(2,ispec2D),njmax_crust_mantle(2,ispec2D)
            iglob=ibool_crust_mantle(i,j,k,ispec)

            vx=veloc_crust_mantle(1,iglob)
            vy=veloc_crust_mantle(2,iglob)
            vz=veloc_crust_mantle(3,iglob)

            nx=normal_xmax_crust_mantle(1,j,k,ispec2D)
            ny=normal_xmax_crust_mantle(2,j,k,ispec2D)
            nz=normal_xmax_crust_mantle(3,j,k,ispec2D)

            vn=vx*nx+vy*ny+vz*nz

            tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
            ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
            tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

            weight=jacobian2D_xmax_crust_mantle(j,k,ispec2D)*wgllwgll_yz(j,k)

            accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
            accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
            accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

            if (SIMULATION_TYPE == 3) then
              b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_xmax_crust_mantle(:,j,k,ispec2D)
            else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
              absorb_xmax_crust_mantle(1,j,k,ispec2D) = tx*weight
              absorb_xmax_crust_mantle(2,j,k,ispec2D) = ty*weight
              absorb_xmax_crust_mantle(3,j,k,ispec2D) = tz*weight
            endif

          enddo
        enddo
      enddo

    else
      ! on GPU
      if( nspec2D_xmax_crust_mantle > 0 ) call compute_stacey_elastic_cuda(Mesh_pointer, &
                                                                absorb_xmax_crust_mantle, &
                                                                1) ! <= xmin
    endif


    ! writes absorbing boundary values
    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_xmax_crust_mantle > 0 ) then
      call write_abs(1,absorb_xmax_crust_mantle,reclen_xmax_crust_mantle,it)
    endif

  endif ! NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB

  !   ymin

  ! reads absorbing boundary values
  if (SIMULATION_TYPE == 3 .and. nspec2D_ymin_crust_mantle > 0)  then
    call read_abs(2,absorb_ymin_crust_mantle, reclen_ymin_crust_mantle,NSTEP-it+1)
  endif

  if( .NOT. GPU_MODE ) then
    ! on CPU
    do ispec2D=1,nspec2D_ymin_crust_mantle

      ispec=ibelm_ymin_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_eta_crust_mantle(1,ispec2D) == 0 .or. nimin_crust_mantle(1,ispec2D) == 0) cycle

      j=1
      do k=nkmin_eta_crust_mantle(1,ispec2D),NGLLZ
        do i=nimin_crust_mantle(1,ispec2D),nimax_crust_mantle(1,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          vx=veloc_crust_mantle(1,iglob)
          vy=veloc_crust_mantle(2,iglob)
          vz=veloc_crust_mantle(3,iglob)

          nx=normal_ymin_crust_mantle(1,i,k,ispec2D)
          ny=normal_ymin_crust_mantle(2,i,k,ispec2D)
          nz=normal_ymin_crust_mantle(3,i,k,ispec2D)

          vn=vx*nx+vy*ny+vz*nz

          tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
          ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
          tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

          weight=jacobian2D_ymin_crust_mantle(i,k,ispec2D)*wgllwgll_xz(i,k)

          accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
          accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
          accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

          if (SIMULATION_TYPE == 3) then
            b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_ymin_crust_mantle(:,i,k,ispec2D)
          else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
            absorb_ymin_crust_mantle(1,i,k,ispec2D) = tx*weight
            absorb_ymin_crust_mantle(2,i,k,ispec2D) = ty*weight
            absorb_ymin_crust_mantle(3,i,k,ispec2D) = tz*weight
          endif

        enddo
      enddo
    enddo

  else
    ! on GPU
    if( nspec2D_ymin_crust_mantle > 0 ) call compute_stacey_elastic_cuda(Mesh_pointer, &
                                                              absorb_ymin_crust_mantle, &
                                                              2) ! <= ymin
  endif


  ! writes absorbing boundary values
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_ymin_crust_mantle > 0 ) then
    call write_abs(2,absorb_ymin_crust_mantle,reclen_ymin_crust_mantle,it)
  endif



  !   ymax

  ! reads absorbing boundary values
  if (SIMULATION_TYPE == 3 .and. nspec2D_ymax_crust_mantle > 0)  then
    call read_abs(3,absorb_ymax_crust_mantle,reclen_ymax_crust_mantle,NSTEP-it+1)
  endif

  if( .NOT. GPU_MODE ) then
    ! on CPU
    do ispec2D=1,nspec2D_ymax_crust_mantle

      ispec=ibelm_ymax_crust_mantle(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_eta_crust_mantle(2,ispec2D) == 0 .or. nimin_crust_mantle(2,ispec2D) == 0) cycle

      j=NGLLY
      do k=nkmin_eta_crust_mantle(2,ispec2D),NGLLZ
        do i=nimin_crust_mantle(2,ispec2D),nimax_crust_mantle(2,ispec2D)
          iglob=ibool_crust_mantle(i,j,k,ispec)

          vx=veloc_crust_mantle(1,iglob)
          vy=veloc_crust_mantle(2,iglob)
          vz=veloc_crust_mantle(3,iglob)

          nx=normal_ymax_crust_mantle(1,i,k,ispec2D)
          ny=normal_ymax_crust_mantle(2,i,k,ispec2D)
          nz=normal_ymax_crust_mantle(3,i,k,ispec2D)

          vn=vx*nx+vy*ny+vz*nz

          tx=rho_vp_crust_mantle(i,j,k,ispec)*vn*nx+rho_vs_crust_mantle(i,j,k,ispec)*(vx-vn*nx)
          ty=rho_vp_crust_mantle(i,j,k,ispec)*vn*ny+rho_vs_crust_mantle(i,j,k,ispec)*(vy-vn*ny)
          tz=rho_vp_crust_mantle(i,j,k,ispec)*vn*nz+rho_vs_crust_mantle(i,j,k,ispec)*(vz-vn*nz)

          weight=jacobian2D_ymax_crust_mantle(i,k,ispec2D)*wgllwgll_xz(i,k)

          accel_crust_mantle(1,iglob)=accel_crust_mantle(1,iglob) - tx*weight
          accel_crust_mantle(2,iglob)=accel_crust_mantle(2,iglob) - ty*weight
          accel_crust_mantle(3,iglob)=accel_crust_mantle(3,iglob) - tz*weight

          if (SIMULATION_TYPE == 3) then
            b_accel_crust_mantle(:,iglob)=b_accel_crust_mantle(:,iglob) - absorb_ymax_crust_mantle(:,i,k,ispec2D)
          else if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD) then
            absorb_ymax_crust_mantle(1,i,k,ispec2D) = tx*weight
            absorb_ymax_crust_mantle(2,i,k,ispec2D) = ty*weight
            absorb_ymax_crust_mantle(3,i,k,ispec2D) = tz*weight
          endif

        enddo
      enddo
    enddo

  else
    ! on GPU
    if( nspec2D_ymax_crust_mantle > 0 ) call compute_stacey_elastic_cuda(Mesh_pointer, &
                                                              absorb_ymax_crust_mantle, &
                                                              3) ! <= ymax
  endif

  ! writes absorbing boundary values
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_ymax_crust_mantle > 0 ) then
    call write_abs(3,absorb_ymax_crust_mantle,reclen_ymax_crust_mantle,it)
  endif

  end subroutine compute_stacey_crust_mantle

