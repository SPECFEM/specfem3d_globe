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

  subroutine compute_stacey_crust_mantle(myrank,ichunk,SIMULATION_TYPE, &
                              NSTEP,it,SAVE_FORWARD,ibool_crust_mantle, &
                              veloc_crust_mantle,accel_crust_mantle,b_accel_crust_mantle, &
                              jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle, &
                              jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle, &
                              wgllwgll_xz,wgllwgll_yz, &
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
                              nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm, &
                              absorb_xmin_crust_mantle,absorb_xmax_crust_mantle, &
                              absorb_ymin_crust_mantle,absorb_ymax_crust_mantle)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer myrank,ichunk,SIMULATION_TYPE
  integer NSTEP,it
  logical SAVE_FORWARD

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE) :: ibool_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE) :: &
    veloc_crust_mantle,accel_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLOB_CRUST_MANTLE_ADJOINT) :: &
    b_accel_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    jacobian2D_xmin_crust_mantle,jacobian2D_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    jacobian2D_ymin_crust_mantle,jacobian2D_ymax_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_CM) :: &
  normal_xmin_crust_mantle,normal_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLY,NSPEC2DMAX_YMIN_YMAX_CM) :: &
  normal_ymin_crust_mantle,normal_ymax_crust_mantle

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE_STACEY) :: &
    rho_vp_crust_mantle,rho_vs_crust_mantle

  integer, dimension(NSPEC2DMAX_XMIN_XMAX_CM) :: ibelm_xmin_crust_mantle,ibelm_xmax_crust_mantle
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_CM) :: ibelm_ymin_crust_mantle,ibelm_ymax_crust_mantle

  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_CM) :: &
    nimin_crust_mantle,nimax_crust_mantle,nkmin_eta_crust_mantle
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_CM) :: &
    njmin_crust_mantle,njmax_crust_mantle,nkmin_xi_crust_mantle

  integer nspec2D_xmin_crust_mantle,nspec2D_xmax_crust_mantle, &
    nspec2D_ymin_crust_mantle,nspec2D_ymax_crust_mantle

  integer reclen_xmin_crust_mantle,reclen_xmax_crust_mantle,&
    reclen_ymin_crust_mantle,reclen_ymax_crust_mantle

  integer nabs_xmin_cm,nabs_xmax_cm,nabs_ymin_cm,nabs_ymax_cm
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nabs_xmin_cm) :: absorb_xmin_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLY,NGLLZ,nabs_xmax_cm) :: absorb_xmax_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,nabs_ymin_cm) :: absorb_ymin_crust_mantle
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ,nabs_ymax_cm) :: absorb_ymax_crust_mantle


  ! local parameters
  real(kind=CUSTOM_REAL) :: weight
  real(kind=CUSTOM_REAL) :: vn,vx,vy,vz,nx,ny,nz,tx,ty,tz
  integer :: i,j,k,ispec,iglob,ispec2D
  integer :: reclen1,reclen2


  ! crust & mantle

  !   xmin
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then

    if (SIMULATION_TYPE == 3 .and. nspec2D_xmin_crust_mantle > 0)  then
      read(51,rec=NSTEP-it+1) reclen1,absorb_xmin_crust_mantle,reclen2
      if (reclen1 /= reclen_xmin_crust_mantle .or. reclen1 /= reclen2)  &
         call exit_MPI(myrank,'Error reading absorbing contribution absorb_xmin')
    endif

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

    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_xmin_crust_mantle > 0 ) &
               write(51,rec=it) reclen_xmin_crust_mantle,absorb_xmin_crust_mantle,reclen_xmin_crust_mantle

  endif

  !   xmax
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then

    if (SIMULATION_TYPE == 3 .and. nspec2D_xmax_crust_mantle > 0)  then
      read(52,rec=NSTEP-it+1) reclen1,absorb_xmax_crust_mantle,reclen2
      if (reclen1 /= reclen_xmax_crust_mantle .or. reclen1 /= reclen2)  &
         call exit_MPI(myrank,'Error reading absorbing contribution absorb_xmax')
    endif

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

    if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_xmax_crust_mantle > 0 ) &
               write(52,rec=it) reclen_xmax_crust_mantle,absorb_xmax_crust_mantle,reclen_xmax_crust_mantle

  endif

  !   ymin

  if (SIMULATION_TYPE == 3 .and. nspec2D_ymin_crust_mantle > 0)  then
    read(53,rec=NSTEP-it+1) reclen1,absorb_ymin_crust_mantle,reclen2
    if (reclen1 /= reclen_ymin_crust_mantle .or. reclen1 /= reclen2)  &
       call exit_MPI(myrank,'Error reading absorbing contribution absorb_ymin')
  endif

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
  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_ymin_crust_mantle > 0 ) &
     write(53,rec=it) reclen_ymin_crust_mantle,absorb_ymin_crust_mantle,reclen_ymin_crust_mantle

  !   ymax

  if (SIMULATION_TYPE == 3 .and. nspec2D_ymax_crust_mantle > 0)  then
    read(54,rec=NSTEP-it+1) reclen1,absorb_ymax_crust_mantle,reclen2
    if (reclen1 /= reclen_ymax_crust_mantle .or. reclen1 /= reclen2)  &
       call exit_MPI(myrank,'Error reading absorbing contribution absorb_ymax')
  endif

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

  if (SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. nspec2D_ymax_crust_mantle > 0 ) &
     write(54,rec=it) reclen_ymax_crust_mantle,absorb_ymax_crust_mantle,reclen_ymax_crust_mantle

  end subroutine compute_stacey_crust_mantle

