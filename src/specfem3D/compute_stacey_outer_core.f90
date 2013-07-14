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

  subroutine compute_stacey_outer_core_forward(ichunk,SAVE_FORWARD, &
                              it,ibool_outer_core, &
                              veloc_outer_core,accel_outer_core, &
                              vp_outer_core,wgllwgll_xz,wgllwgll_yz,wgllwgll_xy, &
                              jacobian2D_bottom_outer_core, &
                              jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core, &
                              jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core, &
                              ibelm_bottom_outer_core, &
                              ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
                              ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
                              nimin_outer_core,nimax_outer_core, &
                              njmin_outer_core,njmax_outer_core, &
                              nkmin_xi_outer_core,nkmin_eta_outer_core, &
                              NSPEC2D_BOTTOM, &
                              nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                              nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                              reclen_zmin, &
                              reclen_xmin_outer_core,reclen_xmax_outer_core, &
                              reclen_ymin_outer_core,reclen_ymax_outer_core, &
                              nabs_zmin_oc, &
                              nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc, &
                              absorb_zmin_outer_core, &
                              absorb_xmin_outer_core,absorb_xmax_outer_core, &
                              absorb_ymin_outer_core,absorb_ymax_outer_core)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer ichunk
  integer it
  logical SAVE_FORWARD

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE) :: &
    veloc_outer_core,accel_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE_STACEY) :: vp_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY) :: wgllwgll_xy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ) :: wgllwgll_yz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,NSPEC2D_BOTTOM_OC) :: &
    jacobian2D_bottom_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,NSPEC2DMAX_XMIN_XMAX_OC) :: &
    jacobian2D_xmin_outer_core,jacobian2D_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,NSPEC2DMAX_YMIN_YMAX_OC) :: &
    jacobian2D_ymin_outer_core,jacobian2D_ymax_outer_core

  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_OC) :: ibelm_xmin_outer_core,ibelm_xmax_outer_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_OC) :: ibelm_ymin_outer_core,ibelm_ymax_outer_core

  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_OC) :: &
    nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_OC) :: &
    njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_BOTTOM
  integer nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core

  integer reclen_zmin,reclen_xmin_outer_core,reclen_xmax_outer_core,&
    reclen_ymin_outer_core,reclen_ymax_outer_core

  integer nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc,nabs_zmin_oc
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nabs_xmin_oc) :: absorb_xmin_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nabs_xmax_oc) :: absorb_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nabs_ymin_oc) :: absorb_ymin_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nabs_ymax_oc) :: absorb_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,nabs_zmin_oc) :: absorb_zmin_outer_core

  ! local parameters
  real(kind=CUSTOM_REAL) :: sn,weight
  !integer :: reclen1,reclen2
  integer :: i,j,k,ispec2D,ispec,iglob

  ! note: we use c functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storages (every 8 step) and/or asynchronious
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  ! outer core
  !   xmin
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then

    do ispec2D=1,nspec2D_xmin_outer_core

      ispec=ibelm_xmin_outer_core(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_xi_outer_core(1,ispec2D) == 0 .or. njmin_outer_core(1,ispec2D) == 0) cycle

      i=1
      do k=nkmin_xi_outer_core(1,ispec2D),NGLLZ
        do j=njmin_outer_core(1,ispec2D),njmax_outer_core(1,ispec2D)
          iglob=ibool_outer_core(i,j,k,ispec)

          sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

          weight = jacobian2D_xmin_outer_core(j,k,ispec2D)*wgllwgll_yz(j,k)

          accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

          if (SAVE_FORWARD) then
            absorb_xmin_outer_core(j,k,ispec2D) = weight*sn
          endif
        enddo
      enddo
    enddo

    ! writes absorbing boundary values
    if (SAVE_FORWARD .and. nspec2D_xmin_outer_core > 0 ) then
      call write_abs(4,absorb_xmin_outer_core,reclen_xmin_outer_core,it)
    endif

  endif

  !   xmax
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then

    do ispec2D=1,nspec2D_xmax_outer_core

      ispec=ibelm_xmax_outer_core(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_xi_outer_core(2,ispec2D) == 0 .or. njmin_outer_core(2,ispec2D) == 0) cycle

      i=NGLLX
      do k=nkmin_xi_outer_core(2,ispec2D),NGLLZ
        do j=njmin_outer_core(2,ispec2D),njmax_outer_core(2,ispec2D)
          iglob=ibool_outer_core(i,j,k,ispec)

          sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

          weight = jacobian2D_xmax_outer_core(j,k,ispec2D)*wgllwgll_yz(j,k)

          accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

          if (SAVE_FORWARD) then
            absorb_xmax_outer_core(j,k,ispec2D) = weight*sn
          endif

        enddo
      enddo
    enddo

    if (SAVE_FORWARD .and. nspec2D_xmax_outer_core > 0 ) then
      call write_abs(5,absorb_xmax_outer_core,reclen_xmax_outer_core,it)
    endif

  endif

  !   ymin

  do ispec2D=1,nspec2D_ymin_outer_core

      ispec=ibelm_ymin_outer_core(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_eta_outer_core(1,ispec2D) == 0 .or. nimin_outer_core(1,ispec2D) == 0) cycle

      j=1
      do k=nkmin_eta_outer_core(1,ispec2D),NGLLZ
        do i=nimin_outer_core(1,ispec2D),nimax_outer_core(1,ispec2D)
          iglob=ibool_outer_core(i,j,k,ispec)

          sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

          weight=jacobian2D_ymin_outer_core(i,k,ispec2D)*wgllwgll_xz(i,k)

          accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

        if (SAVE_FORWARD) then
          absorb_ymin_outer_core(i,k,ispec2D) = weight*sn
        endif

      enddo
    enddo
  enddo

  if (SAVE_FORWARD .and. nspec2D_ymin_outer_core > 0 ) then
    call write_abs(6,absorb_ymin_outer_core,reclen_ymin_outer_core,it)
  endif

  !   ymax
  do ispec2D=1,nspec2D_ymax_outer_core

    ispec=ibelm_ymax_outer_core(ispec2D)

    ! exclude elements that are not on absorbing edges
    if(nkmin_eta_outer_core(2,ispec2D) == 0 .or. nimin_outer_core(2,ispec2D) == 0) cycle

    j=NGLLY
    do k=nkmin_eta_outer_core(2,ispec2D),NGLLZ
      do i=nimin_outer_core(2,ispec2D),nimax_outer_core(2,ispec2D)
        iglob=ibool_outer_core(i,j,k,ispec)

        sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

        weight=jacobian2D_ymax_outer_core(i,k,ispec2D)*wgllwgll_xz(i,k)

        accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

        if (SAVE_FORWARD) then
          absorb_ymax_outer_core(i,k,ispec2D) = weight*sn
        endif

      enddo
    enddo
  enddo

  if (SAVE_FORWARD .and. nspec2D_ymax_outer_core > 0 ) then
    call write_abs(7,absorb_ymax_outer_core,reclen_ymax_outer_core,it)
  endif

  ! for surface elements exactly on the ICB
  do ispec2D = 1,NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

    ispec = ibelm_bottom_outer_core(ispec2D)

    k = 1
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool_outer_core(i,j,k,ispec)

        sn = veloc_outer_core(iglob)/vp_outer_core(i,j,k,ispec)

        weight = jacobian2D_bottom_outer_core(i,j,ispec2D)*wgllwgll_xy(i,j)

        accel_outer_core(iglob) = accel_outer_core(iglob) - weight*sn

        if (SAVE_FORWARD) then
          absorb_zmin_outer_core(i,j,ispec2D) = weight*sn
        endif

      enddo
    enddo
  enddo

  if (SAVE_FORWARD .and. NSPEC2D_BOTTOM(IREGION_OUTER_CORE) > 0 ) then
    call write_abs(8,absorb_zmin_outer_core,reclen_zmin,it)
  endif

  end subroutine compute_stacey_outer_core_forward

!=====================================================================
!=====================================================================

  subroutine compute_stacey_outer_core_backward(ichunk, &
                              NSTEP,it,ibool_outer_core, &
                              b_accel_outer_core, &
                              ibelm_bottom_outer_core, &
                              ibelm_xmin_outer_core,ibelm_xmax_outer_core, &
                              ibelm_ymin_outer_core,ibelm_ymax_outer_core, &
                              nimin_outer_core,nimax_outer_core, &
                              njmin_outer_core,njmax_outer_core, &
                              nkmin_xi_outer_core,nkmin_eta_outer_core, &
                              NSPEC2D_BOTTOM, &
                              nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
                              nspec2D_ymin_outer_core,nspec2D_ymax_outer_core, &
                              reclen_zmin, &
                              reclen_xmin_outer_core,reclen_xmax_outer_core, &
                              reclen_ymin_outer_core,reclen_ymax_outer_core, &
                              nabs_zmin_oc, &
                              nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc, &
                              absorb_zmin_outer_core, &
                              absorb_xmin_outer_core,absorb_xmax_outer_core, &
                              absorb_ymin_outer_core,absorb_ymax_outer_core)

  implicit none

  include "constants.h"
  include "OUTPUT_FILES/values_from_mesher.h"

  integer ichunk
  integer NSTEP,it

  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE) :: ibool_outer_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_OUTER_CORE_ADJOINT) :: &
    b_accel_outer_core

  integer, dimension(NSPEC2D_BOTTOM_OC) :: ibelm_bottom_outer_core
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_OC) :: ibelm_xmin_outer_core,ibelm_xmax_outer_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_OC) :: ibelm_ymin_outer_core,ibelm_ymax_outer_core

  integer, dimension(2,NSPEC2DMAX_YMIN_YMAX_OC) :: &
    nimin_outer_core,nimax_outer_core,nkmin_eta_outer_core
  integer, dimension(2,NSPEC2DMAX_XMIN_XMAX_OC) :: &
    njmin_outer_core,njmax_outer_core,nkmin_xi_outer_core

  integer, dimension(MAX_NUM_REGIONS) :: NSPEC2D_BOTTOM
  integer nspec2D_xmin_outer_core,nspec2D_xmax_outer_core, &
    nspec2D_ymin_outer_core,nspec2D_ymax_outer_core

  integer reclen_zmin,reclen_xmin_outer_core,reclen_xmax_outer_core,&
    reclen_ymin_outer_core,reclen_ymax_outer_core

  integer nabs_xmin_oc,nabs_xmax_oc,nabs_ymin_oc,nabs_ymax_oc,nabs_zmin_oc
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nabs_xmin_oc) :: absorb_xmin_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLY,NGLLZ,nabs_xmax_oc) :: absorb_xmax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nabs_ymin_oc) :: absorb_ymin_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nabs_ymax_oc) :: absorb_ymax_outer_core
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLY,nabs_zmin_oc) :: absorb_zmin_outer_core

  ! local parameters
  !integer :: reclen1,reclen2
  integer :: i,j,k,ispec2D,ispec,iglob

  ! note: we use c functions for I/O as they still have a better performance than
  !           Fortran, unformatted file I/O. however, using -assume byterecl together with Fortran functions
  !           comes very close (only  ~ 4 % slower ).
  !
  !           tests with intermediate storages (every 8 step) and/or asynchronious
  !           file access (by process rank modulo 8) showed that the following,
  !           simple approach is still fastest. (assuming that files are accessed on a local scratch disk)

  !   xmin
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AC) then

    ! reads absorbing boundary values
    if (nspec2D_xmin_outer_core > 0)  then
      ! note: backward/reconstructed wavefields are read in after the Newmark time scheme in the first time loop
      !          this leads to a corresponding boundary condition at time index NSTEP - (it-1) = NSTEP - it + 1
      call read_abs(4,absorb_xmin_outer_core,reclen_xmin_outer_core,NSTEP-it+1)
    endif

    do ispec2D=1,nspec2D_xmin_outer_core

      ispec=ibelm_xmin_outer_core(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_xi_outer_core(1,ispec2D) == 0 .or. njmin_outer_core(1,ispec2D) == 0) cycle

      i=1
      do k=nkmin_xi_outer_core(1,ispec2D),NGLLZ
        do j=njmin_outer_core(1,ispec2D),njmax_outer_core(1,ispec2D)
          iglob=ibool_outer_core(i,j,k,ispec)
          b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_xmin_outer_core(j,k,ispec2D)
        enddo
      enddo
    enddo

  endif

  !   xmax
  ! if two chunks exclude this face for one of them
  if(NCHUNKS_VAL == 1 .or. ichunk == CHUNK_AB) then

    if (nspec2D_xmax_outer_core > 0)  then
      call read_abs(5,absorb_xmax_outer_core,reclen_xmax_outer_core,NSTEP-it+1)
    endif

    do ispec2D=1,nspec2D_xmax_outer_core

      ispec=ibelm_xmax_outer_core(ispec2D)

      ! exclude elements that are not on absorbing edges
      if(nkmin_xi_outer_core(2,ispec2D) == 0 .or. njmin_outer_core(2,ispec2D) == 0) cycle

      i=NGLLX
      do k=nkmin_xi_outer_core(2,ispec2D),NGLLZ
        do j=njmin_outer_core(2,ispec2D),njmax_outer_core(2,ispec2D)
          iglob=ibool_outer_core(i,j,k,ispec)
          b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_xmax_outer_core(j,k,ispec2D)
        enddo
      enddo
    enddo

  endif

  !   ymin
  if (nspec2D_ymin_outer_core > 0)  then
    call read_abs(6,absorb_ymin_outer_core,reclen_ymin_outer_core,NSTEP-it+1)
  endif

  do ispec2D=1,nspec2D_ymin_outer_core

    ispec=ibelm_ymin_outer_core(ispec2D)

    ! exclude elements that are not on absorbing edges
    if(nkmin_eta_outer_core(1,ispec2D) == 0 .or. nimin_outer_core(1,ispec2D) == 0) cycle

    j=1
    do k=nkmin_eta_outer_core(1,ispec2D),NGLLZ
      do i=nimin_outer_core(1,ispec2D),nimax_outer_core(1,ispec2D)
        iglob=ibool_outer_core(i,j,k,ispec)
        b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_ymin_outer_core(i,k,ispec2D)
      enddo
    enddo
  enddo

  !   ymax
  if (nspec2D_ymax_outer_core > 0)  then
    call read_abs(7,absorb_ymax_outer_core,reclen_ymax_outer_core,NSTEP-it+1)
  endif

  do ispec2D=1,nspec2D_ymax_outer_core

    ispec=ibelm_ymax_outer_core(ispec2D)

    ! exclude elements that are not on absorbing edges
    if(nkmin_eta_outer_core(2,ispec2D) == 0 .or. nimin_outer_core(2,ispec2D) == 0) cycle

    j=NGLLY
    do k=nkmin_eta_outer_core(2,ispec2D),NGLLZ
      do i=nimin_outer_core(2,ispec2D),nimax_outer_core(2,ispec2D)
        iglob=ibool_outer_core(i,j,k,ispec)
        b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_ymax_outer_core(i,k,ispec2D)
      enddo
    enddo
  enddo

  ! for surface elements exactly on the ICB
  if (NSPEC2D_BOTTOM(IREGION_OUTER_CORE)> 0)  then
    call read_abs(8,absorb_zmin_outer_core,reclen_zmin,NSTEP-it+1)
  endif

  do ispec2D = 1,NSPEC2D_BOTTOM(IREGION_OUTER_CORE)

    ispec = ibelm_bottom_outer_core(ispec2D)

    k = 1
    do j = 1,NGLLY
      do i = 1,NGLLX
        iglob = ibool_outer_core(i,j,k,ispec)
        b_accel_outer_core(iglob) = b_accel_outer_core(iglob) - absorb_zmin_outer_core(i,j,ispec2D)
      enddo
    enddo
  enddo

  end subroutine compute_stacey_outer_core_backward

