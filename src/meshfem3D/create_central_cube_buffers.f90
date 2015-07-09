!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  7 . 0
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

!
!--- create buffers to assemble with central cube
!

  subroutine create_central_cube_buffers(myrank,iproc_xi,iproc_eta,ichunk, &
                   NPROC_XI,NPROC_ETA,NCHUNKS, &
                   NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                   NSPEC2DMAX_XMIN_XMAX_INNER_CORE,NSPEC2DMAX_YMIN_YMAX_INNER_CORE, &
                   NSPEC2D_BOTTOM_INNER_CORE, &
                   addressing,ibool_inner_core,idoubling_inner_core, &
                   xstore_inner_core,ystore_inner_core,zstore_inner_core, &
                   nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
                   nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
                   ibelm_xmin_inner_core,ibelm_xmax_inner_core, &
                   ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
                   nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices, &
                   receiver_cube_from_slices,sender_from_slices_to_cube,ibool_central_cube, &
                   buffer_slices,buffer_slices2,buffer_all_cube_from_slices)

  use constants

  implicit none

  integer, intent(in) :: myrank,iproc_xi,iproc_eta,ichunk, &
       NPROC_XI,NPROC_ETA,NCHUNKS,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
       NSPEC2DMAX_XMIN_XMAX_INNER_CORE,NSPEC2DMAX_YMIN_YMAX_INNER_CORE,NSPEC2D_BOTTOM_INNER_CORE

  ! for addressing of the slices
  integer, dimension(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1), intent(in) :: addressing

  ! mesh parameters
  integer, dimension(NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), intent(in) :: ibool_inner_core

  ! local to global mapping
  integer, dimension(NSPEC_INNER_CORE), intent(in) :: idoubling_inner_core

  real(kind=CUSTOM_REAL), dimension(NGLOB_INNER_CORE), intent(in) :: xstore_inner_core,ystore_inner_core,zstore_inner_core

  ! boundary parameters locator
  integer, intent(in) :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core,nspec2D_ymin_inner_core,nspec2D_ymax_inner_core
  integer, dimension(NSPEC2DMAX_XMIN_XMAX_INNER_CORE), intent(in) :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(NSPEC2DMAX_YMIN_YMAX_INNER_CORE), intent(in) :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(NSPEC2D_BOTTOM_INNER_CORE), intent(in) :: ibelm_bottom_inner_core

  integer, intent(in) :: nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices

  ! for matching with central cube in inner core
  integer, intent(out) :: receiver_cube_from_slices

  integer, dimension(non_zero_nb_msgs_theor_in_cube), intent(out) :: sender_from_slices_to_cube
  integer, dimension(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices), intent(out) :: ibool_central_cube
  double precision, dimension(npoin2D_cube_from_slices,NDIM), intent(out) :: buffer_slices,buffer_slices2
  double precision, dimension(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM), intent(out) :: &
        buffer_all_cube_from_slices

  ! local variables below
  integer i,j,k,ispec,ispec2D,iglob
  integer sender,receiver,imsg,ipoin,iproc_xi_loop

  double precision x_target,y_target,z_target
  double precision x_current,y_current,z_current

  integer :: nproc_xi_half_floor,nproc_xi_half_ceil

  if (mod(NPROC_XI,2) /= 0) then
    nproc_xi_half_floor = floor(NPROC_XI/2.d0)
    nproc_xi_half_ceil = ceiling(NPROC_XI/2.d0)
  else
    nproc_xi_half_floor = NPROC_XI/2
    nproc_xi_half_ceil = NPROC_XI/2
  endif

  ! check that the number of points in this slice is correct
  if (minval(ibool_inner_core(:,:,:,:)) /= 1 .or. maxval(ibool_inner_core(:,:,:,:)) /= NGLOB_INNER_CORE) &
    call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in inner core')


!--- processor to send information to in cube from slices

! four vertical sides first
  if (ichunk == CHUNK_AC) then
    if (iproc_xi < nproc_xi_half_floor) then
      receiver_cube_from_slices = addressing(CHUNK_AB_ANTIPODE,NPROC_XI-1,iproc_eta)
    else
      receiver_cube_from_slices = addressing(CHUNK_AB,0,iproc_eta)
    endif
  else if (ichunk == CHUNK_BC) then
    if (iproc_xi < nproc_xi_half_floor) then
      receiver_cube_from_slices = addressing(CHUNK_AB_ANTIPODE,NPROC_XI-1-iproc_eta,NPROC_ETA-1)
    else
      receiver_cube_from_slices = addressing(CHUNK_AB,iproc_eta,NPROC_ETA-1)
    endif
  else if (ichunk == CHUNK_AC_ANTIPODE) then
    if (iproc_xi <= ceiling((NPROC_XI/2.d0)-1)) then
      receiver_cube_from_slices = addressing(CHUNK_AB,NPROC_XI-1,iproc_eta)
    else
      receiver_cube_from_slices = addressing(CHUNK_AB_ANTIPODE,0,iproc_eta)
    endif
  else if (ichunk == CHUNK_BC_ANTIPODE) then
    if (iproc_xi < nproc_xi_half_floor) then
      receiver_cube_from_slices = addressing(CHUNK_AB_ANTIPODE,iproc_eta,0)
    else
      receiver_cube_from_slices = addressing(CHUNK_AB,NPROC_XI-1-iproc_eta,0)
    endif
! bottom of cube, direct correspondence but with inverted xi axis
  else if (ichunk == CHUNK_AB_ANTIPODE) then
    receiver_cube_from_slices = addressing(CHUNK_AB,NPROC_XI-1-iproc_xi,iproc_eta)
  else if (ichunk == CHUNK_AB) then
    receiver_cube_from_slices = addressing(CHUNK_AB_ANTIPODE,NPROC_XI-1-iproc_xi,iproc_eta)
  endif


!--- list of processors to receive information from in cube

! only for slices in central cube
  if (ichunk == CHUNK_AB) then
    ! initialize index of sender
    imsg = 0

    ! define sender for xi = xi_min edge
    if (iproc_xi == 0) then
      do iproc_xi_loop = nproc_xi_half_floor,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC,iproc_xi_loop,iproc_eta)
      enddo
    endif

    ! define sender for xi = xi_max edge
    if (iproc_xi == NPROC_XI-1) then
      do iproc_xi_loop = 0, floor((NPROC_XI-1)/2.d0)
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC_ANTIPODE,iproc_xi_loop,iproc_eta)
      enddo
    endif

    ! define sender for eta = eta_min edge
    if (iproc_eta == 0) then
      do iproc_xi_loop = nproc_xi_half_floor,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC_ANTIPODE,iproc_xi_loop,NPROC_ETA-1-iproc_xi)
      enddo
    endif

    ! define sender for eta = eta_max edge
    if (iproc_eta == NPROC_ETA-1) then
      do iproc_xi_loop = nproc_xi_half_floor,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC,iproc_xi_loop,iproc_xi)
      enddo
    endif

    ! define sender for bottom edge
    ! bottom of cube, direct correspondence but with inverted xi axis
    imsg = imsg + 1
    sender_from_slices_to_cube(imsg) = addressing(CHUNK_AB_ANTIPODE,NPROC_XI-1-iproc_xi,iproc_eta)

    ! check that total number of faces found is correct
    if (imsg /= nb_msgs_theor_in_cube) then
      print*,'Error ',myrank,'nb_msgs_theor_in_cube:',nb_msgs_theor_in_cube,imsg
      call exit_MPI(myrank,'wrong number of faces found for central cube')
    endif

  else if (ichunk == CHUNK_AB_ANTIPODE) then
    ! initialize index of sender
    imsg = 0

    ! define sender for xi = xi_min edge
    if (iproc_xi == 0) then
      do iproc_xi_loop = nproc_xi_half_ceil,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC_ANTIPODE,iproc_xi_loop,iproc_eta)
      enddo
    endif

    ! define sender for xi = xi_max edge
    if (iproc_xi == NPROC_XI-1) then
      do iproc_xi_loop = 0, floor((NPROC_XI/2.d0)-1.d0)
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC,iproc_xi_loop,iproc_eta)
      enddo
    endif

    ! define sender for eta = eta_min edge
    if (iproc_eta == 0) then
      do iproc_xi_loop = 0, floor((NPROC_XI/2.d0)-1.d0)
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC_ANTIPODE,iproc_xi_loop,iproc_xi)
      enddo
    endif

    ! define sender for eta = eta_max edge
    if (iproc_eta == NPROC_ETA-1) then
      do iproc_xi_loop = 0, floor((NPROC_XI/2.d0)-1.d0)
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC,iproc_xi_loop,NPROC_ETA-1-iproc_xi)
      enddo
    endif

    ! in case NPROC_XI == 1, the other chunks exchange all bottom points with
    ! CHUNK_AB **and** CHUNK_AB_ANTIPODE
    if (NPROC_XI == 1) then
      ! define sender for xi = xi_min edge
      if (iproc_xi == 0) then
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC_ANTIPODE,0,iproc_eta)
      endif

      ! define sender for xi = xi_max edge
      if (iproc_xi == NPROC_XI-1) then
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC,0,iproc_eta)
      endif

      ! define sender for eta = eta_min edge
      if (iproc_eta == 0) then
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC_ANTIPODE,0,iproc_xi)
      endif

      ! define sender for eta = eta_max edge
      if (iproc_eta == NPROC_ETA-1) then
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC,0,NPROC_ETA-1-iproc_xi)
      endif
    endif

    ! define sender for bottom edge
    ! bottom of cube, direct correspondence but with inverted xi axis
    imsg = imsg + 1
    sender_from_slices_to_cube(imsg) = addressing(CHUNK_AB,NPROC_XI-1-iproc_xi,iproc_eta)

    ! check that total number of faces found is correct
    if (imsg /= nb_msgs_theor_in_cube) then
      print*,'Error ',myrank,'nb_msgs_theor_in_cube:',nb_msgs_theor_in_cube,imsg
      call exit_MPI(myrank,'wrong number of faces found for central cube')
    endif

  else

    ! dummy value in slices
    sender_from_slices_to_cube(1) = -1

  endif


! on chunk AB & AB ANTIPODE, receive all (except bottom) the messages from slices
  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    do imsg = 1,nb_msgs_theor_in_cube-1

    ! receive buffers from slices
    sender = sender_from_slices_to_cube(imsg)
    call recv_dp(buffer_slices,NDIM*npoin2D_cube_from_slices,sender,itag)

    ! copy buffer in 2D array for each slice
    buffer_all_cube_from_slices(imsg,:,:) = buffer_slices(:,:)

    enddo
  endif

  ! send info to central cube from all the slices except those in CHUNK_AB & CHUNK_AB_ANTIPODE
  if (ichunk /= CHUNK_AB .and. ichunk /= CHUNK_AB_ANTIPODE) then
    ! for bottom elements in contact with central cube from the slices side
    ipoin = 0
    do ispec2D = 1,NSPEC2D_BOTTOM_INNER_CORE

      ispec = ibelm_bottom_inner_core(ispec2D)

      ! only for DOFs exactly on surface of central cube (bottom of these elements)
      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          ipoin = ipoin + 1
          iglob = ibool_inner_core(i,j,k,ispec)
          buffer_slices(ipoin,1) = dble(xstore_inner_core(iglob))
          buffer_slices(ipoin,2) = dble(ystore_inner_core(iglob))
          buffer_slices(ipoin,3) = dble(zstore_inner_core(iglob))
        enddo
      enddo
    enddo

    ! send buffer to central cube
    receiver = receiver_cube_from_slices
    call send_dp(buffer_slices,NDIM*npoin2D_cube_from_slices,receiver,itag)

    ! in case NPROC_XI == 1, the other chunks exchange all bottom points with
    ! CHUNK_AB **and** CHUNK_AB_ANTIPODE
    if (NPROC_XI == 1) then
      call send_dp(buffer_slices,NDIM*npoin2D_cube_from_slices,addressing(CHUNK_AB_ANTIPODE,0,iproc_eta),itag)
    endif

  endif  ! end sending info to central cube


  ! exchange of their bottom faces between chunks AB and AB_ANTIPODE
  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then
    ipoin = 0
    do ispec = NSPEC_INNER_CORE, 1, -1
      if (idoubling_inner_core(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE) then
        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            ipoin = ipoin + 1
            iglob = ibool_inner_core(i,j,k,ispec)
            buffer_slices(ipoin,1) = dble(xstore_inner_core(iglob))
            buffer_slices(ipoin,2) = dble(ystore_inner_core(iglob))
            buffer_slices(ipoin,3) = dble(zstore_inner_core(iglob))
          enddo
        enddo
      endif
    enddo
    if (ipoin /= npoin2D_cube_from_slices) then
      print*,'Error',myrank,'bottom points:',npoin2D_cube_from_slices,ipoin
      call exit_MPI(myrank,'wrong number of points found for bottom CC AB or !AB')
    endif

    sender = sender_from_slices_to_cube(nb_msgs_theor_in_cube)

    call sendrecv_dp(buffer_slices,NDIM*npoin2D_cube_from_slices,receiver_cube_from_slices,itag, &
                     buffer_slices2,NDIM*npoin2D_cube_from_slices,sender,itag)

    buffer_all_cube_from_slices(nb_msgs_theor_in_cube,:,:) = buffer_slices2(:,:)

  endif

  !--- now we need to find the points received and create indirect addressing
  ibool_central_cube(:,:) = -1

  if (ichunk == CHUNK_AB .or. ichunk == CHUNK_AB_ANTIPODE) then

   do imsg = 1,nb_msgs_theor_in_cube

    do ipoin = 1,npoin2D_cube_from_slices

      x_target = buffer_all_cube_from_slices(imsg,ipoin,1)
      y_target = buffer_all_cube_from_slices(imsg,ipoin,2)
      z_target = buffer_all_cube_from_slices(imsg,ipoin,3)

      ! x = x_min
      do ispec2D = 1,nspec2D_xmin_inner_core
        ispec = ibelm_xmin_inner_core(ispec2D)
        ! do not loop on elements outside of the central cube
        if (idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
          idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
          idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle
        ! check
        if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE ) stop 'Error xmin ibelm'
        i = 1
        do k = 1,NGLLZ
          do j = 1,NGLLY
           iglob = ibool_inner_core(i,j,k,ispec)
           x_current = dble(xstore_inner_core(iglob))
           y_current = dble(ystore_inner_core(iglob))
           z_current = dble(zstore_inner_core(iglob))
           ! look for matching point
           if (dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
             ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
             goto 100
           endif
          enddo
        enddo
      enddo

      ! x = x_max
      do ispec2D = 1,nspec2D_xmax_inner_core
        ispec = ibelm_xmax_inner_core(ispec2D)
        ! do not loop on elements outside of the central cube
        if (idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
            idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
            idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle
        !check
        if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE ) stop 'Error xmax ibelm'
        i = NGLLX
        do k = 1,NGLLZ
          do j = 1,NGLLY
            iglob = ibool_inner_core(i,j,k,ispec)
            x_current = dble(xstore_inner_core(iglob))
            y_current = dble(ystore_inner_core(iglob))
            z_current = dble(zstore_inner_core(iglob))
            ! look for matching point
            if (dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
              ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
              goto 100
            endif
          enddo
        enddo
      enddo

      ! y = y_min
      do ispec2D = 1,nspec2D_ymin_inner_core
        ispec = ibelm_ymin_inner_core(ispec2D)
        ! do not loop on elements outside of the central cube
        if (idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
            idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
            idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle
        !check
        if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE ) stop 'Error ymin ibelm'
        j = 1
        do k = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            x_current = dble(xstore_inner_core(iglob))
            y_current = dble(ystore_inner_core(iglob))
            z_current = dble(zstore_inner_core(iglob))
            ! look for matching point
            if (dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
              ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
              goto 100
            endif
          enddo
        enddo
      enddo

      ! y = y_max
      do ispec2D = 1,nspec2D_ymax_inner_core
        ispec = ibelm_ymax_inner_core(ispec2D)
        ! do not loop on elements outside of the central cube
        if (idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
            idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
            idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle
        !check
        if (idoubling_inner_core(ispec) == IFLAG_IN_FICTITIOUS_CUBE ) stop 'Error ymax ibelm'
        j = NGLLY
        do k = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            x_current = dble(xstore_inner_core(iglob))
            y_current = dble(ystore_inner_core(iglob))
            z_current = dble(zstore_inner_core(iglob))
            ! look for matching point
            if (dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
              ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
              goto 100
            endif
          enddo
        enddo
      enddo

      ! bottom of cube
      do ispec = 1,NSPEC_INNER_CORE
        ! loop on elements at the bottom of the cube only
        if (idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE) cycle
        k = 1
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool_inner_core(i,j,k,ispec)
            x_current = dble(xstore_inner_core(iglob))
            y_current = dble(ystore_inner_core(iglob))
            z_current = dble(zstore_inner_core(iglob))
            ! look for matching point
            if (dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
              ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
              goto 100
            endif
          enddo
        enddo
      enddo

      ! point not found so far
      if (NPROC_XI == 1) then
        ! ignores point
        ibool_central_cube(imsg,ipoin) = 0
      else
        ! check that a matching point is found in all cases
        call exit_MPI(myrank,'point never found in central cube')
      endif

 100  continue

    enddo ! ipoin

    ! checks ibool array
    if (NPROC_XI == 1) then
      if (minval(ibool_central_cube(imsg,:)) < 0 ) call exit_mpi(myrank,'Error ibool_central_cube point not found')

      ! removes points on bottom surface in antipode chunk for other chunks than its AB sharing chunk
      ! (to avoid adding the same point twice from other chunks)
      if (ichunk == CHUNK_AB_ANTIPODE .and. imsg < nb_msgs_theor_in_cube) then
        do ipoin = 1,npoin2D_cube_from_slices
          x_target = buffer_all_cube_from_slices(imsg,ipoin,1)
          y_target = buffer_all_cube_from_slices(imsg,ipoin,2)
          z_target = buffer_all_cube_from_slices(imsg,ipoin,3)

          ! bottom of cube
          do ispec = 1,NSPEC_INNER_CORE
            ! loop on elements at the bottom of the cube only
            if (idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE) cycle
            k = 1
            do j = 1,NGLLY
              do i = 1,NGLLX
                iglob = ibool_inner_core(i,j,k,ispec)
                x_current = dble(xstore_inner_core(iglob))
                y_current = dble(ystore_inner_core(iglob))
                z_current = dble(zstore_inner_core(iglob))
                ! look for matching point
                if (dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
                  ibool_central_cube(imsg,ipoin) = 0
                  goto 200
                endif
              enddo
            enddo
          enddo

 200      continue

        enddo ! ipoin
      endif

    endif ! NPROC_XI == 1

   enddo ! imsg
  endif

  end subroutine create_central_cube_buffers

!
!----------------------------------
!

  subroutine comp_central_cube_buffer_size(iproc_xi,iproc_eta,ichunk,NPROC_XI,NPROC_ETA,NSPEC2D_BOTTOM_INNER_CORE, &
                nb_msgs_theor_in_cube,npoin2D_cube_from_slices)

!--- compute number of messages to expect in cube as well as their size
!--- take into account vertical sides and bottom side

  use constants

  implicit none

  integer, intent(in) :: iproc_xi,iproc_eta,ichunk,NPROC_XI,NPROC_ETA,NSPEC2D_BOTTOM_INNER_CORE

  integer, intent(out) :: nb_msgs_theor_in_cube,npoin2D_cube_from_slices

  integer :: nproc_xi_half_floor,nproc_xi_half_ceil

  if (mod(NPROC_XI,2) /= 0) then
    nproc_xi_half_floor = floor(NPROC_XI/2.d0)
    nproc_xi_half_ceil = ceiling(NPROC_XI/2.d0)
  else
    nproc_xi_half_floor = NPROC_XI/2
    nproc_xi_half_ceil = NPROC_XI/2
  endif

! only for slices in central cube
  if (ichunk == CHUNK_AB) then
    if (NPROC_XI == 1) then
      ! five sides if only one processor in cube
      nb_msgs_theor_in_cube = 5
    else
      ! case of a corner
      if ((iproc_xi == 0 .or. iproc_xi == NPROC_XI-1).and. &
         (iproc_eta == 0 .or. iproc_eta == NPROC_ETA-1)) then
        ! slices on both "vertical" faces plus one slice at the bottom
        nb_msgs_theor_in_cube = 2*(nproc_xi_half_ceil) + 1
      ! case of an edge
      else if (iproc_xi == 0 .or. iproc_xi == NPROC_XI-1 .or. &
              iproc_eta == 0 .or. iproc_eta == NPROC_ETA-1) then
        ! slices on the "vertical" face plus one slice at the bottom
        nb_msgs_theor_in_cube = nproc_xi_half_ceil + 1
      else
        ! bottom element only
        nb_msgs_theor_in_cube = 1
      endif
    endif
  else if (ichunk == CHUNK_AB_ANTIPODE) then
    if (NPROC_XI == 1) then
      ! five sides if only one processor in cube
      nb_msgs_theor_in_cube = 5
    else
      ! case of a corner
      if ((iproc_xi == 0 .or. iproc_xi == NPROC_XI-1).and. &
         (iproc_eta == 0 .or. iproc_eta == NPROC_ETA-1)) then
        ! slices on both "vertical" faces plus one slice at the bottom
        nb_msgs_theor_in_cube = 2*(nproc_xi_half_floor) + 1
      ! case of an edge
      else if (iproc_xi == 0 .or. iproc_xi == NPROC_XI-1 .or. &
              iproc_eta == 0 .or. iproc_eta == NPROC_ETA-1) then
        ! slices on the "vertical" face plus one slice at the bottom
        nb_msgs_theor_in_cube = nproc_xi_half_floor + 1
      else
        ! bottom element only
        nb_msgs_theor_in_cube = 1
      endif
    endif
  else
    ! not in chunk AB
    nb_msgs_theor_in_cube = 0
  endif

  ! number of points to send or receive (bottom of slices)
  npoin2D_cube_from_slices = NSPEC2D_BOTTOM_INNER_CORE * NGLLX * NGLLY

  end subroutine comp_central_cube_buffer_size

