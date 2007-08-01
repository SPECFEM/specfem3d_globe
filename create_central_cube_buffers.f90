!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!                    and University of Pau, France
! (c) California Institute of Technology and University of Pau, April 2007
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

!
!--- create buffers to assemble with central cube
!

  subroutine create_central_cube_buffers(myrank,iproc_xi,iproc_eta,ichunk, &
       NPROC_XI,NPROC_ETA,NCHUNKS,NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
       NSPEC2DMAX_XMIN_XMAX_INNER_CORE,NSPEC2DMAX_YMIN_YMAX_INNER_CORE,NSPEC2D_BOTTOM_INNER_CORE, &
       addressing,ibool_inner_core,idoubling_inner_core, &
       xstore_inner_core,ystore_inner_core,zstore_inner_core, &
       nspec2D_xmin_inner_core,nspec2D_xmax_inner_core,nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
       ibelm_xmin_inner_core,ibelm_xmax_inner_core,ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
       nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices, &
       receiver_cube_from_slices,sender_from_slices_to_cube,ibool_central_cube,buffer_slices,buffer_all_cube_from_slices)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

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
  double precision, dimension(npoin2D_cube_from_slices,NDIM), intent(out) :: buffer_slices
  double precision, dimension(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM), intent(out) :: &
        buffer_all_cube_from_slices

! local variables below
  integer i,j,k,ispec,ispec2D,iglob,ier
  integer sender,receiver,imsg,ipoin,iproc_xi_loop

  double precision x_target,y_target,z_target
  double precision x_current,y_current,z_current

! MPI status of messages to be received
  integer msg_status(MPI_STATUS_SIZE)

!--- processor to send information to in cube from slices

! four vertical sides first
  if(ichunk == CHUNK_AC) then
    receiver_cube_from_slices = addressing(CHUNK_AB,0,iproc_eta)

  else if(ichunk == CHUNK_BC) then
    receiver_cube_from_slices = addressing(CHUNK_AB,iproc_eta,NPROC_ETA-1)

  else if(ichunk == CHUNK_AC_ANTIPODE) then
    receiver_cube_from_slices = addressing(CHUNK_AB,NPROC_XI-1,iproc_eta)

  else if(ichunk == CHUNK_BC_ANTIPODE) then
    receiver_cube_from_slices = addressing(CHUNK_AB,NPROC_XI-1-iproc_eta,0)

! bottom of cube, direct correspondence but with inverted xi axis
  else if(ichunk == CHUNK_AB_ANTIPODE) then
    receiver_cube_from_slices = addressing(CHUNK_AB,NPROC_XI-1-iproc_xi,iproc_eta)

  else ! case of CHUNK_AB that carries the cube, therefore use a dummy value
    receiver_cube_from_slices = -1
  endif


!--- list of processors to receive information from in cube

! only for slices in central cube
  if(ichunk == CHUNK_AB) then

! initialize index of sender
    imsg = 0

! define sender for bottom edge
! bottom of cube, direct correspondence but with inverted xi axis
    imsg = imsg + 1
    sender_from_slices_to_cube(imsg) = addressing(CHUNK_AB_ANTIPODE,NPROC_XI-1-iproc_xi,iproc_eta)

! define sender for xi = xi_min edge
    if(iproc_xi == 0) then
      do iproc_xi_loop = 0,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC,iproc_xi_loop,iproc_eta)
      enddo
    endif

! define sender for xi = xi_max edge
    if(iproc_xi == NPROC_XI-1) then
      do iproc_xi_loop = 0,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_AC_ANTIPODE,iproc_xi_loop,iproc_eta)
      enddo
    endif

! define sender for eta = eta_min edge
    if(iproc_eta == 0) then
      do iproc_xi_loop = 0,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC_ANTIPODE,iproc_xi_loop,NPROC_ETA-1-iproc_xi)
      enddo
    endif

! define sender for eta = eta_max edge
    if(iproc_eta == NPROC_ETA-1) then
      do iproc_xi_loop = 0,NPROC_XI-1
        imsg = imsg + 1
        sender_from_slices_to_cube(imsg) = addressing(CHUNK_BC,iproc_xi_loop,iproc_xi)
      enddo
    endif

! check that total number of faces found is correct
   if(imsg /= nb_msgs_theor_in_cube) call exit_MPI(myrank,'wrong number of faces found for central cube')

  else

! dummy value in slices
    sender_from_slices_to_cube(1) = -1

  endif

! on chunk AB, receive all the messages from slices
  if(ichunk == CHUNK_AB) then

   do imsg = 1,nb_msgs_theor_in_cube

! receive buffers from slices
  sender = sender_from_slices_to_cube(imsg)
  call MPI_RECV(buffer_slices, &
              NDIM*npoin2D_cube_from_slices,MPI_DOUBLE_PRECISION,sender, &
              itag,MPI_COMM_WORLD,msg_status,ier)

! copy buffer in 2D array for each slice
   buffer_all_cube_from_slices(imsg,:,:) = buffer_slices(:,:)

   enddo
   endif


! send info to central cube from all the slices except those in CHUNK_AB
  if(ichunk /= CHUNK_AB) then

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
    call MPI_SEND(buffer_slices,NDIM*npoin2D_cube_from_slices, &
              MPI_DOUBLE_PRECISION,receiver,itag,MPI_COMM_WORLD,ier)

 endif  ! end sending info to central cube

!--- now we need to find the points received and create indirect addressing

  if(ichunk == CHUNK_AB) then

   do imsg = 1,nb_msgs_theor_in_cube

   do ipoin = 1,npoin2D_cube_from_slices

     x_target = buffer_all_cube_from_slices(imsg,ipoin,1)
     y_target = buffer_all_cube_from_slices(imsg,ipoin,2)
     z_target = buffer_all_cube_from_slices(imsg,ipoin,3)

! x = x_min
  do ispec2D = 1,nspec2D_xmin_inner_core

      ispec = ibelm_xmin_inner_core(ispec2D)

! do not loop on elements outside of the central cube
     if(idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle

     i = 1
     do k = 1,NGLLZ
       do j = 1,NGLLY

         iglob = ibool_inner_core(i,j,k,ispec)
         x_current = dble(xstore_inner_core(iglob))
         y_current = dble(ystore_inner_core(iglob))
         z_current = dble(zstore_inner_core(iglob))

! look for matching point
         if(dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
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
     if(idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle

     i = NGLLX
     do k = 1,NGLLZ
       do j = 1,NGLLY

         iglob = ibool_inner_core(i,j,k,ispec)
         x_current = dble(xstore_inner_core(iglob))
         y_current = dble(ystore_inner_core(iglob))
         z_current = dble(zstore_inner_core(iglob))

! look for matching point
         if(dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
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
     if(idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle

     j = 1
     do k = 1,NGLLZ
       do i = 1,NGLLX

         iglob = ibool_inner_core(i,j,k,ispec)
         x_current = dble(xstore_inner_core(iglob))
         y_current = dble(ystore_inner_core(iglob))
         z_current = dble(zstore_inner_core(iglob))

! look for matching point
         if(dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
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
     if(idoubling_inner_core(ispec) /= IFLAG_MIDDLE_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE .and. &
        idoubling_inner_core(ispec) /= IFLAG_TOP_CENTRAL_CUBE) cycle

     j = NGLLY
     do k = 1,NGLLZ
       do i = 1,NGLLX

         iglob = ibool_inner_core(i,j,k,ispec)
         x_current = dble(xstore_inner_core(iglob))
         y_current = dble(ystore_inner_core(iglob))
         z_current = dble(zstore_inner_core(iglob))

! look for matching point
         if(dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
           ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
           goto 100
         endif

       enddo
     enddo

   enddo

! bottom of cube
  do ispec = 1,NSPEC_INNER_CORE

! loop on elements at the bottom of the cube only
     if(idoubling_inner_core(ispec) /= IFLAG_BOTTOM_CENTRAL_CUBE) cycle

     k = 1
     do j = 1,NGLLY
       do i = 1,NGLLX

         iglob = ibool_inner_core(i,j,k,ispec)
         x_current = dble(xstore_inner_core(iglob))
         y_current = dble(ystore_inner_core(iglob))
         z_current = dble(zstore_inner_core(iglob))

! look for matching point
         if(dsqrt((x_current-x_target)**2 + (y_current-y_target)**2 + (z_current-z_target)**2) < SMALLVALTOL) then
           ibool_central_cube(imsg,ipoin) = ibool_inner_core(i,j,k,ispec)
           goto 100
         endif

       enddo
     enddo

   enddo

! check that a matching point is found in all cases
  call exit_MPI(myrank,'point never found in central cube')

 100 continue

   enddo
   enddo
   endif

  end subroutine create_central_cube_buffers

!
!----------------------------------
!

  subroutine comp_central_cube_buffer_size(iproc_xi,iproc_eta,ichunk,NPROC_XI,NPROC_ETA,NSPEC2D_BOTTOM_INNER_CORE, &
                nb_msgs_theor_in_cube,npoin2D_cube_from_slices)

!--- compute number of messages to expect in cube as well as their size
!--- take into account vertical sides and bottom side

  implicit none

  include "constants.h"

  integer, intent(in) :: iproc_xi,iproc_eta,ichunk,NPROC_XI,NPROC_ETA,NSPEC2D_BOTTOM_INNER_CORE

  integer, intent(out) :: nb_msgs_theor_in_cube,npoin2D_cube_from_slices

! only for slices in central cube
  if(ichunk == CHUNK_AB) then
    if(NPROC_XI == 1) then
! five sides if only one processor in cube
      nb_msgs_theor_in_cube = 5
    else
! case of a corner
      if((iproc_xi == 0 .or. iproc_xi == NPROC_XI-1).and. &
         (iproc_eta == 0 .or. iproc_eta == NPROC_ETA-1)) then
! slices on both "vertical" faces plus one slice at the bottom
        nb_msgs_theor_in_cube = 2*NPROC_XI + 1
! case of an edge
      else if(iproc_xi == 0 .or. iproc_xi == NPROC_XI-1 .or. &
              iproc_eta == 0 .or. iproc_eta == NPROC_ETA-1) then
! slices on the "vertical" face plus one slice at the bottom
        nb_msgs_theor_in_cube = NPROC_XI + 1
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

