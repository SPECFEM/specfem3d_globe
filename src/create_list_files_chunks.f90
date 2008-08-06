!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  4 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory, California Institute of Technology, USA
!             and University of Pau / CNRS / INRIA, France
! (c) California Institute of Technology and University of Pau / CNRS / INRIA
!                            February 2008
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

! subroutine to create the list of messages to assemble between chunks in files if more than one chunk

!! DK DK for merged version: a lot of useless code / useless lines car probably be suppressed
!! DK DK in this new routine below

  subroutine create_list_files_chunks(iregion_code, &
                nglob_ori,NPROC_XI,NPROC_ETA,NPROCTOT,NGLOB1D_RADIAL_CORNER, &
                myrank,addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice,NCHUNKS)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

!! DK DK for the merged version
! include values created by the mesher
  include "OUTPUT_FILES/values_from_mesher.h"

  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NGLOB1D_RADIAL_CORNER

  integer nglob,nglob_ori
  integer NPROC_XI,NPROC_ETA,NPROCTOT,NGLOB1D_RADIAL_my_corner
  integer myrank,NCHUNKS

  character(len=150) OUTPUT_FILES

! pairs generated theoretically
! four sides for each of the three types of messages
  integer, dimension(:), allocatable :: iproc_sender,iproc_receiver,npoin2D_send,npoin2D_receive

! arrays to assemble the corners (3 processors for each corner)
  integer, dimension(:,:), allocatable :: iprocscorners,itypecorner

  integer ichunk_send,iproc_xi_send,iproc_eta_send
  integer ichunk_receive,iproc_xi_receive,iproc_eta_receive
  integer iproc_loop,iproc_xi_loop,iproc_eta_loop
  integer iproc_xi_loop_inv,iproc_eta_loop_inv
  integer imember_corner

  integer iregion_code

  integer iproc_edge_send,iproc_edge_receive
  integer imsg_type,iside,imode_comm,iedge

  integer ier

! current message number
  integer imsg

! for addressing of the slices
  integer ichunk,iproc_xi,iproc_eta,iproc
  integer addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1)
  integer ichunk_slice(0:NPROCTOT-1)
  integer iproc_xi_slice(0:NPROCTOT-1)

  integer iproc_eta_slice(0:NPROCTOT-1)

! this to avoid problem at compile time if less than six chunks
  integer addressing_big(NCHUNKS_MAX,0:NPROC_XI-1,0:NPROC_ETA-1)

! number of faces between chunks
  integer NUM_FACES,NUMMSGS_FACES

! number of corners between chunks
  integer NCORNERSCHUNKS

! number of message types
  integer NUM_MSG_TYPES

  integer NPROC_ONE_DIRECTION

! ************** subroutine starts here **************

! number of corners and faces shared between chunks and number of message types
  if(NCHUNKS == 1 .or. NCHUNKS == 2) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 1
  else if(NCHUNKS == 3) then
    NCORNERSCHUNKS = 1
    NUM_FACES = 1
    NUM_MSG_TYPES = 3
  else if(NCHUNKS == 6) then
    NCORNERSCHUNKS = 8
    NUM_FACES = 4
    NUM_MSG_TYPES = 3
  else
    call exit_MPI(myrank,'number of chunks must be either 1, 2, 3 or 6')
  endif

! if more than one chunk then same number of processors in each direction
  NPROC_ONE_DIRECTION = NPROC_XI

! total number of messages corresponding to these common faces
  NUMMSGS_FACES = NPROC_ONE_DIRECTION*NUM_FACES*NUM_MSG_TYPES

! check that there is more than one chunk, otherwise nothing to do
  if(NCHUNKS == 1) return

! same number of GLL points in each direction for several chunks
  if(NGLLY /= NGLLX) call exit_MPI(myrank,'must have NGLLY = NGLLX for several chunks')

! allocate arrays for faces
  allocate(iproc_sender(NUMMSGS_FACES))
  allocate(iproc_receiver(NUMMSGS_FACES))
  allocate(npoin2D_send(NUMMSGS_FACES))
  allocate(npoin2D_receive(NUMMSGS_FACES))

! allocate array for corners
  allocate(iprocscorners(3,NCORNERSCHUNKS))
  allocate(itypecorner(3,NCORNERSCHUNKS))

! clear arrays allocated
  iproc_sender(:) = 0
  iproc_receiver(:) = 0
  npoin2D_send(:) = 0
  npoin2D_receive(:) = 0
  iprocscorners(:,:) = 0
  itypecorner(:,:) = 0

  if(myrank == 0) then
    write(IMAIN,*) 'There is a total of ',NUMMSGS_FACES,' messages to assemble faces between chunks'
    write(IMAIN,*)
  endif

  imsg = 0

  if(myrank == 0) then

! get the base pathname for output files
    call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! file to store the list of processors for each message for faces
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/list_messages_faces.txt',status='unknown',action='write')

  endif

!!!!!!!!!! DK DK for merged version: beginning of "faces" section here
!!!!!!!!!! DK DK for merged version: beginning of "faces" section here
!!!!!!!!!! DK DK for merged version: beginning of "faces" section here
!!!!!!!!!! DK DK for merged version: beginning of "faces" section here
!!!!!!!!!! DK DK for merged version: beginning of "faces" section here
!!!!!!!!!! DK DK for merged version: beginning of "faces" section here
!!!!!!!!!! DK DK for merged version: beginning of "faces" section here

! create theoretical communication pattern
  do imsg_type = 1,NUM_MSG_TYPES
    do iside = 1,NUM_FACES
      do iproc_loop = 0,NPROC_ONE_DIRECTION-1

! create a new message
! we know there can be no deadlock with this scheme
! because the three types of messages are independent
        imsg = imsg + 1

! check that current message number is correct
        if(imsg > NUMMSGS_FACES) call exit_MPI(myrank,'incorrect message number')

        if(myrank == 0) write(IMAIN,*) 'Generating message ',imsg,' for faces out of ',NUMMSGS_FACES

! we know there is the same number of slices in both directions
        iproc_xi_loop = iproc_loop
        iproc_eta_loop = iproc_loop

! take care of local frame inversions between chunks
        iproc_xi_loop_inv = NPROC_ONE_DIRECTION - iproc_loop - 1
        iproc_eta_loop_inv = NPROC_ONE_DIRECTION - iproc_loop - 1


! define the 12 different messages

! message type M1
        if(imsg_type == 1) then

          if(iside == 1) then
            ichunk_send = CHUNK_AB
            iproc_xi_send = 0
            iproc_eta_send = iproc_eta_loop
            iproc_edge_send = XI_MIN
            ichunk_receive = CHUNK_AC
            iproc_xi_receive = NPROC_XI-1
            iproc_eta_receive = iproc_eta_loop
            iproc_edge_receive = XI_MAX
          endif

          if(iside == 2) then
            ichunk_send = CHUNK_AB
            iproc_xi_send = NPROC_XI-1
            iproc_eta_send = iproc_eta_loop
            iproc_edge_send = XI_MAX
            ichunk_receive = CHUNK_AC_ANTIPODE
            iproc_xi_receive = 0
            iproc_eta_receive = iproc_eta_loop
            iproc_edge_receive = XI_MIN
          endif

          if(iside == 3) then
            ichunk_send = CHUNK_AC_ANTIPODE
            iproc_xi_send = NPROC_XI-1
            iproc_eta_send = iproc_eta_loop
            iproc_edge_send = XI_MAX
            ichunk_receive = CHUNK_AB_ANTIPODE
            iproc_xi_receive = 0
            iproc_eta_receive = iproc_eta_loop
            iproc_edge_receive = XI_MIN
          endif

          if(iside == 4) then
            ichunk_send = CHUNK_AC
            iproc_xi_send = 0
            iproc_eta_send = iproc_eta_loop
            iproc_edge_send = XI_MIN
            ichunk_receive = CHUNK_AB_ANTIPODE
            iproc_xi_receive = NPROC_XI-1
            iproc_eta_receive = iproc_eta_loop
            iproc_edge_receive = XI_MAX
          endif

        endif

! message type M2
        if(imsg_type == 2) then

          if(iside == 1) then
            ichunk_send = CHUNK_AB
            iproc_xi_send = iproc_xi_loop
            iproc_eta_send = NPROC_ETA-1
            iproc_edge_send = ETA_MAX
            ichunk_receive = CHUNK_BC
            iproc_xi_receive = NPROC_XI-1
            iproc_eta_receive = iproc_eta_loop
            iproc_edge_receive = XI_MAX
          endif

          if(iside == 2) then
            ichunk_send = CHUNK_AB
            iproc_xi_send = iproc_xi_loop
            iproc_eta_send = 0
            iproc_edge_send = ETA_MIN
            ichunk_receive = CHUNK_BC_ANTIPODE
            iproc_xi_receive = NPROC_XI-1
            iproc_eta_receive = iproc_eta_loop_inv
            iproc_edge_receive = XI_MAX
          endif

          if(iside == 3) then
            ichunk_send = CHUNK_BC
            iproc_xi_send = 0
            iproc_eta_send = iproc_eta_loop
            iproc_edge_send = XI_MIN
            ichunk_receive = CHUNK_AB_ANTIPODE
            iproc_xi_receive = iproc_xi_loop_inv
            iproc_eta_receive = NPROC_ETA-1
            iproc_edge_receive = ETA_MAX
          endif

          if(iside == 4) then
            ichunk_send = CHUNK_BC_ANTIPODE
            iproc_xi_send = 0
            iproc_eta_send = iproc_eta_loop
            iproc_edge_send = XI_MIN
            ichunk_receive = CHUNK_AB_ANTIPODE
            iproc_xi_receive = iproc_xi_loop
            iproc_eta_receive = 0
            iproc_edge_receive = ETA_MIN
          endif

        endif

! message type M3
        if(imsg_type == 3) then

          if(iside == 1) then
            ichunk_send = CHUNK_AC
            iproc_xi_send = iproc_xi_loop
            iproc_eta_send = NPROC_ETA-1
            iproc_edge_send = ETA_MAX
            ichunk_receive = CHUNK_BC
            iproc_xi_receive = iproc_xi_loop
            iproc_eta_receive = 0
            iproc_edge_receive = ETA_MIN
          endif

          if(iside == 2) then
            ichunk_send = CHUNK_BC
            iproc_xi_send = iproc_xi_loop
            iproc_eta_send = NPROC_ETA-1
            iproc_edge_send = ETA_MAX
            ichunk_receive = CHUNK_AC_ANTIPODE
            iproc_xi_receive = iproc_xi_loop_inv
            iproc_eta_receive = NPROC_ETA-1
            iproc_edge_receive = ETA_MAX
          endif

          if(iside == 3) then
            ichunk_send = CHUNK_AC_ANTIPODE
            iproc_xi_send = iproc_xi_loop
            iproc_eta_send = 0
            iproc_edge_send = ETA_MIN
            ichunk_receive = CHUNK_BC_ANTIPODE
            iproc_xi_receive = iproc_xi_loop_inv
            iproc_eta_receive = 0
            iproc_edge_receive = ETA_MIN
          endif

          if(iside == 4) then
            ichunk_send = CHUNK_AC
            iproc_xi_send = iproc_xi_loop
            iproc_eta_send = 0
            iproc_edge_send = ETA_MIN
            ichunk_receive = CHUNK_BC_ANTIPODE
            iproc_xi_receive = iproc_xi_loop
            iproc_eta_receive = NPROC_ETA-1
            iproc_edge_receive = ETA_MAX
          endif

        endif


! store addressing generated
        iproc_sender(imsg) = addressing(ichunk_send,iproc_xi_send,iproc_eta_send)
        iproc_receiver(imsg) = addressing(ichunk_receive,iproc_xi_receive,iproc_eta_receive)

! check that sender/receiver pair is ordered
        if(iproc_sender(imsg) > iproc_receiver(imsg)) call exit_MPI(myrank,'incorrect order in sender/receiver pair')

! save message type and pair of processors in list of messages
        if(myrank == 0) write(IOUT,*) imsg_type,iproc_sender(imsg),iproc_receiver(imsg)

! loop on sender/receiver (1=sender 2=receiver)
        do imode_comm=1,2

          if(imode_comm == 1) then
            iproc = iproc_sender(imsg)
            iedge = iproc_edge_send

          else if(imode_comm == 2) then
            iproc = iproc_receiver(imsg)
            iedge = iproc_edge_receive

          else
            call exit_MPI(myrank,'incorrect communication mode')
          endif

! only do this if current processor is the right one for MPI version
          if(iproc == myrank) then

! determine chunk number and local slice coordinates using addressing
            ichunk = ichunk_slice(iproc)
            iproc_xi = iproc_xi_slice(iproc)
            iproc_eta = iproc_eta_slice(iproc)

! problem if not on edges
            if(iproc_xi /= 0 .and. iproc_xi /= NPROC_XI-1 .and. &
              iproc_eta /= 0 .and. iproc_eta /= NPROC_ETA-1) call exit_MPI(myrank,'slice not on any edge')

            nglob=nglob_ori
! check that iboolmax=nglob

! end of section done only if right processor for MPI
          endif

! end of loop on sender/receiver
        enddo

! end of loops on all the messages
      enddo
    enddo
  enddo

  if(myrank == 0) close(IOUT)

! check that total number of messages is correct
  if(imsg /= NUMMSGS_FACES) call exit_MPI(myrank,'incorrect total number of messages')

!
!---- check that number of points detected is the same for sender and receiver
!

! synchronize all the processes to make sure all the buffers are ready
  call MPI_BARRIER(MPI_COMM_WORLD,ier)

!!!!!!!!!! DK DK for merged version: beginning of "corner" section here
!!!!!!!!!! DK DK for merged version: beginning of "corner" section here
!!!!!!!!!! DK DK for merged version: beginning of "corner" section here
!!!!!!!!!! DK DK for merged version: beginning of "corner" section here
!!!!!!!!!! DK DK for merged version: beginning of "corner" section here
!!!!!!!!!! DK DK for merged version: beginning of "corner" section here
!!!!!!!!!! DK DK for merged version: beginning of "corner" section here

!
!---- generate the 8 message patterns sharing a corner of valence 3
!

! to avoid problem at compile time, use bigger array with fixed dimension
  addressing_big(:,:,:) = 0
  addressing_big(1:NCHUNKS,:,:) = addressing(1:NCHUNKS,:,:)

  ichunk = 1
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AB,0,NPROC_ETA-1)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_AC,NPROC_XI-1,NPROC_ETA-1)
! this line is ok even for NCHUNKS = 2
  iprocscorners(3,ichunk) = addressing_big(CHUNK_BC,NPROC_XI-1,0)

  itypecorner(1,ichunk) = ILOWERUPPER
  itypecorner(2,ichunk) = IUPPERUPPER
  itypecorner(3,ichunk) = IUPPERLOWER

!! DK DK UGLY in the future, should also assemble second corner when NCHUNKS = 2
!! DK DK UGLY for now we only assemble one corner for simplicity
!! DK DK UGLY formally this is incorrect and should be changed in the future
!! DK DK UGLY in practice this trick works fine

! this only if more than 3 chunks
  if(NCHUNKS > 3) then

  ichunk = 2
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AB,NPROC_XI-1,0)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_AC_ANTIPODE,0,0)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_BC_ANTIPODE,NPROC_XI-1,0)

  itypecorner(1,ichunk) = IUPPERLOWER
  itypecorner(2,ichunk) = ILOWERLOWER
  itypecorner(3,ichunk) = IUPPERLOWER

  ichunk = 3
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AB,0,0)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_AC,NPROC_XI-1,0)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_BC_ANTIPODE,NPROC_XI-1,NPROC_ETA-1)

  itypecorner(1,ichunk) = ILOWERLOWER
  itypecorner(2,ichunk) = IUPPERLOWER
  itypecorner(3,ichunk) = IUPPERUPPER

  ichunk = 4
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AB,NPROC_XI-1,NPROC_ETA-1)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_BC,NPROC_XI-1,NPROC_ETA-1)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_AC_ANTIPODE,0,NPROC_ETA-1)

  itypecorner(1,ichunk) = IUPPERUPPER
  itypecorner(2,ichunk) = IUPPERUPPER
  itypecorner(3,ichunk) = ILOWERUPPER

  ichunk = 5
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AC,0,0)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_BC_ANTIPODE,0,NPROC_ETA-1)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_AB_ANTIPODE,NPROC_XI-1,0)

  itypecorner(1,ichunk) = ILOWERLOWER
  itypecorner(2,ichunk) = ILOWERUPPER
  itypecorner(3,ichunk) = IUPPERLOWER

  ichunk = 6
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AC_ANTIPODE,NPROC_XI-1,0)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_BC_ANTIPODE,0,0)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_AB_ANTIPODE,0,0)

  itypecorner(1,ichunk) = IUPPERLOWER
  itypecorner(2,ichunk) = ILOWERLOWER
  itypecorner(3,ichunk) = ILOWERLOWER

  ichunk = 7
  iprocscorners(1,ichunk) = addressing_big(CHUNK_AC,0,NPROC_ETA-1)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_BC,0,0)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_AB_ANTIPODE,NPROC_XI-1,NPROC_ETA-1)

  itypecorner(1,ichunk) = ILOWERUPPER
  itypecorner(2,ichunk) = ILOWERLOWER
  itypecorner(3,ichunk) = IUPPERUPPER

  ichunk = 8
  iprocscorners(1,ichunk) = addressing_big(CHUNK_BC,0,NPROC_ETA-1)
  iprocscorners(2,ichunk) = addressing_big(CHUNK_AC_ANTIPODE,NPROC_XI-1,NPROC_ETA-1)
  iprocscorners(3,ichunk) = addressing_big(CHUNK_AB_ANTIPODE,0,NPROC_ETA-1)

  itypecorner(1,ichunk) = ILOWERUPPER
  itypecorner(2,ichunk) = IUPPERUPPER
  itypecorner(3,ichunk) = ILOWERUPPER

  endif

! file to store the list of processors for each message for corners
  if(myrank == 0) open(unit=IOUT,file=trim(OUTPUT_FILES)//'/list_messages_corners.txt',status='unknown',action='write')

! loop over all the messages to create the addressing
  do imsg = 1,NCORNERSCHUNKS

  if(myrank == 0) write(IMAIN,*) 'Generating message ',imsg,' for corners out of ',NCORNERSCHUNKS

! save triplet of processors in list of messages
  if(myrank == 0) write(IOUT,*) iprocscorners(1,imsg),iprocscorners(2,imsg),iprocscorners(3,imsg)

! loop on the three processors of a given corner
  do imember_corner = 1,3

! only do this if current processor is the right one for MPI version
! this line is ok even for NCHUNKS = 2
  if(iprocscorners(imember_corner,imsg) == myrank) then

! pick the correct 1D buffer
! this scheme works fine even if NPROC_XI = NPROC_ETA = 1
  if(itypecorner(imember_corner,imsg) == ILOWERLOWER) then
!! DK DK suppressed for merged    filename_in = prname(1:len_trim(prname))//'ibool1D_leftxi_lefteta.txt'
    NGLOB1D_RADIAL_my_corner = NGLOB1D_RADIAL_CORNER(iregion_code,1)
  else if(itypecorner(imember_corner,imsg) == ILOWERUPPER) then
!! DK DK suppressed for merged    filename_in = prname(1:len_trim(prname))//'ibool1D_leftxi_righteta.txt'
    NGLOB1D_RADIAL_my_corner = NGLOB1D_RADIAL_CORNER(iregion_code,4)
  else if(itypecorner(imember_corner,imsg) == IUPPERLOWER) then
!! DK DK suppressed for merged    filename_in = prname(1:len_trim(prname))//'ibool1D_rightxi_lefteta.txt'
    NGLOB1D_RADIAL_my_corner = NGLOB1D_RADIAL_CORNER(iregion_code,2)
  else if(itypecorner(imember_corner,imsg) == IUPPERUPPER) then
!! DK DK suppressed for merged    filename_in = prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt'
    NGLOB1D_RADIAL_my_corner = NGLOB1D_RADIAL_CORNER(iregion_code,3)
  else
    call exit_MPI(myrank,'incorrect corner coordinates')
  endif

! end of section done only if right processor for MPI
  endif

  enddo

  enddo

  if(myrank == 0) close(IOUT)

  end subroutine create_list_files_chunks

