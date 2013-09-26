!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  6 . 0
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and CNRS / INRIA / University of Pau, France
! (c) Princeton University and CNRS / INRIA / University of Pau
!                            August 2013
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

! code to check that all the 2D buffers between chunk faces are okay

  program check_buffers_faces_chunks

  use constants
  use shared_parameters

  implicit none

  integer imsg

  integer npoin2D_sender,npoin2D_receiver
  integer iboolsend,iboolreceive,ipoin2D
  integer iregion_code,iproc

! number of faces between chunks
  integer NUM_FACES,NUMMSGS_FACES

! number of message types
  integer NUM_MSG_TYPES

  double precision xsend,ysend,zsend
  double precision xreceive,yreceive,zreceive
  double precision diff

  integer NPROC_ONE_DIRECTION

! communication pattern for faces between chunks
  integer, dimension(:), allocatable :: iprocfrom_faces,iprocto_faces,imsg_type

  character(len=150) filename,prname

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers between chunk faces'
  print *

! read the parameter file and compute additional parameters
  call read_compute_parameters()

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

  print *
  print *,'There are ',NPROCTOT,' slices numbered from 0 to ',NPROCTOT-1
  print *,'There are ',NCHUNKS,' chunks'
  print *,'There are ',NPROC_XI,' slices along xi in each chunk'
  print *,'There are ',NPROC_ETA,' slices along eta in each chunk'
  print *

! number of corners and faces shared between chunks and number of message types
  if(NCHUNKS == 1 .or. NCHUNKS == 2) then
    NUM_FACES = 1
    NUM_MSG_TYPES = 1
  else if(NCHUNKS == 3) then
    NUM_FACES = 1
    NUM_MSG_TYPES = 3
  else if(NCHUNKS == 6) then
    NUM_FACES = 4
    NUM_MSG_TYPES = 3
  else
    stop 'number of chunks must be either 1, 2, 3 or 6'
  endif

! if more than one chunk then same number of processors in each direction
  NPROC_ONE_DIRECTION = NPROC_XI

! total number of messages corresponding to these common faces
  NUMMSGS_FACES = NPROC_ONE_DIRECTION*NUM_FACES*NUM_MSG_TYPES

  if(NCHUNKS == 1) stop 'only one chunk, nothing to check'

  print *,'There are ',NUMMSGS_FACES,' messages to assemble all the faces'
  print *

! allocate array for messages for faces
  allocate(iprocfrom_faces(NUMMSGS_FACES))
  allocate(iprocto_faces(NUMMSGS_FACES))
  allocate(imsg_type(NUMMSGS_FACES))

! file with the list of processors for each message for faces
  open(unit=IIN,file=trim(OUTPUT_FILES)//'/list_messages_faces.txt',status='old',action='read')
  do imsg = 1,NUMMSGS_FACES
  read(IIN,*) imsg_type(imsg),iprocfrom_faces(imsg),iprocto_faces(imsg)
  if      (iprocfrom_faces(imsg) < 0 &
        .or. iprocto_faces(imsg) < 0 &
        .or. iprocfrom_faces(imsg) > NPROCTOT-1 &
        .or. iprocto_faces(imsg) > NPROCTOT-1) &
    stop 'incorrect chunk faces numbering'
  if (imsg_type(imsg) < 1 .or. imsg_type(imsg) > 3) &
    stop 'incorrect message type labeling'
  enddo
  close(IIN)

! loop over all the regions of the mesh
  do iregion_code = 1,MAX_NUM_REGIONS

  print *
  print *,' ********* checking region ',iregion_code,' *********'
  print *

! loop on all the messages between faces
  do imsg = 1,NUMMSGS_FACES

  print *
  print *,'Checking message ',imsg,' out of ',NUMMSGS_FACES

! read 2-D buffer for the sender and the receiver
  write(filename,"('buffer_faces_chunks_sender_msg',i6.6,'.txt')") imsg
  iproc = iprocfrom_faces(imsg)
  call create_serial_name_database(prname,iproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  open(unit=34,file=prname(1:len_trim(prname))//filename,status='old',action='read')

  write(filename,"('buffer_faces_chunks_receiver_msg',i6.6,'.txt')") imsg
  iproc = iprocto_faces(imsg)
  call create_serial_name_database(prname,iproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  open(unit=35,file=prname(1:len_trim(prname))//filename,status='old',action='read')

  write(*,*) 'reading MPI 2D buffer for sender'
  read(34,*) npoin2D_sender
  read(35,*) npoin2D_receiver

! check that number of points is the same in both buffers
  if(npoin2D_sender /= npoin2D_receiver) &
        stop 'different number of points in the two buffers'

  print *,'this message contains ',npoin2D_sender,' points'

! check all the points based upon their coordinates
  do ipoin2D = 1,npoin2D_sender
    read(34,*) iboolsend,xsend,ysend,zsend
    read(35,*) iboolreceive,xreceive,yreceive,zreceive

    diff = dmax1(dabs(xsend-xreceive),dabs(ysend-yreceive),dabs(zsend-zreceive))
    if(diff > 0.0000001d0) then
      print *,'different : ',ipoin2D,iboolsend,iboolreceive,diff
      print *,'xsend,xreceive = ',xsend,xreceive
      print *,'ysend,yreceive = ',ysend,yreceive
      print *,'zsend,zreceive = ',zsend,zreceive
      stop 'error: different'
    endif

  enddo

  enddo

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_faces_chunks

