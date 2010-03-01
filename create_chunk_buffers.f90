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

! subroutine to create MPI buffers to assemble between chunks

  subroutine create_chunk_buffers(iregion_code,nspec,ibool,idoubling, &
                                  xstore,ystore,zstore, &
                                  nglob_ori, &
                                  NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
                                  NPROC_XI,NPROC_ETA,NPROC,NPROCTOT, &
                                  NGLOB1D_RADIAL_CORNER,NGLOB1D_RADIAL_MAX, &
                                  NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                                  myrank,LOCAL_PATH,addressing, &
                                  ichunk_slice,iproc_xi_slice,iproc_eta_slice,NCHUNKS)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"
  include "precision.h"

  integer, dimension(MAX_NUM_REGIONS,NB_SQUARE_CORNERS) :: NGLOB1D_RADIAL_CORNER

  integer nglob,nglob_ori
  integer NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX
  integer NPROC,NPROC_XI,NPROC_ETA,NPROCTOT,NGLOB1D_RADIAL_MAX,NGLOB1D_RADIAL
  integer NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer nspec
  integer myrank,NCHUNKS

! arrays with the mesh
  double precision xstore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision ystore(NGLLX,NGLLY,NGLLZ,nspec)
  double precision zstore(NGLLX,NGLLY,NGLLZ,nspec)

  character(len=150) OUTPUT_FILES,LOCAL_PATH,ERR_MSG

! array with the local to global mapping per slice
  integer ibool(NGLLX,NGLLY,NGLLZ,nspec)

  integer idoubling(nspec)

! mask for ibool to mark points already found
  logical, dimension(:), allocatable ::  mask_ibool

! array to store points selected for the chunk face buffer
  integer NGLOB2DMAX_XY
  integer, dimension(:), allocatable :: ibool_selected

  double precision, dimension(:), allocatable :: xstore_selected,ystore_selected,zstore_selected

! arrays for sorting routine
  integer, dimension(:), allocatable :: ind,ninseg,iglob,locval,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: work

! pairs generated theoretically
! four sides for each of the three types of messages
  integer, dimension(:), allocatable :: iproc_sender,iproc_receiver,npoin2D_send,npoin2D_receive

! 1D buffers to remove points belonging to corners
  integer ibool1D_leftxi_lefteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D_rightxi_lefteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D_leftxi_righteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D_rightxi_righteta(NGLOB1D_RADIAL_MAX)
  integer ibool1D(NGLOB1D_RADIAL_MAX)
  double precision xread1D(NGLOB1D_RADIAL_MAX)
  double precision yread1D(NGLOB1D_RADIAL_MAX)
  double precision zread1D(NGLOB1D_RADIAL_MAX)
  double precision xdummy,ydummy,zdummy
  integer ipoin1D

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

! boundary parameters per slice
  integer nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, njunk
  integer ibelm_xmin(NSPEC2DMAX_XMIN_XMAX),ibelm_xmax(NSPEC2DMAX_XMIN_XMAX)
  integer ibelm_ymin(NSPEC2DMAX_YMIN_YMAX),ibelm_ymax(NSPEC2DMAX_YMIN_YMAX)

  integer npoin2D,npoin2D_send_local,npoin2D_receive_local

  integer i,j,k,ispec,ispec2D,ipoin2D,ier

! current message number
  integer imsg

! names of the data files for all the processors in MPI
  character(len=150) prname,filename_in,filename_out

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

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '----- creating chunk buffers -----'
    write(IMAIN,*)
    write(IMAIN,*) 'There are ',NPROC_XI,' slices along xi in each chunk'
    write(IMAIN,*) 'There are ',NPROC_ETA,' slices along eta in each chunk'
    write(IMAIN,*) 'There is a total of ',NPROC,' slices in each chunk'
    write(IMAIN,*) 'There are ',NCHUNKS,' chunks'
    write(IMAIN,*) 'There is a total of ',NPROCTOT,' slices in all the chunks'
    write(IMAIN,*)
  endif

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

! define maximum size for message buffers
  NGLOB2DMAX_XY = max(NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX)

! allocate arrays for message buffers with maximum size
  allocate(ibool_selected(NGLOB2DMAX_XY))
  allocate(xstore_selected(NGLOB2DMAX_XY))
  allocate(ystore_selected(NGLOB2DMAX_XY))
  allocate(zstore_selected(NGLOB2DMAX_XY))
  allocate(ind(NGLOB2DMAX_XY))
  allocate(ninseg(NGLOB2DMAX_XY))
  allocate(iglob(NGLOB2DMAX_XY))
  allocate(locval(NGLOB2DMAX_XY))
  allocate(ifseg(NGLOB2DMAX_XY))
  allocate(iwork(NGLOB2DMAX_XY))
  allocate(work(NGLOB2DMAX_XY))


! allocate mask for ibool
  allocate(mask_ibool(nglob_ori))

  imsg = 0

  if(myrank == 0) then

! get the base pathname for output files
    call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! file to store the list of processors for each message for faces
    open(unit=IOUT,file=trim(OUTPUT_FILES)//'/list_messages_faces.txt',status='unknown')

  endif

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
            write(filename_out,"('buffer_faces_chunks_sender_msg',i6.6,'.txt')") imsg
          else if(imode_comm == 2) then
            iproc = iproc_receiver(imsg)
            iedge = iproc_edge_receive
            write(filename_out,"('buffer_faces_chunks_receiver_msg',i6.6,'.txt')") imsg
          else
            call exit_MPI(myrank,'incorrect communication mode')
          endif

! only do this if current processor is the right one for MPI version
          if(iproc == myrank) then

! create the name of the database for each slice
            call create_name_database(prname,iproc,iregion_code,LOCAL_PATH)

! open file for 2D buffer
            open(unit=IOUT_BUFFERS,file=prname(1:len_trim(prname))//filename_out,status='unknown')

! determine chunk number and local slice coordinates using addressing
            ichunk = ichunk_slice(iproc)
            iproc_xi = iproc_xi_slice(iproc)
            iproc_eta = iproc_eta_slice(iproc)

! problem if not on edges
            if(iproc_xi /= 0 .and. iproc_xi /= NPROC_XI-1 .and. &
              iproc_eta /= 0 .and. iproc_eta /= NPROC_ETA-1) call exit_MPI(myrank,'slice not on any edge')

            nglob=nglob_ori
! check that iboolmax=nglob

            if(minval(ibool(:,:,:,1:nspec)) /= 1 .or. maxval(ibool(:,:,:,1:nspec)) /= nglob) &
              call exit_MPI(myrank,ERR_MSG)

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! read boundary parameters

            open(unit=IIN,file=prname(1:len_trim(prname))//'boundary.bin',status='old',action='read',form='unformatted')
            read(IIN) nspec2D_xmin
            read(IIN) nspec2D_xmax
            read(IIN) nspec2D_ymin
            read(IIN) nspec2D_ymax
            read(IIN) njunk
            read(IIN) njunk

            read(IIN) ibelm_xmin
            read(IIN) ibelm_xmax
            read(IIN) ibelm_ymin
            read(IIN) ibelm_ymax
            close(IIN)

! read 1D buffers to remove corner points
            open(unit=IIN,file=prname(1:len_trim(prname))//'ibool1D_leftxi_lefteta.txt',status='old',action='read')
            do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,1)
              read(IIN,*) ibool1D_leftxi_lefteta(ipoin1D),xdummy,ydummy,zdummy
            enddo
            close(IIN)

            open(unit=IIN,file=prname(1:len_trim(prname))//'ibool1D_rightxi_lefteta.txt',status='old',action='read')
            do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,2)
              read(IIN,*) ibool1D_rightxi_lefteta(ipoin1D),xdummy,ydummy,zdummy
            enddo
            close(IIN)

            open(unit=IIN,file=prname(1:len_trim(prname))//'ibool1D_leftxi_righteta.txt',status='old',action='read')
            do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,4)
              read(IIN,*) ibool1D_leftxi_righteta(ipoin1D),xdummy,ydummy,zdummy
            enddo
            close(IIN)

            open(unit=IIN,file=prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt',status='old',action='read')
            do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,3)
              read(IIN,*) ibool1D_rightxi_righteta(ipoin1D),xdummy,ydummy,zdummy
            enddo
            close(IIN)

! erase logical mask
            mask_ibool(:) = .false.

            npoin2D = 0

! create all the points on each face (no duplicates, but not sorted)

! xmin
            if(iedge == XI_MIN) then

! mark corner points to remove them if needed
              if(iproc_eta == 0) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,1)
                  mask_ibool(ibool1D_leftxi_lefteta(ipoin1D)) = .true.
                enddo
              endif

              if(iproc_eta == NPROC_ETA-1) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,4)
                  mask_ibool(ibool1D_leftxi_righteta(ipoin1D)) = .true.
                enddo
              endif

              do ispec2D=1,nspec2D_xmin
                ispec=ibelm_xmin(ispec2D)

! remove central cube for chunk buffers
                if(idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

                i=1
                do k=1,NGLLZ
                  do j=1,NGLLY
                    if(.not. mask_ibool(ibool(i,j,k,ispec))) then
! mask and store points found
                      mask_ibool(ibool(i,j,k,ispec)) = .true.
                      npoin2D = npoin2D + 1
                      if(npoin2D > NGLOB2DMAX_XMIN_XMAX) call exit_MPI(myrank,'incorrect 2D point number in xmin')
                      ibool_selected(npoin2D) = ibool(i,j,k,ispec)

                      xstore_selected(npoin2D) = xstore(i,j,k,ispec)
                      ystore_selected(npoin2D) = ystore(i,j,k,ispec)
                      zstore_selected(npoin2D) = zstore(i,j,k,ispec)
                    endif
                  enddo
                enddo
              enddo

! xmax
            else if(iedge == XI_MAX) then

! mark corner points to remove them if needed

              if(iproc_eta == 0) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,2)
                  mask_ibool(ibool1D_rightxi_lefteta(ipoin1D)) = .true.
                enddo
              endif

              if(iproc_eta == NPROC_ETA-1) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,3)
                  mask_ibool(ibool1D_rightxi_righteta(ipoin1D)) = .true.
                enddo
              endif

              do ispec2D=1,nspec2D_xmax
                ispec=ibelm_xmax(ispec2D)

! remove central cube for chunk buffers
                if(idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

                i=NGLLX
                do k=1,NGLLZ
                  do j=1,NGLLY
                    if(.not. mask_ibool(ibool(i,j,k,ispec))) then
! mask and store points found
                      mask_ibool(ibool(i,j,k,ispec)) = .true.
                      npoin2D = npoin2D + 1
                      if(npoin2D > NGLOB2DMAX_XMIN_XMAX) call exit_MPI(myrank,'incorrect 2D point number in xmax')
                      ibool_selected(npoin2D) = ibool(i,j,k,ispec)

                      xstore_selected(npoin2D) = xstore(i,j,k,ispec)
                      ystore_selected(npoin2D) = ystore(i,j,k,ispec)
                      zstore_selected(npoin2D) = zstore(i,j,k,ispec)
                    endif
                  enddo
                enddo
              enddo

! ymin
            else if(iedge == ETA_MIN) then

! mark corner points to remove them if needed

              if(iproc_xi == 0) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,1)
                  mask_ibool(ibool1D_leftxi_lefteta(ipoin1D)) = .true.
                enddo
              endif

              if(iproc_xi == NPROC_XI-1) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,2)
                  mask_ibool(ibool1D_rightxi_lefteta(ipoin1D)) = .true.
                enddo
              endif

              do ispec2D=1,nspec2D_ymin
                ispec=ibelm_ymin(ispec2D)

! remove central cube for chunk buffers
                if(idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

                j=1
                do k=1,NGLLZ
                  do i=1,NGLLX
                    if(.not. mask_ibool(ibool(i,j,k,ispec))) then
! mask and store points found
                      mask_ibool(ibool(i,j,k,ispec)) = .true.
                      npoin2D = npoin2D + 1
                      if(npoin2D > NGLOB2DMAX_YMIN_YMAX) call exit_MPI(myrank,'incorrect 2D point number in ymin')
                      ibool_selected(npoin2D) = ibool(i,j,k,ispec)

                      xstore_selected(npoin2D) = xstore(i,j,k,ispec)
                      ystore_selected(npoin2D) = ystore(i,j,k,ispec)
                      zstore_selected(npoin2D) = zstore(i,j,k,ispec)
                    endif
                  enddo
                enddo
              enddo

! ymax
            else if(iedge == ETA_MAX) then

! mark corner points to remove them if needed

              if(iproc_xi == 0) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,4)
                  mask_ibool(ibool1D_leftxi_righteta(ipoin1D)) = .true.
                enddo
              endif

              if(iproc_xi == NPROC_XI-1) then
                do ipoin1D = 1,NGLOB1D_RADIAL_CORNER(iregion_code,3)
                  mask_ibool(ibool1D_rightxi_righteta(ipoin1D)) = .true.
                enddo
              endif

              do ispec2D=1,nspec2D_ymax
                ispec=ibelm_ymax(ispec2D)

! remove central cube for chunk buffers
                if(idoubling(ispec) == IFLAG_MIDDLE_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_BOTTOM_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_TOP_CENTRAL_CUBE .or. &
                  idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle

                j=NGLLY
                do k=1,NGLLZ
                  do i=1,NGLLX
                    if(.not. mask_ibool(ibool(i,j,k,ispec))) then
! mask and store points found
                      mask_ibool(ibool(i,j,k,ispec)) = .true.
                      npoin2D = npoin2D + 1
                      if(npoin2D > NGLOB2DMAX_YMIN_YMAX) call exit_MPI(myrank,'incorrect 2D point number in ymax')
                      ibool_selected(npoin2D) = ibool(i,j,k,ispec)

                      xstore_selected(npoin2D) = xstore(i,j,k,ispec)
                      ystore_selected(npoin2D) = ystore(i,j,k,ispec)
                      zstore_selected(npoin2D) = zstore(i,j,k,ispec)
                    endif
                  enddo
                enddo
              enddo

            else

              call exit_MPI(myrank,'incorrect edge code')
            endif

! sort buffer obtained to be conforming with neighbor in other chunk
! sort on x, y and z, the other arrays will be swapped as well

            call sort_array_coordinates(npoin2D,xstore_selected,ystore_selected,zstore_selected, &
              ibool_selected,iglob,locval,ifseg,nglob,ind,ninseg,iwork,work)

! check that no duplicate has been detected
            if(nglob /= npoin2D) call exit_MPI(myrank,'duplicates detected in buffer')

! write list of selected points to output buffer
            write(IOUT_BUFFERS,*) npoin2D
            do ipoin2D = 1,npoin2D
                write(IOUT_BUFFERS,*) ibool_selected(ipoin2D), &
                  xstore_selected(ipoin2D),ystore_selected(ipoin2D),zstore_selected(ipoin2D)
            enddo

            close(IOUT_BUFFERS)

! store result to compare number of points for sender and for receiver
            if(imode_comm == 1) then
              npoin2D_send(imsg) = npoin2D
            else
              npoin2D_receive(imsg) = npoin2D
            endif

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

! gather information about all the messages on all processes
  do imsg = 1,NUMMSGS_FACES

!     gather number of points for sender
      npoin2D_send_local = npoin2D_send(imsg)
      call MPI_BCAST(npoin2D_send_local,1,MPI_INTEGER,iproc_sender(imsg),MPI_COMM_WORLD,ier)
      if(myrank /= iproc_sender(imsg)) npoin2D_send(imsg) = npoin2D_send_local

!     gather number of points for receiver
      npoin2D_receive_local = npoin2D_receive(imsg)
      call MPI_BCAST(npoin2D_receive_local,1,MPI_INTEGER,iproc_receiver(imsg),MPI_COMM_WORLD,ier)
      if(myrank /= iproc_receiver(imsg)) npoin2D_receive(imsg) = npoin2D_receive_local

  enddo

! check the number of points
  do imsg = 1,NUMMSGS_FACES
    if(npoin2D_send(imsg) /= npoin2D_receive(imsg)) &
        call exit_MPI(myrank,'incorrect number of points for sender/receiver pair detected')
  enddo
  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'all the messages for chunk faces have the right size'
    write(IMAIN,*)
  endif

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
  if(myrank == 0) open(unit=IOUT,file=trim(OUTPUT_FILES)//'/list_messages_corners.txt',status='unknown')

! loop over all the messages to create the addressing
  do imsg = 1,NCORNERSCHUNKS

  if(myrank == 0) write(IMAIN,*) 'Generating message ',imsg,' for corners out of ',NCORNERSCHUNKS

! save triplet of processors in list of messages
  if(myrank == 0) write(IOUT,*) iprocscorners(1,imsg),iprocscorners(2,imsg),iprocscorners(3,imsg)

! loop on the three processors of a given corner
  do imember_corner = 1,3

    if(imember_corner == 1) then
      write(filename_out,"('buffer_corners_chunks_master_msg',i6.6,'.txt')") imsg
    else if(imember_corner == 2) then
      write(filename_out,"('buffer_corners_chunks_worker1_msg',i6.6,'.txt')") imsg
    else
      write(filename_out,"('buffer_corners_chunks_worker2_msg',i6.6,'.txt')") imsg
    endif

! only do this if current processor is the right one for MPI version
! this line is ok even for NCHUNKS = 2
  if(iprocscorners(imember_corner,imsg) == myrank) then

! pick the correct 1D buffer
! this scheme works fine even if NPROC_XI = NPROC_ETA = 1
  if(itypecorner(imember_corner,imsg) == ILOWERLOWER) then
    filename_in = prname(1:len_trim(prname))//'ibool1D_leftxi_lefteta.txt'
    NGLOB1D_RADIAL = NGLOB1D_RADIAL_CORNER(iregion_code,1)
  else if(itypecorner(imember_corner,imsg) == ILOWERUPPER) then
    filename_in = prname(1:len_trim(prname))//'ibool1D_leftxi_righteta.txt'
    NGLOB1D_RADIAL = NGLOB1D_RADIAL_CORNER(iregion_code,4)
  else if(itypecorner(imember_corner,imsg) == IUPPERLOWER) then
    filename_in = prname(1:len_trim(prname))//'ibool1D_rightxi_lefteta.txt'
    NGLOB1D_RADIAL = NGLOB1D_RADIAL_CORNER(iregion_code,2)
  else if(itypecorner(imember_corner,imsg) == IUPPERUPPER) then
    filename_in = prname(1:len_trim(prname))//'ibool1D_rightxi_righteta.txt'
    NGLOB1D_RADIAL = NGLOB1D_RADIAL_CORNER(iregion_code,3)
  else
    call exit_MPI(myrank,'incorrect corner coordinates')
  endif

! read 1D buffer for corner
    open(unit=IIN,file=filename_in,status='old',action='read')
    do ipoin1D = 1,NGLOB1D_RADIAL
      read(IIN,*) ibool1D(ipoin1D), &
              xread1D(ipoin1D),yread1D(ipoin1D),zread1D(ipoin1D)
    enddo
    close(IIN)

! sort array read based upon the coordinates of the points
! to ensure conforming matching with other buffers from neighbors
    call sort_array_coordinates(NGLOB1D_RADIAL,xread1D,yread1D,zread1D, &
            ibool1D,iglob,locval,ifseg,nglob,ind,ninseg,iwork,work)

! check that no duplicates have been found
    if(nglob /= NGLOB1D_RADIAL) call exit_MPI(myrank,'duplicates found for corners')

! write file with 1D buffer for corner
    open(unit=IOUT_BUFFERS,file=prname(1:len_trim(prname))//filename_out,status='unknown')
    write(IOUT_BUFFERS,*) NGLOB1D_RADIAL
    do ipoin1D = 1,NGLOB1D_RADIAL
      write(IOUT_BUFFERS,*) ibool1D(ipoin1D), &
              xread1D(ipoin1D),yread1D(ipoin1D),zread1D(ipoin1D)
    enddo
    close(IOUT_BUFFERS)

! end of section done only if right processor for MPI
  endif

  enddo

  enddo

  if(myrank == 0) close(IOUT)

! deallocate arrays
  deallocate(iproc_sender)
  deallocate(iproc_receiver)
  deallocate(npoin2D_send)
  deallocate(npoin2D_receive)

  deallocate(iprocscorners)
  deallocate(itypecorner)

  deallocate(ibool_selected)
  deallocate(xstore_selected)
  deallocate(ystore_selected)
  deallocate(zstore_selected)
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iglob)
  deallocate(locval)
  deallocate(ifseg)
  deallocate(iwork)
  deallocate(work)

  deallocate(mask_ibool)

  end subroutine create_chunk_buffers

