!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  5 . 1
!          --------------------------------------------------
!
!          Main authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!             and University of Pau / CNRS / INRIA, France
! (c) Princeton University / California Institute of Technology and University of Pau / CNRS / INRIA
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

  subroutine read_arrays_buffers_mesher(iregion_code,myrank, &
                               iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
                               npoin2D_xi,npoin2D_eta, &
                               iprocfrom_faces,iprocto_faces,imsg_type, &
                               iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
                               iboolfaces,npoin2D_faces,iboolcorner, &
                               NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB2DMAX_XY,NGLOB1D_RADIAL, &
                               NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT,NPROC_XI,NPROC_ETA,LOCAL_PATH,NCHUNKS)

  implicit none

! standard include of the MPI library
  include 'mpif.h'

  include "constants.h"

  integer iregion_code,myrank,NCHUNKS,ier

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi,npoin2D_eta
  integer NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NGLOB2DMAX_XY,NGLOB1D_RADIAL
  integer NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT,NPROC_XI,NPROC_ETA

  integer npoin2D_faces(NUMFACES_SHARED)

  character(len=150) LOCAL_PATH

  integer, dimension(NGLOB2DMAX_XY,NUMFACES_SHARED) :: iboolfaces
  integer, dimension(NGLOB1D_RADIAL,NUMCORNERS_SHARED) :: iboolcorner
  integer, dimension(NGLOB2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi
  integer, dimension(NGLOB2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta

  integer, dimension(NUMMSGS_FACES) :: iprocfrom_faces,iprocto_faces,imsg_type

! allocate array for messages for corners
  integer, dimension(NCORNERSCHUNKS) :: iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners

  integer npoin2D_xi_mesher,npoin2D_eta_mesher
  integer npoin1D_corner

  integer imsg,icount_faces,icount_corners
  integer ipoin1D,ipoin2D

  double precision xdummy,ydummy,zdummy

! processor identification
  character(len=150) OUTPUT_FILES,prname,filename

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! get the base pathname for output files
  call get_value_string(OUTPUT_FILES, 'OUTPUT_FILES', 'OUTPUT_FILES')

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolleft_xi of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolleft_xi.txt', &
        status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error opening iboolleft_xi file')

  npoin2D_xi(1) = 1
 350  continue
  read(IIN,*) iboolleft_xi(npoin2D_xi(1)),xdummy,ydummy,zdummy
  if(iboolleft_xi(npoin2D_xi(1)) > 0) then
      npoin2D_xi(1) = npoin2D_xi(1) + 1
      goto 350
  endif
! subtract the line that contains the flag after the last point
  npoin2D_xi(1) = npoin2D_xi(1) - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_xi_mesher
  if(npoin2D_xi(1) > NGLOB2DMAX_XMIN_XMAX .or. npoin2D_xi(1) /= npoin2D_xi_mesher) &
      call exit_MPI(myrank,'incorrect iboolleft_xi read')
  close(IIN)

! read iboolright_xi of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolright_xi.txt', &
        status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error opening iboolright_xi file')

  npoin2D_xi(2) = 1
 360  continue
  read(IIN,*) iboolright_xi(npoin2D_xi(2)),xdummy,ydummy,zdummy
  if(iboolright_xi(npoin2D_xi(2)) > 0) then
      npoin2D_xi(2) = npoin2D_xi(2) + 1
      goto 360
  endif
! subtract the line that contains the flag after the last point
  npoin2D_xi(2) = npoin2D_xi(2) - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_xi_mesher
  if(npoin2D_xi(2) > NGLOB2DMAX_XMIN_XMAX .or. npoin2D_xi(2) /= npoin2D_xi_mesher) &
      call exit_MPI(myrank,'incorrect iboolright_xi read')
  close(IIN)

  if(myrank == 0) then
    write(IMAIN,*) '  #max of points in MPI buffers along xi npoin2D_xi = ', &
                                maxval(npoin2D_xi(:))
    write(IMAIN,*) '  #max of array elements transferred npoin2D_xi*NDIM = ', &
                                maxval(npoin2D_xi(:))*NDIM
    write(IMAIN,*)
  endif

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! read 2-D addressing for summation between slices along eta with MPI

! read iboolleft_eta of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolleft_eta.txt', &
        status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error opening iboolleft_eta file')

  npoin2D_eta(1) = 1
 370  continue
  read(IIN,*) iboolleft_eta(npoin2D_eta(1)),xdummy,ydummy,zdummy
  if(iboolleft_eta(npoin2D_eta(1)) > 0) then
      npoin2D_eta(1) = npoin2D_eta(1) + 1
      goto 370
  endif
! subtract the line that contains the flag after the last point
  npoin2D_eta(1) = npoin2D_eta(1) - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_eta_mesher
  if(npoin2D_eta(1) > NGLOB2DMAX_YMIN_YMAX .or. npoin2D_eta(1) /= npoin2D_eta_mesher) &
      call exit_MPI(myrank,'incorrect iboolleft_eta read')
  close(IIN)

! read iboolright_eta of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolright_eta.txt', &
        status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error opening iboolright_eta file')

  npoin2D_eta(2) = 1
 380  continue
  read(IIN,*) iboolright_eta(npoin2D_eta(2)),xdummy,ydummy,zdummy
  if(iboolright_eta(npoin2D_eta(2)) > 0) then
      npoin2D_eta(2) = npoin2D_eta(2) + 1
      goto 380
  endif
! subtract the line that contains the flag after the last point
  npoin2D_eta(2) = npoin2D_eta(2) - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_eta_mesher
  if(npoin2D_eta(2) > NGLOB2DMAX_YMIN_YMAX .or. npoin2D_eta(2) /= npoin2D_eta_mesher) &
      call exit_MPI(myrank,'incorrect iboolright_eta read')
  close(IIN)

  if(myrank == 0) then
    write(IMAIN,*) '  #max of points in MPI buffers along eta npoin2D_eta = ', &
                                maxval(npoin2D_eta(:))
    write(IMAIN,*) '  #max of array elements transferred npoin2D_eta*NDIM = ', &
                                maxval(npoin2D_eta(:))*NDIM
    write(IMAIN,*)
  endif


!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! read chunk messages only if more than one chunk
  if(NCHUNKS /= 1) then

! read messages to assemble between chunks with MPI

  if(myrank == 0) then

    ! file with the list of processors for each message for faces
    open(unit=IIN,file=trim(OUTPUT_FILES)//'/list_messages_faces.txt', &
          status='old',action='read',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening list_messages_faces file')

    do imsg = 1,NUMMSGS_FACES
      read(IIN,*) imsg_type(imsg),iprocfrom_faces(imsg),iprocto_faces(imsg)
      if      (iprocfrom_faces(imsg) < 0 &
          .or. iprocto_faces(imsg) < 0 &
          .or. iprocfrom_faces(imsg) > NPROCTOT-1 &
          .or. iprocto_faces(imsg) > NPROCTOT-1) &
        call exit_MPI(myrank,'incorrect chunk faces numbering')
      if (imsg_type(imsg) < 1 .or. imsg_type(imsg) > 3) &
        call exit_MPI(myrank,'incorrect message type labeling')
    enddo
    close(IIN)

    ! file with the list of processors for each message for corners
    open(unit=IIN,file=trim(OUTPUT_FILES)//'/list_messages_corners.txt', &
          status='old',action='read',iostat=ier)
    if( ier /= 0 ) call exit_MPI(myrank,'error opening list_messages_corners file')

    do imsg = 1,NCORNERSCHUNKS
      read(IIN,*) iproc_master_corners(imsg),iproc_worker1_corners(imsg), &
                            iproc_worker2_corners(imsg)
      if    (iproc_master_corners(imsg) < 0 &
        .or. iproc_worker1_corners(imsg) < 0 &
        .or. iproc_worker2_corners(imsg) < 0 &
        .or. iproc_master_corners(imsg) > NPROCTOT-1 &
        .or. iproc_worker1_corners(imsg) > NPROCTOT-1 &
        .or. iproc_worker2_corners(imsg) > NPROCTOT-1) &
        call exit_MPI(myrank,'incorrect chunk corner numbering')
    enddo
    close(IIN)

  endif

! broadcast the information read on the master to the nodes
  call MPI_BCAST(imsg_type,NUMMSGS_FACES,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iprocfrom_faces,NUMMSGS_FACES,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iprocto_faces,NUMMSGS_FACES,MPI_INTEGER,0,MPI_COMM_WORLD,ier)

  call MPI_BCAST(iproc_master_corners,NCORNERSCHUNKS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iproc_worker1_corners,NCORNERSCHUNKS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(iproc_worker2_corners,NCORNERSCHUNKS,MPI_INTEGER,0,MPI_COMM_WORLD,ier)
  if( ier /= 0 ) call exit_MPI(myrank,'error mpi broadcast')


!---- read indirect addressing for each message for faces of the chunks
!---- a given slice can belong to at most two faces
  icount_faces = 0
  do imsg = 1,NUMMSGS_FACES
    if(myrank == iprocfrom_faces(imsg) .or. myrank == iprocto_faces(imsg)) then
      icount_faces = icount_faces + 1

      if(icount_faces > NUMFACES_SHARED) then
        print*,'error ',myrank,' icount_faces: ',icount_faces,'NUMFACES_SHARED:',NUMFACES_SHARED
        print*,'iregion_code:',iregion_code
        call exit_MPI(myrank,'more than NUMFACES_SHARED faces for this slice')
      endif
      if(icount_faces > 2 .and. (NPROC_XI > 1 .or. NPROC_ETA > 1)) then
        print*,'error ',myrank,' icount_faces: ',icount_faces,'NPROC_XI:',NPROC_XI,'NPROC_ETA:',NPROC_ETA
        print*,'iregion_code:',iregion_code
        call exit_MPI(myrank,'more than two faces for this slice')
      endif

      ! read file with 2D buffer for faces
      if(myrank == iprocfrom_faces(imsg)) then
        write(filename,"('buffer_faces_chunks_sender_msg',i6.6,'.txt')") imsg
      else if(myrank == iprocto_faces(imsg)) then
        write(filename,"('buffer_faces_chunks_receiver_msg',i6.6,'.txt')") imsg
      endif

      open(unit=IIN,file=prname(1:len_trim(prname))//filename,status='old',action='read',iostat=ier)
      if( ier /= 0 ) call exit_MPI(myrank,'error opening buffer_faces file')

      read(IIN,*) npoin2D_faces(icount_faces)
      if(npoin2D_faces(icount_faces) > NGLOB2DMAX_XY) then
        print*,'error ',myrank,' npoin2D_faces: ',npoin2D_faces(icount_faces),icount_faces
        print*,'iregion_code:',iregion_code
        call exit_MPI(myrank,'incorrect nb of points in face buffer')
      endif

      do ipoin2D = 1,npoin2D_faces(icount_faces)
        read(IIN,*) iboolfaces(ipoin2D,icount_faces),xdummy,ydummy,zdummy
      enddo
      close(IIN)

    endif
  enddo


!---- read indirect addressing for each message for corners of the chunks
!---- a given slice can belong to at most one corner
  icount_corners = 0
  do imsg = 1,NCORNERSCHUNKS
    ! if only two chunks then there is no second worker
    if(myrank == iproc_master_corners(imsg) .or. &
         myrank == iproc_worker1_corners(imsg) .or. &
         (NCHUNKS /= 2 .and. myrank == iproc_worker2_corners(imsg))) then

      icount_corners = icount_corners + 1
      if(icount_corners>1 .and. (NPROC_XI > 1 .or. NPROC_ETA > 1)) then
        print*,'error ',myrank,'icount_corners:',icount_corners
        print*,'iregion_code:',iregion_code
        call exit_MPI(myrank,'more than one corner for this slice')
      endif
      if(icount_corners>4) call exit_MPI(myrank,'more than four corners for this slice')

      ! read file with 1D buffer for corner
      if(myrank == iproc_master_corners(imsg)) then
        write(filename,"('buffer_corners_chunks_master_msg',i6.6,'.txt')") imsg
      else if(myrank == iproc_worker1_corners(imsg)) then
        write(filename,"('buffer_corners_chunks_worker1_msg',i6.6,'.txt')") imsg
      else if( NCHUNKS /= 2 .and. myrank == iproc_worker2_corners(imsg)) then
        write(filename,"('buffer_corners_chunks_worker2_msg',i6.6,'.txt')") imsg
      endif

      ! matching codes
      open(unit=IIN,file=prname(1:len_trim(prname))//filename, &
            status='old',action='read',iostat=ier)
      if( ier /= 0 ) call exit_MPI(myrank,'error opening buffer_corners_chunks file')

      read(IIN,*) npoin1D_corner
      if(npoin1D_corner /= NGLOB1D_RADIAL) then
        print*,'error ',myrank,' npoin1D_corner: ',npoin1D_corner,'NGLOB1D_RADIAL:',NGLOB1D_RADIAL
        print*,'iregion_code:',iregion_code
        call exit_MPI(myrank,'incorrect nb of points in corner buffer')
      endif
      do ipoin1D = 1,npoin1D_corner
        read(IIN,*) iboolcorner(ipoin1D,icount_corners),xdummy,ydummy,zdummy
      enddo
      close(IIN)

    endif


  enddo

  endif

  end subroutine read_arrays_buffers_mesher

