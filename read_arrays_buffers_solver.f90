!=====================================================================
!
!          S p e c f e m 3 D  G l o b e  V e r s i o n  3 . 3
!          --------------------------------------------------
!
!                 Dimitri Komatitsch and Jeroen Tromp
!    Seismological Laboratory - California Institute of Technology
!        (c) California Institute of Technology September 2002
!
!    A signed non-commercial agreement is required to use this program.
!   Please check http://www.gps.caltech.edu/research/jtromp for details.
!           Free for non-commercial academic research ONLY.
!      This program is distributed WITHOUT ANY WARRANTY whatsoever.
!      Do not redistribute this program without written permission.
!
!=====================================================================

  subroutine read_arrays_buffers_solver(iregion_code,myrank, &
     iboolleft_xi,iboolright_xi,iboolleft_eta,iboolright_eta, &
     icode_xi,icode_eta,npoin2D_xi,npoin2D_eta, &
     iprocfrom_faces,iprocto_faces,imsg_type, &
     iproc_master_corners,iproc_slave1_corners,iproc_slave2_corners, &
     iboolfaces,icodefaces,npoin2D_faces,iboolcorner,icodecorner, &
     NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY,NPOIN1D_RADIAL, &
     NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT,NPROC_XI,NPROC_ETA,LOCAL_PATH)

  implicit none

  include "constants.h"

  integer iregion_code,myrank

  integer npoin2D_xi,npoin2D_eta
  integer NPOIN2DMAX_XMIN_XMAX,NPOIN2DMAX_YMIN_YMAX,NPOIN2DMAX_XY,NPOIN1D_RADIAL
  integer NUMMSGS_FACES,NCORNERSCHUNKS,NPROCTOT,NPROC_XI,NPROC_ETA

  integer npoin2D_faces(NUMFACES_SHARED)

  character(len=150) LOCAL_PATH

  integer, dimension(NPOIN2DMAX_XY,NUMFACES_SHARED) :: iboolfaces,icodefaces
  integer, dimension(NPOIN1D_RADIAL,NUMCORNERS_SHARED) :: iboolcorner,icodecorner
  integer, dimension(NPOIN2DMAX_XMIN_XMAX) :: iboolleft_xi,iboolright_xi,icode_xi
  integer, dimension(NPOIN2DMAX_YMIN_YMAX) :: iboolleft_eta,iboolright_eta,icode_eta

  integer, dimension(NUMMSGS_FACES) :: iprocfrom_faces,iprocto_faces,imsg_type

! allocate array for messages for corners
  integer, dimension(NCORNERSCHUNKS) :: iproc_master_corners,iproc_slave1_corners,iproc_slave2_corners

  integer npoin2D_xi_mesher,npoin2D_eta_mesher
  integer npoin1D_corner

  integer imsg,icount_faces,icount_corners
  integer ipoin1D,ipoin2D

  integer idummy
  double precision xdummy,ydummy,zdummy

! processor identification
  character(len=150) prname,filename

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! create the name for the database of the current slide and region
  call create_name_database(prname,myrank,iregion_code,LOCAL_PATH)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolleft_xi of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolleft_xi.txt',status='old')
  npoin2D_xi = 1
 350  continue
  read(IIN,*) iboolleft_xi(npoin2D_xi),icode_xi(npoin2D_xi),xdummy,ydummy,zdummy
  if(iboolleft_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 350
  endif
! subtract the line that contains the flag after the last point
  npoin2D_xi = npoin2D_xi - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX .or. npoin2D_xi /= npoin2D_xi_mesher) &
      call exit_MPI(myrank,'incorrect iboolleft_xi read')
  close(IIN)

! read iboolright_xi of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='old')
  npoin2D_xi = 1
 360  continue
! here we do not read code because it is equal to code on the left icode_xi()
  read(IIN,*) iboolright_xi(npoin2D_xi),idummy,xdummy,ydummy,zdummy
  if(iboolright_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 360
  endif
! subtract the line that contains the flag after the last point
  npoin2D_xi = npoin2D_xi - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_xi_mesher
  if(npoin2D_xi > NPOIN2DMAX_XMIN_XMAX .or. npoin2D_xi /= npoin2D_xi_mesher) &
      call exit_MPI(myrank,'incorrect iboolright_xi read')
  close(IIN)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '# of points in MPI buffers along xi npoin2D_xi = ', &
                                npoin2D_xi
    write(IMAIN,*) '# of array elements transferred npoin2D_xi*NDIM = ', &
                                npoin2D_xi*NDIM
    write(IMAIN,*)
  endif

! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! read 2-D addressing for summation between slices along eta with MPI

! read iboolleft_eta of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolleft_eta.txt',status='old')
  npoin2D_eta = 1
 370  continue
  read(IIN,*) iboolleft_eta(npoin2D_eta),icode_eta(npoin2D_eta),xdummy,ydummy,zdummy
  if(iboolleft_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 370
  endif
! subtract the line that contains the flag after the last point
  npoin2D_eta = npoin2D_eta - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX .or. npoin2D_eta /= npoin2D_eta_mesher) &
      call exit_MPI(myrank,'incorrect iboolleft_eta read')
  close(IIN)

! read iboolright_eta of this slice
  open(unit=IIN,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='old')
  npoin2D_eta = 1
 380  continue
! here we do not read code because it is equal to code on the left icode_eta()
  read(IIN,*) iboolright_eta(npoin2D_eta),idummy,xdummy,ydummy,zdummy
  if(iboolright_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 380
  endif
! subtract the line that contains the flag after the last point
  npoin2D_eta = npoin2D_eta - 1
! read nb of points given by the mesher
  read(IIN,*) npoin2D_eta_mesher
  if(npoin2D_eta > NPOIN2DMAX_YMIN_YMAX .or. npoin2D_eta /= npoin2D_eta_mesher) &
      call exit_MPI(myrank,'incorrect iboolright_eta read')
  close(IIN)

  if(myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '# of points in MPI buffers along eta npoin2D_eta = ', &
                                npoin2D_eta
    write(IMAIN,*) '# of array elements transferred npoin2D_eta*NDIM = ', &
                                npoin2D_eta*NDIM
    write(IMAIN,*)
  endif


!! $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! read chunk messages only if more than one chunk
  if(NCHUNKS /= 1) then

! read messages to assemble between chunks with MPI

! file with the list of processors for each message for faces
  open(unit=IIN,file='OUTPUT_FILES/list_messages_faces.txt',status='old')
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
  open(unit=IIN,file='OUTPUT_FILES/list_messages_corners.txt',status='old')
  do imsg = 1,NCORNERSCHUNKS
  read(IIN,*) iproc_master_corners(imsg),iproc_slave1_corners(imsg), &
                          iproc_slave2_corners(imsg)
  if    (iproc_master_corners(imsg) < 0 &
    .or. iproc_slave1_corners(imsg) < 0 &
    .or. iproc_slave2_corners(imsg) < 0 &
    .or. iproc_master_corners(imsg) > NPROCTOT-1 &
    .or. iproc_slave1_corners(imsg) > NPROCTOT-1 &
    .or. iproc_slave2_corners(imsg) > NPROCTOT-1) &
      call exit_MPI(myrank,'incorrect chunk corner numbering')
  enddo
  close(IIN)

!---- read indirect addressing for each message for faces of the chunks
!---- a given slice can belong to at most two faces
  icount_faces = 0
  do imsg = 1,NUMMSGS_FACES
  if(myrank == iprocfrom_faces(imsg) .or. myrank == iprocto_faces(imsg)) then
    icount_faces = icount_faces + 1
    if(icount_faces>NUMFACES_SHARED) call exit_MPI(myrank,'more than NUMFACES_SHARED faces for this slice')
    if(icount_faces>2 .and. (NPROC_XI > 1 .or. NPROC_ETA > 1)) call exit_MPI(myrank,'more than two faces for this slice')

! read file with 2D buffer for faces
    if(myrank == iprocfrom_faces(imsg)) then
      write(filename,500) imsg
    else if(myrank == iprocto_faces(imsg)) then
      write(filename,510) imsg
    endif

    open(unit=IIN,file=prname(1:len_trim(prname))//filename,status='old')
    read(IIN,*) npoin2D_faces(icount_faces)
    if(npoin2D_faces(icount_faces) > NPOIN2DMAX_XY) &
      call exit_MPI(myrank,'incorrect nb of points in face buffer')
    do ipoin2D = 1,npoin2D_faces(icount_faces)
! matching codes
      read(IIN,*) iboolfaces(ipoin2D,icount_faces),icodefaces(ipoin2D,icount_faces),xdummy,ydummy,zdummy
    enddo
    close(IIN)
  endif
  enddo


!---- read indirect addressing for each message for corners of the chunks
!---- a given slice can belong to at most one corner
  icount_corners = 0
  do imsg = 1,NCORNERSCHUNKS
  if(myrank == iproc_master_corners(imsg) .or. &
       myrank == iproc_slave1_corners(imsg) .or. &
       myrank == iproc_slave2_corners(imsg)) then
    icount_corners = icount_corners + 1
    if(icount_corners>1 .and. (NPROC_XI > 1 .or. NPROC_ETA > 1)) &
      call exit_MPI(myrank,'more than one corner for this slice')
    if(icount_corners>4) call exit_MPI(myrank,'more than four corners for this slice')

! read file with 1D buffer for corner
    if(myrank == iproc_master_corners(imsg)) then
      write(filename,600) imsg
    else if(myrank == iproc_slave1_corners(imsg)) then
      write(filename,610) imsg
    else if(myrank == iproc_slave2_corners(imsg)) then
      write(filename,620) imsg
    endif

! matching codes
    open(unit=IIN,file=prname(1:len_trim(prname))//filename,status='old')
    read(IIN,*) npoin1D_corner
    if(npoin1D_corner /= NPOIN1D_RADIAL) &
      call exit_MPI(myrank,'incorrect nb of points in corner buffer')
    do ipoin1D = 1,npoin1D_corner
      read(IIN,*) iboolcorner(ipoin1D,icount_corners),icodecorner(ipoin1D,icount_corners),xdummy,ydummy,zdummy
    enddo
    close(IIN)
  endif
  enddo

  endif

  500 format('buffer_faces_chunks_sender_msg',i4.4,'.txt')
  510 format('buffer_faces_chunks_receiver_msg',i4.4,'.txt')

  600 format('buffer_corners_chunks_master_msg',i4.4,'.txt')
  610 format('buffer_corners_chunks_slave1_msg',i4.4,'.txt')
  620 format('buffer_corners_chunks_slave2_msg',i4.4,'.txt')

  end subroutine read_arrays_buffers_solver

