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

! code to check that all the internal MPI buffers are okay
! inside any given chunk, along both xi and eta
! we compare the coordinates of the points in the buffers

  program check_buffers_2D

  use constants
  use shared_parameters

  implicit none

  integer ithisproc,iotherproc

  integer ipoin

  integer npoin2d_xi_save,npoin2d_xi_mesher,npoin2d_xi
  integer npoin2d_eta_save,npoin2d_eta_mesher,npoin2d_eta

! for addressing of the slices
  integer ichunk,iproc_xi,iproc_eta,iproc
  integer iproc_read,iregion_code
  integer, dimension(:,:,:), allocatable :: addressing

  double precision diff

! 2-D addressing and buffers for summation between slices
  integer, dimension(:), allocatable :: iboolleft_xi,iboolright_xi, &
    iboolleft_eta,iboolright_eta

! coordinates of the points to compare
  double precision, dimension(:), allocatable :: xleft_xi,yleft_xi,zleft_xi, &
     xright_xi,yright_xi,zright_xi,xleft_eta,yleft_eta,zleft_eta, &
     xright_eta,yright_eta,zright_eta

! processor identification
  character(len=150) prname,prname_other

! ************** PROGRAM STARTS HERE **************

  print *
  print *,'Check all MPI buffers along xi and eta inside each chunk'
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

! dynamic memory allocation for arrays
  allocate(addressing(NCHUNKS,0:NPROC_XI-1,0:NPROC_ETA-1))

! open file with global slice number addressing
  print *,'reading slice addressing'
  open(unit=34,file=trim(OUTPUT_FILES)//'/addressing.txt',status='old',action='read')
  do iproc = 0,NPROCTOT-1
      read(34,*) iproc_read,ichunk,iproc_xi,iproc_eta
      if(iproc_read /= iproc) stop 'incorrect slice number read'
      addressing(ichunk,iproc_xi,iproc_eta) = iproc
  enddo
  close(34)

! loop over all the regions of the mesh
  do iregion_code = 1,MAX_NUM_REGIONS

  print *
  print *,' ********* checking region ',iregion_code,' *********'
  print *

! dynamic memory allocation for arrays
  allocate(iboolleft_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(iboolright_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(iboolleft_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(iboolright_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(xleft_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(yleft_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(zleft_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(xright_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(yright_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(zright_xi(NGLOB2DMAX_XMIN_XMAX(iregion_code)))
  allocate(xleft_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(yleft_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(zleft_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(xright_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(yright_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))
  allocate(zright_eta(NGLOB2DMAX_YMIN_YMAX(iregion_code)))

! loop on the chunks
  do ichunk = 1,NCHUNKS

  print *
  print *,'Checking xi in chunk ',ichunk
  print *

! double loop on NPROC_XI and NPROC_ETA
  do iproc_eta=0,NPROC_ETA-1

  print *,'checking row ',iproc_eta

  do iproc_xi=0,NPROC_XI-2

  print *,'checking slice ixi = ',iproc_xi,' in that row'

  ithisproc = addressing(ichunk,iproc_xi,iproc_eta)
  iotherproc = addressing(ichunk,iproc_xi+1,iproc_eta)

! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolright_xi of this slice
  write(*,*) 'reading MPI buffer iboolright_xi slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'iboolright_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 360  continue
  read(34,*) iboolright_xi(npoin2D_xi), &
              xright_xi(npoin2D_xi),yright_xi(npoin2D_xi),zright_xi(npoin2D_xi)
  if(iboolright_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 360
  endif
  npoin2D_xi = npoin2D_xi - 1
  write(*,*) 'found ',npoin2D_xi,' points in iboolright_xi slice ',ithisproc
  read(34,*) npoin2D_xi_mesher
  if(npoin2D_xi > NGLOB2DMAX_XMIN_XMAX(iregion_code) .or. npoin2D_xi /= npoin2D_xi_mesher) then
      stop 'incorrect iboolright_xi read'
  endif
  close(34)

! save to compare to other side
  npoin2D_xi_save = npoin2D_xi

! read iboolleft_xi of other slice
  write(*,*) 'reading MPI buffer iboolleft_xi slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'iboolleft_xi.txt',status='old',action='read')
  npoin2D_xi = 1
 350  continue
  read(34,*) iboolleft_xi(npoin2D_xi), &
              xleft_xi(npoin2D_xi),yleft_xi(npoin2D_xi),zleft_xi(npoin2D_xi)
  if(iboolleft_xi(npoin2D_xi) > 0) then
      npoin2D_xi = npoin2D_xi + 1
      goto 350
  endif
  npoin2D_xi = npoin2D_xi - 1
  write(*,*) 'found ',npoin2D_xi,' points in iboolleft_xi slice ',iotherproc
  read(34,*) npoin2D_xi_mesher
  if(npoin2D_xi > NGLOB2DMAX_XMIN_XMAX(iregion_code) .or. npoin2D_xi /= npoin2D_xi_mesher) then
      stop 'incorrect iboolleft_xi read'
  endif
  close(34)

  if(npoin2D_xi_save == npoin2D_xi) then
      print *,'okay, same size for both buffers'
  else
      stop 'wrong buffer size'
  endif

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin2D_xi
      diff = dmax1(dabs(xleft_xi(ipoin)-xright_xi(ipoin)), &
       dabs(yleft_xi(ipoin)-yright_xi(ipoin)),dabs(zleft_xi(ipoin)-zright_xi(ipoin)))
      if(diff > 0.0000001d0) print *,'different: ',ipoin,iboolleft_xi(ipoin),iboolright_xi(ipoin),diff
  enddo

  enddo
  enddo
  enddo


! loop on the chunks
  do ichunk = 1,NCHUNKS

  print *
  print *,'Checking eta in chunk ',ichunk
  print *

! double loop on NPROC_XI and NPROC_ETA
  do iproc_xi=0,NPROC_XI-1

  print *,'checking row ',iproc_xi

  do iproc_eta=0,NPROC_ETA-2

  print *,'checking slice ieta = ',iproc_eta,' in that row'

  ithisproc = addressing(ichunk,iproc_xi,iproc_eta)
  iotherproc = addressing(ichunk,iproc_xi,iproc_eta+1)

! create the name for the database of the current slide
  call create_serial_name_database(prname,ithisproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)
  call create_serial_name_database(prname_other,iotherproc,iregion_code, &
      LOCAL_PATH,NPROCTOT,OUTPUT_FILES)

! read 2-D addressing for summation between slices along xi with MPI

! read iboolright_eta of this slice
  write(*,*) 'reading MPI buffer iboolright_eta slice ',ithisproc
  open(unit=34,file=prname(1:len_trim(prname))//'iboolright_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 460  continue
  read(34,*) iboolright_eta(npoin2D_eta), &
              xright_eta(npoin2D_eta),yright_eta(npoin2D_eta),zright_eta(npoin2D_eta)
  if(iboolright_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 460
  endif
  npoin2D_eta = npoin2D_eta - 1
  write(*,*) 'found ',npoin2D_eta,' points in iboolright_eta slice ',ithisproc
  read(34,*) npoin2D_eta_mesher
  if(npoin2D_eta > NGLOB2DMAX_YMIN_YMAX(iregion_code) .or. npoin2D_eta /= npoin2D_eta_mesher) then
      stop 'incorrect iboolright_eta read'
  endif
  close(34)

! save to compare to other side
  npoin2D_eta_save = npoin2D_eta

! read iboolleft_eta of other slice
  write(*,*) 'reading MPI buffer iboolleft_eta slice ',iotherproc
  open(unit=34,file=prname_other(1:len_trim(prname_other))//'iboolleft_eta.txt',status='old',action='read')
  npoin2D_eta = 1
 450  continue
  read(34,*) iboolleft_eta(npoin2D_eta), &
              xleft_eta(npoin2D_eta),yleft_eta(npoin2D_eta),zleft_eta(npoin2D_eta)
  if(iboolleft_eta(npoin2D_eta) > 0) then
      npoin2D_eta = npoin2D_eta + 1
      goto 450
  endif
  npoin2D_eta = npoin2D_eta - 1
  write(*,*) 'found ',npoin2D_eta,' points in iboolleft_eta slice ',iotherproc
  read(34,*) npoin2D_eta_mesher
  if(npoin2D_eta > NGLOB2DMAX_YMIN_YMAX(iregion_code) .or. npoin2D_eta /= npoin2D_eta_mesher) then
      stop 'incorrect iboolleft_eta read'
  endif
  close(34)

  if(npoin2D_eta_save == npoin2D_eta) then
      print *,'okay, same size for both buffers'
  else
      stop 'wrong buffer size'
  endif

! check the coordinates of all the points in the buffer
! to see if it is correctly sorted
  do ipoin = 1,npoin2D_eta
      diff = dmax1(dabs(xleft_eta(ipoin)-xright_eta(ipoin)), &
       dabs(yleft_eta(ipoin)-yright_eta(ipoin)),dabs(zleft_eta(ipoin)-zright_eta(ipoin)))
      if(diff > 0.0000001d0) print *,'different: ',ipoin,iboolleft_eta(ipoin),iboolright_eta(ipoin),diff
  enddo

  enddo
  enddo
  enddo

! deallocate arrays
  deallocate(iboolleft_xi)
  deallocate(iboolright_xi)
  deallocate(iboolleft_eta)
  deallocate(iboolright_eta)
  deallocate(xleft_xi)
  deallocate(yleft_xi)
  deallocate(zleft_xi)
  deallocate(xright_xi)
  deallocate(yright_xi)
  deallocate(zright_xi)
  deallocate(xleft_eta)
  deallocate(yleft_eta)
  deallocate(zleft_eta)
  deallocate(xright_eta)
  deallocate(yright_eta)
  deallocate(zright_eta)

  enddo

  print *
  print *,'done'
  print *

  end program check_buffers_2D

