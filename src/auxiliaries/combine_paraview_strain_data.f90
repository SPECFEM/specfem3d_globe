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

program combine_paraview_movie_data

! combines the database files on several slices.
! the local database file needs to have been collected onto the frontend (copy_local_database.pl)

  use constants

  implicit none

  include "OUTPUT_FILES/values_from_mesher.h"

  integer fid,i,ipoint, ios, it,itstart,itstop,dit_movie
  integer iproc, num_node,  npoint_all, nelement_all
  integer np, ne, npoint(1000), nelement(1000), n1, n2, n3, n4, n5, n6, n7, n8

  integer numpoin,nelement_local
!  real(kind=CUSTOM_REAL),dimension(NGLOBMAX_CRUST_MANTLE) :: xstore, ystore, zstore,datstore
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore, ystore, zstore,datstore
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: SEEstore,SNNstore,SZZstore,SNEstore,SNZstore,SEZstore
  real(kind=CUSTOM_REAL) :: x, y, z, dat
  character(len=MAX_STRING_LEN) :: arg(7), prname, dimension_file
  character(len=MAX_STRING_LEN) :: mesh_file, local_element_file, local_data_file
  character(len=3) :: comp
  logical :: MOVIE_COARSE

  do i = 1,6
    call get_command_argument(i,arg(i))
    if (i < 7 .and. len_trim(arg(i)) == 0) then
      print *, ' '
      print *, ' Usage: xcombine_data nnodes dt_movie itstart itstop comp MOVIE_COARSE'
      print *, '   nnodes is the number of slices, dt_movie the time increment of movie snapshots'
      print *, '   itstart and itstop the time increments where to start and stop movie data files'
      print *, '   component can be SEE, SNE,SEZ,SNN,SNZ,SZZ,I1 or I2'
      print *, '   stored in the local directory as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,nspec)  '
      print *, '   MOVIE_COARSE can be either .false. or .true. '
      stop ' Reenter command line options'
    endif
  enddo


  read(arg(1),*) num_node
  read(arg(2),*) dit_movie
  read(arg(3),*) itstart
  read(arg(4),*) itstop
  read(arg(5),*) comp
  read(arg(6),*) MOVIE_COARSE

  if (num_node>1000) stop 'change array sizes for num_node > 1000 and recompile xcombine_paraview_movie_data'

  print *, 'Number of nodes: ',num_node
  print *, ' '
  print *, 'Timeframes every ',dit_movie,'from: ',itstart,' to:',itstop

  ! figure out total number of points
  print *, 'Counting points'
  do iproc = 1, num_node
    ! print *, 'Counting elements: slice ', iproc-1
    write(prname,'(a,i6.6,a)') 'proc',iproc-1,'_'

    dimension_file = 'DATABASES_MPI/' // trim(prname) //'movie3D_info.txt'

    !   print *, 'reading: ',trim(dimension_file)
    open(unit = IIN,file = trim(dimension_file),status='old',action='read', iostat = ios)
    if (ios /= 0) then
      print*, 'Error opening file: ',trim(dimension_file)
      stop 'Error opening file'
    endif
    read(IIN,*) npoint(iproc),nelement(iproc)
    close(IIN)
  enddo

  npoint_all   = sum(npoint(1:num_node))
  nelement_all = sum(nelement(1:num_node))
  print *, 'Total number of points   = ', npoint_all
  print *, 'Total number of elements = ', nelement_all


  do it = itstart, itstop, dit_movie
    print *, '----------- Timeframe ', it, '----------------'

    ! open Paraview output mesh file
    write(mesh_file,'(a,a,a,i6.6,a)')  'DATABASES_MPI/' // 'movie3D_',trim(comp),'_it',it,'.mesh'
    call open_file_fd(trim(mesh_file)//char(0),fid)

    np = 0

    ! write point and scalar information
    print *,'writing point information'
    do iproc = 1, num_node


    !print *, ' '
    !print *, 'Writing points: slice ', iproc-1,'npoints',npoint(iproc)
    write(prname,'(a,i6.6,a)') 'DATABASES_MPI/' // 'proc',iproc-1,'_'

    numpoin = 0


    if (iproc == 1) then
      call write_integer_fd(fid,npoint_all)
    endif

    open(unit = IIN,file = trim(prname)//'movie3D_x.bin',status='old',action='read', iostat = ios,form ='unformatted')
    if (ios /= 0) stop 'Error opening file x.bin'
    if (npoint(iproc)>0) then
      read(IIN) xstore(1:npoint(iproc))
    endif
    close(IIN)

    open(unit = IIN,file = trim(prname)//'movie3D_y.bin',status='old',action='read', iostat = ios,form ='unformatted')
    if (ios /= 0) stop 'Error opening file y.bin'
    if (npoint(iproc)>0) then
      read(IIN) ystore(1:npoint(iproc))
    endif
    close(IIN)

    open(unit = IIN,file = trim(prname)//'movie3D_z.bin',status='old',action='read', iostat = ios,form ='unformatted')
    if (ios /= 0) stop 'Error opening file z.bin'
    if (npoint(iproc)>0) then
      read(IIN) zstore(1:npoint(iproc))
    endif
    close(IIN)

    if ((comp /= 'SI1') .and. (comp /= 'SI2')) then
!comp == 'SEE' .or. comp == 'SNN' .or. comp == 'SZZ' .or. comp == 'SEZ' .or. comp == 'SNZ' .or. comp == 'SNE') then
      write(local_data_file,'(a,a,i6.6,a)') 'movie3D_',comp,it,'.bin'

      !print *,'reading comp:',trim(prname)//trim(local_data_file)

      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
        read(IIN) datstore(1:npoint(iproc))
      endif
      close(IIN)

    else if (comp == 'SI1' .or. comp == 'SI2') then
      write(local_data_file,'(a,i6.6,a)') 'movie3D_SEE',it,'.bin'
      !print *, iproc,'reading from file:'//trim(prname)//trim(local_data_file)
      !print *, 'reading from file:',local_data_file
      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
        read(IIN) SEEstore(1:npoint(iproc))
      endif
      close(IIN)

      write(local_data_file,'(a,i6.6,a)') 'movie3D_SNE',it,'.bin'
      !print *, 'reading from file:',local_data_file
      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
       read(IIN) SNEstore(1:npoint(iproc))
      endif
      close(IIN)

      write(local_data_file,'(a,i6.6,a)') 'movie3D_SEZ',it,'.bin'
      !print *, 'reading from file:',local_data_file
      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
        read(IIN) SEZstore(1:npoint(iproc))
      endif
      close(IIN)

      write(local_data_file,'(a,i6.6,a)') 'movie3D_SNN',it,'.bin'
      !print *, 'reading from file:',local_data_file
      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
        read(IIN) SNNstore(1:npoint(iproc))
      endif
      close(IIN)

      write(local_data_file,'(a,i6.6,a)') 'movie3D_SNZ',it,'.bin'
      !print *, 'reading from file:',local_data_file
      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
        read(IIN) SNZstore(1:npoint(iproc))
      endif
      close(IIN)

      write(local_data_file,'(a,i6.6,a)') 'movie3D_SZZ',it,'.bin'
      !print *, 'reading from file:',local_data_file
      open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ios,form ='unformatted')
      if (ios /= 0) then
        print*,'Error opening file: ',trim(prname)//trim(local_data_file)
        stop 'Error opening file it.bin'
      endif
      if (npoint(iproc)>0) then
        read(IIN) SZZstore(1:npoint(iproc))
      endif
      close(IIN)
    else
      stop 'Error unrecognized component'
    endif !strain or invariant

    datstore=datstore
    do ipoint = 1,npoint(iproc)
       numpoin = numpoin + 1
       x = xstore(ipoint)
       y = ystore(ipoint)
       z = zstore(ipoint)
       dat = datstore(ipoint)
       call write_real_fd(fid,x)
       call write_real_fd(fid,y)
       call write_real_fd(fid,z)
       call write_real_fd(fid,dat)
        !   print *, 'point:',ipoint,x,y,z,dat
    enddo

    if (numpoin /= npoint(iproc)) stop 'Error different number of points'
    np = np + npoint(iproc)

  enddo  ! all slices for points

  if (np /=  npoint_all) stop 'Error: Number of total points are not consistent'
  print *, 'Total number of points: ', np
  print *, ' '

  ne = 0
  ! write element information
  print *, 'Writing element information'
  do iproc = 1, num_node

    ! print *, 'Reading slice ', iproc-1
    write(prname,'(a,i6.6,a)') 'DATABASES_MPI/' // 'proc',iproc-1,'_'

    if (iproc == 1) then
      np = 0
    else
      np = sum(npoint(1:iproc-1))
    endif


    local_element_file = trim(prname) // 'movie3D_elements.bin'
    open(unit = IIN, file = trim(local_element_file), status = 'old', action='read',iostat = ios,form='unformatted')
    if (ios /= 0) stop 'Error opening file'

    !  print *, trim(local_element_file)

    if (iproc == 1) then
      if (MOVIE_COARSE) then
        call write_integer_fd(fid,nelement_all)
      else
        call write_integer_fd(fid,nelement_all*64)
      endif
    endif

    if (MOVIE_COARSE) then
      nelement_local = nelement(iproc)
    else
      nelement_local = nelement(iproc)*64
    endif
    do i = 1, nelement_local
      read(IIN,iostat=ios) n1, n2, n3, n4, n5, n6, n7, n8
      if (ios /= 0) then
        print*,'Error reading file: ',trim(local_element_file)
        stop 'Error reading file movie3D_elements.bin, please check if number of elements is correct (MOVIE_COARSE?)'
      endif
      n1 = n1+np
      n2 = n2+np
      n3 = n3+np
      n4 = n4+np
      n5 = n5+np
      n6 = n6+np
      n7 = n7+np
      n8 = n8+np
      call write_integer_fd(fid,n1)
      call write_integer_fd(fid,n2)
      call write_integer_fd(fid,n3)
      call write_integer_fd(fid,n4)
      call write_integer_fd(fid,n5)
      call write_integer_fd(fid,n6)
      call write_integer_fd(fid,n7)
      call write_integer_fd(fid,n8)
      !write(*,*) n1, n2, n3, n4, n5, n6, n7, n8
    enddo
    close(IIN)

    ne = ne + nelement(iproc)

  enddo ! num_node
  print *, 'Total number of elements: ', ne,' nelement_all',nelement_all
  if (ne /= nelement_all) stop 'Number of total elements are not consistent'

  call close_file_fd(fid)

  print *, 'Done writing '//trim(mesh_file)
  print *, ' '

  enddo ! timesteps
  print *, ' '

end program combine_paraview_movie_data

