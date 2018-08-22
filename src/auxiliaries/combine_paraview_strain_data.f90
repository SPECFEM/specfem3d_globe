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
! the Free Software Foundation; either version 3 of the License, or
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
!
! used to visualize movie files created by MOVIE_VOLUME option in the Par_file

  use constants

  implicit none

  include "OUTPUT_FILES/values_from_mesher.h"

  integer :: fid,i,ipoint,ier,it,itstart,itstop,dit_movie
  integer :: iproc,num_node,npoint_all,nelement_all,nelement_total
  integer :: np, ne
  integer :: n1, n2, n3, n4, n5, n6, n7, n8

  integer,parameter :: MAX_NUM_NODES = 2000
  integer,dimension(MAX_NUM_NODES) :: npoint, nelement

  ! data arrays
  integer :: numpoin,nelement_local
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: xstore, ystore, zstore
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: datastore
  real(kind=CUSTOM_REAL), dimension(NGLOB_CRUST_MANTLE) :: SEEstore,SNNstore,SZZstore,SNEstore,SNZstore,SEZstore
  real(kind=CUSTOM_REAL) :: x, y, z, dat

  character(len=MAX_STRING_LEN) :: prname, dimension_file
  character(len=MAX_STRING_LEN) :: mesh_file, local_element_file, local_data_file, var_name
  character(len=MAX_STRING_LEN) :: arg(7)
  character(len=3) :: comp

  logical :: MOVIE_COARSE

  ! VTK
  logical :: USE_VTK
  character(len=3) :: vtk_format
  ! global point data
  real,dimension(:),allocatable :: total_dat
  real,dimension(:,:),allocatable :: total_dat_xyz
  integer,dimension(:,:),allocatable :: total_dat_con

  ! command-line arguments
  do i = 1,7
    call get_command_argument(i,arg(i))
    if (i < 7 .and. len_trim(arg(i)) == 0) then
      print *, ' '
      print *, ' Usage: xcombine_paraview_strain_data nnodes dt_movie itstart itstop comp MOVIE_COARSE [vtk/vtu]'
      print *, ' with'
      print *, '   nnodes              - is the number of slices'
      print *, '   dt_movie            - the time increment of movie snapshots'
      print *, '   itstart and itstop  - the time increments where to start and stop movie data files'
      print *, '   comp                - component can be SEE, SNE, SEZ, SNN, SNZ, SZZ, VEE, VEN, VEZ, I1 or I2'
      print *, '                         stored in the database directory DATABASES_MPI/'
      print *, '                         (as real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,nspec))'
      print *, '   MOVIE_COARSE        - can be either .false. or .true.'
      print *, '   vtk or vtu          - (optional) to create VTK output files specify "vtk" or "vtu"'
      print *, '                         for legacy .vtk ASCII file or newer .vtu VTK XML binary file, respectively.'
      print *, '                         Default is to output *.mesh files which need mesh2vtu utility tool'

      stop ' Reenter command line options'
    endif
  enddo

  read(arg(1),*) num_node
  read(arg(2),*) dit_movie
  read(arg(3),*) itstart
  read(arg(4),*) itstop
  read(arg(5),*) comp
  read(arg(6),*) MOVIE_COARSE
  ! vtk flag
  USE_VTK  = .false.
  if (trim(arg(7)) == 'vtk' .or. trim(arg(7)) == 'vtu') then
    USE_VTK = .true.
    vtk_format = trim(arg(7))
  endif

  ! variable name
  var_name = comp

  if (num_node > MAX_NUM_NODES) stop 'change array sizes for num_node > 1000 and recompile xcombine_paraview_movie_data'

  ! user output
  print *, 'Number of nodes: ',num_node
  print *, ' '
  print *, 'Timeframes every ',dit_movie,'from: ',itstart,' to:',itstop
  print *, ' '
  print *, 'Component  : ',trim(comp)
  print *, 'Coarse mesh: ',MOVIE_COARSE
  if (USE_VTK) then
    print *, 'Output file format : .',vtk_format
  else
    print *, 'Output file format : .mesh'
  endif
  print *, ' '

  ! figure out total number of points
  print *, 'Counting points'
  do iproc = 1, num_node
    ! print *, 'Counting elements: slice ', iproc-1
    write(prname,'(a,i6.6,a)') 'proc',iproc-1,'_'

    dimension_file = 'DATABASES_MPI/' // trim(prname) //'movie3D_info.txt'

    !   print *, 'reading: ',trim(dimension_file)
    open(unit = IIN,file = trim(dimension_file),status='old',action='read', iostat = ier)
    if (ier /= 0) then
      print *, 'Error opening file: ',trim(dimension_file)
      stop 'Error opening file'
    endif
    ! reads number of points and elements for each process slice
    read(IIN,*) npoint(iproc),nelement(iproc)
    close(IIN)
  enddo

  npoint_all   = sum(npoint(1:num_node))
  nelement_all = sum(nelement(1:num_node))
  print *, 'Total number of points   = ', npoint_all
  print *, 'Total number of elements = ', nelement_all
  print *, ' '

  if (MOVIE_COARSE) then
    nelement_total = nelement_all
  else
    nelement_total = nelement_all * 64
  endif

  ! VTK
  if (USE_VTK) then
    ! creates array to hold point data
    allocate(total_dat(npoint_all),stat=ier)
    if (ier /= 0 ) stop 'Error allocating total_dat array'
    total_dat(:) = 0.0
    allocate(total_dat_xyz(3,npoint_all),stat=ier)
    if (ier /= 0 ) stop 'Error allocating total_dat_xyz array'
    total_dat_xyz(:,:) = 0.0
    allocate(total_dat_con(8,nelement_total),stat=ier)
    if (ier /= 0 ) stop 'Error allocating total_dat_con array'
    total_dat_con(:,:) = 0
  endif

  ! loops over frames
  do it = itstart, itstop, dit_movie
    print *, '----------- Timeframe ', it, '----------------'

    ! open Paraview output mesh file
    if (USE_VTK) then
      ! VTK file format
      ! will collect all data first and then write file at the very end
      continue
    else
      ! MESH file format
      write(mesh_file,'(a,a,a,i6.6,a)')  'DATABASES_MPI/' // 'movie3D_',trim(comp),'_it',it,'.mesh'
      call open_file_fd(trim(mesh_file)//char(0),fid)
    endif

    np = 0

    ! write point and scalar information
    print *,'writing point information'
    do iproc = 1, num_node

      !print *, ' '
      !print *, 'Writing points: slice ', iproc-1,'npoints',npoint(iproc)
      write(prname,'(a,i6.6,a)') 'DATABASES_MPI/' // 'proc',iproc-1,'_'

      if (USE_VTK) then
        ! VTK nothing to write yet
        continue
      else
        ! MESH
        if (iproc == 1) then
          call write_integer_fd(fid,npoint_all)
        endif
      endif

      ! reads in grid point locations
      open(unit = IIN,file = trim(prname)//'movie3D_x.bin',status='old',action='read', iostat = ier,form ='unformatted')
      if (ier /= 0) stop 'Error opening file x.bin'
      if (npoint(iproc) > 0) then
        read(IIN) xstore(1:npoint(iproc))
      endif
      close(IIN)

      open(unit = IIN,file = trim(prname)//'movie3D_y.bin',status='old',action='read', iostat = ier,form ='unformatted')
      if (ier /= 0) stop 'Error opening file y.bin'
      if (npoint(iproc) > 0) then
        read(IIN) ystore(1:npoint(iproc))
      endif
      close(IIN)

      open(unit = IIN,file = trim(prname)//'movie3D_z.bin',status='old',action='read', iostat = ier,form ='unformatted')
      if (ier /= 0) stop 'Error opening file z.bin'
      if (npoint(iproc) > 0) then
        read(IIN) zstore(1:npoint(iproc))
      endif
      close(IIN)

      datastore(:) = 0._CUSTOM_REAL

      ! strain or invariant
      select case (trim(comp))
      case ('SI1','SI2','SI3')
        ! strain invariant options
        ! needs all strain components
        write(local_data_file,'(a,i6.6,a)') 'movie3D_SEE',it,'.bin'
        !print *, iproc,'reading from file:'//trim(prname)//trim(local_data_file)
        !print *, 'reading from file:',local_data_file
        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
          read(IIN) SEEstore(1:npoint(iproc))
        endif
        close(IIN)

        write(local_data_file,'(a,i6.6,a)') 'movie3D_SNE',it,'.bin'
        !print *, 'reading from file:',local_data_file
        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
         read(IIN) SNEstore(1:npoint(iproc))
        endif
        close(IIN)

        write(local_data_file,'(a,i6.6,a)') 'movie3D_SEZ',it,'.bin'
        !print *, 'reading from file:',local_data_file
        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
          read(IIN) SEZstore(1:npoint(iproc))
        endif
        close(IIN)

        write(local_data_file,'(a,i6.6,a)') 'movie3D_SNN',it,'.bin'
        !print *, 'reading from file:',local_data_file
        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
          read(IIN) SNNstore(1:npoint(iproc))
        endif
        close(IIN)

        write(local_data_file,'(a,i6.6,a)') 'movie3D_SNZ',it,'.bin'
        !print *, 'reading from file:',local_data_file
        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
          read(IIN) SNZstore(1:npoint(iproc))
        endif
        close(IIN)

        write(local_data_file,'(a,i6.6,a)') 'movie3D_SZZ',it,'.bin'
        !print *, 'reading from file:',local_data_file
        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
          read(IIN) SZZstore(1:npoint(iproc))
        endif
        close(IIN)

        ! strain tensor invariants
        ! see for example: http://www.continuummechanics.org/principalstrain.html
        if (trim(comp) == 'SI1') then
          ! invariant 1: strain trace (related to hydrostatic component)
          do ipoint = 1,npoint(iproc)
            ! I1 = tr(E) = E11 + E22 + E33
            dat = SEEstore(ipoint) + SNNstore(ipoint) + SZZstore(ipoint)
            datastore(ipoint) = dat
          enddo

        else if (trim(comp) == 'SI2') then
          ! invariant 2: deviatoric component
          do ipoint = 1,npoint(iproc)
            ! I2 = 1/2 [ Ekk**2 - Eij*Eij ] = E11*E22 + E11*E33 + E22*E33 - E12**2 - E13**2 - E23**2
            dat = SNNstore(ipoint)*SEEstore(ipoint) + SNNstore(ipoint)*SZZstore(ipoint) + SEEstore(ipoint)*SZZstore(ipoint) &
                  - SNEstore(ipoint)**2 - SNZstore(ipoint) - SEZstore(ipoint)**2
            datastore(ipoint) = dat
          enddo

        else
          ! invariant 3: determinant (no physical meaning)
          do ipoint = 1,npoint(iproc)
            ! I3 = det(E) = E11*E22*E33 - E11*E23**2 - E22*E13**2 - E33*E12**2 + 2 E12*E13*E23
            dat = SNNstore(ipoint)*SEEstore(ipoint)*SZZstore(ipoint) &
              - SNNstore(ipoint)*SEZstore(ipoint)**2 - SEEstore(ipoint)*SNZstore(ipoint) - SZZstore(ipoint)*SNEstore(ipoint)**2 &
              + 2.0 * SNEstore(ipoint)*SNZstore(ipoint)*SEZstore(ipoint)
            datastore(ipoint) = dat
          enddo
        endif

      case('SEE', 'SNN', 'SZZ', 'SEZ', 'SNZ', 'SNE', &
           'EEE', 'ENN', 'EZZ', 'EEZ', 'ENZ', 'ENE', &
           'PEE', 'PNN', 'PZZ', 'PEZ', 'PNZ', 'PNE', &
           'VEE', 'VEN', 'VEZ', &
           'DIE', 'DIN', 'DIZ')

        ! single component output
        !
        ! movie volume types: see write_movie_volume.f90
        ! for strain movies:
        !movie_prefix='E' ! strain
        !movie_prefix='S' ! time integral of strain
        !movie_prefix='P' ! potency, or integral of strain x \mu
        ! for vector movies:
        !movie_prefix='DI' ! displacement
        !movie_prefix='VE' ! velocity

        write(local_data_file,'(a,a,i6.6,a)') 'movie3D_',comp,it,'.bin'

        !print *,'reading comp:',trim(prname)//trim(local_data_file)

        open(unit = IIN,file = trim(prname)//trim(local_data_file),status='old',action='read', iostat = ier,form ='unformatted')
        if (ier /= 0) then
          print *,'Error opening file: ',trim(prname)//trim(local_data_file)
          stop 'Error opening file it.bin'
        endif
        if (npoint(iproc) > 0) then
          read(IIN) datastore(1:npoint(iproc))
        endif
        close(IIN)

      case default
        print *,'Error invalid component: ',trim(comp)
        stop 'Error unrecognized component'
      end select

      numpoin = 0
      do ipoint = 1,npoint(iproc)
        ! counter
        numpoin = numpoin + 1
        ! point value
        x = xstore(ipoint)
        y = ystore(ipoint)
        z = zstore(ipoint)
        dat = datastore(ipoint)
        !   print *, 'point:',ipoint,x,y,z,dat
        if (USE_VTK) then
          ! VTK
          ! stores all slices into single array
          total_dat(np+numpoin) = dat
          total_dat_xyz(1,np+numpoin) = x
          total_dat_xyz(2,np+numpoin) = y
          total_dat_xyz(3,np+numpoin) = z
        else
          ! MESH format
          call write_real_fd(fid,x)
          call write_real_fd(fid,y)
          call write_real_fd(fid,z)
          call write_real_fd(fid,dat)
        endif
      enddo

      if (numpoin /= npoint(iproc)) stop 'Error different number of points'
      np = np + npoint(iproc)

    enddo  ! all slices for points

    if (np /= npoint_all) stop 'Error: Number of total points are not consistent'
    print *, 'Total number of points: ', np
    print *, ' '

    ne = 0

    ! write element information
    print *, 'Writing element information'
    do iproc = 1, num_node

      ! print *, 'Reading slice ', iproc-1
      write(prname,'(a,i6.6,a)') 'DATABASES_MPI/' // 'proc',iproc-1,'_'

      local_element_file = trim(prname) // 'movie3D_elements.bin'
      open(unit = IIN, file = trim(local_element_file), status = 'old', action='read',iostat = ier,form='unformatted')
      if (ier /= 0) stop 'Error opening file'

      !  print *, trim(local_element_file)

      ! total point counter
      if (iproc == 1) then
        np = 0
      else
        np = sum(npoint(1:iproc-1))
      endif

      ! total number of elements
      if (USE_VTK) then
        ! for VTK, no need to output total
        continue
      else
        ! MESH format
        if (iproc == 1) then
          call write_integer_fd(fid,nelement_total)
        endif
      endif

      if (MOVIE_COARSE) then
        nelement_local = nelement(iproc)
      else
        nelement_local = nelement(iproc)*64
      endif

      ! loops over elements
      do i = 1, nelement_local
        read(IIN,iostat=ier) n1, n2, n3, n4, n5, n6, n7, n8
        if (ier /= 0) then
          print *,'Error reading file: ',trim(local_element_file)
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
        if (USE_VTK) then
          ! VTK
          ! note: indices for VTK start at 0
          total_dat_con(1,i + ne) = n1
          total_dat_con(2,i + ne) = n2
          total_dat_con(3,i + ne) = n3
          total_dat_con(4,i + ne) = n4
          total_dat_con(5,i + ne) = n5
          total_dat_con(6,i + ne) = n6
          total_dat_con(7,i + ne) = n7
          total_dat_con(8,i + ne) = n8
        else
          ! MESH format
          call write_integer_fd(fid,n1)
          call write_integer_fd(fid,n2)
          call write_integer_fd(fid,n3)
          call write_integer_fd(fid,n4)
          call write_integer_fd(fid,n5)
          call write_integer_fd(fid,n6)
          call write_integer_fd(fid,n7)
          call write_integer_fd(fid,n8)
        endif
        !write(*,*) n1, n2, n3, n4, n5, n6, n7, n8
      enddo
      close(IIN)

      ne = ne + nelement_local

    enddo ! num_node

    print *, 'Total number of elements: ', ne,' nelement_total',nelement_total
    if (ne /= nelement_total) stop 'Number of total elements are not consistent'

    ! finishes output file
    if (USE_VTK) then
      select case (vtk_format)
      case ('vtk')
        ! .vtk ascii
        write(mesh_file,'(a,a,a,i6.6,a)')  trim(OUTPUT_FILES_BASE) // '/' // 'movie3D_',trim(comp),'_it',it,'.vtk'
        call write_VTK_movie_data(nelement_total,npoint_all,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
        ! .vtk binary
        !write(mesh_file,'(a,a,a,i6.6,a)')  trim(OUTPUT_FILES_BASE) // '/' // 'movie3D_',trim(comp),'_it',it,'.vtk'
        !call write_VTK_movie_data_binary(nelement_total,npoint_all,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
      case ('vtu')
        ! .vtu ascii
        !write(mesh_file,'(a,a,a,i6.6,a)')  trim(OUTPUT_FILES_BASE) // '/' // 'movie3D_',trim(comp),'_it',it,'.vtu'
        !call write_VTU_movie_data(nelement_total,npoint_all,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
        ! .vtu binary
        write(mesh_file,'(a,a,a,i6.6,a)')  trim(OUTPUT_FILES_BASE) // '/' // 'movie3D_',trim(comp),'_it',it,'.vtu'
        call write_VTU_movie_data_binary(nelement_total,npoint_all,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
      case default
        stop 'Error invalid vtk_format'
      end select
    else
      ! MESH format
      call close_file_fd(fid)
    endif

    print *, 'Done writing '//trim(mesh_file)
    print *, ' '

  enddo ! timesteps

  print *, ' '

end program combine_paraview_movie_data

