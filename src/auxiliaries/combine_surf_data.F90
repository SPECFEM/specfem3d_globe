!=====================================================================
!
!                       S p e c f e m 3 D  G l o b e
!                       ----------------------------
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

program combine_surf_data

  ! combines the database files on several slices.

  use constants, only: &
    CUSTOM_REAL,MAX_STRING_LEN,IIN,NGLLX,NGLLY,NGLLZ

  use shared_parameters, only: &
    LOCAL_PATH

  use constants_solver, only: &
    NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE,NSPEC_OUTER_CORE,NSPEC_INNER_CORE,NPROCTOT_VAL

  implicit none

  character(len=MAX_STRING_LEN) :: arg(20)
  character(len=MAX_STRING_LEN) :: sline,slice_list_name
  character(len=MAX_STRING_LEN) :: filename,surfname,reg_name,belm_name,indir,outdir

  logical :: HIGH_RESOLUTION_MESH,FILE_ARRAY_IS_3D

  integer :: num_node
  integer,dimension(:),allocatable :: node_list,nspec_list,nglob_list

  character(len=MAX_STRING_LEN) :: prname,dimen_name,prname2,nspec2D_file,dimension_file
  character(len=MAX_STRING_LEN) :: ibelm_surf_file,data_file,ibool_file

  integer :: nspec2D_moho_val, nspec2D_400_val, nspec2D_670_val, nspec_surf
  integer :: npoint,nelement, npoint_total, nelement_total
  integer :: np, ne, numpoin, numelem
  integer, allocatable :: ibelm_surf(:)
  real(kind=CUSTOM_REAL), allocatable :: data_2D(:,:,:), data_3D(:,:,:,:)

  integer, dimension(:,:,:,:), allocatable :: ibool
  integer, dimension(:), allocatable :: num_ibool
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: xstore, ystore, zstore
  logical, dimension(:), allocatable :: mask_ibool

  real :: dat, x, y, z
  integer :: ispec, iglob, iglob1, iglob2, iglob3, iglob4, n1, n2, n3, n4, nex
  integer :: i,j,k,ispec_surf,ier
  integer :: it,njunk,ires,idimval,iproc,njunk1,njunk2,njunk3,inx,iny

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
  character(len=MAX_STRING_LEN) :: mesh_file,var_name
  ! global point data
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: total_dat
  real(kind=CUSTOM_REAL),dimension(:,:),allocatable :: total_dat_xyz
  integer,dimension(:,:),allocatable :: total_dat_con
#else
  !!! .mesh specific !!!!!!!!!!!
  integer :: pfd, efd
  character(len=MAX_STRING_LEN) :: command_name
  character(len=MAX_STRING_LEN) :: pt_mesh_file, mesh_file, em_mesh_file
#endif

! ------------------ program starts here -------------------

  do i = 1, 7
    call get_command_argument(i,arg(i))
    if (i < 7 .and. len_trim(arg(i)) == 0) then
      write(*,*) ' '
      write(*,*) ' Usage: xcombine_surf_data slice_list filename surfname input_dir output_dir high/low-resolution 2D/3D'
      write(*,*) ' filename.bin can be either'
      write(*,*) '   real(kind=CUSTOM_REAL) filename(NGLLX,NGLLY,NGLLZ,nspec)'
      write(*,*) '   or ---  filename(NGLLX,NGLLY,NSPEC2D) where'
      write(*,*) '   filename=moho_kernel, d400_kernel, d670_kernel, CMB_kernel, or ICB_kernel'
      write(*,*) ' possible surface names: Moho, 400, 670, CMB, ICB'
      write(*,*) ' files have been collected in input_dir, output mesh file goes to output_dir '
      write(*,*) ' give 0 for low resolution and 1 for high resolution'
      write(*,*) ' give 0 for 2D and 1 for 3D filenames'
      write(*,*) ' region does not have to be specified'
      stop ' Reenter command line options'
    endif
  enddo

  ! slices
  slice_list_name = trim(arg(1))

  ! file to collect
  filename = trim(arg(2))
  surfname = trim(arg(3))

  ! input and output dir
  indir = trim(arg(4))
  outdir = trim(arg(5))

  ! resolution
  read(arg(6),*) ires
  if (ires == 0) then
    HIGH_RESOLUTION_MESH = .false.
    inx = NGLLX-1; iny = NGLLY-1
  else
    HIGH_RESOLUTION_MESH = .true.
    inx = 1; iny = 1
  endif

  ! file dimension
  if (trim(arg(7)) == '') then
    ! default is 2D
    idimval = 0
  else
    read(arg(7),*) idimval
  endif
  if (idimval == 0) then
    FILE_ARRAY_IS_3D = .false.
  else
    FILE_ARRAY_IS_3D = .true.
  endif

  ! user output
  print *, 'xcombine_surf_data'
  print *, '  surface name       : ',trim(surfname)
  print *, '  input directory    : ',trim(indir)
  print *, '  output directory   : ',trim(outdir)
  print *
  print *, '  low/high resolution: ',ires
  print *, '  array is 3D        : ',idimval
  print *
#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
#ifdef USE_VTK_INSTEAD_OF_MESH
  print *, '  using VTK format for file output'
#endif
#ifdef USE_VTU_INSTEAD_OF_MESH
  print *, '  using VTU format for file output'
#endif
#else
  ! .mesh format
  print *, '  using .mesh format for file output'
#endif
  print *

  ! reads mesh parameters
  LOCAL_PATH = indir      ! mesh_parameters.bin file in indir/
  call read_mesh_parameters()

  ! user output
  print *,'mesh parameters (from input directory):'
  print *,'  NSPEC_CRUST_MANTLE = ',NSPEC_CRUST_MANTLE
  print *,'  NSPEC_OUTER_CORE   = ',NSPEC_OUTER_CORE
  print *,'  NSPEC_INNER_CORE   = ',NSPEC_INNER_CORE
  print *

  if (NSPEC_CRUST_MANTLE < NSPEC_OUTER_CORE .or. NSPEC_CRUST_MANTLE < NSPEC_INNER_CORE) &
    stop 'This program needs that NSPEC_CRUST_MANTLE > NSPEC_OUTER_CORE and NSPEC_INNER_CORE'

  ! get slice list
  num_node = 0
  if (trim(slice_list_name) == 'all' .or. trim(slice_list_name) == '-1') then
    ! uses all slices to combine
    num_node = NPROCTOT_VAL
    allocate(node_list(num_node),nspec_list(num_node),nglob_list(num_node), stat=ier)
    if (ier /= 0) stop 'Error allocating list arrays'
    node_list(:) = 0
    nspec_list(:) = 0
    nglob_list(:) = 0
    do i = 1,num_node
      node_list(i) = i-1  ! procs start at 0
    enddo
  else
    ! reads in slice file with specified slices
    open(unit = IIN, file = trim(slice_list_name), status = 'unknown',iostat = ier)
    if (ier /= 0) then
      print *,'no file: ',trim(slice_list_name)
      stop 'Error opening slices file'
    endif

    ! gets number of lines
    num_node = 0
    do while (1 == 1)
      read(IIN,'(a)',iostat=ier) sline
      if (ier /= 0) exit
      read(sline,*,iostat=ier) njunk
      if (ier /= 0) exit
      num_node = num_node + 1
    enddo

    allocate(node_list(num_node),nspec_list(num_node),nglob_list(num_node), stat=ier)
    if (ier /= 0) stop 'Error allocating list arrays'
    node_list(:) = 0
    nspec_list(:) = 0
    nglob_list(:) = 0

    ! reads in list
    rewind(IIN)
    num_node = 0
    do while (1 == 1)
      read(IIN,'(a)',iostat=ier) sline
      if (ier /= 0) exit
      read(sline,*,iostat=ier) njunk
      if (ier /= 0) exit
      num_node = num_node + 1
      node_list(num_node) = njunk
    enddo
    close(IIN)
  endif

  print *, 'Slice list: '
  print *, node_list(1:num_node)
  print *, ' '

  ! checks
  if (num_node < 1) stop 'Error need at least one slice for combining data arrays, please check your input arguments...'

  ! discontinuity surfaces
  if (trim(surfname) == 'Moho' .or. trim(surfname) == '400' .or. trim(surfname) == '670') then
    reg_name = 'reg1_'
    belm_name = trim(reg_name)//'boundary_disc.bin'
  else if (trim(surfname) == 'CMB') then ! assume CMB_top
    reg_name = 'reg1_'
    belm_name = trim(reg_name)//'boundary.bin' ! use reg2_ibelm for CMB_bot
  else if (trim(surfname) == 'ICB') then ! assume ICB_top
    reg_name = 'reg2_'
    belm_name = trim(reg_name)//'boundary.bin' ! use reg3_ibelm for ICB_bot
  else
    stop 'surfname type can only be: Moho, 400, 670, CMB, and ICB'
  endif

  ! figure out the total number of points/elements and allocate arrays
  write(prname,'(a,i6.6,a)') trim(indir)//'/proc',node_list(1),'_'
  nspec2D_file = trim(prname) // trim(belm_name)

  open(IIN,file=trim(nspec2D_file),status='old',form='unformatted')
  if (trim(surfname) == 'CMB' .or. trim(surfname) == 'ICB') then
    read(IIN) njunk
    read(IIN) njunk
    read(IIN) njunk
    read(IIN) njunk
    read(IIN) nspec_surf
    ! debug
    !print *, 'debug: nspec surf = ', nspec_surf
  else
    ! boundary_disc format:
    !#NSPEC2D_MOHO #NSPEC2D_400 #NSPEC2D_670
    read(IIN) nspec2D_moho_val,nspec2D_400_val,nspec2D_670_val
    ! debug
    !print *, 'debug: nspec MOHO / 400 / 670 = ', nspec2D_moho_val,nspec2D_400_val,nspec2D_670_val

    if (trim(surfname) == 'Moho') nspec_surf = nspec2D_moho_val
    if (trim(surfname) == '400') nspec_surf = nspec2D_400_val
    if (trim(surfname) == '670') nspec_surf = nspec2D_670_val
  endif
  close(IIN)

  nex = int(dsqrt(nspec_surf*1.0d0))
  ! debug
  !print *, 'debug: nspec_surf / nex = ', nspec_surf,'/',nex

  if (HIGH_RESOLUTION_MESH) then
    npoint = (nex*(NGLLX-1)+1) * (nex*(NGLLY-1)+1)
    nelement = nspec_surf  * (NGLLX-1) * (NGLLY-1)
  else
    npoint = (nex+1) * (nex+1)
    nelement = nspec_surf
  endif
  npoint_total = npoint * num_node
  nelement_total = nelement * num_node

  ! user output
  print *, 'total number of slices            = ', num_node
  print *, 'total number of spectral elements = ', nspec_surf
  print *
  print *, 'total number of points            = ', npoint_total
  print *, 'total number of elements          = ', nelement_total

  ! ========= write points and elements files ===================
  dimen_name = trim(reg_name)//'solver_data.bin'

  allocate(ibelm_surf(nspec_surf))
  ibelm_surf(:) = 0

  do it = 1, num_node
    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',node_list(it),'_'
    dimension_file = trim(prname) // trim(dimen_name)
    open(unit=IIN,file=trim(dimension_file),status='old',action='read', iostat = ier, form='unformatted')
    if (ier /= 0) stop 'Error opening file'
    read(IIN) nspec_list(it)
    read(IIN) nglob_list(it)
    close(IIN)
  enddo

  if (FILE_ARRAY_IS_3D) then
    allocate(data_3D(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE))
    data_3D(:,:,:,:) = 0.0
  else
    allocate(data_2D(NGLLX,NGLLY,nspec_surf))
    data_2D(:,:,:) = 0.0
  endif

  allocate(ibool(NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
           xstore(NGLOB_CRUST_MANTLE), &
           ystore(NGLOB_CRUST_MANTLE), &
           zstore(NGLOB_CRUST_MANTLE), &
           num_ibool(NGLOB_CRUST_MANTLE), &
           mask_ibool(NGLOB_CRUST_MANTLE),stat=ier)
  if (ier /= 0) stop 'Error allocating mesh arrays'
  ibool(:,:,:,:) = -1
  xstore(:) = 0.0_CUSTOM_REAL; ystore(:) = 0.0_CUSTOM_REAL; zstore(:) = 0.0_CUSTOM_REAL

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU
  ! arrays
  print *
  print *,'VTK initial total points: ',npoint_total
  print *,'VTK initial total elements: ',nelement_total
  print *
  ! creates array to hold point data
  allocate(total_dat(npoint_total),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat array'
  total_dat(:) = 0.0_CUSTOM_REAL
  allocate(total_dat_xyz(3,npoint_total),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_xyz array'
  total_dat_xyz(:,:) = 0.0_CUSTOM_REAL
  allocate(total_dat_con(4,nelement_total),stat=ier)
  if (ier /= 0 ) stop 'Error allocating total_dat_con array'
  total_dat_con(:,:) = 0
#else
  !!! .mesh specific !!!!!!!!!!!
  ! open Paraview output mesh file
  write(mesh_file,'(a,i1,a)') trim(outdir)//'/'//trim(filename)//'.surf'
  write(pt_mesh_file,'(a,i1,a)') trim(outdir)//'/'//trim(filename)//'_point.surf'
  write(em_mesh_file,'(a,i1,a)') trim(outdir)//'/'//trim(filename)//'_element.surf'

  command_name='rm -f '//trim(pt_mesh_file)//' '//trim(em_mesh_file)//' '//trim(mesh_file)
  call system_command(command_name)

  call open_file_fd(trim(pt_mesh_file)//char(0),pfd)
  call open_file_fd(trim(em_mesh_file)//char(0),efd)

  ! .mesh
  call write_integer_fd(pfd,npoint_total)
  call write_integer_fd(efd,nelement_total)
#endif

  np = 0
  ne = 0

  ! loop over slices
  do it = 1, num_node
    ! gets slice id (process id)
    iproc = node_list(it)

    ! output info
    print *, 'Reading slice ', iproc

    write(prname,'(a,i6.6,a)') trim(indir)//'/proc',iproc,'_'
    prname2 = trim(prname)//trim(reg_name)

    ! surface topology file
    ibelm_surf_file = trim(prname) // trim(belm_name)
    print *, trim(ibelm_surf_file)

    open(unit = IIN,file = trim(ibelm_surf_file),status='old', iostat = ier, form='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(ibelm_surf_file); stop
    endif

    if (trim(surfname) == 'Moho' .or. trim(surfname) == '400' .or. trim(surfname) == '670') then
      read(IIN) njunk1,njunk2,njunk3
      if (trim(surfname) == 'Moho') then;
        read(IIN) ibelm_surf  ! moho top
      else if (trim(surfname) == '400' .or. trim(surfname) == '670') then
        read(IIN) njunk       ! moho top
        read(IIN) njunk       ! moho bot
        if (trim(surfname) == '400') then
          read(IIN) ibelm_surf  ! 400 top
        else
          read(IIN) njunk       ! 400 top
          read(IIN) njunk       ! 400 bot
          read(IIN) ibelm_surf  ! 670 top
        endif
      endif
    else ! CMB or ICB
      read(IIN) njunk; read(IIN) njunk; read(IIN) njunk; read(IIN) njunk; read(IIN) njunk; read(IIN) njunk;
      read(IIN) njunk; read(IIN) njunk; read(IIN) njunk; read(IIN) njunk
      read(IIN) ibelm_surf
    endif
    close(IIN)

    ! datafile
    data_file = trim(prname2)//trim(filename)//'.bin'
    print *, trim(data_file)

    open(unit = IIN,file = trim(data_file),status='old', iostat = ier,form ='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(data_file); stop
    endif

    if (FILE_ARRAY_IS_3D) then
      read(IIN) data_3D(:,:,:,1:nspec_list(it))
    else
      read(IIN) data_2D(:,:,1:nspec_surf)
    endif
    close(IIN)

    ! ibool file
    ibool_file = trim(prname2) // 'solver_data.bin'
    print *, trim(ibool_file)

    open(unit = IIN,file = trim(ibool_file),status='old', iostat = ier, form='unformatted')
    if (ier /= 0) then
      print *,'Error opening ',trim(ibool_file); stop
    endif

    read(IIN) nspec_list(it)
    read(IIN) nglob_list(it)
    read(IIN) xstore(1:nglob_list(it))
    read(IIN) ystore(1:nglob_list(it))
    read(IIN) zstore(1:nglob_list(it))
    read(IIN) ibool(:,:,:,1:nspec_list(it))
    close(IIN)

    mask_ibool(:) = .false.
    num_ibool(:) = 0
    numpoin = 0
    k = 1
    do ispec_surf = 1,nspec_surf
      ispec = ibelm_surf(ispec_surf)
      do j = 1, NGLLY, iny
        do i = 1, NGLLX, inx
          iglob = ibool(i,j,k,ispec)
          if (iglob == -1 ) cycle

          if (.not. mask_ibool(iglob)) then
            numpoin = numpoin + 1
            x = xstore(iglob)
            y = ystore(iglob)
            z = zstore(iglob)

            if (FILE_ARRAY_IS_3D) then
              dat = data_3D(i,j,k,ispec)
            else
              dat = data_2D(i,j,ispec_surf)
            endif

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
            ! VTK/VTU
            total_dat(np+numpoin) = dat
            total_dat_xyz(1,np+numpoin) = x
            total_dat_xyz(2,np+numpoin) = y
            total_dat_xyz(3,np+numpoin) = z
#else
            ! .mesh
            call write_real_fd(pfd,x)
            call write_real_fd(pfd,y)
            call write_real_fd(pfd,z)
            call write_real_fd(pfd,dat)
!            call write_real_fd(pfd,real(ispec_surf))
#endif

            mask_ibool(iglob) = .true.
            num_ibool(iglob) = numpoin
          endif

        enddo ! i
      enddo ! j
    enddo !ispec_surf
    if (numpoin /= npoint) stop 'Error: number of points are not consistent'

    ! write element info
    numelem = 0
    do ispec_surf = 1, nspec_surf
      ispec = ibelm_surf(ispec_surf)
      do j = 1, NGLLY-1, iny
        do i = 1, NGLLX-1, inx
          numelem = numelem + 1 ! counts elements

          iglob1 = ibool(i,j,k,ispec)
          iglob2 = ibool(i+inx,j,k,ispec)
          iglob3 = ibool(i+inx,j+iny,k,ispec)
          iglob4 = ibool(i,j+iny,k,ispec)

          n1 = num_ibool(iglob1)+np-1
          n2 = num_ibool(iglob2)+np-1
          n3 = num_ibool(iglob3)+np-1
          n4 = num_ibool(iglob4)+np-1

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
          ! VTK/VTU
          ! note: indices for VTK start at 0
          ! quad4 element using an indexing (left,bottom),(right,bottom),(right,top),(left,top)
          total_dat_con(1,numelem + ne) = n1
          total_dat_con(2,numelem + ne) = n2
          total_dat_con(3,numelem + ne) = n4  ! switch points n3 and n4 to get correct orientation
          total_dat_con(4,numelem + ne) = n3
#else
          ! .mesh
          call write_integer_fd(efd,n1)
          call write_integer_fd(efd,n2)
          call write_integer_fd(efd,n3)
          call write_integer_fd(efd,n4)
#endif

        enddo
      enddo
    enddo
    if (numelem /= nelement) stop 'Error: number of elements are not consistent'

    np = np + numpoin
    ne = ne + numelem

  enddo  ! all slices for points

  if (np /= npoint_total) stop 'Error: Number of total points not consistent'
  if (ne /= nelement_total) stop 'Error: Number of total elements not consistent'

#if defined(USE_VTK_INSTEAD_OF_MESH) || defined(USE_VTU_INSTEAD_OF_MESH)
  ! VTK/VTU format
  var_name = trim(surfname)
#ifdef USE_VTK_INSTEAD_OF_MESH
  ! default VTK
  ! outputs unstructured grid file
  write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // trim(filename) // '.vtk'
  call write_VTK_2Dmovie_data(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
#endif
#ifdef USE_VTU_INSTEAD_OF_MESH
  ! VTU binary format
  write(mesh_file,'(a,i1,a)') trim(outdir)//'/' // trim(filename) // '.vtu'
  call write_VTU_2Dmovie_data_binary(ne,np,total_dat_xyz,total_dat_con,total_dat,mesh_file,var_name)
#endif
  print *
  print *,'written: ',trim(mesh_file)
  print *

  ! free arrays for this region
  !deallocate(total_dat,total_dat_xyz,total_dat_con)
  deallocate(total_dat)
  deallocate(total_dat_xyz)
  deallocate(total_dat_con)
#else
  ! .mesh
  call close_file_fd(pfd)
  call close_file_fd(efd)

  ! cat files
  command_name='cat '//trim(pt_mesh_file)//' '//trim(em_mesh_file)//'>'//trim(mesh_file)
  print *, ' '
  print *, 'cat mesh files ...'
  print *, trim(command_name)
  call system_command(command_name)
#endif

  print *, 'Done writing '//trim(mesh_file)
  print *, ' '

  ! free arrays
  deallocate(node_list,nspec_list,nglob_list)
  deallocate(ibool,xstore,ystore,zstore,num_ibool,mask_ibool)

end program combine_surf_data

