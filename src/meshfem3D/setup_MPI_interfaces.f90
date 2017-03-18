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


  subroutine setup_MPI_interfaces(iregion_code)

  use meshfem3D_par, only: &
    INCLUDE_CENTRAL_CUBE,myrank,NUMFACES_SHARED

  use create_MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  implicit none

  integer,intent(in):: iregion_code

  ! local parameters
  ! assigns initial maximum arrays
  ! for global slices, maximum number of neighbor is around 17 ( 8 horizontal, max of 8 on bottom )
  integer :: MAX_NEIGHBOURS,max_nibool
  integer, dimension(:),allocatable :: my_neighbours,nibool_neighbours
  integer, dimension(:,:),allocatable :: ibool_neighbours
  integer :: ier

  ! allocates temporary arrays for setup routines
  ! estimates a maximum size of needed arrays
  MAX_NEIGHBOURS = 8 + NCORNERSCHUNKS
  if (INCLUDE_CENTRAL_CUBE ) MAX_NEIGHBOURS = MAX_NEIGHBOURS + NUMMSGS_FACES

  allocate(my_neighbours(MAX_NEIGHBOURS), &
          nibool_neighbours(MAX_NEIGHBOURS),stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating my_neighbours array')
  my_neighbours(:) = -1
  nibool_neighbours(:) = 0

  ! estimates initial maximum ibool array
  max_nibool = npoin2D_max_all_CM_IC * NUMFACES_SHARED &
               + non_zero_nb_msgs_theor_in_cube*npoin2D_cube_from_slices

  allocate(ibool_neighbours(max_nibool,MAX_NEIGHBOURS), stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating ibool_neighbours')
  ibool_neighbours(:,:) = 0

  ! sets up MPI interfaces between different processes
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust/mantle
    call setup_MPI_interfaces_cm(MAX_NEIGHBOURS,my_neighbours,nibool_neighbours, &
                                max_nibool,ibool_neighbours)

  case (IREGION_OUTER_CORE)
    ! outer core
    call setup_MPI_interfaces_oc(MAX_NEIGHBOURS,my_neighbours,nibool_neighbours, &
                                max_nibool,ibool_neighbours)

  case (IREGION_INNER_CORE)
    ! inner core
    call setup_MPI_interfaces_ic(MAX_NEIGHBOURS,my_neighbours,nibool_neighbours, &
                                max_nibool,ibool_neighbours)
  end select

  ! frees temporary array
  deallocate(ibool_neighbours)
  deallocate(my_neighbours,nibool_neighbours)

  ! frees arrays not needed any further
  deallocate(iprocfrom_faces,iprocto_faces,imsg_type)
  deallocate(iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners)
  deallocate(buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar)
  deallocate(buffer_send_chunkcorn_vector,buffer_recv_chunkcorn_vector)

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust mantle
    deallocate(iboolcorner_crust_mantle)
    deallocate(iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle)
    deallocate(iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle)
    deallocate(iboolfaces_crust_mantle)
  case (IREGION_OUTER_CORE)
    ! outer core
    deallocate(iboolcorner_outer_core)
    deallocate(iboolleft_xi_outer_core,iboolright_xi_outer_core)
    deallocate(iboolleft_eta_outer_core,iboolright_eta_outer_core)
    deallocate(iboolfaces_outer_core)
  case (IREGION_INNER_CORE)
    ! inner core
    deallocate(iboolcorner_inner_core)
    deallocate(iboolleft_xi_inner_core,iboolright_xi_inner_core)
    deallocate(iboolleft_eta_inner_core,iboolright_eta_inner_core)
    deallocate(iboolfaces_inner_core)
  end select

  ! synchronizes MPI processes
  call synchronize_all()

  end subroutine setup_MPI_interfaces

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_MPI_interfaces_cm(MAX_NEIGHBOURS,my_neighbours,nibool_neighbours, &
                                    max_nibool,ibool_neighbours)

  use meshfem3D_par, only: &
    myrank,iproc_xi,iproc_eta,ichunk,addressing,INCLUDE_CENTRAL_CUBE, &
    NPROC_XI,NPROC_ETA,NPROCTOT, &
    NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NCHUNKS, &
    OUTPUT_FILES,MAX_STRING_LEN

  use meshfem3D_par, only: ibool,is_on_a_slice_edge,xstore_glob,ystore_glob,zstore_glob

  use create_MPI_interfaces_par
  use MPI_crust_mantle_par
  implicit none

  integer,intent(in) :: MAX_NEIGHBOURS,max_nibool
  integer, dimension(MAX_NEIGHBOURS),intent(inout) :: my_neighbours,nibool_neighbours
  integer, dimension(max_nibool,MAX_NEIGHBOURS),intent(inout) :: ibool_neighbours

  ! local parameters
  ! temporary buffers for send and receive between faces of the slices and the chunks
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC) :: &
    buffer_send_faces_scalar,buffer_received_faces_scalar
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: test_flag
  integer,dimension(:),allocatable :: dummy_i
  integer :: i,ier
  !----------------------
  ! debug file output
  logical,parameter :: DEBUG = .false.
  !----------------------
  character(len=MAX_STRING_LEN) :: filename

  ! sets up MPI interfaces
  ! crust mantle region
  if (myrank == 0 ) write(IMAIN,*) 'crust mantle MPI:'

  if (NPROCTOT > 1) then
    allocate(test_flag(NGLOB_CRUST_MANTLE), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating test_flag')

    ! sets flag to rank id (+1 to avoid problems with zero rank)
    test_flag(:) = myrank + 1.0

    ! assembles values
    call assemble_MPI_scalar_block(test_flag,NGLOB_CRUST_MANTLE, &
              iproc_xi,iproc_eta,ichunk,addressing, &
              iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle,iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
              npoin2D_faces_crust_mantle,npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
              iboolfaces_crust_mantle,iboolcorner_crust_mantle, &
              iprocfrom_faces,iprocto_faces,imsg_type, &
              iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
              buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
              buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
              NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
              NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
              NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE),NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
              NGLOB2DMAX_XY,NCHUNKS)

    ! removes own myrank id (+1)
    test_flag(:) = test_flag(:) - ( myrank + 1.0)

    allocate(dummy_i(NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy_i')
    dummy_i(:) = 0

    ! determines neighbor rank for shared faces
    call get_MPI_interfaces(myrank,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                              test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                              num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                              max_nibool,MAX_NEIGHBOURS, &
                              ibool,is_on_a_slice_edge, &
                              IREGION_CRUST_MANTLE,.false.,dummy_i,INCLUDE_CENTRAL_CUBE, &
                              xstore_glob,ystore_glob,zstore_glob,NPROCTOT)

    deallocate(test_flag)
    deallocate(dummy_i)
  else
    ! no interfaces
    num_interfaces_crust_mantle = 0
  endif

  ! stores MPI interfaces information
  allocate(my_neighbours_crust_mantle(num_interfaces_crust_mantle), &
           nibool_interfaces_crust_mantle(num_interfaces_crust_mantle), &
           stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbours_crust_mantle etc.')
  my_neighbours_crust_mantle(:) = -1
  nibool_interfaces_crust_mantle(:) = 0

  ! copies interfaces arrays
  if (num_interfaces_crust_mantle > 0) then
    allocate(ibool_interfaces_crust_mantle(max_nibool_interfaces_cm,num_interfaces_crust_mantle), &
             stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_crust_mantle')
    ibool_interfaces_crust_mantle(:,:) = 0

    ! ranks of neighbour processes
    my_neighbours_crust_mantle(:) = my_neighbours(1:num_interfaces_crust_mantle)
    ! number of global ibool entries on each interface
    nibool_interfaces_crust_mantle(:) = nibool_neighbours(1:num_interfaces_crust_mantle)
    ! global iglob point ids on each interface
    ibool_interfaces_crust_mantle(:,:) = ibool_neighbours(1:max_nibool_interfaces_cm,1:num_interfaces_crust_mantle)
  else
    ! dummy allocation (fortran90 should allow allocate statement with zero array size)
    max_nibool_interfaces_cm = 0
    allocate(ibool_interfaces_crust_mantle(0,0),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy array ibool_interfaces_crust_mantle')
  endif

  ! debug: outputs MPI interface
  if (DEBUG) then
    do i = 1,num_interfaces_crust_mantle
      write(filename,'(a,i6.6,a,i2.2)') trim(OUTPUT_FILES)//'/MPI_points_crust_mantle_proc',myrank, &
                      '_',my_neighbours_crust_mantle(i)
      call write_VTK_data_points(NGLOB_crust_mantle, &
                        xstore_glob,ystore_glob,zstore_glob, &
                        ibool_interfaces_crust_mantle(1:nibool_interfaces_crust_mantle(i),i), &
                        nibool_interfaces_crust_mantle(i),filename)
    enddo
    call synchronize_all()
  endif

  ! checks addressing
  call test_MPI_neighbours(IREGION_CRUST_MANTLE, &
                           num_interfaces_crust_mantle,max_nibool_interfaces_cm, &
                           my_neighbours_crust_mantle,nibool_interfaces_crust_mantle, &
                           ibool_interfaces_crust_mantle)

  ! checks with assembly of test fields
  call test_MPI_cm()

  end subroutine setup_MPI_interfaces_cm

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_MPI_interfaces_oc(MAX_NEIGHBOURS,my_neighbours,nibool_neighbours, &
                                    max_nibool,ibool_neighbours)

  use meshfem3D_par, only: &
    myrank,iproc_xi,iproc_eta,ichunk,addressing,INCLUDE_CENTRAL_CUBE, &
    NPROC_XI,NPROC_ETA,NPROCTOT, &
    NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NCHUNKS, &
    OUTPUT_FILES,MAX_STRING_LEN

  use meshfem3D_par, only: ibool,is_on_a_slice_edge,xstore_glob,ystore_glob,zstore_glob

  use create_MPI_interfaces_par
  use MPI_outer_core_par
  implicit none

  integer :: MAX_NEIGHBOURS,max_nibool
  integer, dimension(MAX_NEIGHBOURS) :: my_neighbours,nibool_neighbours
  integer, dimension(max_nibool,MAX_NEIGHBOURS) :: ibool_neighbours

  ! local parameters
  ! temporary buffers for send and receive between faces of the slices and the chunks
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC) :: &
    buffer_send_faces_scalar,buffer_received_faces_scalar
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: test_flag
  integer,dimension(:),allocatable :: dummy_i
  integer :: i,ier
  !----------------------
  ! debug file output
  logical,parameter :: DEBUG = .false.
  !----------------------
  character(len=MAX_STRING_LEN) :: filename

  ! sets up MPI interfaces
  ! outer core region
  if (myrank == 0 ) write(IMAIN,*) 'outer core MPI:'

  if (NPROCTOT > 1) then
    allocate(test_flag(NGLOB_OUTER_CORE), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating test_flag outer core')

    ! sets flag to rank id (+1 to avoid problems with zero rank)
    test_flag(:) = myrank + 1.0

    ! assembles values
    call assemble_MPI_scalar_block(test_flag,NGLOB_OUTER_CORE, &
                                   iproc_xi,iproc_eta,ichunk,addressing, &
                                   iboolleft_xi_outer_core,iboolright_xi_outer_core, &
                                   iboolleft_eta_outer_core,iboolright_eta_outer_core, &
                                   npoin2D_faces_outer_core,npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
                                   iboolfaces_outer_core,iboolcorner_outer_core, &
                                   iprocfrom_faces,iprocto_faces,imsg_type, &
                                   iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
                                   buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
                                   buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
                                   NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
                                   NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
                                   NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
                                   NGLOB2DMAX_XY,NCHUNKS)


    ! removes own myrank id (+1)
    test_flag(:) = test_flag(:) - ( myrank + 1.0)

    allocate(dummy_i(NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating dummy_i')
    dummy_i(:) = 0

    ! determines neighbor rank for shared faces
    call get_MPI_interfaces(myrank,NGLOB_OUTER_CORE,NSPEC_OUTER_CORE, &
                            test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                            num_interfaces_outer_core,max_nibool_interfaces_oc, &
                            max_nibool,MAX_NEIGHBOURS, &
                            ibool,is_on_a_slice_edge, &
                            IREGION_OUTER_CORE,.false.,dummy_i,INCLUDE_CENTRAL_CUBE, &
                            xstore_glob,ystore_glob,zstore_glob,NPROCTOT)

    deallocate(test_flag)
    deallocate(dummy_i)
  else
    ! no interfaces
    num_interfaces_outer_core = 0
  endif

  ! stores MPI interfaces information
  allocate(my_neighbours_outer_core(num_interfaces_outer_core), &
          nibool_interfaces_outer_core(num_interfaces_outer_core), &
          stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbours_outer_core etc.')
  my_neighbours_outer_core(:) = -1
  nibool_interfaces_outer_core(:) = 0

  ! copies interfaces arrays
  if (num_interfaces_outer_core > 0) then
    allocate(ibool_interfaces_outer_core(max_nibool_interfaces_oc,num_interfaces_outer_core), &
           stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_outer_core')
    ibool_interfaces_outer_core(:,:) = 0

    ! ranks of neighbour processes
    my_neighbours_outer_core(:) = my_neighbours(1:num_interfaces_outer_core)
    ! number of global ibool entries on each interface
    nibool_interfaces_outer_core(:) = nibool_neighbours(1:num_interfaces_outer_core)
    ! global iglob point ids on each interface
    ibool_interfaces_outer_core(:,:) = ibool_neighbours(1:max_nibool_interfaces_oc,1:num_interfaces_outer_core)
  else
    ! dummy allocation (fortran90 should allow allocate statement with zero array size)
    max_nibool_interfaces_oc = 0
    allocate(ibool_interfaces_outer_core(0,0),stat=ier)
  endif

  ! debug: outputs MPI interface
  if (DEBUG) then
    do i = 1,num_interfaces_outer_core
      write(filename,'(a,i6.6,a,i2.2)') trim(OUTPUT_FILES)//'/MPI_points_outer_core_proc',myrank, &
                      '_',my_neighbours_outer_core(i)
      call write_VTK_data_points(NGLOB_OUTER_CORE, &
                        xstore_glob,ystore_glob,zstore_glob, &
                        ibool_interfaces_outer_core(1:nibool_interfaces_outer_core(i),i), &
                        nibool_interfaces_outer_core(i),filename)
    enddo
    call synchronize_all()
  endif

  ! checks addressing
  call test_MPI_neighbours(IREGION_OUTER_CORE, &
                              num_interfaces_outer_core,max_nibool_interfaces_oc, &
                              my_neighbours_outer_core,nibool_interfaces_outer_core, &
                              ibool_interfaces_outer_core)

  ! checks with assembly of test fields
  call test_MPI_oc()

  end subroutine setup_MPI_interfaces_oc

!
!-------------------------------------------------------------------------------------------------
!

  subroutine setup_MPI_interfaces_ic(MAX_NEIGHBOURS,my_neighbours,nibool_neighbours, &
                                    max_nibool,ibool_neighbours)

  use meshfem3D_par, only: &
    myrank,iproc_xi,iproc_eta,ichunk,addressing,INCLUDE_CENTRAL_CUBE, &
    NPROC_XI,NPROC_ETA,NPROCTOT, &
    NGLOB1D_RADIAL,NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX,NCHUNKS, &
    OUTPUT_FILES,IFLAG_IN_FICTITIOUS_CUBE,NGLLX,NGLLY,NGLLZ,NSPEC2D_BOTTOM,MAX_STRING_LEN

  use meshfem3D_par, only: ibool,idoubling,is_on_a_slice_edge,xstore_glob,ystore_glob,zstore_glob

  use create_MPI_interfaces_par
  use MPI_inner_core_par

  implicit none

  integer :: MAX_NEIGHBOURS,max_nibool
  integer, dimension(MAX_NEIGHBOURS) :: my_neighbours,nibool_neighbours
  integer, dimension(max_nibool,MAX_NEIGHBOURS) :: ibool_neighbours

  ! local parameters
  ! temporary buffers for send and receive between faces of the slices and the chunks
  real(kind=CUSTOM_REAL), dimension(npoin2D_max_all_CM_IC) :: &
    buffer_send_faces_scalar,buffer_received_faces_scalar
  real(kind=CUSTOM_REAL),dimension(:),allocatable :: test_flag
  integer :: i,j,k,ispec,iglob,ier
  integer :: ndim_assemble
  !----------------------
  ! debug file output
  logical,parameter :: DEBUG = .false.
  !----------------------
  character(len=MAX_STRING_LEN) :: filename

  ! sets up MPI interfaces
  ! inner core
  if (myrank == 0 ) write(IMAIN,*) 'inner core MPI:'

  if (NPROCTOT > 1) then
    allocate(test_flag(NGLOB_INNER_CORE), &
            stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating test_flag inner core')

    ! sets flag to rank id (+1 to avoid problems with zero rank)
    test_flag(:) = 0.0
    do ispec = 1,NSPEC_INNER_CORE
      ! suppress fictitious elements in central cube
      if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
      ! sets flags
      do k = 1,NGLLZ
        do j = 1,NGLLY
          do i = 1,NGLLX
            iglob = ibool(i,j,k,ispec)
            test_flag(iglob) = myrank + 1.0
          enddo
        enddo
      enddo
    enddo

    ! assembles values
    call assemble_MPI_scalar_block(test_flag,NGLOB_INNER_CORE, &
              iproc_xi,iproc_eta,ichunk,addressing, &
              iboolleft_xi_inner_core,iboolright_xi_inner_core,iboolleft_eta_inner_core,iboolright_eta_inner_core, &
              npoin2D_faces_inner_core,npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
              iboolfaces_inner_core,iboolcorner_inner_core, &
              iprocfrom_faces,iprocto_faces,imsg_type, &
              iproc_master_corners,iproc_worker1_corners,iproc_worker2_corners, &
              buffer_send_faces_scalar,buffer_received_faces_scalar,npoin2D_max_all_CM_IC, &
              buffer_send_chunkcorn_scalar,buffer_recv_chunkcorn_scalar, &
              NUMMSGS_FACES,NUM_MSG_TYPES,NCORNERSCHUNKS, &
              NPROC_XI,NPROC_ETA,NGLOB1D_RADIAL(IREGION_INNER_CORE), &
              NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
              NGLOB2DMAX_XY,NCHUNKS)

    ! debug: idoubling inner core
    if (DEBUG) then
      write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_idoubling_inner_core_proc',myrank
      call write_VTK_data_elem_i(NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                 xstore_glob,ystore_glob,zstore_glob, &
                                 ibool,idoubling,filename)
      call synchronize_all()
    endif

    ! including central cube
    if (INCLUDE_CENTRAL_CUBE) then
      ! user output
      if (myrank == 0 ) write(IMAIN,*) 'inner core with central cube MPI:'

      ! test_flag is a scalar, not a vector
      ndim_assemble = 1

      ! use central cube buffers to assemble the inner core mass matrix with the central cube
      call assemble_MPI_central_cube_block(ichunk,nb_msgs_theor_in_cube, sender_from_slices_to_cube, &
                   npoin2D_cube_from_slices, buffer_all_cube_from_slices, &
                   buffer_slices, buffer_slices2, ibool_central_cube, &
                   receiver_cube_from_slices, ibool, &
                   idoubling, NSPEC_INNER_CORE, &
                   ibelm_bottom_inner_core, NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                   NGLOB_INNER_CORE, &
                   test_flag,ndim_assemble, &
                   iproc_eta,addressing,NCHUNKS,NPROC_XI,NPROC_ETA)

      ! frees array not needed anymore
      deallocate(ibelm_bottom_inner_core)

    endif

    ! removes own myrank id (+1)
    test_flag = test_flag - ( myrank + 1.0)
    where( test_flag < 0.0 ) test_flag = 0.0

    ! debug: in sequential order, for testing purpose
    !do i = 0,NPROCTOT - 1
    !  if (myrank == i) then
    !    ! gets new interfaces for inner_core without central cube yet
    !    ! determines neighbor rank for shared faces
    !    call get_MPI_interfaces(myrank,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
    !                          test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
    !                          num_interfaces_inner_core,max_nibool_interfaces_ic, &
    !                          max_nibool,MAX_NEIGHBOURS, &
    !                          ibool,is_on_a_slice_edge, &
    !                          IREGION_INNER_CORE,.false.,idoubling,INCLUDE_CENTRAL_CUBE, &
    !                          xstore_glob,ystore_glob,zstore_glob,NPROCTOT)
    !  endif
    !  call synchronize_all()
    !enddo

    ! gets new interfaces for inner_core without central cube yet
    ! determines neighbor rank for shared faces
    call get_MPI_interfaces(myrank,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                          test_flag,my_neighbours,nibool_neighbours,ibool_neighbours, &
                          num_interfaces_inner_core,max_nibool_interfaces_ic, &
                          max_nibool,MAX_NEIGHBOURS, &
                          ibool,is_on_a_slice_edge, &
                          IREGION_INNER_CORE,.false.,idoubling,INCLUDE_CENTRAL_CUBE, &
                          xstore_glob,ystore_glob,zstore_glob,NPROCTOT)

    deallocate(test_flag)
  else
    ! no interfaces
    num_interfaces_inner_core = 0
  endif

  ! stores MPI interfaces information
  allocate(my_neighbours_inner_core(num_interfaces_inner_core), &
          nibool_interfaces_inner_core(num_interfaces_inner_core), &
          stat=ier)
  if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array my_neighbours_inner_core etc.')
  my_neighbours_inner_core(:) = -1
  nibool_interfaces_inner_core(:) = 0

  ! copies interfaces arrays
  if (num_interfaces_inner_core > 0) then
    allocate(ibool_interfaces_inner_core(max_nibool_interfaces_ic,num_interfaces_inner_core), &
           stat=ier)
    if (ier /= 0 ) call exit_mpi(myrank,'Error allocating array ibool_interfaces_inner_core')
    ibool_interfaces_inner_core(:,:) = 0

    ! ranks of neighbour processes
    my_neighbours_inner_core(:) = my_neighbours(1:num_interfaces_inner_core)
    ! number of global ibool entries on each interface
    nibool_interfaces_inner_core(:) = nibool_neighbours(1:num_interfaces_inner_core)
    ! global iglob point ids on each interface
    ibool_interfaces_inner_core(:,:) = ibool_neighbours(1:max_nibool_interfaces_ic,1:num_interfaces_inner_core)
  else
    ! dummy allocation (fortran90 should allow allocate statement with zero array size)
    max_nibool_interfaces_ic = 0
    allocate(ibool_interfaces_inner_core(0,0),stat=ier)
  endif

  ! debug: saves MPI interfaces
  if (DEBUG) then
    do i = 1,num_interfaces_inner_core
      write(filename,'(a,i6.6,a,i2.2)') trim(OUTPUT_FILES)//'/MPI_points_inner_core_proc',myrank, &
                      '_',my_neighbours_inner_core(i)
      call write_VTK_data_points(NGLOB_INNER_CORE, &
                        xstore_glob,ystore_glob,zstore_glob, &
                        ibool_interfaces_inner_core(1:nibool_interfaces_inner_core(i),i), &
                        nibool_interfaces_inner_core(i),filename)
    enddo
    call synchronize_all()
  endif

  ! checks addressing
  call test_MPI_neighbours(IREGION_INNER_CORE, &
                              num_interfaces_inner_core,max_nibool_interfaces_ic, &
                              my_neighbours_inner_core,nibool_interfaces_inner_core, &
                              ibool_interfaces_inner_core)

  ! checks with assembly of test fields
  call test_MPI_ic()

  end subroutine setup_MPI_interfaces_ic

