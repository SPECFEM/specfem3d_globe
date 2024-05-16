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

  subroutine get_absorb(prname,iregion,iboun, &
                        nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! Stacey, define flags for absorbing boundaries

  use constants
  use meshfem_par, only: ADIOS_FOR_ARRAYS_SOLVER,nspec

  implicit none

  integer,intent(in) :: iregion

  integer,intent(in) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM

  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX),intent(inout) :: nimin,nimax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX),intent(inout) :: njmin,njmax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX),intent(inout) :: nkmin_xi
  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX),intent(inout) :: nkmin_eta

  logical,intent(in) :: iboun(6,nspec)

  ! global element numbering
  integer :: ispec

  ! counters to keep track of the number of elements on each of the
  ! five absorbing boundaries
  integer :: ispecb1,ispecb2,ispecb3,ispecb4,ispecb5
  integer :: ier

  ! processor identification
  character(len=MAX_STRING_LEN) :: prname

  ! initializes
  ispecb1 = 0
  ispecb2 = 0
  ispecb3 = 0
  ispecb4 = 0
  ispecb5 = 0

  do ispec = 1,nspec

    ! determine if the element falls on an absorbing boundary

    if (iboun(1,ispec)) then

      !   on boundary 1: xmin
      ispecb1 = ispecb1 + 1

      ! this is useful even if it is constant because it can be zero inside the slices
      njmin(1,ispecb1) = 1
      njmax(1,ispecb1) = NGLLY

      !   check for overlap with other boundaries
      nkmin_xi(1,ispecb1) = 1
      if (iboun(5,ispec)) nkmin_xi(1,ispecb1) = 2
    endif

    if (iboun(2,ispec)) then

      !   on boundary 2: xmax
      ispecb2 = ispecb2 + 1

      ! this is useful even if it is constant because it can be zero inside the slices
      njmin(2,ispecb2) = 1
      njmax(2,ispecb2) = NGLLY

      !   check for overlap with other boundaries
      nkmin_xi(2,ispecb2) = 1
      if (iboun(5,ispec)) nkmin_xi(2,ispecb2) = 2
    endif

    if (iboun(3,ispec)) then

      !   on boundary 3: ymin
      ispecb3 = ispecb3 + 1

      !   check for overlap with other boundaries
      nimin(1,ispecb3) = 1
      if (iboun(1,ispec)) nimin(1,ispecb3) = 2
      nimax(1,ispecb3) = NGLLX
      if (iboun(2,ispec)) nimax(1,ispecb3) = NGLLX-1
      nkmin_eta(1,ispecb3) = 1
      if (iboun(5,ispec)) nkmin_eta(1,ispecb3) = 2
    endif

    if (iboun(4,ispec)) then

      !   on boundary 4: ymax
      ispecb4 = ispecb4 + 1

      !   check for overlap with other boundaries
      nimin(2,ispecb4) = 1
      if (iboun(1,ispec)) nimin(2,ispecb4) = 2
      nimax(2,ispecb4) = NGLLX
      if (iboun(2,ispec)) nimax(2,ispecb4) = NGLLX-1
      nkmin_eta(2,ispecb4) = 1
      if (iboun(5,ispec)) nkmin_eta(2,ispecb4) = 2
    endif

    ! on boundary 5: bottom
    if (iboun(5,ispec)) ispecb5 = ispecb5 + 1

  enddo

  ! check theoretical value of elements at the bottom
  if (ispecb5 /= NSPEC2D_BOTTOM) then
    print *,'Error: invalid ispecb5:',ispecb5,NSPEC2D_BOTTOM,'region',iregion,'nspec',nspec
    call exit_MPI(myrank,'ispecb5 should equal NSPEC2D_BOTTOM in absorbing boundary detection')
  endif

  ! save these temporary arrays for the solver for Stacey conditions
  ! old version: left here for reference...
  if (.false.) then
    if (NSPEC2DMAX_YMIN_YMAX > 0) then
      ! This files will be saved with the help of ADIOS if the
      ! ADIOS_FOR_ARRAYS_SOLVER flag is set to true in the Par_file
      if (ADIOS_FOR_ARRAYS_SOLVER) then
        call get_absorb_adios(iregion, &
                              nimin, nimax, njmin, njmax, nkmin_xi, nkmin_eta, &
                              NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)
      else
        open(unit=IOUT,file=prname(1:len_trim(prname))//'stacey.old_format.bin', &
              status='unknown',form='unformatted',action='write',iostat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey.old_format.bin file')
        write(IOUT) nimin
        write(IOUT) nimax
        write(IOUT) njmin
        write(IOUT) njmax
        write(IOUT) nkmin_xi
        write(IOUT) nkmin_eta
        close(IOUT)
      endif
    endif
  endif

  end subroutine get_absorb

!
!-------------------------------------------------------------------------------------------------
!


  subroutine get_absorb_create_Stacey_boundary_arrays(iregion,NSPEC2D_BOTTOM)

! creates absorbing boundary arrays (similar to SPECFEM3D_Cartesian)
! - for future use to replace current looping over xmin/xmax/..

  use constants

  use meshfem_par, only: &
    myrank,NCHUNKS,ichunk

  use regions_mesh_par, only: &
    wxgll,wygll,wzgll

  use regions_mesh_par2, only: &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom, &
    normal_xmin,normal_xmax,normal_ymin,normal_ymax, &
    normal_bottom, &
    jacobian2D_xmin,jacobian2D_xmax,jacobian2D_ymin,jacobian2D_ymax, &
    jacobian2D_bottom, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax, &
    prname

  ! absorb
  use regions_mesh_par2, only: nimin, nimax, njmin, njmax, nkmin_xi,nkmin_eta, &
    num_abs_boundary_faces,abs_boundary_ispec,abs_boundary_npoin, &
    abs_boundary_ijk,abs_boundary_normal,abs_boundary_jacobian2Dw

  ! input parameters
  use shared_parameters, only: REGIONAL_MESH_CUTOFF,ADIOS_FOR_ARRAYS_SOLVER

  ! debug
  use meshfem_par, only: nspec,nglob,ibool,xstore_glob,ystore_glob,zstore_glob

  implicit none

  integer,intent(in) :: iregion
  integer,intent(in) :: NSPEC2D_BOTTOM

  ! local parameters
  integer :: iface,igll,i,j,k,ispec,ispec2D,ier
  real(kind=CUSTOM_REAL) :: nx,ny,nz
  double precision :: weight
  double precision, dimension(NGLLX,NGLLY) :: wgllwgll_xy
  double precision, dimension(NGLLX,NGLLZ) :: wgllwgll_xz
  double precision, dimension(NGLLY,NGLLZ) :: wgllwgll_yz
  logical :: add_bottom_boundary

  ! debugging
  character(len=MAX_STRING_LEN) :: filename
  integer :: ipoints,iglob
  integer, dimension(:), allocatable :: tmp_array_points,tmp_array_elements

  logical, parameter :: DEBUG = .false.

  ! initializes
  num_abs_boundary_faces = 0

  ! checks if anything to do
  ! no Stacey boundary on inner core
  if (iregion == IREGION_INNER_CORE) return

  ! setup GLL weights
  ! weights on surfaces
  do i = 1,NGLLX
    do j = 1,NGLLY
       wgllwgll_xy(i,j) = wxgll(i)*wygll(j)
    enddo
  enddo
  do i = 1,NGLLX
    do k = 1,NGLLZ
       wgllwgll_xz(i,k) = wxgll(i)*wzgll(k)
    enddo
  enddo
  do j = 1,NGLLY
    do k = 1,NGLLZ
       wgllwgll_yz(j,k) = wygll(j)*wzgll(k)
    enddo
  enddo

  ! checks if element arrays available
  if (.not. allocated(nimin)) call exit_MPI(myrank,'Error Stacey element arrays nimin,.. not allocated')

  ! initializes bottom boundary flag
  add_bottom_boundary = .false.

  ! adds bottom surface only for outer core or in case of regional mesh cut-offs for crust-mantle region
  if (iregion == IREGION_OUTER_CORE .or. &
      (REGIONAL_MESH_CUTOFF .and. iregion == IREGION_CRUST_MANTLE)) &
    add_bottom_boundary = .true.

  ! counts number of faces
  ! xmin
  if (NCHUNKS == 1 .or. ichunk == CHUNK_AC) then
    do ispec2D = 1,nspec2D_xmin
      ! exclude elements that are not on absorbing edges
      if (nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle
      num_abs_boundary_faces = num_abs_boundary_faces + 1
    enddo
  endif
  ! xmax
  if (NCHUNKS == 1 .or. ichunk == CHUNK_AB) then
    do ispec2D = 1,nspec2D_xmax
      ! exclude elements that are not on absorbing edges
      if (nkmin_xi(2,ispec2D) == 0 .or. njmin(2,ispec2D) == 0) cycle
      num_abs_boundary_faces = num_abs_boundary_faces + 1
    enddo
  endif
  ! ymin
  do ispec2D = 1,nspec2D_ymin
    ! exclude elements that are not on absorbing edges
    if (nkmin_eta(1,ispec2D) == 0 .or. nimin(1,ispec2D) == 0) cycle
    num_abs_boundary_faces = num_abs_boundary_faces + 1
  enddo
  ! ymax
  do ispec2D = 1,nspec2D_ymax
    ! exclude elements that are not on absorbing edges
    if (nkmin_eta(2,ispec2D) == 0 .or. nimin(2,ispec2D) == 0) cycle
    num_abs_boundary_faces = num_abs_boundary_faces + 1
  enddo
  ! zmin
  if (add_bottom_boundary) then
    do ispec2D = 1,NSPEC2D_BOTTOM
      num_abs_boundary_faces = num_abs_boundary_faces + 1
    enddo
  endif

  ! debug
  if (DEBUG) then
    print *,'debug: save_arrays_Stacey_boundary_to_debug: iregion ',iregion
    print *,'debug: NCHUNKS ',NCHUNKS,' ichunk: ',ichunk, ' - CHUNK_AC = ',CHUNK_AC,' CHUNK_AB = ',CHUNK_AB
    print *,'debug: nspec2D_xmin / xmax = ',nspec2D_xmin,nspec2D_xmax
    print *,'debug: nspec2D_ymin / ymax = ',nspec2D_ymin,nspec2D_ymax
    print *,'debug: nspec2D_zmin = ',NSPEC2D_BOTTOM
    print *,'debug: REGIONAL_MESH_CUTOFF: ',REGIONAL_MESH_CUTOFF,' add_bottom_boundary: ',add_bottom_boundary
    print *,'debug: num_abs_boundary_faces = ',num_abs_boundary_faces
  endif

  ! allocates absorbing boundary arrays
  allocate(abs_boundary_ispec(num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1518')
  allocate(abs_boundary_ijk(3,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1519')
  allocate(abs_boundary_jacobian2Dw(NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1520')
  allocate(abs_boundary_normal(NDIM,NGLLSQUARE,num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1521')
  allocate(abs_boundary_npoin(num_abs_boundary_faces),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1521')
  if (ier /= 0) stop 'Error allocating array abs_boundary_ispec etc.'

  abs_boundary_ispec(:) = 0; abs_boundary_ijk(:,:,:) = 0; abs_boundary_npoin(:) = 0
  abs_boundary_jacobian2Dw(:,:) = 0.0_CUSTOM_REAL; abs_boundary_normal(:,:,:) = 0.0_CUSTOM_REAL

  ! counter for each hexahedral face being an absorbing boundary
  iface = 0

  ! sets absorbing arrays
  ! xmin
  ! if two chunks exclude this face for one of them
  if (NCHUNKS == 1 .or. ichunk == CHUNK_AC) then
    do ispec2D = 1,nspec2D_xmin
      ispec = ibelm_xmin(ispec2D)

      ! exclude elements that are not on absorbing edges
      if (nkmin_xi(1,ispec2D) == 0 .or. njmin(1,ispec2D) == 0) cycle

      ! increase counter
      iface = iface + 1
      abs_boundary_ispec(iface) = ispec
      igll = 0

      i = 1
      do k = nkmin_xi(1,ispec2D),NGLLZ
        do j = njmin(1,ispec2D),njmax(1,ispec2D)

          nx = normal_xmin(1,j,k,ispec2D)
          ny = normal_xmin(2,j,k,ispec2D)
          nz = normal_xmin(3,j,k,ispec2D)

          weight = jacobian2D_xmin(j,k,ispec2D)*wgllwgll_yz(j,k)

          ! stores arrays
          igll = igll + 1
          abs_boundary_ijk(1,igll,iface) = i
          abs_boundary_ijk(2,igll,iface) = j
          abs_boundary_ijk(3,igll,iface) = k
          abs_boundary_normal(1,igll,iface) = nx
          abs_boundary_normal(2,igll,iface) = ny
          abs_boundary_normal(3,igll,iface) = nz
          abs_boundary_jacobian2Dw(igll,iface) = real(weight,kind=CUSTOM_REAL)
        enddo
      enddo
      abs_boundary_npoin(iface) = igll
    enddo
  endif

  ! xmax
  ! if two chunks exclude this face for one of them
  if (NCHUNKS == 1 .or. ichunk == CHUNK_AB) then
    do ispec2D = 1,nspec2D_xmax
      ispec = ibelm_xmax(ispec2D)
      ! exclude elements that are not on absorbing edges
      if (nkmin_xi(2,ispec2D) == 0 .or. njmin(2,ispec2D) == 0) cycle

      ! increase counter
      iface = iface + 1
      abs_boundary_ispec(iface) = ispec
      igll = 0

      i = NGLLX
      do k = nkmin_xi(2,ispec2D),NGLLZ
        do j = njmin(2,ispec2D),njmax(2,ispec2D)
          nx = normal_xmax(1,j,k,ispec2D)
          ny = normal_xmax(2,j,k,ispec2D)
          nz = normal_xmax(3,j,k,ispec2D)
          weight = jacobian2D_xmax(j,k,ispec2D)*wgllwgll_yz(j,k)

          ! stores arrays
          igll = igll + 1
          abs_boundary_ijk(1,igll,iface) = i
          abs_boundary_ijk(2,igll,iface) = j
          abs_boundary_ijk(3,igll,iface) = k
          abs_boundary_normal(1,igll,iface) = nx
          abs_boundary_normal(2,igll,iface) = ny
          abs_boundary_normal(3,igll,iface) = nz
          abs_boundary_jacobian2Dw(igll,iface) = real(weight,kind=CUSTOM_REAL)
        enddo
      enddo
      abs_boundary_npoin(iface) = igll
    enddo
  endif

  ! ymin
  do ispec2D = 1,nspec2D_ymin
    ispec = ibelm_ymin(ispec2D)
    ! exclude elements that are not on absorbing edges
    if (nkmin_eta(1,ispec2D) == 0 .or. nimin(1,ispec2D) == 0) cycle

    ! increase counter
    iface = iface + 1
    abs_boundary_ispec(iface) = ispec
    igll = 0

    j = 1
    do k = nkmin_eta(1,ispec2D),NGLLZ
      do i = nimin(1,ispec2D),nimax(1,ispec2D)
        nx = normal_ymin(1,i,k,ispec2D)
        ny = normal_ymin(2,i,k,ispec2D)
        nz = normal_ymin(3,i,k,ispec2D)
        weight = jacobian2D_ymin(i,k,ispec2D)*wgllwgll_xz(i,k)

        ! stores arrays
        igll = igll + 1
        abs_boundary_ijk(1,igll,iface) = i
        abs_boundary_ijk(2,igll,iface) = j
        abs_boundary_ijk(3,igll,iface) = k
        abs_boundary_normal(1,igll,iface) = nx
        abs_boundary_normal(2,igll,iface) = ny
        abs_boundary_normal(3,igll,iface) = nz
        abs_boundary_jacobian2Dw(igll,iface) = real(weight,kind=CUSTOM_REAL)
      enddo
    enddo
    abs_boundary_npoin(iface) = igll
  enddo

  ! ymax
  do ispec2D = 1,nspec2D_ymax
    ispec = ibelm_ymax(ispec2D)
    ! exclude elements that are not on absorbing edges
    if (nkmin_eta(2,ispec2D) == 0 .or. nimin(2,ispec2D) == 0) cycle

    ! increase counter
    iface = iface + 1
    abs_boundary_ispec(iface) = ispec
    igll = 0

    j = NGLLY
    do k = nkmin_eta(2,ispec2D),NGLLZ
      do i = nimin(2,ispec2D),nimax(2,ispec2D)
        nx = normal_ymax(1,i,k,ispec2D)
        ny = normal_ymax(2,i,k,ispec2D)
        nz = normal_ymax(3,i,k,ispec2D)
        weight = jacobian2D_ymax(i,k,ispec2D)*wgllwgll_xz(i,k)

        ! stores arrays
        igll = igll + 1
        abs_boundary_ijk(1,igll,iface) = i
        abs_boundary_ijk(2,igll,iface) = j
        abs_boundary_ijk(3,igll,iface) = k
        abs_boundary_normal(1,igll,iface) = nx
        abs_boundary_normal(2,igll,iface) = ny
        abs_boundary_normal(3,igll,iface) = nz
        abs_boundary_jacobian2Dw(igll,iface) = real(weight,kind=CUSTOM_REAL)
      enddo
    enddo
    abs_boundary_npoin(iface) = igll
  enddo

  ! zmin
  if (add_bottom_boundary) then
    do ispec2D = 1,NSPEC2D_BOTTOM
      ispec = ibelm_bottom(ispec2D)

      ! increase counter
      iface = iface + 1
      abs_boundary_ispec(iface) = ispec
      igll = 0

      k = 1
      do j = 1,NGLLY
        do i = 1,NGLLX
          nx = normal_bottom(1,i,k,ispec2D)
          ny = normal_bottom(2,i,k,ispec2D)
          nz = normal_bottom(3,i,k,ispec2D)
          weight = jacobian2D_bottom(i,k,ispec2D)*wgllwgll_xy(i,j)

          ! stores arrays
          igll = igll + 1
          abs_boundary_ijk(1,igll,iface) = i
          abs_boundary_ijk(2,igll,iface) = j
          abs_boundary_ijk(3,igll,iface) = k
          abs_boundary_normal(1,igll,iface) = nx
          abs_boundary_normal(2,igll,iface) = ny
          abs_boundary_normal(3,igll,iface) = nz
          abs_boundary_jacobian2Dw(igll,iface) = real(weight,kind=CUSTOM_REAL)
        enddo
      enddo
      abs_boundary_npoin(iface) = igll

    enddo
  endif

  ! check
  if (iface /= num_abs_boundary_faces) then
    print *,'Error: iface = ',iface,' should match num_abs_boundary_faces = ',num_abs_boundary_faces
    call exit_MPI(myrank,'Invalid number of absorbing boundary faces')
  endif

  ! saves arrays to file
  if (ADIOS_FOR_ARRAYS_SOLVER) then
    call get_absorb_stacey_boundary_adios(iregion, num_abs_boundary_faces, &
                                          abs_boundary_ispec,abs_boundary_npoin, &
                                          abs_boundary_ijk,abs_boundary_normal,abs_boundary_jacobian2Dw)
  else
    ! binary format
    open(unit=IOUT,file=prname(1:len_trim(prname))//'stacey.bin', &
         status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey_boundary_arrays.bin file')

    write(IOUT) num_abs_boundary_faces

    if (num_abs_boundary_faces > 0) then
      write(IOUT) abs_boundary_ispec
      write(IOUT) abs_boundary_npoin
      write(IOUT) abs_boundary_ijk
      write(IOUT) abs_boundary_jacobian2Dw
      ! normals only needed for elastic boundary conditions
      if (iregion == IREGION_CRUST_MANTLE) then
        write(IOUT) abs_boundary_normal
      endif
    endif

    close(IOUT)
  endif

  ! debug
  if (DEBUG) then
    ! outputs boundary points for visual inspection
    ! count points
    ipoints = 0
    do iface = 1,num_abs_boundary_faces
      ipoints = ipoints + abs_boundary_npoin(iface)
    enddo
    print *,'debug: number of points ',ipoints

    allocate(tmp_array_points(ipoints),stat=ier)
    if (ier /= 0) stop 'Error allocating tmp_array_points'
    tmp_array_points(:) = 0

    ipoints = 0
    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      do igll = 1,abs_boundary_npoin(iface)
        i = abs_boundary_ijk(1,igll,iface)
        j = abs_boundary_ijk(2,igll,iface)
        k = abs_boundary_ijk(3,igll,iface)
        iglob = ibool(i,j,k,ispec)
        ! sets point
        ipoints = ipoints + 1
        tmp_array_points(ipoints) = iglob
      enddo
    enddo

    filename = trim(prname) // 'stacey_points'
    call write_VTK_data_points(nglob,xstore_glob,ystore_glob,zstore_glob, &
                               tmp_array_points,ipoints,filename)

    deallocate(tmp_array_points)

    ! sets element flags
    allocate(tmp_array_elements(nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating tmp_array_elements'
    tmp_array_elements(:) = 0

    do iface = 1,num_abs_boundary_faces
      ispec = abs_boundary_ispec(iface)
      tmp_array_elements(ispec) = tmp_array_elements(ispec) + 1
    enddo

    filename = trim(prname) // 'stacey_elements'
    call write_VTK_data_elem_i(nspec,nglob, &
                               xstore_glob,ystore_glob,zstore_glob, &
                               ibool,tmp_array_elements,filename)

    deallocate(tmp_array_elements)
  endif

  end subroutine get_absorb_create_Stacey_boundary_arrays
