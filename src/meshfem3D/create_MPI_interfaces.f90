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

  subroutine create_MPI_interfaces(iregion_code)

  implicit none

  integer,intent(in):: iregion_code

  ! sets up arrays
  call cmi_allocate_addressing(iregion_code)

  ! reads "iboolleft_..txt", "iboolright_..txt" (and "list_messages_..txt", "buffer_...txt") files and sets up MPI buffers
  call cmi_get_buffers(iregion_code)

  end subroutine create_MPI_interfaces

!
!-------------------------------------------------------------------------------------------------
!

  subroutine cmi_allocate_addressing(iregion_code)

  use constants, only: CUSTOM_REAL,NUMCORNERS_SHARED,myrank

  use meshfem_par, only: ibool, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
    NSPEC2D_BOTTOM,NSPEC2D_TOP, &
    NSPEC_REGIONS,NGLOB_REGIONS, &
    NGLOB1D_RADIAL

  use MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  use MPI_trinfinite_par
  use MPI_infinite_par

  implicit none

  integer,intent(in):: iregion_code

  ! local parameters
  integer :: ier

  ! parameters from header file
  NGLOB1D_RADIAL_CM = NGLOB1D_RADIAL(IREGION_CRUST_MANTLE)
  NGLOB1D_RADIAL_OC = NGLOB1D_RADIAL(IREGION_OUTER_CORE)
  NGLOB1D_RADIAL_IC = NGLOB1D_RADIAL(IREGION_INNER_CORE)

  ! initializes
  NSPEC_CRUST_MANTLE = 0
  NGLOB_CRUST_MANTLE = 0

  NSPEC_OUTER_CORE = 0
  NGLOB_OUTER_CORE = 0

  NSPEC_INNER_CORE = 0
  NGLOB_INNER_CORE = 0

  NSPEC_TRINFINITE = 0
  NGLOB_TRINFINITE = 0

  NSPEC_INFINITE = 0
  NGLOB_INFINITE = 0

  select case (iregion_code)

  case (IREGION_CRUST_MANTLE)
    NGLOB2DMAX_XMIN_XMAX_CM = NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE)
    NGLOB2DMAX_YMIN_YMAX_CM = NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE)

    NSPEC2DMAX_XMIN_XMAX_CM = NSPEC2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE)
    NSPEC2DMAX_YMIN_YMAX_CM = NSPEC2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE)
    NSPEC2D_BOTTOM_CM = NSPEC2D_BOTTOM(IREGION_CRUST_MANTLE)
    NSPEC2D_TOP_CM = NSPEC2D_TOP(IREGION_CRUST_MANTLE)

    NSPEC_CRUST_MANTLE = NSPEC_REGIONS(IREGION_CRUST_MANTLE)
    NGLOB_CRUST_MANTLE = NGLOB_REGIONS(IREGION_CRUST_MANTLE)

  case (IREGION_OUTER_CORE)
    NGLOB2DMAX_XMIN_XMAX_OC = NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE)
    NGLOB2DMAX_YMIN_YMAX_OC = NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE)

    NSPEC2DMAX_XMIN_XMAX_OC = NSPEC2DMAX_XMIN_XMAX(IREGION_OUTER_CORE)
    NSPEC2DMAX_YMIN_YMAX_OC = NSPEC2DMAX_YMIN_YMAX(IREGION_OUTER_CORE)
    NSPEC2D_BOTTOM_OC = NSPEC2D_BOTTOM(IREGION_OUTER_CORE)
    NSPEC2D_TOP_OC = NSPEC2D_TOP(IREGION_OUTER_CORE)

    NSPEC_OUTER_CORE = NSPEC_REGIONS(IREGION_OUTER_CORE)
    NGLOB_OUTER_CORE = NGLOB_REGIONS(IREGION_OUTER_CORE)

  case (IREGION_INNER_CORE)
    NGLOB2DMAX_XMIN_XMAX_IC = NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE)
    NGLOB2DMAX_YMIN_YMAX_IC = NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE)

    NSPEC2DMAX_XMIN_XMAX_IC = NSPEC2DMAX_XMIN_XMAX(IREGION_INNER_CORE)
    NSPEC2DMAX_YMIN_YMAX_IC = NSPEC2DMAX_YMIN_YMAX(IREGION_INNER_CORE)
    NSPEC2D_BOTTOM_IC = NSPEC2D_BOTTOM(IREGION_INNER_CORE)
    NSPEC2D_TOP_IC = NSPEC2D_TOP(IREGION_INNER_CORE)

    NSPEC_INNER_CORE = NSPEC_REGIONS(IREGION_INNER_CORE)
    NGLOB_INNER_CORE = NGLOB_REGIONS(IREGION_INNER_CORE)

  case(IREGION_TRINFINITE)
    NGLOB2DMAX_XMIN_XMAX_TRINF = NGLOB2DMAX_XMIN_XMAX(IREGION_TRINFINITE)
    NGLOB2DMAX_YMIN_YMAX_TRINF = NGLOB2DMAX_YMIN_YMAX(IREGION_TRINFINITE)

    NSPEC2DMAX_XMIN_XMAX_TRINF = NSPEC2DMAX_XMIN_XMAX(IREGION_TRINFINITE)
    NSPEC2DMAX_YMIN_YMAX_TRINF = NSPEC2DMAX_YMIN_YMAX(IREGION_TRINFINITE)
    NSPEC2D_BOTTOM_TRINF = NSPEC2D_BOTTOM(IREGION_TRINFINITE)
    NSPEC2D_TOP_TRINF = NSPEC2D_TOP(IREGION_TRINFINITE)

    NSPEC_TRINFINITE = NSPEC_REGIONS(IREGION_TRINFINITE)
    NGLOB_TRINFINITE = NGLOB_REGIONS(IREGION_TRINFINITE)

  case(IREGION_INFINITE)
    NGLOB2DMAX_XMIN_XMAX_INF = NGLOB2DMAX_XMIN_XMAX(IREGION_INFINITE)
    NGLOB2DMAX_YMIN_YMAX_INF = NGLOB2DMAX_YMIN_YMAX(IREGION_INFINITE)

    NSPEC2DMAX_XMIN_XMAX_INF = NSPEC2DMAX_XMIN_XMAX(IREGION_INFINITE)
    NSPEC2DMAX_YMIN_YMAX_INF = NSPEC2DMAX_YMIN_YMAX(IREGION_INFINITE)
    NSPEC2D_BOTTOM_INF = NSPEC2D_BOTTOM(IREGION_INFINITE)
    NSPEC2D_TOP_INF = NSPEC2D_TOP(IREGION_INFINITE)

    NSPEC_INFINITE = NSPEC_REGIONS(IREGION_INFINITE)
    NGLOB_INFINITE = NGLOB_REGIONS(IREGION_INFINITE)

  case default
    stop 'Error iregion_code value not recognized'
  end select

  ! checks ibool for mesh
  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! check that the number of points in this slice is correct
    if (NSPEC_CRUST_MANTLE > 0) then
      if (minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_CRUST_MANTLE) &
        call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in crust and mantle')
    endif

  case (IREGION_OUTER_CORE)
    ! check that the number of points in this slice is correct
    if (NSPEC_OUTER_CORE > 0) then
      if (minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_OUTER_CORE) &
        call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in outer core')
    endif

  case (IREGION_INNER_CORE)
    ! check that the number of points in this slice is correct
    if (NSPEC_INNER_CORE > 0) then
      if (minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_INNER_CORE) &
        call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in inner core')
    endif

  case (IREGION_TRINFINITE)
    ! check that the number of points in this slice is correct
    if (NSPEC_TRINFINITE > 0) then
      if (minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_TRINFINITE) &
        call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in transition inifinite region')
    endif

  case (IREGION_INFINITE)
    ! check that the number of points in this slice is correct
    if (NSPEC_INFINITE > 0) then
      if (minval(ibool(:,:,:,:)) /= 1 .or. maxval(ibool(:,:,:,:)) /= NGLOB_INFINITE) &
        call exit_MPI(myrank,'incorrect global numbering: iboolmax does not equal nglob in inifinite region')
    endif

  end select

  ! allocates arrays
  allocate(buffer_send_chunkcorn_scalar(NGLOB1D_RADIAL_CM), &
           buffer_recv_chunkcorn_scalar(NGLOB1D_RADIAL_CM),stat=ier)
  if (ier /= 0) stop 'Error allocating buffer buffer_send_chunkcorn_scalar,.. arrays'
  buffer_send_chunkcorn_scalar(:) = 0.0_CUSTOM_REAL
  buffer_recv_chunkcorn_scalar(:) = 0.0_CUSTOM_REAL

  allocate(buffer_send_chunkcorn_vector(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC), &
           buffer_recv_chunkcorn_vector(NDIM,NGLOB1D_RADIAL_CM + NGLOB1D_RADIAL_IC),stat=ier)
  if (ier /= 0) stop 'Error allocating buffer buffer_send_chunkcorn_vector,.. arrays'
  buffer_send_chunkcorn_vector(:,:) = 0.0_CUSTOM_REAL
  buffer_recv_chunkcorn_vector(:,:) = 0.0_CUSTOM_REAL

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! crust mantle
    allocate(iboolcorner_crust_mantle(NGLOB1D_RADIAL_CM,NUMCORNERS_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolcorner_crust_mantle array'
    iboolcorner_crust_mantle(:,:) = 0

    allocate(iboolleft_xi_crust_mantle(NGLOB2DMAX_XMIN_XMAX_CM), &
             iboolright_xi_crust_mantle(NGLOB2DMAX_XMIN_XMAX_CM),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_xi_crust_mantle array'
    iboolleft_xi_crust_mantle(:) = 0
    iboolright_xi_crust_mantle(:) = 0

    allocate(iboolleft_eta_crust_mantle(NGLOB2DMAX_YMIN_YMAX_CM), &
             iboolright_eta_crust_mantle(NGLOB2DMAX_YMIN_YMAX_CM),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_eta_crust_mantle array'
    iboolleft_eta_crust_mantle(:) = 0
    iboolright_eta_crust_mantle(:) = 0

    allocate(iboolfaces_crust_mantle(NGLOB2DMAX_XY,NUMFACES_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolfaces_crust_mantle array'
    iboolfaces_crust_mantle(:,:) = 0

  case (IREGION_OUTER_CORE)
    ! outer core
    allocate(iboolcorner_outer_core(NGLOB1D_RADIAL_OC,NUMCORNERS_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolcorner_outer_core array'
    iboolcorner_outer_core(:,:) = 0

    allocate(iboolleft_xi_outer_core(NGLOB2DMAX_XMIN_XMAX_OC), &
             iboolright_xi_outer_core(NGLOB2DMAX_XMIN_XMAX_OC),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_xi_outer_core array'
    iboolleft_xi_outer_core(:) = 0
    iboolright_xi_outer_core(:) = 0

    allocate(iboolleft_eta_outer_core(NGLOB2DMAX_YMIN_YMAX_OC), &
             iboolright_eta_outer_core(NGLOB2DMAX_YMIN_YMAX_OC),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_eta_outer_core array'
    iboolleft_eta_outer_core(:) = 0
    iboolright_eta_outer_core(:) = 0

    allocate(iboolfaces_outer_core(NGLOB2DMAX_XY,NUMFACES_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolfaces_outer_core array'
    iboolfaces_outer_core(:,:) = 0

  case (IREGION_INNER_CORE)
    ! inner core
    allocate(iboolcorner_inner_core(NGLOB1D_RADIAL_IC,NUMCORNERS_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolcorner_inner_core array'
    iboolcorner_inner_core(:,:) = 0

    allocate(iboolleft_xi_inner_core(NGLOB2DMAX_XMIN_XMAX_IC), &
             iboolright_xi_inner_core(NGLOB2DMAX_XMIN_XMAX_IC),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_xi_inner_core array'
    iboolleft_xi_inner_core(:) = 0
    iboolright_xi_inner_core(:) = 0

    allocate(iboolleft_eta_inner_core(NGLOB2DMAX_YMIN_YMAX_IC), &
             iboolright_eta_inner_core(NGLOB2DMAX_YMIN_YMAX_IC),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_eta_inner_core array'
    iboolleft_eta_inner_core(:) = 0
    iboolright_eta_inner_core(:) = 0

    allocate(iboolfaces_inner_core(NGLOB2DMAX_XY,NUMFACES_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolfaces_inner_core array'
    iboolfaces_inner_core(:,:) = 0

  case( IREGION_TRINFINITE )
    ! transition infinite
    allocate(iboolcorner_trinfinite(NGLOB1D_RADIAL_TRINF,NUMCORNERS_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolcorner_trinfinite array'
    iboolcorner_trinfinite(:,:) = 0

    allocate(iboolleft_xi_trinfinite(NGLOB2DMAX_XMIN_XMAX_TRINF), &
             iboolright_xi_trinfinite(NGLOB2DMAX_XMIN_XMAX_TRINF),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_xi_trinfinite array'
    iboolleft_xi_trinfinite(:) = 0
    iboolright_xi_trinfinite(:) = 0

    allocate(iboolleft_eta_trinfinite(NGLOB2DMAX_YMIN_YMAX_TRINF), &
             iboolright_eta_trinfinite(NGLOB2DMAX_YMIN_YMAX_TRINF),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_eta_trinfinite array'
    iboolleft_eta_trinfinite(:) = 0
    iboolright_eta_trinfinite(:) = 0

    allocate(iboolfaces_trinfinite(NGLOB2DMAX_XY,NUMFACES_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolfaces_inner_core array'
    iboolfaces_trinfinite(:,:) = 0

  case( IREGION_INFINITE )
    ! infinite
    allocate(iboolcorner_infinite(NGLOB1D_RADIAL_INF,NUMCORNERS_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolcorner_infinite array'
    iboolcorner_infinite(:,:) = 0

    allocate(iboolleft_xi_infinite(NGLOB2DMAX_XMIN_XMAX_INF), &
             iboolright_xi_infinite(NGLOB2DMAX_XMIN_XMAX_INF),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_xi_infinite array'
    iboolleft_xi_infinite(:) = 0
    iboolright_xi_infinite(:) = 0

    allocate(iboolleft_eta_infinite(NGLOB2DMAX_YMIN_YMAX_INF), &
             iboolright_eta_infinite(NGLOB2DMAX_YMIN_YMAX_INF),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolleft_eta_infinite array'
    iboolleft_eta_infinite(:) = 0
    iboolright_eta_infinite(:) = 0

    allocate(iboolfaces_infinite(NGLOB2DMAX_XY,NUMFACES_SHARED),stat=ier)
    if (ier /= 0) stop 'Error allocating iboolfaces_infinite array'
    iboolfaces_infinite(:,:) = 0

  end select

  ! synchronize processes
  call synchronize_all()

  end subroutine cmi_allocate_addressing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine cmi_get_buffers(iregion_code)

  use meshfem_par, only: myrank,MAX_STRING_LEN, &
    NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
    NGLOB1D_RADIAL,NSPEC2D_BOTTOM, &
    NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX, &
    NPROC_XI,NPROC_ETA,NCHUNKS,OUTPUT_FILES,INCLUDE_CENTRAL_CUBE, &
    iproc_xi,iproc_eta,ichunk,addressing, &
    xstore_glob,ystore_glob,zstore_glob

  use meshfem_par, only: &
    ibool,idoubling,is_on_a_slice_edge

  use regions_mesh_par2, only: &
    ibelm_xmin,ibelm_xmax,ibelm_ymin,ibelm_ymax,ibelm_bottom,ibelm_top, &
    nspec2D_xmin,nspec2D_xmax,nspec2D_ymin,nspec2D_ymax

  use MPI_interfaces_par

  use MPI_crust_mantle_par
  use MPI_outer_core_par
  use MPI_inner_core_par

  use MPI_trinfinite_par
  use MPI_infinite_par

  implicit none

  integer,intent(in):: iregion_code

  ! local parameters
  integer :: ier
  ! for central cube buffers
  integer :: nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
            nspec2D_ymin_inner_core,nspec2D_ymax_inner_core
  integer, dimension(:),allocatable :: ibelm_xmin_inner_core,ibelm_xmax_inner_core
  integer, dimension(:),allocatable :: ibelm_ymin_inner_core,ibelm_ymax_inner_core
  integer, dimension(:),allocatable :: ibelm_top_inner_core

  ! debug file output
  character(len=MAX_STRING_LEN) :: filename
  logical,parameter :: DEBUG = .false.

  ! gets 2-D addressing for summation between slices with MPI

  select case (iregion_code)
  case (IREGION_CRUST_MANTLE)
    ! mantle and crust
    if (NSPEC_CRUST_MANTLE > 0) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'crust/mantle region:'
        call flush_IMAIN()
      endif
      call cmi_read_buffer_data(IREGION_CRUST_MANTLE, &
                                NGLOB2DMAX_XMIN_XMAX(IREGION_CRUST_MANTLE), &
                                NGLOB2DMAX_YMIN_YMAX(IREGION_CRUST_MANTLE), &
                                NGLOB1D_RADIAL(IREGION_CRUST_MANTLE), &
                                iboolleft_xi_crust_mantle,iboolright_xi_crust_mantle, &
                                iboolleft_eta_crust_mantle,iboolright_eta_crust_mantle, &
                                npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
                                iboolfaces_crust_mantle,npoin2D_faces_crust_mantle, &
                                iboolcorner_crust_mantle)

      ! note: fix_... routines below update is_on_a_slice_edge_.. arrays:
      !          assign flags for each element which is on a rim of the slice
      !          thus, they include elements on top and bottom not shared with other MPI partitions
      !
      !          we will re-set these flags when setting up inner/outer elements, but will
      !          use these arrays for now as initial guess for the search for elements which share a global point
      !          between different MPI processes
      call fix_non_blocking_slices(is_on_a_slice_edge, &
                                   iboolright_xi_crust_mantle,iboolleft_xi_crust_mantle, &
                                   iboolright_eta_crust_mantle,iboolleft_eta_crust_mantle, &
                                   npoin2D_xi_crust_mantle,npoin2D_eta_crust_mantle, &
                                   ibool, &
                                   NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE,NGLOB2DMAX_XMIN_XMAX_CM,NGLOB2DMAX_YMIN_YMAX_CM)

      ! debug: saves element flags
      if (DEBUG) then
        write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_crust_mantle_proc',myrank
        call write_VTK_data_elem_l(NSPEC_CRUST_MANTLE,NGLOB_CRUST_MANTLE, &
                                   xstore_glob,ystore_glob,zstore_glob, &
                                   ibool,is_on_a_slice_edge,filename)
      endif

      ! added this to reduce the size of the buffers
      ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
      !npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_crust_mantle(:) + npoin2D_xi_inner_core(:)), &
      !                            maxval(npoin2D_eta_crust_mantle(:) + npoin2D_eta_inner_core(:)))
      npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_crust_mantle(:)), &
                                  maxval(npoin2D_eta_crust_mantle(:)))
    endif

  case (IREGION_OUTER_CORE)
    ! outer core
    if (NSPEC_OUTER_CORE > 0) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'outer core region:'
        call flush_IMAIN()
      endif
      call cmi_read_buffer_data(IREGION_OUTER_CORE, &
                                NGLOB2DMAX_XMIN_XMAX(IREGION_OUTER_CORE), &
                                NGLOB2DMAX_YMIN_YMAX(IREGION_OUTER_CORE), &
                                NGLOB1D_RADIAL(IREGION_OUTER_CORE), &
                                iboolleft_xi_outer_core,iboolright_xi_outer_core, &
                                iboolleft_eta_outer_core,iboolright_eta_outer_core, &
                                npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
                                iboolfaces_outer_core,npoin2D_faces_outer_core, &
                                iboolcorner_outer_core)

      ! note: fix_... routines below update is_on_a_slice_edge_.. arrays:
      !          assign flags for each element which is on a rim of the slice
      !          thus, they include elements on top and bottom not shared with other MPI partitions
      !
      !          we will re-set these flags when setting up inner/outer elements, but will
      !          use these arrays for now as initial guess for the search for elements which share a global point
      !          between different MPI processes
      call fix_non_blocking_slices(is_on_a_slice_edge, &
                                   iboolright_xi_outer_core,iboolleft_xi_outer_core, &
                                   iboolright_eta_outer_core,iboolleft_eta_outer_core, &
                                   npoin2D_xi_outer_core,npoin2D_eta_outer_core, &
                                   ibool, &
                                   NSPEC_OUTER_CORE,NGLOB_OUTER_CORE,NGLOB2DMAX_XMIN_XMAX_OC,NGLOB2DMAX_YMIN_YMAX_OC)

      ! debug: saves element flags
      if (DEBUG) then
        write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_outer_core_proc',myrank
        call write_VTK_data_elem_l(NSPEC_OUTER_CORE,NGLOB_OUTER_CORE, &
                                   xstore_glob,ystore_glob,zstore_glob, &
                                   ibool,is_on_a_slice_edge,filename)
      endif

      ! added this to reduce the size of the buffers
      ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
      npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_outer_core(:)), &
                                  maxval(npoin2D_eta_outer_core(:)))
    endif

  case (IREGION_INNER_CORE)
    ! inner core
    if (NSPEC_INNER_CORE > 0) then
      ! user output
      if (myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'inner core region:'
        call flush_IMAIN()
      endif
      call cmi_read_buffer_data(IREGION_INNER_CORE, &
                                NGLOB2DMAX_XMIN_XMAX(IREGION_INNER_CORE), &
                                NGLOB2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
                                NGLOB1D_RADIAL(IREGION_INNER_CORE), &
                                iboolleft_xi_inner_core,iboolright_xi_inner_core, &
                                iboolleft_eta_inner_core,iboolright_eta_inner_core, &
                                npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
                                iboolfaces_inner_core,npoin2D_faces_inner_core, &
                                iboolcorner_inner_core)

      ! central cube buffers
      if (INCLUDE_CENTRAL_CUBE) then

        if (myrank == 0) then
          write(IMAIN,*)
          write(IMAIN,*) 'including central cube'
          call flush_IMAIN()
        endif
        call synchronize_all()

        ! allocates boundary indexing arrays for central cube
        allocate(ibelm_xmin_inner_core(NSPEC2DMAX_XMIN_XMAX_IC), &
                 ibelm_xmax_inner_core(NSPEC2DMAX_XMIN_XMAX_IC), &
                 ibelm_ymin_inner_core(NSPEC2DMAX_YMIN_YMAX_IC), &
                 ibelm_ymax_inner_core(NSPEC2DMAX_YMIN_YMAX_IC), &
                 ibelm_top_inner_core(NSPEC2D_TOP_IC), &
                 ibelm_bottom_inner_core(NSPEC2D_BOTTOM_IC),stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error allocating central cube index arrays')

        ! gets coupling arrays for inner core
        nspec2D_xmin_inner_core = nspec2D_xmin
        nspec2D_xmax_inner_core = nspec2D_xmax
        nspec2D_ymin_inner_core = nspec2D_ymin
        nspec2D_ymax_inner_core = nspec2D_ymax

        ibelm_xmin_inner_core(:) = ibelm_xmin(:)
        ibelm_xmax_inner_core(:) = ibelm_xmax(:)
        ibelm_ymin_inner_core(:) = ibelm_ymin(:)
        ibelm_ymax_inner_core(:) = ibelm_ymax(:)
        ibelm_bottom_inner_core(:) = ibelm_bottom(:)
        ibelm_top_inner_core(:) = ibelm_top(:)

        ! compute number of messages to expect in cube as well as their size
        call comp_central_cube_buffer_size(iproc_xi,iproc_eta,ichunk, &
                                           NPROC_XI,NPROC_ETA,NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                                           nb_msgs_theor_in_cube,npoin2D_cube_from_slices)

        ! this value is used for dynamic memory allocation, therefore make sure it is never zero
        if (nb_msgs_theor_in_cube > 0) then
          non_zero_nb_msgs_theor_in_cube = nb_msgs_theor_in_cube
        else
          non_zero_nb_msgs_theor_in_cube = 1
        endif
        if (myrank == 0) then
          write(IMAIN,*) '  number of messages in cube : ',nb_msgs_theor_in_cube
          write(IMAIN,*) '  number of 2D points in cube: ',npoin2D_cube_from_slices
          call flush_IMAIN()
        endif
        call synchronize_all()

        ! allocate buffers for cube and slices
        allocate(sender_from_slices_to_cube(non_zero_nb_msgs_theor_in_cube), &
                 buffer_all_cube_from_slices(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices,NDIM), &
                 buffer_slices(npoin2D_cube_from_slices,NDIM), &
                 buffer_slices2(npoin2D_cube_from_slices,NDIM), &
                 ibool_central_cube(non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices),stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error allocating cube buffers')
        sender_from_slices_to_cube(:) = -1
        ibool_central_cube(:,:) = -1
        buffer_slices(:,:) = 0.d0
        buffer_slices2(:,:) = 0.d0
        buffer_all_cube_from_slices(:,:,:) = 0.d0

        ! handles the communications with the central cube if it was included in the mesh
        ! create buffers to assemble with the central cube
        call create_central_cube_buffers(iproc_xi,iproc_eta,ichunk, &
                                         NPROC_XI,NPROC_ETA,NCHUNKS, &
                                         NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                         NSPEC2DMAX_XMIN_XMAX(IREGION_INNER_CORE),NSPEC2DMAX_YMIN_YMAX(IREGION_INNER_CORE), &
                                         NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                                         addressing,ibool,idoubling, &
                                         xstore_glob,ystore_glob,zstore_glob, &
                                         nspec2D_xmin_inner_core,nspec2D_xmax_inner_core, &
                                         nspec2D_ymin_inner_core,nspec2D_ymax_inner_core, &
                                         ibelm_xmin_inner_core,ibelm_xmax_inner_core, &
                                         ibelm_ymin_inner_core,ibelm_ymax_inner_core,ibelm_bottom_inner_core, &
                                         nb_msgs_theor_in_cube,non_zero_nb_msgs_theor_in_cube,npoin2D_cube_from_slices, &
                                         receiver_cube_from_slices,sender_from_slices_to_cube,ibool_central_cube, &
                                         buffer_slices,buffer_slices2,buffer_all_cube_from_slices)

        if (myrank == 0) then
          write(IMAIN,*) '  creating central cube done'
          write(IMAIN,*)
          call flush_IMAIN()
        endif
        call synchronize_all()

        ! frees memory
        deallocate(ibelm_xmin_inner_core,ibelm_xmax_inner_core)
        deallocate(ibelm_ymin_inner_core,ibelm_ymax_inner_core)
        deallocate(ibelm_top_inner_core)

      else

        ! allocate fictitious buffers for cube and slices with a dummy size
        ! just to be able to use them as arguments in subroutine calls
        allocate(sender_from_slices_to_cube(1), &
                 buffer_all_cube_from_slices(1,1,1), &
                 buffer_slices(1,1), &
                 buffer_slices2(1,1), &
                 ibool_central_cube(1,1),stat=ier)
        if (ier /= 0 ) call exit_MPI(myrank,'Error allocating dummy buffers')

      endif

      ! note: fix_... routines below update is_on_a_slice_edge_.. arrays:
      !          assign flags for each element which is on a rim of the slice
      !          thus, they include elements on top and bottom not shared with other MPI partitions
      !
      !          we will re-set these flags when setting up inner/outer elements, but will
      !          use these arrays for now as initial guess for the search for elements which share a global point
      !          between different MPI processes
      call fix_non_blocking_slices(is_on_a_slice_edge, &
                                   iboolright_xi_inner_core,iboolleft_xi_inner_core, &
                                   iboolright_eta_inner_core,iboolleft_eta_inner_core, &
                                   npoin2D_xi_inner_core,npoin2D_eta_inner_core, &
                                   ibool, &
                                   NSPEC_INNER_CORE,NGLOB_INNER_CORE,NGLOB2DMAX_XMIN_XMAX_IC,NGLOB2DMAX_YMIN_YMAX_IC)

      if (INCLUDE_CENTRAL_CUBE) then
        ! updates flags for elements on slice boundaries
        call fix_non_blocking_central_cube(is_on_a_slice_edge, &
                                           ibool,NSPEC_INNER_CORE,NGLOB_INNER_CORE,nb_msgs_theor_in_cube,ibelm_bottom_inner_core, &
                                           idoubling,npoin2D_cube_from_slices, &
                                           ibool_central_cube,NSPEC2D_BOTTOM(IREGION_INNER_CORE), &
                                           ichunk,NPROC_XI)
      endif

      ! debug: saves element flags
      if (DEBUG) then
        write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_inner_core_proc',myrank
        call write_VTK_data_elem_l(NSPEC_INNER_CORE,NGLOB_INNER_CORE, &
                                   xstore_glob,ystore_glob,zstore_glob, &
                                   ibool,is_on_a_slice_edge,filename)
      endif

      ! added this to reduce the size of the buffers
      ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
      npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_inner_core(:)), &
                                  maxval(npoin2D_eta_inner_core(:)))
    endif

  case (IREGION_TRINFINITE)
    ! transition infinite layer
    if (NSPEC_TRINFINITE > 0) then
      ! user output
      if(myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'transition infinite region:'
        call flush_IMAIN()
      endif
      call cmi_read_buffer_data(IREGION_TRINFINITE, &
                                NGLOB2DMAX_XMIN_XMAX(IREGION_TRINFINITE), &
                                NGLOB2DMAX_YMIN_YMAX(IREGION_TRINFINITE), &
                                NGLOB1D_RADIAL(IREGION_TRINFINITE), &
                                iboolleft_xi_trinfinite,iboolright_xi_trinfinite, &
                                iboolleft_eta_trinfinite,iboolright_eta_trinfinite, &
                                npoin2D_xi_trinfinite,npoin2D_eta_trinfinite, &
                                iboolfaces_trinfinite,npoin2D_faces_trinfinite, &
                                iboolcorner_trinfinite)

      ! note: fix_... routines below update is_on_a_slice_edge_.. arrays:
      !          assign flags for each element which is on a rim of the slice
      !          thus, they include elements on top and bottom not shared with other MPI partitions
      !
      !          we will re-set these flags when setting up inner/outer elements, but will
      !          use these arrays for now as initial guess for the search for elements which share a global point
      !          between different MPI processes
      call fix_non_blocking_slices(is_on_a_slice_edge, &
                                   iboolright_xi_trinfinite,iboolleft_xi_trinfinite, &
                                   iboolright_eta_trinfinite,iboolleft_eta_trinfinite, &
                                   npoin2D_xi_trinfinite,npoin2D_eta_trinfinite, &
                                   ibool, &
                                   NSPEC_TRINFINITE,NGLOB_TRINFINITE,NGLOB2DMAX_XMIN_XMAX_TRINF,NGLOB2DMAX_YMIN_YMAX_TRINF)

      ! debug: saves element flags
      if (DEBUG) then
        write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_trinfinite_proc',myrank
        call write_VTK_data_elem_l(NSPEC_TRINFINITE,NGLOB_TRINFINITE, &
                                   xstore_glob,ystore_glob,zstore_glob, &
                                   ibool,is_on_a_slice_edge,filename)
      endif

      ! added this to reduce the size of the buffers
      ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
      npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_trinfinite(:)), &
                                  maxval(npoin2D_eta_trinfinite(:)))

    endif

  case (IREGION_INFINITE)
    ! infinite layer
    if (NSPEC_INFINITE > 0) then
      ! user output
      if(myrank == 0) then
        write(IMAIN,*)
        write(IMAIN,*) 'infinite region:'
        call flush_IMAIN()
      endif
      call cmi_read_buffer_data(IREGION_INFINITE, &
                                NGLOB2DMAX_XMIN_XMAX(IREGION_INFINITE), &
                                NGLOB2DMAX_YMIN_YMAX(IREGION_INFINITE), &
                                NGLOB1D_RADIAL(IREGION_INFINITE), &
                                iboolleft_xi_infinite,iboolright_xi_infinite, &
                                iboolleft_eta_infinite,iboolright_eta_infinite, &
                                npoin2D_xi_infinite,npoin2D_eta_infinite, &
                                iboolfaces_infinite,npoin2D_faces_infinite, &
                                iboolcorner_infinite)

      ! note: fix_... routines below update is_on_a_slice_edge_.. arrays:
      !          assign flags for each element which is on a rim of the slice
      !          thus, they include elements on top and bottom not shared with other MPI partitions
      !
      !          we will re-set these flags when setting up inner/outer elements, but will
      !          use these arrays for now as initial guess for the search for elements which share a global point
      !          between different MPI processes
      call fix_non_blocking_slices(is_on_a_slice_edge, &
                                   iboolright_xi_infinite,iboolleft_xi_infinite, &
                                   iboolright_eta_infinite,iboolleft_eta_infinite, &
                                   npoin2D_xi_infinite,npoin2D_eta_infinite, &
                                   ibool, &
                                   NSPEC_INFINITE,NGLOB_INFINITE,NGLOB2DMAX_XMIN_XMAX_INF,NGLOB2DMAX_YMIN_YMAX_INF)

      ! debug: saves element flags
      if (DEBUG) then
        write(filename,'(a,i6.6)') trim(OUTPUT_FILES)//'/MPI_is_on_a_slice_edge_infinite_proc',myrank
        call write_VTK_data_elem_l(NSPEC_INFINITE,NGLOB_INFINITE, &
                                   xstore_glob,ystore_glob,zstore_glob, &
                                   ibool,is_on_a_slice_edge,filename)
      endif

      ! added this to reduce the size of the buffers
      ! size of buffers is the sum of two sizes because we handle two regions in the same MPI call
      npoin2D_max_all_CM_IC = max(maxval(npoin2D_xi_infinite(:)), &
                                  maxval(npoin2D_eta_infinite(:)))

    endif

  end select

  end subroutine cmi_get_buffers


!
!-------------------------------------------------------------------------------------------------
!

  subroutine cmi_read_buffer_data(iregion_code, &
                                  NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX, &
                                  NGLOB1D_RADIAL, &
                                  iboolleft_xi_s,iboolright_xi_s, &
                                  iboolleft_eta_s,iboolright_eta_s, &
                                  npoin2D_xi_s,npoin2D_eta_s, &
                                  iboolfaces_s,npoin2D_faces_s, &
                                  iboolcorner_s)

  use meshfem_par, only: &
    myrank,IMAIN,NDIM,NUMFACES_SHARED,NUMCORNERS_SHARED,NPROC_XI,NPROC_ETA

  use MPI_interfaces_par

  implicit none

  integer :: iregion_code

  integer :: NGLOB2DMAX_XMIN_XMAX,NGLOB2DMAX_YMIN_YMAX
  integer :: NGLOB1D_RADIAL

  integer, dimension(NGLOB2DMAX_XMIN_XMAX) :: iboolleft_xi_s,iboolright_xi_s
  integer, dimension(NGLOB2DMAX_YMIN_YMAX) :: iboolleft_eta_s,iboolright_eta_s

  integer, dimension(NB_SQUARE_EDGES_ONEDIR) :: npoin2D_xi_s,npoin2D_eta_s

  integer, dimension(NGLOB2DMAX_XY,NUMFACES_SHARED) :: iboolfaces_s
  integer, dimension(NUMFACES_SHARED) :: npoin2D_faces_s

  integer, dimension(NGLOB1D_RADIAL,NUMCORNERS_SHARED) :: iboolcorner_s

  ! local parameters
  integer :: icount_faces,imsg

  ! gets 2-D arrays
  npoin2D_xi_s(:) = npoin2D_xi_all(:)
  npoin2D_eta_s(:) = npoin2D_eta_all(:)

  ! gets MPI buffers on sides
  iboolleft_xi_s(:) = iboolleft_xi(:)
  iboolright_xi_s(:) = iboolright_xi(:)
  iboolleft_eta_s(:) = iboolleft_eta(:)
  iboolright_eta_s(:) = iboolright_eta(:)

  ! gets corner info
  iboolcorner_s(:,:) = iboolcorner(:,:)

  ! gets face info
  npoin2D_faces_s(:) = npoin2D_faces(:)
  iboolfaces_s(:,:) = iboolfaces(:,:)

  ! checks indirect addressing for each message for faces of the chunks
  ! a given slice can belong to at most two faces
  if (NPROC_XI > 1 .or. NPROC_ETA > 1) then
    ! checks faces
    icount_faces = 0
    do imsg = 1,NUMMSGS_FACES
      if (myrank == iprocfrom_faces(imsg) .or. myrank == iprocto_faces(imsg)) then
        icount_faces = icount_faces + 1

        if (icount_faces > NUMFACES_SHARED) then
          print *,'Error ',myrank,' icount_faces: ',icount_faces,'NUMFACES_SHARED:',NUMFACES_SHARED
          print *,'iregion_code:',iregion_code
          call exit_MPI(myrank,'more than NUMFACES_SHARED faces for this slice')
        endif
        if (icount_faces > 2 .and. (NPROC_XI > 1 .or. NPROC_ETA > 1)) then
          print *,'Error ',myrank,' icount_faces: ',icount_faces,'NPROC_XI:',NPROC_XI,'NPROC_ETA:',NPROC_ETA
          print *,'iregion_code:',iregion_code
          call exit_MPI(myrank,'more than two faces for this slice')
        endif
      endif
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  #max of points in MPI buffers along xi npoin2D_xi = ', &
                                  maxval(npoin2D_xi_s(:))
      write(IMAIN,*) '  #max of array elements transferred npoin2D_xi*NDIM = ', &
                                  maxval(npoin2D_xi_s(:))*NDIM
      write(IMAIN,*)
      write(IMAIN,*) '  #max of points in MPI buffers along eta npoin2D_eta = ', &
                                  maxval(npoin2D_eta_s(:))
      write(IMAIN,*) '  #max of array elements transferred npoin2D_eta*NDIM = ', &
                                  maxval(npoin2D_eta_s(:))*NDIM
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  else
    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  no MPI buffers needed'
      call flush_IMAIN()
    endif
  endif

  end subroutine cmi_read_buffer_data

