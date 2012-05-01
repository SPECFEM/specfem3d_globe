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

  subroutine setup_color_perm(myrank,iregion_code,nspec,nglob, &
                              ibool,is_on_a_slice_edge,prname, &
                              npoin2D_xi,npoin2D_eta)

  use constants
  use meshfem3D_par,only: NSTEP,DT,NPROC_XI,NPROC_ETA
  implicit none

  ! standard include of the MPI library
  include 'mpif.h'

  integer :: myrank
  integer :: iregion_code

  integer :: nspec,nglob
  integer, dimension(NGLLX,NGLLY,NGLLZ,nspec) :: ibool

  ! this for non blocking MPI
  logical, dimension(nspec) :: is_on_a_slice_edge

  ! name of the database file
  character(len=150) :: prname

  integer :: npoin2D_xi,npoin2D_eta

  ! local parameters
  integer :: nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer
  integer, dimension(:), allocatable :: perm
  integer, dimension(:), allocatable :: first_elem_number_in_this_color
  integer, dimension(:), allocatable :: num_of_elems_in_this_color

  integer :: icolor,ispec_counter
  integer :: nspec_outer_min_global,nspec_outer_max_global
  integer :: ispec,ier

  !!!! David Michea: detection of the edges, coloring and permutation separately
  allocate(perm(nspec))

  ! implement mesh coloring for GPUs if needed, to create subsets of disconnected elements
  ! to remove dependencies and the need for atomic operations in the sum of elemental contributions in the solver
  if(USE_MESH_COLORING_GPU) then

    ! user output
    if(myrank == 0 ) write(IMAIN,*) '  creating mesh coloring'

    allocate(first_elem_number_in_this_color(MAX_NUMBER_OF_COLORS + 1))

    call get_perm_color_faster(is_on_a_slice_edge,ibool,perm,nspec,nglob, &
                              nb_colors_outer_elements,nb_colors_inner_elements,nspec_outer, &
                              first_elem_number_in_this_color,myrank)

    ! for the last color, the next color is fictitious and its first (fictitious) element number is nspec + 1
    first_elem_number_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements + 1) = nspec + 1

    allocate(num_of_elems_in_this_color(nb_colors_outer_elements + nb_colors_inner_elements))

    ! save mesh coloring
    open(unit=99,file=prname(1:len_trim(prname))//'num_of_elems_in_this_color.dat', &
         status='unknown',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening num_of_elems_in_this_color file')

    ! number of colors for outer elements
    write(99,*) nb_colors_outer_elements

    ! number of colors for inner elements
    write(99,*) nb_colors_inner_elements

    ! number of elements in each color
    do icolor = 1, nb_colors_outer_elements + nb_colors_inner_elements
      num_of_elems_in_this_color(icolor) = first_elem_number_in_this_color(icolor+1) &
                                          - first_elem_number_in_this_color(icolor)
      write(99,*) num_of_elems_in_this_color(icolor)
    enddo
    close(99)

    ! check that the sum of all the numbers of elements found in each color is equal
    ! to the total number of elements in the mesh
    if(sum(num_of_elems_in_this_color) /= nspec) then
      print *,'nspec = ',nspec
      print *,'total number of elements in all the colors of the mesh = ',sum(num_of_elems_in_this_color)
      call exit_mpi(myrank,'incorrect total number of elements in all the colors of the mesh')
    endif

    ! check that the sum of all the numbers of elements found in each color for the outer elements is equal
    ! to the total number of outer elements found in the mesh
    if(sum(num_of_elems_in_this_color(1:nb_colors_outer_elements)) /= nspec_outer) then
      print *,'nspec_outer = ',nspec_outer
      print *,'total number of elements in all the colors of the mesh for outer elements = ', &
        sum(num_of_elems_in_this_color)
      call exit_mpi(myrank,'incorrect total number of elements in all the colors of the mesh for outer elements')
    endif

    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

    deallocate(first_elem_number_in_this_color)
    deallocate(num_of_elems_in_this_color)

  else

    !! DK DK for regular C + MPI version for CPUs: do not use colors but nonetheless put all the outer elements
    !! DK DK first in order to be able to overlap non-blocking MPI communications with calculations

    !! DK DK nov 2010, for Rosa Badia / StarSs:
    !! no need for mesh coloring, but need to implement inner/outer subsets for non blocking MPI for StarSs
    ispec_counter = 0
    perm(:) = 0

    ! first generate all the outer elements
    do ispec = 1,nspec
      if(is_on_a_slice_edge(ispec)) then
        ispec_counter = ispec_counter + 1
        perm(ispec) = ispec_counter
      endif
    enddo

    ! make sure we have detected some outer elements
    if(ispec_counter <= 0) stop 'fatal error: no outer elements detected!'

    ! store total number of outer elements
    nspec_outer = ispec_counter

    ! then generate all the inner elements
    do ispec = 1,nspec
      if(.not. is_on_a_slice_edge(ispec)) then
        ispec_counter = ispec_counter + 1
        perm(ispec) = ispec_counter
      endif
    enddo

    ! test that all the elements have been used once and only once
    if(ispec_counter /= nspec) stop 'fatal error: ispec_counter not equal to nspec'

    ! do basic checks
    if(minval(perm) /= 1) stop 'minval(perm) should be 1'
    if(maxval(perm) /= nspec) stop 'maxval(perm) should be nspec'

    call MPI_ALLREDUCE(nspec_outer,nspec_outer_min_global,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ier)
    call MPI_ALLREDUCE(nspec_outer,nspec_outer_max_global,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  endif ! USE_MESH_COLORING_GPU

  !! DK DK and Manh Ha, Nov 2011: added this to use the new mesher in the CUDA or C / StarSs test codes

  if (myrank == 0 .and. iregion_code == IREGION_CRUST_MANTLE) then
    ! write a header file for the Fortran version of the solver
    open(unit=99,file=prname(1:len_trim(prname))//'values_from_mesher_f90.h', &
          status='unknown',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening file values_from_mesher_f90.h')

    write(99,*) 'integer, parameter :: NSPEC = ',nspec
    write(99,*) 'integer, parameter :: NGLOB = ',nglob
    !!! DK DK use 1000 time steps only for the scaling tests
    write(99,*) 'integer, parameter :: NSTEP = 1000 !!!!!!!!!!! ',nstep
    write(99,*) 'real(kind=4), parameter :: deltat = ',DT
    write(99,*)
    write(99,*) 'integer, parameter ::  NGLOB2DMAX_XMIN_XMAX = ',npoin2D_xi
    write(99,*) 'integer, parameter ::  NGLOB2DMAX_YMIN_YMAX = ',npoin2D_eta
    write(99,*) 'integer, parameter ::  NGLOB2DMAX_ALL = ',max(npoin2D_xi,npoin2D_eta)
    write(99,*) 'integer, parameter ::  NPROC_XI = ',NPROC_XI
    write(99,*) 'integer, parameter ::  NPROC_ETA = ',NPROC_ETA
    write(99,*)
    write(99,*) '! element number of the source and of the station'
    write(99,*) '! after permutation of the elements by mesh coloring'
    write(99,*) '! and inner/outer set splitting in the mesher'
    write(99,*) 'integer, parameter :: NSPEC_SOURCE = ',perm(NSPEC/3)
    write(99,*) 'integer, parameter :: RANK_SOURCE = 0'
    write(99,*)
    write(99,*) 'integer, parameter :: RANK_STATION = (NPROC_XI*NPROC_ETA - 1)'
    write(99,*) 'integer, parameter :: NSPEC_STATION = ',perm(2*NSPEC/3)

    ! save coordinates of the seismic source
    !   write(99,*) xstore(2,2,2,10);
    !   write(99,*) ystore(2,2,2,10);
    !   write(99,*) zstore(2,2,2,10);

    ! save coordinates of the seismic station
    !   write(99,*) xstore(2,2,2,nspec-10);
    !   write(99,*) ystore(2,2,2,nspec-10);
    !   write(99,*) zstore(2,2,2,nspec-10);
    close(99)

    !! write a header file for the C version of the solver
    open(unit=99,file=prname(1:len_trim(prname))//'values_from_mesher_C.h', &
          status='unknown',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening file values_from_mesher_C.h')

    write(99,*) '#define NSPEC ',nspec
    write(99,*) '#define NGLOB ',nglob
    !!    write(99,*) '#define NSTEP ',nstep
    !!! DK DK use 1000 time steps only for the scaling tests
    write(99,*) '// #define NSTEP ',nstep
    write(99,*) '#define NSTEP 1000'
    ! put an "f" at the end to force single precision
    write(99,"('#define deltat ',e18.10,'f')") DT
    write(99,*) '#define NGLOB2DMAX_XMIN_XMAX ',npoin2D_xi
    write(99,*) '#define NGLOB2DMAX_YMIN_YMAX ',npoin2D_eta
    write(99,*) '#define NGLOB2DMAX_ALL ',max(npoin2D_xi,npoin2D_eta)
    write(99,*) '#define NPROC_XI ',NPROC_XI
    write(99,*) '#define NPROC_ETA ',NPROC_ETA
    write(99,*)
    write(99,*) '// element and MPI slice number of the source and the station'
    write(99,*) '// after permutation of the elements by mesh coloring'
    write(99,*) '// and inner/outer set splitting in the mesher'
    write(99,*) '#define RANK_SOURCE 0'
    write(99,*) '#define NSPEC_SOURCE ',perm(NSPEC/3)
    write(99,*)
    write(99,*) '#define RANK_STATION (NPROC_XI*NPROC_ETA - 1)'
    write(99,*) '#define NSPEC_STATION ',perm(2*NSPEC/3)
    close(99)

    open(unit=99,file=prname(1:len_trim(prname))//'values_from_mesher_nspec_outer.h', &
          status='unknown',iostat=ier)
    if( ier /= 0 ) call exit_mpi(myrank,'error opening values_from_mesher_nspec_outer.h file')

    write(99,*) '#define NSPEC_OUTER ',nspec_outer_max_global
    write(99,*) '// NSPEC_OUTER_min = ',nspec_outer_min_global
    write(99,*) '// NSPEC_OUTER_max = ',nspec_outer_max_global
    close(99)

  endif

  !! DK DK and Manh Ha, Nov 2011: added this to use the new mesher in the CUDA or C / StarSs test codes

  deallocate(perm)


  end subroutine setup_color_perm
