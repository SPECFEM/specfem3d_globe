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

! we switch between vectorized and non-vectorized version by using pre-processor flag FORCE_VECTORIZATION
! and macros INDEX_IJK, DO_LOOP_IJK, ENDDO_LOOP_IJK defined in config.fh
#include "config.fh"


  subroutine prepare_optimized_arrays()

! optimizes array memory layout to increase computational efficiency

  use constants, only: IMAIN
  use specfem_par, only: myrank

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) "preparing optimized arrays"
    call flush_IMAIN()
  endif

  ! precomputes inverse table of ibool
  call prepare_timerun_ibool_inv_tbl()

  ! prepare fused array for computational kernel
  call prepare_fused_array()

  ! check OpenMP support
  call prepare_openmp()

#ifdef XSMM
  ! prepares LIBXSMM small matrix multiplication functions
  call prepare_xsmm()
#endif

  end subroutine prepare_optimized_arrays

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_timerun_ibool_inv_tbl()

! precomputes inverse table of ibool

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: iphase,ier
  integer :: num_elements

  ! inverse arrays use 1D indexing for better compiler vectorization
  ! only used for Deville routines and FORCE_VECTORIZATION)

  ! checks if used
  if (USE_DEVILLE_PRODUCTS_VAL) then
#ifdef FORCE_VECTORIZATION
    use_inversed_arrays = .true.
#else
    use_inversed_arrays = .false.
#endif
  else
    use_inversed_arrays = .false.
  endif

  ! uses local array to store all element contributions
  if (USE_DEVILLE_PRODUCTS_VAL) then
    ! note: we use allocate for sum_terms arrays rather than defining within subroutine compute_forces_**_Dev() itself
    !       as it will crash when using OpenMP and operating systems with small stack sizes
    !       e.g. see http://stackoverflow.com/questions/22649827/illegal-instruction-error-when-running-openmp-in-gfortran-mac
    allocate(sum_terms_crust_mantle(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE), &
             sum_terms_inner_core(NDIM,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE), &
             sum_terms_outer_core(NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE), &
             stat=ier)
    if (ier /= 0) stop 'Error allocating sum_terms arrays'
    sum_terms_crust_mantle(:,:,:,:,:) = 0._CUSTOM_REAL
    sum_terms_inner_core(:,:,:,:,:) = 0._CUSTOM_REAL
    sum_terms_outer_core(:,:,:,:) = 0._CUSTOM_REAL
  endif

  ! inverse table
  ! this helps to speedup the assembly, especially with OpenMP (or on MIC) threading
  if (use_inversed_arrays) then
    ! allocating arrays
    allocate(ibool_inv_tbl_crust_mantle(NGLLX*NGLLY*NGLLZ*NSPEC_CRUST_MANTLE,2), &
             ibool_inv_tbl_inner_core(NGLLX*NGLLY*NGLLZ*NSPEC_INNER_CORE,2), &
             ibool_inv_tbl_outer_core(NGLLX*NGLLY*NGLLZ*NSPEC_OUTER_CORE,2),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_inv_tbl arrays'

    allocate(ibool_inv_st_crust_mantle(NGLOB_CRUST_MANTLE+1,2), &
             ibool_inv_st_inner_core(NGLOB_INNER_CORE+1,2), &
             ibool_inv_st_outer_core(NGLOB_OUTER_CORE+1,2),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_inv_st arrays'

    allocate(phase_iglob_crust_mantle(NGLOB_CRUST_MANTLE,2), &
             phase_iglob_inner_core(NGLOB_INNER_CORE,2), &
             phase_iglob_outer_core(NGLOB_OUTER_CORE,2),stat=ier)
    if (ier /= 0) stop 'Error allocating phase_iglob arrays'

    ! initializing
    num_globs_crust_mantle(:) = 0
    num_globs_inner_core(:) = 0
    num_globs_outer_core(:) = 0

    ibool_inv_tbl_crust_mantle(:,:) = 0
    ibool_inv_tbl_inner_core(:,:) = 0
    ibool_inv_tbl_outer_core(:,:) = 0

    ibool_inv_st_crust_mantle(:,:) = 0
    ibool_inv_st_inner_core(:,:) = 0
    ibool_inv_st_outer_core(:,:) = 0

    phase_iglob_crust_mantle(:,:) = 0
    phase_iglob_inner_core(:,:) = 0
    phase_iglob_outer_core(:,:) = 0

    !---- make inv. table ----------------------
    ! loops over phases
    ! (1 == outer elements / 2 == inner elements)
    do iphase = 1,2
      ! crust mantle
      if (iphase == 1) then
        ! outer elements (iphase=1)
        num_elements = nspec_outer_crust_mantle
      else
        ! inner elements (iphase=2)
        num_elements = nspec_inner_crust_mantle
      endif
      call make_inv_table(iphase,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                          num_elements,phase_ispec_inner_crust_mantle, &
                          ibool_crust_mantle,phase_iglob_crust_mantle, &
                          ibool_inv_tbl_crust_mantle, ibool_inv_st_crust_mantle, &
                          num_globs_crust_mantle)

      ! inner core
      if (iphase == 1) then
        ! outer elements (iphase=1)
        num_elements = nspec_outer_inner_core
      else
        ! inner elements (iphase=2)
        num_elements = nspec_inner_inner_core
      endif
      call make_inv_table(iphase,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                          num_elements,phase_ispec_inner_inner_core, &
                          ibool_inner_core,phase_iglob_inner_core, &
                          ibool_inv_tbl_inner_core, ibool_inv_st_inner_core, &
                          num_globs_inner_core,idoubling_inner_core)

      ! outer core
      if (iphase == 1) then
        ! outer elements (iphase=1)
        num_elements = nspec_outer_outer_core
      else
        ! inner elements (iphase=2)
        num_elements = nspec_inner_outer_core
      endif
      call make_inv_table(iphase,NGLOB_OUTER_CORE,NSPEC_OUTER_CORE, &
                          num_elements,phase_ispec_inner_outer_core, &
                          ibool_outer_core,phase_iglob_outer_core, &
                          ibool_inv_tbl_outer_core, ibool_inv_st_outer_core, &
                          num_globs_outer_core)
    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)"  inverse table of ibool done"
      call flush_IMAIN()
    endif
  endif

  ! synchronizes processes
  call synchronize_all()

  contains

    subroutine make_inv_table(iphase,nglob,nspec, &
                              phase_nspec,phase_ispec,ibool,phase_iglob, &
                              ibool_inv_tbl,ibool_inv_st,num_globs,idoubling)

    implicit none

    ! arguments
    integer,intent(in) :: iphase
    integer,intent(in) :: nglob
    integer,intent(in) :: nspec
    integer,intent(in) :: phase_nspec
    integer, dimension(:,:),intent(in) :: phase_ispec
    integer, dimension(:,:,:,:),intent(in) :: ibool

    integer, dimension(:,:),intent(inout) :: phase_iglob
    integer, dimension(:,:),intent(inout) :: ibool_inv_tbl
    integer, dimension(:,:),intent(inout) :: ibool_inv_st
    integer, dimension(:),intent(inout) :: num_globs

    integer,dimension(:),optional :: idoubling

    ! local parameters
    integer, dimension(:),   allocatable :: ibool_inv_num
    integer, dimension(:,:), allocatable :: ibool_inv_tbl_tmp
    integer :: num_alloc_ibool_inv_tbl,num_alloc_ibool_inv_tbl_theor
    integer :: num_used_ibool_inv_tbl
    integer :: ip, iglob, ispec_p, ispec, iglob_p, ier
    integer :: inum
#ifdef FORCE_VECTORIZATION
    integer :: ijk
#else
    integer :: i,j,k
#endif
    logical :: is_inner_core

    ! tolerance number of shared degrees per node
    integer, parameter :: N_TOL = 20

    ! checks if anything to do (e.g., no outer elements for single process simulations)
    if (phase_nspec == 0) return

    ! checks if inner core region
    if (present(idoubling)) then
      is_inner_core = .true.
    else
      is_inner_core = .false.
    endif

    ! allocates temporary arrays
    allocate(ibool_inv_num(nglob),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_inv_num array'

    ! gets valence of global degrees of freedom for current phase (inner/outer) elements
    ibool_inv_num(:) = 0
    do ispec_p = 1,phase_nspec
      ispec = phase_ispec(ispec_p,iphase)

      ! exclude fictitious elements in central cube
      if (is_inner_core) then
        if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
      endif

      DO_LOOP_IJK
        iglob = ibool(INDEX_IJK,ispec)
        ! increases valence counter
        ibool_inv_num(iglob) = ibool_inv_num(iglob) + 1
      ENDDO_LOOP_IJK
    enddo

    ! gets maximum valence value
    num_alloc_ibool_inv_tbl = maxval(ibool_inv_num(:))

    ! theoretical number of maximum shared degrees per node
    num_alloc_ibool_inv_tbl_theor = N_TOL*(NGLLX*NGLLY*NGLLZ*nspec/nglob+1)

    ! checks valence
    if (num_alloc_ibool_inv_tbl < 1 .or. num_alloc_ibool_inv_tbl > num_alloc_ibool_inv_tbl_theor) then
      print *,'Error invalid maximum valence:'
      print *,'valence value = ',num_alloc_ibool_inv_tbl,' - theoretical maximum = ',num_alloc_ibool_inv_tbl_theor
      stop 'Error invalid maximum valence value'
    endif
    ! debug
    !print *,myrank,'maximum shared degrees theoretical = ',num_alloc_ibool_inv_tbl_theor ! regional_Greece_small example: 40
    !print *,myrank,'maximum shared degrees from array  = ',maxval(ibool_inv_num(:))      ! regional_Greece_small example: 8 and 16

    allocate(ibool_inv_tbl_tmp(num_alloc_ibool_inv_tbl,nglob),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_inv_tbl_tmp array'

    !---- make temporary array of inv. table : ibool_inv_tbl_tmp
    ibool_inv_tbl_tmp(:,:) = 0
    ibool_inv_num(:) = 0
    do ispec_p = 1,phase_nspec
      ispec = phase_ispec(ispec_p,iphase)

      ! exclude fictitious elements in central cube
      if (is_inner_core) then
        if (idoubling(ispec) == IFLAG_IN_FICTITIOUS_CUBE) cycle
      endif

      DO_LOOP_IJK

        iglob = ibool(INDEX_IJK,ispec)

        ! increases counter
        ibool_inv_num(iglob) = ibool_inv_num(iglob) + 1

        ! inverse table
        ! sets 1D index of local GLL point (between 1 and NGLLCUBE)
#ifdef FORCE_VECTORIZATION
        inum = ijk
#else
        inum = i + (j-1)*NGLLY + (k-1)*NGLLY*NGLLZ
#endif
        ! sets 1D index in local ibool array
        ibool_inv_tbl_tmp(ibool_inv_num(iglob),iglob) = inum + NGLLX*NGLLY*NGLLZ*(ispec-1)

      ENDDO_LOOP_IJK

    enddo

    !---- packing : ibool_inv_tbl_tmp -> ibool_inv_tbl
    ip = 0
    iglob_p = 0
    num_used_ibool_inv_tbl = 0
    do iglob = 1, nglob
      if (ibool_inv_num(iglob) /= 0) then
        iglob_p = iglob_p + 1

        phase_iglob(iglob_p,iphase) = iglob

        ! sets start index of table entry for this global node
        ibool_inv_st(iglob_p,iphase) = ip + 1

        ! sets maximum of used valence
        if ( ibool_inv_num(iglob) > num_used_ibool_inv_tbl ) num_used_ibool_inv_tbl = ibool_inv_num(iglob)

        ! loops over valence
        do inum = 1, ibool_inv_num(iglob)
          ! increases total counter
          ip = ip + 1
          ! maps local 1D index in ibool array
          ibool_inv_tbl(ip,iphase) = ibool_inv_tbl_tmp(inum,iglob)
        enddo
      endif
    enddo
    ! sets last entry in start index table
    ibool_inv_st(iglob_p+1,iphase) = ip + 1

    ! total number global nodes in this phase (inner/outer)
    num_globs(iphase) = iglob_p

    ! checks
    if ( num_used_ibool_inv_tbl > num_alloc_ibool_inv_tbl ) then
      print *,"Error invalid inverse table setting:"
      print *,"  num_alloc_ibool_inv_tbl = ",num_alloc_ibool_inv_tbl
      print *,"  num_used_ibool_inv_tbl  = ",num_used_ibool_inv_tbl
      print *,"invalid value encountered: num_used_ibool_inv_tbl > num_alloc_ibool_inv_tbl"
      print *,"#### Program exits... ##########"
      call exit_MPI(myrank,'Error making inverse table for optimized arrays')
    endif

    ! debug
    !if (myrank == 0) then
    !  print *,'ibool_inv_tbl: '
    !  do iglob_p = 1,200
    !    print *,'  ',iglob_p,'table = ',(ibool_inv_tbl(ip,iphase), &
    !                                     ip = ibool_inv_st(iglob_p,iphase),ibool_inv_st(iglob_p+1,iphase)-1)
    !  enddo
    !endif

    ! frees memory
    deallocate(ibool_inv_num)
    deallocate(ibool_inv_tbl_tmp)

    end subroutine make_inv_table

  end subroutine prepare_timerun_ibool_inv_tbl

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_fused_array()

! prepare fused array for computational kernel
!
! note: fusing separate arrays for xi/eta/gamma into a single one increases efficiency for hardware pre-fetching

  use specfem_par
  use specfem_par_crustmantle
  use specfem_par_innercore
  use specfem_par_outercore

  implicit none

  ! local parameters
  integer :: ispec,ier

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  ! fused array only needed for compute forces in crust/mantle (Deville routine)
  if (USE_DEVILLE_PRODUCTS_VAL) then

    ! crust/mantle
    ! allocates fused array
    allocate(deriv_mapping_crust_mantle(9,NGLLX,NGLLY,NGLLZ,NSPEC_CRUST_MANTLE),stat=ier)
    if (ier /= 0) stop 'Error allocating array deriv_mapping_crust_mantle'

    !---- fused array of mapping matrix
    do ispec = 1,NSPEC_CRUST_MANTLE

      DO_LOOP_IJK

        deriv_mapping_crust_mantle(1,INDEX_IJK,ispec) = xix_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(2,INDEX_IJK,ispec) = xiy_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(3,INDEX_IJK,ispec) = xiz_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(4,INDEX_IJK,ispec) = etax_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(5,INDEX_IJK,ispec) = etay_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(6,INDEX_IJK,ispec) = etaz_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(7,INDEX_IJK,ispec) = gammax_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(8,INDEX_IJK,ispec) = gammay_crust_mantle(INDEX_IJK,ispec)
        deriv_mapping_crust_mantle(9,INDEX_IJK,ispec) = gammaz_crust_mantle(INDEX_IJK,ispec)

      ENDDO_LOOP_IJK

    enddo

    ! inner core
    ! allocates fused array
    allocate(deriv_mapping_inner_core(9,NGLLX,NGLLY,NGLLZ,NSPEC_INNER_CORE),stat=ier)
    if (ier /= 0) stop 'Error allocating array deriv_mapping_inner_core'

    !---- fused array of mapping matrix
    do ispec = 1,NSPEC_INNER_CORE

      DO_LOOP_IJK

        deriv_mapping_inner_core(1,INDEX_IJK,ispec) = xix_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(2,INDEX_IJK,ispec) = xiy_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(3,INDEX_IJK,ispec) = xiz_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(4,INDEX_IJK,ispec) = etax_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(5,INDEX_IJK,ispec) = etay_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(6,INDEX_IJK,ispec) = etaz_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(7,INDEX_IJK,ispec) = gammax_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(8,INDEX_IJK,ispec) = gammay_inner_core(INDEX_IJK,ispec)
        deriv_mapping_inner_core(9,INDEX_IJK,ispec) = gammaz_inner_core(INDEX_IJK,ispec)

      ENDDO_LOOP_IJK

    enddo

    ! outer core
    ! allocates fused array
    allocate(deriv_mapping_outer_core(9,NGLLX,NGLLY,NGLLZ,NSPEC_OUTER_CORE),stat=ier)
    if (ier /= 0) stop 'Error allocating array deriv_mapping_outer_core'

    !---- fused array of mapping matrix
    do ispec = 1,NSPEC_OUTER_CORE

      DO_LOOP_IJK

        deriv_mapping_outer_core(1,INDEX_IJK,ispec) = xix_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(2,INDEX_IJK,ispec) = xiy_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(3,INDEX_IJK,ispec) = xiz_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(4,INDEX_IJK,ispec) = etax_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(5,INDEX_IJK,ispec) = etay_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(6,INDEX_IJK,ispec) = etaz_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(7,INDEX_IJK,ispec) = gammax_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(8,INDEX_IJK,ispec) = gammay_outer_core(INDEX_IJK,ispec)
        deriv_mapping_outer_core(9,INDEX_IJK,ispec) = gammaz_outer_core(INDEX_IJK,ispec)

      ENDDO_LOOP_IJK

    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)"  fused array done"
      call flush_IMAIN()
    endif

    ! synchronizes processes
    call synchronize_all()
  endif

  end subroutine prepare_fused_array

!
!-------------------------------------------------------------------------------------------------
!

#ifdef XSMM
  subroutine prepare_xsmm()

  use constants,only: CUSTOM_REAL,SIZE_DOUBLE,m1,m2,IMAIN

  use specfem_par,only: myrank

  use my_libxsmm,only: libxsmm_init,libxsmm_smm_25_5_5,libxsmm_smm_5_25_5,libxsmm_smm_5_5_5

  implicit none

  ! temporary arrays
  real(kind=CUSTOM_REAL),dimension(5,5) :: A1
  real(kind=CUSTOM_REAL),dimension(5,25) :: B1
  real(kind=CUSTOM_REAL),dimension(5,25) :: C1

  real(kind=CUSTOM_REAL),dimension(25,5) :: A2
  real(kind=CUSTOM_REAL),dimension(5,5) :: B2
  real(kind=CUSTOM_REAL),dimension(25,5) :: C2

  real(kind=CUSTOM_REAL),dimension(5,5,5) :: A3
  real(kind=CUSTOM_REAL),dimension(5,5) :: B3
  real(kind=CUSTOM_REAL),dimension(5,5,5) :: C3

  ! quick check
  if (m1 /= 5) stop 'LibXSMM with invalid m1 constant (must have m1 == 5)'
  if (m2 /= 5*5) stop 'LibXSMM with invalid m2 constant (must have m2 == 5*5)'
  if (CUSTOM_REAL == SIZE_DOUBLE) stop 'LibXSMM optimization only for single precision functions'

  ! initializes LIBXSMM
  call libxsmm_init()

  ! LIBXSMM static functions
  ! use version compilation with: MNK="5 25, 5" ALPHA=1 BETA=0

  ! dummy static calls to check if they work...
  ! (see in compute_forces_**Dev.F90 routines for actual function call)

  ! with A(n1,n2) 5x5-matrix, B(n2,n3) 5x25-matrix and C(n1,n3) 5x25-matrix
  call libxsmm_smm_5_25_5(a=A1, b=B1, c=C1, pa=A1, pb=B1, pc=C1)

  ! with A(n1,n2) 25x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3) 25x5-matrix
  call libxsmm_smm_25_5_5(a=A2, b=B2, c=C2, pa=A2, pb=B2, pc=C2)

  ! with A(n1,n2,n4) 5x5x5-matrix, B(n2,n3) 5x5-matrix and C(n1,n3,n4) 5x5x5-matrix
  call libxsmm_smm_5_5_5(a=A3(1,1,1), b=B3, c=C3(1,1,1),pa=A3(1,1,1), pb=B3, pc=C3(1,1,1))

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) "  LIBXSMM functions ready for small matrix-matrix multiplications"
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_xsmm
#endif
