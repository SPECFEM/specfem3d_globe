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

  use constants,only: IMAIN
  use specfem_par,only: myrank

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
  integer :: iphase

  !---- make inv. table ----------------------
  !---- crust mantle : outer elements (iphase=1)
  iphase = 1
  call make_inv_table(iphase,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                      nspec_outer_crust_mantle,phase_ispec_inner_crust_mantle, &
                      ibool_crust_mantle,phase_iglob_crust_mantle, &
                      ibool_inv_tbl_crust_mantle, ibool_inv_st_crust_mantle, &
                      num_globs_crust_mantle)

  !---- crust mantle : inner elements (iphase=2)
  iphase = 2
  call make_inv_table(iphase,NGLOB_CRUST_MANTLE,NSPEC_CRUST_MANTLE, &
                      nspec_inner_crust_mantle,phase_ispec_inner_crust_mantle, &
                      ibool_crust_mantle,phase_iglob_crust_mantle, &
                      ibool_inv_tbl_crust_mantle, ibool_inv_st_crust_mantle, &
                      num_globs_crust_mantle)

  !---- inner core : outer elements (iphase=1)
  iphase = 1
  call make_inv_table(iphase,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                      nspec_outer_inner_core,phase_ispec_inner_inner_core, &
                      ibool_inner_core,phase_iglob_inner_core, &
                      ibool_inv_tbl_inner_core, ibool_inv_st_inner_core, &
                      num_globs_inner_core)

  !---- inner core : inner elements (iphase=2)
  iphase = 2
  call make_inv_table(iphase,NGLOB_INNER_CORE,NSPEC_INNER_CORE, &
                      nspec_inner_inner_core,phase_ispec_inner_inner_core, &
                      ibool_inner_core,phase_iglob_inner_core, &
                      ibool_inv_tbl_inner_core, ibool_inv_st_inner_core, &
                      num_globs_inner_core)

  !---- outer core : outer elements (iphase=1)
  iphase = 1
  call make_inv_table(iphase,NGLOB_OUTER_CORE,NSPEC_OUTER_CORE, &
                      nspec_outer_outer_core,phase_ispec_inner_outer_core, &
                      ibool_outer_core,phase_iglob_outer_core, &
                      ibool_inv_tbl_outer_core, ibool_inv_st_outer_core, &
                      num_globs_outer_core)

  !---- outer core : inner elements (iphase=2)
  iphase = 2
  call make_inv_table(iphase,NGLOB_OUTER_CORE,NSPEC_OUTER_CORE, &
                      nspec_inner_outer_core,phase_ispec_inner_outer_core, &
                      ibool_outer_core,phase_iglob_outer_core, &
                      ibool_inv_tbl_outer_core, ibool_inv_st_outer_core, &
                      num_globs_outer_core)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  inverse table of ibool done"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  contains

    subroutine make_inv_table(iphase,nglob,nspec,phase_nspec,phase_ispec,ibool,phase_iglob, &
                              ibool_inv_tbl,ibool_inv_st,num_globs)

    implicit none

    ! arguments
    integer,intent(in) :: iphase
    integer,intent(in) :: nglob
    integer,intent(in) :: nspec
    integer :: phase_nspec
    integer, dimension(:,:) :: phase_ispec
    integer, dimension(:,:,:,:) :: ibool
    integer, dimension(:,:) :: phase_iglob
    integer, dimension(:,:),intent(out) :: ibool_inv_tbl
    integer, dimension(:,:),intent(out) :: ibool_inv_st
    integer, dimension(:) :: num_globs

    ! local parameters
    integer, dimension(:),   allocatable :: ibool_inv_num
    integer, dimension(:,:), allocatable :: ibool_inv_tbl_tmp
    integer :: num_alloc_ibool_inv_tbl
    integer :: num_used_ibool_inv_tbl
    integer :: ip, iglob, ispec_p, ispec, iglob_p, ier
    integer :: inum
#ifdef FORCE_VECTORIZATION
    integer :: ijk
#else
    integer :: i,j,k
#endif

    ! tolerance number of shared degrees per node
    integer, parameter :: N_TOL = 20

    ! allocates temporary arrays
    allocate(ibool_inv_num(nglob),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_inv_num array'

    ! number of shared degrees per node
    num_alloc_ibool_inv_tbl = N_TOL*(NGLLX*NGLLY*NGLLZ*nspec/nglob+1)

    allocate(ibool_inv_tbl_tmp(num_alloc_ibool_inv_tbl,nglob),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_inv_tbl_tmp array'


    !---- make temporary array of inv. table : ibool_inv_tbl_tmp
    do iglob = 1, nglob
      ibool_inv_num(iglob) = 0
    enddo

    do ispec_p = 1,phase_nspec
      ispec = phase_ispec(ispec_p,iphase)

      DO_LOOP_IJK

        iglob = ibool(INDEX_IJK,ispec)
        ! increases counter
        ibool_inv_num(iglob) = ibool_inv_num(iglob) + 1
        ! inverse table
#ifdef FORCE_VECTORIZATION
        inum = ijk
#else
        inum = i + (j-1)*NGLLY + (k-1)*NGLLY*NGLLZ
#endif
        ibool_inv_tbl_tmp(ibool_inv_num(iglob),iglob) = inum + NGLLX*NGLLY*NGLLZ*(ispec-1)

      ENDDO_LOOP_IJK

    enddo

    !---- packing : ibool_inv_tbl_tmp -> ibool_inv_tbl
    ip    = 0
    iglob_p = 0
    num_used_ibool_inv_tbl = 0
    do iglob = 1, nglob
      if( ibool_inv_num(iglob) /= 0 ) then
        iglob_p = iglob_p + 1
        phase_iglob(iglob_p,iphase) = iglob
        ibool_inv_st(iglob_p,iphase) = ip + 1
        if( ibool_inv_num(iglob) > num_used_ibool_inv_tbl ) then
           num_used_ibool_inv_tbl = ibool_inv_num(iglob)
        endif
        do inum = 1, ibool_inv_num(iglob)
          ip = ip + 1
          ibool_inv_tbl(ip,iphase) = ibool_inv_tbl_tmp(inum,iglob)
        enddo
      endif
    enddo
    ibool_inv_st(iglob_p+1,iphase) = ip + 1
    num_globs(iphase) = iglob_p

    ! checks
    if( num_used_ibool_inv_tbl > num_alloc_ibool_inv_tbl )  then
      write(IMAIN,*) "Error invalid inverse table setting:"
      write(IMAIN,*) "   num_alloc_ibool_inv_tbl = ",num_alloc_ibool_inv_tbl
      write(IMAIN,*) "   num_used_ibool_inv_tbl  = ",num_used_ibool_inv_tbl
      write(IMAIN,*) "#### num_used_ibool_inv_tbl > num_alloc_ibool_inv_tbl  #####"
      write(IMAIN,*) "#### Program stopped ##########"
      call flush_IMAIN()
    endif

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

  implicit none

  ! local parameters
  integer :: ispec

#ifdef FORCE_VECTORIZATION
  integer :: ijk
#else
  integer :: i,j,k
#endif

  !---- fused array of mapping matrix ----------------------
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

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)"  fused array done"
    call flush_IMAIN()
  endif

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_fused_array

