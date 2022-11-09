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

  subroutine setup_GLL_points()

  use specfem_par

  implicit none

  ! local parameters
  integer :: i,j,k,ier
  integer, dimension(:),allocatable :: ichunk_slice,iproc_xi_slice,iproc_eta_slice

  ! set up GLL points, weights and derivation matrices
  call define_derivation_matrices(xigll,yigll,zigll,wxgll,wygll,wzgll, &
                                  hprime_xx,hprime_yy,hprime_zz, &
                                  hprimewgll_xx,hprimewgll_yy,hprimewgll_zz, &
                                  wgllwgll_xy,wgllwgll_xz,wgllwgll_yz,wgll_cube)

  ! define a 3D extension in order to be able to force vectorization in the compute_forces_**_Dev routines
  do k = 1,NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        wgllwgll_yz_3D(i,j,k) = wgllwgll_yz(j,k)
        wgllwgll_xz_3D(i,j,k) = wgllwgll_xz(i,k)
        wgllwgll_xy_3D(i,j,k) = wgllwgll_xy(i,j)
      enddo
    enddo
  enddo

  if (USE_DEVILLE_PRODUCTS_VAL) then
    ! check that optimized routines from Deville et al. (2002) can be used
    if (NGLLX /= 5 .or. NGLLY /= 5 .or. NGLLZ /= 5) &
      stop 'Deville et al. (2002) routines can only be used if NGLLX = NGLLY = NGLLZ = 5'

    ! define transpose of derivation matrix
    do j = 1,NGLLY
      do i = 1,NGLLX
        hprime_xxT(j,i) = hprime_xx(i,j)
        hprimewgll_xxT(j,i) = hprimewgll_xx(i,j)
      enddo
    enddo
  endif

  ! addressing
  ! 2-D addressing (needed for Stacey boundaries and regular grid kernels)
  allocate(addressing(NCHUNKS_VAL,0:NPROC_XI_VAL-1,0:NPROC_ETA_VAL-1),stat=ier)
  if (ier /= 0) stop 'Error allocating array addressing'
  addressing(:,:,:) = -1

  allocate(ichunk_slice(0:NPROCTOT_VAL-1), &
           iproc_xi_slice(0:NPROCTOT_VAL-1), &
           iproc_eta_slice(0:NPROCTOT_VAL-1),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays ichunk_slice,...'

  ! safety check
  if (NPROC /= NPROC_ETA_VAL * NPROC_XI_VAL) stop 'Error invalid NPROC not matching NPROC_ETA * NPROC_XI value'

  ! creates global slice addressing for solver
  call create_addressing(NCHUNKS_VAL,NPROC,NPROC_ETA_VAL,NPROC_XI_VAL,NPROCTOT_VAL, &
                         addressing,ichunk_slice,iproc_xi_slice,iproc_eta_slice, &
                         OUTPUT_FILES)

  ! determine chunk number and local slice coordinates using addressing
  ! (needed for Stacey conditions)
  ichunk = ichunk_slice(myrank)

  ! frees temporary arrays
  deallocate(ichunk_slice,iproc_xi_slice,iproc_eta_slice)

  end subroutine setup_GLL_points
