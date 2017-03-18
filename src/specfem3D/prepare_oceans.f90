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


  subroutine prepare_oceans()

! prepares arrays for coupling with oceans
!
! note: handling of coupling on GPU needs to avoid using a mutex to update acceleration;
!       tests so far have shown, that with a simple mutex implementation
!       the results differ between successive runs (probably still due to some race conditions?)
!
!       here we now totally avoid mutex usage and still update each global point only once,
!       also facilitating vectorization of the updating loop

  use specfem_par
  use specfem_par_crustmantle

  implicit none
  ! local parameters
  integer :: ipoin,ispec,ispec2D,i,j,k,iglob,ier
  ! flag to mask ocean-bottom degrees of freedom for ocean load
  logical, dimension(:), allocatable :: updated_dof_ocean_load

  ! checks if anything to do
  if (.not. OCEANS_VAL) return

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "preparing oceans arrays"
    call flush_IMAIN()
  endif

  allocate(updated_dof_ocean_load(NGLOB_CRUST_MANTLE_OCEANS),stat=ier)
  if (ier /= 0) stop 'Error allocating arrays updated_dof_ocean_load'
  updated_dof_ocean_load(:) = .false.

  ! counts global points on surface to oceans
  ipoin = 0
  do ispec2D = 1,NSPEC_TOP
    ispec = ibelm_top_crust_mantle(ispec2D)
    k = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ! get global point number
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ! updates once
        if (.not. updated_dof_ocean_load(iglob)) then
          ipoin = ipoin + 1
          updated_dof_ocean_load(iglob) = .true.
        endif
      enddo
    enddo
  enddo
  npoin_oceans = ipoin

  ! user info
  if (myrank == 0) then
    write(IMAIN,*) "  number of global points on oceans = ",npoin_oceans
    call flush_IMAIN()
  endif

  ! allocates arrays with all global points on ocean surface
  allocate(ibool_ocean_load(npoin_oceans), &
           normal_ocean_load(NDIM,npoin_oceans), &
           rmass_ocean_load_selected(npoin_oceans), &
           stat=ier)
  if (ier /= 0 ) stop 'Error allocating oceans arrays'
  ibool_ocean_load(:) = 0
  normal_ocean_load(:,:) = 0._CUSTOM_REAL
  rmass_ocean_load_selected(:) = 0._CUSTOM_REAL

  ! fills arrays for coupling surface at oceans
  updated_dof_ocean_load(:) = .false.
  ipoin = 0
  do ispec2D = 1,NSPEC_TOP
    ispec = ibelm_top_crust_mantle(ispec2D)
    k = NGLLZ
    do j = 1,NGLLY
      do i = 1,NGLLX
        ! get global point number
        iglob = ibool_crust_mantle(i,j,k,ispec)
        ! updates once
        if (.not. updated_dof_ocean_load(iglob)) then
          ipoin = ipoin + 1
          ! fills arrays
          ibool_ocean_load(ipoin) = iglob
          rmass_ocean_load_selected(ipoin) = rmass_ocean_load(iglob)
          normal_ocean_load(:,ipoin) = normal_top_crust_mantle(:,i,j,ispec2D)
          ! masks this global point
          updated_dof_ocean_load(iglob) = .true.
        endif
      enddo
    enddo
  enddo

  ! frees memory
  deallocate(updated_dof_ocean_load)
  deallocate(rmass_ocean_load)

  ! synchronizes processes
  call synchronize_all()

  end subroutine prepare_oceans
