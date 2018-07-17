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

  subroutine get_absorb(prname,iregion,iboun, &
                        nimin,nimax,njmin,njmax,nkmin_xi,nkmin_eta, &
                        NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM)

! Stacey, define flags for absorbing boundaries

  use constants
  use meshfem3D_par, only: ADIOS_FOR_ARRAYS_SOLVER,nspec

  implicit none

  integer,intent(in) :: iregion

  integer,intent(in) :: NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX,NSPEC2D_BOTTOM

  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nimin,nimax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX) :: njmin,njmax
  integer,dimension(2,NSPEC2DMAX_XMIN_XMAX) :: nkmin_xi
  integer,dimension(2,NSPEC2DMAX_YMIN_YMAX) :: nkmin_eta

  logical :: iboun(6,nspec)

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
      ispecb1=ispecb1+1

      ! this is useful even if it is constant because it can be zero inside the slices
      njmin(1,ispecb1) = 1
      njmax(1,ispecb1)=NGLLY

      !   check for overlap with other boundaries
      nkmin_xi(1,ispecb1) = 1
      if (iboun(5,ispec)) nkmin_xi(1,ispecb1) = 2
    endif

    if (iboun(2,ispec)) then

      !   on boundary 2: xmax
      ispecb2=ispecb2+1

      ! this is useful even if it is constant because it can be zero inside the slices
      njmin(2,ispecb2) = 1
      njmax(2,ispecb2)=NGLLY

      !   check for overlap with other boundaries
      nkmin_xi(2,ispecb2) = 1
      if (iboun(5,ispec)) nkmin_xi(2,ispecb2) = 2
    endif

    if (iboun(3,ispec)) then

      !   on boundary 3: ymin
      ispecb3=ispecb3+1

      !   check for overlap with other boundaries
      nimin(1,ispecb3) = 1
      if (iboun(1,ispec)) nimin(1,ispecb3) = 2
      nimax(1,ispecb3)=NGLLX
      if (iboun(2,ispec)) nimax(1,ispecb3)=NGLLX-1
      nkmin_eta(1,ispecb3) = 1
      if (iboun(5,ispec)) nkmin_eta(1,ispecb3) = 2
    endif

    if (iboun(4,ispec)) then

      !   on boundary 4: ymax
      ispecb4=ispecb4+1

      !   check for overlap with other boundaries
      nimin(2,ispecb4) = 1
      if (iboun(1,ispec)) nimin(2,ispecb4) = 2
      nimax(2,ispecb4)=NGLLX
      if (iboun(2,ispec)) nimax(2,ispecb4)=NGLLX-1
      nkmin_eta(2,ispecb4) = 1
      if (iboun(5,ispec)) nkmin_eta(2,ispecb4) = 2
    endif

    ! on boundary 5: bottom
    if (iboun(5,ispec)) ispecb5=ispecb5+1

  enddo

  ! check theoretical value of elements at the bottom
  if (ispecb5 /= NSPEC2D_BOTTOM) &
    call exit_MPI(myrank,'ispecb5 should equal NSPEC2D_BOTTOM in absorbing boundary detection')

  ! save these temporary arrays for the solver for Stacey conditions
  ! This files will be saved with the help of ADIOS if the
  ! ADIOS_FOR_ARRAYS_SOLVER flag is set to true in the Par_file
  if (ADIOS_FOR_ARRAYS_SOLVER) then
    call get_absorb_adios(iregion, &
                          nimin, nimax, njmin, njmax, nkmin_xi, nkmin_eta, &
                          NSPEC2DMAX_XMIN_XMAX,NSPEC2DMAX_YMIN_YMAX)
  else
    open(unit=IOUT,file=prname(1:len_trim(prname))//'stacey.bin', &
          status='unknown',form='unformatted',action='write',iostat=ier)
    if (ier /= 0 ) call exit_MPI(myrank,'Error opening stacey.bin file')
    write(IOUT) nimin
    write(IOUT) nimax
    write(IOUT) njmin
    write(IOUT) njmax
    write(IOUT) nkmin_xi
    write(IOUT) nkmin_eta
    close(IOUT)
  endif

  end subroutine get_absorb

