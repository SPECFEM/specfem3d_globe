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

  subroutine stretching_function(r_top,r_bottom,ner,stretch_tab)

! define stretch_tab which contains r_top and r_bottom for each element layer in the crust for 3D models.
!
! stretch_tab array uses indices index_radius & index_layer :
!   stretch_tab( index_radius (1=top,2=bottom) , index_layer (1=first layer, 2=second layer,..) )

  use shared_parameters, only: MAX_RATIO_CRUST_STRETCHING

  implicit none

  double precision,intent(in) :: r_top, r_bottom
  integer,intent(in) :: ner
  double precision, dimension (2,ner),intent(inout) :: stretch_tab

  ! local parameters
  double precision :: value
  ! for increasing execution speed but have less precision in stretching, increase step
  ! not very effective algorithm, but sufficient : used once per proc for meshing.
  double precision, parameter :: step = 0.001
  integer :: i

  ! safety check
  if (ner <= 1) stop 'Invalid ner value for stretching_function() routine'

  ! initializes array
  ! for example: 2 element layers (ner=2)  for most probable resolutions (NEX < 1000) in the crust
  !                      then stretch_tab(2,1) = 0.5 = stretch_tab(2,2)
  do i = 1,ner
    stretch_tab(2,i)=(1.d0/ner)
  enddo

  ! fill with ratio of the layer one thickness for each element
  do while((stretch_tab(2,1) / stretch_tab(2,ner)) > MAX_RATIO_CRUST_STRETCHING)
    if (modulo(ner,2) /= 0) then
      value = -floor(ner/2.d0)*step
    else
      value = (0.5d0-floor(ner/2.d0))*step
    endif
    do i = 1,ner
      stretch_tab(2,i) = stretch_tab(2,i) + value
      value = value + step
    enddo
  enddo

  ! deduce r_top and r_bottom
  ! r_top
  stretch_tab(1,1) = r_top
  do i=2,ner
    stretch_tab(1,i) = sum(stretch_tab(2,i:ner))*(r_top-r_bottom) + r_bottom
  enddo

  ! r_bottom
  stretch_tab(2,ner) = r_bottom
  do i = 1,ner-1
    stretch_tab(2,i) = stretch_tab(1,i+1)
  enddo

  end subroutine stretching_function

!
!-------------------------------------------------------------------------------------------------
!


  subroutine stretching_function_regional(r_top,r_bottom,ner,stretch_tab)

! define stretch_tab which contains r_top and r_bottom for each element layer in the crust for 3D models.
!
! stretch_tab array uses indices index_radius & index_layer :
!   stretch_tab( index_radius (1=top,2=bottom) , index_layer (1=first layer, 2=second layer,..) )

  use shared_parameters, only: PLANET_TYPE,IPLANET_MARS

  implicit none

  double precision,intent(in) :: r_top, r_bottom !,value
  integer,intent(in) :: ner !,i
  double precision, dimension (2,ner),intent(out) :: stretch_tab

  ! local parameters
  double precision :: R15,R35,R55,R80

!  ! for increasing execution speed but have less precision in stretching, increase step
!  ! not very effective algorithm, but sufficient : used once per proc for meshing.
!  double precision, parameter :: step = 0.001
!
!  ! initializes array
!  ! for example: 2 element layers (ner=2)  for most probable resolutions (NEX < 1000) in the crust
!  !                      then stretch_tab(2,1) = 0.5 = stretch_tab(2,2)
!  do i = 1,ner
!    stretch_tab(2,i)=(1.d0/ner)
!  enddo
!
!  ! fill with ratio of the layer one thickness for each element
!  do while((stretch_tab(2,1) / stretch_tab(2,ner)) > MAX_RATIO_CRUST_STRETCHING)
!    if (modulo(ner,2) /= 0) then
!      value = -floor(ner/2.d0)*step
!    else
!      value = (0.5d0-floor(ner/2.d0))*step
!    endif
!    do i = 1,ner
!      stretch_tab(2,i) = stretch_tab(2,i) + value
!      value = value + step
!    enddo
!  enddo
!
! deduce r_top and r_bottom
!  ! r_top
!  stretch_tab(1,1) = r_top
!  do i=2,ner
!    stretch_tab(1,i) = sum(stretch_tab(2,i:ner))*(r_top-r_bottom) + r_bottom
!  enddo
!
!  ! r_bottom
!  stretch_tab(2,ner) = r_bottom
!  do i = 1,ner-1
!    stretch_tab(2,i) = stretch_tab(1,i+1)
!  enddo

  ! saftey check
  if (ner /= 3 .and. ner /= 5) stop 'Error regional stretching function: ner value'

  ! radii for crustal layering
  if (PLANET_TYPE == IPLANET_MARS) then
    ! Mars
    ! surface radius 3390000.0
    if (ner == 3) then
      R15 = 3370000.d0  ! d = 20km ; 3353 second layer top  doesn't work well, give negative Jacobian
      R35 = 3307000.d0  ! d = 83km
    else
      R15 = 3375000.d0  ! d = 15km
      R35 = 3355000.d0  ! d = 35km
    endif
    R55 = 3335000.d0    ! d = 55km
    R80 = 3310000.d0    ! d = 80km
  else
    ! Earth
    ! surface radius 6371000.0
    R15 = 6356000.d0
    R35 = 6336000.d0
    R55 = 6326000.d0
    R80 = 6306000.d0
  endif

  if (ner == 3) then
    stretch_tab(1,1) = r_top
    stretch_tab(1,2) = R15  ! 15km second layer top
    stretch_tab(1,3) = R35  ! 35km third layer top

    stretch_tab(2,1) = stretch_tab(1,2)  ! bottom first layer
    stretch_tab(2,2) = stretch_tab(1,3)  ! bottom second layer
    stretch_tab(2,3) = r_bottom          ! bottom third layer

  else if (ner == 5) then
    stretch_tab(1,1) = r_top
    stretch_tab(1,2) = R15  ! 2nd layer top
    stretch_tab(1,3) = R35  ! 3rd layer top
    stretch_tab(1,4) = R55  ! 4th layer top
    stretch_tab(1,5) = R80  ! 5th layer top

    stretch_tab(2,1) = stretch_tab(1,2) ! 1st layer bottom
    stretch_tab(2,2) = stretch_tab(1,3) ! 2nd layer bottom
    stretch_tab(2,3) = stretch_tab(1,4) ! 3rd layer bottom
    stretch_tab(2,4) = stretch_tab(1,5)
    stretch_tab(2,5) = r_bottom

    if (r_bottom >= stretch_tab(1,5)) then
      print *, "ERROR: wrong layer in regional stretching table!",stretch_tab(1,5),'bottom',r_bottom
      stop 'Error regional stretching table'
    endif
  endif

  end subroutine stretching_function_regional


