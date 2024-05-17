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


  subroutine interpolate_element_value(xi,eta,gamma,ielem, &
                                       nspec,model_elem, &
                                       val,xigll,yigll,zigll)

! interpolates model value for given location xi/eta/gamma within element ielem using GLL basis functions

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  double precision,intent(in):: xi,eta,gamma
  integer,intent(in):: ielem

  integer,intent(in):: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: model_elem

  real(kind=CUSTOM_REAL),intent(out) :: val

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  ! local parameters
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
                      hgammar(NGLLZ), hpgammar(NGLLZ)
  integer:: i,j,k

  ! checks
  if (ielem < 1 .or. ielem > nspec) stop 'Error invalid ielem index in interpolate_element_value()'

  ! interpolation weights
  call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates value
  val = 0.0
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
          val = val +  real(hxir(i) * hetar(j) * hgammar(k) * model_elem(i,j,k),kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  end subroutine interpolate_element_value

!
!------------------------------------------------------------------------------
!

  subroutine interpolate_limited(xi,eta,gamma,ielem, &
                                 nspec,model_elem, &
                                 val,xigll,yigll,zigll, &
                                 i_selected,j_selected,k_selected)

! interpolates model value for given location xi/eta/gamma within element ielem using GLL basis function;
! checks interpolated value with closest point value and limits value range to 1 % (PERCENTAGE_LIMIT)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLY,NGLLZ

  implicit none

  double precision,intent(in):: xi,eta,gamma
  integer,intent(in):: ielem

  integer,intent(in):: nspec
  real(kind=CUSTOM_REAL),dimension(NGLLX,NGLLY,NGLLZ),intent(in) :: model_elem

  real(kind=CUSTOM_REAL),intent(out) :: val

  ! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX),intent(in) :: xigll
  double precision, dimension(NGLLY),intent(in) :: yigll
  double precision, dimension(NGLLZ),intent(in) :: zigll

  integer,intent(in):: i_selected,j_selected,k_selected

  ! local parameters
  double precision :: hxir(NGLLX), hpxir(NGLLX), hetar(NGLLY), hpetar(NGLLY), &
                      hgammar(NGLLZ), hpgammar(NGLLZ)
  integer:: i,j,k
  real(kind=CUSTOM_REAL) :: val_initial,val_avg,pert,pert_limit

  !--------------------------------------------------------------------

  ! percentage
  real(kind=CUSTOM_REAL), parameter :: PERCENTAGE_LIMIT = 0.01

  !--------------------------------------------------------------------

  ! checks
  if (ielem < 1 .or. ielem > nspec) stop 'Error invalid ielem index in interpolate_limited()'

  ! interpolation weights
  call lagrange_any(xi,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(eta,NGLLY,yigll,hetar,hpetar)
  call lagrange_any(gamma,NGLLZ,zigll,hgammar,hpgammar)

  ! interpolates value
  val = 0.0_CUSTOM_REAL
  do k = 1, NGLLZ
    do j = 1, NGLLY
      do i = 1, NGLLX
        val = val + real(hxir(i) * hetar(j) * hgammar(k) * model_elem(i,j,k),kind=CUSTOM_REAL)
      enddo
    enddo
  enddo

  ! note: interpolation of values close to the surface or 3D moho encounters problems;
  !       this is a fall-back to the closest point value
  !
  ! uses average/closest point value if too far off

  ! closest point value
  val_initial = model_elem(i_selected,j_selected,k_selected)

  ! average value
  val_avg = sum(model_elem(:,:,:)) / NGLLX / NGLLY / NGLLZ

  ! initial difference
  pert = abs(val_initial - val_avg)

  ! upper/lower perturbation bound
  pert_limit = PERCENTAGE_LIMIT * abs(val_avg)
  if (pert > pert_limit) pert_limit = pert

  ! within a certain percentage range
  if (abs(val - val_avg ) > pert_limit) then
    val = val_initial
  endif

  end subroutine interpolate_limited
